"""Self-contained end-to-end smoke tests for the real scoring pipeline.

Unlike the data-backed tests (test_bam_provider / test_frame_rollup / psite),
these generate every fixture in tmp_path (tiny FASTA + bigwig via pyBigWig, tiny
BAM via pysam), so they exercise genome/sequence -> candidate ORF extraction ->
coverage provider -> scoring INSIDE CI, where the gitignored data/ fixtures are
absent. The bigwig test's ORF ends at the contig edge on purpose, so it also
guards the out-of-bounds clamp fix. (The annotation fail-loud guard is covered
separately in test_new_tree_units.)
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import polars as pl
import pytest

from TranslonScorer.workflows import (
    extract_events_workflow,
    score_bams_workflow,
    score_bigwigs_workflow,
)

# A tiny transcript that IS a single ORF: ATG + 40 Ala codons + stop. The ORF
# runs to the contig end, so the padded scoring window overruns the transcript
# length — the exact case the bigwig clamp fix handles.
_ORF_SEQ = "ATG" + "GCT" * 40 + "TAA"


def _read_store(store_dir: str) -> pl.DataFrame:
    files = glob.glob(f"{store_dir}/**/*.parquet", recursive=True)
    assert files, f"no score parquet written under {store_dir}"
    return pl.read_parquet(files)


def _elong_covered(df: pl.DataFrame) -> int:
    """Count elongation events with covered_nt > 0 (coverage actually reached them)."""
    elong = df.filter(pl.col("aspect") == "elongation")
    return sum(1 for e in elong["evidence"].to_list() if json.loads(e).get("covered_nt", 0) > 0)


def test_e2e_bigwig_denovo_orf(tmp_path: Path):
    """genome/sequence -> de-novo ORF extraction -> bigwig scoring, self-contained."""
    pyBigWig = pytest.importorskip("pyBigWig")

    fasta = tmp_path / "tran.fa"
    fasta.write_text(f">tranA\n{_ORF_SEQ}\n")

    bw_path = tmp_path / "cov.bw"
    n = len(_ORF_SEQ)
    bw = pyBigWig.open(str(bw_path), "w")
    bw.addHeader([("tranA", n)])
    # dense coverage over the whole transcript (value 5) -> well above min_reads
    bw.addEntries(["tranA"] * n, list(range(n)), ends=list(range(1, n + 1)), values=[5.0] * n)
    bw.close()

    events_dir = tmp_path / "events"
    extract_events_workflow(str(events_dir), fasta_path=str(fasta), start_codons=["ATG"])

    store_dir = tmp_path / "store"
    score_bigwigs_workflow(str(events_dir), [str(bw_path)], str(store_dir), data_version="v1")
    df = _read_store(str(store_dir))
    assert df.height > 0
    # coverage reached the elongation event(s) despite the ORF ending at the
    # contig edge (regression: bigwig out-of-bounds clamp).
    assert _elong_covered(df) > 0


def test_e2e_bam_gtf(tmp_path: Path):
    """candidate CDS (GTF) -> BAM scoring, self-contained."""
    pysam = pytest.importorskip("pysam")
    from TranslonScorer.model import OffsetParams

    # 300 nt CDS on a small contig (1-based GTF coords).
    gtf = tmp_path / "anno.gtf"
    gtf.write_text('chrT\ttest\tCDS\t10\t309\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n')

    # ~97 unique 30 nt reads tiled in-frame across the CDS.
    unsorted = tmp_path / "reads.unsorted.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chrT", "LN": 1000}]}
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for i, pos in enumerate(range(9, 300, 3)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * 30
            a.flag = 0
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 255
            a.cigartuples = [(0, 30)]  # 30M
            a.query_qualities = pysam.qualitystring_to_array("I" * 30)
            a.set_tag("NH", 1)  # unique
            out.write(a)
    bam = tmp_path / "reads.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))

    events_dir = tmp_path / "events"
    extract_events_workflow(str(events_dir), gtf_path=str(gtf), feature_type="CDS")

    store_dir = tmp_path / "store"
    score_bams_workflow(
        str(events_dir),
        [str(bam)],
        str(store_dir),
        data_version="v1",
        offsets=OffsetParams(method="global", global_offset=12),
    )
    df = _read_store(str(store_dir))
    assert df.height > 0
    assert _elong_covered(df) > 0
