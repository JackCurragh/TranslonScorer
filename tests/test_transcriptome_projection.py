"""T12.1 — transcript→genome projection: pure helpers + BamSetProvider path.

The integration test builds a tiny transcriptome-aligned BAM (references are
transcript ids) and checks that BamSetProvider projects reads to the right
genomic P/A-site, including isoform-multimapper convergence on a shared exon.
"""
from __future__ import annotations

import polars as pl
import pytest

pysam = pytest.importorskip("pysam")


# Two-exon + strand transcript and a - strand transcript, plus a second isoform
# of TX1 that shares the first exon (for multimapper convergence).
def _exon_df() -> pl.DataFrame:
    return pl.DataFrame({
        "tran_id": ["TX1", "TX2", "TX1b"],
        "chr": ["chr1", "chr2", "chr1"],
        "strand": ["+", "-", "+"],
        "start": [[100, 200], [500, 400], [100, 300]],
        "stop": [[130, 230], [520, 420], [130, 330]],
        "tran_start": [[0, 30], [0, 20], [0, 30]],
    })


# ---------------------------------------------------------------------------
# pure projection
# ---------------------------------------------------------------------------

def test_build_and_project_plus_strand():
    from TranslonScorer.coverage.transcriptome import build_exon_index, project_to_genome
    idx = build_exon_index(_exon_df())
    tx1 = idx["TX1"]
    assert tx1.chrom == "chr1" and tx1.strand == 1
    # transcript pos 12 -> exon0 [0,30): 100 + 12
    assert project_to_genome(tx1, 12) == (1, 112)
    # transcript pos 35 -> exon1 (tran_start 30): 200 + (35-30) = 205
    assert project_to_genome(tx1, 35) == (1, 205)
    # out of range
    assert project_to_genome(tx1, 999) is None


def test_project_minus_strand():
    from TranslonScorer.coverage.transcriptome import build_exon_index, project_to_genome
    idx = build_exon_index(_exon_df())
    tx2 = idx["TX2"]
    assert tx2.strand == -1
    # exon0 (ts0): genomic 5' end at stop-1 = 519, descending. tc=5 -> 519-5
    assert project_to_genome(tx2, 5) == (-1, 514)
    # exon1 (ts 20): stop 420 -> 419 descending. tc=25 -> 419 - (25-20)
    assert project_to_genome(tx2, 25) == (-1, 414)


def test_build_exon_index_missing_columns():
    from TranslonScorer.coverage.transcriptome import build_exon_index
    with pytest.raises(ValueError):
        build_exon_index(pl.DataFrame({"tran_id": ["X"]}))


# ---------------------------------------------------------------------------
# BamSetProvider transcriptome integration
# ---------------------------------------------------------------------------

def _write_transcriptome_bam(path, records, refs):
    """records: list of (qname, ref_name, tran_start, length)."""
    header = {"HD": {"VN": "1.0"},
              "SQ": [{"SN": name, "LN": ln} for name, ln in refs]}
    name_to_id = {name: i for i, (name, _) in enumerate(refs)}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for qname, ref, tstart, length in records:
            a = pysam.AlignedSegment(bam.header)
            a.query_name = qname
            a.flag = 0
            a.reference_id = name_to_id[ref]
            a.reference_start = tstart
            a.mapping_quality = 255
            a.cigartuples = [(0, length)]  # length M
            a.query_sequence = "A" * length
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            a.set_tag("NH", 1)
            bam.write(a)
    pysam.index(str(path))


def test_transcriptome_coverage_projects_to_genome(tmp_path):
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import Region, OffsetParams

    bam = tmp_path / "tx.bam"
    # Read r1 on TX1 at transcript pos 0, length 30. P-site offset 12 ->
    # transcript coord 12 -> genomic 112 (+ strand).
    _write_transcriptome_bam(bam, [("r1", "TX1", 0, 30)], [("TX1", 60), ("TX2", 40), ("TX1b", 60)])

    provider = BamSetProvider(
        [str(bam)], exon_df=_exon_df(), transcriptome=True,
        offsets=OffsetParams(method="global", global_offset=12),
    )
    cov = provider.coverage([Region("chr1", 100, 230)], site="P")
    d = dict(zip(cov["pos"].to_list(), cov["count"].to_list()))
    assert d == {112: 1.0}
    # A-site = +3 nt -> 115
    covA = provider.coverage([Region("chr1", 100, 230)], site="A")
    assert dict(zip(covA["pos"].to_list(), covA["count"].to_list())) == {115: 1.0}


def test_transcriptome_isoform_multimapper_convergence(tmp_path):
    """A read hitting TX1 and TX1b (shared first exon) at the same transcript
    offset projects to ONE genomic site -> counted once (unique)."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import Region, OffsetParams

    bam = tmp_path / "tx2.bam"
    # Same qname aligned to both isoforms at pos 0; first exon is shared
    # (chr1:100-130) so both project to genomic 112.
    _write_transcriptome_bam(
        bam,
        [("r1", "TX1", 0, 30), ("r1", "TX1b", 0, 30)],
        [("TX1", 60), ("TX2", 40), ("TX1b", 60)],
    )
    provider = BamSetProvider(
        [str(bam)], exon_df=_exon_df(), transcriptome=True,
        offsets=OffsetParams(method="global", global_offset=12),
    )
    cov = provider.coverage([Region("chr1", 100, 230)], site="P")
    # Convergent => counted once, not twice.
    assert dict(zip(cov["pos"].to_list(), cov["count"].to_list())) == {112: 1.0}


def test_transcriptome_divergent_multimapper_dropped_when_unique(tmp_path):
    """A read hitting two transcripts that project to DIFFERENT genomic sites is
    dropped under multimap='unique'."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import Region, OffsetParams

    bam = tmp_path / "tx3.bam"
    # r1 on TX1 (->chr1:112) and TX2 (->chr2, different site): divergent.
    _write_transcriptome_bam(
        bam,
        [("r1", "TX1", 0, 30), ("r1", "TX2", 0, 30)],
        [("TX1", 60), ("TX2", 40), ("TX1b", 60)],
    )
    provider = BamSetProvider(
        [str(bam)], exon_df=_exon_df(), transcriptome=True, multimap="unique",
        offsets=OffsetParams(method="global", global_offset=12),
    )
    cov = provider.coverage([Region("chr1", 100, 230)], site="P")
    assert cov.is_empty()
