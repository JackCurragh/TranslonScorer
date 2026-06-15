"""Feature-source adapters: GTF/GFF, BED12, FASTA → (blocks, translons) → events.

These prove you can score features without an annotation sqlite DB.
"""

from __future__ import annotations

import polars as pl
import pytest

from TranslonScorer.io import feature_sources as fs

# ---------------------------------------------------------------------------
# canonical schema reshape (+ and - strand, exon order)
# ---------------------------------------------------------------------------


def _gtf(tmp_path) -> str:
    # + strand CDS: exons 101-130 then 201-230 (genomic 100-130, 200-230)
    # - strand CDS: exons 401-430 then 501-530
    lines = [
        'chr1\ts\tCDS\t101\t130\t.\t+\t0\tgene_id "G1"; transcript_id "TXP";',
        'chr1\ts\tCDS\t201\t230\t.\t+\t0\tgene_id "G1"; transcript_id "TXP";',
        'chr1\ts\tCDS\t401\t430\t.\t-\t0\tgene_id "G2"; transcript_id "TXM";',
        'chr1\ts\tCDS\t501\t530\t.\t-\t0\tgene_id "G2"; transcript_id "TXM";',
    ]
    p = tmp_path / "a.gtf"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_from_gtf_blocks_and_translons(tmp_path):
    blocks, translons = fs.from_gtf(_gtf(tmp_path), feature_type="CDS")
    assert set(blocks.columns) == {
        "translon_id",
        "translation_block_rank",
        "bed_chrom",
        "bed_start",
        "bed_end",
        "seq_region_strand",
        "block_length_nt",
    }
    txp = blocks.filter(pl.col("translon_id") == "TXP").sort("translation_block_rank")
    # + strand, transcription order = genomic L->R; 0-based half-open
    assert txp["bed_start"].to_list() == [100, 200]
    assert txp["seq_region_strand"].to_list() == [1, 1]
    assert txp["block_length_nt"].to_list() == [30, 30]
    # - strand transcription order = genomic R->L
    txm = blocks.filter(pl.col("translon_id") == "TXM").sort("translation_block_rank")
    assert txm["bed_start"].to_list() == [500, 400]
    assert txm["seq_region_strand"].to_list() == [-1, -1]
    # translon span = min start .. max stop
    tp = translons.filter(pl.col("translon_id") == "TXP").row(0, named=True)
    assert tp["bed_start"] == 100 and tp["bed_end"] == 230


def test_from_gtf_feeds_extract_events(tmp_path):
    from TranslonScorer.events import extract_events

    blocks, translons = fs.from_gtf(_gtf(tmp_path), feature_type="CDS")
    events, fe, _ = extract_events(blocks, translons)
    assert {"init", "elongation", "term", "junction"} <= set(events["type"])
    assert set(fe["feature_id"]) == {"TXP", "TXM"}


# ---------------------------------------------------------------------------
# BED12
# ---------------------------------------------------------------------------


def test_from_bed12(tmp_path):
    # one 2-block + strand feature: chr1 100-230, blocks 0-30 and 100-130
    bed = tmp_path / "f.bed"
    bed.write_text("chr1\t100\t230\torfA\t0\t+\t100\t230\t0\t2\t30,30\t0,100\n")
    blocks, translons = fs.from_bed12(str(bed))
    b = blocks.filter(pl.col("translon_id") == "orfA").sort("translation_block_rank")
    assert b["bed_start"].to_list() == [100, 200]
    assert b["bed_end"].to_list() == [130, 230]
    assert b["seq_region_strand"].to_list() == [1, 1]


def test_from_bed12_rejects_short_lines(tmp_path):
    bed = tmp_path / "bad.bed"
    bed.write_text("chr1\t100\t230\tx\t0\t+\n")  # only 6 cols
    with pytest.raises(ValueError):
        fs.from_bed12(str(bed))


# ---------------------------------------------------------------------------
# de-novo FASTA (both strands, single-block genomic ORFs)
# ---------------------------------------------------------------------------


def test_from_fasta_denovo_both_strands(tmp_path):
    pytest.importorskip("pyfaidx")
    # Forward ORF: ATG ... TAA. Build a sequence with a clean + strand ORF.
    seq = "ATG" + "AAA" * 5 + "TAA" + "CCCCC"
    fa = tmp_path / "g.fa"
    fa.write_text(f">contig1\n{seq}\n")
    blocks, translons = fs.from_fasta(str(fa), min_len=0, max_len=1000)
    assert not blocks.is_empty()
    assert set(blocks["bed_chrom"].unique()) == {"contig1"}
    # strands present are +/-1; the forward ATG..TAA ORF should appear on +
    assert 1 in set(blocks["seq_region_strand"].to_list())


def test_from_fasta_empty_when_no_orfs(tmp_path):
    pytest.importorskip("pyfaidx")
    fa = tmp_path / "n.fa"
    fa.write_text(">c\nCCCCCCCCCCCC\n")  # no ATG
    blocks, translons = fs.from_fasta(str(fa))
    assert blocks.is_empty()


# ---------------------------------------------------------------------------
# workflow dispatch guard
# ---------------------------------------------------------------------------


def test_workflow_requires_exactly_one_source(tmp_path):
    from TranslonScorer.workflows import extract_events_workflow

    with pytest.raises(ValueError):
        extract_events_workflow(str(tmp_path / "out"))  # none
    with pytest.raises(ValueError):
        extract_events_workflow(str(tmp_path / "out"), gtf_path="a.gtf", bed12_path="b.bed")  # two


def test_extract_events_workflow_from_gtf_end_to_end(tmp_path):
    from TranslonScorer.workflows import extract_events_workflow

    out = tmp_path / "ev"
    summary = extract_events_workflow(str(out), gtf_path=_gtf(tmp_path), feature_type="CDS")
    assert summary["events"] > 0
    assert (out / "events").exists() and (out / "feature_event").exists()
