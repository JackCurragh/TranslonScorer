import json
from pathlib import Path

import polars as pl
import pysam

from TranslonScorer.coverage.profiles import (
    gene_expression_matrix_from_profiles,
    profiles_from_sparse_parquet_matrix,
)


def _write_bam(path: Path) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for query_name, start in [("read_0", 10), ("read_1", 20)]:
            record = pysam.AlignedSegment()
            record.query_name = query_name
            record.query_sequence = "A" * 28
            record.flag = 0
            record.reference_id = 0
            record.reference_start = start
            record.mapping_quality = 255
            record.cigartuples = [(0, 28)]
            record.query_qualities = pysam.qualitystring_to_array("I" * 28)
            bam.write(record)
    pysam.index(str(path))


def _write_sparse_store(root: Path) -> Path:
    bucket = root / "global_counts" / "generation=000001" / "read_bucket=000000"
    bucket.mkdir(parents=True)
    pl.DataFrame(
        {
            "read_bucket": [0, 0, 0],
            "read_id": [0, 0, 1],
            "sample_id": [0, 1, 1],
            "study_id_int": [0, 1, 1],
            "count": [2, 5, 7],
        }
    ).write_parquet(bucket / "part.parquet")
    pl.DataFrame(
        {
            "sample_id": [0, 1],
            "sample_name": ["run_a", "run_b"],
            "study_id": ["study_a", "study_b"],
            "study_id_int": [0, 1],
        }
    ).write_parquet(root / "global_samples.parquet")
    manifest = {
        "matrix_format": "sparse-parquet",
        "read_bucket_size": 100,
        "lookup_tables": {"samples": "global_samples.parquet"},
        "generations": [{"generation": 1, "path": "global_counts/generation=000001"}],
    }
    manifest_path = root / "global_matrix_manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    return manifest_path


def test_sparse_parquet_profiles_preserve_runs_and_gene_matrix(tmp_path):
    bam = tmp_path / "global.bam"
    _write_bam(bam)
    manifest = _write_sparse_store(tmp_path / "matrix")

    exons = pl.DataFrame(
        {
            "tran_id": ["tx1"],
            "chr": ["chr1"],
            "start": [[0]],
            "stop": [[100]],
            "tran_start": [[0]],
            "tran_stop": [[100]],
            "strand": ["+"],
        }
    )
    cds = pl.DataFrame(
        {
            "tran_id": ["tx1"],
            "chr": ["chr1"],
            "start": [0],
            "stop": [90],
            "strand": ["+"],
        }
    )

    profiles, offsets, genomic_counts = profiles_from_sparse_parquet_matrix(
        bam_path=str(bam),
        manifest_path=str(manifest),
        exon_df=exons,
        cds_df=cds,
        default_offset=15,
    )

    assert set(genomic_counts.get_column("sample_id").to_list()) == {"run_a", "run_b"}
    assert offsets
    assert profiles.select("sample_id").unique().height == 2

    offset = offsets[28]
    profile_rows = profiles.sort(["sample_id", "tran_id", "pos"]).to_dicts()
    assert profile_rows == [
        {
            "sample_id": "run_a",
            "sample_index": 0,
            "study_id": "study_a",
            "study_id_int": 0,
            "tran_id": "tx1",
            "pos": 10 + offset,
            "count": 2.0,
        },
        {
            "sample_id": "run_b",
            "sample_index": 1,
            "study_id": "study_b",
            "study_id_int": 1,
            "tran_id": "tx1",
            "pos": 10 + offset,
            "count": 5.0,
        },
        {
            "sample_id": "run_b",
            "sample_index": 1,
            "study_id": "study_b",
            "study_id_int": 1,
            "tran_id": "tx1",
            "pos": 20 + offset,
            "count": 7.0,
        },
    ]

    transcripts = pl.DataFrame({"tran_id": ["tx1"], "gene_id": ["gene1"]})
    long, wide = gene_expression_matrix_from_profiles(
        profiles,
        exon_df=exons,
        transcripts_df=transcripts,
    )

    assert long.sort(["gene_id", "sample_id"]).to_dicts() == [
        {"gene_id": "gene1", "sample_id": "run_a", "count": 2.0},
        {"gene_id": "gene1", "sample_id": "run_b", "count": 12.0},
    ]
    assert wide.to_dicts() == [{"gene_id": "gene1", "run_a": 2.0, "run_b": 12.0}]
