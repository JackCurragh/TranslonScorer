"""Calibration benchmark: do the scores separate coding from non-coding?

Scores annotated CDS (positive) and lncRNA exons (negative) from the SAME
annotation through the SAME pipeline and the same coverage, then reports the
separation.  Two things this catches that unit tests cannot:

  * a silent-absence bug — coverage that stops reaching events shows up as the
    positive set collapsing toward the negative one;
  * a miscalibrated metric — elong_in_frame has a known null (1/3, no
    periodicity), so the negative set landing anywhere other than ~0.33 means
    the metric is not measuring what it claims.

lncRNA is a SOFT negative: some lncRNAs are genuinely translated, which is much
of why translon annotation exists.  Treat the negative-set SUPPORTED fraction as
an upper bound on the false-positive rate, and those events as candidates rather
than errors.

Depth matching matters — lncRNA is lower-coverage than CDS, and that alone can
manufacture a separation.  The report includes a depth-matched comparison; trust
that one.

Usage:
  python scripts/benchmark_cds_vs_lncrna.py \
      --gtf data/genes.gtf --chrom 12 \
      --forward-bigwig fwd.bw --reverse-bigwig rev.bw \
      --out-dir bench_out/
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import polars as pl  # noqa: E402

from TranslonScorer.workflows import pipeline_workflow  # noqa: E402

MIN_READS_FOR_DEPTH_MATCH = 5000


def build_lncrna_gtf(gtf: Path, chrom: str, out: Path) -> int:
    """Extract lncRNA exons for one chromosome.

    Uses awk with an explicit TAB field separator: the default separator splits
    GTF column 9 on spaces and rejoins with tabs, which silently corrupts the
    attributes and yields null transcript_ids downstream.
    """
    prog = (
        f'BEGIN{{OFS="\\t"}} $1=="{chrom}" && $3=="exon" ' f'&& /gene_biotype "lncRNA"/ {{print}}'
    )
    with out.open("w") as fh:
        subprocess.run(["awk", "-F", "\t", prog, str(gtf)], stdout=fh, check=True)
    return sum(1 for _ in out.open())


def score_set(out_dir: Path, gtf: Path, feature_type: str, bigwigs, stranded: bool, chrom: str):
    pipeline_workflow(
        str(out_dir),
        bigwigs=bigwigs,
        stranded=stranded,
        gtf_path=str(gtf),
        feature_type=feature_type,
        chroms=[chrom],
        data_version="bench",
    )
    files = list((out_dir / "scores").rglob("*.parquet"))
    if not files:
        raise SystemExit(f"no scores written under {out_dir}")
    return pl.read_parquet(files)


def elongation(df: pl.DataFrame, label: str) -> pl.DataFrame:
    return df.filter((pl.col("aspect") == "elongation") & (pl.col("n_reads") > 0)).select(
        pl.lit(label).alias("set"), "metric", "n_reads"
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", required=True, type=Path)
    p.add_argument("--chrom", required=True)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--bigwig", action="append", default=[])
    p.add_argument("--forward-bigwig", action="append", default=[])
    p.add_argument("--reverse-bigwig", action="append", default=[])
    args = p.parse_args()

    if len(args.forward_bigwig) != len(args.reverse_bigwig):
        raise SystemExit("--forward-bigwig and --reverse-bigwig must be given in pairs")
    stranded = bool(args.forward_bigwig)
    bigwigs = (
        [{"forward": f, "reverse": r} for f, r in zip(args.forward_bigwig, args.reverse_bigwig)]
        if stranded
        else args.bigwig
    )
    if not bigwigs:
        raise SystemExit("provide --bigwig or --forward-bigwig/--reverse-bigwig")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    neg_gtf = args.out_dir / f"chr{args.chrom}_lncRNA.gtf"
    n_lines = build_lncrna_gtf(args.gtf, args.chrom, neg_gtf)
    print(f"negative set: {n_lines} lncRNA exon lines -> {neg_gtf}")

    pos = score_set(args.out_dir / "positive_cds", args.gtf, "CDS", bigwigs, stranded, args.chrom)
    neg = score_set(
        args.out_dir / "negative_lncrna", neg_gtf, "exon", bigwigs, stranded, args.chrom
    )

    both = pl.concat([elongation(pos, "CDS (positive)"), elongation(neg, "lncRNA (negative)")])
    summary = both.group_by("set").agg(
        pl.len().alias("n"),
        pl.col("metric").median().round(3).alias("median"),
        (pl.col("metric") > 0.5).mean().round(3).alias("frac>0.5"),
    )
    matched = (
        both.filter(pl.col("n_reads") >= MIN_READS_FOR_DEPTH_MATCH)
        .group_by("set")
        .agg(
            pl.len().alias("n"),
            pl.col("metric").median().round(3).alias("median"),
            (pl.col("metric") > 0.5).mean().round(3).alias("frac>0.5"),
        )
    )
    print("\n=== all covered elongation events ===")
    print(summary.sort("set"))
    print(f"\n=== depth-matched (n_reads >= {MIN_READS_FOR_DEPTH_MATCH}) — trust this one ===")
    print(matched.sort("set"))

    print("\n=== elongation calls ===")
    for df, label in ((pos, "CDS"), (neg, "lncRNA")):
        counts = df.filter(pl.col("aspect") == "elongation")["call"].value_counts()
        total = counts["count"].sum()
        print(f"  {label:7s}", {r["call"]: round(r["count"] / total, 3) for r in counts.to_dicts()})

    med = {r["set"]: r["median"] for r in matched.to_dicts()}
    pos_med = med.get("CDS (positive)")
    neg_med = med.get("lncRNA (negative)")
    if pos_med is not None and neg_med is not None:
        print(
            f"\ndepth-matched separation: {pos_med:.3f} vs {neg_med:.3f} (null for no periodicity = 0.333)"
        )
        if neg_med > 0.45:
            print("  WARNING: negative set well above the 1/3 null — metric may be miscalibrated.")
        if pos_med - neg_med < 0.15:
            print("  WARNING: weak separation — check coverage is reaching the positive set.")


if __name__ == "__main__":
    main()
