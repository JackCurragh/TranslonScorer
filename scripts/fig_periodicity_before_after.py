"""Talk figure (slide 9b): per-sample P-site calibration recovers periodicity
that a flat-offset aggregate destroys.

Builds the FrameRollup over a canonical-CDS feature set from the matrix, then
scores frame-0 fraction two ways on the SAME reads:
  - flat offset = 15 (the legacy/broken aggregate path)
  - per-(sample, length) calibrated offsets (the fixed FrameRollup path)

Produces a two-panel figure:
  Left  — frame-0 fraction by read length, flat vs calibrated, 1/3 floor line.
  Right — headline aggregate over dominant RPF lengths (flat ~34% -> calibrated ~87%).

This is the same machinery validated in scripts/step1_validate_frame_rollup.py;
here it is rendered for the deck and generalised to a CDS *set* so the claim is
not a single cherry-picked gene.

Local smoke test (30-sample GAPDH fixture):
  python3 scripts/fig_periodicity_before_after.py \
      --partition-dir data/global_partitioned --gtf data/genes.gtf \
      --transcripts ENST00000229239 --out outputs/fig_periodicity_gapdh.png

HPC (canonical CDS set — recommended for the talk; pass a MANE/principal id list):
  python3 scripts/fig_periodicity_before_after.py \
      --partition-dir /hps/.../global_partitioned \
      --gtf /path/to/annotation.gtf \
      --transcript-file canonical_chr22_mane.txt \
      --out outputs/fig_periodicity_chr22.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")  # headless / HPC
import matplotlib.pyplot as plt  # noqa: E402
import polars as pl  # noqa: E402

from TranslonScorer.io.annotation import build_cds_blocks  # noqa: E402
from TranslonScorer.matrix_rollup import (  # noqa: E402
    build_frame_rollup,
    calibrate_offsets,
    score_frame_rollup,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--partition-dir", required=True, type=Path,
                   help="global_partitioned matrix directory")
    p.add_argument("--gtf", required=True, type=Path, help="annotation GTF")
    p.add_argument("--transcripts", default="ENST00000229239",
                   help="comma-separated transcript ids (default: GAPDH)")
    p.add_argument("--transcript-file", type=Path,
                   help="file with one transcript id per line (overrides --transcripts)")
    p.add_argument("--out", required=True, type=Path, help="output PNG path")
    p.add_argument("--n-workers", type=int, default=4)
    p.add_argument("--dominant", default="28,29,30",
                   help="dominant RPF lengths for the headline aggregate")
    p.add_argument("--min-reads", type=int, default=50,
                   help="min reads per length to plot in the per-length panel")
    return p.parse_args()


def resolve_transcripts(args: argparse.Namespace) -> list[str]:
    if args.transcript_file:
        ids = [ln.strip() for ln in args.transcript_file.read_text().splitlines() if ln.strip()]
    else:
        ids = [t.strip() for t in args.transcripts.split(",") if t.strip()]
    if not ids:
        sys.exit("No transcript ids resolved.")
    return ids


def per_length_table(scored: pl.DataFrame) -> pl.DataFrame:
    """Aggregate per-sample scored rows to one frame-0 fraction per length."""
    return (
        scored.group_by("length")
        .agg(
            pl.col("frame0_count").sum().alias("frame0"),
            pl.col("n_reads").sum().alias("n_reads"),
        )
        .with_columns((pl.col("frame0") / pl.col("n_reads")).alias("f0"))
        .sort("length")
    )


def main() -> None:
    args = parse_args()
    dom = [int(x) for x in args.dominant.split(",")]
    tx_ids = resolve_transcripts(args)

    print(f"Building FrameRollup over {len(tx_ids)} transcript(s)...")
    cds = build_cds_blocks(str(args.gtf))
    feat = cds.filter(pl.col("tran_id").is_in(tx_ids))
    if feat.height == 0:
        sys.exit(f"None of the requested transcripts found in {args.gtf}")
    if feat.height < len(tx_ids):
        print(f"  WARNING: only {feat.height}/{len(tx_ids)} transcripts present in GTF")

    rollup = build_frame_rollup(
        args.partition_dir, feat, multimap_mode="unique", n_workers=args.n_workers
    )
    print(f"  rollup rows={rollup.shape[0]}  samples={rollup['sample_name'].n_unique()}  "
          f"reads={rollup['count'].sum():.0f}")

    agg_for_cal = rollup.group_by(["sample_name", "length", "strand", "phase0"]).agg(
        pl.col("count").sum()
    )
    offsets = calibrate_offsets(agg_for_cal, target_frame=0)

    scored_cal = score_frame_rollup(rollup, offsets, default_offset=12)
    scored_flat = score_frame_rollup(rollup, {}, default_offset=15)

    cal = per_length_table(scored_cal)
    flat = per_length_table(scored_flat)

    merged = cal.join(flat, on="length", how="inner", suffix="_flat").filter(
        (pl.col("n_reads") >= args.min_reads) & pl.col("length").is_in(range(24, 38))
    )

    # headline aggregate over dominant lengths
    cd = scored_cal.filter(pl.col("length").is_in(dom))
    fd = scored_flat.filter(pl.col("length").is_in(dom))
    f0_cal = cd["frame0_count"].sum() / cd["n_reads"].sum()
    f0_flat = fd["frame0_count"].sum() / fd["n_reads"].sum()
    print(f"\nHEADLINE (lengths {dom}):  flat-15 = {f0_flat:.1%}   calibrated = {f0_cal:.1%}")

    # ---- plot ----
    lengths = merged["length"].to_list()
    f0c = merged["f0"].to_list()
    f0f = merged["f0_flat"].to_list()

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2),
                                   gridspec_kw={"width_ratios": [3, 1]})

    x = range(len(lengths))
    w = 0.4
    axL.bar([i - w / 2 for i in x], f0f, width=w, label="flat offset = 15",
            color="#bdbdbd")
    axL.bar([i + w / 2 for i in x], f0c, width=w, label="calibrated per (sample, length)",
            color="#2b8cbe")
    axL.axhline(1 / 3, ls="--", lw=1, color="#999999")
    axL.text(len(lengths) - 0.5, 1 / 3 + 0.01, "random 3-frame floor",
             ha="right", va="bottom", fontsize=8, color="#666666")
    axL.set_xticks(list(x))
    axL.set_xticklabels(lengths)
    axL.set_xlabel("read length (nt)")
    axL.set_ylabel("in-frame (frame-0) fraction")
    axL.set_ylim(0, 1)
    axL.set_title("Frame-0 fraction by read length")
    axL.legend(frameon=False, fontsize=8, loc="upper left")

    bars = axR.bar(["flat-15", "calibrated"], [f0_flat, f0_cal],
                   color=["#bdbdbd", "#2b8cbe"])
    axR.axhline(1 / 3, ls="--", lw=1, color="#999999")
    axR.set_ylim(0, 1)
    axR.set_ylabel("in-frame fraction")
    axR.set_title(f"Canonical CDS\n(lengths {'/'.join(map(str, dom))} nt)")
    for b, v in zip(bars, [f0_flat, f0_cal]):
        axR.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.0%}",
                 ha="center", va="bottom", fontsize=11, fontweight="bold")

    n_samples = rollup["sample_name"].n_unique()
    fig.suptitle(
        f"Per-sample P-site calibration recovers periodicity a flat offset destroys "
        f"({len(tx_ids)} CDS, {n_samples} samples)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
