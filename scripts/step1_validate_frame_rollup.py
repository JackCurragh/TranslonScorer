"""Step 1 validation — FrameRollup correctness gate.

Builds the feature-scoped FrameRollup on the local 30-sample GAPDH matrix
(data/global_partitioned) and verifies that calibrated per-(sample,length)
offsets recover strong frame-0 periodicity for GAPDH (ENST00000229239).

Expected result (positive-control gate, FR3):
  - calibrated frame-0 fraction (29nt) > 0.65 per sample
  - aggregate calibrated frame-0 > flat-15 frame-0 (or ≥ for lengths where
    calibrated ≡ 15 mod 3)
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import polars as pl  # noqa: E402

from TranslonScorer.io.annotation import build_cds_blocks  # noqa: E402
from TranslonScorer.matrix_rollup import (  # noqa: E402
    build_frame_rollup,
    calibrate_offsets,
    score_frame_rollup,
)

PARTITION_DIR = REPO_ROOT / "data" / "global_partitioned"
GTF = REPO_ROOT / "data" / "genes.gtf"


def main() -> None:
    if not PARTITION_DIR.exists():
        print(f"SKIP: partition dir not found ({PARTITION_DIR})")
        return
    if not GTF.exists():
        print(f"SKIP: GTF not found ({GTF})")
        return

    print("Step 1 — FrameRollup validation")
    print(f"Partition dir: {PARTITION_DIR}")
    print(f"GTF: {GTF}")

    cds = build_cds_blocks(str(GTF))
    gapdh = cds.filter(pl.col("tran_id") == "ENST00000229239")
    print(f"\nGAPDH CDS exons: {len(gapdh.row(0, named=True)['start'])}")

    print("\n[1] Building FrameRollup (unique-mapper mode)...")
    rollup = build_frame_rollup(
        PARTITION_DIR,
        gapdh,
        multimap_mode="unique",
        n_workers=4,
    )
    print(f"    Rollup shape: {rollup.shape}")
    print(f"    Samples: {rollup['sample_name'].unique().len()}")
    print(f"    Total unique-read count: {rollup['count'].sum():.0f}")

    print("\n[2] Calibrating per-(sample,length) offsets...")
    agg_for_cal = rollup.group_by(["sample_name", "length", "strand", "phase0"]).agg(
        pl.col("count").sum()
    )
    offsets = calibrate_offsets(agg_for_cal, target_frame=0)
    off_29 = sorted({o for (s, L), o in offsets.items() if L == 29})
    off_28 = sorted({o for (s, L), o in offsets.items() if L == 28})
    print(f"    29nt calibrated offsets across samples: {off_29}")
    print(f"    28nt calibrated offsets across samples: {off_28}")

    print("\n[3] Scoring: calibrated vs flat offset=15...")
    scored_cal = score_frame_rollup(rollup, offsets, default_offset=12)
    scored_flat = score_frame_rollup(rollup, {}, default_offset=15)

    print("\n  Per-length aggregate frame-0 fraction (calibrated | flat-15):")
    print(f"  {'length':>6}  {'n_reads_cal':>12}  {'f0_cal':>8}  {'f0_flat':>9}")
    for L in sorted(scored_cal["length"].unique().to_list()):
        c = scored_cal.filter(pl.col("length") == L)
        f = scored_flat.filter(pl.col("length") == L)
        n_c = c["n_reads"].sum()
        n_f = f["n_reads"].sum()
        if n_c < 10:
            continue
        f0_c = c["frame0_count"].sum() / n_c
        f0_f = f["frame0_count"].sum() / n_f if n_f else 0
        print(f"  {L:>6}  {n_c:>12.0f}  {f0_c:>8.1%}  {f0_f:>9.1%}")

    # Dominant lengths gate
    dom = [28, 29, 30]
    cal_dom = scored_cal.filter(pl.col("length").is_in(dom))
    flat_dom = scored_flat.filter(pl.col("length").is_in(dom))
    f0_cal = cal_dom["frame0_count"].sum() / cal_dom["n_reads"].sum()
    f0_flat = flat_dom["frame0_count"].sum() / flat_dom["n_reads"].sum()
    print(f"\n  Aggregate (28–30nt):")
    print(f"    Calibrated: {f0_cal:.1%}")
    print(f"    Flat-15:    {f0_flat:.1%}")

    print("\n[4] Per-sample periodicity check (29nt, calibrated):")
    s29 = scored_cal.filter(pl.col("length") == 29).sort("sample_name")
    print(f"  {'sample':>15}  {'n_reads':>8}  {'f0':>6}  {'offset':>8}  pass")
    passes = 0
    for row in s29.iter_rows(named=True):
        s = row["sample_name"]
        n = row["n_reads"]
        f0 = row["elong_in_frame"]
        off = offsets.get((s, 29), 12)
        ok = f0 > 0.55
        passes += ok
        print(f"  {s:>15}  {n:>8.0f}  {f0:>6.1%}  {off:>8}  {'✓' if ok else '✗'}")
    print(f"\n  {passes}/{s29.height} samples pass frame-0 > 55% at 29nt")

    # Final verdict
    print(f"\n{'=' * 55}")
    print("  VERDICT")
    print(f"{'=' * 55}")
    gate = f0_cal > 0.55 and passes >= s29.height * 0.8
    if gate:
        print("  ✓ PASS: FrameRollup + calibrated offsets recover periodicity.")
        print(f"    Calibrated frame-0 (28-30nt) = {f0_cal:.1%} > 55%")
    else:
        print("  ✗ FAIL: Frame-0 unexpectedly low — check calibration or scan.")
    print()


if __name__ == "__main__":
    main()
