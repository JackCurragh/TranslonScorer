"""Debug the FrameRollup phase0 convention vs the existing rollup."""
import sys
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import polars as pl
from TranslonScorer.io.annotation import build_cds_blocks
from TranslonScorer.matrix_rollup import (
    build_frame_rollup, build_matrix_rollup, calibrate_offsets, _frame_at
)

cds = build_cds_blocks(str(REPO_ROOT / "data" / "genes.gtf"))
gapdh_cds = cds.filter(pl.col("tran_id") == "ENST00000229239")

PDIR = REPO_ROOT / "data" / "global_partitioned"

# 1. New FrameRollup
print("=== New FrameRollup (tx_pos % 3) ===")
fr = build_frame_rollup(PDIR, gapdh_cds, multimap_mode="unique", n_workers=1)
# Aggregate phase0 distribution for 29nt reads
p0_29 = (fr.filter(pl.col("length") == 29)
            .group_by("phase0").agg(pl.col("count").sum()).sort("phase0"))
print("29nt phase0 distribution:")
print(p0_29)
total_29 = fr.filter(pl.col("length") == 29)["count"].sum()
print(f"total 29nt reads: {total_29}")

# 2. Existing rollup (genomic frame intervals)
print("\n=== Existing build_matrix_rollup (exon phase → phase0) ===")
mr = build_matrix_rollup(PDIR, gapdh_cds, ref_offset=15, n_workers=1, multimap_mode="unique")
p0_29_old = (mr.filter(pl.col("length") == 29)
               .group_by("phase0").agg(pl.col("count").sum()).sort("phase0"))
print("29nt phase0 distribution (existing rollup):")
print(p0_29_old)
print(f"total 29nt reads: {mr.filter(pl.col('length') == 29)['count'].sum()}")

# 3. What does calibrate_offsets do with the new rollup?
print("\n=== calibrate_offsets on new FrameRollup ===")
# Need to aggregate across features first (new rollup has feature_id)
agg_fr = fr.group_by(["sample_name", "length", "strand", "phase0"]).agg(pl.col("count").sum())
offsets_new = calibrate_offsets(agg_fr)
off29 = {s: o for (s,L), o in offsets_new.items() if L==29}
print(f"Calibrated offsets for 29nt (new rollup): {sorted(set(off29.values()))}")

# 4. What does calibrate_offsets return for existing rollup?
print("\n=== calibrate_offsets on existing rollup ===")
offsets_old = calibrate_offsets(mr)
off29_old = {s: o for (s,L), o in offsets_old.items() if L==29}
print(f"Calibrated offsets for 29nt (existing rollup): {sorted(set(off29_old.values()))}")

# 5. Direct check: for a single sample (SRR11005875), 29nt reads, what frame is dominant?
print("\n=== Per-phase0 analysis for SRR11005875, 29nt (new rollup) ===")
s = "SRR11005875"
rows = agg_fr.filter((pl.col("sample_name") == s) & (pl.col("length") == 29))
print(rows)
print("Checking _frame_at for each candidate offset:")
data = [(int(r["phase0"]), float(r["count"])) for r in rows.iter_rows(named=True)]
for o in range(10, 19):
    fc = [0.0, 0.0, 0.0]
    for p0, c in data:
        fc[_frame_at(p0, "+", o)] += c
    total = sum(fc)
    print(f"  offset={o}: f0={fc[0]/total:.1%} f1={fc[1]/total:.1%} f2={fc[2]/total:.1%}")
