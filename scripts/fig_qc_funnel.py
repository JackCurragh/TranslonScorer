"""Talk figure (slide 10) + cohort-number reconciliation tool.

Reads the baked psite_index outputs and:
  1. PRINTS the schema and headline counts for each gate — this resolves the
     2,270 vs 1,387 discrepancy (PIF-usable vs scoring-ready ≥10k reads).
  2. Renders a two-panel figure:
       Left  — the QC funnel: total -> has CDS reads -> QC/PIF-usable -> scoring-ready.
       Right — failure-mode scatter: prop_cds (x) vs in-frame periodicity (y),
               coloured by trimming quality (rpf_28_32_prop). Shows the two
               failure modes separating: low x = rRNA contamination,
               low y = aperiodic/bad processing.

Run on HPC where the parquets live:
  python3 scripts/fig_qc_funnel.py \
      --psite-index /hps/.../translonscorer_6k/psite_index_filtered \
      --out outputs/fig_qc_funnel.png

If a column name differs from the defaults, the script prints every available
column so you can pass the right one via the --*-col flags.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import polars as pl  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--psite-index", required=True, type=Path,
                   help="psite_index_filtered directory")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--min-usable-reads", type=int, default=10_000,
                   help="scoring-ready gate: total usable reads per sample")
    # column-name overrides (defaults match the documented schema)
    p.add_argument("--sample-col", default=None,
                   help="sample id column (auto-detected: sample_id / sample_name)")
    p.add_argument("--propcds-col", default="prop_cds")
    p.add_argument("--period-col", default=None,
                   help="periodicity column (auto: periodicity_score / f0 / f0_pct)")
    p.add_argument("--trim-col", default="rpf_28_32_prop")
    p.add_argument("--reads-col", default=None,
                   help="usable-reads column in usable_sample_lengths (auto: n_reads / count)")
    return p.parse_args()


def pick(df: pl.DataFrame, override: str | None, candidates: list[str], what: str) -> str:
    if override:
        if override not in df.columns:
            raise SystemExit(f"--{what} column '{override}' not in {df.columns}")
        return override
    for c in candidates:
        if c in df.columns:
            return c
    raise SystemExit(f"Could not auto-detect {what}; available columns: {df.columns}")


def load(base: Path, name: str) -> pl.DataFrame | None:
    f = base / name
    if not f.exists():
        print(f"  MISSING: {name}")
        return None
    df = pl.read_parquet(f)
    print(f"  {name}: shape={df.shape}")
    print(f"      columns: {df.columns}")
    return df


def main() -> None:
    args = parse_args()
    base = args.psite_index

    print(f"Reading psite_index outputs from {base}\n")
    qc = load(base, "qc_per_sample.parquet")
    pif = load(base, "usable_samples_pif.parquet")
    ul = load(base, "usable_sample_lengths.parquet")
    _ = load(base, "per_sample_length_frames.parquet")

    if qc is None:
        raise SystemExit("qc_per_sample.parquet is required.")

    sample_col = pick(qc, args.sample_col, ["sample_id", "sample_name"], "sample-col")
    propcds_col = pick(qc, args.propcds_col, ["prop_cds"], "propcds-col")
    period_col = pick(qc, args.period_col,
                      ["periodicity_score", "f0", "f0_pct", "elong_in_frame"], "period-col")

    # ---- funnel counts (this is the reconciliation) ----
    total = qc[sample_col].n_unique()
    has_cds = qc.filter(pl.col(propcds_col) > 0)[sample_col].n_unique()
    pif_usable = pif[pick(pif, args.sample_col, ["sample_id", "sample_name"], "sample-col")
                     ].n_unique() if pif is not None else None

    scoring_ready = None
    if ul is not None:
        ul_sample = pick(ul, args.sample_col, ["sample_id", "sample_name"], "sample-col")
        reads_col = pick(ul, args.reads_col, ["n_reads", "count", "usable_reads"], "reads-col")
        per_sample = ul.group_by(ul_sample).agg(pl.col(reads_col).sum().alias("ur"))
        scoring_ready = per_sample.filter(pl.col("ur") >= args.min_usable_reads).height

    print("\n================ COHORT FUNNEL ================")
    print(f"  total samples            : {total}")
    print(f"  with CDS reads           : {has_cds}")
    print(f"  PIF-usable (QC pass)      : {pif_usable}")
    print(f"  scoring-ready (>={args.min_usable_reads:,} reads): {scoring_ready}")
    print("===============================================")
    print("  ^ 2,270 should be PIF-usable; 1,387 (if that is what appears\n"
          "    on the last line) is the stricter scoring-ready gate.\n")

    # ---- plot ----
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                                   gridspec_kw={"width_ratios": [1.1, 1.4]})

    stages, counts = ["total"], [total]
    stages.append("has CDS reads"); counts.append(has_cds)
    if pif_usable is not None:
        stages.append("QC / PIF usable"); counts.append(pif_usable)
    if scoring_ready is not None:
        stages.append(f"scoring-ready\n(≥{args.min_usable_reads // 1000}k)")
        counts.append(scoring_ready)

    colors = ["#cccccc", "#9ecae1", "#4292c6", "#08519c"][: len(stages)]
    bars = axL.barh(range(len(stages)), counts, color=colors)
    axL.set_yticks(range(len(stages)))
    axL.set_yticklabels(stages)
    axL.invert_yaxis()
    axL.set_xlabel("samples")
    axL.set_title("Cohort triage from the matrix")
    for b, c in zip(bars, counts):
        axL.text(b.get_width(), b.get_y() + b.get_height() / 2, f" {c:,}",
                 va="center", ha="left", fontsize=10, fontweight="bold")

    # scatter: prop_cds vs periodicity, coloured by trimming quality
    sdf = qc.select(
        [c for c in {propcds_col, period_col, args.trim_col} if c in qc.columns]
    ).drop_nulls()
    y = sdf[period_col].to_numpy()
    # normalise periodicity to 0-1 if it is a percentage
    if y.max() > 1.5:
        y = y / 100.0
    sc = axR.scatter(
        sdf[propcds_col].to_numpy(), y,
        c=(sdf[args.trim_col].to_numpy() if args.trim_col in sdf.columns else "#4292c6"),
        cmap="viridis", s=8, alpha=0.5,
    )
    axR.axhline(1 / 3, ls="--", lw=1, color="#999999")
    axR.set_xlabel("prop_cds  (low → rRNA / library contamination)")
    axR.set_ylabel("in-frame fraction  (low → aperiodic / bad processing)")
    axR.set_title("Failure modes separate")
    if args.trim_col in sdf.columns:
        cb = fig.colorbar(sc, ax=axR)
        cb.set_label("rpf_28_32_prop (trim quality)", fontsize=8)

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
