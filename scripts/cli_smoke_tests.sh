#!/usr/bin/env bash
# CLI-level example tests: the go-forward validation method for this scoring
# system (see docs/significance_testing_results.md, "CLI-level example
# tests"). Runs the ACTUAL shipped `translonscorer` command -- not internal
# Python calls -- through two concrete, real-data examples, and checks the
# output shape each one is meant to demonstrate. Re-runnable on any machine:
# the GAPDH example is fully self-contained (repo-committed fixtures); the
# pancreas example needs external data that only exists on this machine and
# is skipped, not failed, when that data isn't present.
#
# Usage:
#   scripts/cli_smoke_tests.sh [OUT_DIR]
#
# OUT_DIR defaults to a fresh temp directory (removed on exit unless
# CLI_SMOKE_KEEP=1 is set). Exit code is 0 iff the GAPDH example (always
# run) passes; the pancreas example's pass/fail is reported but does not
# affect the exit code when skipped for missing data.
#
# Requires the translonscorer package installed in the active environment
# (`pip install -e .` from the repo root, or `source .venv/bin/activate`
# if the repo's own venv already has it) and pyBigWig for the pancreas
# example specifically.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-$(mktemp -d -t translonscorer_cli_smoke)}"
mkdir -p "$OUT_DIR"
KEEP="${CLI_SMOKE_KEEP:-0}"

PANCREAS_ROOT="${PANCREAS_ROOT:-/Users/jackt/projects/all-RiboSeq/hpc_pancreas_local}"

PASS=0
FAIL=0
SKIP=0

pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }
skip() { echo "  SKIP: $1"; SKIP=$((SKIP + 1)); }

cleanup() {
  if [ "$KEEP" != "1" ]; then
    rm -rf "$OUT_DIR"
  else
    echo "Kept output at $OUT_DIR (CLI_SMOKE_KEEP=1)"
  fi
}
trap cleanup EXIT

echo "=== translonscorer CLI smoke tests ==="
echo "repo: $REPO_ROOT"
echo "out:  $OUT_DIR"
command -v translonscorer >/dev/null 2>&1 || {
  echo "FATAL: 'translonscorer' not on PATH. Activate the repo venv first:"
  echo "  source $REPO_ROOT/.venv/bin/activate"
  exit 2
}

# ---------------------------------------------------------------------------
# Example 1: GAPDH BAM + GTF (repo-committed, fully portable)
#
# Demonstrates: event_overlap is wired into score-bams by default (no
# caller opt-in) -- competitor_share is populated in the output evidence
# for any elongation event with a real overlap on disk, and the confidence
# column is populated end to end through report/consequential.
# ---------------------------------------------------------------------------
echo
echo "--- Example 1: GAPDH BAM + GTF (chr12) ---"

G="$OUT_DIR/gapdh"
mkdir -p "$G"
# GTF chrom is "12", BAM uses "chr12" -- rename to match (see
# docs/significance_testing_results.md's GAPDH validation section).
sed 's/^12\t/chr12\t/' "$REPO_ROOT/data/chr12_gapdh_region.gtf" > "$G/region.gtf"

translonscorer extract-events \
  --out-dir "$G/events" \
  --gtf "$G/region.gtf" \
  --feature-type CDS \
  > "$G/extract.log" 2>&1
if [ $? -ne 0 ]; then
  fail "extract-events (GAPDH) -- see $G/extract.log"
else
  pass "extract-events (GAPDH) ran"
fi

translonscorer score-bams \
  --events-dir "$G/events" \
  --bam "$REPO_ROOT/data/gapdh_cohort_genome.bam" \
  --store-dir "$G/scores" \
  --data-version gapdh_cli_smoke \
  --offset-method global \
  > "$G/score.log" 2>&1
if [ $? -ne 0 ]; then
  fail "score-bams (GAPDH) -- see $G/score.log"
else
  pass "score-bams (GAPDH) ran"
fi

translonscorer report \
  --store-dir "$G/scores" \
  --events-dir "$G/events" \
  --out "$G/report.parquet" \
  > "$G/report.log" 2>&1
if [ $? -ne 0 ]; then
  fail "report (GAPDH) -- see $G/report.log"
else
  pass "report (GAPDH) ran"
fi

translonscorer consequential \
  --report "$G/report.parquet" \
  --out "$G/report_consequential.parquet" \
  > "$G/consequential.log" 2>&1
if [ $? -ne 0 ]; then
  fail "consequential (GAPDH) -- see $G/consequential.log"
else
  pass "consequential (GAPDH) ran"
fi

# Check output shape: competitor_share populated (event_overlap wiring),
# confidence column present and non-trivially populated.
python3 - "$G" <<'PYEOF'
import json
import sys

import polars as pl

g = sys.argv[1]
scores = pl.read_parquet(f"{g}/scores/**/*.parquet")
report = pl.read_parquet(f"{g}/report_consequential.parquet")

ok = True

if "confidence" not in scores.columns:
    print("  FAIL: scores parquet has no confidence column")
    ok = False
else:
    n_conf = scores.filter(pl.col("confidence").is_not_null()).height
    print(f"  {'PASS' if n_conf > 0 else 'FAIL'}: {n_conf}/{scores.height} scored events have a confidence value")
    ok &= n_conf > 0

elong = scores.filter(pl.col("aspect") == "elongation")
n_with_competitor = 0
for row in elong.iter_rows(named=True):
    ev = json.loads(row["evidence"])
    if ev.get("competitor_share"):
        n_with_competitor += 1
print(f"  {'PASS' if n_with_competitor > 0 else 'FAIL'}: {n_with_competitor}/{elong.height} elongation events have a non-empty competitor_share (event_overlap wiring)")
ok &= n_with_competitor > 0

if "tier_confidence" not in report.columns or "consequential" not in report.columns:
    print("  FAIL: report parquet missing tier_confidence/consequential columns")
    ok = False
else:
    print(f"  PASS: report has {report.height} translon row(s) with tier_confidence/consequential")

sys.exit(0 if ok else 1)
PYEOF
if [ $? -eq 0 ]; then
  pass "GAPDH output shape checks"
else
  fail "GAPDH output shape checks"
fi

# ---------------------------------------------------------------------------
# Example 2: corrected pancreas bigwig + iRibo pancreas ORF calls
#
# Demonstrates: frame-specific metrics (elongation `metric`, CIF
# significance) land on the correct register once paired with the
# frame-corrected bigwig (merged_bigwigs_a15), unlike the uncorrected
# merged_bigwigs/ pairing used in the original pancreas follow-up (see
# "Finding 2" in docs/significance_testing_results.md).
#
# Needs external data at $PANCREAS_ROOT that is NOT checked into this repo
# -- skipped (not failed) when absent, so this script stays runnable
# anywhere.
# ---------------------------------------------------------------------------
echo
echo "--- Example 2: corrected pancreas bigwig (merged_bigwigs_a15) + iRibo calls ---"

BED12="$PANCREAS_ROOT/pooled_endpoint/annotations/iRibo.Pancreas_pooled.bed12"
FWD="$PANCREAS_ROOT/merged_bigwigs_a15/merged_good_unique_with_junction.merged.forward.bw"
REV="$PANCREAS_ROOT/merged_bigwigs_a15/merged_good_unique_with_junction.merged.reverse.bw"

if [ ! -f "$BED12" ] || [ ! -f "$FWD" ] || [ ! -f "$REV" ]; then
  skip "pancreas data not found under $PANCREAS_ROOT (expected on jackt's machine only)"
else
  P="$OUT_DIR/pancreas"
  mkdir -p "$P"
  # Restrict to the same chr12 GAPDH-region window as the GAPDH example, for
  # a fast, re-runnable smoke check rather than a genome-wide run.
  awk -F'\t' '$1=="chr12" && $2<6700000 && $3>6400000' "$BED12" > "$P/region.bed12"

  translonscorer extract-events \
    --out-dir "$P/events" \
    --bed12 "$P/region.bed12" \
    --chrom chr12 \
    > "$P/extract.log" 2>&1
  if [ $? -ne 0 ]; then
    fail "extract-events (pancreas) -- see $P/extract.log"
  else
    pass "extract-events (pancreas) ran"
  fi

  translonscorer score-bigwig \
    --events-dir "$P/events" \
    --forward-bigwig "$FWD" \
    --reverse-bigwig "$REV" \
    --store-dir "$P/scores" \
    --data-version pancreas_cli_smoke \
    > "$P/score.log" 2>&1
  if [ $? -ne 0 ]; then
    fail "score-bigwig (pancreas) -- see $P/score.log"
  else
    pass "score-bigwig (pancreas) ran"
  fi

  python3 - "$P" <<'PYEOF'
import sys

import polars as pl

p = sys.argv[1]
scores = pl.read_parquet(f"{p}/scores/**/*.parquet")
elong = scores.filter((pl.col("aspect") == "elongation") & (pl.col("n_reads") > 0))

ok = True
if elong.is_empty():
    print("  FAIL: no elongation events with reads")
    ok = False
else:
    # The uncorrected pairing (see Finding 2) put frame-0-specific metrics
    # near/below the 1/3 random floor almost everywhere despite huge depth
    # (elongation SUPPORTED was 1/180). Post-correction, high-depth events
    # should show metric well above 1/3 -- not a hard threshold (real
    # biology varies), just "clearly not still measuring the wrong frame".
    high_depth = elong.filter(pl.col("n_reads") > 500)
    if high_depth.is_empty():
        print("  SKIP: no high-depth (>500 reads) elongation events in this window")
    else:
        mean_metric = high_depth["metric"].mean()
        frac_above_half = (high_depth["metric"] > 0.5).sum() / high_depth.height
        print(f"  high-depth elongation events: {high_depth.height}, mean metric={mean_metric:.3f}, "
              f"fraction with metric>0.5: {frac_above_half:.2f}")
        # Uncorrected run: mean metric was ~0.06 (below the 1/3 random
        # floor). Corrected: should be well above 1/3 on a housekeeping-
        # gene-dense window.
        print(f"  {'PASS' if mean_metric > 0.5 else 'FAIL'}: mean elongation metric ({mean_metric:.3f}) is well above the 1/3 random floor")
        ok &= mean_metric > 0.5

sys.exit(0 if ok else 1)
PYEOF
  if [ $? -eq 0 ]; then
    pass "pancreas output shape checks"
  else
    fail "pancreas output shape checks"
  fi
fi

# ---------------------------------------------------------------------------
echo
echo "=== Summary: $PASS passed, $FAIL failed, $SKIP skipped ==="
if [ "$FAIL" -gt 0 ]; then
  exit 1
fi
exit 0
