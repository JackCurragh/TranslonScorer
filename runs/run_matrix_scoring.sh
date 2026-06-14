#!/usr/bin/env bash
#
# Reproducible event-scoring run on the sparse annotation-scale matrix.
#
# Chains the four new TranslonScorer subcommands end to end:
#   extract-events  annotation sqlite        -> genomic event store
#   score-matrix    events + matrix          -> append-only score store
#   report          scores + events          -> per-translon report + policy
#   consequential   report (strict re-gate)  -> consequential subset
#
# Everything below the CONFIG block is plain CLI calls — no bespoke scripts.
# Edit CONFIG and re-run; outputs land under $OUTDIR, which is git-ignored.
#
# Usage:
#   ./run_matrix_scoring.sh                # defaults (chr12, 16 partitions)
#   CHROM=chr7 N_PARTITIONS=0 ./run_matrix_scoring.sh   # chr7, full cohort
#
set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG (override via environment)
# ---------------------------------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"                      # translonscorer/
ALL_RIBOSEQ="$(cd "$REPO/.." && pwd)"              # all-RiboSeq/

SQLITE="${SQLITE:-$ALL_RIBOSEQ/translon_db_rebuild/translon_db/translons.sqlite}"
MATRIX_DIR="${MATRIX_DIR:-$REPO/data/global_partitioned}"
CHROM="${CHROM:-chr12}"                             # chromosome to score
DATA_VERSION="${DATA_VERSION:-matrix_demo}"         # score-store partition label
ANNOTATION_VERSION="${ANNOTATION_VERSION:-demo}"    # stamped into event ids
# Number of matrix partitions to scan. A subset = a (lower-depth) read sample of
# the cohort and keeps the demo fast; set to 0 to use ALL partitions (full depth).
N_PARTITIONS="${N_PARTITIONS:-16}"
# Strict consequential floor for the re-gate step (1.0 = full chain supported).
STRICT_TIER_CONFIDENCE="${STRICT_TIER_CONFIDENCE:-1.0}"

OUTDIR="${OUTDIR:-$HERE/matrix_${CHROM}}"

# ---------------------------------------------------------------------------
# Resolve the CLI (installed entry point, else module form)
# ---------------------------------------------------------------------------
if command -v translonscorer >/dev/null 2>&1; then
  TS=(translonscorer)
else
  TS=(python3 -m TranslonScorer.cli)
fi

# ---------------------------------------------------------------------------
# Build the --partitions argument list
# ---------------------------------------------------------------------------
ALL_PARTS=()
while IFS= read -r d; do ALL_PARTS+=("$d"); done < <(find "$MATRIX_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
if [[ "${#ALL_PARTS[@]}" -eq 0 ]]; then
  echo "ERROR: no partition directories under $MATRIX_DIR" >&2
  exit 1
fi
if [[ "$N_PARTITIONS" -gt 0 ]]; then
  PARTS=("${ALL_PARTS[@]:0:$N_PARTITIONS}")
else
  PARTS=("${ALL_PARTS[@]}")
fi
PART_ARGS=()
for p in "${PARTS[@]}"; do PART_ARGS+=(--partitions "$p"); done

echo "=================================================================="
echo " TranslonScorer matrix event-scoring run"
echo "   sqlite      : $SQLITE"
echo "   matrix      : $MATRIX_DIR (${#PARTS[@]} of ${#ALL_PARTS[@]} partitions)"
echo "   chromosome  : $CHROM"
echo "   data version: $DATA_VERSION"
echo "   output      : $OUTDIR"
echo "=================================================================="

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# 1. extract-events
# ---------------------------------------------------------------------------
echo ">> [1/4] extract-events ($CHROM)"
"${TS[@]}" extract-events \
  --sqlite "$SQLITE" \
  --out-dir "$OUTDIR/events" \
  --chrom "$CHROM" \
  --annotation-version "$ANNOTATION_VERSION"

# ---------------------------------------------------------------------------
# 2. score-matrix
# ---------------------------------------------------------------------------
echo ">> [2/4] score-matrix"
"${TS[@]}" score-matrix \
  --events-dir "$OUTDIR/events" \
  "${PART_ARGS[@]}" \
  --store-dir "$OUTDIR/scores" \
  --data-version "$DATA_VERSION" \
  --annotation-version "$ANNOTATION_VERSION"

# ---------------------------------------------------------------------------
# 3. report (compose per-translon + apply default policy)
# ---------------------------------------------------------------------------
echo ">> [3/4] report"
"${TS[@]}" report \
  --store-dir "$OUTDIR/scores" \
  --events-dir "$OUTDIR/events" \
  --out "$OUTDIR/report.parquet" \
  --data-version "$DATA_VERSION"

# ---------------------------------------------------------------------------
# 4. consequential (strict re-gate: full translation chain supported)
# ---------------------------------------------------------------------------
echo ">> [4/4] consequential (min-tier-confidence=$STRICT_TIER_CONFIDENCE)"
"${TS[@]}" consequential \
  --report "$OUTDIR/report.parquet" \
  --out "$OUTDIR/consequential_strict.parquet" \
  --min-tier-confidence "$STRICT_TIER_CONFIDENCE"

echo "=================================================================="
echo " DONE. Outputs in $OUTDIR:"
echo "   events/                  genomic event store (events, feature_event)"
echo "   scores/                  append-only fact_event_score store"
echo "   report.parquet           per-translon report + consequentiality_score"
echo "   consequential_strict.parquet   full-chain-supported subset"
echo "=================================================================="
