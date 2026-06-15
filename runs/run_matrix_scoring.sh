#!/usr/bin/env bash
#
# Reproducible event-scoring run on the sparse annotation-scale matrix.
#
# The matrix is sharded by read sequence; ALL partitions under MATRIX_DIR are
# always scanned together (there is no valid single-partition / subset usage).
# Point --matrix-dir at the root and the CLI discovers every partition.
#
# Chains the four subcommands end to end (or use `translonscorer pipeline`):
#   extract-events  annotation sqlite  -> genomic event store
#   score-matrix    events + matrix    -> append-only score store
#   report          scores + events    -> per-translon report + policy
#   consequential   report (re-gate)   -> strict consequential subset
#
# Edit CONFIG and re-run; outputs land under $OUTDIR (git-ignored).
#
# Usage:
#   ./run_matrix_scoring.sh                 # defaults (chr12)
#   CHROM= ./run_matrix_scoring.sh          # whole genome (no chrom restriction)
#
set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG (override via environment)
# ---------------------------------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"                      # translonscorer/
ALL_RIBOSEQ="$(cd "$REPO/.." && pwd)"              # all-RiboSeq/

SQLITE="${SQLITE:-$ALL_RIBOSEQ/translon_db_rebuild/translon_db/translons.sqlite}"
MATRIX_DIR="${MATRIX_DIR:-$REPO/data/global_partitioned}"   # matrix ROOT (all partitions)
CHROM="${CHROM-chr12}"                              # extraction scope; empty = whole genome
DATA_VERSION="${DATA_VERSION:-matrix_demo}"         # score-store partition label
ANNOTATION_VERSION="${ANNOTATION_VERSION:-demo}"    # stamped into event ids
STRICT_TIER_CONFIDENCE="${STRICT_TIER_CONFIDENCE:-1.0}"   # 1.0 = full chain supported

OUTDIR="${OUTDIR:-$HERE/matrix_${CHROM:-genome}}"

# Resolve the CLI (installed entry point, else module form).
if command -v translonscorer >/dev/null 2>&1; then
  TS=(translonscorer)
else
  TS=(python3 -m TranslonScorer.cli)
fi

CHROM_ARGS=()
[[ -n "$CHROM" ]] && CHROM_ARGS=(--chrom "$CHROM")

echo "=================================================================="
echo " TranslonScorer matrix event-scoring run"
echo "   sqlite      : $SQLITE"
echo "   matrix root : $MATRIX_DIR (all partitions)"
echo "   chromosome  : ${CHROM:-<whole genome>}"
echo "   output      : $OUTDIR"
echo "=================================================================="

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo ">> [1/4] extract-events"
"${TS[@]}" extract-events \
  --sqlite "$SQLITE" \
  --out-dir "$OUTDIR/events" \
  "${CHROM_ARGS[@]}" \
  --annotation-version "$ANNOTATION_VERSION"

echo ">> [2/4] score-matrix (all partitions under matrix root)"
"${TS[@]}" score-matrix \
  --events-dir "$OUTDIR/events" \
  --matrix-dir "$MATRIX_DIR" \
  --store-dir "$OUTDIR/scores" \
  --data-version "$DATA_VERSION" \
  --annotation-version "$ANNOTATION_VERSION"

echo ">> [3/4] report"
"${TS[@]}" report \
  --store-dir "$OUTDIR/scores" \
  --events-dir "$OUTDIR/events" \
  --out "$OUTDIR/report.parquet" \
  --data-version "$DATA_VERSION"

echo ">> [4/4] consequential (min-tier-confidence=$STRICT_TIER_CONFIDENCE)"
"${TS[@]}" consequential \
  --report "$OUTDIR/report.parquet" \
  --out "$OUTDIR/consequential_strict.parquet" \
  --min-tier-confidence "$STRICT_TIER_CONFIDENCE"

echo "=================================================================="
echo " DONE. Outputs in $OUTDIR:"
echo "   events/                        genomic event store"
echo "   scores/                        append-only fact_event_score store"
echo "   report.parquet                 per-translon report + consequentiality"
echo "   consequential_strict.parquet   full-chain-supported subset"
echo "=================================================================="
