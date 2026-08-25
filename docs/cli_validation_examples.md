# CLI validation examples

Go-forward validation method for this scoring system (established in the
significance-testing work, see `docs/significance_testing_results.md`):
run the actual shipped `translonscorer` CLI against real data, not internal
Python function calls or one-off validation scripts. Two concrete examples,
codified as a runnable artifact at `scripts/cli_smoke_tests.sh` rather than
a narrative of commands run once.

```
scripts/cli_smoke_tests.sh [OUT_DIR]
```

Requires `translonscorer` on `PATH` (`source .venv/bin/activate` from the
repo root, or `pip install -e .`). `OUT_DIR` defaults to a fresh temp
directory, removed on exit unless `CLI_SMOKE_KEEP=1`. Exit code is 0 iff
the GAPDH example passes; the pancreas example is reported but doesn't
affect the exit code when skipped (see below).

## Example 1: GAPDH BAM + GTF

Fully self-contained — repo-committed fixtures only
(`data/gapdh_cohort_genome.bam`, `data/chr12_gapdh_region.gtf`), so this
half always runs, anywhere.

```bash
# GTF chrom is "12", BAM uses "chr12" -- rename to match.
sed 's/^12\t/chr12\t/' data/chr12_gapdh_region.gtf > /tmp/gapdh_region.gtf

translonscorer extract-events \
  --out-dir out/gapdh/events \
  --gtf /tmp/gapdh_region.gtf \
  --feature-type CDS

translonscorer score-bams \
  --events-dir out/gapdh/events \
  --bam data/gapdh_cohort_genome.bam \
  --store-dir out/gapdh/scores \
  --data-version gapdh_cli_smoke \
  --offset-method global

translonscorer report \
  --store-dir out/gapdh/scores \
  --events-dir out/gapdh/events \
  --out out/gapdh/report.parquet

translonscorer consequential \
  --report out/gapdh/report.parquet \
  --out out/gapdh/report_consequential.parquet
```

**What this demonstrates**: `event_overlap` is wired into `score-bams` by
default (docs/significance_testing_results.md, "event_overlap wired into
production scoring") — no caller opt-in needed. Before that fix, every
elongation event's `competitor_share` was `{}` regardless of real overlaps
on disk.

**What to check in the output** (the script's checks, in order):

1. `out/gapdh/scores/**/*.parquet` has a non-null `confidence` value on
   more than zero rows (§6 hierarchy integration reaching the store).
2. Parsing each elongation row's `evidence` JSON, more than zero rows have
   a non-empty `competitor_share` dict — this is the actual event_overlap-
   wiring check. Last run: 184/336.
3. `out/gapdh/report_consequential.parquet` has `tier_confidence` and
   `consequential` columns and one row per translon.

**How to read a failure**: if (1)/(3) are empty, something broke upstream
of the store or report composition (unrelated to event_overlap). If (2) is
zero specifically, `event_overlap`/`comp_phase` wiring (`workflows.
_overlaps_df_for_scoring`) regressed — that's the one this example exists
to catch.

## Example 2: corrected pancreas bigwig + iRibo pancreas calls

Needs external data at `$PANCREAS_ROOT`
(`/Users/jackt/projects/all-RiboSeq/hpc_pancreas_local` by default,
override with the `PANCREAS_ROOT` env var) that is **not** checked into
this repo — only present on the machine this validation was originally run
on. The script checks for it and **skips** (not fails) this half when
absent, so the script stays runnable anywhere; only Example 1 gates the
exit code.

```bash
BED12="$PANCREAS_ROOT/pooled_endpoint/annotations/iRibo.Pancreas_pooled.bed12"
FWD="$PANCREAS_ROOT/merged_bigwigs_a15/merged_good_unique_with_junction.merged.forward.bw"
REV="$PANCREAS_ROOT/merged_bigwigs_a15/merged_good_unique_with_junction.merged.reverse.bw"

# Restrict to the chr12 GAPDH-region window, same as the GAPDH example, for
# a fast, re-runnable check rather than a genome-wide run.
awk -F'\t' '$1=="chr12" && $2<6700000 && $3>6400000' "$BED12" > /tmp/pancreas_region.bed12

translonscorer extract-events \
  --out-dir out/pancreas/events \
  --bed12 /tmp/pancreas_region.bed12 \
  --chrom chr12

translonscorer score-bigwig \
  --events-dir out/pancreas/events \
  --forward-bigwig "$FWD" \
  --reverse-bigwig "$REV" \
  --store-dir out/pancreas/scores \
  --data-version pancreas_cli_smoke
```

**Why `merged_bigwigs_a15`, not `merged_bigwigs`**: the original pancreas
follow-up found the default `merged_bigwigs/` pairing's P-site placement
disagrees with this codebase's frame convention by exactly one nucleotide,
strand-dependently (Finding 2 in `docs/significance_testing_results.md`).
`merged_bigwigs_a15/` was found (searched for, not built — see that doc's
"Bigwig correction" section) to already be correctly registered: 115/116
checked high-depth events matched `a_e` exactly, mean in-frame share 82.6%.

**What this demonstrates**: with the corrected bigwig, frame-0-specific
metrics (elongation `metric`, and by the same mechanism CIF significance)
land on the right register through the real CLI path, not just in an ad
hoc diagnostic script.

**What to check in the output**: among elongation events with >500 reads
in this window, mean `metric` should be well above the 1/3 random floor
(the uncorrected pairing gave ~0.06–0.08 despite huge depth; the corrected
pairing gives ~0.82, with 99% of high-depth events above 0.5). Last run:
125 high-depth events, mean metric 0.822.

**How to read a failure**: mean metric collapsing back toward ≤1/3 despite
plenty of high-depth events means either the bigwig pairing changed (check
the `merged_bigwigs_a15` path still resolves and hasn't been replaced) or
the frame-registration issue has reappeared for some other reason — treat
it as a real regression worth investigating, not something to route around
by swapping back to a different bigwig without understanding why.

## Extending this

Add a new example the same way: a self-contained block using either a
repo-committed fixture (preferred — keeps the whole script portable) or an
external, existence-checked path (skip cleanly when absent, as Example 2
does), followed by a `python3` block that reads the actual Parquet output
and asserts on something concrete, not just "the command exited 0".
