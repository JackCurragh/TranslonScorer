"""Consequentiality policy — pure transform.

Given per-translon score reports (from ``report.compose_report``) and a
``ConsequentialityPolicy``, derive a continuous consequentiality score and a
boolean label per translon.  No file I/O.

Design (per docs/event_scoring_model.md): consequentiality is *dynamic and
context-aware*, NOT a hard length/biotype gate.  A start-stop ORF or a
translated small RNA stays in; spike-concentrated pileups are down-ranked.  The
score blends:

  tier_confidence       fraction of the confidence aspects (``policy.chain_aspects``,
                        default init/elongation/term) that is SUPPORTED.
  expression_percentile rank of total translon reads across the report —
                        size-factor-normalised expression stands in here when
                        reads are already normalised upstream.

  consequentiality = tier_confidence
                     * (expression_floor + expression_weight * expression_percentile)
                     * context_weight

with the floor/weight/aspect set all carried on the policy, so a clean,
well-supported chain ranks high regardless of depth, but among equally-supported
translons the more highly translated ones rank higher.  Nothing above is fixed
in code — retune it by passing a different policy.

Public API
----------
apply_policy  — report + policy → report + consequentiality_score + consequential
"""

from __future__ import annotations

import polars as pl

from TranslonScorer.model import ConsequentialityPolicy

_SUPPORTED = "SUPPORTED"


def _confidence_call_cols(report: pl.DataFrame, policy: ConsequentialityPolicy) -> list[str]:
    """Report columns that count towards tier_confidence.

    ``policy.chain_aspects`` names the aspects; None means every aspect the
    report carries (any ``{aspect}_call`` column).  Either way, membership comes
    from the data / policy, never a list fixed in this module — a new scorer's
    aspect participates as soon as it is asked for.
    """
    if policy.chain_aspects is None:
        return [c for c in report.columns if c.endswith("_call")]
    return [f"{a}_call" for a in policy.chain_aspects if f"{a}_call" in report.columns]


def apply_policy(
    report: pl.DataFrame,
    policy: ConsequentialityPolicy = ConsequentialityPolicy(),
) -> pl.DataFrame:
    """Apply a consequentiality policy to a per-translon report DataFrame.

    Adds two columns:
      ``consequentiality_score`` (Float64 in [0, 1]) — continuous rank key.
      ``consequential`` (Boolean) — passes the policy's confidence/expression
                                    floors.

    Degrades gracefully: missing chain-call columns → tier_confidence 0;
    missing ``total_reads`` → expression treated as uninformative (percentile 1,
    never penalising) so depth never hard-gates on its own.
    """
    if report.is_empty():
        return report.with_columns(
            pl.lit(0.0).alias("consequentiality_score"),
            pl.lit(False).alias("consequential"),
        )

    # --- tier confidence: mean, across the chain aspects present, of each
    # aspect's own §6 confidence (scoring/attribution.py's compose_confidence,
    # carried through compose_report as ``{aspect}_confidence``) WHEN
    # SUPPORTED, else 0 -- not a flat mean-of-SUPPORTED. An aspect SUPPORTED
    # by a battery where every lens agreed (confidence 1.0) outranks one
    # SUPPORTED on borderline evidence (confidence 0.3); a flat
    # mean-of-SUPPORTED could not express that, since both would just count
    # as "1" once called SUPPORTED.
    #
    # Confidence defaults to 1.0 (full weight, not a penalty) wherever it is
    # unavailable -- the whole ``{aspect}_confidence`` column is missing
    # (older reports, or an aspect with no battery yet, e.g. junction), or a
    # specific translon's value is null (no lens could be evaluated for it,
    # e.g. below a codon floor) -- so the fallback is exactly the pre-§6
    # equal-weight formula, and a translon is never penalised just because
    # the new lens layer couldn't independently corroborate an otherwise-
    # SUPPORTED call (same "never penalise on missing data" convention this
    # function already uses for missing total_reads, below).
    present_calls = _confidence_call_cols(report, policy)
    if present_calls:
        contrib_exprs = []
        present_exprs = []
        for c in present_calls:
            aspect = c[: -len("_call")]
            conf_col = f"{aspect}_confidence"
            conf = pl.col(conf_col).fill_null(1.0) if conf_col in report.columns else pl.lit(1.0)
            is_present = pl.col(c).is_not_null().cast(pl.Float64)
            is_supported = (pl.col(c) == _SUPPORTED).cast(pl.Float64)
            contrib_exprs.append(is_supported * conf * is_present)
            present_exprs.append(is_present)
        n_present = sum(present_exprs[1:], present_exprs[0])
        contrib_sum = sum(contrib_exprs[1:], contrib_exprs[0])
        tier_conf = pl.when(n_present > 0).then(contrib_sum / n_present).otherwise(0.0)
    else:
        tier_conf = pl.lit(0.0)

    df = report.with_columns(tier_conf.alias("tier_confidence"))

    # --- expression percentile: rank of total_reads across the report ---
    if "total_reads" in df.columns and df.height > 1:
        df = df.with_columns(
            (pl.col("total_reads").fill_null(0.0).rank(method="average") - 1.0)
            .truediv(float(df.height - 1))
            .alias("expression_percentile")
        )
    else:
        df = df.with_columns(pl.lit(1.0).alias("expression_percentile"))

    # --- composite score + policy gate ---
    df = df.with_columns(
        (
            pl.col("tier_confidence")
            * (
                float(policy.expression_floor)
                + float(policy.expression_weight) * pl.col("expression_percentile")
            )
            * float(policy.context_weight)
        ).alias("consequentiality_score")
    )

    df = df.with_columns(
        (
            (pl.col("tier_confidence") >= float(policy.min_tier_confidence))
            & (pl.col("expression_percentile") >= float(policy.min_expression_percentile))
        ).alias("consequential")
    )

    return df
