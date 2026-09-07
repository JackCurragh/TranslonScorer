from __future__ import annotations

from typing import Any

import polars as pl


def _fmean(series: pl.Series) -> float:
    """`series.mean()` as a plain float.

    Polars types `Series.mean()` as a union spanning int/float/Decimal/date/
    time/timedelta/str/bytes/ndarray/list/None, because a Series can hold any
    of those. Every arithmetic or `float()` use of it therefore fails type
    checking even on columns that are always numeric here. Narrowing once, in
    one place, is honest about that being a stub limitation rather than
    scattering per-call ignores through the module.

    Callers must guard emptiness themselves (`height`/`is_empty`); on an empty
    Series polars returns None and this raises, which is the right failure.
    """
    return float(series.mean())  # type: ignore[arg-type]


def _fmedian(series: pl.Series) -> float:
    """`series.median()` as a plain float. See `_fmean`."""
    return float(series.median())  # type: ignore[arg-type]


RAW_SCORE_COLUMNS = ["rise_up", "step_down", "hrf", "avg", "nzc", "score"]

FRAME_SCORE_COLUMNS = [
    "frame_posterior_mean",
    "frame_entropy_mean",
    "frame_weighted_count",
    "frame_support_codons",
    "frame_method",
]

AMBIGUITY_COLUMNS = [
    "assignment_weight",
    "assignment_entropy",
    "assignment_model",
    "identifiability_class",
]

PROVENANCE_COLUMNS = [
    "score_mode",
    "input_type",
    "offset_source",
    "annotation_version",
    "method_version",
]

SCORE_SCHEMA_COLUMNS = [
    "orf_id",
    "tran_id",
    "gene_id",
    "start",
    "stop",
    "type",
    *RAW_SCORE_COLUMNS,
    *FRAME_SCORE_COLUMNS,
    *AMBIGUITY_COLUMNS,
    *PROVENANCE_COLUMNS,
]


def _start_col(df: pl.DataFrame) -> str | None:
    for name in ("start", "start_pos_tran", "orf_start", "tran_start"):
        if name in df.columns:
            return name
    return None


def _stop_col(df: pl.DataFrame) -> str | None:
    for name in ("stop", "stop_pos_tran", "orf_stop", "tran_stop"):
        if name in df.columns:
            return name
    return None


def _load_table(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def load_score_table(path: str) -> pl.DataFrame:
    return _load_table(path)


def ensure_orf_id(df: pl.DataFrame) -> pl.DataFrame:
    """Ensure a stable ORF key exists for joining raw/frame score tables."""
    if "orf_id" in df.columns:
        return df

    start = _start_col(df)
    stop = _stop_col(df)
    if "tran_id" in df.columns and start and stop:
        parts = [
            pl.col("tran_id").cast(pl.Utf8),
            pl.col(start).cast(pl.Utf8),
            pl.col(stop).cast(pl.Utf8),
        ]
        if "type" in df.columns:
            parts.append(pl.col("type").cast(pl.Utf8))
        return df.with_columns(pl.concat_str(parts, separator=":").alias("orf_id"))

    return (
        df.with_row_index("_score_row")
        .with_columns(
            pl.concat_str([pl.lit("row"), pl.col("_score_row").cast(pl.Utf8)], separator=":").alias(
                "orf_id"
            )
        )
        .drop("_score_row")
    )


def ensure_score_schema(
    df: pl.DataFrame,
    *,
    score_mode: str = "raw",
    input_type: str | None = None,
    frame_method: str | None = None,
    assignment_model: str = "none",
    identifiability_class: str = "not_evaluated",
    offset_source: str | None = None,
    annotation_version: str | None = None,
    method_version: str | None = None,
) -> pl.DataFrame:
    """Add the score-first schema columns without discarding existing fields."""
    df = ensure_orf_id(df)

    additions: list[pl.Expr] = []

    if "gene_id" not in df.columns:
        additions.append(pl.lit(None).cast(pl.Utf8).alias("gene_id"))

    for col in RAW_SCORE_COLUMNS:
        if col not in df.columns:
            additions.append(pl.lit(None).cast(pl.Float64).alias(col))

    frame_defaults: dict[str, Any] = {
        "frame_posterior_mean": None,
        "frame_entropy_mean": None,
        "frame_weighted_count": None,
        "frame_support_codons": None,
        "frame_method": frame_method,
    }
    for col, value in frame_defaults.items():
        if col not in df.columns:
            dtype = (
                pl.Utf8
                if col == "frame_method"
                else (pl.Int64 if col == "frame_support_codons" else pl.Float64)
            )
            additions.append(pl.lit(value).cast(dtype).alias(col))

    ambiguity_defaults: dict[str, Any] = {
        "assignment_weight": 1.0,
        "assignment_entropy": 0.0,
        "assignment_model": assignment_model,
        "identifiability_class": identifiability_class,
    }
    for col, value in ambiguity_defaults.items():
        if col not in df.columns:
            dtype = pl.Utf8 if col in {"assignment_model", "identifiability_class"} else pl.Float64
            additions.append(pl.lit(value).cast(dtype).alias(col))

    provenance_defaults: dict[str, Any] = {
        "score_mode": score_mode,
        "input_type": input_type,
        "offset_source": offset_source,
        "annotation_version": annotation_version,
        "method_version": method_version,
    }
    for col, value in provenance_defaults.items():
        if col not in df.columns:
            additions.append(pl.lit(value).cast(pl.Utf8).alias(col))

    if additions:
        df = df.with_columns(additions)

    return df


def _collapse_frame_support(frame_support: pl.DataFrame) -> pl.DataFrame:
    if frame_support.is_empty():
        return frame_support
    cols = set(frame_support.columns)
    required = {"tran_id", "codon", "p0", "p1", "p2"}
    if not required.issubset(cols):
        raise ValueError(f"Frame support missing required columns: {sorted(required - cols)}")

    agg_exprs = [
        pl.col("p0").mean().alias("p0"),
        pl.col("p1").mean().alias("p1"),
        pl.col("p2").mean().alias("p2"),
    ]
    if "entropy" in cols:
        agg_exprs.append(pl.col("entropy").mean().alias("entropy"))
    else:
        agg_exprs.append(pl.lit(None).cast(pl.Float64).alias("entropy"))
    if "method" in cols:
        agg_exprs.append(pl.col("method").first().alias("method"))
    else:
        agg_exprs.append(pl.lit(None).cast(pl.Utf8).alias("method"))

    return frame_support.group_by(["tran_id", "codon"]).agg(agg_exprs).sort(["tran_id", "codon"])


def add_frame_score_columns(
    scored: pl.DataFrame,
    frame_support: pl.DataFrame,
    profiles: pl.DataFrame | None = None,
) -> pl.DataFrame:
    """Attach per-ORF frame posterior summaries to a scored ORF table.

    `frame_posterior_mean` is the mean posterior for the ORF's start-frame lane
    across covered codons. `frame_weighted_count` is computed when transcript
    profiles are supplied, using count * posterior for each nucleotide's frame.
    """
    scored = ensure_score_schema(scored)
    if scored.is_empty() or frame_support.is_empty():
        return scored

    start = _start_col(scored)
    stop = _stop_col(scored)
    if not start or not stop or "tran_id" not in scored.columns:
        return scored

    fs = _collapse_frame_support(frame_support)
    fs_by_tx = {
        str(key[0] if isinstance(key, tuple) else key): sub for key, sub in fs.group_by("tran_id")
    }

    prof_by_tx: dict[str, pl.DataFrame] = {}
    if (
        profiles is not None
        and not profiles.is_empty()
        and {"tran_id", "pos", "count"}.issubset(set(profiles.columns))
    ):
        prof_norm = (
            profiles.select(["tran_id", "pos", "count"])
            .group_by(["tran_id", "pos"])
            .agg(pl.col("count").sum().alias("count"))
            .sort(["tran_id", "pos"])
        )
        prof_by_tx = {
            str(key[0] if isinstance(key, tuple) else key): sub
            for key, sub in prof_norm.group_by("tran_id")
        }

    rows: list[dict[str, Any]] = []
    for row in scored.iter_rows(named=True):
        tid = str(row.get("tran_id"))
        tx_fs = fs_by_tx.get(tid)
        start_val = row.get(start)
        stop_val = row.get(stop)
        if tx_fs is None or start_val is None or stop_val is None:
            rows.append(
                {
                    "orf_id": row["orf_id"],
                    "frame_posterior_mean": None,
                    "frame_entropy_mean": None,
                    "frame_weighted_count": None,
                    "frame_support_codons": 0,
                    "frame_method": row.get("frame_method"),
                }
            )
            continue

        s = int(start_val)
        e = int(stop_val)
        if e <= s:
            rows.append(
                {
                    "orf_id": row["orf_id"],
                    "frame_posterior_mean": None,
                    "frame_entropy_mean": None,
                    "frame_weighted_count": None,
                    "frame_support_codons": 0,
                    "frame_method": row.get("frame_method"),
                }
            )
            continue

        frame = row.get("frame")
        if frame is None:
            frame = s % 3
        frame = int(frame) % 3
        p_col = f"p{frame}"
        codon_lo = s // 3
        codon_hi = ((e - 1) // 3) + 1
        sub = tx_fs.filter((pl.col("codon") >= codon_lo) & (pl.col("codon") < codon_hi))

        posterior = _fmean(sub[p_col]) if not sub.is_empty() else None
        entropy = (
            _fmean(sub["entropy"]) if not sub.is_empty() and "entropy" in sub.columns else None
        )
        method = (
            str(sub["method"][0])
            if not sub.is_empty() and "method" in sub.columns and sub["method"][0] is not None
            else row.get("frame_method")
        )

        weighted_count = None
        tx_prof = prof_by_tx.get(tid)
        if tx_prof is not None and not sub.is_empty():
            prof_sub = (
                tx_prof.filter((pl.col("pos") >= s) & (pl.col("pos") < e))
                .with_columns(
                    (pl.col("pos") // 3).alias("codon"), (pl.col("pos") % 3).alias("pos_frame")
                )
                .join(sub.select(["codon", "p0", "p1", "p2"]), on="codon", how="left")
            )
            total = 0.0
            for p_row in prof_sub.iter_rows(named=True):
                pf = int(p_row["pos_frame"])
                p_val = p_row.get(f"p{pf}")
                if p_val is not None:
                    total += float(p_row["count"]) * float(p_val)
            weighted_count = total

        rows.append(
            {
                "orf_id": row["orf_id"],
                "frame_posterior_mean": posterior,
                "frame_entropy_mean": entropy,
                "frame_weighted_count": weighted_count,
                "frame_support_codons": sub.height,
                "frame_method": method,
            }
        )

    summary = pl.from_dicts(rows)
    drop_cols = [c for c in FRAME_SCORE_COLUMNS if c in scored.columns]
    return scored.drop(drop_cols).join(summary, on="orf_id", how="left")


def compare_score_tables(
    raw: pl.DataFrame, frame: pl.DataFrame, *, label_column: str | None = None
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Compare raw and frame-weighted score tables by ORF key."""
    raw_s = ensure_score_schema(raw, score_mode="raw")
    frame_s = ensure_score_schema(frame, score_mode="frame_weighted")

    keep = ["orf_id", "score", "hrf", "avg", "nzc", "rise_up", "step_down"]
    if label_column and label_column in raw_s.columns:
        keep.append(label_column)

    joined = raw_s.select([c for c in keep if c in raw_s.columns]).join(
        frame_s.select([c for c in keep if c in frame_s.columns]),
        on="orf_id",
        how="inner",
        suffix="_frame",
    )

    for col in ("score", "hrf", "avg", "nzc", "rise_up", "step_down"):
        frame_col = f"{col}_frame"
        if col in joined.columns and frame_col in joined.columns:
            joined = joined.with_columns((pl.col(frame_col) - pl.col(col)).alias(f"{col}_delta"))

    summary_values: dict[str, Any] = {
        "n_raw": raw_s.height,
        "n_frame": frame_s.height,
        "n_joined": joined.height,
        "mean_score_raw": (
            _fmean(raw_s["score"]) if "score" in raw_s.columns and raw_s.height else None
        ),
        "mean_score_frame": (
            _fmean(frame_s["score"]) if "score" in frame_s.columns and frame_s.height else None
        ),
        "median_score_delta": (
            _fmedian(joined["score_delta"])
            if "score_delta" in joined.columns and joined.height
            else None
        ),
        "mean_score_delta": (
            _fmean(joined["score_delta"])
            if "score_delta" in joined.columns and joined.height
            else None
        ),
    }

    if label_column and label_column in joined.columns:
        labels = joined.get_column(label_column).drop_nulls().unique().to_list()
        if len(labels) == 2:
            label_a, label_b = labels[0], labels[1]
            a = joined.filter(pl.col(label_column) == label_a)
            b = joined.filter(pl.col(label_column) == label_b)
            summary_values["label_a"] = str(label_a)
            summary_values["label_b"] = str(label_b)
            summary_values["raw_label_score_gap"] = (
                (_fmean(a["score"]) - _fmean(b["score"])) if a.height and b.height else None
            )
            summary_values["frame_label_score_gap"] = (
                (_fmean(a["score_frame"]) - _fmean(b["score_frame"]))
                if a.height and b.height
                else None
            )

    return joined, pl.DataFrame([summary_values])
