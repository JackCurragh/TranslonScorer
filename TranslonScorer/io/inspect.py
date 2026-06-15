from __future__ import annotations

import os
from typing import Iterable, Optional

import click
import polars as pl


def _print_schema(df: pl.DataFrame) -> None:
    click.echo("schema:")
    for c, t in zip(df.columns, df.dtypes):
        click.echo(f"  - {c}: {t}")


def _expand_row(r: dict) -> Iterable[dict]:
    kinds = r.get("comp_kind") or []
    fids = r.get("comp_feature_id") or []
    ss = r.get("comp_slice_start") or []
    ee = r.get("comp_slice_end") or []
    n = max(len(kinds), len(fids), len(ss), len(ee))
    for i in range(n):
        yield {
            "orf_id": r.get("orf_id"),
            "tran_id": r.get("tran_id"),
            "orf_start": r.get("start_pos_tran"),
            "orf_stop": r.get("stop_pos_tran"),
            "kind": kinds[i] if i < len(kinds) else None,
            "feature_id": fids[i] if i < len(fids) else None,
            "slice_start": ss[i] if i < len(ss) else None,
            "slice_end": ee[i] if i < len(ee) else None,
        }


def inspect_parquet(
    path: str, limit: int = 5, expand: bool = True, out_csv: Optional[str] = None
) -> None:
    df = pl.read_parquet(path)
    size = os.path.getsize(path)
    click.echo(f"file: {path}")
    click.echo(f"rows={df.height} cols={len(df.columns)} size={size:,} bytes\n")
    _print_schema(df)

    # Composite columns presence
    comp_cols = {"comp_kind", "comp_feature_id", "comp_slice_start", "comp_slice_end"}
    have_comp = comp_cols.issubset(set(df.columns))

    if not expand:
        click.echo("\npreview (head):")
        head = df.head(limit)
        # Print a subset of columns to keep it tidy
        cols = [
            c
            for c in [
                "orf_id",
                "tran_id",
                "start_pos_tran",
                "stop_pos_tran",
                "locus_id",
                "feature_chain",
                "slice_start",
                "slice_end",
            ]
            if c in head.columns
        ]
        click.echo(head.select(cols))
        return

    if have_comp:
        click.echo("\ncomposite preview:")
        shown = 0
        for r in df.head(limit).iter_rows(named=True):
            shown += 1
            click.echo(
                f"\nORF {r.get('orf_id')} | tran {r.get('tran_id')} | ORF {r.get('start_pos_tran')}→{r.get('stop_pos_tran')}"
            )
            for part in _expand_row(r):
                k = part["kind"]
                fid = part["feature_id"]
                s = part["slice_start"]
                e = part["slice_end"]
                rng = f" [{s}:{e}]" if s is not None and e is not None else ""
                click.echo(f"  - {k}: {fid}{rng}")
        if out_csv:
            rows = []
            for r in df.head(limit).iter_rows(named=True):
                rows.extend(_expand_row(r))
            pl.from_dicts(rows).write_csv(out_csv)
            click.echo(f"\nwrote expanded CSV: {out_csv}")
    else:
        click.echo("\n(no composite columns found; printing basic head)")
        click.echo(df.head(limit))
