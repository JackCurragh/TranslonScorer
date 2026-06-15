from __future__ import annotations

from typing import Dict, List, Tuple

import polars as pl

try:
    import pulp

    HAVE_PULP = True
except Exception:
    HAVE_PULP = False

from ..utils.logging import log_info
import click


def overlap_penalty(orf_i: dict, orf_j: dict) -> float:
    # Simple penalty proportional to genomic overlap length if same frame; else smaller penalty
    same_frame = int(orf_i.get("frame", -1)) == int(orf_j.get("frame", -1))
    # Use transcript positions if available; else length proxy
    Li = int(orf_i.get("length", 0))
    Lj = int(orf_j.get("length", 0))
    base = min(Li, Lj) / 3.0
    return base * (1.0 if same_frame else 0.25)


def greedy_refine(orfs_df: pl.DataFrame, time_limit: int = 10) -> pl.DataFrame:
    # Sort by score; iteratively add if penalty acceptable
    chosen: List[dict] = []
    rows_iter = orfs_df.sort("score", descending=True).iter_rows(named=True)
    iterator = click.progressbar(rows_iter, length=orfs_df.height, label="Greedy assemble")
    for r in iterator:
        penalty = sum(overlap_penalty(r, c) for c in chosen)
        if r["score"] - penalty > 0:
            chosen.append(r)
    return pl.from_dicts(chosen)


def assemble_translome(
    orfs_parquet: str, out_parquet: str, solver: str = "PULP", timeout_sec: int = 60
) -> str:
    orfs = pl.read_parquet(orfs_parquet)
    if orfs.is_empty():
        pl.DataFrame({}).write_parquet(out_parquet)
        return out_parquet

    if solver.upper() == "PULP" and HAVE_PULP:
        prob = pulp.LpProblem("translon_assembly", pulp.LpMaximize)
        xs = {
            i: pulp.LpVariable(f"x_{i}", lowBound=0, upBound=1, cat="Binary")
            for i in range(orfs.height)
        }
        # Objective: sum scores minus pairwise penalties
        score_terms = [float(orfs["score"][i]) * xs[i] for i in range(orfs.height)]
        penalty_terms = []
        # Precompute simple pairwise penalties
        for i in range(orfs.height):
            ri = {c: orfs[c][i] for c in orfs.columns}
            for j in range(i + 1, orfs.height):
                rj = {c: orfs[c][j] for c in orfs.columns}
                pen = overlap_penalty(ri, rj)
                if pen <= 0:
                    continue
                z = pulp.LpVariable(f"z_{i}_{j}", lowBound=0, upBound=1, cat="Binary")
                prob += z >= xs[i] + xs[j] - 1
                penalty_terms.append(pen * z)
        prob += pulp.lpSum(score_terms) - pulp.lpSum(penalty_terms)
        prob.solve(pulp.PULP_CBC_CMD(msg=False, timeLimit=timeout_sec))
        chosen_idx = [i for i in range(orfs.height) if xs[i].value() and xs[i].value() > 0.5]
        out = orfs.take(chosen_idx) if chosen_idx else pl.DataFrame({})
    else:
        out = greedy_refine(orfs)

    out.write_parquet(out_parquet)
    log_info(f"Assembled translome written: {out_parquet}")
    return out_parquet
