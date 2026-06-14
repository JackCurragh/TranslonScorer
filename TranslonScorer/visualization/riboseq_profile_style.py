from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from matplotlib.patches import Patch


FRAME_COLORS = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#54A24B",
}


def phase_colors(positions: Sequence[int]) -> list[str]:
    return [FRAME_COLORS[int(pos) % 3] for pos in positions]


def plot_phase_bars(
    ax: Any,
    positions: Sequence[int],
    counts: Sequence[float],
    *,
    width: float = 0.9,
    alpha: float = 0.95,
) -> None:
    """Plot a standard Ribo-seq profile as phase-coloured vertical bars."""
    if len(positions) == 0:
        return
    pos = np.asarray(positions, dtype=float)
    vals = np.asarray(counts, dtype=float)
    ax.bar(
        pos,
        vals,
        width=width,
        color=phase_colors([int(p) for p in pos]),
        edgecolor="none",
        alpha=alpha,
        align="center",
    )


def style_profile_axis(
    ax: Any,
    *,
    xlabel: str | None = None,
    ylabel: str | None = None,
    show_y_grid: bool = True,
) -> None:
    """Apply a clean track-like style for transcript-position Ribo-seq profiles."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#333333")
    ax.spines["bottom"].set_color("#333333")
    ax.tick_params(axis="both", labelsize=7, colors="#222222", length=3)
    if show_y_grid:
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.7)
    else:
        ax.grid(False)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=8)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7)


def add_phase_legend(ax: Any) -> None:
    handles = [
        Patch(facecolor=FRAME_COLORS[frame], edgecolor="none", label=f"phase {frame}")
        for frame in (0, 1, 2)
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7, ncols=3, loc="upper left")
