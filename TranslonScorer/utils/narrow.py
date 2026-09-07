"""Type narrowing for third-party APIs whose stubs are wider than reality.

Two libraries return values typed far more loosely than they behave, and both
are used constantly here:

* **polars** — ``Series.min()``/``.max()``/``.mean()``/``.median()`` are typed
  as a union spanning ``int | float | Decimal | date | time | timedelta | str |
  bytes | ndarray | list | None``, because a Series can hold any of those. Every
  ``int(series.max())`` on a column that is always numeric therefore fails type
  checking.

* **pysam** — ``AlignedSegment.reference_end`` and ``query_name`` are
  ``Optional``, correct in general (an unmapped record has no reference end)
  but not at call sites already inside a mapped/named-record branch.

Narrowing in one place, with the reason written down, is honest about these
being stub limitations. Scattering ``# type: ignore`` through the modules would
hide the same thing less legibly, and would also silence *real* errors on those
lines later.

These are assertions, not coercions: passing a genuinely unexpected value
raises here rather than propagating a wrong number downstream.
"""

from __future__ import annotations

from typing import Any, Optional

__all__ = ["as_int", "as_float", "opt_float", "require"]


def as_int(value: Any) -> int:
    """A polars scalar (or any numeric) as ``int``.

    Raises ``TypeError`` on None — callers must guard emptiness themselves
    (``height``/``is_empty``), because an aggregate over an empty Series is
    None and silently turning that into 0 would be wrong.
    """
    if value is None:
        raise TypeError("expected a numeric value, got None (empty Series?)")
    return int(value)


def as_float(value: Any) -> float:
    """A polars scalar (or any numeric) as ``float``. See :func:`as_int`."""
    if value is None:
        raise TypeError("expected a numeric value, got None (empty Series?)")
    return float(value)


def opt_float(value: Any) -> Optional[float]:
    """A polars scalar as ``float``, preserving None.

    For aggregates where "no value" is a legitimate result to carry forward
    rather than an error — an absent metric is not a zero one.
    """
    return None if value is None else float(value)


def require(value: Optional[Any], what: str) -> Any:
    """Assert an Optional is present, naming it if not.

    For pysam attributes that are Optional in the stubs but guaranteed by the
    branch the call sits in (a mapped record always has ``reference_end``).
    """
    if value is None:
        raise ValueError(f"{what} is unexpectedly None")
    return value
