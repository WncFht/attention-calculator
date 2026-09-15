"""Dispatch /get_integral_image rendering to the kernel family module."""

import importlib
import re
from fractions import Fraction

from .solve import FAMILY

INT_KEYS = ("m", "n", "au_val", "bu_val", "cu_val", "u_val")
FRAC_KEYS = ("a_val", "b_val", "c_val")

# 前端把有理数显示成 LaTeX \frac{n}{d} 再回传给 /get_integral_image（见
# templates/index.html 的 rationalDisplay/piDisplay），站端接受并原样回显
WIRE_FRAC_RE = re.compile(r"^\\d?frac\{\s*(-?\d+)\s*\}\{\s*(-?\d+)\s*\}$")


def wire_fraction(v: Fraction | str) -> Fraction:
    """Numeric value of a wire field: Fraction passthrough, 'n/d', 'n',
    or the frontend's LaTeX '\\frac{n}{d}'/'\\dfrac{n}{d}' display form."""
    if isinstance(v, Fraction):
        return v
    m = WIRE_FRAC_RE.match(v.strip())
    if m:
        return Fraction(int(m.group(1)), int(m.group(2)))
    return Fraction(v)


def wire_or(v: Fraction | str) -> Fraction | None:
    """Tolerant wire_fraction: None when the text isn't a number.

    /get_integral_image never validates ``coef``/``rational`` — the site
    echoes garbage like ``coef=x`` straight into the LaTeX. Display paths
    (``== 1`` checks) use this; paths that genuinely need the value keep
    wire_fraction so bad input still lands in the site's 500 catch-all.
    """
    try:
        return wire_fraction(v)
    except (ValueError, ZeroDivisionError):
        return None


def rat_tex(v: Fraction | str) -> str:
    """Site-style rational LaTeX: integers bare, fractions as ``\\dfrac``.

    A raw string keeps the request's literal text — '3140/1000' renders
    ``\\dfrac{3140}{1000}`` unreduced, and a LaTeX ``\\frac{n}{d}`` input is
    echoed verbatim (the site re-embeds the macro the user typed).
    """
    if isinstance(v, str):
        num, slash, den = v.partition("/")
        return f"\\dfrac{{{num}}}{{{den}}}" if slash else v
    return (str(v.numerator) if v.denominator == 1
            else f"\\dfrac{{{v.numerator}}}{{{v.denominator}}}")


def coerce_params(query: dict) -> dict:
    """Turn the nine site query fields into the typed params kernels expect.

    m, n, au_val, bu_val, cu_val, u_val become int; a_val, b_val, c_val become
    Fraction (the same shapes the family ``prove`` emits and golden.jsonl
    records). Raises ValueError on malformed input.
    """
    params = {k: int(query[k]) for k in INT_KEYS}
    params.update({k: Fraction(query[k]) for k in FRAC_KEYS})
    return params


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str,
    bound: Fraction | str,
) -> str:
    """Render the proof-equation LaTeX for a solved case.

    ``params`` is the typed dict from ``coerce_params``. ``power``/``bound``
    may be the request's raw text — the site echoes the digits unreduced
    (\\dfrac{314}{100}) — so callers should pass the original strings.
    Returns the site's /get_integral_image ``equation`` string.
    """
    if kind not in FAMILY:
        raise ValueError(f"unsupported type {kind!r}")
    module = importlib.import_module(f"attention_calculator.kernels.{FAMILY[kind]}")
    return module.render_equation(params, kind, power, comp, bound)
