"""Dispatch /get_integral_image rendering to the kernel family module."""

import importlib
import re
from fractions import Fraction
from itertools import pairwise
from math import lcm

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


# coef 位置解析宏时丢掉一切符号：'-\dfrac{6}{4}' 与 '\dfrac{6}{-4}' 都解析成
# (6, 4)（probe: e_q 被积 e^{3x/2}、tan_q 除数 \cos(6/4)）
WIRE_PAIR_MACRO_RE = re.compile(r"\\d?frac\{\s*-?(\d+)\s*\}\{\s*-?(\d+)\s*\}")


def wire_pair(v: Fraction | str) -> tuple[int, int]:
    """Raw (numerator, denominator) the site's render layer derives from a
    wire field — *unreduced* and per-notation quirky:

    - ``\\frac``/``\\dfrac`` macro: unsigned digit groups; a minus inside or
      outside the braces is dropped (e_q '-\\dfrac{6}{4}' -> e^{3x/2}).
    - plain 'n/d' text: signed and unreduced, so ln_q's (q-1) denominator
      treats '4/2' and '2' differently.
    - Fraction passthrough (already reduced — the /calculate path).

    Raises ValueError on unparseable text, like wire_fraction: the site
    echoes garbage into display-only slots but crashes (500) when a kernel
    genuinely needs the value.
    """
    if isinstance(v, Fraction):
        return v.numerator, v.denominator
    m = WIRE_PAIR_MACRO_RE.search(v)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.fullmatch(r"(-?\d+)(?:/(-?\d+))?", v.strip())
    if not m:
        raise ValueError(f"unparseable rational {v!r}")
    return int(m.group(1)), int(m.group(2) or 1)


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
        # raw 'n/1' integerizes to 'n'; 'n/n' with den != 1 stays a \dfrac
        if den == "1":
            return num
        return f"\\dfrac{{{num}}}{{{den}}}" if slash else v
    return (
        str(v.numerator) if v.denominator == 1 else f"\\dfrac{{{v.numerator}}}{{{v.denominator}}}"
    )


def coef_tex(v: Fraction | str) -> str:
    """Multiplier-position text: literal '1' (or Fraction 1) vanishes —
    ``\\pi`` not ``1\\pi`` — while '1/1' still prints ``1\\pi``: the site only
    strips the canonical form (probe: coef '1/1' -> '1C', '1e^{\\pi}').
    '0', '-1', '2/2' and junk all pass through rat_tex's raw echo."""
    return "" if v == "1" or v == Fraction(1) else rat_tex(v)


def join_cdot(tex: list[str]) -> str:
    """Join factor texts the way the site's sympy does: `` \\cdot `` iff the
    left piece ends in '}' and the right is a ``\\left(<digit>...\\right)``
    group — a parenthesized Add carrying no outer exponent."""
    out = tex[0]
    for prev, cur in pairwise(tex):
        cdot = (
            prev.endswith("}")
            and cur.startswith("\\left(")
            and cur[6].isdigit()
            and cur.endswith("\\right)")
        )
        out += " \\cdot " if cdot else " "
        out += cur
    return out


def cdot_tex(tex: str) -> str:
    r"""The site's cdot-spacing fixup on an already-latexed product string:
    `` \\cdot `` before a digit-leading ``\left(`` group when the previous
    token ends in a digit or '}' (trig_pi integrand, beta last-resort block)."""
    return re.sub(r"(?<=[0-9}]) (?=\\left\(\d)", r" \\cdot ", tex)


def lhs_tex(const_tex: str, bound_tex: str, comp: str) -> str:
    """The equation's left-hand side: ``C - r`` for '>', ``r - C`` for '<'."""
    return f"{const_tex} - {bound_tex}" if comp == ">" else f"{bound_tex} - {const_tex}"


def coerce_params(query: dict) -> dict:
    """Turn the nine site query fields into the typed params kernels expect.

    m, n, au_val, bu_val, cu_val, u_val become int; a_val, b_val, c_val become
    Fraction (the same shapes the family ``prove`` emits and golden.jsonl
    records). Raises ValueError on malformed input.
    """
    params = {k: int(query[k]) for k in INT_KEYS}
    params.update({k: Fraction(query[k]) for k in FRAC_KEYS})
    return params


def emit(m: int, n: int, coeffs: list[Fraction], c_val: str | None = None) -> dict:
    """Assemble the site's /calculate payload: the ``parameters`` dict plus the
    ``equations.solution`` string — the inverse direction of coerce_params.

    ``c_val`` overrides the reported c_val where the field is repurposed
    (trig_pi's alpha, artanh/arcoth's reduced q~); ``cu_val`` always carries
    the solved c coefficient.
    """
    a, b = coeffs[0], coeffs[1]
    c = coeffs[2] if len(coeffs) > 2 else Fraction(0)
    u = lcm(a.denominator, b.denominator, c.denominator)
    params = {
        "m": m,
        "n": n,
        "a_val": str(a),
        "b_val": str(b),
        "c_val": str(c) if c_val is None else c_val,
        "au_val": str(a * u),
        "bu_val": str(b * u),
        "cu_val": str(c * u),
        "u_val": str(u),
        "unified_form": {},
    }
    solution = f"a = {a}, b = {b}" + ("" if len(coeffs) <= 2 else f", c= {c}")
    return {"parameters": params, "solution": solution}


def render_equation(
    params: dict,
    kind: str,
    power: Fraction | str,
    comp: str,
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
