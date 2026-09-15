"""Radian trig kernels on [0, 1]: sin_q, cos_q, tan_q, cot_q.

The integrand is x^m (1-x)^n (a + b x + c x^2) sin(q x). With
S_k = int_0^1 x^k sin(qx) dx and C_k = int_0^1 x^k cos(qx) dx, integrating by
parts twice gives

    S_0 = (1 - cos q)/q,
    S_k = -cos q/q + k sin q/q^2 - k(k-1)/q^2 * S_{k-2},

so every basis moment stays in span{sin q, cos q, 1} over QQ(q).

tan/cot first clear the target's denominator (article types 18/19):

    tan q ~ p  <=>  (sin q - p cos q)/cos q ~ 0
    cot q ~ p  <=>  (cos q - p sin q)/sin q ~ 0

and the proof is displayed as ``tan q - p = (int f)/cos q > 0`` with the
constant divisor pulled in front of the integral. The site only admits
0 < q < pi/2 for tan/cot, so the divisor is positive and both directions
reduce to a nonnegative polynomial factor. Ranges on the site:

    sin_q: 0 < q < pi;  cos_q, tan_q, cot_q: 0 < q < pi/2.

Input-format failures are raised as ValueError carrying the site's Chinese
message; the API layer maps it onto the error response.
"""

from __future__ import annotations

import re
from fractions import Fraction
from math import comb, lcm

import sympy as sp

from ..engine import mn_order, search
from ..moment import Moment, combine
from ..render import rat_tex, wire_pair

LIMIT = 10

x = sp.symbols("x")


def sin_moments(kmax: int, q: Fraction) -> list[Moment]:
    """S_k = int_0^1 x^k sin(qx) dx for k = 0..kmax, over {sin_q, cos_q, 1}."""
    s: list[Moment] = [
        {"cos_q": Fraction(-1, q), "1": Fraction(1, q)},
        {"sin_q": Fraction(1, q * q), "cos_q": Fraction(-1, q)},
    ]
    for k in range(2, kmax + 1):
        f = Fraction(-k * (k - 1), q * q)
        s.append(
            {
                "sin_q": Fraction(k, q * q) + f * s[k - 2].get("sin_q", Fraction(0)),
                "cos_q": Fraction(-1, q) + f * s[k - 2].get("cos_q", Fraction(0)),
                "1": f * s[k - 2].get("1", Fraction(0)),
            }
        )
    return [{k: v for k, v in m.items() if v} for m in s]


def basis_moment(s: list[Moment], m: int, n: int, j: int) -> Moment:
    """int_0^1 x^{m+j} (1-x)^n sin(qx) dx = sum_i (-1)^i C(n,i) S_{m+i+j}."""
    return combine(
        [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)],
        [s[m + i + j] for i in range(n + 1)],
    )


def check_input(q: Fraction, bound: Fraction, kind: str) -> None:
    """Reproduce the site's input validation order and messages."""
    if q < 0:
        raise ValueError("左侧系数格式无效")
    if bound < 0:
        raise ValueError("右侧有理数格式无效")
    if kind == "sin_q":
        name, upper_tex, upper = "sin", "π", sp.pi
    else:
        name, upper_tex, upper = kind.removesuffix("_q"), "π/2", sp.pi / 2
    if q <= 0 or not sp.Rational(q) < upper:
        raise ValueError(f"请在{name}后输入一个在(0,{upper_tex})内的数")


def solve(kind: str, q: Fraction, comp: str, bound: Fraction):
    """Find (m, n, a, b, c) whose integrand proves the requested inequality.

    The signed target is s*(A - p*B) with s = +1 for ">", where (A, B) is
    (sin q, 1), (cos q, 1), (sin q, bound*cos q), (cos q, bound*sin q) for
    sin/cos/tan/cot respectively -- so the integrand must be nonnegative.
    """
    check_input(q, bound, kind)
    s = Fraction(1 if comp == ">" else -1)
    if kind == "sin_q":
        target = {"sin_q": s, "1": -s * bound}
    elif kind == "cos_q":
        target = {"cos_q": s, "1": -s * bound}
    elif kind == "tan_q":
        target = {"sin_q": s, "cos_q": -s * bound}
    else:
        target = {"cos_q": s, "sin_q": -s * bound}

    moments = sin_moments(2 * LIMIT + 2, q)
    plans = ((m, n, [basis_moment(moments, m, n, j) for j in range(3)]) for m, n in mn_order(LIMIT))
    return search(plans, target, True)


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """Run the proof search for one trig_q-family request."""
    solved = solve(kind, power, comp, bound)
    a, b, c = solved.coeffs
    u = lcm(a.denominator, b.denominator, c.denominator)
    params = {
        "m": solved.m,
        "n": solved.n,
        "a_val": str(a),
        "b_val": str(b),
        "c_val": str(c),
        "au_val": str(a * u),
        "bu_val": str(b * u),
        "cu_val": str(c * u),
        "u_val": str(u),
        "unified_form": {},
    }
    return {
        "type": kind,
        "parameters": params,
        "solution": f"a = {a}, b = {b}, c= {c}",
    }


def raw_ratio(v: Fraction | str) -> str:
    """Plain 'n/d' text for the tan/cot divisor, e.g. \\cos(6/5).

    A \\frac/\\dfrac macro collapses to its unsigned digit groups joined by
    '/' ('\\dfrac{-6}{4}' -> '6/4'); any other raw string echoes verbatim
    ('-1000...0' -> '\\cos(-1000...0)')."""
    if isinstance(v, str):
        m = re.search(r"\\d?frac\{-?(\d+)\}\{-?(\d+)\}", v)
        return f"{m.group(1)}/{m.group(2)}" if m else v
    return f"{v.numerator}/{v.denominator}" if v.denominator != 1 else str(v.numerator)


def join_product(terms: list[tuple[str, str]]) -> str:
    """Join factor latex the way the site's sympy version does.

    Each entry is (check_tex, shown_tex): check_tex is the term text before
    bracket wrapping, which is what the between-two-numbers rule inspects --
    " \\cdot " goes between a factor ending in a digit (possibly followed by
    '}'/spaces) and a factor starting with a digit or \\frac{d}{d}; everything
    else gets a plain space.
    """
    out = ""
    prev = ""
    for check, shown in terms:
        if out:
            if re.search(r"[0-9][} ]*$", prev) and re.match(r"\d|\\frac{\d+}{\d+}", check):
                out += " \\cdot "
            else:
                out += " "
        out += shown
        prev = shown
    return out


def add_term(expr: sp.Expr) -> tuple[str, str]:
    """(check, shown) for a polynomial factor: wrapped in \\left( \\right)."""
    inner = sp.latex(expr)
    return inner, rf"\left({inner}\right)" if expr.is_Add else inner


def pow_term(expr: sp.Expr, e: int) -> tuple[str, str]:
    """(check, shown) for a basis power: \\left( \\right) even at e == 1.

    At e == 1 the check text stays the bare Add latex (e.g. ``1 - x``): the
    site's sympy applies the between-two-numbers \\cdot rule to the term
    before bracket wrapping, so ``x^{9} \\cdot \\left(1 - x\\right)`` appears.
    """
    inner = sp.latex(expr**e)
    if e == 1:
        return inner, rf"\left({inner}\right)"
    return inner, inner


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Reproduce the site's /get_integral_image LaTeX for this family.

    ``power``/``bound`` may be the raw request strings (echoed unreduced like
    the site does) or already-parsed Fractions.
    """
    m, n = params["m"], params["n"]
    au, bu, cu = (Fraction(params[k]) for k in ("au_val", "bu_val", "cu_val"))
    u = Fraction(params["u_val"])
    # the sin factor takes |q| from the raw wire pair (macro digits unsigned);
    # a plain-signed negative q hoists a '- ' in front of the whole fraction
    qn, qd = wire_pair(power)
    q = Fraction(qn, qd)
    name = kind.removesuffix("_q")

    const = rf"\{name}{rat_tex(power)}"
    lhs = f"{const} - {rat_tex(bound)}" if comp == ">" else f"{rat_tex(bound)} - {const}"

    terms = []
    poly = au + bu * x + cu * x**2
    # sympy absorbs a factor 1; any other constant prints as leading coeff
    if poly.is_Number and poly != 1:
        terms.append((t := sp.latex(poly), t))
    if m:
        terms.append((t := sp.latex(x**m), t))
    if n:
        terms.append(pow_term(1 - x, n))
    if not poly.is_Number:
        terms.append(add_term(poly))
    terms.append((t := sp.latex(sp.sin(sp.Rational(abs(q).numerator, abs(q).denominator) * x)), t))
    num = join_product(terms)
    body = num if u == 1 else rf"\frac{{{num}}}{{{u}}}"
    if q < 0:
        body = "- " + body
    if name == "tan":
        body = rf"\dfrac{{1}}{{\cos({raw_ratio(power)})}} " + body
    elif name == "cot":
        body = rf"\dfrac{{1}}{{\sin({raw_ratio(power)})}} " + body
    return lhs + rf" = \int_0^1 {body} \mathrm{{d}} x > 0"
