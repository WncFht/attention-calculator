"""exp kernels: ``e`` (coef*e), ``e_q`` (e^q), ``e_pi`` (coef*e^pi).

``e``/``e_q`` prove on [0,1] with basis x^m(1-x)^n and factor a+bx;
``e_pi`` proves on [0,pi] with basis sin^m(x)(1-sin x)^n and factor
a+b*sin(x). All moments live in span{const, 1} (kernel-spec.md).
"""

import math
from fractions import Fraction
from math import comb

import sympy as sp

from .. import render
from ..engine import WrongDirection, mn_order, search
from ..moment import Moment, add, combine, scale
from ..render import coef_tex, lhs_tex, rat_tex, wire_pair

LIMIT_E = 30  # e、pi 两类型指数上限 30
LIMIT_OTHER = 10


def exp_x_moments(q: Fraction, sym: str, upto: int) -> list[Moment]:
    """I_k = ∫_0^1 x^k e^{q x} dx for k = 0..upto via IBP I_k = e^q/q - k/q I_{k-1}."""
    out = [{sym: Fraction(1, 1) / q, "1": -Fraction(1, 1) / q}]
    for k in range(1, upto + 1):
        out.append(add({sym: Fraction(1) / q}, scale(out[-1], Fraction(-k) / q)))
    return out


def basis_x_moment(m: int, n: int, j: int, q: Fraction, sym: str) -> Moment:
    """∫_0^1 x^{m+j}(1-x)^n e^{q x} dx = sum_i (-1)^i C(n,i) I_{m+j+i}."""
    imoms = exp_x_moments(q, sym, m + j + n)
    coeffs = [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)]
    return combine(coeffs, imoms[m + j : m + j + n + 1])


def exp_sin_moment(j: int, sym: str = "e_pi") -> Moment:
    """S_j = ∫_0^pi e^x sin^j x dx via multiple-angle expansion of sin^j."""
    if j == 0:
        return {sym: Fraction(1), "1": Fraction(-1)}
    out: Moment = {}
    if j % 2 == 1:  # sin^{2l+1} x = 4^{-l} sum_k (-1)^k C(2l+1,l-k) sin((2k+1)x)
        ell = (j - 1) // 2
        for k in range(ell + 1):
            c = Fraction((-1) ** k * comb(2 * ell + 1, ell - k), 4**ell)
            kk = 2 * k + 1  # ∫ e^x sin(kk x) dx = kk(1+e^pi)/(1+kk^2), kk odd
            out = add(out, {sym: c * kk / (1 + kk * kk), "1": c * kk / (1 + kk * kk)})
    else:  # sin^{2l} x = 4^{-l} [C(2l,l) + 2 sum_k (-1)^k C(2l,l-k) cos(2kx)]
        ell = j // 2
        out = {
            sym: Fraction(comb(2 * ell, ell), 4**ell),
            "1": -Fraction(comb(2 * ell, ell), 4**ell),
        }
        for k in range(1, ell + 1):
            c = Fraction(2 * (-1) ** k * comb(2 * ell, ell - k), 4**ell)
            # ∫ e^x cos(2kx) dx = (e^pi-1)/(1+4k^2)
            out = add(out, {sym: c / (1 + 4 * k * k), "1": -c / (1 + 4 * k * k)})
    return out


def basis_sin_moment(m: int, n: int, j: int) -> Moment:
    """∫_0^pi sin^{m+j} x (1-sin x)^n e^x dx = sum_i (-1)^i C(n,i) S_{m+j+i}."""
    coeffs = [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)]
    return combine(coeffs, [exp_sin_moment(m + j + i) for i in range(n + 1)])


def mul_latex(numer: sp.Expr, u: int) -> str:
    """Site-style LaTeX of the integrand numerator / u.

    The site puts `\\cdot` between a Pow factor and a following parenthesized
    Add factor only when the Add's own text starts with a digit:
    `x^{2} \\cdot \\left(1 - x\\right)` and `(1-x)^{2} \\cdot \\left(8 x + 3\\right)`,
    but `x \\left(1 - x\\right)^{2} \\left(x + 1\\right)` (golden byte-diff).
    Factor order is sympy's as_ordered_factors.
    """
    terms = numer.as_ordered_factors()
    parts = []
    for i, t in enumerate(terms):
        tex = sp.latex(t)
        cdot = terms[i - 1].is_Pow and isinstance(t, sp.Add) and tex[0].isdigit() if i else False
        if isinstance(t, sp.Add):
            tex = rf"\left({tex}\right)"
        if i:
            parts.append(" \\cdot " if cdot else " ")
        parts.append(tex)
    if u == 1:
        return "".join(parts)
    return rf"\frac{{{''.join(parts)}}}{{{u}}}"


def const_tex(kind: str, power: Fraction | str) -> str:
    """Left-side constant text: e.g. e, 2e, \\dfrac{1}{2}e, e^2, e^\\dfrac{1}{2}, e^{\\pi}."""
    if kind == "e":
        return coef_tex(power) + "e"
    if kind == "e_q":
        return "e^" + rat_tex(power)
    return coef_tex(power) + "e^{\\pi}"


def emit(kind: str, m: int, n: int, coeffs: list[Fraction]) -> dict:
    """Assemble the site's parameters dict and equations.solution string.

    Same assembly as render.emit but with the field shape tests/test_exp.py
    pins for this family: au/bu/cu/u_val stay ints and no unified_form key
    (server.calculate's setdefault restores it on the wire).
    """
    res = render.emit(m, n, coeffs)
    params = res["parameters"]
    del params["unified_form"]
    for k in ("au_val", "bu_val", "cu_val", "u_val"):
        params[k] = int(params[k])
    return res


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = False) -> dict:
    """prove(kind, power, comp, bound) -> site /calculate shape.

    Input format is validated one layer up (server.NUM_RE); power==0 for
    e_q reaches the natural 1/q division (site answers 500 there).
    """
    if kind == "e_q":
        # float64 direction pre-check before the kernel: false claims are
        # 方向反了 (e^0 > 2 -> 404), true/equal ones proceed into the 1/q
        # moments which crash for q=0 (-> 500); math.exp overflow -> 500
        c = math.exp(float(power))
        if (float(bound) > c) if comp == ">" else (float(bound) < c):
            raise WrongDirection
        sym, coef, q, limit = "e_q", Fraction(1), power, LIMIT_OTHER
        plans = (
            (m, n, [basis_x_moment(m, n, j, q, sym) for j in (0, 1)]) for m, n in mn_order(limit)
        )
    elif kind == "e":
        sym, coef, q, limit = "e", power, Fraction(1), LIMIT_E
        plans = (
            (m, n, [basis_x_moment(m, n, j, q, sym) for j in (0, 1)]) for m, n in mn_order(limit)
        )
    else:  # e_pi
        sym, coef, limit = "e_pi", power, LIMIT_OTHER
        plans = ((m, n, [basis_sin_moment(m, n, j) for j in (0, 1)]) for m, n in mn_order(limit))
    sign = 1 if comp == ">" else -1
    target = {sym: sign * coef, "1": -sign * bound}
    solved = search(plans, target, True)
    return emit(kind, solved.m, solved.n, solved.coeffs)


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Rebuild the site's get_integral_image LaTeX for solved parameters."""
    x = sp.symbols("x")
    m, n = params["m"], params["n"]
    au, bu, u = params["au_val"], params["bu_val"], params["u_val"]
    if kind == "e_pi":
        s = sp.sin(x)
        integrand = s**m * (1 - s) ** n * (au + bu * s) * sp.exp(x)
    else:
        # e_q's exponent comes from the raw wire pair: a macro parses with
        # unsigned digit groups ('\dfrac{6}{-4}' -> e^{3x/2}, sign dropped)
        q = Fraction(*wire_pair(power)) if kind == "e_q" else sp.Integer(1)
        integrand = x**m * (1 - x) ** n * (au + bu * x) * sp.exp(q * x)
    lhs = lhs_tex(const_tex(kind, power), rat_tex(bound), comp)
    upper = "{\\pi}" if kind == "e_pi" else "1"
    return f"{lhs} = \\int_0^{upper} {mul_latex(integrand, u)} \\mathrm{{d}} x > 0"
