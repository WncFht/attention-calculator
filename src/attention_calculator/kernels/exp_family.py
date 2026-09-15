"""exp kernels: ``e`` (coef*e), ``e_q`` (e^q), ``e_pi`` (coef*e^pi).

``e``/``e_q`` prove on [0,1] with basis x^m(1-x)^n and factor a+bx;
``e_pi`` proves on [0,pi] with basis sin^m(x)(1-sin x)^n and factor
a+b*sin(x). All moments live in span{const, 1} (kernel-spec.md).
"""

from fractions import Fraction
from math import comb, lcm

import sympy as sp

from ..engine import search
from ..moment import Moment, add, combine, scale

LIMIT_E = 30  # e、pi 两类型指数上限 30
LIMIT_OTHER = 10


def mn_order(limit: int):
    """(m,n) in the site's order: m+n asc, |m-n| asc, then m asc.

    engine.mn_order breaks |m-n| ties towards larger m; the live site does
    the opposite (e.g. e_pi power=3/4 picks (1,2) over (2,1), cosh_q picks
    (0,1) over (1,0)). Kept local so family parity does not depend on the
    shared engine ordering.
    """
    for s in range(2 * limit + 1):
        pairs = [(m, s - m) for m in range(s + 1) if s - m <= limit and m <= limit]
        pairs.sort(key=lambda p: (abs(p[0] - p[1]), p[0]))
        yield from pairs


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
        l = (j - 1) // 2
        for k in range(l + 1):
            c = Fraction((-1) ** k * comb(2 * l + 1, l - k), 4**l)
            kk = 2 * k + 1  # ∫ e^x sin(kk x) dx = kk(1+e^pi)/(1+kk^2), kk odd
            out = add(out, {sym: c * kk / (1 + kk * kk), "1": c * kk / (1 + kk * kk)})
    else:  # sin^{2l} x = 4^{-l} [C(2l,l) + 2 sum_k (-1)^k C(2l,l-k) cos(2kx)]
        l = j // 2
        out = {sym: Fraction(comb(2 * l, l), 4**l), "1": -Fraction(comb(2 * l, l), 4**l)}
        for k in range(1, l + 1):
            c = Fraction(2 * (-1) ** k * comb(2 * l, l - k), 4**l)
            # ∫ e^x cos(2kx) dx = (e^pi-1)/(1+4k^2)
            out = add(out, {sym: c / (1 + 4 * k * k), "1": -c / (1 + 4 * k * k)})
    return out


def basis_sin_moment(m: int, n: int, j: int) -> Moment:
    """∫_0^pi sin^{m+j} x (1-sin x)^n e^x dx = sum_i (-1)^i C(n,i) S_{m+j+i}."""
    coeffs = [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)]
    return combine(coeffs, [exp_sin_moment(m + j + i) for i in range(n + 1)])


def frac_latex(r: Fraction) -> str:
    """Site-style rational LaTeX: integers bare, fractions as \\dfrac."""
    if r.denominator == 1:
        return str(r.numerator)
    return f"\\dfrac{{{r.numerator}}}{{{r.denominator}}}"


def const_latex(kind: str, power: Fraction) -> str:
    """Left-side constant text: e.g. e, 2e, \\dfrac{1}{2}e, e^2, e^\\dfrac{1}{2}, e^{\\pi}."""
    if kind == "e":
        return "e" if power == 1 else f"{frac_latex(power)}e"
    if kind == "e_q":
        return "e^" + frac_latex(power)
    return ("e^{\\pi}" if power == 1 else f"{frac_latex(power)}e^{{\\pi}}")


def emit(kind: str, m: int, n: int, coeffs: list[Fraction]) -> dict:
    """Assemble the site's parameters dict and equations.solution string."""
    a, b = coeffs[0], coeffs[1]
    c = coeffs[2] if len(coeffs) > 2 else Fraction(0)
    u = lcm(a.denominator, b.denominator, c.denominator)
    params = {
        "m": m,
        "n": n,
        "a_val": str(a),
        "b_val": str(b),
        "c_val": str(c),
        "au_val": int(a * u),
        "bu_val": int(b * u),
        "cu_val": int(c * u),
        "u_val": int(u),
    }
    solution = f"a = {a}, b = {b}" + ("" if len(coeffs) <= 2 else f", c= {c}")
    return {"parameters": params, "solution": solution}


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """prove(kind, power, comp, bound) -> site /calculate shape."""
    if power <= 0:
        raise ValueError("左侧系数格式无效")
    if bound <= 0:
        raise ValueError("右侧有理数格式无效")
    if kind == "e_q":
        sym, coef, q, limit = "e_q", Fraction(1), power, LIMIT_OTHER
        plans = ((m, n, [basis_x_moment(m, n, j, q, sym) for j in (0, 1)])
                 for m, n in mn_order(limit))
    elif kind == "e":
        sym, coef, q, limit = "e", power, Fraction(1), LIMIT_E
        plans = ((m, n, [basis_x_moment(m, n, j, q, sym) for j in (0, 1)])
                 for m, n in mn_order(limit))
    else:  # e_pi
        sym, coef, limit = "e_pi", power, LIMIT_OTHER
        plans = ((m, n, [basis_sin_moment(m, n, j) for j in (0, 1)])
                 for m, n in mn_order(limit))
    sign = 1 if comp == ">" else -1
    target = {sym: sign * coef, "1": -sign * bound}
    solved = search(plans, target, True)
    return emit(kind, solved.m, solved.n, solved.coeffs)


def render_equation(params: dict, kind: str, power: Fraction, comp: str, bound: Fraction) -> str:
    """Rebuild the site's get_integral_image LaTeX for solved parameters."""
    x = sp.symbols("x")
    m, n = params["m"], params["n"]
    au, bu, u = params["au_val"], params["bu_val"], params["u_val"]
    if kind == "e_pi":
        s = sp.sin(x)
        integrand = s**m * (1 - s) ** n * (au + bu * s) * sp.exp(x)
    else:
        q = power if kind == "e_q" else sp.Integer(1)
        integrand = x**m * (1 - x) ** n * (au + bu * x) * sp.exp(q * x)
    if u != 1:
        integrand = integrand / u
    const = const_latex(kind, power)
    r = frac_latex(bound)
    lhs = f"{const} - {r}" if comp == ">" else f"{r} - {const}"
    upper = "\\pi" if kind == "e_pi" else "1"
    return f"{lhs} = \\int_0^{{{upper}}} {sp.latex(integrand)} \\mathrm{{d}} x > 0"
