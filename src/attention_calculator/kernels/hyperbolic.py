"""hyperbolic kernels: sinh_q, cosh_q, tanh_q, coth_q on [0,1].

Basis x^m(1-x)^n, factor a+bx+cx^2, kernel sinh(qx); moments live in
span{sinh q, cosh q, 1}. tanh/coth are solved as sinh q - p cosh q
(resp. cosh q - p sinh q) = integral, then divided by cosh q (sinh q)
in the rendered equation (kernel-spec.md).
"""

import math
from fractions import Fraction
from math import comb, lcm

import sympy as sp

from ..moment import Moment, add, combine, scale
from ..engine import WrongDirection, mn_order, search
from ..render import rat_tex, wire_fraction
from .exp_family import emit, mul_latex

LIMIT = 10


def sinh_moments(q: Fraction, upto: int) -> tuple[list[Moment], list[Moment]]:
    """J_k = ∫ x^k sinh(qx), C_k = ∫ x^k cosh(qx) on [0,1], k = 0..upto.

    IBP: J_k = cosh(q)/q - (k/q) C_{k-1}, C_k = sinh(q)/q - (k/q) J_{k-1}.
    """
    J: list[Moment] = [{"cosh_q": Fraction(1) / q, "1": Fraction(-1) / q}]
    C: list[Moment] = [{"sinh_q": Fraction(1) / q}]
    for k in range(1, upto + 1):
        J.append(add({"cosh_q": Fraction(1) / q}, scale(C[-1], Fraction(-k) / q)))
        C.append(add({"sinh_q": Fraction(1) / q}, scale(J[-2], Fraction(-k) / q)))
    return J, C


def basis_moment(m: int, n: int, j: int, q: Fraction) -> Moment:
    """∫_0^1 x^{m+j}(1-x)^n sinh(qx) dx = sum_i (-1)^i C(n,i) J_{m+j+i}."""
    J, C = sinh_moments(q, m + j + n)
    coeffs = [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)]
    return combine(coeffs, J[m + j : m + j + n + 1])


def const_latex(kind: str, power: Fraction | str) -> str:
    """\\sinh1 / \\sinh\\dfrac{2}{3} etc."""
    name = {"sinh_q": "sinh", "cosh_q": "cosh", "tanh_q": "tanh", "coth_q": "coth"}[kind]
    return f"\\{name}" + rat_tex(power)


def target_for(kind: str, comp: str, bound: Fraction) -> Moment:
    """Moment vector the integral must equal.

    sinh_q: ±(sinh q - bound); cosh_q: ±(cosh q - bound);
    tanh_q: ±(sinh q - bound*cosh q); coth_q: ±(cosh q - bound*sinh q)
    (sign + for '>', - for '<'; cosh q > 0 and sinh q > 0 for q > 0).
    """
    sign = 1 if comp == ">" else -1
    if kind == "sinh_q":
        return {"sinh_q": Fraction(sign), "1": Fraction(-sign) * bound}
    if kind == "cosh_q":
        return {"cosh_q": Fraction(sign), "1": Fraction(-sign) * bound}
    if kind == "tanh_q":
        return {"sinh_q": Fraction(sign), "cosh_q": Fraction(-sign) * bound}
    return {"cosh_q": Fraction(sign), "sinh_q": Fraction(-sign) * bound}


# the site evaluates the claimed constant in float64 before the kernel runs —
# a strictly-false inequality is 方向反了 without ever touching the moment
# machinery (q=0 would crash it); equality or truth proceeds. math.* overflow
# (sinh 711, cosh 1e15) and the coth 1/tanh(0) division surface as 500.
CONST_F = {"sinh_q": math.sinh, "cosh_q": math.cosh, "tanh_q": math.tanh,
           "coth_q": lambda v: 1 / math.tanh(v)}


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """prove(kind, power, comp, bound) -> site /calculate shape.

    power==0: the float64 pre-check decides first — false claims get 方向反了
    (sinh 0>9, cosh 0<1/2), true/equal ones reach the kernel's 1/q crash
    (sinh 0<0, cosh 0>1/2 -> 500); coth crashes inside the check itself.
    """
    c = CONST_F[kind](float(power))
    if (float(bound) > c) if comp == ">" else (float(bound) < c):
        raise WrongDirection
    plans = ((m, n, [basis_moment(m, n, j, power) for j in (0, 1, 2)])
             for m, n in mn_order(LIMIT))
    solved = search(plans, target_for(kind, comp, bound), True)
    return emit(kind, solved.m, solved.n, solved.coeffs)


def render_equation(params: dict, kind: str, power: Fraction | str,
                    comp: str, bound: Fraction | str) -> str:
    """Rebuild the site's get_integral_image LaTeX for solved parameters."""
    x = sp.symbols("x")
    m, n = params["m"], params["n"]
    au, bu, cu, u = (params["au_val"], params["bu_val"], params["cu_val"], params["u_val"])
    integrand = x**m * (1 - x) ** n * (au + bu * x + cu * x**2) * sp.sinh(wire_fraction(power) * x)
    body = mul_latex(integrand, u)
    if kind == "tanh_q":
        inner = f"\\dfrac{{1}}{{\\cosh({power})}} " + body
    elif kind == "coth_q":
        inner = f"\\dfrac{{1}}{{\\sinh({power})}} " + body
    else:
        inner = body
    const = const_latex(kind, power)
    r = rat_tex(bound)
    lhs = f"{const} - {r}" if comp == ">" else f"{r} - {const}"
    return f"{lhs} = \\int_0^1 {inner} \\mathrm{{d}} x > 0"
