"""hyperbolic kernels: sinh_q, cosh_q, tanh_q, coth_q on [0,1].

Basis x^m(1-x)^n, factor a+bx+cx^2, kernel sinh(qx); moments live in
span{sinh q, cosh q, 1}. tanh/coth are solved as sinh q - p cosh q
(resp. cosh q - p sinh q) = integral, then divided by cosh q (sinh q)
in the rendered equation (kernel-spec.md).
"""

import math
from fractions import Fraction
from math import comb

import sympy as sp

from ..engine import WrongDirection, mn_order, search
from ..moment import Moment, add, combine, scale
from ..render import lhs_tex, rat_tex, wire_fraction
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
    J = sinh_moments(q, m + j + n)[0]
    coeffs = [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)]
    return combine(coeffs, J[m + j : m + j + n + 1])


def const_tex(kind: str, power: Fraction | str) -> str:
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
CONST_F = {
    "sinh_q": math.sinh,
    "cosh_q": math.cosh,
    "tanh_q": math.tanh,
    "coth_q": lambda v: 1 / math.tanh(v),
}


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = False) -> dict:
    """prove(kind, power, comp, bound) -> site /calculate shape.

    power==0: site mode keeps the probed float64-then-1/q-crash behavior;
    exact mode rejects it as a domain error — the kernel collapses there
    (IBP base terms carry 1/q) and the constants are rational (sinh 0 =
    tanh 0 = 0, cosh 0 = 1) or singular (coth 0).

    q<0 (exact only): kernel sinh(qx) is non-positive on [0,1], so a
    non-negative P can only certify the flipped inequality — the site
    path misreports true claims as 方向反了. Reduce by parity instead:
    sinh/tanh are odd so the (|q|, flipped comp, -bound) problem is
    equivalent; cosh is even so only the argument changes. The emitted
    P is negated: (-P)·sinh(qx) >= 0 on [0,1] and its integral is the
    original s·(C - bound) vector. coth_q stays rejected for q<0 — its
    printed 1/sinh(q) prefactor flips the certified sign and would make
    the rendered equation lie.
    """
    if not exact:  # exact mode: direction already certified upstream
        c = CONST_F[kind](float(power))
        if (float(bound) > c) if comp == ">" else (float(bound) < c):
            raise WrongDirection
    q, comp_s, bound_s, negate = power, comp, bound, False
    if exact and power <= 0:
        if kind == "coth_q":  # coth 0 singular; q<0 breaks the 1/sinh q print sign
            raise ValueError("请在coth后输入一个大于0的数")
        if power == 0:
            raise ValueError(f"请在{kind.split('_')[0]}后输入一个非0的数")
        odd = kind in ("sinh_q", "tanh_q")
        q, comp_s = -power, ("<" if comp == ">" else ">") if odd else comp
        bound_s, negate = (-bound if odd else bound), True
    plans = ((m, n, [basis_moment(m, n, j, q) for j in (0, 1, 2)]) for m, n in mn_order(LIMIT))
    solved = search(plans, target_for(kind, comp_s, bound_s), True)
    coeffs = [-v for v in solved.coeffs] if negate else solved.coeffs
    return emit(kind, solved.m, solved.n, coeffs)


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
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
    lhs = lhs_tex(const_tex(kind, power), rat_tex(bound), comp)
    return f"{lhs} = \\int_0^1 {inner} \\mathrm{{d}} x > 0"
