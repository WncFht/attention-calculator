"""invhyp kernel: arsinh_q — arsinh(q) vs a rational bound, exact-mode only.

arsinh(q) = int_0^q dt/sqrt(1+t^2) rescales onto [0, 1] via t = q x:

    arsinh(q) = q * int_0^1 dx / sqrt(1 + q^2 x^2).

The integrand is x^m (1-x)^n (a + b x + c x^2) / sqrt(1 + q^2 x^2). With
M_k = int_0^1 x^k / sqrt(1+q^2 x^2) dx, differentiating x^{k-1} sqrt(1+q^2 x^2)
gives

    sqrt(1+q^2) = (k-1) M_{k-2} + k q^2 M_k        (k >= 2),

so from M_0 = arsinh(q)/q and M_1 = (sqrt(1+q^2) - 1)/q^2 every basis moment
stays in span{arsinh(q), sqrt(1+q^2), 1} over QQ(q): a 3-dim space whose
parasitic sqrt(1+q^2) is eliminated by the target's zero coefficient — the
same pattern as sin_q's cos_q — with 3 unknowns in P = a + b x + c x^2.
The recurrence was verified against mpmath.quad at 200 dps for k <= 8 over
q in {1/2, 1, 3/5, 2, 99/100, 7/3}.

Only q^2 enters the kernel and the recurrence, so negative q works unchanged
(arsinh is odd); q = 0 dies on M_0's 1/q — the same natural crash the
hyperbolic family surfaces for its own 1/q.

arcosh_q is not implemented: arcosh(q) = int_1^q dt/sqrt(t^2-1) rescales via
t = 1 + (q-1) x^2 to sqrt(2(q-1)) * int_0^1 dx/sqrt(1 + (q-1) x^2/2), whose
moments L_k span {arcosh(q)/sqrt(2(q-1)), sqrt((q+1)/2), 1}. The endpoint
branch point at t = 1 forces kernel prefactor sqrt(2(q-1)) to normalize the
arcosh coefficient into QQ, which moves the rational atom onto sqrt(2(q-1))
and the sqrt((q+1)/2) atom onto sqrt(q^2-1): the QQ-moment space
{arcosh_q, sqrt(q^2-1), sqrt(2(q-1))} has no "1" direction, so a nonzero
rational bound can never appear in the identity. A gamma-style composite
integrand could carry one, but its sign-definiteness certificate is not the
poly_nonneg machinery — out of clean scope.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import sympy as sp

from ..engine import mn_order, search
from ..moment import Moment, add, combine, scale
from ..render import emit, lhs_tex, rat_tex, wire_fraction

LIMIT = 10

x = sp.symbols("x")


def arsinh_moments(kmax: int, q: Fraction) -> list[Moment]:
    """M_k = int_0^1 x^k / sqrt(1+q^2 x^2) dx for k = 0..kmax.

    Symbols: "arsinh_q" = arsinh(q), "sqrt_1q2" = sqrt(1+q^2), "1".
    """
    q2 = q * q
    out: list[Moment] = [
        {"arsinh_q": Fraction(1) / q},
        {"sqrt_1q2": Fraction(1) / q2, "1": -Fraction(1) / q2},
    ]
    for k in range(2, kmax + 1):
        prev = scale(out[k - 2], Fraction(-(k - 1)) / (k * q2))
        out.append(add({"sqrt_1q2": Fraction(1) / (k * q2)}, prev))
    return out


def basis_moment(moms: list[Moment], m: int, n: int, j: int) -> Moment:
    """int_0^1 x^{m+j} (1-x)^n / sqrt(1+q^2 x^2) dx = sum_i (-1)^i C(n,i) M_{m+i+j}."""
    return combine(
        [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)],
        [moms[m + i + j] for i in range(n + 1)],
    )


def basis(m: int, n: int, q: Fraction) -> list[Moment]:
    """True moments of x^m (1-x)^n {1, x, x^2} / sqrt(1+q^2 x^2) for (a, b, c)."""
    moms = arsinh_moments(m + n + 2, q)
    return [basis_moment(moms, m, n, j) for j in range(3)]


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """Search (m, n) in the author's order; return the parameter/solution dict.

    ``power`` is the argument q != 0. Raises WrongDirection/NoSolution from
    ..engine on failure; ``exact`` is accepted for signature compatibility —
    this family has no site path.
    """
    moms = arsinh_moments(2 * LIMIT + 2, power)
    plans = ((m, n, [basis_moment(moms, m, n, j) for j in range(3)]) for m, n in mn_order(LIMIT))
    sign = 1 if comp == ">" else -1
    target = {"arsinh_q": Fraction(sign), "1": -sign * bound}
    solved = search(plans, target, True)
    return emit(solved.m, solved.n, solved.coeffs)


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """The proof-equation LaTeX for arsinh_q (no site counterpart — clean form).

    The kernel 1/sqrt(1+q^2 x^2) is printed with cleared denominators:
    q = qn/qd gives qd/sqrt(qn^2 x^2 + qd^2).
    """
    m, n = int(params["m"]), int(params["n"])
    au, bu, cu, u = (Fraction(params[k]) for k in ("au_val", "bu_val", "cu_val", "u_val"))
    q = wire_fraction(power)
    qn, qd = q.numerator, q.denominator
    integrand = (
        x**m * (1 - x) ** n * (au + bu * x + cu * x**2) * qd / (u * sp.sqrt(qn**2 * x**2 + qd**2))
    )
    const = rf"\operatorname{{arsinh}}\left({rat_tex(power)}\right)"
    if wire_fraction(bound) < 0:
        # negative bound: fold the sign instead of printing `- \dfrac{-9}{10}`
        mag = rat_tex(-wire_fraction(bound))
        lhs = f"{const} + {mag}" if comp == ">" else f"-{mag} - {const}"
    else:
        lhs = lhs_tex(const, rat_tex(bound), comp)
    return f"{lhs} = \\int_0^1 {sp.latex(integrand)} \\mathrm{{d}} x > 0"
