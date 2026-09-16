"""gauss-erf kernel family on [0, 1]: gaussint_q, dawson_q, erfiint_q.

Exact-mode-only types — literal ``erf(q) ⋚ r`` is unreachable in this moment
machine (every moment's erf component carries a √π the span cannot absorb;
docs/2026-09-16-w3-research-erf.md §√π障碍), so the family proves the three
scaled constants the same machine does reach:

    gaussint_q: G(q) = ∫₀^q e^{−t²}dt = √π erf(q)/2,  kernel e^{−q²x²}
    dawson_q:   F(q) = e^{−q²}∫₀^q e^{t²}dt,          kernel q·e^{q²(x²−1)}
    erfiint_q:  H(q) = ∫₀^q e^{t²}dt = √π erfi(q)/2,  kernel e^{q²x²}

All three moments come from the same IBP on x^{k−1}e^{±q²x²}: the boundary
term lands on the parasite symbol e^{±q²}, so each chain is two-step,

    gaussint_q: I_0 = G/q, I_1 = t(1 − e^{−q²}), I_k = (k−1)t·I_{k−2} − t·e^{−q²}
    dawson_q:   J'_0 = F,  J'_1 = qt(1 − e^{−q²}), J'_k = qt − (k−1)t·J'_{k−2}
    erfiint_q:  J_0 = H/q, J_1 = t(e^{q²} − 1),    J_k = t·e^{q²} − (k−1)t·J_{k−2}

with t = 1/(2q²); dawson's normalized kernel is what puts F itself (not
e^{q²}F/q) in the span with a rational coefficient. Each moment space is
span{C_q, e^{±q²}, 1} over QQ(q), so P = a + bx + cx² carries the three
unknowns and the parasite cancels to exactly zero in any solved proof.

All three constants are odd in q; the moment machine only sees |q|, so a
negative argument folds into the target sign — ``C(q) ⋚ r`` is proved as
``σ·C(|q|)`` vs ``r`` with σ = sign(q), which displays back as the literal
claim since σ·C(|q|) = C(q). q = 0 degenerates every kernel (I_0 ~ C(q)/q,
the dawson prefactor vanishes) and is rejected.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import sympy as sp

from ..engine import mn_order, search
from ..moment import Moment, add, combine, scale
from ..render import emit, lhs_tex, rat_tex, wire_pair

LIMIT = 10

x = sp.symbols("x")

# LHS constant names in the rendered proof — no site counterparts exist.
NAME = {
    "gaussint_q": "GaussInt",
    "dawson_q": "Dawson",
    "erfiint_q": "ErfiInt",
}


def moments(kind: str, kmax: int, q: Fraction) -> list[Moment]:
    """Power moments I_k = ∫_0^1 x^k K(x) dx for k = 0..kmax, q > 0.

    The moment dicts use the kind name as the constant symbol and
    "e^{-q^2}"/"e^{q^2}" for the parasite exponential.
    """
    t = Fraction(1, 2 * q * q)
    if kind == "gaussint_q":
        eps, c = 1, {"e^{-q^2}": -t}
        out = [{kind: Fraction(1, q)}, {"1": t, "e^{-q^2}": -t}]
    elif kind == "dawson_q":
        eps, c = -1, {"1": q * t}
        out = [{kind: Fraction(1)}, {"1": q * t, "e^{-q^2}": -q * t}]
    else:  # erfiint_q
        eps, c = -1, {"e^{q^2}": t}
        out = [{kind: Fraction(1, q)}, {"e^{q^2}": t, "1": -t}]
    for k in range(2, kmax + 1):
        out.append(add(scale(out[k - 2], eps * (k - 1) * t), c))
    return out


def basis_moment(moms: list[Moment], m: int, n: int, j: int) -> Moment:
    """∫_0^1 x^{m+j} (1-x)^n K(x) dx = sum_i (-1)^i C(n,i) I_{m+j+i}."""
    return combine(
        [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)],
        [moms[m + i + j] for i in range(n + 1)],
    )


def check_input(q: Fraction) -> None:
    """q = 0 collapses every kernel of the family: I_0 = C(q)/q is undefined,
    the recurrences carry 1/q², and the dawson kernel itself vanishes."""
    if q == 0:
        raise ValueError("q不能为0：q=0 时该积分常数退化")


def solve(kind: str, q: Fraction, comp: str, bound: Fraction):
    """Find (m, n, a, b, c) whose integrand proves s·(C(q) − bound) >= 0.

    The target is s·σ·C(|q|) − s·bound over span{C_q, e^{±q²}, 1}: σ = sign(q)
    folds the oddness of C so the |q| moment machine certifies the literal
    claim for negative q as well.
    """
    check_input(q)
    s = Fraction(1 if comp == ">" else -1)
    sigma = Fraction(1 if q > 0 else -1)
    target = {k: v for k, v in {kind: s * sigma, "1": -s * bound}.items() if v}
    moms = moments(kind, 2 * LIMIT + 2, abs(q))
    plans = ((m, n, [basis_moment(moms, m, n, j) for j in range(3)]) for m, n in mn_order(LIMIT))
    return search(plans, target, True)


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("gaussint_q", q, comp, bound) -> /calculate-shaped dict.

    Exact-mode-only family; ``exact`` is accepted for signature compatibility
    (the math is identical — there is no site path to reproduce).
    """
    solved = solve(kind, power, comp, bound)
    result = emit(solved.m, solved.n, solved.coeffs)
    result["type"] = kind
    return result


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX: ``C(q) − r``/``r − C(q)`` = ∫ x^m(1−x)^n P(x) K(x).

    The kernel factor prints ``e^{±q²x²}`` / ``q·e^{q²(x²−1)}`` with q the
    wire value's |q| — q² is sign-blind and dawson's prefactor needs |q|
    since the moments are computed on |q|.
    """
    m, n = int(params["m"]), int(params["n"])
    au, bu, cu = (int(params[k]) for k in ("au_val", "bu_val", "cu_val"))
    u = int(params["u_val"])
    qn, qd = wire_pair(power)
    q2 = sp.Rational(qn * qn, qd * qd)

    const = rf"\mathrm{{{NAME[kind]}}}\left({rat_tex(power)}\right)"
    lhs = lhs_tex(const, rat_tex(bound), comp)

    if kind == "gaussint_q":
        kern = sp.exp(-q2 * x**2)
    elif kind == "dawson_q":
        kern = sp.Rational(abs(qn), qd) * sp.exp(q2 * (x**2 - 1))
    else:
        kern = sp.exp(q2 * x**2)
    body = x**m * (1 - x) ** n * (au + bu * x + cu * x**2) * kern / u
    return f"{lhs} = \\int_0^1 {sp.latex(body)} \\mathrm{{d}} x > 0"
