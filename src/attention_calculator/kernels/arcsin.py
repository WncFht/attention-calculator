"""arcsin_q kernel on [0, 1]: claims ``arcsin q ~ bound`` for rational q in (0,1).

Author's type-10 construction: since ``arcsin q = int_0^1 q dx/sqrt(1-q^2 x^2)``,
take the kernel ``K(x) = 1/sqrt(1-q^2 x^2)`` (the overall q folds into the
moment scale). With ``M_k = int_0^1 x^k K(x) dx``, substituting
``x = sin(theta)/q`` gives ``M_k = q^{-k-1} int_0^{arcsin q} sin^k theta dtheta``
and the classical reduction ``int sin^k = -sin^{k-1} cos / k + (k-1)/k int
sin^{k-2}`` yields

    M_0 = arcsin(q)/q,        M_1 = (1 - w)/q^2,
    M_k = (k-1)/(k q^2) * M_{k-2} - w/(k q^2),      w := sqrt(1 - q^2),

so every basis moment lands in ``span{arcsin q, w, 1}`` over QQ(q). ``w`` is a
parasitic constant -- the target forces its coefficient to 0, exactly like
``cos_q`` in the sin_q family. ``P = a + b x + c x^2`` matches the 3-dim space.

Degenerate case: when ``1 - q^2`` is a rational square (``w = rho`` in QQ, e.g.
the 3-4-5 family q = 3/5, 4/5, 5/13, ...), the true moment space collapses to
``span{arcsin q, 1}``. We then fold w's value into the rational component and
solve the 2-dim system with ``P = a + b x``. The folded family is strictly
stronger than insisting on the symbolic 3x3 solve: q = 4/5 is unprovable in
the symbolic system at any reasonable budget but resolves at (4, 8) folded.

Coverage is intrinsically asymmetric in q: the provable bound window
approaches ``arcsin q`` only through (m, n) with m growing roughly like
1/(1 - q), so q >~ 0.8 needs exponents beyond a modest budget. The site ships
``arcsin_q`` commented out in the frontend -- consistent with the author
having hit the same wall.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd, isqrt

import sympy as sp

from ..engine import NoSolution, mn_order, search
from ..moment import Moment, combine
from ..render import emit, lhs_tex, rat_tex, wire_pair

LIMIT = 30

x = sp.symbols("x")


def arcsin_moments(kmax: int, q: Fraction) -> list[Moment]:
    """M_k = int_0^1 x^k/sqrt(1 - q^2 x^2) dx for k = 0..kmax, over
    {arcsin_q, w, 1} with w = sqrt(1 - q^2) symbolic."""
    q2 = q * q
    m: list[Moment] = [
        {"arcsin_q": Fraction(1, q)},
        {"w": -Fraction(1, q2), "1": Fraction(1, q2)},
    ]
    for k in range(2, kmax + 1):
        f = Fraction(k - 1, k) / q2
        m.append(
            {
                "arcsin_q": f * m[k - 2].get("arcsin_q", Fraction(0)),
                "w": f * m[k - 2].get("w", Fraction(0)) - Fraction(1, k * q2),
                "1": f * m[k - 2].get("1", Fraction(0)),
            }
        )
    return [{k: v for k, v in mom.items() if v} for mom in m]


def sqrt_rational(v: Fraction) -> Fraction | None:
    """sqrt(v) as a Fraction when v is a rational square, else None."""
    rn, rd = isqrt(v.numerator), isqrt(v.denominator)
    if rn * rn == v.numerator and rd * rd == v.denominator:
        return Fraction(rn, rd)
    return None


def fold_w(mom: Moment, rho: Fraction) -> Moment:
    """Substitute the rational value rho for the symbolic w in a moment."""
    out = dict(mom)
    w = out.pop("w", Fraction(0))
    if w:
        out["1"] = out.get("1", Fraction(0)) + w * rho
    return {k: v for k, v in out.items() if v}


def basis_moment(m: list[Moment], mexp: int, nexp: int, j: int) -> Moment:
    """int_0^1 x^{mexp+j} (1-x)^nexp K(x) dx = sum_i (-1)^i C(n,i) M_{m+i+j}."""
    return combine(
        [Fraction((-1) ** i * comb(nexp, i)) for i in range(nexp + 1)],
        [m[mexp + i + j] for i in range(nexp + 1)],
    )


def check_input(q: Fraction) -> None:
    """arcsin_q's domain is the open unit interval (site's disabled UI text)."""
    if q <= 0 or q >= 1:
        raise ValueError("arcsin后的值只能在(0,1)内，请输入一个在(0,1)内的分数")


def solve(q: Fraction, comp: str, bound: Fraction):
    """Find (m, n, a, b, c) whose integrand proves ``s*(arcsin q - bound) >= 0``.

    For rational w = sqrt(1 - q^2) the moment space collapses to
    span{arcsin q, 1} and the search uses P = a + b x (two unknowns);
    otherwise the symbolic-w system needs the full quadratic P.
    """
    check_input(q)
    s = Fraction(1 if comp == ">" else -1)
    target = {"arcsin_q": s, "1": -s * bound}
    rho = sqrt_rational(1 - q * q)
    nunk = 2 if rho is not None else 3
    moms = arcsin_moments(2 * LIMIT + 2, q)
    if rho is not None:
        moms = [fold_w(mom, rho) for mom in moms]
    plans = (
        (mexp, nexp, [basis_moment(moms, mexp, nexp, j) for j in range(nunk)])
        for mexp, nexp in mn_order(LIMIT)
    )
    try:
        return search(plans, target, True)
    except NoSolution:
        if rho is None:
            raise
    # folded 2-dim scan exhausted: retry the symbolic-w 3-unknown system --
    # the checker verifies either emission against the true folded moments
    moms = arcsin_moments(2 * LIMIT + 2, q)
    plans = (
        (mexp, nexp, [basis_moment(moms, mexp, nexp, j) for j in range(3)])
        for mexp, nexp in mn_order(LIMIT)
    )
    return search(plans, target, True)


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """Run the proof search for an arcsin_q request (exact-mode-only type)."""
    solved = solve(power, comp, bound)
    result = emit(solved.m, solved.n, solved.coeffs)
    result["type"] = kind
    return result


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """LaTeX for the emitted proof, in the author's single-fraction style.

    The kernel ``1/sqrt(1 - q^2 x^2) = d/sqrt(d^2 - n^2 x^2)`` (q = n/d) lands
    a factor d in the numerator; reducing d against u keeps the denominator a
    plain integer times the radical -- e.g. the article's own example prints
    ``x(1-x)(1710+37x-236x^2) / (360 sqrt(9-x^2))`` for q = 1/3, u = 1080.
    """
    mexp, nexp = params["m"], params["n"]
    au, bu, cu = (Fraction(params[k]) for k in ("au_val", "bu_val", "cu_val"))
    u = Fraction(params["u_val"])
    qn, qd = wire_pair(power)

    const = rf"\arcsin\left({rat_tex(power)}\right)"
    lhs = lhs_tex(const, rat_tex(bound), comp)

    g = gcd(qd, int(u))
    num = sp.Integer(qd // g) * x**mexp * (1 - x) ** nexp * (au + bu * x + cu * x**2)
    den = sp.Integer(int(u) // g) * sp.sqrt(qd * qd - qn * qn * x**2)
    body = rf"\dfrac{{{sp.latex(num)}}}{{{sp.latex(den)}}}"
    return lhs + rf" = \int_0^1 {body} \mathrm{{d}} x > 0"
