"""li2_q kernel on [0, 1]: ``Li_2(q) ~ bound`` for rational q < 1, q != 0.

Exact-mode-only type (docs/2026-09-16-w3-research-li2.md).  The constant is
``Li_2(q) = int_0^1 phi`` with ``phi(x) = -ln(1-qx)/x = sum_{i>=1} q^i x^{i-1}/i``
(value at x = 0 is the limit q).  Power moments close on
span{Li_2(q), ln(1-q), 1}:

    A_0 = Li_2(q);   A_k = (q^{-k}-1)/k * ln(1-q) + q^{-k}/k * sum_{l=1..k} q^l/l.

Li_2 lives only in A_0, which forces the basis to ``(1-x)^n`` (m == 0 is a
structural corollary, not a search axis) and pins P's constant term to +1:
the comparison sign is absorbed by the kernel ``K = sigma*(phi - L)`` with
sigma = +1 for '>' / -1 for '<', never by flipping P.

The shift polynomial L must be a certified sub-function ('>') resp.
super-function ('<') of phi; ``certified`` tabulates the legal
(family, direction, q-range) pairings — that table IS the sign argument:

    fam      L(x)                                      legal for
    taylor   p_d = sum_{i=1..d} q^i x^{i-1}/i           '>' q>0: tail termwise >= 0
    lift     p_d + q^{d+1}/((d+1)(1-q)) * x^d           '<' q>0: tail <= x^d q^{d+1}/((d+1)(1-q))
    alt      p_d, d odd for '>' / even for '<'          -1 <= q < 0: alternating bound
    chord    q + (-ln(1-q)-q) x                         '>' q<0 or '<' q>0: chord vs phi's convexity
    tan0     q + q^2 x/2  (= p_2)                       '<' q<0: concave tangent (covers q < -1)

phi'' is positive termwise for q > 0 (strictly convex); for q < 0, phi is
strictly concave — x^3 phi''(x) = B(-qx) with B(y) = y^2/(1+y)^2 +
2y/(1+y) - 2 ln(1+y), B(0) = 0 and B'(y) = -2 y^2/(1+y)^3 < 0, so B < 0 for
all u = -q > 0.  For -1 <= q < 0 the alternating-series families supersede
the research doc's lifted shapes (chord + t*x(1-x), tan0 + t*x^2): their
optimal t*(q) is a transcendental univariate minimum with no exact sign
certificate, so they are not implemented.  q < -1 (where the power series
diverges at x > 1/u) keeps only the linear chord/tan0 witnesses — coverage
there is limited to what concavity alone certifies.

Search: n in [0, N_LIMIT] outer, (fam, d) inner, d <= D_LIMIT.  q -> 1^- is
the known weak spot — the gap shrinks like q^d, so q = 99/100 needs d in the
hundreds; honest NoSolution beyond the budget.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import sympy as sp

from ..engine import NoSolution, Solved, WrongDirection, poly_nonneg, poly_nonpos, solve_moment
from ..moment import Moment, add, scale
from ..render import emit, lhs_tex, rat_tex, wire_pair

N_LIMIT = 16
D_LIMIT = 64

x = sp.symbols("x")

# moment symbols over span{Li_2(q), ln(1-q), 1}
LI2 = "li2_q"
LN = "ln_1mq"


def moments(kmax: int, q: Fraction) -> list[Moment]:
    """A_k = int_0^1 x^k phi(x) dx for k = 0..kmax over {li2_q, ln_1mq, 1}."""
    out: list[Moment] = [{LI2: Fraction(1)}]
    for k in range(1, kmax + 1):
        out.append(
            {
                LN: (q**-k - 1) / k,
                "1": q**-k / k * sum(q**e / e for e in range(1, k + 1)),
            }
        )
    return out


def shift_poly(fam: str, d: int, q: Fraction) -> list[Moment]:
    """L's coefficients as Moments indexed by degree (QQ + QQ*ln(1-q) valued)."""
    if fam == "chord":
        return [{"1": q}, {"1": -q, LN: Fraction(-1)}]
    L = [{"1": q**i / i} for i in range(1, d + 1)]
    if fam == "lift":
        L.append({"1": q ** (d + 1) / ((d + 1) * (1 - q))})
    return L


def basis_moment(
    moms: list[Moment], shift: list[Moment], m: int, n: int, j: int, sigma: int
) -> Moment:
    """int_0^1 x^{m+j} (1-x)^n sigma(phi-L) dx
    = sigma sum_t (-1)^t C(n,t) [A_{m+j+t} - sum_i L_i/(m+j+t+i+1)]."""
    out: Moment = {}
    for t in range(n + 1):
        row = dict(moms[m + j + t])
        for i, li in enumerate(shift):
            for sym, v in li.items():
                row[sym] = row.get(sym, Fraction(0)) - v / (m + j + t + i + 1)
        out = add(out, scale({k: v for k, v in row.items() if v}, sigma * (-1) ** t * comb(n, t)))
    return out


def certified(fam: str, d: int, comp: str, q: Fraction) -> bool:
    """Whether (fam, d) is a certified-sign L for this (comp, q) — the kernel
    factor K = sigma(phi - L) is then pointwise >= 0 by the family's lemma.

    taylor: '>' q>0, tail of phi - p_d termwise nonnegative.
    lift:   '<' q>0, tail <= x^d q^{d+1}/((d+1)(1-q)) (x^{i-1} <= x^d and
            1/i <= 1/(d+1) for i > d, then the geometric sum in q).
    alt:    -1 <= q < 0, alternating series with decreasing terms: p_d lies
            below phi for odd d ('>'), above for even d ('<'; d = 0 is the
            bare kernel -phi > 0).
    chord:  '>' q<0 (concave phi above its chord) or '<' q>0 (convex phi
            below its chord); K = |phi - chord| vanishes at both endpoints.
    tan0:   '<' q<0, concave phi below its 0-tangent p_2 — the only '<'
            witness that survives q < -1, where the series diverges.
    """
    if fam == "taylor":
        return comp == ">" and 0 < q < 1 and d >= 0
    if fam == "lift":
        return comp == "<" and 0 < q < 1 and d >= 0
    if fam == "alt":
        return -1 <= q < 0 and d >= 0 and d % 2 == (1 if comp == ">" else 0)
    if fam == "chord":
        return d == 0 and (comp == ">") == (q < 0)
    if fam == "tan0":
        return comp == "<" and q < 0 and d == 2
    return False


def families(q: Fraction, comp: str) -> list[tuple[str, int]]:
    """(fam, d) plan list for this quadrant, in search order."""
    if comp == ">":
        if q > 0:
            return [("taylor", d) for d in range(D_LIMIT + 1)]
        fams = [("chord", 0)]
        if q >= -1:
            fams += [("alt", d) for d in range(1, D_LIMIT + 1, 2)]
        return fams
    if q > 0:
        return [("chord", 0)] + [("lift", d) for d in range(D_LIMIT + 1)]
    if q >= -1:
        return [("alt", d) for d in range(0, D_LIMIT + 1, 2)]
    return [("tan0", 2)]


def check_input(q: Fraction) -> None:
    """q = 0 collapses phi to the zero kernel; q >= 1 puts ln(1-qx)'s branch
    point at or inside the domain (q = 1 also collapses the ln column)."""
    if q == 0:
        raise ValueError("q不能为0：q=0 时核 -ln(1-qx)/x 恒为0")
    if q >= 1:
        raise ValueError("Li2的自变量必须小于1")


def solve(q: Fraction, comp: str, bound: Fraction):
    """Find (n, fam, d, a, b, c) whose integrand proves s*(Li_2(q) - bound) >= 0.

    a solves to +1 identically (the Li_2 column lives only in the j = 0 row);
    a uniformly non-positive solved P certifies the opposite inequality, so
    WrongDirection propagates immediately like engine.search.
    """
    check_input(q)
    s = 1 if comp == ">" else -1
    target = {k: v for k, v in {LI2: Fraction(s), "1": -s * bound}.items() if v}
    moms = moments(N_LIMIT + 2, q)
    fams = families(q, comp)
    for n in range(N_LIMIT + 1):
        for fam, d in fams:
            basis = [basis_moment(moms, shift_poly(fam, d, q), 0, n, j, s) for j in range(3)]
            try:
                coeffs = solve_moment(basis, target)
            except ValueError:
                continue
            if poly_nonneg(coeffs):
                return Solved(0, n, coeffs), fam, d
            if poly_nonpos(coeffs):
                raise WrongDirection
    raise NoSolution


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("li2_q", q, comp, bound) -> /calculate-shaped dict.

    Exact-mode-only type; ``exact`` is accepted for signature compatibility.
    Extra parameters beyond the site's nine: ``fam`` (shift-family tag, the
    checker's sign-lemma key) and ``d`` (family degree parameter).
    """
    solved, fam, d = solve(power, comp, bound)
    result = emit(solved.m, solved.n, solved.coeffs)
    result["parameters"]["fam"] = fam
    result["parameters"]["d"] = d
    result["type"] = kind
    return result


def _shift_sympy(fam: str, d: int, q: sp.Rational) -> sp.Expr:
    """L(x) as a sympy expression, mirroring shift_poly."""
    if fam == "chord":
        return q + (-sp.log(1 - q) - q) * x
    L = sum(q**i / i * x ** (i - 1) for i in range(1, d + 1))
    if fam == "lift":
        L += q ** (d + 1) / ((d + 1) * (1 - q)) * x**d
    return sp.expand(L)


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX: ``Li_2(q) - r``/``r - Li_2(q)`` =
    int (1-x)^n P(x) sigma(phi - L) dx > 0, phi = -ln(1-qx)/x."""
    m, n = int(params["m"]), int(params["n"])
    au, bu, cu = (int(params[k]) for k in ("au_val", "bu_val", "cu_val"))
    u = int(params["u_val"])
    fam, d = str(params["fam"]), int(params["d"])
    qn, qd = wire_pair(power)
    q = sp.Rational(qn, qd)
    s = 1 if comp == ">" else -1

    const = rf"\mathrm{{Li}}_2\left({rat_tex(power)}\right)"
    lhs = lhs_tex(const, rat_tex(bound), comp)

    kern = s * (-sp.log(1 - q * x) / x - _shift_sympy(fam, d, q))
    body = x**m * (1 - x) ** n * (au + bu * x + cu * x**2) * kern / u
    return f"{lhs} = \\int_0^1 {sp.latex(body)} \\mathrm{{d}} x > 0"
