"""psi1_q kernel family on [0, 1]: the trigamma constant psi'(q) vs r.

Exact-mode-only type (the site's dormant type 39 ``psi`` covered only
q in {1/3, 2/3} through a cyclotomic kernel). This family covers all
q in QQ_{>0} through the telescope kernel pair of
docs/2026-09-16-w3-research-trigamma.md:

    '>': K(x) = x^{q-1} (x - 1 - ln x)/(1 - x) >= 0
    '<': K(x) = x^{q+N-2} (1 - x + x ln x)/(1 - x) >= 0,  N = 0 if q > 1 else 1

Both sign lemmas are ln t <= t - 1 (t > 0): x - 1 - ln x is the inequality
at t = x, and 1 - x + x ln x >= 0 is the same inequality at t = 1/x; each
factor vanishes only at the endpoint x = 1. The power moments are
telescoped trigamma tails (recurrence psi'(s+1) = psi'(s) - 1/s^2):

    int x^u K_> dx = psi'(q+u) - 1/(q+u)     = psi'(q) - L_u,
    int x^u K_< dx = 1/(q+N+u-1) - psi'(q+N+u) = U_u - psi'(q),

with L_u = T_u + 1/(q+u), U_u = T_{N+u} + 1/(q+N+u-1),
T_k(q) = sum_{0<=i<k} 1/(q+i)^2. Every moment lives in span{psi'(q), 1};
a (1-x)^n factor cancels the kernel's pole and the psi' row with it
(sum_t (-1)^t C(n,t) = delta_{n0}), so the search is a one-dimensional
m scan at n = 0.

The solved P = a+bx always has a+b = 1, so the site rule (a >= 0 and
a+b >= 0) proves every true bound without a constant tail: '>' hits at the
first m with bound <= L_{m+1} (L_m strictly increases to psi'(q)), '<' at
the first m with bound >= U_{m+1} (U_m strictly decreases to psi'(q));
a larger-than-minimal N only relabels m, since U_u(N+1) = U_{u+1}(N), so
the canonical N is optimal and is not stored in the parameters. A false
claim produces neither a nonneg nor a nonpos P — the scan exhausts to an
honest NoSolution. power carries q itself; the coefficient on psi'(q) is
always 1.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import sympy as sp

from ..engine import search
from ..moment import Moment, combine
from ..render import emit, lhs_tex, rat_tex, wire_pair

# 一维 m 扫描代价极低（每档一个 2x2 有理解），取 256 与 beta 族
# EXACT_LT_LIMIT 同档：q~1 时可证紧界深度 ~1/(2 m^2)，覆盖到 ~1e-5
LIMIT = 256

x = sp.symbols("x")


def shift_n(q: Fraction) -> int:
    """Least N >= 0 with q+N > 1: the '<' kernel needs an integrable
    exponent q+N-2 > -1, and the tail bound psi'(s) < 1/(s-1) needs s > 1."""
    return 0 if q > 1 else 1


def tail_sums(kmax: int, q: Fraction) -> list[Fraction]:
    """T_k(q) = sum_{0<=i<k} 1/(q+i)^2 for k = 0..kmax, built incrementally."""
    out = [Fraction(0)]
    for k in range(kmax):
        out.append(out[-1] + Fraction(1, (q + k) ** 2))
    return out


def moments(kind: str, comp: str, kmax: int, q: Fraction) -> list[Moment]:
    """Power moments int_0^1 x^u K_dir(x) dx over span{psi1_q, 1}, u = 0..kmax.

    '>' kernel: {psi1_q: 1, "1": -L_u}. '<' kernel: {psi1_q: -1, "1": U_u}.
    """
    if comp == ">":
        ts = tail_sums(kmax, q)
        return [{kind: Fraction(1), "1": -(ts[u] + Fraction(1, q + u))} for u in range(kmax + 1)]
    n = shift_n(q)
    ts = tail_sums(kmax + n, q)
    return [
        {kind: Fraction(-1), "1": ts[n + u] + Fraction(1, q + n + u - 1)} for u in range(kmax + 1)
    ]


def basis_moment(moms: list[Moment], m: int, n: int, j: int) -> Moment:
    """int_0^1 x^{m+j} (1-x)^n K(x) dx = sum_i (-1)^i C(n,i) M_{m+j+i}.

    For n >= 1 the psi' row cancels (sum_i (-1)^i C(n,i) = 0): the basis is
    purely rational and never matches a target carrying psi'(q), which is
    why the emitted search scans m at n = 0.
    """
    return combine(
        [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)],
        [moms[m + i + j] for i in range(n + 1)],
    )


def check_input(q: Fraction) -> None:
    """psi' has poles at q = 0, -1, -2, ...; negative non-integers are a
    documented extension (recurrence-shift into (0,1) first) this omits."""
    if q <= 0:
        raise ValueError("psi1_q 仅接受 q>0：q<=0 是 ψ′ 的极点或域外参数")


def solve(kind: str, q: Fraction, comp: str, bound: Fraction):
    """Find (m, a, b) at n = 0 whose integrand proves s*(psi'(q) - bound) >= 0.

    Target is s*psi'(q) - s*bound over span{psi1_q, 1} with s = +1 for '>'.
    """
    check_input(q)
    s = Fraction(1 if comp == ">" else -1)
    target = {k: v for k, v in {kind: s, "1": -s * bound}.items() if v}
    moms = moments(kind, comp, LIMIT + 1, q)
    plans = ((m, 0, [basis_moment(moms, m, 0, j) for j in range(2)]) for m in range(LIMIT + 1))
    return search(plans, target, True)


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("psi1_q", q, comp, bound) -> /calculate-shaped dict.

    Exact-mode-only type; ``exact`` is accepted for signature compatibility
    (there is no site path to reproduce).
    """
    solved = solve(kind, power, comp, bound)
    result = emit(solved.m, solved.n, solved.coeffs)
    result["type"] = kind
    return result


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX: ``psi'(q) - r``/``r - psi'(q)`` = int x^m(1-x)^n P K.

    The kernel factor prints as the merged power x^{m+q-1} or
    x^{m+q+N-2} times the telescope bracket over (1-x), matching the
    direction kernel the params were solved under.
    """
    m, n = int(params["m"]), int(params["n"])
    au, bu, u = (int(params[k]) for k in ("au_val", "bu_val", "u_val"))
    qn, qd = wire_pair(power)
    q = sp.Rational(qn, qd)

    const = rf"\psi_1\left({rat_tex(power)}\right)"
    lhs = lhs_tex(const, rat_tex(bound), comp)

    if comp == ">":
        kern = x ** (q - 1) * (x - 1 - sp.log(x)) / (1 - x)
    else:
        kern = x ** (q + shift_n(Fraction(qn, qd)) - 2) * (1 - x + x * sp.log(x)) / (1 - x)
    body = x**m * (1 - x) ** n * (au + bu * x) * kern / u
    return f"{lhs} = \\int_0^1 {sp.latex(body)} \\mathrm{{d}} x > 0"
