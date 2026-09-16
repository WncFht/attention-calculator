"""Shared machinery for the integral-proof search.

A kernel family module implements ``prove(q, comp, bound, coef)`` which loops
over ``mn_order`` candidates, builds basis moments, calls ``solve_moment`` and
checks sign via ``poly_nonneg``/``poly_nonpos``. The engine only provides the
exact-rational primitives; family-specific assembly stays in the family module
(see docs/kernel-spec.md).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from fractions import Fraction

from .moment import Moment


class WrongDirection(Exception):
    """The identity exists but its integrand is uniformly non-positive."""


class NoSolution(Exception):
    """No candidate within the exponent budget gave a sign-definite proof."""


class InternalError(Exception):
    """Reproduced site-side crash: maps to the site's generic 500.

    Raised for inputs where zhuyidao.net itself errors out (probed behavior,
    e.g. ln_q_square with q in {5, 7} always 500s) rather than for our own
    bugs.
    """


class EqualClaim(ValueError):
    """Bound equals the constant exactly; the site reports 404 '二者相等'.

    Subclasses ValueError so kernel-level callers see the validation-style
    contract, while the server maps it to 404 ahead of the 400 catch.
    """


# 搜索预算: e、pi 两类型指数上限 30, 其余 10 (见 docs/kernel-spec.md 搜索顺序)
# arcsin_q 是 exact-only 型，其核 LIMIT=30——q→1 时可证界要求 m ~ O(1/(1−q))
# pi_sqrt2 是单轴 m 扫描，核 LIMIT=256（沿用 beta 族 exact '<' 档）；
# dixon 族 m 轴 LIMIT=512；li2_q 的 (1-x)^n 轴 N_LIMIT=16；psi1_q 一维扫 256
EXPONENT_LIMIT = {
    "pi": 30,
    "e": 30,
    "arcsin_q": 30,
    "pi_sqrt2": 256,
    "pi3": 512,
    "pi3_u": 512,
    "pi3_a": 512,
    "li2_q": 16,
    "psi1_q": 256,
}


def mn_order(limit: int) -> Iterator[tuple[int, int]]:
    """Yield (m, n) with m+n ascending, |m-n| ascending, then the parity rule.

    Observed site convention: mirror pairs are tried smaller-m first when the
    sum is odd, larger-m first when even — ln²(3/2)<17/100 resolves at (2, 0)
    though (0, 2) is also a valid proof, while pi>8/3 resolves at (0, 1) over
    (1, 0). Equivalent to walking n = ⌈s/2⌉, ⌈s/2⌉-1, ⌈s/2⌉+1, ⌈s/2⌉-2, …
    """
    for s in range(0, 2 * limit + 1):
        for d in range(s % 2, s + 1, 2):
            lo, hi = (s - d) // 2, (s + d) // 2
            pair = (lo, hi) if s % 2 else (hi, lo)
            for m, n in dict.fromkeys((pair, pair[::-1])):
                if m <= limit and n <= limit:
                    yield m, n


def solve_moment(basis: list[Moment], target: Moment) -> list[Fraction]:
    """Solve sum_j u_j * basis[j] == target for u over QQ, uniquely.

    ``basis[j]`` is the moment of the j-th coefficient function of the
    undetermined polynomial factor. Raises ValueError if the linear system has
    no unique solution.
    """
    keys = sorted(set(target) | {k for m in basis for k in m})
    rows = [[m.get(k, Fraction(0)) for m in basis] for k in keys]
    rhs = [target.get(k, Fraction(0)) for k in keys]
    return gauss_solve(rows, rhs)


def gauss_solve(rows: list[list[Fraction]], rhs: list[Fraction]) -> list[Fraction]:
    """Gaussian elimination over Fraction; requires a unique solution."""
    n = len(rhs)
    if len(rows) != n or any(len(r) != n for r in rows):
        raise ValueError("moment system is not square")
    a = [[*r, b] for r, b in zip(rows, rhs, strict=True)]
    for col in range(n):
        piv = next((r for r in range(col, n) if a[r][col] != 0), None)
        if piv is None:
            raise ValueError("singular moment system")
        a[col], a[piv] = a[piv], a[col]
        for r in range(n):
            if r != col and a[r][col] != 0:
                f = a[r][col] / a[col][col]
                a[r] = [x - f * y for x, y in zip(a[r], a[col], strict=True)]
    return [a[i][n] / a[i][i] for i in range(n)]


def poly_nonneg(coeffs: list[Fraction]) -> bool:
    """Sign test for a+bx or a+bx+cx^2 being >= 0 on [0,1] (author's rules)."""
    a = coeffs[0]
    b = coeffs[1] if len(coeffs) > 1 else Fraction(0)
    c = coeffs[2] if len(coeffs) > 2 else Fraction(0)
    if c <= 0:
        return a >= 0 and a + b + c >= 0
    if b < 0 < b + 2 * c:  # vertex inside (0,1)
        return 4 * a * c - b * b >= 0
    return a >= 0 and a + b + c >= 0


def poly_nonpos(coeffs: list[Fraction]) -> bool:
    """Sign test for the polynomial being <= 0 on [0,1]."""
    return poly_nonneg([-c for c in coeffs])


@dataclass
class Solved:
    """A solved candidate: exponent pair and polynomial coefficients."""

    m: int
    n: int
    coeffs: list[Fraction]


def search(
    plans: Iterator[tuple[int, int, list[Moment]]],
    target: Moment,
    nonneg: bool,
    defer: bool = False,
) -> Solved:
    """Iterate (m, n, basis_moments); return first sign-definite solution.

    If a solution's polynomial factor is uniformly non-positive the identity
    proves the opposite inequality: raise WrongDirection, or with ``defer``
    keep searching and raise only when no nonneg candidate appears (the
    zhuyidao trig_pi scan skips non-positive plans -- a later plan can still
    produce the proof, e.g. past its corrupted (1, 8) formula). If it changes
    sign, keep searching.
    """
    saw_nonpos = False
    for m, n, basis in plans:
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            continue
        check = poly_nonneg if nonneg else poly_nonpos
        if check(coeffs):
            return Solved(m, n, coeffs)
        if check([-c for c in coeffs]):
            if not defer:
                raise WrongDirection
            saw_nonpos = True
    if saw_nonpos:
        raise WrongDirection
    raise NoSolution
