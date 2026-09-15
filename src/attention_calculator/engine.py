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


def mn_order(limit: int) -> Iterator[tuple[int, int]]:
    """Yield (m, n) with m+n ascending, then |m-n|, then m ascending.

    Observed site convention: among (m, n) pairs with the same sum and the
    same |m-n|, the smaller-m variant is tried first — e.g. pi>8/3 resolves
    at (0, 1) even though (1, 0) is also a valid proof.
    """
    for s in range(0, 2 * limit + 1):
        pairs = [(m, s - m) for m in range(s + 1) if s - m <= limit and m <= limit]
        pairs.sort(key=lambda p: (abs(p[0] - p[1]), p[0]))
        yield from pairs


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
    a = [r[:] + [b] for r, b in zip(rows, rhs)]
    for col in range(n):
        piv = next((r for r in range(col, n) if a[r][col] != 0), None)
        if piv is None:
            raise ValueError("singular moment system")
        a[col], a[piv] = a[piv], a[col]
        for r in range(n):
            if r != col and a[r][col] != 0:
                f = a[r][col] / a[col][col]
                a[r] = [x - f * y for x, y in zip(a[r], a[col])]
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
    sign: int  # +1 integrand >= 0, -1 integrand <= 0


def search(
    plans: Iterator[tuple[int, int, list[Moment]]],
    target: Moment,
    nonneg: bool,
) -> Solved:
    """Iterate (m, n, basis_moments); return first sign-definite solution.

    If a solution's polynomial factor is uniformly non-positive the identity
    proves the opposite inequality: raise WrongDirection. If it changes sign,
    keep searching.
    """
    for m, n, basis in plans:
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            continue
        check = poly_nonneg if nonneg else poly_nonpos
        if check(coeffs):
            return Solved(m, n, coeffs, +1)
        if check([-c for c in coeffs]):
            raise WrongDirection
    raise NoSolution
