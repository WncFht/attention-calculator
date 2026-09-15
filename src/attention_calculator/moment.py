"""Exact moment vectors: integrals of basis functions times kernels, over QQ.

Every basis integral evaluates to a finite QQ-linear combination of
transcendental constants plus a rational part. A moment is stored as a plain
dict mapping a constant symbol to a fractions.Fraction coefficient, e.g.

    {"1": Fraction(3, 4), "pi": Fraction(-1, 4), "ln2": Fraction(1, 2)}

means the integral equals 3/4 - pi/4 + (ln 2)/2.

Constant symbols are plain strings naming the constant ("pi", "e", "ln2",
"sin_q", "catalan", ...). Which symbols can appear for a given kernel family
is fixed by docs/kernel-spec.md; coefficients are always exact rationals.
"""

from fractions import Fraction

Moment = dict[str, Fraction]


def add(u: Moment, v: Moment) -> Moment:
    """Sum two moment vectors."""
    out = dict(u)
    for k, c in v.items():
        out[k] = out.get(k, Fraction(0)) + c
    return {k: c for k, c in out.items() if c != 0}


def scale(u: Moment, c: Fraction) -> Moment:
    """Multiply a moment vector by a rational scalar."""
    return {k: c * v for k, v in u.items() if c * v != 0}


def combine(coeffs: list[Fraction], moments: list[Moment]) -> Moment:
    """Linear combination of moment vectors with rational coefficients."""
    out: Moment = {}
    for c, m in zip(coeffs, moments, strict=True):
        out = add(out, scale(m, c))
    return out
