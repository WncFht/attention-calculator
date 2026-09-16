"""Exact checker for pi_sqrt2.

params encode ``P = (au + bu·x⁴)/u`` on ``x^{4m+res}(1-x)(1-x⁴)^{-c}`` —
res=2, c=3/4 for '>' and res=3, c=1/4 for '<' (the emitted n is a
direction flag, not read here — comp selects the kernel). The basis
moments are the kernel's own (kernels.pi_sqrt2.basis); the claimed
vector is ``q·S - r`` for '>' resp. ``r - q·S`` for '<' over
{pi_sqrt2, 1} (kernels.pi_sqrt2.target drops the zero keys). nonneg is
the even-kernel rule a >= 0 and a+b >= 0 plus the zero-integrand guard
against the vacuous ``0 = ∫ 0 dx > 0``.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.pi_sqrt2 import basis, target
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m = int(params["m"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params["au_val"]) / u, Fraction(params["bu_val"]) / u]
    tgt = target(comp, power, bound)
    integrand = combine(coeffs, [basis(comp, m, i) for i in (0, 1)])
    return {
        "identity_ok": integrand == tgt,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": tgt,
    }
