"""Exact checker for arcsin_q: claims ``arcsin q ~ bound`` for q in (0,1).

params encode ``P = (au + bu*x + cu*x^2)/u`` multiplying ``x^m (1-x)^n`` and
the kernel ``1/sqrt(1 - q^2 x^2)`` on [0, 1]. The basis is the kernel's own
arcsin_moments/basis_moment; when ``1 - q^2`` is a rational square the
checker reproduces the kernel's fold of w = sqrt(1 - q^2) into the rational
component, so the 2-unknown emission (cu_val = 0) verifies against the
collapsed {arcsin_q, 1} space.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.arcsin import arcsin_moments, basis_moment, fold_w, sqrt_rational
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val")]
    s = Fraction(1 if comp == ">" else -1)
    target = {"arcsin_q": s, "1": -s * bound}
    target = {k: v for k, v in target.items() if v}  # Moments omit zero coeffs
    q = power
    moms = arcsin_moments(m + n + 2, q)
    rho = sqrt_rational(1 - q * q)
    if rho is not None:
        moms = [fold_w(mom, rho) for mom in moms]
    basis = [basis_moment(moms, m, n, j) for j in range(3)]
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
