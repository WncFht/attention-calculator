"""Exact checker for the gauss-erf family: gaussint_q, dawson_q, erfiint_q.

params encode ``P = (au + bu*x + cu*x^2)/u`` multiplying ``x^m (1-x)^n`` and
the family kernel (``e^{−q²x²}`` / ``q·e^{q²(x²−1)}`` / ``e^{q²x²}``) on
[0, 1]. The basis is the kernel's own moments/basis_moment built on |q|;
the claimed vector is s·σ·C(|q|) − s·bound with σ = sign(q), matching the
oddness fold in kernels.gauss_erf.solve. The zero-integrand guard rejects
the vacuous ``0 = ∫ 0 dx > 0``.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.gauss_erf import basis_moment, moments
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val")]
    s = Fraction(1 if comp == ">" else -1)
    sigma = Fraction(1 if power > 0 else -1)
    target = {k: v for k, v in {kind: s * sigma, "1": -s * bound}.items() if v}
    moms = moments(kind, m + n + 2, abs(power))
    basis = [basis_moment(moms, m, n, j) for j in range(3)]
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
