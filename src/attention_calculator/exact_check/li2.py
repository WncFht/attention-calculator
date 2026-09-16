"""Exact checker for li2_q (exact-mode only).

params encode P = (au + bu*x + cu*x^2)/u multiplying (1-x)^n * sigma*(phi - L)
on [0, 1], phi = -ln(1-qx)/x; ``fam``/``d`` select the shift polynomial L and
double as the sign-lemma key — kernels.li2.certified tabulates which
(fam, d, comp, q) combinations make K pointwise nonnegative (termwise tail,
alternating bound, or chord/tangent vs phi's convexity).  The basis moments
are the kernel's own over {li2_q, ln_1mq, 1}; the claimed vector is
s*Li_2(q) - s*bound.  m must be 0 — the Li_2 symbol lives only in A_0, so any
x^m factor kills the constant component and the identity cannot hold anyway.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.li2 import LI2, basis_moment, certified, moments, shift_poly
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val")]
    fam, d = str(params["fam"]), int(params["d"])
    s = 1 if comp == ">" else -1
    target = {k: v for k, v in {LI2: Fraction(s), "1": -s * bound}.items() if v}
    moms = moments(m + n + 2, power)
    basis = [basis_moment(moms, shift_poly(fam, d, power), m, n, j, s) for j in range(3)]
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": certified(fam, d, comp, power) and poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
