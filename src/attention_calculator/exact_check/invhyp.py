"""Exact checker for the invhyp kernel: arsinh_q (exact-mode only).

params encode P = (au + bu*x + cu*x**2)/u multiplying
x^m (1-x)^n / sqrt(1+q^2 x^2) on [0, 1]; the true basis moments are the
kernel's own arsinh_moments/basis_moment over {arsinh_q, sqrt_1q2, 1} —
this family has no site counterpart, so solve and check share it.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.invhyp import basis
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val")]
    sign = 1 if comp == ">" else -1
    raw = {"arsinh_q": Fraction(sign), "1": -sign * bound}
    target = {k: v for k, v in raw.items() if v}
    integrand = combine(coeffs, basis(m, n, power))
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
