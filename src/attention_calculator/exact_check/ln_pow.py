"""Exact checker for ln_q_cube.

params encode P = (au + bu*x + cu*x^2 + du*x^3)/u — ``du_val`` extends the
nine-slot site schema for the cubic's fourth coefficient — over basis
x^m(1-x)^n with denominator (1+cx)^s, s = max(m, n, 1), c = q - 1. The true
moments live in span{ln3, ln2, ln, 1} (kernels.ln_pow.mu), and the printed
claim s*(ln^3 q - bound) is exactly the {"ln3": s, "1": -s*bound} vector
the kernel solved.
"""

from fractions import Fraction

from ..kernels.ln_pow import basis_moment, cubic_nonneg
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val", "du_val")]
    s = max(m, n, 1)
    basis = [basis_moment(power - 1, s, m, n, j) for j in range(4)]
    sign = Fraction(1 if comp == ">" else -1)
    target = {"ln3": sign, "1": -sign * bound}
    target = {k: v for k, v in target.items() if v}  # Moments omit zero coeffs
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": cubic_nonneg(coeffs),
        "integrand": integrand,
        "target": target,
    }
