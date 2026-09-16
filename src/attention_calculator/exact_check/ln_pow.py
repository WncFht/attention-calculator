"""Exact checker for ln_q_cube / ln_q_quad.

params encode P = (au + bu*x + cu*x^2 + du*x^3 [+ eu*x^4])/u — the
``du_val``/``eu_val`` fields extend the nine-slot site schema for the
higher coefficients — over basis x^m(1-x)^n with denominator (1+cx)^s,
s = max(m, n, 1), c = q - 1. The true moments live in
span{ln^r q, ..., ln q, 1} (kernels.ln_pow.mu, r = 3 or 4), and the
printed claim s*(ln^r q - bound) is exactly the {ln^r: s, 1: -s*bound}
vector the kernel solved.
"""

from fractions import Fraction

from ..kernels.ln_pow import basis_moment, cubic_nonneg, quartic_nonneg
from ..moment import combine

# kind -> (kernel log power, P coefficient keys, nonneg rule)
_CFG = {
    "ln_q_cube": (2, ("au_val", "bu_val", "cu_val", "du_val"), cubic_nonneg),
    "ln_q_quad": (3, ("au_val", "bu_val", "cu_val", "du_val", "eu_val"), quartic_nonneg),
}


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    lp, keys, nonneg = _CFG[kind]
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in keys]
    s = max(m, n, 1)
    basis = [basis_moment(power - 1, s, m, n, j, lp) for j in range(lp + 2)]
    sign = Fraction(1 if comp == ">" else -1)
    target = {f"ln{lp + 1}": sign, "1": -sign * bound}
    target = {k: v for k, v in target.items() if v}  # Moments omit zero coeffs
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": nonneg(coeffs),
        "integrand": integrand,
        "target": target,
    }
