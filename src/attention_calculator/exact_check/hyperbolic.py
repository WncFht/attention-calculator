"""Exact checker for hyperbolic kernels: sinh_q, cosh_q, tanh_q, coth_q.

params encode P = (au + bu*x + cu*x^2)/u over basis x^m(1-x)^n with kernel
sinh(qx); moments live in span{sinh q, cosh q, 1} via the kernel's own
basis_moment. This family has no reproduced site bias, so solve and check
share the basis and target builders: tanh_q/coth_q print an overall
1/cosh q (1/sinh q) factor, so the claim tanh q ⋚ bound is solved and
checked as sinh q − bound·cosh q ⋚ 0 (resp. cosh q − bound·sinh q ⋚ 0)
in that space — exactly what target_for returns.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.hyperbolic import basis_moment, target_for
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = params["m"], params["n"]
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val")]
    basis = [basis_moment(m, n, j, power) for j in (0, 1, 2)]
    integrand = combine(coeffs, basis)
    # target_for keeps zero coefficients (bound=0 -> "1": 0); moment dicts
    # compare equal only after dropping them, as add/scale/combine do
    target = {k: c for k, c in target_for(kind, comp, bound).items() if c != 0}
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs),
        "integrand": integrand,
        "target": target,
    }
