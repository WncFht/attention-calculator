"""Exact checker for exp kernels: e, e_q, e_pi.

params encode P = (au + bu*x)/u (e_pi: (au + bu*sin x)/u); the true basis
moments are the kernel's own basis_x_moment / basis_sin_moment — this
family has no reproduced site bias, so solve and check share them.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.exp_family import basis_sin_moment, basis_x_moment
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = params["m"], params["n"]
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params["au_val"]) / u, Fraction(params["bu_val"]) / u]
    sign = Fraction(1 if comp == ">" else -1)
    if kind == "e_pi":
        sym, coef = "e_pi", power
        basis = [basis_sin_moment(m, n, j) for j in (0, 1)]
    elif kind == "e_q":
        sym, coef, q = "e_q", Fraction(1), power
        basis = [basis_x_moment(m, n, j, q, sym) for j in (0, 1)]
    else:
        sym, coef, q = "e", power, Fraction(1)
        basis = [basis_x_moment(m, n, j, q, sym) for j in (0, 1)]
    integrand = combine(coeffs, basis)
    target = {sym: sign * coef, "1": -sign * bound}
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs),
        "integrand": integrand,
        "target": target,
    }
