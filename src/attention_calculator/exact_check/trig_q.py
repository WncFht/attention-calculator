"""Exact checker for radian trig kernels: sin_q, cos_q, tan_q, cot_q.

params encode P = (au + bu*x + cu*x**2)/u multiplying x^m (1-x)^n sin(q x)
on [0, 1]. The basis is the kernel's own sin_moments/basis_moment -- this
family has no reproduced site bias, so solve and check share it.

tan/cot targets keep the kernel's cleared-denominator form: the printed
``tan q ~ bound = (1/cos q) * integral`` means the real content of the
identity is ``sin q - bound*cos q ~ 0`` (cot analogously over sin q), so
the claimed vector lives in {sin_q, cos_q, 1} exactly as solve() builds it.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.trig_q import basis_moment, sin_moments
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val")]
    s = Fraction(1 if comp == ">" else -1)
    if kind == "sin_q":
        raw = {"sin_q": s, "1": -s * bound}
    elif kind == "cos_q":
        raw = {"cos_q": s, "1": -s * bound}
    elif kind == "tan_q":
        raw = {"sin_q": s, "cos_q": -s * bound}
    else:
        raw = {"cos_q": s, "sin_q": -s * bound}
    target = {k: v for k, v in raw.items() if v}
    moments = sin_moments(m + n + 2, power)
    basis = [basis_moment(moments, m, n, j) for j in range(3)]
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs),
        "integrand": integrand,
        "target": target,
    }
