"""Exact checker for log kernels: ln_q, ln_q_square, artanh_q, arcoth_q.

params encode P = (au + bu*x [+ cu*x^2])/u over the basis x^m(1-x)^n with
denominator (1+cx)^s, s = max(m, n, 1) (site convention, log-notes.md).
c = q - 1 for the ln types; for artanh/arcoth the renderer rebuilds the
denominator from the reduced ln argument q~ stashed in c_val, so the check
does the same. The a, b emitted for artanh/arcoth already carry the 1/2
reduction factor — the printed claim s*(artanh q - bound) is exactly the
{"ln": s/2, "1": -s*bound} vector the kernel solved.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.log_family import basis_moment
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    square = kind == "ln_q_square"
    coeffs = [Fraction(params["au_val"]) / u, Fraction(params["bu_val"]) / u]
    if square:
        coeffs.append(Fraction(params["cu_val"]) / u)
    qt = Fraction(params["c_val"]) if kind in ("artanh_q", "arcoth_q") else power
    s = max(m, n, 1)
    basis = [basis_moment(qt - 1, s, m, n, j, square) for j in range(3 if square else 2)]
    sign = Fraction(1 if comp == ">" else -1)
    if square:
        target = {"ln2": sign, "1": -sign * bound}
    else:
        coef = Fraction(1) if kind == "ln_q" else Fraction(1, 2)
        target = {"ln": sign * coef, "1": -sign * bound}
    target = {k: v for k, v in target.items() if v}  # Moments omit zero coeffs
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
