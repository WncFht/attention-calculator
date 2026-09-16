"""Exact checker for beta_even kernels: beta4..beta10 (exact-mode only).

params encode P = (au + bu·x²)/u multiplying x^{2m}(1−x²)^n·ln^r(1/x)/(1+x²);
the true basis moments come from the kernel's own basis() (ln_moment's even
branch, whose symbol coefficient is already β(r+1) — no conversion factor),
and the printed claim s·(power·β(r+1) − bound) is solved directly.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.beta_even import basis
from ..moment import Moment, combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(str(params["u_val"]))
    coeffs = [Fraction(str(params["au_val"])) / u, Fraction(str(params["bu_val"])) / u]

    integrand = combine(coeffs, basis(m, n, kind))
    sign = 1 if comp == ">" else -1
    target: Moment = {k: v for k, v in {kind: sign * power, "1": -sign * bound}.items() if v}
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
