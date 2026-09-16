"""Exact checker for quadlog kernels: pi, pi_n, catalan, zeta3, arctan_q, arccot_q.

params encode P = (au + bu*x^2)/u multiplying x^{2m+odd}(1-x^2)^n·K; the true
basis moments are the kernel's own ln_moment / atan_moment via basis_moments —
this family has no reproduced site bias, and it solves the printed claim
s·(power·C − bound) directly (pi_n: s·(pi^p − bound^q)), so site-mode params
are expected to verify here.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.quadlog import atan_moment, basis_moments, ln_moment, spec
from ..moment import Moment, combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(str(params["u_val"]))
    coeffs = [Fraction(str(params["au_val"])) / u, Fraction(str(params["bu_val"])) / u]
    cfg = spec(kind, power)
    q, r, odd, sym = cfg["q"], cfg["r"], cfg["odd"], cfg["sym"]

    if sym == "arctan":

        def term(k: int) -> tuple[Fraction, Fraction]:
            return atan_moment(k, q)
    else:
        factor = cfg["factor"]

        def term(k: int) -> tuple[Fraction, Fraction]:
            cc, rat = ln_moment(k, r, odd)
            return cc * factor, rat

    integrand = combine(coeffs, basis_moments(m, n, odd, sym, term))
    sign = 1 if comp == ">" else -1
    bound = bound ** cfg.get("pd", 1)  # pi_n p/q: claim pi^p vs bound^q
    target: Moment = {k: v for k, v in {sym: sign * cfg["coef"], "1": -sign * bound}.items() if v}
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
