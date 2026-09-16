"""Exact checker for beta kernels: golden, varpi, gauss.

params encode P = (au + bu*x)/u on x^m(1-x)^n*sqrt(x+4) for golden, and
P = (au + bu*x^4)/u on x^{4m+res}(1-x)/(pi^e*sqrt(1-x^4)) for the
lemniscate pair (n is a direction flag there, not an exponent). The true
basis moments are the kernel's own golden_basis / lemniscate_basis — the
site's transposed '<' solve emits (a, b) satisfying the wrong system, so
those params fail identity_ok here, which is the point.

cu_val != 0 marks the site's fallback templates — cu_val=1 ('>') prints
t*poly*sqrt(1-x^4)/pi + b and cu_val=2 ('<') prints
t*x^res(1-x)*sqrt(1-x^4) + b, with t = au_val/u_val and b = b_val; the
per-kind moment constants live in SQRT_BOUND_MOMENT / LT_BOUND_MOMENT.

An identically-zero integrand (P == 0, or t == b == 0) makes the printed
``... = int 0 dx > 0`` assert 0 > 0 — the equality holds vacuously but the
inequality is false, so nonneg reports False (golden/varpi/gauss 0-vs-0).
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.beta_family import (
    LT_BOUND_MOMENT,
    SQRT_BOUND_MOMENT,
    golden_basis,
    lemniscate_basis,
    lemniscate_target,
)
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params["au_val"]) / u, Fraction(params["bu_val"]) / u]
    cu = int(params["cu_val"])
    if kind != "golden" and cu:
        # fallback template: moment = t*{sym: alpha, "1": beta} + {1: b}
        sym, alpha, beta = {1: SQRT_BOUND_MOMENT, 2: LT_BOUND_MOMENT}[cu][kind]
        t, b = coeffs[0], Fraction(params["b_val"])
        integrand = {k: v for k, v in {sym: t * alpha, "1": t * beta + b}.items() if v}
        target = nonzero(lemniscate_target(kind, comp, power, bound))
        return {
            "identity_ok": integrand == target,
            "nonneg": t >= 0 and b >= 0 and bool(integrand),
            "integrand": integrand,
            "target": target,
        }
    if kind == "golden":
        basis = [golden_basis(int(params["m"]), int(params["n"]), i) for i in (0, 1)]
        sign = Fraction(1 if comp == ">" else -1)
        target = {"phi": sign * power, "1": -sign * bound}
    else:
        basis = [lemniscate_basis(kind, comp, int(params["m"]), i) for i in (0, 1)]
        target = lemniscate_target(kind, comp, power, bound)
    target = nonzero(target)
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }


def nonzero(m: dict) -> dict:
    """Drop zero-valued keys; Moments omit zero coefficients."""
    return {k: v for k, v in m.items() if v}
