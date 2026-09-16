"""Exact checker for the Dixon family: pi3, pi3_u, pi3_a.

params encode ``P = (au + bu·x³)/u`` multiplying ``x^{3m+res}(1−x)`` and the
kernel ``outer·(1−x³)^{−k/3}`` on [0, 1]. (kernel, res, outer) are fully
determined by (kind, comp) — the six configurations are a bijection — and
``n`` is a direction flag, not an exponent (the beta_family/pi_sqrt2
encoding; it carries no claim content, so the check ignores it).

cu_val != 0 marks the ``+b`` fallback: ``t·x^{3m+res}(1−x)·outer·K + b`` with
t = au_val/u_val and b = b_val, verified against the slot-0 basis moment at
the reported m.

P = a + b·x³ is monotone on [0,1], so poly_nonneg's endpoint rule
(a ≥ 0 ∧ a+b ≥ 0) is the exact certificate — the same rule a+bx⁴ uses.
The zero-integrand guard rejects the vacuous ``0 = ∫ 0 dx > 0``.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.dixon import CONFIG, dixon_basis, dixon_target, third_table
from ..moment import add, combine, scale


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    kernel, res, outer = CONFIG[(kind, comp)]
    m = int(params["m"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params["au_val"]) / u, Fraction(params["bu_val"]) / u]
    target = dixon_target(kind, comp, power, bound)
    tab = third_table(kernel, m + 3)
    if int(params["cu_val"]):
        mom = dixon_basis(kernel, res, outer, tab, m, 0)
        t, b = coeffs[0], Fraction(params["b_val"])
        integrand = add(scale(mom, t), {"1": b})
        return {
            "identity_ok": integrand == target,
            "nonneg": t >= 0 and b >= 0 and bool(integrand),
            "integrand": integrand,
            "target": target,
        }
    basis = [dixon_basis(kernel, res, outer, tab, m, i) for i in (0, 1)]
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
