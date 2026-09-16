"""Exact checker for pi/degree trig kernels: sin_pi_q, cos_pi_q,
sin_q_degree, cos_q_degree.

params encode P = (au + bu*sin x)/u multiplying sin^m x (1 - sin x)^n
sin(alpha x) on [0, pi/2]; c_val carries the kernel frequency alpha and
``power`` is the angle argument q itself (degrees for *_q_degree), so the
printed claim is ``sin(pi q) ~ bound`` / ``cos(pi q) ~ bound`` with no
power coefficient.

The basis MUST be the kernel's angle_moment/basis_moment true moments --
never site_basis, whose stored (1, 8) formula carries the BIAS18
corruption. Site params emitted through that corrupted moment fail the
identity here, which is the point.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.trig_pi import angle_moment, basis_moment
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params["au_val"]) / u, Fraction(params["bu_val"]) / u]
    s = Fraction(1 if comp == ">" else -1)
    is_sin = kind.startswith("sin")
    sym = "sin_pi_q" if is_sin else "cos_pi_q"
    target = {k: v for k, v in {sym: s, "1": -s * bound}.items() if v}
    # c_val is the kernel frequency alpha — recompute it from (kind, power)
    # rather than trusting the certificate: a forged alpha prints a different
    # factor, so the encoded integrand is not the claimed one and the
    # identity fails honestly (q_eff = degrees/180)
    q_eff = power / 180 if kind.endswith("_degree") else power
    expected = 1 - 2 * q_eff if is_sin else 2 * q_eff
    alpha = Fraction(params["c_val"])
    if alpha != expected:
        return {"identity_ok": False, "nonneg": False, "integrand": {}, "target": {}}
    t = [angle_moment(j, alpha) for j in range(m + n + 2)]
    basis = [basis_moment(t, m, n, j) for j in range(2)]
    raw = combine(coeffs, basis)
    integrand = {(sym if k == "C" else k): v for k, v in raw.items()}
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
