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
    alpha = Fraction(params["c_val"])
    # cos(pi alpha/2) is the claimed sin(pi q)/cos(pi q) only when c_val
    # matches the request-derived frequency (q_eff = degrees/180); a
    # mismatched c_val makes the integrand a different constant's moment,
    # so it stays keyed "C" and the identity comparison fails honestly.
    q_eff = power / 180 if kind.endswith("_degree") else power
    const_key = sym if alpha == (1 - 2 * q_eff if is_sin else 2 * q_eff) else "C"
    t = [angle_moment(j, alpha) for j in range(m + n + 2)]
    basis = [basis_moment(t, m, n, j) for j in range(2)]
    raw = combine(coeffs, basis)
    integrand = {(const_key if k == "C" else k): v for k, v in raw.items()}
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs),
        "integrand": integrand,
        "target": target,
    }
