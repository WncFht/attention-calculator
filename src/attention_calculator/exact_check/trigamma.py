"""Exact checker for psi1_q — the telescope trigamma kernel pair.

params encode ``P = (au + bu*x)/u`` multiplying ``x^m (1-x)^n`` and the
direction kernel on [0, 1]: '>' uses ``x^{q-1}(x-1-ln x)/(1-x)``, '<' uses
``x^{q+N-2}(1-x+x ln x)/(1-x)`` with N the least integer >= 0 giving
q+N > 1 — N is derived from q, not a stored field (a larger N only
relabels m: U_u(N+1) = U_{u+1}(N)).

Kernel signs are the exact lemma ln t <= t-1 (t > 0): x-1-ln x >= 0 is the
lemma at t = x, and 1-x+x ln x >= 0 is the same lemma at t = 1/x; each
factor is positive on (0,1) and vanishes only at the endpoint x = 1, so
the integrand is sign-definite iff P is — ``poly_nonneg(coeffs)`` plus the
bool(integrand) guard against the vacuous ``0 = int 0 dx > 0`` is the full
certificate. The moment basis is rebuilt from the kernel's own
``moments``/``basis_moment`` over span{psi1_q, 1}; any (1-x)^n factor in
forged params cancels the psi' row and fails the dict equality honestly.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.trigamma import basis_moment, moments
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n = int(params["m"]), int(params["n"])
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val")]
    s = Fraction(1 if comp == ">" else -1)
    target = {k: v for k, v in {kind: s, "1": -s * bound}.items() if v}
    moms = moments(kind, comp, m + n + 1, power)
    basis = [basis_moment(moms, m, n, j) for j in range(2)]
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": poly_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
