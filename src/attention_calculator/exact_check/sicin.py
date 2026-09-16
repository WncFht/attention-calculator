"""Exact checker for the si_q / cin_q ladder family (kernels/sicin.py).

params encode P = (au + bu·x + cu·x² + du·x³)/u multiplying
(1+x)^m (1−x)^n and the level-t_val Taylor-remainder kernel on q̄ = |power|.
The claimed vector is kdir·C(q̄) − s·bound over span{C_q, sin_q, cos_q, 1}
with kdir = s·σ (σ = sign(q) for the odd Si, +1 for the even Cin) — the
same parity fold kernels.sicin.solve uses.

``t_val`` is a trust boundary: t < 0 drops the polynomial tail onto the
*bare* kernel — sin(qx)/x loses sign-definiteness for q > π and would let
a crafted certificate certify a false claim — so it is rejected here, not
trusted from the search side. For t ≥ 0 kernel positivity is the
alternating-truncation lemma of docs/2026-09-16-w3-research-sicin.md
(strict on (0,1], sole zero at x = 0), a theorem rather than a numeric
check. The zero-integrand guard rejects the vacuous 0 = ∫ 0 dx > 0.
"""

from fractions import Fraction

from ..kernels.ln_pow import cubic_nonneg
from ..kernels.sicin import basis_moment, t_moments, u_moments
from ..moment import combine


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    m, n, t = int(params["m"]), int(params["n"]), int(params["t_val"])
    if t < 0:
        raise ValueError("t_val must be >= 0: the bare kernel is not sign-definite")
    u = Fraction(params["u_val"])
    coeffs = [Fraction(params[k]) / u for k in ("au_val", "bu_val", "cu_val", "du_val")]
    s = 1 if comp == ">" else -1
    sigma = (1 if power > 0 else -1) if kind == "si_q" else 1
    kdir = s * sigma
    target = {k: v for k, v in {kind: Fraction(kdir), "1": -s * bound}.items() if v}
    qq = abs(power)
    table = (t_moments if kind == "si_q" else u_moments)(m + n + 3, qq)
    basis = [basis_moment(kind, kdir, t, table, qq, m, n, j) for j in range(4)]
    integrand = combine(coeffs, basis)
    return {
        "identity_ok": integrand == target,
        "nonneg": cubic_nonneg(coeffs) and bool(integrand),
        "integrand": integrand,
        "target": target,
    }
