"""Exact checker for the gamma kernel — the two-part composite integrand.

The printed integrand is ``cf*x^N*K_dir(x) + sub`` with ``cf = power`` — the
request coefficient multiplies every printed term (kernels/gamma
render_equation), so a consistent solve stays consistent for any power.

Main-kernel moments, N = cu_val (verified at 50dps against mpmath, and
derived by splitting x^N = (x^N-1)+1 against the classical
int(1/ln x + 1/(1-x))dx = gamma, int((x^N-1)/ln x)dx = ln(N+1),
int((x^N-1)/(1-x))dx = -H_N):

    '<': int x^N ((2-x)/2 - 1/(1-x) - 1/ln x) = s_N - ln(N+1) - gamma
    '>': int x^N (1/ln x + 1/(1-x) - 1/2)   = gamma + ln(N+1) - r_N

with s_N = H_N + 1/(N+1) - 1/(2(N+2)), r_N = H_N + 1/(2(N+1)). The "ln"
symbol denotes ln(N+1); the sub-fraction's ln moment must cancel the main
term's exactly.

u_val = 0 marks the N=0 degenerate print (cf*K_0 + cf*a_val, au=bu=0).
Otherwise sub is the ln(N+1) sub-proof fraction
(cf/u)*x^m(1-x)^n(au+bu*x)/(1+Nx)^s, s = max(m, n, 1), moments over
{ln(N+1), 1} via log_family.basis_moment.

nonneg: the brackets K_dir >= 0 on (0,1) are the kernel's structural lemma
(the author's proof form); over QQ the check certifies cf > 0 plus the
sub-polynomial resp. tail constant >= 0.
"""

from fractions import Fraction

from ..engine import poly_nonneg
from ..kernels.log_family import basis_moment
from ..moment import Moment, add, combine, scale


def main_moment(comp: str, n: int) -> Moment:
    """Exact moment of x^N*K_dir over {"gamma", "ln"=ln(N+1), "1"}."""
    h = sum(Fraction(1, j) for j in range(1, n + 1))
    if comp == "<":
        s = h + Fraction(1, n + 1) - Fraction(1, 2 * n + 4)  # s_N
        moment: Moment = {"gamma": Fraction(-1), "1": s}
    else:
        r = h + Fraction(1, 2 * n + 2)  # r_N
        moment = {"gamma": Fraction(1), "1": -r}
    if n:  # ln(1) = 0 — the N=0 main term carries no ln part
        moment["ln"] = Fraction(-1 if comp == "<" else 1)
    return moment


def check(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Verify the claimed identity in exact moments; see exact_check."""
    n_k = int(params["cu_val"])
    u = Fraction(params["u_val"])
    a, au, bu = (Fraction(params[k]) for k in ("a_val", "au_val", "bu_val"))
    sign = Fraction(1 if comp == ">" else -1)
    target = {k: v for k, v in {"gamma": sign * power, "1": -sign * bound}.items() if v}
    integrand = scale(main_moment(comp, n_k), power)
    if u == 0:
        if au or bu:  # u=0 with a numerator is the site's zoo print — not an integrand
            return {"identity_ok": False, "nonneg": False, "integrand": integrand, "target": target}
        integrand = add(integrand, {"1": power * a})
        nonneg = power > 0 and a >= 0
    else:
        m, n = int(params["m"]), int(params["n"])
        basis = [basis_moment(Fraction(n_k), max(m, n, 1), m, n, j, False) for j in (0, 1)]
        integrand = add(integrand, combine([power * au / u, power * bu / u], basis))
        # sign rule on the true coefficients au/u, bu/u — params_domain's
        # u_val >= 0 gate makes this equivalent to [au, bu], but the check
        # must stay correct even if called without the gate
        nonneg = power > 0 and poly_nonneg([au / u, bu / u])
    return {
        "identity_ok": integrand == target,
        "nonneg": nonneg,
        "integrand": integrand,
        "target": target,
    }
