"""zeta_odd kernel family: power·ζ(s), s ∈ {5,7,9,11} — exact-mode only.

W3 extension types with no site counterpart (plan doc W3.1: the quadlog ln^r
moment machine already covers them). Same integrand as zeta3 with r = s−1:

    x^{2m+1} (1-x^2)^n (a + b x^2) ln^r(1/x) / (1 + x^2)   on [0,1]

Odd powers keep moments inside span{ζ(r+1), 1}: ln_moment's odd branch is
(-1)^k r!/2^{r+1}·(η(r+1) − partial) with η(s) = (1−2^{1−s})ζ(s), so the
η-coefficient is scaled by 1−2^{−r} — 15/16 for ζ(5), 63/64 for ζ(7) — to
land on the ζ symbol. (Even powers would land on β(r+1) ∈ Q·π^{r+1} and
could not express an odd zeta value.)
"""

from fractions import Fraction

import sympy as sp

from ..engine import mn_order, search
from ..moment import Moment
from ..render import coef_tex, emit, lhs_tex, rat_tex
from .quadlog import basis_moments, factors_tex, ln_moment, numerator

# kind -> (ln exponent r, η(r+1)/ζ(r+1) = 1−2^{−r}): converts ln_moment's
# η-symbol coefficient onto the printed ζ(r+1)
ZETA = {
    "zeta5": (4, Fraction(15, 16)),
    "zeta7": (6, Fraction(63, 64)),
    "zeta9": (8, Fraction(255, 256)),
    "zeta11": (10, Fraction(1023, 1024)),
}

LIMIT = 10


def basis(m: int, n: int, kind: str) -> list[Moment]:
    """True moments of x^{2m+1}(1−x²)^n·{1, x²}·ln^r(1/x)/(1+x²) for (a, b)."""
    r, factor = ZETA[kind]

    def term(k: int) -> tuple[Fraction, Fraction]:
        cc, rat = ln_moment(k, r, True)
        return cc * factor, rat

    return basis_moments(m, n, True, kind, term)


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """Search (m, n) in the author's order; return the parameter/solution dict.

    Raises WrongDirection/NoSolution from ..engine on failure. ``exact`` is
    accepted for signature compatibility — this family has no site path.
    """
    plans = ((m, n, basis(m, n, kind)) for m, n in mn_order(LIMIT))
    sign = 1 if comp == ">" else -1
    target = {kind: sign * power, "1": -sign * bound}
    solved = search(plans, target, nonneg=True)
    return emit(solved.m, solved.n, solved.coeffs)


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """The proof-equation LaTeX, mirroring quadlog's zeta3 layout (q = 1)."""
    r = ZETA[kind][0]
    m, n = int(params["m"]), int(params["n"])
    au, bu, u = int(params["au_val"]), int(params["bu_val"]), int(params["u_val"])

    lhs = lhs_tex(coef_tex(power) + f"\\zeta({r + 1})", rat_tex(bound), comp)

    x = sp.symbols("x")
    num = numerator(m, n, True, au, bu, 1)
    den = u * (x**2 + 1)
    ratio = num / den
    if ratio.as_numer_denom()[1] == 1:
        body = sp.latex(ratio)
    else:
        body = f"\\frac{{{factors_tex(num)}}}{{{sp.latex(den)}}}"
    return f"{lhs} = \\int_0^1 {body}\\ln^{r}(x) \\mathrm{{d}} x > 0"
