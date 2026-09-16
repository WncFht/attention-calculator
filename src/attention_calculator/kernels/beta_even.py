"""beta_even kernel family: power·β(s), s ∈ {4, 6, 8, 10} — exact-mode only.

W7 extension types with no site counterpart. Same 1/(1+x²) kernel as
quadlog's catalan, at the next even-argument rungs of the Dirichlet beta
ladder — β(4) has no known closed form (the even analogue of Catalan
β(2) = C). The integrand uses *even* x-powers with odd ln exponent:

    x^{2m} (1-x^2)^n (a + b x^2) ln^r(1/x) / (1 + x^2)   on [0,1], r = s-1

ln_moment's even branch is (-1)^k r!·(β(r+1) − partial), so the moment's
symbol coefficient lands directly on β(r+1) — factor 1, no η→ζ-style
normalization. (Odd powers at the same r would land on η(r+1) ∈ Q·π^{r+1}
and could not express β(4).) Odd r keeps ln^r(1/x) ≥ 0 on the domain, so
the printed equation must keep the 1/x argument — ln^r(x) would flip the
integrand's sign.
"""

from fractions import Fraction

import sympy as sp

from ..engine import mn_order, search
from ..moment import Moment
from ..render import coef_tex, emit, lhs_tex, rat_tex
from .quadlog import basis_moments, factors_tex, ln_moment, numerator

# kind -> ln exponent r = s-1; ln_moment's even branch already lands on
# β(r+1) with coefficient cc, so no conversion factor is needed
BETA = {
    "beta4": 3,
    "beta6": 5,
    "beta8": 7,
    "beta10": 9,
}

LIMIT = 10


def basis(m: int, n: int, kind: str) -> list[Moment]:
    """True moments of x^{2m}(1−x²)^n·{1, x²}·ln^r(1/x)/(1+x²) for (a, b)."""
    r = BETA[kind]
    return basis_moments(m, n, False, kind, lambda k: ln_moment(k, r, False))


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
    """The proof-equation LaTeX, mirroring zeta_odd's layout (q = 1).

    r is odd, so the weight prints as ln^r(1/x): ln^r(x) would carry the
    opposite sign and contradict the '> 0' tail.
    """
    r = BETA[kind]
    m, n = int(params["m"]), int(params["n"])
    au, bu, u = int(params["au_val"]), int(params["bu_val"]), int(params["u_val"])

    lhs = lhs_tex(coef_tex(power) + f"\\beta({r + 1})", rat_tex(bound), comp)

    x = sp.symbols("x")
    num = numerator(m, n, False, au, bu, 1)
    den = u * (x**2 + 1)
    ratio = num / den
    if ratio.as_numer_denom()[1] == 1:
        body = sp.latex(ratio)
    else:
        body = f"\\frac{{{factors_tex(num)}}}{{{sp.latex(den)}}}"
    return f"{lhs} = \\int_0^1 {body}\\ln^{{{r}}}(1/x) \\mathrm{{d}} x > 0"
