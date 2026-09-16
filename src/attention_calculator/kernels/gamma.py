"""gamma kernel family: prove ``power * gamma ⋚ bound`` for Euler's gamma.

Two-part composite identity on [0,1] — the bracket integrand is

    '<' (gamma < r):  x^N ((2-x)/2 - 1/(1-x) - 1/ln x) + sub
    '>' (gamma > r):  x^N (1/(1-x) + 1/ln x - 1/2)       + sub

with exact main-kernel values (digamma telescoping + a quadratic tail):

    '<':  int main = s_N - ln(N+1) - gamma,  s_N = H_N + 1/(N+1) - 1/(2(N+2))
    '>':  int main = gamma + ln(N+1) - r_N,  r_N = H_N + 1/(2(N+1))

so the remainder is a plain ln(N+1) bound, proved by a second site kernel
``x^{m'}(1-x)^{n'}(a+bx)/(1+Nx)^{s'}`` with ``s' = max(m', n', 1)`` — solved by
delegating to log_family.ln_bound_proof on q = N+1:

    '<':  sub = ln(N+1) - (s_N - r)   needs ln(N+1) > s_N - r
    '>':  sub = (r_N - r) - ln(N+1)   needs ln(N+1) < r_N - r

N selection (site convention, 79/79 probe parity): N = 0 iff the leftover is
already a nonneg constant — r >= s_0 = 3/4 for '<' (constant r - 3/4), r <=
r_0 = 1/2 for '>' (constant 1/2 - r); u_val = 0 then and c_val is 2 resp. 1.
Otherwise the site takes the smallest N >= 1 whose sub-integral is at least
the main-integral — i.e. each part gets at least half the total gap:

    '<':  s_N - ln(N+1) <= (r + gamma)/2        (float64 compare)
    '>':  r_N - ln(N+1) >= (r + gamma)/2

The site reports the sub-proof's (m', n', a, b, u) as its parameters and
stashes N in cu_val (c_val = 0).
"""

from fractions import Fraction
from math import gcd, log

import sympy as sp

from ..engine import NoSolution, WrongDirection
from ..render import coef_tex, lhs_tex, rat_tex, wire_fraction
from ..solve import EULER_F  # site compares in float64
from .log_family import ln_bound_proof, numerator_latex
from .quadlog import factors_tex

N_LIMIT = 600  # H_N denominators stay tractable well past this


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = False) -> dict:
    """prove(kind, power, comp, bound) -> site /calculate shape."""
    # the site's float64 direction pre-check runs before the kernel — a false
    # claim is 方向反了 (gamma 0>1), a true/equal one proceeds into the
    # bound/|power| division which crashes for power=0 (-> 500, gamma 0<1)
    c = float(power) * EULER_F
    if (float(bound) > c) if comp == ">" else (float(bound) < c):
        raise WrongDirection
    # power*gamma ⋚ bound  ≡  gamma ⋚ bound/power (direction flips if power < 0)
    if power < 0:
        comp = ">" if comp == "<" else "<"
    r = bound / abs(power)

    upper = comp == "<"
    n0_const = r - Fraction(3, 4) if upper else Fraction(1, 2) - r
    if n0_const >= 0:
        return {
            "parameters": {
                "m": 0,
                "n": 0,
                "a_val": str(n0_const),
                "b_val": "0",
                "c_val": "2" if upper else "1",
                "au_val": "0",
                "bu_val": "0",
                "cu_val": "0",
                "u_val": "0",
                "unified_form": {},
            },
            "solution": "a = 0, b = 0",
        }

    h = Fraction(0)
    for n in range(1, N_LIMIT + 1):
        h += Fraction(1, n)
        if upper:
            s = h + Fraction(1, n + 1) - Fraction(1, 2 * n + 4)  # s_N
            hit = float(s) - log(n + 1) <= (float(r) + EULER_F) / 2
        else:
            s = h + Fraction(1, 2 * n + 2)  # r_N
            hit = float(s) - log(n + 1) >= (float(r) + EULER_F) / 2
        if not hit:
            continue
        try:
            sub = ln_bound_proof(Fraction(n + 1), ">" if upper else "<", s - r)
        except (NoSolution, WrongDirection):
            continue
        sp = sub["parameters"]
        return {
            "parameters": {
                "m": sp["m"],
                "n": sp["n"],
                "a_val": sp["a_val"],
                "b_val": sp["b_val"],
                "c_val": "0",
                "au_val": sp["au_val"],
                "bu_val": sp["bu_val"],
                "cu_val": str(n),
                "u_val": sp["u_val"],
                "unified_form": {},
            },
            "solution": sub["solution"],
        }
    raise NoSolution


# ------------------------------------------------------------------- rendering


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Rebuild the site's get_integral_image LaTeX for solved parameters."""
    u = int(params["u_val"])
    n = int(params["cu_val"])  # main-kernel exponent and sub denominator coef

    # the coef is a real multiplier: it prefixes the kernel bracket, scales the
    # N=0 constant tail, and folds into the sub-fraction's numerator scale t
    cf = wire_fraction(power)
    ctex = coef_tex(power) + "\\gamma"
    btex = rat_tex(bound)
    lhs = lhs_tex(ctex, btex, comp)

    kern = (
        "\\dfrac{2-x}{2}-\\dfrac{1}{1-x}-\\dfrac{1}{\\ln(x)}"
        if comp == "<"
        else "\\dfrac{1}{1-x}+\\dfrac{1}{\\ln(x)}-\\dfrac{1}{2}"
    )
    x = sp.symbols("x")
    # coef·x^N goes through sympy when N > 0 ('1/2' -> \frac{x^{5}}{2});
    # the bare N = 0 multiplier keeps the request text ('3/2' -> \dfrac{3}{2})
    pre = coef_tex(power) if n == 0 else sp.latex(cf * x**n)
    main = f"{pre}\\left({kern}\\right)"

    if u == 0:
        m, nn = int(params["m"]), int(params["n"])
        au, bu = int(params["au_val"]), int(params["bu_val"])
        if au == 0 and bu == 0:
            # N = 0: leftover is the bare constant a_val * coef
            const = Fraction(params["a_val"]) * cf
            body = main if const == 0 else f"\\left[{main}+{rat_tex(const)}\\right]"
        else:
            # u = 0 denominator collapses the sub-fraction to sympy's
            # zoo * numerator (site shows \tilde{\infty} times the factors)
            s = max(m, nn, 1)
            expr = cf * x**m * (1 - x) ** nn * (au + bu * x) / (u * (1 + n * x) ** s)
            body = f"\\left[{main}+{factors_tex(expr)}\\right]"
    else:
        m, nn = int(params["m"]), int(params["n"])
        au, bu = int(params["au_val"]), int(params["bu_val"])
        s = max(m, nn, 1)
        # the sub-fraction is cf·num/(u·den) printed reduced: the scalar
        # cf/u drops to (cf.n/g) / (cf.d·u/g), so '2' over u=936 shows an
        # unscaled numerator over 468 and '1/2' over u=168 shows 336
        g = gcd(cf.numerator, u)
        tnum = cf.numerator // g
        uden = cf.denominator * (u // g)
        num = numerator_latex(m, nn, au, bu, 0, tnum, Fraction(n), False)
        if s == 1:
            den = sp.latex(uden * n * x + uden)
        else:
            inner = "x + 1" if n == 1 else f"{n} x + 1"
            den = f"\\left({inner}\\right)^{{{s}}}"
            if uden != 1:
                den = f"{uden} {den}"
        body = f"\\left[{main}+\\frac{{{num}}}{{{den}}}\\right]"

    eq = " =\\int_0^1 " if comp == "<" else " = \\int_0^1 "
    return f"{lhs}{eq}{body} \\mathrm{{d}} x > 0"
