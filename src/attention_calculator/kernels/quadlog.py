"""quadlog kernel family: denominator 1+(qx)^2 times an ln^r(1/x) weight.

All six site types share the integrand
    x^{2m+par} (1-x^2)^n (a + b x^2) ln^r(1/x) / (1 + (qx)^2)   on [0,1]
and solve for (a, b) in a 2-dimensional moment space span{C, 1}:

    type       r      parity (x-exponent)   constant        power meaning   limit
    pi         0      even  x^{2m}          pi              coef of pi      30
    pi_n       k-1    even if k odd else odd (x^{2m+1})  pi^k    exponent k      10
    catalan    1      even                Catalan C       coef of C       10
    zeta3      2      odd                 zeta(3)         coef of zeta(3) 10
    arctan_q   0      even                arctan q        argument q      10
    arccot_q   0      even                arccot q = arctan(1/q)          10

Moment closed forms used below (docs/kernel-spec.md):
    ∫x^{2k}   ln^r(1/x)/(1+x^2) = (-1)^k r! (β(r+1) - Σ_{i<k} (-1)^i/(2i+1)^{r+1})
    ∫x^{2k+1} ln^r(1/x)/(1+x^2) = (-1)^k r!/2^{r+1} (η(r+1) - Σ_{l=1..k} (-1)^{l-1}/l^{r+1})
    ∫x^{2k}/(1+(qx)^2)          = (-1)^k q^{-(2k+1)} arctan q + Q(q), via
                                  q^2 T_k + T_{k-1} = 1/(2k-1)
β(1)=pi/4, β(2)=C, β(odd)=Q·pi^odd, η(even)=Q·pi^even — which is why pi_n with
even k needs odd powers (η(k) ~ pi^k) and odd k needs even powers (β(k) ~ pi^k).
"""

import math
from fractions import Fraction
from functools import cache
from itertools import pairwise
from math import comb, factorial, gcd, lcm

import sympy as sp

from ..engine import WrongDirection, mn_order, search
from ..moment import Moment
from ..render import coef_tex, rat_tex, wire_or, wire_pair

# β(k)/pi^k for odd k (Euler numbers: β(2j+1) = (-1)^j E_{2j} pi^{2j+1}/(4^{j+1}(2j)!))
BETA_PI = {1: Fraction(1, 4), 3: Fraction(1, 32), 5: Fraction(5, 1536),
           7: Fraction(61, 184320), 9: Fraction(277, 8257536)}
# η(k)/pi^k for even k (η(k) = (1-2^{1-k})·ζ(k))
ETA_PI = {2: Fraction(1, 12), 4: Fraction(7, 720), 6: Fraction(31, 30240),
          8: Fraction(127, 1209600), 10: Fraction(73, 6842880)}


@cache
def ln_moment(k: int, r: int, odd: bool) -> tuple[Fraction, Fraction]:
    """∫_0^1 x^{2k+odd} ln^r(1/x)/(1+x^2) dx = cc·S + rat.

    S is β(r+1) for even powers, η(r+1) for odd powers; both are returned only
    through their coefficient cc — the caller supplies S's value in the target
    constant symbol.
    """
    if odd:
        cc = Fraction((-1) ** k * factorial(r), 2 ** (r + 1))
        partial = sum(Fraction((-1) ** (i - 1), i ** (r + 1)) for i in range(1, k + 1))
    else:
        cc = Fraction((-1) ** k * factorial(r))
        partial = sum(Fraction((-1) ** i, (2 * i + 1) ** (r + 1)) for i in range(k))
    return cc, -cc * partial


@cache
def atan_moment(k: int, q: Fraction) -> tuple[Fraction, Fraction]:
    """∫_0^1 x^{2k}/(1+(qx)^2) dx = cc·arctan(q) + rat, via q^2 T_k + T_{k-1} = 1/(2k-1)."""
    cc, rat = 1 / q, Fraction(0)
    for j in range(1, k + 1):
        cc, rat = -cc / q**2, (Fraction(1, 2 * j - 1) - rat) / q**2
    return cc, rat


def basis_moments(m: int, n: int, odd: bool, sym: str,
                  term) -> list[Moment]:
    """Moments of x^{2m+odd}(1-x^2)^n·{1, x^2}·K for the unknowns (a, b).

    ``term(k)`` returns (coef_of_symbol, rational) for the x^{2k+odd} power.
    """
    out = []
    for c in (0, 1):
        v: Moment = {}
        for i in range(n + 1):
            w = Fraction((-1) ** i) * comb(n, i)
            cc, rat = term(m + i + c)
            if cc:
                v[sym] = v.get(sym, Fraction(0)) + w * cc
            if rat:
                v["1"] = v.get("1", Fraction(0)) + w * rat
        out.append(v)
    return out


def spec(kind: str, power: Fraction) -> dict:
    """Per-type kernel configuration shared by prove() and render_equation()."""
    if kind == "pi":
        return dict(r=0, odd=False, factor=Fraction(1, 4), sym="pi",
                    coef=power, q=Fraction(1), limit=30)
    if kind == "pi_n":
        # power p/q is reduced to pi^p vs bound^q; the kernel uses ln^{p-1}.
        # no range guard — the site lets the table lookup crash (KeyError ->
        # 500) for exponents outside [1, 10]
        k, pd = power.numerator, power.denominator
        return dict(r=k - 1, odd=k % 2 == 0, sym="pi", coef=Fraction(1),
                    factor=BETA_PI[k] if k % 2 else ETA_PI[k],
                    q=Fraction(1), limit=10, pd=pd)
    if kind == "catalan":
        return dict(r=1, odd=False, factor=Fraction(1), sym="catalan",
                    coef=power, q=Fraction(1), limit=10)
    if kind == "zeta3":
        return dict(r=2, odd=True, factor=Fraction(3, 4), sym="zeta3",
                    coef=power, q=Fraction(1), limit=10)
    q = power if kind == "arctan_q" else 1 / power  # arccot_q -> arctan(1/q)
    return dict(r=0, odd=False, sym="arctan", coef=Fraction(1), q=q, limit=10)


# float64 direction pre-check, same convention as the other kernels: the site
# evaluates the claimed constant before the search and rejects strictly-false
# inequalities with 方向反了. Only the types whose kernel can crash need it —
# arctan/arccot divide by q and pi_n's spec() dies on out-of-range exponents.
PRE_F = {"arctan_q": math.atan, "arccot_q": lambda v: math.atan(1 / v),
         "pi_n": lambda v: math.pi ** v}


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """Search (m, n) in the author's order; return the site's parameter dict."""
    if kind in PRE_F:
        c = PRE_F[kind](float(power))
        if (float(bound) > c) if comp == ">" else (float(bound) < c):
            raise WrongDirection
    cfg = spec(kind, power)
    q, r, odd, sym = cfg["q"], cfg["r"], cfg["odd"], cfg["sym"]

    if sym == "arctan":
        def term(k: int) -> tuple[Fraction, Fraction]:
            return atan_moment(k, q)
    else:
        factor = cfg["factor"]
        def term(k: int) -> tuple[Fraction, Fraction]:
            cc, rat = ln_moment(k, r, odd)
            return cc * factor, rat

    plans = ((m, n, basis_moments(m, n, odd, sym, term))
             for m, n in mn_order(cfg["limit"]))
    sign = 1 if comp == ">" else -1
    bound = bound ** cfg.get("pd", 1)  # pi_n p/q: solve pi^p vs bound^q
    target = {sym: sign * cfg["coef"], "1": -sign * bound}
    solved = search(plans, target, nonneg=True)

    a, b = solved.coeffs
    u = lcm(a.denominator, b.denominator)
    return {
        "parameters": {
            "m": solved.m, "n": solved.n,
            "a_val": str(a), "b_val": str(b), "c_val": "0",
            "au_val": str(a * u), "bu_val": str(b * u), "cu_val": "0",
            "u_val": str(u), "unified_form": {},
        },
        "solution": f"a = {a}, b = {b}",
    }


# ---------------------------------------------------------------- LaTeX render

def const_tex(kind: str, power: Fraction | str) -> str:
    """The target constant as printed on the equation's left-hand side."""
    if kind == "pi":
        return coef_tex(power) + "\\pi"
    if kind == "pi_n":
        # raw (num, den) pair — unreduced: '\frac{-9}{8}' keeps its den 8 and
        # the macro drops its minus inside the pi^... slot (verbatim echo)
        pnum, pden = wire_pair(power)
        if pden != 1:
            return f"\\left(\\pi^{rat_tex(power)}\\right)"
        return f"\\pi^{{{pnum}}}" if pnum >= 10 else f"\\pi^{pnum}"
    if kind == "arctan_q":
        return "\\arctan" + rat_tex(power)
    if kind == "arccot_q":
        return "\\mathrm{arccot}" + rat_tex(power)
    if kind == "catalan":
        return coef_tex(power) + "C"
    return coef_tex(power) + "\\zeta(3)"


def numerator(m: int, n: int, odd: bool, au: int, bu: int, t: int) -> sp.Expr:
    """Numerator t·x^{2m+odd}·(1-x²)^n·(au+bu·x²) as a sympy expression.

    The site binds t into the first surviving factor (x-part, then the basis,
    then P); the Mul then canonicalizes — a monomial or constant P merges into
    the leading coefficient/power, a printed P equal to the basis factor folds
    into its power (so fold needs t·au==1, t·bu==-1 whenever t scales the
    basis), and the Add factors come out in sympy's print order — e.g.
    ``(1-x²)`` before ``(3-3x²)``.
    """
    x = sp.symbols("x")
    pieces = []
    if 2 * m + odd:
        pieces.append(x ** (2 * m + odd))
    if n:
        pieces.append((1 - x ** 2) ** n)
    # the polynomial factor multiplies in even when au == bu == 0, so the
    # degenerate params collapse the whole product to a literal 0 (site
    # renders `0\ln^2(x)` for zeta3 0-vs-0, not x/(x²+1))
    pieces.append(au + bu * x ** 2)
    pieces[0] = t * pieces[0]
    return sp.Mul(*pieces)


def factors_tex(num: sp.Expr) -> str:
    """Per-factor sympy latex joined by the site cdot rule.

    `` \\cdot `` iff the left piece ends in '}' and the right is a
    ``\\left(<digit>...\\right)`` group — a parenthesized Add carrying no
    outer exponent (same rule as beta_family.join_cdot).
    """
    args = num.as_ordered_factors()
    wrap = len(args) > 1
    tex = [f"\\left({sp.latex(a)}\\right)" if a.is_Add and wrap else sp.latex(a)
           for a in args]
    out = tex[0]
    for prev, cur in pairwise(tex):
        cdot = (prev.endswith("}") and cur.startswith("\\left(")
                and cur[6].isdigit() and cur.endswith("\\right)"))
        out += " \\cdot " if cdot else " "
        out += cur
    return out


def render_equation(params: dict, kind: str, power: Fraction | str,
                    comp: str, bound: Fraction | str) -> str:
    """Reproduce the site's /get_integral_image LaTeX string for this family."""
    # coef 原文进 spec：pi/catalan/zeta3 只用显示位，垃圾串也能渲染；
    # arctan/arccot 需要真值，不可解析时 None 自然崩进站端 500；pi_n 的
    # 渲染走 wire_pair 原始 (num, den)——奇偶、外层指数、ln 指数都用未约分值
    pv = wire_or(power)
    if kind == "pi_n":
        pnum, pden = wire_pair(power)
        odd, q = pnum % 2 == 0, Fraction(1)
    else:
        cfg = spec(kind, pv)
        odd, q = cfg["odd"], cfg["q"]
    m, n = int(params["m"]), int(params["n"])
    au, bu, u = int(params["au_val"]), int(params["bu_val"]), int(params["u_val"])

    btex = rat_tex(bound)  # raw 字符串原样回显数位（不约分）
    lhs = (f"{const_tex(kind, power)} - {btex}" if comp == ">"
           else f"{btex} - {const_tex(kind, power)}")
    eq_sep = " = "
    if kind == "pi_n" and pden != 1:
        # fractional power p/q: LHS shows (pi^{p/q})^q - (bound)^q literally;
        # '>' uses the site's tight spacing, '<' the normal one
        c = const_tex(kind, power) + f"^{{{pden}}}"
        b_ = f"\\left({btex}\\right)^{{{pden}}}"
        lhs, eq_sep = (f"{c}- {b_}", "= ") if comp == ">" else (f"{b_} - {c}", " = ")

    # scale t so the denominator u'·(1+q^2 x^2) has integer coefficients
    s = q.denominator
    t = s * s // gcd(s * s, u)

    x = sp.symbols("x")
    qq = sp.Rational(q.numerator, q.denominator)
    num = numerator(m, n, odd, au, bu, t)
    den = u * t * (qq * qq * x ** 2 + 1)

    # the site feeds the assembled fraction through sympy: literally identical
    # num/den factors auto-cancel (pi 0 < 1 -> 1) but a merely proportional
    # numerator stays a fraction (zeta3 0 < 1 keeps x(4x²+4)/(x²+1))
    ratio = num / den
    if ratio.as_numer_denom()[1] == 1:
        body = sp.latex(ratio)
    else:
        body = f"\\frac{{{factors_tex(num)}}}{{{sp.latex(den)}}}"

    if kind == "pi_n":
        r = pnum - 1  # raw numerator, unreduced — '^-4' prints unbraced
        ln = "" if r == 0 else "(\\ln(1/x))" if r == 1 else f"(\\ln(1/x))^{r}"
    elif kind == "catalan":
        ln = "\\ln(1/x)"
    elif kind == "zeta3":
        ln = "\\ln^2(x)"
    else:
        ln = ""
    tail = "  \\mathrm{d} x > 0" if (kind == "pi" and comp == ">") else " \\mathrm{d} x > 0"
    return f"{lhs}{eq_sep}\\int_0^1 {body}{ln}{tail}"
