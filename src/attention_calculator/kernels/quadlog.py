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

from fractions import Fraction
from functools import cache
from math import comb, factorial, gcd, lcm

from ..engine import mn_order, search
from ..moment import Moment

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
        k, pd = power.numerator, power.denominator
        if not 1 <= k <= 10:
            raise ValueError("pi_n requires power numerator in [1, 10]")
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


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """Search (m, n) in the author's order; return the site's parameter dict."""
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

def rat_tex(v: Fraction) -> str:
    """'25' for integers, '\\dfrac{22}{7}' otherwise; sign stays in numerator."""
    return (str(v.numerator) if v.denominator == 1
            else f"\\dfrac{{{v.numerator}}}{{{v.denominator}}}")


def const_tex(kind: str, power: Fraction) -> str:
    """The target constant as printed on the equation's left-hand side."""
    if kind == "pi":
        return "\\pi" if power == 1 else rat_tex(power) + "\\pi"
    if kind == "pi_n":
        if power.denominator != 1:
            return f"\\left(\\pi^\\dfrac{{{power.numerator}}}{{{power.denominator}}}\\right)"
        k = power.numerator
        return f"\\pi^{{{k}}}" if k >= 10 else f"\\pi^{k}"
    if kind == "arctan_q":
        return "\\arctan" + rat_tex(power)
    if kind == "arccot_q":
        return "\\mathrm{arccot}" + rat_tex(power)
    if kind == "catalan":
        return "C" if power == 1 else rat_tex(power) + "C"
    return "\\zeta(3)" if power == 1 else rat_tex(power) + "\\zeta(3)"


def paren_pow(inner: str, e: int) -> str:
    """'\\left(inner\\right)^{e}', exponent omitted for e == 1."""
    return f"\\left({inner}\\right)" if e == 1 else f"\\left({inner}\\right)^{{{e}}}"


def poly_tex(au: int, bu: int) -> str:
    """Inner string of a + b·x^2: positive term first, x^2 term first on ties.

    Leading negative renders as '- 2 x^{2}'; coefficient 1 on x^2 is omitted.
    """
    first_x2 = bu > 0 or au <= 0  # both-nonpositive also lists the x^2 term first
    terms = [(bu, True), (au, False)] if first_x2 else [(au, False), (bu, True)]
    out = ""
    for i, (v, is_x2) in enumerate(terms):
        body = ("x^{2}" if abs(v) == 1 else f"{abs(v)} x^{{2}}") if is_x2 else str(abs(v))
        if i == 0:
            out = ("- " + body) if v < 0 else body
        else:
            out += (" + " if v > 0 else " - ") + body
    return out


def numerator_tex(m: int, n: int, odd: bool, au: int, bu: int, t: int) -> str:
    """Rendered numerator of t·x^{2m+odd}(1-x^2)^n(au+bu·x^2)/u.

    The site canonicalizes monomial factors: a one-term P merges into the x-part
    (a as scalar, b·x^2 as extra degree), and P = 1-x^2 (i.e. au==1, bu==-1)
    merges into the basis power. The scale factor t lands on the first factor.
    Exactly one '\\cdot' is printed; its position is the observed truth table.
    """
    p = 2 * m + odd
    if au == 1 and bu == -1:  # folded: t·x^p·(1-x^2)^{n+1}
        e = n + 1
        xp = f"x^{{{p}}}" if t == 1 else f"{t} x^{{{p}}}"
        b = paren_pow("1 - x^{2}", e) if t == 1 or p else (
            f"\\left({t} - {t} x^{{2}}\\right)" if e == 1
            else f"{t} \\left(1 - x^{{2}}\\right)^{{{e}}}")
        if p == 0:
            return "1 - x^{2}" if (e == 1 and t == 1) else (
                f"{t} - {t} x^{{2}}" if e == 1 else b)
        return xp + (" \\cdot " if e == 1 else " ") + b
    # single-term P merges into the x-part
    mono_b = bu if au == 0 else 0
    p += 2 if mono_b else 0
    c1 = t * (au if bu == 0 else mono_b or 1)  # scalar on the x-part
    pieces = []
    if p:
        pieces.append(f"x^{{{p}}}" if c1 == 1 else f"{c1} x^{{{p}}}")
    if n:
        if p or t == 1:
            pieces.append(paren_pow("1 - x^{2}", n))
        else:
            pieces.append(f"\\left({t} - {t} x^{{2}}\\right)" if n == 1
                          else f"{t} \\left(1 - x^{{2}}\\right)^{{{n}}}")
    two_term = au != 0 and bu != 0
    if two_term:
        if not pieces:
            # lone (possibly t-scaled) polynomial, printed expanded without parens
            return poly_tex(au * t, bu * t)
        pieces.append(f"\\left({poly_tex(au, bu)}\\right)")
    if len(pieces) == 1:
        return pieces[0]
    # joins: cdot before a degree-1 basis factor or before P; plain space otherwise,
    # except a lone unscaled basis factor never takes a cdot.
    if len(pieces) == 2:
        first, second = pieces
        if not two_term:  # [x, B] merged-monomial/constant pair
            sep = " \\cdot " if n == 1 else " "
        elif p:          # [x, P]
            sep = " \\cdot "
        else:            # [B, P]: cdot only when B is a scaled power (t>1, n>1)
            sep = " \\cdot " if t != 1 and n > 1 else " "
        return first + sep + second
    # three pieces [x, B, P]
    x, b, pp = pieces
    j1 = " \\cdot " if n == 1 else " "
    j2 = " " if n == 1 else " \\cdot "
    return x + j1 + b + j2 + pp


def render_equation(params: dict, kind: str, power: Fraction,
                    comp: str, bound: Fraction) -> str:
    """Reproduce the site's /get_integral_image LaTeX string for this family."""
    cfg = spec(kind, power)
    m, n = int(params["m"]), int(params["n"])
    au, bu, u = int(params["au_val"]), int(params["bu_val"]), int(params["u_val"])
    q: Fraction = cfg["q"]

    lhs = (f"{const_tex(kind, power)} - {rat_tex(bound)}" if comp == ">"
           else f"{rat_tex(bound)} - {const_tex(kind, power)}")
    eq_sep = " = "
    if kind == "pi_n" and power.denominator != 1:
        # fractional power p/q: LHS shows (pi^{p/q})^q - (bound)^q literally;
        # '>' uses the site's tight spacing, '<' the normal one
        c = const_tex(kind, power) + f"^{{{power.denominator}}}"
        b_ = f"\\left({rat_tex(bound)}\\right)^{{{power.denominator}}}"
        lhs, eq_sep = (f"{c}- {b_}", "= ") if comp == ">" else (f"{b_} - {c}", " = ")

    # scale t so the denominator u'·(1+q^2 x^2) has integer coefficients
    s = q.denominator
    t = s * s // gcd(s * s, u)
    den = f"x^{{2}} + {u * t}" if q * q * u * t == 1 else f"{u * t * q * q} x^{{2}} + {u * t}"

    num = numerator_tex(m, n, cfg["odd"], au, bu, t)

    if kind == "pi_n":
        r = power.numerator - 1
        ln = "" if r == 0 else "(\\ln(1/x))" if r == 1 else f"(\\ln(1/x))^{r}"
    elif kind == "catalan":
        ln = "\\ln(1/x)"
    elif kind == "zeta3":
        ln = "\\ln^2(x)"
    else:
        ln = ""
    tail = "  \\mathrm{d} x > 0" if (kind == "pi" and comp == ">") else " \\mathrm{d} x > 0"
    return f"{lhs}{eq_sep}\\int_0^1 \\frac{{{num}}}{{{den}}}{ln}{tail}"
