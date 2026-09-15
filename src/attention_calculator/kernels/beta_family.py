"""beta kernel family: ``golden``, ``varpi``, ``gauss``.

- ``golden``: ``x^m (1-x)^n (a+bx) sqrt(x+4)`` on [0,1], moments live in
  span{phi, 1}; (m, n) traversed in engine.mn_order.
- ``varpi`` / ``gauss``: ``x^{4m+res} (1-x) (a+bx^4) / (pi^e sqrt(1-x^4))``.
  The basis polynomial (1-x) is fixed; the returned ``n`` is just a direction
  flag (0 for '>', 1 for '<'), not an exponent. Only m = 0..10 is searched.

Residues and outer factors per (kind, comp):

    varpi '>' : res=1, /pi   -> moments in span{1, varpi^{-1}}
    varpi '<' : res=3        -> moments in span{1, varpi}
    gauss '>' : res=0, /pi   -> moments in span{gauss, 1}
    gauss '<' : res=2        -> moments in span{gauss^{-1}, 1}

Moment ingredients. J_k = B((k+1)/4, 1/2)/4 satisfies J_{k+4}/J_k =
(k+1)/(k+3), giving four interleaved sequences:

    J_{4j}   = (varpi/2) p_j      p_{j+1} = p_j (4j+1)/(4j+3)
    J_{4j+1} = (pi/4) c_j         c_{j+1} = c_j (2j+1)/(2j+2)
    J_{4j+2} = (gauss^{-1}/2) q_j q_{j+1} = q_j (4j+3)/(4j+5)
    J_{4j+3} = d_j / 2            d_{j+1} = d_j (2j+2)/(2j+3)

and dividing a moment by pi maps symbols varpi->gauss, pi->1,
gauss_inv->varpi_inv (since gauss = varpi/pi).

For golden, S_k = int_0^1 x^k sqrt(x+4) dx = A_k + B_k sqrt(5) with
    S_k = sum_j C(k,j)(-4)^{k-j} * 2(5^{j+1} sqrt(5) - 8*4^j)/(2j+3)
mapped to the phi-basis via sqrt(5) = 2 phi - 1.
"""

from fractions import Fraction
from math import comb, lcm

import sympy as sp

from ..engine import mn_order, search
from ..moment import Moment, combine
from ..render import rat_tex, wire_fraction

LIMIT = 10


# --------------------------------------------------------------------- golden


def golden_moment(k: int) -> Moment:
    """Moment of ``x^k sqrt(x+4)`` over {"1", "phi"} (via sqrt(5) = 2 phi - 1)."""
    a = b = Fraction(0)
    for j in range(k + 1):
        w = Fraction(2 * comb(k, j) * (-4) ** (k - j), 2 * j + 3)
        a -= w * 8 * 4 ** j
        b += w * 5 ** (j + 1)
    return {"1": a - b, "phi": 2 * b}


def golden_basis(m: int, n: int, i: int) -> Moment:
    """Moment of ``x^{m+i} (1-x)^n sqrt(x+4)`` for unknown slot i in (a, b)."""
    ms = [golden_moment(m + i + t) for t in range(n + 1)]
    return combine([Fraction((-1) ** t * comb(n, t)) for t in range(n + 1)], ms)


# --------------------------------------------------------------- varpi/gauss


def lemniscate_table(n: int):
    """p, c, q, d sequences of length n+1 for the J_k beta moments."""
    p, c, q, d = [Fraction(1)], [Fraction(1)], [Fraction(1)], [Fraction(1)]
    for j in range(n):
        p.append(p[-1] * Fraction(4 * j + 1, 4 * j + 3))
        c.append(c[-1] * Fraction(2 * j + 1, 2 * j + 2))
        q.append(q[-1] * Fraction(4 * j + 3, 4 * j + 5))
        d.append(d[-1] * Fraction(2 * j + 2, 2 * j + 3))
    return p, c, q, d


def j_moment(k: int, p, c, q, d) -> Moment:
    """J_k = int_0^1 x^k / sqrt(1-x^4) dx over the four symbol classes."""
    j, rem = divmod(k, 4)
    if rem == 0:
        return {"varpi": p[j] / 2}
    if rem == 1:
        return {"pi": c[j] / 4}
    if rem == 2:
        return {"gauss_inv": q[j] / 2}
    return {"1": d[j] / 2}


# dividing a J-moment by pi renames its symbol (varpi/pi = gauss etc.)
DIV_PI = {"varpi": "gauss", "pi": "1", "gauss_inv": "varpi_inv"}


def lemniscate_basis(kind: str, comp: str, m: int, i: int) -> Moment:
    """Moment of ``x^{4m+res+4i} (1-x) / (pi^e sqrt(1-x^4))``, i in (a, b)."""
    res = {"varpi": 1, "gauss": 0}[kind] if comp == ">" else {
        "varpi": 3, "gauss": 2}[kind]
    p, c, q, d = lemniscate_table(m + 2)  # j index reaches m+2 for slot b
    base = combine([Fraction(1), Fraction(-1)],
                   [j_moment(4 * m + res + 4 * i, p, c, q, d),
                    j_moment(4 * m + res + 4 * i + 1, p, c, q, d)])
    if comp == ">":  # the '>' kernels carry an outer 1/pi
        base = {DIV_PI[sym]: v for sym, v in base.items()}
    return base


def lemniscate_target(kind: str, comp: str, power: Fraction,
                      bound: Fraction) -> Moment:
    """Target vector: ``power*C - bound`` arranged per the family's identity."""
    if kind == "varpi":
        return ({"1": power, "varpi_inv": -bound} if comp == ">"
                else {"1": bound, "varpi": -power})
    return ({"gauss": power, "1": -bound} if comp == ">"
            else {"gauss_inv": bound, "1": -power})


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """prove(kind, power, comp, bound) -> site /calculate shape."""
    if kind == "golden":
        plans = ((m, n, [golden_basis(m, n, 0), golden_basis(m, n, 1)])
                 for m, n in mn_order(LIMIT))
        sign = Fraction(1 if comp == ">" else -1)
        target = {"phi": sign * power, "1": -sign * bound}
    else:
        # n is a direction flag here: 0 proves '>', 1 proves '<'
        flag = 0 if comp == ">" else 1
        plans = ((m, flag, [lemniscate_basis(kind, comp, m, 0),
                            lemniscate_basis(kind, comp, m, 1)])
                 for m in range(LIMIT + 1))
        target = lemniscate_target(kind, comp, power, bound)
    solved = search(plans, target, True)

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


# ------------------------------------------------------------------- rendering


def join_cdot(pieces: list[str]) -> str:
    """Join factors the way the site's sympy does: ' \\cdot ' iff the left piece
    ends in '}' and the right is a '\\left(<digit>...\\right)' group."""
    out = pieces[0]
    for prev, cur in zip(pieces, pieces[1:]):
        cdot = (prev.endswith("}") and cur.startswith("\\left(")
                and cur[6].isdigit() and cur.endswith("\\right)"))
        out += " \\cdot " if cdot else " "
        out += cur
    return out


def poly1_tex(au: int, bu: int) -> str:
    """sympy-order text of au + bu*x (positive term first)."""
    return sp.latex(au + bu * sp.symbols("x"))


def golden_numerator(m: int, n: int, au: int, bu: int) -> str:
    """Numerator of x^m (1-x)^n P sqrt(x+4); P sits before sqrt when bu < 0."""
    pieces = []
    if m:
        pieces.append("x" if m == 1 else f"x^{{{m}}}")
    if n:
        pieces.append("\\left(1 - x\\right)" if n == 1
                      else f"\\left(1 - x\\right)^{{{n}}}")
    p = f"\\left({poly1_tex(au, bu)}\\right)"
    if bu < 0:
        pieces += [p, "\\sqrt{x + 4}"]
    else:
        pieces += ["\\sqrt{x + 4}", p]
    return join_cdot(pieces)


def lemniscate_numerator(e: int, au: int, bu: int, u: int) -> str:
    """The x^e (au + bu x^4) / u block preceding (1-x)/sqrt(1-x^4).

    e == 0 falls back to raw sympy (which splits ``(au+bu x^4)/u`` into two
    fractions, matching the site); a monomial P merges into the x-part.
    """
    x = sp.symbols("x")
    if e == 0:
        return sp.latex(sp.Rational(au, u) + sp.Rational(bu, u) * x ** 4)
    if au == 0 or bu == 0:
        coef, deg = (bu, e + 4) if au == 0 else (au, e)
        return sp.latex(sp.Rational(coef, u) * x ** deg)
    xp = "x" if e == 1 else f"x^{{{e}}}"
    num = join_cdot([xp, f"\\left({sp.latex(au + bu * x ** 4)}\\right)"])
    return num if u == 1 else f"\\frac{{{num}}}{{{u}}}"


def const_tex(kind: str, power: Fraction | str) -> str:
    """The constant side as printed: '2\\phi', '\\varpi', '3G', ..."""
    pre = "" if wire_fraction(power) == 1 else rat_tex(power)
    return pre + {"golden": "\\phi", "varpi": "\\varpi", "gauss": "G"}[kind]


def render_equation(params: dict, kind: str, power: Fraction | str,
                    comp: str, bound: Fraction | str) -> str:
    """Rebuild the site's get_integral_image LaTeX for solved parameters."""
    m = int(params["m"])
    au, bu, u = (int(params[k]) for k in ("au_val", "bu_val", "u_val"))
    btex = rat_tex(bound)
    bound_v = wire_fraction(bound)

    if kind == "golden":
        lhs = (f"{const_tex(kind, power)} - {btex}" if comp == ">"
               else f"{btex} - {const_tex(kind, power)}")
        num = golden_numerator(m, int(params["n"]), au, bu)
        return (f"{lhs} = \\int_0^1 \\frac{{{num}}}{{{u}}}"
                " \\mathrm{d} x > 0")

    # varpi '>' shows power - bound*varpi^{-1}; gauss '<' shows bound*G^{-1} -
    # power, both printed with the coefficient on the matching side. gauss
    # drops the spaces around '-'.
    ctex = const_tex(kind, power)
    if kind == "varpi":
        inv = "\\varpi^{-1}" if bound_v == 1 else f"{btex}\\varpi^{{-1}}"
        lhs = (f"{rat_tex(power)} - {inv}" if comp == ">"
               else f"{btex} - {ctex}")
    else:
        ginv = "G^{-1}" if bound_v == 1 else f"{btex}G^{{-1}}"
        lhs = (f"{ctex}-{btex}" if comp == ">" else f"{ginv}-{rat_tex(power)}")

    res = {"varpi": 1, "gauss": 0}[kind] if comp == ">" else {
        "varpi": 3, "gauss": 2}[kind]
    num = lemniscate_numerator(4 * m + res, au, bu, u)
    tail = ("\\dfrac{(1-x)}{\\pi\\sqrt{1-x^4}}" if comp == ">"
            else "\\dfrac{(1-x)}{\\sqrt{1-x^4}}")
    gap = "  " if comp == ">" else " "
    return f"{lhs} = \\int_0^1{gap}{num}{tail} \\mathrm{{d}} x > 0"
