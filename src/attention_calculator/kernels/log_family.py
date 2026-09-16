"""log kernel family: ``ln_q``, ``ln_q_square``, ``artanh_q``, ``arcoth_q``.

All four types share the basis ``x^m(1-x)^n`` on [0,1] and a denominator power
``(1+(q~-1)x)^s`` with ``s = max(m, n, 1)`` — a fixed site convention, never
searched (92/92 golden fits, docs/log-notes.md):

- ``ln_q``:        ``(a+bx) / (1+cx)^s``                 moments in span{ln q, 1}
- ``ln_q_square``: ``(a+bx+cx^2) ln(1+cx) / (1+cx)^s``   span{ln^2 q, ln q, 1}
- ``artanh q = ln(q~)/2`` with ``q~ = (1+q)/(1-q)``, ``q in (0,1)``;
- ``arcoth q = ln(q~)/2`` with ``q~ = (q+1)/(q-1)``, ``q > 1``.

For the artanh/arcoth pair the site solves the reduced ln system with the
ln-coefficient target halved (so the reported a, b already carry the 1/2) and
smuggles ``q~`` to the renderer through ``c_val``.

(m, n) traversal is the site's parity rule: ``m+n`` ascending, ``|m-n|``
ascending, and inside each diagonal pair the smaller-m candidate first when
``m+n`` is odd but the larger-m one first when ``m+n`` is even — e.g. sum 2
yields (1,1), (2,0), (0,2) while sum 3 yields (1,2), (2,1), (0,3), (3,0).
This is engine.mn_order's site-wide rule — required e.g. for
``ln^2(3/2) < 17/100`` where the site picks (2,0) over the also-feasible (0,2).
"""

from fractions import Fraction
from math import comb, gcd

import sympy as sp

from ..engine import InternalError, WrongDirection, mn_order, search
from ..integrand import constant_mpf
from ..moment import Moment, combine
from ..render import emit, join_cdot, lhs_tex, rat_tex, wire_pair

LIMIT = 10


def iota(c: Fraction, s: int, k: int) -> Moment:
    """Moment of ``x^k/(1+cx)^s``: ``{"ln": a, "1": b}`` means a*ln(1+c)+b.

    x^k = c^{-k} sum_r C(k,r)(-1)^{k-r}(1+cx)^r turns the integral into the
    power sum J_p = ((1+c)^{p+1}-1)/(c(p+1)) with J_{-1} = ln(1+c)/c.
    """
    ln = rat = Fraction(0)
    for r in range(k + 1):
        w = Fraction((-1) ** (k - r) * comb(k, r))
        p = r - s
        if p == -1:
            ln += w / c
        else:
            rat += w * ((1 + c) ** (p + 1) - 1) / (c * (p + 1))
    ck = c**k
    return {sym: v / ck for sym, v in (("ln", ln), ("1", rat)) if v}


def kappa(c: Fraction, s: int, k: int) -> Moment:
    """Moment of ``x^k ln(1+cx)/(1+cx)^s`` over span{ln^2, ln, 1} of q = 1+c.

    L_p = (1/c) int_1^q u^p ln u du contributes q^{p+1}*ln q/(c(p+1)) to "ln"
    and -(q^{p+1}-1)/(c(p+1)^2) to "1"; L_{-1} = ln^2(q)/(2c) hits "ln2".
    """
    q = 1 + c
    l2 = l1 = rat = Fraction(0)
    for r in range(k + 1):
        w = Fraction((-1) ** (k - r) * comb(k, r))
        p = r - s
        if p == -1:
            l2 += w / (2 * c)
        else:
            qp = q ** (p + 1)
            l1 += w * qp / (c * (p + 1))
            rat -= w * (qp - 1) / (c * (p + 1) ** 2)
    ck = c**k
    return {sym: v / ck for sym, v in (("ln2", l2), ("ln", l1), ("1", rat)) if v}


def basis_moment(c: Fraction, s: int, m: int, n: int, j: int, square: bool) -> Moment:
    """Moment of ``x^{m+j}(1-x)^n [ln(1+cx)] / (1+cx)^s`` for unknown j."""
    ks = [(kappa if square else iota)(c, s, m + j + t) for t in range(n + 1)]
    return combine([Fraction((-1) ** t * comb(n, t)) for t in range(n + 1)], ks)


def qtilde(kind: str, power: Fraction) -> Fraction:
    """Denominator parameter q~: power itself for ln types, the reduced
    ln argument (1+q)/(1-q) resp. (q+1)/(q-1) for artanh/arcoth."""
    if kind == "artanh_q":
        return (1 + power) / (1 - power)
    if kind == "arcoth_q":
        return (power + 1) / (power - 1)
    return power


def check_input(power: Fraction, bound: Fraction, kind: str) -> None:
    """Reproduce the site's validation: bound format, then left-side format,
    then the per-type domain rule — each with the probed Chinese message."""
    if bound < 0:
        raise ValueError("右侧有理数格式无效")
    if power < 0:
        raise ValueError("左侧系数格式无效")
    if kind in ("ln_q", "ln_q_square") and power <= 1:
        raise ValueError("请在ln后输入一个大于1的数")
    if kind == "artanh_q" and (power.denominator == 1 or not 0 < power < 1):
        raise ValueError("请在输入一个在(0,1)内的分数，本情况不支持整数")
    if kind == "arcoth_q" and power <= 1:
        raise ValueError("请在输入一个大于1的数")


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """prove(kind, power, comp, bound) -> site /calculate shape."""
    check_input(power, bound, kind)
    if kind == "ln_q_square" and power in (Fraction(5), Fraction(7)):
        # 站端 bug：q=5/7 的求解过程必崩（500），但命题为假时仍先走常规
        # 的"方向反了"——golden 16/16 按命题真假分列。方向判定按站端 float64。
        diff = float(constant_mpf(kind, power)) - float(bound)
        if (diff < 0) if comp == ">" else (diff > 0):
            raise WrongDirection
        raise InternalError
    qt = qtilde(kind, power)
    c = qt - 1
    square = kind == "ln_q_square"
    plans = (
        (m, n, [basis_moment(c, max(m, n, 1), m, n, j, square) for j in range(3 if square else 2)])
        for m, n in mn_order(LIMIT)
    )
    sign = Fraction(1 if comp == ">" else -1)
    if square:
        target = {"ln2": sign, "1": -sign * bound}
    else:
        # artanh/arcoth solve ln q~ / 2: the halved target emits a, b halved
        coef = Fraction(1, 2) if kind != "ln_q" else Fraction(1)
        target = {"ln": sign * coef, "1": -sign * bound}
    solved = search(plans, target, True)

    # c_val doubles as the reduced q~ for artanh/arcoth (renderer needs it)
    c_val = None if square or kind == "ln_q" else str(qt)
    return emit(solved.m, solved.n, solved.coeffs, c_val=c_val)


def ln_bound_proof(q: Fraction, comp: str, bound: Fraction) -> dict:
    """Public entry for the gamma family: prove ``ln q ⋚ bound``."""
    return prove("ln_q", q, comp, bound)


# ------------------------------------------------------------------- rendering


def const_tex(kind: str, power: Fraction | str) -> str:
    """Left-side constant text: ``\\ln2``, ``\\ln^2\\dfrac{3}{2}``,
    ``\\mathrm{artanh}\\dfrac{1}{2}``, ``\\mathrm{arcoth}2``."""
    body = rat_tex(power)
    if kind == "ln_q":
        return "\\ln" + body
    if kind == "ln_q_square":
        return "\\ln^2" + body
    return ("\\mathrm{artanh}" if kind == "artanh_q" else "\\mathrm{arcoth}") + body


def poly_latex(au: int, bu: int, cu: int) -> str:
    """sympy-order LaTeX of ``au+bu*x+cu*x^2`` (desc degree, coeff-1 bare)."""
    x = sp.symbols("x")
    return sp.latex(au + bu * x + cu * x**2)


def numerator_latex(
    m: int, n: int, au: int, bu: int, cu: int, t: int, c: Fraction, square: bool
) -> str:
    """Site LaTeX of ``t * x^m (1-x)^n (au+bu*x+cu*x^2/u-free) [log]``.

    The scale ``t`` folds into the first surviving factor (x-part, then the
    (1-x)-part, then P). A monomial P merges into the x-part; P == 1-x merges
    into the basis power; a constant P prints as a bare scalar unless it is 1
    (then it vanishes) and unless it is the only factor. Parentheses wrap an
    Add factor iff the numerator has more than one factor; adjacent factors
    are joined by `` \\cdot `` iff the left ends in ``}`` and the right is a
    ``\\left(`` factor starting with a digit and ending in ``\\right)`` —
    i.e. a parenthesized group carrying no outer exponent.
    """
    p = [au, bu, cu]
    nz = [i for i, v in enumerate(p) if v]
    if nz == [0, 1] and au == 1 and bu == -1:
        n += 1  # P == 1-x folds into the basis power
        nz = []
    mono = nz[0] if len(nz) == 1 else -1  # single-term P merges into x^mono

    pieces = []  # (latex, is_add)
    if m or mono > 0:
        c1 = t * (p[mono] if mono > 0 else 1)
        e = m + max(mono, 0)
        body = "x" if e == 1 else f"x^{{{e}}}"
        pieces.append((body if c1 == 1 else f"{c1} {body}", False))
    if n:
        scaled = not pieces and t != 1  # t lands here only if x-part is absent
        if n == 1:
            pieces.append((f"{t} - {t} x" if scaled else "1 - x", True))
        else:
            base = f"\\left(1 - x\\right)^{{{n}}}"
            pieces.append((f"{t} {base}" if scaled else base, False))
    if len(nz) >= 2:
        pieces.append((poly_latex(*(v * (t if not pieces else 1) for v in p)), True))
    elif mono == 0:  # constant P prints as a bare scalar; 1 vanishes
        v = au * (t if not pieces else 1)
        if v != 1 or not pieces:
            pieces.append((str(v), False))
    if square:
        pieces.append((sp.latex(sp.log(1 + sp.symbols("x") * c)), False))

    wrap = len(pieces) > 1
    tex = [f"\\left({s}\\right)" if is_add and wrap else s for s, is_add in pieces]
    return join_cdot(tex)


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Rebuild the site's get_integral_image LaTeX for solved parameters."""
    m, n = int(params["m"]), int(params["n"])
    au, bu, cu = (int(params[k]) for k in ("au_val", "bu_val", "cu_val"))
    u = int(params["u_val"])
    x = sp.symbols("x")
    if kind in ("artanh_q", "arcoth_q"):
        c = Fraction(params["c_val"]) - 1  # q~ stashed by prove, already reduced
        d, e = c.numerator, c.denominator
    else:
        # ln types derive (q-1) = (d/e) from the raw, unreduced wire pair —
        # '4/2' and '2' build different denominators, and a '\frac' macro
        # parses unsigned (probe: '\frac{-4}{2}' -> (2x+2)^6 like '4/2')
        wn, wd = wire_pair(power)
        d, e = wn - wd, wd
        c = Fraction(d, e)
    s = max(m, n, 1)
    g = gcd(u, e**s)
    t = e**s // g
    up = u // g

    # denominator: s == 1 prints the expanded ``A x + B``; s >= 2 prints
    # ``u' (d x + e)^{s}`` with the u' factor omitted when it is 1; the linear
    # factor goes through sympy's ordering ('4 - 10 x' positive-first)
    if s == 1:
        den = sp.latex(up * d * x + up * e)
    else:
        den = f"\\left({sp.latex(d * x + e)}\\right)^{{{s}}}"
        if up != 1:
            den = f"{up} {den}"

    num = numerator_latex(m, n, au, bu, cu, t, c, kind == "ln_q_square")

    lhs = lhs_tex(const_tex(kind, power), rat_tex(bound), comp)
    return f"{lhs} = \\int_0^1 \\frac{{{num}}}{{{den}}} \\mathrm{{d}} x > 0"
