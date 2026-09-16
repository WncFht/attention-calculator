"""``ln_q_cube`` kernel: prove ``(ln q)^3 ⋚ bound`` for rational q > 1.

Exact-mode-only type — zhuyidao.net has no ln^3 kind, so nothing here is
constrained by site parity; registration into solve.FAMILY / kernels.TYPES
is merge wiring.

Kernel ``ln^2(1+cx)/(1+cx)^s`` with c = q-1 and the ln-family denominator
rule s = max(m, n, 1) over basis ``x^m(1-x)^n``: substituting u = 1+cx turns
each moment into a sum of ``int_1^q u^p ln^2 u du`` terms, where p = -1
gives ``ln^3(q)/3`` — the target constant — and p != -1 terms fill the
lower tower {ln^2 q, ln q, 1}. The kernel carries one power of log below
the target, the same mechanism as ln_q (ln^0 kernel -> ln^1) and
ln_q_square (ln^1 kernel -> ln^2). The moment space is 4-dimensional
{ln3, ln2, ln, 1}, so P = a+bx+cx^2+dx^3 carries four unknowns — a
three-coefficient P would leave the 4-symbol system overdetermined.
Derivation and numeric checks: docs/2026-09-16-ln-cube-derivation.md.

The fourth coefficient extends the nine site param slots with ``d_val`` /
``du_val``; ``u_val`` stays the common denominator of all four.

Nonnegativity of the cubic P on [0,1] is decided exactly in QQ(sqrt D):
interior minima sit at roots x* of P' = 3dx^2+2cx+b, and with
D = c^2-3bd the critical values collapse to ``A ∓ B0*sqrt(D)`` with
A = a + (2c^3-9bcd)/(27d^2) and B0 = 2D/(27d^2) — a quadratic-field sign
test, no Sturm machinery.
"""

from fractions import Fraction
from math import comb, gcd, lcm

import sympy as sp

from ..engine import NoSolution, WrongDirection, mn_order, poly_nonneg, solve_moment
from ..moment import Moment, combine
from ..render import lhs_tex, rat_tex, wire_pair

LIMIT = 10


def mu(c: Fraction, s: int, k: int) -> Moment:
    """Moment of ``x^k ln^2(1+cx)/(1+cx)^s`` over span{ln^3 q, ln^2 q, ln q, 1}.

    u = 1+cx turns x^k into c^{-k} * sum_r C(k,r)(-1)^{k-r} u^r and dx into
    du/c, leaving ``J_p = int_1^q u^p ln^2 u du`` with q = 1+c: J_{-1} =
    ln^3(q)/3 hits "ln3"; for p != -1, IBP on u^p ln^2 u gives
    ``q^{p+1}(ln^2 q/(p+1) - 2 ln q/(p+1)^2 + 2/(p+1)^3) - 2/(p+1)^3``.
    """
    q = 1 + c
    l3 = l2 = l1 = rat = Fraction(0)
    for r in range(k + 1):
        w = Fraction((-1) ** (k - r) * comb(k, r))
        p = r - s
        if p == -1:
            l3 += w / (3 * c)
        else:
            pp = p + 1
            qp = q**pp
            l2 += w * qp / (c * pp)
            l1 -= 2 * w * qp / (c * pp**2)
            rat += w * (2 * qp - 2) / (c * pp**3)
    ck = c**k
    return {sym: v / ck for sym, v in (("ln3", l3), ("ln2", l2), ("ln", l1), ("1", rat)) if v}


def basis_moment(c: Fraction, s: int, m: int, n: int, j: int) -> Moment:
    """Moment of ``x^{m+j}(1-x)^n ln^2(1+cx)/(1+cx)^s`` for unknown j."""
    ks = [mu(c, s, m + j + t) for t in range(n + 1)]
    return combine([Fraction((-1) ** t * comb(n, t)) for t in range(n + 1)], ks)


def qsign(alpha: Fraction, beta: Fraction, disc: Fraction) -> int:
    """Exact sign of ``alpha + beta*sqrt(disc)`` for rational alpha, beta and
    nonnegative rational disc.

    With alpha, beta of opposite sign the comparison reduces to
    alpha^2 vs beta^2*disc — which also returns 0 correctly when disc is a
    perfect square (alpha + beta*sqrt(disc) genuinely vanishes).
    """
    if beta == 0 or disc == 0:
        return (alpha > 0) - (alpha < 0)
    if alpha == 0:
        return (beta > 0) - (beta < 0)
    if (alpha > 0) == (beta > 0):
        return 1 if alpha > 0 else -1
    c = alpha * alpha - beta * beta * disc
    if c == 0:
        return 0
    if alpha > 0:  # beta < 0: value = |alpha| - |beta|*sqrt(disc)
        return 1 if c > 0 else -1
    return 1 if c < 0 else -1  # alpha < 0 < beta: value = |beta|*sqrt(disc) - |alpha|


def cubic_nonneg(coeffs: list[Fraction]) -> bool:
    """Exact ``a+bx+cx^2+dx^3 >= 0 on [0,1]`` test.

    The minimum is attained at an endpoint or an interior root x* of
    P' = 3dx^2+2cx+b. With D = c^2-3bd the critical values are
    ``P(x*±) = A ∓ B0*sqrt(D)`` where A = a + (2c^3-9bcd)/(27d^2) and
    B0 = 2D/(27d^2) — obtained by reducing P mod P'. Membership x* in (0,1)
    and the sign of P(x*) are both qsign calls in QQ(sqrt D).
    """
    a, b, c, d = coeffs
    if d == 0:
        return poly_nonneg([a, b, c])
    if a < 0 or a + b + c + d < 0:
        return False
    disc = c * c - 3 * b * d
    if disc < 0:
        return True  # P' never vanishes: monotone, endpoints already nonneg
    va = a + (2 * c**3 - 9 * b * c * d) / (27 * d**2)
    vb = 2 * disc / (27 * d**2)
    w = 1 if 3 * d > 0 else -1
    for sigma in (Fraction(1), Fraction(-1)):
        # x* = (-c + sigma*sqrt(D)) / (3d); interior iff 0 < x* < 1
        interior = qsign(-c, sigma, disc) * w > 0 and qsign(-c - 3 * d, sigma, disc) * w < 0
        if interior and qsign(va, -sigma * vb, disc) < 0:
            return False
    return True


def cubic_nonpos(coeffs: list[Fraction]) -> bool:
    """Sign test for the cubic being <= 0 on [0,1]."""
    return cubic_nonneg([-v for v in coeffs])


def check_input(power: Fraction, bound: Fraction) -> None:
    """ln-family validation order: bound format, then the q > 1 domain rule."""
    if bound < 0:
        raise ValueError("右侧有理数格式无效")
    if power <= 1:
        raise ValueError("请在ln后输入一个大于1的数")


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("ln_q_cube", q, comp, bound) -> /calculate-shaped dict.

    The scan mirrors engine.search with nonneg=True but swaps the deg<=2
    sign test for cubic_nonneg; a uniformly non-positive solution
    certifies the opposite inequality (WrongDirection), and exhaustion is
    honest NoSolution. ``exact`` is accepted for signature compatibility —
    the type only exists in exact mode, and the math is identical anyway.
    """
    check_input(power, bound)
    c = power - 1
    sign = Fraction(1 if comp == ">" else -1)
    target = {"ln3": sign, "1": -sign * bound}
    for m, n in mn_order(LIMIT):
        basis = [basis_moment(c, max(m, n, 1), m, n, j) for j in range(4)]
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            continue
        if cubic_nonneg(coeffs):
            return emit(m, n, coeffs)
        if cubic_nonpos(coeffs):
            raise WrongDirection
    raise NoSolution


def emit(m: int, n: int, coeffs: list[Fraction]) -> dict:
    """Assemble the /calculate payload for the four-coefficient P.

    Same field semantics as render.emit extended by ``d_val``/``du_val``
    for the x^3 coefficient; u_val = lcm of all four denominators.
    """
    a, b, c, d = coeffs
    u = lcm(a.denominator, b.denominator, c.denominator, d.denominator)
    params = {
        "m": m,
        "n": n,
        "a_val": str(a),
        "b_val": str(b),
        "c_val": str(c),
        "d_val": str(d),
        "au_val": str(a * u),
        "bu_val": str(b * u),
        "cu_val": str(c * u),
        "du_val": str(d * u),
        "u_val": str(u),
        "unified_form": {},
    }
    return {"parameters": params, "solution": f"a = {a}, b = {b}, c= {c}, d= {d}"}


# ------------------------------------------------------------------- rendering


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX in the ln-family layout.

    No site counterpart exists, so the shape follows log_family: LHS
    ``\\ln^3 q ⋚ bound``, integrand ``x^m(1-x)^n P_4(x) ln^2(1+cx) /
    (1+cx)^s`` printed with the site scaling — numerator
    ``t*x^m(1-x)^n(au+bu*x+cu*x^2+du*x^3) ln^2(...)``, denominator
    ``u'(dx+e)^s`` where t/u' = e^s/u after cancelling g = gcd(u, e^s).
    """
    m, n = int(params["m"]), int(params["n"])
    au, bu, cu, du = (int(params[k]) for k in ("au_val", "bu_val", "cu_val", "du_val"))
    u = int(params["u_val"])
    wn, wd = wire_pair(power)
    d, e = wn - wd, wd
    c = Fraction(d, e)
    s = max(m, n, 1)
    g = gcd(u, e**s)
    t, up = e**s // g, u // g
    x = sp.symbols("x")
    num = t * x**m * (1 - x) ** n * (au + bu * x + cu * x**2 + du * x**3) * sp.log(1 + c * x) ** 2
    den = up * (d * x + e) ** s
    frac = sp.latex(num / den).replace("\\log", "\\ln")
    lhs = lhs_tex(f"\\ln^{{3}}{rat_tex(power)}", rat_tex(bound), comp)
    return f"{lhs} = \\int_0^1 {frac} \\mathrm{{d}} x > 0"
