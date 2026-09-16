"""``ln_q_cube`` / ``ln_q_quad`` kernels: prove ``(ln q)^r ⋚ bound`` for
r in {3, 4} and rational q > 1.

Exact-mode-only types — zhuyidao.net has neither ln^3 nor ln^4 kind, so
nothing here is constrained by site parity; registration into solve.FAMILY
/ kernels.EXACT_TYPES is merge wiring.

Kernel ``ln^{r-1}(1+cx)/(1+cx)^s`` with c = q-1 and the ln-family
denominator rule s = max(m, n, 1) over basis ``x^m(1-x)^n``: substituting
u = 1+cx turns each moment into a sum of ``int_1^q u^p ln^{r-1} u du``
terms, where p = -1 gives ``ln^r(q)/r`` — the target constant — and
p != -1 terms fill the lower tower {ln^{r-1} q, ..., ln q, 1}. The kernel
carries one power of log below the target, the same mechanism as ln_q
(ln^0 kernel -> ln^1) and ln_q_square (ln^1 kernel -> ln^2). The moment
space is (r+1)-dimensional, so P carries r+1 unknowns — a degree r-1
polynomial would leave the system overdetermined. Derivations and numeric
checks: docs/2026-09-16-ln-cube-derivation.md (r=3) and
docs/2026-09-16-ln-quad-impl-notes.md (r=4).

The fourth and fifth coefficients extend the nine site param slots with
``d_val``/``du_val`` resp. ``e_val``/``eu_val``; ``u_val`` stays the common
denominator.

Nonnegativity on [0,1] is decided exactly per degree: the cubic P uses a
QQ(sqrt D) critical-value rule (interior minima at roots of P' = 3dx^2+2cx+b
collapse to ``A ∓ B0*sqrt(D)``, a quadratic-field sign test); the quartic P
uses the Sturm machinery in ``quartic_nonneg``.
"""

from fractions import Fraction
from math import comb, factorial, gcd, lcm

import sympy as sp

from ..engine import NoSolution, WrongDirection, mn_order, poly_nonneg, solve_moment
from ..moment import Moment, combine
from ..render import lhs_tex, rat_tex, wire_pair

LIMIT = 10

x_sym = sp.symbols("x")


def ln_sym(j: int) -> str:
    """Moment symbol for ln^j q: "1", "ln", "ln2", ..."""
    return "1" if j == 0 else "ln" if j == 1 else f"ln{j}"


def mu(c: Fraction, s: int, k: int, lp: int) -> Moment:
    """Moment of ``x^k ln^lp(1+cx)/(1+cx)^s`` over span{ln^{lp+1} q, ..., 1}.

    u = 1+cx turns x^k into c^{-k} * sum_r C(k,r)(-1)^{k-r} u^r and dx into
    du/c, leaving ``J_p = int_1^q u^p ln^lp u du`` with q = 1+c: J_{-1} =
    ln^{lp+1}(q)/(lp+1) hits the top symbol; for p != -1, iterated IBP gives
    ``J_p = q^{p+1} * sum_{i=0}^{lp} (-1)^i lp!/(lp-i)! * ln^{lp-i}q/(p+1)^{i+1}
    - (-1)^lp lp!/(p+1)^{lp+1}`` — the u=1 boundary touches only the ln^0
    term, flipping its q^{p+1} factor to q^{p+1}-1.
    """
    q = 1 + c
    acc: Moment = {}
    for r in range(k + 1):
        w = Fraction((-1) ** (k - r) * comb(k, r))
        p = r - s
        if p == -1:
            acc[ln_sym(lp + 1)] = acc.get(ln_sym(lp + 1), Fraction(0)) + w / ((lp + 1) * c)
            continue
        pp = p + 1
        qp = q**pp
        for i in range(lp + 1):
            coef = Fraction((-1) ** i * factorial(lp) // factorial(lp - i), pp ** (i + 1))
            key = ln_sym(lp - i)
            acc[key] = acc.get(key, Fraction(0)) + w * coef * (qp if i < lp else qp - 1) / c
    ck = c**k
    return {sym: v / ck for sym, v in acc.items() if v}


def basis_moment(c: Fraction, s: int, m: int, n: int, j: int, lp: int) -> Moment:
    """Moment of ``x^{m+j}(1-x)^n ln^lp(1+cx)/(1+cx)^s`` for unknown j."""
    ks = [mu(c, s, m + j + t, lp) for t in range(n + 1)]
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


def quartic_nonneg(coeffs: list[Fraction]) -> bool:
    """Exact ``a+bx+cx^2+dx^3+ex^4 >= 0 on [0,1]`` test (Sturm root count).

    P >= 0 on [0,1] iff the lowest-degree nonzero coefficient is positive
    — the sign of P just right of 0, so P(0) >= 0 and P(1) >= 0 come for
    free — and no root in the OPEN interval has odd multiplicity (odd roots
    are exactly the sign changes; even roots only touch). The odd part is
    the product of the odd-power squarefree factors of P after stripping
    the endpoint factors x and x-1, whose roots are admissible regardless
    of parity; with no roots left at 0 or 1, sympy's exact rational Sturm
    ``count_roots(0, 1)`` is unambiguous interior counting. All steps are
    exact over QQ.
    """
    first = next((c for c in coeffs if c != 0), None)
    if first is None:
        return True  # P ≡ 0
    if first < 0:
        return False
    p = sp.Poly([sp.Rational(c.numerator, c.denominator) for c in reversed(coeffs)], x_sym)
    for pt in (0, 1):
        while p.eval(x_sym, pt) == 0:
            p = p.exquo(sp.Poly(x_sym - pt, x_sym))
    odd = sp.Poly(1, x_sym)
    for f, e in p.sqf_list()[1]:
        if e % 2:
            odd *= f
    return odd.count_roots(0, 1) == 0


def quartic_nonpos(coeffs: list[Fraction]) -> bool:
    """Sign test for the quartic being <= 0 on [0,1]."""
    return quartic_nonneg([-v for v in coeffs])


def check_input(power: Fraction, bound: Fraction) -> None:
    """ln-family validation order: bound format, then the q > 1 domain rule."""
    if bound < 0:
        raise ValueError("右侧有理数格式无效")
    if power <= 1:
        raise ValueError("请在ln后输入一个大于1的数")


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("ln_q_cube"/"ln_q_quad", q, comp, bound) -> /calculate dict.

    The scan mirrors engine.search with nonneg=True but swaps the deg<=2
    sign test for the degree-matching exact rule; a uniformly non-positive
    solution certifies the opposite inequality (WrongDirection), and
    exhaustion is honest NoSolution. ``exact`` is accepted for signature
    compatibility — the types only exist in exact mode, and the math is
    identical anyway.
    """
    check_input(power, bound)
    c = power - 1
    sign = Fraction(1 if comp == ">" else -1)
    if kind == "ln_q_quad":
        lp, nonneg, nonpos = 3, quartic_nonneg, quartic_nonpos
    else:
        lp, nonneg, nonpos = 2, cubic_nonneg, cubic_nonpos
    target = {ln_sym(lp + 1): sign, "1": -sign * bound}
    for m, n in mn_order(LIMIT):
        basis = [basis_moment(c, max(m, n, 1), m, n, j, lp) for j in range(lp + 2)]
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            continue
        if nonneg(coeffs):
            return emit(m, n, coeffs)
        if nonpos(coeffs):
            raise WrongDirection
    raise NoSolution


def emit(m: int, n: int, coeffs: list[Fraction]) -> dict:
    """Assemble the /calculate payload for the four/five-coefficient P.

    Same field semantics as render.emit extended by ``d_val``/``du_val``
    (x^3) and ``e_val``/``eu_val`` (x^4); u_val = lcm of all denominators.
    """
    u = lcm(*(c.denominator for c in coeffs))
    params = {
        "m": m,
        "n": n,
        "a_val": str(coeffs[0]),
        "b_val": str(coeffs[1]),
        "c_val": str(coeffs[2]),
        "d_val": str(coeffs[3]),
        "au_val": str(coeffs[0] * u),
        "bu_val": str(coeffs[1] * u),
        "cu_val": str(coeffs[2] * u),
        "du_val": str(coeffs[3] * u),
        "u_val": str(u),
    }
    if len(coeffs) > 4:
        params["e_val"] = str(coeffs[4])
        params["eu_val"] = str(coeffs[4] * u)
    params["unified_form"] = {}
    sol = f"a = {coeffs[0]}, b = {coeffs[1]}, c= {coeffs[2]}, d= {coeffs[3]}"
    return {
        "parameters": params,
        "solution": sol + (f", e= {coeffs[4]}" if len(coeffs) > 4 else ""),
    }


# ------------------------------------------------------------------- rendering


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX in the ln-family layout.

    No site counterpart exists, so the shape follows log_family: LHS
    ``\\ln^r q ⋚ bound``, integrand ``x^m(1-x)^n P(x) ln^{r-1}(1+cx) /
    (1+cx)^s`` printed with the site scaling — numerator
    ``t*x^m(1-x)^n(au+bu*x+cu*x^2+du*x^3[+eu*x^4]) ln^{r-1}(...)``,
    denominator ``u'(dx+e)^s`` where t/u' = e^s/u after cancelling
    g = gcd(u, e^s).
    """
    lp = 3 if kind == "ln_q_quad" else 2
    m, n = int(params["m"]), int(params["n"])
    cus = [int(params[k]) for k in ("au_val", "bu_val", "cu_val", "du_val")]
    if lp == 3:
        cus.append(int(params["eu_val"]))
    u = int(params["u_val"])
    wn, wd = wire_pair(power)
    d, e = wn - wd, wd
    c = Fraction(d, e)
    s = max(m, n, 1)
    g = gcd(u, e**s)
    t, up = e**s // g, u // g
    poly = sum(v * x_sym**i for i, v in enumerate(cus))
    num = t * x_sym**m * (1 - x_sym) ** n * poly * sp.log(1 + c * x_sym) ** lp
    den = up * (d * x_sym + e) ** s
    frac = sp.latex(num / den).replace("\\log", "\\ln")
    lhs = lhs_tex(f"\\ln^{{{lp + 1}}}{rat_tex(power)}", rat_tex(bound), comp)
    return f"{lhs} = \\int_0^1 {frac} \\mathrm{{d}} x > 0"
