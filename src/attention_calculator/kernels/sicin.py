"""si_q / cin_q kernel family on [0, 1]: ``Si(q) ⋚ r`` and ``Cin(q) ⋚ r``.

Exact-mode-only types — zhuyidao.net has neither, so nothing here is
constrained by site parity. One moment machine serves both constants

    Si(q)  = ∫₀^q sin(t)/t dt,     Cin(q) = ∫₀^q (1−cos t)/t dt

through the shifted trig_q sequences T_k = ∫₀¹ x^{k−1} sin(qx) dx and
U_l = ∫₀¹ x^{l−1} (1−cos qx) dx:

    T_0 = Si(q),   T_k = S_{k−1}  for k ≥ 1   (trig_q.sin_moments shift)
    U_0 = Cin(q),  U_l = 1/l − C_l,  C_l = sin q/q − (l−1)/q · T_{l−1}

Every basis moment lands in the 4-symbol span {C_q, sin_q, cos_q, 1}, so
P = a+bx+cx²+dx³ carries four unknowns (a quadratic P would leave the
system overdetermined — same shape as ln_q_cube).

Two structural facts fix the family's form
(docs/2026-09-16-w3-research-sicin.md):

- F(0) pinning: C_q appears only in the l=i=j=0 basis term, so every
  solution has a = ±1. Under a nonnegative kernel a '<' claim would need
  P(0) = −1 — impossible — hence each direction drives its own kernel:
  '>' keeps the special function minus its Taylor lower truncation,
  '<' keeps the complementary upper truncation minus the function.
- The basis must be (1+x)^m(1−x)^n, not x^m(1−x)^n: m ≥ 1 zeroes the
  C_q row of every moment and the system turns singular.

The kernels are Taylor-remainder ladders indexed by t = 0, 1, 2, … —
global one-sided bounds obtained by iterating 1−cos u ≥ 0:

    si  '>': (sin qx − T_{4t+3}(qx))/x      si  '<': (T_{4t+1}(qx) − sin qx)/x
    cin '>': ((1−cos qx) − S_{2t+2}(qx))/x  cin '<': (S_{2t+1}(qx) − (1−cos qx))/x

each strictly positive on (0,1] with x = 0 the only zero (alternating
truncation bounds are strict for u > 0). The polynomial tail contributes
only pure rationals β(m,n,s) = ∫₀¹ x^s (1+x)^m (1−x)^n dx, so the ladder
moments stay in the same span. The bare kernels (sin qx)/x and
(1−cos qx)/x are dominated by t = 0 and unused — sin(qx)/x additionally
loses sign-definiteness for q > π, which would admit fake proofs.

Parity fold: Si is odd, Cin is even. With σ = sign(q) resp. +1, the claim
s·(C(q) − r) > 0 is solved on q̄ = |q| as target {C_q: sσ, 1: −sr} against
the kernel of direction sσ; σ·C(|q|) = C(q) makes the certified identity
the literal claim. q = 0 collapses the recurrences (1/q) and is rejected.

The cubic-P nonnegativity test reuses kernels.ln_pow.cubic_nonneg — the
exact QQ(√D) critical-value rule, no Sturm machinery. F(0) pinning forces
P(0) = +1 on every solution, so WrongDirection is structurally
unreachable: false claims and true-but-over-budget claims both exhaust as
NoSolution, and the mode=exact gate (certified_cmp) decides direction
upstream anyway.
"""

from fractions import Fraction
from math import comb, factorial, lcm

import sympy as sp

from ..engine import NoSolution, mn_order, solve_moment
from ..moment import Moment, add, combine, scale
from ..render import lhs_tex, rat_tex, wire_pair
from .ln_pow import cubic_nonneg
from .trig_q import sin_moments

LIMIT = 10
T_LEVELS = (0, 1, 2)

x = sp.symbols("x")

# LHS constant names in the rendered proof — no site counterparts exist.
NAME = {"si_q": "Si", "cin_q": "Cin"}


def t_moments(kmax: int, q: Fraction) -> list[Moment]:
    """T_k = ∫_0^1 x^{k−1} sin(qx) dx for k = 0..kmax, over {si_q, sin_q, cos_q, 1}.

    T_0 = Si(q) is the target symbol; T_k = S_{k−1} for k ≥ 1 is trig_q's
    sin_moments shifted one index, so only T_0 carries si_q.
    """
    return [{"si_q": Fraction(1)}, *sin_moments(kmax - 1, q)]


def u_moments(kmax: int, q: Fraction) -> list[Moment]:
    """U_l = ∫_0^1 x^{l−1} (1−cos qx) dx for l = 0..kmax, over {cin_q, sin_q, cos_q, 1}.

    U_0 = Cin(q); U_l = 1/l − C_l with C_l = ∫ x^{l−1} cos(qx) dx
    = sin q/q − (l−1)/q·T_{l−1}. At l = 1 the (l−1) factor swallows T_0's
    si_q symbol, leaving C_1 = sin q/q.
    """
    tmoms = t_moments(kmax - 1, q)
    out = [{"cin_q": Fraction(1)}]
    for ell in range(1, kmax + 1):
        cl = add({"sin_q": Fraction(1, q)}, scale(tmoms[ell - 1], Fraction(1 - ell, q)))
        out.append(add({"1": Fraction(1, ell)}, scale(cl, Fraction(-1))))
    return out


def beta(m: int, n: int, s: int) -> Fraction:
    """∫_0^1 x^s (1+x)^m (1−x)^n dx = Σ_{l,i} C(m,l) (−1)^i C(n,i) / (s+l+i+1)."""
    return sum(
        Fraction(comb(m, ell) * (-1) ** i * comb(n, i), s + ell + i + 1)
        for ell in range(m + 1)
        for i in range(n + 1)
    )


def trunc_moment(kind: str, kdir: int, t: int, q: Fraction, m: int, n: int, j: int) -> Fraction:
    """Moment of the ladder kernel's polynomial tail against x^j(1+x)^m(1−x)^n.

    si:  Σ_{r=0}^{R} (−1)^r q^{2r+1}/(2r+1)! · β(m,n,j+2r),   R = 2t+1 or 2t
    cin: Σ_{r=1}^{R} (−1)^{r+1} q^{2r}/(2r)! · β(m,n,j+2r−1), R = 2t+2 or 2t+1

    with R the deeper truncation on the '>' side (kdir > 0).
    """
    if kind == "si_q":
        rmax = 2 * t + (1 if kdir > 0 else 0)
        return sum(
            (-1) ** r * q ** (2 * r + 1) * beta(m, n, j + 2 * r) / factorial(2 * r + 1)
            for r in range(rmax + 1)
        )
    rmax = 2 * t + (2 if kdir > 0 else 1)
    return sum(
        (-1) ** (r + 1) * q ** (2 * r) * beta(m, n, j + 2 * r - 1) / factorial(2 * r)
        for r in range(1, rmax + 1)
    )


def basis_moment(
    kind: str, kdir: int, t: int, table: list[Moment], q: Fraction, m: int, n: int, j: int
) -> Moment:
    """∫_0^1 x^j (1+x)^m (1−x)^n K(x) dx for the level-t kernel of direction kdir.

    K = kdir·(g − tail)/x with g = sin qx (si) or 1−cos qx (cin): the moment
    is kdir·(special − tail_moment), special being the (1+x)^m(1−x)^n
    combination Σ_{l,i} C(m,l)(−1)^i C(n,i)·table[l+i+j] and tail_moment the
    pure-rational trunc_moment.
    """
    special = combine(
        [
            Fraction((-1) ** i * comb(m, ell) * comb(n, i))
            for ell in range(m + 1)
            for i in range(n + 1)
        ],
        [table[ell + i + j] for ell in range(m + 1) for i in range(n + 1)],
    )
    return scale(add(special, {"1": -trunc_moment(kind, kdir, t, q, m, n, j)}), kdir)


def check_input(q: Fraction) -> None:
    """q = 0 collapses the family: the T/U recurrences carry 1/q and
    Si(0) = Cin(0) = 0 is rational — outside the constant types."""
    if q == 0:
        raise ValueError("q不能为0：q=0 时该积分常数退化")


def solve(kind: str, q: Fraction, comp: str, bound: Fraction) -> tuple[int, int, int, list]:
    """Find (m, n, t, coeffs) whose integrand proves s·(C(q) − bound) >= 0.

    The machine runs on q̄ = |q| with kernel direction kdir = s·σ
    (σ = sign(q) for the odd Si, +1 for the even Cin): the solved identity
    kdir·C(q̄) − s·bound = ∫ F·K ≥ 0 is exactly the literal claim. P(0) = +1
    is forced by the C_q row, so a sign-definite solution is nonnegative or
    the candidate is skipped — nonpositivity cannot occur.
    """
    check_input(q)
    s = 1 if comp == ">" else -1
    sigma = (1 if q > 0 else -1) if kind == "si_q" else 1
    kdir = s * sigma
    target = {k: v for k, v in {kind: Fraction(kdir), "1": -s * bound}.items() if v}
    qq = abs(q)
    table = (t_moments if kind == "si_q" else u_moments)(2 * LIMIT + 3, qq)
    for m, n in mn_order(LIMIT):
        for t in T_LEVELS:
            basis = [basis_moment(kind, kdir, t, table, qq, m, n, j) for j in range(4)]
            try:
                coeffs = solve_moment(basis, target)
            except ValueError:
                continue
            if cubic_nonneg(coeffs):
                return m, n, t, coeffs
    raise NoSolution


def emit(m: int, n: int, t: int, coeffs: list[Fraction]) -> dict:
    """Assemble the /calculate payload for the four-coefficient P.

    Same field semantics as kernels.ln_pow.emit extended by ``t_val`` for
    the ladder level; u_val = lcm of all four denominators.
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
        "t_val": str(t),
        "unified_form": {},
    }
    return {"parameters": params, "solution": f"a = {a}, b = {b}, c= {c}, d= {d}"}


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("si_q", q, comp, bound) -> /calculate-shaped dict.

    Exact-mode-only family; ``exact`` is accepted for signature
    compatibility — there is no site path to reproduce.
    """
    m, n, t, coeffs = solve(kind, power, comp, bound)
    result = emit(m, n, t, coeffs)
    result["type"] = kind
    return result


# ------------------------------------------------------------------- rendering


def kernel_expr(kind: str, kdir: int, t: int, qq: sp.Rational) -> sp.Expr:
    """The level-t ladder kernel kdir·(g − tail)/x as a sympy expression."""
    if kind == "si_q":
        rmax = 2 * t + (1 if kdir > 0 else 0)
        tail = sum(
            (-1) ** r * (qq * x) ** (2 * r + 1) / sp.factorial(2 * r + 1) for r in range(rmax + 1)
        )
        g = sp.sin(qq * x)
    else:
        rmax = 2 * t + (2 if kdir > 0 else 1)
        tail = sum(
            (-1) ** (r + 1) * (qq * x) ** (2 * r) / sp.factorial(2 * r) for r in range(1, rmax + 1)
        )
        g = 1 - sp.cos(qq * x)
    return kdir * (g - tail) / x


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX: ``C(q) − r``/``r − C(q)`` = ∫ (1+x)^m(1−x)^n P K.

    The kernel prints on q̄ = |q| — the moments are computed there and the
    parity fold lives in the solved target, not the display.
    """
    m, n, t = int(params["m"]), int(params["n"]), int(params["t_val"])
    au, bu, cu, du = (int(params[k]) for k in ("au_val", "bu_val", "cu_val", "du_val"))
    u = int(params["u_val"])
    qn, qd = wire_pair(power)
    qf = Fraction(qn, qd)
    sigma = (1 if qf > 0 else -1) if kind == "si_q" else 1
    kdir = (1 if comp == ">" else -1) * sigma
    qq = sp.Rational(abs(qf).numerator, abs(qf).denominator)

    const = rf"\mathrm{{{NAME[kind]}}}\left({rat_tex(power)}\right)"
    lhs = lhs_tex(const, rat_tex(bound), comp)
    kern = kernel_expr(kind, kdir, t, qq)
    body = (1 + x) ** m * (1 - x) ** n * (au + bu * x + cu * x**2 + du * x**3) * kern / u
    return f"{lhs} = \\int_0^1 {sp.latex(body)} \\mathrm{{d}} x > 0"
