"""Dixon/Γ(1/3)-lattice kernel family on [0, 1]: ``pi3``, ``pi3_u``, ``pi3_a``.

Exact-mode-only types for the three constants of the equianharmonic lattice
(docs/2026-09-16-w3-research-gamma-third.md):

    pi3:    π₃ = B(1/3,1/3) = √3·Γ(1/3)³/(2π) ≈ 5.2999 — the Dixon period
    pi3_u:  U  = A/π₃ = Γ(2/3)²/Γ(1/3) ≈ 0.6845
    pi3_a:  A  = B(1/3,2/3) = 2π/√3 ≈ 3.6276 — a π×algebraic constant like φ

The machine is the direct cubic analogue of the varpi/gauss lemniscate pair:
the kernels K = (1−x³)^{−k/3} (k = 1, 2) have Beta moments
T_k = (1/3)B((k+1)/3, 1−k/3) that sort by k mod 3 into three rational
sequences,

    k=1 (s=−1/3): T_{3m} = (A/3)·a_m,  T_{3m+1} = U·b_m,   T_{3m+2} = c_m/2
    k=2 (s=−2/3): T_{3m} = (π₃/3)·a'_m, T_{3m+1} = (A/3)·b'_m, T_{3m+2} = c'_m

    a_{m+1} = a_m(3m+1)/(3m+3)     a'_{m+1} = a'_m(3m+1)/(3m+2)
    b_{m+1} = b_m(3m+2)/(3m+4)     b'_{m+1} = b'_m(3m+2)/(3m+3)
    c_{m+1} = c_m(3m+3)/(3m+5)     c'_{m+1} = c'_m(3m+3)/(3m+4)

The basis x^{3m+res}(1−x)(a+bx³) pairs adjacent residue classes into a 2-dim
span, and the outer factor √3/(2π) = 1/A — the exact analogue of varpi's 1/π —
renames the moments' symbols (A→1, U→π₃⁻¹ on kernel 1; π₃→G₃, A→1 on kernel
2, where G₃ = π₃/A = U⁻¹). The six (kind, comp) pairs are a bijection onto
(kernel, res) cells:

    (kind, comp)  kernel  res  √3/2π   span            proves
    pi3     '<'     2      2     —     {1, pi3}         bound − q·π₃
    pi3     '>'     1      0     ✓     {1, pi3_inv}     q − bound·π₃⁻¹
    pi3_u   '>'     1      1     —     {U, 1}           q·U − bound
    pi3_u   '<'     2      0     ✓     {G3, 1}          bound·U⁻¹ − q
    pi3_a   '<'     1      2     —     {1, A}           bound − q·A
    pi3_a   '>'     2      1     —     {A, 1}           q·A − bound

(The '<'-via-reciprocal row mirrors gauss '<': q·U < r ⟺ r·U⁻¹ − q > 0.)

Parameters: ``m`` is the real exponent level, ``n`` is a direction flag
(0 for '>', 1 for '<') — the beta_family/pi_sqrt2 encoding — and
``(au+bu·x³)/u`` is P. Search walks m = 0..LIMIT; the edge ratio of the
one-sided basis closes on the constant like ~1/m, so claims within ~1e-3 of
the constant resolve inside budget — deeper claims report honest NoSolution
(W4 adaptive budget is future work).

Bounds looser than the m=0 interval edge fall to the varpi-style ``+b``
template (cu_val=1): ``t·x^{3m+res}(1−x)·outer·K + b`` with t fixed by the
constant coefficient and b the rational leftover; emitted iff t ≥ 0 and
b > 0 — the scan takes the shallowest m whose moment ratio clears the bound.

power is the coefficient q of q·C ⋚ r (varpi/gauss convention). q = 0 is
rejected: the claim degenerates to a rational comparison no integral proves.
"""

from fractions import Fraction

import sympy as sp

from ..engine import NoSolution, search
from ..moment import Moment, combine
from ..render import coef_tex, emit, lhs_tex, rat_tex, wire_or

LIMIT = 512

x = sp.symbols("x")

SYM_TEX = {"pi3": "\\pi_{3}", "pi3_u": "U", "pi3_a": "A"}

# (kind, comp) -> (kernel exponent k of (1-x^3)^{-k/3}, res, outer √3/(2π) flag)
CONFIG = {
    ("pi3", "<"): (2, 2, False),
    ("pi3", ">"): (1, 0, True),
    ("pi3_u", ">"): (1, 1, False),
    ("pi3_u", "<"): (2, 0, True),
    ("pi3_a", "<"): (1, 2, False),
    ("pi3_a", ">"): (2, 1, False),
}

# multiplying the integrand by √3/(2π) renames each raw symbol; the "1" symbol
# never appears in the res=0 pairings this factor is used with
OUTER_FACTOR = {
    1: {"A": "1", "U": "pi3_inv"},
    2: {"pi3": "G3", "A": "1"},
}


def third_table(kernel: int, n: int):
    """The kernel's three rational moment sequences, each of length n+1."""
    a, b, c = [Fraction(1)], [Fraction(1)], [Fraction(1)]
    for m in range(n):
        if kernel == 1:
            a.append(a[-1] * Fraction(3 * m + 1, 3 * m + 3))
            b.append(b[-1] * Fraction(3 * m + 2, 3 * m + 4))
            c.append(c[-1] * Fraction(3 * m + 3, 3 * m + 5))
        else:
            a.append(a[-1] * Fraction(3 * m + 1, 3 * m + 2))
            b.append(b[-1] * Fraction(3 * m + 2, 3 * m + 3))
            c.append(c[-1] * Fraction(3 * m + 3, 3 * m + 4))
    return a, b, c


def t_moment(kernel: int, k: int, tab) -> Moment:
    """T_k = ∫_0^1 x^k (1-x^3)^{-kernel/3} dx over the three residue classes."""
    a, b, c = tab
    m, rem = divmod(k, 3)
    if kernel == 1:
        if rem == 0:
            return {"A": a[m] / 3}
        if rem == 1:
            return {"U": b[m]}
        return {"1": c[m] / 2}
    if rem == 0:
        return {"pi3": a[m] / 3}
    if rem == 1:
        return {"A": b[m] / 3}
    return {"1": c[m]}


def dixon_basis(kernel: int, res: int, outer: bool, tab, m: int, i: int) -> Moment:
    """Moment of ``x^{3m+res+3i}(1-x) · outer · (1-x^3)^{-kernel/3}``.

    ``i`` indexes the (a, b) unknowns of P = a + b·x³ — slot i multiplies the
    basis monomial by x^{3i}. ``tab`` is third_table(kernel, m+3) or longer.
    """
    k = 3 * m + res + 3 * i
    base = combine(
        [Fraction(1), Fraction(-1)],
        [t_moment(kernel, k, tab), t_moment(kernel, k + 1, tab)],
    )
    if outer:
        base = {OUTER_FACTOR[kernel][sym]: v for sym, v in base.items()}
    return base


def dixon_target(kind: str, comp: str, power: Fraction, bound: Fraction) -> Moment:
    """The claimed vector: q·C ⋚ r as a span-{C-symbols, 1} moment."""
    if kind == "pi3":
        t = {"1": power, "pi3_inv": -bound} if comp == ">" else {"1": bound, "pi3": -power}
    elif kind == "pi3_u":
        t = {"U": power, "1": -bound} if comp == ">" else {"G3": bound, "1": -power}
    else:  # pi3_a
        t = {"A": power, "1": -bound} if comp == ">" else {"1": bound, "A": -power}
    return {k: v for k, v in t.items() if v}


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("pi3", q, comp, bound) -> /calculate-shaped dict.

    Exact-mode-only family; ``exact`` is accepted for signature compatibility
    (there is no site path to reproduce).
    """
    if power == 0:
        raise ValueError("系数不能为0")
    kernel, res, outer = CONFIG[(kind, comp)]
    flag = 0 if comp == ">" else 1
    tab = third_table(kernel, LIMIT + 3)
    plans = (
        (m, flag, [dixon_basis(kernel, res, outer, tab, m, i) for i in (0, 1)])
        for m in range(LIMIT + 1)
    )
    try:
        solved = search(plans, dixon_target(kind, comp, power, bound), True)
    except NoSolution:
        return bound_proof(kind, comp, power, bound)
    result = emit(solved.m, solved.n, solved.coeffs)
    result["type"] = kind
    return result


def bound_proof(kind: str, comp: str, power: Fraction, bound: Fraction) -> dict:
    """The ``+b`` loose-bound fallback: ``∫ t·x^{3m+res}(1-x)·outer·K dx + b``.

    The slot-0 basis moment is {sym: α, "1": β}; t is fixed by the target's
    constant coefficient and b is the rational leftover. Emitted at the
    shallowest m with t ≥ 0 and b > 0 — the moment ratio closes on the
    constant from the loose side, so every bound beyond the 2×2 search's
    interval resolves here (the varpi/gauss cu_val=1 template, with the real
    m and res reported in the m/n slots).
    """
    kernel, res, outer = CONFIG[(kind, comp)]
    flag = 0 if comp == ">" else 1
    target = dixon_target(kind, comp, power, bound)
    sym = next((s for s in target if s != "1"), None)
    if sym is not None:
        tab = third_table(kernel, LIMIT + 3)
        for m in range(LIMIT + 1):
            mom = dixon_basis(kernel, res, outer, tab, m, 0)
            t = target[sym] / mom[sym]
            b = target.get("1", Fraction(0)) - t * mom.get("1", Fraction(0))
            if t >= 0 and b > 0:
                return emit_bound(kind, m, flag, t, b)
    elif target["1"] > 0:
        # bound == 0: the claim is its rational part — t = 0 kills the moment
        return emit_bound(kind, 0, flag, Fraction(0), target["1"])
    raise NoSolution


def emit_bound(kind: str, m: int, flag: int, t: Fraction, b: Fraction) -> dict:
    """The cu_val=1 parameter shape: au/u = t, b_val = b, other slots zero."""
    return {
        "type": kind,
        "parameters": {
            "m": m,
            "n": flag,
            "a_val": "0",
            "b_val": str(b),
            "c_val": "0",
            "au_val": str(t.numerator),
            "bu_val": "0",
            "cu_val": "1",
            "u_val": str(t.denominator),
            "unified_form": {},
        },
        "solution": "a = 0, b = 0",
    }


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX for a solved case (exact types have no site route).

    The reciprocal-direction claims print the reciprocal LHS like the site's
    varpi '>'/gauss '<': ``q − r·π₃^{-1}`` and ``r·U^{-1} − q``.
    """
    m = int(params["m"])
    au, bu, u = (int(params[k]) for k in ("au_val", "bu_val", "u_val"))
    kernel, res, outer = CONFIG[(kind, comp)]
    btex = rat_tex(bound)

    if kind == "pi3" and comp == ">":
        inv = "\\pi_{3}^{-1}" if wire_or(bound) == 1 else f"{btex}\\pi_{{3}}^{{-1}}"
        lhs = f"{rat_tex(power)} - {inv}"
    elif kind == "pi3_u" and comp == "<":
        inv = "U^{-1}" if wire_or(bound) == 1 else f"{btex}U^{{-1}}"
        lhs = f"{inv} - {rat_tex(power)}"
    else:
        lhs = lhs_tex(coef_tex(power) + SYM_TEX[kind], btex, comp)

    kern = (1 - x**3) ** sp.Rational(-kernel, 3)
    fac = sp.sqrt(3) / (2 * sp.pi) if outer else sp.Integer(1)
    if int(params["cu_val"]):
        inner = sp.latex(sp.Rational(au, u) * x ** (3 * m + res) * (1 - x) * fac * kern)
        return f"{lhs} = \\int_0^1 {inner} \\mathrm{{d}} x+{rat_tex(params['b_val'])} > 0"
    body = x ** (3 * m + res) * (1 - x) * (au + bu * x**3) * fac * kern / u
    return f"{lhs} = \\int_0^1 {sp.latex(body)} \\mathrm{{d}} x > 0"
