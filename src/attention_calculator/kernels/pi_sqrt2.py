"""pi_sqrt2 kernel: claims ``q·π√2 ⋚ r`` (exact-mode-only type).

The last clean single-kernel constant on the lemniscate lattice
(docs/2026-09-16-w3-research-gamma-quarter.md §建议核形): with
S := π√2 = Γ(1/4)Γ(3/4) (reflection formula), the two ``(1-x⁴)^{c-1}``
kernels that still pair an S-class moment with the rational class are

    '>' : K(x) = (1-x⁴)^{-3/4}, basis x^{4m+2}(1-x)(a+bx⁴)
    '<' : K(x) = (1-x⁴)^{-1/4}, basis x^{4m+3}(1-x)(a+bx⁴)

The beta moments J'_k = ∫₀¹x^k(1-x⁴)^{-3/4}dx = ¼B((k+1)/4, 1/4) and
J''_k = ¼B((k+1)/4, 3/4) reduce to interleaved rational recurrences:

    J'_{4j+2}  = u_j·S/4     u_j = (3/4)_j/j!,        u_{j+1} = u_j(4j+3)/(4j+4)
    J'_{4j+3}  = v_j         v_j = (1)_j/(5/4)_j,     v_{j+1} = v_j(4j+4)/(4j+5)
    J''_{4j+3} = w_j         w_j = (1/3)(1)_j/(7/4)_j, w_{j+1} = w_j(4j+4)/(4j+7)
    J''_{4j+4} = z_j·S       z_j = (1/16)(5/4)_j/(2)_j, z_{j+1} = z_j(4j+5)/(4j+8)

so the '>' slot-i moment is ``{S: u_{m+i}/4, 1: -v_{m+i}}`` and the '<'
slot-i moment is ``{1: w_{m+i}, S: -z_{m+i}}`` — the same two-symbol
{C, 1} moment shape as the varpi/gauss lemniscate kernels, including the
fixed (1-x) basis factor and the a+bx⁴ nonneg rule (a >= 0 and a+b >= 0,
engine.poly_nonneg). Both directions prove the literal claim directly —
'>' emits ∫P·K = q·S - r, '<' emits ∫P·K = r - q·S; no inverse-symbol
normalization is needed since each kernel's own span already holds S.

power is a pure coefficient — the moments never see q — so q <= 0 needs
no validation: certified_cmp upstream rejects false claims, and a true
claim with non-positive q gives sign-indefinite solves at every m (a
uniformly non-positive P would certify the opposite inequality), ending
in honest NoSolution — the same convention varpi and gamma follow. The
emitted n is a direction flag (0 for '>', 1 for '<'), not an exponent —
same encoding as beta_family.

Search budget: single m-axis, LIMIT = 256 — the lemniscate family's
exact-mode '<' budget (beta_family.EXACT_LT_LIMIT). The four verified
examples need m <= 13; asymptotics u_m/v_m → 4/S and z_m/w_m → 1/S make
every true bound resolvable, with '<' wanting ~3-4x the '>' depth
(gap ~0.057 needs m=13, ~0.007 needs m ~ 116). A full 256-deep miss
costs ~10ms, so over-budget claims report honest NoSolution.
"""

from fractions import Fraction

import sympy as sp

from ..engine import search
from ..moment import Moment
from ..render import coef_tex, emit, lhs_tex, rat_tex

LIMIT = 256

SYM = "pi_sqrt2"

x = sp.symbols("x")


def gt_table(n: int) -> tuple[list[Fraction], list[Fraction]]:
    """u_j = (3/4)_j/j! and v_j = (1)_j/(5/4)_j for j = 0..n ('>' kernel)."""
    u, v = [Fraction(1)], [Fraction(1)]
    for j in range(n):
        u.append(u[-1] * Fraction(4 * j + 3, 4 * j + 4))
        v.append(v[-1] * Fraction(4 * j + 4, 4 * j + 5))
    return u, v


def lt_table(n: int) -> tuple[list[Fraction], list[Fraction]]:
    """w_j = (1/3)(1)_j/(7/4)_j and z_j = (1/16)(5/4)_j/(2)_j, j = 0..n ('<')."""
    w, z = [Fraction(1, 3)], [Fraction(1, 16)]
    for j in range(n):
        w.append(w[-1] * Fraction(4 * j + 4, 4 * j + 7))
        z.append(z[-1] * Fraction(4 * j + 5, 4 * j + 8))
    return w, z


def basis(comp: str, m: int, i: int) -> Moment:
    """Moment of ``x^{4m+res+4i}(1-x)K(x)`` for the i-th P slot (a, b).

    '>' (res=2, K=(1-x⁴)^{-3/4}): ``{S: u_{m+i}/4, 1: -v_{m+i}}``;
    '<' (res=3, K=(1-x⁴)^{-1/4}): ``{1: w_{m+i}, S: -z_{m+i}}``.
    """
    if comp == ">":
        u, v = gt_table(m + 1)
        return {SYM: u[m + i] / 4, "1": -v[m + i]}
    w, z = lt_table(m + 1)
    return {"1": w[m + i], SYM: -z[m + i]}


def target(comp: str, power: Fraction, bound: Fraction) -> Moment:
    """The claim vector: ``q·S - r`` for '>', ``r - q·S`` for '<'."""
    t = {SYM: power, "1": -bound} if comp == ">" else {"1": bound, SYM: -power}
    return {k: v for k, v in t.items() if v}


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("pi_sqrt2", q, comp, bound) -> /calculate-shaped dict.

    Exact-mode-only type; ``exact`` is accepted for signature
    compatibility (there is no site path to reproduce).
    """
    flag = 0 if comp == ">" else 1
    plans = ((m, flag, [basis(comp, m, i) for i in (0, 1)]) for m in range(LIMIT + 1))
    solved = search(plans, target(comp, power, bound), True)
    result = emit(solved.m, solved.n, solved.coeffs)
    result["type"] = kind
    return result


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Proof-equation LaTeX: ``q·π√2 - r`` / ``r - q·π√2`` = ∫x^e(1-x)P·K.

    No site counterpart exists; the integrand prints via sympy with the
    kernel as a denominator power ``(1-x⁴)^{3/4}`` / ``\\sqrt[4]{1-x⁴}``.
    """
    m = int(params["m"])
    au, bu, u = (int(params[k]) for k in ("au_val", "bu_val", "u_val"))
    e, c = (4 * m + 2, sp.Rational(3, 4)) if comp == ">" else (4 * m + 3, sp.Rational(1, 4))
    body = x**e * (1 - x) * (au + bu * x**4) / (u * (1 - x**4) ** c)
    lhs = lhs_tex(coef_tex(power) + "\\pi\\sqrt{2}", rat_tex(bound), comp)
    return f"{lhs} = \\int_0^1 {sp.latex(body)} \\mathrm{{d}} x > 0"
