"""pi-multiple / degree trig kernels on [0, pi/2]: sin_pi_q, cos_pi_q,
sin_q_degree, cos_q_degree (article types 15/17).

The integrand is sin^m x (1 - sin x)^n (a + b sin x) sin(alpha x) with

    alpha = 1 - 2q   for sin(pi q) and sin(q deg),
    alpha = 2 q      for cos(pi q) and cos(q deg),

where q = degrees/180 for the *_q_degree types. The moments
T_j = int_0^{pi/2} sin^j x sin(alpha x) dx are evaluated by expanding sin^j
into multiple angles:

    j odd:  sin^j x = 2^{1-j} sum_l (-1)^{(j-1)/2 - l} C(j,l) sin((j-2l) x)
    j even: sin^j x = 2^{-j} C(j,j/2) + 2^{1-j} sum_l (-1)^{j/2 - l} C(j,l)
                    cos((j-2l) x)

then sin(kx) sin(ax) = [cos((a-k)x) - cos((a+k)x)]/2 and
cos(kx) sin(ax) = [sin((a+k)x) + sin((a-k)x)]/2 integrate termwise into
sin/cos((a+-k)pi/2). Because k always has the parity of j, both reductions
collapse onto +/-cos(pi alpha/2): every moment lands in
span{cos(pi alpha/2), 1}, which is span{sin(pi q), 1} for alpha = 1-2q and
span{cos(pi q), 1} for alpha = 2q.

The site admits a non-integer fraction q in (0, 1/2) for *_pi_q and any
rational in (0, 90) degrees for *_q_degree, answering "type": "sin_pi_q" /
"cos_pi_q" either way. The only rational values reachable are sin(pi/6) =
cos(pi/3) = 1/2 (Niven); hitting them exactly yields "二者相等".
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import sympy as sp
from mpmath import mp

from ..engine import EqualClaim, WrongDirection, mn_order, search
from ..moment import Moment, combine
from ..render import cdot_tex, emit, lhs_tex, rat_tex

LIMIT = 10

x = sp.symbols("x")


def angle_moment(j: int, alpha: Fraction) -> Moment:
    """T_j = int_0^{pi/2} sin^j x sin(alpha x) dx over {C, 1}; C = cos(pi*alpha/2)."""
    out: Moment = {}

    def addto(key: str, v: Fraction) -> None:
        out[key] = out.get(key, Fraction(0)) + v

    if j % 2:
        for i in range((j + 1) // 2):
            k = j - 2 * i
            c = Fraction((-1) ** ((j - 1) // 2 - i) * comb(j, i), 2**j)
            for s, sign in ((-1, 1), (1, -1)):
                # int cos((alpha+s*k)x) -> sin(theta + s*k*pi/2)/(alpha+s*k);
                # odd k: sin(theta + s*k*pi/2) = s*(-1)^{(k-1)/2} cos(theta)
                addto("C", sign * c * s * (-1) ** ((k - 1) // 2) / (alpha + s * k))
    else:
        c0 = Fraction(comb(j, j // 2), 2**j)
        addto("1", c0 / alpha)
        addto("C", -c0 / alpha)
        for i in range(j // 2):
            k = j - 2 * i
            c = Fraction((-1) ** (j // 2 - i) * comb(j, i), 2**j)
            for s in (1, -1):
                # int sin((alpha+s*k)x) -> (1 - cos(theta + s*k*pi/2))/(alpha+s*k);
                # even k: cos(theta + s*k*pi/2) = (-1)^{k/2} cos(theta)
                beta = alpha + s * k
                addto("1", c / beta)
                addto("C", -c * (-1) ** (k // 2) / beta)
    return {k: v for k, v in out.items() if v}


def basis_moment(t: list[Moment], m: int, n: int, j: int) -> Moment:
    """int sin^{m+j} x (1-sin x)^n sin(alpha x) dx = sum_i (-1)^i C(n,i) T_{m+i+j}."""
    return combine(
        [Fraction((-1) ** i * comb(n, i)) for i in range(n + 1)],
        [t[m + i + j] for i in range(n + 1)],
    )


# The site does not integrate; it substitutes into Mathematica-precomputed
# closed forms per (m, n). Its stored j=0 formula at (m, n) = (1, 8) is off by
# delta(alpha) * (C - 1), recovered by rational interpolation over 33 probed
# alpha values (docs/kernel-spec.md). The bias denominator below is exactly the
# natural denominator of the true moment, so the stored formula is the true
# one plus this term. The corrupted moment feeds both the (a, b) solve and the
# sign check: bounds where the biased solve goes nonnegative emit a false
# identity at (1, 8); where it goes nonpositive the plan is skipped.
BIAS18_NUM = (
    8393,
    30548,
    -2295970,
    -8426378,
    249909164,
    928663806,
    -13829413586,
    -52408975814,
    410416416763,
    1608379190294,
    -6230112405308,
    -26062947629208,
    39020672969280,
    192843476410752,
    -20544307623936,
    -376544377036800,
    33710564966400,
)


def bias_18(alpha: Fraction) -> Fraction:
    """delta = (alpha-1)(alpha-2) Q16(alpha) / (256 alpha prod (alpha^2-k^2))."""
    p = Fraction(BIAS18_NUM[0])
    for c in BIAS18_NUM[1:]:
        p = p * alpha + c
    num = (alpha - 1) * (alpha - 2) * p
    den = Fraction(256) * alpha
    for k in range(1, 10):
        den *= alpha * alpha - k * k
    return num / den


def site_basis(t: list[Moment], m: int, n: int, bias: Fraction) -> list[Moment]:
    """Basis moments as the site's stored formulas give them: biased at (1, 8)."""
    basis = [basis_moment(t, m, n, j) for j in range(2)]
    if (m, n) == (1, 8):
        basis[0] = dict(basis[0])
        basis[0]["C"] = basis[0].get("C", Fraction(0)) + bias
        basis[0]["1"] = basis[0].get("1", Fraction(0)) - bias
    return basis


def check_input(q: Fraction, bound: Fraction, kind: str) -> tuple[Fraction, str]:
    """Validate like the site; return (alpha, normalized type)."""
    if q < 0:
        raise ValueError("左侧系数格式无效")
    if bound < 0:
        raise ValueError("右侧有理数格式无效")
    is_sin = kind.startswith("sin")
    if kind.endswith("_degree"):
        if q <= 0 or q >= 90:
            raise ValueError("请在输入一个在(0,90)内的数")
        q_eff = q / 180
    else:
        if not 0 < q < Fraction(1, 2) or q.denominator == 1:
            raise ValueError("请在输入一个在(0,1/2)内的分数，本情况不支持整数")
        q_eff = q
    # Niven: the only rational values are sin(pi/6) = cos(pi/3) = 1/2.
    if q_eff == (Fraction(1, 6) if is_sin else Fraction(1, 3)) and bound == Fraction(1, 2):
        raise EqualClaim("二者相等")
    alpha = 1 - 2 * q_eff if is_sin else 2 * q_eff
    return alpha, "sin_pi_q" if is_sin else "cos_pi_q"


def solve(kind: str, q: Fraction, comp: str, bound: Fraction):
    """Find (m, n, a, b) whose integrand proves the requested inequality."""
    alpha, _ = check_input(q, bound, kind)
    s = Fraction(1 if comp == ">" else -1)
    target = {"C": s, "1": -s * bound}
    # The site sanity-checks the claimed direction against a numerical value
    # of the constant (C = cos(pi alpha/2) equals sin(pi q) resp. cos(pi q))
    # before searching: a wrong-direction request never reaches the scan, so
    # even a corrupted-formula hit cannot emit when bound ~ C fails. Exact
    # equality is impossible here -- Niven cases were caught by check_input.
    with mp.workdps(50):
        const = mp.cos(mp.pi * mp.mpf(alpha.numerator) / mp.mpf(alpha.denominator) / 2)
        gap = mp.mpf(bound.numerator) / mp.mpf(bound.denominator) - const
    if (comp == "<") == (gap <= 0):
        raise WrongDirection
    t = [angle_moment(j, alpha) for j in range(2 * LIMIT + 2)]
    bias = bias_18(alpha)
    plans = ((m, n, site_basis(t, m, n, bias)) for m, n in mn_order(LIMIT))
    return search(plans, target, True, defer=True), alpha


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict:
    """Run the proof search for one trig_pi-family request."""
    solved, alpha = solve(kind, power, comp, bound)
    result = emit(solved.m, solved.n, solved.coeffs, c_val=str(alpha))
    result["type"] = "sin_pi_q" if kind.startswith("sin") else "cos_pi_q"
    return result


def render_equation(
    params: dict, kind: str, power: Fraction | str, comp: str, bound: Fraction | str
) -> str:
    """Reproduce the site's /get_integral_image LaTeX for this family.

    ``power``/``bound`` may be the raw request strings (echoed unreduced like
    the site does) or already-parsed Fractions.
    """
    m, n = params["m"], params["n"]
    au, bu, u = (Fraction(params[k]) for k in ("au_val", "bu_val", "u_val"))
    alpha = Fraction(params["c_val"])
    s = sp.sin(x)

    arg = rat_tex(power) + ("^\\circ" if kind.endswith("_degree") else r"\pi")
    name = "sin" if kind.startswith("sin") else "cos"
    const = rf"\{name}\left({arg}\right)"
    lhs = lhs_tex(const, rat_tex(bound), comp)

    # The site latexes the whole product; its printer puts " \cdot " before a
    # digit-leading parenthesized factor following a number or a power
    # (\d or }), a space elsewhere -- e.g. "(1-s)^8 \cdot (a+b s)" but
    # "... ^{3} \left(-1465926 ..." (a negative group gets no cdot).
    integrand = (au + bu * s) * (1 - s) ** n * sp.sin(alpha * x) * s**m / u
    body = cdot_tex(sp.latex(integrand))
    return lhs + rf" = \int_0^{{\pi/2}} {body} \mathrm{{d}} x > 0"
