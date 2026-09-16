"""mode=exact composite-inequality decomposition: provable-by-construction bounds.

Site counterpart (decompose.py) allocates record-chain bounds chosen for byte
parity with zhuyidao.net — no provability guarantee. This module allocates
rational sub-bounds so that (a) the conjunction certifies the claim and
(b) every sub-claim is provable within the (m,n) search budget, verified by
actually running ``solve.prove(..., exact=True)`` on each.

Model — every sub-claim certifies ``U ⋚ β`` where ``U = constant_mpf(kind,
power)`` is the unit the kernel's certified_cmp evaluates:

- single-atom term ``c·a``: coefficient kinds (pi, e, gamma, golden, catalan,
  zeta3, varpi, gauss, e_pi, plus the exact-mode coefficient kinds — odd
  zetas, even Dirichlet betas, pi_sqrt2, the pi3/Dixon trio, the
  gamma14/34/12 composite kinds)
  fold |c| into the power slot — U = |c|·C and the term is sign(c)·U;
  argument kinds (pi_n, e_q, ln_q, sin_q, and the exact-mode function
  kinds) keep the coefficient outside — U = a^sign(arg), term = c·U^sign(arg).
- product term ``c·Π a_j``: U_j = base value constant_mpf(kind_j, |arg_j|),
  term = c·Π U_j^sign(arg_j); negative exponents are reciprocal factors and
  the sub-proof runs on the base (e_q arg -2 proves e² ⋚ β, contributing 1/β).
  Only the exponent kinds pi_n/e_q read arg<0 as a reciprocal — signed-
  argument kinds (arsinh_q, si_q, cin_q, the gauss-erf family, li2_q, and
  the parity-folding sinh/cosh/tanh/arctan/arccot) pass a negative q to the
  power slot signed, since their kernels prove the literal signed claim.

Allocation (docs/2026-09-16-decompose-math.md): each kernel's reachable set is
a one-sided interval — a bound proves iff its margin |U − β| clears a
per-unit frontier floor ε(kind, power, comp). We measure ε by log-bisection
with the prove oracle, then split the certified slack S ∝ each term's
relative-frontier need ``|v_i|·Σ_j ε_j/U_j`` so every factor's realized
margin is the same multiple K = S/Σneed of its floor (equal safety ratio).
"""

from contextlib import suppress
from dataclasses import dataclass
from fractions import Fraction

import mpmath as mp
import sympy as sp

from . import render, solve
from .decompose import (
    BAD_CHARS,
    CATALAN,
    COEF_KINDS,
    ERR_ATOM,
    ERR_NUMERIC,
    ERR_POW_SPECIAL,
    ERR_RHS_FORM,
    FUNC_KINDS,
    GAMMA,
    GOLDEN,
    LOC,
    PI,
    E,
    atom_label,
    flip,
    frac,
    join_signed,
    merge_atoms,
    parse_atom,
)
from .engine import EqualClaim, NoSolution, WrongDirection
from .integrand import constant_mpf

DPS = 300  # precision for unit values, slack certification, bound snapping
MAX_ROUNDS = 5  # margin-escalation rounds after a NoSolution
UNPROVABLE = Fraction(10) ** 30  # sentinel floor: no bound within reach proves
FLOOR_CACHE: dict = {}  # (kind, power, comp) -> measured frontier floor


def mpf(fr: Fraction):
    """Fraction -> mpf at the ambient precision."""
    return mp.mpf(fr.numerator) / fr.denominator


def unit_value(kind: str, power: Fraction):
    """mpf value of the unit a sub-claim bounds (constant_mpf at ambient dps)."""
    return +constant_mpf(kind, power)


def frac60(x) -> Fraction:
    """Snap an mpf to a Fraction at ~60 significant digits (heuristic use only)."""
    return Fraction(mp.nstr(x, 60))


# ------------------------------------------------------------------ model


@dataclass
class Sub:
    """One atomic sub-claim: ``prove(kind, power, comp, bound)`` certifies U ⋚ bound."""

    kind: str
    power: Fraction  # value handed to solve.prove's power slot
    sigma: int  # +1, or -1 for a reciprocal factor (atom arg < 0)
    arg: Fraction  # original atom arg (sign kept) — display only
    label: str  # atom label (positive-power form)
    coef_disp: Fraction  # displayed coefficient (site's coefficient field)
    comp: str = ">"  # proved direction on U
    bound: Fraction = Fraction(0)  # rational bound on U
    floor: Fraction = Fraction(0)  # measured provability floor on U for comp
    u: object = None  # mpf value of U (filled at DPS)
    proof: dict | None = None
    error: str | None = None

    def t(self):
        """mpf factor value contributed to the term: U**sigma."""
        return self.u if self.sigma > 0 else 1 / self.u

    def phi(self) -> Fraction:
        """Rational factor bound contributed to the term bound: bound**sigma."""
        return self.bound if self.sigma > 0 else 1 / self.bound


@dataclass
class TermAlloc:
    """One lhs term: ``term = outer · Π sub.t``; certified bound ``outer · Π sub.phi``."""

    outer: Fraction
    subs: list
    coef: Fraction  # original signed coefficient (display)
    atoms: list  # parsed atoms (display)
    bound: Fraction = Fraction(0)

    def value(self):
        """mpf term value at ambient precision."""
        v = mpf(self.outer)
        for s in self.subs:
            v *= s.t()
        return v


def term_specs(coef: Fraction, atoms: list) -> tuple[Fraction, list]:
    """(outer, [Sub]) for one parsed term — folds coef per the unit model."""
    if len(atoms) == 1 and atoms[0][0] in EXACT_COEF_KINDS:
        kind, arg = atoms[0]
        sub = Sub(
            kind=kind,
            power=abs(coef),
            sigma=1,
            arg=arg,
            label=atom_label_exact(kind, arg),
            coef_disp=abs(coef),
        )
        return Fraction(-1 if coef < 0 else 1), [sub]
    subs = []
    for kind, arg in atoms:
        # only the exponent kinds read arg<0 as a reciprocal factor; signed
        # function kinds pass q to the power slot with its sign kept
        reciprocal = kind in EXPONENT_KINDS and arg < 0
        subs.append(
            Sub(
                kind=kind,
                power=abs(arg) if kind in EXPONENT_KINDS else arg,
                sigma=-1 if reciprocal else 1,
                arg=arg,
                label=atom_label_exact(kind, arg if kind in SIGNED_ARG_KINDS else abs(arg)),
                coef_disp=Fraction(1),
            )
        )
    return coef, subs


# ------------------------------------------------------------------ exact atoms
# decompose.py's BASE/LOC/COEF_KINDS/FUNC_KINDS are frozen site parity; the
# exact path admits every kind solve.prove(exact=True) can certify, so the
# supersets and the extended parser live here. Coefficient kinds carry the
# constant's multiplier in the power slot (q·C, like pi/e); function kinds
# carry the argument q (like ln_q).

ZETA5 = sp.Symbol("zeta5")
ZETA7 = sp.Symbol("zeta7")
ZETA9 = sp.Symbol("zeta9")
ZETA11 = sp.Symbol("zeta11")
BETA4 = sp.Symbol("beta4")
BETA6 = sp.Symbol("beta6")
BETA8 = sp.Symbol("beta8")
BETA10 = sp.Symbol("beta10")
PI_SQRT2 = sp.Symbol("pi_sqrt2")
PI3 = sp.Symbol("pi3")
PI3_U = sp.Symbol("pi3_u")
PI3_A = sp.Symbol("pi3_a")
GAMMA14 = sp.Symbol("gamma14")
GAMMA34 = sp.Symbol("gamma34")
GAMMA12 = sp.Symbol("gamma12")

EXACT_BASE = {
    ZETA5: ("zeta5", Fraction(1)),
    ZETA7: ("zeta7", Fraction(1)),
    ZETA9: ("zeta9", Fraction(1)),
    ZETA11: ("zeta11", Fraction(1)),
    BETA4: ("beta4", Fraction(1)),
    BETA6: ("beta6", Fraction(1)),
    BETA8: ("beta8", Fraction(1)),
    BETA10: ("beta10", Fraction(1)),
    PI_SQRT2: ("pi_sqrt2", Fraction(1)),
    PI3: ("pi3", Fraction(1)),
    PI3_U: ("pi3_u", Fraction(1)),
    PI3_A: ("pi3_a", Fraction(1)),
    GAMMA14: ("gamma14", Fraction(1)),
    GAMMA34: ("gamma34", Fraction(1)),
    GAMMA12: ("gamma12", Fraction(1)),
}

# Spelling -> sympy Function. Every function name must be an undefined
# Function in locals — sympy auto-evaluates builtins like asin(1/2) -> pi/6,
# which would silently re-parse as a pi coefficient atom. Multiple spellings
# share one Function object so func.__name__ stays the canonical one.
ARCSIN = sp.Function("arcsin")
ARSINH = sp.Function("arsinh")
GAUSSINT = sp.Function("gaussint")
DAWSON = sp.Function("dawson")
ERFIINT = sp.Function("erfiint")
LI2 = sp.Function("li2")
PSI1 = sp.Function("psi1")
SI = sp.Function("Si")
CIN = sp.Function("Cin")
COT = sp.Function("cot")
COSH = sp.Function("cosh")
COTH = sp.Function("coth")
ARCCOT = sp.Function("arccot")
ARTANH = sp.Function("artanh")
ARCOTH = sp.Function("arcoth")
SINPI = sp.Function("sinpi")
COSPI = sp.Function("cospi")
SIND = sp.Function("sind")
COSD = sp.Function("cosd")
GAMMA_F = sp.Function("Gamma")

EXACT_LOC = LOC | {
    "zeta5": ZETA5,
    "zeta7": ZETA7,
    "zeta9": ZETA9,
    "zeta11": ZETA11,
    "beta4": BETA4,
    "beta6": BETA6,
    "beta8": BETA8,
    "beta10": BETA10,
    "pi_sqrt2": PI_SQRT2,
    "pi3": PI3,
    "pi3_u": PI3_U,
    "pi3_a": PI3_A,
    "gamma14": GAMMA14,
    "gamma34": GAMMA34,
    "gamma12": GAMMA12,
    "asin": ARCSIN,
    "arcsin": ARCSIN,
    "asinh": ARSINH,
    "arsinh": ARSINH,
    "gaussint": GAUSSINT,
    "dawson": DAWSON,
    "erfiint": ERFIINT,
    "li2": LI2,
    "psi1": PSI1,
    "trigamma": PSI1,
    "si": SI,
    "Si": SI,
    "cin": CIN,
    "Cin": CIN,
    "cot": COT,
    "cosh": COSH,
    "coth": COTH,
    "arccot": ARCCOT,
    "acot": ARCCOT,
    "artanh": ARTANH,
    "atanh": ARTANH,
    "arcoth": ARCOTH,
    "acoth": ARCOTH,
    "sinpi": SINPI,
    "cospi": COSPI,
    "sind": SIND,
    "cosd": COSD,
    "Gamma": GAMMA_F,
    # "zeta"/"log" stay unmapped on purpose: the builtins keep zeta(2k)
    # auto-reducing to pi powers and log(q) handled by parse_atom's table.
}

# func.__name__ -> kind for the function atoms registered in EXACT_LOC.
EXACT_FUNC_TABLE = {
    "arcsin": "arcsin_q",
    "arsinh": "arsinh_q",
    "gaussint": "gaussint_q",
    "dawson": "dawson_q",
    "erfiint": "erfiint_q",
    "li2": "li2_q",
    "psi1": "psi1_q",
    "Si": "si_q",
    "Cin": "cin_q",
    "cot": "cot_q",
    "cosh": "cosh_q",
    "coth": "coth_q",
    "arccot": "arccot_q",
    "artanh": "artanh_q",
    "arcoth": "arcoth_q",
    "sinpi": "sin_pi_q",
    "cospi": "cos_pi_q",
    "sind": "sin_q_degree",
    "cosd": "cos_q_degree",
}

# ln(q)^e / log(q)^e atoms; e=1 collapses to plain ln before parsing.
LN_POW_KINDS = {2: "ln_q_square", 3: "ln_q_cube", 4: "ln_q_quad"}
# zeta(k) builtin calls map onto odd-zeta coefficient kinds (arg -> bare atom).
ZETA_INT_KINDS = {3: "zeta3", 5: "zeta5", 7: "zeta7", 9: "zeta9", 11: "zeta11"}
# dirichlet_beta(k) builtin calls (the builtin never auto-evaluates here) map
# onto Dirichlet-beta coefficient kinds; beta(2) is Catalan.
BETA_INT_KINDS = {2: "catalan", 4: "beta4", 6: "beta6", 8: "beta8", 10: "beta10"}
# Gamma(1/4 | 1/2 | 3/4) spellings for the composite-gamma coefficient kinds.
GAMMA_FRAC_KINDS = {
    Fraction(1, 4): "gamma14",
    Fraction(1, 2): "gamma12",
    Fraction(3, 4): "gamma34",
}

EXACT_COEF_KINDS = COEF_KINDS | {
    "zeta5",
    "zeta7",
    "zeta9",
    "zeta11",
    "beta4",
    "beta6",
    "beta8",
    "beta10",
    "pi_sqrt2",
    "pi3",
    "pi3_u",
    "pi3_a",
    "gamma14",
    "gamma34",
    "gamma12",
}
EXACT_FUNC_KINDS = FUNC_KINDS | set(EXACT_FUNC_TABLE.values()) | set(LN_POW_KINDS.values())
EXPONENT_KINDS = {"pi_n", "e_q"}  # atom arg is an exponent: arg<0 = reciprocal
# Kernels that accept a signed q (parity fold or a signed moment span): the
# atom arg passes into the power slot with its sign — never a reciprocal.
SIGNED_ARG_KINDS = {
    "arsinh_q",
    "si_q",
    "cin_q",
    "gaussint_q",
    "dawson_q",
    "erfiint_q",
    "li2_q",
    "sinh_q",
    "cosh_q",
    "tanh_q",
    "arctan_q",
    "arccot_q",
}

PI2 = sp.pi / 2

# sympy auto-evaluations that produce builtin constants (zeta(4) -> pi^4/90
# with the real sympy pi, dirichlet_beta(2) -> Catalan, bare "E"/"EulerGamma"/
# "GoldenRatio" spellings) get folded back onto the atom symbols.
SYMPY_CONST = {
    sp.pi: PI,
    sp.E: E,
    sp.EulerGamma: GAMMA,
    sp.Catalan: CATALAN,
    sp.GoldenRatio: GOLDEN,
}


def _lt(q: Fraction, x) -> bool:
    """Exact Fraction < sympy-transcendental comparison (for the pi bounds)."""
    return bool(sp.Rational(q.numerator, q.denominator) < x)


# Function-kind argument domains, mirroring each kernel's own check_input.
# The allocator needs every atom's unit value before any prove call, so
# out-of-domain args must be rejected here — otherwise constant_mpf leaks
# complex/NaN values (log of a negative, asin(2), cin(0)) or crashes.
ATOM_DOMAIN = {
    "ln_q": lambda q: q > 1,
    "ln_q_square": lambda q: q > 1,
    "ln_q_cube": lambda q: q > 1,
    "ln_q_quad": lambda q: q > 1,
    "arcsin_q": lambda q: 0 < q < 1,
    "artanh_q": lambda q: 0 < q < 1,
    "arcoth_q": lambda q: q > 1,
    "psi1_q": lambda q: q > 0,
    "coth_q": lambda q: q > 0,
    "li2_q": lambda q: q != 0 and q < 1,
    "sin_q": lambda q: q > 0 and _lt(q, sp.pi),
    "cos_q": lambda q: q > 0 and _lt(q, PI2),
    "tan_q": lambda q: q > 0 and _lt(q, PI2),
    "cot_q": lambda q: q > 0 and _lt(q, PI2),
    "sin_q_degree": lambda q: 0 < q < 90,
    "cos_q_degree": lambda q: 0 < q < 90,
    "sin_pi_q": lambda q: 0 < q < Fraction(1, 2),
    "cos_pi_q": lambda q: 0 < q < Fraction(1, 2),
    "arctan_q": lambda q: q != 0,
    "arccot_q": lambda q: q != 0,
    "sinh_q": lambda q: q != 0,
    "cosh_q": lambda q: q != 0,
    "tanh_q": lambda q: q != 0,
    "arsinh_q": lambda q: q != 0,
    "gaussint_q": lambda q: q != 0,
    "dawson_q": lambda q: q != 0,
    "erfiint_q": lambda q: q != 0,
    "si_q": lambda q: q != 0,
    "cin_q": lambda q: q != 0,
}


def check_domain(kind: str, arg: Fraction) -> None:
    """Reject a function atom whose argument leaves its kernel's domain."""
    ok = ATOM_DOMAIN.get(kind)
    if ok is not None and not ok(arg):
        raise ValueError(f"{kind} 不接受参数 {arg}：超出核定义域")


EXACT_CONST_LABEL = {
    "zeta5": "\\zeta(5)",
    "zeta7": "\\zeta(7)",
    "zeta9": "\\zeta(9)",
    "zeta11": "\\zeta(11)",
    "beta4": "\\beta(4)",
    "beta6": "\\beta(6)",
    "beta8": "\\beta(8)",
    "beta10": "\\beta(10)",
    "pi_sqrt2": "\\pi\\sqrt{2}",
    "pi3": "\\pi_{3}",
    "pi3_u": "U",
    "pi3_a": "A",
    "gamma14": "\\Gamma\\left(\\dfrac{1}{4}\\right)",
    "gamma34": "\\Gamma\\left(\\dfrac{3}{4}\\right)",
    "gamma12": "\\Gamma\\left(\\dfrac{1}{2}\\right)",
}
# Function-kind label heads matching each kernel's own rendered const name.
EXACT_FUNC_LABEL = {
    "ln_q_square": "\\ln^{2}",
    "ln_q_cube": "\\ln^{3}",
    "ln_q_quad": "\\ln^{4}",
    "arcsin_q": "\\arcsin",
    "arsinh_q": "\\operatorname{arsinh}",
    "gaussint_q": "\\mathrm{GaussInt}",
    "dawson_q": "\\mathrm{Dawson}",
    "erfiint_q": "\\mathrm{ErfiInt}",
    "li2_q": "\\mathrm{Li}_2",
    "psi1_q": "\\psi_1",
    "si_q": "\\mathrm{Si}",
    "cin_q": "\\mathrm{Cin}",
    "cot_q": "\\cot",
    "cosh_q": "\\cosh",
    "coth_q": "\\coth",
    "arccot_q": "\\mathrm{arccot}",
    "artanh_q": "\\operatorname{artanh}",
    "arcoth_q": "\\operatorname{arcoth}",
    "sin_pi_q": "\\sin",
    "cos_pi_q": "\\cos",
    "sin_q_degree": "\\sin",
    "cos_q_degree": "\\cos",
}


def atom_label_exact(kind: str, arg: Fraction) -> str:
    """atom_label over the extended tables; site kinds delegate unchanged."""
    if kind in EXACT_CONST_LABEL:
        return EXACT_CONST_LABEL[kind]
    if kind in EXACT_FUNC_LABEL:
        suffix = (
            "^{\\circ}" if kind.endswith("_degree") else "\\pi" if kind.endswith("_pi_q") else ""
        )
        return f"{EXACT_FUNC_LABEL[kind]}\\left({arg}{suffix}\\right)"
    return atom_label(kind, arg)


def factor_label_exact(kind: str, arg: Fraction) -> str:
    """factor_label over the extended tables; only exponent kinds print 1/base."""
    if arg < 0 and kind in EXPONENT_KINDS:
        return f"\\dfrac{{1}}{{{atom_label_exact(kind, -arg)}}}"
    return atom_label_exact(kind, arg)


def term_piece_exact(coef: Fraction, atoms: list) -> str:
    """term_piece over the extended tables (same layout as the site's)."""
    if len(atoms) == 1:
        kind, arg = atoms[0]
        piece = atom_label_exact(kind, arg)
        if coef == 1:
            return piece
        if coef == -1:
            return "-" + piece
        return frac(coef) + piece
    if coef == 1:
        head = ""
    elif coef == -1:
        head = "-"
    else:
        head = frac(coef) + "\\cdot "
    return head + "\\cdot ".join(factor_label_exact(k, a) for k, a in atoms)


# ------------------------------------------------------------------ parsing


def parse_atom_exact(expr):
    """parse_atom over the extended tables; site shapes delegate unchanged.

    Order matters: exact constant symbols first, then Pow (exact constants to
    a power -> ERR_POW_SPECIAL like the site; ln(q)^e -> the ln-power kinds;
    everything else -> the site rules), then Function (exact table, zeta(k)
    and Gamma(1/n) builtins onto coefficient kinds, else the site rules).
    """
    if expr in EXACT_BASE:
        return EXACT_BASE[expr]
    if expr.is_Pow:
        b, e = expr.base, expr.exp
        if b in EXACT_BASE:
            raise ValueError(ERR_POW_SPECIAL)
        if e.is_Integer and b.is_Function and b.func.__name__ in ("ln", "log"):
            a = b.args[0]
            kind = LN_POW_KINDS.get(int(e)) if a.is_Rational else None
            if kind is None:
                raise ValueError(ERR_ATOM)
            return kind, Fraction(a.p, a.q)
        return parse_atom(expr)
    if expr.is_Function:
        name = expr.func.__name__
        a = expr.args[0] if expr.args else None
        if a is not None and a.is_Rational:
            arg = Fraction(a.p, a.q)
            if name in EXACT_FUNC_TABLE:
                return EXACT_FUNC_TABLE[name], arg
            if name == "zeta" and a.is_Integer and int(a) in ZETA_INT_KINDS:
                return ZETA_INT_KINDS[int(a)], Fraction(1)
            if name == "dirichlet_beta" and a.is_Integer and int(a) in BETA_INT_KINDS:
                return BETA_INT_KINDS[int(a)], Fraction(1)
            if name == "Gamma" and arg in GAMMA_FRAC_KINDS:
                return GAMMA_FRAC_KINDS[arg], Fraction(1)
        return parse_atom(expr)
    return parse_atom(expr)


def single_atom_exact(expr):
    """decompose.single_atom over the extended tables: lone atom -> (kind, arg)."""
    try:
        coef = Fraction(1)
        factors = []
        for f in sp.Mul.make_args(expr):
            if f.is_Rational:
                coef *= Fraction(f.p, f.q)
            else:
                factors.append(f)
        if coef != 1 or len(factors) != 1:
            return None
        atoms = merge_atoms([parse_atom_exact(f) for f in factors])
        if len(atoms) != 1:
            return None
        return atoms[0]
    except (ValueError, TypeError):
        return None


def parse_problem(problem: str):
    """Split ``problem`` reusing decompose.py's parser pieces.

    Returns (comp, terms, rationals, rhs, lhs_expr): terms are (coef, atoms)
    pairs; rhs is Fraction for a numeric right side or a (kind, arg) atom.
    """
    txt = problem.replace(" ", "").replace("\n", "")
    txt = txt.replace("°", "*pi/180").replace("^", "**")
    if txt.count(">") + txt.count("<") != 1:
        raise ValueError("请输入一个只包含一个 > 或 < 的不等式")
    comp = ">" if ">" in txt else "<"
    lt_s, rt_s = txt.split(comp)
    if not lt_s.strip() or not rt_s.strip():
        raise ValueError("表达式不能为空")
    if BAD_CHARS.search(lt_s) or BAD_CHARS.search(rt_s):
        raise ValueError("表达式包含暂不支持的字符")
    lhs = sp.expand(sp.sympify(lt_s, locals=EXACT_LOC)).subs(SYMPY_CONST)
    rhs = sp.sympify(rt_s, locals=EXACT_LOC).subs(SYMPY_CONST)

    terms = []
    rationals = Fraction(0)
    for term in sp.Add.make_args(lhs):
        coef = Fraction(1)
        factors = []
        for f in sp.Mul.make_args(term):
            if f.is_Rational:
                coef *= Fraction(f.p, f.q)
            else:
                factors.append(f)
        if not factors:
            rationals += coef
            continue
        atoms = merge_atoms([parse_atom_exact(f) for f in factors])
        for kind, arg in atoms:
            check_domain(kind, arg)
        if len(atoms) > 1 and any(k in EXACT_FUNC_KINDS for k, _ in atoms):
            raise ValueError(ERR_ATOM)
        terms.append((coef, atoms))

    if rhs.is_Rational or rhs.is_Float:
        r = Fraction(rhs.p, rhs.q) if rhs.is_Rational else Fraction(str(rhs))
    else:
        r = single_atom_exact(rhs)
        if r is None:
            raise ValueError(ERR_RHS_FORM)
        check_domain(*r)
    return comp, terms, rationals, r, lhs, rhs


# ------------------------------------------------------------------ bounds


def need_lower(sub: Sub, dir_sign: int, outer_sign: int) -> bool:
    """Whether the term bound needs a lower bound on this factor's t-value."""
    return dir_sign * outer_sign > 0


def sub_comp(sub: Sub, dir_sign: int, outer_sign: int) -> str:
    """Direction to prove on U: lower bound on t needs U>β for σ=+1, U<β for σ=-1."""
    lower_t = need_lower(sub, dir_sign, outer_sign)
    return ">" if lower_t == (sub.sigma > 0) else "<"


def snap_t(t, margin, lower: bool) -> Fraction:
    """Rational bound on t at margin ∈ [margin/2, margin] on the needed side.

    lower=True picks ceil((t-margin)Q)/Q ∈ [t-margin, t-margin/2];
    lower=False picks floor((t+margin)Q)/Q ∈ [t+margin/2, t+margin].
    """
    m = mpf(margin) if isinstance(margin, Fraction) else margin
    edge = t - m if lower else t + m
    digits = max(2, int(-mp.log10(m)) + 3) if m > 0 else 62
    q = 10 ** min(digits, 90)
    n = mp.ceil(edge * q) if lower else mp.floor(edge * q)
    b = Fraction(int(n), q)
    # quantum safety: bound must stay strictly on its side of t
    if lower and mpf(b) >= t:
        b -= Fraction(1, q)
    if not lower and mpf(b) <= t:
        b += Fraction(1, q)
    return b


def place_bounds(term: TermAlloc, share: Fraction, dir_sign: int) -> None:
    """Choose every sub's comp/bound so term.bound lands within share of value."""
    outer_sign = 1 if term.outer > 0 else -1
    k = len(term.subs)
    v = abs(term.value())
    eta = float(share) / float(v) if v > 0 else 0.5
    eta = min(eta, 0.5)
    # split the relative margin ∝ relative floors, blended with a uniform
    # share so zero-floor factors still get a nonzero margin
    floors = [float(s.floor) / abs(float(s.u)) if s.u else 0.0 for s in term.subs]
    fsum = sum(floors)
    mean = fsum / k
    for s, f in zip(term.subs, floors, strict=True):
        rho = eta * (f + 0.25 * mean) / (1.25 * fsum) if fsum > 0 else eta / k
        rho = min(rho, 0.5)
        s.comp = sub_comp(s, dir_sign, outer_sign)
        t = s.t()
        t_abs = abs(t)
        lower_t = need_lower(s, dir_sign, outer_sign)
        b_t = snap_t(t, t_abs * rho, lower_t)
        if s.sigma > 0:
            s.bound = b_t
        else:
            # bound on U = 1/t: a lower bound on t maps to upper bound on U
            s.bound = 1 / b_t if b_t != 0 else Fraction(10) ** 30
    term.bound = term.outer
    for s in term.subs:
        term.bound *= s.phi()


# ------------------------------------------------------------------ frontier floor


def provable_at(sub: Sub, comp: str, bound: Fraction) -> bool:
    """Prove oracle: exact-mode search certifies ``U comp bound``?"""
    try:
        solve.prove(sub.kind, str(sub.power), comp, str(bound), exact=True)
    except Exception:
        return False
    return True


def margin_floor(sub: Sub, comp: str) -> Fraction:
    """Bisect the provability floor: smallest margin δ with U⋚(U∓δ) provable.

    Returns a Fraction within a factor ~2 of the true floor, or UNPROVABLE
    when no bound proves even at unit scale. Runs ~8 oracle calls.
    """
    u = sub.u
    lower = comp == ">"
    scale = abs(u) if u != 0 else mp.mpf(1)
    hi = scale
    if not provable_at(sub, comp, snap_t(u, hi / 2, lower)):
        hi = scale * 16
        if not provable_at(sub, comp, snap_t(u, hi, lower)):
            return UNPROVABLE
    lo = scale * mp.mpf(10) ** -20
    if provable_at(sub, comp, snap_t(u, lo, lower)):
        return Fraction(0)  # provable below 1e-20·|U|: no measurable floor
    while hi / lo > 8:
        mid = mp.sqrt(lo * hi)
        if provable_at(sub, comp, snap_t(u, mid, lower)):
            hi = mid
        else:
            lo = mid
    return frac60(hi)


# ------------------------------------------------------------------ driver


def certified_slack(terms: list[TermAlloc], comp: str, r: Fraction):
    """sign(Σv - r) certified at escalating precision; returns (sign, |slack|)."""
    want = 1 if comp == ">" else -1
    for dps in (80, DPS, 1200):
        with mp.workdps(dps):
            for s in (s for t in terms for s in t.subs):
                s.u = unit_value(s.kind, s.power)
            tot = mp.mpf(0)
            for t in terms:
                tot += t.value()
            diff = tot - mpf(r)
            guard = mp.mpf(2) ** (30 - dps) * max(1, abs(tot))
            if want * diff > guard:
                return want, abs(diff)
            if want * diff < -guard:
                return -want, abs(diff)
    return 0, mp.mpf(0)


def share_split(terms: list[TermAlloc], slack: Fraction, boost: dict) -> list[Fraction]:
    """Slack shares ∝ per-term frontier need v_i·Σ floor_j/|U_j|, boost-weighted."""
    needs = []
    for i, t in enumerate(terms):
        with mp.workdps(DPS):
            v = abs(t.value())
            need = mp.mpf(0)
            for s in t.subs:
                need += mpf(s.floor) / abs(s.u) if s.u else mp.mpf(0)
            need = v * need * mp.mpf(boost.get(i, 1.0))
        # damp toward equal split so zero-need terms still get slack
        needs.append(need)
    tot = sum(needs)
    if tot == 0:
        needs = [mp.mpf(1)] * len(terms)
        tot = mp.mpf(len(terms))
    damped = [n + tot / len(needs) * mp.mpf("0.02") for n in needs]
    dtot = sum(damped)
    return [slack * frac60(n / dtot) for n in damped]


def decompose_exact(problem: str) -> dict:
    """Decompose ``Σ c_i·Π atoms ⋚ R`` into provable exact-mode sub-claims.

    Returns ``{success, problem, comparison, steps, term_bounds, all_proved,
    certifies, failures, slack, decomposition_latex, normalized_latex}``.
    Every step carries the emitted (kind, power, comparison, bound), the
    prove parameters, achieved depth m+n, and the certified margin.
    """
    comp, terms_raw, rationals, r, lhs_expr, _rhs = parse_problem(problem)
    dir_sign = 1 if comp == ">" else -1

    if not isinstance(r, Fraction):
        return decompose_exact_pair(problem, comp, terms_raw, r)

    R = r - rationals  # lhs rationals fold into the target
    terms = []
    for coef, atoms in terms_raw:
        outer, subs = term_specs(coef, atoms)
        terms.append(TermAlloc(outer=outer, subs=subs, coef=coef, atoms=atoms))

    if not terms:
        ok = (Fraction(0) > R) if comp == ">" else (Fraction(0) < R)
        return {
            "success": ok,
            "problem": problem,
            "comparison": comp,
            "rhs": str(R),
            "steps": [],
            "term_bounds": [],
            "all_proved": ok,
            "certifies": ok,
            "failures": [] if ok else [{"reason": "false_rational_claim"}],
            "slack": str(abs(R)),
        }

    sign, slack_mpf = certified_slack(terms, comp, R)
    if sign != dir_sign:
        return {
            "success": False,
            "problem": problem,
            "comparison": comp,
            "rhs": str(R),
            "steps": [],
            "term_bounds": [],
            "all_proved": False,
            "certifies": False,
            "failures": [{"reason": "claim_certified_false"}],
            "slack": "0",
        }
    # 1% guard below the certified slack — 必须在 Fraction 侧做乘法：
    # slack_mpf * mp.mpf("0.99") 在模块默认 dps=15 下会把预算压回 float64
    slack = frac60(slack_mpf) * Fraction(99, 100)

    # ---- single-term single-atom: prove directly against R (loosest bound)
    if len(terms) == 1 and len(terms[0].subs) == 1 and terms[0].subs[0].sigma > 0:
        t = terms[0]
        s = t.subs[0]
        s.u = unit_value(s.kind, s.power)
        # term = outer·U ⋚ R  →  U ⋚ R/outer (direction flips with outer sign)
        s.comp = comp if t.outer > 0 else flip(comp)
        s.bound = R / t.outer
        run_one(s)
        t.bound = t.outer * s.phi()
        return assemble(problem, comp, R, terms, slack, direct=True, lhs=lhs_expr)

    # ---- frontier floors per distinct sub-claim (process-global cache:
    # a kernel's provability floor is deterministic per (kind,power,comp))
    for t in terms:
        osign = 1 if t.outer > 0 else -1
        for s in t.subs:
            c = sub_comp(s, dir_sign, osign)
            key = (s.kind, s.power, c)
            if key not in FLOOR_CACHE:
                with mp.workdps(DPS):
                    if s.u is None:
                        s.u = unit_value(s.kind, s.power)
                    FLOOR_CACHE[key] = margin_floor(s, c)
            s.floor = FLOOR_CACHE[key]
            s.comp = c

    boost: dict = {}
    proof_cache: dict = {}
    for _ in range(MAX_ROUNDS):
        shares = share_split(terms, slack, boost)
        for t, sh in zip(terms, shares, strict=True):
            place_bounds(t, sh, dir_sign)
        failed_terms = set()
        for i, t in enumerate(terms):
            for s in t.subs:
                key = (s.kind, s.power, s.comp, s.bound)
                if key in proof_cache:
                    s.proof, s.error = proof_cache[key]
                    continue
                run_one(s)
                proof_cache[key] = (s.proof, s.error)
                if s.proof is None:
                    failed_terms.add(i)
        if not failed_terms:
            break
        hard = False
        for i in failed_terms:
            t = terms[i]
            if any(s.error not in ("NoSolution",) for s in t.subs if s.proof is None):
                hard = True  # WrongDirection/ValueError: not fixable by margin
            else:
                boost[i] = boost.get(i, 1.0) * 4.0
        if hard:
            break

    return assemble(problem, comp, R, terms, slack, direct=False, lhs=lhs_expr)


def run_one(s: Sub) -> None:
    """Run one sub-claim through exact prove; record proof or error name."""
    try:
        s.proof = solve.prove(s.kind, str(s.power), s.comp, str(s.bound), exact=True)
        s.error = None
    except (NoSolution, WrongDirection, EqualClaim) as exc:
        s.proof, s.error = None, type(exc).__name__
    except Exception as exc:
        s.proof, s.error = None, f"{type(exc).__name__}: {exc}"


def assemble(problem, comp, R, terms, slack, direct, lhs=None) -> dict:
    """Build the result dict; run the ℚ certificate check Σbounds ⋚ R."""
    total = sum(t.bound for t in terms)
    certifies = (total >= R) if comp == ">" else (total <= R)
    all_proved = all(s.proof is not None for t in terms for s in t.subs)
    failures = []
    for i, t in enumerate(terms):
        for s in t.subs:
            if s.proof is None:
                failures.append(
                    {
                        "term": i,
                        "kind": s.kind,
                        "power": str(s.power),
                        "comparison": s.comp,
                        "bound": str(s.bound),
                        "error": s.error,
                    }
                )
    if not certifies:
        failures.append(
            {"reason": "allocated_bounds_do_not_cover_rhs", "sum": str(total), "rhs": str(R)}
        )
    steps = []
    for t in terms:
        for s in t.subs:
            step = {
                "type": s.kind,
                "power": str(s.power),
                "comparison": s.comp,
                "bound": str(s.bound),
                "bound_latex": frac(s.bound),
                "label": s.label,
                "coefficient": str(s.coef_disp if s.kind in EXACT_COEF_KINDS else t.coef),
            }
            if s.proof is not None:
                params = s.proof.get("parameters")
                if params is not None:
                    step["parameters"] = params
                    step["depth"] = int(params["m"]) + int(params["n"])
                    with suppress(Exception):
                        step["equation"] = render.render_equation(
                            params, s.kind, str(s.power), s.comp, str(s.bound)
                        )
                if s.proof.get("prover") in ("pade", "composite", "agm", "euler_gamma"):
                    step["prover"] = s.proof["prover"]  # non-(m,n) certificates
                step["margin"] = str(abs(frac60(s.u) - s.bound))
            else:
                step["error"] = s.error
            steps.append(step)
    return {
        "success": all_proved and certifies,
        "problem": problem,
        "comparison": comp,
        "rhs": str(R),
        "slack": str(slack),
        "steps": steps,
        "term_bounds": [str(t.bound) for t in terms],
        "bound_sum": str(total),
        "all_proved": all_proved,
        "certifies": certifies,
        "direct": direct,
        "failures": failures,
        "decomposition_latex": decomposition_latex(terms, comp, R, total),
        "normalized_latex": "" if lhs is None else sp.latex(lhs),
    }


def decompose_exact_pair(problem: str, comp: str, terms_raw: list, ratom: tuple) -> dict:
    """Constant-vs-constant ``atom ⋚ atom``: chain two proofs via a rational midpoint."""
    if not (len(terms_raw) == 1 and len(terms_raw[0][1]) == 1 and terms_raw[0][0] == 1):
        raise ValueError(ERR_RHS_FORM)
    lk, la = terms_raw[0][1][0]
    rk, ra = ratom
    # arg<0 reaches here only for exponent kinds (reciprocal — pair mode
    # cannot model 1/U) or signed-arg kinds (a plain negative q, fine)
    if (la < 0 and lk not in SIGNED_ARG_KINDS) or (ra < 0 and rk not in SIGNED_ARG_KINDS):
        raise ValueError(ERR_RHS_FORM)
    lp, rp = la, ra  # atom arg IS the power slot for coef and arg kinds alike
    with mp.workdps(DPS):
        lv = unit_value(lk, lp)
        rv = unit_value(rk, rp)
        mid = frac60((lv + rv) / 2)
        inside = min(lv, rv) < mpf(mid) < max(lv, rv)
        if not inside:
            mid = frac60((lv * 3 + rv) / 4)
    if comp == ">" and not lv > rv:
        return {"success": False, "all_proved": False, "failures": [{"reason": ERR_NUMERIC}]}
    if comp == "<" and not lv < rv:
        return {"success": False, "all_proved": False, "failures": [{"reason": ERR_NUMERIC}]}
    steps = []
    ok = True
    for kind, power, c in ((lk, lp, comp), (rk, rp, flip(comp))):
        s = Sub(
            kind=kind,
            power=power,
            sigma=1,
            arg=power,
            label=atom_label_exact(kind, power),
            coef_disp=Fraction(1),
            comp=c,
            bound=mid,
        )
        run_one(s)
        ok = ok and s.proof is not None
        step = {
            "type": kind,
            "power": str(power),
            "comparison": c,
            "bound": str(mid),
            "label": atom_label_exact(kind, power),
        }
        if s.proof is not None:
            if s.proof.get("parameters") is not None:
                step["parameters"] = s.proof["parameters"]
            if s.proof.get("prover") in ("pade", "composite", "agm", "euler_gamma"):
                step["prover"] = s.proof["prover"]
        else:
            step["error"] = s.error
        steps.append(step)
    return {
        "success": ok,
        "problem": problem,
        "comparison": comp,
        "steps": steps,
        "all_proved": ok,
        "certifies": ok,
        "failures": [] if ok else [{"reason": "subproof_failed"}],
    }


def bound_piece(t: TermAlloc) -> str:
    """Rendered bound for one term: |outer|·Π factor bounds."""
    parts = []
    if abs(t.outer) != 1:
        parts.append(frac(abs(t.outer)))
    for s in t.subs:
        b = frac(abs(s.bound))
        parts.append(b if s.sigma > 0 else f"\\dfrac{{1}}{{{b}}}")
    return "\\cdot ".join(parts) if parts else "0"


def decomposition_latex(terms: list[TermAlloc], comp: str, r: Fraction, total) -> str:
    """``lhs ⋚ bound_sum = total [⋚ R]`` mirroring the site's field."""
    lhs_pieces, bound_pieces = [], []
    for t in terms:
        lhs_pieces.append(term_piece_exact(t.coef, t.atoms))
        piece = bound_piece(t)
        bound_pieces.append("-" + piece if t.bound < 0 else piece)
    dl = f"{join_signed(lhs_pieces)}{comp}{join_signed(bound_pieces)}={frac(total)}"
    if total != r:
        dl += f"{comp}{frac(r)}"
    return dl
