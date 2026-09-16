"""mode=exact composite-inequality decomposition: provable-by-construction bounds.

Site counterpart (decompose.py) allocates record-chain bounds chosen for byte
parity with zhuyidao.net — no provability guarantee. This module allocates
rational sub-bounds so that (a) the conjunction certifies the claim and
(b) every sub-claim is provable within the (m,n) search budget, verified by
actually running ``solve.prove(..., exact=True)`` on each.

Model — every sub-claim certifies ``U ⋚ β`` where ``U = constant_mpf(kind,
power)`` is the unit the kernel's certified_cmp evaluates:

- single-atom term ``c·a``: coefficient kinds (pi, e, gamma, golden, catalan,
  zeta3, varpi, gauss, e_pi) fold |c| into the power slot — U = |c|·C and the
  term is sign(c)·U; argument kinds (pi_n, e_q, ln_q, sin_q, ...) keep the
  coefficient outside — U = a^sign(arg), term = c·U^sign(arg).
- product term ``c·Π a_j``: U_j = base value constant_mpf(kind_j, |arg_j|),
  term = c·Π U_j^sign(arg_j); negative exponents are reciprocal factors and
  the sub-proof runs on the base (e_q arg -2 proves e² ⋚ β, contributing 1/β).

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
    COEF_KINDS,
    ERR_ATOM,
    ERR_NUMERIC,
    ERR_RHS_FORM,
    FUNC_KINDS,
    LOC,
    atom_label,
    flip,
    frac,
    join_signed,
    merge_atoms,
    parse_atom,
    single_atom,
    term_piece,
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
    if len(atoms) == 1 and atoms[0][0] in COEF_KINDS:
        kind, arg = atoms[0]
        sub = Sub(
            kind=kind,
            power=abs(coef),
            sigma=1,
            arg=arg,
            label=atom_label(kind, arg),
            coef_disp=abs(coef),
        )
        return Fraction(-1 if coef < 0 else 1), [sub]
    subs = []
    for kind, arg in atoms:
        subs.append(
            Sub(
                kind=kind,
                power=abs(arg),
                sigma=1 if arg > 0 else -1,
                arg=arg,
                label=atom_label(kind, abs(arg)),
                coef_disp=Fraction(1),
            )
        )
    return coef, subs


# ------------------------------------------------------------------ parsing


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
    lhs = sp.expand(sp.sympify(lt_s, locals=LOC))
    rhs = sp.sympify(rt_s, locals=LOC)

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
        atoms = merge_atoms([parse_atom(f) for f in factors])
        if len(atoms) > 1 and any(k in FUNC_KINDS for k, _ in atoms):
            raise ValueError(ERR_ATOM)
        terms.append((coef, atoms))

    if rhs.is_Rational or rhs.is_Float:
        r = Fraction(rhs.p, rhs.q) if rhs.is_Rational else Fraction(str(rhs))
    else:
        r = single_atom(rhs)
        if r is None:
            raise ValueError(ERR_RHS_FORM)
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
                "coefficient": str(s.coef_disp if s.kind in COEF_KINDS else t.coef),
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
                if s.proof.get("prover") == "pade":
                    step["prover"] = "pade"  # Padé certificate: no (m,n) params
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
    if la < 0 or ra < 0:
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
            label=atom_label(kind, power),
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
            "label": atom_label(kind, power),
        }
        if s.proof is not None:
            if s.proof.get("parameters") is not None:
                step["parameters"] = s.proof["parameters"]
            if s.proof.get("prover") == "pade":
                step["prover"] = "pade"
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
        lhs_pieces.append(term_piece(t.coef, t.atoms))
        piece = bound_piece(t)
        bound_pieces.append("-" + piece if t.bound < 0 else piece)
    dl = f"{join_signed(lhs_pieces)}{comp}{join_signed(bound_pieces)}={frac(total)}"
    if total != r:
        dl += f"{comp}{frac(r)}"
    return dl
