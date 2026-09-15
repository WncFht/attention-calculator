"""Reverse-engineering harness for /decompose_inequality.

For each bench record we know the site's per-step bounds. A model must
reproduce them. Search space:

  * processing order: permutation of the output steps,
  * resid slot: which step takes the leftover instead of a chain bound,
  * resid rule: verbatim resid vs largest record <= resid,
  * product share: how a product term's share is computed,
  * product tail: verbatim quotient vs largest record.

The script scores every combination over all records and reports which
combinations match which records, so the surviving rule set shows itself.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction

import sympy as sp

# ------------------------------------------------------------------- parsing

VARPI = sp.Symbol("varpi")
GAUSS = sp.Symbol("gauss")
LOC = {"varpi": VARPI, "gauss": GAUSS, "ln": sp.log, "e": sp.E, "pi": sp.pi,
       "phi": sp.GoldenRatio, "zeta": sp.zeta, "catalan": sp.Catalan,
       "gamma": sp.EulerGamma, "arctan": sp.Function("arctan"),
       "atan": sp.Function("arctan"), "sinh": sp.sinh,
       "tanh": sp.tanh, "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
       "atan": sp.atan}

BASE = {sp.pi: ("pi", Fraction(1)), sp.E: ("e", Fraction(1)),
        sp.EulerGamma: ("gamma", Fraction(1)),
        sp.GoldenRatio: ("golden", Fraction(1)),
        sp.Catalan: ("catalan", Fraction(1)), sp.zeta(3): ("zeta3", Fraction(1)),
        sp.exp(sp.pi): ("e_pi", Fraction(1))}


def parse_atom(expr) -> tuple[str, Fraction]:
    """Map a sympy factor to (kind, arg)."""
    if expr in BASE:
        return BASE[expr]
    if expr == VARPI:
        return "varpi", Fraction(1)
    if expr == GAUSS:
        return "gauss", Fraction(1)
    if expr.is_Pow:
        b, e = expr.base, expr.exp
        if b in BASE and e.is_Rational:
            return {"pi": "pi_n", "e": "e_q"}[BASE[b][0]], Fraction(e.p, e.q)
    if expr.is_Function:
        name = expr.func.__name__
        a = expr.args[0]
        if name == "exp":
            if a == sp.pi:
                return "e_pi", Fraction(1)
            if a.is_Rational:
                return "e_q", Fraction(a.p, a.q)
        if a.is_Rational:
            arg = Fraction(a.p, a.q)
            table = {"ln": "ln_q", "log": "ln_q", "sin": "sin_q",
                     "cos": "cos_q", "tan": "tan_q", "atan": "arctan_q",
                     "arctan": "arctan_q", "sinh": "sinh_q", "tanh": "tanh_q"}
            if name in table:
                return table[name], arg
    raise ValueError(f"unparsed atom {expr!r}")


def merge_atoms(atoms):
    """pi*pi -> pi_n(2), e*e^2 -> e_q(3): combine same-family powers."""
    fam = {"pi": "pi", "pi_n": "pi", "e": "e", "e_q": "e"}
    out = []
    for kind, arg in atoms:
        if out and fam.get(out[-1][0]) == fam.get(kind) == "pi":
            k, a = out.pop()
            out.append(("pi_n", a + arg))
        elif out and fam.get(out[-1][0]) == fam.get(kind) == "e":
            k, a = out.pop()
            out.append(("e_q", a + arg))
        else:
            out.append((kind, arg))
    # normalize (pi,1)->('pi',1) fine; pi_n(1)->pi
    out = [("pi", a) if k == "pi_n" and a == 1 else (k, a) for k, a in out]
    out = [("e", a) if k == "e_q" and a == 1 else (k, a) for k, a in out]
    return out


def parse_problem(problem: str):
    """-> (terms, comp, R); terms = [(coef, [(kind,arg), ...]), ...]."""
    txt = problem.replace("°", "*pi/180").replace("zeta3", "zeta(3)")
    txt = txt.replace("^", "**")
    for comp in ("<", ">"):
        if comp in txt:
            lt, rt = txt.split(comp)
            lhs = sp.sympify(lt, locals=LOC)
            rhs = sp.sympify(rt, locals=LOC)
            break
    else:
        raise ValueError(f"not relational: {problem}")
    if rhs.is_Rational:
        R = Fraction(rhs.p, rhs.q)
    elif rhs.is_Float:
        R = Fraction(str(rhs))
    elif rhs.is_Number:
        raise ValueError(f"non-numeric rhs: {problem}")
    else:
        # symbolic rhs: fold into lhs (pi>e -> pi-e>0)
        lhs = lhs - rhs
        R = Fraction(0)
    terms = []
    for term in sp.Add.make_args(lhs):
        coef = Fraction(1)
        factors = []
        for f in sp.Mul.make_args(term):
            if f.is_Rational:
                coef *= Fraction(f.p, f.q)
            else:
                factors.append(f)
        terms.append((coef, merge_atoms([parse_atom(f) for f in factors])))
    return terms, comp, R


# --------------------------------------------------------------- true values

VARPI_V = float((sp.gamma(sp.Rational(1, 4)) ** 2
                 / (2 * sp.sqrt(2 * sp.pi))).evalf(40))
GAUSS_V = float((sp.gamma(sp.Rational(1, 4)) ** 2
                 / (2 * sp.pi ** sp.Rational(3, 2))).evalf(40))


def atom_value(kind: str, arg: Fraction) -> float:
    """Float64 value of one atom (bound applies to coef*atom or the atom
    inside a product; powers already folded into arg)."""
    f = float(arg)
    return {
        "pi": lambda: math.pi * f,
        "e": lambda: math.e * f,
        "pi_n": lambda: math.pi ** f,
        "e_q": lambda: math.exp(f),
        "e_pi": lambda: math.e ** math.pi * f,
        "ln_q": lambda: math.log(f),
        "ln_q_square": lambda: math.log(f) ** 2,
        "sin_q": lambda: math.sin(f),
        "cos_q": lambda: math.cos(f),
        "tan_q": lambda: math.tan(f),
        "arctan_q": lambda: math.atan(f),
        "sinh_q": lambda: math.sinh(f),
        "tanh_q": lambda: math.tanh(f),
        "gamma": lambda: f * 0.5772156649015329,
        "golden": lambda: f * (1 + math.sqrt(5)) / 2,
        "catalan": lambda: 0.915965594177219 * f,
        "zeta3": lambda: 1.2020569031595942 * f,
        "varpi": lambda: VARPI_V * f,
        "gauss": lambda: GAUSS_V * f,
    }[kind]()


def term_true(coef: Fraction, atoms) -> float:
    v = float(coef)
    for kind, arg in atoms:
        v *= atom_value(kind, arg)
    return v


# ------------------------------------------------------------- bound chains

class Chain:
    """Upper/lower record bounds of ceil/floor(t*K)/K, K in [k_min, k_max]."""

    def __init__(self, value: float, k_min: int, k_max: int = 3000):
        self.value = value
        up, lo = [], []
        best_u, best_l = math.inf, -math.inf
        for k in range(k_min, k_max + 1):
            n = math.floor(value * k)
            if Fraction(n, k) > Fraction(str(value)):
                n -= 1  # float overshoot guard
            u = Fraction(n + 1, k)
            if u < best_u:
                up.append(u)
                best_u = u
            l = Fraction(n, k)
            if best_l < l < value:
                lo.append(l)
                best_l = l
        self.upper = up
        self.lower = lo

    def largest_upper(self, cap: Fraction) -> Fraction | None:
        ok = [b for b in self.upper if b <= cap]
        return ok[0] if ok else None

    def smallest_lower(self, floor: Fraction) -> Fraction | None:
        ok = [b for b in self.lower if b > floor]
        return ok[-1] if ok else None


K_MIN = {"pi": 2, "e": 1, "gamma": 3, "golden": 3, "ln_q": 4, "sin_q": 1,
         "cos_q": 5, "tan_q": 3, "arctan_q": 5, "zeta3": 9, "catalan": 1,
         "sinh_q": 1, "tanh_q": 1, "varpi": 1, "gauss": 1, "e_q": 1,
         "pi_n": 1, "e_pi": 1}
CHAINS: dict = {}


def chain(kind: str, arg: Fraction) -> Chain:
    key = (kind, arg)
    if key not in CHAINS:
        CHAINS[key] = Chain(atom_value(kind, arg), K_MIN.get(kind, 1))
    return CHAINS[key]


# ---------------------------------------------------------------- brute force

def step_atoms(terms):
    """flat list of (term_idx, atom_idx) in parsed order."""
    return [(i, j) for i, (c, a) in enumerate(terms)
            for j in range(len(a))]


def match_flat(terms, resp_steps):
    """Align site steps (type/coef/label order) to flat atom indices.

    Returns list of flat indices in site order, or None if mismatched.
    """
    flat = step_atoms(terms)
    kinds = [terms[i][1][j][0] for i, j in flat]
    coefs = [terms[i][0] for i, j in flat]
    args = [terms[i][1][j][1] for i, j in flat]
    used = [False] * len(flat)
    order = []
    import re
    for s in resp_steps:
        want = s["type"]
        wantc = Fraction(s["coefficient"])
        # argument hint from label, e.g. \sin\left(2\right) or \ln(3)
        m = re.search(r"\((-?\d+(?:/\d+)?)\s*\\?right?\)|\{(?:\\!)?(-?\d+)\}",
                      s["label"])
        arg_hint = Fraction(m.group(1) or m.group(2)) if m else None
        cand = [i for i in range(len(flat)) if not used[i]
                and kinds[i] == want and coefs[i] == wantc
                and (arg_hint is None or args[i] == arg_hint)]
        if not cand:
            cand = [i for i in range(len(flat)) if not used[i]
                    and kinds[i] == want]
        if not cand:
            return None
        # if several candidates remain (same type+coef), prefer arg order
        # ascending for determinism
        cand.sort(key=lambda i: args[i])
        used[cand[0]] = True
        order.append(cand[0])
    return order


def fit_record(terms, comp, R, site_bounds, resid_rules=("record", "verbatim")):
    """Try every (processing order, resid slot, resid rule). Return the list
    of (order, resid_pos, rule) that reproduce the site bounds exactly.

    Bound semantics:
      '<': each atom bound = largest upper record <= share/rest  (share =
           R - assigned bounds - unprocessed trues, divided by coef);
           resid slot -> verbatim resid or largest record.
      '>': mirrored on lower records; negative-coef terms use upper records.
    """
    flat = step_atoms(terms)
    n = len(flat)
    hits = []

    def termres(tj):
        """current contribution of term tj: product of bounds if all its
        steps are assigned, else its true value."""
        c, atoms = terms[tj]
        idxs = [flat.index((tj, j)) for j in range(len(atoms))]
        if all(i in done for i in idxs):
            v = Fraction(c)
            for i in idxs:
                v *= done[i]
            return v
        return Fraction(str(term_true(c, atoms)))

    for order in itertools.permutations(range(n)):
        for resid_pos in range(n):
            for rr in resid_rules:
                done: dict[int, Fraction] = {}
                okk = True
                for pos, si in enumerate(order):
                    ti, aj = flat[si]
                    coef, atoms = terms[ti]
                    share = R - sum(termres(tj) for tj in range(len(terms))
                                    if tj != ti)
                    kind, arg = atoms[aj]
                    ch = chain(kind, arg)
                    if len(atoms) == 1:
                        cap = share / coef
                    else:
                        rest = Fraction(coef)
                        for k, (k2, a2) in enumerate(atoms):
                            if k == aj:
                                continue
                            # other factors: bound if assigned else true
                            si2 = flat.index((ti, k))
                            rest *= (done[si2] if si2 in done else
                                     Fraction(str(atom_value(k2, a2))))
                        cap = share / rest
                    if pos == resid_pos and rr == "verbatim":
                        b = cap
                    else:
                        # direction-aware record pick
                        neg = (comp == ">") == (coef > 0)
                        if not neg:
                            b = ch.largest_upper(cap)
                        else:
                            b = ch.smallest_lower(cap)
                        if b is None:
                            okk = False
                            break
                    done[si] = b
                if okk and all(done[i] == site_bounds[i] for i in range(n)):
                    hits.append((order, resid_pos, rr))
    return hits


def load():
    out = []
    for line in open("bench/data/decompose.jsonl"):
        r = json.loads(line)
        out.append((r["problem"], r.get("response") or {}))
    return out


def report():
    recs = load()
    hits_n = miss_n = 0
    for prob, resp in recs:
        if not resp.get("success"):
            continue
        try:
            terms, comp, R = parse_problem(prob)
        except ValueError as e:
            print("PARSEFAIL", prob, e)
            continue
        site = [Fraction(s["bound"]) for s in resp["steps"]]
        order = match_flat(terms, resp["steps"])
        flat = step_atoms(terms)
        if order is None:
            print("NOMATCH", prob)
            continue
        # site bounds keyed by flat index
        keyed = [Fraction(0)] * len(flat)
        for si, fi in enumerate(order):
            keyed[fi] = site[si]
        hits = fit_record(terms, comp, R, keyed)
        stypes = [flat_kind(terms, fi) for fi in order]
        if hits:
            hits_n += 1
            flag = "HIT "
        else:
            miss_n += 1
            flag = "MISS"
        print(f"{flag} {prob:32s} site={list(zip(stypes, map(str, site)))}")
        if hits:
            print("     orders:", hits[:4], "..." if len(hits) > 4 else "")
    print(f"hit={hits_n} miss={miss_n}")


def flat_kind(terms, si):
    flat = step_atoms(terms)
    ti, aj = flat[si]
    return terms[ti][1][aj][0]

if __name__ == "__main__":
    report()
