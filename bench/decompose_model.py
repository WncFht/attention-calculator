"""Deterministic model of /decompose_inequality bound allocation.

Term-level semantics
--------------------
Each term has a signed true value v_i = c_i * prod(atom values).
'<' needs an upper bound u_i on v_i; '>' needs a lower bound l_i.
The step's `bound` field is always on |coef|*atom: b_i = |u_i| resp |l_i|.

Processing (sums, i.e. more than one term):
  * resid term: '<' -> last negative-coef term in steps order, else first
    steps-order term; '>' -> first steps-order term.  With product terms
    present, '<' resid is the first NON-product steps-order term instead.
  * non-resid terms are processed wrap-around: steps order after the resid
    first, then steps order before it, resid last.
  * share_i = R - sum(contrib_j), contrib_j = bound if assigned else v_j.

Bound picks ('<'):
  positive single-atom:  u = largest upper record of |v_i| <= share
  negative single-atom:  u = -b, b = smallest lower record of |v_i| >= -share
  product (non-resid):   u = largest product-value record <= share,
                         else verbatim share
  resid: u = R - sum; positive -> largest upper record <= u (product: same
  then factor split), negative -> atom bound = -u verbatim.
'>' mirrors with floor records (strict > for processed positives, <= for
processed negatives); resid verbatim for positive coef, largest upper record
<= -l for negative coef.

Product factor split (term bound B on |c|*A1*A2):
  record factor = lowest FACTOR_RANK atom; cap = B/(|c|*other trues);
  largest record <= cap ('<'); tail factor = B/(|c|*prior bounds) verbatim.
Reciprocal products (a negative-power factor): '>' gives the numerator a
smallest floor record > R*denom_true and the denominator a largest upper
record STRICTLY < num_b/R; both are records.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction

import sympy as sp

# ------------------------------------------------------------------- parsing

VARPI = sp.Symbol("varpi")
GAUSS = sp.Symbol("gauss")
# undefined functions so arctan(1)/ln(1) don't auto-evaluate to pi/4 / 0
LN = sp.Function("ln")
ATAN = sp.Function("arctan")
LOC = {"varpi": VARPI, "gauss": GAUSS, "ln": LN, "e": sp.E, "pi": sp.pi,
       "phi": sp.GoldenRatio, "zeta": sp.zeta, "catalan": sp.Catalan,
       "gamma": sp.EulerGamma, "arctan": ATAN, "atan": ATAN,
       "sinh": sp.sinh, "tanh": sp.tanh, "sin": sp.sin, "cos": sp.cos,
       "tan": sp.tan}

BASE = {sp.pi: ("pi", Fraction(1)), sp.E: ("e", Fraction(1)),
        sp.EulerGamma: ("gamma", Fraction(1)),
        sp.GoldenRatio: ("golden", Fraction(1)),
        sp.Catalan: ("catalan", Fraction(1)), sp.zeta(3): ("zeta3", Fraction(1)),
        sp.exp(sp.pi): ("e_pi", Fraction(1))}


def parse_atom(expr):
    """Map a sympy factor to (kind, arg)."""
    if expr in BASE:
        return BASE[expr]
    if expr == VARPI:
        return "varpi", Fraction(1)
    if expr == GAUSS:
        return "gauss", Fraction(1)
    if expr.is_Pow:
        b, e = expr.base, expr.exp
        if b == sp.E and e == sp.pi:
            return "e_pi", Fraction(1)
        if b in BASE and e.is_Rational:
            return {"pi": "pi_n", "e": "e_q"}[BASE[b][0]], Fraction(e.p, e.q)
        if b in (VARPI, GAUSS):
            raise ValueError("varpi/gauss powers unsupported")
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
        raise ValueError(f"bad arg {expr!r}")
    raise ValueError(f"unparsed atom {expr!r}")


def merge_atoms(atoms):
    """pi*pi -> pi_n(2), e*e^2 -> e_q(3): combine same-family powers."""
    fam = {"pi": "pi", "pi_n": "pi", "e": "e", "e_q": "e"}
    out = []
    for kind, arg in atoms:
        if out and fam.get(out[-1][0]) == fam.get(kind) == "pi" \
                and kind != "e_pi":
            a = out.pop()[1]
            out.append(("pi_n", a + arg))
        elif out and fam.get(out[-1][0]) == fam.get(kind) == "e":
            a = out.pop()[1]
            out.append(("e_q", a + arg))
        else:
            out.append((kind, arg))
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
        if not factors:
            # pure rational term on the lhs folds into R (-pi+8>4 -> -pi>-4)
            R -= coef
            continue
        terms.append((coef, merge_atoms([parse_atom(f) for f in factors])))
    return terms, comp, R


# --------------------------------------------------------------- true values

VARPI_V = float((sp.gamma(sp.Rational(1, 4)) ** 2
                 / (2 * sp.sqrt(2 * sp.pi))).evalf(40))
GAUSS_V = float((sp.gamma(sp.Rational(1, 4)) ** 2
                 / (2 * sp.sqrt(2 * sp.pi ** 3))).evalf(40))


def atom_value(kind: str, arg: Fraction) -> float:
    f = float(arg)
    return {
        "pi": lambda: math.pi * f,
        "e": lambda: math.e * f,
        "pi_n": lambda: math.pi ** f,
        "e_q": lambda: math.exp(f),
        "e_pi": lambda: math.e ** math.pi * f,
        "ln_q": lambda: math.log(f),
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


def term_true(coef, atoms) -> float:
    v = float(coef)
    for kind, arg in atoms:
        v *= atom_value(kind, arg)
    return v


# ------------------------------------------------------------- bound chains

class Chain:
    """Record bounds of a value: upper = running min of ceil(vK)/K,
    lower = running max of floor(vK)/K, K in [k_min, k_max]."""

    def __init__(self, value: float, k_min: int, k_max: int = 3000):
        up, lo = [], []
        best_u, best_l = math.inf, -math.inf
        for k in range(k_min, k_max + 1):
            n = math.floor(value * k)
            # float overshoot guard: n/k must stay below value
            if Fraction(n, k) > Fraction(str(value)):
                n -= 1
            u = Fraction(n + 1, k)
            if u < best_u:
                up.append(u)
                best_u = u
            lo_ = Fraction(n, k)
            if best_l < lo_ < value:
                lo.append(lo_)
                best_l = lo_
        self.upper = up
        self.lower = lo

    def largest_upper(self, cap: Fraction):
        ok = [b for b in self.upper if b <= cap]
        return ok[0] if ok else None

    def largest_upper_strict(self, cap: Fraction):
        ok = [b for b in self.upper if b < cap]
        return ok[0] if ok else None

    def smallest_lower(self, floor: Fraction):
        ok = [b for b in self.lower if b > floor]
        return ok[0] if ok else None

    def smallest_lower_geq(self, floor: Fraction):
        ok = [b for b in self.lower if b >= floor]
        return ok[0] if ok else None


K_MIN_UP = {"pi": 2, "e": 1, "gamma": 3, "golden": 3, "ln_q": 4, "sin_q": 1,
            "cos_q": 5, "tan_q": 3, "arctan_q": 5, "zeta3": 9, "catalan": 1,
            "sinh_q": 1, "tanh_q": 1, "varpi": 1, "gauss": 1, "e_q": 1,
            "pi_n": 1, "e_pi": 1}
# lower-record chains may start earlier (gamma floor 1/2 needs K=2)
K_MIN_LO = dict(K_MIN_UP, gamma=2)
CHAINS: dict = {}


def chain(kind: str, arg: Fraction, side: str) -> Chain:
    key = (kind, arg, side)
    if key not in CHAINS:
        kmin = K_MIN_UP.get(kind, 1) if side == "up" else K_MIN_LO.get(kind, 1)
        CHAINS[key] = Chain(atom_value(kind, arg), kmin)
    return CHAINS[key]


def val_chain(value: float, side: str) -> Chain:
    """Record chain on an arbitrary float (product values, |coef|*atom)."""
    key = ("v", round(value, 15), side)
    if key not in CHAINS:
        CHAINS[key] = Chain(value, 1)
    return CHAINS[key]


# ------------------------------------------------------------------ the model

# record-factor priority inside products: e-family first, then pi, golden, pi_n
FACTOR_RANK = {"e": 0, "e_q": 0, "e_pi": 0, "pi": 1, "golden": 2,
               "gamma": 4, "ln_q": 5, "pi_n": 3, "zeta3": 9, "catalan": 8,
               "sin_q": 6, "cos_q": 6, "tan_q": 7, "arctan_q": 7,
               "sinh_q": 7, "tanh_q": 7, "varpi": 10, "gauss": 10}


def pick_processed(terms, ti, share, comp, done, flat_done):
    """Bound for a processed (non-resid) term. Returns signed term bound."""
    coef, atoms = terms[ti]
    v = term_true(coef, atoms)
    if len(atoms) == 1:
        kind, arg = atoms[0]
        vv = abs(v)
        if coef == 1:
            ch_up = chain(kind, arg, "up")
            ch_lo = chain(kind, arg, "lo")
        else:
            # chain on the |coef|*atom value directly (bound is on c*C)
            ch_up = val_chain(vv, "up")
            ch_lo = val_chain(vv, "lo")
        if comp == "<":
            if coef > 0:
                return ch_up.largest_upper(share)
            b = ch_lo.smallest_lower_geq(-share)
            return -b if b is not None else None
        else:
            if coef > 0:
                return ch_lo.smallest_lower(share)
            b = ch_up.largest_upper(-share)
            return -b if b is not None else None
    # product term
    pch_up = val_chain(abs(v), "up")
    pch_lo = val_chain(abs(v), "lo")
    if comp == "<":
        if coef > 0:
            return pch_up.largest_upper(share) or share
        b = pch_lo.smallest_lower_geq(-share)
        return -b if b is not None else -share
    else:
        if coef > 0:
            return pch_lo.smallest_lower(share) or share
        b = pch_up.largest_upper(-share)
        return -b if b is not None else -share


def resid_bound(terms, ti, resid, comp, has_product):
    """Signed term bound for the resid slot.

    '<': verbatim resid; pure positive-coef sums (no product) snap the resid
    to the largest upper record <= resid instead.
    '>': verbatim for positive coef; negative coef takes the largest upper
    record <= -resid on the atom.
    """
    coef, atoms = terms[ti]
    if comp == "<":
        if coef > 0 and not has_product:
            if len(atoms) == 1 and coef == 1:
                b = chain(atoms[0][0], atoms[0][1], "up").largest_upper(resid)
            else:
                b = val_chain(abs(term_true(coef, atoms)),
                              "up").largest_upper(resid)
            if b is not None:
                return b
        return resid
    else:
        if coef >= 0:
            return resid  # verbatim
        if len(atoms) == 1 and coef == -1:
            b = chain(atoms[0][0], atoms[0][1], "up").largest_upper(-resid)
        else:
            b = val_chain(abs(term_true(coef, atoms)),
                          "up").largest_upper(-resid)
        if b is not None:
            return -b
        return resid


def factor_split(terms, ti, bound_signed, comp):
    """Split a product term's bound into per-atom bounds.

    Record factor = lowest FACTOR_RANK atom; cap = B/(|c|*other trues);
    tail factor = B/(|c|*prior bounds) verbatim.  Returns atom_index -> bound.
    """
    coef, atoms = terms[ti]
    n = len(atoms)
    trues = [atom_value(k, a) for k, a in atoms]
    order = sorted(range(n), key=lambda j: FACTOR_RANK.get(atoms[j][0], 50))
    out = {}
    for pos, j in enumerate(order):
        kind, arg = atoms[j]
        prod_others = Fraction(1)
        for k in range(n):
            if k == j:
                continue
            prod_others *= (out[k] if k in out
                            else Fraction(str(trues[k])))
        cap = bound_signed / (abs(coef) * prod_others)
        if pos < n - 1:
            if comp == "<":
                b = chain(kind, arg, "up").largest_upper(cap)
            else:
                b = chain(kind, arg, "lo").smallest_lower(cap)
            if b is None:
                b = cap
        else:
            b = cap
        out[j] = b
    return out


def reciprocal_split(terms, ti, comp, R, first_j):
    """a * b^{-1} vs R: both factors take records, processed in display order.

    '>' num-first: num = smallest floor > R*den_base_true, then
    denom = largest upper STRICTLY < num_b/R.
    '>' denom-first: denom = largest upper <= num_true/R, then
    num = smallest floor STRICTLY > R*denom_b.
    ('<' mirrors.)
    """
    coef, atoms = terms[ti]
    num_j = next(j for j, (k, a) in enumerate(atoms) if a > 0)
    den_j = next(j for j, (k, a) in enumerate(atoms) if a < 0)
    nk, na = atoms[num_j]
    dk, da = atoms[den_j]
    den_true = Fraction(str(atom_value(dk, -da)))  # bound is on base const
    num_true = Fraction(str(atom_value(nk, na)))
    out = {}
    if comp == ">":
        if first_j == num_j:
            nb = chain(nk, na, "lo").smallest_lower(R * den_true / abs(coef))
            db = chain(dk, -da, "up").largest_upper_strict(
                nb * abs(coef) / R) if nb else None
        else:
            db = chain(dk, -da, "up").largest_upper(num_true * abs(coef) / R)
            nb = chain(nk, na, "lo").smallest_lower(
                R * db / abs(coef)) if db else None
    else:
        if first_j == num_j:
            nb = chain(nk, na, "up").largest_upper(R * den_true / abs(coef))
            db = chain(dk, -da, "lo").smallest_lower(
                nb * abs(coef) / R) if nb else None
        else:
            db = chain(dk, -da, "lo").smallest_lower(num_true * abs(coef) / R)
            nb = chain(nk, na, "up").largest_upper(
                R * db / abs(coef)) if db else None
    if nb is None or db is None:
        return None
    out[num_j], out[den_j] = nb, db
    return out


def is_reciprocal(atoms):
    return any(a < 0 for _, a in atoms)


BASE_OF = {"pi_n": "pi", "e_q": "e"}


def match_atom(kind, arg, coef, site, in_product):
    """Site step (type, coefficient) vs our (kind, arg, term coef).

    Single-atom terms: coefficient = |term coef|.  Product factors:
    coefficient = |arg| (the factor's power); reciprocal factors report
    the BASE type (e^{-1} -> 'e').
    """
    want = site["type"]
    wantc = Fraction(site["coefficient"])
    if in_product:
        if arg < 0:
            return BASE_OF.get(kind, kind) == want and wantc == -arg
        return want == kind and wantc == abs(arg)
    return want == kind and wantc == abs(coef)


def choose_resid(terms, term_order, comp):
    """Index (into term_order list) of the resid term."""
    has_product = any(len(terms[t][1]) > 1 for t in term_order)
    if comp == "<":
        if has_product:
            # first non-product term in steps order
            for pos, t in enumerate(term_order):
                if len(terms[t][1]) == 1:
                    return pos
            return 0
        negs = [pos for pos, t in enumerate(term_order) if terms[t][0] < 0]
        if negs:
            return negs[-1]
        return 0
    return 0


def predict(problem, resp_steps, verbose=False):
    """Predict per-step bounds for a site record.

    Returns (matched: bool, predicted list aligned to resp_steps).
    """
    terms, comp, R = parse_problem(problem)
    # align site steps to (term, atom) flat indices
    flat = [(i, j) for i, (c, a) in enumerate(terms) for j in range(len(a))]
    used = [False] * len(flat)
    order_flat = []
    import re
    for s in resp_steps:
        m = re.search(r"left\((-?\d+(?:/\d+)?)\\right\)|\\l?n?\((\d+)\)",
                      s["label"])
        arg_hint = Fraction(m.group(1) or m.group(2)) if m else None
        cand = [i for i in range(len(flat)) if not used[i]
                and match_atom(terms[flat[i][0]][1][flat[i][1]][0],
                               terms[flat[i][0]][1][flat[i][1]][1],
                               terms[flat[i][0]][0], s,
                               len(terms[flat[i][0]][1]) > 1)
                and (arg_hint is None
                     or abs(terms[flat[i][0]][1][flat[i][1]][1])
                     == arg_hint)]
        if not cand:
            return None, [("NOMATCH", s["type"], s["coefficient"],
                           s["label"])]
        # ambiguous same-type atoms: site lists product factors first
        cand.sort(key=lambda i: (len(terms[flat[i][0]][1]) == 1,
                                 abs(terms[flat[i][0]][1][flat[i][1]][1])))
        used[cand[0]] = True
        order_flat.append(cand[0])
    term_order = list(dict.fromkeys(flat[i][0] for i in order_flat))

    # ---- bound allocation
    n = len(terms)
    if n == 1:
        ti = term_order[0]
        coef, atoms = terms[ti]
        if is_reciprocal(atoms):
            done = {ti: None}
        elif comp == "<":
            done = {ti: R}  # single-term '<': bound = R verbatim
        else:
            done = {ti: resid_bound(terms, ti, R, comp, False)}
    else:
        resid_pos = choose_resid(terms, term_order, comp)
        resid_ti = term_order[resid_pos]
        has_product = any(len(terms[t][1]) > 1 for t in range(n))
        # processing order = steps order with the resid slot removed
        proc = term_order[:resid_pos] + term_order[resid_pos + 1:]
        done = {}

        def contrib(tj):
            if tj in done:
                return done[tj]
            return Fraction(str(term_true(*terms[tj])))

        ok = True
        for ti in proc:
            share = R - sum(contrib(tj) for tj in range(n) if tj != ti)
            b = pick_processed(terms, ti, share, comp, done, None)
            if b is None:
                ok = False
                break
            done[ti] = b
        if ok:
            share_r = R - sum(contrib(tj) for tj in range(n)
                              if tj != resid_ti)
            done[resid_ti] = resid_bound(terms, resid_ti, share_r, comp,
                                         has_product)
        if not ok:
            return False, [("INVALID",)]

    # ---- map term bounds to per-step bounds
    pred = [None] * len(order_flat)
    for ti in range(n):
        coef, atoms = terms[ti]
        if len(atoms) == 1:
            bterm = done[ti]
            b_atom = abs(bterm) if coef != 0 else bterm
            # bound is on |coef|*atom: for negative term b = -signed
            b_atom = bterm if coef > 0 else -bterm
            for idx, fi in enumerate(order_flat):
                if flat[fi][0] == ti:
                    pred[idx] = b_atom
        else:
            if is_reciprocal(atoms):
                first_j = next(flat[fi][1] for fi in order_flat
                               if flat[fi][0] == ti)
                fsplit = reciprocal_split(terms, ti, comp, R, first_j)
            elif done[ti] is None:
                fsplit = None
            else:
                fsplit = factor_split(terms, ti, abs(done[ti]), comp)
            if fsplit is None:
                return False, [("SPLITFAIL",)]
            for idx, fi in enumerate(order_flat):
                t, j = flat[fi]
                if t == ti:
                    pred[idx] = fsplit[j]
    site = [Fraction(s["bound"]) for s in resp_steps]
    matched = pred == site
    if verbose or not matched:
        detail = [(terms[flat[order_flat[i]][0]][1][flat[order_flat[i]][1]][0],
                   str(p), str(s))
                  for i, (p, s) in enumerate(zip(pred, site, strict=False))]
        return matched, detail
    return matched, pred


def load(path="bench/data/decompose.jsonl"):
    with open(path) as fh:
        for line in fh:
            yield json.loads(line)


def report():
    hits = miss = 0
    fails = []
    for r in load():
        resp = r.get("response") or {}
        if not resp.get("success"):
            continue
        try:
            matched, detail = predict(r["problem"], resp["steps"])
        except Exception as e:
            matched, detail = None, [("EXC", repr(e))]
        if matched:
            hits += 1
        else:
            miss += 1
            fails.append((r["problem"], detail))
    print(f"hit={hits} miss={miss}")
    for prob, detail in fails:
        print("MISS", prob, detail)


if __name__ == "__main__":
    report()
