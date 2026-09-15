"""Sweep round-2: more numeric-solve hypotheses for the site's sign check.

Key new mechanisms:
  f64sys: moment entries -> float64 -> Fraction, solve EXACTLY, sign check
          on exact coeffs of the perturbed system. Bound exact or f64.
  mp{15..40}: mpmath lu_solve at dps.
  exactnum: stored Cramer formula (P + Q*B)/D with numerator evaluated
          EXACTLY on exact B, then float(num)/float(D).
  cramer_f: (float(P) + float(Q)*Bf)/float(D)   [from round 1, refuted]
"""

import sys
from fractions import Fraction

from mpmath import mp

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from sim_evals import cramer_terms, system
from sim_float_solve import f_nonneg, f_nonpos
from sim_mechanisms import basis_and_target

from attention_calculator.engine import mn_order, solve_moment


def solve_exact_on(basis, target):
    try:
        return solve_moment(basis, target)
    except ValueError:
        return None


def mech_f64sys_rexact(basis, target, B):
    """Round moments to f64, keep bound exact, solve exactly."""
    fb = [{k: Fraction(float(v)) for k, v in m.items()} for m in basis]
    return solve_exact_on(fb, target)


def mech_f64sys_rf(basis, target, B):
    """Round moments to f64 AND bound to f64, solve exactly."""
    fb = [{k: Fraction(float(v)) for k, v in m.items()} for m in basis]
    t2 = dict(target)
    if "1" in t2:
        t2["1"] = Fraction(float(B)) * (-1 if t2["1"] < 0 else 1)
    return solve_exact_on(fb, t2)


def _forms(basis, target):
    """(P,Q,D) Cramer terms; t1 = d(rhs)/dB = sign(target['1']) on '1'."""
    keys, M, t0, _ = system(basis, target, 1)
    bsign = 1 if target.get("1", 0) > 0 else -1
    t1 = [(Fraction(bsign) if k == "1" else Fraction(0)) for k in keys]
    return cramer_terms(M, t0, t1)


def mech_exactnum(basis, target, B):
    """Stored formula (P + Q*B)/D: numerator exact on exact B, then
    float(num)/float(D)."""
    ct = _forms(basis, target)
    if ct is None:
        return None
    P, Q, D = ct
    return [float(p + q * B) / float(D) for p, q in zip(P, Q, strict=True)]


def mech_exactnum_fden(basis, target, B):
    """Exact rational u then float: float((P + Q*B)/D) — i.e. f64coeff."""
    ct = _forms(basis, target)
    if ct is None:
        return None
    P, Q, D = ct
    return [float((p + q * B) / D) for p, q in zip(P, Q, strict=True)]


def mech_mp(basis, target, B, dps):
    mp.dps = dps
    keys = sorted(set(target) | {k for m in basis for k in m})
    n = len(keys)
    A = mp.matrix(n, n)
    b = mp.matrix(n, 1)
    for i, k in enumerate(keys):
        for j, m in enumerate(basis):
            v = m.get(k, Fraction(0))
            A[i, j] = mp.mpf(v.numerator) / v.denominator
        v = target.get(k, Fraction(0))
        b[i] = mp.mpf(v.numerator) / v.denominator
    try:
        u = mp.lu_solve(A, b)
    except Exception:
        return None
    return [float(u[i]) for i in range(n)]


MECHS = {
    "f64sys_r": mech_f64sys_rexact,
    "f64sys_f": mech_f64sys_rf,
    "exnum": mech_exactnum,
    "exnum2": mech_exactnum_fden,
    "mp16": lambda b, t, B: mech_mp(b, t, B, 16),
    "mp20": lambda b, t, B: mech_mp(b, t, B, 20),
    "mp25": lambda b, t, B: mech_mp(b, t, B, 25),
    "mp30": lambda b, t, B: mech_mp(b, t, B, 30),
    "mp35": lambda b, t, B: mech_mp(b, t, B, 35),
}

CASES = [
    ("pi  >tight", "pi", "1", ">", "6134899525417045/1952799169684491", "WD"),
    ("pi  <float", "pi", "1", "<", "884279719003555/281474976710656", "NS"),
    ("pi  >float", "pi", "1", ">", "884279719003555/281474976710656", "OK(12,20)"),
    ("pi  <tight", "pi", "1", "<", "5706674932067741/1816491048114374", "NS"),
    ("e   >tight", "e", "1", ">", "2124008553358849/781379079653017", "WD"),
    ("e   <float", "e", "1", "<", "6121026514868073/2251799813685248", "NS"),
    ("e   >float", "e", "1", ">", "6121026514868073/2251799813685248", "OK(7,7)"),
    ("e   <tight", "e", "1", "<", "9581113603440437/3524694718233772", "OK(11,12)"),
    ("e_q >tight", "e_q", "1", ">", "2124008553358849/781379079653017", "WD"),
    ("e_q3>tight", "e_q", "3", ">", "9060589063489840/451100167157089", "NS"),
    ("sin >tight", "sin_q", "1", ">", "8199121568317184/9743795943467975", "WD"),
    ("cos >tight", "cos_q", "1", ">", "293104830616638/542483027433479", "NS"),
    ("ln  >tight", "ln_q", "2", ">", "1554903831458736/2243252046704767", "WD"),
    ("lnsq>tight", "ln_q_square", "2", ">", "4646620020445232/9671330777074349", "WD"),
    ("zeta>tight", "zeta3", "1", ">", "5517909911604193/4590389936699688", "NS"),
    ("atan>tight", "arctan_q", "3", ">", "125014145208443/100087721339793", "NS"),
    ("pi_n>tight", "pi_n", "5/2", ">", "4923959516357971/281474976710656", "NS"),
    ("sin >float", "sin_q", "1", ">", "3789648413623927/4503599627370496", "OK"),
]


def predict(mech, kind, pw, comp, rs):
    if isinstance(mech, str):
        mech = MECHS[mech]
    bound = Fraction(rs)
    limit, mk, target, B = basis_and_target(kind, pw, comp, bound)
    for m, n in mn_order(limit):
        basis = mk(m, n)
        try:
            u = mech(basis, target, B)
        except Exception:
            u = None
        if u is None:
            continue
        if f_nonneg(u):
            return f"OK({m},{n})"
        if comp == ">" and f_nonpos(u):
            return f"WD({m},{n})"
    return "NS"


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else None
    names = list(MECHS)
    print(f"{'case':12s} {'site':10s}" + "".join(f"{n:>12s}" for n in names))
    score = dict.fromkeys(names, 0)
    for name, k, p, c, r, want in CASES:
        if which and which not in name:
            continue
        row = f"{name:12s} {want:10s}"
        for mech in names:
            got = predict(mech, k, p, c, r)
            ok = ((want.startswith("OK") and got.startswith("OK"))
                  or (want == "WD" and got.startswith("WD"))
                  or (want == "NS" and got == "NS"))
            score[mech] += ok
            row += f"{got:>12s}"
        print(row, flush=True)
    print("\nscores:", score)
