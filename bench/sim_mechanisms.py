"""Compare candidate numeric mechanisms for the site's sign check.

Mechanisms per candidate (m,n):
  A) float64 Gauss-Jordan solve of the moment system (sim_float_solve)
  B) affine-eval: u(r) = A + B*r exact rationals; sign-check on
     float(A) + float(B)*float(r)
  C) mpmath solve at various dps
  D) exact solve + float64 eval of P at sample points
  E) float64 Cramer for 2x2 / det formulas

For each divergent case print the first sign-definite verdict under each
mechanism and compare with the observed site outcome.
"""

import sys
from fractions import Fraction

from mpmath import mp

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from attention_calculator.engine import mn_order, solve_moment, poly_nonneg
from attention_calculator.kernels import quadlog, exp_family, log_family, trig_q
from sim_float_solve import gauss_float, f_nonneg, f_nonpos


def basis_and_target(kind, power_s, comp, bound_s):
    """Return (limit, mk_basis(m,n), target, r) replicating the kernel setup."""
    power, bound = Fraction(power_s), Fraction(bound_s)
    if kind in ("pi", "pi_n", "catalan", "zeta3", "arctan_q", "arccot_q"):
        cfg = quadlog.spec(kind, power)
        q, rr, odd, sym = cfg["q"], cfg["r"], cfg["odd"], cfg["sym"]
        if sym == "arctan":
            def term(k):
                return quadlog.atan_moment(k, q)
        else:
            factor = cfg["factor"]
            def term(k):
                cc, rat = quadlog.ln_moment(k, rr, odd)
                return cc * factor, rat
        sign = 1 if comp == ">" else -1
        bb = bound ** cfg.get("pd", 1)
        target = {sym: sign * cfg["coef"], "1": -sign * bb}
        return cfg["limit"], (lambda m, n: quadlog.basis_moments(m, n, odd, sym, term)), target, bb
    if kind in ("e", "e_q"):
        if kind == "e_q":
            sym, coef, q, limit = "e_q", Fraction(1), power, 10
        else:
            sym, coef, q, limit = "e", power, Fraction(1), 30
        sign = 1 if comp == ">" else -1
        target = {sym: sign * coef, "1": -sign * bound}
        mk = lambda m, n: [exp_family.basis_x_moment(m, n, j, q, sym) for j in (0, 1)]
        return limit, mk, target, bound
    if kind in ("ln_q", "ln_q_square"):
        qt = log_family.qtilde(kind, power)
        c = qt - 1
        square = kind == "ln_q_square"
        sign = Fraction(1 if comp == ">" else -1)
        target = ({"ln2": sign, "1": -sign * bound} if square
                  else {"ln": sign, "1": -sign * bound})
        mk = lambda m, n: [log_family.basis_moment(c, max(m, n, 1), m, n, j, square)
                           for j in range(3 if square else 2)]
        return 10, mk, target, bound
    if kind in ("sin_q", "cos_q"):
        s = Fraction(1 if comp == ">" else -1)
        target = ({"sin_q": s, "1": -s * bound} if kind == "sin_q"
                  else {"cos_q": s, "1": -s * bound})
        moments = trig_q.sin_moments(22, power)
        mk = lambda m, n: [trig_q.basis_moment(moments, m, n, j) for j in range(3)]
        return 10, mk, target, bound
    raise ValueError(kind)


def mech_gauss(basis, target, rf):
    keys = sorted(set(target) | {k for m in basis for k in m})
    t2 = dict(target)
    rows = [[float(m.get(k, 0.0)) for m in basis] for k in keys]
    rhs = [float(t2.get(k, 0.0)) for k in keys]
    return gauss_float(rows, rhs)


def mech_affine(basis, target, rf):
    """Solve u(r) exactly for r symbolic, then float-eval A + B*rf.

    u = M^{-1} rhs; rhs = t0 + r*t1 where t0 has the symbol part and t1 picks
    the bound coefficient. Compute u0 = solve(rhs with r=0) and du = solve
    for the r-coefficient: rhs_r = coefficient of r in rhs.
    """
    keys = sorted(set(target) | {k for m in basis for k in m})
    rows = [[m.get(k, Fraction(0)) for m in basis] for k in keys]
    # rhs entries: only the "1" component carries r (value -sign*r or bound^pd)
    rhs0, rhsr = [], []
    for k in keys:
        v = target.get(k, Fraction(0))
        if k == "1":
            rhs0.append(Fraction(0))
            rhsr.append(v)  # v = -sign (bound multiplies it)
        else:
            rhs0.append(v)
            rhsr.append(Fraction(0))
    try:
        u0 = solve_moment(basis, dict(zip(keys, rhs0)))
        ur = solve_moment(basis, dict(zip(keys, rhsr)))
    except ValueError:
        return None
    return [float(a) + float(b) * rf for a, b in zip(u0, ur)]


def mech_mp(basis, target, rf, dps):
    mp.dps = dps
    keys = sorted(set(target) | {k for m in basis for k in m})
    n = len(keys)
    t2 = {k: (mp.mpf(v.numerator) / v.denominator if isinstance(v, Fraction) else mp.mpf(v))
          for k, v in target.items()}
    if "1" in t2:  # bound component must use rf (float64) — same value anyway
        pass
    A = mp.matrix(n, n)
    b = mp.matrix(n, 1)
    for i, k in enumerate(keys):
        for j, m in enumerate(basis):
            v = m.get(k, Fraction(0))
            A[i, j] = mp.mpf(v.numerator) / v.denominator
        b[i] = t2.get(k, mp.mpf(0))
    try:
        u = mp.lu_solve(A, b)
    except Exception:
        return None
    return [float(u[i]) for i in range(n)]


def mech_exact_then_sample(basis, target, rf):
    """Exact solve; sign check = P evaluated at float points {0, .25, .5, .75, 1}."""
    try:
        u = solve_moment(basis, target)
    except ValueError:
        return None
    cf = [float(c) for c in u]
    vals = []
    for t in (0.0, 0.25, 0.5, 0.75, 1.0):
        vals.append(cf[0] + (cf[1] * t if len(cf) > 1 else 0)
                    + (cf[2] * t * t if len(cf) > 2 else 0))
    if all(v >= 0 for v in vals):
        return ("nonneg", u)
    if all(v <= 0 for v in vals):
        return ("nonpos", u)
    return ("indef", u)


def mech_exact_f64coeff(basis, target, rf):
    try:
        u = solve_moment(basis, target)
    except ValueError:
        return None
    return [float(c) for c in u]


MECHS = {
    "gauss": mech_gauss,
    "affine": mech_affine,
    "mp15": lambda b, t, r: mech_mp(b, t, r, 15),
    "mp25": lambda b, t, r: mech_mp(b, t, r, 25),
    "mp40": lambda b, t, r: mech_mp(b, t, r, 40),
    "f64coeff": mech_exact_f64coeff,
}

CASES = [
    ("pi  >tight   方向反了", ("pi", "1", ">", "6134899525417045/1952799169684491")),
    ("pi  <float   未找到", ("pi", "1", "<", "884279719003555/281474976710656")),
    ("pi  >float   200@12,20", ("pi", "1", ">", "884279719003555/281474976710656")),
    ("e   >tight   方向反了", ("e", "1", ">", "2124008553358849/781379079653017")),
    ("e   <float   未找到", ("e", "1", "<", "6121026514868073/2251799813685248")),
    ("e   >float   200@7,7", ("e", "1", ">", "6121026514868073/2251799813685248")),
    ("e_q >tight   方向反了", ("e_q", "1", ">", "2124008553358849/781379079653017")),
    ("e_q3>tight   未找到", ("e_q", "3", ">", "9060589063489840/451100167157089")),
    ("sin >tight   方向反了", ("sin_q", "1", ">", "8199121568317184/9743795943467975")),
    ("cos >tight   未找到", ("cos_q", "1", ">", "293104830616638/542483027433479")),
    ("ln  >tight   方向反了", ("ln_q", "2", ">", "1554903831458736/2243252046704767")),
    ("lnsq>tight   方向反了", ("ln_q_square", "2", ">", "4646620020445232/9671330777074349")),
    ("zeta>tight   未找到", ("zeta3", "1", ">", "5517909911604193/4590389936699688")),
    ("atan>tight   未找到", ("arctan_q", "3", ">", "125014145208443/100087721339793")),
]


def verdict_of(u):
    if u is None:
        return "sing"
    if isinstance(u, tuple):
        return u[0]
    if f_nonneg(u):
        return "nonneg"
    if f_nonpos(u):
        return "nonpos"
    return "indef"


def first_hits(kind, pw, comp, rs):
    limit, mk, target, r = basis_and_target(kind, pw, comp, rs)
    rf = float(r)
    out = {}
    for name, mech in MECHS.items():
        res = ("exhaust", None, None)
        for i, (m, n) in enumerate(mn_order(limit)):
            basis = mk(m, n)
            try:
                u = mech(basis, target, rf)
            except Exception:
                u = None
            v = verdict_of(u)
            if v in ("nonneg", "nonpos"):
                res = (v, (m, n), i)
                break
        out[name] = res
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else None
    for name, (k, p, c, r) in CASES:
        if which and which not in name:
            continue
        print(f"### {name}")
        hits = first_hits(k, p, c, r)
        for mech, res in hits.items():
            print(f"   {mech:8s}: {res}")
