"""Trace the (m,n) scan for divergent edge cases: exact vs float64 sign checks.

For each request we rebuild the same plans/target the kernel's prove() builds,
then for every candidate print the solved P coefficients, the exact-rational
verdict, and the verdict a float64-coefficient check would give.
"""

import sys
from fractions import Fraction

sys.path.insert(0, "src")
from attention_calculator.engine import mn_order, poly_nonneg, solve_moment
from attention_calculator.kernels import quadlog, exp_family, log_family, trig_q


def f64_nonneg(coeffs):
    """poly_nonneg on float64(coeff)."""
    return poly_nonneg([float(c) for c in coeffs])


def f64_nonpos(coeffs):
    return poly_nonneg([-float(c) for c in coeffs])


def verdict(coeffs):
    nn = poly_nonneg(coeffs)
    np_ = poly_nonneg([-c for c in coeffs])
    f_nn = f64_nonneg(coeffs)
    f_np = f64_nonpos(coeffs)
    e = "nonneg" if nn else ("nonpos" if np_ else "indef")
    f = "nonneg" if f_nn else ("nonpos" if f_np else "indef")
    return e, f


def plans_for(kind, power, comp):
    """Rebuild the (m,n,basis) plan iterator a kernel's prove() would use."""
    if kind in ("pi", "pi_n", "catalan", "zeta3", "arctan_q", "arccot_q"):
        cfg = quadlog.spec(kind, power)
        q, r, odd, sym = cfg["q"], cfg["r"], cfg["odd"], cfg["sym"]
        if sym == "arctan":
            def term(k):
                return quadlog.atan_moment(k, q)
        else:
            factor = cfg["factor"]
            def term(k):
                cc, rat = quadlog.ln_moment(k, r, odd)
                return cc * factor, rat
        plans = ((m, n, quadlog.basis_moments(m, n, odd, sym, term))
                 for m, n in mn_order(cfg["limit"]))
        sign = 1 if comp == ">" else -1
        bound_sign = sign
        target = {sym: sign * cfg["coef"], "1": -sign * Fraction(0)}  # placeholder
        return plans, cfg, sign
    raise ValueError(kind)


def trace_quadlog(kind, power_s, comp, bound_s, maxplans=400):
    power, bound = Fraction(power_s), Fraction(bound_s)
    cfg = quadlog.spec(kind, power)
    q, r, odd, sym = cfg["q"], cfg["r"], cfg["odd"], cfg["sym"]
    if sym == "arctan":
        def term(k):
            return quadlog.atan_moment(k, q)
    else:
        factor = cfg["factor"]
        def term(k):
            cc, rat = quadlog.ln_moment(k, r, odd)
            return cc * factor, rat
    sign = 1 if comp == ">" else -1
    bb = bound ** cfg.get("pd", 1)
    target = {sym: sign * cfg["coef"], "1": -sign * bb}
    print(f"== {kind} {power_s} {comp} {bound_s}  (limit {cfg['limit']})")
    for i, (m, n) in enumerate(mn_order(cfg["limit"])):
        if i >= maxplans:
            print("  ... truncated")
            break
        basis = quadlog.basis_moments(m, n, odd, sym, term)
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            print(f"  ({m},{n}) singular")
            continue
        e, f = verdict(coeffs)
        mark = "" if e == f else "   <-- DIVERGES"
        if e != "indef" or f != "indef":
            print(f"  ({m},{n}) exact={e:7s} float={f:7s} coeffs={[str(c) for c in coeffs]}{mark}")


def trace_exp(kind, power_s, comp, bound_s, maxplans=2000):
    power, bound = Fraction(power_s), Fraction(bound_s)
    if kind == "e_q":
        sym, coef, q, limit = "e_q", Fraction(1), power, exp_family.LIMIT_OTHER
        mk = lambda m, n: [exp_family.basis_x_moment(m, n, j, q, sym) for j in (0, 1)]
    else:
        sym, coef, q, limit = "e", power, Fraction(1), exp_family.LIMIT_E
        mk = lambda m, n: [exp_family.basis_x_moment(m, n, j, q, sym) for j in (0, 1)]
    sign = 1 if comp == ">" else -1
    target = {sym: sign * coef, "1": -sign * bound}
    print(f"== {kind} {power_s} {comp} {bound_s}  (limit {limit})")
    for i, (m, n) in enumerate(mn_order(limit)):
        if i >= maxplans:
            print("  ... truncated")
            break
        basis = mk(m, n)
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            continue
        e, f = verdict(coeffs)
        if e != "indef" or f != "indef":
            mark = "" if e == f else "   <-- DIVERGES"
            print(f"  ({m},{n}) exact={e:7s} float={f:7s} coeffs={[str(c) for c in coeffs]}{mark}")


def trace_log(kind, power_s, comp, bound_s, maxplans=300):
    power, bound = Fraction(power_s), Fraction(bound_s)
    qt = log_family.qtilde(kind, power)
    c = qt - 1
    square = kind == "ln_q_square"
    sign = Fraction(1 if comp == ">" else -1)
    if square:
        target = {"ln2": sign, "1": -sign * bound}
    else:
        target = {"ln": sign, "1": -sign * bound}
    print(f"== {kind} {power_s} {comp} {bound_s}  (limit 10)")
    for i, (m, n) in enumerate(mn_order(10)):
        basis = [log_family.basis_moment(c, max(m, n, 1), m, n, j, square)
                 for j in range(3 if square else 2)]
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            continue
        e, f = verdict(coeffs)
        if e != "indef" or f != "indef":
            mark = "" if e == f else "   <-- DIVERGES"
            print(f"  ({m},{n}) exact={e:7s} float={f:7s} coeffs={[str(c2) for c2 in coeffs]}{mark}")


def trace_trig(kind, power_s, comp, bound_s):
    q, bound = Fraction(power_s), Fraction(bound_s)
    s = Fraction(1 if comp == ">" else -1)
    if kind == "sin_q":
        target = {"sin_q": s, "1": -s * bound}
    elif kind == "cos_q":
        target = {"cos_q": s, "1": -s * bound}
    moments = trig_q.sin_moments(2 * 10 + 2, q)
    print(f"== {kind} {power_s} {comp} {bound_s}  (limit 10)")
    for m, n in mn_order(10):
        basis = [trig_q.basis_moment(moments, m, n, j) for j in range(3)]
        try:
            coeffs = solve_moment(basis, target)
        except ValueError:
            continue
        e, f = verdict(coeffs)
        if e != "indef" or f != "indef":
            mark = "" if e == f else "   <-- DIVERGES"
            print(f"  ({m},{n}) exact={e:7s} float={f:7s} coeffs={[str(c2) for c2 in coeffs]}{mark}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("pi-false", "all"):
        trace_quadlog("pi", "1", "<", "884279719003555/281474976710656")
    if which in ("pi-tight", "all"):
        trace_quadlog("pi", "1", ">", "6134899525417045/1952799169684491")
    if which in ("e-false", "all"):
        trace_exp("e", "1", "<", "6121026514868073/2251799813685248")
    if which in ("e-tight", "all"):
        trace_exp("e", "1", ">", "2124008553358849/781379079653017")
    if which in ("eq-tight", "all"):
        trace_exp("e_q", "1", ">", "2124008553358849/781379079653017")
    if which in ("sin-tight", "all"):
        trace_trig("sin_q", "1", ">", "8199121568317184/9743795943467975")
    if which in ("ln-tight", "all"):
        trace_log("ln_q", "2", ">", "1554903831458736/2243252046704767")
    if which in ("lnsq-tight", "all"):
        trace_log("ln_q_square", "2", ">", "4646620020445232/9671330777074349")
