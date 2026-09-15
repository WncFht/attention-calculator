"""Simulate the site's hypothesized float64 moment solve + float sign check.

For each (m,n) in site order: build the moment matrix with float64 entries,
solve u * basis = target in float64 (Gaussian elim, partial pivot), then apply
poly_nonneg on the float coefficients. Record the first candidate whose float
verdict is nonneg (site would 200) or nonpos (site would 方向反了, '>' scan),
and compare against the live-site outcome stored in edge-probes.jsonl.
"""

import sys
from fractions import Fraction

sys.path.insert(0, "src")
from attention_calculator.engine import mn_order
from attention_calculator.kernels import exp_family, log_family, quadlog, trig_q


def gauss_float(rows, rhs):
    """Solve A x = rhs in float64; A is list-of-rows (square). None if singular."""
    n = len(rhs)
    a = [[*r[:], b] for r, b in zip(rows, rhs, strict=True)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if a[piv][col] == 0.0:
            return None
        a[col], a[piv] = a[piv], a[col]
        for r in range(n):
            if r != col and a[r][col] != 0.0:
                f = a[r][col] / a[col][col]
                a[r] = [x - f * y for x, y in zip(a[r], a[col], strict=True)]
    return [a[i][n] / a[i][i] for i in range(n)]


def solve_float(basis: list[dict], target: dict) -> list[float] | None:
    """Float64 solve of sum_j u_j basis[j] == target over the symbol space."""
    keys = sorted(set(target) | {k for m in basis for k in m})
    rows = [[float(m.get(k, 0.0)) for m in basis] for k in keys]
    rhs = [float(target.get(k, 0.0)) for k in keys]
    return gauss_float(rows, rhs)


def f_nonneg(c: list[float]) -> bool:
    a = c[0]
    b = c[1] if len(c) > 1 else 0.0
    cc = c[2] if len(c) > 2 else 0.0
    if cc <= 0:
        return a >= 0 and a + b + cc >= 0
    if b < 0 < b + 2 * cc:
        return 4 * a * cc - b * b >= 0
    return a >= 0 and a + b + cc >= 0


def f_nonpos(c):
    return f_nonneg([-v for v in c])


def scan(kind, power_s, comp, bound_s, limit, mk_basis, target, mid_nonpos=True):
    """Walk (m,n); return (event, m, n, coeffs) where event in
    {'nonneg','nonpos','exhaust'}."""
    for i, (m, n) in enumerate(mn_order(limit)):
        basis = mk_basis(m, n)
        u = solve_float(basis, target)
        if u is None:
            continue
        if f_nonneg(u):
            return ("nonneg", m, n, u, i)
        if mid_nonpos and f_nonpos(u):
            return ("nonpos", m, n, u, i)
    return ("exhaust", None, None, None, None)


def case_quadlog(kind, power_s, comp, bound_s):
    power, bound = Fraction(power_s), Fraction(bound_s)
    cfg = quadlog.spec(kind, power)
    q, r_, odd, sym = cfg["q"], cfg["r"], cfg["odd"], cfg["sym"]
    if sym == "arctan":
        def term(k):
            return quadlog.atan_moment(k, q)
    else:
        factor = cfg["factor"]
        def term(k):
            cc, rat = quadlog.ln_moment(k, r_, odd)
            return cc * factor, rat
    sign = 1 if comp == ">" else -1
    bb = bound ** cfg.get("pd", 1)
    target = {sym: sign * cfg["coef"], "1": -sign * bb}
    def mk(m, n):
        return quadlog.basis_moments(m, n, odd, sym, term)
    return scan(kind, power_s, comp, bound_s, cfg["limit"], mk, target)


def case_exp(kind, power_s, comp, bound_s):
    power, bound = Fraction(power_s), Fraction(bound_s)
    if kind == "e_q":
        sym, coef, q, limit = "e_q", Fraction(1), power, 10
    else:
        sym, coef, q, limit = "e", power, Fraction(1), 30
    sign = 1 if comp == ">" else -1
    target = {sym: sign * coef, "1": -sign * bound}
    def mk(m, n):
        return [exp_family.basis_x_moment(m, n, j, q, sym) for j in (0, 1)]
    return scan(kind, power_s, comp, bound_s, limit, mk, target)


def case_log(kind, power_s, comp, bound_s):
    power, bound = Fraction(power_s), Fraction(bound_s)
    qt = log_family.qtilde(kind, power)
    c = qt - 1
    square = kind == "ln_q_square"
    sign = Fraction(1 if comp == ">" else -1)
    target = {"ln2": sign, "1": -sign * bound} if square else {"ln": sign, "1": -sign * bound}
    def mk(m, n):
        return [log_family.basis_moment(c, max(m, n, 1), m, n, j, square)
                for j in range(3 if square else 2)]
    return scan(kind, power_s, comp, bound_s, 10, mk, target)


def case_trig(kind, power_s, comp, bound_s):
    q, bound = Fraction(power_s), Fraction(bound_s)
    s = Fraction(1 if comp == ">" else -1)
    target = {"sin_q": s, "1": -s * bound} if kind == "sin_q" else {"cos_q": s, "1": -s * bound}
    moments = trig_q.sin_moments(22, q)
    def mk(m, n):
        return [trig_q.basis_moment(moments, m, n, j) for j in range(3)]
    return scan(kind, power_s, comp, bound_s, 10, mk, target)


CASES = [
    ("pi-gt-tight  site=方向反了", case_quadlog,
     ("pi", "1", ">", "6134899525417045/1952799169684491")),
    ("pi-float-lt  site=未找到", case_quadlog, ("pi", "1", "<", "884279719003555/281474976710656")),
    ("pi-float-gt  site=200@(12,20)", case_quadlog,
     ("pi", "1", ">", "884279719003555/281474976710656")),
    ("pi-lt-tight  site=未找到", case_quadlog,
     ("pi", "1", "<", "5706674932067741/1816491048114374")),
    ("e-gt-tight   site=方向反了", case_exp, ("e", "1", ">", "2124008553358849/781379079653017")),
    ("e-float-lt   site=未找到", case_exp, ("e", "1", "<", "6121026514868073/2251799813685248")),
    ("e-float-gt   site=200@(7,7)", case_exp, ("e", "1", ">", "6121026514868073/2251799813685248")),
    ("e_q-tight    site=方向反了", case_exp, ("e_q", "1", ">", "2124008553358849/781379079653017")),
    ("e_q3-tight   site=未找到", case_exp, ("e_q", "3", ">", "9060589063489840/451100167157089")),
    ("sin_q-tight  site=方向反了", case_trig,
     ("sin_q", "1", ">", "8199121568317184/9743795943467975")),
    ("cos_q-tight  site=未找到", case_trig, ("cos_q", "1", ">", "293104830616638/542483027433479")),
    ("ln_q-tight   site=方向反了", case_log,
     ("ln_q", "2", ">", "1554903831458736/2243252046704767")),
    ("ln_qsq-tight site=方向反了", case_log,
     ("ln_q_square", "2", ">", "4646620020445232/9671330777074349")),
    ("zeta3-tight  site=未找到", case_quadlog,
     ("zeta3", "1", ">", "5517909911604193/4590389936699688")),
    ("atan3-tight  site=未找到", case_quadlog,
     ("arctan_q", "3", ">", "125014145208443/100087721339793")),
]

if __name__ == "__main__":
    for name, fn, args in CASES:
        try:
            ev, m, n, u, i = fn(*args)
            print(f"{name}: {ev} at ({m},{n}) plan#{i} u={u}")
        except Exception as e:
            print(f"{name}: EXC {e}")
