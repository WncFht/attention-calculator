"""Systematic sweep of float-eval strategies for the site's sign check.

For each candidate (m,n), the solved u(B) is affine in bound B. The site
plausibly stores closed-form coefficients and evaluates them in float64.
We try many natural formulations and see which reproduces ALL observed
site outcomes (200 position / 方向反了 / 未找到).
"""

import sys
from fractions import Fraction

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from sim_float_solve import f_nonneg, f_nonpos, gauss_float
from sim_mechanisms import basis_and_target

from attention_calculator.engine import mn_order, solve_moment


def system(basis, target, bound_sign):
    """Build M (rows=keys, cols=unknowns), t0 (bound-free rhs), t1 (d rhs/d B)."""
    keys = sorted(set(target) | {k for m in basis for k in m})
    M = [[m.get(k, Fraction(0)) for m in basis] for k in keys]
    t0 = [(target[k] if k != "1" else Fraction(0)) for k in keys]
    t1 = [(Fraction(bound_sign) if k == "1" else Fraction(0)) for k in keys]
    return keys, M, t0, t1


def detn(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if n == 3:
        return (
            M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
        )


def cramer_terms(M, t0, t1):
    """Return (P,Q,D): u_j = (P_j + Q_j*B)/D via Cramer."""
    n = len(M)
    D = detn(M)
    if D == 0:
        return None
    P, Q = [], []
    for j in range(n):
        Mj0 = [r[:] for r in M]
        Mj1 = [r[:] for r in M]
        for i in range(n):
            Mj0[i][j] = t0[i]
            Mj1[i][j] = t1[i]
        P.append(detn(Mj0))
        Q.append(detn(Mj1))
    return P, Q, D


def inverse_rows(M):
    """Exact inverse M^{-1} rows for n<=3 via adjugate/det."""
    n = len(M)
    D = detn(M)
    if D == 0:
        return None
    adj = []
    for i in range(n):
        row = []
        for j in range(n):
            # cofactor (j,i): det of minor removing row j col i
            minor = [[M[r2][c2] for c2 in range(n) if c2 != i] for r2 in range(n) if r2 != j]
            sign = -1 if (i + j) % 2 else 1
            row.append(sign * detn(minor) if n > 1 else Fraction(1))
        adj.append(row)
    return [[adj[i][j] / D for j in range(n)] for i in range(n)]  # rows of M^-1


# --- evaluators: each takes (keys, M, t0, t1, Bf) -> list[float] or None ---


def ev_cramer_f64_formula(keys, M, t0, t1, Bf):
    """u_j = (P_j + Q_j*Bf)/D_f : stored (P,Q,D) ints, float eval."""
    ct = cramer_terms(M, t0, t1)
    if ct is None:
        return None
    P, Q, D = ct
    Df = float(D)
    return [(float(p) + float(q) * Bf) / Df for p, q in zip(P, Q, strict=True)]


def ev_affine_f64(keys, M, t0, t1, Bf):
    """u_j = float(P/D) + float(Q/D)*Bf : predivided affine."""
    ct = cramer_terms(M, t0, t1)
    if ct is None:
        return None
    P, Q, D = ct
    return [float(p / D) + float(q / D) * Bf for p, q in zip(P, Q, strict=True)]


def ev_cramer_allfloat(keys, M, t0, t1, Bf):
    """Cramer dets computed in float64 throughout (float moments, float rhs)."""
    n = len(M)
    Mf = [[float(v) for v in row] for row in M]
    tf = [(t0[i] + t1[i] * Fraction(Bf)) for i in range(n)]
    tf = [float(v) for v in tf]

    def detf(A):
        if n == 2:
            return A[0][0] * A[1][1] - A[0][1] * A[1][0]
        if n == 3:
            return (
                A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
            )

    Df = detf(Mf)
    if Df == 0.0:
        return None
    out = []
    for j in range(n):
        Mj = [r[:] for r in Mf]
        for i in range(n):
            Mj[i][j] = tf[i]
        out.append(detf(Mj) / Df)
    return out


def ev_inv_float(keys, M, t0, t1, Bf):
    """u_j = float(i_j0)*t0_0f + ... : stored exact inverse, float dot product."""
    inv = inverse_rows(M)
    if inv is None:
        return None
    tf = [float(t0[i]) + float(t1[i]) * Bf for i in range(len(M))]
    return [sum(float(inv[j][i]) * tf[i] for i in range(len(M))) for j in range(len(M))]


def ev_gauss_f64(keys, M, t0, t1, Bf):
    rhs = [float(t0[i]) + float(t1[i]) * Bf for i in range(len(M))]
    return gauss_float([[float(v) for v in row] for row in M], rhs)


def ev_exact_f64coeff(keys, M, t0, t1, Bf):
    """Exact solve then float64 of coeffs (baseline refuted model)."""
    tgt = {k: t0[i] + t1[i] * Fraction(Bf) for i, k in enumerate(keys)}
    # Bf is float; Fraction(Bf) is the exact dyadic — fine, equals bound intent
    try:
        u = solve_moment(
            [{k: M[i][j] for k in keys for j in [keys.index(k)]} for i in range(0)], tgt
        )  # placeholder
    except Exception:
        return None
    return [float(v) for v in u]


EVALS = {
    "cramer_f": ev_cramer_f64_formula,
    "affine_f": ev_affine_f64,
    "cramer_allf": ev_cramer_allfloat,
    "inv_f": ev_inv_float,
    "gauss_f": ev_gauss_f64,
}


def scan_eval(kind, pw, comp, rs, ev, mid_nonpos):
    """Return (event, (m,n), idx, u). event in nonneg/nonpos/exhaust."""
    bound = Fraction(rs)
    limit, mk, target, B = basis_and_target(kind, pw, comp, bound)
    sign = 1 if comp == ">" else -1
    Bf = float(B)
    for i, (m, n) in enumerate(mn_order(limit)):
        basis = mk(m, n)
        keys, M, t0, t1 = system(basis, target, -sign)
        try:
            u = ev(keys, M, t0, t1, Bf)
        except Exception:
            u = None
        if u is None:
            continue
        if f_nonneg(u):
            return ("nonneg", (m, n), i)
        if mid_nonpos and f_nonpos(u):
            return ("nonpos", (m, n), i)
    return ("exhaust", None, None)


CASES = [
    # name, kind, power, comp, bound, site outcome ('WD','NS','OK@mn')
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
    ("sin >float", "sin_q", "1", ">", "3789648413623927/4503599627370496", "OK(late)"),
]


def predict(evname, kind, pw, comp, rs):
    mid_nonpos = comp == ">"  # only '>' scan aborts on nonpos
    ev = EVALS[evname]
    try:
        res = scan_eval(kind, pw, comp, rs, ev, mid_nonpos)
    except Exception as e:
        return ("EXC", str(e)[:40])
    ev_, mn, _i = res
    if comp == ">":
        return {"nonneg": f"OK{mn}", "nonpos": f"WD@{mn}", "exhaust": "NS"}[ev_]
    else:
        return {"nonneg": f"OK{mn}", "exhaust": "NS", "nonpos": "NS(nonpos-skip)"}[ev_]


if __name__ == "__main__":
    hdr = f"{'case':12s} {'site':10s}" + "".join(f"{n:>14s}" for n in EVALS)
    print(hdr)
    score = {n: 0 for n in EVALS}
    for name, k, p, c, r, want in CASES:
        row = f"{name:12s} {want:10s}"
        for ev in EVALS:
            got = predict(ev, k, p, c, r)
            ok = (
                (want.startswith("OK") and got.startswith("OK"))
                or (want == "WD" and got.startswith("WD"))
                or (want == "NS" and got == "NS")
            )
            if ok:
                score[ev] += 1
            row += f"{got:>14s}"
        print(row)
    print("\nscores:", score)
