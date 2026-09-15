"""Simulate 'stored formula' sign-check mechanisms.

Per candidate (m,n) the solution u(B) is affine in the bound B:
u_j(B) = (P_j + Q_j*B)/D  (Cramer: shared denominator D = det)
       = A_j + B_j*B      (predivided)

Float-eval variants:
  W1: (float(P) + float(Q)*Bf) / float(D)
  W2: float(A) + float(B)*Bf  where A=P/D, B=Q/D exact-divided then float
  W3: float-eval of Horner/other orderings — try (P + Q*B_f)/D with the
      division done once at the end vs per-term
Also GJ-float (from sim_float_solve) for reference.
"""

import sys
from fractions import Fraction

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from sim_float_solve import f_nonneg, f_nonpos
from sim_mechanisms import basis_and_target

from attention_calculator.engine import mn_order, solve_moment


def cramer_forms(basis, target):
    """Return per-unknown (P, Q, D) exact rationals for u(B) = (P + Q*B)/D.

    rhs(B) = rhs0 + B*rhs1 where rhs1 is the bound coefficient vector.
    Solve M u = rhs0 and M u = rhs1 over QQ, combine over common denom.
    """
    keys = sorted(set(target) | {k for m in basis for k in m})
    # The bound enters only through the "1"-component: target["1"] = -sign*B.
    # Other components are bound-independent.
    rhs0 = {k: (Fraction(0) if k == "1" else v) for k, v in target.items()}
    # coefficient of B in each rhs component
    rhs1 = {k: (target["1"] if k == "1" else Fraction(0)) for k in keys}
    # target["1"] = -sign * B  →  derivative wrt B is -sign
    # recover -sign by evaluating: rhs(B) - rhs(0) = B * d; d = -sign on "1"
    # We know target["1"] = -sign*B so d = target["1"]/B -- but we don't have
    # B here; instead solve for the coefficient: since target["1"] = -sign*B,
    # the B-derivative vector is (0,...,-sign,...,0). Get it by solving twice:
    # simpler: return both solves and let caller combine.
    u0 = solve_moment(basis, rhs0)
    # unit-bound derivative: solve with "1"-component = -sign (i.e. per unit B)
    return keys, u0, rhs1


def forms_2x2(basis, target, bound_sign):
    """Direct Cramer for the common 2-unknown case + generic for 3."""
    keys = sorted(set(target) | {k for m in basis for k in m})
    n = len(keys)
    M = [[m.get(k, Fraction(0)) for m in basis] for k in keys]  # rows=keys
    # det
    if n == 2:
        det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    elif n == 3:
        det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
               - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
               + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    if det == 0:
        return None
    # rhs(B) = t0 + B*t1 ; t1 = -sign on "1" (bound enters linearly)
    t0 = [target.get(k, Fraction(0)) if k != "1" else Fraction(0) for k in keys]
    # the B-coefficient: target["1"] = -sign * B -> d = -sign
    t1 = [Fraction(0)] * n
    if "1" in keys:
        t1[keys.index("1")] = Fraction(bound_sign)  # bound_sign = -sign
    # adjugate solve: u_j = det(M with col j replaced by rhs)/det
    P, Q = [], []
    for j in range(n):
        Mj0 = [row[:] for row in M]
        Mj1 = [row[:] for row in M]
        for i in range(n):
            Mj0[i][j] = t0[i]
            Mj1[i][j] = t1[i]
        P.append(detn(Mj0))
        Q.append(detn(Mj1))
    return P, Q, det


def detn(M):
    n = len(M)
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if n == 3:
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    raise ValueError


def scan_formula(kind, pw, comp, rs, evalkind, limit=None):
    """evalkind: 'cramer' -> (P+Q*Bf)/Df ; 'affine' -> float(P/D)+float(Q/D)*Bf."""
    bound = Fraction(rs)
    limit_, mk, target, B = basis_and_target(kind, pw, comp, bound)
    if limit:
        limit_ = limit
    sign = 1 if comp == ">" else -1
    Bf = float(B)
    bound_sign = -sign  # target["1"] = -sign*B
    for i, (m, n) in enumerate(mn_order(limit_)):
        basis = mk(m, n)
        try:
            f = forms_2x2(basis, target, bound_sign)
        except Exception:
            f = None
        if f is None:
            continue
        P, Q, D = f
        u = []
        for p, q in zip(P, Q, strict=True):
            if evalkind == "cramer":
                u.append((float(p) + float(q) * Bf) / float(D))
            elif evalkind == "affine":
                u.append(float(p / D) + float(q / D) * Bf)
            elif evalkind == "cramer_seq":
                # (P + Q*B)/D evaluated strictly left-assoc in float
                u.append((float(p) + float(q) * Bf) / float(D))
        if f_nonneg(u):
            return ("nonneg", (m, n), i, u)
        if f_nonpos(u):
            return ("nonpos", (m, n), i, u)
    return ("exhaust", None, None, None)


CASES = [
    ("pi  >tight   方向反了", "pi", "1", ">", "6134899525417045/1952799169684491"),
    ("pi  <float   未找到", "pi", "1", "<", "884279719003555/281474976710656"),
    ("pi  >float   200@12,20", "pi", "1", ">", "884279719003555/281474976710656"),
    ("pi  <tight   未找到", "pi", "1", "<", "5706674932067741/1816491048114374"),
    ("e   >tight   方向反了", "e", "1", ">", "2124008553358849/781379079653017"),
    ("e   <float   未找到", "e", "1", "<", "6121026514868073/2251799813685248"),
    ("e   >float   200@7,7", "e", "1", ">", "6121026514868073/2251799813685248"),
    ("e   <tight   200@11,12", "e", "1", "<", "9581113603440437/3524694718233772"),
    ("e_q >tight   方向反了", "e_q", "1", ">", "2124008553358849/781379079653017"),
    ("e_q3>tight   未找到", "e_q", "3", ">", "9060589063489840/451100167157089"),
    ("sin >tight   方向反了", "sin_q", "1", ">", "8199121568317184/9743795943467975"),
    ("cos >tight   未找到", "cos_q", "1", ">", "293104830616638/542483027433479"),
    ("ln  >tight   方向反了", "ln_q", "2", ">", "1554903831458736/2243252046704767"),
    ("lnsq>tight   方向反了", "ln_q_square", "2", ">", "4646620020445232/9671330777074349"),
    ("zeta>tight   未找到", "zeta3", "1", ">", "5517909911604193/4590389936699688"),
    ("atan>tight   未找到", "arctan_q", "3", ">", "125014145208443/100087721339793"),
]

if __name__ == "__main__":
    for name, k, p, c, r in CASES:
        out = {}
        for ev in ("cramer", "affine"):
            try:
                out[ev] = scan_formula(k, p, c, r, ev)
            except Exception as e:
                out[ev] = ("EXC", str(e))
        print(f"{name}:  cramer={out['cramer'][:3]}  affine={out['affine'][:3]}")
