# ruff: noqa: RUF001, RUF002, RUF003  # 中文标点属刻意文体
"""批量验证积分恒等式记录（golden.jsonl 或本方 solver 输出）。

每条 success=true 的记录验证两件事：

1. 恒等式：integrand.reconstruct 重建被积函数 f 与积分域，mpmath 50dps
   数值积分得 actual_integral，与恒等式左端 claimed_lhs（由
   type/comparison/rational 确定）比较 → identity_ok。
2. 符号：自适应网格扫描 f 在区间内部不变号 → sign_ok；
   对含多项式因子 P 的类型另做 Sturm 精确实根计数交叉验证 → sign_exact。

结论分类 verdict：
- valid            恒等式成立且被积函数定号（有效证明）
- indefinite-sign  恒等式成立但被积函数变号（证明无效）
- false-identity   积分真值与声称 LHS 不符（本方 solver bug 的证据）
- unresolved-param 自由参数（ln 族分母幂 s、gamma 核指数 k）反解失败
- harness-error    验证器自身异常

对 success=false 的记录核对真实方向是否与站点报错一致
（error_consistent 字段）。

重建不出的自由参数按恒等式数值反解：枚举候选取 |∫f−lhs| 最小者，
粗精度筛不到候选即判 unresolved-param，命中后全精度复核。

用法::

    .venv/bin/python bench/verify.py bench/data/golden.jsonl [-o out.jsonl]
    .venv/bin/python bench/verify.py recs.jsonl --rigorous   # flint.acb 复核

方法学讨论见 docs/verify-notes.md。
"""

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp
from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from attention_calculator.integrand import constant_mpf, lhs_mpf, reconstruct

TOL_REL = mp.mpf("1e-30")     # 恒等式相对容差
SWEEP_TOL = mp.mpf("1e-15")   # 自由参数反解的粗筛容差
GAMMA_K_MAX = 60              # gamma 核指数搜索上限
LN_S_MAX = 15                 # ln 族分母幂搜索上限
GRID_N = 200                  # 符号扫描初始网格点数
ZERO_REL = mp.mpf("1e-20")    # 符号扫描的"视为零"相对阈值

x_sym = sp.Symbol("x")
k_sym = sp.Symbol("k", integer=True, nonnegative=True)
s_sym = sp.Symbol("s", integer=True, nonnegative=True)
t_sym = sp.Symbol("t")

# P 因子的不定元在各类型下的取值域都化归 t∈(0,1)；
# 三项式族用 a+bt+ct²，其余 a+bt。c_val 被复用的类型不在此列 c。
P_QUADRATIC = {"e", "e_q", "ln_q", "ln_q_square", "sin_q", "cos_q", "tan_q",
               "cot_q", "sinh_q", "cosh_q", "tanh_q", "coth_q"}
P_LINEAR = {"golden", "pi", "pi_n", "catalan", "zeta3",
            "arctan_q", "arccot_q", "artanh_q", "arcoth_q",
            "e_pi", "sin_q_degree", "cos_q_degree", "sin_pi_q", "cos_pi_q",
            "varpi", "gauss"}


def mp_func(f: sp.Expr):
    """sympy 表达式 → mpmath 可调用；求值异常时向内微扰重试。"""
    fn = sp.lambdify(x_sym, f, modules="mpmath")

    def g(t):
        try:
            return fn(t)
        except (ZeroDivisionError, ValueError, OverflowError):
            # 端点奇异在有限精度下可能命中（如 1/(1-x)）：向内微扰
            return fn(t * (1 - mp.mpf("1e-35")))

    return g


def quad_integrand(f: sp.Expr, a, b, dps=50):
    """∫_a^b f dx，返回 (值, quad 自带误差估计)。"""
    mp.dps = dps
    g = mp_func(f)
    return mp.quad(g, [mp.mpf(sp.N(a, 30)), mp.mpf(sp.N(b, 30))], error=True)


def sweep_param(f: sp.Expr, a, b, sym, candidates, lhs, dps=30):
    """枚举自由符号取值，返回使 |∫f−lhs| 最小的 (值, 积分, 偏差)。"""
    best = None
    for v in candidates:
        try:
            val, _ = quad_integrand(f.subs(sym, v), a, b, dps=dps)
        except Exception:
            continue
        dev = abs(val - lhs)
        if best is None or dev < best[2]:
            best = (v, val, dev)
    return best


def sign_scan(f: sp.Expr, a, b, dps=30):
    """自适应网格扫描 f 在 (a,b) 内的符号；返回 +1/-1/0（0=变号）。

    端点用近端点采样代替（奇异端点 tanh-sinh 风格避开）；
    在 |f| 最小处自适应加密三轮，捕捉窄下凹。
    """
    mp.dps = dps
    g = mp_func(f)
    am, bm = mp.mpf(sp.N(a, 30)), mp.mpf(sp.N(b, 30))

    def grid(lo, hi, n):
        return [lo + (hi - lo) * mp.mpf(i) / n for i in range(1, n)]

    xs = [*grid(am, bm, GRID_N), am + mp.mpf("1e-25"), bm - mp.mpf("1e-25")]
    vals = [(t, g(t)) for t in xs]
    for _ in range(3):
        mags = [(t, abs(v)) for t, v in vals if mp.isfinite(v)]
        if not mags:
            return 0
        tmin = min(mags, key=lambda p: p[1])[0]
        width = (bm - am) / GRID_N
        vals += [(t, g(t)) for t in grid(tmin - width, tmin + width, 40)
                 if am < t < bm]
    vmax = max(abs(v) for _, v in vals if mp.isfinite(v))
    pos = neg = False
    for _, v in vals:
        if not mp.isfinite(v):
            continue
        if v > ZERO_REL * vmax:
            pos = True
        elif v < -ZERO_REL * vmax:
            neg = True
    return 0 if pos and neg else (1 if pos else -1)


def poly_factor_sign(kind: str, p: dict):
    """P 因子的精确符号判定：QQ 多项式在 t∈(0,1) 上的 Sturm 实根计数。

    返回 +1（严格正）/-1（严格负）/None（区间内有根或端点为零，
    符号判定失效，由数值扫描兜底）。
    """
    if kind not in P_QUADRATIC | P_LINEAR:
        return None
    a = sp.Rational(p["a_val"])
    b = sp.Rational(p["b_val"])
    c = sp.Rational(p["c_val"]) if kind in P_QUADRATIC else sp.Integer(0)
    P = a + b * t_sym + c * t_sym**2
    flip = 1
    # 端点根（t=0/1 对应积分区间端点）不破坏内部定号：剥离后再数内部根。
    # 注意 (t-1) 在 (0,1) 上为负——每剥一次符号翻转；t 为正不影响。
    while P != 0 and P.subs(t_sym, 0) == 0:
        P = sp.div(P, t_sym, t_sym)[0]
    while P != 0 and P.subs(t_sym, 1) == 0:
        P = sp.div(P, t_sym - 1, t_sym)[0]
        flip = -flip
    if P == 0:
        return None
    if sp.Poly(P, t_sym).count_roots(0, 1):
        return None
    v = P.subs(t_sym, sp.Rational(1, 2))
    return int(sp.sign(v)) * flip if v != 0 else None


def verify_record(rec: dict) -> dict:
    """验证单条记录，返回结论 dict（字段供 bench/report.py 聚合）。"""
    kind, comp = rec["type"], rec["comparison"]
    power, rational = Fraction(rec["power"]), Fraction(rec["rational"])
    out = {"type": kind, "power": rec["power"], "comparison": comp,
           "rational": rec["rational"], "success": rec.get("success", False)}
    mp.dps = 50
    lhs = lhs_mpf(kind, comp, power, rational)
    out["claimed_lhs"] = mp.nstr(lhs, 20)
    if not rec.get("success"):
        out["error"] = rec.get("error", "")
        out["error_consistent"] = direction_consistent(
            kind, comp, power, rational, out["error"])
        return out
    f, a, b = reconstruct(kind, comp, power, rec["parameters"])
    for sym, cands, key in ((k_sym, range(GAMMA_K_MAX + 1), "resolved_k"),
                            (s_sym, range(LN_S_MAX + 1), "resolved_s")):
        if sym not in f.free_symbols:
            continue
        hit = sweep_param(f, a, b, sym, cands, lhs)
        if hit is None or hit[2] > SWEEP_TOL:
            out.update(identity_ok=False, sign_ok=None, verdict="unresolved-param")
            return out
        out[key] = hit[0]
        f = f.subs(sym, hit[0])
    val, err = quad_integrand(f, a, b, dps=50)
    out["actual_integral"] = mp.nstr(val, 30)
    out["quad_error_est"] = mp.nstr(err, 5)
    dev = abs(val - lhs)
    out["abs_deviation"] = mp.nstr(dev, 5)
    out["identity_ok"] = bool(dev <= TOL_REL * max(1, abs(lhs)))
    sgn = sign_scan(f, a, b)
    out["sign_scan"] = {1: "nonneg", -1: "nonpos", 0: "changes"}[sgn]
    out["sign_exact"] = poly_factor_sign(kind, rec["parameters"])
    out["sign_ok"] = sgn != 0  # 扫描为准；精确层记在 sign_exact 供对照
    out["verdict"] = (
        "valid" if out["identity_ok"] and out["sign_ok"]
        else "indefinite-sign" if out["identity_ok"]
        else "false-identity")
    return out


def direction_consistent(kind, comp, power, rational, error):
    """报错 case 核对：站点报错类型与真实方向是否自洽。

    "方向反了" 要求真方向与声称相反；"未找到解" 要求真方向与声称一致
    （式子为真但搜索预算内证不出）；输入校验类报错不涉及方向，记 None。
    """
    if "方向反了" in error:
        want = "opposite"
    elif "未找到" in error:
        want = "same"
    else:
        return None
    c = constant_mpf(kind, power)
    r = mp.mpf(rational.numerator) / rational.denominator
    true_dir = ">" if c > r else "<" if c < r else "="
    ok = (true_dir != comp) if want == "opposite" else (true_dir == comp)
    return {"claimed": comp, "true": true_dir, "consistent": ok}


def main():
    ap = argparse.ArgumentParser(description="verify integral identities")
    ap.add_argument("input", help="jsonl of records")
    ap.add_argument("-o", "--out", help="write per-case results jsonl")
    ap.add_argument("--rigorous", action="store_true",
                    help="对非 valid 记录用 flint.acb.integral 复核")
    args = ap.parse_args()
    with open(args.input) as fh:
        recs = [json.loads(line) for line in fh if line.strip()]
    results, stats = [], {"success_cases": 0, "valid": 0, "false-identity": 0,
                          "indefinite-sign": 0, "unresolved-param": 0,
                          "harness-error": 0, "error_cases": 0,
                          "error_consistent": 0}
    for rec in recs:
        try:
            r = verify_record(rec)
        except Exception as e:
            r = {"type": rec.get("type"), "success": rec.get("success"),
                 "comparison": rec.get("comparison"), "rational": rec.get("rational"),
                 "verdict": "harness-error", "detail": repr(e)}
        results.append(r)
        if rec.get("success"):
            stats["success_cases"] += 1
            stats[r["verdict"]] = stats.get(r["verdict"], 0) + 1
        else:
            stats["error_cases"] += 1
            if (r.get("error_consistent") or {}).get("consistent"):
                stats["error_consistent"] += 1
        comp_s, rat_s = str(r.get("comparison", "?")), str(r.get("rational", "?"))
        print(f"{r['type']:>14} {comp_s:>2} {rat_s:>8} "
              f"-> {r.get('verdict', r.get('error', '?'))} "
              f"{r.get('abs_deviation', '')}")
    if args.rigorous:
        rigorous_recheck(recs, results)
    if args.out:
        with open(args.out, "w") as fh:
            for r in results:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    print("\n== summary ==", json.dumps(stats, ensure_ascii=False))


def rigorous_recheck(recs, results):
    """用 flint.acb.integral（认证区间算术）复核所有非 valid 记录。"""
    try:
        import flint
    except ImportError:
        print("python-flint 不可用，跳过 rigorous 复核")
        return
    for rec, r in zip(recs, results, strict=True):
        if r.get("verdict") not in ("false-identity", "indefinite-sign", "harness-error"):
            continue
        f, a, b = reconstruct(rec["type"], rec["comparison"],
                              Fraction(rec["power"]), rec["parameters"])
        if f.free_symbols - {x_sym}:
            continue  # 自由参数未定的记录复核超出本节范围
        # 用 acb 方法重写被积函数：analytic=True 调用会得到复球输入，
        # 非解析点自动返回非有限球（acb.log/sqrt 自带该行为）。
        fa = sp.lambdify(x_sym, f, modules={
            "sin": lambda z: z.sin(), "cos": lambda z: z.cos(),
            "exp": lambda z: z.exp(), "log": lambda z: z.log(),
            "sqrt": lambda z: z.sqrt(), "sinh": lambda z: z.sinh(),
            "cosh": lambda z: z.cosh(), "tanh": lambda z: z.tanh(),
            "pi": flint.acb.pi()})
        lo = flint.acb(mp.nstr(mp.mpf(sp.N(a, 60)), 50))
        hi = flint.acb(mp.nstr(mp.mpf(sp.N(b, 60)), 50))

        def acb_f(t, analytic, fa=fa):
            return fa(t)

        try:
            # 200bit(~60d) 上下文：区间半径须压过 lhs 的 mpf→arb 解析误差
            flint.ctx.prec = 200
            iv = flint.acb.integral(acb_f, lo, hi, abs_tol=flint.arf(1e-35))
            lhs = lhs_mpf(rec["type"], rec["comparison"],
                          Fraction(rec["power"]), Fraction(rec["rational"]))
            contains = iv.real.contains(flint.arb(mp.nstr(lhs, 60)))
            r["acb_interval"], r["acb_contains_lhs"] = str(iv), bool(contains)
            print(f"acb {rec['type']}: {iv} contains_lhs={contains}")
        except Exception as e:
            r["acb_error"] = repr(e)
            print(f"acb {rec['type']}: ERROR {e!r}")


if __name__ == "__main__":
    main()
