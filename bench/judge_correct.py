# 中文标点属刻意文体
"""mode=exact 正确性判官（docs/2026-09-16-math-correctness-plan.md W2）。

对每条 (type, power, comparison, rational) 独立裁决三份证据：

1. ground truth：常数走 cases.true_value（与 solve.certified_cmp 用的
   integrand.constant_mpf 是两套代码路径），协议同为 80→2400 dps 递增 +
   护栏带；不可实值求值（域外、复值等）记 truth=None。golden 自带的
   success/error 站端标签一概不信。
2. exact 路径：solve.prove(exact=True)。emitted 证明必须过
   exact_check.verify（identity_ok 且 nonneg，零容忍）——prove_exact 已
   内联自检，本层是判官侧复核；BUG:*/FLAG:* 前缀的 verdict 都是要追的
   异常，正常类别见 VERDICT_OK。
3. site 路径对照（默认可 --no-site 关掉）：solve.prove(exact=False) 跑同
   一输入，emitted 参数同样送 exact_check——site 假证明在此现形。两侧
   结局不一致时按已知失真簇归因（见 cluster()），落不进已知簇的标
   "unattributed:*"。

语料：positional 输入为 golden 格式 jsonl（只取 4 个输入字段，可多个文件）；
--adversarial 叠加 cases_correct.generate_cases() 的对抗集：power∈{1,2,3,1/2,0,-1}
网格 × 距离 1e-1…1e-12 的夹逼界、连分数收敛子多深度（最深即等值边界）、
Niven/系数零点的等值命题、域边界两侧、零/负有理界。

用法::

    .venv/bin/python bench/judge_correct.py bench/data/golden.jsonl \
        --adversarial -o bench/out/judge_correct.jsonl
    .venv/bin/python bench/judge_correct.py --adversarial --no-site --limit 200
"""

import argparse
import json
import sys
from collections import Counter
from contextlib import nullcontext
from fractions import Fraction
from pathlib import Path

from mpmath import mp, mpf, workdps

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cases import true_value
from cases_correct import generate_cases

from attention_calculator import solve
from attention_calculator.engine import (
    EqualClaim,
    InternalError,
    NoSolution,
    WrongDirection,
)
from attention_calculator.exact_check import verify_response
from attention_calculator.kernels import EXACT_TYPES
from attention_calculator.solve import certified_cmp

TRUTH_DPS = (80, 240, 800, 2400)

# 站端失真簇归因用的型集合
TRIG_PI_KINDS = {"sin_pi_q", "cos_pi_q", "sin_q_degree", "cos_q_degree"}
BETA_KINDS = {"golden", "varpi", "gauss"}
# 正常 verdict 类别（其余以 BUG:/FLAG: 前缀出现）
VERDICT_OK = {
    "proved",  # 真命题 emitted 且过 exact_check
    "unsolved-true",  # 真命题预算内未证出（诚实）
    "rejected-false",  # 假命题被 WrongDirection 认证拒证
    "unsolved-false",  # 假命题耗尽未命中非正解（诚实但信息量低）
    "equal-claim",  # 等值命题报 EqualClaim
    "rejected",  # 域外/格式/不可判输入被拒
}


def truth_sign(kind: str, q: Fraction, r: Fraction) -> int | None:
    """sign(C - r) 的独立认证：+1/-1/0；None = 非实值可判。

    与 solve.certified_cmp 同协议（递增 dps + 2^30-ulp 护栏带），常数实现
    换用 cases.true_value 作交叉验证路径。
    """
    for dps in TRUTH_DPS:
        with workdps(dps):
            try:
                c = true_value(kind, q)
                if not isinstance(c, mp.mpf) or not mp.isfinite(c):
                    return None
                diff = c - mpf(r.numerator) / r.denominator
                guard = mpf(2) ** (30 - dps) * max(1, abs(c))
            except (TypeError, ValueError, OverflowError, ZeroDivisionError):
                return None
            if diff > guard:
                return 1
            if diff < -guard:
                return -1
            if dps == TRUTH_DPS[-1]:
                return 0
    return None


def claim_true(comp: str, truth: int) -> bool:
    """truth = sign(C - r) 下命题 comp 是否为真。"""
    return (truth > 0) == (comp == ">")


def in_domain(kind: str, q: Fraction) -> bool:
    """核内 check_input 域规则的镜像（只为期望裁决分类，不做拦截）。

    系数型（pi/e/gamma/…）与 e_q/e_pi 无显式域校验，返回 True——
    这些型对负/零 power 的处置本身就是被评测对象。pi_n 的域（bb6e35a）：
    分子 >= 1 且分母 <= 64——β/η 生成器让任意正分子可证，分母上界
    防 bound**pd 整数爆炸。
    """
    if kind == "pi_n":
        return q.numerator >= 1 and q.denominator <= 64
    if kind in ("ln_q", "ln_q_square", "ln_q_cube"):
        return q > 1
    if kind == "arcsin_q":
        return 0 < q < 1
    if kind in ("gaussint_q", "dawson_q", "erfiint_q"):
        return q != 0
    # 460c5b1 域规则：e_q/cosh/sinh/tanh 的 q=0 退化拒（e^0=1 等有理点），
    # coth 要求 q>0；q<0 经奇偶归约在域内
    if kind == "e_q" or kind in ("sinh_q", "cosh_q", "tanh_q"):
        return q != 0
    if kind == "coth_q":
        return q > 0
    # 441a7d7 新型域规则：系数型要求 q≠0（q·C 退化有理比较），
    # li2_q 参数域 q<1 且 q≠0，psi1_q 要求 q>0
    if kind in ("pi3", "pi3_u", "pi3_a", "gamma14", "gamma34", "gamma12"):
        return q != 0
    if kind == "li2_q":
        return q != 0 and q < 1
    if kind == "psi1_q":
        return q > 0
    # si_q/cin_q 唯一域检是 q≠0（q<0 走奇偶归约）；zeta/beta 系数型全域
    if kind in ("si_q", "cin_q"):
        return q != 0
    if kind == "artanh_q":
        return 0 < q < 1 and q.denominator != 1
    if kind == "arcoth_q":
        return q > 1
    if kind in ("sin_pi_q", "cos_pi_q"):
        return 0 < q < Fraction(1, 2) and q.denominator != 1
    if kind in ("sin_q_degree", "cos_q_degree"):
        return 0 < q < 90
    with workdps(60):
        x = mpf(q.numerator) / q.denominator
        if kind == "sin_q":
            return 0 < x < mp.pi
        if kind in ("cos_q", "tan_q", "cot_q"):
            return 0 < x < mp.pi / 2
    return True


def run_outcome(kind: str, q: Fraction, comp: str, r: Fraction, exact: bool) -> dict:
    """一次 solve.prove 归约为 outcome 词表；emitted 时附 exact_check 复核。"""
    try:
        resp = solve.prove(kind, str(q), comp, str(r), exact=exact)
    except WrongDirection:
        return {"outcome": "wrong-direction"}
    except NoSolution:
        return {"outcome": "no-solution"}
    except EqualClaim:
        return {"outcome": "equal"}
    except InternalError as e:
        return {"outcome": "internal-error", "detail": str(e)}
    except ZeroDivisionError as e:
        # 判官送入的输入都已过 Fraction 解析，prove 内部除零即内核 bug
        # （如 e_q/cosh_q power=0 的 Fraction(1,0)），不是输入拒绝
        return {"outcome": "crash", "detail": str(e)}
    except ValueError as e:
        return {"outcome": "rejected-input", "detail": str(e)}
    except Exception as e:
        return {"outcome": "crash", "detail": repr(e)}
    out = {"outcome": "emitted", "parameters": resp.get("parameters")}
    if resp.get("prover"):
        out["prover"] = resp["prover"]  # pade/composite 响应无 parameters，证书即证明
    try:
        # verify_response 内部分派 classic/pade/composite 并比对证书所证命题
        chk = verify_response(kind, q, comp, r, resp)
        out["identity_ok"], out["nonneg"] = chk["identity_ok"], chk["nonneg"]
    except Exception as e:
        out["verify_error"] = repr(e)
    return out


def verdict(case: dict, truth: int | None, res: dict) -> str:
    """exact 结局 × ground truth → verdict；BUG:/FLAG: 前缀需人工追。"""
    comp = case["comparison"]
    oc = res["outcome"]
    if oc == "emitted":
        if "verify_error" in res:
            return "BUG:verify-crash"
        if not (res["identity_ok"] and res["nonneg"]):
            # prove_exact 自检失败应抛 InternalError，能走到这是管线漏洞
            return "BUG:false-proof"
        if truth is None:
            return "FLAG:emitted-truth-unknown"
        if truth == 0:
            return "BUG:proved-equal"
        return "proved" if claim_true(comp, truth) else "BUG:proved-false"
    if oc == "wrong-direction":
        if truth is None:
            return "rejected"
        if truth == 0:
            return "FLAG:wd-on-equal"  # certified_cmp 应先报等值
        return "rejected-false" if not claim_true(comp, truth) else "BUG:wd-on-true"
    if oc == "no-solution":
        if truth is None:
            return "rejected"
        if truth == 0:
            return "FLAG:ns-on-equal"
        return "unsolved-true" if claim_true(comp, truth) else "unsolved-false"
    if oc == "equal":
        if truth == 0:
            return "equal-claim"
        if truth is None:
            return "rejected"
        return "BUG:equal-on-unequal"
    if oc == "rejected-input":
        # 负有理界被表层格式校验拦（solve.py 刻意的站端语法对齐），
        # 与域外输入同属设计行为；其余域内真命题被挡才是覆盖缺口
        try:
            if Fraction(case["rational"]) < 0:
                return "rejected"
        except (ValueError, ZeroDivisionError):
            return "rejected"
        if (
            truth is not None
            and truth != 0
            and claim_true(comp, truth)
            and in_domain(case["type"], Fraction(case["power"]))
        ):
            return "FLAG:rejected-true-claim"
        return "rejected"
    if oc == "internal-error":
        return "BUG:internal-error"
    return "BUG:crash"


def cluster(case: dict, truth: int | None, exact_res: dict, site_res: dict) -> str:
    """exact/site 结局分歧的失真簇归因；只在两侧 outcome 不同时调用。"""
    kind, comp = case["type"], case["comparison"]
    q = Fraction(case["power"])
    so, eo = site_res["outcome"], exact_res["outcome"]
    if kind in EXACT_TYPES and so in ("rejected-input", "crash"):
        return "exact-only-type"  # exact 专属型，站端本就 400
    if eo == "rejected-input" and not in_domain(kind, q):
        return "domain-policy"  # exact 域规则拒（如 pi_n 分母 >64），site 可证属预期
    if so == "emitted" and site_res.get("identity_ok") is False:
        # 站端发了假证明，exact 拒绝或改判 —— 按失真簇归因
        if kind in TRIG_PI_KINDS:
            return "trig-pi-bias18"
        if kind in BETA_KINDS and comp == "<":
            return "beta-transposed-lt"
        if q != 1:
            return "power-scaling"
        return "unattributed:false-proof"
    if so == "emitted" and site_res.get("nonneg") is False:
        return "site-sign-indefinite"  # 恒等式成立但被积函数变号 —— 证明无效
    if so == "emitted" and truth == 0:
        # 等值命题上站端发了"严格不等式"证明（恒零被积函数 vacuous 过关），
        # exact 正确报二者相等
        return "site-equality-slippage"
    if so == "internal-error":
        if kind == "ln_q_square" and q in (Fraction(5), Fraction(7)):
            return "ln2-singular-500"
        if q == 0:
            return f"{kind}-power0-crash"
        return "unattributed:site-500"
    if so == "wrong-direction":
        # 站端 float64 预检判反：真命题被误拒（exact 证出/诚实报未找到都算此簇）
        if eo == "emitted" or (
            eo == "no-solution" and truth is not None and truth != 0 and claim_true(comp, truth)
        ):
            return "float64-direction"
        if eo == "no-solution":
            return "direction-vs-exhaustion"  # 双侧都拒证假命题，措辞不同
    if so == "no-solution" and eo == "wrong-direction":
        return "certified-direction"  # exact 认证判反 vs 站端仅搜索耗尽
    if so == "crash" and eo != "crash":
        return "site-crash"  # 站端崩溃的输入 exact 给出干净答案
    if eo == "emitted" and so == "no-solution":
        return "exact-coverage-gain"  # exact 证出站端预算外真命题
    if (
        so == "emitted"
        and site_res.get("identity_ok")
        and site_res.get("nonneg")
        and eo != "emitted"
    ):
        return "exact-coverage-loss"  # 站端证明为真而 exact 没证出 —— 复查
    if "rejected-input" in (so, eo):
        return "domain-policy"  # 输入校验口径差异（负界等）
    return "unattributed:other"


def judge_case(case: dict, site: bool = True) -> dict:
    """裁决一条 (type, power, comparison, rational)，返回记录 dict。"""
    kind, comp = case["type"], case["comparison"]
    rec = {k: case[k] for k in ("type", "power", "comparison", "rational")}
    try:
        q, r = Fraction(case["power"]), Fraction(case["rational"])
    except (ValueError, ZeroDivisionError) as e:
        rec.update(
            truth=None, exact={"outcome": "rejected-input", "detail": str(e)}, verdict="rejected"
        )
        return rec
    truth = truth_sign(kind, q, r)
    try:
        cert = certified_cmp(kind, q, r)
    except Exception:
        cert = "crash"
    ex = run_outcome(kind, q, comp, r, exact=True)
    rec["truth"] = truth
    if cert != truth:  # 两套常数实现同协议交叉核验；不一致本身就是 bug 信号
        rec["truth_mismatch"] = cert
    rec["exact"] = {k: v for k, v in ex.items() if k != "parameters"}
    if ex["outcome"] == "emitted" and ex["parameters"]:
        rec["exact"]["mn"] = [ex["parameters"].get("m"), ex["parameters"].get("n")]
    rec["verdict"] = verdict(case, truth, ex)
    if cert != truth and rec["verdict"] in VERDICT_OK:
        rec["verdict"] = "FLAG:truth-mismatch"
    if site:
        st = run_outcome(kind, q, comp, r, exact=False)
        rec["site"] = {k: v for k, v in st.items() if k != "parameters"}
        if st["outcome"] != ex["outcome"]:
            rec["divergence"] = cluster(case, truth, ex, st)
        elif st["outcome"] == "emitted" and st["parameters"] != ex["parameters"]:
            rec["divergence"] = "params-differ"
    return rec


# --------------------------------------------------------------------- main


def fmt_line(rec: dict) -> str:
    """单行可读输出。"""
    t = {1: ">", -1: "<", 0: "=", None: "?"}[rec["truth"]]
    ex = rec["exact"]
    line = (
        f"{rec['type']:>14} {rec['power']:>10} {rec['comparison']} {rec['rational']:>24}"
        f" | C{t}r | exact={ex['outcome']}"
    )
    if ex["outcome"] == "emitted":
        line += f"(id={ex.get('identity_ok')},nn={ex.get('nonneg')})"
    site = rec.get("site")
    if site:
        line += f" site={site['outcome']}"
        if site["outcome"] == "emitted":
            line += f"(id={site.get('identity_ok')},nn={site.get('nonneg')})"
    line += f" | {rec['verdict']}"
    if "divergence" in rec:
        line += f" ⇄{rec['divergence']}"
    return line


def summarize(recs: list[dict]) -> dict:
    """聚合统计：verdict 分布、每型覆盖率、覆盖缺口与分歧簇。"""
    stats = {
        "cases": len(recs),
        "verdicts": dict(Counter(r["verdict"] for r in recs)),
        "truth": dict(Counter(str(r["truth"]) for r in recs)),
        "divergences": dict(Counter(r["divergence"] for r in recs if "divergence" in r)),
    }
    cov, misses = {}, []
    for r in recs:
        if not (r["truth"] in (1, -1) and claim_true(r["comparison"], r["truth"])):
            continue
        row = cov.setdefault(r["type"], [0, 0])
        row[1] += 1
        if r["verdict"] == "proved":
            row[0] += 1
        else:
            misses.append({k: r[k] for k in ("type", "power", "comparison", "rational", "verdict")})
    stats["coverage_by_type"] = {
        k: f"{p}/{t}" for k, (p, t) in sorted(cov.items(), key=lambda kv: kv[1][0] / kv[1][1])
    }
    stats["coverage_misses"] = misses  # certified-true 但未证出/被误拒
    return stats


def main():
    ap = argparse.ArgumentParser(description="judge mode=exact correctness")
    ap.add_argument("inputs", nargs="*", help="golden 格式 jsonl（只取输入字段，可多个）")
    ap.add_argument("--adversarial", action="store_true", help="叠加对抗语料")
    ap.add_argument("--no-site", action="store_true", help="跳过 site 模式对照")
    ap.add_argument("-o", "--out", help="逐条结果 jsonl 输出路径")
    ap.add_argument("--limit", type=int, help="只跑前 N 条")
    ap.add_argument("--kind", help="只跑指定 type")
    args = ap.parse_args()

    cases = []
    for path in args.inputs:
        with open(path) as fh:
            cases += [
                {k: rec[k] for k in ("type", "power", "comparison", "rational")}
                for rec in map(json.loads, fh)
                if rec.get("type")
            ]
    if args.adversarial:
        cases += generate_cases()
    if args.kind:
        cases = [c for c in cases if c["type"] == args.kind]
    if args.limit:
        cases = cases[: args.limit]
    cases = list(
        {(c["type"], c["power"], c["comparison"], c["rational"]): c for c in cases}.values()
    )

    recs = []
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    # 边跑边写：中途崩溃不丢已完成条目
    with open(args.out, "w") if args.out else nullcontext() as out_fh:
        for case in cases:
            rec = judge_case(case, site=not args.no_site)
            recs.append(rec)
            print(fmt_line(rec))
            if out_fh:
                out_fh.write(json.dumps(rec, ensure_ascii=False, default=str) + "\n")
                out_fh.flush()
    stats = summarize(recs)
    bad = [
        r
        for r in recs
        if not r["verdict"].startswith(("proved", "unsolved", "rejected", "equal"))
        or str(r.get("divergence", "")).startswith("unattributed")
    ]
    if bad:
        stats["flagged"] = [
            {k: r[k] for k in ("type", "power", "comparison", "rational", "verdict") if k in r}
            | ({"divergence": r["divergence"]} if "divergence" in r else {})
            | ({"truth_mismatch": r["truth_mismatch"]} if "truth_mismatch" in r else {})
            for r in bad
        ]
    print("\n== summary ==", json.dumps(stats, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
