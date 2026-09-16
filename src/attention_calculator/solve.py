"""Dispatch a (type, power, comparison, rational) request to its kernel family."""

import importlib
import math
from fractions import Fraction

import mpmath as mp

from .engine import EXPONENT_LIMIT, NoSolution, WrongDirection
from .kernels import EXACT_TYPES, TYPES

# type -> module under attention_calculator.kernels
FAMILY = {
    "pi": "quadlog",
    "pi_n": "quadlog",
    "catalan": "quadlog",
    "zeta3": "quadlog",
    "arctan_q": "quadlog",
    "arccot_q": "quadlog",
    "e": "exp_family",
    "e_q": "exp_family",
    "e_pi": "exp_family",
    "sinh_q": "hyperbolic",
    "cosh_q": "hyperbolic",
    "tanh_q": "hyperbolic",
    "coth_q": "hyperbolic",
    "sin_q": "trig_q",
    "cos_q": "trig_q",
    "tan_q": "trig_q",
    "cot_q": "trig_q",
    "sin_pi_q": "trig_pi",
    "cos_pi_q": "trig_pi",
    "sin_q_degree": "trig_pi",
    "cos_q_degree": "trig_pi",
    "ln_q": "log_family",
    "ln_q_square": "log_family",
    "artanh_q": "log_family",
    "arcoth_q": "log_family",
    "golden": "beta_family",
    "varpi": "beta_family",
    "gauss": "beta_family",
    "gamma": "gamma",
    # EXACT_TYPES — no site counterpart; kernels are exact-path only
    "zeta5": "zeta_odd",
    "zeta7": "zeta_odd",
    "zeta9": "zeta_odd",
    "zeta11": "zeta_odd",
    "beta4": "beta_even",
    "beta6": "beta_even",
    "beta8": "beta_even",
    "beta10": "beta_even",
    "ln_q_cube": "ln_pow",
    "arcsin_q": "arcsin",
    "arsinh_q": "invhyp",
    "gaussint_q": "gauss_erf",
    "dawson_q": "gauss_erf",
    "erfiint_q": "gauss_erf",
    "pi_sqrt2": "pi_sqrt2",
    "pi3": "dixon",
    "pi3_u": "dixon",
    "pi3_a": "dixon",
    "li2_q": "li2",
    "psi1_q": "trigamma",
    "si_q": "sicin",
    "cin_q": "sicin",
    "gamma14": "gamma_special",
    "gamma34": "gamma_special",
    "gamma12": "gamma_special",
}

# 站端在进核前用 float64 求值命题常数 c，再把有理界与 c 做 *精确* 比较
# （Fraction vs float64 的 Python 语义，bound 侧不舍入）：'>' 命题 bound>c、
# '<' 命题 bound<c 即 404 "方向反了"，根本不进搜索。这解释了两类保真分歧：
# float 等值但严格偏上的真命题界被站端误判方向（pi/e/sin/ln/e_q 的 ~1e-32 紧界），
# 以及 '<' 方向下界恰等于 c 时站端越过预检、在搜索耗尽后才报"未找到解"。
# 各型 c 取站端字面表达式的 float64 结果，不一定是正确舍入：
# ln_q_square 是 log(q)**2 的二次舍入（低 1 ulp），e_pi 是正确舍入的 fl(e^π)
# （站端不是 math.exp(math.pi)，那个低 1 ulp）。
with mp.workdps(60):
    CATALAN_F = float(mp.catalan)
    ZETA3_F = float(mp.zeta(3))
    E_PI_F = float(mp.exp(mp.pi))
    VARPI_F = float(mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi)))
    GAUSS_F = float(mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi**3)))
EULER_F = 0.5772156649015329  # 站端 float64 字面量，kernels.gamma / decompose 共用


def trig_in_domain(kind: str, q: Fraction) -> bool:
    """trig_q 的值域判定：sin_q 是 (0, π)，其余 (0, π/2)，与 check_input 同序。"""
    if q <= 0:
        return False
    with mp.workdps(60):
        return mp.mpf(q.numerator) / q.denominator < (mp.pi if kind == "sin_q" else mp.pi / 2)


def direction_f(kind: str, q: Fraction):
    """站端方向预检的 float64 常数 c；None = 该请求没有方向预检。

    域校验先于方向的型（trig_q 值域外、ln_q/ln_q_square q<=1）返回 None，
    交给核内 check_input 报错；ln_q_square q∈{5,7} 由核内 float-diff 特判
    （真→500、假→404），同样不做此预检。artanh_q/arcoth_q/trig_pi 四型的
    方向判定路径未经探测区分，保持核内现状。
    """
    f = float(q)
    if kind == "pi":
        return f * math.pi
    if kind == "e":
        return f * math.e
    if kind == "catalan":
        return f * CATALAN_F
    if kind == "gamma":
        return f * EULER_F
    if kind == "golden":
        return f * (1 + math.sqrt(5)) / 2
    if kind == "varpi":
        return f * VARPI_F
    if kind == "gauss":
        return f * GAUSS_F
    if kind == "e_pi":
        return f * E_PI_F
    if kind == "e_q":
        return math.exp(f)
    if kind == "pi_n":
        return math.pi**f
    if kind == "arctan_q":
        return math.atan(f)
    if kind == "arccot_q":
        return math.atan(1 / f)
    if kind in ("sinh_q", "cosh_q", "tanh_q", "coth_q"):
        return {
            "sinh_q": math.sinh,
            "cosh_q": math.cosh,
            "tanh_q": math.tanh,
            "coth_q": lambda v: 1 / math.tanh(v),
        }[kind](f)
    if kind in ("sin_q", "cos_q", "tan_q", "cot_q"):
        if not trig_in_domain(kind, q):
            return None
        return {
            "sin_q": math.sin,
            "cos_q": math.cos,
            "tan_q": math.tan,
            "cot_q": lambda v: 1 / math.tan(v),
        }[kind](f)
    if kind == "ln_q":
        return math.log(f) if q > 1 else None
    if kind == "ln_q_square":
        return math.log(f) ** 2 if q > 1 and q not in (5, 7) else None
    if kind == "zeta3":
        # 仅 '<'：'>' 的站端判定是 float(bound)>c（阈值在 c 上方 1 ulp 外），
        # 由下方 NoSolution 兜底覆盖，不走精确预检
        return f * ZETA3_F
    return None


def parse_rational(text: str) -> Fraction:
    """Parse '3', '22/7' into a Fraction."""
    return Fraction(text.strip())


def certified_cmp(kind: str, q: Fraction, r: Fraction) -> int | None:
    """Certified sign of C - r for mode=exact direction decisions.

    Evaluates the target constant at escalating mpmath precision and
    accepts the sign only when |C - r| clears a wide guard band around
    the evaluation error (2^30 ulp at the working precision). Returns
    +1/-1/0, or None when the constant is not real-evaluable at this q
    (out-of-domain input — the kernel's own validation then decides).
    """
    from .integrand import constant_mpf

    for dps in (80, 240, 800, 2400):
        with mp.workdps(dps):
            try:
                c = constant_mpf(kind, q)
                if not mp.isfinite(c):
                    return None
                diff = c - mp.mpf(r.numerator) / r.denominator
                guard = mp.mpf(2) ** (30 - dps) * max(1, abs(c))
                if diff > guard:
                    return 1
                if diff < -guard:
                    return -1
            except (TypeError, ValueError, OverflowError, ZeroDivisionError):
                return None
            if abs(diff) <= guard and dps == 2400:
                # |C - r| below the smallest guard band: treat as equality
                # (only reachable for rational C, e.g. Niven points)
                return 0
    raise NoSolution


def failure_text(exc: Exception, kind: str, comp: str) -> str:
    """站端 404 文案：WrongDirection -> 方向反了；NoSolution -> 预算内未找到解。"""
    if isinstance(exc, WrongDirection):
        return "要证明的式子不等号方向反了"
    limit = EXPONENT_LIMIT.get(kind, 10)
    return f"在指数不超过{limit}的范围内未找到{comp}方向的解"


def prove(kind: str, power: str, comp: str, rational: str, exact: bool = False) -> dict:
    """Run the proof search; returns the site's /calculate response shape.

    Raises engine.WrongDirection / engine.NoSolution on failure.

    ``exact`` selects the math-correctness path (plan doc W1): certified
    direction comparison instead of float64, kernels solve the true
    system (no reproduced site bugs), and every emitted proof is
    re-verified by exact_check before returning.
    """
    if kind not in TYPES:
        if exact and kind in EXACT_TYPES:
            module = importlib.import_module(f"attention_calculator.kernels.{FAMILY[kind]}")
            q, r = parse_rational(power), parse_rational(rational)
            return prove_exact(module, kind, q, comp, r)
        raise ValueError(f"unsupported type {kind!r}")
    module = importlib.import_module(f"attention_calculator.kernels.{FAMILY[kind]}")
    q, r = parse_rational(power), parse_rational(rational)
    if exact:
        return prove_exact(module, kind, q, comp, r)
    # 负界在站端被表层格式校验挡掉（右侧有理数格式无效），不进方向预检
    c = None if ((kind == "zeta3" and comp == ">") or r < 0) else direction_f(kind, q)
    if c is not None and ((r > c) if comp == ">" else (r < c)):
        raise WrongDirection
    try:
        return module.prove(kind, q, comp, r)
    except WrongDirection:
        if comp == "<":
            # '<' 扫描途中的非正解在站端是"未找到解"而非"方向反了"——
            # 站端 '<' 的方向判定只在 bound<c 预检发生，扫描里的非正 P 直接耗尽；
            # 未映射的型核内自带方向判定（trig_pi defer 等），WD 原样上报
            if c is not None:
                raise NoSolution from None
            raise
        # '>' 命中的非正 P 在站端同样只是"搜不到"：bound==c 的 cos/golden/sin_q
        # 与 bound<C 的 sin_pi_q 实测均报"未找到>方向的解"而非"方向反了"——
        # 不直接判反，落入下方统一的 float64 真假兜底
    except NoSolution:
        # 已映射 '<' 的"未找到"直接上报；'>' 全型与未映射 '<' 进下方兜底
        if comp == "<" and c is not None:
            raise
    # 搜索耗尽（含 '>' 途中非正）后按 float64 真假区分报错：命题为假报
    # "方向反了"而非"未找到解"（实测 arctan 3 > 5/4 → 方向反了；真命题
    # < 5/4 → 未找到解）。恒等式 ∫f = ±(C−r) 精确成立，真命题不可能搜出
    # 恒≤0 的 P，故仅在耗尽后补判不会误伤已验证路径。
    # 判定精度是 float64：zeta3 '>' 对 float 相等但方向为假的界仍报
    # "未找到解"（float 差为 0 → 放行进入搜索 → 耗尽），故用 float 比较差。
    from .integrand import constant_mpf

    diff = float(constant_mpf(kind, q)) - float(r)
    claim_false = (diff < 0) if comp == ">" else (diff > 0)
    if claim_false:
        raise WrongDirection from None
    raise NoSolution from None


def prove_exact(module, kind: str, q: Fraction, comp: str, r: Fraction) -> dict:
    """mode=exact path: certified direction, true-system solve, self-check.

    Direction is decided by certified_cmp before the search; a mid-scan
    WrongDirection under true moments certifies the opposite inequality,
    so it propagates honestly (no '<'-to-NoSolution remap). The emitted
    parameters are re-verified by exact_check — a failure is our bug,
    reported as InternalError rather than a wrong proof.
    """
    from .certificate import build as build_cert
    from .engine import EqualClaim, InternalError

    sign = certified_cmp(kind, q, r)
    if sign == 0:
        raise EqualClaim("二者相等")
    if sign is not None and (sign < 0) == (comp == ">"):
        raise WrongDirection
    try:
        resp = module.prove(kind, q, comp, r, exact=True)
    except NoSolution:
        # W4 Padé 第二证法：ln_q/arctan_q 的 (m,n) 预算耗尽后用插值
        # Padé 界兜底；证书即证明，无 (m,n,P) 参数可渲染
        if kind in ("ln_q", "arctan_q"):
            from . import pade

            cert = pade.prove(kind, q, comp, r)
            if cert is not None:
                return {"type": kind, "prover": "pade", "certificate": cert}
        # W7 AGM 第二证法：gauss 直出区间包络证书；varpi 走 pi 子证 +
        # AGM 的 composite DAG（cert["prover"]=="composite"）
        if kind in ("varpi", "gauss"):
            from . import agm

            cert = agm.prove(kind, q, comp, r)
            if cert is not None:
                return {"type": kind, "prover": cert["prover"], "certificate": cert}
        raise
    if resp.get("prover") == "composite":
        # W6 Γ-型复合证明：证书是子证明 DAG，children 已在构造时逐个验证
        return resp
    # certificate.build internally runs exact_check.verify — its embedded
    # check block IS the emit self-check, so a failure stays InternalError
    cert = build_cert(kind, q, comp, r, resp["parameters"])
    if not (cert["check"]["identity_ok"] and cert["check"]["nonneg"]):
        raise InternalError("emitted proof failed exact self-check")
    resp["certificate"] = cert
    return resp
