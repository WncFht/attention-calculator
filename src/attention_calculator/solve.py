"""Dispatch a (type, power, comparison, rational) request to its kernel family."""

import importlib
import math
from fractions import Fraction

import mpmath as mp

from .engine import NoSolution, WrongDirection
from .kernels import TYPES

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
    GAUSS_F = float(mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi ** 3)))
EULER_F = 0.5772156649015329  # 与 kernels.gamma 的站端字面量一致


def trig_in_domain(kind: str, q: Fraction) -> bool:
    """trig_q 的值域判定：sin_q 是 (0, π)，其余 (0, π/2)，与 check_input 同序。"""
    if q <= 0:
        return False
    with mp.workdps(60):
        return mp.mpf(q.numerator) / q.denominator < (
            mp.pi if kind == "sin_q" else mp.pi / 2)


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
        return math.pi ** f
    if kind == "arctan_q":
        return math.atan(f)
    if kind == "arccot_q":
        return math.atan(1 / f)
    if kind in ("sinh_q", "cosh_q", "tanh_q", "coth_q"):
        return {"sinh_q": math.sinh, "cosh_q": math.cosh, "tanh_q": math.tanh,
                "coth_q": lambda v: 1 / math.tanh(v)}[kind](f)
    if kind in ("sin_q", "cos_q", "tan_q", "cot_q"):
        if not trig_in_domain(kind, q):
            return None
        return {"sin_q": math.sin, "cos_q": math.cos, "tan_q": math.tan,
                "cot_q": lambda v: 1 / math.tan(v)}[kind](f)
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


def prove(kind: str, power: str, comp: str, rational: str) -> dict:
    """Run the proof search; returns the site's /calculate response shape.

    Raises engine.WrongDirection / engine.NoSolution on failure.
    """
    if kind not in TYPES:
        raise ValueError(f"unsupported type {kind!r}")
    module = importlib.import_module(f"attention_calculator.kernels.{FAMILY[kind]}")
    q, r = parse_rational(power), parse_rational(rational)
    # 负界在站端被表层格式校验挡掉（右侧有理数格式无效），不进方向预检
    c = None if ((kind == "zeta3" and comp == ">") or r < 0) else direction_f(kind, q)
    if c is not None and ((r > c) if comp == ">" else (r < c)):
        raise WrongDirection
    try:
        return module.prove(kind, q, comp, r)
    except WrongDirection:
        # '<' 扫描途中的非正解在站端是"未找到解"而非"方向反了"——
        # 站端 '<' 的方向判定只在 bound<c 预检发生，扫描里的非正 P 直接耗尽
        if comp == "<" and c is not None:
            raise NoSolution from None
        raise
    except NoSolution:
        if c is not None:
            raise
        # 无精确预检的型在搜索耗尽后仍按数值真假区分报错：命题为假时报
        # "方向反了"而非"未找到解"（实测 arctan 3 > 5/4 → 方向反了；真命题
        # < 5/4 → 未找到解）。恒等式 ∫f = ±(C−r) 精确成立，真命题不可能搜出
        # 恒≤0 的 P，故仅在 NoSolution 后补判不会误伤已验证路径。
        # 判定精度是 float64：zeta3 '>' 对 float 相等但方向为假的界仍报
        # "未找到解"（float 差为 0 → 放行进入搜索 → 耗尽），故用 float 比较差。
        from .integrand import constant_mpf
        diff = float(constant_mpf(kind, q)) - float(r)
        claim_false = (diff < 0) if comp == ">" else (diff > 0)
        if claim_false:
            raise WrongDirection from None
        raise
