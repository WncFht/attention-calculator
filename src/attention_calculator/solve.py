"""Dispatch a (type, power, comparison, rational) request to its kernel family."""

import importlib
from fractions import Fraction

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
    try:
        return module.prove(kind, q, comp, r)
    except NoSolution:
        # 站点在搜索耗尽后仍按数值真假区分报错：命题为假时报"方向反了"而非
        # "未找到解"（实测 arctan 3 > 5/4 → 方向反了；真命题 < 5/4 → 未找到解）。
        # 恒等式 ∫f = ±(C−r) 精确成立，真命题不可能搜出恒≤0 的 P，
        # 故仅在 NoSolution 后补判不会误伤已验证路径。
        # 判定精度是 float64：zeta3/gamma 对 float 相等但方向为假的界仍报
        # "未找到解"（float 差为 0 → 放行进入搜索 → 耗尽），故用 float 比较差。
        from .integrand import constant_mpf
        diff = float(constant_mpf(kind, q)) - float(r)
        claim_false = (diff < 0) if comp == ">" else (diff > 0)
        if claim_false:
            raise WrongDirection
        raise
