"""被积函数重建器：站点 parameters + type/power/comparison → sympy 表达式与积分域。

唯一权威实现。验证脚本 bench/verify.py 与本方求解器的一致性检查都走这里。

返回约定::

    reconstruct(kind, comp, power, params) -> (f, a, b)

- f 是 sympy 表达式（变量 x），即站点积分式右端的完整被积函数，
  已含 tan_q 的 1/cos q、varpi 下界的 1/π 等外层因子；
- (a, b) 为积分区间端点（sympy 对象：0、1、pi/2 或 pi）；
- gamma 型例外：f 含自由整数符号 k（γ 核前因子 x^k 的指数，
  同时是 ln 子证明分母 (1+k·x)^n 的系数；站点不返回该值，
  由验证方枚举 k≥0 数值匹配确定）；
- ln_q/ln_q_square/artanh_q/arcoth_q 的分母幂 s 站点同样不回传，
  f 含自由符号 s（实测 s=n 高频出现，但 arcoth 已观测到 s≠n），
  由验证方枚举 s≥0 数值匹配确定。

参数编码规律（线上 API 实测，见 docs/verify-notes.md）：

- 通用：P 系数取 a_val,b_val,c_val（分数），积分域与核按类型定；
- x² 型核（pi/catalan/zeta3/pi_n/arctan_q/arccot_q）：P=a+b·x²；
- 三角 sin-基（sin_pi_q 等四型）：P=a+b·sin x，且 c_val 复用为
  核频率（sin 类 1-2q，cos 类 2q）；
- artanh_q/arcoth_q：归约为 ln_q'，c_val 复用为约化参数 q'；
- ln_q/ln_q_square：分母幂 s 自由（实测常等于 n，但不保证）；
- varpi/gauss：基固定含 (1-x)，n 不出现在被积函数里；
- tan_q/cot_q/tanh_q/coth_q：整体分别乘 1/cos q、1/sin q、
  1/cosh q、1/sinh q。
"""

from fractions import Fraction

import sympy as sp
from mpmath import mp

x = sp.Symbol("x")
k = sp.Symbol("k", integer=True, nonnegative=True)
s = sp.Symbol("s", integer=True, nonnegative=True)

DOMAIN_HALF_PI = (sp.Integer(0), sp.pi / 2)
DOMAIN_PI = (sp.Integer(0), sp.pi)
DOMAIN_UNIT = (sp.Integer(0), sp.Integer(1))


def frac(v) -> Fraction:
    """站点分数字段（"47/120" 或 int）→ Fraction。"""
    return Fraction(str(v))


def poly(p, t):
    """待定多项式因子 a+b·t+c·t²（t 取 x、x²、x⁴ 或 sin x）。"""
    a, b, c = frac(p["a_val"]), frac(p["b_val"]), frac(p["c_val"])
    return sp.Rational(a.numerator, a.denominator) \
        + sp.Rational(b.numerator, b.denominator) * t \
        + sp.Rational(c.numerator, c.denominator) * t**2


def base01(p):
    """x^m(1-x)^n，[0,1] 上最常见的基。"""
    return x ** int(p["m"]) * (1 - x) ** int(p["n"])


def base_quad(p):
    """x^{2m}(1-x²)^n，偶核族基。"""
    return x ** (2 * int(p["m"])) * (1 - x**2) ** int(p["n"])


def base_sin(p):
    """sin^m x·(1-sin x)^n，[0,π/2]/[0,π] 三角基。"""
    return sp.sin(x) ** int(p["m"]) * (1 - sp.sin(x)) ** int(p["n"])


def reconstruct(kind: str, comp: str, power: Fraction, p: dict):
    """重建被积函数与积分域；见模块 docstring 的返回约定。"""
    q = sp.Rational(power.numerator, power.denominator)
    base = base01(p)
    polyx = poly(p, x)
    polyx2 = poly(p, x**2)
    # sin 基与 x⁴ 基的 P 都只有两项（a+b·t）；c_val 在这些类型里被复用或恒零
    a_, b_ = frac(p["a_val"]), frac(p["b_val"])
    polysin = sp.Rational(a_.numerator, a_.denominator) \
        + sp.Rational(b_.numerator, b_.denominator) * sp.sin(x)
    polyx4 = sp.Rational(a_.numerator, a_.denominator) \
        + sp.Rational(b_.numerator, b_.denominator) * x**4

    if kind == "pi":
        return base_quad(p) * polyx2 / (1 + x**2), *DOMAIN_UNIT
    if kind == "pi_n":
        # π^k：矩空间要求指数奇偶与 k 相反（偶 k 用奇指数取 η(k)）
        e = 2 * int(p["m"]) + (1 if int(power) % 2 == 0 else 0)
        f = x**e * (1 - x**2) ** int(p["n"]) * polyx2 \
            * sp.log(1 / x) ** (int(power) - 1) / (1 + x**2)
        return f, *DOMAIN_UNIT
    if kind in ("e", "e_q"):
        return base * polyx * sp.exp(q * x), *DOMAIN_UNIT
    if kind == "e_pi":
        return base_sin(p) * polysin * sp.exp(x), *DOMAIN_PI
    if kind == "ln_q":
        return base * polyx / (1 + (q - 1) * x) ** s, *DOMAIN_UNIT
    if kind == "ln_q_square":
        f = base * polyx * sp.log(1 + (q - 1) * x) / (1 + (q - 1) * x) ** s
        return f, *DOMAIN_UNIT
    if kind in ("sin_q", "cos_q"):
        return base * polyx * sp.sin(q * x), *DOMAIN_UNIT
    if kind in ("tan_q", "cot_q"):
        outer = 1 / sp.cos(q) if kind == "tan_q" else 1 / sp.sin(q)
        return outer * base * polyx * sp.sin(q * x), *DOMAIN_UNIT
    if kind in ("sinh_q", "cosh_q"):
        return base * polyx * sp.sinh(q * x), *DOMAIN_UNIT
    if kind in ("tanh_q", "coth_q"):
        outer = 1 / sp.cosh(q) if kind == "tanh_q" else 1 / sp.sinh(q)
        return outer * base * polyx * sp.sinh(q * x), *DOMAIN_UNIT
    if kind in ("sin_q_degree", "cos_q_degree", "sin_pi_q", "cos_pi_q"):
        # c_val 复用为核频率：sin 型核 sin((1-2q)x)，cos 型核 sin(2qx)
        return base_sin(p) * polysin * sp.sin(frac(p["c_val"]) * x), *DOMAIN_HALF_PI
    if kind in ("arctan_q", "arccot_q"):
        t = q if kind == "arctan_q" else 1 / q
        return base_quad(p) * polyx2 / (1 + (t * x) ** 2), *DOMAIN_UNIT
    if kind in ("artanh_q", "arcoth_q"):
        # c_val 复用为归约后的 ln 参数 q'，不是多项式系数
        qq = frac(p["c_val"])
        f = base * (sp.Rational(frac(p["a_val"]).numerator, frac(p["a_val"]).denominator)
                    + sp.Rational(frac(p["b_val"]).numerator, frac(p["b_val"]).denominator) * x) \
            / (1 + (qq - 1) * x) ** s
        return f, *DOMAIN_UNIT
    if kind == "golden":
        return base * polyx * sp.sqrt(4 + x), *DOMAIN_UNIT
    if kind == "catalan":
        return base_quad(p) * polyx2 * sp.log(1 / x) / (1 + x**2), *DOMAIN_UNIT
    if kind == "zeta3":
        f = x ** (2 * int(p["m"]) + 1) * (1 - x**2) ** int(p["n"]) \
            * polyx2 * sp.log(x) ** 2 / (1 + x**2)
        return f, *DOMAIN_UNIT
    if kind == "varpi":
        e = 4 * int(p["m"]) + (3 if comp == "<" else 1)
        f = x**e * (1 - x) * polyx4 / sp.sqrt(1 - x**4)
        if comp == ">":
            f /= sp.pi
        return f, *DOMAIN_UNIT
    if kind == "gauss":
        e = 4 * int(p["m"]) + (2 if comp == "<" else 0)
        f = x**e * (1 - x) * polyx4 / sp.sqrt(1 - x**4)
        if comp == ">":
            f /= sp.pi
        return f, *DOMAIN_UNIT
    if kind == "gamma":
        return gamma_integrand(comp, p), *DOMAIN_UNIT
    raise ValueError(f"unsupported type {kind!r}")


def gamma_integrand(comp: str, p: dict):
    """gamma 型被积函数，含自由符号 k（见模块 docstring）。

    下界方向核 L(x)=1/(1-x)+1/ln x-1/2，上界方向核
    U(x)=(2-x)/2-1/(1-x)-1/ln x；恒等式
    ∫₀¹ x^{s-1}(1/(1-x)+1/ln x)dx = γ+ln s 使第一段贡献 γ+ln(k+1)，
    剩余 ln(k+1) 由第二段 ln_{k+1} 子证明吸收（u_val=0 时无第二段）。
    """
    m, n, u = int(p["m"]), int(p["n"]), int(p["u_val"])
    a, b = frac(p["a_val"]), frac(p["b_val"])
    if comp == "<":
        first = x**k * ((2 - x) / 2 - 1 / (1 - x) - 1 / sp.log(x))
    else:
        first = x**k * (1 / (1 - x) + 1 / sp.log(x) - sp.Rational(1, 2))
    if u == 0:
        return first
    sub = x**m * (1 - x) ** n \
        * (sp.Rational(a.numerator, a.denominator)
           + sp.Rational(b.numerator, b.denominator) * x) \
        / (u * (1 + k * x) ** n)
    return first + sub


def constant_mpf(kind: str, power: Fraction):
    """目标常数 C 的 mpmath 值（按当前 mp.dps）。"""
    q = mp.mpf(power.numerator) / power.denominator
    const = {
        "pi": mp.pi * q,
        "e": mp.e * q,
        "pi_n": mp.pi ** int(power),
        "e_q": mp.exp(q),
        "ln_q": mp.log(q),
        "ln_q_square": mp.log(q) ** 2,
        "sin_q": mp.sin(q), "cos_q": mp.cos(q),
        "tan_q": mp.tan(q), "cot_q": 1 / mp.tan(q),
        "sin_q_degree": mp.sin(mp.pi * q / 180),
        "cos_q_degree": mp.cos(mp.pi * q / 180),
        "sin_pi_q": mp.sin(mp.pi * q), "cos_pi_q": mp.cos(mp.pi * q),
        "arctan_q": mp.atan(q), "arccot_q": mp.atan(1 / q),
        "sinh_q": mp.sinh(q), "cosh_q": mp.cosh(q),
        "tanh_q": mp.tanh(q), "coth_q": 1 / mp.tanh(q),
        "artanh_q": mp.atanh(q), "arcoth_q": mp.atanh(1 / q),
        "gamma": mp.euler,
        "golden": (1 + mp.sqrt(5)) / 2,
        "catalan": mp.catalan,
        "zeta3": mp.zeta(3),
        "e_pi": mp.exp(mp.pi * q),
        "varpi": mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi)),
        "gauss": mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi ** 3)),
    }[kind]
    return const


def lhs_mpf(kind: str, comp: str, power: Fraction, rational: Fraction):
    """恒等式左端的 mpmath 值：多数类型是 ±(C−r)，varpi/gauss 有倒数形式。

    站点对 ϖ 的下界证 1−r·ϖ⁻¹>0、对 G 的上界证 r·G⁻¹−1>0
    （对应核给出 {ϖ⁻¹,1}/{G⁻¹,1} 矩空间），其余都是 C−r 或 r−C。
    """
    c = constant_mpf(kind, power)
    r = mp.mpf(rational.numerator) / rational.denominator
    if kind == "varpi" and comp == ">":
        return 1 - r / c
    if kind == "gauss" and comp == "<":
        return r / c - 1
    return (c - r) if comp == ">" else (r - c)
