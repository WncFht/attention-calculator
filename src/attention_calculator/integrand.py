# 中文标点属刻意文体
"""被积函数重建器：站点 parameters + type/power/comparison → sympy 表达式与积分域。

唯一权威实现。验证脚本 bench/verify.py 与本方求解器的一致性检查都走这里。

返回约定::

    reconstruct(kind, comp, power, params) -> (f, a, b)

- f 是 sympy 表达式（变量 x），即站点积分式右端的完整被积函数，
  已含 tan_q 的 1/cos q、varpi 下界的 1/π 等外层因子；
- (a, b) 为积分区间端点（sympy 对象：0、1、pi/2 或 pi）；
- gamma 型的核指数 k 即参数 cu_val（渲染式实证 k=cu_val，
  如 cu_val=16→x^16）；不再是自由符号；
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
    return (
        sp.Rational(a.numerator, a.denominator)
        + sp.Rational(b.numerator, b.denominator) * t
        + sp.Rational(c.numerator, c.denominator) * t**2
    )


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
    polyx1 = (
        sp.Rational(a_.numerator, a_.denominator) + sp.Rational(b_.numerator, b_.denominator) * x
    )
    polysin = sp.Rational(a_.numerator, a_.denominator) + sp.Rational(
        b_.numerator, b_.denominator
    ) * sp.sin(x)
    polyx4 = (
        sp.Rational(a_.numerator, a_.denominator) + sp.Rational(b_.numerator, b_.denominator) * x**4
    )

    if kind == "pi":
        return base_quad(p) * polyx2 / (1 + x**2), *DOMAIN_UNIT
    if kind == "pi_n":
        # π^{u/v} 先升到 v 次幂变 π^u 再证：核 ln^{u-1}(1/x)，LHS 为 ±(π^u−r^v)。
        # 矩空间要求基指数奇偶与 u 相反（偶 u 用奇指数取 η(u)）。
        u = power.numerator
        e = 2 * int(p["m"]) + (1 if u % 2 == 0 else 0)
        f = x**e * (1 - x**2) ** int(p["n"]) * polyx2 * sp.log(1 / x) ** (u - 1) / (1 + x**2)
        return f, *DOMAIN_UNIT
    if kind == "e":
        # e 型的 power 是系数（k·e），核恒为 e^x；e_q 的 power 才是指数
        return base * polyx * sp.exp(x), *DOMAIN_UNIT
    if kind == "e_q":
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
        return base * polyx1 / (1 + (qq - 1) * x) ** s, *DOMAIN_UNIT
    if kind == "golden":
        return base * polyx * sp.sqrt(4 + x), *DOMAIN_UNIT
    if kind == "catalan":
        return base_quad(p) * polyx2 * sp.log(1 / x) / (1 + x**2), *DOMAIN_UNIT
    if kind == "zeta3":
        f = (
            x ** (2 * int(p["m"]) + 1)
            * (1 - x**2) ** int(p["n"])
            * polyx2
            * sp.log(x) ** 2
            / (1 + x**2)
        )
        return f, *DOMAIN_UNIT
    if kind == "varpi":
        # '>' 兜底（cu_val=1 标记）：∫ (au/u)·x(1-x)·√(1-x⁴)/π dx + b_val，
        # 而非 1/√ 核矩空间；正常路径即使 a_val=0（如 varpi>0 的 (0,1)）
        # 也走下方矩空间核。
        if comp == ">" and int(p["cu_val"]) != 0:
            f = frac(p["au_val"]) / frac(p["u_val"]) * x * (1 - x) * sp.sqrt(
                1 - x**4
            ) / sp.pi + frac(p["b_val"])
            return f, *DOMAIN_UNIT
        e = 4 * int(p["m"]) + (3 if comp == "<" else 1)
        f = x**e * (1 - x) * polyx4 / sp.sqrt(1 - x**4)
        if comp == ">":
            f /= sp.pi
        return f, *DOMAIN_UNIT
    if kind == "gauss":
        # 同上 '>' 兜底：∫ au·(1-x)·√(1-x⁴)/π dx + b_val；正常路径
        # a_val=0 的记录（如 gauss>3/4 的 (0,6)）cu_val=0，不受影响。
        if comp == ">" and int(p["cu_val"]) != 0:
            f = frac(p["au_val"]) * (1 - x) * sp.sqrt(1 - x**4) / sp.pi + frac(p["b_val"])
            return f, *DOMAIN_UNIT
        e = 4 * int(p["m"]) + (2 if comp == "<" else 0)
        f = x**e * (1 - x) * polyx4 / sp.sqrt(1 - x**4)
        if comp == ">":
            f /= sp.pi
        return f, *DOMAIN_UNIT
    if kind == "gamma":
        return gamma_integrand(comp, p), *DOMAIN_UNIT
    raise ValueError(f"unsupported type {kind!r}")


def gamma_integrand(comp: str, p: dict):
    """gamma 型被积函数（k 即 cu_val，见模块 docstring；无自由符号）。

    第一段 x^{cu_val}·K_dir(x)：'<' 用 L=1/(1-x)+1/ln x-1/2，
    '>' 用 U=(2-x)/2-1/(1-x)-1/ln x；∫x^{s-1}(1/(1-x)+1/ln x)dx=γ+ln s
    使首段贡献 γ±ln(cu_val+1)。第二段是 ln_{cu_val+1} 子证明
    x^m(1-x)^n(au+bu·x)/(u·(1+cu·x)^e)，分母幂 e 实测 '<' 取 n、
    '>' 取 cu_val（后者仅 r=57/100 一单样本）。u_val=0 时无第二段，
    改为加性常数 a_val（a,b 在 u≠0 时由 au/u,bu/u 等价给出）。
    """
    m, n, u = int(p["m"]), int(p["n"]), int(p["u_val"])
    kk = int(p["cu_val"])
    a = frac(p["a_val"])
    au, bu = frac(p["au_val"]), frac(p["bu_val"])
    y = 1 - x  # 合并成单分式消去 1/(1-x) 与 1/ln x 在 x→1 的灾难性对消
    if comp == "<":
        first = x**kk * ((y**2 + y - 2) * sp.log(x) - 2 * y) / (2 * y * sp.log(x))
    else:
        first = x**kk * ((2 - y) * sp.log(x) + 2 * y) / (2 * y * sp.log(x))
    if u == 0:
        return first + sp.Rational(a.numerator, a.denominator)
    e = n if comp == "<" else kk
    sub = (
        x**m
        * (1 - x) ** n
        * (
            sp.Rational(au.numerator, au.denominator)
            + sp.Rational(bu.numerator, bu.denominator) * x
        )
        / (u * (1 + kk * x) ** e)
    )
    return first + sub


def constant_mpf(kind: str, power: Fraction):
    """目标常数 C 的 mpmath 值（按当前 mp.dps）。

    惰性查表：字典里是 zero-arg lambda，只算请求的那一项——q=0 时
    cot/coth/arccot 等倒数常数会 ZeroDivision，惰求值避免无关类型被波及。
    """
    q = mp.mpf(power.numerator) / power.denominator
    const = {
        "pi": lambda: mp.pi * q,
        "e": lambda: mp.e * q,
        "pi_n": lambda: mp.pi**q,
        "e_q": lambda: mp.exp(q),
        "ln_q": lambda: mp.log(q),
        "ln_q_square": lambda: mp.log(q) ** 2,
        "sin_q": lambda: mp.sin(q),
        "cos_q": lambda: mp.cos(q),
        "tan_q": lambda: mp.tan(q),
        "cot_q": lambda: 1 / mp.tan(q),
        "sin_q_degree": lambda: mp.sin(mp.pi * q / 180),
        "cos_q_degree": lambda: mp.cos(mp.pi * q / 180),
        "sin_pi_q": lambda: mp.sin(mp.pi * q),
        "cos_pi_q": lambda: mp.cos(mp.pi * q),
        "arctan_q": lambda: mp.atan(q),
        "arccot_q": lambda: mp.atan(1 / q),
        "sinh_q": lambda: mp.sinh(q),
        "cosh_q": lambda: mp.cosh(q),
        "tanh_q": lambda: mp.tanh(q),
        "coth_q": lambda: 1 / mp.tanh(q),
        "artanh_q": lambda: mp.atanh(q),
        "arcoth_q": lambda: mp.atanh(1 / q),
        "gamma": lambda: q * mp.euler,
        "golden": lambda: q * (1 + mp.sqrt(5)) / 2,
        # catalan/zeta3 的 power 同样是系数（kernel-spec: "coef of C"——
        # 站端 zeta3 0 > 1 判方向反了即 0·ζ(3)>1 为假），漏乘会让 q=0 的
        # 方向兜底把假命题放行为未找到。
        "catalan": lambda: q * mp.catalan,
        "zeta3": lambda: q * mp.zeta(3),
        # exact-mode-only types (kernels/EXACT_TYPES); power is the coefficient
        "zeta5": lambda: q * mp.zeta(5),
        "zeta7": lambda: q * mp.zeta(7),
        "zeta9": lambda: q * mp.zeta(9),
        "zeta11": lambda: q * mp.zeta(11),
        # ln_q_cube 的 power 槽携带 q 本身（与 ln_q 同），常数是 (ln q)^3
        "ln_q_cube": lambda: mp.log(q) ** 3,
        "arcsin_q": lambda: mp.asin(q),
        "arsinh_q": lambda: mp.asinh(q),
        # gauss-erf 族：缩放常数（erf/erfi 本身因 √π 障碍不可达，见
        # docs/2026-09-16-w3-research-erf.md）；power 槽携带 q
        "gaussint_q": lambda: mp.sqrt(mp.pi) * mp.erf(q) / 2,
        # 本环境 mpmath 无 mp.dawson：F(q) = √π/2 · e^{−q²}·erfi(q)
        "dawson_q": lambda: mp.sqrt(mp.pi) * mp.exp(-(q * q)) * mp.erfi(q) / 2,
        "erfiint_q": lambda: mp.sqrt(mp.pi) * mp.erfi(q) / 2,
        # lemniscate 余元常数 S = π√2 = Γ(1/4)Γ(3/4)；power 是系数
        "pi_sqrt2": lambda: q * mp.pi * mp.sqrt(2),
        # Dixon/Γ(1/3) 格点族（kernels/dixon）：power 均是系数 q·C
        # π₃ = B(1/3,1/3) = √3·Γ(1/3)³/(2π)；U = Γ(2/3)²/Γ(1/3)；A = 2π/√3
        "pi3": lambda: q * mp.sqrt(3) * mp.gamma(mp.mpf(1) / 3) ** 3 / (2 * mp.pi),
        "pi3_u": lambda: q * mp.gamma(mp.mpf(2) / 3) ** 2 / mp.gamma(mp.mpf(1) / 3),
        "pi3_a": lambda: q * 2 * mp.pi / mp.sqrt(3),
        # li2_q 的 power 槽携带 q 本身（与 ln_q 同），常数是 Li_2(q)
        "li2_q": lambda: mp.polylog(2, q),
        # psi1_q 同理：power 即参数 q，常数是 ψ′(q)，不乘系数
        "psi1_q": lambda: mp.psi(1, q),
        # si_q/cin_q：power 即参数 q；Cin 是偶函数，核内用 |q| 归约
        "si_q": lambda: mp.si(q),
        "cin_q": lambda: mp.euler + mp.log(mp.fabs(q)) - mp.ci(mp.fabs(q)),
        # Γ 特殊值复合型（kernels/gamma_special）：power 是系数 q·C
        "gamma14": lambda: q * mp.gamma(mp.mpf(1) / 4),
        "gamma34": lambda: q * mp.gamma(mp.mpf(3) / 4),
        "gamma12": lambda: q * mp.sqrt(mp.pi),
        "e_pi": lambda: q * mp.exp(mp.pi),  # power 是 e^π 的系数，非指数
        # varpi/gauss 的 power 是常数倍率（LHS 形如 q·G、q·ϖ），
        # 与 e/golden 一致；此前漏乘导致 q≠1 时真假判定错。
        "varpi": lambda: q * mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi)),
        "gauss": lambda: q * mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi**3)),
    }[kind]()
    return const


def lhs_mpf(kind: str, comp: str, power: Fraction, rational: Fraction):
    """恒等式左端的 mpmath 值：多数类型是 ±(C−r)，varpi/gauss 有倒数形式。

    站点对 ϖ 的下界证 1−r·ϖ⁻¹>0、对 G 的上界证 r·G⁻¹−1>0
    （对应核给出 {ϖ⁻¹,1}/{G⁻¹,1} 矩空间），其余都是 C−r 或 r−C。
    """
    c = constant_mpf(kind, power)
    r = mp.mpf(rational.numerator) / rational.denominator
    if kind == "pi_n" and power.denominator != 1:
        # 分数幂升幂后比较 C^v 与 r^v
        v = power.denominator
        c, r = c**v, r**v
    if kind == "varpi" and comp == ">":
        return 1 - r / c
    if kind == "gauss" and comp == "<":
        return r / c - 1
    return (c - r) if comp == ">" else (r - c)
