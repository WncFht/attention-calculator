"""生成 golden case 列表与组合不等式问题清单。

每个 (type, power) 组合用 mpmath 50 位精度算常数真值 v，取三类有理界：
连分数收敛子（分母 ≤ 10^4）、floor/ceil(v·d)/d 网格（d ∈ {10,100,1000}）、
松散整数界；每个界正反两个方向各一个 case——反方向预期返回"方向反了"
错误，同样是有效样本。

真值一律自己算，不信文档；参数是否合法（如 ln_q 要求 q>1）由服务端判定，
非法参数组合的 case 照样生成，记录服务端报错。
"""

from fractions import Fraction

from mpmath import mp, mpf

mp.dps = 50

# 29 个 type 的 power 取值集（Fraction）。pi/e 的 power 是系数（power=3 → 3π），
# 已实测确认；常数类（gamma/golden/catalan/zeta3/e_pi/varpi/gauss）同为系数，
# 取 {0,1,2,3,1/2}；pi_n 的 5/2 是分数幂探针。
PARAM_SETS = {
    "pi": [Fraction(1), Fraction(3), Fraction(8), Fraction(1, 2)],
    "e": [Fraction(1), Fraction(3), Fraction(8), Fraction(1, 2)],
    "pi_n": [Fraction(n) for n in range(2, 9)] + [Fraction(5, 2)],
    "e_q": [
        Fraction(1, 2),
        Fraction(3, 2),
        Fraction(2),
        Fraction(5, 2),
        Fraction(3),
        Fraction(1, 3),
    ],
    "ln_q": [
        Fraction(3, 2),
        Fraction(4, 3),
        Fraction(5, 4),
        Fraction(2),
        Fraction(3),
        Fraction(4),
        Fraction(5),
        Fraction(6),
        Fraction(7),
        Fraction(8),
        Fraction(10),
        Fraction(11),
        Fraction(16),
        Fraction(22),
        Fraction(64),
        Fraction(100),
        Fraction(257),
        Fraction(1, 2),
    ],
    "ln_q_square": [Fraction(2), Fraction(3), Fraction(3, 2), Fraction(5)],
    "sin_q": [
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(1),
        Fraction(6, 5),
        Fraction(2),
        Fraction(3),
    ],
    "cos_q": [
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(1),
        Fraction(6, 5),
        Fraction(2),
        Fraction(3),
    ],
    "tan_q": [Fraction(1, 2), Fraction(2, 3), Fraction(1), Fraction(6, 5)],
    "cot_q": [Fraction(1, 2), Fraction(2, 3), Fraction(1), Fraction(6, 5)],
    "sin_q_degree": [Fraction(d) for d in (10, 15, 18, 30, 36, 40, 45, 54, 72)],
    "cos_q_degree": [Fraction(d) for d in (10, 15, 18, 30, 36, 40, 45, 54, 72)],
    "sin_pi_q": [
        Fraction(1, 10),
        Fraction(1, 8),
        Fraction(1, 7),
        Fraction(1, 6),
        Fraction(1, 5),
        Fraction(1, 4),
        Fraction(3, 10),
        Fraction(1, 3),
        Fraction(2, 5),
    ],
    "cos_pi_q": [
        Fraction(1, 10),
        Fraction(1, 8),
        Fraction(1, 7),
        Fraction(1, 6),
        Fraction(1, 5),
        Fraction(1, 4),
        Fraction(3, 10),
        Fraction(1, 3),
        Fraction(2, 5),
    ],
    "arctan_q": [
        Fraction(1, 5),
        Fraction(1, 3),
        Fraction(1, 2),
        Fraction(1),
        Fraction(2),
        Fraction(3),
    ],
    "arccot_q": [Fraction(1, 2), Fraction(1), Fraction(2)],
    "sinh_q": [
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(2, 3),
        Fraction(1),
        Fraction(2),
        Fraction(5, 2),
    ],
    "cosh_q": [
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(2, 3),
        Fraction(1),
        Fraction(2),
        Fraction(5, 2),
    ],
    "tanh_q": [
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(2, 3),
        Fraction(1),
        Fraction(2),
        Fraction(5, 2),
    ],
    "coth_q": [
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(2, 3),
        Fraction(1),
        Fraction(2),
        Fraction(5, 2),
    ],
    "artanh_q": [Fraction(1, 5), Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)],
    "arcoth_q": [Fraction(3, 2), Fraction(2), Fraction(3)],
    # 常数类 power 是系数 k（证 k·C vs r），非负有理数；0 走进退化/方向预检路径
    "gamma": [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0)],
    "golden": [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0)],
    "catalan": [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0)],
    "zeta3": [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0)],
    "e_pi": [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0)],
    "varpi": [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0)],
    "gauss": [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0)],
}

GRID_DENOMS = (10, 100, 1000)

# 格式/边界探针：负有理数、小数、畸形、零分母、定义域边界、power=0 等。
PROBE_CASES = [
    {"type": "pi", "power": "1", "comparison": "<", "rational": "-3"},
    {"type": "e", "power": "1", "comparison": ">", "rational": "-5/2"},
    {"type": "pi", "power": "1", "comparison": "<", "rational": "3.5"},
    {"type": "pi", "power": "1", "comparison": "<", "rational": "22/7/2"},
    {"type": "pi", "power": "1", "comparison": "<", "rational": "1/0"},
    {"type": "pi", "power": "1", "comparison": ">", "rational": "0"},
    {"type": "ln_q", "power": "1", "comparison": "<", "rational": "1/2"},
    {"type": "sin_q", "power": "0", "comparison": "<", "rational": "1/2"},
    {"type": "artanh_q", "power": "1", "comparison": "<", "rational": "9/10"},
    {"type": "arcoth_q", "power": "1", "comparison": "<", "rational": "2"},
    {"type": "pi_n", "power": "1", "comparison": ">", "rational": "3"},
    {"type": "pi", "power": "0", "comparison": "<", "rational": "1"},
    {"type": "e_q", "power": "0", "comparison": ">", "rational": "1/2"},
]


def true_value(typ, q):
    """常数真值（mpf，50 dps）。q 是 power 字段的 Fraction。"""
    x = mpf(q.numerator) / q.denominator
    if typ == "pi":
        return mp.pi * x
    if typ == "e":
        return mp.e * x
    if typ == "pi_n":
        return mp.pi**x
    if typ == "e_q":
        return mp.exp(x)
    if typ == "ln_q":
        return mp.log(x)
    if typ == "ln_q_square":
        return mp.log(x) ** 2
    if typ == "ln_q_cube":
        return mp.log(x) ** 3
    if typ == "ln_q_quad":
        return mp.log(x) ** 4
    if typ == "arcsin_q":
        return mp.asin(x)
    if typ == "arsinh_q":
        return mp.asinh(x)
    if typ == "gaussint_q":
        return mp.sqrt(mp.pi) / 2 * mp.erf(x)
    if typ == "dawson_q":
        return mp.sqrt(mp.pi) / 2 * mp.exp(-(x**2)) * mp.erfi(x)
    if typ == "erfiint_q":
        return mp.sqrt(mp.pi) / 2 * mp.erfi(x)
    if typ == "li2_q":
        return mp.polylog(2, x)
    if typ == "psi1_q":
        return mp.psi(1, x)
    if typ == "pi_sqrt2":
        return x * mp.pi * mp.sqrt(2)
    if typ == "pi3":
        return x * mp.sqrt(3) * mp.gamma(mpf(1) / 3) ** 3 / (2 * mp.pi)
    if typ == "pi3_u":
        return x * mp.gamma(mpf(2) / 3) ** 2 / mp.gamma(mpf(1) / 3)
    if typ == "pi3_a":
        return x * 2 * mp.pi / mp.sqrt(3)
    if typ == "gamma14":
        return x * mp.gamma(mpf(1) / 4)
    if typ == "gamma34":
        return x * mp.gamma(mpf(3) / 4)
    if typ == "gamma12":
        return x * mp.sqrt(mp.pi)
    if typ == "sin_q":
        return mp.sin(x)
    if typ == "cos_q":
        return mp.cos(x)
    if typ == "tan_q":
        return mp.tan(x)
    if typ == "cot_q":
        return mp.cot(x)
    if typ == "sin_q_degree":
        return mp.sin(x * mp.pi / 180)
    if typ == "cos_q_degree":
        return mp.cos(x * mp.pi / 180)
    if typ == "sin_pi_q":
        return mp.sin(mp.pi * x)
    if typ == "cos_pi_q":
        return mp.cos(mp.pi * x)
    if typ == "arctan_q":
        return mp.atan(x)
    if typ == "arccot_q":
        return mp.atan(1 / x)
    if typ == "sinh_q":
        return mp.sinh(x)
    if typ == "cosh_q":
        return mp.cosh(x)
    if typ == "tanh_q":
        return mp.tanh(x)
    if typ == "coth_q":
        return mp.coth(x)
    if typ == "artanh_q":
        return mp.atanh(x)
    if typ == "arcoth_q":
        return mp.atanh(1 / x)
    if typ == "gamma":
        return x * mp.euler
    if typ == "golden":
        return x * (1 + mp.sqrt(5)) / 2
    if typ == "catalan":
        return x * mp.catalan
    if typ == "zeta3":
        return x * mp.zeta(3)
    if typ == "zeta5":
        return x * mp.zeta(5)
    if typ == "zeta7":
        return x * mp.zeta(7)
    if typ == "zeta9":
        return x * mp.zeta(9)
    if typ == "zeta11":
        return x * mp.zeta(11)
    if typ == "beta4":
        return x * mp.dirichlet(4, [0, 1, 0, -1])
    if typ == "beta6":
        return x * mp.dirichlet(6, [0, 1, 0, -1])
    if typ == "beta8":
        return x * mp.dirichlet(8, [0, 1, 0, -1])
    if typ == "beta10":
        return x * mp.dirichlet(10, [0, 1, 0, -1])
    if typ == "si_q":
        return mp.si(x)
    if typ == "cin_q":
        return mp.euler + mp.log(mp.fabs(x)) - mp.ci(mp.fabs(x))
    if typ == "e_pi":
        return x * mp.exp(mp.pi)
    if typ == "varpi":
        return x * mp.gamma(mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi))
    if typ == "gauss":
        return x / mp.agm(1, mp.sqrt(2))
    raise ValueError(f"unknown type {typ}")


def convergents(v, max_den=10_000):
    """v(>0 mpf) 的连分数收敛子，分母 ≤ max_den。"""
    out = []
    x = v
    h_im2, h_im1 = 0, 1
    k_im2, k_im1 = 1, 0
    for _ in range(60):
        a = int(mp.floor(x))
        h = a * h_im1 + h_im2
        k = a * k_im1 + k_im2
        if k > max_den:
            break
        out.append(Fraction(h, k))
        h_im2, h_im1 = h_im1, h
        k_im2, k_im1 = k_im1, k
        r = x - a
        if abs(r) < mpf(10) ** -40:
            break
        x = 1 / r
    return out


def bounds_for(v):
    """围绕真值 v 的去重有理界列表；v<0 时给固定的非法格式探测集。"""
    if v <= 0:
        # 服务端不接受负有理数：负界记"格式无效"，非负界记参数定义域错误。
        return [Fraction(int(mp.floor(v * 10)), 10), Fraction(0), Fraction(1, 2), Fraction(1)]
    conv = convergents(v)
    if len(conv) > 5:
        conv = conv[:3] + conv[-2:]
    bounds = list(conv)
    for d in GRID_DENOMS:
        bounds.append(Fraction(int(mp.floor(v * d)), d))
        bounds.append(Fraction(int(mp.ceil(v * d)), d))
    bounds.append(Fraction(int(mp.floor(v)), 1))
    bounds.append(Fraction(int(mp.ceil(v)), 1))
    seen = set()
    out = []
    for f in bounds:
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out


def generate_cases():
    """完整 case 列表：[{type, power, comparison, rational}]，按 (type,q) 顺序去重。"""
    cases = []
    seen = set()
    for typ, powers in PARAM_SETS.items():
        for q in powers:
            for r in bounds_for(true_value(typ, q)):
                for comp in ("<", ">"):
                    key = (typ, str(q), comp, str(r))
                    if key not in seen:
                        seen.add(key)
                        cases.append(
                            {"type": typ, "power": str(q), "comparison": comp, "rational": str(r)}
                        )
    for probe in PROBE_CASES:
        key = (probe["type"], probe["power"], probe["comparison"], probe["rational"])
        if key not in seen:
            seen.add(key)
            cases.append(dict(probe))
    return cases


# /decompose_inequality 问题清单：直接基本型、加减、乘除幂、三项以上混合、
# 括号除法、故意写错与格式边界。
COMBO_PROBLEMS = [
    # 直接基本型
    "pi<22/7",
    "pi>3",
    "e>8/3",
    "e<3",
    "pi^2<10",
    "pi^3>31",
    "ln(2)<7/10",
    "ln2<7/10",
    "sin(1)<17/20",
    "cos(1)>1/2",
    "gamma>1/2",
    "phi<8/5",
    "e^pi>23",
    "phi^2>5/2",
    # 两项加减（含故意方向反）
    "pi+e<6",
    "pi+e>6",
    "e-pi>2/5",
    "pi-e>2/5",
    "pi-e<1/2",
    "pi+phi<5",
    "e+gamma<4",
    "pi+ln(2)<4",
    "e+sin(1)>7/2",
    "pi+gamma>7/2",
    "phi+gamma>2",
    "ln(2)+ln(3)<2",
    "sin(1)+cos(1)>1",
    "pi^2+pi>13",
    "2*pi+e<12",
    # 乘除幂
    "pi*e>17/2",
    "e*pi<9",
    "pi*e<8",
    "pi^2+8*pi>35",
    "e*pi+phi+sin(1)<11",
    "pi*phi>5",
    "e/pi<1",
    "pi/e>11/10",
    "ln(2)*e<2",
    "pi^e>22",
    "pi^e<23",
    "e^2+pi>10",
    "2*pi*e>17",
    "pi^2*e>26",
    "gamma*pi<2",
    "sin(1)*pi>5/2",
    "e^2*pi<24",
    "ln(3)*pi>3",
    # 三项及以上混合
    "pi+e+gamma<7",
    "pi+e+phi>6",
    "pi+e+ln(2)+gamma<8",
    "2*pi+3*e+phi>15",
    "pi+e-gamma<6",
    "pi+pi+pi>9",
    "sin(1)+cos(1)+tan(1)<3",
    "ln(2)+ln(3)+ln(5)<4",
    "e*e+pi>10",
    "pi*e+pi*e<18",
    "pi-e+phi>2",
    "3*pi+2*e<15",
    "pi/2+e/3+gamma>3",
    "e-pi-gamma<0",
    # 除法与括号
    "(pi+e)/2<3",
    "pi/(e+1)<1",
    "1/pi+1/e>1/2",
    "(pi+1)/(e+1)>1",
    "e^2/pi>2",
    "pi^2/e>3",
    "ln(10)/pi<3/4",
    "2/(pi+e)>1/3",
    "2/(pi+e)<7/20",
    # 故意写错 / 格式边界
    "pi<>3",
    "pi",
    "pi>3>2",
    "foo>3",
    "pi+>3",
    ">3",
    "pi>",
    "pi>abc",
    "pi>3.14",
    "PI<22/7",
    "pi < 22/7",
    "π<22/7",
    "pi^2.5>10",
    "",
    "ln(1)>0",
    "sin(1)",
    "pi>=3",
]


def combo_problems():
    """去重后的 /decompose_inequality 问题列表。"""
    return list(dict.fromkeys(COMBO_PROBLEMS))


if __name__ == "__main__":
    cases = generate_cases()
    per_type = {}
    for c in cases:
        per_type[c["type"]] = per_type.get(c["type"], 0) + 1
    print(f"total cases: {len(cases)}, combo problems: {len(combo_problems())}")
    for typ, n in per_type.items():
        combos = len(PARAM_SETS[typ])
        print(f"  {typ:14s} {n:4d} cases over {combos} param values")
