"""mode=exact 正确性判官的对抗语料生成器（docs/2026-09-16-math-correctness-plan.md W2）。

与 cases.py（站端 parity 语料）分工：那边为站端行为采样，这边为数学正确性
施压——紧界、连分数收敛子、等值边界、域边界、power 变体全部围绕"能否骗过
certified_cmp / 骗出假证明 / 逼出 NoSolution"设计。

每型 power 网格 {1, 2, 3, 1/2, 0, -1}（契约指定；参数型里出域的值照样生成，
用于压域判定），并并入 PARAM_SETS 前两值保证域内覆盖。每个 (type, power)
的界集：

- floor/ceil 整数界；
- 距离网格：d ∈ {1e-1 … 1e-12} 各数量级的 floor(C·d)/d 与 ceil(C·d)/d
  双向夹逼（取代表 decade，见 DECADES）；
- 连分数收敛子的多个深度（浅/中/最深），最深收敛子即 best rational
  approximation，天然覆盖等值边界；
- C <= 0 或不可求值时退固定界集（含 0，power=0 的系数型 C=0 即真等值点）。

外加 EQUAL_PROBES 等值点（Niven 点、e^0、power=0 系数型）、EDGE_PROBES
域边界两侧与已知站端崩溃簇、每型零/负有理界各一。
"""

from fractions import Fraction

from cases import PARAM_SETS, convergents, true_value
from mpmath import mp, workdps

# 契约 power 网格；参数型中出域的值用于压域判定/拒证路径
POWERS = [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2), Fraction(0), Fraction(-1)]

# 夹逼界的 decade 分母：1e-1 … 1e-12 的代表档
DECADES = tuple(10**k for k in (1, 2, 3, 5, 7, 9, 11, 12))

# 等值点探针：Niven（sin(π/6)=cos(π/3)=1/2，角度制 30°/60°）、e^0=1、
# pi^0=1，以及系数型 power=0 → C=0 对 bound=0 的真等值。
EQUAL_PROBES = [
    ("sin_pi_q", "1/6", "1/2"),
    ("cos_pi_q", "1/3", "1/2"),
    ("sin_q_degree", "30", "1/2"),
    ("cos_q_degree", "60", "1/2"),
    ("e_q", "0", "1"),
    ("pi_n", "0", "1"),
    ("ln_q_cube", "1", "0"),  # (ln 1)^3 = 0
    *[
        (k, "0", "0")
        for k in (
            "gamma",
            "golden",
            "catalan",
            "zeta3",
            "e_pi",
            "varpi",
            "gauss",
            "pi",
            "e",
            "zeta5",
            "zeta7",
            "zeta9",
            "zeta11",
            "beta4",
            "beta6",
            "beta8",
            "beta10",
        )
    ],
    # si_q/cin_q 的 q=0 是域外拒证（Si(0)=Cin(0)=0 退化有理），不进等值探针
]

# exact-only 型（EXACT_TYPES：无站端对照，site 路径 400）：zeta5/7 是
# 系数倍率；ln_q_cube/arcsin/arsinh/gauss_erf 三型 power 槽携带 q
EXACT_ONLY = {
    "zeta5": POWERS,
    "zeta7": POWERS,
    "ln_q_cube": [
        Fraction(2),
        Fraction(3),
        Fraction(5),
        Fraction(1),
        Fraction(1, 2),
        Fraction(0),
        Fraction(-1),
    ],
    "arcsin_q": [
        Fraction(1, 2),
        Fraction(1, 4),
        Fraction(3, 4),
        Fraction(1),
        Fraction(0),
        Fraction(-1),
        Fraction(3, 2),
    ],
    "arsinh_q": [Fraction(1), Fraction(2), Fraction(1, 2), Fraction(0), Fraction(-1)],
    "gaussint_q": [Fraction(1), Fraction(2), Fraction(1, 2), Fraction(0), Fraction(-1)],
    "dawson_q": [Fraction(1), Fraction(2), Fraction(1, 2), Fraction(0), Fraction(-1)],
    "erfiint_q": [Fraction(1), Fraction(2), Fraction(1, 2), Fraction(0), Fraction(-1)],
    # 441a7d7 批：pi_sqrt2/pi3 系/复合 Γ 是系数型；li2_q/psi1_q 参数型
    "pi_sqrt2": POWERS,
    "pi3": POWERS,
    "pi3_u": POWERS,
    "pi3_a": POWERS,
    "gamma14": POWERS,
    "gamma34": POWERS,
    "gamma12": POWERS,
    "li2_q": [
        Fraction(1, 2),
        Fraction(-1, 2),
        Fraction(1, 4),
        Fraction(3, 4),
        Fraction(-3, 2),
        Fraction(0),
        Fraction(1),
        Fraction(2),
    ],
    "psi1_q": [
        Fraction(1),
        Fraction(1, 2),
        Fraction(2),
        Fraction(5),
        Fraction(1, 3),
        Fraction(0),
        Fraction(-1),
    ],
    # zeta 更高奇阶与 Dirichlet β 偶阶：同为系数型、无域检（zeta_odd/
    # beta_even 两文件均无 check_input）
    "zeta9": POWERS,
    "zeta11": POWERS,
    "beta4": POWERS,
    "beta6": POWERS,
    "beta8": POWERS,
    "beta10": POWERS,
    # si/cin 参数型：唯一域检 q≠0，q<0 走奇偶归约
    "si_q": POWERS,
    "cin_q": POWERS,
}

# pi_n 新域探针（bb6e35a β/η 生成器）：表外整数、分子 >10 的分数幂、
# 分母边界 64/65
PI_N_PROBES = [Fraction(11), Fraction(12), Fraction(25, 7), Fraction(1, 64), Fraction(1, 65)]

# 域边界两侧与已知站端崩溃簇；(kind, power)
EDGE_PROBES = [
    ("ln_q", "1000000001/1000000000"),  # ln(1+1e-9) ≈ 1e-9
    ("ln_q", "999999999/1000000000"),  # q<1 域外
    ("ln_q", "1"),
    ("ln_q_square", "5"),  # 站端 500 簇：奇异矩系统
    ("ln_q_square", "7"),
    ("artanh_q", "999999/1000000"),  # artanh ≈ 7.25
    ("artanh_q", "1000001/1000000"),  # >1 域外
    ("artanh_q", "0"),  # artanh(0) = 0 等值边界
    ("arcoth_q", "1000001/1000000"),  # arcoth ≈ 7.6
    ("arcoth_q", "1"),  # 域边界
    ("sin_q", "103993/33102"),  # π 收敛子（略小），sin ≈ 5.8e-11
    ("sin_q", "104348/33215"),  # > π 域外
    ("cos_q", "104348/33215"),  # > π/2 域外
    ("tan_q", "312689/99532"),  # ≈ π > π/2 域外
    ("tan_q", "355/113"),  # < π/2 但极近边界
    ("cot_q", "103993/33102"),  # > π/2 域外
    ("sin_pi_q", "499999/1000000"),  # 域内边缘
    ("sin_pi_q", "1/2"),  # 域外（含边界）
    ("cos_pi_q", "1/2"),
    ("sin_q_degree", "89"),
    ("sin_q_degree", "90"),  # 域外
    ("tanh_q", "0"),  # 站端 power=0 crash 簇
    ("coth_q", "0"),
    ("e_q", "0"),
    ("gamma", "0"),
]


def grid_bounds(c) -> list[Fraction]:
    """常数 mpf 值附近的夹逼界集；c 非有限正实数退固定集。"""
    if not isinstance(c, mp.mpf) or not mp.isfinite(c):
        return [Fraction(0), Fraction(1), Fraction(-1), Fraction(1, 2)]
    bounds = [Fraction(int(mp.floor(c)), 1), Fraction(int(mp.ceil(c)), 1)]
    bounds.append(Fraction(0))
    if c > 0:
        for d in DECADES:
            bounds.append(Fraction(int(mp.floor(c * d)), d))
            bounds.append(Fraction(int(mp.ceil(c * d)), d))
        conv = convergents(c, max_den=10**15)
        if conv:
            pick = {1, len(conv) // 2, len(conv) - 1}  # 浅/中/深三个深度
            bounds += [conv[i] for i in sorted(pick) if i < len(conv)]
    else:
        for d in DECADES[:4]:
            bounds.append(Fraction(int(mp.floor(c * d)), d))
            bounds.append(Fraction(int(mp.ceil(c * d)), d))
    return bounds


def generate_cases() -> list[dict]:
    """确定性对抗语料；每条 {type, power, comparison, rational}。"""
    cases, seen = [], set()

    def add(kind, power, comp, rat):
        power, rat = Fraction(power), Fraction(rat)
        key = (kind, str(power), comp, str(rat))
        if key not in seen:
            seen.add(key)
            cases.append(
                {"type": kind, "power": str(power), "comparison": comp, "rational": str(rat)}
            )

    def add_bounds(kind, q, bounds):
        for r in dict.fromkeys(bounds):
            add(kind, q, "<", r)
            add(kind, q, ">", r)

    def const_or_none(kind, q):
        try:
            return true_value(kind, q)
        except Exception:
            return None

    with workdps(200):
        for kind, qs in PARAM_SETS.items():
            powers = list(dict.fromkeys([*POWERS, *qs[:2]]))
            for p in powers:
                add_bounds(kind, p, grid_bounds(const_or_none(kind, p)))
        for kind, powers in EXACT_ONLY.items():
            for p in powers:
                add_bounds(kind, p, grid_bounds(const_or_none(kind, p)))
        for p in PI_N_PROBES:
            add_bounds("pi_n", p, grid_bounds(const_or_none("pi_n", p)))
        for kind, qs in EDGE_PROBES:
            add_bounds(kind, qs, grid_bounds(const_or_none(kind, Fraction(qs))))
    # 等值点：两侧方向 + 微扰界
    for kind, qs, rs in EQUAL_PROBES:
        r = Fraction(rs)
        add_bounds(kind, qs, (r, r + Fraction(1, 10**12), r - Fraction(1, 10**12)))
    # 零/负有理界：每型一个代表参数
    for kind, qs in PARAM_SETS.items():
        add_bounds(kind, qs[0], (Fraction(0), Fraction(-1)))
    return cases


if __name__ == "__main__":
    from collections import Counter

    cases = generate_cases()
    print(f"total adversarial cases: {len(cases)}")
    for typ, n in sorted(Counter(c["type"] for c in cases).items()):
        print(f"  {typ:14s} {n:4d}")
