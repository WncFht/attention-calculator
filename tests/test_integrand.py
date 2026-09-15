# ruff: noqa: RUF002  # 中文标点属刻意文体
"""integrand.py 重建器的离线测试：fixture 全部来自 zhuyidao.net 实测返回。

每行 (type, comp, power, rational, m, n, a, b, c, u, s_or_k)：
- s_or_k 为 ln 族分母幂 s / gamma 核指数 k 的实测值（从方程文本读出），
  其余类型为 None。
断言：重建被积函数在域上的积分 = 恒等式左端（mpmath 35dps, 容差 1e-25）。
"""

from fractions import Fraction

import sympy as sp
from mpmath import mp

from attention_calculator.integrand import (
    k,
    lhs_mpf,
    reconstruct,
    s,
    x,
)

CASES = [
    ("arccot_q", ">", "2", "9/20", 0, 1, "1/40", "-3/160", "0", "160", None),
    ("arcoth_q", ">", "2", "6/11", 2, 2, "0", "32/33", "3", "33", 2),
    ("arcoth_q", "<", "3", "5/4", 0, 0, "3/4", "5/4", "2", "4", 1),
    ("arctan_q", ">", "1/2", "9/20", 0, 1, "1/40", "-3/160", "0", "160", None),
    ("artanh_q", ">", "1/2", "6/11", 2, 2, "0", "32/33", "3", "33", 2),
    ("artanh_q", ">", "1/3", "1/3", 1, 1, "0", "1/4", "2", "4", 1),
    ("catalan", ">", "1", "9/10", 1, 1, "1/8", "5/8", "0", "8", None),
    ("cos_pi_q", ">", "1/5", "4/5", 1, 1, "873/2500", "-663/2500", "2/5", "2500", None),
    ("cos_q", ">", "1", "1/2", 0, 0, "1/2", "-1", "1/2", "2", None),
    ("cos_q_degree", ">", "30", "17/20", 1, 1, "95/216", "-1/9", "1/3", "216", None),
    ("cos_q_degree", ">", "36", "4/5", 1, 1, "873/2500", "-663/2500", "2/5", "2500", None),
    ("cosh_q", "<", "1", "8/5", 0, 1, "0", "7/10", "-1/10", "10", None),
    ("cot_q", ">", "1", "3/5", 0, 1, "2/5", "-3/10", "-1/10", "10", None),
    ("coth_q", ">", "1", "13/10", 1, 1, "4/15", "-2/15", "-1/30", "30", None),
    ("e", ">", "1", "8/3", 1, 1, "0", "1/3", "0", "3", None),
    ("e_pi", ">", "1", "23", 2, 5, "182650/112211", "-422240/336633", "0", "336633", None),
    ("e_q", ">", "2", "7", 1, 2, "0", "4", "0", "1", None),
    ("gamma", ">", "1", "1/2", 0, 0, "0", "0", "1", "0", 0),
    ("gamma", ">", "1", "4/7", 3, 4, "11/6", "9925/56", "0", "168", 5),
    ("gamma", ">", "1", "5/9", 1, 2, "0", "4/3", "0", "3", 2),
    ("gamma", "<", "1", "3/5", 2, 4, "937/56", "4615/42", "0", "168", 5),
    ("gamma", "<", "1", "7/12", 5, 6, "382287/2", "-104600375/936", "0", "936", 11),
    ("gauss", ">", "1", "4/5", 2, 0, "6/5", "44/5", "0", "5", None),
    ("gauss", "<", "1", "7/8", 1, 1, "0", "15/4", "0", "4", None),
    ("golden", ">", "1", "8/5", 0, 1, "21/800", "-21/800", "0", "800", None),
    ("ln_q", ">", "2", "2/3", 1, 1, "0", "1/2", "0", "2", 1),
    ("ln_q", ">", "3", "1", 0, 2, "1/3", "4/3", "0", "3", 2),
    ("ln_q", "<", "3", "11/10", 2, 3, "2/15", "8/9", "0", "45", 3),
    ("ln_q", ">", "3/2", "2/5", 1, 1, "1/60", "1/20", "0", "60", 1),
    ("ln_q", ">", "5", "8/5", 2, 3, "48/5", "-128/15", "0", "15", 3),
    ("ln_q_square", "<", "2", "1/2", 0, 2, "1/2", "0", "0", "2", 2),
    ("pi", "<", "1", "22/7", 3, 3, "47/120", "-13/120", "0", "120", None),
    ("pi_n", ">", "2", "49/5", 1, 2, "36/65", "816/65", "0", "65", None),
    ("pi_n", ">", "2", "9", 1, 0, "0", "48", "0", "1", None),
    ("pi_n", ">", "3", "31", 2, 3, "338003449/217788864",
     "-97574279/217788864", "0", "217788864", None),
    ("sin_pi_q", "<", "1/3", "9/10", 0, 0, "3/10", "-4/15", "1/3", "30", None),
    ("sin_pi_q", ">", "1/5", "1/2", 0, 1, "1/6", "637/750", "3/5", "750", None),
    ("sin_q", ">", "1", "4/5", 0, 1, "2/5", "-2/5", "1/5", "5", None),
    ("sin_q", ">", "2", "9/10", 2, 2, "9/35", "1/15", "4/21", "105", None),
    ("sin_q_degree", ">", "36", "1/2", 0, 1, "1/6", "637/750", "3/5", "750", None),
    ("sinh_q", ">", "1", "8/7", 0, 1, "2/7", "-2/7", "1/7", "7", None),
    ("tan_q", "<", "1", "8/5", 1, 1, "7/15", "-1/3", "1/15", "15", None),
    ("tanh_q", "<", "1", "10/13", 1, 1, "8/39", "-4/39", "-1/39", "39", None),
    ("varpi", ">", "1", "5/2", 2, 0, "53/7", "26/7", "0", "7", None),
    ("varpi", "<", "1", "11/4", 1, 1, "33/20", "33/4", "0", "20", None),
    ("varpi", "<", "1", "8/3", 6, 1, "2432199/353600", "663927/70720", "0", "353600", None),
    ("zeta3", ">", "1", "7/6", 1, 0, "0", "16/3", "0", "3", None),
]


def params(m, n, a, b, c, u):
    """打包成站点 parameters 形状。"""
    return {"m": m, "n": n, "a_val": a, "b_val": b, "c_val": c, "u_val": u}


def test_all_identities_numeric():
    """全部实测参数样本：重建后数值积分须等于恒等式左端。"""
    mp.dps = 35
    fails = []
    for kind, comp, pw, rat, m, n, av, bv, cv, uv, free in CASES:
        f, a, b = reconstruct(kind, comp, Fraction(pw),
                              params(m, n, av, bv, cv, uv))
        if s in f.free_symbols:
            f = f.subs(s, free)
        if k in f.free_symbols:
            f = f.subs(k, free)
        lhs = lhs_mpf(kind, comp, Fraction(pw), Fraction(rat))
        fn = sp.lambdify(x, f, modules="mpmath")
        val = mp.quad(fn, [mp.mpf(sp.N(a, 30)), mp.mpf(sp.N(b, 30))])
        if abs(val - lhs) > mp.mpf("1e-25") * max(1, abs(lhs)):
            fails.append((kind, comp, rat, mp.nstr(val - lhs, 4)))
    assert not fails, f"恒等式数值不符: {fails}"


def test_pi_canonical_form():
    """m=3,n=3 → x⁶(1-x²)³(47-13x²)/(120(1+x²))（结构级断言）。"""
    f, a, b = reconstruct("pi", "<", Fraction(1), params(
        3, 3, "47/120", "-13/120", "0", "120"))
    expected = x**6 * (1 - x**2) ** 3 * (47 - 13 * x**2) / (120 * (1 + x**2))
    assert sp.simplify(f - expected) == 0
    assert (a, b) == (0, 1)


def test_trig_pi_uses_cval_as_kernel_frequency():
    """sin_pi_q 的 c_val 是核频率 1-2q 而非多项式系数。"""
    f, a, b = reconstruct("sin_pi_q", ">", Fraction(1, 5), params(
        0, 1, "1/6", "637/750", "3/5", "750"))
    expected = (1 - sp.sin(x)) * (sp.Rational(1, 6) + sp.Rational(637, 750) * sp.sin(x)) \
        * sp.sin(sp.Rational(3, 5) * x)
    assert sp.simplify(f - expected) == 0
    assert (a, b) == (0, sp.pi / 2)


def test_varpi_lower_bound_has_pi_and_inverse_lhs():
    """varpi '>'：核含 1/π，左端形式为 1−r·ϖ⁻¹。"""
    f, a, b = reconstruct("varpi", ">", Fraction(1), params(
        2, 0, "53/7", "26/7", "0", "7"))
    expected = x**9 * (1 - x) * (sp.Rational(53, 7) + sp.Rational(26, 7) * x**4) \
        / (sp.pi * sp.sqrt(1 - x**4))
    assert sp.simplify(f - expected) == 0
    assert (a, b) == (0, 1)
    mp.dps = 30
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi))
    assert abs(lhs_mpf("varpi", ">", Fraction(1), Fraction(5, 2))
               - (1 - mp.mpf(5) / 2 / varpi)) < mp.mpf("1e-25")


def test_free_symbols_for_swept_types():
    """ln 族与 gamma 重建结果含自由符号（s / k），其余类型不含。"""
    f, _, _ = reconstruct("ln_q", "<", Fraction(3), params(2, 3, "2/15", "8/9", "0", "45"))
    assert s in f.free_symbols
    f, _, _ = reconstruct("gamma", "<", Fraction(1), params(2, 4, "937/56", "4615/42", "0", "168"))
    assert k in f.free_symbols
    f, _, _ = reconstruct("e", ">", Fraction(1), params(1, 1, "0", "1/3", "0", "3"))
    assert f.free_symbols == {x}
