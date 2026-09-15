"""Trig family tests: moments vs the article's printed expansions, solver
vs observed zhuyidao.net responses, and numeric validation of the identities."""

import re
from fractions import Fraction as F
from math import comb

import mpmath as mp
import pytest
import sympy as sp

from attention_calculator.engine import WrongDirection
from attention_calculator.kernels import trig_pi, trig_q

q_, a_, b_, c_ = sp.symbols("q a b c")


def symbolic_moment_trig_q(m, n, j):
    """Closed-form moment with symbolic q, mirroring trig_q's recurrence."""
    # S_k over basis {sin q, cos q, 1} as sympy expressions
    S = {0: {"sin": sp.Rational(0), "cos": -1 / q_, "1": 1 / q_},
         1: {"sin": 1 / q_**2, "cos": -1 / q_, "1": sp.Rational(0)}}
    for k in range(2, m + n + j + 1):
        f = -k * (k - 1) / q_**2
        S[k] = {"sin": sp.simplify(k / q_**2 + f * S[k - 2]["sin"]),
                "cos": sp.simplify(-1 / q_ + f * S[k - 2]["cos"]),
                "1": sp.simplify(f * S[k - 2]["1"])}
    out = {k: sp.Rational(0) for k in ("sin", "cos", "1")}
    for i in range(n + 1):
        cf = sp.Rational((-1) ** i * comb(n, i))
        for k in out:
            out[k] += cf * S[m + i + j][k]
    return {k: sp.expand(v) for k, v in out.items()}


def test_trig_q_n1_article_expansion():
    """Article type 14, n=1: printed sin/cos coefficients must match exactly.

    (The article's constant part is printed with flipped signs -- the correct
    value is ((2a-2b)q^2 - 24c)/q^5, verified numerically; we assert ours.)
    """
    m = n = 1
    M = [symbolic_moment_trig_q(m, n, j) for j in range(3)]
    sin_coef = sp.expand(a_ * M[0]["sin"] + b_ * M[1]["sin"] + c_ * M[2]["sin"])
    cos_coef = sp.expand(a_ * M[0]["cos"] + b_ * M[1]["cos"] + c_ * M[2]["cos"])
    one_coef = sp.expand(a_ * M[0]["1"] + b_ * M[1]["1"] + c_ * M[2]["1"])
    art_sin = (-(a_ + b_ + c_) * q_**3 + (6 * b_ + 18 * c_) * q_) / q_**5
    art_cos = (-(2 * a_ + 4 * b_ + 6 * c_) * q_**2 + 24 * c_) / q_**5
    our_one = ((2 * a_ - 2 * b_) * q_**2 - 24 * c_) / q_**5
    assert sp.simplify(sin_coef - art_sin) == 0
    assert sp.simplify(cos_coef - art_cos) == 0
    assert sp.simplify(one_coef - our_one) == 0


def test_trig_q_moments_numeric():
    """S_k recurrence agrees with direct numerical integration."""
    mp.mp.dps = 40
    qf = mp.mpf("0.7")
    s = trig_q.sin_moments(10, F(7, 10))
    for k, mom in enumerate(s):
        exact = (mom.get("sin_q", F(0)) * mp.sin(qf)
                 + mom.get("cos_q", F(0)) * mp.cos(qf) + mom.get("1", F(0)))
        direct = mp.quad(lambda x, kk=k: mp.mpf(x) ** kk * mp.sin(qf * x), [0, 1])
        assert abs(exact - direct) < mp.mpf("1e-30")


def symbolic_T(j):
    """trig_pi.angle_moment with symbolic alpha, keys {C, 1} -> sympy exprs."""
    al = sp.Symbol("al", positive=True)
    out = {"C": sp.Rational(0), "1": sp.Rational(0)}
    if j % 2:
        for i in range((j + 1) // 2):
            k = j - 2 * i
            cf = sp.Rational((-1) ** ((j - 1) // 2 - i) * comb(j, i), 2 ** j)
            for s, sign in ((-1, 1), (1, -1)):
                out["C"] += sign * cf * s * (-1) ** ((k - 1) // 2) / (al + s * k)
    else:
        c0 = sp.Rational(comb(j, j // 2), 2 ** j)
        out["1"] += c0 / al
        out["C"] += -c0 / al
        for i in range(j // 2):
            k = j - 2 * i
            cf = sp.Rational((-1) ** (j // 2 - i) * comb(j, i), 2 ** j)
            for s in (1, -1):
                beta = al + s * k
                out["1"] += cf / beta
                out["C"] += -cf * (-1) ** (k // 2) / beta
    return out


def symbolic_basis_trig_pi(m, n, j):
    out = {"C": sp.Rational(0), "1": sp.Rational(0)}
    for i in range(n + 1):
        cf = sp.Rational((-1) ** i * comb(n, i))
        t = symbolic_T(m + i + j)
        out["C"] += cf * t["C"]
        out["1"] += cf * t["1"]
    return out


def test_trig_pi_n1_article_expansions():
    """Article types 15 and 17 printed n=1 expansions, checked verbatim."""
    a, b = a_, b_
    # type 15: alpha = 1 - 2q, constant sin(pi q)
    M0, M1 = (symbolic_basis_trig_pi(1, 1, j) for j in range(2))
    al = sp.Symbol("al", positive=True)
    M0 = {k: v.subs(al, 1 - 2 * q_) for k, v in M0.items()}
    M1 = {k: v.subs(al, 1 - 2 * q_) for k, v in M1.items()}
    sin_coef = sp.expand(a * M0["C"] + b * M1["C"])
    one_coef = sp.expand(a * M0["1"] + b * M1["1"])
    den = 8 * q_ * (8 * q_**6 - 28 * q_**5 + 14 * q_**4 + 35 * q_**3
                    - 28 * q_**2 - 7 * q_ + 6)
    art_sin = ((-8 * a - 8 * b) * q_**4 + (16 * a + 16 * b) * q_**3
               + (2 * a - 10 * b) * q_**2 + (-10 * a + 2 * b) * q_
               + 12 * a - 9 * b) / den
    art_one = ((16 * b - 16 * a) * q_**4 + (-32 * b + 32 * a) * q_**3
               + (16 * a - 16 * b) * q_**2 + (32 * b - 32 * a) * q_) / den
    assert sp.simplify(sin_coef - art_sin) == 0
    assert sp.simplify(one_coef - art_one) == 0

    # type 17: alpha = 2q, constant cos(pi q)
    M0, M1 = (symbolic_basis_trig_pi(1, 1, j) for j in range(2))
    M0 = {k: v.subs(al, 2 * q_) for k, v in M0.items()}
    M1 = {k: v.subs(al, 2 * q_) for k, v in M1.items()}
    cos_coef = sp.expand(a * M0["C"] + b * M1["C"])
    one_coef = sp.expand(a * M0["1"] + b * M1["1"])
    den17 = 4 * q_ * (16 * q_**6 - 56 * q_**4 + 49 * q_**2 - 9)
    art_cos = ((8 * a + 8 * b) * q_**4 + (-14 * a - 2 * b) * q_**2
               + (-9 * a + 9 * b)) / den17
    art_one17 = ((16 * a - 16 * b) * q_**4 + (-40 * a + 40 * b) * q_**2
                 + (9 * a - 9 * b)) / den17
    assert sp.simplify(cos_coef - art_cos) == 0
    assert sp.simplify(one_coef - art_one17) == 0


def test_trig_pi_moments_numeric():
    """T_j product-to-sum moments agree with direct numerical integration."""
    mp.mp.dps = 40
    alpha = F(3, 8)
    th = mp.pi * mp.mpf(3) / 8 / 2
    for j in range(12):
        d = trig_pi.angle_moment(j, alpha)
        exact = d.get("C", F(0)) * mp.cos(th) + d.get("1", F(0))
        direct = mp.quad(lambda x, jj=j: mp.sin(x) ** jj * mp.sin(mp.mpf(3) / 8 * x),
                         [0, mp.pi / 2])
        assert abs(exact - direct) < mp.mpf("1e-30")


# --- solver parity: parameters observed on zhuyidao.net --------------------

SITE = [
    # (kind, power, comp, bound) -> (m, n, a_val, b_val, c_val, u_val)
    ("sin_q", "1", ">", "5/6", (0, 1, "1/6", "-1/3", "1/6", "6")),
    ("sin_q", "3/5", ">", "1/3", (0, 0, "31/45", "1/25", "4/25", "225")),
    ("sin_q", "1", "<", "6/7", (0, 1, "0", "2/7", "-1/7", "7")),
    ("sin_q", "3", ">", "1/8", (2, 2, "283/272", "-339/136", "759/272", "272")),
    ("sin_q", "1", ">", "0", (0, 0, "2", "-1", "1", "1")),
    ("cos_q", "1", "<", "11/20", (1, 1, "3/40", "1/10", "-1/40", "40")),
    ("cos_q", "1", ">", "8/15", (1, 1, "13/90", "-11/90", "2/45", "90")),
    ("tan_q", "1", "<", "8/5", (1, 1, "7/15", "-1/3", "1/15", "15")),
    ("tan_q", "6/5", ">", "2", (0, 1, "4/5", "24/125", "96/125", "125")),
    ("tan_q", "1", ">", "3/2", (0, 1, "0", "1/4", "1/4", "4")),
    ("cot_q", "1", ">", "3/5", (0, 1, "2/5", "-3/10", "-1/10", "10")),
    ("cot_q", "1", "<", "2/3", (0, 1, "0", "1/6", "1/6", "6")),
    ("cot_q", "6/5", "<", "1/2", (0, 1, "2/5", "12/125", "48/125", "125")),
    ("sin_pi_q", "1/5", ">", "1/2", (0, 1, "1/6", "637/750", "3/5", "750")),
    ("sin_pi_q", "1/3", ">", "6/7", (1, 1, "65/189", "-40/189", "1/3", "189")),
    ("cos_pi_q", "1/7", "<", "10/11", (0, 1, "15/77", "-480/3773", "2/7", "3773")),
    ("cos_pi_q", "2/9", ">", "3/4", (1, 1, "9793/23328", "-4991/23328", "4/9", "23328")),
    ("sin_q_degree", "18", "<", "1/3", (0, 1, "1/5", "-14/125", "4/5", "125")),
    ("sin_q_degree", "45", ">", "7/10", (1, 1, "21/64", "-21/64", "1/2", "64")),
    ("cos_q_degree", "40", ">", "3/4", (1, 1, "9793/23328", "-4991/23328", "4/9", "23328")),
    ("cos_q_degree", "80", "<", "1/5", (0, 1, "7/45", "-26/729", "8/9", "3645")),
]


@pytest.mark.parametrize("kind,power,comp,bound,expected", SITE)
def test_site_params(kind, power, comp, bound, expected):
    res = trig_q.prove(kind, F(power), comp, F(bound)) if kind in (
        "sin_q", "cos_q", "tan_q", "cot_q") else trig_pi.prove(
        kind, F(power), comp, F(bound))
    p = res["parameters"]
    m, n, av, bv, cv, u = expected
    assert (p["m"], p["n"], p["a_val"], p["b_val"], p["c_val"], p["u_val"]) == \
        (m, n, av, bv, cv, u)


def numeric_check_trig_q(kind, q, comp, bound, m, n, au, bu, cu, u):
    """Verify the identity: integral equals the claimed difference, and the
    integrand is sign-definite on [0,1]."""
    mp.mp.dps = 40
    qf = mp.mpf(q.numerator) / q.denominator

    def integ(x):
        return x**m * (1 - x) ** n * (au + bu * x + cu * x * x) * mp.sin(qf * x)
    val = mp.quad(integ, [0, 1]) / u
    sgn = 1 if comp == ">" else -1
    # rendered identity: sgn*(C - bound) = val/divisor, divisor = cos q (tan),
    # sin q (cot) or 1 -- the "p*cos - sin" structure is divided out.
    if kind == "sin_q":
        C, divisor = mp.sin(qf), 1
    elif kind == "cos_q":
        C, divisor = mp.cos(qf), 1
    elif kind == "tan_q":
        C, divisor = mp.tan(qf), mp.cos(qf)
    else:
        C, divisor = 1 / mp.tan(qf), mp.sin(qf)
    expected = sgn * (C - mp.mpf(bound.numerator) / bound.denominator) * divisor
    assert abs(val - expected) < mp.mpf("1e-20")
    # integrand sign check on a grid (poly is continuous, roots isolated)
    xs = [mp.mpf(i) / 200 for i in range(201)]
    assert all(integ(x) >= 0 for x in xs)


@pytest.mark.parametrize("kind,power,comp,bound,expected", SITE)
def test_identity_numeric(kind, power, comp, bound, expected):
    if kind not in ("sin_q", "cos_q", "tan_q", "cot_q"):
        pytest.skip("trig_pi checked separately")
    mod = trig_q
    res = mod.prove(kind, F(power), comp, F(bound))
    p = res["parameters"]
    numeric_check_trig_q(kind, F(power), comp, F(bound), p["m"], p["n"],
                         F(p["au_val"]), F(p["bu_val"]), F(p["cu_val"]),
                         F(p["u_val"]))


def numeric_check_trig_pi(kind, power, comp, bound, m, n, au, bu, u, alpha):
    mp.mp.dps = 40
    alf = mp.mpf(alpha.numerator) / alpha.denominator
    q_deg = F(power)
    q_eff = q_deg / 180 if kind.endswith("_degree") else q_deg
    C = mp.sin(mp.pi * mp.mpf(q_eff.numerator) / q_eff.denominator) if kind.startswith("sin") \
        else mp.cos(mp.pi * mp.mpf(q_eff.numerator) / q_eff.denominator)
    def integ(x):
        return (mp.sin(x) ** m * (1 - mp.sin(x)) ** n
                * (au + bu * mp.sin(x)) * mp.sin(alf * x))
    val = mp.quad(integ, [0, mp.pi / 2]) / u
    sgn = 1 if comp == ">" else -1
    assert abs(val - sgn * (C - mp.mpf(bound.numerator) / bound.denominator)) < mp.mpf("1e-20")
    xs = [mp.pi / 2 * i / 200 for i in range(201)]
    assert all(integ(x) >= 0 for x in xs)


@pytest.mark.parametrize("kind,power,comp,bound,expected", SITE)
def test_identity_numeric_pi(kind, power, comp, bound, expected):
    if kind in ("sin_q", "cos_q", "tan_q", "cot_q"):
        pytest.skip("trig_q checked separately")
    res = trig_pi.prove(kind, F(power), comp, F(bound))
    p = res["parameters"]
    numeric_check_trig_pi(kind, F(power), comp, F(bound), p["m"], p["n"],
                          F(p["au_val"]), F(p["bu_val"]), F(p["u_val"]),
                          F(p["c_val"]))


# --- errors -----------------------------------------------------------------

@pytest.mark.parametrize("kind,power,bound,msg", [
    ("sin_q", "4", "1", "请在sin后输入一个在(0,π)内的数"),
    ("sin_q", "0", "1", "请在sin后输入一个在(0,π)内的数"),
    ("sin_q", "22/7", "1", "请在sin后输入一个在(0,π)内的数"),
    ("cos_q", "4", "1", "请在cos后输入一个在(0,π/2)内的数"),
    ("tan_q", "2", "1", "请在tan后输入一个在(0,π/2)内的数"),
    ("cot_q", "3", "0", "请在cot后输入一个在(0,π/2)内的数"),
    ("sin_q", "-1", "0", "左侧系数格式无效"),
    ("sin_q", "4", "-3/4", "右侧有理数格式无效"),
    ("sin_pi_q", "1/2", "9/10", "请在输入一个在(0,1/2)内的分数，本情况不支持整数"),
    ("sin_pi_q", "0", "1/2", "请在输入一个在(0,1/2)内的分数，本情况不支持整数"),
    ("sin_pi_q", "3/6", "1/2", "请在输入一个在(0,1/2)内的分数，本情况不支持整数"),
    ("cos_pi_q", "1", "1", "请在输入一个在(0,1/2)内的分数，本情况不支持整数"),
    ("sin_q_degree", "90", "1", "请在输入一个在(0,90)内的数"),
    ("cos_q_degree", "200", "0", "请在输入一个在(0,90)内的数"),
    ("sin_pi_q", "1/6", "1/2", "二者相等"),
    ("cos_pi_q", "1/3", "1/2", "二者相等"),
    ("sin_q_degree", "30", "1/2", "二者相等"),
    ("cos_q_degree", "60", "1/2", "二者相等"),
])
def test_input_errors(kind, power, bound, msg):
    mod = trig_q if kind in ("sin_q", "cos_q", "tan_q", "cot_q") else trig_pi
    with pytest.raises(ValueError, match=re.escape(msg)):
        mod.prove(kind, F(power), ">", F(bound))


def test_wrong_direction():
    with pytest.raises(WrongDirection):
        trig_q.prove("sin_q", F(1), "<", F(5, 6))
    with pytest.raises(WrongDirection):
        trig_pi.prove("cos_q_degree", F(45), "<", F(7, 10))


# --- rendering: verbatim strings observed on the site -----------------------

RENDER = [
    ("sin_q", "1", ">", "5/6",
     "\\sin1 - \\dfrac{5}{6} = \\int_0^1 \\frac{\\left(1 - x\\right) "
     "\\left(x^{2} - 2 x + 1\\right) \\sin{\\left(x \\right)}}{6} \\mathrm{d} x > 0"),
    ("sin_q", "3", ">", "1/8",
     "\\sin3 - \\dfrac{1}{8} = \\int_0^1 \\frac{x^{2} \\left(1 - x\\right)^{2} "
     "\\cdot \\left(759 x^{2} - 678 x + 283\\right) \\sin{\\left(3 x \\right)}}"
     "{272} \\mathrm{d} x > 0"),
    ("sin_q", "1", ">", "0",
     "\\sin1 - 0 = \\int_0^1 \\left(x^{2} - x + 2\\right) \\sin{\\left(x \\right)}"
     " \\mathrm{d} x > 0"),
    ("cos_q", "1", "<", "11/20",
     "\\dfrac{11}{20} - \\cos1 = \\int_0^1 \\frac{x \\left(1 - x\\right) "
     "\\left(- x^{2} + 4 x + 3\\right) \\sin{\\left(x \\right)}}{40} \\mathrm{d} x > 0"),
    ("tan_q", "6/5", ">", "2",
     "\\tan\\dfrac{6}{5} - 2 = \\int_0^1 \\dfrac{1}{\\cos(6/5)} "
     "\\frac{\\left(1 - x\\right) \\left(96 x^{2} + 24 x + 100\\right) "
     "\\sin{\\left(\\frac{6 x}{5} \\right)}}{125} \\mathrm{d} x > 0"),
    ("cot_q", "1", "<", "2/3",
     "\\dfrac{2}{3} - \\cot1 = \\int_0^1 \\dfrac{1}{\\sin(1)} \\frac{\\left(1 - x\\right) "
     "\\left(x^{2} + x\\right) \\sin{\\left(x \\right)}}{6} \\mathrm{d} x > 0"),
    ("sin_pi_q", "1/5", ">", "1/2",
     "\\sin\\left(\\dfrac{1}{5}\\pi\\right) - \\dfrac{1}{2} = \\int_0^{\\pi/2} "
     "\\frac{\\left(1 - \\sin{\\left(x \\right)}\\right) \\left(637 \\sin{\\left(x \\right)} "
     "+ 125\\right) \\sin{\\left(\\frac{3 x}{5} \\right)}}{750} \\mathrm{d} x > 0"),
    ("sin_pi_q", "1/7", "<", "4/9",
     "\\dfrac{4}{9} - \\sin\\left(\\dfrac{1}{7}\\pi\\right) = \\int_0^{\\pi/2} "
     "\\frac{\\left(1 - \\sin{\\left(x \\right)}\\right)^{2} \\cdot \\left(25676 - 6656 "
     "\\sin{\\left(x \\right)}\\right) \\sin{\\left(\\frac{5 x}{7} \\right)}}{151263} "
     "\\mathrm{d} x > 0"),
    ("cos_pi_q", "2/9", ">", "3/4",
     "\\cos\\left(\\dfrac{2}{9}\\pi\\right) - \\dfrac{3}{4} = \\int_0^{\\pi/2} "
     "\\frac{\\left(1 - \\sin{\\left(x \\right)}\\right) \\left(9793 - 4991 \\sin{\\left(x "
     "\\right)}\\right) \\sin{\\left(\\frac{4 x}{9} \\right)} \\sin{\\left(x \\right)}}"
     "{23328} \\mathrm{d} x > 0"),
    ("sin_q_degree", "18", "<", "1/3",
     "\\dfrac{1}{3} - \\sin\\left(18^\\circ\\right) = \\int_0^{\\pi/2} "
     "\\frac{\\left(1 - \\sin{\\left(x \\right)}\\right) \\left(25 - 14 \\sin{\\left(x "
     "\\right)}\\right) \\sin{\\left(\\frac{4 x}{5} \\right)}}{125} \\mathrm{d} x > 0"),
    ("cos_q_degree", "40", ">", "3/4",
     "\\cos\\left(40^\\circ\\right) - \\dfrac{3}{4} = \\int_0^{\\pi/2} "
     "\\frac{\\left(1 - \\sin{\\left(x \\right)}\\right) \\left(9793 - 4991 \\sin{\\left(x "
     "\\right)}\\right) \\sin{\\left(\\frac{4 x}{9} \\right)} \\sin{\\left(x \\right)}}"
     "{23328} \\mathrm{d} x > 0"),
]


@pytest.mark.parametrize("kind,power,comp,bound,expected", RENDER)
def test_render_verbatim(kind, power, comp, bound, expected):
    mod = trig_q if kind in ("sin_q", "cos_q", "tan_q", "cot_q") else trig_pi
    res = mod.prove(kind, F(power), comp, F(bound))
    eq = mod.render_equation(res["parameters"], kind, power, comp, bound)
    assert eq == expected
