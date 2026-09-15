"""exp + hyperbolic kernel tests.

Anchored on the article's printed expansions (types 2/4/5/6, 24-27 —
all use the symmetric basis t^n(1-t)^n) and on live-site answers
harvested into bench/data/golden.jsonl.
"""

from fractions import Fraction as F

import pytest
import sympy as sp
from mpmath import mp

from attention_calculator import engine
from attention_calculator.kernels import exp_family, hyperbolic
from attention_calculator.moment import add, scale


def moment_of(*terms):
    """sum_i coef_i * basis-moment — test helper mirroring the P expansion."""
    out = {}
    for coef, mom in terms:
        out = add(out, scale(mom, coef))
    return out


# ---------------------------------------------------------------- article expansions


@pytest.mark.parametrize(
    "n, ce, c1",
    [
        # ∫ x^n(1-x)^n (a+bx) e^x dx = (ce·a + ce·b)e + (c1·a + c1·b)
        (1, (-1, 3), (3, -8)),
        (2, (14, -64), (-38, 174)),  # 2(7a-32b)e + 2(-19a+87b)
        (3, (-426, 2790), (1158, -7584)),  # 6(-71a+465b)e + 6(193a-1264b)
    ],
)
def test_article_e_expansion(n, ce, c1):
    a, b = F(7), F(-3)
    mom = moment_of(
        (a, exp_family.basis_x_moment(n, n, 0, F(1), "e")),
        (b, exp_family.basis_x_moment(n, n, 1, F(1), "e")),
    )
    assert mom == {"e": F(ce[0]) * a + F(ce[1]) * b, "1": F(c1[0]) * a + F(c1[1]) * b}


@pytest.mark.parametrize("q", [F(5, 2), F(3)])
def test_article_e_q_expansion(q):
    """∫ x(1-x)(a+bx) e^{qx} dx =
    [(a+b)q² + (-2a-4b)q + 6b] e^q / q⁴ + [aq² + 2(a-b)q - 6b] / q⁴."""
    a, b = F(2), F(-5)
    mom = moment_of(
        (a, exp_family.basis_x_moment(1, 1, 0, q, "e_q")),
        (b, exp_family.basis_x_moment(1, 1, 1, q, "e_q")),
    )
    assert mom == {
        "e_q": ((a + b) * q**2 + (-2 * a - 4 * b) * q + 6 * b) / q**4,
        "1": (a * q**2 + 2 * (a - b) * q - 6 * b) / q**4,
    }


@pytest.mark.parametrize(
    "n, ce, c1",
    [
        (1, (F(1, 10), F(1, 10)), (F(9, 10), F(-7, 10))),
        (2, (F(7, 85), F(-15, 442)), (F(-109, 85), F(2421, 2210))),
    ],
)
def test_article_e_pi_expansion(n, ce, c1):
    a, b = F(7), F(-3)
    mom = moment_of(
        (a, exp_family.basis_sin_moment(n, n, 0)), (b, exp_family.basis_sin_moment(n, n, 1))
    )
    assert mom == {"e_pi": ce[0] * a + ce[1] * b, "1": c1[0] * a + c1[1] * b}


@pytest.mark.parametrize("q", [F(3, 2), F(2)])
def test_article_hyperbolic_expansion(q):
    """∫ x(1-x)(a+bx+cx²) sinh(qx) dx =
    [(a+b+c)q³ + (6b+18c)q] sinh q / q⁵
    + [-(2a+4b+6c)q² - 24c] cosh q / q⁵ + [(2a-2b)q² + 24c] / q⁵."""
    a, b, c = F(2), F(-5), F(3)
    mom = moment_of(
        (a, hyperbolic.basis_moment(1, 1, 0, q)),
        (b, hyperbolic.basis_moment(1, 1, 1, q)),
        (c, hyperbolic.basis_moment(1, 1, 2, q)),
    )
    assert mom == {
        "sinh_q": ((a + b + c) * q**3 + (6 * b + 18 * c) * q) / q**5,
        "cosh_q": (-(2 * a + 4 * b + 6 * c) * q**2 - 24 * c) / q**5,
        "1": ((2 * a - 2 * b) * q**2 + 24 * c) / q**5,
    }


# ---------------------------------------------------------------- moment <-> numeric integration

mp.dps = 40


def tomf(f):
    """Fraction -> mpf (mpmath won't take Fraction directly)."""
    return mp.mpf(f.numerator) / f.denominator


def mom_val(mom, **consts):
    """Evaluate a Moment dict given mpf values for each constant symbol."""
    return sum(tomf(c) * consts[k] for k, c in mom.items())


@pytest.mark.parametrize(
    "m,n,j,q",
    [
        (0, 0, 0, F(1)),
        (2, 3, 1, F(1)),
        (1, 2, 0, F(5, 2)),
        (3, 1, 1, F(1, 3)),
    ],
)
def test_basis_x_moment_numeric(m, n, j, q):
    qq = tomf(q)
    mom = exp_family.basis_x_moment(m, n, j, q, "e_q")
    val = mom_val(mom, e_q=mp.exp(qq), **{"1": 1})
    want = mp.quad(lambda x: x ** (m + j) * (1 - x) ** n * mp.exp(qq * x), [0, 1])
    assert mp.almosteq(val, want, rel_eps=mp.mpf("1e-30"))


@pytest.mark.parametrize("m,n,j", [(0, 0, 0), (1, 1, 0), (2, 1, 1), (0, 3, 1)])
def test_basis_sin_moment_numeric(m, n, j):
    mom = exp_family.basis_sin_moment(m, n, j)
    val = mom_val(mom, e_pi=mp.exp(mp.pi), **{"1": 1})
    want = mp.quad(lambda x: mp.sin(x) ** (m + j) * (1 - mp.sin(x)) ** n * mp.exp(x), [0, mp.pi])
    assert mp.almosteq(val, want, rel_eps=mp.mpf("1e-30"))


@pytest.mark.parametrize(
    "m,n,j,q",
    [
        (0, 0, 0, F(1)),
        (1, 1, 2, F(1)),
        (2, 0, 1, F(3, 2)),
        (0, 2, 0, F(2, 3)),
    ],
)
def test_basis_sinh_moment_numeric(m, n, j, q):
    qq = tomf(q)
    mom = hyperbolic.basis_moment(m, n, j, q)
    val = mom_val(mom, sinh_q=mp.sinh(qq), cosh_q=mp.cosh(qq), **{"1": 1})
    want = mp.quad(lambda x: x ** (m + j) * (1 - x) ** n * mp.sinh(qq * x), [0, 1])
    assert mp.almosteq(val, want, rel_eps=mp.mpf("1e-30"))


def test_moments_match_sympy_integrate():
    """Small-n cross-validation against live sympy.integrate."""
    x, q = sp.symbols("x q", real=True)
    expr = sp.integrate(x * (1 - x) * sp.exp(q * x), (x, 0, 1))
    got = exp_family.basis_x_moment(1, 1, 0, F(2), "e_q")
    # note: at q=2 the e^q coefficient of this basis moment is exactly 0 —
    # moments may omit zero keys
    assert sp.simplify(expr.subs(q, 2) - (got.get("e_q", F(0)) * sp.E**2 + got.get("1", F(0)))) == 0


# ---------------------------------------------------------------- prove end-to-end (site-verified)


@pytest.mark.parametrize(
    "kind,power,comp,bound,m,n,solution",
    [
        ("e", "1", ">", "8/3", 1, 1, "a = 0, b = 1/3"),
        ("e", "1", "<", "25/9", 1, 1, "a = 1/3, b = -2/9"),
        ("e", "1", ">", "2", 0, 0, "a = 1, b = -1"),
        ("e", "1", "<", "3", 0, 1, "a = 0, b = 1"),
        ("e", "1", ">", "27/10", 1, 2, "a = 1/10, b = 1/10"),
        ("e", "3/4", ">", "2", 1, 1, "a = 0, b = 1/4"),
        ("e_q", "2", "<", "15/2", 2, 2, "a = 0, b = 2"),
        ("e_q", "1/2", ">", "3/2", 0, 0, "a = 1/4, b = -1/4"),
        ("e_q", "3", ">", "20", 2, 4, "a = 33/8, b = -27/8"),
        ("e_q", "1/3", "<", "2", 0, 0, "a = 1/3, b = 1/3"),
        ("e_pi", "1", ">", "23", 2, 5, "a = 182650/112211, b = -422240/336633"),
        ("e_pi", "3/4", ">", "21/2", 1, 2, "a = 165/224, b = 255/28"),
        ("sinh_q", "1", ">", "7/6", 0, 1, "a = 1/6, b = -1/3, c= 1/6"),
        ("cosh_q", "1", "<", "8/5", 0, 1, "a = 0, b = 7/10, c= -1/10"),
        ("tanh_q", "1", ">", "3/4", 0, 1, "a = 0, b = 1/8, c= 1/8"),
        ("coth_q", "1", "<", "2", 0, 0, "a = 2, b = 0, c= -1"),
        ("cosh_q", "2/3", ">", "6/5", 0, 0, "a = 1/5, b = -4/15, c= 2/15"),
        ("tanh_q", "2/3", "<", "3/5", 0, 1, "a = 1/5, b = 0, c= -2/45"),
        ("coth_q", "2/3", ">", "5/3", 0, 1, "a = 1/3, b = 0, c= -2/27"),
    ],
)
def test_prove_site_parity(kind, power, comp, bound, m, n, solution):
    fam = exp_family if kind in ("e", "e_q", "e_pi") else hyperbolic
    res = fam.prove(kind, F(power), comp, F(bound))
    p = res["parameters"]
    assert (p["m"], p["n"]) == (m, n)
    assert res["solution"] == solution
    # site wire shape: 9 parameter fields + solution string
    assert set(p) == {"m", "n", "a_val", "b_val", "c_val", "au_val", "bu_val", "cu_val", "u_val"}
    assert F(p["a_val"]) * p["u_val"] == p["au_val"]
    assert F(p["b_val"]) * p["u_val"] == p["bu_val"]
    assert F(p["c_val"]) * p["u_val"] == p["cu_val"]


# ---------------------------------------------------------------- render_equation byte-parity


@pytest.mark.parametrize(
    "kind,power,comp,bound,equation",
    [
        (
            "e",
            "1",
            ">",
            "8/3",
            r"e - \dfrac{8}{3} = \int_0^1 "
            r"\frac{x^{2} \cdot \left(1 - x\right) e^{x}}{3} \mathrm{d} x > 0",
        ),
        ("e", "1", ">", "2", r"e - 2 = \int_0^1 \left(1 - x\right) e^{x} \mathrm{d} x > 0"),
        (
            "e",
            "1",
            ">",
            "27/10",
            r"e - \dfrac{27}{10} = \int_0^1 \frac{x \left(1 - x\right)^{2} "
            r"\left(x + 1\right) e^{x}}{10} \mathrm{d} x > 0",
        ),
        (
            "e",
            "1",
            "<",
            "68/25",
            r"\dfrac{68}{25} - e = \int_0^1 \frac{x^{2} \left(1 - x\right)^{2} "
            r"\left(x + 1\right) e^{x}}{50} \mathrm{d} x > 0",
        ),
        (
            "e",
            "3/4",
            ">",
            "2",
            r"\dfrac{3}{4}e - 2 = \int_0^1 "
            r"\frac{x^{2} \cdot \left(1 - x\right) e^{x}}{4} \mathrm{d} x > 0",
        ),
        (
            "e_q",
            "2",
            "<",
            "15/2",
            r"\dfrac{15}{2} - e^2 = \int_0^1 2 x^{3} "
            r"\left(1 - x\right)^{2} e^{2 x} \mathrm{d} x > 0",
        ),
        (
            "e_q",
            "1/2",
            ">",
            "3/2",
            r"e^\dfrac{1}{2} - \dfrac{3}{2} = \int_0^1 "
            r"\frac{\left(1 - x\right) e^{\frac{x}{2}}}{4} \mathrm{d} x > 0",
        ),
        (
            "e_q",
            "3",
            ">",
            "20",
            r"e^3 - 20 = \int_0^1 \frac{x^{2} \left(1 - x\right)^{4} "
            r"\cdot \left(33 - 27 x\right) e^{3 x}}{8} \mathrm{d} x > 0",
        ),
        (
            "e_pi",
            "1",
            ">",
            "23",
            r"e^{\pi} - 23 = \int_0^{\pi} \frac{\left(1 - "
            r"\sin{\left(x \right)}\right)^{5} \cdot \left(547950 - "
            r"422240 \sin{\left(x \right)}\right) e^{x} "
            r"\sin^{2}{\left(x \right)}}{336633} \mathrm{d} x > 0",
        ),
        (
            "sinh_q",
            "1",
            ">",
            "7/6",
            r"\sinh1 - \dfrac{7}{6} = \int_0^1 "
            r"\frac{\left(1 - x\right) \left(x^{2} - 2 x + 1\right) "
            r"\sinh{\left(x \right)}}{6} \mathrm{d} x > 0",
        ),
        (
            "cosh_q",
            "2/3",
            ">",
            "6/5",
            r"\cosh\dfrac{2}{3} - \dfrac{6}{5} = \int_0^1 "
            r"\frac{\left(2 x^{2} - 4 x + 3\right) "
            r"\sinh{\left(\frac{2 x}{3} \right)}}{15} \mathrm{d} x > 0",
        ),
        (
            "tanh_q",
            "1",
            ">",
            "3/4",
            r"\tanh1 - \dfrac{3}{4} = \int_0^1 \dfrac{1}{\cosh(1)} "
            r"\frac{\left(1 - x\right) \left(x^{2} + x\right) "
            r"\sinh{\left(x \right)}}{8} \mathrm{d} x > 0",
        ),
        (
            "tanh_q",
            "2/3",
            "<",
            "3/5",
            r"\dfrac{3}{5} - \tanh\dfrac{2}{3} = \int_0^1 "
            r"\dfrac{1}{\cosh(2/3)} \frac{\left(1 - x\right) \left(9 - 2 x^{2}\right) "
            r"\sinh{\left(\frac{2 x}{3} \right)}}{45} \mathrm{d} x > 0",
        ),
        (
            "coth_q",
            "1",
            "<",
            "2",
            r"2 - \coth1 = \int_0^1 \dfrac{1}{\sinh(1)} "
            r"\left(2 - x^{2}\right) \sinh{\left(x \right)} \mathrm{d} x > 0",
        ),
    ],
)
def test_render_equation_bytes(kind, power, comp, bound, equation):
    fam = exp_family if kind in ("e", "e_q", "e_pi") else hyperbolic
    res = fam.prove(kind, F(power), comp, F(bound))
    assert fam.render_equation(res["parameters"], kind, F(power), comp, F(bound)) == equation


# ---------------------------------------------------------------- error paths


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("e", "1", "<", "8/3"),  # e > 8/3 actually
        ("e", "1", ">", "3"),  # e < 3
        ("e_q", "2", "<", "7"),  # e² ≈ 7.389
        ("e_pi", "1", "<", "23"),  # e^π ≈ 23.14
        ("sinh_q", "1", "<", "1"),  # sinh 1 ≈ 1.175
        ("tanh_q", "1", ">", "1"),  # tanh 1 ≈ 0.762
        ("cosh_q", "1", ">", "2"),  # cosh 1 ≈ 1.543
    ],
)
def test_wrong_direction(kind, power, comp, bound):
    fam = exp_family if kind in ("e", "e_q", "e_pi") else hyperbolic
    with pytest.raises(engine.WrongDirection):
        fam.prove(kind, F(power), comp, F(bound))


def test_no_solution_beyond_cap():
    """e² bound tight enough to need (m,n) = (11,11) — site reports
    '在指数不超过10的范围内未找到>方向的解' (probed)."""
    with pytest.raises(engine.NoSolution):
        exp_family.prove("e_q", F(2), ">", F(402325288800, 54448806913))


def test_zero_power_edge_cases():
    """power=0 parity: site answers 方向反了 for e/e_pi/sinh/tanh
    (degenerate solve) and a 500 for e_q/cosh/coth (1/q crash)."""
    with pytest.raises(engine.WrongDirection):
        exp_family.prove("e", F(0), ">", F(1, 2))
    with pytest.raises(engine.WrongDirection):
        exp_family.prove("e_pi", F(0), ">", F(1, 2))
    with pytest.raises(engine.WrongDirection):
        hyperbolic.prove("sinh_q", F(0), ">", F(1, 2))
    with pytest.raises(engine.WrongDirection):
        hyperbolic.prove("tanh_q", F(0), ">", F(1, 2))
    with pytest.raises(ZeroDivisionError):
        exp_family.prove("e_q", F(0), ">", F(1, 2))
    with pytest.raises(ZeroDivisionError):
        hyperbolic.prove("cosh_q", F(0), ">", F(1, 2))
    with pytest.raises(ZeroDivisionError):
        hyperbolic.prove("coth_q", F(0), ">", F(1, 2))


def test_zero_bound_accepted():
    """The site accepts bound=0 (e>0 -> m=n=0, P=1+x)."""
    res = exp_family.prove("e", F(1), ">", F(0))
    assert res["parameters"]["m"] == 0 and res["parameters"]["n"] == 0
    assert res["solution"] == "a = 1, b = 1"
