"""quadlog kernel family tests.

Ground truth comes from three independent sources, compared exactly:
1. the moment expansions printed in the author's article (types 1, 3, 12, 13,
   31, 40) — hand-listed below as Fraction test vectors;
2. live zhuyidao.net responses probed during development (parameters and
   rendered equations, byte-identical);
3. sympy.integrate (exact, r=0 cases) and mpmath quadrature (numeric, r>=1).
"""

from fractions import Fraction

import mpmath as mp
import pytest
import sympy as sp

from attention_calculator import solve
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.kernels import quadlog

F = Fraction
mp.mp.dps = 50


def even_term(k, r):
    """Moment of x^{2k} ln^r(1/x)/(1+x^2) as (const_coef, rat)."""
    return quadlog.ln_moment(k, r, False)


def odd_term(k, r):
    """Moment of x^{2k+1} ln^r(1/x)/(1+x^2) as (const_coef, rat)."""
    return quadlog.ln_moment(k, r, True)


def basis(m, n, r, odd, factor, sym="C"):
    """Moment pair (for a, b) of x^{2m+odd}(1-x^2)^n·{1,x²}·ln^r/(1+x²).

    ``factor`` rescales only the constant coefficient (e.g. eta(3) = 3/4 zeta(3));
    the rational part is untouched, mirroring quadlog.prove's term closure.
    """

    def term(k):
        cc, rat = quadlog.ln_moment(k, r, odd)
        return cc * factor, rat

    return quadlog.basis_moments(m, n, odd, sym, term)


# ---------------------------------------------------------------- article moments

# Type 31 (Catalan): printed expansions for x^{2n}(1-x^2)^n, n=1,2,3.
CATALAN_VECTORS = [
    (1, 1, {"catalan": F(-2), "1": F(17, 9)}, {"catalan": F(2), "1": F(-409, 225)}),
    (2, 2, {"catalan": F(4), "1": F(-40298, 11025)}, {"catalan": F(-4), "1": F(363826, 99225)}),
    (
        3,
        3,
        {"catalan": F(-8), "1": F(12571156, 1715175)},
        {"catalan": F(8), "1": F(-14867117548, 2029052025)},
    ),
]


@pytest.mark.parametrize("m,n,ea,eb", CATALAN_VECTORS)
def test_catalan_moments_article(m, n, ea, eb):
    a_mom, b_mom = basis(m, n, r=1, odd=False, factor=F(1), sym="catalan")
    assert a_mom == ea and b_mom == eb


# Type 40 (zeta odd): printed expansions, site exponent is x^{2m+1}.
ZETA_VECTORS = [
    # (m, n, r, factor, expected a-moment, expected b-moment) — r=2 -> zeta(3)
    (0, 0, 2, F(3, 4), {"z3": F(3, 16)}, {"z3": F(-3, 16), "1": F(1, 4)}),
    (1, 0, 2, F(3, 4), {"z3": F(-3, 16), "1": F(1, 4)}, {"z3": F(3, 16), "1": F(-7, 32)}),
    (0, 1, 2, F(3, 4), {"z3": F(3, 8), "1": F(-1, 4)}, {"z3": F(-3, 8), "1": F(15, 32)}),
    (1, 1, 2, F(3, 4), {"z3": F(-3, 8), "1": F(15, 32)}, {"z3": F(3, 8), "1": F(-193, 432)}),
    # r=4 -> eta(5) = 15/16 zeta(5); article's zeta(5) table
    (0, 0, 4, F(15, 16), {"z3": F(45, 64)}, {"z3": F(-45, 64), "1": F(3, 4)}),
    (1, 0, 4, F(15, 16), {"z3": F(-45, 64), "1": F(3, 4)}, {"z3": F(45, 64), "1": F(-93, 128)}),
    (0, 1, 4, F(15, 16), {"z3": F(45, 32), "1": F(-3, 4)}, {"z3": F(-45, 32), "1": F(189, 128)}),
    (
        1,
        1,
        4,
        F(15, 16),
        {"z3": F(-45, 32), "1": F(189, 128)},
        {"z3": F(45, 32), "1": F(-7549, 5184)},
    ),
]


@pytest.mark.parametrize("m,n,r,factor,ea,eb", ZETA_VECTORS)
def test_zeta_moments_article(m, n, r, factor, ea, eb):
    a_mom, b_mom = basis(m, n, r=r, odd=True, factor=factor, sym="z3")
    assert a_mom == ea and b_mom == eb


# Type 3 (pi^k): article prints x^m(a+bx^2)ln^{k-1}(1/x)/(1+x^2) for the
# site basis x^{2m(+1)}(1-x^2)^0 — i.e. the m_lit exponent, n=0.
PIN_VECTORS = [
    # k, m_lit, expected (a-moment, b-moment) with sym "pi"
    (2, 1, ({"pi": F(1, 48)}, {"pi": F(-1, 48), "1": F(1, 4)})),
    (2, 3, ({"pi": F(-1, 48), "1": F(1, 4)}, {"pi": F(1, 48), "1": F(-3, 16)})),
    (3, 0, ({"pi": F(1, 16)}, {"pi": F(-1, 16), "1": F(2)})),
    (3, 2, ({"pi": F(-1, 16), "1": F(2)}, {"pi": F(1, 16), "1": F(-52, 27)})),
    (4, 1, ({"pi": F(7, 1920)}, {"pi": F(-7, 1920), "1": F(3, 8)})),
    (4, 3, ({"pi": F(-7, 1920), "1": F(3, 8)}, {"pi": F(7, 1920), "1": F(-45, 128)})),
    (5, 0, ({"pi": F(5, 64)}, {"pi": F(-5, 64), "1": F(24)})),
    (5, 2, ({"pi": F(-5, 64), "1": F(24)}, {"pi": F(5, 64), "1": F(-1936, 81)})),
]


@pytest.mark.parametrize("k,m_lit,expected", PIN_VECTORS)
def test_pin_moments_article(k, m_lit, expected):
    odd = k % 2 == 0
    factor = quadlog.ETA_PI[k] if odd else quadlog.BETA_PI[k]
    m = (m_lit - 1) // 2 if odd else m_lit // 2
    a_mom, b_mom = basis(m, 0, r=k - 1, odd=odd, factor=factor, sym="pi")
    assert (a_mom, b_mom) == expected


def test_beta_eta_pi_tables():
    """The hardcoded beta/eta rational tables match sympy-derived values."""
    for j in range(5):  # beta(2j+1)/pi^{2j+1} = (-1)^j E_{2j}/(4^{j+1}(2j)!)
        k = 2 * j + 1
        v = F((-1) ** j * sp.euler(2 * j).p, sp.factorial(2 * j) * 4 ** (j + 1))
        assert quadlog.BETA_PI[k] == v
    for j in range(1, 6):  # eta(2j)/pi^{2j} = (1-2^{1-2j})·(-1)^{j+1} B_{2j} 4^j/(2(2j)!)
        k = 2 * j
        b = sp.bernoulli(2 * j)
        zeta = F((-1) ** (j + 1) * b.p * 4**j, b.q * 2 * sp.factorial(2 * j))
        assert quadlog.ETA_PI[k] == (1 - F(2) ** (1 - k)) * zeta


# ------------------------------------------------------- cross-checks vs sympy/mpmath


def test_even_moments_vs_sympy():
    """r=0 closed forms: sympy.integrate handles 1/(1+x^2) exactly."""
    x = sp.symbols("x")
    for k in range(4):
        v = sp.integrate(x ** (2 * k) / (1 + x**2), (x, 0, 1))
        cc, rat = even_term(k, 0)
        got = cc * sp.pi / 4 + rat
        assert sp.simplify(v - got) == 0


def test_atan_moments_vs_sympy():
    x = sp.symbols("x")
    for q in (sp.Rational(2), sp.Rational(1, 3), sp.Rational(3, 2)):
        for k in range(3):
            v = sp.integrate(x ** (2 * k) / (1 + q**2 * x**2), (x, 0, 1))
            cc, rat = quadlog.atan_moment(k, F(q))
            got = cc * sp.atan(q) + rat
            assert sp.simplify(v - got) == 0


def test_ln_moments_vs_mpmath():
    """Numeric quadrature of the raw integrals vs closed-form moments (n<=3)."""

    def eta(s):
        return (1 - mp.mpf(2) ** (1 - s)) * mp.zeta(s)

    def quad(p, r):
        return mp.quad(lambda x: x**p * mp.log(1 / x) ** r / (1 + x**2), [0, 1])

    # const is the value of beta(r+1) (even) resp. eta(r+1) (odd)
    for odd, r, const in [
        (False, 1, mp.catalan),
        (True, 2, eta(3)),
        (True, 1, eta(2)),
        (False, 2, mp.pi**3 / 32),
    ]:
        term = odd_term if odd else even_term
        for k in range(4):
            cc, rat = term(k, r)
            with mp.workdps(50):
                want = quad(2 * k + odd, r)
                got = mp.mpf(str(cc)) * const + mp.mpf(str(rat))
            assert abs(got - want) < mp.mpf("1e-30")


def test_basis_moments_vs_mpmath():
    """Full (1-x^2)^n-expanded moments vs numeric quadrature, n<=3."""
    for m, n in [(0, 0), (1, 1), (0, 2), (2, 1), (1, 3), (3, 3)]:
        a_mom, b_mom = basis(m, n, r=1, odd=False, factor=F(1))
        for c, mom in enumerate([a_mom, b_mom]):
            with mp.workdps(50):
                want = mp.quad(
                    lambda x, m=m, n=n, c=c: (
                        x ** (2 * m) * (1 - x**2) ** n * x ** (2 * c) * mp.log(1 / x) / (1 + x**2)
                    ),
                    [0, 1],
                )
                got = mp.mpf(str(mom.get("1", 0))) + mp.mpf(str(mom.get("C", 0))) * mp.catalan
            assert abs(got - want) < mp.mpf("1e-28")


# -------------------------------------------------------------- end-to-end vs site

# (kind, power, comp, bound, expected parameters) — live-probed site answers.
SITE_CASES = [
    (
        "pi",
        "1",
        "<",
        "22/7",
        dict(
            m=3,
            n=3,
            a_val="47/120",
            b_val="-13/120",
            au_val="47",
            bu_val="-13",
            cu_val="0",
            u_val="120",
        ),
    ),
    (
        "pi",
        "8",
        ">",
        "25",
        dict(m=2, n=2, a_val="61/8", b_val="-3/8", au_val="61", bu_val="-3", cu_val="0", u_val="8"),
    ),
    (
        "pi",
        "1",
        "<",
        "355/113",
        dict(
            m=5,
            n=8,
            a_val="260341/25919488",
            b_val="-144651/25919488",
            au_val="260341",
            bu_val="-144651",
            cu_val="0",
            u_val="25919488",
        ),
    ),
    (
        "pi",
        "1",
        "<",
        "4",
        dict(m=0, n=0, a_val="0", b_val="4", au_val="0", bu_val="4", cu_val="0", u_val="1"),
    ),
    (
        "pi",
        "1",
        ">",
        "3",
        dict(m=1, n=1, a_val="1/2", b_val="5/2", au_val="1", bu_val="5", cu_val="0", u_val="2"),
    ),
    (
        "pi",
        "1",
        "<",
        "16/5",
        dict(m=1, n=1, a_val="1", b_val="-1", au_val="1", bu_val="-1", cu_val="0", u_val="1"),
    ),
    (
        "pi",
        "1",
        ">",
        "31/10",
        dict(
            m=1, n=2, a_val="5/16", b_val="21/16", au_val="5", bu_val="21", cu_val="0", u_val="16"
        ),
    ),
    (
        "pi",
        "1",
        "<",
        "104348/33215",
        dict(
            m=9,
            n=10,
            a_val="2501801/1219952640",
            b_val="-2263639/1219952640",
            au_val="2501801",
            bu_val="-2263639",
            cu_val="0",
            u_val="1219952640",
        ),
    ),
    ("pi", "1", ">", "3141592653589793/1000000000000000", dict(m=14, n=18)),
    (
        "pi",
        "3/2",
        "<",
        "5",
        dict(m=0, n=1, a_val="0", b_val="3", au_val="0", bu_val="3", cu_val="0", u_val="1"),
    ),
    (
        "pi_n",
        "2",
        "<",
        "10",
        dict(
            m=1,
            n=2,
            a_val="108/13",
            b_val="-48/13",
            au_val="108",
            bu_val="-48",
            cu_val="0",
            u_val="13",
        ),
    ),
    (
        "pi_n",
        "2",
        ">",
        "49/5",
        dict(
            m=1,
            n=2,
            a_val="36/65",
            b_val="816/65",
            au_val="36",
            bu_val="816",
            cu_val="0",
            u_val="65",
        ),
    ),
    (
        "pi_n",
        "3",
        ">",
        "31",
        dict(
            m=2,
            n=3,
            a_val="338003449/217788864",
            b_val="-97574279/217788864",
            au_val="338003449",
            bu_val="-97574279",
            cu_val="0",
            u_val="217788864",
        ),
    ),
    (
        "pi_n",
        "4",
        ">",
        "97",
        dict(
            m=1,
            n=1,
            a_val="1536/455",
            b_val="63936/455",
            au_val="1536",
            bu_val="63936",
            cu_val="0",
            u_val="455",
        ),
    ),
    (
        "pi_n",
        "1",
        "<",
        "22/7",
        dict(
            m=3,
            n=3,
            a_val="47/120",
            b_val="-13/120",
            au_val="47",
            bu_val="-13",
            cu_val="0",
            u_val="120",
        ),
    ),
    (
        "pi_n",
        "5",
        ">",
        "306",
        dict(
            m=1,
            n=2,
            a_val="52763017/902260960",
            b_val="2939998089/902260960",
            au_val="52763017",
            bu_val="2939998089",
            cu_val="0",
            u_val="902260960",
        ),
    ),
    (
        "pi_n",
        "6",
        "<",
        "962",
        dict(
            m=2,
            n=0,
            a_val="25308/155",
            b_val="105948/155",
            au_val="25308",
            bu_val="105948",
            cu_val="0",
            u_val="155",
        ),
    ),
    (
        "pi_n",
        "10",
        ">",
        "93648",
        dict(
            m=3,
            n=0,
            a_val="5827661824/50290065",
            b_val="19130220544/50290065",
            au_val="5827661824",
            bu_val="19130220544",
            cu_val="0",
            u_val="50290065",
        ),
    ),
    (
        "pi_n",
        "5/2",
        ">",
        "17",
        dict(
            m=1,
            n=0,
            a_val="6859/40",
            b_val="7371/40",
            au_val="6859",
            bu_val="7371",
            cu_val="0",
            u_val="40",
        ),
    ),
    (
        "pi_n",
        "5/2",
        "<",
        "35/2",
        dict(
            m=2,
            n=0,
            a_val="1092961/38880",
            b_val="318125/7776",
            au_val="1092961",
            bu_val="1590625",
            cu_val="0",
            u_val="38880",
        ),
    ),
    (
        "pi_n",
        "1/3",
        "<",
        "3/2",
        dict(
            m=0, n=1, a_val="1/16", b_val="33/16", au_val="1", bu_val="33", cu_val="0", u_val="16"
        ),
    ),
    (
        "arctan_q",
        "2",
        "<",
        "6/5",
        dict(
            m=1,
            n=3,
            a_val="1871/500",
            b_val="-177/125",
            au_val="1871",
            bu_val="-708",
            cu_val="0",
            u_val="500",
        ),
    ),
    (
        "arctan_q",
        "1/2",
        "<",
        "1/2",
        dict(m=0, n=0, a_val="0", b_val="1/8", au_val="0", bu_val="1", cu_val="0", u_val="8"),
    ),
    (
        "arctan_q",
        "1/3",
        "<",
        "28/87",
        dict(
            m=1,
            n=1,
            a_val="13/15660",
            b_val="-1/3132",
            au_val="13",
            bu_val="-5",
            cu_val="0",
            u_val="15660",
        ),
    ),
    (
        "arctan_q",
        "2/3",
        "<",
        "3/5",
        dict(
            m=0, n=1, a_val="1/390", b_val="6/65", au_val="1", bu_val="36", cu_val="0", u_val="390"
        ),
    ),
    (
        "arctan_q",
        "2/3",
        "<",
        "59/100",
        dict(
            m=1,
            n=1,
            a_val="407/14040",
            b_val="-97/3510",
            au_val="407",
            bu_val="-388",
            cu_val="0",
            u_val="14040",
        ),
    ),
    (
        "arctan_q",
        "1/3",
        ">",
        "31/97",
        dict(
            m=0,
            n=1,
            a_val="23/5820",
            b_val="-19/5820",
            au_val="23",
            bu_val="-19",
            cu_val="0",
            u_val="5820",
        ),
    ),
    (
        "arctan_q",
        "1",
        "<",
        "4/5",
        dict(m=1, n=1, a_val="1/4", b_val="-1/4", au_val="1", bu_val="-1", cu_val="0", u_val="4"),
    ),
    (
        "arccot_q",
        "1/2",
        "<",
        "6/5",
        dict(
            m=1,
            n=3,
            a_val="1871/500",
            b_val="-177/125",
            au_val="1871",
            bu_val="-708",
            cu_val="0",
            u_val="500",
        ),
    ),
    (
        "arccot_q",
        "2",
        ">",
        "23/50",
        dict(m=1, n=1, a_val="1/40", b_val="1/80", au_val="2", bu_val="1", cu_val="0", u_val="80"),
    ),
    (
        "catalan",
        "1",
        "<",
        "109/119",
        dict(
            m=5,
            n=5,
            a_val="733466923/43978014720",
            b_val="-640846037/43978014720",
            au_val="733466923",
            bu_val="-640846037",
            cu_val="0",
            u_val="43978014720",
        ),
    ),
    (
        "catalan",
        "2",
        "<",
        "2",
        dict(m=0, n=0, a_val="0", b_val="2", au_val="0", bu_val="2", cu_val="0", u_val="1"),
    ),
    (
        "zeta3",
        "1",
        ">",
        "119/99",
        dict(
            m=3,
            n=4,
            a_val="35045/3681689",
            b_val="3786824/11045067",
            au_val="105135",
            bu_val="3786824",
            cu_val="0",
            u_val="11045067",
        ),
    ),
    (
        "zeta3",
        "1",
        ">",
        "1",
        dict(m=1, n=0, a_val="16/3", b_val="32/3", au_val="16", bu_val="32", cu_val="0", u_val="3"),
    ),
    (
        "zeta3",
        "2",
        ">",
        "2",
        dict(m=1, n=0, a_val="32/3", b_val="64/3", au_val="32", bu_val="64", cu_val="0", u_val="3"),
    ),
]


@pytest.mark.parametrize(
    "kind,power,comp,bound,expected",
    SITE_CASES,
    ids=[f"{c[0]}:{c[1]}{c[2]}{c[3]}" for c in SITE_CASES],
)
def test_prove_matches_site(kind, power, comp, bound, expected):
    got = quadlog.prove(kind, F(power), comp, F(bound))["parameters"]
    for key, want in expected.items():
        assert str(got[key]) == str(want), f"{key}: {got[key]} != {want}"


def test_prove_via_dispatcher():
    """solve.prove wiring returns the same parameter set."""
    r = solve.prove("pi", "1", "<", "22/7")
    assert r["parameters"]["m"] == 3 and r["parameters"]["a_val"] == "47/120"
    assert r["solution"] == "a = 47/120, b = -13/120"


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("pi", "1", ">", "22/7"),  # pi < 22/7 is the true direction
        ("pi_n", "3", "<", "30"),  # pi^3 > 30 is the true direction
        ("catalan", "1", "<", "9/10"),  # C > 9/10
        ("zeta3", "1", "<", "1"),  # zeta(3) > 1
        # the kernel search itself finds a <=0 solution (the '>' proof negated):
        ("pi", "1", "<", "3141592653589793/1000000000000000"),
        # no <=0 solution within the limit; the dispatcher's float direction
        # check turns the false claim into WrongDirection (site-verified):
        ("arctan_q", "3", ">", "5/4"),  # arctan 3 < 5/4
        ("catalan", "2", "<", "1"),  # 2C > 1; needs the coefficient
    ],
)
def test_wrong_direction(kind, power, comp, bound):
    """Site 方向反了 cases; asserted through solve.prove, the site contract."""
    with pytest.raises(WrongDirection):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("arctan_q", "3", "<", "5/4"),  # true but needs exponents > 10
        ("pi_n", "2", "<", "986960440109/100000000000"),
    ],
)
def test_no_solution(kind, power, comp, bound):
    """Site 未找到 cases — NoSolution at both the kernel and dispatcher level."""
    with pytest.raises(NoSolution):
        quadlog.prove(kind, F(power), comp, F(bound))
    with pytest.raises(NoSolution):
        solve.prove(kind, power, comp, bound)


def test_pin_power_out_of_range():
    # the site has no range guard: out-of-range exponents die on the moment
    # table lookup -> KeyError -> the generic 500 (probe: pi_n 11>1 -> 500)
    with pytest.raises(KeyError):
        quadlog.prove("pi_n", F(11), ">", F(1))
    with pytest.raises(KeyError):
        quadlog.prove("pi_n", F(11, 2), ">", F(1))


# ------------------------------------------------------------------- rendering

RENDER_CASES = [
    (
        "pi",
        dict(m=3, n=3, au_val="47", bu_val="-13", u_val="120"),
        "<",
        "22/7",
        "1",
        "\\dfrac{22}{7} - \\pi = \\int_0^1 \\frac{x^{6} \\left(1 - x^{2}\\right)^{3} "
        "\\cdot \\left(47 - 13 x^{2}\\right)}{120 x^{2} + 120} \\mathrm{d} x > 0",
    ),
    (
        "pi",
        dict(m=2, n=2, au_val="61", bu_val="-3", u_val="8"),
        ">",
        "25",
        "8",
        "8\\pi - 25 = \\int_0^1 \\frac{x^{4} \\left(1 - x^{2}\\right)^{2} "
        "\\cdot \\left(61 - 3 x^{2}\\right)}{8 x^{2} + 8}  \\mathrm{d} x > 0",
    ),
    (
        "pi",
        dict(m=0, n=0, au_val="0", bu_val="4", u_val="1"),
        "<",
        "4",
        "1",
        "4 - \\pi = \\int_0^1 \\frac{4 x^{2}}{x^{2} + 1} \\mathrm{d} x > 0",
    ),
    (
        "pi",
        dict(m=1, n=1, au_val="1", bu_val="5", u_val="2"),
        ">",
        "3",
        "1",
        "\\pi - 3 = \\int_0^1 \\frac{x^{2} \\cdot \\left(1 - x^{2}\\right) "
        "\\left(5 x^{2} + 1\\right)}{2 x^{2} + 2}  \\mathrm{d} x > 0",
    ),
    (
        "pi",
        dict(m=1, n=1, au_val="1", bu_val="-1", u_val="1"),
        "<",
        "16/5",
        "1",
        "\\dfrac{16}{5} - \\pi = \\int_0^1 \\frac{x^{2} \\left(1 - x^{2}\\right)^{2}}"
        "{x^{2} + 1} \\mathrm{d} x > 0",
    ),
    (
        "pi",
        dict(m=0, n=1, au_val="0", bu_val="3", u_val="1"),
        "<",
        "5",
        "3/2",
        "5 - \\dfrac{3}{2}\\pi = \\int_0^1 \\frac{3 x^{2} \\cdot \\left(1 - x^{2}\\right)}"
        "{x^{2} + 1} \\mathrm{d} x > 0",
    ),
    (
        "pi_n",
        dict(m=1, n=2, au_val="108", bu_val="-48", u_val="13"),
        "<",
        "10",
        "2",
        "10 - \\pi^2 = \\int_0^1 \\frac{x^{3} \\left(1 - x^{2}\\right)^{2} "
        "\\cdot \\left(108 - 48 x^{2}\\right)}{13 x^{2} + 13}(\\ln(1/x)) \\mathrm{d} x > 0",
    ),
    (
        "pi_n",
        dict(m=2, n=3, au_val="338003449", bu_val="-97574279", u_val="217788864"),
        ">",
        "31",
        "3",
        "\\pi^3 - 31 = \\int_0^1 \\frac{x^{4} \\left(1 - x^{2}\\right)^{3} "
        "\\cdot \\left(338003449 - 97574279 x^{2}\\right)}{217788864 x^{2} + 217788864}"
        "(\\ln(1/x))^2 \\mathrm{d} x > 0",
    ),
    (
        "pi_n",
        dict(m=1, n=1, au_val="1536", bu_val="63936", u_val="455"),
        ">",
        "97",
        "4",
        "\\pi^4 - 97 = \\int_0^1 \\frac{x^{3} \\cdot \\left(1 - x^{2}\\right) "
        "\\left(63936 x^{2} + 1536\\right)}{455 x^{2} + 455}(\\ln(1/x))^3 \\mathrm{d} x > 0",
    ),
    (
        "pi_n",
        dict(m=3, n=0, au_val="5827661824", bu_val="19130220544", u_val="50290065"),
        ">",
        "93648",
        "10",
        "\\pi^{10} - 93648 = \\int_0^1 \\frac{x^{7} \\cdot \\left(19130220544 x^{2} + "
        "5827661824\\right)}{50290065 x^{2} + 50290065}(\\ln(1/x))^9 \\mathrm{d} x > 0",
    ),
    (
        "pi_n",
        dict(m=1, n=0, au_val="6859", bu_val="7371", u_val="40"),
        ">",
        "17",
        "5/2",
        "\\left(\\pi^\\dfrac{5}{2}\\right)^{2}- \\left(17\\right)^{2}= \\int_0^1 "
        "\\frac{x^{2} \\cdot \\left(7371 x^{2} + 6859\\right)}{40 x^{2} + 40}"
        "(\\ln(1/x))^4 \\mathrm{d} x > 0",
    ),
    (
        "pi_n",
        dict(m=2, n=0, au_val="1092961", bu_val="1590625", u_val="38880"),
        "<",
        "35/2",
        "5/2",
        "\\left(\\dfrac{35}{2}\\right)^{2} - \\left(\\pi^\\dfrac{5}{2}\\right)^{2} = "
        "\\int_0^1 \\frac{x^{4} \\cdot \\left(1590625 x^{2} + 1092961\\right)}"
        "{38880 x^{2} + 38880}(\\ln(1/x))^4 \\mathrm{d} x > 0",
    ),
    (
        "pi_n",
        dict(m=0, n=1, au_val="1", bu_val="33", u_val="16"),
        "<",
        "3/2",
        "1/3",
        "\\left(\\dfrac{3}{2}\\right)^{3} - \\left(\\pi^\\dfrac{1}{3}\\right)^{3} = "
        "\\int_0^1 \\frac{\\left(1 - x^{2}\\right) \\left(33 x^{2} + 1\\right)}"
        "{16 x^{2} + 16} \\mathrm{d} x > 0",
    ),
    (
        "arctan_q",
        dict(m=1, n=3, au_val="1871", bu_val="-708", u_val="500"),
        "<",
        "6/5",
        "2",
        "\\dfrac{6}{5} - \\arctan2 = \\int_0^1 \\frac{x^{2} \\left(1 - x^{2}\\right)^{3} "
        "\\cdot \\left(1871 - 708 x^{2}\\right)}{2000 x^{2} + 500} \\mathrm{d} x > 0",
    ),
    (
        "arctan_q",
        dict(m=0, n=1, au_val="1", bu_val="36", u_val="390"),
        "<",
        "3/5",
        "2/3",
        "\\dfrac{3}{5} - \\arctan\\dfrac{2}{3} = \\int_0^1 \\frac{\\left(3 - 3 x^{2}\\right) "
        "\\left(36 x^{2} + 1\\right)}{520 x^{2} + 1170} \\mathrm{d} x > 0",
    ),
    (
        "arctan_q",
        dict(m=0, n=1, au_val="23", bu_val="-19", u_val="5820"),
        ">",
        "31/97",
        "1/3",
        "\\arctan\\dfrac{1}{3} - \\dfrac{31}{97} = \\int_0^1 \\frac{\\left(3 - 3 x^{2}\\right) "
        "\\left(23 - 19 x^{2}\\right)}{1940 x^{2} + 17460} \\mathrm{d} x > 0",
    ),
    (
        "arctan_q",
        dict(m=1, n=1, au_val="1", bu_val="-1", u_val="4"),
        "<",
        "4/5",
        "1",
        "\\dfrac{4}{5} - \\arctan1 = \\int_0^1 \\frac{x^{2} \\left(1 - x^{2}\\right)^{2}}"
        "{4 x^{2} + 4} \\mathrm{d} x > 0",
    ),
    (
        "arccot_q",
        dict(m=1, n=3, au_val="1871", bu_val="-708", u_val="500"),
        "<",
        "6/5",
        "1/2",
        "\\dfrac{6}{5} - \\mathrm{arccot}\\dfrac{1}{2} = \\int_0^1 \\frac{x^{2} "
        "\\left(1 - x^{2}\\right)^{3} \\cdot \\left(1871 - 708 x^{2}\\right)}"
        "{2000 x^{2} + 500} \\mathrm{d} x > 0",
    ),
    (
        "arccot_q",
        dict(m=1, n=1, au_val="2", bu_val="1", u_val="80"),
        ">",
        "23/50",
        "2",
        "\\mathrm{arccot}2 - \\dfrac{23}{50} = \\int_0^1 \\frac{x^{2} \\cdot "
        "\\left(1 - x^{2}\\right) \\left(x^{2} + 2\\right)}{20 x^{2} + 80} \\mathrm{d} x > 0",
    ),
    (
        "catalan",
        dict(m=5, n=5, au_val="733466923", bu_val="-640846037", u_val="43978014720"),
        "<",
        "109/119",
        "1",
        "\\dfrac{109}{119} - C = \\int_0^1 \\frac{x^{10} \\left(1 - x^{2}\\right)^{5} "
        "\\cdot \\left(733466923 - 640846037 x^{2}\\right)}{43978014720 x^{2} + 43978014720}"
        "\\ln(1/x) \\mathrm{d} x > 0",
    ),
    (
        "catalan",
        dict(m=0, n=0, au_val="0", bu_val="2", u_val="1"),
        "<",
        "2",
        "2",
        "2 - 2C = \\int_0^1 \\frac{2 x^{2}}{x^{2} + 1}\\ln(1/x) \\mathrm{d} x > 0",
    ),
    (
        "zeta3",
        dict(m=3, n=4, au_val="105135", bu_val="3786824", u_val="11045067"),
        ">",
        "119/99",
        "1",
        "\\zeta(3) - \\dfrac{119}{99} = \\int_0^1 \\frac{x^{7} \\left(1 - x^{2}\\right)^{4} "
        "\\cdot \\left(3786824 x^{2} + 105135\\right)}{11045067 x^{2} + 11045067}"
        "\\ln^2(x) \\mathrm{d} x > 0",
    ),
    (
        "zeta3",
        dict(m=1, n=0, au_val="32", bu_val="64", u_val="3"),
        ">",
        "2",
        "2",
        "2\\zeta(3) - 2 = \\int_0^1 \\frac{x^{3} \\cdot \\left(64 x^{2} + 32\\right)}"
        "{3 x^{2} + 3}\\ln^2(x) \\mathrm{d} x > 0",
    ),
]


@pytest.mark.parametrize("kind,params,comp,bound,power,expected", RENDER_CASES)
def test_render_equation(kind, params, comp, bound, power, expected):
    full = {"cu_val": "0", **params}
    got = quadlog.render_equation(full, kind, F(power), comp, F(bound))
    assert got == expected


# --------------------------------------------------- numeric identity self-check

CONST_VALUE = {
    "pi": lambda coef: mp.pi * coef,
    "pi_n": lambda power: mp.pi ** int(power),
    "catalan": lambda coef: mp.catalan * coef,
    "zeta3": lambda coef: mp.zeta(3) * coef,
    "arctan_q": lambda q: mp.atan(mp.mpf(str(q))),
    "arccot_q": lambda q: mp.atan(mp.mpf(str(1 / q))),
}


@pytest.mark.parametrize("kind,power,comp,bound", [c[:4] for c in SITE_CASES])
def test_identity_holds_numerically(kind, power, comp, bound):
    """Every emitted proof: integral == bound-const and integrand keeps its sign."""
    r = quadlog.prove(kind, F(power), comp, F(bound))["parameters"]
    m, n = int(r["m"]), int(r["n"])
    a, b = F(r["a_val"]), F(r["b_val"])
    q = mp.mpf(str(F(power) if kind == "arctan_q" else (1 / F(power) if kind == "arccot_q" else 1)))
    r_ln = {"pi": 0, "catalan": 1, "zeta3": 2}.get(
        kind, F(power).numerator - 1 if kind == "pi_n" else 0
    )
    odd = (kind == "zeta3") or (kind == "pi_n" and F(power).numerator % 2 == 0)

    def f(x):
        return (
            x ** (2 * m + odd)
            * (1 - x**2) ** n
            * (a + b * x**2)
            * mp.log(1 / x) ** r_ln
            / (1 + q**2 * x**2)
        )

    with mp.workdps(50):
        val = mp.quad(f, [0, 1])

        pw = F(power)
        if kind == "pi_n":
            # the kernel proves pi^p vs bound^q where power = p/q
            target = mp.mpf(str(F(bound) ** pw.denominator)) - mp.pi**pw.numerator
        else:
            const = CONST_VALUE[kind](mp.mpf(str(pw)))
            target = mp.mpf(str(F(bound))) - const
        if comp == ">":
            target = -target
        assert abs(val - target) < mp.mpf("1e-25") * max(1, abs(target))

        xs = [mp.mpf(i) / 200 for i in range(1, 200)]
        assert all(f(x) >= -mp.mpf("1e-30") for x in xs)
