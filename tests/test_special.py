"""beta_family (golden/varpi/gauss) + gamma kernel tests.

Ground truth: live zhuyidao.net responses probed into bench/probes/ and
bench/data/golden.jsonl — parameters, solution strings and rendered equations
are asserted byte-identically.

Two site quirks are pinned down here:

- gauss '<': the correct solve runs only over m<=4; when it exhausts and the
  claim is numerically true, the site solves a fixed m=5 system whose rows are
  the two basis-moment vectors (the transpose of the correct system). The
  emitted "proof" is mathematically false — asserted explicitly below — but
  reproducing it is required for wire parity.
- '>' (varpi/gauss): on exhausted search the site falls back to
  ∫ t·x^k(1-x)·sqrt(1-x^4)/π dx + b (k=0 gauss, k=1 varpi), emitted with
  a_val=0, au/u = t, bu=0, cu=1. That identity IS exact.
"""

from fractions import Fraction

import pytest
import sympy as sp
from mpmath import mp

from attention_calculator import engine, solve
from attention_calculator.integrand import lhs_mpf, reconstruct
from attention_calculator.kernels import beta_family, gamma
from attention_calculator.moment import combine
from attention_calculator.render import coerce_params

F = Fraction
x = sp.symbols("x")


# ------------------------------------------------------------------- moments


def test_lemniscate_seed_moments():
    """J_k = ∫x^k/√(1-x⁴): J0=ϖ/2, J1=π/4, J2=G⁻¹/2, J3=1/2."""
    p, c, q, d = beta_family.lemniscate_table(1)
    assert beta_family.j_moment(0, p, c, q, d) == {"varpi": F(1, 2)}
    assert beta_family.j_moment(1, p, c, q, d) == {"pi": F(1, 4)}
    assert beta_family.j_moment(2, p, c, q, d) == {"gauss_inv": F(1, 2)}
    assert beta_family.j_moment(3, p, c, q, d) == {"1": F(1, 2)}


def test_lemniscate_recurrences():
    """J_{k+4}/J_k = (k+1)/(k+3) for each of the four interleaved sequences."""
    p, c, q, d = beta_family.lemniscate_table(8)
    for j in range(7):
        # p_{j+1}/p_j = (4j+1)/(4j+3), c: (2j+1)/(2j+2),
        # q: (4j+3)/(4j+5), d: (2j+2)/(2j+3)
        assert p[j + 1] == p[j] * F(4 * j + 1, 4 * j + 3)
        assert c[j + 1] == c[j] * F(2 * j + 1, 2 * j + 2)
        assert q[j + 1] == q[j] * F(4 * j + 3, 4 * j + 5)
        assert d[j + 1] == d[j] * F(2 * j + 2, 2 * j + 3)
    assert q[0] == 1 and d[0] == 1 and p[1] == F(1, 3) and c[1] == F(1, 2)


def test_golden_moment_closed_form():
    """S_0 = ∫√(x+4) = 2(5√5-8)/3 -> {"1": -26/3, "phi": 20/3}."""
    assert beta_family.golden_moment(0) == {"1": F(-26, 3), "phi": F(20, 3)}


def test_solved_identity_is_exact():
    """A normal-path proof satisfies a·B0 + b·B1 == target symbol-wise."""
    res = beta_family.prove("gauss", F(1), "<", F(9, 10))
    m = res["parameters"]["m"]
    basis = [beta_family.lemniscate_basis("gauss", "<", m, i) for i in (0, 1)]
    lhs = combine([F(res["parameters"]["a_val"]), F(res["parameters"]["b_val"])], basis)
    assert lhs == {"gauss_inv": F(9, 10), "1": F(-1)}


# ----------------------------------------------------------- site parity: params

# (kind, power, comp, bound, expected subset of parameters, solution string)
SITE_PARAMS = [
    (
        "golden",
        "1",
        ">",
        "0",
        dict(m=0, n=0, a_val="89/160", b_val="13/32", au_val="89", bu_val="65", u_val="160"),
        "a = 89/160, b = 13/32",
    ),
    (
        "golden",
        "1",
        "<",
        "2",
        dict(m=0, n=0, a_val="11/160", b_val="7/32", au_val="11", bu_val="35", u_val="160"),
        "a = 11/160, b = 7/32",
    ),
    (
        "golden",
        "1",
        ">",
        "8/5",
        dict(m=0, n=1, a_val="21/800", b_val="-21/800", u_val="800"),
        "a = 21/800, b = -21/800",
    ),
    (
        "golden",
        "1",
        "<",
        "1619/1000",
        dict(m=1, n=2, a_val="801/128000", b_val="-231/128000", u_val="128000"),
        "a = 801/128000, b = -231/128000",
    ),
    ("varpi", "1", "<", "3", dict(m=0, n=1, a_val="6", b_val="0", u_val="1"), "a = 6, b = 0"),
    (
        "varpi",
        "1",
        "<",
        "14/5",
        dict(m=0, n=1, a_val="0", b_val="42/5", u_val="5"),
        "a = 0, b = 42/5",
    ),
    (
        "varpi",
        "1",
        "<",
        "31/10",
        dict(m=0, n=1, a_val="9", b_val="-21/5", u_val="5"),
        "a = 9, b = -21/5",
    ),
    (
        "varpi",
        "1",
        "<",
        "11/4",
        dict(m=1, n=1, a_val="33/20", b_val="33/4", au_val="33", bu_val="165", u_val="20"),
        "a = 33/20, b = 33/4",
    ),
    (
        "varpi",
        "1",
        "<",
        "53/20",
        dict(
            m=10,
            n=1,
            a_val="3113915036463/1960525168640",
            b_val="36939482878011/1960525168640",
            u_val="1960525168640",
        ),
        "a = 3113915036463/1960525168640, b = 36939482878011/1960525168640",
    ),
    ("varpi", "1", ">", "2", dict(m=0, n=0, a_val="4", b_val="0", u_val="1"), "a = 4, b = 0"),
    (
        "varpi",
        "1",
        ">",
        "5/2",
        dict(m=2, n=0, a_val="53/7", b_val="26/7", u_val="7"),
        "a = 53/7, b = 26/7",
    ),
    # '>' sqrt(1-x^4)-numerator fallback: a_val=0, au/u = t, cu_val=1
    (
        "varpi",
        "1",
        ">",
        "0",
        dict(m=0, n=0, a_val="0", b_val="1", au_val="0", bu_val="0", cu_val="1", u_val="1"),
        "a = 0, b = 0",
    ),
    (
        "varpi",
        "1",
        ">",
        "1/2",
        dict(m=0, n=0, a_val="0", b_val="11/16", au_val="5", bu_val="0", cu_val="1", u_val="2"),
        "a = 0, b = 0",
    ),
    ("gauss", "1", "<", "9/10", dict(m=0, n=1, a_val="0", b_val="3", u_val="1"), "a = 0, b = 3"),
    ("gauss", "1", "<", "1", dict(m=0, n=1, a_val="2", b_val="0", u_val="1"), "a = 2, b = 0"),
    (
        "gauss",
        "1",
        "<",
        "7/8",
        dict(m=1, n=1, a_val="0", b_val="15/4", u_val="4"),
        "a = 0, b = 15/4",
    ),
    (
        "gauss",
        "1",
        "<",
        "171/200",
        dict(m=4, n=1, a_val="7011/3520", b_val="5157/1600", u_val="17600"),
        "a = 7011/3520, b = 5157/1600",
    ),
    # '<' transposed m=5 fallback records (mathematically false, site-emitted)
    (
        "gauss",
        "1",
        "<",
        "21/25",
        dict(m=5, n=1, a_val="765102/209", b_val="2495493/800", u_val="167200"),
        "a = 765102/209, b = 2495493/800",
    ),
    (
        "gauss",
        "1",
        "<",
        "17/20",
        dict(m=5, n=1, a_val="769080/209", b_val="8027019/2560", u_val="535040"),
        "a = 769080/209, b = 8027019/2560",
    ),
    (
        "gauss",
        "1",
        "<",
        "167/200",
        dict(m=5, n=1, a_val="763113/209", b_val="79648569/25600", u_val="5350400"),
        "a = 763113/209, b = 79648569/25600",
    ),
    (
        "gauss",
        "1",
        ">",
        "4/5",
        dict(m=2, n=0, a_val="6/5", b_val="44/5", u_val="5"),
        "a = 6/5, b = 44/5",
    ),
    ("gauss", "1", ">", "1/2", dict(m=0, n=0, a_val="2", b_val="0", u_val="1"), "a = 2, b = 0"),
    ("gauss", "1", ">", "3/8", dict(m=0, n=0, a_val="3", b_val="-3", u_val="1"), "a = 3, b = -3"),
    ("gauss", "1", ">", "3/4", dict(m=0, n=0, a_val="0", b_val="6", u_val="1"), "a = 0, b = 6"),
    # '>' sqrt fallback records
    (
        "gauss",
        "1",
        ">",
        "0",
        dict(m=0, n=0, a_val="0", b_val="3/8", au_val="3", bu_val="0", cu_val="1", u_val="1"),
        "a = 0, b = 0",
    ),
    (
        "gauss",
        "1",
        ">",
        "1/4",
        dict(m=0, n=0, a_val="0", b_val="1/8", au_val="3", bu_val="0", cu_val="1", u_val="1"),
        "a = 0, b = 0",
    ),
    (
        "gauss",
        "2",
        ">",
        "0",
        dict(m=0, n=0, a_val="0", b_val="3/4", au_val="6", bu_val="0", cu_val="1", u_val="1"),
        "a = 0, b = 0",
    ),
    # gamma composite kernel; u_val=0 records print a degenerate solution line
    (
        "gamma",
        "1",
        ">",
        "0",
        dict(m=0, n=0, a_val="1/2", b_val="0", c_val="1", cu_val="0", u_val="0"),
        "a = 0, b = 0",
    ),
    (
        "gamma",
        "1",
        "<",
        "1",
        dict(m=0, n=0, a_val="1/4", b_val="0", c_val="2", cu_val="0", u_val="0"),
        "a = 0, b = 0",
    ),
    (
        "gamma",
        "1",
        ">",
        "1/2",
        dict(m=0, n=0, a_val="0", b_val="0", c_val="1", cu_val="0", u_val="0"),
        "a = 0, b = 0",
    ),
    (
        "gamma",
        "1",
        "<",
        "3/5",
        dict(
            m=2,
            n=4,
            a_val="937/56",
            b_val="4615/42",
            c_val="0",
            cu_val="5",
            au_val="2811",
            bu_val="18460",
            u_val="168",
        ),
        "a = 937/56, b = 4615/42",
    ),
    (
        "gamma",
        "1",
        ">",
        "57/100",
        dict(
            m=4,
            n=2,
            a_val="67",
            b_val="-304/5",
            c_val="0",
            cu_val="4",
            au_val="335",
            bu_val="-304",
            u_val="5",
        ),
        "a = 67, b = -304/5",
    ),
    (
        "gamma",
        "1",
        "<",
        "29/50",
        dict(
            m=7,
            n=8,
            a_val="135299452736339/425425",
            b_val="-122801520422253/425425",
            c_val="0",
            cu_val="16",
            u_val="425425",
        ),
        "a = 135299452736339/425425, b = -122801520422253/425425",
    ),
]


@pytest.mark.parametrize("kind,power,comp,bound,expected,solution", SITE_PARAMS)
def test_prove_site_parity(kind, power, comp, bound, expected, solution):
    res = solve.prove(kind, power, comp, bound)
    p = res["parameters"]
    for k_, want in expected.items():
        assert str(p[k_]) == str(want), f"{k_}: {p[k_]} != {want}"
    assert res["solution"] == solution


# ------------------------------------------------------ render byte-parity

RENDER_CASES = [
    (
        "golden",
        "1",
        ">",
        "0",
        r"\phi - 0 = \int_0^1 \frac{\sqrt{x + 4} \cdot \left(65 x + 89\right)}"
        r"{160} \mathrm{d} x > 0",
    ),
    (
        "golden",
        "1",
        "<",
        "2",
        r"2 - \phi = \int_0^1 \frac{\sqrt{x + 4} \cdot \left(35 x + 11\right)}"
        r"{160} \mathrm{d} x > 0",
    ),
    (
        "golden",
        "1",
        ">",
        "3/2",
        r"\phi - \dfrac{3}{2} = \int_0^1 \frac{\left(7 - 5 x\right) \sqrt{x + 4}}"
        r"{80} \mathrm{d} x > 0",
    ),
    (
        "golden",
        "1",
        "<",
        "17/10",
        r"\dfrac{17}{10} - \phi = \int_0^1 \frac{\left(1 - x\right) \sqrt{x + 4}"
        r" \cdot \left(14 x + 11\right)}{200} \mathrm{d} x > 0",
    ),
    (
        "varpi",
        "1",
        ">",
        "2",
        r"1 - 2\varpi^{-1} = \int_0^1  4 x\dfrac{(1-x)}{\pi\sqrt{1-x^4}}"
        r" \mathrm{d} x > 0",
    ),
    (
        "varpi",
        "1",
        ">",
        "5/2",
        r"1 - \dfrac{5}{2}\varpi^{-1} = \int_0^1  \frac{x^{9} \cdot \left(26 x^{4}"
        r" + 53\right)}{7}\dfrac{(1-x)}{\pi\sqrt{1-x^4}} \mathrm{d} x > 0",
    ),
    (
        "varpi",
        "1",
        "<",
        "3",
        r"3 - \varpi = \int_0^1 6 x^{3}\dfrac{(1-x)}{\sqrt{1-x^4}}"
        r" \mathrm{d} x > 0",
    ),
    (
        "varpi",
        "1",
        "<",
        "14/5",
        r"\dfrac{14}{5} - \varpi = \int_0^1 \frac{42 x^{7}}{5}\dfrac{(1-x)}"
        r"{\sqrt{1-x^4}} \mathrm{d} x > 0",
    ),
    # '>' sqrt fallback renders: integrand + additive constant
    (
        "varpi",
        "1",
        ">",
        "0",
        r"1 - 0\varpi^{-1} = \int_0^1 0\dfrac{\sqrt{1-x^4}}{\pi} \mathrm{d} x"
        r"+1 > 0",
    ),
    (
        "varpi",
        "1",
        ">",
        "1/2",
        r"1 - \dfrac{1}{2}\varpi^{-1} = \int_0^1 \frac{5 x \left(1 - x\right)}{2}"
        r"\dfrac{\sqrt{1-x^4}}{\pi} \mathrm{d} x+\dfrac{11}{16} > 0",
    ),
    (
        "gauss",
        "1",
        ">",
        "0",
        r"G-0 = \int_0^1  \left(3 - 3 x\right)\dfrac{\sqrt{1-x^4}}{\pi}"
        r" \mathrm{d} x+\dfrac{3}{8} > 0",
    ),
    (
        "gauss",
        "1",
        ">",
        "1/4",
        r"G-\dfrac{1}{4} = \int_0^1  \left(3 - 3 x\right)\dfrac{\sqrt{1-x^4}}{\pi}"
        r" \mathrm{d} x+\dfrac{1}{8} > 0",
    ),
    (
        "gauss",
        "1",
        ">",
        "1/2",
        r"G-\dfrac{1}{2} = \int_0^1  2\dfrac{(1-x)}{\pi\sqrt{1-x^4}}"
        r" \mathrm{d} x > 0",
    ),
    (
        "gauss",
        "1",
        ">",
        "3/8",
        r"G-\dfrac{3}{8} = \int_0^1  3 - 3 x^{4}\dfrac{(1-x)}{\pi\sqrt{1-x^4}}"
        r" \mathrm{d} x > 0",
    ),
    (
        "gauss",
        "1",
        ">",
        "4/5",
        r"G-\dfrac{4}{5} = \int_0^1  \frac{x^{8} \cdot \left(44 x^{4} + 6\right)}"
        r"{5}\dfrac{(1-x)}{\pi\sqrt{1-x^4}} \mathrm{d} x > 0",
    ),
    (
        "gauss",
        "1",
        "<",
        "9/10",
        r"\dfrac{9}{10}G^{-1}-1 = \int_0^1 3 x^{6}\dfrac{(1-x)}{\sqrt{1-x^4}}"
        r" \mathrm{d} x > 0",
    ),
    (
        "gauss",
        "1",
        "<",
        "1",
        r"G^{-1}-1 = \int_0^1 2 x^{2}\dfrac{(1-x)}{\sqrt{1-x^4}}"
        r" \mathrm{d} x > 0",
    ),
    (
        "gauss",
        "1",
        "<",
        "21/25",
        r"\dfrac{21}{25}G^{-1}-1 = \int_0^1 \frac{x^{22} \cdot \left(521558037"
        r" x^{4} + 612081600\right)}{167200}\dfrac{(1-x)}{\sqrt{1-x^4}}"
        r" \mathrm{d} x > 0",
    ),
    (
        "gamma",
        "1",
        ">",
        "0",
        r"\gamma - 0 = \int_0^1 \left[\left(\dfrac{1}{1-x}+\dfrac{1}{\ln(x)}"
        r"-\dfrac{1}{2}\right)+\dfrac{1}{2}\right] \mathrm{d} x > 0",
    ),
    (
        "gamma",
        "1",
        "<",
        "1",
        r"1 - \gamma =\int_0^1 \left[\left(\dfrac{2-x}{2}-\dfrac{1}{1-x}"
        r"-\dfrac{1}{\ln(x)}\right)+\dfrac{1}{4}\right] \mathrm{d} x > 0",
    ),
    (
        "gamma",
        "1",
        ">",
        "1/2",
        r"\gamma - \dfrac{1}{2} = \int_0^1 \left(\dfrac{1}{1-x}+\dfrac{1}{\ln(x)}"
        r"-\dfrac{1}{2}\right) \mathrm{d} x > 0",
    ),
    (
        "gamma",
        "1",
        "<",
        "3/5",
        r"\dfrac{3}{5} - \gamma =\int_0^1 \left[x^{5}\left(\dfrac{2-x}{2}"
        r"-\dfrac{1}{1-x}-\dfrac{1}{\ln(x)}\right)+\frac{x^{2} \left(1 - x\right)"
        r"^{4} \cdot \left(18460 x + 2811\right)}{168 \left(5 x + 1\right)^{4}}"
        r"\right] \mathrm{d} x > 0",
    ),
    (
        "gamma",
        "1",
        ">",
        "57/100",
        r"\gamma - \dfrac{57}{100} = \int_0^1 \left[x^{4}\left(\dfrac{1}{1-x}"
        r"+\dfrac{1}{\ln(x)}-\dfrac{1}{2}\right)+\frac{x^{4} \left(1 - x\right)"
        r"^{2} \cdot \left(335 - 304 x\right)}{5 \left(4 x + 1\right)^{4}}"
        r"\right] \mathrm{d} x > 0",
    ),
]


@pytest.mark.parametrize("kind,power,comp,bound,equation", RENDER_CASES)
def test_render_equation_bytes(kind, power, comp, bound, equation):
    res = solve.prove(kind, power, comp, bound)
    typed = coerce_params(
        {
            k_: str(res["parameters"][k_])
            for k_ in ("m", "n", "a_val", "b_val", "c_val", "u_val", "au_val", "bu_val", "cu_val")
        }
    )
    module = beta_family if kind != "gamma" else gamma
    assert module.render_equation(typed, kind, power, comp, bound) == equation


# --------------------------------------------------------------- error paths


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("gauss", "1", "<", "0"),  # 0 < G, claim false
        ("gauss", "1", "<", "4/5"),  # 4/5 < G
        ("gauss", "1", "<", "83/100"),
        ("gauss", "2", "<", "1"),  # 2G > 1, claim false (power is a multiple)
        ("gauss", "1", ">", "1"),  # G < 1
        ("gauss", "1", ">", "9/10"),
        ("gauss", "1", ">", "21/25"),
        ("varpi", "1", "<", "2"),  # 2 < ϖ
        ("varpi", "1", "<", "5/2"),
        ("varpi", "1", ">", "3"),  # 3 > ϖ
        ("golden", "1", "<", "3/2"),  # 3/2 < φ
    ],
)
def test_wrong_direction(kind, power, comp, bound):
    with pytest.raises(engine.WrongDirection):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("gauss", "1", ">", "5/6"),  # true but unproven (sqrt b<0)
        ("gauss", "1", ">", "83/100"),
        ("gauss", "1", ">", "417/500"),
        ("varpi", "1", "<", "263/100"),  # true, needs m > 10
        ("varpi", "1", "<", "21/8"),  # transposed fallback has a<0 -> never fires
        ("varpi", "1", ">", "13/5"),  # sqrt b = -5/8 < 0 -> unproven
        ("varpi", "1", ">", "1311/500"),
    ],
)
def test_no_solution(kind, power, comp, bound):
    with pytest.raises(engine.NoSolution):
        solve.prove(kind, power, comp, bound)


# --------------------------------------------- the m=5 transposed site bug


def test_gauss_lt_transposed_fallback_is_site_bug():
    """gauss '<' m=5 records satisfy basis_i·x = T_i, not a·B0+b·B1 == T.

    The site solves the transpose of the correct system, so the emitted
    equation's identity is false; we reproduce it for wire parity. This test
    pins both facts: our params match the site, and the claimed identity does
    not hold symbolically.
    """
    res = beta_family.prove("gauss", F(1), "<", F(21, 25))
    p = res["parameters"]
    a, b = F(p["a_val"]), F(p["b_val"])
    basis = [beta_family.lemniscate_basis("gauss", "<", p["m"], i) for i in (0, 1)]
    target = {"gauss_inv": F(21, 25), "1": F(-1)}
    # transposed: basis_i · (a,b) == target_i in [gauss_inv, "1"] order
    for bvec, tval in zip(basis, (target["gauss_inv"], target["1"]), strict=True):
        assert a * bvec.get("gauss_inv", F(0)) + b * bvec.get("1", F(0)) == tval
    # while the true coefficient matching does NOT hold
    assert combine([a, b], basis) != target


def test_varpi_lt_transposed_never_fires():
    """varpi '<' transposed solve at the fallback m always gives a < 0, so the
    site never emits it: true-but-hard bounds 404 with 未找到 instead."""
    basis = [beta_family.lemniscate_basis("varpi", "<", 11, i) for i in (0, 1)]
    rows = [[v.get("varpi", F(0)), v.get("1", F(0))] for v in basis]
    coeffs = engine.gauss_solve(rows, [F(-1), F(53, 20)])  # r = 53/20
    assert coeffs[0] < 0


# ------------------------------------------- sqrt fallback is a real identity


@pytest.mark.parametrize(
    "kind,comp,bound,poly",
    [
        ("gauss", ">", "0", 3 * (1 - x) * sp.sqrt(1 - x**4) / sp.pi),
        ("gauss", ">", "1/4", 3 * (1 - x) * sp.sqrt(1 - x**4) / sp.pi),
        ("varpi", ">", "1/2", sp.Rational(5, 2) * x * (1 - x) * sp.sqrt(1 - x**4) / sp.pi),
    ],
)
def test_sqrt_fallback_identity(kind, comp, bound, poly):
    """The '>' fallback record's integrand + b_val really integrates to the
    LHS: ∫ poly·sqrt(1-x^4)/π + b = power·C - bound (exact, mpmath check)."""
    res = solve.prove(kind, "1", comp, bound)
    f, lo, hi = reconstruct(kind, comp, F(1), res["parameters"])
    assert sp.simplify(f - (poly + F(res["parameters"]["b_val"]))) == 0
    mp.dps = 35
    val = mp.quad(sp.lambdify(x, f, modules="mpmath"), [mp.mpf(str(lo)), mp.mpf(str(hi))])
    lhs = lhs_mpf(kind, comp, F(1), F(bound))
    assert abs(val - lhs) < mp.mpf("1e-25")
