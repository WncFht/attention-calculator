"""log kernel family tests (ln_q, ln_q_square, artanh_q, arcoth_q).

Ground truth: live zhuyidao.net responses probed during development
(parameters and rendered equations, byte-identical) plus mpmath quadrature
for the moment closed forms. The denominator power s = max(m, n, 1) and the
(m, n) parity traversal order are both site-verified (docs/log-notes.md).
"""

from fractions import Fraction as F

import pytest
from mpmath import mp

from attention_calculator import engine, solve
from attention_calculator.kernels import log_family

mp.dps = 50


def tomf(f):
    """Fraction -> mpf (mpmath won't take Fraction directly)."""
    return mp.mpf(f.numerator) / f.denominator


# ---------------------------------------------------------------- moment closed forms


@pytest.mark.parametrize(
    "c,s,k,expected",
    [
        (F(1), 1, 0, {"ln": F(1)}),  # ∫1/(1+x) = ln2
        (F(1), 1, 1, {"ln": F(-1), "1": F(1)}),  # ∫x/(1+x) = 1-ln2
        (F(1), 2, 0, {"1": F(1, 2)}),  # ∫1/(1+x)^2 = 1/2
        (F(1), 2, 1, {"ln": F(1), "1": F(-1, 2)}),  # ∫x/(1+x)^2 = ln2-1/2
        (F(2), 1, 0, {"ln": F(1, 2)}),  # ∫1/(1+2x) = (ln3)/2
        (F(1, 2), 2, 1, {"ln": F(4), "1": F(-4, 3)}),  # ∫x/(1+x/2)^2
    ],
)
def test_iota_hand(c, s, k, expected):
    assert log_family.iota(c, s, k) == expected


@pytest.mark.parametrize(
    "c,s,k,expected",
    [
        (F(1), 1, 0, {"ln2": F(1, 2)}),  # ln^2 2 / 2
        (F(1), 1, 1, {"ln2": F(-1, 2), "ln": F(2), "1": F(-1)}),  # ∫x ln(1+x)/(1+x)
        (F(1), 2, 0, {"ln": F(-1, 2), "1": F(1, 2)}),  # ∫ln(1+x)/(1+x)^2
        (F(2), 1, 0, {"ln2": F(1, 4)}),  # ln^2 3 / 4
    ],
)
def test_kappa_hand(c, s, k, expected):
    assert log_family.kappa(c, s, k) == expected


@pytest.mark.parametrize(
    "m,n,j,q,square",
    [
        (0, 0, 0, F(2), False),
        (2, 3, 1, F(3), False),
        (1, 2, 0, F(3, 2), False),
        (3, 0, 1, F(5), False),
        (0, 0, 2, F(2), True),
        (1, 3, 1, F(4), True),
        (2, 0, 0, F(3, 2), True),
    ],
)
def test_basis_moment_numeric(m, n, j, q, square):
    """Closed-form basis moments vs numeric quadrature at s = max(m,n,1)."""
    c = q - 1
    s = max(m, n, 1)
    mom = log_family.basis_moment(c, s, m, n, j, square)
    cq = tomf(c)
    qq = tomf(q)
    val = (
        tomf(mom.get("ln2", F(0))) * mp.log(qq) ** 2
        + tomf(mom.get("ln", F(0))) * mp.log(qq)
        + tomf(mom.get("1", F(0)))
    )

    def f(x):
        base = x ** (m + j) * (1 - x) ** n / (1 + cq * x) ** s
        return base * mp.log(1 + cq * x) if square else base

    want = mp.quad(f, [0, 1])
    assert mp.almosteq(val, want, rel_eps=mp.mpf("1e-30"))


def test_mn_order_log_parity():
    """Site traversal: balanced pairs first; inside a mirror pair the smaller
    m wins on odd m+n, the larger m on even m+n."""
    from attention_calculator import engine

    got = list(engine.mn_order(3))
    assert got[:12] == [
        (0, 0),
        (0, 1),
        (1, 0),
        (1, 1),
        (2, 0),
        (0, 2),
        (1, 2),
        (2, 1),
        (0, 3),
        (3, 0),
        (2, 2),
        (3, 1),
    ]


# ---------------------------------------------------------------- site parity (probed)


@pytest.mark.parametrize(
    "kind,power,comp,bound,expected,solution",
    [
        (
            "ln_q",
            "3",
            "<",
            "10/9",
            dict(
                m=1,
                n=2,
                a_val="0",
                b_val="4/3",
                c_val="0",
                au_val="0",
                bu_val="4",
                cu_val="0",
                u_val="3",
            ),
            "a = 0, b = 4/3",
        ),
        (
            "ln_q",
            "2",
            ">",
            "11/16",
            dict(
                m=0,
                n=2,
                a_val="0",
                b_val="1/8",
                c_val="0",
                au_val="0",
                bu_val="1",
                cu_val="0",
                u_val="8",
            ),
            "a = 0, b = 1/8",
        ),
        (
            "ln_q",
            "2",
            "<",
            "347/500",
            dict(
                m=2,
                n=2,
                a_val="31/500",
                b_val="-2/125",
                c_val="0",
                au_val="31",
                bu_val="-8",
                cu_val="0",
                u_val="500",
            ),
            "a = 31/500, b = -2/125",
        ),
        # article-famous ln5 bound; s = max(6,7) = 7
        (
            "ln_q",
            "5",
            "<",
            "375/233",
            dict(
                m=6,
                n=7,
                a_val="3144/1165",
                b_val="1491584/8155",
                c_val="0",
                au_val="22008",
                bu_val="1491584",
                cu_val="0",
                u_val="8155",
            ),
            "a = 3144/1165, b = 1491584/8155",
        ),
        (
            "ln_q",
            "101/100",
            "<",
            "1/50",
            dict(
                m=0,
                n=0,
                a_val="1/100",
                b_val="1/5000",
                c_val="0",
                au_val="50",
                bu_val="1",
                cu_val="0",
                u_val="5000",
            ),
            "a = 1/100, b = 1/5000",
        ),
        (
            "ln_q",
            "355/113",
            "<",
            "23/20",
            dict(
                m=3,
                n=1,
                a_val="1070835569/489142083",
                b_val="-272755780/163047361",
                c_val="0",
                au_val="1070835569",
                bu_val="-818267340",
                cu_val="0",
                u_val="489142083",
            ),
            "a = 1070835569/489142083, b = -272755780/163047361",
        ),
        # ln^2: const-P case and the parity-order witness (site picks (2,0)
        # over the also-feasible (0,2) — engine.mn_order would get this wrong)
        (
            "ln_q_square",
            "2",
            "<",
            "1/2",
            dict(
                m=0,
                n=2,
                a_val="1/2",
                b_val="0",
                c_val="0",
                au_val="1",
                bu_val="0",
                cu_val="0",
                u_val="2",
            ),
            "a = 1/2, b = 0, c= 0",
        ),
        (
            "ln_q_square",
            "3/2",
            "<",
            "17/100",
            dict(
                m=2,
                n=0,
                a_val="2357/16600",
                b_val="-313/8300",
                c_val="-27/1660",
                au_val="2357",
                bu_val="-626",
                cu_val="-270",
                u_val="16600",
            ),
            "a = 2357/16600, b = -313/8300, c= -27/1660",
        ),
        (
            "ln_q_square",
            "3",
            "<",
            "121/100",
            dict(
                m=1,
                n=4,
                a_val="15343/29700",
                b_val="5581/1980",
                c_val="8168/7425",
                au_val="15343",
                bu_val="83715",
                cu_val="32672",
                u_val="29700",
            ),
            "a = 15343/29700, b = 5581/1980, c= 8168/7425",
        ),
        # artanh solves the halved ln system at q~ = (1+q)/(1-q), reported in c_val
        (
            "artanh_q",
            "1/2",
            "<",
            "11/20",
            dict(
                m=2,
                n=3,
                a_val="1/15",
                b_val="4/9",
                c_val="3",
                au_val="3",
                bu_val="20",
                cu_val="0",
                u_val="45",
            ),
            "a = 1/15, b = 4/9",
        ),
        (
            "artanh_q",
            "1/10",
            "<",
            "11/100",
            dict(
                m=0,
                n=1,
                a_val="89/4950",
                b_val="7/825",
                c_val="11/9",
                au_val="89",
                bu_val="42",
                cu_val="0",
                u_val="4950",
            ),
            "a = 89/4950, b = 7/825",
        ),
        (
            "artanh_q",
            "9/10",
            "<",
            "3/2",
            dict(
                m=4,
                n=5,
                a_val="223830/19",
                b_val="161409348/95",
                c_val="19",
                au_val="1119150",
                bu_val="161409348",
                cu_val="0",
                u_val="95",
            ),
            "a = 223830/19, b = 161409348/95",
        ),
        # arcoth 2 and artanh 1/2 share q~ = 3 -> identical parameters
        (
            "arcoth_q",
            "2",
            "<",
            "11/20",
            dict(
                m=2,
                n=3,
                a_val="1/15",
                b_val="4/9",
                c_val="3",
                au_val="3",
                bu_val="20",
                cu_val="0",
                u_val="45",
            ),
            "a = 1/15, b = 4/9",
        ),
        (
            "arcoth_q",
            "11/10",
            "<",
            "8/5",
            dict(
                m=4,
                n=4,
                a_val="5083468/21",
                b_val="-679520/7",
                c_val="21",
                au_val="5083468",
                bu_val="-2038560",
                cu_val="0",
                u_val="21",
            ),
            "a = 5083468/21, b = -679520/7",
        ),
        (
            "arcoth_q",
            "3/2",
            "<",
            "9/10",
            dict(
                m=1,
                n=2,
                a_val="2",
                b_val="66/5",
                c_val="5",
                au_val="10",
                bu_val="66",
                cu_val="0",
                u_val="5",
            ),
            "a = 2, b = 66/5",
        ),
    ],
)
def test_prove_site_parity(kind, power, comp, bound, expected, solution):
    res = log_family.prove(kind, F(power), comp, F(bound))
    p = res["parameters"]
    for k, want in expected.items():
        assert str(p[k]) == str(want), f"{k}: {p[k]} != {want}"
    assert res["solution"] == solution


def test_prove_via_dispatcher():
    r = solve.prove("ln_q", "3", "<", "10/9")
    assert r["parameters"]["m"] == 1 and r["parameters"]["b_val"] == "4/3"
    assert r["solution"] == "a = 0, b = 4/3"


def test_ln_bound_proof_public():
    """gamma-family reuse entry point."""
    res = log_family.ln_bound_proof(F(5), "<", F(375, 233))
    assert res["parameters"]["m"] == 6 and res["parameters"]["n"] == 7


# ---------------------------------------------------------------- render byte-parity


@pytest.mark.parametrize(
    "kind,power,comp,bound,equation",
    [
        (
            "ln_q",
            "3",
            "<",
            "10/9",
            r"\dfrac{10}{9} - \ln3 = \int_0^1 \frac{4 x^{2} \left(1 - x\right)^{2}}"
            r"{3 \left(2 x + 1\right)^{2}} \mathrm{d} x > 0",
        ),
        (
            "ln_q",
            "2",
            ">",
            "11/16",
            r"\ln2 - \dfrac{11}{16} = \int_0^1 \frac{x \left(1 - x\right)^{2}}"
            r"{8 \left(x + 1\right)^{2}} \mathrm{d} x > 0",
        ),
        (
            "ln_q",
            "2",
            "<",
            "347/500",
            r"\dfrac{347}{500} - \ln2 = \int_0^1 \frac{x^{2} \left(1 - x\right)^{2} "
            r"\cdot \left(31 - 8 x\right)}{500 \left(x + 1\right)^{2}} \mathrm{d} x > 0",
        ),
        (
            "ln_q",
            "3/2",
            ">",
            "2/5",
            r"\ln\dfrac{3}{2} - \dfrac{2}{5} = \int_0^1 \frac{x \left(1 - x\right) "
            r"\left(3 x + 1\right)}{30 x + 60} \mathrm{d} x > 0",
        ),
        (
            "ln_q",
            "101/100",
            "<",
            "1/50",
            r"\dfrac{1}{50} - \ln\dfrac{101}{100} = \int_0^1 \frac{x + 50}"
            r"{50 x + 5000} \mathrm{d} x > 0",
        ),
        (
            "ln_q",
            "3",
            ">",
            "1",
            r"\ln3 - 1 = \int_0^1 \frac{\left(1 - x\right)^{2} \cdot \left(4 x + 1\right)}"
            r"{3 \left(2 x + 1\right)^{2}} \mathrm{d} x > 0",
        ),
        (
            "ln_q",
            "2",
            "<",
            "7/10",
            r"\dfrac{7}{10} - \ln2 = \int_0^1 \frac{x^{2} \cdot \left(1 - x\right)}"
            r"{5 \left(x + 1\right)^{2}} \mathrm{d} x > 0",
        ),
        (
            "ln_q",
            "355/113",
            "<",
            "23/20",
            r"\dfrac{23}{20} - \ln\dfrac{355}{113} = \int_0^1 \frac{x^{3} \cdot "
            r"\left(1 - x\right) \left(1070835569 - 818267340 x\right)}"
            r"{339 \left(242 x + 113\right)^{3}} \mathrm{d} x > 0",
        ),
        (
            "ln_q_square",
            "2",
            "<",
            "1/2",
            r"\dfrac{1}{2} - \ln^22 = \int_0^1 \frac{\left(1 - x\right)^{2} "
            r"\log{\left(x + 1 \right)}}{2 \left(x + 1\right)^{2}} \mathrm{d} x > 0",
        ),
        (
            "ln_q_square",
            "3/2",
            "<",
            "17/100",
            r"\dfrac{17}{100} - \ln^2\dfrac{3}{2} = \int_0^1 \frac{x^{2} \left(- 270 "
            r"x^{2} - 626 x + 2357\right) \log{\left(\frac{x}{2} + 1 \right)}}"
            r"{4150 \left(x + 2\right)^{2}} \mathrm{d} x > 0",
        ),
        (
            "ln_q_square",
            "3/2",
            "<",
            "1",
            r"1 - \ln^2\dfrac{3}{2} = \int_0^1 \frac{\left(4 x^{2} + 10 x + 2\right) "
            r"\log{\left(\frac{x}{2} + 1 \right)}}{x + 2} \mathrm{d} x > 0",
        ),
        (
            "ln_q_square",
            "2",
            ">",
            "12/25",
            r"\ln^22 - \dfrac{12}{25} = \int_0^1 \frac{x \left(1 - x\right)^{3} \cdot "
            r"\left(244 x^{2} + 657 x + 59\right) \log{\left(x + 1 \right)}}"
            r"{3900 \left(x + 1\right)^{3}} \mathrm{d} x > 0",
        ),
        (
            "artanh_q",
            "1/2",
            "<",
            "11/20",
            r"\dfrac{11}{20} - \mathrm{artanh}\dfrac{1}{2} = \int_0^1 \frac{x^{2} "
            r"\left(1 - x\right)^{3} \cdot \left(20 x + 3\right)}{45 \left(2 x + "
            r"1\right)^{3}} \mathrm{d} x > 0",
        ),
        (
            "artanh_q",
            "1/5",
            ">",
            "1/5",
            r"\mathrm{artanh}\dfrac{1}{5} - \dfrac{1}{5} = \int_0^1 \frac{x \left(1 - "
            r"x\right) \left(3 x + 1\right)}{60 x + 120} \mathrm{d} x > 0",
        ),
        (
            "arcoth_q",
            "3/2",
            "<",
            "9/10",
            r"\dfrac{9}{10} - \mathrm{arcoth}\dfrac{3}{2} = \int_0^1 \frac{x \left(1 - "
            r"x\right)^{2} \cdot \left(66 x + 10\right)}{5 \left(4 x + 1\right)^{2}} "
            r"\mathrm{d} x > 0",
        ),
        (
            "arcoth_q",
            "2",
            "<",
            "11/20",
            r"\dfrac{11}{20} - \mathrm{arcoth}2 = \int_0^1 \frac{x^{2} \left(1 - "
            r"x\right)^{3} \cdot \left(20 x + 3\right)}{45 \left(2 x + 1\right)^{3}} "
            r"\mathrm{d} x > 0",
        ),
    ],
)
def test_render_equation_bytes(kind, power, comp, bound, equation):
    res = log_family.prove(kind, F(power), comp, F(bound))
    typed = {
        k: int(res["parameters"][k]) for k in ("m", "n", "au_val", "bu_val", "cu_val", "u_val")
    }
    typed.update({k: F(res["parameters"][k]) for k in ("a_val", "b_val", "c_val")})
    assert log_family.render_equation(typed, kind, F(power), comp, F(bound)) == equation


# ---------------------------------------------------------------- error paths


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("ln_q", "3/2", "<", "0"),  # ln(3/2) > 0
        ("ln_q", "3/2", ">", "1/2"),  # ln(3/2) ≈ 0.405 < 1/2
        ("ln_q", "2", ">", "7/10"),  # ln2 ≈ 0.693 < 7/10
        ("ln_q", "3", ">", "11/10"),  # ln3 ≈ 1.0986 < 1.1
        ("ln_q_square", "2", ">", "1/2"),  # ln^2 2 ≈ 0.480 < 1/2
        ("artanh_q", "1/2", ">", "11/20"),  # artanh(1/2) ≈ 0.549 < 0.55
        ("arcoth_q", "2", ">", "11/20"),  # same q~=3 identity, also false
        # q=5 站端 bug 路径下命题为假（ln²5=2.590290394 > 5709/2204≈2.590290381），
        # 仍先报"反了"而不是 500
        ("ln_q_square", "5", "<", "5709/2204"),
    ],
)
def test_wrong_direction(kind, power, comp, bound):
    with pytest.raises(engine.WrongDirection):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        # 站端 bug：ln_q_square 的 q=5/7 命题为真时求解器必崩 -> 500
        # （5 > 5709/2204 仅真 1.3e-8，亦 500）
        ("ln_q_square", "5", ">", "2"),
        ("ln_q_square", "5", ">", "5709/2204"),
        ("ln_q_square", "7", "<", "4"),
    ],
)
def test_internal_error(kind, power, comp, bound):
    with pytest.raises(engine.InternalError):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("ln_q", "22", "<", "2105/681"),  # true but needs exponents > 10
        ("ln_q", "64", ">", "83/20"),
    ],
)
def test_no_solution(kind, power, comp, bound):
    with pytest.raises(engine.NoSolution):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize(
    "kind,power,bound,msg",
    [
        ("ln_q", "1/2", "0", "请在ln后输入一个大于1的数"),
        ("ln_q", "1", "1/2", "请在ln后输入一个大于1的数"),
        ("ln_q_square", "1/2", "1/2", "请在ln后输入一个大于1的数"),
        ("artanh_q", "1", "1/2", "请在输入一个在(0,1)内的分数，本情况不支持整数"),
        ("artanh_q", "3/2", "1/2", "请在输入一个在(0,1)内的分数，本情况不支持整数"),
        ("arcoth_q", "1", "1/2", "请在输入一个大于1的数"),
        ("arcoth_q", "1/2", "1/2", "请在输入一个大于1的数"),
        ("ln_q", "3", "-1/2", "右侧有理数格式无效"),  # bound checked before domain
        ("ln_q", "1/2", "-7/10", "右侧有理数格式无效"),
        ("ln_q", "-2", "1/2", "左侧系数格式无效"),
    ],
)
def test_domain_and_format_messages(kind, power, bound, msg):
    """Site 404/400 texts, raised as ValueError from the kernel (trig convention)."""
    with pytest.raises(ValueError) as exc:
        solve.prove(kind, power, "<", bound)
    assert str(exc.value) == msg


# --------------------------------------------------- numeric identity self-check

CONST_MPF = {
    "ln_q": mp.log,
    "ln_q_square": lambda q: mp.log(q) ** 2,
    "artanh_q": mp.atanh,
    "arcoth_q": lambda q: mp.atanh(1 / q),
}


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        ("ln_q", "3", "<", "10/9"),
        ("ln_q", "5", "<", "375/233"),
        ("ln_q", "355/113", "<", "23/20"),
        ("ln_q", "2", ">", "11/16"),
        ("ln_q_square", "2", "<", "1/2"),
        ("ln_q_square", "3", "<", "121/100"),
        ("ln_q_square", "3/2", "<", "17/100"),
        ("artanh_q", "1/2", "<", "11/20"),
        ("artanh_q", "9/10", "<", "3/2"),
        ("arcoth_q", "2", "<", "11/20"),
        ("arcoth_q", "3/2", "<", "9/10"),
    ],
)
def test_identity_holds_numerically(kind, power, comp, bound):
    """Every emitted proof: integrand >= 0 and integral == |C - bound|."""
    res = log_family.prove(kind, F(power), comp, F(bound))["parameters"]
    m, n = int(res["m"]), int(res["n"])
    a, b, c2 = F(res["a_val"]), F(res["b_val"]), F(res["c_val"])
    square = kind == "ln_q_square"
    qt = c2 if kind in ("artanh_q", "arcoth_q") else F(power)
    cq = tomf(qt - 1)
    s = max(m, n, 1)
    am, bm = tomf(a), tomf(b)
    cm = tomf(c2) if square else mp.mpf(0)

    def f(x):
        base = x**m * (1 - x) ** n * (am + bm * x + cm * x**2) / (1 + cq * x) ** s
        return base * mp.log(1 + cq * x) if square else base

    val = mp.quad(f, [0, 1])
    const = CONST_MPF[kind](tomf(F(power)))
    target = mp.mpf(str(F(bound))) - const
    if comp == ">":
        target = -target
    assert abs(val - target) < mp.mpf("1e-25") * max(1, abs(target))
    xs = [mp.mpf(i) / 200 for i in range(1, 200)]
    assert all(f(x) >= -mp.mpf("1e-30") for x in xs)
