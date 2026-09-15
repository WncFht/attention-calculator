"""Edge-behavior tests: degenerate proof templates, float64 direction
pre-checks, and raw wire echo in rendering.

Ground truth: live zhuyidao.net probes (bench/data/edge-probes.jsonl k3/k3b
and img/img2 batches) plus bench/data/golden.jsonl — every asserted status,
parameter dict and equation string is a probed site response, byte-identical.
"""

from fractions import Fraction

import pytest

from attention_calculator import engine, solve
from attention_calculator.kernels import (
    beta_family,
    exp_family,
    gamma,
    log_family,
    quadlog,
    trig_q,
    trig_pi,
)
from attention_calculator.render import coef_tex, coerce_params, rat_tex, wire_pair

F = Fraction


def params_of(res):
    """Typed params dict from a prove() result, as the image endpoint sees it."""
    return coerce_params({k: str(res["parameters"][k]) for k in (
        "m", "n", "a_val", "b_val", "c_val", "u_val",
        "au_val", "bu_val", "cu_val")})


# ------------------------------------------------------------- render helpers


@pytest.mark.parametrize("v,expected", [
    ("1/1", "1"),                    # raw 'n/1' integerizes
    ("4/1", "4"),
    ("2/2", "\\dfrac{2}{2}"),        # den != 1 stays an unreduced \dfrac
    ("3140/1000", "\\dfrac{3140}{1000}"),
    ("\\frac{-9}{8}", "\\frac{-9}{8}"),   # macros echo verbatim
    (F(6, 4), "\\dfrac{3}{2}"),
    (F(4), "4"),
])
def test_rat_tex(v, expected):
    assert rat_tex(v) == expected


@pytest.mark.parametrize("v,expected", [
    ("1", ""),                       # canonical 1 vanishes: \pi not 1\pi
    (F(1), ""),
    ("1/1", "1"),                    # ...but '1/1' is not canonical -> 1\pi
    ("-1", "-1"),
    ("0", "0"),
    ("2/2", "\\dfrac{2}{2}"),
])
def test_coef_tex(v, expected):
    assert coef_tex(v) == expected


@pytest.mark.parametrize("v,expected", [
    ("4/2", (4, 2)),                 # plain text: signed, unreduced
    ("-6/4", (-6, 4)),
    ("\\dfrac{-6}{4}", (6, 4)),      # macro: minus inside braces dropped
    ("-\\dfrac{6}{4}", (6, 4)),      # ...and minus outside dropped too
    ("\\dfrac{6}{-4}", (6, 4)),
    ("5", (5, 1)),
    (F(6, 4), (3, 2)),               # Fraction passthrough stays reduced
])
def test_wire_pair(v, expected):
    assert wire_pair(v) == expected


def test_wire_pair_junk():
    with pytest.raises(ValueError):
        wire_pair("x")


# ------------------------------------- '<' degenerate template (cu_val == 2)


@pytest.mark.parametrize("kind,power,bound,expected", [
    ("gauss", "1", "5/4",
     dict(m=0, n=0, a_val="0", b_val="1/24", au_val="25", bu_val="0",
          cu_val="2", u_val="4")),
    ("gauss", "2", "3",
     dict(m=0, n=0, a_val="0", b_val="1/2", au_val="15", bu_val="0",
          cu_val="2", u_val="1")),
    ("gauss", "1/2", "1",
     dict(m=0, n=0, a_val="0", b_val="1/3", au_val="5", bu_val="0",
          cu_val="2", u_val="1")),
    ("gauss", "0", "1",
     dict(m=0, n=0, a_val="0", b_val="5/6", au_val="5", bu_val="0",
          cu_val="2", u_val="1")),
    ("varpi", "1", "15/4",
     dict(m=0, n=0, a_val="0", b_val="1/4", au_val="21", bu_val="0",
          cu_val="2", u_val="1")),
    ("varpi", "1/2", "2",
     dict(m=0, n=0, a_val="0", b_val="1/4", au_val="21", bu_val="0",
          cu_val="2", u_val="2")),
    # t == 0 still fires the template (the gate is t >= 0, b > 0)
    ("varpi", "0", "1/2",
     dict(m=0, n=0, a_val="0", b_val="1/2", au_val="0", bu_val="0",
          cu_val="2", u_val="1")),
])
def test_lt_degenerate_template(kind, power, bound, expected):
    """gauss/varpi '<': when the plain search exhausts, the site emits
    ∫ t·x^res(1-x)√(1-x⁴) + b with cu_val=2 before trying the transposed
    solve — gauss 1<5/4 picks this over an existing (5,1) transposed hit."""
    res = solve.prove(kind, power, "<", bound)
    p = res["parameters"]
    for k_, want in expected.items():
        assert str(p[k_]) == str(want), f"{k_}: {p[k_]} != {want}"
    assert res["solution"] == "a = 0, b = 0"


@pytest.mark.parametrize("kind,power,bound,a,b", [
    ("gauss", "1", "6/5", "6", "-6"),   # b == 0 boundary -> template declines,
    ("varpi", "1", "7/2", "21", "-21"), # falls back to a real solve
])
def test_lt_degenerate_boundary_declines(kind, power, bound, a, b):
    res = solve.prove(kind, power, "<", bound)
    assert res["parameters"]["cu_val"] == "0"
    assert res["solution"] == f"a = {a}, b = {b}"


@pytest.mark.parametrize("params,kind,power,comp,bound,equation", [
    # probed via /get_integral_image directly (img: batch)
    (dict(m=0, n=0, a_val="0", b_val="5/6", c_val="0", u_val="1",
          au_val="5", bu_val="0", cu_val="2"),
     "gauss", "0", "<", "1",
     r"G^{-1}-0 = \int_0^1 \left(5 x^{2} \cdot \left(1 - x\right)\right)"
     r"\sqrt{1-x^4} \mathrm{d} x+\dfrac{5}{6} > 0"),
    (dict(m=0, n=0, a_val="0", b_val="1/4", c_val="0", u_val="1",
          au_val="21", bu_val="0", cu_val="2"),
     "varpi", "1", "<", "15/4",
     r"\dfrac{15}{4} - \varpi = \int_0^1  21 x^{3} \cdot \left(1 - x\right)"
     r"\sqrt{1-x^4} \mathrm{d} x+\dfrac{1}{4} > 0"),
    # golden records
    (dict(m=0, n=0, a_val="0", b_val="1/4", c_val="0", u_val="2",
          au_val="21", bu_val="0", cu_val="2"),
     "varpi", "1/2", "<", "2",
     r"2 - \dfrac{1}{2}\varpi = \int_0^1  \frac{21 x^{3} \cdot \left(1 - x"
     r"\right)}{2}\sqrt{1-x^4} \mathrm{d} x+\dfrac{1}{4} > 0"),
    (dict(m=0, n=0, a_val="0", b_val="1/2", c_val="0", u_val="1",
          au_val="0", bu_val="0", cu_val="2"),
     "varpi", "0", "<", "1/2",
     r"\dfrac{1}{2} - 0\varpi = \int_0^1  0\sqrt{1-x^4} \mathrm{d} x"
     r"+\dfrac{1}{2} > 0"),
    (dict(m=0, n=0, a_val="0", b_val="1/3", c_val="0", u_val="1",
          au_val="5", bu_val="0", cu_val="2"),
     "gauss", "1/2", "<", "1",
     r"G^{-1}-\dfrac{1}{2} = \int_0^1 \left(5 x^{2} \cdot \left(1 - x\right)"
     r"\right)\sqrt{1-x^4} \mathrm{d} x+\dfrac{1}{3} > 0"),
    (dict(m=0, n=0, a_val="0", b_val="5/12", c_val="0", u_val="2",
          au_val="5", bu_val="0", cu_val="2"),
     "gauss", "0", "<", "1/2",
     r"\dfrac{1}{2}G^{-1}-0 = \int_0^1 \left(\frac{5 x^{2} \cdot \left(1 - x"
     r"\right)}{2}\right)\sqrt{1-x^4} \mathrm{d} x+\dfrac{5}{12} > 0"),
])
def test_lt_degenerate_render(params, kind, power, comp, bound, equation):
    typed = coerce_params({k: str(v) for k, v in params.items()})
    assert beta_family.render_equation(typed, kind, power, comp, bound) == equation


# -------------------------------------------- float64 direction pre-checks


@pytest.mark.parametrize("kind,power,comp,bound", [
    ("cosh_q", "0", "<", "1/2"),   # cosh(0)=1, claim 1<1/2 false
    ("e_q", "0", ">", "2"),        # e^0=1 > 2 false
    ("arctan_q", "0", ">", "2"),   # atan(0)=0 > 2 false
    ("pi_n", "0", ">", "2"),       # pi^0=1 > 2 false
    ("gamma", "0", ">", "1"),      # 0*gamma > 1 false — dies before bound/0
    ("gamma", "3", "<", "1"),      # 3*gamma < 1 false
])
def test_float64_precheck_wrong_direction(kind, power, comp, bound):
    """Strictly-false claims in the crashing kernels are rejected by the
    site's float64 pre-check as 方向反了 before the kernel machinery runs."""
    with pytest.raises(engine.WrongDirection):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize("kind,power,comp,bound,exc", [
    # equality/truth proceeds into the kernel, which then crashes (-> 500)
    ("sinh_q", "0", "<", "0", ZeroDivisionError),   # sinh(0)=0, 0<0 not strict
    ("cosh_q", "0", ">", "1/2", ZeroDivisionError), # 1 > 1/2 true -> 1/q crash
    ("coth_q", "0", ">", "1", ZeroDivisionError),   # 1/tanh(0) inside the check
    ("sinh_q", "711", ">", "1", OverflowError),     # math.sinh overflow
    ("cosh_q", "1000000000000000/1", ">", "0", OverflowError),
    ("e_q", "800", ">", "1", OverflowError),        # math.exp overflow
    ("arccot_q", "0", ">", "2", ZeroDivisionError), # atan(1/0)
    ("pi_n", "11", ">", "1", KeyError),             # exponent outside [1,10]
    ("pi_n", "-1", "<", "1", KeyError),             # BETA_PI[-1]
    ("gamma", "0", "<", "1", ZeroDivisionError),    # bound/abs(power)
])
def test_float64_precheck_passthrough_crash(kind, power, comp, bound, exc):
    """True/equal claims pass the pre-check and hit the kernel's own crash —
    the site's 500 surface."""
    with pytest.raises(exc):
        solve.prove(kind, power, comp, bound)


def test_tanh_huge_power_exhausts():
    """tanh(1e15) == 1.0 in float64: claim 1>1 is not strictly false, so the
    search runs and exhausts -> NoSolution (site: 404 未找到)."""
    with pytest.raises(engine.NoSolution):
        solve.prove("tanh_q", "1000000000000000/1", ">", "1")


def test_sinh_700_solves():
    """sinh(700) stays finite in float64 — a real solve comes back."""
    res = solve.prove("sinh_q", "700", ">", "1")
    assert res["parameters"]["m"] == 0 and res["parameters"]["n"] == 0
    assert res["solution"] == "a = 245699/350, b = 488600, c= -489300"


# --------------------------------------------------------- pi_n raw (num, den)


@pytest.mark.parametrize("power,equation", [
    # raw pair drives parity, the outer ^den, and the ln exponent — all
    # unreduced; '\frac' macros parse unsigned
    ("-3/4",
     r"\left(\dfrac{311}{10}\right)^{4} - \left(\pi^\dfrac{-3}{4}\right)^{4}"
     r" = \int_0^1 \frac{x^{4} \cdot \left(1925 x^{2} + 197\right)}"
     r"{108 x^{2} + 108}(\ln(1/x))^-4 \mathrm{d} x > 0"),
    ("5/3",
     r"\left(\dfrac{311}{10}\right)^{3} - \left(\pi^\dfrac{5}{3}\right)^{3}"
     r" = \int_0^1 \frac{x^{4} \cdot \left(1925 x^{2} + 197\right)}"
     r"{108 x^{2} + 108}(\ln(1/x))^4 \mathrm{d} x > 0"),
    ("-5/3",
     r"\left(\dfrac{311}{10}\right)^{3} - \left(\pi^\dfrac{-5}{3}\right)^{3}"
     r" = \int_0^1 \frac{x^{4} \cdot \left(1925 x^{2} + 197\right)}"
     r"{108 x^{2} + 108}(\ln(1/x))^-6 \mathrm{d} x > 0"),
    ("2/4",  # unreduced: num 2 even -> odd powers, den 4 outer, ln^1
     r"\left(\dfrac{311}{10}\right)^{4} - \left(\pi^\dfrac{2}{4}\right)^{4}"
     r" = \int_0^1 \frac{x^{5} \cdot \left(1925 x^{2} + 197\right)}"
     r"{108 x^{2} + 108}(\ln(1/x)) \mathrm{d} x > 0"),
    ("1/1",  # den 1: plain LHS, pi^1 (exponent prints even at 1), no ln tail
     r"\dfrac{311}{10} - \pi^1 = \int_0^1 \frac{x^{4} \cdot "
     r"\left(1925 x^{2} + 197\right)}{108 x^{2} + 108} \mathrm{d} x > 0"),
    ("0",    # pi^0 and a negative ln exponent
     r"\dfrac{311}{10} - \pi^0 = \int_0^1 \frac{x^{5} \cdot "
     r"\left(1925 x^{2} + 197\right)}{108 x^{2} + 108}(\ln(1/x))^-1"
     r" \mathrm{d} x > 0"),
])
def test_pin_raw_pair_render(power, equation):
    params = coerce_params(dict(m="2", n="0", a_val="197/108",
                                b_val="1925/108", c_val="0", u_val="108",
                                au_val="197", bu_val="1925", cu_val="0"))
    assert quadlog.render_equation(params, "pi_n", power, "<", "311/10") == equation


# --------------------------------------------------------------- gamma render


@pytest.mark.parametrize("params,power,comp,bound,equation", [
    # coef folds into the kernel bracket and scales the N=0 constant tail
    (dict(m=0, n=0, a_val="7/3", b_val="0", c_val="1", u_val="0",
          au_val="0", bu_val="0", cu_val="0"), "2", ">", "7/9",
     r"2\gamma - \dfrac{7}{9} = \int_0^1 \left[2\left(\dfrac{1}{1-x}"
     r"+\dfrac{1}{\ln(x)}-\dfrac{1}{2}\right)+\dfrac{14}{3}\right]"
     r" \mathrm{d} x > 0"),
    (dict(m=0, n=0, a_val="7/3", b_val="0", c_val="2", u_val="0",
          au_val="0", bu_val="0", cu_val="0"), "3", "<", "7/9",
     r"\dfrac{7}{9} - 3\gamma =\int_0^1 \left[3\left(\dfrac{2-x}{2}"
     r"-\dfrac{1}{1-x}-\dfrac{1}{\ln(x)}\right)+7\right] \mathrm{d} x > 0"),
    (dict(m=0, n=0, a_val="1/4", b_val="0", c_val="1", u_val="0",
          au_val="0", bu_val="0", cu_val="0"), "3/2", ">", "1/2",
     r"\dfrac{3}{2}\gamma - \dfrac{1}{2} = \int_0^1 \left[\dfrac{3}{2}\left("
     r"\dfrac{1}{1-x}+\dfrac{1}{\ln(x)}-\dfrac{1}{2}\right)+\dfrac{3}{8}"
     r"\right] \mathrm{d} x > 0"),
    # cf·x^N prefix goes through sympy ('1/2' -> \frac{x^{5}}{2}) and the
    # sub-fraction scalar cf/u is printed reduced (24 -> 48)
    (dict(m=0, n=2, a_val="1/8", b_val="7/6", c_val="0", u_val="24",
          au_val="3", bu_val="28", cu_val="2"), "1/2", "<", "1/3",
     r"\dfrac{1}{3} - \dfrac{1}{2}\gamma =\int_0^1 \left[\frac{x^{2}}{2}"
     r"\left(\dfrac{2-x}{2}-\dfrac{1}{1-x}-\dfrac{1}{\ln(x)}\right)"
     r"+\frac{\left(1 - x\right)^{2} \cdot \left(28 x + 3\right)}"
     r"{48 \left(2 x + 1\right)^{2}}\right] \mathrm{d} x > 0"),
    # cf/u reduction cancelling fully (3/3 -> suppressed denominator scalar)
    (dict(m=1, n=2, a_val="0", b_val="4/3", c_val="0", u_val="3",
          au_val="0", bu_val="4", cu_val="2"), "3", ">", "5/3",
     r"3\gamma - \dfrac{5}{3} = \int_0^1 \left[3 x^{2}\left(\dfrac{1}{1-x}"
     r"+\dfrac{1}{\ln(x)}-\dfrac{1}{2}\right)+\frac{4 x^{2} \left(1 - x"
     r"\right)^{2}}{\left(2 x + 1\right)^{2}}\right] \mathrm{d} x > 0"),
    # u_val=0 with a nonzero polynomial: the sub-fraction collapses to
    # sympy's \tilde{\infty} times the numerator factors
    (dict(m=2, n=4, a_val="937/56", b_val="4615/42", c_val="0", u_val="0",
          au_val="2811", bu_val="18460", cu_val="5"), "1", "<", "3/5",
     r"\dfrac{3}{5} - \gamma =\int_0^1 \left[x^{5}\left(\dfrac{2-x}{2}"
     r"-\dfrac{1}{1-x}-\dfrac{1}{\ln(x)}\right)+\tilde{\infty} x^{2} "
     r"\left(1 - x\right)^{4} \cdot \left(18460 x + 2811\right)\right]"
     r" \mathrm{d} x > 0"),
])
def test_gamma_render(params, power, comp, bound, equation):
    typed = coerce_params({k: str(v) for k, v in params.items()})
    assert gamma.render_equation(typed, "gamma", power, comp, bound) == equation


# -------------------------------------------- ln_q raw-pair denominator


@pytest.mark.parametrize("power,equation", [
    # '4/2': (d,e) = (4-2, 2) -> (2x+2)^6; the ln factor echoes unreduced
    ("4/2",
     r"\dfrac{224}{125} - \ln\dfrac{4}{2} = \int_0^1 \frac{64 x^{6} "
     r"\left(1 - x\right)^{4} \cdot \left(32663 - 29785 x\right)}"
     r"{25 \left(2 x + 2\right)^{6}} \mathrm{d} x > 0"),
    # '-6/4': (d,e) = (-10, 4) -> sympy prints (4 - 10x)^6
    ("-6/4",
     r"\dfrac{224}{125} - \ln\dfrac{-6}{4} = \int_0^1 \frac{4096 x^{6} "
     r"\left(1 - x\right)^{4} \cdot \left(32663 - 29785 x\right)}"
     r"{25 \left(4 - 10 x\right)^{6}} \mathrm{d} x > 0"),
    # '4/-2': signed plain text -> (d,e) = (6, -2) -> (6x-2)^6
    ("4/-2",
     r"\dfrac{224}{125} - \ln\dfrac{4}{-2} = \int_0^1 \frac{64 x^{6} "
     r"\left(1 - x\right)^{4} \cdot \left(32663 - 29785 x\right)}"
     r"{25 \left(6 x - 2\right)^{6}} \mathrm{d} x > 0"),
    # '-2': (d,e) = (-3, 1) -> (1-3x)^6, numerator scale 1
    ("-2",
     r"\dfrac{224}{125} - \ln-2 = \int_0^1 \frac{x^{6} "
     r"\left(1 - x\right)^{4} \cdot \left(32663 - 29785 x\right)}"
     r"{25 \left(1 - 3 x\right)^{6}} \mathrm{d} x > 0"),
])
def test_lnq_wire_pair_denominator(power, equation):
    params = coerce_params(dict(m="6", n="4", a_val="32663/25",
                                b_val="-5957/5", c_val="0", u_val="25",
                                au_val="32663", bu_val="-29785", cu_val="0"))
    assert log_family.render_equation(
        params, "ln_q", power, "<", "224/125") == equation


# ------------------------------------------------------- e_q / trig_q raw echo


@pytest.mark.parametrize("power,equation", [
    ("-4",
     r"e^-4 - \dfrac{22}{5} = \int_0^1 \frac{x \left(1 - x\right)^{2} \cdot "
     r"\left(63 x + 12\right) e^{- 4 x}}{80} \mathrm{d} x > 0"),
    # the const echoes the macro verbatim while the kernel exponent takes the
    # macro's unsigned digit groups: '-\dfrac{6}{4}' -> e^{3x/2}
    ("-\\dfrac{6}{4}",
     r"e^-\dfrac{6}{4} - \dfrac{22}{5} = \int_0^1 \frac{x \left(1 - x\right)"
     r"^{2} \cdot \left(63 x + 12\right) e^{\frac{3 x}{2}}}{80}"
     r" \mathrm{d} x > 0"),
])
def test_eq_raw_render(power, equation):
    params = coerce_params(dict(m="1", n="2", a_val="3/20", b_val="63/80",
                                c_val="0", u_val="80", au_val="12",
                                bu_val="63", cu_val="0"))
    assert exp_family.render_equation(params, "e_q", power, ">", "22/5") == equation


def test_eq_junk_coef_crashes():
    params = coerce_params(dict(m="1", n="2", a_val="3/20", b_val="63/80",
                                c_val="0", u_val="80", au_val="12",
                                bu_val="63", cu_val="0"))
    with pytest.raises(ValueError):
        exp_family.render_equation(params, "e_q", "x", ">", "22/5")


@pytest.mark.parametrize("power,equation", [
    # tan divisor takes the macro's unsigned digits: \cos(6/4); the sin
    # factor gets |q| = 3/2
    ("\\dfrac{-6}{4}",
     r"\tan\dfrac{-6}{4} - 1 = \int_0^1 \dfrac{1}{\cos(6/4)} x "
     r"\sin{\left(\frac{3 x}{2} \right)} \mathrm{d} x > 0"),
    ("\\dfrac{6}{4}",
     r"\tan\dfrac{6}{4} - 1 = \int_0^1 \dfrac{1}{\cos(6/4)} x "
     r"\sin{\left(\frac{3 x}{2} \right)} \mathrm{d} x > 0"),
    ("1",
     r"\tan1 - 1 = \int_0^1 \dfrac{1}{\cos(1)} x \sin{\left(x \right)}"
     r" \mathrm{d} x > 0"),
])
def test_tanq_raw_render(power, equation):
    params = coerce_params(dict(m="0", n="0", a_val="0", b_val="1",
                                c_val="0", u_val="1", au_val="0",
                                bu_val="1", cu_val="0"))
    assert trig_q.render_equation(params, "tan_q", power, ">", "1") == equation


def test_tanq_junk_coef_crashes():
    params = coerce_params(dict(m="0", n="0", a_val="0", b_val="1",
                                c_val="0", u_val="1", au_val="0",
                                bu_val="1", cu_val="0"))
    with pytest.raises(ValueError):
        trig_q.render_equation(params, "tan_q", "x", ">", "1")


# --------------------------------------------- coef_tex on solved equations


@pytest.mark.parametrize("kind,comp,power,bound,equation", [
    # '2/2' is not canonical 1 -> \dfrac{2}{2}\pi
    ("pi", "<", "2/2", "22/7",
     r"\dfrac{22}{7} - \dfrac{2}{2}\pi = \int_0^1 \frac{x^{6} "
     r"\left(1 - x^{2}\right)^{3} \cdot \left(47 - 13 x^{2}\right)}"
     r"{120 x^{2} + 120} \mathrm{d} x > 0"),
    ("pi", "<", "-1", "22/7",
     r"\dfrac{22}{7} - -1\pi = \int_0^1 \frac{x^{6} "
     r"\left(1 - x^{2}\right)^{3} \cdot \left(47 - 13 x^{2}\right)}"
     r"{120 x^{2} + 120} \mathrm{d} x > 0"),
    ("pi", "<", "0", "22/7",
     r"\dfrac{22}{7} - 0\pi = \int_0^1 \frac{x^{6} "
     r"\left(1 - x^{2}\right)^{3} \cdot \left(47 - 13 x^{2}\right)}"
     r"{120 x^{2} + 120} \mathrm{d} x > 0"),
    # bound '1/1' integerizes in the constant slot
    ("pi", "<", "1", "1/1",
     r"1 - \pi = \int_0^1 \frac{x^{6} \left(1 - x^{2}\right)^{3} \cdot "
     r"\left(47 - 13 x^{2}\right)}{120 x^{2} + 120} \mathrm{d} x > 0"),
    # varpi '>' shows power - bound·varpi^{-1}: the '4/2' stays unreduced on
    # the inverse side and the '1/1' coef integerizes on the constant side
    ("varpi", ">", "1/1", "4/2",
     r"1 - \dfrac{4}{2}\varpi^{-1} = \int_0^1  \frac{x^{13} \cdot "
     r"\left(47 - 13 x^{4}\right)}{120}\dfrac{(1-x)}{\pi\sqrt{1-x^4}}"
     r" \mathrm{d} x > 0"),
    # gauss '<' shows bound·G^{-1} - power, no spaces around '-'
    ("gauss", "<", "1/1", "4/2",
     r"\dfrac{4}{2}G^{-1}-1 = \int_0^1 \frac{x^{14} \cdot "
     r"\left(47 - 13 x^{4}\right)}{120}\dfrac{(1-x)}{\sqrt{1-x^4}}"
     r" \mathrm{d} x > 0"),
])
def test_coef_and_bound_echo(kind, comp, power, bound, equation):
    """Params are a real (3,3) solve's; only the display slots vary."""
    params = coerce_params(dict(m="3", n="3", a_val="47/120",
                                b_val="-13/120", c_val="0", u_val="120",
                                au_val="47", bu_val="-13", cu_val="0"))
    module = {"pi": quadlog, "varpi": beta_family,
              "gauss": beta_family}[kind]
    assert module.render_equation(params, kind, power, comp, bound) == equation


# ------------------------------------------------ trig_pi cdot on \left(\d…


def test_trig_pi_cdot_rule():
    r"""' \cdot ' appears only before a \left(<digit>…\right) group — a
    negative-leading group keeps a plain space (probed: ' \left(- 1465926').
    A digit-leading group after '}' does take the cdot."""
    params = coerce_params(dict(m="1", n="3", a_val="1905533/18750000",
                                b_val="-244321/3125000", c_val="3/5",
                                u_val="18750000",
                                au_val="-100000000000000000000000",
                                bu_val="-1465926", cu_val="0"))
    eq = trig_pi.render_equation(
        params, "sin_q_degree", "36", ">", "587/1000")
    assert eq == (
        r"\sin\left(36^\circ\right) - \dfrac{587}{1000} = \int_0^{\pi/2} "
        r"\frac{\left(1 - \sin{\left(x \right)}\right)^{3} \left(- 1465926 "
        r"\sin{\left(x \right)} - 100000000000000000000000\right) "
        r"\sin{\left(\frac{3 x}{5} \right)} \sin{\left(x \right)}}"
        r"{18750000} \mathrm{d} x > 0")
