"""Tests for the Padé-interpolation second prover (pade.py, plan doc W4).

The engine proves ``ln q ⋚ r`` and ``arctan q ⋚ r`` by interpolating strict
one-sided Padé bounds whose error functions are manifestly non-negative
rational functions.  Two properties are pinned: every emitted certificate
re-verifies exactly over QQ (verify_cert), and Padé reaches bounds the
site's m,n <= 10 search cannot (the point of the second engine -- the
site-failure assertions here compare against site-mode solve.prove).
"""

from fractions import Fraction
from math import comb

import mpmath as mp
import pytest

from attention_calculator import pade, solve
from attention_calculator.engine import NoSolution, WrongDirection

mp.mp.dps = 60


def rat_bound(c: mp.mpf, digits: int, above: bool) -> Fraction:
    """Floor/ceiling rational approximant of c with denominator 10^digits."""
    D = 10**digits
    k = int(c * D)
    b = Fraction(k, D)
    if above and mp.mpf(b.numerator) / b.denominator <= c:
        b = Fraction(k + 1, D)
    if not above and mp.mpf(b.numerator) / b.denominator >= c:
        b = Fraction(k - 1, D)
    return b


def prove(kind, q, comp, bound, **kw):
    """pade.prove + mandatory verify_cert on every emitted certificate."""
    cert = pade.prove(kind, q, comp, bound, **kw)
    assert cert is not None, f"no certificate for {kind} {q} {comp} {bound}"
    assert pade.verify_cert(cert)
    return cert


# ------------------------------------------------------- article examples


def test_article_ln_lower_interpolation_weights():
    # 进阶教程算例 ln(8/5) > 47/100：a=4499/4500 挂在首个越界项 l_3 上，
    # 与最弱下界 l_1 配对 —— cert 的 b 即文章 a
    cert = prove("ln_q", Fraction(8, 5), ">", Fraction(47, 100))
    assert (cert["n"], cert["m"]) == (1, 3)
    assert cert["b"] == Fraction(4499, 4500)
    assert cert["a"] == Fraction(1, 4500)
    assert cert["resid"] == 0


def test_article_ln_upper():
    # 算例 ln(9/2) < 11/7：u_0 = x 与首个落到界下的 u_2 配对
    cert = prove("ln_q", Fraction(9, 2), "<", Fraction(11, 7))
    assert (cert["n"], cert["m"]) == (0, 2)
    assert cert["a"] + cert["b"] == 1
    assert cert["resid"] == 0


def test_article_arctan_upper():
    # 算例 arctan 2 < 6/5：u_0 = 2 与首个越界 u_2 配对，1749/1918 挂 u_2
    cert = prove("arctan_q", Fraction(2), "<", Fraction(6, 5))
    assert (cert["n"], cert["m"]) == (0, 2)
    assert cert["b"] == Fraction(1749, 1918)


# ------------------------------------------------------------ round trips


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("ln_q", "2", ">", "2/3"),  # equality edge: l_1(1) == 2/3
        ("ln_q", "2", "<", "7/10"),
        ("ln_q", "3", "<", "11/10"),  # golden case, site also proves
        ("ln_q", "3/2", ">", "2/5"),
        ("ln_q", "10", ">", "23/10"),
        ("arctan_q", "1", ">", "157/200"),
        ("arctan_q", "1", "<", "3927/5000"),
        ("arctan_q", "1/2", ">", "11/24"),
        ("arctan_q", "3", "<", "5/4"),
    ],
)
def test_certificates_verify(kind, power, comp, bound):
    prove(kind, Fraction(power), comp, Fraction(bound))


def test_residual_path_single_term():
    # ln 2 > 0: already l_1(1) = 2/3 > 0 -- no below-p approximant exists,
    # so the cert is a single term plus the positive residual l_1 - 0
    cert = prove("ln_q", Fraction(2), ">", Fraction(0))
    assert cert["a"] == 0 and cert["b"] == 1
    assert cert["resid"] == Fraction(2, 3)


# ----------------------------------------------- beyond the site budget
#
# Bounds the site's m,n <= 10 search cannot prove (site-mode solve.prove
# raises); the gaps were measured by sweeping solve.prove against the
# constants.  Padé crosses them at the index recorded in the notes.

SITE_FAIL_CASES = [
    # ln 2 > floor at 1e-20 (gap ~7e-21): site reach is only ~1e-15
    ("ln_q", "2", ">", rat_bound(mp.log(2), 20, False), 14),
    # ln 3 > floor at 1e-18 (gap ~4e-19): site reach ~1e-11
    ("ln_q", "3", ">", rat_bound(mp.log(3), 18, False), 17),
    # ln(3/2) > floor at 1e-25 (gap ~2e-26)
    ("ln_q", "3/2", ">", rat_bound(mp.log(mp.mpf(3) / 2), 25, False), 13),
    # ln(9/2) < 1.5040773968 (gap ~3e-10): site reach ~1e-8 but this
    # particular bound is below it... verified empirically -> NoSolution
    ("ln_q", "9/2", "<", Fraction(940048373, 625000000), 12),
    # arctan 1 = pi/4 < ceil at 1e-20 (gap ~4e-21): site reach ~1e-9
    ("arctan_q", "1", "<", rat_bound(mp.pi / 4, 20, True), 13),
    # arctan 2 < 692/625 = 1.1072 (gap ~5e-5): site reach only ~1e-3 here
    ("arctan_q", "2", "<", Fraction(692, 625), 5),
    # arctan(1/2) > floor at 1e-19 (gap ~2e-20)
    ("arctan_q", "1/2", ">", rat_bound(mp.atan(mp.mpf(1) / 2), 19, False), 8),
]


@pytest.mark.parametrize(("kind", "power", "comp", "bound", "cross"), SITE_FAIL_CASES)
def test_beyond_site_budget(kind, power, comp, bound, cross):
    cert = prove(kind, Fraction(power), comp, bound)
    assert cert["m"] == cross
    with pytest.raises((NoSolution, WrongDirection)):
        solve.prove(kind, power, comp, str(bound))


def test_site_also_proves_loose_bound():
    # sanity: on loose bounds both engines agree
    cert = prove("ln_q", Fraction(2), "<", Fraction(7, 10))
    assert cert["resid"] == 0
    solve.prove("ln_q", "2", "<", "7/10")


# ------------------------------------------------------------ false claims


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("ln_q", "2", ">", "7/10"),  # ln 2 < 7/10: vetoed by u_8
        ("ln_q", "2", "<", "3/5"),  # ln 2 > 3/5
        ("ln_q", "3", "<", "1"),  # ln 3 > 1
        ("arctan_q", "1", "<", "7/10"),  # pi/4 > 7/10
        ("arctan_q", "2", ">", "6/5"),  # arctan 2 < 6/5
    ],
)
def test_false_claim_returns_none(kind, power, comp, bound):
    assert pade.prove(kind, Fraction(power), comp, Fraction(bound)) is None


def test_near_miss_false_claim_exhausts_budget():
    # arctan 2 < (arctan 2 - 2e-8): false by a hair; no approximant crosses
    # within a small budget -> honest None
    bound = rat_bound(mp.atan(2), 7, False)
    assert pade.prove("arctan_q", Fraction(2), "<", bound, max_n=8) is None


# ------------------------------------------------------------------ domain


def test_out_of_domain_returns_none():
    assert pade.prove("ln_q", Fraction(1), ">", Fraction(0)) is None
    assert pade.prove("ln_q", Fraction(1, 2), "<", Fraction(0)) is None
    assert pade.prove("arctan_q", Fraction(-1), "<", Fraction(0)) is None


def test_uncovered_kind_raises():
    with pytest.raises(ValueError):
        pade.prove("zeta3", Fraction(1), ">", Fraction(1))


# -------------------------------------------------------- certificate check


def test_tampered_certificates_fail():
    cert = pade.prove("ln_q", Fraction(8, 5), ">", Fraction(47, 100))

    bad = dict(cert, a=cert["a"] + 1)
    assert not pade.verify_cert(bad)  # a + b != 1
    bad = dict(cert, resid=Fraction(1))
    assert not pade.verify_cert(bad)  # resid identity
    bad = dict(cert, m=cert["m"] + 1)
    assert not pade.verify_cert(bad)  # wrong approximant
    bad = dict(cert, comp="<")
    assert not pade.verify_cert(bad)  # wrong direction
    bad = dict(cert, q=Fraction(1))
    assert not pade.verify_cert(bad)  # wrong point
    bad = dict(cert, serr=[])
    assert not pade.verify_cert(bad)  # missing serr


def test_cert_is_exact_rationals():
    cert = pade.prove("ln_q", Fraction(3), ">", Fraction(1))
    for key in ("q", "p", "a", "b", "resid"):
        assert isinstance(cert[key], Fraction)
    for term in cert["serr"]:
        assert isinstance(term["c"], Fraction)
        assert all(isinstance(v, Fraction) for v in term["den"])


# --------------------------------------------------- error function forms


def test_error_function_closed_form():
    # ln l_n: s_n = x^{2n} / (C(2n,n)^2 (1+x) D_n^2)  -- observed identity
    for n in (1, 2, 5):
        P, Q = pade.approx("ln", n, n)
        c, e = pade.err_form(P, Q, "ln", upper=False)
        assert (c, e) == (Fraction(1, comb(2 * n, n) ** 2), 2 * n)
    # ln u_n: t_n = x^{2n+1} / (C(2n+1,n)^2 (1+x) D~_n^2)
    for n in (0, 1, 4):
        P, Q = pade.approx("ln", n + 1, n)
        c, e = pade.err_form(P, Q, "ln", upper=True)
        assert (c, e) == (Fraction(1, comb(2 * n + 1, n) ** 2), 2 * n + 1)


def test_error_function_exponents_atan():
    # arctan l_n = [2n/2n] -> e = 4n; u_n = [2n+1/2n+1] -> e = 4n+2
    for n in (1, 2, 3):
        P, Q = pade.approx("atan", 2 * n, 2 * n)
        c, e = pade.err_form(P, Q, "atan", upper=False)
        assert c > 0 and e == 4 * n
        P, Q = pade.approx("atan", 2 * n + 1, 2 * n + 1)
        c, e = pade.err_form(P, Q, "atan", upper=True)
        assert c > 0 and e == 4 * n + 2


def test_denominators_pole_free_on_positive_axis():
    # all-nonneg coefficients (with Q(0)=1) => no poles on x > 0; pinned
    # as the structural fact verify_cert relies on
    for n in range(1, 21):
        for func, L, M in (("ln", n, n), ("ln", n + 1, n), ("atan", 2 * n, 2 * n)):
            _, Q = pade.approx(func, L, M)
            assert Q[0] > 0 and all(v >= 0 for v in Q)
