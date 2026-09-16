"""Tests for the Euler--Maclaurin gamma second prover (euler_gamma.py, W7).

The engine proves ``q * gamma ⋚ r`` through a certified rational enclosure
of gamma: EM expansion at N = 2^t plus a two-sided ln-2 child bound.  The
pinned properties: every emitted certificate re-verifies exactly over QQ
(verify_cert replays H_N, Bernoulli numbers, the tail bound, the Bernstein
sign certificate, and both child certs), and the prover reaches gaps the
site kernel's (m, n) budget cannot (kernel fails already at 1e-8 -- those
site-mode failures are asserted alongside).
"""

import json
from fractions import Fraction

import pytest
from mpmath import mp

from attention_calculator import certificate, euler_gamma, solve
from attention_calculator.engine import NoSolution, WrongDirection


@pytest.fixture(autouse=True, scope="module")
def module_dps():
    """本模块数值校验的 mpmath 精度。"""
    with mp.workdps(80):
        yield


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


def prove(power, comp, bound, **kw):
    """euler_gamma.prove + mandatory verify_cert on every emitted cert."""
    cert = euler_gamma.prove("gamma", Fraction(power), comp, Fraction(bound), **kw)
    assert cert is not None, f"no certificate for gamma {power} {comp} {bound}"
    assert euler_gamma.verify_cert(cert)
    return cert


# ------------------------------------------------------------- certificates


@pytest.mark.parametrize(
    ("comp", "digits"),
    [
        (">", 8),  # gap ~5e-9
        ("<", 8),
        (">", 12),  # gap ~5e-13
        ("<", 12),
        (">", 20),  # gap ~6e-21
        ("<", 20),
        (">", 35),  # gap ~1e-36
    ],
)
def test_certificates_verify_beyond_kernel(comp, digits):
    bound = rat_bound(mp.euler, digits, above=(comp == "<"))
    cert = prove(1, comp, bound)
    # certified enclosure actually brackets gamma (numerical sanity)
    lo, hi = Fraction(cert["lo"]), Fraction(cert["hi"])
    gm = mp.mpf(mp.euler)
    assert mp.mpf(lo.numerator) / lo.denominator < gm < mp.mpf(hi.numerator) / hi.denominator
    # the site kernel cannot reach these (fails instantly at gaps <= 1e-8)
    with pytest.raises((NoSolution, WrongDirection)):
        solve.prove("gamma", "1", comp, str(bound))


def test_loose_bounds_both_directions():
    prove(1, ">", Fraction(1, 2))
    prove(1, "<", Fraction(3, 5))
    prove(1, ">", Fraction(0))  # positivity discharge through the same path
    prove(1, "<", Fraction(10))


def test_coefficient_and_sign_flip():
    # 2g < 6/5 normalizes to g < 3/5; -g < -577/1000 normalizes to g > 577/1000
    cert = prove(2, "<", Fraction(6, 5))
    assert cert["q"] == "2" and cert["p"] == "6/5"
    cert = prove(-1, "<", Fraction(-577, 1000))
    assert cert["comp"] == "<"  # literal request kept; flip is internal


def test_enclosure_is_tight():
    # the recorded interval is the EM enclosure -- width comparable to the gap
    bound = rat_bound(mp.euler, 14, above=False)
    cert = prove(1, ">", bound)
    w = Fraction(cert["hi"]) - Fraction(cert["lo"])
    assert 0 < w < Fraction(1, 10**10)


# ------------------------------------------------------------- false claims


@pytest.mark.parametrize(
    ("comp", "bound"),
    [
        (">", "3/5"),  # gamma < 0.578
        ("<", "1/2"),
        (">", "5773/10000"),  # false by ~4e-5
        ("<", "577215664/1000000000"),  # false by ~1e-9
    ],
)
def test_false_claim_returns_none(comp, bound):
    assert euler_gamma.prove("gamma", Fraction(1), comp, Fraction(bound)) is None


def test_zero_coefficient_rejected():
    with pytest.raises(ValueError):
        euler_gamma.prove("gamma", Fraction(0), ">", Fraction(-1))


def test_uncovered_kind_raises():
    with pytest.raises(ValueError):
        euler_gamma.prove("ln_q", Fraction(2), ">", Fraction(0))


def test_honest_none_past_cap():
    # true claim with the ln-2 child budget cut: a 1e-40 slack on ln 2
    # cannot be certified at pade index <= 4 (its reach there is ~1e-6)
    bound = rat_bound(mp.euler, 40, above=False)
    assert euler_gamma.prove("gamma", Fraction(1), ">", bound, max_n=4) is None


# ------------------------------------------------------------- tamper-proof


@pytest.fixture
def cert():
    return euler_gamma.prove("gamma", Fraction(1), ">", Fraction(57721566, 10**8))


def test_tampered_certificates_fail(cert):
    bad = dict(cert, lo=str(Fraction(cert["lo"]) + 1))
    assert not euler_gamma.verify_cert(bad)  # widened interval
    bad = dict(cert, hi=str(Fraction(cert["hi"]) - 1))
    assert not euler_gamma.verify_cert(bad)
    bad = dict(cert, n=cert["n"] * 2)
    assert not euler_gamma.verify_cert(bad)  # wrong N
    bad = dict(cert, j=cert["j"] + 2)
    assert not euler_gamma.verify_cert(bad)  # wrong term count
    bad = dict(cert, tail=str(Fraction(cert["tail"]) * 2))
    assert not euler_gamma.verify_cert(bad)  # wrong tail bound
    bad = dict(cert, comp="<")
    assert not euler_gamma.verify_cert(bad)  # wrong direction
    bad = dict(cert, p="3/5")
    assert not euler_gamma.verify_cert(bad)  # wrong claim
    bad = dict(cert, children=list(reversed(cert["children"])))
    assert not euler_gamma.verify_cert(bad)  # swapped ln children


def test_forged_child_fails(cert):
    # a valid pade cert for a different ln-2 bound does not match expect
    forged = euler_gamma.ln_cert("<", Fraction(8, 10), 30)
    bad = dict(cert, children=[cert["children"][0], forged])
    assert not euler_gamma.verify_cert(bad)
    # expect table disagreeing with the child cert
    bad = dict(cert, expect=[dict(e, bound="1/2") for e in cert["expect"]])
    assert not euler_gamma.verify_cert(bad)
    # a child cert for a different kind entirely
    other = euler_gamma.prove("gamma", Fraction(1), ">", Fraction(1, 2))
    bad = dict(cert, children=[cert["children"][0], other])
    assert not euler_gamma.verify_cert(bad)


def test_missing_and_extra_fields_fail(cert):
    for key in ("n", "j", "tail", "lo", "hi", "expect", "children"):
        bad = dict(cert)
        del bad[key]
        assert not euler_gamma.verify_cert(bad)
    bad = dict(cert, children=cert["children"][:1])
    assert not euler_gamma.verify_cert(bad)


def test_n_must_be_power_of_two(cert):
    bad = dict(cert, n=6)
    assert not euler_gamma.verify_cert(bad)


# ------------------------------------------------------------- wire plumbing


def test_json_round_trip():
    cert = prove(1, ">", Fraction(5772156649, 10**10))
    wire = json.loads(json.dumps(cert))
    assert euler_gamma.verify_cert(euler_gamma.cert_parse(wire))
    assert euler_gamma.verify_cert(wire)  # wire form verifies directly
    assert euler_gamma.cert_jsonable(euler_gamma.cert_parse(wire)) == wire


def test_end_to_end_wired():
    """solve.prove(exact=True) -> certificate.verify_cert round trip.

    The kernel exhausts on this bound, so solve.prove's euler_gamma
    fallback emits the cert; certificate.verify_cert's euler_gamma dispatch
    must accept it, and tampering must fail.
    """
    bound = str(rat_bound(mp.euler, 10, above=False))
    resp = solve.prove("gamma", "1", ">", bound, exact=True)
    assert resp.get("prover") == "euler_gamma" and resp["type"] == "gamma"
    cert = resp["certificate"]
    assert certificate.verify_cert(cert)
    bad = dict(cert, lo=str(Fraction(cert["lo"]) + 1))
    assert not certificate.verify_cert(bad)
