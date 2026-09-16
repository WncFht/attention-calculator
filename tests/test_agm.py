"""Tests for the AGM-iteration second prover (agm.py, plan doc W7).

The lemniscate kernels have the coarsest provable-bound floor of all
families (varpi '>' cannot even prove varpi > 2.6).  agm.py closes it
with rigorous rational AGM enclosures of M = AGM(1, sqrt(2)): a plain
cert for gauss (G = 1/M) and a composite DAG pairing an agm cert with a
pi kernel cert for varpi (varpi = pi*G).  Tests pin the certificate
verification, the floor improvement, tamper rejection, and the wiring:
solve.prove's exact-mode NoSolution fallback plus the certificate.verify_cert
agm/composite dispatch.
"""

import json
from fractions import Fraction

import pytest
from mpmath import mp

from attention_calculator import agm, certificate, solve
from attention_calculator.engine import WrongDirection
from attention_calculator.kernels import gamma_special


@pytest.fixture(autouse=True, scope="module")
def module_dps():
    """本模块数值校验的 mpmath 精度。"""
    with mp.workdps(60):
        yield


with mp.workdps(1200):
    GAUSS = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi**3))
    VARPI = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi))


def rat_bound(c: mp.mpf, digits: int, above: bool) -> Fraction:
    """Floor/ceiling rational approximant of c with denominator 10^digits.

    Everything runs at digits+200 dps so the result is independent of the
    ambient precision other test modules leave behind.
    """
    D = 10**digits
    with mp.workdps(digits + 200):
        k = int(c * D)
        b = Fraction(k, D)
        cb = mp.mpf(b.numerator) / b.denominator
        if above and cb <= c:
            b = Fraction(k + 1, D)
        if not above and cb >= c:
            b = Fraction(k - 1, D)
    return b


def prove(kind, q, comp, bound):
    """agm.prove + mandatory verification of every emitted certificate."""
    cert = agm.prove(kind, Fraction(q), comp, Fraction(bound))
    assert cert is not None, f"no certificate for {kind} {q} {comp} {bound}"
    assert certificate.verify_cert(cert)
    return cert


# ------------------------------------------------------------------- gauss


@pytest.mark.parametrize("digits", (3, 6, 12, 25))
def test_gauss_lower_below_kernel_floor(digits):
    # kernel floor: G > 0.8253 (gap ~9e-3); these sit far below it
    cert = prove("gauss", 1, ">", rat_bound(GAUSS, digits, False))
    assert cert["agm_digits"] == 64


@pytest.mark.parametrize("digits", (3, 6, 12, 25))
def test_gauss_upper_below_kernel_floor(digits):
    # kernel floor: G < 0.8351 (gap ~5e-4)
    cert = prove("gauss", 1, "<", rat_bound(GAUSS, digits, True))
    assert agm.verify_cert(cert)


def test_gauss_coefficient_and_sign_flip():
    # 2G < 7/4 normalizes to G < 7/8; -G < -4/5 normalizes to G > 4/5
    assert certificate.verify_cert(prove("gauss", 2, "<", Fraction(7, 4)))
    assert certificate.verify_cert(prove("gauss", -1, "<", Fraction(-4, 5)))


def test_gauss_negative_bound_trivial():
    # G > -1 is trivially true; the n=0 enclosure already decides it
    cert = prove("gauss", 1, ">", Fraction(-1))
    assert cert["agm_iter"] == 0


# ------------------------------------------------------------------- varpi


@pytest.mark.parametrize(
    "bound",
    ["2.61", "2.622", "2.6220575"],  # kernel '>' floor is ~2.5934
)
def test_varpi_lower_below_kernel_floor(bound):
    cert = prove("varpi", 1, ">", Fraction(bound))
    assert cert["rule"] == "pi_div_agm"
    assert [(e["kind"], e["comp"]) for e in cert["expect"]] == [
        ("pi", ">"),
        ("gauss", ">"),
    ]


@pytest.mark.parametrize(
    "bound",
    ["2.63", "2.6221", "2.62205756"],  # kernel '<' floor gap ~0.008
)
def test_varpi_upper_below_kernel_floor(bound):
    cert = prove("varpi", 1, "<", Fraction(bound))
    assert cert["rule"] == "pi_div_agm"


def test_varpi_coefficient_and_sign_flip():
    # 2ϖ > 5 normalizes to ϖ > 5/2; -ϖ < -25/10 normalizes to ϖ > 5/2
    assert certificate.verify_cert(prove("varpi", 2, ">", Fraction(5)))
    assert certificate.verify_cert(prove("varpi", -1, "<", Fraction(-25, 10)))


def test_varpi_negative_bound_lower():
    # ϖ > -3: trivially true, composite with A=1 still certifies it
    cert = prove("varpi", 1, ">", Fraction(-3))
    assert cert["witness"]["A"] == "1"


def test_varpi_child_claims_extractable():
    # the cert's own claims are readable by gamma_special.child_claim
    cert = agm.prove("varpi", Fraction(1), ">", Fraction(2_622_057, 1_000_000))
    assert gamma_special.child_claim(cert) == (
        "varpi",
        Fraction(1),
        ">",
        Fraction(2_622_057, 1_000_000),
    )
    for e, ch in zip(cert["expect"], cert["children"], strict=True):
        assert gamma_special.child_claim(ch) == (
            e["kind"],
            Fraction(e["power"]),
            e["comp"],
            Fraction(e["bound"]),
        )


# ------------------------------------------------------------ false claims


@pytest.mark.parametrize(
    ("kind", "comp", "bound"),
    [
        ("gauss", ">", "84/100"),  # G < 0.84
        ("gauss", "<", "83/100"),  # G > 0.83
        ("gauss", "<", "0"),  # G > 0
        ("varpi", ">", "2.63"),  # varpi < 2.63
        ("varpi", "<", "2.62"),  # varpi > 2.62
        ("varpi", "<", "0"),  # varpi > 0
    ],
)
def test_false_claim_returns_none(kind, comp, bound):
    assert agm.prove(kind, Fraction(1), comp, Fraction(bound)) is None


def test_deep_gap_needs_deeper_scale():
    # a 1e-100 gap: decided at the 512-digit rung (monkeypatch ladder to
    # prove the escalation works, and that 64 digits cannot decide)
    bound = rat_bound(GAUSS, 100, False)
    old = agm.LADDER
    agm.LADDER = (32, 512)
    try:
        cert = agm.prove("gauss", Fraction(1), ">", bound)
    finally:
        agm.LADDER = old
    assert cert is not None and cert["agm_digits"] == 512
    assert agm.verify_cert(cert)


def test_uncovered_kind_raises():
    with pytest.raises(ValueError):
        agm.prove("zeta3", Fraction(1), ">", Fraction(1))


def test_zero_coefficient_rejected():
    with pytest.raises(ValueError):
        agm.prove("gauss", Fraction(0), ">", Fraction(1))
    with pytest.raises(ValueError):
        agm.prove("varpi", Fraction(0), "<", Fraction(3))


def test_bad_comparison_rejected():
    with pytest.raises(ValueError):
        agm.prove("gauss", Fraction(1), "=", Fraction(1))


# -------------------------------------------------------- tamper rejection


def test_gauss_tampered_certs_fail():
    cert = agm.prove("gauss", Fraction(1), ">", Fraction(417, 500))
    assert agm.verify_cert(cert)
    bad = dict(cert, comparison="<")
    assert not agm.verify_cert(bad)  # wrong direction
    bad = dict(cert, bound=cert["bound"] + 1)
    assert not agm.verify_cert(bad)  # wrong bound
    bad = dict(cert, lo=cert["lo"] - Fraction(1, 10**70))
    assert not agm.verify_cert(bad)  # widened interval end
    bad = dict(cert, hi=cert["hi"] + Fraction(1, 10**70))
    assert not agm.verify_cert(bad)  # widened interval end
    bad = dict(cert, agm_iter=cert["agm_iter"] + 1)
    assert not agm.verify_cert(bad)  # wrong iteration count
    bad = dict(cert, agm_digits=cert["agm_digits"] + 1)
    assert not agm.verify_cert(bad)  # wrong scale
    bad = dict(cert, kind="varpi")
    assert not agm.verify_cert(bad)  # wrong kind
    bad = dict(cert, prover="pade")
    assert not agm.verify_cert(bad)  # wrong prover tag


def test_varpi_tampered_certs_fail():
    cert = agm.prove("varpi", Fraction(1), ">", Fraction(131, 50))
    # forged witness: shrink A below R/B to break the product side
    bad = dict(cert, witness={"A": "1", "B": cert["witness"]["B"]})
    assert not certificate.verify_cert(bad)
    # expect/child mismatch
    bad = dict(cert, expect=[dict(e, bound="1") for e in cert["expect"]])
    assert not certificate.verify_cert(bad)
    # flipped comparison
    bad = dict(cert, comp="<")
    assert not certificate.verify_cert(bad)
    # tampered child bound: widening the agm interval inside the child
    child_g = dict(cert["children"][1], hi=str(Fraction(cert["witness"]["B"]) ** -1 + 1))
    bad = dict(cert, children=[cert["children"][0], child_g])
    assert not certificate.verify_cert(bad)


def test_forged_certs_rejected():
    # plausible-looking but wrong enclosure (true values are n-dependent)
    fake = {
        "prover": "agm",
        "kind": "gauss",
        "power": Fraction(1),
        "comparison": ">",
        "bound": Fraction(0),
        "agm_iter": 3,
        "agm_digits": 64,
        "lo": Fraction(11, 10),  # not the replayed value
        "hi": Fraction(12, 10),
    }
    assert not agm.verify_cert(fake)
    # and a self-consistent forgery attempt on the claim still fails the
    # replay: swap in the n=4 enclosure under n=3's claim
    lo4, hi4 = agm.agm_enclosure(4, 10**64)
    fake["lo"], fake["hi"] = lo4, hi4
    assert not agm.verify_cert(fake)


# ------------------------------------------------------------- wire format


def test_json_round_trip():
    cert = agm.prove("gauss", Fraction(1), ">", Fraction(417, 500))
    wire = json.loads(json.dumps(agm.cert_jsonable(cert)))
    assert agm.verify_cert(agm.cert_parse(wire))
    cert = agm.prove("varpi", Fraction(1), "<", Fraction(263, 100))
    wire = json.loads(json.dumps(cert))  # composite certs are already wire
    assert certificate.verify_cert(wire)


# ------------------------------------------------------------- end to end


def test_solve_fallback_end_to_end():
    """End-to-end through the native wiring: below the kernel floor the
    prove_exact NoSolution fallback routes to agm / composite, and the
    emitted certs pass certificate.verify_cert."""
    # below the kernel floor -> kernel NoSolutions -> agm takes over
    bound = str(rat_bound(GAUSS, 6, False))
    resp = solve.prove("gauss", "1", ">", bound, exact=True)
    assert resp["prover"] == "agm"
    assert certificate.verify_cert(resp["certificate"])

    bound = str(Fraction(2_622_057, 1_000_000))
    resp = solve.prove("varpi", "1", ">", bound, exact=True)
    assert resp["prover"] == "composite"
    assert certificate.verify_cert(resp["certificate"])

    # false claim: certified_cmp raises before the fallback is reached
    with pytest.raises(WrongDirection):
        solve.prove("gauss", "1", ">", "9/10", exact=True)
    with pytest.raises(WrongDirection):
        solve.prove("varpi", "1", "<", "5/2", exact=True)

    # kernel-reachable claim keeps the classic shape (no prover key)
    resp = solve.prove("gauss", "1", ">", "4/5", exact=True)
    assert "prover" not in resp and "parameters" in resp


def test_gamma14_sqrt_mul_uses_agm_varpi():
    """Side benefit: a varpi child beyond the kernel floor still verifies —
    simulate by verifying the composite cert gamma_special would consume."""
    cert = agm.prove("varpi", Fraction(2), ">", Fraction(5_244_115, 1_000_000))
    # child_claim contract used by gamma_special.verify_cert
    assert gamma_special.child_claim(cert) == (
        "varpi",
        Fraction(2),
        ">",
        Fraction(5_244_115, 1_000_000),
    )
    assert certificate.verify_cert(cert)
