"""Composite Γ-value types (kernels/gamma_special.py): gamma14/34/12.

No beta-moment kernel reaches Γ(1/4) — the parity lemma locks it into even
powers — so these prove via a small certificate DAG: rational sub-bounds on
S = 2ϖ (varpi kernel), T = √(2π) (pi kernel), V = π√2 (pi_sqrt2 kernel),
combined by an exact square-root/quotient transfer over QQ.  Tests pin the
transfer rules, the certificate verification, and the failure modes.
"""

from fractions import Fraction

import pytest

from attention_calculator import certificate, solve
from attention_calculator.engine import NoSolution, WrongDirection


def prove(kind, power, comp, bound):
    return solve.prove(kind, str(power), comp, str(bound), exact=True)


# ------------------------------------------------------------------- gamma12


def test_gamma12_lower_bound():
    # √π ≈ 1.7725; 7/4 = 1.75 -> child proves π > 49/16
    resp = prove("gamma12", 1, ">", Fraction(7, 4))
    cert = resp["certificate"]
    assert resp["prover"] == "composite" and cert["rule"] == "sqrt"
    assert certificate.verify_cert(cert)
    assert cert["expect"] == [{"kind": "pi", "power": "1", "comp": ">", "bound": "49/16"}]


def test_gamma12_upper_bound():
    resp = prove("gamma12", 1, "<", Fraction(2))
    assert certificate.verify_cert(resp["certificate"])


def test_gamma12_false_claim_rejected():
    with pytest.raises(WrongDirection):
        prove("gamma12", 1, ">", Fraction(2))


# ------------------------------------------------------------------- gamma14


def test_gamma14_upper_bound():
    # Γ(1/4) ≈ 3.6256; sqrt_mul '<': varpi proves 2ϖ < s1, pi proves π < t2²/2
    resp = prove("gamma14", 1, "<", Fraction(4))
    cert = resp["certificate"]
    assert cert["rule"] == "sqrt_mul" and certificate.verify_cert(cert)
    kinds = [(e["kind"], e["comp"]) for e in cert["expect"]]
    assert kinds == [("varpi", "<"), ("pi", "<")]


def test_gamma14_lower_bound():
    # '>' needs s1·t2 ≥ R² with s1 < S — inside the varpi '>' floor only for
    # modest gaps; R = 3 leaves comfortable room (R² = 9 vs g² ≈ 13.14)
    resp = prove("gamma14", 1, ">", Fraction(3))
    cert = resp["certificate"]
    assert cert["rule"] == "sqrt_mul" and certificate.verify_cert(cert)


def test_gamma14_positivity_transfer():
    # R ≤ 0 '>' is trivially true but still needs a certificate: a child
    # proving g > w > 0 discharges it
    resp = prove("gamma14", 1, ">", Fraction(-1))
    cert = resp["certificate"]
    assert cert["rule"] == "pos_transfer" and certificate.verify_cert(cert)


def test_gamma14_false_claim_rejected():
    with pytest.raises(WrongDirection):
        prove("gamma14", 1, ">", Fraction(4))  # g ≈ 3.6256 < 4


def test_gamma14_negative_bound_lt_rejected():
    with pytest.raises(WrongDirection):
        prove("gamma14", 1, "<", Fraction(-1))  # g > 0 > -1


def test_gamma14_coefficient_and_sign_flip():
    # 2·g < 8 normalizes to g < 4;  -g < -3 normalizes to g > 3
    resp = prove("gamma14", 2, "<", Fraction(8))
    assert certificate.verify_cert(resp["certificate"])
    resp = prove("gamma14", -1, "<", Fraction(-3))
    assert certificate.verify_cert(resp["certificate"])


def test_gamma14_zero_coefficient_rejected():
    with pytest.raises(ValueError):
        prove("gamma14", 0, ">", Fraction(-1))


def test_gamma14_tampered_cert_fails():
    resp = prove("gamma14", 1, ">", Fraction(3))
    cert = resp["certificate"]
    # forged product witness: s1 small enough to violate s1·t2 ≥ R²
    bad = dict(cert, witness={"s1": "1", "t2": cert["witness"]["t2"]})
    assert not certificate.verify_cert(bad)
    # child claim mismatch: expect table disagrees with the child cert
    bad = dict(cert, expect=[dict(e, bound="1") for e in cert["expect"]])
    assert not certificate.verify_cert(bad)
    # wrong rule
    bad = dict(cert, rule="sqrt")
    assert not certificate.verify_cert(bad)


def test_gamma14_honest_no_solution_on_tight_bound():
    # g = 3.62560990…; bound 36257/10000 sits ~9e-5 above it — a gap at the
    # edge of the child floors may honestly exhaust (NoSolution) or just make
    # it through; either outcome must be honest, never a bad cert
    try:
        resp = prove("gamma14", 1, "<", Fraction(36257, 10000))
    except NoSolution:
        return
    assert certificate.verify_cert(resp["certificate"])


# ------------------------------------------------------------------- gamma34
# Γ(3/4) = π√2/Γ(1/4) ≈ 1.2254; needs the pi_sqrt2 kernel (sqrt_div rule).


def test_gamma34_both_directions():
    resp = prove("gamma34", 1, ">", Fraction(1))
    cert = resp["certificate"]
    assert cert["rule"] in ("sqrt_div", "pos_transfer")
    assert certificate.verify_cert(cert)
    resp = prove("gamma34", 1, "<", Fraction(3, 2))
    cert = resp["certificate"]
    assert cert["rule"] == "sqrt_div" and certificate.verify_cert(cert)
    kinds = [(e["kind"], e["comp"]) for e in cert["expect"]]
    assert kinds == [("gamma14", ">"), ("pi_sqrt2", "<")]


def test_gamma34_false_claim_rejected():
    with pytest.raises(WrongDirection):
        prove("gamma34", 1, ">", Fraction(2))  # h ≈ 1.2254 < 2
