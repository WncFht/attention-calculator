"""Semantic forgery cases for verify_cert — audit-exact-math.md M1/M2.

M1 reproduced three end-to-end forgeries that verify_cert accepted for
FALSE claims: a negative u_val flipping the gamma sub-term's true sign,
a negative m accepted as formal moments of a divergent integrand, and a
tampered c_val rebuilding the artanh basis on the wrong reduced argument.
M2 traces all three to the missing params-domain clause — every case
below was internally consistent under the old rules (identity_ok and
nonneg recomputed verbatim True), so rejection must come from the new
entry gate (m,n >= 0, denominators > 0, power in the kind's domain,
kind-derived fields recomputed), not from a moment mismatch.

Each cert is constructed in-code: parameters plus the check block the
forger would have recorded — the values the pre-fix verify recomputed.
"""

from attention_calculator.certificate import verify_cert


def test_gamma_negative_u_val() -> None:
    """M1.1: claim gamma > 1 (false). u_val=-1 makes the true sub-polynomial
    -1-x/2 <= 0, but the old nonneg read the undivided au,bu = 1,1/2 and
    passed."""
    cert = {
        "kind": "gamma",
        "power": "1",
        "comparison": ">",
        "bound": "1",
        "parameters": {
            "m": 0,
            "n": 1,
            "cu_val": "1",
            "u_val": "-1",
            "a_val": "0",
            "au_val": "1",
            "bu_val": "1/2",
        },
        "check": {
            "integrand": {"1": "-1", "gamma": "1"},
            "target": {"1": "-1", "gamma": "1"},
            "identity_ok": True,
            "nonneg": True,
        },
    }
    assert verify_cert(cert) is False


def test_pi_negative_m() -> None:
    """M1.2: claim pi > 4 (false). m=-1 gives formal moments of the
    divergent integrand x^-1(1-x^2)^2*K — the QQ identity held verbatim
    while certifying no nonneg integrand at all."""
    cert = {
        "kind": "pi",
        "power": "1",
        "comparison": ">",
        "bound": "4",
        "parameters": {"m": -1, "n": 2, "au_val": "4", "bu_val": "9", "u_val": "5"},
        "check": {
            "integrand": {"1": "-4", "pi": "1"},
            "target": {"1": "-4", "pi": "1"},
            "identity_ok": True,
            "nonneg": True,
        },
    }
    assert verify_cert(cert) is False


def test_artanh_c_val_tampered() -> None:
    """M1.3: claim artanh(1/2) > 3/5 (false). c_val is attacker input that
    rebuilt the basis on c = c_val-1 = 3; the true reduced argument is
    qtilde(artanh_q, 1/2) = 3, so c_val="4" spanned the wrong basis while
    the recorded moments recomputed verbatim."""
    cert = {
        "kind": "artanh_q",
        "power": "1/2",
        "comparison": ">",
        "bound": "3/5",
        "parameters": {
            "m": 0,
            "n": 2,
            "au_val": "27",
            "bu_val": "162",
            "u_val": "80",
            "c_val": "4",
        },
        "check": {
            "integrand": {"1": "-3/5", "ln": "1/2"},
            "target": {"1": "-3/5", "ln": "1/2"},
            "identity_ok": True,
            "nonneg": True,
        },
    }
    assert verify_cert(cert) is False


def test_zero_denominator() -> None:
    """M2 variant: u_val=0 on a non-gamma kind. The denominator divides
    every numerator field, so zero must be rejected at the entry gate —
    not surface as a ZeroDivisionError deep in the checker."""
    cert = {
        "kind": "pi",
        "power": "1",
        "comparison": ">",
        "bound": "3",
        "parameters": {"m": 1, "n": 1, "au_val": "1", "bu_val": "5", "u_val": "0"},
        "check": {
            "integrand": {"1": "-3", "pi": "1"},
            "target": {"1": "-3", "pi": "1"},
            "identity_ok": True,
            "nonneg": True,
        },
    }
    assert verify_cert(cert) is False


def test_negative_m_variant() -> None:
    """M2 variant: claim pi < 1 (false). m=-2, n=2 solves to the nonneg
    coefficients (0, 1) on the formal basis — old verify accepted; the
    m >= 0 domain clause now rejects before the checker runs."""
    cert = {
        "kind": "pi",
        "power": "1",
        "comparison": "<",
        "bound": "1",
        "parameters": {"m": -2, "n": 2, "au_val": "0", "bu_val": "1", "u_val": "1"},
        "check": {
            "integrand": {"1": "1", "pi": "-1"},
            "target": {"1": "1", "pi": "-1"},
            "identity_ok": True,
            "nonneg": True,
        },
    }
    assert verify_cert(cert) is False


def test_power_out_of_domain() -> None:
    """M2 variant: ln_q with power=-1 is outside the family's q>1 domain —
    ln(-1) is not a real constant, yet the formal {ln,1} basis still solved
    to the nonneg coefficients (0, 8) and old verify accepted."""
    cert = {
        "kind": "ln_q",
        "power": "-1",
        "comparison": "<",
        "bound": "0",
        "parameters": {"m": 0, "n": 1, "au_val": "0", "bu_val": "8", "u_val": "1"},
        "check": {
            "integrand": {"ln": "-1"},
            "target": {"ln": "-1"},
            "identity_ok": True,
            "nonneg": True,
        },
    }
    assert verify_cert(cert) is False
