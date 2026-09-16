"""Tests for the W5 proof certificate: build / verify_cert / cert_tex / CLI.

Every case drives solve.prove(..., exact=True) — prove_exact already gates
emission on exact_check, so each built certificate must verify. Tamper cases
corrupt one field at a time: the recorded check block must match a fresh
recomputation verbatim, so any edit — claim fields, parameters, recorded
moments, recorded verdicts — fails verification.
"""

import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from attention_calculator import solve
from attention_calculator.certificate import build, cert_tex, verify_cert

ROOT = Path(__file__).resolve().parent.parent

# one or more emitted proofs per registered family (cases proven in the
# test_exact_* suites)
CASES = [
    ("pi", "1", ">", "3"),
    ("pi", "1", "<", "22/7"),
    ("pi_n", "2", "<", "10"),
    ("catalan", "1", ">", "9/10"),
    ("zeta3", "1", "<", "5/4"),
    ("arctan_q", "1", ">", "7/9"),
    ("arccot_q", "2", "<", "1/2"),
    ("e", "1", "<", "3"),
    ("e_q", "1/2", "<", "7/4"),
    ("e_pi", "2", "<", "47"),
    ("ln_q", "2", ">", "11/16"),
    ("ln_q_square", "5", ">", "2"),
    ("artanh_q", "1/2", "<", "11/20"),
    ("sin_q", "1/2", ">", "11/23"),
    ("tan_q", "1/2", "<", "3/5"),
    ("sin_pi_q", "1/5", ">", "1/2"),
    ("cos_q_degree", "18", ">", "9/10"),
    ("sinh_q", "1/2", "<", "1"),
    ("tanh_q", "1/2", ">", "2/5"),
    ("coth_q", "1", ">", "13/10"),
    ("golden", "1", ">", "3/2"),
    ("varpi", "1", ">", "5/2"),
    ("gauss", "1", "<", "9/10"),
    ("gamma", "1", ">", "57/100"),
    ("zeta5", "1", "<", "28/27"),
    ("zeta7", "1", ">", "121/120"),
]


def certify(kind: str, power: str, comp: str, bound: str) -> dict:
    """Emit an exact-mode proof and build its certificate."""
    resp = solve.prove(kind, power, comp, bound, exact=True)
    return build(kind, Fraction(power), comp, Fraction(bound), resp["parameters"])


@pytest.mark.parametrize(("kind", "power", "comp", "bound"), CASES)
def test_emitted_cert_verifies(kind, power, comp, bound):
    cert = certify(kind, power, comp, bound)
    assert cert["check"]["identity_ok"]
    assert cert["check"]["nonneg"]
    assert verify_cert(cert)
    # a serialized-then-reparsed cert is the same proof
    assert verify_cert(json.loads(json.dumps(cert, sort_keys=True)))


def test_cert_tex_summary_line():
    tex = cert_tex(certify("pi", "1", "<", "22/7"))
    assert "\\dfrac{22}{7}" in tex and "\\pi" in tex and tex.endswith("> 0")


def _tampers() -> list:
    """Single-field corruptions of a pi>3 certificate."""
    return [
        lambda c: c["parameters"].update(au_val=str(int(c["parameters"]["au_val"]) + 1)),
        lambda c: c["parameters"].update(m=c["parameters"]["m"] + 1),
        lambda c: c["check"]["target"].update({"1": "-4"}),
        lambda c: c["check"]["integrand"].update(pi="2"),
        lambda c: c.update(bound="4"),
        lambda c: c.update(power="2"),
        lambda c: c.update(comparison="<"),
        lambda c: c["check"].update(identity_ok=False),
        lambda c: c["check"].update(nonneg="yes"),
        lambda c: c.update(kind="bogus"),
        lambda c: c.pop("check"),
        lambda c: c.update(power="abc"),
        lambda c: c.update(bound=0.5),
    ]


@pytest.mark.parametrize("tamper", _tampers())
def test_tampered_cert_fails(tamper):
    cert = certify("pi", "1", ">", "3")
    tamper(cert)
    assert not verify_cert(cert)


def test_cli(tmp_path):
    tool = ROOT / "tools" / "verify_cert.py"
    good = tmp_path / "good.json"
    good.write_text(json.dumps(certify("pi", "1", ">", "3")))
    run = subprocess.run(
        [sys.executable, str(tool), str(good)], capture_output=True, text=True, check=False
    )
    assert run.returncode == 0, run.stderr
    assert "VERIFIED" in run.stdout

    bad = tmp_path / "bad.json"
    cert = json.loads(good.read_text())
    cert["bound"] = "4"
    bad.write_text(json.dumps(cert))
    run = subprocess.run(
        [sys.executable, str(tool), str(bad)], capture_output=True, text=True, check=False
    )
    assert run.returncode == 1
    assert "FAILED" in run.stderr

    run = subprocess.run(
        [sys.executable, str(tool)], input="not json", capture_output=True, text=True, check=False
    )
    assert run.returncode == 2
