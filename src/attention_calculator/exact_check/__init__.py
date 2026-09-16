"""Exact proof verifier: recompute the claimed identity over QQ.

The math-correctness layer (docs/2026-09-16-math-correctness-plan.md W0).
Unlike bench/verify.py's 50dps numeric audit, this checks the emitted
parameters in exact rational arithmetic: every integrand moment is a
Moment dict {constant symbol: Fraction}, so ``∫ f dx = target`` is a
plain dict equality — a proof, not a sample.

Family modules (one file per kernels/ family, same module names as
solve.FAMILY) implement::

    def check(kind, power, comp, bound, params) -> dict

returning::

    {
        "identity_ok": bool,   # combine(coeffs, true_basis) == claim vector
        "nonneg": bool,        # P sign-definite on the domain (exact rules)
        "integrand": Moment,   # exact value of the encoded integrand
        "target": Moment,      # claimed s*(power*C - bound) vector
    }

The true-basis requirement matters: kernels that reproduce site bugs
(trig_pi's biased stored formula, beta_family's transposed solve) must
build the check basis from the mathematically correct moments, so site
mode emits params that fail here — that is the point.
"""

import importlib
from fractions import Fraction

from ..solve import FAMILY


def verify(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Dispatch to the family's checker on a response's parameters dict."""
    module = importlib.import_module(f"attention_calculator.exact_check.{FAMILY[kind]}")
    return module.check(kind, power, comp, bound, params)


def verify_response(kind: str, power: Fraction, comp: str, bound: Fraction, resp: dict) -> dict:
    """verify() on a /calculate response body.

    Padé-fallback and composite (proof-DAG) responses carry no
    parameters — the certificate IS the proof; identity_ok reports the
    certificate verifier's verdict (those certs bundle sign and identity
    in one check, so nonneg reports the same bit).
    """
    if resp.get("prover") == "pade":
        from .. import pade
        from ..kernels import gamma_special

        cert = resp["certificate"]
        try:
            match = gamma_special.child_claim(cert) == (kind, power, comp, bound)
        except (KeyError, TypeError, ValueError):
            match = False
        ok = match and pade.verify_cert(pade.cert_parse(cert))
        return {"identity_ok": ok, "nonneg": ok, "integrand": {}, "target": {}}
    if resp.get("prover") == "composite":
        from ..kernels import gamma_special

        cert = resp["certificate"]
        try:
            match = gamma_special.child_claim(cert) == (kind, power, comp, bound)
        except (KeyError, TypeError, ValueError):
            match = False
        ok = match and gamma_special.verify_cert(cert)
        return {"identity_ok": ok, "nonneg": ok, "integrand": {}, "target": {}}
    return verify(kind, power, comp, bound, resp["parameters"])
