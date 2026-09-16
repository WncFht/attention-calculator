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

# Nonnegative-integer fields beyond the universal m/n, keyed by family module:
# gamma's cu_val is the x^N main-kernel exponent; beta_family/dixon's cu_val a
# fallback-template flag; sicin's t_val the Taylor-remainder level; li2's d the
# shift-family degree.  (In every other family *_val fields are free rational
# numerator coefficients and must NOT be sign-checked.)
INT_KEYS = {
    "gamma": ("cu_val",),
    "beta_family": ("cu_val",),
    "dixon": ("cu_val",),
    "sicin": ("t_val",),
    "li2": ("d",),
}


def verify(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Dispatch to the family's checker on a response's parameters dict.

    params_domain runs first: the checkers are moment recomputers that trust
    params to encode a real, convergent, sign-checkable integrand — untrusted
    certificates need the domain gate before any family code runs.
    """
    params_domain(kind, power, params)
    module = importlib.import_module(f"attention_calculator.exact_check.{FAMILY[kind]}")
    return module.check(kind, power, comp, bound, params)


def nonneg_int(v: object) -> bool:
    """True when v encodes a nonnegative integer (int or its decimal text)."""
    if isinstance(v, bool):
        return False
    if isinstance(v, int):
        return v >= 0
    if isinstance(v, str):
        try:
            return int(v) >= 0
        except ValueError:
            return False
    return False


def params_domain(kind: str, power: Fraction, params: dict) -> None:
    """Reject params that cannot encode a legal proof for ``kind``.

    Trust boundary for untrusted certificates (the emitted parameters are
    attacker-controlled JSON, not kernel output):

    - every exponent field (``m``/``n`` plus the INT_KEYS family keys) must be
      a nonnegative integer — a negative exponent makes the encoded factor
      diverge at an endpoint while its formal moment still computes, and a
      negative index silently wraps Python moment tables;
    - the printed denominator ``u_val`` must be positive — P's true
      coefficients are ``au_val/u`` etc., so ``u < 0`` flips every sign a
      nonneg rule certifies (gamma alone admits ``u_val == 0``, its marker
      for the degenerate N=0 print);
    - ``power`` must sit inside the kind's certified-constant domain
      (``power_in_domain``).

    Raises ValueError on violations; certificate.verify_cert maps that to
    "not a proof" without reaching the family checker.
    """
    fam = FAMILY[kind]
    for key in ("m", "n", *INT_KEYS.get(fam, ())):
        if key in params and not nonneg_int(params[key]):
            raise ValueError(f"parameters[{key!r}] must be a nonnegative int, got {params[key]!r}")
    if "u_val" in params:
        u = Fraction(params["u_val"])
        if u < 0 or (u == 0 and fam != "gamma"):
            raise ValueError(f"parameters['u_val'] must be positive, got {params['u_val']!r}")
    if not power_in_domain(kind, power):
        raise ValueError(f"power {power} is outside {kind}'s certified domain")


def power_in_domain(kind: str, q: Fraction) -> bool:
    """Whether ``q`` lies inside ``kind``'s certified-constant domain.

    This is the math boundary a certificate needs, not the site's input-format
    rules: outside it the encoded integrand diverges, the constant is non-real,
    or the family's kernel factor loses sign-definiteness — a formal moment
    equality there certifies nothing.  The bounds mirror each kernel's own
    domain check; kinds where ``power`` is a pure coefficient (the moments
    never see it — quadlog, beta_family, exp_family's e/e_pi, dixon, …)
    accept ``q == 0``, matching the site's own degenerate emissions.
    """
    fam = FAMILY[kind]
    if fam == "log_family":
        return {
            "ln_q": q > 1,
            "ln_q_square": q > 1,
            "artanh_q": 0 < q < 1,
            "arcoth_q": q > 1,
        }[kind]
    if fam == "ln_pow":
        return q > 1
    if fam == "trig_q":
        if q <= 0:
            return False
        import sympy as sp

        upper = sp.pi if kind == "sin_q" else sp.pi / 2
        return bool(sp.Rational(q) < upper)
    if fam == "trig_pi":
        return 0 < (q / 180 if kind.endswith("_degree") else q) < Fraction(1, 2)
    if fam == "arcsin":
        return 0 < q < 1
    if fam == "li2":
        return q != 0 and q < 1
    if fam in ("gamma", "trigamma"):
        return q > 0
    if fam == "hyperbolic":
        return q > 0 if kind == "coth_q" else q != 0
    if fam == "quadlog":
        if kind == "pi_n":
            return q.numerator >= 1 and q.denominator <= 64
        # arctan/arccot bases carry a 1/q division; pi/catalan/zeta3 take
        # power as a pure coefficient (q=0 stays a site-legal emission)
        return kind not in ("arctan_q", "arccot_q") or q != 0
    if fam == "exp_family":
        return kind != "e_q" or q != 0
    if fam in ("sicin", "gauss_erf", "dixon"):
        return q != 0
    return True


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
    if resp.get("prover") == "agm":
        from .. import agm
        from ..kernels import gamma_special

        cert = resp["certificate"]
        try:
            match = gamma_special.child_claim(cert) == (kind, power, comp, bound)
        except (KeyError, TypeError, ValueError):
            match = False
        ok = match and agm.verify_cert(agm.cert_parse(cert))
        return {"identity_ok": ok, "nonneg": ok, "integrand": {}, "target": {}}
    if resp.get("prover") == "euler_gamma":
        from .. import euler_gamma
        from ..kernels import gamma_special

        cert = resp["certificate"]
        try:
            match = gamma_special.child_claim(cert) == (kind, power, comp, bound)
        except (KeyError, TypeError, ValueError):
            match = False
        ok = match and euler_gamma.verify_cert(euler_gamma.cert_parse(cert))
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
