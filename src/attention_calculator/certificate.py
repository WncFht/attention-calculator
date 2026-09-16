"""Machine-checkable certificates for exact-mode proofs (plan doc W5).

A certificate bundles everything an offline verifier needs to recheck an
emitted proof without trusting the pipeline that produced it::

    {
        "kind": "pi",
        "power": "1",
        "comparison": ">",
        "bound": "3",
        "parameters": {"m": 0, "n": 4, "a_val": "0", ..., "u_val": "120"},
        "check": {
            "integrand": {"1": "-3", "pi": "1"},
            "target": {"1": "-3", "pi": "1"},
            "identity_ok": true,
            "nonneg": true
        }
    }

Rationals serialize as ``str(Fraction)`` — "n" when the denominator is 1,
"n/d" otherwise. Moments (moment.Moment: {constant symbol: Fraction})
serialize with symbols sorted, so ``json.dumps(cert, sort_keys=True)`` is
canonical.

``verify_cert`` is the independent recheck: it schema-validates the cert,
recomputes the family's true basis through exact_check.verify (the family
checker IS the math — a QQ dict equality between the integrand's moments
and the claimed target vector, plus the exact sign rule on P), and
requires the recorded check block to match recomputation verbatim, so a
tampered coefficient, field, or verdict fails.

Format spec: docs/2026-09-16-certificate-spec.md.
"""

from fractions import Fraction

from . import exact_check
from .moment import Moment
from .render import rat_tex
from .solve import FAMILY

__all__ = ["build", "cert_tex", "recheck", "verify_cert"]


def build(kind: str, power: Fraction, comp: str, bound: Fraction, params: dict) -> dict:
    """Assemble the certificate for an emitted proof; runs exact_check.verify."""
    res = exact_check.verify(kind, power, comp, bound, params)
    return {
        "kind": kind,
        "power": str(power),
        "comparison": comp,
        "bound": str(bound),
        "parameters": dict(params),
        "check": {
            "integrand": _moment_out(res["integrand"]),
            "target": _moment_out(res["target"]),
            "identity_ok": res["identity_ok"],
            "nonneg": res["nonneg"],
        },
    }


def recheck(cert: dict) -> dict:
    """Recompute a certificate's check block from the kernels.

    Returns the fresh exact_check.verify result (Moments as Fraction
    dicts). Raises ValueError on schema violations; checker exceptions
    from semantically broken parameters propagate.
    """
    kind, power, comp, bound, params, _recorded = _parse(cert)
    return exact_check.verify(kind, power, comp, bound, params)


def verify_cert(cert: dict) -> bool:
    """Independent re-verification: the cert certifies a true proof iff its
    schema is sane, the recorded check block matches a fresh recomputation
    verbatim, and that recomputation proves the claim (identity_ok and
    nonneg both hold).

    Padé certificates (a "serr" key marks the W4 second-prover schema) are
    dispatched to pade.verify_cert — a different exact argument over the
    same wire field name.
    """
    if isinstance(cert, dict) and "serr" in cert:
        from . import pade

        try:
            return pade.verify_cert(pade.cert_parse(cert))
        except (KeyError, TypeError, ValueError, ZeroDivisionError, AttributeError):
            return False
    if isinstance(cert, dict) and cert.get("prover") == "composite":
        from .kernels import gamma_special

        try:
            return gamma_special.verify_cert(cert)
        except (KeyError, TypeError, ValueError, ZeroDivisionError, AttributeError):
            return False
    try:
        kind, power, comp, bound, params, recorded = _parse(cert)
        res = exact_check.verify(kind, power, comp, bound, params)
    except (KeyError, TypeError, ValueError, ZeroDivisionError, AttributeError, ImportError):
        return False  # schema-shaped but unverifiable: not a proof
    return (
        res["identity_ok"]
        and res["nonneg"]
        and recorded["identity_ok"] == res["identity_ok"]
        and recorded["nonneg"] == res["nonneg"]
        and recorded["integrand"] == res["integrand"]
        and recorded["target"] == res["target"]
    )


def cert_tex(cert: dict) -> str:
    r"""Render the certified statement, e.g. '22/7 - \pi = \int f > 0'.

    The target vector already carries the comparison sign (it is
    s·(power·C − bound)), so the integrand's certified inequality is
    always '> 0' — the same shape the site's proof image prints.
    Composite and Padé certs have no check block; print the claim itself.
    """
    if cert.get("prover") == "composite":
        return (
            f"{cert['kind']}({cert['q']}) {cert['comp']} {cert['p']} \\quad (\\mathrm{{composite}})"
        )
    if "serr" in cert:
        return f"{cert['kind']} {cert['comp']} {cert['p']} \\quad (\\mathrm{{Pad\\acute{{e}}}})"
    return f"{_moment_tex(_moment_in(cert['check']['target']))} = \\int f\\,\\mathrm{{d}}x > 0"


def _parse(cert: dict) -> tuple[str, Fraction, str, Fraction, dict, dict]:
    """Schema-validate and decode a certificate; ValueError on any violation."""
    if not isinstance(cert, dict):
        raise ValueError("certificate must be a JSON object")
    kind = cert.get("kind")
    if kind not in FAMILY:
        raise ValueError(f"unknown kind {kind!r}")
    comp = cert.get("comparison")
    if comp not in (">", "<"):
        raise ValueError(f"bad comparison {comp!r}")
    power, bound = _frac(cert.get("power"), "power"), _frac(cert.get("bound"), "bound")
    params = cert.get("parameters")
    if not isinstance(params, dict):
        raise ValueError("parameters must be an object")
    check = cert.get("check")
    if not isinstance(check, dict):
        raise ValueError("check must be an object")
    recorded = {
        "integrand": _moment_in(check.get("integrand")),
        "target": _moment_in(check.get("target")),
        "identity_ok": _flag(check.get("identity_ok"), "identity_ok"),
        "nonneg": _flag(check.get("nonneg"), "nonneg"),
    }
    return kind, power, comp, bound, params, recorded


def _frac(v: object, field: str) -> Fraction:
    """Decode a serialized rational; floats and bools are not exact rationals."""
    if isinstance(v, bool) or not isinstance(v, int | str):
        raise ValueError(f"{field} must be a rational string, got {v!r}")
    try:
        return Fraction(v)
    except ValueError:
        raise ValueError(f"{field} is not a rational: {v!r}") from None


def _flag(v: object, field: str) -> bool:
    """Decode a recorded verdict bit; must be a JSON boolean."""
    if not isinstance(v, bool):
        raise ValueError(f"check.{field} must be a boolean, got {v!r}")
    return v


def _moment_in(d: object) -> Moment:
    """Decode a serialized moment {symbol: rational}; ValueError on bad shape."""
    if not isinstance(d, dict):
        raise ValueError("moment must be an object")
    return {sym: _frac(v, f"moment[{sym!r}]") for sym, v in d.items()}


def _moment_out(m: Moment) -> dict:
    """Serialize a moment deterministically: symbols sorted, coeffs 'n[/d]'."""
    return {sym: str(m[sym]) for sym in sorted(m)}


# constant symbol -> LaTeX for cert_tex; unknown symbols fall back to \mathrm
_SYMBOL_TEX = {
    "pi": "\\pi",
    "e": "e",
    "e_q": "e^{q}",
    "e_pi": "e^{\\pi}",
    "ln": "\\ln",
    "ln2": "\\ln^{2}",
    "ln3": "\\ln^{3}",
    "catalan": "C",
    "zeta3": "\\zeta(3)",
    "zeta5": "\\zeta(5)",
    "zeta7": "\\zeta(7)",
    "zeta9": "\\zeta(9)",
    "zeta11": "\\zeta(11)",
    "beta4": "\\beta(4)",
    "beta6": "\\beta(6)",
    "gamma": "\\gamma",
    "arctan": "\\arctan q",
    "sin_q": "\\sin q",
    "cos_q": "\\cos q",
    "sinh_q": "\\sinh q",
    "cosh_q": "\\cosh q",
    "sin_pi_q": "\\sin(\\pi q)",
    "cos_pi_q": "\\cos(\\pi q)",
    "phi": "\\varphi",
    "varpi": "\\varpi",
    "varpi_inv": "\\varpi^{-1}",
    "gauss": "G",
    "gauss_inv": "G^{-1}",
    "pi_sqrt2": "\\pi\\sqrt{2}",
    "pi3": "\\pi_{3}",
    "pi3_inv": "\\pi_{3}^{-1}",
    "U": "U",
    "G3": "G_{3}",
    "A": "A",
    "li2_q": "\\mathrm{Li}_{2}(q)",
    "ln_1mq": "\\ln(1-q)",
    "psi1_q": "\\psi'(q)",
    "si_q": "\\mathrm{Si}(q)",
    "cin_q": "\\mathrm{Cin}(q)",
    "C": "C",
}


def _moment_tex(m: Moment) -> str:
    """Render a moment as a QQ-linear combination; positive terms lead and the
    '1' term sorts last, matching the site's 'r - C' / 'C - r' print order."""
    out = ""
    for sym in sorted(m, key=lambda s: (m[s] <= 0, s == "1", s)):
        c = m[sym]
        name = _SYMBOL_TEX.get(sym, f"\\mathrm{{{sym}}}")
        term = rat_tex(abs(c)) if sym == "1" else f"{'' if abs(c) == 1 else rat_tex(abs(c))}{name}"
        out += ("-" if c < 0 else "") + term if not out else f" {'-' if c < 0 else '+'} {term}"
    return out or "0"
