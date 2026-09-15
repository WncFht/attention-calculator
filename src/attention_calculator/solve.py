"""Dispatch a (type, power, comparison, rational) request to its kernel family."""

import importlib
from fractions import Fraction

from .kernels import TYPES

# type -> module under attention_calculator.kernels
FAMILY = {
    "pi": "quadlog",
    "pi_n": "quadlog",
    "catalan": "quadlog",
    "zeta3": "quadlog",
    "arctan_q": "quadlog",
    "arccot_q": "quadlog",
    "e": "exp_family",
    "e_q": "exp_family",
    "e_pi": "exp_family",
    "sinh_q": "hyperbolic",
    "cosh_q": "hyperbolic",
    "tanh_q": "hyperbolic",
    "coth_q": "hyperbolic",
    "sin_q": "trig_q",
    "cos_q": "trig_q",
    "tan_q": "trig_q",
    "cot_q": "trig_q",
    "sin_pi_q": "trig_pi",
    "cos_pi_q": "trig_pi",
    "sin_q_degree": "trig_pi",
    "cos_q_degree": "trig_pi",
    "ln_q": "log_family",
    "ln_q_square": "log_family",
    "artanh_q": "log_family",
    "arcoth_q": "log_family",
    "golden": "beta_family",
    "varpi": "beta_family",
    "gauss": "beta_family",
    "gamma": "gamma",
}


def parse_rational(text: str) -> Fraction:
    """Parse '3', '22/7' into a Fraction."""
    return Fraction(text.strip())


def prove(kind: str, power: str, comp: str, rational: str) -> dict:
    """Run the proof search; returns the site's /calculate response shape.

    Raises engine.WrongDirection / engine.NoSolution on failure.
    """
    if kind not in TYPES:
        raise ValueError(f"unsupported type {kind!r}")
    module = importlib.import_module(f"attention_calculator.kernels.{FAMILY[kind]}")
    return module.prove(kind, parse_rational(power), comp, parse_rational(rational))
