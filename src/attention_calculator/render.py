"""Dispatch /get_integral_image rendering to the kernel family module."""

import importlib
from fractions import Fraction

from .solve import FAMILY

INT_KEYS = ("m", "n", "au_val", "bu_val", "cu_val", "u_val")
FRAC_KEYS = ("a_val", "b_val", "c_val")


def coerce_params(query: dict) -> dict:
    """Turn the nine site query fields into the typed params kernels expect.

    m, n, au_val, bu_val, cu_val, u_val become int; a_val, b_val, c_val become
    Fraction (the same shapes the family ``prove`` emits and golden.jsonl
    records). Raises ValueError on malformed input.
    """
    params = {k: int(query[k]) for k in INT_KEYS}
    params.update({k: Fraction(query[k]) for k in FRAC_KEYS})
    return params


def render_equation(
    params: dict, kind: str, power: Fraction, comp: str, bound: Fraction
) -> str:
    """Render the proof-equation LaTeX for a solved case.

    ``params`` is the typed dict from ``coerce_params``. Returns the same
    string the site's /get_integral_image puts under ``equation`` (no
    delimiters).
    """
    if kind not in FAMILY:
        raise ValueError(f"unsupported type {kind!r}")
    module = importlib.import_module(f"attention_calculator.kernels.{FAMILY[kind]}")
    return module.render_equation(params, kind, power, comp, bound)
