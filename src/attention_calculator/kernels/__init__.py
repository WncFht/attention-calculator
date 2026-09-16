"""Kernel family registry: type name -> prove() implementation.

Each family module defines::

    def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict

where ``kind`` is the requested site type (a module may serve several types).
It returns ``{"parameters": {...}, "solution": "a = 47/120, b = -13/120"}``
with ``parameters`` holding the nine site fields (m, n, a_val, b_val, c_val,
au_val, bu_val, cu_val, u_val). Raises WrongDirection or NoSolution from
..engine on failure.

Each family module also defines::

    def render_equation(params: dict, kind: str, power: Fraction | str,
                        comp: str, bound: Fraction | str) -> str

producing the same LaTeX proof string the site's /get_integral_image returns
(e.g. "\\dfrac{22}{7} - \\pi = \\int_0^1 ... \\mathrm{d} x > 0"). Raw
``power``/``bound`` request strings are echoed unreduced like the site does;
``render.wire_fraction`` parses them where the numeric value is needed.
"""

TYPES = [
    "pi",
    "e",
    "pi_n",
    "e_q",
    "ln_q",
    "ln_q_square",
    "sin_q",
    "cos_q",
    "tan_q",
    "cot_q",
    "sin_q_degree",
    "cos_q_degree",
    "sin_pi_q",
    "cos_pi_q",
    "arctan_q",
    "arccot_q",
    "sinh_q",
    "cosh_q",
    "tanh_q",
    "coth_q",
    "artanh_q",
    "arcoth_q",
    "gamma",
    "golden",
    "catalan",
    "zeta3",
    "e_pi",
    "varpi",
    "gauss",
]

# Exact-mode-only types (docs/2026-09-16-math-correctness-plan.md W3):
# no site counterpart — the real site 400s them, so TYPES is untouched and
# /calculate admits them solely under mode=exact.
EXACT_TYPES = [
    "zeta5",
    "zeta7",
    "ln_q_cube",
    "arcsin_q",
    "arsinh_q",
    "gaussint_q",
    "dawson_q",
    "erfiint_q",
]
