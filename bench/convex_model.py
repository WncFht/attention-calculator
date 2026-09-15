"""Offline model of the /convex/prove pipeline.

Mirrors src/attention_calculator/convex.py (the byte-exact port of upstream
ConvexConcaveProver scripts/prove.py) minus the latex/text rendering, so the
sweep can compare status/minimum/proof/reason/errors against captured site
responses. Verified against bench/data/convex-probes.jsonl (333 records):
status 225/225, minimum.x/.value 200/200, tangent_at 77/77,
left_gap_min 77/77, reason 225/225, provided_line 14/14, errors 98/98 exact.

Key semantics (upstream prove.py order):

- normalize_expr: strip -> '^'->'**' -> regex \\bln\\s*\\( -> 'log(' ->
  regex \\be\\s*\\*\\*\\s*x\\b -> 'exp(x)'.
- split_inequality searches ops in order >=, <=, >, <, = and splits on the
  FIRST occurrence of the first op found (not first comparison char).
  '<'/'<=' swap sides. No op -> 'expr >= 0'. '=' is non-strict.
- parse_number: Constant, USub, Div, Pow, sqrt(<number>) only — no
  Add/Sub/Mult/UAdd. collect_terms const-folds Mult/Div (distributive over
  sums); parse_product folds const*atom; atom*atom -> 'products of two
  non-constant atoms are unsupported'. Call: argcount check for ANY name
  ('{name} expects one argument'), then arg must be bare x
  ('{name} only supports argument x'), then dispatch {exp,log,sqrt}.
- combine_terms merges by atom (dict order), drops |coef| <= 1e-9.
- split_positive_negative: c >= 0 stays left; c < 0 negates to right.
- evaluate/derivative: naive sequential accumulation (site ran old CPython;
  Python 3.12+ sum() is Neumaier-compensated and differs in the last ulp).
- curvature: sign(coeff * d2-atom-sign) with |val| <= 1e-9 -> 0; side is
  affine/convex/concave/mixed. Template: left in {convex,affine} and right in
  {concave,affine} and diff in {convex,affine} else inconclusive.
- minimize_convex: lo_eff = max(lo, 1e-8). With hi: d(lo)>=0 -> (lo,f(lo));
  d(hi)<=0 -> (hi,f(hi)); else 160 halvings d(m)<=0 -> a=m else b=m, return
  midpoint. Without hi: d(lo)>=0 -> boundary; else b doubles from 1.0 while
  d(b)<=0 and b<1e8 (reaches 2^27=134217728), still <=0 -> (b,f(b)); else
  160 halvings.
- status: strict ('>'/'<') needs f_min > 1e-8; non-strict f_min >= -1e-8.
- find_line: CF convergents of x_min via float recurrence (max 12 partial
  quotients, stop |frac|<1e-14, fold; keep denom<=10000, dedup, order).
  Skip x0 <= domain[0] (RAW lo, unclamped) or x0 >= domain[1] (only when hi
  given; default None = no upper check). Line = tangent of RIGHT at x0;
  accept first with min(left-line) > 1e-8 (strict) / >= -1e-7 (non-strict).
- --line parsed only after status=proved: affine atoms else
  '--line must be an affine expression in x'; gaps via minimize_convex;
  ok = both gaps >= -1e-8.
- parse_domain runs before inequality parse: split(',', strip parts, exactly
  2 parts; lo=float; hi None when lower() in {inf,infinity,+inf} else float.
  Route layer first: blank -> '请输入一个不等式。', len(strip)>500 ->
  '输入过长，请输入一个较短的单变量不等式'.
"""

import ast
import math
import re
from fractions import Fraction

EPS = 1e-9


class SiteErr(Exception):
    """Error the site reports as a 400 verbatim (route layer or str(exc))."""


def normalize_expr(text):
    text = text.strip()
    text = text.replace("^", "**")
    text = re.sub(r"\bln\s*\(", "log(", text)
    text = re.sub(r"\be\s*\*\*\s*x\b", "exp(x)", text)
    return text


def split_inequality(text):
    text = normalize_expr(text)
    for op in (">=", "<=", ">", "<", "="):
        if op in text:
            lhs, rhs = text.split(op, 1)
            if op in ("<=", "<"):
                lhs, rhs = rhs, lhs
            return lhs.strip(), rhs.strip(), op
    return text, "0", ">="


def dump(node):
    return ast.dump(node, show_empty=True)


def parse_number(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -parse_number(node.operand)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        return parse_number(node.left) / parse_number(node.right)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
        return parse_number(node.left) ** parse_number(node.right)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "sqrt":
        if len(node.args) != 1:
            raise SiteErr("sqrt expects one argument")
        return math.sqrt(parse_number(node.args[0]))
    raise SiteErr(f"expected numeric constant, got {dump(node)}")


def parse_atom(node):
    if isinstance(node, ast.Name) and node.id == "x":
        return ("linear",)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
        name = node.func.id
        if len(node.args) != 1:
            raise SiteErr(f"{name} expects one argument")
        arg = node.args[0]
        if not (isinstance(arg, ast.Name) and arg.id == "x"):
            raise SiteErr(f"{name} only supports argument x")
        if name == "exp":
            return ("exp",)
        if name == "log":
            return ("log",)
        if name == "sqrt":
            return ("power", 0.5)
    if (
        isinstance(node, ast.BinOp)
        and isinstance(node.op, ast.Pow)
        and isinstance(node.left, ast.Name)
        and node.left.id == "x"
    ):
        return ("power", parse_number(node.right))
    raise SiteErr(f"unsupported atom {dump(node)}")


def parse_factor_term(node):
    try:
        return (parse_number(node), ("const",))
    except (SiteErr, ValueError):
        return (1.0, parse_atom(node))


def parse_product(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        left = parse_product(node.left)
        right = parse_product(node.right)
        if left[1][0] == "const":
            return (right[0] * left[0], right[1])
        if right[1][0] == "const":
            return (left[0] * right[0], left[1])
        raise SiteErr("products of two non-constant atoms are unsupported")
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left = parse_product(node.left)
        denom = parse_number(node.right)
        return (left[0] * (1.0 / denom), left[1])
    return parse_factor_term(node)


def collect_terms(node, sign=1.0):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return collect_terms(node.left, sign) + collect_terms(node.right, sign)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Sub):
        return collect_terms(node.left, sign) + collect_terms(node.right, -sign)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        try:
            factor = parse_number(node.left)
            return collect_terms(node.right, sign * factor)
        except (SiteErr, ValueError):
            pass
        try:
            factor = parse_number(node.right)
            return collect_terms(node.left, sign * factor)
        except (SiteErr, ValueError):
            pass
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        try:
            factor = parse_number(node.right)
            return collect_terms(node.left, sign / factor)
        except (SiteErr, ValueError):
            pass
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return collect_terms(node.operand, -sign)
    return [(parse_product(node)[0] * sign, parse_product(node)[1])]


def parse_terms(text):
    tree = ast.parse(normalize_expr(text), mode="eval")
    return collect_terms(tree.body)


def combine_terms(terms):
    buckets = {}
    for c, a in terms:
        buckets[a] = buckets.get(a, 0.0) + c
    return [(c, a) for a, c in buckets.items() if abs(c) > EPS]


def split_positive_negative(terms):
    left, right = [], []
    for c, a in terms:
        if c >= 0:
            left.append((c, a))
        else:
            right.append((-c, a))
    return left, right


def atom_value(atom, x):
    kind = atom[0]
    if kind == "const":
        return 1.0
    if kind == "linear":
        return x
    if kind == "exp":
        return math.exp(x)
    if kind == "log":
        return math.log(x)
    return x ** atom[1]


def atom_d1(atom, x):
    kind = atom[0]
    if kind == "const":
        return 0.0
    if kind == "linear":
        return 1.0
    if kind == "exp":
        return math.exp(x)
    if kind == "log":
        return 1.0 / x
    a = atom[1]
    return a * x ** (a - 1.0)


def atom_d2_sign(atom, coeff):
    kind = atom[0]
    if kind in ("const", "linear"):
        return 0
    if kind == "exp":
        raw = 1.0
    elif kind == "log":
        raw = -1.0
    else:
        a = atom[1]
        raw = a * (a - 1.0)
    val = coeff * raw
    if val > EPS:
        return 1
    if val < -EPS:
        return -1
    return 0


def evaluate(terms, x):
    total = 0
    for c, a in terms:
        total += c * atom_value(a, x)
    return total


def derivative(terms, x):
    total = 0
    for c, a in terms:
        total += c * atom_d1(a, x)
    return total


def curvature(terms):
    signs = {atom_d2_sign(a, c) for c, a in terms}
    signs.discard(0)
    if not signs:
        return "affine"
    if signs == {1}:
        return "convex"
    if signs == {-1}:
        return "concave"
    return "mixed"


def minimize_convex(terms, domain):
    lo, hi = domain
    lo = max(lo, 1e-8)
    if hi is not None:
        if derivative(terms, lo) >= 0:
            return lo, evaluate(terms, lo)
        if derivative(terms, hi) <= 0:
            return hi, evaluate(terms, hi)
        a, b = lo, hi
    else:
        if derivative(terms, lo) >= 0:
            return lo, evaluate(terms, lo)
        a, b = lo, 1.0
        while derivative(terms, b) <= 0 and b < 1e8:
            b *= 2.0
        if b >= 1e8 and derivative(terms, b) <= 0:
            return b, evaluate(terms, b)
    for _ in range(160):
        mid = (a + b) / 2.0
        if derivative(terms, mid) <= 0:
            a = mid
        else:
            b = mid
    x = (a + b) / 2.0
    return x, evaluate(terms, x)


def line_gap_terms(base_terms, m, b, sign):
    line = [(m, ("linear",)), (b, ("const",))]
    if sign == 1:
        return combine_terms(line + [(-c, a) for c, a in base_terms])
    return combine_terms([*base_terms, (-m, ("linear",)), (-b, ("const",))])


def continued_fraction_convergents(x, max_terms=12, max_denominator=10000):
    terms = []
    y = x
    for _ in range(max_terms):
        a = math.floor(y)
        terms.append(a)
        frac_part = y - a
        if abs(frac_part) < 1e-14:
            break
        y = 1.0 / frac_part
    convergents = []
    for i in range(1, len(terms) + 1):
        value = Fraction(terms[i - 1], 1)
        for a in reversed(terms[: i - 1]):
            value = a + Fraction(1, value)
        if value.denominator <= max_denominator and value not in convergents:
            convergents.append(value)
    return convergents


def tangent_line_at(terms, x0):
    m = derivative(terms, x0)
    b = evaluate(terms, x0) - m * x0
    return m, b


def find_line(left, right, domain, x_center, strict):
    """-> (tangent Fraction, x_gap, min_gap) or None."""
    for rational in continued_fraction_convergents(x_center):
        x0 = float(rational)
        if x0 <= domain[0] or (domain[1] is not None and x0 >= domain[1]):
            continue
        m, b = tangent_line_at(right, x0)
        upper_gap = line_gap_terms(left, m, b, -1)
        x_gap, min_gap = minimize_convex(upper_gap, domain)
        line_proves = min_gap > 1e-8 if strict else min_gap >= -1e-7
        if line_proves:
            return rational, x_gap, min_gap
    return None


def parse_domain(text):
    if not text:
        return (1e-8, None)
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 2:
        raise SiteErr("domain must look like '0,inf' or '0,10'")
    lo = float(parts[0])
    hi = None if parts[1].lower() in {"inf", "infinity", "+inf"} else float(parts[1])
    return (lo, hi)


def affine_line_from_expr(text):
    m = 0.0
    b = 0.0
    for c, a in combine_terms(parse_terms(text)):
        if a[0] == "linear":
            m += c
        elif a[0] == "const":
            b += c
        else:
            raise SiteErr("--line must be an affine expression in x")
    return m, b


def verify_line(left, right, domain, m, b):
    left_gap = line_gap_terms(left, m, b, -1)
    right_gap = line_gap_terms(right, m, b, 1)
    x_left, min_left = minimize_convex(left_gap, domain)
    x_right, min_right = minimize_convex(right_gap, domain)
    return {
        "m": m,
        "b": b,
        "left_gap_min_x": x_left,
        "left_gap_min": min_left,
        "right_gap_min_x": x_right,
        "right_gap_min": min_right,
        "ok": min_left >= -1e-8 and min_right >= -1e-8,
    }


REASON_INCONCLUSIVE = (
    "整理后左侧不是凸函数/仿射函数，或右侧不是凹函数/仿射函数，因此当前证明器无法处理。"
)
REASON_FAILED = "数值最小值未达到证明要求；该不等式可能不成立，或超出当前搜索范围。"
REASON_NO_LINE = "不等式数值上已通过，但当前情形没有生成中间直线证明。"
REASON_SEARCH_MISS = "不等式数值上已通过，但内置有限候选搜索没有找到漂亮的有理切点直线。"


def predict(ineq, domain=None, line=None):
    """Route checks + prove() flow -> status/reason/min_x/min_v/tangent/gap/pline."""
    stripped = (ineq or "").strip()
    if not stripped:
        raise SiteErr("请输入一个不等式。")
    if len(stripped) > 500:
        raise SiteErr("输入过长，请输入一个较短的单变量不等式")
    dom = parse_domain(domain or "")
    lhs, rhs, op = split_inequality(ineq)
    f_terms = combine_terms(parse_terms(lhs) + [(-c, a) for c, a in parse_terms(rhs)])
    left, right = split_positive_negative(f_terms)
    diff = combine_terms(left + [(-c, a) for c, a in right])
    left_curv = curvature(left)
    right_curv = curvature(right)
    diff_curv = curvature(diff)
    # normalized/*_latex render Fraction(coeff/exponent) before the template
    # check; non-finite floats overflow there -> 'cannot convert ... ratio'.
    for c, a in diff:
        if not math.isfinite(c):
            raise SiteErr(
                "cannot convert Infinity to integer ratio"
                if math.isinf(c)
                else "cannot convert NaN to integer ratio"
            )
        if a[0] == "power" and not math.isfinite(a[1]):
            raise SiteErr("cannot convert Infinity to integer ratio")
    out = {
        "status": None,
        "reason": None,
        "min_x": None,
        "min_v": None,
        "tangent": None,
        "gap": None,
        "pline": None,
        "curv_left": left_curv,
        "curv_right": right_curv,
        "curv_diff": diff_curv,
    }
    template_ok = left_curv in {"convex", "affine"} and right_curv in {"concave", "affine"}
    if not template_ok or diff_curv not in {"convex", "affine"}:
        out["status"] = "inconclusive"
        out["reason"] = REASON_INCONCLUSIVE
        return out
    x_min, f_min = minimize_convex(diff, dom)
    out["min_x"], out["min_v"] = x_min, f_min
    strict = op in {">", "<"}
    ok = f_min > 1e-8 if strict else f_min >= -1e-8
    if not ok:
        out["status"] = "failed"
        out["reason"] = REASON_FAILED
        return out
    out["status"] = "proved"
    if line:
        m, b = affine_line_from_expr(line)
        out["pline"] = verify_line(left, right, dom, m, b)
    if left_curv != "convex":
        out["reason"] = REASON_NO_LINE
        return out
    found = find_line(left, right, dom, x_min, strict)
    if found is None:
        out["reason"] = REASON_SEARCH_MISS
    else:
        rational, _x_gap, min_gap = found
        out["tangent"] = str(rational) if rational.denominator != 1 else str(rational.numerator)
        out["gap"] = min_gap
        out["reason"] = ""
    return out
