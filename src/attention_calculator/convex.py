"""POST /convex/prove solver — byte-exact clone of zhuyidao.net 凹凸不等式计算器.

Core = lianghuatiaojiushi/ConvexConcaveProver scripts/prove.py verbatim
(float64 term model, derivative-sign bisection with 160 halvings, CF tangent
candidates 12 terms/denominator<=10000) + the site's JSON/latex wrapper:
normalized/*_latex fields, status/reason strings, line/domain params.
"""

import ast
import math
import re
from fractions import Fraction

EPS = 1e-9

REASON_INCONCLUSIVE = (
    "整理后左侧不是凸函数/仿射函数，或右侧不是凹函数/仿射函数，因此当前证明器无法处理。"
)
REASON_FAILED = "数值最小值未达到证明要求；该不等式可能不成立，或超出当前搜索范围。"
REASON_NO_LINE = "不等式数值上已通过，但当前情形没有生成中间直线证明。"
REASON_SEARCH_MISS = "不等式数值上已通过，但内置有限候选搜索没有找到漂亮的有理切点直线。"

# atom = ("const",) | ("linear",) | ("exp",) | ("log",) | ("power", float)
# term = (coeff: float, atom)


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
            raise ValueError("sqrt expects one argument")
        return math.sqrt(parse_number(node.args[0]))
    raise ValueError(f"expected numeric constant, got {dump(node)}")


def parse_atom(node):
    if isinstance(node, ast.Name) and node.id == "x":
        return ("linear",)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
        name = node.func.id
        if len(node.args) != 1:
            raise ValueError(f"{name} expects one argument")
        arg = node.args[0]
        if not (isinstance(arg, ast.Name) and arg.id == "x"):
            raise ValueError(f"{name} only supports argument x")
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
    raise ValueError(f"unsupported atom {dump(node)}")


def parse_factor_term(node):
    try:
        return (parse_number(node), ("const",))
    except ValueError:
        return (1.0, parse_atom(node))


def parse_product(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        left = parse_product(node.left)
        right = parse_product(node.right)
        if left[1][0] == "const":
            return (right[0] * left[0], right[1])
        if right[1][0] == "const":
            return (left[0] * right[0], left[1])
        raise ValueError("products of two non-constant atoms are unsupported")
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
        except ValueError:
            pass
        try:
            factor = parse_number(node.right)
            return collect_terms(node.left, sign * factor)
        except ValueError:
            pass
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        try:
            factor = parse_number(node.right)
            return collect_terms(node.left, sign / factor)
        except ValueError:
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


def evaluate(terms, x):
    # 站端跑旧 CPython（sum 为朴素顺序累加）；3.12+ 的 sum 改 Neumaier 补偿求和，
    # 末位 ulp 会差——手写循环复刻旧行为
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


def fraction_label(frac):
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def fraction_latex(frac):
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"\\frac{{{frac.numerator}}}{{{frac.denominator}}}"


def signed_fraction_label(frac):
    if frac == 0:
        return ""
    sign = "+" if frac > 0 else "-"
    return f" {sign} {fraction_label(abs(frac))}"


def signed_fraction_latex(frac):
    if frac == 0:
        return ""
    sign = "+" if frac > 0 else "-"
    return f" {sign} {fraction_latex(abs(frac))}"


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


def tangent_term_text(coeff, atom, point):
    c = Fraction(coeff).limit_denominator(1000000)
    r = fraction_label(point)
    prefix = "" if c == 1 else "-" if c == -1 else f"{fraction_label(c)}*"
    kind = atom[0]
    if kind == "const":
        return fraction_label(c)
    if kind == "linear":
        return f"{prefix}x"
    if kind == "exp":
        shift = signed_fraction_label(1 - point)
        return f"{prefix}e^({r})(x{shift})"
    if kind == "log":
        reciprocal = fraction_label(1 / point)
        return f"{prefix}({reciprocal}*x - 1 + ln({r}))"
    a = atom[1]
    if abs(a - 0.5) < 1e-12:
        return f"{prefix}((x + {r})/(2*sqrt({r})))"
    return f"{prefix}(({r})^{a:g} + {a:g}*({r})^({a:g}-1)*(x - {r}))"


def tangent_term_latex(coeff, atom, point):
    c = Fraction(coeff).limit_denominator(1000000)
    r = fraction_latex(point)
    prefix = "" if c == 1 else "-" if c == -1 else fraction_latex(c)
    kind = atom[0]
    if kind == "const":
        return fraction_latex(c)
    if kind == "linear":
        return f"{prefix}x"
    if kind == "exp":
        shift = signed_fraction_latex(1 - point)
        return f"{prefix}e^{{{r}}}(x{shift})"
    if kind == "log":
        reciprocal = fraction_latex(1 / point)
        return f"{prefix}({reciprocal}x - 1 + \\ln {r})"
    a = Fraction(atom[1]).limit_denominator(1000000)
    a_tex = fraction_latex(a)
    if abs(atom[1] - 0.5) < 1e-12:
        return f"{prefix}\\frac{{x + {r}}}{{2\\sqrt{{{r}}}}}"
    return f"{prefix}(({r})^{{{a_tex}}} + {a_tex}({r})^{{{a_tex}-1}}(x - {r}))"


def tangent_expr(terms, point, tex):
    fn = tangent_term_latex if tex else tangent_term_text
    return " + ".join(fn(c, a, point) for c, a in terms).replace("+ -", "- ")


def find_line(left, right, domain, x_center, strict, left_tex, right_tex):
    for rational in continued_fraction_convergents(x_center):
        x0 = float(rational)
        if x0 <= domain[0] or (domain[1] is not None and x0 >= domain[1]):
            continue
        m, b = tangent_line_at(right, x0)
        upper_gap = line_gap_terms(left, m, b, -1)
        x_gap, min_gap = minimize_convex(upper_gap, domain)
        line_proves = min_gap > 1e-8 if strict else min_gap >= -1e-7
        if line_proves:
            line_text = tangent_expr(right, rational, False)
            line_latex = tangent_expr(right, rational, True)
            rel = ">" if strict and min_gap > 1e-8 else "\\ge"
            return {
                "formula_latex": f"{left_tex} {rel} {line_latex} \\ge {right_tex}",
                "left_gap_min": min_gap,
                "left_gap_min_x": x_gap,
                "line_latex": line_latex,
                "line_text": line_text,
                "tangent_at": fraction_label(rational),
                "tangent_at_latex": fraction_latex(rational),
            }
    return None


def parse_domain(text):
    if not text:
        return (1e-8, None)
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 2:
        raise ValueError("domain must look like '0,inf' or '0,10'")
    lo = float(parts[0])
    hi = None if parts[1].lower() in {"inf", "infinity", "+inf"} else float(parts[1])
    return (lo, hi)


def atom_text(atom):
    kind = atom[0]
    if kind == "linear":
        return "x"
    if kind == "exp":
        return "e^x"
    if kind == "log":
        return "ln x"
    if kind == "power":
        if abs(atom[1] - 0.5) < 1e-12:
            return "√x"
        return f"x^{atom[1]:g}"
    return ""


def atom_latex(atom):
    kind = atom[0]
    if kind == "linear":
        return "x"
    if kind == "exp":
        return "e^x"
    if kind == "log":
        return "\\ln x"
    if kind == "power":
        if abs(atom[1] - 0.5) < 1e-12:
            return "\\sqrt{x}"
        e = Fraction(atom[1]).limit_denominator(1000000)
        body = (
            str(e.numerator)
            if e.denominator == 1
            else f"\\frac{{{e.numerator}}}{{{e.denominator}}}"
        )
        return f"x^{{{body}}}"
    return ""


def fmt_terms(terms, tex):
    pieces = []
    for coeff_f, atom in terms:
        coeff = Fraction(coeff_f).limit_denominator(1000000)
        if coeff == 0:
            continue
        sign = "-" if coeff < 0 else "+"
        mag = abs(coeff)
        if tex:
            c_str = fraction_latex(mag)
            a_str = atom_latex(atom)
            joiner = ""
        else:
            c_str = fraction_label(mag)
            a_str = atom_text(atom)
            joiner = "*"
        if atom[0] == "const":
            body = c_str
        elif mag == 1:
            body = a_str
        else:
            body = f"{c_str}{joiner}{a_str}"
        pieces.append((sign, body))
    if not pieces:
        return "0"
    first_sign, first_body = pieces[0]
    out = first_body if first_sign == "+" else f"-{first_body}"
    for sign, body in pieces[1:]:
        out += f" {sign} {body}"
    return out


def affine_line_from_expr(text):
    m = 0.0
    b = 0.0
    for c, a in combine_terms(parse_terms(text)):
        if a[0] == "linear":
            m += c
        elif a[0] == "const":
            b += c
        else:
            raise ValueError("--line must be an affine expression in x")
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


def prove(inequality, line=None, domain=None):
    """POST /convex/prove 的 result 字典；异常由路由层映射 400。"""
    dom = parse_domain(domain or "")
    lhs, rhs, op = split_inequality(inequality)
    f_terms = combine_terms(parse_terms(lhs) + [(-c, a) for c, a in parse_terms(rhs)])
    left, right = split_positive_negative(f_terms)
    diff = combine_terms(left + [(-c, a) for c, a in right])
    left_curv = curvature(left)
    right_curv = curvature(right)
    diff_curv = curvature(diff)
    result = {
        "curvature": {
            "difference": diff_curv,
            "left": left_curv,
            "right": right_curv,
        },
        "minimum": None,
        "normalized": {
            "difference": fmt_terms(diff, False),
            "difference_latex": fmt_terms(diff, True),
            "left": fmt_terms(left, False),
            "left_latex": fmt_terms(left, True),
            "right": fmt_terms(right, False),
            "right_latex": fmt_terms(right, True),
        },
        "ok": False,
        "proof": None,
        "provided_line": None,
        "reason": "",
        "status": "",
    }
    template_ok = left_curv in {"convex", "affine"} and right_curv in {
        "concave",
        "affine",
    }
    if not template_ok or diff_curv not in {"convex", "affine"}:
        result["status"] = "inconclusive"
        result["reason"] = REASON_INCONCLUSIVE
        return result
    x_min, f_min = minimize_convex(diff, dom)
    result["minimum"] = {
        "value": f_min,
        "value_text": f"{f_min:.12g}",
        "x": x_min,
        "x_text": f"{x_min:.12g}",
    }
    strict = op in {">", "<"}
    ok = f_min > 1e-8 if strict else f_min >= -1e-8
    if not ok:
        result["status"] = "failed"
        result["reason"] = REASON_FAILED
        return result
    result["status"] = "proved"
    result["ok"] = True
    if line:
        m, b = affine_line_from_expr(line)
        result["provided_line"] = verify_line(left, right, dom, m, b)
    if left_curv != "convex":
        result["reason"] = REASON_NO_LINE
        return result
    proof = find_line(
        left,
        right,
        dom,
        x_min,
        strict,
        result["normalized"]["left_latex"],
        result["normalized"]["right_latex"],
    )
    if proof is None:
        result["reason"] = REASON_SEARCH_MISS
    else:
        result["proof"] = proof
    return result
