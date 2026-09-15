# ruff: noqa: RUF002, RUF003
r"""第二波边缘探针：钉死 /get_integral_image 的校验规则和 /calculate 的校验顺序。

第一波发现站端 image 端点不是笼统 500，而是逐字段 400 校验
（m必须是整数 / m过小 / m过大 / u_val过小 / a_val格式无效 / a_val分母不能为0 /
无效的证明类型 / 无效的不等号方向），且 type/comparison/coef/rational
有缺省值、coef 与 rational 原文不校验不规范化直接回显。

本波要钉死：m/n 的上下界（是否随类型 cap 变化）、u_val 下界是否因类型而异、
a/b/c 的格式域（Fraction() 直解？接受小数/科学计数法？）、各校验的先后次序、
以及 /calculate 中 comparison 与右侧格式的相对次序。

记录格式与 edge_probe.py 相同，追加进同一个 edge-probes.jsonl。
"""
import argparse
import json
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "edge-probes.jsonl"
MIN_INTERVAL = 1.05
TIMEOUT = 100

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (edge probes w2)"

IMG = {"m": "3", "n": "3", "a_val": "47/120", "b_val": "-13/120",
       "c_val": "0", "u_val": "120", "au_val": "47", "bu_val": "-13",
       "cu_val": "0", "type": "pi", "coef": "1", "comparison": "<",
       "rational": "22/7"}

PROBES: list[dict] = []


def calc(pid, note, **form):
    fields = {k: v for k, v in form.items() if v is not None}
    PROBES.append({"id": pid, "endpoint": "/calculate", "method": "POST",
                   "kind": "form", "request": fields, "note": note})


def img(pid, note, **query):
    fields = dict(IMG)
    fields.update(query)
    fields = {k: v for k, v in fields.items() if v is not None}
    PROBES.append({"id": pid, "endpoint": "/get_integral_image",
                   "method": "GET", "kind": "query",
                   "request": fields, "note": note})


def raw(pid, endpoint, method, kind, request, note):
    PROBES.append({"id": pid, "endpoint": endpoint, "method": method,
                   "kind": kind, "request": request, "note": note})


# -------------------------------------------- image: m/n 界与整数字段
img("img2:m30-pi", "pi m=30（cap 边界）", m="30")
img("img2:m31-pi", "pi m=31：过大阈值=30?", m="31")
img("img2:m100-pi", "pi m=100", m="100")
img("img2:m11-atan", "arctan_q m=11：阈值=类型 cap?", type="arctan_q",
    m="11", coef="1", rational="1")
img("img2:m10-atan", "arctan_q m=10", type="arctan_q", m="10", coef="1",
    rational="1")
img("img2:m33-e", "e m=33（若阈值=30 则过大）", type="e", m="33",
    rational="3")
img("img2:m31-e", "e m=31", type="e", m="31", rational="3")
img("img2:n-neg", "n=-1", n="-1")
img("img2:n99-pi", "pi n=99", n="99")
img("img2:n31-pi", "pi n=31", n="31")
img("img2:n11-atan", "arctan_q n=11", type="arctan_q", n="11", coef="1",
    rational="1")
img("img2:m-lead0", "m=03 前导零", m="03")
img("img2:m-plus", "m=+3", m="+3")
img("img2:m-space", "m=' 3 ' 空格", m=" 3 ")
img("img2:m-huge", "m=10^17", m="100000000000000000")
img("img2:n-alpha", "n=abc", n="abc")
img("img2:n-missing", "n 缺失", n=None)

# -------------------------------------------- image: u 与系数字段
img("img2:u0-e", "e 型 u_val=0（gamma 退化可过，e 呢）", type="e", u_val="0",
    rational="3")
img("img2:u0-gamma-plain", "gamma 型 u_val=0 普通参数", type="gamma",
    u_val="0", m="2", n="4", a_val="937/56", b_val="4615/42",
    au_val="2811", bu_val="18460", cu_val="5", coef="1", rational="3/5")
img("img2:u-abc", "u_val=abc", u_val="abc")
img("img2:u-space", "u_val=' 120 '", u_val=" 120 ")
img("img2:u-huge", "u_val=10^20", u_val="100000000000000000000")
img("img2:au-abc", "au_val=abc", au_val="abc")
img("img2:au-float", "au_val=4.5", au_val="4.5")
img("img2:au-missing", "au_val 缺失", au_val=None)
img("img2:bu-missing", "bu_val 缺失", bu_val=None)
img("img2:a-missing", "a_val 缺失", a_val=None)
img("img2:a-dec", "a_val=3.5（Fraction 直解则接受）", a_val="3.5")
img("img2:a-sci", "a_val=1e2（Fraction 接受科学计数）", a_val="1e2")
img("img2:a-space", "a_val=' 47/120 ' 空格", a_val=" 47/120 ")
img("img2:b-alpha", "b_val=abc", b_val="abc")
img("img2:c-alpha", "c_val=abc", c_val="abc")
img("img2:c-dec", "c_val=0.5", c_val="0.5")
img("img2:a-neg", "a_val=-47/120 负", a_val="-47/120")

# -------------------------------------------- image: type/comp/coef/rational
img("img2:type-empty", "type='' 显式空串", type="")
img("img2:comp-empty", "comparison='' 显式空串", comparison="")
img("img2:comp-x", "comparison='x'", comparison="x")
img("img2:coef-11", "coef='1/1'：字符串==1 判等还是解析", coef="1/1")
img("img2:coef-empty", "coef='' 显式空串", coef="")
img("img2:coef-space", "coef=' 1 '", coef=" 1 ")
img("img2:coef-latex", r"coef='\frac{2}{4}' LaTeX 形", coef="\\frac{2}{4}")
img("img2:rat-empty", "rational='' 显式空串", rational="")
img("img2:rat-word", "rational=abc（界是否也免解析）", rational="abc")
img("img2:rat-latex-d", r"rational='\dfrac{22}{7}'", rational="\\dfrac{22}{7}")
img("img2:prec-type-m", "type=foo&m=abc：校验先后", type="foo", m="abc")
img("img2:prec-m-comp", "m=abc&comp='='：m 校验 vs comp", m="abc",
    comparison="=")
img("img2:prec-comp-a", "comp='='&a_val=abc：comp vs a_val", comparison="=",
    a_val="abc")

# -------------------------------------------- /calculate 校验次序
calc("ord:comp-eq-bad-rat", "comp='=' + rational=abc：谁先",
     type="pi", power="1", comparison="=", rational="abc")
calc("ord:comp-eq-bad-pow", "comp='=' + power=x + rational 合法",
     type="pi", power="x", comparison="=", rational="22/7")
calc("ord:comp-eq-no-rat", "comp='=' + rational 缺失", type="pi", power="1",
     comparison="=", rational=None)
calc("ord:type-foo-no-comp", "type=foo + comp 缺失", type="foo", power="1",
     comparison=None, rational="4")
calc("ord:no-comp-bad-rat", "comp 缺失 + rational=abc", type="pi", power="1",
     comparison=None, rational="abc")
calc("ord:bad-type-bad-rat", "type=foo + rational=abc", type="foo", power="1",
     comparison="<", rational="abc")
calc("ord:good-comp-bad-pow", "comp 合法 + power=x + rational=abc",
     type="pi", power="x", comparison="<", rational="abc")

# -------------------------------------------- /calculate 空白与特殊格式
calc("fmt:rat-inner-both", "rational='22 / 7' 两侧空格", type="pi", power="1",
     comparison="<", rational="22 / 7")
calc("fmt:rat-inner-digit", "rational='2 2/7' 数字内空格", type="pi",
     power="1", comparison="<", rational="2 2/7")
calc("fmt:pow-inner", "power='1 / 2' 斜杠两侧空格", type="e_q",
     power="1 / 2", comparison="<", rational="2")
calc("fmt:rat-tab", "rational 前导 tab", type="pi", power="1", comparison="<",
     rational="\t22/7")
calc("fmt:rat-newline", "rational 尾随换行", type="pi", power="1",
     comparison="<", rational="22/7\n")
calc("fmt:rat-2slash", "rational='22/7/1' 双斜杠", type="pi", power="1",
     comparison="<", rational="22/7/1")
calc("fmt:rat-slash-end", "rational='22/'", type="pi", power="1",
     comparison="<", rational="22/")
calc("fmt:rat-slash-start", "rational='/7'", type="pi", power="1",
     comparison="<", rational="/7")
calc("fmt:pow-frac-neg", "power='-1/2' 负分数", type="e_q", power="-1/2",
     comparison="<", rational="2")
calc("fmt:pow-fullwidth-slash", "power='１／２' 全角斜杠", type="e_q",
     power="１／２", comparison="<", rational="2")

# -------------------------------------------- 残余行为
calc("beh:gamma0-gt0", "gamma 0 > 0：0>0 假命题，站端先判方向再除零?",
     type="gamma", power="0", comparison=">", rational="0")
calc("beh:gamma0-lt0b", "gamma 0 < 0：0<0 假", type="gamma", power="0",
     comparison="<", rational="0")
calc("beh:sinh0-lt-big", "sinh 0 < 9999999999999999", type="sinh_q", power="0",
     comparison="<", rational="9999999999999999")
calc("beh:e-neg-coef", "e power=-1 负系数（格式层）", type="e", power="-1",
     comparison="<", rational="1")

raw("mth2:put-calculate", "/calculate", "PUT", "form",
    {"type": "pi", "power": "1", "comparison": "<", "rational": "22/7"},
    "PUT /calculate：站端对其它 method 的行为")
raw("mth2:delete-image", "/get_integral_image", "DELETE", "query", IMG,
    "DELETE /get_integral_image")
raw("mth2:options-calculate", "/calculate", "OPTIONS", "form", {},
    "OPTIONS /calculate")
raw("mth2:get-image-slash", "/get_integral_image/", "GET", "query", IMG,
    "GET /get_integral_image/ 尾斜杠")
raw("mth2:dec-empty", "/decompose_inequality", "POST", "form", {},
    "decompose 缺 problem 字段")
raw("mth2:dec-emptystr", "/decompose_inequality", "POST", "form",
    {"problem": ""}, "decompose problem=''")
raw("mth2:dec-spaces", "/decompose_inequality", "POST", "form",
    {"problem": "   "}, "decompose problem 全空格")
raw("mth2:dec-nosolve", "/decompose_inequality", "POST", "form",
    {"problem": "foo>3"}, "decompose 无法解析的输入")
raw("page2:en-query", "/en?x=1", "GET", "none", {}, "/en 带 query")
raw("page2:health", "/health", "GET", "none", {}, "姊妹页面 /health")
raw("page2:convex", "/convex", "GET", "none", {}, "姊妹页面 /convex")


def site_request(p):
    if p["method"] == "GET":
        resp = SESSION.get(BASE + p["endpoint"], params=p["request"] or None,
                           timeout=TIMEOUT)
    elif p["method"] == "PUT":
        resp = SESSION.put(BASE + p["endpoint"], data=p["request"],
                           timeout=TIMEOUT)
    elif p["method"] == "DELETE":
        resp = SESSION.delete(BASE + p["endpoint"], params=p["request"],
                              timeout=TIMEOUT)
    elif p["method"] == "OPTIONS":
        resp = SESSION.options(BASE + p["endpoint"], timeout=TIMEOUT)
    elif p["kind"] == "query":
        resp = SESSION.post(BASE + p["endpoint"], params=p["request"],
                            timeout=TIMEOUT)
    else:
        resp = SESSION.post(BASE + p["endpoint"], data=p["request"],
                            timeout=TIMEOUT)
    return resp.status_code, resp.text


def ours_request(client, p):
    if p["method"] == "GET":
        resp = client.get(p["endpoint"], query_string=p["request"] or None)
    elif p["method"] == "PUT":
        resp = client.put(p["endpoint"], data=p["request"])
    elif p["method"] == "DELETE":
        resp = client.delete(p["endpoint"], query_string=p["request"])
    elif p["method"] == "OPTIONS":
        resp = client.options(p["endpoint"])
    elif p["kind"] == "query":
        resp = client.post(p["endpoint"], query_string=p["request"])
    else:
        resp = client.post(p["endpoint"], data=p["request"])
    return resp.status_code, resp.get_data(as_text=True)


def bodies_equal(site_raw, ours_raw):
    try:
        s, o = json.loads(site_raw), json.loads(ours_raw)
    except (json.JSONDecodeError, TypeError):
        return site_raw.strip() == ours_raw.strip()
    return s == o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    args = ap.parse_args()

    from attention_calculator.server import app
    client = app.test_client()

    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                try:
                    done.add(json.loads(line)["id"])
                except json.JSONDecodeError:
                    pass
    todo = [p for p in PROBES if p["id"] not in done]
    if args.only:
        todo = [p for p in todo if args.only in p["id"]]
    print(f"{len(done)} done, {len(todo)} to probe", flush=True)

    last = 0.0
    for i, p in enumerate(todo, 1):
        wait = MIN_INTERVAL - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
        t0 = time.monotonic()
        try:
            status, raw_body = site_request(p)
        except (OSError, requests.RequestException) as exc:
            status, raw_body = -1, f"REQUEST_FAILED: {exc}"
        elapsed = time.monotonic() - t0
        last = time.monotonic()
        try:
            ostatus, obody = ours_request(client, p)
        except Exception as exc:
            ostatus, obody = -2, f"LOCAL_CRASH: {type(exc).__name__}: {exc}"
        rec = dict(p)
        rec.update({
            "http_status": status, "raw": raw_body,
            "elapsed_ms": round(elapsed * 1000, 1),
            "ours_status": ostatus, "ours_body": obody,
            "byte_match": raw_body == obody,
            "match": status == ostatus and bodies_equal(raw_body, obody),
        })
        with OUT.open("a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        flag = "" if rec["match"] else "  <<< MISMATCH"
        print(f"[{i}/{len(todo)}] {p['id']}: site={status} ours={ostatus}{flag}",
              flush=True)


if __name__ == "__main__":
    main()
