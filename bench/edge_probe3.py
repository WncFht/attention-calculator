r"""第三波边缘探针：核行为消歧（fuzz/edge 两波遗留的内核分歧点）。

目标：
1. gauss/varpi '<' 退化解（cu_val=2）的触发边界与排序——b=0 是否严格大于才发；
2. 双曲核 q=0 的方向预检模型（math.* 溢出即 500）与其它类型 q=0 假命题；
3. gamma u=0 退化渲染的尾项公式与 coef 乘数形态；
4. pi_n 分数/越界 coef 的图像渲染（ln 指数取 numerator-1 还是 denominator）；
5. ln_q 负/未约分 coef 的分母构造（原始 |num|/den 不约分？）；
6. e_q 负指数的符号、tan_q 负参数的因子排布；
7. coef_tex 规范化规则在各类型常数位的适用面。

记录格式与 edge_probe.py 相同，追加进同一个 edge-probes.jsonl。
"""

import argparse
import contextlib
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
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (edge probes w3)"

IMG = {
    "m": "3",
    "n": "3",
    "a_val": "47/120",
    "b_val": "-13/120",
    "c_val": "0",
    "u_val": "120",
    "au_val": "47",
    "bu_val": "-13",
    "cu_val": "0",
    "type": "pi",
    "coef": "1",
    "comparison": "<",
    "rational": "22/7",
}

PROBES: list[dict] = []


def calc(pid, note, **form):
    fields = {k: v for k, v in form.items() if v is not None}
    PROBES.append(
        {
            "id": pid,
            "endpoint": "/calculate",
            "method": "POST",
            "kind": "form",
            "request": fields,
            "note": note,
        }
    )


def img(pid, note, **query):
    fields = dict(IMG)
    fields.update(query)
    fields = {k: v for k, v in fields.items() if v is not None}
    PROBES.append(
        {
            "id": pid,
            "endpoint": "/get_integral_image",
            "method": "GET",
            "kind": "query",
            "request": fields,
            "note": note,
        }
    )


# -------------------------------------------- 1. gauss/varpi '<' 退化解
calc(
    "k3:gauss-b0-gt",
    "gauss 1<6/5: 退化 b=0 边界（b>0 严格还是 >=0）",
    type="gauss",
    power="1",
    comparison="<",
    rational="6/5",
)
calc(
    "k3:varpi-b0-gt",
    "varpi 1<7/2: 退化 b=0 边界",
    type="varpi",
    power="1",
    comparison="<",
    rational="7/2",
)
calc(
    "k3:gauss-degen-54",
    "gauss 1<5/4: b=1/24>0 应发退化 au=25,u=4",
    type="gauss",
    power="1",
    comparison="<",
    rational="5/4",
)
calc(
    "k3:varpi-degen-154",
    "varpi 1<15/4: b=1/4>0 应发退化 au=21,u=1",
    type="varpi",
    power="1",
    comparison="<",
    rational="15/4",
)
calc(
    "k3:gauss-0-lt-0",
    "gauss 0<0: 退化 b=0 跳过→转置 lhs_mpf 0/0",
    type="gauss",
    power="0",
    comparison="<",
    rational="0",
)
calc(
    "k3:varpi-0-lt-0",
    "varpi 0<0: 退化 b=0 跳过→转置/未找到",
    type="varpi",
    power="0",
    comparison="<",
    rational="0",
)
calc(
    "k3:gauss-gt-14",
    "gauss 1>1/4: '>' 退化 au=3,u=1,b=1/8?",
    type="gauss",
    power="1",
    comparison=">",
    rational="1/4",
)
calc(
    "k3:varpi-gt-32",
    "varpi 1>3/2: '>' 退化 au=15,u=2,b=1/16?",
    type="varpi",
    power="1",
    comparison=">",
    rational="3/2",
)
calc(
    "k3:gauss-degen-2-3",
    "gauss 2<3: b=1/2>0 退化 au=15,u=1",
    type="gauss",
    power="2",
    comparison="<",
    rational="3",
)

# -------------------------------------------- 2. 双曲核 q=0 / 大 q
calc(
    "k3:cosh0-lt-half",
    "cosh 0<1/2 假: 方向预检存在→404 方向反了?",
    type="cosh_q",
    power="0",
    comparison="<",
    rational="1/2",
)
calc(
    "k3:cosh0-gt-half",
    "cosh 0>1/2 真: 过预检→核崩 500?",
    type="cosh_q",
    power="0",
    comparison=">",
    rational="1/2",
)
calc(
    "k3:coth0-gt1",
    "coth 0>1: 1/tanh(0) 除零→500",
    type="coth_q",
    power="0",
    comparison=">",
    rational="1",
)
calc(
    "k3:sinh700-gt1",
    "sinh700≈1.1e308 有限: 真方向→核路径结果",
    type="sinh_q",
    power="700",
    comparison=">",
    rational="1",
)
calc(
    "k3:sinh711-gt1",
    "sinh711 math 溢出→500?",
    type="sinh_q",
    power="711",
    comparison=">",
    rational="1",
)
calc(
    "k3:cosh-big-gt0",
    "cosh 1e15>0: 预检用 math→500; 仅崩后判→404",
    type="cosh_q",
    power="1000000000000000/1",
    comparison=">",
    rational="0",
)
calc(
    "k3:tanh-big-gt1",
    "tanh 1e15>1: 饱和 1.0 相等→核路径",
    type="tanh_q",
    power="1000000000000000/1",
    comparison=">",
    rational="1",
)
calc(
    "k3:eq0-gt2",
    "e_q 0>2 假: 通用方向预检? 404 还是 500",
    type="e_q",
    power="0",
    comparison=">",
    rational="2",
)
calc(
    "k3:atan0-gt2", "arctan 0>2 假: 同上", type="arctan_q", power="0", comparison=">", rational="2"
)
calc(
    "k3:acot0-gt2",
    "arccot 0>2: 1/q 崩与预检谁先",
    type="arccot_q",
    power="0",
    comparison=">",
    rational="2",
)
calc(
    "k3:gamma3-lt1",
    "gamma 3<1 假: 预检→方向反了",
    type="gamma",
    power="3",
    comparison="<",
    rational="1",
)

# -------------------------------------------- 3. gamma 渲染
img(
    "k3:gamma-deg-tail",
    "u=0 尾项消歧: a_val=7/3 coef=2 bound=7/9 "
    "(a*coef=14/3 | coef*r0-bound=2/9 | bound=7/9 | a=7/3)",
    m="0",
    n="0",
    a_val="7/3",
    b_val="0",
    c_val="1",
    u_val="0",
    au_val="0",
    bu_val="0",
    cu_val="0",
    type="gamma",
    coef="2",
    comparison=">",
    rational="7/9",
)
img(
    "k3:gamma-deg-tail-lt",
    "'<' 版: a_val=7/3 c_val=2 coef=3 bound=7/9",
    m="0",
    n="0",
    a_val="7/3",
    b_val="0",
    c_val="2",
    u_val="0",
    au_val="0",
    bu_val="0",
    cu_val="0",
    type="gamma",
    coef="3",
    comparison="<",
    rational="7/9",
)
img(
    "k3:gamma-coef-main",
    "u≠0 正常参数 coef=2: main 的 coef 前缀形态",
    m="4",
    n="2",
    a_val="67",
    b_val="-304/5",
    c_val="0",
    au_val="335",
    bu_val="-304",
    cu_val="4",
    u_val="5",
    type="gamma",
    coef="2",
    comparison=">",
    rational="57/100",
)
img(
    "k3:gamma-coef-frac",
    "退化 coef=3/2: 分数乘数形态与尾项",
    m="0",
    n="0",
    a_val="1/4",
    b_val="0",
    c_val="1",
    u_val="0",
    au_val="0",
    bu_val="0",
    cu_val="0",
    type="gamma",
    coef="3/2",
    comparison=">",
    rational="1/2",
)

# -------------------------------------------- 4. pi_n 越界/分数 coef 渲染
PIN = {
    "type": "pi_n",
    "m": "2",
    "n": "0",
    "a_val": "197/108",
    "b_val": "1925/108",
    "c_val": "0",
    "au_val": "197",
    "bu_val": "1925",
    "cu_val": "0",
    "u_val": "108",
    "comparison": "<",
    "rational": "311/10",
}
img("k3:pin-neg34", "pi_n coef=-3/4: ln 指数 |num|-1=2 vs den=4", coef="-3/4", **PIN)
img("k3:pin-53", "pi_n coef=5/3: num-1=4 vs den=3", coef="5/3", **PIN)
img("k3:pin-neg53", "pi_n coef=-5/3: |num|-1=4 vs den=3 vs num-1=-6", coef="-5/3", **PIN)
img("k3:pin-24", "pi_n coef=2/4 未约分: den=4 还是约分后 den=2", coef="2/4", **PIN)
img("k3:pin-11", "pi_n coef=1/1: 整数化还是 (pi^{1/1})^1", coef="1/1", **PIN)
img("k3:pin-0", "pi_n coef=0: pi^0 渲染", coef="0", **PIN)

# -------------------------------------------- 5. ln_q 负/未约分 coef 分母
LNQ = {
    "type": "ln_q",
    "m": "6",
    "n": "4",
    "a_val": "32663/25",
    "b_val": "-5957/5",
    "c_val": "0",
    "au_val": "32663",
    "bu_val": "-29785",
    "cu_val": "0",
    "u_val": "25",
    "comparison": "<",
    "rational": "224/125",
}
img("k3:lnq-42", "ln_q coef=4/2 正未约分: (2x+2)^6 + t=64?", coef="4/2", **LNQ)
img("k3:lnq-neg64", "ln_q coef=-6/4: (2x+4)^6 + t=4096?", coef="-6/4", **LNQ)
img("k3:lnq-den-neg", "ln_q coef=4/-2 负分母", coef="4/-2", **LNQ)
img("k3:lnq-neg2", "ln_q coef=-2 整数负", coef="-2", **LNQ)

# -------------------------------------------- 6. e_q 负指数 / tan_q 负参数
EQ = {
    "type": "e_q",
    "m": "1",
    "n": "2",
    "a_val": "3/20",
    "b_val": "63/80",
    "c_val": "0",
    "au_val": "12",
    "bu_val": "63",
    "cu_val": "0",
    "u_val": "80",
    "comparison": ">",
    "rational": "22/5",
}
img("k3:eq-neg4", "e_q coef=-4: e^{4x} 还是 e^{-4x}", coef="-4", **EQ)
img("k3:eq-junk", "e_q coef=x 垃圾串", coef="x", **EQ)
img("k3:eq-negdfrac", "e_q coef=-\\dfrac{6}{4}", coef="-\\dfrac{6}{4}", **EQ)
TANQ = {
    "type": "tan_q",
    "m": "0",
    "n": "0",
    "a_val": "0",
    "b_val": "1",
    "c_val": "0",
    "au_val": "0",
    "bu_val": "1",
    "cu_val": "0",
    "u_val": "1",
    "comparison": ">",
    "rational": "1",
}
img(
    "k3:tanq-negdfrac",
    "tan_q coef=\\dfrac{-6}{4}: cos(-6/4)? -x sin(3x/2)?",
    coef="\\dfrac{-6}{4}",
    **TANQ,
)
img(
    "k3:tanq-posdfrac",
    "tan_q coef=\\dfrac{6}{4}: cos(6/4) sin(3x/2)?",
    coef="\\dfrac{6}{4}",
    **TANQ,
)
img("k3:tanq-junk", "tan_q coef=x 垃圾串", coef="x", **TANQ)
img("k3:tanq-1", "tan_q coef=1: \\tan1?", coef="1", **TANQ)

# -------------------------------------------- 7. coef_tex 适用面
img("k3:coef-22", "pi coef=2/2: 整数化 2 还是原文 \\dfrac{2}{2}", coef="2/2")
img("k3:coef-42", "pi coef=4/2", coef="4/2")
img("k3:coef-df42", "pi coef=\\dfrac{4}{2}: 宏解析后整数化?", coef="\\dfrac{4}{2}")
img("k3:coef-neg1", "pi coef=-1", coef="-1")
img("k3:coef-0", "pi coef=0", coef="0")
img("k3:rat-42", "pi rational=4/2: 界位也整数化?", rational="4/2")
img("k3:rat-11", "pi rational=1/1", rational="1/1")
img("k3:coef-e-11", "e coef=1/1", type="e", coef="1/1")
img("k3:coef-cat-11", "catalan coef=1/1", type="catalan", coef="1/1")
img("k3:coef-golden-11", "golden coef=1/1", type="golden", coef="1/1")
img("k3:coef-epi-11", "e_pi coef=1/1", type="e_pi", coef="1/1")
img("k3:coef-zeta3-11", "zeta3 coef=1/1", type="zeta3", coef="1/1")
img(
    "k3:varpi-gt-11-42",
    "varpi '>' coef=1/1 rational=4/2: 独立项位+ϖ^{-1}位",
    type="varpi",
    comparison=">",
    coef="1/1",
    rational="4/2",
)
img(
    "k3:gauss-lt-11-42",
    "gauss '<' coef=1/1 rational=4/2: G^{-1}位-独立项位",
    type="gauss",
    comparison="<",
    coef="1/1",
    rational="4/2",
)
img("k3:eq-1", "e_q coef=1: e^1 还是 e", type="e_q", coef="1")

# -------------------------------------------- 8. 退化解图像渲染
img(
    "k3:img-gauss-cu2",
    "gauss '<' 退化参数 cu_val=2 图像形态",
    type="gauss",
    m="0",
    n="0",
    a_val="0",
    b_val="5/6",
    c_val="0",
    au_val="5",
    bu_val="0",
    cu_val="2",
    u_val="1",
    coef="0",
    comparison="<",
    rational="1",
)
img(
    "k3:img-varpi-cu2",
    "varpi '<' 退化参数 cu_val=2 图像形态",
    type="varpi",
    m="0",
    n="0",
    a_val="0",
    b_val="1/4",
    c_val="0",
    au_val="21",
    bu_val="0",
    cu_val="2",
    u_val="1",
    coef="1",
    comparison="<",
    rational="15/4",
)
img(
    "k3:img-gauss-cu1",
    "gauss '>' 退化参数 cu_val=1 图像形态",
    type="gauss",
    m="0",
    n="0",
    a_val="0",
    b_val="1/8",
    c_val="0",
    au_val="3",
    bu_val="0",
    cu_val="1",
    u_val="1",
    coef="1",
    comparison=">",
    rational="1/4",
)


def site_request(p):
    if p["method"] == "GET":
        resp = SESSION.get(BASE + p["endpoint"], params=p["request"] or None, timeout=TIMEOUT)
    elif p["kind"] == "form":
        resp = SESSION.post(BASE + p["endpoint"], data=p["request"], timeout=TIMEOUT)
    else:
        raise ValueError(p["kind"])
    return resp.status_code, resp.text


def ours_request(client, p):
    if p["method"] == "GET":
        resp = client.get(p["endpoint"], query_string=p["request"] or None)
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
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--only", default=None)
    args = ap.parse_args()

    from attention_calculator.server import app

    client = app.test_client()

    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                with contextlib.suppress(json.JSONDecodeError):
                    done.add(json.loads(line)["id"])
    todo = [p for p in PROBES if p["id"] not in done]
    if args.only:
        todo = [p for p in todo if args.only in p["id"]]
    todo = todo[args.start :]
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
        rec.update(
            {
                "http_status": status,
                "raw": raw_body,
                "elapsed_ms": round(elapsed * 1000, 1),
                "ours_status": ostatus,
                "ours_body": obody,
                "byte_match": raw_body == obody,
                "match": status == ostatus and bodies_equal(raw_body, obody),
            }
        )
        with OUT.open("a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        flag = "" if rec["match"] else "  <<< MISMATCH"
        print(
            f"[{i}/{len(todo)}] {p['id']}: site={status} ours={ostatus} "
            f"({rec['elapsed_ms']:.0f}ms){flag}",
            flush=True,
        )


if __name__ == "__main__":
    main()
