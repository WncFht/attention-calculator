# ruff: noqa: RUF001, RUF002, RUF003
"""随机差分模糊测试：同一请求打 zhuyidao.net 与本机克隆，比对响应找分歧。

/calculate 探针：type 从 29 个合法值里取（偶发非法串），power/rational 覆盖
小分数、真值邻域、大数（含 10^16 上限边界）、零分母、畸形/Unicode 写法，
comparison 偶发非法，偶发缺字段、query 代替 form、GET 方法。
/get_integral_image 探针：取 golden.jsonl 成功记录的 parameters，随机变异
（丢键、值改畸形/负数/巨数、coef/rational 用 \frac{}{} LaTeX 写法等）。

每条探针追加 bench/data/fuzz-probes.jsonl：
{endpoint, request, site_status, site_body, ours_status, ours_body, match, detail}
断点续跑：按 (endpoint, request) 去重跳过已采探针。

限速 ≥1.0s/站端请求，超时 90s；本机调用经线程池提交并在站端等待期间并行跑，
超时 LOCAL_TIMEOUT 记 ours_status=-2。站端 -1（传输失败）与 502/503/504 重试，
500 是站端业务语义不重试。

用法：.venv/bin/python bench/fuzz.py [--seed 1] [--calculate 750] [--image 150]
                                   [--max-minutes 45]
"""
import argparse
import concurrent.futures
import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import requests
from cases import true_value

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.kernels import TYPES
from attention_calculator.server import app

BASE = "https://zhuyidao.net"
DATA_DIR = Path(__file__).parent / "data"
OUT_PATH = DATA_DIR / "fuzz-probes.jsonl"
MIN_INTERVAL = 1.0
TIMEOUT = 90
RETRIES = 2          # 仅 -1 与 502/503/504 重试；500 是站端业务语义
LOCAL_TIMEOUT = 240  # 本机单次调用上限（秒），超时记 ours_status=-2

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-fuzz/0.1 (replication research)"
last_request = 0.0

CLIENT = app.test_client()
EXEC = concurrent.futures.ThreadPoolExecutor(max_workers=8)

# 数字字段畸形/边界写法池。Unicode 数字（٣、３、½）在 Python 里 \d 可匹配、
# int() 可解析，站端正则若按 ASCII 写就是分歧点；1_0 在 int() 下合法同理。
GARBAGE_NUM = [
    "", " ", "abc", "/", "//", "1/2/3", "1/", "/2", "1.5", ".5", "1.", "0.5",
    "1e3", "1E-2", "-1", "-3/2", "1/-2", "-1/2", "+1", "++2", "--1", "1_0",
    "0x10", "007", "01/02", " 3/4 ", "3 /4", "3/ 4", "\t1", "1\n", "1 2",
    "π", "pi", "e", "３", "٣", "٣/٤", "１２/３", "½", "²", "NaN", "nan",
    "inf", "∞", "-inf", "1,000", "1a", "a/1", "50%", "%", "0b1", "2..3",
    "3/2.5", "0/0", "1/0", "-0", "-0/1", "00", "000/000", "9" * 17,
    "1/" + "9" * 17, "9" * 17 + "/1", "10/" + "9" * 16, "1e-3", "\\frac{1}{2}",
    "1//2", "3.14", "x", "None", "null", "0x1F/2", "1 / 2", "  7  ",
]
GARBAGE_COMP = [
    "", " ", "=", "==", "!=", ">=", "<=", "=>", "=<", "≥", "≤", "＞", "＜",
    "><", ">>", "<<", "&gt;", "&lt;", "0", "1", "gt", "lt", "GT", "None",
    ">\t", " >", "> ", "< ", "≥ ", "\n>",
]
GARBAGE_TYPE = [
    "", " ", "PI", "Pi", "pI", " pi", "pi ", " pi ", "π", "pi_n ", "e_q ",
    "sin", "cos", "tan", "ln", "lnq", "ln_q2", "ln_q_squares", "sin_pi",
    "foo", "bar", "0", "1", "null", "None", "true", "pi,e", "pi;e", "Gamma",
    "GAMMA", "zeta", "zeta_3", "arctan", "arcoth", "sinq", "-pi", "pi\n",
    "𝜋", "e_pi ", "golden2", "gauss ",
]
FIELDS = ("type", "power", "comparison", "rational")

# 各类型定义域边界附近的 power 取值：验证站端开闭区间判定与我们是否一致
EDGE_POWER = {
    "ln_q": ["1", "1/1", "2/2", "999/1000", "1001/1000", "0", "1/2"],
    "ln_q_square": ["1", "1/1", "2/2", "999/1000", "1001/1000", "0", "1/2"],
    "sin_q": ["0", "3", "4", "22/7", "355/113", "314/100", "157/50"],
    "cos_q": ["0", "2", "11/7", "157/100", "16/10", "1571/1000"],
    "tan_q": ["0", "2", "11/7", "157/100", "16/10", "1571/1000"],
    "cot_q": ["0", "2", "11/7", "157/100", "16/10", "1571/1000"],
    "sin_q_degree": ["0", "45", "89", "90", "90/1", "899/10", "901/10"],
    "cos_q_degree": ["0", "45", "89", "90", "90/1", "899/10", "901/10"],
    "sin_pi_q": ["0", "1", "2", "1/2", "2/4", "49/100", "51/100", "3/6"],
    "cos_pi_q": ["0", "1", "2", "1/2", "2/4", "49/100", "51/100", "3/6"],
    "artanh_q": ["0", "1", "2", "1/2", "9/10", "10/10", "99/100"],
    "arcoth_q": ["0", "1", "1/2", "9/10", "10/10", "11/10", "2/1"],
}


def log(msg):
    """带时间戳的进度输出。"""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def site_request(method, path, **kwargs):
    """限速 ≥1s 的站端请求，返回 (http_status, body_text, elapsed_ms)。

    -1 传输失败、502/503/504 网关错误重试 ≤RETRIES 次；500 不重试（业务语义）。
    """
    global last_request
    status, body, elapsed = -1, "", 0.0
    for attempt in range(RETRIES + 1):
        wait = MIN_INTERVAL - (time.monotonic() - last_request)
        if wait > 0:
            time.sleep(wait)
        t0 = time.monotonic()
        try:
            resp = SESSION.request(method, BASE + path, timeout=TIMEOUT, **kwargs)
            status, body = resp.status_code, resp.text
        except (OSError, requests.RequestException) as exc:
            status, body = -1, f"REQUEST_FAILED: {exc}"
        elapsed = (time.monotonic() - t0) * 1000
        last_request = time.monotonic()
        if not (status == -1 or status in (502, 503, 504)) or attempt == RETRIES:
            break
        time.sleep(2 ** (attempt + 1))
    return status, body, elapsed


def rand_num(rng):
    """随机数字字段文本：小分数、中位数、10^16 边界巨数、零、畸形。"""
    r = rng.random()
    if r < 0.36:
        n, d = rng.randint(0, 12), rng.randint(1, 12)
        return f"{n}/{d}" if rng.random() < 0.65 else str(n)
    if r < 0.56:
        n, d = rng.randint(0, 10**4), rng.randint(1, 10**4)
        return f"{n}/{d}" if rng.random() < 0.7 else str(n)
    if r < 0.68:
        n = rng.choice([10**15, 10**16 - 1, 10**16, 10**16 + 7,
                        rng.randint(10**16, 10**20), rng.randint(0, 10**17)])
        d = rng.choice([1, 3, 10**16 - 1, 10**16, rng.randint(1, 10**12)])
        return f"{n}/{d}"
    if r < 0.74:
        return rng.choice(["0", "0/1", "0/7"])
    return rng.choice(GARBAGE_NUM)


def near_bound(rng, kind, power):
    """围绕常数真值的随机有理界：floor(v·d)±k / d，d 取杂散分母。"""
    try:
        v = true_value(kind, Fraction(power))
        v = float(v)
        if not (v == v) or v in (float("inf"), float("-inf")) or v <= 0:
            raise ValueError
    except Exception:
        return None
    d = rng.choice([rng.randint(1, 20), rng.randint(1, 500), rng.randint(1, 5000)])
    n = int(v * d) + rng.choice([-1, 0, 0, 0, 1, 1, 2])
    if n < 0:
        n = 0
    return str(n) if rng.random() < 0.15 else f"{n}/{d}"


def gen_calculate(rng):
    """生成一条 /calculate 探针 {method, via, form}；form 可能缺字段。"""
    kind = rng.choice(TYPES) if rng.random() < 0.85 else rng.choice(GARBAGE_TYPE)
    # power：20% 定义域边界、55% 随机数、25% 畸形池
    r = rng.random()
    if r < 0.20 and kind in EDGE_POWER:
        power = rng.choice(EDGE_POWER[kind])
    elif r < 0.75:
        power = rand_num(rng)
    else:
        power = rng.choice(GARBAGE_NUM)
    comp = rng.choice("><") if rng.random() < 0.88 else rng.choice(GARBAGE_COMP)
    # rational：40% 真值邻域、35% 随机数、25% 畸形池
    r = rng.random()
    bound = None
    if r < 0.40:
        bound = near_bound(rng, kind, power)
    if bound is None:
        bound = rand_num(rng) if r < 0.75 else rng.choice(GARBAGE_NUM)
    form = {"type": kind, "power": power, "comparison": comp, "rational": bound}
    # 偶发缺字段（每字段 4%），偶发全空
    for k in FIELDS:
        if rng.random() < 0.04:
            form.pop(k)
    if rng.random() < 0.01:
        form = {}
    via = "query" if rng.random() < 0.03 else "form"
    method = "GET" if rng.random() < 0.015 else "POST"
    return {"method": method, "via": via, "form": form}


def load_golden_pool():
    """golden 成功记录的 image 查询参数池（与 harvest.image_params 同构）。"""
    pool = []
    for line in (DATA_DIR / "golden.jsonl").open():
        if not line.strip():
            continue
        rec = json.loads(line)
        if not rec.get("success"):
            continue
        params = {
            k: json.dumps(v, ensure_ascii=False)
            if isinstance(v, (dict, list)) else str(v)
            for k, v in rec["parameters"].items()
        }
        params.update({"type": rec["type"], "coef": rec["power"],
                       "comparison": rec["comparison"], "rational": rec["rational"]})
        pool.append(params)
    return pool


def mutate_image(rng, params):
    """对 image 查询参数施加 0-3 个变异，返回 (新 dict, 变异描述 list)。"""
    p = dict(params)
    applied = []
    if rng.random() < 0.05:
        # 全随机参数：所有键都给畸形值
        for k in list(p):
            p[k] = rng.choice(GARBAGE_NUM)
        return p, ["scramble"]
    for _ in range(rng.choices([0, 1, 2, 3], weights=[15, 55, 20, 10])[0]):
        op = rng.choice(["drop", "corrupt", "neg", "huge", "zero_u",
                         "latex_coef", "latex_rat", "badtype", "badcomp",
                         "swap", "empty"])
        keys = list(p)
        if op == "drop" and keys:
            k = rng.choice(keys)
            del p[k]
            applied.append(f"drop:{k}")
        elif op == "corrupt":
            k = rng.choice(list(p))
            p[k] = rng.choice(GARBAGE_NUM)
            applied.append(f"corrupt:{k}={p[k]!r}")
        elif op == "neg":
            k = rng.choice(list(p))
            v = p[k]
            p[k] = v[1:] if v.startswith("-") else f"-{v}"
            applied.append(f"neg:{k}={p[k]!r}")
        elif op == "huge":
            k = rng.choice(list(p))
            p[k] = str(rng.choice([-1, 1]) * 10 ** rng.randint(16, 30))
            applied.append(f"huge:{k}")
        elif op == "zero_u":
            p["u_val"] = "0"
            applied.append("zero_u")
        elif op == "latex_coef":
            p["coef"] = f"\\{rng.choice(['frac', 'dfrac'])}" \
                        f"{{{rng.randint(-9, 9)}}}{{{rng.randint(-9, 9)}}}"
            applied.append(f"latex_coef={p['coef']}")
        elif op == "latex_rat":
            p["rational"] = f"\\{rng.choice(['frac', 'dfrac'])}" \
                            f"{{{rng.randint(-9, 9)}}}{{{rng.randint(-9, 9)}}}"
            applied.append(f"latex_rat={p['rational']}")
        elif op == "badtype":
            p["type"] = rng.choice(GARBAGE_TYPE)
            applied.append(f"badtype={p['type']!r}")
        elif op == "badcomp":
            p["comparison"] = rng.choice(GARBAGE_COMP)
            applied.append(f"badcomp={p['comparison']!r}")
        elif op == "swap" and len(keys) >= 2:
            a, b = rng.sample(keys, 2)
            p[a], p[b] = p[b], p[a]
            applied.append(f"swap:{a}<->{b}")
        elif op == "empty":
            k = rng.choice(list(p))
            p[k] = ""
            applied.append(f"empty:{k}")
    return p, applied


def ours_calculate(req):
    """本机 /calculate 调用（线程内执行，返回 status, body）。"""
    if req["method"] == "GET":
        return CLIENT.get("/calculate", query_string=req["form"])
    if req["via"] == "query":
        return CLIENT.post("/calculate", query_string=req["form"])
    return CLIENT.post("/calculate", data=req["form"])


def ours_image(params):
    """本机 /get_integral_image 调用。"""
    return CLIENT.get("/get_integral_image", query_string=params)


def collect(future):
    """取本机调用结果；超时或异常记 -2/-3。"""
    try:
        resp = future.result(timeout=LOCAL_TIMEOUT)
        return resp.status_code, resp.get_data(as_text=True)
    except concurrent.futures.TimeoutError:
        return -2, "OURS_TIMEOUT"
    except Exception as exc:
        return -3, f"OURS_CRASH: {type(exc).__name__}: {exc}"


def parse_body(body):
    """响应体解析为 JSON；失败返回 None。"""
    try:
        return json.loads(body)
    except (json.JSONDecodeError, TypeError):
        return None


def compare(site_status, site_body, ours_status, ours_body):
    """状态码 + 全 JSON 字典比对，返回 (match, detail)。"""
    if site_status != ours_status:
        return False, f"status site={site_status} ours={ours_status}"
    sp, op = parse_body(site_body), parse_body(ours_body)
    if sp is None or op is None:
        if site_body == ours_body:
            return True, "non-json identical"
        return False, f"non-json site={site_body[:80]!r} ours={ours_body[:80]!r}"
    if sp == op:
        return True, "exact"
    diff = {k: [sp.get(k), op.get(k)] for k in set(sp) | set(op)
            if sp.get(k) != op.get(k)}
    return False, f"json diff {json.dumps(diff, ensure_ascii=False)[:300]}"


def seen_keys():
    """已采探针的 (endpoint, request) 键集，用于续跑去重。"""
    keys = set()
    if OUT_PATH.exists():
        for line in OUT_PATH.open():
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            keys.add((rec.get("endpoint"),
                      json.dumps(rec.get("request"), sort_keys=True)))
    return keys


def emit(rec):
    """追加一条探针记录（每次重新 open，防外部文件替换丢写）。"""
    with OUT_PATH.open("a") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def run_probe(endpoint, request, ours_fn, site_args):
    """单条探针：本机调用与站端请求并行（本机借站端限速等待完成），比对落盘。"""
    fut = EXEC.submit(ours_fn)
    status, body, ms = site_request(**site_args)
    ours_status, ours_body = collect(fut)
    match, detail = compare(status, body, ours_status, ours_body)
    rec = {"endpoint": endpoint, "request": request,
           "site_status": status, "site_body": body,
           "ours_status": ours_status, "ours_body": ours_body,
           "match": match, "detail": detail, "site_ms": round(ms, 1)}
    emit(rec)
    return match


def calc_site_args(req):
    """把 {method, via, form} 转成 site_request 的调用参数。"""
    kwargs = {"method": req["method"], "path": "/calculate"}
    kwargs["params" if req["via"] == "query" or req["method"] == "GET"
           else "data"] = req["form"]
    return kwargs


def main():
    """入口：交替跑 /calculate 与 /get_integral_image 探针直到目标数或时限。"""
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--calculate", type=int, default=750)
    ap.add_argument("--image", type=int, default=150)
    ap.add_argument("--max-minutes", type=float, default=45)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    pool = load_golden_pool()
    done = seen_keys()
    deadline = time.time() + args.max_minutes * 60
    stats = {"calc": 0, "img": 0, "div": 0}

    # 预生成全部探针（确定性：同 seed 同序列）；按 (endpoint, request) 去重续跑
    calc_reqs, img_reqs = [], []
    while len(calc_reqs) < args.calculate:
        req = gen_calculate(rng)
        key = ("/calculate", json.dumps(req, sort_keys=True))
        if key not in done:
            done.add(key)
            calc_reqs.append(req)
    while len(img_reqs) < args.image:
        params, muts = mutate_image(rng, rng.choice(pool))
        req = {"params": params, "mutations": muts}
        key = ("/get_integral_image", json.dumps(req, sort_keys=True))
        if key not in done:
            done.add(key)
            img_reqs.append(req)
    log(f"plan: {len(calc_reqs)} calculate + {len(img_reqs)} image, "
        f"deadline {args.max_minutes}min")

    img_i = 0
    for i, req in enumerate(calc_reqs, 1):
        if time.time() > deadline:
            log("time budget exhausted")
            break
        ok = run_probe(
            "/calculate", req,
            lambda req=req: ours_calculate(req),
            calc_site_args(req))
        stats["calc"] += 1
        stats["div"] += not ok
        # 每 5 条 calculate 插 1 条 image 探针
        if i % 5 == 0 and img_i < len(img_reqs) and time.time() <= deadline:
            ireq = img_reqs[img_i]
            img_i += 1
            ok = run_probe(
                "/get_integral_image", ireq,
                lambda ireq=ireq: ours_image(ireq["params"]),
                {"method": "GET", "path": "/get_integral_image",
                 "params": ireq["params"]})
            stats["img"] += 1
            stats["div"] += not ok
        if i % 25 == 0 or i == len(calc_reqs):
            log(f"progress: calc={stats['calc']} img={stats['img']} "
                f"divergences={stats['div']}")

    # 剩余 image 探针补跑
    while img_i < len(img_reqs) and time.time() <= deadline:
        ireq = img_reqs[img_i]
        img_i += 1
        ok = run_probe(
            "/get_integral_image", ireq,
            lambda ireq=ireq: ours_image(ireq["params"]),
            {"method": "GET", "path": "/get_integral_image",
             "params": ireq["params"]})
        stats["img"] += 1
        stats["div"] += not ok
    log(f"done: calc={stats['calc']} img={stats['img']} divergences={stats['div']}")


if __name__ == "__main__":
    sys.exit(main())
