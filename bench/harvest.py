# ruff: noqa: RUF002, RUF003
"""礼貌采集 zhuyidao.net，产出 bench/data/golden.jsonl 与 bench/data/combo.jsonl。

golden 每条 case：POST /calculate → 成功则 GET /get_integral_image 取渲染式；
原始响应体逐字存入 raw_calculate / raw_image，解析字段平铺进记录。
combo 每条 problem：POST /decompose_inequality。

限速 ≥0.7s/请求，超时 90s，网络错误与 5xx 指数退避重试 ≤3 次。
断点续跑：已写入输出文件的 (type,power,comparison,rational) / problem 跳过。

用法：python bench/harvest.py [--only golden|combo] [--limit N]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import requests
from cases import combo_problems, generate_cases

BASE = "https://zhuyidao.net"
DATA_DIR = Path(__file__).parent / "data"
MIN_INTERVAL = 0.7
TIMEOUT = 90
RETRIES = 3

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (replication research)"

last_request = 0.0


def log(msg):
    """带时间戳的进度输出（nohup → harvest.log）。"""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def request(method, path, **kwargs):
    """限速 + 重试的单次请求，返回 (http_status, body_text, elapsed_ms)。

    body_text 原样保留响应体；网络层失败时 body 为 'REQUEST_FAILED: ...'。
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
        except OSError as exc:  # RequestException 是 OSError 子类；certifi 缺失等抛裸 OSError
            status, body = -1, f"REQUEST_FAILED: {exc}"
        elapsed = (time.monotonic() - t0) * 1000
        last_request = time.monotonic()
        if status < 500 or attempt == RETRIES:
            break
        time.sleep(2 ** (attempt + 1))
    return status, body, elapsed


def image_params(case, parameters):
    """/get_integral_image 查询串：parameters 全部字段 + type/coef/comparison/rational。"""
    params = {
        k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else str(v)
        for k, v in parameters.items()
    }
    params.update({
        "type": case["type"], "coef": case["power"],
        "comparison": case["comparison"], "rational": case["rational"],
    })
    return params


def fetch_case(case):
    """单 case：/calculate 记录全量字段；success 时追加渲染 equation。"""
    rec = dict(case)
    status, body, elapsed = request("POST", "/calculate", data=case)
    rec["elapsed_ms"] = round(elapsed, 1)
    rec["http_status"] = status
    rec["raw_calculate"] = body
    try:
        payload = json.loads(body)
    except json.JSONDecodeError:
        rec["success"] = False
        rec["error"] = f"non-json response (http {status})"
        return rec
    rec.update(payload)
    rec["type"] = case["type"]
    rec["success"] = bool(payload.get("success"))  # 错误响应只有 error 字段
    if not rec["success"]:
        return rec
    status2, body2, elapsed2 = request(
        "GET", "/get_integral_image", params=image_params(case, payload["parameters"]))
    rec["image_elapsed_ms"] = round(elapsed2, 1)
    rec["image_http_status"] = status2
    rec["raw_image"] = body2
    try:
        rec["equation"] = json.loads(body2).get("equation")
    except json.JSONDecodeError:
        rec["equation"] = None
        rec["image_error"] = f"non-json response (http {status2})"
    return rec


def load_keys(path, key_fields):
    """读已有 jsonl 的去重键集合；损坏行跳过（上次中断可能留半行）。"""
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            done.add(tuple(rec.get(k) for k in key_fields))
    return done


def harvest_golden(limit=None):
    """采集 /calculate 全量 case，追加写 golden.jsonl。"""
    out_path = DATA_DIR / "golden.jsonl"
    done = load_keys(out_path, ("type", "power", "comparison", "rational"))
    cases = [c for c in generate_cases()
             if (c["type"], c["power"], c["comparison"], c["rational"]) not in done]
    if limit:
        cases = cases[:limit]
    log(f"golden: {len(done)} done, {len(cases)} to fetch")
    for i, case in enumerate(cases, 1):
        rec = fetch_case(case)
        # 每条重新 open：本机 rsync 可能整体替换文件，旧句柄会写进已删 inode
        with out_path.open("a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        if i % 25 == 0 or i == len(cases):
            log(f"golden {i}/{len(cases)} ok={rec['success']}")


def harvest_combo():
    """采集 /decompose_inequality 问题，追加写 combo.jsonl。"""
    out_path = DATA_DIR / "combo.jsonl"
    done = load_keys(out_path, ("problem",))
    problems = [p for p in combo_problems() if (p,) not in done]
    log(f"combo: {len(done)} done, {len(problems)} to fetch")
    for i, problem in enumerate(problems, 1):
        rec = {"problem": problem}
        status, body, elapsed = request(
            "POST", "/decompose_inequality", data={"problem": problem})
        rec["elapsed_ms"] = round(elapsed, 1)
        rec["http_status"] = status
        rec["raw"] = body
        try:
            rec.update(json.loads(body))
        except json.JSONDecodeError:
            rec["error"] = f"non-json response (http {status})"
        rec["success"] = bool(rec.get("success"))  # 错误响应只有 error 字段
        rec["problem"] = problem
        with out_path.open("a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        if i % 10 == 0 or i == len(problems):
            log(f"combo {i}/{len(problems)} ok={rec['success']}")


def main():
    """入口：默认先 golden 后 combo，--only/--limit 可局部跑。"""
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["golden", "combo"], default=None)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    DATA_DIR.mkdir(exist_ok=True)
    if args.only in (None, "golden"):
        harvest_golden(args.limit)
    if args.only in (None, "combo"):
        harvest_combo()


if __name__ == "__main__":
    sys.exit(main())
