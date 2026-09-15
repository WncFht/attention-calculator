"""第三波补探：方向预检的作用域与次序。

已确认 e_q/arctan_q/arccot_q/hyperbolic/gamma 在核前先按 float64 常量值判定
不等号方向（假->404 方向反了，真/相等->进核）。要钉死：
- 预检 vs 定义域检查的先后（domain-first 还是 precheck-first）；
- pi_n / e / e_pi 是否也走这道预检；
- e_q 大指数的 math.exp 溢出。
"""

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
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (edge probes w3b)"

PROBES = [
    dict(
        id="k3b:pin0-gt2",
        note="pi_n 0>2: 共享预检->404方向反了 / 仅核内->500",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "pi_n", "power": "0", "comparison": ">", "rational": "2"},
    ),
    dict(
        id="k3b:pin11-gt1",
        note="pi_n 11>1: pi^11 大真值->过检->spec崩500?",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "pi_n", "power": "11", "comparison": ">", "rational": "1"},
    ),
    dict(
        id="k3b:lnq1-gt2",
        note="ln_q 1>2: 域错 vs 方向反了 谁先",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "ln_q", "power": "1", "comparison": ">", "rational": "2"},
    ),
    dict(
        id="k3b:tanq0-gt1",
        note="tan_q 0>1: 域错 vs 方向反了",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "tan_q", "power": "0", "comparison": ">", "rational": "1"},
    ),
    dict(
        id="k3b:sinq0-gt1",
        note="sin_q 0>1: 域错 vs 方向反了",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "sin_q", "power": "0", "comparison": ">", "rational": "1"},
    ),
    dict(
        id="k3b:atanh0-gt0",
        note="artanh_q 0>0: 域错(整数) vs 方向反了",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "artanh_q", "power": "0", "comparison": ">", "rational": "0"},
    ),
    dict(
        id="k3b:eq800-gt1",
        note="e_q 800>1: math.exp 溢出->500?",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "e_q", "power": "800", "comparison": ">", "rational": "1"},
    ),
    dict(
        id="k3b:e0-gt2",
        note="e 0>2: 假命题(常数=0*e) 方向反了 or 核路径",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "e", "power": "0", "comparison": ">", "rational": "2"},
    ),
    dict(
        id="k3b:epi0-gt2",
        note="e_pi 0>2: 同上",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "e_pi", "power": "0", "comparison": ">", "rational": "2"},
    ),
    dict(
        id="k3b:pi0-gt2",
        note="pi 0>2: 同上",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "pi", "power": "0", "comparison": ">", "rational": "2"},
    ),
    dict(
        id="k3b:sinh0-lt0",
        note="sinh 0<0: 相等->进核崩500?",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "sinh_q", "power": "0", "comparison": "<", "rational": "0"},
    ),
    dict(
        id="k3b:varpi0-gt2",
        note="varpi 0>2: 假->方向反了?",
        endpoint="/calculate",
        method="POST",
        kind="form",
        request={"type": "varpi", "power": "0", "comparison": ">", "rational": "2"},
    ),
]


def main():
    from attention_calculator.server import app

    client = app.test_client()
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["id"])
    todo = [p for p in PROBES if p["id"] not in done]
    print(f"{len(todo)} to probe", flush=True)
    last = 0.0
    for i, p in enumerate(todo, 1):
        wait = MIN_INTERVAL - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
        try:
            resp = SESSION.post(BASE + p["endpoint"], data=p["request"], timeout=TIMEOUT)
            status, raw = resp.status_code, resp.text
        except (OSError, requests.RequestException) as exc:
            status, raw = -1, f"REQUEST_FAILED: {exc}"
        last = time.monotonic()
        try:
            r2 = client.post(p["endpoint"], data=p["request"])
            ostatus, obody = r2.status_code, r2.get_data(as_text=True)
        except Exception as exc:
            ostatus, obody = -2, f"LOCAL_CRASH: {exc}"
        rec = dict(p)
        rec.update(
            {
                "http_status": status,
                "raw": raw,
                "ours_status": ostatus,
                "ours_body": obody,
                "byte_match": raw == obody,
                "match": status == ostatus and raw == obody,
            }
        )
        with OUT.open("a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        flag = "" if rec["match"] else "  <<< MISMATCH"
        print(f"[{i}/{len(todo)}] {p['id']}: site={status} ours={ostatus}{flag}", flush=True)
        print(f"     site: {raw[:140]}", flush=True)


if __name__ == "__main__":
    import json

    main()
