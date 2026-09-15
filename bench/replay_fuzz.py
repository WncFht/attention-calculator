"""离线重放 fuzz-probes.jsonl：只打本方 test_client，不触网。

fuzz.py 采集的 request 形态：
- /calculate ``via=form``  -> POST 表单；``via=query`` -> POST 把字段放 query string
  （站端只读 form，此类一律 400 右侧格式无效）
- /get_integral_image     -> GET，request.params 即 query dict

就地刷新 ours_status/ours_body/match/detail，输出剩余分歧统计。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.server import app

OUT = Path(__file__).parent / "data" / "fuzz-probes.jsonl"


def ours_request(client, rec: dict) -> tuple[int, str]:
    """Replay one stored fuzz request against the current app."""
    req = rec["request"]
    if rec["endpoint"] == "/get_integral_image":
        resp = client.get("/get_integral_image", query_string=req["params"])
    elif req.get("method") == "GET":
        resp = client.get("/calculate", query_string=req["form"])
    elif req.get("via") == "query":
        resp = client.post("/calculate", query_string=req["form"])
    else:
        resp = client.post("/calculate", data=req["form"])
    return resp.status_code, resp.get_data(as_text=True)


def bodies_equal(site_raw: str, ours_raw: str) -> bool:
    """JSON-compare when both parse (byte-level is covered by golden body_match)."""
    try:
        return json.loads(site_raw) == json.loads(ours_raw)
    except (json.JSONDecodeError, TypeError):
        return site_raw.strip() == ours_raw.strip()


def main() -> None:
    """Replay every stored probe and rewrite ours_* fields in place."""
    recs = [json.loads(line) for line in OUT.open() if line.strip()]
    client = app.test_client()
    n_match = 0
    mismatches = []
    for rec in recs:
        status, body = ours_request(client, rec)
        rec["ours_status"], rec["ours_body"] = status, body
        rec["match"] = status == rec["site_status"] and bodies_equal(rec["site_body"], body)
        rec["detail"] = "exact" if rec["match"] else "replayed-diff"
        if rec["match"]:
            n_match += 1
        else:
            mismatches.append(rec)
    with OUT.open("w") as f:
        for rec in recs:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    print(f"replayed {len(recs)}: match={n_match} mismatch={len(mismatches)}")
    for r in mismatches:
        print(f"MISMATCH {r['endpoint']} {json.dumps(r['request'], ensure_ascii=False)[:140]}")
        print(f"  site={r['site_status']} {r['site_body'][:120]}")
        print(f"  ours={r['ours_status']} {r['ours_body'][:120]}")


if __name__ == "__main__":
    main()
