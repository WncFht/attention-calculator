# ruff: noqa: RUF002, RUF003
"""离线重放 edge-probes.jsonl：只打本方 test_client，不触网。

两波探针的 ours_* 字段分别采集自改动前/后的 server.py；
工作树被并行修改后，此脚本把每条已存 request 重放到当前代码，
就地刷新 ours_status/ours_body/byte_match/match。

kind=form-dup 的记录 request 是 list-of-pairs，用原始 urlencoded body 重放。
"""
import json
import sys
from pathlib import Path
from urllib.parse import quote

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

OUT = Path(__file__).parent / "data" / "edge-probes.jsonl"


def ours_request(client, p):
    method, kind, req = p["method"], p["kind"], p["request"]
    if kind == "form-dup":
        body = "&".join(f"{k}={quote(str(v), safe='')}" for k, v in req)
        resp = client.post(
            p["endpoint"], data=body,
            content_type="application/x-www-form-urlencoded")
    elif method == "GET":
        resp = client.get(p["endpoint"], query_string=req or None)
    elif method == "PUT":
        resp = client.put(p["endpoint"], data=req)
    elif method == "DELETE":
        resp = client.delete(p["endpoint"], query_string=req)
    elif method == "OPTIONS":
        resp = client.options(p["endpoint"])
    elif kind == "query":
        resp = client.post(p["endpoint"], query_string=req)
    elif kind == "json":
        resp = client.post(p["endpoint"], json=req)
    else:
        resp = client.post(p["endpoint"], data=req)
    return resp.status_code, resp.get_data(as_text=True)


def bodies_equal(site_raw, ours_raw):
    try:
        s, o = json.loads(site_raw), json.loads(ours_raw)
    except (json.JSONDecodeError, TypeError):
        return site_raw.strip() == ours_raw.strip()
    return s == o


def main():
    from attention_calculator.server import app
    client = app.test_client()

    recs = [json.loads(line) for line in OUT.read_text().splitlines()
            if line.strip()]
    n_match = n_mismatch = n_skip = 0
    for r in recs:
        if r["http_status"] in (-1,):
            n_skip += 1
            continue
        try:
            ostatus, obody = ours_request(client, r)
        except Exception as exc:
            ostatus, obody = -2, f"LOCAL_CRASH: {type(exc).__name__}: {exc}"
        r["ours_status"] = ostatus
        r["ours_body"] = obody
        r["byte_match"] = r["raw"] == obody
        r["match"] = (r["http_status"] == ostatus
                      and bodies_equal(r["raw"], obody))
        if r["match"]:
            n_match += 1
        else:
            n_mismatch += 1
            print(f"MISMATCH {r['id']}: site={r['http_status']} "
                  f"ours={ostatus}", flush=True)

    with OUT.open("w") as f:
        for r in recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"replayed {len(recs)}: match={n_match} mismatch={n_mismatch} "
          f"skipped={n_skip}")


if __name__ == "__main__":
    main()
