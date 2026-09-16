"""离线重放无判官覆盖的采集集：fidelity-probes.jsonl + probes.jsonl。

fidelity 谱系（edge schema）：按 http_status/raw 对比，byte_match 记录逐字节一致。
probes 谱系（旧 schema）：按 http/response（已解析 JSON）对比。
两者均就地刷新 ours_status/ours_body/match 字段，与 replay_edge.py 同约定。
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

DATA = Path(__file__).parent / "data"


def ours_request(client, endpoint, form):
    if endpoint == "/get_integral_image":
        resp = client.get(endpoint, query_string=form)
    else:
        resp = client.post(endpoint, data=form)
    return resp.status_code, resp.get_data(as_text=True)


def bodies_equal(site_raw, ours_raw):
    try:
        s, o = json.loads(site_raw), json.loads(ours_raw)
    except (json.JSONDecodeError, TypeError):
        return site_raw.strip() == ours_raw.strip()
    return s == o


def err_of(body):
    try:
        return repr(json.loads(body).get("error", body))[:120]
    except (json.JSONDecodeError, AttributeError):
        return repr(body[:120])


def replay_edge_schema(client, path):
    recs = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    n_match = n_mismatch = n_skip = n_byte = 0
    for r in recs:
        if r["http_status"] == -1:
            n_skip += 1
            continue
        try:
            ostatus, obody = ours_request(client, r["endpoint"], r["request"])
        except Exception as exc:
            ostatus, obody = -2, f"LOCAL_CRASH: {exc}"
        r["ours_status"] = ostatus
        r["ours_body"] = obody
        r["byte_match"] = r["raw"] == obody
        r["match"] = r["http_status"] == ostatus and bodies_equal(r["raw"], obody)
        n_byte += r["byte_match"]
        if r["match"]:
            n_match += 1
        else:
            n_mismatch += 1
            print(
                f"MISMATCH {r['id']}: site={r['http_status']} {err_of(r['raw'])} "
                f"ours={ostatus} {err_of(obody)} note={r.get('note', '')}",
                flush=True,
            )
    with path.open("w") as f:
        for r in recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(
        f"{path.name}: total={len(recs)} match={n_match} mismatch={n_mismatch} "
        f"skipped={n_skip} byte_exact={n_byte}"
    )


def replay_probes_schema(client, path):
    recs = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    n_match = n_mismatch = n_skip = 0
    # 同一 (endpoint,form) 既有 500 又有非 500 采集 → 该 500 是站端瞬时故障
    # （grid:*:img 七连 500 在 77s 后同参数全 200 实测），跳过时点名
    forms_with_200 = {
        json.dumps(r["form"], sort_keys=True) for r in recs if r["http"] not in (-1, 500)
    }
    for r in recs:
        if r["http"] == -1:
            n_skip += 1
            continue
        if r["http"] == 500 and json.dumps(r["form"], sort_keys=True) in forms_with_200:
            n_skip += 1
            print(f"SKIP-TRANSIENT {r['tag']}: 同参数另有非 500 采集", flush=True)
            continue
        try:
            ostatus, obody = ours_request(client, r["endpoint"], r["form"])
        except Exception as exc:
            ostatus, obody = -2, f"LOCAL_CRASH: {exc}"
        r["ours_status"] = ostatus
        r["ours_body"] = obody
        r["match"] = r["http"] == ostatus and bodies_equal(
            json.dumps(r["response"], ensure_ascii=False), obody
        )
        if r["match"]:
            n_match += 1
        else:
            n_mismatch += 1
            site_err = (
                repr(r["response"].get("error", r["response"]))[:120]
                if isinstance(r["response"], dict)
                else repr(r["response"])[:120]
            )
            print(
                f"MISMATCH {r['tag']}: site={r['http']} {site_err} ours={ostatus} {err_of(obody)}",
                flush=True,
            )
    with path.open("w") as f:
        for r in recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"{path.name}: total={len(recs)} match={n_match} mismatch={n_mismatch} skipped={n_skip}")


def main():
    from attention_calculator.server import app

    client = app.test_client()
    replay_edge_schema(client, DATA / "fidelity-probes.jsonl")
    replay_probes_schema(client, DATA / "probes.jsonl")


if __name__ == "__main__":
    main()
