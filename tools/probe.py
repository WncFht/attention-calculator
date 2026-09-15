"""Probe zhuyidao.net API: /calculate + /get_integral_image for given cases.

Usage: .venv/bin/python tools/probe.py cases.txt outdir
cases.txt lines: type power comparison rational   (comparison as > or <, url-encoded ok)

Rate-limited >=0.6s between requests. Writes one JSON per case: the /calculate
response merged with {"equation": ...} from /get_integral_image when successful.
"""

import json
import sys
import time
import urllib.parse

import requests

BASE = "https://zhuyidao.net"


def calculate(kind: str, power: str, comp: str, rational: str) -> dict:
    """POST /calculate and return parsed JSON."""
    r = requests.post(
        f"{BASE}/calculate",
        data={"type": kind, "power": power, "comparison": comp, "rational": rational},
        timeout=30,
    )
    return r.json()


def integral_image(resp: dict, kind: str, power: str, comp: str, rational: str) -> str:
    """GET /get_integral_image for a successful /calculate response."""
    p = resp["parameters"]
    query = {
        "m": p["m"], "n": p["n"],
        "a_val": p["a_val"], "b_val": p["b_val"], "c_val": p["c_val"],
        "au_val": p["au_val"], "bu_val": p["bu_val"], "cu_val": p["cu_val"],
        "u_val": p["u_val"],
        "comparison": comp, "rational": rational, "coef": power, "type": kind,
    }
    r = requests.get(f"{BASE}/get_integral_image?{urllib.parse.urlencode(query)}", timeout=30)
    return r.json()["equation"]


def main(cases_file: str, outdir: str) -> None:
    """Run all cases; print a compact summary line per case."""
    with open(cases_file) as fh:
        lines = list(fh)
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        kind, power, comp, rational = line.split()
        comp = urllib.parse.unquote(comp)
        time.sleep(0.7)
        try:
            resp = calculate(kind, power, comp, rational)
        except Exception as e:  # 网络边界
            print(f"{kind} {power} {comp} {rational}: HTTP-ERR {e}")
            continue
        out = dict(resp)
        out.update({"type": kind, "power": power, "comparison": comp, "rational": rational})
        if resp.get("success"):
            time.sleep(0.7)
            try:
                out["equation"] = integral_image(resp, kind, power, comp, rational)
            except Exception as e:
                out["equation_error"] = str(e)
        name = "{}_{}_{}_{}.json".format(
            kind, power.replace("/", "d"),
            "lt" if comp == "<" else "gt", rational.replace("/", "d"))
        with open(f"{outdir}/{name}", "w") as f:
            json.dump(out, f, ensure_ascii=False, indent=1)
        p = resp.get("parameters") or {}
        print(f"{kind:>14} {power:>4} {comp} {rational:>8}: ok={resp.get('success')} "
              f"m={p.get('m')} n={p.get('n')} u={p.get('u_val')} "
              f"a={p.get('a_val')} b={p.get('b_val')} c={p.get('c_val')} "
              f"err={resp.get('error', '')}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
