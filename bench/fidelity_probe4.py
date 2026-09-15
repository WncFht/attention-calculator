"""Battery 4: discriminate exact-vs-float compare per type + squeeze windows.

Model variants under test per type:
  EXACT: precheck compares Fraction(r) vs C' (a Fraction, usually
         Fraction(float64(C))). '>' WD iff r>C'; '<' WD iff r<C'.
  FLOAT: precheck compares float64(r) vs float64(C'). '>' WD iff rf>Vf.
         Float-equal bounds always scan.

Discriminating probes:
  - '>' on a bound float-equal to Cf but strictly above Cf:
      EXACT(C'=Cf) -> WD;  FLOAT(Vf=Cf) -> scan.
  - '<' on Cf-1ulp: FLOAT/exact-Cf -> WD;  absent '<' check -> scan.

Also: lnsq dyadic j=-1/+1 (corrected; battery-3 had mp.dps bug -> den 2^54),
and ulp probes for e_pi/gauss/varpi/gamma/catalan '<' side.
"""

import contextlib
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from mpmath import mp

mp.dps = 120
from attention_calculator.integrand import constant_mpf  # noqa: E402

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "fidelity-probes.jsonl"
MIN_INTERVAL = 1.05
TIMEOUT = 120

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (fidelity probes v4)"

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


def dyadic(kind, power, j):
    f = float(constant_mpf(kind, Fraction(power)))
    for _ in range(abs(j)):
        f = math.nextafter(f, math.inf if j > 0 else -math.inf)
    return Fraction(f)


def ulp(tag, kind, power, comp, j):
    r = dyadic(kind, power, j)
    calc(
        f"fid4:{tag}-{comp}{j:+d}",
        f"{kind}({power}){comp} Cf{j:+d}ulp",
        type=kind,
        power=power,
        comparison=comp,
        rational=f"{r.numerator}/{r.denominator}",
    )


# ---- zeta3: '>' on float-equal bounds (squeeze C'_gt / test FLOAT) --------
# FLOAT model -> every floateq bound scans (NS). EXACT C' -> WD iff r>C'.
for tag, r, note in [
    ("a", "461424925/383862797", "r=C+5.76e-19 floateq"),
    ("b", "338896541/281930531", "r=C+1.91e-17 floateq"),
    ("c", "400160733/332896664", "r=C+8.40e-18 floateq"),
    ("d", "61264192/50966133", "r=C-5.05e-17 floateq"),
    ("e", "522689117/434828930", "r=C-5.41e-18 floateq"),
    ("f", "3752663592/3121868509", "r=C-2.58e-19 floateq"),
]:
    calc(
        f"fid4:zeta-gt-fe{tag}",
        f"zeta3> {note}",
        type="zeta3",
        power="1",
        comparison=">",
        rational=r,
    )
ulp("zeta", "zeta3", "1", "<", -1)  # '<' check exists?
ulp("zeta", "zeta3", "1", ">", -1)  # r=Cf-ulp<C: '>' scan -> NS either way

# ---- e_q(3): '>' float-equal above Cf : EXACT->WD, FLOAT->NS(scan ~1.5s) ----
for tag, r, note in [
    ("a", "2424466693/120707089", "r=Cf+3.5e-17 floateq"),
    ("b", "710115019/35354545", "r=Cf+1.44e-15 floateq"),
]:
    calc(f"fid4:eq3-gt-fe{tag}", f"e_q3> {note}", type="e_q", power="3", comparison=">", rational=r)

# ---- arctan_q(3): '>' float-equal above Cf (only one in window) ----
calc(
    "fid4:atan-gt-fea",
    "atan3> r=Cf+2.62e-18 floateq",
    type="arctan_q",
    power="3",
    comparison=">",
    rational="115466635/92443878",
)

# ---- ln_q_square: corrected dyadics + boundary semiconvergents ----
ulp("lnsq", "ln_q_square", "2", ">", -1)  # r=C-3.652e-17 den=2^53 ok
ulp("lnsq", "ln_q_square", "2", ">", 1)  # r=C+7.45e-17  den=2^50 ok
ulp("lnsq", "ln_q_square", "2", "<", -1)
ulp("lnsq", "ln_q_square", "2", "<", 1)
for tag, r, note in [
    ("e", "45714177/95148070", "r=C-3.60e-17"),
]:
    calc(
        f"fid4:lnsq-gt-cf{tag}",
        f"lnsq> {note}",
        type="ln_q_square",
        power="2",
        comparison=">",
        rational=r,
    )
    calc(
        f"fid4:lnsq-lt-cf{tag}",
        f"lnsq< {note}",
        type="ln_q_square",
        power="2",
        comparison="<",
        rational=r,
    )

# ---- e_pi: '<'@Cf fired WD -> threshold above Cf; bracket it --------------
# Candidates: Cf+1ulp (C+4.9e-15) vs 15-digit decimal ~C+3.1e-14 (~+8.7ulp).
for j in (1, 2, 3, 5, 8, 9):
    ulp("epi", "e_pi", "1", "<", j)  # '<' WD iff r < C'_lt
ulp("epi", "e_pi", "1", ">", 1)  # scan iff C'_gt >= r
ulp("epi", "e_pi", "1", ">", 9)  # exact-C'=+3.1e-14 -> WD; float-Vf=+9ulp -> NS
ulp("epi", "e_pi", "1", ">", 10)

# ---- gauss: '<'@Cf fired WD; bracket threshold ----
ulp("gau", "gauss", "1", ">", 1)
ulp("gau", "gauss", "1", "<", 1)
ulp("gau", "gauss", "1", "<", 2)

# ---- varpi: '<'@Cf scanned (NS); check '<'@Cf-1ulp for check existence ----
ulp("var", "varpi", "1", "<", -1)
ulp("var", "varpi", "1", ">", 1)

# ---- catalan / gamma: '<'@Cf-1ulp check existence ----
ulp("cat", "catalan", "1", "<", -1)
ulp("cat", "catalan", "1", ">", 1)
ulp("gam", "gamma", "1", "<", -1)
ulp("gam", "gamma", "1", ">", 1)

# ---- golden: both scanned at Cf; quick ulp bracket ----
ulp("gol", "golden", "1", "<", -1)
ulp("gol", "golden", "1", ">", 1)

# ---- controls: '<'@-1ulp on already-pinned types (regression) ----
ulp("ln3", "ln_q", "3", "<", -1)


def main():
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                with contextlib.suppress(json.JSONDecodeError):
                    done.add(json.loads(line)["id"])
    todo = [p for p in PROBES if p["id"] not in done]
    print(f"{len(done)} done, {len(todo)} to probe", flush=True)
    last = 0.0
    for i, p in enumerate(todo, 1):
        wait = MIN_INTERVAL - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
        t0 = time.monotonic()
        try:
            resp = SESSION.post(BASE + p["endpoint"], data=p["request"], timeout=TIMEOUT)
            status, raw_body = resp.status_code, resp.text
        except (OSError, requests.RequestException) as exc:
            status, raw_body = -1, f"REQUEST_FAILED: {exc}"
        elapsed = time.monotonic() - t0
        last = time.monotonic()
        rec = dict(p)
        rec.update({"http_status": status, "raw": raw_body, "elapsed_ms": round(elapsed * 1000, 1)})
        try:
            body = json.loads(raw_body)
            rec["err"] = body.get("error", "")
            par = body.get("parameters") or {}
            rec["mn"] = [par.get("m"), par.get("n")] if par else None
        except Exception:
            pass
        with OUT.open("a") as fh:
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        short = (rec.get("err") or (f"200 mn={rec.get('mn')}")).strip()
        print(
            f"[{i}/{len(todo)}] {p['id']}: {status} ({rec['elapsed_ms']:.0f}ms) {short}", flush=True
        )


if __name__ == "__main__":
    main()
