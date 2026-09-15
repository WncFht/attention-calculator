"""Battery 3: pin the per-type compare constant C' at ulp granularity.

Model under test: the site decides 方向反了 by an EXACT compare of the parsed
rational bound r against a per-type constant C' (about float64(C)):
  '>' claim -> WrongDirection iff r > C'
  '<' claim -> WrongDirection iff r < C'
  equality  -> exact scan (200 or 未找到)

For each type we probe r = Fraction(nextafter(Cf, dir, j)) — exact dyadics,
always num/den < 1e16 (except ln_q_square j=0 which is substituted by a CF
semiconvergent below). Predictions from the exact solver are in the notes.

Appends to bench/data/fidelity-probes.jsonl; resumable via id skip.
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
from attention_calculator.integrand import constant_mpf

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "fidelity-probes.jsonl"
MIN_INTERVAL = 1.05
TIMEOUT = 120

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (fidelity probes v3)"

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
    """Fraction of the j-th nextafter of float64(C) (j=0 -> Cf itself)."""
    f = float(constant_mpf(kind, Fraction(power)))
    for _ in range(abs(j)):
        f = math.nextafter(f, math.inf if j > 0 else -math.inf)
    return Fraction(f)


def ulp_sweep(tag, kind, power, gt_js, lt_js):
    """'>' probes at gt_js offsets, '<' probes at lt_js offsets from Cf."""
    for j in gt_js:
        r = dyadic(kind, power, j)
        calc(
            f"fid3:{tag}-gt-{j:+d}",
            f"{kind}({power})> Cf{j:+d}ulp r-Cf={j}ulp",
            type=kind,
            power=power,
            comparison=">",
            rational=f"{r.numerator}/{r.denominator}",
        )
    for j in lt_js:
        r = dyadic(kind, power, j)
        calc(
            f"fid3:{tag}-lt-{j:+d}",
            f"{kind}({power})< Cf{j:+d}ulp r-Cf={j}ulp",
            type=kind,
            power=power,
            comparison="<",
            rational=f"{r.numerator}/{r.denominator}",
        )


# ---- e: C' in (Cf-1ulp, C-6.5e-32); test Cf exactly ----
ulp_sweep("e", "e", "1", gt_js=(-1, 0, 1), lt_js=(0, 1))

# ---- e_q(1): same constant e ----
ulp_sweep("eq1", "e_q", "1", gt_js=(0,), lt_js=(0,))

# ---- e_q(3): C' > C+4.7e-28; bracket vs Cf ----
ulp_sweep("eq3", "e_q", "3", gt_js=(-1, 0, 1), lt_js=(-1, 0, 1))

# ---- sin_q(1): C' in (Cf-1ulp, C-9.6e-33) ----
ulp_sweep("sin", "sin_q", "1", gt_js=(-1, 0, 1), lt_js=(0, 1))

# ---- ln_q(2): C' < C-1.3e-32, lower bound unknown ----
ulp_sweep("ln", "ln_q", "2", gt_js=(-1, 0, 1), lt_js=(-1, 0, 1))

# ---- zeta3: C' > C+1.46e-31; bracket vs Cf, Cf+1ulp, Cf+2ulp ----
ulp_sweep("zeta", "zeta3", "1", gt_js=(0, 1, 2), lt_js=(0, 1, 2))

# ---- arctan_q(3): C' > C+5.5e-31; bracket vs Cf ----
ulp_sweep("atan", "arctan_q", "3", gt_js=(0, 1), lt_js=(-1, 0, 1))

# ---- ln_q_square(2): C' < C-4.2e-32; Cf dyadics hit den cap at j=0,+-2 ----
ulp_sweep("lnsq", "ln_q_square", "2", gt_js=(-1, 1), lt_js=(-1, 1))
# CF semiconvergents float-equal to Cf, gaps straddling candidate C' values:
#   gap>0 => r<C ; predicted WD iff r > C'
for tag, gapnote, r in [
    ("g4e-18", "r=C-4.37e-18", "159506815/331992537"),
    ("g1e-18", "r=C-1.03e-18", "432806268/900829541"),
    ("g3e-19", "r=C-2.75e-19", "706105721/1469666545"),
    ("g2e-20", "r=C-2.99e-20", "2664916069/5546673643"),
]:
    calc(
        f"fid3:lnsq-gt-cf{tag}",
        f"lnsq> {gapnote}",
        type="ln_q_square",
        power="2",
        comparison=">",
        rational=r,
    )

# ---- pi_n(5/2): C' < C-4.06e-30; watch for r^v-vs-C^v compare ----
ulp_sweep("pin", "pi_n", "5/2", gt_js=(0, 1), lt_js=(-1, 0))

# ---- generic spot checks on remaining types: '>' and '<' at Cf ----
for tag, kind, power in [
    ("cat", "catalan", "1"),
    ("gam", "gamma", "1"),
    ("gol", "golden", "1"),
    ("epi", "e_pi", "1"),
    ("var", "varpi", "1"),
    ("gau", "gauss", "1"),
    ("sin2", "sin_q", "1/2"),
    ("ln3", "ln_q", "3"),
]:
    ulp_sweep(tag, kind, power, gt_js=(0,), lt_js=(0,))


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
            rec["a_val"] = par.get("a_val")
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
