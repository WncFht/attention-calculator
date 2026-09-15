"""Flip-threshold probes: map where the site switches 200 / 方向反了 / 未找到.

For each type we sweep the bound at sub-ulp granularity around float64(C):
  '>' claims, TRUE side : r = C - 10^-k  (float64-equal for k >= ~16)
  '>' claims, FALSE float-equal : r in (C, Cf+half-ulp)  -- exact nonpos exists;
          方向反了 => mid-scan nonpos check (exact or numeric); 未找到 => none
  '<' claims, FALSE float-equal : r = Cf - 10^-k  -- exact nonneg impossible;
          200 => numerically-faked proof
  '<' claims, TRUE float-equal  : r in (C, Cf+half-ulp)
  1-ulp-off controls both directions.

Appends to bench/data/fidelity-probes.jsonl; resumable via id skip.
"""
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from attention_calculator.integrand import constant_mpf  # noqa: E402

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "fidelity-probes.jsonl"
MIN_INTERVAL = 1.05
TIMEOUT = 100

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (fidelity probes)"

PROBES: list[dict] = []


def calc(pid, note, **form):
    fields = {k: v for k, v in form.items() if v is not None}
    PROBES.append({"id": pid, "endpoint": "/calculate", "method": "POST",
                   "kind": "form", "request": fields, "note": note})


def bound_str(r: Fraction) -> str:
    return f"{r.numerator}/{r.denominator}"


def cf_minus(kind, power, kexp):
    """r = float64(C) - 10^-k  (rounds to Cf; TRUE '>' / FALSE '<')."""
    C = constant_mpf(kind, Fraction(power))
    Cf = float(C)
    return Fraction(Cf) - Fraction(10) ** -kexp


def cf_plus(kind, power, num, den_pow):
    """r = float64(C) + num*10^-den_pow."""
    C = constant_mpf(kind, Fraction(power))
    Cf = float(C)
    return Fraction(Cf) + Fraction(num) * Fraction(10) ** -den_pow


def cf_frac(kind, power):
    return Fraction(float(constant_mpf(kind, Fraction(power))))


def cf_nextafter(kind, power, direction):
    Cf = float(constant_mpf(kind, Fraction(power)))
    return Fraction(math.nextafter(Cf, -math.inf if direction < 0 else math.inf))


def f(kind, power):
    return float(constant_mpf(kind, Fraction(power)))


# ---------------- pi ----------------
K = "pi"; P = "1"
# '>' TRUE float-equal sweep: r = Cf - 10^-k  (gap = (C-Cf) + 10^-k ~ 1.2e-16..)
for k in (17, 20, 24, 28, 30, 32, 34):
    calc(f"fid:pi-gt-true-{k}", f"pi> r=Cf-10^-{k} true floateq", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))
# '>' FALSE float-equal: r = Cf + {1.5e-16, 2e-16} in (pi, Cf+half-ulp)
calc("fid:pi-gt-false-a", "pi> r=Cf+1.5e-16 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 15, 17)))
calc("fid:pi-gt-false-b", "pi> r=Cf+2e-16 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 2, 16)))
# '<' FALSE float-equal sweep: r = Cf - 10^-k (r<C, claim 'pi<r' false)
for k in (17, 25, 32):
    calc(f"fid:pi-lt-false-{k}", f"pi< r=Cf-10^-{k} false floateq", type=K,
         power=P, comparison="<", rational=bound_str(cf_minus(K, P, k)))
# '<' TRUE float-equal: r = Cf+1.5e-16 (>pi)
calc("fid:pi-lt-true-a", "pi< r=Cf+1.5e-16 true floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_plus(K, P, 15, 17)))
# 1-ulp controls
calc("fid:pi-lt-1ulpdn", "pi< r=Cf-1ulp false floatvis", type=K, power=P,
     comparison="<", rational=bound_str(cf_nextafter(K, P, -1)))
calc("fid:pi-gt-1ulpup", "pi> r=Cf+1ulp false floatvis", type=K, power=P,
     comparison=">", rational=bound_str(cf_nextafter(K, P, 1)))
calc("fid:pi-gt-cf", "pi> r=Cf exact true", type=K, power=P,
     comparison=">", rational=bound_str(cf_frac(K, P)))
calc("fid:pi-lt-cf", "pi< r=Cf exact false", type=K, power=P,
     comparison="<", rational=bound_str(cf_frac(K, P)))

# ---------------- e ----------------
K = "e"; P = "1"
for k in (17, 24, 30, 34):
    calc(f"fid:e-gt-true-{k}", f"e> r=Cf-10^-{k}", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))
calc("fid:e-gt-false-a", "e> r=Cf+1.5e-16 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 15, 17)))
calc("fid:e-gt-false-b", "e> r=Cf+2e-16 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 2, 16)))
calc("fid:e-lt-false-30", "e< r=Cf-10^-30 false floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_minus(K, P, 30)))
calc("fid:e-lt-true-a", "e< r=Cf+1.5e-16 true floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_plus(K, P, 15, 17)))
calc("fid:e-lt-1ulpdn", "e< r=Cf-1ulp", type=K, power=P,
     comparison="<", rational=bound_str(cf_nextafter(K, P, -1)))

# ---------------- sin_q ----------------
K = "sin_q"; P = "1"
for k in (17, 25, 32):
    calc(f"fid:sin-gt-true-{k}", f"sin> r=Cf-10^-{k}", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))
# false '>' float-equal: r in (sin1= Cf+1.8e-18, Cf+5.5e-17)
calc("fid:sin-gt-false-a", "sin> r=Cf+2e-17 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 2, 17)))
calc("fid:sin-gt-false-b", "sin> r=Cf+4e-17 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 4, 17)))
calc("fid:sin-lt-false-30", "sin< r=Cf-10^-30 false floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_minus(K, P, 30)))
calc("fid:sin-lt-true-a", "sin< r=Cf+2e-17 true floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_plus(K, P, 2, 17)))
calc("fid:sin-lt-1ulpdn", "sin< r=Cf-1ulp", type=K, power=P,
     comparison="<", rational=bound_str(cf_nextafter(K, P, -1)))

# ---------------- ln_q(2) ----------------
K = "ln_q"; P = "2"
for k in (17, 25, 32):
    calc(f"fid:ln-gt-true-{k}", f"ln> r=Cf-10^-{k}", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))
calc("fid:ln-gt-false-a", "ln> r=Cf+4e-17 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 4, 17)))
calc("fid:ln-lt-false-30", "ln< r=Cf-10^-30 false floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_minus(K, P, 30)))
calc("fid:ln-lt-true-a", "ln< r=Cf+4e-17 true floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_plus(K, P, 4, 17)))

# ---------------- e_q(1) ----------------
K = "e_q"; P = "1"
calc("fid:eq1-gt-false-a", "e_q1> r=Cf+1.5e-16 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 15, 17)))
for k in (17, 30):
    calc(f"fid:eq1-gt-true-{k}", f"e_q1> r=Cf-10^-{k}", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))
calc("fid:eq1-lt-false-30", "e_q1< r=Cf-10^-30 false floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_minus(K, P, 30)))

# ---------------- e_q(3): Cf > C ----------------
K = "e_q"; P = "3"
calc("fid:eq3-gt-false-a", "e_q3> r=Cf-1e-16 false floateq (Cf>C)", type=K,
     power=P, comparison=">", rational=bound_str(cf_minus(K, P, 16)))
for k in (17, 30):
    calc(f"fid:eq3-gt-true-{k}", f"e_q3> r=Cf-10^-{k}... wait Cf>C so r=Cf-1e-17 could be >C", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))

# ---------------- cos_q(1): Cf > C ----------------
K = "cos_q"; P = "1"
calc("fid:cos-gt-cf", "cos> r=Cf exact: FALSE floateq (Cf>C)", type=K, power=P,
     comparison=">", rational=bound_str(cf_frac(K, P)))
calc("fid:cos-lt-cf", "cos< r=Cf exact: TRUE floateq", type=K, power=P,
     comparison="<", rational=bound_str(cf_frac(K, P)))
for k in (17, 30):
    calc(f"fid:cos-gt-true-{k}", f"cos> r=Cf-10^-{k}", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))

# ---------------- zeta3 ----------------
K = "zeta3"; P = "1"
calc("fid:zeta-gt-false-a", "zeta> r=Cf+8e-17 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 8, 17)))
for k in (17, 30):
    calc(f"fid:zeta-gt-true-{k}", f"zeta> r=Cf-10^-{k}", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))

# ---------------- arctan_q(3): Cf > C by 2.2e-18 ----------------
K = "arctan_q"; P = "3"
calc("fid:atan-gt-false-a", "atan3> r=Cf-1e-18 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_minus(K, P, 18)))
calc("fid:atan-gt-true-30", "atan3> r=Cf-10^-30 true", type=K, power=P,
     comparison=">", rational=bound_str(cf_minus(K, P, 30)))

# ---------------- ln_q_square(2): Cf > C ----------------
K = "ln_q_square"; P = "2"
calc("fid:lnsq-gt-false-a", "lnsq> r=Cf-1e-17 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_minus(K, P, 17)))
for k in (18, 30):
    calc(f"fid:lnsq-gt-true-{k}", f"lnsq> r=Cf-10^-{k}", type=K, power=P,
         comparison=">", rational=bound_str(cf_minus(K, P, k)))

# ---------------- pi_n(5/2) ----------------
K = "pi_n"; P = "5/2"
calc("fid:pin-gt-false-a", "pi_n> r=Cf+5e-16 false floateq", type=K, power=P,
     comparison=">", rational=bound_str(cf_plus(K, P, 5, 16)))
calc("fid:pin-gt-true-30", "pi_n> r=Cf-10^-30 true", type=K, power=P,
     comparison=">", rational=bound_str(cf_minus(K, P, 30)))


# ---------------- runner ----------------

def main():
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                try:
                    done.add(json.loads(line)["id"])
                except json.JSONDecodeError:
                    pass
    todo = [p for p in PROBES if p["id"] not in done]
    print(f"{len(done)} done, {len(todo)} to probe", flush=True)
    last = 0.0
    for i, p in enumerate(todo, 1):
        wait = MIN_INTERVAL - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
        t0 = time.monotonic()
        try:
            resp = SESSION.post(BASE + p["endpoint"], data=p["request"],
                                timeout=TIMEOUT)
            status, raw_body = resp.status_code, resp.text
        except (OSError, requests.RequestException) as exc:
            status, raw_body = -1, f"REQUEST_FAILED: {exc}"
        elapsed = time.monotonic() - t0
        last = time.monotonic()
        rec = dict(p)
        rec.update({"http_status": status, "raw": raw_body,
                    "elapsed_ms": round(elapsed * 1000, 1)})
        try:
            rec["err"] = json.loads(raw_body).get("error", "")
            par = json.loads(raw_body).get("parameters", {})
            rec["mn"] = (par.get("m"), par.get("n")) if par else None
        except Exception:
            pass
        with OUT.open("a") as fh:
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        short = raw_body[:110].replace("\n", " ")
        print(f"[{i}/{len(todo)}] {p['id']}: {status} ({rec['elapsed_ms']:.0f}ms) {short}",
              flush=True)


if __name__ == "__main__":
    main()
