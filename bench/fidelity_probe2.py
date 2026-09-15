"""Flip-threshold probes v2: CF-generated bounds (num,den < 1e16).

For each type: '>' TRUE float-equal sweep (gap ladder), '>' FALSE float-equal
(exact nonpos exists / fake-200 hunt), '<' FALSE float-equal (fake-200 hunt),
'<' TRUE float-equal. Appends to bench/data/fidelity-probes.jsonl.
"""
import contextlib
import json
import time
from pathlib import Path

import requests

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "fidelity-probes.jsonl"
MIN_INTERVAL = 1.05
TIMEOUT = 120

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (fidelity probes v2)"

PROBES: list[dict] = []


def calc(pid, note, **form):
    fields = {k: v for k, v in form.items() if v is not None}
    PROBES.append({"id": pid, "endpoint": "/calculate", "method": "POST",
                   "kind": "form", "request": fields, "note": note})


def b(s):
    return s  # already "p/q"


# ---- pi: gap ladder, both directions ----
PI_TRUE = [  # r < pi, '>' TRUE, gap ascending
    ("4.9e-32", "6134899525417045/1952799169684491"),
    ("3.8e-30", "428224593349304/136308121570117"),
    ("8.4e-29", "288469374822515/91822653867264"),
    ("3.1e-28", "148714156295726/47337186164411"),
    ("7.7e-27", "8958937768937/2851718461558"),
    ("3.2e-25", "3587785776203/1142027682075"),
]
PI_FALSE = [  # r > pi
    ("2.3e-31", "5706674932067741/1816491048114374"),
    ("1.4e-30", "4422001152019829/1407566683404023"),
    ("1.97e-30", "3993776558670525/1271258561833906"),
]
for i, (g, r) in enumerate(PI_TRUE):
    calc(f"fid2:pi-gt-T{i}", f"pi> true floateq gap {g}", type="pi", power="1",
         comparison=">", rational=b(r))
for i, (g, r) in enumerate(PI_FALSE):
    calc(f"fid2:pi-gt-F{i}", f"pi> FALSE floateq gap -{g}", type="pi", power="1",
         comparison=">", rational=b(r))
    calc(f"fid2:pi-lt-T{i}", f"pi< TRUE floateq gap -{g}", type="pi", power="1",
         comparison="<", rational=b(r))
for i, (g, r) in enumerate(PI_TRUE[:3]):
    calc(f"fid2:pi-lt-F{i}", f"pi< FALSE floateq gap {g}", type="pi", power="1",
         comparison="<", rational=b(r))

# ---- e ----
E_TRUE = [("6.5e-32", "2124008553358849/781379079653017"),
          ("3.4e-30", "1038929163353808/382200680031313"),
          ("1.06e-29", "992778936702575/365222960440922"),
          ("1.5e-28", "1085079390005041/399178399621704")]
E_FALSE = [("2.98e-31", "9581113603440437/3524694718233772"),
           ("1.02e-30", "3209087943363890/1180557479274721"),
           ("1.5e-28", "46150226651233/16977719590391")]
for i, (g, r) in enumerate(E_TRUE):
    calc(f"fid2:e-gt-T{i}", f"e> true floateq gap {g}", type="e", power="1",
         comparison=">", rational=b(r))
for i, (g, r) in enumerate(E_FALSE):
    calc(f"fid2:e-gt-F{i}", f"e> FALSE floateq gap -{g}", type="e", power="1",
         comparison=">", rational=b(r))
calc("fid2:e-lt-F0", "e< FALSE floateq 6.5e-32", type="e", power="1",
     comparison="<", rational=b(E_TRUE[0][1]))
calc("fid2:e-lt-T0", "e< TRUE floateq 2.98e-31", type="e", power="1",
     comparison="<", rational=b(E_FALSE[0][1]))

# ---- sin_q ----
SIN_TRUE = [("9.6e-33", "8199121568317184/9743795943467975"),
            ("2.15e-32", "7232065840715273/8594551649771164"),
            ("8.93e-32", "4330898657909540/5146818768680731"),
            ("1.38e-31", "3363842930307629/3997574474983920")]
SIN_FALSE = [("7.97e-32", "967055727601911/1149244293696811"),
             ("1.53e-30", "504379980100015/599402699803324")]
for i, (g, r) in enumerate(SIN_TRUE):
    calc(f"fid2:sin-gt-T{i}", f"sin> true floateq gap {g}", type="sin_q",
         power="1", comparison=">", rational=b(r))
for i, (g, r) in enumerate(SIN_FALSE):
    calc(f"fid2:sin-gt-F{i}", f"sin> FALSE floateq gap -{g}", type="sin_q",
         power="1", comparison=">", rational=b(r))
    calc(f"fid2:sin-lt-T{i}", f"sin< TRUE floateq gap -{g}", type="sin_q",
         power="1", comparison="<", rational=b(r))
calc("fid2:sin-lt-F0", "sin< FALSE floateq 9.6e-33", type="sin_q", power="1",
     comparison="<", rational=b(SIN_TRUE[0][1]))

# ---- ln_q(2) ----
LN_TRUE = [("1.3e-32", "1554903831458736/2243252046704767"),
           ("2.57e-31", "1266537308075249/1827227093461019"),
           ("1.36e-30", "689804261308275/995177186973523"),
           ("1.37e-29", "113071214541301/163127280486027")]
LN_FALSE = [("3.43e-32", "6507981849218431/9389033140062816"),
            ("7.78e-32", "3398174186300959/4902529046653282"),
            ("1.06e-30", "288366523383487/416024953243748")]
for i, (g, r) in enumerate(LN_TRUE):
    calc(f"fid2:ln-gt-T{i}", f"ln> true floateq gap {g}", type="ln_q", power="2",
         comparison=">", rational=b(r))
for i, (g, r) in enumerate(LN_FALSE):
    calc(f"fid2:ln-gt-F{i}", f"ln> FALSE floateq gap -{g}", type="ln_q",
         power="2", comparison=">", rational=b(r))
calc("fid2:ln-lt-T0", "ln< TRUE floateq 3.43e-32", type="ln_q", power="2",
     comparison="<", rational=b(LN_FALSE[0][1]))

# ---- e_q(1) : same constant e ----
EQ1_FALSE = E_FALSE[:2]
for i, (g, r) in enumerate(E_TRUE[:3]):
    calc(f"fid2:eq1-gt-T{i}", f"e_q1> true floateq gap {g}", type="e_q",
         power="1", comparison=">", rational=b(r))
for i, (g, r) in enumerate(EQ1_FALSE):
    calc(f"fid2:eq1-gt-F{i}", f"e_q1> FALSE floateq gap -{g}", type="e_q",
         power="1", comparison=">", rational=b(r))

# ---- e_q(3): C_f > C ----
EQ3_TRUE = [("2.19e-30", "9060589063489840/451100167157089"),
            ("2.27e-29", "4342874389265102/216218984131387"),
            ("1.96e-28", "1609176767193097/80116193724221")]
EQ3_FALSE = [("1.67e-29", "2358857337112369/117440591512851"),
             ("4.73e-28", "374840284959636/18662198894315")]
for i, (g, r) in enumerate(EQ3_TRUE[1:]):
    calc(f"fid2:eq3-gt-T{i+1}", f"e_q3> true floateq gap {g}", type="e_q",
         power="3", comparison=">", rational=b(r))
for i, (g, r) in enumerate(EQ3_FALSE):
    calc(f"fid2:eq3-gt-F{i}", f"e_q3> FALSE floateq gap -{g}", type="e_q",
         power="3", comparison=">", rational=b(r))

# ---- cos_q ----
COS_FALSE = [("1.54e-31", "5365928276967407/9931344394959804"),
             ("1.9e-31", "4486613785117493/8303895312659367")]
for i, (g, r) in enumerate(COS_FALSE):
    calc(f"fid2:cos-gt-F{i}", f"cos> FALSE floateq gap -{g}", type="cos_q",
         power="1", comparison=">", rational=b(r))
calc("fid2:cos-lt-F0", "cos< FALSE floateq 3.18e-32", type="cos_q", power="1",
     comparison="<", rational="293104830616638/542483027433479")
calc("fid2:cos-lt-T0", "cos< TRUE floateq 1.54e-31", type="cos_q", power="1",
     comparison="<", rational=COS_FALSE[0][1])

# ---- zeta3 ----
Z_FALSE = [("1.02e-32", "9890949372930056/8228353705163039"),
           ("1.46e-31", "3228169011047533/2685537600227014")]
for i, (g, r) in enumerate(Z_FALSE):
    calc(f"fid2:zeta-gt-F{i}", f"zeta> FALSE floateq gap -{g}", type="zeta3",
         power="1", comparison=">", rational=b(r))
calc("fid2:zeta-gt-T0", "zeta> true floateq 2.45e-31", type="zeta3", power="1",
     comparison=">", rational="1144870450278330/952426168236337")

# ---- arctan_q(3) ----
calc("fid2:atan-gt-F0", "atan3> FALSE floateq 5.5e-31", type="arctan_q",
     power="3", comparison=">", rational="9972443475689970/7984049661000163")
calc("fid2:atan-gt-T0", "atan3> true floateq 4.36e-28", type="arctan_q",
     power="3", comparison=">", rational="28688140985470/22968046183277")

# ---- ln_q_square(2) ----
LNSQ_TRUE = [("1.06e-32", "4646620020445232/9671330777074349"),
             ("4.19e-32", "1585407615904731/3299818233994160"),
             ("1.37e-30", "109602827268961/228123924908131")]
LNSQ_FALSE = [("5.67e-33", "3061212404540501/6371512543080189"),
              ("1.71e-31", "1366201961366809/2843570384177898")]
for i, (g, r) in enumerate(LNSQ_TRUE[1:]):
    calc(f"fid2:lnsq-gt-T{i+1}", f"lnsq> true floateq gap {g}", type="ln_q_square",
         power="2", comparison=">", rational=b(r))
for i, (g, r) in enumerate(LNSQ_FALSE):
    calc(f"fid2:lnsq-gt-F{i}", f"lnsq> FALSE floateq gap -{g}",
         type="ln_q_square", power="2", comparison=">", rational=b(r))

# ---- pi_n(5/2) ----
calc("fid2:pin-gt-F0", "pi_n> FALSE floateq 9.14e-30", type="pi_n", power="5/2",
     comparison=">", rational="2755163350798137/157497139735537")
calc("fid2:pin-gt-T0", "pi_n> true floateq 4.06e-30", type="pi_n", power="5/2",
     comparison=">", rational="8417507566830365/481181402581438")
calc("fid2:pin-gt-T1", "pi_n> true floateq 1.58e-26", type="pi_n", power="5/2",
     comparison=">", rational="133169423484989/7612544386176")


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
        print(f"[{i}/{len(todo)}] {p['id']}: {status} ({rec['elapsed_ms']:.0f}ms) {short}",
              flush=True)


if __name__ == "__main__":
    main()
