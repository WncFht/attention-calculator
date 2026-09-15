"""Probes for the live zhuyidao.net /health app.

POSTs JSON (or raw bodies) to /health/calculate serially (>= 1 req/s) and
appends one JSON line per request to bench/data/health-probes.jsonl. Each
record carries the verbatim response text so parity_health.py can replay
requests against our Flask app and compare status + body bytes.

Usage:
    uv run python bench/probe_health.py --list
    uv run python bench/probe_health.py errors baseline thresholds met
    uv run python bench/probe_health.py --all --delay 1.05
"""

import argparse
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import requests

BASE = "https://zhuyidao.net"
PATH = "/health/calculate"
OUT = Path(__file__).parent / "data" / "health-probes.jsonl"
DELAY = 1.05
LAST_CALL = 0.0

# 必填基线；选填字段一律不携带（与前端留空行为一致：留空 -> 键缺席）
REQ = {"nickname": "Probe", "sex": "male", "age": 35,
       "height_cm": 175, "weight_kg": 70, "pal": 1.55}
FULL = {**REQ, "waist_cm": 82, "hip_cm": 96, "resting_hr": 62,
        "exercise_type": "jogging", "exercise_minutes": 30,
        "bp_context": "clinic", "systolic_bp": 118, "diastolic_bp": 76,
        "sleep_hours": 7.5, "website": ""}


def j(tag, **over):
    """A JSON POST case: REQ updated by ``over``; None removes the key."""
    payload = dict(REQ)
    for k, v in over.items():
        if v is None:
            payload.pop(k, None)
        else:
            payload[k] = v
    return {"tag": tag, "method": "POST", "path": PATH,
            "content_type": "application/json", "json": payload}


def jfull(tag, **over):
    """Same as j() but starts from the fully-filled payload."""
    payload = dict(FULL)
    for k, v in over.items():
        if v is None:
            payload.pop(k, None)
        else:
            payload[k] = v
    return {"tag": tag, "method": "POST", "path": PATH,
            "content_type": "application/json", "json": payload}


def raw(tag, body, content_type="application/json", path=PATH, method="POST"):
    """A raw-body case for transport-level probing."""
    return {"tag": tag, "method": method, "path": path,
            "content_type": content_type, "raw": body}


def send(case):
    """Send one case after the rate-limit gap; append request+response to OUT."""
    global LAST_CALL
    wait = DELAY - (time.monotonic() - LAST_CALL)
    if wait > 0:
        time.sleep(wait)
    t0 = time.monotonic()
    if case["method"] in ("GET", "OPTIONS", "HEAD"):
        r = requests.request(case["method"], BASE + case["path"], timeout=60)
    elif "raw" in case:
        headers = {} if case["content_type"] is None else {
            "Content-Type": case["content_type"]}
        r = requests.post(BASE + case["path"], data=case["raw"],
                          headers=headers, timeout=60)
    else:
        r = requests.post(BASE + case["path"], json=case["json"], timeout=60)
    LAST_CALL = time.monotonic()
    rec = {
        "ts": datetime.now(UTC).isoformat(),
        "tag": case["tag"],
        "method": case["method"],
        "path": case["path"],
        "content_type": case["content_type"],
        "json": case.get("json"),
        "raw": case.get("raw"),
        "http": r.status_code,
        "elapsed_ms": round(1000 * (LAST_CALL - t0), 1),
        "response_text": r.text,
    }
    with OUT.open("a") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return rec


def run_batch(cases, retry_429=0):
    """Send cases serially and print a one-line status per case.

    With ``retry_429`` > 0, a 429 response re-waits ``retry_429`` seconds and
    resends the same case (each attempt still logged verbatim).
    """
    n = 0
    for case in cases:
        attempts = 0
        while True:
            rec = send(case)
            n += 1
            attempts += 1
            if rec["http"] != 429 or attempts > retry_429:
                break
            print(f"[{case['tag']}] -> 429, retrying", flush=True)
            time.sleep(retry_429)
        try:
            body = json.loads(rec["response_text"])
            status = "ok" if body.get("ok") else body.get("error", rec["http"])
        except ValueError:
            status = f"http {rec['http']} (non-json)"
        print(f"[{case['tag']}] -> {status}", flush=True)
    print(f"{n} requests appended to {OUT}")


def pending_cases(names):
    """Cases in ``names`` whose latest record in OUT is still missing/429."""
    latest = {}
    if OUT.exists():
        for line in OUT.open():
            if not line.strip():
                continue
            r = json.loads(line)
            latest[r["tag"]] = r
    cases = []
    for name in names:
        for case in BATCHES[name]():
            r = latest.get(case["tag"])
            if r is None or r["http"] == 429:
                cases.append(case)
    return cases


# --- batch definitions -----------------------------------------------------


def baseline():
    """Item-array shape under each optional-field combination."""
    cases = [
        j("base:req-male"),
        j("base:req-female", sex="female"),
        j("base:waist", waist_cm=82),
        j("base:hip", hip_cm=96),
        j("base:waist-hip", waist_cm=82, hip_cm=96),
        j("base:rhr", resting_hr=62),
        j("base:exercise", exercise_type="jogging", exercise_minutes=30),
        j("base:bp-clinic", bp_context="clinic", systolic_bp=118, diastolic_bp=76),
        j("base:bp-home", bp_context="home", systolic_bp=118, diastolic_bp=76),
        j("base:sleep", sleep_hours=7.5),
        jfull("base:full-male"),
        jfull("base:full-female", sex="female"),
        j("base:pal185", pal=1.85),
        j("base:pal220", pal=2.2),
        j("base:extra-field", unknown_field=123, another="x"),
    ]
    return cases


def transport():
    """Non-JSON bodies, wrong content types, odd JSON top levels, GET."""
    good = json.dumps(REQ)
    return [
        raw("tr:form-urlencoded", "nickname=Probe&sex=male&age=35&height_cm=175&weight_kg=70&pal=1.55",
            content_type="application/x-www-form-urlencoded"),
        raw("tr:bad-json", "{"),
        raw("tr:empty-body", ""),
        raw("tr:null", "null"),
        raw("tr:list", "[]"),
        raw("tr:string", '"hi"'),
        raw("tr:number", "123"),
        raw("tr:bool", "true"),
        raw("tr:json-text-plain", good, content_type="text/plain"),
        raw("tr:json-no-ctype", good, content_type=None),
        raw("tr:json-form-ctype", good,
            content_type="application/x-www-form-urlencoded"),
        raw("tr:get-calculate", "", method="GET"),
        raw("tr:post-slash", good, path=PATH + "/"),
        raw("tr:empty-obj", "{}"),
        raw("tr:unicode-escape-nick", good.replace('"Probe"', '"\\u4f60\\u597d"')),
    ]


def errors():
    """Field-validation catalog: missing/empty/out-of-range/wrong-type."""
    cases = [
        # honeypot
        jfull("err:honeypot", website="x"),
        jfull("err:honeypot-num", website=1),
        j("err:honeypot-bare", website="http://spam"),
        # missing requireds
        j("err:no-nickname", nickname=None),
        j("err:no-sex", sex=None),
        j("err:no-age", age=None),
        j("err:no-height", height_cm=None),
        j("err:no-weight", weight_kg=None),
        j("err:no-pal", pal=None),
        # nickname charset / length
        j("err:nick-empty", nickname=""),
        j("err:nick-spaces", nickname="   "),
        j("err:nick-24", nickname="a" * 24),
        j("err:nick-25", nickname="a" * 25),
        j("err:nick-zh24", nickname="鲸" * 24),
        j("err:nick-zh25", nickname="鲸" * 25),
        j("err:nick-space-in", nickname="a b"),
        j("err:nick-at", nickname="a@b"),
        j("err:nick-dot", nickname="a.b"),
        j("err:nick-newline", nickname="a\nb"),
        j("err:nick-emoji", nickname="蓝鲸\U0001f600"),
        j("err:nick-accent", nickname="café"),
        j("err:nick-num", nickname=123),
        j("err:nick-pad-ok", nickname="  Probe  "),
        j("err:nick-pad-long", nickname="  " + "a" * 24 + "  "),
        # sex enum
        j("err:sex-other", sex="other"),
        j("err:sex-MALE", sex="MALE"),
        j("err:sex-zh", sex="男"),
        j("err:sex-empty", sex=""),
        j("err:sex-num", sex=0),
        # age
        j("err:age-str", age="35"),
        j("err:age-float", age=35.5),
        j("err:age-17", age=17),
        j("err:age-121", age=121),
        j("err:age-0", age=0),
        j("err:age-neg", age=-5),
        j("err:age-alpha", age="abc"),
        j("err:age-17.9", age=17.9),
        j("err:age-120.0", age=120.0),
        # height / weight bounds
        j("err:height-0", height_cm=0),
        j("err:height-0.9", height_cm=0.9),
        j("err:height-250.1", height_cm=250.1),
        j("err:height-str", height_cm="175"),
        j("err:height-neg", height_cm=-5),
        j("err:weight-0", weight_kg=0),
        j("err:weight-0.9", weight_kg=0.9),
        j("err:weight-500.1", weight_kg=500.1),
        j("err:weight-str", weight_kg="70"),
        # pal
        j("err:pal-str", pal="1.55"),
        j("err:pal-1.4", pal=1.4),
        j("err:pal-0", pal=0),
        j("err:pal-null", pal=None),
        j("err:pal-alpha", pal="x"),
        # optional numeric ranges
        j("err:waist-0", waist_cm=0),
        j("err:waist-0.9", waist_cm=0.9),
        j("err:waist-300.1", waist_cm=300.1),
        j("err:waist-neg", waist_cm=-1),
        j("err:waist-str", waist_cm="82"),
        j("err:hip-0", hip_cm=0),
        j("err:hip-300.1", hip_cm=300.1),
        j("err:rhr-19", resting_hr=19),
        j("err:rhr-20", resting_hr=20),
        j("err:rhr-250", resting_hr=250),
        j("err:rhr-251", resting_hr=251),
        j("err:rhr-str", resting_hr="62"),
        j("err:rhr-float", resting_hr=62.5),
        j("err:sys-49", bp_context="clinic", systolic_bp=49, diastolic_bp=76),
        j("err:sys-301", bp_context="clinic", systolic_bp=301, diastolic_bp=76),
        j("err:dia-29", bp_context="clinic", systolic_bp=118, diastolic_bp=29),
        j("err:dia-201", bp_context="clinic", systolic_bp=118, diastolic_bp=201),
        j("err:sleep-neg", sleep_hours=-0.1),
        j("err:sleep-24.1", sleep_hours=24.1),
        j("err:sleep-25", sleep_hours=25),
        j("err:sleep-str", sleep_hours="7.5"),
        j("err:exmin-0", exercise_type="jogging", exercise_minutes=0),
        j("err:exmin-1441", exercise_type="jogging", exercise_minutes=1441),
        j("err:exmin-float", exercise_type="jogging", exercise_minutes=30.5),
        j("err:exmin-str", exercise_type="jogging", exercise_minutes="30"),
        j("err:extype-bogus", exercise_type="bogus", exercise_minutes=30),
        j("err:extype-num", exercise_type=123, exercise_minutes=30),
        j("err:extype-empty", exercise_type="", exercise_minutes=30),
        j("err:bpctx-bogus", bp_context="bogus", systolic_bp=118, diastolic_bp=76),
        j("err:bpctx-empty", bp_context="", systolic_bp=118, diastolic_bp=76),
        j("err:bpctx-num", bp_context=5, systolic_bp=118, diastolic_bp=76),
        j("err:sys-str", bp_context="clinic", systolic_bp="118", diastolic_bp=76),
        # paired-field violations
        j("err:ex-type-only", exercise_type="jogging"),
        j("err:ex-min-only", exercise_minutes=30),
        j("err:bp-ctx-only", bp_context="clinic"),
        j("err:bp-sys-only", systolic_bp=118),
        j("err:bp-dia-only", diastolic_bp=76),
        j("err:bp-sys-dia", systolic_bp=118, diastolic_bp=76),
        j("err:bp-ctx-sys", bp_context="clinic", systolic_bp=118),
        j("err:bp-ctx-dia", bp_context="clinic", diastolic_bp=76),
    ]
    return cases


def error_order():
    """Paired violations to pin down the validation check order."""
    return [
        j("ord:nonick+badage", nickname=None, age=17),
        j("ord:badnick+nosex", nickname="a@b", sex=None),
        j("ord:badsex+badage", sex="other", age=17),
        j("ord:badage+badheight", age=17, height_cm=0),
        j("ord:badheight+badweight", height_cm=0, weight_kg=0),
        j("ord:badweight+badpal", weight_kg=0, pal=1.4),
        j("ord:badpal+honeypot", pal=1.4, website="x"),
        j("ord:honeypot+empty", website="x", nickname=None, sex=None,
          age=None, height_cm=None, weight_kg=None, pal=None),
        j("ord:badwaist+nosex", waist_cm=0, sex=None),
        j("ord:expair+badwaist", exercise_type="jogging", waist_cm=0),
        j("ord:bppair+expair", bp_context="clinic", exercise_type="jogging"),
        j("ord:badhip+badwaist", hip_cm=0, waist_cm=0),
        j("ord:badrhr+badwaist", resting_hr=19, waist_cm=0),
        j("ord:badnick+badwaist", nickname="a@b", waist_cm=0),
        j("ord:badsleep+badbpsys", sleep_hours=25, bp_context="clinic",
          systolic_bp=301, diastolic_bp=76),
        j("ord:badage+honeypot", age=17, website="x"),
    ]


def bmi_bounds():
    """BMI status thresholds: raw-vs-rounded comparisons at the cuts."""
    # h=200 -> bmi = w/4 exactly; probes both sides of each cut
    def c(tag, w, **kw):
        return j(tag, height_cm=200, weight_kg=w, **kw)
    return [
        c("bmi:18.475", 73.9),    # raw 18.475 -> display 18.5, status?
        c("bmi:18.5", 74.0),      # raw exactly 18.5
        c("bmi:18.525", 74.1),    # raw 18.525 -> display 18.5
        c("bmi:23.975", 95.9),    # raw 23.975 -> display 24.0, normal or overweight?
        c("bmi:24.0", 96.0),
        c("bmi:24.025", 96.1),
        c("bmi:27.975", 111.9),   # raw 27.975 -> display 28.0
        c("bmi:28.0", 112.0),
        c("bmi:28.025", 112.1),
        c("bmi:low-15", 60.0),    # deep underweight -> weight_adjustment gain
        c("bmi:high-30", 120.0),  # obese -> weight_adjustment lose
        j("bmi:extreme-60", height_cm=100, weight_kg=60),  # bmi 60 -> extreme
    ]


def thresholds():
    """Status thresholds for waist/whtr/whr/body-fat/hr/bp/sleep."""
    cases = []
    # waist_assessment male 85/90, female 80/85 (h=170 keeps whtr interesting)
    for w in (84.9, 85, 89.9, 90):
        cases.append(j(f"thr:waist-m{w}", waist_cm=w))
    for w in (79.9, 80, 84.9, 85):
        cases.append(j(f"thr:waist-f{w}", sex="female", waist_cm=w))
    # whtr 0.5 cut: waist 85 @ h170 -> exactly 0.5
    cases += [
        j("thr:whtr-0.4994", waist_cm=84.9),
        j("thr:whtr-0.5", waist_cm=85),
        j("thr:whtr-0.5006", waist_cm=85.1),
    ]
    # whr: male 0.90 (hip 100 -> waist 90 = 0.9), female 0.85
    cases += [
        j("thr:whr-m-0.899", waist_cm=89.9, hip_cm=100),
        j("thr:whr-m-0.9", waist_cm=90, hip_cm=100),
        j("thr:whr-f-0.849", sex="female", waist_cm=84.9, hip_cm=100),
        j("thr:whr-f-0.85", sex="female", waist_cm=85, hip_cm=100),
    ]
    # CUN-BAE alert boundaries male ~20/~25, female ~25/~30; plus age 80/81
    # and extreme-BMI 'out of interpretable range'
    cases += [
        j("thr:cunbae-a80", age=80),
        j("thr:cunbae-a81", age=81),
        j("thr:cunbae-f-a81", sex="female", age=81),
        j("thr:bf-m-low", height_cm=160, weight_kg=45, age=25),   # low bf
        j("thr:bf-m-20ish", height_cm=170, weight_kg=66.6, age=35),
        j("thr:bf-m-25ish", height_cm=170, weight_kg=80.9, age=35),
        j("thr:bf-m-high", height_cm=170, weight_kg=95, age=35),
        j("thr:bf-f-25ish", sex="female", height_cm=160, weight_kg=64, age=35),
        j("thr:bf-f-30ish", sex="female", height_cm=160, weight_kg=76.8, age=35),
        j("thr:bmi-extreme-low", height_cm=200, weight_kg=20),    # bmi 5
        j("thr:bmi-extreme-high", height_cm=150, weight_kg=200),  # bmi 88.9
        # resting hr 60/100 cuts
        j("thr:rhr-59", resting_hr=59),
        j("thr:rhr-60", resting_hr=60),
        j("thr:rhr-100", resting_hr=100),
        j("thr:rhr-101", resting_hr=101),
        j("thr:rhr-45", resting_hr=45),
        j("thr:rhr-120", resting_hr=120),
    ]
    # blood pressure grid: clinic thresholds 120/80 ideal, 140/90 screening,
    # grades 1/2/3; home 135/85; low <90/<60
    def bp(tag, s, d, ctx="clinic"):
        return j(tag, bp_context=ctx, systolic_bp=s, diastolic_bp=d)
    cases += [
        bp("bp:119-79", 119, 79),
        bp("bp:120-80", 120, 80),
        bp("bp:139-89", 139, 89),
        bp("bp:140-90", 140, 90),
        bp("bp:159-99", 159, 99),
        bp("bp:160-100", 160, 100),
        bp("bp:179-109", 179, 109),
        bp("bp:180-110", 180, 110),
        bp("bp:89-59", 89, 59),
        bp("bp:90-60", 90, 60),
        bp("bp:160-70", 160, 70),   # grade2 by sys only
        bp("bp:120-95", 120, 95),   # grade1 by dia only
        bp("bp:150-105", 150, 105), # sys g1, dia g2 -> worse wins?
        bp("bp:200-60", 200, 60),   # grade3 sys, low dia
        bp("bp:home-134-84", 134, 84, "home"),
        bp("bp:home-135-85", 135, 85, "home"),
        bp("bp:home-160-100", 160, 100, "home"),
        bp("bp:home-120-80", 120, 80, "home"),
        bp("bp:home-89-59", 89, 59, "home"),
        bp("bp:clinic-134-84", 134, 84, "clinic"),
        bp("bp:clinic-135-85", 135, 85, "clinic"),
    ]
    # sleep: <65 ref 7~8, >=65 ref 6~7; boundaries + direction
    cases += [
        j("slp:6.9-a35", sleep_hours=6.9),
        j("slp:7-a35", sleep_hours=7),
        j("slp:8-a35", sleep_hours=8),
        j("slp:8.1-a35", sleep_hours=8.1),
        j("slp:5.9-a65", sleep_hours=5.9, age=65),
        j("slp:6-a65", sleep_hours=6, age=65),
        j("slp:7-a65", sleep_hours=7, age=65),
        j("slp:7.1-a65", sleep_hours=7.1, age=65),
        j("slp:0-a35", sleep_hours=0),
        j("slp:24-a35", sleep_hours=24),
        j("slp:7.5-a64", sleep_hours=7.5, age=64),
        j("slp:6.5-a64", sleep_hours=6.5, age=64),
    ]
    return cases


def met():
    """Exercise MET table: every type x fixed minutes, plus rounding cases."""
    types = ["walking_slow", "walking_moderate", "walking_brisk", "jogging",
             "running_8kph", "cycling_easy", "cycling_moderate",
             "swimming_leisure", "strength_training", "yoga", "badminton",
             "table_tennis", "stair_climbing", "jump_rope"]
    cases = [j(f"met:{t}", exercise_type=t, exercise_minutes=60) for t in types]
    # minutes=1 -> kcal = met*3.5*70/200 -> x.xxxx precision check on rounding
    cases += [j(f"met:{t}-1min", exercise_type=t, exercise_minutes=1)
              for t in ("walking_slow", "jogging", "jump_rope", "yoga")]
    # minutes boundary values that are accepted
    cases += [
        j("met:min-1", exercise_type="jogging", exercise_minutes=1),
        j("met:min-1440", exercise_type="jogging", exercise_minutes=1440),
        j("met:min-45.7", exercise_type="jogging", exercise_minutes=45.7),
    ]
    return cases


def formulas():
    """Numeric ground truth for every formula across varied inputs."""
    cases = [
        # varied physiques -> regression points for every formula
        j("f:m-180-90-a45", sex="male", age=45, height_cm=180, weight_kg=90,
          waist_cm=95, hip_cm=100, resting_hr=70, pal=1.85, sleep_hours=6.5,
          bp_context="home", systolic_bp=128, diastolic_bp=82),
        j("f:f-160-55-a28", sex="female", age=28, height_cm=160, weight_kg=55,
          waist_cm=70, hip_cm=92, resting_hr=58, pal=2.2, sleep_hours=8.5,
          bp_context="clinic", systolic_bp=105, diastolic_bp=68,
          exercise_type="yoga", exercise_minutes=40),
        j("f:m-165-80-a60", sex="male", age=60, height_cm=165, weight_kg=80,
          waist_cm=98, hip_cm=99, resting_hr=80, pal=1.55, sleep_hours=5.5,
          bp_context="clinic", systolic_bp=145, diastolic_bp=95),
        j("f:f-170-48-a22", sex="female", age=22, height_cm=170, weight_kg=48,
          waist_cm=60, hip_cm=88, resting_hr=55, pal=1.85, sleep_hours=9,
          bp_context="home", systolic_bp=95, diastolic_bp=55),
        j("f:m-190-120-a50", sex="male", age=50, height_cm=190, weight_kg=120,
          waist_cm=115, hip_cm=110, resting_hr=88, pal=2.2, sleep_hours=6,
          bp_context="clinic", systolic_bp=150, diastolic_bp=100),
        j("f:f-155-40-a70", sex="female", age=70, height_cm=155, weight_kg=40,
          waist_cm=58, hip_cm=85, resting_hr=72, pal=1.55, sleep_hours=6.8,
          bp_context="home", systolic_bp=110, diastolic_bp=70),
        j("f:m-175-70-a35-p185", pal=1.85),
        j("f:m-175-70-a35-p220", pal=2.2),
        # decimal inputs: echo types + rounding
        j("f:decimals", height_cm=175.5, weight_kg=70.4, waist_cm=82.3,
          hip_cm=96.7, resting_hr=62.0, sleep_hours=7.55, age=35,
          systolic_bp=118.0, diastolic_bp=76.0, bp_context="clinic",
          exercise_type="jogging", exercise_minutes=33),
        # float age / int-looking floats
        j("f:age-35.0", age=35.0),
        j("f:h-w-int-str", height_cm=175.0, weight_kg=70.0),
        # half-even probes: mifflin male h170 w70 a35 -> 1592.5 exactly
        j("f:mifflin-.5", height_cm=170, weight_kg=70, age=35),
        # harris .5: 88.362+13.397w+4.799h-5.677a = x.5 -> w=70,h=174,a=30:
        # 88.362+937.79+834.606-170.31 = 1690.448 -> not .5; try a=33:
        # 88.362+937.79+834.606-186.741=1674.017; just keep varied points
        j("f:harris-var", height_cm=174, weight_kg=70, age=30),
        # protein range decimals: w=70.5 -> 98.7 / 141.0
        j("f:protein-dec", weight_kg=70.5),
        # exercise kcal .5-ish: met(m)*3.5*w/200*min; pick w so product .5
        j("f:kcal-round", weight_kg=80, exercise_type="jogging",
          exercise_minutes=20),  # 7.5*3.5*80/200*20 = 210.0
        j("f:kcal-round2", weight_kg=71, exercise_type="jogging",
          exercise_minutes=30),  # 7.5*3.5*71/200*30 = 279.5625 -> 280.0?
    ]
    return cases


def tips():
    """Trigger conditions for the tips/tips_en lists."""
    return [
        j("tip:overweight", height_cm=170, weight_kg=80),          # bmi 27.7
        j("tip:obese", height_cm=170, weight_kg=90),               # bmi 31.1
        j("tip:under", height_cm=180, weight_kg=55),               # bmi 17.0
        j("tip:waist-central", waist_cm=95),                       # male >=90
        j("tip:waist-pre", waist_cm=87),                           # male 85-90
        j("tip:whtr-high", waist_cm=95),                           # whtr .559
        j("tip:whr-high", waist_cm=95, hip_cm=100),                # whr .95
        j("tip:bf-alert", height_cm=160, weight_kg=90, age=40),    # high bf
        j("tip:rhr-low", resting_hr=45),
        j("tip:rhr-high", resting_hr=110),
        j("tip:bp-high", bp_context="clinic", systolic_bp=150, diastolic_bp=95),
        j("tip:bp-low", bp_context="clinic", systolic_bp=85, diastolic_bp=55),
        j("tip:bp-home-high", bp_context="home", systolic_bp=140, diastolic_bp=88),
        j("tip:sleep-low", sleep_hours=5),
        j("tip:sleep-high", sleep_hours=10),
        j("tip:exercise", exercise_type="running_8kph", exercise_minutes=45),
        j("tip:everything-bad", height_cm=160, weight_kg=95, waist_cm=105,
          hip_cm=100, resting_hr=105, age=45, sleep_hours=5,
          bp_context="clinic", systolic_bp=170, diastolic_bp=105,
          exercise_type="jogging", exercise_minutes=10),
        j("tip:everything-good-f", sex="female", height_cm=165, weight_kg=55,
          waist_cm=70, hip_cm=95, resting_hr=65, age=30, sleep_hours=7.5,
          bp_context="home", systolic_bp=115, diastolic_bp=75,
          exercise_type="yoga", exercise_minutes=30),
    ]


def misc():
    """Leftover semantics: echoes, absent-vs-null, record_id spacing."""
    return [
        # null == absent? (frontend sends null for blank optional fields!)
        j("misc:nulls", waist_cm=None, hip_cm=None, resting_hr=None,
          exercise_type=None, exercise_minutes=None, bp_context=None,
          systolic_bp=None, diastolic_bp=None, sleep_hours=None,
          website=None),
        # the page literally sends nulls: replicate exact frontend payload
        raw("misc:frontend-payload", json.dumps(
            {"nickname": "Probe", "sex": "male", "age": 35, "height_cm": 175,
             "weight_kg": 70, "waist_cm": None, "hip_cm": None,
             "resting_hr": None, "pal": 1.55, "exercise_type": None,
             "exercise_minutes": None, "bp_context": None, "systolic_bp": None,
             "diastolic_bp": None, "sleep_hours": None, "website": ""})),
        # zero is a real value, not "missing"
        j("misc:sleep-0", sleep_hours=0),
        j("misc:waist-0-with-hip", waist_cm=0, hip_cm=96),
        # nickname echo verbatim
        j("misc:nick-echo", nickname="蓝鲸-07_x"),
        # record_id: three in a row reveal the increment
        j("misc:seq-1"), j("misc:seq-2"), j("misc:seq-3"),
        # extreme but valid inputs
        j("misc:h1-w1", height_cm=1, weight_kg=1),
        j("misc:h250-w500", height_cm=250, weight_kg=500),
        j("misc:age18", age=18), j("misc:age120", age=120),
        j("misc:waist1-hip1", waist_cm=1, hip_cm=1),
        j("misc:waist300-hip300", waist_cm=300, hip_cm=300),
        j("misc:rhr20", resting_hr=20), j("misc:rhr250", resting_hr=250),
        j("misc:sys50-dia30", bp_context="clinic", systolic_bp=50,
          diastolic_bp=30),
        j("misc:sys300-dia200", bp_context="clinic", systolic_bp=300,
          diastolic_bp=200),
        j("misc:exmin1", exercise_type="jogging", exercise_minutes=1),
        j("misc:exmin1440", exercise_type="jogging", exercise_minutes=1440),
    ]


def supplement():
    """Second pass: pin down the ambiguities left by the first probe run."""
    cases = [
        # whtr 0.5 boundary at h=175 (waist/175 exactly 0.5)
        j("s2:whtr-0.4994", waist_cm=87.4),
        j("s2:whtr-0.5", waist_cm=87.5),
        j("s2:whtr-0.5006", waist_cm=87.6),
        # CUN-BAE low-side 'out of range' flag: bf<0 vs bmi bound.
        # a35 male: bf = -27.383+2.472b-0.01865b² -> 0 at bmi ~13.2
        j("s2:cb-w36", height_cm=170, weight_kg=36),    # bmi 12.46, bf -0.62
        j("s2:cb-w37", height_cm=170, weight_kg=37),    # bmi 12.80, bf ~0.0
        j("s2:cb-w38", height_cm=170, weight_kg=38),    # bmi 13.15, bf 0.62
        # a80 male: bf = -4.748+1.572b-0.0092b² -> 0 at bmi ~3.1 (value-based,
        # not bmi-based: same bmi flags at a35 but not a80)
        j("s2:cb-a80-w12", age=80, height_cm=200, weight_kg=12),  # bmi 3.0
        j("s2:cb-a80-w13", age=80, height_cm=200, weight_kg=13),  # bmi 3.25
        # cunbae upper flag? a80 vertex b~85 -> bf ~62.5 (max reachable <=80y)
        j("s2:cb-a80-hi", age=80, height_cm=170, weight_kg=247),  # bmi 85.5
        # deurenberg flag boundary: male deur = 1.2b-8.15, flagged at 87.85
        # and -2.15, normal at 27.85 -> sweep
        j("s2:deur-w110", height_cm=170, weight_kg=110),  # deur 37.5
        j("s2:deur-w116", height_cm=170, weight_kg=116),  # deur 40.0
        j("s2:deur-w120", height_cm=170, weight_kg=120),  # deur 41.7
        j("s2:deur-w125", height_cm=170, weight_kg=125),  # deur 43.7
        j("s2:deur-w130", height_cm=170, weight_kg=130),  # deur 45.8
        j("s2:deur-w140", height_cm=170, weight_kg=140),  # deur 50.0
        j("s2:deur-w160", height_cm=170, weight_kg=160),  # deur 58.3
        j("s2:deur-w200", height_cm=170, weight_kg=200),  # deur 74.9
        # female deur = 1.2b+2.65: same-value sweep to separate deur-bound
        # from bmi-bound
        j("s2:deur-f-69", sex="female", height_cm=160, weight_kg=69),   # 35.0
        j("s2:deur-f-75.4", sex="female", height_cm=160, weight_kg=75.4),  # 38.0
        j("s2:deur-f-79.7", sex="female", height_cm=160, weight_kg=79.7),  # 40.0
        j("s2:deur-f-83.9", sex="female", height_cm=160, weight_kg=83.9),  # 42.0
        j("s2:deur-f-90.3", sex="female", height_cm=160, weight_kg=90.3),  # 45.0
        # deurenberg low side: deur ~0 at bmi 6.8; check flag near 0
        j("s2:deur-lo1", height_cm=200, weight_kg=27),   # bmi 6.75, deur -0.05
        j("s2:deur-lo2", height_cm=200, weight_kg=29),   # bmi 7.25, deur 0.55
        # watson pct>100 flag: a80 h250 male w20/22/26 (tbw 23.5/24.2/25.5)
        j("s2:wat-m-w20", age=80, height_cm=250, weight_kg=20),
        j("s2:wat-m-w22", age=80, height_cm=250, weight_kg=22),
        j("s2:wat-m-w26", age=80, height_cm=250, weight_kg=26),
        j("s2:wat-f-w30", sex="female", age=80, height_cm=250, weight_kg=30),
        j("s2:wat-f-w33", sex="female", age=80, height_cm=250, weight_kg=33),
        # height tip threshold between h=1 (fires) and h=100 (silent)
        *[j(f"s2:htip-{h}", height_cm=h) for h in
          (30, 60, 90, 110, 120, 125, 130, 135, 140, 145, 150)],
        # age>80 + extreme bmi: does 超龄 suppress the 不宜 cascade + tip?
        j("s2:a81-bmi5", age=81, height_cm=200, weight_kg=20),
        # integer-vs-range check order on float inputs
        j("s2:age-119.5", age=119.5), j("s2:age-120.9", age=120.9),
        j("s2:rhr-19.5", resting_hr=19.5), j("s2:rhr-250.9", resting_hr=250.9),
        j("s2:exmin-0.5", exercise_type="jogging", exercise_minutes=0.5),
        j("s2:exmin-1440.9", exercise_type="jogging", exercise_minutes=1440.9),
        j("s2:sys-49.5", bp_context="clinic", systolic_bp=49.5,
          diastolic_bp=76),
        j("s2:sys-118.7", bp_context="clinic", systolic_bp=118.7,
          diastolic_bp=76),
        j("s2:dia-29.5", bp_context="clinic", systolic_bp=118,
          diastolic_bp=29.5),
        j("s2:dia-76.7", bp_context="clinic", systolic_bp=118,
          diastolic_bp=76.7),
        # '必须是数字' for the other numeric fields
        j("s2:h-alpha", height_cm="abc"), j("s2:w-alpha", weight_kg="abc"),
        j("s2:waist-alpha", waist_cm="abc"), j("s2:hip-alpha", hip_cm="abc"),
        j("s2:rhr-alpha", resting_hr="abc"),
        j("s2:sys-alpha", bp_context="clinic", systolic_bp="abc",
          diastolic_bp=76),
        j("s2:dia-alpha", bp_context="clinic", systolic_bp=118,
          diastolic_bp="abc"),
        j("s2:sleep-alpha", sleep_hours="abc"),
        j("s2:exmin-alpha", exercise_type="jogging", exercise_minutes="abc"),
        # empty strings: treated as absent or invalid?
        j("s2:h-empty", height_cm=""), j("s2:waist-empty", waist_cm=""),
        j("s2:sleep-empty", sleep_hours=""),
        # sex normalization edges
        j("s2:sex-Male", sex="Male"), j("s2:sex-pad", sex=" male "),
        j("s2:sex-FEMALE", sex="FEMALE"), j("s2:sex-bool", sex=True),
        # nickname: order length-vs-charset, more charsets, non-str types
        j("s2:nick-25-bad", nickname="a" * 24 + "@"),
        j("s2:nick-jp", nickname="日本語テスト"),
        j("s2:nick-greek", nickname="αβγ"),
        j("s2:nick-pad-25", nickname="  " + "a" * 25 + "  "),
        j("s2:nick-list", nickname=["a"]),
        # pal string variant
        j("s2:pal-str220", pal="2.20"),
        # tdee: raw vs pre-rounded mifflin (1592.5 * pal)
        j("s2:mifflin-p185", height_cm=170, weight_kg=70, age=35, pal=1.85),
        j("s2:mifflin-p220", height_cm=170, weight_kg=70, age=35, pal=2.2),
        # NaN / Infinity raw JSON bodies
        raw("s2:age-nan", '{"nickname":"P","sex":"male","age":NaN,'
            '"height_cm":175,"weight_kg":70,"pal":1.55}'),
        raw("s2:h-inf", '{"nickname":"P","sex":"male","age":35,'
            '"height_cm":Infinity,"weight_kg":70,"pal":1.55}'),
        raw("s2:w-neg-inf", '{"nickname":"P","sex":"male","age":35,'
            '"height_cm":175,"weight_kg":-Infinity,"pal":1.55}'),
        # bool / floaty-string coercion
        j("s2:age-bool", age=True),
        j("s2:h-exp", height_cm="1e2"),
        j("s2:h-pad", height_cm=" 175 "),
        # bp category priority: low-vs-grade conflicts
        j("s2:bp-85-95", bp_context="clinic", systolic_bp=85, diastolic_bp=95),
        j("s2:bp-180-50", bp_context="clinic", systolic_bp=180, diastolic_bp=50),
        j("s2:bp-95-35", bp_context="clinic", systolic_bp=95, diastolic_bp=35),
        j("s2:bp-134-60", bp_context="clinic", systolic_bp=134, diastolic_bp=60),
        # does waist tip fire on whtr/whr alone? waist<85 + whtr>=0.5
        j("s2:whtr-only", height_cm=160, waist_cm=82),   # whtr 0.5125, w<85
        j("s2:whr-only", waist_cm=82, hip_cm=90),        # whr 0.911, w<85
        # bp trio / exercise pair with bogus context only
        j("s2:bpctx-bogus-solo", bp_context="bogus"),
        j("s2:extype-bogus-solo", exercise_type="bogus"),
        # charset=utf-8 suffix on JSON content type
        raw("s2:ct-charset", json.dumps(REQ),
            content_type="application/json; charset=utf-8"),
        # route edges
        raw("s2:get-health", "", path="/health", method="GET"),
        raw("s2:get-health-slash", "", path="/health/", method="GET"),
        raw("s2:get-health-en", "", path="/health/en", method="GET"),
        raw("s2:get-health-en-slash", "", path="/health/en/", method="GET"),
        raw("s2:get-health-static", "", path="/health/static/bg1.png",
            method="GET"),
        raw("s2:options-calc", "", path="/health/calculate", method="OPTIONS"),
        raw("s2:post-health-slash", json.dumps(REQ), path="/health/"),
        raw("s2:get-health-xyz", "", path="/health/xyz", method="GET"),
        raw("s2:post-health-en", json.dumps(REQ), path="/health/en"),
    ]
    return cases


BATCHES = {
    "baseline": baseline,
    "transport": transport,
    "errors": errors,
    "error_order": error_order,
    "bmi_bounds": bmi_bounds,
    "thresholds": thresholds,
    "met": met,
    "formulas": formulas,
    "tips": tips,
    "misc": misc,
    "supplement": supplement,
}


def main():
    global DELAY
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("batches", nargs="*", choices=BATCHES)
    p.add_argument("--all", action="store_true", help="run every batch in order")
    p.add_argument("--list", action="store_true")
    p.add_argument("--delay", type=float, default=DELAY)
    p.add_argument("--retry-429", type=float, default=0, metavar="SECS",
                   help="on HTTP 429 wait SECS and resend the same case")
    p.add_argument("--pending", action="store_true",
                   help="only send cases whose latest record is missing/429")
    args = p.parse_args()
    if args.list:
        total = 0
        for name, f in BATCHES.items():
            print(f"{name}: {len(f())} cases")
            total += len(f())
        print(f"total: {total}")
        return
    DELAY = args.delay
    names = list(BATCHES) if args.all else args.batches
    if not names:
        p.error("give batch names or --all")
    if args.pending:
        cases = pending_cases(names)
        print(f"{len(cases)} pending cases")
        run_batch(cases, retry_429=args.retry_429)
        return
    for name in names:
        run_batch(BATCHES[name](), retry_429=args.retry_429)


if __name__ == "__main__":
    sys.exit(main())
