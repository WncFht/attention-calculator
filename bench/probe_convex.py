"""Replayable probes of the live zhuyidao.net /convex/prove API.

Each batch POSTs /convex/prove cases serially (>= 1 s apart) and appends one
JSON line per HTTP request to bench/data/convex-probes.jsonl. Each record keeps
the full request (method, endpoint, form fields, encoding) plus the verbatim
raw response body and the parsed JSON, keyed by ``tag``.

Re-running a batch re-fetches and appends fresh rows. Usage:

    uv run python bench/probe_convex.py --list
    uv run python bench/probe_convex.py atoms errors
    uv run python bench/probe_convex.py --all --delay 1.1
"""

import argparse
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import requests

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "convex-probes.jsonl"
DELAY = 1.05
LAST_CALL = 0.0


def prove(inequality, tag, line=None, extra=None):
    """One POST /convex/prove case; ``tag`` marks the probing purpose."""
    form = {"inequality": inequality}
    if line is not None:
        form["line"] = line
    if extra:
        form.update(extra)
    return {
        "tag": tag,
        "method": "POST",
        "endpoint": "/convex/prove",
        "encode": "form",
        "form": form,
    }


def raw_case(tag, method, endpoint, form=None, encode="form", body=None, json_body=None):
    """Arbitrary request shape for transport/malformed probing."""
    return {
        "tag": tag,
        "method": method,
        "endpoint": endpoint,
        "encode": encode,
        "form": form or {},
        "body": body,
        "json": json_body,
    }


def send_once(case):
    """Dispatch one case to the server; returns the response object."""
    url = BASE + case["endpoint"]
    if case["method"] == "GET":
        return requests.get(url, params=case["form"], timeout=60)
    if case["method"] == "HEAD":
        return requests.head(url, params=case["form"], timeout=60)
    if case["method"] == "OPTIONS":
        return requests.options(url, timeout=60)
    if case["method"] == "PUT":
        return requests.put(url, data=case["form"], timeout=60)
    if case["encode"] == "multipart":
        files = {k: (None, v) for k, v in case["form"].items()}
        return requests.post(url, files=files, timeout=60)
    if case["encode"] == "json":
        return requests.post(url, json=case["json"], timeout=60)
    if case["encode"] == "raw":
        return requests.post(
            url, data=case["body"], headers={"Content-Type": "text/plain"}, timeout=60
        )
    return requests.post(url, data=case["form"], timeout=60)


def send(case):
    """Send one case after the rate-limit gap; append request+response to OUT."""
    global LAST_CALL
    wait = DELAY - (time.monotonic() - LAST_CALL)
    if wait > 0:
        time.sleep(wait)
    t0 = time.monotonic()
    r = None
    for attempt in range(4):
        try:
            r = send_once(case)
            break
        except requests.RequestException as e:
            print(f"  retry {attempt + 1} after {type(e).__name__}", flush=True)
            time.sleep(5 * (attempt + 1))
    if r is None:
        rec = {
            "ts": datetime.now(UTC).isoformat(),
            "tag": case["tag"],
            "method": case["method"],
            "endpoint": case["endpoint"],
            "encode": case["encode"],
            "form": case["form"],
            "http": None,
            "elapsed_ms": None,
            "raw": None,
            "response": None,
        }
        with OUT.open("a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        return rec
    LAST_CALL = time.monotonic()
    try:
        body = r.json()
    except ValueError:
        body = None
    rec = {
        "ts": datetime.now(UTC).isoformat(),
        "tag": case["tag"],
        "method": case["method"],
        "endpoint": case["endpoint"],
        "encode": case["encode"],
        "form": case["form"],
        "http": r.status_code,
        "elapsed_ms": round(1000 * (LAST_CALL - t0), 1),
        "raw": r.text,
        "response": body,
    }
    with OUT.open("a") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return rec


def done_tags():
    """Tags already present in OUT (for resuming after a crash)."""
    if not OUT.exists():
        return set()
    return {json.loads(line)["tag"] for line in OUT.open()}


def run_batch(cases, skip_done=False):
    """Send cases serially, printing a one-line status per request."""
    done = done_tags() if skip_done else set()
    n = 0
    for case in cases:
        if case["tag"] in done:
            continue
        rec = send(case)
        n += 1
        resp = rec["response"]
        if isinstance(resp, dict):
            status = resp.get("status") or resp.get("error") or resp.get("ok")
        else:
            status = f"non-json:{(rec['raw'] or '')[:60]!r}"
        print(f"[{case['tag']}] http={rec['http']} -> {status}", flush=True)
    print(f"{n} requests appended to {OUT}")


# --- batch definitions -----------------------------------------------------


def atoms():
    """Valid atom inventory: x, powers (int/neg/rational/decimal), exp/log/sqrt,
    constants, coefficients, and input-side syntax variants."""
    ineqs = [
        ("atom:x", "x>0"),
        ("atom:x2", "x^2>0"),
        ("atom:x3", "x^3>0"),
        ("atom:xneg1", "x^-1>0"),
        ("atom:xneg2", "x^-2>0"),
        ("atom:x0", "x^0>0"),
        ("atom:x1", "x^1>0"),
        ("atom:xhalf", "x^(1/2)>0"),
        ("atom:xhalfdec", "x^0.5>0"),
        ("atom:x32", "x^(3/2)>0"),
        ("atom:x15dec", "x^1.5>0"),
        ("atom:x23", "x^(2/3)>0"),
        ("atom:x13", "x^(1/3)>0"),
        ("atom:xneg12", "x^(-1/2)>0"),
        ("atom:x24", "x^(2/4)>0"),
        ("atom:x20", "x^2.0>0"),
        ("atom:xsqrt2", "x^sqrt(2)>0"),
        ("atom:xneg1d2", "x^-1/2>0"),
        ("atom:exp", "exp(x)>0"),
        ("atom:log", "log(x)>0"),
        ("atom:sqrt", "sqrt(x)>0"),
        ("atom:ln", "ln(x)>0"),
        ("atom:estar", "e**x>0"),
        ("atom:ecaret", "e^x>0"),
        ("atom:sqrt4", "sqrt(4)>0"),
        ("atom:int", "2>0"),
        ("atom:frac", "1/2>0"),
        ("atom:dec", "0.5>0"),
        ("atom:coefexp", "2*exp(x)>0"),
        ("atom:x2r", "x*2>0"),
        ("atom:grp3", "(x+1)*3>0"),
        ("atom:xx", "x+x>0"),
        ("atom:1e3", "1e3*x>0"),
        ("atom:x1e2", "x^1e2>0"),
        ("atom:01coef", "0.1*x+1>0"),
        ("atom:longcoef", "0.123456789*x+1>0"),
        ("atom:tinycoef", "1e-10*x+1>0"),
        ("atom:hex", "0x10>x"),
        ("atom:underscore", "1_000>x"),
        ("atom:starstar", "x**2>0"),
        ("atom:caretspace", "x ^ 2>0"),
        ("atom:ws", "  x^2  >  0  "),
    ]
    return [prove(i, t) for t, i in ineqs]


def errors():
    """Malformed atoms and expressions: harvest every distinct error string."""
    ineqs = [
        ("err:foo", "foo(x)>0"),
        ("err:sin", "sin(x)>0"),
        ("err:tan", "tan(x)>0"),
        ("err:y", "x+y>0"),
        ("err:yalone", "y>0"),
        ("err:log2", "log(2)>0"),
        ("err:exp2x", "exp(2*x)>0"),
        ("err:expx1", "exp(x+1)>0"),
        ("err:expneg", "exp(-x)>0"),
        ("err:expX", "exp(X)>0"),
        ("err:EXP", "EXP(x)>0"),
        ("err:Log", "Log(x)>0"),
        ("err:loglog", "log(log(x))>0"),
        ("err:expexp", "exp(exp(x))>0"),
        ("err:sqrtx2", "sqrt(x^2)>0"),
        ("err:sqrt2x", "sqrt(2*x)>0"),
        ("err:abs", "abs(x)>0"),
        ("err:pipe", "|x|>0"),
        ("err:xdx", "x/x>0"),
        ("err:1dx", "1/x>0"),
        ("err:xd0", "x/0>0"),
        ("err:1d0", "1/0>x"),
        ("err:xx_pow", "x^x>0"),
        ("err:2dx", "2^x>0"),
        ("err:exppow", "exp(x)^2>0"),
        ("err:parenpow", "(x+1)^2>0"),
        ("err:negpow", "(-x)^2>0"),
        ("err:xmulx", "x*x>0"),
        ("err:xgrp", "x*(x+1)>0"),
        ("err:explogmul", "exp(x)*log(x)>0"),
        ("err:max", "max(x,1)>0"),
        ("err:log2arg", "log(x,2)>0"),
        ("err:log0arg", "log()>0"),
        ("err:2x", "2x>0"),
        ("err:xgtgt", "x>>0"),
        ("err:xgt", "x>"),
        ("err:gt0", ">0"),
        ("err:eqeq", "x==0"),
        ("err:neq", "x!=0"),
        ("err:chained", "x>0>1"),
        ("err:semi", "x>0;x>1"),
        ("err:unicode_ge", "x^2≥0"),
        ("err:unicode_minus", "x^2−1>0"),
        ("err:times", "2×x>0"),
        ("err:import", "import os"),
        ("err:E", "E^x>0"),
        ("err:eparen", "e^(x)>0"),
        ("err:ln_noparen", "ln x>0"),
        ("err:exp_noparen", "exp x>0"),
        ("err:dot", ".>0"),
        ("err:tuple", "x,x>0"),
        ("err:list", "[x]>0"),
        ("err:dict", "{x:1}>0"),
        ("err:bigint", "x^99999999999999999999>0"),
        ("err:infexp", "x^1e309>0"),
        ("err:nan", "nan>0"),
        ("err:div0pow", "x^(1/0)>0"),
        ("err:sqrtneg", "sqrt(-1)>0"),
        ("err:log-const-in-line", "x^2+1>0"),  # placeholder replaced below
    ]
    cases = [prove(i, t) for t, i in ineqs[:-1]]
    cases.append(prove("x^2+1>0", "err:line-ln-const", line="x+ln(2)"))
    return cases


def directions():
    """Inequality operators, side swaps, constant-only, x on both sides."""
    ineqs = [
        ("dir:gt", "x^2+1>2*x"),
        ("dir:ge", "x^2+1>=2*x"),
        ("dir:lt", "x^2+1<2*x"),
        ("dir:le", "x^2+1<=2*x"),
        ("dir:flip", "log(x)<exp(x)"),
        ("dir:flip2", "0<x"),
        ("dir:eq", "x^2+1=0"),
        ("dir:eqx", "x=0"),
        ("dir:xgtx", "x>x"),
        ("dir:xgex", "x>=x"),
        ("dir:constgt", "2>1"),
        ("dir:constlt", "1>2"),
        ("dir:zerogt", "0>0"),
        ("dir:zeroge", "0>=0"),
        ("dir:xsqgtx", "x^2>x"),
        ("dir:both", "exp(x)+x>log(x)+1"),
        ("dir:explog", "exp(x)>log(x)"),
        ("dir:explog261", "exp(x)-log(x)-261/112>0"),
    ]
    return [prove(i, t) for t, i in ineqs]


def normalization():
    """Term placement, ordering, coefficient formatting in normalized.*."""
    ineqs = [
        ("norm:order1", "1+x^2>0"),
        ("norm:order2", "x+exp(x)>0"),
        ("norm:order3", "exp(x)+x>0"),
        ("norm:order4", "x+x^2>0"),
        ("norm:order5", "x^2+x>0"),
        ("norm:negfirst", "-x+2>0"),
        ("norm:constpos", "x^2+1>log(x)"),
        ("norm:constneg", "x^2-1>log(x)"),
        ("norm:constright", "x^2>log(x)+1"),
        ("norm:constrightneg", "x^2>log(x)-1"),
        ("norm:cancel", "x^2-x^2+1>0"),
        ("norm:negsq", "-x^2+10>0"),
        ("norm:negexp", "-exp(x)+x+10>0"),
        ("norm:neglog", "-log(x)+x>0"),
        ("norm:negsqrt", "-sqrt(x)+x>0"),
        ("norm:negsqrt2", "-sqrt(x)+x+1>0"),
        ("norm:kitchen", "3*x^2+2*exp(x)-4*log(x)-5*sqrt(x)-x+7>0"),
        ("norm:dup", "exp(x)+exp(x)+x+x>0"),
        ("norm:coeffrac", "x^2+1/3*x>0"),
        ("norm:biggroup", "2*(x+exp(x))+3>0"),
    ]
    return [prove(i, t) for t, i in ineqs]


def curvature():
    """x^a convexity by exponent range; combo classification; EPS edges."""
    ineqs = [
        ("curv:neg2", "x^-2>0"),
        ("curv:neg01", "x^-0.1>0"),
        ("curv:01", "x^0.1>0"),
        ("curv:third", "x^(1/3)>0"),
        ("curv:099", "x^0.99>0"),
        ("curv:101", "x^1.01>0"),
        ("curv:eps1", "x^1.000000001>0"),
        ("curv:eps2", "x^1.0000000001>0"),
        ("curv:x10", "x^10>0"),
        ("curv:x100", "x^100>0"),
        ("curv:coeff-eps", "1e-9*exp(x)+1>0"),
        ("curv:coeff-eps2", "1.1e-9*exp(x)>0"),
        ("curv:x2exp", "x^2+exp(x)>0"),
        ("curv:x2log", "x^2+log(x)>0"),
        ("curv:x2sqrt", "x^2+sqrt(x)>0"),
        ("curv:explog", "exp(x)+log(x)>0"),
        ("curv:negsq", "-x^2+10>0"),
        ("curv:logexp", "log(x)-exp(x)>0"),
        ("curv:sqrtlog", "sqrt(x)-log(x)>0"),
        ("curv:twoconc", "-x^0.5-x^0.6+10>0"),
        ("curv:x2x3", "x^2-x^3+10>0"),
    ]
    return [prove(i, t) for t, i in ineqs]


def minimum():
    """minimum.* numerics: boundaries, doubling to 1e8, interior roots."""
    ineqs = [
        ("min:edge", "x>0"),
        ("min:edge2", "x+1>0"),
        ("min:1e8", "x^-1>0"),
        ("min:decr", "100>x"),
        ("min:decrconv", "x^-1+2-x>0"),
        ("min:interior", "x^2-4*x+5>0"),
        ("min:quarter", "x-sqrt(x)>0"),
        ("min:golden", "x^2-x-1>0"),
        ("min:irrat", "x^2-3*x+3>0"),
        ("min:underflow", "x^1000>0"),
        ("min:explinx", "exp(x)-x>0"),
    ]
    return [prove(i, t) for t, i in ineqs]


def tolerance():
    """f_min acceptance boundaries: strict >1e-8 vs non-strict >=-1e-8."""
    ineqs = [
        ("tol:eq0gt", "x^2-2*x+1>0"),
        ("tol:eq0ge", "x^2-2*x+1>=0"),
        ("tol:neg1e9", "x^2-2*x+0.999999999>0"),
        ("tol:neg1e9ge", "x^2-2*x+0.999999999>=0"),
        ("tol:neg1e8ge", "x^2-2*x+0.99999999>=0"),
        ("tol:neg11e8ge", "x^2-2*x+0.999999989>=0"),
        ("tol:pos5e9", "x^2-2*x+1.000000005>0"),
        ("tol:pos2e8", "x^2-2*x+1.00000002>0"),
        ("tol:pos1e8", "x^2-2*x+1.00000001>0"),
    ]
    return [prove(i, t) for t, i in ineqs]


def tangents():
    """Harvest tangent_at across parameter sweeps to pin the candidate set."""
    cases = []
    # article example: exp(x) - 2x - ln x > 1/sqrt(2)
    cases.append(prove("exp(x)-2*x-log(x)-1/sqrt(2)>0", "tan:article"))
    # sweep c in exp(x) - ln x > c ; cliff near f_min ~= 2.33
    for c in [
        "-2",
        "-1",
        "-0.5",
        "0",
        "0.5",
        "1",
        "1.5",
        "1.8",
        "2",
        "2.1",
        "2.2",
        "2.25",
        "2.28",
        "2.29",
        "2.3",
        "2.31",
        "2.32",
        "2.325",
        "2.33",
        "2.3301",
        "2.331",
        "2.332",
        "2.335",
        "2.34",
        "2.35",
        "2.4",
        "3",
    ]:
        cases.append(prove(f"exp(x)-log(x)-{c}>0", f"tan:explog:c={c}"))
    # x^2 + c > ln x ; argmin near 1/sqrt(2) ~= 0.7071
    for c in ["-0.84", "-0.8", "-0.5", "0", "0.5", "1"]:
        cases.append(prove(f"x^2-log(x)+{c}>0", f"tan:x2log:c={c}"))
    # x^2 + b > sqrt(x) ; argmin ~= 0.3969, candidate thresholds computable
    for b in ["0.6", "0.5", "0.479", "0.476", "0.473", "0.4726", "0.47248"]:
        cases.append(prove(f"x^2-sqrt(x)+{b}>0", f"tan:x2sqrt:b={b}"))
    # affine right side: tangent is the right side itself, tangent_at = convergent
    cases += [
        prove("x^2-2*x+2>0", "tan:affine:x0=1"),
        prove("x^2-4*x+5>0", "tan:affine:x0=2"),
        prove("x^2-10*x+40>0", "tan:affine:x0=5"),
        prove("x^2-10*x+40-log(x)>0", "tan:multiright"),
        prove("x^2-x^(2/3)+5>0", "tan:pow23"),
        prove("x^2-x^(1/3)+5>0", "tan:pow13"),
        prove("exp(x)+x^2-log(x)-sqrt(x)-5>0", "tan:bothmulti"),
        prove("exp(x)-log(x)-2>=0", "tan:ge-variant"),
        prove("exp(x)>=x+1", "tan:zero-gap"),
        prove("x+1-log(x)>0", "tan:left-affine"),
        prove("exp(x)-log(x)-5>0", "tan:failmin"),
        prove("x^2+1>0", "tan:noright"),
        prove("x^2>log(x)", "tan:plain"),
        prove("sqrt(x)*2>x", "tan:sqrt-coeff"),
        prove("x^3-x^2+1>0", "tan:cubic-left"),
    ]
    return cases


def provided_line():
    """line= param: valid/invalid lines, ok semantics, no effect on proof."""
    base = "exp(x)-log(x)-261/112>0"
    cases = [
        prove(base, "line:out-syntax", line="30/17*x-1+ln(17/30)"),
        prove(base, "line:float", line="1.7647058823529411*x-1.565105589269355"),
        prove(base, "line:simple", line="2*x-1"),
        prove("x^2+1>0", "line:x", line="x"),
        prove("x^2+1>0", "line:const", line="1"),
        prove("x^2+1>0", "line:bad", line="x+10"),
        prove("x^2+1>0", "line:nonaffine", line="x^2"),
        prove("x^2+1>0", "line:log", line="log(x)"),
        prove("x^2+1>0", "line:exp", line="exp(x)"),
        prove("x^2+1>0", "line:badparse", line="foo("),
        prove("x^2+1>0", "line:badatom", line="foo(x)"),
        prove("x^2+1>0", "line:empty", line=""),
        prove("x^2>x", "line:onfail", line="x+1"),
        prove("log(x)-x+10>0", "line:oninconcl", line="foo("),
        prove("x^2+1>0", "line:neg", line="-x+1"),
        prove("x^2+1>0", "line:frac", line="1/2*x+3/4"),
        prove("x^2+1>0", "line:group", line="2*(x+1)"),
        prove("exp(x)-2*x-log(x)-1/sqrt(2)>0", "line:article", line="41/14*(x-12/161)"),
        prove("x^2+1>0", "line:xx", line="x+x"),
        prove("x^2+1>0", "line:space", line="  x + 1  "),
    ]
    return cases


def http_edges():
    """Transport-level: methods, encodings, stray params, routing."""
    return [
        raw_case("http:get", "GET", "/convex/prove"),
        raw_case("http:get-params", "GET", "/convex/prove", form={"inequality": "x>0"}),
        raw_case("http:options", "OPTIONS", "/convex/prove"),
        raw_case("http:head", "HEAD", "/convex/prove"),
        raw_case("http:put", "PUT", "/convex/prove", form={"inequality": "x>0"}),
        raw_case("http:trailing", "POST", "/convex/prove/", form={"inequality": "x>0"}),
        raw_case("http:case", "POST", "/convex/Prove", form={"inequality": "x>0"}),
        raw_case("http:root-post", "POST", "/convex/", form={"inequality": "x>0"}),
        raw_case(
            "http:json", "POST", "/convex/prove", encode="json", json_body={"inequality": "x>0"}
        ),
        raw_case("http:textplain", "POST", "/convex/prove", encode="raw", body="inequality=x%3E0"),
        raw_case(
            "http:multipart",
            "POST",
            "/convex/prove",
            encode="multipart",
            form={"inequality": "x^2+1>0", "line": "x"},
        ),
        raw_case(
            "http:urlenc-line", "POST", "/convex/prove", form={"inequality": "x^2+1>0", "line": "x"}
        ),
        raw_case("http:line-only", "POST", "/convex/prove", form={"line": "x"}),
        raw_case("http:extra", "POST", "/convex/prove", form={"inequality": "x>0", "foo": "bar"}),
        raw_case(
            "http:domain", "POST", "/convex/prove", form={"inequality": "x^2>x", "domain": "0,10"}
        ),
        raw_case(
            "http:domain-bad", "POST", "/convex/prove", form={"inequality": "x>0", "domain": "abc"}
        ),
        raw_case(
            "http:no_line", "POST", "/convex/prove", form={"inequality": "x^2+1>0", "no_line": "1"}
        ),
        raw_case("http:lang", "POST", "/convex/prove", form={"inequality": "x>0", "lang": "en"}),
        raw_case("http:get-en", "GET", "/convex/en"),
    ]


def followup():
    """Round 2: domain/no_line params, denominator cutoff, cliff edge, misc."""
    cases = [
        prove("100>x", "dom:bind", extra={"domain": "0,10"}),
        prove("x^2>x", "dom:lobound", extra={"domain": "5,10"}),
        prove("x^2-sqrt(x)+0.5>0", "dom:hibound", extra={"domain": "0,0.4"}),
        prove("x^2-4*x+5>0", "dom:neglo", extra={"domain": "-5,inf"}),
        prove("x^2>x", "dom:hi1", extra={"domain": "0,1"}),
        prove("exp(x)-log(x)-2>0", "param:no_line", extra={"no_line": "1"}),
        prove("exp(x)-log(x)-2>0", "param:no-line", extra={"no-line": "1"}),
        prove("exp(x)-log(x)-2>0", "param:noline", extra={"noline": "1"}),
        prove("x^2+10000*x-log(x)>0", "cf:denom1e4"),
        prove("x^2+10002*x-log(x)>0", "cf:denom-over"),
        prove("exp(x)-log(x)-2.3303>0", "tan:cliff1"),
        prove("exp(x)-log(x)-2.33036>0", "tan:cliff2"),
        prove("exp(x)-log(x)-2.330366>0", "tan:cliff3"),
        prove("exp(x)-log(x)-2.330366124762>0", "tan:cliff4"),
        prove("0.000001*x+1>0", "disp:1e-6coef"),
        prove("x^2+1>0", "line:zero", line="0"),
        prove("x^2+1>=x", "tan:ge-affine"),
        prove("x^2-4*x+4>=0", "tan:ge-tangent2"),
        prove("foo(x)>0", "err:precedence", line="x^2"),
        prove("x^2+1>0", "line:ln", line="ln(2)*x"),
        prove("x^2-sqrt(x)+0.5>0", "dom:neglo-crash", extra={"domain": "-1,inf"}),
        prove("x^2+1>=x", "dom:neglo-ta0", extra={"domain": "-1,inf"}),
        prove("x^2-sqrt(x)+0.5>0", "dom:hi-xstar", extra={"domain": "0.5,inf"}),
        raw_case("http:query-post", "POST", "/convex/prove?inequality=x%3E0"),
        prove("x^2-2*x+2>=0", "tan:ge-affine2"),
        prove("x^2-2*x+1.00000001>=0", "tol:pos1e8ge"),
    ]
    return cases


def followup2():
    """Round 3: no-operator default, empty/whitespace input, = on convex,
    malformed domain halves, inverted domain, strict affine-right cert."""
    return [
        prove("x", "dir:noop"),
        prove("", "err:empty"),
        prove("   ", "err:blank"),
        prove("x<0", "dir:lt0"),
        prove("x^2-2*x+1=0", "dir:eqconvex"),
        prove("x>0", "dom:halfbad", extra={"domain": "0,"}),
        prove("x>0", "dom:badlo", extra={"domain": "a,b"}),
        prove("x^2>x", "dom:inverted", extra={"domain": "10,0"}),
        prove("x^2+1>x", "tan:gt-affine"),
        # 3001-term sum: collect_terms recursion depth ~3000 > 1000 -> RecursionError
        # inside prove(); reveals whether the site's except covers non-ValueError.
        prove("x+" * 3000 + "x", "err:recursion"),
    ]


BATCHES = {
    "followup": followup,
    "followup2": followup2,
    "atoms": atoms,
    "errors": errors,
    "directions": directions,
    "normalization": normalization,
    "curvature": curvature,
    "minimum": minimum,
    "tolerance": tolerance,
    "tangents": tangents,
    "provided_line": provided_line,
    "http_edges": http_edges,
}


def main():
    global DELAY
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("batches", nargs="*", choices=BATCHES)
    p.add_argument("--all", action="store_true", help="run every batch in order")
    p.add_argument("--list", action="store_true")
    p.add_argument("--delay", type=float, default=DELAY)
    p.add_argument(
        "--skip-done", action="store_true", help="skip tags already present in the output file"
    )
    args = p.parse_args()
    if args.list:
        total = 0
        for name, f in BATCHES.items():
            n = len(f())
            total += n
            print(f"{name}: {n} cases")
        print(f"total: {total}")
        return
    DELAY = args.delay
    names = list(BATCHES) if args.all else args.batches
    if not names:
        p.error("give batch names or --all")
    for name in names:
        run_batch(BATCHES[name](), skip_done=args.skip_done)


def followup3():
    """Length-guard threshold search + deep-recursion attempt under the cap."""
    cases = [
        prove("x" * 100 + ">0", "len:100"),
        prove("x" * 500 + ">0", "len:500"),
        prove("x" * 1000 + ">0", "len:1000"),
        prove("x" * 1500 + ">0", "len:1500"),
        prove("x" * 2000 + ">0", "len:2000"),
        # 1 char of recursion depth per '-' (UnaryOp USub); needs cap > ~1001
        prove("-" * 1100 + "x", "err:recursion2"),
    ]
    return cases


BATCHES["followup3"] = followup3


def followup4():
    """Pin down the length cap and test length/recursion on the line param."""
    return [
        prove("x" * 198 + ">0", "len:200"),
        prove("x" * 254 + ">0", "len:256"),
        prove("x" * 298 + ">0", "len:300"),
        prove("x" * 398 + ">0", "len:400"),
        prove("x^2+1>0", "line:long", line="x+" * 1500 + "x"),
        prove("x>0", "dom:long", extra={"domain": "0" * 3000}),
    ]


BATCHES["followup4"] = followup4


def followup5():
    """Pin the inequality length cap: total chars in (400, 502]."""
    return [
        prove("x" * 448 + ">0", "len:450"),
        prove("x" * 478 + ">0", "len:480"),
        prove("x" * 498 + ">0", "len:500"),
        prove("x" * 499 + ">0", "len:501"),
    ]


BATCHES["followup5"] = followup5


def followup6():
    """Exact cap: 482 chars ok, 501 rejected. Probe 500 and 490."""
    return [
        prove("x" * 498 + ">0", "len:cap500"),
        prove("x" * 488 + ">0", "len:cap490"),
        prove("x" * 499 + ">0", "len:cap501b"),
    ]


BATCHES["followup6"] = followup6


def followup7():
    """Check order of length-guard vs blank-guard, and line-length guard edge."""
    return [
        prove(" " * 501, "len:blank501"),
        prove(" " * 500, "len:blank500"),
        prove("x^2+1>0", "line:len600", line="x" * 600),
    ]


BATCHES["followup7"] = followup7


def followup8():
    """501 raw chars that strip to 500: raw-length check -> 输入过长."""
    return [prove("x" * 500 + " ", "len:raw501")]


BATCHES["followup8"] = followup8


if __name__ == "__main__":
    sys.exit(main())
