#!/usr/bin/env python3
"""离线复核 exact-mode 证明证书（docs/2026-09-16-certificate-spec.md）。

读取 certificate.build 产出的证书 JSON（文件参数或 stdin），从核函数
重算证明检查并打印判定。退出码：0 = 证书核验通过；1 = 证书格式完好但
核验失败（证明无效或被篡改）；2 = 输入不可读或证书格式错误。

    .venv/bin/python tools/verify_cert.py cert.json
    .venv/bin/python tools/verify_cert.py < cert.json
"""

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.certificate import cert_tex, recheck, verify_cert

# decode failures that mean "malformed" — the same except lists the family
# cert_parse calls sit under inside certificate.verify_cert
PARSE_ERRORS = (KeyError, TypeError, ValueError, ZeroDivisionError, AttributeError)


def second_prover_parse(cert: dict) -> dict:
    """Decode a second-prover cert's wire fields; raises on malformed input.

    Dispatch order mirrors certificate.verify_cert. composite has no
    cert_parse — its schema is the field extraction at the top of
    gamma_special.verify_cert, mirrored here.
    """
    if "serr" in cert:
        from attention_calculator import pade

        return pade.cert_parse(cert)
    if "agm_iter" in cert:
        from attention_calculator import agm

        return agm.cert_parse(cert)
    if cert.get("prover") == "euler_gamma":
        from attention_calculator import euler_gamma

        return euler_gamma.cert_parse(cert)
    return {
        "kind": cert["kind"],
        "comp": cert["comp"],
        "q": Fraction(cert["q"]),
        "p": Fraction(cert["p"]),
        "rule": cert["rule"],
        "witness": {k: Fraction(v) for k, v in cert["witness"].items()},
        "expect": [
            (e["kind"], Fraction(e["power"]), e["comp"], Fraction(e["bound"]))
            for e in cert["expect"]
        ],
        "children": cert["children"],
    }


def main() -> int:
    """CLI entry: parse args, load the cert JSON, recheck, print verdicts."""
    ap = argparse.ArgumentParser(description="verify an exact-mode proof certificate offline")
    ap.add_argument("cert", nargs="?", default="-", help="certificate JSON file ('-' = stdin)")
    args = ap.parse_args()
    try:
        text = sys.stdin.read() if args.cert == "-" else Path(args.cert).read_text()
        cert = json.loads(text)
    except (OSError, json.JSONDecodeError) as e:
        print(f"verify_cert: cannot read certificate: {e}", file=sys.stderr)
        return 2
    if isinstance(cert, dict) and (
        "serr" in cert or "agm_iter" in cert or cert.get("prover") in ("euler_gamma", "composite")
    ):
        # second-prover schemas (Padé/AGM/Euler–Maclaurin/composite):
        # certificate.verify_cert swallows parse failures into False, so the
        # family parse plus claim render run first — a decode failure is a
        # malformed cert (exit 2), not a disproven one
        try:
            second_prover_parse(cert)
            statement = cert_tex(cert)
        except PARSE_ERRORS as e:
            print(f"verify_cert: malformed certificate: {e}", file=sys.stderr)
            return 2
        print(f"statement:   {statement}")
        if verify_cert(cert):
            print("VERIFIED")
            return 0
        print("FAILED: certificate does not certify a valid proof", file=sys.stderr)
        return 1
    try:
        res = recheck(cert)
    except Exception as e:  # untrusted input: any checker/parse failure = malformed cert
        print(f"verify_cert: malformed certificate: {e}", file=sys.stderr)
        return 2
    print(f"identity_ok: {res['identity_ok']}")
    print(f"nonneg:      {res['nonneg']}")
    print(f"statement:   {cert_tex(cert)}")
    if verify_cert(cert):
        print("VERIFIED")
        return 0
    print("FAILED: certificate does not certify a valid proof", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
