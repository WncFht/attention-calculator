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
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.certificate import cert_tex, recheck, verify_cert


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
    if isinstance(cert, dict) and "serr" in cert:
        # Padé second-prover certificate: pade.verify_cert IS the recheck —
        # it re-derives the approximants and error data over QQ
        ok = verify_cert(cert)
        print(f"identity_ok: {ok}")
        print(f"nonneg:      {ok}")
        print(f"statement:   pade {cert.get('kind')} {cert.get('comp')} {cert.get('p')}")
        if ok:
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
