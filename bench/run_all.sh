#!/usr/bin/env bash
# 全量评测（本地跑；本机即 devbox，旧的 rsync→ssh 远端流程已废）。
# 用法: bench/run_all.sh
# golden 数据由 harvest.py 采集；本脚本不重复采，缺数据时先跑 harvest。

set -eo pipefail
cd "$(git rev-parse --show-toplevel)"
PY=.venv/bin/python
mkdir -p bench/out

echo "=== parity ==="
$PY bench/parity.py bench/data/golden.jsonl --out bench/out/parity.jsonl | tee bench/out/parity_summary.txt

echo "=== decompose parity ==="
$PY bench/parity_decompose.py --out bench/out/parity_decompose.jsonl | tee bench/out/parity_decompose_summary.txt

echo "=== edge/fuzz replay ==="
$PY bench/replay_edge.py
$PY bench/replay_fuzz.py

echo "=== capture replay (fidelity + probes) ==="
$PY bench/replay_capture.py

echo "=== sibling parity (health + convex) ==="
$PY bench/parity_health.py
$PY bench/parity_convex.py

echo "=== verify (site 模式恒等式数值审计) ==="
$PY bench/verify.py bench/data/golden.jsonl --out bench/out/verify.jsonl | tee bench/out/verify_summary.txt

echo "=== judge_correct (mode=exact 正确性判官) ==="
$PY bench/judge_correct.py bench/data/golden.jsonl --adversarial -o bench/out/judge_correct.jsonl | tee bench/out/judge_correct_summary.txt

echo "=== report ==="
$PY bench/report.py --parity bench/out/parity.jsonl --verify bench/out/verify.jsonl --out bench/out/report.md
echo "done: bench/out/report.md"
