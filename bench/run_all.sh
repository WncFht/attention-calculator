#!/usr/bin/env bash
# 全量评测：同步代码到 devbox → parity + verify → 聚合报告。
# 用法: bench/run_all.sh [devbox]
# golden 数据由 harvest.py 采集；本脚本不重复采，缺数据时先跑 harvest。

set -eo pipefail
REMOTE=${1:-devbox}
DIR=~/src/attention-calculator
PY=.venv/bin/python

rsync -az --exclude .venv --exclude __pycache__ --exclude .git ./ "$REMOTE:$DIR/"
ssh "$REMOTE" "cd $DIR && $PY -m pip list >/dev/null 2>&1 || true"

echo "=== parity ==="
ssh "$REMOTE" "cd $DIR && $PY bench/parity.py bench/data/golden.jsonl --out bench/out/parity.jsonl" | tee bench/out/parity_summary.txt

echo "=== verify ==="
ssh "$REMOTE" "cd $DIR && $PY bench/verify.py bench/data/golden.jsonl --out bench/out/verify.jsonl" | tee bench/out/verify_summary.txt

echo "=== report ==="
ssh "$REMOTE" "cd $DIR && $PY bench/report.py --parity bench/out/parity.jsonl --verify bench/out/verify.jsonl --out bench/out/report.md"
rsync -az "$REMOTE:$DIR/bench/out/" bench/out/
echo "done: bench/out/report.md"
