#!/usr/bin/env bash
# dev-server.sh — local dev-server lifecycle: restart / stop / smoke.
#
# Why not `pkill -f attention_calculator.server`: -f matches the pattern
# anywhere in a command line, including the *calling* Bash tool's own
# wrapper process — it killed the shell mid-run twice in session history
# (the `serve[r]` bracket trick is fragile). Kill by listening port via ss
# instead.
#
#   bash scripts/dev-server.sh restart [PORT]   # kill :PORT listener, nohup restart, smoke /
#   bash scripts/dev-server.sh smoke  [PORT]    # curl / and POST /calculate sanity
#   bash scripts/dev-server.sh stop   [PORT]
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

CMD=${1:-restart}
PORT=${2:-8080}
LOG=/tmp/attn-demo.log
PY=.venv/bin/python

kill_port() {
  local pids
  pids=$(ss -tlnp "sport = :$PORT" 2>/dev/null | grep -oP 'pid=\K[0-9]+' | sort -u || true)
  [ -n "$pids" ] && kill $pids 2>/dev/null || true
  sleep 1
}

smoke() {
  local code
  code=$(curl -s -o /dev/null -w '%{http_code}' "http://127.0.0.1:$PORT/" || true)
  echo "GET / -> $code"
  curl -s -X POST "http://127.0.0.1:$PORT/calculate" \
    -d 'type=pi&power=1&comparison=%3C&rational=22/7' | head -c 400
  echo
}

case "$CMD" in
restart)
  kill_port
  nohup "$PY" -m attention_calculator.server >"$LOG" 2>&1 &
  echo "pid $! -> $LOG"
  sleep 2
  smoke
  ;;
smoke) smoke ;;
stop)
  kill_port
  echo "stopped :$PORT"
  ;;
*)
  echo "usage: $0 restart|smoke|stop [PORT]" >&2
  exit 2
  ;;
esac
