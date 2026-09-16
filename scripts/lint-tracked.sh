#!/usr/bin/env bash
# lint-tracked.sh — run the repo toolchain over git-tracked files only.
#
# Why not `pre-commit run -a`: pre-commit stashes all unstaged changes for
# the run, and on a shared worktree that stash window races peers' live
# edits (it ate a kernel-spec edit once). Scoping by `git ls-files` never
# touches untracked files and never stashes.
#
# Default is check-only; --fix writes formatting changes (still only to
# tracked files — peers' modified-but-uncommitted tracked files are still
# their business, coordinate before --fix on an active tree).
#
#   bash scripts/lint-tracked.sh          # check everything tracked
#   bash scripts/lint-tracked.sh --fix    # apply formatters
#   bash scripts/lint-tracked.sh py md    # restrict to some extensions
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || exit 1

FIX=0
[ "${1:-}" = "--fix" ] && {
  FIX=1
  shift
}
EXTS=("$@")
want() { [ ${#EXTS[@]} -eq 0 ] || printf '%s\n' "${EXTS[@]}" | grep -qx "$1"; }

PRETTIER=./node_modules/.bin/prettier
MDLINT=./node_modules/.bin/markdownlint-cli2
rc=0
run() { "$@" || rc=1; }

if want py; then
  files=$(git ls-files '*.py')
  if [ -n "$files" ]; then
    if [ $FIX -eq 1 ]; then run xargs ruff format <<<"$files"; else run xargs ruff format --check <<<"$files"; fi
    run xargs ruff check <<<"$files"
  fi
fi

if want sh; then
  files=$(git ls-files '*.sh')
  if [ -n "$files" ]; then
    if [ $FIX -eq 1 ]; then run xargs shfmt -i 2 -w <<<"$files"; else run xargs shfmt -i 2 -d <<<"$files"; fi
    for f in $files; do head -1 "$f" | grep -q zsh || shellcheck -S warning "$f" || rc=1; done
  fi
fi

if want md; then
  # docs/ is out of the chain — verbatim quotes are data, not prose.
  files=$(git ls-files '*.md' | grep -v '^docs/' || true)
  if [ -n "$files" ]; then
    if [ $FIX -eq 1 ]; then run "$MDLINT" --fix $files; else run "$MDLINT" $files; fi
    if [ $FIX -eq 1 ]; then run xargs "$PRETTIER" --write <<<"$files"; else run xargs "$PRETTIER" --check <<<"$files"; fi
  fi
fi

if want yaml; then
  files=$(git ls-files '*.yaml' '*.yml')
  if [ -n "$files" ]; then
    if [ $FIX -eq 1 ]; then run xargs "$PRETTIER" --write <<<"$files"; else run xargs "$PRETTIER" --check <<<"$files"; fi
  fi
fi

if want toml; then
  files=$(git ls-files '*.toml' | grep -v '^uv\.lock$' || true)
  if [ -n "$files" ]; then
    if [ $FIX -eq 1 ]; then run xargs taplo fmt <<<"$files"; else run xargs taplo fmt --check <<<"$files"; fi
  fi
fi

exit $rc
