#!/usr/bin/env bash
# /demo 演示页的前端依赖：站端页面逐字节引用 jsdelivr/hertzen 的
# MathJax/KaTeX/html2canvas，本脚本把它们拉成 static/vendor/ 下的本地副本
# （该目录已 gitignore，属下载产物非源码）。
set -euo pipefail
cd "$(dirname "$0")/.."
VENDOR=src/attention_calculator/static/vendor
mkdir -p "$VENDOR/mathjax" "$VENDOR/katex/fonts"

get() {
  curl -sfL --max-time 60 -o "$VENDOR/$1" "$2" && echo "ok $1"
}

get mathjax/tex-mml-chtml.js https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js
get mathjax/tex-svg.js https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js
get katex/katex.min.css https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css
get katex/katex.min.js https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js
get html2canvas.min.js https://html2canvas.hertzen.com/dist/html2canvas.min.js

# KaTeX css 相对引用 fonts/*.woff2——按 css 实际引用清单全量拉取
grep -o 'fonts/[^)]*\.woff2' "$VENDOR/katex/katex.min.css" | sort -u | while read -r f; do
  get "katex/$f" "https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/$f"
done
