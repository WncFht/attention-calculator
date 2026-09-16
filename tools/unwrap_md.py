#!/usr/bin/env python3
"""Unwrap hard-wrapped markdown paragraphs: one paragraph per line.

Joins soft line breaks inside paragraphs, list-item continuations, and
blockquote continuations. Never touches fenced code blocks, table rows,
headings, HTML lines, setext/thematic-break markers, or a line boundary
that falls inside an open `code span` (odd backtick count -> keep newline,
protects verbatim evidence strings that span lines).
"""

import re
import sys
from pathlib import Path

FENCE = re.compile(r"^\s*(```|~~~)")
HEAD = re.compile(r"^\s*#{1,6}\s")
TABLE = re.compile(r"^\s*\|")
HTML = re.compile(r"^\s*<")
SETEXT = re.compile(r"^\s*(={2,}|-{3,}|\*{3,}|_{3,})\s*$")
ITEM = re.compile(r"^\s*(?:[-*+]|\d+[.)])\s")
QUOTE = re.compile(r"^\s*>\s?")


def unwrap(text: str) -> str:
    out: list[str] = []
    buf: str | None = None
    btype = ""
    fence = False

    def flush():
        nonlocal buf
        if buf is not None:
            out.append(buf)
            buf = None

    for line in text.split("\n"):
        if FENCE.match(line):
            flush()
            out.append(line)
            fence = not fence
            continue
        if fence:
            out.append(line)
            continue
        s = line.strip()
        if not s:
            flush()
            out.append("")
            continue
        if HEAD.match(line) or TABLE.match(line) or HTML.match(line) or SETEXT.match(line):
            flush()
            out.append(line)
            continue
        if QUOTE.match(line):
            if btype == "quote" and buf is not None and buf.count("`") % 2 == 0:
                buf += " " + QUOTE.sub("", line).strip()
            else:
                flush()
                buf = line.lstrip()
                btype = "quote"
            continue
        if ITEM.match(line):
            flush()
            buf = line.rstrip()
            btype = "item"
            continue
        indented_code = line.startswith("    ") and btype == "para"
        if buf is not None and buf.count("`") % 2 == 0 and not indented_code:
            buf += " " + s
        else:
            flush()
            buf = line.rstrip()
            btype = "para"
    flush()
    return "\n".join(out)


def main():
    for arg in sys.argv[1:]:
        p = Path(arg)
        old = p.read_text()
        new = unwrap(old)
        if new == old:
            continue
        p.write_text(new)
        print(f"unwrapped {p}")


if __name__ == "__main__":
    main()
