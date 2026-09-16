"""Byte-level page parity: rendered routes vs captured site pages.

bench/data/site-*.html are the captured ground truth (gitignored corpus);
each case skips when its fixture is absent. /convex page parity lives in
test_convex.py.
"""

from pathlib import Path

import pytest

from attention_calculator import server

DATA = Path("bench/data")

# (route, site capture) — /en is byte-identical to / (client-side language);
# /health/ mirrors /health; /attention/* is the same app under a mount prefix.
PAGES = [
    ("/", "site-root.html"),
    ("/en", "site-root.html"),
    ("/attention", "site-attention.html"),
    ("/attention/", "site-attention_.html"),
    ("/attention/en", "site-attention_en.html"),
    ("/health", "site-health.html"),
    ("/health/", "site-health.html"),
    ("/health/en", "site-health-en.html"),
]


@pytest.fixture
def client():
    return server.app.test_client()


@pytest.mark.parametrize(
    ("route", "fixture"),
    [
        pytest.param(
            route,
            fixture,
            marks=pytest.mark.skipif(
                not (DATA / fixture).exists(), reason="bench/data 语料不入库（rsync 同步）"
            ),
            id=f"{route}->{fixture}",
        )
        for route, fixture in PAGES
    ],
)
def test_page_byte_parity(client, route, fixture):
    """Route renders the captured site page byte-for-byte."""
    expected = (DATA / fixture).read_bytes()
    resp = client.get(route)
    assert resp.status_code == 200
    assert resp.get_data() == expected


def test_demo_route(client):
    """GET /demo serves the index page with CDN URLs rewritten to /static/vendor."""
    resp = client.get("/demo")
    assert resp.status_code == 200
    body = resp.get_data()
    for path in (
        b"/static/vendor/mathjax/tex-mml-chtml.js",
        b"/static/vendor/katex/katex.min.css",
        b"/static/vendor/katex/katex.min.js",
        b"/static/vendor/html2canvas.min.js",
    ):
        assert path in body, path
    assert b"cdn.jsdelivr.net" not in body
    assert b"html2canvas.hertzen.com" not in body
