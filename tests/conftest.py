"""Shared fixtures: mpmath working-precision isolation + the Flask client."""

import pytest
from mpmath import mp

from attention_calculator import server


@pytest.fixture(autouse=True)
def preserve_mp_dps():
    """Snapshot ``mp.dps`` before each test and restore it after, so a bare
    ``mp.dps = N`` in a test body cannot leak into later tests."""
    dps = mp.dps
    yield
    mp.dps = dps


@pytest.fixture
def client():
    """Flask test client over the real app (``server.app``)."""
    return server.app.test_client()
