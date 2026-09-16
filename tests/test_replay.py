"""pytest wrappers for the offline replay judges (replay_edge / replay_fuzz /
replay_capture).

The judge scripts refresh ours_* columns in place; tests redirect their path
constants to tmp copies so the gitignored corpus is not churned by every run.
Each test skips when its dataset is absent (CI checkouts lack bench/data).
"""

import importlib.util
import re
import shutil
from pathlib import Path

import pytest

from attention_calculator import server

BENCH = Path(__file__).resolve().parent.parent / "bench"
DATA = BENCH / "data"


def load_judge(name):
    """Import a bench judge script as a module."""
    spec = importlib.util.spec_from_file_location(name, BENCH / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def assert_zero_mismatch(output: str):
    """Every summary line the judges print must report mismatch=0."""
    counts = re.findall(r"mismatch=(\d+)", output)
    assert counts, "judge printed no summary"
    assert all(c == "0" for c in counts), output[-2000:]


@pytest.mark.skipif(
    not (DATA / "edge-probes.jsonl").exists(), reason="bench/data 语料不入库（rsync 同步）"
)
def test_replay_edge(tmp_path, capsys):
    mod = load_judge("replay_edge")
    mod.OUT = tmp_path / "edge-probes.jsonl"
    shutil.copy(DATA / "edge-probes.jsonl", mod.OUT)
    mod.main()
    assert_zero_mismatch(capsys.readouterr().out)


@pytest.mark.skipif(
    not (DATA / "fuzz-probes.jsonl").exists(), reason="bench/data 语料不入库（rsync 同步）"
)
def test_replay_fuzz(tmp_path, capsys):
    mod = load_judge("replay_fuzz")
    mod.OUT = tmp_path / "fuzz-probes.jsonl"
    shutil.copy(DATA / "fuzz-probes.jsonl", mod.OUT)
    mod.main()
    assert_zero_mismatch(capsys.readouterr().out)


@pytest.mark.skipif(
    not (DATA / "fidelity-probes.jsonl").exists() or not (DATA / "probes.jsonl").exists(),
    reason="bench/data 语料不入库（rsync 同步）",
)
def test_replay_capture(tmp_path, capsys):
    mod = load_judge("replay_capture")
    for name in ("fidelity-probes.jsonl", "probes.jsonl"):
        shutil.copy(DATA / name, tmp_path / name)
    client = server.app.test_client()
    mod.replay_edge_schema(client, tmp_path / "fidelity-probes.jsonl")
    mod.replay_probes_schema(client, tmp_path / "probes.jsonl")
    assert_zero_mismatch(capsys.readouterr().out)
