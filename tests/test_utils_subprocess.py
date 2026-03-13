from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from utils.utils_subprocess import run_command


def test_run_command_forwards_success_result(monkeypatch: pytest.MonkeyPatch) -> None:
    expected = subprocess.CompletedProcess(args=["echo", "ok"], returncode=0)

    def fake_run(
        cmd: list[str],
        check: bool,
        capture_output: bool,
        text: bool,
        cwd: str | None,
    ) -> subprocess.CompletedProcess:
        assert cmd == ["echo", "ok"]
        assert check is True
        assert capture_output is True
        assert text is True
        assert cwd == "/tmp"
        return expected

    monkeypatch.setattr("utils.utils_subprocess.run", fake_run)

    result = run_command(["echo", "ok"], cwd=Path("/tmp"))

    assert result is expected


def test_run_command_raises_system_exit_on_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run(*args, **kwargs):
        raise subprocess.CalledProcessError(
            returncode=2,
            cmd=["bad", "cmd"],
            output="some stdout",
            stderr="some stderr",
        )

    monkeypatch.setattr("utils.utils_subprocess.run", fake_run)

    with pytest.raises(SystemExit, match="Pipeline stopped due to external tool failure"):
        run_command(["bad", "cmd"])
