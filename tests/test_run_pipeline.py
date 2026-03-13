from __future__ import annotations

from pathlib import Path

import pytest

import run_pipeline as rp


def test_run_step_builds_expected_command(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: list[list[str]] = []

    def fake_run(cmd: list[str], check: bool) -> None:
        assert check is True
        captured.append(cmd)

    monkeypatch.setattr(rp.subprocess, "run", fake_run)

    rp.run_step("02-calc_conservation.py", "sample.pdb", extra_args=["--flag", "value"])

    assert captured == [
        [
            rp.sys.executable,
            str(rp.SRC / "02-calc_conservation.py"),
            "sample.pdb",
            "--flag",
            "value",
        ]
    ]


def test_run_step_raises_when_script_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(rp, "SRC", tmp_path)

    with pytest.raises(FileNotFoundError, match="Missing script"):
        rp.run_step("missing.py", "sample.pdb")


def test_run_merge_xy_appends_y_labels_argument(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    captured: list[list[str]] = []
    merge_script = tmp_path / "merge_X_Y.py"
    merge_script.write_text("# placeholder\n", encoding="utf-8")

    def fake_run(cmd: list[str], check: bool) -> None:
        assert check is True
        captured.append(cmd)

    monkeypatch.setattr(rp, "SRC", tmp_path)
    monkeypatch.setattr(rp.subprocess, "run", fake_run)

    rp.run_merge_xy("labels.csv")

    assert captured == [[rp.sys.executable, str(merge_script), "labels.csv"]]


def test_main_runs_expected_steps_without_optional_parts(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    called_steps: list[tuple[str, list[str] | None]] = []
    merge_flags = {"features": False, "xy": False}

    def fake_run_step(
        script_name: str, pdb_file: str, extra_args: list[str] | None = None
    ) -> None:
        assert pdb_file == "input.pdb"
        called_steps.append((script_name, extra_args))

    def fake_merge_features() -> None:
        merge_flags["features"] = True

    def fake_merge_xy(_: str | None) -> None:
        merge_flags["xy"] = True

    monkeypatch.setattr(rp, "run_step", fake_run_step)
    monkeypatch.setattr(rp, "run_merge_features", fake_merge_features)
    monkeypatch.setattr(rp, "run_merge_xy", fake_merge_xy)
    monkeypatch.setattr(
        rp.sys,
        "argv",
        ["run_pipeline.py", "input.pdb", "--skip-merge-features"],
    )

    rp.main()

    assert [name for name, _ in called_steps] == [
        "02-calc_conservation.py",
        "03-get_CDRs.py",
        "04-run_fpocket.py",
        "05-get_all_contacts.py",
        "06-run_pandaprot.py",
        "07-interface_complimentarity.py",
        "08-compute_sasa.py",
        "09-CDR_peptide_features.py",
        "10-build_custom_nb_features.py",
    ]
    assert called_steps[-1][1] is None
    assert merge_flags["features"] is False
    assert merge_flags["xy"] is False
