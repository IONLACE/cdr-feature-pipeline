from pathlib import Path

import pytest

from paths import Paths


def build_config() -> dict:
    return {
        "out_dir": "workspace",
        "inp_pdb": "input_pdbs",
        "meta_subdir": "meta",
        "features_subdir": "features",
        "logs_subdir": "logs",
        "paths": {
            "inp_pdb_dir": "${out_dir}/${inp_pdb}",
            "meta_dir": "${out_dir}/${meta_subdir}",
            "features_dir": "${out_dir}/${features_subdir}",
            "logs_dir": "${out_dir}/${logs_subdir}",
        },
    }


def test_paths_expands_placeholders_and_resolves_relative_paths() -> None:
    root = Path("/tmp/project")
    cfg = build_config()

    paths = Paths(root, cfg)

    assert paths.inp_pdb_dir == root / "workspace" / "input_pdbs"
    assert paths.meta_dir == root / "workspace" / "meta"
    assert paths.features_dir == root / "workspace" / "features"
    assert paths.logs_dir == root / "workspace" / "logs"


def test_paths_keeps_absolute_paths() -> None:
    root = Path("/tmp/project")
    cfg = build_config()
    cfg["paths"]["meta_dir"] = "/data/meta"

    paths = Paths(root, cfg)

    assert paths.meta_dir == Path("/data/meta")


def test_paths_raises_on_unknown_placeholder() -> None:
    root = Path("/tmp/project")
    cfg = build_config()
    cfg["paths"]["meta_dir"] = "${unknown_key}/meta"

    with pytest.raises(KeyError, match="Missing config key for placeholder"):
        Paths(root, cfg)


def test_paths_raises_when_logs_path_key_missing() -> None:
    root = Path("/tmp/project")
    cfg = build_config()
    cfg["paths"].pop("logs_dir")

    with pytest.raises(KeyError, match="Missing required config key"):
        Paths(root, cfg)
