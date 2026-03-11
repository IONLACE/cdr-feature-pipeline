#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"


def run_step(script_name: str, pdb_file: str, extra_args: list[str] | None = None) -> None:
    script_path = SRC / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Missing script: {script_path}")

    cmd = [sys.executable, str(script_path), pdb_file]
    if extra_args:
        cmd.extend(extra_args)

    print("\n=== Running", script_name, "===")
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_merge_features() -> None:
    script_path = SRC / "11-merge_features.py"
    if not script_path.exists():
        raise FileNotFoundError(f"Missing script: {script_path}")

    cmd = [sys.executable, str(script_path)]
    print("\n=== Running 11-merge_features.py ===")
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_merge_xy(y_labels: str | None) -> None:
    script_path = SRC / "merge_X_Y.py"
    if not script_path.exists():
        raise FileNotFoundError(f"Missing script: {script_path}")

    cmd = [sys.executable, str(script_path)]
    if y_labels:
        cmd.append(y_labels)

    print("\n=== Running merge_X_Y.py ===")
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the zff feature pipeline for one complex in a single command. "
            "By default starts at step 02 (assumes MSA already exists)."
        )
    )
    parser.add_argument("pdb_file", help="Input PDB filename or absolute path")
    parser.add_argument(
        "--with-msa",
        action="store_true",
        help="Include step 01-build_msa.py before step 02 (experimental)",
    )
    parser.add_argument(
        "--skip-merge-features",
        action="store_true",
        help="Skip step 11-merge_features.py",
    )
    parser.add_argument(
        "--merge-xy",
        action="store_true",
        help="Run merge_X_Y.py after feature merge",
    )
    parser.add_argument(
        "--y-labels",
        type=str,
        default=None,
        help="Optional y-labels CSV path for merge_X_Y.py",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default=None,
        help="Optional explicit output CSV path for step 10-build_custom_nb_features.py",
    )
    args = parser.parse_args()

    try:
        if args.with_msa:
            run_step("01-build_msa.py", args.pdb_file)

        run_step("02-calc_conservation.py", args.pdb_file)
        run_step("03-get_CDRs.py", args.pdb_file)
        run_step("04-run_fpocket.py", args.pdb_file)
        run_step("05-get_all_contacts.py", args.pdb_file)
        run_step("06-run_pandaprot.py", args.pdb_file)
        run_step("07-interface_complimentarity.py", args.pdb_file)
        run_step("08-compute_sasa.py", args.pdb_file)
        run_step("09-CDR_peptide_features.py", args.pdb_file)

        step10_extra = ["--output_csv", args.output_csv] if args.output_csv else None
        run_step("10-build_custom_nb_features.py", args.pdb_file, extra_args=step10_extra)

        if not args.skip_merge_features:
            run_merge_features()

        if args.merge_xy:
            run_merge_xy(args.y_labels)

        print("\nPipeline completed successfully.")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"Pipeline failed (exit code {e.returncode})")


if __name__ == "__main__":
    main()
