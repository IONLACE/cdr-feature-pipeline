from pathlib import Path
import csv
import shutil
import subprocess
import sys

from paths import PATHS


def extract_chains_from_filename(pdb_file: str):
    stem = Path(pdb_file).stem
    parts = stem.split("_")
    if len(parts) < 3:
        raise ValueError(f"Expected filename like <pdbid>_<nb_chain>_<ag_chains>[_suffix].pdb, got: {Path(pdb_file).name}")

    nb_chain = parts[1]
    ag_chains = parts[2]
    chain_ids = list(ag_chains) + [nb_chain]
    return chain_ids


def run_pandaprot_in_logdir(input_pdb: str, chain_ids):
    """Copy input PDB to a per-structure PandaProt log directory and run PandaProt there.

    Args:
        input_pdb: Path to input PDB file.
        chain_ids: Chain IDs passed to ``--chains``.

    Returns:
        dict with run directory, copied pdb path, command, and log file paths.
    """
    pdb_path = Path(input_pdb)
    if not pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {pdb_path}")

    chain_ids = [str(chain).strip() for chain in chain_ids if str(chain).strip()]
    if not chain_ids:
        raise ValueError("At least one chain ID is required")

    run_dir = PATHS.meta_dir / Path(input_pdb).stem / "pandaprot_run"
    run_dir.mkdir(parents=True, exist_ok=True)

    copied_pdb = run_dir / pdb_path.name
    shutil.copy2(pdb_path, copied_pdb)

    cmd = [
        "pandaprot",
        copied_pdb.name,
        "--chains",
        *chain_ids,
        "--statistics",
        "--residue-summary",
        "--export-vis",
        "--report",
    ]

    stdout_log = run_dir / "pandaprot_stdout.log"
    stderr_log = run_dir / "pandaprot_stderr.log"
    with open(stdout_log, "w") as out_f, open(stderr_log, "w") as err_f:
        subprocess.run(cmd, cwd=run_dir, stdout=out_f, stderr=err_f, check=True)

    return {
        "run_dir": run_dir,
        "copied_pdb": copied_pdb,
        "command": " ".join(cmd),
        "stdout_log": stdout_log,
        "stderr_log": stderr_log,
    }


def _row_atom_endpoint(row, chain_key, residue_key, atom_key):
    return (
        str(row.get(chain_key, "")).strip(),
        str(row.get(residue_key, "")).strip(),
        str(row.get(atom_key, "")).strip(),
    )


def filter_pandaprot_report(report_csv: Path, nb_chain: str, antigen_chains):
    if not report_csv.exists():
        raise FileNotFoundError(f"PandaProt report not found: {report_csv}")

    antigen_set = {str(ch).strip() for ch in antigen_chains if str(ch).strip()}
    nb_chain = str(nb_chain).strip()

    with report_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    kept_rows = []
    seen = set()
    for row in rows:
        chain1 = str(row.get("Chain1", "")).strip()
        chain2 = str(row.get("Chain2", "")).strip()

        is_nb_ag = (
            (chain1 == nb_chain and chain2 in antigen_set)
            or (chain2 == nb_chain and chain1 in antigen_set)
        )
        if not is_nb_ag:
            continue

        endpoint1 = _row_atom_endpoint(row, "Chain1", "Residue1", "Atom1")
        endpoint2 = _row_atom_endpoint(row, "Chain2", "Residue2", "Atom2")
        canonical_pair = tuple(sorted((endpoint1, endpoint2)))
        dedup_key = (str(row.get("Interaction_Type", "")).strip(), canonical_pair)

        if dedup_key in seen:
            continue
        seen.add(dedup_key)
        kept_rows.append(row)

    with report_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(kept_rows)

    return len(rows), len(kept_rows)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python -m src.06-run_pandaprot <pdb_file.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    chain_ids = extract_chains_from_filename(pdb_file)
    pdb_stem = Path(pdb_file).stem
    parts = pdb_stem.split("_")
    nb_chain = parts[1]
    antigen_chains = list(parts[2])

    pdb_id = Path(pdb_file).stem.split("_")[0]

    if not Path(pdb_file).is_absolute():
        pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_file

    print(f"Input PDB file: {pdb_file}, chain IDs: {', '.join(chain_ids)}")

    result = run_pandaprot_in_logdir(pdb_file, chain_ids)
    report_csv = result["run_dir"] / "pandaprot_report.csv"
    before_n, after_n = filter_pandaprot_report(report_csv, nb_chain=nb_chain, antigen_chains=antigen_chains)

    print(f"PandaProt run completed. Logs and outputs are in: {result['run_dir']}")
    print(f"Filtered report saved to {report_csv} ({before_n} -> {after_n} rows)")
