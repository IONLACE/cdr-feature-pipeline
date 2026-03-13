import argparse
import csv
import glob
import os
import shutil
import sys
from pathlib import Path
from paths import PATHS
from utils.utils_subprocess import run_command
import re
import json

DEFAULT_OUTPUT_DIR = PATHS.meta_dir 


def resolve_fpocket_bin() -> str:
    env_bin = os.environ.get("FPOCKET_BIN")
    if env_bin:
        return env_bin

    discovered = shutil.which("fpocket")
    if discovered:
        return discovered

    raise FileNotFoundError(
        "fpocket executable not found. Set FPOCKET_BIN or install fpocket in PATH."
    )

def run_fpocket(input_pdb: Path, nb_ch: str, output_dir: Path = DEFAULT_OUTPUT_DIR):
    input_pdb = Path(input_pdb).resolve()
    output_dir = Path(output_dir).resolve()

    if not input_pdb.exists():
        raise SystemExit(f"Input PDB not found: {input_pdb}")

    fpocket_bin = resolve_fpocket_bin()

    output_dir.mkdir(parents=True, exist_ok=True)
    staged_input = output_dir / input_pdb.name
    shutil.copy2(input_pdb, staged_input)

    cmd = [fpocket_bin, "-f", str(staged_input), "-c", str(nb_ch)]
    return run_command(cmd, cwd=output_dir)


def parse_fpocket_output(pocket_path):
    """
    Parse an fpocket pocket.txt file into a nested dictionary.

    Returns:
        {
            "Pocket 1": {"Score": 0.216, ...},
            "Pocket 2": {...},
            ...
        }
    """
    pockets = {}
    current_pocket = None

    pocket_path = Path(pocket_path)

    with open(pocket_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # Detect new pocket
            pocket_match = re.match(r"(Pocket\s+\d+)\s*:", line)
            if pocket_match:
                pocket_number = pocket_match.group(1)
                current_pocket = f"pocket{pocket_number.split()[1]}"
                pockets[current_pocket] = {}
                continue

            # Parse key-value lines
            if ":" in line and current_pocket is not None:
                key, value = line.split(":", 1)
                key = key.strip()
                value = float(value.strip())

                if key == "Hydrophobicity score" :
                    pockets[current_pocket]['hydrophob_score'] = value
                if key == "Mean local hydrophobic density" :
                    pockets[current_pocket]['hydrophob_density'] = value                
                if key == "Volume score" :
                    pockets[current_pocket]['vol_score'] = value                
                if key == "Cent. of mass - Alpha Sphere max dist" :
                    pockets[current_pocket]['com_asphere_max_dist'] = value
                if key == "Alpha sphere density" :
                    pockets[current_pocket]['asphere_density'] = value

                if key == "Apolar alpha sphere proportion" :
                    pockets[current_pocket]['Apolar_asphere_prop'] = value
                if key == "Polarity score" :
                    pockets[current_pocket]['pol_score'] = value
                if key == "Charge score" :
                    pockets[current_pocket]['charge_score'] = value
                if key == "Proportion of polar atoms" :
                    pockets[current_pocket]['prop_of_polar_atoms'] = value
                if key == "Polar SASA" :
                    pockets[current_pocket]['polar_SASA'] = value
 
    if pocket_path.name.endswith("_info.txt"):
        pdb_stem = pocket_path.stem.replace("_info", "")
        json_out = DEFAULT_OUTPUT_DIR / f"{pdb_stem}"
        json_out.mkdir(parents=True, exist_ok=True)
        with open(json_out / "pocket_features.json", "w") as f:
            json.dump(pockets, f, indent=4)
        print(f"Parsed pocket features saved to {json_out / 'pocket_features.json'}")

    return pockets


def extract_residue_atoms_from_pocket_pdb(pocket_file: Path):
    residue_to_atoms = {}
    with open(pocket_file, "r") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom_name = line[12:16].strip().upper()
            resname = line[17:20].strip().upper()
            chain = line[21].strip() or " "
            resid_raw = line[22:26].strip()
            if not resid_raw:
                continue

            try:
                resid = int(resid_raw)
            except ValueError:
                continue

            residue_key = (chain, resname, resid)
            residue_to_atoms.setdefault(residue_key, set()).add(atom_name)

    return residue_to_atoms


def write_merged_pocket_residues_csv(pdb_stem: str, output_dir: Path = DEFAULT_OUTPUT_DIR) -> Path:
    pockets_pattern = str(
        output_dir / f"{pdb_stem}" / "pocket" / f"{pdb_stem}_out" / "pockets" / "pocket*_atm.pdb"
    )
    pocket_files = sorted(glob.glob(pockets_pattern))

    pocket_residue_atoms = {}
    for pocket_file in pocket_files:
        pocket_name = Path(pocket_file).stem.split("_")[-2]
        residue_atoms = extract_residue_atoms_from_pocket_pdb(Path(pocket_file))
        for residue_key, atoms in residue_atoms.items():
            pocket_key = (residue_key[0], residue_key[1], residue_key[2], pocket_name)
            pocket_residue_atoms.setdefault(pocket_key, set()).update(atoms)

    merged_csv = output_dir / f"{pdb_stem}" / "pocket_residues.csv"
    merged_csv.parent.mkdir(parents=True, exist_ok=True)

    with open(merged_csv, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["chain", "resname", "resid", "atoms", "pocket"])
        writer.writeheader()

        for (chain, resname, resid, pocket_name), atoms in sorted(
            pocket_residue_atoms.items(),
            key=lambda item: (item[0][0], item[0][2], item[0][1], item[0][3]),
        ):
            writer.writerow(
                {
                    "chain": chain,
                    "resname": resname,
                    "resid": resid,
                    "atoms": "|".join(sorted(atoms)),
                    "pocket": pocket_name,
                }
            )

    return merged_csv



if __name__ == "__main__":
    # main()
# STDERR:
#  QH6047 qhull input error: use upper-Delaunay('Qu') or infinity-point('Qz') with Delaunay('d') or Voronoi('v')
    pdb_file = sys.argv[1]
    pdb_stem = Path(pdb_file).stem
    pdb_id = Path(pdb_file).stem.split('_')[0]  # Extract PDB ID from filename
    nb_ch = Path(pdb_file).stem.split('_')[1]

    if not Path(pdb_file).is_absolute():
        pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_file

    print(f"Input PDB file: {pdb_file}, chain ID dropped: {nb_ch}")       

    run_fpocket(pdb_file, nb_ch, DEFAULT_OUTPUT_DIR / f"{pdb_stem}" / "pocket")
    print(f"Pocket detection completed for {pdb_id}. Output in: {DEFAULT_OUTPUT_DIR / f'{pdb_stem}' / 'pocket'}")

    pocket_info_path = DEFAULT_OUTPUT_DIR / f"{pdb_stem}" / "pocket" /f"{pdb_stem}_out" / f"{pdb_stem}_info.txt"
    pocket_feats = parse_fpocket_output (str(pocket_info_path))

    merged_csv = write_merged_pocket_residues_csv(pdb_stem, DEFAULT_OUTPUT_DIR)
    print(f"Merged pocket residues written to {merged_csv}")