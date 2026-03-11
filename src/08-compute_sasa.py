
import freesasa
import tempfile
import os
import csv
import re
import sys
import freesasa
from Bio.PDB import PDBParser
from pathlib import Path

from paths import PATHS


def calculate_residue_sasa(pdb_file: Path):
    structure = freesasa.Structure(str(pdb_file))
    result = freesasa.calc(structure)
    sasa = result.residueAreas()

    residue_sasa = {}
    for chain in sasa:
        for resnum in sasa[chain]:
            area = sasa[chain][resnum].total
            match = re.match(r"(-?\d+)", str(resnum).strip())
            if not match:
                continue
            residue_sasa[(chain, int(match.group(1)))] = area

    return residue_sasa


def get_residue_names(pdb_file: Path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_file))
    model = structure[0]

    residue_names = {}
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            residue_names[(chain.id, int(residue.id[1]))] = residue.get_resname().strip().upper()

    return residue_names


def write_subset_pdb(input_pdb, output_pdb, chains):
    chains_to_keep = set(chains)
    with open(input_pdb, "r") as src, open(output_pdb, "w") as dst:
        for line in src:
            if line.startswith(("ATOM", "HETATM", "ANISOU", "TER")):
                chain_id = line[21] if len(line) > 21 else " "
                if chain_id in chains_to_keep:
                    dst.write(line)
            elif line.startswith(("MODEL", "ENDMDL", "END")):
                dst.write(line)


def compute_dSASA(trimer_pdb: Path, nb_chain: str, ag_chains: str):
    trimer_pdb = Path(trimer_pdb)
    ag_chain_list = list(ag_chains)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        nb_pdb = tmpdir / "nb_only.pdb"
        ag_pdb = tmpdir / "ag_only.pdb"

        write_subset_pdb(trimer_pdb, nb_pdb, [nb_chain])
        write_subset_pdb(trimer_pdb, ag_pdb, ag_chain_list)

        sasa_complex = calculate_residue_sasa(trimer_pdb)
        sasa_nb = calculate_residue_sasa(nb_pdb)
        sasa_ag = calculate_residue_sasa(ag_pdb)

    residue_names = get_residue_names(trimer_pdb)
    rows = []

    for (chain, resid), sasa_unbound in sasa_nb.items():
        if (chain, resid) not in sasa_complex:
            continue
        rows.append(
            {
                "partner": "NB",
                "chain": chain,
                "resname": residue_names.get((chain, resid), "UNK"),
                "resid": resid,
                "sasa_unbound": round(sasa_unbound, 4),
                "sasa_complex": round(sasa_complex[(chain, resid)], 4),
                "dSASA": round(sasa_unbound - sasa_complex[(chain, resid)], 4),
            }
        )

    for (chain, resid), sasa_unbound in sasa_ag.items():
        if (chain, resid) not in sasa_complex:
            continue
        rows.append(
            {
                "partner": "AG",
                "chain": chain,
                "resname": residue_names.get((chain, resid), "UNK"),
                "resid": resid,
                "sasa_unbound": round(sasa_unbound, 4),
                "sasa_complex": round(sasa_complex[(chain, resid)], 4),
                "dSASA": round(sasa_unbound - sasa_complex[(chain, resid)], 4),
            }
        )

    rows.sort(key=lambda r: (r["partner"], r["chain"], r["resid"], r["resname"]))
    return rows


def write_dsasa_csv(rows, output_csv: Path):
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["partner", "chain", "resname", "resid", "sasa_unbound", "sasa_complex", "dSASA"],
        )
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    pdb_file = sys.argv[1]
    pdb_stem = Path(pdb_file).stem
    pdb_id = pdb_stem.split('_')[0]  # Extract PDB ID from filename
    nb_ch = pdb_stem.split('_')[1]  # Extract nanobody chain ID from filename
    ag_ch = pdb_stem.split('_')[2]  # Extract antigen chain ID from filename

    if not Path(pdb_file).is_absolute():
        pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_file

    print(f"Input PDB file: {pdb_file}, nanobody chain: {nb_ch}, antigen chain: {ag_ch}")        


    dsasa_rows = compute_dSASA(pdb_file, nb_ch, ag_ch)
    output_csv = PATHS.meta_dir / f"{Path(pdb_file).stem}" / "sasa" / "dsasa.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    write_dsasa_csv(dsasa_rows, output_csv)
    print(f"dSASA results written to {output_csv}")
    