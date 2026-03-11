#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from anarci import anarci
from paths import PATHS
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import requests
from io import StringIO

def get_cdr_regions_imgt():
    return {
        "CDR1": list(range(27, 39)),
        "CDR2": list(range(56, 66)),
        "CDR3": list(range(105, 118)),
    }

def write_multiple_chains_rcsb(pdb_id: str, ag_chain: str, output_dir: Path):
    pdb_id = pdb_id.upper()
    chains = list(ag_chain.upper())  # "HL" -> ["H", "L"]
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(url)
    response.raise_for_status()
    
    fasta_io = StringIO(response.text)
    
    found = {chain: False for chain in chains}
    
    for record in SeqIO.parse(fasta_io, "fasta"):
        for chain in chains:
            if f"Chain {chain}" in record.description or f"Chains {chain}" in record.description:
                
                out_file = output_dir / f"ag_{chain}_rcsb.fasta"
                SeqIO.write(record, out_file, "fasta")
                found[chain] = True
    
    for chain, status in found.items():
        if not status:
            raise ValueError(f"Chain {chain} not found in PDB {pdb_id}")
    
    return [str(output_dir / f"ag_{c}_rcsb.fasta") for c in chains]


def read_sequence_from_fasta(fasta_path: Path) -> str:
    sequence_lines = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        seen_header = False
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seen_header = True
                if sequence_lines:
                    break
                continue
            sequence_lines.append(line)

    if not seen_header:
        raise ValueError(f"Input file is not FASTA (missing header '>'): {fasta_path}")

    sequence = "".join(sequence_lines).strip().upper()
    if not sequence:
        raise ValueError(f"No sequence found in FASTA: {fasta_path}")
    return sequence


def find_cdrs_from_sequence(sequence: str, cdr_definitions=None):
    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError("Sequence is empty")

    cdr_definitions = cdr_definitions or get_cdr_regions_imgt()
    numbering_results = anarci([("query", sequence)], scheme="imgt", output=False)

    if numbering_results is None or len(numbering_results) == 0:
        raise ValueError("ANARCI could not number this sequence")

    numbered, _, _ = numbering_results
    if numbered[0] is None:
        raise ValueError("No ANARCI numbering result")

    numbering, chain_type, species = numbered[0][0]

    cdrs = {name: [] for name in ["CDR1", "CDR2", "CDR3"]}
    for (imgt_pos, imgt_ins), residue in numbering:
        if residue == "-":
            continue

        for cdr_name, positions in cdr_definitions.items():
            if imgt_pos in positions:
                imgt_full = f"{imgt_pos}{imgt_ins}" if str(imgt_ins).strip() else f"{imgt_pos}"
                cdrs[cdr_name].append(
                    {
                        "cdr_region": cdr_name,
                        "imgt_position": imgt_full,
                        "residue": residue,
                    }
                )

    cdr_sequences = {
        cdr_name: "".join([item["residue"] for item in entries])
        for cdr_name, entries in cdrs.items()
    }

    return cdrs, cdr_sequences, chain_type, species, numbering


def extract_pdb_chain_residues(pdb_file: Path, chain_id: str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_file))
    model = structure[0]

    if chain_id not in model:
        available = ",".join(sorted([chain.id for chain in model]))
        raise ValueError(f"Chain {chain_id} not found in PDB. Available chains: {available}")

    chain = model[chain_id]
    residues = []
    for residue in chain:
        hetflag, resid, icode = residue.id
        if hetflag != " ":
            continue
        if "CA" not in residue:
            continue

        resname = residue.get_resname().strip().upper()
        one_letter = seq1(resname, undef_code="X")
        if one_letter == "X":
            continue

        icode_clean = str(icode).strip() or ""
        residues.append(
            {
                "residue": one_letter,
                "pdb_chain": chain_id,
                "pdb_resid": int(resid),
                "pdb_resname": resname,
            }
        )

    if not residues:
        raise ValueError(f"No protein residues found for chain {chain_id} in {pdb_file}")
    return residues


def build_imgt_to_pdb_map(numbering, pdb_chain_residues):
    imgt_to_pdb = {}
    seq_idx = 0

    for (imgt_pos, imgt_ins), residue in numbering:
        if residue == "-":
            continue
        if seq_idx >= len(pdb_chain_residues):
            break

        imgt_full = f"{imgt_pos}{imgt_ins}" if str(imgt_ins).strip() else f"{imgt_pos}"
        imgt_to_pdb[imgt_full] = pdb_chain_residues[seq_idx]
        seq_idx += 1

    return imgt_to_pdb


def write_cdr_csv(cdrs: dict, cdr_sequences: dict, output_csv: Path, imgt_to_pdb: dict = None):
    imgt_to_pdb = imgt_to_pdb or {}
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "cdr_region",
            "imgt_position",
            "residue",
            "cdr_sequence",
            "pdb_chain",
            "pdb_resid",
            "pdb_resname",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for cdr_name in ["CDR1", "CDR2", "CDR3"]:
            cdr_seq = cdr_sequences.get(cdr_name, "")
            for row in cdrs.get(cdr_name, []):
                pdb_row = imgt_to_pdb.get(row["imgt_position"], {})
                writer.writerow(
                    {
                        "cdr_region": cdr_name,
                        "imgt_position": row["imgt_position"],
                        "residue": row["residue"],
                        "cdr_sequence": cdr_seq,
                        "pdb_chain": pdb_row.get("pdb_chain", ""),
                        "pdb_resid": pdb_row.get("pdb_resid", ""),
                        "pdb_resname": pdb_row.get("pdb_resname", ""),
                    }
                )


def print_cdrs(cdr_sequences: dict):
    print("CDR sequences (IMGT):")
    for cdr_name in ["CDR1", "CDR2", "CDR3"]:
        cdr_seq = cdr_sequences.get(cdr_name, "")
        print(f"{cdr_name}: {cdr_seq} (len={len(cdr_seq)})")


def main():
    parser = argparse.ArgumentParser(
        description="Find CDRs using sequence fetched from RCSB by PDB filename pattern <pdbid>_<chain>_*.pdb"
    )
    parser.add_argument("pdb_file", help="PDB filename like 8q6k_C_BA_Repair.pdb")
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=None,
        help="Optional CSV output path (default: PATHS.meta_dir/<pdb_stem>/cdrs/<pdb_stem>_cdrs.csv)",
    )
    args = parser.parse_args()

    pdb_input = args.pdb_file
    pdb_stem = Path(pdb_input).stem
    stem_parts = pdb_stem.split("_")
    if len(stem_parts) < 2:
        raise ValueError("Expected pdb filename like <pdbid>_<chain>_*.pdb")

    pdb_id = stem_parts[0]
    nb_ch = stem_parts[1]

    pdb_file = Path(pdb_input)
    if not pdb_file.is_absolute():
        pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_input
    print(f"Input PDB file: {pdb_file}, chain ID: {nb_ch}")

    cdr_dir = PATHS.meta_dir / f"{pdb_stem}" / "cdrs"
    fasta_files = write_multiple_chains_rcsb(pdb_id, nb_ch, cdr_dir)
    sequence = read_sequence_from_fasta(Path(fasta_files[0]))

    cdrs, cdr_sequences, chain_type, species, numbering = find_cdrs_from_sequence(sequence)
    pdb_chain_residues = extract_pdb_chain_residues(pdb_file, nb_ch)
    imgt_to_pdb = build_imgt_to_pdb_map(numbering, pdb_chain_residues)
    print(f"Chain type: {chain_type}")
    print(f"Species: {species}")
    print_cdrs(cdr_sequences)

    output_csv = args.out_csv or (cdr_dir / f"{pdb_stem}_cdrs.csv")
    write_cdr_csv(cdrs, cdr_sequences, output_csv, imgt_to_pdb=imgt_to_pdb)
    print(f"CSV written to: {output_csv}")


if __name__ == "__main__":
    main()

