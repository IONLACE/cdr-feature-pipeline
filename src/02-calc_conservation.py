import sys
import csv
from Bio import AlignIO, SeqIO, PDB
from pathlib import Path
from paths import PATHS
from utils.utils_subprocess import run_command

from io import StringIO


def clean_msa(msa_fasta, clean_msa_fasta):
    """Clean the MSA by removing columns with gaps in the query sequence."""
    alignment = AlignIO.read(msa_fasta, "fasta")

    # 1) Clean sequences first
    cleaned_seqs = [record.seq for record in alignment]

    # 2) Determine columns to keep based on query (first sequence)
    query = cleaned_seqs[0]
    keep_cols = [i for i, aa in enumerate(query) if aa != "-"]

    # 3) Build new alignment
    new_records = []
    for record, seq in zip(alignment, cleaned_seqs):
        new_seq = "".join(seq[i] for i in keep_cols)
        record.seq = type(record.seq)(new_seq)
        new_records.append(record)

    # 4) Write cleaned MSA
    SeqIO.write(new_records, clean_msa_fasta, "fasta")

def msa_to_cons_scores(clean_msa_fasta: str, cons_scores_file: str):
    cmd = ["/Users/satyan/Wrk/packages/cons-capra07-main/target/release/cons-capra07", 
    "-i", str(clean_msa_fasta), 
    "-o", str(cons_scores_file), 
    "-w", "va", 
    "-c", "yes", 
    "-t", "yes"]

    result = run_command(cmd)

def write_scores_to_csv(cons_scores_file: str, chain_id: str, csv_out: str = None):
    one_to_three = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }

    scores = {}
    rows = []
    # csv_out = str(Path(cons_scores_file).with_suffix(".csv"))

    with open(cons_scores_file, "r", encoding="utf-8") as f:
        next(f, None)
        for i, line in enumerate(f, 1):
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue

            try:
                resid = int(parts[0])
                score = float(parts[1])
            except ValueError:
                print(f"Skipping invalid line {i}: {line.strip()}")
                continue

            seq_col = parts[2].strip()
            if not seq_col:
                continue

            one_letter = seq_col[0].upper()
            resname = one_to_three.get(one_letter, "UNK")

            scores[resid] = score
            rows.append({
                "chain": chain_id,
                "resname": resname,
                "resid": resid,
                "score": score,
            })

            # print(f"{resname}\t{resid}\t{score}")

    print(f"Parsed {len(scores)} conservation scores")

    if csv_out:
        with open(csv_out, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=["chain", "resname", "resid", "score"])
            writer.writeheader()
            writer.writerows(rows)
        print(f"Residue-wise scores written to {csv_out}")

    return rows


def write_combined_scores_csv(rows, csv_out: str):
    with open(csv_out, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["chain", "resname", "resid", "score"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Combined residue-wise scores written to {csv_out}")

if __name__ == "__main__":

    pdb_file = sys.argv[1] #"8q6k_C_BA_Repair.pdb"
    pdb_stem = Path(pdb_file).stem
    pdb_id = pdb_stem.split('_')[0]  # Extract PDB ID from filename
    ag_ch = pdb_stem.split('_')[2]
    if not Path(pdb_file).is_absolute():
        pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_file
    else:
        pdb_file = pdb_file
    print(f"Input PDB file: {pdb_file}, Antigen chain ID: {ag_ch}")

    msa_dir = PATHS.meta_dir / f"{pdb_stem}" / "msa"
    msa_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    for chain in ag_ch:
        ncbi_aln = msa_dir / f"ag_{chain}_ncbi.fa"
        clean_aln = msa_dir / f"ag_{chain}_ncbi_clean.fasta"

        if not ncbi_aln.exists():
            raise FileNotFoundError(f"Missing input MSA for chain {chain}: {ncbi_aln}")

        clean_msa(msa_fasta=ncbi_aln, clean_msa_fasta=clean_aln)

        cons_scores_file = msa_dir / f"ag_{chain}_cons_scores.tsv"
        msa_to_cons_scores(clean_msa_fasta=clean_aln, cons_scores_file=cons_scores_file)  

        chain_rows = write_scores_to_csv(cons_scores_file, chain_id=chain)
        all_rows.extend(chain_rows)

        # output_pdb = msa_dir / f"ag_{chain}_cons_scores.pdb"
        # map_scores_to_pdb(pdb_file=pdb_file, ag_ch=chain, cons_scores_file=cons_scores_file, output_pdb=output_pdb)

        

    combined_csv_out = PATHS.meta_dir / f"{pdb_stem}" / f"{pdb_stem}_cons_mapped.csv"
    write_combined_scores_csv(all_rows, str(combined_csv_out))

  