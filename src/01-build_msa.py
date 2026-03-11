import requests
from io import StringIO
import sys
from Bio import AlignIO, SeqIO, PDB
from pathlib import Path
from utils.utils_subprocess import run_command
from paths import PATHS

def write_multiple_chains_rcsb(pdb_id: str, ag_chain: str, output_dir: Path) -> None:
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

    for chain in chains:
        path_file = output_dir / f"ag_{chain}_rcsb_fasta_path.txt"
        with path_file.open("w", encoding="utf-8") as handle:
            handle.write(str(output_dir / f"ag_{chain}_rcsb.fasta") + "\n")


def filter_accessions(blast_tsv: Path, accessions_txt: Path, min_id: float, min_qcov: float) -> int:
    accessions = set()
    with blast_tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                accession = parts[0]
                pident = float(parts[1])
                qcovs = float(parts[2])
            except ValueError:
                continue

            if pident >= min_id and qcovs >= min_qcov:
                accessions.add(accession)

    with accessions_txt.open("w", encoding="utf-8") as out:
        for accession in sorted(accessions):
            out.write(f"{accession}\n")

    return len(accessions)


def prepend_query_to_homologs(query_fasta: Path, homologs_fasta: Path, out_fasta: Path) -> None:
    with out_fasta.open("w", encoding="utf-8") as out:
        out.write(query_fasta.read_text(encoding="utf-8").rstrip() + "\n")
        if homologs_fasta.exists() and homologs_fasta.stat().st_size > 0:
            out.write(homologs_fasta.read_text(encoding="utf-8").lstrip())

def build_msa_pipeline(
    query_fasta: Path,
    blastdb_dir: str,
    blast_tsv: Path,
    accessions_txt: Path,
    homologs_fasta: Path,
    homologs_with_query_fasta: Path,
    homologs_nr_fasta: Path,
    aligned_msa_fasta: Path,
    min_id: float,
    min_qcov: float,
    evalue: float,
    max_target_seqs: int,
    cdhit_identity: float,
) -> None:
    run_command(
        [
            "blastp",
            "-query",
            str(query_fasta),
            "-db",
            blastdb_dir,
            "-outfmt",
            "6 sacc pident qcovs evalue bitscore",
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            str(max_target_seqs),
            "-out",
            str(blast_tsv),
        ]
    )

    accession_count = filter_accessions(blast_tsv, accessions_txt, min_id=min_id, min_qcov=min_qcov)
    if accession_count == 0:
        raise ValueError("No accessions passed filtering thresholds. Try lowering --min-id or --min-qcov.")

    run_command(
        [
            "blastdbcmd",
            "-db",
            blastdb_dir,
            "-entry_batch",
            str(accessions_txt),
            "-out",
            str(homologs_fasta),
        ]
    )

    run_command(
        [
            "cd-hit",
            "-i",
            str(homologs_fasta),
            "-o",
            str(homologs_nr_fasta),
            "-c",
            str(cdhit_identity),
        ]
    )
    
    prepend_query_to_homologs(query_fasta, homologs_nr_fasta, homologs_with_query_fasta)
    run_command(
        [
            "clustalo",
            "-i",
            str(homologs_with_query_fasta),
            "-o",
            str(aligned_msa_fasta),
            "--force",
            "--outfmt=fasta",
        ]
    )
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

    write_multiple_chains_rcsb(pdb_id, ag_ch, msa_dir)


    for chain in ag_ch:
        query_fasta = msa_dir / f"ag_{chain}_rcsb.fasta"
        build_msa_pipeline(
            query_fasta=query_fasta,
            blastdb_dir=str(PATHS.blastdb_dir),
            blast_tsv=msa_dir / f"ag_{chain}_blast.tsv",
            accessions_txt=msa_dir / f"ag_{chain}_accessions.txt",
            homologs_fasta=msa_dir / f"ag_{chain}_homologs.fasta",
            homologs_with_query_fasta=msa_dir / f"ag_{chain}_homologs_with_query.fasta",
            homologs_nr_fasta=msa_dir / f"ag_{chain}_homologs_nr.fasta",
            aligned_msa_fasta=msa_dir / f"ag_{chain}_aligned_msa.fasta",
            min_id=30.0,
            min_qcov=50.0,
            evalue=1e-5,
            max_target_seqs=500,
            cdhit_identity=0.9,
        )
