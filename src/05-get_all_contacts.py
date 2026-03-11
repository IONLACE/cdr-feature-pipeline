import csv
import sys
import pandas as pd
import glob
import itertools
import numpy as np
from paths import PATHS
from Bio.PDB import NeighborSearch, PDBParser, Selection
from pathlib import Path

try:
    from pandaprot import PandaProt
except ImportError:
    PandaProt = None


def _parse_chain_ids(chains: str):
    if chains is None:
        return []

    chain_text = str(chains).strip()
    if not chain_text:
        return []

    parsed = list(chain_text)

    unique = []
    seen = set()
    for chain_id in parsed:
        if chain_id not in seen:
            unique.append(chain_id)
            seen.add(chain_id)
    return unique


def _extract_ids_from_stem(pdb_stem: str):
    parts = pdb_stem.split("_")
    if len(parts) < 3:
        raise ValueError(f"Unable to parse pdb_id/nb_chain/ag_chain from file name stem '{pdb_stem}'")

    pdb_id = parts[0]
    nb_chain = parts[1]

    if len(parts) > 3 and parts[-1].lower() in {"repair", "repaired"}:
        antigen_raw = "".join(parts[2:-1]) or parts[2]
    else:
        antigen_raw = parts[2]

    return pdb_id, nb_chain, antigen_raw


def load_cdr_mapping(cdr_csv, nb_chain):
    mapping = {}
    if not cdr_csv:
        return mapping
    with open(cdr_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['cdr_region'] in ('CDR1', 'CDR2', 'CDR3') and row.get('pdb_position'):
                pos = row['pdb_position']
                try:
                    resnum = int(''.join([c for c in pos if c.isdigit()]))
                except Exception:
                    continue
                mapping[(nb_chain, resnum)] = row['cdr_region']
    return mapping


def compute_all_contacts(pdb_file: str, nb_chain: str, ag_chain: str, cutoff: float = 5.0):
    # Write to CSV
    pdb_stem = Path(pdb_file).stem
    csv_dir = PATHS.meta_dir / f"{pdb_stem}" 
    cdr_csv = csv_dir / f"{Path(pdb_file).stem}_cdrs.csv"
    out_csv = csv_dir / "contacts" / f"{Path(pdb_file).stem}_all_and_CDR_contacts.csv"
    # make sure output directory exists
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    parser_pdb = PDBParser(QUIET=True)
    structure = parser_pdb.get_structure('complex', pdb_file)
    model = structure[0]

    ag_chain_list = _parse_chain_ids(ag_chain)
    if not ag_chain_list:
        raise ValueError("No antigen chain IDs were provided")

    chain_nb = model[nb_chain]
    missing_ag_chains = [chain for chain in ag_chain_list if chain not in model]
    if missing_ag_chains:
        raise ValueError(f"Antigen chain(s) {missing_ag_chains} not found in {pdb_file}")

    atoms_nb = [atom for atom in Selection.unfold_entities(chain_nb, 'A') if atom.element != 'H']
    atoms_ag = []
    for chain_id in ag_chain_list:
        chain_ag = model[chain_id]
        atoms_ag.extend([atom for atom in Selection.unfold_entities(chain_ag, 'A') if atom.element != 'H'])

    ns_ag = NeighborSearch(atoms_ag)

    cdr_map = load_cdr_mapping(cdr_csv, nb_chain)

    contacts = []
    for atom_nb in atoms_nb:
        res_nb = atom_nb.get_parent()
        if res_nb.id[0] != ' ':
            continue
        neighbors = ns_ag.search(atom_nb.coord, cutoff)
        for atom_ag in neighbors:
            res_ag = atom_ag.get_parent()
            if res_ag.id[0] != ' ':
                continue
            if res_nb.get_parent().id == res_ag.get_parent().id:
                continue
            nb_resnum = res_nb.id[1]
            ag_resnum = res_ag.id[1]
            cdr_region = cdr_map.get((nb_chain, nb_resnum), '')
            contacts.append({
                'antigen_chain': res_ag.get_parent().id,
                'antigen_resname': res_ag.get_resname(),
                'antigen_resnum': ag_resnum,
                'nanobody_chain': nb_chain,
                'nanobody_resname': res_nb.get_resname(),
                'nanobody_resnum': nb_resnum,
                'cdr_region': cdr_region,
            })

    fieldnames = [
        'antigen_chain',
        'antigen_resname',
        'antigen_resnum',
        'nanobody_chain',
        'nanobody_resname',
        'nanobody_resnum',
        'cdr_region',
    ]
    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in contacts:
            writer.writerow(row)
    print(f"All contacts computed and saved to {out_csv}")
    # return out_csv, len(contacts)
    return

def _load_cdr_region_map(cdr_csv: Path):
    if not cdr_csv.exists():
        return {}
    df = pd.read_csv(cdr_csv)
    if df.empty:
        return {}
    if "cdr_region" not in df.columns or "pdb_resid" not in df.columns:
        return {}

    out = {}
    for _, row in df.iterrows():
        try:
            rn = int(row["pdb_resid"])
        except Exception:
            continue
        out[rn] = str(row["cdr_region"]).strip().upper()
    return out


def build_all_and_cdr_contacts_csv(
    pdb_file: Path,
    nb_chain: str,
    ag_chain: str,
    cdr_csv: Path,
    output_csv: Path,
    cutoff: float = 5.0,
):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(pdb_file))
    model = structure[0]

    ag_chains = _parse_chain_ids(ag_chain)
    if not ag_chains:
        raise ValueError("No antigen chain IDs were provided")
    cdr_map = _load_cdr_region_map(cdr_csv)

    if nb_chain not in model:
        raise ValueError(f"Nanobody chain '{nb_chain}' not found in {pdb_file}")

    rows = []
    nb_obj = model[nb_chain]
    for ag in ag_chains:
        if ag not in model:
            continue
        ag_obj = model[ag]

        for nb_res, ag_res in itertools.product(nb_obj, ag_obj):
            if nb_res.id[0] != " " or ag_res.id[0] != " ":
                continue

            min_dist = min((a1 - a2 for a1 in nb_res for a2 in ag_res), default=np.inf)
            if not np.isfinite(min_dist) or min_dist > cutoff:
                continue

            nb_resid = int(nb_res.get_id()[1])
            ag_resid = int(ag_res.get_id()[1])
            rows.append(
                {
                    "nb_chain": nb_chain,
                    "nb_resnum": nb_resid,
                    "nb_resname": nb_res.get_resname(),
                    "antigen_chain": ag,
                    "antigen_resnum": ag_resid,
                    "antigen_resname": ag_res.get_resname(),
                    "distance": round(float(min_dist), 3),
                    "cdr_region": cdr_map.get(nb_resid, ""),
                }
            )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(rows)
    if not out_df.empty:
        out_df = out_df.drop_duplicates(
            subset=["nb_chain", "nb_resnum", "antigen_chain", "antigen_resnum", "cdr_region"]
        )
        out_df = out_df.sort_values(["antigen_chain", "antigen_resnum", "nb_resnum"])
    out_df.to_csv(output_csv, index=False)
    print(f"[INFO] Unified contacts+CDR mapping saved to {output_csv}")

def get_cdr_pocket_contacts(pdb_file: str):

    pdb_stem = Path(pdb_file).stem

    # /Users/satyan/Wrk/zScore01_stuff/meta/1kxq/contacts/1kxq_G_B_repaired_all_and_CDR_contacts.csv
    contacts_csv = PATHS.meta_dir/ f"{pdb_stem}" / "contacts" / f"{pdb_stem}_all_and_CDR_contacts.csv"
    if not contacts_csv.exists():
       raise FileNotFoundError(f"Contacts CSV not found: {contacts_csv}")

    # Load contacts CSV
    contacts = pd.read_csv(contacts_csv)
    # Only keep contacts where cdr_region is not empty
    contacts = contacts[contacts['cdr_region'].notnull() & (contacts['cdr_region'] != '')]
    
    pockets_pattern = str(PATHS.meta_dir / f"{pdb_stem}" / 'pocket' / f"{pdb_stem}_out" / 'pockets' / 'pocket*_atm.pdb')
    pocket_files = sorted(glob.glob(pockets_pattern))

    if not pocket_files:
        raise FileNotFoundError(f"No pocket files found for pattern: {pockets_pattern}")

    # Build a set of (antigen_chain, antigen_resnum) for each CDR
    cdr_to_ag_res = {}
    for cdr in ['CDR1', 'CDR2', 'CDR3']:
        subset = contacts[contacts['cdr_region'] == cdr]
        cdr_to_ag_res[cdr] = set(zip(subset['antigen_chain'], subset['antigen_resnum']))

    # For each pocket, check if it contains any antigen residues in contact with a CDR
    pocket_cdr_map = {pocket: {cdr: set() for cdr in cdr_to_ag_res} for pocket in pocket_files}
    pocket_features = {}
    for pocket_path in pocket_files:
        pocket_id = Path(pocket_path).name.split('_')[0]  # e.g., pocket1
        # Parse pocket residues
        residues = set()
        with open(pocket_path) as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain = line[21].strip()
                    resnum = int(line[22:26].strip())
                    residues.add((chain, resnum))
        # For each CDR, check if any of its antigen contacts are in this pocket
        for cdr, ag_res_set in cdr_to_ag_res.items():
            intersect = residues & ag_res_set
            if intersect:
                pocket_cdr_map[pocket_path][cdr] = intersect

    # Output: For each pocket, which CDRs are present, and pocket features
    rows = []
    for pocket_path in pocket_files:
        pocket_id = Path(pocket_path).name.split('_')[0]
        cdrs_present = [cdr for cdr, resset in pocket_cdr_map[pocket_path].items() if resset]
        if not cdrs_present:
            continue  # Only keep pockets with CDR contacts
        row = {'pocket_id': pocket_id, 'CDRs_in_pocket': ','.join(cdrs_present)}
        row.update(pocket_features.get(pocket_id, {}))
        for cdr in ['CDR1', 'CDR2', 'CDR3']:
            row[f'{cdr}_nres_in_pocket'] = len(pocket_cdr_map[pocket_path][cdr])
        rows.append(row)

    df = pd.DataFrame(rows)
    out_path = PATHS.meta_dir/ f"{pdb_stem}" / "contacts" / f"{pdb_stem}_CDR_pocket_contacts.csv"
    df.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    pdb_file = sys.argv[1]
    pdb_stem = Path(pdb_file).stem
    pdb_id, nb_ch, ag_ch = _extract_ids_from_stem(pdb_stem)

    if not Path(pdb_file).is_absolute():
        pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_file

    ag_chain_list = _parse_chain_ids(ag_ch)
    print(f"Input PDB file: {pdb_file}, nanobody chain: {nb_ch}, antigen chain(s): {','.join(ag_chain_list)}")        
  

    cdr_csv = PATHS.meta_dir / f"{pdb_stem}" / "cdrs" / f"{pdb_stem}_cdrs.csv"
    merged_contacts_csv = PATHS.meta_dir / f"{pdb_stem}" / "contacts" / f"{pdb_stem}_all_and_CDR_contacts.csv"
    build_all_and_cdr_contacts_csv(
        pdb_file=pdb_file,
        nb_chain=nb_ch,
        ag_chain=ag_ch,
        cdr_csv=cdr_csv,
        output_csv=merged_contacts_csv,
        cutoff=5.0,
    )
    
    get_cdr_pocket_contacts(pdb_file)

