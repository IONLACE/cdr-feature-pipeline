import sys
from pathlib import Path
import peptides
import pandas as pd
from propy import PyPro
import numpy as np

from paths import PATHS, Paths

# Chou–Fasman propensities
HELIX = {
    'A':1.45,'R':1.00,'N':0.67,'D':1.01,'C':0.77,'Q':1.11,'E':1.51,'G':0.57,
    'H':1.00,'I':1.08,'L':1.34,'K':1.07,'M':1.20,'F':1.12,'P':0.59,'S':0.79,
    'T':0.82,'W':1.14,'Y':0.61,'V':1.06
}

SHEET = {
    'A':0.97,'R':0.90,'N':0.89,'D':0.54,'C':1.30,'Q':1.10,'E':0.37,'G':0.75,
    'H':0.87,'I':1.60,'L':1.22,'K':0.74,'M':1.05,'F':1.38,'P':0.55,'S':0.75,
    'T':1.19,'W':1.37,'Y':1.47,'V':1.70
}

TURN = {
    'A':0.66,'R':0.95,'N':1.56,'D':1.46,'C':1.19,'Q':0.98,'E':0.74,'G':1.56,
    'H':0.95,'I':0.47,'L':0.59,'K':1.01,'M':0.60,'F':0.60,'P':1.52,'S':1.43,
    'T':0.96,'W':0.96,'Y':1.14,'V':0.50
}

KIDERA = {
'A':[ -1.56,-1.67,-0.97,-0.27,-0.93,-0.78,-0.20,-0.08, 0.21,-0.48],
'R':[ 0.22, 1.27, 1.37, 1.87,-1.70, 0.46, 0.92,-0.39, 0.23, 0.93],
'N':[ 1.14,-0.07,-0.12, 0.81, 0.18, 0.37,-0.09, 1.23, 1.10,-1.73],
'D':[ 0.58,-0.22,-1.58, 0.81,-0.92, 0.15,-1.52, 0.47, 0.76, 0.70],
'C':[ 0.12,-0.89, 0.45,-1.05,-0.71, 2.41, 1.52,-0.69, 1.13, 1.10],
'Q':[ -0.47, 0.24, 0.07, 1.10, 1.10, 0.59, 0.84,-0.71,-0.03,-2.33],
'E':[ -1.45, 0.19,-1.61, 1.17,-1.31,-0.15,-1.05, 0.15, 0.20, 0.50],
'G':[ 1.46,-1.96,-0.23,-0.16, 0.10,-0.11, 1.32, 2.36,-1.66, 0.46],
'H':[ -0.41, 0.52,-0.28, 0.28, 1.61, 1.01,-1.85, 0.47, 1.13, 1.63],
'I':[ -0.73,-0.16, 1.79,-0.77,-0.54,-1.78,-0.84,-0.51,-0.55,-0.30],
'L':[ -1.04, 0.00, 1.68,-0.77,-0.60,-1.71,-0.98,-0.52,-0.55,-0.16],
'K':[ -0.34, 0.82, 1.23, 1.70, 1.54,-1.62, 1.15,-0.08,-0.48, 0.60],
'M':[ -1.40, 0.18, 1.30,-0.73,-0.26, 0.53,-1.03,-0.54, 0.26, 0.24],
'F':[ -0.21, 0.98, 0.41,-1.52,-0.45, 1.18, 0.31,-1.02,-0.67,-0.44],
'P':[ 2.06,-0.33,-1.15,-0.75, 0.88,-0.45, 0.30,-2.30, 0.74,-0.28],
'S':[ 0.81,-1.08, 0.16, 0.42,-0.21,-0.43,-1.89, 0.52,-0.11,-0.20],
'T':[ 0.26,-0.70, 1.21, 0.63,-0.10, 0.21, 0.24,-1.15,-0.56, 0.19],
'W':[ 0.30, 2.10,-0.72,-1.57,-1.16, 0.57,-0.48,-0.40,-2.30,-0.60],
'Y':[ 1.38, 1.48, 0.80,-0.56,-0.00, 0.68,-0.31, 1.03,-0.05, 0.53],
'V':[ -0.74,-0.71, 2.04,-0.40, 0.50,-1.18,-0.88,-0.31,-0.02,-0.19]
}

KD = {
'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2
}

def avg_propensity(seq, table):
    return sum(table[aa] for aa in seq) / len(seq)



def avg_embedding(seq, table):
    arr = np.array([table[a] for a in seq])
    return arr.mean(axis=0)

def hydrophobic_autocorr(seq, max_lag=5):
    values = np.array([KD[a] for a in seq])
    feats = {}
    L = len(values)

    for lag in range(1, min(max_lag, L-1) + 1):
        corr = np.mean(values[:-lag] * values[lag:])
        feats[f"hydro_autocorr_{lag}"] = corr

    return feats

def cdr_descriptors(seq: str) -> pd.DataFrame:
    seq = seq.upper()
    pep = peptides.Peptide(seq)
    desc = PyPro.GetProDes(seq)

    f = {}

    # ---------------------------
    # BASIC PROPERTIES
    # ---------------------------
    f["length"] = len(seq)
    f["mw"] = pep.molecular_weight()
    f["pI"] = pep.isoelectric_point()
    f["charge_pH7"] = pep.charge(pH=7.0)
    f["hydrophobicity"] = pep.hydrophobicity()
    f["aliphatic_index"] = pep.aliphatic_index()
    f["instability"] = pep.instability_index()
    f["boman"] = pep.boman()

    # ---------------------------
    # AMINO ACID COMPOSITION
    # ---------------------------
    # aa_freq = pep.frequencies()
    # for aa, val in aa_freq.items():
    #     f[f"freq_{aa}"] = val

    # ---------------------------
    # AROMATIC / BINDING HOTSPOTS
    # ---------------------------
    f["aromatic_fraction"] = sum(seq.count(a) for a in "FWYH") / len(seq)
    f["aromatic_FWY_fraction"] = sum(seq.count(a) for a in "FWY") / len(seq)
    f["tyrosine_fraction"] = seq.count("Y") / len(seq)
    f["tryptophan_fraction"] = seq.count("W") / len(seq)
    f["phenylalanine_fraction"] = seq.count("F") / len(seq)

    # ---------------------------
    # FLEXIBILITY / LOOP FEATURES
    # ---------------------------
    f["gly_fraction"] = seq.count("G") / len(seq)
    f["pro_fraction"] = seq.count("P") / len(seq)

    # ---------------------------
    # CHARGE PATCH FEATURES
    # ---------------------------
    f["positive_fraction"] = sum(seq.count(a) for a in "KRH") / len(seq)
    f["negative_fraction"] = sum(seq.count(a) for a in "DE") / len(seq)

    # ---------------------------
    # HYDROPHOBIC PATCH
    # ---------------------------
    f["hydrophobic_fraction"] = sum(seq.count(a) for a in "AILMFWVY") / len(seq)

    # ---------------------------
    # SECONDARY STRUCTURE PROPENSITY
    # ---------------------------
    f["helix_prop"] = avg_propensity(seq, HELIX)
    f["sheet_prop"] = avg_propensity(seq, SHEET)
    f["turn_prop"]  = avg_propensity(seq, TURN)


    # ---------------------------
    # AMINO ACID PROPERTY EMBEDDINGS
    # ---------------------------
# --- Kidera embedding ---
    kidera_vec = avg_embedding(seq, KIDERA)
    for i,val in enumerate(kidera_vec):
        f[f"kidera_{i}"] = val

    # for i, val in enumerate(pep.vhse()):
    #     f[f"vhse_{i}"] = val

    for i, val in enumerate(pep.z_scales()):
        f[f"zscale_{i}"] = val

    # ---------------------------
    # SEQUENCE ORDER FEATURES (PROPY3)
    # propy3 was written for full proteins, not short peptides.
    # For peptides you should use:
    # lag = min(5, len(seq)-1)   
    # This is standard in peptide QSAR papers.
    # ---------------------------

    # --- Hydrophobic autocorrelation (peptide-safe) ---
    f.update(hydrophobic_autocorr(seq))

    # --- Dipeptide composition (safe always) ---
    #f.update(desc.GetDPComp())

    # --- Safe lag for short peptides ---
    lag = max(1, min(5, len(seq) - 1))

    # --- Pseudo AAC (safe) ---
    paac = desc.GetPAAC(lamda=lag)
    f.update(paac)

    # # --- Autocorrelation descriptors (safe) ---
    # auto = desc.GetMoreauBrotoAuto(nlag=lag)
    # f.update(auto)


    return pd.DataFrame([f])

def compute_cdr_features(cdr1_seq, cdr2_seq, cdr3_seq):
    """
    Compute features for all three CDRs and return a single-row DataFrame.
    
    Args:
        cdr1_seq: CDR1 sequence (string)
        cdr2_seq: CDR2 sequence (string)
        cdr3_seq: CDR3 sequence (string)
    
    Returns:
        pd.DataFrame: Single row with all CDR features, prefixed with CDR1_, CDR2_, CDR3_
    """
    results = {}
    
    for cdr_name, seq in [("CDR1", cdr1_seq), ("CDR2", cdr2_seq), ("CDR3", cdr3_seq)]:
        if seq:
            df_feat = cdr_descriptors(seq)
            row = df_feat.iloc[0].to_dict()
            for k, v in row.items():
                results[f"{cdr_name}_{k}"] = v
    
    return pd.DataFrame([results])

# Example usage (for testing):
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python -m src.09-CDR_peptide_features <pdb_file.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    pdb_stem = Path(pdb_file).stem

    # if not Path(pdb_file).is_absolute():
    #     pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_file

    cdr_seq = PATHS.meta_dir / f"{pdb_stem}" / "cdrs" / f"{pdb_stem}_cdrs.csv"
    print(f"Reading CDR sequences from: {cdr_seq}")

    cdr_df = pd.read_csv(cdr_seq)

    cdr1 = ""
    cdr2 = ""
    cdr3 = ""

    if 'cdr_region' in cdr_df.columns and 'cdr_sequence' in cdr_df.columns:
        cdr_map = (
            cdr_df[['cdr_region', 'cdr_sequence']]
            .dropna()
            .drop_duplicates(subset=['cdr_region'])
            .set_index('cdr_region')['cdr_sequence']
            .to_dict()
        )
        cdr1 = cdr_map.get('CDR1', '')
        cdr2 = cdr_map.get('CDR2', '')
        cdr3 = cdr_map.get('CDR3', '')
    elif 'cdr_region' in cdr_df.columns and 'residue' in cdr_df.columns:
        seq_map = cdr_df.groupby('cdr_region')['residue'].apply(lambda s: ''.join(s.astype(str))).to_dict()
        cdr1 = seq_map.get('CDR1', '')
        cdr2 = seq_map.get('CDR2', '')
        cdr3 = seq_map.get('CDR3', '')
    else:
        raise ValueError("CDR CSV must contain either cdr_region+cdr_sequence or cdr_region+residue columns")

    dataset = compute_cdr_features(cdr1, cdr2, cdr3)
    # save to data/features/cdr_pep_features
    out_path = PATHS.features_dir / f"{pdb_stem}" / "cdrs" / f"{pdb_stem}_cdr_peptide_features.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)  
    dataset.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")
