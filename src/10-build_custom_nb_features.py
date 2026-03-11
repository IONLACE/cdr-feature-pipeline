from __future__ import annotations

import argparse
import glob
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay

from paths import PATHS


@dataclass
class AtomRecord:
    chain: str
    resnum: int
    icode: str
    atom: str
    element: str
    coord: np.ndarray


@dataclass
class PocketRecord:
    pocket_id: str
    features: Dict[str, float]
    residues: set[Tuple[str, int]]
    points: np.ndarray


def _parse_complex_tokens(pdb_stem: str) -> Tuple[str, str]:
    toks = pdb_stem.split("_")
    if len(toks) < 3:
        raise ValueError(f"Cannot infer chains from stem: {pdb_stem}")
    return toks[1], toks[2]


def _safe_float(value: object) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _pick_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _parse_resnum_icode(value: object) -> Tuple[Optional[int], str]:
    if pd.isna(value):
        return None, " "
    txt = str(value).strip()
    m = re.match(r"^(-?\d+)([A-Za-z]?)$", txt)
    if not m:
        return None, " "
    return int(m.group(1)), (m.group(2) if m.group(2) else " ")


def _parse_pdb_atoms(pdb_path: Path) -> List[AtomRecord]:
    atoms: List[AtomRecord] = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            chain = line[21].strip() or " "
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            icode = (line[26].strip() or " ")
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = (line[76:78].strip() if len(line) >= 78 else "") or re.sub(r"\d", "", atom_name)[:1]
            atoms.append(
                AtomRecord(
                    chain=chain,
                    resnum=resnum,
                    icode=icode,
                    atom=atom_name,
                    element=element.upper(),
                    coord=np.array([x, y, z], dtype=float),
                )
            )
    return atoms

def _load_cdr3_positions(
    cdr_csv: Path,
    regions={"CDR3"}
) -> set[Tuple[int, str]]:
    if not cdr_csv.exists():
        raise FileNotFoundError(f"Missing CDR file: {cdr_csv}")
    df = pd.read_csv(cdr_csv)
    # cdr_col = _pick_col(df, ["cdr_region", "cdr"])
    cdr_col = "cdr_region"
    if cdr_col is None:
        raise ValueError(f"Missing cdr_region in {cdr_csv}")
    if regions is not None:
        wanted = {r.upper() for r in regions}
        sub = df[df[cdr_col].astype(str).str.upper().isin(wanted)].copy()
    else:
        sub = df.copy()
    if sub.empty:
        return set()

    pos_col = _pick_col(df, ["pdb_resid", "pdb_position"])
    if pos_col is None:
        raise ValueError(f"CDR file has no parseable residue position column: {cdr_csv}")

    out: set[Tuple[int, str]] = set()
    for v in sub[pos_col].tolist():
        rn, ic = _parse_resnum_icode(v)
        if rn is not None:
            out.add((rn, ic))
    return out


def _load_cdr_region_by_resnum(cdr_csv: Path) -> Dict[int, str]:
    if not cdr_csv.exists():
        raise FileNotFoundError(f"Missing CDR file: {cdr_csv}")

    df = pd.read_csv(cdr_csv)
    cdr_col = "cdr_region"
    pos_col = _pick_col(df, ["pdb_resid", "pdb_position"])
    if cdr_col is None or pos_col is None:
        raise ValueError(f"CDR file has no parseable region/position columns: {cdr_csv}")

    out: Dict[int, str] = {}
    for _, row in df.iterrows():
        region = str(row[cdr_col]).strip().upper()
        if region not in {"CDR1", "CDR2", "CDR3"}:
            continue
        rn, _ = _parse_resnum_icode(row[pos_col])
        if rn is None:
            continue
        if rn not in out:
            out[rn] = region
    return out

def _parse_fpocket_features(pockets_dir):
    """
    Parse an fpocket pocket.txt file into a nested dictionary.

    Returns:
        {
            "Pocket 1": {"Score": 0.216, ...},
            "Pocket 2": {...},
            ...
        }
    """
    pocket_feat = {}
    current_pocket = None

    pocket_path = Path(pockets_dir).parent / Path(pockets_dir).parent.name.replace("_out","_info.txt")

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
                pocket_feat[current_pocket] = {}
                continue

            # Parse key-value lines
            if ":" in line and current_pocket is not None:
                key, value = line.split(":", 1)
                key = key.strip()
                value = float(value.strip())

                if key == "Hydrophobicity score" :
                    pocket_feat[current_pocket]['hydrophob_score'] = value
                if key == "Mean local hydrophobic density" :
                    pocket_feat[current_pocket]['hydrophob_density'] = value                
                if key == "Volume score" :
                    pocket_feat[current_pocket]['vol_score'] = value                
                if key == "Cent. of mass - Alpha Sphere max dist" :
                    pocket_feat[current_pocket]['com_asphere_max_dist'] = value
                if key == "Alpha sphere density" :
                    pocket_feat[current_pocket]['asphere_density'] = value

                if key == "Apolar alpha sphere proportion" :
                    pocket_feat[current_pocket]['Apolar_asphere_prop'] = value
                if key == "Polarity score" :
                    pocket_feat[current_pocket]['pol_score'] = value
                if key == "Charge score" :
                    pocket_feat[current_pocket]['charge_score'] = value
                if key == "Proportion of polar atoms" :
                    pocket_feat[current_pocket]['prop_of_polar_atoms'] = value
                if key == "Polar SASA" :
                    pocket_feat[current_pocket]['polar_SASA'] = value

    return pocket_feat



def _load_pockets(pockets_dir: Path) -> Dict[str, PocketRecord]:
    
    
    feats = _parse_fpocket_features(pockets_dir)
    
    pattern = str(pockets_dir / "pocket*_atm.pdb")
    pocket_files = sorted(glob.glob(pattern))
    out: Dict[str, PocketRecord] = {}

    
    for fp in pocket_files:
        p = Path(fp)
        pocket_id = p.name.split("_")[0]
        residues: set[Tuple[str, int]] = set()
        points: List[np.ndarray] = []
        header: List[str] = []

        with open(p, "r") as handle:
            for line in handle:
                if line.startswith("REMARK"):
                    header.append(line)
                    continue
                if line.startswith("ATOM"):
                    chain = line[21].strip() or " "
                    try:
                        resnum = int(line[22:26].strip())
                    except ValueError:
                        continue
                    residues.add((chain, resnum))
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        continue
                    points.append(np.array([x, y, z], dtype=float))

        
        out[pocket_id] = PocketRecord(
            pocket_id=pocket_id,
            features=feats[pocket_id],
            residues=residues,
            points=np.vstack(points) if points else np.zeros((0, 3), dtype=float),
        )

    return out


def _choose_target_pocket(cdr_pocket_csv: Path, pockets: Dict[str, PocketRecord]) -> Optional[str]:
    if cdr_pocket_csv.exists():
        df = pd.read_csv(cdr_pocket_csv)
        if not df.empty and "pocket_id" in df.columns:
            # cdr3_col = _pick_col(df, ["CDR3_nres_in_pocket", "cdr3_nres_in_pocket"])
            cdr3_col = "CDR3_nres_in_pocket"
            if cdr3_col is not None:
                tmp = df[["pocket_id", cdr3_col]].copy()
                tmp[cdr3_col] = pd.to_numeric(tmp[cdr3_col], errors="coerce").fillna(0.0)
                tmp = tmp.sort_values(cdr3_col, ascending=False)
                if len(tmp) > 0 and float(tmp.iloc[0][cdr3_col]) > 0:
                    return str(tmp.iloc[0]["pocket_id"])

            cdr_cols = [c for c in df.columns if c.lower().endswith("_nres_in_pocket")]
            if cdr_cols:
                tmp = df[["pocket_id", *cdr_cols]].copy()
                for c in cdr_cols:
                    tmp[c] = pd.to_numeric(tmp[c], errors="coerce").fillna(0.0)
                tmp["total"] = tmp[cdr_cols].sum(axis=1)
                tmp = tmp.sort_values("total", ascending=False)
                if len(tmp) > 0 and float(tmp.iloc[0]["total"]) > 0:
                    return str(tmp.iloc[0]["pocket_id"])

    if not pockets:
        return None
    # return the pocket_id with the highest pocket_score feature
    return max(pockets, key=lambda pid: pockets[pid].features.get("pocket_score", float("-inf")))


def _pca_axis(points: np.ndarray) -> np.ndarray:
    if points.shape[0] < 3:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    centered = points - points.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    n = np.linalg.norm(axis)
    if n == 0:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    return axis / n


def _angle_deg(v1: np.ndarray, v2: np.ndarray) -> float:
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return float("nan")
    c = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def _line_point_distance(point: np.ndarray, line_point: np.ndarray, line_dir: np.ndarray) -> float:
    n = np.linalg.norm(line_dir)
    if n == 0:
        return float("nan")
    u = line_dir / n
    diff = point - line_point
    proj = np.dot(diff, u) * u
    perp = diff - proj
    return float(np.linalg.norm(perp))


def _is_heavy(atom: AtomRecord) -> bool:
    return atom.element not in {"H", "D"}


def _compute_cdr3_geometry_features(
    cdr3_atoms: List[AtomRecord],
    pocket: PocketRecord,
) -> Dict[str, float]:
    feats: Dict[str, float] = {}
    if not cdr3_atoms or pocket.points.shape[0] < 4:
        return {
            "cdr3_com_to_pocket_centroid_A": float("nan"),
            "cdr3_com_to_pocket_deepest_A": float("nan"),
            "cdr3_insertion_depth_mean_A": float("nan"),
            "cdr3_insertion_depth_max_A": float("nan"),
            "cdr3_insertion_depth_p90_A": float("nan"),
            "cdr3_penetration_fraction": float("nan"),
            "cdr3_axis_to_pocket_axis_deg": float("nan"),
            "cdr3_tip_to_pocket_axis_deg": float("nan"),
            "cdr3_lateral_offset_A": float("nan"),
        }

    cdr3_heavy = np.vstack([a.coord for a in cdr3_atoms if _is_heavy(a)])
    if cdr3_heavy.shape[0] == 0:
        cdr3_heavy = np.vstack([a.coord for a in cdr3_atoms])

    cdr3_com = cdr3_heavy.mean(axis=0)
    pocket_pts = pocket.points
    pocket_centroid = pocket_pts.mean(axis=0)
    pocket_axis = _pca_axis(pocket_pts)

    # Orient axis towards CDR3 COM
    if np.dot(cdr3_com - pocket_centroid, pocket_axis) < 0:
        pocket_axis = -pocket_axis

    proj_pocket = pocket_pts @ pocket_axis
    mouth_proj = float(np.max(proj_pocket))
    deep_proj = float(np.min(proj_pocket))

    proj_cdr3 = cdr3_heavy @ pocket_axis
    insertion_depths = np.maximum(0.0, mouth_proj - proj_cdr3)

    cdr3_axis = _pca_axis(cdr3_heavy)
    tip_idx = int(np.argmin(proj_cdr3))
    tip = cdr3_heavy[tip_idx]
    tip_vec = tip - cdr3_com
    depth_vec = -pocket_axis

    deep_point = pocket_axis * deep_proj
    mouth_point = pocket_axis * mouth_proj

    feats["cdr3_com_to_pocket_centroid_A"] = float(np.linalg.norm(cdr3_com - pocket_centroid))
    feats["cdr3_com_to_pocket_deepest_A"] = float(np.linalg.norm(cdr3_com - deep_point))
    feats["cdr3_insertion_depth_mean_A"] = float(np.mean(insertion_depths))
    feats["cdr3_insertion_depth_max_A"] = float(np.max(insertion_depths))
    feats["cdr3_insertion_depth_p90_A"] = float(np.percentile(insertion_depths, 90))
    feats["cdr3_axis_to_pocket_axis_deg"] = _angle_deg(cdr3_axis, pocket_axis)
    feats["cdr3_tip_to_pocket_axis_deg"] = _angle_deg(tip_vec, depth_vec)
    feats["cdr3_lateral_offset_A"] = _line_point_distance(cdr3_com, mouth_point, depth_vec)

    try:
        hull = Delaunay(pocket_pts)
        inside = hull.find_simplex(cdr3_heavy) >= 0
        feats["cdr3_penetration_fraction"] = float(np.mean(inside))
    except Exception:
        feats["cdr3_penetration_fraction"] = float("nan")

    return feats


def _load_conservation_map(cons_csv: Path, antigen_chain: str) -> Dict[Tuple[str, int], float]:
    if not cons_csv.exists():
        return {}
    antigen_chains = set(antigen_chain)
    if not antigen_chains:
        return {}

    df = pd.read_csv(cons_csv)
    chain_col = "chain"
    rid_col = "resid"
    score_col = "score"
    if chain_col is None or rid_col is None or score_col is None:
        return {}

    out: Dict[Tuple[str, int], float] = {}
    for _, row in df.iterrows():
        chain = str(row[chain_col]).strip()
        if chain not in antigen_chains:
            continue
        rn = _safe_float(row[rid_col])
        sc = _safe_float(row[score_col])
        if math.isnan(rn) or math.isnan(sc):
            continue
        out[(chain, int(rn))] = float(sc)
    return out


def _conservation_features(contacts_csv: Path, cons_map: Dict[Tuple[str, int], float]) -> Dict[str, float]:
    out = {
        "antigen_contact_cons_mean": float("nan"),
        "antigen_contact_cons_max": float("nan"),
        "antigen_contact_cons_top25_mean": float("nan"),
        "cdr_contact_cons_mean": float("nan"),
        "cdr_contact_cons_max": float("nan"),
        "cdr_contact_cons_top25_mean": float("nan"),
        "cdr3_contact_cons_mean": float("nan"),
        "cdr3_contact_cons_max": float("nan"),
        "cdr3_contact_cons_top25_mean": float("nan"),
        "cdr12_contact_cons_mean": float("nan"),
        "cdr12_contact_cons_max": float("nan"),
        "cdr12_contact_cons_top25_mean": float("nan"),
    }
    if not contacts_csv.exists():
        return out

    df = pd.read_csv(contacts_csv)
    # chain_col = _pick_col(df, ["antigen_chain", "ag_chain", "chain_antigen", "Chain"])
    # res_col = _pick_col(df, ["antigen_resnum", "ag_resnum", "ResNum", "resnum"])
    # cdr_col = _pick_col(df, ["cdr_region", "cdr"])
    chain_col = "antigen_chain"
    res_col = "antigen_resnum"
    cdr_col = "cdr_region"
    if chain_col is None or res_col is None:
        return out

    vals_all: List[float] = []
    vals_cdr3: List[float] = []
    vals_cdr12: List[float] = []

    for _, row in df.iterrows():
        chain = str(row[chain_col]).strip()
        rn = _safe_float(row[res_col])
        if math.isnan(rn):
            continue
        key = (chain, int(rn))
        if key not in cons_map:
            continue
        sc = cons_map[key]
        vals_all.append(sc)
        cdr_label = str(row[cdr_col]).strip().upper() if cdr_col in df.columns else ""
        if cdr_label == "CDR3":
            vals_cdr3.append(sc)
        elif cdr_label in {"CDR1", "CDR2"}:
            vals_cdr12.append(sc)

    def _fill(prefix: str, vals: List[float]):
        if not vals:
            return
        arr = np.asarray(vals, dtype=float)
        out[f"{prefix}_mean"] = float(np.mean(arr))
        out[f"{prefix}_max"] = float(np.max(arr))
        top = np.sort(arr)[-max(1, int(math.ceil(0.25 * len(arr)))) :]
        out[f"{prefix}_top25_mean"] = float(np.mean(top))

    _fill("antigen_contact_cons", vals_all)
    vals_cdr = vals_cdr12 + vals_cdr3
    _fill("cdr_contact_cons", vals_cdr)
    _fill("cdr3_contact_cons", vals_cdr3)
    _fill("cdr12_contact_cons", vals_cdr12)
    return out


def _paratope_specificity_features(contacts_csv: Path) -> Dict[str, float]:
    out = {
        "paratope_specificity_ratio": float("nan"),
        "framework_contact_fraction": float("nan"),
        "cdr1_contact_fraction": float("nan"),
        "cdr2_contact_fraction": float("nan"),
        "cdr3_contact_fraction": float("nan"),
    }
    if not contacts_csv.exists():
        return out

    df = pd.read_csv(contacts_csv)
    cdr_col = "cdr_region"
    if cdr_col is None or df.empty:
        return out

    labels = df[cdr_col].astype(str).str.upper().str.strip()
    total = len(labels)
    if total == 0:
        return out

    cdr_mask = labels.isin(["CDR1", "CDR2", "CDR3"])
    cdr_count = int(cdr_mask.sum())
    out["paratope_specificity_ratio"] = cdr_count / total
    out["framework_contact_fraction"] = 1.0 - out["paratope_specificity_ratio"]

    for cdr in ["CDR1", "CDR2", "CDR3"]:
        out[f"{cdr.lower()}_contact_fraction"] = float((labels == cdr).sum() / total)
    return out


def _estimate_dsasa(
    dsasa_interface_csv: Path,
    nb_chain: str,
    cdr3_positions: set[Tuple[int, str]],
    cdr_region_by_nb_resnum: Optional[Dict[int, str]] = None,
) -> Dict[str, float]:
    out = {
        "dsasa_total_A2": float("nan"),
        "dsasa_nanobody_chain_A2": float("nan"),
        "dsasa_cdr1_A2": float("nan"),
        "dsasa_cdr2_A2": float("nan"),
        "dsasa_cdr3_A2": float("nan"),
        "dsasa_framework_A2": float("nan"),
    }
    if not dsasa_interface_csv.exists():
        return out

    df = pd.read_csv(dsasa_interface_csv)
    chain_col = "chain"
    res_col = "resid"
    delta_col = "dSASA"

    required_cols = [chain_col, res_col, delta_col]
    if any(col not in df.columns for col in required_cols):
        return out

    vals = pd.to_numeric(df[delta_col], errors="coerce")

    vals = vals.abs().fillna(0.0)
    out["dsasa_total_A2"] = float(vals.sum())

    chains = df[chain_col].astype(str).str.strip()
    nb_mask = chains == nb_chain
    out["dsasa_nanobody_chain_A2"] = float(vals[nb_mask].sum())

    parsed_rn = df[res_col].apply(_parse_resnum_from_residue_text)
    cdr_region_map = {int(k): str(v).upper() for k, v in (cdr_region_by_nb_resnum or {}).items()}
    cdr3_nums = {x[0] for x in cdr3_positions}

    cdr1_mask = nb_mask & parsed_rn.map(lambda x: x in cdr_region_map and cdr_region_map[x] == "CDR1")
    cdr2_mask = nb_mask & parsed_rn.map(lambda x: x in cdr_region_map and cdr_region_map[x] == "CDR2")
    cdr3_mask = nb_mask & parsed_rn.map(
        lambda x: (x in cdr_region_map and cdr_region_map[x] == "CDR3") or (x in cdr3_nums)
    )

    cdr_any_mask = cdr1_mask | cdr2_mask | cdr3_mask
    framework_mask = nb_mask & (~cdr_any_mask)

    out["dsasa_cdr1_A2"] = float(vals[cdr1_mask].sum())
    out["dsasa_cdr2_A2"] = float(vals[cdr2_mask].sum())
    out["dsasa_cdr3_A2"] = float(vals[cdr3_mask].sum())
    out["dsasa_framework_A2"] = float(vals[framework_mask].sum())

    return out


def _parse_resnum_from_residue_text(value: object) -> Optional[int]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    txt = str(value).strip()
    m = re.search(r"(-?\d+)([A-Za-z]?)$", txt)
    if m:
        try:
            return int(m.group(1))
        except Exception:
            return None
    rn, _ = _parse_resnum_icode(value)
    return rn


def _interaction_bucket(interaction_type: str) -> Optional[str]:
    s = interaction_type.lower().strip()
    if not s:
        return None
    if "hydrogen bond" in s or re.search(r"h[ -]?bond", s):
        return "hbond"
    if "salt bridge" in s or "ionic" in s:
        return "salt_bridge"
    if "cation" in s and "pi" in s:
        return "cation_pi"
    if "pi-pi" in s or "pi pi" in s or "stack" in s:
        return "pi_pi"
    if "hydrophob" in s:
        return "hydrophobic"
    return None


def _extract_interaction_counts(
    pandaprot_csv: Path,
    nb_chain: str,
    antigen_chain: str,
    cdr_region_by_nb_resnum: Optional[Dict[int, str]] = None,
    pocket_antigen_resnums: Optional[set[int]] = None,
) -> Dict[str, float]:
    interaction_types = ["hbond", "salt_bridge", "pi_pi", "cation_pi", "hydrophobic"]
    regions = ["cdr1", "cdr2", "cdr3", "framework"]

    out: Dict[str, float] = {}
    for itype in interaction_types:
        out[f"{itype}_count"] = 0.0
        out[f"{itype}_target_pocket_count"] = 0.0
        for region in regions:
            out[f"{itype}_{region}_count"] = 0.0
            out[f"{itype}_{region}_target_pocket_count"] = 0.0

    if not pandaprot_csv.exists():
        return out

    df = pd.read_csv(pandaprot_csv)
    if df.empty:
        return out

    type_col = _pick_col(df, ["Interaction_Type", "interaction_type", "type"])
    chain1_col = _pick_col(df, ["Chain1", "chain1"])
    chain2_col = _pick_col(df, ["Chain2", "chain2"])
    res1_col = _pick_col(df, ["Residue1", "residue1", "res1"])
    res2_col = _pick_col(df, ["Residue2", "residue2", "res2"])
    atom1_col = _pick_col(df, ["Atom1", "atom1"])
    atom2_col = _pick_col(df, ["Atom2", "atom2"])
    dist_col = _pick_col(df, ["Distance_Å", "Distance_A", "distance", "distance_A"])

    cdr_region_map = {k: v.upper() for k, v in (cdr_region_by_nb_resnum or {}).items()}
    pocket_nums = set(pocket_antigen_resnums or set())

    if type_col is None:
        text_rows = df.fillna("").astype(str).agg(" ".join, axis=1)
        for txt in text_rows.tolist():
            bucket = _interaction_bucket(str(txt))
            if bucket is not None:
                out[f"{bucket}_count"] += 1.0
        return out

    antigen_chains = set(antigen_chain)
    if not antigen_chains:
        return out

    seen_keys: set[Tuple[object, ...]] = set()
    for _, row in df.iterrows():
        bucket = _interaction_bucket(str(row[type_col]))
        if bucket is None:
            continue

        if chain1_col is None or chain2_col is None or res1_col is None or res2_col is None:
            out[f"{bucket}_count"] += 1.0
            continue

        c1 = str(row[chain1_col]).strip()
        c2 = str(row[chain2_col]).strip()
        r1 = _parse_resnum_from_residue_text(row[res1_col])
        r2 = _parse_resnum_from_residue_text(row[res2_col])
        if r1 is None or r2 is None:
            continue

        if c1 == nb_chain and c2 in antigen_chains:
            nb_rn, ag_rn = r1, r2
            ag_chain_row = c2
            nb_atom = str(row[atom1_col]).strip() if atom1_col is not None else ""
            ag_atom = str(row[atom2_col]).strip() if atom2_col is not None else ""
        elif c2 == nb_chain and c1 in antigen_chains:
            nb_rn, ag_rn = r2, r1
            ag_chain_row = c1
            nb_atom = str(row[atom2_col]).strip() if atom2_col is not None else ""
            ag_atom = str(row[atom1_col]).strip() if atom1_col is not None else ""
        else:
            continue

        region_label = cdr_region_map.get(nb_rn, "FRAMEWORK")
        region_key = region_label.lower()
        if region_key not in {"cdr1", "cdr2", "cdr3", "framework"}:
            region_key = "framework"

        dist_val = _safe_float(row[dist_col]) if dist_col is not None else float("nan")
        dist_key = round(float(dist_val), 3) if not math.isnan(dist_val) else "nan"

        endpoint_a = (ag_chain_row, ag_rn, ag_atom)
        endpoint_b = (nb_chain, nb_rn, nb_atom)
        if endpoint_a <= endpoint_b:
            left, right = endpoint_a, endpoint_b
        else:
            left, right = endpoint_b, endpoint_a

        key = (bucket, left, right, dist_key)
        if key in seen_keys:
            continue
        seen_keys.add(key)

        out[f"{bucket}_{region_key}_count"] += 1.0
        if region_key != "framework":
            out[f"{bucket}_count"] += 1.0

        in_target_pocket = (not pocket_nums) or (ag_rn in pocket_nums)
        if in_target_pocket:
            out[f"{bucket}_{region_key}_target_pocket_count"] += 1.0
            if region_key != "framework":
                out[f"{bucket}_target_pocket_count"] += 1.0

    return out


def _parse_interface_complementarity_metrics(interface_csv: Path) -> Dict[str, float]:
    out: Dict[str, float] = {
        "surface_complementarity_sc": float("nan"),
        "contact_density": float("nan"),
        "mean_distance": float("nan"),
        "gap_volume_estimate": float("nan"),
        "packing_uniformity": float("nan"),
        "interface_compactness": float("nan"),
        "chemical_complementarity": float("nan"),
    }
    if not interface_csv.exists():
        return out
    df = pd.read_csv(interface_csv)
    if df.empty:
        return out

    row = df.iloc[0]
    out["surface_complementarity_sc"] = _safe_float(row.get("shape_complementarity_sc_score", float("nan")))
    out["contact_density"] = _safe_float(row.get("packing_statistics_contact_density", float("nan")))
    out["mean_distance"] = _safe_float(row.get("packing_statistics_mean_distance", float("nan")))
    out["gap_volume_estimate"] = _safe_float(row.get("packing_statistics_gap_volume_estimate", float("nan")))
    out["packing_uniformity"] = _safe_float(row.get("packing_statistics_packing_uniformity", float("nan")))
    out["interface_compactness"] = _safe_float(row.get("compactness_compactness_ratio", float("nan")))
    out["chemical_complementarity"] = _safe_float(row.get("chemical_complementarity_chemical_score", float("nan")))
    return out


def _load_selected_cdr_peptide_features(cdr_pep_csv: Path) -> Dict[str, float]:
    selected_metrics = [
        "length",
        "charge_pH7",
        "positive_fraction",
        "negative_fraction",
        "hydrophobic_fraction",
        "aromatic_FWY_fraction",
        "tryptophan_fraction",
        "tyrosine_fraction",
        "gly_fraction",
        "pro_fraction",
        "turn_prop",
        "boman",
        "instability",
    ]
    selected = [f"{cdr}_{metric}" for cdr in ["CDR1", "CDR2", "CDR3"] for metric in selected_metrics]

    out: Dict[str, float] = {}
    if not cdr_pep_csv.exists():
        return out

    df = pd.read_csv(cdr_pep_csv)
    if df.empty:
        return out

    row = df.iloc[0]
    for col in selected:
        if col in df.columns:
            out[col] = _safe_float(row[col])

    cdr1_len = out.get("CDR1_length", float("nan"))
    cdr2_len = out.get("CDR2_length", float("nan"))
    cdr3_len = out.get("CDR3_length", float("nan"))
    if not math.isnan(cdr1_len) and not math.isnan(cdr2_len) and not math.isnan(cdr3_len):
        denom = cdr1_len + cdr2_len
        out["CDR3_len_over_CDR1plus2"] = float(cdr3_len / denom) if denom > 0 else float("nan")

    pos = out.get("CDR3_positive_fraction", float("nan"))
    neg = out.get("CDR3_negative_fraction", float("nan"))
    if not math.isnan(pos) and not math.isnan(neg):
        out["CDR3_net_charge_fraction_proxy"] = float(pos - neg)

    return out


def build_features_for_complex(pdb_file: str, output_csv: Optional[str] = None) -> Path:
    pdb_stem = Path(pdb_file).stem
    # nb_chain, ag_chain = _parse_complex_tokens(pdb_stem)
    _, nb_chain, ag_chain, _ = pdb_stem.split("_")
    
    meta_dir = PATHS.meta_dir / pdb_stem
    feat_dir = PATHS.features_dir

    cdr_csv = meta_dir / "cdrs" / f"{pdb_stem}_cdrs.csv"
    pockets_dir = meta_dir / "pocket" / f"{pdb_stem}_out" / "pockets"
    cdr_pocket_csv = meta_dir / "contacts" / f"{pdb_stem}_CDR_pocket_contacts.csv"
    contacts_csv = meta_dir / "contacts" / f"{pdb_stem}_all_and_CDR_contacts.csv"
    cons_csv = meta_dir / f"{pdb_stem}_cons_mapped.csv"
    dsasa_interface_csv = meta_dir / "sasa" / "dsasa.csv"
    pandaprot_csv = meta_dir / "pandaprot_run" / "pandaprot_report.csv"
    interface_csv = meta_dir / "interface" / f"{pdb_stem}_interface_results.csv"
    cdr_pep_csv = meta_dir / "cdrs" / f"{pdb_stem}_cdr_peptide_features.csv"

    if not Path(pdb_file).is_absolute():
        pdb_file = str(PATHS.inp_pdb_dir / pdb_stem.split("_")[0] / pdb_file)

    cdr3_positions = _load_cdr3_positions(cdr_csv, regions={"CDR3"})
    cdr_region_by_nb_resnum = _load_cdr_region_by_resnum(cdr_csv)
    atoms = _parse_pdb_atoms(Path(pdb_file))
    cdr3_atoms = [a for a in atoms if a.chain == nb_chain and (a.resnum, a.icode) in cdr3_positions]

    pockets = _load_pockets(pockets_dir)
    
    # get the pocket id with max contacts to CDR3
    # else get the pocket id which has max contacts
    # else get the pocket id with max poxket score
    target_pocket = _choose_target_pocket(cdr_pocket_csv, pockets)

    base: Dict[str, float | str] = {
        "ID": pdb_stem,
        "nb_chain": nb_chain,
        "antigen_chain": ag_chain,
        "target_pocket": target_pocket or "",
    }

    if target_pocket and target_pocket in pockets:
        pocket = pockets[target_pocket]
        base.update(pocket.features)
        geom = _compute_cdr3_geometry_features(cdr3_atoms, pocket)
        base.update(geom)
        if "pocket_volume" in pocket.features and "cdr3_penetration_fraction" in geom:
            base["cdr3_pocket_overlap_volume_proxy_A3"] = float(
                pocket.features.get("pocket_volume", float("nan"))
            ) * float(geom.get("cdr3_penetration_fraction", float("nan")))
    else:
        base.update(_compute_cdr3_geometry_features(cdr3_atoms, PocketRecord("", {}, set(), np.zeros((0, 3)))))
        base["cdr3_pocket_overlap_volume_proxy_A3"] = float("nan")

    base.update(_parse_interface_complementarity_metrics(interface_csv))

    cons_map = _load_conservation_map(cons_csv, ag_chain)
    base.update(_conservation_features(contacts_csv, cons_map))
    base.update(_paratope_specificity_features(contacts_csv))

    dsasa = _estimate_dsasa(dsasa_interface_csv, nb_chain, cdr3_positions, cdr_region_by_nb_resnum)
    base.update(dsasa)

    base.update(_load_selected_cdr_peptide_features(cdr_pep_csv))

    pocket_ag_resnums: set[int] = set()
    if target_pocket and target_pocket in pockets:
        antigen_chains = set(ag_chain)
        pocket_ag_resnums = {rn for ch, rn in pockets[target_pocket].residues if ch in antigen_chains}

    interactions = _extract_interaction_counts(
        pandaprot_csv,
        nb_chain=nb_chain,
        antigen_chain=ag_chain,
        cdr_region_by_nb_resnum=cdr_region_by_nb_resnum,
        pocket_antigen_resnums=pocket_ag_resnums,
    )
    base.update(interactions)

    dsasa_total = _safe_float(base.get("dsasa_total_A2", float("nan")))
    if not math.isnan(dsasa_total) and dsasa_total > 0:
        base["hbond_density_per_A2"] = interactions["hbond_count"] / dsasa_total
        base["salt_bridge_density_per_A2"] = interactions["salt_bridge_count"] / dsasa_total
        base["pi_pi_density_per_A2"] = interactions["pi_pi_count"] / dsasa_total
        base["cation_pi_density_per_A2"] = interactions["cation_pi_count"] / dsasa_total
    else:
        base["hbond_density_per_A2"] = float("nan")
        base["salt_bridge_density_per_A2"] = float("nan")
        base["pi_pi_density_per_A2"] = float("nan")
        base["cation_pi_density_per_A2"] = float("nan")

    out_df = pd.DataFrame([base])

    if output_csv:
        out_path = Path(output_csv)
    else:
        out_path = feat_dir / f"{pdb_stem}_custom_features.csv"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build custom nanobody-antigen features from existing outputs.")
    parser.add_argument("pdb_file", type=str, help="PDB filename or absolute path")
    parser.add_argument("--output_csv", type=str, default=None, help="Optional explicit output CSV path")
    args = parser.parse_args()

    out = build_features_for_complex(args.pdb_file, args.output_csv)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
