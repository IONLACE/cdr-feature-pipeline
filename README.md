# cdr-feature-pipeline

Pipeline for building features for nanobody-antigen scoring using ML.

This repository builds structural, interface, conservation, and interaction features from structural models for supervised machine learning model development, model selection, and confidence-based ranking of predicted nanobody-antigen complexes, then merges outputs into ML-ready training/inference tables.

## Project Layout

- `src/`: pipeline scripts
- `utils/`: helpers (subprocess wrappers, etc.)
- `config.yaml`: path configuration
- `paths.py`: resolved path helper based on `config.yaml`
- `run_all_README.md`: example run logs

## Configuration: `config.yaml` and `paths.py`

All scripts resolve input/output locations through `PATHS` from `paths.py`, which is built from `config.yaml`.

Current `config.yaml` keys:

- `out_dir`: root output workspace
- `inp_pdb`: input-PDB subdirectory name
- `meta_subdir`: metadata/intermediate subdirectory name
- `features_subdir`: final feature subdirectory name
- `logs_subdir`: logs subdirectory name
- `paths.*`: expanded path templates using `${...}` placeholders

Effective locations are:

- `PATHS.inp_pdb_dir` -> where input PDBs are read
- `PATHS.meta_dir` -> where intermediate artifacts are written (`msa`, `cdrs`, `contacts`, `interface`, `sasa`, etc.)
- `PATHS.features_dir` -> where per-complex and merged feature CSVs are written
- `PATHS.logs_dir` -> where logs can be written

### Flexibility of input/output locations

You can relocate all pipeline data without changing Python code:

1. Edit `out_dir` to point to a different root (absolute or relative).
2. Optionally rename subfolders via `inp_pdb`, `meta_subdir`, `features_subdir`, `logs_subdir`.
3. Keep `paths.*` templates consistent with your chosen layout.

Example (conceptual):

```yaml
out_dir: /mnt/project_runs/zff_data
inp_pdb: input_pdbs
meta_subdir: meta
features_subdir: features
logs_subdir: logs
paths:
  inp_pdb_dir: ${out_dir}/${inp_pdb}
  meta_dir: ${out_dir}/${meta_subdir}
  features_dir: ${out_dir}/${features_subdir}
  logs_dir: ${out_dir}/${logs_subdir}
```

`paths.py` expands placeholders and supports both absolute and relative paths. If a configured path is relative, it is resolved against the repository root.

## Input Naming Convention

Most scripts expect an input filename pattern like:

`<pdbid>_<nanobody_chain>_<antigen_chain(s)>_<tag>.pdb`

Example:

`8q6k_C_BA_Repair.pdb`

Where:

- `pdbid = 8q6k`
- `nanobody_chain = C`
- `antigen_chain(s) = BA`

## Pipeline Overview

This pipeline is designed to produce robust, ML-consumable feature tables to train and apply model-quality/confidence predictors for nanobody-antigen complex ranking.

### One-command wrapper

You can run the full pipeline for one complex with a single command:

```bash
python run_pipeline.py 8q6k_C_BA_Repair.pdb
```

Default behavior:

- Starts at step 02 (assumes MSA already exists)
- Runs through step 10
- Runs step 11 merge (`custom_features_merged_all.csv`)

Useful options:

```bash
# Include experimental step 01
python run_pipeline.py 8q6k_C_BA_Repair.pdb --with-msa

# Write step-10 output to an explicit path
python run_pipeline.py 8q6k_C_BA_Repair.pdb --output-csv /path/to/out.csv

# Skip step 11 merge
python run_pipeline.py 8q6k_C_BA_Repair.pdb --skip-merge-features

# Run X/Y merge at the end (default ylabels path)
python run_pipeline.py 8q6k_C_BA_Repair.pdb --merge-xy

# Run X/Y merge with explicit y-labels path
python run_pipeline.py 8q6k_C_BA_Repair.pdb --merge-xy --y-labels /path/to/ylabels.csv
```

### 1) Build MSA (experimental)

`01-build_msa.py` is still experimental and needs more testing.

Command:

```bash
python -m src.01-build_msa 8q6k_C_BA_Repair.pdb
```

### 2) Start from conservation if MSA is already computed

If you already have antigen-chain MSA files (for example generated via NCBI `blastp` + COBALT), you can start directly from step 02.

Expected MSA files under:

`<meta>/<complex_id>/msa/`

Format expected by step 02:

- `ag_A_ncbi.fa`
- `ag_B_ncbi.fa`
- ... one per antigen chain

Then run:

```bash
python -m src.02-calc_conservation 8q6k_C_BA_Repair.pdb
python -m src.03-get_CDRs 8q6k_C_BA_Repair.pdb
python -m src.04-run_fpocket 8q6k_C_BA_Repair.pdb
python -m src.05-get_all_contacts 8q6k_C_BA_Repair.pdb
python -m src.06-run_pandaprot 8q6k_C_BA_Repair.pdb
python -m src.07-interface_complimentarity 8q6k_C_BA_Repair.pdb
python -m src.08-compute_sasa 8q6k_C_BA_Repair.pdb
python -m src.09-CDR_peptide_features 8q6k_C_BA_Repair.pdb
python -m src.10-build_custom_nb_features 8q6k_C_BA_Repair.pdb
python -m src.11-merge_features
```

Optional final label merge:

```bash
python -m src.merge_X_Y
# or provide y-labels explicitly
python -m src.merge_X_Y /path/to/ylabels.csv
```

## What Each Python Script Does (`src/`)

- `01-build_msa.py`
  - Fetches antigen sequences from RCSB and builds chain-wise MSA assets (BLAST/CD-HIT/Clustal Omega flow).
  - Status: experimental; validate outputs before production use.

- `02-calc_conservation.py`
  - Cleans per-chain MSA against query-gap columns.
  - Runs `cons-capra07` to compute residue conservation scores.
  - Writes combined conservation mapping CSV: `<meta>/<complex_id>/<complex_id>_cons_mapped.csv`.
  - Fails fast if expected `ag_<chain>_ncbi.fa` input is missing.

- `03-get_CDRs.py`
  - Detects CDR1/CDR2/CDR3 using IMGT numbering (ANARCI) and maps to PDB residues.
  - Writes CDR annotation CSV under `<meta>/<complex_id>/cdrs/`.

- `04-run_fpocket.py`
  - Runs fpocket and parses pocket-level descriptors and pocket residues.

- `05-get_all_contacts.py`
  - Builds nanobody-antigen contact tables and CDR-pocket contact summaries.

- `06-run_pandaprot.py`
  - Runs PandaProt and stores filtered interaction report for downstream counts.

- `07-interface_complimentarity.py`
  - Computes interface complementarity/packing/compactness/chemical metrics.
  - Writes `<complex_id>_interface_results.csv`.

- `08-compute_sasa.py`
  - Computes dSASA at residue level.
  - Canonical output: `<meta>/<complex_id>/sasa/dsasa.csv` with columns:
    - `partner, chain, resname, resid, sasa_unbound, sasa_complex, dSASA`

- `09-CDR_peptide_features.py`
  - Computes peptide-like CDR descriptors (charge, composition, turn propensity, etc.).
  - Writes `<features>/<complex_id>/cdrs/<complex_id>_cdr_peptide_features.csv`.

- `10-build_custom_nb_features.py`
  - Main feature assembly step.
  - Aggregates pocket, CDR geometry, conservation, interface, dSASA, interaction counts, and peptide CDR descriptors into one row per complex.
  - Writes `<features>/<complex_id>_custom_features.csv` (unless `--output_csv` is provided).

- `11-merge_features.py`
  - Recursively finds and merges all `*_custom_features.csv` under features directory.
  - Writes `<features>/custom_features_merged_all.csv`.

- `merge_X_Y.py`
  - Merges merged X features with y-labels CSV by key (`pdb_id` or `ID`, auto-detected).
  - Default y-label path: `<inp_pdb_dir>/ylabels.csv`.
  - Optional CLI override: `sys.argv[1]`.

## Environment Notes

- Several scripts require external tools/packages (`blastp`, `blastdbcmd`, `cd-hit`, `clustalo`, `fpocket`, PandaProt, `cons-capra07`, ANARCI, etc.).
- Use the same environment consistently when running the full pipeline to avoid dependency drift.

## Testing

Use both local tests and GitHub workflows.

### Local tests

Run unit tests with pytest:

```bash
pytest -q
```

You can also run a subset during development:

```bash
pytest -q tests/test_run_pipeline.py
```

### GitHub Actions checks

- `CI` workflow (`.github/workflows/ci.yml`):
  - Runs `pytest -q` on Python 3.10, 3.11, and 3.12.
  - This is the main automated unit-test gate.

- `Docker P1/P2 Build and Push` workflows:
  - Build container images.
  - Run smoke checks only (dependency import check and `run_pipeline.py --help` via container entrypoint).
  - These workflows do not run pytest inside the image.

Recommended pre-PR local check:

```bash
pytest -q
```

## Outputs

Main outputs in `<features_dir>` for machine learning model development and confidence-based model ranking:

- Per-complex: `<complex_id>_custom_features.csv`
- All complexes merged: `custom_features_merged_all.csv`
- Optional X/Y merged table: `merged_features.csv`

## Label And Feature Schema

Use `DockQ` and/or `fnat` from labels as training targets. Keep `dG` as an input feature (not as target) if training a confidence/ranking model.

### Label CSV (`ylabels.csv`)

| Column | Role | Description |
|---|---|---|
| ID | Key | Complex/model identifier for joining with features. |
| DockQ | Label | Primary docking quality target (continuous). |
| iRMSD | Aux Label | Interface RMSD quality signal. |
| LRMSD | Aux Label | Ligand RMSD quality signal. |
| fnat | Label | Fraction of native contacts; useful ranking target. |
| dG | Feature | EvoEF binding energy; recommended as model input feature. |

### Feature CSV (`custom_features_merged_all.csv`)

Total feature columns: **148**

| Column | Description |
|---|---|
| ID | Complex/model identifier used as join key. |
| nb_chain | Nanobody chain ID in the complex. |
| antigen_chain | Antigen chain ID(s) in the complex. |
| target_pocket | Selected fpocket pocket ID used for pocket-centric features. |
| polar_SASA | fpocket-derived pocket physicochemical/shape descriptor. |
| hydrophob_density | fpocket-derived pocket physicochemical/shape descriptor. |
| Apolar_asphere_prop | fpocket-derived pocket physicochemical/shape descriptor. |
| hydrophob_score | fpocket-derived pocket physicochemical/shape descriptor. |
| vol_score | fpocket-derived pocket physicochemical/shape descriptor. |
| pol_score | fpocket-derived pocket physicochemical/shape descriptor. |
| charge_score | fpocket-derived pocket physicochemical/shape descriptor. |
| prop_of_polar_atoms | fpocket-derived pocket physicochemical/shape descriptor. |
| asphere_density | fpocket-derived pocket physicochemical/shape descriptor. |
| com_asphere_max_dist | fpocket-derived pocket physicochemical/shape descriptor. |
| cdr3_com_to_pocket_centroid_A | CDR3-to-pocket geometry or penetration metric. |
| cdr3_com_to_pocket_deepest_A | CDR3-to-pocket geometry or penetration metric. |
| cdr3_insertion_depth_mean_A | CDR3-to-pocket geometry or penetration metric. |
| cdr3_insertion_depth_max_A | CDR3-to-pocket geometry or penetration metric. |
| cdr3_insertion_depth_p90_A | CDR3-to-pocket geometry or penetration metric. |
| cdr3_axis_to_pocket_axis_deg | CDR3-to-pocket geometry or penetration metric. |
| cdr3_tip_to_pocket_axis_deg | CDR3-to-pocket geometry or penetration metric. |
| cdr3_lateral_offset_A | CDR3-to-pocket geometry or penetration metric. |
| cdr3_penetration_fraction | CDR3-to-pocket geometry or penetration metric. |
| surface_complementarity_sc | Interface complementarity/packing metric from interface analysis. |
| contact_density | Interface complementarity/packing metric from interface analysis. |
| mean_distance | Interface complementarity/packing metric from interface analysis. |
| gap_volume_estimate | Interface complementarity/packing metric from interface analysis. |
| packing_uniformity | Interface complementarity/packing metric from interface analysis. |
| interface_compactness | Interface complementarity/packing metric from interface analysis. |
| chemical_complementarity | Interface complementarity/packing metric from interface analysis. |
| antigen_contact_cons_mean | Conservation score summary over contacting residues for the indicated subset. |
| antigen_contact_cons_max | Conservation score summary over contacting residues for the indicated subset. |
| antigen_contact_cons_top25_mean | Conservation score summary over contacting residues for the indicated subset. |
| cdr_contact_cons_mean | Conservation score summary over contacting residues for the indicated subset. |
| cdr_contact_cons_max | Conservation score summary over contacting residues for the indicated subset. |
| cdr_contact_cons_top25_mean | Conservation score summary over contacting residues for the indicated subset. |
| cdr3_contact_cons_mean | Conservation score summary over contacting residues for the indicated subset. |
| cdr3_contact_cons_max | Conservation score summary over contacting residues for the indicated subset. |
| cdr3_contact_cons_top25_mean | Conservation score summary over contacting residues for the indicated subset. |
| cdr12_contact_cons_mean | Conservation score summary over contacting residues for the indicated subset. |
| cdr12_contact_cons_max | Conservation score summary over contacting residues for the indicated subset. |
| cdr12_contact_cons_top25_mean | Conservation score summary over contacting residues for the indicated subset. |
| paratope_specificity_ratio | Fraction of contacts by nanobody region (CDRs/framework). |
| framework_contact_fraction | Fraction of contacts by nanobody region (CDRs/framework). |
| cdr1_contact_fraction | Fraction of contacts by nanobody region (CDRs/framework). |
| cdr2_contact_fraction | Fraction of contacts by nanobody region (CDRs/framework). |
| cdr3_contact_fraction | Fraction of contacts by nanobody region (CDRs/framework). |
| dsasa_total_A2 | dSASA-derived buried surface area feature (A^2). |
| dsasa_nanobody_chain_A2 | dSASA-derived buried surface area feature (A^2). |
| dsasa_cdr1_A2 | dSASA-derived buried surface area feature (A^2). |
| dsasa_cdr2_A2 | dSASA-derived buried surface area feature (A^2). |
| dsasa_cdr3_A2 | dSASA-derived buried surface area feature (A^2). |
| dsasa_framework_A2 | dSASA-derived buried surface area feature (A^2). |
| CDR1_length | CDR peptide length (residues). |
| CDR1_charge_pH7 | Estimated net charge of CDR peptide at pH 7. |
| CDR1_positive_fraction | Fraction of positively charged residues in CDR peptide. |
| CDR1_negative_fraction | Fraction of negatively charged residues in CDR peptide. |
| CDR1_hydrophobic_fraction | Fraction of hydrophobic residues in CDR peptide. |
| CDR1_aromatic_FWY_fraction | Fraction of aromatic residues (F/W/Y) in CDR peptide. |
| CDR1_tryptophan_fraction | Fraction of Trp residues in CDR peptide. |
| CDR1_tyrosine_fraction | Fraction of Tyr residues in CDR peptide. |
| CDR1_gly_fraction | Fraction of Gly residues in CDR peptide. |
| CDR1_pro_fraction | Fraction of Pro residues in CDR peptide. |
| CDR1_turn_prop | Average turn propensity of CDR peptide. |
| CDR1_boman | Boman index proxy for protein-binding potential. |
| CDR1_instability | Instability index of CDR peptide. |
| CDR2_length | CDR peptide length (residues). |
| CDR2_charge_pH7 | Estimated net charge of CDR peptide at pH 7. |
| CDR2_positive_fraction | Fraction of positively charged residues in CDR peptide. |
| CDR2_negative_fraction | Fraction of negatively charged residues in CDR peptide. |
| CDR2_hydrophobic_fraction | Fraction of hydrophobic residues in CDR peptide. |
| CDR2_aromatic_FWY_fraction | Fraction of aromatic residues (F/W/Y) in CDR peptide. |
| CDR2_tryptophan_fraction | Fraction of Trp residues in CDR peptide. |
| CDR2_tyrosine_fraction | Fraction of Tyr residues in CDR peptide. |
| CDR2_gly_fraction | Fraction of Gly residues in CDR peptide. |
| CDR2_pro_fraction | Fraction of Pro residues in CDR peptide. |
| CDR2_turn_prop | Average turn propensity of CDR peptide. |
| CDR2_boman | Boman index proxy for protein-binding potential. |
| CDR2_instability | Instability index of CDR peptide. |
| CDR3_length | CDR peptide length (residues). |
| CDR3_charge_pH7 | Estimated net charge of CDR peptide at pH 7. |
| CDR3_positive_fraction | Fraction of positively charged residues in CDR peptide. |
| CDR3_negative_fraction | Fraction of negatively charged residues in CDR peptide. |
| CDR3_hydrophobic_fraction | Fraction of hydrophobic residues in CDR peptide. |
| CDR3_aromatic_FWY_fraction | Fraction of aromatic residues (F/W/Y) in CDR peptide. |
| CDR3_tryptophan_fraction | Fraction of Trp residues in CDR peptide. |
| CDR3_tyrosine_fraction | Fraction of Tyr residues in CDR peptide. |
| CDR3_gly_fraction | Fraction of Gly residues in CDR peptide. |
| CDR3_pro_fraction | Fraction of Pro residues in CDR peptide. |
| CDR3_turn_prop | Average turn propensity of CDR peptide. |
| CDR3_boman | Boman index proxy for protein-binding potential. |
| CDR3_instability | Instability index of CDR peptide. |
| CDR3_len_over_CDR1plus2 | Derived CDR3 descriptor ratio/proxy from CDR peptide features. |
| CDR3_net_charge_fraction_proxy | Derived CDR3 descriptor ratio/proxy from CDR peptide features. |
| hbond_count | Total hbond interactions involving CDR regions. |
| hbond_target_pocket_count | hbond interactions where antigen residue is in target pocket. |
| hbond_cdr1_count | hbond interactions for nanobody region cdr1. |
| hbond_cdr1_target_pocket_count | hbond interactions for cdr1 with antigen residue in target pocket. |
| hbond_cdr2_count | hbond interactions for nanobody region cdr2. |
| hbond_cdr2_target_pocket_count | hbond interactions for cdr2 with antigen residue in target pocket. |
| hbond_cdr3_count | hbond interactions for nanobody region cdr3. |
| hbond_cdr3_target_pocket_count | hbond interactions for cdr3 with antigen residue in target pocket. |
| hbond_framework_count | hbond interactions for nanobody region framework. |
| hbond_framework_target_pocket_count | hbond interactions for framework with antigen residue in target pocket. |
| salt_bridge_count | Total salt_bridge interactions involving CDR regions. |
| salt_bridge_target_pocket_count | salt_bridge interactions where antigen residue is in target pocket. |
| salt_bridge_cdr1_count | salt_bridge interactions for nanobody region cdr1. |
| salt_bridge_cdr1_target_pocket_count | salt_bridge interactions for cdr1 with antigen residue in target pocket. |
| salt_bridge_cdr2_count | salt_bridge interactions for nanobody region cdr2. |
| salt_bridge_cdr2_target_pocket_count | salt_bridge interactions for cdr2 with antigen residue in target pocket. |
| salt_bridge_cdr3_count | salt_bridge interactions for nanobody region cdr3. |
| salt_bridge_cdr3_target_pocket_count | salt_bridge interactions for cdr3 with antigen residue in target pocket. |
| salt_bridge_framework_count | salt_bridge interactions for nanobody region framework. |
| salt_bridge_framework_target_pocket_count | salt_bridge interactions for framework with antigen residue in target pocket. |
| pi_pi_count | Total pi_pi interactions involving CDR regions. |
| pi_pi_target_pocket_count | pi_pi interactions where antigen residue is in target pocket. |
| pi_pi_cdr1_count | pi_pi interactions for nanobody region cdr1. |
| pi_pi_cdr1_target_pocket_count | pi_pi interactions for cdr1 with antigen residue in target pocket. |
| pi_pi_cdr2_count | pi_pi interactions for nanobody region cdr2. |
| pi_pi_cdr2_target_pocket_count | pi_pi interactions for cdr2 with antigen residue in target pocket. |
| pi_pi_cdr3_count | pi_pi interactions for nanobody region cdr3. |
| pi_pi_cdr3_target_pocket_count | pi_pi interactions for cdr3 with antigen residue in target pocket. |
| pi_pi_framework_count | pi_pi interactions for nanobody region framework. |
| pi_pi_framework_target_pocket_count | pi_pi interactions for framework with antigen residue in target pocket. |
| cation_pi_count | Total cation_pi interactions involving CDR regions. |
| cation_pi_target_pocket_count | cation_pi interactions where antigen residue is in target pocket. |
| cation_pi_cdr1_count | cation_pi interactions for nanobody region cdr1. |
| cation_pi_cdr1_target_pocket_count | cation_pi interactions for cdr1 with antigen residue in target pocket. |
| cation_pi_cdr2_count | cation_pi interactions for nanobody region cdr2. |
| cation_pi_cdr2_target_pocket_count | cation_pi interactions for cdr2 with antigen residue in target pocket. |
| cation_pi_cdr3_count | cation_pi interactions for nanobody region cdr3. |
| cation_pi_cdr3_target_pocket_count | cation_pi interactions for cdr3 with antigen residue in target pocket. |
| cation_pi_framework_count | cation_pi interactions for nanobody region framework. |
| cation_pi_framework_target_pocket_count | cation_pi interactions for framework with antigen residue in target pocket. |
| hydrophobic_count | Total hydrophobic interactions involving CDR regions. |
| hydrophobic_target_pocket_count | hydrophobic interactions where antigen residue is in target pocket. |
| hydrophobic_cdr1_count | hydrophobic interactions for nanobody region cdr1. |
| hydrophobic_cdr1_target_pocket_count | hydrophobic interactions for cdr1 with antigen residue in target pocket. |
| hydrophobic_cdr2_count | hydrophobic interactions for nanobody region cdr2. |
| hydrophobic_cdr2_target_pocket_count | hydrophobic interactions for cdr2 with antigen residue in target pocket. |
| hydrophobic_cdr3_count | hydrophobic interactions for nanobody region cdr3. |
| hydrophobic_cdr3_target_pocket_count | hydrophobic interactions for cdr3 with antigen residue in target pocket. |
| hydrophobic_framework_count | hydrophobic interactions for nanobody region framework. |
| hydrophobic_framework_target_pocket_count | hydrophobic interactions for framework with antigen residue in target pocket. |
| hbond_density_per_A2 | Interaction count normalized by dsasa_total_A2. |
| salt_bridge_density_per_A2 | Interaction count normalized by dsasa_total_A2. |
| pi_pi_density_per_A2 | Interaction count normalized by dsasa_total_A2. |
| cation_pi_density_per_A2 | Interaction count normalized by dsasa_total_A2. |
