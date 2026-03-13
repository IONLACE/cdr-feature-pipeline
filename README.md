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

Detailed column dictionary has been moved to:

- [docs/feature_columns.md](docs/feature_columns.md)
