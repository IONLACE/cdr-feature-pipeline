```bash

python -m src.01-build_msa 8q6k_C_BA_Repair.pdb 

python -m src.02-calc_conservation 8q6k_C_BA_Repair.pdb 
    Input PDB file: /Users/satyan/Wrk/zff_stuff/input_pdbs/8q6k/8q6k_C_BA_Repair.pdb, Antigen chain ID: BA

    ▶ Running command:
    /Users/satyan/Wrk/packages/cons-capra07-main/target/release/cons-capra07 -i /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/msa/ag_B_ncbi_clean.fasta -o /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/msa/ag_B_cons_scores.tsv -w va -c yes -t yes
    ✔ Success
    Parsed 227 conservation scores

    ▶ Running command:
    /Users/satyan/Wrk/packages/cons-capra07-main/target/release/cons-capra07 -i /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/msa/ag_A_ncbi_clean.fasta -o /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/msa/ag_A_cons_scores.tsv -w va -c yes -t yes
    ✔ Success
    Parsed 212 conservation scores
    Combined residue-wise scores written to /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/8q6k_C_BA_Repair_cons_mapped.csv


python -m src.03-get_CDRs 8q6k_C_BA_Repair.pdb 
    Input PDB file: /Users/satyan/Wrk/zff_stuff/input_pdbs/8q6k/8q6k_C_BA_Repair.pdb, chain ID: C
    Chain type: 1
    Species: 128
    CDR sequences (IMGT):
    CDR1: GFTFDDYA (len=8)
    CDR2: IRRSDGST (len=8)
    CDR3: AAAGTPSYYYTEPLSLGTYDY (len=21)
    CSV written to: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/cdrs/8q6k_C_BA_Repair_cdrs.csv


python -m src.04-run_fpocket 8q6k_C_BA_Repair.pdb
    Input PDB file: /Users/satyan/Wrk/zff_stuff/input_pdbs/8q6k/8q6k_C_BA_Repair.pdb, chain ID dropped: C

    ▶ Running command:
    /opt/anaconda3/envs/fpocket/bin/fpocket -f /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/pocket/8q6k_C_BA_Repair.pdb -c C
    ✔ Success
    Pocket detection completed for 8q6k. Output in: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/pocket
    Parsed pocket features saved to /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/pocket_features.json
    Merged pocket residues written to /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/pocket_residues.csv

python -m src.05-get_all_contacts 8q6k_C_BA_Repair.pdb
    Input PDB file: /Users/satyan/Wrk/zff_stuff/input_pdbs/8q6k/8q6k_C_BA_Repair.pdb, nanobody chain: C, antigen chain(s): B,A
    [INFO] Unified contacts+CDR mapping saved to /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/contacts/8q6k_C_BA_Repair_all_and_CDR_contacts.csv
    Wrote /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/contacts/8q6k_C_BA_Repair_CDR_pocket_contacts.csv

python -m src.06-run_pandaprot 8q6k_C_BA_Repair.pdb
    Input PDB file: /Users/satyan/Wrk/zff_stuff/input_pdbs/8q6k/8q6k_C_BA_Repair.pdb, chain IDs: B, A, C
    PandaProt run completed. Logs and outputs are in: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/pandaprot_run
    Filtered report saved to /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/pandaprot_run/pandaprot_report.csv (196 -> 29 rows)

python -m src.07-interface_complimentarity 8q6k_C_BA_Repair.pdb
    Results saved to: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/interface/interface_Sc_summary.txt
    ✓ Analysis complete!
    Results saved to: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/interface/
    Report log: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/interface/interface_Sc_report.log
    Summary report: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/interface/interface_Sc_summary.txt

    ======================================================================
    ANALYSIS RESULTS
    ======================================================================
    Lawrence & Colman Sc: 0.4472
    Overall score: 0.5131
    Quality: Good

    Flat CSV saved to: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/interface/8q6k_C_BA_Repair_interface_results.csv

python -m src.08-compute_sasa 8q6k_C_BA_Repair.pdb 
    Input PDB file: /Users/satyan/Wrk/zff_stuff/input_pdbs/8q6k/8q6k_C_BA_Repair.pdb, nanobody chain: C, antigen chain: BA
    dSASA results written to /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/sasa/dsasa.csv


python -m src.09-CDR_peptide_features 8q6k_C_BA_Repair.pdb
    /opt/anaconda3/envs/testscore/lib/python3.12/site-packages/propy/AAIndex.py:31: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
    import pkg_resources
    Reading CDR sequences from: /Users/satyan/Wrk/zff_stuff/meta/8q6k_C_BA_Repair/cdrs/8q6k_C_BA_Repair_cdrs.csv
    Wrote /Users/satyan/Wrk/zff_stuff/features/8q6k_C_BA_Repair/cdrs/8q6k_C_BA_Repair_cdr_peptide_features.csv


python -m src.build_custom_nb_features  8q6k_C_BA_Repair.pdb   
    Wrote /Users/satyan/Wrk/zff_stuff/features/8q6k_C_BA_Repair_custom_features.csv

```
