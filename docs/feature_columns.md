# Feature Columns

This document contains the detailed column dictionary for `custom_features_merged_all.csv`.

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
