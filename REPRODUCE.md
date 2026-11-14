# REPRODUCE

What of this release regenerates, what does not, and why. Every count is a claims-ledger row printed from its display, with the row id in square brackets.

## 0. Reproducibility levels

- Level one, inspect results: every released table, plotted frame and stored result is in this release with its size, md5, source file and producing script (`RELEASE_MANIFEST.tsv`).
- Level two, regenerate a figure: the figure code in `code/` runs from this release alone (released results, `ledger/claims_ledger_public.tsv`, `code/P31R_INPUT_MAP.tsv`). Only a unit of class `regenerates_from_release` is claimed at this level.
- Level three, recompute from upstream inputs: the analysis scripts are shipped, but the third-party inputs they read (provider summary statistics, reference-panel genotypes) are not redistributed. Each is pointed to by provider in the README; no unit is claimed at this level.

Level two (figure units): `Figure_1`, `Figure_2`, `FigureS1`, `FigureS2`, `FigureS6`, `FigureS6_p2`.
Level one only: `Figure_3`, `Figure_4`, `Figure_5`, `FigureS3`, `FigureS4`, `FigureS5`, `FigureS5_p2`, `annotation`, `colocalization`, `decoy`, `east_asian`, `external_outcomes`, `finemapping`, `kidney`, `literature`, `metabolite`, `prep_stamps`, `proxy`, `replication`, `screen`, `sensitivity_uncited`, `tables`, `workstream_A`, `workstreams`. Each is named with the missing input in section 2.

## 1. Regenerability, per unit

The class of each unit is derived by script from the approved source manifest (the table in section 1; the inputs are named per unit in section 2). Classes:

- `regenerates_from_release`: every input of the unit is in this release.
- `needs_public_download`: an input is a public download; the script is shipped, the input is not.
- `needs_nonredistributable_input`: an input is a third-party or individual-level file that is not redistributed, or an author file that is not in this release (each input line says which); the unit is not claimed as regenerable.
- `does_not_regenerate_engine_output`: an input is a stored engine output whose engine ran on the analysis machine; the unit is not claimed as regenerable.

| unit | kind | class | inputs not in the release, or released only as a copy the code does not read |
|---|---|---|---|
| Figure_1 | figure | regenerates_from_release | none |
| Figure_2 | figure | regenerates_from_release | none |
| Figure_3 | figure | needs_nonredistributable_input | `R2_rewrite/figures/data/Figure_3/umod_signal_ratios_se.tsv`; `V4_local/staging/l3-finemap/UMOD_eGFR_control/conditional/signals.tsv`; `V4_local/staging/l3-finemap/UMOD_eGFR_control/harmonised.tsv.gz` |
| Figure_4 | figure | needs_nonredistributable_input | `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz`; `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi`; `R1_stage1/S1R1_R2_12_chen_instrument/results.tsv` (released as a masked copy that the code does not read); `R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json`; `V4_local/raw/l3-finemap/ARG1_6q23__ref_A1dosage_503xP.i8`; `V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz` |
| Figure_5 | figure | needs_nonredistributable_input | `R2_rewrite/figures/data/Figure_5/eas_partition_members.tsv`; `V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv` |
| FigureS1 | figure | regenerates_from_release | none |
| FigureS2 | figure | regenerates_from_release | none |
| FigureS3 | figure | does_not_regenerate_engine_output | `paper3_ckd_rerun_v3/out/D1/kidney_cystatin_bun_concordance_E_current_20260922.csv` |
| FigureS4 | figure | needs_nonredistributable_input | `R2_rewrite/figures/data/Figure_3/FETCH_RECORD.json`; `R2_rewrite/figures/data/Figure_3/ckdgen_eGFR_EA_chr16_window_grch37.tsv`; `R2_rewrite/figures/data/Figure_3/ensembl_info_assembly.json`; `R2_rewrite/figures/data/Figure_3/ensembl_info_data.json`; `R2_rewrite/figures/data/Figure_3/ensembl_info_software.json`; `R2_rewrite/figures/data/Figure_3/ensembl_overlap.json`; `V4_local/raw/l3-finemap/DPEP1_16q24__ref_A1dosage_503xP.i8`; `V4_local/source_data/l3_metabolite/bridge_class_robustness__ens116pc.csv`; `V4_local/source_data/l3_metabolite/dpep1_locus_annotation.csv`; `V4_local/staging/l3-finemap/DPEP1_16q24/conditional/signals.tsv`; `V4_local/staging/l3-finemap/DPEP1_16q24/harmonised.tsv.gz`; `V4_local/staging/l3-finemap/DPEP1_16q24/traits.tsv` |
| FigureS5 | figure | needs_nonredistributable_input | `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz`; `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi`; `R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv`; `R2_rewrite/figures/data/Figure_4/FETCH_RECORD.json`; `R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json`; `R2_rewrite/figures/data/Figure_S5/arg1_signal_values_W4a.tsv`; `V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz` |
| FigureS5_p2 | figure | needs_nonredistributable_input | `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz`; `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi`; `R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv`; `R2_rewrite/figures/data/Figure_4/FETCH_RECORD.json`; `R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json`; `R2_rewrite/figures/data/Figure_S5/arg1_signal_values_W4a.tsv`; `V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz` |
| FigureS6 | figure | regenerates_from_release | none |
| FigureS6_p2 | figure | regenerates_from_release | none |
| annotation | deliverable | does_not_regenerate_engine_output | none |
| colocalization | deliverable | does_not_regenerate_engine_output | none |
| decoy | deliverable | does_not_regenerate_engine_output | none |
| east_asian | deliverable | does_not_regenerate_engine_output | none |
| external_outcomes | deliverable | does_not_regenerate_engine_output | none |
| finemapping | deliverable | does_not_regenerate_engine_output | none |
| kidney | deliverable | does_not_regenerate_engine_output | none |
| literature | deliverable | does_not_regenerate_engine_output | none |
| metabolite | deliverable | does_not_regenerate_engine_output | none |
| prep_stamps | deliverable | does_not_regenerate_engine_output | none |
| proxy | deliverable | does_not_regenerate_engine_output | none |
| replication | deliverable | does_not_regenerate_engine_output | none |
| screen | deliverable | does_not_regenerate_engine_output | `V4/source_data/tables_from_rds/RDS_CONVERSIONS.csv`; `paper3_deposit_request_v4/scripts/reports/Paper3_Program/run_logs/REGEN_20260817/01_mr_screen.log` |
| sensitivity_uncited | deliverable | does_not_regenerate_engine_output | none |
| tables | deliverable | needs_nonredistributable_input | `R2_rewrite/P31R_REFERENCE_KEYS.tsv`; `R2_rewrite/refs/new_refs_af1.tsv`; `R2_rewrite/refs/new_refs_supp1.tsv`; `R2_rewrite/refs/new_refs_supp2.tsv`; `R2_rewrite/refs/new_refs_table1.tsv`; `R2_rewrite/refs/new_refs_table2.tsv`; `R2_rewrite/scripts/P31R_build.py`; `R2_rewrite/scripts/R2_d1_display_order.tsv`; `R2_rewrite/sections_R2/P31R_Results_R2_draft.md`; `R2_rewrite/supplement/P31R_SM_cards.md`; `R2_rewrite/tables/P31R_Table1.md`; `R2_rewrite/tables/P31R_Table2.md`; `R2_rewrite/tables/P31R_Table2_provenance/devstage_opentargets_20261006.json`; `v3.2/prose/_refs_resolved.tsv` |
| workstream_A | deliverable | does_not_regenerate_engine_output | none |
| workstreams | deliverable | does_not_regenerate_engine_output | none |

A unit is claimed as regenerable only in the class `regenerates_from_release`. The other classes name the input that is missing; the unit is not claimed. For a deliverable unit (a stored result table), the class names the engine output that the released file is a copy of.

## 2. Inputs per unit

A regenerated PDF is compared with the render it reproduces (the figure files of the article and the pages of Additional file 3; this release ships no render) by its text and raster, not by md5: in four regenerations compared byte by byte (`P31R_FigS1.R`; the comparison is recorded in the project's review record, which is not shipped), the files had the same size as the render, and the bytes that differed were the creation date in the document information, the file identifier in the trailer and, in two of the four, the value of one integer object that nothing references.

Each source is shown with where it stands: `out:<path>` if the release holds it (or a file derived from it), else the reason it is not released.

Labels on the input lines: `verbatim copy`, the release holds the input byte for byte (its manifest md5 is the source md5); `masked copy`, the released file holds that input's rows with the cells that restate a provider's per-variant value, or that a data hold withholds, masked and the other cells as in the source (its md5 differs from the source md5); `reduced copy`, the release holds an extract of the input, the rows and columns the generator reads, with provider per-variant values left out; `code bundle; the shipped copy is ...`, a script of the code bundle, byte-identical to the original or derived from it by the portability transform the manifest's transform column names.

### Figure_1 (regenerates_from_release)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

### Figure_2 (regenerates_from_release)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_02_tgfb1/results.tsv -> out:workstreams/W2/results.tsv (masked copy)
  - input: R1_stage1/S1R1_R2_03_extcoloc/results.tsv -> out:workstreams/W3/results.tsv (masked copy)
  - input: R1_stage1/S1R1_R2_04b_arg1sets/results.tsv -> out:workstreams/W4b/results.tsv
  - input: R1_stage1/S1R1_R2_10_integration/corroboration.tsv -> out:workstreams/INT/corroboration.tsv (verbatim copy)
  - input: R1_stage1/S1R1_R2_11_extchecks/results.tsv -> out:external_outcomes/extchecks/results.tsv (masked copy)
  - input: R2_rewrite/review_R2/CLAIM_OUTLINE_STAGE3_2026-10-11_v6.md -> out:code/lift_inputs/Figure_2/CLAIM_OUTLINE_novelty_ledger.md (reduced copy)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

### Figure_3 (needs_nonredistributable_input)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_06a_renal/results.tsv -> out:workstreams/W6a/results.tsv
  - input: R1_stage1/S1R1_R2_06b_renal_d_umod/results.tsv -> out:workstreams/W6b/results.tsv
  - input: R1_stage1/S1R1_R2_1b_verify/S6/adjudication/A_third/a3_umod_q.tsv -> out:kidney/a3_umod_q.tsv
  - input: R2_rewrite/figures/data/Figure_3/umod_signal_ratios_se.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/UMOD_eGFR_control/conditional/signals.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/UMOD_eGFR_control/harmonised.tsv.gz -> non-redistributable (third_party_terms)

### Figure_4 (needs_nonredistributable_input)

  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz -> public download (public)
  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi -> public download (public)
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_04a_arg1grid/results.tsv -> out:workstreams/W4a/results.tsv
  - input: R1_stage1/S1R1_R2_04b_arg1sets/results.tsv -> out:workstreams/W4b/results.tsv
  - input: R1_stage1/S1R1_R2_10_integration/results.tsv -> out:workstreams/INT/results.tsv
  - input: R1_stage1/S1R1_R2_12_chen_instrument/credible_set_status.tsv -> out:workstream_A/credible_set_status.tsv
  - input: R1_stage1/S1R1_R2_12_chen_instrument/results.tsv -> out:workstream_A/results.tsv (masked copy; the figure code does not read it: `code/P31R_INPUT_MAP.tsv` marks the copy neither verbatim nor lifted, so the generator stops on this input)
  - input: R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json -> public download (public)
  - input: V4_local/raw/l3-finemap/ARG1_6q23__ref_A1dosage_503xP.i8 -> non-redistributable (individual_level)
  - input: V4_local/staging/l3-finemap/ARG1_6q23/PREP_STAMP.json -> out:source_data/prep_stamps/ARG1_6q23_PREP_STAMP.json
  - input: V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz -> non-redistributable (third_party_terms)

### Figure_5 (needs_nonredistributable_input)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R2_rewrite/figures/data/Figure_5/eas_partition_members.tsv -> author file not in this release (not_on_approved_list)
  - input: V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv -> out:screen/cis_pqtl_instruments.csv (masked copy)
  - input: V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv -> author file not in this release (not_on_approved_list)

### FigureS1 (regenerates_from_release)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R2_rewrite/posthoc/P31R_posthoc_results.tsv -> out:code/lift_inputs/FigureS1/P31R_posthoc_results.tsv (reduced copy)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

### FigureS2 (regenerates_from_release)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_10_integration/s6d_fitset_readings.tsv -> out:workstreams/INT/s6d_fitset_readings.tsv
  - input: V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv -> out:finemapping/SuppTable_candidate_edge_audit.csv

### FigureS3 (does_not_regenerate_engine_output)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: paper3_ckd_rerun_v3/out/D1/kidney_cystatin_bun_concordance_E_current_20260922.csv -> engine output (engine_output)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/pam_conditional/egfr_mirage_classification.csv -> out:kidney/egfr_mirage_classification.csv

### FigureS4 (needs_nonredistributable_input)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R2_rewrite/figures/data/Figure_3/FETCH_RECORD.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ckdgen_eGFR_EA_chr16_window_grch37.tsv -> non-redistributable (terms_unread_or_hold)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_info_assembly.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_info_data.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_info_software.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_overlap.json -> public download (public)
  - input: V4_local/raw/l3-finemap/DPEP1_16q24__ref_A1dosage_503xP.i8 -> non-redistributable (individual_level)
  - input: V4_local/source_data/l3_metabolite/bridge_class_robustness__ens116pc.csv -> author file not in this release (not_on_approved_list)
  - input: V4_local/source_data/l3_metabolite/dpep1_locus_annotation.csv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/PREP_STAMP.json -> out:source_data/prep_stamps/DPEP1_16q24_PREP_STAMP.json
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/conditional/signals.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/harmonised.tsv.gz -> non-redistributable (third_party_terms)
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/traits.tsv -> author file not in this release (not_on_approved_list)

### FigureS5 (needs_nonredistributable_input)

  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz -> public download (public)
  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi -> public download (public)
  - input: R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv -> author file not in this release (not_on_approved_list)
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_04a_arg1grid/results.tsv -> out:workstreams/W4a/results.tsv
  - input: R1_stage1/S1R1_R2_12_chen_instrument/credible_set_status.tsv -> out:workstream_A/credible_set_status.tsv
  - input: R2_rewrite/figures/data/Figure_4/FETCH_RECORD.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_S5/arg1_signal_values_W4a.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/d2-eas-census/results_census.tsv -> out:code/lift_inputs/FigureS6/results_census.tsv (reduced copy)
  - input: V4_local/staging/d2-eas-census/results_ckbfree.tsv -> out:code/lift_inputs/FigureS6/results_ckbfree.tsv (reduced copy)
  - input: V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz -> non-redistributable (third_party_terms)

### FigureS5_p2 (needs_nonredistributable_input)

  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz -> public download (public)
  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi -> public download (public)
  - input: R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv -> author file not in this release (not_on_approved_list)
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_04a_arg1grid/results.tsv -> out:workstreams/W4a/results.tsv
  - input: R1_stage1/S1R1_R2_12_chen_instrument/credible_set_status.tsv -> out:workstream_A/credible_set_status.tsv
  - input: R2_rewrite/figures/data/Figure_4/FETCH_RECORD.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_S5/arg1_signal_values_W4a.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/d2-eas-census/results_census.tsv -> out:code/lift_inputs/FigureS6/results_census.tsv (reduced copy)
  - input: V4_local/staging/d2-eas-census/results_ckbfree.tsv -> out:code/lift_inputs/FigureS6/results_ckbfree.tsv (reduced copy)
  - input: V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz -> non-redistributable (third_party_terms)

### FigureS6 (regenerates_from_release)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: V4/source_data/supp_tables/SuppTable_EAS_axis.csv -> out:code/lift_inputs/FigureS6/SuppTable_EAS_axis.csv (reduced copy)
  - input: V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv -> out:screen/cis_pqtl_instruments.csv (masked copy)
  - input: V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv -> out:screen/screen_ledger.csv (masked copy)
  - input: V4_local/staging/d2-eas-census/results_census.tsv -> out:code/lift_inputs/FigureS6/results_census.tsv (reduced copy)
  - input: V4_local/staging/d2-eas-census/results_ckbfree.tsv -> out:code/lift_inputs/FigureS6/results_ckbfree.tsv (reduced copy)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

### FigureS6_p2 (regenerates_from_release)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: V4/source_data/supp_tables/SuppTable_EAS_axis.csv -> out:code/lift_inputs/FigureS6/SuppTable_EAS_axis.csv (reduced copy)
  - input: V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv -> out:screen/cis_pqtl_instruments.csv (masked copy)
  - input: V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv -> out:screen/screen_ledger.csv (masked copy)
  - input: V4_local/staging/d2-eas-census/results_census.tsv -> out:code/lift_inputs/FigureS6/results_census.tsv (reduced copy)
  - input: V4_local/staging/d2-eas-census/results_ckbfree.tsv -> out:code/lift_inputs/FigureS6/results_ckbfree.tsv (reduced copy)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

### annotation (does_not_regenerate_engine_output)


### colocalization (does_not_regenerate_engine_output)


### decoy (does_not_regenerate_engine_output)


### east_asian (does_not_regenerate_engine_output)


### external_outcomes (does_not_regenerate_engine_output)


### finemapping (does_not_regenerate_engine_output)


### kidney (does_not_regenerate_engine_output)


### literature (does_not_regenerate_engine_output)


### metabolite (does_not_regenerate_engine_output)


### prep_stamps (does_not_regenerate_engine_output)


### proxy (does_not_regenerate_engine_output)


### replication (does_not_regenerate_engine_output)


### screen (does_not_regenerate_engine_output)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: V4/source_data/tables_from_rds/RDS_CONVERSIONS.csv -> engine output (engine_output)
  - input: V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv -> out:screen/cis_pqtl_instruments.csv (masked copy)
  - input: V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv -> out:screen/screen_ledger.csv (masked copy)
  - input: paper3_deposit_request_v4/scripts/R/00_config.R -> out:code/paper3_deposit_request_v4/scripts/R/00_config.R (code bundle; the shipped copy is derived from the original by the declared portability transform)
  - input: paper3_deposit_request_v4/scripts/reports/Paper3_Program/run_logs/REGEN_20260817/01_mr_screen.log -> author file not in this release (code_bundle_step3)
  - input: paper3_deposit_request_v4/scripts/scripts/paper3.1_ckm/cascade/01_mr_screen.R -> out:code/paper3_deposit_request_v4/scripts/scripts/paper3.1_ckm/cascade/01_mr_screen.R (code bundle; the shipped copy is derived from the original by the declared portability transform)

### sensitivity_uncited (does_not_regenerate_engine_output)


### tables (needs_nonredistributable_input)

  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_06_claims_ledger/ext_table2.py -> out:code/R1_stage1/S1R1_06_claims_ledger/ext_table2.py (code bundle; the shipped copy is byte-identical to the original)
  - input: R2_rewrite/P31R_REFERENCE_KEYS.tsv -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/refs/new_refs_af1.tsv -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/refs/new_refs_supp1.tsv -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/refs/new_refs_supp2.tsv -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/refs/new_refs_table1.tsv -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/refs/new_refs_table2.tsv -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/scripts/P31R_build.py -> author file not in this release (code_bundle_step3)
  - input: R2_rewrite/scripts/P31R_gate.py -> out:code/R2_rewrite/scripts/P31R_gate.py (code bundle; the shipped copy is byte-identical to the original)
  - input: R2_rewrite/scripts/P31R_tablegate.py -> out:code/R2_rewrite/scripts/P31R_tablegate.py (code bundle; the shipped copy is byte-identical to the original)
  - input: R2_rewrite/scripts/R2_d1_display_order.tsv -> author file not in this release (code_bundle_step3)
  - input: R2_rewrite/sections_R2/P31R_Results_R2_draft.md -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/supplement/P31R_SM_cards.md -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/tables/P31R_Table1.md -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/tables/P31R_Table2.md -> author file not in this release (not_on_approved_list)
  - input: R2_rewrite/tables/P31R_Table2_provenance/devstage_opentargets_20261006.json -> public download (public)
  - input: v3.2/prose/_refs_resolved.tsv -> author file not in this release (not_on_approved_list)

### workstream_A (does_not_regenerate_engine_output)


### workstreams (does_not_regenerate_engine_output)


## 3. Commands for the code the bundle ships

Run every command from the archive root (the folder that holds `code/`, `ledger/` and the released tables). The shipped copies of the scripts are portable: the scripts that read the archive configuration resolve their paths from one configuration, the environment variable `P31R_ARCHIVE` (default: the archive root, derived from the script's own place) and `P31R_EXTERNAL` (default: `external_inputs/` beside it, for third-party inputs that are not in the release). `code/P31R_INPUT_MAP.tsv` maps each input's project path to its path in the archive and says whether the released copy is verbatim. `code/P31R_PORTABILITY_LOG.tsv` maps the source md5 of each shipped script, which is what the md5 pins inside the scripts name, to the shipped md5. The figure scripts stop with `input not in the release: <path>` and the reason when asked for an input that is not in the release. The analysis scripts shipped for tracing (the workstream and decoy scripts, the ledger builder) have the same path roots but fail with an ordinary file-not-found error on an input that is not in the release.

`P31R_Fig1.R`, `P31R_Fig2.R`, `P31R_Fig3.R` do not read `P31R_ARCHIVE`: each finds its own folder from the command line and sources `P31R_fig_common.R` from that folder, which derives the figures folder from its own place.

A path printed `<external>/<path>` is a path on a drive of the analysis machine with the drive letter removed. The portable scripts open it below `P31R_EXTERNAL` under a folder named for that drive letter, as their text shows: `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz` is read from `external_inputs/D/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz`; `<external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi` is read from `external_inputs/D/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi`.

A tab-separated frame that holds a line break inside a cell is shipped rewritten (declared transform of the release builder, `normalise_newlines`): each break inside a cell is replaced by ` / `, the quotes that this makes unnecessary are dropped, and the rows end with a line feed alone. Frames so rewritten, whose every row sits on one line: `source_data/figures/Figure_2/fig2b_plotted.tsv`. A fresh render of such a frame differs from the shipped file in those three ways, so it differs in md5; compare the two by cells with a CSV reader, after replacing each break inside a cell by ` / `. The other frames are shipped as rendered: most comma-separated frames end their rows with a carriage return and a line feed, and the ledger and the rewritten frames with a line feed.
Frames shipped as rendered that hold line breaks inside quoted cells (so a row can span lines; read them with a CSV reader): `source_data/figures/FigureS5/FigureS5__per_signal_estimates.csv`, `source_data/figures/FigureS5_p2/FigureS5_p2__credible_set_grid.csv`.

The commands run `python -I` and `Rscript --vanilla`. A command is listed for every active figure unit; only a unit of class `regenerates_from_release` is claimed as regenerable. The figure code reads its printed numbers from `ledger/claims_ledger_public.tsv`; a figure whose numbers include a masked row stops with the reason.

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_Fig1.R
```
Unit(s): Figure_1. Class: regenerates_from_release.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_Fig2.R
```
Unit(s): Figure_2. Class: regenerates_from_release.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_02_tgfb1/results.tsv -> out:workstreams/W2/results.tsv (masked copy)
  - input: R1_stage1/S1R1_R2_03_extcoloc/results.tsv -> out:workstreams/W3/results.tsv (masked copy)
  - input: R1_stage1/S1R1_R2_04b_arg1sets/results.tsv -> out:workstreams/W4b/results.tsv
  - input: R1_stage1/S1R1_R2_10_integration/corroboration.tsv -> out:workstreams/INT/corroboration.tsv (verbatim copy)
  - input: R1_stage1/S1R1_R2_11_extchecks/results.tsv -> out:external_outcomes/extchecks/results.tsv (masked copy)
  - input: R2_rewrite/review_R2/CLAIM_OUTLINE_STAGE3_2026-10-11_v6.md -> out:code/lift_inputs/Figure_2/CLAIM_OUTLINE_novelty_ledger.md (reduced copy)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_Fig3.R
```
Unit(s): Figure_3. Class: needs_nonredistributable_input.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_06a_renal/results.tsv -> out:workstreams/W6a/results.tsv
  - input: R1_stage1/S1R1_R2_06b_renal_d_umod/results.tsv -> out:workstreams/W6b/results.tsv
  - input: R1_stage1/S1R1_R2_1b_verify/S6/adjudication/A_third/a3_umod_q.tsv -> out:kidney/a3_umod_q.tsv
  - input: R2_rewrite/figures/data/Figure_3/umod_signal_ratios_se.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/UMOD_eGFR_control/conditional/signals.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/UMOD_eGFR_control/harmonised.tsv.gz -> non-redistributable (third_party_terms)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_Fig4.R
```
Unit(s): Figure_4. Class: needs_nonredistributable_input.
  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz -> public download (public)
  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi -> public download (public)
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_04a_arg1grid/results.tsv -> out:workstreams/W4a/results.tsv
  - input: R1_stage1/S1R1_R2_04b_arg1sets/results.tsv -> out:workstreams/W4b/results.tsv
  - input: R1_stage1/S1R1_R2_10_integration/results.tsv -> out:workstreams/INT/results.tsv
  - input: R1_stage1/S1R1_R2_12_chen_instrument/credible_set_status.tsv -> out:workstream_A/credible_set_status.tsv
  - input: R1_stage1/S1R1_R2_12_chen_instrument/results.tsv -> out:workstream_A/results.tsv (masked copy; the figure code does not read it: `code/P31R_INPUT_MAP.tsv` marks the copy neither verbatim nor lifted, so the generator stops on this input)
  - input: R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json -> public download (public)
  - input: V4_local/raw/l3-finemap/ARG1_6q23__ref_A1dosage_503xP.i8 -> non-redistributable (individual_level)
  - input: V4_local/staging/l3-finemap/ARG1_6q23/PREP_STAMP.json -> out:source_data/prep_stamps/ARG1_6q23_PREP_STAMP.json
  - input: V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz -> non-redistributable (third_party_terms)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_Fig5.R
```
Unit(s): Figure_5. Class: needs_nonredistributable_input.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R2_rewrite/figures/data/Figure_5/eas_partition_members.tsv -> author file not in this release (not_on_approved_list)
  - input: V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv -> out:screen/cis_pqtl_instruments.csv (masked copy)
  - input: V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv -> author file not in this release (not_on_approved_list)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_FigS1.R
```
Unit(s): FigureS1. Class: regenerates_from_release.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R2_rewrite/posthoc/P31R_posthoc_results.tsv -> out:code/lift_inputs/FigureS1/P31R_posthoc_results.tsv (reduced copy)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_FigS2.R
```
Unit(s): FigureS2. Class: regenerates_from_release.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_10_integration/s6d_fitset_readings.tsv -> out:workstreams/INT/s6d_fitset_readings.tsv
  - input: V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv -> out:finemapping/SuppTable_candidate_edge_audit.csv

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_FigS3.R
```
Unit(s): FigureS3. Class: does_not_regenerate_engine_output.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: paper3_ckd_rerun_v3/out/D1/kidney_cystatin_bun_concordance_E_current_20260922.csv -> engine output (engine_output)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/pam_conditional/egfr_mirage_classification.csv -> out:kidney/egfr_mirage_classification.csv

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_FigS4.R
```
Unit(s): FigureS4. Class: needs_nonredistributable_input.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R2_rewrite/figures/data/Figure_3/FETCH_RECORD.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ckdgen_eGFR_EA_chr16_window_grch37.tsv -> non-redistributable (terms_unread_or_hold)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_info_assembly.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_info_data.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_info_software.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_3/ensembl_overlap.json -> public download (public)
  - input: V4_local/raw/l3-finemap/DPEP1_16q24__ref_A1dosage_503xP.i8 -> non-redistributable (individual_level)
  - input: V4_local/source_data/l3_metabolite/bridge_class_robustness__ens116pc.csv -> author file not in this release (not_on_approved_list)
  - input: V4_local/source_data/l3_metabolite/dpep1_locus_annotation.csv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/PREP_STAMP.json -> out:source_data/prep_stamps/DPEP1_16q24_PREP_STAMP.json
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/conditional/signals.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/harmonised.tsv.gz -> non-redistributable (third_party_terms)
  - input: V4_local/staging/l3-finemap/DPEP1_16q24/traits.tsv -> author file not in this release (not_on_approved_list)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_FigS5.R
```
Unit(s): FigureS5, FigureS5_p2. Class: needs_nonredistributable_input.
  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz -> public download (public)
  - input: <external>/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz.tbi -> public download (public)
  - input: R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv -> author file not in this release (not_on_approved_list)
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: R1_stage1/S1R1_R2_04a_arg1grid/results.tsv -> out:workstreams/W4a/results.tsv
  - input: R1_stage1/S1R1_R2_12_chen_instrument/credible_set_status.tsv -> out:workstream_A/credible_set_status.tsv
  - input: R2_rewrite/figures/data/Figure_4/FETCH_RECORD.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_4/ensembl_overlap_gene_transcript_exon.json -> public download (public)
  - input: R2_rewrite/figures/data/Figure_S5/arg1_signal_values_W4a.tsv -> author file not in this release (not_on_approved_list)
  - input: V4_local/staging/d2-eas-census/results_census.tsv -> out:code/lift_inputs/FigureS6/results_census.tsv (reduced copy)
  - input: V4_local/staging/d2-eas-census/results_ckbfree.tsv -> out:code/lift_inputs/FigureS6/results_ckbfree.tsv (reduced copy)
  - input: V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz -> non-redistributable (third_party_terms)

```
Rscript --vanilla code/R2_rewrite/figures/scripts/P31R_FigS6.R
```
Unit(s): FigureS6, FigureS6_p2. Class: regenerates_from_release.
  - input: R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv -> out:ledger/claims_ledger_public.tsv (public extract; masked rows carry their reason)
  - input: V4/source_data/supp_tables/SuppTable_EAS_axis.csv -> out:code/lift_inputs/FigureS6/SuppTable_EAS_axis.csv (reduced copy)
  - input: V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv -> out:screen/cis_pqtl_instruments.csv (masked copy)
  - input: V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv -> out:screen/screen_ledger.csv (masked copy)
  - input: V4_local/staging/d2-eas-census/results_census.tsv -> out:code/lift_inputs/FigureS6/results_census.tsv (reduced copy)
  - input: V4_local/staging/d2-eas-census/results_ckbfree.tsv -> out:code/lift_inputs/FigureS6/results_ckbfree.tsv (reduced copy)
  - input: v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv -> out:colocalization/table_s_pair_evidence_tiered.csv (masked copy)

Aggregators (shipped for inspection, not classified by the regenerability table and not claimed as regenerable; they read the claims ledger, of which the release holds the public extract only, and the manuscript sections, which this release does not contain; the Additional file 1 builder also imports the manuscript builder `P31R_build.py`, which is not shipped):

- `python -I -X utf8 code/R2_rewrite/figures/af3/P31R_af3_assemble.py`: assembles Additional file 3.
- `python -I -X utf8 code/R2_rewrite/figures/af3/build_legend_AF3_R2.py`: builds the Additional file 3 legend file.
- `python -I -X utf8 code/R2_rewrite/scripts/P31R_build_af1.py`: builds Additional file 1.
- `python -I -X utf8 code/R2_rewrite/figures/render_all_units.py`: runs every active generator once.

Of the aggregators, `build_legend_AF3_R2.py`, `P31R_build_af1.py`, `render_all_units.py` do not read `P31R_ARCHIVE` either: each derives the project root from its own place in the tree.


The release builder (`code/R2_rewrite/release/`) reads the project tree and the withheld files; it is shipped for inspection and is not rerunnable from this repository alone.

## 4. What does not regenerate

- The screen, colocalization, SuSiE and candidate engines ran on the analysis machine. Their stored outputs are released; their scripts are only partly deposited. Stored screen rows: 14,537 [L660].
- Deposit v4 delivered 48 files (44 in scripts/scripts, 4 in scripts/R) [L664]. `release_records/CODE_INCLUDE.tsv` accounts for those files as follows (the counts are computed from that file when this document is built): `47` rows of role `deposit_v4_script`, of which `11` ship (`5` that the v5 delivery re-confirmed against the analysis machine, and `6` R scripts of the primary colocalization run that it did not re-confirm, each with the provenance label of its evidence) and `36` do not (ship = 0, each with a reason, pending a decision); and `01_renal_triad_mr.R`, which ships as reference code under its `paper3_ckd_rerun_v3` path without a v5 re-confirmation, with the label `timestamp evidence only`.
- The v5 delivery (Part A of the request that asked for these scripts, `REQUEST_PROMPT_v5.md`, which is not shipped) supplied the SuSiE, cascade, candidate-chain and East Asian axis scripts. They are in the code bundle with the provenance labels of the README; the labels say how strongly each file is tied to the run that wrote the results, and what they do not establish.
- `l3fm__conditional.py`, the conditioning engine named in Supplementary Sections S-M6 and S-M7, is not in deposit v4. The copy in the project (`V4_local/builders/`) is in the code bundle with its helpers; that it is the version that produced the results is unverified.
- The Stage 1 workstream scripts and the ledger builder with its extension modules are in the code bundle, so that each result can be traced to its script. The ledger builder reads the whole project tree and is not rerunnable from this release.
- The per-cell screen audit files and the full proxy cache stay on the analysis machine; they are not in this release. Rows of the partial proxy record shipped here: 5 [L668].
- The decoy-rate script (`decoy/S1R1_decoy_rates.R`) is released; its inputs are not. Inputs kept out: 6 [L669].
- Individual-level dosage matrices (`.i8`) found: 12 [L665]; none is released.
- Screen counts. `screen/screen_ledger.csv` has one row for each protein-outcome combination of the multiple-testing family, 15,584 [L003]: the tests with an estimate, 14,521 [L004], and the untestable combinations (1,061 with no outcome data; 2 attempted with no estimate [L674]) (the estimates and the untestable combinations add up to the family; checked when this document is built). The stored screen output has 14,537 [L660]: every estimator's rows for the tested pairs, which is the primary row of each tested pair plus the non-primary rows shipped as `screen/screen_nonprimary_estimates.csv` (checked: the rows of that file equal the stored rows minus the tests with an estimate). The coverage statement (86 of 1,948 instrument genes have no stored screen row (the stored screen holds 1,862 genes) [L675]) refers to that stored output, not to `screen/screen_ledger.csv`.

## 5. Provenance caveats

- 3 writer scripts edited after their output was written; 1 hand-relabelled input (A9) [L666] (deposit v4 report, Part B; the report is not part of this release).
- Licence classes are decided per file by content (README, "Licence tiers"). The release does not redistribute any provider's summary-statistics file or any individual-level data. Individual per-variant provider values do appear, in four ways: (a) where the provider's licence permits; (b) where a single variant defines a per-pair estimate or the manuscript itself prints the provider's value at the analysed variant (for the FGFR4-T2D pair, TPMI's effect and standard error at the instrument follow exactly from the shipped Wald ratio and the shipped CKB exposure effect, and TPMI's own p at the instrument ships, as printed in Table 2 and Figure 2b; the CLSA arginine marginal p at the two ARG1 variants ships, as printed in the Results); (c) as values plotted in a figure (some binned or rounded); and (d) quoted in shipped code or audit notes (listed in `release_records/LEAKSCAN_COINCIDENCES.tsv`, by file, line and provider, with the reason `quoted_value`; the values are not repeated there). Separately, Additional file 1 of the article (Table S3) prints, as instrument reporting, the UKB-PPP lead-variant effect, standard error, -log10(p), F statistic and allele frequency of each colocalization-supported pair; this release withholds the same columns in `screen/cis_pqtl_instruments.csv`, which covers every instrument and not only those pairs.

## 6. Environment

Versions as recorded: screen: R 4.5.1, TwoSampleMR 0.7.1, data.table 1.17.8, coloc 5.2.3; coloc log 5.2.3; fine-mapping R 4.4.1, susieR 0.14.2, coloc 5.2.3; susieR 0.14.2 (environment record, after the runs); PLINK v1.90b7.2 (SuSiE LD); Supplementary Note R 4.4.1; ancestry-matched sensitivity coloc 5.2.3 [L447].
Level-two figure rendering (the regenerating units of section 0) needs, as installed on the machine that built this release: R version 4.4.1 (2024-06-14 ucrt); the R packages `GenomicRanges` 1.58.0, `IRanges` 2.40.1, `Rsamtools` 2.22.0, `data.table` 1.15.4, `ggplot2` 4.0.3, `ggrepel` 0.9.8, `ggtext` 0.1.2, `jsonlite` 2.0.0, `patchwork` 1.3.2, `pdftools` 3.8.0, `png` 0.1.9, `scales` 1.4.0, `systemfonts` 1.1.0 (read from the shipped figure code); the base packages `grDevices`, `graphics`, `grid`, `tools`, `utils`; and Python 3.12.4 on the path as `python` with PyMuPDF (`fitz`) 1.27.2.3. `code/R2_rewrite/figures/scripts/P31R_fig_common.R` calls `python` (without `-I`) to write the author metadata into each PDF before the PNG is rasterized from it; without it a generator stops after it has written the frame and the PDF.
The panels of Figure 4 and FigureS5 that read genotype windows also load `Rsamtools`, `GenomicRanges` and `IRanges` and need inputs that are not in this release (section 2).
`software_label` in `RELEASE_MANIFEST.tsv`: `run-time` (the version was written by the run itself), `post hoc` (recorded afterwards in an environment record), `none_recorded` (no version of the producing run is recorded in the claims ledger). A version recorded after the run is not a run-time version.
