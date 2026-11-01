# Release repository for "A genetic evidence map of circulating proteins across cardiovascular, kidney and metabolic outcomes", version `v2.0.5-R2`

Result tables, per-figure plotted frames, tables as data and analysis-tail code for the manuscript. `REPRODUCE.md` says what regenerates and what does not. `RELEASE_MANIFEST.tsv` lists every released file with its size, md5, row and column counts, source file and md5, producing script and md5, software label, citing document and licence class. The manifest cannot hold its own md5, so it and the files `README.md`, `REPRODUCE.md`, `LICENSE`, `CITATION.cff` and `.gitattributes` have no manifest row. `.gitattributes` holds the one line `* -text`, so that git never rewrites line endings and a clone equals the zip byte for byte.

## Version, licence and citation

This is release `v2.0.5-R2`. `LICENSE` gives the terms: the MIT License for code the authors wrote, and `CC BY 4.0` for tables the authors derived (`own_derived` in `RELEASE_MANIFEST.tsv`); files derived from third-party data keep their providers' terms (below). `CITATION.cff` says how to cite the release. The repository is `https://github.com/badboybert/ckm-protein-evidence-map`; its archived versions are on Zenodo under the concept DOI `https://doi.org/10.5281/zenodo.23855497`, which resolves to the latest version.

On Windows, extracting the zip or cloning the repository into a deep folder can fail once a full path exceeds the default Windows path limit (`MAX_PATH`). The longest paths inside this release are `P31R_R2_release_v2.0.5-R2/code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/01_acquire_and_localize_outcomes.py` and `P31R_R2_release_v2.0.5-R2/code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/31_candidate_susie_sensitivity.R` (for a clone, the repository folder name replaces the zip root `P31R_R2_release_v2.0.5-R2/`), so extract or clone close to the drive root, or run `git config core.longpaths true`.

## Licence tiers, third-party data and the TPMI switch

The release ships the authors' derived results and code. Third-party data are pointed to (provider, accession or URL, citation), not redistributed. Each file has a tier, decided by its content:

- First tier, derived per-pair results (MR estimates, colocalization posteriors, counts, classes): ship, with a citation.
- Second tier, per-variant values copied from a provider: ship only where redistribution is confirmed (CC0 or CC BY); otherwise the column or file is withheld, or reduced to the values plotted (regional tracks are binned). The exceptions that remain are listed in the paragraph after these tiers.
- Third tier, no licence (TPMI): the derived per-pair values ship with a citation; the TPMI effect, standard error and allele frequency columns do not.

The release does not redistribute any provider's summary-statistics file or any individual-level data. Individual per-variant provider values do appear, in four ways: (a) where the provider's licence permits; (b) where a single variant defines a per-pair estimate or the manuscript itself prints the provider's value at the analysed variant (for the FGFR4-T2D pair, TPMI's effect and standard error at the instrument follow exactly from the shipped Wald ratio and the shipped CKB exposure effect, and TPMI's own p at the instrument ships, as printed in Table 2 and Figure 2b; the CLSA arginine marginal p at the two ARG1 variants ships, as printed in the Results); (c) as values plotted in a figure (some binned or rounded); and (d) quoted in shipped code or audit notes (listed in `release_records/LEAKSCAN_COINCIDENCES.tsv`, by file, line and provider, with the reason `quoted_value`; the values are not repeated there). Separately, Additional file 1 of the article (Table S3) prints, as instrument reporting, the UKB-PPP lead-variant effect, standard error, -log10(p), F statistic and allele frequency of each colocalization-supported pair; this release withholds the same columns in `screen/cis_pqtl_instruments.csv`, which covers every instrument and not only those pairs.

TPMI is one declared switch (`release_records/TIER_RULES.tsv`, row `SW1`), set to `ship` for this release. The derived per-pair TPMI values ship with a citation of the TPMI PheWeb (`PheCode 250.2`): the Wald ratio, its standard error, p, class and verdict, and the study-level case and control counts. The TPMI effect, standard error and allele frequency columns of the check table are withheld. A single-instrument Wald ratio with its standard error is a derived per-pair result (first tier) and ships. ChEMBL- and DGIdb-derived drug and stage data fields are withheld (ChEMBL is CC BY-SA); the mask covers these data fields only, and drug names in manuscript prose come from the literature. Shipped scripts that quote ChEMBL facts as string literals are declared under Code licence. Open Targets (CC0) fields are native and ship.

Where each third-party source is obtained, and how to cite it. A source whose redistribution is not confirmed contributes derived results; the cases in which an individual per-variant value appears are listed above. A licence is given only where the provider states one. The page that was read for each provider, the wording it states and its status are recorded in `release_records/PROVIDER_TERMS.tsv`. Names, routes and citations are identifiers and are set in code font.

The provider list holds 35 providers, 14 named in the Declarations acknowledgements [L533].

| provider | accession or URL | citation |
|---|---|---|
| `UK Biobank Pharma Proteomics Project` | `Synapse (Sage Bionetworks) project wiki; AWS Open Data registry mirror` | `UK Biobank Pharma Proteomics Project (Sun et al. 2023, Nature 622:329)` |
| `deCODE genetics (Icelandic plasma proteome, SomaScan)` | `decode.com/summarydata download form` | `deCODE genetics, Icelandic SomaScan plasma proteome (Ferkingstad et al. 2021)` |
| `Atherosclerosis Risk in Communities study (SomaScan)` | `nilanjanchatterjeelab.org/pwas (Zhang et al. 2022 data statement)` | `Atherosclerosis Risk in Communities study SomaScan pQTL (Zhang et al. 2022)` |
| `FinnGen R10` | `finngen.fi access form; Google Cloud Storage download` | `FinnGen` |
| `MEGASTROKE consortium` | `megastroke.org download; GCST006906` | `MEGASTROKE consortium (Malik et al. 2018)` |
| `GIGASTROKE consortium` | `GWAS Catalog (GCST90104544 East Asian arm; GCST90104539 European arm used as the East Asian screen's stroke comparator)` | `GIGASTROKE consortium (Mishra et al. 2022); licence CC0 1.0` |
| `CKDGen consortium` | `ckdgen.imbi.uni-freiburg.de Data Sets page (University of Freiburg)` | `CKDGen consortium (Wuttke et al. 2019)` |
| `DIAMANTE consortium (T2D, European, without UK Biobank)` | `diagram-consortium.org downloads (click-through)` | `DIAMANTE consortium (Mahajan et al. 2018 Nat Genet)` |
| `CARDIoGRAMplusC4D consortium` | `cardiogramplusc4d.org data-downloads; IEU OpenGWAS ieu-a-7` | `CARDIoGRAMplusC4D consortium (Nikpay et al. 2015)` |
| `Taiwan Precision Medicine Initiative (PheWeb)` | `TPMI PheWeb, pheweb.ibms.sinica.edu.tw; PheCode 250.2` | `Taiwan Precision Medicine Initiative (Academia Sinica)` |
| `China Kadoorie Biobank (Olink)` | `CKB PheWeb pheweb.ckbiobank.org; encrypted downloads via ckbaccess@ndph.ox.ac.uk` | `China Kadoorie Biobank (Pozarickij et al.; Walters et al. 2023); licence CC BY 4.0` |
| `Biobank Japan (hum0014.v29 AF; hum0197.v3 CKD, ischemic stroke, HF and T2D releases)` | `NBDC Human Database hum0014 / hum0197; GWAS Catalog; pheweb.jp; GCST90018602; hum0014.v29; hum0197.v3` | `Biobank Japan (Sakaue 2021; Ishigaki 2020; others)` |
| `AGEN consortium (T2D, East Asian)` | `GCST010118 (GWAS Catalog); blog.nus.edu.sg/agen` | `AGEN consortium (Spracklen et al. 2020 Nature)` |
| `Canadian Longitudinal Study on Aging (metabolite GWAS)` | `GWAS Catalog GCST90199621-GCST90204063 (Chen et al. 2023, PMID 36635386); individual-level data via clsa-elcv.ca` | `Canadian Longitudinal Study on Aging, metabolite GWAS (Chen et al. 2023)` |
| `Taiwan-Japan eGFR meta-analysis (Taiwan Biobank and Biobank Japan cohorts; TWB2_BBJ_eGFR, PMID 39472450)` | `Figshare 10.6084/m9.figshare.24356587 (article data statement); BBJ part at jenger.riken.jp` | `Taiwan-Japan eGFR meta-analysis (PMID 39472450)` |
| `Type 2 Diabetes Global Genomics Initiative (East Asian T2D)` | `diagram-consortium.org downloads (click-through); GWAS Catalog GCST90492734 (summary statistics not downloadable from the Catalog: full_summary_stats_available = false)` | `Type 2 Diabetes Global Genomics Initiative (Suzuki et al. 2024 Nature)` |
| `1000 Genomes Project (European and East Asian samples, LD reference)` | `internationalgenome.org (IGSR)` | `1000 Genomes Project / IGSR` |
| `Open Targets Platform (development-stage snapshot of 2026-10-06)` | `platform-docs.opentargets.org; Open Targets downloads` | `Open Targets Platform; licence CC0 1.0` |
| `Ensembl (release 116 gene annotation and REST overlap)` | `ensembl.org; REST` | `Ensembl (EMBL-EBI / Wellcome Sanger)` |
| `East Asian stratum (N 282,852) of the multi-ancestry eGFRcrea meta-analysis of Liu H et al., Science 2025 (PMID 39913582); file EGFR_EAS_2P2M, from the Susztak lab GWAS2M page (figshare file 49574424); z-score meta-analysis; includes Biobank Japan participants` | `Susztak lab GWAS2M download page (figshare file 49574424)` | `Liu H et al., Science 2025 (PMID 39913582): East Asian stratum (N 282,852) of the multi-ancestry eGFRcrea meta-analysis; z-score meta-analysis; includes Biobank Japan participants` |
| `pheweb.jp (PheWeb.jp CAD phenotype page, PHEWEB_BBJ_CAD: 32,512 cases and 146,214 controls; Kanai M et al., medRxiv 2021, doi 10.1101/2021.09.03.21262975; not a hum0197 release; also distributes Biobank Japan GWAS)` | `the PheWeb.jp CAD page` | `pheweb.jp (BioBank Japan Project / Okada lab; PheWeb fork)` |
| `Japanese heart failure GWAS (GCST90668009; relation to Biobank Japan unconfirmed)` | `GWAS Catalog; GCST90668009` | `Japanese heart failure GWAS (GCST90668009); licence CC0 1.0` |
| `GIANT consortium (BMI, Locke 2015)` | `giant-consortium.web.broadinstitute.org` | `GIANT consortium (Locke et al. 2015)` |
| `International Consortium of Blood Pressure with UK Biobank (Evangelou 2018, systolic blood pressure)` | `Local file SBP_Evangelou2018_local, a hand-relabelled curated file in deposit v4 with no writer script and no download log. The project's preserved OpenGWAS metadata record for this file (deposit v4 frozen intermediate_gwasinfo) names IEU OpenGWAS id ieu-b-38; the manuscript's Additional file 1 Table S2 cites GWAS Catalog GCST006624 as the source of the sample size only; the Evangelou 2018 data statement offers the UKB+ICBP summary data on request to the corresponding authors and the ICBP summary data on request to the ICBP steering committee; route used: IEU OpenGWAS, dataset id ieu-b-38, the route the project's frozen gwasinfo record supports (the download itself is not logged)` | `International Consortium of Blood Pressure (Evangelou 2018)` |
| `Global Lipids Genetics Consortium (LDL cholesterol, triglycerides, Willer 2013)` | `csg.sph.umich.edu/willer/public/lipids2013` | `Global Lipids Genetics Consortium (Willer 2013)` |
| `MAGIC consortium (HbA1c, Wheeler 2017)` | `magicinvestigators.org` | `MAGIC consortium (Wheeler 2017)` |
| `CHARGE Inflammation Working Group (C-reactive protein, Ligthart 2018)` | `IEU OpenGWAS ieu-b-35` | `CHARGE Inflammation Working Group (Ligthart 2018)` |
| `UK Biobank (contributing cohort of the blood-pressure and cystatin-C eGFR GWAS; not the UKB-PPP)` | `UK Biobank` | `UK Biobank (contributing cohort of BP and cystatin-C GWAS)` |
| `UK10K LD scores provided by the lhcMR authors (Darrous 2021)` | `lhcMR authors (download link in the lhcMR README)` | `UK10K LD scores provided by the lhcMR authors (Darrous 2021)` |
| `ChEMBL (mechanism of action, maximum clinical phase)` | `chembl.gitbook.io; Open Targets` | `ChEMBL (EMBL-EBI); licence CC BY-SA 3.0 Unported` |
| `DGIdb (drug-gene interactions)` | `dgidb.org` | `Drug-Gene Interaction Database (WashU)` |
| `UniProt (membrane-protein annotation)` | `uniprot.org; rest.uniprot.org` | `UniProt Consortium; licence CC BY 4.0` |
| `GWAS Atlas (phenome-wide lookup)` | `atlas.ctglab.nl` | `GWAS ATLAS (Watanabe et al. 2019)` |
| `IEU OpenGWAS (ieu-a-7 CARDIoGRAMplusC4D, ieu-b-35 CRP; phenome-wide lookup)` | `opengwas.io; ieu-a-7; ieu-b-35` | `IEU OpenGWAS (University of Bristol)` |
| `GWAS Catalog (GCST accessions)` | `ebi.ac.uk/gwas` | `GWAS Catalog (NHGRI-EBI)` |

## Folders

- `code/`: the code bundle: analysis-tail and figure code, at its project-relative path (release_records/CODE_INCLUDE.tsv, ship = 1).
- `colocalization/`: stored colocalization results and the tiered pair-evidence table.
  Files: `phase3_coloc__coloc_results_final.csv`, `table_s_pair_evidence_tiered.csv`.
- `decoy/`: the decoy-rate script, its log and its result table.
  Files: `S1R1_decoy_rates.R`, `S1R1_decoy_rates.log`, `S1R1_decoy_rates.tsv`.
- `east_asian/`: East Asian ancestry-matched results (workstream W7 and the stored East Asian scans).
  Files: `ancestry_branches.tsv`, `ancestry_results.tsv`, `phase2_CAD_scan.tsv`, `phase2_Stroke_scan.tsv`, `phase2_T2D_scan.tsv`, `phase2_egfr_scan.tsv`.
- `external_outcomes/`: external outcome checks, overlap and LDSC noise correlations; extchecks/ holds the pre-stated external-check workstream.
  Files: `S1R1_outcome_checks.tsv`, `S1R1_overlap_ldsc.tsv`, `extchecks/branches.tsv`, `extchecks/results.tsv`.
- `finemapping/`: stored SuSiE grid results and the candidate edge audit.
  Files: `SuppTable_candidate_edge_audit.csv`, `phase3_susie__susie_grid_results.csv`.
- `kidney/`: renal classification and the UMOD quantity table.
  Files: `a3_umod_q.tsv`, `egfr_mirage_classification.csv`.
- `ledger/`: the public extract of the claims ledger (every number the shipped figure code prints; a row is masked, with its reason, when it holds a raw per-variant TPMI value or a ChEMBL-derived field or cites a source file withheld for another reason; a row that cites a file withheld only because it prints a TPMI line is not masked, because the derived TPMI values ship; evidence_status is defined below).
  Files: `claims_ledger_public.tsv`.
- `literature/`: prior-art search results and the search log.
  Files: `pairs81_priorart.tsv`, `priorart_branches.tsv`, `priorart_literature_search_log.tsv`, `priorart_results.tsv`.
- `proxy/`: the partial proxy record (the full cache is on the analysis machine).
  Files: `proxy_record.tsv`.
- `release_records/`: the authors' own release records: the tier rules, the provider terms, the code-bundle list, and the declared leak-scan coincidences without their values.
  Files: `CODE_INCLUDE.tsv`, `LEAKSCAN_COINCIDENCES.tsv`, `PROVIDER_TERMS.tsv`, `TIER_RULES.tsv`.
- `replication/`: deCODE and ARIC replication results.
  Files: `decode_replication__decode_mr_results.csv`, `replication_sensitivity__panel1_2.csv`.
- `screen/`: the MR screen: one row per protein x outcome pair, non-primary estimates, reason codes and the instrument list.
  Files: `cis_pqtl_instruments.csv`, `screen_ledger.csv`, `screen_nonprimary_estimates.csv`, `screen_reason_codes.csv`.
- `sensitivity_uncited/`: a sensitivity table the manuscript does not cite.
  Files: `SuppTable_LDaware_IVW_multiinstrument.csv`.
- `source_data/`: per-figure plotted frames (source_data/figures/<figure>/), the metabolite bridge screen and figure preparation stamps.
  Files: `figures/FigureS1/FigureS1__posthoc_28.csv`, `figures/FigureS2/FigureS2__candidates_16.csv`, `figures/FigureS3/FigureS3__benchmark.csv`, `figures/FigureS3/FigureS3__renal_proteins.csv`, `figures/FigureS4/FigureS4__genes.csv`, `figures/FigureS4/FigureS4__marked_variants.csv`, `figures/FigureS4/FigureS4__tracks.csv`, `figures/FigureS5/FigureS5__per_signal_estimates.csv`, `figures/FigureS5/FigureS5__posterior_colocalization.csv`, `figures/FigureS5/FigureS5__whole_locus_reconciliation.csv`, `figures/FigureS5_p2/FigureS5_p2__credible_set_grid.csv`, `figures/FigureS5_p2/FigureS5_p2__overview_genes.csv`, `figures/FigureS5_p2/FigureS5_p2__overview_tracks.csv`, `figures/FigureS6/FigureS6__binary_outcome_pairs.csv`, `figures/FigureS6/FigureS6__strict_T2D.csv`, `figures/FigureS6_p2/FigureS6_p2__ckbfree_composite.csv`, `figures/Figure_1/Figure_1.json`, `figures/Figure_2/fig2_plotted.provenance.json`, `figures/Figure_2/fig2a_plotted.tsv`, `figures/Figure_2/fig2b_plotted.tsv`, `figures/Figure_2/fig2c_plotted.tsv`, `figures/Figure_3/Figure_3__3a_forest.csv`, `figures/Figure_3/Figure_3__3a_stored_results.csv`, `figures/Figure_3/Figure_3__3b_forest.csv`, `figures/Figure_4/Figure_4__4a_forest.csv`, `figures/Figure_4/Figure_4__4b_credible_set_members.csv`, `figures/Figure_4/Figure_4__4b_credible_set_rows.csv`, `figures/Figure_4/Figure_4__4b_exons.csv`, `figures/Figure_4/Figure_4__4b_genes.csv`, `figures/Figure_4/Figure_4__4b_protein_track.csv`, `figures/Figure_4/Figure_4__4b_rs2246012.csv`, `figures/Figure_4/Figure_4__4b_signal_leads.csv`, `figures/Figure_5/Figure_5__5a_boxes.csv`, `figures/Figure_5/Figure_5__5a_labels.csv`, `figures/Figure_5/Figure_5__5a_segments.csv`, `figures/Figure_5/Figure_5__5b_forest.csv`, `figures/Figure_5/Figure_5__5b_text.csv`, `l3_metabolite/bridge_screen_results.csv`, `prep_stamps/ARG1_6q23_PREP_STAMP.json`, `prep_stamps/DPEP1_16q24_PREP_STAMP.json`.
- `tables/`: Table 1 and Table 2 as data, the column definitions and the Table 2 literature search log.
  Files: `Table1_data_sources.tsv`, `Table2_evidence_cards.columns.tsv`, `Table2_evidence_cards.tsv`, `literature_search_log.tsv`.
- `workstream_A/`: the instrument verification workstream (reconciliation, verdict and the second-reader outputs).
  Files: `VERDICT.md`, `credible_set_status.tsv`, `estimates_at_variants.tsv`, `ld_pairs_table.tsv`, `results.tsv`.
- `workstreams/`: Stage 1 workstream result and branch tables, one folder per source id (see the map below).
  Files: `BI/branches_integrated.tsv`, `INT/corroboration.tsv`, `INT/results.tsv`, `INT/s6d_fitset_readings.tsv`, `W1/branches.tsv`, `W1/results.tsv`, `W2/branches.tsv`, `W2/results.tsv`, `W3/branches.tsv`, `W3/results.tsv`, `W4a/branches.tsv`, `W4a/results.tsv`, `W4b/branches.tsv`, `W4b/results.tsv`, `W5/branches.tsv`, `W5/results.tsv`, `W6a/branches.tsv`, `W6a/results.tsv`, `W6a/w6a_posthoc_nopurity.tsv`, `W6b/branches.tsv`, `W6b/results.tsv`, `W6b/umod_signal_se.tsv`.

Workstream folders, files released: 22 files across the workstream folders: W1 2, W2 2, W3 2, W4a 2, W4b 2, W5 2, W6a 3, W6b 3, INT 3, BI 1 [L670]. Branch table: 148 rows: W1 9, W2 8, W3 11, W4a 23, W4b 17, W5 5, W6a 5, W6b 28, W7 42 [L672].

## Table 3 source ids and release folders

The repository-row map of Table 3 (below) and `branches_integrated.tsv` name rows by source id. The ids are not all folder names:

| source id | release folder |
|---|---|
| W1 | `workstreams/W1/` |
| W2 | `workstreams/W2/` |
| W3 | `workstreams/W3/` |
| W4a | `workstreams/W4a/` |
| W4b | `workstreams/W4b/` |
| W5 | `workstreams/W5/` |
| W6a | `workstreams/W6a/` |
| W6b | `workstreams/W6b/` |
| W7 | `east_asian/` |
| INT | `workstreams/INT/` |
| BI | `workstreams/BI/` |
| extchecks | `external_outcomes/extchecks/` |

## Repository rows behind the claims of Table 3

Table 3 of the manuscript prints the values; this list names the result rows that hold them, as code, by folder and row-identifier prefix (rows of workstreams/BI by workstream and row number).

| claim | repository rows |
|---|---|
| `**81 [[L010]] colocalization-supported pairs**` | `--` |
| `*TGFB1*–CAD and *ARG1*–T2D rest on inadmissible reference-LD fits (flagged at first mention under Rule 3(i))` | `workstreams/W2: S2.susie; workstreams/W4b: S4.s` |
| `*RAD23B*–T2D significance depends on the SE convention` | `workstreams/W1: S1.se` |
| `**Pre-declared external outcome checks**` | `--` |
| `Six checks, each classified on the resource's own p` | `external_outcomes/extchecks: C1.<pair>.p_wald1` |
| `*RAD23B*–T2D in DIAMANTE: classification depends on the p source and on rounding` | `workstreams/W1: S1.rad23b.ext` |
| `DIAMANTE regional colocalizations (*PAM*–T2D supported at the default prior; *RAD23B*–T2D shared variant not resolved)` | `workstreams/W3: S3.RAD23B_T2D.pmafdrop, S3.PAM_T2D.pmafdrop; workstreams/BI: W3 rows 4, 6` |
| `*PCSK9*–CAD regional colocalization: not a clean control` | `workstreams/W3: S3.PCSK9_CAD.betase; workstreams/BI: W3 rows 10, 11` |
| `*TGFB1*–CAD protein route is instrument-dependent` | `workstreams/W2: S2.het, S2.p489.het` |
| `*TGFB1*–CAD is sensitive to the reference LD` | `workstreams/W2: S2.coloc_clause; workstreams/BI: W2 row 1` |
| `**Kidney pairs among multi-signal candidates**` | `--` |
| `Multi-signal candidates on inadmissible fits` | `workstreams/W6b: S6.d; workstreams/INT: s6d_fitset_readings.tsv` |
| `Kidney share among the 16 [[L222]] candidates against the other 144 [[L264]] pool pairs (exploratory; sensitive to grouping, the acceptance rule and the LD guard)` | `workstreams/W6a: S6.primary, S6.b.a_pair` |
| `as above` | `workstreams/W6a: S6.b.locus_outcome` |
| `as above` | `workstreams/W6a: S6.c` |
| `as above` | `workstreams/W6b: S6.d; workstreams/BI: W6b row 17` |
| `as above` | `workstreams/W6b: S6.d.kept; workstreams/BI: W6b row 17` |
| `as above` | `ledger/claims_ledger_public.tsv: L538 to L542` |
| `*UMOD* behaved as a control` | `workstreams/W6b: S6.umod` |
| `**ARG1**` | `--` |
| `Arginine: no detectable difference between the signals` | `workstreams/W4a: S4.arg, S4.ld.boot.arg` |
| `Ornithine did not scale with ARG1 in the same proportion` | `workstreams/W4a: S4.orn, S4.ld.boot.orn` |
| `T2D A−B difference, nominal, seen in both outcome resources` | `workstreams/W4a: S4.grid; workstreams/BI: W4a row 16` |
| `Direct T2D-per-arginine comparison: no divergence claim` | `workstreams/W4a: S4.t3; workstreams/BI: W4a row 18` |
| `**DPEP1**` | `--` |
| `No detectable difference in cysteinylglycine per SD of DPEP1` | `workstreams/W5: S5.Q; workstreams/BI: W5 row 1` |
| `No detectable difference in eGFR per SD of DPEP1` | `workstreams/W5: S5.Q; workstreams/BI: W5 row 2` |

## Tables

Rows of Table 1: 15 [L673]. Rows of Table 2: 8 [L671]. Table 2 ships whole except the ChEMBL-derived drug and development-stage cells (a stage set by a drug record, the agent list, and statements or verdicts that name a drug), which read `withheld_pending_terms`; the Open Targets tractability class and the derived TPMI values of the external check ship. `tables/Table2_evidence_cards.columns.tsv` defines every column. The ChEMBL facts behind the masked cells are quoted, under the licence of ChEMBL, in shipped scripts, as declared under Code licence.

## The public claims ledger

`ledger/claims_ledger_public.tsv` carries the id, section, claim, value, display, source, evidence_status and source_md5_before_release_transform of every ledger row. Internal notes and the verification history are dropped. A masked value is the JSON object `{"masked":true,"reason":...}`, never zero or blank, and its display reads `[masked: <reason>]`. A row whose value is a list of word statements (the Table 2 cards) stays unmasked except for the development-stage clauses that carry ChEMBL vocabulary (a stage word, such as approved or in clinical development, with the agents it names): each such clause is replaced inside its list element by `[masked: CHEMBL]`, the clauses before it (Open Targets-native, for example tractability) stay, and the display is the same list with the same replacement. `source` names the file, line and column a value was read from and the md5 of that file AS SHIPPED in this release, with the line numbers of the shipped file. Where the release rewrote the file (masking, reduction, portability rewrite or line-break replacement), `source_md5_before_release_transform` holds the md5 of the file before the rewrite, separated by `;` when a row cites several files; for a cited file that is not in the release both are the original md5, and the column is empty for a row that cites no md5. `evidence_status` is one of:

- `re-derived independently`: a second, separately written computation (or a blinded verifier) reproduced the value from the raw sources.
- `source-reconciled`: the value was read from, or asserted equal to, a named source at build time; one reader, not independently re-derived.
- `not separately checked`: computed by the analysis pipeline and carried as computed; no separate check is recorded.

Masked rows in this release carry the reason `CHEMBL`. A row that cites a source file withheld only because it prints a TPMI line is not masked, because the derived TPMI values ship. No row holds a raw per-variant TPMI value, so none is masked for `TPMI_RAW`.

## Superseded ledger rows

`ledger/claims_ledger_public.tsv` keeps the predecessor of every row that a later ledger batch replaced. The claim of the successor names its predecessor as `[successor of Lxxx: ...]`; read the successor. Predecessor by successor: L048 by L532, L094 by L502, L208 by L534, L210 by L963, L242 by L501, L256 by L959, L261 by L960, L262 by L961, L329 by L535, L403 by L503, L409 by L962, L415 by L504, L419 by L537, L438 by L505, L439 by L506, L446 by L507, L454 by L500, L676 by L533, L967 by L536. Two kinds of text keep wording from before the East Asian eGFR file was identified: the claim of the predecessor L967 ("identity and sample size unconfirmed") and the shipped files that carry the note that its N is unconfirmed, found by scanning the shipped tree: `code/R1_stage1/S1R1_R2_07_ancestry/w7_main.R`, `east_asian/ancestry_branches.tsv`, `east_asian/ancestry_results.tsv`, `workstreams/BI/branches_integrated.tsv`. All of them predate the identification of that file as `East Asian stratum (N 282,852) of the multi-ancestry eGFRcrea meta-analysis of Liu H et al., Science 2025 (PMID 39913582); file EGFR_EAS_2P2M, from the Susztak lab GWAS2M page (figshare file 49574424); z-score meta-analysis; includes Biobank Japan participants`; the successors L535 and L536 carry the identified file. The release does not rewrite any of these texts.

## Legacy labels

The label `BBJ_CAD_hum0197.v5` (also written `BBJ CAD (hum0197.v5)`) in the shipped result tables and frames is the run-time name of the PheWeb.jp CAD file (the `pheweb.jp` resource above: Kanai et al.; it is not a `hum0197` release). The release keeps the name so that a shipped table can be traced to the run output that carries it. Result tables and frames that carry it, found by scanning the shipped files: `code/lift_inputs/FigureS6/results_census.tsv`, `source_data/figures/Figure_5/Figure_5__5b_forest.csv`, `source_data/figures/FigureS6_p2/FigureS6_p2__ckbfree_composite.csv`, `workstreams/W2/results.tsv`. Code that writes or tests the label: `code/R1_stage1/S1R1_R2_02_tgfb1/w2_07_twobytwo.py`, `code/R2_rewrite/scripts/P31R_build_af1.py`, `code/V4_local/builders/eas__census_common.py`. The `hum0197.v3` labels (Biobank Japan T2D, ischemic stroke, CKD and HF) are correct and are unchanged.

## Column notes: the p values of the renal files

- `kidney/egfr_mirage_classification.csv`, columns `crea_p`, `cys_p` and `bun_p`: the two-sided normal p of a single-instrument Wald ratio, `2 * pnorm(-abs(b / se))`, recomputed from the effect and standard error that the GWAS prints to four decimal places (script `code/paper3_ckd_rerun_v3/scripts_sent/01_renal_triad_mr.R`). `class` is concordant, mirage or other, from the cystatin-C result and the sign (the script's header defines the three).
- `source_data/figures/FigureS3/FigureS3__renal_proteins.csv`, columns `cys_p` and `bun_p`: the CKDGen GWAS's own P value for the instrument SNP; for a single SNP it is the first-order Wald p. They set the filled and hollow markers of FigureS3 panels b and c (filled: nominally significant). `egfrcrea_p` is the shipped single-instrument eGFR Wald p. The two definitions agree on every call of nominal significance and differ in magnitude (the recomputation from four-decimal values is coarser); Supplementary Section S-R12 documents both.

## Provenance columns in the plotted frames

Some frame columns are not plotted. They tie each plotted value to its ledger row or printed text:

- `source_data/figures/Figure_4/Figure_4__4a_forest.csv`, column `ledger_p`: ledger row id of the p value.
- `source_data/figures/Figure_4/Figure_4__4a_forest.csv`, column `printed_p`: the p text printed on the figure.
- `source_data/figures/Figure_4/Figure_4__4a_forest.csv`, column `ref_ld_sensitive`: flag: drawn dotted because the estimate depends on the reference LD panel.
- `source_data/figures/Figure_3/Figure_3__3b_forest.csv`, column `quantity_ids`: ids of the quantities behind the row.
- `source_data/figures/Figure_3/Figure_3__3b_forest.csv`, column `ledger_ids`: ledger row ids behind the row.
- `source_data/figures/Figure_3/Figure_3__3b_forest.csv`, column `source_file`: file the row was read from.
- `source_data/figures/Figure_3/Figure_3__3b_forest.csv`, column `printed_text`: the text printed on the figure.
- `source_data/figures/Figure_3/Figure_3__3a_stored_results.csv`, column `ledger_id`: ledger row id of the stored value.
- `source_data/figures/Figure_3/Figure_3__3a_stored_results.csv`, column `source`: file the value was read from.
- `source_data/figures/Figure_4/Figure_4__4b_credible_set_rows.csv`, column `printed_text`: the text printed on the figure.

## Withheld and not shipped

Files whose provider terms are unread, and files under source hold, are kept out of this repository. Individual-level dosage files and third-party summary-statistic windows are not released; `REPRODUCE.md` names them per unit.

## Code bundle

`code/` holds the shipped files of `release_records/CODE_INCLUDE.tsv` at their project-relative paths, and the records the release builder wrote into it (the input map, the portability log and the reduced inputs of the level-two figures). Declared but not shipped (pending a decision): rows of that file with ship = 0, of role `deposit_v4_script`, `figure_prep`, `sm_named_script`, `stage1_workstream_script`. Declared ship = 0 but not pending, because the file ships elsewhere in the release: role `decoy_script` (already shipped byte for byte as decoy/S1R1_decoy_rates.R; the core must not copy it again). The bundle files:

- `code/R2_rewrite/figures/scripts/P31R_Fig1.R`
- `code/R2_rewrite/figures/scripts/P31R_Fig2.R`
- `code/R2_rewrite/figures/scripts/P31R_Fig3.R`
- `code/R2_rewrite/figures/scripts/P31R_Fig4.R`
- `code/R2_rewrite/figures/scripts/P31R_Fig5.R`
- `code/R2_rewrite/figures/scripts/P31R_FigS1.R`
- `code/R2_rewrite/figures/scripts/P31R_FigS2.R`
- `code/R2_rewrite/figures/scripts/P31R_FigS3.R`
- `code/R2_rewrite/figures/scripts/P31R_FigS4.R`
- `code/R2_rewrite/figures/scripts/P31R_FigS5.R`
- `code/R2_rewrite/figures/scripts/P31R_FigS6.R`
- `code/R2_rewrite/figures/scripts/P31R_fig_common.R`
- `code/R2_rewrite/figures/af3/P31R_af3_common.R`
- `code/R2_rewrite/figures/af3/P31R_af3_glyphs.R`
- `code/v3.2/scripts/theme_publication_ckm.R`
- `code/v3.2/scripts/ckm_palette.R`
- `code/R2_rewrite/figures/af3/P31R_af3_assemble.py`
- `code/R2_rewrite/figures/af3/build_legend_AF3_R2.py`
- `code/R2_rewrite/figures/scripts/P31R_figgate.py`
- `code/R2_rewrite/scripts/P31R_gate.py`
- `code/R2_rewrite/scripts/P31R_tablegate.py`
- `code/R2_rewrite/scripts/P31R_af1_devstage.py`
- `code/R2_rewrite/scripts/P31R_build_af1.py`
- `code/R2_rewrite/figures/render_all_units.py`
- `code/R2_rewrite/release/P31R_build_release.py`
- `code/R2_rewrite/release/rel_screen.py`
- `code/R2_rewrite/release/rel_tables.py`
- `code/R2_rewrite/release/rel_docs.py`
- `code/R2_rewrite/release/rel_licence.py`
- `code/R2_rewrite/release/P31R_releasegate.py`
- `code/paper3_deposit_request_v4/scripts/R/00_config.R`
- `code/paper3_deposit_request_v4/scripts/R/_phase1_set.R`
- `code/paper3_deposit_request_v4/scripts/R/coloc.R`
- `code/paper3_deposit_request_v4/scripts/scripts/_shared/localize/_local_outcome_extract.R`
- `code/paper3_deposit_request_v4/scripts/scripts/paper3.1_ckm/cascade/01_mr_screen.R`
- `code/R1_stage1/S1R1_R2_01_rad23b/w1_common.R`
- `code/R1_stage1/S1R1_R2_01_rad23b/w1_compare_seen.R`
- `code/R1_stage1/S1R1_R2_01_rad23b/w1_extract.R`
- `code/R1_stage1/S1R1_R2_01_rad23b/w1_main.R`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_00_prestart.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_01_extract.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_02_locus.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_02b_join_audit.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_02c_fitset7.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_03_cond.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_04_susie_coloc.R`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_04b_fitset7.R`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_05_annot.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_06_coloc_boot.R`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_07_twobytwo.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_08_assemble.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_09_compare_seen.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_10_step2_gate.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_11_gate_mutations.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_common.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_rerun_diff.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_step2_gate_lib.py`
- `code/R1_stage1/S1R1_R2_02_tgfb1/w2_step2_make_expected.py`
- `code/R1_stage1/S1R1_R2_03_extcoloc/w3_00_extract.py`
- `code/R1_stage1/S1R1_R2_03_extcoloc/w3_01_coloc.R`
- `code/R1_stage1/S1R1_R2_03_extcoloc/w3_02_selfcheck.py`
- `code/R1_stage1/S1R1_R2_03_extcoloc/w3_03_compare_seen.py`
- `code/R1_stage1/S1R1_R2_03_extcoloc/w3_snapshot.py`
- `code/R1_stage1/S1R1_R2_04a_arg1grid/w4a_compare_seen.R`
- `code/R1_stage1/S1R1_R2_04a_arg1grid/w4a_final_check.R`
- `code/R1_stage1/S1R1_R2_04a_arg1grid/w4a_funcs.R`
- `code/R1_stage1/S1R1_R2_04a_arg1grid/w4a_main.R`
- `code/R1_stage1/S1R1_R2_04a_arg1grid/w4a_staging_beside.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_00_setup.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_01_data.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_02_susie.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_02b_fitset7.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_03_omnibus.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_04_credsets.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_05_rules.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_06_compare_seen.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_06b_step2_gate.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_06c_gate_mutations.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_07_finish.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_common.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_run_all.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_step2_gate_lib.R`
- `code/R1_stage1/S1R1_R2_04b_arg1sets/w4b_step2_make_expected.R`
- `code/R1_stage1/S1R1_R2_06a_renal/w6a_compare.R`
- `code/R1_stage1/S1R1_R2_06a_renal/w6a_compute.R`
- `code/R1_stage1/S1R1_R2_06a_renal/w6a_mutation_tests.R`
- `code/R1_stage1/S1R1_R2_06a_renal/w6a_run_all.R`
- `code/R1_stage1/S1R1_R2_07_ancestry/w7_main.R`
- `code/R1_stage1/S1R1_R2_07_ancestry/w7_mutations.R`
- `code/R1_stage1/S1R1_R2_10_integration/w9_common.R`
- `code/R1_stage1/S1R1_R2_10_integration/w9_compare_seen.R`
- `code/R1_stage1/S1R1_R2_10_integration/w9_main.R`
- `code/R1_stage1/S1R1_R2_10_integration/w9_mutation.R`
- `code/R1_stage1/S1R1_R2_10_integration/w9_readme.R`
- `code/R1_stage1/S1R1_06_claims_ledger/S1R1_build_ledger.py`
- `code/R1_stage1/S1R1_06_claims_ledger/S1R1_id_alloc.py`
- `code/R1_stage1/S1R1_06_claims_ledger/R2_derived_rows.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_af1.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_r2_release.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_r2_stage1.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_r2_stage5.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_supp.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_supp2.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_table1.py`
- `code/R1_stage1/S1R1_06_claims_ledger/ext_table2.py`
- `code/V4_local/builders/l3fm__conditional.py`
- `code/V4_local/builders/l3fm__common.py`
- `code/V4_local/builders/l3fm__panel2.py`
- `code/V4_local/builders/l3fm__prep.py`
- `code/V4_local/builders/l3__build_bridge_screen_eur.py`
- `code/R1_stage1/S1R1_R2_00_resources/gls_q.R`
- `code/R1_stage1/S1R1_R2_06a_renal/exact_cond_or.R`
- `code/R1_stage1/S1R1_R2_05_dpep1/S5_00_run_all.R`
- `code/R1_stage1/S1R1_R2_05_dpep1/S5_01_prepare.R`
- `code/R1_stage1/S1R1_R2_05_dpep1/S5_02_analyses.R`
- `code/R1_stage1/S1R1_R2_05_dpep1/S5_03_bootstrap.R`
- `code/R1_stage1/S1R1_R2_05_dpep1/S5_04_assemble.R`
- `code/R1_stage1/S1R1_R2_05_dpep1/S5_05_selfcheck.R`
- `code/R1_stage1/S1R1_R2_05_dpep1/S5_common.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/00b_log_reads.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/00_common.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/01_setup.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/02_candidates.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/03_windows.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/04_harmonise.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/05b_fitset7.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/05_susie.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/07_fisher_d.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/08_umod.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/09_collect.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/10b_normalise_inputs.py`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/10_branches.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/11_compare_seen.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/12_end.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/13_mutation_tests.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/14_join_yields.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/run_all.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/step2_diff.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/step2_gate_lib.R`
- `code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/step2_make_expected.R`
- `code/R1_stage1/S1R1_R2_11_extchecks/C1_extchecks.R`
- `code/R1_stage1/S1R1_R2_11_extchecks/C1_compare_seen.R`
- `code/R1_stage1/S1R1_R2_11_extchecks/C1_compare_mutation.R`
- `code/R1_stage1/S1R1_05b_outcome_checks/S1R1_00_instruments.R`
- `code/R1_stage1/S1R1_05b_outcome_checks/S1R1_01_extract.py`
- `code/R1_stage1/S1R1_05b_outcome_checks/S1R1_02_outcome_checks.R`
- `code/R1_stage1/S1R1_R2_12_chen_instrument/scripts/a01_chen_extract.py`
- `code/R1_stage1/S1R1_R2_12_chen_instrument/scripts/b01_ld.R`
- `code/R1_stage1/S1R1_R2_12_chen_instrument/scripts/b02_tables.R`
- `code/R1_stage1/S1R1_R2_12_chen_instrument/scripts/b03_snapshot_check.R`
- `code/R2_rewrite/figures/registry/APPROVED_SOURCE_MANIFEST.tsv`
- `code/R2_rewrite/figures/af3/P31R_ACTIVE_FIGURES.json`
- `code/R2_rewrite/release/rel_ledger.py`
- `code/R2_rewrite/release/rel_portability.py`
- `code/V4_local/builders/_vendor/l3__mr_coloc.py`
- `code/V4_local/builders/_vendor/l3__tabix_remote.py`
- `code/v3.2/prose/scripts/v4_facts.py`
- `code/v3.2/prose/scripts/l2l4_facts.py`
- `code/v3.2/prose/scripts/sidecar_claims.py`
- `code/V4_local/builders/eas__census_common.py`
- `code/v3.2/prose/scripts/addfiles.py`
- `code/v3.2/prose/scripts/renumber_refs.py`
- `code/v3.2/prose/scripts/apply_supplement.py`
- `code/v3.2/prose/scripts/docx_edit.py`
- `code/V4_local/builders/_vendor/eas__census_mr_coloc.py`
- `code/V4_local/builders/_vendor/eas__census_tabix_remote.py`
- `code/R1_stage1/S1R1_R2_06a_renal/fix_F1_diffgate.R`
- `code/R2_rewrite/release/rel_tiers.py`
- `code/paper3_ckd_rerun_v3/scripts_sent/01_renal_triad_mr.R`
- `code/paper3_ckd_rerun_v3/scripts_sent/02_renal_coloc_concordance.R`
- `code/R2_rewrite/release/rel_records.py`
- `code/paper3_deposit_request_v5/scripts/R/_coloc_inputs_key.R`
- `code/paper3_deposit_request_v5/scripts/R/_ld_namespace.R`
- `code/paper3_deposit_request_v5/scripts/R/_outcome_prevalence.R`
- `code/paper3_deposit_request_v5/scripts/R/_prior_bands.R`
- `code/paper3_deposit_request_v5/scripts/R/strand_align.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/aric/01_aric_replication_mr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/decode/01_decode_replication_mr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/01_acquire_and_localize_outcomes.py`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/02_eas_outcome_coverage.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/03_eas_outcome_wald.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/04_eas_pqtl_transferability.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/05_assign_eas_state.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/09_phase1_harmonise_hg38.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/49_susie_grid_corrected.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/50_phase2_sensitivity_offline.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/50c_sens5e8_egfr_offline.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/53_correlated_ivw_recompute.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/61_lead_snp_consequences.py`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/_guard.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/07_susie.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/08_susie_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/09_steiger.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/10_prior_sweep.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/11_prior_sweep_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/12_decode.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/12b_aric.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/16_assert_annotation.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase2/11_phase2s_sensitivity.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase2/15_phase2_sens_5e8.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase2/16_sens_5e8_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/10_phase3s_prior_sweep.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/15_susie_coloc.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/16_susie_grid.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/25_prior_sweep_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/26_susie_coloc_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/28_candidate_edge_enumerate.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/29_candidate_susie.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/29a_palindrome_registry.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/30_candidate_edge_audit.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/31_candidate_susie_sensitivity.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/32_dashed_finalize.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/33_sensitivity_provenance.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/34_coloc_robustness_table.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/39_palindromic_freq_gate.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/_candidate_null_classify.R`
- `code/R2_rewrite/review_R2/S6_classrule_recon.py`
- `code/R1_stage1/S1R1_R2_06a_renal/w6a_posthoc_nopurity.R`
- `code/R1_stage1/S1R1_R2_08_proxy_record/w8_proxy_record.py`
- `code/R1_stage1/S1R1_R2_08_proxy_record/w8_raw_scan.py`
- `code/R1_stage1/S1R1_R2_08_proxy_record/w8_export_ld.R`

The release builder also wrote these files into `code/` (the input map, the portability log and the reduced inputs of the level-two figures):

- `code/P31R_INPUT_MAP.tsv`
- `code/P31R_PORTABILITY_LOG.tsv`
- `code/lift_inputs/FigureS1/P31R_posthoc_results.tsv`
- `code/lift_inputs/FigureS6/SuppTable_EAS_axis.csv`
- `code/lift_inputs/FigureS6/results_census.tsv`
- `code/lift_inputs/FigureS6/results_ckbfree.tsv`
- `code/lift_inputs/Figure_2/CLAIM_OUTLINE_novelty_ledger.md`

## Provenance labels of the engine and deposit scripts

A script that came with a deposit package from the analysis machine (the v4 and v5 deliveries) carries one of three labels, in the `provenance_label` column of `RELEASE_MANIFEST.tsv` and of `release_records/CODE_INCLUDE.tsv`. A label describes the ORIGINAL file as the deposit holds it. It does not say that the copy in `code/` is that file byte for byte. Where the portability transform changed the copy, `source_md5` in the manifest is the original's md5, `md5` is the shipped copy's, and `transform` lists the declared changes; the copy is derived from the matched original and is not byte-identical to it. `evidence_md5` and `provenance_evidence` in `CODE_INCLUDE.tsv` give the recorded hash (or prefix) and where it was recorded, and the build compares them with the original.

| label | what the evidence is | what it does not establish |
|---|---|---|
| `full run-time hash match` | a run recorded the full md5 of the file and it equals the delivered original | proof that all helpers, settings or inputs are identical |
| `partial launch-record match` | only a logged launch-time md5 prefix matches the delivered original | full-file run-version verification |
| `timestamp evidence only` | the file's modification time precedes the outputs of the run, and no run recorded a hash of it (a candidate historical version) | a demonstrated match to the executed version |

A script that evaluates sections of other scripts (the candidate SuSiE script evaluates the function blocks of the core and eGFR SuSiE engines; the sensitivity script patches their text in memory) ships with those parents. A match for one file does not verify the files it evaluates, or the settings and inputs of the run.

Files labelled `full run-time hash match`:

- `code/paper3_deposit_request_v4/scripts/R/00_config.R`
- `code/paper3_deposit_request_v4/scripts/R/coloc.R`
- `code/paper3_deposit_request_v5/scripts/R/strand_align.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/15_susie_coloc.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/16_susie_grid.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/26_susie_coloc_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/28_candidate_edge_enumerate.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/29_candidate_susie.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/30_candidate_edge_audit.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/31_candidate_susie_sensitivity.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/32_dashed_finalize.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/33_sensitivity_provenance.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/_candidate_null_classify.R`

Files labelled `partial launch-record match`:

- `code/paper3_deposit_request_v5/scripts/scripts/_shared/decode/01_decode_replication_mr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/61_lead_snp_consequences.py`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/34_coloc_robustness_table.R`

Files labelled `timestamp evidence only`:

- `code/paper3_deposit_request_v4/scripts/R/_phase1_set.R`
- `code/paper3_deposit_request_v4/scripts/scripts/_shared/localize/_local_outcome_extract.R`
- `code/paper3_deposit_request_v4/scripts/scripts/paper3.1_ckm/cascade/01_mr_screen.R`
- `code/paper3_ckd_rerun_v3/scripts_sent/01_renal_triad_mr.R`
- `code/paper3_ckd_rerun_v3/scripts_sent/02_renal_coloc_concordance.R`
- `code/paper3_deposit_request_v5/scripts/R/_coloc_inputs_key.R`
- `code/paper3_deposit_request_v5/scripts/R/_ld_namespace.R`
- `code/paper3_deposit_request_v5/scripts/R/_outcome_prevalence.R`
- `code/paper3_deposit_request_v5/scripts/R/_prior_bands.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/aric/01_aric_replication_mr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/01_acquire_and_localize_outcomes.py`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/02_eas_outcome_coverage.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/03_eas_outcome_wald.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/04_eas_pqtl_transferability.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/eas_axis/05_assign_eas_state.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/09_phase1_harmonise_hg38.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/49_susie_grid_corrected.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/50_phase2_sensitivity_offline.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/50c_sens5e8_egfr_offline.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/53_correlated_ivw_recompute.R`
- `code/paper3_deposit_request_v5/scripts/scripts/_shared/localize/_guard.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/07_susie.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/08_susie_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/09_steiger.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/10_prior_sweep.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/11_prior_sweep_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/12_decode.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/12b_aric.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/cascade/16_assert_annotation.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase2/11_phase2s_sensitivity.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase2/15_phase2_sens_5e8.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase2/16_sens_5e8_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/10_phase3s_prior_sweep.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/25_prior_sweep_egfr.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/29a_palindrome_registry.R`
- `code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/39_palindromic_freq_gate.R`

## Code licence

The MIT licence (`LICENSE`) applies to all shipped code, which the authors wrote, except a file flagged `third_party:<what>`. The `code_rights` column of `release_records/CODE_INCLUDE.tsv` says `authors` for the engine and deposit scripts that carry a provenance label (an empty `code_rights` does not mark a file as third-party), on the basis of a scan of each file for third-party copyright or licence notices (none was found). A file that contains third-party code is flagged `third_party:<what>`, keeps its own licence and is not relabelled. The R packages these scripts call and external programs such as PLINK are not shipped. A code licence does not cover provider values quoted inside scripts; those follow the cases above and `release_records/LEAKSCAN_COINCIDENCES.tsv`. Some shipped scripts quote ChEMBL facts as string literals: agent names, development-stage statements, approval status and ChEMBL identifiers, chiefly those behind the masked cells of `tables/Table2_evidence_cards.tsv`. Those facts are reproduced under the Creative Commons Attribution-ShareAlike licence of ChEMBL, with attribution to ChEMBL (EMBL-EBI) and the licence version as recorded in `release_records/PROVIDER_TERMS.tsv`; the MIT grant does not cover them, and a work that redistributes them inherits the share-alike term.
