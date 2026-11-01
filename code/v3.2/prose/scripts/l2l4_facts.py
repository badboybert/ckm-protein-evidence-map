#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""l2l4_facts.py -- the L2-L4 update of the prose against V4 (native East Asian census, European bridge screen and its
fine-mapping classes, the UK-Biobank-free metabolite arm): every number and every named pair the new text prints,
DERIVED from the registered files, plus the text itself (texts()) with each number formatted from a fact.

One implementation, imported by pending_corrections.py (INSERTS.md blocks), rebase_v4_edits.py (inherited paragraphs),
apply_blocks.py and apply_supplement.py. The numbers-vs-source GATE for this text (verify_prose.py P24) is a separate,
independent implementation and never imports this module.

Specification: V3/v3.2/prose/L2L4_PROSE_SPEC.md (its AUTHOR RESOLUTIONS override the rest), with the corrections the
2026-09-27 integration found re-deriving it (listed in DEVIATIONS below, each with its evidence).

Sources (P21b: every one inside V3/V4 or V3/V4_local):
  S23N  V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv   (the census, 59 rows)
  STG2  V4_local/staging/d2-eas-census/{results_census,results_ckbfree,rsens_results}.tsv, controls/C10, coloc_inputs/
        (md5-pinned in D2_MANIFEST.tsv; carried by the data layer's --control pass 2)
  BSR   V4_local/source_data/l3_metabolite/bridge_screen_results.csv       (the European bridge screen, 84 rows)
  BCS   V4_local/source_data/l3_metabolite/bridge_class_summary.csv        (the registered class table, 21 rows)
  GRID  V4_local/source_data/l3_metabolite/bridge_class_robustness__ens116pc.csv (the 48-cell grid, Ensembl 116)
  NODE  V4_local/source_data/l3_metabolite/node_census_v3.csv + node_coloc/..._side_by_side.tsv
  S26Q  V4_local/source_data/figures/S26_metabolite_coloc_deflation/*.csv
  F08   V4_local/source_data/figures/Figure_08_arg1_node_locus/Figure_08_arg1_node_locus__source.csv
  F09   V4_local/source_data/figures/Figure_09_three_funnel_spine/Figure_09_three_funnel_spine__source.csv
  V4    V4/source_data/supp_tables/{SuppTable_EAS_axis,SuppTable_coloc_robustness_core,SuppTable_locus_attribution,
        SuppTable_network_edges}.csv; V4/source_data/tables_from_rds/phase6_integration___mr_intermediate_all_fdrp.csv
  STG1  V4_local/staging/d1-bridge-eur/CONTROLS.tsv (C5, C6, C9: the arm-ancestry assertion)
  STGF  V4_local/staging/l3-finemap/bridge_classes__rules_v1.csv (gene context of the SMPD3 and GSTA3 leads, glist)
  LB    V4_local/builders/l3fm__common.py (REF_N), l3fm__robustness.py (the declared grid), eas__census_common.py
        (outcome registry: case/control counts for Table 1)
  FIX round (2026-09-27, the four verifier reports of wf_cf05eb46-f00), additional sources:
  CKBREG V4_local/raw/d2-eas-census/ckb_regions/TGFB1__*.tsv (the CKB region file, REPORTED p; md5 in D2_MANIFEST.tsv)
  EAS3  paper 3.eas/analysis/results/phase2_T2D_scan.tsv + phase2_coloc_ledger.tsv (the East Asian T2D scan FGFR4 was
        selected from; the pilot's FGFR4 row is asserted identical to the scan's) and paper 3.eas/analysis/plink2.log
        (plink2 --ld rs451643 rs351855 on 1000 Genomes East Asians) -- read only, as P23b already reads the pilot outputs
  GLIST V4_local/source_data/l3_metabolite/bridge_class_robustness.csv (the glist-hg19 grid) with the pinned inputs of
        the class summary (CLASS_SUMMARY_STAMP.json): staging/l3-finemap/DPEP1_16q24/harmonised.tsv.gz (GRCh37 -> GRCh38)
        and raw/l3-finemap/ensembl116/ensembl116_genes.tsv (the Ensembl 116 gene table)
  CSB   V4_local/builders/l3fm__class_summary.py (the tier rule's DECLARED wording, printed verbatim in the Methods)
  FIX round 3 (2026-09-28): paper 3.eas/INTERIM_REPORT_v1.md + Paper3.EAS_Manuscript_Interim_v1.md (the unsaved HF
        screen's one candidate, CDKN1A; the eGFR candidates carried for transparency), phase2_egfr_scan.tsv (mr_q,
        in_EUR_eGFR); the DPEP1 positions from the pinned harmonised table and the Ensembl 116 gene table (above)
  FIX round 4 (2026-09-28): DANN V4_local/source_data/l3_metabolite/dpep1_locus_annotation.csv (the DPEP1-locus annotation
        table, Ensembl release 116: GRCh38 positions, VEP consequences on the canonical transcripts, distances to DPEP1 and
        CHMP1A under full gene spans / protein-coding transcripts / canonical transcripts, the transcripts containing each
        variant; built by V3/V4_local/builders/l3fm__dpep1_locus.py, registered in V3/V4_local/MANIFEST.tsv); the V4 MR
        screen's FGFR4 x T2D row (the printed European estimate's outcome GWAS, FinnGen R10) and the paper 3.eas records
        (the scan's European type 2 diabetes comparator, DIAMANTE)
  FIX round 5 (2026-09-29): PSE / PSS / PSC V4_local/source_data/l3_metabolite/bridge_persignal_{evidence,sensitivity,cells}.csv
        (per-signal colocalization evidence and the POST HOC per-signal attribution sensitivity in both directions, any
        colocalizing signal / every one; built by V3/V4_local/builders/l3fm__persignal.py, registered in
        V3/V4_local/MANIFEST.tsv; the registered class table is unchanged and asserted so there)

  python -B l2l4_facts.py            # print every fact and its source
  python -B l2l4_facts.py --texts    # print every generated paragraph
"""
from __future__ import print_function

import ast
import csv
import gzip
import io
import math
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PROSE = os.path.dirname(HERE)
V32 = os.path.dirname(PROSE)
V3 = os.path.dirname(V32)
ROOT = os.path.dirname(V3)
V4SD = os.path.join(V3, "V4", "source_data")
SUPP4 = os.path.join(V4SD, "supp_tables")
TFR4 = os.path.join(V4SD, "tables_from_rds")
LOCAL = os.path.join(V3, "V4_local")
LSD = os.path.join(LOCAL, "source_data")
STG = os.path.join(LOCAL, "staging")
LB = os.path.join(LOCAL, "builders")

S23N = os.path.join(LSD, "figures", "S23_four_quadrant_eur_eas", "S23_native_census__source.csv")
RC = os.path.join(STG, "d2-eas-census", "results_census.tsv")
RCF = os.path.join(STG, "d2-eas-census", "results_ckbfree.tsv")
RSENS = os.path.join(STG, "d2-eas-census", "rsens_results.tsv")
C10 = os.path.join(STG, "d2-eas-census", "controls", "C10_W1_W2_per_region.tsv")
CIN = os.path.join(STG, "d2-eas-census", "coloc_inputs")
BSR = os.path.join(LSD, "l3_metabolite", "bridge_screen_results.csv")
BCS = os.path.join(LSD, "l3_metabolite", "bridge_class_summary.csv")
GRID = os.path.join(LSD, "l3_metabolite", "bridge_class_robustness__ens116pc.csv")
NODE = os.path.join(LSD, "l3_metabolite", "node_census_v3.csv")
NODEC = os.path.join(LSD, "l3_metabolite", "node_coloc", "node_coloc_metabolite_results_ukbfree__side_by_side.tsv")
S26Q = os.path.join(LSD, "figures", "S26_metabolite_coloc_deflation", "S26_metabolite_coloc_deflation__quadrant__source.csv")
S26C = os.path.join(LSD, "figures", "S26_metabolite_coloc_deflation", "S26_metabolite_coloc_deflation__colocgate__source.csv")
F08 = os.path.join(LSD, "figures", "Figure_08_arg1_node_locus", "Figure_08_arg1_node_locus__source.csv")
F09 = os.path.join(LSD, "figures", "Figure_09_three_funnel_spine", "Figure_09_three_funnel_spine__source.csv")
EAS_AXIS = os.path.join(SUPP4, "SuppTable_EAS_axis.csv")
CORE = os.path.join(SUPP4, "SuppTable_coloc_robustness_core.csv")
LOCUS = os.path.join(SUPP4, "SuppTable_locus_attribution.csv")
EDGES = os.path.join(SUPP4, "SuppTable_network_edges.csv")
INTER = os.path.join(TFR4, "phase6_integration___mr_intermediate_all_fdrp.csv")
D1CTL = os.path.join(STG, "d1-bridge-eur", "CONTROLS.tsv")
RULES1 = os.path.join(STG, "l3-finemap", "bridge_classes__rules_v1.csv")
FM_COMMON = os.path.join(LB, "l3fm__common.py")
FM_SUSIE = os.path.join(LB, "l3fm__susie.R")
FM_ROB = os.path.join(LB, "l3fm__robustness.py")
STAGE1 = os.path.join(PROSE, "_stage1_citations.docx")
CKBREG = os.path.join(LOCAL, "raw", "d2-eas-census", "ckb_regions")
EAS3 = os.path.join(ROOT, "paper 3.eas", "analysis")
SCAN = os.path.join(EAS3, "results", "phase2_T2D_scan.tsv")
LEDGER = os.path.join(EAS3, "results", "phase2_coloc_ledger.tsv")
PLINKLOG = os.path.join(EAS3, "plink2.log")
GLIST = os.path.join(LSD, "l3_metabolite", "bridge_class_robustness.csv")
HARM_DPEP1 = os.path.join(STG, "l3-finemap", "DPEP1_16q24", "harmonised.tsv.gz")
ENSGENES = os.path.join(LOCAL, "raw", "l3-finemap", "ensembl116", "ensembl116_genes.tsv")
CSB = os.path.join(LB, "l3fm__class_summary.py")
DANN = os.path.join(LSD, "l3_metabolite", "dpep1_locus_annotation.csv")
PSE = os.path.join(LSD, "l3_metabolite", "bridge_persignal_evidence.csv")        # FIX round 5
PSS = os.path.join(LSD, "l3_metabolite", "bridge_persignal_sensitivity.csv")
PSC = os.path.join(LSD, "l3_metabolite", "bridge_persignal_cells.csv")
MRS4 = os.path.join(TFR4, "phase2_mr_screen__mr_all_results.csv")

GW_Z = 5.4513            # |z| of p = 5e-8, two-sided: the fine-mapping class rules' genome-wide line (A1). NOT the
#                          census's: the census defines genome-wide significance on the REPORTED CKB p (eas__census_common:
#                          sig = pvalue < GWS), and the TGFB1 count below uses that rule (verify:numbers, FIX round).

# Deviations from L2L4_PROSE_SPEC.md, each found re-deriving the spec against the registered files (2026-09-27):
DEVIATIONS = [
    "42a 'Every instrument had F >= 29': the weakest instrument (FURIN) has F = 28.99 (results_census.tsv); S23 rounds it "
    "to 29.0. Printed as 'the weakest instrument had F = 29.0'.",
    "101a DPEP1: 'ratio -0.14 to -0.15 at both DPEP1 signals examined' are the SAME shared signal read on two SNP sets "
    "(signal B on the palindrome-dropped set, signal A on the frequency-resolved set); E351Q proportionality and the "
    "rs58290281-in-CHMP1A clause are not in a registered file (AUTHOR RESOLUTION 4). Replaced by what the registered "
    "grid holds: shared signal and proportionality in 48/48 variants, attribution 16/48 (glist) and 0/48 (Ensembl 116 "
    "gene spans), rival CHMP1A.",
    "101b 'ATXN2 ... PP.H4 <= 0.08': the three values are 0.050/0.071/0.082 (0.082 > 0.08). Printed as a rounded range.",
    "101b/125 'each points to a protein that neither platform measures': false for BRAP -- SH2B3 is a measured neighbour "
    "(registered grid rivals_full) and does not share the signal; ATXN2 is measured. Rewritten per locus.",
    "101b GSTA3 'both of its leads lie in GSTA2': the shared-variant lead is NEAREST GSTA2 (1 kb), not inside it; the "
    "r = 0.9985 GSTA1/GSTA3 correlation is only in ADJUDICATION.md. Printed as 'its lead variants map to GSTA2'.",
    "100 'a broad plateau of about 150 kb': variants within one log10 unit of the T2D peak span 157 kb (the >= 9 span "
    "is 148 kb). Printed with its definition.",
    "53c 'the release used for the study's variant annotation': dropped (CORRECTIONS: VEP 116 is cited for variant "
    "consequences only). AUTHOR RESOLUTION 3 applied: the annotation was fixed after the first results.",
    "124a 'the decisive variant lies inside a CHMP1A transcript': not registered (RESOLUTION 4); replaced.",
    "124a DPEP1 enzyme: Adachi 1989 shows the human enzyme hydrolysing cysteinylglycine ADDUCTS (cystinyl-bis-glycine), "
    "not cysteinylglycine itself; the citation is placed on that clause.",
    "74: TGFB1 sign robustness to the tie rule (RESOLUTION 1): every genome-wide-significant candidate instrument gives a "
    "negative Wald estimate. FIX round (verify:numbers): counted on the census's OWN rule, the REPORTED CKB p < 5e-8 in "
    "the lead's +/-500 kb window with MAF >= 0.01 and present in the outcome GWAS (= results_census n_domain_shared): "
    "43, all negative (-0.419 to -0.344); the first build printed 46, from |z| recomputed on rounded beta/SE.",
    "42a CITE-NEEDED BBJ CAD: pheweb.jp/pheno/CAD (32,512 cases / 146,214 controls = the census registry) cites Kanai M "
    "et al., medRxiv 2021, doi:10.1101/2021.09.03.21262975 (hum0197.v5; not in PubMed).",
    # ---- FIX round (2026-09-27), from the four verifier reports of wf_cf05eb46-f00; each re-derived before changing
    "FIX 42b/75b FGFR4-T2D: disclosed as the one East-Asian-specific survivor of the 331-protein East Asian T2D scan it was "
    "selected from (paper 3.eas phase2_T2D_scan / phase2_coloc_ledger; the pilot row IS the scan row); BBJ is part of "
    "AGEN, so the CKB-free re-run is not independent replication; its instrument rs451643 tags FGFR4 Gly388Arg "
    "(rs351855; plink2 r2 0.968 on 1000G EAS) -- the epitope/altered-function caveat sits where the direction is read.",
    "FIX 53c/124a DPEP1: 'established substrate' -> 'an enzyme with documented activity on the metabolite', the class "
    "summary's declared wording (amendment 2 of its declaration); Kozak & Tate 1982 (reduced Cys-Gly a potent inhibitor of "
    "the rat renal microvillus dipeptidase) is cited where the chemistry is read.",
    "FIX 101a DPEP1: the non-attribution explained from the registered grids -- in the 24 palindrome-drop variants the "
    "shared signal is rs62068372, inside VPS9D1 (Ensembl 116); in the other 24 the lead falls inside CHMP1A under Ensembl "
    "116 (0/24) and is attributed under glist-hg19 in 16, the other 8 blocked by the CHMP1A protein. 'does not hold' -> "
    "'unresolved'.",
    "FIX 74: TGFB1-CAD's European colocalization is prior-fragile (0.753 at p12 = 1e-6), unlike CELSR2-CAD and FGF5-eGFR.",
    "FIX frame: abstract Methods/Conclusions, Background objective, Study Design, Discussion opening, [117] and the main "
    "Conclusions carry the census and the metabolite layer; the ARG1 triangle is described as European (S26 arginine-T2D "
    "EUR-only, heterogeneity p 0.006); the urea-cycle convergence closes the ARG1 discussion.",
    # ---- FIX round 3 (2026-09-28), from the three verifier reports of wf_876cfb67-e41; each re-derived before changing
    "FIX3 75b/112b FGFR4-T2D: Results and Discussion characterise it the same way -- the East-Asian-specific PATTERN (null "
    "European estimate, heterogeneous) is what the scan selected it for; the CKB-free check is stated separately, not as a "
    "selection criterion. Its own Discussion paragraph; the instrument caveat said plainly (both variants at the smallest p "
    "the CKB resource reports; r2 0.97), without the 'rsID-bearing ... reporting floor' pipeline wording.",
    "FIX3 42b Methods: the heart-failure screen's one candidate CDKN1A was shared (paper 3.eas INTERIM_REPORT_v1.md and "
    "Paper3.EAS_Manuscript_Interim_v1.md), so 'neither contributed a candidate' was false; the five eGFR candidates are the "
    "eGFR scan's FDR-significant hits absent from the European eGFR hit table, carried for transparency (records); the "
    "409-word paragraph split (census / pilot).",
    "FIX3 101a DPEP1: the palindrome-kept half's lead (rs58290281, GRCh38 16:89,641,375) is 2,919 bp downstream of DPEP1 "
    "and 559 bp inside CHMP1A's Ensembl 116 span; the drop half's rs62068372 is in VPS9D1, 80,243 bp from DPEP1 with four "
    "protein-coding genes between (not 'the neighbouring VPS9D1'); dropping palindromes removes rs1126464 itself. The "
    "direction is stated once per paragraph (Results F3, Discussion I2); the Discussion gene caveat covers both halves.",
    "FIX3 abstract/Conclusions: 'robust to fine-mapping' restored (framing decision 4) and 'lower T2D risk' made explicit "
    "(350 by Word); the main Conclusions report 7 colocalized and 3 held; numerals throughout one comparison.",
    "FIX3 124b: Piatti 2001 (PMID 11347747) is an oral L-arginine insulin-sensitivity trial, not a vessel-wall study -> "
    "'Clinical studies'; 104 (G2): the missing comma after 'nominations'.",
    # ---- FIX round 4 (2026-09-28), from the two verifier reports of wf_526188f7-dad; each re-derived before changing
    "FIX4 101a DPEP1 (undersell): the registered DPEP1-locus annotation table (Ensembl 116) replaces '559 bp inside the "
    "Ensembl release 116 span of CHMP1A': the shared-variant lead rs58290281 lies between DPEP1 and CHMP1A, 2,919 and "
    "2,862 bp from their protein-coding transcripts; it is inside CHMP1A's gene span only through ENST00000676275 "
    "(protein_coding_CDS_not_defined, not GENCODE basic); the nearer gene is CHMP1A by 57 bp over protein-coding transcripts "
    "and DPEP1 by 118 bp over the canonical ones; DPEP1's signal A peaks at rs409170 in a DPEP1 intron; rs1126464 is a DPEP1 "
    "missense variant (Glu351Gln, canonical MANE transcript). The brief's 'about 2.9 kb from each gene's coding sequence' "
    "is FALSE: the coding sequences lie 3,153 bp (DPEP1) and 4,539 bp (CHMP1A) away, so the text says protein-coding "
    "transcripts. The brief's 'a margin of tens of bp' is 57 bp one way and 118 bp the other; both are printed.",
    "FIX4 53c Methods: gene bodies are the FULL spans of the Ensembl 116 protein-coding genes, non-coding transcripts "
    "included (the class rules read the gene record's start..end); the over-long sentence split after the annotation.",
    "FIX4 75b FGFR4: 'exactly' dropped -- the scan judged the East-Asian-specific pattern against DIAMANTE (paper 3.eas "
    "records; European p 0.21), while the printed European estimate is FinnGen R10 (V4 MR screen; p 0.27); one clause says "
    "so. 112c (Hb): 'the pattern an East Asian scan had selected it for'.",
    "FIX4 Discussion: H ends 'Transfer across ancestry is therefore selective rather than wholesale.' (the triad was stated "
    "a third time); I2 'For DPEP1 the responsible gene also remains open'; F3 'whose strand the allele frequency resolves'.",
    # ---- FIX round 5 (2026-09-29), from the two verifier reports of wf_fce8d104-bec; each re-derived before changing
    "FIX5 101a DPEP1 (undersell): in the palindrome-kept half cysteinylglycine colocalizes with BOTH DPEP1 signals "
    "(conditional PP.H4 0.99994 / 0.99925; the second's top shared variant is the missense rs1126464) and the class rule "
    "reads only the stronger (margin 0.0007). Registered per-signal tables (l3fm__persignal.py, POST HOC, both directions); "
    "F3a states the dual colocalization, the single-signal reading, the any-signal count (8 of 24, the full-panel "
    "not-SuSiE-as-fitted cells) and the every-signal result (neither DPEP1 nor ARG1: the top variant ARG1's lead signal "
    "shares with arginine lies in MED23 [FIX6: round 5 wrote 'ARG1's lead signal peaks in MED23', which is false]). The "
    "brief's 'that choice alone produces 0/24' overstated the counterfactual (reading both signals gives 8/24, not 24/24: "
    "the missense signal is stable in 151 of 200 LD draws, short of 180, and SuSiE as fitted blocks 4 more). Dropped: the "
    "rs409170-intron sentence (a 0.006-z tie with rs423135) and the 57/118-bp transcript-set sentence (the lead's proxy "
    "rs2139455, r2 0.991, leads in 75 of 200 draws, so a base-pair margin is not a property of the signal).",
    "FIX5 53d (new Methods block B5): the per-signal sensitivity defined as post hoc; the UMOD-eGFR renal anchor, whose two "
    "signals both colocalize, is attributed to UMOD in all 48 variants under the any-signal reading and in none under the "
    "registered rule (Ensembl gene spans).",
    "FIX5 124a I2: 'rather than its neighbours' was inaccurate -- the kept half's competitor is the ADJACENT CHMP1A; the "
    "dropped half's rs62068372 lies in VPS9D1, 80,243 bp from DPEP1, whose gene body starts 68,672 bp past DPEP1 with "
    "four protein-coding genes between (CHMP1A, SPATA33, CDK10, SPATA2L) (the brief's 'VPS9D1 is ~80 kb away' conflated the variant with the gene). Now "
    "'the adjacent CHMP1A is not excluded', with the bound 'at most a third of the rule variants' (16 of 48, the largest "
    "count under any annotation or reading). 'attribution to ARG1 held in all 48 rule variants of the fine-mapping grid' "
    "(the declared grid; the every-signal reading attributes neither enzyme and is reported in F3a).",
    # ---- FIX round 6 (2026-09-29), from the round-5 verify:reviewer report (4 MAJOR framing defects) and two brief pre-audits
    "FIX6 101a F3a (M1): 'whose lead signal likewise peaks outside its gene, in the neighbouring MED23' was FALSE -- ARG1's "
    "lead pQTL rs2781668 lies inside ARG1 (in the ARG1/MED23 overlap, 12,089 bp under Ensembl 116, 10,368 bp under "
    "glist-hg19); only the top variant the lead signal shares with arginine after conditioning (rs2749924 kept half / "
    "rs2608954 dropped half) lies in MED23. The brief's '121/200 LD draws' for that placement is boot_n_lead_as_full "
    "(lead identity), not a placement count (the shared variant lands in ARG1 in 77/200); not printed. The clause is gone "
    "from F3a; B5 says 'the top variant that ARG1's lead signal shares with arginine lies in MED23'.",
    "FIX6 101a F3a (M2, m2, m4): finding first; no per-signal COUNT in Results (the any-signal 8 of 24 printed one "
    "direction); the principled contrast from bridge_persignal_evidence.csv boot_n_qualifies_and_A2_ens116pc (missense "
    "signal 151 of 200 < BOOT_MIN 180; the ARG1 signal the rule reads, rs17788484, 198 -- ARG1's SECOND signal, so the text "
    "says 'the ARG1 signal the class rule reads', never 'ARG1's lead signal'); the verdict bounds DPEP1 at 16 of 48 over "
    "every rule, reading (every-signal included) and annotation. Epitope reconciliation (m4): cysteinylglycine is a "
    "Metabolon HD4 UPLC-MS/MS measurement (Chen et al. 2023, [35], Methods 'Metabolite and metabolite ratio data "
    "processing'; the manuscript's own Methods name no platform), so ITS association at Glu351Gln cannot reflect altered "
    "antibody binding; the Olink DPEP1 association there still could -- the clause says both. Cited to [35] only (already "
    "cited in the Methods), not to the epitope reference, which is first cited in the Discussion (P18 order).",
    "FIX6 53d B5 (M2, M3, M4, m5): the every-signal result stated ONCE here with its calibration: under full gene spans it "
    "attributes no pair at all (the registered rule attributes 3 rows, the any reading 5), under glist-hg19 only DPEP1 "
    "rows; UMOD-eGFR is missed by it AND by the registered rule because both need the more strongly colocalizing signal "
    "(B, rs77924615, cond PP.H4 1.000 vs 0.97; the WEAKER pQTL) to meet the position condition, and it lies in PDILT "
    "(the brief's 'reads the stronger signal' was wrong: stronger colocalization, weaker pQTL). UMOD reported under both "
    "named annotations (any 48 full spans / 0 glist; 'conditionally' added: it fails single-variant colocalization). The "
    "UMOD anchor is not called a positive control of gene attribution (it is registered as the multi-signal pipeline "
    "control). Invariance: ARG1-arginine is the only BRIDGE the registered rule and the any-signal reading attribute in "
    "all 48 variants under both annotations (every-biotype identical); the every-signal exception is named. Provenance "
    "exact: 'specified after its results for the DPEP1, ARG1 and UMOD pairs were known' REPLACES 'added after the grid "
    "results had been seen' (l3fm__persignal.py docstring: written AFTER the outcome was known).",
    "FIX6 124a I2 (m3): the DPEP1 sentence leads with the finding (the missense-led shared signal), then the bound and "
    "CHMP1A, then 'so the responsible gene remains open'.",
]


class L2L4Error(Exception):
    pass


# ---------------------------------------------------------------------------------------------------------- readers
def _rows(path, delim=","):
    if not os.path.isfile(path):
        raise L2L4Error("missing source: %s" % path)
    with io.open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def _T(x):
    return str(x).strip().upper() == "TRUE"


def _need(cond, msg):
    if not cond:
        raise L2L4Error(msg)


# ---------------------------------------------------------------------------------------------------------- format
MINUS = u"−"
SUP = str.maketrans("0123456789-", u"⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
         10: "ten", 11: "eleven"}


def sgn(x, nd):
    """signed with the typographic minus: +0.164 / -0.371 -> U+2212."""
    s = ("%." + str(nd) + "f") % x
    return (MINUS + s[1:]) if s.startswith("-") else ("+" + s)


def neg(x, nd):
    s = ("%." + str(nd) + "f") % x
    return (MINUS + s[1:]) if s.startswith("-") else s


def sci(p, nd=0):
    """4.7 x 10^-9 with superscripts; nd = mantissa decimals."""
    m, e = (("%." + str(nd) + "e") % p).split("e")
    return u"%s × 10%s" % (m, str(int(e)).translate(SUP))


def thou(n):
    return "{:,}".format(int(n))


def pair(p, o):
    return u"%s–%s" % (p, o)


def oxford(xs):
    xs = list(xs)
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


# ------------------------------------------------------------------------------------------------------------- facts
_CACHE = {}


def facts():
    if "F" in _CACHE:
        return _CACHE["F"], _CACHE["SRC"]
    F, SRC = {}, {}

    def put(k, v, src):
        F[k] = v
        SRC[k] = src

    # ---- the network (V4 A11 via v4_facts) ---------------------------------------------------------------------------
    sys.path.insert(0, HERE)
    import v4_facts
    vf, _ = v4_facts.facts()
    put("n_pairs", vf["n_pairs"], "v4_facts (A11)")
    put("n_genes", vf["n_genes"], "v4_facts (A11)")

    # ---- D2: the native census (S23N) ---------------------------------------------------------------------------------
    s23 = _rows(S23N)
    cen = [r for r in s23 if r["set"] == "census_network"]
    pil = [r for r in s23 if r["set"] == "pilot_offnetwork"]
    _need(len(s23) == len(cen) + len(pil), "S23N carries a set other than census_network/pilot_offnetwork")
    key = lambda r: (r["protein"], r["outcome"])
    P = lambda r: pair(r["protein"], r["outcome"])
    put("s23n_rows", len(s23), "S23N")
    put("cen_n", len(cen), "S23N set=census_network")
    put("cen_prot", len({r["protein"] for r in cen}), "S23N census proteins")
    nt = [r for r in cen if r["eas_call"].startswith("NOT_TESTABLE")]
    tst = [r for r in cen if not r["eas_call"].startswith("NOT_TESTABLE")]
    _need([P(r) for r in nt] == [pair("APOE", "NAFLD")] and nt[0]["eas_call"] == "NOT_TESTABLE_NO_EAS_GWAS",
          "the untestable census pairs moved: %s" % [P(r) for r in nt])
    put("cen_test", len(tst), "S23N eas_call")
    _need(all(r["eur_call"] == "H4" and r["eur_support"] == "True" for r in cen),
          "a census row is not European-colocalized (the census is defined on network pairs)")
    # instruments
    lead = [r for r in tst if r["eas_instr_is_lead"] == "True"]
    prox = sorted(P(r) for r in tst if r["eas_instr_is_proxy"] == "True")
    _need(len(lead) + len(prox) == len(tst), "lead + proxy != testable")
    put("cen_lead", len(lead), "S23N eas_instr_is_lead")
    put("cen_proxy", prox, "S23N eas_instr_is_proxy")
    rc = {(r["gene"], r["outcome"]): r for r in _rows(RC, "\t")}
    _need(set(rc) == {key(r) for r in cen}, "results_census.tsv and S23N disagree on the census pairs")
    ties = sorted(pair(*k) for k, r in rc.items() if r["status"] == "ok" and int(r["n_tied_top_z"]) > 1)
    put("cen_ties_pz", ties, "results_census.tsv n_tied_top_z > 1 (tied on p AND |z|)")
    fmin = min((float(r["F"]), g) for (g, o), r in rc.items() if r["status"] == "ok")
    put("cen_Fmin", fmin[0], "results_census.tsv F (full precision)")
    put("cen_Fmin_gene", fmin[1], "results_census.tsv F")
    _need(fmin[0] > 10, "a census instrument is weak (F <= 10)")
    # palindromic instruments
    pal = sorted(k for k, r in rc.items() if r["status"] == "ok" and {r["instr_ea"], r["instr_oa"]} in ({"A", "T"}, {"C", "G"}))
    put("cen_pal_pairs", [pair(*k) for k in pal], "results_census.tsv instr_ea/instr_oa")
    put("cen_pal_variants", len({rc[k]["instr"] for k in pal}), "results_census.tsv instr")
    # window
    wins = {int(r["coloc_hi_hg38"]) - int(r["coloc_lo_hg38"]) for r in rc.values() if r["status"] == "ok"}
    _need(wins == {1000000}, "census colocalization windows are not all +/-500 kb: %s" % wins)
    put("cen_window_kb", 500, "results_census.tsv coloc_hi_hg38 - coloc_lo_hg38 = 1,000,000")
    # outcome resources and joins
    acc = {}
    for r in tst:
        acc.setdefault(r["outcome"], set()).add((r["eas_outcome_accession"], r["eas_outcome_resource"]))
    _need(all(len(v) == 1 for v in acc.values()), "a census outcome uses two resources")
    acc = {o: list(v)[0] for o, v in acc.items()}
    put("cen_acc", {o: a for o, (a, _) in acc.items()}, "S23N eas_outcome_accession")
    m = re.search(r"N ([\d,]+)", acc["eGFR"][1])
    _need(m, "the eGFR resource label carries no N")
    put("egfr_eas_n", int(m.group(1).replace(",", "")), "S23N eas_outcome_resource (eGFR)")
    lift = [r for r in tst if "hg38->hg19" in r["eas_join"]]
    direct = [r for r in tst if "hg38 direct" in r["eas_join"]]
    _need(len(lift) + len(direct) == len(tst), "a census join is neither lifted nor direct")
    put("cen_lifted", len(lift), "S23N eas_join")
    # CKB overlap
    ov = Counter(acc[r["outcome"]][0] for r in tst if r["eas_ckb_overlap"] == "True")
    _need(set(ov) == {"GCST010118", "GCST90104544"}, "the CKB-overlapping resources moved: %s" % dict(ov))
    put("ckb_ov_agen", ov["GCST010118"], "S23N eas_ckb_overlap (AGEN)")
    put("ckb_ov_giga", ov["GCST90104544"], "S23N eas_ckb_overlap (GIGASTROKE)")
    put("ckb_ov_resources", len(ov), "S23N eas_ckb_overlap")
    # European-consistent priors (R coloc 5.2.3)
    rs = [r for r in _rows(RSENS, "\t") if r["kind"] == "census"]
    _need(len(rs) == len(tst), "rsens_results.tsv census rows != testable pairs")
    put("rsens_changed", sum(1 for r in rs if r["call_changed"] != "FALSE"), "rsens_results.tsv call_changed")
    _need(all((r["eas_rsens_support"] == r["eas_support"]) for r in tst), "S23N eas_rsens_support differs from eas_support")
    # W1/W2 channel checks
    c10 = _rows(C10, "\t")
    put("w_regions", len(c10), "controls/C10_W1_W2_per_region.tsv")
    put("w2_pass", sum(1 for r in c10 if r["W2"] == "True"), "C10 W2")
    w1f = sorted((r["gene"], float(r["W1_spearman"])) for r in c10 if r["W1"] != "True")
    put("w1_fail", w1f, "C10 W1")
    _need(len(c10) == F["cen_prot"], "C10 regions != census proteins")
    # the calls
    sup = sorted(P(r) for r in tst if r["eas_support"] == "True")
    rob = sorted(P(r) for r in tst if r["eas_support_p12_1e6"] == "True")
    put("cen_coloc", sup, "S23N eas_support")
    put("cen_robust", rob, "S23N eas_support_p12_1e6")
    _need(set(rob) <= set(sup), "a robust pair is not colocalized at the default prior")
    nonc = [r for r in tst if r["eas_support"] != "True"]
    calls = Counter(r["eas_call"] for r in nonc)
    put("cen_nonc", len(nonc), "S23N")
    put("cen_h1", calls["NOT_SUPPORTED_H1"], "S23N eas_call")
    put("cen_h3_pairs", sorted(P(r) for r in nonc if r["eas_call"] == "NOT_SUPPORTED_H3"), "S23N eas_call")
    h4b = [r for r in nonc if r["eas_call"] == "NOT_SUPPORTED_H4_below_rule"]
    _need(len(h4b) == 1 and P(h4b[0]) == pair("DNER", "eGFR"), "the H4-below-rule pair moved")
    put("dner_h4", _f(h4b[0]["eas_coloc_H4"]), "S23N DNER-eGFR eas_coloc_H4")
    _need(sum(calls.values()) == calls["NOT_SUPPORTED_H1"] + len(F["cen_h3_pairs"]) + 1, "an unexpected non-colocalizing call")
    g = {key(r): r for r in s23}
    # per-pair numbers
    ce = g[("CELSR2", "CAD")]
    put("celsr2", dict(h4=_f(ce["eas_coloc_H4"]), h4_6=_f(ce["eas_coloc_H4_p12_1e6"]), eas=_f(ce["eas_beta_ha"]),
                       eur=_f(ce["eur_beta_ha"]), hetp=_f(ce["het_p"])), "S23N CELSR2-CAD")
    fg5 = g[("FGF5", "eGFR")]
    put("fgf5", dict(h4=_f(fg5["eas_coloc_H4"]), h4_6=_f(fg5["eas_coloc_H4_p12_1e6"])), "S23N FGF5-eGFR")
    _need(("%.3f" % F["fgf5"]["h4"], "%.3f" % F["fgf5"]["h4_6"]) == ("1.000", "1.000"), "FGF5-eGFR is not 1.000 at both priors")
    _need(max(_f(r["eas_coloc_H4"]) for r in tst) == F["fgf5"]["h4"], "FGF5-eGFR is not the census's most decisive colocalization")
    tg = g[("TGFB1", "CAD")]
    tgr = rc[("TGFB1", "CAD")]
    _need({tgr["instr_ea"], tgr["instr_oa"]} not in ({"A", "T"}, {"C", "G"}), "the TGFB1 instrument is palindromic")
    loc = {(r["gene"], r["outcome_name"]): r for r in _rows(LOCUS)}
    eur_pos = int(loc[("TGFB1", "CAD")]["instr_pos"])
    eas_pos = int(tg["eas_instr_pos_hg38"].split(":")[1])
    _need(loc[("TGFB1", "CAD")]["instr_chrom"] == tg["eas_instr_pos_hg38"].split(":")[0], "TGFB1 instruments on two chromosomes")
    put("tgfb1", dict(h4=_f(tg["eas_coloc_H4"]), h4_6=_f(tg["eas_coloc_H4_p12_1e6"]), eur=_f(tg["eur_beta_ha"]),
                      eas=_f(tg["eas_beta_ha"]), hetz=_f(tg["het_z"]), maf=_f(tg["eas_instr_maf"]),
                      dist_kb=abs(eas_pos - eur_pos) / 1000.0), "S23N + SuppTable_locus_attribution instr_pos")
    _need(F["tgfb1"]["eur"] > 0 > F["tgfb1"]["eas"], "TGFB1-CAD is no longer opposite in sign (EUR +, EAS -)")
    # the sign's robustness to the instrument tie rule (AUTHOR RESOLUTION 1). FIX round (verify:numbers): genome-wide
    # significance is the census's OWN rule, the REPORTED CKB p < eas__census_common.GWS (sig = pvalue < GWS), not |z|
    # recomputed from rounded beta/SE (which admitted three variants with reported p 5.3-5.8e-8). The candidate instruments
    # are the census's domain -- genome-wide significant, in the lead's +/-500 kb window, MAF >= 0.01 (n_domain) -- that are
    # present and harmonisable in the outcome GWAS (n_domain_shared; exactly the coloc-input rows they join).
    sys.path.insert(0, LB)
    import eas__census_common as ECC
    regf = sorted(f for f in os.listdir(CKBREG) if f.startswith("TGFB1__") and f.endswith(".tsv"))
    _need(len(regf) == 1, "not one TGFB1 CKB region file: %s" % regf)
    lo_, hi_ = int(tgr["coloc_lo_hg38"]), int(tgr["coloc_hi_hg38"])
    dom = [r for r in _rows(os.path.join(CKBREG, regf[0]), "\t")
           if _f(r["pvalue"]) is not None and _f(r["pvalue"]) < ECC.GWS and lo_ <= int(r["position"]) <= hi_
           and _f(r["af"]) is not None and min(_f(r["af"]), 1 - _f(r["af"])) >= ECC.MAF_MIN]
    ci = {r["snp"]: r for r in _rows(os.path.join(CIN, "TGFB1__CAD.tsv"), "\t")}
    shared = [ci[k] for k in ("%s:%s:%s:%s" % (r["chr"], r["position"], r["ref"], r["alt"]) for r in dom) if k in ci]
    _need((len(dom), len(shared)) == (int(tgr["n_domain"]), int(tgr["n_domain_shared"])),
          "TGFB1-CAD candidate instruments %d / %d != the census's n_domain / n_domain_shared %s / %s"
          % (len(dom), len(shared), tgr["n_domain"], tgr["n_domain_shared"]))
    wr = [_f(r["by"]) / _f(r["bx"]) for r in shared]
    put("tgfb1_strong", len(shared), "CKB region file: reported p < GWS, lead +/-500 kb, MAF >= 0.01, joined to "
        "coloc_inputs/TGFB1__CAD.tsv (== results_census n_domain_shared)")
    put("tgfb1_strong_neg", sum(1 for w in wr if w < 0), "coloc_inputs/TGFB1__CAD.tsv sign(by/bx)")
    put("tgfb1_strong_rng", (min(wr), max(wr)), "coloc_inputs/TGFB1__CAD.tsv by/bx")
    il = g[("IL6R", "CAD")]
    core = {(r["gene"], r["outcome"]): r for r in _rows(CORE)}
    put("il6r", dict(h4=_f(il["eas_coloc_H4"]), h4_6=_f(il["eas_coloc_H4_p12_1e6"]),
                     eur_h4_6=_f(core[("IL6R", "CAD")]["prior_pp_h4_1e6"]), eur=_f(il["eur_beta_ha"]),
                     eas=_f(il["eas_beta_ha"])), "S23N IL6R-CAD; SuppTable_coloc_robustness_core prior_pp_h4_1e6")
    _need(F["il6r"]["eur"] < 0 and F["il6r"]["eas"] < 0, "IL6R-CAD is no longer protective in both ancestries")
    # FIX round (verify:accuracy): TGFB1-CAD's EUROPEAN colocalization is prior-fragile, unlike the other two robust edges
    e6 = {k: _f(core[k]["prior_pp_h4_1e6"]) for k in (("CELSR2", "CAD"), ("FGF5", "eGFR"), ("TGFB1", "CAD"))}
    put("eur_1e6", e6, "SuppTable_coloc_robustness_core prior_pp_h4_1e6")
    _need(e6[("TGFB1", "CAD")] < 0.80 <= min(e6[("CELSR2", "CAD")], e6[("FGF5", "eGFR")]) and
          core[("TGFB1", "CAD")]["prior_stable"] == "FALSE",
          "TGFB1-CAD is no longer the one robust census edge whose European colocalization fails at p12 = 1e-6")
    ab = rc[("ABO", "T2D")]
    put("abo_t2d", dict(h4_6=float(ab["PP.H4_p12_1e-6"]), ratio_6=float(ab["PP.H4_p12_1e-6"]) / float(ab["PP.H3_p12_1e-6"])),
        "results_census.tsv ABO-T2D PP.H4/PP.H3 at p12 = 1e-6")
    _need(F["abo_t2d"]["h4_6"] >= 0.80 and F["abo_t2d"]["ratio_6"] < 5, "ABO-T2D no longer fails only the ratio at 1e-6")
    others = sorted(set(sup) - set(rob) - {pair("IL6R", "CAD"), pair("ABO", "T2D")})
    put("cen_fell_below", others, "S23N")
    _need(all(_f(g[tuple(x.split(u"–"))]["eas_coloc_H4_p12_1e6"]) < 0.80 for x in others), "a non-robust pair is >= 0.80")
    # direction of the colocalizing pairs with a strand-resolvable (non-palindromic) instrument
    palset = set(pal)
    sr = sorted(x for x in sup if tuple(x.split(u"–")) not in palset and x != pair("TGFB1", "CAD"))
    same = [x for x in sr if math.copysign(1, _f(g[tuple(x.split(u"–"))]["eur_beta_ha"]))
            == math.copysign(1, _f(g[tuple(x.split(u"–"))]["eas_beta_ha"]))]
    _need(same == sr, "a strand-resolvable colocalizing pair has the opposite sign: %s" % sorted(set(sr) - set(same)))
    put("cen_sr_same", sr, "S23N signs; results_census.tsv alleles")
    put("cen_pal_coloc", sorted(x for x in sup if tuple(x.split(u"–")) in palset), "results_census.tsv")
    # the non-colocalizing pairs
    pc = g[("PCSK9", "CAD")]
    put("pcsk9_hetp", _f(pc["het_p"]), "S23N PCSK9-CAD het_p")
    _need(math.copysign(1, _f(pc["eur_beta_ha"])) == math.copysign(1, _f(pc["eas_beta_ha"])), "PCSK9-CAD native sign flipped")
    ar = g[("ARG1", "T2D")]
    put("arg1_nat", dict(maf=_f(ar["eas_instr_maf"]), hetp=_f(ar["het_p"])), "S23N ARG1-T2D")
    _need(F["arg1_nat"]["maf"] < 0.05, "ARG1's native instrument is no longer low-frequency")
    ax = {(r["gene"], r["outcome"]): r for r in _rows(EAS_AXIS)}
    put("arg1_outcomeside", dict(state=ax[("ARG1", "T2D")]["eas_outcome_state"], beta=_f(ax[("ARG1", "T2D")]["beta_mr_eas"])),
        "SuppTable_EAS_axis ARG1-T2D")
    _need(F["arg1_outcomeside"]["state"] == "EAS_COMPATIBLE" and F["arg1_outcomeside"]["beta"] < 0,
          "ARG1-T2D is no longer outcome-side compatible and protective")
    # CKB-free sensitivity
    ck = [r for r in tst if r["eas_ckb_overlap"] == "True"]
    changed = sorted(P(r) for r in ck if r["eas_ckbfree_support"] != r["eas_support"])
    _need(changed == [pair("APOE", "T2D")], "the CKB-free call changes moved: %s" % changed)
    ap = g[("APOE", "T2D")]
    rcf = {(r["gene"], r["outcome"]): r for r in _rows(RCF, "\t")}
    put("apoe_ckbfree", dict(agen=_f(ap["eas_coloc_H4"]), bbj=_f(ap["eas_ckbfree_coloc_H4"]),
                             bbj_6=float(rcf[("APOE", "T2D")]["PP.H4_p12_1e-6"])), "S23N + results_ckbfree.tsv")
    _need(ap["eas_ckbfree_support_p12_1e6"] == "True", "APOE-T2D is not robust CKB-free")
    free_sup = set(sup) - {P(r) for r in ck if r["eas_ckbfree_support"] != "True"} | {P(r) for r in ck if r["eas_ckbfree_support"] == "True"}
    free_rob = set(rob) - {P(r) for r in ck if r["eas_ckbfree_support_p12_1e6"] != "True"} | \
        {P(r) for r in ck if r["eas_ckbfree_support_p12_1e6"] == "True"}
    put("ckbfree_totals", (len(free_sup), len(free_rob)), "S23N eas_ckbfree_support(_p12_1e6) with the CKB-free rows swapped in")
    _need((rc[("ARG1", "T2D")]["call"], rcf[("ARG1", "T2D")]["call"]) == ("NOT_SUPPORTED_H3", "NOT_SUPPORTED_H1"),
          "ARG1-T2D no longer moves H3 -> H1 CKB-free")
    # overlap with the outcome-side analysis (Figure 3a)
    both = [r for r in cen if key(r) in ax and _f(ax[key(r)]["beta_mr_eas"]) is not None and _f(r["eas_beta_ha"]) is not None]
    sd = sorted(P(r) for r in both if math.copysign(1, _f(ax[key(r)]["beta_mr_eas"])) != math.copysign(1, _f(r["eas_beta_ha"])))
    # both columns on the risk-aligned scale: the outcome side's direction_class agrees with sign(beta_mr_eas) x sign(EUR)
    for r in both:
        dc = ax[key(r)]["direction_class"]
        same_sign = math.copysign(1, _f(ax[key(r)]["beta_mr_eas"])) == math.copysign(1, _f(r["eur_beta_ha"]))
        _need(dc in ("SAME", "OPPOSITE") and (dc == "SAME") == same_sign,
              "outcome-side beta %s is not on the risk-aligned scale of S23N" % P(r))
    put("cen_x_outcomeside", len(both), "S23N x SuppTable_EAS_axis beta_mr_eas")
    put("cen_x_signdiff", sd, "S23N x SuppTable_EAS_axis")
    # the off-network pilot block
    put("pil_n", len(pil), "S23N set=pilot_offnetwork")
    put("pil_prot", len({r["protein"] for r in pil}), "S23N pilot proteins")
    psup = sorted(P(r) for r in pil if r["eas_support"] == "True")
    _need(psup == [pair("FGFR4", "T2D"), pair("UMOD", "eGFR")], "the pilot colocalizations moved: %s" % psup)
    put("pil_coloc", psup, "S23N pilot eas_support")
    _need(all(r["eas_support_p12_1e6"] == "True" for r in pil if r["eas_support"] == "True"), "a pilot colocalization is not robust")
    fg = g[("FGFR4", "T2D")]
    fz = abs(_f(fg["eur_beta_ha"]) / _f(fg["eur_se"]))
    put("fgfr4", dict(h4=_f(fg["eas_coloc_H4"]), h4_6=_f(fg["eas_coloc_H4_p12_1e6"]), eas=_f(fg["eas_beta_ha"]),
                      easp=_f(fg["eas_p"]), eurp=math.erfc(fz / math.sqrt(2)), hetp=math.erfc(abs(_f(fg["het_z"])) / math.sqrt(2)),
                      ckbfree=_f(fg["eas_ckbfree_coloc_H4"])), "S23N FGFR4-T2D (European p and het p from beta/se and het_z)")
    _need(fg["eur_call"] == "" and F["fgfr4"]["eurp"] > 0.05 and fg["eas_ckbfree_support"] == "True",
          "FGFR4-T2D: European colocalization now tested, European MR significant, or CKB-free not supported")
    um = g[("UMOD", "eGFR")]
    put("umod", dict(h4=_f(um["eas_coloc_H4"]), h4_6=_f(um["eas_coloc_H4_p12_1e6"])), "S23N UMOD-eGFR")
    _need(um["eur_call"] == "H3", "UMOD-eGFR is no longer a European coloc.abf H3 pair")
    # FIX round (verify:accuracy): FGFR4-T2D was SELECTED on East Asian data -- the pilot row is the discovery row of the
    # East Asian T2D scan (CKB instruments x AGEN), where it was the one East-Asian-specific survivor
    scan = _rows(SCAN, "\t")
    led = [r for r in _rows(LEDGER, "\t") if not r["protein"].startswith("#")]
    s4 = [r for r in scan if r["protein"] == "FGFR4"]
    _need(len(s4) == 1 and s4[0]["rsid"] == fg["eas_instr"] and abs(float(s4[0]["F"]) - _f(fg["eas_F"])) < 0.05
          and abs(float(s4[0]["eas_mr"]) - _f(fg["eas_beta_ha"])) < 5e-5 and abs(float(s4[0]["eas_p"]) / _f(fg["eas_p"]) - 1) < 1e-6,
          "the pilot FGFR4-T2D row is not the scan's discovery row")
    _need(float(s4[0]["eas_q"]) < 0.05 and float(s4[0]["eur_p"]) > 0.05 and float(s4[0]["het_p"]) < 0.05 and s4[0]["strict"] == "True",
          "FGFR4 no longer passes the scan's East-Asian-specificity filters")
    lt = {r["protein"]: r for r in led if r["outcome"] == "T2D"}
    strict = sorted(r["protein"] for r in scan if r["strict"] == "True")
    _need(set(strict) <= set(lt) and [p for p in strict if _f(lt[p]["coloc_H4_p12_1e6"]) >= 0.80] == ["FGFR4"]
          and sorted(p for p, r in lt.items() if r["verdict"] == "GENUINE_EAS_specific") == ["FGFR4"],
          "FGFR4 is no longer the scan's only East-Asian-specific survivor")
    _need(len({r["protein"] for r in scan}) == len(scan), "a scan protein is duplicated")
    # the text names the COMMON-instrument condition: a protein meeting the significance / null-European / heterogeneity
    # filters on a rare instrument (GNLY, EAS AF 0.045) failed the scan's 'common' filter and was never colocalized
    three = [r for r in scan if float(r["eas_q"]) < 0.05 and float(r["eur_p"]) > 0.05 and float(r["het_p"]) < 0.05]
    _need(all(r["pleio"] == "False" for r in three) and
          sorted(r["protein"] for r in three if r["strict"] == "True") == strict and
          all(r["strict"] == "True" or not (0.05 < float(r["eas_af"]) < 0.95) for r in three),
          "a scan protein meeting the three filters failed 'strict' for a reason other than a rare instrument")
    put("fgfr4_scan_rare", sorted(r["protein"] for r in three if r["strict"] != "True"), "phase2_T2D_scan.tsv (rare-instrument exclusions)")
    put("fgfr4_scan_n", len(scan), "paper 3.eas/analysis/results/phase2_T2D_scan.tsv (proteins tested against AGEN T2D)")
    _need("BBJ" in fg["eas_outcome_resource"] and "AGEN" in fg["eas_outcome_resource"],
          "the FGFR4 outcome resource no longer records that AGEN includes Biobank Japan")
    lg = io.open(PLINKLOG, encoding="utf-8").read()
    _need("--ld rs451643 rs351855" in lg and "--bfile EAS/1000G_EAS" in lg and "504 samples" in lg and fg["eas_instr"] == "rs451643",
          "plink2.log is not the rs451643 / rs351855 LD run on 1000 Genomes East Asians, or rs451643 is not the instrument")
    m = re.search(r"r\^2 = ([0-9.]+)", lg)
    _need(m, "plink2.log carries no r^2")
    put("fgfr4_g388r_r2", float(m.group(1)), "paper 3.eas/analysis/plink2.log (plink2 --ld rs451643 rs351855, 1000G EAS, 504)")
    # FIX round 4 (verify:accuracy-narrative): the scan judged the East-Asian-specific pattern against DIAMANTE (its
    # European T2D comparator), while the European estimate the text prints is the V4 MR screen's, on FinnGen R10 T2D --
    # a different European GWAS. Both read from their own records.
    mr4 = [r for r in _rows(MRS4) if (r["gene"], r["outcome_name"]) == ("FGFR4", "T2D") and r["method"] == "Wald ratio"]
    msr4 = io.open(os.path.join(os.path.dirname(EAS3), "Paper3.EAS_Manuscript_Interim_v1.md"), encoding="utf-8").read()
    irep4 = io.open(os.path.join(os.path.dirname(EAS3), "INTERIM_REPORT_v1.md"), encoding="utf-8").read()
    _need(len(mr4) == 1 and mr4[0]["outcome_id"] == "FinnGen_R10_T2D" and "T2D, DIAMANTE" in msr4
          and "T2D: AGEN N=433k vs DIAMANTE" in irep4 and abs(float(s4[0]["eur_p"]) - F["fgfr4"]["eurp"]) > 0.01,
          "the scan's European T2D comparator (DIAMANTE) and the printed European estimate's outcome (FinnGen R10) are no "
          "longer two different GWAS")
    put("fgfr4_eur_src", dict(printed=mr4[0]["outcome_id"], scan="DIAMANTE", scan_eurp=float(s4[0]["eur_p"])),
        "V4 phase2_mr_screen FGFR4 x T2D outcome_id; paper 3.eas records (the scan's European comparator)")
    # FIX round 2 (verify:numbers/narrative; the brief's pre-audit): the WHOLE East Asian scan FGFR4-T2D was selected
    # from, not its T2D arm alone. Four arms kept per-protein results (phase2_{T2D,CAD,Stroke,egfr}_scan.tsv; one row per
    # protein with both an East Asian and a European outcome row); heart failure and atrial fibrillation were screened on
    # the East Asian side too (inputs data/eas_outcomes/{hf,af}_EAS.tsv; INTERIM_REPORT_v1.md records both screens and
    # their empty yield) but kept no per-protein output and contributed no candidate to the colocalization ledger.
    arms = [("T2D", "phase2_T2D_scan.tsv"), ("CAD", "phase2_CAD_scan.tsv"), ("Stroke", "phase2_Stroke_scan.tsv"),
            ("eGFR", "phase2_egfr_scan.tsv")]
    armrows = {}
    for o_, fn in arms:
        rows_ = _rows(os.path.join(EAS3, "results", fn), "\t")
        _need(len({r["protein"] for r in rows_}) == len(rows_), "a %s scan protein is duplicated" % o_)
        armrows[o_] = rows_
    _need(not [x for x in os.listdir(os.path.join(EAS3, "results")) if re.match(r"(?i)phase2_(hf|af)", x)],
          "an HF/AF scan output exists now: describe its per-protein results")
    _need(all(os.path.isfile(os.path.join(EAS3, "data", "eas_outcomes", x)) for x in ("hf_EAS.tsv", "af_EAS.tsv")),
          "the East Asian HF/AF screening inputs moved")
    irep = io.open(os.path.join(os.path.dirname(EAS3), "INTERIM_REPORT_v1.md"), encoding="utf-8").read()
    _need("HF: Enzan vs HERMES" in irep and "AF: BBJ vs Nielsen" in irep and "HF 0 (" in irep and "AF 0," in irep,
          "the interim report no longer records the HF and AF screens and their empty yield")
    _need({r["outcome"] for r in led} <= {o_ for o_, _ in arms}, "a colocalization candidate comes from an arm without "
          "retained results: %s" % sorted({r["outcome"] for r in led}))
    refined = {r["protein"]: r for r in _rows(os.path.join(EAS3, "results", "phase2_egfr_refined.tsv"), "\t")}
    strict_of = lambda r: (refined[r["protein"]]["strict_EAS_specific"] if r["outcome"] == "eGFR" else
                           [x for x in armrows[r["outcome"]] if x["protein"] == r["protein"]][0]["strict"]) == "True"
    _need(all(r["protein"] in {x["protein"] for x in armrows[r["outcome"]]} for r in led), "a ledger candidate is not a scan row")
    dflt = [r for r in led if _f(r["coloc_H4_p12_1e5"]) >= 0.80]
    cons = [r for r in led if _f(r["coloc_H4_p12_1e6"]) >= 0.80]
    _need([(r["protein"], r["outcome"]) for r in cons] == [("FGFR4", "T2D")] and set(map(id, cons)) <= set(map(id, dflt)),
          "FGFR4-T2D is no longer the one scan candidate that keeps PP.H4 >= 0.80 at the conservative prior")
    pleio = re.findall(r'\("\d+",\d+,\d+\)', io.open(os.path.join(EAS3, "src", "phase2_binary_scan.py"), encoding="utf-8").read())
    _need(len(pleio) >= 1, "the scan's pleiotropic-locus list is not where it was")
    put("fgfr4_scan", dict(arms={o_: len(v) for o_, v in armrows.items()}, tests=sum(len(v) for v in armrows.values()),
                           prot=len({r["protein"] for v in armrows.values() for r in v}), cand=len(led),
                           strict=sum(1 for r in led if strict_of(r)), dflt=len(dflt)),
        "paper 3.eas/analysis/results phase2_{T2D,CAD,Stroke,egfr}_scan.tsv, phase2_egfr_refined.tsv, phase2_coloc_ledger.tsv "
        "(PP.H4 only; the ledger has no PP.H3); HF/AF: data/eas_outcomes + INTERIM_REPORT_v1.md")
    # FIX round 3 (verify:narrative/numbers, wf_876cfb67-e41): (a) the heart-failure screen DID name a candidate, CDKN1A,
    # shared with Europeans (non-significant cross-ancestry heterogeneity) -- 'neither contributed a candidate' was false;
    # it never entered colocalization (not in the ledger). Read from the project's own records, which are the only place
    # the unsaved HF/AF screens are described. (b) WHY the five non-East-Asian-specific eGFR pairs entered colocalization:
    # they are exactly the eGFR scan's false-discovery-significant hits absent from the European eGFR hit table (its
    # 'EUR-/EAS+' quadrant), all failing the strict filter on the refined re-analysis, and the records say all candidates
    # were carried into colocalization for transparency.
    mrep = io.open(os.path.join(os.path.dirname(EAS3), "Paper3.EAS_Manuscript_Interim_v1.md"), encoding="utf-8").read()
    _need(u"HF 0 (CDKN1A → shared, non-significant heterogeneity)" in irep and
          "One heart-failure candidate (CDKN1A) was a shared, not East-Asian-specific, signal (non-significant cross-ancestry "
          "heterogeneity" in mrep and "CDKN1A" not in {r["protein"] for r in led},
          "the records no longer name CDKN1A as the heart-failure screen's one shared (not East-Asian-specific) candidate, "
          "or CDKN1A entered colocalization")
    put("hf_candidate", "CDKN1A", "paper 3.eas INTERIM_REPORT_v1.md ('HF 0 (CDKN1A -> shared, non-significant heterogeneity)') "
        "+ Paper3.EAS_Manuscript_Interim_v1.md; not in phase2_coloc_ledger.tsv")
    eg_scan = armrows["eGFR"]
    eg_hits = sorted(r["protein"] for r in eg_scan if float(r["mr_q"]) < 0.05 and r["in_EUR_eGFR"] == "False")
    eg_led = sorted(r["protein"] for r in led if r["outcome"] == "eGFR")
    _need(eg_hits == eg_led and all(refined[p]["strict_EAS_specific"] == "False" for p in eg_led)
          and sorted(r["protein"] for r in eg_scan if r["quadrant"].startswith("EUR-/EAS+")) == eg_led
          and "all eight were carried into colocalization for transparency" in mrep
          and len(eg_led) + sum(1 for r in led if strict_of(r)) == len(led),
          "the ledger's eGFR candidates are no longer exactly the eGFR scan's false-discovery-significant hits absent from "
          "the European eGFR hit table, all failing the strict filter and carried for transparency")
    put("scan_egfr_carried", len(eg_led), "phase2_egfr_scan.tsv (mr_q < 0.05, in_EUR_eGFR False) == ledger eGFR rows; "
        "phase2_egfr_refined.tsv strict_EAS_specific False; Paper3.EAS_Manuscript_Interim_v1.md 'carried ... for transparency'")
    # the pilot instrument was chosen by file order among the rsID-bearing variants tied at the CKB p floor; the FGFR4
    # missense variant Gly388Arg (rs351855) is one of them (S23 eas_n_tied_top; the pilot region file itself)
    ex = _rows(os.path.join(EAS3, "data", "eas_exposure", "FGFR4_cis_eas.tsv"), "\t")
    rsx = [r for r in ex if r["rsid"].startswith("rs")]
    pmin = min(float(r["pvalue"]) for r in rsx)
    tied = [r["rsid"] for r in rsx if float(r["pvalue"]) == pmin]
    _need(pmin <= 1e-300 and pmin == min(float(r["pvalue"]) for r in ex) and len(tied) == int(fg["eas_n_tied_top"])
          and tied[0] == fg["eas_instr"] and "rs351855" in tied and "file order" in fg["eas_instr_rule"],
          "the FGFR4 instrument is no longer the first of the rsID-bearing variants tied at the CKB p floor, with rs351855 "
          "among them")
    put("fgfr4_tied", len(tied), "paper 3.eas/analysis/data/eas_exposure/FGFR4_cis_eas.tsv (rsID-bearing rows at p = 1e-300) "
        "== S23N eas_n_tied_top")
    # FIX round (verify:underselling): the CKB-free sensitivity kept every colocalization and added one
    lost, gained = sorted(set(sup) - free_sup), sorted(free_sup - set(sup))
    _need(lost == [] and gained == [pair("APOE", "T2D")], "the CKB-free re-run no longer keeps all and adds APOE-T2D")
    put("ckbfree_kept", len(sup), "S23N: every default-resource colocalization is also CKB-free supported")
    put("ckbfree_gained", gained, "S23N eas_ckbfree_support")

    # ---- D3 / D1: the metabolite layer -------------------------------------------------------------------------------
    node = _rows(NODE)
    put("node_n", len(node), "node_census_v3.csv")
    put("node_genes", len({r["gene"] for r in node}), "node_census_v3.csv")
    nn = sorted(u"%s → %s" % (r["gene"], "type 2 diabetes" if r["node_outcome"] == "T2D" else r["node_outcome"])
                for r in node if r["verdict"].startswith("V3_NODE"))
    _need(nn == [u"APOE → type 2 diabetes", u"ARG1 → type 2 diabetes"], "the admitted node pairs moved: %s" % nn)
    put("node_admitted", nn, "node_census_v3.csv verdict")
    ncol = _rows(NODEC, "\t")
    apo = [r for r in ncol if r["gene"] == "APOE" and r["outcome"] == "T2D"]
    _need(len(apo) == 2 and all(_f(r["PP3_N"]) > 0.995 for r in apo), "APOE's node metabolites are no longer distinct-variant")
    _need(sorted(r["metabolite"][:8] for r in apo) == ["1-(1-eny", "sphingom"], "APOE's node metabolites changed")
    put("apoe_node_h3", min(_f(r["PP3_N"]) for r in apo), "node_coloc side_by_side PP3_N")
    ag = [r for r in ncol if (r["gene"], r["metabolite"], r["outcome"]) == ("ARG1", "arginine", "T2D")]
    _need(len(ag) == 1, "no single ARG1 arginine-T2D node row")
    put("arg_t2d", dict(h4=_f(ag[0]["PP4_N"]), h4_5=_f(ag[0]["PP4_p12_5e-6"]), h4_6=_f(ag[0]["PP4_p12_1e-6"]),
                        nsnp=int(ag[0]["nsnp"])), "node_coloc side_by_side (FinnGen R10 T2D)")
    _need(ag[0]["dual_rule_N"] == "TRUE" and ag[0]["dual_1e-6"] == "FALSE", "arginine-T2D: dual rule at default / 1e-6 moved")
    put("arg1_t2d", dict(h4=_f(core[("ARG1", "T2D")]["pp_h4_abf"]), h4_6=_f(core[("ARG1", "T2D")]["prior_pp_h4_1e6"])),
        "SuppTable_coloc_robustness_core ARG1-T2D")
    edges = [r for r in _rows(EDGES) if (r["from"], r["to"]) == ("ARG1", "T2D")]
    _need(len(edges) == 1, "no single ARG1-T2D network edge")
    put("arg1_t2d_mr", dict(b=_f(edges[0]["beta"]), p=_f(edges[0]["pval"])), "SuppTable_network_edges ARG1-T2D")
    bmi = [r for r in _rows(INTER) if (r["gene"], r["outcome_name"]) == ("ARG1", "BMI")]
    _need(len(bmi) == 1 and bmi[0]["fdr_sig"] == "TRUE", "ARG1 -> BMI is not a single FDR-significant row")
    put("arg1_bmi", dict(b=_f(bmi[0]["b"]), p=_f(bmi[0]["pval"])), "phase6 intermediate MR ARG1 -> BMI")
    # the ARG1 locus tracks (Figure 8 source)
    tr = {}
    for r in _rows(F08):
        tr.setdefault(r["track"], []).append((float(r["pos_mb"]), float(r["mlog10p"])))
    tt = {k: max(v, key=lambda x: x[1]) for k, v in tr.items()}
    t2d = [k for k in tr if k.startswith("type 2 diabetes")][0]
    arg = [k for k in tr if k.startswith("arginine")][0]
    prt = [k for k in tr if k.startswith("ARG1")][0]
    win = [p for p, l in tr[t2d] if l >= tt[t2d][1] - 1.0]
    put("arg_locus", dict(arg_prot_bp=abs(tt[arg][0] - tt[prt][0]) * 1e6, t2d_arg_kb=abs(tt[t2d][0] - tt[arg][0]) * 1e3,
                          t2d_span_kb=(max(win) - min(win)) * 1e3), "Figure_08 source (track maxima; T2D within 1 log10 unit)")
    # the bridge screen
    bsr = _rows(BSR)
    put("bsr_n", len(bsr), "bridge_screen_results.csv")
    put("bsr_prot", len({r["protein"] for r in bsr}), "bridge_screen_results.csv")
    put("bsr_metab", len({r["metabolite"] for r in bsr}), "bridge_screen_results.csv")
    _need(all(_f(r["H4"]) is not None for r in bsr) and len({(r["protein"], r["metabolite"]) for r in bsr}) == len(bsr),
          "a screen row is untested or duplicated")
    _need({r["protein_arm"] for r in bsr} == {"UKB-PPP_EUR_discovery"}, "a screen arm is not the UKB-PPP European discovery arm")
    br = [r for r in bsr if r["support"] == "True"]
    put("bridges", len(br), "bridge_screen_results.csv support")
    put("bridge_prot", len({r["protein"] for r in br}), "bridge_screen_results.csv support")
    top = sorted(br, key=lambda r: -_f(r["H4"]))
    put("bridge_top2", [(r["protein"], r["metabolite"], _f(r["H4"])) for r in top[:2]], "bridge_screen_results.csv H4")
    _need([(r["protein"], r["metabolite"]) for r in top[:2]] == [("DPEP1", "cysteinylglycine"), ("ARG1", "arginine")],
          "the two strongest bridges moved")
    rng = lambda prot: (min(_f(r["H4"]) for r in br if r["protein"] == prot), max(_f(r["H4"]) for r in br if r["protein"] == prot))
    put("bridge_rng", {p: rng(p) for p in ("BRAP", "SMPD3", "GSTA3")}, "bridge_screen_results.csv H4")
    orn = [r for r in bsr if (r["protein"], r["metabolite"]) == ("ARG1", "ornithine")]
    _need(len(orn) == 1 and orn[0]["support"] == "False", "ARG1-ornithine is not a single non-bridge screen row")
    put("orn_screen", dict(h4=_f(orn[0]["H4"]), h3=_f(orn[0]["H3"])), "bridge_screen_results.csv ARG1-ornithine")
    # the class table
    bcs = _rows(BCS)
    put("bcs_n", len(bcs), "bridge_class_summary.csv")
    cs = {(r["protein"], r["metabolite"]): r for r in bcs}
    t1 = sorted(k for k, r in cs.items() if r["tier"] == "tier 1")
    t2 = sorted(k for k, r in cs.items() if r["tier"] == "tier 2")
    _need(t1 == [("ARG1", "arginine"), ("DPEP1", "cysteinylglycine")], "bridge tier 1 moved: %s" % t1)
    _need(len(t1) + len(t2) == F["bridges"] and set(t1 + t2) == {(r["protein"], r["metabolite"]) for r in br},
          "the tiered pairs are not the screen bridges")
    put("tier1", t1, "bridge_class_summary.csv tier")
    put("tier2_n", len(t2), "bridge_class_summary.csv tier")
    _need(all(cs[k]["robust"] == "True" and cs[k]["modal_cells"] == "48" for k in t2), "a bridge-tier-2 pair is not one class in 48/48")
    _need({r["n_cells"] for r in bcs} == {"48"}, "the grid is not 48 cells")
    a1 = cs[("ARG1", "arginine")]
    _need((a1["n_gene_attributed"], a1["glist_n_gene_attributed"], a1["ens116all_n_gene_attributed"]) == ("48", "48", "48"),
          "ARG1-arginine is no longer gene-attributed in 48/48 under every annotation")
    put("arg1_ratio", dict(lead=_f(a1["ratio_cis_lead"]), second=_f(a1["ratio_shared_drop"])), "bridge_class_summary.csv ARG1-arginine")
    _need(a1["direction"] == "negative" and a1["shared_signal_drop"] == "B", "ARG1-arginine direction / shared signal moved")
    grid = _rows(GRID)
    ga = [r for r in grid if (r["protein"], r["metabolite"]) == ("ARG1", "arginine")]
    _need({r["signal_A"] for r in ga} == {"rs2781668"}, "ARG1's lead signal (A) is no longer rs2781668")
    put("arg1_lead_rsid", "rs2781668", "bridge_class_robustness__ens116pc.csv signal_A; = SuppTable_EAS_axis instrument_rsID")
    _need(ax[("ARG1", "T2D")]["instrument_rsID"] == "rs2781668", "the ARG1 MR instrument is not the fine-mapping lead")
    on = cs[("ARG1", "ornithine")]
    _need(on["modal_class"] == "shared variant, mediation uncertain" and on["robust"] == "True" and on["direction"] == "positive",
          "ARG1-ornithine is no longer a robust shared-variant, positive-direction pair")
    dp = cs[("DPEP1", "cysteinylglycine")]
    gd = [r for r in grid if (r["protein"], r["metabolite"]) == ("DPEP1", "cysteinylglycine")]
    _need(len(gd) == 48 and all(r["SHARED"] == "True" and r["PROP"] == "True" and r["RATIO_applicable"] == "True" for r in gd),
          "DPEP1-cysteinylglycine: shared signal / proportionality across two signals no longer holds in 48/48")
    _need(all("CHMP1A" in r["rivals_full"] for r in gd), "CHMP1A is no longer the named rival at the DPEP1 locus")
    _need(dp["direction"] == "negative" and dp["modal_class"].startswith("locus-level") and dp["robust"] == "True",
          "DPEP1-cysteinylglycine direction / locus-level class moved")
    put("dpep1_attr", dict(glist=int(dp["glist_n_gene_attributed"]), ens=int(dp["n_gene_attributed"]),
                           ensall=int(dp["ens116all_n_gene_attributed"])), "bridge_class_summary.csv")
    _need(F["dpep1_attr"]["glist"] > F["dpep1_attr"]["ens"], "DPEP1 attribution is no longer annotation-dependent")
    # FIX round (verify:numbers/accuracy): WHY, from the registered grids and the class summary's pinned inputs. In the
    # palindrome-drop half the A2 signal and lead is one variant inside VPS9D1; in the frequency-resolved half the lead
    # falls inside CHMP1A under Ensembl 116 and is attributed under glist-hg19 except where the CHMP1A protein blocks.
    gl = [r for r in _rows(GLIST) if (r["protein"], r["metabolite"]) == ("DPEP1", "cysteinylglycine")]
    hv = {h: [r for r in gd if r["pal"] == h] for h in ("drop", "freq")}
    hg = {h: [r for r in gl if r["pal"] == h] for h in ("drop", "freq")}
    _need(len(gl) == 48 and all(len(v) == 24 for v in list(hv.values()) + list(hg.values())), "the DPEP1 grid halves are not 24 + 24")
    sig = {h: {(r["A2_signal_snp"], r["A2_lead_snp"]) for r in hv[h] + hg[h]} for h in hv}
    _need(all(len(v) == 1 for v in sig.values()), "a half of the DPEP1 grid carries more than one A2 signal / lead")
    with gzip.open(HARM_DPEP1, "rt", encoding="utf-8") as fh:
        harm = {(r["chr"], r["pos37"]): r for r in csv.DictReader(fh, delimiter="\t")}
    egn = [r for r in _rows(ENSGENES, "\t") if r["biotype"] == "protein_coding"]
    inside = lambda ch, p: sorted({r["symbol"] for r in egn if r["chr"] == ch and int(r["start38"]) <= p <= int(r["end38"])})
    where = {}
    for h in ("drop", "freq"):
        snp, lead = list(sig[h])[0]
        ch, p37 = lead.split(":")[:2]
        hr = harm.get((ch, p37))
        _need(hr is not None, "the DPEP1 %s-half lead %s is not in the harmonised table" % (h, lead))
        where[h] = dict(snp=snp, lead_rsid=hr["rsid_ref"], genes=inside(ch, int(hr["pos38"])))
    _need(where["drop"]["genes"] == ["VPS9D1"] and where["drop"]["lead_rsid"] == where["drop"]["snp"],
          "the palindrome-drop half's DPEP1 signal is no longer one variant inside VPS9D1: %s" % where["drop"])
    _need(where["freq"]["genes"] == ["CHMP1A"], "the frequency-resolved half's DPEP1 lead no longer falls inside CHMP1A: %s" % where["freq"])
    ga = lambda rows: sum(1 for r in rows if r["class"] == "gene-attributed bridge")
    blk = [r for r in hg["freq"] if r["class"] != "gene-attributed bridge"]
    _need((ga(hg["drop"]), ga(hv["drop"]), ga(hv["freq"])) == (0, 0, 0) and ga(hg["freq"]) == F["dpep1_attr"]["glist"]
          and all(r["A3"] == "False" and "CHMP1A:" in r["rivals_full"] and "BLOCKS" in r["rivals_full"] for r in blk),
          "the DPEP1 attribution no longer splits as VPS9D1 half / CHMP1A half")
    # FIX round 2 (verify:narrative/numbers): the CAUSE of the split. DPEP1's second forward-selected signal (grid signal_B)
    # differs by half: where frequency-resolvable palindromes are kept it is a palindromic variant inside DPEP1 and the
    # metabolite shares DPEP1's LEAD signal (A); where palindromes are dropped that variant is gone, the second signal is
    # the VPS9D1-located DPEP1 pQTL variant, and the metabolite shares THAT signal (B). The CHMP1A-blocked glist cells are
    # exactly the cells that use SuSiE as fitted (susie == asnow). Registered grids + the class summary's pinned inputs.
    sb = {h: {r["signal_B"] for r in hv[h] + hg[h]} for h in hv}
    sa = {h: {r["signal_A"] for r in hv[h] + hg[h]} for h in hv}
    sh = {h: {r["cond_signal_full"] for r in hv[h] + hg[h]} for h in hv}
    _need(all(len(sb[h]) == 1 and len(sa[h]) == 1 and len(sh[h]) == 1 for h in hv) and sa["drop"] == sa["freq"],
          "the DPEP1 grid halves no longer have one signal A (shared), one signal B and one shared signal each")
    b_drop, b_freq = list(sb["drop"])[0], list(sb["freq"])[0]
    _need(sh == {"drop": {"B"}, "freq": {"A"}} and where["drop"]["snp"] == b_drop and where["freq"]["snp"] == list(sa["freq"])[0],
          "the metabolite no longer shares DPEP1's second signal where palindromes are dropped and its lead signal where "
          "they are kept")
    hb = [r for r in harm.values() if r["rsid_ref"] == b_freq]
    _need(len(hb) == 1 and hb[0]["palindromic"] in ("True", "TRUE", "1") and {hb[0]["a1"], hb[0]["a2"]} in ({"A", "T"}, {"C", "G"})
          and inside(hb[0]["chr"], int(hb[0]["pos38"])) == ["DPEP1"] and b_freq not in {r["signal_A"] for r in hv["drop"]},
          "DPEP1's second signal where palindromes are kept is no longer one palindromic variant inside DPEP1")
    _need(sorted(r["susie"] for r in blk) == ["asnow"] * len(blk) and
          len([r for r in hg["freq"] if r["susie"] == "asnow"]) == len(blk),
          "the CHMP1A-blocked cells are no longer exactly the cells that use SuSiE as fitted")
    put("dpep1_halves", dict(n=24, drop_gene="VPS9D1", freq_gene="CHMP1A", glist_freq=ga(hg["freq"]), blocked=len(blk),
                             b_drop=b_drop, b_freq=b_freq),
        "bridge_class_robustness.csv (glist) + __ens116pc.csv (signal_A/B, cond_signal_full, susie), harmonised pos38 and "
        "palindromic flag, Ensembl 116 protein-coding gene table")
    # FIX round 3 (verify:numbers/narrative): WHERE the two halves' signals sit, so the Ensembl-116 non-attribution is not
    # undersold -- the palindrome-kept half's shared-signal lead lies just past DPEP1's 3' end (DPEP1 is on the + strand)
    # and only just inside CHMP1A's annotated span; the palindrome-drop half's signal B is ~80 kb from DPEP1 in VPS9D1
    # (not DPEP1's neighbour: four protein-coding genes lie between). Dropping palindromic variants removes rs1126464
    # itself (it is palindromic, and it is signal A or B in no palindrome-drop cell).
    gpos = {r["symbol"]: (int(r["start38"]), int(r["end38"]), r["strand"]) for r in egn if r["chr"] == "16"}
    lead_f = harm[tuple(list(sig["freq"])[0][1].split(":")[:2])]
    lp = int(lead_f["pos38"])
    d_start, d_end, d_str = gpos["DPEP1"]
    c_start, c_end, _ = gpos["CHMP1A"]
    b_row = [r for r in harm.values() if r["rsid_ref"] == b_drop]
    _need(len(b_row) == 1, "DPEP1's palindrome-drop signal B is not one row of the harmonised table")
    bp = int(b_row[0]["pos38"])
    between = sorted(s for s, (a_, e_, _) in gpos.items() if a_ > d_end and e_ < bp)
    _need(d_str == "1" and lp > d_end and c_start <= lp <= c_end and inside("16", lp) == ["CHMP1A"] and bp > d_end
          and inside("16", bp) == ["VPS9D1"] and len(between) >= 2
          and b_freq not in {r["signal_A"] for r in hv["drop"] + hg["drop"]} | {r["signal_B"] for r in hv["drop"] + hg["drop"]},
          "the DPEP1 halves no longer sit where the text says: the kept half's lead past DPEP1's 3' end inside CHMP1A's span, "
          "the dropped half's signal in VPS9D1 beyond non-neighbouring genes, rs1126464 absent from the drop half")
    put("dpep1_pos", dict(lead_rsid=lead_f["rsid_ref"], past_dpep1_bp=lp - d_end, in_chmp1a_bp=lp - c_start,
                          b_drop_from_dpep1_bp=bp - d_end, between=between),
        "harmonised.tsv.gz pos38 (the class summary's pinned input) vs raw/l3-finemap/ensembl116/ensembl116_genes.tsv "
        "(protein-coding start38/end38/strand)")
    # FIX round 4 (verify:accuracy-narrative, MAJOR undersell): the registered DPEP1-locus annotation table (Ensembl 116)
    # says what '559 bp inside CHMP1A' hid -- the shared-variant lead lies between the genes, inside CHMP1A's full gene
    # span only through a transcript with no coding sequence, and the nearer gene flips with the transcript set; DPEP1's
    # signal A peaks inside DPEP1; the palindrome-kept half's signal B is a DPEP1 missense variant. Every value is read
    # from the table and cross-checked here against the grid and the harmonised / gene-table derivation above.
    ann = {r["rsid"]: r for r in _rows(DANN)}
    sa_rs = list(sa["freq"])[0]
    _need(sorted(ann) == sorted({lead_f["rsid_ref"], sa_rs, b_freq, b_drop}),
          "the annotation table's variants %s are not the grid's DPEP1 lead / signal A / signal B variants" % sorted(ann))
    _need(all(r["release"] == "116" and r["assembly"] == "GRCh38" for r in ann.values()), "the annotation table is not "
          "Ensembl 116 / GRCh38")
    for rs_, r_ in ann.items():
        hh = [h for h in harm.values() if h["rsid_ref"] == rs_]
        _need(len(hh) == 1 and hh[0]["pos38"] == r_["pos38"], "%s: annotation-table position != harmonised pos38" % rs_)
    la = ann[lead_f["rsid_ref"]]
    ci = lambda r, k: int(r[k])
    _need(ci(la, "DPEP1_gene_bp") == lp - d_end and ci(la, "CHMP1A_gene_bp") == 0 and ci(ann[b_drop], "DPEP1_gene_bp")
          == bp - d_end, "the annotation table disagrees with the gene-table derivation of the DPEP1 positions")
    ctx = [x.strip() for x in la["CHMP1A_containing_tx"].split(";") if x.strip()]
    _need(len(ctx) == 1 and "(protein_coding_CDS_not_defined, not GENCODE basic)" in ctx[0] and not la["DPEP1_containing_tx"],
          "the shared-variant lead is no longer inside CHMP1A through exactly one transcript with no coding sequence: %s" % ctx)
    pc_d, pc_c = ci(la, "DPEP1_pc_tx_bp"), ci(la, "CHMP1A_pc_tx_bp")
    can_d, can_c = ci(la, "DPEP1_canonical_bp"), ci(la, "CHMP1A_canonical_bp")
    _need(pc_c < pc_d and can_d < can_c and round(pc_d / 100.0) == round(pc_c / 100.0) == 29,
          "the nearer gene no longer flips with the transcript set (CHMP1A over protein-coding transcripts, DPEP1 over "
          "the canonical ones), or the lead is no longer about 2.9 kb from both genes' protein-coding transcripts")
    ra = ann[sa_rs]
    _need(ci(ra, "DPEP1_gene_bp") == ci(ra, "DPEP1_canonical_bp") == 0 and
          re.search(r"(^|; )DPEP1 ENST\d+ intron_variant( |;|$)", ra["vep_canonical"]),
          "DPEP1's signal A is no longer intronic in DPEP1's canonical transcript: %s" % ra["vep_canonical"])
    rb_ = ann[b_freq]
    mm = re.search(r"(^|; )DPEP1 ENST\d+ missense_variant p\.([A-Z][a-z]{2}\d+[A-Z][a-z]{2})( MANE)?(;|$)", rb_["vep_canonical"])
    _need(mm and rb_["study_palindromic"] == "True", "the palindrome-kept half's signal B is no longer a palindromic DPEP1 "
          "missense variant: %s" % rb_["vep_canonical"])
    rd = ann[b_drop]
    _need(re.search(r"(^|; )%s ENST\d+ intron_variant" % where["drop"]["genes"][0], rd["vep_canonical"]), "the "
          "palindrome-drop half's signal is no longer intronic in %s: %s" % (where["drop"]["genes"][0], rd["vep_canonical"]))
    put("dpep1_ann", dict(lead=lead_f["rsid_ref"], sa=sa_rs, pc_kb=pc_d / 1000.0, m_pc=pc_d - pc_c, m_can=can_c - can_d,
                          aa=mm.group(2), cds=(ci(la, "DPEP1_cds_bp"), ci(la, "CHMP1A_cds_bp"))),
        "l3_metabolite/dpep1_locus_annotation.csv (Ensembl 116; l3fm__dpep1_locus.py), cross-checked against the grid, "
        "the harmonised pos38 and the class rules' gene table")
    # FIX round 5 (verify:accuracy-narrative, round 4, MAJOR undersell): in the palindrome-kept half cysteinylglycine
    # colocalizes with BOTH DPEP1 signals; the class rule's A2 reads only the stronger (by < 0.001). The registered
    # per-signal tables (l3fm__persignal.py; POST HOC, both directions declared together) give each signal's posterior,
    # top shared variant and position test, and the any-signal / every-signal counts. Every value is cross-checked here
    # against the grid, the annotation table and the registered class table (which the sensitivity must not move).
    pse, pss, psc = _rows(PSE), {(r["protein"], r["metabolite"]): r for r in _rows(PSS)}, _rows(PSC)
    ev = lambda P, M, h, X: [r for r in pse if (r["protein"], r["metabolite"], r["pal"], r["signal"]) == (P, M, h, X)]
    dA, dB = ev("DPEP1", "cysteinylglycine", "freq", "A"), ev("DPEP1", "cysteinylglycine", "freq", "B")
    _need(len(dA) == len(dB) == 1, "the per-signal evidence does not hold one row per DPEP1 kept-half signal")
    dA, dB = dA[0], dB[0]
    hA, hB = float(dA["cond_H4"]), float(dB["cond_H4"])
    _need(min(hA, hB) > 0.999 and dA["read_by_primary_rule"] == "True" and dB["read_by_primary_rule"] == "False"
          and 0 < hA - hB < 0.001, "the DPEP1 kept half no longer colocalizes with both signals above 0.999, the rule no "
          "longer reads the lead signal, or the margin is no longer below 0.001")
    _need(dA["protein_signal_rsid"] == sa_rs and dB["protein_signal_rsid"] == b_freq and dB["lead_rsid"] == b_freq
          and dA["lead_rsid"] == lead_f["rsid_ref"], "the per-signal evidence's DPEP1 signals / leads are not the grid's "
          "(signal A %s, lead %s; signal B %s)" % (sa_rs, lead_f["rsid_ref"], b_freq))
    _need((dA["A2_ens116pc"], dB["A2_ens116pc"]) == ("False", "True") and dA["lead_genes_ens116pc"] == "CHMP1A",
          "the position test no longer fails at DPEP1's lead signal and passes at the missense signal under Ensembl 116")
    dps, aps, ups = pss[("DPEP1", "cysteinylglycine")], pss[("ARG1", "arginine")], pss[("UMOD", "eGFR (Wuttke 2019 EA)")]
    for k, r in pss.items():
        c = cs.get(k)
        _need(c is not None and all(r["primary_n_gene_attributed_%s" % a] == c[pre + "n_gene_attributed"] for a, pre in
                                    (("ens116pc", ""), ("glist", "glist_"), ("ens116all", "ens116all_"))) and r["tier"] == c["tier"],
              "the per-signal table's primary counts / tier for %s are not the registered class table's" % (k,))
    i_ = lambda r, c: int(r[c])
    _need(i_(dps, "any_n_gene_attributed_ens116pc_freq") > 0 and i_(dps, "any_n_gene_attributed_ens116pc_drop") == 0
          and i_(dps, "any_n_gene_attributed_ens116pc") == i_(dps, "any_n_gene_attributed_ens116pc_freq")
          and i_(dps, "all_n_gene_attributed_ens116pc") == 0 and i_(aps, "all_n_gene_attributed_ens116pc") == 0
          and i_(aps, "any_n_gene_attributed_ens116pc") == len(gd),
          "the per-signal counts no longer read: DPEP1 kept half attributed under any, dropped half not, none under every; "
          "ARG1 none under every")
    # the cells the any reading attributes: exactly the kept-half cells that read each criterion once on the full panel
    # and do not use SuSiE as fitted
    got = sorted(r["variant"] for r in psc if (r["protein"], r["metabolite"], r["annotation"]) == ("DPEP1", "cysteinylglycine",
                                                                                                   "ens116pc")
                 and r["class_any"] == "gene-attributed bridge")
    want = sorted(r["variant"] for r in psc if (r["protein"], r["metabolite"], r["annotation"]) == ("DPEP1", "cysteinylglycine",
                                                                                                    "ens116pc")
                  and r["pal"] == "freq" and r["stab"] == "full" and r["susie"] != "asnow")
    _need(got == want and len(got) == i_(dps, "any_n_gene_attributed_ens116pc"),
          "the any-signal attributions of DPEP1-cysteinylglycine are no longer exactly the kept-half, full-panel, "
          "not-SuSiE-as-fitted cells: %s" % got[:3])
    aA = ev("ARG1", "arginine", "freq", "A") + ev("ARG1", "arginine", "drop", "A")
    _need(len(aA) == 2 and all(r["read_by_primary_rule"] == "False" and float(r["cond_H4"]) >= 0.8 and r["A2_ens116pc"] == "False"
                               and r["A2_glist"] == "False" and r["lead_genes_ens116pc"] == "MED23" for r in aA),
          "ARG1's lead signal no longer colocalizes with arginine with its top shared variant in MED23 (outside ARG1)")
    uA = [r for r in pse if (r["protein"], r["metabolite"]) == ("UMOD", "eGFR (Wuttke 2019 EA)")]
    _need(len(uA) == 4 and all(float(r["cond_H4"]) >= 0.8 for r in uA) and i_(ups, "primary_n_gene_attributed_ens116pc") == 0
          and i_(ups, "any_n_gene_attributed_ens116pc") == len(gd), "the UMOD-eGFR anchor no longer has two colocalizing signals, "
          "or is attributed under the registered rule, or is not attributed in all 48 under the any-signal reading")
    mx = max(max(i_(r, "%s_n_gene_attributed_%s" % (m, a)) for m in ("primary", "any") for a in ("ens116pc", "glist", "ens116all"))
             for k, r in pss.items() if k == ("DPEP1", "cysteinylglycine"))
    # 'the adjacent CHMP1A' (Discussion): no protein-coding gene lies between DPEP1's end and CHMP1A's start (Ensembl 116)
    _need(not [s_ for s_, (a_, e_, _) in gpos.items() if a_ > d_end and e_ < c_start] and c_start > d_end,
          "a protein-coding gene now lies between DPEP1 and CHMP1A (CHMP1A is not DPEP1's adjacent gene)")
    put("dpep1_ps", dict(h4_lo="0.999", margin="0.001", n_any=i_(dps, "any_n_gene_attributed_ens116pc_freq"),
                         max_any_primary=mx, arg1_a_gene=aA[0]["lead_genes_ens116pc"],
                         umod_any=i_(ups, "any_n_gene_attributed_ens116pc")),
        "l3_metabolite/bridge_persignal_{evidence,sensitivity,cells}.csv (l3fm__persignal.py), cross-checked against the grid, "
        "the annotation table and the registered class table")
    # FIX round 6 (round-5 verify:reviewer, 4 MAJOR framing defects): the every-signal reading is stated ONCE, in Methods
    # B5, with its calibration -- it attributes no pair at all under full gene spans and only DPEP1 pairs under glist-hg19,
    # and misses the UMOD-eGFR anchor for the same reason as the registered rule (its more strongly colocalizing signal lies
    # in PDILT) -- and the invariance that answers the forking-path worry (ARG1-arginine the only bridge the registered rule
    # and the any-signal reading attribute in all 48 variants under both named annotations). F3a keeps no per-signal
    # COUNT (the any-signal 8 would print one direction in Results); it prints the missense signal's LD-draw stability
    # against the 180 the stability variants require, beside the ARG1 signal the rule reads, and the DPEP1 bound over every
    # rule, reading and annotation.
    bmin = int(re.search(r"^BOOT_MIN = (\d+)", io.open(FM_ROB, encoding="utf-8").read(), re.M).group(1))
    ann3 = ("ens116pc", "glist", "ens116all")
    dmax = max(i_(dps, "%s_n_gene_attributed_%s" % (m, a)) for m in ("primary", "any", "all") for a in ann3)
    _need(3 * dmax <= len(gd) and dmax == max(i_(dps, "primary_n_gene_attributed_%s" % a) for a in ann3),
          "DPEP1-cysteinylglycine's largest count over every rule, reading and annotation is no longer the registered "
          "glist count, or exceeds a third of the grid")
    b_boot = int(dB["boot_n_qualifies_and_A2_ens116pc"])
    aR = [r for h_ in ("freq", "drop") for r in ev("ARG1", "arginine", h_, "A") + ev("ARG1", "arginine", h_, "B")
          if r["read_by_primary_rule"] == "True"]
    a_boot = {int(r["boot_n_qualifies_and_A2_ens116pc"]) for r in aR}
    _need(dB["A2_ens116pc"] == "True" and 0 < b_boot < bmin and len(aR) == 2 and len(a_boot) == 1
          and list(a_boot)[0] >= bmin and all(r["A2_ens116pc"] == "True" and r["signal"] != "A" for r in aR),
          "the missense signal no longer falls short of BOOT_MIN in the LD draws while ARG1's read signal clears it")
    uS = {h_: {x: ev("UMOD", "eGFR (Wuttke 2019 EA)", h_, x)[0] for x in "AB"} for h_ in ("drop", "freq")}
    _need(all(u["B"]["read_by_primary_rule"] == "True" and float(u["B"]["cond_H4"]) > float(u["A"]["cond_H4"]) >= 0.8
              and u["B"]["protein_signal_genes_ens116pc"] == u["B"]["lead_genes_ens116pc"] == "PDILT"
              and u["B"]["A2_ens116pc"] == u["B"]["A2_glist"] == "False" and u["A"]["A2_ens116pc"] == "True"
              and u["A"]["A2_glist"] == "False" for u in uS.values()),
          "UMOD-eGFR's more strongly colocalizing signal no longer lies in PDILT (read by the rule), or its lead signal no "
          "longer meets the position condition under full gene spans only")
    _need(i_(ups, "any_n_gene_attributed_ens116pc") == len(gd) and all(
        i_(ups, "%s_n_gene_attributed_%s" % (m, a)) == 0 for m, a in
        (("any", "glist"), ("primary", "ens116pc"), ("primary", "glist"), ("all", "ens116pc"), ("all", "glist"))),
          "the UMOD-eGFR anchor is no longer any-signal 48 under full spans and 0 under glist, registered / every 0 under both")
    bridges10 = [k for k, r in cs.items() if r["tier"] in ("tier 1", "tier 2")]
    full48 = [k for k in bridges10 if all(i_(pss[k], "%s_n_gene_attributed_%s" % (m, a)) == len(gd)
                                          for m in ("primary", "any") for a in ann3)]
    _need(len(bridges10) == F["bridges"] and full48 == [("ARG1", "arginine")] and all(
        i_(aps, "all_n_gene_attributed_%s" % a) == 0 for a in ann3) and len(aA) == 2
          and all(r["lead_genes_ens116pc"] == "MED23" and r["protein_signal_genes_ens116pc"] != "MED23" for r in aA),
          "ARG1-arginine is no longer the only bridge the registered rule and the any-signal reading attribute in all 48 "
          "variants under every annotation, or the every-signal reading attributes it, or its lead signal's top shared "
          "variant is no longer in MED23")
    ev_all = {a: [k for k, r in pss.items() if i_(r, "all_n_gene_attributed_%s" % a) > 0] for a in ann3}
    _need(ev_all["ens116pc"] == [] and ev_all["ens116all"] == [] and ev_all["glist"]
          and all(k[0] == "DPEP1" for k in ev_all["glist"]),
          "the every-signal reading now attributes a pair under full gene spans, or a non-DPEP1 pair under glist-hg19")
    _need(cs[("DPEP1", "cysteinylglycine")]["gcst"] in {r["metab_gcst"] for r in bsr},
          "DPEP1-cysteinylglycine's metabolite accession is not one of the screen's (CLSA) metabolite accessions")
    put("ps6", dict(dmax=dmax, b_boot=b_boot, a_boot=list(a_boot)[0], bmin=bmin),
        "bridge_persignal_{evidence,sensitivity}.csv (l3fm__persignal.py) + l3fm__robustness.BOOT_MIN; the UMOD, ARG1 and "
        "every-signal statements asserted above")
    # the tier rule, read from the class summary's DECLARATION (so prose and data carry the same words)
    dsrc = re.sub(r"\s+", " ", ast.get_docstring(ast.parse(io.open(CSB, encoding="utf-8").read()), clean=False))
    mt = re.search(r"enzyme-substrate \(protein, metabolite\) is in ENZYME_SUBSTRATE below: the protein is (an enzyme with [^,]+),", dsrc)
    _need(mt and mt.group(1) == "an enzyme with documented activity on the metabolite",
          "the class-summary declaration's tier wording moved: %r" % (mt.group(1) if mt else None))
    put("tier_phrase", mt.group(1), "l3fm__class_summary.py declaration, TIER 'enzyme-substrate' entry")
    # BRAP / ATXN2 / SH2B3; SMPD3 -> SLC7A6; GSTA3 -> GSTA2
    gb = [r for r in grid if r["protein"] == "BRAP" and cs[("BRAP", r["metabolite"])]["tier"] == "tier 2"]
    zb = {_f(r["A1_cis_max_abs_z"]) for r in gb}
    _need(len(zb) == 1 and list(zb)[0] < GW_Z, "BRAP's cis maximum is genome-wide now")
    put("brap_maxz", list(zb)[0], "bridge_class_robustness__ens116pc.csv A1_cis_max_abs_z")
    _need(all("SH2B3: E1 False E2 False" in r["rivals_full"] and "ATXN2: E1 False E2 False" in r["rivals_full"] for r in gb),
          "SH2B3 or ATXN2 now shares a BRAP bridge signal")
    atx = sorted(_f(r["d1_abf_H4"]) for r in bcs if r["protein"] == "ATXN2")
    _need(len(atx) == 3 and max(atx) < 0.8, "ATXN2 x the three 12q24 metabolites moved")
    put("atxn2_h4", atx, "bridge_class_summary.csv d1_abf_H4 (ATXN2)")
    rules = {(r["protein"], r["metabolite"]): r for r in _rows(RULES1)}
    sm = [k for k in t2 if k[0] == "SMPD3"]
    _need(len(sm) == 4 and all("shared-variant lead" in rules[k]["lead_variant_gene_context"] and
                               re.search(r"shared-variant lead \S+ in SLC7A6", rules[k]["lead_variant_gene_context"]) for k in sm),
          "an SMPD3 bridge's shared-variant lead is not inside SLC7A6")
    _need(all("SLC7A6" not in r["rivals_full"] for r in grid if r["protein"] == "SMPD3"), "SLC7A6 is a measured rival now")
    gs = rules[("GSTA3", "1,3,7-trimethylurate")]["lead_variant_gene_context"]
    _need("in GSTA2" in gs and "nearest GSTA2" in gs, "the GSTA3 bridge's leads no longer map to GSTA2")
    _need(all("GSTA2" not in r["rivals_full"] for r in grid if r["protein"] == "GSTA3"), "GSTA2 is a measured rival now")
    _need(cs[("GSTA3", "1,3,7-trimethylurate")]["modal_class"] == "shared variant, mediation uncertain", "GSTA3 class moved")
    # arm ancestry assertion (D1 controls, staging) and fine-mapping settings
    ctl = {r["id"]: r for r in _rows(D1CTL, "\t")}
    _need(all(ctl[c]["pass"] == "True" and "TRIPS" in ctl[c]["mutation_trips"] for c in ("C5", "C6", "C9")),
          "a D1 control failed or its mutation did not trip")
    m5 = re.search(r"median range ([\d.]+)-([\d.]+), shared (\d+)-(\d+)", ctl["C5"]["observed"])
    _need(m5, "C5's observed string changed shape")
    put("anc_af", dict(lo=m5.group(1), hi=m5.group(2), slo=int(m5.group(3)), shi=int(m5.group(4))), "d1 CONTROLS.tsv C5")
    put("anc_extracted", len(re.findall(r"window=True", ctl["C6"]["observed"])), "d1 CONTROLS.tsv C6")
    m9 = re.search(r"distinct (\d+)", ctl["C9"]["observed"])
    _need(m9 and int(m9.group(1)) == F["bsr_metab"], "C9 distinct accessions != screen metabolites")
    src = io.open(FM_COMMON, encoding="utf-8").read()
    m = re.search(r"^REF_N = (\d+)$", src, re.M)
    _need(m, "l3fm__common.REF_N not found")
    put("ld_ref_n", int(m.group(1)), "l3fm__common.py REF_N")
    _need("L = 10" in io.open(FM_SUSIE, encoding="utf-8").read(), "l3fm__susie.R no longer declares L = 10")
    put("susie_L", 10, "l3fm__susie.R settings")
    rob_src = io.open(FM_ROB, encoding="utf-8").read()
    _need("GRID: 2 x 3 x 2 x 2 x 2 = 48 cells" in rob_src and "s < 0.10 for BOTH fits" in rob_src and ">= 180 of 200" in rob_src,
          "the declared grid changed")
    axes = [Counter(r[a] for r in gd) for a in ("pal", "susie", "stab", "ratio", "rival")]
    _need([len(a) for a in axes] == [2, 3, 2, 2, 2], "the registered grid's axes are not 2 x 3 x 2 x 2 x 2")
    put("grid_cells", 48, "bridge_class_robustness__ens116pc.csv axes")
    # S26 (two-ancestry metabolite -> disease)
    q = _rows(S26Q)
    qc = Counter(r["quadrant"] for r in q)
    put("s26_n", len(q), "S26 quadrant source")
    put("s26_q", dict(qc), "S26 quadrant source")
    easassoc = [r for r in q if r["quadrant"] in ("concordant", "EAS-only", "discordant")]
    ev = [r for r in easassoc if r["coloc_evaluated"] == "True"]
    _need(all(r["coloc_evaluated"] != "True" for r in q if r not in easassoc), "a cell without an East Asian association was evaluated")
    col = [r for r in ev if r["colocalizes"] == "True"]
    put("s26_easassoc", len(easassoc), "S26 quadrant: concordant + EAS-only + discordant")
    put("s26_eval", len(ev), "S26 coloc_evaluated")
    put("s26_coloc", len(col), "S26 colocalizes")
    put("s26_coloc_q", dict(Counter(r["quadrant"] for r in col)), "S26")
    _need(all(r["coloc_evaluated"] != "True" for r in q if r["quadrant"] == "discordant"), "the discordant cell was evaluated")
    cg = [r for r in _rows(S26C) if r["tier"] == "1"]
    _need(len(cg) == F["s26_coloc"], "S26 colocgate tier-1 cells != colocalizing cells")
    hub = Counter(r["hub"] for r in cg)
    _need(set(hub) == {"CPS1", "FADS"} and hub["FADS"] == 1, "the colocalizing hubs moved: %s" % dict(hub))
    cps = [r for r in cg if r["hub"] == "CPS1"]
    byo = Counter(r["outcome"] for r in cps)
    _need(set(byo) == {"eGFR", "CAD"}, "CPS1 cells at other outcomes")
    eg = sorted(r["metabolite"] for r in cps if r["outcome"] == "eGFR")
    ca = sorted(r["metabolite"] for r in cps if r["outcome"] == "CAD")
    acyl = lambda xs: sum(1 for x in xs if x.endswith("ylglycine"))
    _need(sorted(x for x in eg if not x.endswith("ylglycine")) == ["creatine", "glycine"] and acyl(eg) == 3,
          "CPS1 eGFR cells are no longer glycine + creatine + three acylglycines: %s" % eg)
    _need(sorted(x for x in ca if not x.endswith("ylglycine")) == ["homoarginine"] and acyl(ca) == 2,
          "CPS1 CAD cells are no longer homoarginine + two acylglycines: %s" % ca)
    # FIX round 2 (verify:narrative, brief pre-audit): the CPS1 cells are gene-level only -- none is at type 2 diabetes and
    # none is a urea-cycle intermediate (the cycle's intermediates, carbamoyl phosphate through urea)
    UREA_CYCLE = {"carbamoyl phosphate", "ornithine", "citrulline", "argininosuccinate", "arginine", "urea"}
    _need("T2D" not in byo and not [x for x in eg + ca if x.lower() in UREA_CYCLE],
          "a CPS1 cell is at type 2 diabetes or is a urea-cycle intermediate: %s" % (eg + ca))
    put("s26_cps1", dict(n=len(cps), egfr=len(eg), cad=len(ca)), "S26 colocgate tier 1 hub CPS1")
    # FIX round: the triangle's arginine -> T2D leg is European (S26: EUR-only, heterogeneous)
    a26 = [r for r in q if (r["metabolite"], r["outcome"]) == ("arginine", "T2D")]
    _need(len(a26) == 1 and a26[0]["quadrant"] == "EUR-only" and _f(a26[0]["p_het"]) < 0.05,
          "arginine -> T2D is no longer a European-only, heterogeneous S26 cell")
    put("arg_t2d_s26", dict(hetp=_f(a26[0]["p_het"]), eur=_f(a26[0]["beta_eur"]), eas=_f(a26[0]["beta_eas"])),
        "S26 quadrant source arginine x T2D")
    fa = [r for r in cg if r["hub"] == "FADS"][0]
    _need((fa["metabolite"].startswith("1_linoleoyl_gpc"), fa["outcome"]) == (True, "T2D"), "the FADS cell moved")
    # Figure 9 funnels
    f9 = {(r["layer"], r["stage"]): r for r in _rows(F09)}
    put("f9", {k: int(v["n"]) for k, v in f9.items()}, "Figure_09 source")
    _need(f9[("L1", "3")]["exemplar"] == "PAM-T2D, PCSK9-CAD, RAD23B-T2D", "the convergent exemplars moved")
    _need((int(f9[("L2", "1")]["n"]), int(f9[("L2", "2")]["n"]), int(f9[("L2", "3")]["n"])) ==
          (F["cen_test"], len(sup), len(rob)), "Figure_09 L2 disagrees with the census")
    _need((int(f9[("L3", "1")]["n"]), int(f9[("L3", "2")]["n"]), int(f9[("L3", "3")]["n"])) ==
          (F["bsr_n"], F["bridges"], len(t1)), "Figure_09 L3 disagrees with the screen / class table")
    _need((int(f9[("L1", "2")]["n"]), int(f9[("L1", "3")]["n"])) == (F["n_pairs"], 3), "Figure_09 L1 disagrees with A11")
    # Table 1: the census outcome resources (the census builder's registry is the only place the counts are recorded)
    sys.path.insert(0, LB)
    import eas__census_common as ECC
    O = ECC.OUTCOMES
    tn = {}
    for k, v in O.items():
        tn[k] = v["ntot"] if v.get("ntot") else v["ncase"] + v["nctrl"]
    put("t1_n", tn, "eas__census_common.OUTCOMES (ncase + nctrl, or ntot)")
    put("t1_acc", {k: v["accession"] for k, v in O.items()}, "eas__census_common.OUTCOMES accession")
    _need(tn["eGFR"] == F["egfr_eas_n"], "the eGFR N in the registry and in S23N disagree")
    _need((O["CAD"]["ncase"], O["CAD"]["nctrl"]) == (32512, 146214), "the BBJ CAD counts moved (they identify its source)")

    _CACHE["F"], _CACHE["SRC"] = F, SRC
    return F, SRC


# ------------------------------------------------------------------------------------------------ references (stage-1)
# Every reference the new text cites, keyed; existing entries by a title fragment of the STAGE-1 list (or of
# renumber_refs.NEW_REFS), new entries by their own NEW_REFS text. Numbers are STAGE-1 numbers: renumber_refs.py puts
# the whole list in first-citation order after the blocks are applied (gate P18).
EXISTING = {
    "UKBPPP": "Plasma proteomic associations with genetics and health in the UK Biobank",
    "SUSIE": "A simple new approach to variable selection in regression",
    "WUTTKE": "A catalog of genetic loci associated with kidney function",
    "WALLACE": "A more accurate method for colocalisation analysis",
    "CKB": "Comparative studies of 2168 plasma proteins",
    "CLSA": "Genomic atlas of the plasma metabolome",
    "CHENG": "Genetic mapping of serum metabolome to chronic diseases among Han Chinese",
    "CALDWELL": "Arginase: a multifaceted enzyme important in",
    "ADACHI1990": "Primary structure of human microsomal dipeptidase",
    "KULIESIUS": "Epitope effect prevalence in affinity-based pQTL studies",
}
NEW_KEYS = ["SPRACKLEN", "MISHRA", "CHEN2024", "KANAI2021", "MIYAZAWA", "ENZAN", "SAKAUE", "SIN", "MCINTYRE",
            "ADACHI1989", "HABIB", "HARTIALA", "WYSS", "MUSUNURU", "NEWTONCHEH", "MALLAT", "GRAINGER", "IL6RGC",
            "KUROSU", "KIR", "HUANG", "GE", "SHEMYAKIN", "PIATTI", "CERVENKA", "BROER", "TASSANEEYAKUL",
            "KOZAK", "ULAGANATHAN"]


def refs():
    """key -> stage-1 reference number."""
    if "R" in _CACHE:
        return _CACHE["R"]
    import docx
    import renumber_refs as RR
    d = docx.Document(STAGE1)
    ps = d.paragraphs
    h = [i for i, p in enumerate(ps) if p.text.strip() == "References"]
    _need(len(h) == 1, "stage-1: one References heading")
    ents = [re.sub(r"\s+", " ", p.text).strip() for p in ps[h[0] + 1:] if p.text.strip()]
    _need(all(re.match(r"^%d\. " % (i + 1), e) for i, e in enumerate(ents)), "stage-1 entries not numbered 1..N")
    allt = ents + ["%d. %s" % (len(ents) + i + 1, t) for i, t in enumerate(RR.NEW_REFS)]
    R = {}
    for k, frag in EXISTING.items():
        hit = [i + 1 for i, e in enumerate(allt) if frag in e]
        _need(len(hit) == 1, "reference %s: fragment %r matches %d entries" % (k, frag, len(hit)))
        R[k] = hit[0]
    _need(len(RR.NEW_REF_KEYS) == len(RR.NEW_REFS), "renumber_refs.NEW_REF_KEYS is not aligned with NEW_REFS")
    for k in NEW_KEYS:
        idx = [j for j, kk in enumerate(RR.NEW_REF_KEYS) if kk == k]
        _need(len(idx) == 1, "new reference %s is not declared once in renumber_refs.NEW_REF_KEYS" % k)
        R[k] = len(ents) + idx[0] + 1
    for k, frag in EXISTING.items():          # a NEW_REFS entry keyed in both places must agree
        if k in RR.NEW_REF_KEYS:
            _need(R[k] == len(ents) + RR.NEW_REF_KEYS.index(k) + 1, "reference %s: fragment and key disagree" % k)
    _CACHE["R"] = R
    return R


def C(*keys):
    """a bracketed stage-1 citation; renumber_refs sorts the numbers after renumbering."""
    R = refs()
    return "[%s]" % ", ".join(str(R[k]) for k in keys)


# ------------------------------------------------------------------------------------------------------------- texts
def texts():
    """block id -> text (INSERTS.md fenced blocks), plus 'edits' for the inherited paragraphs (rebase_v4_edits)."""
    F, _ = facts()
    n, g = F["n_pairs"], F["n_genes"]
    ce, tg, il, fg, um = F["celsr2"], F["tgfb1"], F["il6r"], F["fgfr4"], F["umod"]
    sc = F["fgfr4_scan"]
    _need(sc["cand"] < 10 and sc["strict"] < 10 and sc["dflt"] < 10, "a scan candidate count reached ten: spell-out rule")
    pal_pairs = F["cen_pal_pairs"]
    _need(pair("FGF5", "eGFR") in pal_pairs and pair("GALNT3", "T2D") in pal_pairs, "the named palindromic pairs moved")
    pr = F["cen_proxy"]
    lower = lambda x: x.replace(u"–Stroke", u"–stroke")
    ties = [lower(x) for x in F["cen_ties_pz"]]
    T = {}
    T["A1"] = "Native East Asian Colocalization Census"
    T["A2"] = (
        u"To test whether network edges recur when both sides of the association are measured in East Asians, we "
        u"re-estimated every eligible network pair natively, pairing an East Asian cis-pQTL instrument with an East Asian "
        u"outcome GWAS instead of transporting the European instrument. A pair was eligible when its protein carried a "
        u"genome-wide-significant cis-pQTL (p < 5 × 10⁻⁸ within ±{win} kb of the gene) in the China "
        u"Kadoorie Biobank (CKB) Olink pQTL resource {ckb}, retrieved through its open per-region interface; {cn} of the "
        u"{n} pairs, over {cp} of the {g} proteins, qualified, and {ct} were testable, because no public East Asian NAFLD "
        u"GWAS was available for APOE–NAFLD. The instrument was the cis-pQTL with the smallest p-value that lay within "
        u"±{win} kb of both the gene and the protein's lead variant, had minor-allele frequency ≥ 0.01 and was "
        u"present and harmonisable in the outcome GWAS, with ties broken by the larger |z| and then by position. It was "
        u"the lead variant for {nl} pairs; for {prox} the lead was absent from the outcome GWAS and the next eligible "
        u"variant was used. For {ties}, variants tied on both p and |z|, so position decided the instrument. The weakest "
        u"instrument had F = {fmin}. The matched East Asian outcomes were the Asian Genetic Epidemiology Network (AGEN) "
        u"type 2 diabetes meta-analysis ({a_t2d}) {spr}, the East Asian GIGASTROKE any-stroke analysis ({a_st}) {mis} and "
        u"the Taiwan–Japan eGFR meta-analysis (N = {negfr}; risk-aligned as eGFR*) {chen}. Biobank Japan provided "
        u"coronary artery disease {kan}, atrial fibrillation {miy}, heart failure ({a_hf}) {enz} and chronic kidney "
        u"disease ({a_ckd}) {sak}. Exposure and outcome were joined on position and allele pair, with GRCh38 positions "
        u"lifted to GRCh37 for the {lift} pairs whose outcome file is on that build."
    ).format(win=F["cen_window_kb"], ckb=C("CKB"), cn=F["cen_n"], n=n, cp=F["cen_prot"], g=g, ct=F["cen_test"],
             nl=F["cen_lead"], prox=oxford([lower(x) for x in pr]), ties=oxford(ties), fmin="%.1f" % F["cen_Fmin"],
             a_t2d=F["cen_acc"]["T2D"], spr=C("SPRACKLEN"), a_st=F["cen_acc"]["Stroke"], mis=C("MISHRA"),
             negfr=thou(F["egfr_eas_n"]), chen=C("CHEN2024"), kan=C("KANAI2021"), miy=C("MIYAZAWA"),
             a_hf=F["cen_acc"]["HF"], enz=C("ENZAN"), a_ckd=F["cen_acc"]["CKD"], sak=C("SAKAUE"), lift=F["cen_lifted"])
    T["A3"] = (
        u"Colocalization used the backbone's coloc.abf priors and dual rule (PP.H4 ≥ 0.80 and PP.H4/PP.H3 ≥ 5; "
        u"p1 = p2 = 1 × 10⁻⁴, p12 = 1 × 10⁻⁵) over ±{win} kb around the East Asian lead, "
        u"with p12 = 1 × 10⁻⁶ as the conservative stability check. Wald ratios were compared between "
        u"ancestries by a heterogeneity z-test; the East Asian eGFR GWAS is on a different scale, so eGFR pairs were "
        u"compared by direction only. The instruments for {npal} pairs ({nvar} variants), including FGF5–eGFR and "
        u"GALNT3–T2D, are palindromic; their colocalization posteriors do not depend on allele orientation, but the "
        u"sign of their Wald estimates assumes forward-strand reporting. Three sensitivity analyses tested the East Asian "
        u"data themselves. Because the AGEN type 2 diabetes ({agen} pairs) and GIGASTROKE ({giga} pairs) resources include "
        u"CKB, the exposure cohort, those pairs were repeated against CKB-free Biobank Japan type 2 diabetes and ischaemic "
        u"stroke GWAS {sak}. All {ct} pairs were re-run in R coloc 5.2.3 with the European analysis's effect-size priors "
        u"(case–control outcomes by case fraction; eGFR with a declared standard deviation) in place of the primary "
        u"common prior. Each of the {wreg} CKB regions was checked for rank agreement between |z| and the reported p-value "
        u"(Spearman ρ ≥ 0.95) and for allele-frequency agreement with 1000 Genomes East Asians. Unlike the outcome-side "
        u"and pQTL-portability analyses, the census estimates the full East Asian protein-to-disease association; like "
        u"them, it is reported independently of the European evidence tiers (Supplementary Table S13)."
    ).format(win=F["cen_window_kb"], npal=WORDS[len(pal_pairs)], nvar=WORDS[F["cen_pal_variants"]], agen=F["ckb_ov_agen"],
             giga=F["ckb_ov_giga"], sak=C("SAKAUE"), ct=F["cen_test"], wreg=F["w_regions"])
    # FIX round 3 (verify:narrative): A3 ran to 409 words by Word; the off-network pilot and the scan FGFR4-T2D came from
    # are their own paragraph. The heart-failure screen's one candidate (CDKN1A, shared) is named, and the reader is told
    # why five pairs that failed the East-Asian-specific filter entered colocalization (the records: for transparency).
    _need(sc["cand"] == sc["strict"] + F["scan_egfr_carried"], "the scan's candidates are no longer the strict ones plus "
          "the carried eGFR hits")
    T["A3b"] = (
        u"{pn} off-network pairs over {pp} proteins from an earlier targeted pilot, among them FGFR4–T2D and UMOD–eGFR, "
        u"are reported unchanged as a separate block, under the pilot's instrument rule (the smallest-p cis-pQTL present "
        u"in the outcome GWAS, ties broken by position) and outcome resources. FGFR4–T2D entered that pilot from an "
        u"earlier East Asian scan of CKB-instrumented proteins: {sc_t} protein–outcome tests over {sc_p} proteins against "
        u"type 2 diabetes (AGEN; {sc_t2d} proteins), coronary artery disease ({sc_cad}), stroke ({sc_st}) and eGFR "
        u"({sc_eg}). Heart failure and atrial fibrillation were screened in East Asians too, but their per-protein results "
        u"were not kept; the heart-failure screen's one candidate, {hf}, was shared with Europeans rather than "
        u"East-Asian-specific (no significant cross-ancestry heterogeneity), and neither screen contributed a pair to "
        u"colocalization. {sc_c} of the scan's protein–outcome pairs entered colocalization. {sc_s} were "
        u"East-Asian-specific by its filter: false-discovery-significant in East Asians, null in Europeans, heterogeneous "
        u"between the two and on a common instrument outside known pleiotropic loci. The other {sc_e}, eGFR associations "
        u"false-discovery-significant in East Asians and absent from the European eGFR hits, failed that filter and were "
        u"carried for transparency. "
        u"{sc_d} reached PP.H4 ≥ 0.80 at the default prior, and only FGFR4–T2D kept PP.H4 ≥ 0.80 at the conservative "
        u"prior."
    ).format(pn=WORDS[F["pil_n"]].capitalize(), pp=WORDS[F["pil_prot"]], sc_t=thou(sc["tests"]), sc_p=sc["prot"],
             sc_t2d=sc["arms"]["T2D"], sc_cad=sc["arms"]["CAD"], sc_st=sc["arms"]["Stroke"], sc_eg=sc["arms"]["eGFR"],
             hf=F["hf_candidate"], sc_c=WORDS[sc["cand"]].capitalize(), sc_s=WORDS[sc["strict"]].capitalize(),
             sc_e=WORDS[F["scan_egfr_carried"]], sc_d=WORDS[sc["dflt"]].capitalize())
    T["B2"] = (
        u"To test whether the plasma-protein network extends into the metabolome, we applied the colocalization framework "
        u"of the disease backbone to plasma metabolites in three complementary analyses. Every European disease outcome in "
        u"this layer is the protein arm's own UK-Biobank-free GWAS (Table 1), and the European metabolite QTLs are those of "
        u"the Canadian Longitudinal Study on Aging {clsa}. First, a metabolite node — a protein → metabolite "
        u"→ disease triangle — was admitted only when the protein was a colocalization-supported network protein "
        u"for the same disease that the metabolite arm implicated. The {nn} candidate (gene, outcome) node pairs, over {ng} "
        u"genes, carried from a prior metabolite-QTL analysis were re-evaluated against the {n}-edge backbone, and the "
        u"metabolite–disease edge of each admitted node was colocalized under the dual rule, with p12 = 5 × "
        u"10⁻⁶ and 1 × 10⁻⁶ as prior checks. Second, and independently of that filter, we screened "
        u"every network protein for a colocalized metabolite: each colocalization-supported protein's cis-pQTL was "
        u"colocalized with every plasma metabolite whose cis-acting metabolite quantitative trait locus (mQTL) lead lay "
        u"within ±500 kb of the protein-coding gene "
        u"(coloc.abf; dual rule PP.H4 ≥ 0.80 and PP.H4/PP.H3 ≥ 5; priors p1 = p2 = 1×10⁻⁴ and "
        u"p12 = 1×10⁻⁵), harmonised on GRCh38 by chromosomal position, and all {nb} resulting "
        u"protein–metabolite pairs, over {bp} proteins and {bm} metabolites, were tested. Third, as a check on the "
        u"colocalization gate itself, we compared European {clsa} and East Asian {cheng} metabolite → disease Mendelian "
        u"randomization across five CKM outcomes and retained a cross-ancestry cell only when the metabolite arm "
        u"colocalized with the outcome in East Asians under the same dual rule."
    ).format(clsa=C("CLSA"), nn=F["node_n"], ng=F["node_genes"], n=n, nb=F["bsr_n"], bp=F["bsr_prot"], bm=F["bsr_metab"],
             cheng=C("CHENG"))
    af = F["anc_af"]
    T["B3"] = (
        u"Every protein arm in the screen is the UKB-PPP European-ancestry discovery release {ukb}. Because the release uses "
        u"identical file names across its ancestry folders, each arm's ancestry was asserted rather than read from a file "
        u"name: the {nx} arms extracted from the release archives were taken from the European discovery folder and "
        u"required the discovery sample size, and every arm's allele frequencies matched those of a reference European "
        u"discovery arm (median absolute difference {lo}–{hi} over {slo}–{shi} shared variants). Every metabolite "
        u"accession likewise lay in the European-ancestry block of {clsa}. The five bridge loci were then fine-mapped to ask "
        u"whether each shared signal belongs to the screened protein's gene. SuSiE {susie} (L = {L}) was fitted to every "
        u"trait with a 1000 Genomes European LD reference ({refn} individuals), followed by multi-signal (coloc.susie "
        u"{wal}) and conditional two-signal colocalization. A pair was a gene-attributed bridge when five conditions held: "
        u"protein and metabolite shared a signal; the metabolite tracked the protein in the same ratio at the protein's two "
        u"independent signals, where it had two; the protein's cis signal was genome-wide significant; the protein's signal "
        u"variant and the shared-variant lead both lay in, or nearest to, the bridge gene (the position condition); and no "
        u"other measured protein at the locus shared the signal. A shared signal that was not proportional was classed as a shared variant with "
        u"mediation uncertain, one failing any of the last three conditions as a locus-level colocalization with the gene "
        u"uncertain, and the rest as not supported."
    ).format(ukb=C("UKBPPP"), nx=WORDS[F["anc_extracted"]], lo=af["lo"], hi=af["hi"], slo=thou(af["slo"]), shi=thou(af["shi"]),
             clsa=C("CLSA"), susie=C("SUSIE"), L=F["susie_L"], refn=F["ld_ref_n"], wal=C("WALLACE"))
    T["B4"] = (
        u"So that no single analytic choice could decide a class, every pair was classified under a grid of {cells} rule "
        u"variants, frozen before any class was computed under it. The variants cross five choices: palindromic variants "
        u"dropped or kept where allele frequency resolves their strand; SuSiE evidence used as fitted, used only when its "
        u"LD-mismatch statistic was below 0.10, or not used; each criterion read once on the full LD panel or required in "
        u"at least 180 of 200 LD-reference bootstrap draws; a raw or a standardised proportionality test; and a "
        u"neighbouring protein that shares the signal blocking attribution always, or only when it is genome-wide "
        u"significant at the bridge protein's lead signal. Gene bodies were the full spans of Ensembl release 116 "
        u"protein-coding genes (GRCh38), non-coding transcripts included, with glist-hg19 as a sensitivity annotation. "
        u"Each pair is reported as the number of the {cells} variants giving each class (Supplementary Table S14); the "
        u"class rules, the gene annotation and the tier rule "
        u"below were fixed after the first fine-mapping results had been seen. For presentation, the screen's bridges were "
        u"tiered by chemistry and direction. Bridge tier 1 — named so to keep it distinct from the drug-development "
        u"Tiers 1–4 and from the three evidence tiers — comprises bridges in which the protein is {tier} "
        u"(ARG1 on arginine {c_arg}; DPEP1 on cysteinylglycine {c_mci} and on its S-substituted forms {c_dp}) and the "
        u"shared signal is substrate-consistent, more enzyme with less metabolite; bridge tier 2 comprises every other "
        u"bridge. The bridge tier grades biochemical interpretability, leaves the count of bridges unchanged and is "
        u"reported beside gene attribution rather than implying it."
    ).format(cells=F["grid_cells"], tier=F["tier_phrase"], c_arg=C("SIN", "CALDWELL"), c_mci=C("MCINTYRE"),
             c_dp=C("ADACHI1989", "HABIB"))
    # FIX round 5: the POST HOC per-signal sensitivity, its own Methods paragraph (B4 is 276 words). Declared in both
    # directions in l3fm__persignal.py; the registered classes are not changed. The UMOD-eGFR renal anchor (Results) is the
    # known two-signal case: every one of its four per-signal colocalizations reaches the rule's PP.H4, the registered rule
    # attributes it to UMOD in no variant and the any-signal reading in all 48 (Ensembl gene spans).
    # FIX round 6 (round-5 verify:reviewer): the every-signal reading is stated here ONCE, with its calibration (it
    # attributes no pair under full gene spans; UMOD-eGFR is missed by it and by the registered rule for one stated reason,
    # its more strongly colocalizing signal lying in PDILT), UMOD is reported under BOTH named annotations, and one
    # invariance sentence answers the forking-path worry. The provenance clause is exact (l3fm__persignal.py: 'written
    # AFTER the outcome was known'), REPLACING 'added after the grid results had been seen'.
    # FIX round 6b (round-6 verify:reviewer P1, verify:numbers): the ORDER. Round 6 stated ARG1's every-signal miss before
    # the calibration, so the 'because ... MED23' read as an excuse for the reading rather than as a property of it, and the
    # glist clause ('only DPEP1 pairs') stood near ARG1's 0 of 48, where it can be read as reversing the ARG1/DPEP1 ranking.
    # The order is now: definition premise + provenance -> the two readings -> the grid and S14 -> UMOD under both
    # annotations -> the every-signal calibration, then its ARG1 miss with the MED23 reason -> ENDING on the invariance.
    # Three sentences over 35 words are also split (4 -> 1); no sentence is merged (233 -> 229 words by Word).
    T["B5"] = (
        u"Because the class rules place a shared signal in a gene using only the more strongly colocalizing of a protein's "
        u"two signals, a per-signal sensitivity was specified after its results for the DPEP1, ARG1 and UMOD pairs were "
        u"known. It applied that position condition to each signal "
        u"whose conditional colocalization reached PP.H4 ≥ 0.80, under the same stability requirement. A pair was "
        u"gene-attributed when any such signal met it or, in the opposite direction, only when every one did. The "
        u"{cells} variants, the other conditions and the registered classes were unchanged, and both readings are reported "
        u"for every fine-mapped pair (Supplementary Table S14). The UMOD–eGFR renal anchor (Results) has two conditionally "
        u"colocalizing signals, the more strongly colocalizing of which lies in PDILT. The any-signal reading attributed "
        u"it to UMOD in all {cells} variants under full gene spans and in none under glist-hg19. The registered rule and "
        u"the every-signal reading, which both need that PDILT signal to meet the position condition, attributed it in none "
        u"under either annotation. The every-signal reading attributed no pair at all under full gene spans, and only "
        u"DPEP1 pairs under glist-hg19. It missed ARG1–arginine because the top variant that ARG1's lead signal shares "
        u"with arginine lies in {a1gene}. Under the registered rule and the any-signal reading, ARG1–arginine was the only "
        u"bridge attributed to its gene in all {cells} variants under both annotations."
    ).format(cells=F["grid_cells"], a1gene=F["dpep1_ps"]["arg1_a_gene"])
    _need(F["dpep1_ps"]["umod_any"] == F["grid_cells"], "the any-signal reading no longer attributes UMOD-eGFR in every variant")
    sup = F["cen_coloc"]
    T["E1"] = (
        u"Measured natively, with an East Asian instrument against an East Asian outcome, the network's three most secure "
        u"edges in East Asians were a lipid-locus edge that transferred almost unchanged, a renal edge and a cytokine edge "
        u"whose direction reversed (Figure 3b; Supplementary Table S13). Of the {cn} network pairs whose protein carries a "
        u"genome-wide-significant East Asian cis-pQTL, {ct} were testable and {ns} colocalized: {sup}. {nr} held at the "
        u"conservative prior p12 = 1 × 10⁻⁶. CELSR2–CAD, the edge from the 1p13 LDL-cholesterol locus "
        u"(Methods), colocalized at PP.H4 = {ce_h4} ({ce_h46} at the conservative prior), with an East Asian estimate almost "
        u"identical to the European one ({ce_eas} and {ce_eur} log-odds per standard deviation (SD) of protein; "
        u"heterogeneity p = {ce_p}). "
        u"FGF5–eGFR gave the most decisive colocalization in the census, PP.H4 = {fg_h4} under both priors. "
        u"TGFB1–CAD also held ({tg_h4}; {tg_h46}) but with the opposite sign: genetically higher TGF-β1 "
        u"accompanied higher coronary risk in Europeans ({tg_eur} log-odds per SD) and lower risk in East Asians ({tg_eas}; "
        u"heterogeneity z = {tg_z}), on a common, non-palindromic East Asian instrument (minor-allele frequency {tg_maf}), a "
        u"different variant {tg_kb} kb from the European one; all {tg_ns} genome-wide-significant East Asian cis variants "
        u"that could serve as its instrument (those present in the outcome GWAS) give a negative estimate, so the sign "
        u"does not depend on which of them is chosen. Unlike those of CELSR2–CAD and FGF5–eGFR, however, the European "
        u"colocalization of TGFB1–CAD does not hold at the conservative prior (PP.H4 = {tg_eur6}), so TGFB1–CAD is robust "
        u"on the East Asian side only. The other {nfell} colocalizations did not hold at the conservative prior. The IL6R–CAD positive control colocalized "
        u"natively (PP.H4 = {il_h4}) with the European protective direction but fell to {il_h46}, as its European "
        u"colocalization falls to {il_eur6}; ABO–T2D kept PP.H4 = {ab_h4} but missed the ratio criterion "
        u"(PP.H4/PP.H3 = {ab_r}); and {fell} fell below the rule. Apart from TGFB1–CAD, every colocalizing pair with a "
        u"strand-resolvable East Asian instrument ({sr}) had the European direction of effect."
    ).format(cn=F["cen_n"], ct=F["cen_test"], ns=len(sup), sup=oxford(sup), nr=WORDS[len(F["cen_robust"])].capitalize(),
             ce_h4="%.3f" % ce["h4"], ce_h46="%.3f" % ce["h4_6"], ce_eas=neg(ce["eas"], 3), ce_eur=neg(ce["eur"], 3),
             ce_p="%.2f" % ce["hetp"], fg_h4="%.3f" % F["fgf5"]["h4"], tg_h4="%.3f" % tg["h4"], tg_h46="%.3f" % tg["h4_6"],
             tg_eur=sgn(tg["eur"], 3), tg_eas=neg(tg["eas"], 3), tg_z="%.2f" % tg["hetz"], tg_maf="%.2f" % tg["maf"],
             tg_kb="%d" % round(tg["dist_kb"]), tg_ns=F["tgfb1_strong"], tg_eur6="%.3f" % F["eur_1e6"][("TGFB1", "CAD")],
             nfell=WORDS[len(sup) - len(F["cen_robust"])], il_h4="%.3f" % il["h4"], il_h46="%.3f" % il["h4_6"],
             il_eur6="%.3f" % il["eur_h4_6"], ab_h4="%.3f" % F["abo_t2d"]["h4_6"], ab_r="%.2f" % F["abo_t2d"]["ratio_6"],
             fell=oxford(F["cen_fell_below"]), sr=oxford(F["cen_sr_same"]))
    _need(F["tgfb1_strong_neg"] == F["tgfb1_strong"], "not every strong TGFB1 variant gives a negative estimate")
    ap = F["apoe_ckbfree"]
    w1 = F["w1_fail"]
    T["E2"] = (
        u"Most of the {nc} pairs that did not colocalize natively lacked an East Asian outcome association rather than "
        u"contradicting the European edge. In {h1} the East Asian outcome showed no association at the locus (PP.H1, the "
        u"posterior for an exposure-only signal, dominant); in {nh3} the East Asian data resolved distinct causal variants "
        u"({h3}); and DNER–eGFR reached PP.H4 = {dner}, below the rule. The distinct-variant calls include the "
        u"PCSK9–CAD positive control, whose native estimate nonetheless had the European direction (heterogeneity p = "
        u"{pc_p}), so native non-colocalization cannot by itself be read as non-transfer. They also include ARG1–T2D: "
        u"its European instrument is compatible in East Asian outcomes (Figure 3a), but its native instrument is a "
        u"low-frequency variant (minor-allele frequency {ar_maf}) whose estimate differs from the European one "
        u"(heterogeneity p = {ar_p}). The calls were robust to the effect-size prior: replacing the primary common prior "
        u"with the European analysis's effect-size priors changed none of the {ct}. Replacing the two outcome resources that include CKB, the exposure cohort, changed one "
        u"colocalization call: APOE–T2D, not supported on AGEN (PP.H4 = {ap_a}), colocalized on CKB-free Biobank Japan "
        u"type 2 diabetes ({ap_b}; {ap_b6} at the conservative prior), giving {fs} colocalizing and {fr} robust pairs, while "
        u"ARG1–T2D moved from distinct variants to no association. Allele frequencies agreed with 1000 Genomes East "
        u"Asians in all {wreg} CKB regions, whereas rank agreement between |z| and p fell below 0.95 in two ({w1a}, "
        u"ρ = {w1av}; {w1b}, ρ = {w1bv}), whose calls are reported unadjusted. {nx} census pairs also appear in "
        u"the outcome-side analysis (Figure 3a), which uses the European instrument and, for some endpoints, a different "
        u"East Asian GWAS; for {nsd} of them the East Asian estimates differ in sign between the two analyses."
    ).format(nc=F["cen_nonc"], h1=F["cen_h1"], nh3=len(F["cen_h3_pairs"]), h3=oxford(F["cen_h3_pairs"]),
             dner="%.3f" % F["dner_h4"], pc_p="%.2f" % F["pcsk9_hetp"], ar_maf="%.3f" % F["arg1_nat"]["maf"],
             ar_p="%.3f" % F["arg1_nat"]["hetp"], ct=F["cen_test"], ap_a="%.3f" % ap["agen"], ap_b="%.3f" % ap["bbj"],
             ap_b6="%.3f" % ap["bbj_6"], fs=F["ckbfree_totals"][0], fr=F["ckbfree_totals"][1], wreg=F["w_regions"],
             w1a=w1[0][0], w1av="%.3f" % w1[0][1], w1b=w1[1][0], w1bv="%.3f" % w1[1][1],
             nx=_num_word_cap(F["cen_x_outcomeside"]), nsd=len(F["cen_x_signdiff"]))
    _need(len(w1) == 2 and F["w2_pass"] == F["w_regions"] and F["rsens_changed"] == 0, "the census channel/prior checks moved")
    _need(F["cen_nonc"] - F["cen_h1"] - len(F["cen_h3_pairs"]) == 1 and F["cen_h1"] * 2 > F["cen_nonc"],
          "'most ... lacked an East Asian outcome association' no longer holds")
    T["E3"] = (
        u"Two proteins outside the European network, carried from an earlier pilot (Methods), also colocalized in East "
        u"Asians, both at the conservative prior. FGFR4–T2D colocalized at PP.H4 = {fg_h4} ({fg_h46} at the conservative "
        u"prior), with genetically higher FGFR4 associated with higher type 2 diabetes risk ({fg_b} log-odds per SD; "
        u"p = {fg_p}). Its European estimate is null (p = {fg_eurp}; European colocalization not tested) and differs from "
        u"the East Asian one (heterogeneity p = {fg_het}), the pattern of an East-Asian-specific edge. An earlier East "
        u"Asian scan of CKB-instrumented proteins selected it for that pattern, judged there against a different European "
        u"type 2 diabetes GWAS (Methods): of the {sc_c} candidates the scan carried into colocalization, it alone kept "
        u"PP.H4 ≥ 0.80 at the conservative prior, so its East Asian estimates are selection-stage values. Its "
        u"colocalization held on CKB-free Biobank Japan type 2 diabetes "
        u"(PP.H4 = {fg_free}), but Biobank Japan also contributes to AGEN, so that re-run is not independent replication. "
        u"UMOD–eGFR colocalized at PP.H4 = "
        u"{um_h4} ({um_h46}), at the renal anchor that European coloc.abf resolves into distinct variants and European "
        u"multi-signal SuSiE colocalizes (Supplementary Table S4). The other {rest} pilot pairs did not colocalize."
    ).format(fg_h4="%.3f" % fg["h4"], fg_h46="%.3f" % fg["h4_6"], fg_b=sgn(fg["eas"], 3), fg_p=sci(fg["easp"]),
             fg_eurp="%.2f" % fg["eurp"], fg_het=sci(fg["hetp"]), fg_free="%.3f" % fg["ckbfree"], um_h4="%.3f" % um["h4"],
             um_h46="%.3f" % um["h4_6"], rest=WORDS[F["pil_n"] - len(F["pil_coloc"])], sc_c=WORDS[sc["cand"]])
    _need(fg["eurp"] > 0.05 and fg["hetp"] < 0.05, "FGFR4-T2D: the European estimate is no longer null or no longer differs")
    _need(F["fgfr4_eur_src"]["printed"] == "FinnGen_R10_T2D" and F["fgfr4_eur_src"]["scan"] == "DIAMANTE",
          "'judged there against a different European type 2 diabetes GWAS' no longer holds")
    T["F1"] = "The plasma-protein network extends into the metabolome"
    at, a1t, ab, al = F["arg_t2d"], F["arg1_t2d"], F["arg1_bmi"], F["arg_locus"]
    T["F2"] = (
        u"At ARG1, the network's type 2 diabetes edge meets the enzyme's own substrate: in Europeans, ARG1 → arginine → "
        u"type 2 diabetes is the one protein–metabolite–disease triangle in which all three edges colocalize (Figure 6a, "
        u"c). Of the {nn} candidate (gene, outcome) node pairs, over {ng} genes, two had a network protein colocalized for "
        u"their own outcome, {adm}. APOE did not complete a triangle, because both of its candidate metabolites, a "
        u"plasmalogen phosphatidylethanolamine and a sphingomyelin, resolved into variants distinct from the type 2 "
        u"diabetes signal (PP.H3 = {apo_h3}). ARG1 colocalized with plasma arginine (PP.H4 = {aa}) and with type 2 "
        u"diabetes ({a1}), and arginine colocalized with type 2 diabetes in FinnGen ({at}; {nsnp} variants). Both diabetes "
        u"edges are prior-sensitive ({a1} and {at} at the default prior, {a16} and {at6} at p12 = 1 × 10⁻⁶, below the "
        u"rule), and although the ARG1 and arginine maxima lie {bp930} bp apart, the type 2 diabetes signal is a broad "
        u"plateau, with variants within one log10 unit of its peak spread over about {span} kb, and its maximum lies about "
        u"{kb62} kb from the arginine lead, so the variant it shares with arginine is unresolved. The directions nonetheless "
        u"cohere: genetically higher ARG1 accompanied lower arginine, in the same ratio at both independent ARG1 signals "
        u"(arginine-to-ARG1 z-score ratio {r1} at the lead variant {rsid} and {r2} at the second signal), lower type 2 "
        u"diabetes risk ({mb} log-odds per SD; p = {mp}) and lower BMI ({bb} SD of BMI per SD of ARG1; p = {bp})."
    ).format(nn=F["node_n"], ng=F["node_genes"], adm=oxford(F["node_admitted"]), apo_h3="%.2f" % F["apoe_node_h3"],
             aa="%.3f" % F["bridge_top2"][1][2], a1="%.3f" % a1t["h4"], at="%.3f" % at["h4"], nsnp=thou(at["nsnp"]),
             r1=neg(F["arg1_ratio"]["lead"], 3), rsid=F["arg1_lead_rsid"], r2=neg(F["arg1_ratio"]["second"], 3),
             mb=neg(F["arg1_t2d_mr"]["b"], 3), mp=sci(F["arg1_t2d_mr"]["p"], 1), bb=neg(F["arg1_bmi"]["b"], 3),
             bp=sci(F["arg1_bmi"]["p"], 1), a16="%.3f" % a1t["h4_6"], at6="%.3f" % at["h4_6"],
             bp930="%d" % round(al["arg_prot_bp"]), span="%d" % (round(al["t2d_span_kb"] / 10.0) * 10), kb62="%d" % round(al["t2d_arg_kb"]))
    _need(a1t["h4_6"] < 0.8 and at["h4_6"] < 0.8 and F["arg1_t2d_mr"]["b"] < 0 and F["arg1_bmi"]["b"] < 0
          and F["arg1_ratio"]["lead"] < 0 and F["arg1_ratio"]["second"] < 0, "the ARG1 triangle's premises moved")
    (dpn, dpm, dph), (arn, arm, arh) = F["bridge_top2"]
    da = F["dpep1_attr"]
    dh = F["dpep1_halves"]
    # FIX round 4: the DPEP1 bridge is its own paragraph (F3a). Its gene-attribution account grew by the facts the
    # registered annotation table adds (the missense anchor, the intronic signal A, the transcript-set flip), and F3 would
    # have run to about 380 words across two bridges.
    kw = dict(nb=F["bsr_n"], bp=F["bsr_prot"], nbr=F["bridges"], nbp=WORDS[F["bridge_prot"]], dph="%.3f" % dph, arh="%.3f" % arh,
              cells=F["grid_cells"], oh3="%.2f" % F["orn_screen"]["h3"], mci=C("MCINTYRE"), half=dh["n"], dgene=dh["drop_gene"],
              fgene=dh["freq_gene"], dg=dh["glist_freq"], blk=WORDS[dh["blocked"]], b_freq=dh["b_freq"], b_drop=dh["b_drop"],
              aa=F["dpep1_ann"]["aa"], sa=F["dpep1_ann"]["sa"], lead=F["dpep1_ann"]["lead"],
              pc_kb="%.1f" % F["dpep1_ann"]["pc_kb"], m_pc=F["dpep1_ann"]["m_pc"], m_can=F["dpep1_ann"]["m_can"],
              far="%d" % round(F["dpep1_pos"]["b_drop_from_dpep1_bp"] / 1000.0))
    T["F3"] = (
        u"The bridge screen asked an open question — does any network protein share a causal variant with a nearby "
        u"metabolite? — and its two strongest answers paired enzymes with metabolites they act on (Figure 6b; "
        u"Supplementary Table S14). Of {nb} protein–metabolite pairs over {bp} network proteins, {nbr} colocalized, at "
        u"{nbp} proteins. The two highest posteriors were DPEP1–cysteinylglycine (PP.H4 = {dph}) and ARG1–arginine "
        u"({arh}), and in both the allele that raises the enzyme lowers the metabolite (bridge tier 1). ARG1–arginine, "
        u"the first edge of the triangle above, was attributed to ARG1 in all {cells} fine-mapping rule variants under both "
        u"gene annotations. Ornithine, the product of the arginase reaction, moved the opposite way, rising with ARG1 at a "
        u"signal the two shared in fine-mapping; it is not a bridge, because single-variant colocalization resolved it into "
        u"distinct variants (PP.H3 = {oh3}) and it did not track ARG1 in proportion."
    ).format(**kw)
    # FIX round 5 (verify:accuracy-narrative, round 4, MAJOR undersell): the kept half's dual colocalization and the
    # missense-led second signal are stated; the rule's single-signal reading is named as the reason for the Ensembl-116
    # result; the POST HOC per-signal sensitivity is reported in BOTH directions (any signal / every signal), with the
    # every-signal consequence for ARG1 stated where it bites. The rs409170-intron sentence (a 0.006-z tie) and the
    # 57/118-bp transcript-set flip (a base-pair margin the LD draws do not resolve) are dropped.
    ps = F["dpep1_ps"]
    p6 = F["ps6"]
    kw.update(h4lo=ps["h4_lo"], marg=ps["margin"], clsa=C("CLSA"), b_boot=p6["b_boot"], a_boot=p6["a_boot"],
              bmin=p6["bmin"], dmax=p6["dmax"])
    # FIX round 6 (round-5 verify:reviewer): the finding first (the kept half's dual colocalization at the missense
    # variant), the missense-as-support clause reconciled with the paper's own epitope limitation (the metabolite side is
    # mass spectrometry, the protein side affinity-based), the rule mechanics, the principled contrast (LD-draw stability
    # against BOOT_MIN, beside the ARG1 signal the rule reads), the drop half, and a verdict bounding DPEP1 over every rule,
    # reading and annotation. The every-signal result and ARG1's MED23 variant moved to Methods B5 ('whose lead signal
    # likewise peaks ... in MED23' was false: the lead pQTL rs2781668 lies in the ARG1/MED23 overlap); no per-signal
    # count is printed here (one direction only).
    # FIX round 6b (round-6 verify:numbers, two minors on one sentence): the draw count is the Ensembl-116 full-span column
    # (boot_n_qualifies_and_A2_ens116pc), so it now carries 'under full gene spans' -- the sentence before it names glist
    # too, and no PER-SIGNAL glist draw count is registered (the cells table's glist draw counts are the lead signal's and
    # the readings' unions, a different quantity). The sentence also stated only the shortfall, which reads as a flat
    # failure: the missense signal DOES meet the position condition on the full LD panel (evidence A2_ens116pc, already a
    # gate premise), and that positive half is now stated. The any-signal count 8 of 24 stays UNPRINTED: the per-signal
    # builder's frozen declaration is that neither direction is reported without the other, and the every-signal reading is
    # stated once, in B5 (round-6 M2) -- printing the 8 alone is what FORBID's '8 of 48' exists to stop.
    # PACKAGING PASS A (2026-09-30, round-6b verify:reviewer P1 / verify:numbers): 'under full gene spans' belonged to
    # the DRAW COUNT (boot_n_qualifies_and_A2_ens116pc is an Ensembl-116 full-span column) but round 6b attached it to
    # the full-panel clause, which is true under all three annotations and needs no scope; the count reached the
    # annotation only through 'does so'. The scope now sits on the sentence whose value it qualifies.
    T["F3a"] = (
        u"DPEP1–cysteinylglycine links the membrane dipeptidase to the glutathione-derived dipeptide that a rat renal "
        u"dipeptidase hydrolyses {mci}: in all {cells} rule variants the two shared a signal, and cysteinylglycine tracked "
        u"DPEP1 in the same ratio at DPEP1's two signals. Where palindromic variants whose strand the allele frequency "
        u"resolves are kept ({half} rule variants), cysteinylglycine colocalizes with both DPEP1 signals (conditional "
        u"PP.H4 > {h4lo} each), the second peaking at {b_freq}, a DPEP1 missense variant ({aa}). Because cysteinylglycine "
        u"was measured by mass spectrometry {clsa}, its association at that variant cannot reflect altered antibody "
        u"binding, although DPEP1's own association there could. The class rule places the signal in a gene using only "
        u"the more strongly colocalizing of the two, here DPEP1's lead signal by a PP.H4 margin below {marg}. Its top "
        u"shared variant, {lead}, lies between DPEP1 and {fgene}, about {pc_kb} kb from the protein-coding transcripts of "
        u"each. Full gene spans place it in {fgene}, through a transcript with no coding sequence, and attribute the "
        u"signal to DPEP1 in none of these {half}; glist-hg19 gene bodies attribute it in {dg}, the other {blk} blocked "
        u"because SuSiE, used as fitted, also colocalizes the measured {fgene} protein with cysteinylglycine. In the "
        u"per-signal sensitivity (Methods), the missense signal colocalizes and meets the position condition on the full "
        u"LD panel. Under full gene spans it does so in only {b_boot} of 200 LD-reference draws, short of the {bmin} "
        u"required, against {a_boot} for the ARG1 signal the class rule "
        u"reads. Dropping palindromic variants (the other {half}) removes {b_freq}, and DPEP1's second signal then leads "
        u"at {b_drop}, a DPEP1 pQTL variant inside {dgene}, about {far} kb from DPEP1. Cysteinylglycine shares that "
        u"signal, which no gene annotation assigns to DPEP1. The coding signal points to DPEP1, but no rule, reading or "
        u"annotation attributes the pair to DPEP1 in more than {dmax} of the {cells} variants."
    ).format(**kw)
    _need(da["ens"] == 0 and da["ensall"] == 0 and dh["glist_freq"] == da["glist"] and 2 * dh["n"] == F["grid_cells"],
          "'none of them attributes the signal to DPEP1' / the half counts no longer hold")
    br = F["bridge_rng"]
    ax_ = F["atxn2_h4"]
    T["F3b"] = (
        u"The other {nt2} bridges (bridge tier 2) keep one class in all {cells} variants, and each locates a shared signal "
        u"without naming its gene. BRAP's three metabolites, kynurenine, quinolinate and the uncharacterised X-12100 "
        u"(PP.H4 = {b0}–{b1}), share variants on the 12q24 SH2B3/ATXN2 haplotype, where BRAP's own cis association "
        u"never reaches genome-wide significance (maximum |z| = {bz}); neither measured neighbour explains them, since SH2B3 "
        u"does not share the signal and ATXN2, tested directly, colocalizes with none of them (PP.H4 = {x0}–{x1}). "
        u"SMPD3's four metabolites, lysine, N6,N6,N6-trimethyllysine, glutarylcarnitine and 3-amino-2-piperidone "
        u"({s0}–{s1}), share lead variants inside SLC7A6, which the proteomic panel does not measure. GSTA3's "
        u"1,3,7-trimethylurate ({g0}) shares a variant with GSTA3 but not in proportion, and its lead variants map to the "
        u"unmeasured GSTA2."
    ).format(nt2=WORDS[F["tier2_n"]], cells=F["grid_cells"], b0="%.3f" % br["BRAP"][0], b1="%.3f" % br["BRAP"][1],
             bz="%.2f" % F["brap_maxz"], x0="%.2f" % ax_[0], x1="%.2f" % ax_[-1], s0="%.3f" % br["SMPD3"][0],
             s1="%.3f" % br["SMPD3"][1], g0="%.3f" % br["GSTA3"][0])
    q = F["s26_q"]
    ncq = F["s26_coloc_q"]
    T["F4"] = (
        u"Across ancestries, the metabolite → disease associations that survive colocalization concentrate at the "
        u"urea cycle (Figure 6d). Of {n} metabolite–outcome cells, Mendelian randomization flagged {div} as "
        u"ancestry-divergent: {eur} with an association in Europeans only, {eas} in East Asians only and {dis} with opposite "
        u"signs. Colocalization in East Asians needs an East Asian association and could be run for {ev} of the {ea} cells "
        u"that had one; it held in {nc}, {cc} concordant across ancestries and {ce} East-Asian-only, and the single "
        u"opposite-sign cell was among the {nr} that could not be run. {cps} of the {ncw} sit at the CPS1 locus, which "
        u"encodes the first committed enzyme of the urea cycle and has been associated with plasma glycine and, in women, "
        u"with coronary artery disease {har}: glycine, creatine and three acylglycines at eGFR, and homoarginine and two "
        u"acylglycines at coronary artery disease. Because glycine is a precursor of creatine, and creatine of creatinine "
        u"{wys}, the eGFR cells may partly reflect creatinine production rather than filtration. The ninth is a FADS-locus "
        u"lysophospholipid, 1-linoleoyl-GPC, at type 2 diabetes."
    ).format(n=F["s26_n"], div=q["EUR-only"] + q["EAS-only"] + q["discordant"], eur=q["EUR-only"], eas=q["EAS-only"],
             dis=q["discordant"], ev=F["s26_eval"], ea=F["s26_easassoc"], nc=F["s26_coloc"], cc=ncq.get("concordant", 0),
             ce=ncq.get("EAS-only", 0), nr=F["s26_easassoc"] - F["s26_eval"], cps=WORDS[F["s26_cps1"]["n"]].capitalize(),
             ncw=WORDS[F["s26_coloc"]], har=C("HARTIALA"), wys=C("WYSS"))
    _need(q["discordant"] == 1 and F["s26_coloc"] == 9 and F["s26_cps1"]["n"] == 8, "the S26 prose premises moved")
    T["G1"] = "One colocalization gate distils every layer of the network"
    f9 = F["f9"]
    T["G2"] = (
        u"The same gate distils every layer of the analysis to a few anchors with interpretable biology (Figure 6e). Of "
        u"{l11} European protein–disease Mendelian-randomization tests, {l12} colocalize and {l13} reach the convergent "
        u"tier (PAM–type 2 diabetes, PCSK9–coronary artery disease, RAD23B–type 2 diabetes); of {l21} network "
        u"pairs tested natively in East Asians, {l22} colocalize and {l23} hold at the conservative prior "
        u"(CELSR2–coronary artery disease, FGF5–eGFR and TGFB1–coronary artery disease); and of {l31} "
        u"protein–metabolite pairs, {l32} colocalize and {l33} are enzyme–substrate bridges (ARG1–arginine, "
        u"DPEP1–cysteinylglycine). The final step of each funnel applies a different criterion — evidence tier, "
        u"prior robustness and enzyme–substrate chemistry — but the step that does most of the narrowing, the dual "
        u"rule on coloc.abf posteriors, is the same in every layer, and that gate, calibrated on the European "
        u"protein–disease layer, rather than the raw count of nominations, defines what the network can support."
    ).format(l11=thou(f9[("L1", "1")]), l12=f9[("L1", "2")], l13=f9[("L1", "3")], l21=f9[("L2", "1")], l22=f9[("L2", "2")],
             l23=f9[("L2", "3")], l31=f9[("L3", "1")], l32=f9[("L3", "2")], l33=f9[("L3", "3")])
    _need(F["cen_robust"] == [pair("CELSR2", "CAD"), pair("FGF5", "eGFR"), pair("TGFB1", "CAD")], "the robust three moved")
    T["H"] = (
        u"Measuring both sides of each association in East Asians showed which parts of the network travel. Three edges "
        u"held under the conservative prior, each at a locus with established biology. CELSR2–CAD is the edge from the "
        u"1p13 LDL-cholesterol locus, where a noncoding variant alters hepatic SORT1 expression and, through it, plasma LDL "
        u"cholesterol {mus}; its East Asian estimate matched the European one almost exactly, so the edge transfers as a "
        u"locus. FGF5–eGFR, at an established blood-pressure locus {nch}, gave the census's most decisive "
        u"colocalization, so the renal signal at this locus is shared across ancestries; in Europeans it runs opposite to "
        u"what blood-pressure mediation would require (see the intermediate analysis below). TGFB1–CAD colocalized in "
        u"both ancestries with opposite signs. TGF-β1 is a context-dependent cytokine: it restrains inflammation and "
        u"lesion growth in mouse atherosclerosis {mal} yet promotes fibrosis and neointima formation after vascular injury "
        u"{gra}. Here, however, the two ancestries' instruments are different variants {tg_kb} kb apart, so "
        u"ancestry-specific biology cannot be separated from instruments that tag different regulatory or assay effects; "
        u"either way, the direction of TGF-β1 modulation for coronary disease cannot be carried from one population to "
        u"another. IL6R–CAD, the positive control, was protective in both ancestries, the direction of the IL6R "
        u"Asp358Ala variant, which raises circulating IL-6 receptor and lowers coronary risk {il6}. Transfer across "
        u"ancestry is therefore selective rather than wholesale."
    ).format(mus=C("MUSUNURU"), nch=C("NEWTONCHEH"), mal=C("MALLAT"), gra=C("GRAINGER"), tg_kb="%d" % round(tg["dist_kb"]),
             il6=C("IL6RGC"))
    # FIX round 3 (verify:narrative): FGFR4-T2D gets its own paragraph (the network paragraph ran to 410 words and folded
    # the off-network pair into the network summary). It is characterised as in the Results: the East-Asian-specific
    # PATTERN, which is what the scan selected it for -- never a flat 'East-Asian-specific edge'. The instrument caveat is
    # said plainly: the instrument and Gly388Arg both sit at the smallest p the CKB resource reports and are in near-
    # complete LD (the tie and the r2 are derived above; the rsID count is not printed).
    _need(F["fgfr4_tied"] >= 2 and F["fgfr4_g388r_r2"] >= 0.9, "the FGFR4 instrument is no longer inseparable from Gly388Arg")
    T["Hb"] = (
        u"Off the European map, FGFR4–T2D colocalized in East Asians against a null European estimate, the pattern an "
        u"East Asian scan had selected it for (Results). FGFR4 is the hepatocyte receptor through which the ileal "
        u"hormone FGF19 represses bile-acid synthesis {kur}, and FGF19 also drives insulin-independent hepatic glycogen "
        u"synthesis {kir}. The East Asian direction, more FGFR4 with higher risk, matches the improved glucose tolerance "
        u"and insulin sensitivity of FGFR4-deficient mice made obese on a high-fat diet {ge}, not the glucose intolerance "
        u"and insulin resistance of FGFR4-deficient mice on a normal diet {hu}. Its instrument, however, cannot be "
        u"separated from the missense variant Gly388Arg (rs351855): both reach the smallest p-value the CKB resource "
        u"reports, and the two are in near-complete linkage disequilibrium in 1000 Genomes East Asians (r² = {r2}). "
        u"Gly388Arg alters the receptor's transmembrane segment and its signalling {ula}, so the association may reflect "
        u"an altered receptor or an altered assay epitope {kul} rather than FGFR4 abundance, and it needs orthogonal, "
        u"non–affinity-based confirmation."
    ).format(kur=C("KUROSU"), kir=C("KIR"), ge=C("GE"), hu=C("HUANG"), r2="%.2f" % F["fgfr4_g388r_r2"],
             ula=C("ULAGANATHAN"), kul=C("KULIESIUS"))
    T["I1"] = u"Enzyme–substrate bridges extend the network into the metabolome"
    T["I2"] = (
        u"Extending the network into the metabolome recovered two enzymes paired with metabolites they are documented to "
        u"act on. Arginase 1 hydrolyses arginine to ornithine and urea "
        u"in the last step of the hepatic urea cycle, and its inherited deficiency causes hyperargininaemia {c1}. "
        u"Genetically higher ARG1 accompanied lower plasma arginine in the same proportion at both of its independent "
        u"signals, and attribution to ARG1 held in all {cells} rule variants of the fine-mapping grid. DPEP1, one of the "
        u"two renal anchors (Results), encodes membrane-bound "
        u"dipeptidase {ad90}. A particulate rat renal dipeptidase hydrolyses cysteinylglycine, the dipeptide left when "
        u"γ-glutamyl transpeptidase acts on glutathione {mci}; the human enzyme hydrolyses S-substituted adducts such as "
        u"cystinyl-bis-glycine {ad89}, and mice lacking it excrete more of that adduct {hab}. Which form of the dipeptide "
        u"the enzyme chiefly clears is not settled: in rat renal microvillus membranes it was the principal activity "
        u"hydrolysing those S-substituted derivatives, and reduced cysteinylglycine inhibited it potently {koz}. "
        u"Genetically higher DPEP1 accompanied lower cysteinylglycine, the direction its documented activity predicts. "
        u"Cysteinylglycine also shares the DPEP1 signal that peaks at a missense variant in the enzyme itself, although "
        u"fine-mapping does not exclude the adjacent {fgene} (Results)."
    ).format(c1=C("SIN", "CALDWELL"), ad90=C("ADACHI1990"), mci=C("MCINTYRE"), koz=C("KOZAK"), ad89=C("ADACHI1989"),
             hab=C("HABIB"), cells=F["grid_cells"], fgene=F["dpep1_halves"]["freq_gene"])
    # FIX round 5: 'rather than its neighbours' was inaccurate (the kept half's competitor is the ADJACENT CHMP1A; the
    # dropped half's signal lies in VPS9D1, about 80 kb away with other genes between); 'at most a third' = the largest
    # DPEP1-cysteinylglycine count under any annotation or reading (16) is at most 48 / 3. FIX round 6: the DPEP1 sentence
    # leads with the finding; the bound now covers the every-signal reading too (ps6.dmax)
    _need(3 * F["dpep1_ps"]["max_any_primary"] <= F["grid_cells"] and 3 * F["ps6"]["dmax"] <= F["grid_cells"]
          and F["dpep1_halves"]["freq_gene"] == "CHMP1A",
          "DPEP1's attribution now exceeds a third of the rule variants, or the kept half's competitor is not CHMP1A")
    a26 = F["arg_t2d_s26"]
    T["I2b"] = (
        u"The ARG1 node also completes the network's only protein–metabolite–disease triangle, and its direction "
        u"is informative. Clinical studies cast arginase as harmful: arginase inhibition improves endothelial "
        u"function in patients with coronary artery disease and type 2 diabetes {she}, and oral arginine improved insulin "
        u"sensitivity in a small trial in type 2 diabetes {pia}. If circulating ARG1 acted mainly by depriving nitric oxide "
        u"synthase of arginine {cal}, more ARG1 should accompany higher type 2 diabetes risk. The genetic signal runs the "
        u"other way: more ARG1, less arginine and lower type 2 diabetes risk, together with lower BMI (Results), and the "
        u"protective ARG1 direction recurs in East Asian outcomes with the European instrument (Figure 3a). That pattern "
        u"argues against endothelial arginine scarcity as the link and points, as a hypothesis for direct testing, to "
        u"hepatic arginine turnover or adiposity. The triangle itself is European: its arginine–diabetes leg did not "
        u"recur in East Asians, where genetically predicted arginine showed no association with type 2 diabetes and "
        u"differed from the European estimate (heterogeneity p = {a26p}). In East Asians the metabolite–disease signals "
        u"point to the urea cycle from its other end: CPS1, which catalyses the cycle's first committed "
        u"step {har}, carries {cps} of the {ncw} metabolite–disease cells that colocalized there — {eg} at eGFR and {cad} at "
        u"coronary artery disease, none at type 2 diabetes."
    ).format(she=C("SHEMYAKIN"), pia=C("PIATTI"), cal=C("CALDWELL"), a26p="%.3f" % a26["hetp"], har=C("HARTIALA"),
             cps=WORDS[F["s26_cps1"]["n"]], ncw=WORDS[F["s26_coloc"]], cad=WORDS[F["s26_cps1"]["cad"]],
             eg=WORDS[F["s26_cps1"]["egfr"]])
    T["I3"] = (
        u"The {nt2} bridges of bridge tier 2 locate shared signals without naming their genes. At 12q24, kynurenine, "
        u"quinolinate and X-12100 share variants on the SH2B3/ATXN2 haplotype, an established blood-pressure locus {nch}, "
        u"where BRAP's own cis signal is sub-genome-wide; kynurenine and quinolinate are products of tryptophan catabolism "
        u"along the kynurenine pathway {cer}. At 16q22, lysine and N6,N6,N6-trimethyllysine map inside SLC7A6, which "
        u"encodes the cationic amino-acid exchanger y+LAT2 {bro} and is a more direct candidate for them than SMPD3. At "
        u"6p12, 1,3,7-trimethylurate, a product of caffeine oxidation {tas}, shares a variant with the glutathione "
        u"S-transferase cluster, with its lead variants mapping to GSTA2. The metabolite layer is a targeted colocalization "
        u"extension of the protein network rather than a metabolome-wide causal screen: only network proteins are "
        u"screened, so colocalizing metabolites at non-network genes fall outside it by design. Within that scope, its "
        u"contribution is two enzyme–substrate bridges, one of which completes a protein–metabolite–disease triangle, "
        u"obtained with the same gate, calibrated on the disease backbone, that defines the network's transfer across "
        u"ancestries."
    ).format(nt2=WORDS[F["tier2_n"]], nch=C("NEWTONCHEH"), cer=C("CERVENKA"), bro=C("BROER"), tas=C("TASSANEEYAKUL"))
    T["J"] = (
        u"the pQTL sidecar established exact-variant availability and direction but did not provide a complete native East "
        u"Asian MR and colocalization analysis. The native census supplied that analysis, but only for the {ct} testable "
        u"pairs of the {cn} whose protein carries a genome-wide-significant East Asian cis-pQTL, and most of its "
        u"non-colocalizing pairs ({h1} of {nc}) lacked an East Asian outcome association, which cannot distinguish limited "
        u"power from absent transfer."
    ).format(cn=F["cen_n"], ct=F["cen_test"], h1=F["cen_h1"], nc=F["cen_nonc"])
    rb = F["cen_robust"]
    T["K1"] = (
        u" natively, {ns} of {ct} colocalized, {nr} at a conservative prior ({r0}, {r1}, sign-reversed {r2})."
    ).format(ns=len(sup), ct=F["cen_test"], nr=len(rb), r0=rb[0], r1=rb[1], r2=rb[2])   # FIX round 3: numerals throughout
    # one comparison (the abstract's house style: '21 as concordant, 4 as discordant and 3 as conflicting')
    _need(rb[2] == pair("TGFB1", "CAD") and F["tgfb1"]["eur"] > 0 > F["tgfb1"]["eas"], "the sign-reversed robust edge moved")
    # FIX round 2 (USER decision B): the words freed by moving FGFR4-T2D out of the abstract carry the ARG1 direction.
    # Both edges are derived: ARG1 -> arginine from the class table's z-score ratio at the cis lead, ARG1 -> T2D from the
    # network edge's beta. 'accompanied' (not an arrow chain): no mediation analysis was run.
    _need(F["arg1_ratio"]["lead"] < 0 and F["arg1_t2d_mr"]["b"] < 0, "the ARG1 abstract direction is no longer negative on "
          "both edges (class-table ratio at the cis lead; network-edge beta)")
    # FIX round 3 (verify:narrative): 'robust' had lost its referent (framing decision 4: 'robust to fine-mapping'), and
    # 'lower arginine and T2D risk' was elliptical. 350 words by Word (P11a).
    _need(max(F["dpep1_attr"].values()) < F["grid_cells"], "DPEP1's gene attribution is robust under some annotation now "
          "(ARG1's 48/48 under every annotation is asserted in facts())")
    T["K2"] = (
        u" The two strongest of {nbr} protein–metabolite bridges were enzyme–substrate pairs (ARG1–arginine, "
        u"DPEP1–cysteinylglycine); only ARG1's gene attribution was robust to fine-mapping. Genetically higher ARG1 "
        u"accompanied lower arginine and lower T2D risk."
    ).format(nbr=F["bridges"])
    # T["M5"] (the AGEN abbreviation insert) is RETIRED: PASS A generates the whole list (abbreviations.py, ADD["AGEN"])
    return T


def _num_word_cap(n):
    """a count opening a sentence: spelled out."""
    ones = ["", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]
    tens = ["", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"]
    if n < 10:
        w = ones[n]
    elif n < 20:
        w = ["ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen",
             "nineteen"][n - 10]
    else:
        w = tens[n // 10] + (("-" + ones[n % 10]) if n % 10 else "")
    return w.capitalize()


def inherited_edits():
    """[(anchor, old, new)] for the inherited (stage-1) paragraphs this update touches; old/new are STAGE-1 text
    (citation numbers are stage-1). Run by rebase_v4_edits.main_edits()."""
    F, _ = facts()
    E = []
    from sidecar_claims import facts as _scf
    tt = _scf()["tt"]
    ev, comp = tt["eas_eval"], tt["eas"].get("EAS_COMPATIBLE", 0)
    _need(tt["eas"].get("EAS_DISCORDANT", 0) == 0, "an outcome-side pair is discordant: 'no discordant edge' is false")
    # abstract: the outcome-side/pQTL sentence -> one sentence (the native sentence K1 follows it; spec [main 7])
    E.append(("Results. The network comprised",
              "East Asian outcomes were evaluable for 46 edges (12 compatible, 34 inconclusive, 0 discordant), and 42 edges "
              "had a strong exact-variant EAS pQTL (41 directionally concordant).",
              u"East Asian outcomes showed no discordant edge (%d of %d compatible);" % (comp, ev)))
    # Limitations overlap paragraph: two sentences inserted (spec [main 128]); old is never a substring of new
    R = refs()
    E.append(("Residual participant overlap with UK Biobank is confined by construction",
              "overlaps UKB-PPP. The supplementary cystatin-C",
              u"overlaps UKB-PPP. The European arm of the metabolite layer uses the same eight outcome files, with "
              u"metabolite QTLs from the Canadian Longitudinal Study on Aging [%d]. The supplementary cystatin-C" % R["CLSA"]))
    E.append(("Residual participant overlap with UK Biobank is confined by construction",
              "sample-independent replication. Overlap-induced bias",
              u"sample-independent replication. On the East Asian side, the CKB exposure cohort also contributes to %s "
              u"outcome meta-analyses, and a CKB-free re-run kept all %s colocalizations and added %s (Results). "
              u"Overlap-induced bias"
              % (WORDS[F["ckb_ov_resources"]], WORDS[F["ckbfree_kept"]], oxford(F["ckbfree_gained"]))))
    # List of abbreviations: PP.H1 (first used and glossed in the census Results). PACKAGING PASS A (2026-09-30): this
    # edit is RETIRED -- the whole list is generated by abbreviations.py (REKEY carries the PP.H1 entry), so an anchored
    # edit on the paragraph would no longer find its target.
    # Availability of data and materials: the census outcome resources, the census instruments, the metabolite QTLs
    ac = F["cen_acc"]
    E.append(("All data analysed here are from existing resources",
              "the East Asian outcome-side resources used for the transferability axis (Biobank Japan atrial fibrillation, "
              "coronary artery disease, heart failure and ischaemic stroke; T2DGGI East Asian type 2 diabetes; the East Asian "
              "eGFR meta-analysis; and GCST90018602 chronic kidney disease), all public; Japan COVID-19 Task Force and China "
              "Kadoorie Biobank [%d] pQTL resources for the exact-instrument portability sidecar;" % R["CKB"],
              u"the East Asian outcome-side resources used for the transferability axis (Biobank Japan atrial fibrillation, "
              u"coronary artery disease, heart failure and ischaemic stroke; T2DGGI East Asian type 2 diabetes; the East Asian "
              u"eGFR meta-analysis; and GCST90018602 chronic kidney disease), all public; the East Asian outcomes of the "
              u"native census (the AGEN type 2 diabetes meta-analysis, %s [%d]; the GIGASTROKE East Asian any-stroke analysis, "
              u"%s [%d]; the Taiwan–Japan eGFR meta-analysis [%d]; and the Biobank Japan coronary artery disease [%d], "
              u"atrial fibrillation [%d], heart failure (%s) [%d], chronic kidney disease (%s), type 2 diabetes and "
              u"ischaemic stroke [%d] analyses, from https://pheweb.jp/ and the GWAS Catalog), all public; Japan COVID-19 Task Force "
              u"and China Kadoorie Biobank [%d] pQTL resources for the exact-instrument portability sidecar, the China "
              u"Kadoorie Biobank per-region interface also supplying the native census instruments; plasma metabolite QTLs "
              u"from the Canadian Longitudinal Study on Aging [%d] and the Han Chinese serum metabolome [%d] for the "
              u"metabolite layer; the 1000 Genomes phase 3 European and East Asian reference panels for fine-mapping LD and "
              u"the census allele-frequency check;"
              % (ac["T2D"], R["SPRACKLEN"], ac["Stroke"], R["MISHRA"], R["CHEN2024"], R["KANAI2021"], R["MIYAZAWA"],
                 ac["HF"], R["ENZAN"], ac["CKD"], R["SAKAUE"], R["CKB"], R["CLSA"], R["CHENG"])))
    # ---- FIX round (2026-09-27, verify:underselling): the paper's FRAME carries the two new layers --------------------
    # abstract Methods: the native census and the metabolite screen named (two sentences tightened to hold the 350-word cap)
    E.append(("Methods. We screened",
              "SuSiE, sensitivity analyses and empirical calibration evaluated locus support. deCODE and ARIC assessed exposure "
              "replication and platform robustness; East Asian outcome GWAS and JCTF/CKB pQTL resources assessed outcome and "
              "instrument transportability. All backbone edges were retained in a pre-specified three-tier framework.",
              u"SuSiE, sensitivity analyses and empirical calibration evaluated locus support; deCODE and ARIC assessed "
              u"replication and platform robustness. East Asian outcome, pQTL and native colocalization analyses assessed "
              u"transportability; protein–metabolite colocalization probed the metabolome. Backbone edges were "
              u"graded in three pre-specified tiers."))
    # FIX round 2: the abstract Methods must not state a result (it said the screen 'extended the network into the
    # metabolome'). The permutation z is RESTORED (USER decision B) simply by NOT deleting it: apply_blocks' edit S puts
    # the abstract's null mean and z on the 20,000-permutation run (88.2, z = 5.7), and round 1 then removed the z to hold
    # the word cap. The words freed by moving FGFR4-T2D out of the abstract pay for it, so that deletion is gone. Assert
    # here that the run behind edit S still reads 88.2 / 5.7, since the sidecar claim S1 now prints both.
    sw = _scf()["sweep"]["cell"][("protein", 0.01)]
    _need("%.1f" % float(sw["null_mean"]) == "88.2" and "%.1f" % float(sw["z"]) == "5.7" and int(sw["n_perm"]) == 20000,
          "the abstract's permutation cell moved (%r, z %r over %r permutations): edit S and sidecar claim S1 print it"
          % (sw["null_mean"], sw["z"], sw["n_perm"]))
    # abstract Conclusions: the two layers, finding-led
    E.append(("Conclusions. The evidence-calibrated network",
              "cardiometabolic diseases form an interconnected component, while renal outcomes occupy a distinct axis. Its "
              "graded edges prioritise cross-disease hypotheses while making uncertainty explicit; they do not establish shared "
              "causal proteins or validated therapeutic targets.",
              u"cardiometabolic diseases form an interconnected component and renal outcomes a distinct axis. %s edges held "
              u"natively in East Asians, one sign-reversed; enzyme–substrate bridges reach the metabolome; graded edges "
              u"prioritise cross-disease hypotheses, not validated targets."
              % WORDS[len(F["cen_robust"])].capitalize()))
    # Background objective
    # (P22 requires the old text GONE: every insertion below spans a sentence boundary, so old is never inside new)
    E.append(("Here our primary objective was to construct",
              "then defined how much confidence each edge could carry. Drug tractability",
              u"then defined how much confidence each edge could carry. Two further analyses asked whether network edges "
              u"recur when both sides of each association are measured in East Asians, and whether the network extends, "
              u"through shared genetic signals, into the plasma metabolome. Drug tractability"))
    # Study Design: the pipeline list
    E.append(("The pipeline comprised",
              "(vi) East Asian outcome transportability and an EAS pQTL portability sidecar; (vii) evidence-strength tiering; "
              "(viii) network construction and topology testing; (ix) empirical operating-characteristic calibration; and (x) "
              "translational and intermediate-phenotype annotation.",
              u"(vi) East Asian outcome transportability, an East Asian ancestry (EAS) pQTL portability sidecar and a native "
              u"East Asian colocalization census; (vii) evidence-strength tiering; (viii) network construction and topology "
              u"testing; (ix) "
              u"empirical operating-characteristic calibration; (x) translational and intermediate-phenotype annotation; and "
              u"(xi) metabolite-layer colocalization."))
    # Discussion opening: the two layers after the topology sentence
    rb = F["cen_robust"]
    _need(rb == [pair("CELSR2", "CAD"), pair("FGF5", "eGFR"), pair("TGFB1", "CAD")] and F["tier1"] ==
          [("ARG1", "arginine"), ("DPEP1", "cysteinylglycine")], "the Discussion opening's named edges moved")
    E.append(("We constructed a genetically anchored, evidence-calibrated plasma-protein network spanning eight CKM outcomes",
              "surviving false-discovery correction. The central contribution is therefore",
              u"surviving false-discovery correction. Measured natively in East Asians, %d of the %d testable network edges "
              u"colocalized and %d held under a conservative prior: a lipid-locus edge almost unchanged, a renal edge "
              u"recurring and a TGF-β1 edge reversed in sign. "
              u"The same colocalization gate carried the network into the metabolome through two enzyme–substrate bridges, "
              u"ARG1–arginine and DPEP1–cysteinylglycine, the first completing a European protein–metabolite–disease "
              u"triangle in which genetically higher ARG1 accompanied lower arginine and lower type 2 diabetes risk. "
              u"The central contribution is"
              % (len(F["cen_coloc"]), F["cen_test"], len(rb))))   # FIX round 3: numerals throughout one comparison
    # FIX round 2 (verify:narrative): 'therefore' followed the new metabolome sentence; the list also omitted both layers
    E.append(("We constructed a genetically anchored, evidence-calibrated plasma-protein network spanning eight CKM outcomes",
              "in which network topology and the credibility of each edge are reported together.",
              u"in which network topology, the credibility of each edge, and its reach across ancestry and into the "
              u"metabolome are reported together."))
    # the external-axes paragraph: the census provides the native estimate the sidecar lacks
    E.append(("The external resources calibrate the network because they interrogate different failure modes",
              "but it still does not provide fully native East Asian protein-to-disease MR.",
              u"but it does not estimate the protein–disease association itself, which the native census does, with a native "
              u"instrument, for %d of the %d pairs whose protein carries a genome-wide-significant East Asian cis-pQTL."
              % (F["cen_test"], F["cen_n"])))
    # main Conclusions
    E.append(("We constructed an evidence-calibrated plasma-protein network that provides",
              "rather than a single replication count. The resulting network",
              u"rather than a single replication count. Measured natively in East Asians, %d of the %d testable edges "
              u"colocalized and %d held under a conservative prior, one (TGFB1–CAD) with the opposite sign. The same "
              u"colocalization gate carried the network into the metabolome through two enzyme–substrate bridges, "
              u"ARG1–arginine and DPEP1–cysteinylglycine, the first linking higher ARG1 to lower arginine and lower type 2 "
              u"diabetes risk in Europeans. The resulting network"
              % (len(F["cen_coloc"]), F["cen_test"], len(rb))))   # FIX round 3: the 7 colocalizations are reported too
    return E


# The East Asian census outcome rows appended to Table 1 (Outcome | Dataset identifier | Consortium / Source | N | Type).
def table1_rows():
    F, _ = facts()
    n, a = F["t1_n"], F["t1_acc"]
    return [
        ["T2D (EAS census)", a["T2D"], "AGEN (Spracklen 2020)", thou(n["T2D"]), "Disease"],
        ["Stroke (EAS census)", a["Stroke"], "GIGASTROKE East Asian (Mishra 2022)", thou(n["Stroke"]), "Disease"],
        ["eGFR (EAS census)", a["eGFR"], u"Taiwan–Japan meta-analysis (Chen 2024)", thou(n["eGFR"]), "Continuous trait"],
        ["CAD (EAS census)", a["CAD"], "Biobank Japan (Kanai 2021)", thou(n["CAD"]), "Disease"],
        ["AF (EAS census)", a["AF"], "Biobank Japan (Miyazawa 2023)", thou(n["AF"]), "Disease"],
        ["Heart failure (EAS census)", a["HF"], "Biobank Japan (Enzan 2025)", thou(n["HF"]), "Disease"],
        ["CKD (EAS census)", a["CKD"], "Biobank Japan (Sakaue 2021)", thou(n["CKD"]), "Disease"],
        ["T2D (EAS, CKB-free)", a["T2D_BBJ"], "Biobank Japan (Sakaue 2021)", thou(n["T2D_BBJ"]), "Disease"],
        ["Ischaemic stroke (EAS, CKB-free)", a["Stroke_BBJIS"], "Biobank Japan (Sakaue 2021)", thou(n["Stroke_BBJIS"]), "Disease"],
    ]


# Supplementary Tables S13 and S14 (the supplement's table list; apply_supplement.py prints this, row counts asserted there)
def supp_tables():
    F, _ = facts()
    return (
        u" - Supplementary Table S13 — Native East Asian colocalization census: for each pair, the East Asian "
        u"instrument and outcome resource, European and East Asian Wald estimates with heterogeneity tests, colocalization "
        u"posteriors at the default and conservative priors, and the CKB-free and European-prior sensitivities "
        u"(S23_native_census__source.csv, {rows} rows: the {cn} network pairs, {ct} of them testable, and {pn} off-network "
        u"pilot pairs). - Supplementary Table S14 — Protein–metabolite bridge screen and bridge fine-mapping classes "
        u"(bridge_screen_results.csv, {nb} rows; bridge_class_summary.csv, {nc} rows; bridge_persignal_sensitivity.csv, "
        u"{nps} rows). The screen lists every "
        u"protein–metabolite pair nominated by cis proximity, {nb} pairs over {bp} network proteins and {bm} metabolites, "
        u"each tested on the UKB-PPP European discovery arm; {nbr} meet the colocalization rule. The class table gives, for "
        u"the {nbr} bridges, their follow-up and neighbouring-protein pairs and one pipeline control, the number of the "
        u"{cells} rule variants supporting each class under Ensembl release 116 (glist-hg19 as sensitivity), the direction "
        u"at the shared signal and the bridge tier. The per-signal table gives, for the same pairs, the number of those "
        u"variants in which the pair is gene-attributed when any colocalizing protein signal, or every one, must meet the "
        u"position condition (a post hoc sensitivity; Methods)."
    ).format(rows=F["s23n_rows"], cn=F["cen_n"], ct=F["cen_test"], pn=F["pil_n"], nb=F["bsr_n"], nc=F["bcs_n"],
             bp=F["bsr_prot"], bm=F["bsr_metab"], nbr=F["bridges"], cells=F["grid_cells"], nps=len(_rows(PSS)))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    F, SRC = facts()
    if "--texts" in sys.argv:
        for k, v in texts().items():
            print("### %s\n%s\n" % (k, v))
        for a, o, nw in inherited_edits():
            print("--- EDIT under %r\n  old: %s\n  new: %s\n" % (a[:40], o, nw))
        print(supp_tables())
    else:
        for k in F:
            print("%-22s = %r   [%s]" % (k, F[k], SRC[k]))
