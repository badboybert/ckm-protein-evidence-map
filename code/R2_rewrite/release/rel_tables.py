# -*- coding: utf-8 -*-
"""Stage 5 step 2 -- release module for the two table TSVs (spec section 6, plan S2-D29 card specification).

Writes
  out/tables/Table1_data_sources.tsv        one row per Table 1 body row (shipped)
  out/tables/Table2_evidence_cards.tsv      one row per Table 2 body row, TPMI- and Open Targets-derived data cells replaced by `withheld_pending_terms`
  out/tables/Table2_evidence_cards.columns.tsv   one-line definition of every column of both table files (file, column, definition)
  withheld/tables/Table2_evidence_cards.tsv the full version (never shipped until the provider terms are read)
  withheld/tables/Table2_evidence_cards.withhold_rules.tsv   every rule that turns the full file into the shipped file, with the withheld text

Interface (spec section 2): INPUTS (V3-relative paths the module reads) and build(v3, out, withheld, ledger) -> list[dict].
Module switch: STRICT (default True). The core sets `rel_tables.STRICT = False` under --stage2. Strict refuses a `[HOLD key: ...]` in Table 1 or
Table 2 whose key is not also held in the Results draft (an orphan hold); a hold that Results keeps too is passed through verbatim. Under --stage2
every hold is passed through. Either way the holds found are listed in HOLDS_FOUND and in each returned note.

Rules of this module (spec section 1): no new statistic and no typed value. Every number in a numeric cell is the ledger row's own display token
(the printed precision), asserted to round from the stored leaf; the only literals are ledger row ids, leaf positions and labels. Tag handling is the
manuscript builder's own (P31R_build.TAG, P31R_tablegate.parse_table), imported, not re-implemented. Ledger tags are stripped exactly as the builder
strips them; `[@key]` citations stay inline (their manuscript numbers depend on a full build) and are also listed in `cited_keys`/`literature_pmids`.

Not written here (core): the literature search log copy (tables/literature_search_log.tsv) and the Open Targets snapshot (withheld/tables/); this module only reads them.
Redaction rule of the shipped Table 2 (declared, asserted, reversible; S2-D30(3), amended by S2-D34(1) and (4), 2026-10-18): only DATA columns are withheld.
  * Open Targets (S2-D34(4)): the Open Targets-native fields ship. A cell of development_stage, development_agents, genetic_direction_agreement or development_statement reads
    `withheld_pending_terms` only where it is a ChEMBL-derived drug or stage field: a stage set by a drug record ('approved', 'in clinical development'), a non-empty agent list, a
    statement or direction that names a drug, or a statement of a row that has an agent. The tractability class of a row with no clinical candidate, the one-word genetic-direction
    verdict that names no agent, and the data version and query date ship (rule OT1 lists every withheld cell).
  * TPMI (S2-D34(1)): with the TPMI switch (TIER_RULES.tsv SW1) at 'ship' the derived per-pair TPMI values (ext_* cells and the TPMI note) ship; the raw per-variant TPMI values
    (beta, SE, allele frequency at the variant) are not in any Table 2 column, and a scan asserts that none appears. With SW1 at 'withheld' they read `withheld_pending_terms` (rules TPMI1, TPMI3).
  The full file and the rules are in withheld/. The four prose columns ship verbatim: they equal the built Table 2 after tag resolution and the
withheld card's prose columns, byte for byte (asserted). Table1_data_sources.tsv mirrors the manuscript's Table 1 and ships whole, TPMI case and control counts included.

Coverage (F01/F03/F05): every pair x numeric column is a ledger leaf or a declared `not reported` with a reason code (NR_DECL); every ledger row whose key or claim names the
pair (gene and outcome tokens) and that no card column uses is on the declared NOT_CARRIED list with a reason, or the build fails; the S2-D29 must-carry items are asserted present.

Standalone:  python -I -X utf8 rel_tables.py --build OUT_DIR WITHHELD_DIR     (ledger read from the default path; use a folder outside the project)
             python -I -X utf8 rel_tables.py --selftest
"""
import sys
import ast, csv, hashlib, io, json, os, pathlib, re, shutil, tempfile

PRODUCER = "rel_tables.py"
STRICT = True
HOLDS_FOUND = []
TPMI_SHIP = True          # set by _build from TIER_RULES.tsv row SW1 (target TPMI_DERIVED, handling ship|withheld); read, not typed
TIER_RULES_REL = "R2_rewrite/release/inputs/TIER_RULES.tsv"
WITHHELD = "withheld_pending_terms"
NOT_REPORTED = "not reported"

LEDGER_REL = "R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv"
T1_REL = "R2_rewrite/tables/P31R_Table1.md"
T2_REL = "R2_rewrite/tables/P31R_Table2.md"
CARDS_REL = "R2_rewrite/supplement/P31R_SM_cards.md"   # Stage 7 A2.4c (S2-D51): the full Table 2 cards, Additional file 2 {SM:cards}
PROV_REL = "R2_rewrite/tables/P31R_Table2_provenance"
LOG_REL = PROV_REL + "/literature_search_log.tsv"
SNAP_REL = PROV_REL + "/devstage_opentargets_20261006.json"
OC_REL = "R1_stage1/S1R1_05b_outcome_checks/S1R1_outcome_checks.tsv"
EXT2_REL = "R1_stage1/S1R1_06_claims_ledger/ext_table2.py"
RES_REL = "R2_rewrite/sections_R2/P31R_Results_R2_draft.md"
OLDREFS_REL = "v3.2/prose/_refs_resolved.tsv"
KEYTAB_REL = "R2_rewrite/P31R_REFERENCE_KEYS.tsv"
REFS_RELS = ["R2_rewrite/refs/new_refs_%s.tsv" % n for n in ("af1", "supp1", "supp2", "table1", "table2")]
SCRIPT_RELS = ["R2_rewrite/scripts/P31R_build.py", "R2_rewrite/scripts/P31R_tablegate.py", "R2_rewrite/scripts/P31R_gate.py",
               "R2_rewrite/scripts/R2_d1_display_order.tsv"]

# Declared input list (the core pins and checks these md5s before calling build). The ledger CSV is also passed in as a DataFrame (the core reads it); it is listed here so that
# every returned record carries it in source_paths and source_md5 (F28).
INPUTS = [LEDGER_REL, T1_REL, T2_REL, CARDS_REL, LOG_REL, SNAP_REL, OC_REL, EXT2_REL, RES_REL, OLDREFS_REL, KEYTAB_REL] + REFS_RELS + SCRIPT_RELS

OUT_T1 = "tables/Table1_data_sources.tsv"
OUT_T2 = "tables/Table2_evidence_cards.tsv"
OUT_RULES = "tables/Table2_evidence_cards.withhold_rules.tsv"
OUT_DOC = "tables/Table2_evidence_cards.columns.tsv"

T1_NAMES = ["role", "resource", "ancestry", "trait_and_scale", "cases_controls_or_n", "overlap_with_ukb", "overlap_with_finngen",
            "ldsc_noise_correlation_with_finngen_outcomes"]
T1_EXTRA = ["provider", "licence_class", "cited_keys", "cited_pmids", "ledger_row_ids"]
T1_COLS = T1_NAMES + T1_EXTRA

T2_PROSE = ["established_role_reported_before", "adds", "main_uncertainty", "next_biological_question"]
# Stage 7 A2.4c (S2-D51, review R17): Table 2 prints a summary card; the full cards are Additional file 2's section {SM:cards} (supplement/P31R_SM_cards.md).
# Their four cells ship as full_* columns, read by row label, ledger tags removed (S2-D30(3) applied to the full card).
T2_FULL = ["full_" + _k for _k in T2_PROSE]
CARD_PREFIX = ["*Established role; reported before.* ", "*What this study observes or adds.* ", "*Main uncertainty.* ", "*Next biological question.* "]
# Stage 6.2 header of Table 2 (S2-D37 D2): the card prose columns are these four cells, in this order; the build asserts the header it reads is exactly this
T2_HEADER = ["Pair; evidence tier", "Established role; reported before", "What this study observes or adds", "Main uncertainty", "Next biological question"]
T2_COLS = ["pair_key", "gene", "outcome", "row_group", "table_row_label", "evidence_tier",
           "deCODE_arm", "ARIC_arm", "EAS_state", "card_note",
           "discovery_basis", "discovery_wald", "discovery_p", "pp_h4_default", "pp_h4_conservative",
           "ext_resource", "ext_or", "ext_ci_lo", "ext_ci_hi", "ext_p", "ext_wald", "ext_wald_se", "ext_wald_p", "ext_class", "ext_bonferroni_6",
           "reg_pp_h4_default", "reg_pp_h4_conservative", "reg_pp_h3_default", "reg_pp_h3",
           "eas_basis", "eas_beta", "eas_pp_h1", "eas_pp_h3", "eas_pp_h4", "eas_pp_h4_default", "eas_pp_h4_conservative",
           "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4", "eur_p", "candidate_susie_pp_h4", "candidate_coloc_abf_pp_h4",
           "development_stage", "development_agents", "genetic_direction_agreement", "development_statement", "ot_data_version", "ot_query_date",
           "literature_status", "literature_pmids", "search_date", "search_scope", "cited_keys"] + T2_PROSE + T2_FULL + \
          ["numeric_cell_sources", "not_reported_reasons", "ledger_rows_not_carried", "ledger_row_ids", "licence_class", "source_md5"]
T2_NUMERIC = ["discovery_wald", "discovery_p", "pp_h4_default", "pp_h4_conservative", "ext_or", "ext_ci_lo", "ext_ci_hi", "ext_p", "ext_wald", "ext_wald_se", "ext_wald_p",
              "reg_pp_h4_default", "reg_pp_h4_conservative", "reg_pp_h3_default", "reg_pp_h3", "eas_beta", "eas_pp_h1", "eas_pp_h3", "eas_pp_h4",
              "eas_pp_h4_default", "eas_pp_h4_conservative", "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4", "eur_p",
              "candidate_susie_pp_h4", "candidate_coloc_abf_pp_h4"]
# Data columns withheld in the shipped card (S2-D30(3): the four prose columns are NOT withheld). development_statement is the ledger word statement of stage, agents and direction.
T2_OT_COLS = ["development_stage", "development_agents", "genetic_direction_agreement", "development_statement"]
T2_TPMI_COLS = ["ext_or", "ext_ci_lo", "ext_ci_hi", "ext_p", "ext_wald", "ext_wald_se", "ext_wald_p", "ext_class", "ext_bonferroni_6"]

# (column, ledger row id, leaf).  leaf: int = flat numeric leaf of the stored value (PG._d1_leaf_list order); ("n", j) = element j of the pair's
# row in L223 (nested list, pair in element 0); ("p",) = the pair's leaf in L960 (leaf order = the first six rows of S1R1_outcome_checks.tsv).
# Which stored leaf a column reads is a declared structure; the value, its printed precision and the sign are read and asserted.
NUMSPEC = {
    "PCSK9-CAD": [("discovery_wald", "L350", 0), ("discovery_p", "L351", 0), ("pp_h4_default", "L202", 0), ("pp_h4_conservative", "L352", 0),
                  ("reg_pp_h4_default", "L746", 0), ("reg_pp_h4_conservative", "L747", 0), ("reg_pp_h3", "L745", 0), ("eas_pp_h3", "L361", 0)],
    "TGFB1-CAD": [("discovery_wald", "L085", 0), ("discovery_p", "L353", 0), ("pp_h4_default", "L354", 0), ("pp_h4_conservative", "L214", 2),
                  ("reg_pp_h4_default", "L752", 0), ("reg_pp_h4_conservative", "L753", 0),
                  ("eas_beta", "L085", 1), ("eas_pp_h4_default", "L214", 0), ("eas_pp_h4_conservative", "L214", 1)],
    "PAM-T2D": [("discovery_wald", "L355", 0), ("discovery_p", "L356", 0), ("pp_h4_default", "L357", 0), ("pp_h4_conservative", "L358", 0),
                ("reg_pp_h4_default", "L743", 0), ("reg_pp_h4_conservative", "L744", 0), ("eas_pp_h1", "L362", 0), ("eas_ckbfree_pp_h1", "L363", 0)],
    "RAD23B-T2D": [("discovery_wald", "L359", 0), ("discovery_p", "L204", 0), ("pp_h4_default", "L360", 0), ("pp_h4_conservative", "L204", 1),
                   ("reg_pp_h4_default", "L750", 0), ("reg_pp_h4_conservative", "L751", 0), ("reg_pp_h3_default", "L748", 0), ("reg_pp_h3", "L749", 0)],
    "ARG1-T2D": [("discovery_wald", "L050", 0), ("discovery_p", "L050", 1), ("pp_h4_default", "L051", 0), ("pp_h4_conservative", "L051", 1),
                 ("eas_pp_h3", "L266", 0), ("eas_pp_h4", "L266", 1), ("eas_ckbfree_pp_h1", "L266", 2), ("eas_ckbfree_pp_h3", "L266", 3),
                 ("eas_ckbfree_pp_h4", "L266", 4)],
    "FGFR4-T2D (EAS)": [("discovery_wald", "L086", 0), ("discovery_p", "L086", 1), ("pp_h4_default", "L086", 2), ("pp_h4_conservative", "L215", 0),
                        ("eas_ckbfree_pp_h4", "L086", 3), ("eur_p", "L215", 1)],
    "UMOD-eGFR": [("discovery_wald", "L364", 0), ("discovery_p", "L071", 0), ("eas_pp_h4_default", "L075", 0),
                  ("candidate_susie_pp_h4", "L207", 0), ("candidate_coloc_abf_pp_h4", "L207", 2)],
    "DPEP1-eGFR": [("discovery_wald", "L365", 0), ("discovery_p", "L071", 2), ("candidate_susie_pp_h4", "L207", 1), ("candidate_coloc_abf_pp_h4", "L207", 3)],
}
EXT_SPEC = [("ext_or", "L223", ("n", 1)), ("ext_ci_lo", "L223", ("n", 2)), ("ext_ci_hi", "L223", ("n", 3)), ("ext_p", "L960", ("p",))]
# The pair's external-check row ("<pair>: Wald ratio in <resource> ...", the row that holds the class word) also stores the first-order Wald ratio, its SE and p as numeric leaves 0-2 and the
# discovery estimate of the same quantity as leaf 3 (not printed). The row id is looked up by that claim, never typed; the leaf positions are the declared structure.
EXT_WALD_SPEC = [("ext_wald", 0), ("ext_wald_se", 1), ("ext_wald_p", 2)]
EXT_WALD_DISCOVERY_LEAF = 3
EXT_PAIRS = ["PCSK9-CAD", "TGFB1-CAD", "PAM-T2D", "RAD23B-T2D", "ARG1-T2D", "FGFR4-T2D (EAS)"]   # the six pairs of L223 / L960 (Table 2 rows 1-6)
EXT_DISPLAY = {"DIAMANTE": "DIAMANTE", "CARDIOGRAM": "CARDIoGRAMplusC4D", "TPMI": "TPMI"}   # outcome-checks `resource` -> name Table 2 prints
DISCOVERY_BASIS_COL = "discovery_wald"

# roles of the elements in the L370-L377 word rows, by row group (the element lists are asserted against these lengths)
ROLES = {"pre-selected external-check pair": ["tier", "deCODE_arm", "ARIC_arm", "EAS_state", "stage_dir", "lit"],
         "East Asian, not graded": ["grade", "exposure_side", "ext_check_note", "instrument_note", "stage_dir", "lit"],
         "multi-signal candidate, not graded": ["grade", "EAS_state", "stage_dir", "lit"]}
TIERS = ("Convergent", "Supported")

# ----------------------------------------------------------------------------------------------------------------------------------------------
# declared coverage (F01, F03, F05): nothing defaults to `not reported`, and no ledger row that names a pair goes unaccounted for
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _ids(spec):
    out = []
    for part in spec.split():
        if "-" in part:
            a, b = part.split("-")
            out += ["L%03d" % i for i in range(int(a[1:]), int(b[1:]) + 1)]
        else:
            out.append(part)
    return out


# Every pair x numeric column is a ledger leaf (NUMSPEC / EXT_SPEC / EXT_WALD_SPEC) or sits in exactly one reason group of NR_DECL[pair]. The union must be T2_NUMERIC; a column in
# neither fails the build (RAD23B's regional cells defaulted silently before this).
NR_CODES = ["no_check", "not_run", "no_row", "not_applicable", "not_printed"]
NR_REASON = {
    "no_check": "the pair had no pre-selected external outcome check (only the six pairs of Table 2 rows 1-6 had one), so the external-check and regional-colocalization blocks do not exist for it",
    "not_run": "the analysis was not run for the pair (ledger word row: not tested)",
    "no_row": "no ledger row stores this quantity for the pair; the build checks every ledger row whose key or claim names the pair (NOT_CARRIED lists the rows that exist but are not card columns)",
    "not_printed": "a ledger leaf stores this quantity for the pair, but the manuscript's current Table 2 does not print it, so the shipped card does not carry it (RC1b: the conservative-prior allowance names only values that Table 2 prints)",
    "not_applicable": "the quantity is defined for other pairs only (European p of the East Asian-discovered pair; multi-signal PP.H4 of the two kidney loci) or sits in another column for this pair "
                      "(an East Asian-discovered pair keeps its East Asian colocalization in pp_h4_*; the kidney loci lie outside the 81 map pairs and keep theirs in candidate_*)",
}
_CAND = ["candidate_susie_pp_h4", "candidate_coloc_abf_pp_h4"]
NR_DECL = {
    "PCSK9-CAD": {"no_row": ["reg_pp_h3_default", "eas_beta", "eas_pp_h1", "eas_pp_h4", "eas_pp_h4_default", "eas_pp_h4_conservative",
                             "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4"],
                  "not_applicable": ["eur_p"] + _CAND},
    "TGFB1-CAD": {"no_row": ["reg_pp_h3_default", "reg_pp_h3", "eas_pp_h1", "eas_pp_h3", "eas_pp_h4", "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4"],
                  "not_applicable": ["eur_p"] + _CAND},
    "PAM-T2D": {"no_row": ["reg_pp_h3_default", "reg_pp_h3", "eas_beta", "eas_pp_h3", "eas_pp_h4", "eas_pp_h4_default", "eas_pp_h4_conservative",
                           "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4"],
                "not_applicable": ["eur_p"] + _CAND},
    "RAD23B-T2D": {"no_row": ["eas_beta", "eas_pp_h1", "eas_pp_h3", "eas_pp_h4", "eas_pp_h4_default", "eas_pp_h4_conservative",
                              "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4"],
                   "not_applicable": ["eur_p"] + _CAND},
    "ARG1-T2D": {"no_row": ["reg_pp_h4_default", "reg_pp_h4_conservative", "reg_pp_h3_default", "reg_pp_h3", "eas_beta", "eas_pp_h1",
                            "eas_pp_h4_default", "eas_pp_h4_conservative"],
                 "not_applicable": ["eur_p"] + _CAND},
    "FGFR4-T2D (EAS)": {"no_row": ["reg_pp_h4_default", "reg_pp_h4_conservative", "reg_pp_h3_default", "reg_pp_h3", "eas_pp_h1", "eas_pp_h3",
                                   "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3"],
                        "not_applicable": ["eas_beta", "eas_pp_h4", "eas_pp_h4_default", "eas_pp_h4_conservative"] + _CAND},
    "UMOD-eGFR": {"no_check": ["ext_or", "ext_ci_lo", "ext_ci_hi", "ext_p", "ext_wald", "ext_wald_se", "ext_wald_p",
                               "reg_pp_h4_default", "reg_pp_h4_conservative", "reg_pp_h3_default", "reg_pp_h3"],
                  "no_row": ["eas_beta", "eas_pp_h1", "eas_pp_h3", "eas_pp_h4", "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4"],
                  "not_applicable": ["pp_h4_default", "pp_h4_conservative", "eur_p"],
                  "not_printed": ["eas_pp_h4_conservative"]},
    "DPEP1-eGFR": {"no_check": ["ext_or", "ext_ci_lo", "ext_ci_hi", "ext_p", "ext_wald", "ext_wald_se", "ext_wald_p",
                                "reg_pp_h4_default", "reg_pp_h4_conservative", "reg_pp_h3_default", "reg_pp_h3"],
                   "not_run": ["eas_beta", "eas_pp_h1", "eas_pp_h3", "eas_pp_h4", "eas_pp_h4_default", "eas_pp_h4_conservative",
                               "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4"],
                   "not_applicable": ["pp_h4_default", "pp_h4_conservative", "eur_p"]},
}

# Ledger rows whose key or claim names a card pair (gene and outcome tokens) and that no card column uses. Declared by id with a reason code; the build fails on a row that names a
# pair, is used by no column and is not listed here, and on a listed row that names no pair or is used. Row ids are the only literals.
NC_CODES = ["regional_variant", "duplicate", "signal_level", "ea_outcome_side", "screen_level", "literature_record", "analysis_inputs", "protocol_count"]
NC_REASON = {
    "regional_variant": "alternative handling of a carried regional-colocalization cell (outcome SNPs with a printed frequency of 0 or 1 dropped, marked descriptive; or beta/SE input, marked sensitivity); Table 3 and the supplement",
    "duplicate": "same stored value as a carried leaf (ext_wald_p); asserted equal at build",
    "signal_level": "signal-level or multi-signal test or estimate (conditional ratios, heterogeneity, GLS Q, Holm-adjusted p, credible-set size, attribution grid); Results and Table 3, not a per-pair card quantity",
    "ea_outcome_side": "East Asian outcome-side estimate of the European instrument; the card holds its state in EAS_state and the ancestry-matched colocalization in the eas_* columns",
    "screen_level": "screen- or post hoc-level statistic about the pair (screen counts, East Asian q and heterogeneity p, tied-variant floor, boundary-frequency heterogeneity p)",
    "literature_record": "literature record; the card holds Table 2's wording (literature_status) and the PMIDs",
    "analysis_inputs": "input counts (sample sizes, variants) of the signal analysis",
    "protocol_count": "protocol constant or cross-pair count that names the pair among others",
}
NOT_CARRIED = {
    "regional_variant": _ids("L735-L742"),
    "duplicate": _ids("L904-L908"),
    "signal_level": _ids("L731-L734 L902 L046 L054 L055 L057 L229 L265 L762-L777 L809-L813 L815 L821 L911 L912 L986 L875-L877 L879-L884 L903 L832-L835 L837-L840 L901"),
    "ea_outcome_side": ["L060"],
    "screen_level": _ids("L244 L428 L436 L437"),
    "literature_record": ["L061"],
    "analysis_inputs": ["L443"],
    "protocol_count": _ids("L951 L953 L964 L968"),
}
# the `duplicate` rows: the card pair whose ext_wald_p leaf each row repeats (the row's stored value is asserted equal to that leaf)
DUP_OF = {"L904": "ARG1-T2D", "L905": "FGFR4-T2D (EAS)", "L906": "PAM-T2D", "L907": "PCSK9-CAD", "L908": "TGFB1-CAD"}

# S2-D29 must-carry list (plan l.299): item -> the (pair, column) cells that may not read `not reported` in either card file (a withheld cell counts as carried)
MUST_CARRY = [
    ("discovery Wald ratio and p for PCSK9, TGFB1 and RAD23B", [(p, c) for p in ("PCSK9-CAD", "TGFB1-CAD", "RAD23B-T2D") for c in ("discovery_wald", "discovery_p")]),
    ("UMOD and DPEP1 eGFR* estimates", [("UMOD-eGFR", "discovery_wald"), ("DPEP1-eGFR", "discovery_wald")]),
    ("East Asian ancestry-matched results, PCSK9", [("PCSK9-CAD", "eas_pp_h3")]),
    ("East Asian ancestry-matched results, PAM", [("PAM-T2D", "eas_pp_h1"), ("PAM-T2D", "eas_ckbfree_pp_h1")]),
    ("East Asian ancestry-matched results, ARG1", [("ARG1-T2D", c) for c in ("eas_pp_h3", "eas_pp_h4", "eas_ckbfree_pp_h1", "eas_ckbfree_pp_h3", "eas_ckbfree_pp_h4")]),
    ("DIAMANTE estimate beside the discovery estimate (the former 'larger than discovery' note; the two numbers, no comparison sentence)",
     [("PAM-T2D", "ext_wald"), ("PAM-T2D", "discovery_wald")]),
    ("exposure-side arm results", [(p, c) for p in ("PCSK9-CAD", "TGFB1-CAD", "PAM-T2D", "RAD23B-T2D", "ARG1-T2D", "FGFR4-T2D (EAS)") for c in ("deCODE_arm", "ARIC_arm")]),
]
# the pegzilarginase staging note is the ARG1 ledger word statement itself: ARG1's development_statement must hold this label (and equal the ledger element)
MUST_CARRY_TEXT = [("pegzilarginase staging note", "ARG1-T2D", "development_statement", "pegzilarginase")]

# every element of an L370-L377 word row is held by a card column (substring test), except the literature element, which predates the step 1d Table 2 wording (RAD23B reads
# 'Not reported before' there against Table 2's hedged sentence) and is replaced by literature_status = Table 2's own sentence
ROLE_COLUMN = {"tier": "evidence_tier", "deCODE_arm": "deCODE_arm", "ARIC_arm": "ARIC_arm", "EAS_state": "EAS_state", "stage_dir": "development_statement",
               "grade": "row_group", "exposure_side": "deCODE_arm", "ext_check_note": "card_note", "instrument_note": "card_note", "lit": None}
LEAD_RX = re.compile(r"(?<=[\w\]])\.\s+(?=[A-Z])")      # end of a sentence of the established_role_reported_before cell
REPBEFORE_RX = re.compile(r"Reported before\b")          # start of the cell's statement of what was reported before (Stage 6.2: the cell opens with "Role: ...")

# Table 1 resource cell -> provider (declared; every row must match exactly one rule and every rule must be used)
PROVIDER_RULES = [(r"UKB-PPP", "UKB-PPP"), (r"^deCODE", "deCODE"), (r"\(ARIC\)", "ARIC"), (r"^FinnGen", "FinnGen"), (r"^MEGASTROKE", "MEGASTROKE"),
                  (r"^CKDGen", "CKDGen"), (r"^DIAMANTE", "DIAMANTE"), (r"^CARDIoGRAMplusC4D", "CARDIoGRAMplusC4D"),
                  (r"\(TPMI\)", "TPMI"), (r"^1000 Genomes", "1000Genomes")]   # tokens are space-free (F18): the provider column and the licence class carry the same token
EAS_VOCAB = [(r"\bCKB\b", "CKB"), (r"\bAGEN\b", "AGEN"), (r"\bBBJ\b|Biobank Japan", "BBJ")]
# Providers of the East Asian values of a pair whose ledger text names no East Asian resource (declared, with the file that names it). UMOD-eGFR: the ledger rows (L075) say only
# "East Asian colocalization"; the pilot census row UMOD,eGFR (V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv, eas_outcome_resource
# 'TWB2+BBJ eGFR meta-analysis', accession TWB2_BBJ_eGFR_PMID39472450; instrument from CKB) names TWB, BBJ and CKB, the same class as release input RI057 (phase2_egfr_scan.tsv).
EAS_RESOLVED = {"UMOD-eGFR": ("BBJ", "CKB", "TWB2-BBJ")}

CITE_RX = re.compile(r"\[@([^\]]+)\]")
HOLD_RX = re.compile(r"\[HOLD ([^:\]]+):")
PRINT_RX = re.compile(r"([−+\-]?\d[\d,]*(?:\.\d+)?(?:\s?×\s?10[⁻⁰¹²³⁴-⁹]+)?)\s?\[\[(L\d{3}[a-z]?)\]\]")
BONF_RX = re.compile(r"concordant, (not )?Bonferroni")
SRC_RX = re.compile(r"^(?P<col>\w+)=(?P<rid>L\d{3}[a-z]?)(?:#(?P<flat>\d+)|@(?P<elem>\d+)|\[(?P<pair>[^\]]+)\](?:\.(?P<j>\d+))?)$")
NUMTOK_RX = re.compile(r"^([+-]?)(\d[\d,]*(?:\.\d+)?)(?:\s?×\s?10\^(-?\d+)|[eE]([+-]?\d+))?$")
# the TPMI mentions that are labels or provenance and carry no TPMI value (column -> exact substrings that may remain in a shipped cell)
TPMI_ALLOWED = {"ext_resource": ["TPMI"], "licence_class": ["third_party_derived:TPMI"], "resource": ["Taiwan Precision Medicine Initiative (TPMI) PheWeb", "new:TPMI2025"],
                "provider": ["TPMI"], "cited_keys": ["new:TPMI2025"]}
NOSCAN_COLS = {"source_md5", "ledger_row_ids", "numeric_cell_sources", "search_scope", "literature_pmids", "cited_keys", "cited_pmids", "not_reported_reasons",
               "ledger_rows_not_carried"}
# S2-D30(3): the four Table 2 prose columns ship verbatim (equal to the manuscript's Table 2 and to the withheld card), so the TPMI and Open Targets scans exempt exactly these.
# Table1_data_sources.tsv mirrors Table 1 and keeps its TPMI case and control counts for the same reason; its TPMI mentions are the labels in TPMI_ALLOWED.
PROSE_EXEMPT = set(T2_PROSE) | set(T2_FULL)


class RelTablesError(AssertionError):
    def __init__(self, check, msg):
        AssertionError.__init__(self, "[%s] %s" % (check, msg))
        self.check = check


def _req(cond, check, msg):
    if not cond:
        raise RelTablesError(check, msg)


# definition of every column of both table files (written beside the data; the verifier requires it to list exactly T2_COLS then T1_COLS, in order)
T2F, T1F = "Table2_evidence_cards.tsv", "Table1_data_sources.tsv"
DEF2 = [
    ("pair_key", "Pair as keyed in the ledger and the outcome-check file (hyphen; '(EAS)' marks the East Asian pair)"),
    ("gene", "Protein gene symbol"),
    ("outcome", "Outcome: CAD, T2D or eGFR"),
    ("row_group", "pre-selected external-check pair; East Asian, not graded; multi-signal candidate, not graded"),
    ("table_row_label", "First cell of the Table 2 row, ledger tags and bold markers removed"),
    ("evidence_tier", "Evidence tier from the ledger word row; empty for the three ungraded rows"),
    ("deCODE_arm", "Exposure-side corroboration in deCODE (ledger word row); FGFR4: 'exposure-side not tested'; kidney loci: not reported"),
    ("ARIC_arm", "Exposure-side corroboration in ARIC (ledger word row); FGFR4: 'exposure-side not tested'; kidney loci: not reported"),
    ("EAS_state", "East Asian outcome-side state of the European instrument (graded pairs) or East Asian testing state (kidney loci), from the ledger word row; not the ancestry-matched colocalization of the eas_* columns"),
    ("card_note", "Ledger word statements that no other column holds: the FGFR4 instrument note and TPMI note (the TPMI note reads withheld_pending_terms when the TPMI switch is 'withheld'); not reported for the other rows"),
    ("discovery_basis", "Ledger claim of the row that feeds discovery_wald (names the scale and the resource)"),
    ("discovery_wald", "Discovery Wald ratio, log-odds per SD of protein; kidney loci on the eGFR* scale (-ln eGFR, positive = lower eGFR); FGFR4: East Asian (CKB instrument), as pp_h4_default and pp_h4_conservative are"),
    ("discovery_p", "Discovery Mendelian-randomization p for the same estimate"),
    ("pp_h4_default", "PP.H4 of the pair's discovery colocalization at the default prior p12 = 1e-5 (the map colocalization; FGFR4: the East Asian colocalization of its CKB instrument with AGEN, because its discovery is East Asian); "
                      "not reported where the ledger holds no pair-level value at that prior (kidney loci: see candidate_*)"),
    ("pp_h4_conservative", "PP.H4 at the conservative prior p12 = 1e-6 (same scope as pp_h4_default)"),
    ("ext_resource", "External outcome resource of the pre-selected check"),
    ("ext_or", "Odds ratio per SD of protein from the first-order Wald ratio in the external resource (ledger L223)"),
    ("ext_ci_lo", "Lower 95% limit of ext_or"),
    ("ext_ci_hi", "Upper 95% limit of ext_or"),
    ("ext_p", "The external resource's own reported p (ledger L960); the first-order Wald p of the same estimate is ext_wald_p"),
    ("ext_wald", "First-order Wald ratio in the external resource, log-odds per SD of protein (numeric leaf 0 of the pair's external-check ledger row). The same row stores the discovery estimate it is set beside "
                 "(leaf 3, asserted equal to discovery_wald at the printed precision); no comparison sentence is written here"),
    ("ext_wald_se", "Standard error of ext_wald (leaf 1 of the same ledger row)"),
    ("ext_wald_p", "First-order Wald p of ext_wald (leaf 2 of the same ledger row)"),
    ("ext_class", "Classification of the external check (ledger word of the pair's external-check row)"),
    ("ext_bonferroni_6", "true/false: p below 0.05/6 (Bonferroni over the six pre-selected pairs), from the outcome-check file, asserted against Table 2 and the ledger"),
    ("reg_pp_h4_default", "Regional colocalization of the protein with the external outcome data (primary row): PP.H4 at p12 = 1e-5"),
    ("reg_pp_h4_conservative", "Regional colocalization (primary row): PP.H4 at p12 = 1e-6"),
    ("reg_pp_h3_default", "Regional colocalization (primary row): PP.H3 at p12 = 1e-5; not reported where the ledger stores no primary row at that prior"),
    ("reg_pp_h3", "Regional colocalization (primary row): PP.H3 at the conservative prior p12 = 1e-6; not reported where the ledger stores no primary row at that prior"),
    ("eas_basis", "Verbatim ledger claims of the rows that feed the eas_* columns and eur_p"),
    ("eas_beta", "East Asian log-odds per SD (TGFB1; CKB instrument)"),
    ("eas_pp_h1", "Ancestry-matched East Asian colocalization PP.H1 (resource named in eas_basis); European-discovery pairs only"),
    ("eas_pp_h3", "Ancestry-matched East Asian colocalization PP.H3 (resource named in eas_basis); European-discovery pairs only"),
    ("eas_pp_h4", "Ancestry-matched East Asian colocalization PP.H4, prior not stated in the ledger row (ARG1, AGEN); European-discovery pairs only"),
    ("eas_pp_h4_default", "Ancestry-matched East Asian PP.H4 at p12 = 1e-5; European-discovery pairs only (FGFR4 is East Asian-discovered: its AGEN colocalization is pp_h4_default)"),
    ("eas_pp_h4_conservative", "Ancestry-matched East Asian PP.H4 at p12 = 1e-6; European-discovery pairs only (FGFR4: see pp_h4_conservative)"),
    ("eas_ckbfree_pp_h1", "East Asian CKB-free (Biobank Japan) PP.H1"),
    ("eas_ckbfree_pp_h3", "East Asian CKB-free (Biobank Japan) PP.H3"),
    ("eas_ckbfree_pp_h4", "East Asian CKB-free (Biobank Japan) PP.H4"),
    ("eur_p", "European p (FGFR4 only)"),
    ("candidate_susie_pp_h4", "Multi-signal candidates (the two kidney loci): SuSiE multi-signal PP.H4 (purity 0.5); the ledger row does not state a prior, so it fills no prior-specific column"),
    ("candidate_coloc_abf_pp_h4", "Multi-signal candidates (the two kidney loci): coloc.abf single-variant PP.H4; the ledger row does not state a prior"),
    ("development_stage", "Open Targets development stage (single-protein rule of the ledger builder): 'tractable' is an Open Targets tractability class and ships; 'approved' and 'in clinical development' are set by a drug record (ChEMBL-derived) and read withheld_pending_terms in the shipped file"),
    ("development_agents", "Single-protein agents with their highest stage (ChEMBL-derived drug records); 'none' ships, a list of agents reads withheld_pending_terms in the shipped file"),
    ("genetic_direction_agreement", "Whether the agent's effect agrees with the genetic direction (ledger word row); the verdict ships when it names no drug, and reads withheld_pending_terms when it does"),
    ("development_statement", "The ledger word statement of stage, agents and direction, verbatim (stage_dir element of the L370-L377 row), including staging notes such as an enzyme that does not set the stage; ships for a row with no clinical candidate whose statement names no drug, otherwise reads withheld_pending_terms in the shipped file"),
    ("ot_data_version", "Open Targets Platform data version of the snapshot"),
    ("ot_query_date", "Date the snapshot was queried"),
    ("literature_status", "Table 2's own wording of the prior-evidence status: the sentence of the row's established_role_reported_before cell that begins 'Reported before' (a literal substring of it, asserted). The literature search log's status_in_table column "
                          "predates the Stage 5 step 1d wording and can read stronger than Table 2; it is not repeated here"),
    ("literature_pmids", "PMIDs of the references cited in the row's prose cells (reference tables); a key with no PMID is listed in cited_keys only"),
    ("search_date", "Date(s) of the searches in the log"),
    ("search_scope", "Database, query, hits and screened count of each log row of the pair"),
    ("cited_keys", "[@key] reference keys cited in the row's prose cells, in order of first use"),
    ("established_role_reported_before", "Table 2 column 'Established role; reported before', ledger tags removed, [@key] citations kept inline; identical to the manuscript's Table 2 and to the withheld full card (S2-D30(3))"),
    ("adds", "Table 2 cell, ledger tags removed, [@key] citations kept inline; identical to the manuscript's Table 2 and to the withheld full card (S2-D30(3))"),
    ("main_uncertainty", "Table 2 cell, ledger tags removed, [@key] citations kept inline; identical to the manuscript's Table 2 and to the withheld full card (S2-D30(3))"),
    ("next_biological_question", "Table 2 column 'Next biological question', ledger tags removed, [@key] citations kept inline; identical to the manuscript's Table 2 and to the withheld full card (S2-D30(3))"),
    ("full_established_role_reported_before", "the full card's cell in Additional file 2, section {SM:cards} (supplement/P31R_SM_cards.md), ledger tags removed, [@key] citations kept inline; Table 2 prints a summary of it (S2-D51)"),
    ("full_adds", "the full card's cell in Additional file 2, section {SM:cards} (supplement/P31R_SM_cards.md), ledger tags removed, [@key] citations kept inline; Table 2 prints a summary of it (S2-D51)"),
    ("full_main_uncertainty", "the full card's cell in Additional file 2, section {SM:cards} (supplement/P31R_SM_cards.md), ledger tags removed, [@key] citations kept inline; Table 2 prints a summary of it (S2-D51)"),
    ("full_next_biological_question", "the full card's cell in Additional file 2, section {SM:cards} (supplement/P31R_SM_cards.md), ledger tags removed, [@key] citations kept inline; Table 2 prints a summary of it (S2-D51)"),
    ("numeric_cell_sources", "For every numeric cell, the ledger row and stored leaf it reads (Lxxx#leaf, Lxxx[pair].element, Lxxx@element)"),
    ("not_reported_reasons", "Every numeric column of the row that reads 'not reported', grouped by reason code (code: columns): " + "; ".join("%s = %s" % (c, NR_REASON[c]) for c in NR_CODES)),
    ("ledger_rows_not_carried", "Ledger rows whose key or claim names the pair (gene and outcome) that no card column uses, grouped by reason code (code: ledger ids): " + "; ".join("%s = %s" % (c, NC_REASON[c]) for c in NC_CODES)),
    ("ledger_row_ids", "All ledger rows the card uses: the row's tags in Table 2, the cell sources and the word row"),
    ("licence_class", "own_derived plus a third_party_derived:<provider> class for each provider the row's values derive from (declared rules; the tokens contain no space)"),
    ("source_md5", "md5 of the input files the row was generated from"),
]
DEF1 = [
    ("provider", "Provider of the row's resource (declared resource-to-provider rules; a token without a space)"),
    ("licence_class", "third_party_derived:<provider> of the row's resource. Rule: every Table 1 row is a third-party resource, so the class names its provider. This file mirrors the manuscript's Table 1 verbatim and "
                      "is not masked: a TPMI-derived value that Table 1 prints (the TPMI row's case and control counts) ships under third_party_derived:TPMI because the manuscript prints it (S2-D30(3)); "
                      "the derived TPMI estimates of the external check ship in Table2_evidence_cards.tsv when the TPMI switch is 'ship' (S2-D34(1))"),
    ("cited_keys", "[@key] reference keys cited in the row, in order of first use"),
    ("cited_pmids", "PMIDs of the cited references (reference tables); a key with no PMID is listed in cited_keys only"),
    ("ledger_row_ids", "Ledger rows tagged in the row"),
]


def _column_doc(ctx=None):
    """[(file, column, definition)] for T2_COLS then T1_COLS; the lists must cover the file columns exactly. Keys without a PMID are appended from the build (F17: derived, not typed)."""
    _req([c for c, _ in DEF2] == T2_COLS, "FORMAT", "DEF2 does not list exactly the Table 2 columns, in order")
    d1 = [(c, "Table 1 column of the same name, ledger tags and bold markers removed, [@key] citations kept inline (verbatim mirror of the manuscript's Table 1)") for c in T1_NAMES] + DEF1
    _req([c for c, _ in d1] == T1_COLS, "FORMAT", "DEF1 does not list exactly the Table 1 columns, in order")
    nop = _nopmid(ctx) if ctx is not None else []
    note = ("; keys with no PMID in the reference tables at this build (Stage 6 reference check): " + ", ".join(nop)) if nop else ""
    out = []
    for f, rows in ((T2F, DEF2), (T1F, d1)):
        for c, d in rows:
            out.append((f, c, d + (note if c in ("literature_pmids", "cited_pmids") else "")))
    return out


# ----------------------------------------------------------------------------------------------------------------------------------------------
# dependencies: the manuscript builder's own tag handling and the gate's own number tokenizer, imported (no side effects: build dir redirected)
# ----------------------------------------------------------------------------------------------------------------------------------------------
_DEPS = {}


def _deps(v3):
    key = str(pathlib.Path(v3).resolve())
    if key in _DEPS:
        return _DEPS[key]
    scripts = str(pathlib.Path(v3) / "R2_rewrite" / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    tmp = tempfile.mkdtemp(prefix="rel_tables_mb_")
    old = os.environ.get("P31R_BUILD_OUT")
    os.environ["P31R_BUILD_OUT"] = tmp          # P31R_build makes its output folder at import; point it at a scratch folder
    old_dwb, sys.dont_write_bytecode = sys.dont_write_bytecode, True     # no new .pyc in the project's scripts folder
    try:
        import P31R_gate as PG
        import P31R_tablegate as TG
        import P31R_build as MB
    finally:
        sys.dont_write_bytecode = old_dwb
        if old is None:
            os.environ.pop("P31R_BUILD_OUT", None)
        else:
            os.environ["P31R_BUILD_OUT"] = old
        shutil.rmtree(tmp, ignore_errors=True)
    for obj in (PG, TG, MB):     # the modules must come from this tree: one v3 per process
        _req(str(pathlib.Path(obj.__file__).resolve()).startswith(str(pathlib.Path(v3).resolve())), "DEPS",
             "%s was imported from %s, not from %s" % (obj.__name__, obj.__file__, v3))
    for obj, names in ((PG, ["_d1_tokens", "_d1_leaf_list", "_d1_match1", "_d1_hit_at", "d1_order", "TAG"]), (TG, ["parse_table"]),
                       (MB, ["TAG", "CIT", "NEW_PMID", "ALIAS"])):
        for n in names:
            _req(hasattr(obj, n), "DEPS", "%s.%s is missing: the pinned gate/builder changed shape" % (obj.__name__, n))
    _DEPS[key] = (PG, TG, MB)
    return _DEPS[key]


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read(path):
    with io.open(path, encoding="utf-8", newline="") as fh:
        return fh.read()


def _id_key(rid):
    return (int(re.sub(r"\D", "", rid)), rid)


def norm_num(text):
    """Printed token of a ledger display -> machine form ('3.9 x 10^-12' -> '3.9e-12', '1,234' -> '1234'); the digits are kept as printed."""
    s = text.strip().replace("−", "-")
    m = NUMTOK_RX.match(s)
    _req(m, "NUMCELL", "unparsable printed token %r" % text)
    sign, mant, e1, e2 = m.groups()
    exp = e1 if e1 is not None else e2
    return sign + mant.replace(",", "") + ("e%d" % int(exp) if exp is not None else "")


def cited_keys(text):
    out = []
    for grp in CITE_RX.findall(text):
        for k in grp.split(";"):
            k = k.strip().lstrip("@").strip()
            if k and k not in out:
                out.append(k)
    return out


def _tsv_write(path, cols, rows):
    """Plain TSV: cells are written as they are, with no CSV quoting. A cell may hold a double quote, but no tab, CR or LF, and may not begin with a double quote (a quote-aware reader
    would take it for a quoted field); both are asserted (F18)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(cols)]
    for r in rows:
        cells = []
        for c in cols:
            v = r[c]
            _req(isinstance(v, str) and not re.search(r"[\t\r\n]", v), "FORMAT", "cell %s is not a single-line string without a tab: %r" % (c, v))
            _req(not v.startswith('"'), "FORMAT", "cell %s begins with a double quote: %r" % (c, v[:40]))
            cells.append(v)
        lines.append("\t".join(cells))
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\n".join(lines) + "\n")


def _tsv_read(path):
    text = _read(path)
    _req("\r" not in text and text.endswith("\n"), "FORMAT", "%s: not LF-terminated plain text" % path)
    rows = [ln.split("\t") for ln in text[:-1].split("\n")]
    head, body = rows[0], rows[1:]
    for r in body:
        _req(len(r) == len(head), "FORMAT", "%s: ragged row" % path)
        _req(not any(c.startswith('"') for c in r), "FORMAT", "%s: a cell begins with a double quote" % path)
    return head, [dict(zip(head, r)) for r in body]


# ----------------------------------------------------------------------------------------------------------------------------------------------
# Open Targets helpers: the ledger builder's own functions, extracted from ext_table2.py by AST and run on the saved snapshot (reuse, not re-implement)
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _ext_table2_ns(src, snapshot):
    tree = ast.parse(src)
    want_assign, want_func = {"_LITCELL", "_DS", "_STAGE"}, {"_rows_of", "_single", "_tract", "_stage"}
    nodes, seen = [], set()
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in want_assign:
            nodes.append(node); seen.add(node.targets[0].id)
        elif isinstance(node, ast.FunctionDef) and node.name in want_func:
            nodes.append(node); seen.add(node.name)
    _req(seen == want_assign | want_func, "DEPS", "ext_table2.py no longer defines %s" % sorted((want_assign | want_func) - seen))
    ns = {"_sn": snapshot}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), EXT2_REL, "exec"), ns)
    return ns


# ----------------------------------------------------------------------------------------------------------------------------------------------
# context: everything read from disk / passed in, with the cross-checks that need no output
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _parse_cards(text):
    """Stage 7 A2.4c: the full cards of supplement/P31R_SM_cards.md -> {row label: [four raw cells]}. A card is its label line (the Table 2 row label,
    verbatim) followed by the four labelled paragraphs in CARD_PREFIX order."""
    lines = [x for x in text.replace("\r\n", "\n").split("\n") if x.strip()]
    cards, i = {}, 0
    while i < len(lines):
        if lines[i].startswith("**") and i + 4 < len(lines) and lines[i + 1].startswith(CARD_PREFIX[0]):
            cells = [lines[i + 1 + j] for j in range(4)]
            _req(all(c.startswith(p) for c, p in zip(cells, CARD_PREFIX)), "CARDS", "card %r: its paragraphs are not in CARD_PREFIX order" % lines[i])
            _req(lines[i] not in cards, "CARDS", "card label twice: %r" % lines[i])
            cards[lines[i]] = [c[len(p):] for c, p in zip(cells, CARD_PREFIX)]
            i += 5
        else:
            i += 1
    _req(len(cards) >= 1, "CARDS", "no full card parsed from %s" % CARDS_REL)
    return cards


def _make_context(v3, ledger, t1_text=None, t2_text=None, results_text=None):
    v3 = pathlib.Path(v3)
    PG, TG, MB = _deps(v3)
    df = ledger.fillna("") if hasattr(ledger, "fillna") else ledger
    L = {}
    for r in df.to_dict("records"):
        _req(r["id"] not in L, "LEDGER", "duplicate ledger id %s" % r["id"])
        L[r["id"]] = r
    ctx = {"v3": v3, "PG": PG, "TG": TG, "MB": MB, "L": L, "md5": {}}
    for rel in INPUTS:
        _req((v3 / rel).is_file(), "INPUT", "declared input missing: %s" % rel)
        ctx["md5"][rel] = _md5(v3 / rel)
    ctx["t1_text"] = t1_text if t1_text is not None else _read(v3 / T1_REL)
    ctx["t2_text"] = t2_text if t2_text is not None else _read(v3 / T2_REL)
    ctx["res_text"] = results_text if results_text is not None else _read(v3 / RES_REL)
    ctx["t1"] = TG.parse_table(ctx["t1_text"])
    ctx["t2"] = TG.parse_table(ctx["t2_text"])
    ctx["cards"] = _parse_cards(_read(v3 / CARDS_REL))
    ctx["log"] = list(csv.DictReader(io.open(v3 / LOG_REL, encoding="utf-8", newline=""), delimiter="\t"))
    ctx["oc"] = list(csv.DictReader(io.open(v3 / OC_REL, encoding="utf-8", newline=""), delimiter="\t"))
    ctx["snap"] = json.load(io.open(v3 / SNAP_REL, encoding="utf-8"))
    ctx["ext2_src"] = _read(v3 / EXT2_REL)
    ctx["ns"] = _ext_table2_ns(ctx["ext2_src"], ctx["snap"])
    ctx["oldpmid"] = {r["n"]: r["pmid"] for r in csv.DictReader(io.open(v3 / OLDREFS_REL, encoding="utf-8", newline=""), delimiter="\t")}
    ctx["keytab"] = {r["key"]: r["pmid"] for r in csv.DictReader(io.open(v3 / KEYTAB_REL, encoding="utf-8", newline=""), delimiter="\t")}
    return ctx


# ----------------------------------------------------------------------------------------------------------------------------------------------
# leaves and printed tokens
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _flat_leaves(ctx, rid):
    return ctx["PG"]._d1_leaf_list(ctx["L"][rid]["value"])


def _l223_flat(ctx, pair, j):
    obj = json.loads(ctx["L"]["L223"]["value"])
    idx = 0
    for row in obj:
        nums = [i for i, x in enumerate(row) if isinstance(x, (int, float)) and not isinstance(x, bool)]
        if row[0] == pair:
            _req(j in nums, "NUMCELL", "L223 element %d of %s is not numeric" % (j, pair))
            return idx + nums.index(j), row[j]
        idx += len(nums)
    raise RelTablesError("NUMCELL", "L223 holds no row for %s" % pair)


def _oc_pairs(ctx):
    """Pair names of the first six outcome-check rows: the leaf order of L960 (its source column says 'first six rows')."""
    return [r["pair"] for r in ctx["oc"][:len(EXT_PAIRS)]]


def _l960_flat(ctx, pair):
    order = _oc_pairs(ctx)
    _req(pair in order, "NUMCELL", "L960: %s is not among the first six outcome-check rows %s" % (pair, order))
    i = order.index(pair)
    row = [r for r in ctx["oc"] if r["pair"] == pair][0]
    leaves, _ = _flat_leaves(ctx, "L960")
    want, got = float(row["p_out_file"]), leaves[i]
    _req(abs(want - got) <= 1e-9 * max(abs(want), abs(got)), "NUMCELL", "L960 leaf %d (%s) differs from the outcome-check file p_out_file %r" % (i, pair, row["p_out_file"]))
    return i, got


def _resolve_leaf(ctx, rid, sel, pair):
    """(flat leaf index, leaf value, source ref string suffix) for a NUMSPEC selector."""
    if isinstance(sel, int):
        leaves, _ = _flat_leaves(ctx, rid)
        _req(sel < len(leaves), "NUMCELL", "%s has %d leaves, wanted #%d" % (rid, len(leaves), sel))
        return sel, leaves[sel], "#%d" % sel
    if sel[0] == "n":
        _req(rid == "L223", "NUMCELL", "nested selector only for L223")
        flat, val = _l223_flat(ctx, pair, sel[1])
        return flat, val, "[%s].%d" % (pair, sel[1])
    if sel[0] == "p":
        _req(rid == "L960", "NUMCELL", "pair selector only for L960")
        flat, val = _l960_flat(ctx, pair)
        return flat, val, "[%s]" % pair
    raise RelTablesError("NUMCELL", "bad selector %r" % (sel,))


def _printed_token(ctx, rid, leaf_idx):
    """The ledger display's own printed token for stored leaf #leaf_idx, using the gate's tokenizer and its reading order (declared order if any)."""
    PG, row = ctx["PG"], ctx["L"][rid]
    leaves, ordered = PG._d1_leaf_list(row["value"])
    toks = [t for t in PG._d1_tokens(row["display"]) if t.kind == "num"]
    decl = PG.d1_order().get(rid) if ordered else None
    _req(ordered, "NUMCELL", "%s stores a JSON object: leaf order is not display order" % rid)
    if decl is not None:
        _req(len(decl) == len(toks), "NUMCELL", "%s: declared display order has %d entries for %d tokens" % (rid, len(decl), len(toks)))
    ptr = 0
    for k, t in enumerate(toks):
        if decl is not None:
            li = decl[k]
        else:
            h = PG._d1_hit_at(t, leaves, ptr)
            _req(h is not None, "NUMCELL", "%s: display token %r matches no stored leaf at or after #%d" % (rid, t.text, ptr))
            li = h[0]
            ptr = li + 1
        if li == leaf_idx:
            _req(PG._d1_match1(t, leaves[leaf_idx]) is not None, "NUMCELL", "%s: token %r does not round from leaf #%d" % (rid, t.text, leaf_idx))
            return t.text
    raise RelTablesError("NUMCELL", "%s: its display %r prints no token for stored leaf #%d" % (rid, row["display"][:70], leaf_idx))


def _cell_matches_leaf(ctx, cell, leaf):
    PG = ctx["PG"]
    toks = [t for t in PG._d1_tokens(cell) if t.kind == "num"]
    return len(toks) == 1 and PG._d1_match1(toks[0], leaf) is not None


# ----------------------------------------------------------------------------------------------------------------------------------------------
# Table 1
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _strip(ctx, s):
    return ctx["MB"].TAG.sub("", s).replace("**", "")


def _strip_cites(s):
    return re.sub(r"\s?\[@[^\]]+\]", "", s)


def _pmids_of(ctx, keys):
    """PMIDs of the cited keys. Sources: the program reference table P31R_REFERENCE_KEYS.tsv (the R2 table), the manuscript builder's own table
    (P31R_build.NEW_PMID, which merges refs/new_refs_*.tsv) and the old numbered list; where two give a PMID they must agree. A key unknown to all fails."""
    MB, old, tab = ctx["MB"], ctx["oldpmid"], ctx["keytab"]
    out = []
    for key in keys:
        cand = {}
        if key in tab and tab[key]:
            cand["P31R_REFERENCE_KEYS.tsv"] = tab[key]
        if key in MB.NEW_PMID:
            cand["builder NEW_PMID"] = MB.NEW_PMID[key]
        if key.startswith("old:") and old.get(key[4:]):
            cand["_refs_resolved.tsv"] = old[key[4:]]
        known = key in tab or key in MB.NEW_PMID or (key.startswith("old:") and key[4:] in old) or key in MB.ALIAS
        _req(known, "CITE", "reference key %s is in no reference table" % key)
        _req(len(set(cand.values())) <= 1, "CITE", "reference key %s has conflicting PMIDs %s" % (key, cand))
        p = next(iter(cand.values()), "")
        if p and p not in out:
            out.append(p)
    return out


def _build_table1(ctx):
    n, title, header, body, notes = ctx["t1"]
    _req(len(header) == len(T1_NAMES) and header[0] == "Role" and header[1] == "Resource", "ROWID1", "Table 1 header changed shape: %r" % header[:3])
    used, rows = set(), []
    for cells in body:
        _req(len(cells) == len(header), "ROWID1", "Table 1 ragged row")
        raw_tags = [t for c in cells for t in ctx["PG"].TAG.findall(c)]
        _req(all(t in ctx["L"] for t in raw_tags), "TAGS", "Table 1 tags a ledger id that does not exist")
        r = {name: _strip(ctx, c) for name, c in zip(T1_NAMES, cells)}
        res = _strip_cites(r["resource"])
        hits = [(rx, p) for rx, p in PROVIDER_RULES if re.search(rx, res)]
        _req(len(hits) == 1, "PROVIDER", "Table 1 resource %r matches %d provider rules" % (res, len(hits)))
        used.add(hits[0][0])
        r["provider"] = hits[0][1]
        r["licence_class"] = "third_party_derived:" + hits[0][1]
        _req(not re.search(r"\s", r["provider"]), "LICENCE", "provider token %r contains whitespace" % r["provider"])
        keys = cited_keys(" ".join(cells))
        r["cited_keys"] = "; ".join(keys) if keys else NOT_REPORTED
        pm = _pmids_of(ctx, keys)
        r["cited_pmids"] = "; ".join(pm) if pm else NOT_REPORTED
        ids = sorted({t for c in cells for t in ctx["PG"].TAG.findall(c)}, key=_id_key)
        r["ledger_row_ids"] = "; ".join(ids) if ids else NOT_REPORTED
        rows.append(r)
    _req(used == {rx for rx, _ in PROVIDER_RULES}, "PROVIDER", "unused provider rules: %s" % sorted({rx for rx, _ in PROVIDER_RULES} - used))
    return rows


# ----------------------------------------------------------------------------------------------------------------------------------------------
# Table 2
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _parse_label(ctx, raw):
    """First cell -> (pair_key, row_group_hint, label). '**PCSK9–CAD**; Convergent' / '**FGFR4–T2D**, East Asian; not graded' (builder's tag handling)."""
    label = _strip(ctx, raw)
    m = re.match(r"^([A-Za-z0-9]+)–([A-Za-z0-9]+)(, East Asian)?(?:, [^;]*)?;\s*(.+)$", label)
    _req(m, "ROWID2", "Table 2 first cell has an unexpected shape: %r" % label)
    gene, outcome, eas, rest = m.group(1), m.group(2), m.group(3), m.group(4)
    return gene + "-" + outcome + (" (EAS)" if eas else ""), rest, label


def _ledger_words(ctx, gene, outcome):
    rows = [r for r in ctx["L"].values() if r["section"] == "Table 2 (evidence cards)" and "word statements" in r["claim"]
            and r["claim"].startswith("%s-%s " % (gene, outcome))]
    _req(len(rows) == 1, "WORDS", "%s-%s: %d word rows in the ledger" % (gene, outcome, len(rows)))
    return rows[0]["id"], json.loads(rows[0]["value"])


def _split_top(s, sep=";"):
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == sep and depth == 0:
            out.append(cur.strip()); cur = ""
        else:
            cur += ch
    out.append(cur.strip())
    return out


def _ot_columns(ctx, gene, stage_dir_text):
    ns = ctx["ns"]
    stage = ns["_stage"](gene)
    _req(ns["_DS"][gene] == stage, "OTSTAGE", "%s: snapshot stage %r differs from the ledger builder's declared %r" % (gene, stage, ns["_DS"][gene]))
    first = _split_top(stage_dir_text)[0].lower()
    _req(first.startswith(stage), "OTSTAGE", "%s: ledger words %r do not start with the snapshot stage %r" % (gene, stage_dir_text[:60], stage))
    single = ns["_single"](gene)

    def drug_stage(name):      # the drug record's own maximum stage (the ledger builder asserts approval from it, e.g. tafolecimab)
        for key in (name.lower(), name, name.upper()):
            d = (ctx["snap"].get("drugs") or {}).get(key)
            if d:
                try:
                    return ns["_STAGE"].get(d["data"]["drug"]["maximumClinicalStage"], 0)
                except (KeyError, TypeError):
                    return 0
        return 0
    agents = "; ".join("%s (stage %d)" % (k, max(single[k], drug_stage(k))) for k in sorted(single)) if single else "none"
    dirs = [c for c in _split_top(stage_dir_text)[1:] if re.search(r"\b(agrees?|opposes?)\b", c)]
    return stage, agents, "; ".join(dirs) if dirs else NOT_REPORTED


def _literature(ctx, gene, outcome, words_lit):
    pair = "%s-%s" % (gene, outcome)
    rows = [r for r in ctx["log"] if r["pair"] == pair]
    _req(rows, "LITSTAT", "%s: no row in the literature search log" % pair)
    sts = []
    for r in rows:
        if r["status_in_table"] not in sts:
            sts.append(r["status_in_table"])
    final = [s for s in sts if not re.search(r"\(see .* row\)$", s)]
    _req(len(final) == 1, "LITSTAT", "%s: the log does not settle on one status %r" % (pair, sts))
    _req(words_lit.startswith(final[0]), "LITSTAT", "%s: ledger words %r do not start with the log status %r" % (pair, words_lit[:50], final[0]))
    dates = []
    for r in rows:
        if r["date"] not in dates:
            dates.append(r["date"])
    scope = " || ".join("%s: %s (hits %s, screened %s)" % (r["database"], r["query"], r["hits"], r["screened"]) for r in rows)
    return final[0], "; ".join(dates), scope


def _build_table2_full(ctx, t1_rows):
    PG, MB, L = ctx["PG"], ctx["MB"], ctx["L"]
    n, title, header, body, notes = ctx["t2"]
    _req(len(header) == 5 and len(body) == len(NUMSPEC), "ROWID2", "Table 2 has %d columns and %d rows; the spec is 5 and %d" % (len(header), len(body), len(NUMSPEC)))
    _req([MB.TAG.sub("", h).replace("**", "").strip() for h in header] == T2_HEADER, "COLS2",
         "Table 2 header changed: %r (the cards follow %r)" % ([h.strip() for h in header], T2_HEADER))
    oc = {r["pair"]: r for r in ctx["oc"]}
    cut = json.loads(L["L961"]["value"])[2]
    log_md5 = ctx["md5"][LOG_REL]
    rows, rules_note = [], []
    for cells in body:
        pair_key, rest, label = _parse_label(ctx, cells[0])
        _req(pair_key in NUMSPEC, "ROWID2", "Table 2 row %r is not in the declared pair list" % pair_key)
        gene, tail = pair_key.split("-", 1)
        outcome = tail.split(" ")[0]
        tier_word = next((t for t in TIERS if rest.startswith(t)), "")
        group = ("pre-selected external-check pair" if tier_word else "East Asian, not graded" if "(EAS)" in pair_key else
                 "multi-signal candidate, not graded" if "multi-signal candidate" in rest else None)
        _req(group, "ROWID2", "%s: no row group for label %r" % (pair_key, label))
        _req(bool(tier_word) == ("not graded" not in rest), "TIER", "%s: label %r disagrees on whether the row is graded" % (pair_key, label))
        wid, words = _ledger_words(ctx, gene, outcome)
        roles = ROLES[group]
        _req(len(words) == len(roles), "WORDS", "%s: %s holds %d elements, expected %d" % (pair_key, wid, len(words), len(roles)))
        wd = dict(zip(roles, words))
        if tier_word:
            _req(wd["tier"] == tier_word, "TIER", "%s: label tier %r vs ledger %s tier %r" % (pair_key, tier_word, wid, wd["tier"]))
        r = {"pair_key": pair_key, "gene": gene, "outcome": outcome, "row_group": group, "table_row_label": label, "evidence_tier": wd.get("tier", "")}
        if group == "East Asian, not graded":
            r["deCODE_arm"] = r["ARIC_arm"] = wd["exposure_side"]
            r["EAS_state"] = NOT_REPORTED
            r["card_note"] = " | ".join([wd["instrument_note"], wd["ext_check_note"]])
        else:
            r["deCODE_arm"], r["ARIC_arm"] = wd.get("deCODE_arm", NOT_REPORTED), wd.get("ARIC_arm", NOT_REPORTED)
            r["EAS_state"] = wd["EAS_state"]
            r["card_note"] = NOT_REPORTED
        # numeric cells
        srcs, feeding = [], {}
        for col in T2_NUMERIC:
            r[col] = NOT_REPORTED
        spec = list(NUMSPEC[pair_key]) + (EXT_SPEC if pair_key in EXT_PAIRS else [])
        cl = []
        if pair_key in EXT_PAIRS:      # the pair's external-check row: it holds the class word and the Wald ratio, SE and p
            cl = [x for x in L.values() if x["claim"].startswith(pair_key + ": Wald ratio in")]
            _req(len(cl) == 1, "EXTCLASS", "%s: %d external-check rows in the ledger" % (pair_key, len(cl)))
            spec += [(col, cl[0]["id"], leaf) for col, leaf in EXT_WALD_SPEC]
        for col, rid, sel in spec:
            _req(rid in L, "TAGS", "%s: ledger row %s does not exist" % (col, rid))
            flat, leaf, suffix = _resolve_leaf(ctx, rid, sel, pair_key)
            r[col] = norm_num(_printed_token(ctx, rid, flat))
            _req(_cell_matches_leaf(ctx, r[col], leaf), "NUMCELL", "%s %s: cell %r does not round from the stored leaf %r" % (pair_key, col, r[col], leaf))
            srcs.append("%s=%s%s" % (col, rid, suffix))
            feeding.setdefault(rid, []).append(col)
        if cl:     # the discovery estimate stored in the external-check row is the card's discovery_wald (same quantity, so ext_wald sits beside it without a composed comparison)
            dleaf = _flat_leaves(ctx, cl[0]["id"])[0][EXT_WALD_DISCOVERY_LEAF]
            _req(_cell_matches_leaf(ctx, r["discovery_wald"], dleaf), "WALD", "%s: ledger %s stores the discovery estimate %r, which does not round to discovery_wald %r" %
                 (pair_key, cl[0]["id"], dleaf, r["discovery_wald"]))
        # word cells that come from ledger leaves
        r["ext_resource"], r["ext_class"], r["ext_bonferroni_6"] = NOT_REPORTED, NOT_REPORTED, NOT_REPORTED
        adds = _strip(ctx, cells[2])
        if pair_key in EXT_PAIRS:
            f = oc[pair_key]
            disp = EXT_DISPLAY[f["resource"]]
            _req(disp in adds, "EXTRES", "%s: Table 2 does not name the external resource %s" % (pair_key, disp))
            r["ext_resource"] = disp
            cls = json.loads(cl[0]["value"])[3]
            _req(cls == f["class"] and cls in adds, "EXTCLASS", "%s: class %r (ledger %s) vs file %r vs Table 2" % (pair_key, cls, cl[0]["id"], f["class"]))
            r["ext_class"] = cls
            srcs.append("ext_class=%s@3" % cl[0]["id"])
            feeding.setdefault(cl[0]["id"], []).append("ext_class")
            bm = BONF_RX.findall(adds)
            _req(len(bm) == 1, "BONF", "%s: Table 2 states the Bonferroni class %d times" % (pair_key, len(bm)))
            flag = "false" if bm[0] else "true"
            pval = _l960_flat(ctx, pair_key)[1]
            _req(flag == f["bonferroni_6"].lower() and (pval < cut) == (flag == "true"), "BONF",
                 "%s: Table 2 says %s, file says %s, p %.3g vs cut %.3g" % (pair_key, flag, f["bonferroni_6"], pval, cut))
            r["ext_bonferroni_6"] = flag
            srcs.append("ext_bonferroni_6=%s" % OC_REL.rsplit("/", 1)[1])
        # discovery basis and East Asian basis: the ledger claim of the row feeding the cell, verbatim
        dwr = [s for s in NUMSPEC[pair_key] if s[0] == DISCOVERY_BASIS_COL][0][1]
        r["discovery_basis"] = L[dwr]["claim"]
        er = sorted({rid for col, rid, sel in NUMSPEC[pair_key] if col.startswith("eas_") or col == "eur_p"}, key=_id_key)
        r["eas_basis"] = " | ".join(L[x]["claim"] for x in er) if er else NOT_REPORTED
        # Open Targets columns
        stage, agents, direction = _ot_columns(ctx, gene, wd["stage_dir"])
        r["development_stage"], r["development_agents"], r["genetic_direction_agreement"] = stage, agents, direction
        r["development_statement"] = wd["stage_dir"]      # the ledger word statement verbatim (e.g. ARG1's pegzilarginase note)
        meta = ctx["snap"]["meta"]["data"]["meta"]["dataVersion"]
        r["ot_data_version"], r["ot_query_date"] = "%s.%s" % (meta["year"], meta["month"]), ctx["snap"]["queried"]
        # literature
        _status_log, date, scope = _literature(ctx, gene, outcome, wd["lit"])      # the log's status is checked against the ledger words; it is not carried (it predates the 1d wording)
        prose = [_strip(ctx, c) for c in cells[1:5]]
        for k, c in zip(T2_PROSE, prose):
            r[k] = c
        fc = ctx["cards"].get(cells[0].strip())
        _req(fc is not None, "CARDS", "%s: no full card in %s for the Table 2 row label %r" % (pair_key, CARDS_REL, cells[0].strip()))
        for k, c in zip(T2_FULL, fc):
            r[k] = _strip(ctx, c)
        rb = REPBEFORE_RX.search(prose[0])
        _req(rb, "LITSTAT", "%s: the established_role_reported_before cell has no 'Reported before' statement" % pair_key)
        lead = LEAD_RX.split(prose[0][rb.start():], 1)[0]       # Table 2's own status sentence (Stage 6.2: after the Role sentence)
        _req(lead, "LITSTAT", "%s: the 'Reported before' statement is empty" % pair_key)
        keys = cited_keys(" ".join(cells[1:5]))
        pm = _pmids_of(ctx, keys)
        r["literature_status"], r["search_date"], r["search_scope"] = lead, date, scope
        r["cited_keys"] = "; ".join(keys) if keys else NOT_REPORTED
        r["literature_pmids"] = "; ".join(pm) if pm else NOT_REPORTED
        # ledger ids and provenance
        tags = {t for c in list(cells) + list(fc) for t in PG.TAG.findall(c)}      # Stage 7 A2.4c: the full card's ledger rows stay carried
        _req(all(t in L for t in tags), "TAGS", "%s: Table 2 tags a ledger id that does not exist" % pair_key)
        ids = sorted(tags | set(feeding) | {wid}, key=_id_key)
        r["numeric_cell_sources"] = "; ".join(srcs)
        r["ledger_row_ids"] = "; ".join(ids)
        r["_unfed_tags"] = "; ".join(sorted(tags - set(feeding), key=_id_key))
        md = [("P31R_Table2.md", ctx["md5"][T2_REL]), ("P31R_SM_cards.md", ctx["md5"][CARDS_REL]), ("literature_search_log.tsv", log_md5), ("devstage_opentargets_20261006.json", ctx["md5"][SNAP_REL])]
        if pair_key in EXT_PAIRS:
            md.append(("S1R1_outcome_checks.tsv", ctx["md5"][OC_REL]))
        r["source_md5"] = "; ".join("%s=%s" % x for x in md)
        # licence class
        provs = _providers(ctx, t1_rows, r)
        # R05 row half: the row's class names the providers of ITS values; development-stage prose adds Open Targets / ChEMBL / DGIdb (rel_licence row provider)
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
        import rel_licence as _rl
        r["licence_class"] = _rl.row_class(set(provs) | _rl.table2_row_providers(r["established_role_reported_before"]))
        rows.append(r)
    return rows


def _providers(ctx, t1_rows, r):
    prov = set()
    for t in t1_rows:
        res = _strip_cites(t["resource"])
        if t["role"] == "Exposure":
            prov.add(t["provider"])
        elif t["role"] == "Outcome" and res.endswith(", " + r["outcome"]):
            prov.add(t["provider"])
        elif t["role"].startswith("External outcome check") and r["ext_resource"] != NOT_REPORTED and r["ext_resource"] in res:
            prov.add(t["provider"])
    for arm, name in (("deCODE_arm", "deCODE"), ("ARIC_arm", "ARIC")):
        if r[arm] != NOT_REPORTED and "not tested" not in r[arm]:
            prov.add(name)
    has_eas = any(r[c] != NOT_REPORTED for c in T2_NUMERIC if c.startswith("eas_") or c == "eur_p") or r["row_group"] == "East Asian, not graded"
    if has_eas:
        text = " ".join([r["eas_basis"], r["discovery_basis"]] + [ctx["L"][x].get("note", "") + " " + ctx["L"][x].get("source", "")
                                                                  for x in re.findall(r"L\d{3}[a-z]?", r["numeric_cell_sources"])
                                                                  if x in ("L085", "L086", "L214", "L215", "L266", "L361", "L362", "L363", "L075")]
                        + [r["adds"]])
        found = {name for rx, name in EAS_VOCAB if re.search(rx, text)} | set(EAS_RESOLVED.get(r["pair_key"], ()))
        _req(found, "PROVIDER", "%s: no East Asian provider is named in its ledger text; declare it in EAS_RESOLVED with the file that names it" % r["pair_key"])
        prov |= found
    return prov


# ----------------------------------------------------------------------------------------------------------------------------------------------
# coverage: the declared `not reported` reasons and the declared not-carried ledger rows (F01)
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _named_rows(ctx, gene, outcome):
    """Ledger rows whose key or claim names the pair: the gene token and the outcome token both occur (outcome case-insensitive, so eGFR also meets EGFR)."""
    gx = re.compile(r"(?<![A-Za-z0-9])%s(?![A-Za-z0-9])" % re.escape(gene))
    ox = re.compile(r"(?<![A-Za-z0-9])%s(?![A-Za-z0-9])" % re.escape(outcome), re.I)
    return sorted((rid for rid, row in ctx["L"].items() if gx.search(row["key"] + " " + row["claim"]) and ox.search(row["key"] + " " + row["claim"])), key=_id_key)


def _declared_not_carried():
    _req(set(NC_CODES) == set(NC_REASON) and set(NOT_CARRIED) <= set(NC_REASON), "COVER", "NOT_CARRIED, NC_CODES and NC_REASON disagree")
    dec = {}
    for code in NC_CODES:
        for i in NOT_CARRIED.get(code, []):
            _req(i not in dec, "COVER", "ledger row %s is declared not-carried twice" % i)
            dec[i] = code
    return dec


def _coverage_cells(ctx, r):
    """(not_reported_reasons, ledger_rows_not_carried, named rows) of one card row, computed from the row's own cells. Raises NRDECL on a numeric column that is neither a ledger
    leaf nor declared (or both), COVER on a ledger row that names the pair, is used by no card column and is not declared."""
    pk = r["pair_key"]
    decl = NR_DECL.get(pk)
    _req(decl is not None and set(decl) <= set(NR_CODES) and set(NR_CODES) == set(NR_REASON), "NRDECL", "%s: no (or a malformed) NR_DECL declaration" % pk)
    declared = [c for code in NR_CODES for c in decl.get(code, [])]
    _req(len(declared) == len(set(declared)), "NRDECL", "%s: a numeric column is declared in two reason groups" % pk)
    blank = {c for c in T2_NUMERIC if r[c] == NOT_REPORTED}
    _req(set(declared) == blank, "NRDECL", "%s: the numeric columns that read 'not reported' differ from the declaration (undeclared %s; declared but filled %s)" %
         (pk, sorted(blank - set(declared)), sorted(set(declared) - blank)))
    nr = "; ".join("%s: %s" % (code, ", ".join(sorted(decl[code], key=T2_NUMERIC.index))) for code in NR_CODES if decl.get(code))
    dec = _declared_not_carried()
    used = {x.strip() for x in r["ledger_row_ids"].split(";")}
    named = _named_rows(ctx, r["gene"], r["outcome"])
    miss = [i for i in named if i not in used and i not in dec]
    _req(not miss, "COVER", "%s: %d ledger row(s) name the pair, are used by no card column and are not declared in NOT_CARRIED: %s" % (pk, len(miss), ", ".join(miss)))
    stale = [i for i in named if i in used and i in dec]
    _req(not stale, "COVER", "%s: ledger row(s) declared not-carried are used by the card: %s" % (pk, ", ".join(stale)))
    nc = "; ".join("%s: %s" % (code, ", ".join(ids)) for code in NC_CODES for ids in [[i for i in named if dec.get(i) == code]] if ids)
    return nr, nc or "none", named


def _declare_coverage(ctx, full):
    """Set the two coverage cells of every full card row; check the declared lists themselves (every listed row exists and names a card pair; the duplicates repeat a carried leaf)."""
    dec = _declared_not_carried()
    named_all = set()
    for r in full:
        r["not_reported_reasons"], r["ledger_rows_not_carried"], named = _coverage_cells(ctx, r)
        named_all |= set(named)
    gone = sorted((i for i in dec if i not in ctx["L"]), key=_id_key)
    _req(not gone, "COVER", "declared not-carried ledger row(s) that do not exist: %s" % ", ".join(gone))
    unused = sorted(set(dec) - named_all, key=_id_key)
    _req(not unused, "COVER", "declared not-carried ledger row(s) that name no card pair: %s" % ", ".join(unused))
    _req(set(DUP_OF) == set(NOT_CARRIED["duplicate"]), "COVER", "DUP_OF and the duplicate list disagree")
    for dup, pk in DUP_OF.items():
        r = [x for x in full if x["pair_key"] == pk][0]
        m = [SRC_RX.match(x) for x in r["numeric_cell_sources"].split("; ") if x.startswith("ext_wald_p=")]
        _req(len(m) == 1 and m[0] is not None and m[0].group("flat") is not None, "COVER", "%s: no ext_wald_p source to compare %s with" % (pk, dup))
        leaf = _flat_leaves(ctx, m[0].group("rid"))[0][int(m[0].group("flat"))]
        v = float(ctx["L"][dup]["value"])
        _req(abs(v - leaf) <= 1e-9 * max(abs(v), abs(leaf)), "COVER", "%s is declared a duplicate of the %s Wald p but stores %r against the leaf %r" % (dup, pk, v, leaf))


def _nopmid(ctx):
    """Cited keys of Table 1 and Table 2 with no PMID in the reference tables (F17: recorded in the column definitions, derived at build)."""
    keys = []
    for tab in ("t1", "t2"):
        for cells in ctx[tab][3]:
            for k in cited_keys(" ".join(cells)):
                if k not in keys:
                    keys.append(k)
    return sorted(k for k in keys if not _pmids_of(ctx, [k]))


# ----------------------------------------------------------------------------------------------------------------------------------------------
# withholding: the shipped version from the full version, by declared rules
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _names_in(names, text):
    low = text.lower()
    return any(re.search(r"(?<![a-z0-9])%s(?![a-z0-9])" % re.escape(n), low) for n in names)


def _ot_withheld(col, r, names):
    """Rule OT1 (S2-D34(4)): is this Open Targets-derived cell a ChEMBL-derived drug or stage field? A row has an agent when its agent list is not 'none'; then its stage is set by a drug
    record. The Open Targets-native fields ship: the tractability class of a row with no clinical candidate, a genetic-direction verdict that names no drug, and a statement that names no
    drug (and belongs to a row with no agent). `names` = every drug and agent name in the snapshot."""
    has_agent = r["development_agents"] != "none"
    if col in ("development_stage", "development_agents"):
        return has_agent
    if col == "development_statement":
        return has_agent or _names_in(names, r[col])
    return _names_in(names, r[col])          # genetic_direction_agreement


def _shipped_rows(full, names=frozenset()):
    """(shipped rows, rules). Rules (data columns only, S2-D30(3)): the ChEMBL-derived drug and stage cells of the Open Targets columns (OT1, S2-D34(4)); with the TPMI switch at
    'withheld', the TPMI-derived cells and the TPMI note of card_note for rows whose external resource is TPMI (TPMI1, TPMI3). The four prose columns are not touched."""
    rules, out = [], []
    for r in full:
        s = dict(r)
        pk = r["pair_key"]

        def cell(rid_, col, lic, why):
            rules.append({"rule_id": rid_, "file": "Table2_evidence_cards.tsv", "pair_key": pk, "column": col, "kind": "cell", "licence_class": lic,
                          "reason": why, "shipped_value": WITHHELD, "withheld_value": r[col]})
            s[col] = WITHHELD

        def frag(rid_, col, m_start, new, lic, why):
            _req(m_start is not None, "WITHHOLD", "%s %s: the withholding pattern %s did not match exactly once" % (pk, col, rid_))
            old = r[col][m_start:]
            rules.append({"rule_id": rid_, "file": "Table2_evidence_cards.tsv", "pair_key": pk, "column": col, "kind": "fragment", "licence_class": lic,
                          "reason": why, "shipped_value": new, "withheld_value": old})
            s[col] = r[col][:m_start] + new

        for col in T2_OT_COLS:
            if _ot_withheld(col, r, names):
                cell("OT1", col, "third_party_derived:OpenTargets+ChEMBL",
                     "ChEMBL-derived drug or development-stage field (S2-D34(4)): withheld; the Open Targets-native fields of the card ship")
        if not TPMI_SHIP and r["ext_resource"] == "TPMI":
            for col in T2_TPMI_COLS:
                cell("TPMI1", col, "third_party_derived:TPMI", "TPMI-derived value (S2-D27(8), S2-D29(e)): withheld until the terms are read")
            parts = r["card_note"].split(" | ")
            _req(len(parts) == 2 and "TPMI" in parts[1], "WITHHOLD", "%s: card_note no longer ends with the TPMI note" % pk)
            frag("TPMI3", "card_note", len(parts[0]), " | " + WITHHELD, "third_party_derived:TPMI", "the ledger's TPMI note states the TPMI-derived class")
        out.append(s)
    return out, rules


def _apply_rules(full, rules):
    """The shipped rows implied by the full rows and the rules file (the verifier's own application of the rules)."""
    by = {r["pair_key"]: dict(r) for r in full}
    for ru in rules:
        s = by[ru["pair_key"]]
        cur, wv = s[ru["column"]], ru["withheld_value"]
        if ru["kind"] == "cell":
            _req(cur == wv, "WITHHOLD", "rule %s %s %s: the withheld value is not the full file's cell" % (ru["rule_id"], ru["pair_key"], ru["column"]))
            s[ru["column"]] = ru["shipped_value"]
        else:
            _req(wv and cur.endswith(wv), "WITHHOLD", "rule %s %s %s: the withheld fragment is not the tail of the full cell" % (ru["rule_id"], ru["pair_key"], ru["column"]))
            s[ru["column"]] = cur[:len(cur) - len(wv)] + ru["shipped_value"]
    return [by[r["pair_key"]] for r in full]


# ----------------------------------------------------------------------------------------------------------------------------------------------
# verification: read the written files back and cross the boundary to the ledger and the markdown
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _md_first_cells(text):
    cells = re.findall(r"^\|\s*(?!:?-{3,})([^|]*?)\s*\|", text, flags=re.M)
    return cells[1:]          # drop the header row


def _holds(text, name):
    return [{"file": name, "key": m.group(1).strip(), "text": text[m.start():text.find("]", m.start()) + 1]} for m in HOLD_RX.finditer(text)]


def _snapshot_names(ctx):
    ns, sn = ctx["ns"], ctx["snap"]
    names = set()
    for g in sn["targets"]:
        names |= set(ns["_single"](g))
        names |= {r["drug"]["name"] for r in ns["_rows_of"](g)}
    names |= set(sn.get("drugs", {}))
    names |= {k for k in sn.get("drug_search", {})}
    return {n.lower() for n in names if len(n) >= 4}


def _tpmi_values(ctx):
    f = [r for r in ctx["oc"] if r["pair"] == "FGFR4-T2D (EAS)"]
    _req(len(f) == 1 and f[0]["resource"] == "TPMI", "TPMI", "the outcome-check file no longer has exactly one TPMI row for FGFR4-T2D (EAS)")
    f = f[0]
    if TPMI_SHIP:      # S2-D34(1): the derived values ship; only the raw per-variant effect, SE and frequency must not appear
        return [float(f[c]) for c in ("b_out", "se_out", "eaf_resource_aligned") if f[c] != ""]
    vals = [float(f[c]) for c in ("b_out", "se_out", "p_out_file", "wald", "wald_se", "wald_se_delta", "p", "eaf_resource_aligned", "n_resource") if f[c] != ""]
    obj = json.loads(ctx["L"]["L223"]["value"])
    vals += [float(x) for row in obj if row[0] == "FGFR4-T2D (EAS)" for x in row[1:] if isinstance(x, (int, float))]
    wr = [x for x in ctx["L"].values() if x["claim"].startswith("FGFR4-T2D (EAS): Wald ratio in")]
    _req(len(wr) == 1, "TPMI", "the ledger no longer has exactly one FGFR4-T2D (EAS) external-check row")
    leaves = _flat_leaves(ctx, wr[0]["id"])[0]
    vals += [float(leaves[i]) for _, i in EXT_WALD_SPEC]      # leaves 0-2 only: leaf 3 onward is the East Asian discovery estimate, not TPMI-derived
    return vals


_SCANTOK = re.compile(r"(?<![A-Za-z0-9_.])([+\-−]?\d+(?:\.\d+)?)(\s?×\s?10[⁻⁰¹²³⁴-⁹]+)?")


def _scan_shipped(ctx, t1_ship, t2_ship, names):
    """TPMI and Open Targets scans over every shipped data cell (and the value scan over the TPMI row). The four Table 2 prose columns are exempt (S2-D30(3): they ship verbatim)."""
    PG = ctx["PG"]
    for fname, rows, key in (("Table1", t1_ship, "resource"), ("Table2", t2_ship, "pair_key")):
        for r in rows:
            for col, val in r.items():
                if fname == "Table2" and col in PROSE_EXEMPT:
                    continue
                low = val.lower()
                if "tpmi" in low and not TPMI_SHIP:
                    rest = val
                    for a in TPMI_ALLOWED.get(col, []):
                        rest = rest.replace(a, "")
                    if col == "licence_class":      # the one grammar third_party_derived:A+B+... names TPMI as a provider token (V13)
                        rest = re.sub(r"(?<=[:+])TPMI(?![A-Za-z0-9_.-])", "", rest)
                    _req("tpmi" not in rest.lower(), "TPMI", "%s %s %s: a shipped cell mentions TPMI outside a label/withheld context: %r" % (fname, r[key][:30], col, val[:120]))
                if col not in NOSCAN_COLS:
                    for nm in names:
                        _req(not re.search(r"(?<![a-z0-9])%s(?![a-z0-9])" % re.escape(nm), low), "OT",
                             "%s %s %s: a shipped cell names the Open Targets-derived agent %r" % (fname, r[key][:30], col, nm))
    tv = _tpmi_values(ctx)
    for r in t2_ship:
        if r["ext_resource"] != "TPMI":
            continue
        for col, val in r.items():
            if col in NOSCAN_COLS or col in PROSE_EXEMPT or col.startswith("_"):
                continue
            for m in _SCANTOK.finditer(val):
                digits = re.sub(r"\D", "", m.group(1)).strip("0")
                if len(digits) < (3 if TPMI_SHIP else 2):      # a raw value is identified by three digits; two would collide with the rounded Wald ratio
                    continue
                for t in PG._d1_tokens(m.group(0).strip()):
                    if t.kind == "num":
                        for v in tv:
                            _req(PG._d1_match1(t, v) is None, "TPMI", "%s %s: token %r rounds from a TPMI-derived value (%.6g) of the outcome-check row" % (r["pair_key"], col, m.group(0), v))


def _verify(ctx, out, wh, t1_rows=None):
    """Read the written files back and assert the release-gate conditions for these two tables. Raises RelTablesError(check, ...)."""
    PG, TG, MB, L = ctx["PG"], ctx["TG"], ctx["MB"], ctx["L"]
    out, wh = pathlib.Path(out), pathlib.Path(wh)
    h1, t1 = _tsv_read(out / OUT_T1)
    h2s, t2s = _tsv_read(out / OUT_T2)
    h2f, t2f = _tsv_read(wh / OUT_T2)
    hr, rules = _tsv_read(wh / OUT_RULES)
    _req(h1 == T1_COLS and h2s == T2_COLS and h2f == T2_COLS, "FORMAT", "TSV headers changed")
    dh, dd = _tsv_read(out / OUT_DOC)
    _req(dh == ["file", "column", "definition"] and [(d["file"], d["column"], d["definition"]) for d in dd] == _column_doc(ctx) and all(d["definition"] for d in dd), "FORMAT",
         "the column dictionary does not list exactly the Table 2 then Table 1 columns, in order, each with its definition")
    # ROWID1 / ROWID2: row count and first-column keys, in order, against a second reading of the markdown
    n1, _, hd1, body1, _ = ctx["t1"]
    first1 = _md_first_cells(ctx["t1_text"])
    _req(len(t1) == len(body1) == len(first1), "ROWID1", "Table 1: %d TSV rows, %d parsed rows, %d first-cell matches" % (len(t1), len(body1), len(first1)))
    _req([(r["role"], r["resource"]) for r in t1] == [(_strip(ctx, c[0]), _strip(ctx, c[1])) for c in body1], "ROWID1", "Table 1: (Role, Resource) keys differ from the markdown")
    _req([r["role"] for r in t1] == [_strip(ctx, c) for c in first1], "ROWID1", "Table 1: first-column keys differ from a second reading of the markdown")
    first2 = _md_first_cells(ctx["t2_text"])
    _, _, _, body2, _ = ctx["t2"]
    for rows in (t2s, t2f):
        _req(len(rows) == len(body2) == len(first2), "ROWID2", "Table 2: %d TSV rows, %d parsed rows, %d first-cell matches" % (len(rows), len(body2), len(first2)))
        _req([r["table_row_label"] for r in rows] == [_strip(ctx, c) for c in first2], "ROWID2", "Table 2: first-column keys differ from the markdown")
    # PROSE (S2-D29(b), S2-D30(3)): the four prose columns equal the builder's stripped Table 2 cells in BOTH files, so the shipped prose is identical to the withheld card's
    for r, rs_, cells in zip(t2f, t2s, body2):
        for k, c in zip(T2_PROSE, cells[1:5]):
            want = MB.TAG.sub("", c).replace("**", "")
            _req(r[k] == want, "PROSE", "%s %s: prose cell differs from the builder's stripped cell" % (r["pair_key"], k))
            _req(rs_[k] == want, "PROSE", "%s %s: the shipped prose cell is not identical to the built Table 2 (the prose columns are not withheld)" % (r["pair_key"], k))
    # WITHHOLD: the shipped file is exactly the full file with the rules applied; every rule is read back from the rules file
    _req(hr == ["rule_id", "file", "pair_key", "column", "kind", "licence_class", "reason", "shipped_value", "withheld_value"], "FORMAT", "rules header changed")
    _req(_apply_rules(t2f, rules) == t2s, "WITHHOLD", "the shipped Table 2 is not the full Table 2 with the declared rules applied")
    for ru in rules:
        _req(ru["column"] in T2_COLS and ru["shipped_value"].find(WITHHELD) >= 0, "WITHHOLD", "rule %s does not write the withheld marker" % ru["rule_id"])
    # NUMCELL: every numeric cell equals its ledger leaf at the printed precision (full file); the shipped cell is the same or withheld by a rule
    ruled = {(ru["pair_key"], ru["column"]) for ru in rules if ru["kind"] == "cell"}
    for rf, rs in zip(t2f, t2s):
        srcs = {}
        for part in rf["numeric_cell_sources"].split("; "):
            m = SRC_RX.match(part)
            if m:
                srcs[m.group("col")] = m
            else:
                _req(re.match(r"^ext_bonferroni_6=\S+\.tsv$", part), "NUMCELL", "%s: unparsable source %r" % (rf["pair_key"], part))
        for col in T2_NUMERIC:
            cell = rf[col]
            if cell == NOT_REPORTED:
                _req(col not in srcs, "NUMCELL", "%s %s: has a ledger source but the cell is 'not reported'" % (rf["pair_key"], col))
                _req(rs[col] == NOT_REPORTED, "NUMCELL", "%s %s: shipped cell differs from a 'not reported' cell" % (rf["pair_key"], col))
                continue
            m = srcs.get(col)
            _req(m, "NUMCELL", "%s %s: numeric cell %r has no ledger source" % (rf["pair_key"], col, cell))
            rid = m.group("rid")
            if m.group("flat") is not None:
                leaf = _flat_leaves(ctx, rid)[0][int(m.group("flat"))]
            elif rid == "L223":
                leaf = _flat_leaves(ctx, rid)[0][_l223_flat(ctx, m.group("pair"), int(m.group("j")))[0]]
            else:
                _req(rid == "L960", "NUMCELL", "unexpected nested source %s" % rid)
                leaf = _l960_flat(ctx, m.group("pair"))[1]
            _req(_cell_matches_leaf(ctx, cell, leaf), "NUMCELL", "%s %s: cell %r does not round from ledger %s leaf %r" % (rf["pair_key"], col, cell, rid, leaf))
            _req(rs[col] == cell or ((rf["pair_key"], col) in ruled and rs[col] == WITHHELD), "NUMCELL", "%s %s: shipped cell %r" % (rf["pair_key"], col, rs[col]))
    # PRINTED: where Table 2 prints a stored leaf next to its tag, the cell holds that same number (crosses TSV -> markdown)
    for rf, cells in zip(t2f, body2):
        raw = " ".join(cells)
        for col in T2_NUMERIC:
            m = None
            for part in rf["numeric_cell_sources"].split("; "):
                mm = SRC_RX.match(part)
                if mm and mm.group("col") == col:
                    m = mm
            if m is None:
                continue
            rid = m.group("rid")
            leaf = (_flat_leaves(ctx, rid)[0][int(m.group("flat"))] if m.group("flat") is not None else
                    _flat_leaves(ctx, rid)[0][_l223_flat(ctx, m.group("pair"), int(m.group("j")))[0]] if rid == "L223" else _l960_flat(ctx, m.group("pair"))[1])
            ct = [t for t in PG._d1_tokens(rf[col]) if t.kind == "num"][0]
            for pm in PRINT_RX.finditer(raw):
                if pm.group(2) != rid:
                    continue
                pt = [t for t in PG._d1_tokens(pm.group(1)) if t.kind == "num"]
                if len(pt) == 1 and PG._d1_match1(pt[0], leaf) is not None:
                    _req(abs(pt[0].tv - ct.tv) <= 1e-9 * max(abs(ct.tv), abs(pt[0].tv), 1e-300), "PRINTED",
                         "%s %s: Table 2 prints %r for ledger %s but the cell is %r" % (rf["pair_key"], col, pm.group(1), rid, rf[col]))
    # COVER / NRDECL: the reason cells are recomputed from the written cells and the ledger, in both files
    for rf, rs in zip(t2f, t2s):
        nr, nc, _ = _coverage_cells(ctx, rf)
        for nm, row in (("full", rf), ("shipped", rs)):
            _req(row["not_reported_reasons"] == nr and row["ledger_rows_not_carried"] == nc, "COVER", "%s (%s file): the reason cells differ from the declarations applied to the written cells" % (rf["pair_key"], nm))
    # WORDS: every element of the pair's L370-L377 word row is held by a card column (the literature element is declared superseded)
    for rf in t2f:
        wid, words = _ledger_words(ctx, rf["gene"], rf["outcome"])
        _req(len(words) == len(ROLES[rf["row_group"]]), "WORDS", "%s: %s changed shape" % (rf["pair_key"], wid))
        for role, val in zip(ROLES[rf["row_group"]], words):
            col = ROLE_COLUMN[role]
            _req(col is None or val in rf[col], "WORDS", "%s: the %s element of %s (%r) is held by no card column (expected in %s)" % (rf["pair_key"], role, wid, val[:50], col))
    # MUST (S2-D29 must-carry list): none of these cells reads `not reported`, in either file; the staging note is the ledger word statement itself
    for item, cells_ in MUST_CARRY:
        for rows, nm in ((t2f, "full"), (t2s, "shipped")):
            for pk, col in cells_:
                hit = [x for x in rows if x["pair_key"] == pk]
                _req(len(hit) == 1 and hit[0][col] != NOT_REPORTED, "MUST", "S2-D29 must-carry item %r: %s %s reads 'not reported' in the %s file" % (item, pk, col, nm))
    for item, pk, col, tok in MUST_CARRY_TEXT:
        hit = [x for x in t2f if x["pair_key"] == pk]
        _req(len(hit) == 1 and tok in hit[0][col].lower(), "MUST", "S2-D29 must-carry item %r: %s %s does not hold %r" % (item, pk, col, tok))
        _req([x[col] for x in t2s if x["pair_key"] == pk] == [WITHHELD], "MUST", "%s %s must read %s in the shipped file (Open Targets-derived)" % (pk, col, WITHHELD))
    # LITSTAT: the status begins 'Reported before' and occurs once in Table 2's established_role_reported_before cell in both files (the card never states more than Table 2)
    for rows, nm in ((t2f, "full"), (t2s, "shipped")):
        for x in rows:
            _req(x["literature_status"].startswith("Reported before") and x["established_role_reported_before"].count(x["literature_status"]) == 1, "LITSTAT", "%s (%s file): literature_status is not a 'Reported before' sentence of established_role_reported_before" % (x["pair_key"], nm))
    # LICENCE: tokens contain no whitespace (F18)
    for nm, rows in (("Table1", t1), ("Table2", t2s), ("Table2 full", t2f)):
        for x in rows:
            for tok in x["licence_class"].split("|"):
                _req(re.fullmatch(r"own_derived|third_party_derived:[A-Za-z0-9_.+\-]+", tok), "LICENCE", "%s: licence token %r is not space-free" % (nm, tok))
    # TPMI rows are declared by predicate twice: the external resource and the pair key
    _req({r["pair_key"] for r in t2f if r["ext_resource"] == "TPMI"} == {r["pair_key"] for r in t2f if r["pair_key"].startswith("FGFR4-T2D (EAS)")}, "TPMI",
         "the rows whose external resource is TPMI are not exactly the FGFR4-T2D (EAS) row")
    # LITCELL: each pinned literature clause occurs once in Table 2 and once in its pair's prior-evidence cell (both files)
    lit = ctx["ns"]["_LITCELL"]
    for pair, (cell, toks) in lit.items():
        _req(ctx["t2_text"].count(cell) == 1, "LITCELL", "%s: the pinned clause occurs %d times in Table 2" % (pair, ctx["t2_text"].count(cell)))
        for rows, nm in ((t2f, "full"), (t2s, "shipped")):
            hit = [r for r in rows if r["pair_key"].split(" ")[0] == pair]
            _req(len(hit) == 1 and hit[0]["established_role_reported_before"].count(cell) == 1, "LITCELL", "%s: the pinned clause is not once in the %s file's established_role_reported_before cell" % (pair, nm))
            _req(hit[0]["literature_status"] in (cell, "Reported before: " + cell), "LITCELL", "%s: literature_status is not the pinned clause (alone or after 'Reported before: ') in the %s file" % (pair, nm))
    # HOLDS: strict refuses an orphan hold; every hold is passed through verbatim
    holds = _holds(ctx["t1_text"], "Table1") + _holds(ctx["t2_text"], "Table2")
    res_keys = {h["key"] for h in _holds(ctx["res_text"], "Results")}
    for h in holds:
        _req(not STRICT or h["key"] in res_keys, "HOLD", "strict: %s carries [HOLD %s: ...] that the Results draft does not hold" % (h["file"], h["key"]))
        txt = _strip(ctx, h["text"])
        for rows in ((t2f, t2s) if h["file"] == "Table2" else (t1,)):
            _req(any(txt in v for r in rows for v in r.values()), "HOLD", "%s hold [HOLD %s] was not passed through verbatim" % (h["file"], h["key"]))
    # scans on the shipped files
    _scan_shipped(ctx, t1, t2s, _snapshot_names(ctx))
    return holds


# ----------------------------------------------------------------------------------------------------------------------------------------------
# build
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _set_tpmi_switch(v3):
    """TIER_RULES.tsv row SW1 (kind switch, target TPMI_DERIVED): handling 'ship' (S2-D34(1)) or 'withheld'. Read, never typed here."""
    global TPMI_SHIP
    with open(pathlib.Path(v3) / TIER_RULES_REL, encoding="utf-8", newline="") as f:
        sw = [r for r in csv.DictReader(f, delimiter="	", quoting=csv.QUOTE_NONE) if r["kind"] == "switch" and r["target"] == "TPMI_DERIVED"]
    _req(len(sw) == 1 and sw[0]["handling"].strip() in ("ship", "withheld"), "SWITCH", "TIER_RULES.tsv has no single TPMI_DERIVED switch with handling ship|withheld")
    TPMI_SHIP = sw[0]["handling"].strip() == "ship"


def _build(v3, out, wh, ledger, **texts):
    global HOLDS_FOUND
    v3, out, wh = pathlib.Path(v3), pathlib.Path(out), pathlib.Path(wh)
    _set_tpmi_switch(v3)
    ctx = _make_context(v3, ledger, **texts)
    t1_rows = _build_table1(ctx)
    full = _build_table2_full(ctx, t1_rows)
    unfed = {r["pair_key"]: r.pop("_unfed_tags") for r in full}
    _declare_coverage(ctx, full)
    ship, rules = _shipped_rows(full, _snapshot_names(ctx))
    for s in ship:
        s.pop("_unfed_tags", None)
    _tsv_write(out / OUT_T1, T1_COLS, t1_rows)
    _tsv_write(out / OUT_T2, T2_COLS, ship)
    _tsv_write(wh / OUT_T2, T2_COLS, full)
    rcols = ["rule_id", "file", "pair_key", "column", "kind", "licence_class", "reason", "shipped_value", "withheld_value"]
    _tsv_write(wh / OUT_RULES, rcols, rules)
    _tsv_write(out / OUT_DOC, ["file", "column", "definition"], [{"file": f, "column": c, "definition": d} for f, c, d in _column_doc(ctx)])
    holds = _verify(ctx, out, wh)
    HOLDS_FOUND = [(h["file"], h["text"]) for h in holds]
    ctx["unfed"] = unfed
    ctx["rules_n"] = len(rules)
    return _infos(ctx, t1_rows, full, ship, rules, holds)


def _infos(ctx, t1_rows, full, ship, rules, holds):
    t2_lic = sorted({c for r in full for c in r["licence_class"].split("|")} | {"third_party_derived:OpenTargets"})
    t1_lic = sorted({r["licence_class"] for r in t1_rows})
    hold_txt = "; ".join("%s: [HOLD %s]" % (f, re.match(r"\[HOLD ([^:\]]+)", t).group(1)) for f, t in HOLDS_FOUND) or "none"
    ncols_rules = len(rules)
    t2_src = [LEDGER_REL, T2_REL, LOG_REL, SNAP_REL, OC_REL, EXT2_REL, RES_REL, OLDREFS_REL, KEYTAB_REL] + REFS_RELS
    return [
        {"release_path": OUT_T1, "ship_status": "ship", "source_paths": [LEDGER_REL, T1_REL, OLDREFS_REL, KEYTAB_REL] + REFS_RELS, "producer": PRODUCER,
         "licence_class": "|".join(t1_lic), "rows": len(t1_rows), "columns": len(T1_COLS),
         "note": "Generated from P31R_Table1.md with the manuscript builder's tag handling (ledger tags stripped, [@key] citations kept inline and listed in cited_keys/cited_pmids). "
                 "Provider and licence_class per row from a declared resource->provider rule set (space-free tokens). A verbatim mirror of the manuscript's Table 1, not masked: the TPMI row's case and control counts "
                 "are printed there, so they ship, labelled third_party_derived:TPMI (S2-D30(3); the rule is in the column definitions). Holds: " + hold_txt},
        {"release_path": OUT_T2, "ship_status": "ship", "source_paths": t2_src,
         "producer": PRODUCER, "licence_class": "|".join(t2_lic), "rows": len(ship), "columns": len(T2_COLS),
         "note": "S2-D29 evidence cards, shipped version (S2-D30(3), S2-D34): %d withholding rules replace the ChEMBL-derived drug and stage cells of the Open Targets columns (rule OT1)%s by '%s'; "
                 "the Open Targets-native fields and the derived TPMI values of the external check ship; the four prose columns ship verbatim (equal to the built Table 2 and to the withheld card); the full version is in withheld/%s. Holds passed through verbatim: %s"
                 % (ncols_rules, "" if TPMI_SHIP else " and the TPMI-derived external-check cells and TPMI note (FGFR4-T2D)", WITHHELD, OUT_T2, hold_txt)},
        {"release_path": OUT_DOC, "ship_status": "ship", "source_paths": [LEDGER_REL], "producer": PRODUCER, "licence_class": "own_derived", "rows": len(_column_doc(ctx)), "columns": 3,
         "note": "Definition of every column of Table1_data_sources.tsv and Table2_evidence_cards.tsv (file, column, definition; written by the module; the verifier requires it to list exactly the file columns, in order). "
                 "The reason-code legends and the Table 1 licence rule are in it."},
        {"release_path": OUT_T2, "ship_status": "withheld_pending_terms", "source_paths": t2_src,
         "producer": PRODUCER, "licence_class": "|".join(t2_lic), "rows": len(full), "columns": len(T2_COLS),
         "note": "Full evidence cards (numeric cells = ledger display tokens asserted against stored leaves; prose = builder-stripped Table 2 cells). Contains the TPMI- and Open Targets-derived data values: "
                 "never shipped until the provider terms are read (S2-D27(8), S2-D29(e))."},
        {"release_path": OUT_RULES, "ship_status": "withheld_pending_terms", "source_paths": [LEDGER_REL, T2_REL, SNAP_REL, OC_REL], "producer": PRODUCER,
         "licence_class": "third_party_derived:OpenTargets" + ("" if TPMI_SHIP else "|third_party_derived:TPMI"), "rows": len(rules), "columns": 9,
         "note": "Every cell/fragment that differs between the full and the shipped Table 2, with the withheld text; the release gate re-applies it to the full file and requires equality with the shipped file."},
    ]


def build(v3, out, withheld, ledger):
    """Spec section 2 interface. Writes tables/Table1_data_sources.tsv and tables/Table2_evidence_cards.tsv under out/, the full Table 2 and the withholding rules under withheld/."""
    return _build(v3, out, withheld, ledger)


# ----------------------------------------------------------------------------------------------------------------------------------------------
# self-test: determinism and mutations; each mutation must trip the check it names
# ----------------------------------------------------------------------------------------------------------------------------------------------
def _load_default_ledger(v3):
    import pandas as pd
    return pd.read_csv(pathlib.Path(v3) / LEDGER_REL, dtype=str)


def selftest(v3):
    global STRICT
    v3 = pathlib.Path(v3)
    _set_tpmi_switch(v3)
    ledger = _load_default_ledger(v3)
    root = pathlib.Path(tempfile.mkdtemp(prefix="rel_tables_selftest_"))
    bad = []

    def rep(name, ok, detail=""):
        print("  %-62s %s %s" % (name, "ok" if ok else "MISS", detail))
        if not ok:
            bad.append(name)

    def run(tag, **texts):
        a, w = root / tag / "out", root / tag / "wh"
        return _build(v3, a, w, ledger, **texts), a, w

    def expect(tag, check, **texts):
        try:
            run(tag, **texts)
        except RelTablesError as e:
            rep("build mutation trips %s" % check, e.check == check, "(got %s: %s)" % (e.check, str(e)[:90]))
            return
        rep("build mutation trips %s" % check, False, "(no error)")

    try:
        base, a1, w1 = run("A")
        _, a2, w2 = run("B")
        same = all((a1 / f).read_bytes() == (a2 / f).read_bytes() for f in (OUT_T1, OUT_T2)) and all((w1 / f).read_bytes() == (w2 / f).read_bytes() for f in (OUT_T2, OUT_RULES))
        rep("two builds are byte-identical", same)
        rep("TSVs are UTF-8 with LF only", all(b"\r" not in (a1 / f).read_bytes() for f in (OUT_T1, OUT_T2)))
        ctx = _make_context(v3, ledger)
        t1 = ctx["t1_text"]; t2 = ctx["t2_text"]
        lit = ctx["ns"]["_LITCELL"]
        clause = list(lit.values())[0][0]
        _req(t2.count(clause) == 1 and " were " in clause, "SELFTEST", "mutation target not unique")
        expect("m_lit", "LITCELL", t2_text=t2.replace(clause, clause.replace(" were ", " was ", 1), 1))
        held = t2.replace("| None as a target;", "| [HOLD ZZ9: planted orphan] None as a target;", 1)
        _req(held != t2 and t2.count("| None as a target;") == 1, "SELFTEST", "hold target not unique")
        STRICT = True
        expect("m_hold", "HOLD", t2_text=held)
        STRICT = False
        try:
            run("m_hold_stage2", t2_text=held)
            rep("stage2 passes an orphan hold through", True)
        except RelTablesError as e:
            rep("stage2 passes an orphan hold through", False, str(e)[:80])
        STRICT = True
        # output mutations on copies of build A
        def mutate(tag, check, edit):
            d = root / tag
            shutil.copytree(root / "A", d)
            edit(d / "out", d / "wh")
            try:
                _verify(ctx, d / "out", d / "wh")
            except RelTablesError as e:
                rep("output mutation trips %s" % check, e.check == check, "(got %s: %s)" % (e.check, str(e)[:80]))
                return
            rep("output mutation trips %s" % check, False, "(no error)")

        def rewrite(path, cols, fn):
            head, rows = _tsv_read(path)
            fn(rows)
            _tsv_write(path, head, rows)

        def numcell(o, w):
            def f(rows):
                r = [x for x in rows if x["pair_key"] == "PCSK9-CAD"][0]
                r["discovery_wald"] = r["discovery_wald"].replace("0.182", "0.183")
                _req(r["discovery_wald"].endswith("0.183"), "SELFTEST", "numcell target")
            rewrite(o / OUT_T2, None, f); rewrite(w / OUT_T2, None, f)
        mutate("m_num", "NUMCELL", numcell)

        def prose(o, w):
            def f(rows):
                r = [x for x in rows if x["pair_key"] == "PAM-T2D"][0]
                _req(r["next_biological_question"].count("amidating") == 1, "SELFTEST", "prose target")
                r["next_biological_question"] = r["next_biological_question"].replace("amidating", "amidating,")
            rewrite(o / OUT_T2, None, f); rewrite(w / OUT_T2, None, f)
        mutate("m_prose", "PROSE", prose)

        def drop(o, w):
            f = lambda rows: rows.pop()
            rewrite(o / OUT_T2, None, f); rewrite(w / OUT_T2, None, f)
        mutate("m_row2", "ROWID2", drop)

        def row1(o, w):
            def f(rows):
                _req(rows[3]["resource"].startswith("FinnGen"), "SELFTEST", "t1 target")
                rows[3]["resource"] = rows[3]["resource"].replace("FinnGen", "FinnGen2", 1)
            rewrite(o / OUT_T1, None, f)
        mutate("m_row1", "ROWID1", row1)

        def tpmi(o, w):
            def f(rows):
                full = _tsv_read(w / OUT_T2)[1]
                s = [x for x in rows if x["pair_key"] == "FGFR4-T2D (EAS)"][0]
                if TPMI_SHIP:       # SW1 ship: a derived TPMI cell masked in the shipped card with no declared rule
                    s["ext_p"] = WITHHELD
                else:               # SW1 withheld: the withheld TPMI note restored in the shipped card
                    s["card_note"] = [x for x in full if x["pair_key"] == "FGFR4-T2D (EAS)"][0]["card_note"]
            rewrite(o / OUT_T2, None, f)
        mutate("m_tpmi", "WITHHOLD", tpmi)

        def ot_stage(o, w):      # a ChEMBL-derived stage cell that rule OT1 withholds is restored in the shipped card
            def f(rows):
                full = _tsv_read(w / OUT_T2)[1]
                s = [x for x in rows if x["pair_key"] == "PCSK9-CAD"][0]
                _req(s["development_stage"] == WITHHELD, "SELFTEST", "ot_stage target is not withheld")
                s["development_stage"] = [x for x in full if x["pair_key"] == "PCSK9-CAD"][0]["development_stage"]
            rewrite(o / OUT_T2, None, f)
        mutate("m_ot_stage", "WITHHOLD", ot_stage)

        def both(fn):      # the same edit to the shipped and the full card
            def edit(o, w):
                rewrite(o / OUT_T2, None, fn); rewrite(w / OUT_T2, None, fn)
            return edit

        def ship_prose(o, w):      # S2-D30(3): the shipped prose must equal Table 2 and the full card
            def f(rows):
                r = [x for x in rows if x["pair_key"] == "FGFR4-T2D (EAS)"][0]
                _req(r["adds"].count("TPMI: ") == 1, "SELFTEST", "ship_prose target")
                r["adds"] = r["adds"].replace("TPMI: ", "TPMI: withheld_pending_terms ", 1)
            rewrite(o / OUT_T2, None, f)
        mutate("m_prose_ship", "PROSE", ship_prose)

        def status(rows):
            r = [x for x in rows if x["pair_key"] == "RAD23B-T2D"][0]
            _req(r["literature_status"].startswith("Not found"), "SELFTEST", "status target")
            r["literature_status"] = "Not reported before"
        mutate("m_status", "LITSTAT", both(status))

        def word(rows):
            r = [x for x in rows if x["pair_key"] == "PAM-T2D"][0]
            _req(r["deCODE_arm"] == "same direction", "SELFTEST", "word target")
            r["deCODE_arm"] = "not reported"
        mutate("m_words", "WORDS", both(word))

        def reasons(rows):
            r = [x for x in rows if x["pair_key"] == "RAD23B-T2D"][0]
            r["not_reported_reasons"] = r["not_reported_reasons"].replace("eas_beta, ", "", 1)
        mutate("m_reasons", "COVER", both(reasons))

        def lic(o, w):
            def f(rows):
                r = [x for x in rows if x["provider"] == "1000Genomes"][0]
                r["licence_class"] = "third_party_derived:1000 Genomes"
            rewrite(o / OUT_T1, None, f)
        mutate("m_licence", "LICENCE", lic)

        def fmt(cells):
            try:
                _tsv_write(root / "fmt.tsv", ["a"], [{"a": cells}])
            except RelTablesError as e:
                return e.check
            return None
        rep("TSV writer refuses a tab in a cell (FORMAT)", fmt("x\ty") == "FORMAT")
        rep("TSV writer refuses a cell that begins with a quote (FORMAT)", fmt('"x') == "FORMAT")
        _tsv_write(root / "fmt.tsv", ["a", "b"], [{"a": 'say "hi"', "b": "x"}])
        rep("TSV writer writes an inner double quote unquoted", (root / "fmt.tsv").read_bytes() == b'a\tb\nsay "hi"\tx\n')

        def plant(pair, text):
            def edit(o, w):
                def f(rows):
                    r = [x for x in rows if x["pair_key"] == pair][0]
                    r["discovery_basis"] = r["discovery_basis"] + text
                rewrite(o / OUT_T2, None, f); rewrite(w / OUT_T2, None, f)
            return edit
        _tv = _tpmi_values(ctx) if TPMI_SHIP else []      # [b_out, se_out, eaf] of the withheld outcome-check row, read at run time (no raw value is typed in this file)
        raw_b = " %.3g" % _tv[0] if TPMI_SHIP else " 0.0012"     # SW1 ship: the raw TPMI effect rounded to three significant digits; the Wald ratio rounds to a different token and is allowed
        mutate("m_tpmi_scan", "TPMI", plant("FGFR4-T2D (EAS)", raw_b))
        mutate("m_ot", "OT", plant("PCSK9-CAD", " alirocumab"))

        raw_s = " %.3g" % _tv[1] if TPMI_SHIP else " TPMI estimate"     # SW1 ship: the raw SE at three significant digits; the word TPMI is then allowed
        mutate("m_tpmi_word", "TPMI", plant("FGFR4-T2D (EAS)", raw_s))

        def drop_doc(o, w):
            head, rows = _tsv_read(o / OUT_DOC)
            _tsv_write(o / OUT_DOC, head, rows[:-1])
        mutate("m_doc", "FORMAT", drop_doc)

        def extra_withheld(o, w):
            def f(rows):
                rows[0]["discovery_p"] = WITHHELD
            rewrite(o / OUT_T2, None, f)
        mutate("m_withhold", "WITHHOLD", extra_withheld)
        # source-level mutations: ledger drift, leaf order, printed precision
        import copy
        def ledger_with(rid, field, fn):
            d = ledger.copy()
            i = d.index[d["id"] == rid][0]
            d.loc[i, field] = fn(d.loc[i, field])
            return d
        def expect_ledger(tag, check, df, **texts):
            try:
                _build(v3, root / tag / "out", root / tag / "wh", df, **texts)
            except RelTablesError as e:
                rep("source mutation trips %s (%s)" % (check, tag), e.check == check, "(got %s: %s)" % (e.check, str(e)[:80]))
                return
            rep("source mutation trips %s (%s)" % (check, tag), False, "(no error)")
        v350 = ledger.loc[ledger["id"] == "L350", "value"].iloc[0]
        _req(v350 == "0.181726", "SELFTEST", "L350 value changed: update the mutation")
        expect_ledger("m_ledger_val", "NUMCELL", ledger_with("L350", "value", lambda v: "0.2"))
        expect_ledger("m_ledger_disp", "NUMCELL", ledger_with("L350", "display", lambda v: "+0.183"))
        l960 = json.loads(ledger.loc[ledger["id"] == "L960", "value"].iloc[0])
        l960[1] = 0.039
        expect_ledger("m_l960", "NUMCELL", ledger_with("L960", "value", lambda v: json.dumps(l960)))
        expect_ledger("m_printed", "PRINTED", ledger, t2_text=t2.replace("PP.H4 0.989 [[L357]]", "PP.H4 0.9892 [[L357]]", 1))
        _req(t2.count("PP.H4 0.989 [[L357]]") == 1, "SELFTEST", "printed target not unique")
        import pandas as pd
        row = ledger[ledger["id"] == "L749"].copy()
        _req(len(row) == 1 and "RAD23B_T2D" in row["key"].iloc[0], "SELFTEST", "cover target")
        expect_ledger("m_cover", "COVER", pd.concat([ledger, row.assign(id="L9999")], ignore_index=True))
        expect_ledger("m_dup", "COVER", ledger_with("L907", "value", lambda v: "0.5"))

        def patched(tag, check, mut):      # a declaration edited in memory; the globals are restored
            saved = copy.deepcopy((NUMSPEC, NR_DECL, NOT_CARRIED))
            try:
                mut()
                try:
                    _build(v3, root / tag / "out", root / tag / "wh", ledger)
                except RelTablesError as e:
                    rep("declaration mutation trips %s (%s)" % (check, tag), e.check == check, "(got %s: %s)" % (e.check, str(e)[:80]))
                    return
                rep("declaration mutation trips %s (%s)" % (check, tag), False, "(no error)")
            finally:
                for obj, old in zip((NUMSPEC, NR_DECL, NOT_CARRIED), saved):
                    obj.clear(); obj.update(old)

        def drop_leaf():       # RAD23B's PP.H3 cell silently defaults to 'not reported' again (the F01 failure)
            n = len(NUMSPEC["RAD23B-T2D"])
            NUMSPEC["RAD23B-T2D"] = [x for x in NUMSPEC["RAD23B-T2D"] if x[0] != "reg_pp_h3"]
            _req(len(NUMSPEC["RAD23B-T2D"]) == n - 1, "SELFTEST", "drop_leaf target")
        patched("m_nrdecl", "NRDECL", drop_leaf)

        def drop_must():       # consistent declarations, but a S2-D29 must-carry cell is given up
            NUMSPEC["PCSK9-CAD"] = [x for x in NUMSPEC["PCSK9-CAD"] if x[0] != "eas_pp_h3"]
            NR_DECL["PCSK9-CAD"]["no_row"].append("eas_pp_h3")
            NOT_CARRIED["protocol_count"].append("L361")      # its source row is then used by no column: declared, so that only MUST can trip
        patched("m_must", "MUST", drop_must)

        def declare_unused():  # a declared not-carried row that names no card pair
            NOT_CARRIED["protocol_count"].append("L001")
        patched("m_unused", "COVER", declare_unused)

        def declare_used():    # a declared not-carried row that the card does use
            NOT_CARRIED["protocol_count"].append("L750")
        patched("m_used", "COVER", declare_used)
        ctx_swap = dict(ctx)
        ctx_swap["oc"] = [ctx["oc"][1], ctx["oc"][0]] + ctx["oc"][2:]
        try:
            _build_table2_full(ctx_swap, _build_table1(ctx_swap))
            rep("source mutation trips NUMCELL (outcome-check row order)", False, "(no error)")
        except RelTablesError as e:
            rep("source mutation trips NUMCELL (outcome-check row order)", e.check == "NUMCELL", "(got %s: %s)" % (e.check, str(e)[:80]))
    finally:
        STRICT = True
        shutil.rmtree(root, ignore_errors=True)
    print("rel_tables selftest: %s (%d miss)" % ("PASS" if not bad else "FAIL", len(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    here = pathlib.Path(__file__).resolve()
    v3_root = here.parents[2]
    if "--selftest" in sys.argv[1:]:
        sys.exit(selftest(v3_root))
    if "--build" in sys.argv[1:]:
        i = sys.argv.index("--build")
        o, w = pathlib.Path(sys.argv[i + 1]), pathlib.Path(sys.argv[i + 2])
        for p in (o, w):
            _req(not str(p.resolve()).startswith(str((v3_root / "R2_rewrite" / "release" / "out").resolve())), "SCRATCH", "refusing to write into release/out")
        infos = build(v3_root, o, w, _load_default_ledger(v3_root))
        for d in infos:
            print(json.dumps({k: (v if k != "note" else v[:100]) for k, v in d.items()}, ensure_ascii=False))
        sys.exit(0)
    print(__doc__)
