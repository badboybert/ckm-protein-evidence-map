# -*- coding: utf-8 -*-
"""Stage 3b track T2 -- ledger rows for Table 2 (evidence cards), L350-L399.

Executed by S1R1_build_ledger.py in its own namespace (add, rel, pd, np, V3, CKM, S1, ROWS, LF available).

Layout (revised after the 3b verifiers, finding T2V-09): every numeral Table 2 prints from this file has its OWN
single-valued row (L350-L365), so the number gate can no longer accept a value moved between pairs or between the default
and conservative prior inside one multi-element row. Every value is read from a named source file and asserted against the
rows the main text already uses. The cards' word statements (tier, exposure-side arms, East Asian states, locus facts,
development stage and direction, literature status) are asserted here as well and recorded per pair in L370-L377, so an
upstream change breaks the build instead of silently staling the table (finding T2V-12). Development stage is asserted
against a saved Open Targets snapshot (finding T2V-02) and literature status against the saved search log (T2V-10); both
sit in R2_rewrite/tables/P31R_Table2_provenance/. This file makes no network calls.
"""
import csv as _csv, io as _io, json as _json, os as _os


def _lv(cid):
    """Value of a ledger row already built (the rows the main text tags)."""
    r = [x for x in ROWS if x["id"] == cid]
    assert len(r) == 1, cid
    return _json.loads(r[0]["value"])


_T2 = "Table 2 (evidence cards)"
_PROV = _os.path.join(V3, "R2_rewrite/tables/P31R_Table2_provenance")
_A11 = _os.path.join(V3, "v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv")
_PE = _os.path.join(V3, "v3.2/_engine_root/results/phase4_annotation/table_s_pair_evidence_tiered.csv")
_a = pd.read_csv(_A11).set_index(["gene", "outcome"]).sort_index()
_p = pd.read_csv(_PE).set_index(["gene", "outcome"]).sort_index()
# the ledger's tiered table (phase4) and the A11 copy with East Asian columns must agree on every shared column
assert _a.index.equals(_p.index) and all((_a[c].astype(str) == _p[c].astype(str)).all() for c in _p.columns), "A11 and phase4 tiered tables disagree"
_oc = pd.read_csv(_os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_outcome_checks.tsv"), sep="\t").set_index("pair")
_vep = pd.read_csv(_os.path.join(V3, "V4_local/source_data/annotation/lead_snp_consequences.csv")).set_index(["gene", "outcome"])
_SRC = rel(_A11) + " (= phase4 copy on shared columns, asserted)"

# ---- card word facts read from the tiered table: tier, deCODE arm, ARIC arm, East Asian outcome-side state -------------
_ARM = {"Concordant (deCODE)": "same direction", "Concordant (ARIC)": "same direction", "Not evaluable (weak instrument)": "not evaluable, weak instrument"}
_EAS = {"NOT_EVALUABLE_VARIANT_ABSENT": "not evaluable, instrument absent", "NOT_EVALUABLE_PALINDROMIC_AMBIGUOUS": "not evaluable, palindromic and strand-ambiguous",
        "EAS_SIGNIFICANT_SAME_DIRECTION": "compatible"}
_WANT = {("PCSK9", "CAD"): ("Convergent", "same direction", "same direction", "not evaluable, instrument absent"),
         ("TGFB1", "CAD"): ("Supported", "not evaluable, weak instrument", "not evaluable, weak instrument", "not evaluable, instrument absent"),
         ("PAM", "T2D"): ("Convergent", "same direction", "same direction", "not evaluable, palindromic and strand-ambiguous"),
         ("RAD23B", "T2D"): ("Convergent", "same direction", "same direction", "not evaluable, palindromic and strand-ambiguous"),
         ("ARG1", "T2D"): ("Supported", "same direction", "not evaluable, weak instrument", "compatible")}
_words = {}
for (g, o), want in _WANT.items():
    r = _a.loc[(g, o)]
    got = (r.evidence_tier, _ARM[r.decode_replication_class], _ARM[r.aric_replication_class], _EAS[r.eas_state_reason])
    assert got == want, ((g, o), got, want)
    _words[(g, o)] = list(got)
# the three graded-elsewhere rows are outside the 81 pairs and carry no tier
for _k in [("FGFR4", "T2D"), ("UMOD", "eGFR"), ("DPEP1", "eGFR")]:
    assert _k not in _a.index, ("now among the 81 graded pairs", _k)
assert len(_a) == _lv("L010"), "tiered table no longer holds the 81 pairs"
_genes81 = set(_a.index.get_level_values(0))
assert len(_genes81) == _lv("L011") and not ({"FGFR4", "UMOD"} & _genes81), "FGFR4 or UMOD is now among the 62 proteins of Table S4"

# ---- external outcome checks: class and Bonferroni flag as the cards print them -----------------------------------------
_OC = {"PCSK9-CAD": True, "TGFB1-CAD": True, "PAM-T2D": True, "RAD23B-T2D": False, "ARG1-T2D": False, "FGFR4-T2D (EAS)": True}
for _k, _bf in _OC.items():
    assert _oc.loc[_k, "class"] == "concordant" and bool(_oc.loc[_k, "bonferroni_6"]) == _bf, (_k, _oc.loc[_k, "class"], _oc.loc[_k, "bonferroni_6"])
assert _oc.loc["FGFR4-T2D (EAS)", "resource"] == "TPMI"
_pres = _io.open(_os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), encoding="utf-8").read()
assert "None with CKB, BBJ or KBA" in _pres, "TPMI overlap statement (card: no participants shared with CKB or Biobank Japan) moved"

# ---- L350-L360: discovery estimates and PP.H4, one number per row -------------------------------------------------------
_ONE = [("L350", ("PCSK9", "CAD"), "mr_beta", "discovery Wald ratio (log-odds of CAD per SD of protein)", "%+.3f"),
        ("L351", ("PCSK9", "CAD"), "mr_pval", "discovery MR p", "%.1e"),
        ("L352", ("PCSK9", "CAD"), "PP.H4_p12_1e6", "PP.H4 under the conservative prior", "%.3f"),
        ("L353", ("TGFB1", "CAD"), "mr_pval", "discovery MR p", "%.1e"),
        ("L354", ("TGFB1", "CAD"), "pp_h4", "PP.H4 under the default prior", "%.3f"),
        ("L355", ("PAM", "T2D"), "mr_beta", "discovery Wald ratio (log-odds of T2D per SD of protein)", "%+.3f"),
        ("L356", ("PAM", "T2D"), "mr_pval", "discovery MR p", "%.1e"),
        ("L357", ("PAM", "T2D"), "pp_h4", "PP.H4 under the default prior", "%.3f"),
        ("L358", ("PAM", "T2D"), "PP.H4_p12_1e6", "PP.H4 under the conservative prior", "%.3f"),
        ("L359", ("RAD23B", "T2D"), "mr_beta", "discovery Wald ratio (log-odds of T2D per SD of protein)", "%+.3f"),
        ("L360", ("RAD23B", "T2D"), "pp_h4", "PP.H4 under the default prior", "%.3f")]
_v = {}
for cid, (g, o), col, what, fmt in _ONE:
    _v[cid] = float(_a.loc[(g, o), col])
# cross-checks against the rows the main text already prints
assert abs(_v["L350"] - float(_oc.loc["PCSK9-CAD", "discovery_wald"])) < 1e-6 and abs(_v["L350"] - _lv("L043")[4]) < 1e-6
assert round(float(_a.loc[("PCSK9", "CAD"), "pp_h4"]), 3) == round(_lv("L202")[0], 3)
assert abs(float(_a.loc[("TGFB1", "CAD"), "mr_beta"]) - _lv("L085")[0]) < 1e-3 and round(float(_a.loc[("TGFB1", "CAD"), "PP.H4_p12_1e6"]), 3) == _lv("L214")[2]
assert abs(_v["L355"] - float(_oc.loc["PAM-T2D", "discovery_wald"])) < 1e-6 and abs(_v["L355"] - _lv("L040")[4]) < 1e-6
assert abs(float(_a.loc[("RAD23B", "T2D"), "mr_pval"]) - _lv("L204")[0]) < 1e-12 and round(float(_a.loc[("RAD23B", "T2D"), "PP.H4_p12_1e6"]), 3) == _lv("L204")[1]
assert abs(_v["L359"] - float(_oc.loc["RAD23B-T2D", "discovery_wald"])) < 1e-6
assert abs(float(_a.loc[("ARG1", "T2D"), "mr_beta"]) - _lv("L050")[0]) < 1e-6 and abs(float(_a.loc[("ARG1", "T2D"), "mr_pval"]) - _lv("L050")[1]) < 1e-15
assert [round(float(_a.loc[("ARG1", "T2D"), c]), 3) for c in ("pp_h4", "PP.H4_p12_1e6")] == [round(x, 3) for x in _lv("L051")]
for cid, (g, o), col, what, fmt in _ONE:
    add(cid, _T2, "%s-%s card: %s" % (g, o, what), _v[cid], fmt % _v[cid], _SRC, col, "[D]",
        note="Single-valued row (one printed number per row). Asserted against L043/L202 (PCSK9), L085/L214 (TGFB1), L040 (PAM), L204 (RAD23B).")

# ---- locus facts the cards state in words --------------------------------------------------------------------------------
_pv = _vep.loc[("PCSK9", "CAD")]
assert str(_pv.missense_canonical).upper() == "TRUE" and str(_pv.amino_acids_canonical) == "R/L46", "PCSK9-CAD lead is no longer the missense Arg46Leu"
assert str(_oc.loc["PAM-T2D", "palindrome"]) == "palindromic, AF concordant"
assert str(_oc.loc["RAD23B-T2D", "palindrome"]).startswith("palindromic, AF-uninformative")

# ---- L361-L363: ancestry-matched East Asian results for card pairs the text counts but does not name ---------------------
_cen = pd.read_csv(_os.path.join(V3, "V4_local/staging/d2-eas-census/results_census.tsv"), sep="\t").set_index(["gene", "outcome"])
_ckf = pd.read_csv(_os.path.join(V3, "V4_local/staging/d2-eas-census/results_ckbfree.tsv"), sep="\t").set_index(["gene", "outcome"])
_cp = pd.read_csv(_os.path.join(V3, "V4_local/staging/d2-eas-census/census_pairs.csv")).set_index(["gene", "outcome"])
_s23 = pd.read_csv(_os.path.join(V3, "V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv")).set_index(["protein", "outcome"])
_pc = _cen.loc[("PCSK9", "CAD")]
assert _pc.call == "NOT_SUPPORTED_H3" and not bool(_pc.ckb_overlap) and "BBJ CAD" in _pc.resource_label
_pa, _pb = _cen.loc[("PAM", "T2D")], _ckf.loc[("PAM", "T2D")]
assert _pa.call == "NOT_SUPPORTED_H1" and bool(_pa.ckb_overlap) and _pb.call == "NOT_SUPPORTED_H1" and not bool(_pb.ckb_overlap)
assert float(_pa["PP.H1"]) == max(float(_pa["PP.H%d" % k]) for k in range(5)) and float(_pb["PP.H1"]) == max(float(_pb["PP.H%d" % k]) for k in range(5))
# membership in the counts the text prints (L211: 32 without an East Asian outcome association / 7 distinct-variant pairs)
assert "PCSK9–CAD" in LF["cen_h3_pairs"] and len(LF["cen_h3_pairs"]) == _lv("L211")[1], "PCSK9-CAD is no longer one of the distinct-variant pairs"
assert LF["cen_h1"] == _lv("L211")[0] and _s23.loc[("PAM", "T2D"), "eas_call"] == "NOT_SUPPORTED_H1" and _s23.loc[("PAM", "T2D"), "set"] == "census_network"
assert _cp.loc[("RAD23B", "T2D"), "verdict"] == "NO" and "no cis peak" in _cp.loc[("RAD23B", "T2D"), "reason"]
assert ("RAD23B", "T2D") not in _cen.index
# FGFR4 (earlier screen), UMOD-eGFR (tested separately, L075) and DPEP1-eGFR (not tested) in East Asians
_fg = _s23.loc[("FGFR4", "T2D")]
assert _fg["set"] == "pilot_offnetwork" and str(_fg.eas_instr_rule).startswith("pilot_old: min p first in file order") and int(_fg.eas_n_tied_top) > 1, \
    "FGFR4 instrument is no longer the first of variants tied at the p-value floor under the earlier screen's rule"
assert "CKB-free" in str(_fg.eas_ckbfree_resource) and round(float(_fg.eas_ckbfree_coloc_H4), 3) == round(_lv("L086")[3], 3)
_um = _s23.loc[("UMOD", "eGFR")]
assert _um["set"] == "pilot_offnetwork" and _um.eas_call == "SUPPORTED" and round(float(_um.eas_coloc_H4), 3) == _lv("L075")[0]
assert ("DPEP1", "eGFR") not in _cen.index and ("DPEP1", "eGFR") not in _s23.index and ("UMOD", "eGFR") not in _cen.index
_EASSRC = "V3/V4_local/staging/d2-eas-census/results_census.tsv + results_ckbfree.tsv"
add("L361", _T2, "PCSK9-CAD ancestry-matched (CKB instrument, Biobank Japan CAD): PP.H3 (distinct variants)", float(_pc["PP.H3"]), "%.2f" % float(_pc["PP.H3"]),
    _EASSRC, "PP.H3, call", "[D]", note="One of the seven distinct-variant pairs (L211, asserted). CKB instrument is rs572512, not the European Arg46Leu lead.")
add("L362", _T2, "PAM-T2D ancestry-matched with AGEN (includes CKB): PP.H1, the most probable hypothesis", float(_pa["PP.H1"]), "%.2f" % float(_pa["PP.H1"]),
    _EASSRC, "PP.H1, call", "[D]", note="One of the 32 pairs without an East Asian outcome association (L211, asserted).")
add("L363", _T2, "PAM-T2D ancestry-matched with Biobank Japan (CKB-free): PP.H1, the most probable hypothesis", float(_pb["PP.H1"]), "%.2f" % float(_pb["PP.H1"]),
    _EASSRC, "PP.H1, call", "[D]")

# ---- L364-L365: the two kidney candidates' discovery estimates on the risk-aligned eGFR* scale ----------------------------
_MRF = _os.path.join(V3, "V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv")
_mr = pd.read_csv(_MRF)
_ca = pd.read_csv(_os.path.join(V3, "V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv")).set_index(["Gene", "Endpoint"])
_k = {}
for g in ("UMOD", "DPEP1"):
    r = _mr[(_mr.gene == g) & (_mr.outcome_name == "eGFR") & (_mr.method == "Wald ratio")]
    assert len(r) == 1, g
    _k[g] = (-float(r.b.iloc[0]), float(r.pval.iloc[0]))   # eGFR* = -ln(eGFR) per SD of protein
    assert _k[g][0] > 0 and _ca.loc[(g, "eGFR"), "MR direction"] == "risk" and _ca.loc[(g, "eGFR"), "Status"] == "Drawn", g
_l71 = _lv("L071")
assert abs(_k["UMOD"][1] - _l71[0]) <= 1e-9 * _l71[0] and abs(_k["DPEP1"][1] - _l71[2]) <= 1e-9 * _l71[2], "kidney p values differ from L071"
for cid, g in (("L364", "UMOD"), ("L365", "DPEP1")):
    add(cid, _T2, "%s-eGFR (multi-signal candidate): discovery Wald ratio on the eGFR* scale (-ln(eGFR) per SD of protein; positive = lower eGFR)" % g,
        _k[g][0], "%+.4f" % _k[g][0], rel(_MRF) + " + SuppTable_candidate_edge_audit.csv", "-b (Wald ratio, eGFR); MR direction == risk", "[D]",
        note="p asserted equal to L071; sign flipped to eGFR* as Methods and L078 require for tables.")

# ---- development stage: Open Targets snapshot (finding T2V-02/T2V-03/T2V-04) ---------------------------------------------
_SNAPF = _os.path.join(_PROV, "devstage_opentargets_20261006.json")
_sn = _json.load(_io.open(_SNAPF, encoding="utf-8"))
_STAGE = {"APPROVAL": 4, "PHASE_4": 4, "PHASE_3": 3, "PHASE_2_3": 3, "PHASE_2": 2, "PHASE_1_2": 2, "PHASE_1": 1, "EARLY_PHASE_1": 1}


def _rows_of(g):
    dc = _sn["targets"][g]["data"]["target"]["drugAndClinicalCandidates"]
    return dc["rows"] if dc else []


def _single(g):
    """Agents with a ChEMBL mechanism on the single protein g (nucleic-acid mechanisms excluded), with their highest stage."""
    out = {}
    for r in _rows_of(g):
        for m in ((r["drug"].get("mechanismsOfAction") or {}).get("rows") or []):
            if [t["approvedSymbol"] for t in m["targets"]] == [g] and "RNAI" not in m["actionType"] and "ANTISENSE" not in m["actionType"]:
                out[r["drug"]["name"]] = max(out.get(r["drug"]["name"], 0), _STAGE.get(r["maxClinicalStage"], 0))
    return out


def _tract(g):
    return any(x["value"] and x["modality"] in ("SM", "AB") for x in _sn["targets"][g]["data"]["target"]["tractability"])


def _stage(g):
    s = _single(g)
    if s and max(s.values()) == 4:
        return "approved"
    if s and max(s.values()) >= 1:
        return "in clinical development"
    return "tractable" if _tract(g) else "no known agent or tractability"


_DS = {"PCSK9": "approved", "TGFB1": "in clinical development", "PAM": "tractable", "RAD23B": "tractable", "ARG1": "tractable",
       "FGFR4": "in clinical development", "UMOD": "tractable", "DPEP1": "approved"}
for g, s in _DS.items():
    assert _stage(g) == s, (g, _stage(g), s)
_PCSK9_AB = {"ALIROCUMAB", "EVOLOCUMAB", "TAFOLECIMAB"}   # approved in the drug record (tafolecimab's PCSK9 target row is at phase 3)
assert _PCSK9_AB <= set(_single("PCSK9"))
for _n in _PCSK9_AB:
    _d = _sn["drugs"][_n.lower()]["data"]["drug"]
    assert _d["maximumClinicalStage"] == "APPROVAL" and [[t["approvedSymbol"] for t in m["targets"]] for m in _d["mechanismsOfAction"]["rows"]] == [["PCSK9"]], _n
_tg = _single("TGFB1")
assert set(_tg) == {"FRESOLIMUMAB", "METELIMUMAB", "LY-2382770"}, _tg
_fres = [r for r in _rows_of("TGFB1") if r["drug"]["name"] == "FRESOLIMUMAB"][0]
assert sorted(tuple(t["approvedSymbol"] for t in m["targets"]) for m in _fres["drug"]["mechanismsOfAction"]["rows"]) == [("TGFB1",), ("TGFB2",), ("TGFB3",)], \
    "fresolimumab is no longer listed against all three TGF-beta isoforms"
for _n in ("METELIMUMAB", "LY-2382770"):
    _r = [r for r in _rows_of("TGFB1") if r["drug"]["name"] == _n][0]
    assert [[t["approvedSymbol"] for t in m["targets"]] for m in _r["drug"]["mechanismsOfAction"]["rows"]] == [["TGFB1"]], _n
assert set(_single("FGFR4")) == {"FGF401"}, _single("FGFR4")   # approved pan-FGFR inhibitors act on the family
assert [h["id"] for h in _sn["drug_search"]["roblitinib"]["data"]["search"]["hits"]] == ["CHEMBL3908979"], "roblitinib no longer resolves to FGF401"
assert {r["disease"]["name"] for r in _sn["drugs"]["FGF401"]["data"]["drug"]["indications"]["rows"]} <= {"cancer", "hepatocellular carcinoma"}, "roblitinib is no longer oncology-only"
_ald = [r for r in _rows_of("FGFR4") if r["drug"]["name"] == "ALDAFERMIN"][0]
assert any(m["actionType"] == "ACTIVATOR" and sorted(t["approvedSymbol"] for t in m["targets"]) == ["FGFR4", "KLB"] for m in _ald["drug"]["mechanismsOfAction"]["rows"])
assert set(_single("DPEP1")) == {"CILASTATIN", "CILASTATIN SODIUM"}
for g in ("PAM", "RAD23B", "ARG1", "UMOD"):
    assert not _rows_of(g), (g, "now has a clinical candidate")
_peg = _sn["drugs"]["pegzilarginase"]["data"]["drug"]
assert _peg["name"] == "PEGZILARGINASE" and _peg["drugType"] == "Enzyme" and _peg["maximumClinicalStage"] == "APPROVAL" and not (_peg["mechanismsOfAction"] or {}).get("rows")
assert any(r["maxClinicalStage"] == "APPROVAL" and r["disease"]["name"] == "arginase deficiency" for r in _peg["indications"]["rows"])
# the same stages and agents as Table S4's audit (track D), so the table cannot drift from the supplement
_au = pd.read_csv(_os.path.join(V3, "R2_rewrite/data/P31R_devstage_audit.tsv"), sep="\t", dtype=str).set_index("gene")
assert set(_au.loc["PCSK9", "agents"].split("; ")) == _PCSK9_AB
assert set(_au.loc["TGFB1", "agents"].split("; ")) == set(_tg)
assert _au.loc["DPEP1", "agents"] == "CILASTATIN"
for g in ("PAM", "RAD23B", "ARG1", "FGFR4", "UMOD"):
    assert g not in _au.index, (g, "now in the Table S4 audit: align Table 2 with it")
# genetic directions behind 'agrees with / opposes'
assert _v["L350"] > 0                                                    # PCSK9-CAD risk: inhibitors agree
assert float(_a.loc[("TGFB1", "CAD"), "mr_beta"]) > 0 and _lv("L085")[1] < 0   # inhibition agrees with the European, not the East Asian, sign
assert float(_a.loc[("ARG1", "T2D"), "mr_beta"]) < 0                     # higher ARG1, lower risk: added enzyme agrees
assert _lv("L086")[0] > 0                                                # higher FGFR4, higher risk (East Asian): inhibitors agree, an activator opposes
assert _k["DPEP1"][0] > 0 and float(_a.loc[("DPEP1", "AF"), "mr_beta"]) < 0   # cilastatin agrees with eGFR, opposes AF
_snap_meta = _sn["meta"]["data"]["meta"]["dataVersion"]

# ---- literature status: saved search log (finding T2V-10 / T2-L5) ---------------------------------------------------------
_LOGF = _os.path.join(_PROV, "literature_search_log.tsv")
_log = list(_csv.DictReader(_io.open(_LOGF, encoding="utf-8"), delimiter="\t"))
_LIT = {"PCSK9-CAD": "Established target", "TGFB1-CAD": "Reported before", "PAM-T2D": "Reported before", "RAD23B-T2D": "Not reported before",
        "ARG1-T2D": "Reported before", "FGFR4-T2D": "Reported before", "UMOD-eGFR": "Reported before", "DPEP1-eGFR": "Reported before"}
for _pair, _st in _LIT.items():
    _rr = [x for x in _log if x["pair"] == _pair]
    assert _rr and any(x["status_in_table"] == _st for x in _rr), (_pair, _st)
_nr = [x for x in _log if x["pair"] == "RAD23B-T2D"]
assert sum(x["database"] == "PubMed" for x in _nr) >= 2 and sum(x["database"] == "Europe PMC full text" for x in _nr) >= 1
assert all(x["hits"] == x["screened"] for x in _nr), "a 'not reported before' search was not fully screened"
assert any("36823471" in x["relevant_records"] for x in _log if x["pair"] == "FGFR4-T2D"), "the FGFR4 prior report (Koprulu 2023) left the log"
_refs = {r["key"]: r["pmid"] for r in _csv.DictReader(_io.open(_os.path.join(V3, "R2_rewrite/refs/new_refs_table2.tsv"), encoding="utf-8"), delimiter="\t")}
assert "23747250" not in _refs.values(), "retracted Ezzat 2013 (retraction notice PMID 33027678) must never be cited"
assert _refs.get("new:Koprulu2023FGFR4") == "36823471" and _refs.get("new:Yuan2023ARG1") == "37652020" and _refs.get("new:LiuX2025FGFR4") == "38976909"
assert any("38976909" in x["relevant_records"] for x in _log if x["pair"] == "FGFR4-T2D")

# ---- the shipped literature clauses, pinned: each must appear once in Table 2, and its data-source words must be in the log row --
_T2MD = _os.environ.get("P31R_TABLE2_MD") or _os.path.join(V3, "R2_rewrite/tables/P31R_Table2.md")   # override for mutation tests only
_t2txt = _io.open(_T2MD, encoding="utf-8").read()
# Stage 6.3b (2026-10-20): TGFB1 and PCSK9 clauses re-pinned to the Stage 6.2 Table 2 wording (the "Reported before:" prefix now carries "before"; "established" lower-case after "Reported before:"); claims and log tokens unchanged.
_LITCELL = {"TGFB1-CAD": ("single-variant CAD estimates were published, with the TGFB1 instrument not stated in the main text; whether our external estimate reproduces them is not established [@new:Wang2025TGFB1]", ["PMID 40649979", "the TGFB1 source not stated", "CARDIoGRAMplusC4D 1.3 (1.18-1.44)"]),
            "ARG1-T2D": ("Reported before [@new:Yuan2023ARG1; @new:FawadAliShah2018]; the earlier MR study's T2D data include those of our DIAMANTE check", ["PMID 37652020", "PMID 30297969", "a superset of", "PMID 30233283"]),
            "PCSK9-CAD": ("established target [@new:Sabatine2017FOURIER]; single-variant CAD estimates were published before [@new:Wang2025TGFB1]", ["PMID 28304224", "1.25 (1.13-1.38), p 7.47e-6"])}
# Stage 5 step 1d (S2-D28): clauses moved with the five-column Table 2 and aligned with the R2 Discussion and novelty ledger N1/N2 (TGFB1: reproduction
# not established; ARG1: no SuSiE/single-variant clause, which the R2 text does not carry). Pre-1d pins: ext_table2_preS5S1d.py.
for _pair, (_cell, _toks) in _LITCELL.items():
    assert _t2txt.count(_cell) == 1, ("Table 2 literature clause changed: re-check it against the paper and the log, then update the pin", _pair)
    _lr = " ".join(x["relevant_records"] for x in _log if x["pair"] == _pair)
    for _tk in _toks:
        assert _tk in _lr, ("literature log no longer supports the Table 2 clause", _pair, _tk)

# ---- L370-L377: the word statements of each card, asserted above, recorded per pair --------------------------------------
_W = [("L370", "PCSK9-CAD", _words[("PCSK9", "CAD")] + ["approved: ALIROCUMAB, EVOLOCUMAB, TAFOLECIMAB inhibit PCSK9; agrees", "Established target (PMID 28304224); same CARDIoGRAMplusC4D estimate in PMID 40649979 Figure 2"]),
      ("L371", "TGFB1-CAD", _words[("TGFB1", "CAD")] + ["in clinical development: METELIMUMAB, LY-2382770 (TGF-beta1 only), FRESOLIMUMAB (TGF-beta1, -beta2, -beta3); inhibition agrees with the European sign only",
                                                        "Reported before (PMID 40649979: same direction; UKB-PPP and deCODE exposures; FinnGen coronary-atherosclerosis outcome, not our CAD definition; replicated (P_replication 2.77e-6); MI 1.31; colocalization did not pass, PPH4 0.450; the CARDIoGRAMplusC4D CAD estimate of our external check already published (Figure 2); old:39 pan-TGF-beta neutralisation in mice)"]),
      ("L372", "PAM-T2D", _words[("PAM", "T2D")] + ["tractable, no clinical candidate", "Reported before (PMID 41063246: deCODE/ARIC SomaScan instruments, colocalization, same direction; PMID 30054598)"]),
      ("L373", "RAD23B-T2D", _words[("RAD23B", "T2D")] + ["tractable, no clinical candidate", "Not reported before (3 PubMed + 2 Europe PMC full-text searches, all hits screened)"]),
      ("L374", "ARG1-T2D", _words[("ARG1", "T2D")] + ["tractable, no ARG1 mechanism; pegzilarginase (enzyme, APPROVAL for arginase deficiency, no curated mechanism) does not set the stage",
                                                      "Reported before (PMID 37652020: proteome-wide MR, UKB-PPP and deCODE instruments, DIAGRAM (= DIAMANTE 2018 with UK Biobank, a superset of our external check's data) and FinnGen R8 outcomes, same direction (Figure 2A/B), SuSiE but not single-variant colocalization on deCODE pQTL only (Figure 2C); PMID 30233283)"]),
      ("L375", "FGFR4-T2D (East Asian)", ["not graded", "exposure-side not tested", "TPMI concordant, Bonferroni; no overlap with CKB or BBJ", "earlier-screen instrument, first of tied variants",
                                          "in clinical development: FGF401 (roblitinib), FGFR4-only inhibitor, oncology; aldafermin FGFR4-KLB activator (complex; does not set the stage) opposes",
                                          "Reported before (PMID 36823471: plasma FGFR4 shares the East Asian-specific T2D signal, PP 97%, protein-raising allele raises risk; PMID 38976909: phenome-wide MR, higher FGFR4 raises East Asian T2D risk; both in full text only)"]),
      ("L376", "UMOD-eGFR", ["multi-signal candidate, not graded", "East Asian: tested separately (pilot), colocalized", "tractable (antibody), no clinical candidate",
                             "Reported before (PMIDs 24185693, 29511113)"]),
      ("L377", "DPEP1-eGFR", ["multi-signal candidate, not graded", "East Asian: not tested", "approved: CILASTATIN inhibits DPEP1; agrees with eGFR, opposes AF",
                              "Reported before (PMID 34426578; old:29)"])]
for cid, pr, facts in _W:
    add(cid, _T2, "%s card: word statements (tier or grade; exposure-side arms; East Asian state; development stage and direction; literature status), asserted from sources" % pr,
        facts, "; ".join(facts), rel(_SNAPF) + " (Open Targets %s.%s, queried %s) + " % (_snap_meta["year"], _snap_meta["month"], _sn["queried"]) + rel(_LOGF) + " + " + rel(_A11),
        "evidence_tier, *_replication_class, eas_state_reason; drugAndClinicalCandidates; search log", "[D]",
        note="Words only (no tagged numerals). Every element is asserted above against its source; the development stage follows SM_devstage's single-protein rule.")

# ---- independent checks of the numeric rows (3b verifier, numbers lens) ------------------------------------------------------
for cid in ["L%d" % i for i in range(350, 366)]:
    EXT_CHECKS[cid] = ("MATCH: re-derived blind in fresh Python and R from the raw sources (tiered pair-evidence table; results_census/ckbfree.tsv; "
                       "mr_all_results.csv) by the Stage 3b T2 numbers verifier, 2026-10-06, before the rows were split one value per row")
