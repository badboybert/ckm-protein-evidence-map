# -*- coding: utf-8 -*-
# Stage 3b, track D (Additional file 1 workbook): ledger rows L450-L499.
# Executed by S1R1_build_ledger.py in ITS namespace (add, rel, pd, np, V3, CKM, S1, EXT_CHECKS and the tables read above).
# Every name defined here starts with _af1_ so nothing in the builder's namespace is overwritten.
import io as _af1_io, os as _af1_os, sys as _af1_sys, csv as _af1_csv

_af1_scripts = _af1_os.path.join(V3, "R2_rewrite", "scripts")
if _af1_scripts not in _af1_sys.path:
    _af1_sys.path.insert(0, _af1_scripts)
import P31R_af1_devstage as _af1_dev  # noqa: E402  (one implementation, shared with P31R_build_af1.py)

# ---- L450: development stage after the audit, and mode of action vs the genetically implied direction (Table S4)
_af1_d, _af1_s = _af1_dev.table()
_af1_v = [_af1_s["approved"], _af1_s["clinical"], _af1_s["tractable"], _af1_s["none"],
          len(_af1_s["approved_opposite"]), len(_af1_s["approved_concordant"]), len(_af1_s["approved_unknown"]),
          _af1_s["approved_opposite"], _af1_s["approved_concordant"], _af1_s["approved_unknown"], _af1_s["restaged"], _af1_s["approved_withdrawn"]]
assert _af1_s["restaged"] == ["COL15A1"], _af1_s["restaged"]
add("L450", "Results: map / Table S4", "Development stage of the 62 proteins (direct single-protein agents only, after the audit): approved / in clinical development / tractable / no known agent or tractability; among proteins with an approved direct agent, agent acts opposite / concordant / unknown relative to the genetically implied direction; approved agent listed as withdrawn",
    _af1_v, "%d approved / %d in clinical development / %d tractable / %d none; approved: opposite %d (%s), concordant %d (%s), unknown %d (%s); agent listed as withdrawn: %s" % (
        _af1_v[0], _af1_v[1], _af1_v[2], _af1_v[3], _af1_v[4], ", ".join(_af1_v[7]), _af1_v[5], ", ".join(_af1_v[8]), _af1_v[6], ", ".join(_af1_v[9]), ", ".join(_af1_v[11])),
    "V3/R2_rewrite/scripts/P31R_af1_devstage.py over v3.2 phase4_annotation_master.csv + tiered table + R2_rewrite/data/P31R_devstage_audit.tsv (asserted against data/P31R_devstage_ot_snapshot.json)",
    "drug_tier_direct, audit (direct_single_protein, effect_on_protein, uniprot_location, withdrawal_listed), mr_direction", "[D]",
    note="Stage 3b audit: COL15A1 re-staged (agents act on a collagen group); every direct approved/clinical gene checked against the Open Targets snapshot of 2026-10-06 under one rule (drug-level maximum stage, any jurisdiction); FN1's ocriplasmin row is at WITHDRAWAL and is marked; membrane receptors read 'unknown' (plasma form measured).")
EXT_CHECKS["L450"] = "MATCH: recomputed with separately written code from other sources (UniProt REST locations fetched 2026-10-06; the pipeline's ChEMBL snapshot mechanism_of_action of 2026-09-17; genetic direction from the raw MR table with ln(eGFR) handled): 6 approved / 6 clinical; opposite DPEP1, FN1, PLAU; concordant PCSK9; unknown FLT1, IL6R (V3/R2_rewrite/scripts/P31R_af1_verify_devstage.py, with data/P31R_uniprot_locations.tsv; rerun 2026-10-06 after the snapshot rule added levilimab and tafolecimab and the FN1 withdrawal mark: MATCH). ChEMBL REST itself returned HTTP 500 on 2026-10-06."

# ---- L451: phenome-wide association flags (Methods: every gene flagged on the any-source consensus)
_af1_p = _af1_dev.phewas()
assert _af1_p["consensus_flagged"] == _af1_p["genes"] == 62
add("L451", "Methods / Table S4", "Phenome-wide association flags: genes flagged on the any-source consensus / genes; Open Targets flagged / queryable",
    [_af1_p["consensus_flagged"], _af1_p["genes"], _af1_p["ot_flagged"], _af1_p["ot_queryable"]],
    "%d of %d (any-source consensus); Open Targets %d of %d" % (_af1_p["consensus_flagged"], _af1_p["genes"], _af1_p["ot_flagged"], _af1_p["ot_queryable"]),
    "v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/phase4_annotation_master.csv", "phewas_consensus, phewas_class_ot", "[D]")
_af1_r411 = [r for r in ROWS if r["id"] == "L411"]
if _af1_r411:
    _af1_g = json.loads(_af1_r411[0]["value"])
    assert [_af1_g["any_source"], _af1_g["genes"], _af1_g["ot_flagged"], _af1_g["ot_queryable"]] == [_af1_p["consensus_flagged"], _af1_p["genes"], _af1_p["ot_flagged"], _af1_p["ot_queryable"]], ("L451 vs L411", _af1_g)
EXT_CHECKS["L451"] = "MATCH: S1R1_07_writing_frame/S1R1_DECISIONS_AND_ALTERNATIVES.md section 3 records 62/62 (consensus) and 27/61 (Open Targets) from the earlier manuscript Methods and prioritization_matrix.R"

# ---- L452: the multi-signal candidate audit (Table S6): SuSiE-testable pool / audited / drawn / held / dropped
_af1_ca = pd.read_csv(_af1_os.path.join(V3, "V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv"))
_af1_st = _af1_ca.Status.value_counts().to_dict()
_af1_audp = set(zip(_af1_ca.Gene, _af1_ca.Endpoint))
assert len(_af1_audp) == len(_af1_ca)
assert _af1_audp <= set(zip(_tst.gene, _tst.outcome_name)), "an audited candidate is not in the SuSiE-testable pool (L263)"
assert (_af1_ca["SuSiE PP.H4 (naive)"] >= 0.8).all(), "an audited candidate is below SuSiE PP.H4 0.80 before purity filtering"
assert (_af1_ca["FDR q<0.05"] == "yes").all()
_af1_l222 = [_af1_st.get("Drawn", 0), _af1_st.get("Held (not drawn)", 0), _af1_st.get("Dropped", 0)]
assert _af1_l222 == [16, 1, 9], _af1_l222
assert sum(_af1_l222) == len(_af1_ca)
_af1_kid = int(_af1_ca[_af1_ca.Status == "Drawn"].Endpoint.isin(["eGFR", "CKD"]).sum())
add("L452", "Table S6", "Multi-signal candidate audit: SuSiE-testable pool (second-order-SE FDR, failed coloc.abf, outside the MHC) / audited (SuSiE PP.H4 >= 0.80 before purity filtering) / drawn / not drawn (fit did not converge) / dropped; kidney pairs among the drawn",
    [len(_tst), len(_af1_ca)] + _af1_l222 + [_af1_kid],
    "%d testable; %d audited; %d drawn, %d not drawn (no converged fit), %d dropped; %d kidney among the drawn" % tuple([len(_tst), len(_af1_ca)] + _af1_l222 + [_af1_kid]),
    "V3/V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv; pool as L263", "Status; SuSiE PP.H4 (naive)", "[D]",
    note="The audited set is asserted to lie inside the L263 testable pool; the supplement track's funnel (L403, from the pipeline run logs) is asserted equal when present.")
# cross-derivation: the supplement track derived the same funnel from the run logs (ext_supp.py, L403); assert agreement when that row exists
_af1_r403 = [r for r in ROWS if r["id"] == "L403"]
if _af1_r403:
    _af1_f = json.loads(_af1_r403[0]["value"])
    assert (_af1_f["first_fit_ge080"], _af1_f["first_fit_ge080"] + _af1_f["not_supportive"], _af1_f["kept"], _af1_f["not_converged"],
            _af1_f["below_080"] + _af1_f["no_pure_cs"]) == (len(_af1_ca), len(_tst), _af1_l222[0], _af1_l222[1], _af1_l222[2]), ("L452 vs L403", _af1_f)
    EXT_CHECKS["L452"] = "MATCH: equal to the funnel the supplement track derived independently from the pipeline run logs 29_cand_*.log and 30_candidate_audit.log (L403: 26 of 160 first-fit passes; 16 kept, 1 not converged, 9 dropped), asserted at build"

# ---- L453: STROBE-MR checklist size (Table S1)
_af1_items = list(_af1_csv.DictReader(_af1_io.open(_af1_os.path.join(V3, "R2_rewrite", "data", "P31R_strobe_mr_items.tsv"), encoding="utf-8"), delimiter="\t"))
_af1_ni = sum(1 for r in _af1_items if not r["sub"]); _af1_ns = sum(1 for r in _af1_items if r["sub"])
assert [_af1_ni, _af1_ns] == [20, 30], [_af1_ni, _af1_ns]
add("L453", "Table S1", "STROBE-MR checklist: main items / subitems", [_af1_ni, _af1_ns], "%d items and %d subitems" % (_af1_ni, _af1_ns),
    "PMID 34698778 (where stated); item structure in V3/R2_rewrite/data/P31R_strobe_mr_items.tsv", "abstract: '20 main items and 30 subitems'", "[P]")
EXT_CHECKS["L453"] = "MATCH: the JAMA abstract (PubMed 34698778) states 20 items and 30 subitems; Table 3 of the BMJ explanation paper (PMID 34702754, PMC8546498, Europe PMC XML) parsed to 20 items and 30 subitems, structure identical to P31R_strobe_mr_items.tsv (2026-10-06)"

# ---- L454: date of the Stage 3b development-stage audit (Open Targets GraphQL queries; R2_rewrite/data/P31R_devstage_audit.tsv)
_af1_au = pd.read_csv(_af1_os.path.join(V3, "R2_rewrite", "data", "P31R_devstage_audit.tsv"), sep="	", dtype=str)
_af1_dates = set()
for _af1_s_ in _af1_au.audit_source:
    _af1_dates |= set(__import__("re").findall(r"(20\d\d)-(\d\d)-(\d\d)", _af1_s_))
assert ("2026", "10", "06") in _af1_dates and all(d <= ("2026", "10", "06") for d in _af1_dates), _af1_dates
add("L454", "Table S4 / SM: development stage", "Date of the Stage 3b development-stage audit in Open Targets (every protein with a direct approved or clinical-stage agent)",
    [6, 10, 2026], "6 October 2026", "V3/R2_rewrite/data/P31R_devstage_audit.tsv (audit_source)", "latest query date", "[D]")
