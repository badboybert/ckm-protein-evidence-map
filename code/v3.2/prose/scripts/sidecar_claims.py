#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Numbers the manuscript prints from the files RECOVERED from E: -- each one BUILT from its file.

ONE implementation, imported by every consumer, so the claim list cannot fork:
  * apply_blocks.py      asserts every main-text claim, Table 4 and Table 1 in the manuscript it builds
  * apply_supplement.py  asserts Supplementary Table S12 in the supplement it builds
  * verify_prose.py P20  asserts all of it again against the SHIPPED documents, plus the scope: every recovered
                         file in the engine-root manifest is the source of a claim, or is declared provenance-only

The files came from <E>/20260331-CKM/Paper 3 in two read-only collections (V3/paper3_deposit_request_manual,
2026-09-14; V3/paper3_deposit_request_v2, 2026-09-15) and are deposited into the engine root by
build_engine_root.R. Each claim is the exact substring the manuscript prints, the file(s) it is built from, the
document it is printed in, and the paragraph anchor it must occur in exactly once.

  python sidecar_claims.py                                   # derive every fact, print every claim
  python sidecar_claims.py --side <dir> --against <dump.txt> # check main-text claims against a text dump (pre-build)
  python sidecar_claims.py --scope                           # recovered files that no claim covers
"""
from __future__ import print_function
import csv
import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
V32 = os.path.dirname(os.path.dirname(HERE))
V3 = os.path.dirname(V32)
# SIDECAR_ENGINE points the registry at a COPY of the engine root -- used only by ladder_sidecars.py, which mutates
# that copy and requires each mutation to trip the check it names. verify_prose P20 reports it whenever it is set.
ENGINE_ROOT = os.environ.get("SIDECAR_ENGINE") or os.path.join(V32, "_engine_root")
ENGINE = os.path.join(ENGINE_ROOT, "results")
SIDE = os.path.join(ENGINE, "repro_phase1_hg38")
MANIFEST = os.path.join(ENGINE_ROOT, "ENGINE_ROOT_MANIFEST.csv")
# The colocalization table of the CORRECTED generation (V4 rebase, 2026-09-25). V3/source_data is superseded.
COLOC_CSV = os.path.join(V3, "V4", "source_data", "tables_from_rds", "phase3_coloc__coloc_results_final.csv")
# Files BUILT HERE from the V4 backbone and deposited by build_engine_root.R under results/v4_local/ -- today the
# non-eGFR decoy null with its CKD stratum rebuilt (the recovered A17/A19 are superseded and not deposited).
LOCAL = os.path.join(ENGINE, "v4_local")
SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
RP = "repro_phase1_hg38/"   # claim sources are paths under the engine root's results/, as the manifest lists them
LP = "v4_local/"            # ... and the files built here

# Recovered files the manuscript prints nothing from, each with the reason it is deposited anyway. P20d asserts
# that every recovered file in the manifest is the source of a claim or is listed here -- never both, never neither.
PROVENANCE_ONLY = {
    RP + "phase5_cross_disease/_gene_coords.csv":
        "the gene coordinates the sharing sweep's locus partition is built from; the manuscript prints the "
        "partition's result (the widest locus, claim F3), which is read from the global FWER file",
    RP + "mr_all_results.rds":
        "asserted cell-for-cell identical to the local MR table the engine root is built from "
        "(build_engine_root.R stage 1f), so every MR number re-derived locally is this file's",
    RP + "_fig_results/phase3_coloc/coloc_results_final.rds":
        "asserted cell-for-cell identical to the local colocalization table (build_engine_root.R stage 1f), "
        "which Supplementary Table S12's PP.H4 column is checked against",
}
# The tables asserted cell for cell, and the recovered files they are printed from.
TABLE_SRC = {
    "Table 4": [RP + "table5_prioritisation.csv"],
    "Supplementary Table S12": [RP + "coloc_1mb_sensitivity.csv"],
    "Table 1": ["e_source/R/00_config.R", "e_source/data/intermediate_gwas/frozen/intermediate_gwasinfo.csv"],
}


class ClaimError(Exception):
    pass


def _rows(side, rel):
    p = os.path.join(side, *rel.split("/"))
    if not os.path.exists(p):
        raise ClaimError("deposited sidecar missing: %s -- rebuild the engine root (build_engine_root.R)" % p)
    with io.open(p, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def _true(x):
    return str(x).strip().upper() in ("TRUE", "T", "1")


def _dp(x, nd):
    """Round in DECIMAL, half-up, to nd places. Never binary float formatting: '%.2f' % 0.975 == '0.97'.
    Pass the file's own text where there is one; a computed float is rounded from its shortest repr."""
    from decimal import Decimal, ROUND_HALF_UP
    return str(Decimal(x if isinstance(x, str) else repr(float(x))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def _sci1(p):
    """7.19e-93 -> '7×10⁻⁹³' (one significant digit, as the renal sentence prints it)."""
    m, e = ("%.0e" % p).split("e")
    return "%s×10%s" % (m, str(int(e)).translate(SUP))


def _sci2(p):
    """5e-05 -> '5.00 × 10⁻⁵'."""
    m, e = ("%.2e" % p).split("e")
    return "%s × 10%s" % (m, str(int(e)).translate(SUP))


def _sci1sp(p):
    """6.65e-135 -> '7 × 10⁻¹³⁵' (one significant digit, spaced, as the mechanism sentence prints it)."""
    m, e = ("%.0e" % p).split("e")
    return "%s × 10%s" % (m, str(int(e)).translate(SUP))


def _thou(n):
    return "{:,}".format(int(n))


def facts(side=SIDE, engine=ENGINE):
    f = {}
    tiers = {}
    with io.open(os.path.join(engine, "phase4_annotation", "table_s_pair_evidence_tiered.csv"), encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            tiers[(r["gene"], r["outcome"])] = r["evidence_tier"]
    n_instr = json.load(io.open(os.path.join(engine, "phase1_instruments", "phase1_summary.json"), encoding="utf-8"))["n_proteins_ok"]

    # --- +/-1 Mb colocalization window ---------------------------------------------------------
    c = _rows(side, "coloc_1mb_sensitivity.csv")
    ok = [r for r in c if r["status"] == "ok"]
    sk = [r for r in c if r["status"] == "skipped"]
    if len(ok) + len(sk) != len(c):
        raise ClaimError("coloc_1mb status vocabulary changed: %s" % sorted({r["status"] for r in c}))
    for r in ok:
        h4, h3 = float(r["PP.H4_1mb"]), float(r["PP.H3_1mb"])
        if _true(r["is_coloc_pass_1mb"]) != (h4 >= 0.80 and h4 / max(h3, 1e-300) >= 5):
            raise ClaimError("coloc_1mb pass flag disagrees with the dual rule for %s-%s" % (r["gene"], r["outcome"]))
    if {(r["gene"], r["outcome"]) for r in c} != set(tiers):
        raise ClaimError("the +/-1 Mb pairs are not the %d tiered pairs" % len(tiers))
    # V4 (2026-09-25): the corrected reader re-evaluates every pair at +/-1 Mb. The manuscript's window paragraph
    # has no clause for a pair the wider window cannot rebuild (the V3 clause, HLA-E-CKD, was a reader artifact), so
    # a skipped pair must refuse the build rather than be silently absent from the prose.
    if sk:
        raise ClaimError("%d pair(s) skipped at +/-1 Mb (%s): the window paragraph states that every pair was re-evaluable"
                         % (len(sk), ", ".join("%s-%s" % (r["gene"], r["outcome"]) for r in sk)))
    npass = sum(_true(r["is_coloc_pass_1mb"]) for r in c)
    conv = sorted(k for k, t in tiers.items() if t == "Convergent")
    val = {(r["gene"], r["outcome"]): (float(r["PP.H4_1mb"]) if r["status"] == "ok" else None) for r in c}
    passed = {(r["gene"], r["outcome"]) for r in c if _true(r["is_coloc_pass_1mb"])}
    f["w"] = dict(n=len(c), ok=len(ok), skip=len(sk), passed=npass, flip=len(ok) - npass, notrobust=len(c) - npass,
                  conv=[(g, o, val[(g, o)]) for g, o in conv], conv_all_pass=all(k in passed for k in conv), values=val)

    # --- intermediate -> disease MR-PRESSO -----------------------------------------------------
    p = _rows(side, "phase7_causal_plausibility/phase7_presso.csv")
    if any(r["status"] != "ok" for r in p):
        raise ClaimError("phase7_presso.csv carries non-ok rows; the manuscript counts analysable pairs only")
    f["presso"] = dict(n=len(p), gsig=sum(float(r["presso_global_p"]) < 0.05 for r in p),
                       out=sum(float(r["n_outliers"] or 0) >= 1 for r in p))

    # --- LD re-clumping arm -----------------------------------------------------------------------
    l = _rows(side, "ld_reclump_sensitivity.csv")
    st = {r["status"] for r in l}
    if not st <= {"ok", "partial_ref_coverage", "unavailable_offline"}:
        raise ClaimError("ld_reclump status vocabulary changed: %s" % sorted(st))
    okl = [r for r in l if r["status"] == "ok"]
    if any(str(r["concordant"]).strip().upper() not in ("TRUE", "FALSE") for r in okl):
        raise ClaimError("ld_reclump: an ok row has no TRUE/FALSE concordance")
    if len(l) != n_instr:
        raise ClaimError("ld_reclump covers %d proteins, the engine root instruments %d" % (len(l), n_instr))
    f["clump"] = dict(total=len(l), ok=len(okl), conc=sum(_true(r["concordant"]) for r in okl),
                      chg=sum(not _true(r["concordant"]) for r in okl),
                      partial=sum(r["status"] == "partial_ref_coverage" for r in l),
                      unavail=sum(r["status"] == "unavailable_offline" for r in l))

    # --- sharing: global family-wise p and the locus partition --------------------------------
    g = _rows(side, "phase5_cross_disease/sensitivity_sharing_global_fwer.csv")
    c7 = [r for r in g if r["estimand"] == "collapsed7"]
    ub = [r for r in g if r["estimand"] == "partition_shape" and r.get("unit") == "unbounded_primary"]
    bd = [r for r in g if r["estimand"] == "partition_shape" and r.get("unit") == "bounded_1mb_companion"]
    if len(c7) != 1 or len(ub) != 1 or len(bd) != 1:
        raise ClaimError("global_fwer: expected one collapsed7 row and two partition_shape rows")
    nperm = int(c7[0]["n_perm"])
    fw = float(c7[0]["fwer_p"])
    if abs(fw - 1.0 / nperm) > 1e-12:
        raise ClaimError("family-wise p %r is not the resolution floor 1/%d; the 'p ≤' wording is false" % (fw, nperm))
    if int(float(bd[0]["max_span_bp"])) > int(float(bd[0]["window_bp"])):
        raise ClaimError("the bounded companion partition has a locus wider than its window")
    f["fwer"] = dict(p=fw, nperm=nperm, family=int(c7[0]["family_size"]), excluded=c7[0]["excluded_degenerate"],
                     widest_mb=float(ub[0]["max_span_bp"]) / 1e6)

    # --- renal-bridge enrichment ----------------------------------------------------------------
    rb = {r["threshold"]: float(r["emp_p_prop_uppertail"]) for r in _rows(side, "phase5_cross_disease/sensitivity_renal_bridge_summary.csv")}
    if set(rb) != {"p<0.05", "p<0.01"}:
        raise ClaimError("renal bridge thresholds changed: %s" % sorted(rb))
    f["renal_bridge"] = rb

    # --- sharing permutation nulls: the ORIGINAL E: outputs (2026-09-15 intake) --------------------
    pdir = "_fig_results/phase5_cross_disease/"
    q5, q1 = _rows(side, pdir + "sensitivity_permutation.csv"), _rows(side, pdir + "sensitivity_permutation_p001.csv")
    if len(q5) != 1 or len(q1) != 1:
        raise ClaimError("a permutation null is not single-row")
    q5, q1 = q5[0], q1[0]
    if abs(float(q1["empirical_p"]) - 1.0 / (int(q1["n_perm"]) + 1)) > 1e-15:
        raise ClaimError("the p < 0.01 null is not at its resolution floor 1/(n_perm+1); 'no permutation reached' is false")
    f["perm"] = dict(p05=q5, p01=q1)

    # --- reverse-direction MR -----------------------------------------------------------------------
    rv = _rows(side, "reverse_mr_results.csv")
    f["rev"] = dict(unavail=[(r["gene"], r["disease"]) for r in rv if r["status"] == "unavailable_offline"],
                    ok=[(r["gene"], r["disease"]) for r in rv if r["status"] == "ok"],
                    null=[(r["gene"], r["disease"]) for r in rv if r["status"] == "ok" and not _true(r["reverse_evidence"])],
                    rev=[(r["gene"], r["disease"], float(r["pval"])) for r in rv if r["status"] == "ok" and _true(r["reverse_evidence"])])
    if len(f["rev"]["ok"]) + len(f["rev"]["unavail"]) != len(rv):
        raise ClaimError("reverse MR status vocabulary changed")

    # --- transcript benchmark -------------------------------------------------------------------------
    eq = {r["set"]: r for r in _rows(side, "eqtl_coloc/eqtl_discrimination.csv")}
    if set(eq) != {"strong", "null"}:
        raise ClaimError("eqtl_discrimination sets changed")
    f["eqtl"] = dict(s_genes=int(eq["strong"]["n_genes"]), s_tested=int(eq["strong"]["n_tested"]), s_pass=int(eq["strong"]["n_pass"]),
                     n_genes=int(eq["null"]["n_genes"]), n_tested=int(eq["null"]["n_tested"]), n_pass=int(eq["null"]["n_pass"]),
                     odds=float(eq["strong"]["fisher_or"]), lo=float(eq["strong"]["fisher_ci_lo"]), hi=float(eq["strong"]["fisher_ci_hi"]),
                     p=float(eq["strong"]["fisher_p"]), overlap=int(eq["strong"]["n_null_cis_window_overlap"]))

    # --- PheWAS consensus ---------------------------------------------------------------------------
    pm = _rows(side, "_fig_results/phase4_annotation/phase4_annotation_master.csv")
    cons = [r["phewas_consensus"] for r in pm]
    f["phewas"] = dict(n=len(pm), flagged=sum(1 for x in cons if x.startswith("flag_")),
                       multi=cons.count("flag_multi_source"), single=cons.count("flag_single_source"),
                       ot_clean=sum(1 for r in pm if str(r.get("phewas_class_ot", "")).lower() == "clean"))

    # --- the 20,000-permutation sharing sweep: the abstract, Figure 4d and the Discussion quote it ----
    sw = {}
    for r in _rows(side, pdir + "sensitivity_sharing_sweep.csv"):
        if r["estimand"] == "collapsed7":
            sw[(r["unit"], float(r["threshold"]))] = r
    if not sw:
        raise ClaimError("the sharing sweep has no collapsed7 rows")
    us = sorted(set(u for u, _t in sw))
    ts = sorted(set(t for _u, t in sw))
    if len(sw) != len(us) * len(ts):
        raise ClaimError("the collapsed7 sweep is not a full %d x %d grid" % (len(us), len(ts)))
    if any(_true(r["null_degenerate"]) for r in sw.values()):
        raise ClaimError("a collapsed7 sweep cell has a degenerate null; 'none excluded' is false")
    f["sweep"] = dict(cell=sw, n_units=len(us), n_thr=len(ts), n_cells=len(sw),
                      nperm=int(sw[("protein", 0.01)]["n_perm"]))

    # --- decoy-loci empirical null: two strata, four files -------------------------------------------
    cn = {}
    for tag, root, o_rel, d_rel in (("egfr", side, "coloc_null/null_observed_egfr.csv", "coloc_null/null_draws_egfr.csv"),
                                    ("other", os.path.join(engine, "v4_local"), "coloc_null/null_observed.csv", "coloc_null/null_draws.csv")):
        ob, dr = _rows(root, o_rel), _rows(root, d_rel)
        per, top = {}, {}
        for r in dr:
            k = (r["gene"], r["outcome"])
            per[k] = per.get(k, 0) + 1
            if _true(r["pass"]):
                top[k] = top.get(k, 0) + 1
        if len(set(per.values())) != 1:
            raise ClaimError("the %s decoy null is unbalanced per observed pair: %s" % (tag, sorted(set(per.values()))))
        if set(per) != set((r["gene"], r["outcome"]) for r in ob):
            raise ClaimError("the %s decoy draws do not cover exactly the observed pairs" % tag)
        dp = sum(_true(r["pass"]) for r in dr)
        op = sum(_true(r["pass"]) for r in ob)
        if not op:
            raise ClaimError("the %s stratum has no observed passes, so its FDR is undefined" % tag)
        obs_pass_keys = set((r["gene"], r["outcome"]) for r in ob if _true(r["pass"]))
        mx = max(top.values()) if top else 0
        cn[tag] = dict(n_obs=len(ob), obs_pass=op, n_draw=len(dr), per_pair=list(per.values())[0],
                       draw_pass=dp, rate=100.0 * dp / len(dr), spec=100.0 - 100.0 * dp / len(dr),
                       max_pair=mx, max_pair_genes=sorted(k[0] for k, v in top.items() if v == mx),
                       max_pair_pass=max([v for k, v in top.items() if k in obs_pass_keys] or [0]),
                       fdr=(1.0 * dp / len(dr)) * len(ob) / op * 100.0,
                       by=dict(((r["gene"], r["outcome"]), float(r["H4"])) for r in ob),
                       pass_genes=set(r["gene"] for r in ob if _true(r["pass"])))
    cn["n_pairs"] = cn["egfr"]["n_obs"] + cn["other"]["n_obs"]
    f["cn"] = cn

    # --- the cystatin-C threshold sweep behind the Figure 4e bars -----------------------------------
    f["ms"] = _rows(side, "pam_conditional/egfr_mirage_threshold_sweep.csv")

    # --- phase-6 mechanism layer ---------------------------------------------------------------------
    mm = _rows(side, "_fig_results/phase6_integration/strong_hits_mechanism_summary.csv")
    mi = _rows(side, "_fig_results/phase6_integration/strong_hits_intermediate.csv")
    sig = {}
    for r in mi:
        if r["intermediate_fdr"] and float(r["intermediate_fdr"]) < 0.05:
            sig.setdefault(r["pair_id"], []).append(r)
    cons = set(k for k, rs in sig.items() if any(_true(r["direction_consistent"]) for r in rs))
    dirs = [r for r in mi if r["intermediate_beta"] != "" and r["disease_beta"] != ""]
    worst = min((r for rs in sig.values() for r in rs if not _true(r["direction_consistent"])),
                key=lambda r: float(r["intermediate_fdr"]))
    f["mech"] = dict(n=len(mm), ann=sum(int(r["n_sig_intermediates"]) > 0 for r in mm),
                     unann=sum(int(r["n_sig_intermediates"]) == 0 for r in mm), fdr_any=len(sig),
                     between=sorted(k.replace("~", u"–") for k in set(sig) - cons),
                     conc=sum(_true(r["direction_consistent"]) for r in dirs), tested=len(dirs),
                     worst=(worst["gene"], worst["intermediate"], float(worst["intermediate_fdr"])))

    # --- CELSR2 x LDL-C -------------------------------------------------------------------------------
    ce = _rows(side, "celsr2_ldl_coloc.csv")
    if len(ce) != 1:
        raise ClaimError("celsr2_ldl_coloc.csv is not single-row")
    f["celsr2"] = ce[0]["PP.H4"]

    # --- phase-7 intermediate -> disease census -------------------------------------------------------
    p7a = _rows(side, "phase7_causal_plausibility/intermediate_disease_mr_all.csv")
    p7d = _rows(side, "_fig_results/phase7_causal_plausibility/intermediate_disease_detail.csv")
    if set((r["intermediate"], r["disease"]) for r in p7a) != set((r["intermediate"], r["disease"]) for r in p7d):
        raise ClaimError("the phase-7 census and detail files cover different cells")
    rest = [r for r in p7d if r["pair_flag"] != "tautological"]
    f["p7"] = dict(n=len(p7d), taut=len(p7d) - len(rest), rest=len(rest),
                   bonf=sum(_true(r["sig_bonf"]) for r in rest),
                   nominal=sum(1 for r in rest if _true(r["sig_nominal"]) and not _true(r["sig_bonf"])))

    # --- the tiered evidence table AS DELIVERED (the local build predates the +/-1 Mb and EAS merge) ---
    tt = _rows(side, "_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv")
    if len(tt) != len(tiers):
        raise ClaimError("the delivered tiered table has %d rows, the engine root's %d" % (len(tt), len(tiers)))
    best = {}
    order = {"Convergent": 0, "Supported": 1, "Flagged": 2}
    for r in tt:
        g = r["gene"]
        if g not in best or order[r["evidence_tier"]] < order[best[g]]:
            best[g] = r["evidence_tier"]
    # deCODE-TESTABLE = the protein is in the 2021 release and the variant matched (76 - 27 - 2 = 47); this
    # includes pairs whose instrument is too weak to evaluate, which decode_evaluable_pair excludes (27)
    dev = [r for r in tt if r["decode_status"] not in ("NOT_IN_RELEASE", "NO_EXACT_RSID_MATCH")]
    f["tt"] = dict(n=len(tt), genes=len(best),
                   tiers=dict((t, sum(1 for r in tt if r["evidence_tier"] == t)) for t in order),
                   gene_tiers=dict((t, sum(1 for v in best.values() if v == t)) for t in order),
                   dtest=len(dev), dhigh=sum(_true(r["decode_high_confidence"]) for r in tt),
                   dconc=sum(_true(r["decode_direction_concordant"]) for r in dev),
                   xp=dict((c, sum(1 for r in tt if r["crossplatform_class"] == c))
                           for c in set(r["crossplatform_class"] for r in tt)),
                   eas=dict((s, sum(1 for r in tt if r["eas_outcome_state"] == s))
                            for s in set(r["eas_outcome_state"] for r in tt)),
                   eas_eval=sum(_true(r["eas_outcome_evaluable"]) for r in tt),
                   eas_strong=sum(_true(r["eas_pqtl_any_strong"]) for r in tt),
                   eas_same=sum(1 for r in tt if _true(r["eas_pqtl_any_strong"])
                                and r["eas_pqtl_direction_class"] == "ALL_SAME"),
                   robust=sum(1 for r in tt if float(r["PP.H4_p12_1e6"]) >= 0.80),
                   atten=sum(1 for r in tt if 0.50 <= float(r["PP.H4_p12_1e6"]) < 0.80),
                   fragile=sum(1 for r in tt if float(r["PP.H4_p12_1e6"]) < 0.50),
                   renal=sum(1 for r in tt if r["outcome"] in ("eGFR", "CKD")),
                   renal_genes=[r["gene"] for r in tt if r["outcome"] in ("eGFR", "CKD")])

    # --- transcript colocalization per gene, and its matched control set -----------------------------
    eqg = _rows(side, "eqtl_coloc/eqtl_pqtl_coloc.csv")
    f["eqg"] = dict(n=len(eqg), conc=sum(1 for r in eqg if r["direction"] == "concordant"),
                    disc=sum(1 for r in eqg if r["direction"] == "discordant_FLAG"),
                    controls=len(_rows(side, "eqtl_coloc/eqtl_null_set.csv")))

    # --- locus-attribution census ---------------------------------------------------------------------
    lc = _rows(side, "locus_competitors.csv")
    clean = [r for r in lc if not _true(r["attribution_untestable"])]
    f["lc"] = dict(n=len(lc), contested=sum(_true(r["contested"]) for r in lc),
                   untest=len(lc) - len(clean), a=sum(1 for r in lc if int(r["a"]) >= 1),
                   b=sum(1 for r in lc if int(r["b"]) >= 1), clean=len(clean),
                   clean_empty=sum(1 for r in clean if int(r["n"]) == 0))

    # --- instrument-selection sensitivity -------------------------------------------------------------
    isn = dict((r["param"], r) for r in _rows(side, "phase1_instrument_sensitivity.csv"))
    miss = [k for k in ("p<1e-06 (primary)", "p<5e-08", "F>=30") if k not in isn]
    if miss:
        raise ClaimError("instrument sensitivity params changed; missing %s" % miss)
    f["is"] = dict(prot=int(isn["p<1e-06 (primary)"]["n_proteins"]),
                   instr=int(isn["p<1e-06 (primary)"]["n_instruments"]),
                   p58=int(isn["p<5e-08"]["n_proteins"]), f30=int(isn["F>=30"]["n_proteins"]),
                   f30_lost=int(isn["F>=30"]["n_proteins_lost_vs_primary"]))

    # --- Table 1 inputs: read as TEXT, never sourced ---------------------------------------------------
    cfg = io.open(os.path.join(engine, "e_source", "R", "00_config.R"), encoding="utf-8", errors="replace").read()

    def _vec(name):
        m = re.search(name + r"\s*<-\s*c\((.*?)\)", cfg, re.S)
        if not m:
            raise ClaimError("00_config.R has no %s vector" % name)
        # strip R's integer suffix ONLY from a number: rstrip("L") once ate the last letter of
        # FinnGen_R10_I9_HEARTFAIL, and P20e caught it
        got = dict((k, re.sub(r"^(\d+)L$", r"\1", v.strip().strip('"')))
                   for k, v in re.findall(r'"?([A-Za-z0-9_.]+)"?\s*=\s*"?([^",)\s]+)"?', m.group(1)))
        if not got:
            raise ClaimError("could not parse %s out of 00_config.R" % name)
        return got

    gw = {}
    with io.open(os.path.join(engine, "e_source", "data", "intermediate_gwas", "frozen",
                              "intermediate_gwasinfo.csv"), encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            gw[r["id"]] = r["sample_size"]
    f["t1"] = dict(ids=_vec("OUTCOMES"), n=_vec("OUTCOME_N_TOTAL"), med=gw)

    # --- Table 4 --------------------------------------------------------------------------------------
    with io.open(os.path.join(side, "table5_prioritisation.csv"), encoding="utf-8-sig", newline="") as fh:
        rd = csv.reader(fh)
        f["table4"] = [row for row in rd]

    # --- renal classification: UMOD / DPEP1 and the 52 concordant signals ---------------------------
    mc = _rows(side, "pam_conditional/egfr_mirage_classification.csv")
    by = {r["gene"]: r for r in mc}
    conc = [r["gene"] for r in mc if r["class"] == "concordant"]
    f["renal"] = dict(concordant=len(conc),
                      umod_crea=float(by["UMOD"]["crea_p"]), umod_cys=float(by["UMOD"]["cys_p"]),
                      dpep1_crea=float(by["DPEP1"]["crea_p"]),
                      # retained = concordant signals that also pass the primary filter in the eGFR decoy-null
                      # observed set. That file excludes MHC genes, so a concordant MHC gene counts as not
                      # retained without being tested; today the three such genes also fail the primary rule.
                      retained=sum(1 for g in conc if g in f["cn"]["egfr"]["pass_genes"]),
                      with_cys=sum(1 for g in f["tt"]["renal_genes"] if g in by))
    return f


def claims(f=None, side=SIDE, engine=ENGINE):
    """[{id, anchor, text, src, where, post}] -- src: the recovered file(s) the text is built from, as manifest paths;
    where: 'main' or 'supp'; post=True means the text exists only after this build's own edits."""
    f = f or facts(side, engine)
    out = []

    def add(cid, anchor, text, src, post=False, where="main"):
        if isinstance(src, str):
            src = [src]
        out.append(dict(id=cid, anchor=anchor, text=text, src=list(src), where=where, post=post))

    C1M = RP + "coloc_1mb_sensitivity.csv"
    w = f["w"]
    A = "As a colocalization-window robustness check, the %d colocalization-supported pairs" % w["n"]
    if w["ok"] != w["n"]:
        raise ClaimError("not every pair was re-evaluable at +/-1 Mb, but the paragraph says all were")
    add("W1", A, "the %d colocalization-supported pairs were re-evaluated with the cis window widened to ±1 Mb. "
                 "All %d pairs were re-evaluable" % (w["n"], w["ok"]), C1M)
    add("W2", A, "re-evaluable from the offline cache; %d remained colocalization-supported at ±1 Mb" % w["passed"], C1M)
    add("W3", A, "whereas %d no longer met the threshold under the wider window" % w["flip"], C1M)
    add("W4", A, "These %d window-sensitive pairs indicate" % w["flip"], C1M)
    add("W5", A, "treat the %d/%d concordant pairs as the window-stable core" % (w["passed"], w["n"]), C1M)
    if not w["conv_all_pass"]:
        raise ClaimError("a Convergent-tier pair fails at +/-1 Mb; the window-stable-core sentence is false")
    add("W6", A, "All %d Convergent-tier pairs sit inside this window-stable core (%s at ±1 Mb)"
        % (len(w["conv"]), ", ".join("%s–%s PP.H4 = %.2f" % (g, o, v) for g, o, v in w["conv"])), C1M)
    # W7-W9 RETIRED at V4 (2026-09-25): they carried the clause for the one pair the V3 +/-1 Mb run could not
    # rebuild (HLA-E-CKD). That was the CKDGen reader truncation; the corrected run rebuilds every pair.

    pr = f["presso"]
    add("P1", "the intermediate→disease arm is itself pervaded by heterogeneity",
        "(MR-PRESSO global p < 0.05 in %d of %d analysable pairs, with ≥1 detected outlier in %d)" % (pr["gsig"], pr["n"], pr["out"]),
        RP + "phase7_causal_plausibility/phase7_presso.csv")

    LDR = RP + "ld_reclump_sensitivity.csv"
    k = f["clump"]
    B = "Beyond these five pre-specified analyses"
    add("L1", B, "applied to all %s proteins, single- and multi-instrument alike" % "{:,}".format(k["total"]), LDR)
    add("L2", B, "permitted a clean re-clump for %d of the %s proteins tested; of those, %d returned the same instrument "
                 "set as the distance-based selection and %d did not. The remaining %s had partial reference coverage "
                 "and %d had no offline-mappable clumpable candidate"
        % (k["ok"], "{:,}".format(k["total"]), k["conc"], k["chg"], "{:,}".format(k["partial"]), k["unavail"]), LDR)
    add("L3", "Finally, to guarantee offline reproducibility",
        "could be re-derived for only %d of the %s proteins" % (k["ok"], "{:,}".format(k["total"])), LDR)

    GF = RP + "phase5_cross_disease/sensitivity_sharing_global_fwer.csv"
    fw = f["fwer"]
    add("F1", "Results. The network comprised %d colocalization-supported edges" % w["n"],
        "max-statistic family-wise p ≤ %s" % _sci2(fw["p"]), GF)
    add("F2", "Cross-endpoint sharing was assessed at two pre-specified screens",
        "(%d cells, none excluded for a degenerate null) the max-statistic family-wise p is at the permutation "
        "resolution floor, p ≤ %s (%s permutations)" % (fw["family"], _sci2(fw["p"]), "{:,}".format(fw["nperm"])), GF)
    if fw["excluded"] != "none":
        raise ClaimError("the grid excludes degenerate cells; 'none excluded' is false")
    add("F3", "Analysis units were: individual proteins", "the widest spanning about %.1f Mb" % fw["widest_mb"], GF)

    rb = f["renal_bridge"]
    add("R1", "This multi-endpoint excess is general rather than renal-specific",
        "(p = %.3f at p < 0.05; p = %.2f at p < 0.01)" % (rb["p<0.05"], rb["p<0.01"]),
        RP + "phase5_cross_disease/sensitivity_renal_bridge_summary.csv")

    # The two 1,000-permutation sharing nulls, from the ORIGINAL E: files. Every digit of both Results
    # sentences is built here: a 2026-09-13 local regeneration had rewritten the p < 0.01 clause to 88.4
    # (SD 7.1) while the file the manuscript was written from says 87.955 (SD 6.910), i.e. the inherited
    # 88.0 / 6.9 / 46% / 5.8. Rounded in decimal, half-up, from the file's own text.
    S = "Cross-endpoint sharing was assessed at two pre-specified screens"
    P05 = RP + "_fig_results/phase5_cross_disease/sensitivity_permutation.csv"
    P01 = RP + "_fig_results/phase5_cross_disease/sensitivity_permutation_p001.csv"
    q5, q1 = f["perm"]["p05"], f["perm"]["p01"]
    o5, m5, sd5 = int(q5["observed"]), float(q5["perm_mean"]), float(q5["perm_sd"])
    n5, pe5 = int(q5["n_perm"]), float(q5["empirical_p"])
    se5 = (pe5 * (1 - pe5) / n5) ** 0.5
    add("U1", S, "%d proteins were associated" % o5, P05)
    add("U2", S, "permutation-null mean of %s (SD %s)—a %s%% excess (%s SD)"
        % (_dp(q5["perm_mean"], 1), _dp(q5["perm_sd"], 1), _dp(100.0 * (o5 - m5) / m5, 1), _dp((o5 - m5) / sd5, 1)), P05)
    add("U3", S, "The empirical tail probability at this screen is %s, but at %s permutations its "
                 "Monte-Carlo standard error (MCSE) is %s (95%% interval %s–%s)"     # PASS A: MCSE defined at its first use
        % (_dp(q5["empirical_p"], 3), "{:,}".format(n5), _dp(se5, 3), _dp(pe5 - 1.96 * se5, 3), _dp(pe5 + 1.96 * se5, 3)), P05)
    o1, m1, sd1 = int(q1["observed"]), float(q1["perm_mean"]), float(q1["perm_sd"])
    n1 = int(q1["n_perm"])
    add("U4", S, "%d proteins against a null mean of %s (SD %s)—a %s%% excess (z = %s), with no permutation "
                 "of %s reaching the observed count (empirical p = 1/%s"
        % (o1, _dp(q1["perm_mean"], 1), _dp(q1["perm_sd"], 1), _dp(100.0 * (o1 - m1) / m1, 0),
           _dp((o1 - m1) / sd1, 1), "{:,}".format(n1), "{:,}".format(n1 + 1)), P01)

    REV = RP + "reverse_mr_results.csv"
    rv = f["rev"]
    if len(rv["unavail"]) != 1 or len(rv["rev"]) != 1 or rv["rev"][0][2] >= 0.05:
        raise ClaimError("reverse MR no longer has one untestable pair and one nominal reverse signal")
    add("V1", B, "it was attempted for the single Steiger-reverse-flagged pair (untestable offline—no outcome-side pQTL "
                 "available) and run for %s forward-confirmed control pairs" % WORDS[len(rv["ok"])], REV)
    add("V2", "A formal reverse-direction MR left the testable control pairs forward-consistent",
        "(%s reverse-null; %s–%s showed a nominal reverse signal"
        % (" and ".join("%s–%s" % x for x in rv["null"]), rv["rev"][0][0], rv["rev"][0][1]), REV)
    add("V3", "A formal reverse-direction MR left the testable control pairs forward-consistent",
        "%s–%s could not be reverse-tested" % rv["unavail"][0], REV)

    EQD = RP + "eqtl_coloc/eqtl_discrimination.csv"
    e = f["eqtl"]
    C = "Each colocalization-supported gene was additionally colocalized against whole-blood cis-eQTLs"
    add("E1", C, "Of the %d colocalization-supported genes, %d are present in eQTLGen and were tested, and %d also "
                 "colocalize at the transcript level; %d are absent from eQTLGen"
        % (e["s_genes"], e["s_tested"], e["s_pass"], e["s_genes"] - e["s_tested"]), EQD)
    # PASS A (2026-09-30): 'CI' was used once in the main text and is spelled out; the claim literal follows the text
    add("E2", C, "(%d/%d versus %d/%d; odds ratio %.2f, 95%% confidence interval %.2f–%.2f, p = %.2f"
        % (e["s_pass"], e["s_tested"], e["n_pass"], e["n_tested"], e["odds"], e["lo"], e["hi"], e["p"]), EQD)
    add("E3", C, "%d of those %d control genes share a cis window" % (e["overlap"], e["n_genes"]),
        [EQD, RP + "eqtl_coloc/eqtl_null_set.csv"])

    ph = f["phewas"]
    if ph["flagged"] != ph["n"]:
        raise ClaimError("not every gene is PheWAS-flagged; 'all %d genes flagged' is false" % ph["n"])
    add("H1", "Target genes were annotated across three modules",
        "all %d genes flagged: %d multi-source, %d single-source; Open Targets alone classified %d genes as clean"
        % (ph["n"], ph["multi"], ph["single"], ph["ot_clean"]),
        RP + "_fig_results/phase4_annotation/phase4_annotation_master.csv")

    MIR = RP + "pam_conditional/egfr_mirage_classification.csv"
    rn = f["renal"]
    add("K1", "Two canonical kidney genes, UMOD and DPEP1",
        "(UMOD: eGFR p = %s, cystatin-C-based eGFR p = %s; DPEP1: eGFR p = %s)"
        % (_sci1(rn["umod_crea"]), _sci1(rn["umod_cys"]), _sci1(rn["dpep1_crea"])), MIR)
    add("K2", "proteins with a creatinine-independent eGFR association",
        "of %d proteins with a creatinine-independent eGFR association" % rn["concordant"], MIR)

    # ---- the 20,000-permutation sharing sweep: abstract, Results, Figure 4d, Discussion -------------
    SW = RP + "_fig_results/phase5_cross_disease/sensitivity_sharing_sweep.csv"
    sw = f["sweep"]
    p01, p05 = sw["cell"][("protein", 0.01)], sw["cell"][("protein", 0.05)]
    # FIX round (2026-09-27): the abstract prints the null mean, z and family-wise p of the SAME 20,000-permutation run.
    # Round 1 dropped the z to hold the 350-word cap; fix round 2 restores it (USER decision B) with the words freed by
    # moving FGFR4-T2D out of the abstract. The Results keep the 1,000-permutation z of their own run (claim U-block).
    add("S1", "Proteome-wide effects and network topology identified a CAD–HF–AF–T2D component",
        "versus %s under permutation (z = %s; max-statistic family-wise p" % (_dp(p01["null_mean"], 1), _dp(p01["z"], 1)), SW)
    add("S2", "reached at least two collapsed CKM endpoint-units compared with",
        "compared with %s under the marginal-preserving null (%s permutations)"
        % (_dp(p01["null_mean"], 1), _thou(sw["nperm"])), SW)
    add("S3", "proteins are associated with at least 2 of 7 endpoint-units against a null of",
        "against a null of %s ± %s over %s permutations (z = %s)"
        % (p01["null_mean"], p01["null_sd"], _thou(sw["nperm"]), p01["z"]), SW)
    add("S4", "Cross-endpoint sharing against a coverage-preserving gene-label permutation null",
        "weak at p < 0.05 (z = %s, below the 1.96 reference line)" % p05["z"], SW)
    add("S5", "we report both resolutions rather than the more favourable one",
        "resolves the same cell at %s permutations to p = %s (MCSE %s)"
        % (_thou(sw["nperm"]), _dp(p05["emp_p"], 3), _dp(p05["emp_p_mcse"], 3)), SW)
    add("S6", "The sweep is not uniform across analysis units: at the p < 0.01 screen the excess",
        "protein level (z = %s)" % _dp(p01["z"], 1), SW)
    add("S7", "the pre-specified primary LD control, one representative per single-linkage locus",
        "chained at 1 Mb (z = %s)" % _dp(sw["cell"][("locus_representative", 0.01)]["z"], 1), SW)
    mp = sw["cell"][("locus_minp", 0.01)]
    add("S8", "but not when each locus is represented instead by its minimum-p protein",
        "(z = %s, p = %s)" % (_dp(mp["z"], 2), _dp(mp["emp_p"], 2)), SW)
    add("S9", "none excluded for a degenerate null) the max-statistic family-wise p is at the",
        "%d unit × %d threshold grid (%d cells" % (sw["n_units"], sw["n_thr"], sw["n_cells"]), SW)

    # ---- decoy-loci empirical null -----------------------------------------------------------------
    CN = [LP + "coloc_null/null_observed.csv", LP + "coloc_null/null_draws.csv",
          RP + "coloc_null/null_observed_egfr.csv", RP + "coloc_null/null_draws_egfr.csv"]
    cn = f["cn"]
    eg, ot = cn["egfr"], cn["other"]

    def _ctl(g):
        hits = [v for (gg, o), v in ot["by"].items() if gg == g and "CHD" in o]
        if len(hits) != 1:
            raise ClaimError("expected exactly one CHD control row for %s, found %d" % (g, len(hits)))
        return hits[0]

    add("N1", "The per-test null pass-rate (the filter's false-positive rate under a true null) was",
        "%s%% in the non-eGFR stratum and %s%% for eGFR, i.e. specificity %s%% and %s%%"
        % (_dp(ot["rate"], 2), _dp(eg["rate"], 2), _dp(ot["spec"], 1), _dp(eg["spec"], 1)), CN)
    add("N2", "applied to the observed screen these imply an empirical false-discovery rate of",
        "false-discovery rate of ≈%s%% and ≈%s%%" % (_dp(ot["fdr"], 0), _dp(eg["fdr"], 0)), CN)
    add("N3", "so that side is reported as a distribution rather than by example:",
        "%d of %s decoy pairings (%d pairs × %d donor loci) met the filter (%s%%), the highest for any single "
        "protein being %d of %d" % (eg["draw_pass"], _thou(eg["n_draw"]), eg["n_obs"], eg["per_pair"],
                                    _dp(eg["rate"], 2), eg["max_pair"], eg["per_pair"]), CN)
    add("N4", "so that side is reported as a distribution rather than by example:",
        "against %d of %d (%s%%) of the true pairings"
        % (eg["obs_pass"], eg["n_obs"], _dp(100.0 * eg["obs_pass"] / eg["n_obs"], 1)), CN)
    add("N5", "providing empirical evidence that the colocalization filter is sensitive to at least",
        "PCSK9–CAD (PP.H4 = %s)" % _dp(_ctl("PCSK9"), 3), CN)
    add("N6", "providing empirical evidence that the colocalization filter is sensitive to at least",
        "IL6R–CAD (PP.H4 = %s)" % _dp(_ctl("IL6R"), 3), CN)
    add("N7", "protein–outcome pairs the null covers had its cis-pQTL colocalized against",
        "Each of the %s protein–outcome pairs the null covers" % _thou(cn["n_pairs"]), CN)
    add("N8", "Violins show PP.H4 under the decoy null (light grey) and across all tested pairs",
        "pass-rate was %s%% in the non-eGFR stratum and %s%% for eGFR, implying empirical false-discovery rates "
        "of about %s%% and %s%%" % (_dp(ot["rate"], 2), _dp(eg["rate"], 2), _dp(ot["fdr"], 0), _dp(eg["fdr"], 0)), CN)
    add("N9", "classified Flagged for prior-fragile colocalization (PP.H4 falls from",
        "PP.H4 falls from %s" % _dp(_ctl("FLT1"), 2), CN)

    rn = f["renal"]
    RET = [MIR, RP + "coloc_null/null_observed_egfr.csv"]
    add("N10", "creatinine-independent renal signals. The two strongest of",       # L2-L4: the abstract's next sentence
        "retained %d of %d" % (rn["retained"], rn["concordant"]), RET)
    add("N11", "Supplementary Figure S20), the primary coloc.abf filter retained only",
        "retained only %d (≈%s%%)" % (rn["retained"], _dp(100.0 * rn["retained"] / rn["concordant"], 0)), RET)

    # ---- the cystatin-C loss bars (Figure 4e) -------------------------------------------------------
    MS = RP + "pam_conditional/egfr_mirage_threshold_sweep.csv"
    for i, r in enumerate(f["ms"]):
        add("B%d" % (i + 1), "The bars beside the scatter give the percentage of eGFRcrea-MR-significant proteins",
            "%s%% of n = %s" % (r["pct_mirage"], r["n_classifiable"]), MS)

    # ---- phase-6 mechanism layer -------------------------------------------------------------------
    MM = RP + "_fig_results/phase6_integration/strong_hits_mechanism_summary.csv"
    MI = RP + "_fig_results/phase6_integration/strong_hits_intermediate.csv"
    me = f["mech"]
    add("G1", "The annotation requires both an FDR-significant protein→intermediate association",
        "%d receive a pathway annotation and %d do not" % (me["ann"], me["unann"]), MM)
    add("G2", "own protein→disease effect, and the second condition is doing real work",
        "%d pairs clear FDR on at least one mediator" % me["fdr_any"], MI)
    add("G3", "own protein→disease effect, and the second condition is doing real work",
        "so %d of the unannotated pairs do carry an FDR-significant intermediate association" % len(me["between"]), MI)
    add("G4", "are unannotated because no mediation reading is available for them, not because no signal",
        "(%s)" % ", ".join(me["between"]), MI)
    add("G5", "whose sign is opposite to what mediation through that intermediate would require",
        "(strongest: %s→%s, q = %s)" % (me["worst"][0], me["worst"][1], _sci1sp(me["worst"][2])), MI)
    add("G6", "direction concordance across the testable protein–intermediate–disease comparisons",
        "is %d/%d (%s%%)" % (me["conc"], me["tested"], _dp(100.0 * me["conc"] / me["tested"], 0)), MI)

    add("C1", "the CELSR2 cis-pQTL colocalizes almost completely with LDL cholesterol",
        "(PP.H4 = %s)" % _dp(f["celsr2"], 4), RP + "celsr2_ldl_coloc.csv")

    # ---- phase-7 intermediate -> disease census ----------------------------------------------------
    P7 = [RP + "phase7_causal_plausibility/intermediate_disease_mr_all.csv",
          RP + "_fig_results/phase7_causal_plausibility/intermediate_disease_detail.csv"]
    q = f["p7"]
    add("Q1", "intermediate→disease combinations were estimable; after excluding the",
        "(%d intermediate→disease combinations were estimable" % q["n"], P7)
    add("Q2", "intermediate→disease combinations were estimable; after excluding the",
        "%d remain, of which %d pass the Bonferroni threshold applied over the %d estimable combinations"
        % (q["rest"], q["bonf"], q["n"]), P7)
    add("Q3", "Solid edges are protein → intermediate and dashed edges intermediate → endpoint",
        "%d are Bonferroni-significant over the %d estimable combinations" % (q["bonf"], q["n"]), P7)
    add("Q4", "Solid edges are protein → intermediate and dashed edges intermediate → endpoint",
        "and %d are nominal only" % q["nominal"], P7)

    # ---- the tiered evidence table, as recovered (the local build predates its +/-1 Mb and EAS merge) -
    TT = RP + "_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv"
    tt = f["tt"]
    add("T1", "Proteome-wide effects and network topology identified a CAD–HF–AF–T2D component",
        "%d colocalization-supported edges linking %d proteins" % (tt["n"], tt["genes"]), TT)
    add("T2", "replicated at high confidence; combined deCODE/ARIC evidence classified",
        "Of %d deCODE-testable edges, %d replicated at high confidence" % (tt["dtest"], tt["dhigh"]), TT)
    add("T3", "replicated at high confidence; combined deCODE/ARIC evidence classified",
        "classified %d as concordant, %d as discordant and %d as conflicting"
        % (tt["xp"].get("concordant", -1), tt["xp"].get("discordant", -1), tt["xp"].get("conflicting", -1)), TT)
    # L2-L4 update (2026-09-27): the abstract's outcome-side sentence is "showed no discordant edge (k of n compatible)";
    # a discordant pair would make that sentence false, so the claim text cannot be matched then
    add("T4", "as conflicting. East Asian outcomes showed no discordant edge",
        ("no discordant edge (%d of %d compatible)" % (tt["eas"].get("EAS_COMPATIBLE", -1), tt["eas_eval"]))
        if tt["eas"].get("EAS_DISCORDANT", 0) == 0 else "<an EAS_DISCORDANT pair exists>", TT)
    add("T5", "The exact-instrument EAS pQTL sidecar provided a complementary view",
        "%d had a strong EAS pQTL (F ≥ 10); %d of those %d were directionally concordant"
        % (tt["eas_strong"], tt["eas_same"], tt["eas_strong"]), TT)
    add("T6", "Edge calibration classified",
        "classified %d pairs as Convergent, %d as Supported and %d as Flagged"
        % (tt["tiers"]["Convergent"], tt["tiers"]["Supported"], tt["tiers"]["Flagged"]), TT)
    add("T7", "Classes are assigned at the most conservative prior and the panel prints all three",
        "%d prior-robust (PP.H4 ≥ 0.80), %d prior-attenuated (0.50 ≤ PP.H4 < 0.80) and %d prior-fragile (< 0.50), "
        "the last two together being the %d" % (tt["robust"], tt["atten"], tt["fragile"], tt["atten"] + tt["fragile"]), TT)
    add("T8", "is a gene-level roll-up taking that gene",
        "%d Convergent, %d Supported and %d Flagged over %d genes"
        % (tt["gene_tiers"]["Convergent"], tt["gene_tiers"]["Supported"], tt["gene_tiers"]["Flagged"], tt["genes"]), TT)
    add("T9", "high-confidence means a deCODE instrument with F ≥ 10, a concordant direction and a delta",
        "%d of the %d are direction-concordant" % (tt["dconc"], tt["dtest"]), TT)
    # per PROTEIN (V4): a protein colocalizing with both eGFR and CKD has one instrument and plots as one point
    if rn["with_cys"] != tt["renal"]:
        raise ClaimError("a renal pair lacks a cystatin-C estimate: the Figure 4e sentence says every protein has one")
    add("T10", "The scatter plots the MR effect on cystatin-C-based against creatinine-based eGFR for",
        "for the %d proteins behind the %d colocalization-supported renal pairs" % (len(set(tt["renal_genes"])), tt["renal"]),
        [TT, MIR])

    add("X1", "rather than as negative. Where both effects are estimable,",
        "%d genes agree in direction and %d do not" % (f["eqg"]["conc"], f["eqg"]["disc"]),
        RP + "eqtl_coloc/eqtl_pqtl_coloc.csv")

    # ---- locus-attribution census and instrument-selection sensitivity ------------------------------
    LC = RP + "locus_competitors.csv"
    lc = f["lc"]
    add("L4", "which also pass the colocalization rule for the same endpoint.",
        "%d of the %d pairs are contested in that sense" % (lc["contested"], lc["n"]), LC)
    D = "so for those the neighbours were not excluded by evidence; the two sets overlap"
    add("L5", D, "%d pairs carry at least one untestable rival in the window" % lc["untest"], LC)
    add("L6", D, "no cis instrument (%d pairs)" % lc["a"], LC)
    add("L7", D, "never entered the MR screen (%d pairs)" % lc["b"], LC)
    add("L8", "a rival that failed the MR promotion threshold on a weak instrument is not strong",
        "The remaining %d pairs carry no untestable rival: for %d of them" % (lc["clean"], lc["clean_empty"]), LC)

    IS = RP + "phase1_instrument_sensitivity.csv"
    it = f["is"]
    add("I1", "Multi-SNP outlier diagnostics, however, have no subject among the reported pairs",
        "%s of the %s proteins retain an instrument at the stricter p < 5×10⁻⁸ threshold, and %s at F ≥ 30 (%d lost)"
        % (_thou(it["p58"]), _thou(it["prot"]), _thou(it["f30"]), it["f30_lost"]), IS)
    add("I2", "Circulating protein levels were measured on the Olink Explore 3072 platform",
        "giving %s cis-pQTL instruments for %s plasma proteins" % (_thou(it["instr"]), _thou(it["prot"])), IS)

    ids = [c["id"] for c in out]
    if len(ids) != len(set(ids)):
        raise ClaimError("duplicate claim id(s): %s" % sorted({i for i in ids if ids.count(i) > 1}))
    return out


def s12_values(f=None, side=SIDE, engine=ENGINE):
    """(gene, outcome) -> PP.H4 at +/-1 Mb, or None for the pair the wider window could not evaluate."""
    return (f or facts(side, engine))["w"]["values"]


def table4(f=None, side=SIDE, engine=ENGINE):
    """[header, row, ...] exactly as the deposited table5_prioritisation.csv holds them."""
    return (f or facts(side, engine))["table4"]


def check_text(text_of_anchor, only_pre=False, f=None, side=SIDE, engine=ENGINE, where="main"):
    """text_of_anchor(anchor) -> the paragraph text, in the document `where` names. Returns a list of problems
    (empty = every claim printed in that document holds)."""
    # PASS A stage 2: the claims are written in the legacy 'Supplementary Table/Figure SN' form; the document (from the
    # apply_blocks AF step on, so at the build's own SC assertion and at gate time alike) prints the journal's
    # 'Additional file K: ...' form. journal_form is the one declared transformation and a no-op elsewhere.
    sys.path.insert(0, HERE)
    from addfiles import journal_form
    bad = []
    for c in claims(f, side, engine):
        if c["where"] != where or (only_pre and c["post"]):
            continue
        try:
            t = text_of_anchor(journal_form(c["anchor"]))
        except Exception as ex:
            bad.append("%s: anchor %r: %s" % (c["id"], c["anchor"][:50], ex))
            continue
        n = t.count(journal_form(c["text"]))
        if n != 1:
            bad.append("%s: expected once, found %d: %r" % (c["id"], n, c["text"]))
    return bad


def _cellnorm(s):
    return re.sub(r"\s+", " ", s).strip()


def check_table4(got_rows, f=None, side=SIDE, engine=ENGINE):
    """got_rows: Table 4 as printed, as lists of cell text, header first. [] = cell-for-cell the recovered file."""
    want = [[_cellnorm(c) for c in row] for row in table4(f, side, engine)]
    got = [[_cellnorm(c) for c in row] for row in got_rows]
    if got == want:
        return []
    diffs = [(i, j, x, y) for i, (a, b) in enumerate(zip(got, want)) for j, (x, y) in enumerate(zip(a, b)) if x != y]
    return ["Table 4 differs from table5_prioritisation.csv (%d x %d printed vs %d x %d in the file, %d cells differ%s)"
            % (len(got), len(got[0]) if got else 0, len(want), len(want[0]) if want else 0, len(diffs),
               "; first: row %d col %d prints %r, file %r" % diffs[0] if diffs else "")]


def check_table1(rows, f=None, side=SIDE, engine=ENGINE):
    """rows: Table 1 as printed, as lists of cell text, header first. [] = every dataset identifier and N cell
    equals 00_config.R (the eight outcomes) or the frozen gwasinfo (the six intermediates). Not a claim, because
    a short table row cannot be anchored without quoting the value it asserts."""
    t1 = (f or facts(side, engine))["t1"]
    want = {}
    for trait, ident in t1["ids"].items():
        if trait not in t1["n"]:
            return ["00_config.R gives no N for outcome %s" % trait]
        want[ident] = _thou(t1["n"][trait])
    for ident, n in t1["med"].items():
        want[ident] = _thou(n)
    seen, bad = set(), []
    for r in rows[1:]:
        cells = [_cellnorm(c) for c in r]
        hit = [c for c in cells if c in want]
        if not hit:
            continue
        ident = hit[0]
        seen.add(ident)
        if want[ident] not in cells:
            bad.append("Table 1 row %r prints no N cell equal to %s" % (cells[:2], want[ident]))
    missing = sorted(set(want) - seen)
    if missing:
        bad.append("Table 1 does not print these deposited dataset identifiers: %s" % missing)
    return bad


def _tol(txt):
    """Half a unit in the last printed place: '0.3625' -> 5e-5, '2.99e-04' -> 5e-7."""
    m = txt.lower().split("e")
    dec = len(m[0].split(".")[1]) if "." in m[0] else 0
    return 0.5 * 10 ** ((int(m[1]) if len(m) > 1 else 0) - dec) * 1.0001


def _num(txt):
    try:
        return float(txt.replace(u"−", "-"))   # float() rejects the typographic minus the tables now print
    except ValueError:
        return None


def check_s12(rows, f=None, side=SIDE, engine=ENGINE, coloc_csv=COLOC_CSV):
    """rows: Supplementary Table S12 as lists of cell text, header first. [] = every "PP.H4 (±1 Mb)" cell equals the
    recovered file at the precision printed (blank only where the wider window could not evaluate the pair) and
    every PP.H4 cell equals the deposited colocalization table."""
    if not rows:
        return ["S12: the table has no rows"]
    h = rows[0]
    need = [u"Gene", u"Outcome", u"PP.H4", u"PP.H4 (±1 Mb)"]
    if not all(x in h for x in need):
        return ["S12 header lacks one of %s: %r" % (need, h)]
    gi, oi, pi, wi = [h.index(x) for x in need]
    want = s12_values(f, side, engine)
    coloc = {}
    with io.open(coloc_csv, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            key = (r["gene"], r["outcome_name"])
            if key in coloc:
                return ["duplicate colocalization row for %s" % (key,)]
            coloc[key] = float(r["PP.H4"])
    body = rows[1:]
    if {(c[gi], c[oi]) for c in body} != set(want) or len(body) != len(want):
        return ["S12 pairs are not the recovered file's %d pairs" % len(want)]
    bad = []
    for c in body:
        key = (c[gi], c[oi])
        v, txt = want[key], c[wi].strip()
        if v is None:
            if txt:
                bad.append("%s-%s: +/-1 Mb prints %r, file: not evaluable" % (key[0], key[1], txt))
        elif _num(txt) is None or abs(_num(txt) - v) > _tol(txt):
            bad.append("%s-%s: +/-1 Mb prints %r, file %r" % (key[0], key[1], txt, v))
        cv, pt = coloc.get(key), c[pi].strip()
        if cv is None or _num(pt) is None or abs(_num(pt) - cv) > _tol(pt):
            bad.append("%s-%s: PP.H4 prints %r, table %r" % (key[0], key[1], pt, cv))
    return bad


def scope_problems(f=None, side=SIDE, engine=ENGINE, manifest=MANIFEST):
    """[] = every recovered file in the engine-root manifest is the source of a claim or a checked table, or is
    declared provenance-only -- never both -- and nothing names a file the engine root does not hold."""
    with io.open(manifest, encoding="utf-8", newline="") as fh:
        rec = {r["path"] for r in csv.DictReader(fh)
               if r["kind"] == "output" and (r["path"].startswith(RP) or r["path"].startswith(LP)
                                             or r["path"].startswith("e_source/"))}
    if not rec:
        return ["the engine-root manifest lists no recovered files, so this scope check would be vacuous"]
    used = {}
    for c in claims(f, side, engine):
        for p in c["src"]:
            used.setdefault(p, []).append(c["id"])
    for name, srcs in sorted(TABLE_SRC.items()):
        for p in srcs:
            used.setdefault(p, []).append(name)
    bad = []
    for p in sorted(rec):
        if p not in used and p not in PROVENANCE_ONLY:
            bad.append("recovered file with no claim and no provenance-only declaration: %s" % p)
    for p in sorted(set(used) | set(PROVENANCE_ONLY)):
        if p not in rec:
            bad.append("a claim or declaration names %s, which the engine root does not hold" % p)
    for p in sorted(set(used) & set(PROVENANCE_ONLY)):
        bad.append("%s is declared provenance-only, but %s is built from it" % (p, ", ".join(used[p])))
    for p, why in sorted(PROVENANCE_ONLY.items()):
        if not why.strip():
            bad.append("the provenance-only declaration for %s carries no reason" % p)
    return bad


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    side = sys.argv[sys.argv.index("--side") + 1] if "--side" in sys.argv else SIDE
    f = facts(side)
    cl = claims(f, side)
    for c in cl:
        print("%-3s %-4s %s%s" % (c["id"], c["where"], c["text"], "   [post-edit]" if c["post"] else ""))
    print("Table 4: %d rows x %d columns; S12 +/-1 Mb values: %d"
          % (len(f["table4"]) - 1, len(f["table4"][0]), len(f["w"]["values"])))
    rc = 0
    if "--scope" in sys.argv:
        sp = scope_problems(f, side)
        print("\nscope: %d problem(s)" % len(sp))
        for b in sp:
            print("  " + b)
        rc = rc or (1 if sp else 0)
    if "--against" in sys.argv:
        lines = io.open(sys.argv[sys.argv.index("--against") + 1], encoding="utf-8").read().splitlines()

        def para(anchor):
            hits = [l.split("\t", 1)[1] for l in lines if "\t" in l and anchor in l.split("\t", 1)[1]]
            if len(hits) != 1:
                raise ValueError("anchor matched %d lines" % len(hits))
            return hits[0]

        bad = check_text(para, only_pre=True, f=f, side=side)
        print("\nagainst %s: %d pre-edit main-text claims, %d problem(s)"
              % (os.path.basename(sys.argv[sys.argv.index("--against") + 1]),
                 sum(1 for c in cl if not c["post"] and c["where"] == "main"), len(bad)))
        for b in bad:
            print("  " + b)
        rc = rc or (1 if bad else 0)
    return rc


if __name__ == "__main__":
    sys.exit(main())