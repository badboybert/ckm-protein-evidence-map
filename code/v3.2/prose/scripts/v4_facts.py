#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""v4_facts.py -- every generation-bound number the V4 rebase prints that the claims registry does not already
build, DERIVED from the corrected-generation files. One implementation, imported by the rebase edits
(rebase_v4_edits.py), so no value in the edited prose is typed.

Sources (all V4 or built here from V4; V3/source_data is never read):
  ER   = V3/v3.2/_engine_root/results                   (engine root, rebuilt from V4 + the fourth delivery)
  A11  = ER/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv (77 columns)
  SUPP = V3/V4/source_data/supp_tables/*.csv
  FIGS = V3/V4/source_data/figures/<dump>

  python v4_facts.py            # print every fact as name = value  [source]
"""
from __future__ import print_function
import csv, io, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
V32 = os.path.dirname(os.path.dirname(HERE))
V3 = os.path.dirname(V32)
ER = os.path.join(V32, "_engine_root", "results")
SIDE = os.path.join(ER, "repro_phase1_hg38")
SUPP = os.path.join(V3, "V4", "source_data", "supp_tables")
FIGS = os.path.join(V3, "V4", "source_data", "figures")
TFR = os.path.join(V3, "V4", "source_data", "tables_from_rds")


class FactError(Exception):
    pass


def rows(path):
    if not os.path.isfile(path):
        raise FactError("missing source: %s" % path)
    with io.open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def T(x):
    return str(x).strip().upper() in ("TRUE", "T", "1")


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def facts():
    F, SRC = {}, {}

    def put(name, value, src):
        F[name] = value
        SRC[name] = src

    a11p = os.path.join(SIDE, "_fig_results", "phase4_annotation", "table_s_pair_evidence_tiered.csv")
    a11 = rows(a11p)
    A11 = "A11 tiered table"
    # A11 is the fourth delivery's file, deposited verbatim; its two locus-census columns are BLANK for the 7 new CKD
    # pairs (Part D). The census is read from locus_competitors.csv (A23) only -- make any read of these columns fail.
    for r in a11:
        for c in ("contested", "attribution_untestable"):
            if c not in r:
                raise FactError("A11 no longer carries %r: revisit the Part D guard" % c)
            del r[c]
    pairs = [(r["gene"], r["outcome"]) for r in a11]
    if len(set(pairs)) != len(pairs):
        raise FactError("A11 has duplicate pairs")
    put("n_pairs", len(a11), A11)
    put("n_genes", len({r["gene"] for r in a11}), A11)
    per = {}
    for r in a11:
        per[r["outcome"]] = per.get(r["outcome"], 0) + 1
    put("pairs_per_outcome", per, A11)

    # --- colocalization entry set (Figure 1 legend, Results P61): every suggestive association evaluable? ---------
    col = rows(os.path.join(TFR, "phase3_coloc__coloc_results_final.csv"))
    ev = [r for r in col if fnum(r["PP.H4"]) is not None and fnum(r["PP.H3"]) is not None]
    put("coloc_rows", len(col), "V4 coloc table")
    put("coloc_evaluable", len(ev), "V4 coloc table (rows with both posteriors)")
    put("coloc_not_evaluable", len(col) - len(ev), "V4 coloc table")

    # --- SuSiE classes (P70) -----------------------------------------------------------------------------------------
    sc = {}
    for r in a11:
        sc[r["susie_class"]] = sc.get(r["susie_class"], 0) + 1
    put("susie_counts", sc, A11 + " susie_class")
    put("susie_not_evaluable_pairs", sorted("%s-%s" % p for p, r in zip(pairs, a11) if r["susie_class"] == "not_evaluable"), A11)
    put("susie_discordant_pairs", sorted("%s-%s" % p for p, r in zip(pairs, a11) if r["susie_class"] == "discordant"), A11)
    put("susie_evaluable", len(a11) - sc.get("not_evaluable", 0), A11)

    # --- tiers, core-complete, Bonferroni (Methods P44, Results P70) --------------------------------------------------
    tiers = {}
    for r in a11:
        tiers[r["evidence_tier"]] = tiers.get(r["evidence_tier"], 0) + 1
    put("tiers", tiers, A11)
    cc = [r for r in a11 if T(r["core_complete"])]
    put("core_complete", len(cc), A11 + " core_complete")
    cc_nonconv = [r for r in cc if r["evidence_tier"] != "Convergent"]
    put("core_complete_not_convergent", len(cc_nonconv), A11)
    noarm = [r for r in cc if not T(r["decode_evaluable_pair"]) and not T(r["aric_evaluable_pair"])]
    put("core_complete_no_evaluable_arm", len(noarm), A11 + " decode/aric_evaluable_pair")
    only_arm = [r for r in noarm if not T(r["has_major_weakness"])]
    put("core_complete_noarm_only_gap", len(only_arm), A11 + " has_major_weakness")
    put("core_complete_noarm_with_weakness", sorted("%s-%s" % (r["gene"], r["outcome"]) for r in noarm if T(r["has_major_weakness"])), A11)
    absent_both = [r for r in only_arm if r["decode_status"] == "NOT_IN_RELEASE" and r["aric_status"] == "NOT_ON_ARIC_PANEL"]
    put("noarm_absent_both", len(absent_both), A11)
    put("noarm_absent_both_example", sorted("%s-%s" % (r["gene"], r["outcome"]) for r in absent_both), A11)
    assayed_weak = [r for r in only_arm if r not in absent_both]
    put("noarm_assayed_weak", sorted("%s-%s" % (r["gene"], r["outcome"]) for r in assayed_weak), A11)
    put("noarm_assayed_weak_decode_discordant",
        sum(1 for r in assayed_weak if r["decode_status"] in ("NOT_REPLICATED",) or
            (r["decode_direction_concordant"] != "" and not T(r["decode_direction_concordant"]))), A11)
    bonf = [r for r in a11 if r["mr_sig_class"] == "Bonferroni"]
    put("bonferroni_pairs", len(bonf), A11 + " mr_sig_class")
    nb = [r for r in a11 if r["mr_sig_class"] != "Bonferroni"]
    put("not_bonferroni", len(nb), A11)
    put("not_bonferroni_by_tier", {t: sum(1 for r in nb if r["evidence_tier"] == t) for t in ("Convergent", "Supported", "Flagged")}, A11)
    weakest = max(a11, key=lambda r: fnum(r["mr_pval"]))
    put("weakest_mr_pair", "%s-%s p=%s" % (weakest["gene"], weakest["outcome"], weakest["mr_pval"]), A11)

    # --- Steiger (P69) -------------------------------------------------------------------------------------------------
    st = {}
    for r in a11:
        st[r["steiger_result"]] = st.get(r["steiger_result"], 0) + 1
    put("steiger", st, A11 + " steiger_result")

    # --- instrument threshold (P69, S18a): pairs whose instrument fails p < 5e-8 ------------------------------------
    ins = {}
    for r in rows(os.path.join(TFR, "phase1_instruments__ukb_cis_instruments.csv")):
        ins.setdefault(r["gene"], []).append(float(r["pval"]))
    lost = [p for p in pairs if min(ins[p[0]]) > 5e-8]
    put("pairs_instrument_lost_5e8", len(lost), "V4 instrument table (gene's min p > 5e-8)")
    put("pairs_instrument_kept_5e8", len(pairs) - len(lost), "V4 instrument table")
    Fs = {}
    for r in rows(os.path.join(TFR, "phase1_instruments__ukb_cis_instruments.csv")):
        Fs.setdefault(r["gene"], []).append(float(r["F_stat"]))
    lostF = [max(Fs[g]) for g, o in lost]
    put("lost_5e8_F_range", "%.1f-%.1f" % (min(lostF), max(lostF)) if lostF else "", "V4 instrument table F_stat")

    # --- deCODE / ARIC / cross-platform (Methods P38, Results P83, Discussion P112) -----------------------------------
    ds = {}
    for r in a11:
        ds[r["decode_status"]] = ds.get(r["decode_status"], 0) + 1
    put("decode_status", ds, A11)
    untest = [r for r in a11 if r["decode_status"] in ("NOT_IN_RELEASE", "NO_EXACT_RSID_MATCH")]
    put("decode_not_testable", len(untest), A11)
    put("decode_absent", ds.get("NOT_IN_RELEASE", 0), A11)
    put("decode_no_variant", ds.get("NO_EXACT_RSID_MATCH", 0), A11)
    ab_egfr = [r for r in a11 if r["decode_status"] == "NOT_IN_RELEASE" and r["outcome"] == "eGFR"]
    put("decode_absent_egfr_pairs", len(ab_egfr), A11)
    put("decode_absent_egfr_genes", len({r["gene"] for r in ab_egfr}), A11)
    testable = [r for r in a11 if r not in untest]
    put("decode_testable", len(testable), A11)
    hc = [r for r in testable if T(r["decode_high_confidence"])]
    put("decode_high_conf", len(hc), A11)
    nhc = [r for r in testable if not T(r["decode_high_confidence"])]
    put("decode_not_high_conf", len(nhc), A11)
    weak = [r for r in nhc if not T(r["decode_evaluable_pair"])]
    put("decode_nhc_weak", len(weak), A11 + " decode_evaluable_pair FALSE")
    put("decode_nhc_weak_split", {k: sum(1 for r in weak if r["decode_status"] == k) for k in ("REPLICATED", "NOT_REPLICATED", "ASSAY_DISCORDANT")}, A11)
    evdisc = [r for r in nhc if T(r["decode_evaluable_pair"]) and r["decode_status"] != "ASSAY_DISCORDANT"]
    put("decode_evaluable_discordant", sorted("%s-%s" % (r["gene"], r["outcome"]) for r in evdisc), A11)
    put("decode_evaluable_discordant_aric", {"%s-%s" % (r["gene"], r["outcome"]): r["aric_status"] + "/" + r["aric_direction_concordant"] for r in evdisc}, A11)
    put("decode_evaluable_discordant_tiers", {"%s-%s" % (r["gene"], r["outcome"]): r["evidence_tier"] for r in evdisc}, A11)
    dconc = [r for r in testable if T(r["decode_direction_concordant"])]
    put("decode_dir_concordant", len(dconc), A11)
    dnot = [r for r in testable if not T(r["decode_direction_concordant"])]
    put("decode_dir_not_concordant", sorted("%s->%s|%s|evaluable=%s" % (r["gene"], r["outcome"], r["decode_status"], r["decode_evaluable_pair"]) for r in dnot), A11)
    # Figure 4c draws the single-assay high-confidence pairs: read the drawn set from the panel's own dump
    dd = rows(os.path.join(FIGS, "Figure_04_replication_sensitivity", "replication_sensitivity__panel1_1.csv"))
    drawn = sorted({r["pair_label"].replace(" → ", "->") for r in dd})
    hcset = {"%s->%s" % (r["gene"], r["outcome"]) for r in hc}
    concset = {"%s->%s" % (r["gene"], r["outcome"]) for r in dconc}
    if not set(drawn) <= hcset & concset or len(dd) != 2 * len(drawn):
        raise FactError("Figure 4c draws a pair that is not high-confidence and concordant, or not one row per cohort")
    put("decode_drawn", len(drawn), "V4 Figure 4 panel-(a) dump, distinct pair_label")
    put("decode_hc_not_drawn", sorted(hcset - set(drawn)), "A11 high-confidence minus the drawn set")
    ar_ev = [r for r in a11 if T(r["aric_evaluable_pair"])]
    put("aric_evaluable", len(ar_ev), A11 + " aric_evaluable_pair")
    put("aric_evaluable_concordant", sum(1 for r in ar_ev if T(r["aric_direction_concordant"])), A11)
    put("aric_evaluable_discordant", sum(1 for r in ar_ev if not T(r["aric_direction_concordant"])), A11)
    px = {r["aric_outcome_proxy_used"] for r in a11}
    if not px <= {"proxy", "direct", ""}:
        raise FactError("aric_outcome_proxy_used vocabulary changed: %s" % sorted(px))
    put("aric_proxy_used", sum(1 for r in a11 if r["aric_outcome_proxy_used"] == "proxy"), A11 + " aric_outcome_proxy_used == 'proxy'")
    xp = {}
    for r in a11:
        xp[r["crossplatform_class"]] = xp.get(r["crossplatform_class"], 0) + 1
    put("crossplatform", xp, A11)

    # --- EAS axis (P73, P112) ------------------------------------------------------------------------------------------
    put("eas_exact", sum(1 for r in a11 if T(r["eas_pqtl_any_exact"])), A11 + " eas_pqtl_any_exact")
    strong = [r for r in a11 if T(r["eas_pqtl_any_strong"])]
    put("eas_strong", len(strong), A11)
    dc = {}
    for r in strong:
        dc[r["eas_pqtl_direction_class"]] = dc.get(r["eas_pqtl_direction_class"], 0) + 1
    put("eas_strong_direction", dc, A11)
    nev = {}
    for r in a11:
        if r["eas_outcome_state"] == "EAS_NOT_EVALUABLE":
            nev[r["eas_state_reason"]] = nev.get(r["eas_state_reason"], 0) + 1
    put("eas_not_evaluable_reasons", nev, A11 + " eas_state_reason")
    ev_eas = [r for r in a11 if T(r["eas_outcome_evaluable"])]
    put("eas_eval_binary", sum(1 for r in ev_eas if r["outcome"] != "eGFR"), A11)
    put("eas_eval_egfr", sum(1 for r in ev_eas if r["outcome"] == "eGFR"), A11)
    est = {}
    for r in a11:
        est[r["eas_outcome_state"]] = est.get(r["eas_outcome_state"], 0) + 1
    if set(est) - {"EAS_COMPATIBLE", "EAS_INCONCLUSIVE", "EAS_DISCORDANT", "EAS_NOT_EVALUABLE"}:
        raise FactError("eas_outcome_state vocabulary changed: %s" % sorted(est))
    if len(ev_eas) != len(a11) - est.get("EAS_NOT_EVALUABLE", 0):
        raise FactError("eas_outcome_evaluable disagrees with eas_outcome_state")
    put("eas_states", est, A11 + " eas_outcome_state")

    # --- Figure 2 / 4b legends: the prior bands at the conservative prior, and the binary stability strip --------------
    p6 = [fnum(r["PP.H4_p12_1e6"]) for r in a11]
    put("prior_bands", dict(robust=sum(1 for x in p6 if x >= 0.80), atten=sum(1 for x in p6 if 0.50 <= x < 0.80),
                            fragile=sum(1 for x in p6 if x < 0.50)), A11 + " PP.H4_p12_1e6")
    stable = sum(1 for r in a11 if T(r["stable"]))
    if stable != F["prior_bands"]["robust"]:
        raise FactError("the Stable strip (%d) is no longer the prior-robust band (%d)" % (stable, F["prior_bands"]["robust"]))
    put("prior_stable", stable, A11 + " stable")
    # Figure 2b's second strip: each Flagged pair's single highest-priority weakness (exploratory_reason)
    fr = {}
    for r in a11:
        if r["evidence_tier"] == "Flagged":
            fr[r["exploratory_reason"]] = fr.get(r["exploratory_reason"], 0) + 1
    if sum(fr.values()) != tiers["Flagged"] or "" in fr:
        raise FactError("a Flagged pair carries no exploratory_reason: %s" % fr)
    if any(r["exploratory_reason"] for r in a11 if r["evidence_tier"] != "Flagged"):
        raise FactError("a non-Flagged pair carries an exploratory_reason")
    put("flag_reasons", fr, A11 + " exploratory_reason over Flagged pairs")

    # --- locus census clauses the registry does not assert (P31) ---------------------------------------------------------
    lc = rows(os.path.join(SIDE, "locus_competitors.csv"))
    cont = [r for r in lc if T(r["contested"])]
    put("contested_all_untestable", all(T(r["attribution_untestable"]) for r in cont), "A23 locus_competitors.csv")
    clean = [r for r in lc if not T(r["attribution_untestable"])]
    put("clean_other", len(clean) - sum(1 for r in clean if int(r["n"]) == 0), "A23")

    # --- prior bands and the dual rule at each prior (P82) ----------------------------------------------------------------
    rob = rows(os.path.join(SUPP, "SuppTable_coloc_robustness_core.csv"))
    put("robustness_core_rows", len(rob), "V4 SuppTable_coloc_robustness_core.csv")
    summ = {r["metric"]: r["value"] for r in rows(os.path.join(SUPP, "coloc_robustness_summary.csv"))}
    if int(summ["core_pairs"]) != len(a11) or int(summ["core_prior_stable"]) != sum(1 for r in a11 if T(r["stable"])):
        raise FactError("coloc_robustness_summary.csv disagrees with A11 on the core pairs or prior-stable count")
    put("core_robust", int(summ["core_robust"]), "V4 coloc_robustness_summary.csv core_robust (the strict joint triangulation)")

    # --- candidate (dashed) layer (P35, P81) -------------------------------------------------------------------------------
    cand = rows(os.path.join(SUPP, "SuppTable_candidate_edge_audit.csv"))
    put("candidate_rows", len(cand), "V4 SuppTable_candidate_edge_audit.csv")
    cs = {}
    for r in cand:
        cs[r["Status"]] = cs.get(r["Status"], 0) + 1
    if set(cs) != {"Drawn", "Dropped", "Held (not drawn)"}:
        raise FactError("candidate audit Status vocabulary changed: %s" % sorted(cs))
    put("candidate_status", cs, "V4 SuppTable_candidate_edge_audit.csv Status")
    held = ["%s–%s" % (r["Gene"], r["Endpoint"]) for r in cand if r["Status"] == "Held (not drawn)"]
    put("candidate_held", ", ".join(held), "V4 SuppTable_candidate_edge_audit.csv")
    put("candidate_drawn_abf_negative", sum(1 for r in cand if r["Status"] == "Drawn" and fnum(r["coloc.abf PP.H4"]) is not None
                                           and fnum(r["coloc.abf PP.H4"]) < 0.02), "V4 SuppTable_candidate_edge_audit.csv (drawn, coloc.abf PP.H4 < 0.02)")
    only = sorted({r["Gene"] for r in cand if r["Status"] == "Drawn"} - {g for g, o in pairs})
    put("candidate_only_genes", only, "V4 SuppTable_candidate_edge_audit.csv drawn genes not in A11")
    put("network_gene_nodes", len({g for g, o in pairs}) + len(only), "A11 genes + candidate-only genes")

    # --- renal scrub (P93) ------------------------------------------------------------------------------------------------
    rs = rows(os.path.join(SUPP, "SuppTable_renal_cystatin_BUN_scrub.csv"))
    if {(r["gene"], r["outcome"]) for r in rs} != {p for p in pairs if p[1] in ("eGFR", "CKD")}:
        raise FactError("the renal scrub table does not cover exactly the renal tiered pairs")
    wc = [r for r in rs if fnum(r["cys_b"]) is not None]
    wb = [r for r in rs if fnum(r["bun_b"]) is not None]
    put("renal_scrub", dict(n=len(rs), with_cys=len(wc), cys_concordant=sum(1 for r in wc if T(r["cys_dir_concordant"])),
                            cys_not_concordant=sorted("%s-%s" % (r["gene"], r["outcome"]) for r in wc if not T(r["cys_dir_concordant"])),
                            with_bun=len(wb), bun_concordant=sum(1 for r in wb if T(r["bun_dir_concordant"]))),
        "V4 SuppTable_renal_cystatin_BUN_scrub.csv")
    # the same scrub at PROTEIN level: a protein colocalizing with both eGFR and CKD has one instrument and so one
    # cystatin-C and one BUN estimate; counting its two pairs would count one test twice
    byg = {}
    for r in rs:
        byg.setdefault(r["gene"], []).append(r)
    KEYS = ("SNP", "cys_b", "cys_p", "cys_dir_concordant", "bun_b", "bun_p", "bun_dir_concordant")
    for gname, rr in byg.items():
        if any(tuple(x[k] for k in KEYS) != tuple(rr[0][k] for k in KEYS) for x in rr):
            raise FactError("%s: its renal pairs carry different cystatin-C/BUN rows; protein-level counting is wrong" % gname)
    pr = [rr[0] for rr in byg.values()]
    dn = [r for r in pr if r["gene"] == "DNER"]
    if len(dn) != 1 or fnum(dn[0]["cys_b"]) != 0.0:
        raise FactError("DNER's cystatin-C estimate is no longer the null the renal sentence describes")
    put("dner_cys_p", float(dn[0]["cys_p"]), "V4 SuppTable_renal_cystatin_BUN_scrub.csv DNER cys_p")
    pc = [r for r in pr if fnum(r["cys_b"]) is not None]
    pb = [r for r in pr if fnum(r["bun_b"]) is not None]
    put("renal_scrub_proteins", dict(
        n=len(pr), shared=sorted(g_ for g_, rr in byg.items() if len(rr) > 1),
        with_cys=len(pc), cys_concordant=sum(1 for r in pc if T(r["cys_dir_concordant"])),
        cys_not_concordant=sorted(r["gene"] for r in pc if not T(r["cys_dir_concordant"])),
        with_bun=len(pb), bun_concordant=sum(1 for r in pb if T(r["bun_dir_concordant"]))),
        "V4 SuppTable_renal_cystatin_BUN_scrub.csv, one row per protein (shared rows asserted identical)")
    # the kidney-expression comparison printed on Supplementary Figure S20c (two-sided exact Mann-Whitney; R's
    # wilcox.test on the same dump gives the same W = 50, p = 0.9729)
    from scipy.stats import mannwhitneyu
    kd = rows(os.path.join(FIGS, "S21_sfig_cystatin_scrub", "sfig_cystatin_scrub__panel2_2.csv"))
    cls = sorted({r["class"] for r in kd})
    if len(cls) != 2:
        raise FactError("the kidney-expression dump no longer has two classes: %s" % cls)
    xs = [[float(r["kid_tpm"]) for r in kd if r["class"] == c] for c in cls]
    put("kidney_wilcoxon", dict(p=float(mannwhitneyu(xs[0], xs[1], alternative="two-sided", method="exact").pvalue),
                                n=[len(x) for x in xs]), "V4 S21 dump sfig_cystatin_scrub__panel2_2.csv")

    # --- eGFR-only genes (P92) ---------------------------------------------------------------------------------------------
    by_gene = {}
    for g, o in pairs:
        by_gene.setdefault(g, set()).add(o)
    put("genes_only_egfr", sorted(g for g, s in by_gene.items() if s == {"eGFR"}), A11)
    put("genes_three_plus_endpoints_raw", sorted("%s(%d)" % (g, len(s)) for g, s in by_gene.items() if len(s) >= 3), A11)
    unit = lambda o: "Renal" if o in ("eGFR", "CKD") else o
    put("genes_three_plus_units", sorted("%s(%d)" % (g, len({unit(o) for o in s})) for g, s in by_gene.items()
                                          if len({unit(o) for o in s}) >= 3), A11)
    put("renal_pairs", sum(1 for g, o in pairs if o in ("eGFR", "CKD")), A11)

    # --- Figure 5b / Results: proteins with a targeted-Bonferroni protein -> intermediate association ---------------------
    genes = {g for g, o in pairs}
    alpha = 0.05 / (len(genes) * 6)
    mi = rows(os.path.join(TFR, "phase2b_intermediate__mr_intermediate_all.csv"))
    ints = {r["outcome_name"] for r in mi}
    if len(ints) != 6:
        raise FactError("expected 6 intermediates in the protein->intermediate table, found %d" % len(ints))
    hit = {r["gene"] for r in mi if r["gene"] in genes and fnum(r["pval"]) is not None and float(r["pval"]) < alpha}
    put("targeted_alpha", alpha, "0.05 / (n_genes x 6)")
    put("tripartite_proteins", len(hit), "V4 phase2b_intermediate__mr_intermediate_all.csv, p < 0.05/(n_genes x 6)")
    # Figure 5b as drawn: the dumps are asserted against the derivations above before they are read
    F7 = os.path.join(FIGS, "Figure_07_tripartite_pathway", "tripartite_pathway__panel__layer%s.csv")
    nodes = {r["name"] for r in rows(F7 % "7_GeomPoint")}
    if nodes != hit:
        raise FactError("Figure 5b draws proteins %s, not the targeted-Bonferroni set" % sorted(nodes ^ hit))
    arcs = rows(F7 % "3_GeomCurve")
    if len(arcs) != sum(1 for g_, o in pairs if g_ in hit):
        raise FactError("Figure 5b's context arcs are not the tiered pairs of the drawn proteins")
    edges = rows(F7 % "4_GeomSegment") + rows(F7 % "5_GeomSegment")
    lp = [float(r["lp"]) for r in edges]
    put("f5b", dict(arcs=len(arcs), edges=len(edges), capped=sum(1 for v in lp if v >= 20), cap=max(lp)),
        "V4 Figure_07 dumps (arcs asserted = tiered pairs of the %d drawn proteins)" % len(hit))
    # Figure 5a: the matrix blocks, column-wise, ceil(n/2) per block, in the dumped order
    mx = rows(os.path.join(FIGS, "Figure_06_prioritization_matrix", "prioritization_matrix__panel1_1.csv"))
    gl = [r["gene"] for r in mx]
    if set(gl) != genes or len(gl) != len(genes):
        raise FactError("the matrix dump is not the tiered gene set")
    k = (len(gl) + 1) // 2
    put("f5a_blocks", dict(k=k, n=len(gl), l0=gl[0], l1=gl[k - 1], r0=gl[k], r1=gl[-1]),
        "V4 Figure_06 dump order, ceil(n/2) per block")

    # --- cis windows clipped short of +/-500 kb (Figure 1 legend) --------------------------------------------------------
    gi = {}
    for r in rows(os.path.join(TFR, "phase1_instruments__ukb_cis_instruments.csv")):
        gi.setdefault(r["gene"], r)
    lo = [r for r in gi.values() if int(r["gene_start"]) - int(r["cis_start"]) < 500000]
    hi = [r for r in gi.values() if int(r["cis_end"]) - int(r["gene_end"]) < 500000]
    ends = {}
    for r in gi.values():
        ends[r["gene_chrom"]] = max(ends.get(r["gene_chrom"], 0), int(r["gene_end"]))
    # a start clip is a chromosome end only at position 1; an end clip only where the window stops at one shared
    # coordinate per chromosome that lies beyond every gene on it (the chromosome length) -- anything else would be a
    # neighbouring-gene boundary, which the legend no longer claims
    if any(int(r["cis_start"]) != 1 for r in lo):
        raise FactError("a cis window is clipped at its start somewhere other than position 1")
    per_chr = {}
    for r in hi:
        per_chr.setdefault(r["gene_chrom"], set()).add(int(r["cis_end"]))
    if any(len(v) != 1 or min(v) < ends[c] for c, v in per_chr.items()):
        raise FactError("a cis window is clipped at its end short of the chromosome end")
    put("cis_clipped_chrom_end", len({r["gene"] for r in lo + hi}), "V4 instrument table cis_start/cis_end vs gene bounds")

    # --- instruments per pair (S18a legend) ---------------------------------------------------------------------------
    keys = set(pairs)
    meth = {}
    for r in rows(os.path.join(TFR, "phase2_mr_screen__mr_all_results.csv")):
        k = (r["gene"], r["outcome_name"])
        if k in keys and r["method"] in ("Wald ratio", "Inverse variance weighted"):
            meth[k] = r["method"]
    if set(meth) != keys:
        raise FactError("not every tiered pair has a primary MR row")
    put("single_instrument_pairs", sum(1 for m in meth.values() if m == "Wald ratio"), "V4 MR table, primary method == Wald ratio")

    # --- lead-variant coding consequences (Discussion flagged pairs, Limitations, S12 note) ----------------------------
    vep = rows(os.path.join(V3, "V4_local", "source_data", "annotation", "lead_snp_consequences.csv"))
    if {(r["gene"], r["outcome"]) for r in vep} != keys:
        raise FactError("the lead-SNP consequence table does not cover exactly the tiered pairs")
    mis = [r for r in vep if T(r["missense"])]
    put("missense_led", dict(pairs=len(mis), genes=len({r["gene"] for r in mis}),
                             canonical=sum(1 for r in vep if T(r["missense_canonical"])),
                             truncating=sum(1 for r in vep if T(r["truncating"])), release=vep[0]["ensembl_release"]),
        "V4_local annotation/lead_snp_consequences.csv (build_v4_local_vep.py)")
    # LY96-T2D, the flagged pair the Discussion names in place of IDI2-CKD: every clause of its sentence, asserted
    ly = [r for r in a11 if (r["gene"], r["outcome"]) == ("LY96", "T2D")]
    lv = [r for r in vep if (r["gene"], r["outcome"]) == ("LY96", "T2D")]
    if len(ly) != 1 or len(lv) != 1:
        raise FactError("LY96-T2D is no longer a tiered pair")
    ly, lv = ly[0], lv[0]
    weak = {"prior-fragile": fnum(ly["PP.H4_p12_1e6"]) < 0.50,
            "SuSiE discordant": ly["susie_class"] == "discordant",
            "non-forward Steiger": ly["steiger_result"] != "FORWARD_CONFIRMED",
            "assay discordant": T(ly["flag_assay_discordant"]),
            "replication counter-evidence": ly["decode_status"] in ("NOT_REPLICATED", "ASSAY_DISCORDANT")
                                            or ly["aric_status"] in ("NOT_REPLICATED",)}
    aa = lv["amino_acids_canonical"]
    m = re.fullmatch(r"([A-Z])/([A-Z])(\d+)", aa)
    checks = [ly["evidence_tier"] == "Flagged", [k for k, v in weak.items() if v] == ["prior-fragile"],
              ly["mr_sig_class"] != "Bonferroni", ly["susie_class"] == "inconclusive", ly["lead_snps"] == lv["lead_snp"],
              T(lv["missense_canonical"]), m is not None, len(lv["allele_string"].split("/")) == 3,
              lv["most_severe_site"] == "stop_gained", "stop_gained" not in lv["consequence"],
              "stop_gained" in lv["other_allele_consequences"],
              # "neither cross-platform cohort assays it, and its colocalization does not hold at the ±1 Mb window"
              ly["decode_status"] == "NOT_IN_RELEASE", ly["aric_status"] == "NOT_ON_ARIC_PANEL",
              ly["coloc_1mb_status"] == "ok", not T(ly["is_coloc_pass_1mb"])]
    if not all(checks):
        raise FactError("the LY96-T2D sentence no longer holds: %s" % checks)
    # the analysed allele is the MAJOR allele (fix-round-2 verifier): 'missense' is relative to the GRCh38 reference, so
    # the sentence gives the frequency; three-letter HGVS to match p.Arg46Leu elsewhere in the paper
    aa3 = dict(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu", G="Gly", H="His", I="Ile", L="Leu",
               K="Lys", M="Met", F="Phe", P="Pro", S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val")
    ins = [r for r in rows(os.path.join(TFR, "phase1_instruments__ukb_cis_instruments.csv"))
           if r["SNP"] == lv["lead_snp"] and r["gene"] == "LY96"]
    if len(ins) != 1 or ins[0]["effect_allele"] != lv["analysed_allele"] or lv["grch38_ref"] == lv["analysed_allele"]:
        raise FactError("the LY96 instrument row does not pin the analysed allele as a non-reference allele")
    put("ly96", dict(snp=lv["lead_snp"], change="%s%s%s" % (m.group(1), m.group(3), m.group(2)), w1=float(ly["PP.H4_1mb"]),
                     hgvs="p.%s%s%s" % (aa3[m.group(1)], m.group(3), aa3[m.group(2)]), eaf=float(ins[0]["eaf"])),
        "A11 LY96-T2D row + V4_local lead_snp_consequences.csv + V4 phase1 instrument table (eaf)")
    return F, SRC


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    F, SRC = facts()
    for k in F:
        print("%-38s = %s   [%s]" % (k, json.dumps(F[k], ensure_ascii=False) if not isinstance(F[k], (int, str)) else F[k], SRC[k]))
