"""w2_08_assemble.py -- results.tsv and branches.tsv for section 2, from the component outputs of this folder only.
Computed from the REDACTED pre-declaration (rev 4.2) and amendment 1 before the full text is opened (brief rule 10).

Branch rules implemented (redacted pre-declaration section 2, Principle 0, Rules 3(i), 3(ii), 3(iv), 4):
  s2 (a) after mutual conditioning the ratios have opposite signs in both resources AND heterogeneity p < 0.05 in both
  s2 (b) both conditional ratios have the discovery sign in both resources AND heterogeneity p >= 0.05 in both
  s2 (c) anything else
  heterogeneity p = p_aware_2nd (Rule 5) in FinnGen CHD (primary, family entry) and ieu-a-7, 503 panel.
FIX1 (2026-10-08, Stage 1 fixer; see FIXLOG.md): Rule 3(ii) on the named primary cell now also evaluates that cell's
conditional colocalizations (defect 1); the Rule 3(i) first-mention flag is determined from the existing class (defect 2).
"""
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

COLS = ["quantity_id", "prespec_ref", "estimand", "resource", "panel", "cov", "se_order", "value", "se", "ci_lo", "ci_hi", "p", "n", "unit", "notes"]
SIGS = {"A": "rs73045269", "B": "rs12983775"}
TRN = {"P": "protein", "CHD_FG": "finngen_chd", "CAD_C4D": "ieu_a_7"}


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fmt_p(p):
    return "%.2g" % p


# FIX1 (FIXLOG.md, defect 1): Rule 3(ii) covers "conditional, ratio, GLS and conditional-colocalization results" and is
# evaluated on the named primary cell (TGFB1: FinnGen CHD conditional). The conditional colocalizations at A and B in that
# cell are therefore evaluated beside its ratio test. Conclusion for a colocalization = side of PP.H4 0.80; robust iff the
# same side on 503 and 489 AND >= 90% of the 5,000 503-panel draws on the full-panel side.
def coloc_primary_cell(v, cbs, rs="finngen_chd"):
    out = []
    for lab, sig in SIGS.items():
        h503, h489 = v("S2.coloc_clause.%s.%s.503.PPH4" % (sig, rs)), v("S2.coloc_clause.%s.%s.489.PPH4" % (sig, rs))
        sh = cbs.get((sig, rs), float("nan"))          # missing bootstrap -> NaN -> not robust (never silently robust)
        same = (h503 >= 0.80) == (h489 >= 0.80)
        out.append({"label": lab, "sig": sig, "h503": h503, "h489": h489, "share": sh, "same_panels": same, "robust": same and sh >= 0.90})
    return out


# FIX1 (FIXLOG.md, defect 2): Rule 3(i) "A pair whose class rests on inadmissible evidence is flagged at first mention".
# The existing class and the pipeline fit are READ (raw_input tables; md5 asserted via RESOURCES), never recomputed or rebuilt.
# The LD source of the pipeline SuSiE grid is read from the engine config (role engine; md5 asserted).
def rule3i_existing_class(gene="TGFB1", outcome="CAD", grid_outcome="FinnGen_R10_I9_CHD"):
    tp, gp, cp = W.res("TAB_tiered81"), W.res("TAB_susie_grid"), W.res("ENG_engine_config_00_config.R")
    with open(tp, encoding="utf-8") as fh:
        tier = list(csv.DictReader(fh))
    t = [x for x in tier if x["gene"] == gene and x["outcome"] == outcome]
    assert len(t) == 1, ("TAB_tiered81 rows for %s-%s" % (gene, outcome), len(t))
    t = t[0]
    with open(gp, encoding="utf-8") as fh:
        grid = list(csv.DictReader(fh))
    g = [x for x in grid if x["gene"] == gene]
    assert len(g) == 4 and {x["outcome_id"] for x in g} == {grid_outcome}, ("TAB_susie_grid rows for %s" % gene, [(x["outcome_id"], x["setting"]) for x in g])
    with open(cp, encoding="utf-8") as fh:
        cfg = fh.read().splitlines()
    hdr = [i + 1 for i, line in enumerate(cfg) if "Local 1KG EUR LD Reference (Phase 3-C SuSiE)" in line]
    ldd = [i + 1 for i, line in enumerate(cfg) if line.startswith("LD_EUR_DIR") and "1kg_eur" in line]
    ref_ld = len(hdr) == 1 and len(ldd) == 1
    if not ref_ld:
        W.log("ASSERT_FAILED 00_config.R: Phase 3-C SuSiE 1KG EUR LD-reference lines not found (%s, %s)" % (hdr, ldd))
    # descriptive only (tiers are never rebuilt): pairs matching every other criterion column, without the SuSiE criterion
    keys = [k for k in t if (k.startswith("crit_") or k in ("has_major_weakness", "core_complete", "essential_core_missing", "replication_depth",
                                                             "all_strong_criteria")) and k != "crit_susie_supportive"]
    twins = [x for x in tier if all(x[k] == t[k] for k in keys) and x["crit_susie_supportive"] != "TRUE"]
    return {"tier_row": t, "grid": g, "ref_ld": ref_ld, "cfg_lines": (hdr[0] if hdr else 0, ldd[0] if ldd else 0), "twins": twins,
            "md5": {k: W.res_entry(k)["md5"] for k in ("TAB_tiered81", "TAB_susie_grid", "ENG_engine_config_00_config.R")}}


def rule3i_flag(adm, ev):
    """True when the pair's existing SuSiE class comes from a reference-LD pipeline fit that is inadmissible on the proxy s."""
    t, g = ev["tier_row"], ev["grid"]
    fitted = t["susie_class"] not in ("", "NA", "not_evaluable") and all(x["susie_status"] == "ok" for x in g)
    return (not adm) and fitted and ev["ref_ld"]


def main():
    W.log("==== w2_08_assemble ====")
    rows = []
    for f in ("results_cond.tsv", "results_annot.tsv", "results_twobytwo.tsv"):
        rows += W.read_tsv(W.OUT + "/" + f)
    # het rows: value = p (the id names a p value); the difference d goes to its own .diff row
    out = []
    for r in rows:
        q = r["quantity_id"]
        if ".het." in q and q.split(".")[-1] in ("p_aware_2nd", "p_ign", "p_aware_1st", "p_cojo") and "boot" not in q:
            d, se, lo, hi = r["value"], r["se"], r["ci_lo"], r["ci_hi"]
            r = dict(r, value=r["p"], se="", ci_lo="", ci_hi="",
                     notes=(r["notes"] + "; " if r["notes"] else "") + "A-B difference d = %.6g (se %.6g, 95%% CI %.6g to %.6g)" % (fl(d), fl(se), fl(lo), fl(hi)))
            if q.endswith("p_aware_2nd"):
                out.append(dict(r, quantity_id=q.replace("p_aware_2nd", "diff"), estimand="A-B difference of conditional ratios rho(rs73045269) - rho(rs12983775)",
                                value=d, se=se, ci_lo=lo, ci_hi=hi, p=r["p"], notes="aware, second order (Rule 5 primary SE)"))
        out.append(r)
    rows = out
    # ---------------- SuSiE (Rule 3(i))
    st = W.read_tsv(W.OUT + "/susie_traits.tsv")
    cs = W.read_tsv(W.OUT + "/susie_cs.tsv")
    assert len(st) == 12, len(st)
    smax = {}
    for t in st:
        trait = TRN[t["trait"]]
        css = [c for c in cs if c["set"] == t["set"] and c["trait"] == t["trait"]]
        desc = "; ".join("%s lead %s (pip %.3f, size %s, min|r| %.3f%s%s%s)" % (
            c["cs"], c["lead"], fl(c["lead_pip"]), c["size"], fl(c["purity_min"]), ", has rs73045269" if c["has_A"] == "TRUE" else "",
            ", has rs12983775" if c["has_B"] == "TRUE" else "", ", has rs1800470" if c["has_rs1800470"] == "TRUE" else "") for c in css)
        rows.append({"quantity_id": "S2.susie.%s.%s.s" % (trait, t["set"]), "prespec_ref": "s2; Rule 3(i)",
                     "estimand": "susieR::estimate_s_rss (null-mle) LD-mismatch s", "resource": trait, "panel": t["panel"],
                     "value": t["s"], "n": t["n"], "unit": "s",
                     "notes": "nsnp %s; susie_rss L=10 converged %s (niter %s), warnings: %s; local s = proxy for any pipeline fit" % (
                         t["nsnp"], t["converged"], t["niter"], t["warnings"] or "none")})
        rows.append({"quantity_id": "S2.susie.%s.%s.cs" % (trait, t["set"]), "prespec_ref": "s2; Rule 3(i)",
                     "estimand": "SuSiE-RSS 95% credible sets (L = 10, min_abs_corr 0.5)", "resource": trait, "panel": t["panel"],
                     "value": "%s credible sets" % t["n_cs"], "n": t["nsnp"], "notes": desc or "none"})
        smax.setdefault(trait, []).append(fl(t["s"]))
    # Stage 2 step 2 (user decision 2026-10-09): s_max is over ALL SEVEN distinct fits (all variant sets, both panels): the four
    # fits above + p2np@503, susie@489, p1fr@489 (w2_04b_fitset7.R; p1np@489 = p2np@489). The four-fit value is kept as a sensitivity.
    xs = W.read_tsv(W.OUT + "/susie_s_extra.tsv")
    assert len(xs) == 9, len(xs)
    FIT7_NOTE = {"p2np_503panel": "p2np SNP set on the 503 panel", "susie_489panel": "susie set restricted to SNPs valid on the 489 panel",
                 "p1fr_489panel": "p1fr set restricted to SNPs valid on the 489 panel"}
    smax4 = {tr: max(v) for tr, v in smax.items()}
    for x in xs:
        trait = TRN[x["trait"]]
        rows.append({"quantity_id": "S2.susie.%s.%s.s" % (trait, x["set"]), "prespec_ref": "s2; Rule 3(i); Stage 2 seven-fit set",
                     "estimand": "susieR::estimate_s_rss (null-mle) LD-mismatch s", "resource": trait, "panel": x["panel"],
                     "value": x["s"], "n": x["n"], "unit": "s",
                     "notes": "nsnp %s; %s; s only (no SuSiE fit); added at Stage 2 for the seven-fit s_max; local s = proxy for any pipeline fit" % (x["nsnp"], FIT7_NOTE[x["set"]])})
        smax.setdefault(trait, []).append(fl(x["s"]))
    for trait, v in smax.items():
        assert len(v) == 7, (trait, len(v))
        rows.append({"quantity_id": "S2.susie.%s.s_max" % trait, "prespec_ref": "Rule 3(i)",
                     "estimand": "max s over all seven distinct fits (all variant sets, both panels; p1np on 489 = p2np on 489)",
                     "resource": trait, "panel": "503+489", "value": max(v), "unit": "s"})
        rows.append({"quantity_id": "S2.susie.%s.s_max_4fit" % trait, "prespec_ref": "Rule 3(i) sensitivity",
                     "estimand": "sensitivity: max s over the earlier four-fit set (susie, p1fr, p1np on 503; p2np on 489)",
                     "resource": trait, "panel": "503+489", "value": smax4[trait], "unit": "s"})
    # ---------------- Rule 3(iv) inputs
    cc = W.read_tsv(W.OUT + "/coloc_cond.tsv")
    res = {r["quantity_id"]: r for r in rows}
    for c in cc:
        sig = SIGS[c["signal"]]
        pfx = "S2" if c["panel"] == "503" else "S2.p489"
        zrow = res["%s.cond.%s.%s.beta" % (pfx, c["resource"], sig)]
        z = fl(zrow["value"]) / fl(res["%s.cond.%s.%s.se" % (pfx, c["resource"], sig)]["value"])
        rows.append({"quantity_id": "S2.coloc_clause.%s.%s.%s.z" % (sig, c["resource"], c["panel"]), "prespec_ref": "Rule 3(iv)",
                     "estimand": "CAD conditional association z at %s (given %s)" % (sig, c["conditioned_on"]), "resource": c["resource"],
                     "panel": c["panel"], "value": z, "unit": "z", "notes": "z-space joint fit; r(A,B) from the %s panel" % c["panel"]})
        rows.append({"quantity_id": "S2.coloc_clause.%s.%s.%s.PPH4" % (sig, c["resource"], c["panel"]), "prespec_ref": "Rule 3(iv)",
                     "estimand": "conditional PP.H4 protein x CAD at %s (both conditioned on %s), coloc.abf p1 = p2 = 1e-4, p12 = 1e-5" % (sig, c["conditioned_on"]),
                     "resource": c["resource"], "panel": c["panel"], "value": fl(c["H4"]), "n": c["nsnp"],
                     "notes": "PP.H3 %.4f; top SNP %s; R coloc 5.2.3 (protein quant sdY 1; CAD cc); engine convention (both quant sdY 1) PP.H4 %.4f; window = %s set" % (
                         fl(c["H3"]), c["top_snp_H4"], fl(c["H4_engine_quant_sdY1"]), "p1np" if c["panel"] == "503" else "p2np")})
    cbs = {}
    if os.path.exists(W.OUT + "/coloc_boot_summary.tsv"):
        for b in W.read_tsv(W.OUT + "/coloc_boot_summary.tsv"):
            _, sg, tr1, tr2 = b["quantity"].split("_", 3) if b["quantity"].count("_") == 3 else (None,) * 4
            sig = SIGS[sg]
            rsrc = TRN[tr1 + "_" + tr2]
            cbs[(sig, rsrc)] = fl(b["share_same_side"])
            rows.append({"quantity_id": "S2.coloc_clause.boot.%s.%s.share" % (sig, rsrc), "prespec_ref": "Rule 3(ii) (conditional colocalization)",
                         "estimand": "share of 5,000 503-panel bootstrap draws with conditional PP.H4 on the full-panel side of 0.80", "resource": rsrc,
                         "panel": "503", "value": fl(b["share_same_side"]), "n": b["n_draws"],
                         "notes": "full 503 PP.H4 %.4f (%s); draw PP.H4 2.5/50/97.5%%: %.3f/%.3f/%.3f; MC SE %.4f; same draws as the ratio bootstrap (seed 20261007)" % (
                             fl(b["full_503"]), b["full_side"], fl(b["h4_q025"]), fl(b["h4_q500"]), fl(b["h4_q975"]), fl(b["mc_se"]))})
    else:
        W.log("coloc_boot_summary.tsv missing: conditional-colocalization bootstrap rows NOT RUN")
    # ---------------- Rule 3(i) first-mention flag inputs (FIX1, defect 2): read, not computed
    ev = rule3i_existing_class()
    t, g, m5 = ev["tier_row"], ev["grid"], ev["md5"]
    gs = {x["setting"]: x for x in g}
    rows.append({"quantity_id": "S2.rule3i.class.TGFB1_CAD.susie_class", "prespec_ref": "Rule 3(i) (first-mention flag)",
                 "estimand": "existing SuSiE class of TGFB1-CAD in the 81-pair tiered table (read; classes and tiers are not rebuilt)",
                 "resource": "TAB_tiered81", "value": t["susie_class"],
                 "notes": "read from TAB_tiered81 (raw_input, md5 %s), row gene TGFB1 / outcome CAD: susie_class %s, crit_susie_supportive %s, evidence_tier %s, "
                          "pp_h4 %s; descriptive: %d of 81 pairs match every other criterion column with crit_susie_supportive not TRUE, tiers %s, so the tier "
                          "does not appear to rest on the SuSiE criterion (tier rule source not local; tier not rebuilt)" % (
                              m5["TAB_tiered81"], t["susie_class"], t["crit_susie_supportive"], t["evidence_tier"], t["pp_h4"], len(ev["twins"]),
                              "/".join(sorted({x["evidence_tier"] for x in ev["twins"]})) or "none")})
    rows.append({"quantity_id": "S2.rule3i.pipeline_fit.TGFB1_FinnGen_CHD.PPH4_susie", "prespec_ref": "Rule 3(i) (first-mention flag)",
                 "estimand": "pipeline coloc.susie PP.H4 for TGFB1 x FinnGen R10 I9_CHD, strict setting (read; the fit behind the SuSiE class)",
                 "resource": "TAB_susie_grid", "panel": "pipeline 1KG EUR reference", "value": fl(gs["strict"]["PP.H4_susie"]), "n": gs["strict"]["n_snps"],
                 "notes": "read from TAB_susie_grid (raw_input, md5 %s): %s; classification %s in %d of %d settings; converged exp/out %s" % (
                     m5["TAB_susie_grid"],
                     "; ".join("%s L=%s PP.H4_susie %.3f" % (s, gs[s]["L"], fl(gs[s]["PP.H4_susie"])) for s in ("strict", "relaxed", "exploratory", "robustness")),
                     t["susie_class"], sum(x["classification"] == t["susie_class"] for x in g), len(g),
                     "/".join(sorted({"%s/%s" % (x["converged_exp"], x["converged_out"]) for x in g})))})
    rows.append({"quantity_id": "S2.rule3i.pipeline_fit.ld_source", "prespec_ref": "Rule 3(i) (first-mention flag)",
                 "estimand": "LD source of the pipeline SuSiE grid (Phase 3-C)", "resource": "ENG_engine_config_00_config.R",
                 "value": "1KG EUR reference LD" if ev["ref_ld"] else "NOT FOUND",
                 "notes": "engine config (role engine, read to understand; md5 %s) l.%d header 'Local 1KG EUR LD Reference (Phase 3-C SuSiE)', l.%d LD_EUR_DIR "
                          ".../ld_reference/1kg_eur; the SuSiE engine source (15_susie_coloc.R) is not local (RESOURCES ENG_V4_engine_15_susie_coloc not_found), "
                          "so the reference-LD attribution rests on the config" % ((m5["ENG_engine_config_00_config.R"],) + ev["cfg_lines"])})
    # ---------------- results.tsv
    for r in rows:
        for c in COLS:
            r.setdefault(c, "")
    ids = [r["quantity_id"] for r in rows]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        W.log("ASSERT_FAILED duplicate quantity ids: %s" % dup[:10])
    W.write_tsv(W.OUT + "/results.tsv", COLS, rows)
    R = {r["quantity_id"]: r for r in rows}
    v = lambda q: fl(R[q]["value"])  # noqa: E731
    # ---------------- branches
    br = []
    pf, pc = v("S2.het.finngen_chd.p_aware_2nd"), v("S2.het.ieu_a_7.p_aware_2nd")
    rA = {k: v("S2.ratio.%s.rs73045269.value" % k) for k in ("finngen_chd", "ieu_a_7")}
    rB = {k: v("S2.ratio.%s.rs12983775.value" % k) for k in ("finngen_chd", "ieu_a_7")}
    ds = v("S2.discovery_sign")
    opp = all(math.copysign(1, rA[k]) != math.copysign(1, rB[k]) for k in rA)
    disc = all(math.copysign(1, rA[k]) == ds and math.copysign(1, rB[k]) == ds for k in rA)
    if opp and pf < 0.05 and pc < 0.05:
        branch, tmpl = "(a)", "a reproducible CAD locus whose protein route is instrument-dependent (p %s); TGFB1 leaves the follow-up and candidate sentences" % fmt_p(pf)
    elif disc and pf >= 0.05 and pc >= 0.05:
        branch, tmpl = "(b)", "the marginal sign difference did not persist after mutual conditioning (p %s); TGFB1 stays a candidate, with that caveat" % fmt_p(pf)
    else:
        branch, tmpl = "(c)", "the two-instrument analysis did not resolve the protein route; TGFB1 leaves the candidate lists"
    # Rule 3(ii) on the named primary cell (FinnGen CHD conditional)
    p489 = v("S2.p489.het.finngen_chd.p_aware_2nd")
    d503, d489 = v("S2.het.finngen_chd.diff"), v("S2.p489.het.finngen_chd.diff")
    share = v("S2.het.boot.finngen_chd.share")
    same_panels = (pf < 0.05) == (p489 < 0.05) and math.copysign(1, d503) == math.copysign(1, d489)
    ratio_robust = same_panels and share >= 0.90
    cell = coloc_primary_cell(v, cbs)                                  # FIX1 defect 1: the cell's conditional colocalizations
    robust = ratio_robust and all(c["robust"] for c in cell)
    why = ([] if ratio_robust else ["ratio test not robust (share %.3f)" % share]) + [
        "conditional colocalization at %s %s not robust (PP.H4 503 %.3f / 489 %.3f; share %.3f on the full-panel side of 0.80, threshold 0.90): whether the "
        "Rule 3(iv) clause fires at %s is itself sensitive to the reference LD" % (c["label"], c["sig"], c["h503"], c["h489"], c["share"], c["label"])
        for c in cell if not c["robust"]]
    # Rule 3(iv): fires when |z| >= 1.96 but PP.H4 < 0.80 on either panel
    fires = []
    for sig in SIGS.values():
        for rs in ("finngen_chd", "ieu_a_7"):
            zz = [v("S2.coloc_clause.%s.%s.%s.z" % (sig, rs, pn)) for pn in ("503", "489")]
            hh = [v("S2.coloc_clause.%s.%s.%s.PPH4" % (sig, rs, pn)) for pn in ("503", "489")]
            f = any(abs(z) >= 1.96 for z in zz) and any(h < 0.80 for h in hh)
            fires.append((sig, rs, f, zz, hh))
    # Rule 3(i) admissibility and first-mention flag (computed here so the flag pointer travels with the slot sentence)
    sP, sF, sC = v("S2.susie.protein.s_max"), v("S2.susie.finngen_chd.s_max"), v("S2.susie.ieu_a_7.s_max")
    adm = sP < 0.10 and sF < 0.10
    flag = rule3i_flag(adm, ev)                                        # FIX1 defect 2
    quals = ["nominal: PENDING_HOLM (Rule 4 family member TGFB1-CAD FinnGen CHD; W9 resolves)"] if branch == "(a)" else []
    if not robust:
        quals.append("sensitive to the reference LD (Rule 3(ii) on the named primary cell, FinnGen CHD conditional: %s)" % "; ".join(why))
    for sig, rs, f, zz, hh in fires:
        if f and rs == "finngen_chd":
            quals.append("the CAD association at signal %s did not colocalize with the protein (PP.H4 %.2f), so the ratio there is not anchored to the protein signal" % (sig, min(hh)))
    quals.append("any sentence about rs12983775 adds: it tags a protein-altering allele, so an epitope or assay effect is not excluded")
    quals.append("TGFB1 is never described as an ancestry-specific reversal")
    if flag:
        quals.append("Rule 3(i) first-mention flag (see the Rule 3(i) row): TGFB1-CAD is flagged at first mention, and the flag travels with any statement of "
                     "its SuSiE class or tier: that class rests on inadmissible reference-LD SuSiE, labelled 'proxy s (locally computed; pipeline fit)'")
    inputs = ["S2.ratio.finngen_chd.rs73045269.value", "S2.ratio.finngen_chd.rs12983775.value", "S2.ratio.ieu_a_7.rs73045269.value",
              "S2.ratio.ieu_a_7.rs12983775.value", "S2.het.finngen_chd.p_aware_2nd", "S2.het.ieu_a_7.p_aware_2nd", "S2.discovery_sign"]
    br.append({"rule_ref": "prespec s2 wording rule (a)/(b)/(c)", "inputs": ";".join(inputs), "branch": branch, "filled_template": tmpl,
               "qualifiers": " | ".join(quals) + " || evidence: FinnGen rho_A %+.3f rho_B %+.3f p %s; ieu-a-7 rho_A %+.3f rho_B %+.3f p %s; discovery sign %+d; "
               "opposite signs in both %s; discovery sign in both %s" % (rA["finngen_chd"], rB["finngen_chd"], fmt_p(pf), rA["ieu_a_7"], rB["ieu_a_7"], fmt_p(pc), ds, opp, disc)})
    br.append({"rule_ref": "Rule 3(ii) on the named primary cell (TGFB1: FinnGen CHD conditional; ratio test aware 503, and the cell's conditional colocalizations)",
               "inputs": "S2.het.finngen_chd.p_aware_2nd;S2.het.finngen_chd.diff;S2.p489.het.finngen_chd.p_aware_2nd;S2.p489.het.finngen_chd.diff;S2.het.boot.finngen_chd.share;"
                         + ";".join("S2.coloc_clause.%s.finngen_chd.%s.PPH4" % (c["sig"], pn) for c in cell for pn in ("503", "489")) + ";"
                         + ";".join("S2.coloc_clause.boot.%s.finngen_chd.share" % c["sig"] for c in cell),
               "branch": "robust" if robust else "not robust",
               "filled_template": "no LD qualifier" if robust else "add the qualifier 'sensitive to the reference LD' to the s2 %s sentence" % branch,
               "qualifiers": "ratio test %s (503 p %s (d %+.3f); 489 p %s (d %+.3f); same conclusion on both panels %s; share %.3f); %s; threshold 0.90; "
                             "ieu-a-7 printed beside: ratio share %.4f, conditional colocalization share A %.4f / B %.4f; r(A,B) bootstrap 95%% interval %.3f to %.3f "
                             "(excludes 0: the sign of r is identified); Amendment 1 A1 front-surface choice is W9's" % (
                                 "robust" if ratio_robust else "NOT robust", fmt_p(pf), d503, fmt_p(p489), d489, same_panels, share,
                                 "; ".join("conditional colocalization at %s in the same primary cell %s (PP.H4 503 %.3f / 489 %.3f; share %.3f)" % (
                                     c["label"], "robust" if c["robust"] else "NOT robust", c["h503"], c["h489"], c["share"]) for c in cell),
                                 v("S2.het.boot.ieu_a_7.share"), cbs.get((SIGS["A"], "ieu_a_7"), float("nan")), cbs.get((SIGS["B"], "ieu_a_7"), float("nan")),
                                 v("S2.r.boot.lo"), v("S2.r.boot.hi"))})
    for sig, rs, f, zz, hh in fires:
        br.append({"rule_ref": "Rule 3(iv) colocalization clause (TGFB1 CAD; signal %s %s; %s)" % ("A" if sig == SIGS["A"] else "B", sig, "FinnGen CHD primary" if rs == "finngen_chd" else "ieu-a-7 printed beside"),
                   "inputs": ";".join("S2.coloc_clause.%s.%s.%s.%s" % (sig, rs, pn, k) for pn in ("503", "489") for k in ("z", "PPH4")),
                   "branch": "fires" if f else "does not fire",
                   "filled_template": ("the CAD association at signal %s did not colocalize with the protein (PP.H4 %.2f), so the ratio there is not anchored to the protein signal" % (sig, min(hh)))
                   if f else "no clause (|z| >= 1.96 and PP.H4 >= 0.80 on both panels)" if all(abs(z) >= 1.96 for z in zz) else "no clause (|z| < 1.96: worded by the contrast's own template)",
                   "qualifiers": "z 503 %.2f / 489 %.2f; PP.H4 503 %.4f / 489 %.4f%s" % (zz[0], zz[1], hh[0], hh[1],
                                 "; Rule 3(ii) bootstrap share (PP.H4 side of 0.80) %.4f -> %s" % (cbs[(sig, rs)], "robust" if cbs[(sig, rs)] >= 0.90 and ((hh[0] >= 0.8) == (hh[1] >= 0.8)) else
                                                                                                    "NOT robust: 'sensitive to the reference LD' attaches to the s2 %s sentence (named primary cell; see the Rule 3(ii) row) and to any statement of this colocalization" % branch
                                                                                                    if rs == "finngen_chd" else "NOT robust: 'sensitive to the reference LD' travels with the colocalization statement")
                                 if (sig, rs) in cbs else "")})
    # Rule 3(i)
    br.append({"rule_ref": "Rule 3(i) SuSiE with reference LD admissibility (TGFB1-CAD; s = max over susie/p1fr/p1np on 503 and p2np on 489) and its first-mention flag",
               "inputs": ";".join("S2.susie.%s.%s.s" % (tr, s) for tr in ("protein", "finngen_chd", "ieu_a_7") for s in ("susie", "p1fr", "p1np", "p2np", "p2np_503panel", "susie_489panel", "p1fr_489panel"))
                         + ";" + ";".join("S2.susie.%s.s_max" % tr for tr in ("protein", "finngen_chd", "ieu_a_7"))
                         + ";" + ";".join("S2.susie.%s.s_max_4fit" % tr for tr in ("protein", "finngen_chd", "ieu_a_7"))
                         + ";S2.rule3i.class.TGFB1_CAD.susie_class;S2.rule3i.pipeline_fit.TGFB1_FinnGen_CHD.PPH4_susie;S2.rule3i.pipeline_fit.ld_source",
               "branch": ("admissible" if adm else "inadmissible") + "; first-mention flag " + ("applies" if flag else "does not apply"),
               "filled_template": ("SuSiE PP.H4 and credible sets may be cited" if adm else
                                   "SuSiE (reference-LD) PP.H4 and credible sets reported in the supplement, labelled inadmissible, never cited as support in new text")
               + ("; TGFB1-CAD is flagged at first mention (Rule 3(i)), and the flag travels with any statement of its SuSiE class or tier (Principle 0): its "
                  "existing SuSiE class (%s) rests on reference-LD SuSiE that is inadmissible under Rule 3(i), labelled 'proxy s (locally computed; pipeline fit)'"
                  % t["susie_class"] if flag else ""),
               "qualifiers": "s_max over all seven distinct fits (Stage 2 decision) protein %.4f, CAD FinnGen CHD %.4f (ieu-a-7 beside %.4f); four-fit sensitivity "
                             "s_max protein %.4f, CAD FinnGen CHD %.4f (ieu-a-7 %.4f); threshold 0.10 for BOTH traits; locally computed s is a proxy for the "
                             "pipeline fit (pipeline n_snps %s vs local sets %d-%d); existing class: TAB_tiered81 (md5 %s) TGFB1/CAD susie_class %s, "
                             "crit_susie_supportive %s, evidence_tier %s; pipeline fit: TAB_susie_grid (md5 %s) TGFB1/FinnGen_R10_I9_CHD strict PP.H4_susie %.3f, "
                             "classification %s in %d of 4 settings; LD source: 00_config.R (md5 %s) l.%d-%d, 1KG EUR reference; the tier is not rebuilt and does not "
                             "appear to rest on the SuSiE criterion (%d pairs matching every other criterion without it are %s; descriptive, tier rule source not "
                             "local); Rule 3(i)'s aggregate sensitivity is declared as s6(d) (W6b); W2 computes none" % (
                                 sP, sF, sC, v("S2.susie.protein.s_max_4fit"), v("S2.susie.finngen_chd.s_max_4fit"), v("S2.susie.ieu_a_7.s_max_4fit"), gs["strict"]["n_snps"], min(int(x["nsnp"]) for x in st), max(int(x["nsnp"]) for x in st), m5["TAB_tiered81"], t["susie_class"], t["crit_susie_supportive"], t["evidence_tier"],
                                 m5["TAB_susie_grid"], fl(gs["strict"]["PP.H4_susie"]), t["susie_class"], sum(x["classification"] == t["susie_class"] for x in g),
                                 m5["ENG_engine_config_00_config.R"], ev["cfg_lines"][0], ev["cfg_lines"][1], len(ev["twins"]),
                                 "/".join(sorted({x["evidence_tier"] for x in ev["twins"]})) or "none")})
    # in-every-branch clause
    r2 = v("S2.ld.rs1800470_rs12983775.r2")
    br.append({"rule_ref": "s2 in every branch: rs12983775 epitope clause; no ancestry-specific reversal",
               "inputs": "S2.annot.rs1800470;S2.ld.rs1800470_rs12983775.r2", "branch": "applies (unconditional)",
               "filled_template": "any sentence about rs12983775 adds 'it tags a protein-altering allele, so an epitope or assay effect is not excluded'; TGFB1 is never described as an ancestry-specific reversal",
               "qualifiers": "rs1800470 = TGFB1 missense (P10L, ENST00000221930, Ensembl 116); r2 with rs12983775 = %.3f (503 PLINK)" % r2})
    W.write_tsv(W.OUT + "/branches.tsv", ["rule_ref", "inputs", "branch", "filled_template", "qualifiers"], br)
    W.log("results.tsv: %d rows; branches.tsv: %d rows; s2 branch %s; Rule 3(ii) %s (ratio test %s; conditional colocalization %s); Rule 3(i) %s, first-mention flag %s" % (
        len(rows), len(br), branch, "robust" if robust else "not robust", "robust" if ratio_robust else "not robust",
        ", ".join("%s %s (share %.4f)" % (c["label"], "robust" if c["robust"] else "not robust", c["share"]) for c in cell),
        "admissible" if adm else "inadmissible", "applies" if flag else "does not apply"))
    # completeness: every branch input exists and is non-empty
    miss = [i for b in br for i in b["inputs"].split(";") if i not in R or R[i]["value"] in ("", "NA")]
    if miss:
        W.log("ASSERT_FAILED branch inputs missing or NA: %s" % miss)
    na = [r["quantity_id"] for r in rows if r["value"] in ("", "NA") and ".x22." not in r["quantity_id"]]
    W.log("non-2x2 rows with empty/NA value: %d %s" % (len(na), na[:10]))


if __name__ == "__main__":
    main()
