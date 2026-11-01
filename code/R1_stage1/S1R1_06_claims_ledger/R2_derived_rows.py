# -*- coding: utf-8 -*-
# R2 derived ledger rows (2026-10-14, D1 round). Executed by ext_r2_stage1.py in the builder's namespace, after the protocol constants are
# loaded and BEFORE the Lxxx numbers are assigned. It fills the list _r2d_rows; ext_r2_stage1.py turns every entry into one ledger row
# (id -> Lxxx kept stable through R2_id_to_L_map.tsv, append-only) and runs the gate's own D1 check on it at build time.
# Every name starts with _r2d_ so nothing in the builder's namespace is overwritten.
#
# APPEND-ONLY: add new entries at the END of the list; never reorder, rename or delete (a previous id missing from here stops the build).
#
# Two kinds of entry:
#   SUCCESSOR  a row that repeats an existing row (the old Lxxx stays untouched) with the numbers its DISPLAY prints but its VALUE did not store
#              (the D1 gate found 13 such rows on 2026-10-14). The old row is read from ROWS, its stored leaves are asserted unchanged, and the
#              added leaves are read from other ledger rows or from the source file, never typed.
#   DERIVED    a count or value that is in no results.tsv and no governing document, computed from an existing local table with asserts.
# Nothing here runs a new analysis: every number is a count, a copy or a one-line transform of a value that an existing file already holds.
import json as _r2d_json, os as _r2d_os, re as _r2d_re

_r2d_rows = []


def _r2d_old(lid):
    hit = [r for r in ROWS if r["id"] == lid]
    assert len(hit) == 1, "ledger row %s not found exactly once" % lid
    return hit[0]


def _r2d_val(lid):
    return _r2d_json.loads(_r2d_old(lid)["value"])


def _r2d_const(cid):
    """Value of a protocol constant from R2_protocol_constants.tsv (loaded by ext_r2_stage1.py; its ledger rows do not exist yet at this point)."""
    hit = [r for r in _r2s_consts if r["id"] == cid]
    assert len(hit) == 1, cid
    return _r2d_json.loads(hit[0]["value"])


def _r2d_cell(pkg, qid):
    """Numeric value cell of a Stage 1 package results.tsv row (the same cell the ledger row of that placeholder id reads)."""
    hit = _r2s_BYQ_get(pkg, qid)
    assert len(hit) == 1, (pkg, qid)
    return float(hit[0][1]["value"])


def _r2d_add(**kw):
    need = {"id", "section", "claim", "value", "display", "source", "key", "note", "check"}
    assert set(kw) == need, set(kw) ^ need
    _r2d_rows.append(kw)


_r2d_TIER = _r2d_os.path.join(V3, "v3.2", "_engine_root", "results", "phase4_annotation", "table_s_pair_evidence_tiered.csv")
_r2d_EASAX = _r2d_os.path.join(V3, "V4", "source_data", "supp_tables", "SuppTable_EAS_axis.csv")
_r2d_OCF = _r2d_os.path.join(V3, "R1_stage1", "S1R1_05b_outcome_checks", "S1R1_outcome_checks.tsv")
_r2d_THEME = _r2d_os.path.join(V3, "paper3_deposit_request_v4", "scripts", "R", "theme_paper3.R")
_r2d_UMODQ = _r2d_os.path.join(V3, "R1_stage1", "S1R1_R2_1b_verify", "S6", "adjudication", "A_third", "a3_umod_q.tsv")
_r2d_t = pd.read_csv(_r2d_TIER)
assert len(_r2d_t) == 81 and _r2d_t.groupby(["gene", "outcome"]).ngroups == 81

# ---------------------------------------------------------------------------------------------------- SUCCESSOR of L256 (map significance split)
_r2d_v_o = _r2d_old("L256")
_r2d_v_ov = _r2d_val("L256")
_r2d_v_sc = _r2d_t.mr_sig_class.value_counts().to_dict()
_r2d_v_nom = _r2d_t[_r2d_t.mr_sig_class == "Suggestive"]
assert _r2d_v_ov == [_r2d_v_sc["Bonferroni"], _r2d_v_sc["FDR"], _r2d_v_sc["Suggestive"], "DPEP1-AF", _r2d_v_sc["Bonferroni"] + _r2d_v_sc["FDR"]], (_r2d_v_ov, _r2d_v_sc)
assert list(zip(_r2d_v_nom.gene, _r2d_v_nom.outcome)) == [("DPEP1", "AF")] and set(_r2d_v_sc) == {"Bonferroni", "FDR", "Suggestive"}
assert _r2d_v_o["display"] == "41 / 39 / 1 (DPEP1-AF); 80 of 81 at 5% FDR"
_r2d_v_fdr = _r2d_val("L239")
assert _r2d_v_fdr == 0.05
_r2d_add(id="DER.map_sig", section="Results: map (successor of L256)",
         claim=_r2d_v_o["claim"] + " [successor of L256: the same counts with the denominator 81 and the FDR level 0.05 stored, so the display derives from the value]",
         value=_r2d_v_ov + [len(_r2d_t), _r2d_v_fdr], display=_r2d_v_o["display"],
         source=rel(_r2d_TIER) + " (mr_sig_class); FDR level as L239; the old row L256 is unchanged",
         key="mr_sig_class counts; len(table); L239",
         note="D1 repair. L256's display printed '80 of 81 at 5% FDR' but its value stored only 80, so D1 reported 81 as a MISMATCH against 80 and 5 as UNDERIVED. Counts re-read here from the tiered table and asserted equal to the old leaves; 81 = rows of the table; 0.05 = L239.",
         check="DERIVED: counts recomputed from the tiered table and asserted equal to the stored leaves of L256 (41 Bonferroni, 39 FDR only, 1 Suggestive = DPEP1-AF, 80 = 41 + 39); one reader")

# ---------------------------------------------------------------------------------------------------- SUCCESSOR of L261 (resource p values)
_r2d_v_o = _r2d_old("L261")
_r2d_v_ov = _r2d_val("L261")
_r2d_v_oc = pd.read_csv(_r2d_OCF, sep="\t")
_r2d_v_ord6 = ["PAM-T2D", "RAD23B-T2D", "ARG1-T2D", "PCSK9-CAD", "TGFB1-CAD", "FGFR4-T2D (EAS)"]
_r2d_v_pf = [float(_r2d_v_oc.set_index("pair").loc[k, "p_out_file"]) for k in _r2d_v_ord6]
_r2d_v_a = _r2d_const("PD.alpha")
_r2d_v_m = _r2d_const("PD.ext.bonf_m")
assert _r2d_v_a == 0.05 and _r2d_v_m == 6
assert _r2d_v_ov == _r2d_v_pf + [int(sum(p < _r2d_v_a / _r2d_v_m for p in _r2d_v_pf))] and _r2d_v_ov[-1] == 4, (_r2d_v_ov, _r2d_v_pf)
_r2d_add(id="DER.ext_p", section="Results: external outcome checks (successor of L261)",
         claim=_r2d_v_o["claim"] + " [successor of L261: the same six p values and the count, with the level 0.05 and the denominator 6 stored, so the display derives from the value]",
         value=_r2d_v_ov + [_r2d_v_a, _r2d_v_m], display=_r2d_v_o["display"],
         source=OCS + " (p_out_file, first six rows); level as L928; denominator as L930; the old row L261 is unchanged",
         key="p_out_file (first six rows); L928; L930",
         note="D1 repair. L261's display ends '4 pass 0.05/6' but its value stored only the 4, so D1 reported 0.05 and 6 as UNDERIVED. The six p values are re-read from the outcome-check file and asserted equal to the old leaves; 0.05 = L928 and 6 = L930 (pre-declaration).",
         check="DERIVED: the six resource p values re-read from the outcome-check file and asserted equal to the stored leaves of L261; the count 4 re-derived as p < 0.05/6; one reader")

# ---------------------------------------------------------------------------------------------------- SUCCESSOR of L262 (external-check thresholds)
_r2d_v_o = _r2d_old("L262")
_r2d_v_ov = _r2d_val("L262")
assert _r2d_v_ov == [0.05, 0.05 / 6, 0.05] and _r2d_v_o["display"] == "classification p < 0.05; Bonferroni 0.05/6 = 0.0083; size difference p < 0.05"
_r2d_add(id="DER.ext_thr", section="Methods/Results: external outcome checks (successor of L262)",
         claim=_r2d_v_o["claim"] + " [successor of L262: the same three thresholds, with the denominator 6 stored, so the display derives from the value]",
         value=[_r2d_v_ov[0], _r2d_v_m, _r2d_v_ov[1], _r2d_v_ov[2]], display=_r2d_v_o["display"],
         source=_r2d_v_o["source"] + "; denominator as L930; the old row L262 is unchanged",
         key="pre-specified text; L930",
         note="D1 repair. L262's display prints 'Bonferroni 0.05/6' but its value stored only 0.05 / 6 = 0.0083, so D1 reported 6 as UNDERIVED. Leaves: classification 0.05, the denominator 6 (= L930), the cut 0.05/6, the post hoc size-difference heterogeneity 0.05. " + _r2d_v_o["note"],
         check="DERIVED: leaves of L262 asserted unchanged; the denominator is the value of L930 (pre-declaration line 182); one reader")

# ---------------------------------------------------------------------------------------------------- SUCCESSOR of L409 (grading rules)
_r2d_v_o = _r2d_old("L409")
_r2d_v_ov = _r2d_val("L409")
assert _r2d_v_ov["core_complete"] == {"core": 16, "convergent": 3, "supported_no_arm": 11, "flagged": 2} and _r2d_v_ov["prior_weak_below"] == 0.5 and _r2d_v_ov["susie_discordant_above"] == 0.5
with io.open(_r2d_THEME, encoding="utf-8") as _r2d_fh:
    _r2d_v_th = _r2d_fh.read()
_r2d_v_cc = _r2d_re.search(r"pairs\[, core_complete := \((.*?)\)\]", _r2d_v_th, flags=_r2d_re.S)
assert _r2d_v_cc is not None and _r2d_v_th.count("pairs[, core_complete :=") == 1 and "FOUR, not five." in _r2d_v_th
_r2d_v_crit = _r2d_re.findall(r"crit_\w+", _r2d_v_cc.group(1))
assert len(_r2d_v_crit) == len(set(_r2d_v_crit)) == 4, _r2d_v_crit
_r2d_v_cx = _r2d_t[_r2d_t.core_complete]
assert len(_r2d_v_cx) == 16 and int(((_r2d_v_cx.evidence_tier == "Convergent")).sum()) == 3 and int((_r2d_v_cx.evidence_tier == "Flagged").sum()) == 2
_r2d_v_nv = dict(_r2d_v_ov)
_r2d_v_nv["n_core_criteria"] = len(_r2d_v_crit)
_r2d_add(id="DER.grade", section="SM: grading (successor of L409)",
         claim=_r2d_v_o["claim"] + " [successor of L409: the same leaves plus the number of core criteria, so the display's 'all four internal criteria' derives from the value]",
         value=_r2d_v_nv, display=_r2d_v_o["display"],
         source=_r2d_v_o["source"] + "; core_complete definition read in " + rel(_r2d_THEME) + "; the old row L409 is unchanged",
         key="re-derivation asserted; n_core_criteria = crit_* terms in the core_complete conjunction",
         note="D1 repair. L409's display says 'all four internal criteria' but its value stored no 4 (the nearest leaf is the 3 Convergent pairs), so D1 reported 'four' as a MISMATCH. The count is read from the deposited tier code: core_complete is the conjunction of %s (the code comment reads 'FOUR, not five'). Composition 16 = 3 + 11 + 2 is asserted against the tiered table. " % ", ".join(_r2d_v_crit) + _r2d_v_o["note"],
         check="DERIVED: the four crit_* conjuncts of core_complete counted in the deposited tier code; composition 16 / 3 Convergent / 2 Flagged asserted against the tiered table and the old leaves; one reader")

# ---------------------------------------------------------------------------------------------------- SUCCESSOR of L210 (DPEP1 attribution count)
_r2d_v_o = _r2d_old("L210")
assert _r2d_val("L210") == 16 and _r2d_v_o["display"] == "16 of 48" and _r2d_val("L121") == 48
_r2d_add(id="DER.dpep1_attr", section="Results: kidney (successor of L210)",
         claim=_r2d_v_o["claim"] + " [successor of L210: the count with the number of settings (L121) stored]",
         value=[_r2d_val("L210"), _r2d_val("L121")], display=_r2d_v_o["display"],
         source=_r2d_v_o["source"] + "; settings as L121; the old row L210 is unchanged",
         key=_r2d_v_o["key"] + "; L121",
         note="D1 repair. L210's display prints '16 of 48' but its value stored only the 16, so D1 reported 48 as UNDERIVED. 48 is the value of L121 (analytic settings of the attribution grid).",
         check="DERIVED: leaves of L210 and L121 asserted; no new computation; one reader")

# ---------------------------------------------------------------------------------------------------- DERIVED: Rule 3(i) coverage inside the 81-pair map
_r2d_v_two = _r2d_t[((_r2d_t.gene == "TGFB1") & (_r2d_t.outcome == "CAD")) | ((_r2d_t.gene == "ARG1") & (_r2d_t.outcome == "T2D"))]
_r2d_v_ass = _r2d_const("S9.rule3i.n_map_pairs")
assert len(_r2d_v_two) == _r2d_v_ass == 2 and (_r2d_v_two.susie_class == "supportive").all()
_r2d_v_sup = int((_r2d_t.susie_class == "supportive").sum())
assert _r2d_v_sup == _r2d_val("L017")["supportive"] == 47
_r2d_v_rest = _r2d_t.drop(_r2d_v_two.index)
assert len(_r2d_v_rest) == 79 and int((_r2d_v_rest.susie_class == "supportive").sum()) == 45
assert int((_r2d_v_rest.evidence_tier == "Convergent").sum()) == 3 and int(_r2d_v_rest.core_complete.sum()) == 16
with io.open(_r2d_os.path.join(V3, "R2_rewrite", "review_R2", "CLAIM_OUTLINE_STAGE3_2026-10-11_v6.md"), encoding="utf-8") as _r2d_fh:
    _r2d_v_ol = _r2d_fh.read().splitlines()[15]   # line 16, the outline's H1 row
assert "The other 79 pairs" in _r2d_v_ol and "45 of the 47 SuSiE-supportive" in _r2d_v_ol and "exactly 2" in _r2d_v_ol, "the outline's H1 statement changed"
_r2d_add(id="DER.rule3i_cov", section="Results: map (Rule 3(i) coverage)",
         claim="Rule 3(i) (LD admissibility of reference-LD fine-mapping) was computed for 2 of the 81 map pairs, TGFB1-CAD and ARG1-T2D; the other 79 keep their original tiers and were not assessed; of the 47 pairs with a supportive multi-signal (SuSiE) class, 45 are among the 79",
         value={"n_pairs": len(_r2d_t), "assessed": int(_r2d_v_ass), "not_assessed": len(_r2d_v_rest), "supportive": _r2d_v_sup, "supportive_not_assessed": int((_r2d_v_rest.susie_class == "supportive").sum())},
         display="assessed for LD admissibility: 2 of 81 pairs; not assessed: 79, of which 45 of the 47 with a supportive multi-signal class",
         source=rel(_r2d_TIER) + " (susie_class); assessed pairs as L951",
         key="susie_class == supportive among the pairs other than TGFB1-CAD and ARG1-T2D; L951; L017",
         note="Counts only. The two assessed pairs are asserted to be in the table and both SuSiE-supportive; 47 = L017's supportive count; the 79 contain all 3 Convergent and all 16 core-complete pairs (asserted). Also stated in CLAIM_OUTLINE_STAGE3_2026-10-11_v6.md line 16 (H1: 'The other 79 pairs ... 45 of the 47 SuSiE-supportive'), which is the document that first gave these two counts; they are re-derived here from the table, not copied from it.",
         check="DERIVED: recomputed from the tiered table and compared with the outline's H1 statement (79 and 45 of 47); one reader")

# ---------------------------------------------------------------------------------------------------- DERIVED: UMOD all-signal Q without the palindromic protein signal
_r2d_v_uq = pd.read_csv(_r2d_UMODQ, sep="\t")
_r2d_v_ur = _r2d_v_uq[(_r2d_v_uq.reading.str.startswith("__alt: signal 5")) & (_r2d_v_uq.refine == "exact optimum")]
assert len(_r2d_v_ur) == 1 and int(_r2d_v_ur.df.iloc[0]) == 3
_r2d_v_ua = _r2d_v_uq[(_r2d_v_uq.reading.str.startswith("all 5 signals")) & (_r2d_v_uq.refine == "exact optimum")]
assert len(_r2d_v_ua) == 1 and round(float(_r2d_v_ua.p.iloc[0]), 2) == 0.45 and int(_r2d_v_ua.df.iloc[0]) == 4   # the declared Q of Results (p 0.45, df 4) is the same adjudicator file's primary reading
_r2d_add(id="DER.umod_q_drop5", section="Results: kidney (UMOD check outside the declared analyses)",
         claim="UMOD all-signal GLS Q with the palindromic protein signal 5 (rs62033284, A/T) dropped: Q, df, p [covariance-aware, second-order, 503 panel; the third adjudicator's own implementation, exact optimum]",
         value=[float(_r2d_v_ur.Q.iloc[0]), int(_r2d_v_ur.df.iloc[0]), float(_r2d_v_ur.p.iloc[0])], display="Q 3.182, df 3, p 0.364",
         source=rel(_r2d_UMODQ) + "; the reading '__alt: signal 5 (rs62033284, A/T) dropped', refine 'exact optimum'; same file's primary reading (5 signals): p 0.4473, df 4",
         key="reading=__alt signal 5 dropped; refine=exact optimum",
         note="Read from the verifier's output file, not from the analysis package's results table; no matching test was run for DPEP1. It is not a Stage 1b key and not a declared analysis. The coarse-refine reading of the same file gives p 0.36425 (both print 0.364). Listed as convention question CQ06 reading B in S1R1_R2_10_integration/convention_questions.tsv.",
         check="DERIVED: read from the adjudicator output a3_umod_q.tsv (and the same file's primary reading checked against the declared Q p 0.45); one reader, no second implementation of the dropped-signal Q beyond that file's own")

# ---------------------------------------------------------------------------------------------------- DERIVED: percent bounds of the arginine ratio of ratios
_r2d_v_lo, _r2d_v_hi = _r2d_cell("S1R1_R2_04a_arg1grid", "S4.arg.marg.ror_lo"), _r2d_cell("S1R1_R2_04a_arg1grid", "S4.arg.marg.ror_hi")
_r2d_v_pl, _r2d_v_ph = (_r2d_v_lo - 1) * 100, (_r2d_v_hi - 1) * 100
assert round(_r2d_v_pl) == -31 and round(_r2d_v_ph) == 54
_r2d_add(id="DER.ror_pct", section="Results: ARG1 (percent bounds of the arginine ratio of ratios)",
         claim="Arginine per SD of ARG1, marginal, covariance-aware, second-order SE, 503 panel: the differences between the signals that the 95% CI of the ratio of ratios does not exclude, as percent: lower bound (CI lower end - 1) x 100 and upper bound (CI upper end - 1) x 100",
         value=[round(_r2d_v_pl, 6), round(_r2d_v_ph, 6)], display="-31% to 54%",
         source="V3/R1_stage1/S1R1_R2_04a_arg1grid/results.tsv (S4.arg.marg.ror_lo and S4.arg.marg.ror_hi, as L802 and L801); formula of the pre-declaration (PRESPEC line 392: 'differences from {lo-1}% to {hi-1}% are not excluded')",
         key="(L802 - 1) x 100; (L801 - 1) x 100",
         note="One-line transform of L802 and L801; no new estimate.",
         check="DERIVED: arithmetic on L802 and L801, the pre-declared formula; one reader")

# ---------------------------------------------------------------------------------------------------- DERIVED: the East Asian eGFR file (EGFR_EAS_2P2M) pairs on the outcome side
_r2d_v_ea = pd.read_csv(_r2d_EASAX)
_r2d_v_ev = _r2d_v_ea[_r2d_v_ea.qc_evaluable == True]   # noqa: E712
assert len(_r2d_v_ev) == _r2d_val("L080")[0] == 52
_r2d_v_g = _r2d_v_ev[_r2d_v_ev.eas_source == "EGFR_EAS_2P2M"]
_r2d_v_st = _r2d_v_g.eas_outcome_state.value_counts().to_dict()
_r2d_v_ti = _r2d_v_g.evidence_tier.value_counts().to_dict()
assert len(_r2d_v_g) == 14 and _r2d_v_st == {"EAS_INCONCLUSIVE": 11, "EAS_COMPATIBLE": 3} and _r2d_v_ti == {"Supported": 9, "Flagged": 5}
assert bool(_r2d_v_g.eas_direction_same.all()) and (_r2d_v_g.outcome == "eGFR").all()
_r2d_v_stx = _r2d_v_ea[_r2d_v_ea.eas_outcome_state_strict.isin(["EAS_COMPATIBLE", "EAS_INCONCLUSIVE", "EAS_DISCORDANT"])]
assert len(_r2d_v_stx) == _r2d_val("L081")[0] == 27 and int((_r2d_v_stx.eas_source == "EGFR_EAS_2P2M").sum()) == 14
_r2d_add(id="DER.eas_egfr", section="Results: ancestry (East Asian eGFR file, outcome side)",
         claim="East Asian outcome-side pairs tested in the East Asian eGFR file EGFR_EAS_2P2M (identity and sample size unconfirmed): evaluable pairs, compatible, inconclusive, Supported, Flagged; all with the European sign; all in the strict subset",
         value={"evaluable": len(_r2d_v_g), "compatible": _r2d_v_st["EAS_COMPATIBLE"], "inconclusive": _r2d_v_st["EAS_INCONCLUSIVE"], "supported": _r2d_v_ti["Supported"], "flagged": _r2d_v_ti["Flagged"],
                "same_sign": int(_r2d_v_g.eas_direction_same.sum()), "in_strict_subset": int((_r2d_v_stx.eas_source == "EGFR_EAS_2P2M").sum())},
         display="14 evaluable pairs: 3 compatible, 11 inconclusive; 9 Supported, 5 Flagged; all 14 with the European sign and in the strict subset",
         source=rel(_r2d_EASAX) + " (eas_source, eas_outcome_state, evidence_tier, eas_direction_same, eas_outcome_state_strict)",
         key="eas_source == EGFR_EAS_2P2M among qc_evaluable pairs",
         note="Counts only, from the same table as L080 and L081 (52 evaluable pairs; 27 in the strict subset), which are asserted. The file's identity and N are unconfirmed (L329); no magnitude or scale statement is derived from these pairs.",
         check="DERIVED: counted from SuppTable_EAS_axis.csv with the totals of L080 and L081 asserted; one reader")

# ---------------------------------------------------------------------------------------------------- DERIVED: the (d) estimate under the four-fit / rows-dropped readings (SM Section S-R2)
_r2d_v_fs = pd.read_csv(_r2d_os.path.join(V3, "R1_stage1", "S1R1_R2_10_integration", "s6d_fitset_readings.tsv"), sep="	")
_r2d_v_dd = _r2d_v_fs[_r2d_v_fs.ING1_version.str.startswith("(d) estimate") & ~((_r2d_v_fs.fit_set == "seven") & (_r2d_v_fs.build_shift_rows == "keep"))]
_r2d_v_kk = _r2d_v_fs[_r2d_v_fs.ING1_version.str.startswith("'kept'") & ~((_r2d_v_fs.fit_set == "seven") & (_r2d_v_fs.build_shift_rows == "keep"))]
assert len(_r2d_v_dd) == len(_r2d_v_kk) == 3 and sorted(zip(_r2d_v_dd.fit_set, _r2d_v_dd.build_shift_rows)) == [("four", "drop"), ("four", "keep"), ("seven", "drop")]
for _r2d_c in ("n_reclassified", "or", "lo", "hi", "p"):
    assert _r2d_v_dd[_r2d_c].nunique() == 1 and _r2d_v_kk[_r2d_c].nunique() == 1, _r2d_c   # the three readings give identical values
_r2d_v_d1, _r2d_v_k1 = _r2d_v_dd.iloc[0], _r2d_v_kk.iloc[0]
assert round(float(_r2d_v_d1["or"]), 2) == 2.74 and int(_r2d_v_d1.n_reclassified) == 13 and round(float(_r2d_v_k1["or"]), 2) == 1.36 and float(_r2d_v_k1.p) == 1.0
_r2d_add(id="DER.s6d_alt", section="SM: renal (d) estimate under the four-fit and rows-dropped readings",
         claim="Kidney over-representation among the multi-signal candidates, estimate (d) (candidates on inadmissible reference-LD fits reclassified), under each of the readings that put NOTCH2-T2D below the s limit (four fits with the build-shifted rows kept or dropped; seven fits with them dropped): candidates reclassified, odds ratio, exact 95% CI, p; and the version that keeps ING1-CAD: odds ratio, exact 95% CI, p",
         value={"n_reclassified": int(_r2d_v_d1.n_reclassified), "or": float(_r2d_v_d1["or"]), "lo": float(_r2d_v_d1.lo), "hi": float(_r2d_v_d1.hi), "p": float(_r2d_v_d1.p),
                "kept_or": float(_r2d_v_k1["or"]), "kept_lo": float(_r2d_v_k1.lo), "kept_hi": float(_r2d_v_k1.hi), "kept_p": float(_r2d_v_k1.p)},
         display="13 reclassified: OR 2.74 (exact 95% CI 0.140 to 164.4), p 0.58; ING1-CAD kept: OR 1.36 (exact 95% CI 0.096 to 19.2), p 1",
         source="V3/R1_stage1/S1R1_R2_10_integration/s6d_fitset_readings.tsv (rows four/keep, four/drop, seven/drop, both ING1 versions)",
         key="fit_set, build_shift_rows, ING1_version; n_reclassified, or, lo, hi, p",
         note="The three readings give identical values (asserted). The seven-fit / rows-kept primary reading is a different row (OR not finite) and is not this row; its 'kept' version (OR 2.74, 0.140 to 164.4, p 0.58, 13 reclassified) coincides with this row's numbers but is a different reading (L870-L873). Read from the integration file the Stage 1b adjudicator's table feeds; the file's own columns w9_recount_agrees and adjudicator_levels_agree are TRUE on these rows.",
         check="DERIVED: read from s6d_fitset_readings.tsv with the three-way identity asserted; the file records an adjudicator recount (W9) and agreement of the adjudicator levels on these rows; one reader here")
assert bool(_r2d_v_dd.w9_recount_agrees.all()) and bool(_r2d_v_kk.w9_recount_agrees.all()) and bool(_r2d_v_dd.adjudicator_levels_agree.all())

# ---------------------------------------------------------------------------------------------------- DERIVED: CCND2-T2D Steiger direction test
_r2d_v_sx = pd.read_csv(_r2d_os.path.join(V3, "source_data", "figures", "Figure_04_replication_sensitivity", "replication_sensitivity__panel1_2.csv"))
_r2d_v_sc = _r2d_v_sx[(_r2d_v_sx.gene == "CCND2") & (_r2d_v_sx.outcome_id == "FinnGen_R10_T2D")]
assert len(_r2d_v_sc) == 1 and _r2d_v_sc.steiger_class.iloc[0] == "Non-forward" and bool(_r2d_v_sc.flag_reverse.iloc[0])
_r2d_v_s4 = pd.read_csv(_r2d_os.path.join(V3, "V4", "source_data", "figures", "Figure_04_replication_sensitivity", "replication_sensitivity__panel1_2.csv"))
_r2d_v_s4c = _r2d_v_s4[(_r2d_v_s4.gene == "CCND2") & (_r2d_v_s4.outcome_id == "FinnGen_R10_T2D")]
assert len(_r2d_v_s4c) == 1 and float(_r2d_v_s4c.steiger_z.iloc[0]) == float(_r2d_v_sc.steiger_z.iloc[0]) and float(_r2d_v_s4c.steiger_p.iloc[0]) == float(_r2d_v_sc.steiger_p.iloc[0])   # the V4 copy agrees
assert (_r2d_t[(_r2d_t.gene == "CCND2") & (_r2d_t.outcome == "T2D")].steiger_result == "BORDERLINE").all() and int((_r2d_t.steiger_result == "BORDERLINE").sum()) == 1   # the one non-forward pair of L019
_r2d_add(id="DER.ccnd2_steiger", section="SM: grading (Steiger direction of CCND2-T2D)",
         claim="Steiger directionality test for CCND2-T2D (the one non-forward pair of the 81): z and p [FinnGen R10 T2D, single-variant Wald ratio]",
         value=[float(_r2d_v_sc.steiger_z.iloc[0]), float(_r2d_v_sc.steiger_p.iloc[0])], display="z -3.407, p 7 × 10^-4",
         source="V3/source_data/figures/Figure_04_replication_sensitivity/replication_sensitivity__panel1_2.csv (row CCND2 / FinnGen_R10_T2D: steiger_z, steiger_p, steiger_class Non-forward); the V4 copy of the file agrees; the tiered table marks the pair BORDERLINE (L019)",
         key="gene == CCND2, outcome_id == FinnGen_R10_T2D; steiger_z, steiger_p",
         note="Printed to 4 decimals in the figure source file (z -3.407, p 0.0007); the outline lists 'CCND2 Steiger z -3.41 and p 7e-4' as found in the plan only, and this row finds them in a primary figure source file.",
         check="DERIVED: read from the figure source table (V3 and V4 copies identical on this row) and checked against the tiered table's single BORDERLINE pair; one reader")

# ==================================================================================================== S2-D25 (2026-10-16): the four rows L985-L988
# Carried here from figures/registry/ledger_append_S2D24_B.py, which had appended them to S1R1_CLAIMS_LEDGER.csv by hand. Their Lxxx are
# pre-seeded in R2_id_to_L_map.tsv (DER.corroborated L985, DER.credset_sizes L986, DER.eas_partition L987, DER.r_ab_interval L988), so a rebuild
# gives the same ids. Values are READ from the source files (never typed) and asserted against the neighbouring rows, as in the append script.
# Two display repairs found by the D1 check (P31R_gate.py --d1-all): L987's display printed 'inconclusive 40' but its value stored only the
# 33 + 7 split; L986's display printed the method numbers (3), (4), (5) and the panel sizes 503 / 489 as labels that its value did not hold.
# Both now store what they print (added leaves only; every earlier leaf and key is unchanged, so every idx used by a figure still resolves).
import csv as _r2d_csv, collections as _r2d_col, io as _r2d_io
_r2d_csv.field_size_limit(10 ** 9)


def _r2d_tsv(p, d="\t"):
    with _r2d_io.open(p, newline="", encoding="utf-8") as _r2d_fhx:
        return list(_r2d_csv.DictReader(_r2d_fhx, delimiter=d))


# ---------------- L985: pairs externally corroborated under the rule fixed in advance
_r2d_CO = _r2d_os.path.join(V3, "R1_stage1", "S1R1_R2_10_integration", "corroboration.tsv")
_r2d_w_co = _r2d_tsv(_r2d_CO)
_r2d_w_cor = [r["pair"] for r in _r2d_w_co if r["verdict"] == "CORROBORATED"]
assert _r2d_w_cor == ["PAM-T2D", "ARG1-T2D", "FGFR4-T2D (EAS)"], _r2d_w_cor
assert [r["verdict"] for r in _r2d_w_co if r["pair"] not in _r2d_w_cor and r["verdict"] != "NOT_A_PAIR"] == ["NOT_CORROBORATED"] * 3
_r2d_w_names = [c.replace(" (EAS)", "") for c in _r2d_w_cor]
_r2d_add(id="DER.corroborated", section="Results: external outcome (corroboration rule, Figure 1 and Figure 2 message)",
         claim="Pairs externally corroborated under the corroboration rule fixed in advance: count and pairs",
         value={"n": len(_r2d_w_cor), "pairs": _r2d_w_names},
         display="%d pairs: %s" % (len(_r2d_w_cor), ", ".join(_r2d_w_names)),
         source="V3/R1_stage1/S1R1_R2_10_integration/corroboration.tsv (md5 %s); column verdict" % _r2s_md5(_r2d_CO),
         key="verdict == CORROBORATED (six pre-declared checks; the ARG1 signal-B row is NOT_A_PAIR and is excluded)",
         note="Names agree with the Results sentence ('only PAM-T2D, ARG1-T2D and FGFR4-T2D (in East Asians) are externally corroborated'). "
              "The other three checks (RAD23B-T2D, PCSK9-CAD, TGFB1-CAD) are NOT_CORROBORATED. The verdicts are computed by code from Stage 1b KEY ids "
              "(STAGE1B_REPORT D-4), so this count is PIPE, not a key.",
         check="DERIVED: counted from the verdict column of corroboration.tsv; the six rows and their verdicts were read, not typed; one reader")

# ---------------- L986: credible-set size by method at ARG1 signals A and B
_r2d_CS = _r2d_os.path.join(V3, "R1_stage1", "S1R1_R2_12_chen_instrument", "credible_set_status.tsv")
_r2d_w_tmp = _r2d_col.defaultdict(set)
for _r2d_rx in _r2d_tsv(_r2d_CS):
    _r2d_panx = {"503": "503", "489": "489"}.get(_r2d_rx["panel"], "insample")
    _r2d_mx = _r2d_rx["method"]
    _r2d_kx = ("m1_T2D_503" if _r2d_mx == "1_T2D" else "m1_protein_503" if _r2d_mx == "1_protein" else "m2_insample" if _r2d_mx == "2" else "m%s_%s" % (_r2d_mx, _r2d_panx))
    _r2d_w_tmp[(_r2d_rx["set_signal"], _r2d_kx)].add(int(_r2d_rx["set_size"]))
assert all(len(_r2d_vx) == 1 for _r2d_vx in _r2d_w_tmp.values()), "a method has two sizes"
_r2d_w_sz = _r2d_col.OrderedDict()
for _r2d_sgx in ("A", "B"):
    _r2d_w_sz[_r2d_sgx] = _r2d_col.OrderedDict((_r2d_kx[1], list(_r2d_vx)[0]) for _r2d_kx, _r2d_vx in sorted(_r2d_w_tmp.items(), key=lambda kv: kv[0][1]) if _r2d_kx[0] == _r2d_sgx)
assert _r2d_w_sz["A"]["m2_insample"] == 2 and _r2d_w_sz["A"]["m1_protein_503"] == 2
_r2d_ax = _r2d_w_sz["A"]
_r2d_w_lab = _r2d_col.OrderedDict(methods=[3, 4, 5], panels=[503, 489])   # D1 repair: the numerals the display prints as labels
_r2d_add(id="DER.credset_sizes", section="Results: ARG1 (Figure 4b credible-set tick rows)",
         claim="Credible-set size by method at ARG1 signal A and signal B (variants in the set): SuSiE T2D, SuSiE protein (503 reference), UKB-PPP in-sample set, conditional ABF sets (3)-(5) on the 503 and 489 panels",
         value=_r2d_col.OrderedDict(list(_r2d_w_sz.items()) + [("display_labels", _r2d_w_lab)]),
         display="Signal A: SuSiE T2D %d, SuSiE protein %d, in-sample %d; conditional (3) %d/%d, (4) %d/%d, (5) %d/%d (503/489)" % (
             _r2d_ax["m1_T2D_503"], _r2d_ax["m1_protein_503"], _r2d_ax["m2_insample"], _r2d_ax["m3_503"], _r2d_ax["m3_489"], _r2d_ax["m4_503"], _r2d_ax["m4_489"], _r2d_ax["m5_503"], _r2d_ax["m5_489"]),
         source="V3/R1_stage1/S1R1_R2_12_chen_instrument/credible_set_status.tsv (md5 %s); column set_size, one value per (set_signal, method, panel)" % _r2s_md5(_r2d_CS),
         key="set_signal in {A,B}; method 1_T2D, 1_protein, 2, 3, 4, 5; panel 503, 489 or UKB-PPP in-sample LD; idx = signal then key (e.g. [\"A\",\"m2_insample\"])",
         note="Sizes only; admissibility and robustness are read from the same file's admissibility_rule3i and robustness_rule3ii columns at render, "
              "never inferred for the whole locus. Methods (1) and (2) carry one panel each; (3)-(5) carry both. A size is not a membership: the file lists members per variant. "
              "D1 repair (S2-D25): display_labels stores the method numbers (3, 4, 5) and the panel sizes (503, 489) that the display prints as labels; the sizes are unchanged.",
         check="DERIVED: read from credible_set_status.tsv with the one-size-per-cell assertion; the in-sample set of 2 agrees with the ST16 list in the same file's source_extra (rs2781668, rs2246012); one reader")

# ---------------- L987: East Asian outcome side, sign x significance partitions (canonical tiered table, copy B)
_r2d_TIERB = _r2d_os.path.join(V3, "v3.2", "_engine_root", "results", "repro_phase1_hg38", "_fig_results", "phase4_annotation", "table_s_pair_evidence_tiered.csv")
_r2d_w_ti = _r2d_tsv(_r2d_TIERB, ",")
assert len(_r2d_w_ti) == 81
_r2d_w_ev = [r for r in _r2d_w_ti if r["eas_outcome_state"] != "EAS_NOT_EVALUABLE"]
_r2d_w_ct = _r2d_col.Counter((r["eas_outcome_state"], r["eas_direction_same"]) for r in _r2d_w_ev)
_r2d_w_p = _r2d_col.OrderedDict(
    evaluable=len(_r2d_w_ev), same_sign=sum(1 for r in _r2d_w_ev if r["eas_direction_same"] == "TRUE"), opposite_sign=sum(1 for r in _r2d_w_ev if r["eas_direction_same"] == "FALSE"),
    compatible=_r2d_col.OrderedDict(same=_r2d_w_ct[("EAS_COMPATIBLE", "TRUE")], opposite=_r2d_w_ct[("EAS_COMPATIBLE", "FALSE")]),
    inconclusive=_r2d_col.OrderedDict(same=_r2d_w_ct[("EAS_INCONCLUSIVE", "TRUE")], opposite=_r2d_w_ct[("EAS_INCONCLUSIVE", "FALSE")]),
    discordant=_r2d_col.OrderedDict(same=_r2d_w_ct[("EAS_DISCORDANT", "TRUE")], opposite=_r2d_w_ct[("EAS_DISCORDANT", "FALSE")]),
    not_evaluable=sum(1 for r in _r2d_w_ti if r["eas_outcome_state"] == "EAS_NOT_EVALUABLE"),
    egfr_evaluable=sum(1 for r in _r2d_w_ev if r["outcome"] == "eGFR"), binary_evaluable=sum(1 for r in _r2d_w_ev if r["outcome"] != "eGFR"))
for _r2d_kx in ("compatible", "inconclusive", "discordant"):   # D1 repair: the totals the display prints (compatible 12, inconclusive 40, discordant 0)
    _r2d_w_p[_r2d_kx + "_total"] = _r2d_w_p[_r2d_kx]["same"] + _r2d_w_p[_r2d_kx]["opposite"]
_r2d_w_eg = [r for r in _r2d_rows if r["id"] == "DER.eas_egfr"]
assert len(_r2d_w_eg) == 1
assert [_r2d_w_p["evaluable"], _r2d_w_p["compatible_total"], _r2d_w_p["inconclusive_total"], _r2d_w_p["discordant_total"], _r2d_w_p["not_evaluable"]] == _r2d_val("L080"), "L080 disagrees"
assert [_r2d_w_p["same_sign"], _r2d_w_p["evaluable"]] == _r2d_val("L224")[:2], "L224 disagrees"
assert _r2d_w_p["egfr_evaluable"] == _r2d_w_eg[0]["value"]["evaluable"], "DER.eas_egfr (L967) disagrees"
assert _r2d_w_p["same_sign"] + _r2d_w_p["opposite_sign"] == _r2d_w_p["evaluable"] and _r2d_w_p["egfr_evaluable"] + _r2d_w_p["binary_evaluable"] == _r2d_w_p["evaluable"]
_r2d_add(id="DER.eas_partition", section="Results: ancestry (Figure 5a and AF3 Figure S6: explicit partitions)",
         claim="East Asian outcome-side evaluable pairs: sign versus European (same / opposite) and significance state (compatible / inconclusive / discordant) cross-tabulated; not evaluable; eGFR-file and binary-outcome pairs",
         value=_r2d_w_p,
         display="%d evaluable: %d same sign, %d opposite; compatible %d (all same sign), inconclusive %d (%d same, %d opposite), discordant %d; %d not evaluable; %d eGFR-file pairs, %d binary-outcome pairs" % (
             _r2d_w_p["evaluable"], _r2d_w_p["same_sign"], _r2d_w_p["opposite_sign"], _r2d_w_p["compatible_total"], _r2d_w_p["inconclusive_total"],
             _r2d_w_p["inconclusive"]["same"], _r2d_w_p["inconclusive"]["opposite"], _r2d_w_p["discordant_total"], _r2d_w_p["not_evaluable"],
             _r2d_w_p["egfr_evaluable"], _r2d_w_p["binary_evaluable"]),
         source="V3/v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv (md5 %s); columns eas_outcome_state, eas_direction_same, outcome; memberships in V3/R2_rewrite/figures/data/Figure_5/eas_partition_members.tsv" % _r2s_md5(_r2d_TIERB),
         key="eas_outcome_state != EAS_NOT_EVALUABLE; eas_direction_same TRUE/FALSE; outcome == eGFR vs not",
         note="Counts only. The sign partition (45 / 7) and the significance partition (12 / 40 / 0) overlap and must not be drawn as nested: the 40 inconclusive pairs are 33 same-sign and 7 opposite-sign, so they are not a subset of the 45. "
              "No magnitude is derived: the 14 eGFR pairs are tested in EGFR_EAS_2P2M (HOLD H9, L967, L329). "
              "D1 repair (S2-D25): compatible_total, inconclusive_total and discordant_total store the totals the display prints (12, 40, 0); every earlier key is unchanged.",
         check="DERIVED: tallied from the tiered table with the totals of L080 (52 / 12 / 40 / 0 / 29), L224 (45 of 52) and L967 (14 eGFR pairs) asserted; one reader")

# ---------------- L988: resampling interval of r(A,B) at ARG1
_r2d_W4A = _r2d_os.path.join(V3, "R1_stage1", "S1R1_R2_04a_arg1grid", "results.tsv")
_r2d_w_w = {r["quantity_id"]: r for r in _r2d_tsv(_r2d_W4A)}
_r2d_w_lo, _r2d_w_hi, _r2d_w_r0 = float(_r2d_w_w["S4.ld.boot.r_lo"]["value"]), float(_r2d_w_w["S4.ld.boot.r_hi"]["value"]), float(_r2d_w_w["S4.ld.r.503"]["value"])
assert abs(_r2d_w_r0 - _r2d_val("L058")) < 5e-5 and _r2d_w_lo < _r2d_w_r0 < _r2d_w_hi and _r2d_w_lo < 0 < _r2d_w_hi
_r2d_add(id="DER.r_ab_interval", section="Results: ARG1 (Figure 4b, sign uncertainty of r between the two signals)",
         claim="95% resampling interval for r(A,B), the correlation between ARG1 signals A (rs2781668) and B (rs17788484): lower and upper (503 people resampled 5,000 times; point estimate in L058)",
         value=[_r2d_w_lo, _r2d_w_hi], display="%.2f to %.2f" % (_r2d_w_lo, _r2d_w_hi),
         source="V3/R1_stage1/S1R1_R2_04a_arg1grid/results.tsv (md5 %s); quantity_id=S4.ld.boot.r_lo and S4.ld.boot.r_hi; column value" % _r2s_md5(_r2d_W4A),
         key="quantity_id in {S4.ld.boot.r_lo, S4.ld.boot.r_hi}",
         note="The interval includes 0, so r = -0.04 (L058) is not a securely negative correlation; the sign is uncertain. W4a records this bootstrap as a DEVIATION (bootstrap outside PLINK, section 10).",
         check="DERIVED: read from the W4a results table; the point estimate agrees with L058 (asserted); one reader, bootstrap not independently re-run")
