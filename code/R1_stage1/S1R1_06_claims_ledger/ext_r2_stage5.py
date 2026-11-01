# -*- coding: utf-8 -*-
"""Stage 5 step 1c (written 2026-10-18): ledger rows L500-L520 = the successor rows and the placeholder rows of Stage 5.

Executed by S1R1_build_ledger.py in ITS namespace (add, rel, ROWS, pd, np, V3, CKM, HERE ... are available). Every name defined here starts with s5_.
Specs: STAGE5_SUCCESSOR_SPECS.tsv (the rows with status 'NEW ROW REQUIRED', in file order) then STAGE5_NEW_ROW_SPECS.tsv (S5N01-S5N13).
NOTHING IS TYPED: every value is read from the source named in the row (a file, a line, a column, a quantity id) or from the ledger row that the spec
names, structural anchors are asserted, and each display is built from the value read and then asserted equal to the display the spec prints
(the specs print ASCII 'x' for the multiplication sign; the ledger prints U+00D7, so the comparison folds the two).
The rows L500-L507 are SUCCESSORS: the predecessor row stays (its generator still reads it); the successor stores the leaves the display needs.
Ids are allocated consecutively from L500 and checked against the 'expected' id in the spec.
"""
import csv as s5_csv, io as s5_io, json as s5_json, os as s5_os, re as s5_re, hashlib as s5_hl, datetime as s5_dt

s5_csv.field_size_limit(10 ** 9)
s5_LED = {r["id"]: r for r in ROWS}                       # rows built so far (base rows and earlier extension modules)
s5_val = lambda i: s5_json.loads(s5_LED[i]["value"])
s5_P = lambda *a: s5_os.path.join(V3, *a)
s5_txt = lambda p: s5_io.open(p, encoding="utf-8").read()


def s5_tsv(p):
    with s5_io.open(p, encoding="utf-8", newline="") as fh:
        return list(s5_csv.DictReader(fh, delimiter="\t"))


def s5_md5(p):
    return s5_hl.md5(open(p, "rb").read()).hexdigest()


def s5_norm(t):
    return t.replace("\u00d7", "x")


s5_SPEC = {r["predecessor"]: r for r in s5_tsv(s5_os.path.join(HERE, "STAGE5_SUCCESSOR_SPECS.tsv")) if r["status"] == "NEW ROW REQUIRED"}
s5_ORDER = [r["predecessor"] for r in s5_tsv(s5_os.path.join(HERE, "STAGE5_SUCCESSOR_SPECS.tsv")) if r["status"] == "NEW ROW REQUIRED"]
assert s5_ORDER == ["L454", "L242", "L094", "L403", "L415", "L438", "L439", "L446"], s5_ORDER
s5_NSPEC = [r for r in s5_tsv(s5_os.path.join(HERE, "STAGE5_NEW_ROW_SPECS.tsv"))]
assert [r["placeholder"] for r in s5_NSPEC] == ["S5N%02d" % i for i in range(1, 14)]
s5_next = [500]
s5_idmap = []


def s5_alloc(spec_key, expected=None):
    cid = "L%d" % s5_next[0]
    s5_next[0] += 1
    if expected:
        assert cid == expected, "allocated %s, the spec expects %s (%s)" % (cid, expected, spec_key)
    s5_idmap.append((spec_key, cid))
    return cid


def s5_succ(pred, value, display, section_tail, claim, source, key, note, tag=None):
    """Add the successor of `pred`. The display is checked against the spec's printed display (and, for an unchanged text, the predecessor's)."""
    sp = s5_SPEC[pred]
    exp = sp["successor_id"]
    m = s5_re.search(r"expected (L\d+)", exp)
    sd = sp["successor_display"]
    want = sp["predecessor_display"] if ("identical text" in sd or "unchanged text" in sd) else sd
    want = s5_re.sub(r"\s*\((?:unchanged text|identical text)\)\s*$", "", want)
    assert s5_norm(display) == s5_norm(want), "display of the successor of %s differs from the spec's printed display:\n  built: %r\n  spec : %r" % (pred, display, want)
    if "unchanged text" in sd or "identical text" in sd:
        assert display == s5_LED[pred]["display"], "the successor of %s must print the predecessor's text unchanged" % pred
    cid = s5_alloc(pred, m.group(1))
    add(cid, s5_LED[pred]["section"] + " (successor of %s)" % pred,
        claim + " [successor of %s: %s]" % (pred, note[0]), value, display, source + "; the old row %s is unchanged" % pred, key,
        tag or s5_LED[pred]["tag"], s5_LED[pred]["load_bearing"], note=note[1])
    return cid


# ======================================================================================================== L500: successor of L454 (the audit date as an ISO string)
s5_au = pd.read_csv(s5_P("R2_rewrite", "data", "P31R_devstage_audit.tsv"), sep="\t", dtype=str)
s5_dates = set()
for s5_s in s5_au.audit_source:
    s5_dates |= set(s5_re.findall(r"(20\d\d)-(\d\d)-(\d\d)", s5_s))
assert len(s5_au) == 13 and s5_dates and max(s5_dates) == ("2026", "10", "06") and s5_val("L454") == [6, 10, 2026], (s5_dates, s5_val("L454"))
s5_iso = "%s-%s-%s" % max(s5_dates)
s5_d = s5_dt.date(*[int(x) for x in max(s5_dates)])
s5_disp454 = "%d %s %d" % (s5_d.day, ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"][s5_d.month - 1], s5_d.year)
assert [s5_d.day, s5_d.month, s5_d.year] == s5_val("L454")
s5_succ("L454", s5_iso, s5_disp454, "", "Date of the Stage 3b development-stage audit in Open Targets (every protein with a direct approved or clinical-stage agent)",
        "V3/R2_rewrite/data/P31R_devstage_audit.tsv (audit_source, latest query date; the regex of ext_af1.py l.72-79)", "latest query date (ISO string)",
        ("the same date stored as an ISO string, so that the display date derives from the value",
         "Stage 5 step 1c. D1 repair: L454 stores [6, 10, 2026], and D1 reads a date display only against an ISO date that occurs verbatim in the stored value. The date is the maximum of the query dates in the 13 audit_source cells, "
         "asserted to be 2026-10-06 and to equal the three leaves of L454."))

# ======================================================================================================== L501: successor of L242 (Bonferroni cut at three significant figures)
s5_cut = s5_val("L239") / s5_val("L003")
assert s5_val("L239") == 0.05 and s5_val("L003") == 15584 and s5_cut == s5_val("L242") == 3.2084188911704316e-06
s5_e = int(np.floor(np.log10(s5_cut)))
s5_disp242 = "%.2f \u00d7 10^%d" % (s5_cut / 10.0 ** s5_e, s5_e)
s5_succ("L242", s5_cut, s5_disp242, "", "Bonferroni threshold over the 15,584-test family, printed to three significant figures",
        "0.05 / L003 (the level as L239, the family size as L003); generator S1R1_build_ledger.py", "derived: L239 / L003",
        ("the same value, printed to three significant figures as the R2 Methods, Results, Table 2 and the SM files print it",
         "Stage 5 step 1c. L242 prints 3.2 x 10^-6 (two figures, passes D1); the R2 files print 3.21 x 10^-6. Derived from the ledger values of L239 and L003, not typed; asserted equal to the stored L242 value. "
         "SPEC DISAGREEMENT: the spec prints the display with an ASCII 'x'; the ledger convention (L242, L960) is U+00D7, and the R2 files print U+00D7."))

# ======================================================================================================== L502: successor of L094 (decoy truncation rates computed from the counts)
s5_log = s5_txt(s5_P("R1_stage1", "S1R1_03_decoy_method", "S1R1_decoy_rates.log")).splitlines()
s5_l11 = [l for l in s5_log if l.startswith("CKD draws that truncate the target vector:")]
s5_l10 = [l for l in s5_log if l.startswith("CKD annotated draws:")]
s5_l13 = [l for l in s5_log if l.startswith("CKD decoy passes with the located instrument row cut away:")]
assert len(s5_l11) == len(s5_l10) == len(s5_l13) == 1
s5_m = s5_re.search(r"truncate the target vector: (\d+) \(([\d.]+)%\); that cut a LOCATED instrument row \(idx > k\): (\d+) \(([\d.]+)%\)", s5_l11[0])
s5_trunc, s5_cutrow = int(s5_m.group(1)), int(s5_m.group(3))
s5_draws = int(s5_re.search(r"CKD annotated draws: (\d+);", s5_l10[0]).group(1))
s5_m = s5_re.search(r"cut away: (\d+) of (\d+) passes", s5_l13[0])
s5_pcut, s5_npass = int(s5_m.group(1)), int(s5_m.group(2))
s5_l094, s5_l092 = s5_val("L094"), s5_val("L092")
assert [s5_trunc, s5_cutrow, s5_pcut] == s5_l094 and s5_draws == s5_l092[1] and s5_npass == s5_l092[0]
s5_v094 = [s5_trunc / s5_draws, s5_cutrow / s5_draws, s5_pcut, s5_npass]
s5_pct = s5_re.search(r"truncate the target vector: \d+ \(([\d.]+)%\); that cut a LOCATED instrument row \(idx > k\): \d+ \(([\d.]+)%\)", s5_l11[0]).groups()
assert ("%.1f" % (100 * s5_v094[0]), "%.1f" % (100 * s5_v094[1])) == s5_pct, (s5_v094, s5_pct)
s5_succ("L094", s5_v094, "%.1f%%; %.1f%%; %d of %d" % (100 * s5_v094[0], 100 * s5_v094[1], s5_v094[2], s5_v094[3]), "",
        "CKD decoy draws that truncate the target vector / cut a located instrument row (as shares of all CKD annotated draws); passes with the instrument row cut away, of all CKD passes",
        "V3/R1_stage1/S1R1_03_decoy_method/S1R1_decoy_rates.log (lines 10, 11, 13); counts as L094, draws and passes as L092",
        "truncating draws / L092 draws; located-row draws / L092 draws; passes with the row cut; L092 passes",
        ("the two shares and the pass counts stored as the numbers the display prints",
         "Stage 5 step 1c. D1 repair: L094 stores the raw counts 8160, 1891, 16, whereas the display prints their shares of the 16,200 CKD draws and '16 of 50'. The shares are computed here from the log counts "
         "(8160 / 16200, 1891 / 16200), asserted to round to the log's own printed 50.4% and 11.7%, and to agree with the leaves of L094 and L092; 16 = L094 leaf [2] (log l.13), 50 = L092 leaf [0] (log l.13). "
         "The raw counts are not stored because no display quantity reads them."))

# ======================================================================================================== L503: successor of L403 (SuSiE funnel as a list)
s5_f = s5_val("L403")
s5_ca = pd.read_csv(s5_P("V4", "source_data", "supp_tables", "SuppTable_candidate_edge_audit.csv"))
s5_stc = s5_ca.Status.value_counts().to_dict()
s5_drop = s5_ca[s5_ca.Status == "Dropped"]
s5_held = s5_ca[s5_ca.Status == "Held (not drawn)"]
assert len(s5_ca) == s5_f["first_fit_ge080"] and s5_stc["Drawn"] == s5_f["kept"] and s5_stc["Held (not drawn)"] == s5_f["not_converged"] and s5_stc["Dropped"] == s5_f["below_080"] + s5_f["no_pure_cs"]
assert int(s5_drop.Reason.str.contains("< 0.80 under purity filter").sum()) == s5_f["below_080"] and int(s5_drop.Reason.str.contains("no pure credible set").sum()) == s5_f["no_pure_cs"]
assert int(s5_re.search(r"max_iter=(\d+)", s5_held.Reason.iloc[0]).group(1)) == s5_f["max_iter"] and len(s5_held) == 1
s5_l263, s5_l236 = s5_val("L263"), s5_val("L236")
assert s5_l263[3] == s5_f["testable_egfr"] + s5_f["testable_nonegfr"] == 160 and s5_f["first_fit_ge080"] + s5_f["not_supportive"] == s5_l263[3]
s5_v403 = [s5_f["first_fit_ge080"], s5_l263[3], s5_l236, s5_f["nonegfr_run"], s5_f["egfr_run"], s5_f["not_supportive"], s5_f["kept"], s5_f["not_converged"],
           s5_f["below_080"], s5_l236, s5_f["no_pure_cs"], s5_f["max_iter"]]
s5_hg, s5_he = s5_held.Gene.iloc[0], s5_held.Endpoint.iloc[0]
s5_succ("L403", s5_v403,
        "%d of %d reached PP.H4 >= %.2f on the first fit (%d non-eGFR, %d eGFR); %d did not; kept %d, %d not converged (%s-%s), %d below %.2f, %d no pure credible set; IBSS cap %s iterations"
        % (tuple(s5_v403[:8]) + (s5_hg, s5_he, s5_v403[8], s5_v403[9], s5_v403[10], "{:,}".format(s5_v403[11]))), "",
        "SuSiE candidate funnel over the 160 testable associations: first fit PP.H4 >= 0.80 (of the pool; non-eGFR / eGFR run) / not; purity re-run kept / not converged / below 0.80 / no pure credible set; IBSS iteration cap",
        "V3/V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv + run logs 29_cand_*.log, 30_candidate_audit.log (as L403); pool size as L263 leaf [3]; PP.H4 level as L236",
        "L403 leaves in display reading order; L263[3]; L236 (twice)",
        ("the same funnel stored as a list in display reading order, with the pool size 160 and the 0.80 level (twice) added",
         "Stage 5 step 1c. A list, not an object: D1 reads an object with the any-leaf rule, under which the display token '1' (1 not converged) is matched by 0.8. Leaves from L403 (keys first_fit_ge080, nonegfr_run, egfr_run, "
         "not_supportive, kept, not_converged, below_080, no_pure_cs, max_iter), re-counted here from the audit CSV (Status, Reason, max_iter); 160 = L263 leaf [3] = testable_egfr + testable_nonegfr; 0.8 = L236, used for the "
         "first-fit threshold and the purity threshold."))

# ======================================================================================================== L504: successor of L415 (decoy provenance as a list)
s5_PPH4_DP = 3                                           # decimal places at which the null file stores PP.H4 (the .round(3) of ext_supp.py's xs_eq comparison)
s5_supp = s5_txt(s5_os.path.join(HERE, "ext_supp.py"))
assert 'xs_o["PP.H4"].round(%d)' % s5_PPH4_DP in s5_supp, "ext_supp.py no longer compares at %d dp" % s5_PPH4_DP
s5_rep = s5_txt(s5_P("paper3_deposit_request_v4", "REPORT.md"))
s5_m = s5_re.search(r"\| A17 \| `results/repro_phase1_hg38/coloc_null/null_draws\.csv` \| \d+ \| \d+ \| (\d{4}-\d\d-\d\d) \d\d:\d\d:\d\d", s5_rep)
s5_run = s5_m.group(1)
s5_m = s5_re.search(r"predate the (\d{4}-\d\d-\d\d) engine fix", s5_rep)
s5_fix = s5_m.group(1)
s5_bl = s5_txt(s5_P("V4_local", "source_data", "coloc_null", "BUILD_LOG.txt"))
s5_m = s5_re.search(r"CKD: (\d+) captures, (\d+) MHC-excluded, (\d+) targets, (\d+) draws", s5_bl)
s5_ckdcap, s5_ckdtg = int(s5_m.group(1)), int(s5_m.group(3))
s5_co = pd.read_csv(s5_P("V4", "source_data", "tables_from_rds", "phase3_coloc__coloc_results_final.csv"))
s5_eq = {}
for s5_lab, s5_pth, s5_keep in [("non_ckd", s5_P("v3.2", "_engine_root", "results", "v4_local", "coloc_null", "null_observed.csv"), lambda d: d.outcome != "Wuttke2019_CKD_EA"),
                                ("egfr", s5_P("v3.2", "_engine_root", "results", "repro_phase1_hg38", "coloc_null", "null_observed_egfr.csv"), lambda d: d.outcome == d.outcome)]:
    s5_o = pd.read_csv(s5_pth)
    s5_o = s5_o[s5_keep(s5_o)].merge(s5_co[["gene", "outcome_id", "PP.H4"]], left_on=["gene", "outcome"], right_on=["gene", "outcome_id"], how="left")
    assert s5_o["PP.H4"].notna().all() and s5_o.H4.notna().all(), s5_lab
    s5_eq[s5_lab] = (int((np.abs(s5_o.H4 - s5_o["PP.H4"].round(s5_PPH4_DP)) < 1e-9).sum()), len(s5_o))
s5_rt = s5_txt(s5_P("R1_stage1", "S1R1_03_decoy_method", "S1R1_decoy_rates.R"))
s5_seed = int(s5_re.search(r"set\.seed\((\d+)\); NB <- \d+", s5_rt).group(1))
s5_dm = s5_txt(s5_P("R1_stage1", "S1R1_03_decoy_method", "S1R1_DECOY_METHOD.md"))
s5_m = s5_re.search(r"In (\d+) of (\d+) CKD captures the vector is sorted lexicographically", s5_dm)
s5_srt = (int(s5_m.group(1)), int(s5_m.group(2)))
s5_v415 = [s5_run, s5_fix, s5_ckdcap, s5_ckdtg, s5_PPH4_DP, s5_eq["non_ckd"][0], s5_eq["non_ckd"][1], s5_eq["egfr"][0], s5_eq["egfr"][1], s5_seed, s5_srt[0], s5_srt[1]]
s5_o415 = s5_val("L415")
assert s5_v415 == [s5_o415["run"], s5_o415["fix"], s5_o415["ckd_captures"], s5_o415["ckd_targets"], s5_PPH4_DP] + s5_o415["kept_non_ckd_equal"] + s5_o415["egfr_equal"] + [s5_o415["bootstrap_seed"]] + s5_o415["ckd_sorted"], s5_v415
assert s5_srt[1] == s5_ckdcap


def s5_dmy(iso):
    d = s5_dt.date.fromisoformat(iso)
    return "%d %s %d" % (d.day, ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"][d.month - 1], d.year)


s5_succ("L415", s5_v415,
        "run of %s; reader corrected %s; CKD rebuilt from %d captures (%d targets); kept targets reproducing the corrected PP.H4 at %d dp: %s of %s non-CKD, %d of %d eGFR; bootstrap seed %d; CKD captures alphabetical %d of %d"
        % (s5_dmy(s5_v415[0]), s5_dmy(s5_v415[1]), s5_v415[2], s5_v415[3], s5_v415[4], "{:,}".format(s5_v415[5]), "{:,}".format(s5_v415[6]), s5_v415[7], s5_v415[8], s5_v415[9], s5_v415[10], s5_v415[11]), "",
        "Decoy provenance: date of the run that supplied the six non-CKD non-eGFR outcomes and eGFR; date of the outcome-reader correction; CKD stratum rebuilt from corrected captures (captures, targets); decimal places of the PP.H4 comparison; kept targets whose true-pairing PP.H4 equals the corrected analysis (non-CKD; eGFR); bootstrap seed; CKD captures in alphabetical order of variant identifier (sorted, of all)",
        "V3/paper3_deposit_request_v4/REPORT.md (A17 row; 'predate the ... engine fix') + V3/V4_local/source_data/coloc_null/BUILD_LOG.txt + null_observed{,_egfr}.csv vs V4 coloc table + S1R1_decoy_rates.R + S1R1_DECOY_METHOD.md",
        "mtimes; splice counts; H4 vs round(PP.H4, 3); set.seed; sort count",
        ("the same provenance stored as a list in display reading order, with the decimal places (3) and the second 168 added; the dates are ISO strings",
         "Stage 5 step 1c. Leaves: run date, fix date, 168 captures, 162 targets, 3 dp, 1,170 of 1,170, 321 of 321, seed 20261002, 166, 168 (stored twice, at leaf 2 and leaf 11, because the display reads it twice). "
         "The comparison is repeated here from the null files and the V4 coloc table at the one constant s5_PPH4_DP = 3, asserted to be the .round(3) that ext_supp.py uses; every leaf is asserted equal to the L415 object."))

# ======================================================================================================== L505: successor of L438 (kidney benchmark with the cystatin-C level)
s5_rtp = s5_P("paper3_deposit_request_v4", "scripts", "scripts", "_shared", "renal_triad", "01_renal_triad_mr.R")
s5_rl = s5_txt(s5_rtp).splitlines()
s5_m = s5_re.search(r'fifelse\(cys_p >= (0\.\d+), "mirage"', s5_rl[98])
assert s5_m, "01_renal_triad_mr.R l.99 is not the cystatin-C significance test: %r" % s5_rl[98]
s5_cysp = float(s5_m.group(1))
s5_m = s5_re.search(r"^SWEEP_P <- (0\.\d+)\s*$", s5_rl[43])
assert s5_m and float(s5_m.group(1)) == s5_cysp, "SWEEP_P (l.44) differs from the cystatin-C level (l.99)"
s5_o438 = s5_val("L438")
assert len(s5_o438) == 15 and s5_o438[8:10] == [4, 9]
s5_v438 = s5_o438[:10] + [s5_cysp] + s5_o438[10:]
s5_succ("L438", s5_v438,
        "%d not significant; %d opposite; %d in all; %d of them pass the rule; confirmed: %d of %d pass; palindromes with frequency %.2f-%.2f dropped; Effect/SE at %d dp: %d proteins span p = %g (none passes the rule), confirmed %d-%d, opposite %d-%d (both at z = %d)"
        % tuple(s5_v438), "",
        "Benchmark proteins not confirmed on cystatin C (cystatin p >= the level / opposite sign) / total; of these, passing the colocalization rule; confirmed proteins passing the rule; palindrome band dropped; decimal places of the cystatin-C Effect/StdErr; proteins whose recomputed p spans the level; confirmed count range; opposite-sign count range; z of the opposite-sign proteins",
        "as L438: egfr_mirage_classification.csv (+ UKB-PPP instrument beta); V4 phase3_coloc__coloc_results_final.csv (eGFR); renal_triad/01_renal_triad_mr.R (cystatin level: l.99 cys_p >= 0.05, SWEEP_P l.44)",
        "L438 leaves with the cystatin-C level (01_renal_triad_mr.R l.99) inserted before the confirmed range",
        ("the same leaves with the cystatin-C significance level 0.05 inserted between the 9 and the 47",
         "Stage 5 step 1c. The 0.05 is parsed from 01_renal_triad_mr.R l.99 (cys_p >= 0.05) and asserted equal to SWEEP_P (l.44). The other 15 leaves are those of L438 in stored order. The gene names stay in L438's note."))

# ======================================================================================================== L506: successor of L439 (threshold sweep with the p thresholds)
s5_sw = pd.read_csv(s5_P("v3.2", "_engine_root", "results", "repro_phase1_hg38", "pam_conditional", "egfr_mirage_threshold_sweep.csv"))
assert list(s5_sw.columns) == ["threshold_label", "threshold", "n_classifiable", "n_mirage", "pct_mirage"] and len(s5_sw) == 4
assert s5_sw.threshold_label.iloc[0].startswith("Bonferroni") and list(s5_sw.threshold_label.iloc[1:]) == ["p<1e-4", "p<1e-3", "p<0.05"]
s5_v439 = []
for s5_i in range(4):
    s5_v439 += [float(s5_sw.threshold.iloc[s5_i]), int(s5_sw.n_mirage.iloc[s5_i]), int(s5_sw.n_classifiable.iloc[s5_i])]
s5_o439 = s5_val("L439")
assert [s5_v439[0]] + [s5_v439[i] for i in (1, 2, 4, 5, 7, 8, 10, 11)] == s5_o439, (s5_v439, s5_o439)
assert [float(x) for x in (1e-4, 1e-3, 0.05)] == [s5_v439[3], s5_v439[6], s5_v439[9]]


def s5_thr(x):
    return s5_re.sub(r"e-0?(\d)", r"e-\1", "%.0e" % x) if x < 0.01 else "%g" % x


s5_succ("L439", s5_v439,
        "Bonferroni (%.1e) %d of %d; p < %s %d of %d; p < %s %d of %d; p < %s %d of %d"
        % (s5_v439[0], s5_v439[1], s5_v439[2], s5_thr(s5_v439[3]), s5_v439[4], s5_v439[5], s5_thr(s5_v439[6]), s5_v439[7], s5_v439[8], s5_thr(s5_v439[9]), s5_v439[10], s5_v439[11]), "",
        "Benchmark threshold sweep: at each creatinine-eGFR p threshold (Bonferroni / 1e-4 / 1e-3 / 0.05), proteins with cystatin-C p >= 0.05 of those classifiable (threshold, n_mirage, n_classifiable)",
        "v3.2/_engine_root/results/repro_phase1_hg38/pam_conditional/egfr_mirage_threshold_sweep.csv (threshold, n_mirage, n_classifiable)",
        "threshold, n_mirage, n_classifiable (one triple per sweep row)",
        ("the same sweep stored as a threshold, count and total for each of the four thresholds, so every threshold the display prints derives from the value",
         "Stage 5 step 1c. L439 stored the Bonferroni cut and the counts only; the display also prints p < 1e-4, p < 1e-3 and p < 0.05. All twelve leaves are read from the sweep file's threshold, n_mirage and n_classifiable "
         "columns and asserted equal to the L439 leaves they share. The display prints the cut at two significant figures (3.2e-06) while Methods and the SM print it as 3.21 x 10^-6 [[L242]]: one stored value at two precisions."))

# ======================================================================================================== L507: successor of L446 (Supplementary Note with the thresholds)
s5_tc = pd.read_csv(s5_P("R1_stage1", "S1R1_02_topology_note", "S1R1_topo_supp_table_c_excess.tsv"), sep="\t")
s5_MS = "random SNP sets matched on LD score and MAF"
s5_mrow = lambda thr: s5_tc[(s5_tc.analysis_set == s5_MS) & (s5_tc.threshold == thr) & (s5_tc.n_sets == 30)].iloc[0]
s5_p0 = s5_tc[(s5_tc.analysis_set == "proteins (cis-pQTL MR)") & (s5_tc.unit == "protein") & (s5_tc.threshold == 0.01) & (s5_tc.null_model == "within-unit permutation") & (s5_tc.category == "all")].iloc[0]
s5_lp = s5_tc[s5_tc.analysis_set == "proteins with an LD-scored lead instrument"]
assert int(s5_p0.n_replicates) == 5000 and int(s5_p0.observed_count) == 128 and s5_lp.n_units.nunique() == 1
s5_pp = s5_tc[(s5_tc.analysis_set == "proteins (cis-pQTL MR)") & (s5_tc.unit == "protein") & (s5_tc.threshold == 0.01) & (s5_tc.category == "all") & (s5_tc.observed_count == 128)]
s5_corr = s5_pp[~s5_pp.finngen_noise_correlation.str.startswith(("none", "LDSC intercepts, diagonal only"))]
s5_ind = s5_pp[(s5_pp.null_model == "signal + correlated noise model") & s5_pp.finngen_noise_correlation.str.startswith("none")]
assert len(s5_corr) == 11 and len(s5_ind) == 1 and (s5_corr.n_replicates == 20000).all()
s5_lr = lambda thr: s5_lp[s5_lp.threshold == thr].iloc[0]
s5_o446 = s5_val("L446")
s5_groups = s5_o446[1]                                           # the seven outcome groups: a definition of the count, read from L446 (leaf [1])
s5_a, s5_b = s5_mrow(0.01), s5_mrow(0.005)
s5_v446 = [float(s5_p0.threshold), s5_groups, int(s5_p0.n_replicates), int(s5_lp.n_units.iloc[0]), int(s5_a.n_sets),
           float(s5_a.fold_min), float(s5_a.fold_max), float(s5_a.threshold), float(s5_b.fold_min), float(s5_b.fold_max), float(s5_b.threshold),
           float(s5_p0.null_mean), float(s5_corr.null_mean.min()), float(s5_corr.null_mean.max()), float(s5_ind.null_mean.iloc[0]),
           float(s5_lr(0.005).fold_excess), float(s5_lr(0.005).hit_rate), float(s5_lr(0.001).fold_excess), float(s5_lr(0.001).hit_rate),
           float(s5_a.fold_excess), float(s5_a.hit_rate), float(s5_b.fold_excess), float(s5_b.hit_rate), int(s5_p0.observed_count)]
assert s5_v446[:7] + s5_v446[8:10] + s5_v446[11:] == s5_o446[:7] + s5_o446[7:9] + s5_o446[9:], "successor leaves differ from L446's"
assert len(s5_v446) == 24 and len(s5_o446) == 22
s5_succ("L446", s5_v446,
        "p < %g; %d groups; %d permutations; %d proteins; %d sets; fold range %.2f-%.2f (p < %g), %.2f-%.2f (p < %g); null mean %.1f (permutation), %.1f-%.1f (correlated noise), %.1f (independent noise); "
        "proteins fold %.2f (hit %.1f%%), %.2f (%.1f%%); random fold %.2f (hit %.1f%%), %.2f (%.1f%%); observed %d"
        % tuple(s5_v446[:15] + [s5_v446[15], 100 * s5_v446[16], s5_v446[17], 100 * s5_v446[18], s5_v446[19], 100 * s5_v446[20], s5_v446[21], 100 * s5_v446[22], s5_v446[23]]), "",
        "Multi-outcome count (proteins with p < 0.01 for >= 2 of 7 outcome groups, eGFR and CKD as one): threshold; groups; within-unit permutations; LD-scored proteins; matched random-variant sets; matched-set fold range at p < 0.01 and 0.005 (with the two thresholds); null means; proteins' fold and hit rate at p < 0.005 and 0.001; matched random sets' mean fold and hit rate at p < 0.01 and 0.005; observed count",
        "V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_supp_table_c_excess.tsv", "n_replicates, n_units, n_sets, threshold, fold_min, fold_max, null_mean, fold_excess, hit_rate",
        ("the same Note numbers with the two thresholds that the display prints beside the matched-set fold ranges (0.01 and 0.005) stored in reading order",
         "Stage 5 step 1c. 24 leaves in display reading order: L446's 22 leaves with the second 0.01 inserted after the upper fold of the p < 0.01 range and 0.005 after the upper fold of the p < 0.005 range. "
         "Both are the threshold column of the two matched-random-set rows (30 sets); the first 0.01 is the threshold of the observed-count row; 7 groups as L446 leaf [1]. Re-derived here from the table and asserted equal to L446's leaves. "
         "Flagged for the prose step, not a D1 matter: the Note's hit rates for the proteins (p < 0.005, 0.001) and the matched random sets (p < 0.01, 0.005) are compared without stating the four thresholds."))

# ======================================================================================================== S5N01-S5N13: placeholder rows
s5_cache = {}
s5_N = {}
for s5_r in s5_NSPEC:
    s5_path = s5_P(*s5_r["source_file"].split("/"))
    if s5_path not in s5_cache:
        assert s5_md5(s5_path) == s5_r["source_md5"], "source md5 differs from the spec (%s): %s; do not re-pin, report" % (s5_r["source_file"], s5_md5(s5_path))
        s5_cache[s5_path] = s5_io.open(s5_path, encoding="utf-8", newline="").read().split("\n")
    s5_lines = s5_cache[s5_path]
    s5_hdr = next(s5_csv.reader([s5_lines[0]], delimiter="\t"))
    s5_ln = int(s5_r["source_line"])
    s5_rec = dict(zip(s5_hdr, next(s5_csv.reader([s5_lines[s5_ln - 1]], delimiter="\t"))))
    assert s5_rec["quantity_id"] == s5_r["source_quantity_id"], (s5_r["placeholder"], s5_ln, s5_rec["quantity_id"])
    assert s5_r["source_column"] == "value" and s5_rec["value"] != ""
    s5_x = s5_rec["value"]
    s5_v = int(s5_x) if s5_re.fullmatch(r"-?\d+", s5_x) else float(s5_x)
    s5_disp = ("%d" % s5_v) if isinstance(s5_v, int) else ("%.2f" % s5_v)
    assert s5_disp == s5_r["printed_display"], (s5_r["placeholder"], s5_disp, s5_r["printed_display"])
    if isinstance(s5_v, int):
        assert str(s5_v) == s5_r["value_read_from_source"] and s5_rec["unit"] == "cells"
    else:
        assert float(s5_r["value_read_from_source"]) == s5_v and s5_rec["unit"] == "PP.H4"
    s5_cid = s5_alloc(s5_r["placeholder"])
    s5_N[s5_r["placeholder"]] = s5_cid
    s5_ana = s5_r["placeholder"] != "S5N01"
    add(s5_cid, "SM: " + ("gene attribution (restricted grid)" if s5_ana else "gene attribution (conditional colocalization)") + " (placeholder %s)" % s5_r["placeholder"], s5_r["claim"], s5_v, s5_disp,
        "V3/%s; line %d; quantity id %s; md5 %s" % (s5_r["source_file"], s5_ln, s5_r["source_quantity_id"], s5_r["source_md5"]), "quantity_id=%s; column value" % s5_r["source_quantity_id"], "[R]",
        note=("Stage 5 step 1c, placeholder {{%s}} of %s. " % (s5_r["placeholder"], s5_r["surfaces_to_resolve_in_1c"].split(" (")[0]))
        + ("The branch flags the W4b grid quantities as analyst-computed, not independently re-derived (Rule 7, middle-path DEVIATION); supplement positions only." if s5_ana else
           "Rule 3(iv) clause input; the 489-panel value of the same quantity (0.9985675, S5.coloc_clause.A_rs409170.489.PPH4) is not stored here: it needs its own row if it is ever printed."))

assert s5_next[0] == 521 and len(s5_idmap) == 21, (s5_next, len(s5_idmap))
print("ext_r2_stage5 id map: " + "; ".join("%s=%s" % kv for kv in s5_idmap))
