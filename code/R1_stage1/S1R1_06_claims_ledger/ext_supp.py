# -*- coding: utf-8 -*-
"""Stage 3b track C1 (Supplementary Methods, part 1): ledger rows L400-L424.

Executed by S1R1_build_ledger.py in ITS namespace (add, rel, pd, np, V3, CKM, S1, io, os, json are available).
Every name defined here starts with `xs_` so nothing of the builder is overwritten (the builder's own loop uses
_ef/_lo/_hi/_ep/_n0 after this file runs). Every value is read from a named source file and asserted; no value is typed
unless an assert pins it to the text of its source.
Rows feed: supplement/P31R_SM_{multitest,candidates,decoy,instruments,grading,phewas,devstage}.md.
"""
import re as xs_re
from scipy.stats import norm as xs_norm

xs_DEP = os.path.join(V3, "paper3_deposit_request_v4", "scripts")
xs_LOGS = os.path.join(xs_DEP, "reports", "Paper3_Program", "run_logs", "RERUN_20260915")


def xs_txt(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


def xs_has(text, *needles):
    for n_ in needles:
        assert n_ in text, "source text not found: %r" % n_


# ---------------------------------------------------------------------------------------------------- tables
xs_MR = os.path.join(V3, "V4", "source_data", "tables_from_rds", "phase2_mr_screen__mr_all_results.csv")
xs_mr = pd.read_csv(xs_MR)
xs_prim = xs_mr[xs_mr.method.isin(["Wald ratio", "Inverse variance weighted"])].copy()
xs_fam = 1948 * 8
xs_INS = os.path.join(V3, "V4", "source_data", "tables_from_rds", "phase1_instruments__ukb_cis_instruments.csv")
xs_ins = pd.read_csv(xs_INS)
assert xs_ins.gene.nunique() == 1948 and xs_fam == 15584
xs_TIER = os.path.join(V3, "v3.2", "_engine_root", "results", "phase4_annotation", "table_s_pair_evidence_tiered.csv")
xs_t = pd.read_csv(xs_TIER)
assert len(xs_t) == 81
xs_CA = os.path.join(V3, "V4", "source_data", "supp_tables", "SuppTable_candidate_edge_audit.csv")
xs_ca = pd.read_csv(xs_CA)

# ---------------------------------------------------------------------------------------------------- L400 multiple-testing family
xs_niv = int((xs_prim.method == "Inverse variance weighted").sum())
xs_nwr = int((xs_prim.method == "Wald ratio").sum())
assert xs_prim[["gene", "outcome_name"]].drop_duplicates().shape[0] == len(xs_prim) == xs_niv + xs_nwr
xs_unt = xs_fam - len(xs_prim)
add("L400", "SM: multiple testing", "Multiple-testing family: primary inverse-variance-weighted estimates (one per multi-instrument combination); untestable combinations entered at p = 1",
    {"ivw": xs_niv, "untestable": xs_unt, "pad_p": 1},
    "%d inverse-variance-weighted estimates; %s untestable combinations entered at p = 1" % (xs_niv, "{:,}".format(xs_unt)),
    rel(xs_MR), "method in {Wald ratio, IVW}; 15,584 - testable rows", "[D]",
    note="Padding at p = 1 is the definition in R/00_config.R primary_fdr_set() (asserted below).")
xs_has(xs_txt(os.path.join(xs_DEP, "R", "00_config.R")), "untested cells padded at", "p = 1")

# ---------------------------------------------------------------------------------------------------- L401 stored flag vs family-wide BH
xs_p = np.sort(np.r_[xs_prim.pval.values, np.ones(xs_unt)])
xs_k = np.arange(1, xs_fam + 1)
xs_thr = float(xs_p[(xs_p <= xs_k / xs_fam * 0.05).nonzero()[0].max()])
xs_prim["xs_bh"] = xs_prim.pval <= xs_thr
assert int(xs_prim.xs_bh.sum()) == 292


def xs_bh_within(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return q < 0.05


xs_within = np.zeros(len(xs_prim), bool)
for xs_meth in ["Wald ratio", "Inverse variance weighted"]:
    xs_ix = (xs_prim.method == xs_meth).values
    xs_within[xs_ix] = xs_bh_within(xs_prim.pval.values[xs_ix])
assert (xs_within == xs_prim.fdr_sig.values.astype(bool)).all(), "stored fdr_sig is not BH within each estimator"
assert int(xs_within.sum()) == 304 and int(xs_within[(xs_prim.method == "Inverse variance weighted").values].sum()) == 0
xs_ex = xs_prim[xs_within & ~xs_prim.xs_bh.values]
assert len(xs_ex) == 12 and (xs_ex.method == "Wald ratio").all() and not (xs_prim.xs_bh.values & ~xs_within).any()
xs_pairs81 = set(zip(xs_t.gene, xs_t.outcome))
xs_ex81 = int(sum((g, o) in xs_pairs81 for g, o in zip(xs_ex.gene, xs_ex.outcome_name)))
# second-order-SE pool (L263's 248): recomputed here independently of the builder's code path. The pipeline (28_candidate_pool.log)
# gave second-order SEs to all 14,423 Wald rows; this recomputation gives them only to the 14,402 Wald rows of single-instrument
# genes, because which instrument entered the other 21 is not in the table. Both give 248 (the C1 verifier's R, which uses the
# first-listed instrument for those 21, also gives 248), so the 21 rows decide no membership.
xs_ni = xs_ins.groupby("gene").size()
xs_sg = xs_ins[xs_ins.gene.map(xs_ni) == 1].set_index("gene")
xs_w = ((xs_prim.method == "Wald ratio") & xs_prim.gene.isin(xs_sg.index)).values
assert int(xs_w.sum()) == 14402 and int(((xs_prim.method == "Wald ratio") & ~xs_prim.gene.isin(xs_sg.index)).sum()) == 21
xs_bx = xs_prim.loc[xs_w, "gene"].map(xs_sg["beta"]).values
xs_sx = xs_prim.loc[xs_w, "gene"].map(xs_sg["se"]).values
xs_p2 = xs_prim.pval.values.copy()
xs_b = xs_prim.loc[xs_w, "b"].values; xs_s1 = xs_prim.loc[xs_w, "se"].values
xs_p2[xs_w] = 2 * xs_norm.sf(np.abs(xs_b / np.sqrt(xs_s1 ** 2 + xs_b ** 2 * (xs_sx / xs_bx) ** 2)))
xs_q2 = np.sort(np.r_[xs_p2, np.ones(xs_unt)])
xs_t2 = float(xs_q2[(xs_q2 <= xs_k / xs_fam * 0.05).nonzero()[0].max()])
xs_in2 = xs_p2 <= xs_t2
assert int(xs_in2.sum()) == 248
xs_expool = int((xs_in2 & xs_within & ~xs_prim.xs_bh.values).sum())
xs_tb = xs_t.merge(xs_prim[["gene", "outcome_name", "xs_bh", "fdr_sig"]], left_on=["gene", "outcome"], right_on=["gene", "outcome_name"], how="left")
assert xs_tb.xs_bh.notna().all()
xs_agree = int((xs_tb.xs_bh.astype(bool) == xs_tb.mr_sig_class.isin(["Bonferroni", "FDR"])).sum())
xs_both = int((xs_tb.xs_bh.astype(bool) & xs_tb.fdr_sig.astype(bool)).sum())
assert xs_agree == 81 and xs_both == 80
xs_dp = float(xs_t[(xs_t.gene == "DPEP1") & (xs_t.outcome == "AF")].mr_pval.iloc[0])
assert not xs_tb[(xs_tb.gene == "DPEP1") & (xs_tb.outcome == "AF")].xs_bh.iloc[0]
add("L401", "SM: multiple testing", "Family-wide BH (m = 15,584) vs the stored per-estimator flag fdr_sig: largest p admitted; stored count; extra tests and their p range; extras among the 81 / in the second-order candidate pool; the 81 pairs' FDR class identical under both (agree, pass both); DPEP1-AF p",
    {"max_p_admitted": xs_thr, "stored": int(xs_within.sum()), "extra": len(xs_ex), "extra_p_min": float(xs_ex.pval.min()), "extra_p_max": float(xs_ex.pval.max()),
     "extra_in_81": xs_ex81, "extra_in_pool": xs_expool, "agree_81": xs_agree, "pass_both": xs_both, "dpep1_af_p": xs_dp},
    "largest p admitted %.2e; stored flag %d; %d extra (p %.2e to %.2e); %d among the 81, %d in the candidate pool; %d agree, %d pass both; DPEP1-AF p = %.3f"
    % (xs_thr, int(xs_within.sum()), len(xs_ex), xs_ex.pval.min(), xs_ex.pval.max(), xs_ex81, xs_expool, xs_agree, xs_both, xs_dp),
    rel(xs_MR) + " + " + rel(xs_TIER), "BH(p padded to 15,584) vs p.adjust(BH) within method (01_mr_screen.R fdr_sig line)", "[D]",
    note="Asserted: the stored fdr_sig equals BH within each estimator over testable rows; 12 extra Wald-ratio tests, none colocalization-supported and none in the 248 pool.")
xs_has(xs_txt(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "cascade", "01_mr_screen.R")), 'mr_all[, fdr_sig := p.adjust(pval, method="BH") < FDR_ALPHA, by=method]')

# ---------------------------------------------------------------------------------------------------- L402 second-order SEs (pipeline run log)
xs_cpl = xs_txt(os.path.join(xs_LOGS, "28_candidate_pool.log"))
xs_m = xs_re.search(r"second-order SE on (\d+) single-instrument row\(s\); first-order kept for (\d+) multi-instrument row\(s\)", xs_cpl)
assert xs_m and int(xs_m.group(1)) == xs_nwr == 14423
xs_mo = int(xs_m.group(2))
assert xs_mo == int((~xs_mr.method.isin(["Wald ratio"])).sum()) == 114   # 98 IVW + 8 Egger + 8 weighted median rows
add("L402", "SM: multi-signal candidates", "Candidate-pool FDR: Wald ratios given second-order SEs / multi-instrument rows kept first-order (98 of them the primary IVW estimates)",
    {"wald_second_order": xs_nwr, "multi_rows_first_order": xs_mo, "ivw_primary": xs_niv},
    "%s Wald ratios with second-order SEs; %d multi-instrument rows first-order (%d primary IVW)" % ("{:,}".format(xs_nwr), xs_mo, xs_niv),
    rel(os.path.join(xs_LOGS, "28_candidate_pool.log")), "'second-order SE on ... first-order kept for ...'", "[D]")

# ---------------------------------------------------------------------------------------------------- L403 SuSiE candidate funnel
xs_m = xs_re.search(r"testable\(SuSiE\)=(\d+) \[(\d+) eGFR, (\d+) non-eGFR\] \| MHC=(\d+)", xs_cpl)
xs_test, xs_te, xs_tn, xs_mhc = [int(x) for x in xs_m.groups()]
xs_aud = xs_txt(os.path.join(xs_LOGS, "30_candidate_audit.log"))
xs_resc = int(xs_re.search(r"susie_rescued\s+(\d+)", xs_aud).group(1))
xs_notsup = int(xs_re.search(r"susie_not_supportive\s+(\d+)", xs_aud).group(1))
xs_mhc2 = int(xs_re.search(r"MHC_excluded\s+(\d+)", xs_aud).group(1))
xs_rn = int(xs_re.search(r"29 \[nonegfr\]: \d+ run \| ok=\d+ \| skip_benign=0 \| rescued\(PP.H4_susie>=0.80\)=(\d+)", xs_txt(os.path.join(xs_LOGS, "29_cand_nonegfr.log"))).group(1))
xs_re_ = int(xs_re.search(r"29 \[egfr\]: \d+ run \| ok=\d+ \| skip_benign=0 \| rescued\(PP.H4_susie>=0.80\)=(\d+)", xs_txt(os.path.join(xs_LOGS, "29_cand_egfr.log"))).group(1))
assert xs_resc == len(xs_ca) == xs_rn + xs_re_ == 26 and xs_resc + xs_notsup == xs_test == 160 and xs_mhc == xs_mhc2 == 17
xs_st = xs_ca.Status.value_counts().to_dict()
xs_drop = xs_ca[xs_ca.Status == "Dropped"]
xs_below = int(xs_drop.Reason.str.contains("< 0.80 under purity filter").sum())
xs_nopure = int(xs_drop.Reason.str.contains("no pure credible set").sum())
assert (xs_st["Drawn"], xs_st["Held (not drawn)"], xs_st["Dropped"]) == (16, 1, 9) and xs_below + xs_nopure == 9
xs_held = xs_ca[xs_ca.Status == "Held (not drawn)"]
xs_mi = int(xs_re.search(r"max_iter=(\d+)", xs_held.Reason.iloc[0]).group(1))
assert list(zip(xs_held.Gene, xs_held.Endpoint)) == [("ZFYVE19", "eGFR")]
add("L403", "SM: multi-signal candidates", "SuSiE candidate funnel over the 160 testable associations: first fit PP.H4 >= 0.80 (non-eGFR / eGFR run) / not; purity re-run kept / not converged / below 0.80 / no pure credible set; IBSS iteration cap; testable by run (eGFR / non-eGFR)",
    {"first_fit_ge080": xs_resc, "nonegfr_run": xs_rn, "egfr_run": xs_re_, "not_supportive": xs_notsup, "kept": 16, "not_converged": 1, "below_080": xs_below,
     "no_pure_cs": xs_nopure, "max_iter": xs_mi, "testable_egfr": xs_te, "testable_nonegfr": xs_tn},
    "%d of %d reached PP.H4 >= 0.80 on the first fit (%d non-eGFR, %d eGFR); %d did not; kept 16, 1 not converged (ZFYVE19-eGFR), %d below 0.80, %d no pure credible set; IBSS cap %s iterations"
    % (xs_resc, xs_test, xs_rn, xs_re_, xs_notsup, xs_below, xs_nopure, "{:,}".format(xs_mi)),
    rel(xs_CA) + " + run logs 29_cand_*.log, 30_candidate_audit.log", "Status, Reason; rescued counts", "[D]")

# ---------------------------------------------------------------------------------------------------- L404 MHC bounds and locus window (candidate pool / kidney comparison)
xs_cfg = xs_txt(os.path.join(xs_DEP, "R", "00_config.R"))
xs_m = xs_re.search(r'MHC_CHR <- "(\d+)"; MHC_START <- (\d+)L; MHC_END <- (\d+)L', xs_cfg)
xs_chr, xs_lo, xs_hi = xs_m.group(1), int(xs_m.group(2)), int(xs_m.group(3))
xs_g1 = xs_ins.drop_duplicates("gene").set_index("gene")
xs_poolp = xs_prim[xs_in2]
xs_cp = set(zip(*[pd.read_csv(os.path.join(V3, "V4", "source_data", "tables_from_rds", "phase3_coloc__coloc_results_final.csv")).query("`PP.H4` >= 0.8 and `PP.H4` >= 5 * `PP.H3`")[c] for c in ["gene", "outcome_name"]]))
xs_fail = xs_poolp[[(g, o) not in xs_cp for g, o in zip(xs_poolp.gene, xs_poolp.outcome_name)]].copy()
xs_fail["mhc"] = [str(xs_g1.loc[g, "gene_chrom"]).replace("chr", "") == xs_chr and xs_g1.loc[g, "gene_start"] <= xs_hi and xs_g1.loc[g, "gene_end"] >= xs_lo for g in xs_fail.gene]
assert len(xs_fail) == 177 and int(xs_fail.mhc.sum()) == 17
xs_tst = xs_fail[~xs_fail.mhc][["gene", "outcome_name"]].copy()
xs_dr = set(zip(xs_ca[xs_ca.Status == "Drawn"].Gene, xs_ca[xs_ca.Status == "Drawn"].Endpoint))
xs_tst["grp"] = np.where(xs_tst.outcome_name.isin(["eGFR", "CKD"]), "kidney", xs_tst.outcome_name)
xs_tst["chr"] = [str(xs_g1.loc[g, "gene_chrom"]).replace("chr", "") for g in xs_tst.gene]
xs_tst["pos"] = [float(xs_g1.loc[g, "gene_start"]) for g in xs_tst.gene]
xs_tst["cand"] = [(g, o) in xs_dr for g, o in zip(xs_tst.gene, xs_tst.outcome_name)]
xs_WIN = 1e6
xs_tst = xs_tst.sort_values(["grp", "chr", "pos"]).reset_index(drop=True)
xs_lid, xs_c = [], 0
for xs_i in range(len(xs_tst)):
    if xs_i == 0 or xs_tst.grp[xs_i] != xs_tst.grp[xs_i - 1] or xs_tst.chr[xs_i] != xs_tst.chr[xs_i - 1] or xs_tst.pos[xs_i] - xs_tst.pos[xs_i - 1] > xs_WIN:
        xs_c += 1
    xs_lid.append(xs_c)
xs_tst["locus"] = xs_lid
xs_lc = xs_tst.groupby("locus").agg(kid=("grp", lambda x: x.iloc[0] == "kidney"), cand=("cand", "any"))
xs_loc = (int((xs_lc.cand & xs_lc.kid).sum()), int(xs_lc.cand.sum()), int((~xs_lc.cand & xs_lc.kid).sum()), int((~xs_lc.cand).sum()))
assert len(xs_tst) == 160 and int(xs_tst.cand.sum()) == 16 and xs_loc == (9, 14, 44, 120), xs_loc
add("L404", "SM: multi-signal candidates", "Candidate-pool definitions: extended MHC (GRCh38 chromosome, start Mb, end Mb; gene overlaps) and the locus window for the locus-level kidney comparison (Mb; single linkage of gene starts)",
    {"chrom": int(xs_chr), "mhc_start_mb": xs_lo / 1e6, "mhc_end_mb": xs_hi / 1e6, "locus_window_mb": xs_WIN / 1e6},
    "MHC chromosome %s, %d-%d Mb; loci joined within %d Mb" % (xs_chr, xs_lo / 1e6, xs_hi / 1e6, xs_WIN / 1e6), rel(os.path.join(xs_DEP, "R", "00_config.R")) + " (MHC_CHR/MHC_START/MHC_END)",
    "re-run here: 17 MHC of 177; loci 9/14 vs 44/120 (equal to L263/L264)", "[D]")

# ---------------------------------------------------------------------------------------------------- L405 decoy design constants
xs_nul = xs_txt(os.path.join(xs_DEP, "scripts", "_shared", "coloc_null", "05_run_null.R"))
xs_seed = int(xs_re.search(r"set\.seed\((\d+)\)", xs_nul).group(1))
xs_B = int(xs_re.search(r'B <- as\.integer\(Sys\.getenv\("NPERM", "(\d+)"\)\)', xs_nul).group(1))
xs_kmin = int(xs_re.search(r"if\(k<(\d+)\) return\(NULL\)", xs_nul).group(1))
xs_m = xs_re.search(r'cc\$c=="(\d+)" && cc\$s>=([0-9.]+)e6 && cc\$s<=([0-9.]+)e6', xs_nul)
xs_rates = xs_txt(os.path.join(S1, "S1R1_03_decoy_method", "S1R1_decoy_rates.R"))
xs_nb = int(xs_re.search(r"NB <- (\d+)", xs_rates).group(1))
xs_has(xs_nul, "dr <- sample(elig, B, replace=TRUE)", "abs(cj$s-gci$s)>5e6", "!is_mhc(cj)")
add("L405", "SM: decoy-locus stress test", "Decoy design: donors drawn per target (with replacement); random seed; minimum shared rows; MHC exclusion (chromosome, gene-start window Mb); bootstrap resamples for the intervals",
    {"donors": xs_B, "seed": xs_seed, "min_rows": xs_kmin, "chrom": int(xs_m.group(1)), "mhc_start_mb": float(xs_m.group(2)), "mhc_end_mb": float(xs_m.group(3)), "bootstrap": xs_nb},
    "%d donors per target; seed %d; at least %d rows; MHC chromosome %s, %s-%s Mb (gene start); %s bootstrap resamples" % (xs_B, xs_seed, xs_kmin, xs_m.group(1), xs_m.group(2), xs_m.group(3), "{:,}".format(xs_nb)),
    rel(os.path.join(xs_DEP, "scripts", "_shared", "coloc_null", "05_run_null.R")) + " + V3/R1_stage1/S1R1_03_decoy_method/S1R1_decoy_rates.R", "set.seed, NPERM default, k<, is_mhc, NB", "[D]")

# ---------------------------------------------------------------------------------------------------- L406 decoy strata
xs_drt = pd.read_csv(os.path.join(S1, "S1R1_03_decoy_method", "S1R1_decoy_rates.tsv"), sep="\t").set_index("stratum")
xs_ne, xs_eg = xs_drt.loc["non-eGFR (seven outcomes, CKD included)"], xs_drt.loc["eGFR"]
assert int(xs_ne.draws) == xs_B * int(xs_ne.targets_tested) and int(xs_eg.draws) == xs_B * int(xs_eg.targets_tested)
xs_po = xs_drt.loc[[s for s in xs_drt.index if s.startswith(("FinnGen_", "MEGASTROKE", "Wuttke2019_CKD"))]]
assert len(xs_po) == 7 and int(xs_po.draws.sum()) == int(xs_ne.draws)
xs_sel = "CKD, selection-matched, p < 0.05 at the located instrument (sel_p; 154 contributing targets)"
assert int(xs_drt.loc[xs_sel, "draws"]) == 1535
xs_lo_o, xs_hi_o = xs_po.rate.idxmin(), xs_po.rate.idxmax()
assert xs_lo_o == "FinnGen_R10_NAFLD" and xs_hi_o == "FinnGen_R10_T2D"
add("L406", "SM: decoy-locus stress test", "Decoy targets (non-eGFR / eGFR); per-outcome pass-rate range among the seven non-eGFR outcomes (NAFLD lowest, T2D highest); observed true-pairing passes (non-eGFR / eGFR); two-way bootstrap 95% intervals at full precision (non-eGFR, eGFR, CKD selection-matched)",
    {"targets_non_egfr": int(xs_ne.targets_tested), "targets_egfr": int(xs_eg.targets_tested), "rate_min": float(xs_po.rate.min()), "rate_max": float(xs_po.rate.max()),
     "observed_non_egfr": int(xs_ne.observed_passes), "observed_egfr": int(xs_eg.observed_passes),
     "twoway_non_egfr": [float(xs_ne.twoway_boot_lo), float(xs_ne.twoway_boot_hi)], "twoway_egfr": [float(xs_eg.twoway_boot_lo), float(xs_eg.twoway_boot_hi)],
     "twoway_ckd_selected": [float(xs_drt.loc[xs_sel, "twoway_boot_lo"]), float(xs_drt.loc[xs_sel, "twoway_boot_hi"])]},
    "%s / %d targets; per-outcome %.2f%% (NAFLD) to %.2f%% (T2D); two-way 95%% intervals %.2f-%.2f%% (non-eGFR), %.2f-%.2f%% (eGFR), %.2f-%.2f%% (CKD selection-matched)"
    % ("{:,}".format(int(xs_ne.targets_tested)), int(xs_eg.targets_tested), 100 * xs_po.rate.min(), 100 * xs_po.rate.max(),
       100 * xs_ne.twoway_boot_lo, 100 * xs_ne.twoway_boot_hi, 100 * xs_eg.twoway_boot_lo, 100 * xs_eg.twoway_boot_hi,
       100 * xs_drt.loc[xs_sel, "twoway_boot_lo"], 100 * xs_drt.loc[xs_sel, "twoway_boot_hi"]),
    "V3/R1_stage1/S1R1_03_decoy_method/S1R1_decoy_rates.tsv", "targets_tested, rate by outcome", "[R][S]")

# ---------------------------------------------------------------------------------------------------- L407 instrument clumping, proxies, palindromes
xs_scr = xs_txt(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "cascade", "01_mr_screen.R"))
xs_ldkb = int(xs_re.search(r"^LD_KB\s*<-\s*(\d+)", xs_cfg, xs_re.M).group(1))
xs_has(xs_cfg, "Phase 1 selects by DISTANCE-based greedy clumping in a", "1,000 kb window")
xs_m = xs_re.search(r'"--ld-window-kb","(\d+)","--ld-window","\d+","--ld-window-r2","([0-9.]+)"', xs_scr)
xs_pkb, xs_pr2 = int(xs_m.group(1)), float(xs_m.group(2))
xs_m2 = xs_re.search(r"\.unin <- function\(f\) is\.na\(f\) \|\| \(f > ([0-9.]+) && f < ([0-9.]+)\)", xs_scr)
xs_has(xs_scr, "cands <- pr_av[instrument==g][order(-r2)]", "(dsub$eaf > 0.42 & dsub$eaf < 0.58)")
assert xs_ldkb == 1000 and xs_pkb == 1000 and xs_pr2 == 0.8
add("L407", "SM: instruments", "Instrument selection and harmonisation parameters: distance-clumping window (kb); LD proxy minimum r2 and window (kb), highest r2 first; palindromic frequency treated as uninformative inside (lo, hi)",
    {"clump_kb": xs_ldkb, "proxy_r2": xs_pr2, "proxy_kb": xs_pkb, "pal_lo": float(xs_m2.group(1)), "pal_hi": float(xs_m2.group(2))},
    "distance clumping %s kb; proxies r2 >= %.1f within %s kb; palindromes uninformative at %s-%s" % ("{:,}".format(xs_ldkb), xs_pr2, "{:,}".format(xs_pkb), xs_m2.group(1), xs_m2.group(2)),
    rel(os.path.join(xs_DEP, "R", "00_config.R")) + " + " + rel(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "cascade", "01_mr_screen.R")), "LD_KB; PLINK --r2 args; .unin()", "[D]")

# ---------------------------------------------------------------------------------------------------- L408 instruments lost at p < 5e-8
xs_lg = ["AK1", "BRAP", "CCND2", "FLT1", "KCTD5"]
xs_l = xs_ins[xs_ins.gene.isin(xs_lg)]
assert len(xs_l) == 5 and (xs_l.pval >= 5e-8).all() and (xs_l.pval < 1e-6).all()
xs_lt = xs_t.merge(xs_ins[["gene", "pval"]], on="gene")
assert sorted(xs_lt[xs_lt.pval >= 5e-8].gene.unique()) == xs_lg and int((xs_lt.pval >= 5e-8).sum()) == 11
assert xs_l.set_index("gene").pval.idxmin() == "BRAP" and xs_l.set_index("gene").pval.idxmax() == "CCND2"
add("L408", "SM: instruments", "Instruments of the 11 pairs lost at p < 5e-8: BRAP lead p (lowest), CCND2 lead p (highest), F minimum and maximum",
    {"brap_p": float(xs_l.pval.min()), "ccnd2_p": float(xs_l.pval.max()), "F_min": round(float(xs_l.F_stat.min()), 1), "F_max": round(float(xs_l.F_stat.max()), 1)},
    "lead p %.1e (BRAP) to %.1e (CCND2); F %.1f to %.1f" % (xs_l.pval.min(), xs_l.pval.max(), xs_l.F_stat.min(), xs_l.F_stat.max()),
    rel(xs_INS), "pval, F_stat for AK1/BRAP/CCND2/FLT1/KCTD5", "[D]")

# ---------------------------------------------------------------------------------------------------- L409 grading: thresholds, triggers, the rule re-derived
xs_thm = xs_txt(os.path.join(xs_DEP, "R", "theme_paper3.R"))
xs_rb = xs_txt(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "phase4", "12_rebuild_pair_evidence.R"))
xs_has(xs_thm, "PP.H4_p12_1e6 < 0.5)", 'is_assay_discordant     = (decode_status %in% "ASSAY_DISCORDANT")', "crit_abf_dual)]",
       'pairs[, all_strong_criteria := (core_complete & n_evaluable >= 1L)]', 'fifelse(has_major_weakness, "Flagged",')
xs_has(xs_rb, 'if (h4 >= 0.80) return("supportive")', 'if (!is.na(h3) && h3 > 0.50) return("discordant")')
xs_opp = ((xs_t.decode_status == "NOT_REPLICATED") & xs_t.decode_evaluable_pair) | ((xs_t.aric_status == "NOT_REPLICATED") & xs_t.aric_evaluable_pair)
xs_trig = {"prior": xs_t["PP.H4_p12_1e6"] < 0.5, "susie": xs_t.susie_class == "discordant", "opposing_arm": xs_opp,
           "assay": xs_t.decode_status == "ASSAY_DISCORDANT", "steiger": xs_t.steiger_result != "FORWARD_CONFIRMED"}
xs_weak = np.logical_or.reduce([v.values for v in xs_trig.values()])
xs_core = ((xs_t.mr_pval < 0.05 / xs_fam) & (xs_t["PP.H4_p12_1e6"] >= 0.80) & (xs_t.susie_class == "supportive") & (xs_t.pp_h4 >= 0.80) & (xs_t.pp_h4 >= 5 * xs_t.pp_h3)).values
assert (xs_core == xs_t.core_complete.values).all()
xs_tier = np.where(xs_weak, "Flagged", np.where(xs_core & (xs_t.n_evaluable >= 1).values, "Convergent", "Supported"))
assert (xs_tier == xs_t.evidence_tier.values).all(), "stated grading rules do not reproduce the tiers"
# evaluable arms: the opposite-sign status coincides with NOT_REPLICATED exactly
for xs_a in ["decode", "aric"]:
    xs_e = xs_t[xs_t[xs_a + "_evaluable_pair"]]
    assert ((xs_e[xs_a + "_status"] == "NOT_REPLICATED") == (xs_e[xs_a + "_direction_concordant"] == False)).all()  # noqa: E712
xs_ne_ = xs_t.susie_na_reason.value_counts().to_dict()
xs_cc = {"core": int(xs_core.sum()), "convergent": int((xs_core & (xs_tier == "Convergent")).sum()),
         "supported_no_arm": int((xs_core & (xs_tier == "Supported") & (xs_t.n_evaluable == 0).values).sum()), "flagged": int((xs_core & (xs_tier == "Flagged")).sum())}
assert xs_cc == {"core": 16, "convergent": 3, "supported_no_arm": 11, "flagged": 2}
assert sorted(map(tuple, xs_t[xs_core & (xs_tier == "Flagged")][["gene", "outcome"]].values)) == [("APOE", "T2D"), ("CCND2", "T2D")]
add("L409", "SM: grading", "Grading rules (re-derived from the stated rules; equal to the tiered table for 81 of 81): weakness threshold (PP.H4 at p12 = 1e-6 below); the susie_classify() fallback threshold for discordant (PP.H3 above; used only when no precomputed label exists); weakness triggers; SuSiE not evaluable (non-converged / MHC); pairs meeting all four internal criteria and their tiers",
    {"prior_weak_below": 0.5, "susie_discordant_above": 0.5, "triggers": {k: int(v.sum()) for k, v in xs_trig.items()},
     "susie_not_evaluable": {"nonconverged": int(xs_ne_.get("susie_nonconvergence", 0)), "mhc": int(xs_ne_.get("MHC_excluded", 0))}, "core_complete": xs_cc},
    "weakness: conservative-prior posterior below 0.5 (27); SuSiE discordant, distinct-variant posterior above 0.50 (1); opposing evaluable arm (7); assays disagree (1); Steiger not forward (1); SuSiE not evaluable: 2 non-converged, 1 MHC; all four internal criteria: 16 (3 Convergent, 11 Supported without an arm, 2 Flagged)",
    rel(xs_TIER) + " + theme_paper3.R assign_evidence_tiers + 12_rebuild_pair_evidence.R susie_classify", "re-derivation asserted", "[D]",
    note="susie_classify() returns a precomputed SuSiE-pipeline label when one exists, so the supportive/inconclusive/discordant definition is not fully local (SM prints [DEFINITION TO CONFIRM]).")
assert {k: int(v.sum()) for k, v in xs_trig.items()} == {"prior": 27, "susie": 1, "opposing_arm": 7, "assay": 1, "steiger": 1}
assert (int(xs_ne_.get("susie_nonconvergence", 0)), int(xs_ne_.get("MHC_excluded", 0))) == (2, 1)

# ---------------------------------------------------------------------------------------------------- L410 grading chronology
xs_cmb = xs_txt(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "cascade", "06_combine_coloc.R"))
xs_has(xs_cmb, "2026-07-04 D1")
xs_has(xs_thm, "D2 2026-07-04: Convergent requires study-wide Bonferroni", "EVALUABILITY GATE ADDED 2026-08-10", "SYMMETRIC ACROSS THE TWO PLATFORMS 2026-08-11",
       "CORE / CONVERGENT SPLIT 2026-08-11")
xs_has(xs_rb, "2026-07-02 (editor #3): SuSiE could not be COMPUTED (NA)", "NOT AN ESTIMATE OF ANYTHING (editor ruling 2026-08-08)")
xs_dates = sorted(set(xs_re.findall(r"20\d\d-\d\d-\d\d", xs_thm[xs_thm.index("assign_evidence_tiers <- function("):xs_thm.index("# ── Print summary ──")])))
assert xs_dates[-1] == "2026-08-11", xs_dates   # no dated rule change after 11 August 2026 inside the grading function
xs_drv = xs_txt(os.path.join(xs_LOGS, "driver.log"))
xs_has(xs_drv, "== 14_assign_tiers done 2026-09-17", "== 14_assign_tiers done 2026-09-22")
xs_gd = ["2026-07-02", "2026-07-04", "2026-08-08", "2026-08-10", "2026-08-11", "2026-09-17", "2026-09-22"]   # each pinned by an xs_has() above


def xs_dmy(d):
    import datetime as xs_dt
    x = xs_dt.date.fromisoformat(d)
    return "%d %s %d" % (x.day, x.strftime("%B"), x.year)


add("L410", "SM: grading", "Dates each grading rule was introduced (grading code records) and dates the present analysis was graded: SuSiE not computed = not evaluable; colocalization rule + Convergent needs Bonferroni; non-converged SuSiE = not evaluable; opposing arm counts only if evaluable; symmetric arms + Convergent needs an evaluable arm; graded",
    xs_gd, "; ".join(xs_dmy(d) for d in xs_gd[:5]) + "; graded %s and %s" % (xs_dmy(xs_gd[5]).split()[0], xs_dmy(xs_gd[6])),
    "deposit v4: 12_rebuild_pair_evidence.R, 06_combine_coloc.R, R/theme_paper3.R comments; run_logs/RERUN_20260915/driver.log", "dated comments; 14_assign_tiers done", "[D]",
    note="The latest dated comment inside assign_evidence_tiers() is 2026-08-11 (asserted). Dates are of the code records, written during development on earlier versions.")

# ---------------------------------------------------------------------------------------------------- L411 PheWAS flags
xs_AM = os.path.join(V3, "v3.2", "_engine_root", "results", "repro_phase1_hg38", "_fig_results", "phase4_annotation", "phase4_annotation_master.csv")
xs_am = pd.read_csv(xs_AM)
assert len(xs_am) == 62 and set(xs_am.gene) == set(xs_t.gene)
xs_cs = xs_am.phewas_consensus.value_counts().to_dict()
xs_ot = xs_am.phewas_class_ot.value_counts().to_dict()
xs_otf = xs_ot.get("pleiotropy_flagged", 0) + xs_ot.get("safety_flagged", 0)
xs_otq = len(xs_am) - xs_ot.get("not_queryable", 0)
assert xs_am[xs_am.phewas_class_ot == "not_queryable"].gene.tolist() == ["AGER"]
xs_pw = xs_txt(os.path.join(xs_LOGS, "anno_08_phewas.log"))
xs_has(xs_pw, "Genes carrying at least one flag: %d/%d" % (xs_otf, xs_otq), "[pre-registered stop-loss: >= 90%]", "in95cs=", "passing(PIP>=0.5)")
xs_at = xs_txt(os.path.join(xs_LOGS, "anno_08c_retry.log"))
xs_na = int(xs_re.search(r"derived (\d+) safety-check targets \(high_risk OR not_queryable in OT/IEU", xs_at).group(1))
xs_pm = xs_txt(os.path.join(V3, "V4", "source_data", "scripts", "prioritization_matrix.R"))
xs_m = xs_re.search(r"editorial ruling 2026-08-16.*?fired -- (\d+)/(\d+) genes flagged", xs_pm, xs_re.S)
assert xs_m and xs_m.group(1) == xs_m.group(2)
add("L411", "SM: PheWAS", "PheWAS flags: genes; any-source consensus flagged (two or more sources / one); Open Targets alone flagged of queryable, clean, not queryable; stop-loss threshold (%); Open Targets credible-set level (%) and PIP; GWAS Atlas targets; genes flagged when the rule was first applied (earlier version)",
    {"genes": len(xs_am), "any_source": int(xs_cs.get("flag_multi_source", 0) + xs_cs.get("flag_single_source", 0)), "multi": int(xs_cs.get("flag_multi_source", 0)),
     "single": int(xs_cs.get("flag_single_source", 0)), "ot_flagged": xs_otf, "ot_queryable": xs_otq, "ot_clean": int(xs_ot.get("clean", 0)), "ot_not_queryable": 1,
     "stoploss_pct": 90, "credible_set_pct": 95, "pip": 0.5, "atlas_targets": xs_na, "earlier_flagged": int(xs_m.group(1))},
    "any source %d of %d (%d multi, %d single); Open Targets %d of %d queryable (%d clean; AGER not queryable); stop-loss 90%%; 95%% credible sets, PIP 0.5; GWAS Atlas %d targets; earlier version %s of %s"
    % (xs_cs.get("flag_multi_source", 0) + xs_cs.get("flag_single_source", 0), len(xs_am), xs_cs.get("flag_multi_source", 0), xs_cs.get("flag_single_source", 0),
       xs_otf, xs_otq, xs_ot.get("clean", 0), xs_na, xs_m.group(1), xs_m.group(2)),
    rel(xs_AM) + " + run_logs anno_08_phewas.log, anno_08c_retry.log + V4/source_data/scripts/prioritization_matrix.R", "phewas_consensus, phewas_class_ot", "[D]")
assert (xs_cs.get("flag_multi_source", 0), xs_cs.get("flag_single_source", 0), xs_otf, xs_otq, xs_ot.get("clean", 0)) == (31, 31, 27, 61, 34)

# ---------------------------------------------------------------------------------------------------- L412 PheWAS / development-stage dates
xs_mm = xs_txt(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "phase4", "09_merge_master.py"))
xs_has(xs_mm, "THE PRE-REGISTERED VOCABULARY (phewas_pip_threshold_PREREGISTERED.md", "editor ruling 2026-08-08")
for xs_s in ["anno_03_opentargets", "anno_04_chembl", "anno_05_dgidb", "anno_08_phewas"]:
    assert xs_re.search(r"== %s done 2026-09-17 " % xs_s, xs_drv), xs_s
assert xs_re.search(r"== anno_08b_retry done 2026-09-17 ", xs_drv) and xs_re.search(r"== anno_08c_r2 done 2026-09-17 ", xs_drv)
add("L412", "SM: PheWAS / development stage", "Dates: the PheWAS pre-specification document cited in a recorded ruling; stop-loss first applied (earlier version); annotation queries of the present analysis (Open Targets, ChEMBL, DGIdb, PheWAS)",
    ["2026-08-08", "2026-08-16", "2026-09-17"], "cited %s; first applied %s; queried %s" % (xs_dmy("2026-08-08"), xs_dmy("2026-08-16"), xs_dmy("2026-09-17")),
    "deposit v4 09_merge_master.py comment; V4/source_data/scripts/prioritization_matrix.R comment; run_logs/RERUN_20260915/driver.log", "dated comments; anno_* done", "[D]",
    note="Database releases live on 17 September 2026 are not recorded locally (the raw caches are on the source machine); the SM prints [RELEASE TO CONFIRM].")

# ---------------------------------------------------------------------------------------------------- L413 development-stage phase codes
xs_has(xs_mm, "if chembl_direct_phase is not None and chembl_direct_phase >= 4:\n            drug_tier_direct = 1",
       "elif chembl_direct_phase is not None and chembl_direct_phase >= 1:\n            drug_tier_direct = 2", "# drug_tier_direct: ChEMBL SINGLE PROTEIN only")
add("L413", "SM: development stage", "Development-stage rule: ChEMBL maximum phase for an approved agent (direct single-protein mechanisms)",
    {"approved_phase": 4}, "approved: phase 4",
    rel(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "phase4", "09_merge_master.py")), "drug_tier_direct == 1", "[D]")
add("L414", "SM: development stage", "Development-stage rule: ChEMBL phase range for clinical development (direct single-protein mechanisms; below the approved phase)",
    {"clinical_min": 1, "clinical_max": 3}, "clinical: phase 1 to 3",
    rel(os.path.join(xs_DEP, "scripts", "paper3.1_ckm", "phase4", "09_merge_master.py")), "drug_tier_direct == 2", "[D]",
    note="Split from L413 so that a swap of the approved and clinical phase codes trips the gate (a row holding 1, 3 and 4 would accept either).")

# ---------------------------------------------------------------------------------------------------- L415 decoy provenance (verifier C1V-01, C1V-10)
# The rates splice two runs of 05_run_null.R: the 2026-08-20 run (pre-fix engine) for six non-CKD outcomes and eGFR, and a
# CKD stratum rebuilt from the corrected captures (V3/v3.2/scripts/build_v4_local_null.R). Asserted here: the dates in the
# deposit report, the splice counts in the rebuild log, and that every KEPT target's true-pairing PP.H4 equals the corrected
# coloc table at the three decimals the null file stores (non-CKD and eGFR).
xs_rep = xs_txt(os.path.join(V3, "paper3_deposit_request_v4", "REPORT.md"))
xs_has(xs_rep, "| A17 | `results/repro_phase1_hg38/coloc_null/null_draws.csv` | 131200 | 7466791 | 2026-08-20 15:10:10",
       "| A18 | `results/repro_phase1_hg38/coloc_null/null_draws_egfr.csv` | 32100 | 1778496 | 2026-08-20 15:20:44",
       "Nine of the 33 carry bytes that predate the 2026-09-15 engine fix",
       "A18/A20 hold eGFR rows only, drawn from `05_coloc_egfr.R`'s own captures, whose reader was not the defective one")
xs_bl = xs_txt(os.path.join(V3, "V4_local", "source_data", "coloc_null", "BUILD_LOG.txt"))
xs_m = xs_re.search(r"CKD: (\d+) captures, (\d+) MHC-excluded, (\d+) targets, (\d+) draws", xs_bl)
xs_ckdcap, xs_ckdmhc, xs_ckdtg, xs_ckddr = [int(x) for x in xs_m.groups()]
xs_kept = int(xs_re.search(r"kept verbatim: (\d+) observed \+ (\d+) draw lines", xs_bl).group(1))
assert (xs_ckdcap, xs_ckdtg, xs_ckddr) == (168, 162, 16200) and xs_kept + xs_ckdtg == int(xs_ne.targets_tested) == 1332
xs_has(xs_bl, "spliced: tested 1332 | observed passes 64 | draws 133200 | null passes 501")
xs_co = pd.read_csv(os.path.join(V3, "V4", "source_data", "tables_from_rds", "phase3_coloc__coloc_results_final.csv"))
xs_NO = os.path.join(V3, "v3.2", "_engine_root", "results", "v4_local", "coloc_null", "null_observed.csv")             # read by S1R1_decoy_rates.R
xs_NE = os.path.join(V3, "v3.2", "_engine_root", "results", "repro_phase1_hg38", "coloc_null", "null_observed_egfr.csv")
xs_rates_txt = xs_txt(os.path.join(S1, "S1R1_03_decoy_method", "S1R1_decoy_rates.R"))
xs_has(xs_rates_txt, '"v3.2/_engine_root/results/v4_local/coloc_null/null_observed.csv"',
       '"v3.2/_engine_root/results/repro_phase1_hg38/coloc_null/null_observed_egfr.csv"')
xs_eq = {}
for xs_lab, xs_pth, xs_keep in [("non_ckd", xs_NO, lambda d: d.outcome != "Wuttke2019_CKD_EA"), ("egfr", xs_NE, lambda d: d.outcome == d.outcome)]:
    xs_o = pd.read_csv(xs_pth)
    xs_o = xs_o[xs_keep(xs_o)].merge(xs_co[["gene", "outcome_id", "PP.H4"]], left_on=["gene", "outcome"], right_on=["gene", "outcome_id"], how="left")
    assert xs_o["PP.H4"].notna().all() and xs_o.H4.notna().all(), xs_lab
    xs_eq[xs_lab] = (int((np.abs(xs_o.H4 - xs_o["PP.H4"].round(3)) < 1e-9).sum()), len(xs_o))
assert xs_eq == {"non_ckd": (1170, 1170), "egfr": (321, 321)} and xs_kept == 1170 and int(xs_eg.targets_tested) == 321, xs_eq
xs_bseed = int(xs_re.search(r"set\.seed\((\d+)\); NB <- \d+", xs_rates_txt).group(1))
assert xs_bseed == 20261002
xs_dm = xs_txt(os.path.join(S1, "S1R1_03_decoy_method", "S1R1_DECOY_METHOD.md"))
xs_m = xs_re.search(r"In (\d+) of (\d+) CKD captures the vector is sorted lexicographically", xs_dm)
xs_srt = (int(xs_m.group(1)), int(xs_m.group(2)))
assert xs_srt == (166, 168) and xs_srt[1] == xs_ckdcap
add("L415", "SM: decoy-locus stress test", "Decoy provenance: date of the run that supplied the six non-CKD non-eGFR outcomes and eGFR; date of the outcome-reader correction; CKD stratum rebuilt from corrected captures (captures, targets); kept targets whose true-pairing PP.H4 equals the corrected analysis at 3 dp (non-CKD; eGFR); bootstrap seed of the intervals; CKD captures in alphabetical order of variant identifier",
    {"run": "2026-08-20", "fix": "2026-09-15", "ckd_captures": xs_ckdcap, "ckd_targets": xs_ckdtg, "kept_non_ckd_equal": list(xs_eq["non_ckd"]),
     "egfr_equal": list(xs_eq["egfr"]), "bootstrap_seed": xs_bseed, "ckd_sorted": list(xs_srt)},
    "run of 20 August 2026; reader corrected 15 September 2026; CKD rebuilt from %d captures (%d targets); kept targets reproducing the corrected PP.H4 at 3 dp: %s of %s non-CKD, %d of %d eGFR; bootstrap seed %d; CKD captures alphabetical %d of %d"
    % (xs_ckdcap, xs_ckdtg, "{:,}".format(xs_eq["non_ckd"][0]), "{:,}".format(xs_eq["non_ckd"][1]), xs_eq["egfr"][0], xs_eq["egfr"][1], xs_bseed, xs_srt[0], xs_srt[1]),
    "V3/paper3_deposit_request_v4/REPORT.md (A17/A18 rows) + V3/V4_local/source_data/coloc_null/BUILD_LOG.txt + null_observed{,_egfr}.csv vs V4 coloc table + S1R1_decoy_rates.R + S1R1_DECOY_METHOD.md",
    "mtimes; splice counts; H4 vs round(PP.H4, 3); set.seed; sort count", "[D][S]",
    note="166 of 168 is R1 skeptic 1's count (S1R1_DECOY_METHOD.md); re-derived 2026-10-06 by the C1 builder in R over V3/paper3_ckd_rerun_v3/capture_ckd (scratchpad c1/sortcheck.R): 166 in byte and locale order; exceptions EGFL7, TRAF2. The eGFR equality (321 of 321) is stronger than the deposit report's 'unaffected for the 16 coloc-supported eGFR pairs only'.")

# ---------------------------------------------------------------------------------------------------- L416 eGFR scale factor of the East Asian screen (verifier C1V-04)
xs_EAS = os.path.join(CKM, "paper 3.eas", "analysis")
xs_rf = xs_txt(os.path.join(xs_EAS, "src", "phase2_egfr_refined.py"))
xs_fac = float(xs_re.search(r"^SCALE=([0-9.]+)\s*$", xs_rf, xs_re.M).group(1))
xs_has(xs_rf, "eur_eq=eur_b*SCALE; eur_eq_se=eur_se*SCALE", "hetz=(eas_b-eur_eq)/math.sqrt(eas_se**2+eur_eq_se**2)",
       'r["eas_q"]<0.05 and r["eur_p"]>0.05 and common', 'and r["het_p"]<0.05 and not r["pleiotropic_locus"]')
xs_sj = json.load(io.open(os.path.join(V3, "V4_local", "staging", "d3-metab-ukbfree", "four_quadrant", "EGFR_SCALE_rederivation.json"), encoding="utf-8"))
xs_has(xs_txt(os.path.join(xs_EAS, "src", "scale_calibrate.py")), "if ep<5e-8 and abs(Ub)>1e-6:", "if palin(ea,oa): continue", "w=1.0/ese**2; num+=w*Ub*b; den+=w*Ub**2")
assert xs_fac == 6.16 and xs_sj["reproduces_6_16"] is True and abs(xs_sj["wls_slope"] - xs_fac) < 0.005 and xs_sj["n_loci"] == 2271
xs_rt = pd.read_csv(os.path.join(xs_EAS, "results", "phase2_egfr_refined.tsv"), sep="\t")
xs_cand = xs_rt[xs_rt.eas_q < 0.05].copy()
xs_cl = pd.read_csv(os.path.join(xs_EAS, "results", "phase2_coloc_ledger.tsv"), sep="\t", comment="#")
assert sorted(xs_cand.protein) == sorted(xs_cl[xs_cl.outcome == "eGFR"].protein) == ["AGER", "APOE", "FAP", "PBLD", "TF"]
xs_cand["common"] = xs_cand.eas_af.between(0.05, 0.95, inclusive="neither") & xs_cand.eur_af.between(0.05, 0.95, inclusive="neither")
xs_cand["fail_other"] = ~(xs_cand.eur_p > 0.05) | ~xs_cand.common | xs_cand.pleiotropic_locus.astype(bool)
xs_cand["fail_het"] = ~(xs_cand.het_p < 0.05)
assert not xs_cand.strict_EAS_specific.astype(bool).any()
xs_only_het = sorted(xs_cand[xs_cand.fail_het & ~xs_cand.fail_other].protein)
assert xs_only_het == ["TF"] and int(xs_cand.fail_other.sum()) == 4
add("L416", "SM: eGFR scale", "Earlier East Asian screen, eGFR arm: factor multiplying European ln(eGFR) estimates before the heterogeneity filter; variants behind it (genome-wide significant in both ancestries, non-palindromic; through-origin inverse-variance-weighted slope); eGFR candidates; candidates failing only the heterogeneity filter (TF); failing a scale-free filter",
    {"factor": xs_fac, "wls_slope": round(float(xs_sj["wls_slope"]), 4), "variants": int(xs_sj["n_loci"]), "egfr_candidates": len(xs_cand),
     "only_het": len(xs_only_het), "scale_free_fail": int(xs_cand.fail_other.sum())},
    "factor %.2f (slope %.3f over %s variants); %d eGFR candidates: TF failed only the heterogeneity filter, %d failed a scale-free filter"
    % (xs_fac, xs_sj["wls_slope"], "{:,}".format(xs_sj["n_loci"]), len(xs_cand), int(xs_cand.fail_other.sum())),
    "paper 3.eas/analysis/src/phase2_egfr_refined.py (SCALE) + results/phase2_egfr_refined.tsv + phase2_coloc_ledger.tsv + V3/V4_local/staging/d3-metab-ukbfree/four_quadrant/EGFR_SCALE_rederivation.json (scale_calibrate.py re-run)",
    "SCALE; eas_q < 0.05 rows; filter columns", "[D]",
    note="Used only for the screen's eGFR heterogeneity filter (and the retired metabolite four-quadrant analysis, not in this paper). All five eGFR candidates entered colocalization whatever the filters gave.")

# ---------------------------------------------------------------------------------------------------- L417 BUN scale (verifier C1V-11)
# Wuttke 2019 (old:32, the BUN file used) states only that studies fitted "log(eGFR) and BUN"; Stanzick 2021 (old:34) states that
# the CKDGen BUN results are natural-log BUN. The per-allele SEs implied by the file decide between ln(BUN) (SD about 0.3) and mg/dl
# (SD about 4-5), whose SEs differ about 15-fold at this sample size.
xs_KC = os.path.join(V3, "paper3_ckd_rerun_v3", "out", "D1", "kidney_cystatin_bun_concordance_E_current_20260922.csv")
xs_kc = pd.read_csv(xs_KC).drop_duplicates("gene").merge(xs_ins.drop_duplicates("gene")[["gene", "beta"]], on="gene")
xs_kc = xs_kc[xs_kc.bun_b.notna()]
xs_bse = (np.abs(xs_kc.bun_b * xs_kc.beta) / xs_norm.isf(xs_kc.bun_p / 2)).values
assert len(xs_kc) == 16 and np.isfinite(xs_bse).all() and float(np.median(xs_bse)) < 0.003
xs_has(xs_txt(os.path.join(xs_DEP, "scripts", "_shared", "renal_triad", "01_renal_triad_mr.R")), 'bun <- wald_vs(dt, file.path(KM, "CKDGen_BUN_EA_Wuttke2019.txt.gz")')
add("L417", "SM: eGFR scale", "Blood urea nitrogen file (CKDGen, Wuttke 2019): median per-allele standard error implied by the Wald ratio and its p, over the 16 colocalization-supported kidney proteins with a BUN estimate",
    {"median_se": float(np.median(xs_bse)), "proteins": len(xs_kc)}, "median per-allele SE %.4f over %d proteins" % (np.median(xs_bse), len(xs_kc)),
    rel(xs_KC) + " + UKB-PPP instrument beta", "|b_wald * beta_x| / z(p)", "[D]",
    note="Supports ln(BUN): untransformed mg/dl would give about 0.015 at this sample size. Same check as figures/scripts/P31R_Fig3.R (median < 0.003).")

# ---------------------------------------------------------------------------------------------------- L418 mode of action vs genetic direction, all proteins with an approved or clinical agent (verifier C1V-07)
xs_scripts = os.path.join(V3, "R2_rewrite", "scripts")
import sys as xs_sys
if xs_scripts not in xs_sys.path:
    xs_sys.path.insert(0, xs_scripts)
import P31R_af1_devstage as xs_dev  # noqa: E402  (the one implementation behind Additional file 1: Table S4)
xs_dt, _xs_ds = xs_dev.table()
xs_ag = xs_dt[xs_dt.stage_n <= 2]
xs_vc = {k: sorted(xs_ag[xs_ag.verdict.str.startswith(k)].gene) for k in ("concordant", "opposite", "unknown")}
xs_rec = sorted(xs_ag[xs_ag.verdict.str.contains("membrane receptor")].gene)
assert len(xs_ag) == 12 and sum(map(len, xs_vc.values())) == 12 and xs_rec == xs_vc["unknown"]
assert xs_vc == {"concordant": ["PCSK9", "TGFB1", "TGFB2"], "opposite": ["ADM", "DPEP1", "FN1", "PLAU"], "unknown": ["AGER", "CD27", "CD40", "FLT1", "IL6R"]}, xs_vc
assert xs_ag.set_index("gene").loc["ADM", "stage"] == "in clinical development" and set(_xs_ds["approved_opposite"]) == {"DPEP1", "FN1", "PLAU"}
xs_dvsrc = "V3/R2_rewrite/scripts/P31R_af1_devstage.py table() (the Table S4 builder)"
# one row per count, with no gene names in the display: digits inside gene symbols (TGFB2, IL6R) would otherwise let a
# mutated count match the row (mutation probe, triage round)
add("L418", "SM: development stage", "Proteins with an approved or clinical-stage direct agent (Table S4)", {"proteins": len(xs_ag)},
    "%d proteins with an approved or clinical agent" % len(xs_ag), xs_dvsrc, "stage_n <= 2", "[D]",
    note="Mode of action vs genetic direction: concordant %s; opposite %s; unknown %s (all single-pass type I membrane receptors). ADM is the fourth 'opposite' protein (enibarcimab, an ADM stabiliser, in clinical development)."
    % (", ".join(xs_vc["concordant"]), ", ".join(xs_vc["opposite"]), ", ".join(xs_vc["unknown"])))
add("L420", "SM: development stage", "Of the proteins with an approved or clinical agent: agent mode of action concordant with the genetically implied direction", {"concordant": len(xs_vc["concordant"])},
    "concordant: %d proteins" % len(xs_vc["concordant"]), xs_dvsrc, "verdict == concordant", "[D]")
add("L421", "SM: development stage", "Of the proteins with an approved or clinical agent: agent mode of action opposite to the genetically implied direction", {"opposite": len(xs_vc["opposite"])},
    "opposite: %d proteins" % len(xs_vc["opposite"]), xs_dvsrc, "verdict == opposite", "[D]")
add("L422", "SM: development stage", "Of the proteins with an approved or clinical agent: comparison unknown (agents act in both directions or the protein is a single-pass type I membrane receptor); all unknown ones are receptors",
    {"unknown": len(xs_vc["unknown"]), "unknown_receptor": len(xs_rec)}, "unknown: %d proteins (all membrane receptors)" % len(xs_vc["unknown"]), xs_dvsrc, "verdict startswith unknown", "[D]")

# ---------------------------------------------------------------------------------------------------- L419 SuSiE grid settings of the pair-level analysis (for the Methods purity/L edit)
xs_sg2 = pd.read_csv(os.path.join(V3, "V4", "source_data", "tables_from_rds", "phase3_susie__susie_grid_results.csv"))
xs_gs = xs_sg2.groupby("setting")[["L", "min_abs_corr"]].agg(lambda x: sorted(set(x)))
xs_gsd = {s: (int(xs_gs.loc[s, "L"][0]), float(xs_gs.loc[s, "min_abs_corr"][0])) for s in xs_gs.index}
assert all(len(xs_gs.loc[s, c]) == 1 for s in xs_gs.index for c in ("L", "min_abs_corr"))
assert xs_gsd == {"strict": (10, 0.5), "relaxed": (10, 0.2), "exploratory": (10, 0.1), "robustness": (5, 0.2)}, xs_gsd
add("L419", "Methods / SM: grading", "Pair-level SuSiE grid (deposited results table): maximum signals L and credible-set purity (min_abs_corr) per setting: strict, relaxed, exploratory, robustness",
    {s: {"L": v[0], "purity": v[1]} for s, v in xs_gsd.items()},
    "strict L 10, purity 0.5; relaxed L 10, purity 0.2; exploratory L 10, purity 0.1; robustness L 5, purity 0.2",
    rel(os.path.join(V3, "V4", "source_data", "tables_from_rds", "phase3_susie__susie_grid_results.csv")), "setting, L, min_abs_corr", "[D]",
    note="Which setting the 81 pairs' susie_class comes from is not local (12_rebuild_pair_evidence.R susie_classify returns a precomputed label); the candidate SuSiE's L is not local either (scripts 28-31 not deposited).")

# ---------------------------------------------------------------------------------------------------- independent checks (Stage 3b C1 verifier; different agent, different language and code)
# V3 scratchpad c1v/rederive.R and c1v/grade.R (data.table; BH via p.adjust over the padded family; loci by an all-pairs union-find,
# not the builder's sorted single pass); re-run by the C1 builder on 2026-10-06 and its printed output read against each row.
xs_C1V = "Stage 3b C1 verifier, scratchpad c1v/%s (R, data.table; separately written; re-run 2026-10-06): %s"
EXT_CHECKS["L400"] = xs_C1V % ("rederive.R", "MATCH: 98 IVW (n_snps 2-3), 14,423 Wald, 1,063 untestable")
EXT_CHECKS["L401"] = xs_C1V % ("rederive.R", "MATCH: BH padded to 15,584 admits 292, largest p 9.362e-4; stored fdr_sig 304; 12 extra Wald (p 9.447e-4 to 1.049e-3), 0 in the 81, 0 in the 248 pool; DPEP1-AF p 0.0136")
EXT_CHECKS["L404"] = xs_C1V % ("rederive.R", "MATCH: 177 failing, 17 MHC by gene overlap 25-34 Mb, 160 testable; loci by union-find 9/14 vs 44/120 (p 0.0798), 134 loci")
EXT_CHECKS["L408"] = xs_C1V % ("grade.R", "MATCH: 11 pairs lost at 5e-8; p 5.67e-8 (BRAP) to 8.22e-7 (CCND2); F 24.31 to 29.47")
EXT_CHECKS["L409"] = xs_C1V % ("grade.R", "MATCH: triggers 27/1/7/1/1; SuSiE not evaluable 2 non-converged + 1 MHC; 16 core = 3 Convergent, 11 Supported (no arm), 2 Flagged (CCND2-T2D, APOE-T2D)")

# ---------------------------------------------------------------------------------------------------- provenance note
for xs_r in ROWS:
    if xs_r["id"] in ("L400", "L401", "L403", "L404", "L406", "L408", "L409", "L411"):
        xs_r["note"] = (xs_r["note"] + " " if xs_r["note"] else "") +             "Cross-checked 2026-10-06 in base R by the C1 builder (scratchpad c1/xcheck.R; same agent, different language): all values matched. Not an independent check."
