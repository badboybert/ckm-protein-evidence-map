# 11_compare_seen.R -- COMPARE_SEEN.tsv. Runs only after results.tsv and branches.tsv exist (rule 10: compute first, compare second).
# Only now are the full pre-declaration / amendment and the seen outputs opened. No computed value is changed here.
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "00_common.R"))
say("==== W6b 11_compare_seen ====")
stopifnot(file.exists(file.path(OUT, "results.tsv")), file.exists(file.path(OUT, "branches.tsv")))
say("results.tsv md5 %s and branches.tsv md5 %s exist before the full text is opened", md5(file.path(OUT, "results.tsv")), md5(file.path(OUT, "branches.tsv")))
assert_md5(GOV$prespec_full[1], GOV$prespec_full[2], "GOV:prespec_full (opened now)")
assert_md5(GOV$amend_full[1], GOV$amend_full[2], "GOV:amend_full (opened now)")
R6b <- fread(file.path(OUT, "results.tsv"), colClasses = "character", na.strings = NULL)
v <- function(id) { x <- R6b[quantity_id == id]$value; if (length(x) != 1) stop("no result ", id); x }
C <- list()
cmp <- function(qid, seen, seen_label, s1, s1_label, explanation = "") {
  sn <- suppressWarnings(as.numeric(seen)); s1n <- suppressWarnings(as.numeric(s1))
  d <- if (length(sn) == 1 && length(s1n) == 1 && !is.na(sn) && !is.na(s1n)) abs(sn - s1n) else NA_real_
  C[[length(C) + 1]] <<- data.table(quantity_id = qid, seen_value = as.character(seen), seen_label = seen_label, stage1_value = as.character(s1), stage1_label = s1_label,
                                     abs_diff = if (is.na(d)) "" else formatC(d, digits = 4, format = "g"), explanation = explanation)
}
ST <- file.path(ROOT, "V3/V4_local/staging/l3-finemap/UMOD_eGFR_control")
SEEN_READ <- c()
rd <- function(p, id = NULL) { if (!is.null(id)) res(id, allow_seen = TRUE); SEEN_READ <<- c(SEEN_READ, sprintf("%s md5 %s", p, md5(p))); fread(p) }

# ---- (1) staging two-signal module (cond PP.H4 at A and B; conditional z), 503 = p1np, 489 = p2np ----
for (pan in c("503", "489")) {
  vid <- if (pan == "503") "LOC_UMOD_eGFR_control_two_signal_p1np" else "LOC_UMOD_eGFR_control_two_signal_p2np"
  ts <- rd(res_entry(vid)$path, vid)
  r <- ts[protein == "P__UMOD" & other == "O__eGFR"]
  stopifnot(nrow(r) == 1, r$sig_A == "rs12922822", r$sig_B == "rs77924615")
  lab <- sprintf("staging %s two_signal.tsv (engine l3fm__conditional.py, Python coloc.abf engine)", if (pan == "503") "conditional_p1np" else "conditional_p2np")
  cmp(sprintf("S6.umod.H4.A.%s", pan), r$cond_H4_A, lab, v(sprintf("S6.umod.H4.A.%s", pan)), "Stage 1 fresh R coloc.abf on z-space-conditioned data",
      "same rule (both traits conditioned on B), different coloc implementation (R coloc 5.2.3 vs the vendored D1 engine)")
  cmp(sprintf("S6.umod.H4.B.%s", pan), r$cond_H4_B, lab, v(sprintf("S6.umod.H4.B.%s", pan)), "Stage 1 fresh R coloc.abf on z-space-conditioned data", "as above")
  cmp(sprintf("S6.umod.z.A.%s", pan), r$MA, paste(lab, "MA = z_eGFR(A|B)"), v(sprintf("S6.umod.z.A.%s", pan)), "Stage 1 two-SNP z-space conditioning", "same formula; LD from PLINK vs dosage correlation")
  cmp(sprintf("S6.umod.z.B.%s", pan), r$MB, paste(lab, "MB = z_eGFR(B|A)"), v(sprintf("S6.umod.z.B.%s", pan)), "Stage 1 two-SNP z-space conditioning", "as above")
}

# ---- (2) staging SuSiE fits for UMOD x eGFR: s per trait and set, coloc.susie max PP.H4 ----
SETDIR <- c(susie = "susie", p1fr = "susie_p1fr", p1np = "susie_p1np", p2np = "susie_p2np"); PAN <- c(susie = "503", p1fr = "503", p1np = "503", p2np = "489")
for (sname in names(SETDIR)) {
  tr <- rd(file.path(ST, SETDIR[[sname]], "susie_traits.tsv")); ps <- rd(file.path(ST, SETDIR[[sname]], "pair_summary.tsv"))
  lab <- sprintf("staging UMOD_eGFR_control/%s (engine l3fm__susie.R on the engine harmonised table; LD = cor() of the i8 extract)", SETDIR[[sname]])
  q1 <- sprintf("S6.d.UMOD_eGFR.s.protein.%s.%s", sname, PAN[[sname]]); q2 <- sprintf("S6.d.UMOD_eGFR.s.outcome.%s.%s", sname, PAN[[sname]])
  q3 <- sprintf("S6.d.UMOD_eGFR.coloc_susie_H4.%s.%s", sname, PAN[[sname]])
  cmp(q1, tr[tid == "P__UMOD"]$s_hat, lab, v(q1), "Stage 1 local proxy fit (same SNP set, re-implemented in R; LD from PLINK)",
      "same estimator (estimate_s_rss), same SNPs and n; LD from PLINK vs dosage cor (identical to 1e-16 at the signal SNPs)")
  cmp(q2, tr[tid == "O__eGFR"]$s_hat, lab, v(q2), "Stage 1 local proxy fit", "as above")
  cmp(q3, ps[protein == "P__UMOD" & other == "O__eGFR"]$susie_max_H4, lab, v(q3), "Stage 1 local proxy fit", "same settings (L = 10, purity 0.5, p12 5e-6)")
}

# ---- (3) staging 5 x 3 signal grid (conditional_pairs.tsv signals_grid) ----
cp <- rd(file.path(ST, "conditional", "conditional_pairs.tsv"))
g <- cp[protein == "P__UMOD" & other == "O__eGFR"]$signals_grid
cells <- regmatches(g, gregexpr("P([0-9])\\([^)]*\\)xM([0-9])\\([^)]*\\)=[0-9.eE+-]+", g))[[1]]
for (cl in cells) {
  m <- regmatches(cl, regexec("P([0-9])\\([^)]*\\)xM([0-9])\\([^)]*\\)=([0-9.eE+-]+)", cl))[[1]]
  q <- sprintf("S6.umod.grid.P%sxM%s", m[2], m[3])
  cmp(q, m[4], "staging UMOD_eGFR_control/conditional/conditional_pairs.tsv signals_grid (engine; printed to 3 decimals)", v(q),
      "Stage 1 fresh R coloc.abf on conditioned data (full harmonised set, 503)", "same rule; engine value printed to 3 decimals")
}

# ---- (4) full pre-declaration: section 6 / 5 Seen values and the recorded classifications ----
FT <- readLines(GOV$prespec_full[1], warn = FALSE, encoding = "UTF-8")
grab <- function(pat) { x <- grep(pat, FT, value = TRUE); if (length(x) != 1) stop("full-text line not unique: ", pat); trimws(x) }
l_s <- grab("^  - s: protein 0\\.290")
cmp("S6.d.UMOD_eGFR.s.protein.susie|p1np|p2np", "0.290 / 0.280 / 0.285", paste("full text sec 6 Seen:", l_s),
    sprintf("%s / %s / %s", v("S6.d.UMOD_eGFR.s.protein.susie.503"), v("S6.d.UMOD_eGFR.s.protein.p1np.503"), v("S6.d.UMOD_eGFR.s.protein.p2np.489")),
    "Stage 1 per-set s (susie / p1np / p2np)", "agreement within the printed precision; weak evidence (the UMOD SNP set and estimator are the engine's)")
cmp("S6.d.UMOD_eGFR.s.outcome.susie|p1np|p2np", "0.123 / 0.119 / 0.123", paste("full text sec 6 Seen:", l_s),
    sprintf("%s / %s / %s", v("S6.d.UMOD_eGFR.s.outcome.susie.503"), v("S6.d.UMOD_eGFR.s.outcome.p1np.503"), v("S6.d.UMOD_eGFR.s.outcome.p2np.489")), "Stage 1 per-set s", "as above")
csn <- function(f) { uc <- readRDS(f)$cs; sprintf("protein %d CS (%d single-SNP); eGFR %d CS (%d single-SNP)", sum(uc$trait == "P"), sum(uc$trait == "P" & uc$size == 1),
                                                   sum(uc$trait == "O"), sum(uc$trait == "O" & uc$size == 1)) }
cmp("S6.d.UMOD_eGFR.cs.single_snp.p1np.503", "10 per trait; 8 (protein) and 6 (eGFR) single-SNP on 503", paste("full text sec 6 Seen:", grab("^  - credible sets: 10 per trait")),
    csn(file.path(WORK, "susie", "UMOD_eGFR", "p1np_503.rds")), "Stage 1 p1np/503 fit (work/susie/UMOD_eGFR/p1np_503.rds; counted here, not a results.tsv row)",
    sprintf("agreement: the pre-declaration's label '503' is the p1np set. On the full 'susie' set Stage 1 gives %s, as the staging susie/ fit does", csn(file.path(WORK, "susie", "UMOD_eGFR", "susie_503.rds"))))
others <- R6b[grepl("^S6\\.umod\\.grid\\.", quantity_id) & !quantity_id %in% c("S6.umod.grid.P2xM1", "S6.umod.grid.P1xM2", "S6.umod.grid.P3xM3")]
cmp("S6.umod.grid (named cells)", "P2xM1 1.000, P1xM2 0.906, P3xM3 0.104, other 12 cells ~0.000", paste("full text sec 6 Seen:", grab("^  - signal grid \\(503")),
    sprintf("P2xM1 %s, P1xM2 %s, P3xM3 %s; max of the other 12: %s", v("S6.umod.grid.P2xM1"), v("S6.umod.grid.P1xM2"), v("S6.umod.grid.P3xM3"),
            formatC(max(as.numeric(others$value)), digits = 3, format = "g")),
    "Stage 1 grid", "agreement; weak evidence (engine rule re-implemented on the same inputs)")
l_two <- grab("^  - two_signal: cond_H4_A 0\\.971")
cmp("S6.umod.H4.A/B (503/489)", "0.971 / 0.971 ; 1.000 / 1.000", paste("full text sec 6 Seen:", l_two),
    sprintf("%s / %s ; %s / %s", v("S6.umod.H4.A.503"), v("S6.umod.H4.A.489"), v("S6.umod.H4.B.503"), v("S6.umod.H4.B.489")), "Stage 1 fresh R coloc", "agreement")
rAB <- sub(".*r\\(A,B\\) = ([-0-9.]+).*", "\\1", R6b[quantity_id == "S6.umod.two.503.kA"]$estimand)
cmp("S6.umod.two.503.p_ign | p_aware_2nd", "p 0.21 ign / 0.31-0.32 aware (dev_sd -1.266 / -1.250; r(A,B) +0.58)", paste("full text sec 6 Seen:", l_two),
    sprintf("p_ign %s / p_aware_2nd %s (503); 489: %s / %s; r(A,B) 503 %s", v("S6.umod.two.503.p_ign"), v("S6.umod.two.503.p_aware_2nd"), v("S6.umod.two.489.p_ign"),
            v("S6.umod.two.489.p_aware_2nd"), rAB),
    "Stage 1 Rule 5 two-signal test (beta-ratio, delta method 2nd order)",
    paste("explained: different statistic. Scoping used the engine's z-space deviation dev = MA - kB*PA standardised by its sd (unit-variance conditional z; aware = with the 2 t r_AB term);",
          "Stage 1 uses Rule 5 (eGFR-per-SD-UMOD Wald ratios at A|B and B|A, second-order delta SE, with / without the LD covariance). Same side of 0.05 under both"))
cmp("S6.umod.Q.aware_2nd.p", "none ('no all-signal Q has been computed')", paste("full text sec 6 Seen:", grab("no all-signal Q has been computed")), v("S6.umod.Q.aware_2nd.p"),
    "Stage 1 all-signal GLS (first computation)", "new Stage 1 output; nothing to compare")
cmp("S6.umod.Q.boot.share", "none (no scoping all-signal Q; robustness/UMOD_eGFR_control/p1np/boot_cond.tsv = 200 draws of the two-signal module, not tabulated, not opened)",
    "full text Rule 3(ii): the existing 200-draw outputs for DPEP1 and UMOD were not tabulated", v("S6.umod.Q.boot.share"), "Stage 1 5,000-draw share (R, seed 20261007)",
    "no scoping share exists for this cell; RNG note: Stage 1 uses R Mersenne-Twister seed 20261007; the 200-draw staging file used numpy seed 20260927 for a different estimand")
l_pair <- grab("^- \\*\\*Pair test:\\*\\*")
cmp("S6.d.check_primary.or (lo-hi), p", "OR 3.33 (1.003-12.897), p 0.0329", paste("full text sec 6 Seen:", l_pair),
    sprintf("OR %s (%s-%s), p %s", v("S6.d.check_primary.or"), v("S6.d.check_primary.lo"), v("S6.d.check_primary.hi"), v("S6.d.check_primary.p")),
    "Stage 1 recomputation from the 29-log pool",
    paste("OR, lower bound and p agree; deterministic and visible in the definitions (11/16 vs 57/144); weak evidence.",
          sprintf("The upper bound DIFFERS beyond printed precision (%.3f vs 12.897), EXPLAINED (FIXLOG F1): the seen value is R fisher.test, whose uniroot runs at its default", as.numeric(v("S6.d.check_primary.hi"))),
          "tolerance (about 1.2e-4 on 1/OR); Stage 1 first printed the same 12.897, verification panel L2 flagged it, and the fixer re-solved fisher.test's own equation at tol 1e-13",
          "(agrees with scipy odds_ratio(kind='conditional') to 1e-13 relative). Fixed after the full text had been opened, on the verifier's evidence; it moves away from the seen value."))
l_d <- grab("^- \\*\\*\\(d\\), illustrative:\\*\\*")
nre <- length(strsplit(v("S6.d.reclassified"), ";")[[1]])
cmp("S6.d.or (lo-hi), p", "illustrative: UMOD-eGFR, UMOD-CKD and DPEP1-eGFR reclassified -> p 0.241, OR 2.31 (0.63-9.42); dropped -> p 0.148", paste("full text sec 6 Seen:", l_d),
    sprintf("%d reclassified (%s) -> %s/%d vs %s/%d, p %s, OR %s (%s-%s)", nre, v("S6.d.reclassified"), v("S6.d.a"), as.integer(v("S6.d.a")) + as.integer(v("S6.d.b")), v("S6.d.c"),
            as.integer(v("S6.d.c")) + as.integer(v("S6.d.d")), v("S6.d.p"), v("S6.d.or"), v("S6.d.lo"), v("S6.d.hi")),
    "Stage 1 (d): s computed for all 16 candidates",
    "different reclassification set: scoping illustrated only the three slot pairs; Stage 1 computed s (max over all seven distinct fits, both traits; four-fit sensitivity in S6.d.<cid>.s_max_4fit) for all 16, and the reclassified count is in S6.d.reclassified")
l_cls <- grab("UMOD SuSiE PP\\.H4 is inadmissible \\(s 0\\.29\\)")
cmp("CLASS:Rule3i.UMOD_SuSiE_PPH4", "inadmissible (s 0.29)", paste("full text sec 6 UMOD:", l_cls),
    sprintf("inadmissible (s_max %s over all seven fits; four-fit sensitivity %s)", formatC(as.numeric(v("S6.slot_s.UMOD_eGFR")), digits = 3, format = "f"), formatC(as.numeric(v("S6.d.UMOD_eGFR.s_max_4fit")), digits = 3, format = "f")), "Stage 1 own computed s (Rule 3(i))", "no CLASS_FLIP")
l_dp <- grab("^  - DPEP1 protein: 0\\.4034")
cmp("S6.dpep1egfr.s_protein (susie/p1fr/p1np/p2np)", "0.4034 / 0.4031 / 0.3854 / 0.3829 (DPEP1 locus sets)", paste("full text sec 5 Seen:", l_dp),
    sprintf("%s / %s / %s / %s", v("S6.dpep1egfr.susie.503.s_protein"), v("S6.dpep1egfr.p1fr.503.s_protein"), v("S6.dpep1egfr.p1np.503.s_protein"), v("S6.dpep1egfr.p2np.489.s_protein")),
    "Stage 1 DPEP1-eGFR run (own DPEP1 x eGFR SNP set)",
    "explained: different SNP sets (the locus sets join seven metabolites and two neighbour proteins; the DPEP1-eGFR set joins eGFR only); same classification (>= 0.10)")
cmp("CLASS:Rule3i.DPEP1_SuSiE (sec 5, W5's slot)", "inadmissible for DPEP1 (protein s 0.40)", paste("full text sec 5:", grab("inadmissible for DPEP1 \\(protein s 0\\.40\\)")),
    sprintf("inadmissible (DPEP1-eGFR s_max %s)", formatC(as.numeric(v("S6.slot_s.DPEP1_eGFR")), digits = 3, format = "f")), "Stage 1 own computed s (DPEP1-eGFR run)", "no CLASS_FLIP")
l_iv <- grab("arginine \\(0\\.99 / 1\\.00\\) and UMOD \\(0\\.971 / 1\\.000\\)")
cmp("Rule3iv.UMOD", "do not fire (0.971 / 1.000)", paste("full text Rule 3(iv):", l_iv),
    paste(v("S6.umod.H4.A.503"), v("S6.umod.H4.B.503"), sep = " / "), "Stage 1 branch: does not fire at A or B", "agreement")
l_d6 <- grab("the primary pair OR \\(seen: 1\\.003\\)")
cmp("D6 primary pair OR lower bound", "1.003", paste("full text sec 11:", l_d6), v("S6.d.check_primary.lo"), "Stage 1 recomputation (W6a holds the primary)", "agreement")
AM <- readLines(GOV$amend_full[1], warn = FALSE, encoding = "UTF-8")
cmp("amendment (full text)", "A1/A2 concern front-surface caveats and ARG1 T2D only", paste("full amendment:", length(AM), "lines read"), "no W6b quantity affected", "", "no comparison row needed")

# ---- (5) pipeline values (V4 table; raw_input, not the seen corpus) beside the local proxy, for the reader of the 6(d) proxy ----
pdt <- fread(res("TAB_coloc_dashed"))
cand <- fread(file.path(WORK, "candidates_windows.tsv"))
for (i in seq_len(nrow(cand))) {
  pr <- pdt[gene == cand$gene[i] & outcome == cand$outcome[i]]
  qid <- sprintf("S6.d.%s.coloc_susie_H4.susie.503", cand$cid[i])
  if (nrow(pr) == 1 && qid %in% R6b$quantity_id) {
    loc <- suppressWarnings(as.numeric(v(qid))); smx <- suppressWarnings(as.numeric(v(sprintf("S6.d.%s.s_max", cand$cid[i]))))
    big <- is.na(loc) || abs(loc - pr$pp_h4_susie_strict) > 0.05
    ex <- if (!big) "agrees within 0.05" else if (!is.na(smx) && smx >= 0.10) sprintf("differs; this pair has s_max %.3f >= 0.10, so its reference-LD SuSiE fit is the kind Rule 3(i) declares unreliable; consistent with (not proof of) the explanation", smx) else
      "unexplained: s_max < 0.10 yet the local fit differs from the pipeline fit by > 0.05 (different reference copy and SNP set; the GRCh38 window maps to two GRCh37 segments at NOTCH2)"
    cmp(qid, pr$pp_h4_susie_strict, "pipeline strict SuSiE PP.H4 (TAB_coloc_dashed; raw_input V4 table; pipeline fit on E:)", v(qid), "Stage 1 local proxy fit (susie set, 503)",
        paste("proxy fidelity check, not a reproduction target (different reference copy, 8.55M vs 22.7M variants, and SNP sets; outcomes are the same CKDGen / FinnGen R10 files):", ex))
  }
}

CS <- rbindlist(C)
fwrite(CS, file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", quote = TRUE)
for (x in SEEN_READ) say("SEEN OUTPUT READ (after results/branches): %s", x)
say("COMPARE_SEEN.tsv: %d rows; md5 %s", nrow(CS), md5(file.path(OUT, "COMPARE_SEEN.tsv")))
