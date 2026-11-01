# S5_04_assemble.R -- W5 (§5 DPEP1) step 4: bootstrap summaries, Rule 3(ii) evaluation, results.tsv and branches.tsv.
# Written from the REDACTED pre-declaration only (COMPARE_SEEN.tsv is written afterwards, by hand, from the full text).
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_05_dpep1/S5_common.R"))
say("==== S5_04_assemble start ====")
assert_governing()
P <- readRDS(file.path(WORK, "prep.rds")); A <- readRDS(file.path(WORK, "analyses.rds"))
R2 <- fread(file.path(WORK, "res_02.tsv"), colClasses = "character")
D <- fread(file.path(WORK, "boot_draws.tsv"))
stopifnot(nrow(R2) > 300L, nrow(D) == N_BOOT, setequal(D$draw, seq_len(N_BOOT)))
RES <- list()
addr <- function(quantity_id, prespec_ref, estimand, resource, panel, cov = "", se_order = "", value = NA, se = NA, ci_lo = NA,
                 ci_hi = NA, p = NA, n = NA, unit = "", notes = "")
  RES[[length(RES) + 1]] <<- data.table(quantity_id, prespec_ref, estimand, resource, panel, cov, se_order, value = as.character(value),
                                         se = as.character(se), ci_lo = as.character(ci_lo), ci_hi = as.character(ci_hi),
                                         p = as.character(p), n = as.character(n), unit, notes)
num <- function(x) trimws(formatC(x, digits = 8, format = "g"))   # trimws: formatC pads exact values (e.g. "        1")
g <- A$gls
qv <- function(tn, st, pn, cn) g[[paste(tn, st, pn, cn, sep = ".")]]

# ---------------------------------------------------------------- bootstrap summaries (Rule 3(ii), primary cells)
boot_rows <- list()
for (tn in c("cysgly", "eGFR")) {
  pf <- qv(tn, "all5", "503", "aware_2nd")$p_upper; p489 <- qv(tn, "all5", "489", "aware_2nd")$p_upper
  pb <- D[[paste0("p_Q_", tn)]]
  same <- !is.na(pb) & ((pb < 0.05) == (pf < 0.05))
  share <- sum(same) / N_BOOT
  qs <- quantile(pb, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  panels_agree <- (pf < 0.05) == (p489 < 0.05)
  robust <- panels_agree && share >= 0.90
  est <- sprintf("Rule 3(ii) on the declared primary cell (all-signal GLS, 503, aware, 2nd order), %s: share of %d bootstrap draws (503 people resampled, seed %d) whose conclusion (side of p 0.05) equals the full-panel conclusion (p %.4g, %s)",
                 tn, N_BOOT, SEED, pf, if (pf < 0.05) "p < 0.05" else "p >= 0.05")
  nt <- "DEVIATION: bootstrap outside PLINK (§10); invalid draws (a signal SNP monomorphic) counted as not-same"
  addr(sprintf("S5.Q.boot.%s.share", tn), "Rule 3(ii); §5", est, "1000G EUR 503 (i8 extract, read in place)", "503", "aware", "2nd", num(share), n = N_BOOT, unit = "share",
       notes = sprintf("%s; valid draws %d; share among valid %.4f", nt, sum(!is.na(pb)), mean(same[!is.na(pb)])))
  addr(sprintf("S5.Q.boot.%s.p_lo", tn), "Rule 3(ii); §5", "2.5th percentile of the bootstrap p (aware, 2nd)", "1000G EUR 503", "503", "aware", "2nd", num(qs[1]), n = sum(!is.na(pb)), unit = "p")
  addr(sprintf("S5.Q.boot.%s.p_hi", tn), "Rule 3(ii); §5", "97.5th percentile of the bootstrap p (aware, 2nd)", "1000G EUR 503", "503", "aware", "2nd", num(qs[2]), n = sum(!is.na(pb)), unit = "p")
  addr(sprintf("S5.robust.Q.%s", tn), "Rule 3(ii)", "robust = same conclusion on 503 and 489 AND >= 90% of 5,000 503-panel draws", "", "both", value = robust,
       notes = sprintf("503 p %.4g; 489 p %.4g; panels agree %s; bootstrap share %.4f", pf, p489, panels_agree, share))
  boot_rows[[tn]] <- list(share = share, lo = qs[1], hi = qs[2], robust = robust, p503 = pf, p489 = p489, agree = panels_agree)
  say("RULE 3(ii) Q %s: 503 p %.4g, 489 p %.4g (agree %s); bootstrap share %.4f (p 2.5-97.5%% %.3g-%.3g) -> robust %s", tn, pf, p489, panels_agree, share, qs[1], qs[2], robust)
}
t5 <- A$two[["503"]]$aware_2nd; t9 <- A$two[["489"]]$aware_2nd
same2 <- !is.na(D$p_two) & ((D$p_two < 0.05) == (t5$p < 0.05)) & (sign(D$d_two) == sign(t5$d))
share2 <- sum(same2) / N_BOOT
agree2 <- ((t5$p < 0.05) == (t9$p < 0.05)) && sign(t5$d) == sign(t9$d)
robust2 <- agree2 && share2 >= 0.90
qs2 <- quantile(D$p_two, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
rq <- quantile(D$r_AB, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
addr("S5.two.503.boot_share", "Rule 3(ii); §5 two-signal re-test", sprintf("share of %d draws with the full-panel conclusion: side of p 0.05 (p %.4g) and sign of kA - kB (%+d)", N_BOOT, t5$p, sign(t5$d)),
     "1000G EUR 503 (i8, in place)", "503", "aware", "2nd", num(share2), n = N_BOOT, unit = "share",
     notes = sprintf("DEVIATION: bootstrap outside PLINK (§10); bootstrap p 2.5-97.5%%: %.3g-%.3g; valid %d", qs2[1], qs2[2], sum(!is.na(D$p_two))))
addr("S5.two.489.boot_share", "Rule 3(ii); §5 two-signal re-test", "bootstrap runs on the 503 panel only", "", "489", value = "NOT RUN: Rule 3(ii) bootstraps the 503 panel only")
addr("S5.two.boot.r_AB.ci", "§5 two-signal re-test", "bootstrap 2.5-97.5 percentile interval of r(rs409170, rs62068372), 503", "1000G EUR 503", "503", value = num(A$two[["503"]]$r_AB),
     ci_lo = num(rq[1]), ci_hi = num(rq[2]), n = sum(!is.na(D$r_AB)), unit = "r", notes = sprintf("interval includes 0: %s", rq[1] <= 0 && rq[2] >= 0))
addr("S5.robust.two", "Rule 3(ii)", "two-signal re-test robust = same conclusion (side of p 0.05 and sign) on 503 and 489 AND >= 90% of draws", "", "both", value = robust2,
     notes = sprintf("503 p %.4g sign %+d; 489 p %.4g sign %+d; share %.4f", t5$p, sign(t5$d), t9$p, sign(t9$d), share2))
say("RULE 3(ii) two-signal: 503 p %.4g sign %+d; 489 p %.4g sign %+d; share %.4f -> robust %s; r_AB boot %.3f-%.3f", t5$p, sign(t5$d), t9$p, sign(t9$d), share2, robust2, rq[1], rq[2])
# conditional colocalization at A|B and B|A (Rule 3(ii): side of PP.H4 0.80)
cc <- A$clause
robust_cc <- list()
for (sx in c("A", "B")) {
  h5 <- cc[["503"]][[paste0("H4", sx)]]; h9 <- cc[["489"]][[paste0("H4", sx)]]; hb <- D[[paste0("H4_", sx)]]
  sh <- sum(!is.na(hb) & ((hb >= 0.80) == (h5 >= 0.80))) / N_BOOT
  ag <- (h5 >= 0.80) == (h9 >= 0.80)
  robust_cc[[sx]] <- ag && sh >= 0.90
  sg <- if (sx == "A") SIG_A else SIG_B
  addr(sprintf("S5.coloc_clause.%s_%s.boot_share", sx, sg), "Rule 3(ii) (conditional colocalization: side of PP.H4 0.80)", sprintf("share of %d draws whose conditional PP.H4 at %s is on the full-panel side of 0.80 (full 503 PP.H4 %.4f)", N_BOOT, sg, h5),
       "1000G EUR 503 (i8, in place)", "503", value = num(sh), n = N_BOOT, unit = "share",
       notes = sprintf("bootstrap PP.H4 median %.3f, 2.5-97.5%% %.3f-%.3f; draws with PP.H4 >= 0.80: %d", median(hb, na.rm = TRUE), quantile(hb, 0.025, na.rm = TRUE), quantile(hb, 0.975, na.rm = TRUE), sum(hb >= 0.80, na.rm = TRUE)))
  addr(sprintf("S5.robust.coloc_cond.%s_%s", sx, sg), "Rule 3(ii)", "conditional colocalization robust = same side of 0.80 on 503 and 489 AND >= 90% of draws", "", "both", value = robust_cc[[sx]],
       notes = sprintf("503 PP.H4 %.4f; 489 PP.H4 %.4f; share %.4f", h5, h9, sh))
  say("RULE 3(ii) conditional coloc %s (%s): 503 %.4f, 489 %.4f, share %.4f -> robust %s", sx, sg, h5, h9, sh, robust_cc[[sx]])
}
# Rule 3(iv) clause
fires <- list()
for (sx in c("A", "B")) {
  zz <- c(cc[["503"]][[paste0("z", sx)]], cc[["489"]][[paste0("z", sx)]]); hh <- c(cc[["503"]][[paste0("H4", sx)]], cc[["489"]][[paste0("H4", sx)]])
  fires[[sx]] <- any(abs(zz) >= 1.96) && any(hh < 0.80)
  sg <- if (sx == "A") SIG_A else SIG_B
  addr(sprintf("S5.coloc_clause.%s_%s.fires", sx, sg), "Rule 3(iv)", "clause fires: trait associated at the signal (conditional |z| >= 1.96) but PP.H4 < 0.80 on either panel", "", "both",
       value = fires[[sx]], notes = sprintf("|z| 503 %.3f, 489 %.3f; PP.H4 503 %.4f, 489 %.4f", abs(zz[1]), abs(zz[2]), hh[1], hh[2]))
}

res <- rbind(R2, rbindlist(RES))
stopifnot(!anyDuplicated(res$quantity_id))
num_na <- res[value %in% c("NA", "") | is.na(value)]
say("results: %d rows; rows with empty/NA value: %d (%s)", nrow(res), nrow(num_na), paste(head(num_na$quantity_id, 8), collapse = ", "))
fwrite(res, file.path(OUT, "results.tsv"), sep = "\t", quote = TRUE)

# ---------------------------------------------------------------- branches.tsv (wording rules of §5, Rule 3, Rule 4/4b)
f2 <- function(x) formatC(x, digits = 2, format = "fg", flag = "#")
fp <- function(p) if (p < 0.001) formatC(p, digits = 1, format = "e") else formatC(p, digits = 2, format = "fg", flag = "#")
f3 <- function(x) formatC(x, digits = 3, format = "fg", flag = "#")
B <- list()
addb <- function(rule_ref, inputs, branch, filled_template, qualifiers) B[[length(B) + 1]] <<- data.table(rule_ref, inputs, branch, filled_template, qualifiers)
ldq <- function(rb) if (rb) "" else "sensitive to the reference LD"
for (tn in c("cysgly", "eGFR")) {
  o <- qv(tn, "all5", "503", "aware_2nd"); br <- boot_rows[[tn]]
  lab <- if (tn == "cysgly") "cysteinylglycine" else "eGFR"
  inputs <- sprintf("S5.Q.%s.all5.503.aware_2nd: Q %.3f, df %d, p %.4g (p_ign_2nd %.4g; p_aware_1st %.4g; 489 p %.4g; lower-tail p %.4g); secondary sig123 503 p %.4g; S5.Q.boot.%s.share %.4f",
                    tn, o$Q, o$df, o$p_upper, qv(tn, "all5", "503", "ign_2nd")$p_upper, qv(tn, "all5", "503", "aware_1st")$p_upper, br$p489, o$p_lower,
                    qv(tn, "sig123", "503", "aware_2nd")$p_upper, tn, br$share)
  if (o$p_upper >= 0.05) {
    br_lab <- "p >= 0.05"
    tmpl <- sprintf("no detectable difference in %s per SD of DPEP1 across its signals (Q %s, df %d, p %s)", lab, f2(o$Q), o$df, fp(o$p_upper))
  } else {
    br_lab <- "p < 0.05"
    tmpl <- sprintf("%s per SD of DPEP1 differed across the signals (Q %s, df %d, p %s)", lab, f2(o$Q), o$df, fp(o$p_upper))
  }
  qual <- c("nominal: PENDING_HOLM (Rule 4 family member; Holm by W9)",
            if (o$p_upper >= 0.05) "Rule 4: a criterion met by not rejecting uses the unadjusted p" else NULL,
            if (ldq(br$robust) != "") sprintf("%s (Rule 3(ii): panels agree %s, bootstrap share %.4f)", ldq(br$robust), br$agree, br$share) else sprintf("Rule 3(ii) robust (panels agree; bootstrap share %.4f): no LD qualifier", br$share),
            "Rule 4b: non-rejection worded as 'no detectable difference'; never 'proportional', 'coherent', 'consistent with', 'similar'",
            "Rule 3(iv): all-signal GLS -> per-signal z and PP.H4 in the supplementary table (S5.persig.*); no clause in the sentence",
            {
              ps <- A$persig[trait == tn & panel %in% c("503", "489") & !is.na(PPH4)]
              sprintf("supplementary table content (for the honesty lens): conditional PP.H4 of %s with the protein at signals 1-5 = 503 [%s], 489 [%s]; %d of %d cells < 0.80",
                      lab, paste(sprintf("%.3f", ps[panel == "503"][order(signal)]$PPH4), collapse = ", "), paste(sprintf("%.3f", ps[panel == "489"][order(signal)]$PPH4), collapse = ", "),
                      sum(ps$PPH4 < 0.80), nrow(ps))
            },
            sprintf("printed beside, not part of the §5 template (§5 declares no lower-tail clause): lower-tail p %s", fp(o$p_lower)),
            "placement fixed: one short main-text paragraph after ARG1; detail in the supplement")
  addb(sprintf("§5 wording rule: %s Q", if (tn == "cysgly") "metabolite" else "eGFR"), inputs, br_lab, tmpl, paste(qual, collapse = " | "))
}
# SuSiE admissibility
sm <- A$smax
addb("§5 SuSiE evidence; Rule 3(i)", sprintf("S5.s.protein.max %.4f; S5.s.cysgly.max %.4f (max over susie/p1fr/p1np/p2np x 503/489)", sm[trait == "protein"]$s_max_all, sm[trait == "cysgly"]$s_max_all),
     if (A$admissible) "admissible (s < 0.10 for both traits)" else "inadmissible (s >= 0.10 for at least one trait)",
     if (A$admissible) "SuSiE-based evidence for DPEP1-cysteinylglycine may be cited (s < 0.10 for both traits)" else
       "SuSiE-based evidence (reference-LD PP.H4, credible sets) for DPEP1-cysteinylglycine is reported in the supplement, labelled inadmissible under Rule 3(i), and never cited as support",
     sprintf("locally computed s is a proxy for the pipeline fit (Rule 3(i)); protein s >= 0.10 on every set x panel (%.3f-%.3f); cys-gly s straddles 0.10 across sets (%.3f-%.3f), so the protein alone decides; eGFR s is W6b's",
             min(A$s_tab[trait == "protein"]$s), max(A$s_tab[trait == "protein"]$s), min(A$s_tab[trait == "cysgly"]$s), max(A$s_tab[trait == "cysgly"]$s)))
# Gene attribution
sepv <- max(A$sep[["503"]]); sepd <- sepv < 0.10
gnm <- vapply(A$gn, function(x) x$majority, logical(1))
gin <- sprintf("S5.sep.maxr2 %.4f (503; 489 %.4f) -> separated %s; S5.genename assigned/total: %s", sepv, max(A$sep[["489"]]), sepd,
               paste(sprintf("%s %d/%d (majority %s)", names(A$gn), vapply(A$gn, `[[`, 0, "assigned"), vapply(A$gn, `[[`, 0, "total"), gnm), collapse = "; "))
named <- sepd && all(gnm)
addb("§5 gene attribution; Rule 3(iii)", gin, if (named) "separated AND gene-naming rule met under all three annotations" else
       sprintf("not named: %s", paste(c(if (!sepd) "signals not separated (max r^2 >= 0.10)" else NULL, if (!all(gnm)) sprintf("gene-naming rule not met under %s", paste(names(gnm)[!gnm], collapse = ", ")) else NULL), collapse = "; ")),
     if (named) "DPEP1 is named as the gene" else "the attribution is unresolved (DPEP1 vs CHMP1A)",
     "'unresolved' travels with every attribution sentence (Principle 0); the SM grid prints the restricted counts beside its 48-cell counts")
# SM l.15 two-signal sentence
t5c <- A$two[["503"]]
inp15 <- sprintf("S5.two.503: kA %.4f, kB %.4f, diff %.4f (95%% CI %.4f to %.4f), p_aware_2nd %.4g (p_ign %.4g, p_aware_1st %.4g, p_cojo %.4g); 489 p %.4g; boot_share %.4f; coloc clause A: z 503 %.2f / 489 %.2f, PP.H4 503 %.3f / 489 %.3f; B: z %.2f / %.2f, PP.H4 %.3f / %.3f",
                 t5$kA, t5$kB, t5$d, t5$d - 1.959964 * t5$se_d, t5$d + 1.959964 * t5$se_d, t5$p, t5c$ign_2nd$p, t5c$aware_1st$p, t5c$cojo$p, t9$p, share2,
                 cc[["503"]]$zA, cc[["489"]]$zA, cc[["503"]]$H4A, cc[["489"]]$H4A, cc[["503"]]$zB, cc[["489"]]$zB, cc[["503"]]$H4B, cc[["489"]]$H4B)
tm15 <- if (t5$p >= 0.05)
  sprintf("the two-signal ratio test found no detectable difference in cysteinylglycine per SD of DPEP1 between rs409170 and rs62068372 (%s and %s; difference %s, 95%% CI %s to %s; p %s)",
          f3(t5$kA), f3(t5$kB), f3(t5$d), f3(t5$d - 1.959964 * t5$se_d), f3(t5$d + 1.959964 * t5$se_d), fp(t5$p)) else
  sprintf("cysteinylglycine per SD of DPEP1 differed between rs409170 and rs62068372 (%s vs %s; p %s)", f3(t5$kA), f3(t5$kB), fp(t5$p))
clause_txt <- c()
for (sx in c("A", "B")) if (fires[[sx]]) {
  sg <- if (sx == "A") SIG_A else SIG_B
  hmin <- min(cc[["503"]][[paste0("H4", sx)]], cc[["489"]][[paste0("H4", sx)]])
  clause_txt <- c(clause_txt, sprintf("the cysteinylglycine association at signal %s did not colocalize with the protein (PP.H4 %s), so the ratio there is not anchored to the protein signal", sg, f2(hmin)))
}
q15 <- c(if (robust2) sprintf("Rule 3(ii) two-signal ratio test robust (panels agree; share %.4f)", share2) else sprintf("sensitive to the reference LD (two-signal test; share %.4f)", share2),
         if (length(clause_txt)) paste0("Rule 3(iv) colocalization clause FIRES: ", paste(clause_txt, collapse = "; ")) else "Rule 3(iv) clause does not fire",
         sprintf("conditional colocalization at A %s; at B %s (Rule 3(ii))", if (robust_cc$A) "robust" else "NOT robust: sensitive to the reference LD", if (robust_cc$B) "robust" else "NOT robust: sensitive to the reference LD"),
         "not a Rule 4 family member (unadjusted p)", "Rule 4b: never 'proportional'; the existing ratio_proportional flag is the |dev| <= 3.5 rule, not a test")
st15 <- if (robust_cc$A && robust_cc$B && robust2 && !length(clause_txt)) "SM l.15 stands as updated wording" else
  "SM l.15 must change: its premise that conditional colocalization and the two-signal ratio test 'agreed between panels' holds for the ratio test only"
addb("§5 two-signal re-test; SM attribution l.15 (P31R_SM_attribution.md md5 e7b6e046)", inp15,
     sprintf("%s; %s", if (t5$p >= 0.05) "p >= 0.05 (Rule 4b non-rejection)" else "p < 0.05", st15), tm15, paste(q15, collapse = " | "))
BR <- rbindlist(B)
fwrite(BR, file.path(OUT, "branches.tsv"), sep = "\t", quote = TRUE)
for (i in seq_len(nrow(BR))) say("BRANCH %s -> %s :: %s", BR$rule_ref[i], BR$branch[i], BR$filled_template[i])
say("==== S5_04_assemble done ====")
