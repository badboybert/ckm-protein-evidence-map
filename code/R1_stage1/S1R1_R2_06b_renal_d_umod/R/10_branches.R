# 10_branches.R -- branches.tsv from the redacted texts: section 6(d) classes, D6 (sec 11) + the renal OR rule and qualifier list
# (reads W6a results.tsv; md5 recorded at read and re-checked by 12_end.R), the UMOD control branch, Rule 3(i)/(ii)/(iv), Rule 4.
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
say("==== W6b 10_branches ====")
R6b <- fread(file.path(OUT, "results.tsv"), colClasses = "character", na.strings = NULL)
w6a_path <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_06a_renal/results.tsv")
if (!file.exists(w6a_path)) stop("W6a results.tsv not present: wait for W6a")
w6a_md5 <- md5(w6a_path)
writeLines(c(w6a_path, w6a_md5), file.path(WORK, "w6a_results_read.txt"))
say("W6a results.tsv read: md5 %s", w6a_md5)
R6a <- fread(w6a_path, colClasses = "character", na.strings = NULL)
g <- function(tab, id, col = "value") { x <- tab[quantity_id == id][[col]]; if (length(x) != 1 || is.na(x) || x == "") stop("missing ", id, " ", col); x }
gn <- function(tab, id, col = "value") as.numeric(g(tab, id, col))
f3 <- function(x) trimws(formatC(x, format = "f", digits = 3))   # trimws: an infinite bound (seven-fit (d) table has a zero cell) otherwise prints as " Inf"
f2 <- function(x) formatC(x, format = "f", digits = 2)
fp <- function(p) trimws(if (p < 0.001) formatC(p, format = "e", digits = 1) else formatC(p, format = "g", digits = 2))
B <- list(); addb <- function(rule_ref, inputs, branch, filled, qual = "") B[[length(B) + 1]] <<- data.table(rule_ref = rule_ref, inputs = inputs, branch = branch, filled_template = filled, qualifiers = qual)

# ---- 6(d) per candidate ----
S <- readRDS(file.path(WORK, "d_summary.rds"))
for (i in seq_len(nrow(S$S))) {
  x <- S$S[i]
  addb("6(d); Rule 3(i)", sprintf("%s: s_protein %s; s_outcome %s; s_max %s", x$cid, if (is.na(x$s_protein)) "not computable" else f3(x$s_protein), f3(x$s_outcome), f3(x$s_max)),
       if (x$reclass) "reclassified as not a candidate" else "kept as a candidate",
       sprintf("%s-%s: %s", x$gene, x$outcome, x$why),
       sprintf("s is a locally computed proxy for the pipeline fit (Rule 3(i)); labelled as one; s_max is over all seven distinct fits (Stage 2 decision); four-fit sensitivity s_max %s", f3(x$s_max_4fit)))
}
fd <- S$fd; fk <- S$fk
addb("6(d)", sprintf("%d/%d vs %d/%d; kept version %d/%d vs %d/%d", fd$a, fd$a + fd$b, fd$c, fd$c + fd$d, fk$a, fk$a + fk$b, fk$c, fk$c + fk$d),
     "(d) Fisher test and OR recomputed, both ways",
     sprintf("(d): OR %s (95%% CI %s-%s), p %s; 'kept' version: OR %s (95%% CI %s-%s), p %s", f2(fd$or), f3(fd$lo), f3(fd$hi), fp(fd$p), f2(fk$or), f3(fk$lo), f3(fk$hi), fp(fk$p)),
     "conditional-MLE OR, exact CI (fisher.test's equations solved at tol 1e-13: W6a exact_cond_or.R; FIXLOG F1)")

# ---- D6 (sec 11) and the renal OR wording (sec 6) ----
p_or <- gn(R6a, "S6.primary.or_mle"); p_lo <- gn(R6a, "S6.primary.or_lo"); p_hi <- gn(R6a, "S6.primary.or_hi"); p_p <- gn(R6a, "S6.primary.p")
my_prim <- c(gn(R6b, "S6.d.check_primary.or"), gn(R6b, "S6.d.check_primary.lo"), gn(R6b, "S6.d.check_primary.hi"))
log_assert("S6.primary.W6a_vs_W6b_check", "identical OR and exact CI", sprintf("W6a %.6f (%.6f-%.6f) / W6b %.6f (%.6f-%.6f)", p_or, p_lo, p_hi, my_prim[1], my_prim[2], my_prim[3]),
           max(abs(c(p_or, p_lo, p_hi) - my_prim)) < 1e-6, "W6a S6.primary.or_* vs W6b S6.d.check_primary.* (same 2x2 from independent reads)")
d6_over <- round(p_lo, 3) > 1 && round(fd$lo, 3) > 1
base <- if (d6_over) "over-represented at the pair level" else "proportionally more frequent"
addb("11 D6", sprintf("primary pair OR exact lower bound %s (W6a S6.primary.or_lo); (d) exact lower bound %s (W6b S6.d.lo)", f3(p_lo), f3(fd$lo)),
     if (d6_over) "both lower bounds exceed 1.000 at three decimals -> 'over-represented at the pair level'" else "not both lower bounds exceed 1.000 at three decimals -> 'proportionally more frequent'",
     base, "the sec 6 qualifier rule applies in both cases")
inc1 <- function(lo, hi) lo <= 1 && hi >= 1
lo_lo <- gn(R6a, "S6.b.locus_outcome.or_lo"); lo_hi <- gn(R6a, "S6.b.locus_outcome.or_hi")
c_lo <- gn(R6a, "S6.c.pooled_lo"); c_hi <- gn(R6a, "S6.c.pooled_hi")
a_lo <- gn(R6a, "S6.b.a_pair.or_lo"); a_hi <- gn(R6a, "S6.b.a_pair.or_hi")
lst <- c(if (inc1(lo_lo, lo_hi) || inc1(c_lo, c_hi)) "grouping", if (inc1(a_lo, a_hi)) "the acceptance rule", if (inc1(fd$lo, fd$hi)) "the LD guard")
lst_txt <- if (!length(lst)) "" else if (length(lst) == 1) lst else paste0(paste(lst[-length(lst)], collapse = ", "), " and ", lst[length(lst)])
inputs <- sprintf("pair OR %s (%s-%s), p %s [W6a]; locus-outcome exact CI %s-%s [W6a S6.b.locus_outcome]; (c) pooled CI %s-%s [W6a S6.c]; (a) exact CI %s-%s [W6a S6.b.a_pair]; (d) exact CI %s-%s [W6b S6.d]",
                  f2(p_or), f3(p_lo), f3(p_hi), fp(p_p), f3(lo_lo), f3(lo_hi), f3(c_lo), f3(c_hi), f3(a_lo), f3(a_hi), f3(fd$lo), f3(fd$hi))
if (p_or <= 1) {
  addb("6 renal OR wording", inputs, "pair OR <= 1", sprintf("not more frequent among kidney pairs (OR %s, 95%% CI %s-%s)", f2(p_or), f3(p_lo), f3(p_hi)), "")
} else {
  filled <- sprintf("kidney associations were %s among the multi-signal candidates (pair-level OR %s, 95%% CI %s-%s; p %s)%s", base, f2(p_or), f3(p_lo), f3(p_hi), fp(p_p),
                    if (nzchar(lst_txt)) sprintf(", sensitive to %s", lst_txt) else "")
  addb("6 renal OR wording; 11 D6", inputs, sprintf("pair OR > 1; D6 base phrase '%s'; variants whose CI includes 1: %s", base, if (length(lst)) paste(lst, collapse = "; ") else "none"),
       filled, sprintf("renal 'sensitive to {list}' (Principle 0) = %s", if (nzchar(lst_txt)) lst_txt else "none"))
}

# ---- UMOD ----
U <- readRDS(file.path(WORK, "umod", "umod_summary.rds"))
q5 <- U$G503$aware_2nd
h4 <- c(A503 = gn(R6b, "S6.umod.H4.A.503"), B503 = gn(R6b, "S6.umod.H4.B.503"), A489 = gn(R6b, "S6.umod.H4.A.489"), B489 = gn(R6b, "S6.umod.H4.B.489"))
fail <- c(if (any(h4 < 0.80)) sprintf("conditional PP.H4 below 0.80 at %s", paste(names(h4)[h4 < 0.80], collapse = ", ")),
          if (q5$p_upper < 0.05) sprintf("heterogeneity detected across the 5 protein signals (Q %s, df %d, p %s)", f2(q5$Q), q5$df, fp(q5$p_upper)))
uin <- sprintf("cond PP.H4 A/B 503 %s/%s, 489 %s/%s; all-signal Q %s, df %d, p %s (503, aware, 2nd order, unadjusted); 489 p %s; bootstrap share %s",
               f3(h4[1]), f3(h4[2]), f3(h4[3]), f3(h4[4]), f2(q5$Q), q5$df, fp(q5$p_upper), fp(U$G489$aware_2nd$p_upper), f3(U$share))
qual_u <- c(if (!U$robust) "sensitive to the reference LD (Rule 3(ii))", "nominal: PENDING_HOLM (UMOD-eGFR all-signal Q is a family member; the control criterion itself uses the unadjusted p, Rule 4)",
            "the label means only that no departure was detected (Rule 4b); it does not mean proportionality was shown")
if (!length(fail)) {
  addb("6 UMOD control", uin, "behaved as a control",
       sprintf("both declared signals colocalized with eGFR (PP.H4 %s and %s on the 503 panel; %s and %s on the 489 panel), and no heterogeneity was detected in eGFR per SD of UMOD across its 5 protein signals (Q %s, df %d, p %s)",
               f2(h4[1]), f2(h4[2]), f2(h4[3]), f2(h4[4]), f2(q5$Q), q5$df, fp(q5$p_upper)), paste(qual_u, collapse = "; "))
} else {
  addb("6 UMOD control", uin, "did not behave as a clean control",
       sprintf("UMOD did not behave as a clean control (%s), which limits how far signal-level departures at other loci can be read as evidence against the protein route", paste(fail, collapse = "; ")),
       paste(qual_u, collapse = "; "))
}
addb("Rule 3(ii) (UMOD primary cell: all-signal GLS, 503)", sprintf("503 p %s, 489 p %s, share %s of 5,000 draws", fp(q5$p_upper), fp(U$G489$aware_2nd$p_upper), f3(U$share)),
     if (U$robust) "robust" else "not robust", if (U$robust) "no LD qualifier" else "adds 'sensitive to the reference LD'", "DEVIATION: bootstrap outside PLINK (sec 10)")
# Rule 3(iv): two-signal contrast clause at A and B
for (X in c("A", "B")) {
  z <- c(gn(R6b, sprintf("S6.umod.z.%s.503", X)), gn(R6b, sprintf("S6.umod.z.%s.489", X)))
  hh <- c(h4[[paste0(X, "503")]], h4[[paste0(X, "489")]])
  fires <- any(abs(z) >= 1.96) && any(hh < 0.80)
  addb("Rule 3(iv) (UMOD two-signal contrast)", sprintf("signal %s: eGFR conditional z %s (503) / %s (489); PP.H4 %s / %s", X, f2(z[1]), f2(z[2]), f3(hh[1]), f3(hh[2])),
       if (fires) "fires" else if (all(abs(z) < 1.96)) "does not fire (no association at the signal)" else "does not fire (PP.H4 >= 0.80 on both panels)",
       if (fires) sprintf("the eGFR association at signal %s did not colocalize with the protein (PP.H4 %s), so the ratio there is not anchored to the protein signal", X, f2(min(hh))) else "", "")
}
t5 <- U$TS$`503`$res; t9 <- U$TS$`489`$res
addb("6 UMOD secondary two-signal test; Rule 4b; Rule 5", sprintf("503: kA %.4g kB %.4g p_aware_2nd %s (p_ign %s, p_aware_1st %s); 489: p_aware_2nd %s", t5$kA, t5$kB, fp(t5$p_aware_2nd), fp(t5$p_ign), fp(t5$p_aware_1st), fp(t9$p_aware_2nd)),
     if (t5$p_aware_2nd < 0.05) "p < 0.05 (503)" else "p >= 0.05 (503)",
     if (t5$p_aware_2nd < 0.05) sprintf("eGFR per SD of UMOD differed between signals A and B (p %s)", fp(t5$p_aware_2nd)) else sprintf("no detectable difference in eGFR per SD of UMOD between signals A and B (p %s; not rejected)", fp(t5$p_aware_2nd)),
     "secondary; not in the Rule 4 family; Rule 4b wording")
# Rule 3(i) admissibility of the slot loci's SuSiE evidence
for (sl in c("UMOD_eGFR", "UMOD_CKD", "DPEP1_eGFR")) {
  x <- S$S[cid == sl]
  adm <- !is.na(x$s_protein) && x$s_max < 0.10
  addb("Rule 3(i) (SuSiE with reference LD)", sprintf("%s: s_protein %s, s_outcome %s (max over all seven distinct fits; proxy; four-fit sensitivity: s_protein %s, s_outcome %s)", sl, f3(x$s_protein), f3(x$s_outcome), f3(x$s_protein_4fit), f3(x$s_outcome_4fit)),
       if (adm) "admissible (s < 0.10 for both traits)" else "inadmissible (s >= 0.10 for at least one trait)",
       if (adm) sprintf("%s SuSiE PP.H4 admissible", sl) else sprintf("%s SuSiE PP.H4: reported in the supplement, labelled, never cited as support in new text; flagged at first mention", sl),
       "locally computed s for a pipeline fit is a proxy and is labelled as one")
}
addb("Rule 3(i) scope", "slot loci and the 16 candidates only", "DEVIATION",
     "s not computed (pipeline fit on E:) for every other pair",
     "DEVIATION: Rule 3(i) first-mention flag computed only for the slot loci and the 16 candidates; the declared rule covers every pair whose class rests on reference-LD SuSiE; user decision pending at Stage 2")
BB <- rbindlist(B)
fwrite(BB, file.path(OUT, "branches.tsv"), sep = "\t", quote = TRUE)
say("branches.tsv: %d rows; md5 %s", nrow(BB), md5(file.path(OUT, "branches.tsv")))
