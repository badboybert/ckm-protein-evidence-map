#!/usr/bin/env Rscript
# w4a_compare_seen.R -- Stage 1 W4a, step 2 (Rule 10: compute first, compare second).
# Runs AFTER w4a_main.R has written results.tsv and branches.tsv. Reads the FULL pre-declaration's section 4 "Seen" values
# (transcribed below with their line numbers) and writes COMPARE_SEEN.tsv. It never changes a computed value:
# results.tsv is only read here. Post-comparison DIAGNOSTICS (labelled as such, never results) are computed to explain
# differences and written to work/compare_diagnostics.tsv.
#   Rscript w4a_compare_seen.R
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
suppressPackageStartupMessages({ library(data.table); library(coloc) })
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04a_arg1grid")
WORK <- file.path(OUT, "work")
source(file.path(OUT, "w4a_funcs.R"))
log_init(file.path(OUT, "run.log"))
logf("==== W4a compare-seen step start ====")
FULL <- paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md")
AMEND <- paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md")
assert_md5(FULL, "3319dc39740d9f5b3b8822133d9108b7", "prespec_full (opened for COMPARE_SEEN only)")
assert_md5(AMEND, "377508832ea9bd384133e869c456c909", "amend1_full (opened for COMPARE_SEEN only)")
RES <- fread(file = file.path(OUT, "results.tsv"), sep = "\t", colClasses = "character")
logf("results.tsv md5 %s (%d rows) read for comparison", md5f(file.path(OUT, "results.tsv")), nrow(RES))
fl <- readLines(FULL, encoding = "UTF-8"); al <- readLines(AMEND, encoding = "UTF-8")
chk_line <- function(lines, n, needle) if (!grepl(needle, lines[n], fixed = TRUE)) stop(sprintf("STOP: line %d does not contain '%s'", n, needle))
v <- function(id) { x <- RES$value[RES$quantity_id == id]; if (length(x) != 1) stop(paste("missing", id)); as.numeric(x) }

# ---------------- post-comparison diagnostics (explanations only) ----------------
own <- fread(file = file.path(WORK, "own_harmonised_ARG1win.tsv"), sep = "\t")
iA <- which(own$rsid_ref == "rs2781668"); iB <- which(own$rsid_ref == "rs17788484")
r503 <- v("S4.ld.r.503"); r489 <- v("S4.ld.r.489")
tr <- function(col, fref = FALSE) list(bA = own[[paste0("b_", col)]][iA], sA = own[[paste0("s_", col)]][iA], bB = own[[paste0("b_", col)]][iB],
  sB = own[[paste0("s_", col)]][iB], fA = if (fref) own$af_ref[iA] else own[[paste0("f_", col)]][iA], fB = if (fref) own$af_ref[iB] else own[[paste0("f_", col)]][iB])
D <- list()
dg <- function(id, val, note) { D[[length(D) + 1]] <<- data.frame(diag_id = id, value = sprintf("%.7g", val), note = note); logf("DIAG %s = %.7g (%s)", id, val, note) }
for (fr in c(FALSE, TRUE)) {
  lab <- if (fr) "reffreq" else "ownfreq"
  P <- est_trait(tr("P", fr), r503, "cond", "beta"); Y <- est_trait(tr("fg", fr), r503, "cond", "beta"); M <- est_trait(tr("arg", fr), r503, "cond", "beta")
  a <- ab_test(Y, P, "2nd", TRUE)
  dg(sprintf("diag.cojo_%s.FinnGen.cond.ratioB", lab), a$RB, "beta-space conditioning, T2D per SD ARG1 at B|A")
  dg(sprintf("diag.cojo_%s.FinnGen.cond.p_aware_2nd", lab), a$p, "beta-space, aware, 2nd order")
  dg(sprintf("diag.cojo_%s.FinnGen.cond.p_ign_2nd", lab), ab_test(Y, P, "2nd", FALSE)$p, "beta-space, ign, 2nd order")
  dg(sprintf("diag.cojo_%s.t3.FinnGen.cond.p_aware_2nd", lab), ab_test(Y, M, "2nd", TRUE)$p, "test (iii) beta-space aware 2nd")
  dg(sprintf("diag.cojo_%s.t3.FinnGen.cond.p_aware_1st", lab), ab_test(Y, M, "1st", TRUE)$p, "test (iii) beta-space aware 1st")
}
A <- ab_test(marg_est(own$b_arg[iA], own$s_arg[iA], own$b_arg[iB], own$s_arg[iB], r503), marg_est(own$b_P[iA], own$s_P[iA], own$b_P[iB], own$s_P[iB], r503), "2nd", FALSE)
fi <- fieller(A$RA, A$RB, A$vA, A$vB, 0); ld <- logdelta_ci(A$RA, A$RB, A$vA, A$vB, 0)
dg("diag.arg.marg.ror_ign_fieller_lo", fi$lo, "ign Fieller"); dg("diag.arg.marg.ror_ign_fieller_hi", fi$hi, "ign Fieller")
dg("diag.arg.marg.ror_ign_logdelta_lo", ld[1], "ign log-scale delta"); dg("diag.arg.marg.ror_ign_logdelta_hi", ld[2], "ign log-scale delta")
Aw <- ab_test(marg_est(own$b_arg[iA], own$s_arg[iA], own$b_arg[iB], own$s_arg[iB], r503), marg_est(own$b_P[iA], own$s_P[iA], own$b_P[iB], own$s_P[iB], r503), "2nd", TRUE)
ldw <- logdelta_ci(Aw$RA, Aw$RB, Aw$vA, Aw$vB, Aw$cov)
dg("diag.arg.marg.ror_aware_logdelta_lo", ldw[1], "aware log-scale delta"); dg("diag.arg.marg.ror_aware_logdelta_hi", ldw[2], "aware log-scale delta")
# ornithine z-ratio k (engine two_signal convention: k = conditional z_M / conditional z_P)
for (pn in c("503", "489")) { r <- if (pn == "503") r503 else r489
  P <- cond_z(own$b_P[iA], own$s_P[iA], own$b_P[iB], own$s_P[iB], r); O <- cond_z(own$b_orn[iA], own$s_orn[iA], own$b_orn[iB], own$s_orn[iB], r)
  dg(sprintf("diag.orn.zratio.kA.%s", pn), O$zA / P$zA, "conditional z_M(A|B)/z_P(A|B)"); dg(sprintf("diag.orn.zratio.kB.%s", pn), O$zB / P$zB, "conditional z_M(B|A)/z_P(B|A)") }
# colocalization on ALL 3,067 harmonised SNPs (palindromes kept), 503 dosage LD (dosage r == PLINK r to 5e-7, run.log)
G <- matrix(readBin(paste0(.P31R_V3, "/V4_local/raw/l3-finemap/ARG1_6q23__ref_A1dosage_503xP.i8"), "integer", size = 1, signed = TRUE,
                    n = 503 * 3067), nrow = 503, byrow = TRUE)
rA <- cor(G, G[, iA])[, 1]; rB <- cor(G, G[, iB])[, 1]; rm(G)
NP <- 33455; NM <- 8299; NFG <- 400197; SFG <- 65085 / 400197
afok <- abs(own$f_fg - own$af_ref) <= 0.2
for (set in c("all3067", "p1np")) {
  idx <- if (set == "all3067") seq_len(nrow(own)) else which(own$palindromic == 0)
  kA <- match(iA, idx); kB <- match(iB, idx)
  Pc <- cond_window(own$b_P[idx], own$s_P[idx], rB[idx], kB)
  Mc <- cond_window(own$b_arg[idx], own$s_arg[idx], rB[idx], kB)
  bt <- own$b_fg[idx]; st <- own$s_fg[idx]; bad <- is.na(afok[idx]) | !afok[idx]; bt[bad] <- NA; st[bad] <- NA
  Tc <- cond_window(bt, st, rB[idx], kB)
  ids <- own$snp[idx]
  for (p12 in c(1e-5, 1e-6)) {
    dg(sprintf("diag.coloc.%s.PM.A_given_B.p12_%g", set, p12), coloc_pair(Pc$b, Pc$s, NP, Mc$b, Mc$s, NM, "quant", NA, ids, p12)$H4, "arginine x protein, both conditioned on B")
    dg(sprintf("diag.coloc.%s.PD.A_given_B.cc.p12_%g", set, p12), coloc_pair(Pc$b, Pc$s, NP, Tc$b, Tc$s, NFG, "cc", SFG, ids, p12)$H4, "T2D (cc, W 0.2^2) x protein, both conditioned on B")
    dg(sprintf("diag.coloc.%s.PD.A_given_B.quantW015.p12_%g", set, p12), coloc_pair(Pc$b, Pc$s, NP, Tc$b, Tc$s, NFG, "quant", NA, ids, p12)$H4,
       "T2D treated with prior sd 0.15 (quant sdY 1), both conditioned on B")
    dg(sprintf("diag.coloc.%s.PM.marginal.p12_%g", set, p12), coloc_pair(own$b_P[idx], own$s_P[idx], NP, own$b_arg[idx], own$s_arg[idx], NM, "quant", NA, ids, p12)$H4, "marginal")
    dg(sprintf("diag.coloc.%s.PD.marginal.cc.p12_%g", set, p12), coloc_pair(own$b_P[idx], own$s_P[idx], NP, bt, st, NFG, "cc", SFG, ids, p12)$H4, "marginal, cc")
  }
}
DT <- rbindlist(D); fwrite(DT, file.path(WORK, "compare_diagnostics.tsv"), sep = "\t")
dv <- function(id) as.numeric(DT$value[DT$diag_id == id])

# ---------------- transcribed seen values (full text line numbers) ----------------
S <- list()
sv <- function(qid, seen, seen_label, src, stage1_value, stage1_label, tol = NULL, expl_same = "", expl_diff = "") {
  s <- as.numeric(seen)
  dec <- if (grepl("e", seen, ignore.case = TRUE)) NA else nchar(sub("^[^.]*\\.?", "", seen))
  if (is.null(tol)) tol <- if (is.na(dec)) abs(s) * 0.05 else 0.5 * 10^(-dec) + 1e-12
  d <- abs(as.numeric(stage1_value) - s)
  within <- is.finite(d) && d <= tol
  ex <- if (within) paste0("within printed precision (tol ", signif(tol, 3), ")", if (nzchar(expl_same)) paste0("; ", expl_same) else "") else
    paste0("DIFFERS beyond printed precision (tol ", signif(tol, 3), "): ", if (nzchar(expl_diff)) expl_diff else "unexplained")
  S[[length(S) + 1]] <<- data.frame(quantity_id = qid, seen_value = seen, seen_label = paste0(seen_label, " [", src, "]"),
                                    stage1_value = sprintf("%.7g", as.numeric(stage1_value)), stage1_label = stage1_label,
                                    abs_diff = sprintf("%.3g", d), explanation = ex)
}
nv <- function(qid, seen, seen_label, src, why) S[[length(S) + 1]] <<- data.frame(quantity_id = qid, seen_value = seen, seen_label = paste0(seen_label, " [", src, "]"),
  stage1_value = "NOT COMPUTED", stage1_label = "", abs_diff = "", explanation = why)
F <- "full prespec l."; W <- "weak evidence: value also visible in a definition read before computing (R1_ARG1_rule or outcome-check prespec)"
chk_line(fl, 78, "ign p 0.0336"); chk_line(fl, 336, "aware 0.0365"); chk_line(fl, 366, "0.0445 (503)"); chk_line(fl, 386, "dev_sd 4.011")
sv("S4.grid.FinnGen.cond.p_ign", "0.0336", "Rule 1: WEAKENS frozen on the rule's own statistic (ign p)", paste0(F, 78), v("S4.grid.FinnGen.cond.p_ign"), "FinnGen cond ign 2nd, 503")
sv("S4.ld.r.503", "-0.0445", "r(A,B) 503", paste0(F, 366), v("S4.ld.r.503"), "PLINK --r", expl_same = W)
sv("S4.ld.r.489", "-0.036", "489 r", paste0(F, 366), v("S4.ld.r.489"), "PLINK --r, aligned to a1")
sv("S4.ld.boot.r_lo", "-0.107", "bootstrap 95% lower", paste0(F, 366), v("S4.ld.boot.r_lo"), "R MT seed 20261007, 5,000 draws, type-7 percentile",
   expl_diff = "different bootstrap implementation (scoping RNG, draw count and percentile rule not known to me; seen script not opened); Monte Carlo error of a tail percentile of a skewed r (22 copies of B). Both intervals include 0, so the T2D sign clause fires either way")
sv("S4.ld.boot.r_hi", "0.029", "bootstrap 95% upper", paste0(F, 366), v("S4.ld.boot.r_hi"), "R MT seed 20261007, 5,000 draws, type-7 percentile",
   expl_diff = "as S4.ld.boot.r_lo: implementation/RNG/draw-count difference; not further explained. Same side of 0 (interval includes 0 in both)")
sv("S4.ld.rgrid.crossing.aware", "-0.008", "r-grid crossing, FinnGen cond aware", paste0(F, 368), v("S4.ld.rgrid.crossing.aware"), "grid step 0.0005, linear interpolation")
sv("S4.ld.rgrid.crossing.ign", "-0.006", "r-grid crossing, FinnGen cond ign", paste0(F, 368), v("S4.ld.rgrid.crossing.ign"), "grid step 0.0005, linear interpolation")
sv("S4.arg.marg.ratioA", "-0.7255", "arginine marginal ratio A", paste0(F, 332), v("S4.arg.marg.ratioA"), "b_M/b_P marginal", expl_same = W)
sv("S4.arg.marg.ratioB", "-0.7158", "arginine marginal ratio B", paste0(F, 332), v("S4.arg.marg.ratioB"), "b_M/b_P marginal", expl_same = W)
sv("S4.arg.marg.p_ign", "0.947", "arginine A vs B marginal ign", paste0(F, 333), v("S4.arg.marg.p_ign"), "ign 2nd, 503", expl_same = W)
sv("S4.arg.marg.p_aware_2nd", "0.946", "arginine A vs B marginal aware", paste0(F, 333), v("S4.arg.marg.p_aware_2nd"), "aware 2nd, 503 (primary cell)")
sv("S4.arg.cond.p_ign", "0.952", "arginine A vs B conditional ign", paste0(F, 333), v("S4.arg.cond.p_ign"), "ign 2nd, 503")
sv("S4.arg.cond.p_aware_2nd", "0.953", "arginine A vs B conditional aware", paste0(F, 333), v("S4.arg.cond.p_aware_2nd"), "aware 2nd, 503")
sv("S4.arg.marg.ror", "1.0135", "ratio of ratios A/B", paste0(F, 334), v("S4.arg.marg.ror"), "k_A/k_B marginal")
sv("S4.arg.marg.ror_lo", "0.6831", "ror 95% CI lower (ign)", paste0(F, 334), v("S4.arg.marg.ror_lo"), "Fieller 95%, AWARE 2nd (the pre-declaration's arginine template asks for the aware CI)",
   expl_diff = sprintf("different quantity: seen CI is ign; Stage 1 primary CI is aware Fieller. Post-comparison diagnostics: ign Fieller %.4f-%.4f, ign log-delta %.4f-%.4f, aware log-delta %.4f-%.4f; the seen ign interval is reproduced by %s",
                       dv("diag.arg.marg.ror_ign_fieller_lo"), dv("diag.arg.marg.ror_ign_fieller_hi"), dv("diag.arg.marg.ror_ign_logdelta_lo"), dv("diag.arg.marg.ror_ign_logdelta_hi"),
                       dv("diag.arg.marg.ror_aware_logdelta_lo"), dv("diag.arg.marg.ror_aware_logdelta_hi"),
                       if (abs(dv("diag.arg.marg.ror_ign_logdelta_lo") - 0.6831) < 0.00005) "the ign log-scale delta CI" else if (abs(dv("diag.arg.marg.ror_ign_fieller_lo") - 0.6831) < 0.00005) "the ign Fieller CI" else "neither diagnostic (unexplained)"))
sv("S4.arg.marg.ror_hi", "1.5038", "ror 95% CI upper (ign)", paste0(F, 334), v("S4.arg.marg.ror_hi"), "Fieller 95%, AWARE 2nd",
   expl_diff = "different quantity (ign vs aware CI; see S4.arg.marg.ror_lo). Branch unaffected: p >= 0.05 under both")
sv("S4.grid.FinnGen.cond.ratioA", "-0.1975", "T2D FinnGen cond ratio A", paste0(F, 336), v("S4.grid.FinnGen.cond.ratioA"), "z-space conditional")
sv("S4.grid.FinnGen.cond.ratioB", "-0.0184", "T2D FinnGen cond ratio B", paste0(F, 336), v("S4.grid.FinnGen.cond.ratioB"), "z-space conditional")
sv("S4.grid.FinnGen.cond.p_aware_2nd", "0.0365", "T2D FinnGen cond aware (primary)", paste0(F, 336), v("S4.grid.FinnGen.cond.p_aware_2nd"), "aware 2nd, 503")
sv("S4.grid.FinnGen.cond.p_cojo", "0.0427", "COJO B|A p", paste0(F, 337), v("S4.grid.FinnGen.cond.p_cojo"), "beta-space, own-source allele frequencies, aware 2nd",
   expl_diff = sprintf("convention difference: Stage 1 beta-space uses each source's OWN allele frequency (COJO .ma convention; FinnGen Finnish EAF at B 0.0073 vs EUR 0.022). Post-comparison diagnostic with 503 REFERENCE frequencies: ratio B|A %.4f (seen -0.0285), p aware %.4f, ign %.4f; the seen also lists 0.0462/0.0438 variants. The seen COJO ratio is close to the reference-frequency diagnostic, so the seen COJO most likely used reference frequencies; the residual p difference (0.0427 vs %.4f ign / %.4f aware) is unexplained without the seen script (not opened). Beta-space p is printed beside only",
                       dv("diag.cojo_reffreq.FinnGen.cond.ratioB"), dv("diag.cojo_reffreq.FinnGen.cond.p_aware_2nd"), dv("diag.cojo_reffreq.FinnGen.cond.p_ign_2nd"),
                       dv("diag.cojo_reffreq.FinnGen.cond.p_ign_2nd"), dv("diag.cojo_reffreq.FinnGen.cond.p_aware_2nd")))
sv("S4.grid.FinnGen.marg.ratioA", "-0.1910", "T2D FinnGen marg ratio A", paste0(F, 338), v("S4.grid.FinnGen.marg.ratioA"), "marginal", expl_same = W)
sv("S4.grid.FinnGen.marg.ratioB", "-0.0363", "T2D FinnGen marg ratio B", paste0(F, 338), v("S4.grid.FinnGen.marg.ratioB"), "marginal")
sv("S4.grid.FinnGen.marg.p_ign", "0.0533", "FinnGen marg ign", paste0(F, 338), v("S4.grid.FinnGen.marg.p_ign"), "ign 2nd")
sv("S4.grid.FinnGen.marg.p_aware_2nd", "0.0494", "FinnGen marg aware", paste0(F, 338), v("S4.grid.FinnGen.marg.p_aware_2nd"), "aware 2nd")
sv("S4.grid.DIAMANTE.marg.ratioA", "-0.1138", "DIAMANTE marg ratio A", paste0(F, 340), v("S4.grid.DIAMANTE.marg.ratioA"), "marginal")
sv("S4.grid.DIAMANTE.marg.ratioB", "0.0395", "DIAMANTE marg ratio B", paste0(F, 340), v("S4.grid.DIAMANTE.marg.ratioB"), "marginal")
sv("S4.grid.DIAMANTE.marg.p_ign", "0.0293", "DIAMANTE marg ign", paste0(F, 340), v("S4.grid.DIAMANTE.marg.p_ign"), "ign 2nd")
sv("S4.grid.DIAMANTE.marg.p_aware_2nd", "0.0258", "DIAMANTE marg aware (A2 cell)", paste0(F, 340), v("S4.grid.DIAMANTE.marg.p_aware_2nd"), "aware 2nd")
nv("S4.grid.DIAMANTE.marg.rounding_box", "ign 0.022-0.038; aware 0.019-0.034", "DIAMANTE rounding box", paste0(F, 341), "not declared for W4a in the redacted text (only the SE-from-p variant is); not computed")
sv("S4.grid.DIAMANTE_SEfromP.marg.p_aware_2nd", "0.027", "DIAMANTE SE-from-p variant", paste0(F, 342), v("S4.grid.DIAMANTE_SEfromP.marg.p_aware_2nd"), "SE = |Beta|/z(P), aware 2nd",
   expl_diff = sprintf("cov/SE-order label: the seen 0.027 equals Stage 1's ign value for the same variant (S4.grid.DIAMANTE_SEfromP.marg.p_ign = %.4f); aware is %.4f. Same side of 0.05", v("S4.grid.DIAMANTE_SEfromP.marg.p_ign"), v("S4.grid.DIAMANTE_SEfromP.marg.p_aware_2nd")))
sv("S4.grid.DIAMANTE.cond.ratioA", "-0.1197", "DIAMANTE cond ratio A", paste0(F, 343), v("S4.grid.DIAMANTE.cond.ratioA"), "z-space conditional")
sv("S4.grid.DIAMANTE.cond.ratioB", "0.0478", "DIAMANTE cond ratio B", paste0(F, 343), v("S4.grid.DIAMANTE.cond.ratioB"), "z-space conditional")
sv("S4.grid.DIAMANTE.cond.p_ign", "0.0230", "DIAMANTE cond ign", paste0(F, 343), v("S4.grid.DIAMANTE.cond.p_ign"), "ign 2nd")
sv("S4.grid.DIAMANTE.cond.p_aware_2nd", "0.0261", "DIAMANTE cond aware", paste0(F, 343), v("S4.grid.DIAMANTE.cond.p_aware_2nd"), "aware 2nd")
for (cc in list(c("FinnGen.cond", "-0.179", "-0.347", "-0.011", "0.037", 348), c("FinnGen.marg", "-0.155", "-0.309", "-0.000", "0.049", 349),
                c("DIAMANTE.marg", "-0.153", "-0.288", "-0.019", "0.026", 350), c("DIAMANTE.cond", "-0.168", "-0.315", "-0.020", "0.026", 351))) {
  ln <- paste0(F, cc[6])
  sv(sprintf("S4.grid.%s.diff", cc[1]), cc[2], "A-B difference (aware)", ln, v(sprintf("S4.grid.%s.diff", cc[1])), "aware 2nd")
  sv(sprintf("S4.grid.%s.diff_lo", cc[1]), cc[3], "95% CI lower", ln, v(sprintf("S4.grid.%s.diff_lo", cc[1])), "aware 2nd")
  sv(sprintf("S4.grid.%s.diff_hi", cc[1]), cc[4], "95% CI upper", ln, v(sprintf("S4.grid.%s.diff_hi", cc[1])), "aware 2nd")
  sv(sprintf("S4.grid.%s.p_aware_2nd", cc[1]), cc[5], "p (aware), A-B table", ln, v(sprintf("S4.grid.%s.p_aware_2nd", cc[1])), "aware 2nd")
}
sv("S4.grid.FinnGen.cond.exp_zB", "2.57", "expected |z| at B", paste0(F, 354), abs(v("S4.grid.FinnGen.cond.exp_zB")), "|exp_zB| (Stage 1 prints it signed: negative)")
sv("S4.grid.FinnGen.cond.obs_zB", "-0.24", "observed z at B", paste0(F, 354), v("S4.grid.FinnGen.cond.obs_zB"), "conditional z at B")
sv("S4.grid.FinnGen.marg.exp_zB", "2.63", "expected |z| at B", paste0(F, 355), abs(v("S4.grid.FinnGen.marg.exp_zB")), "|exp_zB| (signed negative)")
sv("S4.grid.FinnGen.marg.obs_zB", "-0.50", "observed z at B", paste0(F, 355), v("S4.grid.FinnGen.marg.obs_zB"), "marginal z at B", expl_same = W)
sv("S4.grid.DIAMANTE.marg.exp_zB", "2.16", "expected |z| at B", paste0(F, 356), abs(v("S4.grid.DIAMANTE.marg.exp_zB")), "|exp_zB| (signed negative)")
sv("S4.grid.DIAMANTE.marg.obs_zB", "0.75", "observed z at B", paste0(F, 356), v("S4.grid.DIAMANTE.marg.obs_zB"), "marginal z at B (2 s.f. Beta/SE)")
sv("S4.grid.FinnGen.cond.power_point", "0.73", "point power", paste0(F, 357), v("S4.grid.FinnGen.cond.power_point"), "Phi(|ez|-1.96)+Phi(-|ez|-1.96)")
sv("S4.grid.DIAMANTE.marg.power_point", "0.58", "point power", paste0(F, 357), v("S4.grid.DIAMANTE.marg.power_point"), "as above")
sv("S4.grid.FinnGen.cond.power_pred", "0.71", "predictive power", paste0(F, 357), v("S4.grid.FinnGen.cond.power_pred"), "Stage 1 definition in results notes")
sv("S4.grid.DIAMANTE.marg.power_pred", "0.56", "predictive power", paste0(F, 357), v("S4.grid.DIAMANTE.marg.power_pred"), "Stage 1 definition in results notes")
sv("S4.t3.FinnGen.cond.503.2nd.p_aware", "0.0462", "test (iii) conditional aware", paste0(F, 359), v("S4.t3.FinnGen.cond.503.2nd.p_aware"), "aware 2nd 503")
sv("S4.t3.FinnGen.marg.503.2nd.p_aware", "0.0620", "test (iii) marginal aware", paste0(F, 359), v("S4.t3.FinnGen.marg.503.2nd.p_aware"), "aware 2nd 503")
sv("S4.t3.FinnGen.cond.489.2nd.p_aware", "0.0497", "test (iii) 489 (conditional)", paste0(F, 360), v("S4.t3.FinnGen.cond.489.2nd.p_aware"), "aware 2nd 489 conditional")
sv("S4.t3.FinnGen.cond.503.1st.p_aware", "0.0378", "test (iii) first-order (conditional)", paste0(F, 360), v("S4.t3.FinnGen.cond.503.1st.p_aware"), "aware 1st 503 conditional")
sv("S4.t3.FinnGen.cond.503.2nd.p_cojo", "0.0597", "test (iii) beta-space 2nd", paste0(F, 361), v("S4.t3.FinnGen.cond.503.2nd.p_cojo"), "beta-space own-source frequencies, aware 2nd",
   expl_diff = sprintf("convention difference (allele-frequency source in the beta-space conditioning). Post-comparison diagnostic with 503 reference frequencies: %.4f (2nd) / %.4f (1st); own-source 1st-order %.4f. %s",
                       dv("diag.cojo_reffreq.t3.FinnGen.cond.p_aware_2nd"), dv("diag.cojo_reffreq.t3.FinnGen.cond.p_aware_1st"), dv("diag.cojo_ownfreq.t3.FinnGen.cond.p_aware_1st"),
                       if (abs(dv("diag.cojo_reffreq.t3.FinnGen.cond.p_aware_2nd") - 0.0597) < 0.00005) "Reference-frequency diagnostic reproduces the seen value" else
                         "The reference-frequency diagnostic comes close to the seen pair (0.0597 / 0.0494), so the difference is mostly the frequency-source convention; the residual (about 0.0002-0.0004) is unexplained without the seen script (not opened). Beta-space is printed beside and is not decisive; the combined-sentence branch uses the z-space p"))
sv("S4.t3.FinnGen.cond.503.2nd.p_cojo", "0.0494", "test (iii) beta-space 1st (compared with Stage 1 2nd-order beta-space; Stage 1 has no 1st-order beta-space id)", paste0(F, 361),
   v("S4.t3.FinnGen.cond.503.2nd.p_cojo"), "beta-space own-source frequencies, aware 2nd", expl_diff = "different SE order and frequency convention; see previous row")
S[[length(S)]]$quantity_id <- "S4.t3.FinnGen.cond.503.1st.p_cojo(seen only)"
sv("S4.t3.DIAMANTE.cond.503.2nd.p_aware", "0.030", "test (iii) DIAMANTE conditional", paste0(F, 362), v("S4.t3.DIAMANTE.cond.503.2nd.p_aware"), "aware 2nd 503")
sv("S4.ld.boot.t3.FinnGen.cond.aware2nd.share_lt05", "0.617", "test (iii) bootstrap share below 0.05, FinnGen (seen as 61.7%)", paste0(F, 363), v("S4.ld.boot.t3.FinnGen.cond.aware2nd.share_lt05"),
   "R MT seed 20261007, 5,000 draws", expl_diff = "bootstrap implementation/RNG differs from scoping (seen 'numpy'/'R' not stated for this share); Monte Carlo SE of a share near 0.62 with 5,000 draws is about 0.007. Same classification: < 0.90 (not robust)")
sv("S4.ld.boot.t3.DIAMANTE.cond.aware2nd.share_lt05", "1.00", "test (iii) bootstrap share, DIAMANTE (seen as 100%)", paste0(F, 363), v("S4.ld.boot.t3.DIAMANTE.cond.aware2nd.share_lt05"), "R MT seed 20261007")
sv("S4.t3.FinnGen.cond.503.2nd.p_ign", "0.043", "test (iii) ign conditional (R1)", paste0(F, 364), v("S4.t3.FinnGen.cond.503.2nd.p_ign"), "ign 2nd 503")
sv("S4.t3.FinnGen.marg.503.2nd.p_ign", "0.066", "test (iii) ign marginal (R1)", paste0(F, 364), v("S4.t3.FinnGen.marg.503.2nd.p_ign"), "ign 2nd 503")
sv("S4.ld.boot.FinnGen.cond.aware2nd.share_lt05", "0.838", "FinnGen cond aware bootstrap share below 0.05 (seen range 83.8-84.3%; lower end compared)", paste0(F, 370),
   v("S4.ld.boot.FinnGen.cond.aware2nd.share_lt05"), "R MT seed 20261007, 5,000 draws",
   expl_same = "inside the seen range 0.838-0.843. Scoping had an R run; its seed/call are not known to me (seen script not opened): if it was R Mersenne-Twister with seed 20261007 and the same sample.int call, agreement is expected (same RNG and seed)",
   tol = 0.0055)
sv("S4.ld.boot.FinnGen.cond.aware2nd.p_lo", "0.021", "FinnGen cond aware bootstrap 2.5th percentile p", paste0(F, 371), v("S4.ld.boot.FinnGen.cond.aware2nd.p_lo"), "R",
   expl_diff = "Monte Carlo/RNG difference between bootstrap implementations (scoping seed and percentile rule not known to me); Stage 1 value 0.02048 sits 0.00002 below the 0.0205 rounding boundary of the printed 0.021")
sv("S4.ld.boot.FinnGen.cond.aware2nd.p_hi", "0.070", "FinnGen cond aware bootstrap 97.5th percentile p", paste0(F, 371), v("S4.ld.boot.FinnGen.cond.aware2nd.p_hi"), "R")
sv("S4.ld.boot.DIAMANTE.cond.aware2nd.share_lt05", "1.00", "DIAMANTE conditional bootstrap share (seen 100%)", paste0(F, 372), v("S4.ld.boot.DIAMANTE.cond.aware2nd.share_lt05"), "R")
sv("S4.check.PM_marginal.PPH4.p1np", "0.996", "R1 P-M marginal PP.H4 p12 1e-5", paste0(F, 327), v("S4.check.PM_marginal.PPH4.p1np"), "coloc 5.2.3 on p1np (2,645 SNPs)", expl_same = W)
sv("S4.coloc_clause.arginine.A.503.PPH4", "0.981", "R1 P-M A|B PP.H4 p12 1e-5 (FinnGen-era R1 coloc)", paste0(F, 327), v("S4.coloc_clause.arginine.A.503.PPH4"), "conditional PP.H4 on p1np, 503",
   expl_diff = sprintf("SNP-set/setup difference. Post-comparison diagnostics: all 3,067 harmonised SNPs (palindromes kept) %.4f; p1np %.4f. %s Rule 3(iv) firing unaffected (>= 0.80 either way)",
                       dv("diag.coloc.all3067.PM.A_given_B.p12_1e-05"), dv("diag.coloc.p1np.PM.A_given_B.p12_1e-05"),
                       if (abs(dv("diag.coloc.all3067.PM.A_given_B.p12_1e-05") - 0.981) < 0.0005 && abs(dv("diag.coloc.all3067.PM.A_given_B.p12_1e-06") - 0.838) < 0.0005)
                         sprintf("EXPLAINED: the all-3,067-SNP diagnostic reproduces the seen R1 values at both priors (%.4f / %.4f vs 0.981 / 0.838), so R1 ran on all harmonised SNPs with palindromes kept; Stage 1 follows the pre-declaration's '503' = p1np (palindromes dropped).",
                                 dv("diag.coloc.all3067.PM.A_given_B.p12_1e-05"), dv("diag.coloc.all3067.PM.A_given_B.p12_1e-06")) else "Not reproduced by either diagnostic: unexplained residual."))
sv("S4.coloc_clause.arginine.B.503.PPH4", "1.000", "R1 P-M B|A PP.H4 p12 1e-5", paste0(F, 327), v("S4.coloc_clause.arginine.B.503.PPH4"), "conditional PP.H4 on p1np, 503")
sv("S4.coloc_clause.T2D.A.503.PPH4", "0.954", "R1 P-D A|B PP.H4 p12 1e-5", paste0(F, 328), v("S4.coloc_clause.T2D.A.503.PPH4"), "conditional PP.H4 on p1np, 503, T2D cc (W 0.2^2)",
   expl_diff = sprintf("setup difference. Post-comparison diagnostics (p12 1e-5): all 3,067 SNPs cc %.4f / prior sd 0.15 %.4f; p1np cc %.4f / prior sd 0.15 %.4f. %s Rule 3(iv) firing unaffected (>= 0.80 on both panels)",
                       dv("diag.coloc.all3067.PD.A_given_B.cc.p12_1e-05"), dv("diag.coloc.all3067.PD.A_given_B.quantW015.p12_1e-05"),
                       dv("diag.coloc.p1np.PD.A_given_B.cc.p12_1e-05"), dv("diag.coloc.p1np.PD.A_given_B.quantW015.p12_1e-05"),
                       if (min(abs(c(dv("diag.coloc.all3067.PD.A_given_B.cc.p12_1e-05"), dv("diag.coloc.all3067.PD.A_given_B.quantW015.p12_1e-05"), dv("diag.coloc.p1np.PD.A_given_B.quantW015.p12_1e-05")) - 0.954)) < 0.0005)
                         sprintf("EXPLAINED: the all-3,067-SNP cc diagnostic reproduces the seen R1 values at both priors (%.4f / %.4f vs 0.954 / 0.674): R1 used all harmonised SNPs (palindromes kept); the T2D prior (cc 0.2 vs 0.15) is immaterial here.", dv("diag.coloc.all3067.PD.A_given_B.cc.p12_1e-05"), dv("diag.coloc.all3067.PD.A_given_B.cc.p12_1e-06")) else "No diagnostic reproduces 0.954 exactly: unexplained residual."))
sv("S4.coloc_clause.T2D.B.503.PPH4", "0.009", "R1 P-D B|A PP.H4 p12 1e-5", paste0(F, 328), v("S4.coloc_clause.T2D.B.503.PPH4"), "conditional PP.H4 on p1np, 503")
nv("S4.coloc.PM.marginal.p12_1e-6", "0.961", "R1 P-M marginal at p12 1e-6", paste0(F, 327), sprintf("not a W4a quantity; diagnostic on p1np %.4f, all 3,067 %.4f", dv("diag.coloc.p1np.PM.marginal.p12_1e-06"), dv("diag.coloc.all3067.PM.marginal.p12_1e-06")))
nv("S4.coloc.PD.marginal", "0.932 / 0.579", "R1 P-D marginal p12 1e-5 / 1e-6", paste0(F, 328), sprintf("not a W4a quantity; diagnostic cc on p1np %.4f / %.4f; all 3,067 %.4f / %.4f",
   dv("diag.coloc.p1np.PD.marginal.cc.p12_1e-05"), dv("diag.coloc.p1np.PD.marginal.cc.p12_1e-06"), dv("diag.coloc.all3067.PD.marginal.cc.p12_1e-05"), dv("diag.coloc.all3067.PD.marginal.cc.p12_1e-06")))
nv("S4.coloc.MD", "marginal 0.930/0.569; A|B 0.962/0.717; B|A 0.009/0.001", "R1 M-D colocalization", paste0(F, 329), "not a W4a quantity (Rule 3(iv) needs trait x protein only); not computed")
nv("S4.R1_joint_GLS", "-0.1668 (SE 0.0302, p 3.4e-8)", "R1 joint two-signal GLS T2D per SD ARG1", paste0(F, 330), "a secondary of the original rule, not in W4a's brief; not computed")
nv("S4.ld.fisher_ci", "-0.131 to +0.043", "Fisher 95% CI of r(A,B)", paste0(F, 366), "not declared for W4a; not computed (B-allele copies on 503 = 2 x 503 x 0.021869 = 22.0, matching the seen '22 copies')")
nv("S4.ukbppp_insample_r2", "0.0044; r 0.49 signal 3 vs A", "UKB-PPP in-sample LD", paste0(F, 367), "not a W4a quantity; not computed")
sv("S4.coloc_clause.ornithine.A.503.PPH4", "0.345", "ornithine cond PP.H4 at A, 503 (two_signal.tsv)", paste0(F, 387), v("S4.coloc_clause.ornithine.A.503.PPH4"), "fresh R coloc, p1np 503")
sv("S4.coloc_clause.ornithine.A.489.PPH4", "0.263", "ornithine cond PP.H4 at A, 489", paste0(F, 387), v("S4.coloc_clause.ornithine.A.489.PPH4"), "fresh R coloc, p2np 489")
sv("S4.coloc_clause.ornithine.B.503.PPH4", "0.9999", "ornithine cond PP.H4 at B, 503", paste0(F, 387), v("S4.coloc_clause.ornithine.B.503.PPH4"), "fresh R coloc")
sv("S4.coloc_clause.ornithine.A.503.z", "-2.70", "ornithine conditional z at A, 503 (Rule 3(iv) seen)", paste0(F, 118), v("S4.coloc_clause.ornithine.A.503.z"), "z-space conditional")
sv("S4.coloc_clause.ornithine.A.489.z", "-2.77", "ornithine conditional z at A, 489", paste0(F, 118), v("S4.coloc_clause.ornithine.A.489.z"), "z-space conditional")
sv("S4.coloc_clause.T2D.B.503.z", "-0.24", "T2D conditional z at B (FinnGen, R1)", paste0(F, 120), v("S4.coloc_clause.T2D.B.503.z"), "z-space conditional")
sv("S4.coloc_clause.arginine.A.503.PPH4", "0.99", "arginine cond PP.H4 (Rule 3(iv) seen, 503)", paste0(F, 121), v("S4.coloc_clause.arginine.A.503.PPH4"), "fresh R coloc")
S[[length(S)]]$quantity_id <- "S4.coloc_clause.arginine.A.503.PPH4 (Rule 3(iv) seen line)"
sv("S4.orn.503.cond.kA", "0.1164", "ornithine k_A, 503 (two_signal.tsv: conditional z-ratio MA/PA)", paste0(F, 385), v("S4.orn.503.cond.kA"), "per-SD beta ratio, conditional (different scale)",
   expl_diff = sprintf("different quantity: seen k is the engine's conditional z-ratio z_M/z_P; Stage 1 k is ornithine units per SD ARG1. Post-comparison diagnostic z-ratio: k_A %.4f (503) / %.4f (489), k_B %.4f / %.4f",
                       dv("diag.orn.zratio.kA.503"), dv("diag.orn.zratio.kA.489"), dv("diag.orn.zratio.kB.503"), dv("diag.orn.zratio.kB.489")))
sv("S4.orn.503.cond.kB", "0.4145", "ornithine k_B, 503 (z-ratio)", paste0(F, 385), v("S4.orn.503.cond.kB"), "per-SD beta ratio, conditional", expl_diff = "different quantity (z-ratio vs beta ratio); see previous row")
sv("S4.orn.503.cond.p_aware_2nd", "5.3e-5", "ornithine p aware 503 (two_signal dev_sd-based)", paste0(F, 386), v("S4.orn.503.cond.p_aware_2nd"), "delta-method ratio test, conditional, aware 2nd",
   expl_diff = "different test statistic: seen is the engine's dev_sd z-space deviation (first-order, made aware); Stage 1 is the Rule 5 delta-method ratio difference. Primary Stage 1 cell (marginal aware 503) p 8.09e-5. Same branch (p < 0.05)")
sv("S4.orn.489.cond.p_aware_2nd", "5.7e-5", "ornithine p aware 489", paste0(F, 386), v("S4.orn.489.cond.p_aware_2nd"), "delta-method ratio test, conditional, aware 2nd", expl_diff = "different statistic; see previous row")
sv("S4.orn.503.cond.p_ign", "6.0e-5", "ornithine p ign 503", paste0(F, 386), v("S4.orn.503.cond.p_ign"), "delta-method, conditional, ign 2nd", expl_diff = "different statistic; see previous rows")
sv("S4.orn.489.cond.p_ign", "6.9e-5", "ornithine p ign 489", paste0(F, 386), v("S4.orn.489.cond.p_ign"), "delta-method, conditional, ign 2nd", expl_diff = "different statistic; see previous rows")
# amendment A2 seen values
chk_line(al, 38, "aware p 0.0365"); chk_line(al, 39, "aware p 0.0258")
sv("S4.grid.FinnGen.cond.p_aware_2nd", "0.0365", "A2 seen: FinnGen conditional aware", "full amendment l.38", v("S4.grid.FinnGen.cond.p_aware_2nd"), "aware 2nd")
S[[length(S)]]$quantity_id <- "S4.grid.FinnGen.cond.p_aware_2nd (A2 seen)"
sv("S4.grid.DIAMANTE.marg.p_aware_2nd", "0.0258", "A2 seen: DIAMANTE marginal aware", "full amendment l.39", v("S4.grid.DIAMANTE.marg.p_aware_2nd"), "aware 2nd")
S[[length(S)]]$quantity_id <- "S4.grid.DIAMANTE.marg.p_aware_2nd (A2 seen)"
S[[length(S) + 1]] <- data.frame(quantity_id = "A2.condition", seen_value = "would be met; both A-B differences negative", seen_label = "A2 seen [full amendment l.40-41]",
  stage1_value = "met", stage1_label = "branches.tsv Amendment 1 A2", abs_diff = "", explanation = "same classification")
S[[length(S) + 1]] <- data.frame(quantity_id = "Rule3ii.T2D.robust", seen_value = "FinnGen 84-85%; DIAMANTE 100%", seen_label = "Why 90% [full prespec l.95]",
  stage1_value = sprintf("FinnGen cond %.4f; DIAMANTE marg %.4f, cond %.4f", v("S4.ld.boot.FinnGen.cond.aware2nd.share_same"), v("S4.ld.boot.DIAMANTE.marg.aware2nd.share_same"), v("S4.ld.boot.DIAMANTE.cond.aware2nd.share_same")),
  stage1_label = "share with the full-panel conclusion", abs_diff = "", explanation = "same classification: FinnGen not robust (< 0.90), DIAMANTE robust")
S[[length(S) + 1]] <- data.frame(quantity_id = "Rule3iv.firing", seen_value = "ornithine A fires; T2D B does not; arginine does not", seen_label = "Rule 3(iv) [full prespec l.117-121]",
  stage1_value = "ornithine A fires; T2D B does not (|z| < 1.96); arginine A/B do not; T2D A does not; ornithine B does not", stage1_label = "branches.tsv Rule 3(iv)", abs_diff = "",
  explanation = "same firing")
CS <- rbindlist(S)
fwrite(CS, file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", quote = FALSE)
nd <- sum(grepl("^DIFFERS", CS$explanation)); nu <- sum(grepl("unexplained", CS$explanation))
logf("COMPARE_SEEN.tsv written: %d rows; %d differ beyond printed precision; %d carry 'unexplained'; md5 %s", nrow(CS), nd, nu, md5f(file.path(OUT, "COMPARE_SEEN.tsv")))
logf("results.tsv md5 after compare step (must equal the pre-compare freeze): %s", md5f(file.path(OUT, "results.tsv")))
logf("==== W4a compare-seen step end ====")
