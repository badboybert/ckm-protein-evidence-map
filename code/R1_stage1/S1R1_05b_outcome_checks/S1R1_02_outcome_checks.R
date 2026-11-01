# External outcome checks, exactly as pre-specified in S1R1_PRESPEC_outcome_checks.md (md5 in S1R1_PRESPEC.md5).
# Inputs: instrument rows extracted by S1R1_01_extract.py (<D>/data/paper31_r1_stage1/extracts, md5 in *.done).
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
suppressPackageStartupMessages(library(data.table))
OD <- paste0(.P31R_V3, "/R1_stage1/S1R1_05b_outcome_checks")
EX <- paste0(.P31R_EXT, "/D/data/paper31_r1_stage1/extracts")
FG <- paste0(.P31R_EXT, "/D/data/paper31_r1_stage1/finngen_ldsc")
L <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); L <<- c(L, s) }
pre <- readLines(file.path(OD, "S1R1_PRESPEC.md5"))[1]
say("prespec md5 recorded %s; now %s", strsplit(pre, " ")[[1]][1], unname(tools::md5sum(file.path(OD, "S1R1_PRESPEC_outcome_checks.md"))))
stopifnot(strsplit(pre, " ")[[1]][1] == unname(tools::md5sum(file.path(OD, "S1R1_PRESPEC_outcome_checks.md"))))

ins <- fread(file.path(OD, "S1R1_eur_instruments.tsv"))
dis <- fread(file.path(OD, "S1R1_eur_discovery_mr.tsv"))
P <- data.table(
  pair = c("PAM-T2D", "RAD23B-T2D", "ARG1-T2D", "PCSK9-CAD", "TGFB1-CAD", "FGFR4-T2D (EAS)", "ARG1-T2D signal B"),
  gene = c("PAM", "RAD23B", "ARG1", "PCSK9", "TGFB1", "FGFR4", "ARG1"),
  rs   = c("rs149802978", "rs7035725", "rs2781668", "rs11591147", "rs73045269", "rs451643", "rs17788484"),
  res  = c("DIAMANTE", "DIAMANTE", "DIAMANTE", "CARDIOGRAM", "CARDIOGRAM", "TPMI", "DIAMANTE"))
# exposure (effect allele, beta, se, eaf) and discovery
exp <- ins[, .(gene, rs = SNP, ea = effect_allele, oa = other_allele, bx = beta, sx = se, eafx = eaf)]
exp <- rbind(exp, data.table(gene = "FGFR4", rs = "rs451643", ea = "T", oa = "G", bx = 0.98, sx = 0.023, eafx = 0.42),
                  data.table(gene = "ARG1", rs = "rs17788484", ea = "T", oa = "C", bx = 0.532254, sx = 0.026935, eafx = 0.0219))
P <- merge(P, exp, by = c("gene", "rs"), sort = FALSE)
dd <- dis[, .(gene, disc_b = b, disc_se = se, disc_p = pval)]
dd <- rbind(dd, data.table(gene = "FGFR4", disc_b = 0.0482, disc_se = 0.0069, disc_p = 3.888478978660555e-12))
P[, k := pair != "ARG1-T2D signal B"]
P <- merge(P, dd, by = "gene", all.x = TRUE, sort = FALSE)
P[pair == "ARG1-T2D signal B", `:=`(disc_b = -0.0193 / 0.532254, disc_se = 0.0387 / 0.532254, disc_p = 0.62)]  # FinnGen at B (triangle log)

# ---- check-resource rows, aligned to the exposure effect allele ----
di <- fread(file.path(EX, "DIAMANTE_instr.tsv"))
ca <- fread(file.path(EX, "CARDIOGRAM_instr.tsv"))
tp <- fread(file.path(EX, "TPMI_250.2_instr.tsv"))
pal <- function(a, b) paste0(sort(c(a, b)), collapse = "") %in% c("AT", "CG")
get_row <- function(r) {
  if (r$res == "DIAMANTE") { x <- di[rsID == r$rs]; if (nrow(x) != 1) return(NULL)
    list(ea = x$EA, oa = x$NEA, b = x$Beta, s = x$SE, eaf = x$EAF, p = x$Pvalue, n = NA, pos = paste0(x$CHR, ":", x$BP)) }
  else if (r$res == "CARDIOGRAM") { x <- ca[rsid == r$rs]; if (nrow(x) != 1) return(NULL)
    list(ea = x$alt, oa = x$ref, b = x$beta_alt, s = x$se, eaf = x$af_alt, p = 10^-x$lp, n = x$ss, pos = paste0(x$chr, ":", x$pos)) }
  else { x <- tp; if (nrow(x) != 1) return(NULL)
    list(ea = x$Allele2, oa = x$Allele1, b = x$BETA, s = x$SE, eaf = x$AF_Allele2, p = x$p.value, n = x$N_case + x$N_ctrl,
         pos = paste0(x$CHR, ":", x$POS), ncase = x$N_case, nctrl = x$N_ctrl) }
}
out <- list()
for (i in seq_len(nrow(P))) {
  r <- P[i]; y <- get_row(r)
  if (is.null(y)) { out[[i]] <- data.table(pair = r$pair, class = "not evaluable (instrument absent)"); next }
  s <- if (toupper(y$ea) == r$ea && toupper(y$oa) == r$oa) 1 else if (toupper(y$ea) == r$oa && toupper(y$oa) == r$ea) -1 else NA
  if (is.na(s)) { out[[i]] <- data.table(pair = r$pair, class = "not evaluable (alleles do not match)"); next }
  by <- s * y$b; eafy <- if (s == 1) y$eaf else 1 - y$eaf
  pflag <- if (pal(r$ea, r$oa)) { maf <- min(r$eafx, 1 - r$eafx)
    if (maf > 0.42) "palindromic, AF-uninformative (forward-strand assumption)" else if (abs(eafy - r$eafx) > 0.10) "palindromic, AF DISCORDANT (possible strand flip)" else "palindromic, AF concordant" } else ""
  wr <- by / r$bx; se1 <- y$s / abs(r$bx); se2 <- abs(wr) * sqrt((y$s / by)^2 + (r$sx / r$bx)^2)
  p <- 2 * pnorm(-abs(wr / se1))
  same <- sign(wr) == sign(r$disc_b)
  cls <- if (same && p < 0.05) "concordant" else if (same) "directionally concordant, not significant" else if (p < 0.05) "discordant" else "inconclusive"
  hz <- (wr - r$disc_b) / sqrt(se1^2 + r$disc_se^2)
  out[[i]] <- data.table(pair = r$pair, resource = r$res, instrument = r$rs, pos_in_resource = y$pos, ea = r$ea,
    eaf_exposure = r$eafx, eaf_resource_aligned = eafy, palindrome = pflag, b_out = by, se_out = y$s, p_out_file = y$p,
    wald = wr, wald_se = se1, wald_se_delta = se2, p = p, discovery_wald = r$disc_b, discovery_se = r$disc_se,
    het_z = hz, het_p = 2 * pnorm(-abs(hz)), class = cls, bonferroni_6 = p < 0.05 / 6,
    n_resource = if (!is.null(y$n)) y$n else NA)
}
R <- rbindlist(out, fill = TRUE)
fwrite(R, file.path(OD, "S1R1_outcome_checks.tsv"), sep = "\t")
say("== Wald ratios (per SD protein; log-odds) in the check resources, classified as pre-specified")
for (i in seq_len(nrow(R))) with(R[i], say("  %-20s %-10s %s EA %s eaf %.3f/%.3f %s | b_out %.4f (SE %.4f) -> Wald %.4f (SE %.4f), p %.2g | discovery %.4f (SE %.4f) | het p %.3f | %s%s",
  pair, resource, instrument, ea, eaf_exposure, eaf_resource_aligned, palindrome, b_out, se_out, wald, wald_se, p, discovery_wald, discovery_se, het_p, class,
  if (isTRUE(bonferroni_6)) " [p < 0.05/6]" else ""))
if ("TPMI" %in% R$resource) say("TPMI N from file: %s", paste(tp$N_case, tp$N_ctrl))

# ---- ARG1 triangle in DIAMANTE (signals A and B), rule of S1R1_ARG1_PRESTATED_RULE.md applied to DIAMANTE ----
a <- R[pair == "ARG1-T2D"]; b <- R[pair == "ARG1-T2D signal B"]
rA <- c(r = a$wald, se = a$wald_se_delta); rB <- c(r = b$wald, se = b$wald_se_delta)
hz <- (rA[["r"]] - rB[["r"]]) / sqrt(rA[["se"]]^2 + rB[["se"]]^2); hp <- 2 * pnorm(-abs(hz))
rAB <- -0.0445
# conditional (reference LD, equal-N approximation) for the T2D z at each signal
zA <- a$b_out / a$se_out; zB <- b$b_out / b$se_out
zAc <- (zA - rAB * zB) / sqrt(1 - rAB^2); zBc <- (zB - rAB * zA) / sqrt(1 - rAB^2)
exp_zB <- rA[["r"]] * 0.532254 / b$se_out
same <- sign(rA[["r"]]) == sign(rB[["r"]])
v <- if (hp < 0.05 || (!same && b$p < 0.05)) "WEAKENS" else if (same && b$p < 0.05) "HOLDS" else if (same) "CONSISTENT, NOT INDEPENDENTLY CONFIRMED" else "INDETERMINATE"
say("== ARG1 in DIAMANTE: A %.4f (SE %.4f), B %.4f (SE %.4f, p %.2g); heterogeneity z %.2f p %.3f; T2D z marginal A %.2f B %.2f, conditional A %.2f B %.2f; expected z at B under A's ratio %.2f",
    rA[["r"]], rA[["se"]], rB[["r"]], rB[["se"]], b$p, hz, hp, zA, zB, zAc, zBc, exp_zB)
say("ARG1 DIAMANTE verdict (same rule): %s   [FinnGen verdict, pre-stated primary: WEAKENS]", v)

# ---- measured overlap with FinnGen: cross-trait LDSC intercepts ----
jk_wls <- function(y, l, w, blk) {
  X1 <- w; Xl <- w * l; Xll <- w * l * l; Y1 <- w * y; Yl <- w * l * y
  S <- rowsum(cbind(X1, Xl, Xll, Y1, Yl), blk); tot <- colSums(S)
  solve2 <- function(s) as.vector(solve(matrix(c(s[1], s[2], s[2], s[3]), 2), c(s[4], s[5])))
  full <- solve2(tot); th <- t(apply(S, 1, function(r) solve2(tot - r))); K <- nrow(S)
  ps <- K * matrix(full, K, 2, byrow = TRUE) - (K - 1) * th
  list(est = full, se = sqrt(apply(ps, 2, var) / K))
}
fit <- function(y, l, b0, wf) { b <- b0; for (it in 1:2) { w <- wf(b, l); X <- cbind(1, l); b <- as.vector(solve(crossprod(X, X * w), crossprod(X, y * w))) }; wf(b, l) }
ld <- fread(paste0(.P31R_EXT, "/D/data/ldscore/LDscores_filtered.csv"), select = c("chr", "pos", "rs", "LDSC"))
ld <- ld[!rs %in% rs[duplicated(rs)]]; setorder(ld, chr, pos); ld[, blk := cut(seq_len(.N), 200, labels = FALSE)]
rdF <- function(tag) { d <- fread(file.path(FG, paste0(tag, "_ldsc.tsv.gz")), header = FALSE, sep = "\t",
  col.names = c("rs", "chrom", "pos", "ref", "alt", "beta", "se", "af", "pval")); d <- d[!rs %in% rs[duplicated(rs)]]
  d[pmin(af, 1 - af) >= 0.01, .(rs, ref, alt, z1 = beta / se)] }
amb <- function(a, b) paste0(pmin(a, b), pmax(a, b)) %in% c("AT", "CG")
pairld <- function(fg, ext, label) {
  W <- merge(fg, ext, by = "rs"); W <- W[!amb(ref, alt)]
  W[, s := fifelse(ea == alt & oa == ref, 1, fifelse(ea == ref & oa == alt, -1, NA_real_))]
  W <- W[!is.na(s) & is.finite(z2)]; W[, z2 := s * z2]
  W <- merge(W, ld[, .(rs, LDSC, blk)], by = "rs")
  k <- W$z1^2 < 80 & W$z2^2 < 80; W <- W[k]
  uni <- function(z) { y <- z^2; wf <- function(b, l) 1 / (pmax(l, 1) * 2 * pmax(b[1] + b[2] * l, 1)^2)
    w <- fit(y, W$LDSC, c(1, mean(y - 1) / mean(W$LDSC)), wf); j <- jk_wls(y, W$LDSC, w, W$blk); c(j$est, j$se) }
  u1 <- uni(W$z1); u2 <- uni(W$z2)
  y <- W$z1 * W$z2
  wf <- function(b, l) 1 / (pmax(l, 1) * (pmax(u1[1] + u1[2] * l, 1) * pmax(u2[1] + u2[2] * l, 1) + (b[1] + b[2] * l)^2))
  w <- fit(y, W$LDSC, c(0, mean(y) / mean(W$LDSC)), wf); j <- jk_wls(y, W$LDSC, w, W$blk)
  data.table(pair = label, n_snps = nrow(W), int_finngen = u1[1], int_resource = u2[1], gcov_int = j$est[1], gcov_int_se = j$se[1],
             null_z_corr = j$est[1] / sqrt(u1[1] * u2[1]), null_z_corr_se = j$se[1] / sqrt(u1[1] * u2[1]), gcov_slope = j$est[2])
}
dl <- fread(file.path(EX, "DIAMANTE_ldsc.tsv"))[!rsID %in% rsID[duplicated(rsID)] & SE > 0, .(rs = rsID, ea = EA, oa = NEA, z2 = Beta / SE)]
cl <- fread(file.path(EX, "CARDIOGRAM_ldsc.tsv"))[!rsid %in% rsid[duplicated(rsid)] & se > 0, .(rs = rsid, ea = alt, oa = ref, z2 = beta_alt / se)]
OV <- rbind(pairld(rdF("T2D"), dl, "FinnGen R10 T2D x DIAMANTE 2018 noUKB"), pairld(rdF("CHD"), cl, "FinnGen R10 I9_CHD x CARDIoGRAMplusC4D 2015"))
fwrite(OV, file.path(OD, "S1R1_overlap_ldsc.tsv"), sep = "\t")
say("== measured overlap (cross-trait LDSC intercept; MAF >= 1%%, chi2 < 80, 200-block jackknife)")
for (i in seq_len(nrow(OV))) with(OV[i], say("  %-46s %d SNPs: univariate intercepts %.3f / %.3f; gcov_int %.4f (SE %.4f) -> null-z correlation %.4f (SE %.4f)",
  pair, n_snps, int_finngen, int_resource, gcov_int, gcov_int_se, null_z_corr, null_z_corr_se))
writeLines(L, file.path(OD, "S1R1_outcome_checks.log"))
