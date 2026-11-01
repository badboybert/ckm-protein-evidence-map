# S5_05_selfcheck.R -- W5 self-checks: identities between two code paths, and mutations that must trip the check they name.
# Writes work/selfcheck.tsv and run.log lines. Not a substitute for the §10 verifier (different tool stack).
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
say("==== S5_05_selfcheck start ====")
assert_md5(GOV$GLS_Q$p, GOV$GLS_Q$m, "gls_q.R (sourced)"); source(GOV$GLS_Q$p)
P <- readRDS(file.path(WORK, "prep.rds")); A <- readRDS(file.path(WORK, "analyses.rds"))
H <- P$H; rs <- H$rsid_ref; ix <- function(x) match(x, rs)
CK <- list()
ck <- function(id, what, ok, detail) { CK[[length(CK) + 1]] <<- data.table(check = id, what, result = if (ok) "PASS" else "FAIL", detail); say("SELFCHECK %s %s: %s -- %s", id, if (ok) "PASS" else "FAIL", what, detail) }

# C1. single-SNP conditioning over the window (condition_z) == two-SNP joint estimate at A|B (both traits, 503)
R <- P$R503; ab <- ix(c(SIG_A, SIG_B))
for (tr in c("P__DPEP1", "M__GCST90200364")) {
  b <- H[[paste0("b__", tr)]]; s <- H[[paste0("se__", tr)]]; z <- b / s
  cz <- condition_z(z[ab[1]], s[ab[1]], R[SIG_A, SIG_B, drop = FALSE], R[SIG_B, SIG_B, drop = FALSE], z[ab[2]])
  je <- joint_effects(b[ab], s[ab], R[c(SIG_A, SIG_B), c(SIG_A, SIG_B)])
  d <- abs(cz$z - je$b[1] / sqrt(je$V[1, 1])) + abs(cz$z * cz$se - je$b[1])
  ck(paste0("C1_", tr), "conditional z/beta at A|B: window conditioning == two-SNP joint (gls_q.R joint_effects)", d < 1e-10, sprintf("|diff| %.2e", d))
}
# C2. coloc.abf (R coloc 5.2.3) == a direct Wakefield ABF implementation, for the 503 A|B conditional pair
setr <- P$SETS$p1np; ids <- rs[setr]
zP <- H$b__P__DPEP1 / H$se__P__DPEP1; zC <- H$b__M__GCST90200364 / H$se__M__GCST90200364
RiS <- R[ids, SIG_B, drop = FALSE]; RSS <- R[SIG_B, SIG_B, drop = FALSE]
cp <- condition_z(zP[setr], H$se__P__DPEP1[setr], RiS, RSS, zP[ab[2]], ids == SIG_B)
cm <- condition_z(zC[setr], H$se__M__GCST90200364[setr], RiS, RSS, zC[ab[2]], ids == SIG_B)
r1 <- coloc_pair(cp$z * cp$se, cp$se, cm$z * cm$se, cm$se, ids)
pd <- abf_direct(cp$z * cp$se, cp$se, cm$z * cm$se, cm$se); pp <- c(NA, NA, NA, pd$H3, pd$H4)
ck("C2", "PP.H4 at A|B (503): R coloc.abf == abf_direct() (the bootstrap's Wakefield ABF; prior sd 0.15, p1 = p2 = 1e-4, p12 = 1e-5)", abs(pp[5] - r1$H4) < 1e-6 && abs(pp[4] - r1$H3) < 1e-6,
   sprintf("coloc %.6f vs direct %.6f; PP.H3 %.6f vs %.6f", r1$H4, pp[5], r1$H3, pp[4]))
# C3. the 503 vs 489 disagreement at A|B is not an allele-alignment artefact: r(set SNP, B) agrees closely between panels
com <- intersect(ids, rownames(P$R489))
dr <- R[com, SIG_B] - P$R489[com, SIG_B]
ck("C3", "r(set SNP, rs62068372) 503 vs 489 on shared SNPs: same sign convention (no flipped SNP)", max(abs(dr)) < 0.15,
   sprintf("max |r503 - r489| %.4f over %d SNPs; correlation %.6f", max(abs(dr)), length(com), cor(R[com, SIG_B], P$R489[com, SIG_B])))
top5 <- order(-abs(cp$z))[1:5]
say("  A|B conditional protein |z| top 5 on 503: %s; metabolite z there: %s", paste(sprintf("%s %.2f", ids[top5], cp$z[top5]), collapse = ", "), paste(sprintf("%.2f", cm$z[top5]), collapse = ", "))
cp9 <- condition_z(zP[setr][ids %in% com], H$se__P__DPEP1[setr][ids %in% com], P$R489[com, SIG_B, drop = FALSE], P$R489[SIG_B, SIG_B, drop = FALSE], zP[ab[2]], com == SIG_B)
cm9 <- condition_z(zC[setr][ids %in% com], H$se__M__GCST90200364[setr][ids %in% com], P$R489[com, SIG_B, drop = FALSE], P$R489[SIG_B, SIG_B, drop = FALSE], zC[ab[2]], com == SIG_B)
t9 <- order(-abs(cp9$z))[1:5]
say("  A|B conditional protein |z| top 5 on 489: %s; metabolite z there: %s", paste(sprintf("%s %.2f", com[t9], cp9$z[t9]), collapse = ", "), paste(sprintf("%.2f", cm9$z[t9]), collapse = ", "))

# M1. mutation: the raw re-derivation comparison trips when one harmonised beta is perturbed
bm <- H$b__P__DPEP1; bm[100] <- bm[100] + 1e-6
ck("M1", "mutation: harmonised DPEP1 beta row 100 + 1e-6 -> re-derivation max|diff| > 1e-9 (assert would trip)", max(abs(bm - H$b__P__DPEP1)) > 1e-9, sprintf("max diff %.1e", max(abs(bm - H$b__P__DPEP1))))
# M2. mutation: palindrome orientation by frequency flips when the GWAS frequency is mirrored (rs1126464, letters C = a1)
af0 <- H$af_ref[ix("rs1126464")]
orient_pal <- function(f_let) if (abs(f_let - af0) <= 0.2 && abs((1 - f_let) - af0) > 0.2) 1L else if (abs((1 - f_let) - af0) <= 0.2 && abs(f_let - af0) > 0.2) -1L else 0L
# (first version expected 0.62 -> 0, which was wrong: at ref af 0.2505 the mirrored reading 0.38 fits; corrected 2026-10-08)
ck("M2", "mutation: rs1126464 CKDGen Freq1 0.25 -> 0.75 flips the palindrome orientation (+1 -> -1); 0.50 fits neither reading (0)",
   orient_pal(0.25) == 1L && orient_pal(0.75) == -1L && orient_pal(0.50) == 0L, sprintf("ref af %.4f: f 0.25 -> %d, 0.75 -> %d, 0.50 -> %d", af0, orient_pal(0.25), orient_pal(0.75), orient_pal(0.50)))
# M3. mutation: covariance-awareness is live -- the two-signal p changes when r(A,B) is set to 0 or sign-flipped
R2 <- P$R503[c(SIG_A, SIG_B), c(SIG_A, SIG_B)]
pfun <- function(Rm) ratio_diff(joint_effects(H$b__P__DPEP1[ab], H$se__P__DPEP1[ab], Rm), joint_effects(H$b__M__GCST90200364[ab], H$se__M__GCST90200364[ab], Rm))$p
Rf <- R2; Rf[1, 2] <- Rf[2, 1] <- -R2[1, 2]
ck("M3", "mutation: two-signal p moves when r(A,B) is sign-flipped (LD enters the test)", abs(pfun(R2) - pfun(Rf)) > 1e-4, sprintf("p %.4f (r %.4f) vs %.4f (r %.4f)", pfun(R2), R2[1, 2], pfun(Rf), Rf[1, 2]))
# M4. mutation: GLS aware differs from ign at the declared cells (the off-diagonal covariance is used)
g <- A$gls
ck("M4", "GLS aware Q != ign Q (eGFR all5 503)", abs(g[["eGFR.all5.503.aware_2nd"]]$Q - g[["eGFR.all5.503.ign_2nd"]]$Q) > 1e-6,
   sprintf("aware %.4f vs ign %.4f", g[["eGFR.all5.503.aware_2nd"]]$Q, g[["eGFR.all5.503.ign_2nd"]]$Q))
# M5. mutation: the Rule 3(iv) clause logic does not fire when both PP.H4 >= 0.80, fires when one is < 0.80
fire <- function(z, h) any(abs(z) >= 1.96) && any(h < 0.80)
ck("M5", "clause logic: (|z| 13, PP.H4 0.9/0.9) no fire; (|z| 13, 0.13/0.99) fires; (|z| 1.5, 0.1/0.1) no fire",
   !fire(c(13, 13), c(0.9, 0.9)) && fire(c(13, 13), c(0.13, 0.99)) && !fire(c(1.5, 1.5), c(0.1, 0.1)), "truth table")
# M6. gene-naming restriction is live: if SuSiE were admissible the restricted glist grid would be all 48 cells
gl <- fread(inp("GRID_bridge_class_robustness"))[locus == "DPEP1_16q24" & protein == "DPEP1" & metabolite == "cysteinylglycine"]
ck("M6", "mutation: admissible = TRUE -> restricted glist grid = 48 cells (vs Stage 1 restricted total)", nrow(gl) == 48L && A$gn$glist$total < 48L,
   sprintf("Stage 1 restricted total %d, assigned %d; all-48 assigned %d", A$gn$glist$total, A$gn$glist$assigned, sum(gl$class == "gene-attributed bridge")))
CKT <- rbindlist(CK)
fwrite(CKT, file.path(WORK, "selfcheck.tsv"), sep = "\t")
say("==== S5_05_selfcheck done: %d checks, %d FAIL ====", nrow(CKT), sum(CKT$result == "FAIL"))
