#!/usr/bin/env Rscript
# w4a_main.R -- Paper 3.1 R2 Stage 1, implementer W4a: ARG1 grid, LD uncertainty, test (iii), ornithine
# (pre-declaration rev 4.2 section 4 (i)(ii)(iii)(vi); amendment 1 A2; Rules 3(ii), 3(iv), 4, 5).
# Computes from the REDACTED pre-declaration only. Writes results.tsv, branches.tsv, run.log and work/ under this folder.
# COMPARE_SEEN.tsv is written afterwards by w4a_compare_seen.R (after the full text is opened).
#   Rscript w4a_main.R
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
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(Rsamtools); library(GenomicRanges); library(coloc)
})
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04a_arg1grid")
WORK <- file.path(OUT, "work")
dir.create(WORK, showWarnings = FALSE)
source(file.path(OUT, "w4a_funcs.R"))
LOGP <- file.path(OUT, "run.log")
if (file.exists(LOGP)) file.rename(LOGP, file.path(WORK, sprintf("run.log.prev_%s", format(Sys.time(), "%Y%m%d_%H%M%S"))))
log_init(LOGP)
T0 <- Sys.time()
logf("==== W4a (S1R1_R2_04a_arg1grid) run start; R %s; pid %d ====", R.version.string, Sys.getpid())
SEED <- 20261007L
NBOOT <- 5000L
RSCRIPT <- file.path(R.home("bin"), "Rscript.exe")

# ===================== S0. governing documents =====================
RESDIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
GOV <- data.frame(
  label = c("briefs", "prespec_redacted", "amend1_redacted", "prespec_full", "amend1_full", "gls_q.R", "RESOURCES.json",
            "outcome_check_prespec", "R1_ARG1_rule"),
  path = c(file.path(RESDIR, "STAGE1_BRIEFS.md"), file.path(RESDIR, "P31R_R2_STAGE1_PRESPEC_REDACTED.md"),
           file.path(RESDIR, "P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"),
           paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"),
           paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"),
           file.path(RESDIR, "gls_q.R"), file.path(RESDIR, "RESOURCES.json"),
           paste0(.P31R_V3, "/R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"),
           paste0(.P31R_V3, "/R1_stage1/S1R1_01_arg1_triangle/S1R1_ARG1_PRESTATED_RULE.md")),
  md5 = c("40d0d52f7055542d774c03f1edf04e00", "a281c3464619b1b67c28e13e7392c41c", "f265d600973c43239bb8377e0cd39d33",
          "3319dc39740d9f5b3b8822133d9108b7", "377508832ea9bd384133e869c456c909", "b7390fdf5b6f2f9c1d30c298d6afdafc",
          "6d21238155a62fb58f71e05a01638dd0", "c9079f2de1ffe75ce31f8a1d597d1fda", "4d2c96a0998563785dadfec1ef9df218"),
  stringsAsFactors = FALSE)
logf("-- S0 governing documents (md5 only for the full texts; they are not read by this script)")
for (i in seq_len(nrow(GOV))) assert_md5(GOV$path[i], GOV$md5[i], GOV$label[i])
logf("tool versions: R %s; data.table %s; coloc %s; Rsamtools %s; jsonlite %s; RNGkind %s",
     R.version.string, packageVersion("data.table"), packageVersion("coloc"), packageVersion("Rsamtools"),
     packageVersion("jsonlite"), paste(RNGkind(), collapse = "/"))
gst <- suppressWarnings(system2(RSCRIPT, c(shQuote(GOV$path[6], type = "cmd"), "--self-test"), stdout = TRUE, stderr = TRUE))
for (l in gst) logf("gls_q.R --self-test | %s", l)
if (!any(grepl("SELF-TEST PASS", gst))) stop("STOP: gls_q.R self-test did not pass")
source(GOV$path[6])   # shared GLS (Rule 6); W4a has no GLS Q, but joint_effects() cross-checks the z-space conditioning below

# ===================== S0b. staging snapshot (start) =====================
STAGING <- paste0(.P31R_V3, "/V4_local/staging/l3-finemap")
snap <- function(tag) {
  f <- list.files(STAGING, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  d <- data.frame(file = f, bytes = file.size(f), md5 = unname(tools::md5sum(f)), stringsAsFactors = FALSE)
  d <- d[order(d$file), ]
  fwrite(d, file.path(WORK, sprintf("staging_snapshot_%s.tsv", tag)), sep = "\t")
  logf("staging snapshot %s: %d files; digest of listing %s", tag, nrow(d), digest_listing(d))
  for (i in seq_len(nrow(d))) cat(sprintf("  [snapshot %s] %s  %s  %s\n", tag, d$md5[i], d$bytes[i], d$file[i]), file = LOGP, append = TRUE)
  d
}
digest_listing <- function(d) { tf <- tempfile(); writeLines(paste(d$file, d$md5), tf); m <- md5f(tf); unlink(tf); m }
SNAP0 <- snap("start")

# ===================== S0c. self-tests of my estimators =====================
logf("-- S0c estimator self-tests")
local({
  # (1) z-space conditioning == gls_q.R joint_effects for two SNPs
  bA <- 0.12; sA <- 0.01; bB <- -0.30; sB <- 0.02; r <- -0.3
  cz <- cond_z(bA, sA, bB, sB, r)
  je <- joint_effects(c(bA, bB), c(sA, sB), matrix(c(1, r, r, 1), 2))
  dv <- max(abs(c(cz$bA, cz$bB) - je$b), abs(c(cz$sA^2, cz$sB^2) - diag(je$V)), abs(cz$c * cz$sA * cz$sB - je$V[1, 2]))
  logf("self-test: cond_z vs gls_q.R joint_effects max abs diff %.3g (must be < 1e-12)", dv)
  if (dv > 1e-12) stop("self-test failed: conditioning")
  # (2) null size of the aware A-vs-B test for conditional estimates under strong LD; ign must be miscalibrated
  set.seed(1)
  r <- 0.6; nsim <- 4000
  Lc <- chol(matrix(c(1, r, r, 1), 2))
  sim <- function() {
    eP <- matrix(rnorm(2 * nsim), nsim) %*% Lc; eY <- matrix(rnorm(2 * nsim), nsim) %*% Lc
    bPt <- c(0.25, 0.5); sP <- c(0.01, 0.01); sY <- c(0.02, 0.02); th <- -0.2
    mP <- as.vector(matrix(c(1, r, r, 1), 2) %*% bPt)            # marginal = R x joint
    bPm <- sweep(eP, 2, sP, "*") + matrix(mP, nsim, 2, byrow = TRUE)
    bYm <- sweep(eY, 2, sY, "*") + matrix(th * mP, nsim, 2, byrow = TRUE)
    list(bPm = bPm, bYm = bYm, sP = sP, sY = sY)
  }
  s <- sim()
  X <- cond_z(s$bPm[, 1], s$sP[1], s$bPm[, 2], s$sP[2], r); Y <- cond_z(s$bYm[, 1], s$sY[1], s$bYm[, 2], s$sY[2], r)
  pa <- ab_test(Y, X, "2nd", TRUE)$p; pi <- ab_test(Y, X, "2nd", FALSE)$p
  logf("self-test: conditional, r = 0.6: size aware %.3f (expect ~0.05), size ign %.3f (must differ from aware)", mean(pa < 0.05), mean(pi < 0.05))
  if (abs(mean(pa < 0.05) - 0.05) > 0.015) stop("self-test failed: aware size")
  if (abs(mean(pa < 0.05) - mean(pi < 0.05)) < 0.01) stop("self-test failed: ign indistinguishable from aware")
  Xm <- marg_est(s$bPm[, 1], s$sP[1], s$bPm[, 2], s$sP[2], r); Ym <- marg_est(s$bYm[, 1], s$sY[1], s$bYm[, 2], s$sY[2], r)
  pm <- ab_test(Ym, Xm, "2nd", TRUE)$p
  logf("self-test: marginal (ratios equal under the null because marginal = R x joint is proportional), size aware %.3f", mean(pm < 0.05))
  if (abs(mean(pm < 0.05) - 0.05) > 0.015) stop("self-test failed: marginal aware size")
  # (3) mutation: a non-proportional B must be detected
  bYm2 <- s$bYm; bYm2[, 2] <- bYm2[, 2] + 0.08
  Y2 <- cond_z(bYm2[, 1], s$sY[1], bYm2[, 2], s$sY[2], r)
  pw <- mean(ab_test(Y2, X, "2nd", TRUE)$p < 0.05)
  logf("self-test: power against a +0.08 shift at B %.3f (must exceed 0.9)", pw)
  if (pw < 0.9) stop("self-test failed: power")
  # (4) Fieller CI excludes 1 iff the difference test rejects
  a <- ab_test(Y2, X, "2nd", TRUE)
  agree <- vapply(1:300, function(i) { f <- fieller(a$RA[i], a$RB[i], a$vA[i], a$vB[i], a$cov[i])
    if (is.na(f$lo)) return(NA); ((f$lo > 1) | (f$hi < 1)) == (a$p[i] < 0.05) }, logical(1))
  logf("self-test: Fieller-CI excludes 1 <=> p < 0.05 in %d of %d bounded cases", sum(agree, na.rm = TRUE), sum(!is.na(agree)))
  if (any(!agree, na.rm = TRUE)) stop("self-test failed: Fieller coherence")
})

# ===================== S1. resource map =====================
RES <- fromJSON(GOV$path[7], simplifyVector = FALSE)
names(RES) <- vapply(RES, function(e) e$id, "")
logf("RESOURCES.json: %d entries", length(RES))
INPUTS <- list()
inp <- function(id) {
  e <- RES[[id]]
  if (is.null(e)) stop(sprintf("STOP: resource id %s not in RESOURCES.json", id))
  if (e$role %in% c("seen_output", "seen_script", "full_text_after_results")) stop(sprintf("STOP: %s has role %s", id, e$role))
  assert_md5(e$path, e$md5, sprintf("%s (%s)", id, e$role))
  INPUTS[[id]] <<- c(id = id, path = e$path, md5 = e$md5, role = e$role)
  e$path
}
inp_extra <- function(label, path, md5, why) {
  assert_md5(path, md5, sprintf("%s [not a RESOURCES id: %s]", label, why))
  INPUTS[[label]] <<- c(id = label, path = path, md5 = md5, role = paste("extra:", why))
  path
}
nums_from <- function(s, pat) as.numeric(gsub(",", "", regmatches(s, regexec(pat, s))[[1]][-1]))

# ===================== results / branches collectors =====================
RR <- list(); BR <- list()
g7 <- function(x) { if (is.null(x) || length(x) == 0) return(""); if (is.character(x)) return(x); ifelse(is.na(x), "NA", sprintf("%.7g", x)) }
addr <- function(qid, ref, estimand, resource, panel, cov = "", se_order = "", value = NA, se = NA, lo = NA, hi = NA, p = NA,
                 n = NA, unit = "", notes = "") {
  RR[[length(RR) + 1]] <<- data.frame(quantity_id = qid, prespec_ref = ref, estimand = estimand, resource = resource, panel = panel,
    cov = cov, se_order = se_order, value = g7(value), se = g7(se), ci_lo = g7(lo), ci_hi = g7(hi), p = g7(p), n = g7(n),
    unit = unit, notes = notes, stringsAsFactors = FALSE)
}
addb <- function(rule_ref, inputs, branch, filled, qualifiers) {
  BR[[length(BR) + 1]] <<- data.frame(rule_ref = rule_ref, inputs = inputs, branch = branch, filled_template = filled,
                                      qualifiers = qualifiers, stringsAsFactors = FALSE)
}
assert_target <- function(id, computed, target, tol, derivation, ref) {
  ok <- is.finite(computed) && abs(computed - target) <= tol
  logf("ASSERT %s: computed %s vs target %s (tol %g) -> %s | derivation: %s", id, g7(computed), g7(target), tol,
       if (ok) "PASS" else "ASSERT_FAILED", derivation)
  addr(paste0("ASSERT.", id), ref, "assertion target", "", "", value = computed,
       notes = sprintf("%s: target %s, tol %g; derivation: %s", if (ok) "ASSERT_PASS" else "ASSERT_FAILED", g7(target), tol, derivation))
  ok
}

# ===================== S2. harmonised table + my own harmonisation =====================
logf("-- S2 harmonisation")
H <- read.delim(gzfile(inp("R1_ARG1_HARMONISED")), stringsAsFactors = FALSE,
                colClasses = c(snp = "character", rsid_ref = "character", rsid_trait = "character", chr = "character",
                               a1 = "character", a2 = "character"))
if (nrow(H) != 3067) stop(sprintf("STOP: harmonised rows %d != 3067", nrow(H)))
logf("R1_ARG1_HARMONISED: %d rows, %d palindromic; null rates: %s", nrow(H), sum(H$palindromic == 1),
     paste(sprintf("%s=%.4f", names(H), colMeans(is.na(H) | H == "")), collapse = " "))
setkey_al <- function(pos, x, y) paste(pos, pmin(x, y), pmax(x, y), sep = ":")
H$key37 <- setkey_al(H$pos37, H$a1, H$a2); H$key38 <- setkey_al(H$pos38, H$a1, H$a2)
if (anyDuplicated(H$key37) || anyDuplicated(H$key38)) stop("STOP: duplicate harmonised keys")
sig <- read.delim(inp("LOC_ARG1_6q23_signals"), stringsAsFactors = FALSE, colClasses = c(snp = "character", rsid = "character"))
RS_A <- "rs2781668"; RS_B <- "rs17788484"     # original rule: A = ARG1 protein signal 1, B = ARG1 protein signal 2
sA <- sig[sig$tid == "P__ARG1" & sig$signal == 1, ]; sB <- sig[sig$tid == "P__ARG1" & sig$signal == 2, ]
if (sA$rsid != RS_A || sB$rsid != RS_B) stop("STOP: signals.tsv protein signals 1/2 differ from the rule's A/B")
iA <- which(H$rsid_ref == RS_A); iB <- which(H$rsid_ref == RS_B)
if (length(iA) != 1 || length(iB) != 1) stop("STOP: A/B not unique in harmonised")
logf("A = %s (%s, pos37 %d, pos38 %d, a1 %s); B = %s (%s, pos37 %d, pos38 %d, a1 %s)", RS_A, H$snp[iA], H$pos37[iA], H$pos38[iA],
     H$a1[iA], RS_B, H$snp[iB], H$pos37[iB], H$pos38[iB], H$a1[iB])

join_report <- function(label, src_keys, n_src) {
  m <- match(H$key_use, src_keys)
  logf("JOIN %s: source rows %d; harmonised matched %d/%d (%.2f%%); p1np matched %d/%d; A %s B %s", label, n_src,
       sum(!is.na(m)), nrow(H), 100 * mean(!is.na(m)), sum(!is.na(m) & H$palindromic == 0), sum(H$palindromic == 0),
       if (is.na(m[iA])) "MISSING" else "ok", if (is.na(m[iB])) "MISSING" else "ok")
  m
}
TR <- data.frame(row = seq_len(nrow(H)))

# --- UKB-PPP ARG1 (ID = GRCh37 chr:pos:ref:alt; pos = GRCh38; effect = ALLELE1)
up <- read.delim(inp("EXP_UKBPPP_EUR_ARG1_OID21445"), stringsAsFactors = FALSE,
                 colClasses = c(ID = "character", ALLELE0 = "character", ALLELE1 = "character"))
if (nrow(up) != 11311) stop("STOP: UKB-PPP ARG1 rows != 11311")
logf("UKB-PPP ARG1 null rates: %s", paste(sprintf("%s=%.4f", names(up), colMeans(is.na(up) | up == "")), collapse = " "))
idp <- do.call(rbind, strsplit(up$ID, ":"))
up$pos37 <- as.integer(idp[, 2])
if (any(idp[, 1] != "6")) stop("STOP: UKB-PPP non-chr6 rows")
k <- setkey_al(up$pos37, up$ALLELE0, up$ALLELE1)
dupk <- k[duplicated(k)]
H$key_use <- H$key37
m <- join_report("UKB-PPP ARG1 (pos37+allele set)", ifelse(k %in% dupk, NA, k), nrow(up))
if (sum(!is.na(m)) != nrow(H)) stop("STOP: UKB-PPP join incomplete for harmonised SNPs")
fl <- ifelse(up$ALLELE1[m] == H$a1, 1, ifelse(up$ALLELE0[m] == H$a1, -1, NA))
TR$b_P <- fl * up$BETA[m]; TR$s_P <- up$SE[m]; TR$f_P <- ifelse(fl == 1, up$A1FREQ[m], 1 - up$A1FREQ[m])
logf("UKB-PPP pos (GRCh38) == harmonised pos38 on %d/%d; max|b - R1 b| %.3g; max|se - R1 se| %.3g", sum(up$pos[m] == H$pos38),
     nrow(H), max(abs(TR$b_P - H$b__P__ARG1)), max(abs(TR$s_P - H$se__P__ARG1)))
if (max(abs(TR$b_P - H$b__P__ARG1)) > 1e-9 || max(abs(TR$s_P - H$se__P__ARG1)) > 1e-9) stop("STOP: own UKB-PPP harmonisation != R1")
N_P <- nums_from(RES[["EXP_UKBPPP_EUR_ARG1_OID21445"]]$n_info, "samplesize ([0-9,]+)")
logf("N for ARG1 (UKB-PPP): %d (RESOURCES n_info 'samplesize at the instrument'; no N column in the window)", N_P)

# --- CLSA arginine / ornithine (GRCh38 harmonised GWAS Catalog windows)
read_clsa <- function(id, lab, rcol) {
  x <- read.delim(gzfile(inp(id)), stringsAsFactors = FALSE,
                  colClasses = c(effect_allele = "character", other_allele = "character", variant_id = "character", rsid = "character"))
  logf("CLSA %s: %d rows; null rates: %s", lab, nrow(x), paste(sprintf("%s=%.4f", names(x), colMeans(is.na(x) | x == "")), collapse = " "))
  if (any(x$chromosome != 6)) stop("STOP: CLSA non-chr6")
  k <- setkey_al(x$base_pair_location, x$effect_allele, x$other_allele); dk <- k[duplicated(k)]
  H$key_use <<- H$key38
  m <- join_report(sprintf("CLSA %s (pos38+allele set)", lab), ifelse(k %in% dk, NA, k), nrow(x))
  fl <- ifelse(x$effect_allele[m] == H$a1, 1, ifelse(x$other_allele[m] == H$a1, -1, NA))
  b <- fl * x$beta[m]; s <- x$standard_error[m]; f <- ifelse(fl == 1, x$effect_allele_frequency[m], 1 - x$effect_allele_frequency[m])
  ok <- !is.na(m)
  d1 <- max(abs(b[ok] - H[[paste0("b__", rcol)]][ok])); d2 <- max(abs(s[ok] - H[[paste0("se__", rcol)]][ok]))
  logf("CLSA %s: own vs R1 on %d matched: max|db| %.3g max|dse| %.3g; R1 rows with no own match: %d", lab, sum(ok), d1, d2, sum(!ok))
  if (d1 > 1e-9 || d2 > 1e-9 || sum(!ok) > 0) stop(sprintf("STOP: own CLSA %s harmonisation != R1", lab))
  list(b = b, s = s, f = f, p = x$p_value[m])
}
cl <- read_clsa("EXP_CLSA_arginine_GCST90200372", "arginine", "M__GCST90200372")
TR$b_arg <- cl$b; TR$s_arg <- cl$s; TR$f_arg <- cl$f; TR$pfile_arg <- cl$p
cl <- read_clsa("EXP_CLSA_ornithine_GCST90200404", "ornithine", "M__GCST90200404")
TR$b_orn <- cl$b; TR$s_orn <- cl$s; TR$f_orn <- cl$f; TR$pfile_orn <- cl$p
N_ARG <- nums_from(RES[["EXP_CLSA_arginine_GCST90200372"]]$n_info, "N ([0-9,]+)")
N_ORN <- nums_from(RES[["EXP_CLSA_ornithine_GCST90200404"]]$n_info, "N ([0-9,]+)")
logf("N for CLSA metabolites: arginine %d, ornithine %d (RESOURCES n_info; the windows have no N column)", N_ARG, N_ORN)

# --- FinnGen R10 T2D (GRCh38; tabix the window)
fgp <- inp("OUT_FG_R10_T2D")
inp_extra("OUT_FG_R10_T2D.tbi", paste0(fgp, ".tbi"), "7c7965b8e6794d20e897be42a8666895", "tabix index named in OUT_FG_R10_T2D notes")
man <- read.delim(inp("OUT_FG_R10_MANIFEST"), stringsAsFactors = FALSE, quote = "")
mt <- man[man$phenocode == "T2D", ]
if (nrow(mt) != 1) stop("STOP: manifest T2D row")
FG_CASE <- mt$num_cases; FG_CTRL <- mt$num_controls; FG_N <- FG_CASE + FG_CTRL
FG_NEFF <- 4 * FG_CASE * FG_CTRL / (FG_CASE + FG_CTRL)
nres <- nums_from(RES[["OUT_FG_R10_T2D"]]$n_info, "cases ([0-9]+) / controls ([0-9]+)")
logf("FinnGen T2D: cases %d controls %d N %d N_eff %.1f (manifest); RESOURCES n_info cases %d controls %d -> %s", FG_CASE, FG_CTRL,
     FG_N, FG_NEFF, nres[1], nres[2], if (all(nres == c(FG_CASE, FG_CTRL))) "same" else "DIFFERENT")
hdr <- readLines(gzfile(fgp), n = 1)
fgcols <- strsplit(sub("^#", "", hdr), "\t")[[1]]
logf("FinnGen header: %s", hdr)
WIN38 <- c(130970832, 132084332)
ps <- jsonlite::fromJSON(inp("LOC_ARG1_6q23_PREP_STAMP"))
if (!all(ps$window_grch38 == WIN38)) stop("STOP: window differs from PREP_STAMP")
t0 <- Sys.time()
fgl <- scanTabix(TabixFile(fgp, index = paste0(fgp, ".tbi")), param = GRanges("6", IRanges(WIN38[1], WIN38[2])))[[1]]
fg <- fread(text = paste(fgl, collapse = "\n"), header = FALSE, sep = "\t", col.names = fgcols,
            colClasses = list(character = c(1, 3, 4, 5, 6)))
logf("FinnGen window 6:%d-%d: %d rows in %.1f s; null rates: %s", WIN38[1], WIN38[2], nrow(fg),
     as.numeric(difftime(Sys.time(), t0, units = "secs")),
     paste(sprintf("%s=%.4f", names(fg), colMeans(is.na(fg) | fg == "")), collapse = " "))
if (nrow(fg) < 1000) stop("STOP: FinnGen window implausibly small")
fwrite(fg, file.path(WORK, "finngen_T2D_ARG1win.tsv"), sep = "\t")
k <- setkey_al(fg$pos, fg$ref, fg$alt); dk <- k[duplicated(k)]
H$key_use <- H$key38
m <- join_report("FinnGen T2D (pos38+allele set)", ifelse(k %in% dk, NA, k), nrow(fg))
fl <- ifelse(fg$alt[m] == H$a1, 1, ifelse(fg$ref[m] == H$a1, -1, NA))
TR$b_fg <- fl * fg$beta[m]; TR$s_fg <- fg$sebeta[m]; TR$f_fg <- ifelse(fl == 1, fg$af_alt[m], 1 - fg$af_alt[m])
afd <- abs(TR$f_fg - H$af_ref)
logf("FinnGen af(a1) vs 503 af_ref: median |diff| %.4f; > 0.2 on %d SNPs (%d non-palindromic); at A %.4f vs %.4f; at B %.4f vs %.4f",
     median(afd, na.rm = TRUE), sum(afd > 0.2, na.rm = TRUE), sum(afd > 0.2 & H$palindromic == 0, na.rm = TRUE),
     TR$f_fg[iA], H$af_ref[iA], TR$f_fg[iB], H$af_ref[iB])

# --- DIAMANTE 2018 noUKB EUR (GRCh37 chr:pos; uncompressed 1.1 GB -> stream-cut, checkpointed)
dp <- inp("OUT_DIAMANTE2018_noUKB_EUR_chrpos")
chain <- inp("REF_LIFTOVER_hg38ToHg19")
blk <- read_chain_blocks(chain, "chr6", "chr6")
lt <- lift_pos(H$pos38, blk)
logf("liftover hg38->hg19 (own chain reader) validated on the harmonised pos38->pos37 pairs: %d/%d identical, %d unmapped",
     sum(lt == H$pos37, na.rm = TRUE), nrow(H), sum(is.na(lt)))
if (sum(lt == H$pos37, na.rm = TRUE) < 0.99 * nrow(H)) stop("STOP: liftover validation failed")
lb <- lift_pos(WIN38, blk)
WIN37 <- c(lb[1], lb[2])
if (any(is.na(WIN37)) || WIN37[1] > min(H$pos37) || WIN37[2] < max(H$pos37)) stop("STOP: lifted window does not cover harmonised pos37")
logf("DIAMANTE window GRCh37 (lifted from GRCh38 %d-%d): %d-%d; covers harmonised pos37 %d-%d", WIN38[1], WIN38[2], WIN37[1], WIN37[2],
     min(H$pos37), max(H$pos37))
dcache <- file.path(WORK, "diamante_chr6_ARG1win.tsv"); dstamp <- file.path(WORK, "diamante_chr6_ARG1win.stamp")
want_stamp <- sprintf("src_md5=%s lo=%d hi=%d", RES[["OUT_DIAMANTE2018_noUKB_EUR_chrpos"]]$md5, WIN37[1], WIN37[2])
if (file.exists(dcache) && file.exists(dstamp) && identical(readLines(dstamp)[1], want_stamp)) {
  logf("DIAMANTE: checkpoint reused (%s; stamp '%s'; line 2 of stamp: %s)", dcache, want_stamp, readLines(dstamp)[2])
} else {
  t0 <- Sys.time()
  con <- file(dp, "r"); dh <- readLines(con, n = 1); ntot <- 0; keep <- character(0)
  repeat {
    x <- readLines(con, n = 1000000L)
    if (!length(x)) break
    ntot <- ntot + length(x)
    y <- x[startsWith(x, "6:")]
    if (length(y)) {
      pos <- as.numeric(sub("^6:([0-9]+)\t.*$", "\\1", y))
      keep <- c(keep, y[pos >= WIN37[1] & pos <= WIN37[2]])
    }
    if (ntot %% 5000000 == 0) logf("DIAMANTE stream: %d rows read, %d kept, %.0f s", ntot, length(keep), as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  close(con)
  writeLines(c(dh, keep), dcache)
  writeLines(c(want_stamp, sprintf("rows_read=%d kept=%d header=%s seconds=%.0f", ntot, length(keep), dh,
                                   as.numeric(difftime(Sys.time(), t0, units = "secs")))), dstamp)
  logf("DIAMANTE stream done: %d data rows read (RESOURCES notes: 21,508,698 data rows by wc -l), %d kept, %.0f s", ntot, length(keep),
       as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
dia <- read.delim(dcache, stringsAsFactors = FALSE, colClasses = c(SNP = "character", EA = "character", NEA = "character"))
logf("DIAMANTE window: %d rows; null rates: %s; md5 of extract %s", nrow(dia),
     paste(sprintf("%s=%.4f", names(dia), colMeans(is.na(dia) | dia == "")), collapse = " "), md5f(dcache))
k <- setkey_al(dia$Pos, dia$EA, dia$NEA); dk <- k[duplicated(k)]
H$key_use <- H$key37
m <- join_report("DIAMANTE (pos37+allele set)", ifelse(k %in% dk, NA, k), nrow(dia))
fl <- ifelse(dia$EA[m] == H$a1, 1, ifelse(dia$NEA[m] == H$a1, -1, NA))
TR$b_dia <- fl * dia$Beta[m]; TR$s_dia <- dia$SE[m]; TR$f_dia <- ifelse(fl == 1, dia$EAF[m], 1 - dia$EAF[m]); TR$p_dia <- dia$Pvalue[m]
DIA_N <- nums_from(RES[["OUT_DIAMANTE2018_noUKB_EUR_chrpos"]]$n_info, "cases ([0-9]+) / controls ([0-9]+)")
DIA_NEFF <- 4 * DIA_N[1] * DIA_N[2] / sum(DIA_N)
logf("DIAMANTE: cases %d controls %d (RESOURCES n_info) N_eff %.1f (computed); at A: b %s se %s p %s EAF(a1) %s; at B: b %s se %s p %s EAF(a1) %s",
     DIA_N[1], DIA_N[2], DIA_NEFF, g7(TR$b_dia[iA]), g7(TR$s_dia[iA]), g7(TR$p_dia[iA]), g7(TR$f_dia[iA]),
     g7(TR$b_dia[iB]), g7(TR$s_dia[iB]), g7(TR$p_dia[iB]), g7(TR$f_dia[iB]))
# placeholder-token check on the parsed numeric columns at A and B
for (cc in c("b_P", "s_P", "b_arg", "s_arg", "b_orn", "s_orn", "b_fg", "s_fg", "b_dia", "s_dia", "p_dia", "f_fg", "f_dia", "f_P", "f_arg", "f_orn"))
  if (any(!is.finite(TR[[cc]][c(iA, iB)]))) stop(sprintf("STOP: non-finite %s at A or B", cc))
fwrite(cbind(H[, c("snp", "rsid_ref", "pos37", "pos38", "a1", "a2", "af_ref", "palindromic")], TR[, -1]),
       file.path(WORK, "own_harmonised_ARG1win.tsv"), sep = "\t")

# ===================== S3. SNP sets and LD (PLINK) =====================
logf("-- S3 SNP sets and PLINK LD")
PLINK <- inp("PLINK19_paper6")
pv <- suppressWarnings(system2(PLINK, "--version", stdout = TRUE, stderr = TRUE)); logf("PLINK --version: %s", paste(pv, collapse = " "))
inp("LD503_bed"); inp("LD503_bim"); inp("LD503_fam")
PREF503 <- sub("\\.bed$", "", RES[["LD503_bed"]]$path)
inp("LD489_chr6_bed"); bim489p <- inp("LD489_chr6_bim")
inp_extra("LD489_chr6_fam", sub("\\.bim$", ".fam", bim489p), RES[["LD489_fam"]]$md5,
          "chr6 .fam of the 489 panel; RESOURCES LD489_fam (chr1) notes all 22 .fam share one md5")
PREF489 <- sub("\\.bim$", "", bim489p)
p1 <- which(H$palindromic == 0)
pst2 <- jsonlite::fromJSON(inp("LOC_ARG1_6q23_PANEL2_STAMP"))
logf("p1np = harmonised minus palindromic: %d SNPs (PANEL2_STAMP definition n_p1np %d)", length(p1), pst2$n_p1np)
if (!(iA %in% p1 && iB %in% p1)) stop("STOP: A/B not in p1np")
q <- function(x) shQuote(x, type = "cmd")
f_p1 <- file.path(WORK, "p1np_ids_503.txt"); writeLines(H$rsid_ref[p1], f_p1)
f_ab <- file.path(WORK, "AB_ids.txt"); writeLines(c(RS_A, RS_B), f_ab)
LDW <- c("--ld-window-r2", "0", "--ld-window", "99999", "--ld-window-kb", "99999")
# 503 freq (A1 record)
run_plink(PLINK, c("--bfile", q(PREF503), "--keep-allele-order", "--extract", q(f_p1), "--freq", "--out", q(file.path(WORK, "ld503_p1np"))), "503 freq")
fr <- read.table(file.path(WORK, "ld503_p1np.frq"), header = TRUE, stringsAsFactors = FALSE, colClasses = c(A1 = "character", A2 = "character"))
if (nrow(fr) != length(p1) || anyDuplicated(fr$SNP)) stop(sprintf("STOP: 503 freq rows %d != p1np %d (duplicate or missing ids)", nrow(fr), length(p1)))
mm <- match(H$rsid_ref[p1], fr$SNP)
A1_503 <- fr$A1[mm]
logf("503 panel: A1 (keep-allele-order) == harmonised a1 on %d/%d p1np SNPs; max |MAF - min(af_ref,1-af_ref)| %.4g", sum(A1_503 == H$a1[p1]),
     length(p1), max(abs(fr$MAF[mm] - pmin(H$af_ref[p1], 1 - H$af_ref[p1]))))
if (any(A1_503 != H$a1[p1])) stop("STOP: 503 A1 != harmonised a1")
# 503 r table for A and B against every p1np SNP
run_plink(PLINK, c("--bfile", q(PREF503), "--keep-allele-order", "--extract", q(f_p1), "--r", "--ld-snp-list", q(f_ab), LDW,
                   "--out", q(file.path(WORK, "ld503_AB_p1np"))), "503 r A,B x p1np")
L503 <- read_ld_table(file.path(WORK, "ld503_AB_p1np.ld"))
assert_ld_complete(L503, c(RS_A, RS_B), H$rsid_ref[p1], "503")
logf("503 LD table: %d rows, all %d x 2 requested pairs returned and finite", nrow(L503), length(p1))
# mutation: one SNP id mutated must trip the assertion
f_abm <- file.path(WORK, "AB_ids_MUTATED.txt"); writeLines(c(RS_A, paste0(RS_B, "_MUT")), f_abm)
mut_ids <- ifelse(H$rsid_ref[p1] == RS_B, paste0(RS_B, "_MUT"), H$rsid_ref[p1])
f_p1m <- file.path(WORK, "p1np_ids_503_MUTATED.txt"); writeLines(mut_ids, f_p1m)
run_plink(PLINK, c("--bfile", q(PREF503), "--keep-allele-order", "--extract", q(f_p1m), "--r", "--ld-snp-list", q(f_abm), LDW,
                   "--out", q(file.path(WORK, "ld503_AB_MUTATED"))), "503 r MUTATION")
mt_msg <- tryCatch({ assert_ld_complete(read_ld_table(file.path(WORK, "ld503_AB_MUTATED.ld")), c(RS_A, paste0(RS_B, "_MUT")),
                                        mut_ids, "503 MUTATION"); "NOT TRIPPED" },
                   error = function(e) conditionMessage(e))
logf("LD assertion mutation test (B id -> %s_MUT): %s", RS_B, if (mt_msg == "NOT TRIPPED") "FAILED TO TRIP" else paste("TRIPPED as expected:", mt_msg))
if (mt_msg == "NOT TRIPPED") stop("STOP: LD assertion mutation did not trip")
r503_col <- function(k) { s <- L503[L503$SNP_A == k, ]; s$R[match(H$rsid_ref[p1], s$SNP_B)] }
R503_A <- r503_col(RS_A); R503_B <- r503_col(RS_B)
r503 <- R503_A[match(iB, p1)]
logf("PLINK r(A,B) 503 = %s (A row) / %s (B row)", g7(r503), g7(R503_B[match(iA, p1)]))
# 489: p2np = p1np present in panel 2 (pos37 + allele set) with panel-2 MAF >= 0.01
b489 <- read.table(bim489p, header = FALSE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "integer", "character", "character"))
logf("489 chr6 bim: %d rows", nrow(b489))
k489 <- setkey_al(b489$V4, b489$V5, b489$V6); dk <- k489[duplicated(k489)]
m489 <- match(H$key37[p1], ifelse(k489 %in% dk, NA, k489))
logf("489 bim match of p1np (pos37+allele set): %d/%d (absent %d; ambiguous duplicated keys in window %d); rsid concordant %d/%d",
     sum(!is.na(m489)), length(p1), sum(is.na(m489)), sum(H$key37[p1] %in% dk), sum(b489$V2[m489] == H$rsid_ref[p1], na.rm = TRUE), sum(!is.na(m489)))
cand <- p1[!is.na(m489)]; cid <- b489$V2[m489[!is.na(m489)]]
if (anyDuplicated(cid)) stop("STOP: duplicated 489 ids among candidates")
f_c489 <- file.path(WORK, "p1np_in489_ids.txt"); writeLines(cid, f_c489)
run_plink(PLINK, c("--bfile", q(PREF489), "--keep-allele-order", "--extract", q(f_c489), "--freq", "--out", q(file.path(WORK, "ld489_cand"))), "489 freq")
fr2 <- read.table(file.path(WORK, "ld489_cand.frq"), header = TRUE, stringsAsFactors = FALSE, colClasses = c(A1 = "character", A2 = "character"))
if (nrow(fr2) != length(cid) || anyDuplicated(fr2$SNP)) stop("STOP: 489 freq rows")
mf <- match(cid, fr2$SNP)
keep2 <- fr2$MAF[mf] >= 0.01
p2 <- cand[keep2]; id489 <- cid[keep2]; A1_489 <- fr2$A1[mf][keep2]
flip489 <- ifelse(A1_489 == H$a1[p2], 1, ifelse(A1_489 == H$a2[p2], -1, NA))
if (any(is.na(flip489))) stop("STOP: 489 A1 not in allele set")
logf("p2np: %d SNPs (dropped MAF<0.01 in 489: %d; PANEL2_STAMP definition n_p2np %d); 489 A1 == harmonised a1 on %d, flipped on %d",
     length(p2), sum(!keep2), pst2$n_p2np, sum(flip489 == 1), sum(flip489 == -1))
if (!(iA %in% p2 && iB %in% p2)) stop("STOP: A/B not in p2np")
idA489 <- id489[p2 == iA]; idB489 <- id489[p2 == iB]
f_p2 <- file.path(WORK, "p2np_ids_489.txt"); writeLines(id489, f_p2)
f_ab2 <- file.path(WORK, "AB_ids_489.txt"); writeLines(c(idA489, idB489), f_ab2)
run_plink(PLINK, c("--bfile", q(PREF489), "--keep-allele-order", "--extract", q(f_p2), "--r", "--ld-snp-list", q(f_ab2), LDW,
                   "--out", q(file.path(WORK, "ld489_AB_p2np"))), "489 r A,B x p2np")
L489 <- read_ld_table(file.path(WORK, "ld489_AB_p2np.ld"))
assert_ld_complete(L489, c(idA489, idB489), id489, "489")
logf("489 LD table: %d rows, all %d x 2 requested pairs returned and finite", nrow(L489), length(p2))
r489_col <- function(k, fk) { s <- L489[L489$SNP_A == k, ]; s$R[match(id489, s$SNP_B)] * flip489 * fk }
R489_A <- r489_col(idA489, flip489[p2 == iA]); R489_B <- r489_col(idB489, flip489[p2 == iB])
r489 <- R489_A[match(iB, p2)]
logf("PLINK r(A,B) 489 (aligned to a1) = %s", g7(r489))
addr("S4.ld.r.503", "4(ii); Rule 3(ii)", "r(A,B), A1 dosage correlation", "1000G EUR g1000_eur (LD503_bed)", "503", value = r503,
     n = 503, unit = "r", notes = sprintf("PLINK 1.9 --r --keep-allele-order; A=%s B=%s; aligned to harmonised a1 (A1 = a1 on both)", RS_A, RS_B))
addr("S4.ld.r.489", "4(ii); Rule 3(ii)", "r(A,B), A1 dosage correlation", "1000G.EUR.QC chr6 (LD489_chr6_bed)", "489", value = r489,
     n = 489, unit = "r", notes = "PLINK 1.9 --r --keep-allele-order; sign aligned to harmonised a1")
assert_target("rule_rAB_503", r503, -0.0445, 0.00005, "PLINK --r on LD503 (A row, B column); target quoted in R1_ARG1_rule (definition; weak evidence)",
              "R1_ARG1_rule")

# ===================== S4. dosage extract (read in place), full-panel check, bootstrap of r =====================
logf("-- S4 dosage extract and bootstrap")
i8p <- inp("R1_ARG1_REF503_I8")
G <- readBin(i8p, what = "integer", size = 1, signed = TRUE, n = 503 * 3067 + 10)
if (length(G) != 503 * 3067) stop("STOP: i8 size")
G <- matrix(G, nrow = 503, byrow = TRUE)
if (!all(G %in% 0:2)) stop("STOP: i8 values outside 0..2")
afd <- max(abs(colMeans(G) / 2 - H$af_ref))
logf("i8: 503 x 3067 read in place; max |mean/2 - af_ref| %.3g", afd)
if (afd > 1e-5) stop("STOP: i8 dosage != af_ref")
gA <- G[, iA]; gB <- G[, iB]
r_dos <- cor(gA, gB)
dRA <- cor(G[, p1], gA)[, 1]; dRB <- cor(G[, p1], gB)[, 1]
logf("dosage r(A,B) %.7f vs PLINK %.7f (|diff| %.2g); window columns: max |dosage r - PLINK r| A %.2g B %.2g", r_dos, r503, abs(r_dos - r503),
     max(abs(dRA - R503_A)), max(abs(dRB - R503_B)))
if (abs(r_dos - r503) > 1e-5 || max(abs(dRA - R503_A)) > 1e-5 || max(abs(dRB - R503_B)) > 1e-5) stop("STOP: dosage r != PLINK r")
rm(G); invisible(gc())
bfile <- file.path(WORK, "boot_r_503.tsv"); bstate <- file.path(WORK, "boot_rng_state.rds"); bstamp <- file.path(WORK, "boot_stamp.txt")
want <- sprintf("seed=%d nboot=%d i8_md5=%s A=%s B=%s RNG=%s", SEED, NBOOT, RES[["R1_ARG1_REF503_I8"]]$md5, RS_A, RS_B, paste(RNGkind(), collapse = "/"))
done <- 0L
if (file.exists(bfile) && file.exists(bstamp) && identical(readLines(bstamp), want) && file.exists(bstate)) {
  done <- nrow(fread(file = bfile))
  if (done %% 500L != 0L) { logf("bootstrap checkpoint has %d rows (not a chunk multiple): restarting", done); done <- 0L }
}
if (done == 0L) { writeLines(want, bstamp); writeLines("draw\tr\tcarriersB", bfile); set.seed(SEED) } else {
  assign(".Random.seed", readRDS(bstate), envir = .GlobalEnv); logf("bootstrap: resuming after %d draws", done) }
t0 <- Sys.time()
while (done < NBOOT) {
  ch <- t(vapply(seq_len(500L), function(i) { idx <- sample.int(503L, 503L, replace = TRUE)
    c(suppressWarnings(cor(gA[idx], gB[idx])), sum(gB[idx] > 0)) }, numeric(2)))
  fwrite(data.table(draw = done + seq_len(500L), r = ch[, 1], carriersB = ch[, 2]), bfile, sep = "\t", append = TRUE)
  done <- done + 500L
  saveRDS(.Random.seed, bstate)
  logf("bootstrap progress: %d/%d draws (%.1f s)", done, NBOOT, as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
BT <- fread(file = bfile)
if (nrow(BT) != NBOOT || anyDuplicated(BT$draw)) stop("STOP: bootstrap rows")
rstar <- BT$r
logf("bootstrap r: %d draws, %d non-finite (zero-variance resample), min carriers of B per draw %d", NBOOT, sum(!is.finite(rstar)), min(BT$carriersB))
rq <- quantile(rstar, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
addr("S4.ld.boot.r_lo", "4(ii); Rule 3(ii)", "bootstrap 2.5th percentile of r(A,B)", "R1_ARG1_REF503_I8 (read in place)", "503", value = rq[1],
     n = sum(is.finite(rstar)), unit = "r", notes = sprintf("DEVIATION: bootstrap outside PLINK (section 10); %d resamples of 503 people, seed %d, R %s; percentile (type 7)", NBOOT, SEED, paste(RNGkind(), collapse = "/")))
addr("S4.ld.boot.r_hi", "4(ii); Rule 3(ii)", "bootstrap 97.5th percentile of r(A,B)", "R1_ARG1_REF503_I8 (read in place)", "503", value = rq[2],
     n = sum(is.finite(rstar)), unit = "r", notes = "DEVIATION: bootstrap outside PLINK (section 10); as S4.ld.boot.r_lo")
R_BOOT_INCL0 <- rq[1] <= 0 && rq[2] >= 0

# ===================== S5. the grid and the A-vs-B tests =====================
logf("-- S5 grid and A-vs-B tests")
tr_ab <- function(col, sefromp = FALSE) {
  b <- TR[[paste0("b_", col)]]; s <- TR[[paste0("s_", col)]]
  if (sefromp) s <- abs(b) / qnorm(TR$p_dia / 2, lower.tail = FALSE)
  list(bA = b[iA], sA = s[iA], bB = b[iB], sB = s[iB], fA = TR[[paste0("f_", col)]][iA], fB = TR[[paste0("f_", col)]][iB])
}
TP <- tr_ab("P"); TARG <- tr_ab("arg"); TORN <- tr_ab("orn"); TFG <- tr_ab("fg"); TDIA <- tr_ab("dia"); TDIAp <- tr_ab("dia", TRUE)
for (nm in c("TP", "TARG", "TORN", "TFG", "TDIA", "TDIAp")) { t <- get(nm)
  logf("inputs %-6s A: b %s se %s f %s | B: b %s se %s f %s", nm, g7(t$bA), g7(t$sA), g7(t$fA), g7(t$bB), g7(t$sB), g7(t$fB)) }
cell <- function(Yt, Xt, r, estimand) {
  Y <- est_trait(Yt, r, estimand, "z"); X <- est_trait(Xt, r, estimand, "z")
  aw2 <- ab_test(Y, X, "2nd", TRUE); ign <- ab_test(Y, X, "2nd", FALSE); aw1 <- ab_test(Y, X, "1st", TRUE)
  cj <- ab_test(est_trait(Yt, r, estimand, "beta"), est_trait(Xt, r, estimand, "beta"), "2nd", TRUE)
  pw <- power_calc(aw2$RA, aw2$vA, X$bB, X$sB, Y$sB)
  list(Y = Y, X = X, aw2 = aw2, ign = ign, aw1 = aw1, cj = cj, pw = pw, obs_zB = Y$bB / Y$sB)
}
RESN <- list(FinnGen = list(t = TFG, lab = "FinnGen R10 T2D (OUT_FG_R10_T2D)", n = FG_N),
             DIAMANTE = list(t = TDIA, lab = "DIAMANTE 2018 noUKB EUR (OUT_DIAMANTE2018_noUKB_EUR_chrpos)", n = sum(DIA_N)),
             DIAMANTE_SEfromP = list(t = TDIAp, lab = "DIAMANTE 2018 noUKB EUR, SE = |Beta|/z(Pvalue) (SE-from-p variant)", n = sum(DIA_N)))
PAN <- list(`503` = r503, `489` = r489)
CELLS <- list()
for (rn in names(RESN)) for (es in c("marg", "cond")) for (pn in names(PAN)) {
  CELLS[[paste(rn, es, pn, sep = "|")]] <- cell(RESN[[rn]]$t, TP, PAN[[pn]], es)
}
est_lab <- c(marg = "marginal", cond = "conditional (z-space, A|B and B|A)")
cojo_note <- function(es) if (es == "marg") "p_cojo: marginal estimand has no conditioning step, so beta-space = z-space (value = p_aware_2nd)" else
  "p_cojo: beta-space (COJO) conditioning b_A|B = (b_A - r sqrt(v_B/v_A) b_B)/(1-r^2), v = 2f(1-f) from each source's own allele frequency"
for (rn in names(RESN)) for (es in c("marg", "cond")) for (pn in names(PAN)) {
  C <- CELLS[[paste(rn, es, pn, sep = "|")]]
  pre <- if (pn == "503") sprintf("S4.grid.%s.%s.", rn, es) else sprintf("S4.grid.%s.%s.489.", rn, es)
  ref <- "4(i) stability grid; Rule 5"
  est <- sprintf("T2D log-odds per SD ARG1, %s; A=%s B=%s", est_lab[es], RS_A, RS_B)
  res <- RESN[[rn]]$lab; nn <- RESN[[rn]]$n
  addr(paste0(pre, "ratioA"), ref, est, res, pn, "aware", "2nd", C$aw2$RA, sqrt(C$aw2$vA), n = nn, unit = "log-odds per SD ARG1")
  addr(paste0(pre, "ratioB"), ref, est, res, pn, "aware", "2nd", C$aw2$RB, sqrt(C$aw2$vB), n = nn, unit = "log-odds per SD ARG1")
  addr(paste0(pre, "diff"), ref, paste(est, "; A-B"), res, pn, "aware", "2nd", C$aw2$d, C$aw2$se, C$aw2$lo, C$aw2$hi, C$aw2$p, nn, "log-odds per SD ARG1")
  addr(paste0(pre, "diff_se"), ref, paste(est, "; A-B"), res, pn, "aware", "2nd", C$aw2$se, n = nn, unit = "log-odds per SD ARG1")
  addr(paste0(pre, "diff_lo"), ref, paste(est, "; A-B"), res, pn, "aware", "2nd", C$aw2$lo, n = nn, unit = "log-odds per SD ARG1", notes = "95% CI, normal")
  addr(paste0(pre, "diff_hi"), ref, paste(est, "; A-B"), res, pn, "aware", "2nd", C$aw2$hi, n = nn, unit = "log-odds per SD ARG1", notes = "95% CI, normal")
  addr(paste0(pre, "p_aware_2nd"), ref, est, res, pn, "aware", "2nd", C$aw2$p, p = C$aw2$p, n = nn, unit = "p", notes = "Rule 5 primary")
  addr(paste0(pre, "p_ign"), ref, est, res, pn, "ign", "2nd", C$ign$p, p = C$ign$p, n = nn, unit = "p")
  addr(paste0(pre, "p_aware_1st"), ref, est, res, pn, "aware", "1st", C$aw1$p, p = C$aw1$p, n = nn, unit = "p")
  addr(paste0(pre, "p_cojo"), ref, est, res, pn, "aware", "2nd", C$cj$p, p = C$cj$p, n = nn, unit = "p", notes = cojo_note(es))
  addr(paste0(pre, "exp_zB"), ref, est, res, pn, "aware", "2nd", C$pw$exp_z, n = nn, unit = "z",
       notes = "expected z at B if B followed A's ratio: ratioA x b_P,B / se_T2D,B (signed; the rule quotes |z|)")
  addr(paste0(pre, "obs_zB"), ref, est, res, pn, "", "", C$obs_zB, n = nn, unit = "z", notes = "observed T2D z at B in this cell")
  addr(paste0(pre, "power_point"), ref, est, res, pn, "aware", "2nd", C$pw$power_point, n = nn, unit = "probability",
       notes = "P(|Z|>1.96) at Z ~ N(exp_zB, 1)")
  addr(paste0(pre, "power_pred"), ref, est, res, pn, "aware", "2nd", C$pw$power_pred, n = nn, unit = "probability",
       notes = sprintf("predictive: Z ~ N(exp_zB, 1 + (b_PB^2 Var(ratioA) + ratioA^2 se_PB^2)/se_YB^2), sd %.3f", C$pw$pred_sd))
}

# arginine (primary: marginal, aware, 503), ornithine (primary: marginal, aware, 503), test (iii) (FinnGen, both estimands)
ror_rows <- function(C, pre, ref, est, res, pn) {
  a <- C$aw2
  fi <- fieller(a$RA, a$RB, a$vA, a$vB, a$cov)
  ld <- logdelta_ci(a$RA, a$RB, a$vA, a$vB, a$cov)
  addr(paste0(pre, "ror"), ref, paste(est, "; ratio of ratios A/B"), res, pn, "aware", "2nd", a$RA / a$RB, lo = fi$lo, hi = fi$hi, p = a$p,
       unit = "ratio", notes = sprintf("CI: Fieller 95%% (%s; coherent with the A-B test); log-scale delta-method CI beside: %s-%s", fi$status, g7(ld[1]), g7(ld[2])))
  addr(paste0(pre, "ror_lo"), ref, paste(est, "; ratio of ratios A/B"), res, pn, "aware", "2nd", fi$lo, unit = "ratio", notes = "Fieller 95%")
  addr(paste0(pre, "ror_hi"), ref, paste(est, "; ratio of ratios A/B"), res, pn, "aware", "2nd", fi$hi, unit = "ratio", notes = "Fieller 95%")
  fi
}
ab_rows <- function(C, pre, ref, est, res, pn, es, unit, kA = "ratioA", kB = "ratioB") {
  addr(paste0(pre, kA), ref, est, res, pn, "aware", "2nd", C$aw2$RA, sqrt(C$aw2$vA), unit = unit)
  addr(paste0(pre, kB), ref, est, res, pn, "aware", "2nd", C$aw2$RB, sqrt(C$aw2$vB), unit = unit)
  addr(paste0(pre, "diff"), ref, paste(est, "; A-B"), res, pn, "aware", "2nd", C$aw2$d, C$aw2$se, C$aw2$lo, C$aw2$hi, C$aw2$p, unit = unit)
  addr(paste0(pre, "p_aware_2nd"), ref, est, res, pn, "aware", "2nd", C$aw2$p, p = C$aw2$p, unit = "p", notes = "Rule 5 primary")
  addr(paste0(pre, "p_ign"), ref, est, res, pn, "ign", "2nd", C$ign$p, p = C$ign$p, unit = "p")
  addr(paste0(pre, "p_aware_1st"), ref, est, res, pn, "aware", "1st", C$aw1$p, p = C$aw1$p, unit = "p")
  addr(paste0(pre, "p_cojo"), ref, est, res, pn, "aware", "2nd", C$cj$p, p = C$cj$p, unit = "p", notes = cojo_note(es))
}
ARGC <- list(); ORNC <- list(); T3C <- list()
for (es in c("marg", "cond")) for (pn in names(PAN)) {
  ARGC[[paste(es, pn)]] <- cell(TARG, TP, PAN[[pn]], es)
  ORNC[[paste(es, pn)]] <- cell(TORN, TP, PAN[[pn]], es)
  for (rn in names(RESN)) T3C[[paste(rn, es, pn)]] <- cell(RESN[[rn]]$t, TARG, PAN[[pn]], es)
}
U_M <- "arginine units per SD ARG1 (CLSA scale; not confirmed locally)"
for (es in c("marg", "cond")) for (pn in names(PAN)) {
  C <- ARGC[[paste(es, pn)]]
  pre <- if (pn == "503") sprintf("S4.arg.%s.", es) else sprintf("S4.arg.%s.489.", es)
  est <- sprintf("arginine per SD ARG1, %s; A=%s B=%s", est_lab[es], RS_A, RS_B)
  ab_rows(C, pre, "4 Arginine; Rule 5", est, "CLSA arginine GCST90200372 / UKB-PPP ARG1 OID21445", pn, es, U_M)
  ror_rows(C, pre, "4 Arginine; Rule 5", est, "CLSA arginine GCST90200372 / UKB-PPP ARG1 OID21445", pn)
}
pA_arg <- 2 * pnorm(-abs(TARG$bA / TARG$sA)); pB_arg <- 2 * pnorm(-abs(TARG$bB / TARG$sB))
addr("S4.arg.pA", "4 Arginine (opening clause)", "arginine marginal association at A", "CLSA arginine GCST90200372", "", value = pA_arg, p = pA_arg,
     n = N_ARG, unit = "p", notes = sprintf("from b/se (z %.3f); file p_value %s", TARG$bA / TARG$sA, g7(TR$pfile_arg[iA])))
addr("S4.arg.pB", "4 Arginine (opening clause)", "arginine marginal association at B", "CLSA arginine GCST90200372", "", value = pB_arg, p = pB_arg,
     n = N_ARG, unit = "p", notes = sprintf("from b/se (z %.3f); file p_value %s", TARG$bB / TARG$sB, g7(TR$pfile_arg[iB])))
U_O <- "ornithine units per SD ARG1 (CLSA scale; not confirmed locally)"
for (es in c("marg", "cond")) for (pn in names(PAN)) {
  C <- ORNC[[paste(es, pn)]]
  pre <- if (es == "marg") sprintf("S4.orn.%s.", pn) else sprintf("S4.orn.%s.cond.", pn)
  est <- sprintf("ornithine per SD ARG1, %s; A=%s B=%s", est_lab[es], RS_A, RS_B)
  ab_rows(C, pre, "4(vi) Ornithine; Rule 5", est, "CLSA ornithine GCST90200404 / UKB-PPP ARG1 OID21445", pn, es, U_O, "kA", "kB")
}
pA_orn <- 2 * pnorm(-abs(TORN$bA / TORN$sA)); pB_orn <- 2 * pnorm(-abs(TORN$bB / TORN$sB))
addr("S4.orn.pA", "4(vi) Ornithine", "ornithine marginal association at A", "CLSA ornithine GCST90200404", "", value = pA_orn, p = pA_orn, n = N_ORN, unit = "p")
addr("S4.orn.pB", "4(vi) Ornithine", "ornithine marginal association at B", "CLSA ornithine GCST90200404", "", value = pB_orn, p = pB_orn, n = N_ORN, unit = "p")
U_T <- "T2D log-odds per arginine unit"
for (rn in names(RESN)) for (es in c("marg", "cond")) for (pn in names(PAN)) {
  C <- T3C[[paste(rn, es, pn)]]
  pre <- sprintf("S4.t3.%s.%s.%s.", rn, es, pn)
  est <- sprintf("test (iii): T2D per unit arginine, %s; A=%s B=%s", est_lab[es], RS_A, RS_B)
  res <- paste(RESN[[rn]]$lab, "/ CLSA arginine")
  addr(paste0(pre, "ratioA"), "4(iii)", est, res, pn, "aware", "2nd", C$aw2$RA, sqrt(C$aw2$vA), unit = U_T)
  addr(paste0(pre, "ratioB"), "4(iii)", est, res, pn, "aware", "2nd", C$aw2$RB, sqrt(C$aw2$vB), unit = U_T)
  addr(paste0(pre, "diff"), "4(iii)", paste(est, "; A-B"), res, pn, "aware", "2nd", C$aw2$d, C$aw2$se, C$aw2$lo, C$aw2$hi, C$aw2$p, unit = U_T)
  addr(paste0(pre, "2nd.p_aware"), "4(iii); Rule 5", est, res, pn, "aware", "2nd", C$aw2$p, p = C$aw2$p, unit = "p",
       notes = if (rn == "FinnGen" && pn == "503") "primary (family member; both estimands required)" else "printed beside")
  addr(paste0(pre, "1st.p_aware"), "4(iii); Rule 5", est, res, pn, "aware", "1st", C$aw1$p, p = C$aw1$p, unit = "p")
  addr(paste0(pre, "2nd.p_ign"), "4(iii); Rule 5", est, res, pn, "ign", "2nd", C$ign$p, p = C$ign$p, unit = "p")
  addr(paste0(pre, "2nd.p_cojo"), "4(iii); Rule 5", est, res, pn, "aware", "2nd", C$cj$p, p = C$cj$p, unit = "p", notes = cojo_note(es))
}

# ---- bootstrap shares (Rule 3(ii)) and the r grid
logf("-- bootstrap shares")
concl <- function(p, d) paste(p < 0.05, sign(d))
boot_eval <- function(tag, Yt, Xt, es, ref) {
  full <- cell(Yt, Xt, r503, es)$aw2; c489 <- cell(Yt, Xt, r489, es)$aw2
  bs <- cell(Yt, Xt, rstar, es)$aw2
  bs$d <- rep_len(bs$d, length(bs$p))
  ok <- is.finite(bs$p)
  sh_lt05 <- mean(bs$p[ok] < 0.05)
  sh_same <- mean(concl(bs$p[ok], bs$d[ok]) == concl(full$p, full$d))
  sh_side <- mean((bs$p[ok] < 0.05) == (full$p < 0.05))
  pq <- quantile(bs$p[ok], c(0.025, 0.975), names = FALSE)
  same489 <- concl(c489$p, c489$d) == concl(full$p, full$d)
  robust <- same489 && sh_same >= 0.90
  nt <- sprintf("cell %s; full 503 p %s (A-B %s), 489 p %s (A-B %s); 489 same conclusion %s; bootstrap share same conclusion (side of 0.05 AND sign of A-B) %.4f, side only %.4f; %d finite draws",
                tag, fmtp(full$p), g7(full$d), fmtp(c489$p), g7(c489$d), same489, sh_same, sh_side, sum(ok))
  addr(sprintf("S4.ld.boot.%s.share_lt05", tag), ref, "share of 5,000 bootstrap draws with p < 0.05", "R1_ARG1_REF503_I8 r(A,B) resamples", "503",
       "aware", "2nd", sh_lt05, n = sum(ok), unit = "share", notes = paste("DEVIATION: bootstrap outside PLINK (section 10);", nt))
  addr(sprintf("S4.ld.boot.%s.p_lo", tag), ref, "2.5th percentile of bootstrap p", "R1_ARG1_REF503_I8 r(A,B) resamples", "503", "aware", "2nd", pq[1],
       n = sum(ok), unit = "p", notes = "DEVIATION: bootstrap outside PLINK (section 10)")
  addr(sprintf("S4.ld.boot.%s.p_hi", tag), ref, "97.5th percentile of bootstrap p", "R1_ARG1_REF503_I8 r(A,B) resamples", "503", "aware", "2nd", pq[2],
       n = sum(ok), unit = "p", notes = "DEVIATION: bootstrap outside PLINK (section 10)")
  addr(sprintf("S4.ld.boot.%s.share_same", tag), ref, "share of draws with the full-panel conclusion (side of p 0.05 and sign of A-B)",
       "R1_ARG1_REF503_I8 r(A,B) resamples", "503", "aware", "2nd", sh_same, n = sum(ok), unit = "share",
       notes = sprintf("Rule 3(ii) robust = 489 same AND share >= 0.90 -> %s; side-only share %.4f", robust, sh_side))
  logf("BOOT %s", nt)
  list(full = full, c489 = c489, share_same = sh_same, share_side = sh_side, share_lt05 = sh_lt05, same489 = same489, robust = robust, plo = pq[1], phi = pq[2])
}
BOOT <- list()
for (rn in names(RESN)) for (es in c("marg", "cond"))
  BOOT[[sprintf("%s.%s.aware2nd", rn, es)]] <- boot_eval(sprintf("%s.%s.aware2nd", rn, es), RESN[[rn]]$t, TP, es,
     if (rn == "FinnGen" && es == "cond") "Rule 3(ii) named primary cell (T2D)" else if (rn == "DIAMANTE" && es == "marg") "Amendment 1 A2: Rule 3(ii) on DIAMANTE marginal aware" else "Rule 3(ii) printed beside")
for (es in c("marg", "cond")) {
  BOOT[[sprintf("arg.%s.aware2nd", es)]] <- boot_eval(sprintf("arg.%s.aware2nd", es), TARG, TP, es,
     if (es == "marg") "Rule 3(ii) named primary cell (arginine)" else "Rule 3(ii) printed beside")
  BOOT[[sprintf("orn.%s.aware2nd", es)]] <- boot_eval(sprintf("orn.%s.aware2nd", es), TORN, TP, es,
     if (es == "marg") "Rule 3(ii) named primary cell (ornithine)" else "Rule 3(ii) printed beside")
  BOOT[[sprintf("t3.FinnGen.%s.aware2nd", es)]] <- boot_eval(sprintf("t3.FinnGen.%s.aware2nd", es), TFG, TARG, es, "Rule 3(ii) named primary cell (test (iii), both estimands)")
  BOOT[[sprintf("t3.DIAMANTE.%s.aware2nd", es)]] <- boot_eval(sprintf("t3.DIAMANTE.%s.aware2nd", es), TDIA, TARG, es, "printed beside")
}
for (es in c("cond", "marg")) {
  b <- BOOT[[sprintf("t3.FinnGen.%s.aware2nd", es)]]
  addr(sprintf("S4.t3.boot.%s.share", es), "4(iii); Rule 3(ii)", sprintf("test (iii) FinnGen %s aware 2nd: share of draws with the full-panel conclusion", est_lab[es]),
       "R1_ARG1_REF503_I8 r(A,B) resamples", "503", "aware", "2nd", b$share_same, n = sum(is.finite(rstar)), unit = "share",
       notes = sprintf("DEVIATION: bootstrap outside PLINK (section 10); side-only share %.4f; share p<0.05 %.4f; 489 same %s; robust %s", b$share_side, b$share_lt05, b$same489, b$robust))
}
logf("-- r grid")
RG <- seq(-0.95, 0.95, by = 0.0005)
crossings <- function(pv) { s <- sign(pv - 0.05); i <- which(s[-1] != s[-length(s)] & is.finite(s[-1]) & is.finite(s[-length(s)]))
  if (!length(i)) return(numeric(0))
  RG[i] + (0.05 - pv[i]) * (RG[i + 1] - RG[i]) / (pv[i + 1] - pv[i]) }
grid_rows <- list()
rgrid_cell <- function(tag, Yt, Xt, es, primary = FALSE) {
  cc <- cell(Yt, Xt, RG, es)
  ca <- crossings(cc$aw2$p); ci <- crossings(cc$ign$p)
  near <- function(x) if (length(x)) x[which.min(abs(x - r503))] else "none in grid [-0.95, 0.95] (p does not cross 0.05)"
  grid_rows[[tag]] <<- data.frame(cell = tag, r = RG, p_aware_2nd = cc$aw2$p, p_ign = cc$ign$p, d = cc$aw2$d)
  for (v in c("aware", "ign")) {
    x <- if (v == "aware") ca else ci
    id <- if (primary) sprintf("S4.ld.rgrid.crossing.%s", v) else sprintf("S4.ld.rgrid.crossing.%s.%s", tag, v)
    addr(id, "4(ii) r grid", sprintf("r(A,B) at which p crosses 0.05, %s", tag), "r grid -0.95..0.95 step 0.0005 (inputs fixed, r varied)", "grid",
         v, "2nd", near(x), unit = "r",
         notes = sprintf("%s; all crossings in grid: %s; observed r 503 %s (p %s), 489 %s; bootstrap 95%% r %s to %s", if (primary) "primary cell FinnGen conditional aware" else "printed beside",
                         if (length(x)) paste(sprintf("%.4f", x), collapse = ",") else "none", g7(r503), fmtp(if (v == "aware") cell(Yt, Xt, r503, es)$aw2$p else cell(Yt, Xt, r503, es)$ign$p),
                         g7(r489), g7(rq[1]), g7(rq[2])))
  }
}
rgrid_cell("FinnGen.cond", TFG, TP, "cond", primary = TRUE)
rgrid_cell("FinnGen.marg", TFG, TP, "marg")
rgrid_cell("DIAMANTE.marg", TDIA, TP, "marg")
rgrid_cell("DIAMANTE.cond", TDIA, TP, "cond")
rgrid_cell("arg.marg", TARG, TP, "marg")
rgrid_cell("orn.marg", TORN, TP, "marg")
rgrid_cell("t3.FinnGen.cond", TFG, TARG, "cond")
rgrid_cell("t3.FinnGen.marg", TFG, TARG, "marg")
fwrite(rbindlist(grid_rows), file.path(WORK, "rgrid.tsv"), sep = "\t")

# ===================== S6. Rule 3(iv) conditional colocalization (fresh R coloc 5.2.3, p12 = 1e-5) =====================
logf("-- S6 conditional colocalization clause")
FG_S <- FG_CASE / FG_N
afok_fg <- abs(TR$f_fg - H$af_ref) <= 0.2
logf("coloc settings: p1 = p2 = 1e-4, p12 = 1e-5; protein quant sdY 1 (per SD NPX) N %d; metabolites quant sdY 1 N %d/%d; T2D cc s %.5f N %d; SNPs with |af(FinnGen) - af_ref| > 0.2 dropped for T2D (%d in p1np)",
     N_P, N_ARG, N_ORN, FG_S, FG_N, sum(!afok_fg[p1], na.rm = TRUE))
TRAITS <- list(arginine = list(col = "arg", type = "quant", N = N_ARG), ornithine = list(col = "orn", type = "quant", N = N_ORN),
               T2D = list(col = "fg", type = "cc", N = FG_N))
COL <- list()
for (pn in c("503", "489")) {
  set <- if (pn == "503") p1 else p2
  rA <- if (pn == "503") R503_A else R489_A; rB <- if (pn == "503") R503_B else R489_B
  kA <- match(iA, set); kB <- match(iB, set)
  ids <- H$snp[set]
  for (sg in c("A", "B")) {
    # signal A: condition on B (only A left); signal B: condition on A
    rk <- if (sg == "A") rB else rA; kk <- if (sg == "A") kB else kA; at <- if (sg == "A") kA else kB
    Pc <- cond_window(TR$b_P[set], TR$s_P[set], rk, kk)
    for (tn in names(TRAITS)) {
      tt <- TRAITS[[tn]]
      b <- TR[[paste0("b_", tt$col)]][set]; s <- TR[[paste0("s_", tt$col)]][set]
      if (tn == "T2D") { bad <- is.na(afok_fg[set]) | !afok_fg[set]; b[bad] <- NA; s[bad] <- NA }
      Yc <- cond_window(b, s, rk, kk)
      cz <- Yc$z[at]
      res <- coloc_pair(Pc$b, Pc$s, N_P, Yc$b, Yc$s, tt$N, tt$type, FG_S, ids, 1e-5)
      res6 <- coloc_pair(Pc$b, Pc$s, N_P, Yc$b, Yc$s, tt$N, tt$type, FG_S, ids, 1e-6)
      COL[[paste(tn, sg, pn)]] <- list(z = cz, H4 = res$H4, H3 = res$H3, nsnp = res$nsnp, H4_6 = res6$H4)
      nt <- sprintf("%s at signal %s conditioned on %s; coloc.abf p1=p2=1e-4 p12=1e-5 on %d SNPs (PP.H3 %.4f; collinear dropped %d); PP.H4 at p12=1e-6 %.4f",
                    tn, sg, if (sg == "A") RS_B else RS_A, res$nsnp, res$H3, Yc$n_collinear, res6$H4)
      addr(sprintf("S4.coloc_clause.%s.%s.%s.z", tn, sg, pn), "Rule 3(iv)", sprintf("%s conditional z at signal %s", tn, sg),
           if (tn == "T2D") "FinnGen R10 T2D" else sprintf("CLSA %s", tn), pn, value = cz, n = res$nsnp, unit = "z", notes = nt)
      addr(sprintf("S4.coloc_clause.%s.%s.%s.PPH4", tn, sg, pn), "Rule 3(iv)", sprintf("conditional PP.H4, ARG1 protein x %s, signal %s", tn, sg),
           if (tn == "T2D") "FinnGen R10 T2D / UKB-PPP ARG1" else sprintf("CLSA %s / UKB-PPP ARG1", tn), pn, value = res$H4, n = res$nsnp, unit = "PP.H4", notes = nt)
      logf("COLOC %s signal %s panel %s: z %.3f PP.H4 %.4f (H3 %.4f, nsnp %d; p12=1e-6 H4 %.4f)", tn, sg, pn, cz, res$H4, res$H3, res$nsnp, res6$H4)
    }
  }
}
# setup check from the original rule (secondary): P-M marginal coloc should reproduce the registered 0.996
pm <- coloc_pair(TR$b_P[p1], TR$s_P[p1], N_P, TR$b_arg[p1], TR$s_arg[p1], N_ARG, "quant", NA, H$snp[p1], 1e-5)
pm_all <- coloc_pair(TR$b_P, TR$s_P, N_P, TR$b_arg, TR$s_arg, N_ARG, "quant", NA, H$snp, 1e-5)
addr("S4.check.PM_marginal.PPH4.p1np", "R1_ARG1_rule secondary (setup control)", "marginal PP.H4 ARG1 protein x arginine", "UKB-PPP ARG1 / CLSA arginine", "503 SNP set (p1np)",
     value = pm$H4, n = pm$nsnp, unit = "PP.H4", notes = sprintf("all 3,067 harmonised SNPs: %.4f (nsnp %d)", pm_all$H4, pm_all$nsnp))
assert_target("rule_PM_marginal_H4", pm_all$H4, 0.996, 0.0005,
              "coloc.abf on all 3,067 harmonised SNPs (P quant sdY 1, arginine quant sdY 1, p12 1e-5); the rule's registered value may come from the D1 SNP set (weak evidence)",
              "R1_ARG1_rule")

# ===================== S7. branches =====================
logf("-- S7 branches")
QH <- "nominal: PENDING_HOLM"
LDQ <- "sensitive to the reference LD"
SIGNQ <- "the sign of the correlation between the signals is not identified"
coloc_fires <- function(tn, sg) {
  a <- COL[[paste(tn, sg, "503")]]; b <- COL[[paste(tn, sg, "489")]]
  assoc503 <- abs(a$z) >= 1.96; assoc489 <- abs(b$z) >= 1.96
  fires <- assoc503 && (a$H4 < 0.80 || b$H4 < 0.80)
  list(fires = fires, assoc503 = assoc503, assoc489 = assoc489, a = a, b = b)
}
for (tn in names(TRAITS)) for (sg in c("A", "B")) {
  f <- coloc_fires(tn, sg)
  inputs <- sprintf("z 503 %.3f / 489 %.3f; PP.H4 503 %.4f / 489 %.4f (p12 1e-5)", f$a$z, f$b$z, f$a$H4, f$b$H4)
  br <- if (!f$assoc503) "does not fire: no association at the signal (|z| < 1.96 on 503)" else if (f$fires) "FIRES: associated (|z| >= 1.96) and PP.H4 < 0.80 on at least one panel" else
    "does not fire: associated and PP.H4 >= 0.80 on both panels"
  if (f$assoc503 != f$assoc489) br <- paste0(br, " [NOTE: |z| >= 1.96 differs between panels: 503 ", f$assoc503, ", 489 ", f$assoc489, "]")
  filled <- if (f$fires) sprintf("the %s association at signal %s did not colocalize with the protein (PP.H4 %s), so the ratio there is not anchored to the protein signal",
                                 tn, sg, if (f$a$H4 < 0.8 && f$b$H4 < 0.8) sprintf("%.2f on 503, %.2f on 489", f$a$H4, f$b$H4) else if (f$a$H4 < 0.8) sprintf("%.2f on 503", f$a$H4) else sprintf("%.2f on 489", f$b$H4)) else ""
  addb(sprintf("Rule 3(iv) %s signal %s", tn, sg), inputs, br, filled, if (tn == "T2D") "resource FinnGen R10" else "")
}
clause_for <- function(tn) { x <- character(0); for (sg in c("A", "B")) { f <- coloc_fires(tn, sg)
  if (f$fires) x <- c(x, sprintf("Rule 3(iv) clause: the %s association at signal %s did not colocalize with the protein (PP.H4 %.2f on 503, %.2f on 489), so the ratio there is not anchored to the protein signal",
                                 tn, sg, f$a$H4, f$b$H4)) }; x }

# Rule 3(ii) per primary cell
r3ii <- function(tag, label) {
  b <- BOOT[[tag]]
  addb(sprintf("Rule 3(ii) %s", label), sprintf("503 p %s (A-B %s); 489 p %s (A-B %s); bootstrap share same conclusion %.4f (side-only %.4f) of %d draws",
                                              fmtp(b$full$p), g7(b$full$d), fmtp(b$c489$p), g7(b$c489$d), b$share_same, b$share_side, sum(is.finite(rstar))),
       if (b$robust) "robust (same conclusion on 489 and in >= 90% of draws)" else "NOT robust",
       if (b$robust) "" else sprintf("add qualifier: %s", LDQ), if (b$robust) "" else LDQ)
  b$robust
}
rob_arg <- r3ii("arg.marg.aware2nd", "arginine A vs B (marginal, aware, 503)")
rob_t2d <- r3ii("FinnGen.cond.aware2nd", "T2D A vs B (FinnGen R10 conditional, aware)")
rob_dia_m <- r3ii("DIAMANTE.marg.aware2nd", "A2: T2D A vs B DIAMANTE marginal aware (conditional share printed beside)")
addb("Rule 3(ii) A2 beside", sprintf("DIAMANTE conditional aware: share same %.4f, robust %s", BOOT[["DIAMANTE.cond.aware2nd"]]$share_same, BOOT[["DIAMANTE.cond.aware2nd"]]$robust),
     "printed beside (not a primary cell)", "", "")
rob_t3c <- r3ii("t3.FinnGen.cond.aware2nd", "test (iii) conditional (FinnGen aware)")
rob_t3m <- r3ii("t3.FinnGen.marg.aware2nd", "test (iii) marginal (FinnGen aware)")
rob_orn <- r3ii("orn.marg.aware2nd", "ornithine A vs B (marginal, aware, 503)")

# Arginine
A5 <- ARGC[["marg 503"]]$aw2
if (pA_arg < 0.05 && pB_arg < 0.05) { op <- sprintf("both signals moved arginine (p_A %s, p_B %s)", fmtp(pA_arg), fmtp(pB_arg)); opb <- "both p < 0.05" } else {
  op <- paste(c(if (pA_arg >= 0.05) sprintf("signal A did not detectably move arginine (p %s)", fmtp(pA_arg)),
                if (pB_arg >= 0.05) sprintf("signal B did not detectably move arginine (p %s)", fmtp(pB_arg))), collapse = "; "); opb <- "a signal with p >= 0.05" }
addb("4 Arginine opening clause", sprintf("p_A %s, p_B %s (marginal arginine association at A, B)", fmtp(pA_arg), fmtp(pB_arg)), opb, op, "")
fi <- fieller(A5$RA, A5$RB, A5$vA, A5$vB, A5$cov)
qa <- c(QH, if (!rob_arg) LDQ, clause_for("arginine"))
if (A5$p >= 0.05) {
  addb("4 Arginine A vs B", sprintf("primary cell marginal aware 503: ratio A %s, ratio B %s, p %s; ratio of ratios %s (Fieller 95%% %s-%s)", g7(A5$RA), g7(A5$RB), fmtp(A5$p),
                                    fmt3(A5$RA / A5$RB), fmt3(fi$lo), fmt3(fi$hi)), "p >= 0.05",
       sprintf("no detectable difference per SD of ARG1 (ratio of ratios A/B %s, 95%% CI %s-%s; differences from %.0f%% to %.0f%% are not excluded)",
               fmt3(A5$RA / A5$RB), fmt3(fi$lo), fmt3(fi$hi), 100 * (fi$lo - 1), 100 * (fi$hi - 1)), paste(qa, collapse = "; "))
} else {
  addb("4 Arginine A vs B", sprintf("primary cell marginal aware 503: p %s; ratio of ratios %s (95%% %s-%s)", fmtp(A5$p), fmt3(A5$RA / A5$RB), fmt3(fi$lo), fmt3(fi$hi)), "p < 0.05",
       sprintf("the arginine associations per SD of ARG1 differed between the signals (ratio of ratios %s, 95%% CI %s-%s, p %s)", fmt3(A5$RA / A5$RB), fmt3(fi$lo), fmt3(fi$hi), fmtp(A5$p)),
       paste(qa, collapse = "; "))
}
# T2D
FC <- CELLS[["FinnGen|cond|503"]]
dirs <- sign(c(CELLS[["FinnGen|marg|503"]]$aw2$d, CELLS[["FinnGen|cond|503"]]$aw2$d, CELLS[["DIAMANTE|marg|503"]]$aw2$d, CELLS[["DIAMANTE|cond|503"]]$aw2$d))
samedir <- length(unique(dirs)) == 1
ins <- sprintf("FinnGen cond aware 503: A-B %s (95%% %s to %s), p %s; obs z_B %s, exp z_B %s; A-B signs in the 4 cells (FG marg, FG cond, DIA marg, DIA cond): %s; r bootstrap 95%% %s to %s",
               fmt4(FC$aw2$d), fmt4(FC$aw2$lo), fmt4(FC$aw2$hi), fmtp(FC$aw2$p), fmt3(FC$obs_zB), fmt3(FC$pw$exp_z), paste(dirs, collapse = ","), fmt4(rq[1]), fmt4(rq[2]))
qt <- c(QH, if (!rob_t2d) LDQ, if (R_BOOT_INCL0) SIGNQ, clause_for("T2D"))
if (FC$aw2$p < 0.05 && samedir) {
  verb <- if (abs(FC$obs_zB) < abs(FC$pw$exp_z)) "fell short of" else "exceeded"
  t2d_txt <- sprintf("signal B's T2D association (z %s) %s the z %s expected if T2D followed ARG1 in proportion (A-B difference %s, 95%% CI %s to %s; p %s)",
                     sprintf("%.2f", FC$obs_zB), verb, sprintf("%.2f", FC$pw$exp_z), fmt3(FC$aw2$d), fmt3(FC$aw2$lo), fmt3(FC$aw2$hi), fmtp(FC$aw2$p))
  t2d_br <- sprintf("primary p < 0.05, same direction in all four cells, |obs| %s |exp|", if (verb == "fell short of") "<" else ">")
} else if (FC$aw2$p < 0.05) {
  t2d_txt <- sprintf("the T2D associations differed between the signals in the primary cell (p %s), but the direction was not consistent across cells", fmtp(FC$aw2$p))
  t2d_br <- "primary p < 0.05, direction differs between cells"
} else {
  t2d_txt <- sprintf("the T2D associations did not detectably differ between the signals (A-B difference %s, 95%% CI %s to %s, p %s)", fmt3(FC$aw2$d), fmt3(FC$aw2$lo), fmt3(FC$aw2$hi), fmtp(FC$aw2$p))
  t2d_br <- "primary p >= 0.05"
}
if (R_BOOT_INCL0) t2d_txt <- paste0(t2d_txt, "; ", SIGNQ)
addb("4 T2D", ins, t2d_br, t2d_txt, paste(qt, collapse = "; "))
# Amendment 1 A2
Dm <- CELLS[["DIAMANTE|marg|503"]]$aw2
a2met <- FC$aw2$p < 0.05 && Dm$p < 0.05 && sign(FC$aw2$d) == sign(Dm$d)
fail_res <- c(if (!rob_t2d) "FinnGen", if (!rob_dia_m) "DIAMANTE")
addb("Amendment 1 A2", sprintf("FinnGen conditional aware p %s (A-B %s); DIAMANTE marginal aware p %s (A-B %s); Rule 3(ii) robust: FinnGen %s, DIAMANTE marginal %s",
                               fmtp(FC$aw2$p), fmt4(FC$aw2$d), fmtp(Dm$p), fmt4(Dm$d), rob_t2d, rob_dia_m),
     if (a2met) "met: both cells p < 0.05, same direction" else "not met: rev-4.2 template applies; DIAMANTE p printed beside",
     if (a2met) sprintf("in both outcome resources (FinnGen p %s; DIAMANTE p %s)", fmtp(FC$aw2$p), fmtp(Dm$p)) else sprintf("(DIAMANTE marginal p %s printed beside)", fmtp(Dm$p)),
     paste(c(if (a2met && length(fail_res)) sprintf("%s in %s", LDQ, paste(fail_res, collapse = " and ")), "never 'replicated'/'replication'",
             "family entry and Holm unchanged (FinnGen primary)"), collapse = "; "))
# Combined (test (iii))
T3c <- T3C[["FinnGen cond 503"]]$aw2; T3m <- T3C[["FinnGen marg 503"]]$aw2
beside <- sprintf("beside: DIAMANTE %s / %s; 489 %s / %s; beta-space %s / %s; first-order %s / %s (conditional / marginal)",
                  fmtp(T3C[["DIAMANTE cond 503"]]$aw2$p), fmtp(T3C[["DIAMANTE marg 503"]]$aw2$p), fmtp(T3C[["FinnGen cond 489"]]$aw2$p), fmtp(T3C[["FinnGen marg 489"]]$aw2$p),
                  fmtp(T3C[["FinnGen cond 503"]]$cj$p), fmtp(T3C[["FinnGen marg 503"]]$cj$p), fmtp(T3C[["FinnGen cond 503"]]$aw1$p), fmtp(T3C[["FinnGen marg 503"]]$aw1$p))
q3 <- c(QH, sprintf("family entry p = larger of the two estimands = %s", fmtp(max(T3c$p, T3m$p))), if (!(rob_t3c && rob_t3m)) LDQ, clause_for("T2D"), clause_for("arginine"))
if (T3c$p < 0.05 && T3m$p < 0.05) {
  addb("4 Combined (test (iii))", sprintf("FinnGen aware 2nd 503: conditional p %s, marginal p %s; %s", fmtp(T3c$p), fmtp(T3m$p), beside), "both estimands p < 0.05",
       sprintf("T2D and arginine did not scale together across the two signals (p %s conditional / %s marginal), which is suggestive evidence against ARG1 abundance acting through arginine in proportion; this is not a formal mediation test",
               fmtp(T3c$p), fmtp(T3m$p)), paste(c(q3, beside), collapse = "; "))
} else {
  addb("4 Combined (test (iii))", sprintf("FinnGen aware 2nd 503: conditional p %s, marginal p %s; %s", fmtp(T3c$p), fmtp(T3m$p), beside), "requirement not met (not both p < 0.05)",
       sprintf("a direct comparison of T2D and arginine across the two signals did not establish that they scale differently (p %s conditional / %s marginal)", fmtp(T3c$p), fmtp(T3m$p)),
       paste(c(q3, beside), collapse = "; "))
}
# Ornithine
O5 <- ORNC[["marg 503"]]$aw2
qo <- c(QH, if (!rob_orn) LDQ, clause_for("ornithine"))
if (O5$p < 0.05) {
  addb("4 Ornithine", sprintf("marginal aware 503: k_A %s, k_B %s, p %s", g7(O5$RA), g7(O5$RB), fmtp(O5$p)), "p < 0.05",
       sprintf("ornithine, the product of the arginase reaction, did not scale with ARG1 in the same proportion at the two signals (ornithine per SD of ARG1: %s at A, %s at B; p %s)",
               fmt3(O5$RA), fmt3(O5$RB), fmtp(O5$p)), paste(qo, collapse = "; "))
} else {
  addb("4 Ornithine", sprintf("marginal aware 503: k_A %s, k_B %s, p %s", g7(O5$RA), g7(O5$RB), fmtp(O5$p)), "p >= 0.05",
       sprintf("no detectable difference in ornithine per SD of ARG1 between the signals (p %s)", fmtp(O5$p)), paste(qo, collapse = "; "))
}
# Rule 4 family entries (for W9)
for (x in list(c("ARG1-arginine A vs B", A5$p), c("ARG1-T2D A vs B (FinnGen conditional aware)", FC$aw2$p),
               c("ARG1 test (iii) (larger of conditional and marginal)", max(T3c$p, T3m$p)), c("ARG1-ornithine A vs B", O5$p)))
  addb("Rule 4 family entry", sprintf("unadjusted p %s", fmtp(as.numeric(x[2]))), x[1], "Holm adjustment by W9", QH)

# ===================== S8. write outputs =====================
RT <- rbindlist(RR)
if (anyDuplicated(RT$quantity_id)) stop(sprintf("STOP: duplicate quantity ids: %s", paste(RT$quantity_id[duplicated(RT$quantity_id)], collapse = ",")))
na_val <- RT[value %in% c("NA", "", "NaN")]
logf("results: %d rows; rows with NA value: %d %s", nrow(RT), nrow(na_val), if (nrow(na_val)) paste(na_val$quantity_id, collapse = ",") else "")
fwrite(RT, file.path(OUT, "results.tsv"), sep = "\t", quote = FALSE)
BT2 <- rbindlist(BR)
fwrite(BT2, file.path(OUT, "branches.tsv"), sep = "\t", quote = FALSE)
logf("branches: %d rows", nrow(BT2))
inl <- do.call(rbind, INPUTS)
fwrite(as.data.table(inl), file.path(WORK, "inputs_md5.tsv"), sep = "\t")
for (i in seq_len(nrow(inl))) logf("INPUT %s | %s | %s | %s", inl[i, "id"], inl[i, "role"], inl[i, "md5"], inl[i, "path"])
logf("seen scripts opened: none")
SNAP1 <- snap("end")
same <- identical(SNAP0$file, SNAP1$file) && identical(SNAP0$md5, SNAP1$md5)
logf("staging snapshot start vs end: %s", if (same) "IDENTICAL (nothing written under V4_local/staging/l3-finemap)" else "DIFFERENT")
if (!same) stop("STOP: staging changed during the run")
logf("==== W4a run end; %.1f min ====", as.numeric(difftime(Sys.time(), T0, units = "mins")))
