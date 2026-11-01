# S5_common.R -- shared helpers for W5 (pre-declaration rev 4.2 §5, DPEP1). Sourced by every S5_*.R script.
# Writes only inside OUT (this folder). Every input is md5-asserted at read time.
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
suppressMessages({library(data.table); library(jsonlite)})

ROOT <- .P31R_CKM
OUT  <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_05_dpep1")
WORK <- file.path(OUT, "work")
LDD  <- file.path(OUT, "ld")
dir.create(WORK, showWarnings = FALSE); dir.create(LDD, showWarnings = FALSE)
LOG  <- file.path(OUT, "run.log")

say <- function(...) {
  s <- sprintf(...)
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG, append = TRUE)
  invisible(s)
}

md5 <- function(p) unname(tools::md5sum(p))
assert_md5 <- function(path, expect, label = basename(path)) {
  if (!file.exists(path)) stop("MISSING input ", label, ": ", path)
  got <- md5(path)
  if (!identical(got, expect)) stop(sprintf("MD5 MISMATCH %s: got %s expected %s (%s)", label, got, expect, path))
  say("md5 OK  %s  %s  (%s bytes)  %s", expect, label, format(file.size(path), big.mark = ","), path)
  invisible(path)
}

# ---- governing documents (md5s from the workflow prompt) ----
GOV <- list(
  BRIEFS        = list(p = file.path(ROOT, "V3/R1_stage1/S1R1_R2_00_resources/STAGE1_BRIEFS.md"), m = "40d0d52f7055542d774c03f1edf04e00"),
  PRESPEC_RED   = list(p = file.path(ROOT, "V3/R1_stage1/S1R1_R2_00_resources/P31R_R2_STAGE1_PRESPEC_REDACTED.md"), m = "a281c3464619b1b67c28e13e7392c41c"),
  AMEND1_RED    = list(p = file.path(ROOT, "V3/R1_stage1/S1R1_R2_00_resources/P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), m = "f265d600973c43239bb8377e0cd39d33"),
  PRESPEC_FULL  = list(p = file.path(ROOT, "V3/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), m = "3319dc39740d9f5b3b8822133d9108b7"),
  AMEND1_FULL   = list(p = file.path(ROOT, "V3/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), m = "377508832ea9bd384133e869c456c909"),
  GLS_Q         = list(p = file.path(ROOT, "V3/R1_stage1/S1R1_R2_00_resources/gls_q.R"), m = "b7390fdf5b6f2f9c1d30c298d6afdafc"),
  RESOURCES     = list(p = file.path(ROOT, "V3/R1_stage1/S1R1_R2_00_resources/RESOURCES.json"), m = "6d21238155a62fb58f71e05a01638dd0"),
  OUTCOME_PRESPEC = list(p = file.path(ROOT, "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), m = "c9079f2de1ffe75ce31f8a1d597d1fda")
)
assert_governing <- function() for (k in names(GOV)) assert_md5(GOV[[k]]$p, GOV[[k]]$m, paste0("GOV_", k))

# ---- inputs: path + md5 taken from RESOURCES.json by id (asserted against the file at read time) ----
.res <- NULL
res <- function(id) {
  if (is.null(.res)) .res <<- fromJSON(GOV$RESOURCES$p, simplifyVector = FALSE)
  e <- Filter(function(x) identical(x$id, id), .res)
  if (length(e) != 1) stop("RESOURCES id not unique / not found: ", id)
  e <- e[[1]]
  if (!identical(e$status, "FOUND")) stop("RESOURCES id not FOUND: ", id)
  list(id = id, path = e$path, md5 = e$md5, role = e$role)
}
inp <- function(id) { r <- res(id); if (identical(r$role, "seen_output")) stop("seen_output used as input: ", id); assert_md5(r$path, r$md5, id); r$path }

# Inputs the resource map lacks (registered in README with path + md5 + reason): the other 7 trait windows at the
# DPEP1 locus, needed only to rebuild the declared p1fr SNP-set rule (every reporting GWAS frequency); md5 = PREP_STAMP.
EXTRA <- list(
  P__CDH15        = list(p = file.path(ROOT, "weight-loss.omics/data/proteome/pqtl_regional/OID21005_chr16_cis.tsv"), m = "2dea4d96f8f858820d0910c5ee52ce53", kind = "protein"),
  M__GCST90199883 = list(p = file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90199883__chr16_89113308_90138456.tsv.gz"), m = "b683648f00480b241fcd52163d26a720", kind = "metabolite"),
  M__GCST90200319 = list(p = file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90200319__chr16_89113308_90138456.tsv.gz"), m = "f838630638ed5f630cec9aa381506aa6", kind = "metabolite"),
  M__GCST90199922 = list(p = file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90199922__chr16_89113308_90138456.tsv.gz"), m = "95133b76803ca4d76eb49a8bbcda7328", kind = "metabolite"),
  M__GCST90200282 = list(p = file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90200282__chr16_89113308_90138456.tsv.gz"), m = "9cdb62380eea55de782a8b73e68afd7a", kind = "metabolite"),
  M__GCST90199961 = list(p = file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90199961__chr16_89113308_90138456.tsv.gz"), m = "a5310ca433e8fa48536c805a7846b40d", kind = "metabolite"),
  M__GCST90200234 = list(p = file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90200234__chr16_89113308_90138456.tsv.gz"), m = "8f238cf2388c5b097e510f3f34f146d9", kind = "metabolite")
)

# ---- constants fixed by the pre-declaration / briefs ----
SEED <- 20261007L
N_BOOT <- 5000L
SIG5 <- c("rs409170", "rs1126464", "rs58290281", "rs117990130", "rs72811595")   # pre-declaration §5 (signals.tsv order)
SIG_A <- "rs409170"; SIG_B <- "rs62068372"                                       # §5 two-signal re-test
CHMP <- c("rs460879", "rs34427386")                                               # §5 attribution check
P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5                                               # Rule 3(iv) prior (confirmed vs COND_STAMP engine)
COLLIN <- 0.9                                                                     # engine's conditioning collinearity guard
WIN38 <- c(89113308L, 90138456L)                                                  # PREP_STAMP window_grch38 (definition)
PLINK <- NULL   # set by S5_01 from RESOURCES PLINK19_paper6

# ---- PLINK wrapper (Windows path with spaces) ----
run_plink <- function(args, tag) {
  exe <- res("PLINK19_paper6")$path
  out <- suppressWarnings(system2(exe, args = args, stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); if (is.null(st)) st <- 0L
  writeLines(out, file.path(LDD, paste0(tag, ".plinklog.txt")))
  say("PLINK %s: exit %d; %s", tag, st, paste(grep("variants|people|Error|Warning|written", out, value = TRUE), collapse = " | "))
  invisible(list(status = st, out = out))
}
q <- function(x) shQuote(x, type = "cmd")

# ---- z-space conditioning (Yang 2012 approximation; = l3fm__conditional.condition and = R1 single-SNP conditioning) ----
# z: trait z over the target SNP set; RiS: r between target SNPs and the conditioning SNPs (p x k); RSS: k x k; zS: trait z at S.
# isS: logical, target SNPs that ARE conditioning SNPs (set NaN). Collinear (multiple r^2 > COLLIN) -> NaN, never inflated.
condition_z <- function(z, se, RiS, RSS, zS, isS = rep(FALSE, length(z))) {
  RiS <- as.matrix(RiS); RSS <- as.matrix(RSS)
  W <- solve(RSS, t(RiS))                       # k x p
  qq <- rowSums(RiS * t(W))                     # r_iS RSS^-1 r_Si
  den <- sqrt(pmax(1 - qq, 1e-12))
  zc <- (z - as.vector(RiS %*% solve(RSS, zS))) / den
  sc <- se / den
  bad <- qq > COLLIN | isS
  zc[bad] <- NA_real_; sc[bad] <- NA_real_
  list(z = zc, se = sc, q = qq, n_collinear = sum(qq > COLLIN & !isS))
}

# ---- coloc.abf (R coloc 5.2.3) on two conditioned datasets; NA rows dropped (engine: >= 20 SNPs required) ----
coloc_pair <- function(b1, s1, b2, s2, snp, ds2_extra = NULL, ds1_extra = NULL) {
  ok <- is.finite(b1) & is.finite(b2) & is.finite(s1) & is.finite(s2)
  if (sum(ok) < 20) return(list(H4 = NA_real_, H3 = NA_real_, nsnp = sum(ok), lead = NA_character_))
  d1 <- c(list(beta = b1[ok], varbeta = s1[ok]^2, snp = snp[ok], type = "quant"), if (is.null(ds1_extra)) list(sdY = 1) else lapply(ds1_extra, function(v) if (length(v) == length(ok)) v[ok] else v))
  d2 <- c(list(beta = b2[ok], varbeta = s2[ok]^2, snp = snp[ok], type = "quant"), if (is.null(ds2_extra)) list(sdY = 1) else lapply(ds2_extra, function(v) if (length(v) == length(ok)) v[ok] else v))
  r <- NULL
  invisible(utils::capture.output(r <- suppressMessages(suppressWarnings(coloc::coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12)))))
  list(H4 = unname(r$summary["PP.H4.abf"]), H3 = unname(r$summary["PP.H3.abf"]), nsnp = sum(ok),
       lead = r$results$snp[which.max(r$results$SNP.PP.H4)])
}

# ---- direct Wakefield ABF colocalization (same model as coloc.abf with type "quant", sdY = 1: prior sd 0.15, p1, p2, p12).
# Used only inside the 5,000-draw bootstrap for speed; S5_05 check C2 asserts it equals coloc::coloc.abf (|diff| < 1e-6).
abf_direct <- function(b1, s1, b2, s2) {
  ok <- is.finite(b1) & is.finite(b2) & is.finite(s1) & is.finite(s2)
  if (sum(ok) < 20) return(list(H4 = NA_real_, H3 = NA_real_, nsnp = sum(ok)))
  labf <- function(z, V, W = 0.15^2) { r <- W / (W + V); 0.5 * (log(1 - r) + r * z^2) }
  lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
  l1 <- labf(b1[ok] / s1[ok], s1[ok]^2); l2 <- labf(b2[ok] / s2[ok], s2[ok]^2)
  a <- lse(l1) + lse(l2); s12 <- lse(l1 + l2)
  lH <- c(0, log(P1) + lse(l1), log(P2) + lse(l2), log(P1) + log(P2) + a + log1p(-exp(s12 - a)), log(P12) + s12)
  pp <- exp(lH - lse(lH))
  list(H4 = pp[5], H3 = pp[4], nsnp = sum(ok))
}

# ---- two-signal A-vs-B test of the protein-scaled ratios (Rule 5): d = kA - kB, k = M/P on two-SNP joint estimates;
# delta method: aware keeps the LD-induced covariance of the joint estimates, ign keeps only their diagonal; 2nd order keeps
# the protein-side term (gradient d/dP), 1st order drops it (R1 implementation: ratio SE = |k| sqrt((sM/M)^2 + (sP/P)^2)).
ratio_diff <- function(jp, jm, cov = "aware", se_order = "2nd") {
  PA <- jp$b[1]; PB <- jp$b[2]; MA <- jm$b[1]; MB <- jm$b[2]
  VM <- jm$V; VP <- jp$V
  if (cov == "ign") { VM <- diag(diag(VM)); VP <- diag(diag(VP)) }
  if (se_order == "1st") VP <- VP * 0
  gM <- c(1 / PA, -1 / PB); gP <- c(-MA / PA^2, MB / PB^2)
  v <- as.numeric(t(gM) %*% VM %*% gM + t(gP) %*% VP %*% gP)
  kA <- MA / PA; kB <- MB / PB
  seA <- sqrt(VM[1, 1] / PA^2 + MA^2 * VP[1, 1] / PA^4); seB <- sqrt(VM[2, 2] / PB^2 + MB^2 * VP[2, 2] / PB^4)
  d <- kA - kB; z <- d / sqrt(v)
  list(kA = kA, kB = kB, seA = seA, seB = seB, d = d, se_d = sqrt(v), z = z, p = 2 * pnorm(-abs(z)), PA = PA, PB = PB, MA = MA, MB = MB)
}
# ---- staging snapshot: md5 of every file under V4_local/staging/l3-finemap (proves nothing was written there)
snap_staging <- function(tag) {
  sd <- file.path(ROOT, "V3/V4_local/staging/l3-finemap")
  f <- sort(list.files(sd, recursive = TRUE, full.names = TRUE, all.files = TRUE))
  m <- md5(f)
  dt <- data.table(file = sub(paste0(sd, "/"), "", f, fixed = TRUE), bytes = file.size(f), md5 = m)
  fwrite(dt, file.path(WORK, sprintf("staging_md5_%s.tsv", tag)), sep = "\t")
  dig <- digest::digest(paste(dt$file, dt$md5, collapse = "\n"), algo = "md5", serialize = FALSE)
  say("STAGING SNAPSHOT %s: %d files under V4_local/staging/l3-finemap; digest of (file,md5) list = %s", tag, nrow(dt), dig)
  invisible(dig)
}
fmt <- function(x, d = 6) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "g"))
