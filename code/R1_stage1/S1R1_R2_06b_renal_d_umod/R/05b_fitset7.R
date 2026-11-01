# 05b_fitset7.R -- Stage 2 step 2 (user decision 2026-10-09): s_max over ALL SEVEN distinct SuSiE fits at every locus.
# 05_susie.R fits four: susie, p1fr, p1np on the 503 panel and p2np on the 489 panel (= p1np on 489). Rule 3(i) reads "all
# variant sets, both panels", so three more fits are needed, per candidate, s only (estimate_s_rss; no SuSiE fit, no credible sets):
#   p2np on 503 (sub-matrix of the PLINK 503 LD), susie on 489, p1fr on 489 (both restricted to SNPs valid on the 489 panel).
# 489 validity = the existing p2np rule applied to EVERY harmonised SNP, palindromes included (04_harmonise.R tests only the
# non-palindromic ones): same GRCh37 position and allele pair (unique), no missing genotype (nchrobs == 2 x 489), panel-2 MAF >=
# 0.01, |AF2 - AF1| <= 0.05. Asserts per candidate: valid & non-palindromic reproduces in_p2np / p2_rsid / p2_sign exactly; and the
# same code path reproduces the stored p1np/503 and p2np/489 s (< 1e-6).
# Checkpoint: work/susie/<cid>/fitset7.rds (finished candidates are skipped). Progress: work/susie/progress.log.
# Usage: Rscript 05b_fitset7.R [n_workers] [cid ...]
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
source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "05_susie.R"))

read_ref7 <- function(pre) {
  bim <- fread(paste0(pre, ".bim"), header = FALSE, col.names = c("chr", "rsid", "cm", "pos37", "a1", "a2"), colClasses = list(character = c(1, 2, 5, 6)))
  frq <- fread(paste0(pre, ".frq.counts"), colClasses = list(character = c("CHR", "SNP", "A1", "A2")))
  stopifnot(nrow(bim) == nrow(frq), all(bim$rsid == frq$SNP), all(bim$a1 == frq$A1), all(bim$a2 == frq$A2))
  bim[, `:=`(af = frq$C1 / (frq$C1 + frq$C2), nchrobs = frq$C1 + frq$C2, a1 = toupper(a1), a2 = toupper(a2))]
  bim[]
}

run_candidate7 <- function(cid) {
  suppressMessages(library(susieR))
  cdir <- file.path(SD, cid); dir.create(cdir, showWarnings = FALSE)
  ck <- file.path(cdir, "fitset7.rds")
  if (file.exists(ck)) return(prog("%s fitset7 checkpoint exists (skip)", cid))
  t00 <- Sys.time()
  H <- fread(file.path(WORK, "harm", paste0(cid, ".tsv.gz")), colClasses = list(character = "chr"))
  st <- fromJSON(file.path(WORK, "harm", paste0(cid, "__stamp.json")))
  ch <- as.character(H$chr[1]); lo37 <- st$window37[1]; hi37 <- st$window37[2]
  # ---- 489 validity for every harmonised SNP ----
  b2 <- read_ref7(file.path(WORK, "ref", paste0("ref489_chr", ch)))[pos37 >= lo37 & pos37 <= hi37]
  j <- merge(H[, .(snp, pos37, a1, a2, af, in_p1np)], b2[, .(pos37, p2_rsid = rsid, p2_a1 = a1, p2_a2 = a2, af2_raw = af, nch2 = nchrobs)], by = "pos37", allow.cartesian = TRUE)
  j <- j[(p2_a1 == a1 & p2_a2 == a2) | (p2_a1 == a2 & p2_a2 == a1)]
  nj <- j[, .N, by = snp]
  j <- j[snp %in% nj[N == 1]$snp]
  j[, p2_sign := ifelse(p2_a1 == a1, 1L, -1L)]
  j[, af2 := ifelse(p2_sign == 1L, af2_raw, 1 - af2_raw)]
  j <- j[nch2 == 2L * 489L & pmin(af2, 1 - af2) >= REF_MAF_MIN & abs(af2 - af) <= P2_AF_DIFF_MAX]
  # the existing p2np set is reproduced exactly
  jn <- j[in_p1np == TRUE]
  hp <- H[in_p2np == TRUE]
  stopifnot(setequal(jn$snp, hp$snp), nrow(jn) == nrow(hp))
  mm <- merge(jn[, .(snp, p2_rsid, p2_sign)], hp[, .(snp, p2_rsid_s = p2_rsid, p2_sign_s = p2_sign)], by = "snp")
  stopifnot(nrow(mm) == nrow(hp), all(mm$p2_rsid == mm$p2_rsid_s), all(mm$p2_sign == mm$p2_sign_s))
  H[, valid489 := snp %in% j$snp]
  hv <- merge(H[valid489 == TRUE, .(snp, pos37)], j[, .(snp, p2_rsid, p2_sign)], by = "snp"); setorder(hv, pos37)
  prog("%s fitset7: %d SNPs; valid-489 %d (non-palindromic %d == in_p2np %d; palindromic %d)", cid, nrow(H), nrow(hv), nrow(jn), nrow(hp), nrow(hv) - nrow(jn))
  # ---- LD ----
  t0 <- Sys.time()
  L503 <- plink_ld(sub("\\.bed$", "", res_entry("LD503_bed")$path), H$rsid_ref, file.path(cdir, "ld503_f7"))
  R503 <- L503$R; dimnames(R503) <- list(H$snp, H$snp)
  pre <- file.path(ROOT, sprintf("lhcmr_eas_refs/1000G_EUR_Phase3_plink/1000G.EUR.QC.%s", ch))
  L489 <- plink_ld(pre, hv$p2_rsid, file.path(cdir, "ld489_f7"))
  R489 <- L489$R * outer(hv$p2_sign, hv$p2_sign); dimnames(R489) <- list(hv$snp, hv$snp)
  prog("%s fitset7 LD 503 %d x %d and 489 %d x %d in %.0fs", cid, nrow(R503), ncol(R503), nrow(R489), ncol(R489), as.numeric(difftime(Sys.time(), t0, units = "secs")))
  # ---- fits ----
  FITS <- list(
    list(fit = "p2np_503", set = "p2np", panel = "503", snps = hp$snp),
    list(fit = "susie_489", set = "susie", panel = "489", snps = hv$snp),
    list(fit = "p1fr_489", set = "p1fr", panel = "489", snps = intersect(H[in_p1fr == TRUE]$snp, hv$snp)),
    list(fit = "chk_p1np_503", set = "p1np", panel = "503", snps = H[in_p1np == TRUE]$snp),
    list(fit = "chk_p2np_489", set = "p2np", panel = "489", snps = hp$snp))
  Hs <- H[, .(snp, z_P, z_O)]
  rows <- list()
  for (F in FITS) {
    if (!length(F$snps)) { prog("%s %s: empty set", cid, F$fit); next }
    R <- if (F$panel == "503") R503[F$snps, F$snps, drop = FALSE] else R489[F$snps, F$snps, drop = FALSE]
    te <- system.time(attr(R, "eigen") <- eigen(R, symmetric = TRUE))[["elapsed"]]
    hh <- Hs[match(F$snps, snp)]
    for (tn in c("P", "O")) {
      z <- if (tn == "P") hh$z_P else hh$z_O; n <- if (tn == "P") st$n_P else st$n_O
      if (all(is.na(z))) { rows[[paste(F$fit, tn)]] <- data.table(fit = F$fit, set = F$set, panel = F$panel, trait = tn, computable = FALSE, n = NA_real_, nsnp = length(z), s_hat = NA_real_); next }
      if (anyNA(z)) stop(cid, " ", F$fit, " ", tn, ": NA z inside a set")
      names(z) <- F$snps
      s_hat <- susieR::estimate_s_rss(z, R, n)
      rows[[paste(F$fit, tn)]] <- data.table(fit = F$fit, set = F$set, panel = F$panel, trait = tn, computable = TRUE, n = n, nsnp = length(z), s_hat = s_hat)
      prog("%s fitset7 %-12s trait %s: %d SNPs n=%.0f s=%.8f (eigen %.0fs)", cid, F$fit, tn, length(z), n, s_hat, te)
    }
  }
  X <- rbindlist(rows)
  # ---- identity checks against the stored four-fit values ----
  for (cc in list(c("chk_p1np_503", "p1np_503"), c("chk_p2np_489", "p2np_489"))) {
    u <- readRDS(file.path(cdir, paste0(cc[2], ".rds")))
    for (tn in c("P", "O")) {
      a <- X[fit == cc[1] & trait == tn]; b <- u$traits[trait == tn]
      if (!nrow(a) || !a$computable) { stopifnot(!nrow(b) || !isTRUE(b$computable)); next }
      d <- abs(a$s_hat - b$s_hat); prog("%s CHECK %s trait %s vs stored: |diff| %.3g", cid, cc[1], tn, d)
      if (!(d < 1e-6)) stop(cid, ": identity check failed for ", cc[1], " ", tn, " (", d, ")")
    }
  }
  saveRDS(list(cid = cid, fits = X[!grepl("^chk_", fit)], idcheck = X[grepl("^chk_", fit)], n_valid489 = nrow(hv), n_valid489_pal = nrow(hv) - nrow(jn),
               ld503_md5 = L503$md5, ld489_md5 = L489$md5, seconds = as.numeric(difftime(Sys.time(), t00, units = "secs"))), paste0(ck, ".part"))
  file.rename(paste0(ck, ".part"), ck)
  unlink(c(file.path(cdir, "ld503_f7.ld.bin"), file.path(cdir, "ld489_f7.ld.bin")))
  prog("%s fitset7 done in %.0fs", cid, as.numeric(difftime(Sys.time(), t00, units = "secs")))
}

if (sys.nframe() == 0L) {
  say("==== W6b 05b_fitset7 driver ====")
  a <- commandArgs(trailingOnly = TRUE)
  nw <- if (length(a)) as.integer(a[1]) else 4L
  only <- if (length(a) > 1) a[-1] else character(0)
  cand <- fread(file.path(WORK, "candidates_windows.tsv"))
  stopifnot(all(file.exists(file.path(WORK, "harm", paste0(cand$cid, "__stamp.json")))))
  sz <- vapply(cand$cid, function(c) fromJSON(file.path(WORK, "harm", paste0(c, "__stamp.json")))$n_susie, 1)
  ord <- cand$cid[order(-sz)]
  if (length(only)) ord <- intersect(ord, only)
  say("05b driver: %d candidates, %d workers, order by size: %s", length(ord), nw, paste(sprintf("%s(%d)", ord, sz[ord]), collapse = " "))
  if (nw <= 1) { for (c in ord) run_candidate7(c) } else {
    library(parallel)
    cl <- makeCluster(nw)
    on.exit(stopCluster(cl))
    clusterCall(cl, function() { source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "05b_fitset7.R")); NULL })
    r <- clusterApplyLB(cl, ord, function(c) tryCatch(run_candidate7(c), error = function(e) { prog("%s ERROR %s", c, conditionMessage(e)); paste("ERROR", c, conditionMessage(e)) }))
    for (x in r) say("05b worker: %s", x)
    if (any(grepl("^ERROR", unlist(r)))) stop("05b: a worker failed")
  }
  miss <- cand$cid[!file.exists(file.path(SD, cand$cid, "fitset7.rds"))]
  if (length(miss)) stop("05b: no fitset7 checkpoint for ", paste(miss, collapse = ","))
  say("05b driver done")
}
