# 05_susie.R -- section 6(d) Rule 3(i) proxy: per candidate, LD from PLINK (--r square bin, --keep-allele-order) on the 503 panel
# (sets susie, p1fr, p1np) and the 489 panel (p2np); one eigendecomposition per (window, set, panel) attached as attr(R, 'eigen');
# per trait estimate_s_rss + susie_rss (engine strict settings: L = 10, max_iter 1000, min_abs_corr 0.5, coverage 0.95,
# estimate_residual_variance FALSE); coloc.susie at p12 = 5e-6 (coloc 5.2.3 default) and 1e-5.
# Checkpoint: work/susie/<cid>/<set>_<panel>.rds; a finished unit is skipped on rerun. Progress: work/susie/progress.log.
# Usage: Rscript 05_susie.R [n_workers]     (driver)   |   sourced by workers (run_candidate)
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
SD <- file.path(WORK, "susie"); dir.create(SD, showWarnings = FALSE)
PROG <- file.path(SD, "progress.log")
prog <- function(fmt, ...) { s <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)); cat(s, "\n", file = PROG, append = TRUE); s }
SETS <- list(susie = c("in_all", "503"), p1fr = c("in_p1fr", "503"), p1np = c("in_p1np", "503"), p2np = c("in_p2np", "489"))

# PLINK LD with the declared assertions: every requested SNP returned, square, finite, unit diagonal, symmetric.
plink_ld <- function(bfile, ids, out_prefix) {
  idf <- paste0(out_prefix, ".ids"); writeLines(ids, idf)
  run_plink(c("--bfile", q(bfile), "--extract", q(idf), "--keep-allele-order", "--r", "square", "bin", "--write-snplist"), out_prefix)
  got <- readLines(paste0(out_prefix, ".snplist"))
  if (length(got) != length(ids) || !setequal(got, ids)) stop(sprintf("LD ASSERT: %d requested, %d returned, %d missing (e.g. %s)", length(ids), length(got),
                                                                     length(setdiff(ids, got)), paste(head(setdiff(ids, got), 3), collapse = ",")))
  p <- length(got); fb <- paste0(out_prefix, ".ld.bin")
  if (file.size(fb) != 8 * p * p) stop("LD ASSERT: .ld.bin size ", file.size(fb), " != 8*p^2")
  R <- matrix(readBin(fb, "double", n = p * p, size = 8), p, p)
  if (any(!is.finite(R))) stop("LD ASSERT: non-finite cells: ", sum(!is.finite(R)))
  if (max(abs(diag(R) - 1)) > 1e-6) stop("LD ASSERT: diagonal != 1")
  if (max(abs(R - t(R))) > 1e-9) stop("LD ASSERT: not symmetric")
  dimnames(R) <- list(got, got)
  o <- match(ids, got)
  R <- R[o, o]
  unlink(fb)                                   # matrix kept in memory only (re-made on rerun); its md5 is logged
  list(R = R, md5 = digest::digest(R, algo = "md5"))
}

fit_trait <- function(z, R, n) {
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(susieR::susie_rss(z = z, R = R, n = n, L = 10, max_iter = 1000, min_abs_corr = 0.5, coverage = 0.95,
                               estimate_residual_variance = FALSE),
             error = function(e) { w <<- c(w, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") })
  if (!is.null(f)) { colnames(f$lbf_variable) <- names(z); names(f$pip) <- names(z) }
  list(fit = f, warn = w)
}

run_candidate <- function(cid) {
  suppressMessages({library(susieR); library(coloc)})
  cdir <- file.path(SD, cid); dir.create(cdir, showWarnings = FALSE)
  H <- fread(file.path(WORK, "harm", paste0(cid, ".tsv.gz")))
  st <- fromJSON(file.path(WORK, "harm", paste0(cid, "__stamp.json")))
  H[, in_all := TRUE]
  todo <- names(SETS)[!file.exists(file.path(cdir, paste0(names(SETS), "_", vapply(SETS, `[`, "", 2), ".rds")))]
  if (!length(todo)) return(prog("%s all sets done (skip)", cid))
  prog("%s start: %d SNPs; sets to do: %s", cid, nrow(H), paste(todo, collapse = ","))
  LD <- list()
  if (any(todo %in% c("susie", "p1fr", "p1np"))) {
    t0 <- Sys.time()
    LD$`503` <- plink_ld(sub("\\.bed$", "", res_entry("LD503_bed")$path), H$rsid_ref, file.path(cdir, "ld503"))
    prog("%s LD 503 %dx%d from PLINK in %.0fs (md5 %s)", cid, nrow(H), nrow(H), as.numeric(difftime(Sys.time(), t0, units = "secs")), LD$`503`$md5)
  }
  if ("p2np" %in% todo && any(H$in_p2np)) {
    t0 <- Sys.time()
    h2 <- H[in_p2np == TRUE]
    pre <- file.path(ROOT, sprintf("lhcmr_eas_refs/1000G_EUR_Phase3_plink/1000G.EUR.QC.%s", h2$chr[1]))
    L2 <- plink_ld(pre, h2$p2_rsid, file.path(cdir, "ld489"))
    sg <- h2$p2_sign
    L2$R <- L2$R * outer(sg, sg)                # panel-2 A1 -> harmonised a1
    dimnames(L2$R) <- list(h2$snp, h2$snp)
    LD$`489` <- L2
    prog("%s LD 489 %dx%d from PLINK in %.0fs (md5 %s)", cid, nrow(h2), nrow(h2), as.numeric(difftime(Sys.time(), t0, units = "secs")), L2$md5)
  }
  for (sname in todo) {
    col <- SETS[[sname]][1]; panel <- SETS[[sname]][2]
    hs <- H[get(col) == TRUE]
    t0 <- Sys.time()
    if (!nrow(hs)) { saveRDS(list(cid = cid, set = sname, panel = panel, nsnp = 0, note = "empty set"), file.path(cdir, paste0(sname, "_", panel, ".rds"))); next }
    R <- if (panel == "503") LD$`503`$R[hs$rsid_ref, hs$rsid_ref, drop = FALSE] else LD$`489`$R[hs$snp, hs$snp, drop = FALSE]
    dimnames(R) <- list(hs$snp, hs$snp)
    ekey <- paste(panel, digest::digest(hs$snp))                 # identical SNP set on the same panel -> same R -> reuse its eigen
    if (!exists("EIG", inherits = FALSE)) EIG <- list()
    if (!is.null(EIG[[ekey]])) { attr(R, "eigen") <- EIG[[ekey]]; te <- 0; prog("%s %s/%s: SNP set identical to an earlier set on this panel; eigen reused", cid, sname, panel) } else {
      te <- system.time(attr(R, "eigen") <- eigen(R, symmetric = TRUE))[["elapsed"]]
      EIG <- list(); EIG[[ekey]] <- attr(R, "eigen")                # keep only the latest decomposition (memory)
    }
    traits <- list(P = list(z = hs$z_P, n = st$n_P), O = list(z = hs$z_O, n = st$n_O))
    trows <- list(); fits <- list(); csr <- list()
    for (tn in names(traits)) {
      z <- traits[[tn]]$z; names(z) <- hs$snp; n <- traits[[tn]]$n
      if (all(is.na(z))) { trows[[tn]] <- data.table(trait = tn, computable = FALSE, note = "no data for this trait"); next }
      if (anyNA(z)) stop(cid, " ", sname, " ", tn, ": NA z inside a set")
      s_hat <- tryCatch(susieR::estimate_s_rss(z, R, n), error = function(e) NA_real_)
      r <- fit_trait(z, R, n); f <- r$fit; fits[[tn]] <- f
      ncs <- if (is.null(f) || is.null(f$sets$cs)) 0L else length(f$sets$cs)
      conv <- !is.null(f) && isTRUE(f$converged)
      trows[[tn]] <- data.table(trait = tn, computable = TRUE, n = n, nsnp = length(z), max_abs_z = max(abs(z)), top_z_snp = names(z)[which.max(abs(z))],
                                s_hat = s_hat, converged = conv, niter = if (is.null(f)) NA_integer_ else f$niter, n_cs = ncs,
                                n_warn = length(r$warn), warnings = paste(unique(r$warn), collapse = " | "), evaluable = conv && !length(r$warn) && ncs >= 1)
      if (ncs) csr[[tn]] <- rbindlist(lapply(seq_len(ncs), function(i) { ix <- f$sets$cs[[i]]; ld <- ix[which.max(f$pip[ix])]
        data.table(trait = tn, cs = names(f$sets$cs)[i], size = length(ix), purity_min = f$sets$purity[i, "min.abs.corr"], lead_snp = names(z)[ld],
                   lead_pip = f$pip[ld], lead_z = z[ld], snps = paste(names(z)[ix], collapse = ";")) }))
    }
    pair <- data.table(n_cs_pairs = 0L, susie_max_H4 = NA_real_, susie_H3_at_max = NA_real_, susie_max_H4_p12_1e5 = NA_real_, hit_P = NA_character_, hit_O = NA_character_)
    if (!is.null(fits$P) && !is.null(fits$O) && length(fits$P$sets$cs) && length(fits$O$sets$cs)) {
      cs5 <- tryCatch(coloc.susie(fits$P, fits$O, p12 = 5e-6)$summary, error = function(e) NULL)
      cs1 <- tryCatch(coloc.susie(fits$P, fits$O, p12 = 1e-5)$summary, error = function(e) NULL)
      if (!is.null(cs5) && nrow(cs5)) { b <- cs5[which.max(cs5$PP.H4.abf)]
        pair <- data.table(n_cs_pairs = nrow(cs5), susie_max_H4 = b$PP.H4.abf, susie_H3_at_max = b$PP.H3.abf,
                           susie_max_H4_p12_1e5 = if (!is.null(cs1) && nrow(cs1)) max(cs1$PP.H4.abf) else NA_real_, hit_P = b$hit1, hit_O = b$hit2) }
    }
    res_u <- list(cid = cid, set = sname, panel = panel, nsnp = nrow(hs), eigen_s = te, ld_md5 = LD[[panel]]$md5,
                  traits = rbindlist(trows, fill = TRUE), cs = rbindlist(csr, fill = TRUE), pair = pair,
                  seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    saveRDS(res_u, file.path(cdir, paste0(sname, "_", panel, ".rds")))
    tr <- res_u$traits
    prog("%s %s/%s: %d SNPs, eigen %.0fs, total %.0fs | s_P %s s_O %s | conv %s | n_cs %s | coloc.susie maxH4 %s",
         cid, sname, panel, nrow(hs), te, res_u$seconds, format(tr[trait == "P"]$s_hat, digits = 4), format(tr[trait == "O"]$s_hat, digits = 4),
         paste(tr$converged, collapse = "/"), paste(tr$n_cs, collapse = "/"), format(pair$susie_max_H4, digits = 4))
  }
  prog("%s done", cid)
}

if (sys.nframe() == 0L) {
  say("==== W6b 05_susie driver ====")
  a <- commandArgs(trailingOnly = TRUE)
  nw <- if (length(a)) as.integer(a[1]) else 6L
  only <- if (length(a) > 1) a[-1] else character(0)
  cand <- fread(file.path(WORK, "candidates_windows.tsv"))
  have <- file.exists(file.path(WORK, "harm", paste0(cand$cid, "__stamp.json")))
  if (any(!have)) say("05 driver: not yet harmonised (skipped this run): %s", paste(cand$cid[!have], collapse = ","))
  cand <- cand[have]
  sz <- vapply(cand$cid, function(c) fromJSON(file.path(WORK, "harm", paste0(c, "__stamp.json")))$n_susie, 1)
  ord <- cand$cid[order(-sz)]
  if (length(only)) ord <- intersect(ord, only)
  say("05 driver: %d candidates, %d workers, order by size: %s", length(ord), nw, paste(sprintf("%s(%d)", ord, sz[ord]), collapse = " "))
  if (nw <= 1) { for (c in ord) run_candidate(c) } else {
    library(parallel)
    cl <- makeCluster(nw)
    on.exit(stopCluster(cl))
    clusterCall(cl, function() { source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "05_susie.R")); NULL })
    r <- clusterApplyLB(cl, ord, function(c) tryCatch(run_candidate(c), error = function(e) { prog("%s ERROR %s", c, conditionMessage(e)); paste("ERROR", c, conditionMessage(e)) }))
    for (x in r) say("05 worker: %s", x)
  }
  say("05 driver done")
}
