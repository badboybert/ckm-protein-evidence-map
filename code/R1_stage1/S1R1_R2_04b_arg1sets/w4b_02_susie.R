# w4b_02_susie.R -- Stage 1 W4b, step 2: SuSiE-RSS (engine settings) and Rule 3(i) s for ARG1 per trait x SNP set x panel.
# DEVIATION: engine path-locked -> the l3fm__susie.R settings and the prep QC are re-implemented here (no engine run).
# Settings (ENG_l3fm__susie.R, read for method): susie_rss L = 10, max_iter = 1000, min_abs_corr = 0.5, coverage = 0.95,
# estimate_residual_variance = FALSE; s = susieR::estimate_s_rss(z, R, n) (null-mle); one eigendecomposition per
# (set, SNP subset) attached as attr(R, "eigen"). Checkpoint per trait x set (susie/<trait>__<set>.rds); resumable.
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets/w4b_common.R"))
suppressMessages(library(susieR))
logf("==== W4b step 2 (SuSiE, s) start; susieR %s", as.character(packageVersion("susieR")))
D <- readRDS(file.path(OUT, "data", "w4b_data.rds")); H <- D$H
dir.create(file.path(OUT, "susie"), showWarnings = FALSE)
TRAITS <- list(P__ARG1 = list(z = "z__P__ARG1", n = D$N$P__ARG1, label = "ARG1 protein (UKB-PPP OID21445)"),
               arginine = list(z = "z__M__GCST90200372", n = D$N$M__GCST90200372, label = "arginine (CLSA GCST90200372)"),
               ornithine = list(z = "z__M__GCST90200404", n = D$N$M__GCST90200404, label = "ornithine (CLSA GCST90200404)"),
               T2D_FG = list(z = "z__FG", n = D$N$FG, label = "T2D (FinnGen R10, N_eff)"))
EIG <- list()
get_R <- function(set, snps) {
  key <- paste(set, length(snps), md5sum_chr(snps))
  Rf <- if (D$PANEL[[set]] == "503") D$R503 else D$R489
  R <- Rf[snps, snps, drop = FALSE]
  if (is.null(EIG[[key]])) {
    t0 <- Sys.time(); EIG[[key]] <<- eigen(R, symmetric = TRUE)
    logf("  eigen(R) set %s (%d SNPs, panel %s) %.0fs; min eigenvalue %.3g", set, length(snps), D$PANEL[[set]],
         as.numeric(difftime(Sys.time(), t0, units = "secs")), min(EIG[[key]]$values))
  }
  attr(R, "eigen") <- EIG[[key]]
  R
}
md5sum_chr <- function(x) { f <- file.path(OUT, "susie", ".snps.tmp"); writeLines(x, f); m <- md5(f); unlink(f); m }

fit_one <- function(z, R, n, L = 10L) {
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(susie_rss(z = z, R = R, n = n, L = L, max_iter = 1000, min_abs_corr = 0.5, coverage = 0.95,
                       estimate_residual_variance = FALSE),
             error = function(e) { w <<- c(w, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") })
  list(fit = f, warn = w)
}

units <- expand.grid(trait = names(TRAITS), set = names(D$SETS), stringsAsFactors = FALSE)
k <- 0L
for (i in seq_len(nrow(units))) {
  tr <- units$trait[i]; st <- units$set[i]; k <- k + 1L
  ck <- file.path(OUT, "susie", sprintf("%s__%s.rds", tr, st))
  if (file.exists(ck)) { logf("[%d/%d] %s x %s: checkpoint exists (skip)", k, nrow(units), tr, st); next }
  t0 <- Sys.time()
  snps <- D$SETS[[st]]
  z <- H[match(snps, snp)][[TRAITS[[tr]]$z]]
  ok <- !is.na(z); snps <- snps[ok]; z <- z[ok]; names(z) <- snps
  R <- get_R(st, snps)
  s_hat <- estimate_s_rss(z, R, TRAITS[[tr]]$n)
  r <- fit_one(z, R, TRAITS[[tr]]$n)
  f <- r$fit
  cs <- NULL
  if (!is.null(f) && length(f$sets$cs)) {
    cs <- rbindlist(lapply(seq_along(f$sets$cs), function(j) {
      ix <- f$sets$cs[[j]]; lead <- ix[which.max(f$pip[ix])]
      data.table(trait = tr, set = st, cs = names(f$sets$cs)[j], size = length(ix), coverage = f$sets$coverage[j],
                 purity_min = f$sets$purity[j, "min.abs.corr"], lead_snp = snps[lead], lead_pip = f$pip[lead], lead_z = z[lead],
                 members = paste(snps[ix], collapse = ";"), pips = paste(signif(f$pip[ix], 4), collapse = ";"))
    }))
  }
  out <- list(trait = tr, set = st, panel = D$PANEL[[st]], nsnp = length(snps), n = TRAITS[[tr]]$n, s = s_hat,
              converged = !is.null(f) && isTRUE(f$converged), niter = if (is.null(f)) NA else f$niter,
              n_cs = if (is.null(cs)) 0L else nrow(cs), warnings = paste(unique(r$warn), collapse = " | "),
              evaluable = !is.null(f) && isTRUE(f$converged) && length(r$warn) == 0 && !is.null(cs),
              cs = cs, pip = if (is.null(f)) NULL else setNames(f$pip, snps), max_abs_z = max(abs(z)), top_z_snp = snps[which.max(abs(z))])
  saveRDS(out, paste0(ck, ".part")); file.rename(paste0(ck, ".part"), ck)
  logf("[%d/%d] %-9s x %-5s (panel %s): %d SNPs n=%d s=%.4f conv=%s niter=%s n_cs=%d warn=%d |z|max %.2f (%s) %.0fs", k, nrow(units),
       tr, st, out$panel, out$nsnp, out$n, s_hat, out$converged, out$niter, out$n_cs, length(r$warn), out$max_abs_z, out$top_z_snp,
       as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
# collect
S <- rbindlist(lapply(seq_len(nrow(units)), function(i) {
  o <- readRDS(file.path(OUT, "susie", sprintf("%s__%s.rds", units$trait[i], units$set[i])))
  data.table(trait = o$trait, set = o$set, panel = o$panel, nsnp = o$nsnp, n = o$n, s = o$s, converged = o$converged,
             niter = o$niter, n_cs = o$n_cs, evaluable = o$evaluable, warnings = o$warnings)
}))
stopifnot(nrow(S) == 16L, !anyNA(S$s))
fwrite(S, file.path(OUT, "susie", "susie_s_table.tsv"), sep = "\t")
CS <- rbindlist(lapply(seq_len(nrow(units)), function(i) readRDS(file.path(OUT, "susie", sprintf("%s__%s.rds", units$trait[i], units$set[i])))$cs), fill = TRUE)
fwrite(CS, file.path(OUT, "susie", "susie_credible_sets.tsv"), sep = "\t")
logf("s table: %s", paste(sprintf("%s/%s=%.4f", S$trait, S$set, S$s), collapse = "; "))
logf("==== W4b step 2 done")
