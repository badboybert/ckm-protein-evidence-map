# w2_04_susie_coloc.R -- TGFB1 SuSiE-RSS (four l3fm-style variant sets, both panels) and the Rule 3(iv) conditional
# colocalization (pre-declaration section 2; Rules 3(i), 3(iv)).
#   Rscript w2_04_susie_coloc.R <OUT dir>
# SuSiE settings = the l3fm engine (l3fm__susie.R): susie_rss(z, R, n, L = 10, max_iter = 1000, min_abs_corr = 0.5,
# coverage = 0.95, estimate_residual_variance = FALSE); s = susieR::estimate_s_rss(z, R, n) (method null-mle, default).
# One eigendecomposition per (set, panel), attached as attr(R, "eigen"); checkpoint per (set, trait) under susie/.
# Sets: susie, p1fr, p1np on the 503 PLINK matrix; p2np on the 489 PLINK matrix (A1-aligned to 503 by w2_02_locus.py).
# n: P = engine n_hat; CHD_FG and CAD_C4D = N_eff from case/control counts (see tgfb1_locus_meta.json).
# Conditional coloc: coloc.abf on the conditioned window (w2_03_cond.py), p1 = p2 = 1e-4, p12 = 1e-5;
# protein type "quant" sdY = 1; CAD type "cc" with the case fraction and total N (primary); engine convention (both
# "quant", sdY = 1) printed beside.
suppressMessages({library(data.table); library(susieR); library(coloc); library(jsonlite)})
a <- commandArgs(trailingOnly = TRUE)
OUT <- a[1]
logf <- file.path(OUT, "run.log")
lg <- function(...) { m <- sprintf("%s  %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
lg("==== w2_04_susie_coloc.R (R %s; susieR %s; coloc %s) ====", R.version.string, packageVersion("susieR"), packageVersion("coloc"))
meta <- fromJSON(file.path(OUT, "tgfb1_locus_meta.json"))
chk <- function(f, m) { got <- unname(tools::md5sum(f)); if (got != m) stop("md5 mismatch ", f); lg("md5 ok %s %s", basename(f), got) }
chk(file.path(OUT, "tgfb1_locus.tsv"), meta$md5$tgfb1_locus.tsv)
chk(file.path(OUT, "ld", "R503_susie.f8"), meta$md5$R503_susie.f8)
chk(file.path(OUT, "ld", "R489_p2np.f8"), meta$md5$R489_p2np.f8)
h <- fread(file = file.path(OUT, "tgfb1_locus.tsv"), sep = "\t")
if (nrow(h) != meta$n$susie) stop("tgfb1_locus.tsv read ", nrow(h), " rows, expected ", meta$n$susie)
p <- nrow(h)
R503 <- matrix(readBin(file.path(OUT, "ld", "R503_susie.f8"), "double", n = p * p, size = 8, endian = "little"), p, p)
i2 <- which(h$in_p2np == 1)
R489 <- matrix(readBin(file.path(OUT, "ld", "R489_p2np.f8"), "double", n = length(i2)^2, size = 8, endian = "little"), length(i2), length(i2))
stopifnot(isSymmetric(R503, tol = 1e-6), isSymmetric(R489, tol = 1e-6), all(abs(diag(R503) - 1) < 1e-6), all(abs(diag(R489) - 1) < 1e-6))
lg("LD read: 503 %d x %d, 489 %d x %d", p, p, length(i2), length(i2))
SIG <- c(A = "rs73045269", B = "rs12983775", ANN = "rs1800470")
traits <- c("P", "CHD_FG", "CAD_C4D")
tn <- unlist(meta$trait_n)
sets <- list(susie = list(idx = seq_len(p), panel = "503"), p1fr = list(idx = which(h$in_p1fr == 1), panel = "503"),
             p1np = list(idx = which(h$in_p1np == 1), panel = "503"), p2np = list(idx = i2, panel = "489"))
dir.create(file.path(OUT, "susie"), showWarnings = FALSE)
trows <- list(); csrows <- list()
for (sn in names(sets)) {
  S <- sets[[sn]]
  ids <- h$rsid503[S$idx]
  todo <- traits[!file.exists(file.path(OUT, "susie", sprintf("%s__%s.rds", sn, traits)))]
  if (length(todo)) {
    R <- if (S$panel == "503") R503[S$idx, S$idx] else R489
    dimnames(R) <- list(ids, ids)
    te <- system.time(attr(R, "eigen") <- eigen(R, symmetric = TRUE))[["elapsed"]]
    lg("[susie %s] %d SNPs, panel %s, eigen %.0fs; min eigenvalue %.3g", sn, length(ids), S$panel, te, min(attr(R, "eigen")$values))
    for (t in todo) {
      z <- h[[paste0("z__", t)]][S$idx]; names(z) <- ids
      n <- tn[[t]]
      s_hat <- tryCatch(estimate_s_rss(z, R, n), error = function(e) NA_real_)
      w <- character(0)
      f <- withCallingHandlers(
        tryCatch(susie_rss(z = z, R = R, n = n, L = 10, max_iter = 1000, min_abs_corr = 0.5, coverage = 0.95,
                           estimate_residual_variance = FALSE),
                 error = function(e) { w <<- c(w, paste("ERROR:", conditionMessage(e))); NULL }),
        warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") })
      cs <- NULL
      if (!is.null(f) && length(f$sets$cs)) {
        cs <- rbindlist(lapply(seq_along(f$sets$cs), function(i) {
          ix <- f$sets$cs[[i]]; lead <- ix[which.max(f$pip[ix])]
          data.table(cs = names(f$sets$cs)[i], size = length(ix), coverage = f$sets$coverage[i],
                     purity_min = f$sets$purity[i, "min.abs.corr"], lead = ids[lead], lead_pip = f$pip[lead], lead_z = z[lead],
                     has_A = SIG[["A"]] %in% ids[ix], has_B = SIG[["B"]] %in% ids[ix], has_rs1800470 = SIG[["ANN"]] %in% ids[ix],
                     pip_A = if (SIG[["A"]] %in% ids) f$pip[match(SIG[["A"]], ids)] else NA_real_,
                     pip_B = if (SIG[["B"]] %in% ids) f$pip[match(SIG[["B"]], ids)] else NA_real_,
                     members = paste(ids[ix], collapse = ";"))
        }))
      }
      res <- list(set = sn, panel = S$panel, trait = t, nsnp = length(ids), n = n, s = s_hat, converged = !is.null(f) && isTRUE(f$converged),
                  niter = if (is.null(f)) NA_integer_ else f$niter, n_cs = if (is.null(cs)) 0L else nrow(cs),
                  warnings = paste(unique(w), collapse = " | "), cs = cs, max_abs_z = max(abs(z)), top_z = ids[which.max(abs(z))])
      saveRDS(res, file.path(OUT, "susie", sprintf("%s__%s.rds", sn, t)))
      lg("[susie %s] %-8s n %.0f  s %.5f  converged %s niter %s  n_cs %d  warn %d  |z|max %.2f (%s)", sn, t, n, s_hat, res$converged,
         res$niter, res$n_cs, length(w), res$max_abs_z, res$top_z)
    }
  } else lg("[susie %s] all traits checkpointed (skip)", sn)
  for (t in traits) {
    res <- readRDS(file.path(OUT, "susie", sprintf("%s__%s.rds", sn, t)))
    trows[[paste(sn, t)]] <- data.table(set = sn, panel = res$panel, trait = t, nsnp = res$nsnp, n = res$n, s = res$s,
                                        converged = res$converged, niter = res$niter, n_cs = res$n_cs, n_warn = nchar(res$warnings) > 0,
                                        warnings = res$warnings, evaluable = res$converged && !nzchar(res$warnings) && res$n_cs >= 1,
                                        max_abs_z = res$max_abs_z, top_z = res$top_z)
    if (!is.null(res$cs)) csrows[[paste(sn, t)]] <- cbind(data.table(set = sn, panel = res$panel, trait = t), res$cs)
  }
}
trows <- rbindlist(trows); csrows <- rbindlist(csrows, fill = TRUE)
fwrite(trows, file.path(OUT, "susie_traits.tsv"), sep = "\t")
fwrite(csrows, file.path(OUT, "susie_cs.tsv"), sep = "\t")
lg("susie: %d trait fits, %d credible sets written", nrow(trows), nrow(csrows))

# ---------------------------------------------------------------- Rule 3(iv): conditional colocalization
cw <- fread(file = file.path(OUT, "cond_window.tsv"), sep = "\t")
if (nrow(cw) < 1000) stop("cond_window.tsv read only ", nrow(cw), " rows")
cf <- meta$case_fraction; nt <- meta$trait_ntot
outn <- c(finngen_chd = "CHD_FG", ieu_a_7 = "CAD_C4D")
crow <- list()
for (pn in c("503", "489")) for (sg in c("A", "B")) for (oc in names(outn)) {
  dP <- cw[panel == pn & signal == sg & trait == "protein"]; dY <- cw[panel == pn & signal == sg & trait == oc]
  stopifnot(nrow(dP) == nrow(dY), all(dP$snp == dY$snp), nrow(dP) > 100)
  d1 <- list(beta = dP$b, varbeta = dP$se^2, snp = dP$snp, type = "quant", sdY = 1, N = tn[["P"]])
  d2 <- list(beta = dY$b, varbeta = dY$se^2, snp = dY$snp, type = "cc", s = cf[[outn[[oc]]]], N = nt[[outn[[oc]]]])
  r1 <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  d2q <- list(beta = dY$b, varbeta = dY$se^2, snp = dY$snp, type = "quant", sdY = 1)
  r2 <- suppressMessages(coloc.abf(d1, d2q, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  top <- r1$results$snp[which.max(r1$results$SNP.PP.H4)]
  crow[[paste(pn, sg, oc)]] <- data.table(panel = pn, signal = sg, conditioned_on = dP$conditioned_on[1], resource = oc, nsnp = nrow(dP),
    H0 = r1$summary[["PP.H0.abf"]], H1 = r1$summary[["PP.H1.abf"]], H2 = r1$summary[["PP.H2.abf"]], H3 = r1$summary[["PP.H3.abf"]],
    H4 = r1$summary[["PP.H4.abf"]], top_snp_H4 = top, H4_engine_quant_sdY1 = r2$summary[["PP.H4.abf"]], H3_engine_quant_sdY1 = r2$summary[["PP.H3.abf"]],
    max_abs_z_protein = max(abs(dP$z)), max_abs_z_outcome = max(abs(dY$z)))
  lg("[cond coloc] panel %s signal %s (| %s) %-11s nsnp %d  H4 %.4f  H3 %.4f  (engine quant/sdY=1: H4 %.4f)  top %s", pn, sg,
     dP$conditioned_on[1], oc, nrow(dP), r1$summary[["PP.H4.abf"]], r1$summary[["PP.H3.abf"]], r2$summary[["PP.H4.abf"]], top)
}
fwrite(rbindlist(crow), file.path(OUT, "coloc_cond.tsv"), sep = "\t")
lg("w2_04 done")
