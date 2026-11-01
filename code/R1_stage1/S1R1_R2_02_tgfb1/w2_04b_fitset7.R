# w2_04b_fitset7.R -- Stage 2 step 2 (user decision 2026-10-09): the three additional Rule 3(i) fits, so that s_max is over
# ALL SEVEN distinct SuSiE fits (all variant sets, both panels). Only s (susieR::estimate_s_rss, null-mle) is needed, so no
# susie_rss fit / credible set is made for the new fits. Existing four fits are untouched (w2_04_susie_coloc.R).
#   p2np_503panel : p2np SNP set on the 503 matrix (sub-matrix of R503_susie.f8)
#   susie_489panel: every valid-489 SNP of the susie set (fitset7_sets.tsv; ld/R489_valid.f8)
#   p1fr_489panel : p1fr SNPs that are valid-489
# (p1np on 489 is p2np on 489, already fitted.) Liveness/identity check: the same code path re-estimates the existing p1np@503
# and p2np@489 s and must reproduce susie_traits.tsv within 1e-6.
#   Rscript w2_04b_fitset7.R <OUT dir>
suppressMessages({library(data.table); library(susieR); library(jsonlite)})
a <- commandArgs(trailingOnly = TRUE)
OUT <- a[1]
logf <- file.path(OUT, "run.log")
lg <- function(...) { m <- sprintf("%s  %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
lg("==== w2_04b_fitset7.R (susieR %s) ====", packageVersion("susieR"))
meta <- fromJSON(file.path(OUT, "tgfb1_locus_meta.json"))
chk <- function(f, m) { got <- unname(tools::md5sum(f)); if (got != m) stop("md5 mismatch ", f); lg("md5 ok %s %s", basename(f), got) }
chk(file.path(OUT, "tgfb1_locus.tsv"), meta$md5$tgfb1_locus.tsv)
chk(file.path(OUT, "ld", "R503_susie.f8"), meta$md5$R503_susie.f8)
chk(file.path(OUT, "ld", "R489_p2np.f8"), meta$md5$R489_p2np.f8)
h <- fread(file = file.path(OUT, "tgfb1_locus.tsv"), sep = "\t")
fs <- fread(file = file.path(OUT, "fitset7_sets.tsv"), sep = "\t")
stopifnot(nrow(fs) == nrow(h), all(fs$i == h$i), all(fs$rsid503 == h$rsid503))
p <- nrow(h)
R503 <- matrix(readBin(file.path(OUT, "ld", "R503_susie.f8"), "double", n = p * p, size = 8, endian = "little"), p, p)
iv <- which(fs$valid489 == 1)
Rv <- matrix(readBin(file.path(OUT, "ld", "R489_valid.f8"), "double", n = length(iv)^2, size = 8, endian = "little"), length(iv), length(iv))
stopifnot(isSymmetric(R503, tol = 1e-6), isSymmetric(Rv, tol = 1e-6), all(abs(diag(Rv) - 1) < 1e-6))
tn <- unlist(meta$trait_n); traits <- c("P", "CHD_FG", "CAD_C4D")
i2 <- which(h$in_p2np == 1)
stopifnot(all(i2 %in% iv))
# fit definitions: idx = rows of tgfb1_locus.tsv; R as a matrix on exactly those rows
fits <- list(
  p2np_503panel  = list(panel = "503", idx = i2,                                R = function(ix) R503[ix, ix]),
  susie_489panel = list(panel = "489", idx = iv,                                R = function(ix) Rv),
  p1fr_489panel  = list(panel = "489", idx = intersect(which(h$in_p1fr == 1), iv), R = function(ix) Rv[match(ix, iv), match(ix, iv)]),
  # identity checks of this code path against the stored four-fit values
  chk_p1np_503   = list(panel = "503", idx = which(h$in_p1np == 1),             R = function(ix) R503[ix, ix]),
  chk_p2np_489   = list(panel = "489", idx = i2,                                R = function(ix) Rv[match(ix, iv), match(ix, iv)]))
old <- fread(file.path(OUT, "susie_traits.tsv"))
rows <- list()
for (fn in names(fits)) {
  F <- fits[[fn]]; ix <- F$idx
  R <- F$R(ix); ids <- h$rsid503[ix]; dimnames(R) <- list(ids, ids)
  stopifnot(nrow(R) == length(ix), all(abs(diag(R) - 1) < 1e-6))
  te <- system.time(attr(R, "eigen") <- eigen(R, symmetric = TRUE))[["elapsed"]]
  for (t in traits) {
    z <- h[[paste0("z__", t)]][ix]; names(z) <- ids
    s_hat <- estimate_s_rss(z, R, tn[[t]])
    lg("[fitset7 %s] %-8s nsnp %d panel %s n %.0f s %.8f (eigen %.0fs)", fn, t, length(ix), F$panel, tn[[t]], s_hat, te)
    rows[[paste(fn, t)]] <- data.table(set = fn, panel = F$panel, trait = t, nsnp = length(ix), n = tn[[t]], s = s_hat)
  }
}
X <- rbindlist(rows)
# identity checks
for (cc in list(c("chk_p1np_503", "p1np"), c("chk_p2np_489", "p2np"))) {
  for (t in traits) {
    d <- abs(X[set == cc[1] & trait == t]$s - old[set == cc[2] & trait == t]$s)
    lg("CHECK %s vs stored %s %s: |diff| %.3g", cc[1], cc[2], t, d)
    stopifnot(d < 1e-6)
  }
}
fwrite(X[!grepl("^chk_", set)], file.path(OUT, "susie_s_extra.tsv"), sep = "\t")
fwrite(X[grepl("^chk_", set)], file.path(OUT, "susie_s_extra_idcheck.tsv"), sep = "\t")
lg("w2_04b done: %d new fit x trait values", nrow(X[!grepl("^chk_", set)]))
