# w2_06_coloc_boot.R -- Rule 3(ii) bootstrap of the conditional colocalization at the TGFB1 signals (section 2).
#   Rscript w2_06_coloc_boot.R <OUT dir>
# Rule 3(ii): conditional-colocalization results are robust when the conclusion (side of PP.H4 0.80) is the same on the
# 503 and 489 panels and in >= 90% of 5,000 503-panel bootstrap draws (seed 20261007). The draws are the SAME person
# resamples as the w2_03_cond.py ratio bootstrap (boot_idx_5000x503.i2, numpy default_rng(20261007)).
# Genotypes: the 503 .bed read IN PLACE with readBin + seek (no copy); SNPs = the p1np set (the "503" label).
# Per draw: r of every p1np SNP with the conditioning SNP from the resampled people; z-space conditioning (r^2 > 0.9
# dropped, as the engine); coloc.abf protein (quant, sdY 1) x CAD (cc, case fraction, total N), p12 = 1e-5.
# A SNP monomorphic in a draw gets r = 0 (counted). Checkpoint every 250 draws (resumable).
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
suppressMessages({library(data.table); library(coloc); library(jsonlite)})
a <- commandArgs(trailingOnly = TRUE); OUT <- a[1]
logf <- file.path(OUT, "run.log")
lg <- function(...) { m <- sprintf("%s  %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
lg("==== w2_06_coloc_boot.R ====")
meta <- fromJSON(file.path(OUT, "tgfb1_locus_meta.json"))
stopifnot(unname(tools::md5sum(file.path(OUT, "tgfb1_locus.tsv"))) == meta$md5$tgfb1_locus.tsv)
h <- fread(file = file.path(OUT, "tgfb1_locus.tsv"), sep = "\t")
stopifnot(nrow(h) == meta$n$susie)
resj <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/RESOURCES.json")
stopifnot(unname(tools::md5sum(resj)) == "546b769ae93a8f566eec2b92f4be1575")
res <- fromJSON(resj)
bedp <- res$path[res$id == "LD503_bed"]; bedmd5 <- res$md5[res$id == "LD503_bed"]
lg("503 bed %s (RESOURCES md5 %s; size asserted, md5 asserted by the Python steps of this run)", bedp, bedmd5)
stopifnot(file.size(bedp) == res$bytes[res$id == "LD503_bed"])
idx <- which(h$in_p1np == 1); m <- length(idx)
nb <- (503 + 3) %/% 4
con <- file(bedp, "rb"); stopifnot(identical(readBin(con, "raw", 3), as.raw(c(0x6c, 0x1b, 0x01))))
lut <- c(2, NA, 1, 0)
G <- matrix(NA_real_, 503, m)
for (j in seq_len(m)) {
  seek(con, 3 + as.numeric(h$bed_idx503[idx[j]]) * nb)
  bits <- as.integer(rawToBits(readBin(con, "raw", nb)))
  code <- bits[c(TRUE, FALSE)] + 2L * bits[c(FALSE, TRUE)]
  G[, j] <- lut[code[1:503] + 1L]
}
close(con)
stopifnot(!anyNA(G))
af <- colMeans(G) / 2
lg("genotypes read in place for %d p1np SNPs; max |A1 freq - af503| = %.2e", m, max(abs(af - h$af503[idx])))
stopifnot(max(abs(af - h$af503[idx])) < 1e-4)
p <- nrow(h)
R503 <- matrix(readBin(file.path(OUT, "ld", "R503_susie.f8"), "double", n = p * p, size = 8, endian = "little"), p, p)
iA <- match("rs73045269", h$rsid503[idx]); iB <- match("rs12983775", h$rsid503[idx])
dA <- max(abs(cor(G, G[, iA]) - R503[idx, idx[iA]])); dB <- max(abs(cor(G, G[, iB]) - R503[idx, idx[iB]]))
lg("readBin genotypes vs PLINK r (columns A and B over p1np): max |diff| %.2e / %.2e", dA, dB)
stopifnot(dA < 1e-5, dB < 1e-5)
tn <- unlist(meta$trait_n); nt <- unlist(meta$trait_ntot); cf <- unlist(meta$case_fraction)
Z <- list(P = h$z__P[idx], CHD_FG = h$z__CHD_FG[idx], CAD_C4D = h$z__CAD_C4D[idx])
SE <- list(P = h$se__P[idx], CHD_FG = h$se__CHD_FG[idx], CAD_C4D = h$se__CAD_C4D[idx])
snp <- h$snp[idx]
one <- function(rows) {
  Gd <- G[rows, , drop = FALSE]
  out <- c()
  nmono <- 0L
  for (sg in c("A", "B")) {
    c0 <- if (sg == "A") iB else iA
    r <- suppressWarnings(as.vector(cor(Gd, Gd[, c0])))
    mono <- is.na(r); nmono <- nmono + sum(mono); r[mono] <- 0
    q <- r^2; keep <- q <= 0.9; keep[c0] <- FALSE
    cz <- function(t) { zc <- (Z[[t]] - r * Z[[t]][c0]) / sqrt(pmax(1 - q, 1e-12)); sc <- SE[[t]] / sqrt(pmax(1 - q, 1e-12)); list(b = (zc * sc)[keep], v = (sc^2)[keep]) }
    P <- cz("P")
    d1 <- list(beta = P$b, varbeta = P$v, snp = snp[keep], type = "quant", sdY = 1, N = tn[["P"]])
    for (oc in c("CHD_FG", "CAD_C4D")) {
      Y <- cz(oc)
      d2 <- list(beta = Y$b, varbeta = Y$v, snp = snp[keep], type = "cc", s = cf[[oc]], N = nt[[oc]])
      out[paste0("H4_", sg, "_", oc)] <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))$summary[["PP.H4.abf"]]
    }
    out[paste0("nsnp_", sg)] <- sum(keep)
  }
  out["n_monomorphic"] <- nmono
  out
}
full <- one(1:503)
cc <- fread(file = file.path(OUT, "coloc_cond.tsv"), sep = "\t")
ref <- c(H4_A_CHD_FG = cc[panel == 503 & signal == "A" & resource == "finngen_chd"]$H4, H4_A_CAD_C4D = cc[panel == 503 & signal == "A" & resource == "ieu_a_7"]$H4,
         H4_B_CHD_FG = cc[panel == 503 & signal == "B" & resource == "finngen_chd"]$H4, H4_B_CAD_C4D = cc[panel == 503 & signal == "B" & resource == "ieu_a_7"]$H4)
lg("identity draw (all 503 people once) vs w2_04 full-panel PP.H4: max |diff| = %.2e (tolerance 1e-3: PLINK prints r to 6 significant digits, the readBin r is exact)", max(abs(full[names(ref)] - ref)))
stopifnot(max(abs(full[names(ref)] - ref)) < 1e-3)
IDX <- matrix(readBin(file.path(OUT, "boot_idx_5000x503.i2"), "integer", n = 5000 * 503, size = 2, endian = "little"), 5000, 503, byrow = TRUE)
stopifnot(min(IDX) >= 0, max(IDX) <= 502)
ckp <- file.path(OUT, "coloc_boot_draws.tsv")
done <- if (file.exists(ckp)) fread(file = ckp, sep = "\t") else NULL
start <- if (is.null(done)) 1L else nrow(done) + 1L
if (!is.null(done)) stopifnot(all(done$draw == seq_len(nrow(done))))
lg("coloc bootstrap: resuming at draw %d of 5000", start)
t0 <- Sys.time()
buf <- list()
for (d in if (start <= 5000L) seq.int(start, 5000L) else integer(0)) {
  buf[[length(buf) + 1L]] <- as.data.table(as.list(c(draw = d, one(IDX[d, ] + 1L))))
  if (d %% 250L == 0L || d == 5000L) {
    fwrite(rbindlist(buf), ckp, sep = "\t", append = file.exists(ckp))
    buf <- list()
    lg("[coloc boot] %d / 5000 draws (%.0fs this session)", d, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
}
B <- fread(file = ckp, sep = "\t")
stopifnot(nrow(B) == 5000, all(B$draw == 1:5000), !anyNA(B))
sm <- rbindlist(lapply(names(ref), function(k) {
  same <- (B[[k]] >= 0.80) == (ref[[k]] >= 0.80)
  data.table(quantity = k, full_503 = ref[[k]], full_side = ifelse(ref[[k]] >= 0.80, ">=0.80", "<0.80"), share_same_side = mean(same),
             mc_se = sqrt(mean(same) * (1 - mean(same)) / 5000), h4_q025 = quantile(B[[k]], 0.025), h4_q500 = quantile(B[[k]], 0.5),
             h4_q975 = quantile(B[[k]], 0.975), n_draws = nrow(B), monomorphic_snp_draws = sum(B$n_monomorphic > 0))
}))
fwrite(sm, file.path(OUT, "coloc_boot_summary.tsv"), sep = "\t")
for (i in seq_len(nrow(sm))) lg("[coloc boot] %s: full %.4f (%s); share same side %.4f (MC SE %.4f); PP.H4 2.5/50/97.5%% %.3f / %.3f / %.3f",
                                sm$quantity[i], sm$full_503[i], sm$full_side[i], sm$share_same_side[i], sm$mc_se[i], sm$h4_q025[i], sm$h4_q500[i], sm$h4_q975[i])
lg("w2_06 done")
