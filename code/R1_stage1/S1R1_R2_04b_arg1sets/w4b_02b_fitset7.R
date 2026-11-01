# w4b_02b_fitset7.R -- Stage 2 step 2 (user decision 2026-10-09): s_max over ALL SEVEN distinct SuSiE fits at every locus.
# Existing four fits (w4b_02_susie.R): susie, p1fr, p1np on 503 and p2np on 489 (= p1np on 489). New (s only, no SuSiE fit):
#   p2np on 503 (sub-matrix of R503), susie on 489, p1fr on 489 (both restricted to SNPs valid on the 489 panel).
# 489 validity is the existing p2np rule, applied to EVERY harmonised SNP (palindromes included; w4b_01_data.R tests only
# non-palindromic ones): same GRCh37 position and allele pair (unique), no missing genotype, panel-2 MAF >= 0.01,
# |A1 freq panel 2 - ref| <= 0.05. Asserts: valid & non-palindromic reproduces SETS$p2np exactly; the p2np sub-matrix of the new
# 489 matrix equals D$R489 (max |diff| < 1e-7); this code path reproduces the stored p1np/503 and p2np/489 s (< 1e-6).
# Writes susie/susie_s_extra.tsv (12 rows = 3 fits x 4 traits) and susie/susie_s_extra_idcheck.tsv. Resumable (skips if present).
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
logf("==== W4b step 2b (seven-fit set: p2np/503, susie/489, p1fr/489) start; susieR %s", as.character(packageVersion("susieR")))
fx <- file.path(OUT, "susie", "susie_s_extra.tsv")
if (file.exists(fx)) { logf("susie_s_extra.tsv exists (skip)"); quit(save = "no") }
D <- readRDS(file.path(OUT, "data", "w4b_data.rds")); H <- D$H
W37 <- D$W37
b489 <- sub("\\.bed$", "", rpath("LD489_chr6_bed")); invisible(rpath("LD489_chr6_bim"))
B2 <- fread(paste0(b489, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "A1", "A2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
B2 <- B2[pos >= W37[1] & pos <= W37[2]]
m2 <- merge(H[, .(row, snp, pos37, a1, a2, af_ref, rsid_ref, palindromic)], B2, by.x = "pos37", by.y = "pos", allow.cartesian = TRUE)
m2 <- m2[(A1 == a1 & A2 == a2) | (A1 == a2 & A2 == a1)]
amb <- m2[, .N, by = row][N > 1]$row
m2 <- m2[!(row %in% amb)]
stopifnot(!anyDuplicated(m2$id))
logf("all harmonised SNPs x 489 bim window: %d of %d present (unique pos+allele pair; %d ambiguous dropped)", nrow(m2), nrow(H), length(amb))
dir.create(file.path(OUT, "plink"), showWarnings = FALSE)
writeLines(m2$id, file.path(OUT, "plink", "p489_cand_all.txt"))
invisible(run_plink(c("--bfile", shQuote(b489), "--chr", "6", "--from-bp", W37[1], "--to-bp", W37[2], "--extract", shQuote(file.path(OUT, "plink", "p489_cand_all.txt")),
            "--keep-allele-order", "--freq", "--out", shQuote(file.path(OUT, "plink", "p489_freq_all"))), "p489_freq_all"))
fq <- fread(file.path(OUT, "plink", "p489_freq_all.frq"))
m2 <- merge(m2, fq[, .(SNP, MAF, NCHROBS, A1f = A1)], by.x = "id", by.y = "SNP")
stopifnot(all(m2$A1f == m2$A1))
m2[, af489 := ifelse(A1 == a1, MAF, 1 - MAF)]
n489 <- 978L   # 2 x 489 samples; asserted against the existing p2np step below
m2[, miss := NCHROBS < n489]
m2 <- m2[!miss & pmin(af489, 1 - af489) >= 0.01 & abs(af489 - af_ref) <= 0.05]
setorder(m2, row)
vp <- m2[palindromic == 0]
stopifnot(setequal(vp$snp, D$SETS$p2np), length(vp$snp) == length(D$SETS$p2np))
logf("valid-489 SNPs: %d (non-palindromic %d == existing p2np %d, set-equal; palindromic %d)", nrow(m2), nrow(vp), length(D$SETS$p2np), sum(m2$palindromic == 1))
L <- plink_r_square(b489, m2$id, 6, W37[1], W37[2], "r489_valid")
sg <- ifelse(L$bim$A1 == m2$a1, 1, -1)
Rv <- L$R * outer(sg, sg); dimnames(Rv) <- list(m2$snp, m2$snp)
dd <- max(abs(Rv[D$SETS$p2np, D$SETS$p2np] - D$R489[D$SETS$p2np, D$SETS$p2np]))
logf("CHECK p2np sub-matrix of the new 489 matrix vs existing R489: max |diff| %.3g", dd)
stopifnot(dd < 1e-7)

TRAITS <- list(P__ARG1 = list(z = "z__P__ARG1", n = D$N$P__ARG1), arginine = list(z = "z__M__GCST90200372", n = D$N$M__GCST90200372),
               ornithine = list(z = "z__M__GCST90200404", n = D$N$M__GCST90200404), T2D_FG = list(z = "z__FG", n = D$N$FG))
FITS <- list(
  p2np_503 = list(set = "p2np", panel = "503", snps = D$SETS$p2np),
  susie_489 = list(set = "susie", panel = "489", snps = m2$snp),
  p1fr_489 = list(set = "p1fr", panel = "489", snps = intersect(D$SETS$p1fr, m2$snp)),
  chk_p1np_503 = list(set = "p1np", panel = "503", snps = D$SETS$p1np),
  chk_p2np_489 = list(set = "p2np", panel = "489", snps = D$SETS$p2np))
old <- fread(file.path(OUT, "susie", "susie_s_table.tsv"))
rows <- list()
for (fn in names(FITS)) {
  F <- FITS[[fn]]
  for (tr in names(TRAITS)) {
    z <- H[match(F$snps, snp)][[TRAITS[[tr]]$z]]
    ok <- !is.na(z); snps <- F$snps[ok]; z <- z[ok]; names(z) <- snps
    R <- if (F$panel == "503") D$R503[snps, snps, drop = FALSE] else Rv[snps, snps, drop = FALSE]
    t0 <- Sys.time()
    s_hat <- estimate_s_rss(z, R, TRAITS[[tr]]$n)
    logf("[fitset7 %s] %-9s panel %s: %d SNPs n=%d s=%.8f %.0fs", fn, tr, F$panel, length(snps), TRAITS[[tr]]$n, s_hat, as.numeric(difftime(Sys.time(), t0, units = "secs")))
    rows[[paste(fn, tr)]] <- data.table(trait = tr, set = F$set, panel = F$panel, nsnp = length(snps), n = TRAITS[[tr]]$n, s = s_hat, fit = fn)
  }
}
X <- rbindlist(rows)
for (fn in c("chk_p1np_503", "chk_p2np_489")) for (tr in names(TRAITS)) {
  o <- old[trait == tr & set == FITS[[fn]]$set & panel == FITS[[fn]]$panel]; stopifnot(nrow(o) == 1)
  d <- abs(X[fit == fn & trait == tr]$s - o$s)
  logf("CHECK %s vs stored s (%s): |diff| %.3g", fn, tr, d); stopifnot(d < 1e-6)
}
fwrite(X[grepl("^chk_", fit)], file.path(OUT, "susie", "susie_s_extra_idcheck.tsv"), sep = "\t")
fwrite(X[!grepl("^chk_", fit)], fx, sep = "\t")
logf("==== W4b step 2b done: %d new fit x trait values", nrow(X[!grepl("^chk_", fit)]))
