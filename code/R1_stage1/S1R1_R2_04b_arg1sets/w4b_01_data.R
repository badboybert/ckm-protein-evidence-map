# w4b_01_data.R -- Stage 1 W4b, step 1: inputs, joins (with yields and per-key null rates), SNP sets, PLINK LD (503, 489).
# Output: data/w4b_data.rds (SNP table with every trait aligned to the 1000G A1 = harmonised a1; SNP sets; LD matrices).
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
suppressMessages(library(Rsamtools))
logf("==== W4b step 1 (data) start")
dir.create(file.path(OUT, "data"), showWarnings = FALSE)
JOIN <- list()   # join audit rows
jaudit <- function(name, n_left, n_right, n_matched, key, notes = "") {
  JOIN[[length(JOIN) + 1]] <<- data.table(join = name, n_left = n_left, n_right = n_right, n_matched = n_matched,
                                          yield_left = n_matched / n_left, join_key = key, notes = notes)
  logf("JOIN %-34s left %6d right %7d matched %6d (yield of left %.4f) key=%s %s", name, n_left, n_right, n_matched,
       n_matched / n_left, key, notes)
  if (n_matched == 0) stop("join yielded 0 rows: ", name)
}
nullrate <- function(dt, cols, tag) for (c in cols) logf("null-rate %-24s %-22s %.5f (n=%d)", tag, c,
                                                       mean(is.na(dt[[c]]) | dt[[c]] %in% c("", "NA", ".", "-", "NULL", "TBD", "#N/A")), nrow(dt))

# ---------------------------------------------------------------- harmonised table (R1-pinned, raw_input)
hp <- rpath("R1_ARG1_HARMONISED")
H <- fread(hp, colClasses = list(character = c("snp", "rsid_ref", "rsid_trait", "chr", "a1", "a2")))
stopifnot(nrow(H) == 3067L)                                   # PREP_STAMP n_snps (definition), RESOURCES notes
logf("harmonised: %d rows, %d cols; palindromic %d; window GRCh37 %d-%d, GRCh38 %d-%d", nrow(H), ncol(H), sum(H$palindromic),
     min(H$pos37), max(H$pos37), min(H$pos38), max(H$pos38))
nullrate(H, c("snp", "rsid_ref", "pos37", "pos38", "a1", "a2", "af_ref", "b__P__ARG1", "b__M__GCST90200372", "b__M__GCST90200404"), "harmonised")
stopifnot(!anyDuplicated(H$snp), all(H$chr == "6"))
logf("harmonised: %d positions carry more than one SNP row", sum(duplicated(H$pos37)))
stamp <- fromJSON(rpath("LOC_ARG1_6q23_PREP_STAMP"))
W37 <- stamp$window_grch37_ref; W38 <- stamp$window_grch38
logf("windows from PREP_STAMP (definition): GRCh37 %d-%d, GRCh38 %d-%d", W37[1], W37[2], W38[1], W38[2])
H[, row := .I]

# ---------------------------------------------------------------- raw exposure cross-checks (my own join vs harmonised)
# UKB-PPP ARG1: ID = GRCh37 chr:pos:ref:alt:imp:v1, pos = GRCh38, effect = ALLELE1 (A1FREQ)
P <- fread(rpath("EXP_UKBPPP_EUR_ARG1_OID21445"))
stopifnot(nrow(P) == 11311L)
nullrate(P, c("ID", "pos", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE"), "UKBPPP_ARG1")
P[, c("c37", "p37", "r37", "a37") := tstrsplit(ID, ":", keep = 1:4)]
P[, p37 := as.integer(p37)]
stopifnot(all(P$r37 == P$ALLELE0), all(P$a37 == P$ALLELE1))
m <- merge(H[, .(row, pos37, pos38, a1, a2, b = b__P__ARG1, se = se__P__ARG1)], P[, .(p37, pos, ALLELE0, ALLELE1, A1FREQ, BETA, SE)],
           by.x = "pos37", by.y = "p37", allow.cartesian = TRUE)
m[, s := align_sign(ALLELE1, ALLELE0, a1, a2)]
m <- m[s != 0]
stopifnot(!anyDuplicated(m$row))
jaudit("harmonised x UKBPPP ARG1 (pos37+alleles)", nrow(H), nrow(P), nrow(m), "GRCh37 pos from ID + allele pair")
stopifnot(nrow(m) == nrow(H), all(m$pos == m$pos38))
d <- max(abs(m$s * m$BETA - m$b)); ds <- max(abs(m$SE - m$se))
logf("UKBPPP ARG1 vs harmonised: max |aligned beta - b| %.3g, max |se diff| %.3g; GRCh38 pos agrees on all rows", d, ds)
stopifnot(d < 1e-6, ds < 1e-6)
H[m$row, af_P := ifelse(m$s == 1, m$A1FREQ, 1 - m$A1FREQ)]

readclsa <- function(id, col) {
  X <- fread(rpath(id))
  nullrate(X, c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency"), id)
  mm <- merge(H[, .(row, pos38, a1, a2, b = get(paste0("b__", col)), se = get(paste0("se__", col)))],
              X[, .(base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency)],
              by.x = "pos38", by.y = "base_pair_location", allow.cartesian = TRUE)
  mm[, s := align_sign(effect_allele, other_allele, a1, a2)]
  mm <- mm[s != 0]
  stopifnot(!anyDuplicated(mm$row))
  jaudit(paste0("harmonised x ", id), nrow(H), nrow(X), nrow(mm), "GRCh38 pos + allele pair")
  stopifnot(nrow(mm) == nrow(H))
  d <- max(abs(mm$s * mm$beta - mm$b)); ds <- max(abs(mm$standard_error - mm$se))
  logf("%s vs harmonised: max |aligned beta - b| %.3g, max |se diff| %.3g", id, d, ds)
  stopifnot(d < 1e-6, ds < 1e-6)
  af <- rep(NA_real_, nrow(H)); af[mm$row] <- ifelse(mm$s == 1, mm$effect_allele_frequency, 1 - mm$effect_allele_frequency)
  af
}
H[, af_M372 := readclsa("EXP_CLSA_arginine_GCST90200372", "M__GCST90200372")]
H[, af_M404 := readclsa("EXP_CLSA_ornithine_GCST90200404", "M__GCST90200404")]

# ---------------------------------------------------------------- FinnGen R10 T2D (GRCh38; tabix the window)
fgp <- rpath("OUT_FG_R10_T2D")
tf <- TabixFile(fgp, index = paste0(fgp, ".tbi"))
logf("FinnGen .tbi: %s md5 %s (RESOURCES notes give 7c7965b8e6794d20e897be42a8666895)", paste0(fgp, ".tbi"), md5(paste0(fgp, ".tbi")))
stopifnot(md5(paste0(fgp, ".tbi")) == "7c7965b8e6794d20e897be42a8666895")
hdr <- headerTabix(tf)$header
lines <- scanTabix(tf, param = GRanges("6", IRanges(W38[1], W38[2])))[[1]]
FG <- fread(text = c(sub("^#", "", hdr), lines), sep = "\t", colClasses = list(character = c("chrom", "ref", "alt", "rsids")))
logf("FinnGen T2D window chr6:%d-%d (GRCh38): %d rows", W38[1], W38[2], nrow(FG))
stopifnot(nrow(FG) > 1000, all(FG$chrom == "6"), all(FG$pos >= W38[1] & FG$pos <= W38[2]))
nullrate(FG, c("pos", "ref", "alt", "rsids", "beta", "sebeta", "af_alt", "pval"), "FinnGen_T2D")
fwrite(FG, file.path(OUT, "data", "finngen_R10_T2D_ARG1window.tsv.gz"), sep = "\t")
m <- merge(H[, .(row, pos38, a1, a2, af_ref)], FG[, .(pos, ref, alt, beta, sebeta, af_alt, pval)], by.x = "pos38", by.y = "pos", allow.cartesian = TRUE)
m[, s := align_sign(alt, ref, a1, a2)]
nmulti <- m[s != 0, .N, by = row][N > 1, .N]
m <- m[s != 0]
m <- m[!(row %in% m[, .N, by = row][N > 1]$row)]          # ambiguous multi-allelic match dropped, never guessed
jaudit("harmonised x FinnGen R10 T2D", nrow(H), nrow(FG), nrow(m), "GRCh38 pos + allele pair", sprintf("ambiguous rows dropped: %d", nmulti))
m[, af_D := ifelse(s == 1, af_alt, 1 - af_alt)]
bad <- abs(m$af_D - m$af_ref) > 0.20
logf("FinnGen T2D: |aligned AF - ref AF| > 0.20 dropped: %d of %d (engine QC)", sum(bad), nrow(m))
m <- m[!bad]
H[m$row, `:=`(b__FG = m$s * m$beta, se__FG = m$sebeta, af_FG = m$af_D, p_FG = m$pval)]
H[, z__FG := b__FG / se__FG]
logf("FinnGen T2D available on %d of %d harmonised SNPs", sum(!is.na(H$b__FG)), nrow(H))

# ---------------------------------------------------------------- DIAMANTE 2018 EUR noUKB (GRCh37; uncompressed 1.1 GB: stream-cut)
dmp <- rpath("OUT_DIAMANTE2018_noUKB_EUR_chrpos")
cut <- file.path(OUT, "data", sprintf("diamante2018_noUKB_EUR__chr6_%d_%d.tsv.gz", W37[1], W37[2]))
if (!file.exists(cut)) {
  t0 <- Sys.time(); con <- file(dmp, "r"); hd <- readLines(con, 1); n <- 0L; keep <- character(0)
  repeat { x <- readLines(con, 2e6); if (!length(x)) break; n <- n + length(x)
    k <- x[startsWith(x, "6:")]
    if (length(k)) { pos <- as.integer(sub("^6:([0-9]+).*$", "\\1", k)); keep <- c(keep, k[pos >= W37[1] & pos <= W37[2]]) }
    if (n %% 4e6 == 0) logf("DIAMANTE stream: %d rows read, %d kept", n, length(keep)) }
  close(con)
  logf("DIAMANTE stream-cut: %d data rows read (RESOURCES notes 21,508,698), %d in window, %.0fs", n, length(keep), as.numeric(difftime(Sys.time(), t0, units = "secs")))
  stopifnot(n == 21508698L, hd == "SNP\tChr\tPos\tEA\tNEA\tEAF\tBeta\tSE\tPvalue")
  gz <- gzfile(paste0(cut, ".part"), "w"); writeLines(c(hd, keep), gz); close(gz); file.rename(paste0(cut, ".part"), cut)
}
logf("DIAMANTE window extract %s md5 %s", cut, md5(cut))
DM <- fread(cut, colClasses = list(character = c("SNP", "Chr", "EA", "NEA")))
stopifnot(nrow(DM) > 1000, all(DM$Chr == "6"))
nullrate(DM, c("SNP", "Pos", "EA", "NEA", "EAF", "Beta", "SE", "Pvalue"), "DIAMANTE")
m <- merge(H[, .(row, pos37, a1, a2, af_ref)], DM[, .(Pos, EA, NEA, EAF, Beta, SE, Pvalue)], by.x = "pos37", by.y = "Pos", allow.cartesian = TRUE)
m[, s := align_sign(EA, NEA, a1, a2)]
nmulti <- m[s != 0, .N, by = row][N > 1, .N]
m <- m[s != 0]; m <- m[!(row %in% m[, .N, by = row][N > 1]$row)]
jaudit("harmonised x DIAMANTE noUKB", nrow(H), nrow(DM), nrow(m), "GRCh37 pos + allele pair", sprintf("ambiguous rows dropped: %d", nmulti))
m[, af_D := ifelse(s == 1, EAF, 1 - EAF)]
bad <- abs(m$af_D - m$af_ref) > 0.20
logf("DIAMANTE: |aligned AF - ref AF| > 0.20 dropped: %d of %d (engine QC)", sum(bad), nrow(m))
m <- m[!bad]
H[m$row, `:=`(b__DM = m$s * m$Beta, se__DM = m$SE, af_DM = m$af_D, p_DM = m$Pvalue)]
H[, z__DM := b__DM / se__DM]
logf("DIAMANTE available on %d of %d harmonised SNPs", sum(!is.na(H$b__DM)), nrow(H))

# ---------------------------------------------------------------- sample sizes (computed, never typed)
nhat <- function(se, f) { ok <- !is.na(f) & pmin(f, 1 - f) > 0.05 & se > 0; v <- 1 / (2 * f[ok] * (1 - f[ok]) * se[ok]^2); stopifnot(sum(ok) >= 50); median(v) }
N <- list(P__ARG1 = round(nhat(H$se__P__ARG1, H$af_P)), M__GCST90200372 = round(nhat(H$se__M__GCST90200372, H$af_M372)),
          M__GCST90200404 = round(nhat(H$se__M__GCST90200404, H$af_M404)))
cc <- function(id) { s <- res(id)$n_info; ca <- as.numeric(sub(".*cases ([0-9]+).*", "\\1", s)); co <- as.numeric(sub(".*controls ([0-9]+).*", "\\1", s)); c(ca = ca, co = co) }
man <- fread(rpath("OUT_FG_R10_MANIFEST"))
fgc <- man[phenocode == "T2D"]; stopifnot(nrow(fgc) == 1)
fgr <- cc("OUT_FG_R10_T2D"); stopifnot(fgr[["ca"]] == fgc$num_cases, fgr[["co"]] == fgc$num_controls)
dmc <- cc("OUT_DIAMANTE2018_noUKB_EUR_chrpos")
neff <- function(x) 4 * x[["ca"]] * x[["co"]] / (x[["ca"]] + x[["co"]])
N$FG <- round(neff(c(ca = fgc$num_cases, co = fgc$num_controls))); N$DM <- round(neff(dmc))
CASEFRAC <- list(FG = fgc$num_cases / (fgc$num_cases + fgc$num_controls), DM = dmc[["ca"]] / (dmc[["ca"]] + dmc[["co"]]))
logf("N used: ARG1 n_hat %d (RESOURCES n_info: l3-finemap n_hat 36979); arginine n_hat %d (8361); ornithine n_hat %d (8301)",
     N$P__ARG1, N$M__GCST90200372, N$M__GCST90200404)
logf("N used: FinnGen T2D N_eff %d from manifest counts %d/%d (RESOURCES n_info N_eff 218000); DIAMANTE N_eff %d from RESOURCES counts %d/%d (n_info 193440)",
     N$FG, fgc$num_cases, fgc$num_controls, N$DM, dmc[["ca"]], dmc[["co"]])

# ---------------------------------------------------------------- 503 panel: bim window and PLINK LD on every harmonised SNP
b503 <- sub("\\.bed$", "", rpath("LD503_bed")); invisible(rpath("LD503_fam"))
logf("LD503_bim (658 MB, md5 asserted once): %s", rpath("LD503_bim"))
invisible(run_plink(c("--bfile", shQuote(b503), "--chr", "6", "--from-bp", W37[1], "--to-bp", W37[2], "--keep-allele-order", "--make-just-bim",
            "--out", shQuote(file.path(OUT, "plink", "win503"))), "win503_bim"))
bw <- fread(file.path(OUT, "plink", "win503.bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "A1", "A2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
logf("503 bim window: %d variants; duplicated ids %d", nrow(bw), sum(duplicated(bw$id)))
mb <- merge(H[, .(row, rsid_ref, pos37, a1, a2)], bw, by.x = "pos37", by.y = "pos", allow.cartesian = TRUE)
mb <- mb[(A1 == a1 & A2 == a2) | (A1 == a2 & A2 == a1)]
stopifnot(!anyDuplicated(mb$row), nrow(mb) == nrow(H), all(mb$id == mb$rsid_ref), !anyDuplicated(bw$id[bw$id %in% H$rsid_ref]))
jaudit("harmonised x 503 bim window", nrow(H), nrow(bw), nrow(mb), "GRCh37 pos + allele pair; id == rsid_ref")
setorder(mb, row)
H[, A1_503 := mb$A1]
logf("503 bim A1 == harmonised a1 on %d of %d SNPs", sum(H$A1_503 == H$a1), nrow(H))
L5 <- plink_r_square(b503, H$rsid_ref, 6, W37[1], W37[2], "r503_all")
sg5 <- ifelse(H$A1_503 == H$a1, 1, -1)
R503 <- L5$R * outer(sg5, sg5); dimnames(R503) <- list(H$snp, H$snp)
logf("R503: %d x %d from PLINK --r square; sign aligned to harmonised a1 (%d SNPs flipped)", nrow(R503), ncol(R503), sum(sg5 < 0))

# cross-check against the R1-pinned 503 dosage extract (read in place; correlation only, no copy written)
gp <- rpath("R1_ARG1_REF503_I8")
G <- matrix(readBin(gp, "integer", n = 503L * 3067L + 1L, size = 1L, signed = TRUE), nrow = 503, byrow = TRUE)
stopifnot(ncol(G) == 3067, all(G %in% 0:2), max(abs(colMeans(G) / 2 - H$af_ref)) < 1e-5)
Rg <- cor(G)
logf("CHECK PLINK R503 vs cor(503 i8 dosage): max |diff| %.3g (expected ~1e-6 or below: same 503 hard calls)", max(abs(Rg - R503)))
stopifnot(max(abs(Rg - R503)) < 1e-4)
rm(Rg)

# LD assertion mutation test: one requested id mutated -> the 'every requested SNP returned' assertion must trip
mut_ids <- H$rsid_ref[1:20]; mut_ids[7] <- paste0(mut_ids[7], "_MUTATED")
tr <- tryCatch({ plink_r_square(b503, mut_ids, 6, W37[1], W37[2], "r503_mutation_test"); "NOT TRIPPED" }, error = function(e) conditionMessage(e))
logf("LD mutation test (id %s): %s", mut_ids[7], tr)
stopifnot(grepl("LD ASSERT FAILED", tr))

# ---------------------------------------------------------------- 489 panel (1000G.EUR.QC.6; nested in 503)
b489 <- sub("\\.bed$", "", rpath("LD489_chr6_bed")); invisible(rpath("LD489_chr6_bim"))
f489 <- xpath("LD489_chr6_fam"); f1 <- rpath("LD489_fam")
stopifnot(md5(f489) == md5(f1))
logf("489 chr6 .fam md5 equals LD489_fam (chr1) md5: identical sample list")
iid503 <- fread(paste0(b503, ".fam"), header = FALSE)$V2; iid489 <- fread(f489, header = FALSE)$V2
logf("489 IIDs in 503: %d of %d", sum(iid489 %in% iid503), length(iid489))
B2 <- fread(paste0(b489, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "A1", "A2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
B2 <- B2[pos >= W37[1] & pos <= W37[2]]
np <- H[palindromic == 0]
m2 <- merge(np[, .(row, pos37, a1, a2, af_ref, rsid_ref)], B2, by.x = "pos37", by.y = "pos", allow.cartesian = TRUE)
m2 <- m2[(A1 == a1 & A2 == a2) | (A1 == a2 & A2 == a1)]
amb <- m2[, .N, by = row][N > 1]$row
m2 <- m2[!(row %in% amb)]
jaudit("p1np x 489 bim window", nrow(np), nrow(B2), nrow(m2), "GRCh37 pos + allele pair (unique)", sprintf("ambiguous %d", length(amb)))
stopifnot(!anyDuplicated(m2$id))
n_cand <- nrow(m2)
# 489 A1 frequency via PLINK --freq (keep-allele-order: frequency of bim A1)
writeLines(m2$id, file.path(OUT, "plink", "p489_cand.txt"))
invisible(run_plink(c("--bfile", shQuote(b489), "--chr", "6", "--from-bp", W37[1], "--to-bp", W37[2], "--extract", shQuote(file.path(OUT, "plink", "p489_cand.txt")),
            "--keep-allele-order", "--freq", "--out", shQuote(file.path(OUT, "plink", "p489_freq"))), "p489_freq"))
fq <- fread(file.path(OUT, "plink", "p489_freq.frq"))
m2 <- merge(m2, fq[, .(SNP, MAF, NCHROBS, A1f = A1)], by.x = "id", by.y = "SNP")
stopifnot(all(m2$A1f == m2$A1))
m2[, af489 := ifelse(A1 == a1, MAF, 1 - MAF)]        # PLINK .frq MAF column = frequency of A1 (with --keep-allele-order)
m2[, miss := NCHROBS < 2 * length(iid489)]
drops <- c(missing = sum(m2$miss), maf = sum(!m2$miss & pmin(m2$af489, 1 - m2$af489) < 0.01),
           afdiff = sum(!m2$miss & pmin(m2$af489, 1 - m2$af489) >= 0.01 & abs(m2$af489 - m2$af_ref) > 0.05))
m2 <- m2[!miss & pmin(af489, 1 - af489) >= 0.01 & abs(af489 - af_ref) <= 0.05]
setorder(m2, row)
logf("p2np: %d SNPs (p1np %d; absent/ambiguous in 489 bim %d; drops: missing %d, MAF<0.01 %d, |dAF|>0.05 %d); max |dAF| %.4f",
     nrow(m2), nrow(np), nrow(np) - n_cand, drops[["missing"]], drops[["maf"]], drops[["afdiff"]], max(abs(m2$af489 - m2$af_ref)))
L4 <- plink_r_square(b489, m2$id, 6, W37[1], W37[2], "r489_p2np")
sg4 <- ifelse(L4$bim$A1 == m2$a1, 1, -1)
R489 <- L4$R * outer(sg4, sg4); snp489 <- H$snp[m2$row]; dimnames(R489) <- list(snp489, snp489)
logf("R489: %d x %d from PLINK --r square on 1000G.EUR.QC.6; sign aligned to harmonised a1 (%d flipped)", nrow(R489), ncol(R489), sum(sg4 < 0))
dd <- max(abs(R489 - R503[snp489, snp489]))
logf("CHECK 489 vs 503 r on p2np: max |diff| %.4f (nested panels; expected small but non-zero)", dd)

# ---------------------------------------------------------------- SNP sets (engine QC; DEVIATION: engine path-locked, re-implemented)
# susie = every harmonised SNP (ref MAF >= 0.01, |dAF| <= 0.20, palindromes only at ref MAF <= 0.40: applied by the R1-pinned prep)
# p1np  = harmonised minus palindromes
# p1fr  = p1np + palindromes with ref MAF <= 0.40 and, for every trait fitted here that reports an AF (UKB-PPP A1FREQ, CLSA EAF,
#         FinnGen af_alt), GWAS MAF <= 0.40 and |GWAS A1 freq - ref A1 freq| <= 0.10 (l3fm__robustness 'freq' rule)
# p2np  = p1np present in the 489 panel (unique pos+alleles, no missing, MAF >= 0.01, |dAF| <= 0.05)
mafr <- pmin(H$af_ref, 1 - H$af_ref)
okf <- function(f) is.na(f) | (pmin(f, 1 - f) <= 0.40 & abs(f - H$af_ref) <= 0.10)
palok <- H$palindromic == 1 & mafr <= 0.40 & okf(H$af_P) & okf(H$af_M372) & okf(H$af_M404) & okf(H$af_FG)
SETS <- list(susie = H$snp, p1fr = H$snp[H$palindromic == 0 | palok], p1np = H$snp[H$palindromic == 0], p2np = snp489)
PANEL <- c(susie = "503", p1fr = "503", p1np = "503", p2np = "489")
for (s in names(SETS)) logf("SNP set %-5s panel %s: %d SNPs (T2D FinnGen available on %d; DIAMANTE on %d)", s, PANEL[[s]], length(SETS[[s]]),
                            sum(!is.na(H[match(SETS[[s]], snp)]$b__FG)), sum(!is.na(H[match(SETS[[s]], snp)]$b__DM)))
logf("p1fr: %d palindromes in harmonised, %d kept by the frequency rule", sum(H$palindromic), sum(palok))

JOIN <- rbindlist(JOIN); fwrite(JOIN, file.path(OUT, "join_audit.tsv"), sep = "\t")
saveRDS(list(H = H, R503 = R503, R489 = R489, SETS = SETS, PANEL = PANEL, N = N, CASEFRAC = CASEFRAC, W37 = W37, W38 = W38,
             A1_489 = data.table(snp = snp489, A1 = L4$bim$A1, af489 = m2$af489)),
        file.path(OUT, "data", "w4b_data.rds"))
logf("==== W4b step 1 done: data/w4b_data.rds")
