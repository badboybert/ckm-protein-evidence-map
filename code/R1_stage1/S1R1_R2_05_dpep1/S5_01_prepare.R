# S5_01_prepare.R -- W5 (§5 DPEP1) step 1: assert governing md5s, snapshot staging, re-derive the harmonised trait
# columns from the raw windows, rebuild the four SNP sets, compute LD with PLINK on both panels, cut eGFR from CKDGen.
# Output: work/prep.rds (+ ld/*.bin matrices, which are LD matrices only; no genotype file is copied).
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_05_dpep1/S5_common.R"))
t0 <- Sys.time()
say("==== S5_01_prepare start ====")
assert_governing()

# ---------------------------------------------------------------- staging snapshot (start)
if (!file.exists(file.path(WORK, "staging_md5_start.tsv"))) snap_staging("start") else say("staging start snapshot exists (resume)")

# ---------------------------------------------------------------- tools
say("R: %s; data.table %s; susieR %s; coloc %s; jsonlite %s; rtracklayer %s; digest %s", R.version.string,
    packageVersion("data.table"), packageVersion("susieR"), packageVersion("coloc"), packageVersion("jsonlite"),
    packageVersion("rtracklayer"), packageVersion("digest"))
say("RNGkind: %s", paste(RNGkind(), collapse = " / "))
pl <- res("PLINK19_paper6"); assert_md5(pl$path, pl$md5, "PLINK19_paper6")
v <- run_plink("--version", "plink_version"); say("PLINK --version: %s", paste(v$out, collapse = " "))
st <- system2(file.path(R.home("bin"), "Rscript.exe"), args = c(q(GOV$GLS_Q$p), "--self-test"), stdout = TRUE, stderr = TRUE)
for (l in st) say("gls_q.R --self-test | %s", l)
if (!any(grepl("SELF-TEST PASS", st))) stop("gls_q.R self-test did not pass")

# ---------------------------------------------------------------- harmonised table (R1-pinned derived input)
hp <- inp("REF_DPEP1_HARMONISED")
H <- fread(hp, colClasses = list(character = c("snp", "rsid_ref", "rsid_trait", "chr", "a1", "a2")))
stopifnot(nrow(H) == 3575L)   # fread on a path with a space can return empty: assert
say("harmonised: %d rows x %d cols; palindromic %d; NA cells in used columns: %d", nrow(H), ncol(H), sum(H$palindromic == 1),
    sum(is.na(H[, .(pos37, pos38, af_ref, b__P__DPEP1, se__P__DPEP1, b__P__CHMP1A, se__P__CHMP1A, b__M__GCST90200364, se__M__GCST90200364)])))
H[, row := .I]

# signals (definition)
sg <- fread(inp("LOC_DPEP1_16q24_signals")); stopifnot(nrow(sg) == 29L)
s5 <- sg[tid == "P__DPEP1"][order(signal)]$rsid
if (!identical(s5, SIG5)) stop("DPEP1 signals in signals.tsv differ from the pre-declaration: ", paste(s5, collapse = ","))
if (!identical(sg[tid == "P__CHMP1A"][order(signal)]$rsid, CHMP)) stop("CHMP1A signals differ")
say("signals.tsv: DPEP1 5 signals = %s (matches pre-declaration §5); CHMP1A = %s; cysteinylglycine = %s",
    paste(s5, collapse = ", "), paste(CHMP, collapse = ", "), paste(sg[tid == "M__GCST90200364"][order(signal)]$rsid, collapse = ", "))
key8 <- c(SIG5, SIG_B, CHMP)
stopifnot(all(key8 %in% H$rsid_ref), !anyDuplicated(H$rsid_ref))
say("key SNPs in harmonised: %s", paste(sprintf("%s(row %d, %s/%s, af %.4f, pal %d)", key8, match(key8, H$rsid_ref),
    H$a1[match(key8, H$rsid_ref)], H$a2[match(key8, H$rsid_ref)], H$af_ref[match(key8, H$rsid_ref)], H$palindromic[match(key8, H$rsid_ref)]), collapse = "; "))

# ---------------------------------------------------------------- 503 dosage extract, read IN PLACE (no copy)
gp <- inp("REF_DPEP1_REF503_I8")
gv <- readBin(gp, what = "integer", n = file.size(gp) + 1L, size = 1L, signed = TRUE)
stopifnot(length(gv) == 503L * 3575L)
G503 <- matrix(gv, nrow = 503L, byrow = TRUE); rm(gv)
stopifnot(all(G503 %in% 0:2))
afd <- max(abs(colMeans(G503) / 2 - H$af_ref))
say("503 i8 read in place: 503 x %d; A1 frequency vs harmonised af_ref max |diff| = %.2e (order check)", ncol(G503), afd)
stopifnot(afd < 1e-5)

# ---------------------------------------------------------------- raw re-derivation of the harmonised trait columns
COMP <- c(A = "T", T = "A", C = "G", G = "C")
ispal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")
read_prot <- function(path) {
  d <- fread(path, colClasses = list(character = c("ID", "ALLELE0", "ALLELE1")))
  d[, `:=`(pos38 = as.integer(pos), ea = toupper(ALLELE1), oa = toupper(ALLELE0), b = BETA, se = SE, eaf = A1FREQ)]
  d[, pos37 := as.integer(tstrsplit(ID, ":", fixed = TRUE)[[2]])]
  d[pos38 >= WIN38[1] & pos38 <= WIN38[2], .(pos38, pos37, ea, oa, b, se, eaf)]
}
read_met <- function(path) {
  d <- fread(path, colClasses = list(character = c("effect_allele", "other_allele", "rsid", "variant_id")))
  d[, .(pos38 = as.integer(base_pair_location), ea = toupper(effect_allele), oa = toupper(other_allele), b = beta,
        se = standard_error, eaf = suppressWarnings(as.numeric(effect_allele_frequency)), rsid)]
}
# orientation of (ea, oa) to the harmonised (a1, a2): +1 / -1 / 0 (exact letters only for palindromes; complement otherwise)
orient <- function(ea, oa, a1, a2) {
  s <- ifelse(ea == a1 & oa == a2, 1L, ifelse(ea == a2 & oa == a1, -1L, 0L))
  np <- !ispal(a1, a2) & s == 0L
  ce <- unname(COMP[ea]); co <- unname(COMP[oa])
  s[np] <- ifelse(ce[np] == a1[np] & co[np] == a2[np], 1L, ifelse(ce[np] == a2[np] & co[np] == a1[np], -1L, 0L))
  s[is.na(s)] <- 0L
  s
}
aset_of <- function(a, b) paste(pmin(a, b), pmax(a, b))
join_trait <- function(d, tid, by = "pos38") {
  d <- copy(d)[, aset := aset_of(ea, oa)]
  dk <- d[, .N, by = c(by, "aset")][N > 1]
  if (nrow(dk)) d <- d[!dk, on = c(by, "aset")]
  hh <- H[, c("row", by, "a1", "a2", "palindromic"), with = FALSE][, aset := aset_of(a1, a2)]
  m1 <- merge(hh, d, by = c(by, "aset"))
  dc <- copy(d)[, `:=`(ea = unname(COMP[ea]), oa = unname(COMP[oa]))][!is.na(ea) & !is.na(oa)][, aset := aset_of(ea, oa)]
  m2 <- merge(hh[palindromic == 0 & !row %in% m1$row], dc, by = c(by, "aset"))
  m <- rbind(m1, m2)
  if (anyDuplicated(m$row)) stop(tid, ": a harmonised row matched twice")
  m[, sgn := orient(ea, oa, a1, a2)]
  m <- merge(data.table(row = H$row), m, by = "row", all.x = TRUE)[order(row)]
  say("join %s -> harmonised on %s + allele set (complement allowed for non-palindromes): raw rows %d (duplicate keys dropped %d), harmonised %d, matched %d (%.2f%%), unoriented %d; raw null rates b %.4f se %.4f eaf %.4f",
      tid, by, nrow(d) + nrow(dk), nrow(dk), nrow(H), sum(!is.na(m$b)), 100 * mean(!is.na(m$b)), sum(m$sgn == 0, na.rm = TRUE),
      mean(is.na(d$b)), mean(is.na(d$se)), mean(is.na(d$eaf)))
  m[sgn == 0, `:=`(b = NA_real_, se = NA_real_)]
  m[, `:=`(b_al = sgn * b, eaf_a1 = ifelse(sgn == 1L, eaf, 1 - eaf))]
  m
}
chk <- list()
for (tid in c("P__DPEP1", "P__CHMP1A", "M__GCST90200364")) {
  path <- switch(tid, P__DPEP1 = inp("EXP_UKBPPP_EUR_DPEP1_OID20897"), P__CHMP1A = inp("EXP_UKBPPP_EUR_CHMP1A_OID20930"),
                 M__GCST90200364 = inp("EXP_CLSA_cysteinylglycine_GCST90200364"))
  d <- if (startsWith(tid, "P__")) read_prot(path) else read_met(path)
  m <- join_trait(d, tid)
  stopifnot(nrow(m) == nrow(H), all(m$row == H$row))
  db <- max(abs(m$b_al - H[[paste0("b__", tid)]]), na.rm = TRUE); ds <- max(abs(m$se - H[[paste0("se__", tid)]]), na.rm = TRUE)
  say("RE-DERIVATION %s: matched %d/%d; max |b_raw*orientation - harmonised b| = %.3g; max |se diff| = %.3g; rows unmatched %d",
      tid, sum(!is.na(m$b_al)), nrow(H), db, ds, sum(is.na(m$b_al)))
  if (sum(!is.na(m$b_al)) != nrow(H) || db > 1e-9 || ds > 1e-9) stop("harmonised ", tid, " does not re-derive from the raw window")
  chk[[tid]] <- m[, .(row, eaf_a1, n_rawrows = nrow(d))]
}
# GWAS frequencies at the key SNPs (palindrome orientation evidence)
for (tid in names(chk)) say("  %s A1 frequency at %s = %s (ref af %s)", tid, paste(key8, collapse = ","),
                           paste(sprintf("%.4f", chk[[tid]]$eaf_a1[match(key8, H$rsid_ref)]), collapse = ","),
                           paste(sprintf("%.4f", H$af_ref[match(key8, H$rsid_ref)]), collapse = ","))

# ---------------------------------------------------------------- SNP sets (rebuilt from the declared rules)
# susie = every harmonised SNP; p1np = non-palindromic; p1fr = p1np + palindromes resolvable by frequency in EVERY reporting
# GWAS at the locus (l3fm__robustness 'pal=freq' rule: ref MAF <= 0.40; GWAS MAF <= 0.40 and |GWAS A1 freq - ref| <= 0.10).
frq_all <- lapply(names(chk), function(t) chk[[t]]$eaf_a1); names(frq_all) <- names(chk)
for (tid in names(EXTRA)) {
  e <- EXTRA[[tid]]; assert_md5(e$p, e$m, tid)
  d <- if (e$kind == "protein") read_prot(e$p) else read_met(e$p)
  m <- join_trait(d, tid)
  pr <- H$palindromic == 1
  if (any(is.na(m$b_al[pr]))) stop(tid, " lacks a record at a palindromic harmonised row (engine refuses)")
  frq_all[[tid]] <- m$eaf_a1
}
pal <- which(H$palindromic == 1)
keep_fr <- vapply(pal, function(j) {
  f0 <- H$af_ref[j]
  if (min(f0, 1 - f0) > 0.40) return(FALSE)
  for (t in names(frq_all)) {
    f <- frq_all[[t]][j]
    if (is.na(f)) next                       # trait does not report a frequency here
    if (min(f, 1 - f) > 0.40 || abs(f - f0) > 0.10) return(FALSE)
  }
  TRUE
}, logical(1))
SETS <- list(susie = H$row, p1np = H$row[H$palindromic == 0], p1fr = sort(c(H$row[H$palindromic == 0], pal[keep_fr])))
say("SNP sets (rebuilt): susie %d; p1np %d; p1fr %d (palindromic %d, kept by the frequency rule %d)", length(SETS$susie),
    length(SETS$p1np), length(SETS$p1fr), length(pal), sum(keep_fr))
for (k in key8) say("  %s in p1np %s, p1fr %s", k, H[rsid_ref == k]$row %in% SETS$p1np, H[rsid_ref == k]$row %in% SETS$p1fr)

# ---------------------------------------------------------------- PLINK LD, panel 503 (paper 6 g1000_eur)
bed503 <- inp("REF_LD503_BED")
pre503 <- sub("\\.bed$", "", bed503)
assert_md5(paste0(pre503, ".bim"), res("LD503_bim")$md5, "LD503_bim"); assert_md5(paste0(pre503, ".fam"), res("LD503_fam")$md5, "LD503_fam")
LO <- min(H$pos37); HI <- max(H$pos37)
wb <- file.path(LDD, "win503")
if (!file.exists(paste0(wb, ".bim")))
  run_plink(c("--bfile", q(pre503), "--chr", "16", "--from-bp", LO, "--to-bp", HI, "--make-just-bim", "--out", q(wb)), "win503_bim")
B503 <- fread(paste0(wb, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "A1", "A2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
stopifnot(nrow(B503) > 3575L)
B503[, aset := paste(pmin(A1, A2), pmax(A1, A2))]
hk <- H[, .(row, pos37, aset = paste(pmin(a1, a2), pmax(a1, a2)), a1, rsid_ref)]
mb <- merge(hk, B503, by.x = c("pos37", "aset"), by.y = c("pos", "aset"))
say("503 window bim: %d variants in [%d, %d]; harmonised matched on pos37+allele set: %d/%d (dups %d); rsID equal %d; A1==a1 %d",
    nrow(B503), LO, HI, nrow(mb), nrow(H), sum(duplicated(mb$row)), sum(mb$rsid_ref == mb$id), sum(mb$A1 == mb$a1))
stopifnot(nrow(mb) == 3575L, !anyDuplicated(mb$row), all(mb$rsid_ref == mb$id), !anyDuplicated(B503$id[B503$id %in% mb$id]))

ld_square <- function(prefix, ids, tag, extra = character(0)) {
  lf <- file.path(LDD, paste0(tag, "_extract.txt")); writeLines(ids, lf)
  ob <- file.path(LDD, tag)
  if (!file.exists(paste0(ob, ".ld.bin")) || !file.exists(paste0(ob, ".snplist")))
    run_plink(c("--bfile", q(prefix), "--chr", "16", "--from-bp", LO, "--to-bp", HI, "--extract", q(lf), "--keep-allele-order",
                "--r", "square", "bin", "--write-snplist", "--freq", extra, "--out", q(ob)), tag)
  sl <- readLines(paste0(ob, ".snplist"))
  n <- length(sl)
  con <- file(paste0(ob, ".ld.bin"), "rb"); x <- readBin(con, "double", n = n * n + 1L, size = 8L); close(con)
  if (length(x) != n * n) stop(tag, ": ld.bin has ", length(x), " values, expected ", n * n)
  M <- matrix(x, n, n, byrow = TRUE); dimnames(M) <- list(sl, sl)
  # ASSERT every requested cell returned and non-NA
  miss <- setdiff(ids, sl)
  ok <- length(miss) == 0 && all(is.finite(M))
  say("PLINK %s: requested %d ids, returned %d, missing %d, non-finite cells %d, max|r| offdiag %.4f, symmetric %s -> %s",
      tag, length(ids), n, length(miss), sum(!is.finite(M)), max(abs(M[upper.tri(M)])), isTRUE(all.equal(M, t(M))), if (ok) "PASS" else "FAIL")
  list(M = M, ok = ok, missing = miss, snplist = sl, frq = fread(paste0(ob, ".frq")))
}
L503 <- ld_square(pre503, mb$id, "ld503_window")
if (!L503$ok) stop("503 LD assertion failed")
o <- match(H$rsid_ref, L503$snplist); stopifnot(!anyNA(o))
R503 <- L503$M[o, o]
sgn503 <- ifelse(mb$A1[match(H$row, mb$row)] == H$a1, 1, -1)
R503 <- R503 * outer(sgn503, sgn503); dimnames(R503) <- list(H$rsid_ref, H$rsid_ref); diag(R503) <- 1
Ri8 <- suppressWarnings(cor(G503))
dR <- max(abs(Ri8 - R503))
say("503 PLINK r (aligned to harmonised a1; sign flips %d) vs Pearson r of the 503 i8 A1 dosages: max |diff| = %.2e over %d cells",
    sum(sgn503 < 0), dR, length(R503))
stopifnot(dR < 1e-6)
rm(Ri8)

# ---------------------------------------------------------------- PLINK LD, panel 489 (1000G.EUR.QC.16)
bed489 <- inp("REF_LD489_chr16_BED")
pre489 <- sub("\\.bed$", "", bed489)
assert_md5(paste0(pre489, ".bim"), res("LD489_chr16_bim")$md5, "LD489_chr16_bim")
assert_md5(paste0(pre489, ".fam"), "fd9673c9469e18b200737b4e5f827c5b", "LD489_chr16_fam (RESOURCES REF_LD489_chr16_BED notes)")
B489 <- fread(paste0(pre489, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "A1", "A2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
stopifnot(nrow(B489) == 316981L)
B489[, aset := paste(pmin(A1, A2), pmax(A1, A2))]
B489w <- B489[pos >= LO & pos <= HI]
m2 <- merge(hk, B489w, by.x = c("pos37", "aset"), by.y = c("pos", "aset"))
dup2 <- m2[, .N, by = row][N > 1]$row
m2 <- m2[!row %in% dup2]
say("489 bim: %d variants on chr16, %d in window; harmonised matched on pos37+allele set: %d/%d (ambiguous %d); rsID equal %d; A1==a1 %d; ids unique in window %s",
    nrow(B489), nrow(B489w), nrow(m2), nrow(H), length(dup2), sum(m2$rsid_ref == m2$id), sum(m2$A1 == m2$a1), !anyDuplicated(B489w$id[B489w$id %in% m2$id]))
stopifnot(!anyDuplicated(B489w$id[B489w$id %in% m2$id]))
L489 <- ld_square(pre489, m2$id, "ld489_window", extra = "--missing")
if (length(L489$missing)) say("489: %d requested ids not returned: %s", length(L489$missing), paste(head(L489$missing, 10), collapse = ","))
lm <- fread(file.path(LDD, "ld489_window.lmiss"))
fr2 <- L489$frq
m2 <- m2[id %in% L489$snplist]
m2[, sgn := ifelse(A1 == a1, 1, -1)]
m2[, f2_A1bim := fr2$MAF[match(id, fr2$SNP)] * ifelse(fr2$A1[match(id, fr2$SNP)] == A1, 1, NA)]
# PLINK --freq with --keep-allele-order reports MAF as the frequency of the A1 column of the bim
stopifnot(all(fr2$A1[match(m2$id, fr2$SNP)] == m2$A1))
m2[, f2_a1 := ifelse(sgn == 1, f2_A1bim, 1 - f2_A1bim)]
m2[, nmiss := lm$N_MISS[match(id, lm$SNP)]]
m2[, af_ref := H$af_ref[row]]
m2[, valid := nmiss == 0 & pmin(f2_a1, 1 - f2_a1) >= 0.01 & abs(f2_a1 - af_ref) <= 0.05 & is.finite(f2_a1)]
say("489 panel-2 validity (p2np rule: no missing genotype, panel-2 MAF >= 0.01, |A1 freq panel2 - panel1| <= 0.05): valid %d of %d matched; drops missing %d, maf %d, af %d",
    sum(m2$valid), nrow(m2), sum(m2$nmiss > 0), sum(pmin(m2$f2_a1, 1 - m2$f2_a1) < 0.01), sum(abs(m2$f2_a1 - m2$af_ref) > 0.05))
v2 <- m2[valid == TRUE][order(row)]
o2 <- match(v2$id, L489$snplist)
R489 <- L489$M[o2, o2] * outer(v2$sgn, v2$sgn); dimnames(R489) <- list(H$rsid_ref[v2$row], H$rsid_ref[v2$row]); diag(R489) <- 1
stopifnot(all(is.finite(R489)))
oa2 <- match(m2[order(row)]$id, L489$snplist); sa2 <- m2[order(row)]$sgn
R489all <- L489$M[oa2, oa2] * outer(sa2, sa2); dimnames(R489all) <- list(H$rsid_ref[m2[order(row)]$row], H$rsid_ref[m2[order(row)]$row]); diag(R489all) <- 1
SETS$p2np <- intersect(SETS$p1np, v2$row)
say("p2np (rebuilt) = p1np SNPs valid in panel 2: %d", length(SETS$p2np))
for (k in key8) say("  %s valid in 489: %s", k, H[rsid_ref == k]$row %in% v2$row)

# ---------------------------------------------------------------- small PLINK runs for the key SNPs + MUTATION test
small_ld <- function(prefix, ids, tag) {
  lf <- file.path(LDD, paste0(tag, "_extract.txt")); writeLines(ids, lf); ob <- file.path(LDD, tag)
  run_plink(c("--bfile", q(prefix), "--extract", q(lf), "--keep-allele-order", "--r", "square", "--write-snplist", "--out", q(ob)), tag)
  sl <- if (file.exists(paste0(ob, ".snplist"))) readLines(paste0(ob, ".snplist")) else character(0)
  M <- if (file.exists(paste0(ob, ".ld")) && length(sl)) as.matrix(read.table(paste0(ob, ".ld"))) else matrix(numeric(0), 0, 0)
  if (length(sl)) dimnames(M) <- list(sl, sl)
  ok <- length(setdiff(ids, sl)) == 0 && length(M) == length(ids)^2 && all(is.finite(M))
  list(M = M, ok = ok, missing = setdiff(ids, sl))
}
S503 <- small_ld(pre503, key8, "ld503_key8")
S489 <- small_ld(pre489, key8, "ld489_key8")
say("ASSERT key-8 PLINK (square text): 503 %s, 489 %s", if (S503$ok) "PASS" else "FAIL", if (S489$ok) "PASS" else "FAIL")
stopifnot(S503$ok, S489$ok)
a8 <- function(S, R, rows, sg) { o <- H$rsid_ref[rows]; M <- S$M[o, o] * outer(sg, sg); max(abs(M - R[o, o])) }
k8rows <- match(key8, H$rsid_ref)
k8in2 <- k8rows[k8rows %in% v2$row]
say("key-8 square (text, 6 d.p.) vs window bin matrix: 503 max|diff| %.2e (8 SNPs); 489 max|diff| %.2e (%d of 8 SNPs valid in panel 2)",
    a8(S503, R503, k8rows, sgn503[k8rows]), a8(S489, R489, k8in2, v2$sgn[match(k8in2, v2$row)]), length(k8in2))
mut <- key8; mut[1] <- paste0(mut[1], "_MUTATED")
SM <- small_ld(pre503, mut, "ld503_key8_MUTATION")
say("MUTATION TEST (id %s -> %s): assertion %s (missing: %s) -- the assertion %s", key8[1], mut[1], if (SM$ok) "PASS" else "FAIL",
    paste(SM$missing, collapse = ","), if (!SM$ok) "TRIPS as required" else "DID NOT TRIP (defect)")
stopifnot(!SM$ok)

# ---------------------------------------------------------------- eGFR: lift the GRCh38 window to GRCh37, stream-cut CKDGen
ch <- inp("REF_LIFTOVER_hg38ToHg19")
chp <- file.path(WORK, "hg38ToHg19.over.chain")
if (!file.exists(chp)) { zz <- gzfile(ch, "rt"); writeLines(readLines(zz), chp); close(zz) }
suppressMessages(library(rtracklayer))
chain <- import.chain(chp)
lift1 <- function(p) { g <- GenomicRanges::GRanges("chr16", IRanges::IRanges(p, width = 1)); r <- unlist(liftOver(g, chain)); if (length(r) != 1) stop("liftover not unique for ", p); start(r) }
win37 <- c(lift1(WIN38[1]), lift1(WIN38[2]))
say("eGFR window: GRCh38 %d-%d lifted (hg38ToHg19 chain) to GRCh37 %d-%d; harmonised pos37 range %d-%d (inside: %s)",
    WIN38[1], WIN38[2], win37[1], win37[2], LO, HI, LO >= win37[1] && HI <= win37[2])
stopifnot(LO >= win37[1], HI <= win37[2])
eg_cache <- file.path(WORK, "egfr_window_raw.tsv")
egp <- inp("OUT_CKDGEN2019_eGFR_EA")
if (!file.exists(eg_cache)) {
  con <- gzfile(egp, "rt"); hdr <- readLines(con, 1); ntot <- 0L; keep <- character(0)
  repeat {
    x <- readLines(con, 500000L); if (!length(x)) break
    ntot <- ntot + length(x)
    y <- x[startsWith(x, "16 ")]
    if (length(y)) { p <- as.integer(vapply(strsplit(y, " ", fixed = TRUE), `[`, "", 2)); keep <- c(keep, y[p >= win37[1] & p <= win37[2]]) }
  }
  close(con)
  say("CKDGen stream: header '%s'; %d data rows read (RESOURCES n_info: 8885712); %d rows in chr16 window", hdr, ntot, length(keep))
  stopifnot(ntot == 8885712L)
  writeLines(c(gsub(" ", "\t", hdr), gsub(" +", "\t", keep)), eg_cache)
}
EG <- fread(eg_cache, colClasses = list(character = c("RSID", "Allele1", "Allele2")))
setnames(EG, c("Chr", "Pos_b37", "RSID", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "n_total_sum"),
         c("chr", "pos37", "rsid", "ea", "oa", "eaf", "b", "se", "p", "n"))
EG[, `:=`(ea = toupper(ea), oa = toupper(oa))]
say("eGFR window rows %d; per-column null rates: b %.4f se %.4f eaf %.4f n %.4f; non-numeric tokens: %d", nrow(EG),
    mean(is.na(EG$b)), mean(is.na(EG$se)), mean(is.na(EG$eaf)), mean(is.na(EG$n)), sum(is.na(suppressWarnings(as.numeric(EG$b)))))
hk2 <- H[, .(row, pos37, a1, a2, af_ref, palindromic, rsid_ref)][, aset := aset_of(a1, a2)]
EG[, aset := aset_of(ea, oa)]
dkE <- EG[, .N, by = .(pos37, aset)][N > 1]
if (nrow(dkE)) EG <- EG[!dkE, on = c("pos37", "aset")]
e1 <- merge(hk2, EG, by = c("pos37", "aset"))
EGc <- copy(EG)[, `:=`(ea = unname(COMP[ea]), oa = unname(COMP[oa]))][!is.na(ea) & !is.na(oa)][, aset := aset_of(ea, oa)]
e2 <- merge(hk2[palindromic == 0 & !row %in% e1$row], EGc, by = c("pos37", "aset"))
E <- rbind(e1, e2)[order(row)]
stopifnot(!anyDuplicated(E$row))
E[, s_let := orient(ea, oa, a1, a2)]                       # letters (complement allowed for non-palindromes)
E[, f_let := ifelse(s_let == 1L, eaf, ifelse(s_let == -1L, 1 - eaf, NA_real_))]   # A1 frequency implied by the letters
# palindromes: the letters cannot tell swap from strand flip, so orientation is decided by frequency (ref MAF <= 0.40 holds
# for every harmonised palindrome; tolerance 0.20 = prep AF_DIFF_MAX); exactly one of the two readings must fit
E[, sgn := s_let]; E[, f_a1 := f_let]
E[palindromic == 1, `:=`(
  sgn = ifelse(abs(f_let - af_ref) <= 0.2 & abs((1 - f_let) - af_ref) > 0.2, s_let,
        ifelse(abs((1 - f_let) - af_ref) <= 0.2 & abs(f_let - af_ref) > 0.2, -s_let, 0L)),
  f_a1 = ifelse(abs(f_let - af_ref) <= 0.2 & abs((1 - f_let) - af_ref) > 0.2, f_let,
         ifelse(abs((1 - f_let) - af_ref) <= 0.2 & abs(f_let - af_ref) > 0.2, 1 - f_let, NA_real_)))]
E[, afbad := sgn != 0L & (!is.finite(f_a1) | abs(f_a1 - af_ref) > 0.2)]
say("eGFR join to harmonised (pos37 + allele set, complement for non-palindromes): window rows %d (duplicate keys dropped %d); matched %d/%d (%.1f%%); unoriented %d (palindromes undecidable by frequency %d); |A1 freq - ref| > 0.20 dropped %d; rsID equal %d of %d",
    nrow(EG) + nrow(dkE), nrow(dkE), nrow(E), nrow(H), 100 * nrow(E) / nrow(H), sum(E$sgn == 0L), sum(E$sgn == 0L & E$palindromic == 1),
    sum(E$afbad), sum(E$rsid == E$rsid_ref), nrow(E))
say("eGFR palindromic orientation: %d palindromes matched; letters kept %d, strand-flipped by frequency %d, undecidable %d",
    sum(E$palindromic == 1), sum(E$palindromic == 1 & E$sgn == E$s_let & E$sgn != 0L), sum(E$palindromic == 1 & E$sgn == -E$s_let & E$sgn != 0L),
    sum(E$palindromic == 1 & E$sgn == 0L))
E <- E[sgn != 0L & !afbad]
E[, `:=`(b_al = sgn * b, z = sgn * b / se)]
EGFR <- data.table(row = H$row)
EGFR[E, on = "row", `:=`(b = i.b_al, se = i.se, n = i.n, f_a1 = i.f_a1)]
EGFR[, z := b / se]
say("eGFR usable at %d/%d harmonised SNPs; p1np %d/%d; at the key SNPs: %s", sum(!is.na(EGFR$b)), nrow(H),
    sum(!is.na(EGFR$b[SETS$p1np])), length(SETS$p1np),
    paste(sprintf("%s b=%.6g se=%.4g n=%d", key8, EGFR$b[match(key8, H$rsid_ref)], EGFR$se[match(key8, H$rsid_ref)], EGFR$n[match(key8, H$rsid_ref)]), collapse = "; "))
stopifnot(all(!is.na(EGFR$b[match(SIG5, H$rsid_ref)])))
N_EGFR_MED <- median(EGFR$n, na.rm = TRUE)
say("eGFR N: per-variant n_total_sum; window median %.0f (min %d, max %d) used as scalar N for coloc sdY estimation; RESOURCES n_info max 567460 (whole file)",
    N_EGFR_MED, min(EGFR$n, na.rm = TRUE), max(EGFR$n, na.rm = TRUE))

# ---------------------------------------------------------------- N for protein and metabolite (logged; never typed into a derived value)
N_PROT <- 33375; N_CYS <- 8299
nhat <- function(se, f) { ok <- is.finite(f) & pmin(f, 1 - f) > 0.05 & se > 0; median(1 / (2 * f[ok] * (1 - f[ok]) * se[ok]^2)) }
NHAT <- c(P__DPEP1 = nhat(H$se__P__DPEP1, chk$P__DPEP1$eaf_a1), M__GCST90200364 = nhat(H$se__M__GCST90200364, chk$M__GCST90200364$eaf_a1),
          P__CHMP1A = nhat(H$se__P__CHMP1A, chk$P__CHMP1A$eaf_a1))
say("N used: DPEP1 protein %d (RESOURCES n_info: UKB-PPP samplesize at the instrument); cysteinylglycine %d (RESOURCES n_info: N 8,299, ledger L326; the CLSA window has no N column). Computed here as a sensitivity, n_hat = median 1/(2f(1-f)se^2) over GWAS MAF > 0.05: DPEP1 %.0f, cysteinylglycine %.0f, CHMP1A %.0f",
    N_PROT, N_CYS, NHAT["P__DPEP1"], NHAT["M__GCST90200364"], NHAT["P__CHMP1A"])

saveRDS(list(H = H, SETS = SETS, R503 = R503, R489 = R489, R489all = R489all, m2 = m2, v2rows = v2$row, sgn489 = v2$sgn, EGFR = EGFR, N_EGFR_MED = N_EGFR_MED,
             N_PROT = N_PROT, N_CYS = N_CYS, NHAT = NHAT, frq = list(DPEP1 = chk$P__DPEP1$eaf_a1, CYS = chk$M__GCST90200364$eaf_a1, CHMP1A = chk$P__CHMP1A$eaf_a1),
             win37 = win37, pal_keep = data.table(row = pal, keep_fr = keep_fr)), file.path(WORK, "prep.rds"))
say("==== S5_01_prepare done in %.0f s ====", as.numeric(difftime(Sys.time(), t0, units = "secs")))
