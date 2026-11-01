# 03_windows.R -- protein windows (UKB-PPP EUR cis files) and outcome window cuts for the 16 candidates.
# GRCh37 bounds = engine rule (l3fm__prep.py): min/max GRCh37 position in the UKB-PPP IDs inside the GRCh38 window, +/- 1000;
# cross-checked against a hg38ToHg19 liftover of the window bounds. Large outcome files are streamed, never loaded whole.
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
suppressMessages({library(rtracklayer); library(GenomicRanges); library(Rsamtools)})
say("==== W6b 03_windows ====")
cand <- fread(file.path(WORK, "candidates.tsv"), colClasses = list(character = "chr"))
EXD <- file.path(WORK, "extracts"); dir.create(EXD, showWarnings = FALSE)
MAN <- file.path(EXD, "EXTRACT_MANIFEST.tsv")
man <- if (file.exists(MAN)) fread(MAN, colClasses = "character") else data.table(file = character(), md5 = character(), rows = character(), source = character(), source_md5 = character(), note = character())
man_add <- function(f, src, src_md5, note) {
  man <<- man[file != basename(f)]
  man <<- rbind(man, data.table(file = basename(f), md5 = md5(f), rows = as.character(nrow(fread(f))), source = src, source_md5 = src_md5, note = note))
  fwrite(man, MAN, sep = "\t")
}

# ---- protein windows: gene -> OID via the UKB-PPP instrument table ----
it <- fread(res("EXP_UKBPPP_INSTRUMENT_TABLE"))
known <- c(UMOD = "EXP_UKBPPP_EUR_UMOD_OID20237", DPEP1 = "EXP_UKBPPP_EUR_DPEP1_OID20897", CELSR2 = "EXP_UKBPPP_EUR_CELSR2_OID30593")
pqdir <- file.path(ROOT, "weight-loss.omics/data/proteome/pqtl_regional")
cand[, oid := it$oid[match(gene, it$gene)]]
stopifnot(!anyNA(cand$oid), all(it[gene %in% cand$gene, .N, by = gene]$N == 1))
cand[, ppath := file.path(pqdir, sprintf("%s_chr%s_cis.tsv", oid, chr))]
cand[, protein_local := file.exists(ppath)]
say("protein windows local: %s; missing: %s", paste(unique(cand[protein_local == TRUE]$gene), collapse = ","), paste(unique(cand[protein_local == FALSE]$gene), collapse = ","))

read_protein <- function(path, ch, lo, hi) {
  d <- fread(path, colClasses = list(character = c("ID", "chrom", "ALLELE0", "ALLELE1")))
  stopifnot(identical(names(d), c("ID", "chrom", "pos", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P")), nrow(d) > 0)
  if (any(d$chrom != ch)) stop(path, " carries another chromosome")
  n0 <- nrow(d)
  d <- d[pos >= lo & pos <= hi]
  idf <- tstrsplit(d$ID, ":", fixed = TRUE)
  if (any(idf[[1]] != ch)) stop("malformed ID chromosome in ", path)
  if (any(toupper(idf[[3]]) != toupper(d$ALLELE0) | toupper(idf[[4]]) != toupper(d$ALLELE1))) stop("ID alleles differ from ALLELE0/1 in ", path)
  out <- data.table(pos38 = d$pos, pos37 = as.integer(idf[[2]]), ea = toupper(d$ALLELE1), oa = toupper(d$ALLELE0),
                    eaf = d$A1FREQ, b = d$BETA, se = d$SE, log10p = d$LOG10P)
  bad <- !is.finite(out$b) | !is.finite(out$se) | out$se <= 0 | !is.finite(out$eaf)
  if (any(bad)) say("  %s: %d rows with non-finite b/se/eaf dropped", basename(path), sum(bad))
  out <- out[!bad]
  out[, akey := paste(pos38, pmin(ea, oa), pmax(ea, oa))]
  dup <- out[duplicated(akey) | duplicated(akey, fromLast = TRUE)]
  out <- out[!akey %in% dup$akey]
  say("  protein %s: %d rows in file, %d in window, %d duplicate-key rows dropped, null rates b/se/eaf 0/0/0 after drop", basename(path), n0, nrow(d), nrow(dup))
  out[, akey := NULL]
  out[]
}

chain_txt <- file.path(WORK, "hg38ToHg19.over.chain")          # decompressed copy of REF_LIFTOVER_hg38ToHg19 (import.chain reads plain text)
if (!file.exists(chain_txt)) { con <- gzfile(res("REF_LIFTOVER_hg38ToHg19")); writeLines(readLines(con), chain_txt); close(con) }
chain19 <- import.chain(chain_txt)
lift19 <- function(ch, pos) {
  g <- GRanges(paste0("chr", ch), IRanges(pos, pos))
  l <- liftOver(g, chain19)
  vapply(l, function(x) if (length(x) == 1) start(x) else NA_integer_, 1L)
}

prot <- list()
for (i in seq_len(nrow(cand))) {
  cc <- cand[i]
  if (!cc$protein_local) { cand[i, `:=`(lo37 = NA_integer_, hi37 = NA_integer_)]; next }
  pth <- if (cc$gene %in% names(known)) res(known[[cc$gene]]) else reg_new(cc$ppath, sprintf("UKB-PPP EUR cis window for %s (%s), section 6(d) candidate %s", cc$gene, cc$oid, cc$cid))
  p <- read_protein(pth, cc$chr, cc$lo38, cc$hi38)
  prot[[cc$cid]] <- p
  fp <- file.path(EXD, paste0(cc$cid, "__protein.tsv.gz"))
  fwrite(p, fp, sep = "\t")
  man_add(fp, cc$ppath, md5c(pth), sprintf("UKB-PPP %s rows in GRCh38 %s:%d-%d", cc$oid, cc$chr, cc$lo38, cc$hi38))
  cand[i, `:=`(lo37 = min(p$pos37) - 1000L, hi37 = max(p$pos37) + 1000L, protein_md5 = md5c(pth), protein_path = pth)]
}
# liftover cross-check of the GRCh37 bounds (and the GRCh37 region for a candidate with no protein window)
cand[, lo37_lift := lift19(chr, lo38)]; cand[, hi37_lift := lift19(chr, hi38)]
for (i in seq_len(nrow(cand))) say("  %-12s GRCh38 %s:%d-%d | GRCh37 from IDs %s-%s | lifted %s-%s", cand$cid[i], cand$chr[i], cand$lo38[i], cand$hi38[i],
                                   cand$lo37[i], cand$hi37[i], cand$lo37_lift[i], cand$hi37_lift[i])
cand[is.na(lo37), `:=`(lo37 = lo37_lift - 1000L, hi37 = hi37_lift + 1000L, bounds37_source = "hg38ToHg19 liftover of the 29-log window bounds +/- 1000 (no protein window)")]
cand[is.na(bounds37_source), bounds37_source := "UKB-PPP ID GRCh37 positions in the window +/- 1000 (engine rule)"]
um <- cand[cid == "UMOD_eGFR"]
ok <- um$lo37 == 19841947L && um$hi37 == 20868706L
log_assert("S6.d.umod_window37", "19841947-20868706 (UMOD PREP_STAMP window_grch37_ref, definition)", paste(um$lo37, um$hi37, sep = "-"), ok,
           "engine rule re-implemented: min/max ID pos37 in GRCh38 window +/- 1000")
dp <- cand[cid == "DPEP1_eGFR"]
log_assert("S6.d.dpep1_window37", "89178787-90189467 (DPEP1 PREP_STAMP window_grch37_ref, definition)", paste(dp$lo37, dp$hi37, sep = "-"),
           dp$lo37 == 89178787L && dp$hi37 == 90189467L, "same rule, DPEP1 window")

# ---- outcome cuts ----
cut_ckdgen <- function(path, src_id, wins) {   # wins: data.table(cid, chr, lo37, hi37); whitespace-separated, GRCh37
  outs <- file.path(EXD, paste0(wins$cid, "__outcome.tsv.gz"))
  if (all(file.exists(outs)) && all(basename(outs) %in% man$file)) { say("  %s cuts exist (skip)", src_id); return(invisible()) }
  con <- gzfile(path, "r"); on.exit(close(con))
  hdr <- strsplit(readLines(con, 1), "\\s+")[[1]]
  want <- c("Chr", "Pos_b37", "RSID", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "n_total_sum")
  if (!identical(hdr, want)) stop(src_id, " header ", paste(hdr, collapse = " "))
  acc <- list(); nread <- 0; t0 <- Sys.time()
  pref <- paste0(unique(wins$chr), " ")
  repeat {
    x <- readLines(con, n = 500000L)
    if (!length(x)) break
    nread <- nread + length(x)
    if (any(grepl("\t", x, fixed = TRUE))) stop(src_id, ": a row contains a tab")
    keep <- Reduce(`|`, lapply(pref, function(p) startsWith(x, p)))
    x <- x[keep]
    if (length(x)) {
      f <- tstrsplit(x, "\\s+")
      d <- data.table(chr = f[[1]], pos = as.integer(f[[2]]), rsid = f[[3]], ea = toupper(f[[4]]), oa = toupper(f[[5]]), eaf = as.numeric(f[[6]]),
                      b = as.numeric(f[[7]]), se = as.numeric(f[[8]]), p = as.numeric(f[[9]]), n = as.numeric(f[[10]]))
      for (w in seq_len(nrow(wins))) {
        dd <- d[chr == wins$chr[w] & pos >= wins$lo37[w] & pos <= wins$hi37[w]]
        if (nrow(dd)) acc[[wins$cid[w]]] <- rbind(acc[[wins$cid[w]]], dd)
      }
    }
    if (nread %% 2000000 == 0) say("  %s: %d lines streamed (%.0fs)", src_id, nread, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  say("  %s: %d data lines streamed in %.0fs", src_id, nread, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  for (w in seq_len(nrow(wins))) {
    d <- acc[[wins$cid[w]]]
    if (is.null(d) || !nrow(d)) stop("empty cut for ", wins$cid[w])
    nn <- sapply(d, function(v) mean(is.na(v)))
    say("  cut %-12s %s rows %d; null rates %s", wins$cid[w], src_id, nrow(d), paste(sprintf("%s=%.3f", names(nn), nn), collapse = " "))
    d[, pos_build := "GRCh37"]
    f <- file.path(EXD, paste0(wins$cid[w], "__outcome.tsv.gz"))
    fwrite(d, f, sep = "\t")
    man_add(f, src_id, md5c(path), sprintf("GRCh37 %s:%d-%d; %d lines streamed", wins$chr[w], wins$lo37[w], wins$hi37[w], nread))
  }
}
fg_parse <- function(x) {
  f <- tstrsplit(x, "\t", fixed = TRUE)
  data.table(chr = f[[1]], pos = as.integer(f[[2]]), oa = toupper(f[[3]]), ea = toupper(f[[4]]), rsid = f[[5]], p = as.numeric(f[[7]]),
             b = as.numeric(f[[9]]), se = as.numeric(f[[10]]), eaf = as.numeric(f[[11]]), n = NA_real_)
}
FGHDR <- c("#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls")
cut_finngen_stream <- function(path, src_id, wins) {   # wins on GRCh38 (chr, lo38, hi38); file sorted by chrom, pos
  outs <- file.path(EXD, paste0(wins$cid, "__outcome.tsv.gz"))
  if (all(file.exists(outs)) && all(basename(outs) %in% man$file)) { say("  %s cuts exist (skip)", src_id); return(invisible()) }
  con <- gzfile(path, "r"); on.exit(close(con))
  hdr <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]]
  if (!identical(hdr, FGHDR)) stop(src_id, " header ", paste(hdr, collapse = " "))
  acc <- list(); nread <- 0; t0 <- Sys.time(); done <- setNames(rep(FALSE, nrow(wins)), wins$cid)
  pref <- paste0(unique(wins$chr), "\t")
  repeat {
    x <- readLines(con, n = 1000000L)
    if (!length(x)) break
    nread <- nread + length(x)
    keep <- Reduce(`|`, lapply(pref, function(p) startsWith(x, p)))
    if (any(keep)) {
      d <- fg_parse(x[keep])
      for (w in seq_len(nrow(wins))) {
        dd <- d[chr == wins$chr[w] & pos >= wins$lo38[w] & pos <= wins$hi38[w]]
        if (nrow(dd)) acc[[wins$cid[w]]] <- rbind(acc[[wins$cid[w]]], dd)
        if (any(d$chr == wins$chr[w] & d$pos > wins$hi38[w])) done[w] <- TRUE
      }
    }
    if (all(done)) { say("  %s: all windows passed after %d lines (sorted file; stop)", src_id, nread); break }
    if (nread %% 5000000 == 0) say("  %s: %d lines streamed (%.0fs)", src_id, nread, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  say("  %s: %d data lines streamed in %.0fs", src_id, nread, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  for (w in seq_len(nrow(wins))) {
    d <- acc[[wins$cid[w]]]
    if (is.null(d) || !nrow(d)) stop("empty cut for ", wins$cid[w])
    nn <- sapply(d, function(v) mean(is.na(v) | (is.character(v) & v == "")))
    say("  cut %-12s %s rows %d; null/empty rates %s", wins$cid[w], src_id, nrow(d), paste(sprintf("%s=%.3f", names(nn), nn), collapse = " "))
    d[, pos_build := "GRCh38"]
    f <- file.path(EXD, paste0(wins$cid[w], "__outcome.tsv.gz"))
    fwrite(d, f, sep = "\t")
    man_add(f, src_id, md5c(path), sprintf("GRCh38 %s:%d-%d; stream", wins$chr[w], wins$lo38[w], wins$hi38[w]))
  }
}
cut_finngen_tabix <- function(path, tbi, src_id, wins) {
  outs <- file.path(EXD, paste0(wins$cid, "__outcome.tsv.gz"))
  if (all(file.exists(outs)) && all(basename(outs) %in% man$file)) { say("  %s cuts exist (skip)", src_id); return(invisible()) }
  tf <- TabixFile(path, index = tbi)
  for (w in seq_len(nrow(wins))) {
    x <- scanTabix(tf, param = GRanges(wins$chr[w], IRanges(wins$lo38[w], wins$hi38[w])))[[1]]
    d <- fg_parse(x)
    d <- d[pos >= wins$lo38[w] & pos <= wins$hi38[w]]
    nn <- sapply(d, function(v) mean(is.na(v) | (is.character(v) & v == "")))
    say("  cut %-12s %s rows %d (tabix); null/empty rates %s", wins$cid[w], src_id, nrow(d), paste(sprintf("%s=%.3f", names(nn), nn), collapse = " "))
    if (!nrow(d)) stop("empty tabix cut ", wins$cid[w])
    d[, pos_build := "GRCh38"]
    f <- file.path(EXD, paste0(wins$cid[w], "__outcome.tsv.gz"))
    fwrite(d, f, sep = "\t")
    man_add(f, src_id, md5c(path), sprintf("GRCh38 %s:%d-%d; tabix", wins$chr[w], wins$lo38[w], wins$hi38[w]))
  }
}

cut_ckdgen(res("OUT_CKDGEN2019_eGFR_EA"), "OUT_CKDGEN2019_eGFR_EA", cand[outcome == "eGFR", .(cid, chr, lo37, hi37)])
cut_ckdgen(res("OUT_CKDGEN2019_CKD_EA"), "OUT_CKDGEN2019_CKD_EA", cand[outcome == "CKD", .(cid, chr, lo37, hi37)])
fgT2D <- res("OUT_FG_R10_T2D")
tbi <- reg_new(paste0(fgT2D, ".tbi"), "tabix index of OUT_FG_R10_T2D (noted in its RESOURCES notes, md5 7c7965b8...)")
stopifnot(md5c(tbi) == "7c7965b8e6794d20e897be42a8666895")
cut_finngen_tabix(fgT2D, tbi, "OUT_FG_R10_T2D", cand[outcome == "T2D", .(cid, chr, lo38, hi38)])
cut_finngen_stream(res("OUT_FG_R10_I9_CHD"), "OUT_FG_R10_I9_CHD", cand[outcome == "CAD", .(cid, chr, lo38, hi38)])
cut_finngen_stream(res("OUT_FG_R10_I9_HEARTFAIL"), "OUT_FG_R10_I9_HEARTFAIL", cand[outcome == "HF", .(cid, chr, lo38, hi38)])

# UMOD eGFR cut must equal the pre-cut engine window (raw_input OUT_CKDGEN2019_eGFR_UMOD_window) row for row
pre <- fread(res("OUT_CKDGEN2019_eGFR_UMOD_window"))
mine <- fread(file.path(EXD, "UMOD_eGFR__outcome.tsv.gz"))
mm <- merge(pre[, .(pos = pos37, rsid, ea = toupper(ea), oa = toupper(oa), eaf, b, se, n)], mine[, .(pos, rsid, ea, oa, eaf2 = eaf, b2 = b, se2 = se, n2 = n)], by = c("pos", "rsid", "ea", "oa"))
okc <- nrow(pre) == nrow(mine) && nrow(mm) == nrow(pre) && max(abs(mm$b - mm$b2)) == 0 && max(abs(mm$se - mm$se2)) == 0 && max(abs(mm$n - mm$n2)) == 0
log_assert("S6.d.umod_egfr_cut_vs_engine_window", sprintf("%d rows identical", nrow(pre)), sprintf("%d mine, %d joined", nrow(mine), nrow(mm)), okc,
           "fresh stream cut of OUT_CKDGEN2019_eGFR_EA vs OUT_CKDGEN2019_eGFR_UMOD_window")

# sample sizes (rule 7: N_eff computed from case/control counts, never typed)
fgm <- fread(res("OUT_FG_R10_MANIFEST"))
nfg <- function(code) { r <- fgm[phenocode == code]; stopifnot(nrow(r) == 1); c(nca = r$num_cases, nco = r$num_controls) }
ckd_n <- res_entry("OUT_CKDGEN2019_CKD_EA")$n_info
ckd_cc <- as.numeric(regmatches(ckd_n, regexec("cases ([0-9]+) / controls ([0-9]+)", ckd_n))[[1]][2:3])
nout <- rbind(
  data.table(outcome = "T2D", nca = nfg("T2D")[["nca"]], nco = nfg("T2D")[["nco"]], source = "R10_manifest phenocode T2D"),
  data.table(outcome = "CAD", nca = nfg("I9_CHD")[["nca"]], nco = nfg("I9_CHD")[["nco"]], source = "R10_manifest phenocode I9_CHD"),
  data.table(outcome = "HF", nca = nfg("I9_HEARTFAIL")[["nca"]], nco = nfg("I9_HEARTFAIL")[["nco"]], source = "R10_manifest phenocode I9_HEARTFAIL"),
  data.table(outcome = "CKD", nca = ckd_cc[1], nco = ckd_cc[2], source = "RESOURCES n_info of OUT_CKDGEN2019_CKD_EA (ledger L311)"))
nout[, neff := neff(nca, nco)]
for (i in seq_len(nrow(nout))) say("N_eff %s = 4*%d*%d/(%d+%d) = %.0f (%s)", nout$outcome[i], nout$nca[i], nout$nco[i], nout$nca[i], nout$nco[i], nout$neff[i], nout$source[i])
fwrite(nout, file.path(WORK, "outcome_n.tsv"), sep = "\t")
fwrite(cand, file.path(WORK, "candidates_windows.tsv"), sep = "\t")
say("03_windows done")
