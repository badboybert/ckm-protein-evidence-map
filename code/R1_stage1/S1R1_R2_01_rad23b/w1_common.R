# W1 (Stage 1, pre-declaration section 1) -- shared helpers: logging, md5 assertions, resource lookup,
# streamed extraction from large summary-statistic files (never loaded whole), PLINK calls.
# Sourced by w1_extract.R and w1_main.R. Writes only under OUT.
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
suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(digest) })

OUT     <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_01_rad23b")
RES_DIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
STAGING <- paste0(.P31R_V3, "/V4_local/staging/l3-finemap")
RSCRIPT <- file.path(R.home("bin"), "Rscript")
for (d in c("extracts", "plink", "ref")) dir.create(file.path(OUT, d), showWarnings = FALSE, recursive = TRUE)

# ---- governing documents (md5s passed in the workflow prompt) ----
GOV <- list(
  briefs           = c(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/STAGE1_BRIEFS.md"), "40d0d52f7055542d774c03f1edf04e00"),
  prespec_redacted = c(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/P31R_R2_STAGE1_PRESPEC_REDACTED.md"), "a281c3464619b1b67c28e13e7392c41c"),
  amend_redacted   = c(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), "f265d600973c43239bb8377e0cd39d33"),
  prespec_full     = c(paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), "3319dc39740d9f5b3b8822133d9108b7"),
  amend_full       = c(paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), "377508832ea9bd384133e869c456c909"),
  gls_q            = c(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/gls_q.R"), "b7390fdf5b6f2f9c1d30c298d6afdafc"),
  resources        = c(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/RESOURCES.json"), "6d21238155a62fb58f71e05a01638dd0"),
  outcome_check_prespec = c(paste0(.P31R_V3, "/R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), "c9079f2de1ffe75ce31f8a1d597d1fda")
)

LOG <- Sys.getenv("W1_LOG", file.path(OUT, "run.log"))
logf <- function(fmt, ...) {
  s <- if (length(list(...))) sprintf(fmt, ...) else fmt
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", s)
  cat(line, "\n", file = LOG, append = TRUE, sep = "")
  cat(line, "\n", sep = ""); flush.console()
  invisible(line)
}

md5 <- function(p) unname(tools::md5sum(p))
assert_md5 <- function(path, expected, label) {
  if (!file.exists(path)) stop(sprintf("MD5 STOP: %s missing: %s", label, path))
  got <- md5(path)
  if (!identical(got, expected)) stop(sprintf("MD5 STOP: %s md5 %s != expected %s (%s)", label, got, expected, path))
  logf("md5 OK  %-28s %s  %s (%s bytes)", label, got, path, format(file.size(path), big.mark = ","))
  invisible(got)
}
assert_governing <- function() for (k in names(GOV)) assert_md5(GOV[[k]][1], GOV[[k]][2], paste0("GOV:", k))

RES <- NULL; RES_IDS <- NULL; MD5_SEEN <- new.env()
load_resources <- function() {
  RES <<- jsonlite::fromJSON(GOV$resources[1], simplifyVector = FALSE)
  RES_IDS <<- vapply(RES, function(e) e$id, "")
  stopifnot(length(RES) == 311L)
}
res_entry <- function(id) { i <- match(id, RES_IDS); if (is.na(i)) stop("resource id not in RESOURCES.json: ", id); RES[[i]] }
# returns the path of a resource after asserting its md5 at read time; role 'seen_output' refused
rpath <- function(id) {
  e <- res_entry(id)
  if (identical(e$role, "seen_output") || identical(e$role, "full_text_after_results")) stop("refusing seen resource as input: ", id)
  if (is.null(e$md5)) stop("resource has no md5: ", id)
  if (is.null(MD5_SEEN[[id]])) {
    assert_md5(e$path, e$md5, paste0("RES:", id))
    logf("        role=%s build=%s", e$role, substr(if (is.null(e$build)) "" else e$build, 1, 120))
    if (!is.null(e$n_info)) logf("        n_info=%s", substr(e$n_info, 1, 200))
    MD5_SEEN[[id]] <- e$md5
  }
  e$path
}
# an input not in the resource map (path + md5 asserted here; listed in README with the reason)
xpath <- function(path, expected, label) {
  if (is.null(MD5_SEEN[[label]])) { assert_md5(path, expected, paste0("EXTRA:", label)); MD5_SEEN[[label]] <- expected }
  path
}

# ---- streamed extraction: keeps only lines with a given prefix or containing a fixed token ----
# Checkpointed: a finished extract is reused only when its spec digest, its source md5 and its own md5 match.
stream_extract <- function(src_id, tag, keep_prefix = character(), keep_fixed = character(), gz = TRUE,
                           vcf = FALSE, chunk = 2000000L) {
  path <- rpath(src_id); src_md5 <- res_entry(src_id)$md5
  out_file <- file.path(OUT, "extracts", paste0(tag, ".tsv")); meta_file <- paste0(out_file, ".meta.json")
  spec <- list(src_id = src_id, src_md5 = src_md5, keep_prefix = sort(keep_prefix), keep_fixed = sort(keep_fixed), vcf = vcf)
  spec_digest <- digest(spec)
  if (file.exists(meta_file) && file.exists(out_file)) {
    m <- jsonlite::fromJSON(meta_file)
    if (identical(m$spec_digest, spec_digest) && isTRUE(m$complete) && identical(md5(out_file), m$out_md5)) {
      logf("extract %s REUSED from checkpoint (spec %s, %s data rows streamed, %d kept, out md5 %s)",
           tag, spec_digest, format(m$n_data_rows, big.mark = ","), m$n_kept, m$out_md5)
      return(m)
    }
  }
  t1 <- Sys.time()
  logf("extract %s: streaming %s (%s) for %d prefixes + %d fixed tokens", tag, src_id, path, length(keep_prefix), length(keep_fixed))
  con <- if (gz) gzfile(path, open = "rt") else file(path, open = "rt")
  header <- NULL; n <- 0; kept <- character(); last <- NA_character_; nchunk <- 0L; n_meta <- 0L
  repeat {
    x <- readLines(con, n = chunk, warn = FALSE)
    if (!length(x)) break
    nchunk <- nchunk + 1L
    if (is.null(header)) {
      if (vcf) { mh <- startsWith(x, "##"); n_meta <- sum(mh); x <- x[!mh] }
      header <- x[1]; x <- x[-1]
    }
    n <- n + length(x)
    k <- logical(length(x))
    for (p in keep_prefix) k <- k | startsWith(x, p)
    for (f in keep_fixed) k <- k | grepl(f, x, fixed = TRUE)
    if (any(k)) kept <- c(kept, x[k])
    if (length(x)) last <- x[length(x)]
    if (nchunk %% 5L == 0L) logf("  %s: %s lines read, %d kept (%.0f s)", tag, format(n, big.mark = ","), length(kept), as.numeric(difftime(Sys.time(), t1, units = "secs")))
  }
  close(con)
  writeLines(c(header, kept), out_file)
  m <- list(src_id = src_id, src_path = path, src_md5 = src_md5, spec_digest = spec_digest, header = header,
            n_meta_lines = n_meta, n_data_rows = n, n_kept = length(kept), last_line = last, complete = TRUE,
            out_md5 = md5(out_file), seconds = round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1),
            keep_prefix = keep_prefix, keep_fixed = keep_fixed, written = format(Sys.time()))
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, pretty = TRUE), meta_file)
  logf("extract %s DONE: %s data rows, %d kept, %.0f s, last line starts '%s'", tag, format(n, big.mark = ","), length(kept), m$seconds, substr(last, 1, 40))
  m
}
read_extract <- function(tag) {
  f <- file.path(OUT, "extracts", paste0(tag, ".tsv"))
  d <- fread(file = f, sep = "\t", colClasses = "character", header = TRUE, fill = TRUE)
  n_lines <- length(readLines(f)) - 1L
  if (nrow(d) != n_lines) stop(sprintf("read_extract %s: fread %d rows vs %d lines", tag, nrow(d), n_lines))
  d
}

# ---- staging snapshot (proves nothing was written under V4_local/staging/l3-finemap/) ----
staging_snapshot <- function(tag) {
  fs <- sort(list.files(STAGING, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE))
  dt <- data.table(file = substring(fs, nchar(STAGING) + 2L), bytes = file.size(fs), md5 = unname(tools::md5sum(fs)))
  f <- file.path(OUT, sprintf("staging_snapshot_%s.tsv", tag)); fwrite(dt, f, sep = "\t")
  dg <- digest(paste(dt$file, dt$md5, collapse = "\n"), algo = "md5", serialize = FALSE)
  logf("staging snapshot [%s]: %d files, %s bytes, digest-of-(file,md5) %s -> %s", tag, nrow(dt), format(sum(dt$bytes), big.mark = ","), dg, basename(f))
  list(dt = dt, digest = dg)
}

# ---- PLINK ----
qa <- function(x) shQuote(x, type = "cmd")
run_plink <- function(plink, args, tag) {
  out <- suppressWarnings(system2(plink, args, stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); if (is.null(st)) st <- 0L
  writeLines(out, file.path(OUT, "plink", paste0(tag, ".stdout.txt")))
  logf("plink %s: exit %d; args: %s", tag, st, paste(args, collapse = " "))
  invisible(list(status = st, out = out))
}

# ---- targets: lead and proxy variants with GRCh37 (503 panel bim) and GRCh38 (liftover) positions ----
PAIRS <- data.table(
  pair         = c("RAD23B_T2D", "ADM_AF", "PECAM1_CAD", "ABO_Stroke", "ABO_T2D"),
  gene         = c("RAD23B", "ADM", "PECAM1", "ABO", "ABO"),
  outcome_name = c("T2D", "AF", "CAD", "Stroke", "T2D"),
  outcome_id   = c("FinnGen_R10_T2D", "FinnGen_R10_I9_AF", "FinnGen_R10_I9_CHD", "MEGASTROKE_EUR_GCST006906", "FinnGen_R10_T2D"),
  out_res      = c("OUT_FG_R10_T2D", "OUT_FG_R10_I9_AF", "OUT_FG_R10_I9_CHD", "OUT_MEGASTROKE_EUR_GCST006906", "OUT_FG_R10_T2D"),
  out_tag      = c("FG_T2D", "FG_AF", "FG_CHD", "MEGA", "FG_T2D"),
  exp_res      = c("EXP_UKBPPP_EUR_RAD23B_OID21166", "EXP_UKBPPP_EUR_ADM_OID21467", NA, "EXP_UKBPPP_EUR_ABO_OID30675", "EXP_UKBPPP_EUR_ABO_OID30675"))
# the external-check instruments (outcome-check pre-specification section 1), for the SE table
EXT <- data.table(
  pair = c("PAM_T2D", "RAD23B_T2D", "ARG1_T2D", "PCSK9_CAD", "TGFB1_CAD"),
  gene = c("PAM", "RAD23B", "ARG1", "PCSK9", "TGFB1"), outcome_name = c("T2D", "T2D", "T2D", "CAD", "CAD"),
  res  = c("DIAM", "DIAM", "DIAM", "CARDIO", "CARDIO"))

liftover_37to38 <- function(chr, pos) {
  suppressPackageStartupMessages({ library(rtracklayer); library(GenomicRanges) })
  gz <- rpath("REF_LIFTOVER_hg19ToHg38")
  ch <- file.path(OUT, "ref", "hg19ToHg38.over.chain")
  if (!file.exists(ch)) R.utils::gunzip(gz, destname = ch, remove = FALSE, overwrite = TRUE)
  chain <- rtracklayer::import.chain(ch)
  gr <- GenomicRanges::GRanges(paste0("chr", chr), IRanges::IRanges(as.integer(pos), width = 1L))
  lo <- rtracklayer::liftOver(gr, chain)
  vapply(seq_along(lo), function(i) { s <- start(lo[[i]]); if (length(s) == 1L) as.integer(s) else NA_integer_ }, 1L)
}

build_targets <- function() {
  f <- file.path(OUT, "extracts", "targets.tsv")
  ins <- fread(file = rpath("TAB_instruments")); stopifnot(nrow(ins) == 1964L)
  fig4 <- fread(file = rpath("FIG04_panel1_2")); stopifnot(nrow(fig4) == 81L)
  rows <- list()
  for (g in unique(PAIRS$gene)) {
    ii <- ins[gene == g]; stopifnot(nrow(ii) == 1L)
    pr <- fig4[gene == g & proxy_used == TRUE]
    m <- regmatches(pr$proxy_snps, regexec("^(rs[0-9]+)<-(rs[0-9]+)\\(r2=([0-9.]+)\\)$", pr$proxy_snps))
    leads <- unique(vapply(m, `[`, "", 2)); proxies <- unique(vapply(m, `[`, "", 3))
    stopifnot(length(leads) == 1L, length(proxies) == 1L, leads == ii$SNP)
    b37 <- as.integer(strsplit(ii$ukb_variant_id, ":")[[1]][2])
    rows[[length(rows) + 1]] <- data.table(gene = g, role = "lead", rsid = ii$SNP, bim_id = ii$rsid_ld_panel,
      chr = as.character(ii$gene_chrom), pos38_src = as.integer(strsplit(ii$chrpos, ":")[[1]][2]), pos37_ukb = b37)
    rows[[length(rows) + 1]] <- data.table(gene = g, role = "proxy", rsid = proxies, bim_id = proxies,
      chr = as.character(ii$gene_chrom), pos38_src = NA_integer_, pos37_ukb = NA_integer_)
  }
  for (g in setdiff(EXT$gene, PAIRS$gene)) {
    ii <- ins[gene == g]; stopifnot(nrow(ii) == 1L)
    rows[[length(rows) + 1]] <- data.table(gene = g, role = "ext_instrument", rsid = ii$SNP, bim_id = ii$rsid_ld_panel,
      chr = as.character(ii$gene_chrom), pos38_src = as.integer(strsplit(ii$chrpos, ":")[[1]][2]),
      pos37_ukb = as.integer(strsplit(ii$ukb_variant_id, ":")[[1]][2]))
  }
  tg <- rbindlist(rows)
  # GRCh37 positions and alleles from the 503 panel, through PLINK (--make-just-bim on the extracted ids)
  plink <- rpath("PLINK19_paper6"); bed <- rpath("LD503_bed"); invisible(rpath("LD503_bim")); invisible(rpath("LD503_fam"))
  bfile <- sub("\\.bed$", "", bed)
  ids <- unique(tg$bim_id); idf <- file.path(OUT, "plink", "targets_ids.txt"); writeLines(ids, idf)
  r <- run_plink(plink, c("--bfile", qa(bfile), "--keep-allele-order", "--extract", qa(idf), "--make-just-bim", "--out", qa(file.path(OUT, "plink", "targets"))), "targets_make_just_bim")
  if (r$status != 0) stop("PLINK --make-just-bim failed")
  bim <- fread(file = file.path(OUT, "plink", "targets.bim"), header = FALSE, col.names = c("chr", "bim_id", "cm", "pos37", "A1", "A2"), colClasses = list(character = 1L))
  if (!setequal(bim$bim_id, ids)) stop(sprintf("targets: PLINK returned %d of %d ids", nrow(bim), length(ids)))
  tg <- merge(tg, bim[, .(bim_id, chr_bim = chr, pos37, A1, A2)], by = "bim_id", all.x = TRUE, sort = FALSE)
  stopifnot(all(tg$chr == tg$chr_bim))
  # leads: the panel GRCh37 position must equal the UKB-PPP variant id's GRCh37 position
  bad <- tg[!is.na(pos37_ukb) & pos37_ukb != pos37]; if (nrow(bad)) stop("targets: bim pos != UKB-PPP GRCh37 pos for ", paste(bad$rsid, collapse = ","))
  tg[, pos38 := liftover_37to38(chr, pos37)]
  tg[, pos38_method := ifelse(is.na(pos38), "UNMAPPED by hg19ToHg38 chain (position in a chain gap); resolved later by rsID + alleles in the GRCh38 outcome file", "liftover hg19ToHg38")]
  bad <- tg[!is.na(pos38_src) & (is.na(pos38) | pos38_src != pos38)]; if (nrow(bad)) stop("targets: liftover != instrument-table GRCh38 pos for ", paste(bad$rsid, collapse = ","))
  if (anyNA(tg[role != "proxy"]$pos38)) stop("targets: a lead or external instrument failed liftover")
  for (i in which(is.na(tg$pos38))) logf("targets: %s (%s) GRCh37 %s:%d does not lift (chain gap); GRCh38 position will be taken from the outcome file by rsID + alleles", tg$rsid[i], tg$gene[i], tg$chr[i], tg$pos37[i])
  fwrite(tg, f, sep = "\t")
  logf("targets: %d variants; GRCh37 from 503 bim via PLINK, GRCh38 by liftover (leads checked against the instrument table: all equal)", nrow(tg))
  tg
}

extraction_jobs <- function(tg) {
  p38 <- function(g, ro) tg[gene == g & role == ro]$pos38
  p37 <- function(g, ro) tg[gene == g & role == ro]$pos37
  list(
    FG_AF  = list(src = "OUT_FG_R10_I9_AF", gz = TRUE,
                  prefix = sprintf("11\t%d\t", na.omit(c(p38("ADM", "lead"), p38("ADM", "proxy")))),
                  fixed = c("rs2403310", "rs113773056", "rs2015148")),
    FG_CHD = list(src = "OUT_FG_R10_I9_CHD", gz = TRUE,
                  prefix = sprintf("17\t%d\t", na.omit(c(p38("PECAM1", "lead"), p38("PECAM1", "proxy")))),
                  fixed = c("rs1962147", "rs8070460")),
    MEGA   = list(src = "OUT_MEGASTROKE_EUR_GCST006906", gz = TRUE, prefix = character(),
                  fixed = c(sprintf("\t%d\t", na.omit(c(p38("ABO", "lead"), p38("ABO", "proxy")))), "rs505922", "rs576123")),
    DIAM   = list(src = "OUT_DIAMANTE2018_noUKB_EUR_chrpos", gz = FALSE,
                  prefix = sprintf("%s:%d\t", c("9", "9", "5", "6"), c(p37("RAD23B", "lead"), p37("RAD23B", "proxy"), p37("PAM", "ext_instrument"), p37("ARG1", "ext_instrument"))),
                  fixed = character()),
    DIAMRS = list(src = "OUT_DIAMANTE2018_noUKB_EUR_rsid", gz = FALSE,
                  prefix = sprintf("%s\t", c("rs7035725", "rs7856112", "rs149802978", "rs2781668")),
                  fixed = sprintf("\t%s:%d\t", c("9", "9", "5", "6"), c(p37("RAD23B", "lead"), p37("RAD23B", "proxy"), p37("PAM", "ext_instrument"), p37("ARG1", "ext_instrument")))),
    CARDIO = list(src = "OUT_CARDIOGRAMplusC4D_ieu-a-7", gz = TRUE, vcf = TRUE,
                  prefix = c(sprintf("1\t%d\t", p37("PCSK9", "ext_instrument")), sprintf("19\t%d\t", p37("TGFB1", "ext_instrument"))),
                  fixed = c("rs11591147", "rs73045269")),
    BIM503 = list(src = "LD503_bim", gz = FALSE, prefix = character(), lineno = TRUE,
                  fixed = sprintf("\t%s\t", unique(tg[role != "ext_instrument"]$bim_id)))
  )
}

run_job <- function(name, jobs) {
  j <- jobs[[name]]
  if (isTRUE(j$lineno)) return(stream_bim_lines(j$src, name, j$fixed))
  stream_extract(j$src, name, keep_prefix = j$prefix, keep_fixed = j$fixed, gz = j$gz, vcf = isTRUE(j$vcf))
}

# the 503 bim streamed with line numbers (variant index for the in-place .bed read in w1_main.R)
stream_bim_lines <- function(src_id, tag, fixed, chunk = 2000000L) {
  path <- rpath(src_id); src_md5 <- res_entry(src_id)$md5
  out_file <- file.path(OUT, "extracts", paste0(tag, ".tsv")); meta_file <- paste0(out_file, ".meta.json")
  spec_digest <- digest(list(src_id = src_id, src_md5 = src_md5, fixed = sort(fixed), lineno = TRUE))
  if (file.exists(meta_file) && file.exists(out_file)) {
    m <- jsonlite::fromJSON(meta_file)
    if (identical(m$spec_digest, spec_digest) && isTRUE(m$complete) && identical(md5(out_file), m$out_md5)) {
      logf("extract %s REUSED from checkpoint (%s lines)", tag, format(m$n_data_rows, big.mark = ",")); return(m) }
  }
  t1 <- Sys.time(); con <- file(path, open = "rt"); n <- 0; kept <- character()
  repeat {
    x <- readLines(con, n = chunk, warn = FALSE); if (!length(x)) break
    k <- logical(length(x)); for (f in fixed) k <- k | grepl(f, x, fixed = TRUE)
    if (any(k)) kept <- c(kept, paste0(n + which(k), "\t", x[k]))
    n <- n + length(x)
  }
  close(con)
  writeLines(c("lineno\tchr\tbim_id\tcm\tpos37\tA1\tA2", kept), out_file)
  m <- list(src_id = src_id, src_md5 = src_md5, spec_digest = spec_digest, n_data_rows = n, n_kept = length(kept),
            complete = TRUE, out_md5 = md5(out_file), seconds = round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1))
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, pretty = TRUE), meta_file)
  logf("extract %s DONE: %s bim lines, %d kept, %.0f s", tag, format(n, big.mark = ","), length(kept), m$seconds)
  m
}
