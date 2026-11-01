#!/usr/bin/env Rscript
## C1 -- EXTERNAL CHECKS FROM RAW (Paper 3.1 R2, Stage 1b, "middle path").
## Governing: STAGE1B_BRIEFS_v3.md "C1"; STAGE1_BRIEFS.md v3 "Governing documents" + section 0;
##            outcome-check pre-specification (OCP) sections 1, 2, 4, 5; pre-declaration (redacted) Principle 0 and section 1.
## Rerun everything:  Rscript "<this folder>/C1_extchecks.R"  then  Rscript "<this folder>/C1_compare_seen.R"
## Resumable: each resource scan writes extract_<tag>.tsv + extract_<tag>.tsv.done; a rerun skips a scan whose
## .done records the same source md5.  Progress lines are flushed to run.log every 2M rows.
## This script never opens a seen_output except S1R1_overlap_ldsc.tsv, read for the r0 label only (READ_FROM_SEEN),
## as the C1 brief directs.  COMPARE_SEEN is written by the second script, after results.tsv exists.

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
suppressPackageStartupMessages({
  library(data.table); library(jsonlite)
  library(GenomicRanges); library(rtracklayer)
})
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_11_extchecks")
RES <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
V3  <- .P31R_V3
LOG <- file.path(OUT, "run.log")
T0  <- Sys.time()
logf <- function(...) {
  s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(s, "\n", sep = "", file = LOG, append = TRUE); cat(s, "\n", sep = ""); flush.console(); invisible(s)
}
md5 <- function(p) unname(tools::md5sum(p))
`%||%` <- function(a, b) if (is.null(a)) b else a
fmt <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.logical(x)) return(ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE")))
  if (is.numeric(x)) return(ifelse(is.na(x), "NA", trimws(formatC(x, digits = 12, format = "g"))))   # formatC pads to width; trimmed
  as.character(x)
}

## ---------------------------------------------------------------- assertion ledger (asserts never stop the run)
ASSERTS <- list()
chk <- function(id, cond, msg) {
  ok <- isTRUE(cond)
  ASSERTS[[length(ASSERTS) + 1]] <<- data.table(assert_id = id, pass = ok, msg = msg)
  logf(sprintf("ASSERT %-46s %s  %s", id, if (ok) "PASS  " else "FAILED", msg))
  ok
}

logf("==================== C1 RUN START ====================")
logf("tools: ", R.version.string, "; data.table ", as.character(packageVersion("data.table")),
     "; jsonlite ", as.character(packageVersion("jsonlite")), "; rtracklayer ", as.character(packageVersion("rtracklayer")),
     "; GenomicRanges ", as.character(packageVersion("GenomicRanges")), "; R.utils ", as.character(packageVersion("R.utils")))
logf("seeds: none used (C1 has no bootstrap or other random step)")
logf("PLINK: not used (C1 computes no LD); gls_q.R: not sourced (C1 runs no GLS Q) -- section 0 rules 3, 4, 6 not applicable")

## ---------------------------------------------------------------- 1. governing md5s (STOP on mismatch, section 0 rule 1)
## Re-pinned 2026-10-08 (resume) to the md5s frozen in the resuming workflow prompt. STAGE1B_VERIFIER_BRIEF 09b81ec2... supersedes
## 71dc91af... and 4a50dcb3... (additive; tolerance-class table unchanged); RESOURCES.json is now 546b769a... (was 968b814e...).
GOV <- rbindlist(list(
  list("STAGE1B_BRIEFS_v3",      file.path(RES, "STAGE1B_BRIEFS_v3.md"), "e65e601cd3f466a02039ef05b964fef9", "workflow prompt"),
  list("STAGE1B_VERIFIER_BRIEF", file.path(RES, "STAGE1B_VERIFIER_BRIEF.md"), "09b81ec22e9cd737fc42e3139f229b35", "workflow prompt (resume; supersedes 71dc91af.../4a50dcb3...)"),
  list("STAGE1_BRIEFS_v3",       file.path(RES, "STAGE1_BRIEFS.md"), "40d0d52f7055542d774c03f1edf04e00", "workflow prompt"),
  list("PRESPEC_REDACTED",       file.path(RES, "P31R_R2_STAGE1_PRESPEC_REDACTED.md"), "a281c3464619b1b67c28e13e7392c41c", "workflow prompt"),
  list("AMEND1_REDACTED",        file.path(RES, "P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), "f265d600973c43239bb8377e0cd39d33", "workflow prompt"),
  list("OUTCOME_CHECK_PRESPEC",  file.path(V3, "R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), "c9079f2de1ffe75ce31f8a1d597d1fda", "workflow prompt"),
  list("RESOURCES_JSON",         file.path(RES, "RESOURCES.json"), "546b769ae93a8f566eec2b92f4be1575", "workflow prompt (resume)"),
  list("PRESPEC_FULL",           file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), "3319dc39740d9f5b3b8822133d9108b7",
       "workflow prompt (md5 only here; opened only by C1_compare_seen.R after results.tsv)"),
  list("AMEND1_FULL",            file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), "377508832ea9bd384133e869c456c909",
       "workflow prompt (md5 only here; opened only by C1_compare_seen.R after results.tsv)"),
  list("GLS_Q_R",                file.path(RES, "gls_q.R"), "b7390fdf5b6f2f9c1d30c298d6afdafc",
       "workflow prompt (md5 only; not sourced: C1 runs no GLS Q)")
))
setnames(GOV, c("id", "path", "md5_expected", "source"))
GOV[, md5_actual := vapply(path, md5, "")]
gov_mismatch <- function(G) G$id[G$md5_actual != G$md5_expected]   # the gate; M6 mutates its input
for (i in seq_len(nrow(GOV))) logf(sprintf("governing %-24s %s  expected %s  actual %s  [%s]", GOV$id[i],
                                           if (GOV$md5_actual[i] == GOV$md5_expected[i]) "OK      " else "MISMATCH",
                                           GOV$md5_expected[i], GOV$md5_actual[i], GOV$source[i]))
if (length(gov_mismatch(GOV))) { logf("STOP: governing md5 mismatch: ", paste(gov_mismatch(GOV), collapse = ", ")); quit(status = 2) }

## ---------------------------------------------------------------- 2. resource map; inputs by id with md5 at read time
RJ <- fromJSON(file.path(RES, "RESOURCES.json"), simplifyVector = FALSE)
logf("RESOURCES.json entries: ", length(RJ), " (Stage 1 brief says 311; logged as a difference, not used)")
rget <- function(id) { e <- Filter(function(z) identical(z$id, id), RJ); if (length(e) != 1) stop("resource id not unique/found: ", id); e[[1]] }
INPUTS <- list()
inp <- function(id) {
  e <- rget(id)
  if (!(e$role %in% c("raw_input", "definition", "tool"))) stop("role not usable as input before COMPARE_SEEN: ", id, " = ", e$role)
  m <- md5(e$path)
  logf(sprintf("input %-36s role=%-10s md5 %s expected %s  %s", id, e$role, m, e$md5, e$path))
  if (m != e$md5) { logf("STOP: input md5 mismatch for ", id); quit(status = 3) }
  INPUTS[[id]] <<- list(path = e$path, md5 = m, role = e$role, bytes = file.size(e$path))
  e$path
}

## ---------------------------------------------------------------- 3. staging snapshot (start)
STG <- file.path(V3, "V4_local/staging/l3-finemap")
snap <- function() {
  f <- sort(list.files(STG, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE))
  data.table(file = sub(paste0(STG, "/"), "", f, fixed = TRUE), bytes = file.size(f), md5 = unname(tools::md5sum(f)))
}
SNAP0 <- snap()
agg_md5 <- function(sn) { tf <- tempfile(tmpdir = OUT, fileext = ".tmp"); writeLines(paste(sn$file, sn$md5), tf); m <- md5(tf); file.remove(tf); m }
logf("staging snapshot START: ", nrow(SNAP0), " files under ", STG, "; aggregate md5 of the (file md5) list = ", agg_md5(SNAP0))
for (i in seq_len(nrow(SNAP0))) cat(sprintf("  STAGING_START	%s	%s	%s
", SNAP0$md5[i], SNAP0$bytes[i], SNAP0$file[i]), file = LOG, append = TRUE)

## ---------------------------------------------------------------- 4. definitions from the OCP section 1 table (asserts only)
## Transcribed from the governing OCP (md5 c9079f2d...). Used ONLY to assert that what this script rebuilds from raw
## matches the declared instruments / discovery signs. No computed value is taken from here.
PS <- data.table(
  pair   = c("PAM_T2D", "RAD23B_T2D", "ARG1_T2D", "PCSK9_CAD", "TGFB1_CAD", "FGFR4_T2D"),
  gene   = c("PAM", "RAD23B", "ARG1", "PCSK9", "TGFB1", "FGFR4"),
  chr    = c("5", "9", "6", "1", "19", "5"),
  rsid   = c("rs149802978", "rs7035725", "rs2781668", "rs11591147", "rs73045269", "rs451643"),
  ea     = c("G", "G", "T", "T", "T", "T"),
  oa     = c("C", "C", "C", "G", "C", "G"),
  pos38  = c(103011162L, 107287633L, 131576138L, 55039974L, 41319286L, 177084431L),
  pos37  = c(102346866L, 110049914L, 131897278L, 55505647L, 41825191L, NA_integer_),
  disc_wald_printed = c(-0.0900, -0.4804, -0.1910, 0.1817, 0.1643, 0.0482),
  ## OCP section 1 exposure columns as printed, with the half-unit of the last printed digit as the assert tolerance
  bexp_printed = c(-1.0994, -0.0641, -0.2372, -1.1533, 0.2880, 0.98),  bexp_tol = c(rep(5e-5, 5), 5e-3),
  sexp_printed = c(0.0165, 0.0076, 0.0099, 0.0278, 0.0099, 0.023),     sexp_tol = c(rep(5e-5, 5), 5e-4),
  eaf_printed  = c(0.054, 0.550, 0.167, 0.018, 0.171, 0.42),           eaf_tol  = c(rep(5e-4, 5), 5e-3),
  disc_outcome = c("FinnGen_R10_T2D", "FinnGen_R10_T2D", "FinnGen_R10_T2D", "FinnGen_R10_I9_CHD", "FinnGen_R10_I9_CHD", "AGEN_T2D_Spracklen2020"),
  resource = c(rep("OUT_DIAMANTE2018_noUKB_EUR_rsid", 3), rep("OUT_CARDIOGRAMplusC4D_ieu-a-7", 2), "OUT_TPMI_T2D_250_2"),
  resource_label = c(rep("DIAMANTE 2018 EUR T2D without UKB (Mahajan 2018)", 3), rep("CARDIoGRAMplusC4D 2015 1000G (ieu-a-7)", 2),
                     "TPMI PheWeb SAIGE PheCode 250.2"),
  ncase  = c(rep(55005L, 3), rep(60801L, 2), 59289L),   # OCP section 2 (governing); DIAMANTE has no N column
  nctrl  = c(rep(400308L, 3), rep(123504L, 2), 269125L)
)
PS[, n_eff := 4 * as.numeric(ncase) * nctrl / (as.numeric(ncase) + nctrl)]
for (i in seq_len(nrow(PS))) {
  e <- rget(PS$resource[i])
  logf(sprintf("N_eff %s: %d/%d -> %.0f (computed); RESOURCES n_info: %s", PS$pair[i], PS$ncase[i], PS$nctrl[i], PS$n_eff[i], e$n_info %||% "(none)"))
}

## ---------------------------------------------------------------- pure helpers (also exercised by the mutation tests)
NUMRE <- "^[-+]?([0-9]*\\.?[0-9]+|[0-9]+\\.)([eE][-+]?[0-9]+)?$"
is_num_token <- function(v) !is.na(v) & grepl(NUMRE, v)
num <- function(v, what) {
  ok <- is_num_token(v)
  chk(paste0("numeric_tokens:", what), all(ok), sprintf("%d/%d parse as numbers; non-numeric=[%s]", sum(ok), length(v), paste(unique(v[!ok]), collapse = ",")))
  suppressWarnings(as.numeric(ifelse(ok, v, NA)))
}
allele_set_ok <- function(a1, a2, ea, oa) !is.na(a1) && !is.na(a2) && a1 != a2 && setequal(c(a1, a2), c(ea, oa))
align <- function(eff, oth, b, f, ea, oa) {
  if (identical(eff, ea) && identical(oth, oa)) return(list(b = b, f = f, flipped = FALSE))
  if (identical(eff, oa) && identical(oth, ea)) return(list(b = -b, f = 1 - f, flipped = TRUE))
  NULL
}
is_palin <- function(a, b) paste0(a, b) %in% c("AT", "TA", "CG", "GC")
strand_status <- function(palin, eaf_exp, eaf_out) {
  if (!palin) return("not palindromic")
  maf <- min(eaf_exp, 1 - eaf_exp)
  if (maf > 0.42) return("forward-strand assumption")
  if (abs(eaf_exp - eaf_out) > 0.10) return("possibly strand-flipped")
  "palindromic; frequencies agree (<=0.10)"
}
classify <- function(wald, disc, p, evaluable) {
  if (!isTRUE(evaluable)) return("not evaluable")
  if (anyNA(c(wald, disc, p))) stop("classify: missing input with evaluable=TRUE")
  if (sign(wald) == 0 || sign(disc) == 0) stop("classify: zero sign")
  same <- sign(wald) == sign(disc)
  if (same && p < 0.05) "concordant" else if (same) "directionally concordant, not significant" else if (p < 0.05) "discordant" else "inconclusive"
}
filt_alle <- function(dt, c1, c2, ea, oa) {
  if (!nrow(dt)) return(dt)
  dt[vapply(seq_len(nrow(dt)), function(j) allele_set_ok(dt[[c1]][j], dt[[c2]][j], ea, oa), logical(1))]
}
build_at_neutral <- function(rows, chr, p37, p38, alle) {
  ## rows: data.table with CHR, POS, Allele1, Allele2 (character). Decides the build from the neutral SNP only.
  hit <- function(p) nrow(filt_alle(rows[rows$CHR == chr & rows$POS == as.character(p)], "Allele1", "Allele2", alle[1], alle[2])) > 0
  h37 <- hit(p37); h38 <- hit(p38)
  if (h38 && !h37) "GRCh38" else if (h37 && !h38) "GRCh37" else "undetermined"
}
yield_ok <- function(dt, col, key) nrow(dt[dt[[col]] == key]) == 1L
diam_pos_ok <- function(r, e) nrow(r) == 1L && r$CHR == e$chr && as.integer(r$BP) == e$pos37 && r$SNP == paste0(e$chr, ":", e$pos37)
diam_identity_ok <- function(o, r) {
  nrow(o) == 1L && nrow(r) == 1L && o$EA == r$EA && o$NEA == r$NEA &&
    isTRUE(all.equal(as.numeric(c(o$EAF, o$Beta, o$SE, o$Pvalue)), as.numeric(c(r$EAF, r$Beta, r$SE, r$Pvalue)), tolerance = 1e-12))
}
snap_same <- function(a, b) nrow(a) == nrow(b) && all(a$file == b$file) && all(a$md5 == b$md5)
no_pad <- function(v) all(!grepl("^\\s|\\s$", v))
class_invariant_ok <- function(evaluable, p, cls) {
  if (!isTRUE(evaluable)) return(cls == "not evaluable")
  if (p < 0.05) cls %in% c("concordant", "discordant") else cls %in% c("directionally concordant, not significant", "inconclusive")
}

## ---------------------------------------------------------------- 5. exposure instruments (TAB_instruments; FGFR4 per OCP)
p_ins <- inp("TAB_instruments")
INS <- fread(p_ins, colClasses = "character", sep = ",")
chk("rows:TAB_instruments", nrow(INS) == 1964L, sprintf("rows=%d (RESOURCES: 1965 lines incl. header)", nrow(INS)))
EXPO <- list()
for (g in PS$gene[1:5]) {
  r <- INS[gene == g]
  ps <- PS[gene == g]
  chk(paste0("join_yield:TAB_instruments:", g), nrow(r) == 1L, sprintf("rows for gene %s = %d", g, nrow(r)))
  r <- r[1]
  cp <- strsplit(r$chrpos, ":", fixed = TRUE)[[1]]
  uv <- strsplit(r$ukb_variant_id, ":", fixed = TRUE)[[1]]
  pos37 <- as.integer(uv[2])
  chk(paste0("instrument_identity:", g), r$SNP == ps$rsid && r$effect_allele == ps$ea && r$other_allele == ps$oa &&
        cp[1] == ps$chr && as.integer(cp[2]) == ps$pos38 && uv[1] == ps$chr && pos37 == ps$pos37,
      sprintf("TAB_instruments %s %s %s/%s GRCh38 %s, GRCh37 (ukb_variant_id) %s vs OCP %s %s:%d (%d) %s/%s",
              g, r$SNP, r$effect_allele, r$other_allele, r$chrpos, pos37, ps$rsid, ps$chr, ps$pos38, ps$pos37, ps$ea, ps$oa))
  EXPO[[g]] <- list(rsid = r$SNP, ea = r$effect_allele, oa = r$other_allele, chr = cp[1], pos38 = as.integer(cp[2]), pos37 = pos37,
                    beta = num(r$beta, paste0("beta_exp:", g)), se = num(r$se, paste0("se_exp:", g)), eaf = num(r$eaf, paste0("eaf_exp:", g)),
                    p = num(r$pval, paste0("p_exp:", g)), n = as.integer(r$samplesize), p_note = sprintf("pval_underflow=%s; log10p=%s", r$pval_underflow, r$log10p),
                    source = "UKB-PPP cis-pQTL lead (TAB_instruments)", oid = r$oid)
}
p_ckb <- inp("CKB_FGFR4_cis_eas")
CKB <- fread(p_ckb, colClasses = "character", sep = "\t")
chk("rows:CKB_FGFR4_cis_eas", nrow(CKB) == 3057L, sprintf("rows=%d (RESOURCES: 3057 data rows)", nrow(CKB)))
logf("CKB_FGFR4_cis_eas per-key null rate: rsid empty ", sprintf("%.4f", mean(CKB$rsid == "" | is.na(CKB$rsid))),
     "; position empty ", sprintf("%.4f", mean(CKB$position == "" | is.na(CKB$position))))
r <- CKB[rsid == "rs451643"]
chk("join_yield:CKB_FGFR4_cis_eas:rs451643", nrow(r) == 1L, sprintf("rows=%d", nrow(r)))
r <- r[1]
chk("instrument_identity:FGFR4", r$chr == "5" && as.integer(r$position) == 177084431L && r$alt == "T" && r$ref == "G",
    sprintf("CKB file %s %s:%s alt %s ref %s vs OCP rs451643 5:177084431 (GRCh38) T/G", r$rsid, r$chr, r$position, r$alt, r$ref))
ni_ckb <- rget("CKB_FGFR4_cis_eas")$n_info %||% ""
n_ckb <- suppressWarnings(as.integer(gsub(",", "", sub("^.*\\bN ([0-9,]+).*$", "\\1", ni_ckb, perl = TRUE))))
chk("n_from_RESOURCES:CKB_FGFR4_cis_eas", !is.na(n_ckb) && n_ckb > 0, sprintf("N parsed from RESOURCES n_info '%s' = %s (file has no N column)", ni_ckb, n_ckb))
EXPO[["FGFR4"]] <- list(rsid = r$rsid, ea = r$alt, oa = r$ref, chr = r$chr, pos38 = as.integer(r$position), pos37 = NA_integer_,
                        beta = num(r$beta, "beta_exp:FGFR4"), se = num(r$sebeta, "se_exp:FGFR4"), eaf = num(r$af, "eaf_exp:FGFR4"),
                        p = num(r$pvalue, "p_exp:FGFR4"), n = n_ckb,
                        p_note = sprintf("CKB file pvalue string '%s'; N from RESOURCES n_info (manifest)", r$pvalue),
                        source = "CKB pQTL (CKB_FGFR4_cis_eas), per OCP section 1", oid = "CKB FGFR4")
## rebuilt exposure vs the OCP section 1 table (asserts only; a sign or allele error here would flip the class)
exp_ok <- function(e, ps) isTRUE(abs(e$beta - ps$bexp_printed) <= ps$bexp_tol + 1e-12 && abs(e$se - ps$sexp_printed) <= ps$sexp_tol + 1e-12 &&
                                   abs(e$eaf - ps$eaf_printed) <= ps$eaf_tol + 1e-12)
for (i in seq_len(nrow(PS))) {
  ps <- PS[i]; e <- EXPO[[ps$gene]]
  chk(paste0("exposure_matches_OCP:", ps$gene), exp_ok(e, ps),
      sprintf("rebuilt beta %s se %s eaf %s vs OCP printed %s (%s), EAF %s", fmt(e$beta), fmt(e$se), fmt(e$eaf), ps$bexp_printed, ps$sexp_printed, ps$eaf_printed))
}

## ---------------------------------------------------------------- 6. liftover (validated on the five EUR instruments)
dir.create(file.path(OUT, "_ref"), showWarnings = FALSE)
gunz <- function(src, dst) { if (!file.exists(dst)) R.utils::gunzip(src, destname = dst, remove = FALSE, overwrite = TRUE); dst }
CH38to19 <- import.chain(gunz(inp("REF_LIFTOVER_hg38ToHg19"), file.path(OUT, "_ref", "hg38ToHg19.over.chain")))
CH19to38 <- import.chain(gunz(inp("REF_LIFTOVER_hg19ToHg38"), file.path(OUT, "_ref", "hg19ToHg38.over.chain")))
lift1 <- function(chr, pos, chain) {
  g <- GRanges(paste0("chr", chr), IRanges(pos, pos))
  r <- unlist(liftOver(g, chain))
  if (length(r) != 1) NA_integer_ else as.integer(start(r))
}
for (g in PS$gene[1:5]) {
  e <- EXPO[[g]]
  l37 <- lift1(e$chr, e$pos38, CH38to19); l38 <- lift1(e$chr, e$pos37, CH19to38)
  chk(paste0("liftover_validation:", g), identical(l37, e$pos37) && identical(l38, e$pos38),
      sprintf("hg38->hg19 %d -> %s (UKB id %d); hg19->hg38 %d -> %s", e$pos38, l37, e$pos37, e$pos37, l38))
}
EXPO$FGFR4$pos37 <- lift1("5", EXPO$FGFR4$pos38, CH38to19)
logf("FGFR4 rs451643 GRCh38 5:", EXPO$FGFR4$pos38, " -> GRCh37 5:", EXPO$FGFR4$pos37, " (hg38ToHg19 chain)")

## ---------------------------------------------------------------- 7. discovery estimates (signs for classification)
p_mr <- inp("TAB_mr_all_results")
MR <- fread(p_mr, colClasses = "character", sep = ",")
chk("rows:TAB_mr_all_results", nrow(MR) == 14537L, sprintf("rows=%d (RESOURCES: 14538 lines incl. header)", nrow(MR)))
DISC <- list()
for (i in 1:5) {
  ps <- PS[i]
  r <- MR[gene == ps$gene & outcome_id == ps$disc_outcome & method == "Wald ratio"]
  chk(paste0("join_yield:TAB_mr_all_results:", ps$pair), nrow(r) == 1L, sprintf("rows=%d", nrow(r)))
  r <- r[1]
  b <- num(r$b, paste0("disc_b:", ps$pair)); s <- num(r$se, paste0("disc_se:", ps$pair))
  chk(paste0("discovery_matches_OCP:", ps$pair), sign(b) == sign(ps$disc_wald_printed) && abs(b - ps$disc_wald_printed) <= 6e-5,
      sprintf("TAB_mr_all_results %s %s b=%.6f vs OCP printed %.4f (n_proxy=%s)", ps$gene, ps$disc_outcome, b, ps$disc_wald_printed, r$n_proxy))
  DISC[[ps$pair]] <- list(wald = b, se = s, src = sprintf("TAB_mr_all_results %s/%s Wald ratio (n_proxy=%s; SE as stored = pipeline convention)", ps$gene, ps$disc_outcome, r$n_proxy))
}

## ---------------------------------------------------------------- 8. resource scans (streamed, read-only, resumable)
scan_resource <- function(tag, path, src_md5, gz, keep, vcf = FALSE, chunk = 1000000L) {
  ext <- file.path(OUT, sprintf("extract_%s.tsv", tag)); done <- paste0(ext, ".done")
  if (file.exists(done) && file.exists(ext)) {
    d <- readLines(done)
    if (paste0("src_md5\t", src_md5) %in% d) { logf("scan ", tag, ": RESUME -- .done present with matching src md5; skipped"); return(list(ext = ext, done = d)) }
  }
  con <- if (gz) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  t1 <- Sys.time(); n <- 0; hdr <- NULL; kept <- character(0); nmeta <- 0; nchunk <- 0
  repeat {
    x <- readLines(con, n = chunk, warn = FALSE)
    if (!length(x)) break
    nchunk <- nchunk + 1
    x <- sub("\r$", "", x)
    if (vcf) { m <- startsWith(x, "##"); nmeta <- nmeta + sum(m); x <- x[!m] }
    if (is.null(hdr) && length(x)) { hdr <- x[1]; x <- x[-1] }
    n <- n + length(x)
    if (length(x)) kept <- c(kept, x[keep(x)])
    if (nchunk %% 2 == 0) logf(sprintf("scan %s: %.0f rows read, %d kept, %.0f s", tag, n, length(kept), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  }
  secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  writeLines(c(hdr, kept), ext)
  d <- c(paste0("src_md5\t", src_md5), paste0("src_path\t", path), sprintf("rows_scanned\t%.0f", n), paste0("rows_kept\t", length(kept)),
         sprintf("meta_lines\t%.0f", nmeta), sprintf("seconds\t%.1f", secs), paste0("finished\t", format(Sys.time())))
  writeLines(d, done)
  logf(sprintf("scan %s: DONE %.0f rows, %d kept, %.1f s", tag, n, length(kept), secs))
  list(ext = ext, done = d)
}
read_extract <- function(sc) {
  x <- readLines(sc$ext)
  kept <- as.integer(sub("rows_kept\t", "", grep("^rows_kept\t", sc$done, value = TRUE)))
  dt <- if (length(x) > 1) fread(text = x, colClasses = "character", sep = "\t", header = TRUE) else
    setnames(data.table(matrix(character(0), 0, length(strsplit(x[1], "\t")[[1]]))), strsplit(x[1], "\t")[[1]])
  chk(paste0("extract_rowcount:", basename(sc$ext)), nrow(dt) == kept, sprintf("parsed %d rows; .done says %d", nrow(dt), kept))
  dt
}
starts_any <- function(pfx) function(x) { k <- logical(length(x)); for (p in pfx) k <- k | startsWith(x, p); k }
grepl_any  <- function(pat) function(x) { k <- logical(length(x)); for (p in pat) k <- k | grepl(p, x, fixed = TRUE); k }

## 8a DIAMANTE (rsID file = OCP section 2 primary; original chr:pos file = identity check)
p_dr <- inp("OUT_DIAMANTE2018_noUKB_EUR_rsid")
eur_t2d <- PS[resource == "OUT_DIAMANTE2018_noUKB_EUR_rsid"]
SC_DR <- scan_resource("DIAMANTE_rsid", p_dr, INPUTS[["OUT_DIAMANTE2018_noUKB_EUR_rsid"]]$md5, FALSE, starts_any(paste0(eur_t2d$rsid, "\t")))
p_dc <- inp("OUT_DIAMANTE2018_noUKB_EUR_chrpos")
SC_DC <- scan_resource("DIAMANTE_chrpos", p_dc, INPUTS[["OUT_DIAMANTE2018_noUKB_EUR_chrpos"]]$md5, FALSE,
                       starts_any(paste0(eur_t2d$chr, ":", eur_t2d$pos37, "\t")))
## 8b CARDIoGRAMplusC4D VCF (two instruments + the neutral SNP rs7412 for the TPMI build check)
p_cg <- inp("OUT_CARDIOGRAMplusC4D_ieu-a-7")
SC_CG <- scan_resource("CARDIoGRAM_ieu-a-7", p_cg, INPUTS[["OUT_CARDIOGRAMplusC4D_ieu-a-7"]]$md5, TRUE,
                       grepl_any(paste0("\t", c(PS[resource == "OUT_CARDIOGRAMplusC4D_ieu-a-7"]$rsid, "rs7412"), "\t")), vcf = TRUE)
CG <- read_extract(SC_CG)
n7412 <- CG[ID == "rs7412"]
chk("join_yield:ieu-a-7:rs7412(neutral)", nrow(n7412) == 1L, sprintf("rows=%d", nrow(n7412)))
R7412 <- list(chr = n7412$`#CHROM`[1], p37 = as.integer(n7412$POS[1]), a = c(n7412$REF[1], n7412$ALT[1]))
R7412$p38 <- lift1(R7412$chr, R7412$p37, CH19to38)
logf(sprintf("neutral SNP rs7412: GRCh37 %s:%d (ieu-a-7) -> GRCh38 %s:%s (hg19ToHg38); alleles %s/%s",
             R7412$chr, R7412$p37, R7412$chr, R7412$p38, R7412$a[1], R7412$a[2]))
## 8c TPMI (one pass; keeps rs7412 at both candidate builds and FGFR4 at both; build decided at rs7412 ONLY)
p_tp <- inp("OUT_TPMI_T2D_250_2")
SC_TP <- scan_resource("TPMI_250_2", p_tp, INPUTS[["OUT_TPMI_T2D_250_2"]]$md5, TRUE,
                       starts_any(c(paste0(R7412$chr, "\t", c(R7412$p37, R7412$p38), "\t"),
                                    paste0("5\t", c(EXPO$FGFR4$pos37, EXPO$FGFR4$pos38), "\t"))))
## 8d AGEN raw (FGFR4 discovery, GRCh37)
p_ag <- inp("EAS_AGEN_T2D_Spracklen2020")
SC_AG <- scan_resource("AGEN_T2D", p_ag, INPUTS[["EAS_AGEN_T2D_Spracklen2020"]]$md5, TRUE,
                       grepl_any(paste0("\t5\t", EXPO$FGFR4$pos37, "\t")))
p_aw <- inp("EAS_AGEN_T2D_FGFR4_window")

## ---------------------------------------------------------------- 9. per-pair harmonisation and estimation
DR <- read_extract(SC_DR); DC <- read_extract(SC_DC); TP <- read_extract(SC_TP); AG <- read_extract(SC_AG)
logf(sprintf("per-key null rates in extracts: DIAMANTE_rsid rsID empty %.3f; DIAMANTE_chrpos SNP empty %.3f; ieu-a-7 ID empty %.3f; TPMI MarkerID empty %.3f; AGEN Chrpos empty %.3f",
             mean(DR$rsID == ""), mean(DC$SNP == ""), mean(CG$ID == ""), if (nrow(TP)) mean(TP$MarkerID == "") else NA, if (nrow(AG)) mean(AG$Chrpos == "") else NA))

## FGFR4 discovery (AGEN raw, aligned to T, divided by the CKB beta) with the window file as a cross-check
{
  e <- EXPO$FGFR4
  a <- AG[AG$Chr == "5" & AG$Pos == as.character(e$pos37)]
  a[, EAu := toupper(EA)][, NEAu := toupper(NEA)]
  a <- filt_alle(a, "EAu", "NEAu", e$ea, e$oa)
  if (!chk("join_yield:AGEN_T2D:FGFR4", nrow(a) == 1L, sprintf("rows at 5:%d with allele set {%s,%s} = %d", e$pos37, e$ea, e$oa, nrow(a))))
    stop("FGFR4 discovery row not found in AGEN raw; the discovery sign cannot be derived")
  al <- align(a$EAu[1], a$NEAu[1], num(a$Beta[1], "AGEN_beta"), num(a$EAF[1], "AGEN_eaf"), e$ea, e$oa)
  sAG <- num(a$SE[1], "AGEN_se")
  W <- fread(p_aw, colClasses = "character", sep = "\t")
  w <- W[rsid == e$rsid]
  chk("crosscheck:AGEN_window_vs_raw:FGFR4", nrow(w) == 1L && w$alt == e$ea && abs(as.numeric(w$beta) - al$b) <= 1e-4 && abs(as.numeric(w$sebeta) - sAG) <= 1e-4,
      sprintf("window file %s alt %s beta %s se %s vs raw aligned beta %.4f se %.4f", w$rsid, w$alt, w$beta, w$sebeta, al$b, sAG))
  dw <- al$b / e$beta
  chk("discovery_matches_OCP:FGFR4_T2D", sign(dw) == 1 && abs(dw - 0.0482) <= 6e-5, sprintf("AGEN %.4f / CKB %.3f = %.6f vs OCP printed +0.0482", al$b, e$beta, dw))
  DISC[["FGFR4_T2D"]] <- list(wald = dw, se = sAG / abs(e$beta), src = "AGEN Spracklen 2020 raw (GRCh37) / CKB pQTL beta; first-order SE")
}

ROWS <- list(); BR <- list(); STATE <- list()
addr <- function(qid, ref, est, res, value, se = NA, lo = NA, hi = NA, p = NA, n = NA, unit = "", notes = "", panel = "NA", cov = "NA", se_order = "NA")
  ROWS[[length(ROWS) + 1]] <<- data.table(quantity_id = qid, prespec_ref = ref, estimand = est, resource = res, panel = panel, cov = cov, se_order = se_order,
                                          value = fmt(value), se = fmt(se), ci_lo = fmt(lo), ci_hi = fmt(hi), p = fmt(p), n = fmt(n), unit = unit, notes = notes)

for (i in seq_len(nrow(PS))) {
  ps <- PS[i]; e <- EXPO[[ps$gene]]; pr <- ps$pair; notes_out <- character(0)
  out <- NULL
  if (ps$resource == "OUT_DIAMANTE2018_noUKB_EUR_rsid") {
    r <- DR[DR$rsID == e$rsid]
    chk(paste0("join_yield:DIAMANTE_rsid:", pr), yield_ok(DR, "rsID", e$rsid), sprintf("rows with rsID %s = %d", e$rsid, nrow(r)))
    if (nrow(r) >= 1) {
      r <- filt_alle(r, "EA", "NEA", e$ea, e$oa)
      pos_ok <- chk(paste0("position_alleles:DIAMANTE_rsid:", pr), diam_pos_ok(r, e),
                    sprintf("rows with matching allele set = %d; CHR %s BP %s SNP %s vs GRCh37 %s:%d", nrow(r), r$CHR[1], r$BP[1], r$SNP[1], e$chr, e$pos37))
      if (nrow(r) == 1L) {
        al <- align(r$EA, r$NEA, num(r$Beta, paste0("Beta:", pr)), num(r$EAF, paste0("EAF:", pr)), e$ea, e$oa)
        out <- list(b = al$b, se = num(r$SE, paste0("SE:", pr)), p = num(r$Pvalue, paste0("Pvalue:", pr)), eaf = al$f, flipped = al$flipped,
                    pcol = "Pvalue (DIAMANTE rsID file)", raw = sprintf("EA %s NEA %s EAF %s Beta %s SE %s P %s", r$EA, r$NEA, r$EAF, r$Beta, r$SE, r$Pvalue),
                    n = ps$ncase + ps$nctrl, nnote = "N from OCP section 2 counts (file has no N column)")
        ## identity with the original chr:pos file
        o <- filt_alle(DC[DC$SNP == paste0(e$chr, ":", e$pos37)], "EA", "NEA", e$ea, e$oa)
        chk(paste0("identity:DIAMANTE_rsid_vs_original:", pr), diam_identity_ok(o, r),
            sprintf("original rows=%d; original EA %s NEA %s EAF %s Beta %s SE %s P %s | rsID file %s", nrow(o), o$EA[1], o$NEA[1], o$EAF[1], o$Beta[1], o$SE[1], o$Pvalue[1], out$raw))
        notes_out <- c(notes_out, "DIAMANTE Beta/SE printed to 2 s.f., so Wald and p_wald1 are rounded-Wald quantities")
      }
    }
  } else if (ps$resource == "OUT_CARDIOGRAMplusC4D_ieu-a-7") {
    r <- CG[CG$ID == e$rsid]
    chk(paste0("join_yield:ieu-a-7:", pr), yield_ok(CG, "ID", e$rsid), sprintf("rows with ID %s = %d", e$rsid, nrow(r)))
    if (nrow(r) >= 1) {
      r <- filt_alle(r, "REF", "ALT", e$ea, e$oa)
      chk(paste0("position_alleles:ieu-a-7:", pr), nrow(r) == 1L && r$`#CHROM` == e$chr && as.integer(r$POS) == e$pos37,
          sprintf("rows with matching allele set = %d; %s:%s REF %s ALT %s vs GRCh37 %s:%d", nrow(r), r$`#CHROM`[1], r$POS[1], r$REF[1], r$ALT[1], e$chr, e$pos37))
      if (nrow(r) == 1L) {
        keys <- strsplit(r$FORMAT, ":", fixed = TRUE)[[1]]; vals <- strsplit(r[["ieu-a-7"]], ":", fixed = TRUE)[[1]]
        names(vals) <- keys
        chk(paste0("vcf_format:", pr), all(c("ES", "SE", "LP", "AF", "SS") %in% keys), paste("FORMAT", r$FORMAT))
        al <- align(r$ALT, r$REF, num(vals[["ES"]], paste0("ES:", pr)), num(vals[["AF"]], paste0("AF:", pr)), e$ea, e$oa)
        out <- list(b = al$b, se = num(vals[["SE"]], paste0("SE:", pr)), p = 10^(-num(vals[["LP"]], paste0("LP:", pr))), eaf = al$f, flipped = al$flipped,
                    pcol = "10^-LP (ieu-a-7 VCF)", raw = sprintf("REF %s ALT %s %s=%s", r$REF, r$ALT, r$FORMAT, r[["ieu-a-7"]]),
                    n = as.integer(vals[["SS"]]), nnote = sprintf("VCF SS; cases/controls %d/%d from OCP section 2", ps$ncase, ps$nctrl))
      }
    }
  } else {
    bld <- build_at_neutral(TP, R7412$chr, R7412$p37, R7412$p38, R7412$a)
    chk("build_check:TPMI_at_rs7412", bld %in% c("GRCh37", "GRCh38"), sprintf("TPMI build at neutral rs7412 = %s (GRCh37 %d / GRCh38 %s)", bld, R7412$p37, R7412$p38))
    STATE$tpmi_build <- bld
    fpos <- if (bld == "GRCh38") e$pos38 else if (bld == "GRCh37") e$pos37 else NA
    r <- TP[TP$CHR == "5" & TP$POS == as.character(fpos)]
    chk(paste0("join_yield:TPMI:", pr), nrow(r) >= 1L, sprintf("rows at 5:%s (%s) = %d", fpos, bld, nrow(r)))
    r <- filt_alle(r, "Allele1", "Allele2", e$ea, e$oa)
    chk(paste0("position_alleles:TPMI:", pr), nrow(r) == 1L, sprintf("rows with allele set {%s,%s} = %d; MarkerID %s", e$ea, e$oa, nrow(r), paste(r$MarkerID, collapse = ",")))
    if (nrow(r) == 1L) {
      al <- align(r$Allele2, r$Allele1, num(r$BETA, paste0("BETA:", pr)), num(r$AF_Allele2, paste0("AF_Allele2:", pr)), e$ea, e$oa)
      chk(paste0("tpmi_N_matches_OCP:", pr), as.integer(r$N_case) == ps$ncase && as.integer(r$N_ctrl) == ps$nctrl,
          sprintf("file N_case %s N_ctrl %s vs OCP %d/%d", r$N_case, r$N_ctrl, ps$ncase, ps$nctrl))
      out <- list(b = al$b, se = num(r$SE, paste0("SE:", pr)), p = num(r$p.value, paste0("p.value:", pr)), eaf = al$f, flipped = al$flipped,
                  pcol = "p.value (TPMI SAIGE)", raw = sprintf("MarkerID %s Allele1 %s Allele2 %s AF_Allele2 %s BETA %s SE %s p.value %s Is.SPA %s imputationInfo %s QC %s",
                                                            r$MarkerID, r$Allele1, r$Allele2, r$AF_Allele2, r$BETA, r$SE, r$p.value, r$Is.SPA, r$imputationInfo, r$QC),
                  n = as.integer(r$N_case) + as.integer(r$N_ctrl), nnote = sprintf("N_case %s + N_ctrl %s from the file", r$N_case, r$N_ctrl), qc = r$QC)
      if (!identical(r$QC, "PASS")) notes_out <- c(notes_out, sprintf("TPMI QC column = %s at the instrument (OCP sets no QC rule; class not altered) [IMPLEMENTER CHOICE]", r$QC))
    }
  }
  evaluable <- !is.null(out)
  palin <- is_palin(e$ea, e$oa)
  d <- DISC[[pr]]
  if (evaluable) {
    wald <- out$b / e$beta
    se1 <- out$se / abs(e$beta)
    sed <- sqrt(out$se^2 / e$beta^2 + out$b^2 * e$se^2 / e$beta^4)
    pw1 <- 2 * pnorm(-abs(wald / se1)); pw2 <- 2 * pnorm(-abs(wald / sed))
    cls <- classify(wald, d$wald, out$p, TRUE)
    sst <- strand_status(palin, e$eaf, out$eaf)
    hz <- (d$wald - wald) / sqrt(d$se^2 + se1^2); hp <- 2 * pnorm(-abs(hz))
  } else {
    wald <- se1 <- sed <- pw1 <- pw2 <- hz <- hp <- NA_real_
    cls <- classify(NA, d$wald, NA, FALSE); sst <- if (palin) "not evaluable" else "not palindromic"
    out <- list(b = NA, se = NA, p = NA, eaf = NA, pcol = "", raw = "instrument absent or alleles do not match", n = NA, nnote = "", flipped = NA)
  }
  bonf <- if (evaluable) out$p < 0.05 / 6 else NA
  cc1 <- if (cls == "concordant") "met" else "not met"
  ## asserted invariants of the classifier on this pair
  chk(paste0("class_invariant:", pr), class_invariant_ok(evaluable, out$p, cls), sprintf("class '%s' consistent with the p side", cls))
  if (ps$gene == "RAD23B") chk("rad23b_forward_strand_label_derived", sst == "forward-strand assumption", paste("derived strand status:", sst))
  r0 <- NA_real_
  STATE[[pr]] <- list(wald = wald, se1 = se1, sed = sed, pw1 = pw1, pw2 = pw2, cls = cls, sst = sst, bonf = bonf, cc1 = cc1, out = out, d = d, hz = hz, hp = hp,
                      evaluable = evaluable, palin = palin, notes = notes_out)
}

## ---------------------------------------------------------------- 10. r0 label (READ_FROM_SEEN, label only, per the C1 brief)
p_r0 <- file.path(V3, "R1_stage1/S1R1_05b_outcome_checks/S1R1_overlap_ldsc.tsv")
m_r0 <- md5(p_r0)
logf("READ_FROM_SEEN (label only): ", p_r0, " md5 ", m_r0, " (not in RESOURCES; registered in README)")
OV <- fread(text = readLines(p_r0), sep = "\t", colClasses = "character")
chk("rows:S1R1_overlap_ldsc", nrow(OV) == 2L && "null_z_corr" %in% names(OV), sprintf("rows=%d; columns: %s", nrow(OV), paste(names(OV), collapse = ",")))
R0 <- list(DIAMANTE = as.numeric(OV[grepl("DIAMANTE", pair, fixed = TRUE)]$null_z_corr),
           CARDIoGRAM = as.numeric(OV[grepl("CARDIoGRAMplusC4D", pair, fixed = TRUE)]$null_z_corr))
chk("r0_label_values_found", length(R0$DIAMANTE) == 1 && length(R0$CARDIoGRAM) == 1, "one null_z_corr row per European check resource")

## ---------------------------------------------------------------- 11. results.tsv and branches.tsv
OCP4 <- "OCP section 4"; OCP5 <- "OCP section 5"
for (i in seq_len(nrow(PS))) {
  ps <- PS[i]; e <- EXPO[[ps$gene]]; pr <- ps$pair; s <- STATE[[pr]]; o <- s$out
  q <- function(k) paste0("C1.", pr, ".", k)
  res_out <- sprintf("%s [%s]", ps$resource_label, ps$resource)
  instr <- sprintf("%s (%s:%d GRCh38; EA %s / OA %s)", e$rsid, e$chr, e$pos38, e$ea, e$oa)
  addr(q("beta_exp"), "OCP section 1; C1 brief (TAB_instruments; FGFR4 per OCP)", paste("per-allele protein effect of EA at", instr), e$source,
       e$beta, se = e$se, p = e$p, n = e$n, unit = "SD protein per EA", notes = e$p_note)
  addr(q("se_exp"), "OCP section 1", paste("SE of beta_exp at", e$rsid), e$source, e$se, unit = "SD protein")
  addr(q("eaf_exp"), "OCP section 1 / section 4", paste("EA frequency, exposure resource, at", e$rsid), e$source, e$eaf, unit = "frequency")
  addr(q("beta_out"), "OCP section 4 (harmonisation by allele label)", paste("per-EA log-odds in the check resource at", instr, "aligned to the exposure EA"), res_out,
       o$b, se = o$se, p = o$p, n = o$n, unit = "log-odds per EA",
       notes = paste0(if (isTRUE(o$flipped)) "resource allele flipped to exposure EA; " else if (isFALSE(o$flipped)) "resource effect allele = exposure EA; " else "", o$raw, "; ", o$nnote))
  addr(q("se_out"), "OCP section 4", paste("SE of beta_out at", e$rsid), res_out, o$se, unit = "log-odds")
  addr(q("eaf_out"), "OCP section 4 (palindromic: EAF in both resources)", paste("exposure-EA frequency in the check resource at", e$rsid), res_out, o$eaf, unit = "frequency")
  addr(q("palindromic"), "OCP section 4", paste("instrument alleles", e$ea, "/", e$oa, "palindromic (A/T or C/G)"), "definition", s$palin, unit = "logical",
       notes = paste("strand status:", s$sst))
  addr(q("wald"), OCP4, "Wald ratio beta_out/beta_exp (marginal, single instrument)", res_out, s$wald, se = s$se1, lo = s$wald - qnorm(0.975) * s$se1,
       hi = s$wald + qnorm(0.975) * s$se1, p = s$pw1, n = o$n, unit = "log-odds per SD protein", se_order = "first-order",
       notes = paste(c("CI and p from the first-order SE", s$notes), collapse = "; "))
  addr(q("se1"), OCP4, "first-order Wald SE = se_out/|beta_exp|", res_out, s$se1, unit = "log-odds per SD protein", se_order = "first-order")
  addr(q("se_delta"), OCP4, "delta-method Wald SE = sqrt(se_out^2/beta_exp^2 + beta_out^2*se_exp^2/beta_exp^4)", res_out, s$sed,
       unit = "log-odds per SD protein", se_order = "delta (second-order)")
  addr(q("p_resource"), "OCP section 5; pre-declaration section 1 wording rule ('resource's own p'); C1 brief", paste("the resource's own p at", e$rsid), res_out,
       o$p, p = o$p, unit = "p", notes = paste("column:", o$pcol))
  addr(q("p_wald1"), "C1 brief (printed beside the class)", "two-sided normal p of wald/se1", res_out, s$pw1, p = s$pw1, unit = "p", se_order = "first-order",
       notes = if (ps$resource == "OUT_DIAMANTE2018_noUKB_EUR_rsid") "rounded-Wald p (DIAMANTE Beta/SE at 2 s.f.)" else "")
  addr(q("class"), OCP5, "classification on p_resource with the discovery sign", res_out, s$cls, unit = "category",
       notes = sprintf("discovery Wald %s (%s); sign(check Wald) %s", formatC(s$d$wald, digits = 6, format = "g"), s$d$src, if (is.na(s$wald)) "NA" else if (s$wald > 0) "+" else "-"))
  addr(q("bonf6"), OCP5, "p_resource < 0.05/6", res_out, s$bonf, unit = "logical",
       notes = sprintf("threshold 0.05/6 = %.6f; on p_wald1 the result would be %s (shown, not used)", 0.05 / 6, fmt(s$pw1 < 0.05 / 6)))
  addr(paste0("C1.corroboration_condition1.", pr), "pre-declaration Principle 0 ('Corroborated', condition 1)",
       "pre-declared external check concordant (resource p < 0.05 with the discovery sign)", res_out, s$cc1, unit = "category",
       notes = "conditions 2 (section 3 branch) and 3 (section 2 branch, TGFB1) are not evaluated by C1; corroboration itself is left to W9")
  ## extras (declared in the OCP, not in the C1 key list)
  addr(q("p_wald2"), "pre-declaration section 1 wording rule (second-order p beside)", "two-sided normal p of wald/se_delta", res_out, s$pw2, p = s$pw2, unit = "p",
       se_order = "delta (second-order)", notes = "EXTRA: not in the C1 key list")
  addr(q("disc_wald"), "OCP section 1 / section 4 comparison", paste("discovery Wald ratio,", ps$disc_outcome), s$d$src, s$d$wald, se = s$d$se,
       unit = "log-odds per SD protein", notes = "EXTRA: not in the C1 key list; supplies the discovery sign for the class")
  addr(q("het_z"), "OCP section 4 comparison (discovery vs check, z-test; anticonservative under overlap)", "(disc_wald - wald)/sqrt(disc_se^2 + se1^2)", res_out,
       s$hz, p = s$hp, unit = "z", se_order = "check: first-order; discovery: as stored",
       notes = "EXTRA: not in the C1 key list. [IMPLEMENTER CHOICE] discovery SE taken as stored in TAB_mr_all_results (pipeline convention); FGFR4 discovery SE first-order")
  lab <- if (ps$gene == "FGFR4") "independent East Asian check" else
    sprintf("UKB-free external consortium check; FinnGen overlap measured at r0 = %.3f", if (grepl("DIAMANTE", ps$resource)) R0$DIAMANTE else R0$CARDIoGRAM)
  addr(q("label"), "OCP section 5 label", "check label (r0 is READ_FROM_SEEN, not computed, not a key)", res_out, lab, unit = "label",
       notes = if (ps$gene == "FGFR4") "no r0 for the East Asian check" else sprintf("READ_FROM_SEEN: r0 = null_z_corr from S1R1_overlap_ldsc.tsv (md5 %s); printed at 3 dp [IMPLEMENTER CHOICE]", m_r0))
  ## branches
  quals <- c(if (ps$gene != "FGFR4") "not independent replication (Principle 0 list)",
             if (s$sst == "forward-strand assumption") "forward-strand assumption (OCP section 4)",
             if (s$sst == "possibly strand-flipped") "possibly strand-flipped (OCP section 4)",
             if (ps$gene == "FGFR4") "AGEN-TPMI participant overlap not checked (OCP section 2)",
             if (ps$gene == "RAD23B") "section 1 p-source/rounding clause and its 0.05/6 result: W1 owns the rounding intervals",
             if (ps$gene == "PCSK9" || ps$gene == "TGFB1") "resource 77% European / 13% South Asian / 6% East Asian (OCP section 2)",
             s$notes,
             sprintf("0.05/6: %s", if (isTRUE(s$bonf)) "p_resource < 0.0083" else if (isFALSE(s$bonf)) "p_resource >= 0.0083" else "NA"))
  BR[[length(BR) + 1]] <- data.table(rule_ref = paste0(OCP5, " classification [", q("class"), "]"),
                                     inputs = sprintf("sign(wald)=%s; sign(disc)=%s; p_resource=%s; evaluable=%s", if (is.na(s$wald)) "NA" else if (s$wald > 0) "+" else "-",
                                                      if (s$d$wald > 0) "+" else "-", fmt(o$p), s$evaluable),
                                     branch = s$cls, filled_template = lab, qualifiers = paste(quals, collapse = " | "))
  BR[[length(BR) + 1]] <- data.table(rule_ref = paste0("Principle 0 'Corroborated' condition 1 [C1.corroboration_condition1.", pr, "]"),
                                     inputs = paste("class =", s$cls), branch = s$cc1, filled_template = "",
                                     qualifiers = "conditions 2-3 and the corroboration verdict are W9's")
  if (s$palin) BR[[length(BR) + 1]] <- data.table(rule_ref = paste0("OCP section 4 palindromic orientation [", q("palindromic"), "]"),
                                                  inputs = sprintf("MAF_exp=%s; eaf_exp=%s; eaf_out=%s", fmt(min(e$eaf, 1 - e$eaf)), fmt(e$eaf), fmt(o$eaf)),
                                                  branch = s$sst, filled_template = "", qualifiers = "")
}
if (!is.null(STATE$tpmi_build)) logf("TPMI build (decided at rs7412 only): ", STATE$tpmi_build)

## ---------------------------------------------------------------- 12. mutation tests (each must trip the check it names)
MT <- list()
mt <- function(id, check, mutation, tripped) {
  MT[[length(MT) + 1]] <<- data.table(mutation_id = id, check_named = check, mutation = mutation, tripped = tripped, verdict = if (isTRUE(tripped)) "PASS (tripped)" else "FAIL (did not trip)")
}
{
  ## M1 join yield: mutate an instrument id
  mt("M1", "join_yield:DIAMANTE_rsid", "instrument rsID rs2781668 -> rs2781668X", yield_ok(DR, "rsID", "rs2781668") && !yield_ok(DR, "rsID", "rs2781668X"))
  ## M2 allele set: mutate the resource EA
  r <- DR[DR$rsID == "rs2781668"]; mt("M2", "position_alleles / align (alleles do not match -> not evaluable)", paste("EA", r$EA, "-> A"),
                                      !allele_set_ok("A", r$NEA, "T", "C") && is.null(align("A", r$NEA, 1, 0.5, "T", "C")) && classify(NA, -1, NA, FALSE) == "not evaluable")
  ## M3 placeholder token
  mt("M3", "numeric_tokens", "Pvalue -> 'TBD'", !all(is_num_token(c(r$Pvalue, "TBD"))))
  ## M4 build check at the neutral SNP: shift POS by +1
  if (nrow(TP)) { TPm <- copy(TP); TPm[TPm$CHR == R7412$chr, POS := as.character(as.integer(POS) + 1L)]
    mt("M4", "build_check:TPMI_at_rs7412", "rs7412 rows POS + 1", build_at_neutral(TPm, R7412$chr, R7412$p37, R7412$p38, R7412$a) == "undetermined") }
  ## M5 classifier boundaries and invariants
  mt("M5a", "classify", "p 0.049 same sign", classify(1, 1, 0.049, TRUE) == "concordant")
  mt("M5b", "classify", "p 0.051 same sign", classify(1, 1, 0.051, TRUE) == "directionally concordant, not significant")
  mt("M5c", "classify", "p 0.01 opposite sign", classify(-1, 1, 0.01, TRUE) == "discordant")
  mt("M5d", "classify", "p 0.5 opposite sign", classify(-1, 1, 0.5, TRUE) == "inconclusive")
  mt("M5e", "classify (missing p raises)", "p NA with evaluable=TRUE", inherits(try(classify(1, 1, NA, TRUE), silent = TRUE), "try-error"))
  ## M6 md5 gate: the gate function itself, on a copy of GOV with one expected md5 altered (target asserted unique)
  G6 <- copy(GOV); k6 <- which(G6$id == "STAGE1_BRIEFS_v3"); G6$md5_expected[k6] <- sub("^4", "5", G6$md5_expected[k6])
  mt("M6", "governing md5 gate (gov_mismatch)", "expected md5 of STAGE1_BRIEFS_v3 40d0... -> 50d0...",
     length(k6) == 1L && length(gov_mismatch(GOV)) == 0L && identical(gov_mismatch(G6), "STAGE1_BRIEFS_v3"))
  ## M7 strand flag for PAM (palindromic, MAF <= 0.42): replace eaf_out by 1 - eaf_out
  s <- STATE$PAM_T2D
  if (isTRUE(s$evaluable)) mt("M7", "strand_status (possibly strand-flipped)", "PAM eaf_out -> 1 - eaf_out",
                              strand_status(TRUE, EXPO$PAM$eaf, 1 - s$out$eaf) == "possibly strand-flipped" && s$sst != "possibly strand-flipped")
  ## M8 identity check: mutate the original-file Beta string
  o <- filt_alle(DC[DC$SNP == paste0(EXPO$ARG1$chr, ":", EXPO$ARG1$pos37)], "EA", "NEA", "T", "C"); o2 <- copy(o); o2$Beta <- paste0(o2$Beta, "1")
  mt("M8", "identity:DIAMANTE_rsid_vs_original", "original-file Beta string -> Beta followed by '1'", diam_identity_ok(o, r) && !diam_identity_ok(o2, r))
  ## M9 liftover validation: mutate pos38 by +1e6
  mt("M9", "liftover_validation", "ARG1 pos38 + 1e6", !identical(lift1("6", EXPO$ARG1$pos38 + 1000000L, CH38to19), EXPO$ARG1$pos37))
  ## M10 position check: mutate BP
  r2 <- copy(r); r2$BP <- as.character(as.integer(r2$BP) + 1L)
  mt("M10", "position_alleles:DIAMANTE_rsid", "BP + 1", diam_pos_ok(r, EXPO$ARG1) && !diam_pos_ok(r2, EXPO$ARG1))
  ## M11 class invariant: a significant result labelled 'inconclusive' must be caught
  mt("M11", "class_invariant", "label 'inconclusive' with p 0.01; 'concordant' with p 0.2; 'concordant' when not evaluable",
     !class_invariant_ok(TRUE, 0.01, "inconclusive") && !class_invariant_ok(TRUE, 0.2, "concordant") && !class_invariant_ok(FALSE, NA, "concordant") &&
       class_invariant_ok(TRUE, 0.01, "concordant"))
  ## M12 staging snapshot comparison: one file's md5 altered in a copy of the start snapshot
  S12 <- copy(SNAP0); if (nrow(S12)) S12$md5[1] <- paste0(substr(S12$md5[1], 1, 31), if (substr(S12$md5[1], 32, 32) == "0") "1" else "0")
  mt("M12", "staging_unchanged (snap_same)", "first staging file md5, last hex digit changed", nrow(S12) > 0 && snap_same(SNAP0, SNAP0) && !snap_same(SNAP0, S12))
  ## M13 exposure vs OCP section 1: flip the sign of the rebuilt ARG1 beta_exp (an allele/sign error in the instrument table)
  e13 <- EXPO$ARG1; e13$beta <- -e13$beta
  mt("M13", "exposure_matches_OCP (exp_ok)", "ARG1 beta_exp sign flipped", exp_ok(EXPO$ARG1, PS[gene == "ARG1"]) && !exp_ok(e13, PS[gene == "ARG1"]))
  ## M14 padded cells: the untrimmed formatC string (the pre-fix fmt) must trip; the current fmt must not
  mt("M14", "no_padded_numeric_cells (no_pad)", "formatC(4.6e-15, digits = 12) without trimws",
     !no_pad(formatC(4.6e-15, digits = 12, format = "g")) && no_pad(fmt(c(4.6e-15, 0.0538925, -0.136440539213))))
}
MTD <- rbindlist(MT)
writeLines(c("# C1 mutation tests: each mutation must trip the check it names (TRUE = tripped).",
             paste(names(MTD), collapse = "\t"), apply(MTD, 1, paste, collapse = "\t")), file.path(OUT, "mutation_tests.txt"))
logf(sprintf("mutation tests: %d of %d tripped the check they name", sum(MTD$tripped), nrow(MTD)))

## ---------------------------------------------------------------- 13. staging snapshot (end; taken before the outputs so a failure reaches results.tsv)
SNAP1 <- snap()
same <- snap_same(SNAP0, SNAP1)
chk("staging_unchanged", same, sprintf("start %d files, end %d files; differing: %s", nrow(SNAP0), nrow(SNAP1),
                                       paste(setdiff(paste(SNAP1$file, SNAP1$md5), paste(SNAP0$file, SNAP0$md5)), collapse = ";")))
chk("no_padded_numeric_cells", no_pad(unlist(lapply(ROWS, function(d) unlist(d[, .(value, se, ci_lo, ci_hi, p, n)])))),
    "no results value/se/ci/p/n cell begins or ends with whitespace")

## ---------------------------------------------------------------- 14. write outputs
AS <- rbindlist(ASSERTS)
RESU <- rbindlist(ROWS)
if (any(!AS$pass)) RESU <- rbind(RESU, AS[pass == FALSE, .(quantity_id = paste0("ASSERT_FAILED.", assert_id), prespec_ref = "section 0 rule 9", estimand = "assertion",
                                                         resource = "", panel = "NA", cov = "NA", se_order = "NA", value = "ASSERT_FAILED", se = "", ci_lo = "", ci_hi = "",
                                                         p = "", n = "", unit = "", notes = msg)])
fwrite(RESU, file.path(OUT, "results.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(BR), file.path(OUT, "branches.tsv"), sep = "\t", quote = FALSE)
fwrite(AS, file.path(OUT, "asserts.tsv"), sep = "\t", quote = FALSE)
invisible(chk("results_written", file.exists(file.path(OUT, "results.tsv")) && nrow(fread(text = readLines(file.path(OUT, "results.tsv")), sep = "\t")) == nrow(RESU),
    sprintf("results.tsv rows = %d", nrow(RESU))))
AS <- rbindlist(ASSERTS); fwrite(AS, file.path(OUT, "asserts.tsv"), sep = "\t", quote = FALSE)   # rewritten so it holds every assert
logf(sprintf("asserts: %d pass, %d failed (all asserts, incl. results_written)", sum(AS$pass), sum(!AS$pass)))
logf("results.tsv md5 ", md5(file.path(OUT, "results.tsv")), "; branches.tsv md5 ", md5(file.path(OUT, "branches.tsv")))
logf(sprintf("C1 main script finished in %.1f s", as.numeric(difftime(Sys.time(), T0, units = "secs"))))
