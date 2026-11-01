# w4b_common.R -- shared helpers for Stage 1 W4b (pre-declaration section 4(iv)(v): ARG1 omnibus, credible sets, s,
# effector, gene naming). Sourced by every W4b step script. Writes only inside OUT.
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
suppressMessages({library(data.table); library(jsonlite)})
options(stringsAsFactors = FALSE, scipen = 0, digits = 10, warn = 1)

ROOT <- .P31R_CKM
V3 <- file.path(ROOT, "V3")
OUT <- Sys.getenv("W4B_OUT", file.path(V3, "R1_stage1/S1R1_R2_04b_arg1sets"))   # override only for a verifier re-run into another folder
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RES_DIR <- file.path(V3, "R1_stage1/S1R1_R2_00_resources")
LOG <- file.path(OUT, "run.log")
SEED <- 20261007L
NBOOT <- 5000L

# ---- governing documents (md5s passed in the workflow prompt) ----
GOV <- list(
  BRIEFS = list(path = file.path(RES_DIR, "STAGE1_BRIEFS.md"), md5 = "40d0d52f7055542d774c03f1edf04e00"),
  PRESPEC_REDACTED = list(path = file.path(RES_DIR, "P31R_R2_STAGE1_PRESPEC_REDACTED.md"), md5 = "a281c3464619b1b67c28e13e7392c41c"),
  AMEND1_REDACTED = list(path = file.path(RES_DIR, "P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), md5 = "f265d600973c43239bb8377e0cd39d33"),
  PRESPEC_FULL = list(path = file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), md5 = "3319dc39740d9f5b3b8822133d9108b7"),
  AMEND1_FULL = list(path = file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), md5 = "377508832ea9bd384133e869c456c909"),
  GLS_Q_R = list(path = file.path(RES_DIR, "gls_q.R"), md5 = "b7390fdf5b6f2f9c1d30c298d6afdafc"),
  RESOURCES = list(path = file.path(RES_DIR, "RESOURCES.json"), md5 = "546b769ae93a8f566eec2b92f4be1575"))

# ---- logging ----
logf <- function(...) {
  msg <- sprintf(...)
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG, append = TRUE)
  invisible(line)
}

md5 <- function(p) unname(tools::md5sum(p))
assert_md5 <- function(p, want, label = basename(p)) {
  if (!file.exists(p)) stop(sprintf("MISSING input %s (%s)", label, p))
  got <- md5(p)
  if (!identical(got, want)) stop(sprintf("MD5 MISMATCH %s: got %s want %s (%s)", label, got, want, p))
  logf("md5 OK  %-44s %s  %s", label, got, p)
  invisible(got)
}

# ---- resource map ----
.RES <- NULL
res <- function(id) {
  if (is.null(.RES)) {
    assert_md5(GOV$RESOURCES$path, GOV$RESOURCES$md5, "RESOURCES.json")
    .RES <<- fromJSON(GOV$RESOURCES$path, simplifyVector = FALSE)
  }
  hit <- Filter(function(e) identical(e$id, id), .RES)
  if (length(hit) != 1) stop("resource id not unique / not found: ", id)
  hit[[1]]
}
# path of a resource id, asserting its md5 at read time (and that it is not a seen_output)
rpath <- function(id, allow_seen = FALSE) {
  e <- res(id)
  if (!allow_seen && identical(e$role, "seen_output")) stop("resource ", id, " is seen_output; not an input before COMPARE_SEEN")
  assert_md5(e$path, e$md5, id)
  e$path
}
# inputs not in RESOURCES (registered in README with path, md5 and reason)
EXTRA <- list(
  GLIST_HG19 = list(path = file.path(ROOT, "paper 6/analysis/data/annotation/glist-hg19"), md5 = "d72643c547803fba23c7f70490fbe269"),
  LD489_chr6_fam = list(path = file.path(ROOT, "lhcmr_eas_refs/1000G_EUR_Phase3_plink/1000G.EUR.QC.6.fam"), md5 = NA_character_))
xpath <- function(id) {
  e <- EXTRA[[id]]
  if (is.na(e$md5)) { logf("extra input %s md5 %s (registered in README)  %s", id, md5(e$path), e$path); return(e$path) }
  assert_md5(e$path, e$md5, id)
  e$path
}

# ---- results / branches accumulators (each step writes its own part file; 05 assembles) ----
RCOLS <- c("quantity_id", "prespec_ref", "estimand", "resource", "panel", "cov", "se_order", "value", "se", "ci_lo", "ci_hi",
           "p", "n", "unit", "notes")
new_rows <- function() list()
add_row <- function(acc, quantity_id, prespec_ref, estimand, resource = "", panel = "", cov = "", se_order = "", value = NA,
                    se = NA, ci_lo = NA, ci_hi = NA, p = NA, n = NA, unit = "", notes = "") {
  fmt <- function(x) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) "" else if (is.numeric(x)) format(signif(x, 10), scientific = NA, trim = TRUE) else as.character(x)
  acc[[length(acc) + 1]] <- data.table(quantity_id = quantity_id, prespec_ref = prespec_ref, estimand = estimand,
                                       resource = resource, panel = as.character(panel), cov = cov, se_order = se_order,
                                       value = fmt(value), se = fmt(se), ci_lo = fmt(ci_lo), ci_hi = fmt(ci_hi), p = fmt(p),
                                       n = fmt(n), unit = unit, notes = notes)
  acc
}
write_part <- function(acc, name) {
  dt <- rbindlist(acc); setcolorder(dt, RCOLS)
  if (anyDuplicated(dt$quantity_id)) stop("duplicate quantity ids in part ", name, ": ", paste(head(dt$quantity_id[duplicated(dt$quantity_id)]), collapse = ", "))
  f <- file.path(OUT, "parts", paste0("results_", name, ".tsv"))
  dir.create(dirname(f), showWarnings = FALSE)
  fwrite(dt, f, sep = "\t", quote = FALSE)
  logf("wrote %s (%d rows)", f, nrow(dt))
  invisible(dt)
}

# ---- allele helpers ----
COMP <- c(A = "T", T = "A", C = "G", G = "C")
is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")
# sign to put an effect reported on (ea, oa) onto (a1, a2): +1 / -1 / 0 (no match). Exact matching only for palindromes.
align_sign <- function(ea, oa, a1, a2) {
  s <- integer(length(ea))
  same <- ea == a1 & oa == a2; swap <- ea == a2 & oa == a1
  pal <- is_pal(a1, a2)
  cea <- unname(COMP[ea]); coa <- unname(COMP[oa])
  csame <- !pal & !is.na(cea) & cea == a1 & coa == a2; cswap <- !pal & !is.na(cea) & cea == a2 & coa == a1
  s[same | csame] <- 1L; s[swap | cswap] <- -1L
  s
}

# ---- PLINK ----
PLINK_ID <- "PLINK19_paper6"
plink_exe <- function() rpath(PLINK_ID)
run_plink <- function(args, tag) {
  exe <- plink_exe()
  t0 <- Sys.time()
  out <- system2(exe, args, stdout = TRUE, stderr = TRUE)
  st <- attr(out, "status"); if (is.null(st)) st <- 0L
  writeLines(out, file.path(OUT, "plink", paste0(tag, ".stdout.txt")))
  logf("PLINK %s: exit %s in %.1fs; args: %s", tag, st, as.numeric(difftime(Sys.time(), t0, units = "secs")), paste(args, collapse = " "))
  if (st != 0) stop("PLINK failed: ", tag, "\n", paste(tail(out, 20), collapse = "\n"))
  out
}
# square r matrix from PLINK for the SNP ids in `ids` (bim ids), on bfile prefix `bfile`, restricted to chr/window.
# Returns r with dimnames = ids, in the order of `ids`, sign as PLINK's A1 (bim column 5) with --keep-allele-order.
# ASSERTS every requested id is returned, unique, and every cell non-NA.
plink_r_square <- function(bfile, ids, chr, from, to, tag) {
  dir.create(file.path(OUT, "plink"), showWarnings = FALSE)
  lst <- file.path(OUT, "plink", paste0(tag, ".extract.txt"))
  writeLines(ids, lst)
  pre <- file.path(OUT, "plink", tag)
  run_plink(c("--bfile", shQuote(bfile), "--chr", chr, "--from-bp", from, "--to-bp", to, "--extract", shQuote(lst),
              "--keep-allele-order", "--r", "square", "bin", "--threads", "2", "--out", shQuote(pre)), tag)
  # PLINK writes the matrix in .bim order of the extracted variants; recover that order with --write-snplist
  run_plink(c("--bfile", shQuote(bfile), "--chr", chr, "--from-bp", from, "--to-bp", to, "--extract", shQuote(lst),
              "--keep-allele-order", "--make-just-bim", "--out", shQuote(paste0(pre, "_order"))), paste0(tag, "_order"))
  bo <- fread(paste0(pre, "_order.bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "A1", "A2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
  k <- nrow(bo)
  binf <- paste0(pre, ".ld.bin")
  v <- readBin(binf, what = "double", n = k * k + 1L, size = 8L)
  if (length(v) != k * k) stop(sprintf("LD bin %s has %d values, expected %d", binf, length(v), k * k))
  R <- matrix(v, k, k, byrow = TRUE); dimnames(R) <- list(bo$id, bo$id)
  # assertions: every requested id returned exactly once, no NA cell
  miss <- setdiff(ids, bo$id)
  if (length(miss)) stop(sprintf("LD ASSERT FAILED (%s): %d requested SNP(s) not returned by PLINK: %s", tag, length(miss), paste(head(miss, 5), collapse = ",")))
  if (anyDuplicated(bo$id)) stop(sprintf("LD ASSERT FAILED (%s): duplicated ids returned", tag))
  if (k != length(ids)) stop(sprintf("LD ASSERT FAILED (%s): %d returned vs %d requested", tag, k, length(ids)))
  R <- R[ids, ids, drop = FALSE]
  if (anyNA(R) || any(!is.finite(R))) stop(sprintf("LD ASSERT FAILED (%s): %d NA/non-finite cells", tag, sum(!is.finite(R))))
  list(R = R, bim = bo[match(ids, bo$id)])
}

# ---- gene annotation (glist semantics: genes whose body contains pos, else the single nearest gene) ----
genes_at <- function(gt, pos) {
  ins <- gt$symbol[gt$start <= pos & pos <= gt$end]
  if (length(ins)) return(list(genes = sort(unique(ins)), dist = 0))
  d <- pmin(abs(pos - gt$start), abs(pos - gt$end))
  j <- order(d, gt$symbol)[1]          # ties broken by symbol, as Python min((d, name)) in the engine
  list(genes = gt$symbol[j], dist = d[j])
}
