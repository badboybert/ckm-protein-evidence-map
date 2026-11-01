# W6b common helpers: paths, logging, md5 assertions, resource lookup by id, PLINK runner, harmonisation units.
# Sourced by every W6b script. Writes only under OUT.
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
suppressMessages({library(data.table); library(jsonlite)})
setDTthreads(1)
options(stringsAsFactors = FALSE, scipen = 0, warn = 1)

ROOT   <- .P31R_CKM
OUT    <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_06b_renal_d_umod")
WORK   <- file.path(OUT, "work")
RESDIR <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_00_resources")
LOGF   <- file.path(OUT, "run.log")
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

say <- function(fmt, ...) {
  s <- if (length(list(...))) sprintf(fmt, ...) else fmt
  s <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s)
  cat(s, "\n", sep = "")
  cat(s, "\n", sep = "", file = LOGF, append = TRUE)
}
md5 <- function(p) unname(tools::md5sum(p))
.md5_cache <- new.env()
md5c <- function(p) {           # md5 cached within one R session (large files are hashed once per run)
  k <- normalizePath(p, winslash = "/", mustWork = TRUE)
  if (!exists(k, envir = .md5_cache, inherits = FALSE)) assign(k, md5(p), envir = .md5_cache)
  get(k, envir = .md5_cache)
}
assert_md5 <- function(p, expect, what = basename(p)) {
  if (!file.exists(p)) stop("MISSING input ", what, ": ", p)
  got <- md5c(p)
  if (!identical(got, expect)) stop(sprintf("MD5 MISMATCH %s: got %s expected %s (%s)", what, got, expect, p))
  say("md5 OK  %-44s %s  %s", what, got, p)
  invisible(p)
}

# ---- governing documents (md5s passed in the workflow prompt) ----
GOV <- list(
  briefs        = c(file.path(RESDIR, "STAGE1_BRIEFS.md"), "40d0d52f7055542d774c03f1edf04e00"),
  prespec_red   = c(file.path(RESDIR, "P31R_R2_STAGE1_PRESPEC_REDACTED.md"), "a281c3464619b1b67c28e13e7392c41c"),
  amend_red     = c(file.path(RESDIR, "P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), "f265d600973c43239bb8377e0cd39d33"),
  prespec_full  = c(file.path(ROOT, "V3/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), "3319dc39740d9f5b3b8822133d9108b7"),
  amend_full    = c(file.path(ROOT, "V3/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), "377508832ea9bd384133e869c456c909"),
  gls_q         = c(file.path(RESDIR, "gls_q.R"), "b7390fdf5b6f2f9c1d30c298d6afdafc"),
  resources     = c(file.path(RESDIR, "RESOURCES.json"), "546b769ae93a8f566eec2b92f4be1575"))
assert_governing <- function() for (k in names(GOV)) assert_md5(GOV[[k]][1], GOV[[k]][2], paste0("GOV:", k))

# ---- resource map lookup by id; md5 asserted at read time; seen roles refused unless allow_seen ----
RES <- fromJSON(GOV$resources[1], simplifyVector = FALSE)
RES_IDS <- vapply(RES, function(e) e$id, "")
res_entry <- function(id) {
  i <- which(RES_IDS == id)
  if (length(i) != 1) stop("resource id not unique/absent: ", id)
  RES[[i]]
}
res <- function(id, allow_seen = FALSE) {
  e <- res_entry(id)
  if (!identical(e$status, "FOUND")) stop("resource ", id, " status ", e$status)
  if (!allow_seen && e$role %in% c("seen_output", "seen_script", "full_text_after_results"))
    stop("resource ", id, " has role ", e$role, ": not an input before COMPARE_SEEN")
  assert_md5(e$path, e$md5, id)
  e$path
}

# ---- inputs the map lacks: registered here (path, bytes, md5, reason); asserted unchanged on every later read ----
NEWREG <- file.path(OUT, "new_inputs.tsv")
reg_new <- function(path, reason) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  m <- md5c(path); b <- file.size(path)
  old <- if (file.exists(NEWREG)) fread(NEWREG, colClasses = "character") else data.table(path = character(), bytes = character(), md5 = character(), reason = character())
  if (path %in% old$path) {
    if (old$md5[old$path == path] != m) stop("NEW input changed since registration: ", path)
  } else {
    old <- rbind(old, data.table(path = path, bytes = as.character(b), md5 = m, reason = reason))
    fwrite(old, NEWREG, sep = "\t")
  }
  say("md5 OK  %-44s %s  %s", paste0("NEW:", basename(path)), m, path)
  path
}

# ---- PLINK (one executable: PLINK19_paper6) ----
q <- function(x) shQuote(x, type = "cmd")
plink_exe <- function() res("PLINK19_paper6")
run_plink <- function(args, out_prefix, tries = 3L, timeout = 3600) {
  # A PLINK process was seen to finish its work (log 'End time') yet not exit, deadlocking the caller (2026-10-08, chr2 panel-2
  # --freq), and a 300 s timeout was too short on the shared, saturated machine. So the process is polled (processx): it succeeds
  # when its .log reaches 'End time' and it exits 0; if it lingers > 60 s after 'End time' it is killed and the run is accepted
  # only because the log is complete (logged); a hard timeout and a retry remain.
  exe <- plink_exe()
  args <- sub('^"(.*)"$', "\\1", args)          # processx passes arguments verbatim: drop the shell quotes callers add with q()
  for (k in seq_len(tries)) {
    so <- paste0(out_prefix, ".stdout"); se <- paste0(out_prefix, ".stderr"); lg <- paste0(out_prefix, ".log")
    unlink(c(so, se, lg))
    pr <- processx::process$new(exe, c(args, "--memory", "3000", "--threads", "1", "--out", out_prefix), stdout = so, stderr = se)
    t0 <- Sys.time(); t_end <- NA; status <- NA
    repeat {
      if (!pr$is_alive()) { status <- pr$get_exit_status(); break }
      done <- file.exists(lg) && any(grepl("^End time", readLines(lg, warn = FALSE)))
      if (done && is.na(t_end)) t_end <- Sys.time()
      if (!is.na(t_end) && difftime(Sys.time(), t_end, units = "secs") > 60) { pr$kill(); status <- "killed_after_end"; break }
      if (difftime(Sys.time(), t0, units = "secs") > timeout) { pr$kill(); status <- "timeout"; break }
      Sys.sleep(1)
    }
    done <- file.exists(lg) && any(grepl("^End time", readLines(lg, warn = FALSE)))
    if (done && (identical(status, 0L) || identical(status, "killed_after_end"))) {
      if (identical(status, "killed_after_end")) say("PLINK %s: log complete ('End time') but the process did not exit within 60 s; killed and accepted", basename(out_prefix))
      return(invisible(so))
    }
    say("PLINK attempt %d/%d for %s: status %s, log complete %s -> %s", k, tries, basename(out_prefix), status, done, if (k < tries) "retry" else "fail")
  }
  stop("PLINK failed after ", tries, " attempts: see ", out_prefix, ".log")
}

# ---- harmonisation units (re-implemented from l3fm__common.py; engine read, not run) ----
COMP <- c(A = "T", T = "A", C = "G", G = "C")
is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")
# +1 when (ea, oa) == (a1, a2); -1 when swapped; complement allowed only when NOT palindromic; 0 otherwise. Vectorised.
align_to_ref <- function(ea, oa, a1, a2) {
  pal <- is_pal(ea, oa)
  s <- integer(length(ea))
  s[ea == a1 & oa == a2] <- 1L
  s[s == 0L & ea == a2 & oa == a1] <- -1L
  cea <- unname(COMP[ea]); coa <- unname(COMP[oa])
  k <- s == 0L & !pal & !is.na(cea) & !is.na(coa)
  s[k & cea == a1 & coa == a2] <- 1L
  s[k & cea == a2 & coa == a1 & s == 0L] <- -1L
  s
}
REF_MAF_MIN <- 0.01; AF_DIFF_MAX <- 0.20; PALIN_MAF_MAX <- 0.40; P2_AF_DIFF_MAX <- 0.05; P1FR_AF_DIFF_MAX <- 0.10
n_hat <- function(se, f) {      # engine rule: median over 0.05 < MAF of 1 / (2 f (1 - f) se^2)
  k <- is.finite(f) & pmin(f, 1 - f) > 0.05 & se > 0
  if (sum(k) < 50) stop("n_hat: only ", sum(k), " usable SNPs")
  median(1 / (2 * f[k] * (1 - f[k]) * se[k]^2))
}
neff <- function(nca, nco) 4 * nca * nco / (nca + nco)

# ---- results accumulation (one TSV per script under work/rows/, merged by 09_collect.R) ----
ROWDIR <- file.path(WORK, "rows"); dir.create(ROWDIR, showWarnings = FALSE, recursive = TRUE)
RCOLS <- c("quantity_id", "prespec_ref", "estimand", "resource", "panel", "cov", "se_order", "value", "se", "ci_lo", "ci_hi", "p", "n", "unit", "notes")
rrow <- function(quantity_id, prespec_ref, estimand, resource = "", panel = "", cov = "", se_order = "", value = NA, se = NA,
                 ci_lo = NA, ci_hi = NA, p = NA, n = NA, unit = "", notes = "") {
  f <- function(x) if (is.numeric(x)) { if (is.na(x)) "" else formatC(x, digits = 10, format = "g") } else as.character(x)
  data.table(quantity_id = quantity_id, prespec_ref = prespec_ref, estimand = estimand, resource = resource, panel = panel,
             cov = cov, se_order = se_order, value = f(value), se = f(se), ci_lo = f(ci_lo), ci_hi = f(ci_hi), p = f(p),
             n = f(n), unit = unit, notes = notes)
}
write_rows <- function(dt, name) { fwrite(dt, file.path(ROWDIR, paste0(name, ".tsv")), sep = "\t"); say("rows written: %s (%d)", name, nrow(dt)) }
ASSERTS <- file.path(WORK, "asserts.tsv")
log_assert <- function(id, target, value, ok, derivation) {
  a <- data.table(assert_id = id, target = as.character(target), value = as.character(value), status = if (isTRUE(ok)) "PASS" else "ASSERT_FAILED",
                  derivation = derivation, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  fwrite(a, ASSERTS, sep = "\t", append = file.exists(ASSERTS))
  say("ASSERT %s %s: target %s, value %s (%s)", a$status, id, a$target, a$value, derivation)
  invisible(isTRUE(ok))
}

# Stage 2 step 2: with seven fits the (d) table has a zero cell (b = 0), so the conditional-MLE OR and its exact upper bound are +Inf by
# definition. +Inf is admitted ONLY for S6.d.or / S6.d.hi (and the kept equivalents) of a table whose b cell is "0"; any other non-finite
# numeric value is returned as bad. Used by 09_collect.R and mutation-tested in 13_mutation_tests.R.
nonfinite_bad <- function(dt, num = suppressWarnings(as.numeric(dt$value))) {
  inf_ok <- function(pre) { b <- dt[quantity_id == paste0(pre, ".b")]$value; length(b) == 1L && b == "0" }
  ok_inf <- !is.na(num) & is.infinite(num) & num > 0 & ((dt$quantity_id %in% c("S6.d.or", "S6.d.hi") & inf_ok("S6.d")) |
                                                       (dt$quantity_id %in% c("S6.d.kept.or", "S6.d.kept.hi") & inf_ok("S6.d.kept")))
  dt[!is.na(num) & !is.finite(num) & !ok_inf]
}
