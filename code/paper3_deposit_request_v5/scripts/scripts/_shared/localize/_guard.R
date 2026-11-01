# ============================================================================
# 防呆 / fail-loud guards — kill SILENT failures in the rerun pipeline.
# (User requirement 2026-06-30: "最不想看到靜默問題出現".)
#
# Principle — distinguish two kinds of empty:
#   * BENIGN empty  : one protein x outcome pair with no shared SNP. Normal,
#                     happens thousands of times in a proteome-wide screen.
#                     -> caller should `next`/skip. DO NOT guard these.
#   * UNEXPECTED empty: a WHOLE outcome / source / mediator returns NULL or 0
#                     rows, or a whole outcome is missing from the result set
#                     (file gone, wrong columns, API failed, no JWT).
#                     -> that is a SILENT BUG. ABORT loudly.
#
# stop() halts the run with a red error; we also append the reason to
# RERUN_ABORTED.log so the failure is never lost even if console scrolls.
# Source this AFTER 00_config.R (uses BASE_DIR if present).
# ============================================================================

.guard_log <- function() {
  bd <- if (exists("BASE_DIR")) get("BASE_DIR") else getwd()
  file.path(bd, "RERUN_ABORTED.log")
}

# Loud, persistent, terminal abort.
rerun_abort <- function(msg) {
  line <- sprintf("[ABORT %s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  # The log write is BEST EFFORT and the abort message now says so when it failed. It used to be
  # wrapped in try(silent=TRUE) beside a comment promising a persistent record, so an unwritable
  # RERUN_ABORTED.log left the reader believing a file existed that did not. The abort itself must
  # never depend on the log succeeding -- that is why it is still a try().
  .w <- try(cat(line, "\n", file = .guard_log(), append = TRUE), silent = TRUE)
  if (inherits(.w, "try-error")) line <- paste0(line, " [NOT logged: could not write ", .guard_log(), "]")
  stop(line, call. = FALSE)
}

# A whole outcome/source/mediator must not be silently NULL/empty.
# min_rows defaults to 1 (just non-empty) — deliberately NOT a tight floor, so it
# never false-aborts on legitimate variation; raise it only where you are sure.
guard_nonempty <- function(x, label, min_rows = 1L) {
  if (is.null(x)) rerun_abort(sprintf("'%s' returned NULL (file missing / wrong columns / API failed / no JWT?)", label))
  # dim(), not tryCatch(nrow()). nrow() on a plain vector RETURNS NULL rather than raising, so the
  # tryCatch never fired, `is.na(NULL)` is logical(0), and `if (logical(0))` threw R's generic
  # "argument is of length zero" -- bypassing rerun_abort() and its log entirely. The guard whose
  # whole job is to fail loudly and persistently was the one place that failed quietly.
  n <- if (is.null(dim(x))) NA_integer_ else nrow(x)
  if (is.na(n))        rerun_abort(sprintf("'%s' has no rows attribute (not a table?)", label))
  if (n < min_rows)    rerun_abort(sprintf("'%s' returned %d rows (< expected min %d) — likely a silent failure", label, n, min_rows))
  invisible(x)
}

# Every expected key must appear in `got`; a missing one = a whole outcome
# silently vanished from the results.
guard_coverage <- function(got, expected, what = "outcomes") {
  got      <- unique(as.character(got))
  expected <- unique(as.character(expected))
  miss <- setdiff(expected, got)
  if (length(miss)) rerun_abort(sprintf("%s missing from results: %s  (silent drop — expected %d, got %d)",
                                        what, paste(miss, collapse = ", "), length(expected), length(intersect(expected, got))))
  invisible(TRUE)
}

# Result outcome-id set must EXACTLY match the canonical expected set:
#   missing -> a whole outcome silently dropped (NULL/empty extract);
#   extra   -> a stale/old extract leaked in (mixed old+new directory).
guard_exact_set <- function(got, expected, what = "outcome ids") {
  got <- unique(as.character(got)); expected <- unique(as.character(expected))
  miss <- setdiff(expected, got); extra <- setdiff(got, expected)
  if (length(miss) || length(extra))
    rerun_abort(sprintf("%s set mismatch -> missing:[%s] | unexpected/stale:[%s]",
                what, paste(miss, collapse = ", "), paste(extra, collapse = ", ")))
  invisible(TRUE)
}
