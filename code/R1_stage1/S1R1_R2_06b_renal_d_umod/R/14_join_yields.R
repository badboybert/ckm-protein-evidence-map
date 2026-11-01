# 14_join_yields.R -- join yields per candidate (rule 7), from the extracts and the harmonisation stamps (nothing re-computed),
# plus a copy of the assertion log at the top level. Writes join_yields.tsv and asserts.tsv.
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
say("==== W6b 14_join_yields ====")
cand <- fread(file.path(WORK, "candidates_windows.tsv"))
Y <- rbindlist(lapply(cand$cid, function(cid) {
  st <- fromJSON(file.path(WORK, "harm", paste0(cid, "__stamp.json")))
  pf <- file.path(WORK, "extracts", paste0(cid, "__protein.tsv.gz"))
  np <- if (file.exists(pf)) nrow(fread(pf)) else NA_integer_
  dr <- st$drops; g <- function(k) if (is.null(dr[[k]])) 0L else as.integer(dr[[k]])
  joined <- if (is.na(np)) st$n_outcome_rows - g("dup_key_outcome") else np - g("nojoin_outcome") - g("ambiguous_outcome")
  refm <- joined - g("ref_absent") - g("ref_ambiguous")
  data.table(cid = cid, protein_rows_in_window = np, outcome_rows_in_cut = st$n_outcome_rows,
             joined_protein_outcome = if (is.na(np)) NA_integer_ else joined,
             join_yield_vs_protein = if (is.na(np)) NA_real_ else round(joined / np, 4),
             reference_matched = refm, after_qc_susie_set = st$n_susie, p1fr = st$n_p1fr, p1np = st$n_p1np, p2np = st$n_p2np,
             qc_drops = paste(names(dr), unlist(dr), sep = "=", collapse = "; "), panel2_drops = paste(names(st$drops_p2), unlist(st$drops_p2), sep = "=", collapse = "; "),
             rsid_concordance_outcome = st$rsid_concordance_outcome, n_P = st$n_P, n_O = st$n_O, n_O_source = st$n_O_source)
}))
fwrite(Y, file.path(OUT, "join_yields.tsv"), sep = "\t")
for (i in seq_len(nrow(Y))) say("join yield %-12s protein %s -> joined %s (%.1f%%) -> ref %d -> QC'd %d (susie) / %d p1fr / %d p1np / %d p2np; rsID concordance %s", Y$cid[i],
                                Y$protein_rows_in_window[i], Y$joined_protein_outcome[i], 100 * Y$join_yield_vs_protein[i], Y$reference_matched[i], Y$after_qc_susie_set[i],
                                Y$p1fr[i], Y$p1np[i], Y$p2np[i], Y$rsid_concordance_outcome[i])
file.copy(ASSERTS, file.path(OUT, "asserts.tsv"), overwrite = TRUE)
a <- fread(file.path(OUT, "asserts.tsv")); last <- a[, .SD[.N], by = assert_id]
say("asserts.tsv: %d log rows, %d distinct checks; final status: %d PASS, %d ASSERT_FAILED (%s)", nrow(a), nrow(last), sum(last$status == "PASS"), sum(last$status != "PASS"),
    paste(last[status != "PASS"]$assert_id, collapse = ", "))
