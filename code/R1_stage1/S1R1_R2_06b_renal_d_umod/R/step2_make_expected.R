# HONESTY NOTE (2026-10-10): These tables pin the md5 of what the pipeline emitted at step 2; they are REGRESSION guards (reject any later
# change) and do NOT establish that the pinned values are correct. Independent anchors: W6b s_max - A_third + L1/L2 within 2.1e-8; W2/W4b
# s_max, per-fit s and text cells - none.
# step2_make_expected.R -- ONE-OFF generator of step2_expected_changes.tsv (2026-10-10). NOT part of run_all.R: the table is the recorded,
# reviewed expected-change list, so regenerating it from a new output would defeat the gate. Run only to re-record a deliberate change.
# It records the post-step-2 output (the one the Stage 2 step 2 report describes) and assigns one reason per cell by id pattern;
# every cell must get a reason, or the script stops.
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod/R/step2_gate_lib.R"))
D <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"); PRE <- file.path(D, "_pre_step2_20261010")
reason_res <- function(id, col, kind) {
  if (kind == "added") {
    if (grepl("_4fit$", id)) return("new row: the four-fit value kept as the labelled sensitivity")
    if (grepl("[.]s[.](protein|outcome)[.](p2np[.]503|susie[.]489|p1fr[.]489)$", id)) return("new row: per-fit s of one of the three added SuSiE fits")
    return(NA_character_)
  }
  if (grepl("^S6[.]d[.][A-Za-z0-9]+_[A-Za-z0-9]+[.](s_protein|s_outcome|s_max)$", id))
    return(switch(col, value = "s is now the maximum over all seven distinct SuSiE fits (was four)", estimand = "estimand text names the seven-fit set", notes = "notes name the seven-fit argmax set", NA_character_))
  if (grepl("^S6[.]d[.][A-Za-z0-9]+_[A-Za-z0-9]+[.]class_d$", id) && col %in% c("notes", "value")) return("class_d text/value follows the seven-fit s_max (NOTCH2-T2D crosses 0.10)")
  if (grepl("^S6[.]d[.](kept[.])?(a|b|c|d|or|lo|hi|p|reclassified)$", id)) return("NOTCH2-T2D reclassified (seven-fit s_max 0.1010 >= 0.10): the (d) table moves 2/1/66/91 to 2/0/66/92")
  if (grepl("^S6[.]slot_s[.]", id)) return("slot s is now the seven-fit maximum")
  NA_character_
}
reason_br <- function(col) switch(col, inputs = "ids of the _4fit sensitivity rows appended to inputs", filled_template = "template prints the seven-fit s_max and its argmax set",
                                  qualifiers = "qualifier text names both maxima (seven-fit, four-fit)", branch = "NOTCH2-T2D reclassified as not a candidate", NA_character_)
rows <- list()
for (f in c("results.tsv", "branches.tsv")) {
  nw <- g_rd(file.path(D, f)); ol <- g_rd(file.path(PRE, f)); kn <- g_keys(f, nw); ko <- g_keys(f, ol)
  add <- setdiff(kn, ko)
  if (length(add)) { r <- vapply(add, reason_res, "", col = "*", kind = "added"); if (f != "results.tsv" || anyNA(r)) stop("unreasoned added rows: ", paste(add[is.na(r)], collapse = ","))
    rows[[length(rows) + 1]] <- as.data.table(list(file = f, key = add, column = "*", new_md5 = g_rowmd5(nw[match(add, kn)]), kind = "added", reason = unname(r))) }
  common <- intersect(ko, kn); io <- match(common, ko); inn <- match(common, kn)
  for (cc in names(ol)) { if (f == "results.tsv" && cc == "quantity_id") next
    ch <- which(ol[[cc]][io] != nw[[cc]][inn]); if (!length(ch)) next
    rs <- if (f == "results.tsv") mapply(reason_res, common[ch], MoreArgs = list(col = cc, kind = "changed"), USE.NAMES = FALSE) else rep(reason_br(cc), length(ch))
    if (anyNA(rs)) stop("unreasoned changed cells in ", f, " ", cc, ": ", paste(common[ch][is.na(rs)], collapse = ","))
    rows[[length(rows) + 1]] <- as.data.table(list(file = f, key = common[ch], column = cc, new_md5 = g_md5(nw[[cc]][inn][ch]), kind = "changed", reason = unname(rs))) } }
X <- rbindlist(rows)
fwrite(X, file.path(D, "step2_expected_changes.tsv"), sep = "\t", quote = TRUE)
cat(sprintf("step2_expected_changes.tsv: %d entries (%d added rows, %d changed cells)\n", nrow(X), sum(X$kind == "added"), sum(X$kind == "changed")))
