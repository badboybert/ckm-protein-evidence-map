# 09_collect.R -- merge work/rows/*.tsv into results.tsv (declared columns); assert unique ids and no empty value where one is expected.
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
say("==== W6b 09_collect ====")
f <- sort(list.files(ROWDIR, pattern = "\\.tsv$", full.names = TRUE))
dt <- rbindlist(lapply(f, function(x) fread(x, colClasses = "character", na.strings = NULL)), fill = TRUE)
for (k in RCOLS) if (!k %in% names(dt)) dt[, (k) := ""]
dt <- dt[, ..RCOLS]
for (k in RCOLS) dt[is.na(get(k)), (k) := ""]
dup <- dt[duplicated(quantity_id)]$quantity_id
if (length(dup)) stop("duplicate quantity_id: ", paste(head(dup, 10), collapse = ", "))
empty <- dt[value == "" & !grepl("NOT RUN", notes)]
log_assert("results.no_empty_values", 0, nrow(empty), nrow(empty) == 0, paste("rows with an empty value and no NOT RUN note:", paste(head(empty$quantity_id, 10), collapse = ",")))
num <- suppressWarnings(as.numeric(dt$value))
bad <- nonfinite_bad(dt, num)
ok_inf <- !is.na(num) & is.infinite(num) & num > 0 & !(dt$quantity_id %in% bad$quantity_id)
log_assert("results.no_nonfinite_values", 0, nrow(bad), nrow(bad) == 0, sprintf("numeric values that are Inf/NaN, other than +Inf in S6.d(.kept).or/.hi of a table with b = 0 (%d such values admitted)", sum(ok_inf)))
fwrite(dt, file.path(OUT, "results.tsv"), sep = "\t", quote = FALSE)
say("results.tsv: %d rows from %d row files; md5 %s", nrow(dt), length(f), md5(file.path(OUT, "results.tsv")))
