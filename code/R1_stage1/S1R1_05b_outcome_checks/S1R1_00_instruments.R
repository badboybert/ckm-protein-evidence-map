# Exposure-side and discovery-side inputs for the six predeclared pairs (no check-resource row is read here).
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
suppressPackageStartupMessages(library(data.table))
ER <- paste0(.P31R_V3, "/v3.2/_engine_root/results")
OD <- paste0(.P31R_V3, "/R1_stage1/S1R1_05b_outcome_checks")
ri <- as.data.table(readRDS(file.path(ER, "phase1_instruments/ukb_cis_instruments.rds")))
mr <- as.data.table(readRDS(file.path(ER, "phase2_mr_screen/mr_all_results.rds")))
cat("instrument cols:", paste(names(ri), collapse = ","), "\n")
cat("mr cols:", paste(names(mr), collapse = ","), " rows", nrow(mr), "\n")
pairs <- data.table(gene = c("PAM", "RAD23B", "ARG1", "PCSK9", "TGFB1"), outcome = c("T2D", "T2D", "T2D", "CAD", "CAD"))
oc <- intersect(c("outcome", "outcome_name", "disease"), names(mr))[1]
x <- merge(pairs, ri, by = "gene")
y <- mr[paste(gene, get(oc)) %in% paste(pairs$gene, pairs$outcome)]
print(x); print(y)
fwrite(x, file.path(OD, "S1R1_eur_instruments.tsv"), sep = "\t")
fwrite(y, file.path(OD, "S1R1_eur_discovery_mr.tsv"), sep = "\t")
cat("md5 instruments rds", unname(tools::md5sum(file.path(ER, "phase1_instruments/ukb_cis_instruments.rds"))), "\n")
cat("md5 mr_all rds", unname(tools::md5sum(file.path(ER, "phase2_mr_screen/mr_all_results.rds"))), "\n")
