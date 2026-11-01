# w8_export_ld.R (Stage 7 A2.2, S2-D47): export the LD rows for the five proxy pairs from the saved LD-proxy table.
# Reads  : the saved PLINK in-phase proxy table (RDS, 9 columns) + its .log (panel sample size) + the instruments file (alias column rsid_ld_panel)
#          + the replication/Steiger panel file (proxy_used, proxy_snps = 'lead<-proxy(r2=...)').
# Writes : w8_ld_rows.tsv (one row per proxy pair, full-precision r2 and the phased haplotype alleles) and w8_ld_panel.txt (n people from the .log).
# NOTHING IS TYPED: leads, proxies and panel ids come from the files above; w8_proxy_record.py pins the md5 of each and asserts them.
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
suppressMessages(library(data.table))
V3  <- .P31R_V3
LDD <- file.path(V3, "paper3_deposit_request_v5/files/data/ld_reference/1kg_eur")
RDS <- file.path(LDD, "instr_ld_proxies_phase1_hg38.rds")
LOG <- file.path(LDD, "instr_ld_proxies_phase1_hg38.log")
INS <- file.path(V3, "V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv")
PAN <- file.path(V3, "R2_rewrite/release/out/replication/replication_sensitivity__panel1_2.csv")
OUT <- file.path(V3, "R1_stage1/S1R1_R2_08_proxy_record")
ld  <- readRDS(RDS)
stopifnot(is.data.frame(ld), nrow(ld) == 42538L,
          identical(names(ld), c("instrument","proxy","r2","MAF_A","MAF_B","a_hap1","b_hap1","a_hap2","b_hap2")))
ins <- fread(INS)
stopifnot(all(c("SNP","gene","rsid_ld_panel") %in% names(ins)))
pan <- fread(PAN)
stopifnot(all(c("gene","outcome_id","outcome_name","proxy_used","proxy_snps") %in% names(pan)))
pp <- pan[proxy_used == TRUE]
m  <- regmatches(pp$proxy_snps, regexec("^(rs[0-9]+)<-(rs[0-9]+)[(]r2=([0-9.]+)[)]$", pp$proxy_snps))
stopifnot(all(lengths(m) == 4L), nrow(pp) == 5L)
res <- list()
for (k in seq_len(nrow(pp))) {
  lead <- m[[k]][2]; prox <- m[[k]][3]
  ii <- ins[SNP == lead & gene == pp$gene[k]]
  stopifnot(nrow(ii) == 1L)
  pid <- ii$rsid_ld_panel                      # the id the panel stores the lead under (alias), by the instruments file
  rr  <- ld[instrument == pid & proxy == prox]
  stopifnot(nrow(rr) == 1L)
  res[[k]] <- data.table(pair = paste0(pp$gene[k], "-", pp$outcome_name[k]), outcome_id = pp$outcome_id[k], lead_rsid = lead, panel_instrument_id = pid,
                         proxy_rsid = prox, r2 = sprintf("%.17g", rr$r2), MAF_A = sprintf("%.17g", rr$MAF_A), MAF_B = sprintf("%.17g", rr$MAF_B),
                         a_hap1 = rr$a_hap1, b_hap1 = rr$b_hap1, a_hap2 = rr$a_hap2, b_hap2 = rr$b_hap2, proxy_snps_text = pp$proxy_snps[k])
}
res <- rbindlist(res)
fwrite(res, file.path(OUT, "w8_ld_rows.tsv"), sep = "\t")
lg <- readLines(LOG, warn = FALSE)
nl <- grep("[0-9]+ people [(]", lg, value = TRUE)
stopifnot(length(nl) == 1L)
writeLines(c(paste0("log_line\t", nl), paste0("log_file\t", LOG)), file.path(OUT, "w8_ld_panel.txt"))
cat("w8_export_ld.R: wrote", nrow(res), "pair rows; panel line:", nl, "\n")
