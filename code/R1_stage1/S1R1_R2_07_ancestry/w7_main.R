#!/usr/bin/env Rscript
# W7 -- Paper 3.1 R2 Stage 1: pre-declaration section 7 (ancestry numbers and scales) and section 8 (post hoc numbers).
# Computed from the REDACTED pre-declaration (md5 a281c346...) and its redacted amendment; the full texts are only md5-checked here.
# Rerun everything:  Rscript "<archive>/code/R1_stage1/S1R1_R2_07_ancestry/w7_main.R"
# Writes ONLY into OUT. Reads every input in place. Streams DIAMANTE (1.3 GB) and CARDIoGRAMplusC4D (277 MB gz) in chunks,
# caching the instrument rows in OUT/extracts/ (resumable: a .done marker records the source md5 and line counts).
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
options(stringsAsFactors = FALSE, warn = 1, width = 200)
suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(digest) })

OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_07_ancestry")
RES_DIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
STAGING <- paste0(.P31R_V3, "/V4_local/staging/l3-finemap")
EXTRACTS <- file.path(OUT, "extracts")
dir.create(EXTRACTS, showWarnings = FALSE)
# Mutation tests (w7_mutations.R): W7_MUTATE=<id> injects one in-memory defect and redirects every output to OUT/mutations/<id>/;
# the cached outcome extracts are still read from OUT/extracts. Unset in the real run.
MUT <- Sys.getenv("W7_MUTATE", "")
if (nzchar(MUT)) { OUT <- file.path(OUT, "mutations", MUT); dir.create(OUT, recursive = TRUE, showWarnings = FALSE) }
LOGF <- file.path(OUT, "run.log")
T0 <- Sys.time()
cat("", file = LOGF)
say <- function(...) {
  s <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOGF, append = TRUE)
}
md5 <- function(p) unname(tools::md5sum(p))
f10 <- function(x) ifelse(is.na(x), "NA", sprintf("%.10g", x))
g3 <- function(x) as.character(signif(x, 3))
pl <- function(...) { a <- list(...); if (!length(a[[1]])) character(0) else do.call(paste0, a) }   # paste0 that stays empty on zero rows

say("W7 start. R", R.version.string, "| data.table", as.character(packageVersion("data.table")),
    "| jsonlite", as.character(packageVersion("jsonlite")), "| digest", as.character(packageVersion("digest")))
say("Seeds: none used. W7 draws no random numbers (no bootstrap is declared in sections 7-8); seed 20261007 is therefore not invoked.")
say("LD: none computed. Sections 7-8 need no LD, so PLINK19_paper6 is not called (rule 3 not applicable).")

# ------------------------------------------------------------------ governing documents (rule 1)
GOV <- list(
  briefs          = c(file.path(RES_DIR, "STAGE1_BRIEFS.md"), "40d0d52f7055542d774c03f1edf04e00"),
  prespec_redact  = c(file.path(RES_DIR, "P31R_R2_STAGE1_PRESPEC_REDACTED.md"), "a281c3464619b1b67c28e13e7392c41c"),
  amend1_redact   = c(file.path(RES_DIR, "P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), "f265d600973c43239bb8377e0cd39d33"),
  prespec_full    = c(paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), "3319dc39740d9f5b3b8822133d9108b7"),
  amend1_full     = c(paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), "377508832ea9bd384133e869c456c909"),
  gls_q           = c(file.path(RES_DIR, "gls_q.R"), "b7390fdf5b6f2f9c1d30c298d6afdafc"),
  resources       = c(file.path(RES_DIR, "RESOURCES.json"), "6d21238155a62fb58f71e05a01638dd0"),
  outcome_prespec = c(paste0(.P31R_V3, "/R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), "c9079f2de1ffe75ce31f8a1d597d1fda")
)
for (k in names(GOV)) {
  m <- md5(GOV[[k]][1]); say("GOV", k, m, GOV[[k]][1])
  if (!identical(m, GOV[[k]][2])) stop("GOVERNING md5 mismatch for ", k, ": ", m, " != ", GOV[[k]][2])
}
say("GOV all md5 asserted equal to the workflow prompt (full texts md5-checked only; not opened before results.tsv).")

# gls_q.R is not used by W7 (no GLS Q is declared in sections 7-8); its self-test is still logged (rule 6).
if (nzchar(MUT)) say("MUTATION RUN:", MUT, "-- outputs go to", OUT)
st <- if (nzchar(MUT)) "skipped in mutation runs" else
  tryCatch(system2(file.path(R.home("bin"), "Rscript"), c(shQuote(GOV$gls_q[1]), "--self-test"), stdout = TRUE, stderr = TRUE),
           error = function(e) paste("self-test failed to run:", conditionMessage(e)))
for (l in st) say("gls_q.R --self-test |", l)

# ------------------------------------------------------------------ resource map
RES <- fromJSON(GOV$resources[1], simplifyVector = FALSE)
rid <- vapply(RES, function(e) as.character(e$id), "")
stopifnot(length(RES) == 311, !anyDuplicated(rid))
res_get <- function(id) { i <- match(id, rid); if (is.na(i)) stop("RESOURCES lacks id ", id); RES[[i]] }
INPUTS <- list()
inp <- function(id) {
  e <- res_get(id); p <- e$path; m <- md5(p)
  say("INPUT", id, "| role", e$role, "| md5", m, "| expected", e$md5, "|", p)
  if (!identical(m, e$md5)) stop("INPUT md5 mismatch for ", id)
  if (!(e$role %in% c("raw_input", "definition"))) stop("W7 reads only raw_input/definition files as inputs; ", id, " is ", e$role)
  INPUTS[[id]] <<- c(p, m); p
}

# ------------------------------------------------------------------ staging snapshot and seen-output md5 at start (rule 2)
snap <- function(tag) {
  fs <- list.files(STAGING, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  d <- data.table(file = substring(fs, nchar(STAGING) + 2L), bytes = file.size(fs), md5 = unname(tools::md5sum(fs)))
  setorder(d, file)
  fwrite(d, file.path(OUT, paste0("staging_snapshot_", tag, ".tsv")), sep = "\t")
  dg <- digest(paste0(d$md5, "  ", d$file, collapse = "\n"), algo = "md5", serialize = FALSE)
  say("STAGING snapshot", tag, "| files", nrow(d), "| bytes", sum(d$bytes), "| digest of sorted md5 list", dg)
  d
}
SNAP0 <- snap("start")
POSTHOC_IDS <- c("POSTHOC_results", "POSTHOC_log")
POSTHOC_MD5_0 <- vapply(POSTHOC_IDS, function(i) { e <- res_get(i); m <- md5(e$path)
  say("SEEN-OUTPUT md5 at start (not read):", i, m, "| map", e$md5); if (!identical(m, e$md5)) stop("seen output changed before start: ", i); m }, "")
for (i in c("POSTHOC_spec", "POSTHOC_erratum", "POSTHOC_check_py", "ENG_eas__census_c6_census.R")) {
  e <- res_get(i); m <- md5(e$path); say("DEFINITION/ENGINE md5 (read for method):", i, e$role, m, "| map", e$md5)
  if (!identical(m, e$md5)) stop("md5 mismatch ", i)
}
say("SEEN SCRIPTS opened for method questions (all before any computation; reasons):")
say("  SEEN_wpE_mh_R ef0a393ae18e9bc4ff241717d255c3f2 -- which tiers enter the 'Supported-vs-Flagged' contrast and which exact stratified test (R mantelhaen.test exact=TRUE on tier x same-sign x outcome; Fisher for the unstratified p).")
say("  SEEN_wpE_scale_ckb_py 9f49411291f92793244eed81a8ed113a -- variant filter (0.05<af<0.95, sebeta>0) and what lo/hi are (10th/90th percentiles) for the implied-SD median.")
say("  SEEN_wpE_scale_eas_egfr_py 7397f786a987bc5d3ff57796c78fe5a4 -- the derived route the brief requires (se_y = se_mr_eas*|bx|, f = CKB af at the GRCh38 instrument position, n = n_total_variant).")
for (i in c("SEEN_wpE_mh_R", "SEEN_wpE_scale_ckb_py", "SEEN_wpE_scale_eas_egfr_py")) {
  e <- res_get(i); m <- md5(e$path); if (!identical(m, e$md5)) stop("seen script md5 mismatch ", i)
}

# ------------------------------------------------------------------ collectors
RESULTS <- list(); BRANCHES <- list()
addr <- function(quantity_id, prespec_ref, estimand, resource, value = NA, se = NA, ci_lo = NA, ci_hi = NA, p = NA, n = NA,
                 unit = "", notes = "", panel = "NA", cov = "NA", se_order = "NA") {
  v <- if (is.numeric(value)) f10(value) else as.character(value)
  RESULTS[[length(RESULTS) + 1L]] <<- data.table(quantity_id = quantity_id, prespec_ref = prespec_ref, estimand = estimand,
    resource = resource, panel = panel, cov = cov, se_order = se_order, value = v, se = f10(se), ci_lo = f10(ci_lo),
    ci_hi = f10(ci_hi), p = f10(p), n = ifelse(is.na(n), "NA", as.character(n)), unit = unit, notes = notes)
}
addb <- function(rule_ref, inputs, branch, filled_template, qualifiers) {
  BRANCHES[[length(BRANCHES) + 1L]] <<- data.table(rule_ref = rule_ref, inputs = inputs, branch = branch,
    filled_template = filled_template, qualifiers = qualifiers)
}
ASSERTS <- character(0)
assert_target <- function(id, computed, target, derivation, prespec_ref) {
  ok <- isTRUE(all.equal(as.numeric(computed), as.numeric(target)))
  tag <- if (ok) "ASSERT_PASSED" else "ASSERT_FAILED"
  if (!ok) ASSERTS <<- c(ASSERTS, paste0(id, ": computed ", computed, " vs target ", target))
  addr(paste0("ASSERT.", id), prespec_ref, paste("assertion target printed in the redacted pre-declaration:", derivation),
       "derived", value = computed, notes = paste0(tag, "; target ", target,
       " is visible in the redacted copy, so agreement is weak evidence (computed first, compared second)"))
  say(tag, id, "computed", computed, "target", target, "|", derivation)
}
count_lines <- function(p) length(readLines(p, warn = FALSE))
INSIDE <- c(0.80, 1.14)
scale_branch <- function(m) if (is.na(m)) "NOT_RUN" else if (m >= INSIDE[1] && m <= INSIDE[2]) "INSIDE" else "OUTSIDE"
scale_text <- function(b, m) switch(b,
  INSIDE  = sprintf("consistent with per-SD (median implied SD %.3f; acceptance interval 0.80-1.14)", m),
  OUTSIDE = sprintf("scale not established; magnitudes are not compared (median implied SD %.3f, outside 0.80-1.14)", m),
  NOT_RUN = "NOT RUN")

# =================================================================== SECTION 7(a): CKB pQTL empirical scale QC
say("=== S7 scale QC: CKB pQTL (implied SD = sebeta*sqrt(2*N*f*(1-f)); 0.05<af<0.95, sebeta>0; per-protein N from the manifest)")
ckb_e <- res_get("CKB_REGIONS_FOLDER"); ckb_dir <- ckb_e$path
tsvs <- sort(list.files(ckb_dir, pattern = "\\.tsv$"), method = "radix")
stopifnot(length(tsvs) == 38)
m5 <- md5(file.path(ckb_dir, tsvs))
folder_dig <- digest(paste0(m5, " *", tsvs, "\n", collapse = ""), algo = "md5", serialize = FALSE)
say("INPUT CKB_REGIONS_FOLDER | folder digest (md5 of 'md5 *name' lines, C-sorted, 38 .tsv) ", folder_dig, "| expected", ckb_e$md5,
    "| bytes", sum(file.size(file.path(ckb_dir, tsvs))), "| expected", ckb_e$bytes)
if (MUT == "folder_digest") m5[5] <- "00000000000000000000000000000000"   # mutation: one region file's md5 differs
folder_dig <- digest(paste0(m5, " *", tsvs, "\n", collapse = ""), algo = "md5", serialize = FALSE)
if (!identical(folder_dig, ckb_e$md5)) stop("CKB_REGIONS_FOLDER digest mismatch")
for (i in seq_along(tsvs)) say("  CKB file md5", m5[i], tsvs[i])
INPUTS[["CKB_REGIONS_FOLDER"]] <- c(ckb_dir, folder_dig)
p_man <- inp("CKB_PHENOTYPE_MANIFEST")
man <- fromJSON(paste(readLines(p_man, warn = FALSE), collapse = "\n"), simplifyVector = FALSE)
stopifnot(length(man) == 3149)
man_ph <- vapply(man, function(r) as.character(r$phenocode), "")
man_n <- vapply(man, function(r) if (is.null(r$num_samples)) NA_real_ else as.numeric(r$num_samples), 0)
stopifnot(!anyDuplicated(man_ph))
gene_files <- setdiff(tsvs, "FETCH_LOG.tsv"); stopifnot(length(gene_files) == 37)

scale_one <- function(path, N, label) {
  nl <- count_lines(path) - 1L
  d <- fread(path, sep = "\t", header = TRUE, colClasses = list(character = c("rsid", "id", "ref", "alt", "nearest_genes")))
  if (nrow(d) != nl || nl <= 0) stop("row count mismatch for ", path, ": fread ", nrow(d), " vs lines ", nl)
  if (!is.numeric(d$af) || !is.numeric(d$sebeta)) stop("non-numeric af/sebeta (placeholder token?) in ", path)
  k <- d[!is.na(af) & !is.na(sebeta) & af > 0.05 & af < 0.95 & sebeta > 0]
  s <- k$sebeta * sqrt(2 * N * k$af * (1 - k$af))
  data.table(file = label, N = N, n_total = nrow(d), na_af = mean(is.na(d$af)), na_se = mean(is.na(d$sebeta)), n_rows = nrow(k),
             median = median(s), lo = quantile(s, 0.10, names = FALSE, type = 7), hi = quantile(s, 0.90, names = FALSE, type = 7),
             min = min(s), max = max(s), median_N3971 = median(k$sebeta * sqrt(2 * 3971 * k$af * (1 - k$af))))
}
SC <- rbindlist(lapply(gene_files, function(f) {
  g <- sub("__.*", "", f); N <- man_n[match(g, man_ph)]
  if (is.na(N)) stop("no manifest N for ", g)
  cbind(gene = g, scale_one(file.path(ckb_dir, f), N, f))
}))
p_fgfr4 <- inp("CKB_FGFR4_cis_eas")
nfg <- man_n[match("FGFR4", man_ph)]; stopifnot(!is.na(nfg))
SC <- rbind(SC, cbind(gene = "FGFR4_pilot", scale_one(p_fgfr4, nfg, basename(p_fgfr4))))
SC[, branch := vapply(median, scale_branch, "")]
fwrite(SC, file.path(OUT, "s7_scale_ckb.tsv"), sep = "\t")
say("CKB per-file implied SD (gene, N, n_rows, median [p10, p90]):")
for (i in seq_len(nrow(SC))) say(sprintf("  %-12s N %d n_total %d n_rows %d median %.4f [%.4f, %.4f] min %.3f max %.3f | N3971 median %.4f | %s",
  SC$gene[i], SC$N[i], SC$n_total[i], SC$n_rows[i], SC$median[i], SC$lo[i], SC$hi[i], SC$min[i], SC$max[i], SC$median_N3971[i], SC$branch[i]))
for (i in seq_len(nrow(SC))) {
  q <- paste0("S7.scale.CKB_", SC$gene[i])
  src <- if (SC$gene[i] == "FGFR4_pilot") "CKB_FGFR4_cis_eas" else "CKB_REGIONS_FOLDER"
  est <- "median over variants (0.05<af<0.95, sebeta>0) of implied SD = sebeta*sqrt(2*N*af*(1-af)), N = manifest num_samples"
  addr(paste0(q, ".median"), "s7 scale QC", est, paste0(src, ":", SC$file[i]), value = SC$median[i], n = SC$n_rows[i], unit = "implied SD",
       notes = sprintf("N %d (manifest); n_total %d; branch %s; median with N=3971 (scoping's single N) %.4f", SC$N[i], SC$n_total[i], SC$branch[i], SC$median_N3971[i]))
  addr(paste0(q, ".lo"), "s7 scale QC", "10th percentile of implied SD (same variants)", paste0(src, ":", SC$file[i]), value = SC$lo[i], n = SC$n_rows[i], unit = "implied SD")
  addr(paste0(q, ".hi"), "s7 scale QC", "90th percentile of implied SD (same variants)", paste0(src, ":", SC$file[i]), value = SC$hi[i], n = SC$n_rows[i], unit = "implied SD")
  addr(paste0(q, ".n_rows"), "s7 scale QC", "variants entering the median (0.05<af<0.95, sebeta>0)", paste0(src, ":", SC$file[i]), value = SC$n_rows[i], n = SC$n_total[i], unit = "variants",
       notes = sprintf("rows in file %d; af NA rate %.4f; sebeta NA rate %.4f", SC$n_total[i], SC$na_af[i], SC$na_se[i]))
}
reg <- SC[gene != "FGFR4_pilot"]
mm <- median(reg$median); n_in <- sum(reg$branch == "INSIDE")
addr("S7.scale.CKB_all.median_of_medians", "s7 scale QC", "median across the 37 CKB region files of the per-file median implied SD", "CKB_REGIONS_FOLDER",
     value = mm, n = nrow(reg), unit = "implied SD", notes = sprintf("min %.4f max %.4f; files inside 0.80-1.14: %d of %d", min(reg$median), max(reg$median), n_in, nrow(reg)))
addr("S7.scale.CKB_all.min", "s7 scale QC", "minimum per-file median implied SD (37 region files)", "CKB_REGIONS_FOLDER", value = min(reg$median), n = nrow(reg), unit = "implied SD")
addr("S7.scale.CKB_all.max", "s7 scale QC", "maximum per-file median implied SD (37 region files)", "CKB_REGIONS_FOLDER", value = max(reg$median), n = nrow(reg), unit = "implied SD")
addr("S7.scale.CKB_all.n_inside", "s7 scale QC", "region files whose median implied SD lies in 0.80-1.14", "CKB_REGIONS_FOLDER", value = n_in, n = nrow(reg), unit = "files")
say(sprintf("CKB across 37 region files: median of medians %.4f (min %.4f, max %.4f); inside 0.80-1.14: %d of %d; FGFR4 pilot %.4f",
            mm, min(reg$median), max(reg$median), n_in, nrow(reg), SC[gene == "FGFR4_pilot", median]))
b_all <- scale_branch(mm)
addb("s7 wording rule (scale QC; acceptance interval 0.80-1.14)", sprintf("CKB median of 37 per-file medians %.4f; files inside %d/37; range %.3f-%.3f", mm, n_in, min(reg$median), max(reg$median)),
     paste0("CKB_pQTL_overall:", b_all), scale_text(b_all, mm),
     if (n_in < nrow(reg)) sprintf("per-file branch differs for %d file(s): %s", nrow(reg) - n_in, paste(reg[branch != "INSIDE", sprintf("%s %.3f", gene, median)], collapse = "; ")) else "every region file inside the interval")
for (i in seq_len(nrow(SC))) addb("s7 wording rule (scale QC), per file", sprintf("%s median %.4f, n_rows %d, N %d", SC$gene[i], SC$median[i], SC$n_rows[i], SC$N[i]),
                                  paste0("CKB_", SC$gene[i], ":", SC$branch[i]), scale_text(SC$branch[i], SC$median[i]), "")

# =================================================================== SECTION 7(a'): East Asian eGFR, derived route
say("=== S7 scale QC: outcome-side East Asian eGFR (EGFR_EAS_2P2M). Raw file NOT FOUND locally (RESOURCES EAS_eGFR_2P2M status ",
    res_get("EAS_eGFR_2P2M")$status, "); derived route.")
DEV_EGFR <- paste("DEVIATION: derived route; reproduces the scoping value by construction; not new evidence; N unconfirmed;",
                  "assumed SE order of se_mr_eas = first order (se_mr_eas = se_y/|bx|, so se_y = se_mr_eas*|bx|), as the scoping route assumed;",
                  "the order cannot be checked because the raw file is absent")
p_ax <- inp("TAB_EAS_axis")
stopifnot(identical(md5(p_ax), res_get("EAS_AXIS_SUPPTABLE")$md5))
ax <- fread(p_ax)
if (nrow(ax) != count_lines(p_ax) - 1L || nrow(ax) != 81 || ncol(ax) != 41) stop("EAS axis table shape")
if (MUT == "eas_sign") ax[gene == "ARG1" & outcome == "T2D", beta_mr_eas := -beta_mr_eas]   # mutation: label and estimate disagree
stopifnot(is.logical(ax$qc_evaluable), !anyNA(ax$qc_evaluable))
p_ins <- inp("TAB_instruments")
stopifnot(identical(md5(p_ins), res_get("EXP_UKBPPP_INSTRUMENT_TABLE")$md5))
ins <- fread(p_ins, colClasses = list(character = c("chrpos", "ukb_variant_id")))
if (nrow(ins) != 1964) stop("instrument table rows ", nrow(ins))
stopifnot(is.numeric(ins$beta), is.numeric(ins$se), !anyNA(ins$beta), !anyNA(ins$se))
ins_first <- ins[!duplicated(gene)]   # the scoping route: first instrument row per gene (pandas drop_duplicates("gene"))
eg <- ax[qc_evaluable == TRUE & eas_source == "EGFR_EAS_2P2M"]
say("EAS eGFR evaluable rows:", nrow(eg), "of", ax[eas_source == "EGFR_EAS_2P2M", .N], "EGFR_EAS_2P2M rows")
EGR <- rbindlist(lapply(seq_len(nrow(eg)), function(i) {
  r <- eg[i]; j <- match(r$gene, ins_first$gene); if (is.na(j)) stop("no instrument for ", r$gene)
  bx <- abs(ins_first$beta[j]); pos <- as.integer(sub(".*:", "", ins_first$chrpos[j]))
  same_snp <- identical(ins_first$SNP[j], r$instrument_rsID)
  n_gene_rows <- ins[gene == r$gene, .N]
  fl <- list.files(ckb_dir, pattern = paste0("^", r$gene, "__chr.*\\.tsv$"), full.names = TRUE)
  f <- NA_real_; n_hit <- 0L; allele_ok <- NA
  if (length(fl)) {
    cc <- fread(fl[1], sep = "\t", colClasses = list(character = c("rsid", "id", "ref", "alt", "nearest_genes")))
    h <- cc[position == pos]; n_hit <- nrow(h)
    if (n_hit) { f <- h$af[1]; allele_ok <- setequal(toupper(c(h$ref[1], h$alt[1])), toupper(c(ins_first$effect_allele[j], ins_first$other_allele[j]))) }
  }
  n <- r$n_total_variant; se_y <- r$se_mr_eas * bx
  sd <- if (is.na(n) || is.na(f)) NA_real_ else se_y * sqrt(2 * n * f * (1 - f))
  data.table(gene = r$gene, instrument_rsID = r$instrument_rsID, scoping_row_snp = ins_first$SNP[j], same_snp = same_snp,
             n_instrument_rows_gene = n_gene_rows, bx_abs = bx, pos_hg38 = pos, se_mr_eas = r$se_mr_eas, se_y = se_y,
             n_total_variant = n, ckb_file = if (length(fl)) basename(fl[1]) else NA_character_, ckb_rows_at_pos = n_hit,
             ckb_alleles_match_instrument = allele_ok, f_eas = f, sdY_hat = sd)
}))
fwrite(EGR, file.path(OUT, "s7_scale_eas_egfr_derived.tsv"), sep = "\t")
for (i in seq_len(nrow(EGR))) say(sprintf("  %-8s %-11s same_snp %s bx %.4f n %s f %s ckb_rows %d alleles_ok %s se_y %.5f sdY %s",
  EGR$gene[i], EGR$instrument_rsID[i], EGR$same_snp[i], EGR$bx_abs[i], EGR$n_total_variant[i], f10(EGR$f_eas[i]), EGR$ckb_rows_at_pos[i],
  EGR$ckb_alleles_match_instrument[i], EGR$se_y[i], f10(EGR$sdY_hat[i])))
ok <- EGR[!is.na(sdY_hat)]
say("EAS eGFR join yield: CKB frequency found for", nrow(ok), "of", nrow(EGR), "evaluable eGFR rows; instrument row = pair instrument for",
    sum(EGR$same_snp), "of", nrow(EGR))
eg_med <- median(ok$sdY_hat)
addr("S7.scale.EAS_eGFR_2P2M.median", "s7 scale QC", "median over evaluable eGFR pairs of implied SD = se_mr_eas*|bx|*sqrt(2*n_total_variant*f*(1-f)), f = CKB af at the instrument",
     "TAB_EAS_axis (EGFR_EAS_2P2M rows) + TAB_instruments + CKB_REGIONS_FOLDER", value = eg_med, n = nrow(ok), unit = "implied SD",
     notes = paste(DEV_EGFR, sprintf("| min %.4f max %.4f | rows with f found %d of %d", min(ok$sdY_hat), max(ok$sdY_hat), nrow(ok), nrow(EGR))))
addr("S7.scale.EAS_eGFR_2P2M.lo", "s7 scale QC", "10th percentile of the derived implied SD", "derived", value = quantile(ok$sdY_hat, 0.1, names = FALSE), n = nrow(ok), unit = "implied SD", notes = DEV_EGFR)
addr("S7.scale.EAS_eGFR_2P2M.hi", "s7 scale QC", "90th percentile of the derived implied SD", "derived", value = quantile(ok$sdY_hat, 0.9, names = FALSE), n = nrow(ok), unit = "implied SD", notes = DEV_EGFR)
addr("S7.scale.EAS_eGFR_2P2M.min", "s7 scale QC", "minimum derived implied SD", "derived", value = min(ok$sdY_hat), n = nrow(ok), unit = "implied SD", notes = DEV_EGFR)
addr("S7.scale.EAS_eGFR_2P2M.max", "s7 scale QC", "maximum derived implied SD", "derived", value = max(ok$sdY_hat), n = nrow(ok), unit = "implied SD", notes = DEV_EGFR)
addr("S7.scale.EAS_eGFR_2P2M.n_rows", "s7 scale QC", "evaluable eGFR pairs with a CKB frequency at the instrument", "derived", value = nrow(ok), n = nrow(EGR), unit = "pairs", notes = DEV_EGFR)
b_eg <- scale_branch(eg_med)
addb("s7 wording rule (scale QC; acceptance interval 0.80-1.14)", sprintf("EAS eGFR derived median %.4f over %d pairs (min %.3f, max %.3f)", eg_med, nrow(ok), min(ok$sdY_hat), max(ok$sdY_hat)),
     paste0("EAS_eGFR_2P2M:", b_eg), scale_text(b_eg, eg_med), DEV_EGFR)

# =================================================================== SECTION 7(b): outcome-stratified exact test (D12-A) and opposite-sign counts
say("=== S7 outcome-stratified exact test of the Supported-vs-Flagged sign contrast (East Asian evaluable pairs)")
p_t81 <- inp("TAB_tiered81")
t81 <- fread(p_t81)
if (nrow(t81) != 81 || nrow(t81) != count_lines(p_t81) - 1L) stop("tiered81 rows")
k_ax <- merge(ax, t81[, .(gene, outcome, t81_tier = evidence_tier, mr_beta)], by = c("gene", "outcome"), all.x = TRUE)
say("JOIN EAS axis x tiered81 on gene+outcome: matched", sum(!is.na(k_ax$mr_beta)), "of", nrow(ax), "| key null rates gene", mean(is.na(ax$gene)), "outcome", mean(is.na(ax$outcome)))
stopifnot(sum(!is.na(k_ax$mr_beta)) == 81, all(k_ax$t81_tier == k_ax$evidence_tier))
ev <- k_ax[qc_evaluable == TRUE]
stopifnot(all(ev$direction_class %in% c("SAME", "OPPOSITE")), !anyNA(ev$beta_mr_eas))
ev[, same_recomputed := sign(beta_mr_eas) == sign(mr_beta)]
n_dir_disagree <- ev[same_recomputed != (direction_class == "SAME"), .N]
say("direction_class recomputed from sign(beta_mr_eas) vs sign(discovery mr_beta): disagreements", n_dir_disagree, "of", nrow(ev))
stopifnot(n_dir_disagree == 0, all(ev$eas_direction_same == (ev$direction_class == "SAME")))
n_conv <- ev[evidence_tier == "Convergent", .N]
say("evaluable by tier:", paste(ev[, .N, by = evidence_tier][, paste(evidence_tier, N)], collapse = "; "), "| Convergent evaluable:", n_conv)
e2 <- ev[evidence_tier %in% c("Supported", "Flagged")]
e2[, tierf := factor(evidence_tier, levels = c("Supported", "Flagged"))]
e2[, samef := factor(direction_class == "SAME", levels = c(TRUE, FALSE), labels = c("same", "opposite"))]
e2[, outf := factor(outcome)]
tab2 <- table(e2$tierf, e2$samef)
ft <- fisher.test(tab2)
tab3 <- table(e2$tierf, e2$samef, e2$outf)
mh <- mantelhaen.test(tab3, exact = TRUE)
# check: only strata with both tiers AND both sign classes carry information; the exact conditional p must be identical
inf <- dimnames(tab3)[[3]][apply(tab3, 3, function(m) all(rowSums(m) > 0) && all(colSums(m) > 0))]
mh_inf <- mantelhaen.test(tab3[, , inf, drop = FALSE], exact = TRUE)
say("informative strata:", paste(inf, collapse = ","), "| p all strata", f10(mh$p.value), "| p informative strata", f10(mh_inf$p.value))
stopifnot(isTRUE(all.equal(mh$p.value, mh_inf$p.value)))
st_tab <- as.data.table(as.data.frame(tab3)); setnames(st_tab, c("tier", "sign", "outcome", "n"))
fwrite(st_tab, file.path(OUT, "s7_strat_table.tsv"), sep = "\t")
for (o in dimnames(tab3)[[3]]) say(sprintf("  stratum %-7s Supported same %d opp %d | Flagged same %d opp %d", o, tab3["Supported", "same", o], tab3["Supported", "opposite", o],
                                           tab3["Flagged", "same", o], tab3["Flagged", "opposite", o]))
say(sprintf("pooled: Supported same %d opp %d | Flagged same %d opp %d | Fisher p %s OR %s | stratified exact p %s, conditional MLE OR %s (95%% CI %s-%s)",
  tab2["Supported", "same"], tab2["Supported", "opposite"], tab2["Flagged", "same"], tab2["Flagged", "opposite"], f10(ft$p.value), f10(ft$estimate),
  f10(mh$p.value), f10(mh$estimate), f10(mh$conf.int[1]), f10(mh$conf.int[2])))
nS <- sum(tab2["Supported", ]); nF <- sum(tab2["Flagged", ])
addr("S7.strat.p_exact", "s7; s11 D12-A", "exact conditional test of a common OR = 1 (tier Supported vs Flagged x same vs opposite sign), stratified by outcome; East Asian qc-evaluable pairs",
     "TAB_EAS_axis", value = mh$p.value, p = mh$p.value, ci_lo = mh$conf.int[1], ci_hi = mh$conf.int[2], n = nrow(e2), unit = "p",
     notes = sprintf("R mantelhaen.test(exact=TRUE); conditional MLE OR (same sign, Supported vs Flagged) %s; strata %s; informative strata %s (identical p); Convergent evaluable pairs: %d (none enter); post hoc (D12)",
                     f10(mh$estimate), paste(dimnames(tab3)[[3]], collapse = ","), paste(inf, collapse = ","), n_conv))
addr("S7.strat.p_unstrat", "s7; s11 D12-A (printed for completeness; D12-A forbids reporting it alone)", "Fisher exact test, pooled 2x2 (tier Supported vs Flagged x same vs opposite sign); East Asian qc-evaluable pairs",
     "TAB_EAS_axis", value = ft$p.value, p = ft$p.value, ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2], n = nrow(e2), unit = "p",
     notes = sprintf("conditional MLE OR %s; never reported unstratified (s7 wording rule)", f10(ft$estimate)))
addr("S7.strat.Supported.same", "s7", "Supported East Asian evaluable pairs with the discovery sign", "TAB_EAS_axis", value = tab2["Supported", "same"], n = nS, unit = "pairs")
addr("S7.strat.Flagged.same", "s7", "Flagged East Asian evaluable pairs with the discovery sign", "TAB_EAS_axis", value = tab2["Flagged", "same"], n = nF, unit = "pairs")
addb("s11 D12 (D12-A)", sprintf("Supported %d/%d same sign; Flagged %d/%d; stratified exact p %s; unstratified p %s (not reported)",
       tab2["Supported", "same"], nS, tab2["Flagged", "same"], nF, f10(mh$p.value), f10(ft$p.value)), "D12-A",
     sprintf("D12-A fill: Supported %d/%d same sign vs Flagged %d/%d same sign; outcome-stratified exact p %s; label 'post hoc'",
             tab2["Supported", "same"], nS, tab2["Flagged", "same"], nF, g3(mh$p.value)),
     "post hoc; reported only with the outcome-stratified exact p (never unstratified)")

# opposite-sign counts (section 8 'Also derived': 7 opposite-sign, 43 distinct lead variants)
n_opp <- ev[direction_class == "OPPOSITE", .N]; n_lead <- uniqueN(ev$instrument_rsID)
addr("S7.opp.n_opposite", "s7; s8 Also derived", "East Asian qc-evaluable pairs whose Wald ratio has the opposite sign to discovery", "TAB_EAS_axis", value = n_opp, n = nrow(ev), unit = "pairs",
     notes = paste("pairs:", paste(ev[direction_class == "OPPOSITE", pl(gene, "-", outcome, "(", evidence_tier, ")")], collapse = ", ")))
addr("S7.opp.n_evaluable", "s7", "East Asian qc-evaluable pairs", "TAB_EAS_axis", value = nrow(ev), n = nrow(ax), unit = "pairs")
addr("S7.opp.n_distinct_leads", "s7; s8 Also derived", "distinct instrument (lead) variants among East Asian qc-evaluable pairs", "TAB_EAS_axis", value = n_lead, n = nrow(ev), unit = "variants")
for (t in c("Convergent", "Supported", "Flagged")) {
  d <- ev[evidence_tier == t]
  addr(paste0("S7.opp.", t), "s7", paste("opposite-sign East Asian evaluable pairs,", t), "TAB_EAS_axis", value = d[direction_class == "OPPOSITE", .N], n = nrow(d), unit = "pairs")
}
for (o in sort(unique(ev$outcome))) {
  d <- ev[outcome == o]
  addr(paste0("S7.opp.outcome.", o), "s7", paste("opposite-sign East Asian evaluable pairs, outcome", o), "TAB_EAS_axis", value = d[direction_class == "OPPOSITE", .N], n = nrow(d), unit = "pairs")
}
assert_target("S8.eas_opposite_sign", n_opp, 7, "count of direction_class OPPOSITE among the 52 qc-evaluable East Asian pairs", "s8 Also derived")
assert_target("S8.eas_distinct_leads", n_lead, 43, "distinct instrument_rsID among the qc-evaluable East Asian pairs", "s8 Also derived")

# =================================================================== SECTION 7(c): CKB-free composition
say("=== S7 CKB-free composition (census rows whose outcome resource includes CKB, replaced by their CKB-free rows)")
p_cen <- inp("TAB_results_census"); stopifnot(identical(md5(p_cen), res_get("EAS_CENSUS_RESULTS")$md5))
p_ckf <- inp("TAB_results_ckbfree")
cen <- fread(p_cen, sep = "\t"); ckf <- fread(p_ckf, sep = "\t")
if (nrow(cen) != count_lines(p_cen) - 1L || nrow(cen) != 48) stop("census rows")
if (nrow(ckf) != count_lines(p_ckf) - 1L || nrow(ckf) != 13) stop("ckbfree rows")
tb <- function(x) as.character(x) %in% c("TRUE", "True", "true")
cen_ok <- cen[status == "ok"]
say("census rows", nrow(cen), "| status ok", nrow(cen_ok), "| not ok:", paste(cen[status != "ok", paste(gene, outcome, status)], collapse = "; "))
stopifnot(!anyNA(cen_ok$ckb_overlap))
ov <- cen_ok[tb(ckb_overlap)]; nov <- cen_ok[!tb(ckb_overlap)]
stopifnot(all(!tb(ckf$ckb_overlap)), all(ckf$status == "ok"))
if (MUT == "ckbfree_drop") ckf <- ckf[-1]   # mutation: one replacement row missing
jj <- merge(ov[, .(gene, outcome, tier, cen_resource = resource_label, cen_support = tb(support), cen_support6 = tb(`support_p12_1e-6`), cen_mr_beta = mr_beta, cen_mr_p = mr_p, cen_H4 = PP.H4)],
            ckf[, .(gene, outcome, ckf_resource = resource_label, ckf_support = tb(support), ckf_support6 = tb(`support_p12_1e-6`), ckf_mr_beta = mr_beta, ckf_mr_p = mr_p, ckf_H4 = PP.H4)],
            by = c("gene", "outcome"), all = TRUE)
say("JOIN census CKB-overlap rows x ckbfree rows on gene+outcome: overlap rows", nrow(ov), "| ckbfree rows", nrow(ckf), "| matched", jj[!is.na(cen_H4) & !is.na(ckf_H4), .N])
stopifnot(nrow(jj) == nrow(ov), nrow(jj) == nrow(ckf), !anyNA(jj$cen_H4), !anyNA(jj$ckf_H4))
fwrite(jj, file.path(OUT, "s7_ckbfree_pairs.tsv"), sep = "\t")
comp <- rbind(nov[, .(gene, outcome, tier, resource_label, support = tb(support), support6 = tb(`support_p12_1e-6`), ckb_free_source = "census (CKB-free resource)")],
              merge(ckf, ov[, .(gene, outcome, tier)], by = c("gene", "outcome"))[, .(gene, outcome, tier, resource_label, support = tb(support), support6 = tb(`support_p12_1e-6`), ckb_free_source = "CKB-free replacement")])
stopifnot(nrow(comp) == nrow(cen_ok))
fwrite(comp, file.path(OUT, "s7_ckbfree_composite.tsv"), sep = "\t")
cs <- cen_ok[tb(support)]; ks <- comp[support == TRUE]
addr("S7.ckbfree.n_census_tested", "s7 CKB-free composition", "East Asian census pairs with status ok", "TAB_results_census", value = nrow(cen_ok), n = nrow(cen), unit = "pairs")
addr("S7.ckbfree.n_ckb_overlap", "s7 CKB-free composition", "census pairs whose outcome resource includes CKB (ckb_overlap TRUE)", "TAB_results_census", value = nrow(ov), n = nrow(cen_ok), unit = "pairs",
     notes = paste("resources:", paste(unique(ov$resource_label), collapse = " | "), "| pairs:", paste(ov[, pl(gene, "-", outcome)], collapse = ", ")))
addr("S7.ckbfree.n_replaced", "s7 CKB-free composition", "CKB-overlap pairs with a CKB-free replacement row", "TAB_results_ckbfree", value = nrow(jj), n = nrow(ov), unit = "pairs",
     notes = paste("replacement resources:", paste(unique(ckf$resource_label), collapse = " | ")))
addr("S7.ckbfree.census_supported", "s7 CKB-free composition", "census pairs colocalization-supported (support TRUE, p12 1e-5)", "TAB_results_census", value = nrow(cs), n = nrow(cen_ok), unit = "pairs",
     notes = paste(cs[, pl(gene, "-", outcome, "(", tier, ifelse(tb(ckb_overlap), ", CKB-overlap", ""), ")")], collapse = ", "))
addr("S7.ckbfree.census_supported_ckbfree_resources", "s7 CKB-free composition", "census-supported pairs whose resource is already CKB-free", "TAB_results_census", value = nrow(cs[!tb(ckb_overlap)]), n = nrow(cs), unit = "pairs")
addr("S7.ckbfree.composite_supported", "s7 CKB-free composition", "supported pairs (p12 1e-5) when every CKB-overlap resource is replaced by its CKB-free row", "TAB_results_census + TAB_results_ckbfree",
     value = nrow(ks), n = nrow(comp), unit = "pairs", notes = paste(ks[, pl(gene, "-", outcome, "(", tier, "; ", ckb_free_source, ")")], collapse = ", "))
addr("S7.ckbfree.composite_supported_p12_1e6", "s7 CKB-free composition", "supported pairs at p12 1e-6 in the CKB-free composite", "TAB_results_census + TAB_results_ckbfree",
     value = comp[support6 == TRUE, .N], n = nrow(comp), unit = "pairs", notes = paste(comp[support6 == TRUE, pl(gene, "-", outcome)], collapse = ", "))
addr("S7.ckbfree.census_supported_p12_1e6", "s7 CKB-free composition", "census pairs supported at p12 1e-6", "TAB_results_census",
     value = cen_ok[tb(`support_p12_1e-6`), .N], n = nrow(cen_ok), unit = "pairs", notes = paste(cen_ok[tb(`support_p12_1e-6`), pl(gene, "-", outcome)], collapse = ", "))
gained <- jj[ckf_support & !cen_support, pl(gene, "-", outcome)]; lost <- jj[cen_support & !ckf_support, pl(gene, "-", outcome)]
addr("S7.ckbfree.replaced_gained", "s7 CKB-free composition", "CKB-overlap pairs supported only in the CKB-free replacement", "TAB_results_ckbfree", value = length(gained), n = nrow(jj), unit = "pairs", notes = paste(gained, collapse = ", "))
addr("S7.ckbfree.replaced_lost", "s7 CKB-free composition", "CKB-overlap pairs supported only with the CKB-including resource", "TAB_results_ckbfree", value = length(lost), n = nrow(jj), unit = "pairs", notes = paste(lost, collapse = ", "))
addr("S7.ckbfree.replaced_sign_flips", "s7 CKB-free composition", "replaced pairs whose MR beta changes sign (CKB-including vs CKB-free resource)", "TAB_results_census + TAB_results_ckbfree",
     value = jj[sign(cen_mr_beta) != sign(ckf_mr_beta), .N], n = nrow(jj), unit = "pairs", notes = paste(jj[sign(cen_mr_beta) != sign(ckf_mr_beta), pl(gene, "-", outcome)], collapse = ", "))
for (t in c("Convergent", "Supported", "Flagged")) addr(paste0("S7.ckbfree.composite_supported.", t), "s7 CKB-free composition", paste("CKB-free composite supported pairs, tier", t),
     "TAB_results_census + TAB_results_ckbfree", value = ks[tier == t, .N], n = comp[tier == t, .N], unit = "pairs")
say("CKB-free composition: census supported", nrow(cs), "| composite supported", nrow(ks), "| gained", paste(gained, collapse = ","), "| lost", paste(lost, collapse = ","))

# =================================================================== SECTION 8: post hoc numbers (re-derived in R; the engine is NOT executed)
say("=== S8 post hoc external outcome check of the 28 colocalization-supported T2D/CAD pairs (re-derivation; P31R_posthoc_check.py not executed)")
p_mr <- inp("TAB_mr_all_results")
mr <- fread(p_mr)
if (nrow(mr) != count_lines(p_mr) - 1L || nrow(mr) != 14537) stop("mr_all_results rows")
pairs <- t81[outcome %in% c("T2D", "CAD")]
stopifnot(nrow(pairs) == 28, !any(grepl("[;, ]", pairs$lead_snps)))
pairs <- merge(pairs, ins[, .(gene, SNP, effect_allele, other_allele, bx = beta, sx = se, eaf, ukb_variant_id, chrpos)],
               by.x = c("gene", "lead_snps"), by.y = c("gene", "SNP"), all.x = TRUE, sort = FALSE)
say("JOIN 28 pairs x instruments on gene+lead SNP: matched", sum(!is.na(pairs$bx)), "of 28")
stopifnot(nrow(pairs) == 28, !anyNA(pairs$bx))
mrw <- mr[method %in% c("Wald ratio", "Inverse variance weighted"), .(gene, outcome = outcome_name, disc_b = b, disc_se = se, method, n_snps, any_proxy)]
pairs <- merge(pairs, mrw, by = c("gene", "outcome"), all.x = TRUE, sort = FALSE)
say("JOIN pairs x mr_all_results (Wald ratio / IVW) on gene+outcome: rows", nrow(pairs), "| matched", sum(!is.na(pairs$disc_b)),
    "| methods", paste(names(table(pairs$method)), table(pairs$method), collapse = "; "))
stopifnot(nrow(pairs) == 28, !anyNA(pairs$disc_b), !anyNA(pairs$disc_se))
say("discovery b (mr_all_results) vs tiered81 mr_beta (printed to 6 decimals): max abs diff", f10(max(abs(pairs$disc_b - pairs$mr_beta))))
stopifnot(max(abs(pairs$disc_b - pairs$mr_beta)) < 1e-6)   # tiered81 rounds mr_beta to 6 decimals (max rounding error 5e-7)
outc <- readLines(GOV$outcome_prespec[1], warn = FALSE, encoding = "UTF-8")
predecl <- c("PAM|T2D", "RAD23B|T2D", "ARG1|T2D", "PCSK9|CAD", "TGFB1|CAD")   # outcome-check pre-specification s1 (five European pairs; FGFR4-T2D is East Asian)
for (pp in predecl) {   # each pre-declared pair must be a row of that document's s1 table ("| GENE<en dash>OUTCOME |")
  g <- sub("\\|.*", "", pp); o <- sub(".*\\|", "", pp)
  if (!any(startsWith(outc, paste0("| ", g)) & grepl(paste0(o, " |"), outc, fixed = TRUE))) stop("pre-declared pair not found in the outcome-check prespec: ", pp)
}
pairs[, predeclared := paste(gene, outcome, sep = "|") %in% predecl]
pairs[, prior_robust := `PP.H4_p12_1e6` >= 0.80]
pairs[, pos37 := as.integer(sapply(strsplit(ukb_variant_id, ":"), `[`, 2))]
rs <- unique(pairs$lead_snps)
say("pairs 28 | T2D", pairs[outcome == "T2D", .N], "CAD", pairs[outcome == "CAD", .N], "| predeclared", sum(pairs$predeclared),
    "| distinct lead variants", length(rs), "| tiers", paste(names(table(pairs$evidence_tier)), table(pairs$evidence_tier), collapse = "; "))

# ---- streaming extracts (checkpointed; one unit = one source file)
extract_done <- function(tag, src_md5) {
  dn <- file.path(EXTRACTS, paste0(tag, ".done"))
  if (!file.exists(dn)) return(NULL)
  meta <- readLines(dn, warn = FALSE)
  if (!any(meta == paste("src_md5", src_md5)) || !any(meta == paste("rs", paste(sort(rs), collapse = ",")))) return(NULL)
  meta
}
stream_diamante <- function(src, src_md5) {
  tag <- "DIAMANTE_w7"; out <- file.path(EXTRACTS, paste0(tag, ".tsv"))
  meta <- extract_done(tag, src_md5)
  if (!is.null(meta)) { say("CHECKPOINT reuse", tag, paste(meta, collapse = " | ")); return(fread(out, sep = "\t", colClasses = "character")) }
  t1 <- Sys.time(); con <- file(src, open = "r"); on.exit(close(con))
  hdr <- readLines(con, n = 1L)
  stopifnot(identical(hdr, "rsID\tBP\tSNP\tCHR\tEA\tNEA\tEAF\tBeta\tSE\tPvalue"))
  n <- 1L; ch <- 0L; hits <- character(0)
  repeat {
    x <- readLines(con, n = 1000000L); if (!length(x)) break
    n <- n + length(x); ch <- ch + 1L
    f1 <- substr(x, 1L, regexpr("\t", x, fixed = TRUE) - 1L)
    hits <- c(hits, x[f1 %in% rs])
    if (ch %% 5L == 0L) say(sprintf("  DIAMANTE progress: %d lines read, %d hits, %.0f s", n, length(hits), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  }
  sp <- strsplit(hits, "\t", fixed = TRUE); stopifnot(all(lengths(sp) == 10L))
  d <- as.data.table(do.call(rbind, sp)); setnames(d, strsplit(hdr, "\t")[[1]])
  fwrite(d, out, sep = "\t")
  writeLines(c(paste("src", src), paste("src_md5", src_md5), paste("rs", paste(sort(rs), collapse = ",")), paste("lines_incl_header", n),
               paste("hits", nrow(d)), paste("seconds", round(as.numeric(difftime(Sys.time(), t1, units = "secs"))))), file.path(EXTRACTS, paste0(tag, ".done")))
  say("DIAMANTE streamed:", n, "lines incl header;", nrow(d), "rows at the", length(rs), "instruments;", round(as.numeric(difftime(Sys.time(), t1, units = "secs"))), "s")
  d
}
stream_vcf <- function(src, src_md5) {
  tag <- "CARDIOGRAM_w7"; out <- file.path(EXTRACTS, paste0(tag, ".tsv"))
  meta <- extract_done(tag, src_md5)
  if (!is.null(meta)) { say("CHECKPOINT reuse", tag, paste(meta, collapse = " | ")); return(fread(out, sep = "\t", colClasses = "character")) }
  t1 <- Sys.time(); con <- gzfile(src, open = "r"); on.exit(close(con))
  nrec <- 0L; nhdr <- 0L; ch <- 0L; hits <- character(0)
  repeat {
    x <- readLines(con, n = 1000000L); if (!length(x)) break
    ch <- ch + 1L
    h <- startsWith(x, "#"); nhdr <- nhdr + sum(h); x <- x[!h]; nrec <- nrec + length(x)
    id <- sub("^[^\t]*\t[^\t]*\t([^\t]*)\t.*$", "\\1", x, perl = TRUE)
    hits <- c(hits, x[id %in% rs])
    if (ch %% 3L == 0L) say(sprintf("  CARDIoGRAM progress: %d records, %d hits, %.0f s", nrec, length(hits), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  }
  sp <- strsplit(hits, "\t", fixed = TRUE); stopifnot(all(lengths(sp) == 10L))
  d <- as.data.table(do.call(rbind, sp)); setnames(d, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"))
  fwrite(d, out, sep = "\t")
  writeLines(c(paste("src", src), paste("src_md5", src_md5), paste("rs", paste(sort(rs), collapse = ",")), paste("records", nrec), paste("header_lines", nhdr),
               paste("hits", nrow(d)), paste("seconds", round(as.numeric(difftime(Sys.time(), t1, units = "secs"))))), file.path(EXTRACTS, paste0(tag, ".done")))
  say("CARDIoGRAM streamed:", nrec, "records,", nhdr, "header lines;", nrow(d), "rows at the instruments;", round(as.numeric(difftime(Sys.time(), t1, units = "secs"))), "s")
  d
}
p_dia <- inp("OUT_DIAMANTE2018_noUKB_EUR_rsid"); dia <- stream_diamante(p_dia, INPUTS[["OUT_DIAMANTE2018_noUKB_EUR_rsid"]][2])
p_car <- inp("OUT_CARDIOGRAMplusC4D_ieu-a-7"); car <- stream_vcf(p_car, INPUTS[["OUT_CARDIOGRAMplusC4D_ieu-a-7"]][2])
dmeta <- readLines(file.path(EXTRACTS, "DIAMANTE_w7.done")); cmeta <- readLines(file.path(EXTRACTS, "CARDIOGRAM_w7.done"))
n_dlines <- as.integer(sub("lines_incl_header ", "", grep("^lines_incl_header", dmeta, value = TRUE)))
n_crec <- as.integer(sub("records ", "", grep("^records", cmeta, value = TRUE)))
say("TRUNCATION guards: DIAMANTE lines incl header", n_dlines, "= data rows", n_dlines - 1L, "(RESOURCES: rows=21,508,158; wc -l on the CRLF file gives 21,508,159 lines incl header) | CARDIoGRAM records", n_crec, "(RESOURCES: 8,479,693 records)")
stopifnot(n_dlines - 1L == 21508158L, n_crec == 8479693L)
# outcome rows in a common shape: rsid, chr, pos37, ea, oa, b, se, eaf, p
dia2 <- dia[, .(rsid = rsID, chr = CHR, pos = as.integer(BP), ea = toupper(EA), oa = toupper(NEA), b = as.numeric(Beta), se = as.numeric(SE), eaf = as.numeric(EAF), p = as.numeric(Pvalue))]
if (MUT == "dia_drop_arg1") dia2 <- dia2[rsid != "rs2781668"]   # mutation: ARG1 instrument absent from DIAMANTE
if (MUT == "dia_flip_arg1") dia2[rsid == "rs2781668", b := -b]   # mutation: ARG1 check estimate sign flipped
fmt <- function(FORMAT, SAMPLE, key) mapply(function(f, s) { v <- strsplit(s, ":", fixed = TRUE)[[1]]; k <- strsplit(f, ":", fixed = TRUE)[[1]]; if (key %in% k) v[match(key, k)] else NA_character_ }, FORMAT, SAMPLE, USE.NAMES = FALSE)
car2 <- car[, .(rsid = ID, chr = CHROM, pos = as.integer(POS), ea = toupper(ALT), oa = toupper(REF), b = as.numeric(fmt(FORMAT, SAMPLE, "ES")),
                se = as.numeric(fmt(FORMAT, SAMPLE, "SE")), eaf = as.numeric(fmt(FORMAT, SAMPLE, "AF")), p = 10^-as.numeric(fmt(FORMAT, SAMPLE, "LP")))]
for (nm in c("dia2", "car2")) { d <- get(nm); stopifnot(!anyNA(d$b), !anyNA(d$se), !anyNA(d$p)); say(nm, "rows", nrow(d), "| distinct rsid", uniqueN(d$rsid), "| NA eaf", sum(is.na(d$eaf))) }

# ---- per-pair estimates (engine logic re-implemented, plus SE-order variants)
pal <- function(a, b) paste(sort(c(a, b)), collapse = "") %in% c("AT", "CG")
P <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  r <- pairs[i]; src <- if (r$outcome == "T2D") dia2 else car2
  y <- src[rsid == r$lead_snps]
  base <- data.table(gene = r$gene, outcome = r$outcome, tier = r$evidence_tier, prior_robust = r$prior_robust, instrument = r$lead_snps,
                     predeclared = r$predeclared, disc_b = r$disc_b, disc_se = r$disc_se, bx = r$bx, sx = r$sx, eaf_exp = r$eaf,
                     n_resource_rows = nrow(y), pos37_ukb = r$pos37)
  if (nrow(y) != 1L) return(cbind(base, cls = "not evaluable (instrument absent)"))
  s <- if (y$ea == r$effect_allele && y$oa == r$other_allele) 1 else if (y$ea == r$other_allele && y$oa == r$effect_allele) -1 else 0
  if (s == 0) return(cbind(base, cls = "not evaluable (alleles do not match)"))
  by <- s * y$b; eafy <- if (s == 1) y$eaf else 1 - y$eaf
  flag <- ""
  if (pal(r$effect_allele, r$other_allele)) { maf <- min(r$eaf, 1 - r$eaf)
    flag <- if (maf > 0.42) "palindromic, AF-uninformative" else if (abs(eafy - r$eaf) > 0.10) "palindromic, AF DISCORDANT" else "palindromic, AF concordant" }
  w <- by / r$bx; se1 <- y$se / abs(r$bx); se2 <- sqrt(y$se^2 / r$bx^2 + by^2 * r$sx^2 / r$bx^4)
  disc_se2 <- sqrt(r$disc_se^2 + r$disc_b^2 * r$sx^2 / r$bx^2)
  same <- (w > 0) == (r$disc_b > 0); p <- y$p
  cls <- if (same && p < 0.05) "concordant" else if (same) "directionally concordant" else if (p < 0.05) "discordant" else "inconclusive"
  hp <- function(sc, sd) 2 * pnorm(-abs((w - r$disc_b) / sqrt(sc^2 + sd^2)))
  cbind(base, cls = cls, same_sign = same, wald = w, se = se1, se_2nd = se2, p = p, p_wald_1st = 2 * pnorm(-abs(w / se1)), p_wald_2nd = 2 * pnorm(-abs(w / se2)),
        OR = exp(w), OR_lo = exp(w - 1.96 * se1), OR_hi = exp(w + 1.96 * se1), disc_se_2nd = disc_se2,
        het_p = hp(se1, r$disc_se), het_p_disc2 = hp(se1, disc_se2), het_p_check2 = hp(se2, r$disc_se), het_p_both2 = hp(se2, disc_se2),
        palindrome = flag, eaf_res = eafy, pos37_res = y$pos, pos_match = (y$pos == r$pos37), res_chr = y$chr)
}), fill = TRUE)
fwrite(P, file.path(OUT, "s8_pairs.tsv"), sep = "\t")
ev8 <- P[!startsWith(cls, "not evaluable")]
say("evaluable", nrow(ev8), "of", nrow(P), "| not evaluable:", paste(P[startsWith(cls, "not evaluable"), paste(gene, outcome, cls)], collapse = "; "))
say("position check (resource GRCh37 position vs UKB-PPP variant id GRCh37):", sum(ev8$pos_match), "of", nrow(ev8), "match")
for (i in seq_len(nrow(P))) say(sprintf("  %-8s %-4s %-10s pre %-5s robust %-5s cls %-26s w %s se %s p %s het %s %s",
  P$gene[i], P$outcome[i], P$tier[i], P$predeclared[i], P$prior_robust[i], P$cls[i], f10(P$wald[i]), f10(P$se[i]), f10(P$p[i]), f10(P$het_p[i]), P$palindrome[i]))
for (i in seq_len(nrow(P))) {
  r <- P[i]
  addr(paste0("S8.pair.", r$gene, ".", r$outcome), "s8; posthoc spec item 1", "Wald ratio beta_out/beta_exp at the discovery lead (no proxy), allele-label harmonised; class on the resource's own p",
       if (r$outcome == "T2D") "OUT_DIAMANTE2018_noUKB_EUR_rsid" else "OUT_CARDIOGRAMplusC4D_ieu-a-7", value = r$wald, se = r$se,
       ci_lo = r$wald - 1.96 * r$se, ci_hi = r$wald + 1.96 * r$se, p = r$p, n = 1, unit = "log-odds per SD of protein", se_order = "1st (SE_out/|bx|)",
       notes = sprintf("class %s; tier %s; predeclared %s; prior_robust %s; disc_b %s (se %s); het_p %s; p_wald_1st %s; p_wald_2nd %s; %s",
                       r$cls, r$tier, r$predeclared, r$prior_robust, f10(r$disc_b), f10(r$disc_se), f10(r$het_p), f10(r$p_wald_1st), f10(r$p_wald_2nd), r$palindrome))
}

# ---- discovery SE convention: verify against FinnGen R10 T2D (tabix, read in place)
say("--- discovery SE order check against FinnGen R10 T2D at the T2D leads (scanTabix on the existing .tbi; GRCh38 positions from the instrument table)")
p_fg <- inp("OUT_FG_R10_T2D"); tbi <- paste0(p_fg, ".tbi"); m_tbi <- md5(tbi)
say("INPUT OUT_FG_R10_T2D index", tbi, m_tbi, "| expected 7c7965b8e6794d20e897be42a8666895 (RESOURCES notes)")
stopifnot(identical(m_tbi, "7c7965b8e6794d20e897be42a8666895"))
suppressPackageStartupMessages(library(Rsamtools))
say("Rsamtools", as.character(packageVersion("Rsamtools")))
tf <- TabixFile(p_fg, index = tbi)
t2d <- pairs[outcome == "T2D"]
FG <- rbindlist(lapply(seq_len(nrow(t2d)), function(i) {
  r <- t2d[i]; ch <- sub(":.*", "", r$chrpos); ps <- as.integer(sub(".*:", "", r$chrpos))
  x <- scanTabix(tf, param = GenomicRanges::GRanges(ch, IRanges::IRanges(ps, ps)))[[1]]
  if (!length(x)) return(data.table(gene = r$gene, fg_rows = 0L))
  sp <- do.call(rbind, strsplit(x, "\t", fixed = TRUE))
  k <- which((sp[, 3] == r$other_allele & sp[, 4] == r$effect_allele) | (sp[, 3] == r$effect_allele & sp[, 4] == r$other_allele))
  if (length(k) != 1L) return(data.table(gene = r$gene, fg_rows = length(x), fg_allele_rows = length(k)))
  s <- if (sp[k, 4] == r$effect_allele) 1 else -1
  by <- s * as.numeric(sp[k, 9]); sy <- as.numeric(sp[k, 10])
  data.table(gene = r$gene, fg_rows = length(x), fg_allele_rows = 1L, b_fg = by / r$bx, se1_fg = sy / abs(r$bx),
             se2_fg = sqrt(sy^2 / r$bx^2 + by^2 * r$sx^2 / r$bx^4), disc_b = r$disc_b, disc_se = r$disc_se, any_proxy = r$any_proxy)
}), fill = TRUE)
FG[, b_match := abs(b_fg - disc_b) < 1e-6 * abs(disc_b)]
FG[, se_is_1st := abs(se1_fg - disc_se) < 1e-6 * disc_se]
FG[, se_is_2nd := abs(se2_fg - disc_se) < 1e-6 * disc_se]
fwrite(FG, file.path(OUT, "s8_discovery_se_order_check.tsv"), sep = "\t")
for (i in seq_len(nrow(FG))) say(sprintf("  %-8s rows %s b_match %s se_is_1st %s se_is_2nd %s proxy %s", FG$gene[i], FG$fg_rows[i], FG$b_match[i], FG$se_is_1st[i], FG$se_is_2nd[i], FG$any_proxy[i]))
nb <- FG[b_match == TRUE, .N]; n1 <- FG[b_match == TRUE & se_is_1st == TRUE, .N]; n2 <- FG[b_match == TRUE & se_is_2nd == TRUE, .N]
say("stored mr_all_results SE vs FinnGen at matching-b T2D pairs:", n1, "of", nb, "equal the FIRST-order SE;", n2, "equal the second-order SE")
DISC_SE_ORDER <- if (nb > 0 && n1 == nb) "1st" else if (nb > 0 && n2 == nb) "2nd" else "mixed"
addr("S8.disc_se_order_check", "s8 (SE convention of disc_se)", "SE order of the stored discovery Wald SE (mr_all_results 'se'), checked against FinnGen R10 T2D at T2D leads where b matches",
     "TAB_mr_all_results + OUT_FG_R10_T2D", value = DISC_SE_ORDER, n = nb, unit = "order",
     notes = sprintf("first-order equal at %d of %d; second-order equal at %d of %d; b mismatch (proxy or other) at %d of %d T2D pairs: %s; CAD rows not checked (I9_CHD has no local .tbi), the stored column is one pipeline column",
                     n1, nb, n2, nb, nrow(FG) - nb, nrow(FG), paste(FG[b_match != TRUE | is.na(b_match), gene], collapse = ",")))
if (DISC_SE_ORDER != "1st") say("WARNING: stored discovery SE is not uniformly first order; the 'disc 2nd' variants below assume it is first order")

# ---- summary statistics (spec items 1-4)
bt <- function(k, n) binom.test(k, n, 0.5)$p.value
summ <- function(d, lab) {
  k <- sum(d$same_sign); n <- nrow(d)
  list(lab = lab, n = n, same = k, binom_p = bt(k, n), conc = sum(d$cls == "concordant"), p05 = sum(d$p < 0.05),
       het = sum(d$het_p < 0.05), het_d2 = sum(d$het_p_disc2 < 0.05), het_c2 = sum(d$het_p_check2 < 0.05), het_b2 = sum(d$het_p_both2 < 0.05),
       conc_w1 = sum(d$same_sign & d$p_wald_1st < 0.05), conc_w2 = sum(d$same_sign & d$p_wald_2nd < 0.05),
       classes = paste(names(table(d$cls)), table(d$cls), sep = " ", collapse = "; "))
}
S_all <- summ(ev8, "ALL"); S_np <- summ(ev8[predeclared == FALSE], "NOT PRE-DECLARED")
for (S in list(S_all, S_np)) say(sprintf("%s: n %d | classes %s | same sign %d of %d (binomial p %s) | concordant %d | resource p<0.05 %d | het p<0.05: disc1/check1 %d, disc2/check1 %d, disc1/check2 %d, disc2/check2 %d | concordant by Wald p 1st %d, 2nd %d",
  S$lab, S$n, S$classes, S$same, S$n, f10(S$binom_p), S$conc, S$p05, S$het, S$het_d2, S$het_c2, S$het_b2, S$conc_w1, S$conc_w2))
nf <- ev8[tier != "Flagged"]; fl <- ev8[tier == "Flagged"]
tabT <- matrix(c(sum(nf$same_sign), nrow(nf) - sum(nf$same_sign), sum(fl$same_sign), nrow(fl) - sum(fl$same_sign)), 2, byrow = TRUE)
rb <- ev8[prior_robust == TRUE]; nr <- ev8[prior_robust == FALSE]
tabP <- matrix(c(sum(rb$same_sign), nrow(rb) - sum(rb$same_sign), sum(nr$same_sign), nrow(nr) - sum(nr$same_sign)), 2, byrow = TRUE)
fT <- fisher.test(tabT); fP <- fisher.test(tabP)
say("exploratory Fisher (not Flagged vs Flagged, same sign):", paste(tabT, collapse = ","), "(row-major: notFlagged same, opp; Flagged same, opp) p", f10(fT$p.value))
say("exploratory Fisher (prior-robust vs not, same sign):", paste(t(tabP), collapse = ","), "p", f10(fP$p.value))
addr("S8.binom.all.p", "s8; posthoc spec item 2", "exact binomial test (two-sided, H0 0.5) of same-sign count among evaluable pairs; descriptive (pairs at shared loci are not independent)",
     "OUT_DIAMANTE2018_noUKB_EUR_rsid + OUT_CARDIOGRAMplusC4D_ieu-a-7", value = S_all$binom_p, p = S_all$binom_p, n = S_all$n, unit = "p",
     notes = sprintf("%d of %d same sign; post hoc", S_all$same, S_all$n))
addr("S8.binom.notpredeclared.p", "s8; posthoc spec item 2 + erratum (23 not 22)", "exact binomial test of same-sign count among evaluable pairs not pre-declared",
     "OUT_DIAMANTE2018_noUKB_EUR_rsid + OUT_CARDIOGRAMplusC4D_ieu-a-7", value = S_np$binom_p, p = S_np$binom_p, n = S_np$n, unit = "p",
     notes = sprintf("%d of %d same sign; post hoc", S_np$same, S_np$n))
addr("S8.fisher.tier.p", "s8; posthoc spec item 3", "Fisher exact test, same sign: not Flagged (Convergent+Supported) vs Flagged; exploratory", "derived",
     value = fT$p.value, p = fT$p.value, n = nrow(ev8), unit = "p", notes = sprintf("table [[%d,%d],[%d,%d]]; OR %s", tabT[1, 1], tabT[1, 2], tabT[2, 1], tabT[2, 2], f10(fT$estimate)))
addr("S8.fisher.prior.p", "s8; posthoc spec item 3", "Fisher exact test, same sign: prior-robust (PP.H4 at p12 1e-6 >= 0.80) vs not; exploratory", "derived",
     value = fP$p.value, p = fP$p.value, n = nrow(ev8), unit = "p", notes = sprintf("table [[%d,%d],[%d,%d]]; OR %s", tabP[1, 1], tabP[1, 2], tabP[2, 1], tabP[2, 2], f10(fP$estimate)))
for (t in c("Convergent", "Supported", "Flagged")) { d <- ev8[tier == t]
  addr(paste0("S8.tier.", t, ".same"), "s8; posthoc spec item 3", paste("same-sign pairs, tier", t), "derived", value = sum(d$same_sign), n = nrow(d), unit = "pairs")
  addr(paste0("S8.tier.", t, ".concordant"), "s8; posthoc spec item 3", paste("concordant pairs (same sign, resource p<0.05), tier", t), "derived", value = sum(d$cls == "concordant"), n = nrow(d), unit = "pairs") }
for (S in list(S_all, S_np)) {
  tg <- if (S$lab == "ALL") "all" else "notpredeclared"
  addr(paste0("S8.", tg, ".n_evaluable"), "s8", paste("evaluable pairs,", S$lab), "derived", value = S$n, unit = "pairs", notes = S$classes)
  addr(paste0("S8.", tg, ".same_sign"), "s8", paste("pairs with the discovery sign,", S$lab), "derived", value = S$same, n = S$n, unit = "pairs")
  addr(paste0("S8.", tg, ".concordant"), "s8", paste("concordant pairs (same sign, resource p<0.05),", S$lab), "derived", value = S$conc, n = S$n, unit = "pairs")
  addr(paste0("S8.", tg, ".p_lt_005"), "s8", paste("pairs with resource p<0.05 (any sign),", S$lab), "derived", value = S$p05, n = S$n, unit = "pairs")
  addr(paste0("S8.", tg, ".het.disc1_check1.n_lt_005"), "s8; posthoc spec item 4", paste("discovery-vs-check heterogeneity z-test p<0.05 (disc SE as stored =", DISC_SE_ORDER, "order; check SE 1st order) -- the engine's convention,", S$lab),
       "derived", value = S$het, n = S$n, unit = "pairs", se_order = paste0("disc ", DISC_SE_ORDER, " (stored) / check 1st"))
  addr(paste0("S8.", tg, ".het.disc2_check1.n_lt_005"), "s8", paste("heterogeneity p<0.05 with discovery SE second order (stored SE + b^2 sx^2/bx^2), check SE 1st order,", S$lab),
       "derived", value = S$het_d2, n = S$n, unit = "pairs", se_order = "disc 2nd / check 1st")
  addr(paste0("S8.", tg, ".het.disc1_check2.n_lt_005"), "s8", paste("heterogeneity p<0.05 with discovery SE as stored, check SE second order,", S$lab),
       "derived", value = S$het_c2, n = S$n, unit = "pairs", se_order = "disc stored / check 2nd")
  addr(paste0("S8.", tg, ".het.disc2_check2.n_lt_005"), "s8", paste("heterogeneity p<0.05 with both SEs second order,", S$lab),
       "derived", value = S$het_b2, n = S$n, unit = "pairs", se_order = "disc 2nd / check 2nd")
  addr(paste0("S8.", tg, ".concordant_by_waldp_1st"), "s8", paste("same sign and Wald p (1st-order SE) < 0.05,", S$lab), "derived", value = S$conc_w1, n = S$n, unit = "pairs", se_order = "check 1st")
  addr(paste0("S8.", tg, ".concordant_by_waldp_2nd"), "s8", paste("same sign and Wald p (2nd-order SE) < 0.05,", S$lab), "derived", value = S$conc_w2, n = S$n, unit = "pairs", se_order = "check 2nd")
}
addr("S8.n_pairs", "s8 Also derived", "colocalization-supported T2D/CAD pairs in the tiered table", "TAB_tiered81", value = nrow(P), unit = "pairs",
     notes = sprintf("T2D %d, CAD %d; tiers %s", P[outcome == "T2D", .N], P[outcome == "CAD", .N], paste(names(table(P$tier)), table(P$tier), collapse = ", ")))
addr("S8.n_predeclared", "s8 Also derived", "pairs among the 28 pre-declared in the outcome-check pre-specification", "outcome-check prespec + TAB_tiered81", value = sum(P$predeclared), n = nrow(P), unit = "pairs")
addr("S8.n_notpredeclared", "s8 Also derived", "pairs among the 28 not pre-declared", "derived", value = sum(!P$predeclared), n = nrow(P), unit = "pairs")
addr("S8.n_leads", "s8 Also derived", "distinct lead (instrument) variants among the 28 pairs", "TAB_tiered81", value = uniqueN(P$instrument), n = nrow(P), unit = "variants")

# assertion targets visible in the redacted section 8 'Also derived' list (rule 9: logged, never stop the run)
assert_target("S8.all28_predeclared", sum(P$predeclared), 5, "pre-declared pairs among the 28 (outcome-check prespec s1, European five)", "s8 Also derived")
assert_target("S8.all28_notpredeclared", sum(!P$predeclared), 23, "28 minus pre-declared", "s8 Also derived")
assert_target("S8.all28_same_sign", S_all$same, 27, "same-sign pairs among the evaluable (all 28 if all evaluable)", "s8 Also derived")
assert_target("S8.all28_p_lt_005", S_all$p05, 20, "pairs with resource p<0.05 among evaluable (any sign)", "s8 Also derived")
assert_target("S8.all28_concordant", S_all$conc, 20, "concordant pairs (same sign and resource p<0.05) among evaluable -- second reading of '20 at p<0.05'", "s8 Also derived")
assert_target("S8.all28_leads", uniqueN(P$instrument), 27, "distinct lead variants among the 28", "s8 Also derived")
say("NOTE '16 -> 12 dependence on the SE convention': the redacted text does not define the quantity. Candidates computed before reading the full text:",
    sprintf("het count ALL disc1/check1 %d, disc2/check1 %d, disc1/check2 %d, disc2/check2 %d; NOT-PREDECLARED %d, %d, %d, %d; concordant by Wald p 1st %d / 2nd %d (ALL), %d / %d (23)",
            S_all$het, S_all$het_d2, S_all$het_c2, S_all$het_b2, S_np$het, S_np$het_d2, S_np$het_c2, S_np$het_b2, S_all$conc_w1, S_all$conc_w2, S_np$conc_w1, S_np$conc_w2))
cand <- list(c("het disc1/check1 -> disc2/check1 (ALL)", S_all$het, S_all$het_d2), c("het disc1/check1 -> disc1/check2 (ALL)", S_all$het, S_all$het_c2),
             c("het disc1/check1 -> disc2/check2 (ALL)", S_all$het, S_all$het_b2), c("concordant Wald p 1st -> 2nd (ALL)", S_all$conc_w1, S_all$conc_w2),
             c("het disc1/check1 -> disc2/check1 (23)", S_np$het, S_np$het_d2), c("het disc1/check1 -> disc1/check2 (23)", S_np$het, S_np$het_c2),
             c("het disc1/check1 -> disc2/check2 (23)", S_np$het, S_np$het_b2), c("concordant Wald p 1st -> 2nd (23)", S_np$conc_w1, S_np$conc_w2))
# candidate round 2 (added after round 1 matched nothing, still BEFORE the full text was opened): discovery-side significance of the
# 28 pairs under the stored (first-order) vs second-order discovery SE, against the family Bonferroni threshold 0.05/15,584 (Rule 1)
BONF <- 0.05 / 15584
P[, disc_p1 := 2 * pnorm(-abs(disc_b / disc_se))]; P[, disc_p2 := 2 * pnorm(-abs(disc_b / disc_se_2nd))]
for (sub_lab in c("ALL", "23", "T2D", "CAD")) {
  d <- switch(sub_lab, ALL = P, `23` = P[predeclared == FALSE], T2D = P[outcome == "T2D"], CAD = P[outcome == "CAD"])
  cand[[length(cand) + 1L]] <- c(sprintf("discovery Bonferroni (p<0.05/15584) stored 1st -> 2nd order SE (%s, n %d)", sub_lab, nrow(d)),
                                 sum(d$disc_p1 < BONF), sum(d$disc_p2 < BONF))
  addr(paste0("S8.disc_bonf.", sub_lab), "s8 (16 -> 12 candidate, round 2)", paste("pairs with discovery p < 0.05/15584 under the stored first-order SE; second-order count in notes;", sub_lab),
       "TAB_mr_all_results + TAB_instruments", value = sum(d$disc_p1 < BONF), n = nrow(d), unit = "pairs", se_order = "disc 1st (stored)",
       notes = sprintf("second-order disc SE (stored SE^2 + b^2 sx^2/bx^2): %d; movers: %s", sum(d$disc_p2 < BONF), paste(d[disc_p1 < BONF & disc_p2 >= BONF, pl(gene, "-", outcome)], collapse = ", ")))
}
say("16 -> 12 candidate round 2 (discovery Bonferroni by SE order):", paste(vapply(cand[-(1:8)], function(z) paste0(z[1], " = ", z[2], " -> ", z[3]), ""), collapse = "; "))
hit <- Filter(function(z) z[2] == "16" && z[3] == "12", cand)
addr("ASSERT.S8.se_dependence_16_to_12", "s8 Also derived", "the redacted text prints '16 -> 12 dependence on the SE convention' without defining the count; candidate counts computed first",
     "derived", value = if (length(hit)) paste(vapply(hit, `[`, "", 1), collapse = " | ") else "no candidate gives 16 -> 12",
     notes = paste0(if (length(hit)) "ASSERT_PASSED (definition inferred, weak evidence): " else "ASSERT_FAILED/UNDETERMINED: ",
                    paste(vapply(cand, function(z) paste0(z[1], " = ", z[2], " -> ", z[3]), ""), collapse = "; ")))
if (!length(hit)) ASSERTS <- c(ASSERTS, "S8.se_dependence_16_to_12: no candidate definition reproduces 16 -> 12")
say("16 -> 12 candidates matching:", if (length(hit)) paste(vapply(hit, `[`, "", 1), collapse = " | ") else "none")
addb("s8 (printed numbers; no wording branch)", sprintf("binomial p all %s (%d/%d); binomial p 23 %s (%d/%d); Fisher tier p %s; Fisher prior p %s",
       f10(S_all$binom_p), S_all$same, S_all$n, f10(S_np$binom_p), S_np$same, S_np$n, f10(fT$p.value), f10(fP$p.value)), "PRINTED",
     sprintf("s8 fill: not pre-declared %d evaluable, %d same sign (binomial p %s), %d concordant; all %d evaluable, %d same sign (binomial p %s); Fisher tier p %s; Fisher prior-band p %s",
             S_np$n, S_np$same, g3(S_np$binom_p), S_np$conc, S_all$n, S_all$same, g3(S_all$binom_p),
             g3(fT$p.value), g3(fP$p.value)),
     "post hoc; UKB-free external consortium check, not independent replication (FinnGen overlap); descriptive binomial (shared loci not independent); tier/prior Fisher exploratory")

# =================================================================== write outputs, end-of-run checks
R <- rbindlist(RESULTS); B <- rbindlist(BRANCHES)
stopifnot(!anyDuplicated(R$quantity_id))
fwrite(R, file.path(OUT, "results.tsv"), sep = "\t", quote = FALSE)
fwrite(B, file.path(OUT, "branches.tsv"), sep = "\t", quote = FALSE)
say("results.tsv rows", nrow(R), "| branches.tsv rows", nrow(B), "| ASSERT_FAILED:", length(ASSERTS), if (length(ASSERTS)) paste(ASSERTS, collapse = " ; ") else "")
SNAP1 <- snap("end")
same_snap <- identical(SNAP0$file, SNAP1$file) && identical(SNAP0$md5, SNAP1$md5)
say("STAGING unchanged start vs end:", same_snap)
POSTHOC_MD5_1 <- vapply(POSTHOC_IDS, function(i) md5(res_get(i)$path), "")
say("POSTHOC seen outputs unchanged at exit:", identical(POSTHOC_MD5_0, POSTHOC_MD5_1), paste(POSTHOC_IDS, POSTHOC_MD5_1, collapse = "; "))
stopifnot(same_snap, identical(POSTHOC_MD5_0, POSTHOC_MD5_1))
for (f in c("results.tsv", "branches.tsv", "s7_scale_ckb.tsv", "s7_scale_eas_egfr_derived.tsv", "s7_strat_table.tsv", "s7_ckbfree_pairs.tsv", "s7_ckbfree_composite.tsv",
            "s8_pairs.tsv", "s8_discovery_se_order_check.tsv", "extracts/DIAMANTE_w7.tsv", "extracts/CARDIOGRAM_w7.tsv"))
  say("OUTPUT md5", md5(file.path(OUT, f)), f)
say("W7 done in", round(as.numeric(difftime(Sys.time(), T0, units = "secs"))), "s")
