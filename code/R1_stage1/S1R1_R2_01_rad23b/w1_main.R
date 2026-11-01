# W1 -- Stage 1, pre-declaration section 1: RAD23B, the five proxy-marked pairs, the SE convention.
# One command reruns everything:  Rscript w1_main.R   (extraction checkpoints in extracts/ are reused when valid)
# Computed from the REDACTED pre-declaration only. COMPARE_SEEN.tsv is written afterwards by w1_compare_seen.R.
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
t_start <- Sys.time()
OUTDIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_01_rad23b")
Sys.setenv(W1_LOG = file.path(OUTDIR, "run.log"))
source(file.path(OUTDIR, "w1_common.R"))
cat("", file = LOG)
SEED <- 20261007L; set.seed(SEED)
logf("W1 S1R1_R2_01_rad23b start. Pre-declaration section 1 (proxy ledger, direct-lead, exposure orientation, SE convention, RAD23B external check).")
logf("seed %d set (W1 draws no random numbers: no bootstrap is declared for section 1)", SEED)
logf("seen scripts opened before COMPARE_SEEN: NONE. Seen outputs read: NONE.")
assert_governing(); load_resources()
snap0 <- staging_snapshot("start")

# ---------------------------------------------------------------- tools
pkgv <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "absent")
logf("R: %s | data.table %s | jsonlite %s | digest %s | Rsamtools %s | rtracklayer %s | BEDMatrix %s",
     R.version.string, pkgv("data.table"), pkgv("jsonlite"), pkgv("digest"), pkgv("Rsamtools"), pkgv("rtracklayer"), pkgv("BEDMatrix"))
assert_md5(RSCRIPT, res_entry("R_4.4.1")$md5, "TOOL:Rscript")
PLINK <- rpath("PLINK19_paper6")
pv <- system2(PLINK, "--version", stdout = TRUE, stderr = TRUE); logf("PLINK --version: %s", paste(pv, collapse = " "))
st <- suppressWarnings(system2(RSCRIPT, c(qa(GOV$gls_q[1]), "--self-test"), stdout = TRUE, stderr = TRUE))
for (l in st) logf("gls_q.R --self-test: %s", l)
if (!any(grepl("SELF-TEST PASS", st))) stop("gls_q.R self-test did not pass")
logf("gls_q.R is not called by W1 (section 1 has no GLS Q); md5 asserted and self-test logged per rule 6.")

# ---------------------------------------------------------------- result helpers
RESULTS <- list(); BR <- list(); AFAIL <- character()
fmtv <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.logical(x)) return(ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE")))
  if (is.numeric(x)) return(ifelse(is.na(x), "NA", trimws(formatC(x, digits = 8, format = "g"))))
  gsub("[\t\r\n]+", " ", as.character(x))
}
addr <- function(id, ref, estimand, resource = "", panel = "", cov = "", se_order = "", value = "", se = NA, ci_lo = NA,
                 ci_hi = NA, p = NA, n = NA, unit = "", notes = "") {
  RESULTS[[length(RESULTS) + 1L]] <<- data.table(quantity_id = id, prespec_ref = ref, estimand = estimand, resource = resource,
    panel = panel, cov = cov, se_order = se_order, value = fmtv(value), se = fmtv(se), ci_lo = fmtv(ci_lo), ci_hi = fmtv(ci_hi),
    p = fmtv(p), n = fmtv(n), unit = unit, notes = fmtv(notes))
}
assert_target <- function(id, ref, computed, target, derivation) {
  ok <- identical(as.character(computed), as.character(target))
  logf("ASSERT %s: computed %s | target %s | %s | derivation: %s", id, paste(computed, collapse = ","), paste(target, collapse = ","),
       if (ok) "PASS" else "ASSERT_FAILED", derivation)
  addr(paste0("S1.assert.", id), ref, "pre-declared assertion target (computed first, then compared)",
       value = if (ok) "PASS" else "ASSERT_FAILED",
       notes = sprintf("computed %s; target %s; derivation: %s", paste(computed, collapse = ","), paste(target, collapse = ","), derivation))
  if (!ok) AFAIL <<- c(AFAIL, sprintf("%s: computed %s vs target %s", id, paste(computed, collapse = ","), paste(target, collapse = ",")))
  invisible(ok)
}
addb <- function(rule_ref, inputs, branch, filled, qualifiers) BR[[length(BR) + 1L]] <<- data.table(rule_ref = rule_ref,
  inputs = fmtv(inputs), branch = fmtv(branch), filled_template = fmtv(filled), qualifiers = fmtv(qualifiers))
p2s <- function(z) 2 * pnorm(-abs(z))
null_rates <- function(d, cols, label) for (cc in cols) logf("null rate %s.%s = %.4f (n=%d)", label, cc, mean(is.na(d[[cc]]) | d[[cc]] %in% c("", "NA")), nrow(d))
COMP <- c(A = "T", T = "A", C = "G", G = "C")
is_pal <- function(a, b) toupper(a) %in% names(COMP) && unname(COMP[toupper(a)]) == toupper(b)

# ---------------------------------------------------------------- inputs (small tables)
mr   <- fread(file = rpath("TAB_mr_all_results"));  if (nrow(mr) != 14537L) stop("mr_all_results rows ", nrow(mr))
ins  <- fread(file = rpath("TAB_instruments"), colClasses = list(character = c("chrpos", "ukb_variant_id", "effect_allele", "other_allele")))
if (nrow(ins) != 1964L) stop("instrument rows ", nrow(ins))
tier <- fread(file = rpath("TAB_tiered81"));        if (nrow(tier) != 81L) stop("tiered81 rows")
fig4 <- fread(file = rpath("FIG04_panel1_2"));      if (nrow(fig4) != 81L) stop("fig4 rows")
dec  <- fread(file = rpath("TAB_decode_mr_results")); if (nrow(dec) < 1L) stop("decode rows")
log28 <- readLines(rpath("LOG_28_candidate_pool"))
man  <- fread(file = rpath("OUT_FG_R10_MANIFEST"), sep = "\t")
null_rates(mr, c("gene", "outcome_id", "method", "b", "se", "pval"), "mr_all")
null_rates(ins, c("gene", "SNP", "effect_allele", "other_allele", "beta", "se", "eaf", "F_stat"), "instruments")
null_rates(tier, c("gene", "outcome", "mr_pval", "evidence_tier", "crit_bonferroni", "crit_prior_stable", "crit_susie_supportive", "has_major_weakness", "n_evaluable"), "tiered81")

# N_eff (computed from case/control counts; compared with RESOURCES n_info; logged only -- no W1 estimate uses N)
neff <- function(ca, co) 4 * ca * co / (ca + co)
for (k in list(c("OUT_FG_R10_T2D", "T2D"), c("OUT_FG_R10_I9_AF", "I9_AF"), c("OUT_FG_R10_I9_CHD", "I9_CHD"))) {
  r <- man[phenocode == k[2]]; stopifnot(nrow(r) == 1L)
  ne <- neff(r$num_cases, r$num_controls); ni <- res_entry(k[1])$n_info
  stated <- as.numeric(sub(".*= ([0-9]+) \\(computed\\).*", "\\1", ni))
  logf("N_eff %s: cases %d controls %d -> %.0f (manifest); n_info states %s; diff %.1f", k[1], r$num_cases, r$num_controls, ne, stated, ne - stated)
}
for (k in c("OUT_MEGASTROKE_EUR_GCST006906", "OUT_DIAMANTE2018_noUKB_EUR_chrpos", "OUT_CARDIOGRAMplusC4D_ieu-a-7")) {
  ni <- res_entry(k)$n_info; ca <- as.numeric(sub("^cases ([0-9]+) .*", "\\1", ni)); co <- as.numeric(sub("^cases [0-9]+ / controls ([0-9]+).*", "\\1", ni))
  stated <- as.numeric(sub(".*= ([0-9]+) \\(computed\\).*", "\\1", ni))
  logf("N_eff %s: cases %.0f controls %.0f (from n_info; file has no N column) -> %.0f; n_info states %s; diff %.1f", k, ca, co, neff(ca, co), stated, neff(ca, co) - stated)
}

# targets (lead/proxy positions) and extraction checkpoints
tg <- if (file.exists(file.path(OUT, "extracts", "targets.tsv"))) fread(file = file.path(OUT, "extracts", "targets.tsv"), colClasses = list(character = c("chr", "chr_bim"))) else build_targets()
tg_chk <- build_targets()   # always rebuilt and compared, so a stale targets file cannot drive the extracts
if (!isTRUE(all.equal(as.data.frame(tg[, .(bim_id, pos37, pos38, A1, A2)]), as.data.frame(tg_chk[, .(bim_id, pos37, pos38, A1, A2)]), check.attributes = FALSE))) stop("targets.tsv is stale")
tg <- tg_chk
jobs <- extraction_jobs(tg)
for (j in names(jobs)) run_job(j, jobs)
for (j in names(jobs)) { m <- jsonlite::fromJSON(file.path(OUT, "extracts", paste0(j, ".tsv.meta.json"))); logf("extract %s: %s data rows streamed in %.0f s; %d kept; complete=%s; out md5 %s", j, format(m$n_data_rows, big.mark = ","), m$seconds, m$n_kept, m$complete, m$out_md5) }
for (lf in sort(list.files(file.path(OUT, "extracts"), pattern = "^extract_.*[.]log$", full.names = TRUE))) {   # the parallel pre-fill runs, copied into run.log
  for (l in readLines(lf)) if (grepl("DONE|streaming|targets:|staging snapshot", l)) logf("[%s] %s", basename(lf), l) }
# completeness of the streams: rows against the resource map where it records them; FinnGen files end on chromosome 23
chk_rows <- function(j, expect) { m <- jsonlite::fromJSON(file.path(OUT, "extracts", paste0(j, ".tsv.meta.json"))); ok <- m$n_data_rows == expect
  logf("stream completeness %s: %s rows vs %s recorded in RESOURCES -> %s", j, format(m$n_data_rows, big.mark = ","), format(expect, big.mark = ","), if (ok) "EQUAL" else "DIFFERENT"); ok }
stopifnot(chk_rows("DIAM", 21508698), chk_rows("DIAMRS", 21508158), chk_rows("CARDIO", 8479693), chk_rows("BIM503", 22665064))
for (j in c("FG_AF", "FG_CHD")) { m <- jsonlite::fromJSON(file.path(OUT, "extracts", paste0(j, ".tsv.meta.json"))); if (!startsWith(m$last_line, "23\t")) stop(j, " stream did not reach chromosome 23"); logf("stream completeness %s: last line on chromosome 23 (%s rows)", j, format(m$n_data_rows, big.mark = ",")) }

# FinnGen T2D through its tabix index (read-only), at the RAD23B and ABO leads and proxies
suppressPackageStartupMessages(library(Rsamtools))
fg_t2d <- rpath("OUT_FG_R10_T2D"); xpath(paste0(fg_t2d, ".tbi"), "7c7965b8e6794d20e897be42a8666895", "FG_T2D_tbi")
tf <- TabixFile(fg_t2d, index = paste0(fg_t2d, ".tbi"))
fg_hdr <- strsplit(headerTabix(tf)$header, "\t")[[1]]
reg <- tg[gene %in% c("RAD23B", "ABO") & role == "lead", .(chr, lo = pos38 - 100000L, hi = pos38 + 100000L)]
tx <- unlist(lapply(seq_len(nrow(reg)), function(i) scanTabix(tf, param = GenomicRanges::GRanges(sprintf("%s:%d-%d", reg$chr[i], reg$lo[i], reg$hi[i])))[[1]]))
fgT2D <- fread(text = c(paste(fg_hdr, collapse = "\t"), tx), sep = "\t", colClasses = "character")
logf("FinnGen T2D tabix: %d rows in the two +/-100 kb windows (RAD23B, ABO leads)", nrow(fgT2D))
fwrite(fgT2D[pos %in% as.character(na.omit(tg$pos38)) | grepl("rs576123|rs505922|rs7035725|rs7856112", rsids)], file.path(OUT, "extracts", "FG_T2D_targets.tsv"), sep = "\t")

ex <- list(FG_T2D = fgT2D, FG_AF = read_extract("FG_AF"), FG_CHD = read_extract("FG_CHD"), MEGA = read_extract("MEGA"),
           DIAM = read_extract("DIAM"), DIAMRS = read_extract("DIAMRS"), CARDIO = read_extract("CARDIO"), BIM503 = read_extract("BIM503"))
for (k in c("FG_T2D", "FG_AF", "FG_CHD")) setnames(ex[[k]], "#chrom", "chrom", skip_absent = TRUE)

# ---------------------------------------------------------------- outcome-row lookup and harmonisation
# returns the outcome row for a variant: by GRCh38 position + allele set (forward or complement), else by rsID + allele set
out_row <- function(tag, chr, pos38, rsid, alleles) {
  d <- ex[[tag]]; al <- toupper(alleles); alc <- unname(COMP[al])
  if (startsWith(tag, "FG")) {
    d <- d[chrom == chr]
    p38c <- if (is.na(pos38)) "__none__" else as.character(pos38)
    cand <- d[pos == p38c | grepl(paste0("(^|,)", rsid, "(,|$)"), rsids)]
    cand <- cand[(toupper(ref) %in% al & toupper(alt) %in% al) | (toupper(ref) %in% alc & toupper(alt) %in% alc)]
    if (!nrow(cand)) return(NULL)
    if (nrow(cand) > 1L) stop("ambiguous outcome rows for ", rsid, " in ", tag)
    how <- if (!is.na(pos38) && cand$pos == as.character(pos38)) "GRCh38 position + alleles" else "rsID + alleles (position from the file)"
    return(list(tag = tag, ea = toupper(cand$alt), oa = toupper(cand$ref), beta = as.numeric(cand$beta), se = as.numeric(cand$sebeta),
                p = as.numeric(cand$pval), p_str = cand$pval, eaf = as.numeric(cand$af_alt), rsid_file = cand$rsids,
                pos_file = as.integer(cand$pos), how = how))
  }
  if (tag == "MEGA") {
    p38c <- if (is.na(pos38)) "__none__" else as.character(pos38)
    cand <- d[(hm_chrom == chr & hm_pos == p38c) | hm_rsid == rsid | variant_id == rsid]
    cand <- cand[toupper(hm_effect_allele) %in% c(al, alc) & toupper(hm_other_allele) %in% c(al, alc)]
    if (!nrow(cand)) return(NULL)
    if (nrow(cand) > 1L) stop("ambiguous MEGASTROKE rows for ", rsid)
    how <- if (!is.na(pos38) && cand$hm_pos == as.character(pos38)) "GRCh38 position + alleles" else "rsID + alleles (position from the file)"
    return(list(tag = tag, ea = toupper(cand$hm_effect_allele), oa = toupper(cand$hm_other_allele), beta = as.numeric(cand$hm_beta),
                se = as.numeric(cand$standard_error), p = as.numeric(cand$p_value), p_str = cand$p_value,
                eaf = as.numeric(cand$hm_effect_allele_frequency), rsid_file = cand$hm_rsid, pos_file = as.integer(cand$hm_pos), how = how))
  }
  stop("out_row: unknown tag ", tag)
}
# align an outcome row to target effect allele X (other allele Y). Palindromic pairs: by label only (forward-strand assumption).
align_row <- function(r, X, Y) {
  X <- toupper(X); Y <- toupper(Y)
  if (r$ea == X && r$oa == Y) return(list(beta = r$beta, eaf = r$eaf, how = "effect allele = target (by label)"))
  if (r$ea == Y && r$oa == X) return(list(beta = -r$beta, eaf = 1 - r$eaf, how = "swapped to target (by label)"))
  if (is_pal(X, Y)) stop("palindromic pair not matched by label")
  if (r$ea == COMP[[X]] && r$oa == COMP[[Y]]) return(list(beta = r$beta, eaf = r$eaf, how = "complement strand, effect allele = target"))
  if (r$ea == COMP[[Y]] && r$oa == COMP[[X]]) return(list(beta = -r$beta, eaf = 1 - r$eaf, how = "complement strand, swapped to target"))
  stop("alleles do not match")
}

# ---------------------------------------------------------------- the five proxy-marked pairs: identity checks
fig4p <- fig4[proxy_used == TRUE]
assert_target("proxy_pairs_n", "s1 (five proxy-marked pairs)", nrow(fig4p), 5L, "FIG04_panel1_2 rows with proxy_used == TRUE")
got_pairs <- sort(paste(fig4p$gene, fig4p$outcome_name, sep = "_"))
assert_target("proxy_pairs_set", "s1 (five proxy-marked pairs)", paste(got_pairs, collapse = ";"), paste(sort(PAIRS$pair), collapse = ";"),
              "FIG04_panel1_2 proxy_used == TRUE, gene_outcome; compared with the five pairs the pre-declaration names")
jy <- merge(PAIRS[, .(pair, gene, outcome_id)], mr, by = c("gene", "outcome_id"))
logf("join yield: 5 declared pairs -> mr_all_results rows: %d/5 (methods: %s)", nrow(jy), paste(unique(jy$method), collapse = ","))
stopifnot(nrow(jy) == 5L, all(jy$method == "Wald ratio"), all(jy$any_proxy == TRUE))
jt <- merge(PAIRS[, .(pair, gene, outcome = outcome_name)], tier, by = c("gene", "outcome"))
logf("join yield: 5 declared pairs -> tiered81: %d/5", nrow(jt)); stopifnot(nrow(jt) == 5L)

# ---------------------------------------------------------------- PLINK r (503 panel), square and table modes; mutation tests
bfile <- sub("\\.bed$", "", rpath("LD503_bed")); invisible(rpath("LD503_bim")); invisible(rpath("LD503_fam"))
fam_n <- length(readLines(rpath("LD503_fam"))); stopifnot(fam_n == 503L); logf("503 panel: %d samples in .fam", fam_n)
read_square <- function(prefix, ids_req) {
  f <- paste0(prefix, ".ld"); if (!file.exists(f)) stop("PLINK square: no .ld for ", basename(prefix))
  M <- as.matrix(fread(file = f, header = FALSE)); dimnames(M) <- NULL
  if (!all(dim(M) == length(ids_req))) stop(sprintf("PLINK square %s: %dx%d matrix for %d requested ids -- a requested SNP was not returned", basename(prefix), nrow(M), ncol(M), length(ids_req)))
  if (any(!is.finite(M))) stop(sprintf("PLINK square %s: non-finite cell", basename(prefix)))
  if (any(abs(diag(M) - 1) > 1e-6)) stop("PLINK square: diagonal not 1")
  M
}
plink_square <- function(gene, ids, tag) {
  idf <- file.path(OUT, "plink", paste0(tag, "_ids.txt")); writeLines(ids, idf)
  r <- run_plink(PLINK, c("--bfile", qa(bfile), "--keep-allele-order", "--extract", qa(idf), "--r", "square", "--out", qa(file.path(OUT, "plink", tag))), tag)
  read_square(file.path(OUT, "plink", tag), ids)
}
LD <- list()
for (g in unique(PAIRS$gene)) {
  L <- tg[gene == g & role == "lead"]; P <- tg[gene == g & role == "proxy"]
  ord <- if (L$pos37 < P$pos37) c("lead", "proxy") else c("proxy", "lead")   # PLINK writes the matrix in .bim order
  M <- plink_square(g, c(L$bim_id, P$bim_id), paste0("sq_", g))
  r <- M[1, 2]
  LD[[g]] <- list(r_A1A1 = r, r2 = r^2, order = ord, lead_A1 = L$A1, lead_A2 = L$A2, prox_A1 = P$A1, prox_A2 = P$A2)
  logf("PLINK --r square %s: r(lead %s A1=%s, proxy %s A1=%s) = %.6f, r2 = %.6f (order in .ld: %s)", g, L$bim_id, L$A1, P$bim_id, P$A1, r, r^2, paste(ord, collapse = ","))
}
# mutation test: a mutated SNP id must trip the 'every requested cell returned' assertion
mut <- tryCatch({ plink_square("RAD23B", c(tg[gene == "RAD23B" & role == "lead"]$bim_id, "rs7856112_MUTATED"), "sq_MUTATION_RAD23B"); "NOT TRIPPED" },
                error = function(e) paste("TRIPPED:", conditionMessage(e)))
logf("MUTATION TEST (square): %s", mut); if (!startsWith(mut, "TRIPPED")) stop("LD assertion is vacuous: mutated id did not trip it")
# table mode (second representation): r with in-phase alleles
all_ids <- unique(tg[role != "ext_instrument"]$bim_id); lead_ids <- tg[role == "lead"]$bim_id
writeLines(all_ids, file.path(OUT, "plink", "tab_ids.txt")); writeLines(lead_ids, file.path(OUT, "plink", "tab_leads.txt"))
run_plink(PLINK, c("--bfile", qa(bfile), "--keep-allele-order", "--extract", qa(file.path(OUT, "plink", "tab_ids.txt")), "--r", "in-phase", "with-freqs",
                   "--ld-snp-list", qa(file.path(OUT, "plink", "tab_leads.txt")), "--ld-window-r2", "0", "--ld-window", "99999", "--ld-window-kb", "99999",
                   "--out", qa(file.path(OUT, "plink", "tab"))), "tab_in_phase")
read_tab <- function(f, need) {
  tb <- fread(file = f, header = TRUE)
  for (k in seq_len(nrow(need))) {
    hit <- tb[(SNP_A == need$a[k] & SNP_B == need$b[k]) | (SNP_A == need$b[k] & SNP_B == need$a[k])]
    if (nrow(hit) != 1L || is.na(hit$R)) stop(sprintf("PLINK table: pair %s-%s not returned (or NA)", need$a[k], need$b[k]))
  }
  tb
}
need <- tg[role == "lead", .(gene, a = bim_id)][tg[role == "proxy", .(gene, b = bim_id)], on = "gene"]
TAB <- read_tab(file.path(OUT, "plink", "tab.ld"), need)
mut2 <- tryCatch({ read_tab(file.path(OUT, "plink", "tab.ld"), rbind(need, data.table(gene = "X", a = "rs7035725", b = "rs7856112_MUTATED"))); "NOT TRIPPED" },
                 error = function(e) paste("TRIPPED:", conditionMessage(e)))
logf("MUTATION TEST (table): %s", mut2); if (!startsWith(mut2, "TRIPPED")) stop("table-mode assertion is vacuous")
run_plink(PLINK, c("--bfile", qa(bfile), "--keep-allele-order", "--extract", qa(file.path(OUT, "plink", "tab_ids.txt")), "--freq", "--out", qa(file.path(OUT, "plink", "frq"))), "freq")
FRQ <- fread(file = file.path(OUT, "plink", "frq.frq")); stopifnot(nrow(FRQ) == length(all_ids), !anyNA(FRQ$MAF))
for (g in names(LD)) {
  L <- tg[gene == g & role == "lead"]; P <- tg[gene == g & role == "proxy"]
  h <- TAB[(SNP_A == L$bim_id & SNP_B == P$bim_id) | (SNP_A == P$bim_id & SNP_B == L$bim_id)]
  stopifnot(abs(abs(h$R) - abs(LD[[g]]$r_A1A1)) < 1e-4)
  # PHASE like "CG/GT": allele of SNP_A followed by its in-phase allele of SNP_B
  ph <- strsplit(h$PHASE, "/")[[1]]; a_first <- h$SNP_A == L$bim_id
  pairs <- t(vapply(ph, function(s) if (a_first) c(substr(s, 1, 1), substr(s, 2, 2)) else c(substr(s, 2, 2), substr(s, 1, 1)), c("", "")))
  LD[[g]]$phase_str <- h$PHASE; LD[[g]]$phase_pairs <- pairs   # rows: lead allele, in-phase proxy allele
  LD[[g]]$r_table <- h$R; LD[[g]]$frq <- FRQ[SNP %in% c(L$bim_id, P$bim_id)]
  logf("PLINK table %s: R=%.6f PHASE=%s (SNP_A=%s); square r=%.6f", g, h$R, h$PHASE, h$SNP_A, LD[[g]]$r_A1A1)
}
# in-place .bed read (BEDMatrix, memory-mapped, nothing copied): a second reader for the sign of r and the counted allele
suppressPackageStartupMessages(library(BEDMatrix))
X <- BEDMatrix(paste0(bfile, ".bed"), n = 503L, p = 22665064L)
bh <- ex$BIM503; bh[, lineno := as.integer(lineno)]; stopifnot(all(all_ids %in% bh$bim_id), anyDuplicated(bh[bim_id %in% all_ids]$bim_id) == 0L)
for (g in names(LD)) {
  L <- tg[gene == g & role == "lead"]; P <- tg[gene == g & role == "proxy"]
  il <- bh[bim_id == L$bim_id]$lineno; ip <- bh[bim_id == P$bim_id]$lineno
  G <- X[, c(il, ip)]
  rb <- cor(G[, 1], G[, 2], use = "complete.obs")
  f_cnt <- mean(G[, 1], na.rm = TRUE) / 2; f_A1 <- FRQ[SNP == L$bim_id]$MAF
  counted <- if (abs(f_cnt - f_A1) < 6e-4) "A1" else if (abs(f_cnt - (1 - f_A1)) < 6e-4) "A2" else "UNKNOWN"   # .frq prints 4 decimals
  if (counted == "UNKNOWN") stop(sprintf("BEDMatrix counted allele not identified for %s: mean/2 %.5f vs PLINK A1 freq %.4f", g, f_cnt, f_A1))
  rb_A1 <- rb   # both columns count the same allele position (A1 or A2), so r(A1,A1) = r(A2,A2)
  LD[[g]]$r_bedmatrix <- rb_A1; LD[[g]]$bed_counted <- counted
  logf("BEDMatrix in place %s: r = %.6f (counts %s; PLINK square r = %.6f; |diff| %.2e); missing genotypes lead %d proxy %d",
       g, rb_A1, counted, LD[[g]]$r_A1A1, abs(rb_A1 - LD[[g]]$r_A1A1), sum(is.na(G[, 1])), sum(is.na(G[, 2])))
  if (abs(rb_A1 - LD[[g]]$r_A1A1) > 1e-4) stop("PLINK and BEDMatrix disagree on r for ", g)
}
rm(X); invisible(gc())

# ---------------------------------------------------------------- UKB-PPP cis windows (exposure side at the proxies)
read_cis <- function(id) { d <- fread(file = rpath(id), colClasses = list(character = c("ID", "ALLELE0", "ALLELE1")))
  logf("cis window %s: %d rows", id, nrow(d)); d[, pos37 := as.integer(sapply(strsplit(ID, ":"), `[`, 2))]; d }
cis <- list()
for (id in na.omit(unique(PAIRS$exp_res))) cis[[id]] <- read_cis(id)
cis_row <- function(id, ppos, alleles) {
  d <- cis[[id]][pos37 == ppos]
  d <- d[toupper(ALLELE0) %in% toupper(alleles) & toupper(ALLELE1) %in% toupper(alleles)]
  if (nrow(d) > 1L) stop("ambiguous cis row"); if (!nrow(d)) return(NULL); d
}

# ---------------------------------------------------------------- section 1: proxy/allele ledger, direct-lead sensitivity, orientation
BONF_GRID <- uniqueN(ins$gene) * uniqueN(mr$outcome_id)
assert_target("family_grid", "Rule 1 (family over 15,584 tests)", BONF_GRID, 15584L, "unique instrument genes x unique outcome_id in mr_all_results")
BONF <- 0.05 / BONF_GRID
LEDGER <- list(); DIRECT <- list()
inphase_allele <- function(g, EA) {
  L <- tg[gene == g & role == "lead"]; P <- tg[gene == g & role == "proxy"]; r <- LD[[g]]$r_A1A1
  if (!(toupper(EA) %in% c(L$A1, L$A2))) stop("exposure EA not a panel allele for ", g)
  s <- if (toupper(EA) == L$A1) 1 else -1
  x_sq <- if (s * r > 0) P$A1 else P$A2
  pp <- LD[[g]]$phase_pairs; x_tab <- pp[pp[, 1] == toupper(EA), 2]
  if (length(x_tab) != 1L || x_tab != x_sq) stop(sprintf("in-phase allele disagreement for %s: square-sign %s vs PLINK in-phase %s", g, x_sq, paste(x_tab, collapse = ",")))
  list(X = x_sq, Y = if (x_sq == P$A1) P$A2 else P$A1, r_EA = s * r)
}
for (i in seq_len(nrow(PAIRS))) {
  pr <- PAIRS[i]; g <- pr$gene
  ii <- ins[gene == g]; EA <- toupper(ii$effect_allele); OA <- toupper(ii$other_allele); bx <- ii$beta; sx <- ii$se
  L <- tg[gene == g & role == "lead"]; P <- tg[gene == g & role == "proxy"]
  ip <- inphase_allele(g, EA)
  mrr <- mr[gene == g & outcome_id == pr$outcome_id]; stopifnot(nrow(mrr) == 1L)
  f4 <- fig4[gene == g & outcome_name == pr$outcome_name]
  # outcome rows
  ro_p <- out_row(pr$out_tag, L$chr, P$pos38, P$rsid, c(P$A1, P$A2))
  ro_l <- out_row(pr$out_tag, L$chr, L$pos38, L$rsid, c(EA, OA))
  if (is.null(ro_p)) stop("proxy row not found for ", pr$pair)
  ap <- align_row(ro_p, ip$X, ip$Y)
  b_w <- ap$beta / bx; se1 <- ro_p$se / abs(bx); se2 <- sqrt(se1^2 + b_w^2 * sx^2 / bx^2)
  p1 <- p2s(b_w / se1); p2 <- p2s(b_w / se2)
  rel <- function(a, b) abs(a - b) / max(abs(b), 1e-12)
  match_pipe <- rel(b_w, mrr$b) < 1e-4 && rel(se1, mrr$se) < 1e-4 && rel(p1, mrr$pval) < 1e-3
  # exposure side at the proxy (UKB-PPP cis window), oriented to the in-phase allele by label
  ce <- if (!is.na(pr$exp_res)) cis_row(pr$exp_res, P$pos37, c(P$A1, P$A2)) else NULL
  ce_l <- if (!is.na(pr$exp_res)) cis_row(pr$exp_res, L$pos37, c(EA, OA)) else NULL
  bxp <- sxp <- pxp <- fxp <- NA_real_
  if (!is.null(ce)) { s_ <- if (toupper(ce$ALLELE1) == ip$X) 1 else -1; bxp <- s_ * ce$BETA; sxp <- ce$SE; pxp <- 10^(-ce$LOG10P); fxp <- if (s_ == 1) ce$A1FREQ else 1 - ce$A1FREQ }
  if (!is.null(ce_l)) { s_ <- if (toupper(ce_l$ALLELE1) == EA) 1 else -1
    logf("%s: UKB-PPP cis row at the lead reproduces the instrument table: beta %.6f vs %.6f, se %.6f vs %.6f", g, s_ * ce_l$BETA, bx, ce_l$SE, sx)
    stopifnot(abs(s_ * ce_l$BETA - bx) < 1e-5, abs(ce_l$SE - sx) < 1e-6) }
  # frequencies
  frq <- LD[[g]]$frq; fA1l <- frq[SNP == L$bim_id]$MAF; fA1p <- frq[SNP == P$bim_id]$MAF
  f503_EA <- if (EA == L$A1) fA1l else 1 - fA1l; f503_X <- if (ip$X == P$A1) fA1p else 1 - fA1p
  # the reason the proxy was used (computed from the data, with the pipeline's harmonisation rule in 01_mr_screen.R)
  lead_pal <- is_pal(EA, OA); unin <- function(f) is.na(f) || (f > 0.42 && f < 0.58)
  al_l <- if (!is.null(ro_l)) align_row(ro_l, EA, OA) else NULL
  rsid_ok <- !is.null(ro_l) && grepl(paste0("(^|,)", L$rsid, "(,|$)"), ro_l$rsid_file)
  reason <- if (is.null(ro_l)) {
    sprintf("lead %s absent from %s: no row at GRCh38 %s:%d with alleles %s/%s and no row with rsID %s (whole file streamed)", L$rsid, pr$out_res, L$chr, L$pos38, EA, OA, L$rsid)
  } else if (!rsid_ok) {
    sprintf("lead present at GRCh38 %s:%d with alleles %s/%s but labelled %s in %s; the pipeline matches outcome rows by rsID, so the lead had no direct row and the proxy path ran", L$chr, L$pos38, ro_l$oa, ro_l$ea, ro_l$rsid_file, pr$out_res)
  } else if (lead_pal && (unin(ii$eaf) || unin(al_l$eaf))) {
    sprintf("lead present; palindromic %s/%s with EAF %.3f (UKB-PPP) and %.3f (%s), at least one inside (0.42, 0.58): strand not resolvable by frequency, so the pipeline routes the instrument to the proxy path (01_mr_screen.R read_local_outcome, the TwoSampleMR action=2 rule)", EA, OA, ii$eaf, al_l$eaf, pr$out_res)
  } else "UNEXPLAINED: the lead is present and harmonisable"
  wald_formula <- sprintf("Wald = beta_out(%s, in-phase allele %s; alleles %s/%s) / beta_exp(%s, EA %s); SE1 = SE_out/|beta_exp|; SE2 = sqrt(SE_out^2/beta_exp^2 + beta_out^2*SE_exp^2/beta_exp^4)",
                          P$rsid, ip$X, ip$X, ip$Y, L$rsid, EA)
  LEDGER[[pr$pair]] <- data.table(pair = pr$pair, gene = g, outcome = pr$outcome_name, outcome_resource = pr$out_res,
    lead = L$rsid, lead_panel_id = L$bim_id, lead_GRCh37 = sprintf("%s:%d", L$chr, L$pos37), lead_GRCh38 = sprintf("%s:%d", L$chr, L$pos38),
    queried = P$rsid, queried_GRCh37 = sprintf("%s:%d", P$chr, P$pos37),
    queried_GRCh38 = sprintf("%s:%s", P$chr, if (is.na(P$pos38)) paste0(ro_p$pos_file, " (from the outcome file; the GRCh37 position does not lift)") else P$pos38),
    lead_alleles_exposure = sprintf("EA %s / OA %s (UKB-PPP instrument table)", EA, OA), lead_palindromic = lead_pal,
    queried_alleles_panel = sprintf("A1 %s / A2 %s (503 bim)", P$A1, P$A2),
    queried_alleles_outcome_before = sprintf("effect %s / other %s (%s)", ro_p$ea, ro_p$oa, pr$out_res),
    alleles_after_harmonisation = sprintf("outcome beta re-expressed for in-phase allele %s (travels with lead EA %s; %s)", ip$X, EA, ap$how),
    eaf_exposure_lead = ii$eaf, eaf_ukbppp_queried = fxp, eaf_outcome_lead = if (!is.null(al_l)) al_l$eaf else NA_real_,
    eaf_outcome_queried = ap$eaf, eaf_503_lead_EA = f503_EA, eaf_503_queried_inphase = f503_X,
    r_plink_EA_inphase = abs(ip$r_EA), r_plink_A1A1 = LD[[g]]$r_A1A1, r2_plink = LD[[g]]$r2, phase_plink = LD[[g]]$phase_str,
    r2_pipeline_fig4 = f4$proxy_r2_min,
    beta_exp_lead = bx, se_exp_lead = sx, p_exp_lead = ii$pval, beta_exp_queried = bxp, se_exp_queried = sxp, p_exp_queried = pxp,
    beta_out_queried = ap$beta, se_out_queried = ro_p$se, p_out_queried = ro_p$p,
    beta_out_lead = if (!is.null(al_l)) al_l$beta else NA_real_, se_out_lead = if (!is.null(ro_l)) ro_l$se else NA_real_,
    p_out_lead = if (!is.null(ro_l)) ro_l$p else NA_real_, outcome_lead_rsid_in_file = if (!is.null(ro_l)) ro_l$rsid_file else NA_character_,
    proxy_reason = reason, wald_formula = wald_formula, wald_b = b_w, wald_se1 = se1, wald_p1 = p1, wald_se2 = se2, wald_p2 = p2,
    pipeline_b = mrr$b, pipeline_se = mrr$se, pipeline_p = mrr$pval, recompute_matches_pipeline = match_pipe,
    proxy_found_by = ro_p$how, lead_found_by = if (!is.null(ro_l)) ro_l$how else "not found")
  logf("ledger %s: proxy Wald %.6f (SE1 %.6f, p1 %.4g; SE2 %.6f, p2 %.4g) vs pipeline %.6f (%.6f, %.4g): match=%s; reason: %s",
       pr$pair, b_w, se1, p1, se2, p2, mrr$b, mrr$se, mrr$pval, match_pipe, reason)
  # direct-lead sensitivity (forward strand: outcome allele matched to the exposure EA by label)
  if (is.null(ro_l)) {
    DIRECT[[pr$pair]] <- data.table(pair = pr$pair, status = sprintf("NOT RUN: %s", reason))
  } else {
    bl <- al_l$beta; wl <- bl / bx; sl1 <- ro_l$se / abs(bx); sl2 <- sqrt(sl1^2 + wl^2 * sx^2 / bx^2)
    DIRECT[[pr$pair]] <- data.table(pair = pr$pair, status = "RUN", lead = L$rsid, found_by = ro_l$how, rsid_in_file = ro_l$rsid_file,
      beta = bl, se = ro_l$se, p = ro_l$p, wald = wl, wald_se = sl1, wald_p1 = p2s(wl / sl1), wald_se2 = sl2, wald_p2 = p2s(wl / sl2),
      strand = if (lead_pal) "forward-strand assumption (palindromic lead; matched by label)" else "non-palindromic; matched by label")
  }
}
LED <- rbindlist(LEDGER, fill = TRUE); DIR <- rbindlist(DIRECT, fill = TRUE)
stopifnot(nrow(LED) == 5L)
fwrite(LED, file.path(OUT, "ledger.tsv"), sep = "\t"); fwrite(DIR, file.path(OUT, "directlead.tsv"), sep = "\t")
for (i in seq_len(nrow(LED))) {
  x <- LED[i]; pid <- x$pair; res <- x$outcome_resource
  for (f in setdiff(names(x), c("pair", "gene", "outcome", "outcome_resource"))) {
    addr(sprintf("S1.ledger.%s.%s", pid, f), "s1 Analyses: proxy/allele ledger", sprintf("%s-%s ledger field %s", x$gene, x$outcome, f),
         resource = res, panel = if (grepl("503|plink", f)) "503 (PLINK --r, --keep-allele-order)" else "",
         se_order = if (grepl("se1|p1", f)) "1st" else if (grepl("se2|p2", f)) "2nd" else "", value = x[[f]])
  }
}
for (i in seq_len(nrow(DIR))) {
  x <- DIR[i]; pid <- x$pair; pr <- PAIRS[pair == pid]
  if (x$status != "RUN") {
    for (f in c("beta", "se", "p", "wald", "wald_se")) addr(sprintf("S1.directlead.%s.%s", pid, f), "s1 Analyses: direct-lead sensitivity", "outcome association at the lead and its Wald ratio", resource = pr$out_res, value = x$status)
    next
  }
  est <- sprintf("%s at lead %s, aligned to the exposure EA by label (%s)", pr$outcome_name, x$lead, x$strand)
  addr(sprintf("S1.directlead.%s.beta", pid), "s1 Analyses: direct-lead sensitivity", paste("outcome log-odds per EA;", est), resource = pr$out_res, value = x$beta, se = x$se, p = x$p,
       notes = sprintf("found by %s; rsID in file %s", x$found_by, x$rsid_in_file))
  addr(sprintf("S1.directlead.%s.se", pid), "s1 Analyses: direct-lead sensitivity", paste("SE of the outcome association;", est), resource = pr$out_res, value = x$se)
  addr(sprintf("S1.directlead.%s.p", pid), "s1 Analyses: direct-lead sensitivity", paste("outcome p (the resource's own p);", est), resource = pr$out_res, value = x$p)
  addr(sprintf("S1.directlead.%s.wald", pid), "s1 Analyses: direct-lead sensitivity", paste("Wald ratio per SD protein;", est), resource = pr$out_res, se_order = "1st",
       value = x$wald, se = x$wald_se, p = x$wald_p1, ci_lo = x$wald - qnorm(0.975) * x$wald_se, ci_hi = x$wald + qnorm(0.975) * x$wald_se,
       notes = sprintf("second-order SE %.6g, p %.4g; meets Bonferroni 0.05/%d: 1st %s, 2nd %s", x$wald_se2, x$wald_p2, BONF_GRID, x$wald_p1 < BONF, x$wald_p2 < BONF))
  addr(sprintf("S1.directlead.%s.wald_se", pid), "s1 Analyses: direct-lead sensitivity", paste("first-order Wald SE;", est), resource = pr$out_res, se_order = "1st", value = x$wald_se,
       notes = sprintf("second-order SE %.6g", x$wald_se2))
}

# exposure-side orientation at non-palindromic in-phase proxies (UKB-PPP cis windows)
for (i in seq_len(nrow(PAIRS))) {
  pr <- PAIRS[i]; g <- pr$gene; x <- LED[pair == pr$pair]; P <- tg[gene == g & role == "proxy"]
  ref <- "s1 Analyses: exposure-side orientation"
  if (is.na(pr$exp_res)) { addr(sprintf("S1.orient.%s", pr$pair), ref, "UKB-PPP beta at the proxy vs the lead", value = "NOT RUN: PECAM1 UKB-PPP cis window not local (EXP_UKBPPP_EUR_PECAM1_OID21131 not_found)"); next }
  if (is_pal(P$A1, P$A2)) { addr(sprintf("S1.orient.%s", pr$pair), ref, "UKB-PPP beta at the proxy vs the lead", value = "NOT RUN: proxy palindromic"); next }
  if (is.na(x$beta_exp_queried)) {
    addr(sprintf("S1.orient.%s", pr$pair), ref, "UKB-PPP beta at the proxy vs the lead", resource = pr$exp_res,
         value = sprintf("NOT RUN: proxy %s absent from the UKB-PPP window %s (no row at GRCh37 %s)", x$queried, pr$exp_res, x$queried_GRCh37),
         notes = if (g == "ABO") "ABO is not in the pre-declaration's parenthetical list (RAD23B, ADM); attempted because RESOURCES lists the ABO window for section 1 orientation" else "")
    next
  }
  ok <- sign(x$beta_exp_queried) == sign(x$beta_exp_lead)
  # phase-free proxy Wald: outcome and exposure both at the proxy, both on the in-phase allele by label
  wpf <- x$beta_out_queried / x$beta_exp_queried; spf <- x$se_out_queried / abs(x$beta_exp_queried)
  addr(sprintf("S1.orient.%s", pr$pair), ref, sprintf("sign of the UKB-PPP beta at proxy %s (in-phase allele) vs the lead %s beta (EA)", x$queried, x$lead),
       resource = pr$exp_res, panel = "503 (phase)", value = if (ok) "CONSISTENT" else "INCONSISTENT",
       notes = sprintf("UKB-PPP beta at proxy %.6f (SE %.6f, p %.3g, freq %.4f) vs lead %.6f (SE %.6f, EAF %.4f); r(EA, in-phase allele) %.4f; phase-free proxy Wald (both at the proxy) %.6f (SE1 %.6f, p %.4g) vs phase-based %.6f",
                       x$beta_exp_queried, x$se_exp_queried, x$p_exp_queried, x$eaf_ukbppp_queried, x$beta_exp_lead, x$se_exp_lead, x$eaf_exposure_lead,
                       x$r_plink_EA_inphase, wpf, spf, p2s(wpf / spf), x$wald_b))
}

# ---------------------------------------------------------------- external checks (outcome-check pre-specification), incl. RAD23B
half_unit <- function(s) {
  s <- sub("^[-+]", "", trimws(as.character(s)))
  if (grepl("[eE]", s)) { mant <- sub("[eE].*$", "", s); ex_ <- as.integer(sub("^.*[eE]", "", s))
    d <- if (grepl("\\.", mant)) nchar(sub("^[^.]*\\.", "", mant)) else 0L; return(0.5 * 10^(ex_ - d)) }
  d <- if (grepl("\\.", s)) nchar(sub("^[^.]*\\.", "", s)) else 0L
  0.5 * 10^(-d)
}
stopifnot(half_unit("-0.017") == 0.0005, abs(half_unit("0.0084") - 5e-5) < 1e-18, abs(half_unit("3.8e-02") - 5e-4) < 1e-18, half_unit("0.038") == 0.0005)
p_range <- function(B, hB, S, hS, f) {   # p over the rounding box; |z| is monotone in |B| and in S, so the corners bound it
  cand <- expand.grid(b = c(B - hB, B + hB), s = c(S - hS, S + hS))
  ps <- mapply(f, cand$b, cand$s); if ((B - hB) * (B + hB) <= 0) ps <- c(ps, 1)
  range(ps)
}
classify <- function(same_sign, p) if (same_sign && p < 0.05) "concordant" else if (same_sign) "directionally concordant, not significant" else if (p < 0.05) "discordant" else "inconclusive"
EXTRES <- list()
for (k in seq_len(nrow(EXT))) {
  e <- EXT[k]; ii <- ins[gene == e$gene]; xEA <- toupper(ii$effect_allele); xOA <- toupper(ii$other_allele); bx <- ii$beta; sx <- ii$se
  t <- tg[gene == e$gene & role %in% c("lead", "ext_instrument")]; stopifnot(nrow(t) == 1L)
  disc <- mr[gene == e$gene & outcome_name == e$outcome_name & method %in% c("Wald ratio", "Inverse variance weighted")]; stopifnot(nrow(disc) == 1L)
  if (e$res == "DIAM") {
    d <- ex$DIAM[SNP == sprintf("%s:%d", t$chr, t$pos37) & toupper(EA) %in% c(xEA, xOA) & toupper(NEA) %in% c(xEA, xOA)]
    drs <- ex$DIAMRS[rsID == ii$SNP]
    if (nrow(d) != 1L || nrow(drs) != 1L) stop("DIAMANTE row not unique for ", e$pair)
    same <- drs$SNP == d$SNP && drs$EA == d$EA && drs$NEA == d$NEA && as.numeric(drs$EAF) == as.numeric(d$EAF) && as.numeric(drs$Beta) == as.numeric(d$Beta) &&
            as.numeric(drs$SE) == as.numeric(d$SE) && abs(as.numeric(drs$Pvalue) - as.numeric(d$Pvalue)) < 1e-15
    logf("DIAMANTE %s: chrpos file %s %s/%s Beta %s SE %s P %s | rsID file Beta %s SE %s P %s -> numerically identical: %s", e$pair, d$SNP, d$EA, d$NEA, d$Beta, d$SE, d$Pvalue, drs$Beta, drs$SE, drs$Pvalue, same)
    if (!same) stop("the two DIAMANTE files disagree at ", e$pair)
    B_str <- d$Beta; S_str <- d$SE; P_str <- d$Pvalue; eaf <- as.numeric(d$EAF); ea <- toupper(d$EA); oa <- toupper(d$NEA)
    p_res <- as.numeric(P_str); res_lab <- "DIAMANTE 2018 European T2D without UKB (OUT_DIAMANTE2018_noUKB_EUR_chrpos; rsID file identical)"
  } else {
    d <- ex$CARDIO[`#CHROM` == t$chr & POS == as.character(t$pos37) & ID == ii$SNP]; if (nrow(d) != 1L) stop("CARDIoGRAMplusC4D row not unique for ", e$pair)
    fmt <- strsplit(d$FORMAT, ":")[[1]]; val <- strsplit(d[["ieu-a-7"]], ":")[[1]]; names(val) <- fmt
    B_str <- val[["ES"]]; S_str <- val[["SE"]]; P_str <- val[["LP"]]; eaf <- as.numeric(val[["AF"]]); ea <- toupper(d$ALT); oa <- toupper(d$REF)
    p_res <- 10^(-as.numeric(P_str)); res_lab <- "CARDIoGRAMplusC4D 1000G (OUT_CARDIOGRAMplusC4D_ieu-a-7; p = 10^-LP)"
  }
  B <- as.numeric(B_str); S <- as.numeric(S_str)
  if (ea == xEA && oa == xOA) { sg <- 1; eaf_al <- eaf } else if (ea == xOA && oa == xEA) { sg <- -1; eaf_al <- 1 - eaf } else stop("alleles do not match for ", e$pair)
  Bal <- sg * B; pal <- is_pal(xEA, xOA); maf_exp <- min(ii$eaf, 1 - ii$eaf)
  strand <- if (!pal) "non-palindromic" else if (maf_exp > 0.42) sprintf("forward-strand assumption (palindromic, MAF %.2f > 0.42: not orientable by frequency)", maf_exp) else if (abs(eaf_al - ii$eaf) > 0.10) "POSSIBLY STRAND-FLIPPED (palindromic; frequencies disagree by > 0.10)" else sprintf("palindromic; frequencies agree after alignment (resource %.3f vs UKB-PPP %.3f)", eaf_al, ii$eaf)
  w <- Bal / bx; s1 <- S / abs(bx); s2 <- sqrt(s1^2 + w^2 * sx^2 / bx^2); pw1 <- p2s(w / s1); pw2 <- p2s(w / s2)
  hB <- half_unit(B_str); hS <- half_unit(S_str)
  f1 <- function(b, s) p2s(b / s); f2 <- function(b, s) p2s(b / sqrt(s^2 + b^2 * sx^2 / bx^2))
  r1 <- p_range(Bal, hB, S, hS, f1); r2 <- p_range(Bal, hB, S, hS, f2)
  pr_res <- if (e$res == "DIAM") c(p_res - half_unit(P_str), p_res + half_unit(P_str)) else c(10^(-(as.numeric(P_str) + half_unit(P_str))), 10^(-(as.numeric(P_str) - half_unit(P_str))))
  ss <- sign(w) == sign(disc$b)
  EXTRES[[e$pair]] <- data.table(pair = e$pair, gene = e$gene, outcome = e$outcome_name, resource = res_lab, instrument = ii$SNP,
    resource_EA = ea, resource_OA = oa, resource_EAF = eaf, exposure_EA = xEA, eaf_aligned = eaf_al, eaf_ukbppp = ii$eaf, strand = strand,
    beta_printed = B_str, se_printed = S_str, p_printed = P_str, beta_aligned = Bal, se = S, p_resource = p_res, p_resource_lo = pr_res[1], p_resource_hi = pr_res[2],
    wald = w, wald_se1 = s1, wald_p1 = pw1, wald_p1_lo = r1[1], wald_p1_hi = r1[2], wald_se2 = s2, wald_p2 = pw2, wald_p2_lo = r2[1], wald_p2_hi = r2[2],
    discovery_b = disc$b, same_sign = ss, class_resource_p = classify(ss, p_res), class_1st = classify(ss, pw1), class_2nd = classify(ss, pw2),
    bonf6_resource = p_res < 0.05 / 6, bonf6_1st = pw1 < 0.05 / 6, bonf6_2nd = pw2 < 0.05 / 6)
  logf("external check %s: aligned beta %.4g (SE %s) -> Wald %.4f; resource p %.4g; Wald p 1st %.4g [%.4g, %.4g], 2nd %.4g [%.4g, %.4g]; class (resource p) %s; %s",
       e$pair, Bal, S_str, w, p_res, pw1, r1[1], r1[2], pw2, r2[1], r2[2], classify(ss, p_res), strand)
}
EXTD <- rbindlist(EXTRES); fwrite(EXTD, file.path(OUT, "ext_checks.tsv"), sep = "\t")
logf("external check FGFR4-T2D (East Asian, TPMI): NOT RUN by W1 (CKB exposure and TPMI outcome are outside W1 inputs; FGFR4-T2D is not one of the 81 European pairs)")

# RAD23B external-check quantities
xr <- EXTD[pair == "RAD23B_T2D"]
ref1 <- "s1 Wording rule: external check (classify on the resource p at the declared instrument)"
addr("S1.rad23b.ext.resource_p_lead", ref1, "DIAMANTE noUKB T2D p at the declared instrument rs7035725 (printed)", resource = xr$resource, value = xr$p_resource,
     ci_lo = xr$p_resource_lo, ci_hi = xr$p_resource_hi, p = xr$p_resource,
     notes = sprintf("printed '%s' (rounding interval in ci_lo..ci_hi); EA %s/%s EAF %.2f -> aligned to UKB-PPP EA %s: beta %+.3g (printed '%s'), SE '%s'; Wald %.4f vs discovery %.4f, same sign %s; class on resource p: %s; %s; label: UKB-free external consortium check, not independent replication",
                     xr$p_printed, xr$resource_EA, xr$resource_OA, xr$resource_EAF, xr$exposure_EA, xr$beta_aligned, xr$beta_printed, xr$se_printed, xr$wald, xr$discovery_b, xr$same_sign, xr$class_resource_p, xr$strand))
dpx <- ex$DIAM[SNP == sprintf("9:%d", tg[gene == "RAD23B" & role == "proxy"]$pos37)]; stopifnot(nrow(dpx) == 1L)
ipR <- inphase_allele("RAD23B", toupper(ins[gene == "RAD23B"]$effect_allele))
bpx <- if (toupper(dpx$EA) == ipR$X) as.numeric(dpx$Beta) else if (toupper(dpx$NEA) == ipR$X) -as.numeric(dpx$Beta) else NA
addr("S1.rad23b.ext.resource_p_proxy", ref1, "DIAMANTE noUKB T2D p at the discovery proxy rs7856112 (printed)", resource = xr$resource, value = as.numeric(dpx$Pvalue),
     p = as.numeric(dpx$Pvalue), ci_lo = as.numeric(dpx$Pvalue) - half_unit(dpx$Pvalue), ci_hi = as.numeric(dpx$Pvalue) + half_unit(dpx$Pvalue),
     notes = sprintf("printed '%s'; EA %s/%s Beta '%s' SE '%s'; in-phase allele with lead EA G is %s (PLINK, 503): aligned beta %+.3g, phase-based Wald %.4f (SE1 %.4f)",
                     dpx$Pvalue, dpx$EA, dpx$NEA, dpx$Beta, dpx$SE, ipR$X, bpx, bpx / ins[gene == "RAD23B"]$beta, as.numeric(dpx$SE) / abs(ins[gene == "RAD23B"]$beta)))
wnote <- "rounding interval over the printed Beta and SE (half a unit in the last printed digit; corners of the box)"
addr("S1.rad23b.ext.wald1_p", ref1, "first-order rounded-Wald p at rs7035725 (DIAMANTE Beta/SE as printed)", resource = xr$resource, se_order = "1st",
     value = xr$wald_p1, p = xr$wald_p1, ci_lo = xr$wald_p1_lo, ci_hi = xr$wald_p1_hi, notes = sprintf("Wald %.5f, SE1 %.5f; %s", xr$wald, xr$wald_se1, wnote))
addr("S1.rad23b.ext.wald1_lo", ref1, "lower end of the first-order rounded-Wald p interval", resource = xr$resource, se_order = "1st", value = xr$wald_p1_lo, notes = wnote)
addr("S1.rad23b.ext.wald1_hi", ref1, "upper end of the first-order rounded-Wald p interval", resource = xr$resource, se_order = "1st", value = xr$wald_p1_hi, notes = wnote)
addr("S1.rad23b.ext.wald2_p", ref1, "second-order (delta) rounded-Wald p at rs7035725", resource = xr$resource, se_order = "2nd",
     value = xr$wald_p2, p = xr$wald_p2, ci_lo = xr$wald_p2_lo, ci_hi = xr$wald_p2_hi, notes = sprintf("Wald %.5f, SE2 %.5f (exposure beta %.6f, SE %.6f exact); %s", xr$wald, xr$wald_se2, ins[gene == "RAD23B"]$beta, ins[gene == "RAD23B"]$se, wnote))
addr("S1.rad23b.ext.wald2_lo", ref1, "lower end of the second-order rounded-Wald p interval", resource = xr$resource, se_order = "2nd", value = xr$wald_p2_lo, notes = wnote)
addr("S1.rad23b.ext.wald2_hi", ref1, "upper end of the second-order rounded-Wald p interval", resource = xr$resource, se_order = "2nd", value = xr$wald_p2_hi, notes = wnote)
addr("S1.rad23b.ext.bonf6", ref1, "resource p at the declared instrument < 0.05/6 (six external checks)", resource = xr$resource, value = xr$bonf6_resource,
     notes = sprintf("0.05/6 = %.6f; resource p %.4g -> %s; first-order Wald p %.4g [%.4g, %.4g] -> %s; second-order %.4g [%.4g, %.4g] -> %s", 0.05 / 6, xr$p_resource, xr$bonf6_resource,
                     xr$wald_p1, xr$wald_p1_lo, xr$wald_p1_hi, xr$bonf6_1st, xr$wald_p2, xr$wald_p2_lo, xr$wald_p2_hi, xr$bonf6_2nd))

# ---------------------------------------------------------------- SE-convention table (sections 1 and Rule 1)
n_w <- sum(mr$method == "Wald ratio"); n_ivw <- sum(mr$method == "Inverse variance weighted"); n_eg <- sum(mr$method == "MR Egger"); n_wm <- sum(mr$method == "Weighted median")
l28 <- grep("second-order SE on", log28, value = TRUE); stopifnot(length(l28) == 1L)
n_w28 <- as.integer(sub(".*second-order SE on ([0-9]+) single-instrument.*", "\\1", l28)); n_m28 <- as.integer(sub(".*first-order kept for ([0-9]+) multi-instrument.*", "\\1", l28))
fdr28 <- as.integer(sub(".*MR FDR-significant pairs = ([0-9]+).*", "\\1", grep("MR FDR-significant pairs", log28, value = TRUE)))
logf("28_candidate_pool.log: second-order SE on %d single-instrument rows; first-order kept on %d multi-instrument rows; 'MR FDR-significant pairs = %d'", n_w28, n_m28, fdr28)
assert_target("se.wald_rows", "s1 SE-convention table (14,423 Wald rows)", n_w, 14423L, "mr_all_results rows with method == 'Wald ratio'")
assert_target("se.multi_rows", "s1 SE-convention table (114 = 98 IVW + 8 Egger + 8 weighted median)", sprintf("%d=%d+%d+%d", n_ivw + n_eg + n_wm, n_ivw, n_eg, n_wm), "114=98+8+8", "mr_all_results rows by method")
prim <- mr[method %in% c("Wald ratio", "Inverse variance weighted")]
if (anyDuplicated(prim[, .(gene, outcome_id)])) stop("more than one primary-method row per pair")
n_an <- nrow(prim); logf("primary-method rows (Wald or IVW, one per pair): %d; grid %d; untested cells padded at p = 1: %d", n_an, BONF_GRID, BONF_GRID - n_an)
ninst <- ins[, .N, by = gene]
one <- ins[gene %in% ninst[N == 1L]$gene, .(gene, bx = beta, sx = se, Fst = F_stat)]
logf("instrument table: F_stat vs (beta/se)^2, max relative difference %.2e (single-instrument genes)", one[, max(abs(Fst - (bx / sx)^2) / Fst)])
prim <- merge(prim, one, by = "gene", all.x = TRUE, sort = FALSE)
amb <- prim[method == "Wald ratio" & is.na(Fst)]
logf("join yield: Wald rows with a single-instrument gene %d/%d; Wald rows from multi-instrument genes (which instrument was usable is not recorded) %d: %s",
     sum(prim$method == "Wald ratio" & !is.na(prim$Fst)), n_w, nrow(amb), paste(unique(amb$gene), collapse = ","))
Frng <- ins[gene %in% amb$gene, .(Fmin = min(F_stat), Fmax = max(F_stat)), by = gene]
prim <- merge(prim, Frng, by = "gene", all.x = TRUE, sort = FALSE)
logf("Wald p consistency: max |pval - 2*Phi(-|b/se|)| / pval over Wald rows = %.2e", prim[method == "Wald ratio", max(abs(pval - p2s(b / se)) / pmax(pval, 1e-300))])
prim[, `:=`(F_a = fifelse(!is.na(Fst), Fst, Fmax), F_b = fifelse(!is.na(Fst), Fst, Fmin))]
prim[, `:=`(se2 = fifelse(method == "Wald ratio", sqrt(se^2 + b^2 / F_a), se), se2_b = fifelse(method == "Wald ratio", sqrt(se^2 + b^2 / F_b), se))]
prim[, `:=`(p1 = pval, p2 = p2s(b / se2), p2_b = p2s(b / se2_b))]
bhq <- function(p) p.adjust(c(p, rep(1, BONF_GRID - length(p))), method = "BH")[seq_along(p)]
prim[, `:=`(q1 = bhq(p1), q2 = bhq(p2), q2_b = bhq(p2_b))]
prim[, `:=`(bonf1 = p1 < BONF, bonf2 = p2 < BONF, bonf2_b = p2_b < BONF, bh1 = q1 < 0.05, bh2 = q2 < 0.05, bh2_b = q2_b < 0.05)]
amb_inert <- identical(prim$bonf2, prim$bonf2_b) && identical(prim$bh2, prim$bh2_b)
logf("multi-instrument Wald rows bracketed over their candidate instruments (largest vs smallest F): Bonferroni and BH sets identical under both = %s", amb_inert)
# reproduction of the table's own flags (logged only)
mr2 <- copy(mr); mr2[, fdr_bymethod := p.adjust(pval, "BH") < 0.05, by = method]
logf("reproduction: csv bonferroni_sig TRUE %d vs p < 0.05/%d on all rows %d (agree on %d/%d rows); csv fdr_sig TRUE %d vs BH within method family %d (agree %d/%d)",
     sum(mr$bonferroni_sig), BONF_GRID, sum(mr$pval < BONF), sum(mr$bonferroni_sig == (mr$pval < BONF)), nrow(mr), sum(mr$fdr_sig), sum(mr2$fdr_bymethod), sum(mr2$fdr_bymethod == mr$fdr_sig), nrow(mr))
fam_note <- sprintf("primary row per pair (Wald ratio or IVW; %d pairs); Bonferroni p < 0.05/%d = %.4g; BH over the %d-cell grid with %d untested cells padded at p = 1 (00_config.R primary_fdr_set); second order = delta method on the %d Wald rows, first order kept on the multi-instrument rows (28_candidate_pool.log); %d Wald rows from multi-instrument genes bracketed over candidate instruments, sets identical: %s",
                    n_an, BONF_GRID, BONF, BONF_GRID, BONF_GRID - n_an, n_w, nrow(amb), amb_inert)
mres <- "mr_all_results (FinnGen R10 CAD/HF/T2D/NAFLD/AF; MEGASTROKE; Wuttke 2019 eGFR/CKD)"
for (o in c("1st", "2nd")) {
  bf <- if (o == "1st") prim$bonf1 else prim$bonf2; bhh <- if (o == "1st") prim$bh1 else prim$bh2
  addr(paste0("S1.se.family_bonf.", o), "s1 SE-convention table; Rule 1 family", sprintf("pairs meeting study-wide Bonferroni, %s-order Wald SE", o), resource = mres, se_order = o, value = sum(bf), n = n_an, notes = fam_note)
  addr(paste0("S1.se.family_bh.", o), "s1 SE-convention table; Rule 1 family", sprintf("pairs in the family BH set (q < 0.05), %s-order Wald SE", o), resource = mres, se_order = o, value = sum(bhh), n = n_an, notes = fam_note)
}
assert_target("se.family_bonf_2nd", "s1 SE-convention table (assert 72 Bonferroni)", sum(prim$bonf2), 72L, "pipeline convention (second-order on Wald rows, first-order on multi-instrument rows), primary row per pair, p < 0.05/15,584")
assert_target("se.family_bh_2nd", "s1 SE-convention table (assert 248 BH)", sum(prim$bh2), 248L, "pipeline convention, BH over the 15,584-cell grid padded at p = 1")
logf("family counts: Bonferroni 1st %d / 2nd %d; BH 1st %d / 2nd %d", sum(prim$bonf1), sum(prim$bonf2), sum(prim$bh1), sum(prim$bh2))

# the 81 pairs: p, class and tier under each convention
t81 <- merge(tier, prim[, .(gene, outcome_name, method, b, se, se2, p1, p2, p2_b, q1, q2, bonf1, bonf2, bonf2_b, bh1, bh2, bh2_b, F_a)],
             by.x = c("gene", "outcome"), by.y = c("gene", "outcome_name"), all.x = TRUE, sort = FALSE)
logf("join yield: tiered81 -> primary rows %d/81 (null p1: %d)", sum(!is.na(t81$p1)), sum(is.na(t81$p1)))
stopifnot(nrow(t81) == 81L, !anyNA(t81$p1))
logf("tiered81 mr_pval vs first-order p: max relative difference %.2e", t81[, max(abs(mr_pval - p1) / p1)])
# sign convention of the tiered table (hygiene: asserted over every shared row): identity for binary outcomes, exact negation for eGFR (risk-aligned)
sg_id <- t81[outcome != "eGFR", max(abs(mr_beta - b))]; sg_neg <- t81[outcome == "eGFR", max(abs(mr_beta + b))]
logf("tiered81 mr_beta vs mr_all b: non-eGFR rows (%d) identical to 6 dp (max |diff| %.1e); eGFR rows (%d) exact negation (risk-aligned; max |mr_beta + b| %.1e)",
     t81[outcome != "eGFR", .N], sg_id, t81[outcome == "eGFR", .N], sg_neg)
stopifnot(sg_id < 1e-6, sg_neg < 1e-6)
cls <- function(bf, bhh) ifelse(bf, "Bonferroni", ifelse(bhh, "FDR", "Suggestive"))
tier_rule <- function(d, bonf) {   # P31R_SM_grading.md: any weakness -> Flagged; all four internal criteria + an evaluable arm -> Convergent; else Supported
  weak <- as.logical(d$has_major_weakness); all4 <- bonf & as.logical(d$crit_prior_stable) & as.logical(d$crit_susie_supportive)
  ifelse(weak, "Flagged", ifelse(all4 & d$n_evaluable >= 1, "Convergent", "Supported"))
}
t81[, class1 := cls(bonf1, bh1)]; t81[, class2 := cls(bonf2, bh2)]
rep_tier <- tier_rule(t81, as.logical(t81$crit_bonferroni))
logf("tier rule reproduces the declared evidence_tier from the table's own crit_bonferroni on %d/81 pairs", sum(rep_tier == t81$evidence_tier))
if (any(rep_tier != t81$evidence_tier)) stop("tier rule does not reproduce the declared tiers: ", paste(t81[rep_tier != evidence_tier, paste(gene, outcome)], collapse = ", "))
logf("crit_bonferroni (table) vs first-order Bonferroni here: agree on %d/81; mr_sig_class (table) vs first-order class here: agree on %d/81 (%s)",
     sum(as.logical(t81$crit_bonferroni) == t81$bonf1), sum(t81$mr_sig_class == t81$class1),
     paste(t81[mr_sig_class != class1, sprintf("%s-%s table %s vs %s", gene, outcome, mr_sig_class, class1)], collapse = "; "))
t81[, tier1 := tier_rule(t81, bonf1)]; t81[, tier2 := tier_rule(t81, bonf2)]
t81[, `:=`(mover_bonf = bonf1 != bonf2, mover_fdr = bh1 != bh2, mover_tier = tier1 != tier2)]
t81 <- merge(t81, EXTD[, .(gene, outcome, ext_resource = resource, ext_p_resource = p_resource, ext_p_wald1 = wald_p1, ext_p_wald2 = wald_p2,
                           ext_class_resource_p = class_resource_p, ext_class_1st = class_1st, ext_class_2nd = class_2nd)], by = c("gene", "outcome"), all.x = TRUE, sort = FALSE)
SE81 <- t81[, .(gene, outcome, method, b, se1 = se, p1, se2, p2, bh_q1 = q1, bh_q2 = q2, bonf1, bonf2, bh1, bh2, class1, class2,
                tier_declared = evidence_tier, tier1, tier2, mover_bonf, mover_fdr, mover_tier, ext_resource, ext_p_resource, ext_p_wald1, ext_p_wald2,
                ext_class_resource_p, ext_class_1st, ext_class_2nd)]
setorder(SE81, p1); fwrite(SE81, file.path(OUT, "se_convention_table.tsv"), sep = "\t")
nm <- function(d, col) paste(d[, sprintf("%s-%s (%s p %.3g -> %s p %.3g)", gene, outcome, if (col == "bonf") ifelse(bonf1, "Bonf", "notBonf") else ifelse(bh1, "BH", "notBH"),
                                          p1, if (col == "bonf") ifelse(bonf2, "Bonf", "notBonf") else ifelse(bh2, "BH", "notBH"), p2)], collapse = "; ")
mb <- t81[mover_bonf == TRUE][order(p1)]; mf <- t81[mover_fdr == TRUE][order(p1)]; mt <- t81[mover_tier == TRUE]
for (o in c("1st", "2nd")) {
  addr(paste0("S1.se.core_bonf.", o), "s1 SE-convention table (within the 81)", sprintf("of the 81 graded pairs, number meeting Bonferroni, %s-order SE", o), resource = mres, se_order = o,
       value = sum(if (o == "1st") t81$bonf1 else t81$bonf2), n = 81L)
  addr(paste0("S1.se.core_fdr.", o), "s1 SE-convention table (within the 81)", sprintf("of the 81 graded pairs, number in the family BH set, %s-order SE (Bonferroni pairs included)", o), resource = mres, se_order = o,
       value = sum(if (o == "1st") t81$bh1 else t81$bh2), n = 81L,
       notes = sprintf("classes %s: %s", o, paste(names(table(if (o == "1st") t81$class1 else t81$class2)), table(if (o == "1st") t81$class1 else t81$class2), sep = " ", collapse = ", ")))
}
addr("S1.se.movers_bonf", "s1 SE-convention table (named Bonferroni movers within the 81)", "pairs whose Bonferroni status differs between first- and second-order SE", resource = mres, se_order = "1st vs 2nd",
     value = nrow(mb), n = 81L, notes = nm(mb, "bonf"))
addr("S1.se.movers_fdr", "s1 SE-convention table (named FDR movers within the 81)", "pairs whose family-BH status differs between first- and second-order SE", resource = mres, se_order = "1st vs 2nd",
     value = nrow(mf), n = 81L, notes = nm(mf, "fdr"))
addr("S1.se.movers_tier", "s1 SE-convention table (tier under each convention; descriptive, Rule 2)", "pairs whose tier by the SM grading rule differs between conventions", resource = mres, se_order = "1st vs 2nd",
     value = nrow(mt), n = 81L, notes = paste(mt[, sprintf("%s-%s %s -> %s", gene, outcome, tier1, tier2)], collapse = "; "))
assert_target("se.movers_bonf_n", "s1 SE-convention table (11 Bonferroni movers)", nrow(mb), 11L, "within the 81: bonf1 != bonf2")
assert_target("se.movers_fdr_n", "s1 SE-convention table (9 FDR movers)", nrow(mf), 9L, "within the 81: bh1 != bh2 (family BH)")
logf("Bonferroni movers (%d): %s", nrow(mb), nm(mb, "bonf")); logf("FDR movers (%d): %s", nrow(mf), nm(mf, "fdr"))

# ---------------------------------------------------------------- RAD23B: r2, tier, deCODE, the cut at lead / proxy / second order
LR <- LED[pair == "RAD23B_T2D"]; DR <- DIR[pair == "RAD23B_T2D"]; T81R <- t81[gene == "RAD23B" & outcome == "T2D"]
addr("S1.rad23b.r2_plink", "s1 Wording rule: first mention (Stage 1 PLINK r2)", "r2 between lead rs7035725 and proxy rs7856112", resource = "1000 Genomes EUR (LD503)",
     panel = "503 (PLINK --r square, --keep-allele-order)", value = LR$r2_plink,
     notes = sprintf("signed r(A1 C, A1 G) %.6f; PLINK in-phase %s (lead EA G travels with proxy %s); BEDMatrix in-place r %.6f; pipeline (FIG04) r2 %s", LR$r_plink_A1A1, LR$phase_plink,
                     ipR$X, LD$RAD23B$r_bedmatrix, LR$r2_pipeline_fig4))
addr("S1.rad23b.tier", "s1 Wording rule (RAD23B keeps its declared Convergent tier)", "declared evidence tier of RAD23B-T2D", resource = "TAB_tiered81", value = T81R$evidence_tier,
     notes = sprintf("SM grading rule, descriptive: first-order %s (class %s), second-order %s (class %s); Rule 2: no tier is rebuilt", T81R$tier1, T81R$class1, T81R$tier2, T81R$class2))
dr <- dec[gene == "RAD23B" & outcome_name == "T2D" & status == "MR_OK"]; stopifnot(nrow(dr) == 1L)
F_back <- dr$decode_b^2 / (dr$delta_se^2 - dr$decode_se^2); dse_chk <- sqrt(dr$decode_se^2 + dr$decode_b^2 / dr$decode_F_stat)
addr("S1.rad23b.decode_F", "s1 (deCODE F)", "deCODE exposure-side F for RAD23B (SomaScan 12522_6)", resource = "TAB_decode_mr_results", value = dr$decode_F_stat,
     notes = sprintf("table value; back-derived from decode_b, decode_se and delta_se: %.4f; weak_instrument=%s (F < 10 rule); palindromic=%s", F_back, dr$weak_instrument, dr$palindromic))
p2dec <- p2s(dr$decode_b / dr$delta_se)
addr("S1.rad23b.decode_p2", "s1 (deCODE second-order delta p, recomputed)", "deCODE-instrumented RAD23B-T2D Wald ratio, second-order (delta) p", resource = "TAB_decode_mr_results", se_order = "2nd",
     value = p2dec, p = p2dec, se = dr$delta_se,
     notes = sprintf("2*Phi(-|decode_b %.6f / delta_se %.6f|); table delta_pval %.6g (|diff| %.1e); delta_se recomputed from decode_se and F: %.6f; first-order p %.4g", dr$decode_b, dr$delta_se, dr$delta_pval, abs(p2dec - dr$delta_pval), dse_chk, p2s(dr$decode_b / dr$decode_se)))
cut_lead <- DR$wald_p1 < BONF; cut_proxy <- LR$wald_p1 < BONF; cut_2nd <- LR$wald_p2 < BONF
addr("S1.rad23b.cut.lead", "s1 Wording rule: first mention (meets the cut at the lead)", "FinnGen T2D p at the lead rs7035725 < 0.05/15,584", resource = "OUT_FG_R10_T2D", se_order = "1st",
     value = cut_lead, p = DR$wald_p1, notes = sprintf("lead p %.4g (forward-strand assumption); threshold %.4g", DR$wald_p1, BONF))
addr("S1.rad23b.cut.proxy", "s1 Wording rule: first mention (meets the cut at the proxy)", "FinnGen T2D p at the proxy rs7856112 < 0.05/15,584", resource = "OUT_FG_R10_T2D", se_order = "1st",
     value = cut_proxy, p = LR$wald_p1, notes = sprintf("proxy p %.4g (= pipeline p %.4g); threshold %.4g", LR$wald_p1, LR$pipeline_p, BONF))
addr("S1.rad23b.cut.second_order", "s1 Wording rule: first mention (meets the cut under second-order SE)", "proxy-based Wald p with second-order SE < 0.05/15,584", resource = "OUT_FG_R10_T2D", se_order = "2nd",
     value = cut_2nd, p = LR$wald_p2, notes = sprintf("second-order p %.4g at the proxy; at the lead %.4g; family BH q (2nd) %.4g", LR$wald_p2, DR$wald_p2, T81R$q2))

# ---------------------------------------------------------------- branches (section 1 wording rule; Rule 2 first-mention disclosures)
fp <- function(p) formatC(p, digits = 2, format = "g")
addb("s1 Wording rule: tier", sprintf("evidence_tier=%s (TAB_tiered81); grading rule 1st=%s, 2nd=%s", T81R$evidence_tier, T81R$tier1, T81R$tier2),
     "RAD23B keeps its declared Convergent tier (fixed by the rule, whatever the sensitivities show)",
     "RAD23B-T2D keeps its declared Convergent tier.",
     sprintf("Rule 2: sensitivities are descriptive; under second-order SE the pair's class would be %s and its grade %s (disclosed at first mention)", T81R$class2, T81R$tier2))
fm_branch <- sprintf("meets the cut at the lead: %s; at the proxy: %s; under second-order SE: %s", ifelse(cut_lead, "yes", "no"), ifelse(cut_proxy, "yes", "no"), ifelse(cut_2nd, "yes", "no"))
addb("s1 Wording rule: first mention",
     sprintf("p_lead=%.4g (S1.directlead.RAD23B_T2D.p); p_proxy=%.4g (S1.ledger.RAD23B_T2D.p_out_queried); r2=%.4f (S1.rad23b.r2_plink); p_2nd=%.4g; threshold 0.05/15584=%.4g", DR$wald_p1, LR$wald_p1, LR$r2_plink, LR$wald_p2, BONF),
     fm_branch,
     sprintf("RAD23B-T2D: lead rs7035725 p %s, proxy rs7856112 p %s (Stage 1 PLINK r2 %.3f, 1000 Genomes European panel); the pair %s the study-wide cut (3.2 x 10^-6) at the lead, %s it at the proxy, and %s it under second-order SE (p %s).",
             fp(DR$wald_p1), fp(LR$wald_p1), LR$r2_plink, ifelse(cut_lead, "meets", "does not meet"), ifelse(cut_proxy, "meets", "does not meet"), ifelse(cut_2nd, "meets", "does not meet"), fp(LR$wald_p2)),
     sprintf("'proxy' label at first mention (Principle 0 list); the lead estimate rests on the forward-strand assumption (palindromic G/C, EAF %.2f); Rule 2: class %s under second-order SE", ins[gene == "RAD23B"]$eaf, T81R$class2))
any_ge <- any(c(xr$p_resource, xr$wald_p1, xr$wald_p2) >= 0.05)
straddle <- any(c(xr$p_resource_lo, xr$wald_p1_lo, xr$wald_p2_lo) < 0.05 & c(xr$p_resource_hi, xr$wald_p1_hi, xr$wald_p2_hi) >= 0.05)
clause <- any_ge || straddle
addb("s1 Wording rule: external check",
     sprintf("resource p %.4g [%.4g, %.4g]; 1st-order rounded-Wald p %.4g [%.4g, %.4g]; 2nd-order %.4g [%.4g, %.4g]; Wald sign %s vs discovery %s; any printed p >= 0.05: %s; any interval straddles 0.05: %s",
             xr$p_resource, xr$p_resource_lo, xr$p_resource_hi, xr$wald_p1, xr$wald_p1_lo, xr$wald_p1_hi, xr$wald_p2, xr$wald_p2_lo, xr$wald_p2_hi, sign(xr$wald), sign(xr$discovery_b), any_ge, straddle),
     sprintf("classified on the resource p at the declared instrument: %s; p-source/rounding clause %s; 0.05/6: %s", xr$class_resource_p, if (clause) "ADDED" else "not added", if (xr$bonf6_resource) "met" else "not met"),
     sprintf("In DIAMANTE 2018 (European, without UK Biobank) the RAD23B instrument rs7035725 gave %s (p %s; first-order rounded-Wald p %s, rounding interval %s-%s; second-order %s, %s-%s)%s; the check %s the six-check threshold 0.05/6 (0.0083).",
             if (xr$class_resource_p == "concordant") "a concordant estimate" else xr$class_resource_p, fp(xr$p_resource), fp(xr$wald_p1), fp(xr$wald_p1_lo), fp(xr$wald_p1_hi), fp(xr$wald_p2), fp(xr$wald_p2_lo), fp(xr$wald_p2_hi),
             if (clause) "; the classification depends on the p source and on rounding of the printed DIAMANTE estimates" else "", if (xr$bonf6_resource) "meets" else "does not meet"),
     "labels that travel (outcome-check pre-specification): 'forward-strand assumption' (palindromic, MAF 0.45); 'UKB-free external consortium check, not independent replication' with the FinnGen overlap r0 (not computed by W1); the RAD23B p-source/rounding clause and its 0.05/6 result (Principle 0)")
for (i in seq_len(nrow(PAIRS))) {
  pr <- PAIRS[i]; x <- LED[pair == pr$pair]; dl <- DIR[pair == pr$pair]; tt <- t81[gene == pr$gene & outcome == pr$outcome_name]
  addb("Principle 0 / Rules 2-3(i): first-mention 'proxy' label",
       sprintf("proxy %s for lead %s (r2 %.4f, 503); reason: %s; direct-lead %s", x$queried, x$lead, x$r2_plink, x$proxy_reason,
               if (dl$status == "RUN") sprintf("Wald %.4g p %.3g", dl$wald, dl$wald_p1) else dl$status),
       "proxy-marked pair: the 'proxy' label travels with the pair",
       sprintf("%s-%s is estimated at proxy %s (r2 %.3f with lead %s)%s.", pr$gene, pr$outcome_name, x$queried, x$r2_plink, x$lead,
               if (dl$status == "RUN") sprintf("; at the lead itself the estimate is %.3g (p %s)", dl$wald, fp(dl$wald_p1)) else "; the lead is absent from the outcome resource"),
       sprintf("'proxy' label; class under second-order SE: %s (first-order %s)", tt$class2, tt$class1))
}
if (nrow(mb) || nrow(mf)) addb("Rule 2: a pair whose class differs under a sensitivity is disclosed at first mention",
     sprintf("Bonferroni movers %d; FDR movers %d (S1.se.movers_bonf, S1.se.movers_fdr)", nrow(mb), nrow(mf)),
     "class differs between first- and second-order SE: disclose at first mention",
     paste(unique(rbind(mb, mf)[, sprintf("%s-%s: %s under first-order SE, %s under second-order SE", gene, outcome, class1, class2)]), collapse = " | "),
     "descriptive only (Rule 2): no tier, class or candidate set is rebuilt")

# ---------------------------------------------------------------- write outputs and verify them
RD <- rbindlist(RESULTS)
if (anyDuplicated(RD$quantity_id)) stop("duplicate quantity ids: ", paste(RD$quantity_id[duplicated(RD$quantity_id)], collapse = ","))
fwrite(RD, file.path(OUT, "results.tsv"), sep = "\t", quote = FALSE)
BD <- rbindlist(BR); fwrite(BD, file.path(OUT, "branches.tsv"), sep = "\t", quote = FALSE)
chk <- fread(file = file.path(OUT, "results.tsv"), sep = "\t", quote = "", colClasses = "character")
stopifnot(nrow(chk) == nrow(RD), ncol(chk) == 15L)
na_vals <- chk[value %in% c("", "NA")]
logf("results.tsv: %d rows, %d columns; rows with empty/NA value: %d (%s)", nrow(chk), ncol(chk), nrow(na_vals), paste(na_vals$quantity_id, collapse = ", "))
need_ids <- c(sprintf("S1.ledger.%s.%s", rep(PAIRS$pair, each = 3), c("r2_plink", "wald_b", "proxy_reason")),
              sprintf("S1.directlead.%s.%s", rep(PAIRS$pair, each = 5), c("beta", "se", "p", "wald", "wald_se")), sprintf("S1.orient.%s", PAIRS$pair),
              paste0("S1.se.", c("family_bonf.1st", "family_bonf.2nd", "family_bh.1st", "family_bh.2nd", "core_bonf.1st", "core_bonf.2nd", "core_fdr.1st", "core_fdr.2nd", "movers_bonf", "movers_fdr")),
              paste0("S1.rad23b.ext.", c("resource_p_lead", "resource_p_proxy", "wald1_p", "wald1_lo", "wald1_hi", "wald2_p", "wald2_lo", "wald2_hi", "bonf6")),
              paste0("S1.rad23b.", c("r2_plink", "tier", "decode_F", "decode_p2")))
miss <- setdiff(need_ids, chk$quantity_id); if (length(miss)) stop("briefed quantity ids missing: ", paste(miss, collapse = ","))
logf("all %d briefed quantity ids present; branches.tsv %d rows; se_convention_table.tsv %d rows; ledger.tsv %d rows; ext_checks.tsv %d rows", length(need_ids), nrow(BD), nrow(SE81), nrow(LED), nrow(EXTD))
logf("ASSERT summary: %d assertion targets; failed: %s", sum(startsWith(RD$quantity_id, "S1.assert.")), if (length(AFAIL)) paste(AFAIL, collapse = " | ") else "none")
snap1 <- staging_snapshot("end")
pre <- fread(file = file.path(OUT, "staging_snapshot_pre_extract.tsv"))
same_all <- identical(snap0$digest, snap1$digest) && identical(digest(paste(pre$file, pre$md5, collapse = "\n"), algo = "md5", serialize = FALSE), snap1$digest)
logf("staging snapshot pre-extract / start / end identical: %s", same_all); if (!same_all) stop("files under V4_local/staging/l3-finemap changed during the run")
logf("output md5: results.tsv %s; branches.tsv %s; ledger.tsv %s; se_convention_table.tsv %s; ext_checks.tsv %s; directlead.tsv %s",
     md5(file.path(OUT, "results.tsv")), md5(file.path(OUT, "branches.tsv")), md5(file.path(OUT, "ledger.tsv")), md5(file.path(OUT, "se_convention_table.tsv")),
     md5(file.path(OUT, "ext_checks.tsv")), md5(file.path(OUT, "directlead.tsv")))
logf("W1 finished in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
