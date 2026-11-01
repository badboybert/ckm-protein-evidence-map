#!/usr/bin/env Rscript
# ============================================================================
# 39_palindromic_freq_gate.R -- MANDATORY, COMPREHENSIVE, BLOCKING frequency gate for palindromic fit-SNPs.
#
# WHY (cold re-review High, 2026-07-17). The strand-align fix (R/strand_align.R) leaves A/T and C/G palindromic
# SNPs UNCHANGED in the SuSiE fits, because allele matching cannot decide their strand. That is defensible ONLY
# if a post-run check verifies, for every palindromic SNP that entered a fit, that its LD row (oriented to the
# .bim A1) is physically the pQTL EFFECT allele -- i.e. that leaving it unchanged was correct. A palindrome
# whose .bim A1 is physically the OTHER allele is silently anti-aligned (z is ea-oriented, R is oa-oriented),
# corrupting that fit. This script is that check, and it BLOCKS (non-zero exit) on any violation.
#
# HOW. During the re-run each engine's strand_align_R() call appends its palindromic fit-SNPs to a sink dir
# (env CKM_PALINDROME_SINK=<dir>), from the single choke point, so coverage is comprehensive across
# 15/26/16/29/31 with no per-engine plumbing. Here we:
#   1. read the sink (label, rsid, ea_pqtl, oa_pqtl, bim_a1, z_pqtl);
#   2. look up the pQTL EAF = A1FREQ = freq(ea_pqtl) per (gene, rsid) from the UKB-PPP sumstats (get_pqtl_region);
#   3. look up the 1000G EUR freq of .bim A1 (EUR_freqcnt.frq.counts, made with --keep-allele-order so its A1 IS
#      the LD/.bim A1) = C1/(C1+C2);
#   4. classify: freq(bim_A1) ~ EAF  -> ALIGNED  (leaving unchanged is correct);
#               freq(bim_A1) ~ 1-EAF -> SWAPPED  (should have been flipped -> VIOLATION);
#               EAF ~ 0.5            -> AMBIGUOUS (freq cannot resolve; reported, cannot be auto-cleared);
#               neither within tol / no EAF / no 1000G freq -> UNRESOLVED (reported).
#
# EXIT (editor ruling 2026-07-18): 0 UNLESS there is a genuine SWAPPED palindrome in an EVALUABLE (non-MHC-
# excluded) pair -- the only case that corrupts a USED SuSiE fit. near-0.5 AMBIGUOUS (frequency cannot decide
# strand) and UNRESOLVED (rsID-mapping gap / allele-label technicals) are RECORDED and DISCLOSED, not blocked;
# a SWAPPED in an MHC-excluded pair is immaterial (its SuSiE number is not_evaluable and unused). The PASS marker
# records the full audit (aligned rate, named swaps) as the source for the strand/QC disclosure. READ-ONLY.
# ============================================================================

# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"))
suppressMessages({ library(data.table) })

B  <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
# PHASE1_SET, NOT PHASE1_TAG. `_phase1_set.R` (sourced above) defines PHASE1_SET; PHASE1_TAG does
# not exist until R/00_config.R runs, and that is 18 lines below this one. So this line resolved
# only when the script happened to be sourced into a session that had already loaded the config --
# i.e. this gate, which BLOCKS the SuSiE certification, could never be run as a standalone Rscript.
# 00_config.R:1089 is `PHASE1_TAG <- PHASE1_SET`, so the two are the same string by construction and
# this substitution changes no path.
RC <- file.path(B, file.path("results", paste0("repro_", PHASE1_SET)))
D  <- file.path(RC, "phase3_susie")

# tunables (mirror the reverted in-engine constants; single source here)
FREQ_TOL   <- 0.20   # nearer of {EAF, 1-EAF} to the 1000G freq must be within this, else UNRESOLVED
FREQ_AMBIG <- 0.05   # |EAF-0.5| below this -> AMBIGUOUS (freq cannot distinguish aligned from swapped)
FREQ_SEP   <- 0.10   # aligned vs swapped distances must differ by MORE than this; a near-tie is uninformative -> AMBIGUOUS
Z_AMBIG    <- 5      # |z| threshold to flag HIGH-signal AMBIGUOUS/UNRESOLVED palindromes in the audit (disclosed, not blocked -- editor ruling 2026-07-18)

# sink dir MUST match what strand_align_R() writes to: env CKM_PALINDROME_SINK (single, env-only mechanism).
# Read the SAME env var (High: selector mismatch), falling back to the canonical default only if unset.
SINK <- if (exists("PALIN_SINK")) PALIN_SINK else { .e <- Sys.getenv("CKM_PALINDROME_SINK", ""); if (nzchar(.e)) .e else file.path(D, "palindrome_sink") }
FRQ  <- file.path(B, "data/ld_reference/1kg_eur/EUR_freqcnt.frq.counts")
INSTR<- if (exists("PALIN_INSTR")) PALIN_INSTR else PHASE1_INSTRUMENTS
OUT  <- file.path(D, "palindromic_freq_gate.csv")
PASS_MARK <- file.path(D, "palindromic_freq_gate.PASS")   # written ONLY on a clean pass; 33 requires it to certify the dashed layer
if (file.exists(PASS_MARK)) file.remove(PASS_MARK)          # clear any stale marker up front so a prior pass can never linger

source(file.path(B, "R/00_config.R"))
COLOC_SOURCE_ONLY <- TRUE; suppressMessages(source(file.path(B, "R/coloc.R"))); rm(COLOC_SOURCE_ONLY)

fail <- function(msg) { cat(sprintf("\n[39 GATE FAIL] %s\n", msg)); quit(status = 1L, save = "no") }   # note: never writes PASS_MARK

# ── 1. read the sink ─────────────────────────────────────────────────────────
if (!dir.exists(SINK)) fail(sprintf("sink dir not found (%s) -- was CKM_PALINDROME_SINK set during the re-run? Cannot verify palindromes that were never captured.", SINK))
RUN_ID <- Sys.getenv("CKM_RUN_ID", "")
if (!nzchar(RUN_ID)) fail("CKM_RUN_ID not set -- refusing to certify palindromes without a run identity binding the sink and PASS marker to THIS run. Export CKM_RUN_ID before the capture run (see RERUN_RUNBOOK.md).")

# SEEN markers: strand_align_R writes one per fit-call (even at 0 palindromes), so their presence is POSITIVE
# proof capture ran for every fit -- distinguishing a true zero-palindrome run from capture-off, and letting the
# coverage check below detect a fit that was skipped (resumed checkpoint / cache hit) (cold re-review High/Medium).
seen_files <- list.files(SINK, pattern = "^seen_.*\\.csv$", full.names = TRUE)
if (!length(seen_files)) fail(sprintf("no seen_*.csv in %s -- palindrome capture did NOT run (CKM_PALINDROME_SINK unset during the fits, or wrong dir). Refusing to certify.", SINK))
seen_dt <- rbindlist(lapply(seen_files, fread), use.names = TRUE, fill = TRUE)
if (any(!c("label", "run_id", "n_palindromic") %in% names(seen_dt))) fail("seen markers malformed (need label/run_id/n_palindromic) -- older helper; clear the sink and re-capture.")
.sbad <- unique(as.character(seen_dt$run_id[as.character(seen_dt$run_id) != RUN_ID]))
if (length(.sbad)) fail(sprintf("seen markers carry other run_id(s) [%s] != current [%s] -- stale/mixed sink; clear it and re-run capture.", paste(head(.sbad, 3), collapse = ","), RUN_ID))
seen_labels <- unique(seen_dt$label)
cat(sprintf("[39] capture ran: %d fit-call(s), %d distinct label(s), %d palindromic SNP(s) recorded\n",
            nrow(seen_dt), length(seen_labels), sum(seen_dt$n_palindromic)))

# COVERAGE (HARD): every SHIPPED pair that was actually FIT -- >=10 LD SNPs after strand-align drops, non-MHC,
# which is EXACTLY the set strand_align_R records a seen marker for -- must appear in seen_labels. Filtering the
# output CSVs on n_snps_ld>=10 mirrors the sink's fit-gate, so this cannot false-fail on MHC/not-evaluable or
# <10-SNP skips, and DOES catch a fit skipped by a resumed checkpoint / cache hit (cold re-review High).
seen_total <- sum(seen_dt$n_palindromic)
# EVERY return path that used to yield character(0) now FAILS instead. Those silent empties were the
# gate's own escape hatch: `expected` is the set of fits the capture is obliged to cover, and it is
# built from these same CSVs, so a missing file, a renamed column or a truncated table shrank the
# obligation to match whatever was captured -- and an empty `expected` made the coverage check pass
# vacuously while the PASS marker said nothing about how many fits were expected. A gate whose
# obligation is derived from a file must fail when that file cannot define the obligation.
.fit_labels <- function(f) {
  if (!file.exists(f))
    fail(sprintf("cannot establish the capture obligation: %s is absent. `expected` is built from this file, so proceeding would shrink the coverage requirement to whatever happened to be captured.", basename(f)))
  d <- fread(f)
  if (!all(c("gene", "outcome_name") %in% names(d)))
    fail(sprintf("cannot establish the capture obligation: %s lacks gene/outcome_name (has: %s) -- the pair labels cannot be built.", basename(f), paste(head(names(d), 12), collapse = ", ")))
  nl <- if ("n_snps_ld" %in% names(d)) "n_snps_ld" else if ("strict_n_snps_ld" %in% names(d)) "strict_n_snps_ld" else NA_character_
  if (is.na(nl))
    fail(sprintf("cannot establish the capture obligation: %s has neither n_snps_ld nor strict_n_snps_ld, so the >=10-SNP fit gate cannot be mirrored and every row would drop out of `expected`.", basename(f)))
  if ("mhc_excluded" %in% names(d)) d <- d[!(mhc_excluded %in% c(TRUE, "TRUE"))]
  d <- d[!is.na(get(nl)) & get(nl) >= 10]
  # paste0() treats a ZERO-LENGTH argument as "" whenever another argument has length 1, so once the
  # filters emptied the table this returned the phantom label "~" rather than nothing -- coverage
  # then reported a missing pair that never existed, and the empty-obligation guard below could never
  # fire because the obligation was never empty. Return early on an emptied table.
  if (!nrow(d)) return(character(0))
  oc <- if ("outcome_id" %in% names(d)) as.character(d$outcome_id) else rep(NA_character_, nrow(d))
  unique(paste0(d$gene, "~", ifelse(d$outcome_name == "eGFR", "eGFR", oc)))
}
# include the strict sensitivity (31) and grid (16) fits explicitly; they refit the same pairs/SNPs as 29 (so
# their palindromes are a subset already captured), but requiring their labels makes the completeness explicit.
expected <- unique(c(.fit_labels(file.path(D, "susie_coloc_results.csv")),
                     .fit_labels(file.path(D, "candidate_susie_egfr.csv")),
                     .fit_labels(file.path(D, "candidate_susie_nonegfr.csv")),
                     .fit_labels(file.path(D, "candidate_susie_sensitivity_egfr.csv")),
                     .fit_labels(file.path(D, "candidate_susie_sensitivity_nonegfr.csv"))))
# An empty obligation is never a legitimate state here: this gate only runs after the SuSiE chain,
# which by construction fitted at least one pair.
if (!length(expected))
  fail("the capture obligation is EMPTY -- no fitted pair was derived from the five SuSiE output CSVs. Coverage would pass without checking anything.")
.miss_cov <- setdiff(expected, seen_labels)
if (length(.miss_cov))
  fail(sprintf("INCOMPLETE capture: %d fitted pair(s) (>=10 SNPs, non-MHC) have NO seen marker this run -- a fit was skipped (resume/cache), so its palindromes are unverified. The capture run must force-fresh (clear cache/checkpoint). Missing: %s", length(.miss_cov), paste(head(.miss_cov, 8), collapse = "; ")))
if (length(expected)) cat(sprintf("[39] coverage OK: all %d fitted shipped pair(s) captured\n", length(expected)))

files <- list.files(SINK, pattern = "^palin_.*\\.csv$", full.names = TRUE)
sink_raw <- if (length(files)) rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE) else data.table()
if (nrow(sink_raw) && !"label" %in% names(sink_raw)) fail("palin rows lack a label column -- older helper; clear the sink and re-capture.")
# RECONCILE PER LABEL: seen markers count palindromes per fit; the palin rows for EACH label must total exactly
# the seen n_palindromic for that label. A global total can hide an extra row for label A cancelling a missing
# row for label B; per-label catches that so no palindrome goes unclassified (cold re-review Medium).
.seen_by  <- seen_dt[, .(seen_n = sum(n_palindromic)), by = label]
.palin_by <- if (nrow(sink_raw)) sink_raw[, .(palin_n = .N), by = label] else data.table(label = character(0), palin_n = integer(0))
.rec <- merge(.seen_by, .palin_by, by = "label", all = TRUE)
.rec[is.na(seen_n), seen_n := 0L]; .rec[is.na(palin_n), palin_n := 0L]
.rec_bad <- .rec[seen_n != palin_n]
if (nrow(.rec_bad))
  fail(sprintf("sink inconsistent (per-label): %d label(s) where seen n_palindromic != palin rows -- a palin_*.csv is missing/partial/duplicated: %s. Clear the sink and re-capture.",
               nrow(.rec_bad), paste(head(sprintf("%s(seen=%d,palin=%d)", .rec_bad$label, .rec_bad$seen_n, .rec_bad$palin_n), 5), collapse = "; ")))

if (nrow(sink_raw) == 0L) {   # zero palindromes AND per-label reconciled: capture provably ran, nothing to verify -> PASS
  writeLines(c(sprintf("run_id=%s", RUN_ID), sprintf("time=%s", format(Sys.time())),
               "n_palindromic=0", sprintf("n_fit_calls=%d", nrow(seen_dt)),
               # the obligation this PASS was measured against, so a later reader can tell a genuine
               # all-clear from a coverage check that had nothing to check
               sprintf("n_expected_fits=%d", length(expected)),
               "verdicts=NONE(zero palindromic fit-SNPs)"), PASS_MARK)
  cat(sprintf("\n[39 GATE PASS] zero palindromic fit-SNPs across %d captured fit-call(s) -- nothing to verify.\n  PASS marker: %s\n", nrow(seen_dt), PASS_MARK))
  quit(status = 0L, save = "no")
}
sink_dt <- sink_raw
req <- c("label", "rsid", "ea_pqtl", "oa_pqtl", "bim_a1", "z_pqtl", "z_out", "run_id")
if (any(!req %in% names(sink_dt))) fail(sprintf("palindrome rows missing columns: %s", paste(setdiff(req, names(sink_dt)), collapse = ", ")))
.bad <- unique(as.character(sink_dt$run_id[as.character(sink_dt$run_id) != RUN_ID]))
if (length(.bad)) fail(sprintf("palindrome rows carry other run_id(s) [%s] != current [%s] -- stale/mixed sink; clear it and re-run capture.", paste(head(.bad, 3), collapse = ","), RUN_ID))
sink_dt[, `:=`(ea_pqtl = toupper(ea_pqtl), oa_pqtl = toupper(oa_pqtl), bim_a1 = toupper(bim_a1),
               gene = sub("~.*", "", label))]
# severity = the WORSE of the exposure and outcome fit z at this SNP: the same (mis-signed) R feeds BOTH SuSiE
# fits, so a low-z-exposure / high-z-outcome palindrome still corrupts the outcome fit (cold re-review High).
sink_dt[, absz := pmax(abs(z_pqtl), abs(z_out), na.rm = TRUE)]
sink_dt <- sink_dt[order(-absz)][, .SD[1], by = .(label, rsid)]
cat(sprintf("[39] %d palindromic fit-SNP record(s) across %d fit(s), %d gene(s)\n",
            nrow(sink_dt), uniqueN(sink_dt$label), uniqueN(sink_dt$gene)))

# ── 2. pQTL EAF per (gene, rsid) via get_pqtl_region ─────────────────────────
if (!file.exists(INSTR)) fail(sprintf("instruments file not found (%s) -- needed for gene coords", INSTR))
instr <- readRDS(INSTR); setDT(instr)
rsid_map <- load_rsid_map()
if (!exists("CIS_WINDOW_KB")) fail("CIS_WINDOW_KB undefined (config) -- cannot reconstruct the pQTL region")

eaf_lu <- list()   # gene -> data.table(rsid, A1FREQ, ea_a1)
for (g in unique(sink_dt$gene)) {
  gi <- instr[gene == g][1]
  if (nrow(gi) == 0 || is.na(gi$gene_chrom) || is.na(gi$gene_start) || is.na(gi$gene_end)) {
    cat(sprintf("  [warn] no instrument coords for %s -- its palindromes will be UNRESOLVED\n", g)); next
  }
  rs <- max(1L, gi$gene_start - CIS_WINDOW_KB * 1000L); re <- gi$gene_end + CIS_WINDOW_KB * 1000L
  pd <- tryCatch(get_pqtl_region(g, as.character(gi$gene_chrom), rs, re, rsid_map),
                 error = function(e) { cat(sprintf("  [warn] get_pqtl_region(%s) failed: %s\n", g, conditionMessage(e))); NULL })
  if (is.null(pd) || nrow(pd) == 0 || !all(c("rsid", "A1FREQ", "ALLELE1") %in% names(pd))) next
  pd <- pd[!is.na(rsid) & !is.na(A1FREQ)]
  eaf_lu[[g]] <- unique(pd[, .(rsid, A1FREQ = as.numeric(A1FREQ), ea_a1 = toupper(ALLELE1))], by = "rsid")
}

get_eaf <- function(g, r) { t <- eaf_lu[[g]]; if (is.null(t)) return(list(NA_real_, NA_character_))
  h <- t[rsid == r]; if (!nrow(h)) return(list(NA_real_, NA_character_)); list(h$A1FREQ[1], h$ea_a1[1]) }

# ── 3. 1000G EUR freq of .bim A1 ─────────────────────────────────────────────
if (!file.exists(FRQ)) fail(sprintf("1000G freq-counts file not found (%s) -- run `plink --freq counts --keep-allele-order`", FRQ))
frq <- fread(FRQ, select = c("SNP", "A1", "A2", "C1", "C2"))
frq <- frq[SNP %in% unique(sink_dt$rsid)]
frq[, `:=`(A1 = toupper(A1), A2 = toupper(A2), freq_a1 = C1 / (C1 + C2))]
setkey(frq, SNP)

# ── 4. classify each palindromic fit-SNP ─────────────────────────────────────
verdict1 <- function(row) {
  ee <- get_eaf(row$gene, row$rsid); eaf <- ee[[1]]; ea_a1 <- ee[[2]]
  fr <- frq[SNP == row$rsid]
  if (is.na(eaf))                       return(list("UNRESOLVED", NA_real_, NA_real_, "no pQTL A1FREQ"))
  if (!nrow(fr))                        return(list("UNRESOLVED", eaf, NA_real_, "not in 1000G freq panel"))
  # the freq-panel A1 must be the sink's .bim A1 (both --keep-allele-order from the same panel); else orientation unknown
  if (!identical(fr$A1[1], row$bim_a1)) return(list("UNRESOLVED", eaf, NA_real_, sprintf("freq-panel A1 (%s) != sink bim_a1 (%s)", fr$A1[1], row$bim_a1)))
  # sink ea_pqtl must equal the sumstat ALLELE1 that A1FREQ refers to (consistency of the EAF's meaning)
  if (!is.na(ea_a1) && !identical(ea_a1, row$ea_pqtl)) return(list("UNRESOLVED", eaf, fr$freq_a1[1], sprintf("sumstat ALLELE1 (%s) != sink ea_pqtl (%s)", ea_a1, row$ea_pqtl)))
  f1000 <- fr$freq_a1[1]
  if (abs(eaf - 0.5) < FREQ_AMBIG)      return(list("AMBIGUOUS", eaf, f1000, "EAF within FREQ_AMBIG of 0.5"))
  d_al <- abs(f1000 - eaf); d_sw <- abs(f1000 - (1 - eaf))
  if (abs(d_al - d_sw) < FREQ_SEP)      return(list("AMBIGUOUS", eaf, f1000, sprintf("aligned/swapped near-tie (|f-EAF|=%.3f vs |f-(1-EAF)|=%.3f)", d_al, d_sw)))
  if (d_al < d_sw && d_al <= FREQ_TOL)  return(list("ALIGNED",  eaf, f1000, sprintf("|f-EAF|=%.3f", d_al)))
  if (d_sw < d_al && d_sw <= FREQ_TOL)  return(list("SWAPPED",  eaf, f1000, sprintf("|f-(1-EAF)|=%.3f", d_sw)))
  list("UNRESOLVED", eaf, f1000, sprintf("neither within tol (|f-EAF|=%.3f, |f-(1-EAF)|=%.3f)", d_al, d_sw))
}

# per-row classify via lapply+rbindlist (explicit columns; avoids the data.table `res[, -1]` numeric-j pitfall
# that would bind -1 instead of dropping a column -- cold re-review Medium)
res <- rbindlist(lapply(seq_len(nrow(sink_dt)), function(i) {
  v <- verdict1(sink_dt[i])
  data.table(verdict = v[[1]], pqtl_eaf = v[[2]], eur_freq_bimA1 = v[[3]], note = v[[4]])
}))
out <- cbind(sink_dt[, .(label, gene, rsid, ea_pqtl, oa_pqtl, bim_a1, z_pqtl, z_out, absz)], res)
setorder(out, verdict, -absz)
fwrite(out, OUT)

# ── 5. summarise, classify by MATERIALITY, BLOCK only on a MATERIAL swap ──────
# Editor ruling (journal-editor-review, 2026-07-18): the gate PRODUCES the audit and blocks ONLY on a genuine
# SWAPPED palindrome in a pair whose SuSiE result is actually USED (evaluable = NOT MHC-excluded). Palindromes
# that are near-0.5 (frequency cannot decide strand -- a documented limitation, NOT evidence of error) or
# UNRESOLVED (rsID-mapping gap / allele-label edge cases -- technical) are DISCLOSED, not blocked; and a SWAPPED
# in an MHC-excluded pair is immaterial (its SuSiE number is not_evaluable and unused). The measured audit here
# (aligned rate + named swaps) is what the strand/QC disclosure reports; a separate drop-palindrome sensitivity
# re-fit on the SuSiE-tier pairs hardens the "keep + disclose" decision.
tab <- out[, .N, by = verdict][order(-N)]
cat("\n[39] palindromic freq-gate verdicts:\n"); for (i in seq_len(nrow(tab))) cat(sprintf("  %-10s %d\n", tab$verdict[i], tab$N[i]))
cat(sprintf("[39] report: %s\n", OUT))

# headline QC: alignment rate over the frequency-RESOLVABLE palindromes (ALIGNED + SWAPPED)
n_al  <- out[verdict == "ALIGNED", .N]; n_res <- out[verdict %in% c("ALIGNED", "SWAPPED"), .N]
aln_pct <- 100 * n_al / max(1L, n_res)
cat(sprintf("[39] frequency-resolvable palindromes: %d ALIGNED / %d resolvable = %.4f%%\n", n_al, n_res, aln_pct))

# MHC-excluded pairs: their SuSiE result is not_evaluable and unused -> a swap there is IMMATERIAL (editor ruling)
mhc_lab <- character(0)
if (file.exists(file.path(D, "susie_coloc_results.csv"))) {
  .sc <- fread(file.path(D, "susie_coloc_results.csv"))
  if ("mhc_excluded" %in% names(.sc)) {
    .m <- .sc[mhc_excluded %in% c(TRUE, "TRUE")]
    mhc_lab <- unique(paste0(.m$gene, "~", ifelse(.m$outcome_name == "eGFR", "eGFR", as.character(.m$outcome_id))))
  }
}
swp <- out[verdict == "SWAPPED"]; if (nrow(swp)) swp[, material := !(label %in% mhc_lab)]
n_amb_hz <- out[verdict == "AMBIGUOUS" & absz >= Z_AMBIG, .N]
n_unres  <- out[verdict == "UNRESOLVED", .N]

if (nrow(swp)) { cat("\nSWAPPED palindromes (genuine mis-orientation; `material`=in a USED, non-MHC pair):\n")
                 print(swp[, .(label, rsid, z_pqtl, z_out, pqtl_eaf, eur_freq_bimA1, material)]) }
if (n_amb_hz)  cat(sprintf("\n[disclose] %d high-|z| AMBIGUOUS palindromes: EAF~0.5, frequency cannot decide strand (documented limitation, not mis-orientation).\n", n_amb_hz))
if (n_unres) { cat(sprintf("\n[disclose] %d UNRESOLVED palindromes (technical -- rsID-mapping gap / allele-label), by reason:\n", n_unres))
               print(out[verdict == "UNRESOLVED", .N, by = note][order(-N)]) }

n_material_swap <- if (nrow(swp)) swp[material == TRUE, .N] else 0L
if (n_material_swap)
  fail(sprintf("%d genuine SWAPPED palindrome(s) in an EVALUABLE (non-MHC) pair -- these enter a USED SuSiE fit with mis-signed LD and MUST be re-oriented/dropped and the pair re-fit: %s",
               n_material_swap, paste(head(swp[material == TRUE, sprintf("%s:%s", label, rsid)], 6), collapse = "; ")))

# PASS: no material swap. Record the full audit in the marker (source for the strand/QC disclosure numbers).
writeLines(c(sprintf("run_id=%s", RUN_ID), sprintf("time=%s", format(Sys.time())),
             sprintf("n_palindromic=%d", nrow(out)),
             sprintf("n_expected_fits=%d n_fit_calls=%d", length(expected), nrow(seen_dt)),
             sprintf("aligned=%d resolvable=%d aligned_pct=%.4f", n_al, n_res, aln_pct),
             sprintf("swapped_total=%d swapped_material=%d swapped_cleared_mhc=%d", nrow(swp), n_material_swap, nrow(swp) - n_material_swap),
             sprintf("high_z_ambiguous=%d unresolved=%d", n_amb_hz, n_unres),
             sprintf("verdicts=%s", paste(sprintf("%s:%d", tab$verdict, tab$N), collapse = ";"))), PASS_MARK)
cat(sprintf("\n[39 GATE PASS] %d/%d frequency-resolvable palindromes ALIGNED (%.4f%%); %d SWAPPED, all in MHC-excluded (unused) pairs -> 0 MATERIAL; %d high-|z| AMBIGUOUS + %d UNRESOLVED are disclosed limitations (editor ruling 2026-07-18).\n  PASS marker: %s\n",
            n_al, n_res, aln_pct, nrow(swp), n_amb_hz, n_unres, PASS_MARK))
