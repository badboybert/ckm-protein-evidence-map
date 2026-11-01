# Phase 3-S: Coloc Prior Sensitivity (p12 sweep)
# =========================================================
# For the strong coloc hits (dual rule PP.H4>=0.80 & PP.H4/PP.H3>=5), re-run coloc.abf
# with p12 in {1e-6, 5e-6, 1e-5} and compare PP.H4 stability.
#
# Method:
#   Sources coloc.R (COLOC_SOURCE_ONLY=TRUE, no main() call)
#   Overrides run_coloc_abf() to accept a custom p12 value
#   Re-runs run_coloc_pair() for each strong hit × 3 p12 values
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/scripts/phase3s_prior_sweep.R"

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
suppressPackageStartupMessages(library(data.table))

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))   # for is_coloc_pass() (explicit direct source; coloc also loads it)
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_prior_bands.R"))   # is_prior_stable()/prior_band(): ONE definition, deliberately NOT in 00_config.R (coloc inputs key)
# ── Source coloc.R without running main() ──────────────────────────────────
COLOC_SOURCE_ONLY <- TRUE
suppressMessages(source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R")))
rm(COLOC_SOURCE_ONLY)

cat("\n=== Phase 3-S: Coloc Prior Sensitivity Sweep ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

if (!exists("PSW_OUT")) PSW_OUT <- file.path(RESULTS_DIR, "phase3_sensitivity")   # override for corrected cascade
out_dir <- PSW_OUT
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# WRITING THE NON-eGFR TABLE OVER A TABLE THAT ALREADY HAS eGFR ROWS ERASES THEM. The documented
# order is 10 -> 11, and 11 appends its 16 eGFR rows to these same canonical filenames. Nothing
# enforced the order: run both again concurrently and, if 11 finishes its 16 pairs first, the save
# at the end replaces the 76-row table with 60 non-eGFR rows. Both processes report success and the
# result is a fresh, valid-looking file with no record that 16 pairs were erased. eGFR rows here
# mean 11 has already run (this is out of order) or is running now -- both are stops.
#
# CHECKED BEFORE THE SWEEP, NOT BEFORE THE SAVE. Placed at the save it still aborted, but only after
# computing all 60 pairs -- an hour spent to reach a stop that is knowable in milliseconds.
# It does not close the concurrent-write window, which needs a lock; it stops the ordering mistake.
local({
  .w <- file.path(out_dir, "prior_sweep_wide.rds")
  if (file.exists(.w)) {
    .old <- tryCatch(as.data.table(readRDS(.w)), error = function(e) NULL)
    if (is.null(.old))
      stop(sprintf("10: %s exists but cannot be read -- it is half-written or corrupt. Delete it for a clean 10 -> 11 sequence rather than sweeping into an unknown state.", basename(.w)))
    if ("outcome_name" %in% names(.old) && any(.old$outcome_name == "eGFR"))
      stop(sprintf("10: %s already holds %d eGFR row(s). This script writes the non-eGFR half only, so its save would erase them. Run 10 BEFORE 11, or delete the file for a clean 10 -> 11 sequence.",
                   basename(.w), sum(.old$outcome_name == "eGFR")))
  }
})

# ── Override run_coloc_abf to use global CURRENT_P12 ─────────────────────────
# After source(), run_coloc_pair() calls run_coloc_abf() and looks it up
# in the global env. Re-assigning here overrides the function for all callers.
CURRENT_P12 <- 1e-5  # default; modified in loop below

run_coloc_abf <- function(pqtl_beta, pqtl_varbeta, n_pqtl,
                           out_beta, out_varbeta, meta, snp_ids) {
  if (length(pqtl_beta) < 20) { cat("    [skip] SNP < 20\n"); return(NULL) }
  ds1 <- list(beta = pqtl_beta, varbeta = pqtl_varbeta,
              N = n_pqtl, sdY = 1, type = "quant", snp = snp_ids)
  ds2 <- list(beta = out_beta, varbeta = out_varbeta,
              N = meta$n, type = meta$type, snp = snp_ids)
  if (meta$type == "cc") ds2$s <- meta$ncase / meta$n else ds2$sdY <- meta$sdY   # the outcome's own phenotype SD (GWAS_META), not a literal 1: coloc's quant prior is 0.15*sdY and eGFR is on ln scale
  tryCatch(
    coloc::coloc.abf(dataset1 = ds1, dataset2 = ds2, p12 = CURRENT_P12),
    error = function(e) fail_loud_data(sprintf("prior-sweep coloc.abf FAILED (data problem, NOT a legit skip): %s", conditionMessage(e)))  # 2026-07-01 (closure B1): fail-loud (classed) so the outer per-pair catch can't demote it to a NA cell
  )
}
cat("run_coloc_abf overridden to use CURRENT_P12\n\n")

# ── Load the strong hits ───────────────────────────────────────────────────────
if (!exists("PSW_COLOC")) PSW_COLOC <- file.path(PHASE3_DIR, "coloc_results_v7_stage1.rds")   # override for corrected (v6)
if (!exists("PSW_INSTR")) PSW_INSTR <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")
coloc_orig <- readRDS(PSW_COLOC)
instr      <- readRDS(PSW_INSTR)
setDT(coloc_orig); setDT(instr)

strong <- coloc_orig[is_coloc_pass(PP.H4, PP.H3)][order(-PP.H4)]   # dual coloc rule (D1 2026-07-04; was H4-only -> 53)
.egfr_id <- if (exists("OUTCOMES") && !is.null(OUTCOMES[["eGFR"]])) OUTCOMES[["eGFR"]] else { warning("[eGFR-id] OUTCOMES$eGFR missing; legacy fallback ebi-a-GCST003372"); "ebi-a-GCST003372" }
if (exists("PSW_NONEGFR_ONLY") && isTRUE(PSW_NONEGFR_ONLY)) strong <- strong[outcome_id != .egfr_id]  # eGFR handled by 25 (id from config)
if (exists("PSW_LIMIT")) strong <- strong[seq_len(min(PSW_LIMIT, nrow(strong)))]   # smoke-test cap (scope only, not results)
cat(sprintf("Strong hits to sweep: %d pairs\n\n", nrow(strong)))

# --- optional pair-chunk parallelism (PSW_CHUNK_N>0 = worker; PSW_MERGE=1 = merge) ---
.PCN <- as.integer(Sys.getenv("PSW_CHUNK_N","0")); .PCI <- as.integer(Sys.getenv("PSW_CHUNK_I","0"))
.PMERGE <- Sys.getenv("PSW_MERGE","") %in% c("1","true","TRUE")
.PCHUNK <- .PCN > 0 && !.PMERGE
if (.PCHUNK) { .psz <- ceiling(nrow(strong)/.PCN); .plo <- .PCI*.psz + 1; .phi <- min((.PCI+1)*.psz, nrow(strong))
  strong <- if (.plo <= nrow(strong)) strong[.plo:.phi] else strong[0]
  cat(sprintf("PSW CHUNK %d/%d: %d pairs\n", .PCI, .PCN, nrow(strong))) }
if (.PMERGE) strong <- strong[0]   # merge mode: skip the p12 loop; chunk results loaded below

rsid_map <- load_rsid_map()

# ── p12 sweep ─────────────────────────────────────────────────────────────────
p12_values <- c(1e-6, 5e-6, 1e-5)
all_sweep  <- list()

for (p12_val in p12_values) {
  CURRENT_P12 <<- p12_val
  cat(sprintf("\n=== p12 = %.0e ===\n", p12_val))
  p12_results <- list()

  for (i in seq_len(nrow(strong))) {
    g   <- strong$gene[i]
    oid <- strong$outcome_id[i]

    gene_info <- instr[gene == g][1]
    if (nrow(gene_info) == 0) {
      cat(sprintf("  [skip] %s: no instrument\n", g)); next
    }

    cat(sprintf("  [%d/%d] %s ~ %s\n", i, nrow(strong), g, strong$outcome_name[i]))
    flush(stdout())

    res <- tryCatch(
      run_coloc_pair(g, oid,
                     gene_info$gene_chrom, gene_info$gene_start, gene_info$gene_end,
                     rsid_map),
      error = function(e) {
        # 2026-07-01 (closure B1): ANY error from a strong-pair sweep call is effectively fatal (a strong
        # pair must yield a complete grid; 99 D2 rejects any NA anyway) -> flag it, whether or not it came
        # via fail_loud_data(). Legit-empty skips return NULL WITHOUT erroring, so they are unaffected.
        message(sprintf("    [ERROR] %s ~ %s: %s", g, oid, e$message))
        assign(".CKM_DATA_FAILURE", sprintf("%s ~ %s: %s", g, oid, conditionMessage(e)), envir = .GlobalEnv); NULL
      }
    )
    if (exists(".CKM_DATA_FAILURE", envir = .GlobalEnv))  # abort now, do NOT continue with a NA row
      stop(paste("10 ABORT -- data-foundation failure during prior sweep:", get(".CKM_DATA_FAILURE", envir = .GlobalEnv)))

    if (!is.null(res)) {
      res[, p12 := p12_val]
      p12_results[[paste(g, oid)]] <- res
    } else {
      # P3-F2 (cold audit): run_coloc_pair returned NULL for a STRONG pair -> keep the tracked NA row but emit a
      # LOUD marker so it cannot silently demote/undercount the pair (99_completeness_audit also asserts no NA cells).
      cat(sprintf("  [P3-F2 WARN] prior-sweep run_coloc_pair NULL for STRONG pair %s ~ %s at p12=%.1e -> NA row (investigate; not a silent skip)\n", g, oid, p12_val))
      p12_results[[paste(g, oid)]] <- data.table(
        gene = g, outcome_id = oid, p12 = p12_val,
        PP.H4 = NA_real_, PP.H3 = NA_real_
      )
    }
    Sys.sleep(0.1)
  }

  sweep_dt <- rbindlist(p12_results, fill = TRUE)
  all_sweep[[as.character(p12_val)]] <- sweep_dt
  cat(sprintf("  p12=%.0e: %d/%d pairs completed, PP.H4>=0.80: %d\n",
              p12_val,
              sum(!is.na(sweep_dt$PP.H4)),
              nrow(strong),
              sum(sweep_dt$PP.H4 >= 0.80, na.rm = TRUE)))
}

# ── Combine and build comparison table ────────────────────────────────────────
sweep_all <- rbindlist(all_sweep, fill = TRUE)
if (.PCHUNK) {
  saveRDS(sweep_all, file.path(out_dir, sprintf("prior_sweep_long_chunk%02d.rds", .PCI)))
  cat(sprintf("PSW CHUNK %d done: %d rows -> chunk%02d\n", .PCI, nrow(sweep_all), .PCI)); quit(save = "no")
}
if (.PMERGE) {
  cf <- sort(Sys.glob(file.path(out_dir, "prior_sweep_long_chunk*.rds")))
  sweep_all <- rbindlist(lapply(cf, readRDS), fill = TRUE)
  cat(sprintf("PSW MERGE: %d chunk files -> %d rows\n", length(cf), nrow(sweep_all)))
}
sweep_all[, outcome_name := names(OUTCOMES)[match(outcome_id, OUTCOMES)]]

# Pivot wider
sweep_wide <- dcast(
  # FULL PRECISION. These posteriors were rounded to four decimals here, before dcast -- so
  # `stable` and `sensitive` below were computed from the rounded values, and the wide table
  # SAVED them rounded while the eGFR engine writes the same columns unrounded. The combined
  # table therefore carried outcome-dependent precision (FGF5~AF's baseline 0.996730757478863
  # was stored as 0.9967), and a raw value inside the rounding interval of 0.80 or 0.20 would
  # flip a classification. Rounding belongs in the printed summary, not in the values a
  # threshold is applied to or the file downstream reads.
  sweep_all[, .(gene, outcome_id, outcome_name, p12, PP.H4, PP.H3)],
  gene + outcome_id + outcome_name ~ p12,
  value.var = c("PP.H4", "PP.H3")
)
setnames(sweep_wide,
  c("PP.H4_1e-06", "PP.H4_5e-06", "PP.H4_1e-05"),
  c("PP.H4_p12_1e6", "PP.H4_p12_5e6", "PP.H4_p12_1e5"),
  skip_absent = TRUE
)
# Add original result
orig_sub <- coloc_orig[is_coloc_pass(PP.H4, PP.H3),
                        .(gene, outcome_id, PP.H4_original = PP.H4)]   # dual coloc rule (D1 2026-07-04)
sweep_wide <- merge(sweep_wide, orig_sub, by = c("gene","outcome_id"), all.x = TRUE)

# Stability classification
if ("PP.H4_p12_1e6" %in% names(sweep_wide) &&
    "PP.H4_p12_1e5" %in% names(sweep_wide)) {
  # SINGLE SOURCE OF TRUTH (00_config.R is_prior_stable / prior_band; editor ruling 2026-08-01).
  # The old local rule ANDed p12 = 1e-5, which is coloc.abf's default and therefore vacuous.
  # Verified equivalent on this generation before the change: 23 pairs either way.
  sweep_wide[, stable := is_prior_stable(PP.H4_p12_1e6)]
  sweep_wide[, prior_band := prior_band(PP.H4_p12_1e6)]
  # RENAMED from `sensitive`: that name collided with the figure's "Prior-sensitive" BAND while
  # meaning something else (any non-stable pair, or a >0.20 shift between priors). Nothing consumes
  # this column; it stays as a diagnostic under a name that cannot be mistaken for the band.
  sweep_wide[, prior_shift_gt_020 := !stable | abs(PP.H4_p12_1e6 - PP.H4_p12_1e5) > 0.20]
}

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(sweep_all,  file.path(out_dir, "prior_sweep_long.rds"))
saveRDS(sweep_wide, file.path(out_dir, "prior_sweep_wide.rds"))
fwrite(sweep_wide,  file.path(out_dir, "prior_sweep_wide.csv"))
cat(sprintf("\nSweep results saved to: %s\n", out_dir))

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== Prior Sensitivity Summary ===\n")
cat(sprintf("Pairs tested:             %d\n", nrow(strong)))
cat(sprintf("Prior-robust (PP.H4>=0.80 at p12=1e-6; is_prior_stable): %d\n",
    sum(sweep_wide$stable, na.rm = TRUE)))
# THE LABEL NOW STATES THE DEFINITION IT COUNTS. `sensitive` is the OR of two conditions -- not
# stable at either prior, OR a change above 0.20 -- but this line read "PP.H4 changes >0.20", which
# is only the second. On the current 60 non-eGFR pairs it prints 42 while 39 actually move by more
# than 0.20: ABO~T2D (0.1870), FN1~CAD (0.1764) and IL6R~CAD (0.1787) are counted because they fall
# below 0.80, not because they moved. The tier-relevant `stable` flag is unaffected; what was wrong
# was the stated meaning of a published count. Both components are printed so neither can be read
# as the other.
cat(sprintf("Sensitive (NOT stable OR |dPP.H4| >0.20):   %d\n",
    sum(sweep_wide$prior_shift_gt_020, na.rm = TRUE)))
cat(sprintf("  of which |dPP.H4| >0.20:                  %d\n",
    sum(abs(sweep_wide$PP.H4_p12_1e6 - sweep_wide$PP.H4_p12_1e5) > 0.20, na.rm = TRUE)))
cat(sprintf("  of which not stable but |dPP.H4| <=0.20:  %d\n",
    sum(!sweep_wide$stable & abs(sweep_wide$PP.H4_p12_1e6 - sweep_wide$PP.H4_p12_1e5) <= 0.20,
        na.rm = TRUE)))

cat("\nFull comparison table:\n")
print_cols <- intersect(c("gene","outcome_name","PP.H4_original",
                           "PP.H4_p12_1e6","PP.H4_p12_5e6","PP.H4_p12_1e5",
                           "stable","prior_band","prior_shift_gt_020"),
                        names(sweep_wide))
print(sweep_wide[order(-sweep_wide$PP.H4_original), ..print_cols])
