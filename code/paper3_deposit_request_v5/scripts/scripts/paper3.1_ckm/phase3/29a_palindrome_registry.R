# ============================================================================
# 29a: palindrome exclusion registry -- pre-fit application + post-gate append.
# ----------------------------------------------------------------------------
# WHY. 39_palindromic_freq_gate.R adjudicates every USED palindromic fit-SNP
# post hoc; a verdict=SWAPPED there costs a re-fit of the affected arm. That
# tuition must be paid AT MOST ONCE per SNP (PI order 2026-08-19): the verdict
# is persisted to the DATA-LAYER registry data/palindrome_exclusions/registry.csv
# (an input, not a result -- it survives regeneration), and every candidate-chain
# fit launched THROUGH this script applies it up front via the reviewed drop-hook
# (strand_align.R `ckm.palindrome_drop`, label-exact match), so a from-scratch
# rebuild flows straight into 29 with the exclusion already in force and 39
# passes first time. 39 stays BLOCKING: its job is now to catch NEW swaps only;
# each one is appended here (mode=append) and becomes permanently excluded.
# The registry is machine-derived from the gate's verdict artifact -- never a
# hand-typed list (editor consult 2026-08-19: APPROVE-WITH-MODIFICATIONS;
# thresholds stay single-sourced in 39, this script copies none of them).
#
# MODES (env CKM_29A_MODE):
#   append   derive verdict=="SWAPPED" rows from phase3_susie/palindromic_freq_gate.csv,
#            append any new (label,rsid) to the registry with provenance columns.
#            No fitting. Run after a 39 BLOCK. Idempotent (duplicates skipped).
#   launch   (default) read the registry, set the drop option, then source the
#            target fit script IN THIS R SESSION (options() cannot cross an
#            Rscript process boundary -- that is why this is a launcher, not a
#            side-file). Target: env CKM_29A_TARGET = 29 (default) | 31.
#            CAND_SUSIE_MODE / CKM_RUN_ID / CKM_PALINDROME_SINK pass through
#            untouched. Empty or absent registry -> plain launch (this wrapper
#            is the STANDARD candidate-chain path, not a remediation special case).
# Set CKM_29A_DRY=1 to print what would be appended/applied and exit before any
# write or fit (smoke hook).
#
# FAIL-CLOSED GUARDS (launch mode):
#   * a registry label NOT ending in "~eGFR" -> STOP. The non-eGFR engine
#     (15_susie_coloc.R) has no palindrome_dropped_rsids self-declaration column
#     yet (editor blocker #1, 2026-08-19); pre-applying there would be a silent
#     drop. Add the column (next-generation batch) before admitting such rows.
#   * a registry label present in the MAIN SuSiE table -> STOP. This script
#     launches candidate-chain fits only; a main-table hit means cascade/07(/08)
#     must re-run with the drop -- a PI decision, never an automatic action.
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
B    <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
D    <- file.path(B, "results", paste0("repro_", PHASE1_SET), "phase3_susie")
REG  <- file.path(B, "data/palindrome_exclusions/registry.csv")
MODE <- tolower(Sys.getenv("CKM_29A_MODE", "launch"))
DRY  <- nzchar(Sys.getenv("CKM_29A_DRY", ""))

read_reg <- function() if (file.exists(REG)) read.csv(REG, stringsAsFactors = FALSE) else NULL

if (MODE == "append") {
  gvf <- file.path(D, "palindromic_freq_gate.csv")
  if (!file.exists(gvf)) stop("29a append: gate verdict file absent: ", gvf)
  gv <- read.csv(gvf, stringsAsFactors = FALSE)
  keep_cols <- c("label", "gene", "rsid", "verdict", "ea_pqtl", "oa_pqtl", "bim_a1",
                 "pqtl_eaf", "eur_freq_bimA1", "note")
  miss <- setdiff(keep_cols, names(gv))
  if (length(miss)) stop("29a append: verdict file lacks column(s): ", paste(miss, collapse = ", "))
  sw  <- gv[gv$verdict == "SWAPPED", keep_cols, drop = FALSE]
  reg <- read_reg()
  new <- if (is.null(reg)) sw else sw[!paste(sw$label, sw$rsid) %in% paste(reg$label, reg$rsid), , drop = FALSE]
  if (!nrow(new)) { cat("29a append: no new SWAPPED verdicts -- registry already current.\n"); quit(save = "no", status = 0) }
  new$adjudicated <- format(Sys.Date())
  new$phase1_set  <- PHASE1_SET
  new$ckm_run_id  <- Sys.getenv("CKM_RUN_ID", "")
  new$source_file <- "phase3_susie/palindromic_freq_gate.csv"
  new$freq_panel  <- "1000G EUR (plink --freq counts --keep-allele-order)"
  cat("29a append: new SWAPPED verdict row(s) to persist:\n")
  print(new[, c("label", "rsid", "verdict", "pqtl_eaf", "eur_freq_bimA1", "note")], row.names = FALSE)
  if (DRY) { cat(sprintf("29a append DRY: would append %d row(s); registry untouched.\n", nrow(new))); quit(save = "no", status = 0) }
  dir.create(dirname(REG), recursive = TRUE, showWarnings = FALSE)
  out <- if (is.null(reg)) new else rbind(reg[, names(new), drop = FALSE], new)
  # atomic replacement (Task-1 F1): write the full file to a sibling temp path, then swap it in.
  # An interrupted write can no longer leave a truncated-but-readable registry that launch mode
  # would silently read as "no exclusions"; the worst interruption leaves the OLD registry intact.
  tmpf <- paste0(REG, ".tmp")
  write.csv(out, tmpf, row.names = FALSE)
  if (file.exists(REG) && !file.remove(REG)) stop("29a append: cannot replace registry: ", REG)
  if (!file.rename(tmpf, REG)) stop("29a append: rename failed: ", tmpf, " -> ", REG)
  cat(sprintf("29a append: %d row(s) appended -> %s (%d total)\n", nrow(new), REG, nrow(out)))
  quit(save = "no", status = 0)
}

if (MODE != "launch") stop("29a: CKM_29A_MODE must be 'append' or 'launch' (got '", MODE, "')")

TGT <- Sys.getenv("CKM_29A_TARGET", "29")
tgt_file <- if (TGT == "29") {
  file.path(B, "scripts/paper3.1_ckm/phase3/29_candidate_susie.R")
} else if (TGT == "31") {
  file.path(B, "scripts/paper3.1_ckm/phase3/31_candidate_susie_sensitivity.R")
} else stop("29a: CKM_29A_TARGET must be '29' or '31' (got '", TGT, "')")

reg <- read_reg()
if (!is.null(reg) && nrow(reg)) {
  if (!all(c("label", "rsid") %in% names(reg))) stop("29a: registry lacks label/rsid columns: ", REG)
  bad <- unique(reg$label[!grepl("~eGFR$", reg$label)])
  if (length(bad))
    stop("29a: registry has non-eGFR label(s) [", paste(bad, collapse = ", "),
         "] -- 15_susie_coloc.R has no palindrome_dropped_rsids self-declaration column yet ",
         "(editor blocker #1, 2026-08-19); add it before pre-applying non-eGFR exclusions.")
  mt <- file.path(D, "susie_coloc_results.csv")
  if (file.exists(mt)) {
    d  <- read.csv(mt, stringsAsFactors = FALSE)
    on <- if ("outcome_name" %in% names(d)) d$outcome_name else d$outcome
    lb <- paste0(d$gene, "~", ifelse(on == "eGFR", "eGFR", as.character(d$outcome_id)))
    hit <- intersect(unique(reg$label), unique(lb))
    if (length(hit))
      stop("29a: registry label(s) present in the MAIN SuSiE table: ", paste(hit, collapse = ", "),
           " -- re-fitting the main chain (cascade/07/08) with the drop is a PI decision; ",
           "this launcher covers the candidate chain only.")
  }
  options(ckm.palindrome_drop = reg[, c("label", "rsid")])
  cat("29a: registry exclusions in force (label-exact, via strand_align drop-hook):\n")
  print(reg[, c("label", "rsid", "adjudicated")], row.names = FALSE)
} else {
  cat("29a: registry empty/absent -- plain launch, no exclusions.\n")
}

if (DRY) {
  cat(sprintf("29a launch DRY: would source %s (CAND_SUSIE_MODE='%s'); stopping before the fit.\n",
              tgt_file, Sys.getenv("CAND_SUSIE_MODE", "")))
  quit(save = "no", status = 0)
}
source(tgt_file)
