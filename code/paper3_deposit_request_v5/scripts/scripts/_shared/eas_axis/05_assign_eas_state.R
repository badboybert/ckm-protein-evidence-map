# ============================================================================
# EAS AXIS 05: the four states (contract SS9 ONLY; SS10 item 5)
# ----------------------------------------------------------------------------
# Derives eas_outcome_state + the mutually exclusive state_reason from 03's
# estimates. This file must contain NO tier fields, NO F_lo, NO replication
# depth, and NO pQTL fields in any decision expression -- 07 verifies that by
# static scan. state_reason uses the contract SS9.1 vocabulary VERBATIM: the
# seven INCONCLUSIVE labels, the two significant labels, and
# NOT_EVALUABLE_<SS7.3 code>. No other label may appear.
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
suppressPackageStartupMessages(library(data.table))
B <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"); setwd(B)
EAS_DIR <- file.path("results", paste0("repro_", PHASE1_SET), "eas_axis")
est <- fread(file.path(EAS_DIR, "eas_outcome_estimates.csv"))
stopifnot(nrow(est) > 0L)

# SS7.3 codes are the only legal suffixes for NOT_EVALUABLE rows.
QC_CODES <- c("OUTCOME_UNAVAILABLE", "PHENOTYPE_NOT_ACCEPTABLE", "VARIANT_ABSENT",
              "BUILD_UNRESOLVED", "POSITION_CONFLICT", "ALLELE_MISMATCH",
              "MULTIALLELIC_OR_DUPLICATE", "PALINDROMIC_AMBIGUOUS", "LOW_INFO",
              "LOW_MAC", "LOW_STRATUM_COUNT", "SAMPLE_SIZE_UNRESOLVED", "INVALID_EFFECT")
stopifnot("05 needs 02/03's qc_reason to name NOT_EVALUABLE rows (SS9.1)" =
            "qc_reason" %in% names(est),
          all(est[qc_evaluable == FALSE, qc_reason] %in% QC_CODES))

# [DECISION-SCOPE-BEGIN] 07's static scan covers from here to the END marker:
# everything that builds a decision input (direction_class) or decides a state
# or a reason. The contract's excluded-field identifiers (SS9: the tier system,
# the winner's-curse bound, the replication and transferability sidecars) must
# not appear inside this fence; 07 names them on its side.
# EUR harmful-direction MR beta: KEYED join (not positional -- merge() re-sorts).
t4 <- fread(file.path("results", paste0("repro_", PHASE1_SET), "table_s_pair_evidence_tiered.csv"))
est[t4, on = .(gene, outcome), mr_beta_eur := i.mr_beta]
stopifnot(sum(is.na(est$mr_beta_eur)) == 0L)
est[, dir_prod := betaMR_EAS * mr_beta_eur]
est[, direction_class := fifelse(!is.finite(dir_prod), "NULL_OR_UNDEFINED",
                          fifelse(dir_prod > 0, "SAME",
                          fifelse(dir_prod < 0, "OPPOSITE", "NULL_OR_UNDEFINED")))]
est[, eas_direction_same := fifelse(direction_class == "SAME", TRUE,
                             fifelse(direction_class == "OPPOSITE", FALSE, NA))]

SIG <- 0.05
# SS9.1: fixed priority NOT_EVALUABLE -> q_eas<0.05 -> INCONCLUSIVE, and inside
# INCONCLUSIVE the seven-row table in contract order. The q<0.05 &
# NULL_OR_UNDEFINED cell is unreachable for a Wald estimate (beta 0 or NA never
# yields q<0.05); it falls to NULL_OR_UNDEFINED_DIRECTION rather than to an
# invented label.
assign1 <- function(evaluable, qc_code, dirc, q, qd, mt) {
  if (!isTRUE(evaluable)) return(c("EAS_NOT_EVALUABLE", paste0("NOT_EVALUABLE_", qc_code)))
  sig <- is.finite(q) && q < SIG
  if (sig && dirc == "SAME")     return(c("EAS_COMPATIBLE", "EAS_SIGNIFICANT_SAME_DIRECTION"))
  if (sig && dirc == "OPPOSITE") return(c("EAS_DISCORDANT", "EAS_SIGNIFICANT_OPPOSITE_DIRECTION"))
  r <- if (dirc == "NULL_OR_UNDEFINED") "NULL_OR_UNDEFINED_DIRECTION"
  else if (dirc == "SAME") {
    if (!isTRUE(mt))                  "SAME_DIR_EAS_NS_MAG_UNTESTABLE"
    else if (is.finite(qd) && qd < SIG) "SAME_DIR_EAS_NS_MAG_DIFFERENT"
    else                              "SAME_DIR_EAS_NS_NO_MAG_DIFFERENCE_DETECTED"
  } else {
    if (!isTRUE(mt))                  "OPPOSITE_DIR_EAS_NS_MAG_UNTESTABLE"
    else if (is.finite(qd) && qd < SIG) "OPPOSITE_DIR_EAS_NS_MAG_DIFFERENT"
    else                              "OPPOSITE_DIR_EAS_NS_NO_MAG_DIFFERENCE_DETECTED"
  }
  c("EAS_INCONCLUSIVE", r)
}
sr <- t(mapply(assign1, est$qc_evaluable, est$qc_reason, est$direction_class,
               est$q_eas, est$q_diff, est$magnitude_testable))
est[, eas_outcome_state := sr[, 1]]
est[, state_reason      := sr[, 2]]
est[, eas_state_short := c(EAS_COMPATIBLE = "EA-C", EAS_DISCORDANT = "EA-D",
                           EAS_INCONCLUSIVE = "EA-I", EAS_NOT_EVALUABLE = "EA-NE")[eas_outcome_state]]

# SS9.2 strict VARIANT_LEVEL sensitivity: the SAME partition function over
# q_eas_strict, for the counterfactual "only the strict family was admitted".
# Rows evaluable under primary but excluded by strictness get the explicit
# marker STRICT_EXCLUDED (pre-spec SS5) -- they are neither a state nor missing.
est[, strict_family := qc_evaluable == TRUE & sample_count_basis == "VARIANT_LEVEL"]
srS <- t(mapply(function(fam, evaluable, qc_code, dirc, q, qd, mt) {
  if (isTRUE(fam)) return(assign1(TRUE, qc_code, dirc, q, qd, mt))
  if (isTRUE(evaluable)) return(c("STRICT_EXCLUDED", "STRICT_EXCLUDED"))
  assign1(FALSE, qc_code, dirc, q, qd, mt)
}, est$strict_family, est$qc_evaluable, est$qc_reason, est$direction_class,
   est$q_eas_strict, est$q_diff, est$magnitude_testable))
est[, eas_outcome_state_strict := srS[, 1]]
est[, state_reason_strict      := srS[, 2]]
# [DECISION-SCOPE-END]

INCONCLUSIVE_REASONS <- c("NULL_OR_UNDEFINED_DIRECTION",
                          "SAME_DIR_EAS_NS_MAG_UNTESTABLE",
                          "SAME_DIR_EAS_NS_MAG_DIFFERENT",
                          "SAME_DIR_EAS_NS_NO_MAG_DIFFERENCE_DETECTED",
                          "OPPOSITE_DIR_EAS_NS_MAG_UNTESTABLE",
                          "OPPOSITE_DIR_EAS_NS_MAG_DIFFERENT",
                          "OPPOSITE_DIR_EAS_NS_NO_MAG_DIFFERENCE_DETECTED")

# partition checks (07 re-asserts; failing here writes nothing)
stopifnot(sum(is.na(est$eas_outcome_state)) == 0L,
          nrow(est[, .N, by = .(gene, outcome)][N > 1]) == 0L,
          all(est[eas_outcome_state == "EAS_NOT_EVALUABLE",
                  sub("^NOT_EVALUABLE_", "", state_reason)] %in% QC_CODES),
          all(est[eas_outcome_state == "EAS_COMPATIBLE",  state_reason] == "EAS_SIGNIFICANT_SAME_DIRECTION"),
          all(est[eas_outcome_state == "EAS_DISCORDANT",  state_reason] == "EAS_SIGNIFICANT_OPPOSITE_DIRECTION"),
          all(est[eas_outcome_state == "EAS_INCONCLUSIVE", state_reason] %in% INCONCLUSIVE_REASONS))

# SS10.1: the authoritative profile must carry the audit fields; a reduced table
# must not be written under the authoritative name (fail loud, not project away).
REQUIRED_PROFILE_COLS <- c(
  "gene", "outcome", "instrument_rsID", "eas_source", "eas_phenotype", "build",
  "snp_present", "allele_match", "info", "eaf", "source_qc_status",
  "n_total_study", "n_cases_study", "n_controls_study",
  "n_total_variant", "n_cases_variant", "n_controls_variant", "n_eff_variant",
  "n_total_mac_basis", "n_cases_mac_basis", "n_controls_mac_basis",
  "sample_count_basis",
  "mac_total_approx", "mac_cases_approx", "mac_controls_approx",
  "qc_evaluable", "qc_reason",
  "betaX_UKB", "seX_UKB", "F_original", "F_lo",
  "betaY_EAS_raw", "betaY_EAS_harmful", "seY_EAS",
  "betaMR_EAS", "seMR_EAS", "ci95_low_EAS", "ci95_high_EAS", "p_eas", "q_eas",
  "direction_class", "magnitude_testable", "magnitude_untestable_reason",
  "z_diff", "p_diff", "q_diff",
  "eas_outcome_state", "state_reason")
miss <- setdiff(REQUIRED_PROFILE_COLS, names(est))
if (length(miss))
  stop("05: eas_outcome_estimates.csv lacks SS10.1 profile columns: ",
       paste(miss, collapse = ", "),
       " -- fix 01/02/03 to carry them; the authoritative profile is not written reduced.")

stopifnot(all(est[strict_family == TRUE, eas_outcome_state_strict] %in%
                c("EAS_COMPATIBLE", "EAS_DISCORDANT", "EAS_INCONCLUSIVE")),
          all(est[qc_evaluable == TRUE & strict_family == FALSE,
                  eas_outcome_state_strict] == "STRICT_EXCLUDED"),
          all(est[qc_evaluable == FALSE, eas_outcome_state_strict] == "EAS_NOT_EVALUABLE"))

fwrite(est, file.path(EAS_DIR, "eas_outcome_profile.csv"))
cat("05 states:\n"); print(est[, .N, by = eas_outcome_state][order(-N)])
cat("strict sensitivity states (SS9.2):\n")
print(est[, .N, by = eas_outcome_state_strict][order(-N)])
chg <- est[strict_family == TRUE & eas_outcome_state != eas_outcome_state_strict,
           .(gene, outcome, eas_outcome_state, eas_outcome_state_strict)]
cat(sprintf("primary-vs-strict changed pairs: %d\n", nrow(chg))); if (nrow(chg)) print(chg)
cat("INCONCLUSIVE sub-reasons (contract SS9.1 order):\n")
print(est[eas_outcome_state == "EAS_INCONCLUSIVE",
          .N, by = factor(state_reason, levels = INCONCLUSIVE_REASONS)])
cat(sprintf("sign flips visible: %d OPPOSITE among evaluable (opposite-direction subtotal, SS9.2)\n",
            nrow(est[direction_class == "OPPOSITE" & qc_evaluable == TRUE])))
