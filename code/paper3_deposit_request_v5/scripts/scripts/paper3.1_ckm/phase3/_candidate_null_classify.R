# ============================================================================
# Candidate-mode NULL reason classifier (Stage 0 dashed layer).
# SOURCED BY (kept current -- a stale list here taught r45 that 31/32 were not consumers):
#   29_candidate_susie.R              -- classify_null() in the per-pair loop
#   30_candidate_edge_audit.R         -- CAND_BENIGN_TAGS, freeze-time gate
#   31_candidate_susie_sensitivity.R  -- classify_null() in the WORKERS; CAND_BENIGN_TAGS in the MAIN
#                                        process at the write boundary. BOTH sources are required: a
#                                        worker-only source left the main process without the symbol and
#                                        killed the script after every edge had run (r23, and again r44).
#   32_dashed_finalize.R              -- CAND_BENIGN_TAGS: `dropped:collapsed` asserts a BENIGN collapse,
#                                        so the consumer checks the reason against this list (r43).
#   _test_candidate_null_classify.R   -- unit test
# phase3/33 hashes this file and lists it as a dependency of every artifact the four above produce.
#
# The reused SuSiE engines (phase3/15, phase3/26) return R NULL for BOTH:
#   (benign) biological insufficiency  -- too few overlapping SNPs / no computable LD / empty
#            credible set / empty coloc.susie summary (a legitimate "no rescue possible here"); AND
#   (FATAL)  genuine computation/data failure -- coloc.susie() threw; pQTL region missing required
#            columns; an MR-passing candidate whose pQTL/outcome region reads as no-data.
# Only the emitted log text distinguishes them. This file is the SINGLE SOURCE for that decision,
# with DEFAULT-FATAL on any UNRECOGNISED reason, so string drift can never silently demote a real
# failure to a benign "uncomputable" skip (editor middle-path; Codex round-1 High #1).
#
# classify_null(log_lines) -> list(status = "benign" | "FATAL", tag = <clean tag>, detail = <text>)
#   "benign" -> 29 records run_status="skip_benign", error_msg=tag (a true no-rescue).
#   "FATAL"  -> 29 ABORTS (checkpoint-excludes the pair -> retried on resume). NEVER frozen.
# CAND_BENIGN_TAGS is the SINGLE SOURCE for "which collapses are benign", and therefore for whether a
# `dropped:collapsed` edge is a biological negative or a computation failure wearing that label. Every
# consumer sources this file; nobody copies the list.
# ============================================================================

# type-(b) computation/data failures -> must ABORT. The two engines emit DIFFERENT coloc strings:
#   15 -> message("    [coloc.susie error] ..."), 26 -> message("    [coloc.susie] ..."); both also
#   cat("[skip] coloc.susie failed"). The bracketed prefix "[coloc.susie" matches BOTH messages and
#   NOT the benign "[skip] coloc.susie returned empty summary" (whose bracket precedes "skip").
CAND_FATAL_SIG <- "\\[coloc\\.susie|coloc\\.susie failed|pQTL missing columns|no pQTL data"

# benign no-rescue signatures -> clean tags (no ambiguous substrings that could trip 30's gate).
cand_benign_tag <- function(log_str) {
  # STRAND-ALIGN DROP IS THE SAME GATE, ONE STEP LATER. strand_align_R() runs AFTER the LD
  # alignment and removes strand-inconsistent SNPs from the same set, so "< 10 after strand-align
  # drop" means exactly what "< 10 SNPs after LD" means: too few SNPs at this locus to fit. It was
  # simply not listed, so it fell to the default-FATAL arm and would abort 29 on a legitimate
  # low-SNP locus. All three engines emit it (15:409, 16:370, 26:456), in two wordings.
  if (grepl("< 5 rsIDs|< 10 rsIDs|< 10 SNPs after LD|< 10 after LD|< 10 SNPs after strand-align drop|< 10 after strand-align drop|Only [0-9]+ SNPs matched|Only [0-9]+ matched LD", log_str)) return("few_snps_for_ld")
  if (grepl("\\[LD local\\] Failed|LD matrix failed|\\[skip\\] LD failed", log_str)) return("ld_unavailable")
  if (grepl("pQTL only", log_str)) return("few_pqtl_snps")
  if (grepl("< 20 after merge|< 20 merged|< 20 after harmonise", log_str)) return("few_overlap_snps")
  if (grepl("susie_rss failed for|susie_rss exp failed|susie_rss out failed", log_str)) return("no_credible_set")
  if (grepl("empty summary", log_str)) return("empty_coloc_summary")
  NA_character_
}
CAND_BENIGN_TAGS <- c("few_snps_for_ld", "ld_unavailable", "few_pqtl_snps",
                      "few_overlap_snps", "no_credible_set", "empty_coloc_summary")

classify_null <- function(log_lines) {
  nz         <- trimws(log_lines[nzchar(trimws(log_lines))])
  detail_all <- if (length(nz)) paste(nz, collapse = " | ") else "(no engine output captured)"
  log_str    <- paste(log_lines, collapse = " ")
  # FATAL takes priority: a fatal signature anywhere => abort, regardless of any benign line present.
  if (grepl(CAND_FATAL_SIG, log_str))
    return(list(status = "FATAL", tag = "computation_or_data_failure", detail = detail_all))
  bt <- cand_benign_tag(log_str)
  if (!is.na(bt)) return(list(status = "benign", tag = bt, detail = bt))
  # DEFAULT-FATAL on unknown: an unrecognised NULL reason must never be frozen as a benign skip.
  list(status = "FATAL", tag = "unrecognized_null", detail = detail_all)
}
