# ============================================================================
# 33 (Stage 0 · provenance): write a self-documenting sidecar for the B1 sensitivity re-run + finalization,
# recording engine/script/input/output md5s, the in-memory patch spec (asserted 2/2 sites on the FROZEN
# engines), settings, and R/package versions (editor Q6 m-3).
#
# WHAT IT DOES, STATED HONESTLY (r12/r13): it RECORDS the md5 of the code and artifacts present WHEN IT
# RUNS, and CHECKS that each artifact postdates its LISTED code and LISTED artifact dependencies. It does NOT
# prove which code produced the bytes -- an md5 of the current file plus an mtime comparison cannot do that, and the
# gate below says so where it is defined. The original header claimed this sidecar "pins exactly WHICH
# engine version + patch produced them"; it does not, and 29's checkpoint hole (found r12) was a live
# demonstration of the gap. Read it as: the code present at certification time, plus a staleness tripwire
# -- not a proof of lineage. See also the "what this sidecar does NOT pin" section at the end.
# NON-DESTRUCTIVE: reads only; writes ONE new file candidate_susie_sensitivity_provenance.md.
# Run: cd "<E>/20260331-CKM/Paper 3" && Rscript scripts/paper3.1_ckm/phase3/33_sensitivity_provenance.R
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
Sys.setenv(CKM_OFFLINE = "1")
# base-R only (readLines / tools::md5sum / vapply / writeLines) -- no data.table needed. Package versions are
# recorded below via packageVersion(), which does not require the packages to be attached (Codex review Low).
B  <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
P  <- file.path(B, "scripts/paper3.1_ckm/phase3")
# _phase1_set.R defines PHASE1_SET, not PHASE1_TAG; the old reference aborted every launch
# before the CKM_RUN_ID check this script exists to enforce.
RC <- file.path(B, file.path("results", paste0("repro_", PHASE1_SET))); D <- file.path(RC, "phase3_susie")
OUT <- file.path(D, "candidate_susie_sensitivity_provenance.md")

md5 <- function(f) if (file.exists(f)) unname(tools::md5sum(f)) else NA_character_
ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "NA")

# assert the A" 2-site patch still applies to the FROZEN engines, so this provenance is not stale
pat <- "min_abs_corr = 0, verbose = FALSE, check_LD = FALSE"
sites <- function(f) if (file.exists(f)) sum(grepl(pat, readLines(f), fixed = TRUE)) else NA_integer_
n15 <- sites(file.path(P, "15_susie_coloc.R")); n26 <- sites(file.path(P, "26_susie_coloc_egfr.R"))
if (is.na(n15) || is.na(n26) || n15 != 2L || n26 != 2L)
  stop(sprintf("33: engine patch-site count != 2 (15=%s, 26=%s) -- the frozen engines changed; the sensitivity provenance would be stale. Re-run 31/32 against the current engines before recording provenance.", n15, n26))

# 2026-07-15: 28, 29 and 30 were ADDED (with R/00_config.R and R/coloc.R below). Every one of them can
# change an artifact this sidecar certifies -- 28 builds the pool and now applies the primary FDR rule,
# 29 produces the candidate SuSiE runs, 30 produces the audit that 32's dashed set is selected from.
# A sidecar that hashes the artifacts while omitting the code that produced them certifies the right
# bytes under an unknown rule, which is the exact failure this whole fix is about.
# 33 hashes ITSELF: the sidecar must identify the generator + gate version that wrote it, or a later
# edit to these gates leaves an old record indistinguishable from one the current gates passed (r11).
scripts <- c("33_sensitivity_provenance.R",
             "15_susie_coloc.R", "26_susie_coloc_egfr.R", "_candidate_null_classify.R",
             "28_candidate_edge_enumerate.R",   # builds the pool (and now applies the primary FDR rule)
             "29_candidate_susie.R",            # produces candidate_susie_{egfr,nonegfr}.csv
             "30_candidate_edge_audit.R",       # produces candidate_edge_audit.csv -> 32's dashed set
             "31_candidate_susie_sensitivity.R", "32_dashed_finalize.R")
# Logic that lives outside P but can change what this sidecar certifies. 00_config carries
# primary_fdr_set() (the FDR rule 28/32 apply); R/coloc.R is sourced by 29 and 31 on their NON-eGFR
# branch only (the eGFR path evals 26's head, which sources 00_config alone) and defines the engine they
# run. Both are hashed unconditionally here -- the sidecar records what the run could reach; the
# dependency gate above is where the eGFR/non-eGFR split matters. Hashing the artifacts while omitting the code that produced them certifies the
# right bytes under an unknown rule -- which is the exact failure this fix exists to remove.
EXTRA_LOGIC <- file.path(B, c("R/00_config.R", "R/coloc.R", "R/strand_align.R"))   # strand_align.R is sourced by both engine heads (LD sign-alignment); hash it too (cold re-review Low)
# 2026-07-15: `mr_significant_fdr.rds` was REPLACED by `mr_all_results.rds`. That artifact is the saved
# by-method one (see the three-surface block in R/00_config.R -- do not restate it here); 28/32 no longer read it
# -- they derive the manuscript's PRIMARY FDR from mr_all_results.rds via primary_fdr_set(). A sidecar
# hashing an input the scripts do not read certifies the wrong lineage, which is worse than none.
# 28's own raw inputs were also ADDED: hashing only the candidate_pool.csv it PRODUCED proves what the
# pool is, not what produced it -- so a changed coloc/Bonferroni/instrument input would leave this
# sidecar certifying a pool it can no longer account for.
inputs  <- c(file.path(D, c("candidate_edge_audit.csv", "candidate_pool.csv",
                            "candidate_susie_egfr.csv", "candidate_susie_nonegfr.csv")),
             file.path(RC, c("mr_all_results.rds", "coloc_results_v7.rds",
                             "mr_significant_bonferroni.rds")),
             PHASE1_INSTRUMENTS)
# candidate_pool_no_coloc.csv is 28's exclusion-audit surface: which FDR-significant pairs had no
# coloc row at all. Narrowing the FDR gate changes who reaches it, so leaving it unpinned would let that
# surface drift unnoticed. 28 now always writes it (empty if none), so it is safe to require here.
outputs <- file.path(D, c("candidate_susie_sensitivity_egfr.csv", "candidate_susie_sensitivity_nonegfr.csv",
                          "candidate_pool_no_coloc.csv", "dashed_edges_final.csv"))
for (.f in c(file.path(P, scripts), EXTRA_LOGIC, inputs, outputs))
  if (!file.exists(.f)) stop(sprintf("33: expected file missing -- run the chain first: 28 -> 29 (both modes) -> 30 -> 31 (both modes) -> 32, then this. Missing: %s", .f))

# PALINDROME GATE (39) must have PASSED for this run before the dashed layer is certified.
#
# WHAT PASS MEANS, precisely, because this comment used to overstate it. Palindromic SNPs are left
# UNCHANGED in the SuSiE fits (strand_align.R). 39 audits every one of them by frequency and BLOCKS on
# a genuine SWAPPED palindrome in an evaluable (non-MHC) pair. It deliberately does NOT block on
# AMBIGUOUS (EAF near 0.5, where frequency cannot decide strand even in principle) or on UNRESOLVED
# (rsID-mapping / allele-label gaps): those are DISCLOSED, per the journal-editor ruling of
# 2026-07-18, with the drop-palindrome sensitivity re-fit (41) supporting the keep-and-disclose call.
#
# So PASS certifies, precisely: no SWAP was DETECTED among the palindromes whose strand frequency
# could resolve, in a used fit, and the ones it could not resolve are on the record. It does NOT
# certify that every palindrome is same-strand -- for the AMBIGUOUS ones no method here can say. The stronger reading was written here
# and is wrong: DPEP1~AF's fit window carries rs154656 at pQTL EAF 0.420 / EUR 0.476 with |z| 94.3,
# an AMBIGUOUS SNP that stays in the fit and does not stop PASS. The PP.H4 for such a pair is computed
# with a strand-undecidable SNP, which is a stated limitation of the method, not a certification of it.
.pass_mark <- file.path(D, "palindromic_freq_gate.PASS")
if (!file.exists(.pass_mark))
  stop("33: palindromic_freq_gate.PASS missing -- 39_palindromic_freq_gate.R has not passed for this run. Run the SuSiE/dashed chain with the capture sink set (CKM_PALINDROME_SINK + CKM_RUN_ID), then 39, before certifying the dashed layer.")
.cur_rid <- Sys.getenv("CKM_RUN_ID", "")
if (!nzchar(.cur_rid))
  stop("33: CKM_RUN_ID not set -- cannot verify the palindrome gate PASS marker is from THIS run; a stale marker could otherwise satisfy this check. Export CKM_RUN_ID for the certified run.")
.pm_rid <- sub("^run_id=", "", grep("^run_id=", readLines(.pass_mark, warn = FALSE), value = TRUE)[1])
if (is.na(.pm_rid) || !identical(.pm_rid, .cur_rid))   # NA when the marker has no run_id= line (malformed)
  stop(sprintf("33: palindromic_freq_gate.PASS run_id (%s) != current CKM_RUN_ID (%s) -- the gate passed on a DIFFERENT run (or the marker is malformed); re-run 39 for this run.", .pm_rid, .cur_rid))
# The gate must have run AFTER the fits it certifies: a PASS marker OLDER than any SuSiE fit output means a
# fit (29/31) was re-run after the gate, leaving those fits unverified (cold re-review High).
.fit_out <- file.path(D, c("susie_coloc_results.csv",   # the 52-solid tier fits (15/26) are gate-covered too
                           "candidate_susie_egfr.csv", "candidate_susie_nonegfr.csv",
                           "candidate_susie_sensitivity_egfr.csv", "candidate_susie_sensitivity_nonegfr.csv"))
.fit_out <- .fit_out[file.exists(.fit_out)]
if (length(.fit_out) && file.mtime(.pass_mark) < max(file.mtime(.fit_out)))
  stop("33: palindromic_freq_gate.PASS is OLDER than a SuSiE fit output (susie_coloc_results.csv / candidate_susie_*.csv) -- a fit was re-run after the gate; re-run 39 so the gate covers the current fits before certifying.")

# ── STALENESS GATE: existence is NOT enough (r8/r9/r10 Blockers, 2026-07-16) ─────────────────
# Every artifact this sidecar certifies must POSTDATE every LISTED dependency -- the CODE that produced
# it and the ARTIFACTS it derived from. "Listed" is load-bearing: runtime reference data (rsid map,
# outcome GWAS, UKB-PPP tars, LD panels) is in neither map and is neither hashed nor staleness-checked. Without this, 33 hashes the CURRENT source alongside
# outputs an OLDER source produced and emits a record that is precisely false: it would state the
# full-grid PRIMARY rule (see the drawn-rule line below) while certifying bytes built under the
# saved by-method rule. A provenance file that certifies the wrong lineage is worse than none: it is believed.
#
# This gate took three rounds to get right, and each gap is worth stating so it is not reintroduced:
#   r8: existence was the only check -- a pre-fix artifact sailed through.
#   r9: checking artifact-vs-PRODUCER-SCRIPT missed a PARTIAL re-run: re-run 28 and 32 alone and their
#       outputs refresh while 29/30's still postdate their own untouched scripts, so every check passed
#       while this file hashed a NEW pool beside an OLD audit.
#   r10: checking 28-32 missed the SHARED code -- R/00_config.R, R/coloc.R, the 15/26 engines, the null
#       classifier. Those are in `scripts`/`EXTRA_LOGIC` below precisely BECAUSE they can change an
#       artifact; gating against only some of the code they hash was the same hole from a third side.
#
# So ONE map, not two: for each output, everything it depends on. The producer script is just its first
# code dependency -- keeping a separate PRODUCES map would let the two drift apart, which is the failure
# mode this whole exercise is about.
# The code lists are read off each script's ACTUAL source()/eval() structure, not assumed. Two things
# that are easy to get wrong and were: 29/31 source R/coloc.R ONLY on the non-eGFR branch
# (`if (MODE == "nonegfr")`) -- the eGFR path evals 26's head, and 26 sources 00_config alone -- so .COL
# hangs off the non-eGFR outputs only (r11: a first version attached it to both while claiming to be
# derived from what the scripts source; over-gating, and the claim was false).
# AND: this comment used to say "30 does NOT source 00_config". That was true when written and was made
# FALSE by r33, which gave 30 the truth guards it needed -- a fix that created a dependency and left the
# map, and the prose, describing the old world (r34 Blocker). If you change what a script sources, this
# map is part of the change.
#
# mtime is a WEAK signal: touching a file moves it, a fast enough edit-then-run could tie, and it says
# nothing about CONTENT. It is a tripwire against the stale case, not proof of lineage -- the md5s below
# are that, and only for what the "does NOT pin" section admits they cover.
# It fires today by design: 28/32 and R/00_config.R were edited 2026-07-16; the artifacts are 07-12.
.CFG <- file.path(B, "R/00_config.R")
.COL <- file.path(B, "R/coloc.R")
.CLS <- file.path(P, "_candidate_null_classify.R")
.E_E <- file.path(P, "26_susie_coloc_egfr.R")     # eGFR engine
.E_N <- file.path(P, "15_susie_coloc.R")          # non-eGFR engine
.SA  <- file.path(B, "R/strand_align.R")          # LD sign-alignment helper, sourced by BOTH engine heads (2026-07-17); run-defining for any output that runs on .E_E/.E_N
.dep <- list()
# output                                              <- code it runs on           + artifacts it derives from
.dep[[file.path(D, "candidate_pool.csv")]] <- c(
  file.path(P, "28_candidate_edge_enumerate.R"), .CFG,
  file.path(RC, c("mr_all_results.rds", "coloc_results_v7.rds", "mr_significant_bonferroni.rds")),
  PHASE1_INSTRUMENTS)
.dep[[file.path(D, "candidate_pool_no_coloc.csv")]] <- .dep[[file.path(D, "candidate_pool.csv")]]
.dep[[file.path(D, "candidate_susie_egfr.csv")]] <- c(          # no .COL: eGFR branch never sources coloc.R
  file.path(P, "29_candidate_susie.R"), .CLS, .CFG, .E_E, .SA, file.path(D, "candidate_pool.csv"))
.dep[[file.path(D, "candidate_susie_nonegfr.csv")]] <- c(
  file.path(P, "29_candidate_susie.R"), .CLS, .CFG, .COL, .E_N, .SA, file.path(D, "candidate_pool.csv"))
.dep[[file.path(D, "candidate_edge_audit.csv")]] <- c(
  # .CFG: since r33, 30 sources config and its classification depends on derive_runnability() and
  # assert_coloc_status_truth(). A config edit after the audit was produced would otherwise not trip this
  # gate, while 33 hashed the current config beside an audit built under the old one (r34 Blocker).
  file.path(P, "30_candidate_edge_audit.R"), .CLS, .CFG,
  file.path(D, c("candidate_pool.csv", "candidate_susie_egfr.csv", "candidate_susie_nonegfr.csv")))
# BOTH loose SuSiE files, for BOTH modes: since r30, each 31 run re-derives the dashed set from the
# complete loose output, so a mode's sensitivity output depends on the OTHER mode's file too (r31 High).
.dep[[file.path(D, "candidate_susie_sensitivity_egfr.csv")]] <- c(   # no .COL: eGFR branch never sources it
  file.path(P, "31_candidate_susie_sensitivity.R"), .CLS, .CFG, .E_E, .SA,
  file.path(D, c("candidate_edge_audit.csv", "candidate_susie_egfr.csv", "candidate_susie_nonegfr.csv",
                 "candidate_pool.csv")))
.dep[[file.path(D, "candidate_susie_sensitivity_nonegfr.csv")]] <- c(
  file.path(P, "31_candidate_susie_sensitivity.R"), .CLS, .CFG, .COL, .E_N, .SA,
  # candidate_pool.csv: 31 reads it directly, not only through the audit. The transitive audit edge would
  # usually catch a stale pool, but this map claims to list what each output depends on (r17 Low).
  # Both loose SuSiE files: see the eGFR entry above (r31 High).
  file.path(D, c("candidate_edge_audit.csv", "candidate_susie_egfr.csv", "candidate_susie_nonegfr.csv",
                 "candidate_pool.csv")))
.dep[[file.path(D, "dashed_edges_final.csv")]] <- c(
  # .CLS: since r43, 32 sources _candidate_null_classify.R for CAND_BENIGN_TAGS -- the single source for
  # which collapses are benign, and therefore for whether a `dropped:collapsed` edge is a biological
  # negative or a computation failure. A helper edit after this file was produced would not have tripped
  # the gate. Identical to r34's Blocker, one round after that lesson was written down (r44).
  file.path(P, "32_dashed_finalize.R"), .CFG, .CLS,
  file.path(D, c("candidate_edge_audit.csv", "candidate_susie_sensitivity_egfr.csv",
                 "candidate_susie_sensitivity_nonegfr.csv")),
  file.path(RC, "mr_all_results.rds"))          # 32 re-derives the primary FDR set from it
.viol <- character(0)
for (.o in names(.dep)) {
  .omt <- file.mtime(.o)
  for (.i in .dep[[.o]]) {
    .imt <- file.mtime(.i)
    if (is.na(.omt) || is.na(.imt))
      stop(sprintf("33: cannot read mtime for %s or %s -- refusing to certify what cannot be dated",
                   basename(.o), basename(.i)))
    if (.omt < .imt)
      .viol <- c(.viol, sprintf("%-42s %s  <  %s  %s", basename(.o), format(.omt, "%Y-%m-%d %H:%M"),
                                format(.imt, "%Y-%m-%d %H:%M"), basename(.i)))
  }
}
if (length(.viol))
  stop(sprintf("33: %d dependency violation(s) -- an artifact PREDATES a LISTED code or artifact dependency, so this sidecar would certify a stale or mixed-lineage set (a FALSE provenance record). Re-run the chain in order: 28 -> 29 (both modes) -> 30 -> 31 (both modes) -> 32.\n  %-42s %-17s    %-17s %s\n  %s",
               length(.viol), "artifact", "artifact mtime", "dependency mtime", "depends on",
               paste(.viol, collapse = "\n  ")))

kv <- function(paths, base_of = P) unname(vapply(paths, function(f) sprintf("- `%s` : %s",
            f, md5(if (identical(base_of, "")) f else file.path(base_of, f))), ""))

ln <- c(
  "# B1 sensitivity re-run + finalization — provenance sidecar",
  "",
  sprintf("- generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  # LABELLED, not fixed: these are read from THIS session, i.e. when the sidecar was generated -- not when
  # 29/31/32 produced the artifacts. The dependency gate has no package dimension, so a package upgrade
  # between the run and this call would be recorded here as though it were the producing environment
  # (r13 Medium). Saying so is the honest option; making it true needs 29/31 to persist session metadata
  # for 33 to read back, which is tracked in REBUILD_LEDGER §2.1c(d).
  sprintf("- R (present when THIS SIDECAR was generated, not necessarily when the artifacts were produced): %s",
          R.version.string),
  sprintf("- packages (same caveat -- present at sidecar generation): data.table %s | susieR %s | coloc %s | ieugwasr %s",
          ver("data.table"), ver("susieR"), ver("coloc"), ver("ieugwasr")),
  "",
  "## in-memory engine patch (A\")",
  sprintf("- patch: `min_abs_corr = 0 -> 0.5` (fixed-string, in memory; frozen 15/26 unchanged on disk)"),
  sprintf("- the same fixed string also rewrites `check_LD = FALSE -> TRUE`, which is INERT: `check_LD` is not a"),
  sprintf("  parameter of susieR::susie_rss (absent from its formals; swallowed by `...`), and susie_rss returns an"),
  sprintf("  identical object either way. It is patched only because the frozen engine's text contains it. The only"),
  sprintf("  setting susieR acts on here is min_abs_corr -- a post-hoc credible-set purity filter, applied after the"),
  sprintf("  IBSS loop, so it cannot affect convergence. This sensitivity tests PURITY, not external-LD consistency."),
  sprintf("- patch sites found: 15=%d, 26=%d (asserted == 2 each)", n15, n26),
  # "FDR q<0.05" alone is the exact ambiguity this repair exists to remove -- it reads as either
  # definition, and this sidecar ships with the submission package (r5 Low).
  "- SuSiE: L=10, min_abs_corr=0.5, max_iter=1000 (check_LD is passed but inert -- see the patch note above)",
  "  (max_iter=1000 is inherited from the frozen engines 15/26, which this sensitivity patches in memory; it",
  "   is a run-defining setting and is recorded here even though it is not part of the strict/loose diff.)",
  paste0("- drawn rule = full-grid PRIMARY FDR q<0.05 & strict PP.H4>=0.80 & n_warn==0. \"Primary FDR\" ",
         "= primary_fdr_set() (R/00_config.R): BH once over N_PROTEINS_CANONICAL x length(OUTCOMES) with ",
         "untested cells padded at p=1. NOT mr_significant_fdr.rds, which is the saved by-method artifact ",
         "that 28/32 read until 2026-07-15 -- a third quantity again, distinct from the sensitivity value ",
         "the manuscript prints; see the three-surface block above primary_fdr_set() in R/00_config.R."),
  "",
  "## scripts / engines (md5)",
  kv(scripts, P),
  # These live outside P so they need the absolute-path form. They are listed with the scripts, not
  # the inputs, because they ARE logic: primary_fdr_set() is the rule 28/32 apply, and R/coloc.R is the
  # engine 29/31 source.
  kv(EXTRA_LOGIC, ""),
  "",
  "## inputs (md5)",
  kv(inputs, ""),
  "",
  "## outputs (md5)",
  kv(outputs, ""),
  "",
  "## what this sidecar does NOT pin",
  "Runtime data reached indirectly through `R/00_config.R` and `R/coloc.R` is not hashed here: gene",
  "coordinates, the cis-rsid cache, the LD reference panels, the outcome GWAS files, and the UKB-PPP",
  "tars. They are GB-scale and pinning them is a separate decision, not something this sidecar settled.",
  "So, exactly: this RECORDS the code present when this sidecar was generated, HASHES the frozen artifacts",
  "it certifies, and CHECKS that each artifact postdates the LISTED code and LISTED artifact dependencies",
  "-- runtime reference data is neither hashed nor staleness-checked, see above. It is NOT",
  "proof of producing lineage -- an md5 of a file as it stands today plus an mtime comparison cannot",
  "establish that those bytes produced those artifacts. Nor does it certify that the reference data",
  "underneath was the same. Stated explicitly because a provenance file that is silent about its own",
  "boundary is read as claiming there isn't one."
)
writeLines(ln, OUT)
cat(sprintf("33: wrote %s\n\n", OUT)); cat(paste(ln, collapse = "\n"), "\n")
