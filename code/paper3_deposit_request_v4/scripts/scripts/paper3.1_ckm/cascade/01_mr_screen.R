# ============================================================================
# HASHED AS THE FIRST STATEMENT, before anything this script does can change the file, so the
# artifacts below can name the code that produced them and not merely the inputs they were made
# from. Same pattern and same purpose as cascade/03b's SELF_MD5.
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
SELF_MD5 <- local({
  .a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .p <- if (length(.a)) normalizePath(sub("^--file=", "", .a[1]), mustWork = FALSE) else ""
  if (nzchar(.p) && file.exists(.p)) unname(tools::md5sum(.p)) else NA_character_
})
# NOT OPTIONAL. A stamp that silently becomes NA records "some version of this script" and every
# comparison against it then passes, which is the opposite of what it is for.
CKM_RUN_ID <- sprintf("%s-%d-%s", format(Sys.time(), "%Y%m%dT%H%M%OS3"), Sys.getpid(),
                      substr(digest::digest(list(Sys.time(), Sys.getpid(), tempfile())), 1, 12))
if (is.na(SELF_MD5))
  stop("01_mr_screen: cannot determine its own path, so the artifacts it writes cannot record which version of this script produced them. Run it with Rscript.")
# PHASE 2: the corrected proteome-wide MR screen, on whichever Phase-1 generation is selected
# (R/_phase1_set.R). Self-contained: it builds its own LD-proxy cache, extracts all 8 local
# outcomes, runs offline MR and flags Bonferroni/FDR/suggestive. localize/02 and 03 did the first
# two of those for the original pipeline and are SUPERSEDED by the code here; their outputs have no
# reader.
# Every artifact it writes is named after the selected generation (config: LOCAL_EXTRACT_DIR,
# LD_PROXY_PREFIX, INSTR_SNPLIST, mr_artifact()), so a rebuild never overwrites the previous
# generation's files and no name can outlive what it names.
# NOTE the denominators differ and that is deliberate: BONFERRONI uses the full
# N_PROTEINS_CANONICAL x 8 grid; `fdr_sig` is BH WITHIN each method family over ANALYSABLE rows
# only, BEFORE the primary dedup. It is NOT the primary FDR, and it is NOT the
# sensitivity value the manuscript prints either -- that one is computed fresh by
# the generator as pooled BH after the dedup. Three distinct surfaces; they are
# named once, in the block above primary_fdr_set() in R/00_config.R. See also the
# caution at the fdr_sig line below. The drift block near the end is a COMPARISON, not a
# reproduction proof: it does not expect near-zero drift (see its own note).
# Writes repro_corrected/mr_all_results_<generation>.rds (+ _suggestive/_bonferroni/_fdr).
# ============================================================================
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"))
options(width = 200)
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
REPRO <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
# CREATED NOW, not discovered at the first write. Nothing wrote here until after PLINK, all
# eight outcome reads and every MR fit, so on a clean tree the run spent hours and then failed
# on a missing directory -- having already published the proxy cache and the outcome extracts.
dir.create(REPRO, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(REPRO)) stop("01_mr_screen: cannot create the output directory ", REPRO)
t_start <- proc.time()["elapsed"]

# ---- (0) obtain the SELECTED instrument set ----
# 2026-07-01 (gate finding 2): Phase 1 is REUSED/frozen. If the canonical, assay-deduped file
# EXISTS, read it DIRECTLY (do NOT rebuild from the mutable 1635 + delta53 sources -> no silent
# value/OID/allele drift). The rebuild is only the first-time fallback.
# NAMED FOR WHAT IT IS, not for a gene count. `.p1_file` described one generation and kept the
# name after Phase 1 was rebuilt to 1,948 genes, so every message printed here named a file the
# variable no longer held.
.p1_file <- PHASE1_INSTRUMENTS
if (file.exists(.p1_file)) {
  # read_phase1_instruments(), not a second readRDS of the path: config read this file once,
  # checked those bytes against the builder's certificate and derived N_PROTEINS_CANONICAL and
  # BONFERRONI_ALPHA from them. Re-opening the path here would let this script analyse one
  # generation while the threshold it applies describes another.
  inst <- as.data.table(read_phase1_instruments())
  # WHICH instrument set this is, said out loud. The frozen branch is taken whenever the file at
  # PHASE1_INSTRUMENTS exists, and that path defaults to results/repro_phase1_1682/ -- so a rebuilt
  # Phase 1 written anywhere else is ignored and the cascade continues on the previous instruments
  # while printing "reusing frozen foundation", which reads like a deliberate choice. Naming the path
  # and the counts makes a silent substitution visible in the log of every run.
  cat(sprintf("01_mr_screen: Phase-1 foundation = %s (%d genes / %d rows / %d distinct SNP)\n",
              .p1_file, uniqueN(inst$gene), nrow(inst), uniqueN(inst$SNP)))
  # The gene count is checked against a constant fixed for the PREVIOUS instrument set. A rebuild
  # legitimately changes it, and then the correct wiring aborts here while stale wiring passes -- so
  # the mismatch must name both numbers and the file, rather than failing as a bare stopifnot with no
  # indication of which of the two is out of date. Whether the constant should move is an analysis
  # decision and is deliberately NOT made here.
  if (uniqueN(inst$gene) != N_PROTEINS_CANONICAL)
    stop(sprintf(paste0("01_mr_screen: %s holds %d genes but N_PROTEINS_CANONICAL is %d. If this is the ",
                        "rebuilt instrument set, the canonical count and the Bonferroni denominator ",
                        "derived from it (R/00_config.R) have not been updated to match it; if it is ",
                        "not, PHASE1_INSTRUMENTS is pointing at the wrong table."),
                 .p1_file, uniqueN(inst$gene), N_PROTEINS_CANONICAL))
  stopifnot(inst[, uniqueN(oid), by = gene][V1 > 1, .N] == 0)   # single canonical assay per gene
  cat("01_mr_screen: reusing the frozen Phase-1 foundation named above (single-OID/gene)\n")
} else {
  # UNREACHABLE AS THE PIPELINE NOW STANDS, and said so rather than left looking like a live option:
  # R/00_config.R reads and certifies the selected instrument RDS while it is being sourced, above,
  # so a missing file stops the run there and execution never arrives here. The branch is kept
  # because it is the record of how the 1,682-gene generation was assembled from canon + delta; it
  # cannot run, and it cannot affect any number.
  # THE OTHER BRANCH BUILDS THE FILE THE FIRST BRANCH THEN REUSES FOREVER. Once installed, changed
  # canon/delta inputs are ignored on every later run, and deleting the file makes this same code
  # produce a different foundation -- so the result depends on filesystem history rather than on the
  # inputs. Said in the log because the branch is otherwise invisible in the output.
  # ...and that installation happens under whatever name PHASE1_INSTRUMENTS carries. When that name
  # denotes a REBUILT generation that has not been produced yet, this branch reconstructs the
  # PREVIOUS one from canon + delta and installs it there, after which every later run reuses it as
  # though it were the rebuild. Nothing downstream can tell: the reconstruction is internally
  # consistent and its gene count is exactly the constant it is checked against. So the branch has
  # to be asked for rather than reached by the absence of a file.
  # THIS BRANCH NO LONGER RECONSTRUCTS ANYTHING, and the reasons are cumulative.
  # It cannot be reached in the current pipeline: R/00_config.R reads and certifies the selected
  # instrument RDS while it is being sourced above, so an absent file stops the run there. If it
  # WERE reached -- the file deleted between config and here -- its recovery does not exist: it
  # would not re-read the missing table but rebuild the PREVIOUS 1,682-gene generation from canon
  # plus delta and install it under whatever name PHASE1_INSTRUMENTS carries. During a rebuild
  # that name denotes the NEW generation, so every later run would reuse the old instruments as
  # though they were the rebuild, and nothing downstream could tell: the reconstruction is
  # internally consistent and its gene count is exactly the constant it is checked against.
  # That is the worst failure available to this script, and it was guarded by an environment
  # variable test that accepted "0" and "false" as permission, followed by a file.copy() whose
  # result was not checked. Fixing three defects in code that cannot run and must not run is not
  # the same as removing the hazard.
  # The reconstruction it performed is recorded in REBUILD_LEDGER.md; the 1,682-gene generation
  # it produced is results/repro_phase1_1682/, which still exists and is still readable.
  stop(sprintf(paste0("01_mr_screen: %s does not exist. This script does NOT reconstruct a ",
                      "Phase-1 foundation: reconstructing from canon plus delta yields the ",
                      "PREVIOUS 1,682-gene generation, and installing it under this path would ",
                      "make every later run analyse the old instruments while reporting the ",
                      "rebuild. Point PHASE1_INSTRUMENTS at the instrument table you mean to ",
                      "analyse, or rebuild it with the Phase-1 scripts."),
               .p1_file))
}
cat(sprintf("01_mr_screen: %s instrument set = %d rows / %d distinct SNP / %d gene\n",
            PHASE1_TAG, nrow(inst), uniqueN(inst$SNP), uniqueN(inst$gene)))
# 2026-07-01 fix: gene-set helpers usable in BOTH branches (canon/delta objects exist only when building).
# THE canon/delta SPLIT IS A PROPERTY OF ONE GENERATION, not of the analysis. It records how the
# 1,682-gene set was assembled (1,635 canonical + a 47-gene local recovery). The GRCh38 rebuild is a
# single extraction and has no delta, so on any other generation these two vectors would slice the
# instrument set by a partition that does not apply to it: `.canon_genes` would silently mean "every
# gene except 47 names from another generation", and the two report blocks at the end would print
# subsets that look meaningful and are not.
.has_delta   <- ckm_is_known_previous(.p1_file) &&
                file.exists("results/repro_phase1/instruments_delta53.rds")
.delta_genes <- if (.has_delta) {
                  unique(as.data.table(readRDS("results/repro_phase1/instruments_delta53.rds"))$gene)
                } else character(0)
.canon_genes <- if (.has_delta) setdiff(unique(inst$gene), .delta_genes) else unique(inst$gene)
if (!.has_delta)
  cat("01_mr_screen: this generation has no canon/delta split -- the drift comparison runs over all genes and
",
      "    the delta-hits listing is skipped.
", sep = "")

# ---- (A) combined LD-proxy table for the selected instrument set (PLINK) ----
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_ld_namespace.R"))
# Panel-space, for the same reason as localize/03: this builds a SECOND proxy cache with the same
# PLINK parameters, and the two must agree about which variants can be proxied.
# NA means the reference panel does not contain this variant, so it is dropped rather than replaced
# by its outcome-space name: asking PLINK for a name the panel was already known not to hold either
# finds nothing, or finds a record the position/allele matching had rejected. On a legacy table the
# panel key IS SNP and nothing is dropped. ld_query_keys() applies both rules; the `.ldkey` column
# this line used to add to `inst` is gone with the inline filter that was its only reader.
# THE SHARED DEFINITION, not a second copy of it. This line used to filter to rs-format inline,
# which silently DROPS an instrument whose rsid_ld_panel is present but malformed -- a value that
# claims the panel holds the variant while naming something PLINK cannot be asked for. The same
# rows are refused by ld_query_keys(), but only when ld_query_index_md5() reaches it while the
# provenance record is built -- which is AFTER PLINK has run without them and after the proxy
# cache has been written. Calling the helper here makes the refusal happen before either.
snps <- ld_query_keys(inst)
# THE NAMES FOLLOW THE GENERATION (config: INSTR_SNPLIST / LD_PROXY_PREFIX). They were literal
# "instr_1682.*", which after the rebuild would have written a 1,948-gene cache into the file the
# previous generation's consumers read -- and the previous generation's cache is what the submitted
# numbers were computed from.
snpfile <- INSTR_SNPLIST; writeLines(snps, snpfile)
out <- LD_PROXY_PREFIX
if (file.exists(paste0(out, ".rds"))) { file.remove(paste0(out, ".rds"))  # 2026-07-01 gate: always rebuild from the current snplist
  if (file.exists(paste0(out, ".rds"))) stop(sprintf("01_mr_screen: could NOT remove stale proxy cache %s -- a failed remove must not silently reuse old proxies", paste0(basename(out), ".rds"))) }  # Codex-phase2 #1
if (!file.exists(paste0(out, ".rds"))) {
  cat(sprintf("PLINK --r2 in-phase for %d SNPs ...\n", length(snps))); flush(stdout())
  if (file.exists(paste0(out, ".ld"))) { file.remove(paste0(out, ".ld"))   # 2026-06-30 gate obj #3: no stale .ld fallback
    if (file.exists(paste0(out, ".ld"))) stop("01_mr_screen: could NOT remove stale .ld -- refusing stale fallback") }  # Codex-phase2 #1
  # THE PANEL IS AUTHENTICATED BEFORE IT IS USED, not recorded after. Its BED/BIM/FAM hashes were
  # written into the provenance sidecar below -- after PLINK had already run on them -- and never
  # compared with anything, so replacing EUR.bed or EUR.fam while leaving the BIM identifiers intact
  # changes r2, phase, every proxy choice and every MR estimate downstream, and the build reports
  # success. The BIM hash certified by Phase 1 is already carried on the instrument table
  # (attribute ckm_ld_panel_md5, written by 09_phase1_harmonise_hg38.R), so this compares what is
  # about to be read against what Phase 1 recorded reading.
  # ALL THREE PANEL FILES, not just the BIM. The BIM names the variants; the BED holds the genotypes
  # and the FAM the samples, so replacing either changes r2, phase, which proxy is chosen and every
  # MR estimate that follows, while the variant names -- and therefore a BIM-only check -- stay
  # identical. They were merely RECORDED into the provenance sidecar AFTER PLINK had read them,
  # which documents a swap instead of preventing it.
  # EUR_MANIFEST.tsv is the panel's own certificate (role=panel rows), so this compares against a
  # recorded truth rather than a hash this script generated itself.
  .pman <- tryCatch(fread(file.path(LD_EUR_DIR, "EUR_MANIFEST.tsv")), error = function(e) NULL)
  if (is.null(.pman) || !all(c("role", "file", "md5") %in% names(.pman)))
    stop("01_mr_screen: no usable EUR_MANIFEST.tsv beside the LD panel, so the genotypes PLINK is ",
         "about to read cannot be checked against anything.")
  { .dup <- .pman[role == "panel"][, .N, by = file][N > 1L]
    if (nrow(.dup))
      stop("01_mr_screen: EUR_MANIFEST.tsv certifies ", paste(.dup$file, collapse = ", "),
           " more than once. Taking the first row makes two disagreeing certificates look like one.") }
  for (.pf in c("EUR.bed", "EUR.bim", "EUR.fam")) {
    .exp <- .pman[role == "panel" & file == .pf, md5][1]
    if (is.na(.exp) || !nzchar(.exp))
      stop("01_mr_screen: EUR_MANIFEST.tsv does not certify ", .pf, ".")
    .got <- unname(tools::md5sum(file.path(LD_EUR_DIR, .pf)))
    if (!identical(tolower(.got), tolower(.exp)))
      stop(sprintf(paste0("01_mr_screen: %s hashes to %s but its manifest certifies %s. LD, phase ",
                          "and proxy choice -- and every MR estimate downstream -- would come from ",
                          "a panel other than the certified one."),
                   .pf, substr(.got, 1, 12), substr(.exp, 1, 12)))
  }
  # ...and the BIM must ALSO be the one Phase 1 defined its instruments against.
  # THE WHOLE PANEL, not only the BIM. Replacing EUR.bed or EUR.fam and updating its manifest row
  # passes the loop above, and the BIM check below cannot see it -- yet PLINK derives LD, phase and
  # proxy choice from the genotypes and samples. Phase 1 recorded which manifest it pinned; that is
  # the identity to compare.
  .pman_expect <- attr(read_phase1_instruments(), "ckm_ld_panel_manifest_md5")
  if (!is.null(.pman_expect) && !is.na(.pman_expect) &&
      !identical(tolower(unname(tools::md5sum(file.path(LD_EUR_DIR, "EUR_MANIFEST.tsv")))),
                 tolower(.pman_expect)))
    stop("01_mr_screen: EUR_MANIFEST.tsv is not the panel manifest Phase 1 pinned (recorded ",
         substr(.pman_expect, 1, 12), "). The proxies would come from a panel other than the one ",
         "the instruments were defined against.")
  .bim_expect <- attr(read_phase1_instruments(), "ckm_ld_panel_md5")
  .bim_now    <- unname(tools::md5sum(file.path(LD_EUR_DIR, "EUR.bim")))
  if (!is.null(.bim_expect) && !is.na(.bim_expect) &&
      !identical(tolower(.bim_now), tolower(.bim_expect)))
    stop(sprintf(paste0("01_mr_screen: the LD panel BIM at %s hashes to %s, but the Phase-1 instrument table ",
                        "was built against %s. Proxies, r2 and phase would come from a different ",
                        "panel than the instruments were defined on."),
                 file.path(LD_EUR_DIR, "EUR.bim"), substr(.bim_now, 1, 12), substr(.bim_expect, 1, 12)))
  # Kept, so the same bytes can be re-checked after PLINK has read them: verifying a file and
  # then handing its PATH to another program leaves the interval between the two unguarded, and
  # the sidecar below was recording the after-hashes without ever comparing them.
  .panel_before <- vapply(c("EUR.bed", "EUR.bim", "EUR.fam"),
                          function(.f) unname(tools::md5sum(file.path(LD_EUR_DIR, .f))), character(1))
  .plink_args <- c("--bfile", shQuote(file.path(LD_EUR_DIR,"EUR")), "--r2","in-phase","with-freqs",
    "--ld-snp-list", shQuote(snpfile), "--ld-window-kb","1000","--ld-window","99999","--ld-window-r2","0.8",
    "--out", shQuote(out))
  .pst <- system2(PLINK_BIN, .plink_args, stdout=TRUE, stderr=TRUE)
  if (!is.null(attr(.pst,"status")) && attr(.pst,"status")!=0)
  # PLINK's own output was captured and thrown away, so "check PLINK/bfile" was all anyone got.
  stop("01_mr_screen: PLINK --r2 failed (exit ", attr(.pst, "status"), ") — LD proxy NOT built. ",
       "PLINK said:
  ", paste(utils::tail(.pst, 12), collapse = "
  "))
  if (!file.exists(paste0(out, ".ld"))) stop("01_mr_screen: PLINK produced no .ld — refusing to build proxy from stale/absent data")
  ld <- fread(paste0(out, ".ld"))
# EVERY requested key must come back. PLINK emits a self-pair for each identifier it recognises, so
# a key it never saw -- a stale or mistyped rsid_ld_panel, or one from a previous BIM -- simply does
# not appear, the remaining instruments still produce a nonempty .ld, and that instrument later has
# no proxy row and disappears without trace. Only for tables that CARRY rsid_ld_panel: the legacy
# snplist has 202 identifiers absent from the current BIM and that behaviour must not change.
.unseen <- setdiff(snps, unique(ld$SNP_A))
if (length(.unseen)) {
  if ("rsid_ld_panel" %in% names(inst))
    stop(sprintf("%s: PLINK did not recognise %d of %d requested panel identifier(s) (e.g. %s). ",
                 "cascade/01_mr_screen", length(.unseen), length(snps),
                 paste(utils::head(.unseen, 4), collapse = ", ")),
         "Those instruments would have no proxy and would be dropped silently.")
  # A legacy table has no panel name to give, so its outcome-space identifier goes to PLINK and 202
  # of them are not in this BIM: no proxy row, and the instrument vanishes from every outcome that
  # lacks its direct name. Aborting would only make the legacy path unrunnable -- the real repair is
  # the rebuilt table's rsid_ld_panel. What must not continue is the SILENCE: FST is queried as
  # rs4370233 while the panel holds the same A/C variant as rs184943135, 34 proxies at r2>=0.8 exist,
  # 26 are orientable in Stroke, CKD and eGFR, and FST is absent from all three stored extracts.
  .uf <- file.path(LD_EUR_DIR, paste0(basename(LD_PROXY_PREFIX), "_UNSEEN_KEYS.txt"))
  writeLines(sort(.unseen), .uf)
  cat(sprintf(paste0("01_mr_screen: WARNING -- the reference panel does not hold %d of %d requested identifier(s), ",
                     "covering %d gene(s). Those instruments get NO proxy and will be missing from every ",
                     "outcome that lacks their direct name. This is the legacy table's namespace gap, ",
                     "resolved by rebuilding Phase 1 so the table carries rsid_ld_panel. Full list: %s\n"),
              length(.unseen), length(snps), uniqueN(inst$gene[inst$SNP %in% .unseen]), .uf))
}
  ok <- grepl("^[ACGT][ACGT]/[ACGT][ACGT]$", ld$PHASE)
  prox <- ld[ok & SNP_A != SNP_B]
  prox[, `:=`(a_hap1=substr(PHASE,1,1), b_hap1=substr(PHASE,2,2), a_hap2=substr(PHASE,4,4), b_hap2=substr(PHASE,5,5))]
  setnames(prox, c("SNP_A","SNP_B","R2"), c("instrument","proxy","r2"), skip_absent=TRUE)
  prox <- prox[, .(instrument, proxy, r2, MAF_A, MAF_B, a_hap1, b_hap1, a_hap2, b_hap2)]
  # PER INDEX, not globally. The guard below asks only whether ANY proxy for ANY instrument
  # survived, so an index variant whose every non-self edge was dropped as malformed PHASE passes
  # unnoticed: it gets no proxy row, disappears from every outcome that lacks its direct name, and is
  # then classified as an ordinary no_outcome_data cell. The raw .ld already says which indices had
  # edges to lose, so comparing the two costs one setdiff and no new inputs.
  .raw_idx  <- unique(ld[SNP_A != SNP_B, SNP_A])
  .lost_idx <- setdiff(.raw_idx, unique(prox$instrument))
  if (length(.lost_idx))
    stop(sprintf(paste0("01_mr_screen: %d index variant(s) had LD partners in the raw PLINK output ",
                        "but none survived phase parsing (e.g. %s). Their instruments would get no ",
                        "proxy and would be recorded as ordinary missing outcome data. Inspect %s."),
                 length(.lost_idx), paste(utils::head(.lost_idx, 5), collapse = ", "),
                 paste0(out, ".ld")))
  if (nrow(prox) == 0L || uniqueN(prox$instrument) == 0L)   # P2-F01: a present-but-degenerate .ld (blank / non-ACGT PHASE, all self-pairs) must NOT silently save an empty proxy table -> keep_set collapses to direct-only and MR can regress without a proxy-specific abort
    stop(sprintf("01_mr_screen: LD proxy parse produced %d rows / %d instruments from %d raw .ld pairs -- refusing to save an empty proxy table (would silently drop every proxy substitution)", nrow(prox), uniqueN(prox$instrument), nrow(ld)))
  .panel_after <- vapply(c("EUR.bed", "EUR.bim", "EUR.fam"),
                         function(.f) unname(tools::md5sum(file.path(LD_EUR_DIR, .f))), character(1))
  if (!identical(tolower(.panel_before), tolower(.panel_after)))
    fail_loud_data(sprintf(paste0("01_mr_screen: the LD panel changed while PLINK was reading it ",
                                  "(%s). The proxies just computed do not come from the panel that ",
                                  "was certified."),
                           paste(names(.panel_before)[tolower(.panel_before) != tolower(.panel_after)],
                                 collapse = ", ")))
  setorder(prox, instrument, -r2); saveRDS(prox, paste0(out, ".rds"))
  # The same sidecar localize/03 writes, because this is a SECOND cache built with the same
  # parameters and localize/41_3 consumes it. Without one, 41_3 cannot show the cache was queried in the panel
  # namespace against this instrument set, and its only remaining choice is to refuse to run.
  saveRDS(list(built = format(Sys.time(), tz = "UTC"),
               n_index = length(snps),
               index_md5 = ld_query_index_md5(inst),
               cache_md5 = unname(tools::md5sum(paste0(out, ".rds"))),
               bfile_bim_md5 = unname(tools::md5sum(file.path(LD_EUR_DIR, "EUR.bim"))),
               bfile_bed_md5 = unname(tools::md5sum(file.path(LD_EUR_DIR, "EUR.bed"))),
               bfile_fam_md5 = unname(tools::md5sum(file.path(LD_EUR_DIR, "EUR.fam"))),
               plink_bin = PLINK_BIN,
               plink_md5 = if (file.exists(PLINK_BIN)) unname(tools::md5sum(PLINK_BIN)) else NA_character_,
               plink_args = .plink_args,
               ld_params = "r2>=0.8, +/-1Mb, in-phase with-freqs",
               namespace = if ("rsid_ld_panel" %in% names(inst)) "panel" else "legacy_snp"),
          paste0(out, ".provenance.rds"))
}
prox <- readRDS(paste0(out, ".rds")); setDT(prox)
# The cache is written by PLINK in PANEL space; every join below is keyed on SNP. BOTH columns are
# translated, not just the index: a proxy is normally an arbitrary panel variant and passes through
# unchanged, but 27 proxy variants in the current cache are themselves instruments (34 proxy edges),
# and for those the mapping is known and must be applied or the lookup searches the outcome for a
# name it does not carry. On a legacy table the map is the identity and nothing changes.
.ldmap <- ld_to_snp_map(inst)
prox[, instrument := ld_names_to_snp(instrument, .ldmap)]
prox[, proxy      := ld_names_to_snp(proxy,      .ldmap)]
# A row whose instrument or proxy translated to NA is UNUSABLE, so it goes. Left in, the NA reaches
# keep_set, and `NA %in% c(NA)` is TRUE in R -- so outcome rows with a missing rsID (1,489,577 of
# them in FinnGen T2D) matched it, and one of them could be attached to an instrument as its proxy.
prox <- prox[!is.na(instrument) & !is.na(proxy)]
cat(sprintf("combined proxy table: %d pairs, %d instruments with >=1 proxy\n", nrow(prox), uniqueN(prox$instrument)))

# ---- (B) outcome extraction (mirror localize/02 exactly) ----
# SNP is the OUTCOME-space name; rs format is a property of the LD panel, not of it (.ldkey above is
# where that requirement belongs, and it is applied there). Requiring it here dropped instruments
# before any outcome was read: no direct match, and no proxy either.
# eaf_exp is the LEAST informative exposure frequency among the proteins sharing this SNP. Collapsing
# to the first protein's value let one protein's informative frequency keep a direct palindromic row
# that harmonise_data(action=2) then drops for another -- after the proxy path had already closed.
inst_al <- inst[!is.na(SNP) & nzchar(SNP),
                .(ie=toupper(effect_allele)[1], io=toupper(other_allele)[1],
                  eaf_exp = eaf[which.min(abs(eaf - 0.5))][1]), by=SNP]
setkey(inst_al, SNP)
instr_rsids <- inst_al$SNP
keep_set <- unique(c(instr_rsids, unique(prox$proxy)))
# 2026-07-01: N/ncase are READ from config's SINGLE SOURCE (OUTCOME_N_TOTAL / OUTCOME_N_CASE, guarded by
# config's stopifnot) — NOT hardcoded here, so this can never drift from config/coloc/09_steiger. id ->
# trait-name via config OUTCOMES. eGFR stays quant (per-SNP n_total_sum); CKD keeps n_col but its nominal N
# comes from config. (config values == the previous literals here, so results are unchanged.)
.id2nm <- setNames(names(OUTCOMES), unname(OUTCOMES))              # "FinnGen_R10_T2D" -> "T2D", ...
.oN  <- function(oid) { nm <- .id2nm[[oid]]; v <- if (!is.null(nm)) OUTCOME_N_TOTAL[[nm]] else NA
                        if (is.null(v) || is.na(v)) stop(sprintf("01_mr: no config OUTCOME_N_TOTAL for outcome %s", oid)); as.integer(v) }
.oNC <- function(oid) { nm <- .id2nm[[oid]]; v <- if (!is.null(nm)) OUTCOME_N_CASE[[nm]] else NA
                        if (is.null(v) || is.na(v)) stop(sprintf("01_mr: OUTCOME_N_CASE is missing or NA for case-control outcome %s -- config missing ncase (only cc outcomes call .oNC)", oid)); as.integer(v) }  # 2026-07-01 Codex-phase2 #4
fg <- function(file, oid) list(file=file, sep="\t", rsid="rsids", ea="alt", oa="ref", beta="beta", se="sebeta", eaf="af_alt", pval="pval", n_fixed=.oN(oid), ncase=.oNC(oid), type="cc", multi_rsid=TRUE)
LOCAL_META <- list(
  "FinnGen_R10_I9_CHD"=fg("finngen_R10_I9_CHD.gz","FinnGen_R10_I9_CHD"), "FinnGen_R10_I9_HEARTFAIL"=fg("finngen_R10_I9_HEARTFAIL.gz","FinnGen_R10_I9_HEARTFAIL"),
  "FinnGen_R10_T2D"=fg("finngen_R10_T2D.gz","FinnGen_R10_T2D"), "FinnGen_R10_NAFLD"=fg("finngen_R10_NAFLD.gz","FinnGen_R10_NAFLD"),
  "FinnGen_R10_I9_AF"=fg("finngen_R10_I9_AF.gz","FinnGen_R10_I9_AF"),
  "MEGASTROKE_EUR_GCST006906"=list(file="MEGASTROKE_EUR_harmonised_hg38.h.tsv.gz",sep="\t",rsid="hm_rsid",rsid_fallback="variant_id",ea="hm_effect_allele",oa="hm_other_allele",beta="hm_beta",se="standard_error",eaf="hm_effect_allele_frequency",pval="p_value",
    # All-or-nothing per row: an unharmonisable row has every hm_ field NA while the submitted
    # columns are complete, so substituting only the identifier left a row with no alleles and no
    # beta -- dropped, and the instrument then looked absent (localize/02 carries this).
    ea_fallback="effect_allele",oa_fallback="other_allele",beta_fallback="beta",eaf_fallback="effect_allele_frequency",
    n_fixed=.oN("MEGASTROKE_EUR_GCST006906"),ncase=.oNC("MEGASTROKE_EUR_GCST006906"),type="cc"),
  "Wuttke2019_CKD_EA"=list(file="CKDGen_CKD_Wuttke2019_EA.txt.gz",sep=" ",rsid="RSID",ea="Allele1",oa="Allele2",beta="Effect",se="StdErr",eaf="Freq1",pval="P-value",n_col="n_total_sum",n_fixed=.oN("Wuttke2019_CKD_EA"),ncase=.oNC("Wuttke2019_CKD_EA"),type="cc"),
  "Wuttke2019_eGFR_EA"=list(file="CKDGen_eGFR_EA_hg37.txt.gz",sep=" ",rsid="RSID",ea="Allele1",oa="Allele2",beta="Effect",se="StdErr",eaf="Freq1",pval="P-value",n_col="n_total_sum",type="quant"))

# The run's single record of what each outcome file was. Computed on first use and reused, so every
# artifact this run writes names the same bytes even if the file changes underneath.
.OUTCOME_MD5 <- new.env(parent = emptyenv())
# The hash taken BEFORE the extraction, for the sources that publish no checksum of their own. It is
# not provenance -- it cannot say the bytes are the ones the provider released -- but it is what the
# post-read comparison needs so that a file replaced mid-run is noticed rather than stamped with the
# replacement's hash.
.OUTCOME_PRE <- new.env(parent = emptyenv())
# The SOURCE's own checksums, where the provider publishes them. This is the difference between
# DETECTING CHANGE and ESTABLISHING IDENTITY: a hash this script computes says only "the same bytes
# as the last run", while a hash the data provider published says "these are the bytes they
# released". Fetched once by HEAD request (no re-download of 800 MB) into
# OUTCOME_SOURCE_CHECKSUMS.tsv and used offline thereafter, per the local-first rule.
# An outcome ABSENT from the manifest means no such certificate exists for it -- reported as that,
# never silently treated as a pass.
.SRC_CHK <- local({
  .f <- file.path(OUTCOME_DIR, "OUTCOME_SOURCE_CHECKSUMS.tsv")
  if (file.exists(.f)) fread(.f, colClasses = "character") else NULL
})
# THE SOURCES THAT PUBLISH NO CHECKSUM, named. Without this list a MISSING certificate and a source
# that never had one are the same event, so deleting a row -- or the whole manifest -- turns every
# check off without a word. Membership here is a fact about the provider, not a threshold.
.NO_PUBLISHED_CHECKSUM <- c("CKDGen")   # 2026-08-05: MEGASTROKE removed -- EBI DOES publish an md5 for
# the harmonised file, in <file>-meta.yaml. Claiming otherwise kept the strongest available check off
# for the one outcome that turned out to be the wrong accession.
.verify_outcome_source <- function(fname, oid) {
  .exempt <- any(vapply(.NO_PUBLISHED_CHECKSUM, function(.k) grepl(.k, fname, ignore.case = TRUE),
                        logical(1)))
  if (!is.null(.SRC_CHK) && !all(c("file", "outcome_id") %in% names(.SRC_CHK)))
    fail_loud_data(paste0("01_mr_screen: OUTCOME_SOURCE_CHECKSUMS.tsv has no outcome_id column, so a ",
                          "checksum can certify a FILE without binding it to a phenotype -- which is ",
                          "the substitution this check exists to prevent."))
  if (is.null(.SRC_CHK) || !("file" %in% names(.SRC_CHK))) {
    if (.exempt) return("provider publishes no checksum")
    fail_loud_data(sprintf(paste0("01_mr_screen: OUTCOME_SOURCE_CHECKSUMS.tsv is absent, so %s ",
                                  "cannot be checked against what its provider released."), fname))
  }
  .r <- .SRC_CHK[file == fname]
  # A SIZE IS NOT A CHECKSUM, and is recorded as neither more nor less than it is. The three sources
  # that publish no checksum do publish a size, and each local file was matched against its release
  # (Stroke additionally byte-for-byte on its first and last megabyte). Weaker than a checksum, but
  # it is an established identity rather than none, and the artifacts must not call it "verified".
  if (.exempt) {
    # The SAME structural rules as a checksummed source: one row, and that row must name this
    # phenotype. Returning early skipped both, so deleting or duplicating a row silently disabled the
    # only automated release check these three have.
    if (!nrow(.r))
      fail_loud_data(sprintf(paste0("01_mr_screen: %s has no row in OUTCOME_SOURCE_CHECKSUMS.tsv. ",
                                    "Its provider publishes no checksum, so that row (release URL ",
                                    "and size) is its only recorded identity."), fname))
    if (nrow(.r) > 1L)
      fail_loud_data(sprintf("01_mr_screen: %s has %d rows in OUTCOME_SOURCE_CHECKSUMS.tsv; one file, one certificate.",
                             fname, nrow(.r)))
    if ("outcome_id" %in% names(.SRC_CHK) && !identical(.r$outcome_id[1], oid))
      fail_loud_data(sprintf(paste0("01_mr_screen: %s is the recorded file for outcome %s, but it is ",
                                    "being read as %s."), fname, .r$outcome_id[1], oid))
    if (!("published_size" %in% names(.SRC_CHK)) || !nzchar(.r$published_size[1]))
      fail_loud_data(sprintf(paste0("01_mr_screen: %s has a row but no published_size, so nothing ",
                                    "about it can be checked against its release."), fname))
    {
      .sz <- file.info(file.path(OUTCOME_DIR, fname))$size
      if (!identical(as.numeric(.sz), as.numeric(.r$published_size[1])))
        fail_loud_data(sprintf(paste0("01_mr_screen: %s is %.0f bytes but its release is %s bytes ",
                                      "(%s). Its provider publishes no checksum, so size is the only ",
                                      "identity available and it does not match."),
                               fname, .sz, .r$published_size[1], .r$source[1]))
      return("size matches the release; provider publishes no checksum")
    }
    return("provider publishes no checksum")
  }
  if (!nrow(.r)) {
    if (.exempt) return("provider publishes no checksum")
    fail_loud_data(sprintf(paste0("01_mr_screen: %s has no row in OUTCOME_SOURCE_CHECKSUMS.tsv. Its ",
                                  "provider does publish checksums, so a missing row means the ",
                                  "certificate was lost, not that none exists."), fname))
  }
  # TWO CERTIFICATES FOR ONE FILE IS NOT A CERTIFICATE. Taking [1] silently picks whichever was
  # written first, so a second row disagreeing on outcome_id or published_md5 decides nothing.
  if (nrow(.r) > 1L)
    fail_loud_data(sprintf("01_mr_screen: %s has %d rows in OUTCOME_SOURCE_CHECKSUMS.tsv; one file, one certificate.",
                           fname, nrow(.r)))
  # THE FILE *AND* THE PHENOTYPE IT IS. Verifying the checksum alone proves only that these bytes
  # are the ones the provider released -- not that they are the released file for the outcome this
  # loop is about to label them with. Point the T2D metadata at the certified CHD filename and the
  # checksum passes, after which CHD associations are published as T2D.
  if ("outcome_id" %in% names(.SRC_CHK) && !identical(.r$outcome_id[1], oid))
    fail_loud_data(sprintf(paste0("01_mr_screen: %s is the certified file for outcome %s, but it ",
                                  "is being read as %s."), fname, .r$outcome_id[1], oid))
  .got <- unname(tools::md5sum(file.path(OUTCOME_DIR, fname)))
  if (!identical(tolower(.got), tolower(.r$published_md5[1])))
    fail_loud_data(sprintf(paste0("01_mr_screen: %s does not match the checksum its provider ",
                                  "published (%s on disk, %s from %s). This is not the released file."),
                           fname, substr(.got, 1, 12), substr(.r$published_md5[1], 1, 12), .r$source[1]))
  "verified"
}


.outcome_md5 <- function(fname) {
  if (!exists(fname, envir = .OUTCOME_MD5, inherits = FALSE)) {
    .h <- unname(tools::md5sum(file.path(OUTCOME_DIR, fname)))
    # THE BYTES THAT WERE READ, NOT THE BYTES THAT WERE CHECKED. .verify_outcome_source() hashes the
    # file BEFORE the extraction; this hash is taken after it, and is the one stamped on the artifact.
    # Comparing it with the published checksum here is what makes the two hashes describe one file:
    # a replacement between the check and the read fails at this line instead of producing an extract
    # from uncertified bytes stamped with the replacement's own hash.
    if (!is.null(.SRC_CHK) && "file" %in% names(.SRC_CHK)) {
      .c <- .SRC_CHK[file == fname]
      if ((nrow(.c) == 0L || !nzchar(.c$published_md5[1])) &&
          exists(fname, envir = .OUTCOME_PRE, inherits = FALSE) &&
          !identical(tolower(.h), tolower(get(fname, envir = .OUTCOME_PRE, inherits = FALSE))))
        fail_loud_data(sprintf(paste0("01_mr_screen: %s changed while it was being read (%s before, ",
                                      "%s after). Its provider publishes no checksum, so this is the ",
                                      "only check available, and the extract just built is not from ",
                                      "the bytes that were read."),
                               fname, substr(get(fname, envir = .OUTCOME_PRE, inherits = FALSE), 1, 12),
                               substr(.h, 1, 12)))
      if (nrow(.c) == 1L && nzchar(.c$published_md5[1]) &&
          !identical(tolower(.h), tolower(.c$published_md5[1])))
        fail_loud_data(sprintf(paste0("01_mr_screen: %s no longer matches its published checksum ",
                                      "after being read (%s on disk, %s published). It was replaced ",
                                      "during the run; the extract just built is not from the ",
                                      "certified file."),
                               fname, substr(.h, 1, 12), substr(.c$published_md5[1], 1, 12)))
    }
    assign(fname, .h, envir = .OUTCOME_MD5)
  }
  get(fname, envir = .OUTCOME_MD5, inherits = FALSE)
}

build_outcome_lookup <- function(m) {
  f <- file.path(OUTCOME_DIR, m$file)
  if (!file.exists(f))
    fail_loud_data(sprintf("01_mr_screen: outcome %s expects %s, which is absent", oid, f))
  want <- c(m$rsid, m$ea, m$oa, m$beta, m$se, m$pval)
  if (!is.null(m$eaf)) want <- c(want, m$eaf); if (!is.null(m$n_col)) want <- c(want, m$n_col)
  for (.fb in c("rsid_fallback","ea_fallback","oa_fallback","beta_fallback","eaf_fallback")) if (!is.null(m[[.fb]])) want <- c(want, m[[.fb]])
  dt <- fread(f, sep=m$sep, select=unique(want), showProgress=FALSE)
  rid <- as.character(dt[[m$rsid]])
  if (!is.null(m$rsid_fallback)) { fb <- as.character(dt[[m$rsid_fallback]])
    # ANY unusable harmonised field triggers the row-level fallback, not just a missing hm_rsid.
    # 23,760 stroke rows carry a present hm_rsid beside an unusable hm_beta while the submitted
    # association is complete; 88 of them are in the current lookup set, and keying the fallback on
    # the identifier alone left those rows to be dropped by the beta filter and replaced by proxies.
    # B4GAT1 rs10896113 took proxy rs11227464 (beta -0.0007) over its own direct beta 0; NID2
    # rs2516600 took rs941622 (-0.0022) over 0; MSLN rs2235503 vanished from stroke entirely. The
    # stored MR table carries the consequences (B4GAT1 b=0.002218, NID2 b=-0.003820).
    .hb <- suppressWarnings(as.numeric(dt[[m$beta]]))
    # ANY of them, which is what the comment above says and what the replacement loop below assumes:
    # the harmonised and submitted fields are swapped as a set, so a row whose harmonised beta is
    # fine but whose harmonised ALLELES are unusable, or whose harmonised EAF is missing, was left
    # with a half-usable harmonised record and no way to reach its complete submitted one. SE and
    # p-value are deliberately not included -- LOCAL_META defines no fallback column for them, so
    # they are shared rather than part of the alternative set.
    na <- is.na(rid) | rid == "" | rid == "NA" | !is.finite(.hb)
    for (.p in list(c("ea","ea_fallback"), c("oa","oa_fallback"))) {
      if (is.null(m[[.p[1]]]) || is.null(m[[.p[2]]])) next
      .a <- toupper(as.character(dt[[m[[.p[1]]]]]))
      na <- na | is.na(.a) | !grepl("^[ACGT]+$", .a)
    }
    if (!is.null(m$eaf) && !is.null(m$eaf_fallback)) {
      .e <- suppressWarnings(as.numeric(dt[[m$eaf]]))
      na <- na | !is.finite(.e) | .e < 0 | .e > 1
    }
    rid[na] <- fb[na]
    for (.p in list(c("ea","ea_fallback"),c("oa","oa_fallback"),c("beta","beta_fallback"),c("eaf","eaf_fallback")))
      if (!is.null(m[[.p[1]]]) && !is.null(m[[.p[2]]])) { .v <- as.character(dt[[m[[.p[1]]]]]); .w <- as.character(dt[[m[[.p[2]]]]]); .v[na] <- .w[na]; set(dt, j=m[[.p[1]]], value=.v) } }
  if (isTRUE(m$multi_rsid)) { multi <- grepl(",", rid, fixed=TRUE)
    if (any(multi)) { idx <- which(multi)
      expl <- rbindlist(lapply(idx, function(i){ toks <- trimws(strsplit(rid[i],",")[[1]]); toks <- toks[toks %in% keep_set]; if(length(toks)) data.table(row=i, rs=toks) else NULL }))
      if (nrow(expl)) { add <- dt[expl$row]; dt <- rbind(dt, add); rid <- c(rid, expl$rs) } } }
  # THE FILE MUST BE IN THE NAMESPACE WE ARE QUERYING IT WITH. Everything downstream matches on
  # rsID. If this outcome's identifier column is in another namespace -- chr:pos, a different build's
  # rsIDs, an internal accession -- the match returns nearly nothing, every gene-outcome cell becomes
  # an ordinary `no_outcome_data`, the grid accounting completes, and an almost-empty outcome is
  # published as a normal result. guard_nonempty() needs one row and coverage needs one MR result, so
  # a single accidental match satisfies both.
  # This is not a coverage threshold and invents no number: it asks whether the column CONTAINS THE
  # KIND OF IDENTIFIER being looked up. A file of rsIDs answers yes whatever its overlap with our
  # instruments; a file of chr:pos answers no however many rows it has.
  .rs_like <- grepl("^rs[0-9]+$", rid)
  if (!any(.rs_like))
    fail_loud_data(sprintf(paste0("01_mr_screen: not one of %d identifiers in column '%s' of %s is ",
                                  "an rsID. The extraction matches on rsID, so this file is in a ",
                                  "different identifier namespace and would yield an almost empty ",
                                  "extract that every guard accepts."),
                           length(rid), m$rsid, m$file))
  k <- rid %in% keep_set; dt <- dt[k]; rid <- rid[k]
  lk <- data.table(rs=rid, ea=toupper(as.character(dt[[m$ea]])), oa=toupper(as.character(dt[[m$oa]])),
    beta=as.numeric(dt[[m$beta]]), se=as.numeric(dt[[m$se]]), pval=as.numeric(dt[[m$pval]]),
    eaf=if(!is.null(m$eaf)) suppressWarnings(as.numeric(dt[[m$eaf]])) else NA_real_,
    n=if(!is.null(m$n_col)) suppressWarnings(as.numeric(dt[[m$n_col]])) else m$n_fixed)
  # USABLE, not merely non-missing. `!is.na(beta) & !is.na(se) & se > 0` admits Inf and NaN (is.na()
  # is FALSE for Inf), a p-value outside [0,1], an EAF outside [0,1] and an allele that is not a
  # sequence. Every one of those survives into harmonisation and into an MR estimate: an infinite
  # beta or SE propagates through the IVW weights, an out-of-range EAF changes palindrome handling,
  # and a non-sequence allele silently fails to match the instrument's pair. These are properties of
  # the ROW, checked where the row is read, with no threshold invented -- each condition is simply
  # what the quantity must be to mean anything.
  # WHAT AN MR ESTIMATE ACTUALLY NEEDS is the association: a finite beta, a finite positive SE, and
  # alleles it can be harmonised on. The p-value, the effect-allele frequency and the per-SNP N are
  # ancillary -- IVW and Wald use none of them -- yet a row was DELETED when one of them was
  # malformed, taking a perfectly usable association with it. Malformed ancillary values are now set
  # to missing instead, so they cannot mislead a later consumer while the association survives.
  # (Ordering below still puts valid p-values before missing ones, so deduplication is unchanged.)
  .n0 <- nrow(lk)
  set(lk, i = which(!is.na(lk$pval) & !(is.finite(lk$pval) & lk$pval >= 0 & lk$pval <= 1)),
      j = "pval", value = NA_real_)
  set(lk, i = which(!is.na(lk$eaf)  & !(is.finite(lk$eaf)  & lk$eaf  >= 0 & lk$eaf  <= 1)),
      j = "eaf",  value = NA_real_)
  if (!is.null(m$n_col))
    set(lk, i = which(!is.na(lk$n) & !(is.finite(lk$n) & lk$n > 0)), j = "n", value = NA_real_)
  lk <- lk[is.finite(beta) & is.finite(se) & se > 0 &
           grepl("^[ACGT]+$", ea) & grepl("^[ACGT]+$", oa)]
  if (nrow(lk) < .n0)
    cat(sprintf("  [%s] dropped %d unusable row(s) (non-finite beta/se or a non-ACGT allele)
",
                m$file, .n0 - nrow(lk)))
  # ONE ROW PER (rsID, unordered allele pair) -- not one row per rsID. Deduplicating on the rsID alone
  # happens BEFORE anything knows which allele pair is wanted, so at a multiallelic site the row with
  # the smallest p-value won and the record describing the instrument's actual variant was discarded.
  # Measured on the shipped extracts: 43 direct rows across 16 proteins carry an allele pair
  # incompatible with their instrument, and for 23 of them the exact pair was present in the same raw
  # file and thrown away only for being less significant (C3/CHD rs11569466: A/G at p=0.118843 kept,
  # G/T at p=0.883646 discarded). Those cells were then counted as "directly present", which also
  # suppressed the proxy fallback, so they produced no MR result at all.
  lk[, .apair := paste(pmin(toupper(ea), toupper(oa)), pmax(toupper(ea), toupper(oa)))]
  { .k <- paste(lk$rs, lk$.apair)
    # Oriented onto one allele order first: a row written with ea/oa swapped carries -beta and is the
    # SAME association, not a contradiction. Comparing raw values would abort on a file that simply
    # lists both orientations.
    # A COLUMN, not an external vector. Inside j with by=, an external full-length vector is seen
    # whole by every group -- the same recycling error that aborted 03b's chromosome ranking.
    lk[, `:=`(.ckey = .k, .obeta = fifelse(toupper(ea) <= toupper(oa), beta, -beta))]
    .conf <- lk[, .(n = uniqueN(paste(signif(.obeta, 12), signif(se, 12)))), by = .ckey][n > 1L]
    setnames(.conf, ".ckey", "k")
    lk[, c(".ckey", ".obeta") := NULL]
    if (nrow(.conf))
      fail_loud_data(sprintf(paste0("01_mr_screen: %s has %d variant/allele-pair key(s) carrying ",
                                    "DIFFERENT beta/se (e.g. %s). Keeping the smallest p-value picks ",
                                    "the more extreme of two records that contradict each other."),
                             m$file, nrow(.conf), paste(head(.conf$k, 3), collapse = "; "))) }
  lk <- lk[order(pval)][!duplicated(paste(rs, .apair))]; setkey(lk, rs); lk }

# Hoisted out of read_local_outcome() so the proxy path can use the same map the direct path does.
.cmp <- c(A = "T", T = "A", C = "G", G = "C")
orient_proxy <- function(ie, pr, po) {
  p_inphase <- if (ie==pr$a_hap1) pr$b_hap1 else if (ie==pr$a_hap2) pr$b_hap2 else NA_character_
  if (is.na(p_inphase)) return(NULL)
  # BOTH panel alleles must be the outcome row's pair (see 02_extract_outcomes_local.R): matching only
  # the in-phase allele carries an A-versus-T beta as if it were A-versus-G, for the wrong allele
  # contrast and possibly with the wrong sign.
  p_other <- if (ie==pr$a_hap1) pr$b_hap2 else pr$b_hap1
  if (is.na(p_other)) return(NULL)
  # THE SAME VARIANT ON THE OTHER STRAND IS STILL THAT VARIANT -- the direct path above already
  # accepts the reverse-complement pair, and this one did not, so a valid opposite-strand proxy was
  # filed as "no compatible variant" or replaced by a lower-ranked one. The expected effect allele
  # has to be complemented too: relaxing only the pair test would compare the outcome's effect
  # allele with an uncomplemented in-phase allele and can carry the wrong sign.
  .pu <- toupper(c(p_inphase, p_other)); .ou <- toupper(c(po$ea, po$oa))
  # A PALINDROME CANNOT BE ORIENTED FROM ITS ALLELES ALONE, but it can be from FREQUENCIES -- which
  # is what harmonise_data(action=2) does for the direct path. Refusing all of them made the proxy
  # path strictly stricter than the direct one for the same class of variant (4,473 of 30,063 SNV
  # pairs in the panel cache are palindromic). The rule below is TwoSampleMR's, not a new one:
  # missing or strictly inside (0.42, 0.58) is uninformative; otherwise the two frequencies agreeing
  # on which side of 0.5 they fall means the literal mapping, disagreeing means the complement.
  #
  # PHASE ALONE IS NOT ENOUGH: it says which panel B allele travels with the chosen panel A allele,
  # not whether the outcome's A/T labels are the panel's strand or its reverse complement -- and
  # those give opposite signs. b_hap1 is B1 under coupling and B2 under repulsion, but this cache is
  # built at r2 >= 0.8 and an informative MAF is at most 0.42, so repulsion caps r2 at
  # 0.42/0.58 = 0.724 < 0.8: any palindromic row this rule can accept is in coupling, and b_hap1
  # therefore carries frequency MAF_B.
  .unin <- function(f) is.na(f) || (f > 0.42 && f < 0.58)
  .pal_pair <- function(a, b) a %in% names(.cmp) && identical(unname(.cmp[[a]]), b)
  if (.pal_pair(.pu[1], .pu[2])) {
    # STILL THE SAME VARIANT. The frequency rule decides STRAND; it does not decide whether this
    # outcome row is the proxy at all. For a palindrome the literal and complemented pairs are the
    # same set, so one setequal() covers both -- but it has to be made.
    if (!setequal(.pu, .ou)) return(NULL)
    .fb <- suppressWarnings(as.numeric(pr$MAF_B))
    if (.unin(.fb) || .unin(suppressWarnings(as.numeric(po$eaf)))) return(NULL)
    .f_ea <- if (identical(toupper(po$ea), toupper(pr$b_hap1))) .fb else 1 - .fb
    .expect <- if ((.f_ea < 0.5) == (as.numeric(po$eaf) < 0.5)) toupper(p_inphase)
               else unname(.cmp[[toupper(p_inphase)]])
    return(if (toupper(po$ea) == .expect) list(beta=po$beta,eaf=po$eaf,se=po$se,n=po$n)
           else list(beta=-po$beta,eaf=if(is.na(po$eaf)) NA_real_ else 1-po$eaf,se=po$se,n=po$n))
  }
  .expect <- toupper(p_inphase)
  if (!setequal(.pu, .ou)) {
    if (!all(.pu %in% names(.cmp)) || !setequal(unname(.cmp[.pu]), .ou)) return(NULL)
    .expect <- unname(.cmp[toupper(p_inphase)])
  }
  if (toupper(po$ea)==.expect) list(beta=po$beta,eaf=po$eaf,se=po$se,n=po$n) else list(beta=-po$beta,eaf=if(is.na(po$eaf)) NA_real_ else 1-po$eaf,se=po$se,n=po$n) }

read_local_outcome <- function(oid) {
  m <- LOCAL_META[[oid]]; lk <- build_outcome_lookup(m)
  # A DIRECT match is the same VARIANT, not merely the same name. The outcome row's unordered allele
  # pair has to be the instrument's; otherwise the beta describes a different contrast at that
  # position. An rsID present under some other pair is NOT directly available, so it now falls
  # through to the proxy path instead of being silently accepted and blocking it.
  .ia <- inst_al[.(intersect(instr_rsids, unique(lk$rs)))]
  .ia <- .ia[!is.na(ie) & !is.na(io)]
  .ia[, .apair := paste(pmin(toupper(ie), toupper(io)), pmax(toupper(ie), toupper(io)))]
  # The same variant on the other strand is still that variant (localize/02 carries this too).
  # .cmp is defined once, above orient_proxy(), so the direct and proxy paths cannot disagree on it.
  .rc <- .ia[grepl("^[ACGT]$", toupper(ie)) & grepl("^[ACGT]$", toupper(io))]
  if (nrow(.rc)) .rc[, .apair := paste(pmin(.cmp[toupper(ie)], .cmp[toupper(io)]), pmax(.cmp[toupper(ie)], .cmp[toupper(io)]))]
  dsub <- merge(lk, unique(rbind(.ia[, .(rs = SNP, .apair)], .rc[, .(rs = SNP, .apair)])), by = c("rs", ".apair"), sort = FALSE)
  dsub <- dsub[order(pval)][!duplicated(rs)]
  # A DIRECT ROW THAT HARMONISATION WILL DISCARD IS NOT A DIRECT ROW. A palindromic pair whose EAF
  # cannot resolve the strand is dropped later by harmonise_data(action=2) as ambiguous -- but by
  # then the instrument has been removed from `missing`, so no proxy was tried and the cell yields
  # nothing at all. Classifying it as not-directly-available lets the proxy path run; if no proxy
  # orients either, the cell is reported as having no compatible variant, which is what it is.
  .pe <- toupper(dsub$ea); .po <- toupper(dsub$oa)
  .pal <- .pe %chin% names(.cmp) & .po %chin% names(.cmp) & .po == unname(.cmp[.pe])
  # EITHER SIDE. harmonise_data(action=2) calls a palindrome ambiguous when the exposure OR the
  # outcome frequency is uninformative, so testing the outcome alone left rows that it would drop
  # classified as directly available -- and the proxy path is skipped for exactly those.
  .ex <- inst_al[.(dsub$rs)]$eaf_exp
  .amb <- .pal & (is.na(dsub$eaf) | (dsub$eaf > 0.42 & dsub$eaf < 0.58) |
                  is.na(.ex)      | (.ex      > 0.42 & .ex      < 0.58))
  if (any(.amb)) {
    cat(sprintf("  [%s] %d direct palindromic row(s) whose EAF cannot resolve strand -> proxy path
",
                m$file, sum(.amb)))
    dsub <- dsub[!.amb]
  }
  d_ids <- dsub$rs; in_out <- unique(lk$rs)
  direct <- if (length(d_ids)) data.table(SNP=d_ids, effect_allele.outcome=dsub$ea, other_allele.outcome=dsub$oa,
    beta.outcome=dsub$beta, se.outcome=dsub$se, pval.outcome=dsub$pval, eaf.outcome=dsub$eaf,
    samplesize.outcome=dsub$n, ncase.outcome=if(!is.null(m$ncase)) m$ncase else NA_integer_,
    proxy.outcome=FALSE, proxy_snp.outcome=NA_character_, proxy_r2.outcome=NA_real_) else NULL
  # `missing` is now every instrument WITHOUT a compatible direct row, which includes those present
  # under another allele pair -- exactly the cells the old dedupe silently closed off.
  missing <- setdiff(instr_rsids, d_ids); pr_av <- prox[instrument %in% missing & proxy %in% in_out]; rows <- list()
  for (g in unique(pr_av$instrument)) { cands <- pr_av[instrument==g][order(-r2)]; ie <- inst_al[.(g)]$ie; io <- inst_al[.(g)]$io
    for (j in seq_len(nrow(cands))) { .pos <- lk[.(cands$proxy[j])]
      # The proxy's rsID can now hold several allele pairs too; orient_proxy() takes ONE row and
      # decides whether that record's pair is the panel pair, so each is offered in turn.
      # po must be the accepted record: .pos[1] paired one contrast's p-value with another's beta.
      or <- NULL; po <- .pos
      for (.k in seq_len(nrow(.pos))) { or <- orient_proxy(ie, cands[j], .pos[.k]); if (!is.null(or)) { po <- .pos[.k]; break } }
      if (!is.null(or)) { rows[[length(rows)+1]] <- data.table(SNP=g, effect_allele.outcome=ie, other_allele.outcome=io,
        beta.outcome=or$beta, se.outcome=or$se, pval.outcome=po$pval, eaf.outcome=or$eaf, samplesize.outcome=or$n,
        ncase.outcome=if(!is.null(m$ncase)) m$ncase else NA_integer_, proxy.outcome=TRUE, proxy_snp.outcome=cands$proxy[j], proxy_r2.outcome=cands$r2[j]); break } } }
  proxied <- if (length(rows)) rbindlist(rows) else NULL
  out <- rbind(direct, proxied); if (is.null(out)||!nrow(out)) return(NULL)
  # WHERE EACH INSTRUMENT WENT, REPORTED -- NOT ASSERTED.
  # This was a fatal guard on the premise that an instrument identifier present in `lk` and
  # absent from `out` had been lost by code. The premise is false, and measurably so: `lk` is
  # keyed on the identifier, while a variant is the identifier AND its allele pair. A multiallelic
  # rsID can be in the outcome file under a different pair, which is a different variant and
  # carries a different contrast. On the rebuilt CHD extraction six instruments are exactly that
  # -- EIF4G3 rs11444926 (instrument GA/G, FinnGen GAA/G), VEGFC rs142109874, SERPINB5
  # rs148729762, PBXIP1 rs35704440, IL17F rs66834286, SEZ6 rs778315867 -- all with
  # rsid_ld_panel = NA, so no proxy can stand in. The guard reported them as code loss and stopped
  # the run on the FIRST outcome, so the rebuild could not start.
  #
  # And the repaired version would have been worse than useless: `direct` is built from every row
  # of `dsub`, so "in dsub, absent from out" is empty by construction. Any version of this test at
  # this point compares the output of a step with its own input. What actually separates the three
  # outcomes is a property of the DATA -- whether the file holds this variant, whether a proxy for
  # it survived -- so it is counted and printed, and the run continues.
  # The failure this was reaching for -- an almost-empty outcome published as a complete screen --
  # is not detectable here at all: it needs the eight outcomes compared against each other, which
  # is where the per-outcome asymmetry check below does it.
  .n_direct  <- length(d_ids)
  .n_proxied <- if (!is.null(proxied)) uniqueN(proxied$SNP) else 0L
  .n_absent  <- length(setdiff(instr_rsids, unique(as.character(out$SNP))))
  cat(sprintf("[%s direct / %s proxy / %s no compatible variant] ",
              format(.n_direct, big.mark = ","), format(.n_proxied, big.mark = ","),
              format(.n_absent, big.mark = ",")))
  out[, `:=`(outcome=oid, id.outcome=oid, mr_keep.outcome=TRUE)]; out[] }

# ---- (C) offline MR (mirror localize/04) ----
.SEC_MR_FAIL <- data.table()   # 2026-07-01 Codex-phase2 #6: audit trail for secondary-MR (Egger/WM/Q) ERRORS -- NOT a silent NULL
.sec_mr_note <- function(gene, oid, what, msg) assign(".SEC_MR_FAIL", rbind(get(".SEC_MR_FAIL", envir=.GlobalEnv), data.table(gene=gene, outcome_id=oid, what=what, msg=substr(msg,1,200))), envir=.GlobalEnv)
mr_one <- function(exp_fmt, out_dat, gene, oid) {
  # DIRECT AND PROXY ROWS ARE HARMONISED SEPARATELY, because they carry different kinds of evidence
  # about which allele is the effect allele. orient_proxy() above has already established a proxy
  # row's orientation from LD phase and verified both of its alleles -- but it then attached the
  # PROXY variant's frequency while relabelling the alleles as the instrument's. harmonise_data()
  # with action=2 reads that as a palindromic frequency conflict and drops the instrument: measured
  # on the rebuilt set, 7 instrument-outcome cells are discarded this way (CNTNAP2-Stroke
  # rs2373282; CXCL16-Stroke/CKD/eGFR rs78086290; WFDC12-Stroke/CKD/eGFR rs75002084), each one
  # changing the beta, SE and p-value written to the MR table for that pair.
  # eaf.outcome is set to NA rather than to a substitute frequency -- no frequency for the
  # instrument's own alleles exists in the outcome study -- and those rows go through action=1,
  # which still refuses incompatible alleles but does not adjudicate by frequency. action=2 for the
  # direct rows is unchanged. Phase-orienting the proxy frequency (0 of 7 retained) and borrowing
  # the 1000G index frequency (1 of 7) were both measured and rejected; this retains 7 of 7.
  # ON THE PREVIOUS GENERATION this changes five cells and grows the MR table from 11,538 to 11,542
  # rows (recovering CXCL16-CKD/eGFR/Stroke and GCHFR-Stroke, and moving CNTNAP2-Stroke to its
  # three-SNP estimate). No significance set moves: Bonferroni 75 -> 75, by-method FDR 234 -> 234
  # with identical members, suggestive 1,368 -> 1,368. It is a change to a published artifact and is
  # to be disclosed as one.
  .is_prox <- if ("proxy.outcome" %in% names(out_dat)) !is.na(out_dat$proxy.outcome) & out_dat$proxy.outcome else rep(FALSE, nrow(out_dat))
  .hz <- function(d, act) if (!nrow(d)) NULL else
    tryCatch(harmonise_data(exp_fmt, d, action=act), error=function(e) fail_loud_data(sprintf("01_mr_screen MR harmonise_data THREW (%s ~ %s, action=%d): %s -- malformed data, not a legit no-overlap", gene, oid, act, conditionMessage(e))))
  .direct <- .hz(out_dat[!.is_prox, , drop=FALSE], 2L)
  .pd <- out_dat[.is_prox, , drop=FALSE]
  if (nrow(.pd)) .pd$eaf.outcome <- NA_real_
  .proxy <- .hz(.pd, 1L)
  harm <- if (is.null(.direct)) .proxy else if (is.null(.proxy)) .direct else
            as.data.frame(data.table::rbindlist(list(.direct, .proxy), fill=TRUE))
  if (is.null(harm)||nrow(harm)==0||sum(harm$mr_keep)==0) return(NULL)
  n_snps <- sum(harm$mr_keep); pm <- if (n_snps==1) "mr_wald_ratio" else "mr_ivw"
  r <- tryCatch(mr(harm, method_list=pm), error=function(e) fail_loud_data(sprintf("01_mr_screen MR mr() THREW (%s ~ %s): %s -- computation failure", gene, oid, conditionMessage(e)))); if (is.null(r)||nrow(r)==0) return(NULL)  # 2026-07-01 (closure sweep): thrown mr() = fatal
  # A RETURNED ROW IS NOT AN ESTIMATE. mr() returning something non-empty was the whole test, so a
  # primary row carrying NA, NaN or Inf for b, se or pval counted as a filled cell, satisfied the
  # coverage audit and was published. Fatal rather than dropped: an unusable number here is a
  # computation that went wrong, not a gene-outcome pair with no data.
  { .pr <- as.data.table(r)[method %in% c("Inverse variance weighted", "Wald ratio")]
    # EXACTLY ONE. "if any primary row exists" let a result whose method label is neither of these
    # two bypass the check entirely and still be saved as a filled cell.
    if (nrow(.pr) != 1L || !(.pr$method[1] %in% c("Inverse variance weighted", "Wald ratio")))
      fail_loud_data(sprintf(paste0("01_mr_screen: mr() returned %d primary row(s) for %s ~ %s ",
                                    "(methods: %s). One IVW or Wald-ratio row is what every count ",
                                    "downstream assumes."),
                             nrow(.pr), gene, oid, paste(unique(as.data.table(r)$method), collapse = ", ")))
    if (nrow(.pr) && !all(is.finite(.pr$b) & is.finite(.pr$se) & .pr$se > 0 &
                          is.finite(.pr$pval) & .pr$pval >= 0 & .pr$pval <= 1))
      fail_loud_data(sprintf(paste0("01_mr_screen: mr() returned a non-finite primary estimate for ",
                                    "%s ~ %s (b=%s, se=%s, p=%s). A cell like this counts as a ",
                                    "result everywhere downstream."),
                             gene, oid, paste(.pr$b, collapse = ","), paste(.pr$se, collapse = ","),
                             paste(.pr$pval, collapse = ","))) }
  hk <- as.data.table(harm)[mr_keep==TRUE]; n_prox <- if ("proxy.outcome" %in% names(hk)) sum(hk$proxy.outcome, na.rm=TRUE) else NA_integer_
  # OR ONLY WHERE AN ODDS RATIO EXISTS. exp(b) was written into columns named OR for every outcome,
# including the continuous ones (eGFR, and any other type != "cc"), where it is an exponentiated
# SD-scale effect and not an odds ratio at all. b, se and pval were always correct; the three OR
# columns were mislabelled. They are NA for continuous outcomes rather than silently wrong.
  .cc <- identical(LOCAL_META[[oid]]$type, "cc")
  .fin <- function(x) fifelse(is.finite(x), x, NA_real_)   # exp() overflows to Inf, which is not an odds ratio
  res <- as.data.table(r)[, .(gene, outcome_id=oid, n_snps, method, b, se, pval, OR=if (.cc) .fin(exp(b)) else NA_real_, OR_lo95=if (.cc) .fin(exp(b-1.96*se)) else NA_real_, OR_up95=if (.cc) .fin(exp(b+1.96*se)) else NA_real_, n_proxy=n_prox, any_proxy=isTRUE(n_prox>0))]
  if (n_snps>=3) { ex <- tryCatch(mr(harm, method_list=c("mr_egger_regression","mr_weighted_median")), error=function(e){ .sec_mr_note(gene, oid, "egger/wm", conditionMessage(e)); NULL })   # 2026-07-01 Codex-phase2 #6: audit (not silent) secondary-MR errors
    # A RETURNED SECONDARY ROW IS NOT A SECONDARY RESULT EITHER. mr() has non-throwing null paths --
    # Egger returns an all-missing row on collinear input -- so .SEC_MR_FAIL could report no error
    # while an unusable row was appended to the published table. Recorded as the failure it is.
    if (!is.null(ex) && nrow(ex) > 0) {
      .bad <- as.data.table(ex)[!(is.finite(b) & is.finite(se) & se > 0 &
                                  is.finite(pval) & pval >= 0 & pval <= 1)]
      if (nrow(.bad)) { .sec_mr_note(gene, oid, paste(.bad$method, collapse = "/"),
                                     "returned a non-finite or out-of-range estimate")
                        ex <- as.data.table(ex)[!(method %in% .bad$method)] } }
    if (!is.null(ex)&&nrow(ex)>0) res <- rbind(res, as.data.table(ex)[, .(gene, outcome_id=oid, n_snps, method, b, se, pval, OR=if (.cc) .fin(exp(b)) else NA_real_, OR_lo95=if (.cc) .fin(exp(b-1.96*se)) else NA_real_, OR_up95=if (.cc) .fin(exp(b+1.96*se)) else NA_real_, n_proxy=n_prox, any_proxy=isTRUE(n_prox>0))], fill=TRUE)
    { .want <- c("MR Egger", "Weighted median")
      .got  <- if (!is.null(ex) && nrow(ex) > 0) as.data.table(ex)$method else character(0)
      .miss <- setdiff(.want, .got)
      # Requested and not returned is a failure too. mr() has non-throwing paths that return fewer
      # rows than methods asked for, and the audit only recorded thrown errors.
      if (length(.miss)) .sec_mr_note(gene, oid, paste(.miss, collapse = "/"), "requested but not returned") }
    het <- tryCatch(mr_heterogeneity(harm, method_list="mr_ivw"), error=function(e){ .sec_mr_note(gene, oid, "heterogeneity", conditionMessage(e)); NULL })
    if (!is.null(het)&&nrow(het)>0) res[method=="Inverse variance weighted", `:=`(Q=het$Q, Q_df=het$Q_df, Q_pval=het$Q_pval)] }
  res }

if (anyDuplicated(unname(OUTCOMES)))
  stop("01_mr_screen: OUTCOMES lists a duplicate outcome id (", 
       paste(unname(OUTCOMES)[duplicated(unname(OUTCOMES))], collapse = ", "),
       "). The Bonferroni denominator counts entries, so a duplicate makes every threshold stricter ",
       "than the number of tests actually performed.")
# THE OUTCOMES TESTED AND THE OUTCOMES THE THRESHOLD COUNTS MUST BE THE SAME SET. Extraction and
# coverage run over names(LOCAL_META); BONFERRONI_ALPHA is 0.05 / (N_PROTEINS_CANONICAL *
# length(OUTCOMES)). Adding an outcome to config while LOCAL_META still lists eight leaves the screen
# testing eight and the threshold dividing by nine -- every p-value judged against a bar for an
# analysis that was not run.
# Compared as OUTCOME IDS on both sides. LOCAL_META is keyed by id (FinnGen_R10_I9_CHD) and OUTCOMES
# by name (CAD), so comparing names(LOCAL_META) with names(OUTCOMES) is always FALSE -- the same
# two-namespace mistake as the proxy/instrument comparison above, caught here only because the check
# was run against the real config before it was sent anywhere.
if (!setequal(names(LOCAL_META), unname(OUTCOMES)))
  stop("01_mr_screen: this script extracts ", length(LOCAL_META), " outcome(s) (",
       paste(sort(setdiff(names(LOCAL_META), unname(OUTCOMES))), collapse = ", "),
       " not in config; ", paste(sort(setdiff(unname(OUTCOMES), names(LOCAL_META))), collapse = ", "),
       " in config but not extracted), while BONFERRONI_ALPHA divides by length(OUTCOMES) = ",
       length(OUTCOMES), ". The tested grid and the threshold would describe different analyses.")

# CLEARED AT THE START OF THE RUN. .CKM_DATA_FAILURE is a session-global set by a caught MR
# failure; nothing ever removed it, so in a reused R session a flag left by an earlier run made
# the first otherwise-successful cell abort with a stale message about a different pair.
if (exists(".CKM_DATA_FAILURE", envir = .GlobalEnv)) rm(".CKM_DATA_FAILURE", envir = .GlobalEnv)

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_guard.R"))   # 防呆 fail-loud guards (no silent NULL/empty/missing)
cat(sprintf("\nExtracting %d outcomes for %d instruments (%d genes) ...\n",
            length(LOCAL_META), nrow(inst), uniqueN(inst$gene))); flush(stdout())
# PUBLICATION MODEL, stated so the next person does not build a transaction here.
# Each extract is written to its final path as it is produced, and the MR table plus its three
# significance subsets are written after every guard has passed. That is FAIL-CLOSED PUBLICATION,
# not an atomic transaction: an interrupted run can leave five new extracts beside three from an
# earlier attempt, or a new MR table beside older subsets.
# What makes that safe is downstream, not here. Every artifact carries ckm_phase1_md5 (and now
# ckm_outcome_md5), and every consumer that combines them compares the stamps and refuses a mismatch
# -- 53 refuses unstamped or foreign extracts, 02_promote checks the MR table, 41_3 checks the proxy
# cache. So an interrupted run costs AVAILABILITY: the run must be repeated. It cannot silently mix
# generations, because nothing downstream will accept the mixture.
# The same question was put to review for cascade/03b, whose answer applies here: there is no smaller
# correct transactional repair inside a multi-file publication step, and rollback machinery adds
# destructive states without improving integrity. DO NOT add it.
extract_dir <- LOCAL_EXTRACT_DIR; dir.create(extract_dir, showWarnings=FALSE, recursive=TRUE)
# Checked here rather than discovered at the first saveRDS, which is after PLINK and after the
# whole first outcome has been read.
if (!dir.exists(extract_dir))
  stop("01_mr_screen: cannot create the outcome extract directory ", extract_dir)
out_list <- list()
for (oid in names(LOCAL_META)) { cat(" ", oid, " "); flush(stdout())
  .vs <- .verify_outcome_source(LOCAL_META[[oid]]$file, oid)
  cat(if (identical(.vs, "verified")) "[src-md5 ok] " else sprintf("[src-md5 %s] ", .vs))
  assign(LOCAL_META[[oid]]$file,
         unname(tools::md5sum(file.path(OUTCOME_DIR, LOCAL_META[[oid]]$file))),
         envir = .OUTCOME_PRE)
  o <- read_local_outcome(oid)
  guard_nonempty(o, paste0("outcome:", oid))   # a whole outcome NULL/empty = silent bug -> ABORT (before saving)
  # WHICH GWAS THIS IS, by content. Only the FILENAME was checked, and the extract carried nothing
  # but the instrument hash -- so a different or edited sumstats file placed at the expected path is
  # assigned the hardcoded outcome id, passes every guard, and its estimates are stamped as
  # canonical. The file that was actually read is named and hashed here.
  # THE PRODUCER, TOO. The stamps below name the inputs; they do not name what turned those
  # inputs into this file. Two runs of DIFFERENT extraction or proxy code over the same Phase 1
  # and the same outcome files therefore carry identical stamps, and an extract left by the
  # earlier one is accepted beside artifacts from the later one.
  attr(o, "ckm_entry_script_md5") <- SELF_MD5
  attr(o, "ckm_run_id") <- CKM_RUN_ID
  attr(o, "ckm_phase1_md5")   <- unname(INSTRUMENTS_MD5)
  attr(o, "ckm_phase1_file")  <- .p1_file
  # WHAT THIS HASH DOES AND DOES NOT ESTABLISH, said plainly because the distinction matters.
  # It is computed here, from whatever is at the expected path. So it DETECTS CHANGE -- the same
  # outcome file hashing differently in a later run, or an extract and an MR table disagreeing about
  # which bytes they came from. It does NOT establish IDENTITY: nothing local certifies that
  # finngen_R10_I9_CHD.gz is the CHD phenotype. The FinnGen manifest is a phenotype catalogue with
  # no checksums, these files carry no per-row N to compare against config's expected sample sizes,
  # and only shah2020_HF (not in the current eight) ships an .md5sum.txt. A file swapped before the
  # first run of a generation would be hashed and stamped as canonical.
  # Closing that needs a certificate this project does not have -- recorded in REBUILD_LEDGER rather
  # than papered over with a check that would only re-assert what this line already assumes.
  attr(o, "ckm_outcome_file") <- LOCAL_META[[oid]]$file
  # ONE HASH PER FILE PER RUN, taken when the file is read and reused wherever it is recorded. It
  # used to be computed separately here and again for the MR table, so a file replaced between the
  # two received two different identities: the extract said version A, the MR table built from that
  # very extract said version B, and every downstream stamp comparison would then be satisfied by
  # whichever it happened to consult.
  attr(o, "ckm_outcome_md5")  <- .outcome_md5(LOCAL_META[[oid]]$file)
  # WHETHER THIS OUTCOME WAS CERTIFIED BY ITS PROVIDER, recorded rather than only printed. Three of
  # the eight sources publish no checksum, so refusing them would abort the run on legitimate data;
  # what must not happen is an uncertified outcome being indistinguishable downstream from a
  # certified one, which is the case while the status exists only in the console.
  attr(o, "ckm_outcome_certified") <- .vs
  saveRDS(o, file.path(extract_dir, sprintf("%s.rds", oid))); out_list[[oid]] <- o
  cat(sprintf("rows=%d (proxy=%d)\n", nrow(o), sum(o$proxy.outcome, na.rm=TRUE))); flush(stdout()) }
out_all <- rbindlist(out_list, fill=TRUE)
guard_exact_set(out_all$id.outcome, names(LOCAL_META), paste(PHASE1_TAG, "rerun outcome ids"))   # all 8 present, no stale extras, or ABORT
cat(sprintf("total outcome rows: %d (%d proxied)\n", nrow(out_all), sum(out_all$proxy.outcome, na.rm=TRUE)))

set.seed(2026)   # NEW-4 fix: deterministic mr_weighted_median bootstrap SE
genes <- inst[, unique(gene)]; oids <- unique(out_all$id.outcome); all_res <- list(); attempted_cells <- list(); no_data_cells <- list(); null_exp_genes <- character()   # gate obj #5 + P1-F01: cell-level (gene x outcome) audit incl no-outcome-data + null-exposure genes
for (i in seq_along(genes)) { g <- genes[i]; pdt <- inst[gene==g]
  exp_fmt <- tryCatch(format_data(as.data.frame(pdt), type="exposure", snp_col="SNP", beta_col="beta", se_col="se",
    effect_allele_col="effect_allele", other_allele_col="other_allele", eaf_col="eaf", pval_col="pval", samplesize_col="samplesize", phenotype_col="exposure"), error=function(e) fail_loud_data(sprintf("01_mr_screen MR format_data(exposure) THREW for gene %s: %s -- malformed instrument, fatal", g, conditionMessage(e))))  # 2026-07-01 (closure sweep)
  if (is.null(exp_fmt)) { null_exp_genes <- c(null_exp_genes, g); next }   # P1-F01: track null-exposure genes (else their 8 cells vanish from grid accounting)
  for (oid in oids) { od <- out_all[id.outcome==oid & SNP %in% pdt$SNP]
    if (nrow(od)==0) { no_data_cells[[length(no_data_cells)+1]] <- data.table(gene=g, outcome_id=oid); next }   # P1-F01 (cold audit): record no-outcome-data cells (were silently skipped -> full-grid accounting was incomplete)
    attempted_cells[[length(attempted_cells)+1]] <- data.table(gene=g, outcome_id=oid)   # gate obj #5: cell attempted (>=1 instrument SNP present in outcome)
    r <- tryCatch(mr_one(exp_fmt, as.data.frame(od), g, oid), error=function(e) { assign(".CKM_DATA_FAILURE", sprintf("01_mr_screen MR %s ~ %s: %s", g, oid, conditionMessage(e)), envir=.GlobalEnv); NULL })  # 2026-07-01 (closure sweep): flag any thrown MR failure (incl fail_loud_data from harmonise/mr)
    if (exists(".CKM_DATA_FAILURE", envir=.GlobalEnv)) stop(paste("01_mr_screen ABORT -- MR computation failure (thrown, not a legit no-overlap skip):", get(".CKM_DATA_FAILURE", envir=.GlobalEnv)))
    if (!is.null(r)&&nrow(r)>0) all_res[[length(all_res)+1]] <- r }
  if (i %% 200 == 0) cat(sprintf("  [%d/%d] %s | rows=%d | %.0fs\n", i, length(genes), g, length(all_res), proc.time()["elapsed"]-t_start)) }

mr_all <- rbindlist(all_res, fill=TRUE)
# If NO cell produced a result, mr_all has no columns at all, and the next line dies on a missing
# `outcome_id` before guard_coverage can say what happened. That case is the most severe version of
# exactly what guard_coverage exists to catch, so it must be named here.
if (!nrow(mr_all))
  fail_loud_data(sprintf(paste0("01_mr_screen: MR produced no result for ANY gene-outcome cell ",
                                "(%d genes x %d outcomes attempted). This is total harmonisation or ",
                                "extraction failure, not a coverage gap."),
                         length(genes), length(oids)))
mr_all[, outcome_name := names(OUTCOMES)[match(outcome_id, OUTCOMES)]]
# 2026-06-30 (gate obj #1): a whole outcome producing ZERO MR results = systematic harmonization/source
# failure. The cascade audit starts from the suggestive subset and would NOT see it -> assert here.
guard_coverage(unique(mr_all$outcome_id), names(LOCAL_META), paste(PHASE1_TAG, "MR: every outcome has >=1 result"))
N_PROTEINS <- uniqueN(inst$gene); BONF <- 0.05/(N_PROTEINS*length(OUTCOMES))
mr_all[, bonferroni_sig := pval < BONF]
# `fdr_sig` is BH WITHIN each method family, over ANALYSABLE rows only -- deliberately NOT the
# manuscript's primary FDR, which is BH ONCE over the full N_PROTEINS x length(OUTCOMES) grid with
# untested cells padded at p=1 (make_draft_manuscript_redesigned.py:84-91; recomputed locally at
# make_analysis_reports.R:58). Alongside it the reports print a "僅可分析 sensitivity" value -- but that
# one is computed FRESH by the generator (pooled BH, after the primary dedup), NOT read from the
# artifact this line writes. The two are a documented pair, established by cert r7 #2 on 2026-07-05
# -- see make_analysis_reports.R:50-51.
#
# CAUTION (2026-07-16): the artifact THIS line writes is NOT the same computation as the "234" the
# manuscript prints. This line is BH within each method family, BEFORE the primary dedup. The generator
# computes its sensitivity value fresh as POOLED BH over the 11,534 analysable rows, AFTER the dedup.
# The two select the IDENTICAL SET on the frozen data (setequal, verified 2026-07-16) -- an observed fact,
# not a guarantee: 0 pairs carry both Wald and IVW, and none falls in the 1.9% gap between the two
# denominators. Do not describe this artifact as "the sensitivity value the manuscript reports".
# DO NOT "fix" this line to match the primary definition: doing so collapses this artifact onto the
# primary set, and make_analysis_reports.R:174 would then print the same number twice while calling
# one of them a sensitivity analysis. The real defect is downstream -- phase3/28, phase3/32 and
# cascade/14_assign_tiers.R consumed THIS artifact while describing it as "MR FDR q<0.05" -- using the
# saved by-method definition [2] and calling it the primary one [1]. Fix them, not this.
# STATUS 2026-07-16: 28 and 32 are FIXED -- they read mr_all_results.rds and derive the primary rule via
# primary_fdr_set() (R/00_config.R). 14_assign_tiers.R still reads this artifact, deliberately: its use
# is provably inert on the frozen 52 and fixing it means changing assign_evidence_tiers()'s signature in
# theme_paper3.R, which carries a byte contract. Tracked in REBUILD_LEDGER §2.1c.
mr_all[, fdr_sig := p.adjust(pval, method="BH") < FDR_ALPHA, by=method]
# Stamped with the instrument set it was computed from, for the same reason the outcome extracts are:
# 53 reads this file's recorded p-values alongside instruments and LD it recomputes itself, and
# without a stamp an older MR run is indistinguishable from a current one.
# WHICH OUTCOMES WERE CERTIFIED, carried onto the table people actually read. It is attached to each
# extract, but rbindlist() drops table-level attributes, so the promoted tables disclosed nothing --
# and 02 does not read the extracts.
attr(mr_all, "ckm_outcome_certified") <- vapply(names(LOCAL_META), function(.o) {
  .a <- attr(out_list[[.o]], "ckm_outcome_certified", exact = TRUE)
  if (is.null(.a)) NA_character_ else as.character(.a) }, character(1))
attr(mr_all, "ckm_entry_script_md5") <- SELF_MD5
attr(mr_all, "ckm_run_id") <- CKM_RUN_ID
attr(mr_all, "ckm_declared_partial_env") <- paste(
  c(paste0("R=", getRversion()),
    vapply(c("TwoSampleMR", "data.table", "coloc"),
           function(.p) tryCatch(paste0(.p, "=", as.character(utils::packageVersion(.p))),
                                 error = function(e) paste0(.p, "=absent")), character(1))),
  collapse = "; ")
attr(mr_all, "ckm_phase1_md5")  <- unname(INSTRUMENTS_MD5)
attr(mr_all, "ckm_phase1_file") <- PHASE1_INSTRUMENTS
# ...and the eight outcome files these estimates came from, by content, for the same reason.
attr(mr_all, "ckm_outcome_md5") <- vapply(names(LOCAL_META),
  function(.o) .outcome_md5(LOCAL_META[[.o]]$file), character(1))
# NOT SAVED HERE. The cell-accounting and cell-loss guards below can still stop this run, and
# writing the full MR table before them meant a refused run left a new all-results file on disk
# beside the PREVIOUS run's three significance subsets -- a mixed set that looks like one
# analysis. It is written together with those subsets, after every guard has passed.
# P1-T025 fix: genes that entered MR but produced no result (silent-skip audit trail).
# COMPUTED here because the guards below need these sets; WRITTEN after them, for the reason given
# at the saveRDS block: a run refused by a guard used to leave new audit CSVs beside the previous
# run's MR tables, and nothing in the CSVs says which run they belong to.
.no_res_genes <- data.table(gene = sort(setdiff(genes, unique(mr_all$gene))))
# 2026-06-30 gate obj #5: cell-level audit — (gene,outcome) attempted (>=1 instrument SNP present in outcome) but NO MR result
.att <- unique(rbindlist(attempted_cells)); .att[, .k := paste(gene, outcome_id)]
.gotc <- unique(mr_all[, paste(gene, outcome_id)])
.no_res_cells <- .att[!.k %in% .gotc, .(gene, outcome_id)]
# P1-F01 (cold audit): also persist no-OUTCOME-DATA cells (gene had an instrument but the outcome GWAS had no
# row) so the full N_PROTEINS_CANONICAL x 8 grid accounting is complete: with-result + attempted-no-result + no-outcome-data = grid.
.nd <- if (length(no_data_cells)) unique(rbindlist(no_data_cells)) else data.table(gene=character(), outcome_id=character())
# (.nd itself is written with the others, after the guards.)
.grid <- uniqueN(inst$gene) * length(oids)
.acc <- length(.gotc) + nrow(.att[!.k %in% .gotc]) + nrow(.nd) + length(unique(null_exp_genes)) * length(oids)
cat(sprintf("MR cells: %d with-result, %d attempted-no-result, %d no-outcome-data, %d null-exposure-gene | full-grid %d/%d %s\n",
            length(.gotc), nrow(.att[!.k %in% .gotc]), nrow(.nd), length(unique(null_exp_genes)) * length(oids), .acc, .grid,
            if (.acc == .grid) "COMPLETE" else "INCOMPLETE"))
# P1-F01 / P5-NEW-01 (cold audit): every grid cell must be classified (with-result | attempted-no-result |
# no-outcome-data | null-exposure) -> fail loud rather than leave a partial no-result provenance surface.
if (.acc != .grid) stop(sprintf("01: MR cell accounting INCOMPLETE (%d/%d) -- %d cells unclassified; refusing to leave a partial no-result provenance surface", .acc, .grid, .grid - .acc))
# 2026-07-01 Codex-phase2 #6: secondary-MR (Egger/WM/Q) ERRORS are audited (not silently dropped). n_snps<3
# legitimately has no Egger/WM and is NOT recorded here; only thrown computation errors are.
if (nrow(.SEC_MR_FAIL) > 0) {
  cat(sprintf("NOTE: %d secondary-MR (Egger/WM/Q) computation error(s) -> the secondary-failures CSV (primary IVW/Wald unaffected)\n", nrow(.SEC_MR_FAIL))) } else cat("secondary-MR (Egger/WM/Q): 0 computation errors\n")
# 2026-06-30 (editor enhancement #1): per-outcome cell-loss breakdown so differential loss (one outcome
# losing an anomalous share -> biased Phase-5 comparative network) is quantitative, not eyeballed.
.att[, .got := .k %in% .gotc]
cat("  per-outcome MR cell loss (attempted / with-result / no-result):\n")
print(.att[, .(attempted = .N, with_result = sum(.got), no_result = sum(!.got)), by = outcome_id][order(outcome_id)])
# 2026-07-01 (gate #13 finding 3): a single outcome with anomalous cell-loss = outcome-specific systematic
# failure that the overall bound can miss -> per-outcome bound (disease GWAS cover ~82-91% of instruments,
# so normal per-outcome loss is well under 50%).
# HOW MUCH OF THE INSTRUMENT SET EACH OUTCOME ACTUALLY COVERED. The loss rates below are computed
# over ATTEMPTED cells, so an outcome that overlapped almost nothing has a tiny denominator and a
# near-zero loss rate -- it passes as a clean outcome. The eight outcomes are queried with the
# same instrument set, so they are comparable to each other, and that comparison needs no
# threshold invented for it: the same outlier rule already applied to loss below is applied to
# coverage here.
.cov <- .att[, .(n_attempted = uniqueN(gene)), by = outcome_id]
cat("  per-outcome instrument coverage (genes with >=1 instrument found):\n"); print(.cov[order(n_attempted)])
# The abort that stood here compared each outcome with half the median of the others. Removed, not
# replaced: "half the median" is an invented number, and no statistic computed within one run can
# see a shortfall that every outcome shares -- the median falls with them. It could, however, abort
# a legitimately sparse outcome. Detecting shared undercoverage needs an expected coverage per
# outcome that the study's authors supply; it is not something this script can derive. The table
# above is printed so the shortfall is visible.
.olr <- .att[, .(loss = mean(!.got)), by = outcome_id]
if (any(.olr$loss > 0.50)) stop(sprintf("01_mr_screen: outcome(s) with >50%% MR cell-loss (systematic per-outcome failure?): %s",
                                         paste(.olr[loss > 0.50, sprintf("%s=%.0f%%", outcome_id, 100*loss)], collapse = ", ")))
# 2026-07-01 (editor REQUIRE): the coarse 50% bound misses a SUBTLE per-outcome differential, which would
# bias the Phase-5 cross-disease comparability/network (one disease appearing to share fewer proteins than
# it truly does). Flag any outcome whose loss is an OUTLIER vs the median (asymmetry), not just gross loss.
.medl <- median(.olr$loss)
.outl <- .olr[loss > pmax(0.40, .medl + 0.15)]
if (nrow(.outl) > 0) stop(sprintf("01_mr_screen: per-outcome MR cell-loss ASYMMETRY (median %.1f%%; outlier(s): %s) -- a differential loss biases the cross-disease network; investigate before promoting",
                                  100*.medl, paste(.outl[, sprintf("%s=%.0f%%", outcome_id, 100*loss)], collapse = ", ")))
# 2026-07-01 (gate): catastrophic MR-cell loss = systematic harmonization/extraction failure, not normal
# variation -> automated bound (per-cell loss is variable, so a generous 60% catches gross loss only).
.loss_rate <- nrow(.att[!.k %in% .gotc]) / nrow(.att)
# The arguments used to be in the wrong order -- the artifact NAME landed in the %.0f slot, so the
# guard threw "invalid format" instead of reporting the loss rate and where to look.
if (.loss_rate > 0.60)
  stop(sprintf(paste0("01_mr_screen: MR cell-loss %.0f%% > 60%% -- systematic failure likely; ",
                      "investigate before promoting (see %s)"),
               100 * .loss_rate, mr_artifact("no_result_cells", "csv")))
# ...and NOW, with every guard passed, the full table, its three subsets and the audit CSVs that
# describe them are written together.
.secf <- file.path(REPRO, mr_artifact("secondary_failures", "csv"))
# STAMPED LIKE THE TABLES THEY DESCRIBE. These are plain CSVs, so the run identity has to be a
# column: without it a directory can hold one run's audits beside another run's MR tables and
# nothing distinguishes them.
for (.a in list(quote(.no_res_genes), quote(.no_res_cells), quote(.nd), quote(.SEC_MR_FAIL))) {
  .d <- get(as.character(.a))
  if (is.data.frame(.d) && nrow(.d)) set(.d, j = "ckm_run_id", value = CKM_RUN_ID)
}
fwrite(.no_res_genes, file.path(REPRO, mr_artifact("no_result_genes", "csv")))
fwrite(.no_res_cells, file.path(REPRO, mr_artifact("no_result_cells", "csv")))
fwrite(.nd,           file.path(REPRO, mr_artifact("no_outcome_data_cells", "csv")))
# WRITTEN EVERY RUN, empty when there were no failures. Removing it only on a clean run made the
# file's ABSENCE the signal, which meant a previous run's report had to be deleted, which needed a
# writability preflight, which covered the removal path but not the write path. An always-written
# file is always this run's, and the three audit CSVs beside it are published the same way.
fwrite(.SEC_MR_FAIL, .secf)
saveRDS(mr_all, file.path(REPRO, mr_artifact("all_results")))
prim <- mr_all[method %in% c("Inverse variance weighted","Wald ratio")]
# STAMPED LIKE THE TABLE THEY COME FROM, explicitly. These are four separate writes, so an
# interruption can leave a new all-results table beside an earlier run's subsets, and cascade/02
# checks pair MEMBERSHIP rather than identity -- a stale subset whose pairs are all still present in
# the new table is promoted and reported as current. The generation stamps are what distinguishes
# them, and MEASURED (data.table 1.17.8): a row subset `d[i]` does carry them over, while a
# column-selecting `d[, .(...)]` and rbindlist() drop them. Row subsetting is what these three lines
# use, so today they would be stamped anyway -- but that is undocumented behaviour of the subsetting
# operator, and the provenance of a published artifact should not rest on it.
.stamped <- function(d) {
  for (.a in c("ckm_run_id", "ckm_entry_script_md5", "ckm_phase1_md5", "ckm_phase1_file", "ckm_outcome_md5"))
    attr(d, .a) <- attr(mr_all, .a)
  d
}
saveRDS(.stamped(prim[bonferroni_sig==TRUE]), file.path(REPRO, mr_artifact("significant_bonferroni")))
saveRDS(.stamped(prim[fdr_sig==TRUE]),        file.path(REPRO, mr_artifact("significant_fdr")))
saveRDS(.stamped(prim[pval<SUGGESTIVE_P]),    file.path(REPRO, mr_artifact("suggestive")))

cat(sprintf("\n=== %s corrected MR ===\nN_PROTEINS=%d, Bonf=%.4g | rows=%d genes=%d | primary=%d\n", PHASE1_TAG,
            N_PROTEINS, BONF, nrow(mr_all), uniqueN(mr_all$gene), nrow(prim)))
cat(sprintf("Bonferroni=%d | FDR(by=method)=%d | suggestive=%d\n",
            sum(prim$bonferroni_sig,na.rm=TRUE), sum(prim$fdr_sig,na.rm=TRUE), prim[pval<SUGGESTIVE_P,.N]))

# ---- (D) reproduction check: canonical genes vs existing corrected MR ----
# 2026-06-30: this 1635-subset drift check is MOOT after the UKB-free source change (AF Nielsen->R10,
# CKD Pattaro->Wuttke, ncase -> 1635 subset will NOT match old-source results) AND the old canonical
# file is intentionally isolated for the clean rerun. Guard it so a clean rerun does not crash here.
.oldf <- file.path(REPRO, "mr_all_results.rds")
# Everything above is already on disk. This block only COMPARES against an older run, so a corrupt
# or schema-incompatible old file used to end the process non-zero over a complete run's outputs --
# a finished analysis reported as a failure. Its failure is reported as what it is.
if (file.exists(.oldf)) tryCatch({
  old <- readRDS(.oldf); setDT(old)
  oldp <- old[method %in% c("Inverse variance weighted","Wald ratio"), .(gene, outcome_id, method, b_old=b, p_old=pval)]
  newp <- prim[gene %in% .canon_genes, .(gene, outcome_id, method, b_new=b, p_new=pval)]
  chk <- merge(oldp, newp, by=c("gene","outcome_id","method"))
  # "1635-subset" named a set this generation does not have: with no delta file, .canon_genes is
  # every gene in the instrument set. The label has to be derived from what was compared.
  cat(sprintf("\n[repro check] pairs among the %d canonical genes matched old: %d (old primary=%d)\n",
              length(.canon_genes), nrow(chk), nrow(oldp)))
  # max() of an empty vector is -Inf, which prints as a number and reads like a measured agreement.
  if (nrow(chk))
    cat(sprintf("  max |b_new-b_old|=%.3g | max |p_new-p_old|=%.3g  (NOTE: NOT ~0 expected — AF/CKD sources changed)\n",
                max(abs(chk$b_new-chk$b_old),na.rm=TRUE), max(abs(chk$p_new-chk$p_old),na.rm=TRUE)))
  else cat("  no shared gene-outcome-method pair — nothing compared\n")
}, error = function(e)
  cat(sprintf("\n[repro check] SKIPPED — the old %s could not be compared (%s). The new artifacts above are
  unaffected: this block reads nothing they depend on.\n", .oldf, conditionMessage(e)))) else
  cat("\n[repro check] skipped — old mr_all_results.rds isolated/absent (expected on clean UKB-free rerun)\n")

if (.has_delta) {
  cat(sprintf("\n[delta hits] primary Bonferroni/suggestive among the %d recovered genes:\n",
              length(.delta_genes)))
  # head(), not [1:15]: data.table pads an out-of-range index with NA rows, so a run with fewer
  # than 15 delta hits printed blank pseudo-hits underneath the real ones.
  print(utils::head(prim[gene %in% .delta_genes & pval<SUGGESTIVE_P][order(pval)], 15)[, .(gene, outcome_name, n_snps, b=round(b,3), pval=signif(pval,3), bonf=bonferroni_sig, any_proxy)])
}
cat(sprintf(paste0("\nelapsed %.0fs. Wrote *_%s.rds and the audit CSVs to %s;\n",
                   "  outcome extracts to %s\n",
                   "  LD-proxy artifacts to %s\n"),
            proc.time()["elapsed"]-t_start, PHASE1_TAG, REPRO, LOCAL_EXTRACT_DIR,
            dirname(LD_PROXY_PREFIX)))
