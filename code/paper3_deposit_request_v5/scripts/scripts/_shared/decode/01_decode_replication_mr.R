# deCODE Replication MR
# ============================================================================
# Independent replication of the strong coloc hits (dual rule) using deCODE pQTL
# (Iceland, N≈35,287). Same SNP instruments from UKB Phase 1, but with
# independent effect estimates from deCODE.
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/scripts/_shared/decode/01_decode_replication_mr.R"

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
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
})

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))

# ── Local mode (offline): use LOCAL outcome extraction instead of OpenGWAS API ──
# Set DECODE_LOCAL <- TRUE (e.g. via localize/11 wrapper) to read outcomes from
# local sumstats; otherwise falls back to extract_outcome_data() (API).
if (!exists("DECODE_LOCAL")) DECODE_LOCAL <- FALSE
if (isTRUE(DECODE_LOCAL)) source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_local_outcome_extract.R"))

if (!exists("DECODE_OUT_DIR")) DECODE_OUT_DIR <- file.path(RESULTS_DIR, "decode_replication")
OUT_DIR   <- DECODE_OUT_DIR
CACHE_DIR <- file.path(OUT_DIR, "cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== deCODE Replication MR ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── 1. Load inputs ──────────────────────────────────────────────────────────

# Phase 1 instruments (override via DECODE_INSTR)
if (!exists("DECODE_INSTR")) DECODE_INSTR <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")
instruments <- readRDS(DECODE_INSTR)
cat(sprintf("Instruments: %d SNPs, %d genes\n", nrow(instruments), uniqueN(instruments$gene)))

# Phase 3 strong hits (override via DECODE_COLOC)
if (!exists("DECODE_COLOC")) DECODE_COLOC <- file.path(PHASE3_DIR, "coloc_results_final.rds")
coloc_all <- readRDS(DECODE_COLOC)
# THE STRONG SET MUST COME FROM A COMBINED FILE THAT MATCHES THE STAGE ARTIFACTS ON DISK
# (objective review round 2, finding 6). `12c` and `28` check this; deCODE did not, so a skipped or
# interrupted 06_combine_coloc would give a complete-looking summary over the PREVIOUS strong set,
# and the row-count assertion at the end compares only against that same stale input.
local({
  # The stage artifacts are located from the GENERATION, not from whatever directory the coloc
  # input happens to sit in -- deriving them from dirname(DECODE_COLOC) means a redirected input
  # takes its provenance check with it, which is the opposite of a check. 12c and 28 use the
  # generation constants; so does this.
  .rc <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
  .stage <- c(coloc_inputs_key_noneGFR = file.path(.rc, "coloc_results_v7_stage1.rds"),
              coloc_inputs_key_eGFR    = file.path(.rc, "coloc_egfr_results.rds"))
  for (.a in names(.stage)) {
    .v <- attr(coloc_all, .a, exact = TRUE)
    if (is.null(.v) || !nzchar(as.character(.v)[1]))
      stop(sprintf("deCODE: %s carries no %s -- it predates the stamping or came from another generation. Re-run cascade/06_combine_coloc.R.", basename(DECODE_COLOC), .a))
    if (!file.exists(.stage[[.a]]))
      stop(sprintf("deCODE: cannot verify %s -- the stage artifact %s is absent.", .a, basename(.stage[[.a]])))
    .sv <- attr(readRDS(.stage[[.a]]), "coloc_inputs_key", exact = TRUE)
    if (is.null(.sv) || !identical(as.character(.sv)[1], as.character(.v)[1]))
      stop(sprintf("deCODE: %s does not match %s -- the combined file was not built from the stage artifact now on disk. Re-run cascade/06_combine_coloc.R.", .a, basename(.stage[[.a]])))
  }
  cat("deCODE: coloc provenance verified against both stage artifacts\n")
})
strong_hits <- coloc_all[is_coloc_pass(PP.H4, PP.H3), .(gene, outcome_id, outcome_name, PP.H4)]   # dual coloc rule (D1 2026-07-04; was H4-only -> 53)
cat(sprintf("Strong hits: %d pairs, %d genes\n", nrow(strong_hits), uniqueN(strong_hits$gene)))

# Phase 2 UKB MR results (for direction comparison; override via DECODE_UKBMR)
if (!exists("DECODE_UKBMR")) DECODE_UKBMR <- file.path(PHASE2_DIR, "mr_all_results.rds")
ukb_mr <- readRDS(DECODE_UKBMR)
# Keep primary method only (Wald ratio or IVW)
ukb_mr_primary <- ukb_mr[method %in% c("Wald ratio", "Inverse variance weighted")]

# ── 2. Map deCODE files ─────────────────────────────────────────────────────

decode_pqtl_dir <- file.path(DECODE_DIR, "pqtl")
decode_files <- list.files(decode_pqtl_dir, pattern = "[.]txt[.]gz$", full.names = TRUE)

# ── ONE RELEASE, and files from any other are excluded by name rather than mis-parsed ────────────
# The deCODE download folder carries three families: the 2021 Ferkingstad release
# ({SomaID}_{ver}_{GENE}_{protein}.txt.gz), a later re-release in TWO equal variants
# (Proteomics_PC0_* and Proteomics_SMP_PC0_*, 5,284 files each), and GBR_UKB*, which is UK Biobank
# OLINK -- this study's DISCOVERY platform, so reading it here would make the replication circular.
# The manuscript cites the 2021 release, so this analysis reads only that one and the replication
# set rests on a single, describable release.
#
# The previous parser accepted the re-release too, and its else branch caught anything it did not
# recognise: Proteomics_SMP_PC0_10000_28_CRYBB2_CRBB2 has parts[2] == "SMP", so it fell through and
# took gene = parts[3] = "PC0" -- a silent gene misassignment, over 5,284 files in the release.
# Excluding by name makes that branch unreachable instead of repairing it in place.
#
# It also set source = "deCODE_2023_SomaScan", a year inferred from the filename and never checked
# against the release documentation; the field was computed and then written to no output at all.
# Removed rather than propagated. Provenance that IS kept: filename and soma_id, per row.
.other_release <- grepl("^(Proteomics_|GBR_UKB|00README)", basename(decode_files))
if (any(.other_release)) {
  cat(sprintf("deCODE: %d file(s) EXCLUDED as not the 2021 release (analysis uses one release):
",
              sum(.other_release)))
  cat(paste0("    ", basename(decode_files)[.other_release], collapse = "
"), "
")
  decode_files <- decode_files[!.other_release]
}
# Parse gene and SomaID from the ONE accepted form: {SomaID}_{ver}_{GENE}_{protein}.txt.gz
parse_decode_filename <- function(f) {
  bn <- basename(f)
  parts <- strsplit(gsub("[.]txt[.]gz$", "", bn), "_")[[1]]
  if (length(parts) < 3)
    stop("deCODE: cannot parse ", bn, " as {SomaID}_{ver}_{GENE}_{protein}")
  data.table(file = f, soma_id = parts[1], gene = parts[3], filename = bn)
}
stopifnot(identical(parse_decode_filename("x/10000_28_CRYBB2_CRBB2.txt.gz")$gene, "CRYBB2"))
decode_map <- rbindlist(lapply(decode_files, parse_decode_filename))
cat(sprintf("deCODE files: %d files, %d unique genes (2021 Ferkingstad release only)
",
            nrow(decode_map), uniqueN(decode_map$gene)))

# CLEAR THE PREVIOUS RUN'S OUTPUTS NOW, not at the saves 400 lines below. Both artifacts were only
# ever replaced by the final fwrite/saveRDS, so any failure inside the per-pair loop -- an
# unreadable gz, a harmonisation error, memory -- left the PREVIOUS summary in place, and
# cascade/13 checks only that the file exists. A partial or abandoned deCODE run would therefore be
# consumed as a complete one, and pairs that this run would have found replicated would silently
# read as absent. Removing them here makes an abandoned run fail the next step loudly instead.
# MOVED ABOVE THE COVERAGE GATE (objective review 2026-08-06, finding 11). It used to sit after
# it, on the rationale that 'a run that was never going to produce anything stops above and leaves
# the previous results alone'. That rationale is wrong: the previous results belong to a previous
# run, and a hard abort here leaves a complete-looking summary that cascade/13 accepts on
# existence alone. Aborting must never leave a readable artifact that is not this run's.
# testable) stops above and leaves the previous results alone.
local({
  .stale <- c(file.path(OUT_DIR, "decode_replication_summary.csv"),
              file.path(OUT_DIR, "decode_mr_results.rds"))
  .present <- .stale[file.exists(.stale)]
  if (length(.present)) {
    unlink(.present, force = TRUE)
    # unlink() can no-op on a locked Windows file and returns 0 either way, so the removal is
    # verified rather than assumed -- an unremovable stale file is the exact state this prevents.
    .left <- .present[file.exists(.present)]
    if (length(.left))
      stop(sprintf(paste0("deCODE: could not remove the previous run's output(s) (%s). They would ",
                          "survive a failure in this run and be read as this run's results."),
                   paste(basename(.left), collapse = ", ")))
    cat(sprintf("deCODE: cleared %d previous output(s) before the run\n", length(.present)))
  }
})

# ── THE MIRROR MUST BE COMPLETE FOR THIS STRONG SET, judged against the RELEASE MANIFEST ────────
# "Which genes have a file on disk" is not the same question as "which proteins deCODE measured",
# and answering the second with the first is the defect this whole layer was found to have: a gene
# whose file has simply not been downloaded is indistinguishable from one deCODE never assayed, and
# both end up labelled 'no deCODE pQTL data' in the manuscript. Measured 2026-08-03 on the rebuilt
# strong set: the release holds 49 assays over 41 of the 61 strong genes, the mirror held 19, and
# running this script in that state would have reported 55 pairs as NO_DECODE_DATA when 29 of them
# were merely absent from disk.
#
# So the manifest decides, and a mirror that does not cover the release's offering stops the run.
# Size equality is the check here rather than a checksum: the downloader verifies each object's MD5
# against the manifest ETag as it writes it, and re-hashing ~46 GB on every run would not.
.mani_f <- file.path(DECODE_DIR, "folder_manifest.json")
.mani <- jsonlite::fromJSON(.mani_f)$files
# "IS THIS PROTEIN IN A deCODE RELEASE" AND "IS IT IN THE ONE WE ANALYSE" ARE DIFFERENT QUESTIONS,
# AND THE STATES ANSWER THE SECOND. The manifest is read for every family so both questions can be
# answered, but the coverage STATES are computed against the 2021 family alone, because that is the
# only release this analysis reads. A protein deCODE published only in the re-release is therefore
# NOT_IN_RELEASE -- true of the release named -- and the per-gene counts keep the other question
# visible: n_assays_any_release says deCODE published something, n_assays_in_2021 says whether this
# analysis could read it. Naming both after "release" was what let one row assert both at once.
# Symbol tokens sit at different positions per family -- taking token 3 for a Proteomics_* file is
# how 5,284 files once got gene = "PC0" -- so each family is parsed by its own index, mirroring
# 00_coverage_audit.R:105-111. GBR_UKB* is excluded: it is UK Biobank Olink, this study's discovery
# platform, and counting it as deCODE coverage would be circular.
.RE2021 <- "^[0-9]+_[0-9]+_[^_]+_.+[.]txt[.]gz$"
.mani <- data.table::as.data.table(.mani)[, .(Key, Size, .bn = basename(Key))]
.mani[, family := data.table::fifelse(grepl(.RE2021, .bn), "ferkingstad_2021",
              data.table::fifelse(grepl("^Proteomics_SMP_PC0_", .bn), "rerelease_smp",
              data.table::fifelse(grepl("^Proteomics_PC0_", .bn), "rerelease_pc0", NA_character_)))]
.mani <- .mani[!is.na(family)]
.tok <- function(f, i) vapply(strsplit(sub("[.]txt[.]gz$", "", f), "_", fixed = TRUE), `[`, character(1), i)
.mani[family == "ferkingstad_2021", `:=`(gene = .tok(.bn, 3L), seqid = paste(.tok(.bn, 1L), .tok(.bn, 2L), sep = "-"))]
.mani[family == "rerelease_pc0",    `:=`(gene = .tok(.bn, 5L), seqid = paste(.tok(.bn, 3L), .tok(.bn, 4L), sep = "-"))]
.mani[family == "rerelease_smp",    `:=`(gene = .tok(.bn, 6L), seqid = paste(.tok(.bn, 4L), .tok(.bn, 5L), sep = "-"))]
.mani <- as.data.frame(.mani[!is.na(gene) & nzchar(gene), .(Key, Size, gene, family, seqid)])

strong_genes <- unique(strong_hits$gene)

# THE KEY IS UniProt, NOT THE FILENAME'S SYMBOL TOKEN (objective review 2026-08-06, finding 1).
# `.mani$gene` above is the third `_`-separated token of the filename, and matching it to
# `strong_genes` by string equality is the exact failure `_shared/decode/00_coverage_audit.R` exists
# to expose: deCODE names NECTIN2's file for its legacy symbol PVRL2, so the symbol route puts it in
# NOT_IN_RELEASE by a key failure of ours rather than by the release. (CD40 and RTBDN also sit in
# NOT_IN_RELEASE, but correctly: they exist only in the re-release, which this analysis does not
# read.) Measured effect of switching the key on this generation at the time: NOT_IN_RELEASE /
# NOT_RETRIEVED went from 28/29 to 25/32; the release-basis change below moved it again.
#
# Resolution mirrors 00_coverage_audit.R exactly (same file, same normalisation, symbol -> accession
# SET, membership tested from the instrument table's accession) so the two cannot drift: joining the
# other way multiplies rows and lets a TrEMBL entry decide a label.
# THE DECOMPRESSION'S EXIT STATUS IS OBSERVED (round 5, finding 1). `fread(cmd = "gzip -dc ...")`
# cannot see the producer, so a read error or truncation part-way through the idmapping delivers a
# SHORT table: `.res` stays non-empty so the run continues, and every gene whose mapping lay past
# the failure comes out NOT_IN_RELEASE with assay counts 0/0/0/0 -- indistinguishable from deCODE
# never having assayed it. Same pattern as the eQTL engine's reads, inline because these are the
# only two call sites and a shared copy in a third place is what drifts.
.idm_gz  <- file.path(BASE_DIR, "data/annotation/HUMAN_9606_idmapping.dat.gz")
.idm_out <- tempfile(fileext = ".dat")
# NO SHELL IN THE PATH OF THIS CHECK. This ran `bash -c 'gzip -dc ... ; echo $? > status'` and read
# the status file back. On 2026-09-09 a WSL launcher appeared at <C>/Windows/System32/bash.exe, and
# Windows resolves System32 ahead of PATH, so the bare name `bash` reached an uninstalled WSL:
# nothing ran, no status file was written, and the guard below aborted (measured 2026-09-16).
# system2() returns the child's exit status when stdout is a file, so the status is still observed.
.idm_bin <- Sys.which("gzip")
if (!nzchar(.idm_bin))
  stop("deCODE: no gzip executable on PATH -- cannot decompress the idmapping, and a partial mapping would label genes NOT_IN_RELEASE as though deCODE had never assayed them")
.idm_rc <- suppressWarnings(system2(.idm_bin,
  c("-dc", shQuote(normalizePath(.idm_gz, winslash = "/", mustWork = TRUE))),
  stdout = .idm_out, stderr = FALSE))
.idm_rc <- suppressWarnings(as.integer(.idm_rc[1]))
if (!identical(.idm_rc, 0L))
  stop(sprintf("deCODE: gzip exited %s while decompressing %s -- any rows read would be a PARTIAL mapping, and every gene past the failure would come out NOT_IN_RELEASE as though deCODE had never assayed it",
               .idm_rc, basename(.idm_gz)))
.idm <- data.table::fread(file = .idm_out, header = FALSE, sep = "\t",
                          col.names = c("acc", "type", "value"), quote = "")
unlink(.idm_out)
.idm <- .idm[type %chin% c("Gene_Name", "Gene_Synonym")]
.idm[, acc_base := sub("-.*$", "", acc)]
.sym2acc <- unique(.idm[, .(symbol = toupper(value), acc = acc_base)])
.gene_acc <- unique(data.table::as.data.table(
  data.table::fread(file = file.path(PHASE1_DIR, "ukb_cis_instruments.csv")))[gene %chin% strong_genes & nzchar(uniprot), .(gene, uniprot)])
.mdt <- data.table::as.data.table(.mani)[, .(Key, Size, seqid, symbol_uc = toupper(gene))]
.res <- unique(merge(merge(.mdt, .sym2acc, by.x = "symbol_uc", by.y = "symbol", allow.cartesian = TRUE),
                     .gene_acc, by.x = "acc", by.y = "uniprot", allow.cartesian = TRUE)[, .(gene, Key, Size, seqid)])
if (!nrow(.res))
  stop("deCODE: the UniProt route resolved no manifest file for any strong gene -- refusing to fall back to symbol equality, which is the bug this replaces.")
# gene identity now comes from the UniProt resolution, not from the filename token
.res  <- merge(.res, unique(data.table::as.data.table(.mani)[, .(Key, family)]), by = "Key")
# `.req` -- what this ANALYSIS may read -- is the 2021 family only, unchanged.
# `.genes_in_rel` below is taken from `.res`, i.e. every family. It is NOT a state basis and NOT the
# display basis; both use `.genes_2021`, so every reported decomposition is of the analysed release.
.req  <- as.data.frame(.res[family == "ferkingstad_2021", .(Key, Size, gene)])
if (!nrow(.req))
  stop("deCODE: the UniProt route resolved no 2021-release file for any strong gene.")
.want <- setNames(.req$Size, basename(.req$Key))
.got  <- file.size(file.path(decode_pqtl_dir, names(.want)))
.bad  <- names(.want)[is.na(.got) | .got != .want]

# THREE STATES, NOT TWO, AND THE RUN NO LONGER STOPS ON THE THIRD (2026-08-05, PI decision).
# Stopping was right while the alternative was to mislabel; it is wrong now that the two absences can
# be told apart in the output. download.decode.is terminates transfers after ~30 s and does not
# honour Range, so "retrieve the rest, then run" has no timeline. What must never happen is the two
# absences sharing one label: 'deCODE never assayed this protein' is a checkable claim about a third
# party's release, and reporting it for a protein deCODE did measure is false. So:
#
#   NOT_IN_RELEASE  -- the ANALYSED release has no assay for this gene. Permanent; no download helps.
#   NOT_RETRIEVED   -- the analysed release has it, this mirror does not (absent or wrong size).
#                      Dated, ours, and fixable by downloading.
#   testable        -- file present at the manifest's size; the criterion applies.
#
# "THE RELEASE" IS THE ONE THIS ANALYSIS READS. Both labels used to be computed against ANY deCODE
# release, on the reasoning that it made "deCODE never assayed this protein" literally true. But this
# analysis reads only the 2021 Ferkingstad release (see the exclusion above: the manuscript cites
# that release and the replication set must rest on a single describable one), so a protein assayed
# only in the later Proteomics_PC0 re-release was NOT assayed in the release being analysed -- and
# calling it NOT_RETRIEVED told the reader we had merely failed to fetch it, i.e. that another
# download would fix it. It would not: there is nothing in the analysed release to fetch. On this
# generation that mislabelled CD40 and RTBDN, both of which exist only in the re-release. Computing
# both states against the 2021 family makes each label true of the release it names.
#
# Counts are printed every run so a shrinking mirror is visible rather than inferred.
.ok_files      <- setdiff(names(.want), .bad)
.genes_in_rel  <- unique(.res$gene)     # ANY deCODE release; not a state basis
.genes_2021    <- unique(.req$gene)     # the ANALYSED release; the basis for both absence states
.genes_ok      <- unique(.res$gene[basename(.res$Key) %chin% .ok_files])   # resolved gene, not the filename token
DECODE_NOT_IN_RELEASE <- sort(setdiff(strong_genes, .genes_2021))
DECODE_NOT_RETRIEVED  <- sort(setdiff(.genes_2021,  .genes_ok))

cat(sprintf(paste0("deCODE coverage of this strong set: %d gene(s) testable / %d in the ANALYSED ",
                   "release but NOT RETRIEVED / %d NOT IN THE ANALYSED RELEASE  (of %d strong gene(s))\n"),
            length(.genes_ok), length(DECODE_NOT_RETRIEVED), length(DECODE_NOT_IN_RELEASE),
            length(strong_genes)))
if (length(DECODE_NOT_RETRIEVED))
  cat(sprintf("  not retrieved (%d): %s\n", length(DECODE_NOT_RETRIEVED),
              paste(DECODE_NOT_RETRIEVED, collapse = ", ")))
if (length(DECODE_NOT_IN_RELEASE))
  cat(sprintf("  not in release (%d): %s\n", length(DECODE_NOT_IN_RELEASE),
              paste(DECODE_NOT_IN_RELEASE, collapse = ", ")))
# A run with nothing testable produces a table in which the criterion never applied. That is not a
# result and must not be written as one.
if (!length(.genes_ok))
  stop(sprintf(paste0("deCODE: not one assay file for this strong set is present at its manifest ",
                      "size in %s, so every pair would be an absence. Nothing to replicate against."),
               decode_pqtl_dir))

# ONLY the files the gate validated. decode_map was built from everything in the directory that
# parses, so a file that is not a manifest entry for a strong gene -- an unfamiliar release family, a
# renamed copy -- became an extra assay for that gene and entered the multi-assay mean, the minimum
# p and the replication status. Narrowing to the validated files removes that path rather than
# testing for it, and also drops the wrong-size ones, which `.req` alone would have kept.
# HOW MANY ASSAYS THIS GENE SHOULD HAVE HAD, NOT JUST HOW MANY IT HAS (round 3, finding 4).
# Validation is per FILE but coverage is judged per GENE, so a gene keeps "testable" status while
# one of its assays is absent or wrong-sized -- and the summary then shows n_assays = 1 with no
# record that a second existed. On this generation that is not hypothetical: FLT1 has two assays,
# and losing the weak one (F = 2.36, delta p = 0.1468) would leave the strong one (F = 10.25,
# delta p = 0.0096) and flip the pair from non-evaluable to high-confidence, with an output shape
# identical to a legitimate single-assay evaluation.
# THE COUNTS AND THE STATUS MUST BE OVER THE SAME RELEASE SET (round 4, finding 1) -- and they are
# again, from the other side. That round found CD40~Stroke and RTBDN~T2D coming out NOT_RETRIEVED
# beside n_assays_in_release = 0, and resolved it by counting every family so the count matched the
# then any-release status. The states are now computed against the 2021 family, which put the same
# two pairs at NOT_IN_RELEASE beside a count of 1 -- the same contradiction mirrored. The count that
# sits beside the status is now n_assays_in_2021; n_assays_any_release keeps the wider fact under a
# name that cannot be read as the status's release.
.res_dt    <- data.table::as.data.table(.res)
# COUNT ASSAYS, NOT RELEASE-FILE COPIES (round 6, finding 1). deCODE ships the same SomaScan
# aptamer in more than one release family, so counting manifest ROWS made FLT1's two assays
# (16315-105 and 8231-122, each present in 2021, PC0 and SMP-PC0) read as six, PCSK9's one as three
# and CD40's one as two. The MR results were unaffected -- only the two 2021 files ever enter
# decode_map -- but a column named and documented as an assay count told the reader there were six
# independent FLT1 assays. SeqId is what identifies an assay, so it is what gets counted.
.assay_cov <- .res_dt[, .(n_assays_any_release = uniqueN(seqid)), by = gene]
.assay_2021 <- .res_dt[family == "ferkingstad_2021", .(n_assays_in_2021 = uniqueN(seqid)), by = gene]
.assay_ok  <- .res_dt[basename(Key) %chin% .ok_files, .(n_assays_usable = uniqueN(seqid)), by = gene]
.assay_cov <- merge(merge(.assay_cov, .assay_2021, by = "gene", all.x = TRUE),
                    .assay_ok, by = "gene", all.x = TRUE)
.assay_cov[is.na(n_assays_in_2021), n_assays_in_2021 := 0L]
.assay_cov[is.na(n_assays_usable), n_assays_usable := 0L]
# MISSING IS RELATIVE TO THE FAMILY THIS ANALYSIS READS, not to every family deCODE published.
# in_release - usable counted the re-release copies as "missing" -- FLT1 came out 6/2/2/4 when
# nothing of its is actually absent -- which is the same defect as the one this block just fixed:
# a count whose name does not match what it counts. A re-release assay is DELIBERATELY not read;
# only a 2021 assay that is absent or wrong-sized is missing.
.assay_cov[, n_assays_missing := n_assays_in_2021 - n_assays_usable]
if (any(.assay_cov$n_assays_missing > 0))
  cat(sprintf("deCODE: %d gene(s) have a 2021 assay that is absent or wrong-sized while another passed: %s%s",
              nrow(.assay_cov[n_assays_missing > 0 & n_assays_usable > 0]),
              paste(.assay_cov[n_assays_missing > 0 & n_assays_usable > 0,
                               sprintf("%s(%d/%d)", gene, n_assays_usable, n_assays_in_2021)], collapse = ", "),
              "\n"))

decode_map <- decode_map[filename %chin% .ok_files]

# AND ITS `gene` COLUMN IS REPLACED BY THE RESOLVED ONE. `parse_decode_filename()` takes gene from
# the filename's third token, so `has_decode` below (strong_genes %in% decode_map$gene) would still
# miss a protein filed under a legacy symbol even once the coverage gate above resolves it. The two
# must use the same identity or the gate and the loop disagree about which genes exist.
decode_map <- merge(decode_map[, .(file, soma_id, filename)],
                    unique(data.table(filename = basename(.res$Key), gene = .res$gene)),
                    by = "filename")

# Genes the RELEASE does not cover. Not "absent from disk": those were the same question until the
# gate above, and conflating them is what let 'not downloaded' be reported as 'no deCODE data'.
has_decode <- strong_genes[strong_genes %in% decode_map$gene]
# THE ANALYSED RELEASE, so this line and the states above decompose the same 63 genes. Built from
# .genes_in_rel (any release) it printed 42 testable beside only 19 missing -- 42 + 19 does not reach
# 63, because the two genes assayed only in the re-release fell out of both halves of the display.
no_decode  <- setdiff(strong_genes, .genes_2021)   # release coverage by the UniProt route, not the filename token
cat(sprintf("Testable genes: %d / %d (missing: %s)\n\n",
            length(has_decode), length(strong_genes), paste(no_decode, collapse=", ")))

# ── 3. Read deCODE + exact rsID match ───────────────────────────────────────

read_decode_for_snps <- function(decode_file, target_rsids) {
  # 2026-07-02 (OOM fix): STREAM + filter to target rsids -- do NOT fread the whole file. deCODE assay files are
  # ~1GB gzipped (AGER = 953MB); fread-ing the full file spikes R memory to several GB and the process is
  # OOM-killed (killed twice at the first big assay). Mirror the coloc extraction: a Python pass over the gzip
  # keeps only rows whose comma-separated `rsids` column contains a target rsID, then fread the tiny result.
  if (length(target_rsids) == 0) return(NULL)
  .tfr <- tempfile(fileext = ".txt"); writeLines(unique(as.character(target_rsids)), .tfr)
  .py  <- tempfile(fileext = ".py")
  writeLines(sprintf('
import gzip, sys
targets = set(l.strip() for l in open(r"%s") if l.strip())
with gzip.open(r"%s", "rt") as f:
    hdr = f.readline(); sys.stdout.write(hdr)
    cols = hdr.rstrip("\\n").split("\\t")
    ri = cols.index("rsids")
    for line in f:
        c = line.rstrip("\\n").split("\\t")
        if len(c) > ri and any(r in targets for r in c[ri].split(",")):
            sys.stdout.write(line)
', gsub("\\\\", "/", .tfr), gsub("\\\\", "/", decode_file)), .py)
  # stderr KEPT -- same edit and same reason as R/coloc.R: a discarded stderr turned a whole failed
  # run into a bare exit code (2026-08-02, exit 322 on all 1,380 pairs; six hours to establish
  # that python itself was fine and the machine had simply run out of commit).
  .pyerr <- tempfile(fileext = ".stderr")
  .out <- system2(PYTHON_BIN, args = .py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(.pyerr)
  unlink(c(.py, .tfr))
  if (!is.null(attr(.out, "status")) && attr(.out, "status") != 0)
    fail_loud_data(sprintf("decode STREAM extraction FAILED (python exit %s) for %s -- not a legit empty region. python stderr: %s", attr(.out, "status"), basename(decode_file), if (nzchar(.emsg)) .emsg else "<empty: no diagnostic from the process, which points at process creation rather than the script>"))
  if (length(.out) <= 1) return(NULL)
  dt <- fread(text = paste(.out, collapse = "\n"), sep = "\t")

  # Exact token match: split rsids by comma, check exact %in%
  dt[, rsid_list := strsplit(as.character(rsids), ",")]
  matched <- dt[sapply(rsid_list, function(x) any(x %in% target_rsids))]

  if (nrow(matched) == 0) return(NULL)

  # Extract which rsID matched (exact)
  matched[, matched_rsid := sapply(rsid_list, function(x) {
    hit <- x[x %in% target_rsids]
    if (length(hit) > 0) hit[1] else NA_character_
  })]

  matched[, rsid_list := NULL]
  matched[!is.na(matched_rsid)]
}

# ── 3b. Vendor effect-allele frequencies, keyed by the vendor's own join column ─────────────────
#
# BUILT ONCE, LAZILY, BY ONE STREAMING PASS. The annotation is 27,178,728 rows; holding all of it
# would cost more than it is worth, so the pass keeps only rows whose rsids field contains one of
# this generation's instrument rsIDs. That is a SELECTION filter; the JOIN is still on `Name`, as
# the vendor documents, so a multiallelic rsID cannot attach the wrong allele pair's frequency.
#
# ⚠ THIS INPUT HAS NEVER BEEN BYTE-VERIFIED. The size gate at :113-116 covers only the 2021 summary
# family's filename pattern, so this file has never been checked against a manifest and there is no
# manifest entry for it. That is the exact shape of the RTBDN failure (33.6 MB short on disk, the
# analysis ran anyway and produced beta = -15.408, F = 0.069). A truncated annotation here would not
# corrupt a number -- it would silently return NA for the missing tail, which reverts those variants
# to the old palindrome-dropping behaviour. Loud, per-run reporting of the miss rate below is what
# makes that visible; it is not a substitute for a manifest entry, which must be recorded as absent.
.DC_EAF <- NULL
.dc_eaf_load <- function() {
  if (!is.null(.DC_EAF)) return(invisible())
  f <- file.path(DECODE_DIR, "assocvariants.annotated.txt.gz")
  if (!file.exists(f)) fail_loud_data(sprintf("deCODE annotation absent: %s -- proteomics_readme.txt names it as the source of effectAlleleFreq, and without it every palindromic variant is discarded by harmonisation", f))
  .want <- unique(as.character(instruments$SNP))
  .tf <- tempfile(fileext = ".txt"); writeLines(.want, .tf)
  .py <- tempfile(fileext = ".py")
  writeLines(sprintf('
import gzip, io, sys
want = set(x.strip() for x in io.open(r"%s", encoding="utf-8") if x.strip())
sys.stdout.write("Name\\teffectAlleleFreq\\n")
with gzip.open(r"%s", "rt", encoding="utf-8", errors="replace") as fh:
    hdr = fh.readline().rstrip("\\n").split("\\t")
    ni, ri, fi = hdr.index("Name"), hdr.index("rsids"), hdr.index("effectAlleleFreq")
    for line in fh:
        p = line.rstrip("\\n").split("\\t")
        if any(t in want for t in p[ri].split(",")):
            sys.stdout.write(p[ni] + "\\t" + p[fi] + "\\n")
', .tf, f), .py)
  .err <- tempfile(fileext = ".stderr")
  .out <- system2(PYTHON_BIN, args = .py, stdout = TRUE, stderr = .err)
  .msg <- if (file.exists(.err)) paste(readLines(.err, warn = FALSE), collapse = " | ") else ""
  unlink(c(.tf, .py, .err))
  if (!is.null(attr(.out, "status")) && attr(.out, "status") != 0)
    fail_loud_data(sprintf("deCODE annotation stream FAILED (python exit %s): %s", attr(.out, "status"), .msg))
  m <- fread(text = paste(.out, collapse = "\n"), sep = "\t", colClasses = list(character = "Name"))
  m[, effectAlleleFreq := suppressWarnings(as.numeric(effectAlleleFreq))]
  # ONE NAME, ONE FREQUENCY. A duplicate would mean the vendor's stated unique key is not unique,
  # and picking one silently is how the wrong allele pair's frequency gets attached.
  .dup <- m[, .N, by = Name][N > 1]
  if (nrow(.dup)) fail_loud_data(sprintf("deCODE annotation: %d Name(s) appear more than once (e.g. %s) -- the documented join key is not unique, so no frequency can be attached without a choice", nrow(.dup), .dup$Name[1]))
  setkey(m, Name)
  .DC_EAF <<- m
  cat(sprintf("[decode] annotation: %d variant frequencies loaded for %d instrument rsIDs\n", nrow(m), length(.want)))
  invisible()
}
# COUNTED OVER DISTINCT VARIANTS, NOT OVER LOOKUPS. Variants are shared three ways here: one gene
# against several outcomes (PCSK9's chr1:55039974:T:G for both AF and CAD), two aptamers of one gene
# (FLT1's chr13:28446046:A:C), and two different genes (KCTD5 and PRSS27 both at
# chr16:2717893:T:G). Counting lookups inflates numerator and denominator equally: this generation's
# 20 lookups are 17 distinct variants. The quantity this exists to expose is missing rows in a
# possibly truncated annotation, and that is a property of variants, not of how often they were
# asked for.
.DC_EAF_SEEN <- character(0); .DC_EAF_GOT <- character(0)
.dc_eaf_for <- function(nm) {
  .dc_eaf_load()
  .nm <- as.character(nm)
  v <- .DC_EAF[.(.nm), effectAlleleFreq]
  .DC_EAF_SEEN <<- union(.DC_EAF_SEEN, .nm)
  .DC_EAF_GOT  <<- union(.DC_EAF_GOT,  .nm[is.finite(v)])
  v
}

# ── 4. Main MR loop ─────────────────────────────────────────────────────────

assay_results <- list()   # assay-level MR results
summary_rows  <- list()   # gene-outcome summary

for (i in seq_len(nrow(strong_hits))) {
  gene_i    <- strong_hits$gene[i]
  outcome_i <- strong_hits$outcome_id[i]
  outname_i <- strong_hits$outcome_name[i]

  cat(sprintf("[%d/%d] %s ~ %s (%s)\n", i, nrow(strong_hits), gene_i, outname_i, outcome_i))

  # ── Check deCODE availability
  if (!gene_i %in% decode_map$gene) {
    # WHICH absence. NO_DECODE_DATA said "deCODE has no data for this protein", which is a claim about
    # the release and is false for a gene we merely failed to download. The gate above has already
    # separated the two; a gene in neither set is an absence nothing explains, so it stops.
    .st <- if (gene_i %chin% DECODE_NOT_IN_RELEASE) "NOT_IN_RELEASE"
           else if (gene_i %chin% DECODE_NOT_RETRIEVED) "NOT_RETRIEVED"
           else stop(sprintf(paste0("deCODE: %s has no assay in decode_map, yet the manifest gate ",
                                    "classified it as neither absent-from-release nor not-retrieved. ",
                                    "That is an unexplained absence, not a data gap."), gene_i))
    summary_rows[[i]] <- data.table(
      gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
      status = .st, decode_b = NA_real_, decode_se = NA_real_,
      decode_b_min = NA_real_, decode_b_max = NA_real_,
      decode_pval_iut = NA_real_, delta_se = NA_real_, delta_pval_iut = NA_real_,
      decode_F_min = NA_real_,
      assay_Q = NA_real_, assay_Q_df = NA_real_, assay_Q_pval = NA_real_,
      ukb_b = NA_real_, direction_concordant = NA,
      # No assay was ever attempted for this pair, so this is NOT "did not replicate".
      weak_instrument = NA, high_confidence = FALSE,
      n_assays = 0L, n_assays_evaluable = 0L, evaluable_class = "none",
      assay_concordance = NA_character_,
      harmonise_palindromic = NA, harmonise_ambiguous = NA,
      harmonise_mr_keep = NA, harmonise_remove = NA
    )
    cat(sprintf("  -> %s\n", .st))
    next
  }

  # ── Get instrument SNPs for this gene
  gene_instr <- instruments[gene == gene_i]
  target_rsids <- gene_instr$SNP[grepl("^rs[0-9]+$", gene_instr$SNP)]

  if (length(target_rsids) == 0) {
    summary_rows[[i]] <- data.table(
      gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
      status = "NO_EXACT_RSID_MATCH", decode_b = NA_real_, decode_se = NA_real_,
      decode_b_min = NA_real_, decode_b_max = NA_real_,
      decode_pval_iut = NA_real_, delta_se = NA_real_, delta_pval_iut = NA_real_,
      decode_F_min = NA_real_,
      assay_Q = NA_real_, assay_Q_df = NA_real_, assay_Q_pval = NA_real_,
      ukb_b = NA_real_, direction_concordant = NA,
      # No assay was ever attempted for this pair, so this is NOT "did not replicate".
      weak_instrument = NA, high_confidence = FALSE,
      n_assays = 0L, n_assays_evaluable = 0L, evaluable_class = "none",
      assay_concordance = NA_character_,
      harmonise_palindromic = NA, harmonise_ambiguous = NA,
      harmonise_mr_keep = NA, harmonise_remove = NA
    )
    cat("  -> NO_EXACT_RSID_MATCH (instrument has no rsID)\n")
    next
  }

  # ── Read all deCODE assays for this gene
  gene_decode_files <- decode_map[gene == gene_i]
  assay_mr_results <- list()

  for (j in seq_len(nrow(gene_decode_files))) {
    soma_id_j  <- gene_decode_files$soma_id[j]
    file_j     <- gene_decode_files$file[j]
    filename_j <- gene_decode_files$filename[j]

    cat(sprintf("  Assay %d/%d: %s\n", j, nrow(gene_decode_files), filename_j))

    # Read & match
    matched <- read_decode_for_snps(file_j, target_rsids)
    if (is.null(matched) || nrow(matched) == 0) {
      cat("    -> No rsID match in deCODE file\n")
      assay_mr_results[[j]] <- data.table(
        gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
        soma_id = soma_id_j, filename = filename_j,
        status = "NO_EXACT_RSID_MATCH"
      )
      next
    }

    # Compute decode F-stat
    matched[, decode_F := Beta^2 / SE^2]
    cat(sprintf("    Matched %d SNP(s), deCODE F-stat: %s\n",
                nrow(matched), paste(round(matched$decode_F, 1), collapse=", ")))

    # weak_flag is decided at the fit, from the contrast the fit kept (see F_val below).
    # Deciding it here, over every row read out of the assay, let a DISCARDED duplicate
    # allele contrast set it.

    # ── Format deCODE as exposure, WITH the vendor-supplied effect-allele frequency
    #
    # THE FREQUENCY WAS AVAILABLE ALL ALONG AND WAS NEVER READ. proteomics_readme.txt:18-19 states
    # that deCODE provides assocvariants.annotated.txt.gz specifically because "the effectAllele is
    # not always the minor allele, and therefore the ImpMAF is not always the frequency of the
    # effect allele" -- i.e. the very objection recorded in the comment this replaces. The file has
    # been on disk since 2026-04-27 (441,697,192 bytes).
    #
    # WITHOUT IT, harmonise_data(action = 2) DISCARDS EVERY PALINDROMIC VARIANT, because strand for
    # an A/T or C/G variant can only be settled by frequency. That is a systematic loss over a whole
    # class of variants, not a property of any one pair: 8 of the 63 strong-set instruments are
    # palindromic, carrying 10 pairs. On this generation's retrievable files it changes exactly one
    # reported pair (PAM~T2D) because the other nine are already blocked upstream -- that narrowness
    # is a consequence of the download failure, not of the rule.
    #
    # JOINED ON `Name`, NOT ON rsID, because the vendor says so ("This file can be mapped to the
    # other files using the Name column") and because rsID is neither unique nor always present in
    # the annotation. A Name resolving to more than one annotation row is a stop, not a pick-one.
    .af <- .dc_eaf_for(matched$Name)
    exp_df <- data.frame(
      SNP            = matched$matched_rsid,
      beta           = matched$Beta,
      se             = matched$SE,
      pval           = matched$Pval,
      effect_allele  = matched$effectAllele,
      other_allele   = matched$otherAllele,
      eaf            = .af,
      samplesize     = matched$N,
      exposure       = paste0(gene_i, "_deCODE_", soma_id_j),
      stringsAsFactors = FALSE
    )

    exp_fmt <- tryCatch(
      format_data(
        dat               = exp_df,
        type              = "exposure",
        snp_col           = "SNP",
        beta_col          = "beta",
        se_col            = "se",
        effect_allele_col = "effect_allele",
        other_allele_col  = "other_allele",
        pval_col          = "pval",
        samplesize_col    = "samplesize",
        phenotype_col     = "exposure",
        # eaf_col is the annotation's effectAlleleFreq, NOT ImpMAF. ImpMAF is a minor-allele
        # frequency and would be wrong for every variant whose effect allele is the major one.
        eaf_col           = "eaf"
      ),
      error = function(e) fail_loud_data(sprintf("deCODE format_data(exposure) THREW for %s: %s -- malformed deCODE instrument, fatal", gene_i, conditionMessage(e)))  # 2026-07-01 (closure sweep): thrown = genuine failure feeding tiering, not a legit non-result
    )
    if (is.null(exp_fmt)) {
      assay_mr_results[[j]] <- data.table(
        gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
        soma_id = soma_id_j, filename = filename_j,
        status = "FORMAT_ERROR"
      )
      next
    }

    # ── Get outcome data (with cache)
    cache_file <- file.path(CACHE_DIR, sprintf("outcome_%s_%s.rds",
                                                gsub("[^a-zA-Z0-9]", "_", outcome_i),
                                                paste(target_rsids, collapse = "_")))
    # 2026-07-01 Codex-p3rest #2: in DECODE_LOCAL mode NEVER read this shared-key (outcome+rsids only) cache --
    # an older run may have cached OpenGWAS/proxy or pre-UKB-free outcome data under the same key; local
    # extraction is offline + cheap, so always re-extract. Cache is only for the (slow) API path.
    if (!isTRUE(DECODE_LOCAL) && file.exists(cache_file)) {
      out_dat <- readRDS(cache_file)
    } else {
      out_dat <- tryCatch(
        if (isTRUE(DECODE_LOCAL))
          # proxies = TRUE MATCHES THE API BRANCH BELOW AND THE DISCOVERY ARM. The local extractor
          # defaults proxies = FALSE so each caller opts in; the API branch already passes TRUE, so
          # the two branches of this same call were reading different outcome data, and the offline
          # branch is the one this project actually runs.
          # palindromic_fallthrough = TRUE for the same reason as the ARIC arm: a palindromic direct
          # row that harmonisation will discard is not a direct row, and leaving it in blocks the
          # proxy. This arm's one HARMONISE_PALINDROMIC_DROP is PAM~T2D, whose direction verdict has
          # therefore never been computed; opting in produces it for the first time.
          local_extract_outcome_data(snps = target_rsids, outcomes = outcome_i,
                                     proxies = TRUE, palindromic_fallthrough = TRUE)
        else
          extract_outcome_data(
            snps          = target_rsids,
            outcomes      = outcome_i,
            proxies       = TRUE,
            rsq           = 0.8,
            align_alleles = 1,
            palindromes   = 1,
            maf_threshold = 0.3
          ),
        error = function(e) fail_loud_data(sprintf("deCODE outcome extraction THREW for %s ~ %s: %s -- fatal (legit no-data is the 0-row OUTCOME_MISSING below)", gene_i, outcome_i, conditionMessage(e)))  # 2026-07-01 (closure sweep)
      )
      if (!isTRUE(DECODE_LOCAL) && !is.null(out_dat)) saveRDS(out_dat, cache_file)   # 2026-07-01 Codex-p3rest #2: only the API path uses this cache; don't leave a local cache a later API run could misread
    }

    if (is.null(out_dat) || nrow(out_dat) == 0) {
      assay_mr_results[[j]] <- data.table(
        gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
        soma_id = soma_id_j, filename = filename_j,
        status = "OUTCOME_MISSING"
      )
      cat("    -> OUTCOME_MISSING\n")
      next
    }

    # ── Harmonise
    # DIRECT AND PROXY ROWS SEPARATELY, proxy rows at action = 1. Identical to the ARIC arm and for
    # the identical reason: a proxy row's strand was already resolved inside the orienter by the
    # same frequency rule, under the coupling guarantee r2 >= 0.8 supplies. Re-testing it here
    # re-adjudicates a settled question, and for a palindrome the answer is "ambiguous" -- which
    # discards the very proxy that was fetched to rescue an unresolvable direct palindrome.
    # `%in% TRUE` not `== TRUE`: NA in a row subscript yields a row of NAs instead of dropping it.
    .dc_prox <- if ("proxy.outcome" %in% names(out_dat)) out_dat$proxy.outcome %in% TRUE
                else rep(FALSE, nrow(out_dat))
    .dc_harm <- function(d, act) {
      if (is.null(d) || !nrow(d)) return(NULL)
      tryCatch(
        harmonise_data(exposure_dat = exp_fmt, outcome_dat = d, action = act),
        error = function(e) fail_loud_data(sprintf("deCODE harmonise_data THREW for %s ~ %s (action=%d): %s -- fatal (legit palindromic/mr_keep drop is the HARMONISE_DROPPED below)", gene_i, outcome_i, act, conditionMessage(e)))
      )
    }
    harm <- rbind(.dc_harm(out_dat[!.dc_prox, , drop = FALSE], 2L),
                  .dc_harm(out_dat[ .dc_prox, , drop = FALSE], 1L))

    if (is.null(harm) || nrow(harm) == 0 || sum(harm$mr_keep) == 0) {
      # Check reason
      drop_reason <- "HARMONISE_DROPPED"
      if (!is.null(harm) && nrow(harm) > 0) {
        if (any(harm$palindromic & harm$ambiguous, na.rm = TRUE)) {
          drop_reason <- "HARMONISE_PALINDROMIC_DROP"
        }
      }
      assay_mr_results[[j]] <- data.table(
        gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
        soma_id = soma_id_j, filename = filename_j,
        status = drop_reason
      )
      cat(sprintf("    -> %s\n", drop_reason))
      next
    }

    # Record harmonise diagnostics
    harm_diag <- data.table(
      palindromic = harm$palindromic[1],
      ambiguous   = if ("ambiguous" %in% names(harm)) harm$ambiguous[1] else NA,
      mr_keep     = harm$mr_keep[1],
      remove      = if ("remove" %in% names(harm)) harm$remove[1] else NA
    )

    harm_keep <- harm[harm$mr_keep, ]
    n_snps <- nrow(harm_keep)

    # ── Run MR
    primary_method <- if (n_snps == 1) "mr_wald_ratio" else "mr_ivw"
    mr_res <- tryCatch(mr(harm, method_list = primary_method), error = function(e) fail_loud_data(sprintf("deCODE mr() THREW for %s ~ %s: %s -- computation failure feeding tiering, fatal", gene_i, outcome_i, conditionMessage(e))))  # 2026-07-01 (closure sweep)

    if (is.null(mr_res) || nrow(mr_res) == 0)
      # 2026-07-01 (closure sweep, Codex): mr() returning EMPTY after a SUCCESSFUL harmonise (>=1 mr_keep SNP)
      # is a genuine computation failure, not a legit non-result. It feeds tiering (decode_status) and 99
      # treats deCODE as report-only, so fail loud rather than record a silent MR_FAILED row.
      fail_loud_data(sprintf("deCODE mr() returned EMPTY/NULL after successful harmonise for %s ~ %s -- computation failure feeding tiering, fatal", gene_i, outcome_i))

    mr_dt <- as.data.table(mr_res)

    # Delta-method SE: properly incorporates exposure SE uncertainty
    # For Wald ratio: SE_delta = sqrt(SE_wald^2 + b_MR^2 / F)
    # For IVW: TwoSampleMR already uses proper multi-SNP SE; delta approximation
    #   is negligible for strong instruments (F>>10)
    b_mr  <- mr_dt$b[1]
    se_mr <- mr_dt$se[1]
    # F COMES FROM THE CONTRAST THE FIT KEPT, not from every row read out of the assay. deCODE
    # can report two allele contrasts for one rsID; format_data keeps the first instance, so a
    # min over `matched` could take the F of the row the fit DISCARDED -- measured on
    # TINAGL1 ~ CAD, rs945214: [value withheld, release rule F15] (discarded T/A) in place of [value withheld, release rule F15] (kept G/A), and that
    # propagated into delta_se and delta_pval. beta/se here are the exposure values harmonise
    # handed the estimator, so this is the same Beta^2/SE^2 as the read-level diagnostic above,
    # evaluated on the rows actually used. (2026-06-30 gate obj #12 made this a MIN rather than
    # [1]; the MIN is kept, its domain is what changes.)
    F_val     <- min((harm_keep$beta.exposure / harm_keep$se.exposure)^2, na.rm = TRUE)
    weak_flag <- F_val < 10
    delta_se_val   <- sqrt(se_mr^2 + b_mr^2 / F_val)
    delta_pval_val <- 2 * pnorm(-abs(b_mr / delta_se_val))

    assay_mr_results[[j]] <- data.table(
      gene         = gene_i,
      outcome_id   = outcome_i,
      outcome_name = outname_i,
      soma_id      = soma_id_j,
      filename     = filename_j,
      n_snps       = n_snps,
      method       = mr_dt$method[1],
      decode_b     = b_mr,
      decode_se    = se_mr,
      decode_pval  = mr_dt$pval[1],
      delta_se     = delta_se_val,
      delta_pval   = delta_pval_val,
      decode_F_stat = F_val,
      weak_instrument = weak_flag,
      palindromic     = harm_diag$palindromic,
      ambiguous       = harm_diag$ambiguous,
      mr_keep         = harm_diag$mr_keep,
      remove          = harm_diag$remove,
      status          = "MR_OK"
    )
    cat(sprintf("    -> MR_OK: b=%.4f, p=%.2e (delta_p=%.2e), F=%.1f%s\n",
                b_mr, mr_dt$pval[1], delta_pval_val, F_val,
                if (weak_flag) " [WEAK]" else ""))
  }

  # ── Combine assay results for this gene-outcome pair
  assay_dt <- rbindlist(assay_mr_results, fill = TRUE)

  # Add to global assay-level results
  assay_results[[i]] <- assay_dt

  # ── Build gene-outcome summary row
  ok_assays <- assay_dt[status == "MR_OK"]

  # ══════════════════════════════════════════════════════════════════════════════════════════════
  # ASSAY BOOKKEEPING (statistical + adversarial editor, 2026-08-07, concurring)
  #
  # `n_assays` COUNTS ASSAYS ATTEMPTED, NOT ASSAYS THAT WORKED. It used to be nrow(ok_assays) in the
  # success branches and nrow(assay_dt) in the failure branch, so a pair whose second aptamer failed
  # to harmonise reported n_assays = 1 and became indistinguishable from a genuinely single-aptamer
  # pair -- exactly the collapse that lets a partially-evaluated pair be read as fully evaluated.
  # The two counts are now separate columns and both are emitted on every row.
  #
  # THE BOOLEAN `evaluable` IS GONE. It could not distinguish "no assay cleared the gate" from
  # "some did", and every downstream gate that read it silently treated partial evidence as
  # complete. `evaluable_class` is three-valued and every gate below tests == "all"; a FALSE that
  # arises from evaluable_class != "all" means NOT EVALUABLE, and must never be reported as
  # "did not replicate".
  # ⚠️ UNSIGNED, DELIBERATELY ABSENT: the multi-assay FLAG (conflict 2) is the one item the two
  # editors have not converged on -- adversarial wants `multi_assay_estimate = (n_assays > 1)`,
  # statistical wants `assay_discordant = (sign disagreement) OR (evaluable_class == "some")`.
  # Neither is written here. Adding whichever I preferred is exactly the act that put an invented
  # rule ("least favourable assay = largest delta p") through seven review rounds and into a
  # production run before both editors vetoed it. It goes in once, after they sign.
  n_att <- nrow(assay_dt)
  n_ev  <- nrow(ok_assays)
  ev_class <- if (n_ev == 0L) "none" else if (n_ev == n_att) "all" else "some"

  # COCHRAN'S Q ACROSS EVALUABLE ASSAYS. Reported as a LOWER BOUND on heterogeneity, not an
  # unbiased test: the assays are estimated against the SAME outcome GWAS, so their sampling errors
  # are positively correlated and the usual independence assumption understates the true spread.
  # A large Q is therefore informative; a small Q is not evidence of agreement.
  # WEIGHTED BY THE DELTA SE, NOT THE FIRST-ORDER MR SE (review finding 4, 2026-08-07). The
  # first-order SE treats the exposure estimate as known, so exposure uncertainty was charged to
  # between-assay heterogeneity. This pipeline computes a second-order (delta) SE precisely because
  # that assumption is false. Measured on FLT1~CAD, whose weaker assay has F = 2.36: first-order
  # gives Q = 4.624, df = 1, p = 0.0315; the delta SE gives Q = 0.569, p = 0.451. The p-value is
  # UNCALIBRATED either way -- the assays share one outcome GWAS, so their errors are positively
  # correlated and Q is a LOWER BOUND on heterogeneity, never a test.
  .cochran_q <- function(b, se) {
    ok <- is.finite(b) & is.finite(se) & se > 0
    if (sum(ok) < 2L) return(c(Q = NA_real_, df = NA_real_, p = NA_real_))
    w <- 1 / se[ok]^2; bb <- sum(w * b[ok]) / sum(w)
    Q <- sum(w * (b[ok] - bb)^2); df <- sum(ok) - 1L
    c(Q = Q, df = df, p = stats::pchisq(Q, df, lower.tail = FALSE))
  }

  if (nrow(ok_assays) == 0) {
    # Use the most informative failure status
    fail_status <- assay_dt$status[1]
    summary_rows[[i]] <- data.table(
      gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
      status = fail_status, decode_b = NA_real_, decode_se = NA_real_,
      decode_b_min = NA_real_, decode_b_max = NA_real_,
      decode_pval_iut = NA_real_, delta_se = NA_real_, delta_pval_iut = NA_real_,
      decode_F_min = NA_real_,
      assay_Q = NA_real_, assay_Q_df = NA_real_, assay_Q_pval = NA_real_,
      ukb_b = NA_real_, direction_concordant = NA,
      # NOT "did not replicate": nothing was evaluable, so nothing was tested. high_confidence is
      # FALSE rather than NA because every downstream gate treats it as a boolean, and an NA there
      # propagates into gates that then read as "unknown, therefore pass" in `isTRUE`-free code.
      weak_instrument = NA, high_confidence = FALSE,
      n_assays = n_att, n_assays_evaluable = 0L, evaluable_class = "none",
      assay_concordance = NA_character_,
      harmonise_palindromic = NA, harmonise_ambiguous = NA,
      harmonise_mr_keep = NA, harmonise_remove = NA
    )
    next
  }

  # Get UKB MR direction
  ukb_row <- ukb_mr_primary[gene == gene_i & outcome_id == outcome_i]
  # 2026-07-01 (editor REQUIRE): a deCODE-tested (strong) pair should BY CONSTRUCTION have a UKB primary MR
  # row; a missing one is a cross-file JOIN INCONSISTENCY (not deCODE partiality), and silently recording
  # NO_UKB_COMPARISON would demote it to Moderate. Fail loud so the inconsistency is investigated.
  if (nrow(ukb_row) == 0)
    fail_loud_data(sprintf("deCODE %s ~ %s: NO matching UKB primary MR row for a tested pair -- cross-file join inconsistency", gene_i, outcome_i))
  ukb_b_val <- ukb_row$b[1]

  # For single-assay genes: straightforward
  if (nrow(ok_assays) == 1) {
    dir_conc <- if (!is.na(ukb_b_val)) sign(ok_assays$decode_b[1]) == sign(ukb_b_val) else NA

    final_status <- if (is.na(dir_conc)) {
      "NO_UKB_COMPARISON"
    } else if (ok_assays$decode_pval[1] < 0.05 && dir_conc) {
      "REPLICATED"
    } else if (dir_conc) {
      "DIRECTION_ONLY"
    } else {
      "NOT_REPLICATED"
    }

    # FORCED FALSE WHEN NOT EVERY ATTEMPTED ASSAY IS EVALUABLE. One evaluable assay out of two is
    # partial evidence; calling it high-confidence would let an aptamer that could not be tested
    # count as one that agreed.
    hi_conf <- ev_class == "all" && !is.na(dir_conc) && dir_conc &&
               ok_assays$decode_F_stat[1] >= 10 &&
               ok_assays$delta_pval[1] < 0.05

    summary_rows[[i]] <- data.table(
      gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
      status = final_status,
      decode_b = ok_assays$decode_b[1], decode_se = ok_assays$decode_se[1],
      decode_b_min = ok_assays$decode_b[1], decode_b_max = ok_assays$decode_b[1],
      decode_pval_iut = ok_assays$decode_pval[1],
      delta_se = ok_assays$delta_se[1], delta_pval_iut = ok_assays$delta_pval[1],
      decode_F_min = ok_assays$decode_F_stat[1],
      assay_Q = NA_real_, assay_Q_df = NA_real_, assay_Q_pval = NA_real_,
      ukb_b = ukb_b_val, direction_concordant = dir_conc,
      weak_instrument = ok_assays$weak_instrument[1],
      high_confidence = hi_conf,
      n_assays = n_att, n_assays_evaluable = n_ev, evaluable_class = ev_class,
      assay_concordance = NA_character_,
      harmonise_palindromic = ok_assays$palindromic[1],
      harmonise_ambiguous = ok_assays$ambiguous[1],
      harmonise_mr_keep = ok_assays$mr_keep[1],
      harmonise_remove = ok_assays$remove[1]
    )
  } else {
    # Multi-assay genes (NOTCH2, FLT1): both assays, check concordance
    directions <- sign(ok_assays$decode_b)
    assay_conc <- if (length(unique(directions)) == 1) "ROBUST" else "ASSAY_DISCORDANT"

    # Use mean of assay betas for summary direction check
    # THE REPORTED EFFECT IS AN ASSAY'S, NOT AN AVERAGE OF ASSAYS (objective review 2026-08-06,
    # finding 3). The arithmetic mean emits a number belonging to no assay: for FLT1~CAD the two
    # betas are 2.8649 and 1.3261 with SEs 0.6494 and 0.3006, and the mean pair 2.0955 / 0.4750 is
    # neither, while downstream forest intervals are drawn from these fields. Under the
    # intersection-union rule every statistic takes its LEAST FAVOURABLE value, which for an effect
    # size is the smallest |beta| and for its uncertainty the largest SE -- the same convention the
    # ARIC arm uses, so the two arms remain comparable.
    # NO SINGLE ESTIMATE IS EMITTED FOR A MULTI-ASSAY PAIR (statistical + adversarial editor,
    # 2026-08-07, both vetoing the previous rule).
    #
    # Round 1 took min|beta| and max SE INDEPENDENTLY, producing a beta from one assay beside an SE
    # from another -- a pair belonging to no estimator, hence no sampling distribution and no valid
    # interval, while reporting/decode_forest.R:34 draws decode_b +/- 1.96*decode_se from exactly
    # these fields. Round 2 replaced that with "the least favourable assay, chosen by largest delta
    # p, supplies beta and SE together". That was MY OWN INVENTION and both editors rejected it:
    # delta p and the reported beta/se are monotone functions of the same z, so selecting the assay
    # BY its delta p and then reporting that assay's interval is result-dependent selection -- the
    # interval's coverage is destroyed, and it is destroyed in the direction that makes the pair
    # look worse, which is not conservatism but a different estimand.
    #
    # The resolution is to emit NOTHING that pretends to be "the" estimate. Multi-assay pairs carry
    # NA in decode_b / decode_se / delta_se, the observed RANGE in decode_b_min / decode_b_max, and
    # Cochran's Q for heterogeneity. Consumers that need a plottable interval must draw one line
    # PER ASSAY from decode_mr_results.rds, which persists per-assay statistics; the four consumers
    # are required to ERROR on n_assays > 1 rather than silently take row 1.
    lf_b   <- NA_real_
    lf_se  <- NA_real_
    lf_dse <- NA_real_
    b_rng  <- range(ok_assays$decode_b, na.rm = TRUE)
    .q     <- .cochran_q(ok_assays$decode_b, ok_assays$delta_se)

    # INTERSECTION-UNION ACROSS APTAMERS, matching the ARIC arm (statistical + adversarial editor,
    # 2026-08-06, concurring). This previously took min(p) across assays while taking min(F): a
    # union statistic for the effect combined with an intersection statistic for instrument
    # strength. That composite corresponds to no null hypothesis, and it lets the reported F come
    # from a different aptamer than the reported p. Under the IUT every assay must replicate, each
    # statistic takes its LEAST favourable value, and the test is level-alpha with no multiplicity
    # adjustment needed. It also matters scientifically here: disagreement between aptamers is an
    # epitope/binding signal in this project, and min(p) systematically washes that signal out.
    # min(F) was already the least favourable value and is unchanged.
    #
    # NA HANDLING IS NOT SYMMETRIC UNDER THE SWITCH. min(x, na.rm=TRUE) over all-NA returns +Inf,
    # which failed the "< 0.05" test and was therefore safe; max() returns -Inf, which would PASS
    # it. The direction of the trap reversed with the rule, so the all-NA case is made explicit.
    max_p <- if (all(is.na(ok_assays$decode_pval))) NA_real_ else max(ok_assays$decode_pval, na.rm = TRUE)
    # DIRECTION COMES FROM THE PER-ASSAY SIGNS, NEVER FROM decode_b. That field is now NA for every
    # multi-assay pair by design, so `sign(lf_b)` silently returned NA and the pair was recorded as
    # NO_UKB_COMPARISON -- i.e. "we had nothing to compare against" -- when in fact both assays had
    # estimates and both agreed. Caught by the smoke test on FLT1~CAD, which is exactly the pair the
    # minimal trigger set exists to exercise.
    .dirs <- sign(ok_assays$decode_b)
    .dirs <- .dirs[is.finite(.dirs) & .dirs != 0]
    dir_conc <- if (!is.na(ukb_b_val) && length(.dirs) && length(unique(.dirs)) == 1L)
                  .dirs[1] == sign(ukb_b_val) else NA

    final_status <- if (assay_conc == "ASSAY_DISCORDANT") {
      "ASSAY_DISCORDANT"
    } else if (is.na(dir_conc)) {
      "NO_UKB_COMPARISON"
    } else if (isTRUE(max_p < 0.05) && dir_conc) {
      "REPLICATED"
    } else if (dir_conc) {
      "DIRECTION_ONLY"
    } else {
      "NOT_REPLICATED"
    }

    # Delta for multi-assay: least favourable across assays (IUT; see the note above)
    max_delta_p <- if (all(is.na(ok_assays$delta_pval))) NA_real_ else max(ok_assays$delta_pval, na.rm = TRUE)
    min_F       <- min(ok_assays$decode_F_stat)
    # INTERSECTION ACROSS ASSAYS, AND FORCED FALSE WHEN NOT EVERY ATTEMPTED ASSAY IS EVALUABLE.
    # Direction is taken from the PER-ASSAY signs, never from decode_b -- which is now NA here, so
    # a consumer that reconstructed direction from it would silently get NA rather than a sign.
    hi_conf <- ev_class == "all" && assay_conc == "ROBUST" && !is.na(dir_conc) && dir_conc &&
               min_F >= 10 && isTRUE(max_delta_p < 0.05)

    summary_rows[[i]] <- data.table(
      gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
      status = final_status,
      # NA BY DESIGN for multi-assay: no single number represents these assays. See the note above.
      decode_b = lf_b, decode_se = lf_se,
      decode_b_min = b_rng[1], decode_b_max = b_rng[2],
      decode_pval_iut = max_p,
      delta_se = lf_dse, delta_pval_iut = max_delta_p,
      decode_F_min = min_F,
      # A LOWER BOUND ON HETEROGENEITY, not an unbiased test: the assays share one outcome GWAS, so
      # their errors are positively correlated and Q understates the true spread.
      assay_Q = unname(.q["Q"]), assay_Q_df = unname(.q["df"]), assay_Q_pval = unname(.q["p"]),
      ukb_b = ukb_b_val, direction_concordant = dir_conc,
      weak_instrument = any(ok_assays$weak_instrument),
      high_confidence = hi_conf,
      n_assays = n_att, n_assays_evaluable = n_ev, evaluable_class = ev_class,
      assay_concordance = assay_conc,
      harmonise_palindromic = ok_assays$palindromic[1],
      harmonise_ambiguous = ok_assays$ambiguous[1],
      harmonise_mr_keep = ok_assays$mr_keep[1],
      harmonise_remove = ok_assays$remove[1]
    )
  }

  cat(sprintf("  => %s\n", summary_rows[[i]]$status))
}

# ── 5. Combine & save ───────────────────────────────────────────────────────

all_assay   <- rbindlist(assay_results, fill = TRUE)
all_summary <- rbindlist(summary_rows, fill = TRUE)

# ── 6. Assertions ────────────────────────────────────────────────────────────

# INVARIANTS ON THE ASSAY BOOKKEEPING (both editors, 2026-08-07, verbatim requirement).
#
# DERIVATION CONVENTIONS, stated once so four consumers cannot each invent their own:
#   sign disagreement across assays  <=>  decode_b_min * decode_b_max < 0
#   gate-level conflict              <=>  evaluable_class == "some"
# Neither is emitted as a column: both were proposed (multi_assay_estimate, assay_discordant) and
# both editors rejected both, because each is a pure function of columns already on the row and the
# OR of the two loses the distinction between an estimate-level conflict and a gate-level one.
#
# "some" REQUIRES AT LEAST TWO ATTEMPTED ASSAYS, by the definition 0 < n_assays_evaluable <
# n_assays. A single attempted assay that fails the gate gives n_assays_evaluable = 0 = all of them,
# hence "none". If "some" ever appears on a one-assay row the derivation is broken, and it breaks in
# the direction this class exists to prevent -- partial evidence reading as complete.
.bad_some <- all_summary[evaluable_class == "some" & n_assays <= 1L]
if (nrow(.bad_some))
  fail_loud_data(sprintf("evaluable_class == 'some' on %d row(s) with n_assays <= 1: %s",
                         nrow(.bad_some), paste(.bad_some$gene, collapse = ", ")))

# decode_b_min/max ARE THE UNIVERSAL ACCESSOR and must be populated wherever an assay was
# evaluable -- degenerate to the single value when there is one. decode_b is the fail-closed field
# (NA for multi-assay by design); min/max must NOT be, or every consumer needs an NA branch.
.bad_rng <- all_summary[n_assays_evaluable >= 1L & (is.na(decode_b_min) | is.na(decode_b_max))]
if (nrow(.bad_rng))
  fail_loud_data(sprintf("decode_b_min/max NA on %d evaluable row(s): %s",
                         nrow(.bad_rng), paste(.bad_rng$gene, collapse = ", ")))

cat("\n=== Assertions ===\n")

# A1: one summary row per strong-hit pair (invariant; count is whatever the input set is)
stopifnot("Summary must match strong_hits" = nrow(all_summary) == nrow(strong_hits))
cat(sprintf("PASS: %d summary rows\n", nrow(all_summary)))

# A2/A3: testable vs the two absences — counts are DATA-DEPENDENT on the strong-hit set.
# (the strong-hit set comes from 06's combined v7 table and includes eGFR; the absences follow from
#  what the release holds and what this mirror has.) Do NOT hardcode the counts — assert only the
# structural invariant that every pair is exactly one of the three, and report.
#
# 2026-08-05: this counted `status != "NO_DECODE_DATA"` as testable. That token is no longer emitted
# for a coverage absence -- the two absences are NOT_IN_RELEASE and NOT_RETRIEVED -- so on the real
# pair set it would have reported every pair testable (74 rather than 19) and no absences at all,
# putting the replication rate over a denominator four times too large. The CSV statuses were right;
# only this summary was wrong, which is the shape that gets read and quoted rather than recomputed.
.ABSENT_STATUS <- c("NOT_IN_RELEASE", "NOT_RETRIEVED")
n_testable       <- sum(!all_summary$status %chin% .ABSENT_STATUS)
n_not_in_release <- sum(all_summary$status == "NOT_IN_RELEASE")
n_not_retrieved  <- sum(all_summary$status == "NOT_RETRIEVED")
stopifnot("testable + NOT_IN_RELEASE + NOT_RETRIEVED must equal total pairs" =
            n_testable + n_not_in_release + n_not_retrieved == nrow(all_summary))
cat(sprintf("INFO: %d NON-ABSENT pair(s); %d NOT_RETRIEVED, %d NOT_IN_RELEASE (of %d total)\n",
            n_testable, n_not_retrieved, n_not_in_release, nrow(all_summary)))
# THE MANUSCRIPT'S "testable" IS A DIFFERENT PREDICATE FROM `n_testable`, and printing the latter
# under that word invited 49 to be quoted as the denominator of the 24 high-confidence replications
# (p3_r2c F2). `n_testable` is structurally the complement of the two coverage absences -- the
# stopifnot above depends on that -- so it is left alone and the manuscript's rule is emitted as its
# own named quantity: a pair-level estimate exists, OR the pair has more than one assay. On the
# current set that is 47, not 49; the two ABO pairs with NO_EXACT_RSID_MATCH have neither an estimate
# nor a second assay, so they are non-absent without being testable in the manuscript's sense.
n_mr_testable <- sum(!is.na(all_summary$decode_b) | all_summary$n_assays > 1L)
cat(sprintf("INFO: %d pair(s) TESTABLE in the manuscript's sense (pair-level decode_b present, or n_assays > 1)\n",
            n_mr_testable))
# (The line that used to sit here announced the rate as being over TESTABLE pairs. It was left
#  behind when the rates moved to the evaluable denominator, so the script printed two
#  incompatible meanings for its own headline and invited 19 to be cited instead of 12.)

# EVALUABLE IS NARROWER THAN TESTABLE (objective review 2026-08-06, finding 2). `n_testable` removes
# only the two COVERAGE absences, so it still counts (i) pairs whose MR never produced an estimate --
# a harmonisation drop is not a replication failure, it is an untested pair -- and (ii) pairs whose
# deCODE instrument is too weak to carry the test at all (F < 10), where weak-instrument bias pulls
# the estimate back toward the confounded observational association. Reporting "replicated / testable"
# over that denominator states a replication FAILURE for pairs that were never evaluable.
# The column is DERIVED and additive: `status` is untouched, so nothing downstream changes meaning;
# what changes is which denominator is printed and that the reasons are counted separately.
# Per-gene assay coverage travels with the result, so a single-assay row can be told apart from a
# gene that should have had two (finding 4).
all_summary <- merge(all_summary, .assay_cov[, .(gene, n_assays_any_release, n_assays_in_2021, n_assays_usable, n_assays_missing)],
                     by = "gene", all.x = TRUE, sort = FALSE)
all_summary[is.na(n_assays_any_release), `:=`(n_assays_any_release = 0L, n_assays_in_2021 = 0L, n_assays_usable = 0L, n_assays_missing = 0L)]
# THE BOOLEAN `evaluable` IS REPLACED BY A THREE-VALUED CLASS (statistical + adversarial editor,
# 2026-08-07, concurring). A boolean cannot distinguish "no attempted assay cleared the gate" from
# "some did and some did not", so every gate reading it treated partial evidence as complete. The
# per-row evaluable_class is set where the row is built (none/some/all across ATTEMPTED assays);
# what is added here is the pair-level gate, which additionally requires a usable direction and a
# deCODE instrument strong enough to carry the test.
#
# EVERY GATE TESTS == "all". A pair whose second aptamer could not be evaluated is NOT evaluable,
# and a FALSE arising that way must be reported as "非可評估", never as "未複製" -- the two are
# different claims and only one of them is about replication.
all_summary[, evaluable_pair := evaluable_class == "all" &
                                !(status %chin% .ABSENT_STATUS) &
                                !is.na(direction_concordant) &
                                is.finite(decode_F_min) & decode_F_min >= 10]
# FIT ABSENCE AND ASSAY DISAGREEMENT ARE DIFFERENT EXCLUSIONS. `direction_concordant` is NA when
# no MR estimate exists, AND when the assays of one gene all fit but disagree in sign
# (ASSAY_DISCORDANT), AND when there is no UKB beta to compare against. Reading NA as "no MR
# estimate" put ADAM15 ~ eGFR (two assays, exposure betas 0.0456 and -0.0013, F_min 0.0127)
# into the no-MR count and left it out of the weak-instrument count. What records the presence
# of an estimate is decode_pval_iut: every no-fit row above writes it as NA.
# The three counts are mutually exclusive, in this order: no estimate -> weak instrument ->
# estimate with a strong instrument but no usable direction.
.ne_nomr  <- sum(!(all_summary$status %chin% .ABSENT_STATUS) & is.na(all_summary$decode_pval_iut))
.ne_weak  <- sum(!(all_summary$status %chin% .ABSENT_STATUS) & !is.na(all_summary$decode_pval_iut) &
                 (!is.finite(all_summary$decode_F_min) | all_summary$decode_F_min < 10))
.ne_dirna <- sum(!(all_summary$status %chin% .ABSENT_STATUS) & !is.na(all_summary$decode_pval_iut) &
                 is.finite(all_summary$decode_F_min) & all_summary$decode_F_min >= 10 &
                 is.na(all_summary$direction_concordant))
# COUNTED SEPARATELY, because it is the one exclusion that is NOT a statement about replication:
# these pairs had an evaluable assay AND a strong instrument AND a direction, and were excluded
# only because a SIBLING assay could not be evaluated.
.ne_partial <- sum(!(all_summary$status %chin% .ABSENT_STATUS) & all_summary$evaluable_class == "some" &
                   !is.na(all_summary$direction_concordant) &
                   is.finite(all_summary$decode_F_min) & all_summary$decode_F_min >= 10)
n_evaluable  <- sum(all_summary$evaluable_pair)
n_repl_eval  <- sum(all_summary$evaluable_pair & all_summary$status == "REPLICATED")
cat(sprintf(paste0("INFO: of those %d testable, %d are EVALUABLE (MR estimated AND deCODE F >= 10); ",
                   "%d had no MR estimate (harmonisation/fit), %d had a weak deCODE instrument, ",
                   "%d had an estimate and a strong instrument but no usable direction ",
                   "(assay disagreement or no UKB comparison).\n"),
            n_testable, n_evaluable, .ne_nomr, .ne_weak, .ne_dirna))
cat(sprintf(paste0("INFO: %d further pair(s) excluded ONLY because a SIBLING assay was not evaluable ",
                   "(evaluable_class == 'some'). 非可評估, NOT 未複製.
"), .ne_partial))
cat(sprintf("INFO: replicated among EVALUABLE pairs: %d of %d. The 'of %d testable' figure above must NOT be quoted as a replication rate.\n",
            n_repl_eval, n_evaluable, n_testable))

# A4: multi-assay genes (e.g. NOTCH2, FLT1) — report assay-row multiplicity (informational;
# not asserted to a fixed count, which would be stale if the strong-hit set changes)
for (g in intersect(c("NOTCH2", "FLT1"), all_assay$gene)) {
  cnt <- all_assay[gene == g, .N, by = outcome_id]
  cat(sprintf("INFO: %s assay rows per outcome: %s\n", g, paste(cnt$N, collapse = ", ")))
}

# A5: PCSK9 and IL6R exact rsID match + allele alignment
pcsk9_assay <- all_assay[gene == "PCSK9" & status == "MR_OK"]
il6r_assay  <- all_assay[gene == "IL6R" & status == "MR_OK"]
# 2026-06-30 (gate obj #11): positive controls — warn loudly if absent, do NOT hard-abort. Their
# absence under the corrected UKB-free set (e.g. protein not in deCODE SomaScan) is informative, not
# a pipeline error; a hard stop here would conflate deCODE non-coverage with a method failure.
if (nrow(pcsk9_assay) == 0) cat("WARNING: PCSK9 has NO deCODE MR result (positive control missing — investigate; not auto-fatal)\n")
if (nrow(il6r_assay)  == 0) cat("WARNING: IL6R has NO deCODE MR result (positive control missing — investigate; not auto-fatal)\n")
if (nrow(pcsk9_assay) > 0 && nrow(il6r_assay) > 0) cat("PASS: PCSK9 & IL6R have deCODE MR results\n")

# A7: Positive controls
pcsk9_row <- all_summary[gene == "PCSK9" & outcome_name == "CAD"]
il6r_row  <- all_summary[gene == "IL6R" & outcome_name == "CAD"]
cat(sprintf("Positive control PCSK9~CAD: status=%s, p=%s\n",
            pcsk9_row$status, format(pcsk9_row$decode_pval_iut, digits=3)))
cat(sprintf("Positive control IL6R~CAD: status=%s, p=%s\n",
            il6r_row$status, format(il6r_row$decode_pval_iut, digits=3)))

# ── 7. Save results ─────────────────────────────────────────────────────────

saveRDS(all_assay, file.path(OUT_DIR, "decode_mr_results.rds"))
# WRITE THEN RENAME, so a partial file never appears under the final name (round 6, finding 2).
# fwrite() wrote straight to the final path, so an interrupted or failing write left a parseable
# header plus a prefix of the 76 rows -- and cascade/13 checks only that the file exists, while
# phase4/12 turns every absent pair into decode_status = 'NO_DATA' / 'Not testable'. Lost rows
# would therefore read as ordinary missing deCODE evidence. Rename on one volume is atomic: the
# consumer sees the complete file or no file.
.sum_final <- file.path(OUT_DIR, "decode_replication_summary.csv")
.sum_tmp   <- paste0(.sum_final, ".partial")
fwrite(all_summary, .sum_tmp)
if (!file.rename(.sum_tmp, .sum_final)) {
  unlink(.sum_tmp)
  stop(sprintf("deCODE: could not move the completed summary into place (%s) -- refusing to leave a partial or absent artifact that cascade/13 would accept on existence alone", basename(.sum_final)))
}

cat("\n=== Final Summary ===\n")
cat(sprintf("Total pairs: %d\n", nrow(all_summary)))
cat(sprintf("NOT_RETRIEVED (in the 2021 release, absent from this mirror): %d\n", n_not_retrieved))
cat(sprintf("NOT_IN_RELEASE (not assayed in the 2021 release this analysis reads): %d\n", n_not_in_release))
cat("\nTestable pairs breakdown:\n")
print(all_summary[!status %chin% .ABSENT_STATUS, .N, by = status][order(-N)])

n_replicated <- sum(all_summary$status == "REPLICATED")
n_dir_only   <- sum(all_summary$status == "DIRECTION_ONLY")
n_weak       <- sum(all_summary$weak_instrument == TRUE, na.rm = TRUE)
n_hi_conf    <- sum(all_summary$high_confidence == TRUE, na.rm = TRUE)
# ONE DENOMINATOR, NOT THREE. The `evaluable` block above printed 11/12 and then this section
# printed the same numerator over `testable` (15/19), so the script contradicted itself and the
# larger, wrong figure was the one that read like the headline. Every RATE here is now over
# EVALUABLE pairs; the testable count is still reported, but as coverage, never as a
# denominator for replication.
# FREQUENCY COVERAGE IS REPORTED, NOT ASSUMED. A row that got no frequency reverts to the old
# behaviour for that row alone -- silently, if nobody prints this. The annotation is not
# byte-verified against any manifest (see 3b), so a truncated file would show up here as a miss
# rate and nowhere else.
if (length(.DC_EAF_SEEN) > 0L)
  cat(sprintf("\nVendor effect-allele frequency attached to %d of %d distinct deCODE variants (%.1f%%); variants without one keep the pre-2026-08-09 behaviour, i.e. palindromic ones among them are still discarded by harmonisation.\n",
              length(.DC_EAF_GOT), length(.DC_EAF_SEEN), 100 * length(.DC_EAF_GOT) / length(.DC_EAF_SEEN)))

# THE DECOMPOSITION IS EXHAUSTIVE OR IT IS A FALSE CLAIM. Every non-absent pair falls in exactly
# one bucket: no estimate / weak instrument / estimate and strong instrument but no usable
# direction / evaluable / excluded only because a SIBLING assay was not evaluable. This line
# named three of the five, so the printed sum could not reach the non-absent total.
cat(sprintf(paste0("\nEVALUABLE pairs (MR estimated AND deCODE F>=10): %d  [non-absent %d = ",
                   "evaluable %d + no-MR %d + weak-F %d + no-direction %d + sibling-not-evaluable %d]\n"),
            n_evaluable, n_testable, n_evaluable, .ne_nomr, .ne_weak, .ne_dirna, .ne_partial))
cat(sprintf("REPLICATED among evaluable: %d / %d\n", n_repl_eval, n_evaluable))
cat(sprintf("High-confidence among evaluable: %d / %d\n",
            sum(all_summary$evaluable_pair & all_summary$high_confidence %in% TRUE), n_evaluable))
cat(sprintf("DIRECTION_ONLY among evaluable: %d / %d\n",
            sum(all_summary$evaluable_pair & all_summary$status == "DIRECTION_ONLY"), n_evaluable))
cat(sprintf("(whole-table counts, NOT rates: REPLICATED %d, DIRECTION_ONLY %d over all %d pairs)\n",
            n_replicated, n_dir_only, nrow(all_summary)))
cat(sprintf("Weak deCODE instrument (F<10): %d\n", n_weak))

cat("\n=== REPLICATED but NOT high-confidence (weak instrument) ===\n")
weak_rep <- all_summary[status == "REPLICATED" & high_confidence == FALSE]
if (nrow(weak_rep) > 0) {
  print(weak_rep[, .(gene, outcome_name, decode_F_min, decode_pval_iut, delta_pval_iut)])
} else {
  cat("(none)\n")
}

cat("\nSaved to:\n")
cat(sprintf("  %s\n", file.path(OUT_DIR, "decode_mr_results.rds")))
cat(sprintf("  %s\n", file.path(OUT_DIR, "decode_replication_summary.csv")))
cat("\nDone.\n")
