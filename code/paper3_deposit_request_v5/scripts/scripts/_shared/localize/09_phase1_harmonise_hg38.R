# =============================================================================================
# Phase 1 (GRCh38), step 2 of 2: turn the raw extraction into the analysis input.
#
# WHY THIS EXISTS
# ---------------
# 08_phase1_extract_hg38.py emits the UKB-PPP columns as they come out of the archives
# (CHROM/GENPOS/ALLELE0/ALLELE1/A1FREQ/N/BETA/SE/LOG10P) plus provenance. Everything downstream of
# Phase 1 reads a HARMONISED table instead -- TwoSampleMR-style names, an rsID in SNP, an exposure
# label, and r2_exposure. In the shipped pipeline that gap was closed by a chain of one-off scripts
# (phase1/04, 05, 06, 12) that were written once and then frozen; 06 in particular now begins with
#     stop("Phase 1 frozen: 06_add_rsid_to_instruments.R must not be re-run.")
# and its own comment says to remove that only when rebuilding Phase 1 from scratch. That is exactly
# what is happening, but rather than unfreeze a historical script and have it overwrite the artifact
# that documents the submitted analysis, the same transformation is done here, forward-only.
#
# THE rsID PROBLEM, AND WHY IT IS NOT A METADATA PROBLEM
# ------------------------------------------------------
# Nothing here fetches. The downstream extractors keep only instruments whose SNP matches
# ^rs[0-9]+$ (02_extract_outcomes_local.R:26, cascade/01_mr_screen.R:90 and :143), so an instrument
# without an rsID is DISCARDED and its protein can leave the analysis entirely.
#
# Section 2 therefore uses TWO local sources, both keyed on a VARIANT rather than a coordinate:
#   - ukb_rsid_map_by_id.tsv -- UKB-PPP's own per-chromosome annotation, keyed on the full variant ID
#     (chr:pos:ref:alt:imp:v1) this extraction already carries. Built once by
#     scripts/_shared/phase1/04_download_rsid_map.py; REQUIRED and pinned to an expected md5 here.
#   - the 1000G EUR .bim, matched through the GRCh37 coordinate embedded in the UKB-PPP ID.
# A THIRD source -- the previous instrument table, read by literal path from
# results/repro_phase1_1682/ -- was removed 2026-08-02 so that a from-scratch reproduction needs no
# retired generation on disk. See the block at Section 2 for what that cost and what it did not.
#
# The UKB annotation was added 2026-07-31 and IS part of the contract, not an optional extra: it is
# the only rsID evidence for 66 instruments covering 66 proteins, so with it this extraction yields
# 1,964 instruments over 1,948 proteins and without it 1,898 over 1,882. Its absence is an error, not
# a fallback -- otherwise the same certified extraction would give two answers 66 proteins apart
# depending on the state of the disk.
#
# The SUPERSEDED data/ukb_ppp/ukb_rsid_map.tsv is NOT one of them: it is position-only, its builder
# discarded the alleles, and it holds 1,587 positions chosen for the OLD GRCh37-derived windows.
# Because both sources identify a variant rather than a coordinate, a disagreement between them is
# two NAMES for one variant, not two variants -- so both names are kept rather than discarding both
# and losing the protein over a naming difference. SNP takes the UKB annotation's name
# and rsid_ld_panel carries the reference panel's, because the two
# downstream consumers of SNP work in different identifier spaces and one column cannot serve both;
# every case is written out. The run does NOT stop on a coverage figure: the 90% / 2% gates that
# used to be here were invented in this file, and a threshold on a result cannot make an unnamed
# variant named -- it can only refuse to write, and the only way to pass it is to move the number.
# What IS refused is a broken join (not one instrument received an rsID), which is a failure of this
# script rather than a property of the data. Every instrument without a usable rsID is listed with
# its reason, and the genes that lose all of theirs are named.
#
# (The 95.4% the shipped run achieved was not the map alone: 49 of those rsIDs came from an
# allele-aware .bim backfill in 15_backfill53_build_instruments.R. The map by itself gives 92.6%.)
#
# LOCAL-FIRST: reads only local files. No network.
# Run: Rscript --no-init-file scripts/_shared/localize/09_phase1_harmonise_hg38.R
# =============================================================================================
# THIS file, wherever it was invoked from, hashed before it does anything. The certificate used to
# hash a hard-coded canonical path at the end of the run, which records the canonical file even when
# a modified copy is what executed, and misses an edit made while the run was in progress.
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
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  # It must be THIS file. --file= names whatever Rscript was pointed at, so a wrapper that sources
  # this script recorded the wrapper's hash under builder = "09_phase1_harmonise_hg38.R", and an
  # invocation with no --file recorded NA and installed the output anyway.
  if (!length(f) || !file.exists(f[1]))
    stop("09: cannot identify the file being executed (no --file= argument), so the certificate ",
         "cannot record which version of this script produced the output. Run it with Rscript.")
  if (!identical(basename(normalizePath(f[1], winslash = "/")), "09_phase1_harmonise_hg38.R"))
    stop("09: this script is being sourced from ", basename(f[1]), ", so builder_md5 would name ",
         "that file while the certificate names 09_phase1_harmonise_hg38.R. Run it directly.")
  unname(tools::md5sum(normalizePath(f[1])))
})

# NOT _phase1_set.R. This script is the PRODUCER of the new instrument generation, not a consumer of
# a selected one: it uses neither PHASE1_SET nor PHASE1_INSTRUMENTS. Sourcing the selector pulled in its rebuild guard,
# which exists to stop CONSUMERS running on a stale generation -- and it stopped this file, the one
# whose whole purpose is to produce the generation that guard is waiting for. The guard was right
# about the previous table (382 (gene, chrpos) keys and 31 genes absent from the rebuild); it simply
# had no business being asked here.
suppressPackageStartupMessages(library(data.table))
setwd(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"))

IN    <- "results/phase1_hg38/ukb_cis_instruments.csv"
MANI  <- "results/phase1_hg38/phase1_OUTPUTS.tsv"
RSMAP <- "data/ukb_ppp/ukb_rsid_map.tsv"
OUTD  <- "results/phase1_hg38"
OUT   <- file.path(OUTD, "ukb_cis_instruments.rds")


# ---- 0. the lock comes FIRST -----------------------------------------------------------------
# Before ANY destructive step. It used to be taken just before the install, hundreds of lines after
# this script had already deleted the previous outputs and written its sidecars -- so two harmonisers
# started together both deleted, both computed, and only then did one of them discover it could not
# proceed. By that point the other process's outputs were already gone.
.lock <- file.path(OUTD, ".phase1_run.lock")
if (!dir.exists(OUTD))
  stop(OUTD, " does not exist -- run 08_phase1_extract_hg38.py first")
if (!dir.create(.lock, showWarnings = FALSE)) {
  .who <- paste(list.files(.lock), collapse = ", ")
  stop("another Phase 1 process holds ", .lock, " (", .who, "). 08 finalising while this script ",
       "runs would leave an RDS derived from the previous extraction beside a newer certified ",
       "one. Wait for it, or remove the directory if no such process is alive.")
}
# REGISTERED BEFORE the PID file is written. With the order reversed, a failure to write inside the
# lock directory -- a full disk, a permission change -- left the directory created and no finalizer
# yet attached, so the lock outlived the run and blocked every later one with a PID that never held it.
# reg.finalizer, NOT on.exit(): this script runs directly under Rscript, so there is no function
# frame for on.exit() to attach to and the call is a silent no-op -- the lock would then outlive
# every run and block the next one. (The same mistake was found and fixed in 03b.)
reg.finalizer(environment(), function(e) unlink(.lock, recursive = TRUE), onexit = TRUE)
writeLines(paste0(Sys.getpid(), "\n", format(Sys.time())),
           file.path(.lock, paste0("pid_09_", Sys.getpid())))

# ---- 1. refuse anything the extractor did not certify ---------------------------------------------
# 08 writes phase1_OUTPUTS.tsv last, and refuses to write it at all if any archive could not be read.
# Its absence therefore means the extraction did not finish or was not complete, and harmonising such
# a table would launder an uncertified instrument set into the analysis input.
if (!file.exists(IN))   stop("missing ", IN, " -- run 08_phase1_extract_hg38.py first")
if (!file.exists(MANI)) stop(MANI, " is absent, so the extraction is NOT certified complete. ",
                             "Harmonising it would present an incomplete instrument set as the ",
                             "analysis input. Re-run the extraction.")
# READ ONCE, then hashed and parsed from the SAME bytes -- the pattern this script already uses for
# the instrument CSV and the .bim. Parsing the path here and hashing it eighty lines later left a
# window in which the certificate could be replaced between the two, so the output would record the
# hash of one manifest while every decision above came from another.
.mani_bytes <- readBin(MANI, "raw", file.info(MANI)$size)
MANI_MD5 <- digest::digest(.mani_bytes, algo = "md5", serialize = FALSE)
mani <- fread(text = rawToChar(.mani_bytes), colClasses = "character")
# A smoke extraction writes phase1_SMOKE_OUTPUTS.tsv rather than this file, so it cannot normally
# reach here at all. This refuses the marker anyway, in the manifest or in the table itself: a file
# renamed into place by hand should still not be able to enter the analysis as an exposure set.
if ("SMOKE_TEST_NOT_AN_EXPOSURE_SET" %in% mani$file)
  stop(basename(MANI), " carries a SMOKE_TEST marker: it certifies a deliberately limited test ",
       "extraction, not an instrument set. Run 08 without P1_LIMIT.")
row  <- mani[file == basename(IN)]
if (nrow(row) != 1L) stop(basename(MANI), " does not describe ", basename(IN))
if (!identical(tolower(unname(tools::md5sum(IN))), tolower(row$md5[1])))
  stop(basename(IN), " does not match its certificate; it changed after the extraction wrote it.")

# Any exit from here on leaves the previous outputs gone rather than plausibly current. No consumer
# inspects this script's certificate, so a stale RDS beside a failed rerun reads as canonical;
# removing it first makes a failure legible as a missing input instead of a silent old answer.
# The removal is CHECKED -- an unlink that fails and is ignored leaves exactly the state this is
# meant to prevent, while the comment above says it cannot happen.
.stale <- c(OUT, paste0(OUT, ".manifest.tsv"), file.path(OUTD, "instruments_without_rsid.csv"),
            file.path(OUTD, "rsid_name_disagreements.csv"))
.stale <- .stale[file.exists(.stale)]
if (length(.stale)) {
  unlink(.stale)
  .left <- .stale[file.exists(.stale)]
  if (length(.left))
    stop("could not remove the previous output(s): ", paste(basename(.left), collapse = ", "),
         ". If this run then failed, they would be indistinguishable from a fresh result. ",
         "Close whatever holds them and rerun.")
}

# Read ONCE, and hashed from the bytes that were read. Hashing the path above and calling fread()
# on it again is two reads: 08 replaces this CSV before writing its new manifest, so a concurrent
# extraction could have its old table validated here and its new one parsed -- and if that
# extraction then failed without certifying, this script would harmonise it under the old hash.
.in_bytes <- readBin(IN, "raw", file.info(IN)$size)
# The hash OF THOSE BYTES, kept for the provenance the output carries. Hashing the path again later
# records whatever is there then, so a CSV replaced after this read leaves the table holding results
# from the original bytes beside an attribute naming the replacement.
IN_MD5 <- digest::digest(.in_bytes, algo = "md5", serialize = FALSE)
if (!identical(tolower(IN_MD5), tolower(row$md5[1])))
  stop(basename(IN), " changed between its certificate check and this read.")
# The certificate declares a ROW COUNT as well as a hash, and it was never compared -- a manifest
# saying rows=1 beside the right md5 for a 1,700-row table passed. It also lists the summary that
# records WHICH extractor and which coordinate table produced this CSV, and that was neither
# required nor hashed, so the output could not establish either.
if ("rows" %in% names(mani) && !is.na(row$rows[1]) && nzchar(row$rows[1])) {
  # A WHOLE NUMBER, spelled as one. as.integer() reads "663.9" as 663 and "6.63e2" as 663 too, so a
  # certificate that does not actually declare an integer row count passed while looking checked.
  if (!grepl("^[0-9]+$", trimws(row$rows[1])))
    stop(basename(MANI), " declares ", sQuote(row$rows[1]), " rows for ", basename(IN),
         ", which is not a whole number.")
  .declared <- suppressWarnings(as.integer(trimws(row$rows[1])))
  if (is.na(.declared))
    stop(basename(MANI), " declares a non-numeric row count for ", basename(IN), ".")
} else {
  stop(basename(MANI), " does not declare a row count for ", basename(IN), ", so the certificate ",
       "cannot say whether the table is complete.")
}
# ...and COMPARED. Requiring the field and then not checking it against the parsed table left a
# manifest declaring rows=1 beside the correct md5 for a 1,700-row CSV passing exactly as before.
# The comparison happens after the parse, below.
SUMMARY_F <- file.path(OUTD, "phase1_summary.json")
if (!basename(SUMMARY_F) %in% mani$file)
  stop(basename(MANI), " does not certify ", basename(SUMMARY_F), ", which is the only record of ",
       "which extractor version and which coordinate table produced this CSV.")
if (!file.exists(SUMMARY_F))
  stop(SUMMARY_F, " is certified by the manifest but absent.")
# Read once, hashed and parsed from the same bytes: hashing the path and re-parsing it below left a
# window in which the summary could be replaced between the two, so every extractor identity read
# out of it would come from a file other than the one certified here.
.sum_bytes <- readBin(SUMMARY_F, "raw", file.info(SUMMARY_F)$size)
SUMMARY_MD5 <- digest::digest(.sum_bytes, algo = "md5", serialize = FALSE)
.srow <- mani[file == basename(SUMMARY_F)]
if (nrow(.srow) != 1L || !identical(tolower(SUMMARY_MD5), tolower(.srow$md5[1])))
  stop(basename(SUMMARY_F), " does not match its certificate.")
# The extractor's own identifiers, carried into this script's provenance so the RDS can name the
# run that produced its input rather than only the bytes.
.sj <- tryCatch(jsonlite::fromJSON(rawToChar(.sum_bytes)), error = function(e) NULL)
# 08 puts code_md5 inside its RUN SIGNATURE payload, which is hashed rather than written out; the
# summary carries run_signature and run_signature_inputs_except_archive_fingerprints. Reading a key
# that is not there recorded NA for every genuine 08 output while looking like a real check, so the
# fields that DO exist are read and their absence is fatal rather than silent.
EXTRACTOR_MD5 <- {
  # NESTED, under run_signature_inputs_except_archive_fingerprints. Read at the top level it was
  # NA for every genuine 08 output, so the field looked like a check and recorded nothing.
  v <- if (is.null(.sj)) NULL else
       (.sj[["code_md5"]] %||% .sj[["run_signature_inputs_except_archive_fingerprints"]][["code_md5"]])
  if (is.null(v)) NA_character_ else as.character(v)
}
if (is.na(EXTRACTOR_MD5) || !nzchar(EXTRACTOR_MD5))
  stop(basename(SUMMARY_F), " records no extractor code hash, so the version of 08 that produced ",
       basename(IN), " cannot be identified.")
# ...AND IT MUST BE THE EXTRACTOR ON DISK. Requiring the field to be merely present identified
# nothing: "" satisfied it, and so did a summary written by any other version of 08 -- one with
# different QC thresholds, a different clumping rule or a different panel decision, all of which
# produce a self-consistent CSV, summary and manifest that this script would accept as the Phase-1
# instrument set. The comparison is what makes the recorded hash a check rather than a note.
EXTRACTOR_F <- "scripts/_shared/localize/08_phase1_extract_hg38.py"
if (!file.exists(EXTRACTOR_F))
  stop(EXTRACTOR_F, " is not on disk, so the extraction cannot be tied to the code that produced it.")
if (!identical(tolower(EXTRACTOR_MD5), tolower(unname(tools::md5sum(EXTRACTOR_F)))))
  stop(basename(IN), " was produced by a different extractor than the one on disk (summary ",
       substr(EXTRACTOR_MD5, 1, 12), ", ", basename(EXTRACTOR_F), " ",
       substr(unname(tools::md5sum(EXTRACTOR_F)), 1, 12), "). The thresholds, the clumping rule and ",
       "the panel decision behind these instruments are therefore not the ones in the current ",
       "script. Re-run 08, or check out the extractor that produced this table.")
EXTRACT_RUN_SIG <- {
  v <- if (!is.null(.sj)) (.sj[["run_signature"]] %||% .sj[["run_sig"]]) else NULL
  if (is.null(v)) NA_character_ else as.character(v)
}
# The run signature is what identifies the extraction, and 08 always writes it. Its absence means
# the summary was not written by this extractor.
if (is.na(EXTRACT_RUN_SIG))
  stop(basename(SUMMARY_F), " carries no run signature, so the extraction that produced ",
       basename(IN), " cannot be identified. It was not written by 08_phase1_extract_hg38.py.")
# The liftover chain is an input to the coordinates this table inherits, and 08 records it.
EXTRACT_CHAIN_MD5 <- {
  v <- if (!is.null(.sj)) .sj[["liftover_chain_md5"]] else NULL
  if (is.null(v)) NA_character_ else as.character(v)
}
# REQUIRED, not optional. Every GRCh37<->GRCh38 relationship in this pipeline -- the ID-to-GENPOS
# check, the coordinate table, the windows downstream queries -- rests on one chain file, and a
# summary that omits its hash produces an output whose provenance cannot name it.
if (is.na(EXTRACT_CHAIN_MD5))
  stop(basename(SUMMARY_F), " records no liftover_chain_md5, so the chain behind every ",
       "GRCh37/GRCh38 relationship in ", basename(IN), " cannot be identified.")

# ── the gene coordinates are GRCh38 because they were CHECKED, not because it is written here ──
# This script stamped gene_coord_build = "GRCh38" unconditionally while inspecting nothing about
# gene_chrom/start/end. A certified CSV carrying A1BG's historical GRCh37 interval
# 19:58,856,544-58,864,865 instead of the GRCh38 19:58,345,178-58,353,502 passes every ID/GENPOS
# build test -- those look at variants, not at gene bodies -- and the output then presents a GRCh37
# window as GRCh38, so colocalisation queries the wrong region.
EXTRACT_BUILD  <- { v <- if (!is.null(.sj)) .sj[["genome_build"]] else NULL
                    if (is.null(v)) NA_character_ else as.character(v) }
EXTRACT_COORDS <- { v <- if (!is.null(.sj)) .sj[["coords_file"]] else NULL
                    if (is.null(v)) NA_character_ else as.character(v) }
EXTRACT_COORDS_MD5 <- { v <- if (!is.null(.sj)) .sj[["coords_md5"]] else NULL
                        if (is.null(v)) NA_character_ else as.character(v) }
# ANCHORED. A substring test passes on "not GRCh38", "GRCh38-like" and "lifted from GRCh38", each of
# which is a statement that the coordinates are NOT plainly GRCh38. The field 08 writes is
# "GRCh38 (GENCODE v50)", so the build token has to start the string and end at a word boundary --
# which pins the build without pinning the annotation release into this script.
# ...and \b is not the anchor this needs. In R, "-" is a non-word character, so "^GRCh38\\b" matches
# "GRCh38-like" -- a boundary exists between "8" and "-". The label must END there or continue with
# something that cannot make it a different build name.
if (is.na(EXTRACT_BUILD) || !grepl("^GRCh38($|[^A-Za-z0-9_-])", EXTRACT_BUILD))
  stop(basename(SUMMARY_F), " records genome_build = ",
       if (is.na(EXTRACT_BUILD)) "<absent>" else sQuote(EXTRACT_BUILD),
       ", not GRCh38. This script labels its gene coordinates GRCh38 and everything downstream ",
       "believes that label.")
if (is.na(EXTRACT_COORDS_MD5))
  stop(basename(SUMMARY_F), " records no coords_md5, so the coordinate table behind ", basename(IN),
       " cannot be identified.")
COORDS_TBL <- file.path("data/annotation", basename(EXTRACT_COORDS))
if (!file.exists(COORDS_TBL))
  stop("the extraction was built against ", basename(EXTRACT_COORDS), ", which is not on disk, so ",
       "its gene coordinates cannot be verified.")
# Same read-once rule. The coordinate table was hashed from its path and then re-parsed from that
# path, so the gene intervals every window in the output rests on could come from a file other than
# the one whose hash was checked -- and the coordinate builder is not covered by the Phase-1 lock.
.coords_bytes <- readBin(COORDS_TBL, "raw", file.info(COORDS_TBL)$size)
COORDS_MD5 <- digest::digest(.coords_bytes, algo = "md5", serialize = FALSE)
if (!identical(tolower(COORDS_MD5), tolower(EXTRACT_COORDS_MD5)))
  stop(basename(COORDS_TBL), " is not the table the extraction used (summary ",
       substr(EXTRACT_COORDS_MD5, 1, 12), ", on disk ", substr(COORDS_MD5, 1, 12), ").")
raw <- fread(text = rawToChar(.in_bytes), colClasses = list(character = c("CHROM", "chrpos")))
.coords_ref <- fread(text = rawToChar(.coords_bytes), colClasses = list(character = "chrom"))
# KEYED ON oid, THE ASSAY, not on the gene symbol. Two Olink assays can measure the same protein at
# the same variant with different beta/SE -- CXCL8 has OID20631 and OID21430 -- so a row whose oid is
# wrong but whose symbol and interval are right passed every check here, and resolve_ukbppp_tar()
# downstream then opened the OTHER assay's archive: MR on one assay, colocalisation on another.
# The identity fields are compared for the same reason: they are copied into the output and used
# downstream, and nothing was checking them at all.
# assay_file is the one that actually settles it. Keying on oid alone cannot catch a mislabelled
# CXCL8 row, because OID20631 and OID21430 are both genuinely CXCL8 at the same locus and every
# other field agrees -- but the ARCHIVE the numbers were read out of does not, and the reference
# assigns exactly one file to each oid.
.idc <- c(assay_file = "file", ensg = "ensg", uniprot = "acc",
          assignment_basis = "assignment_basis", alt_locus_ensg = "alt_locus_ensg",
          alt_locus_start = "alt_locus_start", alt_locus_end = "alt_locus_end")
# REQUIRED, not intersected away. Reducing the identity fields to whatever both tables happen to
# carry meant a coordinate table missing `file`, `acc` or `ensg` silently stopped checking that
# field -- and those are the fields that decide WHICH archive a row's numbers came from, so a row
# carrying another assay's beta, se, frequency and sample size under this protein's identity would
# pass. The extracted CSV's columns are required a few dozen lines below; the reference's were not
# required anywhere.
.miss_ref <- setdiff(unname(.idc), names(.coords_ref))
.miss_raw <- setdiff(names(.idc), names(raw))
if (length(.miss_ref) || length(.miss_raw))
  stop("the identity cross-check cannot run: ", basename(COORDS_TBL), " lacks ",
       if (length(.miss_ref)) paste(.miss_ref, collapse = ", ") else "nothing",
       " and the extraction lacks ",
       if (length(.miss_raw)) paste(.miss_raw, collapse = ", ") else "nothing",
       ". These fields decide which archive each row's numbers were read from; without them a row ",
       "carrying another assay's estimates under this protein's identity is undetectable.")
.blank <- function(x) { x <- trimws(as.character(x)); x[is.na(x) | x %in% c("", '""')] <- ""; x }
.L <- unique(raw[, c("oid", "gene", "gene_chrom", "gene_start", "gene_end", names(.idc)), with = FALSE])
.R <- unique(.coords_ref[, c("oid", "gene", "chrom", "start", "end", unname(.idc)), with = FALSE])
setnames(.R, c("oid", "r_gene", "r_chrom", "r_start", "r_end", paste0("r_", names(.idc))))
.chk <- merge(.L, .R, by = "oid", all.x = TRUE)
.bad <- .chk[is.na(r_start) | .blank(gene) != .blank(r_gene) |
             .blank(gene_chrom) != .blank(r_chrom) |
             as.integer(gene_start) != r_start | as.integer(gene_end) != r_end]
for (.k in names(.idc))
  .bad <- rbind(.bad, .chk[!is.na(r_start) & .blank(get(.k)) != .blank(get(paste0("r_", .k)))])
.bad <- unique(.bad)
if (nrow(.bad))
  stop(nrow(.bad), " assay(s) carry a gene symbol, interval or identity that does not match ",
       basename(COORDS_TBL), " -- the table the extraction says it used. The build label this ",
       "script writes would then be a claim about a window that came from somewhere else, and the ",
       "identity fields decide which archive downstream opens. First: ",
       paste(utils::head(.bad[, paste0(oid, " ", gene, " ", gene_chrom, ":", gene_start, "-",
                                       gene_end, " vs ", r_gene, " ", r_chrom, ":", r_start, "-",
                                       r_end)], 3), collapse = "; "), ". Nothing was written.")
cat(sprintf("gene coordinates and assay identity: all %d assay(s) match %s (%d identity field(s))\n",
            nrow(.chk), basename(COORDS_TBL), length(.idc)))

# chrpos IS CHROM:GENPOS, and is asserted to be. It independently keys the historical rsID map
# below, and a conflict there is deliberately resolved in that map's favour -- so a chrpos naming a
# different locus than the row's own CHROM/GENPOS silently renames the variant. On the A1BG row,
# setting chrpos to 19:54301977 while leaving CHROM/GENPOS/ID correct yields LILRA3's rs431001.
# NA COUNTS AS A VIOLATION. which() drops NA comparisons, so a row with a missing chrpos, CHROM or
# GENPOS passed this test silently and could still take an rsID from its ID -- and no later assertion
# required chrpos to be present either.
.cp <- which(raw$chrpos != paste0(raw$CHROM, ":", raw$GENPOS) |
             is.na(raw$chrpos) | is.na(raw$CHROM) | is.na(raw$GENPOS))
if (length(.cp))
  stop(length(.cp), " row(s) have a chrpos that is not CHROM:GENPOS (e.g. ",
       paste(utils::head(raw[.cp, paste0(gene, " ", chrpos, " vs ", CHROM, ":", GENPOS)], 3),
             collapse = "; "), "). chrpos keys the historical rsID map, and that map wins a ",
       "disagreement, so the row would be named after whatever variant sits at the other position.")

# THE ID's OWN FIELDS. ID is chr:pos37:ref:alt:imp:v1 and only the first two were ever read, so its
# chromosome and alleles were never compared with the row's. An ID naming a different chromosome
# resolves against the panel to that chromosome's variant, and swapping ALLELE0/ALLELE1 without
# reversing BETA still matched through the unordered key -- emitting the opposite effect allele
# against an unchanged effect size. All 663 rows of the current partial extraction satisfy both.
.idf <- tstrsplit(as.character(raw$ID), ":", fixed = TRUE)
if (length(.idf) < 4L)
  stop("the ID column is not chr:pos:ref:alt (first value ", sQuote(raw$ID[1]), ").")
# PER ROW, and NA counts as a violation. The width test above is global: in a table where one ID
# carries only two fields, tstrsplit() still returns four vectors and that row's allele comparisons
# are NA -- which which() drops, so the malformed row passed the check written to catch it.
.ib <- which(.idf[[1]] != raw$CHROM | toupper(.idf[[3]]) != toupper(raw$ALLELE0) |
             toupper(.idf[[4]]) != toupper(raw$ALLELE1) |
             is.na(.idf[[1]]) | is.na(.idf[[2]]) | is.na(.idf[[3]]) | is.na(.idf[[4]]))
if (length(.ib))
  stop(length(.ib), " row(s) have an ID whose chromosome or ordered alleles disagree with the ",
       "row's own CHROM/ALLELE0/ALLELE1 (e.g. ",
       paste(utils::head(raw[.ib, paste0(gene, " ", ID, " vs ", CHROM, " ", ALLELE0, "/", ALLELE1)],
                         3), collapse = "; "),
       "). The ID resolves the rsID, so a wrong chromosome names another variant, and a swapped ",
       "allele pair emits the opposite effect allele against an unchanged BETA.")

if (nrow(raw) != .declared)
  stop(basename(MANI), " declares ", .declared, " row(s) for ", basename(IN), " but it parses to ",
       nrow(raw), ". The certificate and the table describe different extractions.")
# AND AGAINST THE SUMMARY, which is a second, independently certified statement of the same two
# numbers. Checking only the manifest left one editable path open: delete a mapped single-instrument
# protein's row from the CSV and update that one file's rows+md5 in phase1_OUTPUTS.tsv. Every check
# above then agrees, the deleted protein is absent from every denominator the coverage gates use,
# and the run emits an apparently complete RDS that is silently missing it. The summary is hashed by
# the same manifest, so an edit now has to be made consistently in three places instead of one.
.sn <- suppressWarnings(as.integer(.sj[["n_snps_total"]] %||% NA))
.sg <- suppressWarnings(as.integer(.sj[["n_genes_ok"]] %||% NA))
if (is.na(.sn) || is.na(.sg))
  stop(basename(SUMMARY_F), " records no n_snps_total/n_genes_ok, so the table cannot be checked ",
       "against the extraction's own account of what it produced.")
if (nrow(raw) != .sn || uniqueN(raw$ensg) != .sg)
  stop(basename(SUMMARY_F), " says the extraction produced ", .sn, " row(s) over ", .sg,
       " gene(s); ", basename(IN), " holds ", nrow(raw), " row(s) over ", uniqueN(raw$ensg),
       ". The certified summary and the certified table describe different extractions.")
if ("SMOKE_TEST_NOT_AN_EXPOSURE_SET" %in% names(raw))
  stop(basename(IN), " carries a SMOKE_TEST_NOT_AN_EXPOSURE_SET column: it came from a limited ",
       "test extraction, not a full one. Run 08 without P1_LIMIT.")
cat(sprintf("read %s certified instruments (%s genes)\n",
            format(nrow(raw), big.mark = ","), format(uniqueN(raw$ensg), big.mark = ",")))

need <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P",
          "F_stat", "pval", "chrpos", "gene", "ensg", "oid", "uniprot", "gene_chrom", "gene_start",
          "gene_end", "assignment_basis", "id_build", "pval_underflow",
          "alt_locus_ensg", "alt_locus_start", "alt_locus_end", "cis_start", "cis_end",
          # REQUIRED. It was optional, so when it was absent the identity cross-check below silently
          # dropped the one field that settles which archive a row's numbers were read out of --
          # OID20631 and OID21430 are both genuinely CXCL8 at the same locus, and every other field
          # agrees between them.
          "assay_file")
miss <- setdiff(need, names(raw))
if (length(miss)) stop("the extraction is missing columns: ", paste(miss, collapse = ", "))

# ONE ASSAY PER PROTEIN, asserted here rather than assumed from 08's panel rule. 08 fixes the panel
# in advance and records the decision, but nothing downstream of that re-checks it: two valid panels
# of one protein contributing distinct lead variants keep (ensg, SNP) unique and pass every other
# assertion in this file, while downstream groups both rows under one `gene`/`id.exposure` and
# resolve_ukbppp_tar() (R/00_config.R) opens `.oid[1]`. MR would then combine two assays' estimands
# while colocalisation read only the first assay's archive. The pre-final partial extraction really
# did carry four CXCL8 OIDs across three distinct variants, so this is the observed shape, not a
# hypothetical one.
.asy <- raw[, .(n_oid = uniqueN(oid), n_file = uniqueN(assay_file)), by = ensg][n_oid > 1L | n_file > 1L]
if (nrow(.asy))
  stop(nrow(.asy), " protein(s) carry instruments from more than one assay (e.g. ",
       paste(utils::head(raw[ensg %in% .asy$ensg,
                            paste0(gene, ": ", paste(sort(unique(oid)), collapse = "/")), by = ensg]$V1,
                         3), collapse = "; "),
       "). 08 fixes one panel per protein in advance; a table with two means the panel decision did ",
       "not survive into it, and MR and colocalisation can then read different assays for one gene.")

# THE WINDOW THE INSTRUMENTS ACTUALLY CAME FROM, carried rather than recomputed. 08 truncates the
# window at the chromosome end (A1BG: 57,845,178-58,617,616, while gene_end + 500 kb would be
# 58,853,502), so a consumer that rebuilds `gene +/- CIS_KB` for its own regional query does not get
# the interval the exposure was defined on. Required rather than optional, and with no fallback that
# recomputes it: an absent column has to stop the run, or the truncation silently disappears again.
.cis_kb <- suppressWarnings(as.integer(
  .sj$run_signature_inputs_except_archive_fingerprints$cis_kb))
if (length(.cis_kb) != 1L || is.na(.cis_kb) || .cis_kb <= 0L)
  stop(basename(SUMMARY_F), " records no usable cis_kb, so the window width behind cis_start/",
       "cis_end cannot be checked against the certified extraction.")
.cw <- raw[, .(nwin   = uniqueN(paste(cis_start, cis_end)),
               # THE CHROMOSOME TOO. The window test is numeric, so a chromosome-18 variant whose
               # position happens to fall inside a chromosome-19 protein's interval satisfied it and
               # entered the table as a cis instrument for that protein. Nothing else in this script
               # related the variant's CHROM to the protein's gene_chrom -- the coordinate cross-check
               # only establishes that gene_chrom itself is right.
               off_chrom = sum(as.character(CHROM) != as.character(gene_chrom)),
               outside = sum(GENPOS < cis_start | GENPOS > cis_end),
               bad_lo  = sum(cis_start != pmax(1L, gene_start - .cis_kb * 1000L)),
               # Only the UPPER bound is truncated by a chromosome length this script does not hold,
               # so it is bounded rather than equated; the lower bound is exact. STATED LIMIT: for a
               # truncated assay -- 4 of the 655 in the current extraction (A1BG, ARSA, DLL1, ECHS1) --
               # an upper bound that is wrong but still inside [gene_end, gene_end + CIS_KB] passes
               # here. Duplicating GRCh38 chromosome lengths into this script to close that would
               # create a second source of truth for the quantity 08 already applied; the CSV's own
               # md5 is checked against the extraction certificate above, which is what makes the
               # carried value trustworthy. This check is consistency evidence, not anti-tampering.
               bad_hi  = sum(cis_end < gene_end | cis_end > gene_end + .cis_kb * 1000L)), by = oid]
if (any(.cw$nwin != 1L) || sum(.cw$off_chrom, .cw$outside, .cw$bad_lo, .cw$bad_hi) > 0L)
  stop(sum(.cw$nwin != 1L), " assay(s) carry more than one cis window, ", sum(.cw$off_chrom),
       " row(s) sit on a chromosome other than their own protein's, ", sum(.cw$outside),
       " row(s) sit outside their own window, and ", sum(.cw$bad_lo, .cw$bad_hi), " row(s) have a ",
       "window that is not gene +/- ", .cis_kb, " kb truncated at the chromosome ends. The window ",
       "is a property of the assay and decides every regional query downstream.")
cat(sprintf("cis windows: %d assay(s), one window each, every GENPOS inside it\n", nrow(.cw)))

# ---- 2. rsID, from two local sources: UKB's own annotation and the LD panel ----
# This matters more than a metadata field should. The downstream extractors keep only instruments
# whose SNP matches ^rs[0-9]+$ (02_extract_outcomes_local.R:26, cascade/01_mr_screen.R:90 and :143), so an
# instrument that falls back to a chr:pos string is not "joined by position" as an earlier version of
# this comment claimed -- it is DISCARDED, and its protein leaves the analysis. Coverage is therefore
# part of the result, not a diagnostic.
#
# Two local sources, neither requiring the network, BOTH keyed on a variant rather than a position:
#   A  ukb_rsid_map_by_id.tsv, UKB-PPP's own annotation, joined on the full variant ID this
#      extraction already carries. (Not ukb_rsid_map.tsv -- see below.)
#   B  the 1000G EUR .bim, matched through the GRCh37 coordinate embedded in the UKB-PPP ID and
#      keyed on the unordered allele pair. The .bim is GRCh37 and so is that ID, which is why they
#      match at all; measured coverage of QC-passing candidates is ~85%.
# Where both name a variant and the names DIFFER, BOTH names are kept: each has matched the
# same variant identity, so the disagreement is about naming (merged or retired rsIDs), not about
# which variant is meant. SNP takes A's name
# and rsid_ld_panel carries B's, because the two downstream consumers of SNP work in different
# identifier spaces (02_extract_outcomes_local.R matches outcome rows by SNP; 03_build_ld_proxies.R
# hands SNP to PLINK against this panel) and one column cannot serve both. Every such case is written
# to rsid_name_disagreements.csv.
# NOT ukb_rsid_map.tsv. That file looks like a position annotation but is not one: its builder
# (scripts/_shared/phase1/04_download_rsid_map.py) matched each OLD instrument's full ID to the UKB
# annotation and then DISCARDED the ID and the alleles, keeping only chrpos -> rsid. So a row means
# "the rsID of whichever variant the old extraction selected at this position", not "the rsID at this
# position". The new extraction can select a different allele pair there, and nothing in the file can
# say so -- against the closest local raw extraction, 195 rows would take an rsID on that basis.
# Checking that the NEW table has only one allele pair at the position does not help either: it
# establishes that the new table is unambiguous, not that its variant is the one the map meant.
#
# THE PREVIOUS INSTRUMENT TABLE WAS REMOVED AS AN rsID SOURCE, 2026-08-02.
# It supplied SNP's first choice, read by literal path from results/repro_phase1_1682/. That made a
# retired generation a required input of the Phase-1 rebuild, so a from-scratch reproduction could
# not be run without it. Measured on the current extraction before removing it: the UKB annotation
# named 1,958 of 1,964 instruments (99.7%), it named ALL 1,302 rows the previous table covered, and
# on every one of those 1,302 the two agreed on the rsID -- so removing the source changes no SNP
# value here. What it does remove is a hedge: the previous table's names were the ones the eight
# outcome GWAS were known to resolve, so had the two ever disagreed, its alias was the one that kept
# the instrument in the analysis. That disagreement has never occurred; if it does, the evidence is
# rsid_name_disagreements.csv and the adjudication belongs there, not in a pinned historical file.
# The UKB annotation (ukb_rsid_map_by_id.tsv, below) is now SNP's first source; the 1000G EUR .bim
# remains the second, and rsid_ld_panel still carries the .bim name for PLINK.

BIM <- "data/ld_reference/1kg_eur/EUR.bim"
if (!file.exists(BIM)) stop("missing ", BIM, " -- needed as the allele-aware rsID source")
# PINNED, not merely present. This file names variants that become SNP for every row the previous
# table does not cover, and rsid_ld_panel for every row it does -- so a different plausible .bim at
# this path silently redefines the downstream join key and the LD proxies computed from it, while
# every check in this script still passes. The manifest states which panel these bytes are and where
# they came from; the trio was verified byte-identical to the members of the source archive it names.
BIM_MANI <- "data/ld_reference/1kg_eur/EUR_MANIFEST.tsv"
if (!file.exists(BIM_MANI))
  stop("missing ", BIM_MANI, ": the LD reference panel is uncertified, so this script cannot say ",
       "which panel supplied the rsIDs and LD proxies in its output.")
# Read once, hashed and parsed from the same bytes -- the rule this script applies to every other
# certificate. Written the other way round (hash the path, then fread the path) it would carry the
# same replacement window the manifest and summary reads were just fixed for.
.bm_bytes <- readBin(BIM_MANI, "raw", file.info(BIM_MANI)$size)
BIM_MANI_MD5 <- digest::digest(.bm_bytes, algo = "md5", serialize = FALSE)
.bm <- fread(text = rawToChar(.bm_bytes), colClasses = "character")
if (!all(c("role", "file", "md5", "source") %in% names(.bm)))
  stop(basename(BIM_MANI), " is not a panel certificate (needs role/file/md5/source).")
.bm_row <- .bm[file == basename(BIM)]
if (nrow(.bm_row) != 1L)
  stop(basename(BIM_MANI), " records ", nrow(.bm_row), " row(s) for ", basename(BIM),
       ", not exactly one.")
# The .bim is 242 MB, so it is read once as text and hashed from that same text rather than being
# visited twice; the hash then describes the rows the join actually used.
.bim_txt <- readChar(BIM, file.info(BIM)$size, useBytes = TRUE)
BIM_MD5  <- digest::digest(.bim_txt, algo = "md5", serialize = FALSE)
# THE CERTIFICATE ITSELF IS PINNED, the way the UKB annotation is. Verifying the panel against a manifest that can
# be rewritten alongside it closes nothing: replacing an rsID at an existing position and updating
# the manifest md5 passes every coordinate and allele check downstream. The BIM supplies
# rsid_ld_panel for 1,762 rows and is the only source of SNP for six instruments, so a substitution
# here changes join keys and variant identity silently.
BIM_MANI_EXPECT_MD5 <- "64d30e854da52bcfe5cdb6cd0153789d"
if (!identical(tolower(BIM_MANI_MD5), BIM_MANI_EXPECT_MD5))
  stop(basename(BIM_MANI), " is not the pinned panel certificate (expected ",
       substr(BIM_MANI_EXPECT_MD5, 1, 12), ", found ", substr(BIM_MANI_MD5, 1, 12), "). The panel ",
       "it certifies decides rsid_ld_panel for most instruments; a certificate that can be rewritten ",
       "with the panel it describes certifies nothing. If the panel was legitimately rebuilt, update ",
       "this constant deliberately.")
# Hashed from the SAME text the join is built from, then compared with the certificate -- so the
# bytes that were verified and the bytes that were used cannot be two different files.
if (!identical(tolower(BIM_MD5), tolower(.bm_row$md5[1])))
  stop(basename(BIM), " is not the panel ", basename(BIM_MANI), " certifies (certificate ",
       substr(.bm_row$md5[1], 1, 12), ", on disk ", substr(BIM_MD5, 1, 12), "). It supplies rsIDs ",
       "and every LD proxy downstream, so a different panel here changes the join key silently.")
cat(sprintf("LD panel: %s certified as %s (%s)\n", basename(BIM), .bm_row$source[1],
            substr(BIM_MD5, 1, 12)))
bim <- fread(text = .bim_txt, header = FALSE,
             col.names = c("chr", "rsid_bim", "cm", "pos", "a1", "a2"),
             colClasses = list(character = 1L))
rm(.bim_txt)
bim[, key37 := paste0(chr, ":", pos)]
# The .bim is not exclusively rsIDs: 10,253 of its entries are ss-prefixed submitted-SNP accessions.
# Treating those as rsIDs made two things go wrong at once -- they became "the allele-aware answer",
# and where the UKB map supplied a real rsID for the same variant they were scored as a CONFLICT, so
# the valid rsID was thrown away and the row fell back to chr:pos for no reason (13 instruments in the
# closest raw extraction, e.g. C1R ss1388031242 vs the valid rs10849546).
bim <- bim[grepl("^rs[0-9]+$", rsid_bim)]

# The GRCh37 position and alleles carried by the UKB-PPP ID field. 08 has already asserted that this
# ID's chromosome and alleles agree with the row it sits on, so only its POSITION differs from GENPOS.
# The whole .bim join rests on ID carrying a GRCh37 position while GENPOS carries a GRCh38 one. 08
# declares that in id_build, but 08 ASSIGNS that label unconditionally, so reading it back proves
# only that a string survived a CSV round-trip -- it is not evidence about the coordinates. The
# premise is tested here against the reference panel instead, which can actually distinguish the two:
# the .bim is GRCh37, so if the premise holds, ID positions hit it far more often than GENPOS
# positions do. If the two rates were ever comparable, the join would be selecting rows by chance and
# every rsID it produced would be suspect.
if (!all(raw$id_build == "GRCh37"))
  stop("the extraction does not declare ID as GRCh37 (id_build = ",
       paste(unique(raw$id_build), collapse = "/"), "). The 1000G .bim is GRCh37, so matching ID ",
       "against it would pair instruments with the wrong variant.")
.bimpos <- unique(bim[, .(key = paste0(chr, ":", pos))])$key
.r_id   <- mean(paste0(tstrsplit(raw$ID, ":", fixed = TRUE)[[1]], ":",
                       tstrsplit(raw$ID, ":", fixed = TRUE)[[2]]) %in% .bimpos)
.r_gp   <- mean(raw$chrpos %in% .bimpos)
cat(sprintf("build check: ID positions in the GRCh37 panel %.1f%%, GENPOS positions %.1f%%\n",
            100 * .r_id, 100 * .r_gp))
if (.r_id < 0.50 || .r_id <= .r_gp * 1.5)
  stop(sprintf("ID does not behave like a GRCh37 coordinate: %.1f%% of ID positions are in the ",
               100 * .r_id),
       sprintf("GRCh37 reference panel against %.1f%% of GENPOS positions. ", 100 * .r_gp),
       "The rsID join reads ID as GRCh37; if that is wrong it pairs each instrument with whatever ",
       "variant happens to sit at the same number in the other build, silently and undetectably. ",
       "Nothing was written.")
raw[, key37 := paste0(tstrsplit(ID, ":", fixed = TRUE)[[1]], ":",
                      tstrsplit(ID, ":", fixed = TRUE)[[2]])]
# PER ROW, because the rate above is an average and an average passes with a minority of rows in the
# wrong build: at the measured 83.9% / 1.44%, a table in which only 60% of IDs are correct still
# scores 0.51 and clears the 0.50 gate. The rate says the COLUMN behaves like GRCh37; it says nothing
# about any particular row, and every row is then trusted.
#
# The dangerous row is a specific one: its ID position does NOT carry its allele pair in the panel
# while its GENPOS position DOES. Such a row would take an rsID from the panel entry sitting at its
# GRCh38 coordinate -- a different variant that happens to share the number -- and the allele-aware
# join cannot tell, because the alleles match too. Rows that simply fail to match are harmless: they
# get no rsID and are reported as such.
raw[, .akey_chk := paste(pmin(toupper(ALLELE0), toupper(ALLELE1)),
                         pmax(toupper(ALLELE0), toupper(ALLELE1)))]
.panel_va <- unique(bim[, paste0(chr, ":", pos, "|",
                                 paste(pmin(toupper(a1), toupper(a2)),
                                       pmax(toupper(a1), toupper(a2))))])
# The condition below only means something where the two coordinates DIFFER. Written without that
# restriction it read "the ID position is not in the panel AND the GENPOS position is" -- and the
# substitution it exists to catch is precisely the one where the ID has been rewritten with GENPOS,
# so key37 == chrpos and the test reduces to !X & X, false for every row. It could never fire on its
# own target.
.same_pos <- raw$key37 == raw$chrpos
.wrong <- raw[!.same_pos &
              !paste0(key37, "|", .akey_chk) %in% .panel_va &
               paste0(chrpos, "|", .akey_chk) %in% .panel_va]
if (nrow(.wrong))
  stop(nrow(.wrong), " row(s) carry an ID whose position is absent from the GRCh37 panel while ",
       "their GRCh38 position IS present with the same alleles -- so the ID on those rows behaves ",
       "like a GRCh38 coordinate, and the join would name them after whichever variant sits at that ",
       "number in the other build. The aggregate rate above cannot see this. First: ",
       paste(utils::head(.wrong[, paste0(gene, " ", ID, " (GENPOS ", chrpos, ")")], 5),
             collapse = "; "), ". Nothing was written.")
# And the case the test above structurally cannot reach: ID and GENPOS identical. The UKB-PPP ID
# carries a GRCh37 coordinate and GENPOS a GRCh38 one, so on a correctly built table they agree only
# where the two builds happen to agree for that variant -- which is common enough that this is a
# rate, not an absolute. A table in which they agree for MOST rows is one whose ID column has been
# overwritten with GENPOS, and no comparison against the panel can distinguish those rows.
.frac_same <- mean(.same_pos)
if (.frac_same > 0.90)
  stop(sprintf("%.1f%% of rows have ID and GENPOS at the same position. ", 100 * .frac_same),
       "The ID field is supposed to carry a GRCh37 coordinate and GENPOS a GRCh38 one; agreement on ",
       "almost every row means the ID column has been overwritten with the GRCh38 position, and ",
       "every rsID this script assigns would then come from whatever variant sits at that number in ",
       "the GRCh37 panel. Nothing was written.")
cat(sprintf("build check: ID and GENPOS coincide on %.1f%% of rows (both builds agree there)\n",
            100 * .frac_same))
# PER ROW for the coinciding ones too. Excluding them from .wrong and then testing only whether they
# are a MAJORITY leaves a handful of rows whose ID was replaced by GENPOS passing everything: the
# measured rates are r_id 83.9% / r_gp 1.4% / frac_same 1.3%, so one such row moves nothing.
#
# A row where the two coincide is only legitimate if the two builds really do place that variant at
# the same coordinate, which the liftover chain can decide. Rows that fail it are rows whose ID is
# not a GRCh37 coordinate at all, and they are the ones that would take an rsID from whatever the
# panel happens to hold at their GRCh38 position.
# EVERY ROW, not only the coinciding ones. The tests above are indirect: the aggregate rate cannot
# see a single row, and .wrong only fires when the ID's position is ABSENT from the panel while the
# GRCh38 position is present. The dangerous case is the opposite one -- an ID whose position exists
# in the panel carrying the same allele pair as the real variant. On the real A1BG row, ID
# 19:58858676:A:G resolves to rs1265957, while 19:69984 is another genuine A/G key in the same panel
# resolving to rs372156287; a row carrying 19:69984:A:G:imp:v1 passes the chromosome check, the
# allele check, the rate check, .wrong and the coinciding branch, and A1BG has no previous-table
# match, so SNP would simply become the wrong variant. The chain settles it directly: the ID's
# GRCh37 position must lift, uniquely and on the row's own chromosome, onto that row's GENPOS. This
# is what 08 asserts over its retained instruments; 09 does not authenticate which 08 wrote the CSV,
# so it asks the chain itself rather than inheriting the claim. ~1,700 lifts, seconds of work.
if (TRUE) {
  # A row where ID and GENPOS coincide is legitimate only if the two builds really do place that
  # variant at the same coordinate. Excluding those rows from .wrong and then asking only whether
  # they are a MAJORITY left a handful of overwritten IDs passing everything -- at the measured
  # rates (ID in panel 83.9%, GENPOS in panel 1.4%, coinciding 1.3%) a few rows move nothing.
  #
  # The chain decides it, and this script asks the chain itself rather than trusting that the
  # extractor did: 09 does not authenticate which version of 08 wrote the CSV.
  # PYTHON_BIN comes from R/00_config.R, which this script does not source -- it deliberately
  # loads only _phase1_set.R. Resolved here, with the config as the source of truth and an explicit
  # failure rather than a call to a name that does not exist.
  .pybin <- if (exists("PYTHON_BIN", envir = globalenv())) get("PYTHON_BIN", envir = globalenv())
            else {
              .cfg <- readLines(file.path("R", "00_config.R"), warn = FALSE)
              .ln <- grep("^\\s*PYTHON_BIN\\s*<-", .cfg, value = TRUE)
              if (!length(.ln))
                stop("09: PYTHON_BIN is not defined and R/00_config.R does not declare it, so the ",
                     "coinciding-position build check cannot run.")
              eval(parse(text = sub("^\\s*", "", .ln[1])))
            }
  if (!nzchar(.pybin) || !file.exists(.pybin))
    stop("09: PYTHON_BIN resolves to ", .pybin, ", which is not on disk. The coinciding-position ",
         "build check needs it.")
  .chain <- "data/ld_reference/hg19ToHg38.over.chain.gz"
  # THE CHAIN THE EXTRACTION WAS CERTIFIED WITH, not merely a chain. The summary's hash was required
  # to exist and then never compared, so a chain replaced after extraction would decide every
  # per-variant verdict here while the certificate named the other one: a malformed ID that fails
  # under the certified chain but happens to map under the replacement passes, and SNP then names
  # the wrong panel variant. The opposite mismatch is no better -- it aborts after the previous RDS
  # has already been removed above.
  if (file.exists(.chain) &&
      !identical(tolower(unname(tools::md5sum(.chain))), tolower(EXTRACT_CHAIN_MD5)))
    stop(basename(.chain), " is not the chain the extraction was certified with (summary ",
         substr(EXTRACT_CHAIN_MD5, 1, 12), ", on disk ",
         substr(unname(tools::md5sum(.chain)), 1, 12), "). Every GRCh37/GRCh38 verdict in this ",
         "script and in the extraction would then rest on two different chains.")
  if (!file.exists(.chain))
    stop("missing ", .chain, " -- needed to decide whether a row whose ID and GENPOS coincide does ",
         "so because both builds agree there, or because the ID was overwritten.")
  # Keyed on (chromosome, ID position, GENPOS): the same GRCh37 position can legitimately appear for
  # more than one row, and what is being asserted is per row -- that THIS row's ID lifts onto THIS
  # row's GENPOS.
  .cand <- unique(raw[, .(chrom  = tstrsplit(ID, ":", fixed = TRUE)[[1]],
                          pos37  = as.integer(tstrsplit(ID, ":", fixed = TRUE)[[2]]),
                          genpos = as.integer(GENPOS))])
  .inf <- tempfile(fileext = ".tsv"); .outf <- tempfile(fileext = ".tsv")
  .pyf <- tempfile(fileext = ".py")
  # DEDUPLICATED before it is written: .cand is unique on (chrom, pos37, genpos), so one GRCh37
  # position shared by two rows with different GENPOS would otherwise be lifted twice and the merge
  # below would multiply. The per-row assertion is restored by that merge, not by the input file.
  fwrite(unique(.cand[, .(chrom, pos37)]), .inf, sep = "\t", col.names = FALSE)
  # Written to a FILE, not passed as -c: the script contains quotes, tabs and newlines, and every
  # previous attempt to inline something like this through a shell mangled one of them.
  writeLines(c(
    "import sys",
    "from pyliftover import LiftOver",
    "lo = LiftOver(sys.argv[1])",
    "with open(sys.argv[2]) as fi, open(sys.argv[3], 'w') as fo:",
    "    for line in fi:",
    # DOUBLE-ESCAPED, because R resolves the escape before writeLines() ever sees it: "'\n'" puts a
    # real newline inside the Python string literal, so the generated script is
    #     c, p = line.rstrip('
    #     ').split('  ')
    # which is a SyntaxError. system2() then returns nonzero and the run stops -- after the previous
    # RDS has already been removed above. Eight instruments in the reference extraction (ACP1, AZU1,
    # BSG, CD274, CDHR5, CFD, CLUL1, COLEC12) reach this branch, all of them legitimate rows, so this
    # fires on real data rather than on a malformed table.
    "        c, p = line.rstrip('\\n').split('\\t')",
    # EVERY mapping and the chromosome it lands on, not h[0][1]. Taking the first hit accepted a
    # position the chain sends to several places -- which says nothing about the build -- and
    # discarding the returned chromosome accepted one that lands on a DIFFERENT chromosome at the
    # same number. Both are exactly the substitution this check exists to catch.
    "        h = lo.convert_coordinate('chr' + c, int(p) - 1) or []",
    "        fo.write('\\t'.join([c, p, str(len(h)),",
    "                            h[0][0][3:] if h else 'NA',",
    "                            str(h[0][1] + 1) if h else 'NA']) + '\\n')"),
    .pyf)
  .st <- suppressWarnings(system2(.pybin, shQuote(c(.pyf, .chain, .inf, .outf)),
                                  stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(.st, "status")) && attr(.st, "status") != 0)
    stop("09: the liftover of coinciding-position rows failed: ", paste(.st, collapse = " "))
  # colClasses is positional here: col.names is applied AFTER the columns are typed, so naming the
  # column in colClasses matched nothing and the chromosome came back as an integer -- which loses
  # X and Y.
  .lf <- fread(.outf, header = FALSE,
               colClasses = c("character", "integer", "integer", "character", "integer"),
               col.names = c("chrom", "pos37", "n_hit", "chrom38", "pos38"))
  # UNIQUE, ON THE SAME CHROMOSOME, AT THE SAME COORDINATE. A position the chain sends to several
  # places establishes nothing, and one that lands on another chromosome at the same number is the
  # substitution being looked for, not evidence against it.
  .lf <- merge(.lf, .cand, by = c("chrom", "pos37"), all.x = TRUE, sort = FALSE)
  .bad_same <- .lf[is.na(pos38) | n_hit != 1L | is.na(chrom38) | chrom38 != chrom | pos38 != genpos]
  if (nrow(.bad_same))
    stop(nrow(.bad_same), " instrument position(s) whose ID does not lift, uniquely and on its own ",
         "chromosome, onto that row's GENPOS -- so the ID is not the GRCh37 coordinate of the ",
         "variant the row describes, and the panel join would name the row after whatever variant ",
         "sits at the ID's number instead. First: ",
         paste(utils::head(.bad_same[, paste0(chrom, ":", pos37, " -> ",
                                              fifelse(is.na(pos38), "unmappable",
                                                      paste0(chrom38, ":", pos38, " (", n_hit,
                                                             " mapping(s))")))], 5),
               collapse = "; "), ". Nothing was written.")
  cat(sprintf("build check: all %d instrument position(s) lift from their ID onto their own GENPOS\n",
              nrow(.lf)))
}
cat(sprintf("build check: 0 of %d row(s) have an ID that behaves like GRCh38\n", nrow(raw)))
raw[, .akey_chk := NULL]
# Collapsed to ONE row per (position, unordered allele pair) BEFORE the join, so the join is 1:1 by
# construction. The previous order -- join on position, multiply, then resolve -- meant a position
# with several .bim entries produced extra rows and tripped the row-multiplication check below before
# the ambiguity rule could apply, making that rule unreachable. Where a single (position, allele
# pair) carries more than one rsID, neither is usable and both are dropped here.
bim[, akey := paste(pmin(toupper(a1), toupper(a2)), pmax(toupper(a1), toupper(a2)))]
# SNVs ONLY. There is no sound ordered key against a .bim: its a1/a2 carry no REF/ALT semantics, so
# no reconstruction of the UKB orientation is possible and an unordered indel match is exactly the
# collapse this key was just fixed to avoid -- with no second source able to adjudicate it. A
# REF/ALT-aware source (a normalised dbSNP VCF keyed on build/chrom/pos/REF/ALT) would be needed to
# name an indel from the panel, and that is a different input, not a smaller change. Refusing costs
# nothing measurable here: the panel matched 0 of the 61 indel rows in the reference extraction.
bl <- unique(bim[grepl("^[ACGT]$", toupper(a1)) & grepl("^[ACGT]$", toupper(a2)),
                 .(key37, akey, rsid_bim)])
# `ambig` is kept, not just the NA: a key the panel refused because it names the same variant more
# than once is a different situation from a key the panel does not hold at all, and the loss sidecar
# below reported both as "the panel has this position but not this allele pair".
bl <- bl[, .(rsid_bim = if (uniqueN(rsid_bim) == 1L) rsid_bim[1] else NA_character_,
             ambig = uniqueN(rsid_bim) > 1L),
         by = .(key37, akey)]
raw[, akey := paste(pmin(toupper(ALLELE0), toupper(ALLELE1)),
                    pmax(toupper(ALLELE0), toupper(ALLELE1)))]
raw <- merge(raw, bl, by = c("key37", "akey"), all.x = TRUE, sort = FALSE)
raw[, akey_ambig := !is.na(ambig) & ambig]
raw[, ambig := NULL]
# The unordered allele key above IS the allele test: a row only matches a .bim entry whose two
# alleles are the same two, in either orientation. Rejected-for-alleles is therefore indistinguishable
# from absent here, which the loss sidecar below accounts for explicitly.

# 08 guarantees one row per (gene, variant), so this join must not have multiplied anything; if it
# did, the later nrow(out) == nrow(raw) check could not see it, because it compares against the
# deduplicated object.
.n_before <- nrow(raw)
raw <- unique(raw, by = c("chrpos", "ensg"))
if (nrow(raw) != .n_before)
  stop(sprintf("the .bim join multiplied %d certified row(s); a position matched several bim entries ",
               .n_before - nrow(raw)),
       "whose alleles both matched, so which rsID is kept would depend on row order.")

# The previous instrument table's rsID join stood here and was REMOVED 2026-08-02 (see the header):
# a reproduction run must not need results/repro_phase1_1682/. `akey` was set for the .bim join above
# and is dropped here, where the removed join used to drop it.
raw[, akey := NULL]

# ── THE AUTHORITATIVE UKB ANNOTATION, joined on the FULL variant ID ──────────────────────────────
# THIS IS A THIRD rsID SOURCE, and the file's stated contract now has three, not two. It is declared
# here rather than left as an undocumented addition because it is not decoration: it names SNP for
# 1,958 rows and is the ONLY source for 66 of them (with it, 89 instruments and 88 proteins are
# removed for lack of an rsID; without it, 155 and 154). Unlike the UKB annotation, which is pinned to an expected
# md5, this map has no pinned content hash -- it is pinned instead to the EXTRACTION it was built
# for, via its sidecar's instruments_md5, and its own bytes are hashed once and recorded in the RDS.
# Whether it should outrank the pinned table and the panel is an analysis decision; what is fixed
# here is that its identity and its effect are both recorded rather than implicit.
# The BIM is indirect: it is matched through the GRCh37 coordinate embedded in the UKB-PPP ID, and
# cannot name a variant this rebuild found and it did not have, and the other is a reference panel
# that holds only what it holds. UKB-PPP publishes its own per-chromosome rsID annotation keyed on
# exactly the identifier this extraction carries (chr:pos:ref:alt:imp:v1), so it needs no coordinate
# projection and no allele key -- the ID already fixes position, REF and ALT, which is what made the
# earlier position-only map unusable at multiallelic sites and reciprocal indel pairs.
#
# It is joined FIRST because it is the authority: where UKB names a variant, that is its name.
#
# REQUIRED, NOT OPTIONAL. It used to be used only if the file happened to exist, which made the
# instrument universe depend on the state of the disk rather than on the data: with the map the RDS
# holds 1,964 instruments over 1,948 proteins, without it 1,898 over 1,882 -- the same certified
# extraction giving two answers 66 proteins apart. A source that decides 66 proteins is either part
# of the contract or it is not; it cannot be a bonus that appears when a file happens to be present.
#
# ...AND PINNED to an expected hash, like the .bim. The sidecar alone is self-referential: it supplies
# both the expected map hash and the claim about which extraction the map belongs to, so replacing
# the map and rewriting its sidecar passes every check -- while this map decides SNP for 1,958 rows
# and is the sole evidence for 66. A pinned constant is the one thing a rewritten sidecar cannot
# satisfy. If the map is legitimately rebuilt, update this constant deliberately.
.UKBMAP <- "data/ukb_ppp/ukb_rsid_map_by_id.tsv"
.UKBMAP_EXPECT_MD5 <- "e300734eb695938ef360ec69cec6dcbf"
{
  if (!file.exists(.UKBMAP))
    stop(basename(.UKBMAP), " is absent. It is a REQUIRED rsID source, not an optional one: without ",
         "it this extraction yields 66 fewer proteins, and which of those two answers you get would ",
         "depend only on whether a file is on disk. Run scripts/_shared/phase1/04_download_rsid_map.py.")
  .um_prov <- paste0(.UKBMAP, ".provenance.tsv")
  if (!file.exists(.um_prov))
    stop(basename(.UKBMAP), " has no provenance sidecar, so it cannot be shown to have been built ",
         "for this instrument set. Rerun scripts/_shared/phase1/04_download_rsid_map.py.")
  .um_p <- fread(.um_prov, sep = "\t")
  .um_want <- unname(tools::md5sum(IN))
  .um_got  <- .um_p[field == "instruments_md5", value][1]
  if (!identical(tolower(.um_got %||% ""), tolower(.um_want)))
    stop(sprintf(paste0("%s was built from a different extraction than this run reads (sidecar %s, ",
                        "%s is %s). Its rsIDs would name variants from another generation. Rerun ",
                        "04_download_rsid_map.py against %s."),
                 basename(.UKBMAP), substr(.um_got %||% "<none>", 1, 12), basename(IN),
                 substr(.um_want, 1, 12), IN))
  # ...and the sidecar must describe THIS table. 04 installs the two files with separate replaces, so
  # a crash or a concurrent run can leave a sidecar beside a table it does not describe -- detectable
  # only because the sidecar carries the table's hash, and only if someone compares them.
  # ONE READ, and everything below is about THOSE bytes. Hashing the path, then reopening it to
  # parse, then hashing the path again for the RDS attribute leaves two windows in which the file can
  # change -- and it DID change during a review of this file, from f535d869 (328,565 bytes) to
  # e300734e (340,259 bytes), because 04 was being re-run at the same time. The two versions happened
  # to carry the same 1,949 identifiers; nothing in the code required that. rsid_ukb has first
  # priority over both stated sources, so a swap here changes variant identity in the shipped RDS.
  .um_bytes <- readBin(.UKBMAP, "raw", file.info(.UKBMAP)$size)
  .UKBMAP_MD5 <- digest::digest(.um_bytes, algo = "md5", serialize = FALSE)
  if (!identical(tolower(.UKBMAP_MD5), .UKBMAP_EXPECT_MD5))
    stop(sprintf(paste0("%s is not the pinned map (expected %s, on disk %s). Its own sidecar cannot ",
                        "settle this, because the sidecar is rewritten together with the map. This ",
                        "file decides SNP for most instruments and is the only evidence for 66 of ",
                        "them, so a changed map is a deliberate act -- update .UKBMAP_EXPECT_MD5 if ",
                        "it was one."),
                 basename(.UKBMAP), substr(.UKBMAP_EXPECT_MD5, 1, 12), substr(.UKBMAP_MD5, 1, 12)))
  .um_omd5 <- .um_p[field == "output_md5", value][1]
  if (!identical(tolower(.um_omd5 %||% ""), tolower(.UKBMAP_MD5)))
    stop(sprintf(paste0("%s is not the table its provenance sidecar describes (sidecar output_md5 ",
                        "%s, bytes read %s). The two are installed by separate writes and one of them ",
                        "did not complete, or the file changed under this run. Rerun ",
                        "04_download_rsid_map.py."),
                 basename(.UKBMAP), substr(.um_omd5 %||% "<none>", 1, 12), substr(.UKBMAP_MD5, 1, 12)))
  .um <- fread(text = rawToChar(.um_bytes), sep = "\t", select = c("ukb_variant_id", "rsid"))
  setnames(.um, c("ID", "rsid_ukb"))
  # ONE rsID PER ID, or refuse. unique(by="ID") keeps whichever row came first, so an ID annotated
  # with two different rsIDs would be resolved by file order and an arbitrary name would become the
  # join key this instrument is searched for in every outcome.
  .um <- .um[grepl("^rs[0-9]+$", rsid_ukb)]
  .um_c <- .um[, .N, by = .(ID, rsid_ukb)][, .N, by = ID][N > 1L]
  if (nrow(.um_c))
    stop(nrow(.um_c), " variant ID(s) in ", basename(.UKBMAP), " carry more than one rsID (e.g. ",
         paste(utils::head(.um_c$ID, 3), collapse = ", "), "). Row order would decide which name ",
         "this instrument is searched for.")
  .um <- unique(.um, by = "ID")
  raw <- merge(raw, .um, by = "ID", all.x = TRUE, sort = FALSE)
  .UKBMAP_USED <- .UKBMAP; .UKBMAP_N <- sum(!is.na(raw$rsid_ukb))
  cat(sprintf("  UKB annotation: %d of %d instrument row(s) named (%.1f%%)\n",
              .UKBMAP_N, nrow(raw), 100 * mean(!is.na(raw$rsid_ukb))))
}
# Both sources are now allele-aware, and that changes what a disagreement MEANS. When the map was
# keyed on position alone, the two could be describing different variants, and discarding both was
# the only safe answer. Now they have each matched the same position AND the same allele pair, so a
# disagreement is two NAMES for one variant -- merged or retired rsIDs, most often -- not two
# variants. Discarding both would drop real proteins over a naming difference: on the closest raw
# table seven single-instrument proteins (ADM, CD84, CEP170, FST, GUSB, MIA, PGD) hit this rule, and
# because seven is under the loss gate the script would have succeeded while silently removing them.
#
# BOTH names are kept, because SNP is asked to serve two identifier spaces that do not agree:
#
#   - 03_build_ld_proxies.R:32 selects the panel key and :58 hands it to PLINK --ld-snp-list against
#     which knows the .bim name.
#   - 02_extract_outcomes_local.R:26 builds its lookup from inst$SNP and keeps exact string matches
#     only, against outcome files that carry the UKB annotation's name.
#
# One column cannot satisfy both, and the earlier rule -- ".bim wins, it is the identifier space
# every downstream LD operation works in" -- decided it for the LD side and silently cost the other.
# Checked on disk: the seven conflicts (ADM, CD84, CEP170, FST, GUSB, MIA, PGD) appear in FinnGen
# T2D, Wuttke CKD, Wuttke eGFR and the harmonised MEGASTROKE under the PREVIOUS name, and under the
# panel name in none of them. With the panel name in SNP those variants read as absent from the
# outcome and are proxy-substituted or dropped, although the exact variant is present.
#
# So SNP takes the outcome-resolvable name and rsid_ld_panel carries the panel's, and the choice is
# asymmetric for a reason: the proxy machinery exists to rescue variants the OUTCOME lacks, so a
# name that resolves in the outcome removes the need for a proxy, whereas a name that resolves only
# in the panel produces a proxy for a variant that was there all along.
# ANY disagreement between the panel and the UKB annotation. An earlier revision compared the panel
# only against the previous instrument table, so a UKB/BIM disagreement that table had no opinion
# about -- LRRC37A2 -- never reached the sidecar the disagreements are itemised in.
# EVERY disagreement between the two remaining sources, not just those involving the panel. The
# previous form required rsid_bim to be present, so a UKB-versus-pinned disagreement at a variant the
# panel does not hold went unrecorded -- and since rsid_ukb wins, that is precisely the case where
# the name the OUTCOME files resolve (the pinned one) would be silently replaced and the instrument
# could vanish from outcome extraction. Zero such rows exist in the current data; the point is that
# nothing was looking.
raw[, conflict := (!is.na(rsid_bim) & !is.na(rsid_ukb) & rsid_bim != rsid_ukb)]
# THE UKB ANNOTATION FIRST. It is the authority on what a variant IS, keyed on the full variant ID
# rather than a coordinate, and it names 99.7% of the instruments on its own. The previous
# instrument table used to take precedence over it; that source was removed (see the block above),
# and on the data it covered the two never disagreed.
raw[, SNP := fifelse(!is.na(rsid_ukb) & nzchar(rsid_ukb), rsid_ukb,
             fifelse(!is.na(rsid_bim), rsid_bim, chrpos))]
# Never NA and never empty: SNP is a join key in two downstream stages, and an empty string there
# matches nothing while looking like a value.
if (anyNA(raw$SNP) || !all(nzchar(raw$SNP)))
  stop(sum(is.na(raw$SNP) | !nzchar(raw$SNP)), " row(s) have no SNP identifier at all.")
raw[, rsid_ld_panel := fifelse(!is.na(rsid_bim), rsid_bim, NA_character_)]
if (any(raw$conflict)) {
  # rsid_ukb IS one of the columns. It decides SNP whenever it is present, so omitting it from the
  # disagreement listing left the reader unable to see which name won or why.
  fwrite(raw[conflict == TRUE, .(gene, chrpos, ALLELE0, ALLELE1, rsid_ukb, rsid_bim,
                                 used_for_outcome = SNP, used_for_ld = rsid_ld_panel)],
         file.path(OUTD, "rsid_name_disagreements.csv"))
  cat(sprintf("  %d variant(s) are named differently by the reference panel and the UKB annotation; ",
              sum(raw$conflict)),
      "SNP carries the UKB annotation's name and ",
      "rsid_ld_panel carries the panel's. Every case is listed in rsid_name_disagreements.csv. ",
      "03_build_ld_proxies.R resolves the panel name through ld_panel_key() (03:27), so these rows ",
      "are already looked up under rsid_ld_panel rather than under SNP.\n", sep = "")
}

n_rs   <- sum(grepl("^rs[0-9]+$", raw$SNP))
# THE DENOMINATOR THE COVERAGE WAS MEASURED ON, captured before the unmapped rows are
# removed below. Reported against nrow(out) it is necessarily 100% whenever an RDS is written
# at all -- 619/663 became "100%" -- so the certificate and the closing message concealed
# exactly the loss they exist to report.
n_rows_certified <- nrow(raw)
n_conf <- sum(raw$conflict)
# WHICH SOURCE ACTUALLY SUPPLIED THE NAME, counted the way SNP is assigned -- the branch that WON,
# not a membership test, so the two source counts sum to n_rs.
cat(sprintf(paste0("rsID: %d / %d (%.1f%%) -- UKB annotation %d, LD panel %d; ",
                   "named differently by two sources %d; none %d\n"),
            n_rs, nrow(raw), 100 * n_rs / nrow(raw),
            sum(!is.na(raw$rsid_ukb) & nzchar(raw$rsid_ukb)),
            sum((is.na(raw$rsid_ukb) | !nzchar(raw$rsid_ukb)) & !is.na(raw$rsid_bim)),
            # NOT nrow - n_rs - n_conf: a conflict row keeps the panel's rsID, so it is already
            # counted in n_rs. Subtracting it again understated the genuinely unmapped rows by
            # exactly the number of conflicts.
            n_conf, nrow(raw) - n_rs))

# Defined unconditionally. It was assigned only inside `if (nrow(lost))`, and the manifest below
# always reads it -- so a run in which EVERY instrument got an rsID died with "object 'gone' not
# found", after the RDS had already been installed and before its certificate was written.
gone <- character(0)
lost <- raw[!grepl("^rs[0-9]+$", SNP)]
if (nrow(lost)) {
  fwrite(lost[, .(gene, chrpos, ID, ALLELE0, ALLELE1, oid,
                  # Distinguishes "nothing was found" from "something was found and refused", which
                  # the single previous label conflated -- the sidecar exists to show why a protein
                  # left the analysis, and those are different answers.
                  # Both source columns hold either a valid rsID or NA, so a lost row always has
                  # both missing -- "found but not usable" was unreachable. What IS distinguishable
                  # is WHY the panel gave nothing: the exact allele pair carried several rsIDs and
                  # was therefore refused, the position exists with different alleles, or neither.
                  reason = fifelse(akey_ambig,
                                   "the reference panel names this exact variant more than once",
                            fifelse(key37 %in% bl$key37,
                                    "the reference panel has this position but not this allele pair",
                                    "no rsID in either source")))][order(gene)],
         file.path(OUTD, "instruments_without_rsid.csv"))
  # ENSG, not the symbol. 08 keys proteins by ENSG because symbols collide and get renamed; counting
  # losses by symbol would let one protein's surviving instrument hide another's loss when two ENSGs
  # share a name. No collision exists in the current certified table, so this is a latent difference
  # rather than a changed number today.
  gone <- setdiff(unique(raw$ensg), unique(raw[grepl("^rs[0-9]+$", SNP)]$ensg))
  cat(sprintf("  %d instrument(s) have no usable rsID; %d gene(s) lose EVERY instrument as a result\n",
              nrow(lost), length(gone)))
  cat("  listed in ", file.path(OUTD, "instruments_without_rsid.csv"), "\n", sep = "")
  if (length(gone))
    cat("  genes lost entirely: ", paste(head(sort(gone), 15), collapse = ", "),
        if (length(gone) > 15) " ..." else "", "\n", sep = "")
  # HOW MANY VARIANTS HAVE AN rsID IS A PROPERTY OF THE DATA, NOT A QUALITY BAR TO PASS.
  #
  # This block used to abort unless >=90% of instruments carried an rsID and <=2% of proteins lost
  # all of theirs. Both numbers were invented here -- no reviewer, editor or reference set them --
  # and a threshold on a result is the wrong instrument entirely: it cannot make a variant that
  # dbSNP does not name become named, it can only refuse to write the table. Worse, the obvious way
  # to make a run pass is to move the number, which is choosing a criterion to fit the result.
  #
  # The 2% gate also measured the wrong thing. It asked how many of the genes THIS extraction found
  # are unusable, so a rebuild that finds MORE proteins is penalised for the extra ones it cannot
  # name: 08 found 2,036 genes against the previous generation's 1,682, of which 1,948 end up
  # analysable -- a net gain of 266 proteins -- and the gate reported that as a 4.3% failure. The
  # 88 it counted as lost were never in the previous set at all.
  #
  # So the loss is REPORTED, not gated: every instrument without a usable rsID is itemised in
  # instruments_without_rsid.csv, the genes that lose all of theirs are named in the log above, and
  # the counts go into the certificate. A reader can see exactly what is missing and decide what it
  # means; this script does not decide for them.
  #
  # What IS still refused is a broken join, which is not a data property but a failure of this
  # script: if NOT ONE instrument received an rsID, the map did not load, the key did not match, or
  # the merge silently produced nothing. Zero is not a tunable threshold.
  if (n_rs == 0L)
    stop("not one of ", nrow(raw), " instruments received an rsID. That is not a property of the ",
         "data -- it means the rsID sources did not load or did not join. Check the panel .bim ",
         "and ", basename(.UKBMAP), " before reading anything into this run.")
}

# ...and REMOVED, so the table this script ships IS the set that can be analysed. They used to stay
# in the RDS with SNP = chrpos, and 02_extract_outcomes_local.R:26 then discarded every non-rs row --
# so the RDS asserted proteins into the instrument universe that no outcome was ever searched for, a
# multi-instrument protein could silently be analysed on fewer variants than the table held, and any
# protein count taken from the RDS exceeded the analysable count. The gates above still measure the
# loss on the FULL table, it is reported in the log and itemised in instruments_without_rsid.csv, and
# the certificate below records it: what changes here is only that the shipped table stops
# disagreeing with what downstream actually uses.
.n_unmapped <- nrow(lost)
raw <- raw[grepl("^rs[0-9]+$", SNP)]
if (.n_unmapped)
  cat(sprintf("  removed %d instrument(s) with no usable rsID; the table now holds %d row(s) over %d gene(s)\n",
              .n_unmapped, nrow(raw), uniqueN(raw$ensg)))
if (!nrow(raw))
  stop("every instrument was removed for lack of an rsID; there is no exposure set to write.")

# ---- 3. the harmonised contract -------------------------------------------------------------------
# r2 for a single SNP: F/(F+N-2). Kept identical to phase1/12_phase1_r2.R so the column means the
# same thing it did in the submitted analysis.
# The shipped analysis input has exactly 20 columns and `exposure` is "GENE (ACCESSION)", not the bare
# symbol -- id.exposure is the symbol. Both were verified against
# the previous Phase 1 set (R/_phase1_set.R). Three further columns are added deliberately and
# named here so the schema change is declared rather than discovered: ensg (a stable key that does not
# move when a symbol is renamed), assignment_basis (which assays entered under a declared convention),
# and ukb_variant_id (UKB-PPP's own name for the variant, chr:pos:ref:alt:imp:v1). It was called
# id_grch37, which is accurate about the position embedded in the string and badly misleading about
# the table: every coordinate in this file is GRCh38, and a column announcing GRCh37 makes a reader
# stop and check which build the table is in -- in an analysis whose original defect was exactly a
# build mix-up. It is an IDENTIFIER, not a coordinate: never use it as a position in either build.
# and pval_underflow. That is FOUR added columns,
# not three; an earlier version of this comment said three and had already added the fourth.
# The twenty original columns keep their original names AND their original order, so a positional
# consumer of the old table still reads the same fields; the four additions come after them. That
# invariant is asserted below rather than left to the ordering of this expression.
# THE SHIPPED ORDER, verbatim. A previous version claimed the original twenty keep their names and
# their order with the additions appended, and then wrote a different order entirely -- beta and se
# before the alleles, oid and uniprot swapped, ensg inserted in the middle. Every positional consumer
# of the old table would have read the wrong field, silently. The order is asserted below against a
# literal copy of the shipped header, so this cannot drift again.
out <- raw[, .(
  SNP,
  chrpos,
  effect_allele = ALLELE1,
  other_allele  = ALLELE0,
  beta          = BETA,
  se            = SE,
  eaf           = A1FREQ,
  pval,
  samplesize    = N,
  F_stat,
  log10p        = LOG10P,
  gene,
  uniprot,
  oid,
  gene_chrom, gene_start, gene_end,
  r2_exposure   = F_stat / (F_stat + N - 2),
  exposure      = paste0(gene, " (", uniprot, ")"),
  id.exposure   = gene,
  # ---- the four additions, AFTER the original twenty ----
  ensg,
  assignment_basis,          # carried so a declared-convention assay stays identifiable downstream
  ukb_variant_id = ID,
  # 08 floors p-values that underflow binary64 (495 of 1,723 rows in the shipped reference are
  # affected). Without this flag a floored p-value is indistinguishable from an exact one, and log10p
  # is the only place the information survives. Anything ranking instruments must use log10p.
  pval_underflow,
  # The reference panel's name for this variant, when it differs from the one in SNP. SNP has to be
  # the name the OUTCOME files resolve; the LD panel knows a different one for seven of these
  # variants, and PLINK is given SNP, so the panel's name has to travel too or those rows simply
  # cannot be looked up in the panel. NA only when the panel does not name the variant at all: it is
  # filled for EVERY .bim match, including the ones where the two names agree. (The comment used to
  # say "NA when the two agree", which the assignment a few lines below has never done -- and
  # ld_panel_key() in R/_ld_namespace.R relies on it being populated for the agreeing rows too.)
  rsid_ld_panel,
  alt_locus_ensg, alt_locus_start, alt_locus_end,
  # APPENDED, never inserted: the shipped twenty keep their positions, so a positional consumer of
  # the old table still reads the same fields. These two are the interval the instruments were
  # actually extracted from, and every regional query downstream has to use them instead of
  # recomputing gene +/- CIS_WINDOW -- which loses 08's chromosome-end truncation.
  cis_start, cis_end
)]

# The shipped 20-column contract, spelled out rather than described. results/repro_phase1_1682's RDS
# has exactly these names in exactly this order.
SHIPPED_20 <- c("SNP", "chrpos", "effect_allele", "other_allele", "beta", "se", "eaf", "pval",
                "samplesize", "F_stat", "log10p", "gene", "uniprot", "oid", "gene_chrom",
                "gene_start", "gene_end", "r2_exposure", "exposure", "id.exposure")
# The additions, in shipped order. The name no longer states a count: it said "four" while listing
# eight, and a name that has to be revised every time the list grows is a comment that will be wrong
# again. The three alt_locus_* columns are here because 08 carries them specifically so a
# declared-convention assignment stays reconstructible -- keeping only assignment_basis said that an
# alternative exists without saying what it is, which is the half the outstanding alternative-window
# sensitivity analysis actually needs. cis_start/cis_end are the extracted interval (see above).
ADDED <- c("ensg", "assignment_basis", "ukb_variant_id", "pval_underflow", "rsid_ld_panel",
           "alt_locus_ensg", "alt_locus_start", "alt_locus_end", "cis_start", "cis_end")
# id.exposure IS the gene symbol, while uniqueness is checked on (ensg, SNP), so two loci sharing a
# symbol pass every check and are then grouped as one exposure by everything downstream -- combining
# instruments from two different genes. Zero collisions in the current 2,915-row coordinate table;
# asserted rather than assumed, because the check that would catch it does not look at this.
.gmap <- unique(out[, .(gene, ensg)])
if (anyDuplicated(.gmap$gene))
  stop(sum(duplicated(.gmap$gene)), " gene symbol(s) stand for more than one ENSG (",
       paste(utils::head(unique(.gmap$gene[duplicated(.gmap$gene)]), 3), collapse = ", "),
       "). id.exposure is the symbol, so downstream code would pool instruments from two loci into ",
       "one exposure.")
if (!identical(names(out), c(SHIPPED_20, ADDED)))
  stop("the output schema does not match the contract.\n  expected: ",
       paste(c(SHIPPED_20, ADDED), collapse = ", "), "\n  got:      ",
       paste(names(out), collapse = ", "))

stopifnot(
  # .declared MINUS the rows removed for lack of an rsID, not nrow(raw): `out` is a column projection
  # of `raw`, so comparing the two could never fail. The certified count is still the only external
  # number to check against, and the removal is the one accounted-for difference from it.
  nrow(out) == .declared - .n_unmapped,
  !anyNA(out$SNP), all(nzchar(out$SNP)),
  !anyNA(out$beta), !anyNA(out$se), all(out$se > 0),
  !anyNA(out$eaf), all(out$eaf > 0 & out$eaf < 1),
  !anyNA(out$samplesize), all(out$samplesize > 0),
  all(out$r2_exposure > 0 & out$r2_exposure < 1),
  # THE REST OF THE SHIPPED COLUMNS, which were never required to be present. A mapped row with a
  # missing pval, allele, coordinate or underflow flag reached the RDS and was dropped much later by
  # formatting or harmonisation -- silently, and by then attributable to neither this script nor the
  # extraction. Every one of these is read downstream.
  # F_stat IS (BETA/SE)^2, recomputed rather than trusted: r2_exposure is derived straight from the
  # supplied value, so a positive but wrong F_stat put a wrong F_stat AND a wrong r2 into the RDS
  # while satisfying every range test. The tolerance is relative because the extraction writes the
  # value as float64 text.
  all(abs(out$F_stat - (out$beta / out$se)^2) <= 1e-6 * pmax(1, (out$beta / out$se)^2)),
  !anyNA(out$pval), !anyNA(out$log10p), !anyNA(out$pval_underflow),
  !anyNA(out$effect_allele), all(nzchar(out$effect_allele)),
  !anyNA(out$other_allele),  all(nzchar(out$other_allele)),
  !anyNA(out$chrpos),        all(nzchar(out$chrpos)),
  !anyNA(out$gene_chrom),    all(nzchar(out$gene_chrom)),
  !anyNA(out$gene_start), !anyNA(out$gene_end),
  !anyNA(out$cis_start),  !anyNA(out$cis_end),
  all(out$assignment_basis %in% c("hgnc_route", "declared_convention")),
  # SNP is the downstream join key, so a duplicate would make the join ambiguous
  !any(duplicated(out[, .(ensg, SNP)])),
  # ... and uniqueness WITHIN a gene is not enough. The same variant legitimately instruments several
  # proteins (ten rsIDs do), but every such row must carry the same coordinate and the same alleles.
  # An rsID meaning two different variants in two genes is exactly what a position-only rsID source
  # can produce, and 02_extract_outcomes_local.R keeps the first allele pair it sees for each SNP --
  # so the second protein would be harmonised against the wrong variant with nothing to show for it.
  out[grepl("^rs", SNP), uniqueN(paste(chrpos, effect_allele, other_allele)), by = SNP][, all(V1 == 1L)]
)

n_conv <- uniqueN(out[assignment_basis == "declared_convention"]$gene)
if (n_conv)
  cat(sprintf("[!] %d gene(s) enter the analysis under a declared convention rather than a settled\n",
              n_conv),
      "    assignment; assignment_basis carries that through this table too.\n", sep = "")

# Re-checked at COMMIT time. Reading the CSV once and hashing those bytes stops this run from
# parsing a half-replaced file, but it says nothing about whether that extraction is still the
# certified one when the RDS is installed: 08 can start a new extraction meanwhile, and this script
# would then leave an RDS derived from the previous one beside a newly certified CSV.
if (!file.exists(MANI))
  stop("the extraction certificate disappeared while this script was running: 08 is extracting ",
       "again. Nothing was installed; rerun 09 when that extraction has certified.")
.mani_now <- fread(MANI, colClasses = "character")
.row_now  <- .mani_now[file == basename(IN)]
if (nrow(.row_now) != 1L || !identical(tolower(.row_now$md5[1]), tolower(row$md5[1])))
  stop("the certified extraction changed while this script was running (was ",
       substr(row$md5[1], 1, 12), ", now ",
       if (nrow(.row_now) == 1L) substr(.row_now$md5[1], 1, 12) else "absent",
       "). The table this run harmonised is no longer the current one; nothing was installed.")

# The lock was taken at the top of this script and is still held, so 08 cannot have finalised while
# this run was computing. The re-read above remains as the check that the input is still the one
# that was harmonised; together they close the window rather than narrowing it.

tmp <- sprintf("%s.%d.tmp", OUT, Sys.getpid())
# The table carries its OWN provenance, as attributes, in addition to the sidecar certificate.
# Installing the table and its certificate takes two renames, and no in-process handler can cover a
# hard kill between them -- the canonical path would then hold an uncertified table that every
# consumer reads normally, because none of them inspects the sidecar. A table that certifies itself
# is checkable no matter how the run ended, and by anything that loads it.
# The bytes `raw` was PARSED FROM, not a fresh read of the path. Re-reading records the hash of
# whatever is at that path now, so a CSV replaced after the read leaves the object holding results
# from the original bytes and an attribute naming the replacement.
attr(out, "ckm_source_md5")  <- IN_MD5
attr(out, "ckm_source_file") <- basename(IN)
attr(out, "ckm_builder")     <- "09_phase1_harmonise_hg38.R"
attr(out, "ckm_built")       <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
# The build of each coordinate space, IN the table. gene_start/gene_end kept their names when they
# changed from GRCh37 to GRCh38, so nothing in the object itself contradicts a consumer that applies
# the old contract and lifts them a second time.
# The two inputs that DECIDE SNP. Without them the object identifies the association data but not
# either rsID source, and changing either changes every downstream join key. They lived only in the
# sidecar, which a kill between the two renames leaves absent while the RDS is already canonical.
attr(out, "ckm_ld_panel")        <- BIM
attr(out, "ckm_ld_panel_md5")    <- BIM_MD5
# ...AND THE THIRD ONE, which is the one that decides most of them. The UKB annotation names SNP for
# 1,958 rows and is the ONLY source for 66 of them: with it, 89 instruments and 88 proteins are
# removed for lack of an rsID; without it, 155 and 154. Whether that file is present therefore
# changes the instrument universe by 66 proteins -- and it was recorded nowhere, so two RDS files
# built from identical inputs but a different answer to "was the map on disk" were indistinguishable.
attr(out, "ckm_ukb_rsid_map")     <- if (nzchar(.UKBMAP_USED)) .UKBMAP_USED else NA_character_
# The hash of the BYTES THAT SUPPLIED THE rsIDs, carried down from the single read above -- not a
# fresh md5sum() of the path, which would describe whatever is there now rather than what was used.
attr(out, "ckm_ukb_rsid_map_md5") <- if (nzchar(.UKBMAP_USED)) .UKBMAP_MD5 else NA_character_
attr(out, "ckm_ukb_rsid_named")   <- .UKBMAP_N
attr(out, "ckm_ld_panel_manifest_md5") <- BIM_MANI_MD5
# THE THREE CERTIFICATES THIS TABLE RESTS ON, in the object rather than only in the sidecar written
# afterwards. A kill between the RDS rename and the manifest rename leaves a complete table at the
# canonical path with no certificate beside it, and nothing downstream inspects the sidecar -- so
# without these the object could not say which coordinate table, which extraction summary or which
# completion manifest produced it. The run signature is opaque: it cannot reveal the coordinate hash
# it contains without the summary, which is exactly the file that may be missing.
attr(out, "ckm_coords")               <- COORDS_TBL
attr(out, "ckm_coords_md5")           <- EXTRACT_COORDS_MD5
attr(out, "ckm_extract_summary_md5")  <- SUMMARY_MD5
attr(out, "ckm_extract_manifest_md5") <- MANI_MD5
# How many certified instruments this table does NOT contain, and why. Without it the row count
# cannot be reconciled with the extraction's own certified n_snps_total.
attr(out, "ckm_unmapped_removed")     <- .n_unmapped
# The fields that decide how this table may be read, carried IN it. Installing the RDS and its
# sidecar takes two renames; a kill in between leaves the canonical table with no sidecar, and the
# downstream readers load the table without requiring one. What matters there is not every hash but
# whether a gate was bypassed and what was lost.
attr(out, "ckm_builder_md5")      <- SELF_MD5
attr(out, "ckm_extract_run_sig")  <- EXTRACT_RUN_SIG
attr(out, "ckm_extract_chain_md5") <- EXTRACT_CHAIN_MD5   # the chain every 37<->38 step here used
attr(out, "ckm_extractor_md5")    <- EXTRACTOR_MD5
attr(out, "ckm_proteins_lost")    <- sort(gone)
attr(out, "ckm_gene_coord_build") <- "GRCh38"
attr(out, "ckm_chrpos_build")     <- "GRCh38"
attr(out, "ckm_id_build")         <- "GRCh37"
# The extracted interval is a THIRD coordinate space in this table, and the one that regional queries
# are built from. It cannot catch a GRCh38 window handed to a GRCh37 resource -- that has to be
# asserted where the query is made -- but it does let a consumer state which build it is holding.
attr(out, "ckm_cis_coord_build")  <- "GRCh38"
saveRDS(out, tmp)
# Compares CONTENT, not just dim: a same-shape corruption that still parses would otherwise promote,
# and the manifest would then certify the wrong bytes with a perfectly valid hash (phase1audit F-06).
# NOT identical(): a data.table's .internal.selfref is a pointer that necessarily differs after a
# round trip, so identical() is FALSE for every correct write (measured, R 4.5.1 / data.table 1.17.8)
# and the first honest rerun would abort. all.equal() ignores that pointer; tolerance = 0 keeps the
# numeric comparison exact. Measured to catch row-order changes and dropped attributes.
if (!isTRUE(all.equal(readRDS(tmp), out, tolerance = 0))) {
  unlink(tmp); stop("staged RDS did not read back identically")
}
# file.rename() DOES replace an existing destination on this platform (verified, R 4.5.1 on Windows:
# it returns TRUE and the destination holds the new bytes), so one rename installs the table without
# any interval in which the path is missing.
# The certificate is BUILT AND STAGED FIRST, from the staged RDS, and both are installed together.
# Installing the table and only then constructing the manifest left a window -- a kill, a hashing
# failure, a write failure -- in which the canonical path held a table with no certificate at all.
# No consumer inspects that certificate, so such a table reads as current; the whole point of
# removing the previous outputs up front was that a failure should be legible as a missing input.
.tmp_md5 <- unname(tools::md5sum(tmp))
.mtmp <- sprintf("%s.manifest.tsv.%d.tmp", OUT, Sys.getpid())
fwrite(data.table(file = basename(OUT), rows = nrow(out),
                  md5 = .tmp_md5,
                  source = basename(IN), source_md5 = row$md5[1],
                  # SNP is decided by these two as much as by the extraction itself; either can change
                  # the downstream join key while source_md5 stays identical.
                  # The map that actually determines SNP. RSMAP is the position-only file this
                  # script REFUSES to use; hashing it recorded nothing about the rsIDs in this table
                  # and would have changed the certificate when an unused file changed.
                  # The PATH, not the basename. Both `file` and this field read
                  # "ukb_cis_instruments.rds", which hid that the allele map came from
                  # results/repro_phase1_1682/ while the output goes to results/phase1_hg38/ --
                  # two different tables with one name, in a certificate whose purpose is to say
                  # which inputs produced this output.
                  bim_source = BIM,
                  bim_md5 = BIM_MD5,
                  # WHICH VERSION of this script. `builder` gave a filename, which identifies the
                  # code no more than the basename above identified the map.
                  # THIS script, wherever it was run from, hashed at startup. A hard-coded
                  # canonical path records the canonical file even when a modified copy is what
                  # executed, and re-reading at the end misses an edit made during the run.
                  builder_md5 = SELF_MD5,
                  extraction_manifest_md5 = MANI_MD5,
                  extraction_summary_md5  = SUMMARY_MD5,
                  extractor_code_md5      = EXTRACTOR_MD5,
                  extraction_run_signature = EXTRACT_RUN_SIG,
                  # 5 Medium: the chain decides every GRCh37<->GRCh38 statement this table rests on,
                  # in 08 and again in this script's own check, and it was recorded by neither.
                  extraction_chain_md5 = EXTRACT_CHAIN_MD5,
                  gene_coords_source = basename(COORDS_TBL),
                  gene_coords_md5    = EXTRACT_COORDS_MD5,
                  liftover_chain_md5 = tryCatch(
                    unname(tools::md5sum("data/ld_reference/hg19ToHg38.over.chain.gz")),
                    error = function(e) NA_character_),
                  # The coordinate columns changed BUILD without changing name: gene_start and
                  # gene_end are GRCh38 here and were GRCh37 in the pinned 20-column table (A1BG
                  # 58,345,178 vs 58,856,544). A consumer applying the old contract lifts them a
                  # second time, and nothing in the artifact contradicts it -- the retired
                  # cascade/03_gene_coords.R called them pos37 for exactly that reason.
                  gene_coord_build = "GRCh38",
                  chrpos_build = "GRCh38",
                  id_build = "GRCh37",
                  disagreements_md5 = {
                    .d <- file.path(OUTD, "rsid_name_disagreements.csv")
                    if (file.exists(.d)) unname(tools::md5sum(.d)) else ""
                  },
                  rsid_pct = round(100 * n_rs / n_rows_certified, 2),
                  # WHICH proteins have no usable rsID, not merely how many. There is no longer a
                  # gate to bypass -- the count is a measured property of the release, so what the
                  # certificate has to carry is the identity of what is missing, not a flag saying
                  # someone allowed it.
                  n_proteins_lost = length(gone),
                  proteins_lost   = paste(sort(gone), collapse = "|"),
                  lost_sidecar_md5 = {
                    .f <- file.path(OUTD, "instruments_without_rsid.csv")
                    if (file.exists(.f)) unname(tools::md5sum(.f)) else ""
                  },
                  builder = "09_phase1_harmonise_hg38.R",
                  built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
       .mtmp, sep = "\t")
if (!file.exists(.mtmp)) { unlink(c(tmp, .mtmp)); stop("could not stage the certificate at ", .mtmp) }
# Both staged and verified; install them. Two renames is the minimum. If the second fails the table
# is present but uncertified, which is the state this script's own gate at the top refuses -- so the
# next run stops on a missing certificate rather than reading an uncertified table as canonical.
if (!file.rename(tmp, OUT)) {
  unlink(c(tmp, .mtmp))
  stop("cannot install ", OUT, "; the staged copies have been removed.")
}
if (!file.rename(.mtmp, paste0(OUT, ".manifest.tsv"))) {
  # The RDS goes too. Leaving it was justified by "the next run stops on a missing certificate",
  # which is false: the gate at the top of this script checks 08's phase1_OUTPUTS.tsv, not this
  # certificate -- so an uncertified RDS at the canonical path reads as a normal analysis input.
  unlink(c(.mtmp, OUT))
  # CHECKED. On Windows a lock can make that deletion fail, and the message would then claim a
  # removal that did not happen -- leaving the canonical path holding exactly the uncertified table
  # it says it removed, which everything downstream reads normally.
  if (file.exists(OUT))
    stop(basename(OUT), "'s certificate could not be written AND the table could not be removed. ",
         "It is still at ", OUT, ", uncertified and indistinguishable from a certified one to ",
         "everything downstream. DELETE IT BY HAND before running anything else.")
  stop(basename(OUT), "'s certificate could not be written, so the table has been removed too: an ",
       "uncertified RDS at that path is indistinguishable from a certified one to everything ",
       "downstream. Rerun this script.")
}
cat(sprintf("\nwrote %s\n  %s instruments | %s genes | rsID %.1f%%\n",
            OUT, format(nrow(out), big.mark = ","), format(uniqueN(out$gene), big.mark = ","),
            100 * n_rs / n_rows_certified))
cat("[NOTE] phase1/06_add_rsid_to_instruments.R and phase1/12_phase1_r2.R are the frozen originals of\n",
    "       this transformation. They are NOT re-run: 06 carries a freeze guard and both write over\n",
    "       the artifact that documents the submitted analysis.\n", sep = "")
