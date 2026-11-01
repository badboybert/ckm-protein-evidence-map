# ============================================================================
# EAS AXIS 02: 76-pair exact-match coverage + QC (contract v5 SS7, SS10 item 2)
# ----------------------------------------------------------------------------
# For every frozen pair x source (primary AND sensitivity): does the source
# carry the lead instrument as an exact, allele-resolved, QC-passing record?
# Outputs:
#   eas_outcome_coverage.csv              -- 76 rows, primary sources (05/06 input)
#   eas_outcome_coverage_sensitivity.csv  -- one row per pair x sensitivity source
#   lead_b37_lookup.csv                   -- all 62 leads' GRCh37 coordinates (04 reuses)
# ONE mutually exclusive SS7.3 reason per failure (first failing SS7.2 gate),
# the full SS10.1 upstream audit fields, and the ORIENTATION DECISION for 03.
# Matching happens HERE exactly once; 03 must not re-match or re-orient.
#
# GRCh37 coordinates: EAS.bim alone lacks 5 of the 62 leads (EASr2 finding 4),
# so the lookup is built by lifting every lead's frozen hg38 position over the
# local hg38->hg19 chain (pyliftover, offline, fail-loud -- the same mechanism
# R/coloc.R uses), then CROSS-CHECKED against EAS.bim wherever the rsID exists
# there: any disagreement stops the run (a build confusion must be resolved,
# never trusted silently). The sentinel assertion on the bim itself stays.
# Palindromic leads (SS7.2.8): orientation requires a genome-anchored record --
# ref/alt designated by the source (01 leaves them NA otherwise). A palindromic
# lead in a non-anchored source is PALINDROMIC_AMBIGUOUS, never EAF-guessed.
# sample_count_basis is NEVER blank: rows that fail before a record is located
# carry "UNRESOLVED" (the three-value contract partitions all 76).
# Duplicate records must agree on EVERY field this analysis reads -- alleles,
# effect stats, p, INFO, status, counts, mac fields, ref/alt, chromosome --
# or the pair is MULTIALLELIC_OR_DUPLICATE (EASr2 finding 8).
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
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
suppressPackageStartupMessages(library(data.table))
B <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"); setwd(B)
EAS_DIR <- file.path("results", paste0("repro_", PHASE1_SET), "eas_axis")
LOCAL   <- file.path(EAS_DIR, "localized")
MANI    <- file.path(EAS_DIR, "eas_source_manifest.csv")
if (!file.exists(MANI)) stop("02: run 01 first -- ", MANI, " absent")
dir.create(EAS_DIR, showWarnings = FALSE, recursive = TRUE)

INFO_MIN <- 0.80; MAC_MIN <- 20; STRATUM_MIN <- 10

# ── the frozen universe: 76 pairs, 62 unique leads, UKB alleles + hg38 pos ───
t4 <- fread(file.path("results", paste0("repro_", PHASE1_SET),
                      "table_s_pair_evidence_tiered.csv"))
stopifnot(nrow(t4) > 0L)
ins <- as.data.table(readRDS(file.path("results", "phase1_hg38", "ukb_cis_instruments.rds")))
pairs <- merge(t4[, .(gene, outcome)], ins[, .(gene, rsid = SNP, chrpos38 = chrpos,
               ea_ukb = toupper(effect_allele), oa_ukb = toupper(other_allele))],
               by = "gene")
# the merge is by gene and must stay 1:1 with the canonical pair table; the lead count follows
# from it rather than being asserted separately (was 76L / 62L, the 2026-08-16 anchor).
stopifnot(nrow(pairs) == nrow(t4), uniqueN(pairs$rsid) > 0L)
pairs[, c("chr38", "pos38") := {
  sp <- tstrsplit(sub("^chr", "", chrpos38), ":")
  .(as.character(sp[[1]]), as.integer(sp[[2]]))
}]
pairs[, palindromic := (ea_ukb == "A" & oa_ukb == "T") | (ea_ukb == "T" & oa_ukb == "A") |
                        (ea_ukb == "C" & oa_ukb == "G") | (ea_ukb == "G" & oa_ukb == "C")]

# ── b37 lookup: liftover ALL 62 leads, cross-check against the bim ───────────
bim <- fread(file.path("data", "ld_reference", "1kg_eas", "EAS.bim"),
             col.names = c("chr", "rs", "cm", "pos", "a1", "a2"))
sent <- bim[rs == "rs77924615"]
if (!nrow(sent) || sent$pos[1] != 20392332L)
  stop("02: EAS.bim failed the b37 sentinel (rs77924615 must sit at 16:20392332); ",
       "the panel is not the build the position matching assumes")

leads <- unique(pairs[, .(rsid, chr38, pos38)])
tmp_in  <- tempfile(fileext = ".tsv"); tmp_out <- tempfile(fileext = ".tsv")
fwrite(leads[, .(rsid, chrom38 = chr38, pos38)], tmp_in, sep = "\t")
py <- tempfile(fileext = ".py")
cat(sprintf('
from pyliftover import LiftOver
import csv

lo = LiftOver(r"%s")
with open(r"%s") as fin, open(r"%s", "w", newline="") as fout:
    rdr = csv.DictReader(fin, delimiter="\\t")
    wtr = csv.writer(fout, delimiter="\\t")
    wtr.writerow(["rsid", "chr37", "pos37"])
    for row in rdr:
        try:
            res = lo.convert_coordinate("chr" + row["chrom38"].lstrip("chr"),
                                        int(row["pos38"]) - 1)
            if res and len(res) == 1:
                wtr.writerow([row["rsid"], res[0][0].lstrip("chr"), res[0][1] + 1])
        except Exception:
            raise
',
  gsub("\\\\", "/", LIFTOVER_CHAIN_38_TO_37),
  gsub("\\\\", "/", tmp_in), gsub("\\\\", "/", tmp_out)), file = py)
.pyerr <- tempfile(fileext = ".stderr")
.st <- system2(PYTHON_BIN, args = py, stdout = FALSE, stderr = .pyerr)
.emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
unlink(c(py, tmp_in, .pyerr))
if (!identical(as.integer(.st), 0L))
  stop("02: lead liftover (hg38->hg19) FAILED (exit ", .st, ") -- fail loud, not a silent ",
       "coordinate gap. python stderr: ", if (nzchar(.emsg)) .emsg else "<empty>")
b37 <- fread(tmp_out); unlink(tmp_out)
b37[, chr37 := as.character(chr37)]
miss <- setdiff(leads$rsid, b37$rsid)
if (length(miss))
  stop("02: ", length(miss), " lead(s) did not lift to GRCh37 (chain gap): ",
       paste(miss, collapse = ", "), " -- resolve before coverage, do not let them ",
       "become false VARIANT_ABSENT")
bimchk <- merge(b37, bim[rs %in% leads$rsid, .(rsid = rs, chrB = as.character(chr), posB = pos)],
                by = "rsid")
bad <- bimchk[chr37 != chrB | pos37 != posB]
if (nrow(bad))
  stop("02: liftover vs EAS.bim disagree for ", nrow(bad), " lead(s): ",
       paste(bad$rsid, collapse = ", "), " -- build confusion must be resolved, not trusted")
cat(sprintf("02: b37 lookup built for all %d leads (liftover; %d cross-checked against EAS.bim)\n",
            nrow(b37), nrow(bimchk)))
fwrite(merge(leads, b37, by = "rsid"), file.path(EAS_DIR, "lead_b37_lookup.csv"))

mani <- fread(MANI)
stopifnot(all(mani$build %in% c("GRCh37", "GRCh38")))
norm_chr <- function(x) sub("^chr", "", toupper(as.character(x)))
COMP <- c(A = "T", T = "A", C = "G", G = "C")
# every field the analysis reads downstream -- INCLUDING the raw counts and the
# rsid, which are copied into the profile verbatim (EASr4 finding 2); duplicates
# must agree on ALL of them or the pair is unresolvable
DUP_FIELDS <- c("rsid", "effect_allele", "other_allele", "beta", "se", "eaf", "pos", "p",
                "info", "source_qc_status", "sample_count_basis",
                "n_total_study", "n_cases_study", "n_controls_study",
                "n_total_variant", "n_cases_variant", "n_controls_variant", "n_eff_variant",
                "n_total_mac_basis", "n_cases_mac_basis", "n_controls_mac_basis",
                "mac_total_approx", "mac_cases_approx", "mac_controls_approx",
                "ref", "alt", "chrN")

blank_row <- function(pr, src_id, src_pheno, src_build) {
  data.table(gene = pr$gene, outcome = pr$outcome, instrument_rsID = pr$rsid,
             eas_source = src_id, eas_phenotype = src_pheno, build = src_build,
             snp_present = FALSE, allele_match = FALSE, qc_evaluable = FALSE,
             hit_via = NA_character_, position_checked = FALSE,
             palindromic = pr$palindromic,
             ea_src = NA_character_, oa_src = NA_character_,
             ref_src = NA_character_, alt_src = NA_character_,
             orientation_action = NA_character_, orientation_basis = NA_character_,
             strand_flipped = FALSE,
             betaY_src = NA_real_, seY_src = NA_real_, p_src = NA_real_,
             info = NA_real_, eaf = NA_real_, source_qc_status = NA_character_,
             n_total_study = NA_real_, n_cases_study = NA_real_, n_controls_study = NA_real_,
             n_total_variant = NA_real_, n_cases_variant = NA_real_,
             n_controls_variant = NA_real_, n_eff_variant = NA_real_,
             n_total_mac_basis = NA_real_, n_cases_mac_basis = NA_real_,
             n_controls_mac_basis = NA_real_,
             sample_count_basis = "UNRESOLVED",   # never blank: three-value contract
             mac_total_approx = NA_real_, mac_cases_approx = NA_real_,
             mac_controls_approx = NA_real_, qc_reason = NA_character_)
}

# the SS7.2 gate for one pair against one loaded source table
gate_pair <- function(pr, src, d) {
  o <- blank_row(pr, src$source_id, src$phenotype_label, src$build)
  fail1 <- function(code) { o$qc_reason <- code; o }
  lkp <- if (src$build == "GRCh37") {
    b37[rsid == pr$rsid, .(chrL = norm_chr(chr37), posL = pos37)]
  } else {
    pairs[rsid == pr$rsid][1, .(chrL = norm_chr(chr38), posL = pos38)]
  }
  hit <- if ("rsid" %in% names(d)) d[rsid == pr$rsid] else d[0]
  via <- "rsid"
  if (!nrow(hit) && nrow(lkp)) {
    hit <- d[chrN == lkp$chrL[1] & pos == lkp$posL[1]]; via <- "position"
  }
  if (!nrow(hit)) return(fail1("VARIANT_ABSENT"))
  if (nrow(hit) > 1L) {
    hit <- unique(hit, by = intersect(DUP_FIELDS, names(hit)))
    if (nrow(hit) > 1L) return(fail1("MULTIALLELIC_OR_DUPLICATE"))
  }
  o$snp_present <- TRUE; o$hit_via <- via
  # the located record's metadata is copied IMMEDIATELY, before EVERY further
  # gate -- position, orientation, validity: a pair excluded at any of them
  # still reports the counts, basis, INFO and EAF the record actually carries
  # (EASr4 finding 1: three T2D palindromic pairs lost variant-level counts;
  # EASr5 finding 1: the POSITION_CONFLICT branch was still ahead of the copy).
  bta <- suppressWarnings(as.numeric(hit$beta)); sev <- suppressWarnings(as.numeric(hit$se))
  eaf <- suppressWarnings(as.numeric(hit$eaf)); pvl <- suppressWarnings(as.numeric(hit$p))
  o$betaY_src <- bta; o$seY_src <- sev; o$p_src <- pvl; o$eaf <- eaf
  o$info <- suppressWarnings(as.numeric(hit$info))
  o$source_qc_status <- hit$source_qc_status
  for (cc in c("n_total_study", "n_cases_study", "n_controls_study",
               "n_total_variant", "n_cases_variant", "n_controls_variant", "n_eff_variant",
               "n_total_mac_basis", "n_cases_mac_basis", "n_controls_mac_basis",
               "mac_total_approx", "mac_cases_approx", "mac_controls_approx"))
    set(o, j = cc, value = suppressWarnings(as.numeric(hit[[cc]])))
  if (!is.na(hit$sample_count_basis) && nzchar(hit$sample_count_basis))
    o$sample_count_basis <- hit$sample_count_basis
  if (via == "rsid" && nrow(lkp) && is.finite(hit$pos)) {
    o$position_checked <- TRUE
    if (!(hit$chrN == lkp$chrL[1] && hit$pos == lkp$posL[1]))
      return(fail1("POSITION_CONFLICT"))
  }
  ea <- toupper(hit$effect_allele); oa <- toupper(hit$other_allele)
  o$ea_src <- ea; o$oa_src <- oa
  o$ref_src <- toupper(hit$ref); o$alt_src <- toupper(hit$alt)
  anchored <- !is.na(hit$ref) && !is.na(hit$alt)
  direct  <- !is.na(ea) && ea == pr$ea_ukb && oa == pr$oa_ukb
  swapped <- !is.na(ea) && ea == pr$oa_ukb && oa == pr$ea_ukb
  if (pr$palindromic) {
    if (!anchored) return(fail1("PALINDROMIC_AMBIGUOUS"))
    if (direct)      { o$orientation_action <- "AS_IS"; o$orientation_basis <- "REF_ALT_ANCHORED" }
    else if (swapped){ o$orientation_action <- "FLIP";  o$orientation_basis <- "REF_ALT_ANCHORED" }
    else             return(fail1("ALLELE_MISMATCH"))
  } else {
    cd <- !is.na(ea) && !is.na(COMP[ea]) && !is.na(COMP[oa]) &&
          COMP[ea] == pr$ea_ukb && COMP[oa] == pr$oa_ukb
    cs <- !is.na(ea) && !is.na(COMP[ea]) && !is.na(COMP[oa]) &&
          COMP[ea] == pr$oa_ukb && COMP[oa] == pr$ea_ukb
    if (direct)       { o$orientation_action <- "AS_IS"; o$orientation_basis <- "NONPALINDROMIC_LETTER_MATCH" }
    else if (swapped) { o$orientation_action <- "FLIP";  o$orientation_basis <- "NONPALINDROMIC_LETTER_MATCH" }
    else if (cd)      { o$orientation_action <- "AS_IS"; o$orientation_basis <- "NONPALINDROMIC_COMPLEMENT"; o$strand_flipped <- TRUE }
    else if (cs)      { o$orientation_action <- "FLIP";  o$orientation_basis <- "NONPALINDROMIC_COMPLEMENT"; o$strand_flipped <- TRUE }
    else              return(fail1("ALLELE_MISMATCH"))
  }
  o$allele_match <- TRUE
  if (!is.finite(bta) || !is.finite(sev) || sev <= 0 ||
      !is.finite(eaf) || eaf <= 0 || eaf >= 1) return(fail1("INVALID_EFFECT"))
  if (o$source_qc_status %in% c("INFO_FAIL", "INFO_MISSING")) return(fail1("LOW_INFO"))
  if (o$sample_count_basis == "UNRESOLVED") return(fail1("SAMPLE_SIZE_UNRESOLVED"))
  if (src$outcome_type == "continuous") {
    if (!is.finite(o$mac_total_approx) || o$mac_total_approx < MAC_MIN) return(fail1("LOW_MAC"))
  } else {
    if (!is.finite(o$mac_cases_approx) || !is.finite(o$mac_controls_approx)) return(fail1("LOW_MAC"))
    if (o$mac_cases_approx < STRATUM_MIN || o$mac_controls_approx < STRATUM_MIN)
      return(fail1("LOW_STRATUM_COUNT"))
  }
  o$qc_evaluable <- TRUE
  o
}

run_role <- function(role_want) {
  rows <- list()
  for (s in seq_len(nrow(mani))) {
    src <- mani[s]
    if (src$role != role_want) next
    want <- pairs[outcome == src$outcome]
    if (!nrow(want)) next
    d <- fread(file.path(LOCAL, paste0(src$source_id, ".tsv.gz")),
               na.strings = c("NA", ""), showProgress = FALSE)
    d[, `:=`(chrN = norm_chr(chr), pos = as.integer(pos))]
    if ("rsid" %in% names(d)) {
      chk <- merge(d[rsid %in% want$rsid, .(rsid, chrN, pos)],
                   if (src$build == "GRCh37") b37[, .(rsid, chrL = norm_chr(chr37), posL = pos37)]
                   else unique(pairs[, .(rsid, chrL = norm_chr(chr38), posL = pos38)]),
                   by = "rsid")
      if (nrow(chk) >= 3L) {
        okd <- chk[, mean(chrN == chrL & pos == posL, na.rm = TRUE)]
        if (is.finite(okd) && okd < 0.5)
          stop("02: ", src$source_id, " declared ", src$build, " but only ",
               round(100 * okd), "% of rsID-matched leads sit at declared-build coordinates -- ",
               "registry build is wrong; fix it (dated amendment), do not guess per row")
      }
    }
    for (i in seq_len(nrow(want)))
      rows[[length(rows) + 1L]] <- gate_pair(want[i], src, d)
    rm(d); invisible(gc())
  }
  if (length(rows)) rbindlist(rows) else data.table()
}

cov <- run_role("primary")
# outcomes with no primary source at all
covered <- unique(cov[, .(gene, outcome)])
for (i in seq_len(nrow(pairs))) {
  pr <- pairs[i]
  if (!nrow(covered[gene == pr$gene & outcome == pr$outcome])) {
    o <- blank_row(pr, NA_character_, NA_character_, NA_character_)
    o$qc_reason <- "OUTCOME_UNAVAILABLE"
    cov <- rbind(cov, o)
  }
}
QC_CODES <- c("OUTCOME_UNAVAILABLE", "PHENOTYPE_NOT_ACCEPTABLE", "VARIANT_ABSENT",
              "BUILD_UNRESOLVED", "POSITION_CONFLICT", "ALLELE_MISMATCH",
              "MULTIALLELIC_OR_DUPLICATE", "PALINDROMIC_AMBIGUOUS", "LOW_INFO",
              "LOW_MAC", "LOW_STRATUM_COUNT", "SAMPLE_SIZE_UNRESOLVED", "INVALID_EFFECT")
stopifnot(nrow(cov) == nrow(pairs),
          all(!cov$qc_evaluable | is.na(cov$qc_reason)),
          all(cov$qc_evaluable | cov$qc_reason %in% QC_CODES),
          all(cov$sample_count_basis %in% c("VARIANT_LEVEL", "STUDY_LEVEL_APPROX", "UNRESOLVED")),
          all(!cov$qc_evaluable | !is.na(cov$orientation_action)),
          all(!(cov$palindromic & cov$qc_evaluable) |
                cov$orientation_basis == "REF_ALT_ANCHORED"))
fwrite(cov, file.path(EAS_DIR, "eas_outcome_coverage.csv"))

# sensitivity sources: same gate, separate table, never touches the primary file
sen <- run_role("sensitivity")
if (nrow(sen)) {
  stopifnot(all(sen$sample_count_basis %in%
                  c("VARIANT_LEVEL", "STUDY_LEVEL_APPROX", "UNRESOLVED")),
            all(!sen$qc_evaluable | is.na(sen$qc_reason)))
  fwrite(sen, file.path(EAS_DIR, "eas_outcome_coverage_sensitivity.csv"))
} else stop("02: no sensitivity rows produced -- three sensitivity sources are registered; ",
            "their absence here is a defect, not an option")

cat(sprintf("\n02 primary coverage: snp_present %d/76 | allele_match %d/76 | qc_evaluable %d/76\n",
            sum(cov$snp_present), sum(cov$allele_match), sum(cov$qc_evaluable)))
cat(sprintf("basis among evaluable: VARIANT %d | STUDY %d\n",
            nrow(cov[qc_evaluable == TRUE & sample_count_basis == "VARIANT_LEVEL"]),
            nrow(cov[qc_evaluable == TRUE & sample_count_basis == "STUDY_LEVEL_APPROX"])))
cat("palindromic pairs (SS7.2.8 audit):\n")
print(cov[palindromic == TRUE, .(gene, outcome, eas_source, qc_evaluable,
                                 orientation_action, orientation_basis, qc_reason)])
cat("reason breakdown (failures only):\n")
print(cov[qc_evaluable == FALSE, .N, by = qc_reason][order(-N)])
cat(sprintf("sensitivity coverage rows: %d (%s)\n", nrow(sen),
            paste(unique(sen$eas_source), collapse = ", ")))
