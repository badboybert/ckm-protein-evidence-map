# Phase 3-C Extension: SuSiE-coloc for 14 eGFR Strong Hits
# ============================================================================
# The original 15_susie_coloc.R excludes eGFR (rsid_only in GWAS_META).
# This script uses the local CKDGen file to run SuSiE fine-mapping coloc
# for eGFR strong hits, then merges with existing results.
#
# Usage: cd "<E>/20260331-CKM/Paper 3" && Rscript scripts/paper3.1_ckm/phase3/26_susie_coloc_egfr.R
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
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(ieugwasr)
})

if (!requireNamespace("susieR", quietly = TRUE))
  stop("Package 'susieR' not installed. Run: install.packages('susieR')")
library(susieR)

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
COLOC_SOURCE_ONLY <- TRUE   # define functions only; R/coloc.R:1213 runs main() on a bare source()
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R"))   # liftover_to_hg37(): ONE definition, shared
rm(COLOC_SOURCE_ONLY)
# Root-cause LD sign-alignment helper (2026-07-17, editor ruling 2b). In the head so the A" diagnostics/27
# (which eval this engine's pre-main-loop block) load strand_align_R too, not just the full-source cascade.
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/strand_align.R"))

cat("=== Phase 3-C Extension: SuSiE-coloc for eGFR ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── Settings ─────────────────────────────────────────────────────────────────
EGFR_FILE  <- file.path(OUTCOME_DIR, "CKDGen_eGFR_EA_hg37.txt.gz")
# eGFR id from config (corrected cascade = Wuttke2019_eGFR_EA); legacy fallback + warn if absent.
EGFR_ID    <- if (exists("OUTCOMES") && !is.null(OUTCOMES[["eGFR"]])) OUTCOMES[["eGFR"]] else { warning("[eGFR-id] OUTCOMES$eGFR missing in config; legacy fallback ebi-a-GCST003372"); "ebi-a-GCST003372" }
# 2026-06-30 (gate obj #6): derive eGFR N from the data (median n_total_sum), not a stale literal,
# consistent with the coloc/prior regional-N path. Fail loud if the column is absent.
# LAZY. This reads the whole n_total_sum column (8,885,712 rows out of a 192.6 MB gzip) and is only
# a FALLBACK: the fit uses the regional median carried through the merge, so on current data the
# value is never consulted. It was paid on every evaluation of this head -- by 27/29/31/35/36/37/41,
# which slice and eval it, and by every PSOCK worker. delayedAssign() runs the same code, only when
# something actually asks for EGFR_N. The eager cat() that used to print it here is what forced it.
delayedAssign("EGFR_N", {
  .nc <- tryCatch(data.table::fread(EGFR_FILE, select = "n_total_sum", showProgress = FALSE)[[1]], error = function(e) NULL)
  if (is.null(.nc) || !length(.nc)) stop("26 eGFR SuSiE: cannot read n_total_sum from eGFR file — refusing stale N literal (was 567460)")
  as.integer(median(.nc, na.rm = TRUE))
})
MAF_MIN    <- 0.01

stopifnot(file.exists(EGFR_FILE))

# ── SuSiE iteration ceiling (override via SUSIE_MAX_ITER) ─────────────────────
# DEFAULT 1000, matching 15_susie_coloc.R. No eGFR strong pair failed to converge at 1000, so
# this engine is not the reason the ceiling exists; it carries the parameter so that "the strong
# set ran at one ceiling" stays true across both of its halves.
if (!exists("SUSIE_MAX_ITER")) SUSIE_MAX_ITER <- 1000L

if (!exists("SUSIE_EGFR_OUT")) SUSIE_EGFR_OUT <- file.path(RESULTS_DIR, "phase3_susie")   # override for corrected
out_dir <- SUSIE_EGFR_OUT
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ────────────────────────────────────────────────────────────────
if (!exists("SUSIE_EGFR_COLOC")) SUSIE_EGFR_COLOC <- file.path(PHASE3_DIR, "coloc_results_final.rds")   # override for corrected (v6)
if (!exists("SUSIE_EGFR_INSTR")) SUSIE_EGFR_INSTR <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")
coloc_final <- readRDS(SUSIE_EGFR_COLOC)
# THE COMBINED TABLE MUST NAME ITS INPUTS. This read whatever was at the path and ran SuSiE on it, so
# a schema-compatible table from another generation -- or one produced before 06 started stamping --
# would be mixed with the current instruments and outcomes without anything noticing. 06 records the
# non-eGFR and eGFR input keys on the union it writes; requiring them here is what makes that useful.
# Only their PRESENCE is required, not a recomputation: 06 already verified each stratum against the
# live inputs, and recomputing here would re-read the archives for no additional guarantee.
for (.a in c("coloc_inputs_key_noneGFR", "coloc_inputs_key_eGFR")) {
  .v <- attr(coloc_final, .a, exact = TRUE)
  if (is.null(.v) || !nzchar(as.character(.v)[1]))
    stop("26 eGFR SuSiE: ", basename(SUSIE_EGFR_COLOC), " carries no ", .a, ", so the inputs its ",
         "posteriors were computed from cannot be established. It predates the stamping, or it came ",
         "from another generation. Re-run 06_combine_coloc.R.")
}

cat(sprintf("[provenance] v7 inputs keys: non-eGFR %s | eGFR %s\n",
            substr(attr(coloc_final, "coloc_inputs_key_noneGFR", exact = TRUE), 1, 12),
            substr(attr(coloc_final, "coloc_inputs_key_eGFR", exact = TRUE), 1, 12)))
instruments <- readRDS(SUSIE_EGFR_INSTR)
setDT(coloc_final); setDT(instruments)

strong_egfr <- coloc_final[is_coloc_pass(PP.H4, PP.H3) & outcome_id == EGFR_ID][order(-PP.H4)]   # dual coloc rule (D1; centralized 2026-07-05)
if (exists("SUSIE_LIMIT")) strong_egfr <- strong_egfr[seq_len(min(SUSIE_LIMIT, nrow(strong_egfr)))]   # smoke cap (scope only)

# A BYTE-EXACT INPUT KEY WAS TRIED HERE AND WITHDRAWN (2026-08-03, PI decision). Round 6 found that
# the attribute check above proves only that SOMETHING stamped the file, not that it was this
# generation, so the stored key was compared against one recomputed from the live inputs. On real
# data it blocked immediately: coloc_results_v7.rds keys 4a86fc8d43b8, the tree now keys
# 7ee365b481ad, because four covered files were edited after v7 was written -- R/coloc.R and
# R/00_config.R at 20:56/21:02 on 2026-08-02, plus localize/22 and /23. Those edits were two cat()
# labels (Nielsen -> outcome_id) and comment cleanups; the data was verified unchanged at the time
# (REBUILD_LEDGER 3528, 3542).
#
# A byte fingerprint cannot tell a log-label edit from a computation edit, so enforcing it here
# would demand a six-hour 04->05->06 re-run to certify a change of wording. This project has
# already been bitten by exactly that and wrote it down: _coloc_inputs_key.R fingerprints
# mr_suggestive on its FIELDS rather than its bytes because re-saving an identical object with
# different RDS compression 'altered the byte hash and made every worker reject unchanged
# scientific work'. The same argument applies to the engine source.
#
# So the presence check above stands, and the gap it leaves -- a foreign but stamped artifact
# passes -- is accepted and recorded rather than closed with a check that cannot distinguish the
# two cases. Closing it properly needs a key over what the inputs MEAN, not over their bytes.
cat(sprintf("eGFR strong hits: %d pairs\n\n", nrow(strong_egfr)))

# ── rsID map ─────────────────────────────────────────────────────────────────
load_rsid_map <- function() {
  f <- file.path(BASE_DIR, "cache/cis_rsid_map.tsv")
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, colClasses = list(character = "chrom"))
  # KEYED ON THE ALLELES TOO -- the same edit as R/coloc.R and 05_coloc_egfr.R, which this script
  # carries its own copy of. A position-only key collapses two variants at one coordinate to one
  # rsID: at 1:2102322 the map holds T/G -> rs201135834 and T/TG -> rs34905415, both present in
  # NADK's current archive, and after the MAF/INFO filter only the indel survives -- so it took the
  # SNV's rsID, matched the wrong outcome row, and the position fallback could not rescue it.
  dt[, lookup_key := paste(chrom, pos, toupper(allele0), toupper(allele1), sep = ":")]
  dt <- dt[order(as.integer(chrom), pos)]
  # ONE RSID -> ONE POSITION, BUT EVERY ALLELE PAIR AT THAT POSITION -- same edit and same reason
  # as R/coloc.R. `.SD[1]` collapsed the legacy-multiple-POSITIONS case it was written for AND the
  # distinct-variants case it was not: 1:1251285 carries G/A and G/T under one rsID (rs6657544).
  dt <- dt[, .SD[chrom == chrom[1L] & pos == pos[1L]], by = rsid]
  dt <- dt[, .SD[1], by = lookup_key]
  cat(sprintf("[rsID map] %s entries\n\n", formatC(nrow(dt), big.mark = ",")))
  dt
}
rsid_map <- load_rsid_map()

# ── Helper functions (same as 21_rerun_egfr_coloc.R) ────────────────────────
complement_allele <- function(a) chartr("ACGTacgt", "TGCAtgca", toupper(a))
is_palindromic    <- function(a1, a2) complement_allele(a1) == toupper(a2)

harmonize_alleles <- function(dt, ea_p, oa_p, ea_o, oa_o, beta_o) {
  ep <- toupper(dt[[ea_p]]); op <- toupper(dt[[oa_p]])
  eo <- toupper(dt[[ea_o]]); oo <- toupper(dt[[oa_o]])
  comp_ep <- complement_allele(ep); comp_op <- complement_allele(op)
  same      <- eo == ep  & oo == op
  flip_beta <- eo == op  & oo == ep
  same_str  <- eo == comp_ep & oo == comp_op
  flip_str  <- eo == comp_op & oo == comp_ep
  palin      <- is_palindromic(ep, op)
  # SNVs ONLY in the strand branches -- the same edit as R/coloc.R and 05_coloc_egfr.R, whose
  # harmonisation this script copies. complement_allele() is chartr(), which complements without
  # reversing ("GGA" -> "CCT", not "TCC"), so a genuinely reverse-complemented indel is dropped while
  # an unrelated multi-base pair can match by coincidence. Reversing the string is not the fix: a
  # point liftover does not re-anchor an indel on the target reference. Exact and swapped matches
  # still accept multi-base alleles; only the strand inference is restricted.
  strand_snv <- nchar(ep) == 1L & nchar(op) == 1L & nchar(eo) == 1L & nchar(oo) == 1L
  keep       <- allele_compatible(ep, op, eo, oo)   # engine predicate (R/coloc.R, sourced above): the rsID merge below asks the same question, and two copies would drift
  needs_flip <- flip_beta | (strand_snv & !palin & flip_str)
  dt <- dt[keep]
  if (nrow(dt) > 0 && any(needs_flip[keep]))
    dt[needs_flip[keep], (beta_o) := -get(beta_o)]
  dt
}

get_pqtl_region <- function(gene, chrom, region_start, region_end) {
  tar_file <- resolve_ukbppp_tar(gene)   # 2026-07-01 (gate finding 1): canonical assay by instrument OID, not [1]
  if (is.null(tar_file)) return(NULL)

  py <- tempfile(fileext = ".py")
  writeLines(sprintf('
import tarfile, gzip, pandas as pd, sys
with tarfile.open(r"%s", "r:") as tf:
    for m in tf.getmembers():
        if "chr%s_" in m.name and m.name.endswith(".gz"):
            f = tf.extractfile(m)
            with gzip.open(f) as gz:
                df = pd.read_csv(gz, sep=" ")
            df = df[(df["GENPOS"] >= %d) & (df["GENPOS"] <= %d)]
            df.to_csv(sys.stdout, index=False)
            break
', gsub("\\\\", "/", tar_file), chrom, region_start, region_end), py)

  # stderr KEPT -- same edit and same reason as R/coloc.R: a discarded stderr turned a whole failed
  # run into a bare exit code (2026-08-02, exit 322 on all 1,380 pairs; six hours to establish
  # that python itself was fine and the machine had simply run out of commit).
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(.pyerr)
  unlink(py)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0)  # 2026-07-01 (closure B1): extraction failure -> fail loud, not a NA sensitivity cell
    fail_loud_data(sprintf("26 eGFR SuSiE extraction FAILED (python exit %s) -- not a legit skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: no diagnostic from the process, which points at process creation rather than the script>"))
  if (length(out) <= 1) return(NULL)

  dt <- suppressWarnings(tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("26 eGFR region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e)))))  # 2026-07-01 (closure sweep): parallel of 03/06
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("26 eGFR region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")

  dt[, chrpos := paste0(chrom, ":", GENPOS)]
  if ("A1FREQ" %in% names(dt))
    dt <- dt[!is.na(A1FREQ) & A1FREQ >= MAF_MIN & A1FREQ <= 1 - MAF_MIN]

  if (!is.null(rsid_map)) {
    dt[, lookup_key := paste(chrpos, toupper(ALLELE0), toupper(ALLELE1), sep = ":")]
    dt <- merge(dt, rsid_map[, .(lookup_key, rsid)],
                by = "lookup_key", all.x = TRUE)
  } else {
    dt[, rsid := NA_character_]
  }
  dt
}

get_egfr_region <- function(chrom, region_start, region_end) {
  py <- tempfile(fileext = ".py")
  writeLines(sprintf('
import gzip, sys
chrom_target = "%s"
start = %d
end   = %d
with gzip.open(r"%s", "rt") as f:
    hdr = f.readline().strip()
    cols = hdr.split()
    ci = cols.index("Chr")
    pi = cols.index("Pos_b37")
    print("\\t".join(cols))
    for line in f:
        row = line.strip().split()
        if len(row) < len(cols):
            raise ValueError("short GWAS association row")
        try:
            chr_val = row[ci].lstrip("chr")
            if chr_val == chrom_target and start <= int(row[pi]) <= end:
                print("\\t".join(row))
        except (ValueError, IndexError):
            raise
', as.character(chrom), region_start, region_end,
   gsub("\\\\", "/", EGFR_FILE)), py)

  # stderr KEPT -- same edit and same reason as R/coloc.R: a discarded stderr turned a whole failed
  # run into a bare exit code (2026-08-02, exit 322 on all 1,380 pairs; six hours to establish
  # that python itself was fine and the machine had simply run out of commit).
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(.pyerr)
  unlink(py)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0)  # 2026-07-01 (closure B1): extraction failure -> fail loud, not a NA sensitivity cell
    fail_loud_data(sprintf("26 eGFR SuSiE extraction FAILED (python exit %s) -- not a legit skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: no diagnostic from the process, which points at process creation rather than the script>"))
  if (length(out) <= 1) return(NULL)

  dt <- suppressWarnings(tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("26 eGFR region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e)))))  # 2026-07-01 (closure sweep): parallel of 03/06
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("26 eGFR region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")

  setnames(dt, trimws(names(dt)))
  if ("Chr" %in% names(dt)) dt[, Chr := gsub("^chr", "", as.character(Chr))]
  dt[, rsid_out := RSID]
  dt[is.na(rsid_out) | !grepl("^rs", rsid_out), rsid_out := NA_character_]
  dt[, chrpos_out := paste0(Chr, ":", Pos_b37)]
  dt
}

rsid_pos_merge <- function(pqtl_sub, out_sub) {
  result <- list()
  matched_pqtl_var <- character(0)   # (chrpos:EA:OA) of the pQTL rows the primary merge already used

  has_p_rsid <- any(!is.na(pqtl_sub$rsid))
  has_o_rsid <- any(!is.na(out_sub$rsid_out))

  if (has_p_rsid && has_o_rsid) {
    # THE COORDINATE HAS TO AGREE TOO -- same edit and same reason as R/coloc.R and 05_coloc_egfr.R,
    # both of which this script copies. rsID alone accepted an outcome record with that rsID at a
    # different position, and marked the pQTL position matched so the right record was unreachable.
    p_r <- pqtl_sub[!is.na(rsid) & !is.na(chrpos_pqtl), .(rsid, chrpos_pqtl, pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl)]
    # out_n CARRIED THROUGH: this select() is what decides which outcome columns reach `common`,
    # and the regional-N block before susie_rss reads it there. Without it common_ld$out_n is NULL,
    # all(is.na(NULL)) is TRUE, and the fit silently falls back to the file-wide median.
    o_r <- out_sub[!is.na(rsid_out) & !is.na(chrpos_out), .(rsid_out, chrpos_out, out_beta, out_varbeta, ea_out, oa_out, out_n)]
    rsid_m <- merge(p_r, o_r, by.x = c("rsid", "chrpos_pqtl"), by.y = c("rsid_out", "chrpos_out"))
    if (nrow(rsid_m) > 0) {
      rsid_m[, snp_id := rsid]
      # ONLY A HARMONISABLE rsID HIT CLOSES A COORDINATE -- the same edit as R/coloc.R, whose matching
      # this script copies. An outcome record can carry the pQTL's rsID at the right position with
      # alleles that cannot describe the same variant; harmonisation drops it downstream, but the
      # coordinate had already been marked matched, so the compatible record at that same coordinate
      # was never reached and the pair silently lost the SNP.
      .ok <- allele_compatible(rsid_m$ea_pqtl, rsid_m$oa_pqtl, rsid_m$ea_out, rsid_m$oa_out)
      # (position, alleles), not position -- see R/coloc.R. Closing the coordinate also removed every
      # other pQTL variant at it, although those are different variants the fallback could match.
      matched_pqtl_var <- rsid_m[!is.na(chrpos_pqtl) & !is.na(.ok) & .ok,
                                 paste(chrpos_pqtl, toupper(ea_pqtl), toupper(oa_pqtl), sep = ":")]
      result[["rsid"]] <- rsid_m
    }
  }

  has_p_pos <- any(!is.na(pqtl_sub$chrpos_pqtl))
  has_o_pos <- any(!is.na(out_sub$chrpos_out))

  if (has_p_pos && has_o_pos) {
    p_fb <- pqtl_sub[!is.na(chrpos_pqtl) &
                     !(paste(chrpos_pqtl, toupper(ea_pqtl), toupper(oa_pqtl), sep = ":") %in% matched_pqtl_var),
                     .(chrpos_pqtl, pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl)]
    o_fb <- out_sub[!is.na(chrpos_out), .(chrpos_out, rsid_out, out_beta, out_varbeta, ea_out, oa_out, out_n)]  # 2026-07-02: carry rsid_out; 2026-08-03: carry out_n (see o_r above)
    pos_m <- merge(p_fb, o_fb, by.x = "chrpos_pqtl", by.y = "chrpos_out")
    if (nrow(pos_m) > 0) {
      # 2026-07-02 (editor-endorsed ACP1 recovery): a position-matched SNP is the SAME variant on both sides, so the
      # OUTCOME's dbSNP rsID is a valid rsID for LD (the local EUR panel is rsID-keyed). Use it when valid; fall back
      # to chrpos only when the outcome has none. Recovers pairs whose cis_rsid_map lacks coverage (e.g. ACP1~eGFR at
      # the chr2 start: pQTL side had 0 mapped rsIDs, but CKDGen supplies valid rsIDs that hit EUR.bim).
      pos_m[, snp_id := fifelse(!is.na(rsid_out) & grepl("^rs[0-9]+$", rsid_out), rsid_out, chrpos_pqtl)]
      result[["pos"]] <- pos_m
    }
  }

  if (length(result) == 0) return(data.table())
  rbindlist(result, use.names = TRUE, fill = TRUE)
}

# ── LD matrix helper (local plink preferred, API fallback) ───────────────────
get_ld_matrix_auto <- function(rsids, gene_chrom, max_snps_api = 400) {
  rsids <- unique(rsids[grepl("^rs[0-9]+$", rsids)])
  if (length(rsids) < 5) {
    cat("    [LD] < 5 rsIDs available\n"); return(NULL)
  }

  use_local <- exists("PLINK_BIN") && !is.null(PLINK_BIN) && file.exists(PLINK_BIN) &&
               exists("LD_EUR_DIR") && !is.null(LD_EUR_DIR)
  if (use_local) {
    bfile    <- file.path(LD_EUR_DIR, "EUR")
    bim_file <- paste0(bfile, ".bim")
    if (file.exists(bim_file)) {
      cat(sprintf("    [LD local] %d SNPs, chr%s\n", length(rsids), gene_chrom))
      ld <- tryCatch(
        ieugwasr::ld_matrix(rsids, bfile = bfile, plink_bin = PLINK_BIN),
        error = function(e) fail_loud_data(sprintf("26 eGFR SuSiE local LD (PLINK) FAILED: %s -- LD infrastructure failure, not a biological skip", conditionMessage(e)))  # 2026-07-01 (closure B1, Codex)
      )
      if (!is.null(ld)) {
        cat(sprintf("    [LD local] %d x %d matrix\n", nrow(ld), ncol(ld)))
        return(ld)
      }
      cat("    [LD local] Failed (local-only; NO API fallback by design)\n"); return(NULL)
    }
    fail_loud_data(sprintf("26 eGFR SuSiE: LD reference EUR.bim NOT FOUND (%s) -- infrastructure failure; cannot compute LD for ANY pair, not a per-pair skip", bim_file))   # 2026-07-01 Codex-p3rest #4
  }
  # No API fallback (localized 2026-06-05, mirrors phase3/15): require local PLINK + 1KG EUR.
  stop("eGFR SuSiE LD requires local PLINK_BIN + LD_EUR_DIR/EUR.bim (no OpenGWAS API).")
}

# ── Run SuSiE-coloc for one eGFR pair ───────────────────────────────────────
run_susie_egfr_pair <- function(gene, gene_chrom, gene_start, gene_end, pp4_abf) {
  region_start <- max(1L, gene_start - CIS_WINDOW_KB * 1000L)
  region_end   <- gene_end + CIS_WINDOW_KB * 1000L

  cat(sprintf("  susie: %s ~ eGFR | chr%s:%d-%d\n", gene, gene_chrom, region_start, region_end))
  flush(stdout())

  # 1. pQTL data
  pqtl_dt <- get_pqtl_region(gene, as.character(gene_chrom), region_start, region_end)
  if (is.null(pqtl_dt) || nrow(pqtl_dt) == 0) {
    # A DATA-FOUNDATION FAILURE, NOT A LOW-SNP SKIP. R/coloc.R:963-969 stops on an empty extraction
    # and on missing required columns; these copies returned NULL, which the caller renders as an NA
    # row or an omitted pair -- so a broken archive or a schema change looked like an uncomputable
    # biological result. fail_loud_data() already exists here and the main loop turns it into an abort.
    fail_loud_data(sprintf("26: gene %s pQTL extraction returned NO data (tar resolved but empty) -- data-foundation failure; mirrors R/coloc.R:963", gene))
  }
  req_cols <- c("ALLELE0", "ALLELE1", "BETA", "SE", "INFO", "N")
  if (any(!req_cols %in% names(pqtl_dt))) {
    fail_loud_data(sprintf("26: gene %s pQTL missing required columns: %s -- data-format failure; mirrors R/coloc.R:969", gene, paste(setdiff(req_cols, names(pqtl_dt)), collapse = ", ")))
  }
  pqtl_dt <- pqtl_dt[!is.na(INFO) & INFO >= INFO_MIN & !is.na(BETA) & !is.na(SE)]
  if (nrow(pqtl_dt) < 30) {
    cat(sprintf("    [skip] pQTL only %d SNPs\n", nrow(pqtl_dt))); return(NULL)
  }
  pqtl_dt[, pqtl_varbeta := SE^2]
  n_pqtl <- as.integer(pqtl_dt[["N"]][1])

  # PROJECTED TO GRCh37. pqtl_dt$chrpos is chrom:GENPOS from the UKB-PPP archive and GENPOS is
  # GRCh38, while this script reads CKDGen_eGFR_EA_hg37.txt.gz (Pos_b37). Comparing them directly
  # does not merely fail to match -- it matches the WRONG variants wherever the two builds happen to
  # agree numerically, which on the delivered eGFR analysis produced 25 such pairings.
  pqtl_37 <- liftover_to_hg37(pqtl_dt)
  .psrc <- if (!is.null(pqtl_37)) pqtl_37 else pqtl_dt
  .cp37 <- if (!is.null(pqtl_37)) pqtl_37$chrpos37 else rep(NA_character_, nrow(pqtl_dt))
  # EVERY FIELD FROM THE SAME TABLE. rsid and chrpos_pqtl were taken from the lifted table while
  # beta, varbeta and the alleles were taken from pqtl_dt -- and liftover_to_hg37() ends in
  # merge(..., by.x = "chrpos"), which REORDERS: chrpos is a character key, so "100190" sorts before
  # "10572" and the two tables share no row order at all. Every row therefore carried one variant's
  # identifier and position with another variant's effect and alleles. For ACP1, a current eGFR
  # suggestive pair, all 2,786 post-filter ranks differ between the two orders. 05_coloc_egfr.R takes
  # all six from .psrc and was never affected; this copy is where they were mixed.
  pqtl_sub <- data.table(
    rsid = .psrc$rsid, chrpos_pqtl = .cp37,
    pqtl_beta = .psrc$BETA, pqtl_varbeta = .psrc$pqtl_varbeta,
    ea_pqtl = toupper(.psrc$ALLELE1), oa_pqtl = toupper(.psrc$ALLELE0)
  )

  # 2. eGFR data
  # THE WINDOW IS PROJECTED TOO. region_start/region_end are GRCh38, from the GRCh38 gene body, while
  # get_egfr_region() filters Pos_b37 -- the comment above says exactly that about the variants and
  # then handed the bounds over unconverted. Same edit as 05_coloc_egfr.R, which this script copies.
  # The endpoints are not lifted (they can land in different chain blocks or out of order); the
  # bounds come from the projected positions of the pQTLs actually being colocalised.
  # AN OFF-CHROMOSOME LIFT IS A PROPERTY OF ONE VARIANT, NOT OF THE PAIR -- the engine's policy
  # (R/coloc.R:1049), which this script had diverged from. Aborting because a single pQTL landed on
  # another chromosome treats a legitimate chain outcome as a data-foundation failure: that variant
  # has no GRCh37 counterpart usable against this outcome, so it is dropped and the pair proceeds on
  # the rest. GUSB has two such variants among ~2,860, and the old condition stopped the whole run
  # there. What IS fatal is having NO variant left on this chromosome -- then there is no window, and
  # falling back to the GRCh38 bounds would query the wrong region silently. The coordinate is
  # cleared rather than merely excluded from the bounds, so the variant also leaves the rsID merge
  # and the position fallback.
  .chr37  <- if (!is.null(pqtl_37)) sub(":.*", "", pqtl_37$chrpos37) else NA_character_
  .same37 <- .chr37 %in% as.character(gene_chrom)
  .off37  <- !is.na(.chr37) & !.same37
  if (!any(.same37))
    fail_loud_data(sprintf(paste0("SuSiE eGFR coloc %s: the pQTL liftover to GRCh37 leaves NO variant ",
                                  "on chr%s (mapped=%d, off-chromosome=%d, unmapped=%d), so there is ",
                                  "no window. Falling back to the GRCh38 bounds would query the wrong ",
                                  "region silently."),
                           gene, gene_chrom, sum(.same37), sum(.off37), sum(is.na(.chr37))))
  if (any(.off37)) {
    cat(sprintf("    [liftover] %d pQTL(s) lifted off chr%s -- dropped for this outcome\n",
                sum(.off37), gene_chrom))
    pqtl_37$chrpos37[.off37] <- NA_character_
  }
  .pos37 <- as.integer(sub(".*:", "", pqtl_37$chrpos37[.same37]))
  out_dt <- get_egfr_region(as.character(gene_chrom), min(.pos37), max(.pos37))
  if (is.null(out_dt) || nrow(out_dt) < 20) {
    cat("    [skip] eGFR no data\n"); return(NULL)
  }

  # n_total_sum CARRIED THROUGH, because susie_rss needs the N of the variants it is actually
  # fitting. See the regional-N block before susie_rss(out) below.
  out_sub <- out_dt[, .(rsid_out, chrpos_out,
    ea_out = toupper(Allele1), oa_out = toupper(Allele2),
    out_beta = as.numeric(Effect), out_varbeta = as.numeric(StdErr)^2,
    out_n = if ("n_total_sum" %in% names(out_dt)) as.numeric(n_total_sum) else NA_real_)]
  out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

  # 3. Merge + harmonise
  common <- rsid_pos_merge(pqtl_sub, out_sub)
  cat(sprintf("    merged: %d SNPs\n", nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] < 20 merged\n"); return(NULL) }

  # HARMONISE FIRST, THEN DROP WHAT IS STILL AMBIGUOUS -- same edit and same reason as R/coloc.R.
  # A position matching two outcome records had both deleted before the alleles were looked at.
  common <- harmonize_alleles(common, "ea_pqtl", "oa_pqtl", "ea_out", "oa_out", "out_beta")
  # A REUSED RSID IS NOT AN AMBIGUOUS VARIANT -- same edit and same reason as R/coloc.R. Only the
  # colliding rows are re-keyed, so unambiguous rsIDs keep their names for the rsID-based LD lookup.
  common[duplicated(snp_id) | duplicated(snp_id, fromLast = TRUE),
         snp_id := paste(chrpos_pqtl, ea_pqtl, oa_pqtl, sep = ":")]
  common <- common[, if (.N == 1L) .SD, by = snp_id]
  cat(sprintf("    after harmonise: %d SNPs\n", nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] < 20 after harmonise\n"); return(NULL) }

  # 4. LD matrix
  common[, z_pqtl_abs := abs(pqtl_beta / sqrt(pqtl_varbeta))]
  common <- common[order(-z_pqtl_abs)]

  snps_for_ld <- common$snp_id[grepl("^rs[0-9]+$", common$snp_id)]
  cat(sprintf("    rsIDs for LD: %d / %d\n", length(snps_for_ld), nrow(common)))
  if (length(snps_for_ld) < 10) {
    cat("    [skip] < 10 rsIDs\n"); return(NULL)
  }

  ld_mat <- get_ld_matrix_auto(snps_for_ld, gene_chrom)
  if (is.null(ld_mat)) { cat("    [skip] LD failed\n"); return(NULL) }

  # 5. Align to LD matrix
  ld_snps_raw <- rownames(ld_mat)
  ld_snps_rs  <- sub("_[A-Za-z]+_[A-Za-z]+$", "", ld_snps_raw)
  idx_in_ld   <- match(ld_snps_rs, common$snp_id)
  keep_ld_pos <- which(!is.na(idx_in_ld))

  if (length(keep_ld_pos) < 10) {
    cat(sprintf("    [skip] Only %d matched LD\n", length(keep_ld_pos))); return(NULL)
  }

  common_ld  <- common[idx_in_ld[keep_ld_pos]]
  ld_mat_f   <- ld_mat[keep_ld_pos, keep_ld_pos, drop = FALSE]
  rownames(ld_mat_f) <- colnames(ld_mat_f) <- ld_snps_rs[keep_ld_pos]

  cat(sprintf("    LD-aligned: %d SNPs\n", nrow(common_ld)))
  if (nrow(common_ld) < 10) { cat("    [skip] < 10 after LD\n"); return(NULL) }

  # 5b. STRAND-ALIGN R to the pQTL effect allele (root-cause fix 2026-07-17, editor 2b). Re-orients R where
  # the .bim A1 (from the allele-suffixed LD rownames) is the pQTL OTHER allele; no-op for aligned SNPs;
  # palindromes left as-is (verified by frequency); allele_mismatch dropped. (eGFR outcome is fixed.)
  .sa <- strand_align_R(common_ld, ld_snps_raw[keep_ld_pos], ld_mat_f, label = paste0(gene, "~eGFR"))
  common_ld <- .sa$common_ld; ld_mat_f <- .sa$R
  # WHICH VARIANTS THIS FIT DROPPED, CARRIED INTO THE RESULT ROW. strand_align_R returns the rsIDs it
  # removed, and this engine used to discard that log -- so a fit that silently lost a SNP looked
  # identical to one that did not, and the only trace was stdout. That matters now because
  # 39_palindromic_freq_gate can require a SWAPPED palindrome to be removed via
  # options(ckm.palindrome_drop = ...), i.e. a deliberate, operator-specified removal in a production
  # run; without this column the published fit table would carry no record that it happened.
  # "" means nothing was dropped -- NOT NA, which would be indistinguishable from "this engine does
  # not record it" (which is what the non-eGFR half written by 15_susie_coloc.R still looks like).
  .sa_dropped <- paste(as.character(.sa$log$dropped_rsids), collapse = ",")
  if (nrow(common_ld) < 10) { cat("    [skip] < 10 after strand-align drop\n"); return(NULL) }   # re-check gate after drops (R2 Medium)

  # 6. SuSiE
  z_exp <- common_ld$pqtl_beta / sqrt(common_ld$pqtl_varbeta)
  z_out <- common_ld$out_beta  / sqrt(common_ld$out_varbeta)

  # THE N OF THE VARIANTS BEING FITTED, not of the whole file. EGFR_N is the median n_total_sum over
  # the entire CKDGen file (525,153); the fitted loci differ from it and from each other -- SDCCAG8
  # 518,470, ADAM15 524,509, GSTA3 540,514, HCLS1 565,862 -- so the SuSiE and ABF paths were using
  # different sample sizes for the same pair. N is inert for coloc.abf given beta, varbeta and sdY,
  # but susie_rss uses it. The repository's own diagnostic (35_dashed_ld_diagnostic.R:34) already
  # calls the whole-genome scalar misspecified. Falls back to the file-wide value only if the column
  # is absent, which get_egfr_region() would have to stop reporting for that to happen.
  .n_out_region <- if (all(is.na(common_ld$out_n))) EGFR_N else
                     as.integer(median(common_ld$out_n, na.rm = TRUE))
  if (.n_out_region != EGFR_N)
    cat(sprintf("    outcome N: %d (regional median over %d fitted SNPs; file-wide %d)\n",
                .n_out_region, sum(!is.na(common_ld$out_n)), EGFR_N))

  S1 <- tryCatch(
    # max_iter = SUSIE_MAX_ITER (default 1000), kept in step with the non-eGFR engine
    # (15_susie_coloc.R) -- two engines at
    # different ceilings would make "the ceiling" a per-outcome property, which it is not.
    # CEILING HISTORY: susieR's default is 100, a package default this study never chose and never
    # stated. Raised to 1000 on 2026-07-16, to 10000 on 2026-08-08.
    #
    # WHY THE CEILING IS THIS STUDY'S TO SET. coloc::runsusie, the wrapper the package provides for
    # exactly this task, defaults to repeat_until_convergence = TRUE, because its authors treat
    # non-convergence at the default limit as a computational artifact rather than a result. This
    # engine bypasses that wrapper and calls susie_rss directly, so the stopping rule is ours.
    #
    # IT IS FORWARDED, NOT SWALLOWED -- and that had to be checked, because `check_LD` on this same
    # line is NOT: max_iter is absent from formals(susie_rss) exactly as check_LD is, but
    # susie_suff_stat (which susie_rss forwards `...` to) HAS max_iter as a formal and check_LD is a
    # formal of neither. Verified empirically 2026-07-16.
    #
    # RAISING IT IS FREE FOR FITS THAT ALREADY CONVERGE. IBSS stops at its own niter, not at the
    # ceiling -- measured 2026-08-08, RTBDN-T2D exposure converged at niter = 13 under a ceiling of
    # 100. No eGFR pair failed to converge at 1000, so this raise is not expected to move any eGFR
    # result; it is here for engine parity, not for an outcome.
    susieR::susie_rss(z = z_exp, R = ld_mat_f, n = n_pqtl,
                       L = 10, min_abs_corr = 0, verbose = FALSE, check_LD = FALSE, max_iter = SUSIE_MAX_ITER),
    error = function(e) fail_loud_data(sprintf("26 eGFR SuSiE susie_rss(exp) THREW an error (malformed data/LD, not a legit no-CS result): %s", conditionMessage(e)))  # 2026-07-01 (closure B1, Codex)
  )
  if (is.null(S1)) { cat("    [skip] susie_rss exp failed\n"); return(NULL) }

  S2 <- tryCatch(
    susieR::susie_rss(z = z_out, R = ld_mat_f, n = .n_out_region,
                       L = 10, min_abs_corr = 0, verbose = FALSE, check_LD = FALSE, max_iter = SUSIE_MAX_ITER),
    error = function(e) fail_loud_data(sprintf("26 eGFR SuSiE susie_rss(out) THREW an error (malformed data/LD, not a legit no-CS result): %s", conditionMessage(e)))  # 2026-07-01 (closure B1, Codex)
  )
  if (is.null(S2)) { cat("    [skip] susie_rss out failed\n"); return(NULL) }

  n_cs_exp <- length(S1$sets$cs)
  n_cs_out <- length(S2$sets$cs)
  cat(sprintf("    SuSiE CS: exp=%d, out=%d\n", n_cs_exp, n_cs_out))

  # 7. coloc.susie
  cs_res <- tryCatch(
    coloc::coloc.susie(S1, S2),
    error = function(e) { message("    [coloc.susie] ", conditionMessage(e)); NULL }
  )
  # THE FITS RAN, SO THEIR CONVERGENCE IS KNOWN AND MUST SURVIVE THIS EXIT. coloc.susie failing (or
  # returning an empty summary) is a failure of the COLOCALIZATION, not of the SuSiE fits -- S1 and S2 are
  # complete by this point and their `$converged` is the answer this change exists to record. Returning a
  # bare NULL here sends the pair to the caller's no-fit fallback, whose row has no convergence columns, so
  # a fit whose convergence we KNOW would be written as if we did not. That path is not hypothetical: it is
  # exactly the `dropped:collapsed` class (4 of the 17 dashed candidates). rbindlist(fill=TRUE) then hides
  # it by filling NA. Return the row with the posteriors NA and the convergence carried.
  if (is.null(cs_res)) {
    # THE CONSOLE TEXT IS A MACHINE-READABLE CONTRACT, NOT A LABEL. 29 captures these lines and
    # _candidate_null_classify.R matches them: CAND_FATAL_SIG on "coloc.susie failed" (:36) and the
    # benign tag on "empty summary" (:45), case-sensitively. Rewording them to say what the row now
    # carries broke BOTH matches, turning an expected benign collapse into an unrecognised NULL: 29
    # aborts, 31 records run_status=error, 32 aborts, and the dashed layer produces nothing. The
    # wording is restored verbatim; the distinguishing information lives in the DATA below.
    # TWO DIFFERENT EVENTS MUST NOT WRITE THE SAME ROW. This branch is coloc.susie THROWING; the
    # branch below is coloc.susie succeeding and returning an empty summary. Both wrote NA for
    # n_signal_pairs, so the artifact could not tell a computation failure from a genuine zero --
    # and both are real, expected outcomes here (the dropped:collapsed class). No column is added:
    # NA and 0 are simply the honest values. Nothing was counted on this path, so it stays NA.
    cat("    [skip] coloc.susie failed\n")   # EXACT STRING -- _candidate_null_classify.R:36 matches it
    return(data.table(
      gene = gene, outcome_name = "eGFR", outcome_id = EGFR_ID,
      PP.H4_abf = as.numeric(pp4_abf), PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      concordant_susie = NA,
      n_snps_merged = nrow(common), n_snps_ld = nrow(common_ld),
      n_signals_exp = n_cs_exp, n_signals_out = n_cs_out, n_signal_pairs = NA_integer_,
      converged_exp = isTRUE(S1$converged), converged_out = isTRUE(S2$converged),
      max_iter_used = SUSIE_MAX_ITER,
      palindrome_dropped_rsids = .sa_dropped))
  }

  # THE PACKAGE'S NO-RESULT SENTINEL IS A BARE data.table(nsnps = NA), WITH NO `summary` COMPONENT.
  # Verified against the installed coloc 5.2.3: coloc.susie() has `return(data.table(nsnps = NA))` on
  # its no-credible-set path. A guard that treats any missing `summary` as fatal therefore aborts a
  # legitimate result -- the defect this guard was written to prevent, committed inside the guard.
  # So the sentinel is recognised first and takes the empty-result path; only a shape that is neither
  # the sentinel nor a summary-bearing object is a schema change and stops.
  .cs_sentinel <- is.data.frame(cs_res) && identical(names(cs_res), "nsnps")
  if (!.cs_sentinel && !("summary" %in% names(cs_res)))
    fail_loud_data(sprintf("26 eGFR SuSiE %s: coloc.susie returned an object that is neither its no-CS sentinel nor a summary-bearing result (class %s, names %s) -- the result schema changed; treating it as an empty colocalisation would discard whatever it does contain.", gene, class(cs_res)[1], paste(names(cs_res), collapse = ", ")))
  summary_dt <- as.data.table(cs_res$summary)
  # A MISSING COLUMN IS NOT A ZERO. This arm also accepted a NON-empty summary that simply lacks
  # PP.H4.abf -- a coloc package or schema change -- and recorded it with n_signal_pairs = 0 and
  # absent posteriors, identical to a genuine empty result, while discarding real signal-pair
  # rows that may have included a supportive H4. The two are separated: zero rows is a count,
  # a summary without the column is a data-format failure and stops. The message deliberately
  # avoids the strings _candidate_null_classify.R matches ("empty summary", "coloc.susie failed"),
  # because this is neither of those events.
  if (nrow(summary_dt) > 0 && !("PP.H4.abf" %in% names(summary_dt)))
    fail_loud_data(sprintf("26 eGFR SuSiE %s: coloc.susie returned %d summary row(s) with no PP.H4.abf column (present: %s) -- the summary schema changed; treating this as a null colocalisation would discard real signal pairs.",
                           gene, nrow(summary_dt), paste(names(summary_dt), collapse = ", ")))
  if (nrow(summary_dt) == 0) {
    # coloc.susie RAN and found no signal pair. That is a count of zero, not an unknown -- see the
    # note in the throw branch above.
    cat("    [skip] coloc.susie returned empty summary\n")   # EXACT STRING -- :45 matches "empty summary", case-sensitively
    return(data.table(
      gene = gene, outcome_name = "eGFR", outcome_id = EGFR_ID,
      PP.H4_abf = as.numeric(pp4_abf), PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      concordant_susie = NA,
      n_snps_merged = nrow(common), n_snps_ld = nrow(common_ld),
      n_signals_exp = n_cs_exp, n_signals_out = n_cs_out, n_signal_pairs = 0L,
      converged_exp = isTRUE(S1$converged), converged_out = isTRUE(S2$converged),
      max_iter_used = SUSIE_MAX_ITER,
      palindrome_dropped_rsids = .sa_dropped))
  }

  best <- summary_dt[which.max(PP.H4.abf)]
  pp4_susie <- best$PP.H4.abf
  pp3_susie <- best$PP.H3.abf

  cat(sprintf("    coloc.susie: %d pairs | PP.H4=%.3f PP.H3=%.3f\n",
              nrow(summary_dt), pp4_susie, pp3_susie))

    # Convergence, PER FIT AND PER SIDE. susie_rss returns `$converged` directly, so this needs no warning
    # capture. It has to reach the ARTIFACT, not just the row: the first version of this change added the
    # columns to the returned table and the col_order restriction below silently dropped them again --
    # applied but inert, which is the whole failure mode this recording exists to end.
  # NOT delegated: 26 has its own S1/S2 and its own return table. Without these two lines the eGFR fits --
  # which are the ONLY fits with a recorded non-convergence in this study, and the reason this change was
  # made -- would silently carry no convergence record.
  data.table(
    gene = gene, outcome_name = "eGFR", outcome_id = EGFR_ID,
    # STORED UNROUNDED. concordant_susie below is computed from the full-precision value, but the
    # downstream candidate layer re-applies `PP.H4_susie >= 0.80` to the STORED column -- so a raw
    # 0.79996 was written as 0.8000 and called positive downstream while this row said FALSE. Same
    # class as the abf rounding in R/coloc.R; rounding for display belongs to the reporting layer.
    PP.H4_abf = as.numeric(pp4_abf),
    PP.H4_susie = as.numeric(pp4_susie),
    PP.H3_susie = as.numeric(pp3_susie),
    # CONVERGENCE GATES THE FLAG, same change as 15_susie_coloc.R. The two `converged_*` columns below
    # were recorded and read by nothing here, so a nonconverged fit's posterior was thresholded as an
    # ordinary result. No current eGFR pair fails to converge at the 1000 ceiling (08's header records
    # that), so this changes no current eGFR value -- it removes the path by which a future one would
    # be published as a negative. NA, not FALSE: this column is already NA on the two skip paths above.
    concordant_susie = if (isTRUE(S1$converged) && isTRUE(S2$converged)) pp4_susie >= 0.80 else NA,
    n_snps_merged = nrow(common),
    n_snps_ld = nrow(common_ld),
    n_signals_exp = n_cs_exp,
    n_signals_out = n_cs_out,
    n_signal_pairs = nrow(summary_dt),
    max_iter_used = SUSIE_MAX_ITER,   # see 15_susie_coloc.R: the ceiling this fit ran under
    converged_exp = isTRUE(S1$converged),
    converged_out = isTRUE(S2$converged),
    palindrome_dropped_rsids = .sa_dropped
  )
}

# ── Main loop ────────────────────────────────────────────────────────────────
results <- list()
for (i in seq_len(nrow(strong_egfr))) {
  g     <- strong_egfr$gene[i]
  pp4   <- strong_egfr$PP.H4[i]

  # THE ROW COUNT IS TAKEN BEFORE , NOT AFTER, as in R/coloc.R:1321.  on an
  # absent gene gives zero rows, but  on a zero-row data.table returns ONE row of NA -- so
  # nrow(gene_info) was never 0 and the skip below was unreachable. The pair went on with NA
  # chrom/start/end. Every current suggestive gene is in the instrument table, so this is latent.
  .gi_all   <- instruments[gene == g]
  gene_info <- if (nrow(.gi_all)) .gi_all[1] else .gi_all
  if (nrow(gene_info) == 0) {
    cat(sprintf("[%d/%d] %s: no instrument\n", i, nrow(strong_egfr), g))
    results[[g]] <- data.table(
      gene = g, outcome_name = "eGFR", outcome_id = EGFR_ID,
      # unrounded on EVERY path: the same column was rounded only on the fallback rows.
      # THE TWO ASSIGNMENTS BELOW WERE INSIDE THAT TRAILING COMMENT, so both fallback rows lacked
      # PP.H4_susie and PP.H3_susie entirely and rbindlist(fill = TRUE) filled them in silently.
      PP.H4_abf = as.numeric(pp4),
      PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      concordant_susie = NA, n_snps_merged = NA_integer_, n_snps_ld = NA_integer_,
      n_signals_exp = NA_integer_, n_signals_out = NA_integer_, n_signal_pairs = NA_integer_,
      # NA, not FALSE -- no fit ran on this path; see 15's identical fallback.
      converged_exp = NA, converged_out = NA, max_iter_used = SUSIE_MAX_ITER
    )
    next
  }

  cat(sprintf("\n[%d/%d] %s ~ eGFR (PP.H4_abf=%.3f)\n", i, nrow(strong_egfr), g, pp4))
  flush(stdout())

  res <- tryCatch(
    run_susie_egfr_pair(g, gene_info$gene_chrom, gene_info$gene_start,
                        gene_info$gene_end, pp4),
    error = function(e) {  # 2026-07-01 (closure B1): ANY error on a strong pair is fatal (99 rejects NA); flag all, not only fail_loud_data() paths. Legit-empty skips return NULL without erroring.
      message(sprintf("  [ERROR] %s: %s", g, conditionMessage(e)))
      assign(".CKM_DATA_FAILURE", sprintf("%s ~ eGFR SuSiE: %s", g, conditionMessage(e)), envir = .GlobalEnv); NULL
    }
  )
  if (exists(".CKM_DATA_FAILURE", envir = .GlobalEnv))  # abort, do NOT record a NA SuSiE row
    stop(paste("26 ABORT -- data-foundation failure during eGFR SuSiE:", get(".CKM_DATA_FAILURE", envir = .GlobalEnv)))

  if (!is.null(res)) {
    results[[g]] <- res
  } else {
    results[[g]] <- data.table(
      gene = g, outcome_name = "eGFR", outcome_id = EGFR_ID,
      # unrounded on EVERY path: the same column was rounded only on the fallback rows.
      # THE TWO ASSIGNMENTS BELOW WERE INSIDE THAT TRAILING COMMENT, so both fallback rows lacked
      # PP.H4_susie and PP.H3_susie entirely and rbindlist(fill = TRUE) filled them in silently.
      PP.H4_abf = as.numeric(pp4),
      PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      concordant_susie = NA, n_snps_merged = NA_integer_, n_snps_ld = NA_integer_,
      n_signals_exp = NA_integer_, n_signals_out = NA_integer_, n_signal_pairs = NA_integer_,
      # 26 has TWO fallbacks; the other one (no instruments) already carried these and this one did not,
      # so the schema differed between two rows of the same table and rbindlist(fill=TRUE) filled the gap
      # silently. NA, not FALSE: no fit ran on this path.
      converged_exp = NA, converged_out = NA, max_iter_used = SUSIE_MAX_ITER
    )
  }
}

# ── Combine eGFR results ────────────────────────────────────────────────────
egfr_dt <- rbindlist(results, fill = TRUE, use.names = TRUE)
col_order <- c("gene", "outcome_name", "outcome_id",
               "PP.H4_abf", "PP.H4_susie", "PP.H3_susie",
               "concordant_susie", "n_snps_merged", "n_snps_ld",
               "n_signals_exp", "n_signals_out", "n_signal_pairs",
               "converged_exp", "converged_out",
               "max_iter_used",
               "palindrome_dropped_rsids")   # see the return table above; omit any of them and it is inert
col_order <- intersect(col_order, names(egfr_dt))
egfr_dt <- egfr_dt[, ..col_order][order(-PP.H4_abf)]

cat(sprintf("\neGFR SuSiE results: %d pairs\n", nrow(egfr_dt)))
cat(sprintf("With result: %d\n", sum(!is.na(egfr_dt$PP.H4_susie))))
cat(sprintf("Concordant (PP.H4_susie >= 0.80): %d\n",
    sum(egfr_dt$concordant_susie, na.rm = TRUE)))

if (exists("SUSIE_LIMIT")) {
  cat("\n[SMOKE] SUSIE_LIMIT set: eGFR SuSiE ran OK; skipping merge/save.\n")
  print(egfr_dt[, intersect(c("gene","PP.H4_abf","PP.H4_susie","concordant_susie"), names(egfr_dt)), with = FALSE])
  q(save = "no", status = 0)
}

# ── Merge with existing non-eGFR results ─────────────────────────────────────
existing <- readRDS(file.path(out_dir, "susie_coloc_results.rds"))
setDT(existing)

# THE MERGE CANNOT TELL A STALE SCHEMA FROM A MISSING VALUE, SO CHECK BEFORE IT. `rbind(..., fill = TRUE)`
# below back-fills any column `existing` lacks with NA. If 07_susie (engine 15) has not been re-run since
# convergence recording was added, every non-eGFR row silently acquires `converged_* = NA` -- asserting
# "we do not know whether these fits converged" about fits that ran and whose answer we would have had.
# Half the artifact would then be new and half stale, with no error anywhere. The execution order in
# REBUILD_LEDGER 2.5a (07 -> 08) is what prevents this; nothing enforced it until now.
# max_iter_used JOINS THE PRESENCE CHECK (review finding, 2026-08-08). Without it this script
# will merge a non-eGFR half produced under a different iteration ceiling and the combined table
# records no trace of the mix: 60 rows at one ceiling beside 16 at another, every column valid.
.need_conv <- c("converged_exp", "converged_out", "max_iter_used")
.miss_conv <- setdiff(.need_conv, names(existing))
if (length(.miss_conv))
  # THE CAUSE IS NAMED FROM WHAT IS ACTUALLY MISSING (2026-08-08). This guard now fires for two
  # different reasons and the message asserted only the older one: a half that has both convergence
  # columns and lacks only max_iter_used does NOT predate the convergence recording, and the damage
  # is not "convergence unknown". Naming the wrong cause sends reconciliation to look in the wrong
  # place -- the same failure 03b_gene_coords_hg38.R:731 records for its own four-cause refusal.
  # The remedy is identical either way, so only the diagnosis was wrong, not the instruction.
  stop(sprintf("26: susie_coloc_results.rds is missing %s. %s rbind(fill=TRUE) would write NA into %s for every existing non-eGFR row, so the merged table could not certify that both halves share a convergence record and an iteration ceiling. Re-run 07_susie FIRST, then this script (REBUILD_LEDGER 2.5a gives the order).",
               paste(.miss_conv, collapse = ", "),
               if (all(c("converged_exp", "converged_out") %in% .miss_conv))
                 "It predates the per-fit convergence recording, so 07_susie (engine 15) has not been re-run since."
               else if (identical(.miss_conv, "max_iter_used"))
                 "It carries the convergence columns but predates the iteration-ceiling recording, so it was produced by an earlier 07_susie run under an unrecorded ceiling."
               else
                 "It was produced by an earlier 07_susie run than this one.",
               paste(.miss_conv, collapse = "/")))
# PRESENCE IS NOT CORRECTNESS. A file that was ALREADY tainted by an earlier fill=TRUE HAS both columns --
# every value NA -- and a presence check waves it through. Check the premise instead: a row that shows
# evidence the fit RAN (a posterior, or an LD SNP count) cannot honestly say "convergence unknown". This is
# the failure mode this layer has a name for -- verifying a label exists rather than that its premise holds
# -- committed inside the guard written to fix that very class.
# `!is.na(PP.H4_susie)` was the wrong evidence and misses the rows this recording was extended for: a
# post-fit collapse legitimately has NA posterior while its fits ran and their convergence is known.
# susie_fit_ran() (R/00_config.R) is the single definition; see its header for why.
assert_convergence_honest(existing, "26: susie_coloc_results.rds")

# THE EXISTING NON-eGFR HALF MUST BE THIS GENERATION'S. The convergence check above only proves
# the file was written by a version that records convergence -- an earlier same-schema file passes.
# The realistic trigger is an interrupted rerun: 07 fails before replacing the previous
# susie_coloc_results.rds, then 08/26 is resumed and appends this generation's eGFR rows to the
# previous generation's non-eGFR rows. Every row then looks ordinary and nothing records the
# mixture.
#
# Checked SEMANTICALLY, not by byte fingerprint: the non-eGFR keys must be exactly this
# generation's strong non-eGFR pairs, and their PP.H4_abf must be the values coloc_final carries.
# That is what the withdrawn input-key comparison was reaching for, without being defeated by a
# cat() label edit (see the note at the provenance check above).
.ne_have <- unique(existing[outcome_id != EGFR_ID, .(gene, outcome_id)])
.ne_want <- unique(as.data.table(coloc_final)[is_coloc_pass(PP.H4, PP.H3) & outcome_id != EGFR_ID,
                                              .(gene, outcome_id)])
.ne_missing <- fsetdiff(.ne_want, .ne_have)
.ne_extra   <- fsetdiff(.ne_have, .ne_want)
if (nrow(.ne_missing) || nrow(.ne_extra))
  stop(sprintf("26 eGFR SuSiE: the non-eGFR half of %s is not this generation's strong set (%d expected pair(s) absent, %d unexpected present).%s  absent: %s%s  unexpected: %s%sThe two halves would be from different runs. Re-run 07_susie before this script (REBUILD_LEDGER 2.5a gives the order).",
               basename(file.path(out_dir, "susie_coloc_results.rds")), nrow(.ne_missing), nrow(.ne_extra), "\n",
               paste(utils::head(.ne_missing[, paste(gene, outcome_id)], 6), collapse = ", "), "\n",
               paste(utils::head(.ne_extra[, paste(gene, outcome_id)], 6), collapse = ", "), "\n"))
.ne_chk <- merge(existing[outcome_id != EGFR_ID, .(gene, outcome_id, PP.H4_abf)],
                 as.data.table(coloc_final)[, .(gene, outcome_id, .pp = PP.H4)],
                 by = c("gene", "outcome_id"))
.ne_bad <- .ne_chk[!is.na(PP.H4_abf) & !is.na(.pp) & abs(PP.H4_abf - .pp) > 1e-6]
if (nrow(.ne_bad))
  stop(sprintf("26 eGFR SuSiE: %d non-eGFR row(s) in %s carry a PP.H4_abf that does not match coloc_results_v7 (e.g. %s: %.6f vs %.6f) -- the SuSiE table was computed against a different coloc result than the one being read now.",
               nrow(.ne_bad), "susie_coloc_results.rds", .ne_bad$gene[1], .ne_bad$PP.H4_abf[1], .ne_bad$.pp[1]))
cat(sprintf("[provenance] non-eGFR half reconciles: %d pair(s), PP.H4_abf matches coloc_results_v7\n", nrow(.ne_have)))
saveRDS(existing, file.path(out_dir, "susie_coloc_results.rds.bak_pre_egfr"))

# AN eGFR ROW THAT IS NO LONGER A STRONG PAIR MUST GO, not merely be superseded. The de-dup below
# replaces keys that appear in BOTH tables, so a rerun updates every pair the current strong set
# still contains -- and leaves behind, untouched and fully populated, any eGFR pair the previous
# run had and this one does not. The saved artifact then reports a numeric SuSiE result for a pair
# that is not a strong hit of this generation, in the same schema as the real rows and with nothing
# marking it. Only THIS script's own outcome is dropped: non-eGFR rows belong to 15 and are not
# this script's to judge.
.eg_old <- existing[outcome_id %chin% EGFR_ID]
.eg_new_keys <- unique(egfr_dt$gene)
.eg_drop <- setdiff(unique(.eg_old$gene), .eg_new_keys)
if (length(.eg_drop)) {
  cat(sprintf("[merge] dropping %d eGFR row(s) from the previous run that are not strong pairs of this generation: %s
",
              length(.eg_drop), paste(.eg_drop, collapse = ", ")))
  existing <- existing[!(outcome_id %chin% EGFR_ID & gene %chin% .eg_drop)]
}

combined <- rbind(existing, egfr_dt, fill = TRUE, use.names = TRUE)
# 2026-07-01 Codex-p3rest #1: idempotent -- if `existing` already holds eGFR rows (a rerun without clearing),
# blind rbind would DUPLICATE them; keep one row per identifying key, preferring the freshly-appended eGFR.
.keycols <- intersect(c("gene","outcome_id","outcome_name"), names(combined))
# 2026-07-01: use data.table's native `by=` for de-dup. The prior `combined[, ..keycols]` was a runtime bug --
# `..keycols` resolves to a variable named `keycols` (no dot), but the var is `.keycols` -> "not found in calling
# scope" abort (static review missed it; it only shows when the merge path runs).
if (length(.keycols)) combined <- combined[!duplicated(combined, by = .keycols, fromLast = TRUE)]
combined <- combined[order(-PP.H4_abf)]
# WRITTEN BESIDE, THEN PROMOTED (review finding, 2026-08-08). These two were saved straight to their
# final paths, RDS first and CSV second, so an interruption between them left a NEW rds beside the
# PREVIOUS csv -- and the two are read by different consumers: 26_susie_coloc_egfr.R takes the rds,
# while 13_rebuild_pair_evidence.R and 99_completeness_audit.R take the csv. Both files parse, both
# look complete, and nothing compares them, so a split-ceiling artifact is indistinguishable from a
# correct one. Writing both temporaries first means a death before the first rename leaves the two
# final files as they were -- consistent with each other -- rather than half-updated.
.rds_f <- file.path(out_dir, "susie_coloc_results.rds")
.csv_f <- file.path(out_dir, "susie_coloc_results.csv")
saveRDS(combined, paste0(.rds_f, ".tmp"))
fwrite(combined,  paste0(.csv_f, ".tmp"))
if (!all(file.rename(c(paste0(.rds_f, ".tmp"), paste0(.csv_f, ".tmp")), c(.rds_f, .csv_f))))
  stop("26 eGFR SuSiE: could not promote both result files; the pair on disk may now disagree -- inspect ",
       out_dir, " before using either.")

cat(sprintf("\nCombined: %d rows (%d non-eGFR + %d eGFR)\n",
    nrow(combined), nrow(existing), nrow(egfr_dt)))

# ── Summary ──────────────────────────────────────────────────────────────────
cat(sprintf("\n=== SuSiE-coloc Summary (All %d pairs) ===\n", nrow(combined)))
cat(sprintf("Total pairs: %d\n", nrow(combined)))
cat(sprintf("With SuSiE result: %d\n", sum(!is.na(combined$PP.H4_susie))))
cat(sprintf("Concordant (susie >= 0.80): %d\n",
    sum(combined$concordant_susie, na.rm = TRUE)))
# THE DISCORDANT COUNT EXCLUDES NONCONVERGED FITS, like concordant_susie above. It counted every
# numeric PP.H4_susie < 0.80 regardless of convergence, so NCAN-CAD (neither side converged,
# PP.H4_susie = 1.777e-48) and ADM-AF (exposure did not converge) were reported as evidence AGAINST
# colocalisation. A fit that did not converge is not a negative result, and folding the two into the
# negative count is the same error the NA above was introduced to stop. The counts now partition the
# pairs that have a fit: concordant + discordant + not-evaluable.
.conv_ok <- combined$converged_exp %in% TRUE & combined$converged_out %in% TRUE
cat(sprintf("Discordant (abf >= 0.80, susie < 0.80, both fits converged): %d\n",
    sum(.conv_ok & !is.na(combined$PP.H4_susie) & combined$PP.H4_abf >= 0.80 &
        combined$PP.H4_susie < 0.80, na.rm = TRUE)))
cat(sprintf("Not evaluable (a fit did not converge): %d\n",
    sum(!.conv_ok & !is.na(combined$PP.H4_susie), na.rm = TRUE)))

cat("\neGFR pairs:\n")
print(egfr_dt[, .(gene, PP.H4_abf, PP.H4_susie, n_signals_exp, n_signals_out, concordant_susie)])

cat("\n=== Done ===\n")
