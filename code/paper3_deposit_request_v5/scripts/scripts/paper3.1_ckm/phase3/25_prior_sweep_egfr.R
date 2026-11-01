# Phase 3-S Extension: Prior Sweep for 14 eGFR Strong Hits
# ============================================================================
# The original prior sweep (10_phase3s_prior_sweep.R) used coloc_results_v7_stage1.rds
# which has no eGFR data. This script runs the p12 sweep for the 14 eGFR
# strong hits using the local CKDGen file, then merges with existing results.
#
# Usage: cd "<E>/20260331-CKM/Paper 3" && Rscript scripts/paper3.1_ckm/phase3/25_prior_sweep_egfr.R
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
})

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_prior_bands.R"))   # is_prior_stable()/prior_band(): ONE definition, deliberately NOT in 00_config.R (coloc inputs key)
COLOC_SOURCE_ONLY <- TRUE   # define functions only; a bare source() of R/coloc.R runs main()
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R"))   # liftover_to_hg37(): ONE definition, shared
rm(COLOC_SOURCE_ONLY)

cat("=== Phase 3-S Extension: eGFR Prior Sweep ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── Settings ─────────────────────────────────────────────────────────────────
EGFR_FILE  <- file.path(OUTCOME_DIR, "CKDGen_eGFR_EA_hg37.txt.gz")
# eGFR id from config (corrected cascade = Wuttke2019_eGFR_EA); legacy fallback + warn if absent.
EGFR_ID    <- if (exists("OUTCOMES") && !is.null(OUTCOMES[["eGFR"]])) OUTCOMES[["eGFR"]] else { warning("[eGFR-id] OUTCOMES$eGFR missing in config; legacy fallback ebi-a-GCST003372"); "ebi-a-GCST003372" }
MAF_MIN    <- 0.01
PP4_THRESH <- 0.80

stopifnot(file.exists(EGFR_FILE))

if (!exists("PSW_OUT")) PSW_OUT <- file.path(RESULTS_DIR, "phase3_sensitivity")   # override for corrected cascade
out_dir <- PSW_OUT
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load rsID map ────────────────────────────────────────────────────────────
load_rsid_map <- function() {
  f <- file.path(BASE_DIR, "cache/cis_rsid_map.tsv")
  if (!file.exists(f)) { cat("[warn] cis_rsid_map.tsv not found\n"); return(NULL) }
  dt <- fread(f, colClasses = list(character = "chrom"))
  # KEYED ON THE ALLELES TOO -- the same edit as R/coloc.R, 05 and 26, whose map this copies. A
  # position-only key collapses two variants at one coordinate to one rsID: at 1:2102322 the map
  # holds T/G -> rs201135834 and T/TG -> rs34905415, and after the MAF filter only one survives, so
  # it took the other's rsID and matched the wrong outcome row.
  dt[, lookup_key := paste(chrom, pos, toupper(allele0), toupper(allele1), sep = ":")]
  dt <- dt[order(as.integer(chrom), pos)]
  # ONE RSID -> ONE POSITION, BUT EVERY ALLELE PAIR AT IT. `.SD[1], by = rsid` collapsed the
  # legacy-multiple-POSITIONS case it was written for AND the distinct-variants case it was not:
  # 1:1251285 carries G/A and G/T under one rsID (rs6657544).
  dt <- dt[, .SD[chrom == chrom[1L] & pos == pos[1L]], by = rsid]
  dt <- dt[, .SD[1], by = lookup_key]
  cat(sprintf("[rsID map] %s entries\n\n", formatC(nrow(dt), big.mark = ",")))
  dt
}
rsid_map <- load_rsid_map()

# ── Load 14 eGFR strong hits ────────────────────────────────────────────────
if (!exists("PSW_COLOC")) PSW_COLOC <- file.path(PHASE3_DIR, "coloc_results_final.rds")   # override for corrected (v6)
if (!exists("PSW_INSTR")) PSW_INSTR <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")
coloc_final <- readRDS(PSW_COLOC)
instruments <- readRDS(PSW_INSTR)
setDT(coloc_final); setDT(instruments)

strong_egfr <- coloc_final[is_coloc_pass(PP.H4, PP.H3) & outcome_id == EGFR_ID][order(-PP.H4)]   # P3-F04: centralized dual-rule is_coloc_pass() (R/00_config.R), not a re-inlined predicate (single source of truth; eGFR pairs all pass the ratio)
if (exists("PSW_LIMIT")) strong_egfr <- strong_egfr[seq_len(min(PSW_LIMIT, nrow(strong_egfr)))]   # smoke cap (scope only)
cat(sprintf("eGFR strong hits: %d pairs\n\n", nrow(strong_egfr)))

# --- optional pair-chunk parallelism (PSW_CHUNK_N>0 = worker; PSW_MERGE=1 = merge) ---
# Code version for the chunk protocol, the same three files 16 hashes: a chunk written by different
# code has the right filename, row count and columns, and passes every count-based assertion.
.psw_code_ver <- tryCatch(paste(vapply(c("R/strand_align.R", "R/coloc.R",
                                         "scripts/paper3.1_ckm/phase3/25_prior_sweep_egfr.R"),
                                       function(.f) unname(tools::md5sum(file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"), .f))),
                                       character(1)), collapse = ":"),
                          error = function(e) NA_character_)
.PCN <- as.integer(Sys.getenv("PSW_CHUNK_N","0")); .PCI <- as.integer(Sys.getenv("PSW_CHUNK_I","0"))
.PMERGE <- Sys.getenv("PSW_MERGE","") %in% c("1","true","TRUE")
.PCHUNK <- .PCN > 0 && !.PMERGE
if (.PCHUNK) { .psz <- ceiling(nrow(strong_egfr)/.PCN); .plo <- .PCI*.psz + 1; .phi <- min((.PCI+1)*.psz, nrow(strong_egfr))
  strong_egfr <- if (.plo <= nrow(strong_egfr)) strong_egfr[.plo:.phi] else strong_egfr[0]
  cat(sprintf("PSW-eGFR CHUNK %d/%d: %d pairs\n", .PCI, .PCN, nrow(strong_egfr))) }
# in MERGE mode keep strong_egfr FULL (orig_sub needs it); the loop below is skipped via its bound.

# ── Helper: allele harmonization ─────────────────────────────────────────────
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
  # SNVs ONLY in the strand branches, and the predicate comes from the engine (R/coloc.R,
  # sourced above) rather than a fourth copy. complement_allele() is chartr(), which complements
  # without reversing ("GGA" -> "CCT", not "TCC"), so a genuinely reverse-complemented indel was
  # dropped while an unrelated multi-base pair could match by coincidence.
  strand_snv <- nchar(ep) == 1L & nchar(op) == 1L & nchar(eo) == 1L & nchar(oo) == 1L
  keep       <- allele_compatible(ep, op, eo, oo)
  needs_flip <- flip_beta | (strand_snv & !palin & flip_str)
  dt <- dt[keep]
  if (nrow(dt) > 0 && any(needs_flip[keep]))
    dt[needs_flip[keep], (beta_o) := -get(beta_o)]
  dt
}

# ── Helper: pQTL region extraction ───────────────────────────────────────────
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
    fail_loud_data(sprintf("25 eGFR prior-sweep extraction FAILED (python exit %s) -- not a legit skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: no diagnostic from the process, which points at process creation rather than the script>"))
  if (length(out) <= 1) return(NULL)

  dt <- suppressWarnings(tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("25 eGFR region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e)))))  # 2026-07-01 (closure sweep): parallel of 03/06
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("25 eGFR region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")

  dt[, chrpos := paste0(chrom, ":", GENPOS)]
  if ("A1FREQ" %in% names(dt))
    dt <- dt[!is.na(A1FREQ) & A1FREQ >= MAF_MIN & A1FREQ <= 1 - MAF_MIN]

  if (!is.null(rsid_map)) {
    # Position AND alleles, written exactly as 05 and 26 write it: joining on chrpos alone would
    # reinstate the collapse the map key was widened to prevent, and a copy that computes the same
    # string a different way is where the next drift starts.
    dt[, lookup_key := paste(chrpos, toupper(ALLELE0), toupper(ALLELE1), sep = ":")]
    dt <- merge(dt, rsid_map[, .(lookup_key, rsid)],
                by = "lookup_key", all.x = TRUE)
  } else {
    dt[, rsid := NA_character_]
  }
  dt
}

# ── Helper: eGFR region extraction ───────────────────────────────────────────
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
            # RAISE, DO NOT SKIP -- the same policy as 26 and R/coloc.R. A short or unparseable row is a
            # format failure in the GWAS file, not a variant to pass over: skipping it silently thins
            # the extracted region while the R code around it advertises fail-loud extraction.
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
    fail_loud_data(sprintf("25 eGFR prior-sweep extraction FAILED (python exit %s) -- not a legit skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: no diagnostic from the process, which points at process creation rather than the script>"))
  if (length(out) <= 1) return(NULL)

  dt <- suppressWarnings(tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("25 eGFR region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e)))))  # 2026-07-01 (closure sweep): parallel of 03/06
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("25 eGFR region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")

  setnames(dt, trimws(names(dt)))
  if ("Chr" %in% names(dt)) dt[, Chr := gsub("^chr", "", as.character(Chr))]
  dt[, rsid_out := RSID]
  dt[is.na(rsid_out) | !grepl("^rs", rsid_out), rsid_out := NA_character_]
  dt[, chrpos_out := paste0(Chr, ":", Pos_b37)]
  dt
}

# ── Helper: merge by rsid + chrpos fallback ──────────────────────────────────
rsid_pos_merge <- function(pqtl_sub, out_sub) {
  result <- list()
  matched_pqtl_var <- character(0)   # (chrpos:EA:OA) of the pQTL rows the primary merge already used

  has_p_rsid <- any(!is.na(pqtl_sub$rsid))
  has_o_rsid <- any(!is.na(out_sub$rsid_out))

  if (has_p_rsid && has_o_rsid) {
    # KEYED ON (rsID, POSITION), like the engine (R/coloc.R) and 26 -- this copy keyed on the rsID
    # alone, and o_r did not even carry chrpos_out, so a variant whose rsID appears at a different
    # coordinate in the outcome was accepted as the same variant. That is precisely what the dual
    # key exists to reject, and it is reachable here because the eGFR file and the lifted pQTL
    # positions are both GRCh37 and so are directly comparable.
    p_r <- pqtl_sub[!is.na(rsid) & !is.na(chrpos_pqtl),
                    .(rsid, chrpos_pqtl, pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl)]
    o_r <- out_sub[!is.na(rsid_out) & !is.na(chrpos_out),
                   .(rsid_out, chrpos_out, out_beta, out_varbeta, ea_out, oa_out)]
    rsid_m <- merge(p_r, o_r, by.x = c("rsid", "chrpos_pqtl"), by.y = c("rsid_out", "chrpos_out"))
    if (nrow(rsid_m) > 0) {
      rsid_m[, snp_id := rsid]
      # (position, alleles), not position -- see R/coloc.R. An rsID hit that cannot harmonise used to
      # close the coordinate for the compatible record at it, and a harmonisable one closed it for every
      # OTHER pQTL variant sitting there.
      .ok <- allele_compatible(rsid_m$ea_pqtl, rsid_m$oa_pqtl, rsid_m$ea_out, rsid_m$oa_out)
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
    o_fb <- out_sub[!is.na(chrpos_out), .(chrpos_out, out_beta, out_varbeta, ea_out, oa_out)]
    pos_m <- merge(p_fb, o_fb, by.x = "chrpos_pqtl", by.y = "chrpos_out")
    if (nrow(pos_m) > 0) {
      pos_m[, snp_id := chrpos_pqtl]
      result[["pos"]] <- pos_m
    }
  }

  if (length(result) == 0) return(data.table())
  rbindlist(result, use.names = TRUE, fill = TRUE)
}

# ── Run coloc with custom p12 for one eGFR pair ─────────────────────────────
run_egfr_coloc_p12 <- function(gene, chrom, gene_start, gene_end, p12_val) {
  region_start <- max(1L, gene_start - CIS_WINDOW_KB * 1000L)
  region_end   <- gene_end + CIS_WINDOW_KB * 1000L

  pqtl_dt <- get_pqtl_region(gene, chrom, region_start, region_end)
  if (is.null(pqtl_dt) || nrow(pqtl_dt) == 0)
    fail_loud_data(sprintf("25: gene %s pQTL extraction returned NO data (tar resolved but empty) -- data-foundation failure; mirrors R/coloc.R:963", gene))

  req_cols <- c("ALLELE0", "ALLELE1", "BETA", "SE", "INFO", "N")
  # A DATA-FOUNDATION FAILURE, NOT A LOW-SNP SKIP. R/coloc.R:963-969 stops on an empty extraction
  # and on missing required columns; these copies returned NULL, which the caller renders as an NA
  # row or an omitted pair -- so a broken archive or a schema change looked like an uncomputable
  # biological result. fail_loud_data() already exists here and the main loop turns it into an abort.
  if (length(setdiff(req_cols, names(pqtl_dt))) > 0)
    fail_loud_data(sprintf("25: gene %s pQTL missing required columns: %s -- data-format failure; mirrors R/coloc.R:969", gene, paste(setdiff(req_cols, names(pqtl_dt)), collapse = ", ")))

  pqtl_dt <- pqtl_dt[!is.na(INFO) & INFO >= INFO_MIN & !is.na(BETA) & !is.na(SE)]
  if (nrow(pqtl_dt) < 30) return(NULL)

  pqtl_dt[, pqtl_varbeta := SE^2]
  n_pqtl <- as.integer(pqtl_dt[["N"]][1])

  # PROJECTED TO GRCh37. pqtl_dt$chrpos is chrom:GENPOS from the UKB-PPP archive and GENPOS is
  # GRCh38, while this script reads CKDGen_eGFR_EA_hg37.txt.gz (Pos_b37). Comparing them directly
  # does not merely fail to match -- it matches the WRONG variants wherever the two builds happen to
  # agree numerically, which on the delivered eGFR analysis produced 25 such pairings.
  pqtl_37 <- liftover_to_hg37(pqtl_dt)
  .psrc <- if (!is.null(pqtl_37)) pqtl_37 else pqtl_dt
  .cp37 <- if (!is.null(pqtl_37)) pqtl_37$chrpos37 else rep(NA_character_, nrow(pqtl_dt))
  # EVERY FIELD FROM THE SAME TABLE. rsid and chrpos_pqtl came from the lifted table while beta,
  # varbeta and the alleles came from pqtl_dt -- and liftover_to_hg37() ends in
  # merge(..., by.x = "chrpos"), which REORDERS: chrpos is a character key, so "100190" sorts before
  # "10572" and the two tables share no row order. Every row then carried one variant's identifier
  # and position with another variant's effect and alleles. For ACP1, a current eGFR suggestive pair,
  # 17 row positions change, and a targeted reconstruction produced 14 hybrid candidates of which 7
  # passed allele compatibility -- i.e. the wrong beta reached coloc. Same edit as 26.
  pqtl_sub <- data.table(
    rsid = .psrc$rsid, chrpos_pqtl = .cp37,
    pqtl_beta = .psrc$BETA, pqtl_varbeta = .psrc$pqtl_varbeta,
    ea_pqtl = toupper(.psrc$ALLELE1), oa_pqtl = toupper(.psrc$ALLELE0)
  )

  # THE WINDOW IS PROJECTED TOO, and an off-chromosome lift drops one variant rather than the pair.
  # region_start/region_end are GRCh38, from the GRCh38 gene body, while get_egfr_region() filters
  # Pos_b37 -- this script lifted the variants and then handed the bounds over unconverted, which is
  # the same build-mixing the rebuild exists to remove. For DNER the two windows overlap in only
  # 1,401 of 4,653 / 4,837 eGFR rows. The endpoints are not lifted (they can land in different chain
  # blocks or out of order); the bounds come from the projected positions of the pQTLs actually being
  # colocalised. Policy for an off-chromosome lift follows the engine (R/coloc.R:1049): clear that
  # variant's coordinate so it leaves the bounds, the rsID merge and the position fallback alike, and
  # fail only when nothing is left on this chromosome.
  .chr37  <- if (!is.null(pqtl_37)) sub(":.*", "", pqtl_37$chrpos37) else NA_character_
  .same37 <- .chr37 %in% as.character(chrom)
  .off37  <- !is.na(.chr37) & !.same37
  if (!any(.same37))
    fail_loud_data(sprintf(paste0("25 eGFR prior sweep %s: the pQTL liftover to GRCh37 leaves NO ",
                                  "variant on chr%s (mapped=%d, off-chromosome=%d, unmapped=%d), so ",
                                  "there is no window. Falling back to the GRCh38 bounds would query ",
                                  "the wrong region silently."),
                           gene, chrom, sum(.same37), sum(.off37), sum(is.na(.chr37))))
  if (any(.off37)) {
    cat(sprintf("    [liftover] %d pQTL(s) lifted off chr%s -- dropped for this outcome\n",
                sum(.off37), chrom))
    pqtl_37$chrpos37[.off37] <- NA_character_
    pqtl_sub[.off37, chrpos_pqtl := NA_character_]
  }
  .pos37 <- as.integer(sub(".*:", "", pqtl_37$chrpos37[.same37]))
  out_dt <- get_egfr_region(chrom, min(.pos37), max(.pos37))
  if (is.null(out_dt) || nrow(out_dt) < 20) return(NULL)

  out_sub <- out_dt[, .(rsid_out, chrpos_out,
    ea_out = toupper(Allele1), oa_out = toupper(Allele2),
    out_beta = as.numeric(Effect), out_varbeta = as.numeric(StdErr)^2)]
  out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

  common <- rsid_pos_merge(pqtl_sub, out_sub)
  if (nrow(common) < 20) return(NULL)

  # HARMONISE FIRST, THEN DROP WHAT IS STILL AMBIGUOUS -- the engine's order (R/coloc.R:1142), which
  # this copy had reversed. Run the other way round, a position matching two outcome records (the
  # pQTL's A/G against an A/G and an A/T, or any multiallelic site) had BOTH rows deleted before the
  # alleles were looked at, even though harmonisation would have kept the compatible one. Measured:
  # ACP1 and SMPD3 each lose one variant relative to harmonise-first. The ambiguity drop stays --
  # a group still multiple AFTER harmonisation is genuinely unresolved.
  common <- harmonize_alleles(common, "ea_pqtl", "oa_pqtl", "ea_out", "oa_out", "out_beta")
  # A REUSED RSID IS NOT AN AMBIGUOUS VARIANT. snp_id is the rsID, and two distinct variants can
  # share one (1:1251285 G/A and G/T are both rs6657544); without this they arrive as a duplicated
  # group and are both deleted although harmonisation resolved each to its own outcome row. Only the
  # colliding ones are re-keyed, so unambiguous rsIDs keep their name.
  common[duplicated(snp_id) | duplicated(snp_id, fromLast = TRUE),
         snp_id := paste(chrpos_pqtl, ea_pqtl, oa_pqtl, sep = ":")]
  common <- common[, if (.N == 1L) .SD, by = snp_id]
  if (nrow(common) < 20) return(NULL)

  if ("n_total_sum" %in% names(out_dt)) {
    out_n <- as.integer(median(out_dt$n_total_sum, na.rm = TRUE))
  } else {
    stop("25 eGFR prior sweep: n_total_sum absent in Wuttke eGFR — refusing stale global-N fallback (was 567460); check the GWAS file")  # 2026-06-30 fail-loud (gate obj #7)
  }

  ds1 <- list(beta = common$pqtl_beta, varbeta = common$pqtl_varbeta,
              N = n_pqtl, sdY = 1, type = "quant", snp = common$snp_id)
  ds2 <- list(beta = common$out_beta, varbeta = common$out_varbeta,
              N = out_n, sdY = GWAS_META[[EGFR_ID]]$sdY, type = "quant", snp = common$snp_id)   # phenotype SD from GWAS_META (R/coloc.R), not a literal 1: eGFR effects are on the ln scale and coloc's quant prior is 0.15*sdY

  cr <- tryCatch(
    coloc::coloc.abf(dataset1 = ds1, dataset2 = ds2, p12 = p12_val),
    error = function(e) fail_loud_data(sprintf("25 eGFR prior-sweep coloc.abf FAILED on valid merged data (p12=%.0e): %s -- computation failure, not a legit non-result", p12_val, conditionMessage(e)))  # 2026-07-01 (closure B1, Codex)
  )
  if (is.null(cr)) return(NULL)

  pp <- cr$summary
  data.table(
    gene = gene, outcome_id = EGFR_ID, outcome_name = "eGFR",
    p12 = p12_val,
    # NOT ROUNDED. The stability classification below thresholds these at 0.80 and 0.20, so a raw
    # 0.79996 rounded to 0.8000 changes the class. The engine returns full precision; this copy did not.
    PP.H4 = unname(pp["PP.H4.abf"]),
    PP.H3 = unname(pp["PP.H3.abf"])
  )
}

# ── p12 sweep ────────────────────────────────────────────────────────────────
p12_values <- c(1e-6, 5e-6, 1e-5)
all_results <- list()

for (i in if (.PMERGE) integer(0) else seq_len(nrow(strong_egfr))) {
  g <- strong_egfr$gene[i]
  # THE ROW COUNT IS TAKEN BEFORE , NOT AFTER, as in R/coloc.R:1321.  on an
  # absent gene gives zero rows, but  on a zero-row data.table returns ONE row of NA -- so
  # nrow(gene_info) was never 0 and the skip below was unreachable. The pair went on with NA
  # chrom/start/end. Every current suggestive gene is in the instrument table, so this is latent.
  .gi_all   <- instruments[gene == g]
  gene_info <- if (nrow(.gi_all)) .gi_all[1] else .gi_all
  if (nrow(gene_info) == 0) {
    stop(sprintf("25 eGFR prior sweep: %s is a strong eGFR pair in coloc_results_v7.rds but has NO row in the Phase-1 instrument table -- the two inputs of this generation disagree, which is not a biological skip. 16_susie_grid treats the same condition as fatal; 15 and 26 record a traceable NA row instead, but this sweep table is gene x p12 and has no such row to record, so omitting the gene here is silent.", g))
  }

  for (p12_val in p12_values) {
    cat(sprintf("[%d/%d] %s ~ eGFR, p12=%.0e ... ", i, nrow(strong_egfr), g, p12_val))
    flush(stdout())

    res <- tryCatch(
      run_egfr_coloc_p12(g, gene_info$gene_chrom, gene_info$gene_start,
                         gene_info$gene_end, p12_val),
      error = function(e) {  # 2026-07-01 (closure B1): ANY error on a strong pair is fatal (99 rejects NA); flag all, not only fail_loud_data() paths. Legit-empty skips return NULL without erroring.
        message(e$message); assign(".CKM_DATA_FAILURE", sprintf("%s ~ eGFR: %s", g, conditionMessage(e)), envir = .GlobalEnv); NULL
      }
    )
    if (exists(".CKM_DATA_FAILURE", envir = .GlobalEnv))  # abort, do NOT record a skip
      stop(paste("25 ABORT -- data-foundation failure during eGFR prior sweep:", get(".CKM_DATA_FAILURE", envir = .GlobalEnv)))

    if (!is.null(res)) {
      all_results[[paste(g, p12_val)]] <- res
      cat(sprintf("PP.H4=%.3f\n", res$PP.H4))
    } else {
      cat(sprintf("[P3-F2 WARN] eGFR prior-sweep run_coloc_pair NULL for STRONG pair %s at p12=%.1e -> NA row (investigate; not a silent skip)\n", g, p12_val))
      all_results[[paste(g, p12_val)]] <- data.table(
        gene = g, outcome_id = EGFR_ID, outcome_name = "eGFR",
        p12 = p12_val, PP.H4 = NA_real_, PP.H3 = NA_real_
      )
    }
  }
}

# ── Combine eGFR sweep ──────────────────────────────────────────────────────
egfr_sweep <- rbindlist(all_results, fill = TRUE)
if (.PCHUNK) {
  attr(egfr_sweep, "psw_code_ver") <- .psw_code_ver
  saveRDS(egfr_sweep, file.path(out_dir, sprintf("egfr_sweep_chunk%02d.rds", .PCI)))
  cat(sprintf("PSW-eGFR CHUNK %d done: %d rows\n", .PCI, nrow(egfr_sweep))); quit(save = "no")
}
if (.PMERGE) {
  cf <- sort(Sys.glob(file.path(out_dir, "egfr_sweep_chunk*.rds")))
  # EVERY CHUNK, OR NONE. This took whatever files happened to be on disk, so a partial or stale set
  # became the final sweep silently. The merge is a SEPARATE invocation, so PSW_CHUNK_N has to be
  # required here rather than merely consulted -- otherwise a merge launched without it skips the
  # check and the hole is still open. PSW_CHUNK_N is the worker count the run declared, not a
  # threshold invented here.
  if (.PCN <= 0L)
    stop("25 eGFR prior sweep MERGE: PSW_MERGE=1 was set without PSW_CHUNK_N. The merge cannot tell a complete chunk set from a partial one without knowing how many chunks the workers were split into; re-run the merge with the same PSW_CHUNK_N the workers used.")
  if (length(cf) != .PCN)
    stop(sprintf("25 eGFR prior sweep MERGE: found %d chunk file(s) in %s but the run declared PSW_CHUNK_N=%d. Refusing to merge a partial sweep.", length(cf), out_dir, .PCN))
  # THE EXACT FILES, NOT THAT-MANY FILES. A count check passes on a stale or misnumbered set, so
  # the names are required to be chunk00..chunk(N-1), each chunk must name the code that wrote it,
  # and the merged table must carry exactly one row per (gene, p12) that this run expects.
  .want <- file.path(out_dir, sprintf("egfr_sweep_chunk%02d.rds", seq_len(.PCN) - 1L))
  .miss <- setdiff(basename(.want), basename(cf))
  if (length(.miss))
    stop(sprintf("25 eGFR prior sweep MERGE: missing chunk file(s) %s in %s; found %s.",
                 paste(.miss, collapse = ", "), out_dir, paste(basename(cf), collapse = ", ")))
  .parts <- lapply(.want, readRDS)
  for (.i in seq_along(.want)) {
    .cv <- attr(.parts[[.i]], "psw_code_ver", exact = TRUE)
    if (is.null(.cv) || !identical(.cv, .psw_code_ver))
      stop(sprintf("25 eGFR prior sweep MERGE: %s carries code version %s but this run is %s -- it was computed by different code. Delete and regenerate that chunk.",
                   basename(.want[.i]), if (is.null(.cv)) "<unstamped>" else substr(.cv, 1, 12), substr(.psw_code_ver, 1, 12)))
  }
  egfr_sweep <- rbindlist(.parts, fill = TRUE)
  .exp <- CJ(gene = sort(unique(strong_egfr$gene)), p12 = p12_values)
  .got <- unique(as.data.table(egfr_sweep)[, .(gene, p12)])
  .gap <- fsetdiff(.exp, .got)
  if (nrow(.gap))
    stop(sprintf("25 eGFR prior sweep MERGE: %d expected (gene, p12) cell(s) absent from the merged sweep, e.g. %s -- the chunks do not cover this run's pairs.",
                 nrow(.gap), paste(utils::head(.gap[, paste(gene, p12)], 5), collapse = "; ")))
  .dup <- as.data.table(egfr_sweep)[, .N, by = .(gene, p12)][N > 1L]
  if (nrow(.dup))
    stop(sprintf("25 eGFR prior sweep MERGE: %d (gene, p12) cell(s) appear more than once, e.g. %s -- overlapping or duplicated chunks.",
                 nrow(.dup), paste(utils::head(.dup[, paste(gene, p12, N)], 5), collapse = "; ")))
  cat(sprintf("PSW-eGFR MERGE: %d chunk files -> %d rows\n", length(cf), nrow(egfr_sweep)))
}

# Pivot wider
egfr_wide <- dcast(
  egfr_sweep[, .(gene, outcome_id, outcome_name, p12, PP.H4, PP.H3)],   # unrounded: see the return above
  gene + outcome_id + outcome_name ~ p12,
  value.var = c("PP.H4", "PP.H3")
)
setnames(egfr_wide,
  c("PP.H4_1e-06", "PP.H4_5e-06", "PP.H4_1e-05"),
  c("PP.H4_p12_1e6", "PP.H4_p12_5e6", "PP.H4_p12_1e5"),
  skip_absent = TRUE
)

# Add original PP.H4
orig_sub <- strong_egfr[, .(gene, outcome_id, PP.H4_original = PP.H4)]
egfr_wide <- merge(egfr_wide, orig_sub, by = c("gene", "outcome_id"), all.x = TRUE)

# Stability classification
if ("PP.H4_p12_1e6" %in% names(egfr_wide) && "PP.H4_p12_1e5" %in% names(egfr_wide)) {
  # SINGLE SOURCE OF TRUTH (00_config.R; editor ruling 2026-08-01) -- see 10_phase3s_prior_sweep.R.
  egfr_wide[, stable := is_prior_stable(PP.H4_p12_1e6)]
  egfr_wide[, prior_band := prior_band(PP.H4_p12_1e6)]
  egfr_wide[, prior_shift_gt_020 := !stable | abs(PP.H4_p12_1e6 - PP.H4_p12_1e5) > 0.20]
}

cat(sprintf("\neGFR prior sweep: %d pairs\n", nrow(egfr_wide)))
cat(sprintf("Stable: %d | Sensitive: %d\n",
    sum(egfr_wide$stable, na.rm = TRUE), sum(egfr_wide$prior_shift_gt_020, na.rm = TRUE)))

if (exists("PSW_LIMIT")) {
  cat("\n[SMOKE] PSW_LIMIT set: eGFR coloc ran OK; skipping merge/save.\n")
  print(egfr_wide[, intersect(c("gene","PP.H4_original","PP.H4_p12_1e6","PP.H4_p12_1e5"), names(egfr_wide)), with = FALSE])
  q(save = "no", status = 0)
}

# ── Merge with existing non-eGFR results ─────────────────────────────────────
existing_wide <- readRDS(file.path(out_dir, "prior_sweep_wide.rds"))
setDT(existing_wide)

# Backup
saveRDS(existing_wide, file.path(out_dir, "prior_sweep_wide.rds.bak_pre_egfr"))

# 2026-07-01 Codex-p3rest #1: idempotent -- drop any pre-existing eGFR rows before appending fresh eGFR (a
# rerun where `existing` already holds eGFR would otherwise DUPLICATE them). Works for wide AND multi-row long.
.drop_egfr <- function(dt) { setDT(dt); if ("outcome_id" %in% names(dt)) dt[outcome_id != "Wuttke2019_eGFR_EA"] else if ("outcome_name" %in% names(dt)) dt[outcome_name != "eGFR"] else if ("outcome" %in% names(dt)) dt[outcome != "eGFR"] else dt }
existing_wide <- .drop_egfr(existing_wide)
combined_wide <- rbind(existing_wide, egfr_wide, fill = TRUE)
cat(sprintf("\nCombined: %d rows (%d non-eGFR + %d eGFR)\n",
    nrow(combined_wide), nrow(existing_wide), nrow(egfr_wide)))

# Save combined
saveRDS(combined_wide, file.path(out_dir, "prior_sweep_wide.rds"))
fwrite(combined_wide,  file.path(out_dir, "prior_sweep_wide.csv"))

# Also save eGFR-only long format
# NB: use egfr_sweep (long-format eGFR data, valid in all modes) not all_results,
# which is EMPTY in merge mode (loop skipped) -> would drop the eGFR long rows.
egfr_long <- egfr_sweep
existing_long <- readRDS(file.path(out_dir, "prior_sweep_long.rds"))
saveRDS(existing_long, file.path(out_dir, "prior_sweep_long.rds.bak_pre_egfr"))
existing_long <- .drop_egfr(existing_long)   # 2026-07-01 Codex-p3rest #1: idempotent (remove old eGFR before re-append)
combined_long <- rbind(existing_long, egfr_long, fill = TRUE)
saveRDS(combined_long, file.path(out_dir, "prior_sweep_long.rds"))

# ── Summary ──────────────────────────────────────────────────────────────────
cat(sprintf("\n=== Prior Sweep Summary (Combined %d pairs) ===\n", nrow(combined_wide)))
cat(sprintf("Total pairs: %d\n", nrow(combined_wide)))
cat(sprintf("Prior-robust (PP.H4>=0.80 at p12=1e-6; is_prior_stable): %d\n",
    sum(combined_wide$stable, na.rm = TRUE)))
cat(sprintf("Sensitive: %d\n",
    sum(combined_wide$prior_shift_gt_020, na.rm = TRUE)))

cat("\neGFR pairs detail:\n")
print_cols <- intersect(c("gene", "outcome_name", "PP.H4_original",
                          "PP.H4_p12_1e6", "PP.H4_p12_5e6", "PP.H4_p12_1e5",
                          "stable", "prior_band", "prior_shift_gt_020"),
                       names(egfr_wide))
print(egfr_wide[order(-PP.H4_original), ..print_cols])

cat("\n=== Done ===\n")
