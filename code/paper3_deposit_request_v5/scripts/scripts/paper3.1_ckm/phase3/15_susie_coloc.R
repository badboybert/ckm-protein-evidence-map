# Phase 3-C: SuSiE Fine-mapping Colocalization for the strong coloc hits (dual rule PP.H4>=0.80 & PP.H4/PP.H3>=5)
# =========================================================
# Method:
#   Re-extracts regional summary stats (same pipeline as coloc.R),
#   then adds LD matrix via ieugwasr::ld_matrix() (1KG EUR), runs
#   susie_rss() on exposure and outcome, then coloc.susie().
#
# NOTE (P3-F4 cold audit): runs fully OFFLINE (local PLINK + local coords); NO IEU OpenGWAS JWT required —
#   the JWT check was removed in the corrected cascade and the wrapper cascade/07_susie.R sets offline mode.
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase3/15_susie_coloc.R"

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

# susieR check
if (!requireNamespace("susieR", quietly = TRUE))
  stop("Package 'susieR' not installed. Run: install.packages('susieR')")
library(susieR)

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))   # for is_coloc_pass() (explicit direct source; coloc also loads it)
# ── Source coloc.R without running main() ──────────────────────────────────
COLOC_SOURCE_ONLY <- TRUE
suppressMessages(source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R")))
rm(COLOC_SOURCE_ONLY)

cat("\n=== Phase 3-C: SuSiE Fine-mapping Colocalization ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── JWT check REMOVED (localized 2026-06-04) ──────────────────────────────────
# SuSiE now uses LOCAL plink LD + LOCAL gene coords (via the fixed R/coloc.R);
# no OpenGWAS API is called, so no token is needed. (ieugwasr::user() was a live
# API preflight that aborts offline — removed per the local-first audit.)
cat("JWT check skipped: SuSiE runs fully local (plink LD + local gene coords).\n\n")

# ── Run control ───────────────────────────────────────────────────────────────
# Set FORCE_RERUN <- TRUE to discard checkpoint and re-run all pairs from scratch
# (needed when switching from API LD to local plink LD)
FORCE_RERUN <- exists("SUSIE_OUT_DIR") || nzchar(Sys.getenv("CKM_PALINDROME_SINK", ""))  # 2026-07-01 (gate finding 3): corrected cascade (sets SUSIE_OUT_DIR) forces a FRESH SuSiE run. Also force fresh when palindrome capture is on, so EVERY fit calls strand_align_R and the sink is complete (cold re-review High: a resumed checkpoint would skip capture)

# ── Output directory (override via SUSIE_OUT_DIR for the corrected cascade) ────
if (!exists("SUSIE_OUT_DIR")) SUSIE_OUT_DIR <- file.path(RESULTS_DIR, "phase3_susie")
out_dir <- SUSIE_OUT_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load strong hits (override inputs via SUSIE_COLOC_INPUT / SUSIE_INSTR) ──────
if (!exists("SUSIE_COLOC_INPUT")) SUSIE_COLOC_INPUT <- file.path(PHASE3_DIR, "coloc_results_final.rds")
if (!exists("SUSIE_INSTR"))       SUSIE_INSTR       <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")
coloc_final <- readRDS(SUSIE_COLOC_INPUT)
instr       <- readRDS(SUSIE_INSTR)
setDT(coloc_final); setDT(instr)

strong <- coloc_final[is_coloc_pass(PP.H4, PP.H3)][order(-PP.H4)]   # dual coloc rule (D1 2026-07-04)
# eGFR (rsid_only) has no local GWAS file → skip from the start
.egfr_id <- if (exists("OUTCOMES") && !is.null(OUTCOMES[["eGFR"]])) OUTCOMES[["eGFR"]] else { warning("[eGFR-id] OUTCOMES$eGFR missing; legacy fallback ebi-a-GCST003372"); "ebi-a-GCST003372" }
strong <- strong[outcome_id != .egfr_id]
# P3-F1 (cold audit): flag MHC-region genes -- the frozen method excludes MHC from SuSiE, but the raw artifact
# retained computed PP.H4 for MHC pairs (AGER~HF/TRIM40~T2D). Flag them so a DIRECT consumer of this artifact
# sees the method-mandated exclusion (downstream 12_rebuild already overrides MHC rows to not_evaluable/MHC_excluded).
strong[, mhc_excluded := if (exists("is_mhc_gene")) vapply(gene, function(x) isTRUE(is_mhc_gene(x)), logical(1)) else FALSE]
if (any(strong$mhc_excluded)) cat(sprintf("  [P3-F1] %d MHC-region strong pair(s) flagged mhc_excluded (SuSiE method excludes MHC): %s\n", sum(strong$mhc_excluded), paste(strong[mhc_excluded == TRUE, paste(gene, outcome_id)], collapse = ", ")))
cat(sprintf("Strong hits to analyse: %d pairs (eGFR excluded)\n\n", nrow(strong)))

rsid_map <- load_rsid_map()

# ── Pre-cache hg38 gene coordinates (avoid rate-limiting during main loop) ─────
# Outcomes needing hg38 coords: FinnGen (hg38) and GWAS Catalog (hg38_hm)
# hg37 outcomes (Nielsen AF) don't call get_gene_hg38_coords
needs_hg38 <- unique(strong$gene[
  strong$outcome_id %in% names(GWAS_META)[
    sapply(GWAS_META, function(m) !is.null(m$genome) && m$genome %in% c("hg38","hg38_hm"))
  ]
])
cat(sprintf("Pre-fetching hg38 coords for %d genes...\n", length(needs_hg38)))
gene_hg38_cache <- list()
for (g in needs_hg38) {
  for (attempt in seq_len(3)) {
    coords <- tryCatch(get_gene_hg38_coords(g), error = function(e) NULL)
    if (!is.null(coords)) {
      gene_hg38_cache[[g]] <- coords
      cat(sprintf("  OK  %s : chr%s:%d-%d\n", g, coords$chrom, coords$start, coords$end))
      break
    }
    cat(sprintf("  [retry %d] %s\n", attempt, g))
    Sys.sleep(3)
  }
  if (is.null(gene_hg38_cache[[g]])) cat(sprintf("  [FAIL] %s\n", g))
  Sys.sleep(1)  # 1 s between genes
}
cat(sprintf("Pre-cache done: %d / %d genes resolved\n\n",
            sum(!sapply(gene_hg38_cache, is.null)), length(needs_hg38)))

# 2026-07-01 Codex-p3rest #3: every needed gene MUST resolve from the LOCAL coord table (gene_hg38_coords.tsv,
# built by 03_gene_coords). A missing one is a provenance break in the local table -> fail loud; do NOT paper
# it over with a hardcoded coordinate (the old PCSK9 hardcoded fallback + its "direct API call" note removed).
.unresolved <- needs_hg38[vapply(needs_hg38, function(g) is.null(gene_hg38_cache[[g]]), logical(1))]
if (length(.unresolved)) stop(sprintf("15 SuSiE: %d strong-pair gene(s) missing from the local hg38 coord table (%s) -- rebuild gene_hg38_coords.tsv via cascade/03b_gene_coords_hg38.R (NOT 03_gene_coords.R, which is the retired liftover version); refusing a hardcoded/API coordinate fallback", length(.unresolved), paste(head(.unresolved, 10), collapse = ",")))

# Override get_gene_hg38_coords with cache lookup
base_get_gene_hg38_coords <- get_gene_hg38_coords  # keep original
get_gene_hg38_coords <- function(gene_name) {
  if (gene_name %in% names(gene_hg38_cache)) return(gene_hg38_cache[[gene_name]])
  # Not in cache: try original (with extra delay), then cache result
  Sys.sleep(1)
  coords <- tryCatch(base_get_gene_hg38_coords(gene_name), error = function(e) NULL)
  if (!is.null(coords)) gene_hg38_cache[[gene_name]] <<- coords
  coords
}

# ── Helper: LD matrix — local plink preferred, API fallback ────────────────────
# When local 1KG EUR files exist (PLINK_BIN + LD_EUR_DIR in config), uses all
# available rsID SNPs without the 400-SNP API limit.
# Root-cause LD sign-alignment helper (2026-07-17, editor ruling 2b). Sourced INSIDE this block (not at file
# top) so the A" consumers (29/31/35/37), which eval only the sliced helper-through-fitting-loop region, also
# define strand_align_R, not just the full-source cascade run. (This comment must NOT contain the literal
# slice-marker phrases the A" greps search for, or it would shift their slice boundary -- learned the hard way.)
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/strand_align.R"))

# THE ITERATION CEILING LIVES HERE, NOT AT THE FILE HEAD, for the same reason strand_align.R does:
# 29/31/35/37 evaluate only this sliced region, so anything they need must be defined inside it. A
# definition at the top would be invisible to them and every candidate fit would fail on an unknown
# object. DEFAULT 1000 is deliberate -- those consumers must keep the ceiling their frozen artefacts
# were produced under, without either script being edited. Only the strong-set trunk (cascade/07,
# cascade/08) sets SUSIE_MAX_ITER, and it does so before sourcing, so this line does not override it.
# 2026-08-08: raising it globally was tried and reverted -- at 10000 the candidate fit OMG~T2D
# terminated with susieR's own error for an unreasonably large estimated prior variance.
if (!exists("SUSIE_MAX_ITER")) SUSIE_MAX_ITER <- 1000L

get_ld_matrix_auto <- function(rsids, gene_chrom, max_snps_api = 400) {
  rsids <- unique(rsids[grepl("^rs[0-9]+$", rsids)])
  if (length(rsids) < 5) {
    cat("    [LD] < 5 rsIDs available\n"); return(NULL)
  }

  # ── Local plink path ────────────────────────────────────────────────────────
  use_local <- exists("PLINK_BIN") && !is.null(PLINK_BIN) && file.exists(PLINK_BIN) &&
               exists("LD_EUR_DIR") && !is.null(LD_EUR_DIR)
  if (use_local) {
    # 1KG EUR tgz 是整合的單一 EUR.bed/bim/fam（非 per-chromosome）
    bfile    <- file.path(LD_EUR_DIR, "EUR")
    bim_file <- paste0(bfile, ".bim")
    if (file.exists(bim_file)) {
      cat(sprintf("    [LD local] EUR plink: %d SNPs, chr%s (no limit)\n",
                  length(rsids), gene_chrom))
      ld <- tryCatch(
        ieugwasr::ld_matrix(rsids, bfile = bfile, plink_bin = PLINK_BIN),
        error = function(e) fail_loud_data(sprintf("15 SuSiE local LD (PLINK) FAILED: %s -- LD infrastructure failure, not a biological skip", conditionMessage(e)))  # 2026-07-01 (closure B1, Codex): mirror of 26:226
      )
      if (!is.null(ld)) {
        cat(sprintf("    [LD local] Returned %d x %d matrix\n", nrow(ld), ncol(ld)))
        return(ld)
      }
      cat("    [LD local] Failed (local-only; NO API fallback by design)\n"); return(NULL)
    } else {
      fail_loud_data(sprintf("15 SuSiE: LD reference EUR.bim NOT FOUND (%s) -- infrastructure failure; cannot compute LD for ANY pair, not a per-pair skip", bim_file))   # 2026-07-01 Codex-p3rest #4
    }
  }

  # ── No API fallback (localized 2026-06-04): require local PLINK + 1KG EUR ────
  stop("SuSiE LD requires local PLINK_BIN + LD_EUR_DIR/EUR.bim (no OpenGWAS API fallback).")
}

# ── Helper: run SuSiE coloc for one (gene, outcome_id) pair ────────────────────
run_susie_pair <- function(gene, outcome_id, gene_chrom, gene_start, gene_end,
                            rsid_map, pp4_abf) {
  meta   <- GWAS_META[[outcome_id]]
  genome <- meta$genome

  region_start <- max(1L, gene_start - CIS_WINDOW_KB * 1000L)
  region_end   <- gene_end + CIS_WINDOW_KB * 1000L

  cat(sprintf("  susie: %s ~ %s | chr%s:%d-%d\n",
              gene, outcome_id, gene_chrom, region_start, region_end))
  flush(stdout())

  # ── 1. pQTL regional data (hg37) ────────────────────────────────────────────
  pqtl_dt <- get_pqtl_region(gene, as.character(gene_chrom),
                               region_start, region_end, rsid_map)
  if (is.null(pqtl_dt) || nrow(pqtl_dt) == 0) {
    # A DATA-FOUNDATION FAILURE, NOT A LOW-SNP SKIP. R/coloc.R:963-969 stops on an empty extraction
    # and on missing required columns; these copies returned NULL, which the caller renders as an NA
    # row or an omitted pair -- so a broken archive or a schema change looked like an uncomputable
    # biological result. fail_loud_data() already exists here and the main loop turns it into an abort.
    fail_loud_data(sprintf("15: gene %s pQTL extraction returned NO data (tar resolved but empty) -- data-foundation failure; mirrors R/coloc.R:963", gene))
  }
  req_cols <- c("ALLELE0", "ALLELE1", "BETA", "SE", "INFO", "N")
  if (any(!req_cols %in% names(pqtl_dt))) {
    fail_loud_data(sprintf("15: gene %s pQTL missing required columns: %s -- data-format failure; mirrors R/coloc.R:969", gene, paste(setdiff(req_cols, names(pqtl_dt)), collapse = ", ")))
  }
  pqtl_dt <- pqtl_dt[!is.na(INFO) & INFO >= INFO_MIN & !is.na(BETA) & !is.na(SE)]
  if (nrow(pqtl_dt) < 30) {
    cat(sprintf("    [skip] pQTL only %d SNPs after filter\n", nrow(pqtl_dt)))
    return(NULL)
  }
  pqtl_dt[, pqtl_varbeta := SE^2]
  pqtl_dt[, ea_pqtl      := toupper(ALLELE1)]
  pqtl_dt[, oa_pqtl      := toupper(ALLELE0)]
  n_pqtl <- as.integer(pqtl_dt[["N"]][1])

  # ── 2. Outcome data + pQTL_sub construction (mirrors run_coloc_pair) ─────────
  if (genome == "hg38") {
    # NO LIFTOVER: pqtl_dt$chrpos is chrom:GENPOS from the UKB-PPP archive and GENPOS is
    # already GRCh38 (verified against the archive itself, not against a comment).
    # liftover_to_hg38() has been DELETED -- it lifted GRCh38 a second time, displacing
    # positions by up to 864,716 bp, and every position-fallback match from it paired
    # two different variants.
    out_dt    <- get_finngen_region(outcome_id, gene)
    if (is.null(out_dt) || nrow(out_dt) < 20) {
      cat("    [skip] FinnGen no data\n"); return(NULL)
    }
    chrpos38_col <- pqtl_dt$chrpos
    pqtl_sub <- data.table(
      rsid         = pqtl_dt$rsid,
      chrpos_pqtl  = chrpos38_col,
      pqtl_beta    = pqtl_dt$BETA,
      pqtl_varbeta = pqtl_dt$pqtl_varbeta,
      ea_pqtl      = pqtl_dt$ea_pqtl,
      oa_pqtl      = pqtl_dt$oa_pqtl
    )
    out_sub <- out_dt[, .(
      rsid_out, chrpos_out,
      ea_out      = toupper(get(meta$ea_col)),
      oa_out      = toupper(get(meta$oa_col)),
      out_beta    = as.numeric(get(meta$beta_col)),
      out_varbeta = as.numeric(get(meta$se_col))^2
    )]
    out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

  } else if (genome == "hg37") {
    # PROJECTED TO GRCh37. The outcome here reports Pos_b37 while pqtl_dt$chrpos is GRCh38, so the
    # two cannot be compared as they stand -- and a direct comparison does not merely fail to match,
    # it matches the WRONG variants wherever the two builds happen to agree numerically.
    pqtl_37 <- liftover_to_hg37(pqtl_dt)
    .psrc <- if (!is.null(pqtl_37)) pqtl_37 else pqtl_dt

    # THE WINDOW IS PROJECTED TOO, and the extraction now happens AFTER the lift rather than before.
    # region_start/region_end are GRCh38, built from the GRCh38 gene body, while get_nielsen_region()
    # filters Pos_b37: this branch lifted the variants and handed the bounds over unconverted, which
    # is the build-mixing the rebuild exists to remove. Measured on DNER (chr2): the GRCh38 bounds
    # 228,857,519-230,214,586 versus the lifted 229,722,334-231,078,645 -- an offset of +864,815 bp,
    # so the two intervals barely overlap. The hg37 outcome on this path is Wuttke CKD, so every CKD
    # pair here was affected. The endpoints are not lifted (they can land in different chain blocks
    # or out of order); the bounds come from the projected positions of the pQTLs being colocalised.
    # An off-chromosome lift drops that one variant, per the engine (R/coloc.R:1049); only having
    # nothing left on this chromosome is fatal.
    .chr37  <- if (!is.null(pqtl_37)) sub(":.*", "", pqtl_37$chrpos37) else NA_character_
    .same37 <- .chr37 %in% as.character(gene_chrom)
    .off37  <- !is.na(.chr37) & !.same37
    if (!any(.same37))
      fail_loud_data(sprintf(paste0("15 SuSiE %s ~ %s: the pQTL liftover to GRCh37 leaves NO variant ",
                                    "on chr%s (mapped=%d, off-chromosome=%d, unmapped=%d), so there ",
                                    "is no window. Falling back to the GRCh38 bounds would query the ",
                                    "wrong region silently."),
                             gene, outcome_id, gene_chrom, sum(.same37), sum(.off37),
                             sum(is.na(.chr37))))
    if (any(.off37)) {
      cat(sprintf("    [liftover] %d pQTL(s) lifted off chr%s -- dropped for this outcome\n",
                  sum(.off37), gene_chrom))
      pqtl_37$chrpos37[.off37] <- NA_character_
    }
    .pos37 <- as.integer(sub(".*:", "", pqtl_37$chrpos37[.same37]))
    out_dt <- get_nielsen_region(outcome_id, as.character(gene_chrom), min(.pos37), max(.pos37))
    if (is.null(out_dt) || nrow(out_dt) < 20) {
      cat("    [skip] Nielsen no data\n"); return(NULL)
    }
    pqtl_sub <- data.table(
      rsid         = .psrc$rsid,
      chrpos_pqtl  = if (!is.null(pqtl_37)) pqtl_37$chrpos37 else rep(NA_character_, nrow(pqtl_dt)),
      pqtl_beta    = .psrc$BETA,
      pqtl_varbeta = .psrc$pqtl_varbeta,
      ea_pqtl      = .psrc$ea_pqtl,
      oa_pqtl      = .psrc$oa_pqtl
    )
    out_sub <- out_dt[, .(
      rsid_out, chrpos_out,
      ea_out      = toupper(get(meta$ea_col)),
      oa_out      = toupper(get(meta$oa_col)),
      out_beta    = as.numeric(get(meta$beta_col)),
      out_varbeta = as.numeric(get(meta$se_col))^2
    )]
    out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

  } else if (genome == "hg38_hm") {
    # NO LIFTOVER: pqtl_dt$chrpos is chrom:GENPOS from the UKB-PPP archive and GENPOS is
    # already GRCh38 (verified against the archive itself, not against a comment).
    # liftover_to_hg38() has been DELETED -- it lifted GRCh38 a second time, displacing
    # positions by up to 864,716 bp, and every position-fallback match from it paired
    # two different variants.
    out_dt    <- get_harmonized_hg38_region(outcome_id, gene)
    if (is.null(out_dt) || nrow(out_dt) < 20) {
      cat("    [skip] harmonized no data\n"); return(NULL)
    }
    pqtl_src <- pqtl_dt
    pqtl_sub <- data.table(
      rsid         = pqtl_src$rsid,
      chrpos_pqtl  = pqtl_dt$chrpos,
      pqtl_beta    = pqtl_src$BETA,
      pqtl_varbeta = pqtl_src$pqtl_varbeta,
      ea_pqtl      = pqtl_src$ea_pqtl,
      oa_pqtl      = pqtl_src$oa_pqtl
    )
    out_sub <- out_dt[, .(
      rsid_out, chrpos_out,
      ea_out      = toupper(get(meta$ea_col)),
      oa_out      = toupper(get(meta$oa_col)),
      out_beta    = as.numeric(get(meta$beta_col)),
      out_varbeta = as.numeric(get(meta$se_col))^2
    )]
    out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

  } else {
    cat(sprintf("    [skip] unknown genome: %s\n", genome)); return(NULL)
  }

  # ── 3. Merge + harmonise (same as run_coloc_pair) ───────────────────────────
  common <- rsid_pos_merge(pqtl_sub, out_sub)
  # THE OUTCOME NAMED THIS VARIANT, SO ITS rsID IS USABLE FOR LD. The engine's position fallback
  # (R/coloc.R:282) drops rsid_out and sets snp_id to the coordinate, and the LD step below keeps
  # only `^rs[0-9]+$` -- so a position-matched variant the outcome named was excluded from LD even
  # though the EUR panel is rsID-keyed and would have found it. 26:271 already does this (editor-
  # endorsed, the ACP1 case); this copy did not. Done here and not in rsid_pos_merge() because that
  # function is shared with 04/05/06, where snp_id is also the de-duplication key.
  #
  # KEYED ON (COORDINATE, EA, OA), NOT ON THE COORDINATE. "Exactly one rsID somewhere at this
  # position" does not establish that the outcome row which actually merged owns it: at a
  # multiallelic coordinate a compatible row with no rsID would inherit an incompatible row's sole
  # rsID. `common` carries the merged outcome alleles, so the triple names one outcome record --
  # the same identity the allele-aware rsID map and the engine's matched_pqtl_var use.
  if (nrow(common) && all(c("chrpos_pqtl", "ea_out", "oa_out") %in% names(common)) &&
      "rsid_out" %in% names(out_sub)) {
    .ok <- out_sub[!is.na(chrpos_out) & !is.na(rsid_out) & grepl("^rs[0-9]+$", rsid_out)]
    .ok[, .key := paste(chrpos_out, toupper(ea_out), toupper(oa_out), sep = ":")]
    .uid <- .ok[, .(.n_id = uniqueN(rsid_out), .id = rsid_out[1L]), by = .key]
    .ck  <- paste(common$chrpos_pqtl, toupper(common$ea_out), toupper(common$oa_out), sep = ":")
    .hit <- .uid[match(.ck, .uid$.key)]
    .rec <- !grepl("^rs[0-9]+$", common$snp_id) & !is.na(.hit$.id) & .hit$.n_id == 1L
    if (any(.rec)) {
      common[.rec, snp_id := .hit$.id[.rec]]
      cat(sprintf("    rsID recovered for %d position-matched SNP(s)\n", sum(.rec)))
    }
  }

  cat(sprintf("    merged: %d SNPs\n", nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] < 20 after merge\n"); return(NULL) }

  # HARMONISE BEFORE DE-DUPLICATING, then re-key only the rsIDs that still collide -- the engine
  # order (R/coloc.R:1148-1163). Dropping ambiguous groups first deleted variants that
  # harmonisation would have resolved to distinct outcome rows: 1:1251285 G/A and G/T are both
  # rs6657544, and both were being thrown away.
  common <- harmonize_alleles(common, "ea_pqtl", "oa_pqtl", "ea_out", "oa_out", "out_beta")
  cat(sprintf("    after allele align: %d SNPs\n", nrow(common)))
  common[duplicated(snp_id) | duplicated(snp_id, fromLast = TRUE),
         snp_id := paste(chrpos_pqtl, ea_pqtl, oa_pqtl, sep = ":")]
  n_before <- nrow(common)
  common <- common[, if (.N == 1L) .SD, by = snp_id]
  if (nrow(common) < n_before)
    cat(sprintf("    dedup: %d -> %d\n", n_before, nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] < 20 after harmonise\n"); return(NULL) }

  # ── 4. LD matrix (1KG EUR via ieugwasr) ─────────────────────────────────────
  # Sort by pQTL |z| descending so the pQTL lead SNP is always in the top 400
  common[, z_pqtl_abs := abs(pqtl_beta / sqrt(pqtl_varbeta))]
  common <- common[order(-z_pqtl_abs)]

  snps_for_ld <- common$snp_id[grepl("^rs[0-9]+$", common$snp_id)]
  cat(sprintf("    SNPs with rsID for LD: %d / %d (sorted by pQTL |z|)\n",
              length(snps_for_ld), nrow(common)))

  if (length(snps_for_ld) < 10) {
    cat("    [skip] < 10 rsIDs for LD computation\n"); return(NULL)
  }

  ld_mat <- get_ld_matrix_auto(snps_for_ld, gene_chrom)
  if (is.null(ld_mat)) {
    cat("    [skip] LD matrix failed\n"); return(NULL)
  }

  # ── 5. Align common to LD matrix ─────────────────────────────────────────────
  # ieugwasr::ld_matrix() returns rownames as "rsID_A_B" (allele-suffixed).
  # NOTE (verified 2026-06-05): plink runs with --keep-allele-order and our local
  # 1KG EUR bim is effect-allele-aligned to UKB-PPP, so ld_matrix A1 == pQTL effect
  # allele for ALL SNPs (checked: 0 sign-flips across PCSK9/MINDY1/CCND2 regions).
  # Hence matching by plain rsID is sign-correct here; no LD sign harmonisation is
  # needed. (The IBSS non-convergence comes from external-LD-vs-GWAS-sample mismatch,
  # not allele misalignment.) Re-verify the flip count if the LD panel ever changes.
  ld_snps_raw <- rownames(ld_mat)
  ld_snps_rs  <- sub("_[A-Za-z]+_[A-Za-z]+$", "", ld_snps_raw)  # rs12345_T_G → rs12345

  # Match common$snp_id (plain rsID) to ld_snps_rs
  idx_in_ld   <- match(ld_snps_rs, common$snp_id)
  keep_ld_pos <- which(!is.na(idx_in_ld))

  if (length(keep_ld_pos) < 10) {
    cat(sprintf("    [skip] Only %d SNPs matched between common and LD matrix\n",
                length(keep_ld_pos)))
    return(NULL)
  }

  common_ld  <- common[idx_in_ld[keep_ld_pos]]
  ld_mat_f   <- ld_mat[keep_ld_pos, keep_ld_pos, drop = FALSE]
  # Keep plain rsID rownames for clarity
  rownames(ld_mat_f) <- colnames(ld_mat_f) <- ld_snps_rs[keep_ld_pos]

  cat(sprintf("    LD-aligned: %d SNPs (from %d in LD matrix)\n",
              nrow(common_ld), length(ld_snps_raw)))
  if (nrow(common_ld) < 10) {
    cat("    [skip] < 10 SNPs after LD alignment\n"); return(NULL)
  }

  # ── 5b. STRAND-ALIGN R to the pQTL effect allele (root-cause fix 2026-07-17, editor 2b) ────────────────
  # Re-orient R for SNPs where the .bim A1 (recovered from the allele-suffixed LD rownames) is the pQTL OTHER
  # allele, so z (ea_pqtl-oriented) and R agree. No-op for the already-aligned majority. Palindromes are left
  # as-is (verified separately by allele frequency); allele_mismatch SNPs are dropped.
  .sa <- strand_align_R(common_ld, ld_snps_raw[keep_ld_pos], ld_mat_f, label = paste0(gene, "~", outcome_id))
  common_ld <- .sa$common_ld; ld_mat_f <- .sa$R
  if (nrow(common_ld) < 10) { cat("    [skip] < 10 SNPs after strand-align drop\n"); return(NULL) }   # re-check gate after drops (R2 Medium)

  # ── 6. susie_rss ─────────────────────────────────────────────────────────────
  z_exp <- common_ld$pqtl_beta / sqrt(common_ld$pqtl_varbeta)
  z_out <- common_ld$out_beta  / sqrt(common_ld$out_varbeta)
  # THE N OF THE VARIANTS BEING FITTED, when the outcome reports one per variant. meta$n is a single
  # declared scalar; the Wuttke CKD file carries n_total_sum per variant, median 402,682 over the
  # file and about 403,000 at the two strong CKD loci, while the declared 480,698 is ABOVE that
  # file''s maximum of 444,971 -- so it is not a coarser version of the truth, it is a value no
  # variant has. susie_rss uses n to build its sufficient statistics, unlike coloc.abf where n is
  # inert once beta, varbeta and sdY are given. Outcomes with no per-variant column keep meta$n.
  n_out <- {
    .nc <- if ("n_total_sum" %in% names(out_dt)) out_dt[match(common_ld$chrpos_pqtl, chrpos_out), n_total_sum] else NULL
    if (is.null(.nc) || all(is.na(.nc))) as.integer(meta$n) else as.integer(median(.nc, na.rm = TRUE))
  }

    # max_iter = SUSIE_MAX_ITER (default 1000; the strong-set trunk passes 10000).
    # CEILING HISTORY: susieR's default is 100 -- a package default, never a choice
    # this study made and never stated anywhere in it. Raised to 1000 on 2026-07-16, to 10000 on
    # 2026-08-08. Both raises are recorded below; the number on the calls is the current ceiling.
    #
    # WHY THE CEILING IS THIS STUDY'S TO SET. coloc::runsusie, the wrapper the package provides for
    # exactly this task, defaults to repeat_until_convergence = TRUE, because its authors treat
    # non-convergence at the default limit as a computational artifact rather than a result. This
    # engine bypasses that wrapper and calls susie_rss directly, so the stopping rule is ours and has
    # to be stated rather than inherited.
    #
    # IT IS FORWARDED, NOT SWALLOWED -- and that had to be checked, because `check_LD` on this same
    # line is NOT: max_iter is absent from formals(susie_rss) exactly as check_LD is, but
    # susie_suff_stat (which susie_rss forwards `...` to) HAS max_iter as a formal and check_LD is a
    # formal of neither. Verified empirically 2026-07-16: max_iter=1 -> converged FALSE,
    # max_iter=1000 -> converged TRUE.
    #
    # RAISING IT IS FREE FOR FITS THAT ALREADY CONVERGE. IBSS stops at its own niter, not at the
    # ceiling -- measured 2026-08-08, RTBDN-T2D exposure converged at niter = 13 under a ceiling of
    # 100 -- and 47 of 52 pairs were byte-identical across the 100 -> 1000 raise. Only fits that do
    # not converge run longer.
    #
    # THE 2026-08-08 RAISE follows a rule registered before the result was known: any fit reporting
    # converged == FALSE is refit once at the higher ceiling. Two strong pairs qualify, NCAN-CAD
    # (both sides) and ADM-AF (exposure), and their ELBO traces at the 1000 ceiling are not alike.
    # ADM's exposure step size decayed 4.3 orders of magnitude (230.3 -> 0.0103): it sits on the
    # optimum and merely missed the tolerance. NCAN's did not decay at all -- its whole-run mean step
    # equals its last-50 median to within 0.3% (208.21 vs 208.8), so that ELBO is linear from the
    # first iteration, which is what an unattainable optimum looks like rather than a slow one.
  S1 <- tryCatch(
    susieR::susie_rss(z = z_exp, R = ld_mat_f, n = n_pqtl,
                       L = 10, min_abs_corr = 0, verbose = FALSE, check_LD = FALSE, max_iter = SUSIE_MAX_ITER),
    error = function(e) fail_loud_data(sprintf("15 SuSiE susie_rss(exp) THREW an error (malformed data/LD, not a legit no-CS result): %s", conditionMessage(e)))  # 2026-07-01 (closure B1, Codex): mirror of 26:328
  )
  if (is.null(S1)) { cat("    [skip] susie_rss failed for exposure\n"); return(NULL) }

  S2 <- tryCatch(
    susieR::susie_rss(z = z_out, R = ld_mat_f, n = n_out,
                       L = 10, min_abs_corr = 0, verbose = FALSE, check_LD = FALSE, max_iter = SUSIE_MAX_ITER),
    error = function(e) fail_loud_data(sprintf("15 SuSiE susie_rss(out) THREW an error (malformed data/LD, not a legit no-CS result): %s", conditionMessage(e)))  # 2026-07-01 (closure B1, Codex): mirror of 26:335
  )
  if (is.null(S2)) { cat("    [skip] susie_rss failed for outcome\n"); return(NULL) }

  n_cs_exp <- length(S1$sets$cs)
  n_cs_out <- length(S2$sets$cs)
  cat(sprintf("    SuSiE credible sets: exp=%d, out=%d\n", n_cs_exp, n_cs_out))

  # ── 7. coloc.susie ───────────────────────────────────────────────────────────
  cs_res <- tryCatch(
    coloc::coloc.susie(S1, S2),
    error = function(e) {
      message("    [coloc.susie error] ", conditionMessage(e)); NULL
    }
  )
  # THE FITS RAN, SO THEIR CONVERGENCE IS KNOWN AND MUST SURVIVE THIS EXIT. coloc.susie failing (or
  # returning an empty summary) is a failure of the COLOCALIZATION, not of the SuSiE fits -- S1 and S2 are
  # complete by this point and their `$converged` is the answer this change exists to record. Returning a
  # bare NULL here sends the pair to the caller's no-fit fallback, whose row has no convergence columns, so
  # a fit whose convergence we KNOW would be written as if we did not. That path is not hypothetical: it is
  # exactly the `dropped:collapsed` class (4 of the 17 dashed candidates). rbindlist(fill=TRUE) then hides
  # it by filling NA. Return the row with the posteriors NA and the convergence carried.
  if (is.null(cs_res)) {
    cat("    [skip] coloc.susie failed\n")
    return(data.table(
      gene = gene, outcome_id = outcome_id,
      n_snps_merged = nrow(common), n_snps_ld = nrow(common_ld),
      n_signals_exp = n_cs_exp, n_signals_out = n_cs_out, n_signal_pairs = NA_integer_,
      PP.H4_abf = as.numeric(pp4_abf), PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      concordant_susie = NA,
      converged_exp = isTRUE(S1$converged), converged_out = isTRUE(S2$converged),
      max_iter_used = SUSIE_MAX_ITER))
  }

  summary_dt <- as.data.table(cs_res$summary)
  if (nrow(summary_dt) == 0 || !"PP.H4.abf" %in% names(summary_dt)) {
    cat("    [skip] coloc.susie returned empty summary\n")
    return(data.table(
      gene = gene, outcome_id = outcome_id,
      n_snps_merged = nrow(common), n_snps_ld = nrow(common_ld),
      n_signals_exp = n_cs_exp, n_signals_out = n_cs_out, n_signal_pairs = NA_integer_,
      PP.H4_abf = as.numeric(pp4_abf), PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      concordant_susie = NA,
      converged_exp = isTRUE(S1$converged), converged_out = isTRUE(S2$converged),
      max_iter_used = SUSIE_MAX_ITER))
  }

  # Best signal pair = row with max PP.H4
  best <- summary_dt[which.max(PP.H4.abf)]
  pp4_susie <- best$PP.H4.abf
  pp3_susie <- best$PP.H3.abf
  n_pairs   <- nrow(summary_dt)

  cat(sprintf("    coloc.susie: %d signal pairs | best PP.H4=%.3f PP.H3=%.3f\n",
              n_pairs, pp4_susie, pp3_susie))

  data.table(
    gene             = gene,
    outcome_id       = outcome_id,
    n_snps_merged    = nrow(common),
    n_snps_ld        = nrow(common_ld),
    n_signals_exp    = n_cs_exp,
    n_signals_out    = n_cs_out,
    n_signal_pairs   = n_pairs,
    # NOT ROUNDED. concordant_susie below is computed from the full-precision value, so storing a
    # rounded one let a row contradict its own flag; 26 writes unrounded eGFR rows into the same
    # combined artifact, so the precision was also outcome-dependent. Same change as 25.
    PP.H4_abf        = as.numeric(pp4_abf),
    PP.H4_susie      = as.numeric(pp4_susie),
    PP.H3_susie      = as.numeric(pp3_susie),
    # CONVERGENCE GATES THE CLASSIFICATION, not merely the record. Both fits' `converged` were being
    # stored below and read by nothing here, so a pair whose fit never converged was classified from
    # its posterior as an ordinary result: NCAN-CAD (both sides FALSE) came out DISCORDANT on
    # PP.H4_susie = 1.777416e-48, and ADM-AF (exposure FALSE) came out inconclusive. The downstream
    # pair-evidence code does check convergence and reclassifies both as not_evaluable, so the
    # authoritative 43/28/1/4 split is unaffected -- but this sidecar is read DIRECTLY by the grid and
    # by supplementary tables, where the classification was wrong.
    # NA, not FALSE: "the fit did not converge" is not "the pair does not colocalise", and NA is
    # already this column's value on the coloc.susie-failed path above, so no consumer sees a new
    # value class.
    concordant_susie = if (isTRUE(S1$converged) && isTRUE(S2$converged)) pp4_susie >= 0.80 else NA,
    # Convergence, RECORDED PER FIT AND PER SIDE (2026-07-16). Until now nothing in this study recorded it
    # except 31, and only per EDGE across both calls: the canonical logs preserve unattributed counts (6
    # over the non-eGFR pairs, 3 over the eGFR pairs) and the 268 candidate fits had no record at all. So
    # "did the primary fits converge?" was unanswerable, and the attribution of the three held edges'
    # non-convergence to the eGFR outcome rather than the UKB-PPP exposure was never established -- one
    # `n_warn` per edge cannot tell those apart. susie_rss returns `$converged` directly, so this needs no
    # warning capture: read the fit's own answer, per side.
    max_iter_used    = SUSIE_MAX_ITER,   # the ceiling this fit ACTUALLY ran under; without it a
                                         # non-converged row cannot be told apart from one that
                                         # silently ran at the default (REBUILD_LEDGER §9bb)
    converged_exp    = isTRUE(S1$converged),
    converged_out    = isTRUE(S2$converged)
  )
}

# ── Main loop with checkpoint ──────────────────────────────────────────────────
checkpoint_file <- file.path(out_dir, "susie_checkpoint.rds")
results    <- list()
done_pairs <- character(0)

if (FORCE_RERUN && file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
  cat("[FORCE_RERUN] Checkpoint cleared — fresh start\n\n")
}

# Strand-align version guard (cold re-review Medium): a DIRECT run leaves FORCE_RERUN=FALSE, so a checkpoint
# written before the LD sign-alignment fix (unaligned R) could be resumed with its already-done pairs skipped.
# Stamp checkpoints with the strand_align.R md5; a missing/mismatched stamp means "different engine" -> discard.
# THE TOKEN MUST COVER THE CODE THAT PRODUCED THE CACHED VALUES, not one helper. A cache carries the
# LD matrix, the matched SNP set and n_out, all of which depend on this script and on R/coloc.R
# (liftover, merge, harmonisation), not only on the strand helper. Any of the three unreadable ->
# token NA -> rebuild everything (fail safe).
.susie_code_ver <- tryCatch(paste(vapply(c("R/strand_align.R", "R/coloc.R", "scripts/paper3.1_ckm/phase3/15_susie_coloc.R"),
                            function(.f) unname(tools::md5sum(file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"), .f))),
                            character(1)), collapse = ":"),
                           error = function(e) NA_character_)
if (file.exists(checkpoint_file)) {
  .cpv <- tryCatch(readRDS(checkpoint_file)$code_ver, error = function(e) NULL)
  if (is.null(.cpv) || !identical(.cpv, .susie_code_ver)) {
    file.remove(checkpoint_file)
    cat("[strand-align] checkpoint predates the LD sign-alignment fix (no/mismatched code_ver) — discarded, fresh start\n\n")
  }
}

# A CHECKPOINT FROM A DIFFERENT ENGINE IS NOT A RESUME, IT IS A SILENT HALF-RUN.
# The rows it carries were computed by whatever this file said at the time; resuming skips them, so they
# reach the artifact unchanged while every other row comes from the current engine. `rbindlist(fill = TRUE)`
# at the combine then papers over any schema difference with NA, and nothing errors.
# Not hypothetical (2026-07-16): a checkpoint dated 07-09 held **38 completed pairs — every non-eGFR pair
# in the study** — with no `converged_*` columns. Resuming it would have recorded "convergence unknown" for
# 38 fits that ran and whose answer was available, i.e. the per-fit convergence recording would have been
# inert for the entire 52 while appearing to work.
# 29 does not have this hole because its checkpoint is keyed on an md5 of the engine source (`code_md5`), so
# any engine edit invalidates it. That is the more general guard and 15 should eventually have it (§3.1).
# This is the narrow version: the schema the current engine promises must be present, or refuse to resume.
if (file.exists(checkpoint_file)) {
  .cp <- readRDS(checkpoint_file)
  .cr <- .cp$results
  if (length(.cr)) {
    # EVERY row, and the VALUES -- not names(.cr[[1]]). Checking the first element's column names is a
    # sample, not a check: a checkpoint half-written across an engine change has both shapes in it, and the
    # first one may be the good one. And presence is not correctness -- a checkpoint already tainted by an
    # earlier fill has the columns and NA values. assert_convergence_honest() does both.
    assert_convergence_honest(data.table::rbindlist(.cr, fill = TRUE, use.names = TRUE),
                              sprintf("15: %s (checkpoint, %d completed pair(s))", basename(checkpoint_file), length(.cr)))
  }
}

if (file.exists(checkpoint_file)) {
  cp <- readRDS(checkpoint_file)
  results    <- cp$results
  done_pairs <- cp$done_pairs
  cat(sprintf("Resuming from checkpoint: %d done, %d with results\n\n",
              length(done_pairs), length(results)))
} else {
  cat("Fresh start\n\n")
}

for (i in seq_len(nrow(strong))) {
  g     <- strong$gene[i]
  oid   <- strong$outcome_id[i]
  oname <- strong$outcome_name[i]
  pp4   <- strong$PP.H4[i]
  key   <- paste(g, oid, sep = "_")

  if (key %in% done_pairs) {
    cat(sprintf("[%d/%d] SKIP (done): %s ~ %s\n", i, nrow(strong), g, oname))
    next
  }

  # 2026-07-19 (Task123 D1 fix): MHC pairs are method-excluded from SuSiE, so do NOT compute PP.H4 for
  # them -- the raw susie_coloc_results artifact must not carry numeric SuSiE for AGER-HF/TRIM40-T2D.
  # Record a not_evaluable NA row (identical schema to the pre-fit NULL fallback below).
  if (isTRUE(strong$mhc_excluded[i])) {
    cat(sprintf("[%d/%d] MHC-excluded (method-mandated, no SuSiE): %s ~ %s\n", i, nrow(strong), g, oname))
    results[[key]] <- data.table(
      gene = g, outcome_id = oid, outcome_name = oname,
      PP.H4_abf = as.numeric(pp4), PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      n_snps_merged = NA_integer_, n_snps_ld = NA_integer_,
      n_signals_exp = NA_integer_, n_signals_out = NA_integer_,
      n_signal_pairs = NA_integer_, concordant_susie = NA,
      converged_exp = NA, converged_out = NA, max_iter_used = SUSIE_MAX_ITER
    )
    done_pairs <- c(done_pairs, key); next
  }

  cat(sprintf("\n[%d/%d] %s ~ %s  (PP.H4_abf=%.3f)\n",
              i, nrow(strong), g, oname, pp4))
  flush(stdout())

  # THE ROW COUNT IS TAKEN BEFORE , NOT AFTER, as in R/coloc.R:1321.  on an
  # absent gene gives zero rows, but  on a zero-row data.table returns ONE row of NA -- so
  # nrow(gene_info) was never 0 and the skip below was unreachable. The pair went on with NA
  # chrom/start/end. Every current suggestive gene is in the instrument table, so this is latent.
  .gi_all   <- instr[gene == g]
  gene_info <- if (nrow(.gi_all)) .gi_all[1] else .gi_all
  if (nrow(gene_info) == 0) {
    # A TRACEABLE ROW, NOT A SILENT DROP. This marked the pair done and omitted it, so the pair
    # simply disappeared from the output with nothing recording why. 26 writes a row for the same
    # condition; the schema below is 15's own not-evaluable fallback, used verbatim.
    cat(sprintf("[%d/%d] no instrument row (not evaluable): %s ~ %s\n", i, nrow(strong), g, oname))
    results[[key]] <- data.table(
      gene = g, outcome_id = oid, outcome_name = oname,
      PP.H4_abf = as.numeric(pp4), PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      n_snps_merged = NA_integer_, n_snps_ld = NA_integer_,
      n_signals_exp = NA_integer_, n_signals_out = NA_integer_,
      n_signal_pairs = NA_integer_, concordant_susie = NA,
      converged_exp = NA, converged_out = NA, max_iter_used = SUSIE_MAX_ITER
    )
    done_pairs <- c(done_pairs, key); next
  }

  res <- tryCatch(
    run_susie_pair(g, oid,
                    gene_info$gene_chrom, gene_info$gene_start, gene_info$gene_end,
                    rsid_map, pp4),
    error = function(e) {  # 2026-07-01 (closure B1): ANY error on a strong pair is fatal (99 D3 rejects NA); flag it. Legit-empty skips return NULL without erroring.
      message(sprintf("  [ERROR] %s ~ %s: %s", g, oid, conditionMessage(e)))
      assign(".CKM_DATA_FAILURE", sprintf("%s ~ %s SuSiE: %s", g, oid, conditionMessage(e)), envir = .GlobalEnv); NULL
    }
  )
  if (exists(".CKM_DATA_FAILURE", envir = .GlobalEnv))  # abort, do NOT record a NA SuSiE row
    stop(paste("15 ABORT -- data-foundation failure during SuSiE:", get(".CKM_DATA_FAILURE", envir = .GlobalEnv)))

  if (!is.null(res)) {
    res[, outcome_name := oname]
    results[[key]] <- res
  } else {
    results[[key]] <- data.table(
      gene = g, outcome_id = oid, outcome_name = oname,
      PP.H4_abf = as.numeric(pp4), PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
      n_snps_merged = NA_integer_, n_snps_ld = NA_integer_,
      n_signals_exp = NA_integer_, n_signals_out = NA_integer_,
      n_signal_pairs = NA_integer_, concordant_susie = NA,
      # NA, not FALSE: this fallback is reached when run_susie_pair() returned NULL BEFORE fitting (no pQTL
      # data, too few SNPs, no LD). No fit ran, so convergence is genuinely unknown. `FALSE` here would
      # assert a fit that failed to converge, which is a different and false claim. The post-fit exits
      # inside run_susie_pair() carry the real value and never reach this row.
      converged_exp = NA, converged_out = NA, max_iter_used = SUSIE_MAX_ITER
    )
  }
  done_pairs <- c(done_pairs, key)

  # Checkpoint every 5 pairs
  if (i %% 5 == 0 || i == nrow(strong)) {
    saveRDS(list(results = results, done_pairs = done_pairs, code_ver = .susie_code_ver), checkpoint_file)
    cat(sprintf("  [checkpoint saved: %d done, %d with results]\n",
                length(done_pairs), length(results)))
  }
  flush(stdout())
}

# ── Combine & save ─────────────────────────────────────────────────────────────
final_dt <- rbindlist(results, fill = TRUE, use.names = TRUE)
# P3-F1 (cold audit): expose the MHC-exclusion flag ON THE SAVED ARTIFACT so a DIRECT consumer of
# susie_coloc_results sees the method-mandated MHC exclusion (previously only on the in-memory `strong` table
# and dropped by the col_order restriction below). Downstream 12_rebuild independently sets not_evaluable/MHC_excluded.
final_dt[, mhc_excluded := if (exists("is_mhc_gene")) vapply(gene, function(x) isTRUE(is_mhc_gene(x)), logical(1)) else FALSE]
col_order <- c("gene", "outcome_name", "outcome_id",
               "PP.H4_abf", "PP.H4_susie", "PP.H3_susie",
               "concordant_susie", "mhc_excluded", "n_snps_merged", "n_snps_ld",
               "n_signals_exp", "n_signals_out", "n_signal_pairs",
               # Without these two the convergence recording added above is INERT: run_susie_pair() returns
               # them and this restriction throws them away, so susie_coloc_results would carry no
               # convergence record at all -- exactly the gap this change exists to close.
               # max_iter_used is here for the same reason and would be dropped the same way. A row
               # saying converged_exp = FALSE is not interpretable without the ceiling it ran under:
               # 1000 and 10000 produce an identical-looking row (review finding, 2026-08-08).
               "converged_exp", "converged_out", "max_iter_used")
col_order <- intersect(col_order, names(final_dt))
final_dt  <- final_dt[, ..col_order][order(-PP.H4_abf)]

# WRITTEN BESIDE, THEN PROMOTED (review finding, 2026-08-08). These two were saved straight to their
# final paths, RDS first and CSV second, so an interruption between them left a NEW rds beside the
# PREVIOUS csv -- and the two are read by different consumers: 26_susie_coloc_egfr.R takes the rds,
# while 13_rebuild_pair_evidence.R and 99_completeness_audit.R take the csv. Both files parse, both
# look complete, and nothing compares them, so a split-ceiling artifact is indistinguishable from a
# correct one. Writing both temporaries first means a death before the first rename leaves the two
# final files as they were -- consistent with each other -- rather than half-updated.
.rds_f <- file.path(out_dir, "susie_coloc_results.rds")
.csv_f <- file.path(out_dir, "susie_coloc_results.csv")
saveRDS(final_dt, paste0(.rds_f, ".tmp"))
fwrite(final_dt,  paste0(.csv_f, ".tmp"))
if (!all(file.rename(c(paste0(.rds_f, ".tmp"), paste0(.csv_f, ".tmp")), c(.rds_f, .csv_f))))
  stop("15 SuSiE: could not promote both result files; the pair on disk may now disagree -- inspect ",
       out_dir, " before using either.")

# ── Summary ────────────────────────────────────────────────────────────────────
cat("\n=== Phase 3-C SuSiE-coloc Summary ===\n")
cat(sprintf("Pairs attempted (non-eGFR):       %d\n", nrow(strong)))
cat(sprintf("Pairs with SuSiE result:           %d\n", sum(!is.na(final_dt$PP.H4_susie))))
cat(sprintf("Concordant (PP.H4_susie >= 0.80):  %d\n",
    sum(final_dt$concordant_susie, na.rm = TRUE)))
# THE DISCORDANT COUNT EXCLUDES NONCONVERGED FITS, like concordant_susie above. It counted every
# numeric PP.H4_susie < 0.80 regardless of convergence, so NCAN-CAD (neither side converged,
# PP.H4_susie = 1.777e-48) and ADM-AF (exposure did not converge) were reported as evidence AGAINST
# colocalisation. A fit that did not converge is not a negative result, and folding the two into the
# negative count is the same error the NA above was introduced to stop. The counts now partition the
# pairs that have a fit: concordant + discordant + not-evaluable.
.conv_ok <- final_dt$converged_exp %in% TRUE & final_dt$converged_out %in% TRUE
cat(sprintf("Discordant (PP.H4_abf>=0.8, PP.H4_susie<0.8, both fits converged): %d\n",
    sum(.conv_ok & !is.na(final_dt$PP.H4_susie) & final_dt$PP.H4_abf >= 0.80 &
        final_dt$PP.H4_susie < 0.80, na.rm = TRUE)))
cat(sprintf("Not evaluable (a fit did not converge):        %d\n",
    sum(!.conv_ok & !is.na(final_dt$PP.H4_susie), na.rm = TRUE)))

cat("\n--- Comparison table (top 20 by PP.H4_abf) ---\n")
print_cols <- intersect(c("gene","outcome_name","PP.H4_abf","PP.H4_susie",
                           "n_signals_exp","n_signals_out","concordant_susie"),
                        names(final_dt))
print(final_dt[seq_len(min(20, nrow(final_dt))), ..print_cols], digits = 3)

cat(sprintf("\nResults saved to: %s\n", out_dir))
