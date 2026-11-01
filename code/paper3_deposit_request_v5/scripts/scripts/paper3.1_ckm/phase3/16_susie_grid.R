# Phase 3-C: SuSiE Grid Sensitivity Analysis
# =============================================
# Runs susie_rss() ONCE per pair (L=10 and L=5), then re-filters credible sets
# with different min_abs_corr values via susie_get_cs(). This avoids recomputing
# LD matrices for each setting.
#
# Grid settings (per GPT / Wallace 2021 guidance):
#   strict:      L=10, min_abs_corr=0.5  (SuSiE default purity)
#   relaxed:     L=10, min_abs_corr=0.2
#   exploratory: L=10, min_abs_corr=0.1
#   robustness:  L=5,  min_abs_corr=0.2
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase3/16_susie_grid.R"

# CHUNK MODE IS REFUSED, not fingerprinted harder. The chunk protocol writes per-worker files and
# merges them, so the merge has to certify that every chunk came from the same run. Three rounds of
# review found the certificate incomplete each time -- first it covered only code, then it read a
# per-process scope so no two workers could agree, then it still omitted EUR.bed, the SG_COLOC
# artifact itself, and the susieR/ieugwasr/PLINK versions. That target has no boundary: every
# further input is one more thing the stamp does not cover, and a stamp that is trusted while
# incomplete is worse than none, because a mixed grid then looks complete.
#
# The canonical driver (_shared/localize/49_susie_grid_corrected.R) never sets SG_CHUNK_*; only this
# file's header ever documented them. The sibling made the same call already: cascade/11 clears
# PSW_CHUNK_N/I and PSW_MERGE outright and runs 25 serially. This grid is 56 pairs x 4 settings at
# roughly 40-140 s per pair, which is hours, not days, run serially. (PI decision 2026-08-03.)
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
if (nzchar(Sys.getenv("SG_CHUNK_N")) || nzchar(Sys.getenv("SG_CHUNK_I")) || nzchar(Sys.getenv("SG_MERGE")))
  stop("16_susie_grid: chunk mode is withdrawn. SG_CHUNK_N / SG_CHUNK_I / SG_MERGE are set, but the 
merge cannot certify that separately-computed chunks share one set of inputs -- three review rounds 
found the fingerprint incomplete, and a chunked grid that merges anyway looks complete while holding 
mixed computations. Run serially: unset those variables. The canonical driver 
(_shared/localize/49_susie_grid_corrected.R) already does not set them.")
# PLACED BEFORE ANY I/O. Sitting after the SG_COLOC read, this gate was unreachable on a bare run:
# the script died on the missing default input first, so an operator asking for chunk mode was told
# about a missing file instead of about the withdrawn mode.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(ieugwasr)
  library(susieR)
})

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))   # for is_coloc_pass() (explicit direct source; coloc also loads it)
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/strand_align.R"))   # root-cause LD sign-alignment (2026-07-17, editor ruling 2b)
# ── Source coloc.R (data extraction functions) ───────────────────────────────
COLOC_SOURCE_ONLY <- TRUE
suppressMessages(source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R")))
rm(COLOC_SOURCE_ONLY)

cat("\n=== Phase 3-C: SuSiE Grid Sensitivity ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── JWT check (skipped offline: local plink LD + local coords need no API) ───────
# ONE PREDICATE FOR ONE FLAG. This used nzchar() while the block refusing the online LD
# fallthrough tested for the exact string "1", so CKM_OFFLINE=TRUE skipped the preflight yet
# still permitted an online call. ckm_env_yes() (R/_phase1_set.R:20) is what R/00_config.R:1264
# itself uses, and nzchar() is TRUE even for "0" and "false".
if (!ckm_env_yes("CKM_OFFLINE")) {
  cat("Checking JWT validity...\n")
  jwt_ok <- tryCatch({
    resp <- ieugwasr::user()
    cat(sprintf("  JWT valid: user = %s\n\n", resp[["email"]]))
    TRUE
  }, error = function(e) { FALSE })
  if (!jwt_ok) stop("JWT invalid — update OPENGWAS_JWT in R/00_config.R")
} else cat("CKM_OFFLINE: skipping JWT preflight (local LD/coords).\n\n")

# ── Grid settings ──────────────────────────────────────────────────────────────
GRID <- list(
  strict      = list(L = 10, mac = 0.5),
  relaxed     = list(L = 10, mac = 0.2),
  exploratory = list(L = 10, mac = 0.1),
  robustness  = list(L = 5,  mac = 0.2)
)

# ── Load data ──────────────────────────────────────────────────────────────────
# A WRAPPER-SUPPLIED OUTPUT DIR MEANS A FRESH RUN, as in 15:42. Computed BEFORE the default below,
# so it distinguishes "a wrapper set this" from "this script defaulted it". The cache stores the LD
# matrix, the matched SNP set and n_out; its version token covers the three source files but not the
# rsID map, coordinate table, outcome GWAS or LD reference. The complete key already exists as
# coloc_inputs_key() in R/_coloc_inputs_key.R, but that file sources R/coloc.R, and by this point
# COLOC_SOURCE_ONLY has been rm()d -- calling it here would run the engine main(). Forcing fresh on
# the wrapper path removes the reliance instead of widening the token.
SG_FORCE_FRESH <- exists("SG_OUT")
if (!exists("SG_OUT")) SG_OUT <- file.path(RESULTS_DIR, "phase3_susie")   # A6 corrected wrapper overrides
out_dir   <- SG_OUT
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
# Cache version token = md5 of the strand-align helper (root-cause fix 2026-07-17, editor 2b, gap #2). The
# pair-prep cache stores LD keyed on gene_outcome only, with no version key; a cache written BEFORE this fix
# carries UNALIGNED R and, lacking this token, is treated as stale and rebuilt -- it can never bypass the
# alignment on re-run. Mirrors 29's code_md5-keyed checkpoint. Helper unreadable -> token NA -> rebuild all (fail safe).
# THE TOKEN MUST COVER THE CODE THAT PRODUCED THE CACHED VALUES, not one helper. A cache carries the
# LD matrix, the matched SNP set and n_out, all of which depend on this script and on R/coloc.R
# (liftover, merge, harmonisation), not only on the strand helper. Any of the three unreadable ->
# token NA -> rebuild everything (fail safe).
.sg_cache_ver <- tryCatch(paste(vapply(c("R/strand_align.R", "R/coloc.R", "scripts/paper3.1_ckm/phase3/16_susie_grid.R"),
                          function(.f) unname(tools::md5sum(file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"), .f))),
                          character(1)), collapse = ":"),
                         error = function(e) NA_character_)

# The chunk provenance token (.sg_chunk_ver) stood here and was DELETED with chunk mode. It existed
# only so a merge could certify that separately-computed chunks shared one set of inputs; with no
# merge there is nothing for it to certify. .sg_cache_ver below still guards the per-pair cache.

if (!exists("SG_COLOC")) SG_COLOC <- file.path(PHASE3_DIR, "coloc_results_final.rds")  # A6 override (corrected v6)
if (!exists("SG_INSTR")) SG_INSTR <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")  # A6 override (corrected 1682 (previous generation))
coloc_final <- readRDS(SG_COLOC)
instr       <- readRDS(SG_INSTR)
setDT(coloc_final); setDT(instr)

strong <- coloc_final[is_coloc_pass(PP.H4, PP.H3)][order(-PP.H4)]   # dual coloc rule (D1 2026-07-04)
# eGFR handling: default EXCLUDE by outcome_name (preserves live behavior, robust to id);
# corrected wrapper sets SG_KEEP_EGFR=TRUE to include eGFR (corrected susie = 52 incl eGFR).
if (!(exists("SG_KEEP_EGFR") && isTRUE(SG_KEEP_EGFR))) strong <- strong[outcome_name != "eGFR"]
# Codex backfill-review fix: honour config's PRE-SPECIFIED MHC SuSiE not_evaluable policy (00_config.R:142).
# is_mhc_gene() returns NA if gene coords unavailable -> fail loud (never silently treat an undetermined pair as non-MHC).
if (exists("is_mhc_gene")) {
  .mhc <- vapply(strong$gene, function(.g) is_mhc_gene(.g), logical(1))
  if (anyNA(.mhc)) stop("16_susie_grid: is_mhc_gene() NA (gene coords unavailable) for ", paste(unique(strong$gene[is.na(.mhc)]), collapse = ", "), " -- cannot apply MHC exclusion")
  if (any(.mhc)) { cat(sprintf("MHC not_evaluable (pre-specified, excluded from grid): %s\n", paste(unique(strong$gene[.mhc]), collapse = ", "))); strong <- strong[!.mhc] }
} else stop("16_susie_grid: is_mhc_gene() not defined (config not sourced?) -- cannot apply pre-specified MHC exclusion")
if (exists("SG_SMOKE") && isTRUE(SG_SMOKE))            # smoke: 1 non-eGFR + 1 eGFR pair
  strong <- rbind(head(strong[outcome_name != "eGFR"], 1), head(strong[outcome_name == "eGFR"], 1))
cat(sprintf("Strong hits: %d pairs\n\n", nrow(strong)))
.EXPECTED_KEYS <- strong[, paste(gene, outcome_id, sep = "_")]   # Codex backfill-review fix: capture the full non-MHC expected set BEFORE chunk/merge subsetting zeroes `strong` (merge-mode coverage validation)
# THE SAME SCOPE IN EVERY PROCESS. .sg_chunk_ver below fingerprints the run, so it must not depend
# on which slice of the run this process happens to hold. Reading `strong` inside that promise did:
# a worker forces it after `strong` has been cut to its chunk, and the merge forces it after `strong`
# has been reduced to zero rows -- so no two processes could ever agree and every chunk would be
# rejected. Captured HERE, next to .EXPECTED_KEYS, for the same reason: this is the last point at
# which every process still sees the whole post-MHC set.
.sg_scope <- list(genes = sort(unique(strong$gene)), outcome_ids = sort(unique(strong$outcome_id)))

rsid_map <- load_rsid_map()

# ── Pre-fetch hg38 coords ─────────────────────────────────────────────────────
needs_hg38 <- unique(strong$gene[
  strong$outcome_id %in% names(GWAS_META)[
    sapply(GWAS_META, function(m) !is.null(m$genome) && m$genome %in% c("hg38","hg38_hm"))
  ]
])
cat(sprintf("Pre-fetching hg38 coords for %d genes...\n", length(needs_hg38)))
gene_hg38_cache <- list()
for (g in needs_hg38) {
  for (attempt in 1:3) {
    coords <- tryCatch(get_gene_hg38_coords(g), error = function(e) NULL)
    if (!is.null(coords)) {
      gene_hg38_cache[[g]] <- coords
      cat(sprintf("  OK  %s\n", g)); break
    }
    Sys.sleep(3)
  }
  Sys.sleep(1)
}
cat(sprintf("Pre-cache done: %d / %d\n\n", length(gene_hg38_cache), length(needs_hg38)))

# Override hg38 lookup with cache
base_get_gene_hg38_coords <- get_gene_hg38_coords
get_gene_hg38_coords <- function(gene_name) {
  if (gene_name %in% names(gene_hg38_cache)) return(gene_hg38_cache[[gene_name]])
  Sys.sleep(1)
  coords <- tryCatch(base_get_gene_hg38_coords(gene_name), error = function(e) NULL)
  if (!is.null(coords)) gene_hg38_cache[[gene_name]] <<- coords
  coords
}

# ── LD matrix helper (copied from 15_susie_coloc.R) ──────────────────────────
get_ld_matrix_auto <- function(rsids, gene_chrom, max_snps_api = 400) {
  rsids <- unique(rsids[grepl("^rs[0-9]+$", rsids)])
  if (length(rsids) < 5) { cat("    [LD] < 5 rsIDs\n"); return(NULL) }
  use_local <- exists("PLINK_BIN") && !is.null(PLINK_BIN) && file.exists(PLINK_BIN) &&
               exists("LD_EUR_DIR") && !is.null(LD_EUR_DIR)
  if (use_local) {
    bfile <- file.path(LD_EUR_DIR, "EUR")
    if (file.exists(paste0(bfile, ".bim"))) {
      cat(sprintf("    [LD local] %d SNPs, chr%s\n", length(rsids), gene_chrom))
      # A LOCAL PLINK FAILURE IS INFRASTRUCTURE, NOT BIOLOGY -- same policy as 15:140 and 26:226,
      # which this helper claims to be copied from. Swallowing it here fell through to online LD,
      # i.e. a different LD reference and a 400-SNP cap, without the result saying so.
      ld <- tryCatch(ieugwasr::ld_matrix(rsids, bfile = bfile, plink_bin = PLINK_BIN),
                     error = function(e) fail_loud_data(sprintf("16_susie_grid local LD (PLINK) FAILED: %s -- LD infrastructure failure, not a biological skip", conditionMessage(e))))
      if (!is.null(ld)) { cat(sprintf("    [LD] %d x %d\n", nrow(ld), ncol(ld))); return(ld) }
    }
  }
  # A NULL LOCAL RESULT IS ALSO A LOCAL FAILURE. Making only the thrown exception fatal left the
  # other half of the same event falling through to online OpenGWAS LD -- a different reference and
  # a 400-SNP cap -- whenever offline mode was off, with nothing in the result saying the LD came
  # from elsewhere. 15 and 26 never switch reference. Reaching here with use_local TRUE means the
  # local panel was present but produced nothing.
  if (use_local)
    fail_loud_data(sprintf("16_susie_grid get_ld_matrix_auto: the local EUR panel returned no LD matrix for %d rsID(s) on chr%s -- LD infrastructure failure, not a biological skip; refusing to switch to online LD", length(rsids), gene_chrom))
  # Codex backfill-review fix: OFFLINE must NOT silently fall through to online OpenGWAS LD.
  if (ckm_env_yes("CKM_OFFLINE"))
    stop(sprintf("16_susie_grid get_ld_matrix_auto: CKM_OFFLINE=1 but local LD unavailable/failed (use_local=%s, %d rsIDs, chr%s) -- refusing online LD fallthrough", use_local, length(rsids), gene_chrom))
  if (length(rsids) > max_snps_api) rsids <- rsids[seq_len(max_snps_api)]
  Sys.sleep(1)
  tryCatch(ieugwasr::ld_matrix(rsids, pop = "EUR"), error = function(e) NULL)
}

# ── Phase 1: Prepare data + LD for each pair (cached) ────────────────────────
prepare_pair_data <- function(gene, outcome_id, gene_chrom, gene_start, gene_end,
                               rsid_map) {
  cache_file <- file.path(cache_dir, sprintf("%s_%s.rds", gene, gsub("[^a-zA-Z0-9]", "_", outcome_id)))
  if (file.exists(cache_file) && !isTRUE(SG_FORCE_FRESH) && !nzchar(Sys.getenv("CKM_PALINDROME_SINK", ""))) {   # skip cache when capturing palindromes: must rebuild so strand_align_R runs & writes the sink (cold re-review High: a cache hit skips capture)
    .cached <- readRDS(cache_file)
    if (is.list(.cached) && !is.null(.cached$cache_ver) && identical(.cached$cache_ver, .sg_cache_ver)) {
      cat(sprintf("  [cache] Loading %s ~ %s\n", gene, outcome_id))
      return(.cached)
    }
    cat(sprintf("  [cache] STALE (pre-strand-align) %s ~ %s -- rebuilding LD\n", gene, outcome_id))
  }

  meta   <- GWAS_META[[outcome_id]]
  genome <- meta$genome
  region_start <- max(1L, gene_start - CIS_WINDOW_KB * 1000L)
  region_end   <- gene_end + CIS_WINDOW_KB * 1000L

  cat(sprintf("  [prep] %s ~ %s | chr%s:%d-%d\n", gene, outcome_id, gene_chrom, region_start, region_end))

  # ── pQTL data ──
  pqtl_dt <- get_pqtl_region(gene, as.character(gene_chrom), region_start, region_end, rsid_map)
  if (is.null(pqtl_dt) || nrow(pqtl_dt) == 0) return(NULL)
  req_cols <- c("ALLELE0", "ALLELE1", "BETA", "SE", "INFO", "N")
  if (any(!req_cols %in% names(pqtl_dt))) return(NULL)
  pqtl_dt <- pqtl_dt[!is.na(INFO) & INFO >= INFO_MIN & !is.na(BETA) & !is.na(SE)]
  if (nrow(pqtl_dt) < 30) return(NULL)
  pqtl_dt[, pqtl_varbeta := SE^2]
  pqtl_dt[, ea_pqtl      := toupper(ALLELE1)]
  pqtl_dt[, oa_pqtl      := toupper(ALLELE0)]
  n_pqtl <- as.integer(pqtl_dt[["N"]][1])

  # ── Outcome data ──
  if (genome == "hg38") {
    # NO LIFTOVER: pqtl_dt$chrpos is chrom:GENPOS from the UKB-PPP archive and GENPOS is
    # already GRCh38 (verified against the archive itself, not against a comment).
    # liftover_to_hg38() has been DELETED -- it lifted GRCh38 a second time, displacing
    # positions by up to 864,716 bp, and every position-fallback match from it paired
    # two different variants.
    out_dt    <- get_finngen_region(outcome_id, gene)
    if (is.null(out_dt) || nrow(out_dt) < 20) return(NULL)
    chrpos38_col <- pqtl_dt$chrpos
    pqtl_sub <- data.table(rsid = pqtl_dt$rsid,
                           chrpos_pqtl = chrpos38_col,
                           pqtl_beta = pqtl_dt$BETA,
                           pqtl_varbeta = pqtl_dt$pqtl_varbeta,
                           ea_pqtl = pqtl_dt$ea_pqtl,
                           oa_pqtl = pqtl_dt$oa_pqtl)
    out_sub <- out_dt[, .(rsid_out, chrpos_out, ea_out = toupper(get(meta$ea_col)),
                          oa_out = toupper(get(meta$oa_col)),
                          out_beta = as.numeric(get(meta$beta_col)),
                          out_varbeta = as.numeric(get(meta$se_col))^2)]
  } else if (genome == "hg37") {
    # PROJECTED TO GRCh37: the outcome reports Pos_b37, pqtl_dt$chrpos is GRCh38 (see the hg38
    # branches above). Comparing them directly matches the WRONG variants where the builds coincide.
    pqtl_37 <- liftover_to_hg37(pqtl_dt)
    .psrc <- if (!is.null(pqtl_37)) pqtl_37 else pqtl_dt
    # LIFT FIRST, THEN EXTRACT. The window bounds are GRCh38 (they come from the gene coordinates
    # this build is keyed to) and the outcome file is GRCh37, so passing them straight through
    # queried the wrong interval -- on DNER the two differ by +864,815 bp and barely overlap.
    # Off-chromosome lifts are dropped for this outcome rather than aborting, which is the engine's
    # policy at R/coloc.R:1049; a whole-gene failure to lift is still fatal.
    .chr37  <- if (!is.null(pqtl_37)) sub(":.*", "", pqtl_37$chrpos37) else NA_character_
    .same37 <- .chr37 %in% as.character(gene_chrom)
    .off37  <- !is.na(.chr37) & !.same37
    if (is.null(pqtl_37) || !any(.same37)) return(NULL)
    if (any(.off37)) {
      cat(sprintf("    [liftover] %d pQTL(s) lifted off chr%s -- dropped for this outcome\n", sum(.off37), gene_chrom))
      pqtl_37$chrpos37[.off37] <- NA_character_
    }
    .pos37 <- as.integer(sub(".*:", "", pqtl_37$chrpos37[.same37]))
    out_dt <- get_nielsen_region(outcome_id, as.character(gene_chrom), min(.pos37), max(.pos37))
    if (is.null(out_dt) || nrow(out_dt) < 20) return(NULL)
    pqtl_sub <- data.table(rsid = .psrc$rsid,
                           chrpos_pqtl = pqtl_37$chrpos37,
                           pqtl_beta = .psrc$BETA, pqtl_varbeta = .psrc$pqtl_varbeta,
                           ea_pqtl = .psrc$ea_pqtl, oa_pqtl = .psrc$oa_pqtl)
    out_sub <- out_dt[, .(rsid_out, chrpos_out, ea_out = toupper(get(meta$ea_col)),
                          oa_out = toupper(get(meta$oa_col)),
                          out_beta = as.numeric(get(meta$beta_col)),
                          out_varbeta = as.numeric(get(meta$se_col))^2)]
  } else if (genome == "hg38_hm") {
    # NO LIFTOVER: pqtl_dt$chrpos is chrom:GENPOS from the UKB-PPP archive and GENPOS is
    # already GRCh38 (verified against the archive itself, not against a comment).
    # liftover_to_hg38() has been DELETED -- it lifted GRCh38 a second time, displacing
    # positions by up to 864,716 bp, and every position-fallback match from it paired
    # two different variants.
    out_dt    <- get_harmonized_hg38_region(outcome_id, gene)
    if (is.null(out_dt) || nrow(out_dt) < 20) return(NULL)
    pqtl_src <- pqtl_dt
    pqtl_sub <- data.table(rsid = pqtl_src$rsid,
                           chrpos_pqtl = pqtl_dt$chrpos,
                           pqtl_beta = pqtl_src$BETA, pqtl_varbeta = pqtl_src$pqtl_varbeta,
                           ea_pqtl = pqtl_src$ea_pqtl, oa_pqtl = pqtl_src$oa_pqtl)
    out_sub <- out_dt[, .(rsid_out, chrpos_out, ea_out = toupper(get(meta$ea_col)),
                          oa_out = toupper(get(meta$oa_col)),
                          out_beta = as.numeric(get(meta$beta_col)),
                          out_varbeta = as.numeric(get(meta$se_col))^2)]
  } else { return(NULL) }
  out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

  # ── Merge + harmonise ──
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

  if (nrow(common) < 20) return(NULL)
  # HARMONISE BEFORE DE-DUPLICATING, then re-key only the rsIDs that still collide -- the engine's
  # order (R/coloc.R:1148-1163). Dropping ambiguous groups first deleted variants that harmonisation
  # would have resolved to distinct outcome rows (1:1251285 G/A and G/T are both rs6657544).
  n_before <- nrow(common)
  common <- harmonize_alleles(common, "ea_pqtl", "oa_pqtl", "ea_out", "oa_out", "out_beta")
  common[duplicated(snp_id) | duplicated(snp_id, fromLast = TRUE),
         snp_id := paste(chrpos_pqtl, ea_pqtl, oa_pqtl, sep = ":")]
  common <- common[, if (.N == 1L) .SD, by = snp_id]
  if (nrow(common) < 20) return(NULL)

  # ── LD matrix ──
  common[, z_pqtl_abs := abs(pqtl_beta / sqrt(pqtl_varbeta))]
  common <- common[order(-z_pqtl_abs)]
  snps_for_ld <- common$snp_id[grepl("^rs[0-9]+$", common$snp_id)]
  if (length(snps_for_ld) < 10) return(NULL)

  ld_mat <- get_ld_matrix_auto(snps_for_ld, gene_chrom)
  if (is.null(ld_mat)) return(NULL)

  # ── Align to LD ──
  ld_snps_raw <- rownames(ld_mat)
  ld_snps_rs  <- sub("_[A-Za-z]+_[A-Za-z]+$", "", ld_snps_raw)
  idx_in_ld   <- match(ld_snps_rs, common$snp_id)
  keep_ld_pos <- which(!is.na(idx_in_ld))
  if (length(keep_ld_pos) < 10) return(NULL)

  common_ld <- common[idx_in_ld[keep_ld_pos]]
  ld_mat_f  <- ld_mat[keep_ld_pos, keep_ld_pos, drop = FALSE]
  rownames(ld_mat_f) <- colnames(ld_mat_f) <- ld_snps_rs[keep_ld_pos]

  # STRAND-ALIGN R to the pQTL effect allele (root-cause fix 2026-07-17, editor 2b). Must run BEFORE the cache
  # save below so the cached LD carries the corrected orientation. No-op for aligned SNPs; palindromes kept
  # (verified by frequency); allele_mismatch dropped.
  .sa <- strand_align_R(common_ld, ld_snps_raw[keep_ld_pos], ld_mat_f, label = paste0(gene, "~", outcome_id))
  common_ld <- .sa$common_ld; ld_mat_f <- .sa$R
  if (nrow(common_ld) < 10) { cat("    [skip] < 10 after strand-align drop\n"); return(NULL) }   # re-check gate after drops (R2 Medium)

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

  result <- list(z_exp = z_exp, z_out = z_out, ld_mat = ld_mat_f,
                 n_pqtl = n_pqtl, n_out = n_out, n_snps = nrow(common_ld),
                 cache_ver = .sg_cache_ver)   # stamp for stale-cache invalidation (editor 2b, gap #2)
  saveRDS(result, cache_file)
  cat(sprintf("    cached: %d SNPs\n", nrow(common_ld)))
  result
}

# ── Phase 2: Run SuSiE grid on cached data ──────────────────────────────────
run_susie_grid_pair <- function(dat, gene, outcome_id, pp4_abf) {
  rows <- list()

  # Run susie_rss ONCE per L value, with min_abs_corr=0 (keep all CS)
  fits <- list()
  for (L_val in c(10, 5)) {
    key <- paste0("L", L_val)
    # max_iter = 1000 (2026-07-16). The default is 100 -- a susieR default, never a choice this study made
    # and never stated anywhere in it. coloc::runsusie, the wrapper the package provides for exactly this
    # task, defaults to repeat_until_convergence = TRUE because its authors treat non-convergence at the
    # default limit as a computational artifact, not a result; this engine bypasses that wrapper and calls
    # susie_rss directly, at the default, with no retry. Raising it removes an unspecified stopping rule.
    # IT IS FORWARDED, NOT SWALLOWED -- and that had to be checked, because `check_LD` on this same line is
    # NOT: max_iter is absent from formals(susie_rss) exactly as check_LD is, but susie_suff_stat (which
    # susie_rss forwards `...` to) HAS max_iter as a formal and check_LD is a formal of neither. Verified
    # empirically 2026-07-16, which is the only proof that counts: max_iter=1 -> converged FALSE,
    # max_iter=1000 -> converged TRUE. And a fit that already converged is UNCHANGED by a higher ceiling
    # (identical() on the returned object), so this can only move fits that did not converge at 100.
    S1 <- tryCatch(susie_rss(z = dat$z_exp, R = dat$ld_mat, n = dat$n_pqtl,
                              L = L_val, min_abs_corr = 0, verbose = FALSE, check_LD = FALSE, max_iter = 1000),
                   error = function(e) NULL)
    S2 <- tryCatch(susie_rss(z = dat$z_out, R = dat$ld_mat, n = dat$n_out,
                              L = L_val, min_abs_corr = 0, verbose = FALSE, check_LD = FALSE, max_iter = 1000),
                   error = function(e) NULL)
    fits[[key]] <- list(S1 = S1, S2 = S2)
  }

  # For each grid setting, re-filter CS and run coloc.susie
  for (setting_name in names(GRID)) {
    setting <- GRID[[setting_name]]
    L_key <- paste0("L", setting$L)
    S1_raw <- fits[[L_key]]$S1
    S2_raw <- fits[[L_key]]$S2

    # Codex backfill-review fix: a SuSiE fit EXCEPTION on an authoritative non-MHC strong pair must throw, not become
    # an all-NA "fit_failed" row classified downstream as weak/unstable evidence.
    if (is.null(S1_raw) || is.null(S2_raw))
      stop(sprintf("16_susie_grid: susie_rss() FAILED for %s ~ %s [%s] (exp=%s, out=%s) -- refusing all-NA rows for a strong hit",
                   gene, outcome_id, setting_name, if (is.null(S1_raw)) "NULL" else "ok", if (is.null(S2_raw)) "NULL" else "ok"))

    # Re-filter credible sets with this setting's min_abs_corr
    S1 <- S1_raw; S2 <- S2_raw
    S1$sets <- susie_get_cs(S1, Xcorr = dat$ld_mat, min_abs_corr = setting$mac)
    S2$sets <- susie_get_cs(S2, Xcorr = dat$ld_mat, min_abs_corr = setting$mac)

    n_cs_exp <- length(S1$sets$cs)
    n_cs_out <- length(S2$sets$cs)

    if (n_cs_exp == 0 || n_cs_out == 0) {
      rows[[setting_name]] <- data.table(
        gene = gene, outcome_id = outcome_id, setting = setting_name,
        L = setting$L, min_abs_corr = setting$mac,
        PP.H4_abf = as.numeric(pp4_abf),
        PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
        PP.H1_susie = NA_real_, PP.H2_susie = NA_real_,
        n_snps = dat$n_snps, n_cs_exp = n_cs_exp, n_cs_out = n_cs_out,
        n_signal_pairs = 0L,
        susie_status = ifelse(n_cs_out == 0, "no_cs_outcome", "no_cs_exposure"),
        # Both fits RAN here (a fit with no credible set is not a failed fit), so their convergence is
        # meaningful and is recorded on every branch -- a row that omits it would read as "unknown" when
        # the answer is known.
        converged_exp = isTRUE(S1$converged), converged_out = isTRUE(S2$converged)
      )
      next
    }

    cs_res <- tryCatch(coloc::coloc.susie(S1, S2), error = function(e) NULL)
    # Codex backfill-review fix: a coloc.susie EXCEPTION on a strong pair must throw; an empty-but-SUCCESSFUL summary
    # (ran, no signal pair) stays a legitimate "coloc_empty" status.
    if (is.null(cs_res))
      stop(sprintf("16_susie_grid: coloc.susie() FAILED for %s ~ %s [%s] -- refusing all-NA rows for a strong hit", gene, outcome_id, setting_name))
    # THE PACKAGE'S NO-RESULT SENTINEL IS A BARE data.table(nsnps = NA), WITH NO `summary` COMPONENT.
    # Verified against the installed coloc 5.2.3: coloc.susie() has `return(data.table(nsnps = NA))` on
    # its no-credible-set path. A guard that treats any missing `summary` as fatal therefore aborts a
    # legitimate result -- the defect this guard was written to prevent, committed inside the guard.
    # So the sentinel is recognised first and takes the empty-result path; only a shape that is neither
    # the sentinel nor a summary-bearing object is a schema change and stops.
    .cs_sentinel <- is.data.frame(cs_res) && identical(names(cs_res), "nsnps")
    if (!.cs_sentinel && !("summary" %in% names(cs_res)))
      stop(sprintf("16_susie_grid: coloc.susie returned an object that is neither its no-CS sentinel nor a summary-bearing result for %s ~ %s [%s] (class %s, names %s) -- the result schema changed; treating it as an empty colocalisation would discard whatever it does contain.", gene, outcome_id, setting_name, class(cs_res)[1], paste(names(cs_res), collapse = ", ")))
    if (nrow(as.data.table(cs_res$summary)) == 0) {
      rows[[setting_name]] <- data.table(
        gene = gene, outcome_id = outcome_id, setting = setting_name,
        L = setting$L, min_abs_corr = setting$mac,
        PP.H4_abf = as.numeric(pp4_abf),
        PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
        PP.H1_susie = NA_real_, PP.H2_susie = NA_real_,
        n_snps = dat$n_snps, n_cs_exp = n_cs_exp, n_cs_out = n_cs_out,
        n_signal_pairs = 0L, susie_status = "coloc_empty",
        converged_exp = isTRUE(S1$converged), converged_out = isTRUE(S2$converged)
      )
      next
    }

    sdt <- as.data.table(cs_res$summary)
    best <- sdt[which.max(PP.H4.abf)]

    rows[[setting_name]] <- data.table(
      gene = gene, outcome_id = outcome_id, setting = setting_name,
      L = setting$L, min_abs_corr = setting$mac,
      # NOT ROUNDED. The supportive/discordant classification below thresholds these directly, so
      # rounding could change the class; and 26 writes unrounded rows into the same combined
      # artifact. Same change as 25 and 15.
      PP.H4_abf = as.numeric(pp4_abf),
      PP.H4_susie = as.numeric(best$PP.H4.abf),
      PP.H3_susie = as.numeric(best$PP.H3.abf),
      PP.H1_susie = as.numeric(best$PP.H1.abf),
      PP.H2_susie = as.numeric(best$PP.H2.abf),
      n_snps = dat$n_snps, n_cs_exp = n_cs_exp, n_cs_out = n_cs_out,
      n_signal_pairs = nrow(sdt), susie_status = "ok",
      converged_exp = isTRUE(S1$converged), converged_out = isTRUE(S2$converged)
    )
  }
  rbindlist(rows)
}

# ── Main loop ──────────────────────────────────────────────────────────────────
all_results <- list()

for (i in seq_len(nrow(strong))) {
  g     <- strong$gene[i]
  oid   <- strong$outcome_id[i]
  oname <- strong$outcome_name[i]
  pp4   <- strong$PP.H4[i]
  key   <- paste(g, oid, sep = "_")

  cat(sprintf("\n[%d/%d] %s ~ %s (PP.H4_abf=%.3f)\n", i, nrow(strong), g, oname, pp4))

  # THE ROW COUNT IS TAKEN BEFORE , NOT AFTER, as in R/coloc.R:1321.  on an
  # absent gene gives zero rows, but  on a zero-row data.table returns ONE row of NA -- so
  # nrow(gene_info) was never 0 and the skip below was unreachable. The pair went on with NA
  # chrom/start/end. Every current suggestive gene is in the instrument table, so this is latent.
  .gi_all   <- instr[gene == g]
  gene_info <- if (nrow(.gi_all)) .gi_all[1] else .gi_all
  # Codex backfill-review fix: a strong pair MUST have an instrument -- fail loud, never silently skip an authoritative hit
  if (nrow(gene_info) == 0) stop(sprintf("16_susie_grid: strong pair %s ~ %s has NO instrument row", g, oid))

  dat <- tryCatch(
    prepare_pair_data(g, oid, gene_info$gene_chrom, gene_info$gene_start, gene_info$gene_end, rsid_map),
    error = function(e) { message("  [ERROR] ", conditionMessage(e)); NULL }
  )

  # Codex backfill-review fix: MHC pairs are already excluded upstream (pre-specified not_evaluable). A NON-MHC strong
  # hit that cannot be prepped (missing pQTL/outcome/LD or low coverage) is a genuine coverage failure -> fail loud,
  # rather than emit 4 all-NA "data_prep_failed" rows that classify downstream as weak/unstable evidence for a strong hit.
  if (is.null(dat)) stop(sprintf("16_susie_grid: data prep FAILED for non-MHC strong pair %s ~ %s (missing pQTL/outcome/LD or <min SNP coverage) -- refusing to emit all-NA rows for an authoritative strong hit", g, oid))

  cat(sprintf("  Running SuSiE grid (4 settings)...\n"))
  # Codex backfill-review fix: do NOT swallow grid errors to NULL -- internal fit/coloc exceptions must propagate as fail-loud stops
  grid_res <- run_susie_grid_pair(dat, g, oid, pp4)

  # Codex backfill-review fix: a grid error must NOT silently omit a strong pair (all_results[[key]]<-NULL would drop it)
  if (is.null(grid_res)) stop(sprintf("16_susie_grid: SuSiE grid FAILED for %s ~ %s -- refusing to silently omit a strong pair", g, oid))
  grid_res[, outcome_name := oname]
  for (j in seq_len(nrow(grid_res))) {
    sn <- grid_res$setting[j]
    cat(sprintf("    %s: cs(%d,%d) H4=%.3f H3=%.3f [%s]\n",
                sn,
                ifelse(is.na(grid_res$n_cs_exp[j]), 0L, grid_res$n_cs_exp[j]),
                ifelse(is.na(grid_res$n_cs_out[j]), 0L, grid_res$n_cs_out[j]),
                ifelse(is.na(grid_res$PP.H4_susie[j]), 0, grid_res$PP.H4_susie[j]),
                ifelse(is.na(grid_res$PP.H3_susie[j]), 0, grid_res$PP.H3_susie[j]),
                grid_res$susie_status[j]))
  }
  all_results[[key]] <- grid_res
  flush(stdout())
}

# Codex backfill-review fix: assert every strong pair produced a full grid (no silent omission via next/NULL)
.exp_keys <- strong[, paste(gene, outcome_id, sep = "_")]
if (!setequal(names(all_results), .exp_keys)) stop(sprintf("16_susie_grid: grid coverage mismatch -- missing pair(s): %s", paste(setdiff(.exp_keys, names(all_results)), collapse = ", ")))
.badn <- names(all_results)[vapply(all_results, function(x) nrow(x) != length(GRID), logical(1))]
if (length(.badn)) stop(sprintf("16_susie_grid: pair(s) without %d grid rows: %s", length(GRID), paste(.badn, collapse = ", ")))

# ── Combine & classify ────────────────────────────────────────────────────────
final <- rbindlist(all_results, fill = TRUE, use.names = TRUE)
# The chunk-write and chunk-merge blocks stood here and were DELETED with chunk mode itself; see the
# note at the SG_CHUNK_* refusal above. `final` is now always the serial rbindlist, and the coverage
# assertions the merge used to perform are the per-pair ones immediately above it.

# Classification per pair: based on strict + relaxed settings
classify_pair <- function(dt) {
  # A NONCONVERGED FIT'S POSTERIOR IS NOT EVIDENCE, and this grid was thresholding it as if it were:
  # NCAN-CAD (both sides converged = FALSE) came out "discordant" under all four settings on
  # PP.H3_susie = 1.000, and ADM-AF (exposure FALSE) came out inconclusive. phase4/12:545 already
  # refuses such a pair BEFORE reading this label, so the authoritative 43/28/1/4 split never used
  # these values -- but this CSV is read directly by reporting/susie_grid.R and by the supplement.
  # NA_character_ rather than a new "not_evaluable" level: 12:548 already tests `!is.na(classification)`
  # and falls through to its own convergence gate, whereas a new level would land in
  # reporting/susie_grid.R's unlabelled bucket and be shown as if it were inconclusive.
  #
  # GATED ON THE TWO SETTINGS THIS FUNCTION READS, not by filtering `dt` first. Filtering was measured
  # to mislabel ADM-AF: only its `robustness` fit converged, so dropping the nonconverged rows left
  # `strict` and `relaxed` empty, every threshold fell through, and the pair came out "unstable" --
  # which asserts that the settings disagreed, a different and false claim. Both of the pairs at issue
  # now return NA, matching what phase4/12:545 independently decides from the relaxed row.
  strict  <- dt[setting == "strict"]
  relaxed <- dt[setting == "relaxed"]
  .converged <- function(x) nrow(x) > 0 && isTRUE(x$converged_exp[1]) && isTRUE(x$converged_out[1])
  if (!.converged(strict) && !.converged(relaxed)) return(NA_character_)

  # Supportive: H4 >= 0.80 in strict OR relaxed
  if (nrow(strict) && !is.na(strict$PP.H4_susie) && strict$PP.H4_susie >= 0.80) return("supportive")
  if (nrow(relaxed) && !is.na(relaxed$PP.H4_susie) && relaxed$PP.H4_susie >= 0.80) return("supportive")

  # Discordant: H3 > 0.50 in strict or relaxed with credible CS
  if (nrow(strict) && !is.na(strict$PP.H3_susie) && strict$PP.H3_susie > 0.50) return("discordant")
  if (nrow(relaxed) && !is.na(relaxed$PP.H3_susie) && relaxed$PP.H3_susie > 0.50) return("discordant")

  # Inconclusive: no CS in outcome under strict AND relaxed
  strict_nocs  <- nrow(strict)  && strict$susie_status  == "no_cs_outcome"
  relaxed_nocs <- nrow(relaxed) && relaxed$susie_status == "no_cs_outcome"
  if (strict_nocs && relaxed_nocs) return("inconclusive")

  # Unstable: has result but classification differs across settings
  return("unstable")
}

pair_class <- final[, .(classification = classify_pair(.SD)),
                    by = .(gene, outcome_id)]
final <- merge(final, pair_class, by = c("gene", "outcome_id"), all.x = TRUE)

# ── Save ───────────────────────────────────────────────────────────────────────
# A SMOKE RUN MUST NOT REPLACE THE REAL GRID. SG_SMOKE cuts `strong` to one pair before the
# coverage assertions, so every assertion passes against the reduced set and the four resulting
# rows were written to the canonical paths with nothing marking them as a smoke. Run after a full
# grid, that silently replaced 224 rows with 4, and downstream counts, classifications and figures
# would follow the four. 26 and 25 already take this shape for their own caps: compute, print,
# do not save.
if (exists("SG_SMOKE") && isTRUE(SG_SMOKE)) {
  cat(sprintf("\n[SMOKE] SG_SMOKE set: grid ran on %d pair(s) x %d setting(s) = %d row(s); NOT saving over %s\n",
              uniqueN(final[, paste(gene, outcome_id)]), length(GRID), nrow(final), out_dir))
  print(final[, intersect(c("gene", "outcome_id", "setting", "PP.H4_susie", "susie_status"), names(final)), with = FALSE])
  q(save = "no", status = 0)
}
fwrite(final, file.path(out_dir, "susie_grid_results.csv"))
saveRDS(final, file.path(out_dir, "susie_grid_results.rds"))

# ── Summary ────────────────────────────────────────────────────────────────────
cat("\n\n=== Classification Summary ===\n")
class_tbl <- pair_class[, .N, by = classification][order(-N)]
for (j in seq_len(nrow(class_tbl))) {
  cat(sprintf("  %s: %d\n", class_tbl$classification[j], class_tbl$N[j]))
}

cat("\n--- Per-setting overview ---\n")
for (sn in names(GRID)) {
  sub <- final[setting == sn]
  n_ok   <- sum(sub$susie_status == "ok", na.rm = TRUE)
  n_nocs <- sum(sub$susie_status == "no_cs_outcome", na.rm = TRUE)
  n_h4   <- sum(!is.na(sub$PP.H4_susie) & sub$PP.H4_susie >= 0.80, na.rm = TRUE)
  cat(sprintf("  %s (L=%d, mac=%.1f): %d ok, %d no_cs_out, %d H4>=0.80\n",
              sn, GRID[[sn]]$L, GRID[[sn]]$mac, n_ok, n_nocs, n_h4))
}

cat(sprintf("\nResults: %s\n", file.path(out_dir, "susie_grid_results.csv")))
cat("Done.\n")
