# Phase 3: Colocalization Validation (v5)
# =========================================================
# v5 vs v4 主要改動：
#   1. 輸入: mr_suggestive.rds（807 pairs）替換 mr_significant_fdr.rds（38 pairs）
#   2. rsID 作為主要合併鍵（primary merge by rsID, fallback by position）
#   3. MAF filter（pQTL A1FREQ 0.01-0.99）
#   4. 新增 MEGASTROKE（hg38 harmonized）、CKDGen（hg38 harmonized）outcome 類型
#   5. FinnGen 和 Nielsen AF 亦提取 rsID（rsids / rs_dbSNP147 欄位）
#   6. cis_rsid_map.tsv 由 scripts/_shared/phase1/05_build_cis_rsid_map.py 預先建立
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/R/coloc.R"

# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
suppressPackageStartupMessages(library(data.table))
check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' not installed.", pkg))
  library(pkg, character.only = TRUE)
}
check_pkg("coloc")

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))

cat("=== Phase 3: Colocalization (v5) ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
flush(stdout())

# ── 全域設定 ──────────────────────────────────────────────────────────────────
CACHE_DIR <- file.path(BASE_DIR, "cache")
MAF_MIN   <- 0.01   # pQTL MAF filter

# ── GWAS 登錄表 ─────────────────────────────────────────────────────────────
# genome 類型：
#   "hg38"    = FinnGen（含 rsids 欄）
#   "hg37"    = Nielsen AF（含 rs_dbSNP147 欄）
#   "hg38_hm" = GWAS Catalog harmonized，hg38，含 hm_rsid 欄
#   "rsid_only" = 無本地完整檔，skip
GWAS_META <- list(
  "FinnGen_R10_I9_CHD" = list(
    file = "finngen_R10_I9_CHD.gz", genome = "hg38",
    chrom_col = "#chrom", pos_col = "pos", rsid_col = "rsids",
    ea_col = "alt", oa_col = "ref",
    beta_col = "beta", se_col = "sebeta",
    sep = "\t", type = "cc", n = 412181L, ncase = 46959L
  ),
  "FinnGen_R10_I9_HEARTFAIL" = list(
    file = "finngen_R10_I9_HEARTFAIL.gz", genome = "hg38",
    chrom_col = "#chrom", pos_col = "pos", rsid_col = "rsids",
    ea_col = "alt", oa_col = "ref",
    beta_col = "beta", se_col = "sebeta",
    sep = "\t", type = "cc", n = 412181L, ncase = 29672L
  ),
  "FinnGen_R10_T2D" = list(
    file = "finngen_R10_T2D.gz", genome = "hg38",
    chrom_col = "#chrom", pos_col = "pos", rsid_col = "rsids",
    ea_col = "alt", oa_col = "ref",
    beta_col = "beta", se_col = "sebeta",
    sep = "\t", type = "cc", n = 400197L, ncase = 65085L
  ),
  "FinnGen_R10_NAFLD" = list(
    file = "finngen_R10_NAFLD.gz", genome = "hg38",
    chrom_col = "#chrom", pos_col = "pos", rsid_col = "rsids",
    ea_col = "alt", oa_col = "ref",
    beta_col = "beta", se_col = "sebeta",
    sep = "\t", type = "cc", n = 412181L, ncase = 2568L
  ),
  "FinnGen_R10_I9_AF" = list(   # AF: Nielsen(UKB)->FinnGen R10 I9_AF(no UKB) 2026-06-30; key legacy pending rename
    file = "finngen_R10_I9_AF.gz", genome = "hg38",
    chrom_col = "#chrom", pos_col = "pos", rsid_col = "rsids",
    ea_col = "alt", oa_col = "ref",
    beta_col = "beta", se_col = "sebeta",
    sep = "\t", type = "cc", n = 261395L, ncase = 50743L
  ),
  "MEGASTROKE_EUR_GCST006906" = list(   # MEGASTROKE 歐洲裔（hg38 harmonized）
    file = "MEGASTROKE_EUR_harmonised_hg38.h.tsv.gz", genome = "hg38_hm",
    chrom_col = "hm_chrom", pos_col = "hm_pos", rsid_col = "hm_rsid",
    ea_col = "hm_effect_allele", oa_col = "hm_other_allele",
    beta_col = "hm_beta", se_col = "standard_error",
    sep = "\t", type = "cc", n = 446696L, ncase = 40585L
  ),
  "Wuttke2019_CKD_EA" = list(         # CKD: Pattaro(mixed,hg38)->Wuttke 2019 EA(EUR,hg37) 2026-06-30; key legacy pending rename
    file = "CKDGen_CKD_Wuttke2019_EA.txt.gz", genome = "hg37",
    chrom_col = "Chr", pos_col = "Pos_b37", rsid_col = "RSID",
    ea_col = "Allele1", oa_col = "Allele2",
    beta_col = "Effect", se_col = "StdErr",
    sep = " ", type = "cc", n = 480698L, ncase = 41395L
  ),
  # LEGACY/live eGFR entry. The CORRECTED cascade handles eGFR via localize/06 (Wuttke2019_eGFR_EA + CKDGen file),
  # NOT via this GWAS_META entry. It is UNUSED in the corrected coloc: the non-eGFR path (localize/21) excludes
  # eGFR by outcome_id, and 06 uses its own file. rsid_only -> skip even if ever hit. Kept for live back-compat.
  "ebi-a-GCST003372" = list(   # eGFR (Pattaro 2016) — LEGACY stub, UNUSED in corrected coloc (eGFR runs via localize/06)
    file = NA, genome = "rsid_only",
    type = "quant", n = 133413L
  ),
  "Wuttke2019_eGFR_EA" = list( # eGFR Wuttke 2019 EA (corrected; mirrors 02 + both Steiger; consistency entry) 2026-06-30
    file = "CKDGen_eGFR_EA_hg37.txt.gz", genome = "hg37",
    chrom_col = "Chr", pos_col = "Pos_b37", rsid_col = "RSID",
    ea_col = "Allele1", oa_col = "Allele2",
    beta_col = "Effect", se_col = "StdErr",
    sep = " ", type = "quant", n = 567460L,
    # THE PHENOTYPE SD, because coloc's quantitative prior is 0.15 * sdY and N plays no part in it
    # once beta and varbeta are given (coloc 5.2.3 approx.bf.estimates). Every quantitative dataset
    # here used to declare sdY = 1, which asserts a variance-one phenotype. These are CKDGen effects
    # on natural-log eGFR, not a standardised trait, so that assertion was wrong by a factor of five
    # and every eGFR SNP's Bayes factor with it -- for rs4970384 the log ABF is -3.790 at sdY = 1 and
    # -2.093 at the true scale. Estimated over the whole file (8,864,521 variants with finite SE, N
    # and frequency) by coloc's own sdY.est relationship, using each SNP's n_total_sum rather than a
    # single median N (which gives 0.18079 and is the wrong quantity). It is stable: SD 0.0000157
    # across 100 deterministic partitions, 0.18290-0.18394 across chromosomes, and unchanged to four
    # decimals when rare variants are excluded. It belongs to THIS file -- a different eGFR release,
    # or a differently transformed phenotype, needs its own estimate.
    sdY = 0.183122888
  )
)

# 2026-07-01: N/ncase for the 8 config outcomes are AUTHORITATIVELY (re)set from config's SINGLE SOURCE
# (OUTCOME_N_TOTAL / OUTCOME_N_CASE, guarded by config's stopifnot), overwriting the inline literals above so
# GWAS_META can never drift from config / 01_mr / 09_steiger. coloc.abf uses one outcome-level N (not per-SNP).
# Entries NOT in config OUTCOMES (the legacy Pattaro 'ebi-a-GCST003372' eGFR stub, UNUSED in the corrected
# cascade) keep their own literal.
.gm_id2nm <- setNames(names(OUTCOMES), unname(OUTCOMES))
for (.gm_oid in names(GWAS_META)) { .gm_nm <- unname(.gm_id2nm[.gm_oid]); if (is.na(.gm_nm)) next   # `[` not `[[`: non-config ids (Pattaro stub) -> NA, skip (not a throw)
  GWAS_META[[.gm_oid]]$n <- as.integer(OUTCOME_N_TOTAL[[.gm_nm]])
  if (isTRUE(GWAS_META[[.gm_oid]]$type == "cc")) {
    if (is.na(OUTCOME_N_CASE[[.gm_nm]])) stop(sprintf("03 coloc: config OUTCOME_N_CASE is NA for case-control outcome %s", .gm_oid))
    GWAS_META[[.gm_oid]]$ncase <- as.integer(OUTCOME_N_CASE[[.gm_nm]]) } }

# ── 載入 cis rsID map ─────────────────────────────────────────────────────────
load_rsid_map <- function() {
  f <- file.path(CACHE_DIR, "cis_rsid_map.tsv")
  if (!file.exists(f)) {
    cat("[warn] cis_rsid_map.tsv not found — rsID primary matching disabled\n")
    cat("       Run scripts/_shared/phase1/05_build_cis_rsid_map.py first\n\n")
    return(NULL)
  }
  dt <- fread(f, colClasses = list(character = "chrom"))
  # KEYED ON THE ALLELES TOO. The map file carries allele0/allele1 and this threw them away, so two
  # variants at one position collapsed to one rsID. At 1:2102322 the map holds T/G -> rs201135834 and
  # T/TG -> rs34905415, and both occur in current cis members: the indel was joined to the SNV's
  # rsID, its position was then marked as matched so the position fallback could not rescue it, and
  # harmonisation rejected T/TG against T/G. The row was lost, not merely mislabelled.
  dt[, lookup_key := paste(chrom, pos, toupper(allele0), toupper(allele1), sep = ":")]

  # 去重：同一 rsid 可能對應多個位置（dbSNP legacy artifact）
  # 排序後再 dedup → 選擇可重現（非隨機）
  # rule: 優先 chrom 最小、pos 最小（canonical primary mapping）
  n_raw <- nrow(dt)
  dt <- dt[order(as.integer(chrom), pos)]
  # ONE RSID -> ONE POSITION, BUT EVERY ALLELE PAIR AT THAT POSITION. `.SD[1]` kept a single row per
  # rsID, which collapsed the legacy-multiple-POSITIONS case it was written for AND the distinct
  # variants case it was not: 1:1251285 carries G/A and G/T under one rsID (rs6657544), and both are
  # in NADK's current archive with both outcome rows present in FinnGen T2D. One of them was silently
  # dropped before the allele-aware lookup_key below could distinguish them -- so making that key
  # allele-aware, as was done earlier, changed nothing while this line ran first.
  dt <- dt[, .SD[chrom == chrom[1L] & pos == pos[1L]], by = rsid]   # one rsid -> one position, all its alleles
  dt <- dt[, .SD[1], by = lookup_key]  # one position+alleles → one rsid
  n_dedup <- nrow(dt)
  if (n_dedup < n_raw)
    cat(sprintf("[rsID map] Deduplicated %s → %s entries (removed %d duplicates)\n",
                formatC(n_raw, big.mark = ","),
                formatC(n_dedup, big.mark = ","),
                n_raw - n_dedup))

  cat(sprintf("[rsID map] Loaded %s entries\n\n", formatC(nrow(dt), big.mark = ",")))
  dt
}

# ── 載入 Phase 2 suggestive pairs（812）──────────────────────────────────────
load_phase2_suggestive <- function() {
  sug_file <- file.path(PHASE2_DIR, "mr_suggestive.rds")
  if (!file.exists(sug_file))
    stop("Phase 2 suggestive results not found: ", sug_file)
  pairs <- readRDS(sug_file)
  cat(sprintf("Phase 2 suggestive pairs: %d (%d unique proteins)\n",
              nrow(pairs), uniqueN(pairs$gene)))
  pairs
}

# ── 互補鹼基 ────────────────────────────────────────────────────────────────
complement_allele <- function(a) chartr("ACGTacgt", "TGCAtgca", toupper(a))
is_palindromic    <- function(a1, a2) complement_allele(a1) == toupper(a2)

# ── Allele Harmonization（含 strand flip）────────────────────────────────────
# CAN THESE TWO RECORDS DESCRIBE THE SAME VARIANT? Factored out of harmonize_alleles() because the
# rsID merge has to ask the same question BEFORE it declares a coordinate matched: an rsID hit whose
# alleles cannot harmonise used to block the position fallback, so the compatible record at the same
# coordinate became unreachable and the pair lost SNPs. Two copies of the predicate would drift, and
# harmonize_alleles() cannot be reused there because it also flips beta and so is not idempotent.
allele_compatible <- function(ea_p, oa_p, ea_o, oa_o) {
  ep <- toupper(ea_p); op <- toupper(oa_p); eo <- toupper(ea_o); oo <- toupper(oa_o)
  strand_snv <- nchar(ep) == 1L & nchar(op) == 1L & nchar(eo) == 1L & nchar(oo) == 1L
  (eo == ep & oo == op) | (eo == op & oo == ep) |
    (strand_snv & !is_palindromic(ep, op) &
       ((eo == complement_allele(ep) & oo == complement_allele(op)) |
        (eo == complement_allele(op) & oo == complement_allele(ep))))
}

harmonize_alleles <- function(common, ea_p, oa_p, ea_o, oa_o, beta_o) {
  ep <- toupper(common[[ea_p]]); op <- toupper(common[[oa_p]])
  eo <- toupper(common[[ea_o]]); oo <- toupper(common[[oa_o]])

  comp_ep <- complement_allele(ep); comp_op <- complement_allele(op)

  same      <- eo == ep  & oo == op
  flip_beta <- eo == op  & oo == ep
  same_str  <- eo == comp_ep & oo == comp_op
  flip_str  <- eo == comp_op & oo == comp_ep

  palin      <- is_palindromic(ep, op)
  # THE STRAND BRANCHES ARE FOR SNVs ONLY. complement_allele() is chartr(), which complements each
  # base without reversing: "GGA" becomes "CCT", whereas the reverse complement is "TCC". So a
  # genuinely reverse-complemented indel was dropped, while an unrelated pair could match by
  # character-complement alone -- PDZK1 currently accepts pQTL TAAA/T against outcome ATTT/A through
  # same_str on exactly that coincidence. Reversing the string is NOT the fix: a point liftover does
  # not re-anchor or renormalise an indel on the target reference, so no string operation establishes
  # that two multi-base representations describe the same variant. Exact and swapped matches still
  # accept multi-base alleles; only the strand inference is restricted. This can only remove rows,
  # never manufacture them (PDZK1-T2D 447 -> 446 harmonised rows).
  strand_snv <- nchar(ep) == 1L & nchar(op) == 1L & nchar(eo) == 1L & nchar(oo) == 1L
  keep       <- allele_compatible(ep, op, eo, oo)
  needs_flip <- flip_beta | (strand_snv & !palin & flip_str)

  common <- common[keep]
  if (nrow(common) > 0 && any(needs_flip[keep])) {
    common[needs_flip[keep], (beta_o) := -get(beta_o)]
  }
  common
}

# ── rsID 主要合併 + 位置後備合併 ─────────────────────────────────────────────
# 輸入（標準化欄位名，避免衝突）：
#   pqtl_sub: rsid, chrpos_pqtl (可 NA), pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl
#   out_sub:  rsid_out, chrpos_out (可 NA), out_beta, out_varbeta, ea_out, oa_out
# 回傳: data.table with snp_id + needed columns for coloc/harmonize
rsid_pos_merge <- function(pqtl_sub, out_sub) {
  result <- list()
  matched_pqtl_var <- character(0)   # (chrpos:EA:OA) of the pQTL rows the primary merge already used

  # ── Primary: rsID merge ──────────────────────────────────────────────────
  has_p_rsid <- "rsid"     %in% names(pqtl_sub) && any(!is.na(pqtl_sub$rsid))
  has_o_rsid <- "rsid_out" %in% names(out_sub)  && any(!is.na(out_sub$rsid_out))

  if (has_p_rsid && has_o_rsid) {
    # Select minimal columns (no chrpos conflict)
    # THE COORDINATE HAS TO AGREE TOO. Merging on rsID alone accepted an outcome record carrying the
    # same rsID at a DIFFERENT position -- and because the pQTL position was then added to
    # matched_pqtl_pos, the correct record at the right coordinate could no longer be reached by the
    # fallback. One wrong row displaced one right one. Allowing the fallback to run anyway is not the
    # fix: that keeps both, double-counting the pQTL. This is only valid because both sides are now in
    # one build -- GRCh38 outcomes against archive GENPOS, GRCh37 outcomes against projected pQTL
    # positions -- which was not true before this week.
    p_r <- pqtl_sub[!is.na(rsid) & !is.na(chrpos_pqtl),
                     .(rsid, chrpos_pqtl, pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl)]
    o_r <- out_sub[!is.na(rsid_out) & !is.na(chrpos_out),
                    .(rsid_out, chrpos_out, out_beta, out_varbeta, ea_out, oa_out)]
    rsid_m <- merge(p_r, o_r, by.x = c("rsid", "chrpos_pqtl"), by.y = c("rsid_out", "chrpos_out"))
    if (nrow(rsid_m) > 0) {
      rsid_m[, snp_id := rsid]
      # WHAT THE PRIMARY MERGE CONSUMES IS A VARIANT, NOT A COORDINATE. Two things were wrong with
      # excluding by position. (a) An outcome record can carry the pQTL's rsID at the right position
      # with alleles that cannot describe the same variant; harmonisation drops it downstream, but the
      # coordinate had already been closed, so the compatible record at that coordinate was never
      # reached. (b) Even when the rsID hit IS harmonisable, closing its coordinate also removes every
      # OTHER pQTL variant sitting there -- a multiallelic site, or an SNV beside an indel -- although
      # those are different variants the fallback could still match. Keying on (position, alleles)
      # excludes exactly the pQTL row the primary merge already used, and nothing else.
      .ok <- allele_compatible(rsid_m$ea_pqtl, rsid_m$oa_pqtl, rsid_m$ea_out, rsid_m$oa_out)
      matched_pqtl_var <- rsid_m[!is.na(chrpos_pqtl) & !is.na(.ok) & .ok,
                                 paste(chrpos_pqtl, toupper(ea_pqtl), toupper(oa_pqtl), sep = ":")]
      result[["rsid"]] <- rsid_m
      cat(sprintf("    rsID merge: %d SNPs\n", nrow(rsid_m)))
    }
  }

  # ── Fallback: position merge（跳過已 rsID 配對的）──────────────────────
  has_p_pos <- "chrpos_pqtl" %in% names(pqtl_sub) && any(!is.na(pqtl_sub$chrpos_pqtl))
  has_o_pos <- "chrpos_out"  %in% names(out_sub)  && any(!is.na(out_sub$chrpos_out))

  if (has_p_pos && has_o_pos) {
    p_fb <- pqtl_sub[!is.na(chrpos_pqtl) &
                      !(paste(chrpos_pqtl, toupper(ea_pqtl), toupper(oa_pqtl), sep = ":") %in% matched_pqtl_var),
                      .(chrpos_pqtl, pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl)]
    o_fb <- out_sub[!is.na(chrpos_out),
                     .(chrpos_out, out_beta, out_varbeta, ea_out, oa_out)]
    pos_m <- merge(p_fb, o_fb, by.x = "chrpos_pqtl", by.y = "chrpos_out")
    if (nrow(pos_m) > 0) {
      pos_m[, snp_id := chrpos_pqtl]
      result[["pos"]] <- pos_m
      cat(sprintf("    position fallback: +%d SNPs\n", nrow(pos_m)))
    }
  }

  if (length(result) == 0) {
    cat("    rsID merge: 0 SNPs\n")
    return(data.table())
  }
  rbindlist(result, use.names = TRUE, fill = TRUE)
}

# ── 從 UKB-PPP tar 取出 pQTL cis 區域（v5：MAF filter + rsID）──────────────
get_pqtl_region <- function(gene, chrom, region_start, region_end, rsid_map) {
  tar_file <- resolve_ukbppp_tar(gene)   # 2026-07-01 (gate finding 1): canonical assay by instrument OID, not [1]
  if (is.null(tar_file)) return(NULL)

  py <- tempfile(fileext = ".py")
  cat(sprintf('
import tarfile, gzip, pandas as pd, sys

tar_path = r"%s"
chrom = "%s"
region_start = %d
region_end = %d

with tarfile.open(tar_path, "r:") as tf:
    for m in tf.getmembers():
        if f"chr{chrom}_" in m.name and m.name.endswith(".gz"):
            f = tf.extractfile(m)
            with gzip.open(f) as gz:
                df = pd.read_csv(gz, sep=" ")
            df = df[(df["GENPOS"] >= region_start) & (df["GENPOS"] <= region_end)]
            df.to_csv(sys.stdout, index=False)
            break
',
    gsub("\\\\", "/", tar_file), chrom, region_start, region_end),
  file = py)

  # stderr KEPT. It used to be discarded, so a failure reached the log as a bare exit code and
  # nothing else. On 2026-08-02 all 1,380 pairs returned "exit 322" and establishing that python
  # was in fact fine -- the machine had run out of commit and could not create the process --
  # took six hours and a temporary sitecustomize.py hook. Whatever python has to say belongs in
  # the message that reports the failure.
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(c(py, .pyerr))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop(sprintf("region extraction python failed (exit %s) -- NOT a legit empty-region skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: the process produced no diagnostic, which points at process creation rather than the script>"))
  if (length(out) <= 1) return(NULL)

  pqtl_dt <- suppressWarnings(fread(text = paste(out, collapse = "\n")))
  pqtl_dt[, chrpos := paste0(chrom, ":", GENPOS)]

  # MAF filter
  if ("A1FREQ" %in% names(pqtl_dt)) {
    pqtl_dt <- pqtl_dt[!is.na(A1FREQ) & A1FREQ >= MAF_MIN & A1FREQ <= 1 - MAF_MIN]
  }

  # Add rsID from cis_rsid_map (by chrom:pos:alleles key -- see load_rsid_map)
  if (!is.null(rsid_map)) {
    pqtl_dt[, lookup_key := paste(chrpos, toupper(ALLELE0), toupper(ALLELE1), sep = ":")]
    pqtl_dt <- merge(pqtl_dt, rsid_map[, .(lookup_key, rsid)],
                     by = "lookup_key", all.x = TRUE)
  } else {
    pqtl_dt[, rsid := NA_character_]
  }

  pqtl_dt
}

# ── Liftover pQTL GRCh38→GRCh37 ────────────────────────────────────────────
# GRCh38 -> GRCh37, for the two CKDGen outcomes. Their positions are `Pos_b37`, the pQTL positions
# are GRCh38, and comparing the two directly is the same error this whole rebuild exists to fix: it
# does not merely fail to match, it MATCHES THE WRONG VARIANTS. Measured on the delivered eGFR
# analysis, 25 SNPs across 8 genes were two different variants merged into one row (DNER's
# rs10201872 paired with eGFR's rs116361129, 864,716 bp apart). Projecting the pQTL side down to
# GRCh37 makes the comparison meaningful; a variant with no GRCh37 counterpart gets NA and is simply
# not available for the position fallback, which is the correct outcome for it.
liftover_to_hg37 <- function(pqtl_dt) {
  tmp_in  <- tempfile(fileext = ".tsv"); tmp_out <- tempfile(fileext = ".tsv")
  fwrite(pqtl_dt[, .(chrpos,
                     chrom38 = sub(":.*", "", chrpos),
                     pos38   = as.integer(sub(".*:", "", chrpos)))], tmp_in, sep = "\t")
  py <- tempfile(fileext = ".py")
  # `raise`, not `pass`, in the handler below. A chain GAP already produces an empty result and
  # therefore no output row, which becomes chrpos37 = NA in the all.x merge -- the correct answer for
  # an unmappable variant. `pass` gave a MALFORMED row that same answer, and since the subprocess
  # still exited 0 the fail-loud check below never saw it. The two are now distinguishable: no
  # mapping is NA, a malformed row is a non-zero exit and a stopped run.
  cat(sprintf('
from pyliftover import LiftOver
import csv

lo = LiftOver(r"%s")
with open(r"%s") as fin, open(r"%s", "w", newline="") as fout:
    rdr = csv.DictReader(fin, delimiter="\\t")
    wtr = csv.writer(fout, delimiter="\\t")
    wtr.writerow(["chrpos38", "chrpos37"])
    for row in rdr:
        try:
            chrom = "chr" + row["chrom38"].lstrip("chr")
            pos38 = int(row["pos38"])
            res = lo.convert_coordinate(chrom, pos38 - 1)
            if res and len(res) == 1:
                wtr.writerow([row["chrpos"], "%%s:%%d" %% (res[0][0].lstrip("chr"), res[0][1] + 1)])
        except Exception:
            raise
',
    gsub("\\\\", "/", LIFTOVER_CHAIN_38_TO_37),
    gsub("\\\\", "/", tmp_in), gsub("\\\\", "/", tmp_out)), file = py)
  # stderr KEPT -- same reason as the extraction sites above.
  .pyerr <- tempfile(fileext = ".stderr")
  .st <- system2(PYTHON_BIN, args = py, stdout = FALSE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(c(py, tmp_in, .pyerr))
  if (!identical(as.integer(.st), 0L))
    stop(sprintf("coloc pQTL liftover (hg38->hg37) FAILED (exit %s) -- fail loud, not silent SNP loss. python stderr: %s", .st, if (nzchar(.emsg)) .emsg else "<empty>"))
  if (!file.exists(tmp_out) || file.info(tmp_out)$size == 0) { unlink(tmp_out); return(NULL) }
  lmap <- fread(tmp_out, col.names = c("chrpos38", "chrpos37")); unlink(tmp_out)
  # A GRCh38 position that lifts to more than one GRCh37 place is ambiguous; the Python above already
  # emits only unique mappings, and this keeps the join 1:1 on what survives.
  lmap <- lmap[, .SD[1], by = chrpos38]
  merge(pqtl_dt, lmap, by.x = "chrpos", by.y = "chrpos38", all.x = TRUE)
}


# ── Ensembl GRCh38 gene 座標 ─────────────────────────────────────────────────
# Refuses to let colocalisation run against a coordinate table whose provenance cannot be shown.
# 03_gene_coords.R writes to this same path and does not remove the sidecar, so "a certificate exists"
# certifies nothing on its own -- the hash and the builder name both have to match.
.assert_coords_certified <- function() {
  cert <- paste0(GENE_HG38_COORDS, ".manifest.tsv")
  bad  <- NULL
  # ONE read. Hashing the path and then re-reading the path lets a concurrent rebuild slip
  # unvalidated bytes into the cache between the two, so the bytes that are checked are the bytes
  # that are parsed.
  # Checked HERE too, not only by get_gene_hg38_coords(). is_mhc_gene() calls this function
  # directly, and without this a missing table on that path raised a bare file-open error instead of
  # the message naming the script that rebuilds it.
  if (!file.exists(GENE_HG38_COORDS))
    stop("missing ", basename(GENE_HG38_COORDS), " -- run ",
         "scripts/paper3.1_ckm/cascade/03b_gene_coords_hg38.R (03_gene_coords.R is the retired ",
         "liftover version and must not be used to rebuild it)")
  bytes <- readBin(GENE_HG38_COORDS, "raw", file.size(GENE_HG38_COORDS))
  have  <- digest::digest(bytes, algo = "md5", serialize = FALSE)
  if (!file.exists(cert)) {
    bad <- paste0("no certificate (", basename(cert), ") beside it")
  } else {
    cm <- tryCatch(data.table::fread(cert), error = function(e) NULL)
    # REQUIRED, not checked-if-present. Each of the four build/window/manifest fields below was
    # guarded by `"x" %in% names(cm) &&`, so a certificate that simply omitted them passed the very
    # check whose error message says the table's build was verified.
    if (is.null(cm) || !all(c("file", "md5", "builder", "phase1_gene_coord_build",
                              "phase1_chrpos_build", "cis_window_kb",
                              "phase1_manifest_md5") %in% names(cm)) || nrow(cm) != 1L)
      bad <- "its certificate is unreadable or malformed"
    # `file` and `rows` were required to EXIST and then never compared, so a certificate could
    # name a different table and state a different row count while still being accepted. The md5
    # prevents that from changing any coordinate, but a certificate that misstates what it
    # describes is not provenance.
    # ...and `rows`, which the comment above says was compared and which the code never looked at.
    # Counted from the bytes already in hand, not parsed: this function runs on EVERY
    # get_gene_hg38_coords() call, and a full fread per call to check one integer is a cost the
    # md5 above does not pay. The table has no embedded newlines, so lines-minus-header is the
    # row count exactly.
    # Digits only, then compared. as.integer() alone turns "2920.9", "+2920" and "2.920e3" into
    # 2920, so a certificate that does not state the row count is accepted as though it did.
    else if (!("rows" %in% names(cm)) ||
             !grepl("^[[:space:]]*[0-9]+[[:space:]]*$", as.character(cm$rows[1])) ||
             !identical(as.integer(cm$rows[1]), sum(bytes == as.raw(10L)) - 1L))
      bad <- "its certificate states a row count that the table does not have"
    else if (!identical(cm$file[1], basename(GENE_HG38_COORDS)))
      bad <- paste0("its certificate describes '", cm$file[1], "', not ",
                    basename(GENE_HG38_COORDS))
    else if (!identical(tolower(have), tolower(cm$md5[1])))
      bad <- paste0("its certificate describes a DIFFERENT file (md5 mismatch); the table was most ",
                    "likely rebuilt by the retired 03_gene_coords.R afterwards")
    else if (!identical(cm$builder[1], "03b_gene_coords_hg38.R"))
      bad <- paste0("its certificate names builder '", cm$builder[1], "', not 03b_gene_coords_hg38.R")
    else if (!"builder_md5" %in% names(cm) || is.na(cm$builder_md5[1]) || !nzchar(cm$builder_md5[1]))
      bad <- paste0("its certificate does not record WHICH VERSION of 03b_gene_coords_hg38.R wrote ",
                    "it, so an output from any earlier 03b would be accepted here")
    else {
      # The builder's VERSION, symmetrically with the gate 03b applies to ITS upstream builder. A
      # coordinate table carries whatever policy the script that wrote it had -- which genes get an
      # alias row, which are withheld -- so accepting a table from an unknown 03b means accepting an
      # unknown policy.
      .b3 <- file.path(BASE_DIR, "scripts/paper3.1_ckm/cascade/03b_gene_coords_hg38.R")
      .b3m <- if (file.exists(.b3)) unname(tools::md5sum(.b3)) else NA_character_
      if (is.na(.b3m))
        bad <- paste0("cannot hash ", basename(.b3), " to check the table against the builder ",
                      "currently on disk")
      else if (!identical(tolower(cm$builder_md5[1]), tolower(.b3m)))
        bad <- paste0("it was built by a DIFFERENT version of 03b_gene_coords_hg38.R (certificate ",
                      substr(cm$builder_md5[1], 1, 12), ", current ", substr(.b3m, 1, 12),
                      "); its alias and withholding policy may predate the current rules")
      else {
        # And the UPSTREAM table it was built FROM. 03b records source/source_md5 for exactly this,
        # and nothing was checking them: rebuild gene_coords_hg38_ukbppp.tsv, then fail to rerun 03b,
        # and the older installed coordinates stay fully certified. Every gene whose GENCODE
        # boundaries moved would then define a stale coloc window, silently -- which is the failure
        # this certificate exists to make impossible.
        .src <- if ("source" %in% names(cm) && nzchar(cm$source[1]))
                  file.path(BASE_DIR, "data/annotation", basename(cm$source[1])) else NA_character_
        if (!"source_md5" %in% names(cm) || is.na(cm$source_md5[1]) || !nzchar(cm$source_md5[1]))
          bad <- "its certificate does not record which upstream coordinate table it was built from"
        else if (is.na(.src) || !file.exists(.src))
          bad <- paste0("the upstream table it names (", cm$source[1], ") is not on disk, so it ",
                        "cannot be shown to be current")
        else if (!identical(tolower(unname(tools::md5sum(.src))), tolower(cm$source_md5[1])))
          bad <- paste0("the upstream table ", cm$source[1], " has been REBUILT since this one was ",
                        "derived from it (certificate ", substr(cm$source_md5[1], 1, 12), ", current ",
                        substr(unname(tools::md5sum(.src)), 1, 12), "). Rerun ",
                        "cascade/03b_gene_coords_hg38.R")
        # AND THE PHASE-1 SET THE ATTESTATION WAS MADE AGAINST. 03b installs this table only after
        # checking that every instrumented gene is covered by it AND that the boundaries agree with
        # the ones the instruments were selected in; it records instruments_md5 for exactly that
        # reason, and nothing was reading the field. Replace Phase 1 without rerunning 03b and the
        # table stays fully certified while its attestation describes the PREVIOUS generation: a
        # newly instrumented symbol then fails its lookup at coloc time, and a gene whose
        # boundaries differ makes the instrument window and the coloc region two different
        # intervals -- the GRCh37-against-GRCh38 defect this rebuild exists to remove.
        else if (!"instruments_md5" %in% names(cm) || is.na(cm$instruments_md5[1]) ||
                 !nzchar(cm$instruments_md5[1]))
          bad <- paste0("its certificate does not record which Phase-1 instrument set its coverage ",
                        "and boundary check was made against, so a Phase-1 rebuild would go unnoticed")
        # Not skipped when the path is unknown or absent: this is the same check, and a check that
        # cannot be made is not a check that passed.
        else if (!exists("PHASE1_INSTRUMENTS") || !file.exists(PHASE1_INSTRUMENTS))
          bad <- paste0("the Phase-1 instrument set its coverage check was made against is not ",
                        "available here, so the table cannot be shown to describe the active ",
                        "instrument generation")
        # THE GENERATION THIS SESSION IS RUNNING, not whatever is at the path now. INSTRUMENTS_MD5 is
        # pinned when 00_config.R is sourced, and every threshold, the gene universe and the outcome
        # map are derived from it; comparing the certificate with a freshly hashed path accepts a
        # table certified against a generation this session is not using.
        # THE BUILD AND THE WINDOW THE TABLE WAS CERTIFIED AGAINST. 03b records them; nothing read
        # them, so a Phase-1 manifest relabelled to GRCh37 afterwards, or a CIS_WINDOW_KB changed in
        # config, left every hash intact and the table fully accepted.
        else if (!all(c(cm$phase1_gene_coord_build[1], cm$phase1_chrpos_build[1]) == "GRCh38"))
          bad <- "its certificate records a Phase-1 build other than GRCh38"
        # Digits only, for the same reason as `rows` above: as.integer("500.9") is 500, so a
        # certificate stating a window this run does not use was accepted as though it matched.
        else if (exists("CIS_WINDOW_KB") &&
                 (!grepl("^[[:space:]]*[0-9]+[[:space:]]*$", as.character(cm$cis_window_kb[1])) ||
                  !identical(as.integer(cm$cis_window_kb[1]), as.integer(CIS_WINDOW_KB))))
          bad <- paste0("it was certified for a ", cm$cis_window_kb[1], " kb cis window and this run ",
                        "uses ", CIS_WINDOW_KB, " kb")
        else if (!identical(tolower(cm$phase1_manifest_md5[1]),
                            tolower(unname(tools::md5sum(paste0(PHASE1_INSTRUMENTS, ".manifest.tsv"))))))
          bad <- "the Phase-1 manifest has changed since this table was certified against it"
        else if (!identical(tolower(if (exists("INSTRUMENTS_MD5")) INSTRUMENTS_MD5
                                    else unname(tools::md5sum(PHASE1_INSTRUMENTS))),
                            tolower(cm$instruments_md5[1])))
          bad <- paste0("the Phase-1 instrument set has been REBUILT since this table's coverage and ",
                        "boundary check (certificate ", substr(cm$instruments_md5[1], 1, 12),
                        ", current ", substr(unname(tools::md5sum(PHASE1_INSTRUMENTS)), 1, 12),
                        "). Rerun cascade/03b_gene_coords_hg38.R")
        else {
          # The whole upstream chain, not just the table 03b read. 03b refuses to run unless the
          # coordinate builder AND the GENCODE annotation match what is on disk, but those can both
          # change afterwards without changing the coordinate TSV's bytes -- and then this table stays
          # accepted under rules that no longer apply. It matters most for the alias rows, which are
          # decided FROM the GTF: a different release can make an existing alias the current symbol of
          # some other locus while every lookup here goes on succeeding.
          # Hashed FROM DISK. A previous revision compared the certificate's copy of these hashes
          # against the upstream MANIFEST's copy -- both written from the same source at the same
          # moment -- so the comparison could not fail and proved nothing. What has to be checked is
          # the builder script and the annotation as they are NOW.
          # Keyed on the table AND the certificate AND the cheap upstream hash -- not the table
          # alone, which let a changed builder or a replaced certificate reuse an earlier pass. The
          # GTF is the one input whose hash is not recomputed per call: 76 MB times a thousand gene
          # lookups is 76 GB, so it is keyed on its size and mtime instead, and a same-size
          # replacement with a preserved timestamp is the one substitution this cannot see.
          # Size and mtime in the KEY, hashes only on a miss. Hashing 39 MB to decide whether the
          # cache is valid, and only then consulting it, meant ~39 GB of reading across a coloc run
          # -- about 100 seconds -- to answer a question the cache already knew. Same trade as the
          # GTF, and the same stated limit: a same-size replacement with a preserved timestamp is
          # not seen.
          .oth <- list(idmapping_md5 = file.path(BASE_DIR, "data/annotation",
                                                 "HUMAN_9606_idmapping.dat.gz"),
                       release_listing_md5 = file.path(BASE_DIR, "data/ukb_ppp",
                                                       "eur_folder_contents.json"))
          .oi <- lapply(.oth, function(f) {
            i <- file.info(f)
            paste(i$size, format(as.numeric(i$mtime), digits = 22), sep = "|")
          })
          .b93k <- file.path(BASE_DIR, "scripts/_shared/localize/93_build_gene_coords_hgnc.R")
          .gtfk <- file.path(BASE_DIR, "data/annotation/gencode.v50.basic.annotation.gtf.gz")
          .gi <- file.info(.gtfk)
          .key <- paste(have,
                        if (file.exists(cert)) unname(tools::md5sum(cert)) else "nocert",
                        if (file.exists(.b93k)) unname(tools::md5sum(.b93k)) else "no93",
                        # The full stored value, not a rounded one. format() gives whole seconds and
                        # %.6f still collapses distinct sub-microsecond timestamps, so each rounding
                        # made the accepted exception wider than the sentence above admits.
                        .gi$size, format(as.numeric(.gi$mtime), digits = 22),
                        # The two inputs checked on a cache MISS have to be in the key, or the first
                        # successful lookup makes every later change to them invisible. Both are
                        # small enough to hash outright.
                        .oi$idmapping_md5, .oi$release_listing_md5, sep = "|")
          if (identical(get0(".GENE_HG38_PROV_OK", envir = .GlobalEnv), .key)) {
            NULL          # this exact combination already passed the upstream chain this session
          } else if (!all(c("upstream_builder_md5", "gtf_md5") %in% names(cm)))
            bad <- paste0("its certificate does not record the upstream builder and annotation it ",
                          "was derived under, so a change to either would go unnoticed")
          else {
            .b93 <- file.path(BASE_DIR, "scripts/_shared/localize/93_build_gene_coords_hgnc.R")
            .gtf <- file.path(BASE_DIR, "data/annotation/gencode.v50.basic.annotation.gtf.gz")
            # All four inputs 03b validates, not two of them. The UniProt mapping decides which gene
            # each assay resolves to and the release listing decides the assay universe; either can
            # change without touching the coordinate table's bytes, and the coordinates would then be
            # derived under a mapping that no longer exists while every hash here still matched.
            .oth <- list(idmapping_md5 = file.path(BASE_DIR, "data/annotation",
                                                   "HUMAN_9606_idmapping.dat.gz"),
                         release_listing_md5 = file.path(BASE_DIR, "data/ukb_ppp",
                                                         "eur_folder_contents.json"))
            .ub  <- if (file.exists(.b93)) unname(tools::md5sum(.b93)) else NA_character_
            .ug  <- if (file.exists(.gtf)) unname(tools::md5sum(.gtf)) else NA_character_
            if (is.na(.ub) || is.na(.ug))
              bad <- paste0("cannot hash ", if (is.na(.ub)) basename(.b93) else basename(.gtf),
                            ", so this table cannot be checked against the inputs it came from")
            else if (!identical(tolower(cm$upstream_builder_md5[1]), tolower(.ub)))
              bad <- paste0("the coordinate BUILDER on disk differs from the one this table was ",
                            "derived under (certificate ", substr(cm$upstream_builder_md5[1], 1, 12),
                            ", current ", substr(.ub, 1, 12), "). Rerun ",
                            "scripts/_shared/localize/93_build_gene_coords_hgnc.R, then ",
                            "cascade/03b_gene_coords_hg38.R")
            else if (!identical(tolower(cm$gtf_md5[1]), tolower(.ug)))
              bad <- paste0("the GENCODE annotation on disk differs from the one this table was ",
                            "derived under (certificate ", substr(cm$gtf_md5[1], 1, 12), ", current ",
                            substr(.ug, 1, 12), "); its alias rows may now point at other genes. ",
                            "Rebuild both, then rerun cascade/03b_gene_coords_hg38.R")
            else for (.k in names(.oth)) {
              if (!.k %in% names(cm) || is.na(cm[[.k]][1]) || !nzchar(cm[[.k]][1])) {
                bad <- paste0("its certificate does not record ", basename(.oth[[.k]]),
                              ", one of the inputs the coordinates were derived from")
                break
              }
              .cur2 <- if (file.exists(.oth[[.k]])) unname(tools::md5sum(.oth[[.k]])) else NA_character_
              if (is.na(.cur2) || !identical(tolower(cm[[.k]][1]), tolower(.cur2))) {
                bad <- paste0(basename(.oth[[.k]]), " differs from the copy this table was derived ",
                              "under (certificate ", substr(cm[[.k]][1], 1, 12), ", current ",
                              substr(.cur2, 1, 12), "). Rerun ",
                              "scripts/_shared/localize/93_build_gene_coords_hgnc.R, then ",
                              "cascade/03b_gene_coords_hg38.R")
                break
              }
            }
          }
        }
      }
    }
  }
  if (is.null(bad)) {
    # The expensive half of this check -- hashing the 80 MB GTF and both builder scripts -- is
    # memoised against the coordinate table's own hash. It runs on every gene lookup, and roughly a
    # thousand lookups per coloc run meant ~80 GB of GTF reads for an answer that cannot change while
    # the table does not. The CHEAP half (hashing this small table and comparing it to its
    # certificate) still runs every time, so a table replaced mid-session is still caught.
    if (exists(".key", inherits = FALSE)) assign(".GENE_HG38_PROV_OK", .key, envir = .GlobalEnv)
    return(list(md5 = have, bytes = bytes))
  }
  msg <- paste0(basename(GENE_HG38_COORDS), ": ", bad, ". It may be the retired liftover table from ",
                "03_gene_coords.R rather than the GRCh38 table from ",
                "scripts/paper3.1_ckm/cascade/03b_gene_coords_hg38.R. Rebuild it, or set ",
                "CKM_ALLOW_UNCERTIFIED_COORDS=1 to proceed deliberately (historical reruns).")
  # FATAL by default -- this is deliberate, and the header above says "stops loudly" for that reason.
  # A warning is useless here: cascade/04_coloc.R redirects each parallel worker's
  # stderr to a file and reports success from exit codes alone, so nothing would reach anyone.
  # Compared to "1", not merely non-empty, so CKM_ALLOW_UNCERTIFIED_COORDS=0 cannot switch it off.
  if (!identical(Sys.getenv("CKM_ALLOW_UNCERTIFIED_COORDS"), "1")) stop(msg)
  # Once per DISTINCT failure -- not once per gene, and not once per session. Per gene prints the
  # same line hundreds of times in a loop until nobody reads it; per session means a table replaced
  # mid-run, or a second and different failure, passes in silence after the first warning.
  seen <- get0(".GENE_HG38_WARNED", envir = .GlobalEnv, ifnotfound = character())
  key  <- paste(have, bad)
  if (!key %in% seen) {
    warning(msg, call. = FALSE, immediate. = TRUE)
    assign(".GENE_HG38_WARNED", c(seen, key), envir = .GlobalEnv)
  }
  list(md5 = have, bytes = bytes)
}

get_gene_hg38_coords <- function(gene_name) {
  # LOCAL, no API: read the vendored gene -> GRCh38 coords table
  # (data/ld_reference/gene_hg38_coords.tsv, built by cascade/03b_gene_coords_hg38.R from
  # the certified GRCh38 coordinate table. localize/26 and 03_gene_coords.R are the retired liftover
  # route -- and 03_gene_coords.R STILL WRITES TO THIS PATH if anyone runs it, which is the whole
  # reason the certificate below exists: the file being present proves nothing about who made it).
  # Replaces the former Ensembl REST call whose transient failures silently dropped
  # 41 coloc pairs incl the PCSK9->CAD positive control. No network; NO Ensembl fallback
  # by design -- a missing gene STOPS (it was once a silent skip, which dropped 41 pairs), and an
  # uncertified or missing TABLE stops too.
  # Validated on EVERY call, not only when the cache is cold. A long-lived session that loaded the
  # table before this check existed would otherwise keep using it forever, which is exactly the
  # situation the check is for. The work is one md5 of a small TSV; the table itself stays cached.
  if (!file.exists(GENE_HG38_COORDS))
    stop("missing local gene_hg38_coords.tsv -- run scripts/paper3.1_ckm/cascade/03b_gene_coords_hg38.R (03_gene_coords.R ",
         "is the retired liftover version and must not be used to rebuild it)")
  cur <- .assert_coords_certified()
  # Set on the FIRST certified read, by whoever gets there -- not only by main(). A script that
  # sources this file and loops over genes was previously unpinned, so a table replaced between two
  # iterations was accepted and the results silently spanned two versions.
  .pin <- get0(".COLOC_RUN_COORD_MD5", envir = .GlobalEnv)
  if (is.null(.pin)) {
    assign(".COLOC_RUN_COORD_MD5", cur$md5, envir = .GlobalEnv)
    .pin <- cur$md5
  }
  if (!identical(.pin, cur$md5))
    stop("the coordinate table changed while this run was in progress (pinned ",
         substr(.pin, 1, 12), ", now ", substr(cur$md5, 1, 12), "). Pairs completed before the ",
         "change used the old cis windows and pairs after it would use new ones, in one result set. ",
         "Nothing further was computed; delete the checkpoint and start again.")
  # Keyed to the validated hash, not merely to "a table has been loaded". Otherwise a process that
  # cached the table before a rebuild keeps serving the OLD coordinates while the certificate check
  # passes on every call -- validating the file on disk proves nothing about the object in memory.
  if (!exists(".GENE_HG38_TBL", envir = .GlobalEnv) ||
      !identical(get0(".GENE_HG38_MD5", envir = .GlobalEnv), cur$md5)) {
    t <- fread(text = rawToChar(cur$bytes)); setkey(t, gene)
    assign(".GENE_HG38_TBL", t,       envir = .GlobalEnv)
    assign(".GENE_HG38_MD5", cur$md5, envir = .GlobalEnv)
  }
  r <- get(".GENE_HG38_TBL", envir = .GlobalEnv)[gene == gene_name]
  if (nrow(r) == 0) stop(sprintf("gene not in local hg38 coords table: %s -- a suggestive hg38 coloc pair would be SILENTLY dropped; add it via cascade/03b_gene_coords_hg38.R or exclude it explicitly", gene_name))  # 2026-07-01 gate: fatal, not silent skip
  list(chrom = as.character(r$chrom[1]), start = as.integer(r$start[1]), end = as.integer(r$end[1]))
}

# ── FinnGen（hg38）local 過濾 + rsids 欄位 ─────────────────────────────────
get_finngen_region <- function(outcome_id, gene_name) {
  meta       <- GWAS_META[[outcome_id]]
  local_file <- file.path(OUTCOME_DIR, meta$file)
  if (!file.exists(local_file))   # 2026-07-01 Codex-phase3 #2: missing outcome FILE is systematic (drops the whole outcome), NOT a legit per-pair skip
    stop(sprintf("03 coloc: outcome GWAS file MISSING: %s (outcome %s) -- data-foundation failure; every configured outcome must have its file", basename(local_file), outcome_id))

  hg38 <- get_gene_hg38_coords(gene_name)
  if (is.null(hg38)) {
    cat("    [error] 無法取得 GRCh38 座標:", gene_name, "\n"); return(NULL)
  }
  q_chrom <- hg38$chrom
  q_start <- max(1L, hg38$start - COLOC_REGION_KB * 1000L)
  q_end   <- hg38$end   + COLOC_REGION_KB * 1000L

  py <- tempfile(fileext = ".py")
  cat(sprintf('
import gzip, sys
with gzip.open(r"%s", "rt") as f:
    hdr = next(f)
    sys.stdout.write(hdr.lstrip("#"))
    for line in f:
        p = line.split("\\t")
        if p[0] == "%s":
            try:
                if %d <= int(p[1]) <= %d:
                    sys.stdout.write(line)
            except ValueError:
                raise
',
    gsub("\\\\", "/", local_file), q_chrom, q_start, q_end),
  file = py)

  # stderr KEPT. It used to be discarded, so a failure reached the log as a bare exit code and
  # nothing else. On 2026-08-02 all 1,380 pairs returned "exit 322" and establishing that python
  # was in fact fine -- the machine had run out of commit and could not create the process --
  # took six hours and a temporary sitecustomize.py hook. Whatever python has to say belongs in
  # the message that reports the failure.
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(c(py, .pyerr))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop(sprintf("region extraction python failed (exit %s) -- NOT a legit empty-region skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: the process produced no diagnostic, which points at process creation rather than the script>"))
  if (length(out) <= 1) { cat("    [空] FinnGen 區域無資料\n"); return(NULL) }

  dt <- tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("coloc outcome-region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e))))  # 2026-07-01 (closure sweep)
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("coloc outcome-region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")
  # fread stops at the first malformed line with only a warning, which used to be suppressed here: CKDGen
  # regions were cut short silently. Every extracted row must come back.
  if (nrow(dt) != length(out) - 1L) fail_loud_data(sprintf("coloc outcome-region fread kept %d of %d extracted rows -- the parse stopped early, so the region is truncated", nrow(dt), length(out) - 1L))

  setnames(dt, names(dt), gsub("^#", "", names(dt)))

  # rsID from rsids column (take first if multiple comma-separated)
  if ("rsids" %in% names(dt)) {
    dt[, rsid_out := sub(",.*", "", rsids)]
    dt[is.na(rsid_out) | rsid_out == "" | rsid_out == ".", rsid_out := NA_character_]
  } else {
    dt[, rsid_out := NA_character_]
  }

  dt[, chrpos_out := paste0(chrom, ":", pos)]
  dt
}

# ── Nielsen AF（hg37）local 過濾 + rs_dbSNP147 ─────────────────────────────
get_nielsen_region <- function(outcome_id, chrom, region_start, region_end) {
  meta       <- GWAS_META[[outcome_id]]
  local_file <- file.path(OUTCOME_DIR, meta$file)
  if (!file.exists(local_file))   # 2026-07-01 Codex-phase3 #2: missing outcome FILE is systematic (drops the whole outcome), NOT a legit per-pair skip
    stop(sprintf("03 coloc: outcome GWAS file MISSING: %s (outcome %s) -- data-foundation failure; every configured outcome must have its file", basename(local_file), outcome_id))

  py <- tempfile(fileext = ".py")
  py_sep <- if (identical(meta$sep, "\t")) "\\t" else meta$sep   # 2026-06-30 generalized delimiter (Wuttke CKD = space)
  cat(sprintf('
import gzip, csv, sys
chrom_target = "%s"
start = %d
end   = %d
with gzip.open(r"%s", "rt") as f:
    # space-delimited files (Wuttke CKDGen) pad some fields with extra spaces: csv.reader turns
    # each extra space into an empty field and shifts the row, so split on runs of whitespace instead
    rdr = (l.split() for l in f) if "%s" == " " else csv.reader(f, delimiter="%s")
    hdr = next(rdr)
    try:
        ci = hdr.index("%s")
        pi = hdr.index("%s")
    except ValueError as e:
        sys.exit(1)
    print("\\t".join(hdr))
    for row in rdr:
        try:
            chr_val = row[ci].strip().lstrip("chr")
            if chr_val == chrom_target and start <= int(row[pi]) <= end:
                print("\\t".join(row))
        except (ValueError, IndexError):
            raise
',
    as.character(chrom), region_start, region_end,
    gsub("\\\\", "/", local_file), py_sep, py_sep,
    meta$chrom_col, meta$pos_col),
  file = py)

  # stderr KEPT. It used to be discarded, so a failure reached the log as a bare exit code and
  # nothing else. On 2026-08-02 all 1,380 pairs returned "exit 322" and establishing that python
  # was in fact fine -- the machine had run out of commit and could not create the process --
  # took six hours and a temporary sitecustomize.py hook. Whatever python has to say belongs in
  # the message that reports the failure.
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(c(py, .pyerr))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop(sprintf("region extraction python failed (exit %s) -- NOT a legit empty-region skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: the process produced no diagnostic, which points at process creation rather than the script>"))
  if (length(out) <= 1) { cat(sprintf("    [空] %s 區域無資料\n", outcome_id)); return(NULL) }

  dt <- tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("coloc outcome-region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e))))  # 2026-07-01 (closure sweep)
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("coloc outcome-region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")
  # fread stops at the first malformed line with only a warning, which used to be suppressed here: CKDGen
  # regions were cut short silently. Every extracted row must come back.
  if (nrow(dt) != length(out) - 1L) fail_loud_data(sprintf("coloc outcome-region fread kept %d of %d extracted rows -- the parse stopped early, so the region is truncated", nrow(dt), length(out) - 1L))

  dt[, (meta$chrom_col) := gsub("^chr", "", get(meta$chrom_col))]

  # rsID from meta$rsid_col (generalized 2026-06-30: Nielsen rs_dbSNP147 / Wuttke RSID / FinnGen rsids)
  rcol <- meta$rsid_col
  if (!is.null(rcol) && rcol %in% names(dt)) {
    dt[, rsid_out := as.character(get(rcol))]
    dt[is.na(rsid_out) | !grepl("^rs", rsid_out), rsid_out := NA_character_]
  } else {
    dt[, rsid_out := NA_character_]
  }

  dt[, chrpos_out := paste0(get(meta$chrom_col), ":", get(meta$pos_col))]
  dt
}

# ── GWAS Catalog harmonized hg38（MEGASTROKE / CKDGen）region 過濾 ─────────
get_harmonized_hg38_region <- function(outcome_id, gene_name) {
  meta       <- GWAS_META[[outcome_id]]
  local_file <- file.path(OUTCOME_DIR, meta$file)
  if (!file.exists(local_file))   # 2026-07-01 Codex-phase3 #2: missing outcome FILE is systematic (drops the whole outcome), NOT a legit per-pair skip
    stop(sprintf("03 coloc: outcome GWAS file MISSING: %s (outcome %s) -- data-foundation failure; every configured outcome must have its file", basename(local_file), outcome_id))

  hg38 <- get_gene_hg38_coords(gene_name)
  if (is.null(hg38)) {
    cat("    [error] 無法取得 GRCh38 座標:", gene_name, "\n"); return(NULL)
  }
  q_chrom <- hg38$chrom
  q_start <- max(1L, hg38$start - COLOC_REGION_KB * 1000L)
  q_end   <- hg38$end   + COLOC_REGION_KB * 1000L

  chrom_col <- meta$chrom_col   # "hm_chrom"
  pos_col   <- meta$pos_col     # "hm_pos"

  py <- tempfile(fileext = ".py")
  cat(sprintf('
import gzip, csv, sys
target_chrom = "%s"
start = %d
end   = %d
with gzip.open(r"%s", "rt") as f:
    rdr = csv.reader(f, delimiter="\\t")
    hdr = next(rdr)
    try:
        ci = hdr.index("%s")
        pi = hdr.index("%s")
    except ValueError as e:
        sys.exit(1)
    print("\\t".join(hdr))
    for row in rdr:
        try:
            chr_val = row[ci].lstrip("chr")
            if chr_val == target_chrom and start <= int(row[pi]) <= end:
                print("\\t".join(row))
        except (ValueError, IndexError):
            raise
',
    q_chrom, q_start, q_end,
    gsub("\\\\", "/", local_file),
    chrom_col, pos_col),
  file = py)

  # stderr KEPT. It used to be discarded, so a failure reached the log as a bare exit code and
  # nothing else. On 2026-08-02 all 1,380 pairs returned "exit 322" and establishing that python
  # was in fact fine -- the machine had run out of commit and could not create the process --
  # took six hours and a temporary sitecustomize.py hook. Whatever python has to say belongs in
  # the message that reports the failure.
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(c(py, .pyerr))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop(sprintf("region extraction python failed (exit %s) -- NOT a legit empty-region skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: the process produced no diagnostic, which points at process creation rather than the script>"))
  if (length(out) <= 1) { cat("    [空] harmonized 區域無資料\n"); return(NULL) }

  dt <- tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("coloc outcome-region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e))))  # 2026-07-01 (closure sweep)
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("coloc outcome-region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")
  # fread stops at the first malformed line with only a warning, which used to be suppressed here: CKDGen
  # regions were cut short silently. Every extracted row must come back.
  if (nrow(dt) != length(out) - 1L) fail_loud_data(sprintf("coloc outcome-region fread kept %d of %d extracted rows -- the parse stopped early, so the region is truncated", nrow(dt), length(out) - 1L))

  dt[, (chrom_col) := gsub("^chr", "", get(chrom_col))]

  # rsID from hm_rsid (or meta$rsid_col)
  rsid_col <- meta$rsid_col
  if (rsid_col %in% names(dt)) {
    dt[, rsid_out := get(rsid_col)]
    dt[is.na(rsid_out) | !grepl("^rs", rsid_out), rsid_out := NA_character_]
  } else {
    dt[, rsid_out := NA_character_]
  }

  dt[, chrpos_out := paste0(get(chrom_col), ":", get(pos_col))]
  dt
}

# ── coloc.abf wrapper ────────────────────────────────────────────────────────
run_coloc_abf <- function(pqtl_beta, pqtl_varbeta, n_pqtl, out_beta, out_varbeta, meta, snp_ids) {
  if (length(pqtl_beta) < 20) { cat("    [skip] SNP < 20\n"); return(NULL) }

  ds1 <- list(beta = pqtl_beta, varbeta = pqtl_varbeta,
               N = n_pqtl, sdY = 1, type = "quant", snp = snp_ids)

  ds2 <- list(beta = out_beta, varbeta = out_varbeta,
               N = meta$n, type = meta$type, snp = snp_ids)
  # sdY comes from the outcome's own entry, not a literal 1. If a quantitative outcome is added
  # without one, ds2$sdY is NULL and coloc refuses the dataset -- which is the right failure, rather
  # than silently assuming a standardised phenotype.
  if (meta$type == "cc") ds2$s <- meta$ncase / meta$n else ds2$sdY <- meta$sdY

  tryCatch(
    coloc::coloc.abf(dataset1 = ds1, dataset2 = ds2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),  # explicit standard priors (== coloc defaults; de-ambiguate per Plan)
    # 2026-07-01 (gate): a coloc.abf FAILURE is not a legit skip — propagate so the per-pair loop
    # records status "error" (-> end-of-loop abort), instead of silently degrading to NULL -> "skip".
    error = function(e) stop(sprintf("coloc.abf failed (data problem, NOT a legit skip): %s", conditionMessage(e)))
  )
}

# ── 主 coloc 函式（v5）──────────────────────────────────────────────────────
run_coloc_pair <- function(gene, outcome_id, chrom, gene_start, gene_end, rsid_map) {
  region_start <- max(1L, gene_start - COLOC_REGION_KB * 1000L)
  region_end   <- gene_end   + COLOC_REGION_KB * 1000L

  cat(sprintf("  coloc: %s ~ %s | chr%s:%d-%d\n",
              gene, outcome_id, chrom, region_start, region_end))
  flush(stdout())

  meta   <- GWAS_META[[outcome_id]]
  genome <- meta$genome

  if (genome == "rsid_only") {
    cat("    [skip] rsid_only GWAS（無完整本地檔）\n")
    return(NULL)
  }

  # ── pQTL 資料（含 MAF filter + rsID）────────────────────────────────────
  pqtl_dt <- get_pqtl_region(gene, chrom, region_start, region_end, rsid_map)
  # 2026-07-01 (gate #13 finding 2): a resolved tar returning NO pQTL data, or missing required columns,
  # is a DATA-FOUNDATION failure (not a legit low-SNP skip) -> propagate as error, don't silently skip.
  if (is.null(pqtl_dt) || nrow(pqtl_dt) == 0)
    stop(sprintf("coloc: gene %s pQTL extraction returned NO data (tar resolved but empty) -- data-foundation failure", gene))

  req_cols <- c("ALLELE0", "ALLELE1", "BETA", "SE", "INFO", "N")
  missing  <- setdiff(req_cols, names(pqtl_dt))
  if (length(missing) > 0)
    stop(sprintf("coloc: gene %s pQTL missing required columns: %s -- data-format failure", gene, paste(missing, collapse = ", ")))

  pqtl_dt <- pqtl_dt[!is.na(INFO) & INFO >= INFO_MIN & !is.na(BETA) & !is.na(SE)]
  if (nrow(pqtl_dt) < 30) {
    cat(sprintf("    [skip] pQTL 過濾後 %d SNPs\n", nrow(pqtl_dt))); return(NULL)
  }

  pqtl_dt[, pqtl_varbeta := SE^2]
  pqtl_dt[, ea_pqtl      := toupper(ALLELE1)]
  pqtl_dt[, oa_pqtl      := toupper(ALLELE0)]
  n_pqtl <- as.integer(pqtl_dt[["N"]][1])

  n_with_rsid <- sum(!is.na(pqtl_dt$rsid))
  cat(sprintf("    pQTL SNPs: %d (rsID覆蓋: %d)\n", nrow(pqtl_dt), n_with_rsid))
  flush(stdout())

  # ══════════════════════════════════════════════════════════════════════════
  # 策略 A: hg38（FinnGen）
  # ══════════════════════════════════════════════════════════════════════════
  if (genome == "hg38") {
    # NO LIFTOVER. pqtl_dt$chrpos is chrom:GENPOS read straight from the UKB-PPP archive, and GENPOS
    # is ALREADY GRCh38 -- verified against the archive itself, not against a comment: for A1BG's
    # chr19 member, 27 GENPOS values fall inside A1BG's GRCh38 span and none inside its GRCh37 span,
    # and lifting the GRCh37 position embedded in the variant ID (chr19:58,856,694) forward returns
    # exactly that row's GENPOS (58,345,328). The ID string carries a legacy GRCh37 position; GENPOS
    # does not. The call that used to stand here -- a function then named liftover_to_hg38(), which
    # no longer exists -- therefore lifted GRCh38 coordinates a SECOND time,
    # displacing them by 0 to 864,716 bp depending on the locus, and every position-fallback match
    # made from the result paired two different variants.

    # FinnGen 區域資料（含 rsid_out, chrpos_out）
    out_dt <- get_finngen_region(outcome_id, gene)
    if (is.null(out_dt) || nrow(out_dt) < 20) {
      cat(sprintf("    [skip] FinnGen SNP=%d\n", if (is.null(out_dt)) 0L else nrow(out_dt)))
      return(NULL)
    }

    # Standardise pQTL columns -- chrpos is GRCh38 and the outcome is GRCh38, so they are comparable.
    pqtl_sub <- data.table(
      rsid         = pqtl_dt$rsid,
      chrpos_pqtl  = pqtl_dt$chrpos,
      pqtl_beta    = pqtl_dt$BETA,
      pqtl_varbeta = pqtl_dt$pqtl_varbeta,
      ea_pqtl      = pqtl_dt$ea_pqtl,
      oa_pqtl      = pqtl_dt$oa_pqtl
    )

    # Standardise outcome columns
    out_sub <- out_dt[, .(
      rsid_out,
      chrpos_out,
      ea_out      = toupper(get(meta$ea_col)),
      oa_out      = toupper(get(meta$oa_col)),
      out_beta    = as.numeric(get(meta$beta_col)),
      out_varbeta = as.numeric(get(meta$se_col))^2
    )]
    out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

    common <- rsid_pos_merge(pqtl_sub, out_sub)

  # ══════════════════════════════════════════════════════════════════════════
  # 策略 B: hg37（Nielsen AF）
  # ══════════════════════════════════════════════════════════════════════════
  } else if (genome == "hg37") {
    # PROJECTED DOWN, not assumed. The comment here used to read "already hg37"; the archive says
    # otherwise, and that one false comment was the only thing making this branch look safe. The
    # outcome's positions are Pos_b37, so the pQTL positions have to be brought to GRCh37 before they
    # can be compared at all. Rows with no GRCh37 counterpart get NA and are excluded from the
    # position fallback, which is the correct answer for them.
    pqtl_37 <- liftover_to_hg37(pqtl_dt)
    .n37 <- if (!is.null(pqtl_37)) sum(!is.na(pqtl_37$chrpos37)) else 0L
    cat(sprintf("    liftover hg38->hg37: %d/%d\n", .n37, nrow(pqtl_dt)))
    # ...AND SO IS THE WINDOW. region_start/region_end are GRCh38, built from the GRCh38 gene body;
    # get_nielsen_region() filters the outcome's Pos_b37. Handing GRCh38 bounds to a GRCh37 file is
    # the same defect this rebuild exists to remove, with the builds swapped -- for DNER it queried
    # 2:228857519-230214586 instead of roughly 2:229722235-231079301, losing most of the true window
    # and adding an equally large wrong flank. The two ENDPOINTS are not lifted: a window's endpoints
    # can land in different chain blocks or in reverse order (PDZK1 does), which would invent a
    # multi-megabase region. The bounds are taken from the projected positions of the pQTLs actually
    # being colocalised, which is exactly the set of outcome rows this pair needs.
    # AN OFF-CHROMOSOME LIFT IS A PROPERTY OF ONE VARIANT, NOT OF THE PAIR. Aborting the whole run
    # because a single pQTL landed on another chromosome treated a legitimate chain outcome as a
    # data-foundation failure: that variant simply has no GRCh37 counterpart usable against this
    # outcome, so it is dropped, and the pair proceeds on the rest. What IS fatal is having NO
    # variant left on this chromosome -- then there is no window, and falling back to the GRCh38
    # bounds would query the wrong region silently. Dropping is done by clearing the coordinate, so
    # the variant is excluded from BOTH the rsID merge (keyed on rsid+chrpos) and the position
    # fallback, rather than only from the window bounds.
    .chr37  <- if (!is.null(pqtl_37)) sub(":.*", "", pqtl_37$chrpos37) else NA_character_
    .same37 <- .chr37 %in% as.character(chrom)
    .off37  <- !is.na(.chr37) & !.same37
    if (!any(.same37))
      fail_loud_data(sprintf(paste0("coloc %s ~ %s: the pQTL liftover to GRCh37 leaves NO variant on ",
                                    "chr%s (mapped=%d, off-chromosome=%d, unmapped=%d), so there is no ",
                                    "window. Falling back to the GRCh38 bounds would query the wrong ",
                                    "region silently."),
                             gene, outcome_id, chrom, sum(.same37), sum(.off37), sum(is.na(.chr37))))
    if (any(.off37)) {
      cat(sprintf("    [liftover] %d pQTL(s) lifted off chr%s -- dropped for this outcome\n",
                  sum(.off37), chrom))
      pqtl_37$chrpos37[.off37] <- NA_character_
    }
    .pos37 <- as.integer(sub(".*:", "", pqtl_37$chrpos37[.same37]))
    out_dt <- get_nielsen_region(outcome_id, chrom, min(.pos37), max(.pos37))
    if (is.null(out_dt) || nrow(out_dt) < 20) {
      cat(sprintf("    [skip] %s SNP=%d\n", outcome_id, if (is.null(out_dt)) 0L else nrow(out_dt)))
      return(NULL)
    }

    pqtl_src <- if (!is.null(pqtl_37)) pqtl_37 else pqtl_dt
    pqtl_sub <- data.table(
      rsid         = pqtl_src$rsid,
      chrpos_pqtl  = if (!is.null(pqtl_37)) pqtl_37$chrpos37 else rep(NA_character_, nrow(pqtl_dt)),
      pqtl_beta    = pqtl_src$BETA,
      pqtl_varbeta = pqtl_src$pqtl_varbeta,
      ea_pqtl      = pqtl_src$ea_pqtl,
      oa_pqtl      = pqtl_src$oa_pqtl
    )

    out_sub <- out_dt[, .(
      rsid_out,
      chrpos_out,
      ea_out      = toupper(get(meta$ea_col)),
      oa_out      = toupper(get(meta$oa_col)),
      out_beta    = as.numeric(get(meta$beta_col)),
      out_varbeta = as.numeric(get(meta$se_col))^2
    )]
    out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

    common <- rsid_pos_merge(pqtl_sub, out_sub)

  # ══════════════════════════════════════════════════════════════════════════
  # 策略 C: hg38 harmonized（MEGASTROKE / CKDGen）
  # ══════════════════════════════════════════════════════════════════════════
  } else if (genome == "hg38_hm") {
    # NO LIFTOVER, for the same reason as the hg38 branch above: GENPOS is already GRCh38 and the
    # harmonised outcome is GRCh38.
    out_dt <- get_harmonized_hg38_region(outcome_id, gene)
    if (is.null(out_dt) || nrow(out_dt) < 20) {
      cat(sprintf("    [skip] harmonized SNP=%d\n", if (is.null(out_dt)) 0L else nrow(out_dt)))
      return(NULL)
    }

    pqtl_src <- pqtl_dt
    pqtl_sub <- data.table(
      rsid         = pqtl_src$rsid,
      chrpos_pqtl  = pqtl_src$chrpos,
      pqtl_beta    = pqtl_src$BETA,
      pqtl_varbeta = pqtl_src$pqtl_varbeta,
      ea_pqtl      = pqtl_src$ea_pqtl,
      oa_pqtl      = pqtl_src$oa_pqtl
    )

    out_sub <- out_dt[, .(
      rsid_out,
      chrpos_out,
      ea_out      = toupper(get(meta$ea_col)),
      oa_out      = toupper(get(meta$oa_col)),
      out_beta    = as.numeric(get(meta$beta_col)),
      out_varbeta = as.numeric(get(meta$se_col))^2
    )]
    out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

    common <- rsid_pos_merge(pqtl_sub, out_sub)

  } else {
    cat(sprintf("    [skip] 未知 genome: %s\n", genome)); return(NULL)
  }

  # ── 共同後處理 ──────────────────────────────────────────────────────────
  cat(sprintf("    總合併: %d SNPs\n", nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] 合併後不足 20\n"); return(NULL) }

  # HARMONISE FIRST, THEN DROP WHAT IS STILL AMBIGUOUS. These two ran the other way round, so a
  # position that matched two outcome records -- say the pQTL's A/G against an A/G and an A/T, or any
  # multiallelic site -- had BOTH rows deleted before the alleles were ever looked at, even though
  # harmonisation would have kept the one compatible row and rejected the other. Order is the whole
  # fix; the ambiguity drop stays, because a group that is still multiple AFTER harmonisation is
  # genuinely unresolved and must not enter coloc twice under one identifier.
  common <- harmonize_alleles(common, "ea_pqtl", "oa_pqtl", "ea_out", "oa_out", "out_beta")
  cat(sprintf("    allele align: %d SNPs\n", nrow(common)))

  # A REUSED RSID IS NOT AN AMBIGUOUS VARIANT. snp_id is the rsID, so two distinct variants sharing
  # one rsID -- 1:1251285 G/A and G/T are both rs6657544 -- arrived here as a duplicated group and
  # were both deleted by the safety net below, even though harmonisation had just resolved each to
  # its own outcome row. Only the COLLIDING ones are re-keyed: an unambiguous rsID keeps its name, so
  # the SuSiE copy's rsID-based LD lookup is untouched, while a collided pair gets an allele-aware id
  # that its `^rs[0-9]+$` filter correctly excludes -- an LD panel keyed by rsID cannot represent two
  # variables with one name, and picking either arbitrarily would misalign z against R.
  common[duplicated(snp_id) | duplicated(snp_id, fromLast = TRUE),
         snp_id := paste(chrpos_pqtl, ea_pqtl, oa_pqtl, sep = ":")]

  # 去重（安全網）：有歧義的 snp_id（多列）整組丟掉，不隨機取一列
  n_before <- nrow(common)
  common <- common[, if (.N == 1L) .SD, by = snp_id]
  if (nrow(common) < n_before)
    cat(sprintf("    dedup snp_id: %d → %d SNPs (dropped %d ambiguous)\n",
                n_before, nrow(common), n_before - nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] 對齊後不足 20\n"); return(NULL) }

  # Non-destructive capture hook (only when COLOC_CAPTURE_DIR set): save the exact
  # coloc-ready vectors so an empirical-null (decoy-loci) can reuse run_coloc_abf identically.
  .capdir <- Sys.getenv("COLOC_CAPTURE_DIR", "")
  if (nzchar(.capdir)) {
    dir.create(.capdir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(gene = gene, outcome_id = outcome_id, n_pqtl = n_pqtl, meta = meta,
                 common = common[, .(snp_id, pqtl_beta, pqtl_varbeta, out_beta, out_varbeta)]),
            file = file.path(.capdir, sprintf("cap__%s__%s.rds", gene, gsub("[^A-Za-z0-9]+", "_", outcome_id))))
  }
  cr <- run_coloc_abf(common$pqtl_beta, common$pqtl_varbeta, n_pqtl,
                      common$out_beta, common$out_varbeta, meta,
                      snp_ids = common$snp_id)
  if (is.null(cr)) return(NULL)

  pp <- cr$summary
  # 2026-07-01 (editor REQUIRE): a NON-FINITE PP on a pair that REACHED coloc.abf is a degenerate result that
  # would otherwise be recorded status=="ok" and then silently filtered from the strong set (!is.na & >=0.80) --
  # an undetectable completeness hole (+ possible outcome-differential bias to Phase 5). Fail loud, don't skip.
  if (!is.finite(pp["PP.H4.abf"]) || !is.finite(pp["PP.H3.abf"]))
    fail_loud_data(sprintf("coloc %s ~ %s: coloc.abf returned NON-FINITE PP (H4=%s H3=%s) on a reached pair -- degenerate, refusing a silent status=ok", gene, outcome_id, pp["PP.H4.abf"], pp["PP.H3.abf"]))
  cat(sprintf("  → PP.H4=%.3f  PP.H3=%.3f  nSNPs=%d\n",
              pp["PP.H4.abf"], pp["PP.H3.abf"], as.integer(pp["nsnps"])))
  flush(stdout())

  data.table(
    gene               = gene,
    outcome_id         = outcome_id,
    n_snps_coloc       = as.integer(pp["nsnps"]),
    # ALL FIVE STORED UNROUNDED, not only H3/H4. Rounding H0-H2 to four decimals meant the published
    # table did not contain what coloc.abf returned and the five posteriors need not sum to one, so a
    # reader recomputing any quantity from them gets a different answer from the engine. Whatever is
    # displayed is a formatting choice at the point of display, not of storage.
    PP.H0              = as.numeric(pp["PP.H0.abf"]),
    PP.H1              = as.numeric(pp["PP.H1.abf"]),
    PP.H2              = as.numeric(pp["PP.H2.abf"]),
    # coloc_support below is computed from the full-precision posteriors, but 06,
    # the parallel summary and the phase-5/reporting scripts all RECOMPUTE the same rule from these
    # stored values -- so a raw H4 of 0.79996, which the engine calls unsupported, was stored as
    # 0.8000 and called supported downstream (23_merge:214 already detected the contradiction and
    # published the rounded call anyway). Fixing the producer makes every existing recomputation
    # agree; rounding for display is the reporting layer's job, not the engine's. H3 as well as H4,
    # because both enter the H4/H3 ratio arm of the rule.
    PP.H3              = as.numeric(pp["PP.H3.abf"]),
    PP.H4              = as.numeric(pp["PP.H4.abf"]),
    coloc_support      = is_coloc_pass(pp["PP.H4.abf"], pp["PP.H3.abf"]),   # D1 dual rule; centralized (PP finite: guaranteed by fail_loud above)
    independent_signal = pp["PP.H3.abf"] >= 0.80
  )
}

# ── Main ─────────────────────────────────────────────────────────────────────
main <- function() {
  rsid_map    <- load_rsid_map()
  pairs       <- load_phase2_suggestive()
  # THE PINNED SNAPSHOT. Reopening the path here loaded whatever was on disk at that moment, while
  # the certificate comparison and the checkpoint key both used the hash config pinned when it was
  # sourced -- so a rebuild between the two produced results computed from generation B and filed
  # under generation A.
  instruments <- read_phase1_instruments()

  out_dir <- PHASE3_DIR
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  pairs_unique <- unique(pairs[method %in% c("Inverse variance weighted", "Wald ratio"),
                               .(gene, outcome_id, outcome_name)])
  cat(sprintf("\nUnique pairs: %d\n\n", nrow(pairs_unique)))
  flush(stdout())

  checkpoint_file <- file.path(out_dir, "coloc_checkpoint_v7.rds")
  all_coloc  <- list()
  done_pairs <- character(0)
  pair_status <- list()   # 2026-06-30 gate obj #2: per-pair status (ok / skip:reason / error) -> no silent drop

  # 2026-07-01 Codex-phase3 #1: checkpoint PROVENANCE key = code version + config N/ncase. A pre-2026-07-01
  # checkpoint, or one built under a different outcome N, is discarded (its ABF used different inputs) even if
  # the pair set matches -- the prior done-pairs check below only catches a different PAIR set.
  # The COORDINATE TABLE is part of the provenance. Without it, rebuilding the table and resuming
  # gives completed pairs their OLD regions and new pairs the new ones, in one results object, with
  # nothing to say so -- and moving a cis window is exactly what a coordinate rebuild does.
  # From the VALIDATED read, not a separate hash of the live path. Hashing the path here and
  # validating it a few lines later are two reads: a 03b install landing between them pins table A
  # while the run then validates and uses table B, and the checkpoint records A.
  # ONE validated read supplies both the pin and the provenance key. Validating once for the pin and
  # again a few lines later meant a 03b install landing between them pinned "absent" (or table A)
  # while the run then used table B -- and every gene lookup failed against a pin nothing could match,
  # recording pair errors instead of stopping cleanly at startup.
  .cert0 <- .assert_coords_certified()
  .COORD_MD5 <- if (is.null(.cert0)) "absent" else .cert0$md5
  assign(".COLOC_RUN_COORD_MD5", .COORD_MD5, envir = .GlobalEnv)
  # PINNED for the whole run. get_gene_hg38_coords() reloads whenever the table changes underneath
  # it, which is right for a long-lived session but wrong inside one analysis: earlier pairs would
  # use table A and later pairs table B while the checkpoint recorded only A. Any change now stops
  # the run instead of quietly splitting it.
  assign(".COLOC_RUN_COORD_MD5", .COORD_MD5, envir = .GlobalEnv)
  # THE FILE'S OWN HASH, not a date someone has to remember to change. The literal stayed at
  # 2026-07-01 through every edit since, so a computation-changing change could resume a checkpoint
  # built before it.
  .CODE_MD5 <- tryCatch(unname(tools::md5sum(file.path(BASE_DIR, "R/coloc.R"))),
                        error = function(e) "unhashable")
  .PROV_KEY <- paste0("coloc|", substr(.CODE_MD5, 1, 12), "|N=", paste(OUTCOME_N_TOTAL, collapse=","), "|NC=",
                      paste(OUTCOME_N_CASE, collapse=","), "|COORDS=", .COORD_MD5,
                      # ...AND THE INSTRUMENT GENERATION AND THE WINDOW. 03b's output is coordinate
                      # driven, so its bytes can be identical across two Phase-1 generations while
                      # its certificate names different ones; and CIS_WINDOW_KB changes the region
                      # every pair is computed over. Without both, a checkpoint from one generation
                      # is resumed under another wherever the pair names overlap.
                      "|INSTR=", if (exists("INSTRUMENTS_MD5")) INSTRUMENTS_MD5 else "unpinned",
                      "|WIN=", if (exists("COLOC_REGION_KB")) COLOC_REGION_KB else "unset")
  # Checked ONCE, before the checkpoint is read and before any pair runs. It used to surface only
  # inside the per-pair tryCatch, which recorded every pair as an "error", let the run continue to
  # the end, and then wrote those errors into the checkpoint -- so restoring the certificate was not
  # enough, the checkpoint had to be deleted by hand as well.
  invisible(.assert_coords_certified())

  if (file.exists(checkpoint_file)) {
    cp <- readRDS(checkpoint_file)
    if (is.null(cp$prov_key) || !identical(cp$prov_key, .PROV_KEY)) {
      cat("[checkpoint] PROVENANCE mismatch (stale code/N, e.g. a pre-2026-07-01 checkpoint) -> discarding, fresh start\n"); cp <- list()
    }
    # 2026-07-01 (gate finding 3): validate the checkpoint belongs to THIS run — any done pair outside
    # the current suggestive set means stale/different inputs -> discard, don't reuse.
    # ...AND THE STATUSES, NOT ONLY THE DONE PAIRS. An errored pair is deliberately kept out of
    # done_pairs so it will be retried, but it stays in pair_status -- so a checkpoint whose only
    # foreign entry is an error passed this ownership check, its stale error was restored, and the
    # run aborted at completion on a pair the current universe does not contain.
    .cur_keys <- paste(pairs_unique$gene, pairs_unique$outcome_id, sep = "_")
    if (length(setdiff(c(cp$done_pairs, names(cp$pair_status)), .cur_keys)) > 0) {
      cat("[checkpoint] STALE (done pairs or statuses not in current run) -> discarding, fresh start\n"); cp <- list()
    }
    if (!is.null(cp$results) && length(cp$results) > 0) {
      all_coloc  <- cp$results
      done_pairs <- cp$done_pairs
      if (!is.null(cp$pair_status)) pair_status <- cp$pair_status
      cat(sprintf("Resuming: %d done, %d with results\n\n",
                  length(done_pairs), length(all_coloc)))
    } else {
      cat("Fresh start (empty checkpoint)\n\n")
    }
  } else {
    cat("Fresh start (no checkpoint)\n\n")
  }
  flush(stdout())

  for (i in seq_len(nrow(pairs_unique))) {
    g   <- pairs_unique$gene[i]
    oid <- pairs_unique$outcome_id[i]
    key <- paste(g, oid, sep = "_")

    if (key %in% done_pairs) {
      cat(sprintf("[%d/%d] SKIP: %s ~ %s\n", i, nrow(pairs_unique), g, oid))
      next
    }

    # THE ROW COUNT IS TAKEN BEFORE `[1]`, NOT AFTER. `dt[gene == g]` on an absent gene gives zero
    # rows, but `[1]` on a zero-row data.table returns ONE row of NA -- so nrow(gene_info) was never
    # 0 and this skip was unreachable. The pair went on with NA chrom/start/end and surfaced later as
    # an error, which aborts the whole run at the end-of-loop check instead of being recorded as the
    # skip it is. Every current suggestive gene is in the instrument table, so this is latent today.
    .gi_all <- instruments[gene == g]
    gene_info <- .gi_all[1]
    if (nrow(.gi_all) == 0) {
      cat(sprintf("[%d/%d] %s: 無 instrument 資料\n", i, nrow(pairs_unique), g))
      pair_status[[key]] <- "skip:no_instrument"; done_pairs <- c(done_pairs, key); next
    }

    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, nrow(pairs_unique), g, oid))
    flush(stdout())

    .cerr <- FALSE
    res <- tryCatch(
      run_coloc_pair(g, oid,
                     gene_info$gene_chrom, gene_info$gene_start, gene_info$gene_end,
                     rsid_map),
      error = function(e) {
        message(sprintf("  [ERROR] %s ~ %s: %s", g, oid, e$message)); .cerr <<- TRUE; NULL
      }
    )

    if (!is.null(res)) { all_coloc[[key]] <- res; pair_status[[key]] <- "ok" } else
      pair_status[[key]] <- if (.cerr) "error" else "skip"   # gate obj #2: error (bad) vs legit skip (low-SNP/no-region)
    # An ERRORED pair is not done. Recording it here meant the resume skipped it as completed while
    # the end-of-loop error check aborted on its status, so a transient failure poisoned every resume
    # until the checkpoint was deleted by hand. Retrying is not silent: the status is still written,
    # and a pair that fails again records "error" again and aborts again.
    if (!.cerr) done_pairs <- c(done_pairs, key)

    if (i %% 5 == 0 || i == nrow(pairs_unique)) {
      saveRDS(list(results = all_coloc, done_pairs = done_pairs, pair_status = pair_status, prov_key = .PROV_KEY), checkpoint_file)
      cat(sprintf("  [checkpoint: %d done, %d results]\n",
                  length(done_pairs), length(all_coloc)))
    }
    flush(stdout())
    Sys.sleep(0.2)
  }

  saveRDS(list(results = all_coloc, done_pairs = done_pairs, pair_status = pair_status, prov_key = .PROV_KEY), checkpoint_file)   # 2026-07-01 Codex-phase3 reverify: final save also needs prov_key (loader discards a checkpoint without it)
  # 2026-06-30 gate obj #2: persist per-pair status + FAIL LOUD if any suggestive pair ERRORED (vs legit skip).
  .ps <- data.table(pair = names(pair_status), status = unlist(pair_status))
  fwrite(.ps, file.path(out_dir, "coloc_pair_status_v7.csv"))
  .errp <- .ps[status == "error"]
  if (nrow(.errp) > 0) stop(sprintf("03_coloc: %d pair(s) ERRORED (not a legit skip): %s -- investigate; do not ship silent gaps",
                                    nrow(.errp), paste(.errp$pair, collapse = ", ")))

  if (length(all_coloc) == 0) {
    cat("\n結果：無任何 coloc 成功\n")
    return(invisible(NULL))
  }

  coloc_dt <- rbindlist(all_coloc, fill = TRUE)
  coloc_dt[, outcome_name := names(OUTCOMES)[match(outcome_id, OUTCOMES)]]
  saveRDS(coloc_dt, file.path(out_dir, "coloc_results_v7_stage1.rds"))

  cat("\n=== Phase 3 (v5) Summary ===\n")
  # "pairs with a posterior", not "pairs tested": only status=ok pairs produce a row, while attempted
  # is ok + skip. The archived generation attempted 1,368 and stored 1,084, so the old wording
  # understated the tested universe by 284 in the one place a reader would look for it.
  cat(sprintf("Pairs with a posterior: %d\n", nrow(coloc_dt)))
  .dual <- is_coloc_pass(coloc_dt$PP.H4, coloc_dt$PP.H3)   # D1 dual rule (centralized 2026-07-05)
  cat(sprintf("Dual rule (PP.H4>=0.80 & H4/H3>=5): %d pairs, %d proteins\n",
              sum(.dual, na.rm = TRUE),
              length(unique(coloc_dt$gene[.dual]))))
  cat(sprintf("PP.H3 >= 0.80: %d pairs\n",
              sum(coloc_dt$PP.H3 >= 0.80, na.rm = TRUE)))

  cat("\n--- Top Results (by PP.H4) ---\n")
  ord <- order(-coloc_dt$PP.H4)
  print(coloc_dt[ord[1:min(20, nrow(coloc_dt))],
                  c("gene", "outcome_id", "n_snps_coloc", "PP.H4", "PP.H3")])

  invisible(coloc_dt)
}

if (!exists("COLOC_SOURCE_ONLY") || !isTRUE(COLOC_SOURCE_ONLY)) coloc_results <- main()
