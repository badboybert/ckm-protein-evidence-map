# ============================================================================
# EAS AXIS 04: EAS pQTL instrument-transferability sidecar (contract SS5, SS10 item 4)
# ----------------------------------------------------------------------------
# Independent of the outcome axis: for each of the 62 frozen lead instruments,
# does JCTF (Japanese Olink) or CKB (Chinese Olink, PheWeb caches) carry an
# exact, allele-resolved SNP-protein association, with F_EAS=(beta/se)^2 and
# direction vs UKB? Writes the FULL per-source/per-assay sidecar
# (eas_pqtl_profile.csv, SS10.1 schema) and the deterministic any-source L4
# summary (eas_pqtl_l4_summary.csv, SS5.2: counts + flags + four-value
# direction class only -- NO representative source, assay, or F).
# BOTH registered sources are REQUIRED to be present and readable: a missing
# source stops the run rather than being summarized as "no evidence"
# (EASr2 finding 9). Duplicates must agree on effect stats AND alleles AND
# assay identity (finding 10). Nothing here reads outcome states, tiers, or
# the outcome profile.
# JCTF format (probed 2026-08-15 against the installed file): GTEx-style
#   variant_id_hg19  "1:26963434:G:A"  = chr:pos:REF:ALT on GRCh37
#   variant_id_hg38  "chr1:26636943:G:A"
#   rsid (often EMPTY), gene_id, gene_name, ..., slope, slope_se, ...
# slope is per ALT allele (GTEx convention); REF/ALT designate the genome
# reference -> records are strand-anchored. Matching: rsid where present,
# else GRCh37 chr:pos from variant_id_hg19 against 02's lead_b37_lookup.csv
# (02 must have run first; its absence stops the run).
# CKB caches lacking source URL / query date / trait ID provenance are refused
# (pre-spec) -- fail-closed, not silently skipped.
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
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
`%||%` <- function(a, b) if (is.null(a)) b else a
B <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"); setwd(B)
EAS_DIR <- file.path("results", paste0("repro_", PHASE1_SET), "eas_axis")
dir.create(EAS_DIR, showWarnings = FALSE, recursive = TRUE)
F_MIN <- 10

t4  <- fread(file.path("results", paste0("repro_", PHASE1_SET), "table_s_pair_evidence_tiered.csv"))
ins <- as.data.table(readRDS(file.path("results", "phase1_hg38", "ukb_cis_instruments.rds")))
lead <- unique(merge(t4[, .(gene)], ins[, .(gene, rsid = SNP, beta_ukb = beta,
                    ea_ukb = toupper(effect_allele), oa_ukb = toupper(other_allele))],
                    by = "gene"))
stopifnot(uniqueN(lead$rsid) > 0L)
LK37 <- file.path(EAS_DIR, "lead_b37_lookup.csv")
if (!file.exists(LK37))
  stop("04: ", LK37, " absent -- run 02 first; JCTF matching needs the leads' GRCh37 coordinates")
lk37 <- fread(LK37)[, .(rsid, chr37 = as.character(chr37), pos37 = as.integer(pos37))]
lead <- merge(lead, lk37, by = "rsid", all.x = TRUE)
stopifnot(sum(is.na(lead$pos37)) == 0L)
skey <- function(g, r) paste(g, r, sep = "|")

rows <- list()
row1 <- function(li, source, assay, build, map_status, exact, am, b, s, keep_reason) {
  Fe <- if (is.finite(b) && is.finite(s) && s > 0) (b / s)^2 else NA_real_
  data.table(
    gene = li$gene, instrument_rsID = li$rsid, eas_pqtl_source = source,
    assay_id = assay, build = build, protein_mapping_status = map_status,
    exact_match = exact, allele_match = am,
    betaX_EAS = b, seX_EAS = s, F_EAS = Fe,
    direction_same_as_UKB = ifelse(is.finite(b), sign(b) == sign(li$beta_ukb), NA),
    eligible_strong_row = isTRUE(am) && is.finite(Fe) && Fe >= F_MIN,
    eas_pqtl_sidecar_key = skey(li$gene, li$rsid),
    eas_pqtl_reason = if (!exact) "VARIANT_ABSENT" else if (!am) "ALLELE_MISMATCH"
                      else if (!is.finite(Fe)) "INVALID_EFFECT"
                      else if (Fe < F_MIN) "WEAK_F" else keep_reason)
}

# ── CKB PheWeb JSON caches: REQUIRED source ──────────────────────────────────
CKB <- file.path("data", "eas_pqtl", "ckb_json")
ckb_files <- list.files(CKB, pattern = "\\.json$", full.names = TRUE)
if (!dir.exists(CKB) || !length(ckb_files))
  stop("04: CKB cache directory absent or empty (", CKB, ") -- a registered source that ",
       "is unavailable stops the run; absence must never be summarized as no evidence")
# PER-LEAD, NOT JUST PER-DIRECTORY. The check above is global: with 63 caches present, one gene's cache
# going missing leaves the directory non-empty, every guard passes, and that gene receives only its
# JCTF row -- so a query that never ran is written out looking exactly like a queried-and-absent
# negative. The loop below iterates over the FILES, so a lead with no file is simply never visited and
# nothing notices. This asserts the other direction. A protein that CKB genuinely does not carry needs
# a cache recording that (the provenance fields are already mandatory below); it must not be
# represented by the absence of a file, because absence cannot be told apart from a failed extraction.
.ckb_have <- toupper(sub("\\.json$", "", basename(ckb_files)))
.ckb_want <- unique(toupper(lead$gene))
.ckb_miss <- setdiff(.ckb_want, .ckb_have)
if (length(.ckb_miss))
  stop(sprintf(paste0("04: %d of %d lead gene(s) have no CKB cache file in %s: %s. A missing cache is ",
                      "not a negative result -- either extract it, or write a cache that records the ",
                      "query returning nothing. Continuing would give these genes a JCTF-only profile ",
                      "indistinguishable from one where CKB was queried and had no signal."),
               length(.ckb_miss), length(.ckb_want), CKB,
               paste(utils::head(sort(.ckb_miss), 12), collapse = ", ")))
cat(sprintf("CKB caches: %d file(s) covering all %d lead gene(s)\n", length(ckb_files), length(.ckb_want)))
for (f in ckb_files) {
  # NAMED `.gsym`, NOT `gene`. Inside `lead[...]` the bare name `gene` resolves to the TABLE COLUMN,
  # not to this loop scalar, so `lead[toupper(lead$gene) == gene]` asked whether each row's uppercased
  # symbol equalled its own symbol -- true for every row whose symbol is already uppercase. Every one
  # of the 63 CKB caches was therefore applied to all 63 leads: 3,969 CKB profile rows instead of 63,
  # and 3,962 VARIANT_ABSENT instead of 56. The published EAS counts are unaffected only because the
  # summary de-duplicates by source and key downstream, and because no cross-applied cache happens to
  # contain another gene's lead on the current data.
  .gsym <- toupper(sub("\\.json$", "", basename(f)))
  if (!.gsym %in% toupper(lead$gene)) next
  j <- tryCatch(fromJSON(f, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(j)) stop("04: unparseable CKB cache: ", f)
  meta_ok <- all(c("source_url", "query_date", "trait_id") %in% names(j))
  if (!meta_ok)
    stop("04: CKB cache ", basename(f), " lacks provenance (source_url/query_date/trait_id) -- ",
         "the pre-spec refuses caches whose origin is unstated; re-query with provenance")
  v <- as.data.table(j$variants %||% j$data %||% list())
  if (nrow(v) && !any(c("rsid", "SNP") %in% names(v)))
    stop("04: CKB cache ", basename(f), " has no rsid/SNP column -- cannot match; ",
         "re-query with a variant identifier, do not guess")
  vid <- if ("rsid" %in% names(v)) "rsid" else "SNP"
  gset <- lead[toupper(lead$gene) == .gsym]
  for (k in seq_len(nrow(gset))) {
    li <- gset[k]
    hit <- if (nrow(v)) v[get(vid) == li$rsid] else v
    if (nrow(hit) > 1L) {
      # agreement must cover stats AND alleles (EASr2 finding 10); assay is the
      # cache's single trait_id, identical within one file by construction
      hb <- unique(data.table(b = as.character(hit$beta),
                              s = as.character(hit$sebeta %||% hit$se),
                              ea = toupper(as.character(hit$alt %||% hit$effect_allele)),
                              oa = toupper(as.character(hit$ref %||% hit$other_allele))))
      if (nrow(hb) > 1L)
        stop("04: CKB duplicate rows disagree for ", li$gene, "/", li$rsid, " -- refusing (SS5.2)")
      hit <- hit[1]
    }
    exact <- nrow(hit) == 1L
    b <- if (exact) suppressWarnings(as.numeric(hit$beta[1])) else NA_real_
    s <- if (exact) suppressWarnings(as.numeric(hit$sebeta[1] %||% hit$se[1])) else NA_real_
    am <- if (exact) {
      ea <- toupper(hit$alt[1] %||% hit$effect_allele[1] %||% "")
      oa <- toupper(hit$ref[1] %||% hit$other_allele[1] %||% "")
      al <- (ea == li$ea_ukb & oa == li$oa_ukb) | (ea == li$oa_ukb & oa == li$ea_ukb)
      if (isTRUE(al) && ea == li$oa_ukb) b <- -b
      isTRUE(al)
    } else FALSE
    # build read from the cache's own provenance (the 2026-08-16 re-query is
    # GRCh38 -- verified via the panel's ABO top-hit coordinate), never assumed
    rows[[length(rows) + 1L]] <- row1(li, "CKB", as.character(j$trait_id %||% NA_character_),
                                      as.character(j$build %||% NA_character_),
                                      "MAPPED_BY_CACHE_NAME", exact, am, b, s, "OK")
  }
}

# ── JCTF cis-pQTL sumstats: REQUIRED source ──────────────────────────────────
JCTF <- file.path("data", "eas_pqtl", "jctf",
                  "pqtl_sumstats_wang_qs_et_al_2024_japan_covid19_taskforce.tsv.gz")
if (!file.exists(JCTF))
  stop("04: JCTF file absent (", JCTF, ") -- a registered source that is unavailable ",
       "stops the run; absence must never be summarized as no evidence")
jd <- fread(JCTF, showProgress = FALSE)
need <- c("variant_id_hg19", "rsid", "gene_name", "slope", "slope_se")
if (!all(need %in% names(jd)))
  stop("04: JCTF header mismatch -- expected [", paste(need, collapse = ","),
       "], got [", paste(head(names(jd), 15), collapse = ","), "]; re-probe and amend")
jd[, c("chr37j", "pos37j", "ref_j", "alt_j") := {
  sp <- tstrsplit(variant_id_hg19, ":")
  .(as.character(sp[[1]]), suppressWarnings(as.integer(sp[[2]])),
    toupper(sp[[3]]), toupper(sp[[4]]))
}]
if (jd[, mean(!is.na(pos37j) & !is.na(ref_j) & !is.na(alt_j))] < 0.99)
  stop("04: JCTF variant_id_hg19 does not parse as chr:pos:REF:ALT on >=99% of rows -- ",
       "the format assumption is wrong; re-probe and amend")
for (k in seq_len(nrow(lead))) {
  li <- lead[k]
  gsub_rows <- jd[toupper(gene_name) == toupper(li$gene)]
  hit <- gsub_rows[!is.na(rsid) & nzchar(rsid) & rsid == li$rsid]
  if (!nrow(hit)) hit <- gsub_rows[chr37j == li$chr37 & pos37j == li$pos37]
  if (nrow(hit) > 1L) {
    # agreement over stats AND alleles AND assay identity (gene_id)
    hb <- unique(hit[, .(slope, slope_se, ref_j, alt_j, gene_id)])
    if (nrow(hb) > 1L)
      stop("04: JCTF duplicate rows disagree for ", li$gene, "/", li$rsid, " -- refusing (SS5.2)")
    hit <- hit[1]
  }
  exact <- nrow(hit) == 1L
  b <- if (exact) as.numeric(hit$slope) else NA_real_        # slope per ALT allele
  s <- if (exact) as.numeric(hit$slope_se) else NA_real_
  am <- if (exact) {
    ea <- hit$alt_j; oa <- hit$ref_j
    al <- (ea == li$ea_ukb & oa == li$oa_ukb) | (ea == li$oa_ukb & oa == li$ea_ukb)
    if (isTRUE(al) && ea == li$oa_ukb) b <- -b
    isTRUE(al)
  } else FALSE
  rows[[length(rows) + 1L]] <- row1(li, "JCTF",
                                    if (exact) as.character(hit$gene_id) else NA_character_,
                                    "GRCh37", "MAPPED", exact, am, b, s, "OK")
}
rm(jd); invisible(gc())

prof <- rbindlist(rows, fill = TRUE)
PROFILE_COLS <- c("gene", "instrument_rsID", "eas_pqtl_source", "assay_id", "build",
                  "protein_mapping_status", "exact_match", "allele_match",
                  "betaX_EAS", "seX_EAS", "F_EAS", "direction_same_as_UKB",
                  "eligible_strong_row", "eas_pqtl_sidecar_key", "eas_pqtl_reason")
stopifnot(nrow(prof) > 0L, identical(sort(names(prof)), sort(PROFILE_COLS)),
          all(c("CKB", "JCTF") %in% prof$eas_pqtl_source))
setcolorder(prof, PROFILE_COLS)
fwrite(prof, file.path(EAS_DIR, "eas_pqtl_profile.csv"))

# ── deterministic any-source L4 summary (SS5.2): counts + flags only ─────────
strong <- prof[eligible_strong_row == TRUE]
summ <- unique(lead[, .(gene, instrument_rsID = rsid)])
agg <- strong[, .(eas_pqtl_n_sources_strong = uniqueN(eas_pqtl_source),
                  eas_pqtl_n_assays_strong = uniqueN(paste(eas_pqtl_source, assay_id)),
                  n_same = sum(direction_same_as_UKB == TRUE, na.rm = TRUE),
                  n_opp  = sum(direction_same_as_UKB == FALSE, na.rm = TRUE)),
              by = .(gene, instrument_rsID)]
summ <- merge(summ, agg, by = c("gene", "instrument_rsID"), all.x = TRUE)
anyx <- prof[exact_match == TRUE, .(eas_pqtl_any_exact = TRUE), by = .(gene, instrument_rsID)]
summ <- merge(summ, anyx, by = c("gene", "instrument_rsID"), all.x = TRUE)
summ[is.na(eas_pqtl_any_exact), eas_pqtl_any_exact := FALSE]
summ[is.na(eas_pqtl_n_sources_strong), `:=`(eas_pqtl_n_sources_strong = 0L,
                                            eas_pqtl_n_assays_strong = 0L,
                                            n_same = 0L, n_opp = 0L)]
summ[, eas_pqtl_any_strong := eas_pqtl_n_sources_strong > 0L]
# contract SS5.2 vocabulary VERBATIM (EASr3 finding 1: round 2 wrongly renamed it)
summ[, eas_pqtl_direction_class := fifelse(!eas_pqtl_any_strong |
                                             (n_same == 0L & n_opp == 0L), "NO_STRONG_EAS_PQTL",
                                    fifelse(n_same > 0L & n_opp > 0L, "MIXED",
                                    fifelse(n_opp == 0L, "ALL_SAME", "ALL_OPPOSITE")))]
summ[, eas_pqtl_sidecar_key := skey(gene, instrument_rsID)]
summ[, c("n_same", "n_opp") := NULL]
stopifnot(anyDuplicated(summ, by = c("gene", "instrument_rsID")) == 0L)
fwrite(summ, file.path(EAS_DIR, "eas_pqtl_l4_summary.csv"))
cat(sprintf("04: %d profile rows; %d/%d instruments any_exact; %d any_strong\n",
            nrow(prof), sum(summ$eas_pqtl_any_exact), nrow(summ), sum(summ$eas_pqtl_any_strong)))
print(summ[, .N, by = eas_pqtl_direction_class])
