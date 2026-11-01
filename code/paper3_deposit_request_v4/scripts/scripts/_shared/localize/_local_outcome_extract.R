# ============================================================================
# Shared LOCAL outcome extractor — returns, for a set of rsIDs and one outcome id,
# a data.frame in TwoSampleMR::extract_outcome_data() format, read from the LOCAL
# full sumstats (same files coloc uses). Direct match by DEFAULT; OPTIONAL LD-proxy substitution via proxies=/
# rsq= args (still no API). P7-F4: Phase 7's direct-only behaviour is because its wrapper does not forward those
# args, NOT a helper limitation — the proxy path + provenance fields are implemented below.
# Reuses the LOCAL_META + lookup logic from localize/02. Source this, then call
# local_extract_outcome_data(snps, outcome_id). Requires OUTCOME_DIR (from config).
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
suppressPackageStartupMessages(library(data.table))

# 2026-06-30 FIX (UKB-free overhaul, ref DATA_SOURCE_FINAL_TABLE.md + REMEDIATION_PLAN B1/B3):
#   AF  : nielsen2018 (UKB) -> finngen_R10_I9_AF.gz (UKB-free); ncase 50,743; N 261,395
#   CKD : Pattaro CKDGen_CKD_harmonised_hg38 (deleted, mixed-anc) -> CKDGen_CKD_Wuttke2019_EA.txt.gz (EUR); ncase 41,395; N 480,698
#   FinnGen ncase -> R10 manifest authoritative (CHD 46,959/HF 29,672/T2D 65,085/NAFLD 2,568); T2D N 400,197
# NOTE: this Phase7 helper duplicates outcome metadata; central-config consolidation (C2) is a later round.
.LO_fg <- function(file, ncase, n = 412181L) list(file = file, sep = "\t", rsid = "rsids", ea = "alt", oa = "ref",
    beta = "beta", se = "sebeta", eaf = "af_alt", pval = "pval", n_fixed = n, ncase = ncase,
    type = "cc", multi_rsid = TRUE)
.LO_META <- list(
  "FinnGen_R10_I9_CHD"       = .LO_fg("finngen_R10_I9_CHD.gz",       46959L),
  "FinnGen_R10_I9_HEARTFAIL" = .LO_fg("finngen_R10_I9_HEARTFAIL.gz", 29672L),
  "FinnGen_R10_T2D"          = .LO_fg("finngen_R10_T2D.gz",          65085L, n = 400197L),
  "FinnGen_R10_NAFLD"        = .LO_fg("finngen_R10_NAFLD.gz",         2568L),
  "FinnGen_R10_I9_AF"        = .LO_fg("finngen_R10_I9_AF.gz",        50743L, n = 261395L),
  "MEGASTROKE_EUR_GCST006906" = list(file = "MEGASTROKE_EUR_harmonised_hg38.h.tsv.gz", sep = "\t",
        rsid = "hm_rsid", rsid_fallback = "variant_id", ea = "hm_effect_allele", oa = "hm_other_allele",
        beta = "hm_beta", se = "standard_error", eaf = "hm_effect_allele_frequency", pval = "p_value",
        # All-or-nothing per row: an unharmonisable row has every hm_ field NA while the submitted
        # columns are complete, so substituting only the identifier left the row with no alleles and
        # no beta -- it was dropped, and a measured instrument looked absent (localize/02 too).
        ea_fallback = "effect_allele", oa_fallback = "other_allele",
        beta_fallback = "beta", eaf_fallback = "effect_allele_frequency",
        n_fixed = 446696L, ncase = 40585L, type = "cc"),
  "Wuttke2019_CKD_EA"          = list(file = "CKDGen_CKD_Wuttke2019_EA.txt.gz", sep = " ", rsid = "RSID",
        ea = "Allele1", oa = "Allele2", beta = "Effect", se = "StdErr", eaf = "Freq1", pval = "P-value",
        n_col = "n_total_sum", n_fixed = 480698L, ncase = 41395L, type = "cc"),
  "Wuttke2019_eGFR_EA"    = list(file = "CKDGen_eGFR_EA_hg37.txt.gz", sep = " ", rsid = "RSID",
        ea = "Allele1", oa = "Allele2", beta = "Effect", se = "StdErr", eaf = "Freq1", pval = "P-value",
        n_col = "n_total_sum", type = "quant")
)

# 2026-07-01: N/ncase READ from config's SINGLE SOURCE (OUTCOME_N_TOTAL/OUTCOME_N_CASE) -- overwrites the
# inline literals above so this Phase7/decode helper never drifts from config/01_mr/coloc (this IS the
# "central-config consolidation (C2)" the header noted). Callers MUST source config first (it needs OUTCOME_DIR).
if (!exists("OUTCOME_N_TOTAL") || !exists("OUTCOMES")) stop("_local_outcome_extract: config not sourced first (needs OUTCOME_N_TOTAL/OUTCOMES)")
.lo_id2nm <- setNames(names(OUTCOMES), unname(OUTCOMES))
for (.lo_oid in names(.LO_META)) { .lo_nm <- unname(.lo_id2nm[.lo_oid]); if (is.na(.lo_nm)) next   # `[` not `[[`: non-config id -> NA, skip (not a throw)
  .LO_META[[.lo_oid]]$n_fixed <- as.integer(OUTCOME_N_TOTAL[[.lo_nm]])
  if (isTRUE(.LO_META[[.lo_oid]]$type == "cc")) {
    if (is.na(OUTCOME_N_CASE[[.lo_nm]])) stop(sprintf("_local_outcome_extract: config OUTCOME_N_CASE is NA for cc outcome %s", .lo_oid))
    .LO_META[[.lo_oid]]$ncase <- as.integer(OUTCOME_N_CASE[[.lo_nm]]) } }

# read outcome file -> lookup keyed by rsID, restricted to `keep` rsIDs (no API)
.LO_lookup <- function(m, keep) {
  f <- file.path(OUTCOME_DIR, m$file); if (!file.exists(f)) stop("_local_outcome_extract: outcome sumstats file not found: ", f)
  want <- c(m$rsid, m$ea, m$oa, m$beta, m$se, m$pval)
  if (!is.null(m$eaf)) want <- c(want, m$eaf)
  if (!is.null(m$n_col)) want <- c(want, m$n_col)
  for (.fb in c("rsid_fallback", "ea_fallback", "oa_fallback", "beta_fallback", "eaf_fallback"))
    if (!is.null(m[[.fb]])) want <- c(want, m[[.fb]])
  dt <- fread(f, sep = m$sep, select = unique(want), showProgress = FALSE)
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
    na <- is.na(rid) | rid == "" | rid == "NA" | !is.finite(.hb)
    # ALIGNED WITH THE DISCOVERY SCREEN (phase2audit F-04, PI ruling 2026-08-17: align, do not merely
    # disclose). cascade/01_mr_screen.R:514-521 also falls back when the harmonised ALLELES or EAF are
    # unusable, not only the rsID or beta -- the harmonised and submitted fields are swapped as a set,
    # so a row with a fine harmonised beta but unusable harmonised alleles was left half-usable here
    # and complete there. Measured consequence of the narrower rule: NADK and CNTNAP2 produced 10
    # sensitivity estimates that differed from discovery on an unchanged instrument set, which makes a
    # threshold-sensitivity arm answer two questions at once.
    for (.p in list(c("ea", "ea_fallback"), c("oa", "oa_fallback"))) {
      if (is.null(m[[.p[1]]]) || is.null(m[[.p[2]]])) next
      .a <- toupper(as.character(dt[[m[[.p[1]]]]]))
      na <- na | is.na(.a) | !grepl("^[ACGT]+$", .a)
    }
    if (!is.null(m$eaf) && !is.null(m$eaf_fallback)) {
      .e <- suppressWarnings(as.numeric(dt[[m$eaf]]))
      na <- na | !is.finite(.e) | .e < 0 | .e > 1
    }
    rid[na] <- fb[na]
    for (.p in list(c("ea", "ea_fallback"), c("oa", "oa_fallback"), c("beta", "beta_fallback"), c("eaf", "eaf_fallback")))
      if (!is.null(m[[.p[1]]]) && !is.null(m[[.p[2]]])) { .v <- as.character(dt[[m[[.p[1]]]]]); .w <- as.character(dt[[m[[.p[2]]]]]); .v[na] <- .w[na]; set(dt, j = m[[.p[1]]], value = .v) } }
  if (isTRUE(m$multi_rsid)) {
    multi <- grepl(",", rid, fixed = TRUE)
    if (any(multi)) { idx <- which(multi)
      expl <- rbindlist(lapply(idx, function(i){ toks <- trimws(strsplit(rid[i], ",")[[1]]); toks <- toks[toks %in% keep]; if (length(toks)) data.table(row = i, rs = toks) else NULL }))
      if (nrow(expl)) { add <- dt[expl$row]; dt <- rbind(dt, add); rid <- c(rid, expl$rs) } }
  }
  k <- rid %in% keep; dt <- dt[k]; rid <- rid[k]
  lk <- data.table(rs = rid, ea = toupper(as.character(dt[[m$ea]])), oa = toupper(as.character(dt[[m$oa]])),
    beta = as.numeric(dt[[m$beta]]), se = as.numeric(dt[[m$se]]), pval = as.numeric(dt[[m$pval]]),
    eaf = if (!is.null(m$eaf)) suppressWarnings(as.numeric(dt[[m$eaf]])) else NA_real_,
    n   = if (!is.null(m$n_col)) suppressWarnings(as.numeric(dt[[m$n_col]])) else m$n_fixed)
  lk <- lk[!is.na(beta) & !is.na(se) & se > 0]
  if (!is.null(m$n_col)) { .n0 <- nrow(lk); lk <- lk[is.finite(n) & n > 0]   # 2026-07-01 (same silent-NA bug as 01_mr Codex-phase2 #3): n_col per-SNP N must be finite/positive
    if (.n0 > 0 && nrow(lk) < 0.95 * .n0) stop(sprintf("_local_outcome_extract: outcome %s n_col='%s' -- %d/%d SNP non-finite/non-positive N (malformed sample-size column)", m$file, m$n_col, .n0 - nrow(lk), .n0)) }
  # ONE ROW PER (rsID, unordered allele pair) -- see cascade/01_mr_screen.R. Deduplicating on the
  # rsID alone runs before anything knows which allele pair is wanted, so at a multiallelic site the
  # most significant row won and the instrument's actual variant was discarded, then counted as
  # "directly present", which also suppressed the proxy fallback.
  lk[, .apair := paste(pmin(toupper(ea), toupper(oa)), pmax(toupper(ea), toupper(oa)))]
  lk <- lk[order(pval)][!duplicated(paste(rs, .apair))]; setkey(lk, rs); lk
}

# ── LOCAL LD-proxy substitution (Codex backfill-review fix #4a: honour proxies=TRUE OFFLINE) ──────────
# Mirrors 02_extract_outcomes_local.R's orient_proxy + substitution using the 03_build_ld_proxies table, so the
# offline extractor recovers instruments missing from an outcome via an r2>=0.8 in-phase LD proxy. Same logic +
# same prox table + same sumstats -> consistent with the main corrected cascade (16_backfill53's read_local_outcome,
# which is proxy-included). Loaded lazily (only when proxies=TRUE) so direct-only callers pay no cost.
.LO_prox <- NULL; .LO_inst_al <- NULL
# The instrument alleles are needed by the DIRECT-match gate too, which runs whether or not proxies
# were asked for. Loading them only inside .lo_load_proxy() left the default proxies=FALSE path
# dereferencing a NULL data.table -- "could not find function \".\"" -- so direct-only extraction
# could not return a single row.
.lo_load_inst_al <- function() {
  if (!is.null(.LO_inst_al)) return(invisible())
  if (!exists("INSTRUMENTS_FILE") || !file.exists(INSTRUMENTS_FILE))
    stop("_local_outcome_extract: INSTRUMENTS_FILE is needed for the direct-match allele gate ",
         "(an outcome row is only the same variant if its allele pair is the instrument's)")
  .inst <- as.data.table(readRDS(INSTRUMENTS_FILE))
  # Not rs-format: SNP is the outcome-space name (see local_extract_outcome_data). An instrument
  # with no row here has no allele pair, so it can be neither direct-matched nor proxied.
  # eaf_exp CARRIES THE EXPOSURE-SIDE FREQUENCY for the palindromic ambiguity test, defined exactly
  # as the discovery screen defines it (cascade/01_mr_screen.R:332-334): the EAF closest to 0.5,
  # i.e. the least informative one, so a variant that is ambiguous under ANY of its records is
  # treated as ambiguous. NA when the instrument table carries no eaf, which makes the test fire --
  # the conservative direction, since an untestable palindrome is one harmonisation would drop.
  # ie/io are unchanged: `unique(by="SNP")` and `[1]` per group both take the first record.
  .LO_inst_al <<- .inst[!is.na(SNP) & nzchar(SNP),
                        .(ie = toupper(effect_allele)[1], io = toupper(other_allele)[1],
                          eaf_exp = if ("eaf" %in% names(.inst)) eaf[which.min(abs(eaf - 0.5))][1] else NA_real_),
                        by = SNP]
  setkey(.LO_inst_al, SNP)
}
.lo_load_proxy <- function() {
  if (!is.null(.LO_prox)) return(invisible())
  pf <- paste0(LD_PROXY_PREFIX, ".rds")   # corrected 1682 (previous generation) cascade rebuilds the _1682 table (not the legacy instr_ld_proxies.rds)
  if (!file.exists(pf)) stop("_local_outcome_extract: proxies=TRUE but LD-proxy table absent: ", pf, " (run cascade/01_mr_screen.R, which writes the _1682 cache)")
  if (!exists("INSTRUMENTS_FILE") || !file.exists(INSTRUMENTS_FILE)) stop("_local_outcome_extract: proxies=TRUE needs INSTRUMENTS_FILE (instrument alleles)")
  .inst <- as.data.table(readRDS(INSTRUMENTS_FILE))
  # THE CACHE IS IN PANEL SPACE, and every join below is on outcome-space SNP. Neither column was
  # translated here, so an instrument whose panel name differs from its outcome name (ADM is
  # rs113773056 in the panel, rs2403310 in the outcomes) was looked up under the wrong namespace,
  # and a panel identifier could be asked of a GWAS that does not use it. The cache is also not
  # tied to the instrument set in use: built for an earlier generation it simply has no row for a
  # new instrument, which then disappears from every outcome lacking its direct name -- silently.
  source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_ld_namespace.R"))
  # The producer names the sidecar from the stem, not from the .rds path: cascade/01 and 03 both
  # write paste0(out, ".provenance.rds"). Building it from `pf` asked for
  # instr_ld_proxies_1682.rds.provenance.rds, which nothing writes, so a correctly produced cache
  # would be rejected as having no provenance the moment the instrument table carried rsid_ld_panel.
  .pv <- paste0(sub("\\.rds$", "", pf), ".provenance.rds")
  if ("rsid_ld_panel" %in% names(.inst)) {
    if (!file.exists(.pv))
      stop("_local_outcome_extract: the instrument table carries rsid_ld_panel but ", basename(pf),
           " has no provenance sidecar, so it cannot be shown to cover this instrument set in the ",
           "panel namespace. Rebuild via cascade/01_mr_screen.R, which is what writes the _1682 cache (03 writes the unversioned instr_ld_proxies.*).")
    .pp <- readRDS(.pv)
    if (!identical(tolower(.pp$cache_md5 %||% ""), tolower(unname(tools::md5sum(pf)))))
      stop("_local_outcome_extract: ", basename(pf), " is not the cache its provenance describes.")
    if (!identical(.pp$index_md5 %||% "", ld_query_index_md5(.inst)))
      stop("_local_outcome_extract: the proxy cache was built for a different instrument set; ",
           "instruments it never queried have no proxy row and vanish silently. Rebuild via cascade/01_mr_screen.R, which is what writes the _1682 cache.")
  }
  # THE PANEL AND THE PARAMETERS, not just the identity of the cache and its index. r2, phase and
  # the proxy's identity are all properties of the genotypes it was computed from; a cache built
  # against a different BED, a different PLINK, or a looser threshold answers with a variant that
  # does not stand for the instrument, and both hashes above still match. Checked whenever a sidecar
  # exists, so a legacy cache without one keeps working exactly as before.
  if (file.exists(.pv)) {
    .pp <- readRDS(.pv)
    for (.k in c("bfile_bim_md5", "bfile_bed_md5", "bfile_fam_md5")) {
      .file <- file.path(LD_EUR_DIR, sub("^bfile_(...)_md5$", "EUR.\\1", .k))
      if (!is.null(.pp[[.k]]) && !identical(tolower(.pp[[.k]]), tolower(unname(tools::md5sum(.file)))))
        stop("_local_outcome_extract: the proxy cache was computed against a different ",
             basename(.file), " than the one configured; every r2 in it belongs to other genotypes.")
    }
    if (!is.null(.pp$plink_md5) && !is.na(.pp$plink_md5) && file.exists(PLINK_BIN) &&
        !identical(tolower(.pp$plink_md5), tolower(unname(tools::md5sum(PLINK_BIN)))))
      stop("_local_outcome_extract: the proxy cache was built by a different PLINK binary than the ",
           "configured one; defaults and tie-breaking are not guaranteed to agree.")
    .want_ns <- if ("rsid_ld_panel" %in% names(.inst)) "panel" else "legacy_snp"
    if (!is.null(.pp$namespace) && !identical(.pp$namespace, .want_ns))
      stop("_local_outcome_extract: the proxy cache was queried in the '", .pp$namespace,
           "' namespace but this instrument table is '", .want_ns, "'.")
  }
  .LO_prox <<- as.data.table(readRDS(pf))
  # Panel space -> outcome space, on BOTH columns, before any join.
  .map <- ld_to_snp_map(.inst)
  .LO_prox[, `:=`(instrument = ld_names_to_snp(instrument, .map),
                  proxy      = ld_names_to_snp(proxy, .map))]
  # NA is not a name -- the same removal 02 and 41b already do, missing only here. Left in, an NA
  # proxy reaches `keep`, and `NA %in% c(NA)` is TRUE in R, so an outcome row with a missing rsID
  # matches it: FinnGen T2D alone has 1,489,577 such rows, any one of which could be attached to an
  # instrument as its proxy on an allele-pair match and nothing else.
  .LO_prox <<- .LO_prox[!is.na(instrument) & !is.na(proxy)]
  setkey(.LO_prox, instrument)
  .lo_load_inst_al()
}
.LO_CMP <- c(A = "T", T = "A", C = "G", G = "C")

# Orient a proxy outcome row (po) to represent instrument effect allele ie (returns list or NULL) -- copy of 02's orient_proxy
#
# TWO ORIENTERS, ONE GATED. `full = FALSE` is the historical body, byte-for-byte, so every caller
# that does not opt in behaves exactly as before. `full = TRUE` is the discovery screen's
# orient_proxy() (cascade/01_mr_screen.R:606-653) COPIED VERBATIM -- it carries a coupling bound
# that must not be re-derived. Two things the historical body lacks and the discovery one has:
#   (1) a palindrome branch that orients from frequencies instead of refusing outright, and
#   (2) a reverse-complement fallback on the pair test.
# Without (1) the proxy path is strictly stricter than the direct path for the same class of
# variant (4,473 of 30,063 SNV pairs in the panel cache are palindromic); without (2) a valid
# opposite-strand proxy is filed as "no compatible variant".
.lo_orient_proxy <- function(ie, pr, po, full = FALSE) {
  p_inphase <- if (ie == pr$a_hap1) pr$b_hap1 else if (ie == pr$a_hap2) pr$b_hap2 else NA_character_
  if (is.na(p_inphase)) return(NULL)                  # instrument allele not on reference strand -> skip (conservative)
  # BOTH panel alleles must be the outcome row's pair (see 02_extract_outcomes_local.R): matching only
  # the in-phase allele accepts a panel A/G proxy against an outcome row recorded as A/T, carrying an
  # A-versus-T beta as if it were A-versus-G.
  p_other <- if (ie == pr$a_hap1) pr$b_hap2 else pr$b_hap1
  if (is.na(p_other)) return(NULL)

  if (!isTRUE(full)) {
    if (!setequal(toupper(c(p_inphase, p_other)), toupper(c(po$ea, po$oa)))) return(NULL)
    return(if (po$ea == p_inphase) list(beta = po$beta, eaf = po$eaf, se = po$se, n = po$n)
           else                    list(beta = -po$beta, eaf = if (is.na(po$eaf)) NA_real_ else 1 - po$eaf, se = po$se, n = po$n))
  }

  # ── verbatim from cascade/01_mr_screen.R:619-653 ──────────────────────────────────────────────
  .cmp <- .LO_CMP
  .pu <- toupper(c(p_inphase, p_other)); .ou <- toupper(c(po$ea, po$oa))
  # A PALINDROME CANNOT BE ORIENTED FROM ITS ALLELES ALONE, but it can be from FREQUENCIES -- which
  # is what harmonise_data(action=2) does for the direct path. The rule below is TwoSampleMR's, not
  # a new one: missing or strictly inside (0.42, 0.58) is uninformative; otherwise the two
  # frequencies agreeing on which side of 0.5 they fall means the literal mapping, disagreeing means
  # the complement.
  #
  # PHASE ALONE IS NOT ENOUGH: b_hap1 is B1 under coupling and B2 under repulsion, but this cache is
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
    return(if (toupper(po$ea) == .expect) list(beta = po$beta, eaf = po$eaf, se = po$se, n = po$n)
           else list(beta = -po$beta, eaf = if (is.na(po$eaf)) NA_real_ else 1 - po$eaf, se = po$se, n = po$n))
  }
  .expect <- toupper(p_inphase)
  if (!setequal(.pu, .ou)) {
    if (!all(.pu %in% names(.cmp)) || !setequal(unname(.cmp[.pu]), .ou)) return(NULL)
    .expect <- unname(.cmp[toupper(p_inphase)])
  }
  if (toupper(po$ea) == .expect) list(beta = po$beta, eaf = po$eaf, se = po$se, n = po$n)
  else list(beta = -po$beta, eaf = if (is.na(po$eaf)) NA_real_ else 1 - po$eaf, se = po$se, n = po$n)
}

# extract_outcome_data()-compatible LOCAL replacement: DIRECT matches + optional r2>=0.8 LD proxies.
# proxies defaults FALSE so each caller OPTS IN explicitly (avoids a silent proxy blast-radius across shared consumers
# e.g. phase7/deCODE); the phase2 5e-8 sensitivity callers pass proxies=TRUE to stay consistent with the proxy-included
# main cascade. Proxy provenance (proxy.outcome/proxy_snp.outcome/proxy_r2.outcome) is recorded on every row.
# palindromic_fallthrough OPTS INTO THE DISCOVERY SCREEN'S TWO EXTRA RULES, and defaults FALSE for
# the same reason `proxies` does: every existing caller stays byte-identical. Without it, a direct
# palindromic row whose strand cannot be resolved is counted as a direct match here, removed from
# `missing`, and so blocks the proxy path -- and then harmonise_data(action = 2) discards it in the
# caller. The discovery screen (cascade/01_mr_screen.R:681-687) excludes exactly those rows before
# computing `missing`, on the stated ground that "a direct row that harmonisation will discard is
# not a direct row". Measured consequence with this flag off: ARIC's PECAM1-CAD (rs1962147, G/C)
# and RAD23B-T2D (rs7035725, G/C) reach NO_OVERLAP_AFTER_HARMONISE while the discovery arm used a
# proxy for both -- i.e. the two arms were NOT reading the same outcome data, which is the claim
# both replication scripts make in their own headers.
local_extract_outcome_data <- function(snps, outcomes, proxies = FALSE,
                                       palindromic_fallthrough = FALSE, ...) {
  oid <- outcomes[1]; m <- .LO_META[[oid]]
  if (is.null(m)) stop("_local_outcome_extract: unknown outcome id (not in .LO_META): ", oid)
  # These are OUTCOME-space names. rs format is a property of the LD panel, not of them, and
  # requiring it here discarded instruments carrying a chr:pos-style identifier before any outcome
  # was read -- no direct match and no proxy, for a variant the outcome does contain.
  snps <- unique(snps[!is.na(snps) & nzchar(snps)]); if (!length(snps)) return(NULL)
  keep <- snps
  .lo_load_inst_al()                                   # the direct gate below needs these either way
  if (isTRUE(proxies)) { .lo_load_proxy(); keep <- unique(c(snps, .LO_prox[.(snps), proxy, nomatch = 0L])) }
  lk <- .LO_lookup(m, keep); if (is.null(lk) || !nrow(lk)) return(NULL)
  mk <- function(SNP, ea, oa, beta, se, pval, eaf, n, is_proxy = FALSE, proxy_snp = NA_character_, proxy_r2 = NA_real_) data.frame(SNP = SNP,
    beta.outcome = beta, se.outcome = se, effect_allele.outcome = ea, other_allele.outcome = oa,
    eaf.outcome = eaf, pval.outcome = pval, samplesize.outcome = n,
    ncase.outcome = if (!is.null(m$ncase)) m$ncase else NA_integer_,
    outcome = oid, id.outcome = oid, mr_keep.outcome = TRUE,
    proxy.outcome = is_proxy, proxy_snp.outcome = proxy_snp, proxy_r2.outcome = proxy_r2, stringsAsFactors = FALSE)
  # A direct match is the same VARIANT: the outcome row's unordered allele pair must be the
  # instrument's, on either strand. This gate governed only proxy eligibility and was never applied
  # to `direct` itself, so an rsID carried under a different allele pair was published under the
  # instrument's identity -- a different contrast reported as that instrument's association.
  .ia <- .LO_inst_al[.(intersect(snps, unique(lk$rs)))][!is.na(ie) & !is.na(io)]
  .ia[, .apair := paste(pmin(toupper(ie), toupper(io)), pmax(toupper(ie), toupper(io)))]
  .cm <- c(A = "T", T = "A", C = "G", G = "C")
  .rc <- .ia[grepl("^[ACGT]$", toupper(ie)) & grepl("^[ACGT]$", toupper(io))]
  if (nrow(.rc)) .rc[, .apair := paste(pmin(.cm[toupper(ie)], .cm[toupper(io)]), pmax(.cm[toupper(ie)], .cm[toupper(io)]))]
  .ds <- merge(lk, unique(rbind(.ia[, .(rs = SNP, .apair)], .rc[, .(rs = SNP, .apair)])), by = c("rs", ".apair"), sort = FALSE)
  .ds <- .ds[order(pval)][!duplicated(rs)]
  # ── verbatim from cascade/01_mr_screen.R:676-687, gated ───────────────────────────────────────
  # EITHER SIDE. harmonise_data(action=2) calls a palindrome ambiguous when the exposure OR the
  # outcome frequency is uninformative, so testing the outcome alone leaves rows that it would drop
  # classified as directly available -- and the proxy path is skipped for exactly those.
  if (isTRUE(palindromic_fallthrough) && nrow(.ds)) {
    .pe <- toupper(.ds$ea); .po <- toupper(.ds$oa)
    .pal <- .pe %chin% names(.LO_CMP) & .po %chin% names(.LO_CMP) & .po == unname(.LO_CMP[.pe])
    .ex  <- .LO_inst_al[.(.ds$rs)]$eaf_exp
    .amb <- .pal & (is.na(.ds$eaf) | (.ds$eaf > 0.42 & .ds$eaf < 0.58) |
                    is.na(.ex)     | (.ex     > 0.42 & .ex     < 0.58))
    if (any(.amb)) {
      cat(sprintf("  [%s] %d direct palindromic row(s) whose EAF cannot resolve strand -> proxy path\n",
                  basename(m$file), sum(.amb)))
      .ds <- .ds[!.amb]
    }
  }
  direct <- if (nrow(.ds)) mk(.ds$rs, .ds$ea, .ds$oa, .ds$beta, .ds$se, .ds$pval, .ds$eaf, .ds$n) else NULL
  proxied <- NULL; n_prox <- 0L
  if (isTRUE(proxies)) {
    rsq_thr <- { .dd <- list(...); if (!is.null(.dd$rsq)) as.numeric(.dd$rsq) else 0.8 }   # honour caller rsq (table is r2>=0.8)
    missing <- setdiff(snps, .ds$rs)                             # instrument SNPs with NO compatible direct match
    pr_av <- .LO_prox[instrument %in% missing & proxy %in% lk$rs & r2 >= rsq_thr]
    rows <- list()
    for (g in unique(pr_av$instrument)) {
      cands <- pr_av[instrument == g][order(-r2)]; ia <- .LO_inst_al[.(g)]
      if (is.na(ia$ie)) next
      for (j in seq_len(nrow(cands))) {
        .pos <- lk[.(cands$proxy[j])]
        # The proxy rsID can now carry several allele pairs; the orienter judges ONE record.
        # po must be the record the orienter accepted, not the first one, or the p-value comes from
        # a different allele contrast than the beta and SE beside it.
        or <- NULL; po <- .pos
        for (.k in seq_len(nrow(.pos))) { or <- .lo_orient_proxy(ia$ie, cands[j], .pos[.k], full = isTRUE(palindromic_fallthrough)); if (!is.null(or)) { po <- .pos[.k]; break } }
        if (!is.null(or)) { rows[[length(rows) + 1]] <- mk(g, ia$ie, ia$io, or$beta, or$se, po$pval, or$eaf, or$n, is_proxy = TRUE, proxy_snp = cands$proxy[j], proxy_r2 = cands$r2[j]); break }
      }
    }
    if (length(rows)) { proxied <- do.call(rbind, rows); n_prox <- nrow(proxied) }
  }
  out <- rbind(direct, proxied)
  # Codex backfill-review fix: coverage visibility (never silently drop instruments) -- report direct/proxy/missing
  cat(sprintf("    [local-extract %s] requested=%d direct=%d proxy=%d missing=%d\n",
      oid, length(snps), if (is.null(direct)) 0L else nrow(direct), n_prox,
      length(snps) - (if (is.null(out)) 0L else nrow(out))))
  if (is.null(out) || !nrow(out)) return(NULL)
  out
}
