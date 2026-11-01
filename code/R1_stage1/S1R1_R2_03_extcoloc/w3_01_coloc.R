# W3 step 1 (pre-declaration §3): external regional colocalization with R coloc 5.2.3.
# Pairs: PAM-T2D, RAD23B-T2D (DIAMANTE 2018 EUR noUKB), TGFB1-CAD, PCSK9-CAD (CARDIoGRAMplusC4D ieu-a-7).
# Exposure: UKB-PPP EUR cis window, beta/varbeta, quantitative, sdY = 1 (per SD of inverse-rank-normalised NPX).
# Outcome, primary: p values + MAF + N with the resource's case fraction (type "cc"); sensitivity: printed beta/SE.
# Join: GRCh37 chr:pos parsed from the UKB-PPP ID (the UKB-PPP `pos` column is GRCh38) + unordered allele set.
# Palindromic SNPs dropped if EITHER source's frequency is in [0.42, 0.58]. p1 = p2 = 1e-4; p12 = 1e-5 and 1e-6.
# Reads only: governing docs (md5), RESOURCES.json, the four exposure files, this folder's extracts (md5 from
# extract_manifest.json, written by w3_00_extract.py), the 489-panel chr5/chr9 .bim (rsID labels only), and the
# manuscript Methods (definition of the colocalization rule, Rule 1).
# No LD is used (coloc.abf assumes one causal variant per trait), so PLINK, the LD bootstrap and Rule 3(ii) do not apply.

# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
suppressPackageStartupMessages({ library(data.table); library(coloc); library(jsonlite) })
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_03_extcoloc")
RESDIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
LOGF <- file.path(OUT, "run.log")
t_start <- Sys.time()
log_ <- function(fmt, ...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " [coloc] ", sprintf(fmt, ...))
  cat(line, "\n", sep = ""); cat(line, "\n", sep = "", file = LOGF, append = TRUE)
}
md5 <- function(p) unname(tools::md5sum(p))
assert_md5 <- function(p, m, label) {
  got <- md5(p)
  if (is.na(got) || got != m) stop(sprintf("MD5 MISMATCH %s: %s expected %s got %s", label, p, m, got))
  log_("md5 OK %s %s %s", label, m, p)
}

GOV <- setNames(c("40d0d52f7055542d774c03f1edf04e00", "a281c3464619b1b67c28e13e7392c41c", "f265d600973c43239bb8377e0cd39d33", "3319dc39740d9f5b3b8822133d9108b7", "377508832ea9bd384133e869c456c909", "b7390fdf5b6f2f9c1d30c298d6afdafc", "6d21238155a62fb58f71e05a01638dd0"), c(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/STAGE1_BRIEFS.md"), paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/P31R_R2_STAGE1_PRESPEC_REDACTED.md"), paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/gls_q.R"), paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/RESOURCES.json")))
# Definition of the colocalization rule (Rule 1 primary), not in RESOURCES.json; registered in README.md.
METHODS <- setNames(c("75bb35c9f48255cbc4c6a604b07797ad"), c(paste0(.P31R_V3, "/R2_rewrite/sections/P31R_Methods.md")))

log_("R %s; coloc %s; data.table %s; jsonlite %s; no RNG used by this step",
     R.version.string, as.character(packageVersion("coloc")), as.character(packageVersion("data.table")),
     as.character(packageVersion("jsonlite")))
stopifnot(as.character(packageVersion("coloc")) == "5.2.3")
for (p in names(GOV)) assert_md5(p, GOV[[p]], "governing")
assert_md5(names(METHODS), METHODS[[1]], "definition: manuscript Methods (colocalization rule)")
rule_line <- grep("A pair was colocalization-supported when", readLines(names(METHODS), encoding = "UTF-8", warn = FALSE), value = TRUE)
stopifnot(length(rule_line) == 1, grepl("at least 0.80", rule_line), grepl("at least 5 \\[\\[L237\\]\\] times the posterior for distinct variants", rule_line))
log_("colocalization rule (Methods): PP.H4 >= 0.80 AND PP.H4 >= 5 x PP.H3, priors p1 = p2 = 1e-4, p12 = 1e-5")
log_("LD: none used. coloc.abf assumes one causal variant per trait, so PLINK19_paper6 is not run, no LD panel is read, and Rule 3(ii) (503/489 + 5,000-draw bootstrap) and Rule 3(iv) do not apply to §3; no RNG is used, so seed 20261007 is not needed")
log_("other files read (method/definition only, md5 at read): manuscript Methods %s (definition of the colocalization rule; not in RESOURCES.json); engine V3/V4_local/builders/_vendor/l3__mr_coloc.py fa6db396645e411504c4825864f75dd4 (RESOURCES role engine: read to see how the screen fed coloc.abf, beta/se with prior sd 0.15; no code copied); R2_rewrite/supplement/P31R_SM_candidates.md 4ecdf93d473efc9dee320a4724d20643 (DEF_SM_CANDIDATES, definition; read while locating the rule); outcome-check pre-specification c9079f2de1ffe75ce31f8a1d597d1fda (governing; instrument table and N). Seen scripts opened: none. Seen outputs opened before COMPARE_SEEN: none", names(METHODS))

# gls_q.R is not used by §3 (no GLS Q); its md5 is asserted above and its self-test is run for the log (rule 6).
if (Sys.getenv("W3_DEV") == "1") { log_("DEV RUN: gls_q.R self-test skipped (development iteration only; the final run does not set W3_DEV)") } else {
  st <- system2(file.path(R.home("bin"), "Rscript.exe"), c(shQuote(names(GOV)[6]), "--self-test"), stdout = TRUE, stderr = TRUE)
  log_("gls_q.R --self-test (not used by W3): exit %s; last line: %s", as.character(attr(st, "status") %||% 0L), tail(st, 1))
  if (!identical(tail(st, 1), "SELF-TEST PASS")) stop("gls_q.R self-test did not pass")
}

RES <- fromJSON(names(GOV)[7], simplifyVector = FALSE)
RESL <- setNames(RES, vapply(RES, function(e) e$id, ""))
res_get <- function(id) { e <- RESL[[id]]; if (is.null(e)) stop("RESOURCES id not found: ", id); e }

MAN <- fromJSON(file.path(OUT, "extract_manifest.json"), simplifyVector = FALSE)
ENG <- setNames(c("fa6db396645e411504c4825864f75dd4", "4ecdf93d473efc9dee320a4724d20643", "c9079f2de1ffe75ce31f8a1d597d1fda"), c(paste0(.P31R_V3, "/V4_local/builders/_vendor/l3__mr_coloc.py"), paste0(.P31R_V3, "/R2_rewrite/supplement/P31R_SM_candidates.md"), paste0(.P31R_V3, "/R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md")))
for (q in names(ENG)) assert_md5(q, ENG[[q]], "read for method/definition only")

PAIRS <- data.table(
  pair = c("PAM_T2D", "RAD23B_T2D", "TGFB1_CAD", "PCSK9_CAD"),
  protein = c("PAM", "RAD23B", "TGFB1", "PCSK9"),
  exp_id = c("EXP_UKBPPP_EUR_PAM_OID20291", "EXP_UKBPPP_EUR_RAD23B_OID21166", "EXP_UKBPPP_EUR_TGFB1_OID20621", "EXP_UKBPPP_EUR_PCSK9_OID20235"),
  out_id = c("OUT_DIAMANTE2018_noUKB_EUR_chrpos", "OUT_DIAMANTE2018_noUKB_EUR_chrpos", "OUT_CARDIOGRAMplusC4D_ieu-a-7", "OUT_CARDIOGRAMplusC4D_ieu-a-7"),
  bim_id = c("LD489_chr5_bim", "LD489_chr9_bim", NA, NA))
P1 <- 1e-4; P2 <- 1e-4; P12S <- c(1e-5, 1e-6)
PAL_LO <- 0.42; PAL_HI <- 0.58

# ---------------- N and case fraction from RESOURCES n_info (computed, never typed) ----------------
parse_n <- function(id) {
  ni <- res_get(id)$n_info
  m <- regmatches(ni, regexec("cases ([0-9]+) / controls ([0-9]+) \\(N ([0-9]+)\\); N_eff = [^=]+= ([0-9]+)", ni))[[1]]
  if (length(m) != 5) stop("cannot parse n_info for ", id, ": ", ni)
  ca <- as.numeric(m[2]); co <- as.numeric(m[3]); Nprint <- as.numeric(m[4]); Neff_print <- as.numeric(m[5])
  N <- ca + co; s <- ca / N; Neff <- 4 * ca * co / (ca + co)
  if (N != Nprint) stop("N != cases + controls for ", id)
  # coloc's case-control variance uses N s (1 - s), which equals N_eff / 4: same information as N_eff with s = 0.5
  stopifnot(abs(N * s * (1 - s) - Neff / 4) < 1e-6)
  log_("n_info %s: cases %d, controls %d, N %d, case fraction s = %.6f; N_eff computed %.1f (RESOURCES prints %d; diff %.1f)",
       id, ca, co, N, s, Neff, Neff_print, Neff - Neff_print)
  list(cases = ca, controls = co, N = N, s = s, Neff = Neff, Neff_print = Neff_print, n_info = ni)
}
NINFO <- list(OUT_DIAMANTE2018_noUKB_EUR_chrpos = parse_n("OUT_DIAMANTE2018_noUKB_EUR_chrpos"),
              `OUT_CARDIOGRAMplusC4D_ieu-a-7` = parse_n("OUT_CARDIOGRAMplusC4D_ieu-a-7"))

# ---------------- helpers ----------------
comp <- c(A = "T", T = "A", C = "G", G = "C")
is_pal <- function(a, b) nchar(a) == 1 & nchar(b) == 1 & ((a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C"))
aset <- function(a, b) ifelse(a < b, paste(a, b, sep = "/"), paste(b, a, sep = "/"))
flipset <- function(a, b) {
  ok <- grepl("^[ACGT]+$", a) & grepl("^[ACGT]+$", b)
  fa <- vapply(strsplit(a, ""), function(x) paste(comp[x], collapse = ""), ""); fb <- vapply(strsplit(b, ""), function(x) paste(comp[x], collapse = ""), "")
  ifelse(ok, aset(fa, fb), NA_character_)
}
null_rates <- function(dt, cols) paste(sprintf("%s=%.4f", cols, vapply(cols, function(c) mean(is.na(dt[[c]]) | dt[[c]] == ""), 0)), collapse = ", ")

read_exposure <- function(id, n_expected) {
  e <- res_get(id); stopifnot(e$role == "raw_input"); assert_md5(e$path, e$md5, id)
  d <- fread(file = e$path, sep = "\t", colClasses = list(character = c("ID", "ALLELE0", "ALLELE1")))
  if (nrow(d) != n_expected) stop(sprintf("exposure %s rows %d != %d (blind reader?)", id, nrow(d), n_expected))
  sp <- tstrsplit(d$ID, ":", fixed = TRUE)
  d[, `:=`(chr = sp[[1]], pos37 = as.integer(sp[[2]]), ref37 = sp[[3]], alt37 = sp[[4]])]
  stopifnot(all(d$ref37 == d$ALLELE0), all(d$alt37 == d$ALLELE1), all(d$chr == as.character(d$chrom)))
  for (c in c("A1FREQ", "BETA", "SE", "LOG10P")) stopifnot(is.numeric(d[[c]]), all(is.finite(d[[c]])))
  d
}

read_outcome <- function(pair, out_id) {
  x <- MAN$extracts[[pair]]
  stopifnot(x$source == out_id, MAN$source_md5[[out_id]] == res_get(out_id)$md5)
  p <- file.path(OUT, x$file); assert_md5(p, x$md5, sprintf("extract %s (stream-cut of %s)", pair, out_id))
  parse_outcome(p, out_id, x$rows)
}
parse_outcome <- function(p, out_id, n_expected) {
  if (out_id == "OUT_DIAMANTE2018_noUKB_EUR_chrpos") {
    d <- fread(file = p, sep = "\t", colClasses = list(character = c("SNP", "Chr", "EA", "NEA", "EAF", "Beta", "SE", "Pvalue")))
    o <- data.table(chr = d$Chr, pos37 = as.integer(d$Pos), ea = d$EA, nea = d$NEA, eaf_txt = d$EAF,
                    eaf = as.numeric(d$EAF), beta = as.numeric(d$Beta), se = as.numeric(d$SE), p = as.numeric(d$Pvalue),
                    p_txt = d$Pvalue, rsid = NA_character_, ss = NA_real_)
  } else {
    d <- fread(file = p, sep = "\t", colClasses = list(character = c("CHROM", "ID", "REF", "ALT")))
    stopifnot(all(d$multiallelic == 0))
    o <- data.table(chr = d$CHROM, pos37 = as.integer(d$POS), ea = d$ALT, nea = d$REF, eaf_txt = as.character(d$AF),
                    eaf = as.numeric(d$AF), beta = as.numeric(d$ES), se = as.numeric(d$SE), p = 10^(-as.numeric(d$LP)),
                    p_txt = paste0("LP=", d$LP), rsid = d$ID, ss = as.numeric(d$SS))
  }
  if (nrow(o) != n_expected) stop(sprintf("extract %s rows %d != manifest %d (blind reader?)", p, nrow(o), n_expected))
  stopifnot(!anyNA(o[, .(chr, pos37, ea, nea, eaf, beta, se, p)]))
  o
}

# ---------------- join + QC for one pair ----------------
# MAF for the outcome p-value input (decision recorded before seeing its results; see README "MAF handling"):
#   the outcome's own printed EAF; where DIAMANTE prints EAF as 0 or 1 (MAF < 0.005 not representable at 2 s.f.),
#   the UKB-PPP European A1FREQ of the same variant is substituted ("maf_src = exposure"). The variant set therefore
#   stays the declared UKB-PPP cis window. A descriptive variant that drops those SNPs instead ("pmafdrop") is printed.
harmonise_pair <- function(pr, mutate_shift = 0L, mutate_chr = FALSE, quiet = FALSE) {
  w <- MAN$windows[[pr$pair]]
  ex <- read_exposure(pr$exp_id, w$n_exposure_rows)
  if (mutate_shift != 0L) ex[, pos37 := pos37 + mutate_shift]
  if (mutate_chr) ex[, chr := paste0("chr", chr)]
  oc <- read_outcome(pr$pair, pr$out_id)
  ex[, key := paste(chr, pos37, aset(ALLELE0, ALLELE1), sep = ":")]
  oc[, `:=`(key = paste(chr, pos37, aset(ea, nea), sep = ":"), keyflip = paste(chr, pos37, flipset(ea, nea), sep = ":"))]
  n_out_bulk <- oc[pos37 >= w$start37 & pos37 <= w$end37, .N]
  # a key carried by two exposure rows (e.g. G:GGGTT and GGGTT:G at one position) cannot be matched unambiguously: drop both
  dupk_ex <- unique(ex$key[duplicated(ex$key)]); dupk_oc <- unique(oc$key[duplicated(oc$key)])
  n_dup_ex_rows <- sum(ex$key %in% dupk_ex); n_dup_oc_rows <- sum(oc$key %in% dupk_oc)
  ex_u <- ex[!key %in% dupk_ex]; oc_u <- oc[!key %in% dupk_oc]
  pos_match <- length(intersect(unique(paste(ex$chr, ex$pos37)), unique(paste(oc$chr, oc$pos37))))
  # exact allele-set match first; then strand-complement match for non-palindromic variants only
  m1 <- merge(ex_u, oc_u[, !"keyflip"], by = "key", suffixes = c(".x", ".y"))
  m1[, match_type := rep("exact", .N)]
  rest_ex <- ex_u[!key %in% m1$key & !is_pal(ALLELE0, ALLELE1)]
  oc2 <- oc_u[!key %in% m1$key & !is.na(keyflip) & !is_pal(ea, nea)]
  oc2[, key := keyflip]
  m2 <- merge(rest_ex, oc2[, !"keyflip"], by = "key", suffixes = c(".x", ".y"))
  m2[, match_type := rep("strand_complement", .N)]
  m <- rbind(m1, m2, fill = TRUE)
  if (nrow(m) == 0) stop(sprintf("JOIN YIELD 0 for %s: a join yielding 0 is a bug, not a result", pr$pair))
  stopifnot(!anyDuplicated(m$ID), !anyDuplicated(m$key))
  # outcome effect allele on the exposure strand (for the beta/SE sensitivity; coloc.abf uses beta^2, so sign is immaterial)
  m[, ea_on_exp := ifelse(match_type == "exact", ea, vapply(strsplit(ea, ""), function(x) paste(comp[x], collapse = ""), ""))]
  stopifnot(all(m$ea_on_exp == m$ALLELE1 | m$ea_on_exp == m$ALLELE0))
  m[, `:=`(beta_out_al = ifelse(ea_on_exp == ALLELE1, beta, -beta), eaf_out_al = ifelse(ea_on_exp == ALLELE1, eaf, 1 - eaf))]
  m[, pal := is_pal(ALLELE0, ALLELE1)]
  m[, pal_drop := pal & ((A1FREQ >= PAL_LO & A1FREQ <= PAL_HI) | (eaf >= PAL_LO & eaf <= PAL_HI))]
  m[, maf_exp := pmin(A1FREQ, 1 - A1FREQ)]
  m[, maf_out_own := pmin(eaf, 1 - eaf)]
  m[, maf0 := !pal_drop & maf_out_own <= 0]
  m[, maf_out := ifelse(maf0, maf_exp, maf_out_own)]
  m[, maf_src := ifelse(maf0, "exposure", "outcome")]
  m[, bad_drop := !pal_drop & (!is.finite(p) | p <= 0 | p > 1 | se <= 0 | SE <= 0 | maf_out <= 0 | maf_exp <= 0)]
  keep <- m[!pal_drop & !bad_drop]
  setorder(keep, pos37.x)
  yield_exp <- nrow(m) / nrow(ex); yield_out <- nrow(m) / n_out_bulk; yield_min <- nrow(m) / min(nrow(ex), n_out_bulk)
  sub <- keep[maf_src == "exposure"]
  j <- list(pair = pr$pair, n_exp = nrow(ex), n_out_window = nrow(oc), n_out_bulk = n_out_bulk, pos_match = pos_match,
            matched = nrow(m), matched_exact = nrow(m1), matched_complement = nrow(m2),
            yield_exp = yield_exp, yield_out = yield_out, yield_min = yield_min,
            dup_key_rows_exp = n_dup_ex_rows, dup_key_rows_out = n_dup_oc_rows, n_pal = sum(m$pal), pal_dropped = sum(m$pal_drop),
            maf_substituted = nrow(sub), bad_dropped = sum(m$bad_drop), nsnp_in = nrow(keep),
            nulls_exp = null_rates(ex, c("chr", "pos37", "ALLELE0", "ALLELE1")), nulls_out = null_rates(oc, c("chr", "pos37", "ea", "nea")),
            sub_maf_exp_max = if (nrow(sub)) max(sub$maf_exp) else NA_real_,
            sub_maf_exp_gt_0.01 = if (nrow(sub)) sum(sub$maf_exp > 0.01) else 0L,
            sub_maxlp_exp = if (nrow(sub)) max(sub$LOG10P) else NA_real_,
            sub_minp_out = if (nrow(sub)) min(sub$p) else NA_real_,
            maxlp_exp_pal = if (any(m$pal_drop)) max(m$LOG10P[m$pal_drop]) else NA_real_,
            minp_out_pal = if (any(m$pal_drop)) min(m$p[m$pal_drop]) else NA_real_,
            ss_range = if (all(is.na(keep$ss))) "no per-variant N column" else paste(range(keep$ss), collapse = "-"))
  if (!quiet) log_(paste0("join %s: exposure %d rows, outcome window %d rows (%d in the GRCh37 bulk window); positions shared %d; ",
                          "ambiguous duplicate-key rows dropped exp %d / out %d; allele-matched %d (exact %d, strand-complement %d); ",
                          "yield %.3f of exposure, %.3f of outcome bulk, %.3f of min; per-key null rates exposure [%s] outcome [%s]; ",
                          "palindromic %d, dropped by the 0.42-0.58 rule %d (max exposure -log10p among them %.2f, min outcome p %s); ",
                          "outcome EAF printed 0 or 1 -> UKB-PPP MAF substituted on %d SNPs (their UKB-PPP MAF max %.4f, %d above 0.01; ",
                          "max exposure -log10p %.2f; min outcome p %s); other invalid dropped %d; nsnp into coloc %d; outcome per-variant N %s"),
                   pr$pair, j$n_exp, j$n_out_window, j$n_out_bulk, j$pos_match, j$dup_key_rows_exp, j$dup_key_rows_out,
                   j$matched, j$matched_exact, j$matched_complement, j$yield_exp, j$yield_out, j$yield_min, j$nulls_exp, j$nulls_out,
                   j$n_pal, j$pal_dropped, as.numeric(j$maxlp_exp_pal), format(j$minp_out_pal), j$maf_substituted,
                   as.numeric(j$sub_maf_exp_max), j$sub_maf_exp_gt_0.01, as.numeric(j$sub_maxlp_exp), format(j$sub_minp_out),
                   j$bad_dropped, j$nsnp_in, j$ss_range)
  list(keep = keep, all = m, join = j)
}
YIELD_MIN <- 0.5
check_yield <- function(j) if (j$yield_min < YIELD_MIN) stop(sprintf("join yield %s %.3f < %.2f of min(exposure, outcome bulk)", j$pair, j$yield_min, YIELD_MIN))

run_coloc <- function(keep, ninfo, input, p12) {
  if (input == "pmafdrop") keep <- keep[maf_src == "outcome"]
  D1 <- list(snp = keep$ID, beta = keep$BETA, varbeta = keep$SE^2, type = "quant", sdY = 1, MAF = keep$maf_exp)
  if (input %in% c("pmaf", "pmafdrop")) {
    D2 <- list(snp = keep$ID, pvalues = keep$p, MAF = keep$maf_out, N = ninfo$N, s = ninfo$s, type = "cc")
  } else {
    D2 <- list(snp = keep$ID, beta = keep$beta_out_al, varbeta = keep$se^2, type = "cc", s = ninfo$s, N = ninfo$N, MAF = keep$maf_out)
  }
  # coloc's informational "minimum p value is" warnings are logged, not printed to stderr
  wmsg <- character(0)
  r <- withCallingHandlers({
    check_dataset(D1); check_dataset(D2)
    out <- NULL; invisible(capture.output(out <- coloc.abf(D1, D2, p1 = P1, p2 = P2, p12 = p12))); out
  }, warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") })
  attr(r, "warnings") <- unique(sub("\n.*", "", wmsg))
  attr(r, "nsnp_expected") <- nrow(keep)
  r
}

# ---------------- branch logic (§3 wording rule; PCSK9 control criterion) ----------------
rule_met <- function(h4, h3) h4 >= 0.80 & h4 >= 5 * h3
s3_branch <- function(h4_5, h3_5, h4_6, h3_6) {
  m5 <- rule_met(h4_5, h3_5); m6 <- rule_met(h4_6, h3_6)
  if (m5 && m6) return(list(branch = 1L, addon = FALSE))
  if (m5 && !m6) return(list(branch = 2L, addon = h3_6 >= 0.5))
  if (!m5 && (h3_5 >= 0.5 || h3_6 >= 0.5)) return(list(branch = 3L, addon = FALSE))
  list(branch = 4L, addon = FALSE)
}
pcsk9_control <- function(b) {
  if (b$branch == 1L || (b$branch == 2L && !b$addon)) return("behaved as a control")
  if (b$branch == 3L) return("did not behave as a clean control: favoured distinct variants")
  "did not behave as a clean control: did not resolve a shared variant"
}
# unit tests (boundaries included); each case names the branch it must produce
ut <- list(list(c(0.90, 0.05, 0.85, 0.10), 1L, FALSE), list(c(0.80, 0.16, 0.80, 0.16), 1L, FALSE),
           list(c(0.80, 0.161, 0.80, 0.16), 4L, FALSE),   # ratio 4.97 < 5 at 1e-5 -> rule not met; PP.H3 < 0.5 at both -> branch 4
           list(c(0.90, 0.05, 0.40, 0.55), 2L, TRUE), list(c(0.90, 0.05, 0.60, 0.30), 2L, FALSE),
           list(c(0.30, 0.60, 0.04, 0.90), 3L, FALSE), list(c(0.40, 0.30, 0.10, 0.50), 3L, FALSE),
           list(c(0.40, 0.30, 0.10, 0.49), 4L, FALSE), list(c(0.79, 0.01, 0.30, 0.01), 4L, FALSE))
for (k in seq_along(ut)) {
  v <- ut[[k]][[1]]; b <- s3_branch(v[1], v[2], v[3], v[4])
  if (b$branch != ut[[k]][[2]] || b$addon != ut[[k]][[3]]) stop(sprintf("branch unit test %d failed: got %d/%s", k, b$branch, b$addon))
}
stopifnot(pcsk9_control(list(branch = 1L, addon = FALSE)) == "behaved as a control",
          pcsk9_control(list(branch = 2L, addon = TRUE)) != "behaved as a control",
          pcsk9_control(list(branch = 2L, addon = FALSE)) == "behaved as a control",
          grepl("distinct", pcsk9_control(list(branch = 3L, addon = FALSE))),
          grepl("resolve", pcsk9_control(list(branch = 4L, addon = FALSE))))
log_("branch-function unit tests: %d of %d passed (boundaries PP.H4 = 0.80, ratio 5, PP.H3 = 0.5 included); PCSK9 mapping 5/5", length(ut), length(ut))

# ---------------- main runs ----------------
results <- list(); branches <- list(); asserts <- list()
addr <- function(qid, value, notes = "", resource = "", estimand = "", n = NA, unit = "", p = NA, prespec_ref = "§3") {
  results[[length(results) + 1]] <<- data.table(quantity_id = qid, prespec_ref = prespec_ref, estimand = estimand, resource = resource,
    panel = "none (coloc.abf uses no LD)", cov = "NA", se_order = "NA", value = as.character(value), se = NA, ci_lo = NA, ci_hi = NA,
    p = p, n = n, unit = unit, notes = notes)
}
fmtpp <- function(x) trimws(formatC(x, digits = 6, format = "g"))
rsid_lookup <- function(bim_id, q_chr, q_pos, q_a0, q_a1) {
  if (is.na(bim_id)) return(NA_character_)
  e <- res_get(bim_id); if (is.null(BIMS[[bim_id]])) { assert_md5(e$path, e$md5, bim_id)
    BIMS[[bim_id]] <<- fread(file = e$path, header = FALSE, col.names = c("bchr", "rsid", "cm", "bpos", "ba1", "ba2"), colClasses = list(character = c(1, 2, 5, 6)))
    stopifnot(nrow(BIMS[[bim_id]]) > 100000) }
  qs <- aset(q_a0, q_a1); qf <- flipset(q_a0, q_a1)
  b <- BIMS[[bim_id]][bchr == q_chr & bpos == q_pos]
  b <- b[aset(ba1, ba2) == qs | aset(ba1, ba2) == qf]
  if (nrow(b) == 1) b$rsid else if (nrow(b) == 0) "not in 489 bim" else paste(b$rsid, collapse = "|")
}
BIMS <- list()

ALLRUNS <- list(); JOINS <- list()
INPUT_LAB <- c(pmaf = "p + MAF + N [primary]", betase = "printed beta/SE [sensitivity]",
               pmafdrop = "p + MAF + N, SNPs with outcome EAF printed 0/1 dropped instead of MAF-substituted [descriptive]")
for (i in seq_len(nrow(PAIRS))) {
  pr <- PAIRS[i]
  ninfo <- NINFO[[pr$out_id]]
  t0 <- Sys.time()
  h <- harmonise_pair(pr)
  JOINS[[pr$pair]] <- h$join
  keep <- h$keep; j <- h$join
  # the analysed variant set (summary statistics only), for the verifier and w3_02_selfcheck.py
  fwrite(keep[, .(ID, key, chr = chr.x, pos37 = pos37.x, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P, match_type, ea, nea, eaf_txt, eaf,
                  beta_out = beta, se_out = se, p_out = p, p_txt, beta_out_al, maf_exp, maf_out, maf_src, pal, rsid)],
         file.path(OUT, sprintf("coloc_input_%s.tsv", pr$pair)), sep = "\t")
  res_lab <- if (pr$out_id == "OUT_DIAMANTE2018_noUKB_EUR_chrpos") "DIAMANTE 2018 EUR noUKB (OUT_DIAMANTE2018_noUKB_EUR_chrpos)" else "CARDIoGRAMplusC4D 2015 (OUT_CARDIOGRAMplusC4D_ieu-a-7)"
  ychk <- tryCatch({ check_yield(j); "pass" }, error = function(e) conditionMessage(e))
  if (ychk != "pass") {
    asserts[[length(asserts) + 1]] <- sprintf("ASSERT_FAILED S3.%s.join.yield_min: %s", pr$pair, ychk)
    log_("ASSERT_FAILED %s (run continues)", ychk)
  } else log_("join-yield check %s: %.3f >= %.2f (pass)", pr$pair, j$yield_min, YIELD_MIN)
  # exposure / outcome leads in the analysed set
  le <- keep[which.max(LOG10P)]; lo <- keep[which.min(p)]
  le_rs <- if (!is.na(pr$bim_id)) rsid_lookup(pr$bim_id, le$chr.x, le$pos37.x, le$ALLELE0, le$ALLELE1) else le$rsid
  lo_rs <- if (!is.na(pr$bim_id)) rsid_lookup(pr$bim_id, lo$chr.x, lo$pos37.x, lo$ALLELE0, lo$ALLELE1) else lo$rsid
  jn <- sprintf("exposure rows %d; outcome window rows %d (GRCh37 bulk %d); ambiguous duplicate-key exposure rows dropped %d; allele-matched %d (exact %d, complement %d); palindromic dropped %d; outcome MAF substituted from UKB-PPP %d; invalid dropped %d",
                j$n_exp, j$n_out_window, j$n_out_bulk, j$dup_key_rows_exp, j$matched, j$matched_exact, j$matched_complement, j$pal_dropped, j$maf_substituted, j$bad_dropped)
  addr(sprintf("S3.%s.join.matched", pr$pair), j$matched, jn, res_lab, "join on GRCh37 chr:pos (UKB-PPP ID) + allele set", unit = "SNPs")
  addr(sprintf("S3.%s.join.yield_exp", pr$pair), sprintf("%.4f", j$yield_exp), "allele-matched / exposure rows", res_lab, "join yield", unit = "fraction")
  addr(sprintf("S3.%s.join.yield_out", pr$pair), sprintf("%.4f", j$yield_out), "allele-matched / outcome rows in the GRCh37 bulk window", res_lab, "join yield", unit = "fraction")
  addr(sprintf("S3.%s.join.pal_dropped", pr$pair), j$pal_dropped, sprintf("palindromic SNPs with either source's frequency in [0.42, 0.58]; of %d palindromic matched; max exposure -log10p among dropped %s; min outcome p among dropped %s", j$n_pal, format(j$maxlp_exp_pal), format(j$minp_out_pal)), res_lab, "palindrome rule (§3)", unit = "SNPs")
  addr(sprintf("S3.%s.join.maf_substituted", pr$pair), j$maf_substituted, sprintf("outcome EAF printed as 0 or 1 (MAF < 0.005 not representable at 2 s.f.): UKB-PPP MAF substituted in the p+MAF+N input (their UKB-PPP MAF max %s, %d above 0.01; max exposure -log10p %s; min outcome p %s); the pmafdrop rows drop them instead", format(j$sub_maf_exp_max), j$sub_maf_exp_gt_0.01, format(j$sub_maxlp_exp), format(j$sub_minp_out)), res_lab, "MAF handling", unit = "SNPs")
  addr(sprintf("S3.%s.lead_exposure", pr$pair), le$ID, sprintf("rsid %s; -log10p %.2f; A1FREQ %.4f; outcome p at this SNP %s", le_rs, le$LOG10P, le$A1FREQ, format(le$p)), res_lab, "exposure lead in the analysed set")
  addr(sprintf("S3.%s.lead_outcome", pr$pair), lo$ID, sprintf("rsid %s; outcome p %s (printed %s); exposure -log10p at this SNP %.2f", lo_rs, format(lo$p), lo$p_txt, lo$LOG10P), res_lab, "outcome lead in the analysed set")
  nexp <- as.numeric(sub(".*samplesize ([0-9]+) at the instrument.*", "\\1", res_get(pr$exp_id)$n_info))
  sdy_est <- coloc:::sdY.est(keep$SE^2, keep$maf_exp, nexp)
  log_("scale check %s: coloc sdY.est from varbeta, MAF and N = %d (n_info samplesize at the instrument) = %.4f; sdY = 1 is used (BETA per SD of inverse-rank-normalised NPX)", pr$pair, nexp, sdy_est)
  addr(sprintf("S3.%s.exposure_sdY_est", pr$pair), sprintf("%.4f", sdy_est), sprintf("coloc:::sdY.est(varbeta, MAF, N = %d); diagnostic only: the runs use sdY = 1", nexp), sprintf("UKB-PPP EUR (%s)", pr$exp_id), "exposure trait SD implied by the summary statistics", n = nexp, unit = "SD")
  pp <- list()
  inputs <- c("pmaf", "betase", if (j$maf_substituted > 0) "pmafdrop")
  for (input in inputs) for (p12 in P12S) {
    r <- run_coloc(keep, ninfo, input, p12)
    s <- r$summary
    stopifnot(abs(sum(s[2:6]) - 1) < 1e-9, s[["nsnps"]] == attr(r, "nsnp_expected"), all(is.finite(s)))
    top <- r$results[which.max(r$results$SNP.PP.H4), ]
    tk <- keep[ID == top$snp]
    stopifnot(nrow(tk) == 1)
    top_rs <- if (!is.na(pr$bim_id)) rsid_lookup(pr$bim_id, tk$chr.x, tk$pos37.x, tk$ALLELE0, tk$ALLELE1) else tk$rsid
    p12s <- formatC(p12, format = "e", digits = 0)
    est <- sprintf("coloc.abf %s; p1 = p2 = 1e-4, p12 = %s; exposure UKB-PPP beta/varbeta (quant, sdY = 1); outcome %s (cc, N %d, s %.5f)", pr$pair, p12s,
                   INPUT_LAB[[input]], ninfo$N, ninfo$s)
    for (hh in c("H0", "H1", "H2", "H3", "H4"))
      addr(sprintf("S3.%s.%s.%s.%s", pr$pair, input, p12s, hh), fmtpp(s[[paste0("PP.", hh, ".abf")]]), "", res_lab, est, n = ninfo$N, unit = "posterior probability")
    addr(sprintf("S3.%s.%s.%s.nsnp", pr$pair, input, p12s), s[["nsnps"]], if (input == "pmafdrop") paste(jn, "; minus the MAF-substituted SNPs") else jn, res_lab, est, n = ninfo$N, unit = "SNPs")
    addr(sprintf("S3.%s.%s.%s.top", pr$pair, input, p12s), top$snp,
         sprintf("max SNP.PP.H4 = %s; rsid %s; exposure -log10p %.2f; outcome p %s; outcome MAF source %s", fmtpp(top$SNP.PP.H4), top_rs, tk$LOG10P, format(tk$p), tk$maf_src), res_lab, est, n = ninfo$N, unit = "variant (GRCh37 chr:pos:ref:alt)")
    met <- rule_met(s[["PP.H4.abf"]], s[["PP.H3.abf"]])
    addr(sprintf("S3.%s.%s.%s.rule_met", pr$pair, input, p12s), as.integer(met),
         sprintf("colocalization rule (Methods): PP.H4 >= 0.80 and PP.H4 >= 5 x PP.H3; PP.H4/PP.H3 = %s", fmtpp(s[["PP.H4.abf"]] / s[["PP.H3.abf"]])), res_lab, est, n = ninfo$N, unit = "1 = met")
    pp[[paste(input, p12s)]] <- s
    ALLRUNS[[length(ALLRUNS) + 1]] <- data.table(pair = pr$pair, input = input, p12 = p12s, nsnp = s[["nsnps"]],
      H0 = s[["PP.H0.abf"]], H1 = s[["PP.H1.abf"]], H2 = s[["PP.H2.abf"]], H3 = s[["PP.H3.abf"]], H4 = s[["PP.H4.abf"]],
      rule_met = met, top = top$snp, top_rsid = top_rs, top_snppph4 = top$SNP.PP.H4)
    log_("coloc %s %s p12=%s: nsnp %d; PP.H0 %s H1 %s H2 %s H3 %s H4 %s; rule met %s; top %s (%s) SNP.PP.H4 %s; coloc messages: %s",
         pr$pair, input, p12s, s[["nsnps"]], fmtpp(s[["PP.H0.abf"]]), fmtpp(s[["PP.H1.abf"]]), fmtpp(s[["PP.H2.abf"]]), fmtpp(s[["PP.H3.abf"]]),
         fmtpp(s[["PP.H4.abf"]]), met, top$snp, top_rs, fmtpp(top$SNP.PP.H4), paste(attr(r, "warnings"), collapse = " | "))
  }
  # branches: the primary input decides; the other inputs are printed beside it
  bb <- list()
  for (input in inputs) {
    a <- pp[[paste(input, "1e-05")]]; b <- pp[[paste(input, "1e-06")]]
    bb[[input]] <- s3_branch(a[["PP.H4.abf"]], a[["PP.H3.abf"]], b[["PP.H4.abf"]], b[["PP.H3.abf"]])
  }
  brl <- function(x) sprintf("%d%s", x$branch, if (x$branch == 2L) if (x$addon) " (with add-on)" else " (no add-on)" else "")
  a <- pp[["pmaf 1e-05"]]; b <- pp[["pmaf 1e-06"]]
  inputs_txt <- sprintf("primary input (p + MAF + N): PP.H4 %s / PP.H3 %s at p12 = 1e-5; PP.H4 %s / PP.H3 %s at p12 = 1e-6; nsnp %d",
                        fmtpp(a[["PP.H4.abf"]]), fmtpp(a[["PP.H3.abf"]]), fmtpp(b[["PP.H4.abf"]]), fmtpp(b[["PP.H3.abf"]]), a[["nsnps"]])
  pn <- sub("_", "–", pr$pair)
  v2 <- function(x) sprintf("%.2f", x)
  tmpl <- switch(as.character(bb$pmaf$branch),
    "1" = sprintf("%s: supported in external outcome data, robust to the prior (PP.H4 %s at p12 = 10⁻⁵; %s at 10⁻⁶)", pn, v2(a[["PP.H4.abf"]]), v2(b[["PP.H4.abf"]])),
    "2" = paste0(sprintf("%s: supported at the default prior (PP.H4 %s at p12 = 10⁻⁵; %s at 10⁻⁶)", pn, v2(a[["PP.H4.abf"]]), v2(b[["PP.H4.abf"]])),
                 if (bb$pmaf$addon) sprintf("; at the conservative prior the data favour distinct variants (PP.H3 %s at 10⁻⁶)", v2(b[["PP.H3.abf"]])) else ""),
    "3" = sprintf("%s: the external data favour distinct variants (PP.H3 %s at p12 = 10⁻⁵; %s at 10⁻⁶; PP.H4 %s at 10⁻⁵)", pn, v2(a[["PP.H3.abf"]]), v2(b[["PP.H3.abf"]]), v2(a[["PP.H4.abf"]])),
    "4" = sprintf("%s: the external data do not resolve a shared variant (PP.H4 %s, PP.H3 %s at p12 = 10⁻⁵; PP.H4 %s, PP.H3 %s at 10⁻⁶)", pn, v2(a[["PP.H4.abf"]]), v2(a[["PP.H3.abf"]]), v2(b[["PP.H4.abf"]]), v2(b[["PP.H3.abf"]])))
  qual <- c(sprintf("printed beta/SE sensitivity gives branch %s", brl(bb$betase)),
            if (!is.null(bb$pmafdrop)) sprintf("descriptive MAF-drop variant gives branch %s", brl(bb$pmafdrop)),
            if (bb$pmaf$branch == 2L && bb$pmaf$addon) "the §3 conservative-prior add-on travels with the claim to every surface (Principle 0)",
            if (pr$protein %in% c("TGFB1", "PCSK9")) "single-SNP estimates published before (Wang 2025, Fig 2); the regional test is the new element (§3 prior publication)",
            if (pr$out_id != "OUT_DIAMANTE2018_noUKB_EUR_chrpos") "CARDIoGRAMplusC4D 2015 is 77% European / 13% South Asian / 6% East Asian (RESOURCES notes)",
            "Rule 2 named exception: §3 describes the pair, never changes its tier; under branch 3 the pair leaves the follow-up and candidate sentences",
            "Rule 4 'nominal' not applicable (§3 is not a heterogeneity test in the family); Rule 3(ii) not applicable (no LD, not conditional)",
            "corroboration: left to W9")
  branches[[length(branches) + 1]] <- data.table(rule_ref = "§3 wording rule (branches 1-4); Rule 1 colocalization rule = PP.H4 >= 0.80 and PP.H4 >= 5 x PP.H3",
    inputs = paste(pr$pair, inputs_txt), branch = brl(bb$pmaf), filled_template = tmpl, qualifiers = paste(qual, collapse = "; "))
  for (input in setdiff(inputs, "pmaf")) {
    a2 <- pp[[paste(input, "1e-05")]]; b2 <- pp[[paste(input, "1e-06")]]
    branches[[length(branches) + 1]] <- data.table(rule_ref = sprintf("§3 wording rule applied to the %s input, descriptive only", input),
      inputs = sprintf("%s %s: PP.H4 %s / PP.H3 %s at 1e-5; PP.H4 %s / PP.H3 %s at 1e-6; nsnp %d", pr$pair, input, fmtpp(a2[["PP.H4.abf"]]), fmtpp(a2[["PP.H3.abf"]]), fmtpp(b2[["PP.H4.abf"]]), fmtpp(b2[["PP.H3.abf"]]), a2[["nsnps"]]),
      branch = brl(bb[[input]]), filled_template = "not a wording input (the primary input decides)", qualifiers = INPUT_LAB[[input]])
  }
  if (pr$protein == "PCSK9") {
    ctl <- pcsk9_control(bb$pmaf)
    ctl_t <- if (ctl == "behaved as a control") sprintf("PCSK9–CAD behaved as a control (branch %s)", brl(bb$pmaf)) else if (bb$pmaf$branch == 3L)
      "PCSK9 did not behave as a clean control: the regional test favoured distinct variants at a known drug target [stated as a calibration of the method, with what it means for the other pairs' regional results: Stage 2 text]" else
      "PCSK9 did not behave as a clean control: the regional test did not resolve a shared variant at a known drug target [with the same statement of what this means for the other pairs: Stage 2 text]"
    branches[[length(branches) + 1]] <- data.table(rule_ref = "§3 PCSK9 control criterion (matches UMOD's; the slot stays under Principle 0)",
      inputs = sprintf("PCSK9_CAD primary-input branch %s", brl(bb$pmaf)), branch = ctl, filled_template = ctl_t,
      qualifiers = paste0("consequences written after the East Asian PP.H3 was seen (§3); printed beta/SE sensitivity would give: ", pcsk9_control(bb$betase)))
  }
  log_("pair %s done in %.1f s: primary branch %s; %s", pr$pair, as.numeric(difftime(Sys.time(), t0, units = "secs")),
       brl(bb$pmaf), paste(sprintf("%s branch %s", setdiff(inputs, "pmaf"), vapply(setdiff(inputs, "pmaf"), function(x) brl(bb[[x]]), "")), collapse = "; "))
}

# N_eff rows
for (oid in names(NINFO)) {
  ni <- NINFO[[oid]]
  addr(sprintf("S3.N_eff.%s", oid), sprintf("%.1f", ni$Neff), sprintf("computed 4*Nca*Nco/(Nca+Nco) from RESOURCES n_info (cases %d, controls %d); RESOURCES prints %d; coloc input uses N %d with s %.6f (N s (1-s) = N_eff/4)", ni$cases, ni$controls, ni$Neff_print, ni$N, ni$s),
       oid, "effective sample size", n = ni$N, unit = "individuals")
}

# ---------------- mutation tests: each check must trip on its own defect ----------------
mut <- list()
tr <- function(name, expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  ok <- !is.na(msg) && grepl(pattern, msg)
  mut[[length(mut) + 1]] <<- data.table(test = name, expected_error = pattern, error = ifelse(is.na(msg), "", substr(gsub("[[:cntrl:]]", " ", msg), 1, 200)), tripped = ok)
  log_("MUTATION %s: %s (error: %s)", name, if (ok) "tripped the named check" else "DID NOT TRIP THE NAMED CHECK", msg)
}
tr("md5 assert on a wrong expected md5", assert_md5(names(GOV)[1], "00000000000000000000000000000000", "mutation"), "MD5 MISMATCH")
tr("join-yield threshold check when exposure GRCh37 positions are shifted by +1", check_yield(harmonise_pair(PAIRS[1], mutate_shift = 1L, quiet = TRUE)$join), "^join yield PAM_T2D")
tr("join-yield-0 hard stop when the exposure chromosome is renamed (chr5 -> chrchr5)", harmonise_pair(PAIRS[1], mutate_chr = TRUE, quiet = TRUE), "JOIN YIELD 0")
tr("join-yield threshold check on a copy of the real PAM join audit with yield_min set to 0.30", check_yield(modifyList(JOINS$PAM_T2D, list(yield_min = 0.30))), "^join yield PAM_T2D 0.300")
stopifnot(YIELD_MIN == 0.5, is.null(check_yield(JOINS$PAM_T2D)))
tr("N parse assert when n_info N is altered", { RESL$OUT_DIAMANTE2018_noUKB_EUR_chrpos$n_info <- sub("N 455313", "N 455314", RESL$OUT_DIAMANTE2018_noUKB_EUR_chrpos$n_info); parse_n("OUT_DIAMANTE2018_noUKB_EUR_chrpos") }, "N != cases [+] controls")
RESL <- setNames(RES, vapply(RES, function(e) e$id, ""))  # restore the unmutated map
stopifnot(grepl("N 455313", RESL$OUT_DIAMANTE2018_noUKB_EUR_chrpos$n_info))
tr("branch unit test with a wrong expectation", { v <- c(0.90, 0.05, 0.85, 0.10); if (s3_branch(v[1], v[2], v[3], v[4])$branch != 2L) stop("branch mismatch") }, "branch mismatch")
# palindrome rule: synthetic SNPs, each must land on the named side of the rule
pal_test <- data.table(A0 = c("A", "C", "A", "G", "A"), A1 = c("T", "G", "T", "C", "G"), fe = c(0.50, 0.30, 0.30, 0.42, 0.50), fo = c(0.30, 0.30, 0.58, 0.20, 0.50))
pal_test[, drop := is_pal(A0, A1) & ((fe >= PAL_LO & fe <= PAL_HI) | (fo >= PAL_LO & fo <= PAL_HI))]
stopifnot(identical(pal_test$drop, c(TRUE, FALSE, TRUE, TRUE, FALSE)))
tr("palindrome rule with a wrong expectation", stopifnot(identical(pal_test$drop, c(FALSE, FALSE, TRUE, TRUE, FALSE))), "pal_test")
log_("palindrome-rule synthetic test: 5/5 on the expected side (exposure-band, outcome-band, boundary 0.42 and 0.58, non-palindromic at 0.50 kept)")
# coloc.abf sign invariance (justifies not orienting palindromic/complement SNPs for the coloc itself)
h <- harmonise_pair(PAIRS[3], quiet = TRUE); k2 <- copy(h$keep); r0 <- run_coloc(k2, NINFO[[PAIRS$out_id[3]]], "betase", 1e-5)
k2[, beta_out_al := -beta_out_al]; r1 <- run_coloc(k2, NINFO[[PAIRS$out_id[3]]], "betase", 1e-5)
dsign <- max(abs(r0$summary[2:6] - r1$summary[2:6]))
log_("coloc.abf sign invariance (TGFB1_CAD, beta/SE input, all outcome signs flipped): max |dPP| = %.3g", dsign)
stopifnot(dsign < 1e-12)
# placeholder guard in the extract reader: a copy of an extract with one 'NA' token must be refused
tmpd <- tempfile("w3mut_"); dir.create(tmpd)
x <- readLines(file.path(OUT, MAN$extracts$PAM_T2D$file)); x[2] <- sub("\t[^\t]+$", "\tNA", x[2]); writeLines(x, file.path(tmpd, "m.tsv"))
tr("non-numeric outcome p token refused by the real outcome parser (copy of the PAM extract, one Pvalue set to NA)",
   suppressWarnings(parse_outcome(file.path(tmpd, "m.tsv"), "OUT_DIAMANTE2018_noUKB_EUR_chrpos", MAN$extracts$PAM_T2D$rows)), "anyNA")
stopifnot(nrow(parse_outcome(file.path(OUT, MAN$extracts$PAM_T2D$file), "OUT_DIAMANTE2018_noUKB_EUR_chrpos", MAN$extracts$PAM_T2D$rows)) == MAN$extracts$PAM_T2D$rows)
unlink(tmpd, recursive = TRUE)
MUT <- rbindlist(mut)
fwrite(MUT, file.path(OUT, "mutation_tests.tsv"), sep = "\t")
if (!all(MUT$tripped)) stop("a mutation did not trip its check")

# ---------------- write outputs ----------------
R <- rbindlist(results, fill = TRUE)
for (a in asserts) R <- rbind(R, data.table(quantity_id = "ASSERT_FAILED", prespec_ref = "§3", estimand = "", resource = "", panel = "", cov = "", se_order = "",
  value = "", se = NA, ci_lo = NA, ci_hi = NA, p = NA, n = NA, unit = "", notes = a), fill = TRUE)
stopifnot(!anyDuplicated(R$quantity_id[R$quantity_id != "ASSERT_FAILED"]))
need <- R[grepl("\\.(H[0-4]|nsnp|top)$", quantity_id)]
stopifnot(nrow(need) == 7 * length(ALLRUNS), length(ALLRUNS) == 2 * (2 * 2 + 2 * 3), !anyNA(need$value), all(need$value != ""), !any(need$value %in% c("NA", "NaN", "Inf")))
B <- rbindlist(branches)
notab <- function(dt) for (cc in names(dt)) if (any(grepl("[\t\n\r]", as.character(dt[[cc]])))) stop("tab/newline inside field of column ", cc)
notab(R); notab(B)
fwrite(R, file.path(OUT, "results.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(B, file.path(OUT, "branches.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(ALLRUNS), file.path(OUT, "coloc_runs_wide.tsv"), sep = "\t")
fwrite(rbindlist(lapply(JOINS, as.data.table)), file.path(OUT, "join_audit.tsv"), sep = "\t")
log_("wrote results.tsv (%d rows; %d coloc-quantity rows, none NA), branches.tsv (%d rows), coloc_runs_wide.tsv, join_audit.tsv, mutation_tests.tsv (%d/%d tripped)",
     nrow(R), nrow(need), nrow(B), sum(MUT$tripped), nrow(MUT))
log_("coloc step done in %.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs")))
