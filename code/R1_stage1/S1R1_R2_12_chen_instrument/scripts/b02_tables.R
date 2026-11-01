# b02_tables.R -- Workstream A: credible-set annotations (existing outputs only), existing estimates at the three variants,
# the locus comparison table, and the assembled results.tsv.   Run after a01_chen_extract.py and b01_ld.R.
# Reads existing files only. No new fine-mapping, no new association estimates, no MR, no download.
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
options(stringsAsFactors = FALSE, digits = 15, scipen = 0, warn = 1)
ROOT <- .P31R_CKM
PK   <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_12_chen_instrument")
W4B  <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_04b_arg1sets")
L3   <- file.path(ROOT, "V3/V4_local/staging/l3-finemap/ARG1_6q23")
logf <- function(...) { l <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(...)); cat(l, "\n", sep = ""); cat(l, "\n", sep = "", file = file.path(PK, "run.log"), append = TRUE) }
fmt <- function(x) vapply(x, function(z) if (is.na(z)) "" else if (is.numeric(z)) sprintf("%.17g", z) else as.character(z), "", USE.NAMES = FALSE)
O <- readRDS(file.path(PK, "ld/ld_objects.rds")); VARS <- O$VARS; hv <- O$hv      # names C, A, B
stopifnot(identical(names(VARS), c("C", "A", "B")))
rd <- function(p, ...) { x <- read.delim(p, stringsAsFactors = FALSE, ...); stopifnot(nrow(x) > 0); x }   # (not fread: spaces in path)
md5 <- function(p) unname(tools::md5sum(p))
`%+%` <- function(a, b) paste0(a, b)

# ================================================================== 1. credible-set annotations (existing outputs) =========
CS <- rd(file.path(W4B, "credible_sets_by_method.tsv"), quote = "")
stopifnot(nrow(CS) == 12, md5(file.path(W4B, "credible_sets_by_method.tsv")) == md5(file.path(W4B, "credible_sets_by_method.tsv")))
D  <- readRDS(file.path(W4B, "data/w4b_data.rds")); Hh <- D$H
inset <- function(setname, rs) { snp <- Hh$snp[match(rs, Hh$rsid_ref)]; snp %in% D$SETS[[setname]] }
hasz  <- function(col, rs) !is.na(Hh[[col]][match(rs, Hh$rsid_ref)])
METH <- c(`1_T2D` = "(1) SuSiE-RSS, 503 reference (p1np), FinnGen R10 T2D", `1_protein` = "(1) SuSiE-RSS, 503 reference (p1np), UKB-PPP ARG1 protein",
          `2` = "(2) UKB-PPP in-sample SuSiE (Sun 2023 ST16), ARG1 protein", `3` = "(3) shared ARG1-T2D(FinnGen) 95% ABF set, conditioned on the other signal, z-space, 503 LD",
          `4` = "(4) T2D(FinnGen) 95% ABF set, conditioned on the other signal, z-space, 503 LD", `5` = "(5) shared ARG1-T2D(DIAMANTE no-UKB) 95% set, conditioned on the other signal, z-space, 503 LD")
OUTCOME <- c(`3` = "z__FG", `4` = "z__FG", `5` = "z__DM")
adm_of <- function(m, adm_col) {
  if (m %in% c("1_T2D", "1_protein")) return(sub("^(admissible|inadmissible).*", "\\1", adm_col))
  if (m == "2") return("admissible")
  "unresolved (not a SuSiE fit, so the Rule 3(i) wording does not cover it; Rule 3(ii) qualifier shown in the robustness column)" }
rows <- list()
for (i in seq_len(nrow(CS))) {
  m <- CS$method[i]; X <- CS$signal[i]; mem503 <- strsplit(CS$members[i], ";")[[1]]
  variants <- list(`503` = mem503)
  if (m %in% c("3", "4", "5")) { e <- sub("^.*489 set: ", "", CS$extra[i]); e <- sub(" \\(.*\\)\\s*$", "", e); variants[["489"]] <- strsplit(e, ";")[[1]] }
  for (pn in names(variants)) {
    mem <- variants[[pn]]
    nocs <- CS$size[i] == 0 || grepl("NO CREDIBLE SET", CS$members[i])
    for (k in names(VARS)) {
      rs <- VARS[[k]]
      st <- if (nocs) "no available credible set/output (no set returned at this signal)" else if (rs %in% mem) "included" else if (m == "2") "not in the listed set; whether the variant was in the fit's input is not recorded (external ST16 list)" else {
        sn <- if (pn == "503") "p1np" else "p2np"
        ok <- inset(sn, rs) && (if (m %in% names(OUTCOME)) hasz(OUTCOME[[m]], rs) else TRUE)
        if (ok) "tested but not in the set" else "not in the fit's input variants" }
      rows[[length(rows) + 1]] <- data.frame(set_signal = X, method = m, method_label = METH[[m]], panel = if (m == "2") "UKB-PPP in-sample LD" else pn, variant = rs, variant_role = k,
        status = st, set_size = if (nocs) 0L else length(mem), admissibility_rule3i = adm_of(m, CS$admissibility[i]), admissibility_text_in_source = CS$admissibility[i], robustness_rule3ii = CS$robustness[i], source_extra = CS$extra[i])
    }
  }
}
CSS <- do.call(rbind, rows); write.table(CSS, file.path(PK, "credible_set_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
logf("credible-set status rows: %d (12 stored sets; methods 3-5 also have the 489 version)", nrow(CSS))
# cross-check against the stored results.tsv membership rows (independent read of a second W4b file)
RS <- rd(file.path(W4B, "results.tsv"), quote = "\"")
chk <- 0L
for (i in seq_len(nrow(CS))) { id <- sprintf("S4.cs.%s.%s.members", CS$method[i], CS$signal[i]); r <- RS[RS$quantity_id == id, ]; stopifnot(nrow(r) == 1)
  if (!grepl("NO CREDIBLE", r$value)) stopifnot(identical(r$value, CS$members[i])); chk <- chk + 1L }
logf("cross-check: members string identical in credible_sets_by_method.tsv and results.tsv for %d sets", chk)
# set-level S1R1_R2_04b brief expectation (the 4 sets named in AMEND2 section 4 contain rs2246012): confirm from the data
named <- c("S4.cs.2.A.members", "S4.cs.1_protein.A.members", "S4.cs.1_T2D.A.members", "S4.cs.4.A.members")
for (id in named) { r <- RS[RS$quantity_id == id, ]; logf("AMEND2 s4 check %s contains rs2246012: %s", id, grepl("rs2246012", r$value)) }

# ================================================================== 2. existing estimates at the three variants =============
pos38 <- setNames(hv$pos38, hv$rsid_ref); a1 <- setNames(hv$a1, hv$rsid_ref); a2 <- setNames(hv$a2, hv$rsid_ref)
est <- list()
addE <- function(rs, source, trait, ea, nea, eaf, beta, se, p, scale, prov, extra = "") {
  flip <- if (is.na(ea)) NA else if (ea == a1[[rs]]) FALSE else if (ea == a2[[rs]]) TRUE else NA
  est[[length(est) + 1]] <<- data.frame(rsid = rs, source = source, trait = trait, effect_allele = ea, other_allele = nea, eaf_source = eaf, beta_or_logOR_source = beta, se_source = se, p_source = p,
    scale_as_reported = scale, panel_a1 = a1[[rs]], beta_aligned_to_panel_a1 = if (isTRUE(flip)) -as.numeric(beta) else if (isFALSE(flip)) as.numeric(beta) else NA,
    eaf_aligned_to_panel_a1 = if (isTRUE(flip)) 1 - as.numeric(eaf) else if (isFALSE(flip)) as.numeric(eaf) else NA,
    provenance = prov, note = extra, stringsAsFactors = FALSE) }
# UKB-PPP OID21445 (raw regional file, GRCh38 position, md5 a6d1f14e)
u <- rd(file.path(ROOT, "weight-loss.omics/data/proteome/pqtl_regional/OID21445_chr6_cis.tsv"))
for (rs in VARS) { r <- u[u$pos == pos38[[rs]], ]; stopifnot(nrow(r) == 1)
  addE(rs, "UKB-PPP OID21445 (ARG1 protein), regional file", "ARG1 plasma protein", r$ALLELE1, r$ALLELE0, r$A1FREQ, r$BETA, r$SE, 10^(-r$LOG10P), "per-allele beta on the file's protein scale (as reported; not rescaled); p from LOG10P",
       "existing file, md5 a6d1f14e0f8f171911495efc54b5aebf; row matched on GRCh38 position; ID " %+% r$ID) }
# CLSA arginine GCST90200372 window (Chen metabolite GWAS release; md5 dbd68a8c)
cl <- read.delim(gzfile(file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90200372__chr6_130970832_132084332.tsv.gz")), stringsAsFactors = FALSE); stopifnot(nrow(cl) > 0)
for (rs in VARS) { r <- cl[cl$rsid == rs, ]; stopifnot(nrow(r) == 1, r$base_pair_location == pos38[[rs]])
  addE(rs, "CLSA arginine GCST90200372 (GWAS Catalog window)", "arginine", r$effect_allele, r$other_allele, r$effect_allele_frequency, r$beta, r$standard_error, r$p_value, "per-allele beta, metabolite scale as reported by the source (not rescaled)",
       "existing window file, md5 dbd68a8cc11b8735353f24994408169b; matched on rsid and GRCh38 position") }
# FinnGen R10 T2D window and DIAMANTE no-UKB window (W4b extracts)
fg <- read.delim(gzfile(file.path(W4B, "data/finngen_R10_T2D_ARG1window.tsv.gz")), stringsAsFactors = FALSE); stopifnot(nrow(fg) > 0)
for (rs in VARS) { r <- fg[fg$pos == pos38[[rs]], ]; stopifnot(nrow(r) == 1, grepl(rs, r$rsids))
  addE(rs, "FinnGen R10 T2D (W4b window extract)", "type 2 diabetes", r$alt, r$ref, r$af_alt, r$beta, r$sebeta, r$pval, "per-allele log-odds (as reported)",
       "existing W4b extract data/finngen_R10_T2D_ARG1window.tsv.gz; matched on GRCh38 position; alt allele = effect allele") }
dm <- read.delim(gzfile(file.path(W4B, "data/diamante2018_noUKB_EUR__chr6_131291089_132406428.tsv.gz")), stringsAsFactors = FALSE); stopifnot(nrow(dm) > 0)
for (rs in VARS) { r <- dm[dm$Pos == hv$pos37[hv$rsid_ref == rs], ]; stopifnot(nrow(r) == 1)
  addE(rs, "DIAMANTE 2018 no-UKB EUR (W4b window extract)", "type 2 diabetes", r$EA, r$NEA, r$EAF, r$Beta, r$SE, r$Pvalue, "per-allele log-odds (as reported; values rounded to 2-3 sig. digits by the source file)",
       "existing W4b extract; matched on GRCh37 position (the file is GRCh37); 55,005 cases / 400,308 controls") }
# Chen 2023 (publisher supplement; rows read by a01_chen_extract.py)
ce <- read.delim(file.path(PK, "chen_estimates.tsv"), header = FALSE, skip = 1, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
g <- function(v, sh, pat) { r <- ce[ce$V1 == v & ce$V2 == sh & grepl(pat, ce$V4), ]; stopifnot(nrow(r) == 1); r }
for (v in c("rs2246012", "rs17788484")) { r <- g(v, "ST5", "arginine GWAS")
  addE(v, "Chen 2023 ST5 (CLSA arginine GWAS row)", "arginine", r$V7, r$V8, as.numeric(r$V9), as.numeric(r$V10), as.numeric(r$V11), as.numeric(r$V12), "per-allele beta, metabolite scale as reported in ST5 (not rescaled)",
       sprintf("publisher xlsx md5 9a245ed1...; ST5 row %s; N %s", r$V3, r$V13), "Chen's GWAS-catalogue-release row; the same CLSA association as the GCST90200372 file above") }
r <- g("rs2246012", "ST17", "T2D block")
addE("rs2246012", "Chen 2023 ST17 T2D block (instrument-outcome row)", "type 2 diabetes (Chen outcome; Mahajan 2018 total per AMEND2 s4)", r$V5, r$V6, as.numeric(r$V7), as.numeric(r$V8), as.numeric(r$V9), as.numeric(r$V10), "per-allele log-odds as reported (ST17 columns 68-73)",
     sprintf("publisher xlsx; ST17 row %s; block located by position, validated: exp(0.037/0.173255)=1.2381 equals ST19 OR", r$V3), "outcome cohort is NOT our DIAMANTE file (74,124/824,006 per AMEND2 s4; provenance equivalence not claimed)")
r <- g("rs2246012", "ST19", "Wald ratio")
addE("rs2246012", "Chen 2023 ST19 Wald ratio arginine->T2D", "arginine -> T2D (Chen MR)", NA, NA, NA, NA, NA, as.numeric(r$V11), "odds ratio per Chen arginine unit",
     sprintf("publisher xlsx; ST19 row %s; OR %s (95%% CI %s to %s); P.adj %s; nSNP %s. Not re-derived or rescaled here", r$V3, r$V8, r$V9, r$V10, r$V12, r$V7), "Chen's Wald ratio, one-SNP; reported as published")
r <- g("rs2246012", "ST24", "coloc")
addE("rs2246012", "Chen 2023 ST24 colocalization arginine/T2D", "arginine / T2D (Chen coloc.abf)", NA, NA, NA, NA, NA, NA, "posterior probabilities",
     sprintf("publisher xlsx; ST24 row %s; nSNP %s; PP.H0 %s PP.H1 %s PP.H2 %s PP.H3 %s PP.H4 %s", r$V3, r$V7, r$V8, r$V9, r$V10, r$V11, r$V12), "reported as published; no re-computation")
E <- do.call(rbind, est); write.table(E, file.path(PK, "estimates_at_variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
crl <- readLines(file.path(PK, "chen_rows.tsv"), encoding = "UTF-8"); stopifnot(length(crl) > 5, !any(grepl("rs2781668", crl)))
logf("Chen supplement: rs2781668 appears in none of ST5/ST14/ST17/ST19/ST24 (searched by rsID; see chen_rows.tsv)")
# W4b-vs-raw consistency (aligned betas): harmonised b__P__ARG1 / arginine vs the raw rows above
for (rs in VARS) { h <- Hh[match(rs, Hh$rsid_ref), ]
  e1 <- E[E$rsid == rs & grepl("^UKB-PPP", E$source), ]; e2 <- E[E$rsid == rs & grepl("^CLSA", E$source), ]
  stopifnot(abs(h$b__P__ARG1 - e1$beta_aligned_to_panel_a1) < 1e-9, abs(h$b__M__GCST90200372 - e2$beta_aligned_to_panel_a1) < 1e-9) }
logf("harmonised aligned betas (protein, arginine) equal the raw-file betas aligned to panel A1 at all three variants")

# ================================================================== 3. locus comparison table ===============================
R <- O$RES; gv <- function(rn, key) { x <- read.delim(file.path(PK, "ld_results_part.tsv"), stringsAsFactors = FALSE); x[x$id == key & x$panel == rn, ] }
LDP <- read.delim(file.path(PK, "ld_results_part.tsv"), stringsAsFactors = FALSE)
val <- function(id, pn) { r <- LDP[LDP$id == id & LDP$panel == pn, ]; stopifnot(nrow(r) == 1); r }
tab <- list()
for (k in names(VARS)) {
  rs <- VARS[[k]]; ee <- E[E$rsid == rs, ]
  cs_sum <- function(sig) { x <- CSS[CSS$variant == rs & CSS$set_signal == sig, ]; paste(sprintf("%s[%s;%s;%s]", x$method, x$panel, sub(" \\(.*$", "", x$status), sub(" \\(.*$", "", x$admissibility_rule3i)), collapse = " | ") }
  rA <- if (k == "C") c(val("r2.C_A.503", "503"), val("r2.C_A.489", "489")) else NULL
  tab[[k]] <- data.frame(role = c(C = "Chen 2023 instrument (arginine; effector gene ARG1)", A = "ARG1 protein signal A lead", B = "ARG1 protein signal B lead")[[k]], rsid = rs,
    panel_A1_A2_GRCh37 = sprintf("%s/%s @ %d", hv$a1[hv$rsid_ref == rs], hv$a2[hv$rsid_ref == rs], hv$pos37[hv$rsid_ref == rs]), GRCh38_position = hv$pos38[hv$rsid_ref == rs],
    freq_panel503_A1 = fmt(as.numeric(val(sprintf("freq.A1.%s.503", k), "503")$value)), freq_panel489_A1 = fmt(as.numeric(val(sprintf("freq.A1.%s.489", k), "489")$value)),
    freq_by_source = { fe <- ee[!is.na(ee$eaf_source), ]; paste(sprintf("%s: %s=%s", sub(" \\(.*$", "", fe$source), fe$effect_allele, signif(as.numeric(fe$eaf_source), 4)), collapse = "; ") },
    r2_with_rs2246012_503 = if (k == "C") "self" else fmt(as.numeric(val(sprintf("r2.C_%s.503", k), "503")$value)), n_503 = if (k == "C") "" else val(sprintf("r2.C_%s.503", k), "503")$n,
    r2_with_rs2246012_489 = if (k == "C") "self" else fmt(as.numeric(val(sprintf("r2.C_%s.489", k), "489")$value)), n_489 = if (k == "C") "" else val(sprintf("r2.C_%s.489", k), "489")$n,
    credible_sets_signal_A = cs_sum("A"), credible_sets_signal_B = cs_sum("B"),
    estimates = paste(sprintf("%s [%s; EA %s; beta %s; SE %s; scale: %s]", sub(" \\(.*$", "", ee$source), ee$trait, ee$effect_allele, fmt(as.numeric(ee$beta_or_logOR_source)), fmt(as.numeric(ee$se_source)), ee$scale_as_reported), collapse = " || "), stringsAsFactors = FALSE) }
T <- do.call(rbind, tab); write.table(T, file.path(PK, "locus_comparison_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pairs <- data.frame(pair = c("rs2246012 vs rs2781668 (A)", "rs2246012 vs rs17788484 (B)", "rs2781668 (A) vs rs17788484 (B)"), key = c("C_A", "C_B", "A_B"), stringsAsFactors = FALSE)
P <- do.call(rbind, lapply(seq_len(3), function(i) data.frame(pair = pairs$pair[i], r2_503 = val(paste0("r2.", pairs$key[i], ".503"), "503")$value, n_503 = val(paste0("r2.", pairs$key[i], ".503"), "503")$n,
   r2_489 = val(paste0("r2.", pairs$key[i], ".489"), "489")$value, n_489 = val(paste0("r2.", pairs$key[i], ".489"), "489")$n, plink_503 = val(paste0("r2.plink.", pairs$key[i], ".503"), "503")$value, plink_489 = val(paste0("r2.plink.", pairs$key[i], ".489"), "489")$value)))
write.table(P, file.path(PK, "ld_pairs_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ================================================================== 4. assemble results.tsv =================================
res <- LDP
cs_rows <- data.frame(id = sprintf("cs.%s.%s.%s.%s", CSS$method, CSS$set_signal, CSS$panel, CSS$variant_role), value = CSS$status, panel = CSS$panel, n = CSS$set_size,
   notes = sprintf("%s; signal %s; variant %s; Rule 3(i) admissibility: %s; robustness (Rule 3(ii)): %s; source: credible_sets_by_method.tsv (existing W4b output, no new fit)", CSS$method_label, CSS$set_signal, CSS$variant, CSS$admissibility_rule3i, CSS$robustness_rule3ii), stringsAsFactors = FALSE)
est_rows <- do.call(rbind, lapply(seq_len(nrow(E)), function(i) { e <- E[i, ]; tag <- sprintf("est.%s.%s", gsub("[^A-Za-z0-9]+", "_", sub(" \\(.*$", "", e$source)), e$rsid)
  data.frame(id = paste0(tag, c(".beta_source", ".se_source", ".eaf_source", ".p_source")), value = c(fmt(as.numeric(e$beta_or_logOR_source)), fmt(as.numeric(e$se_source)), fmt(as.numeric(e$eaf_source)), fmt(as.numeric(e$p_source))),
   panel = "", n = "", notes = sprintf("%s; %s; effect allele %s (other %s); scale: %s; %s", e$source, e$trait, e$effect_allele, e$other_allele, e$scale_as_reported, e$provenance), stringsAsFactors = FALSE) }))
est_rows <- est_rows[est_rows$value != "", ]
chen_text <- E[grepl("ST19|ST24", E$source), ]; est_rows <- rbind(est_rows, data.frame(id = c("est.Chen2023_ST19.rs2246012.text", "est.Chen2023_ST24.rs2246012.text"), value = chen_text$provenance, panel = "", n = "", notes = chen_text$source, stringsAsFactors = FALSE))
all <- rbind(res, cs_rows, est_rows); write.table(all, file.path(PK, "results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
logf("results.tsv: %d rows (LD/verdict %d, credible-set %d, estimates %d)", nrow(all), nrow(res), nrow(cs_rows), nrow(est_rows))
