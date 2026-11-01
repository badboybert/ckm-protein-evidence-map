# b01_ld.R -- Workstream A (AMEND2 v3, sections 2, 3, 6): variant identity, panel availability, LD (r2), verdict.
# Builder implementation: R + PLINK 1.9. Run:  Rscript b01_ld.R        (about 2-3 min; reads inputs in place, writes only inside PK)
# No download, no new fine-mapping, no MR. Inputs are md5-pinned; the pre-declaration md5 is asserted first (STOP if it differs).
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
OUTD <- file.path(PK, "ld"); dir.create(OUTD, showWarnings = FALSE)
LOG  <- file.path(PK, "run.log")
logf <- function(...) { l <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(...)); cat(l, "\n", sep = ""); cat(l, "\n", sep = "", file = LOG, append = TRUE) }
md5 <- function(p) unname(tools::md5sum(p))
fmt <- function(x) if (length(x) == 0 || is.na(x)) "" else if (is.numeric(x)) sprintf("%.17g", x) else as.character(x)

# ---- 0. governing document -------------------------------------------------------------------------------------------
GOV <- file.path(ROOT, "V3/R2_rewrite/review_R2/P31R_R2_STAGE2_AMEND2_v3.md")
if (!identical(md5(GOV), "7daeec61dde5bd9a98956498c44a9f34")) stop("STOP: pre-declaration md5 differs: ", md5(GOV))
logf("pre-declaration md5 OK 7daeec61dde5bd9a98956498c44a9f34")

# ---- 1. inputs: md5 snapshot (start) ----------------------------------------------------------------------------------
L3 <- file.path(ROOT, "V3/V4_local/staging/l3-finemap/ARG1_6q23")
W4B <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_04b_arg1sets")
P503 <- file.path(ROOT, "paper 6/analysis/g1000_eur/g1000_eur")
P489 <- file.path(ROOT, "lhcmr_eas_refs/1000G_EUR_Phase3_plink/1000G.EUR.QC.6")
PLINK <- file.path(ROOT, "paper 6/analysis/tools/plink19/plink.exe")
INPUTS <- c(GOV = GOV, signals = file.path(L3, "conditional/signals.tsv"), harmonised = file.path(L3, "harmonised.tsv.gz"),
            prep_stamp = file.path(L3, "PREP_STAMP.json"), panel2_stamp = file.path(L3, "snpsets/PANEL2_STAMP.json"),
            p1np = file.path(L3, "snpsets/p1np.tsv"), p2np = file.path(L3, "snpsets/p2np.tsv"),
            ld503_bed = paste0(P503, ".bed"), ld503_bim = paste0(P503, ".bim"), ld503_fam = paste0(P503, ".fam"),
            ld489_bed = paste0(P489, ".bed"), ld489_bim = paste0(P489, ".bim"), ld489_fam = paste0(P489, ".fam"),
            plink_exe = PLINK, w4b_results = file.path(W4B, "results.tsv"), w4b_credsets = file.path(W4B, "credible_sets_by_method.tsv"),
            w4b_data_rds = file.path(W4B, "data/w4b_data.rds"),
            chen_xlsx = file.path(ROOT, "V3/R1_stage1/S1R1_R2_09_priorart_downloads_20261010/chen2023_41588_2022_1270_MOESM4_ESM.xlsx"),
            ukbppp_arg1 = file.path(ROOT, "weight-loss.omics/data/proteome/pqtl_regional/OID21445_chr6_cis.tsv"),
            clsa_arginine_window = file.path(ROOT, "V3/V4_local/raw/d1-bridge-eur/chen_windows/GCST90200372__chr6_130970832_132084332.tsv.gz"),
            finngen_window = file.path(W4B, "data/finngen_R10_T2D_ARG1window.tsv.gz"),
            diamante_window = file.path(W4B, "data/diamante2018_noUKB_EUR__chr6_131291089_132406428.tsv.gz"))
EXPECT <- c(signals = "b1de3285132999f1ae78cf119b8572e9", harmonised = "c394b6a1c0acbe32bfabc65d136ac974",
            ld503_bed = "543f92fba6ffe987dabf848b169c968d", ld503_bim = "48fbc06e2ec736fbcc58989cb1727322", ld503_fam = "6935b96c7e6f2a8f6122a049969bd0b8",
            ld489_bed = "4efd6cc853bc41b7daeb8a46c3d2fcd4", ld489_bim = "07a98ae97b76dc4df5672ca5ba2bcf19", ld489_fam = "fd9673c9469e18b200737b4e5f827c5b",
            plink_exe = "974d486f5c803047af9e05ba6fff9a2d", chen_xlsx = "9a245ed156dd900b4e53c9ff4c5a7586",
            ukbppp_arg1 = "a6d1f14e0f8f171911495efc54b5aebf", clsa_arginine_window = "dbd68a8cc11b8735353f24994408169b")
snap <- function(tag) {
  d <- data.frame(id = names(INPUTS), path = unname(INPUTS), stringsAsFactors = FALSE)
  d$bytes <- file.size(d$path); d$md5 <- vapply(d$path, md5, ""); d$expected_md5 <- unname(EXPECT[d$id]); d$expected_md5[is.na(d$expected_md5)] <- ""
  stopifnot(!anyNA(d$bytes), all(d$expected_md5 == "" | d$expected_md5 == d$md5))
  write.table(d, file.path(PK, sprintf("staging_snapshot_%s.tsv", tag)), sep = "\t", quote = FALSE, row.names = FALSE); d }
S0 <- snap("start"); logf("staging snapshot (start) written: %d inputs, all expected md5s match", nrow(S0))

# ---- 2. variant identity and orientation, from source ----------------------------------------------------------------
sg <- read.delim(INPUTS[["signals"]], stringsAsFactors = FALSE)
rsA <- sg$rsid[sg$tid == "P__ARG1" & sg$signal == 1]; rsB <- sg$rsid[sg$tid == "P__ARG1" & sg$signal == 2]
rsC <- sg$rsid[sg$tid == "M__GCST90200372" & sg$signal == 1]          # arginine signal 1 in our CLSA analysis
logf("signals.tsv: protein signal 1 (A) = %s ; protein signal 2 (B) = %s ; arginine signal 1 = %s", rsA, rsB, rsC)
logf("the brief's [UNVERIFIED] names (A rs2781668, B rs17788484) %s", if (identical(c(rsA, rsB), c("rs2781668", "rs17788484"))) "AGREE with signals.tsv" else "DISAGREE with signals.tsv: using signals.tsv")
Chen_instr <- "rs2246012"   # Chen 2023 ST17 row for arginine (read from the xlsx by a01_chen_extract.py; asserted below)
stopifnot(identical(rsC, Chen_instr))
VARS <- c(C = Chen_instr, A = rsA, B = rsB)
H <- read.delim(gzfile(INPUTS[["harmonised"]]), stringsAsFactors = FALSE); stopifnot(nrow(H) == 3067)
hv <- H[match(VARS, H$rsid_ref), c("rsid_ref", "snp", "chr", "pos37", "pos38", "a1", "a2", "palindromic", "af_ref")]
stopifnot(!anyNA(hv$rsid_ref), all(hv$palindromic == 0))
logf("harmonised.tsv.gz (c394b6a1...): %s", paste(sprintf("%s %s:%s/%s pos37=%d pos38=%d", hv$rsid_ref, hv$chr, hv$a1, hv$a2, hv$pos37, hv$pos38), collapse = " | "))
ce <- read.delim(file.path(PK, "chen_estimates.tsv"), header = FALSE, skip = 1, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
st17 <- ce[ce$V1 == "rs2246012" & ce$V2 == "ST17" & grepl("^instrument", ce$V4), ]
stopifnot(nrow(st17) == 1, st17$V6 == 131577068, st17$V7 == 131898208)    # pos38, pos37 from Chen ST17
stopifnot(hv$pos38[1] == st17$V6, hv$pos37[1] == st17$V7)
logf("Chen ST17: rs2246012 POS_hg38 %s, POS_hg37 %s = harmonised pos38/pos37 (no lift-over performed; no chain file used)", st17$V6, st17$V7)

# ---- 3. PLINK per panel -----------------------------------------------------------------------------------------------
plink <- function(args, tag) {
  out <- system2(PLINK, args, stdout = TRUE, stderr = TRUE); st <- attr(out, "status"); if (is.null(st)) st <- 0L
  logf("PLINK %s: exit %d ; args: %s", tag, st, paste(args, collapse = " "))
  if (st != 0) stop("PLINK failed: ", tag, "\n", paste(tail(out, 10), collapse = "\n")); invisible(out) }
logf("PLINK version: %s", system2(PLINK, "--version", stdout = TRUE)[1])
Q <- function(x) shQuote(x, type = "cmd")
snplist <- file.path(OUTD, "three_snps.txt"); writeLines(unname(VARS), snplist)
PANELS <- list(`503` = P503, `489` = P489)
fam <- lapply(PANELS, function(p) read.table(paste0(p, ".fam"), stringsAsFactors = FALSE))
iids489 <- fam[["489"]]$V2; stopifnot(all(iids489 %in% fam[["503"]]$V2), nrow(fam[["503"]]) == 503, length(iids489) == 489)
keep489 <- file.path(OUTD, "keep489.txt"); write.table(fam[["489"]][, 1:2], keep489, quote = FALSE, row.names = FALSE, col.names = FALSE)
RUNS <- list(`503` = list(p = P503, keep = NULL), `489` = list(p = P489, keep = NULL), `503in489ids` = list(p = P503, keep = keep489))  # third = diagnostic only
RES <- list(); GT <- list()
for (rn in names(RUNS)) {
  pre <- file.path(OUTD, paste0("p", rn, "_3snp")); keep <- if (is.null(RUNS[[rn]]$keep)) character() else c("--keep", Q(RUNS[[rn]]$keep))
  # extract the three variants (id-matched within the locus window so a duplicated rsID elsewhere cannot enter); allele order kept
  plink(c("--bfile", Q(RUNS[[rn]]$p), "--chr", "6", "--from-bp", "131800000", "--to-bp", "132000000", "--extract", Q(snplist),
          keep, "--keep-allele-order", "--make-bed", "--out", Q(pre)), paste0(rn, " extract"))
  bim <- read.table(paste0(pre, ".bim"), stringsAsFactors = FALSE, col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
  avail <- VARS %in% bim$id
  names(avail) <- names(VARS)
  logf("panel %s availability: %s", rn, paste(sprintf("%s=%s", VARS, avail), collapse = " "))
  if (!all(avail)) { RES[[rn]] <- list(avail = avail); next }
  stopifnot(nrow(bim) == 3, !anyDuplicated(bim$id))
  # PLINK frequencies (allele counts), missingness, and unphased r2 with all pairs retained
  plink(c("--bfile", Q(pre), "--keep-allele-order", "--freq", "counts", "--out", Q(paste0(pre, "_freq"))), paste0(rn, " freq"))
  plink(c("--bfile", Q(pre), "--keep-allele-order", "--missing", "--out", Q(paste0(pre, "_miss"))), paste0(rn, " missing"))
  plink(c("--bfile", Q(pre), "--keep-allele-order", "--r2", "--ld-window-r2", "0", "--ld-window", "99999", "--ld-window-kb", "1000", "--out", Q(paste0(pre, "_ld"))), paste0(rn, " r2"))
  plink(c("--bfile", Q(pre), "--keep-allele-order", "--recode", "A", "--out", Q(pre)), paste0(rn, " recode A"))
  fq <- read.table(paste0(pre, "_freq.frq.counts"), header = TRUE, stringsAsFactors = FALSE)
  ms <- read.table(paste0(pre, "_miss.lmiss"), header = TRUE, stringsAsFactors = FALSE)
  ld <- read.table(paste0(pre, "_ld.ld"), header = TRUE, stringsAsFactors = FALSE)
  stopifnot(nrow(ld) == 3)       # all three pairs retained, r2 < 0.20 included (window-r2 0)
  rw <- read.table(paste0(pre, ".raw"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  G <- as.matrix(rw[, 7:9]); colnames(G) <- sub("_.*$", "", colnames(G)); G <- G[, VARS, drop = FALSE]
  a1 <- sub("^.*_", "", colnames(rw)[7:9]); names(a1) <- sub("_.*$", "", colnames(rw)[7:9])
  stopifnot(nrow(G) == (if (rn == "489" || rn == "503in489ids") 489 else 503))
  RES[[rn]] <- list(avail = avail, bim = bim, fq = fq, ms = ms, ld = ld, G = G, a1 = a1)
}

# ---- 4. R recomputation (squared Pearson correlation of 0/1/2 counts; pairwise complete; no mean imputation) ------------
pair <- function(G, x, y) { ok <- complete.cases(G[, x], G[, y]); n <- sum(ok)
  v <- c(var(G[ok, x]), var(G[ok, y])); r <- if (n > 2 && all(v > 0)) cor(G[ok, x], G[ok, y]) else NA_real_
  list(n = n, r = r, r2 = r^2, mono = any(v == 0)) }
PAIRS <- list(c("C", "A"), c("C", "B"), c("A", "B"))
rows <- list(); add <- function(id, value, panel = "", n = NA, notes = "") rows[[length(rows) + 1]] <<- data.frame(id = id, value = fmt(value), panel = panel, n = fmt(n), notes = notes, stringsAsFactors = FALSE)
CAT <- list()
add("input.preregistration.md5", "7daeec61dde5bd9a98956498c44a9f34", "", NA, "AMEND2 v3, asserted before any computation")
add("variant.A.rsid", rsA, "", NA, "signals.tsv P__ARG1 signal 1 (z_marginal withheld_pending_terms; pos37 131897278)")
add("variant.B.rsid", rsB, "", NA, "signals.tsv P__ARG1 signal 2 (pos37 131894358)")
add("variant.C.rsid", rsC, "", NA, "signals.tsv M__GCST90200372 (arginine, CLSA) signal 1 = Chen ST17 instrument for arginine")
for (k in seq_len(3)) { add(sprintf("variant.%s.pos37_pos38_a1_a2", names(VARS)[k]), sprintf("%d_%d_%s_%s", hv$pos37[k], hv$pos38[k], hv$a1[k], hv$a2[k]), "", NA, "harmonised.tsv.gz (c394b6a1); a1/a2 = reference-panel allele order; non-palindromic") }
for (rn in names(RES)) {
  R <- RES[[rn]]
  add(sprintf("avail.all3.%s", rn), all(R$avail), rn, NA, paste(sprintf("%s=%s", VARS, R$avail), collapse = "; "))
  if (!all(R$avail)) next
  n_samp <- nrow(R$G)
  add(sprintf("panel.n_samples.%s", rn), n_samp, rn, n_samp, "individuals in the extracted panel (all founders: PLINK --r2 uses founders; fam parent ids all 0)")
  for (k in seq_len(3)) {
    id <- VARS[k]; fr <- R$fq[R$fq$SNP == id, ]; ms <- R$ms[R$ms$SNP == id, ]; b <- R$bim[R$bim$id == id, ]
    ncall <- n_samp - ms$N_MISS; nal <- fr$C1 + fr$C2
    add(sprintf("bim.%s.%s", names(VARS)[k], rn), sprintf("%d:%s/%s", b$pos, b$a1, b$a2), rn, NA, "panel .bim chr6 position (GRCh37), A1/A2")
    add(sprintf("count.A1.%s.%s", names(VARS)[k], rn), fr$C1, rn, nal, sprintf("PLINK --freq counts: A1=%s count %d, A2 count %d, called alleles %d (=2 x %d called individuals), missing individuals %d", fr$A1, fr$C1, fr$C2, nal, ncall, ms$N_MISS))
    add(sprintf("freq.A1.%s.%s", names(VARS)[k], rn), fr$C1 / nal, rn, ncall, sprintf("A1=%s (A1 frequency in this panel)", fr$A1))
    add(sprintf("monomorphic.%s.%s", names(VARS)[k], rn), fr$C1 == 0 || fr$C2 == 0, rn, ncall, "")
  }
  for (pr in PAIRS) {
    a <- pair(R$G, VARS[[pr[1]]], VARS[[pr[2]]]); lk <- R$ld[(R$ld$SNP_A == VARS[[pr[1]]] & R$ld$SNP_B == VARS[[pr[2]]]) | (R$ld$SNP_B == VARS[[pr[1]]] & R$ld$SNP_A == VARS[[pr[2]]]), ]
    stopifnot(nrow(lk) == 1)
    d <- abs(lk$R2 - a$r2); stopifnot(d < 5e-6)   # PLINK prints 6 significant digits
    add(sprintf("r2.%s_%s.%s", pr[1], pr[2], rn), a$r2, rn, a$n, sprintf("%s vs %s: squared Pearson r of 0/1/2 A1-allele counts, pairwise complete (R, unrounded); PLINK 1.9 --r2 %s; |R - PLINK| %.2g; signed r %.6f", VARS[[pr[1]]], VARS[[pr[2]]], fmt(lk$R2), d, a$r))
    add(sprintf("r2.plink.%s_%s.%s", pr[1], pr[2], rn), lk$R2, rn, a$n, "PLINK 1.9 --r2 --ld-window-r2 0 output (6 significant digits as printed)")
    RES[[rn]]$r2[[paste(pr, collapse = "_")]] <- a
  }
}
# ---- 5. verdict (AMEND2 v3 section 3.4), unrounded values -------------------------------------------------------------
categorize <- function(rA, rB, evaluable) {
  if (!evaluable || is.na(rA) || is.na(rB)) return(c("UNRESOLVED (a required comparison is not evaluable)", "not evaluable"))
  if (rA >= 0.80 && rB < 0.20) return(c("STRONGLY TAGS A", sprintf("r2 with A %.17g >= 0.80 and with B %.17g < 0.20", rA, rB)))
  if (rB >= 0.80 && rA < 0.20) return(c("STRONGLY TAGS B", sprintf("r2 with B %.17g >= 0.80 and with A %.17g < 0.20", rB, rA)))
  if (rA < 0.20 && rB < 0.20) return(c("NOT STRONGLY TAGGED", sprintf("r2 below 0.20 with both leads (A %.17g, B %.17g): no strong pairwise LD with either lead", rA, rB)))
  c("PARTIAL / UNRESOLVED", sprintf("intermediate or high LD with both leads (A %.17g, B %.17g)", rA, rB)) }
for (rn in c("503", "489", "503in489ids")) {
  R <- RES[[rn]]; ev <- all(R$avail)
  if (ev) { ra <- R$r2[["C_A"]]; rb <- R$r2[["C_B"]]; ev <- !(ra$mono || rb$mono) && !is.na(ra$r2) && !is.na(rb$r2); ct <- categorize(ra$r2, rb$r2, ev) }
  else ct <- categorize(NA, NA, FALSE)
  CAT[[rn]] <- ct
  add(sprintf("category.%s", rn), ct[1], rn, NA, paste0(ct[2], if (rn == "503in489ids") " [DIAGNOSTIC ONLY: 503 panel restricted to the 489 individuals; not part of the verdict]" else ""))
}
if (CAT[["503"]][1] == CAT[["489"]][1] && !grepl("^UNRESOLVED", CAT[["503"]][1])) { fin <- CAT[["503"]][1]; why <- "503 (primary) and 489 (nested subset) give the same category" } else if (grepl("^UNRESOLVED", CAT[["503"]][1]) || grepl("^UNRESOLVED", CAT[["489"]][1])) { fin <- "UNRESOLVED (a required comparison is not evaluable)"; why <- "503: %s; 489: %s" } else { fin <- "PARTIAL / UNRESOLVED"; why <- "different categories in the two panels (503: %s; 489: %s)" }
why <- sprintf(why, CAT[["503"]][1], CAT[["489"]][1])
add("category.final", fin, "503 primary; 489 must agree", NA, why)
write.table(do.call(rbind, rows), file.path(PK, "ld_results_part.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(RES = RES, VARS = VARS, hv = hv, CAT = CAT, fin = fin), file.path(OUTD, "ld_objects.rds"))
logf("VERDICT 503: %s | 489: %s | final: %s", CAT[["503"]][1], CAT[["489"]][1], fin)
S1 <- snap("end"); stopifnot(identical(S0[, c("id", "md5")], S1[, c("id", "md5")])); logf("staging snapshot (end) equals start")
logf("done; versions: %s ; PLINK md5 %s", R.version.string, md5(PLINK))
