# 04_harmonise.R -- DEVIATION: engine path-locked. Re-implements in R the l3fm__prep.py QC and the l3fm__panel2.py /
# l3fm__robustness.py SNP sets (read, not run) for the 16 candidate windows:
#   join   protein x outcome on position (GRCh37 from the UKB-PPP ID for CKDGen; GRCh38 for FinnGen) + compatible alleles
#          (unique match; complement only when not palindromic); reference = 503 panel bim at GRCh37 position + compatible alleles.
#   QC     ref MAF >= 0.01; monomorphic dropped; palindromic kept only at ref MAF <= 0.40; |aligned GWAS AF - ref AF| <= 0.20
#          for every trait reporting a frequency.
#   sets   susie = all QC'd SNPs (503); p1np = non-palindromic (503); p1fr = p1np + palindromes whose every reporting trait has
#          GWAS MAF <= 0.40 and |GWAS freq of ref A1 - ref A1 freq| <= 0.10 (503); p2np = p1np SNPs in the 489 panel at the same
#          GRCh37 position and allele pair (unique), no missing genotype, panel-2 MAF >= 0.01, |AF2 - AF1| <= 0.05 (489).
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
source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "00_common.R"))
say("==== W6b 04_harmonise ====")
cand <- fread(file.path(WORK, "candidates_windows.tsv"), colClasses = list(character = "chr"))
nout <- fread(file.path(WORK, "outcome_n.tsv"))
EXD <- file.path(WORK, "extracts"); REFD <- file.path(WORK, "ref"); HD <- file.path(WORK, "harm")
dir.create(REFD, showWarnings = FALSE); dir.create(HD, showWarnings = FALSE)
man <- fread(file.path(EXD, "EXTRACT_MANIFEST.tsv"), colClasses = "character")
chk_extract <- function(f) { r <- man[file == basename(f)]; if (nrow(r) != 1 || md5(f) != r$md5) stop("extract md5 mismatch ", f); f }

# ---- reference panel files (registered: md5 asserted) ----
ref503 <- sub("\\.bed$", "", res("LD503_bed")); res("LD503_bim"); res("LD503_fam")
p2pre <- function(ch) {
  id_bed <- sprintf("LD489_chr%s_bed", ch); id_bim <- sprintf("LD489_chr%s_bim", ch)
  pre <- file.path(ROOT, sprintf("lhcmr_eas_refs/1000G_EUR_Phase3_plink/1000G.EUR.QC.%s", ch))
  if (id_bed %in% RES_IDS) { res(id_bed); res(id_bim) } else {
    reg_new(paste0(pre, ".bed"), sprintf("489 panel chr%s .bed (absent from RESOURCES; brief rule 3)", ch))
    reg_new(paste0(pre, ".bim"), sprintf("489 panel chr%s .bim (absent from RESOURCES; brief rule 3)", ch))
  }
  fam <- reg_new(paste0(pre, ".fam"), sprintf("489 panel chr%s .fam", ch))
  stopifnot(md5c(fam) == res_entry("LD489_fam")$md5)       # the 22 per-chromosome .fam files share one md5 (RESOURCES note)
  pre
}

# ---- one PLINK pass over the 503 panel for all windows: bim rows + A1 frequencies ----
rng <- unique(cand[, .(chr, lo37, hi37)])
rng[, label := paste0("w", seq_len(.N))]
rf <- file.path(REFD, "ranges37.txt"); fwrite(rng[, .(chr, lo37, hi37, label)], rf, sep = "\t", col.names = FALSE)
if (!file.exists(file.path(REFD, "ref503.frq.counts"))) {
  say("PLINK 503: --extract range (%d windows) --freq --make-just-bim", nrow(rng))
  run_plink(c("--bfile", q(ref503), "--extract", "range", q(rf), "--keep-allele-order", "--freq", "counts", "--make-just-bim"), file.path(REFD, "ref503"))
}
read_ref <- function(pre) {
  bim <- fread(paste0(pre, ".bim"), header = FALSE, col.names = c("chr", "rsid", "cm", "pos37", "a1", "a2"), colClasses = list(character = c(1, 2, 5, 6)))
  frq <- fread(paste0(pre, ".frq.counts"), colClasses = list(character = c("CHR", "SNP", "A1", "A2")))
  stopifnot(nrow(bim) == nrow(frq), all(bim$rsid == frq$SNP), all(bim$a1 == frq$A1), all(bim$a2 == frq$A2))
  # --keep-allele-order --freq counts: C1 = A1 allele count, C2 = A2 count, G0 = missing genotypes -> exact A1 frequency
  bim[, `:=`(af = frq$C1 / (frq$C1 + frq$C2), nchrobs = frq$C1 + frq$C2, a1 = toupper(a1), a2 = toupper(a2))]
  bim[]
}
R503 <- read_ref(file.path(REFD, "ref503"))
stopifnot(!anyDuplicated(R503$rsid))
say("503 panel rows in windows: %d (IDs unique)", nrow(R503))

# ---- 489 panel per chromosome ----
R489 <- list()
for (ch in unique(cand$chr)) {
  pre <- p2pre(ch)
  out <- file.path(REFD, paste0("ref489_chr", ch))
  r <- rng[chr == ch]
  rf2 <- file.path(REFD, paste0("ranges37_chr", ch, ".txt")); fwrite(r[, .(chr, lo37, hi37, label)], rf2, sep = "\t", col.names = FALSE)
  if (!file.exists(paste0(out, ".frq.counts"))) run_plink(c("--bfile", q(pre), "--extract", "range", q(rf2), "--keep-allele-order", "--freq", "counts", "--make-just-bim"), out)
  R489[[ch]] <- read_ref(out)
  stopifnot(!anyDuplicated(R489[[ch]]$rsid))
}

# ---- per candidate ----
harmonise <- function(cc) {
  O <- fread(chk_extract(file.path(EXD, paste0(cc$cid, "__outcome.tsv.gz"))), colClasses = list(character = c("chr", "rsid", "ea", "oa")))
  bad <- !is.finite(O$b) | !is.finite(O$se) | O$se <= 0
  drops <- list(outcome_nonfinite = sum(bad)); O <- O[!bad]
  O[, akey := paste(pos, pmin(ea, oa), pmax(ea, oa))]
  dupO <- O[duplicated(akey) | duplicated(akey, fromLast = TRUE), unique(akey)]
  drops$dup_key_outcome <- O[akey %in% dupO, .N]; O <- O[!akey %in% dupO]
  ref <- R503[chr == cc$chr & pos37 >= cc$lo37 & pos37 <= cc$hi37]
  has_p <- isTRUE(cc$protein_local)
  if (has_p) {
    P <- fread(chk_extract(file.path(EXD, paste0(cc$cid, "__protein.tsv.gz"))), colClasses = list(character = c("ea", "oa")))
    jpos <- if (O$pos_build[1] == "GRCh37") "pos37" else "pos38"
    P[, row := .I]
    # outcome candidates at the protein position with compatible alleles (unique)
    m <- merge(P[, .(row, jp = get(jpos), ea, oa)], O[, .(jp = pos, oea = ea, ooa = oa, ob = b, ose = se, oeaf = eaf, orsid = rsid, on = n)],
               by = "jp", allow.cartesian = TRUE)
    m[, s := align_to_ref(oea, ooa, ea, oa)]
    m <- m[s != 0L]
    nm <- m[, .N, by = row]
    drops$nojoin_outcome <- nrow(P) - nrow(nm); drops$ambiguous_outcome <- nm[N > 1, .N]
    m <- m[row %in% nm[N == 1]$row]
    X <- merge(P, m[, .(row, s_o = s, ob, ose, oeaf, orsid, on)], by = "row")
    # reference at pos37, compatible alleles (unique)
    r <- merge(X[, .(row, pos37, ea, oa)], ref[, .(pos37, rsid_ref = rsid, a1, a2, af, nchrobs)], by = "pos37", allow.cartesian = TRUE)
    r[, sref := align_to_ref(ea, oa, a1, a2)]
    r <- r[sref != 0L]
    nr <- r[, .N, by = row]
    drops$ref_absent <- nrow(X) - nrow(nr); drops$ref_ambiguous <- nr[N > 1, .N]
    r <- r[row %in% nr[N == 1]$row]
    X <- merge(X, r[, .(row, rsid_ref, a1, a2, af, nchrobs, sref)], by = "row")
    X[, `:=`(eaf_P_a1 = ifelse(sref == 1L, eaf, 1 - eaf),
             eaf_O_a1 = { fo <- ifelse(s_o == 1L, oeaf, 1 - oeaf); ifelse(sref == 1L, fo, 1 - fo) },
             b_P = sref * b, se_P = se, b_O = s_o * sref * ob, se_O = ose)]
  } else {
    # no protein window: outcome-only set, outcome joined to the reference by rsID (build-agnostic) + compatible alleles
    O[, row := .I]
    ors <- O[, .(rs = unlist(strsplit(rsid, ",", fixed = TRUE))), by = row]
    r <- merge(ors, ref[, .(rs = rsid, pos37, rsid_ref = rsid, a1, a2, af, nchrobs)], by = "rs")
    r <- merge(r, O[, .(row, ea, oa)], by = "row")
    r[, sref := align_to_ref(ea, oa, a1, a2)]
    r <- r[sref != 0L]
    nr <- r[, .N, by = row]; nref <- r[, .N, by = rsid_ref]
    drops$ref_absent <- nrow(O) - nrow(nr); drops$ref_ambiguous <- nr[N > 1, .N] + nref[N > 1, .N]
    r <- r[row %in% nr[N == 1]$row & rsid_ref %in% nref[N == 1]$rsid_ref]
    X <- merge(O[, .(row, pos38 = pos, ea, oa, ob = b, ose = se, oeaf = eaf, orsid = rsid, on = n)], r[, .(row, pos37, rsid_ref, a1, a2, af, nchrobs, sref)], by = "row")
    X[, `:=`(eaf_O_a1 = ifelse(sref == 1L, oeaf, 1 - oeaf), b_O = sref * ob, se_O = ose, eaf_P_a1 = NA_real_, b_P = NA_real_, se_P = NA_real_, eaf = NA_real_)]
  }
  # QC on the 503 reference frequency (A1)
  X[, maf_ref := pmin(af, 1 - af)]
  X[, pal := is_pal(a1, a2)]
  X[, drop := NA_character_]
  X[is.na(drop) & (af == 0 | af == 1), drop := "monomorphic_ref"]
  X[is.na(drop) & maf_ref < REF_MAF_MIN, drop := "ref_maf_lt_0.01"]
  X[is.na(drop) & pal & maf_ref > PALIN_MAF_MAX, drop := "palindromic_maf_gt_0.40"]
  if (has_p) X[is.na(drop) & is.finite(eaf_P_a1) & abs(eaf_P_a1 - af) > AF_DIFF_MAX, drop := "af_discordant_protein"]
  X[is.na(drop) & is.finite(eaf_O_a1) & abs(eaf_O_a1 - af) > AF_DIFF_MAX, drop := "af_discordant_outcome"]
  for (k in unique(na.omit(X$drop))) drops[[k]] <- X[drop == k, .N]
  H <- X[is.na(drop)]
  setorder(H, pos37)
  H[, snp := sprintf("%s:%d:%s:%s", cc$chr, pos37, a1, a2)]
  stopifnot(!anyDuplicated(H$snp), !anyDuplicated(H$rsid_ref))
  H[, `:=`(z_P = b_P / se_P, z_O = b_O / se_O)]
  # p1np, p1fr
  H[, in_p1np := !pal]
  fr_ok <- function(f) is.na(f) | (pmin(f, 1 - f) <= PALIN_MAF_MAX)
  H[, in_p1fr := in_p1np | (pal & maf_ref <= PALIN_MAF_MAX &
                             (is.na(eaf_P_a1) | (fr_ok(eaf_P_a1) & abs(eaf_P_a1 - af) <= P1FR_AF_DIFF_MAX)) &
                             (is.na(eaf_O_a1) | (fr_ok(eaf_O_a1) & abs(eaf_O_a1 - af) <= P1FR_AF_DIFF_MAX)))]
  # p2np
  b2 <- R489[[cc$chr]][pos37 >= cc$lo37 & pos37 <= cc$hi37]
  j <- merge(H[in_p1np == TRUE, .(snp, pos37, a1, a2, af)], b2[, .(pos37, p2_rsid = rsid, p2_a1 = a1, p2_a2 = a2, af2_raw = af, nch2 = nchrobs)], by = "pos37", allow.cartesian = TRUE)
  j <- j[(p2_a1 == a1 & p2_a2 == a2) | (p2_a1 == a2 & p2_a2 == a1)]
  nj <- j[, .N, by = snp]
  d2 <- list(panel2_absent = H[in_p1np == TRUE, .N] - nrow(nj), panel2_ambiguous = nj[N > 1, .N])
  j <- j[snp %in% nj[N == 1]$snp]
  j[, p2_sign := ifelse(p2_a1 == a1, 1L, -1L)]
  j[, af2 := ifelse(p2_sign == 1L, af2_raw, 1 - af2_raw)]
  j[, drop2 := NA_character_]
  j[is.na(drop2) & nch2 != 2L * 489L, drop2 := "panel2_missing_genotype"]
  j[is.na(drop2) & (pmin(af2, 1 - af2) < REF_MAF_MIN), drop2 := "panel2_maf_lt_0.01"]
  j[is.na(drop2) & abs(af2 - af) > P2_AF_DIFF_MAX, drop2 := "panel2_af_discordant"]
  for (k in unique(na.omit(j$drop2))) d2[[k]] <- j[drop2 == k, .N]
  j <- j[is.na(drop2)]
  H <- merge(H, j[, .(snp, p2_rsid, p2_sign, af2)], by = "snp", all.x = TRUE, sort = FALSE)
  setorder(H, pos37)
  H[, in_p2np := !is.na(p2_sign)]
  # sample sizes
  n_O <- if (cc$outcome %in% c("eGFR")) median(O$n, na.rm = TRUE) else nout[outcome == cc$outcome]$neff
  n_O_src <- if (cc$outcome == "eGFR") "median n_total_sum over the cut region (engine rule for Wuttke eGFR)" else sprintf("N_eff from case/control counts (%s)", nout[outcome == cc$outcome]$source)
  n_P <- if (has_p) n_hat(H$se_P, H$eaf) else NA_real_
  st <- list(cid = cc$cid, gene = cc$gene, outcome = cc$outcome, chr = cc$chr, window38 = c(cc$lo38, cc$hi38), window37 = c(cc$lo37, cc$hi37),
             bounds37_source = cc$bounds37_source, n_outcome_rows = nrow(O), protein = has_p,
             n_susie = nrow(H), n_p1fr = H[in_p1fr == TRUE, .N], n_p1np = H[in_p1np == TRUE, .N], n_p2np = H[in_p2np == TRUE, .N],
             drops = drops, drops_p2 = d2, n_P = n_P, n_P_source = "n_hat = median 1/(2f(1-f)se^2), GWAS A1FREQ, MAF>0.05, over the QC'd SNPs (engine rule)",
             n_O = n_O, n_O_source = n_O_src,
             rsid_concordance_outcome = if ("orsid" %in% names(H)) sprintf("%d/%d", H[, sum(mapply(function(a, b) a %in% strsplit(b, ",", fixed = TRUE)[[1]], rsid_ref, orsid))], H[grepl("^rs", orsid), .N]) else "")
  keep <- c("snp", "rsid_ref", "chr", "pos37", "pos38", "a1", "a2", "af", "maf_ref", "pal", "b_P", "se_P", "z_P", "eaf_P_a1", "b_O", "se_O", "z_O", "eaf_O_a1",
            "in_p1fr", "in_p1np", "in_p2np", "p2_rsid", "p2_sign", "af2")
  H[, chr := cc$chr]
  if (!"pos38" %in% names(H)) H[, pos38 := NA_integer_]
  fwrite(H[, ..keep], file.path(HD, paste0(cc$cid, ".tsv.gz")), sep = "\t")
  writeLines(toJSON(st, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(HD, paste0(cc$cid, "__stamp.json")))
  say("harm %-12s susie %d | p1fr %d | p1np %d | p2np %d | n_P %s n_O %.0f | drops %s | p2 drops %s | rsID %s", cc$cid, st$n_susie, st$n_p1fr, st$n_p1np, st$n_p2np,
      if (has_p) sprintf("%.0f", n_P) else "NA", n_O, paste(names(drops), unlist(drops), sep = "=", collapse = " "),
      paste(names(d2), unlist(d2), sep = "=", collapse = " "), st$rsid_concordance_outcome)
  st
}
stamps <- lapply(seq_len(nrow(cand)), function(i) harmonise(cand[i]))
names(stamps) <- cand$cid

# ---- validation of the re-implementation against the engine's UMOD harmonised table (raw_input) and stamps (definitions) ----
eng <- fread(res("LOC_UMOD_eGFR_control_harmonised"))
mine <- fread(file.path(HD, "UMOD_eGFR.tsv.gz"))
mm <- merge(eng[, .(snp, bP = b__P__UMOD, sP = se__P__UMOD, bO = b__O__eGFR, sO = se__O__eGFR, af_ref)], mine[, .(snp, b_P, se_P, b_O, se_O, af)], by = "snp")
okh <- nrow(eng) == nrow(mine) && nrow(mm) == nrow(eng) && max(abs(mm$bP - mm$b_P)) < 1e-12 && max(abs(mm$bO - mm$b_O)) < 1e-12 &&
  max(abs(mm$sP - mm$se_P)) < 1e-12 && max(abs(mm$sO - mm$se_O)) < 1e-12 && max(abs(mm$af_ref - mm$af)) < 5e-4
log_assert("S6.d.reimpl_vs_engine_UMOD_harmonised", sprintf("%d SNPs, identical betas/SEs", nrow(eng)),
           sprintf("%d mine, %d joined, max|db_P| %.2g, max|db_O| %.2g, max|daf| %.2g", nrow(mine), nrow(mm), max(abs(mm$bP - mm$b_P)), max(abs(mm$bO - mm$b_O)), max(abs(mm$af_ref - mm$af))),
           okh, "R re-implementation of l3fm__prep.py QC vs LOC_UMOD_eGFR_control_harmonised")
ps <- fromJSON(res("LOC_UMOD_eGFR_control_PREP_STAMP")); p2s <- fromJSON(res("LOC_UMOD_eGFR_control_PANEL2_STAMP"))
su <- stamps[["UMOD_eGFR"]]
log_assert("S6.d.reimpl_UMOD_set_sizes", sprintf("susie %d / p1np %d / p2np %d (PREP/PANEL2 stamps, definitions)", ps$n_snps, p2s$n_p1np, p2s$n_p2np),
           sprintf("susie %d / p1np %d / p2np %d", su$n_susie, su$n_p1np, su$n_p2np),
           su$n_susie == ps$n_snps && su$n_p1np == p2s$n_p1np && su$n_p2np == p2s$n_p2np, "set sizes of the re-implementation; weak evidence (stamp values visible)")
saveRDS(stamps, file.path(HD, "stamps.rds"))
say("04_harmonise done")
