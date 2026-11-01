# w4b_03_omnibus.R -- Stage 1 W4b, step 3: section 4(iv) all-signal omnibus (GLS Q via the shared gls_q.R) and the
# per-signal supplementary table (Rule 3(iv): z and PP.H4 at each signal, each trait conditioned on the OTHER signals).
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets/w4b_common.R"))
suppressMessages({library(coloc); library(readxl)})
logf("==== W4b step 3 (omnibus, per-signal) start")
assert_md5(GOV$GLS_Q_R$path, GOV$GLS_Q_R$md5, "gls_q.R (sourced)")
source(GOV$GLS_Q_R$path)
D <- readRDS(file.path(OUT, "data", "w4b_data.rds")); H <- D$H
acc <- new_rows()

# ---------------------------------------------------------------- signal sets
sig <- fread(rpath("LOC_ARG1_6q23_signals"), colClasses = list(character = c("snp", "rsid")))
S5 <- sig[tid == "P__ARG1"][order(signal)]
stopifnot(nrow(S5) == 5L, all(S5$snp %in% H$snp))
logf("5-signal set (LOC_ARG1_6q23_signals, P__ARG1): %s", paste(sprintf("%d=%s(%s)", S5$signal, S5$rsid, S5$snp), collapse = ", "))
SIGA <- "6:131897278:T:C"; SIGB <- "6:131894358:T:C"      # R1_ARG1_rule: A = rs2781668, B = rs17788484
stopifnot(S5$snp[1] == SIGA, S5$snp[2] == SIGB)

st16 <- read_excel(rpath("EXP_SUN2023_SUPP_TABLES_ST16"), sheet = "ST16", skip = 4)
names(st16)[1] <- "protein_id"
x6 <- as.data.table(st16)[protein_id == "ARG1:P05089:OID21445:v1" & `Cis/trans` == "cis"]
stopifnot(nrow(x6) == 6L)
x6[, c("c37", "p37", "ref", "alt") := tstrsplit(`Variant ID`, ":", keep = 1:4)][, p37 := as.integer(p37)]
x6 <- merge(x6, H[, .(pos37, snp, a1, a2, rsid_ref)], by.x = "p37", by.y = "pos37", all.x = TRUE)
x6 <- x6[is.na(snp) | (alt == a1 & ref == a2) | (alt == a2 & ref == a1)]
stopifnot(nrow(x6) == 6L, !anyNA(x6$snp), all(x6$rsID == x6$rsid_ref))
setorder(x6, -`-log10(P)(cond)`)
logf("6-signal set (ST16 ARG1 cis rows, UKB-PPP in-sample): %s", paste(sprintf("%s(%s, PIP %.3f)", x6$rsID, x6$snp, as.numeric(x6$PIP)), collapse = ", "))
S6 <- x6$snp

# ---------------------------------------------------------------- GLS omnibus
OMNI_TRAITS <- list(arginine = list(b = "b__M__GCST90200372", se = "se__M__GCST90200372", res = "CLSA GCST90200372"),
                    T2D_FG = list(b = "b__FG", se = "se__FG", res = "FinnGen R10 T2D"),
                    T2D_DM = list(b = "b__DM", se = "se__DM", res = "DIAMANTE 2018 EUR noUKB T2D"))
run_omni <- function(tr, setname, snps, note = "") {
  h <- H[match(snps, snp)]
  bY <- h[[OMNI_TRAITS[[tr]]$b]]; sY <- h[[OMNI_TRAITS[[tr]]$se]]
  if (anyNA(bY)) stop("missing trait value at a signal: ", tr, " ", setname)
  R <- D$R503[snps, snps]
  jp <- joint_effects(h$b__P__ARG1, h$se__P__ARG1, R); jy <- joint_effects(bY, sY, R)
  cells <- list(aware = c("aware", "2nd"), ign = c("ign", "2nd"), aware_1st = c("aware", "1st"))
  out <- list()
  for (cn in names(cells)) {
    q <- profile_Q(jy$b, jy$V, jp$b, jp$V, cov = cells[[cn]][1], se_order = cells[[cn]][2])
    if (q$at_bound) logf("WARNING theta at the search bound: %s %s %s", tr, setname, cn)
    out[[cn]] <- q
    est <- sprintf("GLS Q of %s-per-SD-ARG1 joint ratios across the %d-signal set (%s); joint effects from marginal z + 503 LD (PLINK)",
                   tr, length(snps), setname)
    base <- sprintf("S4.omni.%s.%s.%s", tr, setname, cn)
    nt <- paste0(sprintf("cov=%s se_order=%s; signals %s; theta at bound: %s", cells[[cn]][1], cells[[cn]][2],
                         paste(h$rsid_ref, collapse = ","), q$at_bound), if (nzchar(note)) paste0("; ", note) else "")
    for (f in c("Q", "df", "p_upper", "p_lower", "theta")) {
      acc <<- add_row(acc, paste0(base, ".", f), "4(iv)", est, OMNI_TRAITS[[tr]]$res, "503", cells[[cn]][1], cells[[cn]][2],
                      value = q[[f]], p = if (f %in% c("p_upper", "p_lower")) q[[f]] else NA, n = length(snps),
                      unit = switch(f, Q = "chi-square", df = "df", p_upper = "p", p_lower = "p", theta = "trait units per SD ARG1"), notes = nt)
    }
    logf("OMNI %-8s %-14s %-9s Q=%.4f df=%d p_upper=%.4g p_lower=%.4g theta=%.5f", tr, setname, cn, q$Q, q$df, q$p_upper, q$p_lower, q$theta)
  }
  # Rule 5 ids for a GLS Q
  ids <- c(p_aware_2nd = out$aware$p_upper, p_ign_2nd = out$ign$p_upper, p_aware_1st = out$aware_1st$p_upper)
  for (k in names(ids)) acc <<- add_row(acc, sprintf("S4.omni.%s.%s.%s", tr, setname, k), "4(iv); Rule 5", "upper-tail p of the GLS Q",
                                         OMNI_TRAITS[[tr]]$res, "503", sub("^p_([a-z]+)_.*$", "\\1", k), sub("^.*_", "", k),
                                         value = ids[[k]], p = ids[[k]], n = length(snps), unit = "p",
                                         notes = if (k == "p_aware_2nd") "primary (Rule 5)" else "printed beside the primary (Rule 5)")
  acc <<- add_row(acc, sprintf("S4.omni.%s.%s.p_cojo", tr, setname), "4(iv); Rule 5", "beta-space (COJO) p", OMNI_TRAITS[[tr]]$res,
                  "503", "", "", value = "NOT RUN for GLS Q", notes = "Rule 5: p_cojo applies to two-signal tests only")
  invisible(out)
}
OM <- list()
for (tr in names(OMNI_TRAITS)) OM[[paste(tr, "5sig")]] <- run_omni(tr, "5sig", S5$snp)
# six-signal sensitivity (UKB-PPP in-sample set)
for (tr in c("arginine", "T2D_DM")) OM[[paste(tr, "6sig")]] <- run_omni(tr, "6sig", S6,
   note = "rs2698107 (G/C palindrome, ref MAF 0.379) kept: its allele frequency agrees in every source used (p1fr frequency rule)")
miss6 <- S6[is.na(H[match(S6, snp)]$b__FG)]
logf("6-signal set: FinnGen R10 T2D lacks %s", paste(H[match(miss6, snp)]$rsid_ref, collapse = ","))
acc <- add_row(acc, "S4.omni.T2D_FG.6sig.status", "4(iv)", "six-signal FinnGen T2D GLS", "FinnGen R10 T2D", "503",
               value = sprintf("DEVIATION: %s (UKB-PPP signal %s) is absent from FinnGen R10 (no row at its GRCh38 position); run with a proxy (6sigproxy) and with the signal dropped (6sigdrop5)",
                               paste(H[match(miss6, snp)]$rsid_ref, collapse = ","), paste(miss6, collapse = ",")))
if (length(miss6)) {
  # proxy = the highest-PIP member of that signal's UKB-PPP credible set that is in the harmonised table and in FinnGen
  r6 <- x6[snp %in% miss6]
  cs_ids <- strsplit(r6$`CS variant IDs`, ", ")[[1]]; cs_pip <- as.numeric(strsplit(r6$`CS PIPs`, ", ")[[1]])
  cs <- data.table(id = cs_ids, pip = cs_pip)[, c("p37", "ref", "alt") := tstrsplit(id, ":", keep = 2:4)][, p37 := as.integer(p37)]
  cs <- merge(cs, H[!is.na(b__FG), .(pos37, snp, a1, a2, rsid_ref)], by.x = "p37", by.y = "pos37")
  cs <- cs[(alt == a1 & ref == a2) | (alt == a2 & ref == a1)][order(-pip)]
  stopifnot(nrow(cs) >= 1)
  px <- cs$snp[1]
  r2px <- D$R503[px, miss6]^2
  logf("FinnGen 6-signal proxy for %s: %s (%s), UKB-PPP CS PIP %.3f, r2 with the signal on 503 = %.3f", miss6, cs$rsid_ref[1], px, cs$pip[1], r2px)
  S6p <- replace(S6, S6 == miss6, px)
  OM[["T2D_FG 6sigproxy"]] <- run_omni("T2D_FG", "6sigproxy", S6p,
     note = sprintf("DEVIATION: %s absent from FinnGen R10; proxy %s (highest-PIP member of its UKB-PPP credible set present in FinnGen; r2 %.3f on 503); protein side uses the same proxy", H[match(miss6, snp)]$rsid_ref, cs$rsid_ref[1], r2px))
  OM[["T2D_FG 6sigdrop5"]] <- run_omni("T2D_FG", "6sigdrop5", setdiff(S6, miss6),
     note = sprintf("DEVIATION: %s absent from FinnGen R10; the remaining five UKB-PPP signals", H[match(miss6, snp)]$rsid_ref))
}

# ---------------------------------------------------------------- conditioning (z-space; Yang 2012 approximation; engine collinearity 0.9)
condition <- function(z, se, R, S) {
  RSS <- R[S, S, drop = FALSE]; RiS <- R[, S, drop = FALSE]
  W <- solve(RSS, t(RiS)); q <- rowSums(RiS * t(W))
  zc <- (z - RiS %*% solve(RSS, z[S])) / sqrt(pmax(1 - q, 1e-12)); sec <- se / sqrt(pmax(1 - q, 1e-12))
  bad <- q > 0.9; bad[S] <- TRUE
  zc[bad] <- NA; sec[bad] <- NA
  list(z = as.vector(zc), se = as.vector(sec))
}
coloc_cond <- function(snps, R, X, others, trait_b, trait_se, trait_type, casefrac) {
  h <- H[match(snps, snp)]
  S <- match(others, snps); stopifnot(!anyNA(S))
  cp <- condition(h$z__P__ARG1, h$se__P__ARG1, R, S)
  zt <- h[[trait_b]] / h[[trait_se]]
  ct <- condition(zt, h[[trait_se]], R, S)
  ix <- which(!is.na(cp$z) & !is.na(ct$z))
  d1 <- list(beta = cp$z[ix] * cp$se[ix], varbeta = cp$se[ix]^2, type = "quant", sdY = 1, snp = snps[ix])
  d2 <- if (trait_type == "cc") list(beta = ct$z[ix] * ct$se[ix], varbeta = ct$se[ix]^2, type = "cc", s = casefrac, snp = snps[ix]) else
    list(beta = ct$z[ix] * ct$se[ix], varbeta = ct$se[ix]^2, type = "quant", sdY = 1, snp = snps[ix])
  ab <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  jx <- match(X, snps)
  list(z_trait = ct$z[jx], z_prot = cp$z[jx], PPH4 = unname(ab$summary["PP.H4.abf"]), PPH3 = unname(ab$summary["PP.H3.abf"]),
       nsnp = length(ix), top = ab$results$snp[which.max(ab$results$SNP.PP.H4)])
}
PS_TRAITS <- list(arginine = list(b = "b__M__GCST90200372", se = "se__M__GCST90200372", type = "quant", cf = NA, res = "CLSA GCST90200372"),
                  T2D_FG = list(b = "b__FG", se = "se__FG", type = "cc", cf = D$CASEFRAC$FG, res = "FinnGen R10 T2D"),
                  T2D_DM = list(b = "b__DM", se = "se__DM", type = "cc", cf = D$CASEFRAC$DM, res = "DIAMANTE 2018 EUR noUKB T2D"))
PSROWS <- list()
for (pn in c("503", "489")) {
  setn <- if (pn == "503") "p1np" else "p2np"
  R <- if (pn == "503") D$R503 else D$R489
  for (tr in names(PS_TRAITS)) {
    snps <- D$SETS[[setn]]; snps <- snps[!is.na(H[match(snps, snp)][[PS_TRAITS[[tr]]$b]])]
    stopifnot(all(S5$snp %in% snps))
    Rs <- R[snps, snps]
    for (k in seq_len(nrow(S5))) {
      X <- S5$snp[k]; oth <- setdiff(S5$snp, X)
      r <- coloc_cond(snps, Rs, X, oth, PS_TRAITS[[tr]]$b, PS_TRAITS[[tr]]$se, PS_TRAITS[[tr]]$type, PS_TRAITS[[tr]]$cf)
      PSROWS[[length(PSROWS) + 1]] <- data.table(trait = tr, panel = pn, set = setn, signal = S5$signal[k], rsid = S5$rsid[k], snp = X,
                                                 z_trait_cond = r$z_trait, z_protein_cond = r$z_prot, PPH4 = r$PPH4, PPH3 = r$PPH3,
                                                 nsnp = r$nsnp, top_shared = H[match(r$top, snp)]$rsid_ref)
      est <- sprintf("signal %d (%s): %s and ARG1 each conditioned on the other four protein signals (z-space, %s LD); coloc.abf p1=p2=1e-4, p12=1e-5",
                     S5$signal[k], S5$rsid[k], tr, pn)
      acc <- add_row(acc, sprintf("S4.persig.%s.%d.%s.z", tr, S5$signal[k], pn), "4(iv); Rule 3(iv) supplementary table", est,
                     PS_TRAITS[[tr]]$res, pn, value = r$z_trait, n = r$nsnp, unit = "z (conditional)",
                     notes = sprintf("SNP set %s; protein conditional z %.3f", setn, r$z_prot))
      acc <- add_row(acc, sprintf("S4.persig.%s.%d.%s.PPH4", tr, S5$signal[k], pn), "4(iv); Rule 3(iv) supplementary table", est,
                     PS_TRAITS[[tr]]$res, pn, value = r$PPH4, n = r$nsnp, unit = "posterior",
                     notes = sprintf("SNP set %s; PP.H3 %.4f; top shared SNP %s", setn, r$PPH3, H[match(r$top, snp)]$rsid_ref))
      logf("PERSIG %-8s panel %s signal %d %-11s z_trait|others %7.3f  z_prot|others %7.3f  PP.H4 %.4f PP.H3 %.4f nsnp %d top %s", tr, pn,
           S5$signal[k], S5$rsid[k], r$z_trait, r$z_prot, r$PPH4, r$PPH3, r$nsnp, H[match(r$top, snp)]$rsid_ref)
    }
  }
}
PS <- rbindlist(PSROWS); fwrite(PS, file.path(OUT, "persignal_table.tsv"), sep = "\t")
stopifnot(nrow(PS) == 2 * 3 * 5, !anyNA(PS$z_trait_cond), !anyNA(PS$PPH4))
write_part(acc, "03_omnibus")
saveRDS(list(S5 = S5, S6 = S6, OM = OM), file.path(OUT, "data", "w4b_omnibus.rds"))
logf("==== W4b step 3 done")
