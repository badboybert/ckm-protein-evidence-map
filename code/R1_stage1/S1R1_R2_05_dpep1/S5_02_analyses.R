# S5_02_analyses.R -- W5 (§5 DPEP1) step 2: GLS Q (cysteinylglycine, eGFR), the two-signal re-test (Rule 5 ids), the
# Rule 3(iv) colocalization-clause inputs, the per-signal table, s over four SNP sets x two panels, the separation r^2 and
# the Rule 3(iii) gene-naming counts. Bootstrap is in S5_03. Output: work/analyses.rds, work/res_02.tsv.
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_05_dpep1/S5_common.R"))
t0 <- Sys.time()
say("==== S5_02_analyses start ====")
assert_governing()
assert_md5(GOV$GLS_Q$p, GOV$GLS_Q$m, "gls_q.R (sourced)")
source(GOV$GLS_Q$p)   # joint_effects(), profile_Q() -- the one shared implementation
P <- readRDS(file.path(WORK, "prep.rds"))
H <- P$H; R503 <- P$R503; R489 <- P$R489; EGFR <- P$EGFR
rs <- H$rsid_ref
ix <- function(x) match(x, rs)
RES <- list()
addr <- function(quantity_id, prespec_ref, estimand, resource, panel, cov = "", se_order = "", value = NA, se = NA, ci_lo = NA,
                 ci_hi = NA, p = NA, n = NA, unit = "", notes = "") {
  RES[[length(RES) + 1]] <<- data.table(quantity_id, prespec_ref, estimand, resource, panel, cov, se_order, value = as.character(value),
                                         se = as.character(se), ci_lo = as.character(ci_lo), ci_hi = as.character(ci_hi),
                                         p = as.character(p), n = as.character(n), unit, notes)
}
num <- function(x) if (is.numeric(x)) formatC(x, digits = 8, format = "g") else as.character(x)

TR <- list(
  cysgly = list(b = H$b__M__GCST90200364, se = H$se__M__GCST90200364, label = "cysteinylglycine (CLSA GCST90200364)", unit = "metabolite units per allele"),
  eGFR   = list(b = EGFR$b, se = EGFR$se, label = "eGFR (CKDGen 2019 EA, ln eGFR)", unit = "ln(eGFR) per allele"))
PRT <- list(b = H$b__P__DPEP1, se = H$se__P__DPEP1)
PANELS <- list(`503` = R503, `489` = R489)
sig <- ix(SIG5)
stopifnot(!anyNA(sig), all(SIG5 %in% rownames(R489)))

# ================================================================ A. GLS Q (pre-declaration §5; Rule 5 GLS ids)
COVS <- list(aware_2nd = c("aware", "2nd"), ign_2nd = c("ign", "2nd"), aware_1st = c("aware", "1st"))
gls <- list()
for (tn in names(TR)) for (pn in names(PANELS)) {
  Rs <- PANELS[[pn]][SIG5, SIG5]
  jp <- joint_effects(PRT$b[sig], PRT$se[sig], Rs)
  jy <- joint_effects(TR[[tn]]$b[sig], TR[[tn]]$se[sig], Rs)
  # 3-SNP-only joint model for signals 1-3 (a labelled sensitivity of the secondary set; the declared secondary uses idx 1:3
  # of the five-signal joint model, i.e. signals 1-3 conditioned on all five)
  jp3 <- joint_effects(PRT$b[sig[1:3]], PRT$se[sig[1:3]], Rs[1:3, 1:3]); jy3 <- joint_effects(TR[[tn]]$b[sig[1:3]], TR[[tn]]$se[sig[1:3]], Rs[1:3, 1:3])
  for (st in c("all5", "sig123", "sig123_joint3only")) for (cn in names(COVS)) {
    o <- switch(st,
      all5 = profile_Q(jy$b, jy$V, jp$b, jp$V, idx = 1:5, cov = COVS[[cn]][1], se_order = COVS[[cn]][2]),
      sig123 = profile_Q(jy$b, jy$V, jp$b, jp$V, idx = 1:3, cov = COVS[[cn]][1], se_order = COVS[[cn]][2]),
      sig123_joint3only = profile_Q(jy3$b, jy3$V, jp3$b, jp3$V, idx = 1:3, cov = COVS[[cn]][1], se_order = COVS[[cn]][2]))
    gls[[paste(tn, st, pn, cn, sep = ".")]] <- o
    id <- sprintf("S5.Q.%s.%s.%s.%s", tn, st, pn, cn)
    est <- sprintf("GLS proportionality Q of %s-per-DPEP1 joint ratios across %s; joint (conditional) per-allele effects from marginal z with %s-panel PLINK r among the 5 declared signals (zJ = R^-1 z, bJ = zJ*se); %s",
                   TR[[tn]]$label, switch(st, all5 = "all 5 declared DPEP1 signals (primary set)", sig123 = "signals 1-3 (secondary set; idx 1:3 of the five-signal joint model, i.e. conditioned on all five)",
                                          sig123_joint3only = "signals 1-3 from a 3-SNP joint model (labelled sensitivity, not declared)"),
                   pn, if (cn == "aware_2nd") "Rule 5 PRIMARY (aware, 2nd order)" else cn)
    prim <- st == "all5" && pn == "503" && cn == "aware_2nd"
    nt <- sprintf("%s%s; theta at_bound=%s", if (prim) "Rule 3(ii) primary cell; Rule 4 family member (Holm in W9); " else "", if (st == "sig123_joint3only") "DEVIATION-free extra sensitivity" else "", o$at_bound)
    res_lab <- sprintf("UKB-PPP OID20897 x %s", if (tn == "cysgly") "CLSA GCST90200364" else "CKDGen 2019 EA eGFR (a995e47a)")
    addr(paste0(id, ".Q"), "§5 Analyses; Rule 5", est, res_lab, pn, COVS[[cn]][1], COVS[[cn]][2], num(o$Q), unit = "chi-square", notes = nt)
    addr(paste0(id, ".df"), "§5 Analyses; Rule 5", est, res_lab, pn, COVS[[cn]][1], COVS[[cn]][2], o$df, unit = "df")
    addr(paste0(id, ".p"), "§5 Analyses; Rule 5", est, res_lab, pn, COVS[[cn]][1], COVS[[cn]][2], num(o$p_upper), p = num(o$p_upper), unit = "p (upper tail)",
         notes = if (cn == "aware_2nd") "p_aware_2nd (primary)" else if (cn == "ign_2nd") "p_ign_2nd" else "p_aware_1st")
    addr(paste0(id, ".p_lower"), "§5 Analyses (printed beside)", est, res_lab, pn, COVS[[cn]][1], COVS[[cn]][2], num(o$p_lower), p = num(o$p_lower), unit = "p (lower tail)")
    addr(paste0(id, ".theta"), "§5 Analyses", est, res_lab, pn, COVS[[cn]][1], COVS[[cn]][2], num(o$theta), unit = paste0(TR[[tn]]$unit, " per SD DPEP1 (common ratio)"))
    if (cn == "aware_2nd") addr(paste0(sprintf("S5.Q.%s.%s.%s", tn, st, pn), ".p_cojo"), "Rule 5", est, res_lab, pn, "", "", "NOT RUN for GLS Q",
                                notes = "Rule 5: p_cojo NOT RUN for GLS Q (gls_q.R has no beta-space option)")
  }
  say("GLS %s %s: all5 aware_2nd Q %.3f df 4 p %.4g (lower %.3g) theta %.5g | ign_2nd p %.4g | aware_1st p %.4g | sig123 aware_2nd p %.4g",
      tn, pn, gls[[paste(tn, "all5", pn, "aware_2nd", sep = ".")]]$Q, gls[[paste(tn, "all5", pn, "aware_2nd", sep = ".")]]$p_upper,
      gls[[paste(tn, "all5", pn, "aware_2nd", sep = ".")]]$p_lower, gls[[paste(tn, "all5", pn, "aware_2nd", sep = ".")]]$theta,
      gls[[paste(tn, "all5", pn, "ign_2nd", sep = ".")]]$p_upper, gls[[paste(tn, "all5", pn, "aware_1st", sep = ".")]]$p_upper,
      gls[[paste(tn, "sig123", pn, "aware_2nd", sep = ".")]]$p_upper)
}

# ================================================================ B. two-signal re-test, rs409170 (A) vs rs62068372 (B), cysteinylglycine
ab <- ix(c(SIG_A, SIG_B))
# ratio_diff() (the two-signal A-vs-B test, Rule 5 ids) lives in S5_common.R so S5_03 uses the same implementation
# beta-space (COJO, Yang et al. 2012) joint estimates for the two SNPs: D = diag(2pq N), B = D^1/2 R D^1/2,
# bJ = B^-1 D b, sigma2_J = ((N-1) Vp - bJ' D b) / (N - k), V = sigma2_J B^-1; Vp = median over window SNPs of 2pq((N-1)se^2 + b^2).
cojo_joint <- function(b, se, f, R, N, bw, sew, fw) {
  h <- 2 * f * (1 - f); D <- diag(h * N, nrow = length(b))
  okw <- is.finite(fw) & pmin(fw, 1 - fw) > 0.01
  Vp <- median(2 * fw[okw] * (1 - fw[okw]) * ((N - 1) * sew[okw]^2 + bw[okw]^2))
  Bm <- sqrt(D) %*% R %*% sqrt(D)
  bJ <- as.vector(solve(Bm, D %*% b))
  s2 <- ((N - 1) * Vp - sum(bJ * (D %*% b))) / (N - length(b))
  list(b = bJ, V = s2 * solve(Bm), Vp = Vp)
}
two <- list()
for (pn in names(PANELS)) {
  R2 <- PANELS[[pn]][c(SIG_A, SIG_B), c(SIG_A, SIG_B)]
  jp <- joint_effects(PRT$b[ab], PRT$se[ab], R2); jm <- joint_effects(TR$cysgly$b[ab], TR$cysgly$se[ab], R2)
  a2 <- ratio_diff(jp, jm, "aware", "2nd"); ig <- ratio_diff(jp, jm, "ign", "2nd"); a1 <- ratio_diff(jp, jm, "aware", "1st")
  cp <- cojo_joint(PRT$b[ab], PRT$se[ab], P$frq$DPEP1[ab], R2, P$N_PROT, PRT$b, PRT$se, P$frq$DPEP1)
  cm <- cojo_joint(TR$cysgly$b[ab], TR$cysgly$se[ab], P$frq$CYS[ab], R2, P$N_CYS, TR$cysgly$b, TR$cysgly$se, P$frq$CYS)
  cj <- ratio_diff(cp, cm, "aware", "2nd")
  two[[pn]] <- list(r_AB = R2[1, 2], aware_2nd = a2, ign_2nd = ig, aware_1st = a1, cojo = cj)
  est <- sprintf("cysteinylglycine-per-DPEP1 ratio at signal A (rs409170 | rs62068372) vs signal B (rs62068372 | rs409170); two-SNP joint (= single-SNP conditional) estimates, %s-panel PLINK r(A,B) = %.5f; test of kA - kB, delta method", pn, R2[1, 2])
  rl <- "UKB-PPP OID20897 x CLSA GCST90200364"
  addr(sprintf("S5.two.%s.r_AB", pn), "§5 two-signal re-test", "PLINK r(rs409170, rs62068372), aligned to harmonised A1", paste0("1000G EUR ", pn), pn, value = num(R2[1, 2]), unit = "r")
  addr(sprintf("S5.two.%s.kA", pn), "§5 two-signal re-test", est, rl, pn, "", "2nd", num(a2$kA), se = num(a2$seA), unit = "cys-gly units per SD DPEP1")
  addr(sprintf("S5.two.%s.kB", pn), "§5 two-signal re-test", est, rl, pn, "", "2nd", num(a2$kB), se = num(a2$seB), unit = "cys-gly units per SD DPEP1")
  addr(sprintf("S5.two.%s.diff", pn), "§5 two-signal re-test", est, rl, pn, "aware", "2nd", num(a2$d), se = num(a2$se_d), ci_lo = num(a2$d - 1.959964 * a2$se_d), ci_hi = num(a2$d + 1.959964 * a2$se_d), p = num(a2$p), unit = "kA - kB")
  addr(sprintf("S5.two.%s.p_aware_2nd", pn), "§5; Rule 5", est, rl, pn, "aware", "2nd", num(a2$p), se = num(a2$se_d), p = num(a2$p), unit = "p", notes = sprintf("PRIMARY; z %.4f; sign(kA-kB) %+d", a2$z, sign(a2$d)))
  addr(sprintf("S5.two.%s.p_ign", pn), "§5; Rule 5", est, rl, pn, "ign", "2nd", num(ig$p), se = num(ig$se_d), p = num(ig$p), unit = "p", notes = sprintf("z %.4f", ig$z))
  addr(sprintf("S5.two.%s.p_aware_1st", pn), "§5; Rule 5", est, rl, pn, "aware", "1st", num(a1$p), se = num(a1$se_d), p = num(a1$p), unit = "p", notes = sprintf("z %.4f", a1$z))
  addr(sprintf("S5.two.%s.p_cojo", pn), "§5; Rule 5", paste0(est, "; joint estimates in beta space (COJO: D = 2pqN with GWAS A1 frequencies, Vp = window median)"), rl, pn, "aware", "2nd",
       num(cj$p), se = num(cj$se_d), p = num(cj$p), unit = "p", notes = sprintf("kA %.5f kB %.5f z %.4f; N protein %d, cys-gly %d", cj$kA, cj$kB, cj$z, P$N_PROT, P$N_CYS))
  say("TWO-SIGNAL %s: r_AB %.5f; PA %.4f PB %.4f MA %.4f MB %.4f; kA %.5f (se %.5f) kB %.5f (se %.5f); diff %.5f se %.5f; p_aware_2nd %.4g, p_ign %.4g, p_aware_1st %.4g, p_cojo %.4g",
      pn, R2[1, 2], a2$PA, a2$PB, a2$MA, a2$MB, a2$kA, a2$seA, a2$kB, a2$seB, a2$d, a2$se_d, a2$p, ig$p, a1$p, cj$p)
}

# ================================================================ C. conditional colocalization (Rule 3(iv)) and per-signal table
set_rows <- list(`503` = P$SETS$p1np, `489` = P$SETS$p2np, `503full` = P$SETS$susie)
set_R <- list(`503` = R503, `489` = R489, `503full` = R503)
set_lab <- list(`503` = "p1np (palindromes dropped), 503-panel LD", `489` = "p2np, 489-panel LD", `503full` = "all harmonised SNPs (palindromes resolvable by frequency kept), 503-panel LD [sensitivity]")
maf_ref <- pmin(H$af_ref, 1 - H$af_ref)
cond_set <- function(trait_z, trait_se, setr, Rm, condS) {
  ids <- rs[setr]; ids <- ids[ids %in% rownames(Rm)]
  RiS <- Rm[ids, condS, drop = FALSE]; RSS <- Rm[condS, condS, drop = FALSE]
  zS <- trait_z[ix(condS)]
  rr <- ix(ids)
  cz <- condition_z(trait_z[rr], trait_se[rr], RiS, RSS, zS, isS = ids %in% condS)
  list(ids = ids, rr = rr, z = cz$z, se = cz$se, ncol = cz$n_collinear)
}
zP <- PRT$b / PRT$se
TRz <- lapply(TR, function(t) t$b / t$se)
coloc_cond <- function(tn, setk, condS) {
  Rm <- set_R[[setk]]
  cp <- cond_set(zP, PRT$se, set_rows[[setk]], Rm, condS)
  ct <- cond_set(TRz[[tn]], TR[[tn]]$se, set_rows[[setk]], Rm, condS)
  stopifnot(identical(cp$ids, ct$ids))
  extra <- if (tn == "eGFR") list(N = P$N_EGFR_MED, MAF = maf_ref[ct$rr]) else NULL
  r <- coloc_pair(cp$z * cp$se, cp$se, ct$z * ct$se, ct$se, cp$ids, ds2_extra = extra)
  r$n_collinear <- cp$ncol; r
}
# Rule 3(iv) clause for the two-signal re-test: trait conditional z at A|B and B|A, conditional PP.H4 (protein|other vs trait|other)
clause <- list()
for (pn in c("503", "489")) {
  R2 <- PANELS[[pn]][c(SIG_A, SIG_B), c(SIG_A, SIG_B)]
  jm <- joint_effects(TR$cysgly$b[ab], TR$cysgly$se[ab], R2)
  zc <- jm$b / sqrt(diag(jm$V))
  cA <- coloc_cond("cysgly", pn, SIG_B); cB <- coloc_cond("cysgly", pn, SIG_A)
  clause[[pn]] <- list(zA = zc[1], zB = zc[2], H4A = cA$H4, H4B = cB$H4, H3A = cA$H3, H3B = cB$H3, nA = cA$nsnp, nB = cB$nsnp, leadA = cA$lead, leadB = cB$lead)
  for (sx in c("A", "B")) {
    sg <- if (sx == "A") SIG_A else SIG_B; oth <- if (sx == "A") SIG_B else SIG_A
    cc <- if (sx == "A") cA else cB
    zz <- if (sx == "A") zc[1] else zc[2]
    addr(sprintf("S5.coloc_clause.%s_%s.%s.z", sx, sg, pn), "Rule 3(iv); §5", sprintf("cysteinylglycine conditional z at %s given %s (two-SNP joint, %s-panel r)", sg, oth, pn),
         "CLSA GCST90200364", pn, value = num(zz), unit = "z")
    addr(sprintf("S5.coloc_clause.%s_%s.%s.PPH4", sx, sg, pn), "Rule 3(iv); §5", sprintf("coloc.abf PP.H4, DPEP1 protein | %s vs cysteinylglycine | %s, over %s; p1 = p2 = 1e-4, p12 = 1e-5, prior sd 0.15 (sdY = 1)", oth, oth, set_lab[[pn]]),
         "UKB-PPP OID20897 x CLSA GCST90200364", pn, value = num(cc$H4), n = cc$nsnp, unit = "PP.H4",
         notes = sprintf("PP.H3 %.4g; lead SNP %s; %d SNPs; %d collinear (multiple r^2 > 0.9) dropped", cc$H3, cc$lead, cc$nsnp, cc$n_collinear))
  }
  say("COLOC CLAUSE %s: z_cysgly(A|B) %.3f PP.H4_A %.4f (H3 %.3g, lead %s, n %d) | z_cysgly(B|A) %.3f PP.H4_B %.4f (H3 %.3g, lead %s, n %d)",
      pn, zc[1], cA$H4, cA$H3, cA$lead, cA$nsnp, zc[2], cB$H4, cB$H3, cB$lead, cB$nsnp)
}
# per-signal table for the all-signal GLS (Rule 3(iv): supplementary table at every locus alike)
persig <- list()
for (pn in c("503", "489")) {
  Rs <- PANELS[[pn]][SIG5, SIG5]
  jp <- joint_effects(PRT$b[sig], PRT$se[sig], Rs)
  for (tn in names(TR)) {
    jy <- joint_effects(TR[[tn]]$b[sig], TR[[tn]]$se[sig], Rs)
    for (k in seq_along(SIG5)) {
      th <- jy$b[k] / jp$b[k]; thse <- sqrt(jy$V[k, k] / jp$b[k]^2 + jy$b[k]^2 * jp$V[k, k] / jp$b[k]^4)
      id <- sprintf("S5.persig.%s.s%d_%s.%s", tn, k, SIG5[k], pn)
      addr(paste0(id, ".z_marg_trait"), "Rule 3(iv) per-signal table", sprintf("%s marginal z at %s", TR[[tn]]$label, SIG5[k]), TR[[tn]]$label, pn, value = num(TRz[[tn]][sig[k]]), unit = "z")
      addr(paste0(id, ".z_joint_trait"), "Rule 3(iv) per-signal table", sprintf("%s z at %s conditioned on the other four signals (five-SNP joint, %s-panel r)", TR[[tn]]$label, SIG5[k], pn), TR[[tn]]$label, pn, value = num(jy$b[k] / sqrt(jy$V[k, k])), unit = "z")
      addr(paste0(id, ".ratio"), "Rule 3(iv) per-signal table", sprintf("joint %s-per-DPEP1 Wald ratio at %s (five-SNP joint), 2nd-order delta SE", TR[[tn]]$label, SIG5[k]), TR[[tn]]$label, pn, "", "2nd",
           num(th), se = num(thse), ci_lo = num(th - 1.959964 * thse), ci_hi = num(th + 1.959964 * thse), unit = paste(TR[[tn]]$unit, "per SD DPEP1"))
      if (tn == "cysgly") addr(sprintf("S5.persig.protein.s%d_%s.%s.z_joint", k, SIG5[k], pn), "Rule 3(iv) per-signal table", sprintf("DPEP1 protein z at %s conditioned on the other four signals", SIG5[k]), "UKB-PPP OID20897", pn, value = num(jp$b[k] / sqrt(jp$V[k, k])), unit = "z")
      persig[[paste(tn, k, pn)]] <- data.table(trait = tn, signal = k, rsid = SIG5[k], panel = pn, z_marg = TRz[[tn]][sig[k]], z_joint = jy$b[k] / sqrt(jy$V[k, k]),
                                               zP_joint = jp$b[k] / sqrt(jp$V[k, k]), ratio = th, ratio_se = thse)
    }
  }
}
for (setk in c("503", "489", "503full")) for (tn in names(TR)) for (k in seq_along(SIG5)) {
  cc <- coloc_cond(tn, setk, SIG5[-k])
  inset <- sig[k] %in% set_rows[[setk]]
  id <- sprintf("S5.persig.%s.s%d_%s.%s.PPH4", tn, k, SIG5[k], setk)
  addr(id, "Rule 3(iv) per-signal table", sprintf("coloc.abf PP.H4, DPEP1 protein vs %s, both conditioned on the other four signals, over %s; p12 = 1e-5", TR[[tn]]$label, set_lab[[setk]]),
       sprintf("UKB-PPP OID20897 x %s", TR[[tn]]$label), sub("full", "", setk), value = num(cc$H4), n = cc$nsnp, unit = "PP.H4",
       notes = sprintf("PP.H3 %.4g; lead %s; signal SNP itself in the set: %s; %d collinear dropped%s", cc$H3, cc$lead, inset, cc$n_collinear,
                       if (tn == "eGFR") sprintf("; eGFR sdY estimated by coloc from N %.0f (window median n_total_sum) and reference MAF", P$N_EGFR_MED) else ""))
  persig[[paste("coloc", tn, k, setk)]] <- data.table(trait = tn, signal = k, rsid = SIG5[k], panel = setk, PPH4 = cc$H4, PPH3 = cc$H3, lead = cc$lead, nsnp = cc$nsnp, signal_in_set = inset)
  say("PER-SIGNAL coloc %s %s s%d %s: PP.H4 %.4f PP.H3 %.4f lead %s (n %d, signal in set %s)", setk, tn, k, SIG5[k], cc$H4, cc$H3, cc$lead, cc$nsnp, inset)
}

# ================================================================ D. s (Rule 3(i)) for protein and cysteinylglycine, 4 sets x 2 panels
suppressMessages(library(susieR))
s_tab <- list()
valid2 <- P$v2rows
for (setn in c("susie", "p1fr", "p1np", "p2np")) for (pn in c("503", "489")) {
  rows <- P$SETS[[setn]]
  if (pn == "489") rows <- rows[rows %in% valid2]
  ids <- rs[rows]
  Rm <- if (pn == "503") R503[ids, ids] else R489[ids, ids]
  stopifnot(all(is.finite(Rm)))
  te <- system.time(attr(Rm, "eigen") <- eigen(Rm, symmetric = TRUE))[["elapsed"]]
  for (tn in c("P__DPEP1", "M__GCST90200364")) {
    z <- H[[paste0("b__", tn)]][rows] / H[[paste0("se__", tn)]][rows]
    nn <- if (tn == "P__DPEP1") P$N_PROT else P$N_CYS
    s1 <- estimate_s_rss(z, Rm, nn)
    s2 <- estimate_s_rss(z, Rm, P$NHAT[[tn]])
    pipeline_run <- (setn %in% c("susie", "p1fr", "p1np") && pn == "503") || (setn == "p2np" && pn == "489")
    tl <- if (tn == "P__DPEP1") "protein" else "cysgly"
    addr(sprintf("S5.s.%s.%s.%s", tl, setn, pn), "Rule 3(i); §5", sprintf("susieR::estimate_s_rss (null-mle) for %s z over SNP set %s on the %s panel (PLINK r), n = %d",
         tl, setn, pn, nn), if (tn == "P__DPEP1") "UKB-PPP OID20897" else "CLSA GCST90200364", pn, value = num(s1), n = length(rows), unit = "s",
         notes = sprintf("%s; n_hat sensitivity (n = %.0f): s = %.5f; eigen %.0fs; locally computed s is a proxy for the pipeline fit",
                         if (pipeline_run) "set x panel on which the l3-finemap pipeline ran a SuSiE fit" else "Stage 1 extension: no pipeline fit ran on this set x panel", P$NHAT[[tn]], s2, te))
    s_tab[[paste(tn, setn, pn)]] <- data.table(trait = tl, set = setn, panel = pn, nsnp = length(rows), s = s1, s_nhat = s2, pipeline_run = pipeline_run)
    say("s %s %s %s: %d SNPs, s = %.5f (n %d) | n_hat %.0f -> %.5f | pipeline fit set x panel: %s", tl, setn, pn, length(rows), s1, nn, P$NHAT[[tn]], s2, pipeline_run)
  }
}
s_tab <- rbindlist(s_tab)
smax <- s_tab[, .(s_max_all = max(s), s_max_pipeline = max(s[pipeline_run]), s_max_nhat = max(s_nhat)), by = trait]
for (i in seq_len(nrow(smax))) {
  addr(sprintf("S5.s.%s.max", smax$trait[i]), "Rule 3(i)", "max s over the four sets x two panels (Rule 3(i) s for the pair's trait)", "", "both", value = num(smax$s_max_all[i]), unit = "s",
       notes = sprintf("max over the 4 set x panel combinations the pipeline actually fitted = %.5f; with n_hat = %.5f", smax$s_max_pipeline[i], smax$s_max_nhat[i]))
}
admissible <- all(smax$s_max_all < 0.10)
addr("S5.s.pair.admissible", "Rule 3(i)", "SuSiE with reference LD admissible for DPEP1-cysteinylglycine (s < 0.10 for both traits, s = max over every set and panel)", "", "both",
     value = admissible, notes = sprintf("protein max s %.4f; cys-gly max s %.4f", smax[trait == "protein"]$s_max_all, smax[trait == "cysgly"]$s_max_all))

# ================================================================ E. attribution separation (PLINK r^2, 503; 489 beside)
sep <- list()
for (pn in c("503", "489")) {
  Rm <- if (pn == "503") R503 else P$R489all
  r2 <- Rm[SIG5, CHMP]^2
  for (a in SIG5) for (b in CHMP) addr(sprintf("S5.sep.r2.%s.%s.%s", a, b, pn), "§5 attribution check", sprintf("PLINK r^2 between DPEP1 signal %s and CHMP1A signal %s", a, b),
                                       paste0("1000G EUR ", pn), pn, value = num(r2[a, b]), unit = "r^2")
  sep[[pn]] <- r2
  say("SEPARATION %s: max r^2 %.5f at %s; r^2 matrix rows %s x cols %s: %s", pn, max(r2), paste(which(r2 == max(r2), arr.ind = TRUE)[1, ], collapse = ","),
      paste(SIG5, collapse = ","), paste(CHMP, collapse = ","), paste(sprintf("%.4f", r2), collapse = " "))
}
addr("S5.sep.maxr2", "§5 attribution check", "max r^2 (503, PLINK --r) between every DPEP1 signal and each CHMP1A signal", "1000G EUR 503", "503",
     value = num(max(sep[["503"]])), unit = "r^2", notes = sprintf("separated iff < 0.10: %s; 489 beside: %.5f", max(sep[["503"]]) < 0.10, max(sep[["489"]])))
addr("S5.sep.maxr2.489", "§5 attribution check (beside)", "as S5.sep.maxr2 on the 489 panel", "1000G EUR 489", "489", value = num(max(sep[["489"]])), unit = "r^2")

# ================================================================ F. Rule 3(iii) gene naming on the restricted grid
GR <- list(ens116pc = inp("GRID_bridge_class_robustness__ens116pc"), glist = inp("GRID_bridge_class_robustness"), ens116all = inp("GRID_bridge_class_robustness__ens116all"))
gn <- list()
for (an in names(GR)) {
  g <- fread(GR[[an]]); stopifnot(nrow(g) == 1008L)
  g <- g[locus == "DPEP1_16q24" & protein == "DPEP1" & metabolite == "cysteinylglycine"]
  stopifnot(nrow(g) == 48L, uniqueN(g$variant) == 48L)
  gate_closed <- g$susie == "sguard" & (suppressWarnings(as.numeric(g$s_protein)) >= 0.10 | suppressWarnings(as.numeric(g$s_metabolite)) >= 0.10)
  comp <- g$susie == "none" | (g$susie %in% c("asnow", "sguard") & admissible) | gate_closed
  assigned <- sum(comp & g$class == "gene-attributed bridge"); tot <- sum(comp)
  strict_tot <- sum(g$susie == "none"); strict_ass <- sum(g$susie == "none" & g$class == "gene-attributed bridge")
  gn[[an]] <- list(assigned = assigned, total = tot, majority = tot > 0 && assigned > tot / 2, all48 = sum(g$class == "gene-attributed bridge"),
                   strict = c(strict_ass, strict_tot), gate_closed = sum(gate_closed))
  est <- sprintf("settings of the 48-cell SM attribution grid (%s annotation) that meet Rule 3(i) [susie=none; or sguard with its own s gate closed (grid s_protein or s_metabolite >= 0.10); asnow/sguard only if the pair's Stage 1 s < 0.10 for both traits (here %s)] assigning DPEP1 ('gene-attributed bridge')", an, admissible)
  addr(sprintf("S5.genename.%s.assigned", an), "Rule 3(iii); §5", est, basename(GR[[an]]), "503", value = assigned, n = tot, unit = "settings",
       notes = sprintf("majority (> half) = %s; unrestricted 48-cell count %d/48; strict restriction (susie=none only) %d/%d; sguard cells with gate closed %d", gn[[an]]$majority, gn[[an]]$all48, strict_ass, strict_tot, sum(gate_closed)))
  addr(sprintf("S5.genename.%s.total", an), "Rule 3(iii); §5", est, basename(GR[[an]]), "503", value = tot, unit = "settings")
  say("GENE NAMING %s: restricted grid %d settings, DPEP1 assigned in %d (majority %s); all 48: %d; strict none-only %d/%d; sguard gates closed %d",
      an, tot, assigned, gn[[an]]$majority, gn[[an]]$all48, strict_ass, strict_tot, sum(gate_closed))
}

saveRDS(list(gls = gls, two = two, clause = clause, persig = rbindlist(persig, fill = TRUE), s_tab = s_tab, smax = smax, admissible = admissible,
             sep = sep, gn = gn), file.path(WORK, "analyses.rds"))
fwrite(rbindlist(RES), file.path(WORK, "res_02.tsv"), sep = "\t")
say("==== S5_02_analyses done in %.0f s (%d result rows) ====", as.numeric(difftime(Sys.time(), t0, units = "secs")), length(RES))
