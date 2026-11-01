# 08_umod.R -- section 6 UMOD control. Inputs: LOC_UMOD_eGFR_control_harmonised (raw_input), LOC_UMOD_eGFR_control_signals
# (definition), REF_UMOD_REF503_I8 (read in place, bootstrap only), the 503/489 panels via PLINK, and this folder's p2np rebuild.
#   primary   all-signal GLS Q of eGFR-per-UMOD ratios over the 5 protein signals, 503 (gls_q.R; aware 2nd order); 489 beside
#   Rule 3(ii) 5,000 503-panel bootstrap draws, seed 20261007 (DEVIATION: bootstrap outside PLINK (sec 10))
#   secondary two-signal A (rs12922822) vs B (rs77924615) test, both panels; Rule 5 ids p_aware_2nd / p_ign / p_aware_1st / p_cojo
#   Rule 3(iv) cond PP.H4 at A and B: fresh R coloc.abf on z-space-conditioned data (engine two-signal module), p12 = 1e-5
#   grid      5 protein signals x 3 eGFR signals, each conditioned on its other signals (full harmonised set, 503)
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
source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "05_susie.R"))   # plink_ld + common
suppressMessages(library(coloc))
say("==== W6b 08_umod ====")
UD <- file.path(WORK, "umod"); dir.create(UD, showWarnings = FALSE)
source(res("DEF_GLS_Q_R"))                                    # shared estimator (md5 asserted by res())
stopifnot(md5c(GOV$gls_q[1]) == GOV$gls_q[2])

H <- fread(res("LOC_UMOD_eGFR_control_harmonised"))
sig <- fread(res("LOC_UMOD_eGFR_control_signals"))
stopifnot(nrow(H) == 2823, nrow(sig) == 8, !anyDuplicated(H$snp))
nn <- sapply(H[, .(b__P__UMOD, se__P__UMOD, b__O__eGFR, se__O__eGFR, af_ref)], function(v) sum(!is.finite(v)))
say("UMOD harmonised: %d SNPs; non-finite counts %s", nrow(H), paste(names(nn), nn, sep = "=", collapse = " "))
stopifnot(all(nn == 0))
H[, `:=`(bP = b__P__UMOD, sP = se__P__UMOD, bY = b__O__eGFR, sY = se__O__eGFR)]
psig <- sig[tid == "P__UMOD"][order(signal)]; ysig <- sig[tid == "O__eGFR"][order(signal)]
stopifnot(nrow(psig) == 5, nrow(ysig) == 3, all(psig$snp %in% H$snp), all(ysig$snp %in% H$snp))
A <- "rs12922822"; B <- "rs77924615"
ok_decl <- psig$rsid[1] == A && psig$rsid[2] == B && ysig$rsid[1] == B
log_assert("S6.umod.declared_signals", "A = rs12922822 (protein signal 1), B = rs77924615 (protein signal 2 = eGFR signal 1)",
           sprintf("P1 %s, P2 %s, eGFR1 %s", psig$rsid[1], psig$rsid[2], ysig$rsid[1]), ok_decl, "LOC_UMOD_eGFR_control_signals")
iA <- match(psig$snp[1], H$snp); iB <- match(psig$snp[2], H$snp); i5 <- match(psig$snp, H$snp)

# ---- LD: PLINK (503 full harmonised set; 489 p2np set from this folder's rebuild; 489 for the 5 signal SNPs) ----
ld503f <- file.path(UD, "R503_full.rds")
if (!file.exists(ld503f)) { L <- plink_ld(sub("\\.bed$", "", res("LD503_bed")), H$rsid_ref, file.path(UD, "ld503_full")); saveRDS(L, ld503f) }
L503 <- readRDS(ld503f); R503 <- L503$R; dimnames(R503) <- list(H$snp, H$snp)
say("UMOD LD 503: %dx%d (md5 %s)", nrow(R503), ncol(R503), L503$md5)
mine <- fread(file.path(WORK, "harm", "UMOD_eGFR.tsv.gz"))
stopifnot(all(mine$snp %in% H$snp), nrow(mine) == nrow(H))
p2 <- mine[in_p2np == TRUE]
pre2 <- file.path(ROOT, "lhcmr_eas_refs/1000G_EUR_Phase3_plink/1000G.EUR.QC.16"); res("LD489_chr16_bed"); res("LD489_chr16_bim")
ld489f <- file.path(UD, "R489_p2np.rds")
if (!file.exists(ld489f)) { L <- plink_ld(pre2, p2$p2_rsid, file.path(UD, "ld489_p2np")); L$R <- L$R * outer(p2$p2_sign, p2$p2_sign); saveRDS(L, ld489f) }
L489 <- readRDS(ld489f); R489 <- L489$R; dimnames(R489) <- list(p2$snp, p2$snp)
say("UMOD LD 489 (p2np): %dx%d (md5 %s)", nrow(R489), ncol(R489), L489$md5)
# the 5 signal SNPs on panel 2 (signal 5 rs62033284 is palindromic, so it is not in p2np; matched here by position + exact allele pair,
# orientation checked by frequency because the 489 are a subset of the 503)
b2 <- fread(file.path(WORK, "ref", "ref489_chr16.bim"), header = FALSE, col.names = c("chr", "rsid", "cm", "pos37", "a1", "a2"))
f2 <- fread(file.path(WORK, "ref", "ref489_chr16.frq.counts"))
stopifnot(all(f2$SNP == b2$rsid))
b2[, af2 := f2$C1 / (f2$C1 + f2$C2)]
s5 <- H[i5, .(snp, rsid_ref, pos37, a1, a2, af_ref)]
s5 <- merge(s5, b2[, .(pos37, p2_rsid = rsid, p2_a1 = a1, p2_a2 = a2, af2)], by = "pos37", all.x = TRUE, sort = FALSE)
s5 <- s5[(p2_a1 == a1 & p2_a2 == a2) | (p2_a1 == a2 & p2_a2 == a1)]
s5 <- s5[match(H$snp[i5], snp)]
stopifnot(nrow(s5) == 5, !anyNA(s5$p2_rsid))
s5[, p2_sign := ifelse(p2_a1 == a1, 1L, -1L)]
s5[, af2a := ifelse(p2_sign == 1L, af2, 1 - af2)]
stopifnot(max(abs(s5$af2a - s5$af_ref)) <= P2_AF_DIFF_MAX)
L5 <- plink_ld(pre2, s5$p2_rsid, file.path(UD, "ld489_sig5"))
R5_489 <- L5$R * outer(s5$p2_sign, s5$p2_sign); dimnames(R5_489) <- list(s5$snp, s5$snp)
R5_503 <- R503[i5, i5]

# ---- LD assertion mutation test (rule 3): one SNP id mutated must trip the 'every requested SNP returned' assertion ----
mut <- tryCatch({ plink_ld(pre2, c(s5$p2_rsid[1:4], paste0(s5$p2_rsid[5], "_MUT")), file.path(UD, "ld_mutation")); "NOT TRIPPED" },
                error = function(e) paste("TRIPPED:", conditionMessage(e)))
log_assert("MUT.ld_every_pair_returned", "TRIPPED", mut, startsWith(mut, "TRIPPED"), "plink_ld with rs62033284 id mutated to *_MUT")

# ---- cross-representation check: PLINK r vs r from the i8 dosage read in place (503 x 2823, SNP order = harmonised rows) ----
gfile <- res("REF_UMOD_REF503_I8")
gv <- readBin(gfile, what = "integer", n = file.size(gfile) + 1L, size = 1L, signed = TRUE)
stopifnot(length(gv) == 503L * nrow(H))
G <- matrix(gv, nrow = 503L, byrow = TRUE); rm(gv)
stopifnot(all(G %in% 0:2), max(abs(colMeans(G) / 2 - H$af_ref)) < 1e-5)
G5 <- G[, i5]
dR <- max(abs(cor(G5) - R5_503))
log_assert("S6.umod.ld_plink_vs_dosage", "max |r_PLINK - r_dosage| < 1e-6 over the 5 signal SNPs", sprintf("%.3g", dR), dR < 1e-6, "PLINK --r square vs cor() of REF_UMOD_REF503_I8 columns")
say("r(A,B): 503 %.5f | 489 %.5f", R503[iA, iB], R5_489[1, 2])

rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- rrow(...)
RS <- "UKB-PPP OID20237 x CKDGen 2019 EA eGFR (LOC_UMOD_eGFR_control_harmonised)"

# ---- primary: all-signal GLS Q ----
gls_cell <- function(R) {
  jP <- joint_effects(H$bP[i5], H$sP[i5], R); jY <- joint_effects(H$bY[i5], H$sY[i5], R)
  out <- list()
  for (cv in c("aware", "ign")) for (so in c("2nd", "1st")) {
    o <- profile_Q(jY$b, jY$V, jP$b, jP$V, cov = cv, se_order = so)
    out[[paste(cv, so, sep = "_")]] <- o
  }
  out
}
G503 <- gls_cell(R5_503); G489 <- gls_cell(R5_489)
for (pan in c("503", "489")) {
  GG <- if (pan == "503") G503 else G489
  pre <- if (pan == "503") "S6.umod.Q" else "S6.umod.Q.p489"
  for (k in c("aware_2nd", "ign_2nd", "aware_1st")) {
    o <- GG[[k]]; cv <- sub("_.*", "", k); so <- sub(".*_", "", k)
    rid <- c(aware_2nd = "p_aware_2nd (primary)", ign_2nd = "p_ign_2nd", aware_1st = "p_aware_1st")[[k]]
    est <- sprintf("GLS Q of eGFR-per-SD-UMOD ratios across the 5 conditionally distinct UMOD protein signals (joint effects from marginal z with %s LD); Rule 5 id %s", pan, rid)
    add(paste0(pre, ".", k, ".Q"), "6 UMOD primary; Rule 5", est, RS, pan, cv, so, value = o$Q, unit = "chi-square", notes = if (o$at_bound) "theta at grid bound" else "")
    add(paste0(pre, ".", k, ".df"), "6 UMOD primary; Rule 5", est, RS, pan, cv, so, value = o$df, unit = "df")
    add(paste0(pre, ".", k, ".p"), "6 UMOD primary; Rule 5", est, RS, pan, cv, so, value = o$p_upper, p = o$p_upper, unit = "p (upper tail)")
    add(paste0(pre, ".", k, ".p_lower"), "6 UMOD primary (lower tail, printed)", est, RS, pan, cv, so, value = o$p_lower, unit = "p (lower tail)")
    add(paste0(pre, ".", k, ".theta"), "6 UMOD primary", paste(est, "- profiled common ratio"), RS, pan, cv, so, value = o$theta, unit = "ln(eGFR) per SD UMOD")
  }
  add(paste0(pre, ".p_cojo"), "Rule 5", "beta-space (COJO) p for a GLS Q", RS, pan, value = NA, notes = "NOT RUN for GLS Q")
}
say("GLS 503: Q %.3f df %d p_aware_2nd %.4g | ign %.4g | 1st %.4g ; 489: p_aware_2nd %.4g", G503$aware_2nd$Q, G503$aware_2nd$df, G503$aware_2nd$p_upper,
    G503$ign_2nd$p_upper, G503$aware_1st$p_upper, G489$aware_2nd$p_upper)
# check of the shared estimator's optimum: gls_q.R profile_Q refines with optimize() at its default tolerance (~1.2e-4 in theta),
# which is coarse when theta ~ 0.01 (eGFR in ln units per SD). The declared value is the gls_q.R value; the exact optimum
# (optimize tol 1e-12 on +/-0.02 around it, cross-checked by a 40,001-point grid) is printed beside it, for every cell.
exact_opt <- function(R, cv, so) {
  jP <- joint_effects(H$bP[i5], H$sP[i5], R); jY <- joint_effects(H$bY[i5], H$sY[i5], R)
  VY <- jY$V; VP <- jP$V
  if (cv == "ign") { VY <- diag(diag(VY)); VP <- diag(diag(VP)) }
  if (so == "1st") VP <- VP * 0
  Qf <- function(t) { e <- jY$b - t * jP$b; V <- VY + t^2 * VP; as.numeric(t(e) %*% solve(V, e)) }
  o <- profile_Q(jY$b, jY$V, jP$b, jP$V, cov = cv, se_order = so)
  ox <- optimize(Qf, c(o$theta - 0.02, o$theta + 0.02), tol = 1e-12)
  th <- o$theta + seq(-0.02, 0.02, length.out = 40001); gmin <- min(vapply(th, Qf, 1))
  list(Q = o$Q, Qx = min(ox$objective, gmin), p = o$p_upper, px = pchisq(min(ox$objective, gmin), o$df, lower.tail = FALSE), theta_x = ox$minimum)
}
EXO <- list()
for (pan in c("503", "489")) for (k in c("aware_2nd", "ign_2nd", "aware_1st")) {
  e <- exact_opt(if (pan == "503") R5_503 else R5_489, sub("_.*", "", k), sub(".*_", "", k)); EXO[[paste(pan, k)]] <- e
  say("GLS optimum %s %s: gls_q.R Q %.6f p %.6f | exact Q %.6f p %.6f (dQ %.2g)", pan, k, e$Q, e$p, e$Qx, e$px, e$Q - e$Qx)
}
for (pan in c("503", "489")) for (k in c("aware_2nd", "ign_2nd", "aware_1st")) {
  e <- EXO[[paste(pan, k)]]; pre <- if (pan == "503") "S6.umod.Q" else "S6.umod.Q.p489"
  add(paste0(pre, ".", k, ".p_exact_optimum"), "6 UMOD (estimator check, printed beside)", "same GLS cell, Q minimised exactly (optimize tol 1e-12 + fine grid) instead of gls_q.R's default-tolerance optimize",
      RS, pan, sub("_.*", "", k), sub(".*_", "", k), value = e$px, p = e$px, unit = "p (upper tail)",
      notes = sprintf("Q exact %.6f vs gls_q.R %.6f (dQ %.2g); theta exact %.8g; the declared value is the gls_q.R row", e$Qx, e$Q, e$Q - e$Qx, e$theta_x))
}
e0 <- EXO[["503 aware_2nd"]]
log_assert("S6.umod.Q.profile_min_strict", "gls_q.R Q within 1e-8 of the exact minimum (informative)", sprintf("dQ %.3g", e0$Q - e0$Qx), (e0$Q - e0$Qx) < 1e-8,
           "optimize tol 1e-12 + 40,001-point grid on +/-0.02 around the gls_q.R theta")
matl <- all(vapply(EXO, function(e) abs(e$p - e$px) < 1e-3 && ((e$p < 0.05) == (e$px < 0.05)), TRUE))
log_assert("S6.umod.Q.profile_min_material", "every cell: |p_gls_q - p_exact| < 1e-3 and the same side of 0.05", matl, matl, "6 cells (2 panels x 3 Rule 5 variants)")

# ---- Rule 3(ii): 5,000 503-panel bootstrap draws of the primary cell ----
# Draw d resamples people with the d-th sample.int() call of one seed-20261007 stream (identical to a sequential loop), so the
# draws can be shared among workers without changing any result. Checkpoints: boot_gls.tsv (first sequential run) + boot/w*.tsv.
BOOT <- file.path(UD, "boot_gls.tsv"); BOOTD <- file.path(UD, "boot"); dir.create(BOOTD, showWarnings = FALSE)
set.seed(20261007L)
IDX <- t(vapply(1:5000, function(d) sample.int(503L, 503L, replace = TRUE), integer(503)))
read_done <- function() {
  f <- c(if (file.exists(BOOT)) BOOT, list.files(BOOTD, pattern = "^w[0-9]+\\.tsv$", full.names = TRUE))
  if (!length(f)) return(data.table(draw = integer(), p = numeric(), Q = numeric(), theta = numeric(), n_mono = integer()))
  unique(rbindlist(lapply(f, fread)), by = "draw")
}
boot_one <- function(d) {
  Gd <- G5[IDX[d, ], , drop = FALSE]
  sdv <- apply(Gd, 2, sd)
  Rd <- suppressWarnings(cor(Gd)); mono <- which(sdv == 0)
  if (length(mono)) { Rd[mono, ] <- 0; Rd[, mono] <- 0; diag(Rd) <- 1 }
  jP <- joint_effects(H$bP[i5], H$sP[i5], Rd); jY <- joint_effects(H$bY[i5], H$sY[i5], Rd)
  o <- profile_Q(jY$b, jY$V, jP$b, jP$V, cov = "aware", se_order = "2nd")
  data.table(draw = d, p = o$p_upper, Q = o$Q, theta = o$theta, n_mono = length(mono))
}
done <- read_done()
todo <- setdiff(1:5000, done$draw)
say("UMOD GLS bootstrap: %d draws already checkpointed, %d to do", nrow(done), length(todo))
if (length(todo)) {
  nw <- 4L; t0 <- Sys.time()
  chunks <- split(todo, rep_len(seq_len(nw), length(todo)))
  library(parallel); cl <- makeCluster(nw)
  clusterExport(cl, c("IDX", "G5", "H", "i5", "BOOTD", "joint_effects", "profile_Q", "boot_one"), envir = environment())
  invisible(clusterEvalQ(cl, suppressMessages(library(data.table))))
  parLapply(cl, seq_len(nw), function(k, chunks) {
    f <- file.path(BOOTD, sprintf("w%d.tsv", k)); buf <- list(); dd <- chunks[[k]]
    for (j in seq_along(dd)) {
      buf[[length(buf) + 1]] <- boot_one(dd[j])
      if (length(buf) == 100L || j == length(dd)) { fwrite(rbindlist(buf), f, sep = "\t", append = file.exists(f)); buf <- list() }
    }
    NULL
  }, chunks = chunks)
  stopCluster(cl)
  say("UMOD GLS bootstrap: %d draws computed by %d workers in %.0fs", length(todo), nw, as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
bt <- read_done(); setorder(bt, draw)
# identity of the parallel scheme: re-compute 3 already-checkpointed draws in this process and compare
chk_d <- bt$draw[c(1, 2500, 5000)]
re <- rbindlist(lapply(chk_d, boot_one))
log_assert("S6.umod.boot.recompute", "3 checkpointed draws re-computed identically", sprintf("max|dp| %.3g", max(abs(re$p - bt[match(chk_d, draw)]$p))),
           max(abs(re$p - bt[match(chk_d, draw)]$p)) < 1e-12, "draws 1st, 2500th, 5000th of the checkpoint vs a fresh computation")
stopifnot(nrow(bt) == 5000, !anyDuplicated(bt$draw), !anyNA(bt$p))
share <- mean((bt$p < 0.05) == (G503$aware_2nd$p_upper < 0.05))
pq <- quantile(bt$p, c(0.025, 0.975), names = FALSE)
estb <- "Rule 3(ii) bootstrap of the primary cell (all-signal GLS, aware, 2nd order, 503): 5,000 draws of the 503 people with replacement, seed 20261007; share = draws on the full-panel side of p 0.05"
add("S6.umod.Q.boot.share", "Rule 3(ii); 6 UMOD", estb, RS, "503 bootstrap", "aware", "2nd", value = share, n = 5000, unit = "share",
    notes = sprintf("DEVIATION: bootstrap outside PLINK (sec 10): r from cor() of REF_UMOD_REF503_I8 read in place; draws with a monomorphic signal SNP: %d (r set to 0)", sum(bt$n_mono > 0)))
add("S6.umod.Q.boot.p_lo", "Rule 3(ii); 6 UMOD", paste(estb, "- 2.5th percentile of draw p"), RS, "503 bootstrap", "aware", "2nd", value = pq[1], n = 5000, unit = "p")
add("S6.umod.Q.boot.p_hi", "Rule 3(ii); 6 UMOD", paste(estb, "- 97.5th percentile of draw p"), RS, "503 bootstrap", "aware", "2nd", value = pq[2], n = 5000, unit = "p")
say("UMOD GLS bootstrap: share %.4f, p 2.5-97.5%% %.4g-%.4g", share, pq[1], pq[2])

# ---- two-signal A vs B (secondary) ----
two_signal <- function(r, fA, fB, mode = c("conditional", "marginal", "cojo")) {
  mode <- match.arg(mode)
  b <- list(P = H$bP[c(iA, iB)], Y = H$bY[c(iA, iB)]); s <- list(P = H$sP[c(iA, iB)], Y = H$sY[c(iA, iB)])
  R2 <- matrix(c(1, r, r, 1), 2)
  if (mode == "conditional") { jP <- joint_effects(b$P, s$P, R2); jY <- joint_effects(b$Y, s$Y, R2) }
  if (mode == "marginal") { jP <- list(b = b$P, V = diag(s$P) %*% R2 %*% diag(s$P)); jY <- list(b = b$Y, V = diag(s$Y) %*% R2 %*% diag(s$Y)) }
  if (mode == "cojo") {
    v <- 2 * c(fA, fB) * (1 - c(fA, fB)); Sx <- matrix(c(v[1], r * sqrt(v[1] * v[2]), r * sqrt(v[1] * v[2]), v[2]), 2); Si <- solve(Sx)
    cj <- function(bb, ss) list(b = as.vector(Si %*% (v * bb)), V = mean(v * ss^2) * Si)
    jP <- cj(b$P, s$P); jY <- cj(b$Y, s$Y)
  }
  k <- jY$b / jP$b
  v2 <- (jY$V[cbind(1:2, 1:2)] / jP$b^2) + (jY$b^2 * jP$V[cbind(1:2, 1:2)] / jP$b^4)
  c2 <- jY$V[1, 2] / (jP$b[1] * jP$b[2]) + jY$b[1] * jY$b[2] / (jP$b[1]^2 * jP$b[2]^2) * jP$V[1, 2]
  v1 <- jY$V[cbind(1:2, 1:2)] / jP$b^2; c1 <- jY$V[1, 2] / (jP$b[1] * jP$b[2])
  d <- k[1] - k[2]
  se_a2 <- sqrt(v2[1] + v2[2] - 2 * c2); se_i2 <- sqrt(v2[1] + v2[2]); se_a1 <- sqrt(v1[1] + v1[2] - 2 * c1)
  list(kA = k[1], kB = k[2], d = d, se_aware_2nd = se_a2, se_ign = se_i2, se_aware_1st = se_a1,
       p_aware_2nd = 2 * pnorm(-abs(d / se_a2)), p_ign = 2 * pnorm(-abs(d / se_i2)), p_aware_1st = 2 * pnorm(-abs(d / se_a1)),
       zP = jP$b / sqrt(jP$V[cbind(1:2, 1:2)]), zY = jY$b / sqrt(jY$V[cbind(1:2, 1:2)]))
}
afA2 <- s5$af2a[1]; afB2 <- s5$af2a[2]
TS <- list(`503` = list(r = R503[iA, iB], fA = H$af_ref[iA], fB = H$af_ref[iB]), `489` = list(r = R5_489[1, 2], fA = afA2, fB = afB2))
for (pan in names(TS)) {
  x <- TS[[pan]]
  cnd <- two_signal(x$r, x$fA, x$fB, "conditional"); mar <- two_signal(x$r, x$fA, x$fB, "marginal"); cj <- two_signal(x$r, x$fA, x$fB, "cojo")
  est <- sprintf("two-signal test, eGFR-per-SD-UMOD ratio at A=rs12922822 (conditioned on B) vs at B=rs77924615 (conditioned on A), z-space conditioning with %s-panel r(A,B) = %.5f", pan, x$r)
  pp <- sprintf("S6.umod.two.%s", pan)
  add(paste0(pp, ".kA"), "6 UMOD secondary", est, RS, pan, value = cnd$kA, unit = "ln(eGFR) per SD UMOD")
  add(paste0(pp, ".kB"), "6 UMOD secondary", est, RS, pan, value = cnd$kB, unit = "ln(eGFR) per SD UMOD")
  add(paste0(pp, ".diff"), "6 UMOD secondary", paste(est, "- kA - kB"), RS, pan, "aware", "2nd", value = cnd$d, se = cnd$se_aware_2nd, ci_lo = cnd$d - 1.959964 * cnd$se_aware_2nd, ci_hi = cnd$d + 1.959964 * cnd$se_aware_2nd, unit = "ln(eGFR) per SD UMOD")
  add(paste0(pp, ".p_aware_2nd"), "6 UMOD secondary; Rule 5 primary", est, RS, pan, "aware", "2nd", value = cnd$p_aware_2nd, se = cnd$se_aware_2nd, p = cnd$p_aware_2nd, unit = "p")
  add(paste0(pp, ".p_ign"), "6 UMOD secondary; Rule 5", est, RS, pan, "ign", "2nd", value = cnd$p_ign, se = cnd$se_ign, p = cnd$p_ign, unit = "p")
  add(paste0(pp, ".p_aware_1st"), "6 UMOD secondary; Rule 5", est, RS, pan, "aware", "1st", value = cnd$p_aware_1st, se = cnd$se_aware_1st, p = cnd$p_aware_1st, unit = "p")
  add(paste0(pp, ".p_cojo"), "6 UMOD secondary; Rule 5", sub("z-space conditioning", "beta-space (COJO) joint estimates from reference A1 frequencies, equal-N, residual variance = mean(2f(1-f)se^2) over A and B,", est),
      RS, pan, "aware", "2nd", value = cj$p_aware_2nd, se = cj$se_aware_2nd, p = cj$p_aware_2nd, unit = "p", notes = sprintf("kA %.6g kB %.6g", cj$kA, cj$kB))
  add(paste0(pp, ".marginal.p_aware_2nd"), "6 UMOD secondary (marginal estimand, printed beside)", sub("\\(conditioned on [AB]\\)", "(marginal)", est), RS, pan, "aware", "2nd",
      value = mar$p_aware_2nd, se = mar$se_aware_2nd, p = mar$p_aware_2nd, unit = "p", notes = sprintf("marginal kA %.6g kB %.6g; p_ign %.4g; p_aware_1st %.4g", mar$kA, mar$kB, mar$p_ign, mar$p_aware_1st))
  say("two-signal %s: kA %.5g kB %.5g | p_aware_2nd %.4g p_ign %.4g p_aware_1st %.4g p_cojo %.4g | marginal p_aware_2nd %.4g", pan, cnd$kA, cnd$kB, cnd$p_aware_2nd,
      cnd$p_ign, cnd$p_aware_1st, cj$p_aware_2nd, mar$p_aware_2nd)
  TS[[pan]]$res <- cnd
}

# ---- conditioning (engine l3fm__conditional.condition: single or multi-SNP S; SNPs with r_iS R_SS^-1 r_Si > 0.9 dropped) ----
condition <- function(z, se, R, S) {
  RSS <- R[S, S, drop = FALSE]; RiS <- R[, S, drop = FALSE]
  W <- solve(RSS, t(RiS)); qv <- rowSums(RiS * t(W))
  zc <- (z - as.vector(RiS %*% solve(RSS, z[S]))) / sqrt(pmax(1 - qv, 1e-12)); sc <- se / sqrt(pmax(1 - qv, 1e-12))
  bad <- qv > 0.9; bad[S] <- TRUE
  zc[bad] <- NA; sc[bad] <- NA
  list(z = zc, se = sc)
}
abf <- function(cp, cy, snps, N = NULL, maf = NULL) {
  ok <- !is.na(cp$z) & !is.na(cy$z)
  if (sum(ok) < 20) return(NULL)
  d1 <- list(beta = (cp$z * cp$se)[ok], varbeta = cp$se[ok]^2, snp = snps[ok], type = "quant")
  d2 <- list(beta = (cy$z * cy$se)[ok], varbeta = cy$se[ok]^2, snp = snps[ok], type = "quant")
  if (is.null(N)) { d1$sdY <- 1; d2$sdY <- 1 } else { d1$N <- N[1]; d2$N <- N[2]; d1$MAF <- maf[ok]; d2$MAF <- maf[ok] }
  r <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  c(H4 = unname(r$summary["PP.H4.abf"]), H3 = unname(r$summary["PP.H3.abf"]), nsnp = sum(ok))
}
st <- fromJSON(file.path(WORK, "harm", "UMOD_eGFR__stamp.json"))
for (pan in c("503", "489")) {
  if (pan == "503") { hs <- H[palindromic == 0]; R <- R503[hs$snp, hs$snp] } else { hs <- H[match(p2$snp, snp)]; R <- R489[hs$snp, hs$snp] }
  ja <- match(psig$snp[1], hs$snp); jb <- match(psig$snp[2], hs$snp)
  stopifnot(!is.na(ja), !is.na(jb))
  zP <- hs$bP / hs$sP; zY <- hs$bY / hs$sY
  pB <- condition(zP, hs$sP, R, jb); yB <- condition(zY, hs$sY, R, jb)          # only A left
  pA <- condition(zP, hs$sP, R, ja); yA <- condition(zY, hs$sY, R, ja)          # only B left
  hA <- abf(pB, yB, hs$snp); hB <- abf(pA, yA, hs$snp)
  maf <- pmin(hs$af_ref, 1 - hs$af_ref)
  hA2 <- abf(pB, yB, hs$snp, c(st$n_P, st$n_O), maf); hB2 <- abf(pA, yA, hs$snp, c(st$n_P, st$n_O), maf)
  set_lab <- if (pan == "503") "p1np (503)" else "p2np (489)"
  estc <- "conditional PP.H4 of eGFR with UMOD at signal %s (both traits conditioned on signal %s; coloc.abf p1 = p2 = 1e-4, p12 = 1e-5, prior sd 0.15 = engine convention sdY = 1)"
  add(sprintf("S6.umod.H4.A.%s", pan), "6 UMOD colocalization; Rule 3(iv)", sprintf(estc, "A", "B"), RS, pan, value = hA[["H4"]], n = hA[["nsnp"]], unit = "PP.H4",
      notes = sprintf("set %s; PP.H3 %.4g; sensitivity with sdY estimated from N and MAF: PP.H4 %.4g", set_lab, hA[["H3"]], hA2[["H4"]]))
  add(sprintf("S6.umod.H4.B.%s", pan), "6 UMOD colocalization; Rule 3(iv)", sprintf(estc, "B", "A"), RS, pan, value = hB[["H4"]], n = hB[["nsnp"]], unit = "PP.H4",
      notes = sprintf("set %s; PP.H3 %.4g; sensitivity with sdY estimated from N and MAF: PP.H4 %.4g", set_lab, hB[["H3"]], hB2[["H4"]]))
  cz <- TS[[pan]]$res
  add(sprintf("S6.umod.z.A.%s", pan), "Rule 3(iv)", "conditional eGFR z at A given B (two-SNP z-space conditioning)", RS, pan, value = cz$zY[1], unit = "z", notes = sprintf("UMOD protein z(A|B) %.3f", cz$zP[1]))
  add(sprintf("S6.umod.z.B.%s", pan), "Rule 3(iv)", "conditional eGFR z at B given A (two-SNP z-space conditioning)", RS, pan, value = cz$zY[2], unit = "z", notes = sprintf("UMOD protein z(B|A) %.3f", cz$zP[2]))
  say("cond coloc %s: H4_A %.4f (H3 %.4f) H4_B %.4f (H3 %.4f) | eGFR z(A|B) %.3f z(B|A) %.3f", pan, hA[["H4"]], hA[["H3"]], hB[["H4"]], hB[["H3"]], cz$zY[1], cz$zY[2])
  TS[[pan]]$H4 <- c(A = hA[["H4"]], B = hB[["H4"]])
}

# ---- 5 x 3 signal grid (each signal conditioned on its trait's OTHER signals; full harmonised set, 503) ----
zP <- H$bP / H$sP; zY <- H$bY / H$sY
iy <- match(ysig$snp, H$snp)
for (kp in 1:5) {
  cp <- condition(zP, H$sP, R503, setdiff(i5, i5[kp]))
  jt <- which.max(abs(cp$z))                    # diagnostic: the protein's top residual SNP after conditioning on its other signals
  diag_p <- sprintf("protein top residual after conditioning: %s (z %.2f; r with signal %d = %.3f); z at signal %d itself %.2f", H$rsid_ref[jt], cp$z[jt], kp,
                    R503[jt, i5[kp]], kp, cp$z[i5[kp]])
  for (km in 1:3) {
    cy <- condition(zY, H$sY, R503, setdiff(iy, iy[km]))
    h <- abf(cp, cy, H$snp)
    add(sprintf("S6.umod.grid.P%dxM%d", kp, km), "6 UMOD colocalization grid", sprintf("coloc.abf PP.H4, UMOD signal %d (%s) | other UMOD signals x eGFR signal %d (%s) | other eGFR signals; full harmonised set, 503 LD, p12 = 1e-5",
        kp, psig$rsid[kp], km, ysig$rsid[km]), RS, "503", value = h[["H4"]], n = h[["nsnp"]], unit = "PP.H4", notes = sprintf("PP.H3 %.4g; %s", h[["H3"]], diag_p))
  }
}

# ---- control criterion + Rule 3(ii) + admissibility inputs ----
crit_h4 <- all(c(TS$`503`$H4, TS$`489`$H4) >= 0.80)
crit_q <- G503$aware_2nd$p_upper >= 0.05
add("S6.umod.control.criterion", "6 UMOD 'behaved as a control'", "cond_H4_A and cond_H4_B >= 0.80 on both panels AND all-signal Q p (503, aware, 2nd order, unadjusted) >= 0.05",
    RS, "503+489", value = as.character(crit_h4 && crit_q), notes = sprintf("H4 criterion %s; Q criterion %s", crit_h4, crit_q))
robust <- ((G503$aware_2nd$p_upper < 0.05) == (G489$aware_2nd$p_upper < 0.05)) && share >= 0.90
add("S6.umod.Q.rule3ii.robust", "Rule 3(ii)", "conclusion (side of p 0.05) same on 503 and 489 AND in >= 90% of 5,000 503-panel draws", RS, "503+489", "aware", "2nd",
    value = as.character(robust), notes = sprintf("503 p %.4g, 489 p %.4g, share %.4f", G503$aware_2nd$p_upper, G489$aware_2nd$p_upper, share))
# ---- S2-D25 item 2: per-signal ratio SEs for Figure 3a (NEW output only; nothing above or in umod_summary.rds changes) ----
# Per-signal ratio k_i = bY_i / bP_i from the joint effects of the primary cell (jP, jY as in gls_cell). SE = square root of the diagonal of the
# delta-method covariance used by the two-signal test (two_signal(), v2): V_YY[i,i] / bP_i^2 + bY_i^2 V_PP[i,i] / bP_i^4 ('2nd'; the primary
# se_order); the protein-side term dropped is the '1st' column (two_signal() v1). The diagonal does not depend on aware/ign (only off-diagonals do).
signal_se <- rbindlist(lapply(c("503", "489"), function(pan) {
  Rm <- if (pan == "503") R5_503 else R5_489
  jP <- joint_effects(H$bP[i5], H$sP[i5], Rm); jY <- joint_effects(H$bY[i5], H$sY[i5], Rm)
  GG <- if (pan == "503") G503 else G489
  k <- jY$b / jP$b
  stopifnot(isTRUE(all.equal(k, GG$aware_2nd$ratios, tolerance = 1e-12)))     # the ratio vector is the one stored by profile_Q
  vY <- diag(jY$V); vP <- diag(jP$V)
  data.table(panel = pan, signal = 1:5, rsid = psig$rsid, snp = psig$snp, ratio = k, bY_joint = jY$b, bP_joint = jP$b,
             se_2nd = sqrt(vY / jP$b^2 + jY$b^2 * vP / jP$b^4), se_1st = sqrt(vY / jP$b^2))
}))
stopifnot(all(is.finite(signal_se$se_2nd)), all(signal_se$se_2nd > 0), all(signal_se$se_2nd >= signal_se$se_1st))
fwrite(signal_se, file.path(UD, "umod_signal_se.tsv"), sep = "\t")
say("per-signal ratio SEs written: work/umod/umod_signal_se.tsv (%d rows)", nrow(signal_se))
saveRDS(list(G503 = G503, G489 = G489, share = share, pq = pq, TS = TS, crit_h4 = crit_h4, crit_q = crit_q, robust = robust), file.path(UD, "umod_summary.rds"))
write_rows(rbindlist(rows), "08_umod")
say("08_umod done")
