# S5_03_bootstrap.R -- W5 (§5 DPEP1) step 3: Rule 3(ii) LD bootstrap, 5,000 draws of the 503 reference people (seed 20261007).
# DEVIATION: bootstrap outside PLINK (§10) -- PLINK cannot resample individuals; the full-panel r comes from PLINK (S5_01),
# each draw's r is the Pearson r of the 503 A1 dosages read IN PLACE from the i8 extract (identical to PLINK r on the full panel,
# max |diff| 4.4e-16, logged in S5_01).
# Per draw: (1) GLS Q, all five signals, aware, 2nd order (the declared primary cell) for cysteinylglycine and eGFR;
# (2) the two-signal re-test rs409170 vs rs62068372 (aware, 2nd order): p and sign(kA - kB);
# (3) the conditional PP.H4 at A|B and B|A on p1np (Rule 3(ii) 'side of PP.H4 0.80' for conditional colocalization, printed),
#     computed with abf_direct() (S5_common.R; equal to coloc::coloc.abf, asserted by S5_05 check C2 and below on draw 1).
# Modes:
#   Rscript S5_03_bootstrap.R               main: merge every checkpoint/part file (overlaps must agree), re-verify 24
#                                           checkpointed draws by recomputation, compute the rest, write work/boot_draws.tsv
#   Rscript S5_03_bootstrap.R --range a:b   worker: draws a..b not yet in any checkpoint -> work/boot_draws.part_a_b.tsv
# The draws are generated in full from the seed before any is used, so every process uses the same draw b.
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
args <- commandArgs(trailingOnly = TRUE)
RANGE <- if (length(args) >= 2 && args[1] == "--range") as.integer(strsplit(args[2], ":", fixed = TRUE)[[1]]) else NULL
TAG <- if (is.null(RANGE)) "main" else sprintf("worker %d:%d", RANGE[1], RANGE[2])
t0 <- Sys.time()
say("==== S5_03_bootstrap start (%s) ====", TAG)
assert_governing()
assert_md5(GOV$GLS_Q$p, GOV$GLS_Q$m, "gls_q.R (sourced)")
source(GOV$GLS_Q$p)
P <- readRDS(file.path(WORK, "prep.rds"))
H <- P$H; rs <- H$rsid_ref; ix <- function(x) match(x, rs)
gp <- inp("REF_DPEP1_REF503_I8")
gv <- readBin(gp, what = "integer", n = file.size(gp) + 1L, size = 1L, signed = TRUE)
stopifnot(length(gv) == 503L * 3575L)
G <- matrix(gv, nrow = 503L, byrow = TRUE); rm(gv)
storage.mode(G) <- "double"

sig <- ix(SIG5); ab <- ix(c(SIG_A, SIG_B)); setr <- P$SETS$p1np
stopifnot(all(ab %in% setr))
bP <- H$b__P__DPEP1; sP <- H$se__P__DPEP1; zP <- bP / sP
bC <- H$b__M__GCST90200364; sC <- H$se__M__GCST90200364; zC <- bC / sC
bE <- P$EGFR$b; sE <- P$EGFR$se

# draws: fixed by the seed, generated in full before any is used
set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
IDX <- t(vapply(seq_len(N_BOOT), function(b) sample.int(503L, 503L, replace = TRUE), integer(503)))
idx_md5 <- digest::digest(IDX, algo = "md5")
say("bootstrap (%s): %d draws x 503 people, seed %d (Mersenne-Twister / Inversion / Rejection); draw-matrix md5 %s; first draw starts %s",
    TAG, N_BOOT, SEED, idx_md5, paste(IDX[1, 1:6], collapse = ","))
ckm <- file.path(WORK, "boot_draws.idxmd5")
if (file.exists(ckm)) { if (readLines(ckm) != idx_md5) stop("checkpoints were made with different draw indices; delete work/boot_draws*.tsv") } else writeLines(idx_md5, ckm)

# keep only the columns a draw needs (signals, A/B, p1np), re-indexed; no change to any value
keepc <- sort(unique(c(sig, ab, setr)))
Gk <- G[, keepc, drop = FALSE]; rm(G)
kx <- function(j) match(j, keepc)
sig_k <- kx(sig); ab_k <- kx(ab); set_k <- kx(setr)
COLS <- c("draw", "valid_sig", "p_Q_cysgly", "p_Q_eGFR", "valid_ab", "r_AB", "p_two", "d_two", "H4_A", "H4_B", "n_mono_set")
one_draw <- function(b, use_coloc_pkg = FALSE) {
  g <- Gk[IDX[b, ], , drop = FALSE]
  gs <- g[, sig_k]
  sdv <- apply(gs, 2, sd)
  out <- list(draw = b, valid_sig = all(sdv > 0))
  if (out$valid_sig) {
    R5 <- cor(gs)
    jp <- joint_effects(bP[sig], sP[sig], R5)
    jc <- joint_effects(bC[sig], sC[sig], R5); je <- joint_effects(bE[sig], sE[sig], R5)
    out$p_Q_cysgly <- profile_Q(jc$b, jc$V, jp$b, jp$V, cov = "aware", se_order = "2nd")$p_upper
    out$p_Q_eGFR <- profile_Q(je$b, je$V, jp$b, jp$V, cov = "aware", se_order = "2nd")$p_upper
  } else { out$p_Q_cysgly <- NA_real_; out$p_Q_eGFR <- NA_real_ }
  gab <- g[, ab_k]
  out$valid_ab <- all(apply(gab, 2, sd) > 0)
  if (out$valid_ab) {
    r2 <- cor(gab); out$r_AB <- r2[1, 2]
    tw <- ratio_diff(joint_effects(bP[ab], sP[ab], r2), joint_effects(bC[ab], sC[ab], r2), "aware", "2nd")
    out$p_two <- tw$p; out$d_two <- tw$d
    # conditional colocalization on p1np: r of every set SNP with A and with B in this draw (monomorphic -> r = 0, counted)
    gsSet <- g[, set_k, drop = FALSE]
    rr <- suppressWarnings(cor(gsSet, gab)); nmono <- sum(!is.finite(rr[, 1])); rr[!is.finite(rr)] <- 0
    for (sx in c("A", "B")) {
      cS <- if (sx == "A") 2L else 1L                  # A|B conditions on B; B|A conditions on A
      cidx <- if (sx == "A") ab[2] else ab[1]
      RiS <- matrix(rr[, cS], ncol = 1); RSS <- matrix(1, 1, 1)
      isS <- setr == cidx
      cp <- condition_z(zP[setr], sP[setr], RiS, RSS, zP[cidx], isS)
      cm <- condition_z(zC[setr], sC[setr], RiS, RSS, zC[cidx], isS)
      out[[paste0("H4_", sx)]] <- if (use_coloc_pkg) coloc_pair(cp$z * cp$se, cp$se, cm$z * cm$se, cm$se, rs[setr])$H4 else
        abf_direct(cp$z * cp$se, cp$se, cm$z * cm$se, cm$se)$H4
    }
    out$n_mono_set <- nmono
  } else { out$r_AB <- NA_real_; out$p_two <- NA_real_; out$d_two <- NA_real_; out$H4_A <- NA_real_; out$H4_B <- NA_real_; out$n_mono_set <- NA_integer_ }
  as.data.table(out)[, ..COLS]
}

# ---------------------------------------------------------------- read every checkpoint (main + parts + the 2026-10-08 helper file)
ck_main <- file.path(WORK, "boot_draws.tsv")
PART_PAT <- "^boot_draws\\.(part_[0-9a-z_]+|worker)\\.tsv$"   # unmerged part files only (merged ones end .merged_into_main.tsv)
read_all <- function() {
  fs <- c(ck_main, list.files(WORK, pattern = PART_PAT, full.names = TRUE))
  fs <- fs[file.exists(fs) & file.size(fs) > 0]
  if (!length(fs)) return(data.table())
  d <- rbindlist(lapply(fs, function(f) fread(f)[, src := basename(f)]), fill = TRUE)
  d
}
D0 <- read_all()
if (nrow(D0)) {
  dupd <- D0[duplicated(draw) | duplicated(draw, fromLast = TRUE)]
  if (nrow(dupd)) {
    num_cols <- setdiff(COLS, "draw")
    dd <- dupd[, lapply(.SD, function(x) diff(range(as.numeric(x), na.rm = TRUE))), by = draw, .SDcols = num_cols]
    worst <- max(unlist(dd[, ..num_cols]), na.rm = TRUE)
    say("checkpoint overlap (%s): %d draws present in more than one file; max abs disagreement over all columns %.3g", TAG, uniqueN(dupd$draw), worst)
    if (worst > 1e-6) stop("checkpoint files disagree on a shared draw")
  }
  D0 <- unique(D0, by = "draw")
}
say("checkpoints (%s): %d distinct draws already computed (%s)", TAG, nrow(D0), if (nrow(D0)) paste(D0[, .N, by = src][, sprintf("%s %d", src, N)], collapse = "; ") else "none")

if (is.null(RANGE)) {
  # re-verify a spread of checkpointed draws by recomputation with the current code (incl. coloc.abf vs abf_direct on draw 1)
  if (nrow(D0)) {
    vd <- unique(c(1L, D0$draw[unique(round(seq(1, nrow(D0), length.out = 24)))]))
    vd <- vd[vd %in% D0$draw]
    re <- rbindlist(lapply(vd, one_draw))
    old <- D0[match(vd, draw)]
    num_cols <- setdiff(COLS, c("draw", "valid_sig", "valid_ab", "n_mono_set"))
    dmax <- max(abs(as.matrix(re[, ..num_cols]) - as.matrix(old[, ..num_cols])), na.rm = TRUE)
    pk <- one_draw(1L, use_coloc_pkg = TRUE)
    dpk <- max(abs(c(pk$H4_A - re$H4_A[re$draw == 1L], pk$H4_B - re$H4_B[re$draw == 1L])))
    say("RE-VERIFY (%s): %d checkpointed draws recomputed with the current code: max abs difference %.3g over %s; draw 1 coloc::coloc.abf vs abf_direct max |diff| %.3g",
        TAG, length(vd), dmax, paste(num_cols, collapse = ","), dpk)
    if (dmax > 1e-6 || dpk > 1e-6) stop("re-verification failed")
  }
  todo <- setdiff(seq_len(N_BOOT), D0$draw)
  ck <- file.path(WORK, "boot_draws.part_main.tsv")
} else {
  todo <- setdiff(seq(RANGE[1], RANGE[2]), D0$draw)
  ck <- file.path(WORK, sprintf("boot_draws.part_%d_%d.tsv", RANGE[1], RANGE[2]))
}
say("bootstrap (%s): %d draws to compute -> %s", TAG, length(todo), basename(ck))
buf <- list(); tl <- Sys.time(); ndone <- 0L
for (b in todo) {
  buf[[length(buf) + 1]] <- one_draw(b); ndone <- ndone + 1L
  if (length(buf) == 50L || b == tail(todo, 1)) {
    fwrite(rbindlist(buf), ck, sep = "\t", append = file.exists(ck) && file.size(ck) > 0)
    buf <- list()
    if (ndone %% 250L == 0L || b == tail(todo, 1))
      say("bootstrap progress (%s): %d/%d of this run's draws done (%.2f s/draw)", TAG, ndone, length(todo),
          as.numeric(difftime(Sys.time(), tl, units = "secs")) / ndone)
  }
}
if (!is.null(RANGE)) { say("==== S5_03_bootstrap worker %d:%d done in %.0f s ====", RANGE[1], RANGE[2], as.numeric(difftime(Sys.time(), t0, units = "secs"))); quit(save = "no") }

# ---------------------------------------------------------------- main: merge everything into work/boot_draws.tsv
D <- read_all()
D <- unique(D, by = "draw")[order(draw)]
stopifnot(nrow(D) == N_BOOT, !anyDuplicated(D$draw), setequal(D$draw, seq_len(N_BOOT)))
prov <- D[, .N, by = src]
fwrite(D[, ..COLS], file.path(WORK, "boot_draws.merged.tmp"), sep = "\t")
file.rename(file.path(WORK, "boot_draws.merged.tmp"), ck_main)
for (f in list.files(WORK, pattern = PART_PAT, full.names = TRUE)) file.rename(f, sub("\\.tsv$", ".merged_into_main.tsv", f))
say("bootstrap complete: %d draws (sources before merge: %s); invalid (a signal SNP monomorphic in the draw) %d; NA p (Q cysgly %d, Q eGFR %d, two %d); monomorphic p1np SNPs per draw median %s",
    nrow(D), paste(prov[, sprintf("%s %d", src, N)], collapse = "; "), sum(!D$valid_sig | !D$valid_ab), sum(is.na(D$p_Q_cysgly)), sum(is.na(D$p_Q_eGFR)),
    sum(is.na(D$p_two)), median(D$n_mono_set, na.rm = TRUE))
say("==== S5_03_bootstrap done in %.0f s ====", as.numeric(difftime(Sys.time(), t0, units = "secs")))
