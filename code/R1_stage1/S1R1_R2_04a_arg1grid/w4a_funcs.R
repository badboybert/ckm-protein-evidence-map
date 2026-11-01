# w4a_funcs.R -- functions for Stage 1 W4a (pre-declaration rev 4.2 section 4 (i)(ii)(iii)(vi); amendment 1 A2).
# Pure functions plus logging/md5 helpers. Sourced by w4a_main.R. No side effects on load.
#
# Method (implemented here from the pre-declaration, the original rule S1R1_ARG1_PRESTATED_RULE.md and the
# l3fm__conditional.py docstring (engine, read only); the R1 triangle script (seen_script) was NOT opened):
#   z-space conditioning (Yang 2012 z-approximation, equal N, standardised genotypes), two SNPs A and B with LD r:
#     z_A|B = (z_A - r z_B)/sqrt(1 - r^2),  se_A|B = se_A/sqrt(1 - r^2),  b_A|B = z_A|B * se_A|B   (and symmetrically)
#     corr(b_A|B, b_B|A) = -r  (joint estimates: V = D R^-1 D, as gls_q.R joint_effects)
#   marginal estimates at A and B from one GWAS sample: corr(b_A, b_B) = r
#   beta-space (COJO) conditioning: b_A|B = (b_A - r sqrt(v_B/v_A) b_B)/(1 - r^2), v = 2f(1-f) from the source's own
#     effect-allele frequency (COJO .ma convention), se_A|B = se_A/sqrt(1 - r^2), corr = -r
#   ratio R_X = b_Y,X / b_P,X ; delta-method variance, second order (Rule 5 primary):
#     Var(R_X) = se_Y^2/b_P^2 + b_Y^2 se_P^2/b_P^4        ("1st" order drops the protein-side term)
#     Cov(R_A,R_B) = c_Y se_YA se_YB/(b_PA b_PB) + R_A R_B c_P se_PA se_PB/(b_PA b_PB)   (aware; 0 when "ign")
#     the protein (or arginine) and the outcome come from different samples, so they carry no cross-covariance.
#   A-vs-B test: d = R_A - R_B, z = d / sqrt(Var(R_A) + Var(R_B) - 2 Cov), two-sided normal p.

.W4A_LOG <- NULL
log_init <- function(path) { .W4A_LOG <<- path; invisible(path) }
logf <- function(fmt, ...) {
  s <- if (length(list(...))) sprintf(fmt, ...) else fmt
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", s)
  cat(line, "\n", sep = "")
  if (!is.null(.W4A_LOG)) cat(line, "\n", sep = "", file = .W4A_LOG, append = TRUE)
  invisible(s)
}
md5f <- function(p) unname(as.character(tools::md5sum(p)))
assert_md5 <- function(p, want, label) {
  if (!file.exists(p)) stop(sprintf("MISSING %s: %s", label, p))
  got <- md5f(p)
  ok <- identical(got, want)
  logf("md5 %s %s  bytes=%s  %s  [%s]", if (ok) "OK      " else "MISMATCH", got, format(file.size(p), scientific = FALSE), p, label)
  if (!ok) stop(sprintf("STOP: md5 mismatch for %s: got %s want %s", label, got, want))
  invisible(got)
}

# ---------------- estimators ----------------
cond_z <- function(bA, sA, bB, sB, r) {
  zA <- bA / sA; zB <- bB / sB; k <- sqrt(1 - r^2)
  zAc <- (zA - r * zB) / k; zBc <- (zB - r * zA) / k
  sAc <- sA / k; sBc <- sB / k
  list(bA = zAc * sAc, sA = sAc, bB = zBc * sBc, sB = sBc, c = -r, zA = zAc, zB = zBc)
}
cond_beta <- function(bA, sA, bB, sB, r, fA, fB) {
  vA <- 2 * fA * (1 - fA); vB <- 2 * fB * (1 - fB)
  bAc <- (bA - r * sqrt(vB / vA) * bB) / (1 - r^2)
  bBc <- (bB - r * sqrt(vA / vB) * bA) / (1 - r^2)
  sAc <- sA / sqrt(1 - r^2); sBc <- sB / sqrt(1 - r^2)
  list(bA = bAc, sA = sAc, bB = bBc, sB = sBc, c = -r, zA = bAc / sAc, zB = bBc / sBc)
}
marg_est <- function(bA, sA, bB, sB, r) list(bA = bA, sA = sA, bB = bB, sB = sB, c = r, zA = bA / sA, zB = bB / sB)

# estimates for one trait under an estimand; vectorised over r
est_trait <- function(t, r, estimand, space = "z") {
  if (estimand == "marg") return(marg_est(t$bA, t$sA, t$bB, t$sB, r))
  if (space == "z") return(cond_z(t$bA, t$sA, t$bB, t$sB, r))
  if (space == "beta") return(cond_beta(t$bA, t$sA, t$bB, t$sB, r, t$fA, t$fB))
  stop("space")
}

ab_test <- function(Y, X, order = "2nd", aware = TRUE) {
  # Y, X: lists with bA, sA, bB, sB, c (correlation of the A and B estimates within that sample)
  RA <- Y$bA / X$bA; RB <- Y$bB / X$bB
  o2 <- as.numeric(order == "2nd")
  vA <- Y$sA^2 / X$bA^2 + o2 * Y$bA^2 * X$sA^2 / X$bA^4
  vB <- Y$sB^2 / X$bB^2 + o2 * Y$bB^2 * X$sB^2 / X$bB^4
  cv <- if (aware) Y$c * Y$sA * Y$sB / (X$bA * X$bB) + o2 * RA * RB * X$c * X$sA * X$sB / (X$bA * X$bB) else 0 * RA
  d <- RA - RB; vd <- vA + vB - 2 * cv
  se <- sqrt(vd); z <- d / se
  q <- qnorm(0.975)
  list(RA = RA, RB = RB, vA = vA, vB = vB, cov = cv, d = d, se = se, lo = d - q * se, hi = d + q * se, z = z,
       p = 2 * pnorm(-abs(z)))
}

fieller <- function(a, b, v11, v22, v12, level = 0.95) {
  t2 <- qnorm(1 - (1 - level) / 2)^2
  A <- b^2 - t2 * v22; B <- -2 * (a * b - t2 * v12); C <- a^2 - t2 * v11
  disc <- B^2 - 4 * A * C
  if (!(A > 0) || disc < 0) return(list(lo = NA_real_, hi = NA_real_, status = "unbounded (denominator not significantly non-zero)"))
  rt <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
  list(lo = rt[1], hi = rt[2], status = "bounded")
}
logdelta_ci <- function(a, b, v11, v22, v12) {
  ror <- a / b
  if (ror <= 0) return(c(NA_real_, NA_real_))
  v <- v11 / a^2 + v22 / b^2 - 2 * v12 / (a * b)
  exp(log(ror) + c(-1, 1) * qnorm(0.975) * sqrt(v))
}

power_calc <- function(RA, vRA, xB, sxB, syB) {
  # expected z at B if B followed A's ratio, given the outcome's actual SE at B (original rule)
  q <- qnorm(0.975)
  ez <- RA * xB / syB
  pp <- pnorm(abs(ez) - q) + pnorm(-abs(ez) - q)
  # predictive: Z_B ~ N(ez, 1 + (xB^2 Var(R_A) + R_A^2 se_PB^2)/se_YB^2)
  sdv <- sqrt(1 + (xB^2 * vRA + RA^2 * sxB^2) / syB^2)
  pr <- 1 - pnorm((q - ez) / sdv) + pnorm((-q - ez) / sdv)
  list(exp_z = ez, power_point = pp, power_pred = pr, pred_sd = sdv)
}

# conditioning a whole window on ONE SNP k (engine convention: drop k and SNPs with r^2 > 0.9 to k)
cond_window <- function(b, s, r_to_k, k, collin = 0.9) {
  z <- b / s; zk <- z[k]
  q <- r_to_k^2
  zc <- (z - r_to_k * zk) / sqrt(pmax(1 - q, 1e-12))
  sc <- s / sqrt(pmax(1 - q, 1e-12))
  bad <- q > collin; bad[k] <- TRUE
  zc[bad] <- NA; sc[bad] <- NA
  list(b = zc * sc, s = sc, z = zc, n_collinear = sum(q > collin) - 1)
}

# ---------------- liftover (UCSC chain, single positions, + strand chains only) ----------------
read_chain_blocks <- function(path, tname, qname) {
  L <- readLines(gzfile(path))
  hdr <- which(startsWith(L, "chain"))
  hdr_end <- c(hdr[-1] - 1, length(L))
  out <- list()
  for (h in seq_along(hdr)) {
    f <- strsplit(L[hdr[h]], " ")[[1]]
    if (f[3] != tname || f[8] != qname) next
    if (f[5] != "+" || f[10] != "+") next
    score <- as.numeric(f[2]); tS <- as.numeric(f[6]); qS <- as.numeric(f[11])
    body <- L[(hdr[h] + 1):hdr_end[h]]
    body <- body[nzchar(body)]
    m <- do.call(rbind, lapply(strsplit(body, "[ \t]"), function(x) { x <- as.numeric(x); c(x, rep(0, 3 - length(x))) }))
    tt <- tS; qq <- qS
    tb <- numeric(nrow(m)); qb <- numeric(nrow(m))
    for (i in seq_len(nrow(m))) { tb[i] <- tt; qb[i] <- qq; tt <- tt + m[i, 1] + m[i, 2]; qq <- qq + m[i, 1] + m[i, 3] }
    out[[length(out) + 1]] <- data.frame(score = score, t0 = tb, t1 = tb + m[, 1], q0 = qb)
  }
  b <- do.call(rbind, out)
  b[order(-b$score, b$t0), ]
}
lift_pos <- function(pos1, blocks) {
  vapply(pos1, function(p) {
    p0 <- p - 1
    i <- which(blocks$t0 <= p0 & p0 < blocks$t1)
    if (!length(i)) return(NA_real_)
    i <- i[1]
    blocks$q0[i] + (p0 - blocks$t0[i]) + 1
  }, numeric(1))
}

# ---------------- PLINK ----------------
run_plink <- function(exe, args, label) {
  io <- which(args == "--out")
  if (length(io)) { pre <- gsub('^"|"$', "", args[io + 1]); unlink(paste0(pre, c(".ld", ".frq", ".log", ".nosex"))) }  # no stale outputs
  t0 <- Sys.time()
  out <- suppressWarnings(system2(exe, args, stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); st <- if (is.null(st)) 0L else st
  logf("PLINK [%s] exit %d in %.1f s: %s", label, st, as.numeric(difftime(Sys.time(), t0, units = "secs")),
       paste(args, collapse = " "))
  w <- grep("Warning|Error", out, value = TRUE)
  if (length(w)) logf("PLINK [%s] messages: %s", label, paste(w, collapse = " | "))
  list(status = st, out = out)
}
read_ld_table <- function(path) {
  if (!file.exists(path)) stop(sprintf("LD table missing: %s", path))
  x <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  x$R <- suppressWarnings(as.numeric(x$R))
  x
}
assert_ld_complete <- function(tab, list_ids, set_ids, label) {
  for (k in list_ids) {
    s <- tab[tab$SNP_A == k, ]
    miss <- setdiff(set_ids, s$SNP_B)
    if (length(miss)) stop(sprintf("LD ASSERT [%s]: %d requested pairs with %s not returned (e.g. %s)", label, length(miss), k,
                                   paste(head(miss, 3), collapse = ",")))
    if (anyDuplicated(s$SNP_B)) stop(sprintf("LD ASSERT [%s]: duplicated pairs for %s", label, k))
    s <- s[s$SNP_B %in% set_ids, ]
    if (any(!is.finite(s$R))) stop(sprintf("LD ASSERT [%s]: %d NA/NaN r values for %s", label, sum(!is.finite(s$R)), k))
  }
  invisible(TRUE)
}

# ---------------- coloc ----------------
coloc_pair <- function(bP, sP, nP, bY, sY, nY, typeY, sfrac = NA, ids, p12 = 1e-5, sdY_trait = 1) {
  ok <- is.finite(bP) & is.finite(sP) & is.finite(bY) & is.finite(sY)
  d1 <- list(beta = bP[ok], varbeta = sP[ok]^2, type = "quant", sdY = 1, N = nP, snp = ids[ok])
  d2 <- if (typeY == "quant") list(beta = bY[ok], varbeta = sY[ok]^2, type = "quant", sdY = sdY_trait, N = nY, snp = ids[ok]) else
    list(beta = bY[ok], varbeta = sY[ok]^2, type = "cc", s = sfrac, N = nY, snp = ids[ok])
  res <- NULL
  invisible(capture.output(res <- suppressMessages(coloc::coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12))))
  s <- res$summary
  list(nsnp = unname(s["nsnps"]), H0 = unname(s["PP.H0.abf"]), H1 = unname(s["PP.H1.abf"]), H2 = unname(s["PP.H2.abf"]),
       H3 = unname(s["PP.H3.abf"]), H4 = unname(s["PP.H4.abf"]))
}

# ---------------- formatting ----------------
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, sprintf("%.2e", p), sprintf("%.3f", p)))
fmt3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
fmt4 <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))
