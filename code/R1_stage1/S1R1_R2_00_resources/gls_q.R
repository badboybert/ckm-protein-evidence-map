# Shared GLS proportionality (Q) test for Stage 1 (pre-declaration §4(iv), §5, §6 UMOD).
# One implementation, sourced by W4b, W5 and W6b with an md5 assertion; verifiers re-implement it in Python.
# Method (as in the scoping script seen/wpD/wpD_03_multisignal.py, re-implemented here):
#   joint per-allele effects from marginal z with reference LD R:  zJ = R^-1 z,  bJ = zJ * s,  VJ = D R^-1 D  (D = diag(s))
#   for protein P and trait Y (independent samples): e(theta) = bY - theta * bP ;  V(theta) = VY + theta^2 VP
#   Q(theta) = e' V(theta)^-1 e, profiled over theta; Q_min ~ chi-square on k - 1 df under proportionality.
#   cov = "aware" keeps the LD-induced off-diagonal covariance; cov = "ign" uses only the diagonal of VY and VP.
# Run `Rscript gls_q.R --self-test` to check it.

joint_effects <- function(b, s, R) {
  stopifnot(length(b) == length(s), nrow(R) == length(b), ncol(R) == length(b), all(is.finite(c(b, s, R))), all(s > 0))
  Ri <- solve(R)
  z <- b / s
  D <- diag(s, nrow = length(s))
  list(b = as.vector(Ri %*% z) * s, V = D %*% Ri %*% D)
}

profile_Q <- function(bY, VY, bP, VP, idx = seq_along(bY), cov = c("aware", "ign"), se_order = c("2nd", "1st"),
                      bounds = c(-10, 10)) {
  # se_order "2nd" (Rule 5 primary): V(theta) = VY + theta^2 VP; "1st": the protein-side term is dropped (VP = 0).
  # There is no beta-space (COJO) option for the GLS Q: p_cojo is NOT RUN for Q tests (it applies to two-signal tests).
  cov <- match.arg(cov); se_order <- match.arg(se_order)
  bY <- bY[idx]; bP <- bP[idx]; VY <- VY[idx, idx, drop = FALSE]; VP <- VP[idx, idx, drop = FALSE]
  if (cov == "ign") { VY <- diag(diag(VY), nrow = length(idx)); VP <- diag(diag(VP), nrow = length(idx)) }
  if (se_order == "1st") VP <- VP * 0
  k <- length(idx); stopifnot(k >= 2)
  Qf <- function(th) { e <- bY - th * bP; V <- VY + th^2 * VP; as.numeric(t(e) %*% solve(V, e)) }
  # global search on a grid, then refine, so a boundary or a second minimum is not missed
  grid <- seq(bounds[1], bounds[2], length.out = 4001)
  qs <- vapply(grid, Qf, numeric(1))
  i <- which.min(qs)
  lo <- grid[max(1, i - 1)]; hi <- grid[min(length(grid), i + 1)]
  o <- optimize(Qf, c(lo, hi))
  at_bound <- (i == 1 || i == length(grid))
  list(theta = o$minimum, Q = o$objective, df = k - 1,
       p_upper = pchisq(o$objective, k - 1, lower.tail = FALSE),
       p_lower = pchisq(o$objective, k - 1, lower.tail = TRUE),
       ratios = bY / bP, cov = cov, se_order = se_order, at_bound = at_bound)
}

.self_test <- function() {
  set.seed(1)
  k <- 5
  R <- 0.2^abs(outer(1:k, 1:k, "-"))
  sP <- rep(0.01, k); sY <- rep(0.02, k)
  bPt <- c(-0.24, 0.53, 0.10, -0.15, 0.30)
  L <- chol(R)
  draw <- function(beta, s) { as.vector(R %*% beta) + as.vector(t(L) %*% rnorm(k)) * s }  # marginal = R * joint + noise
  pv <- replicate(400, {
    bP <- draw(bPt, sP); bY <- draw(-0.7 * bPt, sY)
    jp <- joint_effects(bP, sP, R); jy <- joint_effects(bY, sY, R)
    profile_Q(jy$b, jy$V, jp$b, jp$V)$p_upper
  })
  size <- mean(pv < 0.05)
  cat(sprintf("self-test: null size at 0.05 = %.3f (expect about 0.05)\n", size))
  stopifnot(size > 0.01, size < 0.10)
  # mutation: a non-proportional signal must be detected
  bYt <- -0.7 * bPt; bYt[1] <- bYt[1] + 0.5
  pw <- replicate(200, {
    bP <- draw(bPt, sP); bY <- draw(bYt, sY)
    jp <- joint_effects(bP, sP, R); jy <- joint_effects(bY, sY, R)
    profile_Q(jy$b, jy$V, jp$b, jp$V)$p_upper
  })
  cat(sprintf("self-test: power against a +0.5 deviation at signal 1 = %.3f (expect near 1)\n", mean(pw < 0.05)))
  stopifnot(mean(pw < 0.05) > 0.9)
  # mutation: the ign option must differ from aware when LD is strong
  R2 <- 0.8^abs(outer(1:k, 1:k, "-"))
  jp <- joint_effects(bPt, sP, R2); jy <- joint_effects(-0.7 * bPt + c(0.05, 0, 0, 0, 0), sY, R2)
  a <- profile_Q(jy$b, jy$V, jp$b, jp$V, cov = "aware")$Q; g <- profile_Q(jy$b, jy$V, jp$b, jp$V, cov = "ign")$Q
  cat(sprintf("self-test: aware Q %.3f vs ign Q %.3f under strong LD (must differ)\n", a, g))
  stopifnot(abs(a - g) > 1e-6)
  # mutation: first order (VP dropped) must give a larger Q than second order when the protein is noisy
  jpn <- joint_effects(bPt, rep(0.2, k), R); jyn <- joint_effects(-0.7 * bPt + c(0.2, 0, 0, 0, 0), sY, R)
  q2 <- profile_Q(jyn$b, jyn$V, jpn$b, jpn$V, se_order = "2nd")$Q; q1 <- profile_Q(jyn$b, jyn$V, jpn$b, jpn$V, se_order = "1st")$Q
  cat(sprintf("self-test: 1st-order Q %.3f vs 2nd-order Q %.3f with a noisy protein (1st must be larger)\n", q1, q2))
  stopifnot(q1 > q2)
  cat("SELF-TEST PASS\n")
}

if (sys.nframe() == 0L && "--self-test" %in% commandArgs(trailingOnly = TRUE)) .self_test()
