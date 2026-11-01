# exact_cond_or.R -- conditional-MLE odds ratio and its exact two-sided CI for a 2 x 2 table, solved to a tight tolerance.
# Added by the S6 FIXER (defect F1, 2026-10-08; see FIXLOG.md). One implementation, sourced by W6a (w6a_compute.R) and by
# W6b (R/07_fisher_d.R, md5 asserted there).
#
# The estimand and its defining equations are those of stats::fisher.test (R 4.4.1), unchanged. Given all margins, X = tab[1, 1]
# follows the noncentral hypergeometric distribution with odds ratio psi, and
#   conditional MLE:                E[X; psi_hat] = x
#   exact lower bound (two-sided):  P(X >= x; psi_lo) = (1 - conf.level) / 2
#   exact upper bound (two-sided):  P(X <= x; psi_hi) = (1 - conf.level) / 2
# fisher.test solves each equation with uniroot() at its default tolerance (.Machine$double.eps^0.25, about 1.2e-4) on psi or on
# 1/psi in (0, 1]. That leaves the upper bound wrong in the third decimal (11/16 vs 57/144: 12.8969 against 12.9014). Here each
# root is found on theta = log(psi) with uniroot(tol = tol), default 1e-13. The residual of every defining equation is returned,
# so callers can assert it. Edge cases follow fisher.test: x at the bottom of the support gives MLE 0 and lower bound 0; x at the
# top gives MLE Inf and upper bound Inf. The p value is not computed here (fisher.test's p involves no root finding).

.ecor_setup <- function(tab) {
  if (!is.matrix(tab) || !all(dim(tab) == 2L) || anyNA(tab) || any(tab < 0) || any(tab != round(tab)))
    stop("exact_cond_or: tab must be a 2 x 2 matrix of non-negative integer counts")
  m <- sum(tab[, 1]); n <- sum(tab[, 2]); k <- sum(tab[1, ]); x <- tab[1, 1]
  lo <- max(0, k - n); hi <- min(k, m)
  if (lo == hi) stop("exact_cond_or: degenerate support (a margin is zero)")
  support <- lo:hi
  logdc <- dhyper(support, m, n, k, log = TRUE)
  dn <- function(theta) { d <- logdc + theta * support; d <- exp(d - max(d)); d / sum(d) }
  list(x = x, lo = lo, hi = hi,
       mean_x = function(theta) sum(support * dn(theta)),
       p_ge   = function(theta) sum(dn(theta)[support >= x]),
       p_le   = function(theta) sum(dn(theta)[support <= x]))
}

# Residuals of the three defining equations at given values (est, lo, hi): E[X] - x, P(X >= x) - alpha, P(X <= x) - alpha.
# An edge value (0 or Inf) has residual 0 only where fisher.test would return it; anywhere else it is Inf.
exact_cond_or_resid <- function(tab, est, lo, hi, conf.level = 0.95) {
  s <- .ecor_setup(tab); alpha <- (1 - conf.level) / 2
  r1 <- function(v, f, target, edge_ok) if (is.finite(v) && v > 0) f(log(v)) - target else if (edge_ok) 0 else Inf
  c(mle = r1(est, s$mean_x, s$x, (s$x == s$lo && est == 0) || (s$x == s$hi && est == Inf)),
    lo  = r1(lo, s$p_ge, alpha, s$x == s$lo && lo == 0),
    hi  = r1(hi, s$p_le, alpha, s$x == s$hi && hi == Inf))
}

exact_cond_or <- function(tab, conf.level = 0.95, tol = 1e-13) {
  s <- .ecor_setup(tab); alpha <- (1 - conf.level) / 2
  th0 <- log((tab[1, 1] + .5) * (tab[2, 2] + .5) / ((tab[1, 2] + .5) * (tab[2, 1] + .5)))   # Haldane log OR: bracket centre
  solve <- function(f, target, dir) uniroot(function(t) f(t) - target, th0 + c(-1, 1), extendInt = dir, tol = tol, maxiter = 10000L)$root
  est <- if (s$x == s$lo) 0 else if (s$x == s$hi) Inf else exp(solve(s$mean_x, s$x, "upX"))   # E[X] increases with theta
  lo  <- if (s$x == s$lo) 0 else exp(solve(s$p_ge, alpha, "upX"))                            # P(X >= x) increases with theta
  hi  <- if (s$x == s$hi) Inf else exp(solve(s$p_le, alpha, "downX"))                        # P(X <= x) decreases with theta
  list(estimate = est, conf.int = c(lo, hi), conf.level = conf.level, tol = tol,
       resid = exact_cond_or_resid(tab, est, lo, hi, conf.level))
}
