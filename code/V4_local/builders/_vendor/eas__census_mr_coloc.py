"""Core MR + colocalization functions — Python reimplementation (no R).
Wald ratio, IVW (random-effects), Steiger directionality, and coloc.abf.
coloc.abf follows Giambartolomei et al. 2014 (Wakefield ABF). Validate against
R/coloc on a known locus before trusting at scale.
"""
import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def fstat(beta, se):
    return (np.asarray(beta) / np.asarray(se)) ** 2


def harmonise(ea_x, oa_x, ea_y, oa_y, by, byse, palindromic_drop=True):
    """Align an outcome effect to the exposure effect allele.
    Returns (by_aligned, byse) or None if undecidable. Alleles uppercase single-char.
    """
    ea_x, oa_x, ea_y, oa_y = (a.upper() for a in (ea_x, oa_x, ea_y, oa_y))
    palindromic = ea_x == COMPLEMENT.get(oa_x)
    if ea_x == ea_y and oa_x == oa_y:
        return by, byse
    if ea_x == oa_y and oa_x == ea_y:                       # simple allele swap
        return -by, byse
    if ea_x == COMPLEMENT.get(ea_y) and oa_x == COMPLEMENT.get(oa_y):   # strand flip
        return None if (palindromic and palindromic_drop) else (by, byse)
    if ea_x == COMPLEMENT.get(oa_y) and oa_x == COMPLEMENT.get(ea_y):   # flip + swap
        return None if (palindromic and palindromic_drop) else (-by, byse)
    return None                                             # incompatible alleles


def wald_ratio(bx, by, byse):
    beta = by / bx
    se = abs(byse / bx)                                     # first-order delta method
    return beta, se, 2 * norm.sf(abs(beta / se))


def ivw(bx, bxse, by, byse):
    """Inverse-variance weighted MR with multiplicative random-effects SE."""
    bx, bxse, by, byse = (np.asarray(v, float) for v in (bx, bxse, by, byse))
    w = 1.0 / byse ** 2
    beta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = np.sqrt(1.0 / np.sum(w * bx ** 2))
    resid = by - beta * bx
    q = np.sum(w * resid ** 2)
    dof = max(len(bx) - 1, 1)
    se *= max(1.0, np.sqrt(q / dof))                        # inflate for heterogeneity
    return beta, se, 2 * norm.sf(abs(beta / se))


def steiger(bx, bxse, n_x, by, byse, n_y):
    """Directional check via variance explained (z^2/(z^2+n)). Approximate for binary
    outcomes (use as a flag, not a hard filter). Returns (correct_dir, r2_exp, r2_out)."""
    zx, zy = bx / bxse, by / byse
    r2x = zx ** 2 / (zx ** 2 + n_x)
    r2y = zy ** 2 / (zy ** 2 + n_y)
    return bool(r2x > r2y), float(r2x), float(r2y)


def _labf(beta, se, W):
    V = np.asarray(se, float) ** 2
    z = np.asarray(beta, float) / np.asarray(se, float)
    r = W / (V + W)
    return 0.5 * (np.log(1 - r) + r * z ** 2)


def coloc_abf(beta1, se1, beta2, se2, sd1=0.15, sd2=0.15,
              p1=1e-4, p2=1e-4, p12=1e-5):
    """coloc.abf over one region; arrays aligned SNP-by-SNP. Returns PP.H0..H4."""
    l1 = _labf(beta1, se1, sd1 ** 2)
    l2 = _labf(beta2, se2, sd2 ** 2)
    lh1, lh2 = logsumexp(l1), logsumexp(l2)
    lh4 = logsumexp(l1 + l2)
    # H3 (different causal variants) = log(exp(lh1+lh2) - exp(lh4))
    lh3 = (lh1 + lh2) + np.log1p(-np.exp(lh4 - (lh1 + lh2)))
    parts = np.array([0.0,
                      np.log(p1) + lh1,
                      np.log(p2) + lh2,
                      np.log(p1) + np.log(p2) + lh3,
                      np.log(p12) + lh4])
    pp = np.exp(parts - logsumexp(parts))
    return dict(zip(["H0", "H1", "H2", "H3", "H4"], pp.tolist()))
