#!/usr/bin/env python
"""W3 step 2: implementer self-check (NOT the §10 verifier). Recomputes every coloc.abf posterior in results.tsv
from coloc_input_<pair>.tsv with a from-scratch numpy implementation of Wakefield ABFs, so an error in how the R
datasets were assembled (wrong MAF, N, s, varbeta, prior sd, SNP set) would show as a difference.

Approximations as in R coloc 5.2.3 (documented there): quantitative trait with beta/varbeta and sdY -> prior sd 0.15*sdY;
case-control with p + MAF + N + s -> z = qnorm(p/2), V = 1 / (2 N f (1-f) s (1-s)), prior sd 0.2;
case-control with beta/varbeta -> V = varbeta, prior sd 0.2.
"""
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import csv, math, os, sys, time
import numpy as np
from scipy.stats import norm
from scipy.special import logsumexp

OUT = (_P31R_V3 + "/R1_stage1/S1R1_R2_03_extcoloc")
NINFO = {}  # parsed from results.tsv N_eff notes, not typed


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S") + " [selfcheck] " + msg
    print(line, flush=True)
    with open(os.path.join(OUT, "run.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def labf(z, V, W):
    r = W / (V + W)
    return 0.5 * (np.log(1 - r) + r * z ** 2)


def pps(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5):
    lh1, lh2, lh4 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    lsum = lh1 + lh2
    lh3 = lsum + np.log1p(-np.exp(lh4 - lsum))
    parts = np.array([0.0, math.log(p1) + lh1, math.log(p2) + lh2, math.log(p1) + math.log(p2) + lh3, math.log(p12) + lh4])
    return np.exp(parts - logsumexp(parts))


def main():
    import scipy
    log("python %s; numpy %s; scipy %s" % (sys.version.split()[0], np.__version__, scipy.__version__))
    res = {}
    with open(os.path.join(OUT, "results.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            res[r["quantity_id"]] = r
    for k, r in res.items():
        if k.startswith("S3.N_eff."):
            oid = k[len("S3.N_eff."):]
            import re
            m = re.search(r"cases (\d+), controls (\d+)", r["notes"])
            ca, co = int(m.group(1)), int(m.group(2))
            NINFO[oid] = (ca + co, ca / (ca + co))
    pairs = {"PAM_T2D": "OUT_DIAMANTE2018_noUKB_EUR_chrpos", "RAD23B_T2D": "OUT_DIAMANTE2018_noUKB_EUR_chrpos",
             "TGFB1_CAD": "OUT_CARDIOGRAMplusC4D_ieu-a-7", "PCSK9_CAD": "OUT_CARDIOGRAMplusC4D_ieu-a-7"}
    worst, n_cmp = 0.0, 0
    for pair, oid in pairs.items():
        N, s = NINFO[oid]
        rows = list(csv.DictReader(open(os.path.join(OUT, "coloc_input_%s.tsv" % pair), encoding="utf-8"), delimiter="\t", quoting=csv.QUOTE_NONE))
        assert len(rows) > 1000, (pair, len(rows))
        f = lambda c: np.array([float(x[c]) for x in rows])
        bx, sx = f("BETA"), f("SE")
        l1 = labf(bx / sx, sx ** 2, 0.15 ** 2)
        maf_out, p_out = f("maf_out"), f("p_out")
        by, sy = f("beta_out_al"), f("se_out")
        src = np.array([x["maf_src"] for x in rows])
        for inp in ("pmaf", "betase", "pmafdrop"):
            keep = np.ones(len(rows), bool) if inp != "pmafdrop" else (src == "outcome")
            if inp == "pmafdrop" and keep.all():
                continue
            if inp in ("pmaf", "pmafdrop"):
                z = norm.isf(p_out / 2)
                V = 1.0 / (2 * N * maf_out * (1 - maf_out) * s * (1 - s))
            else:
                z, V = by / sy, sy ** 2
            l2 = labf(z, V, 0.2 ** 2)
            for p12, tag in ((1e-5, "1e-05"), (1e-6, "1e-06")):
                pp = pps(l1[keep], l2[keep], p12=p12)
                for h, v in zip(("H0", "H1", "H2", "H3", "H4"), pp):
                    qid = "S3.%s.%s.%s.%s" % (pair, inp, tag, h)
                    rv = float(res[qid]["value"])
                    d = abs(rv - v)
                    worst = max(worst, d)
                    n_cmp += 1
                    if d > 1e-5:
                        log("SELF-CHECK DIFF %s: R %g vs numpy %g" % (qid, rv, v))
                nq = "S3.%s.%s.%s.nsnp" % (pair, inp, tag)
                assert int(res[nq]["value"]) == int(keep.sum()), (nq, res[nq]["value"], keep.sum())
    log("numpy re-implementation vs R coloc 5.2.3: %d posteriors compared, max |diff| = %.3g (values printed to 6 s.f.); nsnp equal on every run"
        % (n_cmp, worst))
    if worst > 1e-5:
        sys.exit(4)


if __name__ == "__main__":
    main()
