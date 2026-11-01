"""w2_02c_fitset7.py -- Stage 2 step 2 (user decision 2026-10-09: s_max over ALL SEVEN distinct SuSiE fits at every locus).

The four-fit set of w2_02_locus.py / w2_04_susie_coloc.R is susie, p1fr, p1np on the 503 panel and p2np on the 489 panel.
Rule 3(i) reads "all variant sets, both panels": seven distinct fits = those four + p2np@503 + susie@489 + p1fr@489
(p1np@489 is p2np@489, the same SNP set on the same panel).

This step adds only what the three new fits need and changes nothing already computed:
  * panel-2 (489) validity for EVERY SNP of the susie set, palindromes included (w2_02_locus.py checks only the
    non-palindromic ones, because p2np = p1np-in-489). Same rule: same GRCh37 position and allele pair (unique), no
    missing genotype, panel-2 MAF >= 0.01, |A1 freq panel 2 - panel 1| <= 0.05;
  * PLINK --r square on all valid-489 SNPs (ld/R489_valid.f8, rows/cols in tgfb1_locus.tsv order, aligned to the 503 A1).
Asserts: valid489 AND in_p1np == in_p2np for every row (the existing p2np set is reproduced exactly), and the p2np
sub-matrix of R489_valid equals the existing R489_p2np.f8 (max |diff| < 1e-7).
Writes ld/R489_valid.f8, fitset7_sets.tsv (i, rsid503, rsid489, valid489, s489, af489).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402
import w2_02_locus as L  # noqa: E402


def main():
    W.log("==== w2_02c_fitset7 (TGFB1: panel-2 validity for every susie-set SNP; LD on the valid set) ====")
    tab = W.read_tsv(W.OUT + "/tgfb1_locus.tsv")
    bim489 = L.read_bim489()
    by37b, rc2 = {}, {}
    for b in bim489:
        by37b.setdefault(b["pos37"], []).append(b)
        rc2[b["rsid"]] = rc2.get(b["rsid"], 0) + 1
    cand = {}
    drops = {"panel2_absent": 0, "panel2_ambiguous": 0}
    for r in tab:
        a1, a2 = r["a1"], r["a2"]
        hits = [b for b in by37b.get(int(r["pos37"]), []) if {b["a1"], b["a2"]} == {a1, a2}]
        if len(hits) != 1 or rc2[hits[0]["rsid"]] != 1:
            drops["panel2_absent" if not hits else "panel2_ambiguous"] += 1
            continue
        cand[int(r["i"])] = hits[0]
    bfile489 = W.res("LD489_chr19_bed").replace(".bed", "")
    fr2 = {x["snp"]: x for x in L.plink_freq(bfile489, [b["rsid"] for b in cand.values()], L.LDD + "/p489_cand_all")}
    valid, s489, af489 = {}, {}, {}
    for r in tab:
        i = int(r["i"])
        if i not in cand:
            valid[i] = False
            continue
        b = cand[i]
        f = fr2[b["rsid"]]
        af2 = f["maf"] if f["a1"] == b["a1"] else 1 - f["maf"]
        sg = 1 if b["a1"] == r["a1"] else -1
        a503 = float(r["af503"])
        af2_on503 = af2 if sg == 1 else 1 - af2
        why = None
        if f["nchrobs"] != 978:
            why = "panel2_missing_genotype"
        elif min(af2, 1 - af2) < L.REF_MAF_MIN:
            why = "panel2_maf_lt_0.01"
        elif abs(af2_on503 - a503) > L.P2_AF_DIFF_MAX:
            why = "panel2_af_discordant"
        valid[i] = why is None
        s489[i], af489[i] = sg, af2_on503
        if why:
            drops[why] = drops.get(why, 0) + 1
    # the existing p2np set must be reproduced exactly
    for r in tab:
        i = int(r["i"])
        assert (valid[i] and r["in_p1np"] == "1") == (r["in_p2np"] == "1"), ("p2np not reproduced at row", i, r["rsid503"])
        if r["in_p2np"] == "1":
            assert cand[i]["rsid"] == r["rsid489"] and s489[i] == int(r["s489"]), ("489 rsid/sign differ at row", i)
    vi = [int(r["i"]) for r in tab if valid[int(r["i"])]]
    npal = sum(1 for r in tab if valid[int(r["i"])] and r["palindromic"] == "1")
    W.log("valid-489 SNPs (any susie-set SNP): %d of %d (palindromic %d); p2np reproduced exactly (%d); drops %s" % (
        len(vi), len(tab), npal, sum(1 for r in tab if r["in_p2np"] == "1"), drops))
    ids = [cand[i]["rsid"] for i in vi]
    order, Rraw, _ = L.plink_ld_square(bfile489, ids, L.LDD + "/ld489_valid")
    pos = {s: k for k, s in enumerate(order)}
    ix = [pos[s] for s in ids]
    sg = np.array([s489[i] for i in vi], dtype=float)
    R = Rraw[np.ix_(ix, ix)] * np.outer(sg, sg)
    assert np.isfinite(R).all() and np.allclose(np.diag(R), 1.0, atol=1e-6)
    # cross-check against the existing p2np matrix
    p = len(tab)
    i2 = [int(r["i"]) for r in tab if r["in_p2np"] == "1"]
    R2 = np.fromfile(W.OUT + "/ld/R489_p2np.f8", dtype="<f8").reshape(len(i2), len(i2))
    posv = {i: k for k, i in enumerate(vi)}
    sub = R[np.ix_([posv[i] for i in i2], [posv[i] for i in i2])]
    dmax = float(np.max(np.abs(sub - R2)))
    W.log("CHECK p2np sub-matrix of R489_valid vs existing R489_p2np.f8: max |diff| %.3g (%d x %d)" % (dmax, len(i2), len(i2)))
    assert dmax < 1e-7, dmax
    R.astype("<f8").tofile(L.LDD + "/R489_valid.f8")
    out = [{"i": int(r["i"]), "rsid503": r["rsid503"], "rsid489": cand[int(r["i"])]["rsid"] if int(r["i"]) in cand else "",
            "valid489": int(valid[int(r["i"])]), "s489": s489.get(int(r["i"]), ""), "af489": af489.get(int(r["i"]), float("nan"))} for r in tab]
    W.write_tsv(W.OUT + "/fitset7_sets.tsv", ["i", "rsid503", "rsid489", "valid489", "s489", "af489"], out)
    W.log("w2_02c done: ld/R489_valid.f8 md5 %s (%d x %d); fitset7_sets.tsv md5 %s" % (
        W.md5_file(L.LDD + "/R489_valid.f8"), R.shape[0], R.shape[1], W.md5_file(W.OUT + "/fitset7_sets.tsv")))


if __name__ == "__main__":
    main()
