"""w2_03_cond.py -- TGFB1 two-instrument analysis (pre-declaration section 2; Rule 5; Rule 3(ii)).

Signals: A = rs73045269 (discovery instrument), B = rs12983775 (CKB instrument). All betas on the 503 bim A1 allele
(rs73045269:T, rs12983775:G); r(A, B) from PLINK --r on those alleles.

Joint / conditional (z-space; Yang et al. 2012 z-approximation, as in the l3fm engine): with z = (zA, zB), s = (sA, sB),
  bJ = s * R^-1 z,  V(bJ) = D R^-1 D (D = diag s).  For two SNPs bJ_A = b(A | B) and bJ_B = b(B | A) exactly.
Ratio per signal X: rho_X = bJ_Y,X / bJ_P,X (CAD log-odds per SD of TGFB1), protein and CAD from independent samples.
A-vs-B test of d = rho_A - rho_B, delta method:
  p_aware_2nd (primary, Rule 5): Var(rho_X) = V_Y,XX / bP_X^2 + bY_X^2 V_P,XX / bP_X^4,
                                 Cov(rho_A, rho_B) = V_Y,AB / (bP_A bP_B) + bY_A bY_B V_P,AB / (bP_A^2 bP_B^2)
  p_ign        : the same with the LD-induced covariances V_.,AB set to 0 (diagonal only), second order
  p_aware_1st  : protein-side terms dropped (V_P = 0), covariance kept
  p_cojo       : beta-space (COJO-style) joint estimates -- bJ = (H^1/2 R H^1/2)^-1 H b with H = diag(2 f (1 - f)) from
                 each resource's own reported allele frequency, V(bJ) = c (H^1/2 R H^1/2)^-1 with c = the regional median
                 of h_j s_j^2 (QC'd window SNPs, MAF > 0.05) -- then the aware, second-order test.
Bootstrap (Rule 3(ii)): 5,000 draws of the 503 people with replacement (numpy default_rng(20261007), PCG64), the .bed
read in place (DEVIATION: bootstrap outside PLINK (section 10)); per draw r(A, B), the joint fits, the ratios and the
aware second-order p; conclusion = side of p 0.05 and sign of d; share = fraction of draws with the full-panel
conclusion. A draw where A or B is monomorphic gets r = 0 (counted).
Also writes the conditioned window z/se (each signal isolated by conditioning on the other; SNPs with r^2 > 0.9 to the
conditioning SNP dropped, as in the engine) on 503 (p1np) and 489 (p2np) for the Rule 3(iv) colocalization in R.
"""
import json
import math
import os
import sys

import numpy as np
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

A_RS, B_RS = "rs73045269", "rs12983775"
TR = {"protein": "P", "finngen_chd": "CHD_FG", "ieu_a_7": "CAD_C4D"}
OUTC = ["finngen_chd", "ieu_a_7"]
NBOOT = 5000
COLLIN = 0.9


def p2(z):
    return float(2 * norm.sf(abs(z)))


def joint_z(b, s, r):
    R = np.array([[1.0, r], [r, 1.0]])
    Ri = np.linalg.inv(R)
    z = b / s
    bJ = (Ri @ z) * s
    V = np.diag(s) @ Ri @ np.diag(s)
    return bJ, V


def joint_cojo(b, s, f, r, c):
    h = 2 * f * (1 - f)
    Hh = np.diag(np.sqrt(h))
    M = Hh @ np.array([[1.0, r], [r, 1.0]]) @ Hh
    Mi = np.linalg.inv(M)
    bJ = Mi @ (h * b)
    return bJ, c * Mi


def ratio_test(bY, VY, bP, VP, cov=True, second=True):
    rho = bY / bP
    gY = 1.0 / bP                       # d rho_X / d bY_X
    gP = -bY / bP ** 2                  # d rho_X / d bP_X
    VYm = VY.copy() if cov else np.diag(np.diag(VY))
    VPm = VP.copy() if cov else np.diag(np.diag(VP))
    if not second:
        VPm = VPm * 0
    C = np.outer(gY, gY) * VYm + np.outer(gP, gP) * VPm     # covariance of (rho_A, rho_B)
    d = rho[0] - rho[1]
    vd = C[0, 0] + C[1, 1] - 2 * C[0, 1]
    z = d / math.sqrt(vd)
    return {"rhoA": rho[0], "rhoB": rho[1], "seA": math.sqrt(C[0, 0]), "seB": math.sqrt(C[1, 1]), "cov": C[0, 1],
            "d": d, "se_d": math.sqrt(vd), "z": z, "p": p2(z)}


def bed_rows(bed, n, idxs):
    nb = (n + 3) // 4
    lut = np.array([2.0, np.nan, 1.0, 0.0])
    G = np.empty((n, len(idxs)))
    with open(bed, "rb") as fh:
        assert fh.read(3) == b"\x6c\x1b\x01"
        for j, i in enumerate(idxs):
            fh.seek(3 + i * nb)
            by = np.frombuffer(fh.read(nb), dtype=np.uint8)
            G[:, j] = lut[np.stack([(by >> s) & 3 for s in (0, 2, 4, 6)], axis=1).reshape(-1)[:n]]
    assert not np.isnan(G).any()
    return G


def main():
    W.log("==== w2_03_cond (TGFB1 two-instrument analysis) ====")
    meta = json.load(open(W.OUT + "/tgfb1_locus_meta.json"))
    W.assert_md5(W.OUT + "/tgfb1_locus.tsv", meta["md5"]["tgfb1_locus.tsv"], "tgfb1_locus.tsv (own)")
    rows = W.read_tsv(W.OUT + "/tgfb1_locus.tsv")
    nS = len(rows)
    R503 = np.fromfile(W.OUT + "/ld/R503_susie.f8", dtype="<f8").reshape(nS, nS)
    p2rows = [i for i, r in enumerate(rows) if r["in_p2np"] == "1"]
    R489 = np.fromfile(W.OUT + "/ld/R489_p2np.f8", dtype="<f8").reshape(len(p2rows), len(p2rows))
    rs = [r["rsid503"] for r in rows]
    iA, iB = rs.index(A_RS), rs.index(B_RS)
    jA, jB = p2rows.index(iA), p2rows.index(iB)
    r503, r489 = float(R503[iA, iB]), float(R489[jA, jB])
    W.log("A = %s (%s, A1 %s), B = %s (%s, A1 %s); r(A,B) PLINK 503 = %.6f, 489 = %.6f (A1 coding)" % (
        A_RS, rows[iA]["snp"], rows[iA]["a1"], B_RS, rows[iB]["snp"], rows[iB]["a1"], r503, r489))
    # ---------------- per-resource marginal values at A and B (on the 503 A1)
    marg = {}
    for res_name, t in TR.items():
        b = np.array([float(rows[iA]["b__" + t]), float(rows[iB]["b__" + t])])
        s = np.array([float(rows[iA]["se__" + t]), float(rows[iB]["se__" + t])])
        f = np.array([float(rows[iA]["fA1__" + t]), float(rows[iB]["fA1__" + t])])
        if np.isnan(f).any():
            f = np.array([float(rows[iA]["af503"]), float(rows[iB]["af503"])])
        marg[res_name] = (b, s, f)
    # COJO residual-variance constant per resource (regional median of h s^2)
    cconst = {}
    for res_name, t in TR.items():
        v = []
        for r in rows:
            f = float(r["fA1__" + t])
            if f != f:
                f = float(r["af503"])
            if 0.05 < min(f, 1 - f):
                v.append(2 * f * (1 - f) * float(r["se__" + t]) ** 2)
        cconst[res_name] = float(np.median(v))
    out = []

    def add(qid, estimand, resource, panel, value, se="", p="", cov="", se_order="", unit="", notes="", n="", ci=("", "")):
        out.append({"quantity_id": qid, "prespec_ref": "s2", "estimand": estimand, "resource": resource, "panel": panel,
                    "cov": cov, "se_order": se_order, "value": value, "se": se, "ci_lo": ci[0], "ci_hi": ci[1], "p": p,
                    "n": n, "unit": unit, "notes": notes})
    add("S2.r.503", "r(rs73045269:T, rs12983775:G) A1-allele dosage correlation", "LD503 (1000G EUR, g1000_eur)", "503", r503,
        n=503, unit="r", notes="PLINK19_paper6 --r square --keep-allele-order; on the protein-raising alleles (T, A) r = %+.6f" % (-r503))
    add("S2.r.489", "r(rs73045269:T, rs12983775:G) A1-allele dosage correlation", "LD489 (1000G.EUR.QC chr19)", "489", r489,
        n=489, unit="r", notes="PLINK19_paper6 --r square --keep-allele-order; 489 A1 equals 503 A1 at both SNPs")
    disc = {}
    full = {}
    for panel, r in (("503", r503), ("489", r489)):
        J = {k: joint_z(*marg[k][:2], r) for k in TR}
        Jc = {k: joint_cojo(marg[k][0], marg[k][1], marg[k][2], r, cconst[k]) for k in TR}
        pre = "S2" if panel == "503" else "S2.p489"
        for k in TR:
            bJ, V = J[k]
            for x, inst in ((0, A_RS), (1, B_RS)):
                other = B_RS if x == 0 else A_RS
                allele = rows[iA if x == 0 else iB]["a1"]
                add("%s.cond.%s.%s.beta" % (pre, k, inst), "joint (= conditional) effect of %s:%s given %s (z-space, r from %s)" % (
                    inst, allele, other, panel), k, panel, float(bJ[x]), se=math.sqrt(V[x, x]), p=p2(bJ[x] / math.sqrt(V[x, x])),
                    unit="per allele (%s)" % ("SD protein" if k == "protein" else "log-odds"),
                    notes="z_cond = %.4f; marginal beta %.6g se %.6g" % (bJ[x] / math.sqrt(V[x, x]), marg[k][0][x], marg[k][1][x]))
                add("%s.cond.%s.%s.se" % (pre, k, inst), "SE of the joint effect of %s:%s given %s" % (inst, allele, other), k, panel,
                    math.sqrt(V[x, x]))
        for k in OUTC:
            bY, VY = J[k]
            bP, VP = J["protein"]
            tests = {"p_aware_2nd": ratio_test(bY, VY, bP, VP, True, True), "p_ign": ratio_test(bY, VY, bP, VP, False, True),
                     "p_aware_1st": ratio_test(bY, VY, bP, VP, True, False)}
            bYc, VYc = Jc[k]
            bPc, VPc = Jc["protein"]
            tests["p_cojo"] = ratio_test(bYc, VYc, bPc, VPc, True, True)
            t0 = tests["p_aware_2nd"]
            for x, inst in ((0, A_RS), (1, B_RS)):
                add("%s.ratio.%s.%s.value" % (pre, k, inst), "conditional ratio: CAD log-odds per SD TGFB1 at %s (mutually conditioned)" % inst,
                    k, panel, t0["rhoA"] if x == 0 else t0["rhoB"], se=t0["seA"] if x == 0 else t0["seB"], cov="aware",
                    se_order="2nd", unit="log-odds per SD protein",
                    ci=((t0["rhoA"] if x == 0 else t0["rhoB"]) - 1.959964 * (t0["seA"] if x == 0 else t0["seB"]),
                        (t0["rhoA"] if x == 0 else t0["rhoB"]) + 1.959964 * (t0["seA"] if x == 0 else t0["seB"])),
                    p=p2((t0["rhoA"] if x == 0 else t0["rhoB"]) / (t0["seA"] if x == 0 else t0["seB"])))
                add("%s.ratio.%s.%s.se" % (pre, k, inst), "SE (delta method, second order) of the conditional ratio at %s" % inst, k, panel,
                    t0["seA"] if x == 0 else t0["seB"], se_order="2nd")
            for pk, t in tests.items():
                add("%s.het.%s.%s" % (pre, k, pk), "A-vs-B difference of conditional ratios rho(rs73045269) - rho(rs12983775)", k, panel,
                    t["d"], se=t["se_d"], p=t["p"], cov="aware" if pk in ("p_aware_2nd", "p_aware_1st", "p_cojo") else "ign",
                    se_order="1st" if pk == "p_aware_1st" else "2nd",
                    ci=(t["d"] - 1.959964 * t["se_d"], t["d"] + 1.959964 * t["se_d"]),
                    unit="log-odds per SD protein", notes=("PRIMARY (Rule 5)" if pk == "p_aware_2nd" else "") + (
                        "; beta-space (COJO-style) conditioning, c = %.4g (protein) %.4g (outcome)" % (cconst["protein"], cconst[k]) if pk == "p_cojo" else "")
                    + ("; family entry (Rule 4)" if (k == "finngen_chd" and pk == "p_aware_2nd" and panel == "503") else ""))
            full[(panel, k)] = tests
        # discovery sign: marginal Wald ratio at the discovery instrument rs73045269 (both resources)
        if panel == "503":
            for k in OUTC:
                disc[k] = marg[k][0][0] / marg["protein"][0][0]
                b0, s0 = marg[k][0], marg[k][1]
                bp, sp = marg["protein"][0], marg["protein"][1]
                for x, inst in ((0, A_RS), (1, B_RS)):
                    wr = b0[x] / bp[x]
                    wse = math.sqrt(s0[x] ** 2 / bp[x] ** 2 + b0[x] ** 2 * sp[x] ** 2 / bp[x] ** 4)
                    add("S2.marg.ratio.%s.%s.value" % (k, inst), "marginal Wald ratio at %s (descriptive; not conditioned)" % inst, k, panel,
                        wr, se=wse, p=p2(wr / wse), se_order="2nd", unit="log-odds per SD protein",
                        notes="marginal z protein %.3f, outcome %.3f" % (bp[x] / sp[x], b0[x] / s0[x]))
    # ---------------- discovery sign cross-check with the pipeline MR table (raw_input)
    import csv
    mr = [r for r in csv.DictReader(open(W.res("TAB_mr_all_results"), encoding="utf-8")) if r["gene"] == "TGFB1" and r["outcome_id"] == "FinnGen_R10_I9_CHD"]
    assert len(mr) == 1
    dsign = 1 if float(mr[0]["b"]) > 0 else -1
    W.log("discovery TGFB1-CAD estimate (TAB_mr_all_results, FinnGen_R10_I9_CHD, Wald ratio): b = %s -> discovery sign %+d; "
          "recomputed marginal ratio at rs73045269: FinnGen %.6f, ieu-a-7 %.6f" % (mr[0]["b"], dsign, disc["finngen_chd"], disc["ieu_a_7"]))
    if (disc["finngen_chd"] > 0) != (dsign > 0):
        W.log("ASSERT_FAILED discovery sign differs between TAB_mr_all_results and the recomputed marginal ratio")
    add("S2.discovery_sign", "sign of the discovery TGFB1-CAD estimate (Wald ratio at rs73045269)", "TAB_mr_all_results (FinnGen_R10_I9_CHD)", "",
        dsign, notes="pipeline b = %s; recomputed marginal Wald ratio FinnGen %.6f" % (mr[0]["b"], disc["finngen_chd"]))
    # ---------------- bootstrap (5,000 draws of the 503 people)
    bed = W.res("LD503_bed")
    G = bed_rows(bed, 503, [int(rows[iA]["bed_idx503"]), int(rows[iB]["bed_idx503"])])
    r_full_np = float(np.corrcoef(G, rowvar=False)[0, 1])
    W.log("bed read in place (no copy): r(A,B) full 503 from the numpy reader = %.6f vs PLINK %.6f (diff %.2e)" % (r_full_np, r503, r_full_np - r503))
    if abs(r_full_np - r503) > 1e-5:
        W.log("ASSERT_FAILED numpy-reader r differs from PLINK r")
    rng = np.random.default_rng(W.SEED)
    draws = []
    nmono = 0
    IDX = np.empty((NBOOT, 503), dtype="<i2")      # person indices of every draw (0-based), reused by the R coloc bootstrap
    for d in range(NBOOT):
        ix = rng.integers(0, 503, 503)
        IDX[d] = ix
        g = G[ix]
        if g[:, 0].std() == 0 or g[:, 1].std() == 0:
            rb = 0.0
            nmono += 1
        else:
            rb = float(np.corrcoef(g, rowvar=False)[0, 1])
        rec = {"draw": d + 1, "r": rb}
        JP = joint_z(*marg["protein"][:2], rb)
        for k in OUTC:
            JY = joint_z(*marg[k][:2], rb)
            t = ratio_test(JY[0], JY[1], JP[0], JP[1], True, True)
            rec["d_" + k], rec["p_" + k] = t["d"], t["p"]
            rec["rhoA_" + k], rec["rhoB_" + k] = t["rhoA"], t["rhoB"]
        draws.append(rec)
        if (d + 1) % 1000 == 0:
            W.log("[boot] %d / %d draws" % (d + 1, NBOOT))
    W.write_tsv(W.OUT + "/boot_draws.tsv", list(draws[0].keys()), draws)
    IDX.tofile(W.OUT + "/boot_idx_5000x503.i2")
    W.log("bootstrap person indices written: boot_idx_5000x503.i2 (int16 LE, 0-based; md5 %s)" % W.md5_file(W.OUT + "/boot_idx_5000x503.i2"))
    rv = np.array([x["r"] for x in draws])
    lo, hi = np.percentile(rv, [2.5, 97.5])
    add("S2.r.boot.lo", "2.5th percentile of r(rs73045269:T, rs12983775:G) over 5,000 bootstrap draws of the 503 people", "LD503", "503",
        float(lo), n=NBOOT, notes="seed 20261007, numpy default_rng (PCG64); DEVIATION: bootstrap outside PLINK (section 10); %d monomorphic draws" % nmono)
    add("S2.r.boot.hi", "97.5th percentile of r(A,B) over 5,000 bootstrap draws", "LD503", "503", float(hi), n=NBOOT,
        notes="seed 20261007, numpy default_rng (PCG64); DEVIATION: bootstrap outside PLINK (section 10)")
    boot = {}
    for k in OUTC:
        t = full[("503", k)]["p_aware_2nd"]
        concl = (t["p"] < 0.05, np.sign(t["d"]))
        same = np.array([(x["p_" + k] < 0.05, np.sign(x["d_" + k])) == concl for x in draws])
        share = float(same.mean())
        pv = np.array([x["p_" + k] for x in draws])
        boot[k] = share
        add("S2.het.boot.%s.share" % k, "share of 5,000 503-panel bootstrap draws with the full-panel conclusion (side of p 0.05 and sign of d) "
            "of the aware second-order A-vs-B test", k, "503", share, n=NBOOT, cov="aware", se_order="2nd",
            notes="full-panel conclusion: p %s 0.05, d %s 0; bootstrap p 2.5-97.5%%: %.3g-%.3g; MC SE of share %.4f%s" % (
                "<" if concl[0] else ">=", "+" if concl[1] > 0 else "-", np.percentile(pv, 2.5), np.percentile(pv, 97.5),
                math.sqrt(share * (1 - share) / NBOOT), "; NAMED PRIMARY CELL for Rule 3(ii)" if k == "finngen_chd" else ""))
        add("S2.het.boot.%s.p_lo" % k, "2.5th percentile of the aware 2nd-order p over the draws", k, "503", float(np.percentile(pv, 2.5)), n=NBOOT)
        add("S2.het.boot.%s.p_hi" % k, "97.5th percentile of the aware 2nd-order p over the draws", k, "503", float(np.percentile(pv, 97.5)), n=NBOOT)
    W.log("bootstrap: r 95%% interval %.4f to %.4f; share same conclusion FinnGen %.4f, ieu-a-7 %.4f" % (lo, hi, boot["finngen_chd"], boot["ieu_a_7"]))
    # ---------------- conditioned window data for the Rule 3(iv) colocalization (R coloc)
    cw = []
    for panel, idx, R in (("503", [i for i, r in enumerate(rows) if r["in_p1np"] == "1"], R503), ("489", p2rows, R489)):
        if panel == "503":
            Rp = R[np.ix_(idx, idx)]
        else:
            Rp = R
        loc = {i: j for j, i in enumerate(idx)}
        for sig, cond_on in (("A", iB), ("B", iA)):
            c = loc[cond_on]
            rc = Rp[:, c]
            q = rc ** 2
            for res_name, t in TR.items():
                z = np.array([float(rows[i]["z__" + t]) for i in idx])
                se = np.array([float(rows[i]["se__" + t]) for i in idx])
                zc = (z - rc * z[c]) / np.sqrt(np.clip(1 - q, 1e-12, None))
                sc = se / np.sqrt(np.clip(1 - q, 1e-12, None))
                for j, i in enumerate(idx):
                    if q[j] > COLLIN or j == c:
                        continue
                    cw.append({"panel": panel, "signal": sig, "conditioned_on": rows[cond_on]["rsid503"], "trait": res_name,
                               "snp": rows[i]["snp"], "rsid": rows[i]["rsid503"], "b": zc[j] * sc[j], "se": sc[j], "z": zc[j]})
    W.write_tsv(W.OUT + "/cond_window.tsv", ["panel", "signal", "conditioned_on", "trait", "snp", "rsid", "b", "se", "z"], cw)
    W.log("conditioned window data written: %d rows (503 = p1np, 489 = p2np; r^2 > %.1f to the conditioning SNP dropped)" % (len(cw), COLLIN))
    W.write_tsv(W.OUT + "/results_cond.tsv", list(out[0].keys()), out)
    for o in out:
        if ".het." in o["quantity_id"] and "boot" not in o["quantity_id"]:
            W.log("  %-40s panel %s d %+.4f se %.4f p %.3g" % (o["quantity_id"], o["panel"], o["value"], o["se"], o["p"]))
        if ".ratio." in o["quantity_id"] and o["quantity_id"].endswith("value"):
            W.log("  %-50s panel %s %+.4f (se %.4f)" % (o["quantity_id"], o["panel"], o["value"], o["se"]))


if __name__ == "__main__":
    main()
