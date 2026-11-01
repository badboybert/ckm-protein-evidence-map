#!/usr/bin/env python
"""l3fm__conditional.py -- summary-statistic CONDITIONAL colocalization (COJO-style z-approximation with the SAME external
1000G EUR LD as the SuSiE step) at each fine-mapping locus, plus the D1 H4-lead re-derivation. OFFLINE; reads only;
writes staging/l3-finemap/<locus>/conditional/ (or --outroot for mutation runs).

  python -B V3/V4_local/builders/l3fm__conditional.py [--locus L ...] [--outroot DIR] [--mutation none|ld_shuffle]

Conditioning (Yang et al. 2012 z-approximation; equal-N, standardised genotypes):
  z_i|S = (z_i - r_iS R_SS^-1 z_S) / sqrt(1 - r_iS R_SS^-1 r_Si),  se_i|S = se_i / sqrt(1 - r_iS R_SS^-1 r_Si),
  SNPs with multiple r^2 to S > 0.9 are dropped from the conditioned dataset (collinearity), never inflated.
Signals: forward stepwise selection (lead first; add the top conditional SNP while p_cond < 5e-8, candidates with multiple
  r^2 > 0.9 to the selected set excluded; at most 5; out-of-sample-LD guard: a secondary must have marginal p < 1e-5 and
  the same sign conditionally as marginally -- see select_signals).
Per protein x metabolite pair (coloc.abf = the D1 vendored engine, p1 = p2 = 1e-4, p12 = 1e-5):
  base            unconditioned, on the fine-mapping SNP set
  P|leadP         protein conditioned on its own lead -> is there a SECOND protein signal that colocalizes?
  M|leadP         metabolite conditioned on the protein lead -> does the metabolite signal vanish (max |z| before/after)?
  M|allP          metabolite conditioned on ALL the protein's selected signals (max |z| after) -> is it fully explained?
  P|leadM         protein conditioned on the metabolite lead (max |z| before/after)
  signals         coloc.signals-style: each protein signal x each metabolite signal, each conditioned on its OTHER signals
  M|leadQ         for every other protein Q at the locus with a genome-wide lead: metabolite conditioned on Q's lead, then
                  coloc with P (does P's colocalization survive removing Q's signal?) and vice versa.

TWO-SIGNAL MODULE (--variant p1np | p2np; added 2026-09-26, promoting the independent verifier's decisive methods; the
rules are POST HOC -- adopted after the verifier's results at DPEP1/ARG1 had been seen):
  SNP sets   p1np = harmonised SNPs minus palindromic ones, panel-1 LD (503); p2np = the p1np SNPs present in panel 2
             (1000G.EUR.QC, 489 -- all among panel 1's 503), panel-2 LD. Built by l3fm__panel2.py.
  signals    A = protein lead, B = the protein's second forward-selected signal (select_signals, same GW threshold and
             out-of-sample-LD guard; max 2), selected ONCE on p1np with panel-1 LD; p2np conditions on the SAME SNPs
             (a SNP absent from p2np -> that pair's p2np row is NA). Palindromic SNPs can therefore never be A or B.
  cond-both  both traits conditioned on the SAME set (the verifier's 'coloc DPEP1 x M|B-hap'): H4_A = coloc.abf(P|B, M|B)
             (only signal A left in both), H4_B = coloc.abf(P|A, M|A); cond_H4 = max(H4_A, H4_B). A single-signal
             protein (no GW B) gets the unconditioned coloc.abf on the set (cond_signal = 'single').
  ratio      the verifier's two-signal proportionality rule: with PA = z_P(A|B), PB = z_P(B|A), MA = z_M(A|B),
             MB = z_M(B|A), k_A = MA/PA, k_B = MB/PB, dev = MA - k_B * PA; proportional iff |dev| <= DEV_MAX = 3.5 z
             (the verifier's declared conditional-z error budget for reference LD at |z| ~ 100). dev_sd = the dev
             standardised by its first-order sd, sqrt(1 + k_B^2 + PA^2 (1 + k_B^2) / PB^2) (unit-variance conditional z),
             reported only. Applicable only when the protein has a GW second signal.
  reproduction  rows with signal_rule = 'verifier_handpicked' recompute the module at the verifier's hand-picked A/B
             (DPEP1 rs409170/rs164749, ARG1 rs2781668/rs17788484) so its numbers can be compared; NOT used for classes.
"""
import argparse
import csv
import gzip
import json
import math
import os
import sys
import time

import numpy as np
from scipy.stats import norm

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import l3fm__common as C  # noqa: E402

B = C.B
E = B.load_engine()
GW = 5e-8
Z_GW = float(norm.isf(GW / 2))
COLLIN = 0.9
Z_MARG = float(norm.isf(1e-5 / 2))
DEV_MAX = 3.5
VERIFIER_AB = {"DPEP1_16q24": ("P__DPEP1", "rs409170", "rs164749"), "ARG1_6q23": ("P__ARG1", "rs2781668", "rs17788484")}
VARIANTS = {"primary": (None, "p1"), "p1np": ("p1np", "p1"), "p2np": ("p2np", "p2")}


def load_locus(ldir, gfile, mutation="none", seed=20260926, snpset=None, nref=503):
    with gzip.open(os.path.join(ldir, "harmonised.tsv.gz"), "rt", encoding="utf-8") as fh:
        h = list(csv.DictReader(fh, delimiter="\t"))
    tr = list(csv.DictReader(open(os.path.join(ldir, "traits.tsv"), encoding="utf-8"), delimiter="\t"))
    p_full = len(h)
    G = np.fromfile(gfile, dtype=np.int8)
    if G.size % nref:
        C.refuse("%s: genotype extract has %d values, not a multiple of %d" % (ldir, G.size, nref))
    G = G.reshape(nref, G.size // nref).astype(np.float64)
    af_exp = np.array([float(r["af_ref"]) for r in h])
    if snpset:
        ss = list(csv.DictReader(open(snpset, encoding="utf-8"), delimiter="\t"))
        rows = [int(r["row"]) for r in ss]
        if any(h[j]["snp"] != r["snp"] for j, r in zip(rows, ss)):
            C.refuse("%s: snpset ids differ from the harmonised rows" % snpset)
        if G.shape[1] == p_full and len(rows) != p_full:
            G = G[:, rows]
        h = [h[j] for j in rows]
        af_exp = np.array([float(r["af"]) for r in ss])
    p = len(h)
    if G.shape[1] != p:
        C.refuse("%s: genotype extract has %d SNP columns, expected %d" % (ldir, G.shape[1], p))
    af = G.mean(axis=0) / 2
    if np.max(np.abs(af - af_exp)) > 1e-5:
        C.refuse("%s: genotype extract != expected A1 frequency" % ldir)
    if mutation == "ld_shuffle":
        G = G[:, np.random.default_rng(seed).permutation(p)]
    R = np.corrcoef(G, rowvar=False)
    data = {}
    for t in tr:
        tid = t["tid"]
        b = np.array([float(r["b__" + tid]) for r in h]); se = np.array([float(r["se__" + tid]) for r in h])
        data[tid] = {"b": b, "se": se, "z": b / se, "kind": t["kind"], "label": t["label"], "n": float(t["n"])}
    snps = [r["snp"] for r in h]
    rs = [r["rsid_ref"] for r in h]
    pos = np.array([int(r["pos37"]) for r in h])
    return h, tr, R, data, snps, rs, pos


def condition(z, se, R, S):
    """conditional z and se given the selected index set S; NaN where collinear (multiple r^2 > COLLIN) or in S."""
    S = list(S)
    if not S:
        return z.copy(), se.copy()
    RSS = R[np.ix_(S, S)]
    RiS = R[:, S]
    W = np.linalg.solve(RSS, RiS.T)                     # |S| x p
    q = np.einsum("ij,ji->i", RiS, W)                   # r_iS R_SS^-1 r_Si
    zc = (z - RiS @ np.linalg.solve(RSS, z[S])) / np.sqrt(np.clip(1 - q, 1e-12, None))
    sec = se / np.sqrt(np.clip(1 - q, 1e-12, None))
    bad = q > COLLIN
    bad[S] = True
    zc[bad] = np.nan
    sec[bad] = np.nan
    return zc, sec


def select_signals(z, se, R, max_k=5):
    """forward selection. OUT-OF-SAMPLE-LD GUARD (declared deviation from plain COJO): a secondary candidate must also be
    marginally associated (|z_marginal| >= Z_MARG, p < 1e-5) and keep the SIGN of its marginal z; a candidate whose
    conditional z is created or sign-flipped by conditioning is an LD-mismatch artefact at strong leads (seen at
    DPEP2/UMOD/GSTA3 without the guard) and is skipped, not selected. Returns (selected, z_at_selection, n_rejected)."""
    S = [int(np.nanargmax(np.abs(z)))]
    if abs(z[S[0]]) < Z_GW:
        return S, [float(z[S[0]])], 0
    zs = [float(z[S[0]])]
    rejected = 0
    while len(S) < max_k:
        zc, _ = condition(z, se, R, S)
        ok = (~np.isnan(zc)) & (np.abs(zc) >= Z_GW)
        if not ok.any():
            break
        guard = ok & (np.abs(z) >= Z_MARG) & (np.sign(zc) == np.sign(z))
        if not guard.any():
            rejected += int(ok.sum())
            break
        j = int(np.nanargmax(np.where(guard, np.abs(zc), -1)))
        rejected += int((ok & ~guard & (np.abs(zc) > abs(zc[j]))).sum())     # candidates the guard passed over
        S.append(j); zs.append(float(zc[j]))
    return S, zs, rejected


def abf(b1, se1, b2, se2):
    ok = ~(np.isnan(b1) | np.isnan(b2) | np.isnan(se1) | np.isnan(se2))
    if ok.sum() < 20:
        return None, int(ok.sum())
    return E.coloc_abf(b1[ok], se1[ok], b2[ok], se2[ok], p1=1e-4, p2=1e-4, p12=1e-5), int(ok.sum())


def d1_h4_lead(protein, metabolite):
    """re-derive the D1 H4-lead (argmax SNP.PP.H4 = argmax of lABF1 + lABF2) on the D1 SNP set, and its protein |z|."""
    units = [u for u in B.derive_units() if u["protein"] == protein and u["metabolite"] == metabolite]
    if len(units) != 1:
        C.refuse("D1 unit %s/%s not unique" % (protein, metabolite))
    u = units[0]
    prot, _, _ = B.load_eur_protein(B.eur_window_path(u["protein"], u["oid"], u["chr"]), u["chr"], u["lo"], u["hi"])
    met, _, _, _ = B.load_metab_window(B.metab_window_path(u))
    keys, b1, s1, b2, s2 = [], [], [], [], []
    for k, pp in prot.items():
        m = met.get(k)
        if not m:
            continue
        if (pp["ea"], pp["oa"]) == (m["ea"], m["oa"]): mb = m["b"]
        elif (pp["ea"], pp["oa"]) == (m["oa"], m["ea"]): mb = -m["b"]
        elif (pp["ea"], pp["oa"]) == (B.COMP.get(m["ea"]), B.COMP.get(m["oa"])): mb = m["b"]
        elif (pp["ea"], pp["oa"]) == (B.COMP.get(m["oa"]), B.COMP.get(m["ea"])): mb = -m["b"]
        else: continue
        keys.append(k); b1.append(pp["b"]); s1.append(pp["se"]); b2.append(mb); s2.append(m["se"])
    b1, s1, b2, s2 = map(np.array, (b1, s1, b2, s2))
    l1 = E._labf(b1, s1, 0.15 ** 2); l2 = E._labf(b2, s2, 0.15 ** 2)
    j = int(np.argmax(l1 + l2))
    pp = E.coloc_abf(b1, s1, b2, s2, p1=1e-4, p2=1e-4, p12=1e-5)
    w = np.exp((l1 + l2) - np.max(l1 + l2)); w /= w.sum()
    return {"nsnp": len(keys), "H4": pp["H4"], "lead_pos38": keys[j][1], "lead_snp_pp_h4": float(w[j]),
            "z_protein_at_lead": float(b1[j] / s1[j]), "z_metab_at_lead": float(b2[j] / s2[j])}


def run_locus(name, outroot=None, mutation="none", log=print):
    ldir = os.path.join(C.STAGE_FM, name)
    gfile = os.path.join(C.RAW_FM, "%s__ref_A1dosage_503xP.i8" % name)
    h, tr, R, D, snps, rs, pos = load_locus(ldir, gfile, mutation)
    odir = outroot or os.path.join(ldir, "conditional")
    os.makedirs(odir, exist_ok=True)
    p38 = {int(r["pos38"]) if r["pos38"] else -1: i for i, r in enumerate(h)}
    prots = [t["tid"] for t in tr if t["kind"] == "protein"]
    others = [t["tid"] for t in tr if t["kind"] != "protein"]
    sig = {}
    srows = []
    for tid, d in D.items():
        S, zs, nrej = select_signals(d["z"], d["se"], R)
        sig[tid] = S
        for k, (j, zz) in enumerate(zip(S, zs)):
            srows.append({"tid": tid, "label": d["label"], "signal": k + 1, "snp": snps[j], "rsid": rs[j], "pos37": int(pos[j]),
                          "z_marginal": "%.3f" % d["z"][j], "z_conditional_at_selection": "%.3f" % zz,
                          "gw_significant": abs(zz) >= Z_GW, "guard_rejected_stronger_candidates": nrej})
    rows = []
    for pt in prots:
        P = D[pt]
        lp = sig[pt][0]
        for mt in others:
            M = D[mt]
            lm = sig[mt][0]
            r = {"protein": pt, "other": mt, "nsnp": len(snps), "lead_protein": rs[lp], "lead_other": rs[lm],
                 "r2_leads": "%.4f" % R[lp, lm] ** 2, "z_protein_lead": "%.2f" % P["z"][lp], "z_other_lead": "%.2f" % M["z"][lm],
                 "z_protein_at_other_lead": "%.2f" % P["z"][lm], "z_other_at_protein_lead": "%.2f" % M["z"][lp],
                 "protein_gw": abs(P["z"][lp]) >= Z_GW, "n_signals_protein": len(sig[pt]), "n_signals_other": len(sig[mt])}
            base, n0 = abf(P["b"], P["se"], M["b"], M["se"])
            r["abf_H4_base"], r["abf_H3_base"] = base["H4"], base["H3"]
            zc, sc = condition(P["z"], P["se"], R, [lp])
            res, n1 = abf(zc * sc, sc, M["b"], M["se"])
            r["abf_H4_P_given_leadP"] = res["H4"] if res else ""
            r["maxabsz_P_given_leadP"] = "%.2f" % np.nanmax(np.abs(zc))
            zc, sc = condition(M["z"], M["se"], R, [lp])
            res, _ = abf(P["b"], P["se"], zc * sc, sc)
            r["abf_H4_M_given_leadP"] = res["H4"] if res else ""
            r["maxabsz_M_before"] = "%.2f" % np.max(np.abs(M["z"]))
            r["maxabsz_M_given_leadP"] = "%.2f" % np.nanmax(np.abs(zc))
            zc_all, _ = condition(M["z"], M["se"], R, sig[pt])
            r["maxabsz_M_given_allP"] = "%.2f" % np.nanmax(np.abs(zc_all))
            zc, sc = condition(P["z"], P["se"], R, [lm])
            r["maxabsz_P_before"] = "%.2f" % np.max(np.abs(P["z"]))
            r["maxabsz_P_given_leadM"] = "%.2f" % np.nanmax(np.abs(zc))
            # signal-wise
            best = (None, -1)
            cells = []
            for kp, jp in enumerate(sig[pt]):
                zp, sp = condition(P["z"], P["se"], R, [x for x in sig[pt] if x != jp])
                for km, jm in enumerate(sig[mt]):
                    zm, sm = condition(M["z"], M["se"], R, [x for x in sig[mt] if x != jm])
                    res, _ = abf(zp * sp, sp, zm * sm, sm)
                    h4 = res["H4"] if res else float("nan")
                    cells.append("P%d(%s)xM%d(%s)=%.3f" % (kp + 1, rs[jp], km + 1, rs[jm], h4))
                    if h4 == h4 and h4 > best[1]:
                        best = ((kp + 1, rs[jp], km + 1, rs[jm]), h4)
            r["signals_max_H4"] = "%.4f" % best[1] if best[0] else ""
            r["signals_argmax"] = "P%d(%s) x M%d(%s)" % best[0] if best[0] else ""
            r["signals_grid"] = "; ".join(cells)
            # other proteins' leads
            qtxt = []
            for qt in prots:
                if qt == pt:
                    continue
                Q = D[qt]
                lq = sig[qt][0]
                if abs(Q["z"][lq]) < Z_GW:
                    continue
                zm, sm = condition(M["z"], M["se"], R, [lq])
                res_pq, _ = abf(P["b"], P["se"], zm * sm, sm)
                zm2, sm2 = condition(M["z"], M["se"], R, [lp])
                res_qp, _ = abf(Q["b"], Q["se"], zm2 * sm2, sm2)
                res_q, _ = abf(Q["b"], Q["se"], M["b"], M["se"])
                qtxt.append("%s: H4(Q,M)=%.3f; H4(P,M|leadQ %s)=%s; H4(Q,M|leadP)=%s; r2(leadP,leadQ)=%.3f; max|zM|leadQ|=%.2f" % (
                    qt, res_q["H4"], rs[lq], "%.3f" % res_pq["H4"] if res_pq else "NA", "%.3f" % res_qp["H4"] if res_qp else "NA",
                    R[lp, lq] ** 2, np.nanmax(np.abs(zm))))
            r["other_protein_conditioning"] = " || ".join(qtxt)
            rows.append(r)
            log("  %-10s x %-16s base=%.3f P|leadP=%s M|leadP=%s  max|zM| %s -> %s ; signals %s %s" % (
                pt, mt, base["H4"], _f(r["abf_H4_P_given_leadP"]), _f(r["abf_H4_M_given_leadP"]), r["maxabsz_M_before"],
                r["maxabsz_M_given_leadP"], r["signals_max_H4"], r["signals_argmax"]))
    C.write_tsv(os.path.join(odir, "conditional_pairs.tsv"), list(rows[0].keys()), rows)
    C.write_tsv(os.path.join(odir, "signals.tsv"), list(srows[0].keys()), srows)
    # D1 H4-lead re-derivation (bridge protein only; skip for the UMOD control, which has no D1 unit)
    d1rows = []
    if mutation == "none" and not name.startswith("UMOD"):
        bp = [t for t in tr if t["role"] == "bridge_protein"][0]
        for t in tr:
            if t["kind"] != "metabolite":
                continue
            d = d1_h4_lead(bp["gene"], t["label"])
            i = p38.get(d["lead_pos38"])
            d1rows.append({"protein": bp["gene"], "metabolite": t["label"], "d1_nsnp": d["nsnp"], "d1_H4": "%.4f" % d["H4"],
                           "d1_H4_lead_pos38": d["lead_pos38"], "d1_H4_lead_snp_pp": "%.3f" % d["lead_snp_pp_h4"],
                           "z_protein_at_d1_lead": "%.2f" % d["z_protein_at_lead"], "z_metab_at_d1_lead": "%.2f" % d["z_metab_at_lead"],
                           "d1_lead_in_finemap_set": i is not None, "d1_lead_rsid": rs[i] if i is not None else "",
                           "d1_lead_pos37": int(pos[i]) if i is not None else ""})
        C.write_tsv(os.path.join(odir, "d1_h4_leads.tsv"), list(d1rows[0].keys()), d1rows)
    C.write_json(os.path.join(odir, "COND_STAMP.json"), {
        "locus": name, "mutation": mutation, "nsnp": len(snps), "z_gw": Z_GW, "collinearity_r2": COLLIN,
        "builder_md5": C.md5_file(os.path.abspath(__file__)), "common_md5": C.md5_file(os.path.join(HERE, "l3fm__common.py")),
        "harmonised_md5": C.md5_file(os.path.join(ldir, "harmonised.tsv.gz")), "geno_md5": C.md5_file(gfile),
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S")})
    return rows, srows, d1rows, R


def _f(x):
    return "%.3f" % x if isinstance(x, float) else str(x)


def ratio_rule(MA, MB, PA, PB):
    """the verifier's two-signal proportionality rule -> (proportional, dev, k_A, k_B, dev_sd)."""
    kA, kB = MA / PA, MB / PB
    dev = MA - kB * PA
    sd = math.sqrt(1 + kB * kB + PA * PA * (1 + kB * kB) / (PB * PB))
    return abs(dev) <= DEV_MAX, dev, kA, kB, dev / sd


def abf_lead(bp, sp, bm, sm, snps):
    """coloc.abf on the non-NaN SNPs plus the SNP with the largest SNP.PP.H4 (lABF1 + lABF2)."""
    ok = ~(np.isnan(bp) | np.isnan(bm) | np.isnan(sp) | np.isnan(sm))
    if ok.sum() < 20:
        return None, "", int(ok.sum())
    pp = E.coloc_abf(bp[ok], sp[ok], bm[ok], sm[ok], p1=1e-4, p2=1e-4, p12=1e-5)
    l = E._labf(bp[ok], sp[ok], 0.15 ** 2) + E._labf(bm[ok], sm[ok], 0.15 ** 2)
    return pp, np.array(snps)[ok][int(np.argmax(l))], int(ok.sum())


def two_signal_rows(name, variant, D, R, snps, rs, pos, prots, others, sigAB, rule):
    """one row per protein x other trait: cond-both coloc + two-signal ratio (see the module docstring)."""
    ix = {x: i for i, x in enumerate(rs)}
    out = []
    for pt in prots:
        P = D[pt]
        ab = sigAB.get(pt)
        for mt in others:
            M = D[mt]
            r = {"locus": name, "variant": variant, "signal_rule": rule, "protein": pt, "other": mt, "nsnp": len(snps),
                 "sig_A": "", "sig_B": "", "cond_signal": "", "cond_H4": "", "cond_H3": "", "cond_H4_A": "", "cond_H4_B": "",
                 "cond_lead_A": "", "cond_lead_B": "", "PA": "", "PB": "", "MA": "", "MB": "", "k_A": "", "k_B": "",
                 "ratio_dev": "", "ratio_dev_sd": "", "ratio_proportional": "", "maxabsz_M_given_AB": "", "note": ""}
            if ab is None:
                r["note"] = "protein signals not selected on p1np"
                out.append(r)
                continue
            a_rs, b_rs = ab
            if a_rs not in ix or (b_rs and b_rs not in ix):
                r["note"] = "signal SNP absent from this SNP set (%s/%s)" % (a_rs, b_rs)
                out.append(r)
                continue
            r["sig_A"], r["sig_B"] = a_rs, b_rs
            a = ix[a_rs]
            if not b_rs:
                pp, lead, n = abf_lead(P["b"], P["se"], M["b"], M["se"], snps)
                r.update({"cond_signal": "single", "cond_H4": pp["H4"], "cond_H3": pp["H3"], "cond_lead_A": lead,
                          "note": "protein has one GW signal: unconditioned coloc.abf; ratio test not applicable"})
                zc, _ = condition(M["z"], M["se"], R, [a])
                r["maxabsz_M_given_AB"] = "%.2f" % np.nanmax(np.abs(zc))
                out.append(r)
                continue
            b = ix[b_rs]
            zpB, spB = condition(P["z"], P["se"], R, [b])      # only A left
            zmB, smB = condition(M["z"], M["se"], R, [b])
            zpA, spA = condition(P["z"], P["se"], R, [a])      # only B left
            zmA, smA = condition(M["z"], M["se"], R, [a])
            ppA, leadA, _ = abf_lead(zpB * spB, spB, zmB * smB, smB, snps)
            ppB, leadB, _ = abf_lead(zpA * spA, spA, zmA * smA, smA, snps)
            h4a = ppA["H4"] if ppA else float("nan")
            h4b = ppB["H4"] if ppB else float("nan")
            best = "A" if (h4a >= h4b or h4b != h4b) else "B"
            PA, PB, MA, MB = float(zpB[a]), float(zpA[b]), float(zmB[a]), float(zmA[b])
            prop, dev, kA, kB, dsd = ratio_rule(MA, MB, PA, PB)
            zab, _ = condition(M["z"], M["se"], R, [a, b])
            r.update({"cond_signal": best, "cond_H4": max(h4a, h4b), "cond_H3": (ppA if best == "A" else ppB)["H3"],
                      "cond_H4_A": h4a, "cond_H4_B": h4b, "cond_lead_A": leadA, "cond_lead_B": leadB,
                      "PA": "%.3f" % PA, "PB": "%.3f" % PB, "MA": "%.3f" % MA, "MB": "%.3f" % MB, "k_A": "%.4f" % kA, "k_B": "%.4f" % kB,
                      "ratio_dev": "%.3f" % dev, "ratio_dev_sd": "%.3f" % dsd, "ratio_proportional": prop,
                      "maxabsz_M_given_AB": "%.2f" % np.nanmax(np.abs(zab))})
            out.append(r)
    return out


def run_two_signal(name, variant, outroot=None, log=print):
    """the two-signal module on SNP set p1np or p2np (signals selected once, on p1np)."""
    sset, panel = VARIANTS[variant]
    ldir = os.path.join(C.STAGE_FM, name)
    sfile = os.path.join(ldir, "snpsets", sset + ".tsv")
    gfile = os.path.join(C.RAW_FM, ("%s__ref_A1dosage_503xP.i8" if panel == "p1" else "%s__ref2_A1dosage_489xP.i8") % name)
    h, tr, R, D, snps, rs, pos = load_locus(ldir, gfile, snpset=sfile, nref=503 if panel == "p1" else C.REF2_N)
    odir = outroot or os.path.join(ldir, "conditional_" + variant)
    os.makedirs(odir, exist_ok=True)
    prots = [t["tid"] for t in tr if t["kind"] == "protein"]
    others = [t["tid"] for t in tr if t["kind"] != "protein"]
    sigfile = os.path.join(ldir, "conditional_p1np", "signals_AB.tsv")
    if variant == "p1np":
        srows = []
        for tid in prots:
            S, zs, nrej = select_signals(D[tid]["z"], D[tid]["se"], R, max_k=2)
            gwB = len(S) > 1 and abs(zs[1]) >= Z_GW
            srows.append({"tid": tid, "sig_A": rs[S[0]], "z_A": "%.3f" % D[tid]["z"][S[0]], "A_gw": abs(zs[0]) >= Z_GW,
                          "sig_B": rs[S[1]] if gwB else "", "z_B_given_A_at_selection": "%.3f" % zs[1] if gwB else "",
                          "guard_rejected_stronger_candidates": nrej})
        C.write_tsv(os.path.join(odir, "signals_AB.tsv"), list(srows[0].keys()), srows)
        sigfile = os.path.join(odir, "signals_AB.tsv")
    if not os.path.exists(sigfile):
        C.refuse("%s: %s missing (run --variant p1np first)" % (name, sigfile))
    sigAB = {r["tid"]: (r["sig_A"], r["sig_B"]) for r in csv.DictReader(open(sigfile, encoding="utf-8"), delimiter="\t")}
    rows = two_signal_rows(name, variant, D, R, snps, rs, pos, prots, others, sigAB, "declared_AB")
    if name in VERIFIER_AB:
        pt, a_rs, b_rs = VERIFIER_AB[name]
        rows += two_signal_rows(name, variant, D, R, snps, rs, pos, [pt], others, {pt: (a_rs, b_rs)}, "verifier_handpicked")
    C.write_tsv(os.path.join(odir, "two_signal.tsv"), list(rows[0].keys()), rows)
    C.write_json(os.path.join(odir, "TWO_SIGNAL_STAMP.json"), {
        "locus": name, "variant": variant, "nsnp": len(snps), "dev_max": DEV_MAX, "z_gw": Z_GW,
        "snpset_md5": C.md5_file(sfile), "signals_md5": C.md5_file(sigfile), "geno_md5": C.md5_file(gfile),
        "harmonised_md5": C.md5_file(os.path.join(ldir, "harmonised.tsv.gz")),
        "builder_md5": C.md5_file(os.path.abspath(__file__)), "common_md5": C.md5_file(os.path.join(HERE, "l3fm__common.py")),
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S")})
    for r in rows:
        if r["cond_H4"] != "":
            log("  [%s %s] %-10s x %-16s A %s B %s cond %s H4 %.3f dev %s (%s)" % (
                variant, r["signal_rule"][:8], r["protein"], r["other"], r["sig_A"], r["sig_B"] or "-", r["cond_signal"],
                float(r["cond_H4"]), r["ratio_dev"] or "n/a", r["ratio_proportional"]))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", nargs="*")
    ap.add_argument("--outroot")
    ap.add_argument("--mutation", default="none", choices=["none", "ld_shuffle"])
    ap.add_argument("--variant", default="primary", choices=list(VARIANTS))
    a = ap.parse_args(argv)
    names = a.locus or list(C.LOCI) + ["UMOD_eGFR_control"]
    for n in names:
        print("[cond %s] %s" % (a.variant, n), flush=True)
        if a.variant == "primary":
            run_locus(n, os.path.join(a.outroot, n) if a.outroot else None, a.mutation, lambda m: print(m, flush=True))
        else:
            run_two_signal(n, a.variant, os.path.join(a.outroot, n) if a.outroot else None, lambda m: print(m, flush=True))


if __name__ == "__main__":
    main()
