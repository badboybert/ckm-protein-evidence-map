"""w2_02_locus.py -- the TGFB1 locus table, the four l3fm-style variant sets and the PLINK LD matrices (section 2).

Window   = the UKB-PPP cis window clipped to +/-500 kb of the gene = instrument-table cis_start..cis_end (GRCh38;
           cis = gene_start - 500 kb .. gene_end + 500 kb in the pipeline's own gene coordinates).
Traits   = P (UKB-PPP OID20621, effect = ALLELE1), CHD_FG (FinnGen R10 I9_CHD, effect = alt, GRCh38),
           CAD_C4D (CARDIoGRAMplusC4D ieu-a-7, effect = ALT, GRCh37).
Joins    = P x FinnGen on GRCh38 pos + compatible alleles; P x ieu-a-7 and P x reference on the GRCh37 position embedded
           in the UKB-PPP ID + compatible alleles (complement only for non-palindromic SNPs; a non-unique match is dropped).
Engine QC (DEVIATION: engine path-locked -- re-implemented from l3fm__prep.py / l3fm__panel2.py / l3fm__robustness.py):
  susie = joined SNPs with ref (503) MAF >= 0.01, not monomorphic, palindromic only at ref MAF <= 0.40, and
          |aligned GWAS A1 frequency - ref A1 frequency| <= 0.20 for every trait reporting a frequency;
  p1np  = susie minus palindromes;
  p1fr  = p1np + susie palindromes with every reporting GWAS MAF <= 0.40 and |GWAS A1 freq - ref A1 freq| <= 0.10;
  p2np  = p1np SNPs present in panel 2 (489, chr19 LDSC fileset) at the same GRCh37 position with the same allele pair
          (unique), no missing genotype, panel-2 MAF >= 0.01 and |A1 freq panel 2 - panel 1| <= 0.05.
LD       = PLINK19_paper6 --r square --keep-allele-order --extract (503 on the susie set; 489 on p2np); every cell
           asserted returned and non-NA; a one-id mutation must trip the assertion. A1 per SNP recorded and every beta
           expressed on the 503 A1 (489 rows/cols sign-flipped where its A1 differs).
"""
import csv
import gzip
import json
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

CH = "19"
A_RS, B_RS, ANN_RS = "rs73045269", "rs12983775", "rs1800470"
TRAITS = ["P", "CHD_FG", "CAD_C4D"]
REF_MAF_MIN, AF_DIFF_MAX, PAL_MAF_MAX, PAL_FREQ_TOL, P2_AF_DIFF_MAX = 0.01, 0.20, 0.40, 0.10, 0.05
LDD = W.OUT + "/ld"


def plink():
    return W.res("PLINK19_paper6", check=False)


def run_plink(args, tag):
    cmd = [plink()] + args
    r = subprocess.run(cmd, capture_output=True, text=True)
    W.log("[plink %s] %s -> rc %d" % (tag, " ".join(a if " " not in a else '"%s"' % a for a in args), r.returncode), also_print=False)
    if r.returncode != 0:
        raise RuntimeError("plink failed (%s): %s" % (tag, r.stdout[-2000:] + r.stderr[-2000:]))
    return r


def plink_freq(bfile, ids, out):
    with open(out + ".snps", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    run_plink(["--bfile", bfile, "--extract", out + ".snps", "--keep-allele-order", "--freq", "--out", out], "freq")
    rows = []
    with open(out + ".frq") as fh:
        next(fh)
        for line in fh:
            f = line.split()
            rows.append({"snp": f[1], "a1": f[2], "a2": f[3], "maf": float(f[4]), "nchrobs": int(f[5])})
    return rows


def plink_ld_square(bfile, ids, out):
    """r on the A1 allele (--keep-allele-order) for ids, in the PLINK (bim) order; asserts every requested id is
    returned exactly once and every cell is finite. Returns (order, R)."""
    with open(out + ".snps", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    run_plink(["--bfile", bfile, "--extract", out + ".snps", "--keep-allele-order", "--r", "square", "--out", out], "r square")
    fr = plink_freq(bfile, ids, out + "_order")
    order = [r["snp"] for r in fr]
    R = np.loadtxt(out + ".ld", dtype=np.float64, ndmin=2) if len(ids) > 1 else np.ones((1, 1))
    if sorted(order) != sorted(ids) or len(set(order)) != len(order):
        raise AssertionError("PLINK returned %d SNPs for %d requested (missing: %s)" % (
            len(order), len(ids), sorted(set(ids) - set(order))[:5]))
    if R.shape != (len(ids), len(ids)):
        raise AssertionError("PLINK matrix %s for %d requested SNPs" % (R.shape, len(ids)))
    if not np.isfinite(R).all():
        raise AssertionError("PLINK matrix has %d non-finite cells" % int((~np.isfinite(R)).sum()))
    return order, R, {r["snp"]: r for r in fr}


def scan_bim503(lo37, hi37):
    bim = W.res("LD503_bim")
    out = []
    seen19 = False
    with open(bim, "rb") as fh:
        for i, line in enumerate(fh):
            if not line.startswith(b"19\t"):
                if seen19:
                    break
                continue
            seen19 = True
            f = line.split()
            p = int(f[3])
            if lo37 <= p <= hi37:
                out.append({"idx": i, "rsid": f[1].decode(), "pos37": p, "a1": f[4].decode(), "a2": f[5].decode()})
    return out


def read_bim489():
    out = []
    with open(W.res("LD489_chr19_bim"), encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            f = line.split()
            out.append({"idx": i, "rsid": f[1], "pos37": int(f[3]), "a1": f[4], "a2": f[5]})
    return out


def main():
    os.makedirs(LDD, exist_ok=True)
    W.log("==== w2_02_locus (TGFB1) ====")
    ins = [r for r in csv.DictReader(open(W.res("EXP_UKBPPP_INSTRUMENT_TABLE"), encoding="utf-8")) if r["gene"] == "TGFB1"]
    assert len(ins) == 1
    ins = ins[0]
    lo38, hi38 = int(ins["cis_start"]), int(ins["cis_end"])
    assert lo38 == int(ins["gene_start"]) - 500000 and hi38 == int(ins["gene_end"]) + 500000
    W.log("window GRCh38 %s:%d-%d = instrument-table cis window (gene %s-%s +/- 500 kb); UKB-PPP file spans 40319325-42319263"
          % (CH, lo38, hi38, ins["gene_start"], ins["gene_end"]))
    # ---------------- protein
    prot = {}
    ndup = 0
    nrow = 0
    with open(W.res("EXP_UKBPPP_EUR_TGFB1_OID20621"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            nrow += 1
            c, p37, ref, alt = r["ID"].split(":")[:4]
            assert c == CH and r["chrom"] == CH and ref == r["ALLELE0"] and alt == r["ALLELE1"], r["ID"]
            p38 = int(r["pos"])
            if not lo38 <= p38 <= hi38:
                continue
            k = (p38, frozenset((ref, alt)))
            if k in prot:
                prot[k] = None
                ndup += 1
                continue
            prot[k] = {"pos38": p38, "pos37": int(p37), "ea": alt, "oa": ref, "b": float(r["BETA"]), "se": float(r["SE"]),
                       "eaf": float(r["A1FREQ"]), "log10p": float(r["LOG10P"]), "id": r["ID"]}
    prot = {k: v for k, v in prot.items() if v}
    W.log("protein: %d rows in file, %d in window, %d duplicate keys dropped" % (nrow, len(prot), ndup))
    # ---------------- FinnGen CHD (GRCh38) and ieu-a-7 (GRCh37) from the stream-cut extracts
    ext = W.OUT + "/extracts/"
    st = json.load(open(ext + "OUT_FG_R10_I9_CHD.stamp.json"))
    assert st["source_md5"] == W.res_entry("OUT_FG_R10_I9_CHD")["md5"]
    fg = {}
    p = ext + "OUT_FG_R10_I9_CHD__TGFB1.tsv"
    W.assert_md5(p, st["file_md5"]["OUT_FG_R10_I9_CHD__TGFB1.tsv"], "extract FG CHD TGFB1")
    nfg = 0
    for r in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"):
        nfg += 1
        fg.setdefault(int(r["pos"]), []).append({"ea": r["alt"], "oa": r["ref"], "b": float(r["beta"]), "se": float(r["sebeta"]),
                                                 "eaf": float(r["af_alt"]), "p": float(r["pval"]), "rsid": r["rsids"]})
    st = json.load(open(ext + "OUT_CARDIOGRAMplusC4D_ieu-a-7.stamp.json"))
    assert st["source_md5"] == W.res_entry("OUT_CARDIOGRAMplusC4D_ieu-a-7")["md5"]
    p = ext + "OUT_CARDIOGRAMplusC4D_ieu-a-7__TGFB1.tsv"
    W.assert_md5(p, st["file_md5"]["OUT_CARDIOGRAMplusC4D_ieu-a-7__TGFB1.tsv"], "extract ieu-a-7 TGFB1")
    c4 = {}
    nc4 = 0
    with open(p, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d = dict(zip(hdr, f))
            fmt = dict(zip(d["FORMAT"].split(":"), d["ieu-a-7"].split(":")))
            nc4 += 1
            c4.setdefault(int(d["POS"]), []).append({"ea": d["ALT"], "oa": d["REF"], "b": float(fmt["ES"]), "se": float(fmt["SE"]),
                                                     "eaf": float(fmt["AF"]) if fmt.get("AF", ".") not in (".", "") else float("nan"),
                                                     "p": 10 ** (-float(fmt["LP"])), "rsid": d["ID"], "ss": fmt.get("SS", "")})
    W.log("extract rows: FinnGen CHD %d, ieu-a-7 %d" % (nfg, nc4))
    # ---------------- reference panels
    lo37 = min(v["pos37"] for v in prot.values()) - 1000
    hi37 = max(v["pos37"] for v in prot.values()) + 1000
    bim503 = scan_bim503(lo37, hi37)
    W.log("503 bim: %d variants on chr19 %d-%d (GRCh37 extent of the protein rows +/- 1 kb)" % (len(bim503), lo37, hi37))
    by37 = {}
    for b in bim503:
        by37.setdefault(b["pos37"], []).append(b)
    rscount = {}
    for b in bim503:
        rscount[b["rsid"]] = rscount.get(b["rsid"], 0) + 1

    def match(p, cands):
        hits = [(c, W.align(c["ea"], c["oa"], p["ea"], p["oa"])) for c in cands]
        hits = [h for h in hits if h[1]]
        if len(hits) == 1:
            return hits[0], ""
        return (None, 0), ("absent" if not hits else "ambiguous")

    drops = {}

    def drop(k):
        drops[k] = drops.get(k, 0) + 1
    rows = []
    for k, p in prot.items():
        (cf, sf), why = match(p, fg.get(p["pos38"], []))
        if not cf:
            drop("finngen_" + why); continue
        (cc, sc), why = match(p, c4.get(p["pos37"], []))
        if not cc:
            drop("ieu-a-7_" + why); continue
        hits = [(r, W.align(p["ea"], p["oa"], r["a1"], r["a2"])) for r in by37.get(p["pos37"], [])]
        hits = [h for h in hits if h[1]]
        if len(hits) != 1:
            drop("ref503_absent" if not hits else "ref503_ambiguous"); continue
        ref, sref = hits[0]
        if rscount[ref["rsid"]] != 1 or ref["rsid"] == ".":
            drop("ref503_rsid_not_unique"); continue
        rows.append({"p": p, "fg": cf, "sf": sf, "c4": cc, "sc": sc, "ref": ref, "sref": sref})
    W.log("join P x FinnGen x ieu-a-7 x 503: %d of %d protein rows (yield %.3f); drops %s" % (
        len(rows), len(prot), len(rows) / len(prot), drops))
    # yield by class (calibration, hygiene section 1): ieu-a-7 carries no indels and few MAF < 0.01 variants, so the
    # gate is on QC-eligible SNVs (protein MAF >= 0.01), where a broken key would show as a collapse
    joined = {(x["p"]["pos38"], x["p"]["ea"], x["p"]["oa"]) for x in rows}
    cls = {}
    for p in prot.values():
        maf = min(p["eaf"], 1 - p["eaf"])
        c = ("snv" if len(p["ea"]) == len(p["oa"]) == 1 else "indel") + ("_maf_ge_0.01" if maf >= 0.01 else "_maf_lt_0.01")
        n, j = cls.get(c, (0, 0))
        cls[c] = (n + 1, j + ((p["pos38"], p["ea"], p["oa"]) in joined))
    W.log("join yield by class: " + "; ".join("%s %d/%d = %.3f" % (c, j, n, j / n) for c, (n, j) in sorted(cls.items())))
    n, j = cls.get("snv_maf_ge_0.01", (0, 0))
    if not n or j / n < 0.70:
        W.log("ASSERT_FAILED join yield on QC-eligible SNVs %d/%d below 0.70" % (j, n))
    else:
        W.log("join-yield gate (QC-eligible SNVs >= 0.70): pass, %d/%d = %.3f" % (j, n, j / n))
    rows.sort(key=lambda x: x["ref"]["idx"])
    ids = [x["ref"]["rsid"] for x in rows]
    fr = {r["snp"]: r for r in plink_freq(W.res("LD503_bed").replace(".bed", ""), ids, LDD + "/p503_joined")}
    assert set(fr) == set(ids), "503 freq: missing ids"
    # ---------------- engine QC -> susie set
    qc_drops = {}

    def qdrop(k):
        qc_drops[k] = qc_drops.get(k, 0) + 1
    keep = []
    for x in rows:
        f = fr[x["ref"]["rsid"]]
        assert f["a1"] in (x["ref"]["a1"], x["ref"]["a2"]) and f["nchrobs"] == 1006, f
        af1 = f["maf"] if f["a1"] == x["ref"]["a1"] else 1 - f["maf"]       # frequency of the bim A1
        x["af503"] = af1
        maf = min(af1, 1 - af1)
        x["pal"] = W.is_pal(x["ref"]["a1"], x["ref"]["a2"])
        if maf == 0:
            qdrop("monomorphic_ref"); continue
        if maf < REF_MAF_MIN:
            qdrop("ref_maf_lt_0.01"); continue
        if x["pal"] and maf > PAL_MAF_MAX:
            qdrop("palindromic_maf_gt_0.40"); continue
        # GWAS frequencies on the bim A1
        p = x["p"]
        x["fA1"] = {}
        x["fA1"]["P"] = p["eaf"] if x["sref"] == 1 else 1 - p["eaf"]
        fpe = x["fg"]["eaf"] if x["sf"] == 1 else 1 - x["fg"]["eaf"]
        x["fA1"]["CHD_FG"] = fpe if x["sref"] == 1 else 1 - fpe
        if x["c4"]["eaf"] == x["c4"]["eaf"]:
            fpe = x["c4"]["eaf"] if x["sc"] == 1 else 1 - x["c4"]["eaf"]
            x["fA1"]["CAD_C4D"] = fpe if x["sref"] == 1 else 1 - fpe
        else:
            x["fA1"]["CAD_C4D"] = float("nan")
        bad = [t for t, v in x["fA1"].items() if v == v and abs(v - af1) > AF_DIFF_MAX]
        if bad:
            qdrop("af_discordant_" + bad[0]); continue
        keep.append(x)
    W.log("engine QC on the 503 panel: %d -> susie set %d; drops %s" % (len(rows), len(keep), qc_drops))
    rows = keep
    # p1fr palindromes
    for x in rows:
        x["in_p1np"] = not x["pal"]
        if x["pal"]:
            ok = all(v == v and min(v, 1 - v) <= PAL_MAF_MAX and abs(v - x["af503"]) <= PAL_FREQ_TOL for v in x["fA1"].values() if v == v)
            x["in_p1fr"] = ok
        else:
            x["in_p1fr"] = True
    # ---------------- panel 2 (489)
    bim489 = read_bim489()
    by37b = {}
    for b in bim489:
        by37b.setdefault(b["pos37"], []).append(b)
    rc2 = {}
    for b in bim489:
        rc2[b["rsid"]] = rc2.get(b["rsid"], 0) + 1
    p2drops = {}
    cand = []
    for x in rows:
        if not x["in_p1np"]:
            continue
        hits = [b for b in by37b.get(x["ref"]["pos37"], []) if {b["a1"], b["a2"]} == {x["ref"]["a1"], x["ref"]["a2"]}]
        if len(hits) != 1 or rc2[hits[0]["rsid"]] != 1:
            k = "panel2_absent" if not hits else "panel2_ambiguous"
            p2drops[k] = p2drops.get(k, 0) + 1
            x["in_p2np"] = False
            continue
        x["ref2"] = hits[0]
        x["s489"] = 1 if hits[0]["a1"] == x["ref"]["a1"] else -1
        cand.append(x)
    bfile489 = W.res("LD489_chr19_bed").replace(".bed", "")
    fr2 = {r["snp"]: r for r in plink_freq(bfile489, [x["ref2"]["rsid"] for x in cand], LDD + "/p489_cand")}
    for x in cand:
        f = fr2[x["ref2"]["rsid"]]
        af2 = f["maf"] if f["a1"] == x["ref2"]["a1"] else 1 - f["maf"]       # freq of the 489 bim A1
        af2_on503 = af2 if x["s489"] == 1 else 1 - af2
        x["af489"] = af2_on503
        why = None
        if f["nchrobs"] != 978:
            why = "panel2_missing_genotype"
        elif min(af2, 1 - af2) < REF_MAF_MIN:
            why = "panel2_maf_lt_0.01"
        elif abs(af2_on503 - x["af503"]) > P2_AF_DIFF_MAX:
            why = "panel2_af_discordant"
        x["in_p2np"] = why is None
        if why:
            p2drops[why] = p2drops.get(why, 0) + 1
    n_sets = {"susie": len(rows), "p1fr": sum(x["in_p1fr"] for x in rows), "p1np": sum(x["in_p1np"] for x in rows),
              "p2np": sum(x.get("in_p2np", False) for x in rows)}
    W.log("variant sets: %s; panel-2 drops %s; palindromes in susie set %d" % (n_sets, p2drops, sum(x["pal"] for x in rows)))
    # ---------------- PLINK LD matrices + mutation test
    bfile503 = W.res("LD503_bed").replace(".bed", "")
    ids503 = [x["ref"]["rsid"] for x in rows]
    order503, R503, frq503 = plink_ld_square(bfile503, ids503, LDD + "/ld503_susie")
    assert order503 == ids503, "PLINK order differs from bim-index order"
    for x in rows:
        assert frq503[x["ref"]["rsid"]]["a1"] in (x["ref"]["a1"], x["ref"]["a2"])
    p2 = [x for x in rows if x.get("in_p2np")]
    ids489 = [x["ref2"]["rsid"] for x in p2]
    order489, R489raw, _ = plink_ld_square(bfile489, ids489, LDD + "/ld489_p2np")
    pos = {s: i for i, s in enumerate(order489)}
    ix = [pos[s] for s in ids489]
    R489 = R489raw[np.ix_(ix, ix)]
    sg = np.array([x["s489"] for x in p2], dtype=float)
    R489 = R489 * np.outer(sg, sg)
    W.log("LD: 503 matrix %s on the susie set; 489 matrix %s on p2np; all cells finite; 489 A1 flips %d" % (
        R503.shape, R489.shape, int((sg < 0).sum())))
    # mutation: one id altered -> the assertion must trip
    mut = list(ids503[:50])
    mi = mut.index(A_RS) if A_RS in mut else 0
    mut[mi] = mut[mi] + "_MUT"
    try:
        plink_ld_square(bfile503, mut, LDD + "/mutation_test")
        W.log("ASSERT_FAILED LD mutation test: a mutated SNP id did NOT trip the returned/non-NA assertion")
    except (AssertionError, RuntimeError) as e:
        W.log("LD mutation test OK: id %s -> assertion tripped (%s)" % (mut[mi], str(e)[:120]))
    # ---------------- write the locus table and LD binaries
    trait_n = {}
    v = [1.0 / (2 * x["p"]["eaf"] * (1 - x["p"]["eaf"]) * x["p"]["se"] ** 2) for x in rows if 0.05 < min(x["p"]["eaf"], 1 - x["p"]["eaf"])]
    trait_n["P"] = float(np.median(v))
    man = {r["phenocode"]: r for r in csv.DictReader(open(W.res("OUT_FG_R10_MANIFEST"), encoding="utf-8"), delimiter="\t")}
    nca, nco = int(man["I9_CHD"]["num_cases"]), int(man["I9_CHD"]["num_controls"])
    trait_n["CHD_FG"] = 4.0 * nca * nco / (nca + nco)
    alt = [r for r in csv.DictReader(open(W.res("D_ALT_OUTCOMES_SUMMARY"), encoding="utf-8")) if r["gwas_id"] == "ieu-a-7"]
    assert len(alt) == 1
    ca, co = int(alt[0]["n_cases"]), int(alt[0]["n_controls"])
    trait_n["CAD_C4D"] = 4.0 * ca * co / (ca + co)
    casefrac = {"CHD_FG": nca / (nca + nco), "CAD_C4D": ca / (ca + co)}
    ntot = {"P": trait_n["P"], "CHD_FG": nca + nco, "CAD_C4D": ca + co}
    W.log("N: P n_hat (engine rule, median 1/(2f(1-f)se^2), MAF > 0.05, %d SNPs) = %.0f (instrument-table samplesize 33988); "
          "CHD_FG cases %d controls %d -> N_eff %.0f (RESOURCES n_info 166436), case fraction %.5f; "
          "CAD_C4D cases %d controls %d -> N_eff %.0f (RESOURCES n_info 162973), case fraction %.5f" % (
              len(v), trait_n["P"], nca, nco, trait_n["CHD_FG"], casefrac["CHD_FG"], ca, co, trait_n["CAD_C4D"], casefrac["CAD_C4D"]))
    out = []
    for i, x in enumerate(rows):
        p, r = x["p"], x["ref"]
        s = x["sref"]
        o = {"i": i, "snp": "%s:%d:%s:%s" % (CH, r["pos37"], r["a1"], r["a2"]), "rsid503": r["rsid"], "bed_idx503": r["idx"],
             "rsid489": x["ref2"]["rsid"] if x.get("in_p2np") else "", "bed_idx489": x["ref2"]["idx"] if x.get("in_p2np") else "",
             "pos37": r["pos37"], "pos38": p["pos38"], "a1": r["a1"], "a2": r["a2"], "ukb_id": p["id"],
             "af503": x["af503"], "af489": x.get("af489", float("nan")), "palindromic": int(x["pal"]),
             "in_susie": 1, "in_p1fr": int(x["in_p1fr"]), "in_p1np": int(x["in_p1np"]), "in_p2np": int(bool(x.get("in_p2np"))),
             "s489": x.get("s489", "")}
        for t, b, se, sg2, lp in (("P", p["b"], p["se"], s, p["log10p"]),
                                  ("CHD_FG", x["fg"]["b"], x["fg"]["se"], s * x["sf"], -math.log10(max(x["fg"]["p"], 1e-300))),
                                  ("CAD_C4D", x["c4"]["b"], x["c4"]["se"], s * x["sc"], -math.log10(max(x["c4"]["p"], 1e-300)))):
            o["b__" + t] = sg2 * b
            o["se__" + t] = se
            o["z__" + t] = sg2 * b / se
            o["log10p__" + t] = lp
            o["fA1__" + t] = x["fA1"][t]
        o["rsid_fg"] = x["fg"]["rsid"]
        o["rsid_c4d"] = x["c4"]["rsid"]
        out.append(o)
    cols = list(out[0].keys())
    W.write_tsv(W.OUT + "/tgfb1_locus.tsv", cols, out)
    R503.astype("<f8").tofile(LDD + "/R503_susie.f8")
    R489.astype("<f8").tofile(LDD + "/R489_p2np.f8")
    # rsID concordance (diagnostic: joins are positional)
    agree = sum(1 for o in out if o["rsid_fg"] and o["rsid503"] in o["rsid_fg"].split(","))
    agree2 = sum(1 for o in out if o["rsid_c4d"] == o["rsid503"])
    meta = {"window_grch38": [lo38, hi38], "chr": CH, "n": n_sets, "trait_n": trait_n, "trait_ntot": ntot, "case_fraction": casefrac,
            "join_drops": drops, "qc_drops": qc_drops, "panel2_drops": p2drops,
            "rsid_concordance_finngen": "%d/%d" % (agree, len(out)), "rsid_concordance_ieu": "%d/%d" % (agree2, len(out)),
            "R503": "ld/R503_susie.f8 (float64 LE, %d x %d, order = tgfb1_locus.tsv rows)" % R503.shape,
            "R489": "ld/R489_p2np.f8 (float64 LE, %d x %d, order = tgfb1_locus.tsv rows with in_p2np == 1; aligned to the 503 A1)" % R489.shape,
            "md5": {"tgfb1_locus.tsv": W.md5_file(W.OUT + "/tgfb1_locus.tsv"), "R503_susie.f8": W.md5_file(LDD + "/R503_susie.f8"),
                    "R489_p2np.f8": W.md5_file(LDD + "/R489_p2np.f8")}}
    json.dump(meta, open(W.OUT + "/tgfb1_locus_meta.json", "w"), indent=1)
    W.log("rsID concordance (diagnostic): FinnGen %s, ieu-a-7 %s" % (meta["rsid_concordance_finngen"], meta["rsid_concordance_ieu"]))
    for rs in (A_RS, B_RS, ANN_RS):
        hit = [o for o in out if o["rsid503"] == rs]
        if hit:
            o = hit[0]
            W.log("  %s: %s in sets susie/p1fr/p1np/p2np = 1/%d/%d/%d; z P %.2f CHD_FG %.2f CAD_C4D %.2f" % (
                rs, o["snp"], o["in_p1fr"], o["in_p1np"], o["in_p2np"], o["z__P"], o["z__CHD_FG"], o["z__CAD_C4D"]))
        else:
            W.log("  %s: NOT in the susie set" % rs + ("" if rs == ANN_RS else " -- ASSERT_FAILED signal SNP missing"))


if __name__ == "__main__":
    main()
