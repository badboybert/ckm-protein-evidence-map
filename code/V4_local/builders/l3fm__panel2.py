#!/usr/bin/env python
"""l3fm__panel2.py -- the SNP-set variants every sensitivity stage reads, and the SECOND LD panel's genotype extract
(promoting the independent verifier's two decisive methods: a second 1000G EUR panel, and dropping palindromic SNPs).
OFFLINE; reads the locus harmonised table + the panel-2 per-chromosome PLINK fileset; writes staging/l3-finemap/<locus>/
snpsets/, raw/l3-finemap/<locus>__ref2_A1dosage_489xP.i8 and its RAW_MANIFEST row.

  python -B V3/V4_local/builders/l3fm__panel2.py [--locus L ...] [--force]

Per locus (declared rules; resumable: skipped when PANEL2_STAMP.json still matches the harmonised table, the panel-2
.bim/.bed sizes and this script):
  p1np  = every harmonised SNP that is NOT palindromic (A/T or C/G), panel 1 (paper 6 g1000_eur, 503) LD.
  p2np  = the p1np SNPs present in the panel-2 .bim (1000G.EUR.QC.<chr>, 489) at the same GRCh37 position with the same
          allele pair (unique match; a panel-2 A1 equal to the harmonised A2 is recoded 2 - dosage), no missing genotype,
          panel-2 MAF >= 0.01 and |A1 frequency panel 2 - panel 1| <= P2_AF_DIFF_MAX (0.05; the panels share 489 samples,
          so a larger difference can only be an allele-coding error). Every drop is counted.
  Output snpsets/<set>.tsv: snp, row (0-based row of harmonised.tsv.gz), af (A1 = harmonised a1 frequency in that panel).
CAUTION (measured here and stamped): panel 2's 489 IIDs are all among panel 1's 503; the second panel is a perturbation of
the first (14 fewer people, different variant QC), not an independent population sample.
"""
import argparse
import csv
import gzip
import json
import os
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import l3fm__common as C  # noqa: E402
import l3fm__prep as P  # noqa: E402

LOCI = list(C.LOCI) + ["UMOD_eGFR_control"]


def geno2_path(loc):
    return os.path.join(C.RAW_FM, "%s__ref2_A1dosage_489xP.i8" % loc)


def build(loc, log=print):
    ldir = os.path.join(C.STAGE_FM, loc)
    hpath = os.path.join(ldir, "harmonised.tsv.gz")
    with gzip.open(hpath, "rt", encoding="utf-8") as fh:
        h = list(csv.DictReader(fh, delimiter="\t"))
    ch = h[0]["chr"]
    pre2 = C.REF2_PREFIX_FMT % ch
    if C.fam_n(pre2) != C.REF2_N:
        C.refuse("panel-2 .fam for chr%s has %d samples, expected %d" % (ch, C.fam_n(pre2), C.REF2_N))
    iid1, iid2 = set(C.fam_iids(C.REF_PREFIX)), C.fam_iids(pre2)
    overlap = sum(1 for i in iid2 if i in iid1)
    drops = {}

    def drop(k):
        drops[k] = drops.get(k, 0) + 1
    p1np = []
    for j, r in enumerate(h):
        if r["palindromic"] == "1":
            drop("palindromic")
            continue
        p1np.append((j, r))
    lo37, hi37 = min(int(r["pos37"]) for r in h), max(int(r["pos37"]) for r in h)
    bim2 = C.ref2_bim_region(ch, lo37, hi37)
    bypos = {}
    for b in bim2:
        bypos.setdefault(b["pos37"], []).append(b)
    cand = []
    for j, r in p1np:
        hits = [b for b in bypos.get(int(r["pos37"]), []) if {b["a1"], b["a2"]} == {r["a1"], r["a2"]}]
        if len(hits) != 1:
            drop("panel2_absent" if not hits else "panel2_ambiguous")
            continue
        b = hits[0]
        cand.append((j, r, b, 1 if b["a1"] == r["a1"] else -1))
    G = C.read_bed_dosage([c[2]["idx"] for c in cand], n=C.REF2_N, prefix=pre2)
    for k, c in enumerate(cand):
        if c[3] == -1:
            G[:, k] = 2 - G[:, k]
    keep = []
    rs_same = 0
    for k, (j, r, b, s) in enumerate(cand):
        g = G[:, k]
        if np.isnan(g).any():
            drop("panel2_missing_genotype"); continue
        f2 = g.mean() / 2
        if min(f2, 1 - f2) < C.REF_MAF_MIN or g.std() == 0:
            drop("panel2_maf_lt_0.01"); continue
        if abs(f2 - float(r["af_ref"])) > C.P2_AF_DIFF_MAX:
            drop("panel2_af_discordant"); continue
        keep.append(k)
        rs_same += b["rsid"] == r["rsid_ref"]
    G2 = G[:, keep]
    sdir = os.path.join(ldir, "snpsets")
    os.makedirs(sdir, exist_ok=True)
    C.write_tsv(os.path.join(sdir, "p1np.tsv"), ["snp", "row", "af"], [{"snp": r["snp"], "row": j, "af": r["af_ref"]} for j, r in p1np])
    C.write_tsv(os.path.join(sdir, "p2np.tsv"), ["snp", "row", "af"],
                [{"snp": cand[k][1]["snp"], "row": cand[k][0], "af": "%.6f" % (G2[:, i].mean() / 2)} for i, k in enumerate(keep)])
    gp = geno2_path(loc)
    G2.astype(np.int8).tofile(gp + ".part")
    os.replace(gp + ".part", gp)
    P.register([P.entry(gp, "ref2_extract", "1000G.EUR.QC chr%s (LDSC, 489 = subset of panel-1 503) A1 dosage, int8, 489 x %d, "
                                            "SNP order = %s/snpsets/p2np.tsv" % (ch, len(keep), loc))])
    af1 = np.array([float(cand[k][1]["af_ref"]) for k in keep])
    af2 = G2.mean(axis=0) / 2
    stamp = {"locus": loc, "chr": ch, "n_harmonised": len(h), "n_p1np": len(p1np), "n_p2np": len(keep), "drops": drops,
             "panel2_prefix": C.rel(pre2), "panel2_n": C.REF2_N, "panel2_iids_in_panel1": "%d/%d" % (overlap, len(iid2)),
             "rsid_concordance_p2np": "%d/%d" % (rs_same, len(keep)), "max_abs_af_diff": float(np.max(np.abs(af2 - af1))),
             "harmonised_md5": C.md5_file(hpath), "script_md5": C.md5_file(os.path.abspath(__file__)),
             "common_md5": C.md5_file(os.path.join(HERE, "l3fm__common.py")),
             "panel2_bed_bytes": os.path.getsize(pre2 + ".bed"), "panel2_bim_bytes": os.path.getsize(pre2 + ".bim"),
             "outputs": {"snpsets/p1np.tsv": C.md5_file(os.path.join(sdir, "p1np.tsv")), "snpsets/p2np.tsv": C.md5_file(os.path.join(sdir, "p2np.tsv")),
                         C.rel(gp): C.md5_file(gp)},
             "built_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    C.write_json(os.path.join(sdir, "PANEL2_STAMP.json"), stamp)
    log("[panel2] %-18s harmonised %d -> p1np %d -> p2np %d; IIDs %s; rsID %s; max|dAF| %.4f; drops %s" % (
        loc, len(h), len(p1np), len(keep), stamp["panel2_iids_in_panel1"], stamp["rsid_concordance_p2np"], stamp["max_abs_af_diff"], drops))
    return stamp


def fresh(loc):
    sp = os.path.join(C.STAGE_FM, loc, "snpsets", "PANEL2_STAMP.json")
    if not os.path.exists(sp):
        return False
    s = json.load(open(sp, encoding="utf-8"))
    pre2 = C.REF2_PREFIX_FMT % s["chr"]
    outs_ok = all(os.path.exists(os.path.join(C.STAGE_FM, loc, k) if not k.startswith("V3/") else os.path.join(C.ROOT, k)) and
                  C.md5_file(os.path.join(C.STAGE_FM, loc, k) if not k.startswith("V3/") else os.path.join(C.ROOT, k)) == v
                  for k, v in s["outputs"].items())
    return (outs_ok and s["harmonised_md5"] == C.md5_file(os.path.join(C.STAGE_FM, loc, "harmonised.tsv.gz"))
            and s["script_md5"] == C.md5_file(os.path.abspath(__file__)) and s["common_md5"] == C.md5_file(os.path.join(HERE, "l3fm__common.py"))
            and s["panel2_bed_bytes"] == os.path.getsize(pre2 + ".bed") and s["panel2_bim_bytes"] == os.path.getsize(pre2 + ".bim"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", nargs="*")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    for loc in a.locus or LOCI:
        if not a.force and fresh(loc):
            print("[panel2] %-18s up to date (skip)" % loc, flush=True)
            continue
        build(loc, lambda m: print(m, flush=True))


if __name__ == "__main__":
    main()
