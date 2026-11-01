#!/usr/bin/env python
"""l3__build_bridge_screen_eur.py -- L3 protein<->metabolite bridge screen on UKB-PPP EUROPEAN DISCOVERY protein arms
(user decision D1, 2026-09-26). OFFLINE (no network), no side effects on import, never writes outside --outroot.

Why this exists. The shipped screen (V3/source_data/l3_metabolite/build_bridge_screen.py -> S14) took 20 of its 23
protein arms from the UKB-PPP AFRICAN folder (syn51365304, members 'combined_chr*', N = 888) and ABO/APOE/ARG1 from a
deCODE cache. This builder re-runs the SAME engine -- coloc() below is copied verbatim from the shipped builder
(allele or strand matching, n >= 20, coloc.abf p1 = p2 = 1e-4, p12 = 1e-5, dual rule H4 >= 0.80 and H4/H3 >= 5) -- on:
  * protein arms, --arm eur (default): UKB-PPP European discovery (syn51365303), ea = ALLELE1, oa = ALLELE0, restricted
      to gene +-500 kb, keyed (chr, pos) last-row-wins exactly as shipped:
        17 proteins  weight-loss.omics/data/proteome/pqtl_regional/<OID>_chr<c>_cis.tsv   (read-only; lead +-1 Mb)
         6 proteins  V3/V4_local/raw/d1-bridge-eur/ukbppp_eur_windows/<OID>_chr<c>_cis.tsv (l3__fetch_bridge_eur_inputs.py)
  * protein arms, --arm shipped (control C1): V3/source_data/l3_metabolite/protein_regional/<GENE>.tsv, else
      protein_regional_cache/<GENE>.tsv -- the shipped load_protein() precedence.
  * metabolite arms (both modes): Chen-CLSA EUROPEAN deposit = lowest accession, asserted inside
      [GCST90199621, GCST90201020]; gene +-500 kb; pre-fetched by l3__fetch_bridge_eur_inputs.py to
      V3/V4_local/raw/d1-bridge-eur/chen_windows/ and parsed exactly as the shipped stream_metab().
  * the 20 worklist keys marked UNMAPPED are mapped by recover_unmapped_bridges.norm() (copied verbatim) to the lowest
      accession of the normalised name (mapping = normalised_name).

  python -B V3/V4_local/builders/l3__build_bridge_screen_eur.py --arm eur     --outroot V3/V4_local/staging/d1-bridge-eur
  python -B V3/V4_local/builders/l3__build_bridge_screen_eur.py --arm shipped --outroot V3/V4_local/staging/d1-bridge-eur/control_shipped

Outputs in --outroot: bridge_screen_results.csv (84 rows, worklist order; shipped columns first, then oid, protein_arm,
mapping, prot_window_md5, metab_window_md5, pos_collisions), bridge_screen_pp_full.csv (all five PPs unrounded + flags),
BUILD_STAMP.json. REFUSES (no silent STREAM_FAIL / PENDING): a missing window, an md5 differing from RAW_MANIFEST.tsv, a
non-European GCST, vendored-engine drift, a key set that is not the 84 worklist keys.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sys
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))          # .../CKM papers
WORKLIST = os.path.join(ROOT, "V3", "V4_local", "source_data", "l3_metabolite", "bridge_worklist.csv")
CHEN_MAP = os.path.join(ROOT, "weight-loss.omics", "data", "metabolome", "chen_clsa_metabolite_gcst_map.tsv")
V4_INSTR = os.path.join(ROOT, "V3", "V4", "source_data", "tables_from_rds", "phase1_instruments__ukb_cis_instruments_with_r2.csv")
TIERED = os.path.join(ROOT, "V3", "v3.2", "_engine_root", "results", "repro_phase1_hg38", "_fig_results",
                      "phase4_annotation", "table_s_pair_evidence_tiered.csv")
WL_PQTL = os.path.join(ROOT, "weight-loss.omics", "data", "proteome", "pqtl_regional")
RAW = os.path.join(ROOT, "V3", "V4_local", "raw", "d1-bridge-eur")
RAW_MANIFEST = os.path.join(RAW, "RAW_MANIFEST.tsv")
EUR_WIN_DIR = os.path.join(RAW, "ukbppp_eur_windows")
CHEN_WIN_DIR = os.path.join(RAW, "chen_windows")
SHIPPED_FETCHED = os.path.join(ROOT, "V3", "source_data", "l3_metabolite", "protein_regional")
SHIPPED_CACHE = os.path.join(ROOT, "V3", "source_data", "l3_metabolite", "protein_regional_cache")
VENDOR = os.path.join(HERE, "_vendor")
# vendored engine copies must stay byte-identical to the modules the shipped builder imported
VENDOR_MD5 = {"l3__mr_coloc.py": "fa6db396645e411504c4825864f75dd4",       # paper 3.eas/analysis/engine_ref/mr_coloc.py
              "l3__tabix_remote.py": "391fcc8c5385c72a8a095f50b2f69210"}   # paper 3.metaboQTL/src/metabolome/tabix_remote.py

WIN = 500_000                          # bridge window = gene +- 500 kb (shipped)
EUR_BLOCK = (90199621, 90201020)       # Chen-CLSA European deposits (GWAS Catalog initialSampleSize; verifier 2026-09-26)
DL_GENES = {"AK1", "ANXA4", "BRAP", "CCND2", "KCTD5", "SMPD3"}   # no local European window; extracted from syn51365303 tars
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
OUT_COLS = ["protein", "metabolite", "metab_gcst", "nsnp", "H4", "H3", "support", "protein_source",
            "oid", "protein_arm", "mapping", "prot_window_md5", "metab_window_md5", "pos_collisions"]


class Refuse(SystemExit):
    pass


def refuse(msg):
    raise Refuse("REFUSED: " + msg)


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 22), b""):
            h.update(ch)
    return h.hexdigest()


def load_engine():
    """import the vendored coloc engine after asserting it is byte-identical to the shipped engine."""
    for f, want in VENDOR_MD5.items():
        got = md5_file(os.path.join(VENDOR, f))
        if got != want:
            refuse("vendored %s md5 %s != %s (engine drift)" % (f, got, want))
    if VENDOR not in sys.path:
        sys.path.insert(0, VENDOR)
    import l3__mr_coloc as E
    return E


# ---------------------------------------------------------------- units (single implementation; fetcher + checker import it)
def norm(s):   # verbatim from V3/source_data/l3_metabolite/recover_unmapped_bridges.py
    s=s.lower().replace(" levels",""); s=re.sub(r"\*+","",s); return re.sub(r"[^a-z0-9]+"," ",s).strip()


def gcst_num(g):
    if not re.fullmatch(r"GCST\d{8}", g or ""):
        refuse("malformed GCST %r" % g)
    return int(g[4:])


def assert_european_gcst(g):
    """C9: every metabolite arm must be a Chen-CLSA European deposit."""
    n = gcst_num(g)
    if not (EUR_BLOCK[0] <= n <= EUR_BLOCK[1]):
        refuse("metabolite accession %s lies outside the European block GCST%d-GCST%d" % (g, EUR_BLOCK[0], EUR_BLOCK[1]))
    return g


def gcst_url(g):   # verbatim logic from the shipped builder
    n=int(g.replace("GCST","")); lo=((n-1)//1000)*1000+1
    return f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{lo:08d}-GCST{lo+999:08d}/{g}/harmonised/{g}.h.tsv.gz"


def load_oids():
    """gene -> OID from the V4 instruments table (one OID per gene, asserted)."""
    m = defaultdict(set)
    for r in csv.DictReader(io.open(V4_INSTR, encoding="utf-8")):
        m[r["gene"]].add(r["oid"])
    return m


def derive_units():
    """the 84 worklist keys, in worklist order, with GCST (worklist or normalised lowest accession) and window."""
    wl = list(csv.DictReader(io.open(WORKLIST, encoding="utf-8-sig")))
    idx = defaultdict(set)
    for r in csv.DictReader(io.open(CHEN_MAP, encoding="utf-8"), delimiter="\t"):
        idx[norm(r["trait"])].add(r["gcst"])
    oids = load_oids()
    units, seen = [], set()
    for r in wl:
        key = (r["protein"], r["metabolite"])
        if key in seen:
            refuse("worklist key duplicated: %r" % (key,))
        seen.add(key)
        g = r["protein"]
        if len(oids.get(g, ())) != 1:
            refuse("gene %s has %d OIDs in the V4 instruments table" % (g, len(oids.get(g, ()))))
        if r["metab_gcst"] != "UNMAPPED":
            gc, mapping = r["metab_gcst"], "worklist"
        else:
            accs = sorted(idx.get(norm(r["metabolite"]), []))
            if not accs:
                refuse("UNMAPPED key %r has no normalised-name accession" % (key,))
            gc, mapping = accs[0], "normalised_name"
        assert_european_gcst(gc)
        ch = str(r["gene_chr"])
        gs, ge = int(r["gene_start"]), int(r["gene_end"])
        units.append(dict(protein=g, metabolite=r["metabolite"], gcst=gc, mapping=mapping, chr=ch,
                          gene_start=gs, gene_end=ge, lo=gs - WIN, hi=ge + WIN, oid=next(iter(oids[g]))))
    if len(units) != 84 or len({u["protein"] for u in units}) != 23:
        refuse("worklist gives %d keys over %d proteins, expected 84 over 23" % (len(units), len({u["protein"] for u in units})))
    return units


def metab_window_path(u):
    return os.path.join(CHEN_WIN_DIR, "%s__chr%s_%d_%d.tsv.gz" % (u["gcst"], u["chr"], u["lo"], u["hi"]))


def eur_window_path(gene, oid, ch):
    if gene in DL_GENES:
        return os.path.join(EUR_WIN_DIR, "%s_chr%s_cis.tsv" % (oid, ch))
    return os.path.join(WL_PQTL, "%s_chr%s_cis.tsv" % (oid, ch))


def read_manifest(path=RAW_MANIFEST):
    if not os.path.exists(path):
        return {}
    return {r["path"]: r for r in csv.DictReader(io.open(path, encoding="utf-8"), delimiter="\t")}


def rel(p):
    """path relative to the CKM papers root; absolute (forward slashes) when outside it or on another drive."""
    a = os.path.abspath(p)
    try:
        r = os.path.relpath(a, ROOT)
    except ValueError:                      # different drive (e.g. <D>/data/...)
        return a.replace("\\", "/")
    return a.replace("\\", "/") if r.startswith("..") else r.replace("\\", "/")


def checked_md5(path, manifest):
    """md5 of an input, refusing when it is missing or differs from RAW_MANIFEST (staleness / tamper guard)."""
    if not os.path.exists(path):
        refuse("missing input window %s" % path)
    got = md5_file(path)
    row = manifest.get(rel(path))
    if row is None:
        refuse("input %s is not registered in RAW_MANIFEST.tsv" % rel(path))
    if row["md5"] != got:
        refuse("input %s md5 %s != RAW_MANIFEST %s" % (rel(path), got, row["md5"]))
    return got


# ---------------------------------------------------------------- protein arms
def load_eur_protein(path, ch, lo, hi):
    """UKB-PPP 9-column window -> {(chr,pos): allele/beta/se}; ea = ALLELE1, oa = ALLELE0; last row wins (as shipped)."""
    d, coll, skipped = {}, 0, 0
    with io.open(path, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"ID", "chrom", "pos", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"}
        if not need <= set(rd.fieldnames or []):
            refuse("%s lacks the 9-column UKB-PPP schema (has %s)" % (path, rd.fieldnames))
        for r in rd:
            if str(r["chrom"]) != ch:
                refuse("%s carries chromosome %s, expected %s" % (path, r["chrom"], ch))
            p = int(r["pos"])
            if p < lo or p > hi:
                continue
            try:
                v = {"ea": r["ALLELE1"].upper(), "oa": r["ALLELE0"].upper(), "b": float(r["BETA"]), "se": float(r["SE"])}
            except ValueError:
                skipped += 1
                continue
            if (ch, p) in d:
                coll += 1
            d[(ch, p)] = v
    return d, coll, skipped


def load_shipped_protein(gene):
    """the shipped load_protein() precedence: protein_regional/<GENE>.tsv ('fetched') then protein_regional_cache/ ('cache').
    The deCODE genome-wide fallback is not reproduced: every shipped protein resolved to one of these two files."""
    fet = os.path.join(SHIPPED_FETCHED, gene + ".tsv")
    cache = os.path.join(SHIPPED_CACHE, gene + ".tsv")
    coll = 0
    if os.path.exists(fet):
        d = {}
        for r in csv.DictReader(io.open(fet, encoding="utf-8"), delimiter="\t"):
            try:
                k = (str(r["chr"]), int(r["pos"]))
                v = {"ea": r["ea"].upper(), "oa": r["oa"].upper(), "b": float(r["beta"]), "se": float(r["se"])}
            except (KeyError, ValueError):
                continue
            coll += k in d
            d[k] = v
        return d, "fetched", fet, coll
    if os.path.exists(cache):
        d = {}
        for r in csv.DictReader(io.open(cache, encoding="utf-8"), delimiter="\t"):
            k = (r["chr"], int(r["pos"]))
            coll += k in d
            d[k] = {"ea": r["ea"], "oa": r["oa"], "b": float(r["beta"]), "se": float(r["se"])}
        return d, "cache", cache, coll
    refuse("shipped arm for %s: neither %s nor %s exists" % (gene, fet, cache))


# ---------------------------------------------------------------- metabolite arms
def load_metab_window(path):
    """parse a pre-fetched Chen window exactly as the shipped stream_metab() parsed the streamed rows."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    hdr = lines[0].rstrip("\r").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    out, coll, skipped, nrow = {}, 0, 0, 0
    for line in lines[1:]:
        if not line:
            continue
        nrow += 1
        f = line.split("\t")
        try:
            c = f[ci["chromosome"]].replace("chr", ""); p = int(f[ci["base_pair_location"]])
            b = float(f[ci["beta"]]); se = float(f[ci["standard_error"]])
        except (KeyError, ValueError, IndexError):
            skipped += 1
            continue
        coll += (c, p) in out
        out[(c, p)] = {"ea": f[ci["effect_allele"]].upper(), "oa": f[ci["other_allele"]].upper(), "b": b, "se": se}
    return out, coll, skipped, nrow


# ---------------------------------------------------------------- engine
def coloc(prot, metab, E):
    # verbatim from V3/source_data/l3_metabolite/build_bridge_screen.py (engine module passed in instead of a global)
    import numpy as np
    b1=[];se1=[];b2=[];se2=[]
    for k,p in prot.items():
        m=metab.get(k)
        if not m: continue
        if (p["ea"],p["oa"])==(m["ea"],m["oa"]): mb=m["b"]
        elif (p["ea"],p["oa"])==(m["oa"],m["ea"]): mb=-m["b"]
        elif (p["ea"],p["oa"])==(COMP.get(m["ea"]),COMP.get(m["oa"])): mb=m["b"]
        elif (p["ea"],p["oa"])==(COMP.get(m["oa"]),COMP.get(m["ea"])): mb=-m["b"]
        else: continue
        b1.append(p["b"]);se1.append(p["se"]);b2.append(mb);se2.append(m["se"])
    n=len(b1)
    if n<20: return {"nsnp":n,"H4":"","H3":"","support":"TOO_FEW"}
    pp=E.coloc_abf(np.array(b1),np.array(se1),np.array(b2),np.array(se2),p1=1e-4,p2=1e-4,p12=1e-5)
    h4,h3=pp["H4"],pp["H3"]
    return {"nsnp":n,"H4":round(h4,3),"H3":round(h3,3),"support":bool(h4>=0.80 and (h4/h3 if h3>0 else 9e9)>=5)}


def pp_full(prot, metab, E):
    """all five posteriors on the same alignment as coloc(); asserted consistent with coloc() for every row."""
    import numpy as np
    b1, se1, b2, se2 = [], [], [], []
    for k, p in prot.items():
        m = metab.get(k)
        if not m:
            continue
        if (p["ea"], p["oa"]) == (m["ea"], m["oa"]): mb = m["b"]
        elif (p["ea"], p["oa"]) == (m["oa"], m["ea"]): mb = -m["b"]
        elif (p["ea"], p["oa"]) == (COMP.get(m["ea"]), COMP.get(m["oa"])): mb = m["b"]
        elif (p["ea"], p["oa"]) == (COMP.get(m["oa"]), COMP.get(m["ea"])): mb = -m["b"]
        else: continue
        b1.append(p["b"]); se1.append(p["se"]); b2.append(mb); se2.append(m["se"])
    if len(b1) < 20:
        return len(b1), None
    return len(b1), E.coloc_abf(np.array(b1), np.array(se1), np.array(b2), np.array(se2), p1=1e-4, p2=1e-4, p12=1e-5)


# ---------------------------------------------------------------- flags (diagnostic; not used by the verdict)
def mhc_genes():
    out = set()
    for r in csv.DictReader(io.open(TIERED, encoding="utf-8")):
        if str(r.get("abf_mhc_caution", "")).upper() == "TRUE":
            out.add(r["gene"])
    return out


def shared_window_map(units):
    """protein pairs whose bridge windows overlap AND that share a metabolite key (multi-signal / shared-locus caution)."""
    win = {}
    for u in units:
        win[u["protein"]] = (u["chr"], u["lo"], u["hi"])
    mets = defaultdict(set)
    for u in units:
        mets[u["protein"]].add(u["gcst"])
    share = defaultdict(set)
    ps = sorted(win)
    for i, a in enumerate(ps):
        for b in ps[i + 1:]:
            if win[a][0] == win[b][0] and win[a][1] <= win[b][2] and win[b][1] <= win[a][2]:
                share[a].add(b); share[b].add(a)
    return share, mets


# ---------------------------------------------------------------- build
def build(arm, outroot, prot_override=None, units=None, log=print):
    """run the screen; prot_override maps gene -> window path (used ONLY by the checker's mutation ladder)."""
    E = load_engine()
    units = units if units is not None else derive_units()
    for u in units:
        assert_european_gcst(u["gcst"])
    manifest = read_manifest()
    mhc = mhc_genes()
    share, mets = shared_window_map(units)
    pcache = {}
    rows, full = [], []
    t0 = time.time()
    for i, u in enumerate(units, 1):
        g = u["protein"]
        if g not in pcache:
            if arm == "eur" or (prot_override and g in prot_override):
                path = (prot_override or {}).get(g) or eur_window_path(g, u["oid"], u["chr"])
                pm = md5_file(path) if (prot_override and g in prot_override) else checked_md5(path, manifest)
                d, coll, skipped = load_eur_protein(path, u["chr"], u["lo"], u["hi"])
                src = "eur_tar_extract" if g in DL_GENES else "wlomics_pqtl_regional"
                parm = "UKB-PPP_EUR_discovery"
                if prot_override and g in prot_override:
                    src, parm = "override", "OVERRIDE:" + os.path.basename(path)
            else:
                d, src, path, coll = load_shipped_protein(g)
                skipped = 0
                pm = md5_file(path)
                parm = "shipped_deCODE" if src == "cache" else "shipped_UKB-PPP_AFR"
            if not d:
                refuse("protein arm %s loaded 0 rows from %s" % (g, path))
            pcache[g] = (d, src, pm, parm, coll, skipped, rel(path))
        d, src, pm, parm, pcoll, pskip, ppath = pcache[g]
        mpath = metab_window_path(u)
        mm = checked_md5(mpath, manifest)
        metab, mcoll, mskip, mrows = load_metab_window(mpath)
        if not metab:
            refuse("metabolite window %s parsed to 0 rows" % mpath)
        r = coloc(d, metab, E)
        n2, pp = pp_full(d, metab, E)
        if n2 != r["nsnp"] or (pp is not None and (round(pp["H4"], 3) != r["H4"] or round(pp["H3"], 3) != r["H3"])):
            refuse("pp_full disagrees with coloc() for %s/%s" % (g, u["metabolite"]))
        rows.append({"protein": g, "metabolite": u["metabolite"], "metab_gcst": u["gcst"], **r, "protein_source": src,
                     "oid": u["oid"], "protein_arm": parm, "mapping": u["mapping"], "prot_window_md5": pm,
                     "metab_window_md5": mm, "pos_collisions": "prot=%d;metab=%d" % (pcoll, mcoll)})
        sw = sorted(p for p in share.get(g, ()) if u["gcst"] in mets[p])
        full.append({"protein": g, "metabolite": u["metabolite"], "metab_gcst": u["gcst"], "nsnp": r["nsnp"],
                     **({k: "%.6g" % pp[k] for k in ("H0", "H1", "H2", "H3", "H4")} if pp else {k: "" for k in ("H0", "H1", "H2", "H3", "H4")}),
                     "dominant": (max(pp, key=pp.get) if pp else ""), "support": r["support"],
                     "mhc_caution": g in mhc, "shared_window_same_metabolite": ";".join(sw),
                     "prot_window": ppath, "n_prot_window": len(d), "n_metab_rows": mrows, "n_metab_keys": len(metab),
                     "prot_skipped": pskip, "metab_skipped": mskip})
        log("[%2d/%d] %-7s x %-40s %s nsnp=%-5s H4=%-5s H3=%-5s support=%s%s" % (
            i, len(units), g, u["metabolite"][:40], u["gcst"], r["nsnp"], r["H4"], r["H3"], r["support"],
            "  <== BRIDGE" if r["support"] is True else ""))
    # C4 inline: the 84 keys once each, no transient state
    keys = [(r["protein"], r["metabolite"]) for r in rows]
    if len(keys) != len(set(keys)) or len(keys) != len(units):
        refuse("output keys are not the worklist keys once each")
    bad = [k for k, r in zip(keys, rows) if r["support"] not in (True, False)]
    if bad:
        refuse("non-verdict rows (TOO_FEW or transient): %s" % bad)
    os.makedirs(outroot, exist_ok=True)
    out = os.path.join(outroot, "bridge_screen_results.csv")
    _atomic_csv(out, OUT_COLS, rows)
    fout = os.path.join(outroot, "bridge_screen_pp_full.csv")
    _atomic_csv(fout, list(full[0].keys()), full)
    stamp = {"arm": arm, "built_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "seconds": round(time.time() - t0, 1),
             "builder": rel(os.path.abspath(__file__)), "builder_md5": md5_file(os.path.abspath(__file__)),
             "vendored_md5": VENDOR_MD5, "worklist_md5": md5_file(WORKLIST), "chen_map_md5": md5_file(CHEN_MAP),
             "raw_manifest_md5": md5_file(RAW_MANIFEST) if os.path.exists(RAW_MANIFEST) else None,
             "prot_override": {k: rel(v) for k, v in (prot_override or {}).items()},
             "outputs": {os.path.basename(out): md5_file(out), os.path.basename(fout): md5_file(fout)},
             "n_rows": len(rows), "n_bridges": sum(r["support"] is True for r in rows)}
    with open(os.path.join(outroot, "BUILD_STAMP.json.part"), "w", encoding="utf-8") as fh:
        json.dump(stamp, fh, indent=1)
    os.replace(os.path.join(outroot, "BUILD_STAMP.json.part"), os.path.join(outroot, "BUILD_STAMP.json"))
    log("DONE arm=%s rows=%d bridges=%d -> %s (%.1fs)" % (arm, len(rows), stamp["n_bridges"], out, time.time() - t0))
    return rows, full


def _atomic_csv(path, cols, rows):
    tmp = path + ".part"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--arm", choices=["eur", "shipped"], default="eur")
    ap.add_argument("--outroot", required=True)
    a = ap.parse_args(argv)
    build(a.arm, os.path.abspath(a.outroot), log=lambda s: print(s, flush=True))


if __name__ == "__main__":
    main()
