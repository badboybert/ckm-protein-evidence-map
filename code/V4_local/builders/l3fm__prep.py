#!/usr/bin/env python
"""l3fm__prep.py -- harmonise every trait at each fine-mapping locus onto ONE SNP set matched to the external 1000G EUR
reference, and extract that reference's genotypes for the LD matrix. OFFLINE; reads only; writes under raw/l3-finemap
(reference extracts, the Wuttke region extract, RAW_MANIFEST.tsv) and staging/l3-finemap/<locus>/ (harmonised tables).

  python -B V3/V4_local/builders/l3fm__prep.py [--locus BRAP_12q24 ...] [--force]

Per locus (resumable: a locus whose PREP_STAMP.json input md5s still match is skipped):
  traits  = the bridge protein (D1 window, md5 asserted against the D1 RAW_MANIFEST), every LOCAL neighbour protein from the
            derived Olink inventory (md5 registered here), every D1 worklist metabolite of the bridge protein (D1 window).
  join    = (chr, GRCh38 pos) + compatible alleles (same / swapped; complement only when not palindromic); an ambiguous
            multi-allelic match is dropped, never guessed.
  ref     = GRCh37 position from the UKB-PPP ID + compatible alleles against g1000_eur.bim; A1 = counted allele.
  QC      = ref MAF >= 0.01; |aligned AF - ref AF| <= 0.20 for the protein AND the metabolite; palindromic kept only at
            ref MAF <= 0.40; monomorphic dropped. Every drop is counted in the stamp.
  control = UMOD_eGFR_control: UMOD (OID20237, European discovery) x Wuttke 2019 EA eGFR in the V4 engine's logged window.
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from collections import defaultdict

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import l3fm__common as C  # noqa: E402

B = C.B


def manifest_rows():
    p = os.path.join(C.RAW_FM, "RAW_MANIFEST.tsv")
    if not os.path.exists(p):
        return {}
    return {r["path"]: r for r in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")}


def register(entries):
    """add/refresh rows in raw/l3-finemap/RAW_MANIFEST.tsv (path, bytes, md5, kind, note)."""
    m = manifest_rows()
    for e in entries:
        m[e["path"]] = e
    os.makedirs(C.RAW_FM, exist_ok=True)
    C.write_tsv(os.path.join(C.RAW_FM, "RAW_MANIFEST.tsv"), ["path", "bytes", "md5", "kind", "note"],
                [m[k] for k in sorted(m)])


def entry(path, kind, note=""):
    return {"path": C.rel(path), "bytes": os.path.getsize(path), "md5": C.md5_file(path), "kind": kind, "note": note}


def index_by_pos(d):
    out = defaultdict(list)
    for k, v in d.items():
        out[(k[0], k[1])].append(v)
    return out


def match_one(p, cands):
    """the unique candidate whose alleles are compatible with p (+1 same orientation, -1 swapped)."""
    hits = []
    for c in cands:
        s = C.align_to_ref(c["ea"], c["oa"], p["ea"], p["oa"])
        if s:
            hits.append((c, s))
    return hits[0] if len(hits) == 1 else (None, 0 if not hits else "ambiguous")


def build_locus(name, spec, bimrows, log):
    """spec: dict(chr, lo, hi, traits=[dict(tid, kind, label, path, loader, md5, extra)]). Returns stamp dict."""
    ch = spec["chr"]
    primary = spec["traits"][0]
    assert primary["kind"] == "protein"
    loaded = {}
    drops = defaultdict(int)
    for t in spec["traits"]:
        if t["kind"] == "protein":
            d, ndup = C.read_protein_full(t["path"], ch, spec["lo"], spec["hi"])
        elif t["kind"] == "metabolite":
            d, ndup = C.read_metab_full(t["path"])
        elif t["kind"] == "outcome_gwas":
            d, ndup = t["preloaded"], t.get("ndup", 0)
        else:
            C.refuse("unknown trait kind %s" % t["kind"])
        if not d:
            C.refuse("%s/%s loaded 0 rows" % (name, t["tid"]))
        loaded[t["tid"]] = d
        drops["dup_key_" + t["tid"]] = ndup
    # join every trait to the primary protein
    pos_key = "pos37" if spec.get("join_on") == "pos37" else "pos38"
    idx = {}
    for t in spec["traits"][1:]:
        if pos_key == "pos37":
            ix = defaultdict(list)
            for k, v in loaded[t["tid"]].items():
                ix[(k[0], v["pos37"])].append(v)
            idx[t["tid"]] = ix
        else:
            idx[t["tid"]] = index_by_pos(loaded[t["tid"]])
    refidx = defaultdict(list)
    for r in bimrows:
        refidx[r["pos37"]].append(r)
    rows = []
    for k, p in loaded[primary["tid"]].items():
        rec = {"pos38": p.get("pos38", ""), "pos37": p["pos37"], "traits": {primary["tid"]: (p, 1)}}
        ok = True
        for t in spec["traits"][1:]:
            c, s = match_one(p, idx[t["tid"]].get((ch, p[pos_key]), []))
            if not c:
                drops["nojoin_" + t["tid"] if s == 0 else "ambiguous_" + t["tid"]] += 1
                ok = False
                break
            rec["traits"][t["tid"]] = (c, s)
        if not ok:
            continue
        # reference match
        hits = [(r, C.align_to_ref(p["ea"], p["oa"], r["a1"], r["a2"])) for r in refidx.get(p["pos37"], [])]
        hits = [h for h in hits if h[1]]
        if len(hits) != 1:
            drops["ref_absent" if not hits else "ref_ambiguous"] += 1
            continue
        r, sref = hits[0]
        rec["ref"], rec["sref"] = r, sref
        rows.append(rec)
    if not rows:
        C.refuse("%s: no SNP survives the join" % name)
    # genotypes for QC + LD
    rows.sort(key=lambda x: x["ref"]["idx"])
    G = C.read_bed_dosage([x["ref"]["idx"] for x in rows])
    if np.isnan(G).any():
        C.refuse("%s: missing genotypes in the reference extract" % name)
    af = G.mean(axis=0) / 2.0                           # A1 frequency in 503 EUR
    sd = G.std(axis=0)
    keep = []
    for j, x in enumerate(rows):
        maf = min(af[j], 1 - af[j])
        pal = C.is_palindromic(x["ref"]["a1"], x["ref"]["a2"])
        if sd[j] == 0:
            drops["monomorphic_ref"] += 1; continue
        if maf < C.REF_MAF_MIN:
            drops["ref_maf_lt_0.01"] += 1; continue
        if pal and maf > C.PALIN_MAF_MAX:
            drops["palindromic_maf_gt_0.40"] += 1; continue
        bad = False
        for tid, (v, s) in x["traits"].items():
            f = v.get("eaf", float("nan"))
            if f != f:
                continue
            # effect allele of trait tid expressed on the primary protein's ea, then on the reference A1
            f_on_pea = f if s == 1 else 1 - f
            f_a1 = f_on_pea if x["sref"] == 1 else 1 - f_on_pea
            if abs(f_a1 - af[j]) > C.AF_DIFF_MAX:
                drops["af_discordant_" + tid] += 1
                bad = True
                break
        if bad:
            continue
        keep.append(j)
    rows = [rows[j] for j in keep]
    G = G[:, keep]
    af = af[keep]
    # harmonised table: every beta on the reference counted allele A1
    tids = [t["tid"] for t in spec["traits"]]
    out = []
    for j, x in enumerate(rows):
        r = x["ref"]
        o = {"snp": "%s:%d:%s:%s" % (ch, r["pos37"], r["a1"], r["a2"]), "rsid_ref": r["rsid"], "bim_idx": r["idx"],
             "chr": ch, "pos37": r["pos37"], "pos38": x["pos38"], "a1": r["a1"], "a2": r["a2"], "af_ref": "%.6f" % af[j],
             "palindromic": int(C.is_palindromic(r["a1"], r["a2"]))}
        rsid_trait = ""
        for tid in tids:
            v, s = x["traits"][tid]
            sign = s * x["sref"]
            o["b__" + tid] = repr(sign * v["b"])
            o["se__" + tid] = repr(v["se"])
            o["z__" + tid] = repr(sign * v["b"] / v["se"])
            if v.get("rsid"):
                rsid_trait = v["rsid"]
        o["rsid_trait"] = rsid_trait
        out.append(o)
    # rsID concordance (diagnostic: the join is positional; rsIDs are an independent key)
    both = [o for o in out if o["rsid_trait"].startswith("rs") and o["rsid_ref"].startswith("rs")]
    rs_agree = sum(o["rsid_trait"] == o["rsid_ref"] for o in both)
    # sample sizes
    ninfo = {}
    for t in spec["traits"]:
        if t.get("n_fixed"):
            ninfo[t["tid"]] = {"n": t["n_fixed"], "n_source": t["n_source"]}
            continue
        vals = []
        for x in rows:
            v, s = x["traits"][t["tid"]]
            vals.append((v["b"], v["se"], v.get("eaf", float("nan"))))
        nh = C.n_hat(vals)
        ninfo[t["tid"]] = {"n": round(nh), "n_source": "n_hat = median 1/(2f(1-f)se^2), MAF>0.05, %d SNPs" % len(rows)}
        if t.get("n_meta"):
            ninfo[t["tid"]]["n_meta"] = t["n_meta"]
    odir = os.path.join(C.STAGE_FM, name)
    os.makedirs(odir, exist_ok=True)
    cols = ["snp", "rsid_ref", "rsid_trait", "bim_idx", "chr", "pos37", "pos38", "a1", "a2", "af_ref", "palindromic"]
    for tid in tids:
        cols += ["b__" + tid, "se__" + tid, "z__" + tid]
    C.write_tsv(os.path.join(odir, "harmonised.tsv.gz"), cols, out, gz=True)
    gpath = os.path.join(C.RAW_FM, "%s__ref_A1dosage_503xP.i8" % name)
    tmp = gpath + ".part"
    G.astype(np.int8).tofile(tmp)
    os.replace(tmp, gpath)
    traits_rows = []
    for t in spec["traits"]:
        traits_rows.append({"tid": t["tid"], "kind": t["kind"], "label": t["label"], "gene": t.get("gene", ""),
                            "oid": t.get("oid", ""), "gcst": t.get("gcst", ""), "role": t.get("role", ""),
                            "source": C.rel(t["path"]), "source_md5": t["md5"], "n": ninfo[t["tid"]]["n"],
                            "n_source": ninfo[t["tid"]]["n_source"], "n_meta": ninfo[t["tid"]].get("n_meta", ""),
                            "d1_H4": t.get("d1_H4", ""), "d1_H3": t.get("d1_H3", ""), "d1_nsnp": t.get("d1_nsnp", ""),
                            "d1_support": t.get("d1_support", "")})
    C.write_tsv(os.path.join(odir, "traits.tsv"), list(traits_rows[0].keys()), traits_rows)
    stamp = {"locus": name, "chr": ch, "window_grch38": [spec["lo"], spec["hi"]], "join_on": pos_key,
             "window_grch37_ref": [spec["lo37"], spec["hi37"]], "n_snps": len(out), "n_primary_rows": len(loaded[primary["tid"]]),
             "drops": dict(drops), "rsid_concordance": "%d/%d" % (rs_agree, len(both)),
             "ref_panel": C.rel(C.REF_PREFIX), "ref_n": int(G.shape[0]), "ref_maf_min": C.REF_MAF_MIN,
             "af_diff_max": C.AF_DIFF_MAX, "palin_maf_max": C.PALIN_MAF_MAX,
             "inputs": {t["tid"]: {"path": C.rel(t["path"]), "md5": t["md5"]} for t in spec["traits"]},
             "outputs": {"harmonised.tsv.gz": C.md5_file(os.path.join(odir, "harmonised.tsv.gz")),
                         "traits.tsv": C.md5_file(os.path.join(odir, "traits.tsv")),
                         C.rel(gpath): C.md5_file(gpath)},
             "builder_md5": C.md5_file(os.path.abspath(__file__)), "common_md5": C.md5_file(os.path.join(HERE, "l3fm__common.py")),
             "built_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    C.write_json(os.path.join(odir, "PREP_STAMP.json"), stamp)
    register([entry(gpath, "ref_extract", "g1000_eur A1 dosage, int8, 503 x %d, SNP order = %s/harmonised.tsv.gz" % (len(out), name))])
    log("[prep] %-18s snps=%d  rsid %s  drops=%s" % (name, len(out), stamp["rsid_concordance"], dict(drops)))
    return stamp


def locus_specs(which):
    loci = C.locus_units()
    inv = C.neighbour_inventory(loci)
    os.makedirs(C.STAGE_FM, exist_ok=True)
    C.write_tsv(os.path.join(C.STAGE_FM, "neighbour_inventory.tsv"),
                ["locus", "gene", "oid", "gene_start", "gene_end", "named_by_user", "on_olink", "status", "local_window",
                 "synapse_id", "synapse_tar"], inv)
    d1 = {(r["protein"], r["metabolite"]): r for r in csv.DictReader(open(os.path.join(C.D1_STAGE, "bridge_screen_results.csv"), encoding="utf-8"))}
    d1man = B.read_manifest()
    specs = {}
    reg = []
    for loc, L in loci.items():
        if which and loc not in which:
            continue
        ppath = B.eur_window_path(L["bridge"], L["oid"], L["chr"])
        meta = ppath.replace(".tsv", ".meta.json")
        nmeta = json.load(open(meta))["N_max"] if os.path.exists(meta) else ""
        traits = [{"tid": "P__" + L["bridge"], "kind": "protein", "label": L["bridge"], "gene": L["bridge"], "oid": L["oid"],
                   "path": ppath, "md5": B.checked_md5(ppath, d1man), "role": "bridge_protein", "n_meta": nmeta}]
        for r in inv:
            if r["locus"] == loc and r["status"].startswith("LOCAL"):
                p = os.path.join(C.ROOT, r["local_window"])
                if os.path.abspath(p).startswith(os.path.abspath(C.NEIGHBOUR_WIN_DIR)):
                    # extracted here by l3fm__extract_neighbours.py: already registered -> md5 must match, never re-registered
                    md5 = B.checked_md5(p, manifest_rows())
                else:
                    md5 = C.md5_file(p)
                    reg.append(entry(p, "input_readonly_protein_window", "neighbour %s at %s (European discovery)" % (r["gene"], loc)))
                traits.append({"tid": "P__" + r["gene"], "kind": "protein", "label": r["gene"], "gene": r["gene"], "oid": r["oid"],
                               "path": p, "md5": md5, "role": "neighbour_protein"})
        for u in L["units"]:
            mp = B.metab_window_path(u)
            dr = d1[(u["protein"], u["metabolite"])]
            traits.append({"tid": "M__" + u["gcst"], "kind": "metabolite", "label": u["metabolite"], "gcst": u["gcst"],
                           "path": mp, "md5": B.checked_md5(mp, d1man), "role": "metabolite",
                           "d1_H4": dr["H4"], "d1_H3": dr["H3"], "d1_nsnp": dr["nsnp"], "d1_support": dr["support"]})
        specs[loc] = {"chr": L["chr"], "lo": L["lo"], "hi": L["hi"], "traits": traits}
    if not which or "UMOD_eGFR_control" in which:
        ch, lo, hi = C.engine_umod_window()
        om = C.olink_map()["UMOD"][0]
        ppath = C.local_protein_window(om["oid"], ch)
        if not ppath:
            C.refuse("UMOD European window not local")
        reg.append(entry(ppath, "input_readonly_protein_window", "UMOD positive control (European discovery)"))
        specs["UMOD_eGFR_control"] = {"chr": ch, "lo": lo, "hi": hi, "join_on": "pos37",
                                      "traits": [{"tid": "P__UMOD", "kind": "protein", "label": "UMOD", "gene": "UMOD", "oid": om["oid"],
                                                  "path": ppath, "md5": C.md5_file(ppath), "role": "control_protein"},
                                                 {"tid": "O__eGFR", "kind": "outcome_gwas", "label": "eGFR (Wuttke 2019 EA)",
                                                  "path": C.WUTTKE_EA, "md5": "", "role": "control_outcome"}]}
    if reg:
        register(reg)
    return specs


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", nargs="*")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    t0 = time.time()
    if C.fam_n() != C.REF_N:
        C.refuse("reference .fam has %d samples, expected %d" % (C.fam_n(), C.REF_N))
    specs = locus_specs(a.locus)
    # GRCh37 extent of each locus = the primary protein rows inside the GRCh38 window (from their IDs)
    todo = {}
    for name, s in specs.items():
        stp = os.path.join(C.STAGE_FM, name, "PREP_STAMP.json")
        p = s["traits"][0]
        d, _ = C.read_protein_full(p["path"], s["chr"], s["lo"], s["hi"])
        s["lo37"], s["hi37"] = min(v["pos37"] for v in d.values()) - 1000, max(v["pos37"] for v in d.values()) + 1000
        if s.get("join_on") == "pos37":
            wout = os.path.join(C.RAW_FM, "wuttke_EA_eGFR__chr%s_%d_%d.tsv.gz" % (s["chr"], s["lo37"], s["hi37"]))
            if not os.path.exists(wout):
                wd, wdup, wn = C.read_wuttke_region(C.WUTTKE_EA, s["chr"], s["lo37"], s["hi37"])
                cols = ["pos37", "rsid", "ea", "oa", "eaf", "b", "se", "n"]
                C.write_tsv(wout, cols, sorted(wd.values(), key=lambda v: v["pos37"]), gz=True)
                register([entry(wout, "outcome_region_extract", "Wuttke 2019 EA eGFR, %d rows read, %d dup keys dropped; source %s (%d B)"
                                % (wn, wdup, C.WUTTKE_EA, os.path.getsize(C.WUTTKE_EA)))])
            wd = {}
            with __import__("gzip").open(wout, "rt", encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    v = {"pos37": int(r["pos37"]), "rsid": r["rsid"], "ea": r["ea"], "oa": r["oa"], "eaf": float(r["eaf"]),
                         "b": float(r["b"]), "se": float(r["se"]), "n": float(r["n"])}
                    wd[(s["chr"], v["pos37"], frozenset((v["ea"], v["oa"])))] = v
            s["traits"][1]["preloaded"] = wd
            s["traits"][1]["path"] = wout
            s["traits"][1]["md5"] = C.md5_file(wout)
            s["traits"][1]["n_fixed"] = int(np.median([v["n"] for v in wd.values()]))
            s["traits"][1]["n_source"] = "median n_total_sum over the extracted region (engine rule: regional median)"
        if os.path.exists(stp) and not a.force:
            old = json.load(open(stp))
            if all(old["inputs"].get(t["tid"], {}).get("md5") == t["md5"] for t in s["traits"]) and \
                    old.get("builder_md5") == C.md5_file(os.path.abspath(__file__)) and \
                    old.get("common_md5") == C.md5_file(os.path.join(HERE, "l3fm__common.py")):
                print("[prep] %-18s up to date (skip)" % name, flush=True)
                continue
        todo[name] = s
    if todo:
        regions = {(s["chr"], s["lo37"], s["hi37"]): n for n, s in todo.items()}
        t1 = time.time()
        bim = C.ref_bim_regions(list(regions))
        print("[prep] bim scan %.1fs: %s" % (time.time() - t1, {regions[r]: len(v) for r, v in bim.items()}), flush=True)
        for r, name in regions.items():
            build_locus(name, todo[name], bim[r], lambda m: print(m, flush=True))
    print("[prep] done in %.1fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
