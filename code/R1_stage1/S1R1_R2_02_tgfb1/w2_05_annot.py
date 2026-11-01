"""w2_05_annot.py -- annotation of rs1800470 and the TGFB1 neighbouring genes (pre-declaration section 2).

rs1800470: Ensembl REST GET /variation and /vep (JSON cached under ensembl_rest/; a cached response is reused).
Neighbouring genes: GRID_ens116_genes (the declared source) has NO TGFB1 locus (its loci are ARG1, BRAP, DPEP1, GSTA3,
SMPD3, UMOD) -> DEVIATION: genes from Ensembl REST /overlap/region over the cis window (GRCh38), release recorded.
r^2(rs1800470, rs12983775): PLINK19_paper6 --r2 on 503 (and 489 beside).
"""
import csv
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

REST = "https://rest.ensembl.org"
CACHE = W.OUT + "/ensembl_rest"
SIG = {"rs73045269": "A (discovery instrument)", "rs12983775": "B (CKB instrument)", "rs1800470": "annotated missense"}


def get(name, path):
    p = "%s/%s.json" % (CACHE, name)
    if os.path.exists(p):
        return json.load(open(p, encoding="utf-8")), "cache"
    import requests
    url = REST + path
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
    r.raise_for_status()
    d = {"url": url, "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "response": r.json()}
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=1)
    W.log("ensembl REST GET %s -> cached %s (md5 %s)" % (url, p, W.md5_file(p)))
    return d, "fetched"


def plink_r2(bfile, a, b, out):
    with open(out + ".snps", "w") as fh:
        fh.write("%s\n%s\n" % (a, b))
    plink = W.res("PLINK19_paper6", check=False)
    subprocess.run([plink, "--bfile", bfile, "--extract", out + ".snps", "--keep-allele-order", "--r", "square", "--out", out],
                   capture_output=True, text=True, check=True)
    import numpy as np
    R = np.loadtxt(out + ".ld")
    assert R.shape == (2, 2) and np.isfinite(R).all(), R
    return float(R[0, 1])


def main():
    os.makedirs(CACHE, exist_ok=True)
    W.log("==== w2_05_annot ====")
    out = []

    def add(qid, estimand, resource, value, notes="", panel="", unit=""):
        out.append({"quantity_id": qid, "prespec_ref": "s2", "estimand": estimand, "resource": resource, "panel": panel, "cov": "",
                    "se_order": "", "value": value, "se": "", "ci_lo": "", "ci_hi": "", "p": "", "n": "", "unit": unit, "notes": notes})
    try:
        info, _ = get("info_data", "/info/data")
        sw, _ = get("info_software", "/info/software")
        asm, _ = get("info_assembly", "/info/assembly/homo_sapiens")
        rel = "Ensembl release(s) %s; REST API %s; assembly %s" % (info["response"].get("releases"), sw["response"].get("release"),
                                                                 asm["response"].get("assembly_name"))
        var, _ = get("variation_rs1800470", "/variation/human/rs1800470")
        vep, _ = get("vep_rs1800470", "/vep/human/id/rs1800470?canonical=1")
        gene, _ = get("lookup_TGFB1", "/lookup/symbol/homo_sapiens/TGFB1")
    except Exception as e:  # network unavailable
        W.log("Ensembl REST unavailable: %s -> S2.annot.* NOT RUN" % e)
        add("S2.annot.rs1800470", "rs1800470 annotation (Ensembl REST)", "Ensembl REST", "NOT RUN: Ensembl REST unreachable (%s)" % str(e)[:80])
        add("S2.annot.genes", "neighbouring genes", "Ensembl REST", "NOT RUN: Ensembl REST unreachable")
        W.write_tsv(W.OUT + "/results_annot.tsv", list(out[0].keys()), out)
        return
    W.log("Ensembl: " + rel)
    v = var["response"]
    maps = [m for m in v.get("mappings", []) if m.get("assembly_name") == "GRCh38"]
    m = maps[0] if maps else {}
    vr = vep["response"][0]
    tcs = [t for t in vr.get("transcript_consequences", []) if t.get("gene_symbol") == "TGFB1"]
    can = [t for t in tcs if t.get("canonical")] or tcs
    t = can[0] if can else {}
    ann = ("rs1800470 GRCh38 %s:%s alleles %s; ancestral %s; minor allele %s (global MAF %s); most severe consequence %s; "
           "TGFB1 %s %s: %s, codon %s, amino acids %s at protein position %s; Ensembl %s" % (
               m.get("seq_region_name"), m.get("start"), m.get("allele_string"), m.get("ancestral_allele"), v.get("minor_allele"), v.get("MAF"),
               vr.get("most_severe_consequence"), t.get("transcript_id"), "(canonical)" if t.get("canonical") else "",
               ",".join(t.get("consequence_terms", [])), t.get("codons"), t.get("amino_acids"), t.get("protein_start"), rel))
    W.log("rs1800470: " + ann)
    add("S2.annot.rs1800470", "rs1800470 position, alleles and VEP consequence in TGFB1", "Ensembl REST (cached JSON ensembl_rest/)", ann,
        notes="variation + vep JSON md5 %s / %s" % (W.md5_file(CACHE + "/variation_rs1800470.json"), W.md5_file(CACHE + "/vep_rs1800470.json")))
    g = gene["response"]
    W.log("TGFB1 gene (Ensembl lookup): %s %s:%s-%s strand %s; instrument table gene_start/end 41288203-41353961" % (
        g.get("id"), g.get("seq_region_name"), g.get("start"), g.get("end"), g.get("strand")))
    # neighbouring genes over the cis window (GRCh38)
    ins = [r for r in csv.DictReader(open(W.res("EXP_UKBPPP_INSTRUMENT_TABLE"), encoding="utf-8")) if r["gene"] == "TGFB1"][0]
    lo, hi = int(ins["cis_start"]), int(ins["cis_end"])
    # the 1.07 Mb query returned HTTP 500 from the server; query 4 contiguous chunks (retry once each), de-duplicate by id
    step = (hi - lo) // 4 + 1
    gd = {}
    for k in range(4):
        a, b = lo + k * step, min(hi, lo + (k + 1) * step - 1)
        for attempt in range(3):
            try:
                ov, _ = get("overlap_TGFB1_cis_window_part%d" % (k + 1), "/overlap/region/human/19:%d-%d?feature=gene" % (a, b))
                break
            except Exception as e:
                W.log("overlap chunk %d attempt %d failed: %s" % (k + 1, attempt + 1, str(e)[:100]))
                time.sleep(5)
        else:
            raise SystemExit("Ensembl overlap failed for chunk %d" % (k + 1))
        for x in ov["response"]:
            gd[x["id"]] = x
    genes = sorted(gd.values(), key=lambda x: x["start"])
    pc = [x for x in genes if x.get("biotype") == "protein_coding"]
    W.log("genes overlapping the cis window 19:%d-%d (GRCh38): %d (%d protein-coding)" % (lo, hi, len(genes), len(pc)))
    pos38 = {}
    for r in W.read_tsv(W.OUT + "/tgfb1_locus.tsv"):
        if r["rsid503"] in SIG:
            pos38[r["rsid503"]] = int(r["pos38"])
    near = []
    for rs, p in sorted(pos38.items(), key=lambda kv: kv[1]):
        inside = [x.get("external_name") or x["id"] for x in pc if x["start"] <= p <= x["end"]]
        dist = sorted(((min(abs(p - x["start"]), abs(p - x["end"])) if not (x["start"] <= p <= x["end"]) else 0, x.get("external_name") or x["id"]) for x in pc))
        w50 = [nm for d, nm in dist if d <= 50000]
        near.append("%s (%s, GRCh38 19:%d): inside %s; protein-coding genes within 50 kb: %s; nearest %s" % (
            rs, SIG[rs], p, ",".join(inside) or "none", ",".join(w50) or "none", dist[0][1] if dist else "NA"))
    genes_txt = "protein-coding genes in the cis window (GRCh38 19:%d-%d, %d): %s || %s" % (
        lo, hi, len(pc), ", ".join("%s(%d-%d)" % (x.get("external_name") or x["id"], x["start"], x["end"]) for x in pc), " || ".join(near))
    for s in near:
        W.log("  " + s)
    add("S2.annot.genes", "neighbouring genes of the TGFB1 signals", "Ensembl REST /overlap/region (cached JSON)", genes_txt,
        notes="DEVIATION: GRID_ens116_genes has no TGFB1 locus; genes from Ensembl REST (%s); all biotypes in the cache file" % rel)
    # LD r^2 rs1800470 - rs12983775
    b503 = W.res("LD503_bed").replace(".bed", "")
    b489 = W.res("LD489_chr19_bed").replace(".bed", "")
    r503 = plink_r2(b503, "rs12983775", "rs1800470", W.OUT + "/ld/annot503")
    r489 = plink_r2(b489, "rs12983775", "rs1800470", W.OUT + "/ld/annot489")
    add("S2.ld.rs1800470_rs12983775.r2", "r^2(rs1800470, rs12983775)", "LD503", r503 ** 2, panel="503", unit="r2",
        notes="PLINK19_paper6 --r square --keep-allele-order; signed r on A1 (rs12983775:G, rs1800470:G) = %.6f" % r503)
    add("S2.p489.ld.rs1800470_rs12983775.r2", "r^2(rs1800470, rs12983775)", "LD489", r489 ** 2, panel="489", unit="r2",
        notes="PLINK19_paper6; signed r on A1 = %.6f" % r489)
    W.log("r2(rs1800470, rs12983775): 503 %.4f (r %.4f); 489 %.4f" % (r503 ** 2, r503, r489 ** 2))
    W.write_tsv(W.OUT + "/results_annot.tsv", list(out[0].keys()), out)


if __name__ == "__main__":
    main()
