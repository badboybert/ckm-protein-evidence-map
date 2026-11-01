"""w2_07_twobytwo.py -- the descriptive 2 x 2 table (pre-declaration section 2), regenerated from raw files by one script.

Cells (per EUROPEAN protein-raising allele; SNP-outcome association on the outcome's own scale):
  EURinstr_EURout, CKBinstr_EURout, CKBinstr_EASout, EURinstr_EASout.
Sentinels: TGFB1, FGFR4, CELSR2, IL6R, ABO, GALNT3, SPON1, FGF5, ARG1 signal A.
EUR instrument = the UKB-PPP cis instrument (instrument table; ARG1 = signal A rs2781668; TGFB1 = rs73045269); its
  protein-raising allele from the UKB-PPP cis file BETA (asserted equal to the instrument table).
CKB instrument = the CKB pQTL cis lead: the variant with the smallest printed p among p < 5e-8 inside the gene +/- 500 kb
  window (instrument-table cis_start..cis_end, GRCh38); ties (e.g. p floored at 1e-300) broken by the larger |beta/sebeta|,
  then position (tie counts logged). Beta is printed to 2 significant figures, so |beta/se| alone is too coarse to rank.
  TGFB1 must be rs12983775 (GRCh38 19:41351534; pre-declared) -- asserted. Its EUROPEAN protein-raising allele comes from
  the UKB-PPP cis file at that variant; if the variant is absent there, the CKB-raising allele is used and flagged.
  The instrument is held fixed across outcomes (no proxies; an absent variant is reported as absent).
EUR outcome: FinnGen R10 (CHD for TGFB1/CELSR2/IL6R; T2D for ABO/GALNT3/ARG1/FGFR4; AF for SPON1), CKDGen 2019 EA eGFR
  for FGF5, DIAMANTE 2018 noUKB EUR as FGFR4's second EUR cell.
EAS outcome: BBJ CAD (TGFB1, CELSR2, IL6R); AGEN T2D (ABO, GALNT3, ARG1, FGFR4); BBJ AF (SPON1); TWB2+BBJ eGFR (FGF5);
  BBJ T2D as FGFR4's second EAS cell. FGFR4's AGEN cell is re-cut from the raw AGEN file (not the pre-lifted window).
Joins: GRCh38 files on pos38 + alleles; GRCh37 files on the GRCh37 position from the UKB-PPP variant ID (or, for a CKB
  variant absent from UKB-PPP, its hg38 -> hg19 liftover with REF_LIFTOVER_hg38ToHg19) + alleles; complement only for
  non-palindromic variants; palindromic variants matched on exact labels and checked against the population frequency
  (EUR outcomes vs UKB-PPP A1FREQ, EAS outcomes vs CKB af): |diff| > 0.10 or MAF > 0.42 -> orientation not verified.
census_value = EAS_CENSUS_RESULTS 'by' (CKB-instrument x EAS-outcome cell only) re-oriented to the same allele, printed as
  a cross-check; it is a pipeline output and never an input to a cell.
"""
import csv
import gzip
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

EXT = W.OUT + "/extracts/"
CKB_DIR = W.ROOT + "/V3/V4_local/raw/d2-eas-census/ckb_regions/"
D2_MANIFEST = (W.ROOT + "/V3/V4_local/staging/d2-eas-census/D2_MANIFEST.tsv", "99ebeff1d7df613c62dcedac2264fb1a")

# resource id -> (label, build, parser)
RES = {
    "OUT_FG_R10_I9_CHD": ("FinnGen R10 I9_CHD", "38", "finngen"),
    "OUT_FG_R10_T2D": ("FinnGen R10 T2D", "38", "finngen"),
    "OUT_FG_R10_I9_AF": ("FinnGen R10 I9_AF", "38", "finngen"),
    "OUT_DIAMANTE2018_noUKB_EUR_chrpos": ("DIAMANTE 2018 EUR noUKB T2D", "37", "diamante"),
    "OUT_CKDGEN2019_eGFR_EA": ("CKDGen 2019 EA eGFR (ln eGFR)", "37", "ckdgen"),
    "EAS_BBJ_CAD_hum0197v5": ("BBJ CAD (hum0197.v5)", "37", "saige"),
    "EAS_AGEN_T2D_Spracklen2020": ("AGEN T2D Spracklen 2020 (includes CKB)", "37", "agen"),
    "EAS_BBJ_AF_2022": ("BBJ AF 2022", "37", "bbjaf"),
    "EAS_TWB2_BBJ_eGFR": ("TWB2+BBJ eGFR (per SD RINT)", "37", "twb2"),
    "EAS_BBJ_T2D_hum0197v3": ("BBJ T2D (hum0197.v3)", "37", "saige"),
}
LOCI = {
    "TGFB1": {"exp": "EXP_UKBPPP_EUR_TGFB1_OID20621", "eur": ["OUT_FG_R10_I9_CHD"], "eas": ["EAS_BBJ_CAD_hum0197v5"], "census": "CAD", "ckb": "TGFB1__chr19_40751534_41951534.tsv"},
    "FGFR4": {"exp": "EXP_UKBPPP_EUR_FGFR4_OID30646", "eur": ["OUT_FG_R10_T2D", "OUT_DIAMANTE2018_noUKB_EUR_chrpos"],
              "eas": ["EAS_AGEN_T2D_Spracklen2020", "EAS_BBJ_T2D_hum0197v3"], "census": None, "ckb": "CKB_FGFR4_cis_eas"},
    "CELSR2": {"exp": "EXP_UKBPPP_EUR_CELSR2_OID30593", "eur": ["OUT_FG_R10_I9_CHD"], "eas": ["EAS_BBJ_CAD_hum0197v5"], "census": "CAD", "ckb": "CELSR2__chr1_108679521_109879521.tsv"},
    "IL6R": {"exp": "EXP_UKBPPP_EUR_IL6R_OID20385", "eur": ["OUT_FG_R10_I9_CHD"], "eas": ["EAS_BBJ_CAD_hum0197v5"], "census": "CAD", "ckb": "IL6R__chr1_153844459_155044459.tsv"},
    "ABO": {"exp": "EXP_UKBPPP_EUR_ABO_OID30675", "eur": ["OUT_FG_R10_T2D"], "eas": ["EAS_AGEN_T2D_Spracklen2020"], "census": "T2D", "ckb": "ABO__chr9_132657521_133857521.tsv"},
    "GALNT3": {"exp": "EXP_UKBPPP_EUR_GALNT3_OID20471", "eur": ["OUT_FG_R10_T2D"], "eas": ["EAS_AGEN_T2D_Spracklen2020"], "census": "T2D", "ckb": "GALNT3__chr2_165125521_166346201.tsv"},
    "SPON1": {"exp": "EXP_UKBPPP_EUR_SPON1_OID20759", "eur": ["OUT_FG_R10_I9_AF"], "eas": ["EAS_BBJ_AF_2022"], "census": "AF", "ckb": "SPON1__chr11_13417074_14768133.tsv"},
    "FGF5": {"exp": "EXP_UKBPPP_EUR_FGF5_OID20490", "eur": ["OUT_CKDGEN2019_eGFR_EA"], "eas": ["EAS_TWB2_BBJ_eGFR"], "census": "eGFR", "ckb": "FGF5__chr4_79663187_80863187.tsv"},
    "ARG1": {"exp": "EXP_UKBPPP_EUR_ARG1_OID21445", "eur": ["OUT_FG_R10_T2D"], "eas": ["EAS_AGEN_T2D_Spracklen2020"], "census": "T2D", "ckb": "ARG1__chr6_130970832_132173339.tsv"},
}
TAG = {"OUT_FG_R10_T2D": "FinnGen", "OUT_DIAMANTE2018_noUKB_EUR_chrpos": "DIAMANTE", "EAS_AGEN_T2D_Spracklen2020": "AGEN", "EAS_BBJ_T2D_hum0197v3": "BBJ"}
PLACEHOLDERS = {"", "NA", "na", "NaN", "nan", ".", "-", "NULL", "#N/A", "TBD", "Inf", "-Inf"}


def fnum(x):
    if x in PLACEHOLDERS:
        return float("nan")
    return float(x)


def parse_row(kind, f):
    """-> dict(pos, ea, oa, b, se, p, eaf) with effect allele ea; alleles upper-cased."""
    if kind == "finngen":
        return dict(pos=int(f["pos"]), ea=f["alt"], oa=f["ref"], b=fnum(f["beta"]), se=fnum(f["sebeta"]), p=fnum(f["pval"]), eaf=fnum(f["af_alt"]), id=f["rsids"])
    if kind == "diamante":
        return dict(pos=int(f["Pos"]), ea=f["EA"].upper(), oa=f["NEA"].upper(), b=fnum(f["Beta"]), se=fnum(f["SE"]), p=fnum(f["Pvalue"]), eaf=fnum(f["EAF"]), id=f["SNP"])
    if kind == "ckdgen":
        return dict(pos=int(f["Pos_b37"]), ea=f["Allele1"].upper(), oa=f["Allele2"].upper(), b=fnum(f["Effect"]), se=fnum(f["StdErr"]), p=fnum(f["P-value"]), eaf=fnum(f["Freq1"]), id=f["RSID"])
    if kind == "saige":
        pos = f["POS"]
        pos = int(pos) if pos.isdigit() else int(float(pos))
        return dict(pos=pos, ea=f["Allele2"].upper(), oa=f["Allele1"].upper(), b=fnum(f["BETA"]), se=fnum(f["SE"]), p=fnum(f["p.value"]), eaf=fnum(f["AF_Allele2"]), id=f["SNPID"])
    if kind == "agen":
        pos = f["Pos"]
        pos = int(pos) if pos.isdigit() else int(f["Chrpos"].split(":")[1]) if f["Chrpos"].split(":")[1].isdigit() else int(float(pos))
        return dict(pos=pos, ea=f["EA"].upper(), oa=f["NEA"].upper(), b=fnum(f["Beta"]), se=fnum(f["SE"]), p=fnum(f["P"]), eaf=fnum(f["EAF"]), id=f["MarkerName"])
    if kind == "bbjaf":
        return dict(pos=int(f["POS"]), ea=f["ALT"].upper(), oa=f["REF"].upper(), b=fnum(f["BETA"]), se=fnum(f["SE"]), p=fnum(f["P"]), eaf=fnum(f["AAF"]), id="")
    if kind == "twb2":
        return dict(pos=int(f["BP"]), ea=f["Allele1"].upper(), oa=f["Allele2"].upper(), b=fnum(f["Effect"]), se=fnum(f["StdErr"]), p=fnum(f["P"]), eaf=float("nan"), id=f["SNP"])
    raise ValueError(kind)


def load_extract(rid, locus):
    st = json.load(open(EXT + rid + ".stamp.json", encoding="utf-8"))
    assert st["source_md5"] == W.res_entry(rid)["md5"], rid
    fn = "%s__%s.tsv" % (rid, locus)
    W.assert_md5(EXT + fn, st["file_md5"][fn], "extract " + fn)
    kind = RES[rid][2]
    ws = kind in ("ckdgen", "twb2")
    rows, nulls, n = {}, {}, 0
    with open(EXT + fn, encoding="utf-8") as fh:
        hdr = fh.readline().split() if ws else fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.split() if ws else line.rstrip("\n").split("\t")))
            n += 1
            try:
                r = parse_row(kind, f)
            except ValueError:
                nulls["unparsable_pos"] = nulls.get("unparsable_pos", 0) + 1
                continue
            for k in ("ea", "oa"):
                if r[k] in PLACEHOLDERS:
                    nulls[k] = nulls.get(k, 0) + 1
            for k in ("b", "se", "p", "eaf"):
                if r[k] != r[k]:
                    nulls[k] = nulls.get(k, 0) + 1
            rows.setdefault(r["pos"], []).append(r)
    return rows, n, {k: round(v / n, 5) for k, v in nulls.items()}


def match(cands, a_eff, a_oth):
    """the unique candidate with allele set {a_eff, a_oth} (complement only when not palindromic) -> (row, sign, how)."""
    hits = []
    for c in cands:
        s = W.align(c["ea"], c["oa"], a_eff, a_oth)
        if s:
            hits.append((c, s))
    if len(hits) == 1:
        c, s = hits[0]
        how = "exact" if {c["ea"], c["oa"]} == {a_eff, a_oth} else "complement"
        return c, s, how
    return None, 0, ("absent" if not hits else "ambiguous (%d)" % len(hits))


def read_ukb(rid):
    out = {}
    with open(W.res(rid), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c, p37, ref, alt = r["ID"].split(":")[:4]
            out.setdefault(int(r["pos"]), []).append({"pos38": int(r["pos"]), "pos37": int(p37), "ea": alt, "oa": ref, "b": float(r["BETA"]),
                                                      "se": float(r["SE"]), "eaf": float(r["A1FREQ"]), "id": r["ID"]})
    return out


def read_ckb(locus, fname, manifest):
    if fname == "CKB_FGFR4_cis_eas":
        p = W.res("CKB_FGFR4_cis_eas")
    else:
        p = CKB_DIR + fname
        key = "V3/V4_local/raw/d2-eas-census/ckb_regions/" + fname
        W.assert_md5(p, manifest[key], "CKB region " + fname)
    rows = []
    with open(p, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({"pos38": int(r["position"]), "ref": r["ref"], "alt": r["alt"], "b": fnum(r["beta"]), "se": fnum(r["sebeta"]),
                         "p": fnum(r["pvalue"]), "af": fnum(r["af"]), "rsid": r["rsid"], "id": r["id"]})
    return rows, p


def main():
    W.log("==== w2_07_twobytwo ====")
    W.assert_md5(*D2_MANIFEST, "D2_MANIFEST.tsv (CKB region md5 source)")
    manifest = {r["path"]: r["md5"] for r in csv.DictReader(open(D2_MANIFEST[0], encoding="utf-8"), delimiter="\t")}
    ins = {r["gene"]: r for r in csv.DictReader(open(W.res("EXP_UKBPPP_INSTRUMENT_TABLE"), encoding="utf-8"))}
    census = {(r["gene"], r["outcome"]): r for r in csv.DictReader(open(W.res("EAS_CENSUS_RESULTS"), encoding="utf-8"), delimiter="\t")}
    from pyliftover import LiftOver
    lo = LiftOver(W.res("REF_LIFTOVER_hg38ToHg19"))
    out = []
    joinlog = []

    def add(qid, estimand, resource, value, se="", p="", notes="", unit=""):
        out.append({"quantity_id": qid, "prespec_ref": "s2 (2 x 2 table, descriptive)", "estimand": estimand, "resource": resource, "panel": "",
                    "cov": "", "se_order": "", "value": value, "se": se, "ci_lo": "", "ci_hi": "", "p": p, "n": "", "unit": unit, "notes": notes})
    ext_cache = {}
    for locus, L in LOCI.items():
        ir = ins[locus]
        lo38, hi38 = int(ir["cis_start"]), int(ir["cis_end"])
        ukb = read_ukb(L["exp"])
        # ---- EUR instrument
        pos37e, ref_e, alt_e = int(ir["ukb_variant_id"].split(":")[1]), ir["ukb_variant_id"].split(":")[2], ir["ukb_variant_id"].split(":")[3]
        cand = [u for u in ukb.get(int(ir["chrpos"].split(":")[1]), []) if u["id"].startswith(ir["ukb_variant_id"])]
        assert len(cand) == 1, (locus, ir["ukb_variant_id"])
        ue = cand[0]
        s_tab = 1 if (ir["effect_allele"], ir["other_allele"]) == (ue["ea"], ue["oa"]) else -1
        assert abs(s_tab * float(ir["beta"]) - ue["b"]) < 1e-5 * max(1, abs(ue["b"])), (locus, ir["beta"], ue["b"])
        if locus == "ARG1":
            assert ir["SNP"] == "rs2781668"
        if locus == "TGFB1":
            assert ir["SNP"] == "rs73045269"
        e_raise, e_other = (ue["ea"], ue["oa"]) if ue["b"] > 0 else (ue["oa"], ue["ea"])
        eur_instr = {"rsid": ir["SNP"], "pos38": ue["pos38"], "pos37": ue["pos37"], "raise": e_raise, "other": e_other,
                     "bx_ukb": abs(ue["b"]), "se_ukb": ue["se"], "f_raise_ukb": ue["eaf"] if e_raise == ue["ea"] else 1 - ue["eaf"],
                     "pal": W.is_pal(e_raise, e_other)}
        # ---- CKB instrument
        ckb, ckbp = read_ckb(locus, L["ckb"], manifest)
        inwin = [r for r in ckb if lo38 <= r["pos38"] <= hi38 and r["p"] == r["p"] and r["p"] < 5e-8 and r["se"] > 0]
        for r in inwin:
            r["absz"] = abs(r["b"] / r["se"])
        inwin.sort(key=lambda r: (r["p"], -r["absz"], r["pos38"]))
        top = inwin[0]
        ties_p = sum(1 for r in inwin if r["p"] == top["p"])
        ties = sum(1 for r in inwin if r["p"] == top["p"] and abs(r["absz"] - top["absz"]) < 1e-9)
        if locus == "TGFB1":
            assert top["pos38"] == 41351534, "TGFB1 CKB lead is not rs12983775 (GRCh38 41351534): %s" % top
        uk = [u for u in ukb.get(top["pos38"], []) if W.align(u["ea"], u["oa"], top["alt"], top["ref"])]
        if len(uk) == 1:
            u = uk[0]
            k_raise, k_other = (u["ea"], u["oa"]) if u["b"] > 0 else (u["oa"], u["ea"])
            pos37k, src37 = u["pos37"], "UKB-PPP variant ID"
            bx_ukb_k, se_ukb_k = abs(u["b"]), u["se"]
            f_raise_ukb_k = u["eaf"] if k_raise == u["ea"] else 1 - u["eaf"]
            orient = "EUROPEAN protein-raising allele (UKB-PPP BETA %+.4f per %s)" % (u["b"], u["ea"])
        else:
            k_raise, k_other = (top["alt"], top["ref"]) if top["b"] > 0 else (top["ref"], top["alt"])
            r37 = lo.convert_coordinate("chr" + ir["gene_chrom"], top["pos38"] - 1)
            pos37k, src37 = (r37[0][1] + 1 if r37 else None), "hg38->hg19 liftover"
            bx_ukb_k, se_ukb_k, f_raise_ukb_k = float("nan"), float("nan"), float("nan")
            orient = "FLAG: variant absent from the UKB-PPP cis file -> CKB protein-raising allele used"
        bx_ckb_k = top["b"] if k_raise == top["alt"] else -top["b"]
        f_raise_ckb = top["af"] if k_raise == top["alt"] else 1 - top["af"]
        ckb_instr = {"id": top["id"], "rsid": top["rsid"], "pos38": top["pos38"], "pos37": pos37k, "raise": k_raise, "other": k_other,
                     "bx_ckb": bx_ckb_k, "se_ckb": top["se"], "p_ckb": top["p"], "f_raise_ckb": f_raise_ckb, "bx_ukb": bx_ukb_k,
                     "se_ukb": se_ukb_k, "f_raise_ukb": f_raise_ukb_k, "pal": W.is_pal(k_raise, k_other)}
        # EUR protein beta and CKB protein beta at the EUR instrument (descriptive)
        ck_at_e = [r for r in ckb if r["pos38"] == eur_instr["pos38"] and W.align(r["alt"], r["ref"], e_raise, e_other)]
        bx_ckb_e = (W.align(ck_at_e[0]["alt"], ck_at_e[0]["ref"], e_raise, e_other) * ck_at_e[0]["b"]) if len(ck_at_e) == 1 else float("nan")
        se_ckb_e = ck_at_e[0]["se"] if len(ck_at_e) == 1 else float("nan")
        f_raise_ckb_e = (ck_at_e[0]["af"] if ck_at_e[0]["alt"] == e_raise else 1 - ck_at_e[0]["af"]) if len(ck_at_e) == 1 else float("nan")
        eur_instr.update({"bx_ckb": bx_ckb_e, "se_ckb": se_ckb_e, "f_raise_ckb": f_raise_ckb_e})
        add("S2.x22.%s.EURinstr.snp" % locus, "EUR (UKB-PPP) instrument and its European protein-raising allele", L["exp"],
            "%s GRCh38 %s:%d GRCh37 %d raising %s other %s" % (eur_instr["rsid"], ir["gene_chrom"], eur_instr["pos38"], eur_instr["pos37"], e_raise, e_other),
            notes="UKB-PPP beta per raising allele %.4f (se %.4f), raising-allele freq %.3f; CKB beta per same allele %s (se %s); palindromic %s" % (
                eur_instr["bx_ukb"], eur_instr["se_ukb"], eur_instr["f_raise_ukb"], "%.4f" % bx_ckb_e if bx_ckb_e == bx_ckb_e else "absent",
                "%.4f" % se_ckb_e if se_ckb_e == se_ckb_e else "NA", eur_instr["pal"]))
        add("S2.x22.%s.CKBinstr.snp" % locus, "CKB instrument (CKB cis lead) and its European protein-raising allele", ckbp,
            "%s %s GRCh38 %s:%d GRCh37 %s (%s) raising %s other %s" % (ckb_instr["id"], ckb_instr["rsid"] or "(no rsid)", ir["gene_chrom"], ckb_instr["pos38"],
                                                                     ckb_instr["pos37"], src37, k_raise, k_other),
            notes="%s; CKB beta per raising allele %.4f (se %.4f, p %.2g, raising-allele freq %.3f); UKB-PPP beta per same allele %s; "
                  "CKB lead rule: smallest p, then largest |beta/se|, among p < 5e-8 in %s:%d-%d (%d candidates; %d tied at the top p, %d also tied on |z|); palindromic %s" % (
                      orient, bx_ckb_k, top["se"], top["p"], f_raise_ckb, "%.4f" % bx_ukb_k if bx_ukb_k == bx_ukb_k else "absent",
                      ir["gene_chrom"], lo38, hi38, len(inwin), ties_p, ties, ckb_instr["pal"]))
        W.log("[2x2 %s] EUR instr %s raise %s (bx UKB %.3f, CKB %s); CKB instr %s %s raise %s (bx CKB %.3f, UKB %s); ties %d" % (
            locus, eur_instr["rsid"], e_raise, eur_instr["bx_ukb"], "%.3f" % bx_ckb_e if bx_ckb_e == bx_ckb_e else "absent", ckb_instr["id"],
            ckb_instr["rsid"], k_raise, bx_ckb_k, "%.3f" % bx_ukb_k if bx_ukb_k == bx_ukb_k else "absent", ties))
        # ---- cells
        multi = len(L["eur"]) > 1
        for pop, rids in (("EUR", L["eur"]), ("EAS", L["eas"])):
            for rid in rids:
                key = (rid, locus)
                if key not in ext_cache:
                    ext_cache[key] = load_extract(rid, locus)
                    rows, n, nulls = ext_cache[key]
                    W.log("[2x2 %s] extract %s: %d rows; per-key null/placeholder rates %s" % (locus, rid, n, nulls or "{} (none)"))
                    joinlog.append({"locus": locus, "resource": rid, "rows": n, "null_rates": json.dumps(nulls)})
                rows, n, nulls = ext_cache[key]
                label, build, _ = RES[rid]
                for iname, I in (("EURinstr", eur_instr), ("CKBinstr", ckb_instr)):
                    cell = "%s_%sout" % (iname, pop) + ("_" + TAG[rid] if multi else "")
                    pos = I["pos38"] if build == "38" else I["pos37"]
                    c, s, how = match(rows.get(pos, []) if pos else [], I["raise"], I["other"])
                    note = []
                    if c is None:
                        W.log("[2x2 %s] %s x %s: %s at %s:%s" % (locus, cell, rid, how, build, pos))
                        add("S2.x22.%s.%s.beta" % (locus, cell), "SNP-outcome beta per European protein-raising allele", label, "NA",
                            notes="instrument %s %s in %s (GRCh%s %s): %s; no proxy (instrument held fixed)" % (
                                I.get("rsid") or I.get("id"), I["raise"], label, build, pos, how))
                        for fld in ("se", "p"):
                            add("S2.x22.%s.%s.%s" % (locus, cell, fld), "", label, "NA", notes="variant %s" % how)
                        joinlog.append({"locus": locus, "resource": rid, "cell": cell, "result": how})
                    else:
                        b = s * c["b"]
                        f_out = float("nan") if c["eaf"] != c["eaf"] else (c["eaf"] if s == 1 else 1 - c["eaf"])
                        fref = I["f_raise_ukb"] if pop == "EUR" else I.get("f_raise_ckb", float("nan"))
                        if I["pal"]:
                            if f_out != f_out or fref != fref:
                                note.append("PALINDROMIC: orientation by allele labels only (no frequency to check)")
                            elif abs(f_out - fref) > 0.10 or min(f_out, 1 - f_out) > 0.42:
                                note.append("PALINDROMIC: orientation NOT verified (outcome freq %.3f vs %s freq %.3f)" % (f_out, "UKB-PPP" if pop == "EUR" else "CKB", fref))
                            else:
                                note.append("palindromic: orientation consistent with frequency (outcome %.3f vs %.3f)" % (f_out, fref))
                        elif f_out == f_out and fref == fref and abs(f_out - fref) > 0.20:
                            note.append("FLAG: raising-allele frequency outcome %.3f vs reference %.3f" % (f_out, fref))
                        note.append("join %s on GRCh%s %d; outcome row %s %s/%s; raising-allele freq in outcome %s" % (
                            how, build, pos, c["id"], c["ea"], c["oa"], "%.4f" % f_out if f_out == f_out else "NA"))
                        add("S2.x22.%s.%s.beta" % (locus, cell), "SNP-outcome beta per European protein-raising allele (%s %s)" % (
                            I.get("rsid") or I.get("id"), I["raise"]), label, b, se=c["se"], p=c["p"], notes="; ".join(note),
                            unit="log-odds" if "eGFR" not in label else ("ln eGFR" if "CKDGen" in label else "SD of RINT eGFR"))
                        add("S2.x22.%s.%s.se" % (locus, cell), "SE", label, c["se"])
                        add("S2.x22.%s.%s.p" % (locus, cell), "p (as printed by the resource)", label, c["p"])
                        joinlog.append({"locus": locus, "resource": rid, "cell": cell, "result": how, "pos": pos})
                    add("S2.x22.%s.%s.resource" % (locus, cell), "outcome resource of the cell", label, "%s [%s]" % (label, rid))
                    # census cross-check: only the CKB-instrument x EAS-outcome cell
                    cv, cn = "NA", "no census value for this cell (the census covers CKB instrument x EAS outcome only)"
                    if iname == "CKBinstr" and pop == "EAS" and L["census"]:
                        cr = census.get((locus, L["census"]))
                        if cr is None or not cr["instr"]:
                            cn = "census row absent"
                        else:
                            cpos = int(cr["instr_pos_hg38"])
                            if cpos == ckb_instr["pos38"] and {cr["instr_ea"], cr["instr_oa"]} == {ckb_instr["raise"], ckb_instr["other"]}:
                                cv = float(cr["by"]) * (1 if cr["instr_ea"] == ckb_instr["raise"] else -1)
                                cn = "census by %s per %s (resource %s) re-oriented to %s; census instrument = this cell's instrument" % (
                                    cr["by"], cr["instr_ea"], cr["resource_label"], ckb_instr["raise"])
                                if L["census"] == "eGFR":
                                    cn += ("; the census eGFR value is risk-aligned eGFR* = -eGFR (eas__census_run.py l.77 flip=(outcome == "
                                           "'eGFR'); eas__census_common.py l.487 by = -beta), so on the raw eGFR scale of this cell it is %+.4f" % (-cv))
                            else:
                                cn = "census instrument %s (%s, GRCh38 %s) differs from this cell's instrument; census by %s per %s" % (
                                    cr["instr"], cr["instr_label"], cr["instr_pos_hg38"], cr["by"], cr["instr_ea"])
                    elif iname == "CKBinstr" and pop == "EAS" and not L["census"]:
                        cn = "no census row for FGFR4"
                    add("S2.x22.%s.%s.census_value" % (locus, cell), "census-derived cross-check (EAS_CENSUS_RESULTS by)", "EAS_CENSUS_RESULTS", cv, notes=cn)
    W.write_tsv(W.OUT + "/results_twobytwo.tsv", list(out[0].keys()), out)
    W.write_tsv(W.OUT + "/twobytwo_joinlog.tsv", ["locus", "resource", "cell", "result", "pos", "rows", "null_rates"], joinlog)
    # compact human-readable table
    comp = []
    for locus in LOCI:
        for o in out:
            q = o["quantity_id"]
            if q.startswith("S2.x22.%s." % locus) and q.endswith(".beta"):
                cell = q.split(".")[3]
                cvs = [x for x in out if x["quantity_id"] == "S2.x22.%s.%s.census_value" % (locus, cell)]
                comp.append({"locus": locus, "cell": cell, "resource": o["resource"], "beta": o["value"], "se": o["se"], "p": o["p"],
                             "census_value": cvs[0]["value"] if cvs else "NA", "notes": o["notes"]})
    W.write_tsv(W.OUT + "/twobytwo_table.tsv", ["locus", "cell", "resource", "beta", "se", "p", "census_value", "notes"], comp)
    for c in comp:
        b = c["beta"]
        W.log("  %-7s %-26s %-38s beta %s  p %s  census %s" % (c["locus"], c["cell"], c["resource"][:38],
                                                             "%+.4f" % b if isinstance(b, float) else b,
                                                             "%.2g" % c["p"] if isinstance(c["p"], float) else c["p"],
                                                             "%+.4f" % c["census_value"] if isinstance(c["census_value"], float) else "NA"))


if __name__ == "__main__":
    main()
