#!/usr/bin/env python
"""l3fm__common.py -- shared units for the L3 bridge FINE-MAPPING step (user decision 2026-09-26: fine-map the
ambiguous bridge loci before deciding how the paper presents the 10 European bridges). Import-only; no side effects.

ONE implementation of every loader: the protein/metabolite windows and the coloc.abf engine are the D1 builder's
(l3__build_bridge_screen_eur.py, imported, never copied); the only new readers are (a) the extra fields the D1 loaders
drop (the GRCh37 position inside the UKB-PPP ID, the Chen rsid/eaf), cross-checked value-for-value against the D1
loaders, and (b) the external LD reference panel.

LD reference (declared; the project's own panel <E>/20260331-CKM/Paper 3/data/ld_reference/1kg_eur/EUR is NOT mounted
on this machine -- E: absent 2026-09-26): 'paper 6/analysis/g1000_eur/g1000_eur' = 1000 Genomes phase 3, the 503 EUR
samples (the same 503 as the engine's panel, whose log reads '503 people'), bi-allelic SNPs, GRCh37, MAGMA build
(README in that folder). The engine's panel carried 8,550,156 variants against this one's 22,665,064, i.e. it was
frequency-filtered; we apply reference MAF >= 0.01 (declared, not re-derived from the engine, whose script is not local).

SECOND LD PANEL (added 2026-09-26, promoting the independent verifier's decisive method): 'lhcmr_eas_refs/1000G_EUR_Phase3_
plink/1000G.EUR.QC.<chr>' = the LDSC 1000G phase 3 EUR per-chromosome fileset, 489 samples. CAUTION, measured on the .fam
files: all 489 are among panel 1's 503 (489/489 IID overlap), so panel 2 is NOT an independent reference sample -- it is
panel 1 minus 14 people with a different variant QC. Agreement across the two panels therefore tests robustness to a
small perturbation of the LD estimate (sample and SNP set), not robustness to a different population. l3fm__panel2.py
builds its extracts; P2_AF_DIFF_MAX guards allele coding (the two panels share samples, so A1 frequencies must agree).
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
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
import re
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))            # .../CKM papers


def _load(name):
    s = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


B = _load("l3__build_bridge_screen_eur")          # D1 builder: windows, manifests, coloc engine (import has no side effects)

RAW_FM = os.path.join(ROOT, "V3", "V4_local", "raw", "l3-finemap")
STAGE_FM = os.path.join(ROOT, "V3", "V4_local", "staging", "l3-finemap")
D1_STAGE = os.path.join(ROOT, "V3", "V4_local", "staging", "d1-bridge-eur")
REF_PREFIX = os.path.join(ROOT, "paper 6", "analysis", "g1000_eur", "g1000_eur")
REF_N = 503
REF2_PREFIX_FMT = os.path.join(ROOT, "lhcmr_eas_refs", "1000G_EUR_Phase3_plink", "1000G.EUR.QC.%s")   # per chromosome
REF2_N = 489
P2_AF_DIFF_MAX = 0.05       # |A1 frequency panel 2 - panel 1| above this -> SNP dropped from the panel-2 set (coding guard)
NEIGHBOUR_WIN_DIR = os.path.join(ROOT, "V3", "V4_local", "raw", "l3-finemap", "ukbppp_eur_windows")   # l3fm__extract_neighbours.py
OLINK_MAP = os.path.join(ROOT, "weight-loss.omics", "data", "proteome", "olink_protein_map_3k_v1.tsv")
SYN_MAP = os.path.join(ROOT, "weight-loss.omics", "data", "proteome", "synapse_discovery_map.tsv")
WUTTKE_EA = (_P31R_EXT + "/D/data/ckdgen/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz")
ENGINE_LOG_29 = os.path.join(ROOT, "V3", "paper3_deposit_request_v4", "scripts", "reports", "Paper3_Program", "run_logs",
                             "RERUN_20260915", "29_cand_egfr.log")
DASHED = os.path.join(ROOT, "V3", "V4", "supp_tables", "SuppTable_coloc_robustness_dashed.csv")

REF_MAF_MIN = 0.01          # reference-panel MAF filter (declared)
AF_DIFF_MAX = 0.20          # |aligned GWAS AF - reference AF| above this -> dropped (allele/strand error guard)
PALIN_MAF_MAX = 0.40        # palindromic SNPs kept only when unambiguous by frequency
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# the four ambiguous loci + the ARG1 control. Metabolite lists are NOT typed: every D1 worklist key of the bridge protein.
# 'named_neighbours' are the genes the user named; the full neighbour inventory is DERIVED from the Olink map.
LOCI = {
    "BRAP_12q24": {"bridge": "BRAP", "named_neighbours": ["SH2B3", "ATXN2"], "role": "ambiguous"},
    "SMPD3_16q22": {"bridge": "SMPD3", "named_neighbours": ["SLC7A6"], "role": "ambiguous"},
    "GSTA3_6p12": {"bridge": "GSTA3", "named_neighbours": ["GSTA1", "GSTA2"], "role": "ambiguous"},
    "DPEP1_16q24": {"bridge": "DPEP1", "named_neighbours": ["CHMP1A"], "role": "ambiguous"},
    "ARG1_6q23": {"bridge": "ARG1", "named_neighbours": [], "role": "control_clean_bridge"},
}


class Refuse(SystemExit):
    pass


def refuse(msg):
    raise Refuse("REFUSED: " + msg)


def md5_file(path):
    return B.md5_file(path)


def rel(p):
    return B.rel(p)


# ---------------------------------------------------------------- inventories
def olink_map():
    """gene -> list of Olink assays (OID, chr, start, end; GRCh38 coordinates as in the UKB-PPP map)."""
    out = defaultdict(list)
    for r in csv.DictReader(io.open(OLINK_MAP, encoding="utf-8"), delimiter="\t"):
        out[r["HGNC.symbol"]].append({"oid": r["OlinkID"], "chr": r["chr"], "start": int(r["gene_start"]),
                                      "end": int(r["gene_end"]), "assay": r["Assay"], "panel": r["Panel"]})
    return out


def syn_map():
    return {r["OID"]: (r["synid"], r["tarname"]) for r in csv.DictReader(io.open(SYN_MAP, encoding="utf-8"), delimiter="\t")}


def local_protein_window(oid, ch):
    """European-discovery window on disk for an OID, or None. Three locations only (all syn51365303-derived): the D1
    extracts, the weight-loss.omics windows, and the neighbour extracts of l3fm__extract_neighbours.py."""
    for p in (os.path.join(B.EUR_WIN_DIR, "%s_chr%s_cis.tsv" % (oid, ch)), os.path.join(B.WL_PQTL, "%s_chr%s_cis.tsv" % (oid, ch)),
              os.path.join(NEIGHBOUR_WIN_DIR, "%s_chr%s_cis.tsv" % (oid, ch))):
        if os.path.exists(p):
            return p
    return None


def locus_units(units=None):
    """locus -> dict(bridge unit list, window) derived from the D1 worklist units."""
    units = units if units is not None else B.derive_units()
    out = {}
    for loc, spec in LOCI.items():
        us = [u for u in units if u["protein"] == spec["bridge"]]
        if not us:
            refuse("locus %s: bridge protein %s has no worklist keys" % (loc, spec["bridge"]))
        wins = {(u["chr"], u["lo"], u["hi"]) for u in us}
        if len(wins) != 1:
            refuse("locus %s: bridge keys carry %d different windows" % (loc, len(wins)))
        ch, lo, hi = wins.pop()
        out[loc] = {"units": us, "chr": ch, "lo": lo, "hi": hi, "oid": us[0]["oid"], **spec}
    return out


def neighbour_inventory(loci):
    """every Olink protein whose gene overlaps the locus window (GRCh38), plus the user-named genes whether or not assayed."""
    om, sm = olink_map(), syn_map()
    rows = []
    for loc, L in loci.items():
        seen = set()
        for g, assays in sorted(om.items()):
            for a in assays:
                if a["chr"] == L["chr"] and a["start"] <= L["hi"] and a["end"] >= L["lo"] and g != L["bridge"]:
                    p = local_protein_window(a["oid"], L["chr"])
                    syn = sm.get(a["oid"], ("", ""))
                    rows.append({"locus": loc, "gene": g, "oid": a["oid"], "gene_start": a["start"], "gene_end": a["end"],
                                 "named_by_user": g in L["named_neighbours"], "on_olink": True,
                                 "local_window": rel(p) if p else "", "synapse_id": syn[0], "synapse_tar": syn[1],
                                 "status": "LOCAL" if p else "NOT_LOCAL (request: Synapse %s, parent syn51365303)" % (syn[0] or "?")})
                    seen.add(g)
        for g in L["named_neighbours"]:
            if g in seen:
                continue
            if g in om:
                a = om[g][0]
                p = local_protein_window(a["oid"], a["chr"])
                syn = sm.get(a["oid"], ("", ""))
                rows.append({"locus": loc, "gene": g, "oid": a["oid"], "gene_start": a["start"], "gene_end": a["end"],
                             "named_by_user": True, "on_olink": True, "local_window": rel(p) if p else "",
                             "synapse_id": syn[0], "synapse_tar": syn[1],
                             "status": ("LOCAL" if p else "NOT_LOCAL (request: Synapse %s)" % syn[0]) + " [gene outside window]"})
            else:
                rows.append({"locus": loc, "gene": g, "oid": "", "gene_start": "", "gene_end": "", "named_by_user": True,
                             "on_olink": False, "local_window": "", "synapse_id": "", "synapse_tar": "",
                             "status": "NOT_ASSAYED (no Olink Explore 3072 assay -> no UKB-PPP protein exists to request)"})
    return rows


# ---------------------------------------------------------------- readers that keep what the D1 loaders drop
def read_protein_full(path, ch, lo, hi):
    """UKB-PPP 9-column window: rows in [lo, hi] (GRCh38) keyed (chr, pos38, alleles); keeps pos37 from the ID.
    Cross-checked against B.load_eur_protein for every key the D1 loader keeps."""
    out, dup = {}, set()
    with io.open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["chrom"] != ch:
                refuse("%s carries chromosome %s" % (path, r["chrom"]))
            p38 = int(r["pos"])
            if p38 < lo or p38 > hi:
                continue
            idf = r["ID"].split(":")
            if len(idf) < 4 or idf[0] != ch:
                refuse("%s: malformed ID %r" % (path, r["ID"]))
            if (idf[2].upper(), idf[3].upper()) != (r["ALLELE0"].upper(), r["ALLELE1"].upper()):
                refuse("%s: ID alleles differ from ALLELE0/ALLELE1 at %s" % (path, r["ID"]))
            try:
                v = {"pos38": p38, "pos37": int(idf[1]), "ea": r["ALLELE1"].upper(), "oa": r["ALLELE0"].upper(),
                     "b": float(r["BETA"]), "se": float(r["SE"]), "eaf": float(r["A1FREQ"])}
            except ValueError:
                continue
            k = (ch, p38, frozenset((v["ea"], v["oa"])))
            if k in out:
                dup.add(k)
            out[k] = v
    for k in dup:
        out.pop(k, None)          # an exact duplicate key is ambiguous -> dropped (counted by the caller)
    d1, _, _ = B.load_eur_protein(path, ch, lo, hi)
    bad = 0
    for k, v in out.items():
        w = d1.get((k[0], k[1]))
        if w is not None and (w["ea"], w["oa"]) == (v["ea"], v["oa"]) and (w["b"] != v["b"] or w["se"] != v["se"]):
            bad += 1
    if bad:
        refuse("%s: %d keys disagree with the D1 loader" % (path, bad))
    return out, len(dup)


def read_metab_full(path):
    """Chen harmonised window keyed (chr, pos38, alleles), keeping rsid + eaf; cross-checked against B.load_metab_window."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    hdr = lines[0].rstrip("\r").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    out, dup = {}, set()
    for line in lines[1:]:
        if not line:
            continue
        f = line.split("\t")
        try:
            c = f[ci["chromosome"]].replace("chr", ""); p = int(f[ci["base_pair_location"]])
            v = {"pos38": p, "ea": f[ci["effect_allele"]].upper(), "oa": f[ci["other_allele"]].upper(),
                 "b": float(f[ci["beta"]]), "se": float(f[ci["standard_error"]]),
                 "eaf": float(f[ci["effect_allele_frequency"]]) if f[ci["effect_allele_frequency"]] not in ("", "NA") else float("nan"),
                 "rsid": f[ci["rsid"]] if "rsid" in ci else ""}
        except (KeyError, ValueError, IndexError):
            continue
        k = (c, p, frozenset((v["ea"], v["oa"])))
        if k in out:
            dup.add(k)
        out[k] = v
    for k in dup:
        out.pop(k, None)
    d1, _, _, _ = B.load_metab_window(path)
    bad = sum(1 for k, v in out.items() if (k[0], k[1]) in d1 and (d1[(k[0], k[1])]["ea"], d1[(k[0], k[1])]["oa"]) == (v["ea"], v["oa"])
              and (d1[(k[0], k[1])]["b"] != v["b"] or d1[(k[0], k[1])]["se"] != v["se"]))
    if bad:
        refuse("%s: %d keys disagree with the D1 loader" % (path, bad))
    return out, len(dup)


def read_wuttke_region(path, ch, lo37, hi37):
    """Wuttke 2019 EA eGFR (GRCh37, whitespace-separated; split on whitespace runs -- the padded-Freq1 lesson)."""
    out, dup, n = {}, set(), 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        hdr = fh.readline().split()
        want = ["Chr", "Pos_b37", "RSID", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "n_total_sum"]
        if hdr != want:
            refuse("Wuttke header %r != %r" % (hdr, want))
        for line in fh:
            if "\t" in line:
                refuse("Wuttke row contains a tab (reader assumes whitespace-separated)")
            f = line.split()
            if f[0] != ch:
                continue
            p = int(f[1])
            if p < lo37 or p > hi37:
                continue
            n += 1
            v = {"pos37": p, "rsid": f[2], "ea": f[3].upper(), "oa": f[4].upper(), "eaf": float(f[5]),
                 "b": float(f[6]), "se": float(f[7]), "n": float(f[9])}
            k = (ch, p, frozenset((v["ea"], v["oa"])))
            if k in out:
                dup.add(k)
            out[k] = v
    for k in dup:
        out.pop(k, None)
    return out, len(dup), n


# ---------------------------------------------------------------- reference panel
def ref_bim_region(ch, lo37, hi37):
    """rows of the reference .bim in [lo37, hi37] with their variant index (bed row)."""
    rows = []
    idx = -1
    with open(REF_PREFIX + ".bim", "r", encoding="ascii") as fh:
        for line in fh:
            idx += 1
            if not line.startswith(ch + "\t"):
                continue
            f = line.split("\t")
            p = int(f[3])
            if lo37 <= p <= hi37:
                rows.append({"idx": idx, "rsid": f[1], "pos37": p, "a1": f[4].upper(), "a2": f[5].strip().upper()})
    return rows


def ref_bim_regions(regions):
    """one pass over the bim for several (chr, lo37, hi37) regions -> {region: rows}."""
    out = {r: [] for r in regions}
    bych = defaultdict(list)
    for r in regions:
        bych[r[0]].append(r)
    idx = -1
    with open(REF_PREFIX + ".bim", "r", encoding="ascii") as fh:
        for line in fh:
            idx += 1
            t = line.find("\t")
            ch = line[:t]
            if ch not in bych:
                continue
            f = line.split("\t")
            p = int(f[3])
            for r in bych[ch]:
                if r[1] <= p <= r[2]:
                    out[r].append({"idx": idx, "rsid": f[1], "pos37": p, "a1": f[4].upper(), "a2": f[5].strip().upper()})
    return out


def read_bed_dosage(indices, n=REF_N, prefix=REF_PREFIX):
    """A1-allele dosage (0/1/2, NaN missing) for the given bed variant indices, shape (n, len(indices)).
    prefix/n default to panel 1; panel 2 passes REF2_PREFIX_FMT % chr and REF2_N."""
    nb = (n + 3) // 4
    with open(prefix + ".bed", "rb") as fh:
        magic = fh.read(3)
        if magic != b"\x6c\x1b\x01":
            refuse("reference .bed is not a variant-major PLINK1 file")
        G = np.empty((n, len(indices)), dtype=np.float64)
        lut = np.array([2.0, np.nan, 1.0, 0.0])        # 00 hom A1, 01 missing, 10 het, 11 hom A2  -> A1 dosage
        for j, i in enumerate(indices):
            fh.seek(3 + i * nb)
            by = np.frombuffer(fh.read(nb), dtype=np.uint8)
            codes = np.stack([(by >> s) & 3 for s in (0, 2, 4, 6)], axis=1).reshape(-1)[:n]
            G[:, j] = lut[codes]
    return G


def fam_n(prefix=REF_PREFIX):
    with open(prefix + ".fam", "r", encoding="ascii") as fh:
        return sum(1 for _ in fh)


def fam_iids(prefix):
    with open(prefix + ".fam", "r", encoding="ascii") as fh:
        return [ln.split()[1] for ln in fh if ln.strip()]


def ref2_bim_region(ch, lo37, hi37):
    """rows of the panel-2 per-chromosome .bim in [lo37, hi37] with their variant index (bed row of THAT chromosome file)."""
    rows = []
    with open((REF2_PREFIX_FMT % ch) + ".bim", "r", encoding="ascii") as fh:
        for idx, line in enumerate(fh):
            f = line.split()
            if f[0] != ch:
                refuse("panel-2 bim for chr%s carries chromosome %s" % (ch, f[0]))
            p = int(f[3])
            if lo37 <= p <= hi37:
                rows.append({"idx": idx, "rsid": f[1], "pos37": p, "a1": f[4].upper(), "a2": f[5].upper()})
    return rows


# ---------------------------------------------------------------- harmonisation to the reference counted allele (A1)
def align_to_ref(ea, oa, a1, a2):
    """+1 when ea==A1 (effect already on the counted allele), -1 when ea==A2, 0 when not matchable.
    Exact allele matching only for palindromic SNPs (strand-ambiguous); complement allowed otherwise."""
    pal = {ea, oa} in ({"A", "T"}, {"C", "G"})
    if (ea, oa) == (a1, a2):
        return 1
    if (ea, oa) == (a2, a1):
        return -1
    if not pal:
        if (COMP.get(ea), COMP.get(oa)) == (a1, a2):
            return 1
        if (COMP.get(ea), COMP.get(oa)) == (a2, a1):
            return -1
    return 0


def is_palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def n_hat(b_se_eaf):
    """effective N of a standardised trait from its per-SNP SEs: median over 0.05 < MAF of 1 / (2 f (1-f) se^2)."""
    v = [1.0 / (2 * f * (1 - f) * se * se) for _, se, f in b_se_eaf if f == f and 0.05 < min(f, 1 - f) and se > 0]
    if len(v) < 50:
        refuse("n_hat: only %d usable SNPs" % len(v))
    return float(np.median(v))


def write_tsv(path, cols, rows, gz=False):
    tmp = path + ".part"
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in cols})
    data = buf.getvalue().encode("utf-8")
    with open(tmp, "wb") as raw:
        if gz:       # gzip header mtime 0: identical content -> identical md5, so downstream stamps stay fresh on a rebuild
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as g:
                g.write(data)
        else:
            raw.write(data)
    os.replace(tmp, path)


def write_json(path, obj):
    with open(path + ".part", "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(path + ".part", path)


def engine_umod_window():
    """UMOD~eGFR SuSiE window as the V4 engine logged it (read from the log, not typed)."""
    with io.open(ENGINE_LOG_29, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.search(r"susie: UMOD ~ eGFR \| chr(\d+):(\d+)-(\d+)", line)
            if m:
                return m.group(1), int(m.group(2)), int(m.group(3))
    refuse("UMOD~eGFR window line not found in %s" % ENGINE_LOG_29)


def engine_umod_expected():
    """the shipped strict SuSiE UMOD~eGFR row (V4 supp table) -- the pipeline positive control's target."""
    for r in csv.DictReader(io.open(DASHED, encoding="utf-8")):
        if r["gene"] == "UMOD" and r["outcome"] == "eGFR":
            return {k: r[k] for k in ("pp_h4_abf", "pp_h4_susie_orig", "pp_h4_susie_strict", "pp_h3_susie_strict", "run_status", "n_warn")}
    refuse("UMOD/eGFR row absent from %s" % DASHED)
