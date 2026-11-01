#!/usr/bin/env python
"""eas__census_common.py -- shared code for TRACK D2, the V4 native East Asian census.

ONE implementation of every piece the census scripts share (paths, write guard, md5, the pinned
engine, liftover, the outcome-resource registry, readers, the pair runner). The fetcher, slicer,
runner, controls and S23 builder all import this module; nothing here is copied into them.

Run every census script with  python -B  (or PYTHONDONTWRITEBYTECODE=1) so no __pycache__ is written.

WRITE RULE: every write goes through atomic_write()/open_out(), which refuse any path outside
V3/V4_local/raw/d2-eas-census, V3/V4_local/staging/d2-eas-census, <D>/data/paper31_v4_raw/d2-eas-census
or the track scratch folder. A partial file never replaces a good one (<name>.part then os.replace).
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
import io
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True

ROOT = _P31R_CKM
V4L = ROOT + "/V3/V4_local"
BUILDERS = V4L + "/builders"
VENDOR = BUILDERS + "/_vendor"
RAW = V4L + "/raw/d2-eas-census"
STAGE = V4L + "/staging/d2-eas-census"
DRAW = (_P31R_EXT + "/D/data/paper31_v4_raw/d2-eas-census")
SCRATCH = ((_P31R_EXT + "/scratch"
))
PREP = RAW + "/prep"   # the D2 prep folder: the author's 2026-09-27 copy of the session scratchpad's prep/eas-census
                       # (raw/d2-eas-census/prep/PREP_MD5.txt); the two tables the census step reads are pinned here
PREP_MD5 = {"fetch_regions_head.csv": "4e80426bbe74e77ad53a7f908bc26822",
            "fetch_regions_head_local7.csv": "6f88edc6776a983e9ed46532a396dbbf"}
WRITABLE = [RAW, STAGE, DRAW, SCRATCH]

# ---------------------------------------------------------------- inputs (read-only, md5-pinned)
TIER = ROOT + "/V3/v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv"
TIER_MD5 = "31576d124c8eca551fc8e357acbc5520"
INSTR = ROOT + "/V3/V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv"
INSTR_MD5 = "7170afca7086d63918093b445a8fe6ff"
EUR_MR = ROOT + "/V3/V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv"
EUR_COLOC = ROOT + "/V3/V4/source_data/tables_from_rds/phase3_coloc__coloc_results_final.csv"
PHEN = ROOT + "/paper 3.eas/analysis/data/ckb_meta/ckb_phenotypes.json"
PHEN_MD5 = "254bfd16ae7e0ab6617f5582d79f5e8a"
TOPHITS_PREP = PREP + "/ckb_top_hits.tsv"
TOPHITS_MD5 = "bfb8252fa625e4a68039cedbc6820130"
TOPHITS_URL = "https://pheweb.ckbiobank.org/download/top_hits.tsv"
TOPHITS_ETAG = "1772686289.3223631-1297534-2320963119"
TOPHITS_BYTES = 1297534
EAS = ROOT + "/paper 3.eas/analysis"
PILOT_EXP = EAS + "/data/eas_exposure"
PILOT_OUT = EAS + "/data/eas_outcomes"
PILOT_RES = EAS + "/results"
PFKFB2_LOCAL = EAS + "/data/verify_pfkfb2/PFKFB2_ckb_region.tsv"
S23_CURRENT = V4L + "/source_data/figures/S23_four_quadrant_eur_eas/S23_four_quadrant__source.csv"
S23_CURRENT_MD5 = "1d60355284f03058dc9c459cc3c20f04"
S23_BUILDER = BUILDERS + "/eas__build_four_quadrant.py"
S23_BUILDER_MD5 = "73e5f5f805836e48e777910d1529c541"
T2D_UPGRADE = ROOT + "/paper 3 integrated/analysis/t2d_upgrade_results.csv"
T2D_UPGRADE_MD5 = "96f1ed66c2cab8ea3df4bdb5b9a7fa5f"
CHAIN = (_P31R_EXT + "/.pyliftover/hg38ToHg19.over.chain.gz")
CHAIN_MD5 = "ff3031d93792f4cbb86af44055efd903"
G1K_BIM = ROOT + "/paper 4/independent_build/data/g1000_eas_merged.bim"
G1K_BED = ROOT + "/paper 4/independent_build/data/g1000_eas_merged.bed"
G1K_FAM = ROOT + "/paper 4/independent_build/data/g1000_eas_merged.fam"
ENGINE_SRC = EAS + "/engine_ref/mr_coloc.py"
ENGINE_MD5 = "fa6db396645e411504c4825864f75dd4"
ENGINE_VENDOR = VENDOR + "/eas__census_mr_coloc.py"
TABIX_VENDOR = VENDOR + "/eas__census_tabix_remote.py"
TABIX_MD5 = "391fcc8c5385c72a8a095f50b2f69210"
GS544_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/"
             "GCST90104544/harmonised/GCST90104544.h.tsv.gz")
GS544_TBI_LOCAL = EAS + "/data/GCST90104544.h.tsv.gz.tbi"
GS544_TBI_BYTES = 1490689

GWS = 5e-8
WIN = 500_000          # coloc window: lead +/- 500 kb ; instrument domain also restricted to lead +/- 500 kb
MAF_MIN = 0.01
MIN_COLOC_SNPS = 20
P1 = P2 = 1e-4
P12 = 1e-5
P12_STRICT = 1e-6


# ---------------------------------------------------------------- write guard + io
def _norm(p):
    return os.path.normcase(os.path.abspath(p)).replace("\\", "/")


def assert_writable(path):
    n = _norm(path)
    for w in WRITABLE:
        if n.startswith(_norm(w) + "/") or n == _norm(w):
            return
    raise SystemExit(f"WRITE REFUSED (outside the D2 track areas): {path}")


def ensure_dir(d):
    assert_writable(d)
    os.makedirs(d, exist_ok=True)


def atomic_write_bytes(path, data):
    assert_writable(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".part"
    with open(tmp, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def atomic_write_text(path, text):
    atomic_write_bytes(path, text.encode("utf-8"))


def write_table(path, rows, fields=None, sep="\t"):
    """Write rows (list of dicts) atomically. Returns (n_rows, md5)."""
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, delimiter=sep, lineterminator="\n", extrasaction="raise")
    w.writeheader()
    for r in rows:
        w.writerow({k: fmt(r.get(k)) for k in fields})
    atomic_write_text(path, buf.getvalue())
    return len(rows), md5(path)


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):          # includes numpy float64 (a float subclass whose repr is 'np.float64(x)')
        v = float(v)
        if math.isnan(v):
            return "NA"
        return repr(v)
    if type(v).__module__ == "numpy":
        return fmt(v.item())
    return str(v)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_md5(path, want, what=""):
    got = md5(path)
    if got != want:
        raise SystemExit(f"MD5 MISMATCH {what}: {path} is {got}, expected {want}")
    return got


def read_csv(path, sep=","):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=sep))


class Log:
    """Append-only, flushed progress log (one line per unit)."""

    def __init__(self, path):
        assert_writable(path)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.path = path

    def __call__(self, msg):
        line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}"
        print(line, flush=True)
        with open(self.path, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


# ---------------------------------------------------------------- pinned engine
def engine():
    """Import the vendored, md5-pinned pilot engine (byte copy of paper 3.eas engine_ref/mr_coloc.py)."""
    assert_md5(ENGINE_VENDOR, ENGINE_MD5, "vendored engine")
    if os.path.exists(ENGINE_SRC):
        assert_md5(ENGINE_SRC, ENGINE_MD5, "engine source")
    import importlib.util
    spec = importlib.util.spec_from_file_location("eas__census_mr_coloc", ENGINE_VENDOR)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tabix():
    assert_md5(TABIX_VENDOR, TABIX_MD5, "vendored tabix_remote")
    import importlib.util
    spec = importlib.util.spec_from_file_location("eas__census_tabix_remote", TABIX_VENDOR)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- alleles + liftover
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def aset(a, b):
    return frozenset((a.upper(), b.upper()))


def is_palindromic(a, b):
    a, b = a.upper(), b.upper()
    return len(a) == 1 and len(b) == 1 and COMP.get(a) == b


class Lifter:
    """hg38 -> hg19 with strand and multiplicity handled explicitly.
    lift(chrom, pos38, ref, alt) -> (pos19, ref19, alt19, status) ; status in
    ok | minus_complemented | unlifted | multimap | other_chrom | minus_palindrome | minus_indel."""

    def __init__(self):
        from pyliftover import LiftOver
        assert_md5(CHAIN, CHAIN_MD5, "hg38ToHg19 chain")
        self.lo = LiftOver(CHAIN)
        self.cache = {}

    def lift_pos(self, chrom, pos38):
        k = (chrom, pos38)
        if k not in self.cache:
            self.cache[k] = self.lo.convert_coordinate("chr" + str(chrom), int(pos38) - 1)
        return self.cache[k]

    def lift(self, chrom, pos38, ref, alt):
        c = self.lift_pos(chrom, pos38)
        if not c:
            return None, ref, alt, "unlifted"
        if len(c) > 1:
            return None, ref, alt, "multimap"
        ch, p0, strand = c[0][0], c[0][1], c[0][2]
        if ch != "chr" + str(chrom):
            return None, ref, alt, "other_chrom"
        if strand == "+":
            return p0 + 1, ref, alt, "ok"
        if len(ref) != 1 or len(alt) != 1:
            return None, ref, alt, "minus_indel"
        if is_palindromic(ref, alt):
            return None, ref, alt, "minus_palindrome"
        return p0 + 1, COMP[ref.upper()], COMP[alt.upper()], "minus_complemented"


# ---------------------------------------------------------------- outcome-resource registry
# effect-allele column recorded per resource (spec C): AGEN EA; BBJ SAIGE Allele2; TWB2 Allele1;
# CKD effect_allele; BBJ_AF_2022 ALT; Enzan effect_allele; GIGASTROKE-544 harmonised effect_allele.
OUTCOMES = {
    "T2D": dict(label="AGEN Spracklen 2020 T2D (includes CKB, KBA and BBJ)", accession="GCST010118",
                path=EAS + "/data/eas_outcomes/AGEN_T2D_full.txt.gz", md5="10abea38901654f8d96bc89bd10dbd44",
                build="hg19", sep="\t", chr="Chr", pos="Pos", ea="EA", oa="NEA", beta="Beta", se="SE", p="P",
                eaf="EAF", n="Neff", rsid=None, info=None, direction=None, binary=True,
                ncase=77418, nctrl=356122, ckb_overlap=True, expect_rows=11825584),
    "CAD": dict(label="BBJ CAD (hum0197.v5, pheweb.jp CAD)", accession="BBJ_CAD_hum0197.v5",
                path=EAS + "/data/eas_outcomes/BBJ_CAD.txt.gz", md5="93fec56d62ab2e371e0fade6a71b4fff",
                build="hg19", sep="\t", chr="CHR", pos="POS", ea="Allele2", oa="Allele1", beta="BETA", se="SE",
                p="p.value", eaf="AF_Allele2", n="N", rsid="SNPID", info="imputationInfo", direction=None,
                binary=True, ncase=32512, nctrl=146214, ckb_overlap=False, expect_rows=None),
    "eGFR": dict(label="TWB2+BBJ eGFR meta-analysis (N 244,952)", accession="TWB2_BBJ_eGFR_PMID39472450",
                 path=EAS + "/data/eas_outcomes/TWB2_BBJ_eGFR_hg19.gz", md5="30981f43003261a3c08ada6b2a9c4dc6",
                 build="hg19", sep=None, chr="CHR", pos="BP", ea="Allele1", oa="Allele2", beta="Effect",
                 se="StdErr", p="P", eaf=None, n=None, rsid="SNP", info=None, direction="Direction",
                 binary=False, ncase=None, nctrl=None, ntot=244952, ckb_overlap=False, expect_rows=4512904),
    "CKD": dict(label="BBJ chronic renal failure (Sakaue 2021)", accession="GCST90018602",
                path=ROOT + "/paper 6/analysis/data/outcomes/eas/CKD_GCST90018602.gz",
                md5="b91c07c37abded1045ac5cdc9932584c", build="hg19", sep="\t", chr="chromosome",
                pos="base_pair_location", ea="effect_allele", oa="other_allele", beta="beta", se="standard_error",
                p="p_value", eaf="effect_allele_frequency", n=None, rsid=None, info=None, direction=None,
                binary=True, ncase=2117, nctrl=174345, ckb_overlap=False, expect_rows=None),
    "AF": dict(label="BBJ atrial fibrillation 2022 (Miyazawa 2023)", accession="BBJ_AF_2022_hum0014.v29",
               path=EAS + "/data/eas_outcomes/BBJ_AF_2022.sumstats.gz", md5="496ce37cef199e4b1a9b91b0eebd90c0",
               build="hg19", sep="\t", chr="CHR", pos="POS", ea="ALT", oa="REF", beta="BETA", se="SE", p="P",
               eaf="AAF", n="N", rsid=None, info="INFO", direction=None, binary=True, ncase=9826, nctrl=140446,
               ckb_overlap=False, expect_rows=None),
    "HF": dict(label="Enzan 2025 heart failure (GWAS Catalog harmonised)", accession="GCST90668009",
               path=EAS + "/data/eas_outcomes/GCST90668009.h.tsv.gz", md5="f0d97595f5ba16f4cb921d9013e56758",
               build="hg38", sep="\t", chr="chromosome", pos="base_pair_location", ea="effect_allele",
               oa="other_allele", beta="beta", se="standard_error", p="p_value", eaf="effect_allele_frequency",
               n="n", rsid="rsid", info=None, direction=None, binary=True, ncase=16251, nctrl=197577,
               ckb_overlap=False, expect_rows=7786701),
    "Stroke": dict(label="GIGASTROKE East Asian any stroke (Mishra 2022; includes CKB)", accession="GCST90104544",
                   path=None, remote=GS544_URL, md5=None, build="hg38", sep="\t", chr="chromosome",
                   pos="base_pair_location", ea="effect_allele", oa="other_allele", beta="beta",
                   se="standard_error", p="p_value", eaf="effect_allele_frequency", n=None, rsid="rsid",
                   info=None, direction=None, binary=True, ncase=27413, nctrl=237242, ckb_overlap=True,
                   expect_rows=None),
    # ---- CKB-free sensitivity resources (local) ----
    "T2D_BBJ": dict(label="BBJ T2D (Sakaue-Kanai 2020; CKB-free)", accession="BBJ_T2D_hum0197.v3",
                    path=(_P31R_EXT + "/D/data/bbj/hum0197.v3.BBJ.T2D.v1/hum0197.v3.BBJ.T2D.v1/GWASsummary_T2D_Japanese_SakaueKanai2020.auto.txt.gz"),
                    md5=None, build="hg19", sep="\t", chr="CHR", pos="POS", ea="Allele2", oa="Allele1", beta="BETA",
                    se="SE", p="p.value", eaf="AF_Allele2", n="N", rsid="SNPID", info="imputationInfo",
                    direction=None, binary=True, ncase=45383, nctrl=132032, ckb_overlap=False, expect_rows=None),
    "Stroke_BBJIS": dict(label="BBJ ischaemic stroke (Sakaue-Kanai 2020; CKB-free)", accession="BBJ_IS_hum0197.v3",
                         path=(_P31R_EXT + "/D/data/bbj/hum0197.v3.BBJ.IS.v1/hum0197.v3.BBJ.IS.v1/GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz"),
                         md5=None, build="hg19", sep="\t", chr="CHR", pos="POS", ea="Allele2", oa="Allele1",
                         beta="BETA", se="SE", p="p.value", eaf="AF_Allele2", n="N", rsid="SNPID",
                         info="imputationInfo", direction=None, binary=True, ncase=22664, nctrl=152022,
                         ckb_overlap=False, expect_rows=None),
}
CKBFREE_FOR = {"T2D": "T2D_BBJ", "Stroke": "Stroke_BBJIS"}   # every other default resource is already CKB-free
NOT_TESTABLE = {"NAFLD": "NOT_TESTABLE_NO_EAS_GWAS"}


# ---------------------------------------------------------------- CKB region files
REGION_COLS = ["af", "alt", "beta", "chr", "end", "id", "nearest_genes", "position", "pvalue", "ref", "rsid", "sebeta"]


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def read_region(path, require_rsid=False):
    """Read a CKB region TSV (pilot schema). Returns (rows, n_bad). Numerics parsed; rows keep file order."""
    rows, bad = [], 0
    with open(path, newline="", encoding="utf-8") as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            b, s, p = fnum(r.get("beta")), fnum(r.get("sebeta")), fnum(r.get("pvalue"))
            try:
                pos = int(r["position"])
            except (TypeError, ValueError, KeyError):
                bad += 1
                continue
            if b is None or s is None or p is None or s <= 0:
                bad += 1
                continue
            rs = (r.get("rsid") or "").strip()
            if rs in (".", "NA"):
                rs = ""
            if require_rsid and not rs:
                continue
            rows.append(dict(chr=str(r["chr"]).replace("chr", ""), position=pos, ref=r["ref"].upper(),
                             alt=r["alt"].upper(), beta=b, sebeta=s, pvalue=p, af=fnum(r.get("af")), rsid=rs,
                             id=r.get("id", ""), idx=i))
    return rows, bad


def maf(af):
    return None if af is None else min(af, 1.0 - af)


def zabs(r):
    return abs(r["beta"] / r["sebeta"])


def key_new(r):
    """User-decided tie rule: reported p ascending, |beta/sebeta| descending, position ascending."""
    return (r["pvalue"], -zabs(r), r["position"])


def key_old(r):
    """Pilot rule: min(p), first in file order."""
    return (r["pvalue"], r["idx"])


def key_maxz(r):
    """The LITERAL 'largest |z| first' reading (NOT the decided rule; used only by the C4 mutation)."""
    return (-zabs(r), r["pvalue"], r["position"])


def vid(r):
    return f"{r['chr']}:{r['position']}:{r['ref']}:{r['alt']}"


# ---------------------------------------------------------------- slices (normalised outcome rows)
SLICE_COLS = ["chr", "pos", "ea", "oa", "beta", "se", "p", "eaf", "n", "rsid", "info", "direction"]


def read_slice(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            b, s, p = fnum(r["beta"]), fnum(r["se"]), fnum(r["p"])
            if b is None or s is None or s <= 0:
                continue
            rows.append(dict(chr=r["chr"], pos=int(r["pos"]), ea=r["ea"].upper(), oa=r["oa"].upper(), beta=b,
                             se=s, p=p, eaf=fnum(r["eaf"]), n=r["n"], rsid=r["rsid"], info=r["info"],
                             direction=r["direction"]))
    return rows


def read_pilot_outcome(path):
    """Pilot outcome slice (pilot schema: alt = effect allele). Returns normalised rows + file-order index."""
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            b, s, p = fnum(r.get("beta")), fnum(r.get("sebeta")), fnum(r.get("pvalue"))
            if b is None or s is None or p is None:
                continue
            try:
                pos = int(float(r["position"]))
            except (TypeError, ValueError, KeyError):
                continue
            rs = (r.get("rsid") or "").strip()
            if rs in (".", "NA"):
                rs = ""
            rows.append(dict(chr=str(r["chr"]).replace("chr", ""), pos=pos, ea=r["alt"].upper(), oa=r["ref"].upper(),
                             beta=b, se=s, p=p, eaf=fnum(r.get("af")), n="", rsid=rs, info="", direction="", idx=i))
    return rows


# ---------------------------------------------------------------- the pair runner (one code path for all modes)
def run_pair(E, exp, out, cfg, lifter=None):
    """Estimate one protein -> outcome pair.

    exp  : CKB rows from read_region (hg38)
    out  : normalised outcome rows (read_slice / read_pilot_outcome)
    cfg  : dict with
        join          'rsid' | 'posallele'
        rule          'old' | 'new'
        domain        'file'  (pilot: every row of the exposure file) |
                      'census' (GWS in gene window AND lead +/- 500 kb, MAF >= 0.01)
        coloc_domain  'file' | 'lead500'
        require_rsid  drop exposure rows without an rsID (pilot load())
        flip          negate the outcome beta (eGFR -> eGFR*)
        lift          True when the outcome is hg19 (exposure lifted hg38 -> hg19)
        gene_window   (chr, lo, hi) hg38 -- census only
        fixed_instr   optional variant id 'chr:pos:ref:alt' (hold the instrument, C5)
        p12s          list of p12 values for coloc (first is the primary)
    Returns a flat result dict.
    """
    import numpy as np
    rule = {"new": key_new, "old": key_old, "maxz": key_maxz}[cfg["rule"]]
    if cfg.get("require_rsid"):
        exp = [r for r in exp if r["rsid"]]
    res = dict(join=cfg["join"], rule=cfg["rule"], domain=cfg["domain"], coloc_domain=cfg["coloc_domain"])

    # ---- outcome index
    if cfg["join"] == "rsid":
        oidx = {}
        for o in out:
            if o["rsid"]:
                oidx[o["rsid"]] = o            # pilot semantics: a later duplicate rsID overwrites
        res["out_dup_keys_dropped"] = 0
    else:
        oidx, dups = {}, set()
        for o in out:
            k = (o["chr"], o["pos"], aset(o["ea"], o["oa"]))
            if k in oidx:
                dups.add(k)
            oidx[k] = o
        for k in dups:
            del oidx[k]                        # ambiguous: two outcome rows for one variant
        res["out_dup_keys_dropped"] = len(dups)

    lift_counts = {}

    def okey(r):
        if cfg["join"] == "rsid":
            return r["rsid"] or None, r["ref"], r["alt"]
        if cfg.get("lift"):
            p19, ref, alt, st = lifter.lift(r["chr"], r["position"], r["ref"], r["alt"])
            lift_counts[st] = lift_counts.get(st, 0) + 1
            if p19 is None:
                return None, ref, alt
            return (r["chr"], p19, aset(ref, alt)), ref, alt
        return (r["chr"], r["position"], aset(r["ref"], r["alt"])), r["ref"], r["alt"]

    okc = {}

    def okey_full(r):
        if id(r) not in okc:
            okc[id(r)] = okey(r)
        return okc[id(r)]

    def join(r):
        k, ref, alt = okey_full(r)
        if k is None:
            return None
        o = oidx.get(k)
        if o is None:
            return None
        by = -o["beta"] if cfg["flip"] else o["beta"]
        h = E.harmonise(alt, ref, o["ea"], o["oa"], by, o["se"])
        if h is None:
            return None
        return o, h[0], h[1]

    # ---- absolute lead and windows
    if cfg["domain"] == "census":
        gchr, glo, ghi = cfg["gene_window"]
        exp = [r for r in exp if r["chr"] == gchr]
        in_gene = [r for r in exp if glo <= r["position"] <= ghi]
        sig = [r for r in in_gene if r["pvalue"] < GWS]
        if not sig:
            res.update(status="NO_GWS_EXPOSURE_IN_WINDOW")
            return res
        abs_lead = min(sig, key=rule)
        lo, hi = abs_lead["position"] - WIN, abs_lead["position"] + WIN
        pool = [r for r in sig if lo <= r["position"] <= hi and (maf(r["af"]) or 0) >= MAF_MIN]
        cwin = [r for r in exp if lo <= r["position"] <= hi]
        res.update(coloc_lo_hg38=lo, coloc_hi_hg38=hi, n_exp_gene_window=len(in_gene), n_exp_gws=len(sig),
                   n_exp_coloc_window=len(cwin))
    else:
        sig = [r for r in exp if r["pvalue"] < GWS]
        abs_lead = min(sig if sig else exp, key=rule)
        pool = sig
        cwin = exp
        res.update(n_exp_file=len(exp), n_exp_gws=len(sig))
    res.update(abs_lead=vid(abs_lead), abs_lead_rsid=abs_lead["rsid"], abs_lead_p=abs_lead["pvalue"],
               abs_lead_maf=maf(abs_lead["af"]), n_out_rows=len(out))

    # ---- join every candidate row once (pool and coloc window), cache by object id
    cache = {}

    def jc(r):
        if id(r) not in cache:
            cache[id(r)] = join(r)
        return cache[id(r)]

    def present(r):
        """Pilot semantics (harmonise_in_domain False): the rsID/key merely exists in the outcome.
        Census semantics (True): the variant is present AND harmonisable."""
        if cfg.get("harmonise_in_domain", True):
            return jc(r) is not None
        k = okey_full(r)[0]
        return k is not None and k in oidx

    # ---- instrument
    shared = [r for r in pool if present(r)]
    if cfg.get("fixed_instr"):
        shared = [r for r in shared if vid(r) == cfg["fixed_instr"]]
        if len(shared) != 1:
            res.update(status="FIXED_INSTRUMENT_NOT_SHARED")
    if shared:
        lead = min(shared, key=rule)
        topp = [r for r in shared if r["pvalue"] == lead["pvalue"]]
        topz = [r for r in topp if zabs(r) == zabs(lead)]
        res.update(instr=vid(lead), instr_rsid_ckb=lead["rsid"], instr_pos_hg38=lead["position"],
                   instr_ea=lead["alt"], instr_oa=lead["ref"], instr_p=lead["pvalue"], instr_maf=maf(lead["af"]),
                   bx=lead["beta"], bx_se=lead["sebeta"], F=(lead["beta"] / lead["sebeta"]) ** 2,
                   is_abs_lead=vid(lead) == vid(abs_lead), n_domain=len(pool), n_domain_shared=len(shared),
                   n_tied_top_p=len(topp), n_tied_top_z=len(topz))
        j = jc(lead)
        if j is None:
            res.update(mr_status="HARMONISE_FAIL (palindromic/incompatible)")
        else:
            o, by, byse = j
            b, s, p = E.wald_ratio(lead["beta"], by, byse)
            res.update(mr_status="ok", instr_out_pos=o["pos"], instr_out_rsid=o["rsid"], by=by, by_se=byse,
                       by_p=o["p"], mr_beta=b, mr_se=s, mr_p=p)
    else:
        lead = abs_lead
        res.update(mr_status="NO_GWS_SHARED_INSTRUMENT", instr=vid(abs_lead), instr_rsid_ckb=abs_lead["rsid"],
                   is_abs_lead=True, n_domain=len(pool), n_domain_shared=0, n_tied_top_p=0, n_tied_top_z=0,
                   F=(abs_lead["beta"] / abs_lead["sebeta"]) ** 2)
    res["instr_is_proxy"] = not res["is_abs_lead"]

    # ---- coloc
    bx, bxse, byy, byse, snps, extra = [], [], [], [], [], []
    for r in cwin:
        j = jc(r)
        if j is None:
            continue
        o, b_y, s_y = j
        bx.append(r["beta"]); bxse.append(r["sebeta"]); byy.append(b_y); byse.append(s_y); snps.append(vid(r))
        extra.append((maf(r["af"]), o.get("n", ""), o.get("direction", "")))
    res["coloc_nsnp"] = len(bx)
    res["_coloc_input"] = (snps, bx, bxse, byy, byse, extra)
    if lift_counts:
        res["lift_counts"] = json.dumps(lift_counts, sort_keys=True)
    p12s = cfg.get("p12s", [P12, P12_STRICT])
    if len(bx) >= MIN_COLOC_SNPS:
        for p12 in p12s:
            pp = E.coloc_abf(np.array(bx), np.array(bxse), np.array(byy), np.array(byse), p1=P1, p2=P2, p12=p12)
            tag = "" if p12 == p12s[0] else f"_p12_{p12:.0e}".replace("-0", "-")
            for h, v in pp.items():
                res[f"PP.{h}{tag}"] = v
            res[f"support{tag}"] = dual_rule(pp["H4"], pp["H3"])
        res["coloc_status"] = "ok"
    else:
        res["coloc_status"] = "TOO_FEW_SHARED_SNPS"
        res["support"] = "TOO_FEW_SHARED_SNPS"
    if "status" not in res:
        res["status"] = "ok"
    return res


def dual_rule(h4, h3):
    return bool(h4 >= 0.80 and (h4 / h3 if h3 > 0 else 9e9) >= 5)


def call_of(res, tag=""):
    """Human call: SUPPORTED / H3_distinct / H1_power / H0_H2 / TOO_FEW."""
    if res.get("coloc_status") != "ok":
        return "TOO_FEW_SHARED_SNPS"
    if res.get(f"support{tag}"):
        return "SUPPORTED"
    pp = {h: res.get(f"PP.{h}{tag}") for h in ("H0", "H1", "H2", "H3", "H4")}
    top = max(pp, key=lambda h: pp[h])
    return {"H3": "NOT_SUPPORTED_H3", "H1": "NOT_SUPPORTED_H1", "H4": "NOT_SUPPORTED_H4_below_rule",
            "H0": "NOT_SUPPORTED_H0", "H2": "NOT_SUPPORTED_H2"}[top]
