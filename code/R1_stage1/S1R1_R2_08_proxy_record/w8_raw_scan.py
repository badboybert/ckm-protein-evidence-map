# -*- coding: utf-8 -*-
"""w8_raw_scan.py (Stage 7 A2.2, S2-D47): read the outcome rows at the five pairs' leads and proxies from the RAW local outcome files.

Why: the saved extracts (review_R2/seen/wpb/fg_*.txt, mega.txt) hold the proxy rows, but not every lead (the ABO lead is absent from fg_T2D.txt). The raw files are on D: (already downloaded;
no download here), so the lead and proxy rows are read from them by POSITION, which also answers 'is the lead present under another label'. Extracts are cross-checked against these rows in
w8_proxy_record.py.  Internal output (provider values): w8_outcome_rows_raw.tsv is NOT a public file.
Positions are NOT typed: lead positions = instruments file chrpos (GRCh38); proxy positions = the position column of the saved extract row whose rsID is the proxy id (from w8_ld_rows.tsv).
Resumable: one block per raw file appended to the output with a marker; a re-run skips finished files.  python -I -X utf8 w8_raw_scan.py
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
import csv, gzip, hashlib, io, os, sys

V3 = _P31R_V3
HERE = V3 + "/R1_stage1/S1R1_R2_08_proxy_record"
RAW = (_P31R_EXT + "/D/data/paper31_v4_raw/ukbfree_outcomes")
WPB = V3 + "/R2_rewrite/review_R2/seen/wpb"
OUT = HERE + "/w8_outcome_rows_raw.tsv"
FILES = {"FinnGen_R10_T2D": "finngen_R10_T2D.gz", "FinnGen_R10_I9_CHD": "finngen_R10_I9_CHD.gz", "FinnGen_R10_I9_AF": "finngen_R10_I9_AF.gz",
         "MEGASTROKE_EUR_GCST006906": "29531354-GCST006906-EFO_0000712.h.tsv.gz"}
EXTRACT = {"FinnGen_R10_T2D": "fg_T2D.txt", "FinnGen_R10_I9_CHD": "fg_CHD.txt", "FinnGen_R10_I9_AF": "fg_AF.txt"}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


ins = {}
with io.open(V3 + "/V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv", encoding="utf-8", newline="") as fh:
    for r in csv.DictReader(fh):
        ins[(r["SNP"], r["gene"])] = r
ld = list(csv.DictReader(io.open(HERE + "/w8_ld_rows.tsv", encoding="utf-8", newline=""), delimiter="\t"))
assert len(ld) == 5
lead_pos = {}                                     # (chrom, pos) -> label
for r in ld:
    chrom, pos = ins[(r["lead_rsid"], r["pair"].split("-")[0])]["chrpos"].split(":")
    lead_pos[(chrom, int(pos))] = "lead:" + r["lead_rsid"]
prox_pos = {}
for oid, fn in EXTRACT.items():
    for line in io.open(WPB + "/" + fn, encoding="utf-8"):
        a = line.rstrip("\n").split("\t")
        for r in ld:
            if a[4] == r["proxy_rsid"]:
                prox_pos[(a[0], int(a[1]))] = "proxy:" + r["proxy_rsid"]
mega_hdr = open(WPB + "/mega.txt", encoding="utf-8").readline().rstrip("\n").split("\t")
for line in io.open(WPB + "/mega.txt", encoding="utf-8").read().splitlines()[1:]:
    a = dict(zip(mega_hdr, line.split("\t")))
    prox_pos[(a["hm_chrom"], int(a["hm_pos"]))] = "proxy:" + a["hm_rsid"]
want = dict(lead_pos); want.update(prox_pos)
print("positions wanted:", len(want), sorted(want.items()), flush=True)

done = set()
if os.path.exists(OUT):
    for line in io.open(OUT, encoding="utf-8"):
        if line.startswith("#DONE\t"):
            done.add(line.split("\t")[1].strip())
cols = ["outcome_id", "label", "chrom", "pos", "ref_or_other", "alt_or_effect", "rsid_label", "p", "beta", "se", "eaf", "extra"]
if not os.path.exists(OUT):
    io.open(OUT, "w", encoding="utf-8").write("\t".join(cols) + "\n")
for oid, fn in FILES.items():
    if oid in done:
        print("skip (done)", oid, flush=True)
        continue
    p = RAW + "/" + fn
    rows = []
    n = 0
    with gzip.open(p, "rt", encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if oid.startswith("FinnGen"):
            assert hdr[:13] == ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"], hdr
            ci, pi = hdr.index("#chrom"), hdr.index("pos")
        else:
            for need in ("hm_chrom", "hm_pos", "hm_rsid", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency"):
                assert need in hdr, (need, hdr)
            ci, pi = hdr.index("hm_chrom"), hdr.index("hm_pos")
        for line in f:
            n += 1
            a = line.rstrip("\n").split("\t")
            try:
                k = (a[ci], int(a[pi]))
            except ValueError:
                continue
            if k in want:
                d = dict(zip(hdr, a))
                if oid.startswith("FinnGen"):
                    rows.append([oid, want[k], k[0], k[1], d["ref"], d["alt"], d["rsids"], d["pval"], d["beta"], d["sebeta"], d["af_alt"], ""])
                else:
                    rows.append([oid, want[k], k[0], k[1], d["hm_other_allele"], d["hm_effect_allele"], d["hm_rsid"], d["p_value"], d["hm_beta"], d["standard_error"],
                                 d["hm_effect_allele_frequency"], "eff=%s;oth=%s;rs=%s" % (d["effect_allele"], d["other_allele"], d["variant_id"])])
    with io.open(OUT, "a", encoding="utf-8", newline="") as o:
        for r in rows:
            o.write("\t".join(str(x) for x in r) + "\n")
        o.write("#DONE\t%s\t%s\tsize=%d\tmd5=%s\tdata_rows=%d\n" % (oid, fn, os.path.getsize(p), md5(p), n))
    print("done", oid, "rows hit:", len(rows), "data rows:", n, flush=True)
print("ALL DONE", flush=True)
