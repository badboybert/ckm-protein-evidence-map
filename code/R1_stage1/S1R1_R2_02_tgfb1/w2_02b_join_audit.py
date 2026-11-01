"""w2_02b_join_audit.py -- hygiene pre-flight (skill gwas-sumstats-hygiene, scripts/join_audit.py) on the two TGFB1 locus
joins: UKB-PPP protein x FinnGen CHD (GRCh38 chr:pos:alleles) and UKB-PPP protein x ieu-a-7 (GRCh37 position from the
UKB-PPP ID + alleles). Writes keyed audit inputs to join_audit/ and logs the audit verbatim."""
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
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

AUD = os.path.expanduser((_P31R_EXT + "/tools/gwas-sumstats-hygiene/scripts/join_audit.py"))


def key(c, p, a, b):
    return "%s:%s:%s" % (c, p, ":".join(sorted((a.upper(), b.upper()))))


def main():
    d = W.OUT + "/join_audit"
    os.makedirs(d, exist_ok=True)
    W.log("==== w2_02b_join_audit ====")
    ins = [r for r in csv.DictReader(open(W.res("EXP_UKBPPP_INSTRUMENT_TABLE"), encoding="utf-8")) if r["gene"] == "TGFB1"][0]
    lo, hi = int(ins["cis_start"]), int(ins["cis_end"])
    P38, P37 = [], []
    for r in csv.DictReader(open(W.res("EXP_UKBPPP_EUR_TGFB1_OID20621"), encoding="utf-8"), delimiter="\t"):
        if not lo <= int(r["pos"]) <= hi:
            continue
        c, p37, ref, alt = r["ID"].split(":")[:4]
        base = {"chr": c, "ref": ref, "alt": alt, "beta": r["BETA"], "se": r["SE"], "eaf": r["A1FREQ"]}
        P38.append(dict(base, key=key(c, r["pos"], ref, alt), pos=r["pos"]))
        P37.append(dict(base, key=key(c, p37, ref, alt), pos=p37))
    FG = [{"key": key("19", r["pos"], r["ref"], r["alt"]), "chr": "19", "pos": r["pos"], "ref": r["ref"], "alt": r["alt"],
           "beta": r["beta"], "se": r["sebeta"], "pval": r["pval"], "eaf": r["af_alt"]}
          for r in csv.DictReader(open(W.OUT + "/extracts/OUT_FG_R10_I9_CHD__TGFB1.tsv", encoding="utf-8"), delimiter="\t")]
    C4 = []
    with open(W.OUT + "/extracts/OUT_CARDIOGRAMplusC4D_ieu-a-7__TGFB1.tsv", encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            x = dict(zip(hdr, line.rstrip("\n").split("\t")))
            fm = dict(zip(x["FORMAT"].split(":"), x["ieu-a-7"].split(":")))
            C4.append({"key": key("19", x["POS"], x["REF"], x["ALT"]), "chr": "19", "pos": x["POS"], "ref": x["REF"], "alt": x["ALT"],
                       "beta": fm["ES"], "se": fm["SE"], "eaf": fm.get("AF", "")})
    cols = ["key", "chr", "pos", "ref", "alt", "beta", "se", "eaf"]
    W.write_tsv(d + "/protein_grch38.tsv", cols, P38)
    W.write_tsv(d + "/protein_grch37.tsv", cols, P37)
    W.write_tsv(d + "/finngen_chd.tsv", cols + ["pval"], FG)
    W.write_tsv(d + "/ieu_a_7.tsv", cols, C4)
    ok = True
    for lab, e, o in (("P x FinnGen CHD (GRCh38)", "protein_grch38.tsv", "finngen_chd.tsv"), ("P x ieu-a-7 (GRCh37)", "protein_grch37.tsv", "ieu_a_7.tsv")):
        r = subprocess.run([sys.executable, "-X", "utf8", AUD, "--exposure", d + "/" + e, "--outcome", d + "/" + o, "--key", "key"],
                           capture_output=True, text=True)
        W.log("join_audit %s -> exit %d" % (lab, r.returncode))
        for line in (r.stdout + r.stderr).strip().splitlines():
            W.log("  join_audit | " + line)
        ok &= r.returncode == 0
    st = subprocess.run([sys.executable, "-X", "utf8", AUD, "--self-test"], capture_output=True, text=True)
    W.log("join_audit --self-test exit %d: %s" % (st.returncode, (st.stdout + st.stderr).strip().splitlines()[-1] if (st.stdout + st.stderr).strip() else ""))
    if not ok:
        W.log("join_audit reported a FAIL (see lines above); interpreted in README")


if __name__ == "__main__":
    main()
