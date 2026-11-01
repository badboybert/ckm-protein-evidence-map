"""w2_common.py -- shared helpers for Stage 1 W2 (pre-declaration section 2: TGFB1 and the 2 x 2 table).

Every input path is taken from RESOURCES.json by id and its md5 is asserted at read time (brief section 0, rule 1).
Writes only inside OUT (rule 2). Every script appends timestamped lines to OUT/run.log.
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
import hashlib
import json
import os
import sys
import time

ROOT = _P31R_CKM
OUT = os.environ.get("W2_OUT", ROOT + "/V3/R1_stage1/S1R1_R2_02_tgfb1").replace("\\", "/")   # W2_OUT: re-run elsewhere (verifier L3)
RES_DIR = ROOT + "/V3/R1_stage1/S1R1_R2_00_resources"
STAGING_FM = ROOT + "/V3/V4_local/staging/l3-finemap"

GOVERNING = {
    "STAGE1_BRIEFS.md": (RES_DIR + "/STAGE1_BRIEFS.md", "40d0d52f7055542d774c03f1edf04e00"),
    "PRESPEC_REDACTED": (RES_DIR + "/P31R_R2_STAGE1_PRESPEC_REDACTED.md", "a281c3464619b1b67c28e13e7392c41c"),
    "PRESPEC_AMEND1_REDACTED": (RES_DIR + "/P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md", "f265d600973c43239bb8377e0cd39d33"),
    "PRESPEC_FULL": (ROOT + "/V3/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md", "3319dc39740d9f5b3b8822133d9108b7"),
    "PRESPEC_AMEND1_FULL": (ROOT + "/V3/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md", "377508832ea9bd384133e869c456c909"),
    "gls_q.R": (RES_DIR + "/gls_q.R", "b7390fdf5b6f2f9c1d30c298d6afdafc"),
    "RESOURCES.json": (RES_DIR + "/RESOURCES.json", "546b769ae93a8f566eec2b92f4be1575"),
}

SEED = 20261007
_MD5_CACHE = {}


def log(msg, also_print=True):
    line = "%s  %s" % (time.strftime("%Y-%m-%dT%H:%M:%S"), msg)
    with open(OUT + "/run.log", "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    if also_print:
        print(line, flush=True)


def md5_file(path, chunk=1 << 22):
    key = (path, os.path.getsize(path), os.path.getmtime(path))
    if key in _MD5_CACHE:
        return _MD5_CACHE[key]
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    _MD5_CACHE[key] = h.hexdigest()
    return _MD5_CACHE[key]


def assert_md5(path, expected, label):
    got = md5_file(path)
    if got != expected:
        log("MD5 MISMATCH %s: %s expected %s got %s -- STOP" % (label, path, expected, got))
        raise SystemExit("md5 mismatch for %s" % label)
    log("md5 ok %-40s %s  %s" % (label, got, path), also_print=False)
    return got


def assert_governing():
    for k, (p, m) in GOVERNING.items():
        assert_md5(p, m, k)


def resources():
    assert_md5(*GOVERNING["RESOURCES.json"], "RESOURCES.json")
    with open(GOVERNING["RESOURCES.json"][0], encoding="utf-8") as fh:
        d = json.load(fh)
    return {e["id"]: e for e in d}


_RES = None


def res(rid, check=True):
    """path of resource rid, md5 asserted at read time; refuses seen_output / not_found roles."""
    global _RES
    if _RES is None:
        _RES = resources()
    e = _RES[rid]
    if e.get("role") in ("seen_output", "full_text_after_results", "not_found"):
        raise SystemExit("resource %s has role %s: not an input before COMPARE_SEEN" % (rid, e.get("role")))
    p = e["path"]
    if check and e.get("md5") and os.path.isfile(p):
        assert_md5(p, e["md5"], rid)
    return p


def res_entry(rid):
    global _RES
    if _RES is None:
        _RES = resources()
    return _RES[rid]


def write_tsv(path, cols, rows):
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r.get(c, "")) for c in cols) + "\n")
    os.replace(tmp, path)


def _fmt(v):
    try:
        import numpy as _np
        if isinstance(v, _np.integer):
            return str(int(v))
        if isinstance(v, _np.bool_):
            return str(bool(v))
    except ImportError:
        pass
    if isinstance(v, float):                 # includes numpy float64 (a float subclass; its numpy-2 repr is not numeric)
        v = float(v)
        if v != v:
            return "NA"
        return repr(v)
    return str(v)


def read_tsv(path):
    import csv
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def is_pal(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def align(ea, oa, a1, a2):
    """+1 if (ea, oa) == (a1, a2); -1 if swapped; complement allowed only for non-palindromic SNPs; 0 otherwise."""
    if (ea, oa) == (a1, a2):
        return 1
    if (ea, oa) == (a2, a1):
        return -1
    if not is_pal(ea, oa):
        if (COMP.get(ea), COMP.get(oa)) == (a1, a2):
            return 1
        if (COMP.get(ea), COMP.get(oa)) == (a2, a1):
            return -1
    return 0


if __name__ == "__main__":
    sys.exit("helper module")
