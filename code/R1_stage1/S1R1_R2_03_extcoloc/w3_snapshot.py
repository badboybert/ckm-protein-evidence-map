#!/usr/bin/env python
"""W3 staging snapshot (brief §0 rule 2): md5 of every file under V4_local/staging/l3-finemap/.

usage: python -I w3_snapshot.py start|end
  start -> writes staging_snapshot_start.tsv
  end   -> writes staging_snapshot_end.tsv and asserts it is identical to the start snapshot
Both append a summary line to run.log.
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
import hashlib, os, sys, time

OUT = (_P31R_V3 + "/R1_stage1/S1R1_R2_03_extcoloc")
STAGING = (_P31R_V3 + "/V4_local/staging/l3-finemap")


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S") + " [snapshot] " + msg
    print(line, flush=True)
    with open(os.path.join(OUT, "run.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def snapshot():
    rows = []
    for root, _dirs, files in os.walk(STAGING):
        for f in files:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, STAGING).replace("\\", "/")
            rows.append((rel, os.path.getsize(p), md5(p)))
    rows.sort()
    return rows


def write(rows, path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("relpath\tbytes\tmd5\n")
        for r in rows:
            fh.write("%s\t%d\t%s\n" % r)


def digest(rows):
    return hashlib.md5("".join("%s|%d|%s\n" % r for r in rows).encode()).hexdigest()


if __name__ == "__main__":
    mode = sys.argv[1]
    assert mode in ("start", "end"), mode
    rows = snapshot()
    assert len(rows) > 0, "empty staging snapshot: reader is blind"
    path = os.path.join(OUT, "staging_snapshot_%s.tsv" % mode)
    write(rows, path)
    log("%s: %d files under %s; digest of (relpath|bytes|md5) list = %s" % (mode, len(rows), STAGING, digest(rows)))
    if mode == "end":
        start_path = os.path.join(OUT, "staging_snapshot_start.tsv")
        with open(start_path, encoding="utf-8") as fh:
            start = [tuple(l.rstrip("\n").split("\t")) for l in fh.readlines()[1:]]
        start = [(a, int(b), c) for a, b, c in start]
        if start != rows:
            diff = sorted(set(start) ^ set(rows))
            log("STAGING CHANGED between start and end: %d differing rows, first: %s" % (len(diff), diff[:5]))
            sys.exit(2)
        log("end snapshot IDENTICAL to start snapshot (%d files): nothing written under l3-finemap" % len(rows))
