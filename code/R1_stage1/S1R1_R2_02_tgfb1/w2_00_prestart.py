"""w2_00_prestart.py -- W2 start-of-run checks (brief section 0, rules 1, 2, 3, 6):
governing md5s, the staging snapshot (start|end), tool versions, PLINK --version, gls_q.R --self-test.
  python w2_00_prestart.py start|end
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
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

RSCRIPT = _P31R_RSCRIPT


def snapshot(tag):
    rows = []
    for d, _, fs in os.walk(W.STAGING_FM):
        for f in fs:
            p = os.path.join(d, f).replace("\\", "/")
            rows.append((p, os.path.getsize(p), W.md5_file(p)))
    rows.sort()
    out = W.OUT + "/staging_snapshot_%s.tsv" % tag
    W.write_tsv(out, ["path", "bytes", "md5"], [{"path": a, "bytes": b, "md5": c} for a, b, c in rows])
    import hashlib
    dig = hashlib.md5("".join("%s\t%s\n" % (a, c) for a, _, c in rows).encode()).hexdigest()
    W.log("STAGING SNAPSHOT %s: %d files under %s, digest %s (%s)" % (tag, len(rows), W.STAGING_FM, dig, out))
    return rows, dig


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "start"
    W.log("==== w2_00_prestart %s ====" % tag)
    W.assert_governing()
    W.log("governing md5s asserted: " + ", ".join("%s=%s" % (k, m) for k, (_, m) in W.GOVERNING.items()))
    rows, dig = snapshot(tag)
    if tag == "end":
        a = W.read_tsv(W.OUT + "/staging_snapshot_start.tsv")
        b = W.read_tsv(W.OUT + "/staging_snapshot_end.tsv")
        same = [(r["path"], r["md5"]) for r in a] == [(r["path"], r["md5"]) for r in b]
        W.log("STAGING SNAPSHOT start == end: %s (%d vs %d files)" % (same, len(a), len(b)))
        if not same:
            W.log("ASSERT_FAILED staging snapshot changed during the run")
        return
    plink = W.res("PLINK19_paper6")
    v = subprocess.run([plink, "--version"], capture_output=True, text=True).stdout.strip()
    W.log("PLINK19_paper6 --version: %s (%s)" % (v, plink))
    import numpy, scipy, pandas
    W.log("Python %s; numpy %s; scipy %s; pandas %s" % (sys.version.split()[0], numpy.__version__, scipy.__version__, pandas.__version__))
    rv = subprocess.run([RSCRIPT, "-e", 'cat(R.version.string, "; coloc", as.character(packageVersion("coloc")), "; susieR", '
                         'as.character(packageVersion("susieR")), "; data.table", as.character(packageVersion("data.table")), "\\n")'],
                        capture_output=True, text=True)
    W.log("R: %s" % rv.stdout.strip())
    gls = W.GOVERNING["gls_q.R"][0]
    st = subprocess.run([RSCRIPT, gls, "--self-test"], capture_output=True, text=True)
    for line in (st.stdout + st.stderr).strip().splitlines():
        W.log("gls_q.R --self-test | " + line)
    if "SELF-TEST PASS" not in st.stdout:
        W.log("ASSERT_FAILED gls_q.R self-test did not pass")
    W.log("gls_q.R is not used by W2 (section 2 has no GLS Q: the TGFB1 contrast is a two-signal A-vs-B test); self-test run for rule 6 compliance")


if __name__ == "__main__":
    main()
