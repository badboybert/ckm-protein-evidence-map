"""rel_records.py -- the authors' own release records, shipped under out/release_records/ (RC2, S2-D35(4); audit step 6, minor 'files a stranger is pointed to').

The README and REPRODUCE point to these records, so they ship:
  * TIER_RULES.tsv         verbatim copy of inputs/TIER_RULES.tsv (the tier rules and the reasons);
  * PROVIDER_TERMS.tsv     verbatim copy of inputs/PROVIDER_TERMS.tsv (the page read for each provider and the wording it states);
  * CODE_INCLUDE.tsv       verbatim copy of inputs/CODE_INCLUDE.tsv (the code-bundle list);
  * LEAKSCAN_COINCIDENCES.tsv  the declared leak-scan coincidences REDACTED: file, line(s), provider and reason; never the value (the value column of
                           inputs/LEAKSCAN_COINCIDENCES.tsv is a withheld provider value for some rows, and the release gate searches the whole release for it).
Asserted here: the redacted file has no column named value, no decimal token of >= 4 significant digits in any cell, and exactly one row per row of the input.

Module interface (as rel_ledger): INPUTS, build(v3, out, withheld, ledger) -> records.  Runs after rel_portability and rel_ledger, so the line numbers refer to the files as shipped.
"""
import csv
import re
from pathlib import Path

REL_DIR = "R2_rewrite/release"
RECORDS_DIR = "release_records"
SRC = {"TIER_RULES.tsv": REL_DIR + "/inputs/TIER_RULES.tsv", "PROVIDER_TERMS.tsv": REL_DIR + "/inputs/PROVIDER_TERMS.tsv", "CODE_INCLUDE.tsv": REL_DIR + "/inputs/CODE_INCLUDE.tsv"}
LEAK_SRC = REL_DIR + "/inputs/LEAKSCAN_COINCIDENCES.tsv"
INPUTS = list(SRC.values()) + [LEAK_SRC]
RED_COLS = ["file", "line", "provider", "reason"]
DEC_RX = re.compile(r"(?<![\w.])(?:\d+\.\d+|\.\d+|\d+(?:\.\d+)?[eE][-+]?\d+)(?![\w])")
OT_RX = re.compile(r"(?i)opentargets|devstage|annotation")


class RecordsError(AssertionError):
    pass


def _req(c, m):
    if not c:
        raise RecordsError(m)


def _tsv(p):
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def _sig(tok):
    m = re.split(r"[eE]", tok)[0]
    if "." in m:
        m = m.rstrip("0")
    return len(m.replace(".", "").lstrip("0"))


def _same(a, b):
    """The same number, exactly or (both >= 6 significant digits) at 6 digits (the release gate's R14 comparison)."""
    fa, fb = float(a), float(b)
    return fa == fb or (_sig(a) >= 6 and _sig(b) >= 6 and float("%.6g" % fa) == float("%.6g" % fb))


def _providers(v3, srcs):
    """Providers of the withheld source files, read from TIER_RULES (file rows) and RELEASE_INPUTS (licence_class of the release path); a ledger row is a ChEMBL-masked row."""
    d = Path(v3) / REL_DIR / "inputs"
    by_t = {r["target"]: r["providers"] for r in _tsv(d / "TIER_RULES.tsv") if r["kind"] == "file" and r["providers"]}
    by_ri = {}
    for r in _tsv(d / "RELEASE_INPUTS.tsv"):
        lc = r["licence_class"]
        if r["release_path"] and lc.startswith("third_party_derived:"):
            by_ri.setdefault(r["release_path"], set()).update(x for x in lc.split(":", 1)[1].split("+") if x)
    res = set()
    for s in srcs:
        if s.startswith("ledger:"):
            res.add("ChEMBL (ledger row)")
        elif s in by_t:
            res.update(x for x in by_t[s].split(";") if x)
        elif OT_RX.search(s):
            res.update(("OpenTargets", "ChEMBL", "DGIdb"))
        elif s in by_ri:
            res.update(by_ri[s])
    return sorted(res)


def build(v3, out, withheld, ledger):
    v3, out = Path(v3), Path(out)
    d = out / RECORDS_DIR
    d.mkdir(parents=True, exist_ok=True)
    recs = []
    for name, src in SRC.items():
        _req((v3 / src).is_file(), "release record source %s is missing" % src)
        (d / name).write_bytes((v3 / src).read_bytes())
        recs.append(dict(release_path="%s/%s" % (RECORDS_DIR, name), ship_status="ship", source_paths=[src], producer="rel_records.py", licence_class="own_derived",
                         rows=len(_tsv(d / name)), columns=len(_tsv(d / name)[0]), note="release record: verbatim copy of the authors' build record %s" % src))
    rows = _tsv(v3 / LEAK_SRC)
    red = []
    for r in rows:
        f = r["out_path"]
        lines = []
        pth = out / f
        if pth.is_file() and not f.startswith(RECORDS_DIR + "/"):
            for i, ln in enumerate(pth.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
                if any(_same(x, r["value"]) for x in DEC_RX.findall(ln)):
                    lines.append(i)
        srcs = [s for s in r["withheld_source"].split(";") if s]
        reason = r["reason"]
        _req(not [x for x in DEC_RX.findall(reason) if _sig(x) >= 4], "a leak-scan reason carries a number: %s" % reason)
        red.append(dict(file=f, line=";".join(str(x) for x in lines[:8]) or "not found now", provider=";".join(_providers(v3, srcs)) or "not recorded", reason=reason))
    red.sort(key=lambda x: (x["file"], x["line"], x["reason"]))
    path = d / "LEAKSCAN_COINCIDENCES.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(RED_COLS)
        for x in red:
            w.writerow([x[c] for c in RED_COLS])
    got = _tsv(path)
    _req(len(got) == len(rows) and list(got[0].keys()) == RED_COLS and "value" not in got[0], "redacted leak-scan record: shape differs from the input")
    _req(not any(_sig(t) >= 4 for g in got for c in g.values() for t in DEC_RX.findall(c)), "redacted leak-scan record carries a number of >= 4 significant digits")
    recs.append(dict(release_path="%s/LEAKSCAN_COINCIDENCES.tsv" % RECORDS_DIR, ship_status="ship", source_paths=[LEAK_SRC], producer="rel_records.py", licence_class="own_derived",
                     rows=len(got), columns=len(RED_COLS), note="release record: the declared leak-scan coincidences by file, line, provider and reason, without the values"))
    return recs


def selftest(v3):
    import shutil
    import tempfile
    tmp = Path(tempfile.mkdtemp(prefix="rel_records_selftest_"))
    here = Path(v3)
    shutil.copytree(here / REL_DIR / "out", tmp / "out")
    recs = build(v3, tmp / "out", None, None)
    bad = []

    def rep(name, ok):
        print("  %-70s %s" % (name, "ok" if ok else "MISS"))
        if not ok:
            bad.append(name)
    red = _tsv(tmp / "out" / RECORDS_DIR / "LEAKSCAN_COINCIDENCES.tsv")
    vals = [r["value"] for r in _tsv(here / LEAK_SRC)]
    txt = (tmp / "out" / RECORDS_DIR / "LEAKSCAN_COINCIDENCES.tsv").read_text(encoding="utf-8")
    rep("four records, one row per declared coincidence", len(recs) == 4 and len(red) == len(vals))
    rep("no declared value is printed in the redacted record", not any(v in txt for v in vals if len(v) >= 6))
    rep("every row names file, line and provider", all(r["file"] and r["line"] and r["provider"] for r in red))
    print("rel_records selftest: %s" % ("PASS" if not bad else "FAIL %s" % bad))
    return not bad


if __name__ == "__main__":
    import sys
    sys.exit(0 if selftest(sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[2])) else 1)
