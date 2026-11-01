#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 rewrite (Stage 3b) -- gate over the TABLE sources (tables/P31R_Table*.md) and the SUPPLEMENTARY METHODS
sources (supplement/P31R_SM_*.md, supplement/P31R_SNote.md). Curated supplementary tables written by hand (e.g. tables/P31R_TableS2.md, "## Table S2")
are gated the same way; data-export tables (S3-S8) are derived by scripts/P31R_build_af1.py directly from source files.

It applies the prose gate's own checks (P31R_gate.check: N1 untagged numbers, N2 value vs ledger at printed precision
with sign, N2b rounding boundary, N3 unknown tag, N4 spelled numbers, V1 banned vocabulary) -- one implementation,
imported, not forked -- to:
  * every table CELL, title line and footnote paragraph (each cell is checked as its own paragraph, so a tag in one
    cell can never vouch for a number in another);
  * every supplement paragraph (lists are allowed in the supplement; markdown tables are not -- use a table file).
Table structure checks:
  T1  the file starts with "## Table N", then one line "**Table N. <title>.**", then a markdown table, then optional
      footnote paragraphs; every row has the header's column count;
  T2  no cell is empty (write "--" or "not evaluable", never leave a blank a reader must interpret).
Run:  python P31R_tablegate.py [--selftest]      exit 1 on any failure.
D1 (added 2026-10-14; implemented in P31R_gate.py, imported): every ledger row named by a tag in a checked table / supplement file must have a
  display equal to its stored value rounded at the display's precision; failures print as "FAIL D1 Lxxx [MISMATCH|UNDERIVED|DATE] ...".
  --ledger FILE.csv reads another ledger (self-test: a COPY with one value changed).
Stage 5 step 1a (2026-10-16): the mutation-target scanner is P31R_gate.scan_pairs (one implementation, refuses a four-digit tag loudly); a tag in a table
  or supplement file is a citation for the D1 uncited-row exemption (P31R_gate.d1_tags(where=file)); selftest_s5 below exercises both.
R2 mode (added 2026-10-14; the default run above is unchanged and reads ../tables and ../supplement as before):
      python P31R_tablegate.py --sections DIR [--tables A.md,B.md] [--supp C.md,D.md] [--exemptions FILE]
      python P31R_tablegate.py --r2     = --sections ../sections_R2, table P31R_Table3_R2_draft.md, supplements P31R_SM_MethodsDetail_R2.md and
                                          P31R_SM_ResultsDetail_R2.md, and the extra exemptions in R2_numeral_exemptions.tsv
  Same checks, same rules. Observed on the earlier R2 Table 3 draft (header made conformant 2026-10-14): its first line was '## Table 3. <title>' (one line), not the two-line
  '## Table N' + '**Table N. <title>.**' form that T1 requires; its section rows still have empty cells (T2); the cell-level N/V checks do not
  run on a file that fails T1. That is a format question for the convert step, not something this option changes.
"""
import io, os, re, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TAB = os.path.join(ROOT, "tables")
SUP = os.path.join(ROOT, "supplement")
sys.path.insert(0, HERE)
import P31R_gate as PG  # noqa: E402

R2_TABLES = ["P31R_Table3_R2_draft.md"]
R2_SUPP = ["P31R_SM_MethodsDetail_R2.md", "P31R_SM_ResultsDetail_R2.md"]

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def _check_text(text, L, allow_lists=False):
    """Run P31R_gate.check on a block of markdown via a temporary file (no budget applies)."""
    fd, p = tempfile.mkstemp(suffix=".md", prefix="tg_")
    os.close(fd)
    try:
        io.open(p, "w", encoding="utf-8").write(text)
        fails, _ = PG.check(p, L)
    finally:
        os.remove(p)
    if allow_lists:
        fails = [f for f in fails if not f.startswith("S1 ")]
    return fails


def parse_table(raw):
    lines = raw.splitlines()
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    head = lines[i].strip() if i < len(lines) else ""
    m = re.fullmatch(r"## Table (S?\d+)", head)
    if not m:
        raise ValueError("T1 first line must be '## Table N', got %r" % head[:40])
    n = m.group(1)
    i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    title = lines[i].strip() if i < len(lines) else ""
    if not re.match(r"\*\*Table %s\. .+\.\*\*" % n, title):
        raise ValueError("T1 second line must be '**Table %s. <title>.**' (optionally followed by caption text), got %r" % (n, title[:60]))
    i += 1
    rows = []
    while i < len(lines) and not lines[i].strip().startswith("|"):
        i += 1
    while i < len(lines) and lines[i].strip().startswith("|"):
        cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
        rows.append(cells)
        i += 1
    if len(rows) < 3 or not all(re.fullmatch(r":?-{3,}:?", c) for c in rows[1]):
        raise ValueError("T1 a markdown table with a header row and a --- separator is required")
    header, body = rows[0], rows[2:]
    notes = "\n".join(lines[i:]).strip()
    return n, title, header, body, notes


def check_table(path, L):
    raw = io.open(path, encoding="utf-8").read()
    try:
        n, title, header, body, notes = parse_table(raw)
    except ValueError as e:
        return [str(e)], 0
    fails = []
    fails += ["title: " + f for f in _check_text(title, L)]
    for j, h in enumerate(header):
        fails += ["header col %d: %s" % (j + 1, f) for f in _check_text(h, L)]
    for r, row in enumerate(body, 1):
        if len(row) != len(header):
            fails.append("T1 row %d has %d cells, header has %d" % (r, len(row), len(header)))
        for c, cell in enumerate(row, 1):
            if not cell.strip():
                fails.append("T2 row %d col %d is empty" % (r, c))
                continue
            fails += ["row %d col %d (%s): %s" % (r, c, header[c - 1][:20] if c - 1 < len(header) else "?", f) for f in _check_text(cell, L)]
    if notes:
        fails += ["notes: " + f for f in _check_text(notes, L)]
    return fails, len(body)


def check_supp(path, L):
    raw = io.open(path, encoding="utf-8").read()
    fails = []
    if re.search(r"^\s*\|.*\|\s*$", raw, flags=re.M):
        fails.append("markdown table in a supplement section: put tables in a table file")
    fails += _check_text(raw, L, allow_lists=True)
    return fails


def _parse(argv):
    a = {"sections": None, "tables": None, "supp": None, "r2": False, "exemptions": None, "ledger": None}
    i = 0
    while i < len(argv):
        x = argv[i]
        if x in ("--sections", "--tables", "--supp", "--exemptions", "--ledger"):
            assert i + 1 < len(argv), "%s needs a value" % x
            a[x[2:]] = [f for f in argv[i + 1].split(",") if f] if x in ("--tables", "--supp") else os.path.abspath(argv[i + 1])
            i += 2
        elif x == "--r2":
            a["r2"] = True
            i += 1
        elif x == "--selftest":
            i += 1
        else:
            raise SystemExit("unknown option %s (see the docstring)" % x)
    return a


def main(argv=None):
    a = _parse(sys.argv[1:] if argv is None else argv)
    r2 = a["r2"]
    custom = r2 or any(a[k] for k in ("sections", "tables", "supp", "exemptions"))
    tabdir = supdir = a["sections"] or (PG.R2_SEC if r2 else None)
    tabdir, supdir = tabdir or TAB, supdir or SUP
    ex = a["exemptions"] or (PG.R2_EXEMPT_FILE if r2 else None)
    PG.EXTRA_EXEMPT[:] = PG.load_exemptions(ex) if ex else []
    old_ledger = PG.LEDGER
    if a["ledger"]:
        PG.LEDGER = a["ledger"]
    try:
        L = PG.ledger()
        total = 0
        tabs = a["tables"] or (R2_TABLES if r2 else (sorted(f for f in os.listdir(tabdir) if re.fullmatch(r"P31R_TableS?\d+\.md", f)) if os.path.isdir(tabdir) else []))
        sups = a["supp"] or (R2_SUPP if r2 else (sorted(f for f in os.listdir(supdir) if re.fullmatch(r"P31R_(SM_[a-z0-9_-]+|SNote)\.md", f)) if os.path.isdir(supdir) else []))
        if custom:
            print("MODE: dir=%s tables=%d supplements=%d exemption patterns=%d" % (tabdir, len(tabs), len(sups), len(PG.EXTRA_EXEMPT)))
        for f in tabs:
            fails, nrow = check_table(os.path.join(tabdir, f), L)
            fails = fails + PG.d1_tags(L, PG.d1_file_tags(os.path.join(tabdir, f)), where=f)
            print("%s: %d body rows, %d failure(s)" % (f, nrow, len(fails)))
            for x in fails:
                print("  FAIL " + x)
            total += len(fails)
        for f in sups:
            fails = check_supp(os.path.join(supdir, f), L) + PG.d1_tags(L, PG.d1_file_tags(os.path.join(supdir, f)), where=f)
            raw_ = io.open(os.path.join(supdir, f), encoding="utf-8").read()
            words = len(re.findall(r"\S+", PG.TAG.sub("", raw_)))
            print("%s: %d words, %d failure(s)" % (f, words, len(fails)))
            for x in fails:
                print("  FAIL " + x)
            if r2:
                npl = len(re.findall(r"\{\{[A-Za-z0-9_.]+\}\}", raw_))
                if npl:
                    print("  NOTE %d unresolved {{id}} placeholder(s) in %s (not a gate failure; convert them before this file counts as gated)" % (npl, f))
            total += len(fails)
        print("TABLE/SUPPLEMENT GATE %s: %d failure(s) over %d table(s) and %d supplement file(s)" % ("RED" if total else "GREEN", total, len(tabs), len(sups)))
        return 1 if total else 0
    finally:
        PG.EXTRA_EXEMPT[:] = []
        PG.LEDGER = old_ledger


def _cli(argv):
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = main(argv)
    out = buf.getvalue()
    return rc, [l[len("  FAIL "):] for l in out.splitlines() if l.startswith("  FAIL ")], out


def selftest_r2(L):
    """R2 mode through the command-line path: a mutation in a copy of the R2 SM Results draft must trip N2, N1 and N3 (supplement checks), the R2
    Table 3 draft must be opened from the directory it is pointed at, and --r2 must open all three R2 table/supplement files."""
    import collections, shutil
    bad = 0

    def rep(lab, ok, detail=""):
        nonlocal bad
        print("  R2 %-40s %s %s" % (lab, "ok" if ok else "MISS", detail))
        bad += not ok

    name = "P31R_SM_ResultsDetail_R2.md"
    d = tempfile.mkdtemp()
    for f in (name, R2_TABLES[0]):
        if not os.path.isfile(os.path.join(PG.R2_SEC, f)):
            rep("R2 draft present", False, f)
            return 1
        shutil.copy(os.path.join(PG.R2_SEC, f), os.path.join(d, f))
    cp = os.path.join(d, name)

    def fails_of(sec):
        return _cli(["--sections", sec, "--supp", name, "--exemptions", PG.R2_EXEMPT_FILE])

    rc0, base, _ = fails_of(d)
    rc1, real, _ = fails_of(PG.R2_SEC)
    rep("copy == real file (--sections --supp)", base == real and rc0 == rc1, "(%d baseline failure lines)" % len(base))
    text = io.open(cp, encoding="utf-8", newline="").read()
    cands = []
    for tok, tag, lid in PG.scan_pairs(text):
        pat = tok + " " + tag
        if text.count(pat) == 1 and not any(("'%s'" % tok) in x and lid in x for x in base):
            cands.append((tok, tag, lid, pat))
    rep("mutation targets exist", len(cands) >= 3, "(%d unique '<number> [[Lxxx]]' pairs)" % len(cands))

    def new_lines(mutated):
        io.open(cp, "w", encoding="utf-8", newline="").write(mutated)
        try:
            rc, now, _ = fails_of(d)
        finally:
            io.open(cp, "w", encoding="utf-8", newline="").write(text)
        return rc, list((collections.Counter(now) - collections.Counter(base)).elements())

    hit = {"N2": None, "N1": None, "N3": None}
    for tok, tag, lid, pat in cands[:60]:
        if hit["N2"] is None:
            mt = "7" + tok
            rc, nl = new_lines(text.replace(pat, mt + " " + tag, 1))
            if rc == 1 and any((x.startswith("N2 ") and ("'%s'" % mt) in x) for x in nl):
                hit["N2"] = pat
        if hit["N1"] is None:
            rc, nl = new_lines(text.replace(pat, tok, 1))
            if rc == 1 and any((x.startswith("N1 ") and ("'%s'" % tok) in x) for x in nl):
                hit["N1"] = pat
        if hit["N3"] is None:
            rc, nl = new_lines(text.replace(pat, tok + " [[L000]]", 1))
            if rc == 1 and any((x.startswith("N3 unknown tag L000")) for x in nl):
                hit["N3"] = pat
        if all(hit.values()):
            break
    for code, why in (("N2", "a wrong number next to a tag"), ("N1", "a tag removed"), ("N3", "an unknown tag")):
        rep("mutation in R2 copy trips %s" % code, hit[code] is not None, "(%s; target %s)" % (why, hit[code] or "none"))
    rc_t, f_t, out_t = _cli(["--sections", d, "--tables", R2_TABLES[0], "--exemptions", PG.R2_EXEMPT_FILE])
    rep("Table 3 R2 copy is opened from --sections", ("%s:" % R2_TABLES[0]) in out_t and "MODE:" in out_t)
    rc_r, f_r, out_r = _cli(["--r2"])
    rep("--r2 opens table + both SM drafts", all(("%s:" % f) in out_r for f in R2_TABLES + R2_SUPP))
    probe = os.path.join(d, "probe.md")
    io.open(probe, "w", encoding="utf-8").write("Rule 3(ii) uses 5,000 draws in Stage 1.\n")
    toks = lambda fs: sorted(re.findall(r"untagged number '([^']+)'", " ".join(fs)))
    wi = toks(_cli(["--sections", d, "--supp", "probe.md", "--exemptions", PG.R2_EXEMPT_FILE])[1])
    wo = toks(_cli(["--sections", d, "--supp", "probe.md"])[1])
    rep("exemptions reach the supplement checks", wi == ["5,000"] and len(wo) > 1, "(with: %s; without: %s)" % (wi, wo))
    return 1 if bad else 0


def selftest_d1():
    """D1 reaches table cells and supplement paragraphs: a ledger COPY with only L004's value changed must trip D1 on L004 in a table and in a
    supplement file; the unchanged copy must not."""
    import csv
    bad = 0
    d = tempfile.mkdtemp()
    with io.open(PG.LEDGER, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh)
        fields, rows = rd.fieldnames, [dict(r) for r in rd]

    def ledger_copy(name, value):
        p = os.path.join(d, name)
        with io.open(p, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
            w.writeheader()
            for r in rows:
                if r["id"] == "L004" and value is not None:
                    r = dict(r, value=value)
                w.writerow(r)
        return p

    io.open(os.path.join(d, "P31R_Table9.md"), "w", encoding="utf-8").write(
        "## Table 9\n\n**Table 9. Test table.**\n\n| Resource | N |\n|---|---|\n| Screen | 14,521 [[L004]] |\n")
    io.open(os.path.join(d, "P31R_SM_probe.md"), "w", encoding="utf-8").write("The screen tested 14,521 [[L004]] pairs.\n")
    same, mut = ledger_copy("ledger_same.csv", None), ledger_copy("ledger_mut.csv", "14522")
    for lab, args in (("table", ["--tables", "P31R_Table9.md"]), ("supplement", ["--supp", "P31R_SM_probe.md"])):
        rc0, f0, _ = _cli(["--sections", d, "--ledger", same] + args)
        rc1, f1, _ = _cli(["--sections", d, "--ledger", mut] + args)
        ok = rc0 == 0 and not f0 and rc1 == 1 and any(x.startswith("D1 L004 [MISMATCH]") for x in f1) and PG.LEDGER.endswith("S1R1_CLAIMS_LEDGER.csv")
        print("  D1 %-40s %s %s" % ("%s: value-only change trips, copy passes" % lab, "ok" if ok else "MISS", (f1[:1] or f0[:1] or [""])[0][:60]))
        bad += not ok
    return bad


def selftest_s5():
    """Stage 5 step 1a through the table / supplement command-line path: a four-digit tag in a table cell or a supplement paragraph fails loudly (N3),
    the shared scanner refuses it, and a tag of a D1-exempt uncited row (L032) in a table cell or supplement paragraph fails as [EXEMPT-CITED]."""
    bad = 0
    d = tempfile.mkdtemp()
    base = "## Table 9\n\n**Table 9. Test table.**\n\n| Resource | N |\n|---|---|\n| Screen | 14,521 [[%s]] |\n"
    cases = []
    for tag, want in (("L1000", "N3 four-digit ledger tag L1000"), ("L032", "D1 L032 [EXEMPT-CITED]"), ("L004", None)):
        io.open(os.path.join(d, "P31R_Table9.md"), "w", encoding="utf-8").write(base % tag)
        io.open(os.path.join(d, "P31R_SM_probe.md"), "w", encoding="utf-8").write("The screen tested 14,521 [[%s]] pairs.\n" % tag)
        for lab, args in (("table", ["--tables", "P31R_Table9.md"]), ("supplement", ["--supp", "P31R_SM_probe.md"])):
            rc, f, _ = _cli(["--sections", d] + args)
            ok = (rc == 0 and not f) if want is None else (rc == 1 and any(want in x for x in f))
            print("  S5 %-12s %-11s %-34s %s %s" % (lab, tag, want or "clean", "ok" if ok else "MISS", (f[:1] or [""])[0][:60]))
            bad += not ok
    try:
        PG.scan_pairs("so 7 [[L1000]]")
        loud = False
    except AssertionError:
        loud = True
    print("  S5 %-12s %-11s %-34s %s" % ("scan_pairs", "L1000", "refuses loudly", "ok" if loud else "MISS"))
    bad += not loud
    return bad


def selftest():
    L = PG.ledger()
    d = tempfile.mkdtemp()
    def tab(body):
        p = os.path.join(d, "P31R_Table9.md")
        io.open(p, "w", encoding="utf-8").write(body)
        return check_table(p, L)[0]
    ok = "## Table 9\n\n**Table 9. Test table.**\n\n| Resource | N |\n|---|---|\n| Screen | 14,521 [[L004]] |\n\nFootnote: tests 14,521 [[L004]].\n"
    cases = [("ok", tab(ok), None),
             ("untagged cell", tab(ok.replace("14,521 [[L004]] |", "14,521 |")), "N1"),
             ("wrong cell", tab(ok.replace("14,521 [[L004]] |", "14,522 [[L004]] |")), "N2"),
             ("tag in other cell does not vouch", tab(ok.replace("| Screen | 14,521 [[L004]] |", "| Screen 14,521 | [[L004]] |")), "N1"),
             ("empty cell", tab(ok.replace("| Screen | 14,521 [[L004]] |", "| Screen |  |")), "T2"),
             ("banned word", tab(ok.replace("Screen", "Screen endpoint")), "V1"),
             ("bad header", tab(ok.replace("## Table 9", "## Tab 9")), "T1"),
             ("ragged row", tab(ok.replace("| Screen | 14,521 [[L004]] |", "| Screen | 14,521 [[L004]] | x |")), "T1")]
    bad = 0
    for lab, out, want in cases:
        good = (not out) if want is None else any(want in x.split(":")[0] or x.startswith(want) or (" " + want + " ") in x or (": " + want) in x for x in out)
        print("  %-32s %s %s" % (lab, "ok" if good else "MISS", out[:1]))
        bad += not good
    print("SELFTEST %s (%d cases)" % ("PASS" if not bad else "FAIL", len(cases)))
    r2bad = selftest_r2(L)
    d1bad = selftest_d1()
    s5bad = selftest_s5()
    return 1 if (bad or r2bad or d1bad or s5bad) else 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv[1:] else main())
