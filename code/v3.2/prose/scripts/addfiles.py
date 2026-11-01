#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""addfiles.py -- the journal's "Additional file N" naming: the declared file map, the citation form, the tables
workbook (Additional file 1) and the listing section (PACKAGING PASS A, stage 2, 2026-09-30).

Cardiovascular Diabetology (BMC): additional files are named "Additional file 1" and so on, cited by that name in the
body, cited in sequence, and listed in a separate section of the manuscript text (file name, format, title,
description). The supplement docx does NOT hold every supplementary table: Tables S1-S11, S13 and S14 are separate
comma-separated files (its own inventory says so; 15 files, widths up to 54 columns), so ONE additional file cannot
carry "Table SN" for every N. The declared map:

    Additional file 1  the supplementary-tables workbook (.xlsx) built HERE from the canonical CSVs: Tables S1-S11,
                       S13 and S14, one worksheet per SOURCE FILE (13 tables in 15 worksheets: S14 has three), every
                       cell the source file's text; the Contents worksheet lists table, title, source file, rows and
                       the source's md5. Its creator and last-modified-by are the standing author (global CLAUDE.md
                       section 8; PASS A FIX ROUND -- openpyxl wrote 'openpyxl').
    Additional file 2  the supplement docx: Figures S1-S20 and Table S12 (reproduced in full there).

Numbering follows the journal's "cited in sequence": the first supplementary item the text cites is Table S1 (Study
Design, the STROBE-MR checklist), the first figure comes a paragraph later, so the tables file is 1. The scout brief
called the docx "Additional file 1" on the premise that it held every table; the premise was wrong (pre-audit), and
the docx is file 2.

The citation form in the MAIN text is "Additional file K: Table SN" / "Additional file K: Figure SN" (ranges and lists
kept as written: "Tables S5 and S6", "Figures S2-S3"); the S-numbers are unchanged (USER). The SOURCES -- INSERTS.md
blocks, the legend files, the generators, the rebase edits, the claim literals -- keep the legacy "Supplementary Table
SN" wording, and `journal_form` is applied ONCE, to the built main document (apply_blocks, before the sidecar
assertion) and at every compare site that matches a source literal against the built text. A cross-file table range
("Supplementary Tables S1-S14") has no journal form and REFUSES; the one such phrase (Availability) is rewritten by a
declared rebase edit before the transformation runs. The supplement's own internal wording (its headings, inventory
and cross-references) is untouched except its title and the sentence that said the tables were separate CSV files.

  python addfiles.py --self-test
  python addfiles.py --build <main.docx> <supplement.docx>    # workbook + listing section (build step 4c)
  python addfiles.py --check <main.docx> <supplement.docx>    # report only
"""
from __future__ import print_function
import csv
import hashlib
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROSE = os.path.dirname(HERE)
V32 = os.path.dirname(PROSE)
V3 = os.path.dirname(V32)
ROOT = os.path.dirname(V3)

WORKBOOK_DIR = os.path.join(PROSE, "additional_files")
WORKBOOK = os.path.join(WORKBOOK_DIR, "Paper3.1_Additional_file_1_Supplementary_Tables.xlsx")
# a gate may point the reader at a COPY (the mutation ladder); the build always writes WORKBOOK
def workbook_path():
    return os.environ.get("VP_ADDFILE1") or WORKBOOK


class AddfileError(Exception):
    pass


# ---------------------------------------------------------------- the declared map ------------------------------------
TABLES_IN_DOCX = {12}                     # reproduced in full at the end of the supplement docx
TABLE_NUMBERS = list(range(1, 15))        # S1..S14
FIGURE_NUMBERS = list(range(1, 21))       # S1..S20
FILES = {
    1: {"kind": "tables", "format": ".xlsx (Microsoft Excel workbook)",
        "title": u"Supplementary Tables S1–S11, S13 and S14",
        "tables": [n for n in TABLE_NUMBERS if n not in TABLES_IN_DOCX], "figures": []},
    2: {"kind": "docx", "format": ".docx (Microsoft Word)",
        "title": u"Supplementary Figures S1–S20 and Supplementary Table S12",
        "tables": sorted(TABLES_IN_DOCX), "figures": FIGURE_NUMBERS},
}
# the supplement's title paragraph (its p0, Heading 1) and the inventory sentence about the separate CSV files
SUPP_TITLE_OLD = u"Supplementary Figures and Tables"
SUPP_TITLE = u"Additional file 2 — Supplementary Figures S1–S20 and Supplementary Table S12"
SUPP_TABLES_OLD = (u"Supplementary Tables S1–S11, S13 and S14 are supplied as separate comma-separated files")
# PASS A FIX ROUND: 'one worksheet per table' was inexact -- 13 tables sit in 15 worksheets (S14 has three source files)
SUPP_TABLES_NEW = (u"Supplementary Tables S1–S11, S13 and S14 are supplied as the worksheets of Additional file 1, "
                   u"one worksheet per source file, each a verbatim copy of the comma-separated source file named below")


def file_of(kind, n):
    """kind 'Table' / 'Figure', n the S-number -> the additional file that holds it (or AddfileError)."""
    key = "tables" if kind == "Table" else "figures"
    hits = [k for k, f in FILES.items() if n in f[key]]
    if len(hits) != 1:
        raise AddfileError("%s S%d is held by %d additional files" % (kind, n, len(hits)))
    return hits[0]


# ---------------------------------------------------------------- the citation form -----------------------------------
_ITEM = r"S\d+[a-z]?"
_LIST = _ITEM + r"(?:\s*[–—-]\s*S?\d+[a-z]?)?(?:(?:,\s*|\s+and\s+)" + _ITEM + r")*"
LEGACY_RX = re.compile(r"Supplementary (Table|Figure)(s?) (" + _LIST + r")")
JOURNAL_RX = re.compile(r"Additional file (\d+): (Table|Figure)(s?) (" + _LIST + r")")
PREFIX_RX = re.compile(r"Additional files? \d+(?: and \d+)?: ")


def items_of(lst):
    """'S5 and S6' -> [5, 6]; 'S9-S12' -> [9, 10, 11, 12]; 'S18b' -> [18]."""
    m = re.match(r"S(\d+)[a-z]?\s*[–—-]\s*S?(\d+)", lst)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        if b < a:
            raise AddfileError("descending range %r" % lst)
        return list(range(a, b + 1))
    return [int(x) for x in re.findall(r"S(\d+)", lst)]


def journal_form(text):
    """Every legacy citation of `text` in the journal form; the rest of the text byte-identical. Refuses a citation whose
    items are not all in one additional file."""
    def sub(m):
        kind, plural, lst = m.group(1), m.group(2), m.group(3)
        files = sorted({file_of(kind, n) for n in items_of(lst)})
        if len(files) != 1:
            raise AddfileError("%r spans additional files %s: no journal form" % (m.group(0), files))
        return u"Additional file %d: %s%s %s" % (files[0], kind, plural, lst)
    return LEGACY_RX.sub(sub, text)


def cross_file_form(text, here):
    """Every legacy citation in `text` of items held by an additional file OTHER than `here`, in the journal form; a
    citation of `here`'s own items, and the rest of the text, byte-identical. PASS A FIX ROUND (2026-10-01): the
    supplement (Additional file 2) cites Tables S10 and S11, which Additional file 1 holds, as 'Supplementary Table SN'.
    Refuses a citation that spans files, as journal_form does."""
    def sub(m):
        kind, lst = m.group(1), m.group(3)
        files = sorted({file_of(kind, n) for n in items_of(lst)})
        if files == [here]:
            return m.group(0)
        if len(files) != 1:
            raise AddfileError("%r spans additional files %s: no journal form" % (m.group(0), files))
        return journal_form(m.group(0))
    return LEGACY_RX.sub(sub, text)


def legacy_citations(text):
    return [m.group(0) for m in LEGACY_RX.finditer(text)]


def cited(text):
    """[(file, kind, S-number)] of every journal-form citation in `text`, ranges and lists expanded, in order."""
    out = []
    for m in JOURNAL_RX.finditer(text):
        for n in items_of(m.group(4)):
            out.append((int(m.group(1)), m.group(2), n))
    return out


def mask_prefix(text):
    """'Additional file 2: Figure S6' -> 'Figure S6': for numeral scans (the file number is a citation, not a datum)."""
    return PREFIX_RX.sub("", text)


# ---------------------------------------------------------------- the workbook (Additional file 1) ---------------------
SUPPTABLE_DIR = os.path.join(V3, "V4", "source_data", "supp_tables")       # the canonical SuppTable_*.csv (V4)
SUPPTABLE_MANIFEST = os.path.join(SUPPTABLE_DIR, "report_source_manifest.csv")
V4LOCAL_MANIFEST = os.path.join(V3, "V4_local", "MANIFEST.tsv")
CONTENTS = "Contents"
WORKBOOK_AUTHOR = "Bertrand Chin-Ming Tan"      # the checker's literal (the writer uses finalise_docx.AUTHOR)
WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
INV_RX = re.compile(r"^Supplementary Table S(\d+) — (.+)$")
PAREN_RX = re.compile(r"\(([^()]*)\)")
FILE_RX = re.compile(r"([A-Za-z0-9_]+\.csv), ([\d,]+) rows")


def _md5(path):
    return hashlib.md5(io.open(path, "rb").read()).hexdigest()


def _registered_md5(basename):
    """The md5 the data layer registers for a source CSV: report_source_manifest.csv for the V4 SuppTable files,
    V4_local/MANIFEST.tsv for the S13/S14 files. Refuses an unregistered file."""
    for r in csv.DictReader(io.open(SUPPTABLE_MANIFEST, encoding="utf-8-sig", newline="")):
        if r["table"] == basename:
            return r["md5"], SUPPTABLE_MANIFEST
    for r in csv.DictReader(io.open(V4LOCAL_MANIFEST, encoding="utf-8-sig", newline=""), delimiter="\t"):
        if os.path.basename(r["path"]) == basename:
            return r["md5"], V4LOCAL_MANIFEST
    raise AddfileError("%s is registered in neither manifest" % basename)


def _source_path(basename):
    p = os.path.join(SUPPTABLE_DIR, basename)
    if os.path.isfile(p):
        return p
    sys.path.insert(0, HERE)
    import apply_supplement
    rel = apply_supplement.NEW_TABLE_FILES.get(basename)
    if rel:
        return os.path.join(ROOT, rel.replace("/", os.sep))
    raise AddfileError("no canonical copy of %s" % basename)


def supplement_inventory(supp_docx):
    """[(S-number, title, [(basename, rows), ...])] from the built supplement's table inventory paragraphs -- the ONE
    place the table titles and file names are declared (apply_supplement asserts the row counts against the files)."""
    import docx
    d = docx.Document(supp_docx)
    out = []
    for p in d.paragraphs[:40]:
        t = re.sub(r"\s+", " ", p.text).strip()
        m = INV_RX.match(t)
        if not m:
            continue
        # the title runs up to the first parenthetical that gives row counts (a title may carry its own parenthesis:
        # 'Transcript (eQTLGen whole-blood) colocalization annotation (SuppTable_transcript_coloc.csv, 62 rows)')
        pm = next((x for x in PAREN_RX.finditer(m.group(2)) if "rows" in x.group(1)), None)
        if pm is None:
            raise AddfileError("S%s: no '(..., N rows)' parenthetical in %r" % (m.group(1), m.group(2)[:80]))
        title = m.group(2)[:pm.start()].strip()
        files = [(b, int(n.replace(",", ""))) for b, n in FILE_RX.findall(pm.group(1))]
        if not files and int(m.group(1)) not in TABLES_IN_DOCX:
            raise AddfileError("S%s: no '<file>.csv, N rows' in %r" % (m.group(1), pm.group(1)[:80]))
        out.append((int(m.group(1)), title, files))
    if [n for n, _, _ in out] != TABLE_NUMBERS:
        raise AddfileError("the supplement inventory lists tables %s, expected S1..S14" % [n for n, _, _ in out])
    return out


def _cell(s):
    """A CSV string -> the cell value. EVERY value is stored as text: a number stored as a number comes back through
    Excel's 15-significant-digit serialisation ('2.2737367544350304e-15' -> '2.27373675443503e-15') and an integral
    float comes back as an integer ('88109.0' -> 88109), both measured on this data, so a numeric cell cannot be a
    verbatim copy of its source. Text is."""
    return s


def cell_text(v):
    """The canonical text of a cell value read back (the inverse of _cell)."""
    if v is None:
        return u""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    return u"%s" % v


def sheet_name(n, k, nfiles):
    if nfiles == 1:
        return "Table S%d" % n
    return "Table S%d (%s)" % (n, ("screen", "classes", "per-signal")[k])


def workbook_rows(supp_docx):
    """[(sheet, S-number, title, basename, path, rows, md5, manifest)] -- every source found, its md5 equal to the
    registered one, its row count equal to the inventory's."""
    rows = []
    for n, title, files in supp_docx if isinstance(supp_docx, list) else supplement_inventory(supp_docx):
        if n not in FILES[1]["tables"]:
            continue
        for k, (base, nrows) in enumerate(files):
            path = _source_path(base)
            md5, man = _registered_md5(base)
            got = _md5(path)
            if got != md5:
                raise AddfileError("%s is %s on disk, registered %s in %s" % (base, got[:8], md5[:8], os.path.basename(man)))
            with io.open(path, encoding="utf-8-sig", newline="") as fh:
                actual = sum(1 for _ in csv.reader(fh)) - 1
            if actual != nrows:
                raise AddfileError("%s: %d rows on disk, the inventory prints %d" % (base, actual, nrows))
            rows.append((sheet_name(n, k, len(files)), n, title, base, path, nrows, md5, os.path.basename(man)))
    return rows


def build_workbook(supp_docx, out=None):
    """Write Additional file 1. Returns the rows written."""
    from openpyxl import Workbook
    rows = workbook_rows(supp_docx)
    out = out or WORKBOOK
    if not os.path.isdir(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))
    wb = Workbook()
    ws = wb.active
    ws.title = CONTENTS
    ws.append(["Worksheet", "Table", "Title", "Source file", "Rows", "Source md5", "Registered in"])
    for sheet, n, title, base, path, nrows, md5, man in rows:
        ws.append([sheet, "Supplementary Table S%d" % n, title, base, nrows, md5, man])
    ws.append([])
    ws.append(["Each worksheet is a verbatim copy of its comma-separated source file (header row first). Every value "
               "is stored as text exactly as the source prints it, numbers included, so that no value is reformatted "
               "or rounded; convert a column with Excel's VALUE() where numeric operations are needed."])
    for sheet, n, title, base, path, nrows, md5, man in rows:
        w = wb.create_sheet(sheet)
        with io.open(path, encoding="utf-8-sig", newline="") as fh:
            for i, r in enumerate(csv.reader(fh)):
                w.append(r if i == 0 else [_cell(c) for c in r])
    # PASS A FIX ROUND (m1): the standing author, not the tool default 'openpyxl' (global CLAUDE.md section 8; the two
    # docx files carry it since 2026-09-26). check_workbook pins the literal, so the writer's constant cannot drift.
    from finalise_docx import AUTHOR
    wb.properties.creator = AUTHOR
    wb.properties.lastModifiedBy = AUTHOR
    wb.save(out)
    return rows


def check_workbook(supp_docx, path=None):
    """[problems]: the workbook on disk holds exactly the declared sheets and, cell for cell, the source files' text."""
    from openpyxl import load_workbook
    path = path or workbook_path()
    bad = []
    if not os.path.isfile(path):
        return ["Additional file 1 is not on disk: %s" % path]
    rows = workbook_rows(supp_docx)
    wb = load_workbook(path, read_only=True)
    # PASS A FIX ROUND (m1): the file's own metadata names the standing author (it read 'openpyxl')
    for field in ("creator", "lastModifiedBy"):
        got = getattr(wb.properties, field, None)
        if got != WORKBOOK_AUTHOR:
            bad.append("docProps %s is %r, not %r" % (field, got, WORKBOOK_AUTHOR))
    want = [CONTENTS] + [r[0] for r in rows]
    if wb.sheetnames != want:
        bad.append("worksheets %s, expected %s" % (wb.sheetnames, want))
    for sheet, n, title, base, spath, nrows, md5, man in rows:
        if sheet not in wb.sheetnames:
            continue
        with io.open(spath, encoding="utf-8-sig", newline="") as fh:
            src = list(csv.reader(fh))
        got = [[cell_text(v) for v in r] for r in wb[sheet].iter_rows(values_only=True)]
        # openpyxl drops trailing empty cells of a row and trailing empty rows
        src_n = [[c for c in r] for r in src]
        while got and all(c == "" for c in got[-1]):
            got.pop()
        if len(got) != len(src_n):
            bad.append("%s: %d rows in the workbook, %d in %s" % (sheet, len(got), len(src_n), base))
            continue
        for i, (a, b) in enumerate(zip(got, src_n)):
            a = a + [""] * (len(b) - len(a))
            if a[:len(b)] != b or any(c != "" for c in a[len(b):]):
                j = next((k for k in range(max(len(a), len(b))) if (a + [""])[k] != (b + [""])[k]), 0)
                bad.append("%s row %d col %d: workbook %r, %s %r" % (sheet, i, j, a[j] if j < len(a) else None, base,
                                                                      b[j] if j < len(b) else None))
                break
    # the Contents sheet names every source and its registered md5
    if CONTENTS in wb.sheetnames:
        ct = [[cell_text(v) for v in r] for r in wb[CONTENTS].iter_rows(values_only=True)]
        listed = {(r[3], r[5]) for r in ct[1:] if len(r) >= 6 and r[3].endswith(".csv")}
        for sheet, n, title, base, spath, nrows, md5, man in rows:
            if (base, md5) not in listed:
                bad.append("Contents does not list %s with md5 %s" % (base, md5[:8]))
    wb.close()
    return bad


# ---------------------------------------------------------------- the listing section (main text) -----------------------
LISTING_HEADING = u"Additional files"


def listing_paragraphs(supp_docx, workbook=None):
    """The journal's four fields per file, one paragraph per additional file, derived from the built supplement (its
    inventory, headings and S12 caption) and the workbook rows -- never typed."""
    import docx
    rows = workbook_rows(supp_docx)
    inv = supplement_inventory(supp_docx)
    titles = {n: t for n, t, _ in inv}
    d = docx.Document(supp_docx)
    figs = [m.group(1) for p in d.paragraphs
            for m in [re.match(r"Supplementary Figure (S\d+[a-c]?)\.", re.sub(r"\s+", " ", p.text).strip())] if m]
    nums = sorted({int(re.sub(r"\D", "", f)) for f in figs})
    if nums != FIGURE_NUMBERS:
        raise AddfileError("the supplement's figure headings run %s, expected S1..S20" % nums)
    s12 = [re.sub(r"\s+", " ", p.text).strip() for p in d.paragraphs if p.text.strip().startswith("Supplementary Table S12.")]
    if len(s12) != 1:
        raise AddfileError("expected one Supplementary Table S12 caption in the supplement, found %d" % len(s12))
    m12 = re.search(r"(\d+) (?:colocalization-supported )?protein–outcome pairs", s12[0]) or \
        re.search(r"\((\d+) rows", titles.get(12, ""))
    n12 = m12.group(1) if m12 else re.search(r"(\d+)", titles[12]).group(1)
    by_n = {}
    for sheet, n, title, base, path, nrows, md5, man in rows:
        by_n.setdefault(n, []).append((base, nrows))
    entries = []
    for n in sorted(by_n):
        fl = by_n[n]
        if len(fl) == 1:
            entries.append(u"Table S%d — %s (%s rows)" % (n, titles[n], "{:,}".format(fl[0][1])))
        else:
            entries.append(u"Table S%d — %s, in %s worksheets (%s)" % (
                n, titles[n], WORDS.get(len(fl), str(len(fl))),
                "; ".join("%s, %s rows" % (b, "{:,}".format(r)) for b, r in fl)))
    f1 = FILES[1]
    # PASS A FIX ROUND (m10): 'N worksheets, one per table' was inexact (13 tables in 15 worksheets); both counts derived
    p1 = (u"Additional file 1. File format: %s. Title: %s. Description: %d tables in %d worksheets, one per source "
          u"file, each a verbatim copy of that comma-separated file; the Contents worksheet names every table's source "
          u"file, row count and file checksum. %s." % (f1["format"], f1["title"], len(by_n), len(rows), "; ".join(entries)))
    f2 = FILES[2]
    split = sorted({f for f in figs if f[-1].isalpha()})
    stem = sorted({f[:-1] for f in split})
    panels = (u" (%d panels; %s as %s)" % (len(figs), " and ".join(stem), " and ".join(split))) if split else u""
    p2 = (u"Additional file 2. File format: %s. Title: %s. Description: the %d supplementary figures%s with their "
          u"legends, the inventory of all fourteen supplementary tables with their row counts and source file names, "
          u"and Supplementary Table S12 — %s, reproduced in full."
          % (f2["format"], f2["title"], len(nums), panels, titles[12]))
    return [p1, p2]


def insert_listing(doc, supp_docx):
    """Append the 'Additional files' section (Heading 1 + one paragraph per file) after the reference list. Refuses if a
    heading of that name already exists."""
    from docx_edit import para_text
    import copy
    paras = doc.paragraphs
    if any(re.sub(r"\s+", " ", para_text(p)).strip() == LISTING_HEADING and p.style.name.startswith("Heading")
           for p in paras):
        raise AddfileError("an '%s' heading is already present" % LISTING_HEADING)
    heads = [p for p in paras if re.sub(r"\s+", " ", para_text(p)).strip() == "References" and p.style.name == "Heading 1"]
    if len(heads) != 1:
        raise AddfileError("expected one 'References' Heading 1, found %d" % len(heads))
    body_paras = [p for p in paras if re.sub(r"\s+", " ", para_text(p)).strip()]
    last = body_paras[-1]
    if not re.match(r"^\d+\.\s", re.sub(r"\s+", " ", para_text(last)).strip()):
        raise AddfileError("the last paragraph is not a reference entry: %r" % para_text(last)[:60])
    # a plain body paragraph to clone for the entries: the first 'First Paragraph'/'Body Text' paragraph after Declarations
    body_like = next((p for p in paras if p.style.name in ("Body Text", "First Paragraph")), None)
    if body_like is None:
        raise AddfileError("no body paragraph to clone the listing style from")
    from docx.text.paragraph import Paragraph
    h_el = copy.deepcopy(heads[0]._p)
    last._p.addnext(h_el)
    h = Paragraph(h_el, last._parent)
    _set_text(h, LISTING_HEADING)
    anchor = h_el
    for txt in listing_paragraphs(supp_docx):
        el = copy.deepcopy(body_like._p)
        anchor.addnext(el)
        anchor = el
        _set_text(Paragraph(el, last._parent), txt)
    return 2


def _set_text(p, text):
    """Replace a cloned paragraph's text keeping its first run's formatting and dropping every other run and any
    drawing / section property it carried."""
    from docx.oxml.ns import qn
    pPr = p._p.find(qn("w:pPr"))
    if pPr is not None:
        for tag in ("w:sectPr", "w:pageBreakBefore"):
            for el in pPr.findall(qn(tag)):
                pPr.remove(el)
    runs = p.runs
    for r in runs[1:]:
        r._r.getparent().remove(r._r)
    if runs:
        runs[0].text = text
        for el in list(runs[0]._r):
            if el.tag not in (qn("w:rPr"), qn("w:t")):
                runs[0]._r.remove(el)
    else:
        p.add_run(text)


def listing_section(doc):
    """(heading paragraph index, [entry texts]) of the built main document, or (None, [])."""
    from docx_edit import para_text
    paras = doc.paragraphs
    hi = [i for i, p in enumerate(paras) if re.sub(r"\s+", " ", para_text(p)).strip() == LISTING_HEADING
          and p.style.name == "Heading 1"]
    if len(hi) != 1:
        return None, []
    out = []
    for p in paras[hi[0] + 1:]:
        t = re.sub(r"\s+", " ", para_text(p)).strip()
        if t:
            out.append(t)
    return hi[0], out


# ---------------------------------------------------------------- self-test --------------------------------------------
def self_test():
    ok = True

    def chk(label, cond):
        nonlocal ok
        print("  [%s] %s" % ("OK  " if cond else "FAIL", label))
        ok = ok and bool(cond)
    chk("one table -> file 1", journal_form(u"is undefined (Supplementary Table S2).") == u"is undefined (Additional file 1: Table S2).")
    chk("S12 -> file 2", journal_form(u"in Supplementary Table S12 and") == u"in Additional file 2: Table S12 and")
    chk("figure with panel letter", journal_form(u"shown in Supplementary Figure S18b);") == u"shown in Additional file 2: Figure S18b);")
    chk("figure range", journal_form(u"(Supplementary Figures S9–S12).") == u"(Additional file 2: Figures S9–S12).")
    chk("table list", journal_form(u"Supplementary Tables S5 and S6") == u"Additional file 1: Tables S5 and S6")
    chk("comma after the item is not a list", journal_form(u"(Supplementary Table S3, cystatin-C/BUN)") == u"(Additional file 1: Table S3, cystatin-C/BUN)")
    chk("lower-case generic word untouched", journal_form(u"a supplementary annotation") == u"a supplementary annotation")
    chk("already journal form is a fixed point", journal_form(u"Additional file 1: Table S3") == u"Additional file 1: Table S3")
    try:
        journal_form(u"Supplementary Tables S1–S14")
        chk("cross-file range REFUSES", False)
    except AddfileError:
        chk("cross-file range REFUSES", True)
    chk("cited() expands", cited(u"x (Additional file 2: Figures S13–S16) y Additional file 1: Tables S5 and S6") ==
        [(2, "Figure", 13), (2, "Figure", 14), (2, "Figure", 15), (2, "Figure", 16), (1, "Table", 5), (1, "Table", 6)])
    chk("legacy_citations finds the legacy form only",
        legacy_citations(u"Supplementary Table S1 and Additional file 1: Table S2") == [u"Supplementary Table S1"])
    chk("mask_prefix", mask_prefix(u"see Additional file 2: Figure S6 and Additional files 1 and 2: Tables S1") == u"see Figure S6 and Tables S1")
    chk("cross_file_form maps a file-1 table cited inside file 2",
        cross_file_form(u"edges in Supplementary Table S11. Edges", 2) == u"edges in Additional file 1: Table S11. Edges")
    chk("cross_file_form leaves file 2's own table and figures alone",
        cross_file_form(u"in Supplementary Table S12 and Supplementary Figure S6", 2)
        == u"in Supplementary Table S12 and Supplementary Figure S6")
    chk("_cell keeps every value as text", all(_cell(x) == x for x in ("0.0290", "0.029", "12", "88109.0", "NA", "")))
    chk("cell_text inverts _cell", all(cell_text(_cell(s)) == s for s in ("0.0290", "0.029", "12", "-3", "1e-05", "NA", "", "True")))
    print("addfiles self-test: %s" % ("OK" if ok else "FAIL"))
    return ok


def main():
    if "--self-test" in sys.argv:
        return 0 if self_test() else 1
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 2:
        raise SystemExit(__doc__)
    main_path, supp_path = args
    import docx
    if "--build" in sys.argv:
        rows = build_workbook(supp_path)
        print("Additional file 1 written: %s (%d worksheets + %s)" % (WORKBOOK, len(rows), CONTENTS))
        bad = check_workbook(supp_path, WORKBOOK)
        if bad:
            raise SystemExit("REFUSED: the written workbook does not equal its sources: %s" % bad[:3])
        print("  read back: every worksheet equals its source file cell for cell; every source md5 is the registered one")
        d = docx.Document(main_path)
        n = insert_listing(d, supp_path)
        d.save(main_path)
        print("listing section appended: '%s' + %d paragraphs" % (LISTING_HEADING, n))
        d = docx.Document(main_path)
        hi, entries = listing_section(d)
        if hi is None or len(entries) != len(FILES):
            raise SystemExit("REFUSED: the listing section did not land (%r, %d entries)" % (hi, len(entries)))
        for e in entries:
            print("  " + e[:110] + ("..." if len(e) > 110 else ""))
        return 0
    d = docx.Document(main_path)
    from docx_edit import para_text
    legacy = [c for p in d.paragraphs for c in legacy_citations(para_text(p))]
    cites = [c for p in d.paragraphs for c in cited(para_text(p))]
    hi, entries = listing_section(d)
    print("legacy citations: %d; journal-form items cited: %d; listing section: %s (%d entries)"
          % (len(legacy), len(cites), "present" if hi is not None else "ABSENT", len(entries)))
    print("workbook check: %s" % (check_workbook(supp_path) or "clean"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
