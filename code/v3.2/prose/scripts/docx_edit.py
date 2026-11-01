#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run-safe docx editing primitives for the Route-B prose phase.

Everything here works on `docx.text.paragraph.Paragraph` run lists, so character
formatting (bold/italic/superscript) survives a splice. Nothing here knows about
figure numbers; the mapping lives in CITATION_MAP.tsv and the callers.

Primitives
----------
para_text(p)                       raw concatenation of run text (what offsets index)
splice(p, edits)                   edits = [(start, end, new_text)], applied right-to-left
replace_paragraph_text(p, text)    keep the paragraph and its style, replace all its text
insert_paragraph_after(p, text, style)
insert_paragraph_before(p, text, style)
delete_paragraph(p)
copy_image_paragraph(src_p, dst_p) not used; see set_image
add_picture_paragraph(doc, after_p, image_path, width_emu, style)

Self-test: python docx_edit.py --self-test   (round-trips a real project docx)
"""
from __future__ import print_function
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import copy
import sys
import os

import docx
from docx.oxml.ns import qn
from docx.shared import Emu


class EditError(Exception):
    pass


def para_text(p):
    return u"".join(r.text or u"" for r in p.runs)


def _run_offsets(p):
    offs = []
    at = 0
    for r in p.runs:
        t = r.text or u""
        offs.append((at, at + len(t)))
        at += len(t)
    return offs


def splice(p, edits):
    """edits: iterable of (start, end, new_text) on para_text(p). Right-to-left."""
    edits = sorted(edits, key=lambda e: e[0], reverse=True)
    for (s, e, new) in edits:
        txt = para_text(p)
        if not (0 <= s <= e <= len(txt)):
            raise EditError("span (%d,%d) outside paragraph of len %d" % (s, e, len(txt)))
        offs = _run_offsets(p)
        if not offs:
            raise EditError("paragraph has no runs; cannot splice")
        i = j = None
        for k, (a, b) in enumerate(offs):
            if i is None and a <= s < b:
                i = k
            if a < e <= b:
                j = k
        if s == len(txt):
            i = len(offs) - 1
        if i is None:
            i = len(offs) - 1
        if j is None:
            j = len(offs) - 1
        runs = p.runs
        if i == j:
            a = offs[i][0]
            runs[i].text = (runs[i].text or u"")[: s - a] + new + (runs[i].text or u"")[e - a:]
        else:
            a_i = offs[i][0]
            a_j = offs[j][0]
            tail = (runs[j].text or u"")[e - a_j:]
            runs[i].text = (runs[i].text or u"")[: s - a_i] + new
            for k in range(i + 1, j):
                runs[k].text = u""
            runs[j].text = tail
    return p


def replace_paragraph_text(p, text):
    """Replace all text, keeping the paragraph, its style and the first run's font."""
    runs = p.runs
    if not runs:
        p.add_run(text)
        return p
    runs[0].text = text
    for r in runs[1:]:
        r.text = u""
    return p


def _new_p_like(p, text, style=None):
    new = copy.deepcopy(p._p)
    # strip every run and any drawing from the copy
    for child in list(new):
        if child.tag in (qn("w:r"), qn("w:hyperlink"), qn("w:bookmarkStart"), qn("w:bookmarkEnd")):
            new.remove(child)
    p._p.addnext(new)
    from docx.text.paragraph import Paragraph
    np = Paragraph(new, p._parent)
    if style is not None:
        np.style = style
    np.add_run(text)
    return np


def insert_paragraph_after(p, text, style=None):
    return _new_p_like(p, text, style)


def insert_paragraph_before(p, text, style=None):
    from docx.text.paragraph import Paragraph
    new = copy.deepcopy(p._p)
    for child in list(new):
        if child.tag in (qn("w:r"), qn("w:hyperlink"), qn("w:bookmarkStart"), qn("w:bookmarkEnd")):
            new.remove(child)
    p._p.addprevious(new)
    np = Paragraph(new, p._parent)
    if style is not None:
        np.style = style
    np.add_run(text)
    return np


def body_children(doc):
    """Body-level elements in document order: paragraphs AND tables.

    python-docx's `doc.paragraphs` hides tables, which is how a caption came to be
    separated from its own table by five printed pages: the insert was anchored on
    the caption PARAGRAPH, and the table object that follows it is invisible at
    that level.
    """
    from docx.table import Table
    from docx.text.paragraph import Paragraph
    out = []
    for ch in doc.element.body.iterchildren():
        if ch.tag == qn("w:p"):
            out.append(("p", Paragraph(ch, doc)))
        elif ch.tag == qn("w:tbl"):
            out.append(("tbl", Table(ch, doc)))
    return out


def table_after(doc, para, within=3):
    """The first table element following `para`, or None.

    `within` bounds how many body elements may sit between them, so a caption whose
    table is far away does not silently match a different table.
    """
    kids = list(doc.element.body.iterchildren())
    try:
        i = kids.index(para._p)
    except ValueError:
        raise EditError("paragraph is not a body child")
    for j in range(i + 1, min(i + 1 + within, len(kids))):
        if kids[j].tag == qn("w:tbl"):
            from docx.table import Table
            return Table(kids[j], doc)
    return None


def insert_paragraph_after_element(el, like_para, text, style=None):
    """Insert a new paragraph directly after an arbitrary body element.

    `like_para` supplies the paragraph properties to copy, as _new_p_like does.
    """
    from docx.text.paragraph import Paragraph
    new = copy.deepcopy(like_para._p)
    for child in list(new):
        if child.tag in (qn("w:r"), qn("w:hyperlink"), qn("w:bookmarkStart"),
                         qn("w:bookmarkEnd")):
            new.remove(child)
    el.addnext(new)
    np = Paragraph(new, like_para._parent)
    if style is not None:
        np.style = style
    np.add_run(text)
    return np


def delete_paragraph(p):
    el = p._p
    el.getparent().remove(el)
    el._p = el._element = None


def has_drawing(p):
    return len(p._p.findall(".//" + qn("a:blip"))) > 0


def set_image(doc, p, image_path, width_emu):
    """Replace the drawing in paragraph p with image_path at width_emu.

    Implemented as: clear the paragraph's runs, then add a fresh picture run.
    Keeps the paragraph object, its style and its position.
    """
    for r in list(p.runs):
        r._r.getparent().remove(r._r)
    run = p.add_run()
    run.add_picture(image_path, width=Emu(width_emu))
    return p


# ---------------------------------------------------------------- self-test

def _self_test():
    """Exercise every primitive on a REAL project docx, not a fixture we wrote."""
    root = _P31R_V3
    src = os.path.join(root, "Paper3.1_Draft_Manuscript_codex.docx")
    if not os.path.exists(src):
        print("SELF-TEST SKIP: %s not found" % src)
        return 1
    import tempfile
    ok = 0
    fail = 0

    def chk(name, cond, detail=""):
        nonlocal ok, fail
        if cond:
            ok += 1
            print("  PASS %s" % name)
        else:
            fail += 1
            print("  FAIL %s  %s" % (name, detail))

    d = docx.Document(src)
    # This project's docx is single-run-per-paragraph (histogram {0: 2, 1: 170}),
    # so T1 is the real path and T2 constructs a multi-run case for the arithmetic.
    ps = [p for p in d.paragraphs if len(p.runs) == 1 and len(para_text(p)) > 200]
    chk("T0 found a long paragraph to test on", len(ps) > 0)
    if not ps:
        return 1
    p = ps[0]
    before = para_text(p)

    # T1 splice inside one run
    i = before.find(u"the")
    splice(p, [(i, i + 3, u"THE")])
    after = para_text(p)
    chk("T1 single-run splice changes exactly those 3 chars",
        after[:i] == before[:i] and after[i:i + 3] == u"THE" and after[i + 3:] == before[i + 3:],
        repr(after[max(0, i - 10):i + 13]))

    # T1b two splices in one paragraph must not corrupt each other's offsets
    d1b = docx.Document(src)
    p1b = [q for q in d1b.paragraphs if para_text(q).count(u"the") >= 2 and len(para_text(q)) > 200][0]
    t1b = para_text(p1b)
    a1 = t1b.find(u"the")
    a2 = t1b.find(u"the", a1 + 3)
    splice(p1b, [(a1, a1 + 3, u"AAAAAAA"), (a2, a2 + 3, u"B")])
    want = t1b[:a1] + u"AAAAAAA" + t1b[a1 + 3:a2] + u"B" + t1b[a2 + 3:]
    chk("T1b two splices of different lengths in one paragraph", para_text(p1b) == want,
        repr(para_text(p1b)[:80]))

    # T2 splice spanning a run boundary (multi-run paragraph built from a real one)
    d2 = docx.Document(src)
    p2 = [q for q in d2.paragraphs if len(q.runs) == 1 and len(para_text(q)) > 200][0]
    whole = para_text(p2)
    half = len(whole) // 2
    p2.runs[0].text = whole[:half]
    p2.add_run(whole[half:])
    chk("T2a constructed a 2-run paragraph", len(p2.runs) == 2 and para_text(p2) == whole)
    splice(p2, [(half - 1, half + 1, u"##")])
    chk("T2b cross-run splice keeps length arithmetic",
        para_text(p2) == whole[: half - 1] + u"##" + whole[half + 1:],
        repr(para_text(p2)[half - 8:half + 8]))

    # T3 a splice that should be REFUSED
    d3 = docx.Document(src)
    p3 = [q for q in d3.paragraphs if len(para_text(q)) > 50][0]
    try:
        splice(p3, [(0, len(para_text(p3)) + 5, u"x")])
        chk("T3 out-of-range splice refuses", False, "it did not raise")
    except EditError:
        chk("T3 out-of-range splice refuses", True)

    # T4 insert / delete round-trip through a saved file
    d4 = docx.Document(src)
    n0 = len(d4.paragraphs)
    anchor = d4.paragraphs[60]
    np = insert_paragraph_after(anchor, u"ZZMARKERZZ", anchor.style)
    chk("T4a insert adds one paragraph", len(d4.paragraphs) == n0 + 1)
    tmp = os.path.join(tempfile.gettempdir(), "_docx_edit_selftest.docx")
    d4.save(tmp)
    d5 = docx.Document(tmp)
    found = [q for q in d5.paragraphs if para_text(q) == u"ZZMARKERZZ"]
    chk("T4b marker survives save and reload", len(found) == 1)
    if found:
        delete_paragraph(found[0])
        chk("T4c delete removes it", len(d5.paragraphs) == n0)
        d5.save(tmp)
        d6 = docx.Document(tmp)
        chk("T4d deletion survives save", not any(para_text(q) == u"ZZMARKERZZ" for q in d6.paragraphs))

    # T5 replace_paragraph_text keeps the style
    d7 = docx.Document(src)
    cap = [q for q in d7.paragraphs if q.style.name == "Caption"][0]
    st = cap.style.name
    replace_paragraph_text(cap, u"NEW CAPTION TEXT")
    chk("T5 replace keeps style and text", para_text(cap) == u"NEW CAPTION TEXT" and cap.style.name == st)

    # T6 MUTATION: a splice must actually change what a reader sees, not just the object
    d8 = docx.Document(src)
    p8 = [q for q in d8.paragraphs if u"colocalization" in para_text(q)][0]
    k = para_text(p8).find(u"colocalization")
    splice(p8, [(k, k + len(u"colocalization"), u"QQMUTATEDQQ")])
    tmp2 = os.path.join(tempfile.gettempdir(), "_docx_edit_selftest2.docx")
    d8.save(tmp2)
    d9 = docx.Document(tmp2)
    chk("T6 mutation is visible after a save+reload",
        any(u"QQMUTATEDQQ" in para_text(q) for q in d9.paragraphs))

    # T7 has_drawing finds the figure paragraphs and only those
    d10 = docx.Document(src)
    drawn = [i for i, q in enumerate(d10.paragraphs) if has_drawing(q)]
    chk("T7 has_drawing finds exactly the 7 figure paragraphs", len(drawn) == 7, str(drawn))

    print("\nself-test: %d passed, %d failed" % (ok, fail))
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(_self_test())
    print(__doc__)
