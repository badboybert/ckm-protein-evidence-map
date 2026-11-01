#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 rewrite (R2) -- assemble Additional file 3 (the supplementary figures) into ONE PDF.

Scope is the ACTIVE registry (figures/af3/P31R_ACTIVE_FIGURES.json), never a glob: the AF3 figures in registry order, each as its page unit(s)
(FigureS5, FigureS5_p2 ... continuation pages are allowed and are never shrunk to fit) followed by its legend. Layout: A4, 20-mm margins, a
figure page is placed at its own size (170 mm wide, no scaling), 9-pt Arial legend under the last page of the figure, on a continuation page
when it does not fit. Title page and running footer use the 7/9/11-pt house sizes only.

Refuses (exit 1) unless the figure gate is GREEN for every AF3 unit: G0-G9 and G11 (blocking), G12 source manifest, G13 staleness and
panels. NOT G1b (decision F9: the 225-mm figure-plus-legend rule is a main-text print rule). The pointer gate (G14) is a release check over the
whole text and is run with `python scripts/P31R_figgate.py --pointers`, not here.

After saving, the PDF is READ BACK with two readers other than the builder's own layout code (PyMuPDF text layer via pdftext, and pdfminer):
characters per page inside a two-sided band and the two readers within 15%; every legend title recovered; for every unit, text lines known to be
in the unit PDF recovered on the assembled page by pdfminer and drawn at unit scale (width ratio 1.000 +/- 0.005, so nothing was shrunk);
author metadata. The row it appends to af3/AF3_assembly_log.tsv pins the assembled md5 and every input md5 (unit PDFs, legends, active
registry, title page, this script, figgate).

    python -I P31R_af3_assemble.py            build (refuses on a red gate), read back, log
    python -I P31R_af3_assemble.py --verify   re-check the shipped PDF against the log and the files on disk now (stale input = FAIL)
    python -I P31R_af3_assemble.py --selftest sandbox build with fixture units, then mutations that must each trip THEIR check
"""
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import contextlib, csv, hashlib, io, json, os, re, shutil, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(FIGS, "scripts"))
import P31R_figgate as FG  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

AUTHOR = "Bertrand Chin-Ming Tan"
A4 = (595.276, 841.89)
MARGIN = 56.693                      # 20 mm
FOOT_Y = A4[1] - 32.0
FONTDIR = _P31R_FONTDIR
REGULAR, BOLD = os.path.join(FONTDIR, "arial.ttf"), os.path.join(FONTDIR, "arialbd.ttf")
BODY, TAG, FOOT = 9.0, 11.0, 7.0     # the three house sizes
PITCH = 1.28 * BODY
LOG_COLS = ["built", "pdf", "pdf_md5", "pages", "chars_pymupdf", "chars_pdfminer", "units", "legends", "active_md5", "title_md5",
            "assembler_md5", "figgate_md5", "gate"]
CHAR_BAND_PER_PAGE = (60, 9000)      # two-sided; calibrated against the second reader, re-derive when either reader changes
DEFAULT_CFG = dict(out=os.path.join(HERE, "out", "Additional_file_3.pdf"), log=os.path.join(HERE, "AF3_assembly_log.tsv"),
                   title_page=os.path.join(os.path.dirname(FIGS), "sections", "P31R_TitlePage.md"))


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


class Refuse(RuntimeError):
    pass


# ---------------------------------------------------------------- inputs
def af3_plan(active):
    """[(figure entry, [unit ids in page order])] for the AF3 figures, in registry order."""
    units = FG.active_units(active)
    return [(f, [u for u, e, _ in units if e["id"] == f["id"]]) for f in active["figures"] if f["kind"] == "af3"]


def manuscript_title(path):
    """Read, never typed: the working title on the title page (a hand-copied title drifts)."""
    m = re.search(r"\*\*Title[^*]*\*\*\s*(.+)", io.open(path, encoding="utf-8").read())
    if not m:
        raise Refuse("manuscript title not found in %s" % path)
    return m.group(1).strip()


def parse_legend(txt):
    """(lead, body): the bold opener '**Figure S1. Title.**' and the rest, ledger tags, emphasis and heading removed."""
    t = re.sub(r"\s?\[\[[^\]]*\]\]", "", txt)
    t = "\n".join(l for l in t.splitlines() if not l.startswith("#")).strip()
    m = re.match(r"\*\*(Figure S\d+\.[^*]*)\*\*\s*(.*)", t, flags=re.S)
    if not m:
        raise Refuse("legend does not open with a bold '**Figure S<n>. Title.**'")
    return re.sub(r"\s+", " ", m.group(1)).strip(), re.sub(r"\s+", " ", m.group(2).replace("**", "")).strip()


# ---------------------------------------------------------------- layout
class Pager:
    def __init__(self, fitz, doc):
        self.fitz, self.doc = fitz, doc
        self.font = {REGULAR: fitz.Font(fontfile=REGULAR), BOLD: fitz.Font(fontfile=BOLD)}
        self.page = None
        self.y = MARGIN

    def new_page(self):
        self.page = self.doc.new_page(width=A4[0], height=A4[1])
        for f in (REGULAR, BOLD):
            self.page.insert_font(fontname="ar" if f == REGULAR else "ab", fontfile=f)
        self.y = MARGIN
        return self.page

    def width(self, text, size, fontfile=REGULAR):
        return self.font[fontfile].text_length(text, fontsize=size)

    def wrap(self, text, size, width, fontfile=REGULAR):
        lines, cur = [], ""
        for w in text.split(" "):
            t = (cur + " " + w).strip()
            if cur and self.width(t, size, fontfile) > width:
                lines.append(cur)
                cur = w
            else:
                cur = t
        if cur:
            lines.append(cur)
        return lines

    def line(self, text, size=BODY, bold=False, x=MARGIN):
        self.page.insert_text((x, self.y + size), text, fontname="ab" if bold else "ar", fontsize=size)
        self.y += PITCH if size <= BODY else 1.3 * size

    def room(self, need):
        return self.y + need <= A4[1] - MARGIN - 14


def title_page(P, title, plan, fig_pages, legend_text):
    P.new_page()
    P.line("Additional file 3", TAG, True)
    P.y += 4
    for l in P.wrap("Supplementary figures to: " + title, BODY, A4[0] - 2 * MARGIN):
        P.line(l)
    P.y += 8
    P.line("Contents", BODY, True)
    for f, units in plan:
        label = parse_legend(legend_text[f["id"]])[0]      # the legend's own title: the contents cannot drift from it
        pg = fig_pages[f["id"]]
        lines = P.wrap(label, BODY, A4[0] - 2 * MARGIN - 40)
        for i, l in enumerate(lines):
            P.line(l)
            if i == 0:
                P.page.insert_text((A4[0] - MARGIN - 20, P.y - PITCH + BODY), str(pg), fontname="ar", fontsize=BODY)


def assemble(plan, title, unit_pdf, legend_text, out_path):
    import fitz
    # two passes: the first fixes where each figure starts so the contents page can name page numbers; the layout is deterministic
    fig_pages = {f["id"]: 0 for f, _ in plan}
    for _ in range(2):
        doc = fitz.open()
        P = Pager(fitz, doc)
        title_page(P, title, plan, fig_pages, legend_text)
        starts = {}
        for f, units in plan:
            lead, body = parse_legend(legend_text[f["id"]])
            for k, u in enumerate(units):
                src = fitz.open(unit_pdf[u])
                w, h = src[0].rect.width, src[0].rect.height
                P.new_page()
                starts.setdefault(f["id"], doc.page_count)
                if k > 0:
                    P.line("%s (continued, page %d of %d)" % (f["id"].replace("FigureS", "Figure S"), k + 1, len(units)), BODY, True)
                    P.y += 2
                P.page.show_pdf_page(fitz.Rect(MARGIN, P.y, MARGIN + w, P.y + h), src, 0)
                P.y += h + 10
                src.close()
            # legend under the last page of the figure, on a continuation page when it does not fit
            width = A4[0] - 2 * MARGIN
            lines_l = P.wrap(lead, BODY, width, BOLD)
            lines_b = P.wrap(body, BODY, width)
            first = True
            for ls, bold in [(lines_l, True), (lines_b, False)]:
                for l in ls:
                    if not P.room(PITCH):
                        P.new_page()
                        P.line("%s (legend, continued)" % f["id"].replace("FigureS", "Figure S"), BODY, True)
                        P.y += 2
                    P.line(l, BODY, bold)
        fig_pages = starts
        n = doc.page_count
        for i, pg in enumerate(doc):
            if i:
                pg.insert_text((MARGIN, FOOT_Y), "Additional file 3  |  page %d of %d" % (i + 1, n), fontname="ar", fontsize=FOOT)
        if _ == 0:
            doc.close()
    doc.set_metadata({"author": AUTHOR, "creator": AUTHOR, "title": "Additional file 3", "subject": "", "keywords": "", "producer": doc.metadata.get("producer", "")})
    try:
        doc.subset_fonts()
    except Exception:
        pass
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp = out_path + ".tmp"
    doc.save(tmp, garbage=4, deflate=True)
    doc.close()
    os.replace(tmp, out_path)
    return fig_pages


# ---------------------------------------------------------------- read back
def _norm(s):
    return re.sub(r"\s+", "", s)


def readback(out_path, plan, unit_pdf, legend_text, expect_pages=None):
    """Independent reading of the shipped PDF. Returns (failures, stats). Never returns an empty reader result as clean."""
    import fitz
    from pdfminer.high_level import extract_text
    import pdftext
    f = []
    doc = fitz.open(out_path)
    n = doc.page_count
    if expect_pages is not None and n != expect_pages:
        f.append("R1 %d pages, the log says %d" % (n, expect_pages))
    md = doc.metadata
    if md.get("author") != AUTHOR or md.get("creator") != AUTHOR:
        f.append("R2 author/creator metadata is %r / %r, not %r" % (md.get("author"), md.get("creator"), AUTHOR))
    t1 = pdftext.pdf_text(out_path)                                   # raises when nothing is recovered
    pages_pm = [extract_text(out_path, page_numbers=[i]) for i in range(n)]
    c1, c2 = len(_norm(t1)), sum(len(_norm(p)) for p in pages_pm)
    if c2 == 0 or abs(c1 - c2) > 0.15 * max(c1, c2):
        f.append("R3 readers disagree: PyMuPDF/pdftext %d chars vs pdfminer %d" % (c1, c2))
    lo, hi = CHAR_BAND_PER_PAGE
    for i, p in enumerate(pages_pm):
        c = len(_norm(p))
        if not (lo <= c <= hi):
            f.append("R4 page %d: %d characters outside the plausible band %d-%d" % (i + 1, c, lo, hi))
    allpm = _norm("".join(pages_pm))
    for fe, units in plan:
        lead, _ = parse_legend(legend_text[fe["id"]])
        if _norm(lead) not in allpm:
            f.append("R5 legend title of %s not recovered by pdfminer: %r" % (fe["id"], lead[:60]))
        for u in units:
            ex = FG.extract(unit_pdf[u])
            cand = sorted((l["text"] for l in ex["lines"] if len(l["text"]) >= 12 and re.fullmatch(r"[A-Za-z0-9 \u00a0,.;:()/\u2013-]+", l["text"])), key=len, reverse=True)[:3]
            if not cand:
                f.append("R6 %s: no plain text line to probe the assembled page with" % u)
            for c in cand:
                page_hit = [i for i, p in enumerate(pages_pm) if _norm(c) in _norm(p)]
                if not page_hit:
                    f.append("R6 %s: line %r of the unit PDF not recovered on any assembled page" % (u, c[:50]))
                    continue
                # drawn at unit scale: width of the same string in the unit and on the assembled page
                src = fitz.open(unit_pdf[u])
                ru = src[0].search_for(c)
                ra = doc[page_hit[0]].search_for(c)
                src.close()
                if ru and ra:
                    ratio = ra[0].width / ru[0].width
                    if abs(ratio - 1.0) > 0.005:
                        f.append("R7 %s: line %r drawn at %.3f x its size in the unit PDF (a figure page must not be scaled)" % (u, c[:40], ratio))
                else:
                    f.append("R7 %s: line %r not located by the PyMuPDF search on both PDFs" % (u, c[:40]))
    doc.close()
    return f, dict(pages=n, chars_pymupdf=c1, chars_pdfminer=c2)


# ---------------------------------------------------------------- build / verify
def gate_units(unit_ids):
    """The figure gate over exactly these units. Returns (rc, captured output)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = FG.run(unit_ids)
    return rc, buf.getvalue()


def inputs(cfg):
    active = FG.load_active()
    errs = FG.active_errors(active)
    if errs:
        raise Refuse("active registry invalid: " + "; ".join(errs))
    plan = af3_plan(active)
    if not plan:
        raise Refuse("no AF3 figure in the active registry")
    unit_pdf = {u: os.path.join(FG.OUT, u + ".pdf") for _, us in plan for u in us}
    legend_text = {}
    for f, _ in plan:
        lp = FG.v3_path(f["legend"])
        if not os.path.exists(lp):
            raise Refuse("legend missing: %s" % f["legend"])
        legend_text[f["id"]] = io.open(lp, encoding="utf-8").read()
    return active, plan, unit_pdf, legend_text


def log_row(cfg, active_path, plan, unit_pdf, legend_text, stats, gate_out):
    return dict(built=time.strftime("%Y-%m-%d %H:%M:%S"), pdf=os.path.basename(cfg["out"]), pdf_md5=md5(cfg["out"]), pages=stats["pages"],
                chars_pymupdf=stats["chars_pymupdf"], chars_pdfminer=stats["chars_pdfminer"],
                units=";".join("%s=%s" % (u, md5(p)) for u, p in unit_pdf.items()),
                legends=";".join("%s=%s" % (f["id"], md5(FG.v3_path(f["legend"]))) for f, _ in plan),
                active_md5=md5(active_path), title_md5=md5(cfg["title_page"]), assembler_md5=md5(os.path.abspath(__file__)),
                figgate_md5=md5(FG.__file__), gate="GREEN: %s" % ", ".join(u for u in unit_pdf))


def build(cfg=None):
    cfg = dict(DEFAULT_CFG, **(cfg or {}))
    active, plan, unit_pdf, legend_text = inputs(cfg)
    rc, out = gate_units(list(unit_pdf))
    if rc != 0:
        print(out)
        raise Refuse("the figure gate is not green for the AF3 units: not assembling")
    title = manuscript_title(cfg["title_page"])
    assemble(plan, title, unit_pdf, legend_text, cfg["out"])
    fails, stats = readback(cfg["out"], plan, unit_pdf, legend_text)
    if fails:
        raise Refuse("read-back of the assembled PDF failed:\n  " + "\n  ".join(fails))
    row = log_row(cfg, FG.ACTIVE_JSON, plan, unit_pdf, legend_text, stats, out)
    new = not os.path.exists(cfg["log"])
    with io.open(cfg["log"], "a", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, LOG_COLS, delimiter="\t", lineterminator="\n")
        if new:
            w.writeheader()
        w.writerow(row)
    print("assembled %s: %d pages, md5 %s; read-back clean (%d / %d chars)" % (cfg["out"], stats["pages"], row["pdf_md5"], stats["chars_pymupdf"], stats["chars_pdfminer"]))
    return row


def verify(cfg=None):
    cfg = dict(DEFAULT_CFG, **(cfg or {}))
    f = []
    if not os.path.exists(cfg["log"]) or not os.path.exists(cfg["out"]):
        return ["V0 no assembled PDF or no assembly log: run the build"]
    rows = list(csv.DictReader(io.open(cfg["log"], encoding="utf-8", newline=""), delimiter="\t"))
    row = rows[-1]
    active, plan, unit_pdf, legend_text = inputs(cfg)
    if md5(cfg["out"]) != row["pdf_md5"]:
        f.append("V1 the shipped PDF's md5 is not the logged one")
    if {u: m for u, m in (x.split("=") for x in row["units"].split(";"))} != {u: md5(p) for u, p in unit_pdf.items()}:
        f.append("V2 the unit set or a unit PDF changed since the assembly (units: %s)" % ", ".join(unit_pdf))
    if row["legends"] != ";".join("%s=%s" % (fe["id"], md5(FG.v3_path(fe["legend"]))) for fe, _ in plan):
        f.append("V3 a legend (or the AF3 figure list) changed since the assembly")
    for k, p in (("active_md5", FG.ACTIVE_JSON), ("title_md5", cfg["title_page"]), ("assembler_md5", os.path.abspath(__file__)), ("figgate_md5", FG.__file__)):
        if md5(p) != row[k]:
            f.append("V4 %s changed since the assembly: %s" % (k, p))
    t = os.path.getmtime(cfg["out"])
    late = [u for u, p in unit_pdf.items() if os.path.getmtime(p) > t + FG.STALE_TOL_S]
    if late:
        f.append("V5 unit PDF(s) newer than the assembled PDF: %s" % ", ".join(late))
    rb, _ = readback(cfg["out"], plan, unit_pdf, legend_text, int(row["pages"]))
    return f + rb


# ---------------------------------------------------------------- self-test
def _fixture_unit(path, tag_letters, lines, w_mm=170.0, h_mm=120.0):
    """A cairo-like single-page unit: Arial text at 7/9/11 pt, panel tags 11 pt bold (flag 16). No digits, so G7 has nothing to claim."""
    import fitz
    d = fitz.open()
    pg = d.new_page(width=w_mm / 25.4 * 72, height=h_mm / 25.4 * 72)
    pg.insert_font(fontname="ar", fontfile=REGULAR)
    pg.insert_font(fontname="ab", fontfile=BOLD)
    y = 20
    for t in tag_letters:
        pg.insert_text((10, y), t, fontname="ab", fontsize=11)
        y += 18
    for l in lines:
        pg.insert_text((40, y), l, fontname="ar", fontsize=9)
        y += 14
    pg.draw_line((10, 100), (200, 100), width=0.6)
    d.save(path)
    d.close()


def selftest():
    import copy, fitz
    tmp = tempfile.mkdtemp(prefix="af3_selftest_")
    keep = {k: getattr(FG, k) for k in ("OUT", "REG", "LOGD", "ACTIVE_JSON", "MANIFEST", "load_active", "render_log")}
    results = []

    def rep(label, ok, detail=""):
        results.append(bool(ok))
        print("  %-52s %s %s" % (label, "ok" if ok else "MISS", detail))

    try:
        out, reg, logd, leg = (os.path.join(tmp, d) for d in ("out", "registry", "render_log", "legends"))
        for d in (out, reg, logd, leg):
            os.makedirs(d)
        gen = os.path.join(tmp, "P31R_FigS1.R")
        common = os.path.join(FG.HERE, "P31R_fig_common.R")
        io.open(gen, "w").write("# fixture generator\n")
        src_csv = os.path.join(tmp, "src_a.csv")
        io.open(src_csv, "w").write("a,b\n")
        figs = [dict(id="FigureS1", kind="af3", order=1, pages=1, panels=["a", "b"], title="First fixture figure", illustrates=["S-R1"],
                     generator=gen, legend=os.path.join(leg, "P31R_Legend_FigS1.md")),
                dict(id="FigureS2", kind="af3", order=2, pages=2, panels=["a"], title="Second fixture figure with a continuation page", illustrates=["S-R2"],
                     generator=gen, legend=os.path.join(leg, "P31R_Legend_FigS2.md"))]
        active = dict(schema=1, figures=figs, retired=[], text_files=["x"], text_ignore=[], text_dirs=[])
        # G13's one-generator-per-figure rule is not exercised here (units are run explicitly); both fixtures share a file
        ajson = os.path.join(tmp, "active.json")
        json.dump(active, io.open(ajson, "w"))
        manifest = os.path.join(tmp, "manifest.tsv")
        def write_manifest(disp):
            rows = [["AF3_S%d" % k, "-", "fixture series", src_csv.replace("\\", "/"), md5(src_csv), "x", "PIPE", disp, "data", "fixture"] for k in (1, 2)]
            io.open(manifest, "w", encoding="utf-8", newline="").write("\n".join("\t".join(r) for r in [FG.MANIFEST_COLS] + rows) + "\n")
        write_manifest("APPROVED")
        FG.OUT, FG.REG, FG.LOGD, FG.ACTIVE_JSON, FG.MANIFEST = out, reg, logd, ajson, manifest
        FG.load_active = lambda path=None: json.load(io.open(ajson, encoding="utf-8"))
        title = os.path.join(tmp, "TitlePage.md")
        io.open(title, "w", encoding="utf-8").write("**Title (working):** A fixture title for the self-test\n")
        cfg = dict(out=os.path.join(tmp, "af3out", "Additional_file_3.pdf"), log=os.path.join(tmp, "AF3_assembly_log.tsv"), title_page=title)

        def unit(u, fig, tags, lines):
            _fixture_unit(os.path.join(out, u + ".pdf"), tags, lines)
            row = dict(figure=u, w_mm="170", h_mm="120", generator=gen.replace("\\", "/"), generator_md5=md5(gen), common_md5=md5(common), deps="",
                       data="%s=%s" % (src_csv.replace("\\", "/"), md5(src_csv)), pdf_md5=md5(os.path.join(out, u + ".pdf")), png_md5="0" * 32, rendered="x")
            io.open(os.path.join(logd, u + ".tsv"), "w", encoding="utf-8", newline="").write(
                "\t".join(row) + "\n" + "\t".join(row.values()) + "\n")
            json.dump(dict(panels=tags, ledger=[], derived=[], literals=[], sources=[dict(path=src_csv.replace("\\", "/"), panel="")]),
                      io.open(os.path.join(reg, u + ".json"), "w"))
            # a PNG twin for G8
            png = os.path.join(out, u + ".png")
            io.open(png, "wb").write(b"png" + u.encode())
            row["png_md5"] = md5(png)
            io.open(os.path.join(logd, u + ".tsv"), "w", encoding="utf-8", newline="").write("\t".join(row) + "\n" + "\t".join(row.values()) + "\n")
        unit("FigureS1", "FigureS1", ["a", "b"], ["Estimates per signal and outcome", "Interval bars are pointwise ranges"])
        unit("FigureS2", "FigureS2", ["a"], ["Continuation fixture first page text", "Second line of the first page"])
        unit("FigureS2_p2", "FigureS2", [], ["Continuation page carries more detail", "Final line of the continuation"])
        for fe, body in (("S1", "**a**, first panel; **b**, second panel. " + "Words fill the legend. " * 30),
                         ("S2", "**a**, the only panel. " + "More words to wrap over several lines of the legend. " * 90)):
            io.open(os.path.join(leg, "P31R_Legend_Fig%s.md" % fe), "w", encoding="utf-8").write(
                "## Figure %s\n\n**Figure %s. Fixture title for %s.** %s\n" % (fe, fe, fe, body))
        # the active registry names generator by absolute path in this sandbox (v3_path accepts absolute paths); G13 whole-set checks are not run
        # clean build
        try:
            row = build(cfg)
            rep("clean sandbox build assembles and reads back", True, "(%s pages)" % row["pages"])
        except Refuse as e:
            rep("clean sandbox build assembles and reads back", False, str(e)[:300])
            return 1
        v = verify(cfg)
        rep("verify of the fresh PDF is clean", not v, str(v[:2]))
        # the continuation page and legend overflow were exercised
        n = fitz.open(cfg["out"]).page_count
        rep("continuation + legend overflow produce extra pages", n > 1 + 1 + 2, "(%d pages: title + 3 unit pages + overflow)" % n)
        # ---- mutations: each must trip the check it NAMES
        def refuses(label, needle, fn_mut, restore):
            fn_mut()
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    build(dict(cfg, out=os.path.join(tmp, "af3out", "mut.pdf"), log=os.path.join(tmp, "mut_log.tsv")))
                rep(label, False, "built although it should refuse")
            except Refuse as e:
                rep(label, needle in str(e), "(%s)" % str(e).splitlines()[0][:70] if needle in str(e) else "(wrong reason: %s)" % str(e)[:200])
            finally:
                restore()
        p1 = os.path.join(out, "FigureS1.pdf")
        b1 = open(p1, "rb").read()
        refuses("a unit PDF edited after its render", "figure gate is not green", lambda: open(p1, "ab").write(b"\n%edit\n"), lambda: open(p1, "wb").write(b1))
        lg = os.path.join(leg, "P31R_Legend_FigS1.md")
        lb = open(lg, "rb").read()
        refuses("a legend naming the wrong panels", "figure gate is not green", lambda: open(lg, "w", encoding="utf-8").write("## Figure S1\n\n**Figure S1. T.** **a**, only a.\n"), lambda: open(lg, "wb").write(lb))
        refuses("a legend with no bold opener", "figure gate is not green", lambda: open(lg, "w", encoding="utf-8").write("## Figure S1\n\nTitle. **a**, x; **b**, y.\n"), lambda: open(lg, "wb").write(lb))
        rg = os.path.join(reg, "FigureS1.json")
        rb_ = open(rg, "rb").read()
        refuses("a source whose manifest row is SUPERSEDED", "figure gate is not green", lambda: write_manifest("SUPERSEDED"), lambda: write_manifest("APPROVED"))
        refuses("a source whose manifest row is BLOCKED", "figure gate is not green", lambda: write_manifest("BLOCKED"), lambda: write_manifest("APPROVED"))
        # verify catches an input changed after the build
        io.open(title, "a", encoding="utf-8").write("changed\n")
        v = verify(cfg)
        rep("verify: title page changed after the build", any(x.startswith("V4 title_md5") for x in v), str(v[:1]))
        io.open(title, "w", encoding="utf-8").write("**Title (working):** A fixture title for the self-test\n")
        b = open(cfg["out"], "rb").read()
        open(cfg["out"], "ab").write(b"\n%tamper\n")
        v = verify(cfg)
        rep("verify: shipped PDF edited after the build", any(x.startswith("V1") for x in v), str(v[:1]))
        open(cfg["out"], "wb").write(b)
        # a blind readback: the checker must raise/fail on an empty text layer
        try:
            blank = fitz.open()
            blank.new_page()
            bp = os.path.join(tmp, "blank.pdf")
            blank.save(bp)
            plan, up, lt = af3_plan(active), {u: os.path.join(out, u + ".pdf") for u in ("FigureS1", "FigureS2", "FigureS2_p2")}, {"FigureS1": open(lg).read(), "FigureS2": open(os.path.join(leg, "P31R_Legend_FigS2.md")).read()}
            try:
                fl, _ = readback(bp, plan, up, lt)
                rep("read-back of a blank PDF is not clean", bool(fl), str(fl[:1]))
            except Exception as e:
                rep("read-back of a blank PDF raises", True, "(%s)" % type(e).__name__)
        except Exception as e:
            rep("read-back of a blank PDF", False, repr(e))
        # a scaled figure page must be caught (R7): rebuild with the unit drawn at 0.9
        orig_show = fitz.Page.show_pdf_page
        def shrink(self, rect, *a, **k):
            r = fitz.Rect(rect.x0, rect.y0, rect.x0 + rect.width * 0.9, rect.y0 + rect.height * 0.9)
            return orig_show(self, r, *a, **k)
        fitz.Page.show_pdf_page = shrink
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                build(dict(cfg, out=os.path.join(tmp, "af3out", "mut2.pdf"), log=os.path.join(tmp, "mut2_log.tsv")))
            rep("a shrunken figure page is refused", False, "built")
        except Refuse as e:
            rep("a shrunken figure page is refused", "R7" in str(e), "(%s)" % str(e).splitlines()[1][:80] if "R7" in str(e) else str(e)[:160])
        finally:
            fitz.Page.show_pdf_page = orig_show
    finally:
        for k, v in keep.items():
            setattr(FG, k, v)
        shutil.rmtree(tmp, ignore_errors=True)
    print("AF3 ASSEMBLER SELFTEST %s (%d checks)" % ("PASS" if all(results) else "FAIL", len(results)))
    return 0 if all(results) else 1


if __name__ == "__main__":
    a = sys.argv[1:]
    try:
        if "--selftest" in a:
            sys.exit(selftest())
        if "--verify" in a:
            fl = verify()
            for x in fl:
                print("  FAIL " + x)
            print("AF3 VERIFY %s" % ("RED: %d failure(s)" % len(fl) if fl else "GREEN"))
            sys.exit(1 if fl else 0)
        build()
    except Refuse as e:
        print("REFUSED: %s" % e)
        sys.exit(1)
