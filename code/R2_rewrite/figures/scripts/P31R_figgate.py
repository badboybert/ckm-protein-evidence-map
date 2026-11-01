#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 rewrite (R2) -- gate over the SHIPPED figure PDFs (figures/out/Figure_N.pdf).

Reads the PDF (PyMuPDF rawdict, char level), never the PNG, and never the R source alone.
  G0  liveness: every expected figure exists; the reader recovers its panel tags; characters per page sit in a
      plausible band AND agree with a second, independent reader (pdfminer) within 15%  [a blind or flooded reader
      passes every absence check]
  G1  page: 170 or 85 mm wide (+/-0.5), at most 225 mm high
  G1b figure AND legend fit the journal's print height: image height at print width + image-caption gap + caption lines
      <= 225 mm. Cardiovascular Diabetology: "maximum height of 225 mm for figure and legend" (submission guidelines,
      fetched 2026-10-07). Caption geometry is calibrated on published CD PDFs (8-pt captions, 3.56 mm line pitch, ~153
      characters per 170.6-mm line, 3.6 mm gap; R2 review scout G, DOIs 10.1186/s12933-025-02851-0, -02738-0,
      10.1186/s12933-026-03209-w, -03203-2) and must be re-checked on a Word render at Stage 4. An image-only 225 mm check
      (G1) passed every round while figure + legend overran: that is the defect this closes.
  G11 line widths: every VISIBLE line is >= 0.25 pt at print scale ("All lines should be wider than 0.25 pt"): stroked
      lines by their stroke width, differently coloured outlines likewise, and elongated filled shapes (bars, exons; long
      side > 2 pt) by their short side plus any same-colour outline. Point markers never trip it.
  G1b and G11 are WARN-ONLY until the R2 figure stage (Stage 4) sets HEIGHT_BLOCKING (or env P31R_HEIGHT_GATE=block):
  blocking now would freeze every build (P31R_build.py refuses on a red figure gate) before the figures are redesigned.
  G2  type sizes: every span is 7 pt (annotation), 9 pt (body) or 11 pt (panel tag); only short numeric superscripts
      (exponents) may be smaller, never below 4.5 pt
  G3  glyph coverage: every drawn char exists in the INSTALLED font file its span names (cairo boxes a missing glyph in
      the PDF while the PNG twin falls back silently); unknown fonts fail
  G4  nothing drawn off the page
  G5  no two text lines overlap (bbox intersection > 15% of the smaller line)
  G6  vocabulary sheet: no banned term (shared list from P31R_gate.py)
  G7  numbers vs source: every printed number OCCURRENCE is claimed by exactly one entry of figures/registry/<fig>.json
      ("count", default 1). A ledger entry names ONE element of its row ("idx", 1-based like p31r_ledger; required when
      the row holds more than one number) and is compared at its printed precision, sign included; a derived entry is a
      ratio of two ledger elements; a literal needs why = "<axis_tick|coordinate|scale|other>: <reason>". Axis ticks name
      their "axis" (and "transform": linear/log10/sqrt/asinh/custom); a linear or log10 axis must be evenly stepped, with a
      round linear step. A non-tick literal may not equal a ledger value (any value of a cited row, any integer of any
      row) unless it says "not_ledger": "<why>". An entry whose text another entry shares needs an "anchor" (a regex with
      {} for the number, matched against printed lines), so two roles of one number cannot be swapped. "%" is part of
      the token. Digits inside gene/variant names, cytobands,
      figure pointers, chr coordinates, 'Type 2' etc. are masked with the prose gate's EXEMPT list. G7m: a printed
      negative must use U+2212, not a hyphen-minus. Exponents are read only from a smaller, RAISED span after "10".
  G8  provenance: the render-log row's generator / common / deps / data md5s equal the files on disk now, and the
      PDF and PNG md5s equal the shipped files [an edit after the render passes every other gate]
  G9  legend coupling: sections/P31R_Legend_FigN.md exists and names exactly the panels the PDF carries
  G10 text coupling: every "Figure Nx" cited in the manuscript sections names a panel that exists; uncited panels
      are listed as a note
AF3 PIPELINE (added 2026-10-16, decision S2-D24; G0-G11 above are unchanged). The scope is the ACTIVE registry
figures/af3/P31R_ACTIVE_FIGURES.json (main Figure_1-5, AF3 FigureS<n>, continuation pages FigureS<n>_p2 ...), not a glob:
  G9  (AF3)  legend = af3/legends/P31R_Legend_FigS<n>.md: '## Figure S<n>' heading, '**Figure S<n>. Title.**' opener, bold panel letters
      equal to the union of the panels of the figure's pages
  G11 (AF3)  blocking (a main figure's G11 stays warn-only until Stage 4); G1b and G10 are NOT applied to AF3 (decision F9; G14 replaces G10)
  G12 sources    every plotted series names ONE source in the approved manifest (af3/P31R_SOURCE_MANIFEST.tsv, env P31R_SOURCE_MANIFEST):
                 the file the render read is listed at the md5 it was read at, status approved, for this figure and series; a superseded
                 version, an unlisted file, a file under a hold (status held, e.g. H9) or an approved file edited since approval fails;
                 acceptance is by manifest row (path + md5), never by directory name
  G13 scope      active registry schema; every active generator exists; a scripts/P31R_Fig*.R that is neither active nor under 'retired'
                 fails; planned panels equal the panels the pages carry; a shipped PDF older than its generator, shim, data or ledger fails
  G14 pointers   over the enumerated text/SM files (text_files; any other .md in text_dirs must match text_ignore): every active figure is
                 cited, every explicit panel citation resolves, a whole-figure citation covers its panels (panels cited by no letter, in a
                 figure cited only by letters, are a NOTE)
KNOWN LIMITS of the AF3 checks:
  - G12 proves the files read are approved versions and each series names one source; it does not prove a plotted value came from that file
    (the generator's source frame af3/source_data/<unit>__<name>.csv and the independent visual check cover that);
  - G13 staleness compares modification times: a copy that preserves mtime (shutil.copy2) can hide a stale file; the md5 pins (G8) cover content;
  - G14 reads 'Figure N' / 'Fig. N' as OUR figure: a citation of another paper's figure in the enumerated text would be misread (none today);
  - G10 still reads the superseded R1-era ../sections/ text for the five main figures and is left as it was.
KNOWN LIMITS (3a-iii integration audit; the independent figure verifiers cover them, the gate does not):
  - a role swap between two numbers that each have a unique printed text and no anchor passes (G7 checks values and
    counts, not position); anchors are required only where texts collide;
  - axis ticks are pooled per (text, anchor) across the figure, not bound to their axis geometry, so a missing tick
    offset by a lone label with the same text, or two ticks swapped in place, passes;
  - a line-anchored pattern (^...$) counts one match per line, so two identical tokens on one merged line undercount
    (fails closed);
  - G5 compares text with text only: a line or point crossing a label is not seen.
Run:  python P31R_figgate.py [--selftest] [--pointers | --active] [Figure_3 FigureS1 FigureS5_p2 ...]    exit 1 on any failure.
      no unit given = every active unit + the global checks (G13 scope, G12 manifest integrity, G14 pointers); --pointers = G14 only with
      the citation table; --active = the global checks only.
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
import csv, hashlib, io, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.dirname(HERE)
ROOT = os.path.dirname(FIGS)
V3 = os.path.dirname(ROOT)
OUT = os.path.join(FIGS, "out")
REG = os.path.join(FIGS, "registry")
SEC = os.path.join(ROOT, "sections")
LOGD = os.path.join(FIGS, "render_log")   # one <figure>.tsv per figure (written by p31r_save)
EXPECTED = ["Figure_%d" % i for i in range(1, 6)]
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import P31R_gate as PG  # noqa: E402  (shared ledger parsing, value matching and the banned-vocabulary list)

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

SIZES = (7.0, 9.0, 11.0)
FONTDIR = _P31R_FONTDIR
FONTFILE = {("arial", False, False): "arial.ttf", ("arial", True, False): "arialbd.ttf",
            ("arial", False, True): "ariali.ttf", ("arial", True, True): "arialbi.ttf"}
# a bare power (10^-6) is tried FIRST, or the "10" is consumed as a plain number and the exponent lost (3a gate adversary);
# a thousands group is exactly 3 digits, so "1,948," and "3,45,33" are not read as one number
NUMRX = re.compile(r"(?<![\w.^_])10\^([−\-+]?\d+)(?!\d)|(?<![\w.^_])([+−\-]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?%?)(?!\d)(?!,\d)(\s?×\s?10\^([−\-+]?\d+)(?!\d))?")
WHY_KINDS = ("axis_tick", "coordinate", "scale", "other")
SUPMAP = {ord(k): v for k, v in zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")}
# G1b / G11 (R2 review section 21): print geometry of Cardiovascular Diabetology, see the docstring for the calibration
PAGE_MAX_MM = 225.0
PRINT_W_MM = 170.0            # stricter than the ~166 mm observed in print: a full-width image is not enlarged
CAPTION_GAP_MM = 3.6
CAPTION_LINE_MM = 3.56
CAPTION_CHARS_PER_LINE = 153
MIN_STROKE_PT = 0.25
HEIGHT_BLOCKING = os.environ.get("P31R_HEIGHT_GATE", "").lower() == "block"   # Stage 4 flips this to True


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


# ---------------------------------------------------------------- extraction
def extract(pdf):
    import fitz
    doc = fitz.open(pdf)
    if doc.page_count != 1:
        raise RuntimeError("%s: %d pages (want 1)" % (pdf, doc.page_count))
    pg = doc[0]
    rect = pg.rect
    lines, spans = [], []
    for b in pg.get_text("rawdict")["blocks"]:
        for ln in b.get("lines", []):
            ls = []
            for s in ln["spans"]:
                t = "".join(c["c"] for c in s["chars"])
                if not t:
                    continue
                sp = dict(text=t, size=s["size"], font=s["font"], flags=s["flags"], bbox=s["bbox"],
                          chars=[(c["c"], c["bbox"]) for c in s["chars"]])
                if t.strip():
                    spans.append(sp)
                ls.append(sp)
            if any(x["text"].strip() for x in ls):
                horiz = abs(ln.get("dir", (1, 0))[1]) < 0.01
                lines.append(dict(text=assemble(ls, horiz), bbox=ln["bbox"], spans=[x for x in ls if x["text"].strip()], horiz=horiz))
    strokes = [dict(type=d.get("type"), width=d.get("width"), color=d.get("color"), fill=d.get("fill"),
                    rw=d["rect"].width, rh=d["rect"].height) for d in pg.get_drawings()]
    out = dict(w_mm=rect.width / 72 * 25.4, h_mm=rect.height / 72 * 25.4, rect=tuple(rect), lines=lines, spans=spans, strokes=strokes)
    doc.close()
    return out


def _script(sp, prev):
    """'^' for a superscript, '_' for a subscript: a short numeric span smaller than its neighbour (< 0.85x) and shifted
    by at least 0.15 em from the neighbour's baseline (raised = superscript, lowered = subscript)."""
    if prev is None or not re.fullmatch(r"\s*[−\-+]?\d{1,3}\s*", sp["text"]) or sp["size"] >= 0.85 * prev["size"]:
        return ""
    shift = prev["bbox"][3] - sp["bbox"][3]
    if shift >= 0.15 * prev["size"]:
        return "^"
    if -shift >= 0.15 * prev["size"] or sp["bbox"][1] > prev["bbox"][1] + 0.25 * prev["size"]:
        return "_"
    return ""


def assemble(ls, horiz=True):
    """Join a line's spans. A superscript is marked '^' (an exponent when it follows '10'; otherwise not a data number),
    a subscript '_' (log10, p12: never a data number). Between two ordinary spans that differ in size or font a space
    is inserted, so a 7-pt count printed right after 9-pt text stays its own token (cairo leaves only ~0.4 pt)."""
    parts, prev = [], None
    for sp in ls:
        t = sp["text"]
        if not t.strip():
            parts.append(" ")
            continue
        mark = _script(sp, prev) if horiz else ""   # bbox y says nothing about scripts in rotated text
        if not mark and prev is not None and parts and not parts[-1].endswith(" ") and not t.startswith(" "):
            if abs(sp["size"] - prev["size"]) > 0.3 or sp["font"] != prev["font"] or sp["bbox"][0] - prev["bbox"][2] > 0.2 * sp["size"]:
                parts.append(" ")
        parts.append(mark + t.strip() if mark else t)
        if not mark:
            prev = sp
    return "".join(parts)


def masked(text):
    """Blank the figure exemptions shared with the prose gate (P31R_gate.EXEMPT_FIG: gene and variant names, cytobands,
    figure pointers, 'Type 2', 'PP.H4', chr coordinates, IL-6 ...), so their digits are never read as numbers."""
    for pat in PG.EXEMPT_FIG:
        text = re.sub(pat, lambda m: "§" * len(m.group(0)), text)
    return text


def second_reader_chars(pdf):
    from pdfminer.high_level import extract_text
    return len(re.sub(r"\s", "", extract_text(pdf)))


def numbers(text):
    return [m.group(0).strip() for m in NUMRX.finditer(masked(text))]


def unread_digits(text):
    """Digits the tokenizer did not consume and no exemption masked: '14,5210', '3,45,33', '1e-06' leave such residues."""
    t = masked(text)
    t = NUMRX.sub(lambda m: "#" * len(m.group(0)), t)
    t = re.sub(r"[\^_][−\-+]?\d+", "", t)     # superscripts/subscripts outside a power (r^2, log_10)
    return re.findall(r"\S*\d\S*", t)


def norm(tok):
    """One canonical spelling for a printed token or a registry text: U+2212 -> '-', Unicode superscripts -> '^digits',
    '× 10' always followed by '^', no spaces, no trailing '%'."""
    t = tok.strip()
    t = re.sub(r"([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)", lambda m: "^" + m.group(1).translate(SUPMAP), t)
    t = t.replace("−", "-").replace(" ", "")
    t = re.sub(r"×10(?=[-+]?\d)", "×10^", t)
    return t


def to_gate_form(tok):
    """Canonical token -> the spelling P31R_gate.parse_text_num understands: 'a × 10e' for a mantissa with a power, and
    '1 × 10e' for a bare power (positive or negative)."""
    t = norm(tok).lstrip("+").rstrip("%")
    m = re.fullmatch(r"(-?\d[\d,]*(?:\.\d+)?)×10\^([-+]?\d+)", t)
    if m:
        return "%s × 10%s" % (m.group(1), m.group(2).lstrip("+"))
    m = re.fullmatch(r"10\^([-+]?\d+)", t)
    if m:
        return "1 × 10%s" % m.group(1).lstrip("+")
    return t


def panel_tags(ex):
    return sorted({s["text"].strip() for s in ex["spans"]
                   if abs(s["size"] - 11.0) < 0.15 and (s["flags"] & 16) and re.fullmatch(r"[a-h]", s["text"].strip())})


# ---------------------------------------------------------------- checks (each takes extracted data; pure, so the self-test can feed defects)
def g1_page(ex):
    f = []
    if not any(abs(ex["w_mm"] - w) <= 0.5 for w in (170.0, 85.0)):
        f.append("G1 width %.1f mm (want 170 or 85)" % ex["w_mm"])
    if ex["h_mm"] > 225.5:
        f.append("G1 height %.1f mm > 225" % ex["h_mm"])
    return f


def legend_caption_text(md):
    """The legend as it prints: ledger tags, markdown emphasis and the '## Figure N' heading line removed."""
    t = re.sub(r"\[\[[^\]]*\]\]", "", md)
    t = "\n".join(l for l in t.splitlines() if not l.startswith("#"))
    t = t.replace("**", "")
    return re.sub(r"\s+", " ", t).strip()


def g1b_height(ex, caption):
    """Image at print width + gap + caption lines must fit the 225 mm print height (figure AND legend)."""
    if not caption:
        return ["G1b no legend text to measure"]
    lines = -(-len(caption) // CAPTION_CHARS_PER_LINE)
    img = ex["h_mm"] * PRINT_W_MM / ex["w_mm"] if ex["w_mm"] > 100 else ex["h_mm"]   # an 85-mm figure is not enlarged
    total = img + CAPTION_GAP_MM + lines * CAPTION_LINE_MM
    if total > PAGE_MAX_MM + 0.05:
        return ["G1b figure + legend %.1f mm > %.0f (image %.1f + gap %.1f + %d caption lines %.1f; %d legend characters): "
                "shorten the legend or the image by %.1f mm" % (total, PAGE_MAX_MM, img, CAPTION_GAP_MM, lines,
                                                               lines * CAPTION_LINE_MM, len(caption), total - PAGE_MAX_MM)]
    return []


def _same_colour(a, b):
    return a is not None and b is not None and len(a) == len(b) and all(abs(x - y) < 1e-3 for x, y in zip(a, b))


def g11_strokes(ex):
    """VISIBLE line width at print scale must be >= 0.25 pt (R2 Stage 0 verifier: measure what prints, not the stroke
    operand). A stroked line prints its stroke width. A filled shape whose outline has the fill colour prints its short side
    plus the stroke, so it is checked only when elongated (long side > 2 pt): a thin bar, exon or rule. A filled shape with a
    differently coloured outline prints that outline as a line. A filled shape without outline prints its short side.
    Point markers (both sides >= 0.25 pt) never trip the check."""
    scale = PRINT_W_MM / ex["w_mm"] if ex.get("w_mm", 0) > 100 else 1.0
    lines, bars = [], []
    for d in ex.get("strokes", []):
        t, w = d.get("type"), d.get("width") or 0.0
        short, long_ = min(d.get("rw", 0), d.get("rh", 0)), max(d.get("rw", 0), d.get("rh", 0))
        if t == "s":
            if w * scale < MIN_STROKE_PT:
                lines.append(w)
        elif t == "fs" and not _same_colour(d.get("color"), d.get("fill")):
            if w * scale < MIN_STROKE_PT:
                lines.append(w)
        elif t in ("fs", "f"):
            vis = short + (w if t == "fs" else 0.0)
            if long_ > 2.0 and vis * scale < MIN_STROKE_PT:
                bars.append(vis)
    out = []
    if lines:
        out.append("G11 %d stroked line(s) thinner than %.2f pt at print scale (thinnest %.3f pt)" % (len(lines), MIN_STROKE_PT, min(lines)))
    if bars:
        out.append("G11 %d filled bar(s) or exon(s) narrower than %.2f pt at print scale (narrowest %.3f pt)" % (len(bars), MIN_STROKE_PT, min(bars)))
    return out


def route_geometry(geo, blocking):
    """G1b/G11 results: failures when blocking, warnings otherwise. (fails, notes)"""
    return (list(geo), []) if blocking else ([], ["  WARN (not yet blocking) " + x for x in geo])


def g2_sizes(ex):
    """Every span is 7, 9 or 11 pt, except a superscript or subscript (smaller than and shifted from its neighbour,
    numeric, at most 3 digits), which may be smaller but not below 4.5 pt."""
    f = []
    groups = [(l["spans"], l.get("horiz", True)) for l in ex["lines"] if l.get("spans")] or [([sp], True) for sp in ex["spans"]]
    for grp, horiz in groups:
        prev = None
        for sp in grp:
            sz = round(sp["size"], 1)
            ok = any(abs(sz - a) <= 0.15 for a in SIZES)
            if not ok and sz >= 4.5 and horiz and _script(sp, prev):
                ok = True
            if not ok:
                f.append("G2 %.1f pt span %r" % (sz, sp["text"][:40]))
            if not _script(sp, prev):
                prev = sp
    return f


_fontcache = {}


def _font_for(name):
    n = name.lower()
    fam = "arial" if "arial" in n else None
    if fam is None:
        return None
    bold = "bold" in n
    ital = "italic" in n or "oblique" in n
    fn = FONTFILE.get((fam, bold, ital))
    if fn not in _fontcache:
        import fitz
        _fontcache[fn] = fitz.Font(fontfile=os.path.join(FONTDIR, fn))
    return _fontcache[fn]


def g3_glyphs(ex):
    f = []
    for s in ex["spans"]:
        if "symbol" in s["font"].lower():
            continue
        ft = _font_for(s["font"])
        if ft is None:
            f.append("G3 unknown font %r for %r" % (s["font"], s["text"][:30]))
            continue
        for c, _ in s["chars"]:
            if c.strip() and not ft.has_glyph(ord(c)):
                f.append("G3 glyph U+%04X %r missing from %s (span %r)" % (ord(c), c, s["font"], s["text"][:30]))
    return f


def g4_offpage(ex):
    f = []
    x0, y0, x1, y1 = ex["rect"]
    for s in ex["spans"]:
        for c, bb in s["chars"]:
            if not c.strip():
                continue
            if bb[0] < x0 - 0.5 or bb[1] < y0 - 0.5 or bb[2] > x1 + 0.5 or bb[3] > y1 + 0.5:
                f.append("G4 char %r of %r off the page" % (c, s["text"][:30]))
                break
    return f


def _inter(a, b):
    w = min(a[2], b[2]) - max(a[0], b[0])
    h = min(a[3], b[3]) - max(a[1], b[1])
    return w * h if w > 0 and h > 0 else 0.0


def _area(a):
    return max(0.0, a[2] - a[0]) * max(0.0, a[3] - a[1])


def g5_overlap(ex):
    f = []
    L = ex["lines"]
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            a, b = L[i]["bbox"], L[j]["bbox"]
            ov = _inter(a, b)
            if ov and ov > 0.15 * min(_area(a), _area(b)):
                f.append("G5 overlapping text %r / %r" % (L[i]["text"][:30], L[j]["text"][:30]))
    return f


def g6_vocab(ex):
    f = []
    txt = " ".join(l["text"] for l in ex["lines"])
    for pat in PG.BANNED:
        for m in re.finditer(pat, txt, flags=re.I):
            f.append("G6 banned term %r ... %s ..." % (m.group(0), txt[max(0, m.start() - 30):m.end() + 30]))
    for m in re.finditer(r"replicat\w*", txt, flags=re.I):
        f.append("G6 banned term %r" % m.group(0))
    return f


def _flat_get(v, idx):
    """Element of a ledger value by index path: an integer is a 1-BASED list position (as p31r_ledger(id, idx) in R), a
    string is a dict key; a list of these walks a nested value. Integers on a dict are refused (use the key)."""
    for k in (idx if isinstance(idx, list) else [idx]):
        if isinstance(k, bool) or not isinstance(k, (int, str)):
            raise TypeError("idx step %r: use an integer (list) or a string (dict key)" % (k,))
        if isinstance(k, int):
            if not isinstance(v, list):
                raise TypeError("integer idx %d on a %s value: use the dict key" % (k, type(v).__name__))
            if not 1 <= k <= len(v):
                raise IndexError("idx %d outside 1..%d" % (k, len(v)))
            v = v[k - 1]
        else:
            if not isinstance(v, dict):
                raise TypeError("key %r on a %s value" % (k, type(v).__name__))
            v = v[k]
    return v


def _elem(L, lid, idx=None):
    v = json.loads(L[lid]["value"].strip())
    x = _flat_get(v, idx) if idx is not None else v
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        raise TypeError("%s%s is not a number: %r" % (lid, "" if idx is None else "[%r]" % (idx,), x))
    return float(x)


def _numeric_leaves(v):
    if isinstance(v, bool):
        return []
    if isinstance(v, (int, float)):
        return [float(v)]
    if isinstance(v, list):
        return [x for e in v for x in _numeric_leaves(e)]
    if isinstance(v, dict):
        return [x for e in v.values() for x in _numeric_leaves(e)]
    return []


_LEAVES = {}


def _all_leaves(L):
    if not _LEAVES:
        for lid, r in L.items():
            try:
                _LEAVES[lid] = _numeric_leaves(json.loads(r["value"]))
            except ValueError:
                _LEAVES[lid] = []
    return _LEAVES


def _round(x):
    """A tick step is a 'nice' number: 1, 2, 2.5, 3, 4 or 5 times a power of ten (the steps axis-break algorithms use)."""
    import math
    if x <= 0:
        return False
    m = x / 10 ** math.floor(math.log10(x) + 1e-12)
    return any(abs(m - k) < 1e-6 for k in (1, 2, 2.5, 3, 4, 5, 10))


TRANSFORMS = ("linear", "log10", "sqrt", "asinh", "custom")


def _ticks_ok(ax, tl, transform):
    """At least two ticks; on a linear or log10 axis the (transformed) steps are constant and the linear step is round."""
    import math
    out = []
    if len(tl) < 2:
        return ["G7 axis %r declares only %d tick: an axis carries at least two" % (ax, len(tl))]
    try:
        v = sorted({PG.parse_text_num(to_gate_form(t))[0] for t in tl})
    except ValueError:
        return ["G7 axis %r has an unparsable tick" % ax]
    if transform not in TRANSFORMS:
        return ["G7 axis %r: transform %r not in %s" % (ax, transform, "/".join(TRANSFORMS))]
    if transform in ("linear", "log10") and len(v) >= 3:
        w = [math.log10(x) for x in v] if transform == "log10" and all(x > 0 for x in v) else v
        st = [b - a for a, b in zip(w, w[1:])]
        if max(st) - min(st) > 1e-6 * max(1.0, abs(max(st))):
            out.append("G7 axis %r: ticks %s are not evenly spaced on a %s scale" % (ax, v, transform))
        elif transform == "linear" and not _round(st[0]):
            out.append("G7 axis %r: tick step %g is not a round number" % (ax, st[0]))
    return out


def _schema(kind, e):
    errs = []
    if not isinstance(e, dict) or not isinstance(e.get("text"), str) or not e.get("text").strip():
        return ["entry without a text string: %r" % (e,)]
    c = e.get("count", 1)
    if isinstance(c, bool) or not isinstance(c, int) or c < 1:
        errs.append("count must be an integer >= 1")
    if kind == "ledger" and not isinstance(e.get("id"), str):
        errs.append("ledger entry needs an id")
    if kind == "derived" and not (isinstance(e.get("num"), list) and isinstance(e.get("den"), list)):
        errs.append("derived entry needs num and den")
    if kind == "literal":
        why = str(e.get("why", ""))
        k, _, reason = why.partition(":")
        if k.strip() not in WHY_KINDS or not reason.strip():
            errs.append("literal needs why = '<kind>: <reason>' with kind in %s and a reason" % "/".join(WHY_KINDS))
        if k.strip() == "axis_tick" and not isinstance(e.get("axis"), str):
            errs.append("axis_tick literal needs an 'axis' name")
        if k.strip() == "axis_tick" and e.get("transform", "linear") not in TRANSFORMS:
            errs.append("axis_tick transform must be one of %s" % "/".join(TRANSFORMS))
    return ["G7 registry %r: %s" % (e.get("text") if isinstance(e, dict) else e, x) for x in errs]


FORMAT_RULES = [(r"\d[eE][+\-−]?\d", "e-notation (write × 10 with a superscript)"),
                (r"\d\s?[xX]\s?10\b", "a letter x for × in a power of ten"),
                (r"\bNA\b|\bNaN\b|\bnan\b|\bInf\b", "a missing or infinite value printed"),
                (r"(?<![\d\w.])\.\d", "a decimal without its leading zero"),
                (r"(?:^|(?<=[\s(\[=<>≤≥,:;/^±]))–(?=\d)", "an en-dash used as a minus sign")]


def g7_numbers(ex, reg, L):
    """Every printed number OCCURRENCE is claimed by exactly one registry entry ("count", default 1). A ledger entry is
    compared with ONE element of its row: "idx" (1-based list position or dict key), or, without idx, the row's single
    numeric value. A literal needs why = "<axis_tick|coordinate|scale|other>: <reason>"; axis ticks need an "axis" name,
    must be round and come at least three per axis; any other literal that equals a value of ANY ledger row fails unless
    it carries "not_ledger": "<why the coincidence is not a claim>". Digits the tokenizer cannot read fail."""
    from collections import Counter, defaultdict
    f = []
    printed = Counter()
    for l in ex["lines"]:
        for t in numbers(l["text"]):
            printed[norm(t)] += 1
        mk = masked(l["text"])
        for m in re.finditer(r"(?:^|(?<=[A-Za-z])(?=-\d+\.\d)|(?<=[\s(\[=<>≤≥,:;/^±]))-(?=\d)", mk):
            f.append("G7m hyphen-minus before a number (use U+2212): %r" % l["text"][:50])
        for rx, what in FORMAT_RULES:
            if re.search(rx, mk):
                f.append("G7f %s: %r" % (what, l["text"][:50]))
        for u in unread_digits(l["text"]):
            f.append("G7 digits the gate cannot read as one number: %r in %r" % (u, l["text"][:50]))
    declared = Counter()
    ticks = defaultdict(list)
    tfm = defaultdict(set)
    anchor_groups = defaultdict(list)
    texts = Counter(norm(e["text"]) for kind in ("ledger", "derived", "literals") for e in reg.get(kind, [])
                    if isinstance(e, dict) and isinstance(e.get("text"), str))
    entries = [("ledger", e) for e in reg.get("ledger", [])] + [("derived", e) for e in reg.get("derived", [])] + \
              [("literal", e) for e in reg.get("literals", [])]
    for kind, e in entries:
        se = _schema(kind, e)
        if se:
            f += se
            continue
        declared[norm(e["text"])] += int(e.get("count", 1))
        # role anchors: an entry whose text another entry shares must say WHERE it is printed ("anchor": a regex with {}
        # standing for the number, matched against whole printed lines); any entry may carry one, and it is then checked
        if texts[norm(e["text"])] > 1 and not e.get("anchor"):
            f.append("G7 %r is declared by %d entries: each needs an 'anchor' naming its line" % (e["text"], texts[norm(e["text"])]))
        if e.get("anchor"):
            # anchors are checked per (text, anchor) GROUP after the loop: the summed counts of every entry sharing a text
            # and an anchor must equal the printed occurrences exactly (ticks such as "0" on several axes share both)
            anchor_groups[(norm(e["text"]), e["anchor"])].append(e)
        if kind == "literal":
            k = str(e["why"]).partition(":")[0].strip()
            if k == "axis_tick":
                ticks[e["axis"]].append(e["text"])
                tfm[e["axis"]].add(e.get("transform", "linear"))
                continue
            try:
                tv, tol, sci = PG.parse_text_num(to_gate_form(e["text"]))
            except ValueError:
                f.append("G7 unparsable literal %r" % e["text"])
                continue
            if e.get("not_ledger"):
                continue
            # a literal may not stand in for a ledger value: compared with every value of the rows this registry cites,
            # and with the INTEGER values of every row (continuous leaves of uncited rows match almost any short decimal)
            cited = {x.get("id") for x in reg.get("ledger", [])}
            hit = [lid for lid, lv in _all_leaves(L).items()
                   if PG.matches(tv, tol, sci, lv if lid in cited else [x for x in lv if float(x).is_integer()], False) is True]
            if hit and not (k == "scale" and float(tv).is_integer() and abs(tv) <= 12):
                f.append("G7 literal %r equals a value of %s: declare it as a ledger entry, or add \"not_ledger\": \"<why>\"" % (e["text"], ", ".join(hit[:4])))
        elif kind == "derived":
            try:
                nv, dv = (_elem(L, *e["num"]), _elem(L, *e["den"]))
                val = nv / dv * (100 if e.get("pct") else 1)
                tv, tol, _ = PG.parse_text_num(to_gate_form(e["text"].rstrip("%")))
            except Exception as ex_:
                f.append("G7 derived entry %r unusable: %s" % (e["text"], ex_))
                continue
            if not abs(tv - val) < tol * (1 - 1e-9) + 1e-15:
                f.append("G7 derived %r != %s/%s = %.6g" % (e["text"], e["num"], e["den"], val))
        else:
            if e["id"] not in L:
                f.append("G7 unknown ledger id %s for %r" % (e["id"], e["text"]))
                continue
            try:
                tv, tol, sci = PG.parse_text_num(to_gate_form(e["text"]))
            except ValueError:
                f.append("G7 unparsable registry number %r" % e["text"])
                continue
            try:
                if "idx" in e:
                    cands = [_elem(L, e["id"], e["idx"])]
                else:
                    leaves = sorted(set(_all_leaves(L)[e["id"]]))
                    if len(leaves) != 1:
                        f.append("G7 %r: ledger row %s holds %d numeric values; name the element with idx" % (e["text"], e["id"], len(leaves)))
                        continue
                    cands = leaves
            except Exception as ex_:
                f.append("G7 %r: %s idx %r unusable (%s)" % (e["text"], e["id"], e.get("idx"), ex_))
                continue
            mm = PG.matches(tv, tol, sci, cands, bool(e.get("pct")))
            if mm == "BOUNDARY":
                f.append("G7 %r sits on a rounding boundary of %s (print from full precision)" % (e["text"], e["id"]))
            elif not mm:
                f.append("G7 %r does not match ledger %s%s (%s)" % (e["text"], e["id"], "[%s]" % (e["idx"],) if "idx" in e else "",
                                                                   L[e["id"]]["display"][:60]))
    for (nt, anc), es in anchor_groups.items():
        txt = es[0]["text"]
        num = re.escape(txt).replace("−", "[−-]").replace(r"\-", "[−-]").replace(r"\ ", r"\s?")
        try:
            arx = re.compile(anc.replace("{}", "(?:%s)" % num))
        except re.error as ex_:
            f.append("G7 anchor %r unusable: %s" % (anc, ex_))
            continue
        hits = sum(len(arx.findall(l["text"].replace("^", ""))) for l in ex["lines"])
        want = sum(int(e.get("count", 1)) for e in es)
        if hits != want:
            f.append("G7 %r: anchor %r matches %d printed occurrence(s), the entries using it declare %d" % (txt, anc, hits, want))
    for ax, tl in ticks.items():
        if len(tfm[ax]) > 1:
            f.append("G7 axis %r declares conflicting transforms %s" % (ax, sorted(tfm[ax])))
        f += _ticks_ok(ax, tl, next(iter(tfm[ax])))
    for t in sorted(set(printed) | set(declared)):
        if printed[t] != declared[t]:
            if declared[t] == 0:
                f.append("G7 printed number %r (x%d) not declared in the registry" % (t, printed[t]))
            elif printed[t] == 0:
                f.append("G7 registry entry %r is not printed on the figure (stale registry?)" % t)
            else:
                f.append("G7 %r printed %d time(s) but declared %d" % (t, printed[t], declared[t]))
    return f


def g8_provenance(fig, row):
    f = []
    if row is None:
        return ["G8 no render-log row for %s" % fig]
    def chk(rel, want, what):
        p = rel if os.path.isabs(rel) or re.match(r"^[A-Za-z]:", rel) else os.path.join(V3, rel)
        if not os.path.exists(p):
            f.append("G8 %s missing: %s" % (what, p))
        elif md5(p) != want:
            f.append("G8 %s changed since the render: %s" % (what, rel))
    chk(row["generator"], row["generator_md5"], "generator")
    chk(os.path.join(HERE, "P31R_fig_common.R"), row["common_md5"], "common helper")
    for pin in filter(None, row["deps"].split(";")):
        p, _, h = pin.rpartition("=")
        chk(p, h, "house dependency")
    for pin in filter(None, row["data"].split(";")):
        p, _, h = pin.rpartition("=")
        chk(p, h, "data file")
    if not row["data"]:
        f.append("G8 %s read no registered data file (a figure with no data pins cannot be re-derived)" % fig)
    chk(os.path.join(OUT, fig + ".pdf"), row["pdf_md5"], "shipped PDF")
    chk(os.path.join(OUT, fig + ".png"), row["png_md5"], "PNG raster")
    return f


def legend_panels(txt):
    body = re.sub(r"\[\[[^\]]*\]\]", "", txt)
    got = set(re.findall(r"\*\*\(?([a-h])\)?\*\*", body))
    got |= set(re.findall(r"(?:^|[\s(])\(([a-h])\)", body))
    return sorted(got)


def main_legend_path(fig, active=None):
    try:
        a = active if active is not None else load_active()
    except Exception:
        a = None
    for e in (a or {}).get("figures", []):
        if e.get("id") == fig and e.get("kind") == "main" and e.get("legend"):
            return v3_path(e["legend"])
    return os.path.join(SEC, "P31R_Legend_Fig%s.md" % fig.split("_")[1])


def g9_legend(fig, ex, active=None):
    n = fig.split("_")[1]
    p = main_legend_path(fig, active)
    if not os.path.exists(p):
        return ["G9 legend file missing: %s" % os.path.basename(p)]
    tags = panel_tags(ex)
    lp = legend_panels(io.open(p, encoding="utf-8").read())
    if len(tags) <= 1 and not lp:
        return []
    if tags != lp:
        return ["G9 legend panels %s != figure panels %s" % (lp, tags)]
    return []


def cited_panels():
    cites = {}
    for f in os.listdir(SEC):
        if not (f.startswith("P31R_") and f.endswith(".md")) or f.startswith("P31R_Legend_"):
            continue
        t = io.open(os.path.join(SEC, f), encoding="utf-8").read()
        for m in re.finditer(r"Figures?\s+(\d)([a-h](?:[–-][a-h])?(?:,\s?[a-h])*)?", t):
            n = int(m.group(1))
            letters = set()
            if m.group(2):
                for part in m.group(2).split(","):
                    part = part.strip()
                    if re.fullmatch(r"[a-h][–-][a-h]", part):
                        letters |= {chr(c) for c in range(ord(part[0]), ord(part[-1]) + 1)}
                    else:
                        letters.add(part)
            cites.setdefault(n, set()).update(letters or {""})
    return cites


def g10_text(fig, ex, cites):
    f, notes = [], []
    n = int(fig.split("_")[1])
    tags = set(panel_tags(ex))
    got = cites.get(n, set())
    if not got:
        f.append("G10 %s is never cited in the manuscript sections" % fig)
    for l in sorted(x for x in got if x):
        if l not in tags:
            f.append("G10 text cites %s%s but the PDF has panels %s" % (fig.replace("_", " "), l, sorted(tags)))
    unc = sorted(tags - got) if "" not in got else []
    if unc and len(tags) > 1:
        notes.append("  G10 note: %s panels never cited by letter: %s" % (fig, unc))
    return f, notes


# ---------------------------------------------------------------- active registry, source manifest, pointers (AF3 pipeline; S2-D24)
# Added 2026-10-16 (plan FIGURE_MESSAGES_R2 section 6; decision S2-D24; reviewer brief "Reuse the verification system; do not expand it").
# NOTHING above this line changed meaning: G0-G11 are the same functions with the same rules. These are NEW checks:
#   G12 source manifest   every plotted series names ONE source; the file the render read is an approved version in the
#                         manifest (by path AND md5, not by directory); a superseded or unlisted file fails
#   G13 active registry   the figure list in af3/P31R_ACTIVE_FIGURES.json is the scope: every active generator exists, a generator on
#                         disk that is neither active nor retired fails (a retired generator is never picked up by a glob),
#                         planned panels equal rendered panels, an output is not older than its generator or data
#   G14 pointers          every active figure is cited in the enumerated text/SM files, every explicit panel citation resolves,
#                         a whole-figure citation covers its panels (panels never cited by letter are a NOTE, not a failure)
# AF3 units (FigureS<n>, continuation pages FigureS<n>_p2 ...) get G0-G9 and G11 (blocking) and the new checks; NOT G1b (a main-text
# print rule, decision F9) and not G10 (superseded by G14, which reads the R2 text).
AF3 = os.path.join(FIGS, "af3")
ACTIVE_JSON = os.path.join(AF3, "P31R_ACTIVE_FIGURES.json")
MANIFEST = os.environ.get("P31R_SOURCE_MANIFEST") or os.path.join(REG, "APPROVED_SOURCE_MANIFEST.tsv")
MANIFEST_COLS = ["figure", "panel", "series", "source_file", "md5", "columns_or_ids", "status", "disposition", "kind", "notes"]
MAIN_ID = re.compile(r"Figure_[1-9]\d*")
AF3_ID = re.compile(r"FigureS[1-9]\d*")
STALE_TOL_S = 2.0       # an output written by a generator is at least as new as the files it read, within filesystem granularity


def rel_v3(p):
    """Path relative to V3 with forward slashes (the spelling the render log pins)."""
    p = os.path.abspath(p).replace("\\", "/")
    v = V3.replace("\\", "/") + "/"
    return p[len(v):] if p.startswith(v) else p


def active_errors(a):
    """Schema of the active figure registry. Pure; returns a list of messages (empty = valid)."""
    e = []
    if not isinstance(a, dict) or a.get("schema") != 1:
        return ["active registry: 'schema' must be 1"]
    figs = a.get("figures")
    if not isinstance(figs, list) or not figs:
        return ["active registry: 'figures' must be a non-empty list"]
    seen, last = set(), 0
    for f in figs:
        if not isinstance(f, dict):
            e.append("figure entry is not an object: %r" % (f,))
            continue
        i, kind = f.get("id"), f.get("kind")
        if kind not in ("main", "af3"):
            e.append("%s: kind must be 'main' or 'af3'" % i)
        elif not (MAIN_ID if kind == "main" else AF3_ID).fullmatch(str(i)):
            e.append("%s: id does not fit kind %s (Figure_<n> for main, FigureS<n> for af3)" % (i, kind))
        if i in seen:
            e.append("%s: listed twice" % i)
        seen.add(i)
        o = f.get("order")
        if isinstance(o, bool) or not isinstance(o, int) or o <= last:
            e.append("%s: 'order' must be an integer larger than the previous figure's" % i)
        else:
            last = o
        pg = f.get("pages", 1)
        if isinstance(pg, bool) or not isinstance(pg, int) or pg < 1 or (kind == "main" and pg != 1):
            e.append("%s: 'pages' must be an integer >= 1 (a main figure is exactly 1)" % i)
        pn, opt = f.get("panels"), f.get("optional_panels", [])
        if pn is None:
            if kind != "af3":
                e.append("%s: a main figure must declare its planned panels (an empty list for a single-panel figure)" % i)
        elif not isinstance(pn, list) or any(not (isinstance(x, str) and re.fullmatch(r"[a-h]", x)) for x in pn) or len(set(pn)) != len(pn):
            e.append("%s: 'panels' must be a list of distinct letters a-h" % i)
        if not isinstance(opt, list) or any(not (isinstance(x, str) and re.fullmatch(r"[a-h]", x)) for x in opt) or (isinstance(pn, list) and set(opt) & set(pn)):
            e.append("%s: 'optional_panels' must be letters a-h not already in 'panels'" % i)
        for k in ("generator", "legend", "title"):
            if not isinstance(f.get(k), str) or not f.get(k).strip():
                e.append("%s: '%s' must be a non-empty string" % (i, k))
        if kind == "af3" and not (isinstance(f.get("illustrates"), list) and f["illustrates"] and all(isinstance(x, str) for x in f["illustrates"])):
            e.append("%s: an af3 figure must say which SM sections it illustrates ('illustrates')" % i)
    mains = [f.get("id") for f in figs if isinstance(f, dict) and f.get("kind") == "main"]
    afs = [f.get("id") for f in figs if isinstance(f, dict) and f.get("kind") == "af3"]
    if mains != ["Figure_%d" % k for k in range(1, len(mains) + 1)]:
        e.append("main figures must run Figure_1..Figure_n without a gap, in order: %s" % mains)
    if afs != ["FigureS%d" % k for k in range(1, len(afs) + 1)]:
        e.append("af3 figures must run FigureS1..FigureSn without a gap, in order: %s" % afs)
    kinds = [f.get("kind") for f in figs if isinstance(f, dict)]
    if "af3" in kinds and "main" in kinds[kinds.index("af3"):]:
        e.append("every main figure must come before the first af3 figure")
    if not (isinstance(a.get("retired"), list) and all(isinstance(x, str) for x in a["retired"])):
        e.append("'retired' must be a list of generator paths (empty if none)")
    e += holds_errors(a)
    for k in ("text_files", "text_ignore", "text_dirs"):
        v = a.get(k)
        if not (isinstance(v, list) and all(isinstance(x, str) and x for x in v)) or (k == "text_files" and not v):
            e.append("'%s' must be a list of strings%s" % (k, " (non-empty)" if k == "text_files" else ""))
    return e


def load_active(path=None):
    p = path or ACTIVE_JSON
    if not os.path.exists(p):
        raise RuntimeError("active figure registry missing: %s" % p)
    return json.load(io.open(p, encoding="utf-8"))


def active_units(a):
    """[(unit, figure entry, page)]: page 1 is the figure's own id, continuation pages are <id>_p2, <id>_p3 ..."""
    return [(f["id"] if k == 1 else "%s_p%d" % (f["id"], k), f, k) for f in a["figures"] for k in range(1, int(f.get("pages", 1)) + 1)]


def g13_active(a, on_disk, exists):
    """Scope of the figure set. `on_disk`: V3-relative paths of scripts/P31R_Fig*.R that exist; `exists(rel)`: file test."""
    f = ["G13 " + x for x in active_errors(a)]
    if f:
        return f
    gens = {}
    for x in a["figures"]:
        gens.setdefault(x["generator"], []).append(x["id"])
    ret = set(a["retired"])
    for g, ids in sorted(gens.items()):
        if len(ids) > 1:
            f.append("G13 generator %s is listed for %d figures (%s): one generator per figure" % (g, len(ids), ", ".join(ids)))
        if g in ret:
            f.append("G13 generator %s is both active and retired" % g)
        elif not exists(g):
            f.append("G13 active generator missing: %s" % g)
    for g in sorted(set(on_disk) - set(gens) - ret):
        f.append("G13 generator %s is on disk but neither active nor retired: list it in the active registry or under 'retired'" % g)
    return f


def g13_panels(entry, union_tags):
    """Planned panels (active registry) against the panels the rendered pages carry. (fails, notes)"""
    pl, opt, got = entry.get("panels"), set(entry.get("optional_panels", [])), set(union_tags)
    if pl is None:
        return [], ["  G13 note: %s panels not declared in the active registry (rendered: %s)" % (entry["id"], "".join(sorted(got)) or "-")]
    f = []
    if set(pl) - got:
        f.append("G13 %s: planned panels %s missing from the rendered pages (carry %s)" % (entry["id"], sorted(set(pl) - got), sorted(got)))
    if got - set(pl) - opt:
        f.append("G13 %s: rendered panels %s are not in the active registry (planned %s, optional %s)" % (entry["id"], sorted(got - set(pl) - opt), pl, sorted(opt)))
    return f, []


def g13_stale(unit, pdf_mtime, inputs):
    """An output must not predate its generator, shim or data. `inputs`: {label: mtime}."""
    late = sorted((lab for lab, t in inputs.items() if t > pdf_mtime + STALE_TOL_S))
    if not late:
        return []
    return ["G13 %s: shipped PDF is older than %d input(s) it was built from (%s%s): re-render" % (unit, len(late), ", ".join(late[:3]), ", ..." if len(late) > 3 else "")]


# ---- G12: the approved source manifest (registry/APPROVED_SOURCE_MANIFEST.tsv, written by S2-D24 task B) -----------
# One row per plotted series or printed-number source per figure panel: figure (Figure_N or AF3_S<n>; '(tiered table)', '(ledger)', '(plan)'
# for global rows), panel, series, source_file (under V3/ or absolute), md5, columns_or_ids, status (PIPE|KEY), disposition
# (APPROVED|BLOCKED|SUPERSEDED), kind (data|ledger|text|metadata|frame), notes. A version is accepted because a row APPROVES it (path + md5), whatever
# directory it sits in; a SUPERSEDED row blocks it by path + md5. A ledger row is pinned by id, not by md5 (the ledger is append-only).
def manifest_names(fid):
    """The manifest's spellings of an active figure id."""
    m = re.fullmatch(r"FigureS(\d+)", fid)
    return {fid, "AF3_S%s" % m.group(1), "Figure_S%s" % m.group(1)} if m else {fid}


def read_manifest(path=None):
    """(rows, errors). A missing file is an error, never an empty manifest."""
    p = path or MANIFEST
    if not os.path.exists(p):
        return [], ["G12 approved source manifest missing: %s" % p]
    with io.open(p, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames != MANIFEST_COLS:
            return [], ["G12 manifest columns %s != %s" % (rd.fieldnames, MANIFEST_COLS)]
        return [dict(r) for r in rd], []


def manifest_errors(rows, md5_of, figure_ids):
    """Integrity of the manifest itself. `md5_of(rel)` returns the md5 of V3-relative `rel` or None if absent."""
    f, seen, appr, sup = [], {}, {}, {}
    names = set().union(*[manifest_names(i) for i in figure_ids]) if figure_ids else set()
    for k, r in enumerate(rows, 2):
        tag = "row %d (%s %s)" % (k, r.get("figure"), (r.get("panel") or "")[:12])
        if r.get("disposition") not in ("APPROVED", "BLOCKED", "SUPERSEDED"):
            f.append("G12 manifest %s: disposition must be APPROVED, BLOCKED or SUPERSEDED, got %r" % (tag, r.get("disposition")))
        if r.get("status") not in ("PIPE", "KEY"):
            f.append("G12 manifest %s: status must be PIPE or KEY, got %r" % (tag, r.get("status")))
        if r.get("kind") not in ("data", "ledger", "text", "metadata", "frame"):
            f.append("G12 manifest %s: kind %r is not data, ledger, text, metadata or frame" % (tag, r.get("kind")))
        if not re.fullmatch(r"[0-9a-f]{32}", r.get("md5", "")):
            f.append("G12 manifest %s: md5 is not 32 hex digits" % tag)
        if not r.get("source_file"):
            f.append("G12 manifest %s: source_file missing" % tag)
        fg = r.get("figure", "")
        if fg not in names and not fg.startswith("("):
            f.append("G12 manifest %s: figure %r is not an active figure (nor a '(...)' global row)" % (tag, fg))
        if r.get("disposition") == "APPROVED":
            now = md5_of(r.get("source_file", ""))
            if now is None:
                f.append("G12 manifest %s: approved file missing: %s" % (tag, r.get("source_file")))
            elif r.get("kind") != "ledger" and now != r.get("md5"):
                f.append("G12 manifest %s: approved file changed on disk since approval: %s" % (tag, r.get("source_file")))
            key = (fg, r.get("panel"), r.get("series"))
            if key in seen:
                f.append("G12 manifest %s: series %r of %s panel %s is approved from two sources (rows %d and %d): ONE named source per series" % (tag, (r.get("series") or "")[:40], fg, r.get("panel"), seen[key], k))
            seen[key] = k
            appr.setdefault((r.get("source_file", "").replace("\\", "/"), r.get("md5")), k)
        elif r.get("disposition") == "SUPERSEDED":
            sup.setdefault((r.get("source_file", "").replace("\\", "/"), r.get("md5")), k)
    for key in sorted(set(appr) & set(sup)):
        f.append("G12 manifest: %s at %s is both APPROVED (row %d) and SUPERSEDED (row %d)" % (key[0], (key[1] or "")[:8], appr[key], sup[key]))
    return f


def same_file(a, b):
    """Two spellings of one file (a pin is the path as the generator saw it, which resolves junctions and symlinks)."""
    try:
        return os.path.realpath(v3_path(a)) == os.path.realpath(v3_path(b))
    except Exception:
        return a == b


def parse_pins(s):
    out = []
    for pin in filter(None, (s or "").split(";")):
        p, _, h = pin.rpartition("=")
        out.append((p.replace("\\", "/"), h))
    return out


def _panel_ok(row_panel, want):
    """A declared panel matches a manifest row of that panel letter ('c (optional)' counts as c) or a figure-level row ('-', 'flow'); '' = any."""
    if not want:
        return True
    rp = (row_panel or "").split()[0] if (row_panel or "").split() else ""
    return rp in (want, "-", "flow")


def g12_sources(unit, fig_id, reg, row, rows, ledger_path=None, require_sources=False):
    """Every data file a render read is an APPROVED version for this figure (path + md5). `reg`: the unit's registry JSON; `row`: its render-log row;
    `rows`: the manifest. A generator that opened its sources through af3_read() declares them as reg['sources'] = [{path, panel}]: each must be
    approved for this figure and panel, each must have been read, and nothing else may have been (require_sources: an AF3 unit must declare them)."""
    f = []
    names = manifest_names(fig_id)
    mine = [r for r in rows if r["figure"] in names]
    pins = parse_pins(row["data"]) if row else []
    unlisted, read_paths = [], []
    for p, h in pins:
        if ledger_path and same_file(p, ledger_path):
            continue                                   # claims ledger: pinned by G8, read through p31r_ledger, not a plotted series
        read_paths.append(p)
        sup = [r for r in rows if r["disposition"] == "SUPERSEDED" and r["md5"] == h and same_file(r["source_file"], p)]
        if sup:
            f.append("G12 %s: %s was read at a SUPERSEDED version (%s): blocked by the manifest" % (unit, p, (sup[0].get("notes") or "").strip()[:60]))
            continue
        here = [r for r in mine if same_file(r["source_file"], p)]
        if any(r["disposition"] == "APPROVED" and r["md5"] == h for r in here):
            continue
        appr = [r for r in here if r["disposition"] == "APPROVED"]
        if appr:
            f.append("G12 %s: %s was read at md5 %s, which is not the approved version (%s)" % (unit, p, h[:8], ", ".join(sorted({r["md5"][:8] for r in appr}))))
        elif any(r["disposition"] == "BLOCKED" for r in here):
            f.append("G12 %s: %s is BLOCKED for this figure (%s)" % (unit, p, "; ".join((r.get("notes") or "").strip()[:60] for r in here if r["disposition"] == "BLOCKED")))
        elif any(r["disposition"] == "APPROVED" and r["md5"] == h and same_file(r["source_file"], p) for r in rows):
            f.append("G12 %s: %s is approved for other figures only (%s)" % (unit, p, ", ".join(sorted({r["figure"] for r in rows if same_file(r["source_file"], p) and r["disposition"] == "APPROVED"}))))
        else:
            unlisted.append(p)
    if unlisted:   # one line per unit: a reconciliation works from the list, not from 30 lines of the same sentence
        f.append("G12 %s: %d data file(s) read but not in the approved source manifest: %s%s" % (unit, len(unlisted), "; ".join(unlisted[:3]), "; ..." if len(unlisted) > 3 else ""))
    decl = reg.get("sources") if isinstance(reg, dict) else None
    if decl is None:
        if require_sources:
            f.append("G12 %s: the registry declares no 'sources' (an AF3 unit opens its data through af3_read())" % unit)
        return f
    if not isinstance(decl, list) or not all(isinstance(d, dict) and d.get("path") for d in decl):
        return f + ["G12 %s: registry 'sources' must be a list of {path, panel}" % unit]
    for d in decl:
        ok = [r for r in mine if r["disposition"] == "APPROVED" and same_file(r["source_file"], d["path"]) and _panel_ok(r["panel"], d.get("panel", ""))]
        if not ok:
            f.append("G12 %s: declared source %s (panel %r) has no APPROVED manifest row for this figure" % (unit, d["path"], d.get("panel", "")))
    dp = [d["path"] for d in decl]
    for p in read_paths:
        if not any(same_file(p, q) for q in dp):
            f.append("G12 %s: %s was read but is not declared in the registry's sources (open data through af3_read)" % (unit, p))
    for q in dp:
        if not any(same_file(p, q) for p in read_paths):
            f.append("G12 %s: declared source %s is not among the files the render read" % (unit, q))
    return f


def _blank(v):
    return v is None or str(v).strip().lower() in ("", "na", "nan", "null", "none")


def holds_errors(a):
    e = []
    hs = a.get("holds", [])
    if not isinstance(hs, list):
        return ["'holds' must be a list"]
    for h in hs:
        if not (isinstance(h, dict) and isinstance(h.get("id"), str) and h.get("id") and isinstance(h.get("source_columns"), list) and h["source_columns"]
                and isinstance(h.get("source_value"), str) and h["source_value"] and isinstance(h.get("magnitude_columns"), list) and h["magnitude_columns"]):
            e.append("hold %r needs id, source_columns, source_value and magnitude_columns" % (h.get("id") if isinstance(h, dict) else h))
    return e


def g12_holds(holds, unit, frame_tables):
    """A source hold (H9: no magnitude from EGFR_EAS_2P2M on any published surface, main or AF3) applied to the plotted frames a unit wrote
    (`frame_tables`: {frame name: list of row dicts}). Rows from the held source must have blank magnitude cells; a frame that carries magnitude
    columns but no source column cannot be checked and fails."""
    f = []
    for h in holds:
        for name, rows in sorted(frame_tables.items()):
            if not rows:
                continue
            cols = set(rows[0])
            mag = [c for c in h["magnitude_columns"] if c in cols]
            if not mag:
                continue
            src = [c for c in h["source_columns"] if c in cols]
            if not src:
                f.append("G12 %s: frame %r carries magnitude column(s) %s but no source column %s: %s cannot be checked" % (unit, name, mag[:3], h["source_columns"], h["id"]))
                continue
            held = [r for r in rows if any(r.get(c) == h["source_value"] for c in src)]
            bad = [(c, r) for r in held for c in mag if not _blank(r.get(c))]
            if bad:
                f.append("G12 %s: frame %r holds %d magnitude value(s) (%s) on %d row(s) from the held source %s (%s)" % (unit, name, len(bad), bad[0][0], len(held), h["source_value"], h["id"]))
    return f


# ---- G14: pointers ----------------------------------------------------------------------------------------------
_HEAD = re.compile(r"(?<![A-Za-z])Fig(?:ures|ure|s\.|\.)\s+(?=S?\d)")
_NUM = re.compile(r"(S?)(\d{1,2})(?!\d)")
_LET = re.compile(r"([a-z])(?:[–-]([a-z]))?(?![a-z0-9])")
_SEP = re.compile(r"\s?,\s?(?:and\s)?|\s(?:and|&)\s")


def _fig_id(s, n):
    return ("FigureS%s" if s else "Figure_%s") % n


def parse_citations(text):
    """Every figure citation in `text`: [{fig, panels (frozenset of letters, or None for the whole figure), raw, pos, bad}].
    Case-sensitive 'Figure'/'Figures'/'Fig.'/'Figs.' only (so 'significant figures' is not a citation). Handles 'Figure 2a', '2a, b',
    '2a–c', 'Figure S1', 'Additional file 3: Figure S1', and (plural heads only) 'Figures 2d–e, 3, 6b–d' and 'Figures S1–S3'. A bare
    letter after a comma counts only when it ends a list ('2a, b)', '2a, b and c'), so 'Figure 2b, a forest' cites panel b alone."""
    out = []
    for h in _HEAD.finditer(text):
        plural = h.group(0).startswith(("Figures", "Figs"))
        pos = h.end()
        cur = None
        while True:
            m = _NUM.match(text, pos)
            if not m:
                break
            s, n = m.group(1), m.group(2)
            pos = m.end()
            ids = [_fig_id(s, n)]
            if plural:
                rg = re.match(r"[–-](S?)(\d{1,2})(?!\d)", text[pos:])
                if rg and bool(rg.group(1)) == bool(s) and int(rg.group(2)) > int(n):
                    ids = [_fig_id(s, k) for k in range(int(n), int(rg.group(2)) + 1)]
                    pos += rg.end()
            letters, bad = None, False
            ml = _LET.match(text, pos)
            if ml and len(ids) == 1:
                a_, b_ = ml.group(1), ml.group(2)
                letters = set(chr(c) for c in range(ord(a_), ord(b_ or a_) + 1)) if (not b_ or b_ >= a_) else set()
                bad = bool(b_ and b_ < a_)
                pos = ml.end()
            elif len(ids) == 1 and re.match(r"[a-z]", text[pos:pos + 1] or " "):
                bad = True      # letters glued to the number but not a clean panel token ('2ab')
            for i in ids:
                out.append(dict(fig=i, panels=frozenset(letters) if letters else None, raw=text[h.start():pos], pos=h.start(), bad=bad))
            cur = out[-1] if letters and len(ids) == 1 else None
            # continuation: bare letters after a letter item, or a further figure number after a plural head
            while True:
                sp = _SEP.match(text, pos)
                if not sp:
                    break
                nxt = sp.end()
                mb = re.match(r"([a-z])(?:[–-]([a-z]))?(?=\s?[,;)\].]|\s(?:and|&)\s|\s?$)", text[nxt:])
                if cur is not None and mb and not _NUM.match(text, nxt):
                    a_, b_ = mb.group(1), mb.group(2)
                    more = set(chr(c) for c in range(ord(a_), ord(b_ or a_) + 1))
                    cur["panels"] = frozenset(set(cur["panels"]) | more)
                    pos = nxt + mb.end()
                    continue
                break
            sp = _SEP.match(text, pos)
            nm = _NUM.match(text, sp.end()) if sp else None
            if plural and nm and (nm.group(1) or len(nm.group(2)) == 1):
                pos = sp.end()
                continue
            break
    return out


def pointer_text_scope(a, listing):
    """`listing`: V2... paths (relative to R2_rewrite) of every *.md in a['text_dirs']. Every one is either enumerated in text_files or
    matches text_ignore; an unclassified file or a missing enumerated file is a failure (a pointer gate that opens a subset of the text is
    worse than none)."""
    import fnmatch
    f = []
    listed = set(a["text_files"])
    for p in sorted(listed - set(listing)):
        f.append("G14 enumerated text file missing: %s" % p)
    for p in sorted(set(listing) - listed):
        if not any(fnmatch.fnmatch(os.path.basename(p), pat) for pat in a["text_ignore"]):
            f.append("G14 text file %s is neither enumerated in text_files nor matched by text_ignore: classify it (a gate scope is the gate)" % p)
    return f


def g14_pointers(a, texts):
    """texts: {relative name: text} for EVERY enumerated file. (fails, notes, table)"""
    figs = {f["id"]: f for f in a["figures"]}
    whole, byp, bad_f, notes = {}, {}, [], []
    n_cit = 0
    for name, t in sorted(texts.items()):
        if len(t.strip()) < 200:
            bad_f.append("G14 %s: only %d characters read (a blind reader passes every pointer check)" % (name, len(t.strip())))
        for c in parse_citations(t):
            n_cit += 1
            fe = figs.get(c["fig"])
            where = "%s: %r" % (name, c["raw"])
            if fe is None:
                bad_f.append("G14 %s cites %s, which is not an active figure" % (where, c["fig"]))
                continue
            if c["bad"]:
                bad_f.append("G14 %s is not a clean panel citation" % where)
                continue
            if c["panels"] is None:
                whole.setdefault(c["fig"], []).append(name)
                continue
            pl = fe.get("panels")
            ok = set(pl or []) | set(fe.get("optional_panels", []))
            for p in sorted(c["panels"]):
                if pl is None:
                    bad_f.append("G14 %s cites panel %s but the panels of %s are not declared yet" % (where, p, c["fig"]))
                elif p not in ok:
                    bad_f.append("G14 %s: %s has no panel %s (panels %s)" % (where, c["fig"], p, sorted(ok)))
                else:
                    byp.setdefault(c["fig"], {}).setdefault(p, []).append(name)
    if n_cit == 0:
        bad_f.append("G14 no figure citation recovered from %d text file(s): the parser or the file list is blind" % len(texts))
    table = []
    for i, fe in figs.items():
        if i not in whole and i not in byp:
            bad_f.append("G14 %s is never cited in the %d enumerated text/SM file(s)" % (i, len(texts)))
        elif i not in whole:
            un = sorted((set(fe.get("panels") or []) - set(byp.get(i, {}))))
            if un:
                notes.append("  G14 note: %s is cited only by panel letter; panels never cited: %s" % (i, un))
        table.append((i, sorted(set(whole.get(i, []))), {p: sorted(set(v)) for p, v in sorted(byp.get(i, {}).items())}))
    return bad_f, notes, table


def list_text_files(a):
    out = []
    for d in a["text_dirs"]:
        p = os.path.join(ROOT, d)
        if os.path.isdir(p):
            out += ["%s/%s" % (d, x) for x in sorted(os.listdir(p)) if x.endswith(".md") and os.path.isfile(os.path.join(p, x))]
    return out


def read_texts(a):
    return {p: io.open(os.path.join(ROOT, p), encoding="utf-8").read() for p in a["text_files"] if os.path.exists(os.path.join(ROOT, p))}


def legend_format_errors(fig_id, txt):
    """AF3 legend format (figures/af3/legends/README.md): a '## Figure S<n>' heading, then a paragraph opening '**Figure S<n>. ...'."""
    n = fig_id.replace("FigureS", "")
    f = []
    if not re.search(r"(?m)^##\s+Figure S%s\s*$" % n, txt):
        f.append("G9 legend of %s has no '## Figure S%s' heading line" % (fig_id, n))
    if not re.search(r"(?m)^\*\*Figure S%s\.\s+[^*]+\*\*" % n, txt):
        f.append("G9 legend of %s does not open with '**Figure S%s. Title.**'" % (fig_id, n))
    return f


def g9_af3_legend(entry, union_tags, txt):
    f = legend_format_errors(entry["id"], txt)
    lp = legend_panels(txt)
    tags = sorted(union_tags)
    if not (len(tags) <= 1 and not lp) and tags != lp:
        f.append("G9 legend panels %s != figure panels %s (%s, all pages)" % (lp, tags, entry["id"]))
    return f


def v3_path(rel):
    return rel if re.match(r"^[A-Za-z]:", rel) or os.path.isabs(rel) else os.path.join(V3, rel)


def md5_rel(rel):
    p = v3_path(rel)
    return md5(p) if os.path.isfile(p) else None


def active_context():
    """Everything the global checks need from disk: (active registry or None, [failures])."""
    try:
        a = load_active()
    except Exception as e:
        return None, ["G13 %s" % e]
    return a, []


def global_checks(a, rows_m, m_errs):
    """G13 scope, G12 manifest integrity, G14 pointers. (fails, notes, pointer table)"""
    gdir = os.path.join(FIGS, "scripts")
    on_disk = ["R2_rewrite/figures/scripts/" + x for x in sorted(os.listdir(gdir)) if re.fullmatch(r"P31R_Fig.*\.R", x)]
    f = g13_active(a, on_disk, lambda rel: os.path.isfile(v3_path(rel)))
    f += list(m_errs)
    f += manifest_errors(rows_m, md5_rel, {x["id"] for x in a["figures"]}) if not m_errs else []
    notes = []
    stems = {os.path.splitext(x)[0] for d in (OUT, REG, LOGD) if os.path.isdir(d) for x in os.listdir(d)}
    act_units = {u for u, _, _ in active_units(a)}
    orphan = sorted(s for s in stems if re.fullmatch(r"Figure(?:_\d+|S\d+(?:_p\d+)?)", s) and s not in act_units)
    if orphan:
        notes.append("  G13 note: outputs/registries/render logs not in the active registry (ignored, not shipped): %s" % ", ".join(orphan))
    f += pointer_text_scope(a, list_text_files(a))
    pf, pn, table = g14_pointers(a, read_texts(a))
    return f + pf, notes + pn, table


# ---------------------------------------------------------------- driver
def render_log():
    out = {}
    if not os.path.isdir(LOGD):
        return out
    for f in os.listdir(LOGD):
        if f.endswith(".tsv"):
            with io.open(os.path.join(LOGD, f), encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            if len(rows) != 1 or rows[0]["figure"] + ".tsv" != f:
                raise RuntimeError("render log %s must hold exactly its own figure's row" % f)
            out[rows[0]["figure"]] = rows[0]
    return out


def run(figs=None):
    """figs=None: every unit of the active registry (main Figure_1-5, AF3 FigureS<n> and continuation pages) plus the global checks
    (G13 scope, G12 manifest integrity, G14 pointers). Explicit unit ids run only those units (the global checks are skipped)."""
    L = PG.ledger()
    log = render_log()
    cites = cited_panels()
    total = 0
    active, actf = active_context()
    if active is not None:
        ae = active_errors(active)
        if ae:
            actf += ["G13 " + x for x in ae]
            active = None
    units = {u: (e, k) for u, e, k in active_units(active)} if active else {}
    full = figs is None
    if full:
        figs = list(units) if active else list(EXPECTED)
    if actf:
        print("ACTIVE REGISTRY: %d failure(s)" % len(actf))
        for x in actf:
            print("  FAIL " + x)
        total += len(actf)
    rows_m, _m_errs = read_manifest()
    ledger_path = PG.LEDGER
    tags_by_fig, done_units = {}, set()
    for fig in figs:
        entry, page = units.get(fig, (None, 1))
        is_af3 = entry is not None and entry["kind"] == "af3"
        pdf = os.path.join(OUT, fig + ".pdf")
        fails, notes = [], []
        if not os.path.exists(pdf):
            print("%s: FAIL G0 missing %s" % (fig, pdf))
            total += 1
            continue
        try:
            ex = extract(pdf)
        except Exception as e:  # the reader must raise, never return empty
            print("%s: FAIL G0 reader error %s" % (fig, e))
            total += 1
            continue
        nch = sum(len(re.sub(r"\s", "", s["text"])) for s in ex["spans"])
        n2 = second_reader_chars(pdf)
        if not (40 <= nch <= 12000):
            fails.append("G0 implausible character count %d (band 40-12,000)" % nch)
        if n2 == 0 or abs(nch - n2) > 0.15 * max(nch, n2):
            fails.append("G0 readers disagree: PyMuPDF %d chars vs pdfminer %d" % (nch, n2))
        tags = panel_tags(ex)
        regp = os.path.join(REG, fig + ".json")
        reg = json.load(io.open(regp, encoding="utf-8")) if os.path.exists(regp) else None
        if reg is None:
            fails.append("G7 registry missing: %s" % os.path.relpath(regp, ROOT))
            reg = {}
        want_tags = reg.get("panels")
        if want_tags is not None and tags != sorted(want_tags):
            fails.append("G0 reader recovered panel tags %s, registry says %s" % (tags, sorted(want_tags)))
        fails += g1_page(ex) + g2_sizes(ex) + g3_glyphs(ex) + g4_offpage(ex) + g5_overlap(ex) + g6_vocab(ex)
        try:
            fails += g7_numbers(ex, reg, L)
        except Exception as e_:   # a malformed registry must not hide the other figures' results
            fails.append("G7 registry could not be checked: %s: %s" % (type(e_).__name__, e_))
        fails += g8_provenance(fig, log.get(fig))
        if is_af3:
            # AF3 (decision F9): G11 blocking, no G1b (a main-text print rule), no G10 (G14 reads the R2 text); the legend is checked per figure below
            gf, gn = route_geometry(g11_strokes(ex), True)
            fails += gf
            notes += gn
        else:
            fails += g9_legend(fig, ex, active)
            lgp = main_legend_path(fig, active)
            geo = g1b_height(ex, legend_caption_text(io.open(lgp, encoding="utf-8").read()) if os.path.exists(lgp) else "") + g11_strokes(ex)
            gf, gn = route_geometry(geo, HEIGHT_BLOCKING)
            fails += gf
            notes += gn
            if active is None:   # active registry loaded: G14 replaces G10
                f10, n10 = g10_text(fig, ex, cites)
                fails += f10
                notes += n10
        if entry is not None:
            tags_by_fig.setdefault(entry["id"], set()).update(tags)
            done_units.add(fig)
            row = log.get(fig)
            if row is not None:
                fails += g12_sources(fig, entry["id"], reg, row, rows_m, ledger_path, require_sources=is_af3)
                frames = {}
                for nm in reg.get("frames", []) if isinstance(reg.get("frames"), list) else []:
                    fp = os.path.join(AF3, "source_data", "%s__%s.csv" % (fig, nm))
                    if os.path.exists(fp):
                        with io.open(fp, encoding="utf-8-sig", newline="") as fh:
                            frames[nm] = list(csv.DictReader(fh))
                    else:
                        fails.append("G12 %s: plotted frame %s is missing" % (fig, os.path.relpath(fp, ROOT)))
                fails += g12_holds(active.get("holds", []), fig, frames)
                inputs = {}
                for lab, rel in [("generator", row["generator"]), ("common helper", os.path.join(HERE, "P31R_fig_common.R"))] + \
                                [("house dependency " + os.path.basename(p), p) for p, _ in parse_pins(row["deps"])] + \
                                [("data file " + os.path.basename(p), p) for p, _ in parse_pins(row["data"])]:
                    q = v3_path(rel)
                    if os.path.exists(q):
                        inputs[lab] = os.path.getmtime(q)
                fails += g13_stale(fig, os.path.getmtime(pdf), inputs)
        print("%s: %.0f x %.0f mm, %d chars (pdfminer %d), panels %s, %d failure(s)"
              % (fig, ex["w_mm"], ex["h_mm"], nch, n2, "".join(tags) or "-", len(fails)))
        for n in notes:
            print(n)
        for x in fails:
            print("  FAIL " + x)
        total += len(fails)
    # figure level (all pages of a figure): planned panels vs rendered panels; an AF3 legend names the union of its pages' panels
    if active is not None:
        for e in active["figures"]:
            mine = [u for u, ee, _ in active_units(active) if ee["id"] == e["id"]]
            if not mine or not all(u in done_units for u in mine):
                continue
            ff, fn = g13_panels(e, tags_by_fig.get(e["id"], set()))
            if e["kind"] == "af3":
                lp = v3_path(e["legend"])
                if os.path.exists(lp):
                    ff += g9_af3_legend(e, tags_by_fig.get(e["id"], set()), io.open(lp, encoding="utf-8").read())
                else:
                    ff.append("G9 legend file missing: %s" % e["legend"])
            print("%s (figure level, %d page(s)): %d failure(s)" % (e["id"], len(mine), len(ff)))
            for n in fn:
                print(n)
            for x in ff:
                print("  FAIL " + x)
            total += len(ff)
    if full and active is not None:
        gf, gn, _t = global_checks(active, *read_manifest())
        print("ACTIVE REGISTRY, SOURCE MANIFEST, POINTERS: %d failure(s)" % len(gf))
        for n in gn:
            print(n)
        for x in gf:
            print("  FAIL " + x)
        total += len(gf)
    print("FIGURE GATE %s: %d failure(s) over %d figure(s)" % ("RED" if total else "GREEN", total, len(figs)))
    return 1 if total else 0


def cli_global(pointers_only):
    """--active: G13 scope + G12 manifest integrity + G14 pointers.  --pointers: G14 only, with the citation table."""
    a, f0 = active_context()
    if a is not None and active_errors(a):
        f0, a = ["G13 " + x for x in active_errors(a)], None
    if a is None:
        for x in f0:
            print("  FAIL " + x)
        print("ACTIVE REGISTRY RED: %d failure(s)" % len(f0))
        return 1
    rows, me = read_manifest()
    if pointers_only:
        fails = pointer_text_scope(a, list_text_files(a))
        pf, notes, table = g14_pointers(a, read_texts(a))
        fails += pf
    else:
        fails, notes, table = global_checks(a, rows, me)
    print("%-10s %-4s %-46s %s" % ("figure", "kind", "cited whole in", "cited by panel"))
    kinds = {x["id"]: x["kind"] for x in a["figures"]}
    for i, w, p in table:
        print("%-10s %-4s %-46s %s" % (i, kinds[i], ", ".join(os.path.basename(x) for x in w) or "-",
                                      "; ".join("%s: %s" % (k, ",".join(os.path.basename(x) for x in v)) for k, v in p.items()) or "-"))
    for n in notes:
        print(n)
    for x in fails:
        print("  FAIL " + x)
    print("%s %s: %d failure(s) over %d enumerated text file(s)" % ("POINTERS" if pointers_only else "ACTIVE REGISTRY", "RED" if fails else "GREEN", len(fails), len(a["text_files"])))
    return 1 if fails else 0


def af3_cases():
    """Self-test cases for G9 (AF3 legend), G12, G13 and G14. Each defect is injected into a COPY of a clean fixture and must be named by the
    intended rule (the needle) with the intended prefix; each clean fixture must pass. Returns [(label, output, wanted prefix or None)]."""
    import copy
    nd = lambda out, needle: [x for x in out if needle in x]
    fig = lambda i, k, o, panels, **kw: dict(dict(id=i, kind=k, order=o, pages=1, panels=panels, title="t", generator="g/P31R_%s.R" % i, legend="l/%s.md" % i,
                                                  **({"illustrates": ["S-R11"]} if k == "af3" else {})), **kw)
    A = dict(schema=1, figures=[fig("Figure_1", "main", 1, []), fig("Figure_2", "main", 2, ["a", "b", "c"]),
                                fig("FigureS1", "af3", 3, None), fig("FigureS2", "af3", 4, ["a", "b"])],
             retired=["g/P31R_Fig9.R"], text_files=["sections_R2/R.md"], text_ignore=["*_pre*.md"], text_dirs=["sections_R2"])
    gens = [x["generator"] for x in A["figures"]]

    def mut(fn):
        b = copy.deepcopy(A)
        fn(b)
        return b
    allexist = lambda rel: True
    cs = []
    # ---- G13: the active list is the scope
    cs += [("G13 ok", g13_active(A, gens + ["g/P31R_Fig9.R"], allexist), None),
           ("G13 dup id", nd(g13_active(mut(lambda b: b["figures"].append(fig("Figure_2", "main", 9, ["a"]))), gens, allexist), "listed twice"), "G13"),
           ("G13 kind/id", nd(g13_active(mut(lambda b: b["figures"][1].update(kind="af3")), gens, allexist), "does not fit kind"), "G13"),
           ("G13 main 2 pages", nd(g13_active(mut(lambda b: b["figures"][1].update(pages=2)), gens, allexist), "exactly 1"), "G13"),
           ("G13 main no panels", nd(g13_active(mut(lambda b: b["figures"][1].update(panels=None)), gens, allexist), "must declare"), "G13"),
           ("G13 af3 no illustrates", nd(g13_active(mut(lambda b: b["figures"][2].pop("illustrates")), gens, allexist), "illustrates"), "G13"),
           ("G13 gap", nd(g13_active(mut(lambda b: b["figures"].pop(0)), gens, allexist), "without a gap"), "G13"),
           ("G13 order", nd(g13_active(mut(lambda b: b["figures"][1].update(order=1)), gens, allexist), "larger than the previous"), "G13"),
           ("G13 main after af3", nd(g13_active(mut(lambda b: b["figures"].__setitem__(slice(None), [b["figures"][i] for i in (0, 2, 1, 3)])), gens, allexist), "before the first af3"), "G13"),
           ("G13 optional overlaps", nd(g13_active(mut(lambda b: b["figures"][1].update(optional_panels=["a"])), gens, allexist), "optional_panels"), "G13"),
           ("G13 generator missing", nd(g13_active(A, gens, lambda rel: rel != "g/P31R_Figure_2.R"), "active generator missing"), "G13"),
           ("G13 unlisted generator", nd(g13_active(A, gens + ["g/P31R_Fig7.R"], allexist), "neither active nor retired"), "G13"),
           ("G13 retired not picked up", g13_active(A, gens + ["g/P31R_Fig9.R"], allexist), None),
           ("G13 active and retired", nd(g13_active(mut(lambda b: b["retired"].append("g/P31R_Figure_1.R")), gens, allexist), "both active and retired"), "G13"),
           ("G13 shared generator", nd(g13_active(mut(lambda b: b["figures"][1].update(generator=b["figures"][0]["generator"])), gens, allexist), "one generator per figure"), "G13"),
           ("G13 panels ok", g13_panels(A["figures"][1], {"a", "b", "c"})[0], None),
           ("G13 panel missing", nd(g13_panels(A["figures"][1], {"a", "b"})[0], "missing from the rendered"), "G13"),
           ("G13 panel extra", nd(g13_panels(A["figures"][1], {"a", "b", "c", "d"})[0], "not in the active registry"), "G13"),
           ("G13 optional ok", g13_panels(dict(A["figures"][1], panels=["a", "b"], optional_panels=["c"]), {"a", "b", "c"})[0], None),
           ("G13 undeclared is a note", g13_panels(A["figures"][2], {"a"})[0], None),
           ("G13 stale", nd(g13_stale("FigureS1", 1000.0, {"generator": 2000.0, "data file x": 900.0}), "older than 1 input"), "G13"),
           ("G13 stale ok", g13_stale("FigureS1", 1000.0, {"generator": 999.0, "data file x": 1001.0}), None)]
    # ---- G12: the approved source manifest (registry/APPROVED_SOURCE_MANIFEST.tsv schema)
    LED = "R1/ledger.csv"
    M = lambda **k: dict(dict(figure="Figure_2", panel="a", series="s1", source_file="d/a.csv", md5="a" * 32, columns_or_ids="x", status="PIPE",
                              disposition="APPROVED", kind="data", notes=""), **k)
    pins = lambda *x: dict(data=";".join("%s=%s" % (p, h) for p, h in x))
    ROW = pins((LED, "f" * 32), ("d/a.csv", "a" * 32))
    g = lambda reg=None, row=ROW, rows=None, fid="Figure_2", req=False: g12_sources(fid, fid, reg or {}, row, rows if rows is not None else [M()], LED, req)
    ids = ["Figure_2", "Figure_3", "FigureS1"]
    md5_ok = lambda rel: "a" * 32
    SRC = lambda *p: {"sources": [dict(path=x, panel=y) for x, y in p]}
    cs += [("G12 ok", g(), None),
           ("G12 unlisted pin", nd(g(row=pins((LED, "f" * 32), ("d/a.csv", "a" * 32), ("d/b.csv", "b" * 32))), "not in the approved source manifest"), "G12"),
           ("G12 unapproved version", nd(g(row=pins((LED, "f" * 32), ("d/a.csv", "c" * 32))), "not the approved version"), "G12"),
           ("G12 superseded row", nd(g(rows=[M(disposition="SUPERSEDED")]), "SUPERSEDED"), "G12"),
           ("G12 superseded by a global row", nd(g(rows=[M(), M(figure="(tiered table)", panel="-", disposition="SUPERSEDED")]), "SUPERSEDED"), "G12"),
           ("G12 blocked source", nd(g(rows=[M(disposition="BLOCKED", notes="no per-signal SE")]), "BLOCKED for this figure"), "G12"),
           ("G12 blocked beside approved is accepted", g(rows=[M(disposition="BLOCKED", panel="c"), M()]), None),
           ("G12 approved for another figure only", nd(g(rows=[M(figure="Figure_3")]), "other figures only"), "G12"),
           ("G12 AF3 manifest spelling", g12_sources("FigureS1_p2", "FigureS1", {}, ROW, [M(figure="AF3_S1", panel="-")], LED), None),
           ("G12 ledger not a series", g(row=pins((LED, "e" * 32), ("d/a.csv", "a" * 32))), None),
           ("G12 declared sources ok", g(reg=SRC(("d/a.csv", "a")), req=True), None),
           ("G12 AF3 must declare", nd(g(req=True), "declares no 'sources'"), "G12"),
           ("G12 declared not approved", nd(g(reg=SRC(("d/a.csv", "a"), ("d/z.csv", "a")), req=True), "has no APPROVED manifest row"), "G12"),
           ("G12 declared wrong panel", nd(g(reg=SRC(("d/a.csv", "b")), req=True), "has no APPROVED manifest row"), "G12"),
           ("G12 figure-level row serves a panel", g(reg=SRC(("d/a.csv", "b")), rows=[M(panel="-")], req=True), None),
           ("G12 'c (optional)' serves panel c", g(reg=SRC(("d/a.csv", "c")), rows=[M(panel="c (optional)")], req=True), None),
           ("G12 read but not declared", nd(g(reg=SRC(("d/a.csv", "a")), rows=[M(), M(source_file="d/b.csv", series="s2")], row=pins((LED, "f" * 32), ("d/a.csv", "a" * 32), ("d/b.csv", "b" * 32)) | {}), "not declared"), "G12"),
           ("G12 declared not read", nd(g(reg=SRC(("d/a.csv", "a")), row=pins((LED, "f" * 32))), "not among the files the render read"), "G12"),
           ("G12 manifest ok", manifest_errors([M(), M(figure="(ledger)", panel="-", kind="ledger", source_file="R1/ledger.csv", md5="9" * 32)], md5_ok, ids), None),
           ("G12 manifest dup series", nd(manifest_errors([M(), M(source_file="d/b.csv")], md5_ok, ids), "ONE named source per series"), "G12"),
           ("G12 manifest same series other panel ok", manifest_errors([M(), M(panel="b", source_file="d/b.csv")], md5_ok, ids), None),
           ("G12 manifest changed", nd(manifest_errors([M()], lambda rel: "b" * 32, ids), "changed on disk"), "G12"),
           ("G12 manifest ledger md5 moves", manifest_errors([M(kind="ledger")], lambda rel: "b" * 32, ids), None),
           ("G12 manifest missing file", nd(manifest_errors([M()], lambda rel: None, ids), "approved file missing"), "G12"),
           ("G12 manifest bad disposition", nd(manifest_errors([M(disposition="OK")], md5_ok, ids), "disposition must be"), "G12"),
           ("G12 manifest bad status", nd(manifest_errors([M(status="draft")], md5_ok, ids), "status must be"), "G12"),
           ("G12 manifest bad kind", nd(manifest_errors([M(kind="png")], md5_ok, ids), "kind"), "G12"),
           ("G12 manifest bad md5", nd(manifest_errors([M(md5="xyz")], md5_ok, ids), "md5 is not 32 hex"), "G12"),
           ("G12 manifest unknown figure", nd(manifest_errors([M(figure="Figure_9")], md5_ok, ids), "not an active figure"), "G12"),
           ("G12 manifest AF3 spelling", manifest_errors([M(figure="AF3_S1", panel="-")], md5_ok, ids), None),
           ("G12 manifest approved and superseded", nd(manifest_errors([M(), M(figure="(plan)", disposition="SUPERSEDED")], md5_ok, ids), "both APPROVED"), "G12"),
           ("G12 manifest superseded needs no file", manifest_errors([M(disposition="SUPERSEDED")], lambda rel: None, ids), None),
           ("G12 manifest file missing", nd(read_manifest(os.path.join(HERE, "no_such_manifest.tsv"))[1], "manifest missing"), "G12")]
    H9 = dict(id="H9", source_columns=["eas_outcome_source", "eas_source"], source_value="EGFR_EAS_2P2M", magnitude_columns=["beta_mr_eas", "p_eas"], note="x")
    fr = lambda *rows: {"f": [dict(eas_outcome_source=s, beta_mr_eas=b, p_eas=p) for s, b, p in rows]}
    cs += [("G12 hold ok (held rows blank)", g12_holds([H9], "FigureS6", fr(("EGFR_EAS_2P2M", "", "NA"), ("T2DGGI_EAS", "0.2", "0.01"))), None),
           ("G12 hold violated", nd(g12_holds([H9], "FigureS6", fr(("EGFR_EAS_2P2M", "0.2", ""), ("T2DGGI_EAS", "0.2", "0.01"))), "held source"), "G12"),
           ("G12 hold, no source column", nd(g12_holds([H9], "FigureS6", {"f": [dict(beta_mr_eas="0.2")]}), "cannot be checked"), "G12"),
           ("G12 hold, frame without magnitudes", g12_holds([H9], "FigureS6", {"f": [dict(eas_outcome_source="EGFR_EAS_2P2M", sign="same")]}), None),
           ("G12 hold, other source column name", nd(g12_holds([H9], "FigureS6", {"f": [dict(eas_source="EGFR_EAS_2P2M", beta_mr_eas="0.2")]}), "held source"), "G12"),
           ("G12 hold schema", nd(["G12 " + x for x in holds_errors({"holds": [dict(id="H9")]})], "needs id"), "G12")]
    # ---- G14: pointers
    PAD = "x" * 220 + " "
    T = lambda s: {"sections_R2/R.md": PAD + s}
    P = lambda s, a=A: g14_pointers(a, T(s))
    ok_txt = "Figure 1; Figure 2a, b and c; Additional file 3: Figure S1; Figure S2."
    CIT = lambda s: sorted((c["fig"], tuple(sorted(c["panels"] or ()))) for c in parse_citations(s) if not c["bad"])
    cs += [("G14 ok", P(ok_txt)[0], None),
           ("G14 unresolved panel", nd(P(ok_txt + " Figure 2e.")[0], "has no panel e"), "G14"),
           ("G14 non-active main", nd(P(ok_txt + " Figure 7.")[0], "not an active figure"), "G14"),
           ("G14 non-active af3", nd(P(ok_txt + " Figure S9.")[0], "not an active figure"), "G14"),
           ("G14 uncited figure", nd(P("Figure 1; Figure 2a, b, c; Figure S1.")[0], "FigureS2 is never cited"), "G14"),
           ("G14 mutation: citation deleted", nd(P(ok_txt.replace("Additional file 3: Figure S1;", ""))[0], "FigureS1 is never cited"), "G14"),
           ("G14 mutation: panel renamed", nd(P(ok_txt.replace("2a", "2z"))[0], "has no panel z"), "G14"),
           ("G14 whole covers panels", [x for x in P("Figure 1; Figure 2; Figure S1; Figure S2")[0]] + P("Figure 1; Figure 2; Figure S1; Figure S2")[1], None),
           ("G14 letters only is a note", [] if P("Figure 1; Figure 2a; Figure S1; Figure S2")[0] == [] and nd(P("Figure 1; Figure 2a; Figure S1; Figure S2")[1], "never cited: ['b', 'c']") else ["G14 note missing"], None),
           ("G14 panel of undeclared af3", nd(P(ok_txt + " Figure S1a.")[0], "not declared yet"), "G14"),
           ("G14 declared af3 panel ok", P(ok_txt + " Figure S2b.")[0], None),
           ("G14 malformed", nd(P(ok_txt + " Figure 2ab.")[0], "not a clean panel citation"), "G14"),
           ("G14 short text", nd(g14_pointers(A, {"sections_R2/R.md": "Figure 1"})[0], "characters read"), "G14"),
           ("G14 blind parser", nd(g14_pointers(A, T("no pointers here"))[0], "blind"), "G14"),
           ("G14 parse range", [] if CIT("Figure 2a–c") == [("Figure_2", ("a", "b", "c"))] else ["G14 parse"], None),
           ("G14 parse list", [] if CIT("Figure 2a, b)") == [("Figure_2", ("a", "b"))] else ["G14 parse"], None),
           ("G14 parse and", [] if CIT("Figure 2a and b") == [("Figure_2", ("a", "b"))] else ["G14 parse"], None),
           ("G14 parse article", [] if CIT("Figure 2b, a forest") == [("Figure_2", ("b",))] else ["G14 parse"], None),
           ("G14 parse plural", [] if CIT("Figures 2d–e, 3, 6b–d") == [("Figure_2", ("d", "e")), ("Figure_3", ()), ("Figure_6", ("b", "c", "d"))] else ["G14 parse %r" % CIT("Figures 2d–e, 3, 6b–d")], None),
           ("G14 parse S range", [] if CIT("Figures S1–S3") == [("FigureS1", ()), ("FigureS2", ()), ("FigureS3", ())] else ["G14 parse"], None),
           ("G14 parse AF3 prefix", [] if CIT("(Additional file 3: Figure S1; Additional file 1: Table S5)") == [("FigureS1", ())] else ["G14 parse"], None),
           ("G14 parse not a citation", [] if CIT("three significant figures; the figure 2 and 3; Table 2") == [] else ["G14 parse"], None),
           ("G14 parse singular no list", [] if CIT("Figure 3, 45 of 52") == [("Figure_3", ())] else ["G14 parse"], None),
           ("G14 parse Fig.", [] if CIT("(Fig. 4b)") == [("Figure_4", ("b",))] else ["G14 parse"], None)]
    lst = ["sections_R2/R.md", "sections_R2/R_pre1.md", "sections_R2/New.md"]
    cs += [("G14 scope ok", pointer_text_scope(A, ["sections_R2/R.md", "sections_R2/R_pretagfix.md"]), None),
           ("G14 scope unclassified", nd(pointer_text_scope(A, lst), "New.md"), "G14"),
           ("G14 scope missing", nd(pointer_text_scope(A, ["sections_R2/R_pre1.md"]), "enumerated text file missing"), "G14")]
    # ---- G9 for an AF3 legend
    LG = lambda n, body: "## Figure S%d\n\n%s\n" % (n, body)
    cs += [("G9 af3 ok", g9_af3_legend(A["figures"][3], {"a", "b"}, LG(2, "**Figure S2. Title.** **a**, x; **b**, y. Abbreviations.")), None),
           ("G9 af3 heading", nd(g9_af3_legend(A["figures"][3], {"a", "b"}, "**Figure S2. Title.** **a**, x; **b**, y."), "heading"), "G9"),
           ("G9 af3 opener", nd(g9_af3_legend(A["figures"][3], {"a", "b"}, LG(2, "Title. **a**, x; **b**, y.")), "does not open"), "G9"),
           ("G9 af3 panels", nd(g9_af3_legend(A["figures"][3], {"a", "b"}, LG(2, "**Figure S2. Title.** **a**, x.")), "legend panels"), "G9"),
           ("G9 af3 single panel", g9_af3_legend(A["figures"][2], set(), LG(1, "**Figure S1. Title.** Text.")), None)]
    # ---- the real files: the readers must recover something (a blind reader passes every absence check) and a mutation of the REAL text must trip
    real = load_active()
    assert not active_errors(real), active_errors(real)
    texts = read_texts(real)
    assert len(texts) == len(real["text_files"]), "an enumerated text file could not be opened: %s" % sorted(set(real["text_files"]) - set(texts))
    found = [(n, c) for n, t in sorted(texts.items()) for c in parse_citations(t)]
    assert found, "the pointer parser recovered no figure citation from the real R2 text"
    n0, c0 = found[0]
    t0 = texts[n0]
    assert t0[c0["pos"]:c0["pos"] + len(c0["raw"])] == c0["raw"], "mutation target not at its recorded position"
    bad = dict(texts)
    bad[n0] = t0[:c0["pos"]] + "Figure 9" + t0[c0["pos"] + len(c0["raw"]):]
    cs.append(("G14 real-text mutation: Figure 9", nd(g14_pointers(real, bad)[0], "Figure_9"), "G14"))
    fig0 = c0["fig"]
    kill = {n: t for n, t in texts.items()}
    for n, t in texts.items():
        for c in sorted((c for c in parse_citations(t) if c["fig"] == fig0), key=lambda c: -c["pos"]):
            kill[n] = kill[n][:c["pos"]] + "[x]" + kill[n][c["pos"] + len(c["raw"]):]
    cs.append(("G14 real-text mutation: %s uncited" % fig0, nd(g14_pointers(real, kill)[0], "%s is never cited" % fig0), "G14"))
    return cs


def selftest():
    """Every check is fed a defect it must name (the failure must come from THAT check, and for G7 from the named rule),
    and a clean case it must pass. Synthetic spans mimic cairo geometry: ~0.4 pt between spans, scripts shifted 0.3 em."""
    L = PG.ledger()
    R = (0, 0, 481.89, 600.0)  # 170 x 211.7 mm

    def ex(lines=(), spans=(), w=170.0, h=200.0):
        return dict(w_mm=w, h_mm=h, rect=R, lines=list(lines), spans=list(spans))

    def sp(t, size=9.0, font="ArialMT", flags=0, bbox=(10, 10, 60, 20)):
        n = max(1, len(t))
        cw = (bbox[2] - bbox[0]) / n
        return dict(text=t, size=size, font=font, flags=flags, bbox=bbox,
                    chars=[(c, (bbox[0] + i * cw, bbox[1], bbox[0] + (i + 1) * cw, bbox[3])) for i, c in enumerate(t)])

    def ln(t, bbox=(10, 10, 60, 20)):
        return dict(text=t, bbox=bbox)

    def line(*spans):   # a real line: assembled text + its spans
        return dict(text=assemble(list(spans)), bbox=(spans[0]["bbox"][0], 10, spans[-1]["bbox"][2], 20), spans=list(spans))

    def g7(lines, reg, rule):
        return [x for x in g7_numbers(ex(lines=lines), reg, L) if x.startswith(rule)]

    def only(lines, reg, needle):
        return [x for x in g7_numbers(ex(lines=lines), reg, L) if needle in x]

    T = lambda txt: [ln(txt)]
    KC, RC = (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)   # black, red (not R: R is the page rect above)

    def S(t, w, color=KC, fill=None, rw=10.0, rh=0.0):   # a drawing record as extract() writes it
        return dict(type=t, width=w, color=color if t != "f" else None, fill=fill, rw=rw, rh=rh)
    tick = lambda t: {"text": t, "why": "axis_tick: x", "axis": "x"}
    sup9 = line(sp("p = 2.9 × 10", 7.0, bbox=(10, 10, 60, 20)), sp("−15", 5.6, bbox=(60.4, 8, 70, 17.9)))
    sub = line(sp("−log", 9.0, bbox=(10, 10, 30, 20)), sp("10", 7.2, bbox=(30.4, 13, 38, 23)), sp("(p)", 9.0, bbox=(38.4, 10, 50, 20)))
    small_same = line(sp("Tests", 9.0, bbox=(10, 10, 40, 20)), sp("12", 7.0, bbox=(40.4, 12, 47, 20)))
    cases = [
        ("G1", g1_page(ex(w=183.0)), "G1"), ("G1 ok", g1_page(ex()), None), ("G1 tall", g1_page(ex(h=240.0)), "G1"),
        # G1b: a 224-mm image with a 250-word legend overruns; a 150-mm image with the same legend fits; a 210-mm image
        # fits with a 20-word legend but not with 250 words (the legend, not the image, decides)
        ("G1b over", g1b_height(ex(h=224.0), "w " * 750), "G1b"), ("G1b ok", g1b_height(ex(h=150.0), "w " * 750), None),
        ("G1b legend decides", g1b_height(ex(h=210.0), "x" * 1500), "G1b"), ("G1b short ok", g1b_height(ex(h=210.0), "x" * 120), None),
        ("G1b no legend", g1b_height(ex(h=150.0), ""), "G1b"),
        # G11: visible width at print scale. KC = black, RC = red; rw/rh = path bounding box (pt)
        ("G11 thin line", g11_strokes(dict(w_mm=170, strokes=[S("s", 0.213), S("s", 0.5)])), "G11"),
        ("G11 hairline", g11_strokes(dict(w_mm=170, strokes=[S("s", 0.0)])), "G11"),
        ("G11 thin fill bar", g11_strokes(dict(w_mm=170, strokes=[S("f", None, fill=KC, rw=0.09, rh=6.0)])), "G11"),
        ("G11 outline line", g11_strokes(dict(w_mm=170, strokes=[S("fs", 0.2, color=RC, fill=KC, rw=3.0, rh=3.0)])), "G11"),
        ("G11 thin exon", g11_strokes(dict(w_mm=170, strokes=[S("fs", 0.1, color=KC, fill=KC, rw=0.06, rh=5.0)])), "G11"),
        ("G11 exon ok", g11_strokes(dict(w_mm=170, strokes=[S("fs", 0.2134, color=KC, fill=KC, rw=0.06, rh=5.0)])), None),
        ("G11 markers ok", g11_strokes(dict(w_mm=170, strokes=[S("fs", 0.0075, color=KC, fill=KC, rw=2.0, rh=2.0), S("f", None, fill=KC, rw=2.0, rh=2.0),
                                                              S("s", 0.25)])), None),
        ("G11 route warn", route_geometry(["G1b x"], False)[0], None), ("G11 route block", route_geometry(["G1b x"], True)[0], "G1b"),
        ("G2", g2_sizes(ex(spans=[sp("PCSK9", size=8.0)])), "G2"), ("G2 ok", g2_sizes(ex(spans=[sp("x", 7.0), sp("y", 9.0), sp("a", 11.0)])), None),
        ("G2 sup ok", g2_sizes(ex(lines=[sup9])), None), ("G2 sub ok", g2_sizes(ex(lines=[sub])), None),
        ("G2 lone small", g2_sizes(ex(spans=[sp("−15", size=5.6)])), "G2"),
        ("G2 tiny script", g2_sizes(ex(lines=[line(sp("× 10", 7.0, bbox=(10, 10, 30, 20)), sp("−9", 4.0, bbox=(30.4, 8, 36, 16)))])), "G2"),
        ("G3", g3_glyphs(ex(spans=[sp("10\u207b\u00b3")])), "G3"), ("G3 ok", g3_glyphs(ex(spans=[sp("\u2212 \u2265 \u00d7 \u03b2")])), None),
        ("G3 font", g3_glyphs(ex(spans=[sp("x", font="Helvetica")])), "G3"),
        ("G4", g4_offpage(ex(spans=[sp("edge", bbox=(470, 10, 500, 20))])), "G4"), ("G4 ok", g4_offpage(ex(spans=[sp("in")])), None),
        ("G5", g5_overlap(ex(lines=[ln("PAM", (10, 10, 40, 20)), ln("PCSK9", (12, 12, 50, 22))])), "G5"),
        ("G5 ok", g5_overlap(ex(lines=[ln("PAM", (10, 10, 40, 20)), ln("PCSK9", (10, 30, 50, 40))])), None),
        ("G6", g6_vocab(ex(lines=T("CKD endpoint"))), "G6"), ("G6 rep", g6_vocab(ex(lines=T("deCODE replication"))), "G6"),
        ("G6 ok", g6_vocab(ex(lines=T("CKD outcome"))), None),
        # G7: declaration, counting, roles
        ("G7 undeclared", only(T("14,521 tests"), {}, "not declared"), "G7"),
        ("G7 ok", g7_numbers(ex(lines=T("14,521 tests")), {"ledger": [{"text": "14,521", "id": "L004"}]}, L), None),
        ("G7 wrong", only(T("14,522"), {"ledger": [{"text": "14,522", "id": "L004"}]}, "does not match"), "G7"),
        ("G7 stale", only([], {"ledger": [{"text": "14,521", "id": "L004"}]}, "not printed"), "G7"),
        ("G7 count", only(T("45 Supported; 45 of 52"), {"ledger": [{"text": "45", "id": "L224", "idx": 1}, {"text": "52", "id": "L224", "idx": 2}]}, "printed 2"), "G7"),
        ("G7 count ok", g7_numbers(ex(lines=T("45 Supported; 45 of 52")), {"ledger": [{"text": "45", "id": "L014", "idx": "Supported", "anchor": "{} Supported"}, {"text": "45", "id": "L224", "idx": 1, "anchor": "{} of 52"}, {"text": "52", "id": "L224", "idx": 2}]}, L), None),
        ("G7 anchor", only(T("45 of 52"), {"ledger": [{"text": "45", "id": "L224", "idx": 1, "anchor": "{} Supported"}, {"text": "52", "id": "L224", "idx": 2}]}, "anchor"), "G7"),
        ("G7 anchor exact", only(T("7") + T("Distinct 7"), {"ledger": [{"text": "7", "id": "L082", "idx": 3, "anchor": "{}"},
                                                                    {"text": "7", "id": "L211", "idx": 2, "anchor": "Distinct {}"}]}, "anchor"), "G7"),
        ("G7 anchor needed", only(T("45 Supported; 45 of 52"), {"ledger": [{"text": "45", "id": "L014", "idx": "Supported"}, {"text": "45", "id": "L224", "idx": 1}, {"text": "52", "id": "L224", "idx": 2}]}, "needs an 'anchor'"), "G7"),
        ("G7 pct", only(T("5"), {"ledger": [{"text": "5%", "id": "L239", "pct": True}]}, "not"), "G7"),
        ("G7 pct ok", g7_numbers(ex(lines=T("at 5% FDR")), {"ledger": [{"text": "5%", "id": "L239", "pct": True}]}, L), None),
        ("G7 role", only(T("52 of 45"), {"ledger": [{"text": "52", "id": "L224", "idx": 1}, {"text": "45", "id": "L224", "idx": 2}]}, "does not match"), "G7"),
        ("G7 dup", only(T("1,948 and 1,948"), {"ledger": [{"text": "1,948", "id": "L999"}, {"text": "1,948", "id": "L001"}]}, "unknown ledger id"), "G7"),
        ("G7 needidx", only(T("45"), {"ledger": [{"text": "45", "id": "L224"}]}, "name the element"), "G7"),
        ("G7 display digits", only(T("48"), {"ledger": [{"text": "48", "id": "L032"}]}, "does not match"), "G7"),
        ("G7 int idx on dict", only(T("45"), {"ledger": [{"text": "45", "id": "L014", "idx": 1}]}, "use the dict key"), "G7"),
        ("G7 sign", only(T("+0.191"), {"ledger": [{"text": "+0.191", "id": "L050", "idx": 1}]}, "does not match"), "G7"),
        ("G7 sign ok", g7_numbers(ex(lines=T("−0.191")), {"ledger": [{"text": "−0.191", "id": "L050", "idx": 1}]}, L), None),
        ("G7 sci", g7_numbers(ex(lines=[line(sp("p = 4.7 × 10", 7.0, bbox=(10, 10, 60, 20)), sp("−9", 5.6, bbox=(60.4, 8, 66, 17.9)))]),
                              {"ledger": [{"text": "4.7 × 10⁻⁹", "id": "L050", "idx": 2}]}, L), None),
        ("G7 sci wrong", only([line(sp("p = 4.6 × 10", 7.0, bbox=(10, 10, 60, 20)), sp("−14", 5.6, bbox=(60.4, 8, 70, 17.9)))],
                              {"ledger": [{"text": "4.6 × 10⁻¹⁴", "id": "L261", "idx": 1}]}, "does not match"), "G7"),
        ("G7 power", g7_numbers(ex(lines=T("p < 10^−6")), {"ledger": [{"text": "10⁻⁶", "id": "L230"}]}, L), None),
        ("G7 power sign", only(T("p < 10^6"), {"ledger": [{"text": "10⁶", "id": "L230"}]}, "does not match"), "G7"),
        ("G7 small same line", only([small_same], {}, "'12'"), "G7"),
        ("G7 subscript ignored", g7_numbers(ex(lines=[sub]), {}, L), None),
        ("G7 extra digit", only(T("14,5210 tests"), {"ledger": [{"text": "14,521", "id": "L004"}]}, "cannot read"), "G7"),
        ("G7 derived ok", g7_numbers(ex(lines=T("69%")), {"derived": [{"text": "69%", "num": ["L225", 1], "den": ["L225", 2], "pct": True}]}, L), None),
        ("G7 derived bad", only(T("70%"), {"derived": [{"text": "70%", "num": ["L225", 1], "den": ["L225", 2], "pct": True}]}, "derived"), "G7"),
        # literals
        ("G7 noreason", only(T("0.5"), {"literals": [{"text": "0.5"}]}, "needs why"), "G7"),
        ("G7 kind no reason", only(T("0.5"), {"literals": [{"text": "0.5", "why": "axis_tick", "axis": "x"}]}, "needs why"), "G7"),
        ("G7 tick ok", g7_numbers(ex(lines=T("0 0.5 1")), {"literals": [tick("0"), tick("0.5"), tick("1")]}, L), None),
        ("G7 tick uneven", only(T("0 0.5 14,522"), {"literals": [tick("0"), tick("0.5"), tick("14,522")]}, "evenly spaced"), "G7"),
        ("G7 tick step", only(T("0 0.37 0.74"), {"literals": [tick("0"), tick("0.37"), tick("0.74")]}, "not a round number"), "G7"),
        ("G7 tick coords ok", g7_numbers(ex(lines=T("131.570 131.575 131.580")), {"literals": [tick("131.570"), tick("131.575"), tick("131.580")]}, L), None),
        ("G7 two ticks ok", g7_numbers(ex(lines=T("0 5")), {"literals": [tick("0"), tick("5")]}, L), None),
        ("G7 lone tick", only(T("0.5"), {"literals": [tick("0.5")]}, "only 1 tick"), "G7"),
        ("G7 shadow", only(T("0.87 0.84"), {"ledger": [{"text": "0.84", "id": "L223", "idx": [1, 3]}], "literals": [{"text": "0.87", "why": "other: label"}]}, "equals a value of"), "G7"),
        ("G7 shadow int", only(T("52"), {"literals": [{"text": "52", "why": "other: label"}]}, "equals a value of"), "G7"),
        ("G7 shadow ack", g7_numbers(ex(lines=T("0.87")), {"literals": [{"text": "0.87", "why": "other: label", "not_ledger": "test"}]}, L), None),
        ("G7 schema", only(T("3"), {"literals": [{"text": "3", "why": "other: x", "count": "two"}]}, "count must be"), "G7"),
        # format and minus
        ("G7m", g7(T("−0.5 -0.25 0"), {"literals": [tick("−0.5"), tick("−0.25"), tick("0")]}, "G7m"), "G7m"),
        ("G7m exponent", g7(T("p = 4.7 × 10^-9"), {"ledger": [{"text": "4.7 × 10⁻⁹", "id": "L050", "idx": 2}]}, "G7m"), "G7m"),
        ("G7m after comma", g7(T("(0.1,-0.2)"), {"literals": [{"text": "0.1", "why": "other: t", "not_ledger": "t"}, {"text": "-0.2", "why": "other: t", "not_ledger": "t"}]}, "G7m"), "G7m"),
        ("G7m range ok", g7(T("0.68–1.50 and 0.68-1.50"), {"literals": [{"text": "0.68", "why": "other: t", "count": 2, "not_ledger": "t"}, {"text": "1.50", "why": "other: t", "count": 2, "not_ledger": "t"}]}, "G7m"), None),
        ("G7f e-notation", g7(T("1e-06"), {}, "G7f"), "G7f"),
        ("G7f NA", g7(T("p = NA"), {}, "G7f"), "G7f"),
        ("G7f leading dot", g7(T("p < .05"), {}, "G7f"), "G7f"),
        ("G7f en-dash minus", g7(T("beta –0.19"), {}, "G7f"), "G7f"),
    ]
    cases += af3_cases()   # AF3 pipeline (2026-10-16): G9 af3 legend, G12 sources, G13 scope, G14 pointers
    assert numbers("n = 1,948, p") == ["1,948"], numbers("n = 1,948, p")
    assert numbers("IL-6 and Type 2 diabetes at chr6:131.5 (SLC16A1, 1000 Genomes; Figure 2b; rs2781668; 12q24)") == [], \
        numbers("IL-6 and Type 2 diabetes at chr6:131.5 (SLC16A1, 1000 Genomes; Figure 2b; rs2781668; 12q24)")
    assert numbers("OR 1.2 [1, 3]") == ["1.2", "1", "3"], numbers("OR 1.2 [1, 3]")
    assert [x for x in g7_numbers(ex(lines=T("CAD-0.19")), {"literals": [{"text": "0.19", "why": "other: t", "not_ledger": "t"}]}, L) if x.startswith("G7m")], "a hyphen glued to a label before a decimal is a minus sign"
    assert norm("4.7 × 10⁻⁹") == norm("4.7×10^−9") == "4.7×10^-9" and norm("10⁶") == "10^6" != norm("106")
    assert PG.parse_text_num(to_gate_form("10^−6"))[0] == 1e-6 and PG.parse_text_num(to_gate_form("10^6"))[0] == 1e6
    assert legend_panels("**Figure 2. Map.** **a**, the map. **b**, forest. (c) inset") == ["a", "b", "c"]
    # liveness on a REAL cairo PDF: the 3a gate adversary's fixture (9-pt 'Tests' then 7-pt '12', 0.4 pt apart)
    fx = os.path.join(HERE, "fixtures", "smallint_after.pdf")
    if os.path.exists(fx):
        got = [t for l in extract(fx)["lines"] for t in numbers(l["text"])]
        assert "12" in got, got
    bad = 0
    for lab, out, want in cases:
        ok = (not out) if want is None else (bool(out) and all(x.startswith(want) for x in out))
        print("  %-20s %s %s" % (lab, "ok" if ok else "MISS", out[:1]))
        bad += not ok
    print("SELFTEST %s (%d cases + real-PDF liveness %s)" % ("PASS" if not bad else "FAIL", len(cases), "on" if os.path.exists(fx) else "MISSING"))
    return 1 if bad or not os.path.exists(fx) else 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--selftest" in a:
        sys.exit(selftest())
    if "--pointers" in a or "--active" in a:
        sys.exit(cli_global("--pointers" in a))
    sys.exit(run([x for x in a if x.startswith("Figure")] or None))
