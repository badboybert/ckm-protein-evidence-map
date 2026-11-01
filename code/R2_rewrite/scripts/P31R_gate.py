#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 rewrite -- source gate for the section markdown (Stage 2).

Checks, per section file in ../sections/P31R_*.md:
  N1  every numeral in the prose carries a ledger tag [[Lxxx]] (or sits in an exempt context: figure/table/file labels,
      reference brackets, gene/variant/build names);
  N2  every tagged numeral matches a number in that ledger row (value or display), to the printed precision;
  N3  every tag names a ledger row that exists;
  V1  no banned vocabulary (writing frame section 3);
  W1  word count within +/-10% of the section budget; W2 share of sentences over 35 words (reported; flagged > 20%).
Run:  python P31R_gate.py [--selftest]      exit 1 if any N1/N2/N3/V1/W1 failure.
R2 mode (added 2026-10-14; the default run above is unchanged and reads ../sections as before):
      python P31R_gate.py --sections DIR [--files A.md,B.md]          the same checks on another directory / file list
      python P31R_gate.py --r2                                         = --sections ../sections_R2, the three R2 prose drafts
                                                                         (Methods, Results, Discussion) and --exemptions R2_numeral_exemptions.tsv
      python P31R_gate.py --sections DIR --exemptions FILE            extra non-quantity exemption patterns (class, regex, justification, example)
  The checks N1-N4, S1, S2, V1, W1, W2 are the same functions with the same rules; only the directory, the file list and, when a
  file is given, extra EXEMPT patterns differ. R2 mode also prints a NOTE (never a failure) for every unresolved {{id}} placeholder.
  The table / supplement R2 drafts (Table 3, SM Methods, SM Results) are run by P31R_tablegate.py --r2.
  --selftest also exercises R2 mode on a copy of the R2 Results draft (a mutation there must trip N2, N1 and N3).
D1 (added 2026-10-14): for every ledger row named by a tag in a checked file, the row's `display` must equal its stored `value` rounded at
  the display's own precision (see the D1 header below). Runs in every mode, after N1-N4/S/V/W, as "FAIL D1 Lxxx [MISMATCH|UNDERIVED|DATE] ...".
      python P31R_gate.py --d1-all                                     D1 over EVERY ledger row (no prose files); exit 1 if any row is inconsistent
      python P31R_gate.py ... --ledger FILE.csv                        read another ledger (the self-test uses a COPY with one value changed)
Stage 5 step 1a (2026-10-16, S2-D27) -- D1 changes, all with a self-test and a mutation that trips the check it names (selftest_s5):
      POSITION    a display number must match the stored leaf in ITS OWN POSITION: leaves of a list value are read in stored order and the
                  display's quantities must follow them in order. Step 1d (S2-D28) made the pointer STRICT: it moves past a matched leaf, so a
                  leaf is read by one display quantity and a changed leaf cannot hide behind the next equal one (L430: HF 6 -> 7 was silent).
                  A row whose display reads one leaf twice (p < 0.05 / p >= 0.05), or reads the leaves in another order ("9 of 10" stored as
                  [10, 9]), declares the leaf each display quantity reads in R2_d1_display_order.tsv; the declared leaf is then the only one
                  that may match. Values holding a JSON object are not positional (key order is not display order) and keep the any-leaf rule.
      sci floor   a scientific-notation token compares at a RELATIVE floor (1e-9 x |value|), not the absolute 1e-15 (7.2e-93 matched 9.8e-93).
      uncited     D1_UNCITED: rows exempt from D1 only while no R2 surface cites them; a citation (a tag in a checked file, or the bare id in
                  an operative surface) fails as [EXEMPT-CITED]. D1_UNCITED_PENDING is a declared, temporary allowance (empty since step 1c).
                  Step 1d: L409 (the row-specific label exemption, now removed) is exempt on this list instead; its citations moved to L962.
      four-digit  a [[L1000]] tag fails loudly as N3 (the id space is L000-L999); scan_pairs() is the one mutation-target scanner.
      RESERVED_TEST_IDS (L971-L984, L999) are the ids the self-tests use: no real row may take them.
      python P31R_gate.py --d1-propose-order                          print display-order declarations for rows that fail only on POSITION
      python P31R_gate.py ... --d1-order FILE --surfaces A,B           another order file / another citation-surface list (self-tests)
"""
import csv, glob, io, json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SEC = os.path.join(ROOT, "sections")
LEDGER = os.path.join(os.path.dirname(ROOT), "R1_stage1", "S1R1_06_claims_ledger", "S1R1_CLAIMS_LEDGER.csv")
BUDGET = {"P31R_Methods.md": 2100, "P31R_Results.md": 2400, "P31R_Background.md": 550, "P31R_Discussion.md": 1450, "P31R_Conclusions.md": 100,
          "P31R_Abstract.md": 350, "P31R_ResearchInsights.md": 170}
BUDGET.update({"P31R_Legend_Fig%d.md" % i: 250 for i in range(1, 6)})
HARD = {"P31R_Legend_Fig%d.md" % i for i in range(1, 6)}   # figure legends: <= 250 words is a journal-frame ceiling, no +10%
SKIP = {"P31R_TitlePage.md"}   # non-claim front matter (ORCID, postcode, phone)
TAG = re.compile(r"\s?\[\[(L\d{3}[a-z]?(?:[A-Za-z0-9]*)?)\]\]")
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
BANNED = [r"calibrat\w*", r"architecture", r"\baxis\b", r"\bcomponents?\b", r"\bgates?\b", r"sidecar", r"\bnative\b", r"\bcensus\b",
          r"\banchor\w*", r"\bheld\b", r"\bfrozen\b", r"\bregistered\b", r"mirage", r"in this generation", r"by construction",
          r"not chosen post hoc", r"re-derived", r"\bcache\b", r"\bedges?\b", r"\bbackbone\b", r"\bendpoints?\b", r"\bscrub\w*",
          r"false[- ]negative", r"\bspecificity\b", r"two-sided", r"\bdecisive\b", r"sign-reversed", r"ancestry-divergent", r"rule variants"]
REPLICAT_OK = [r"not (?:an? )?independent replication", r"rather than (?:an? )?independent replication", r"not replication"]
EXEMPT = [r"\b\d{3}-\d{4}-[A-Z]-\d{3}-\d{3}-MY\d\b",   # Stage 10 step 4 (2026-10-31): NSTC grant identifiers in the Funding paragraph (a label, not a quantity); FIRST, because the masking runs in list order and the upper-case token pattern below would take "B-182-011-MY3" first
          r"\{SM:[a-z0-9_-]+\}", r"\[@[^\]]+\]", r"\b(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\d+(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter)\b", r"(?:Figures?|Fig\.)\s+\d+[a-z]?(?:[–-]\d+[a-z]?)?", r"Tables?\s+S?\d+", r"Additional file\s+\d", r"\bSM\d+(?:[–-]SM?\d+)?",
          r"\[\d+(?:[,–-]\s?\d+)*\]", r"1000 Genomes", r"FinnGen R10", r"GRCh3[78]", r"\brs\d+", r"\b\d+[pq]\d+(?:\.\d+)?\b",
          r"\b[A-Z][A-Z0-9-]*\d[A-Z0-9-]*\b", r"\bp1(?:2)?\b", r"\bPP\.H[34]\b", r"\bp2\b", r"\bH[34]\b", r"Glu351Gln", r"Gly388Arg",
          r"[αβγδκ]\d", r"[Tt]ype [12] [Dd]iabetes", r"(?:DIAMANTE|CARDIoGRAMplusC4D|GIGASTROKE|MEGASTROKE)\s(?:19|20)\d\d", r"PheCode\s[\d.]+",
          r"95%\s(?:interval|confidence interval|CI)", r"Olink Explore 3072", r"coloc 5\.2\.3", r"\bL\s?=\s?10\b", r"ln\(eGFR\)", r"\b(?:two|one)-sample\b", r"SomaScan v\d", r"\bR\s\d\.\d\.\d", r"CC BY 4\.0", r"\bdoi\.org/10\.\d{4,9}/[A-Za-z0-9._-]+"]   # Stage 10 step 3a: the licence name in the Software paragraph (S2-D66 C2); G5 (2026-11-01): the archive DOI URL in the Declarations (a locator, not a quantity); APPENDED, never inserted: _PROSE_ONLY indexes EXEMPT
# Figure subset of EXEMPT, shared with figures/scripts/P31R_figgate.py (one list, no fork): drop the prose-only patterns
# ({SM:key}, [@key], SMn, reference brackets, which on a figure would hide real numbers such as "[1, 3]"), and use an
# upper-case token pattern WITHOUT hyphens, so "CAD-0.19" is not masked as a gene name; IL-6 and chr coordinates are added.
_PROSE_ONLY = {EXEMPT[1], EXEMPT[2], EXEMPT[7], EXEMPT[8], EXEMPT[14]}   # Stage 10 step 4 (2026-10-31): indices shifted by one for the grant pattern at index 0
EXEMPT_FIG = [x for x in EXEMPT if x not in _PROSE_ONLY] + [r"\b[A-Z][A-Z0-9]*\d[A-Z0-9]*\b", r"\bIL-\d+[A-Z]?\b",
                                                               r"\bchr\d+[:\d.,–\-]*", r"[Tt]ype [12]\b",
                                                               r"\bCARDIoGRAMplusC4D\b", r"\bChromosome\s\d+\b",
                                                               r"Figures?\s+\d+[a-z]?(?:(?:,\s?|\s(?:and|&)\s|[–-])\d*[a-z]?)*"]
EXTRA_EXEMPT = []   # R2 mode only: extra non-quantity patterns from R2_numeral_exemptions.tsv (load_exemptions); empty in the default run
R2_SEC = os.path.join(ROOT, "sections_R2")
R2_PROSE = ["P31R_Methods_R2_draft.md", "P31R_Results_R2_draft.md", "P31R_Discussion_R2_draft.md"]
R2_EXEMPT_FILE = os.path.join(HERE, "R2_numeral_exemptions.tsv")
# An exemption must carry its label words. A pattern that matches any of these bare quantities is refused (it would hide a number).
EXEMPT_PROBES = ["0.10", "0.05", "0.90", "0.80", "503", "489", "5,000", "14", "52", "p 0.364", "3.21 × 10⁻⁶", "90%", "97%", "2026", "2023", "1.96",
                 "0.993", "−31% to 54%", "14 of the 52", "the 503 panel", "489 panel", "p12 = 10⁻⁵", "20,000 draws", "OR 2.74", "of 81"]
NUM = re.compile(r"(?<![A-Za-z0-9_.\-])([−-]?\d[\d,]*(?:\.\d+)?)(\s?×\s?10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)?|(?<![A-Za-z0-9_.])(10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)")


def ledger():
    L = {}
    with io.open(LEDGER, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            L[r["id"]] = r
    return L


def nums_in(row):
    out = []
    txt = row["value"] + " " + row["display"]
    txt = txt.replace("−", "-").replace("^", "")
    for m in re.finditer(r"([-+]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?)(?:\s?×\s?10(-?\d+))?", txt):
        try:
            v = float(m.group(1).replace(",", ""))
        except ValueError:
            continue
        if m.group(2):
            v = v * 10 ** int(m.group(2))
        out.append(v)
    # bare 10^-k written as "10-6" after stripping ^
    for m in re.finditer(r"10(-\d+)", txt):
        out.append(10 ** int(m.group(1)))
    return out


def parse_text_num(s):
    s = s.replace("−", "-").translate(SUP)
    m = re.match(r"\s*([-]?\d[\d,]*(?:\.\d+)?)\s?(?:×\s?10(-?\d+))?\s*$", s)
    if m:
        base = m.group(1).replace(",", "")
        v = float(base)
        dec = len(base.split(".")[1]) if "." in base else 0
        if m.group(2):
            e = int(m.group(2))
            return v * 10 ** e, 0.5 * 10 ** (-dec) * 10 ** e, True
        return v, 0.5 * 10 ** (-dec), False
    m = re.match(r"\s*10(-?\d+)\s*$", s)
    if m:
        v = 10 ** int(m.group(1))
        return v, v * 1e-9, True
    raise ValueError(s)


def matches(tv, tol, sci, cands, pct):
    """True if some candidate rounds to the printed value. A candidate exactly on a rounding boundary counts only if
    it is itself no more precise than the print (it may be a rounded copy: double rounding is how +0.0395 printed
    as +0.040 when the full value was +0.03945). Sign must match: there is no sign-free fallback."""
    for c in cands:
        for cc in ([c, c * 100] if pct else [c]):
            d = abs(tv - cc)
            if sci:
                if cc != 0 and d <= max(tol, abs(cc) * 1e-9) * 1.0001:
                    return True
            elif d < tol * (1 - 1e-9) + 1e-15:
                return True
            elif abs(d - tol) <= tol * 1e-6 + 1e-15:
                return "BOUNDARY"
    return False


# ---------------------------------------------------------------------------------------------------------------------
# D1 (added 2026-10-14): the ledger's DISPLAY string must be the STORED VALUE rounded at the display's own precision.
# N2 compares prose to a row using value AND display together, so a row whose display drifted from its value (or whose value was
# changed after the display was written) still passes N2. D1 closes that: every number printed in `display` must be recoverable from
# `value` alone. Rules (one implementation; called per checked file by main(), by P31R_tablegate.py, and by --d1-all):
#   * value is parsed as JSON (scalar / list / dict) else as raw text; numeric leaves are collected; number-bearing string leaves
#     (e.g. "24.3-29.5", "PMID 39870808") are tokenised by the SAME tokenizer as the display;
#   * each display number must round-trip to some leaf at the display's printed precision (decimals, thousands commas, signs,
#     x 10^k, 1.2e-6, 10^-6, "%" = leaf x 100, "a-b%" = both ends percent, ranges "0.12-0.32", words one..nine);
#   * a signed token (+/-) must match the leaf's sign; an unsigned token matches the leaf's magnitude;
#   * ISO dates (and "D Month YYYY" / "D and D Month YYYY") are exact strings: the ISO date must occur verbatim in the stored value;
#   * labels (EXEMPT patterns, D1_LABELS, letter-adjacent numerals such as T2D, PP.H3, p12, rs123, "arginase-1") are not quantities; a bare year after a
#     capitalised word is a label ONLY when it matches no leaf;
#   * a failing token is classed MISMATCH (a stored number within one unit of the printed last place exists: the display drifted from
#     the value) or UNDERIVED (nothing near it is stored: the display carries a number its value cannot produce). Both fail D1.
# ---------------------------------------------------------------------------------------------------------------------
MONTHS = {m: i for i, m in enumerate("January February March April May June July August September October November December".split(), 1)}
D1_WORDS = {w: i for i, w in enumerate("one two three four five six seven eight nine".split(), 1)}
D1_TOK = re.compile(
    r"(?P<iso>(?<![\dA-Za-z])\d{4}-\d{2}-\d{2}(?![\dA-Za-z]))"
    r"|(?P<tdate>(?<![A-Za-z0-9_.])(?P<days>\d{1,2}(?:(?:,\s?|\s+and\s+)\d{1,2})*)\s+(?P<mon>" + "|".join(MONTHS) + r")\s+(?P<yr>\d{4})(?!\d))"
    r"|(?:(?<![A-Za-z0-9_.])(?<![A-Za-z]-)(?P<ssign>[-+])?(?P<sm>\d+(?:\.\d+)?)\s?×\s?10\^?(?P<se>-?\d+))"
    r"|(?:(?<![A-Za-z0-9_.])(?<![A-Za-z]-)(?P<esign>[-+])?(?P<em>\d+(?:\.\d+)?)[eE](?P<ee>[-+]?\d+)(?![\dA-Za-z]))"
    r"|(?:(?<![A-Za-z0-9_.])10\^(?P<be>-?\d+))"
    r"|(?:(?<![A-Za-z0-9_.])(?<![A-Za-z]-)(?P<psign>[-+])?(?P<pi>\d{1,3}(?:,\d{3})+|\d+)(?P<pf>\.\d+)?)"
    r"|(?:\b(?P<w>" + "|".join(D1_WORDS) + r")\b(?!-))", re.I)


# Label-only numerals that carry their label words (D1 only; the prose gate's EXEMPT list is unchanged): "1000G" (1000 Genomes), "6:131573339" (chromosome:position).
D1_LABELS = [r"\b1000G\b", r"\b(?:\d{1,2}|X|Y):(?=\d{4,})"]
# Rows exempt from D1 as "not cited by R2" (S2-D27): thresholds or denominators in the display that the value does not store. The exemption holds only
# while no operative R2 surface cites the row (a tag in a checked file, or its bare id in a surface of r2_surface_files()); a citation fails as
# [EXEMPT-CITED] (so a listed row that is cited is NOT exempt: its full D1 follows, and 1c's successor retags make it uncited again; the list is
# therefore also the home of the rows that await a successor, and nothing needs deleting when a citation is found). D1_UNCITED_PENDING is a declared, temporary allowance
# (a file that may still name an exempt row until a planned repoint); an entry fails as [EXEMPT-STALE] once the mention is gone. Step 1c (2026-10-18) repointed the
# manifest's L210 mention to L963, so the allowance is empty; L256, L261, L262 and L454 joined the list after their surfaces moved to L959, L960, L961 and L500.
# L409 (step 1d, S2-D28): its display says "all four internal criteria" but its value stores no 4; the successor L962 stores n_core_criteria = 4 and carries every
# citation, so L409 is exempt here (and the one-row label mechanism that used to pass it is gone).
D1_UNCITED = ["L032", "L094", "L403", "L415", "L438", "L439", "L446", "L210", "L409", "L256", "L261", "L262", "L454"]
D1_UNCITED_PENDING = {}
# Ids the self-tests use (synthetic rows L971-L984, the unknown id L999). No real row may take them; the ledger builder asserts this too
# (S1R1_id_alloc.RESERVED_TEST_IDS, which must equal this set).
RESERVED_TEST_IDS = frozenset(["L%d" % n for n in range(971, 985)] + ["L999"])
D1_ORDER_FILE = os.path.join(HERE, "R2_d1_display_order.tsv")
SURFACES_OVERRIDE = []   # --surfaces: replaces r2_surface_files() (self-tests only)


class _Tok:
    __slots__ = ("kind", "text", "tv", "tol", "sci", "bare", "signed", "pct", "iso", "start")


def _d1_tokens(s, mask=True):
    """Quantity tokens of a display / leaf string (see the D1 header). Returns a list of _Tok."""
    s = s.replace("−", "-")
    if mask:
        for pat in EXEMPT + EXTRA_EXEMPT + D1_LABELS:
            s = re.sub(pat, lambda m: "§" * len(m.group(0)), s)
    s = re.sub(r"10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)", lambda m: "10^" + m.group(1).translate(SUP), s)
    out = []
    for m in D1_TOK.finditer(s):
        t = _Tok()
        t.text, t.start, t.pct, t.sci, t.bare, t.signed, t.iso = m.group(0), m.start(), False, False, False, False, None
        if m.group("iso"):
            t.kind, t.iso = "date", [m.group("iso")]
        elif m.group("tdate"):
            t.kind = "date"
            mon = MONTHS[[k for k in MONTHS if k.lower() == m.group("mon").lower()][0]]
            t.iso = ["%s-%02d-%02d" % (m.group("yr"), mon, int(d)) for d in re.findall(r"\d+", m.group("days"))]
        elif m.group("sm") is not None:
            t.kind, t.sci, t.signed = "num", True, bool(m.group("ssign"))
            dec = len(m.group("sm").split(".")[1]) if "." in m.group("sm") else 0
            e = int(m.group("se"))
            t.tv = float(m.group("sm")) * 10 ** e * (-1 if m.group("ssign") == "-" else 1)
            t.tol = 0.5 * 10 ** (-dec) * 10 ** e
        elif m.group("em") is not None:
            t.kind, t.sci, t.signed = "num", True, bool(m.group("esign"))
            dec = len(m.group("em").split(".")[1]) if "." in m.group("em") else 0
            e = int(m.group("ee"))
            t.tv = float(m.group("em")) * 10 ** e * (-1 if m.group("esign") == "-" else 1)
            t.tol = 0.5 * 10 ** (-dec) * 10 ** e
        elif m.group("be") is not None:
            t.kind, t.sci, t.bare = "num", True, True
            t.tv = 10.0 ** int(m.group("be"))
            t.tol = t.tv * 1e-6
        elif m.group("pi") is not None:
            t.kind, t.signed = "num", bool(m.group("psign"))
            dec = len(m.group("pf")) - 1 if m.group("pf") else 0
            t.tv = float(m.group("pi").replace(",", "") + (m.group("pf") or "")) * (-1 if m.group("psign") == "-" else 1)
            t.tol = 0.5 * 10 ** (-dec)
        else:
            t.kind, t.tv, t.tol = "num", float(D1_WORDS[m.group("w").lower()]), 0.5
        if t.kind == "num":
            tail = s[m.end():m.end() + 24]
            t.pct = bool(re.match(r"(?:\s?[-–]\s?[-+]?[\d,]*\.?\d+)?\s?%", tail))
        out.append(t)
    return out


def _d1_leaf_list(value):
    """(numeric leaves in stored order, ordered). ordered is False when a JSON object is anywhere in the value: its key order says nothing about the
    order a display reads it in, so such a row keeps the any-leaf rule."""
    try:
        obj = json.loads(value)
    except ValueError:
        obj = value
    nums = []
    has_obj = []

    def walk(o):
        if o is None or isinstance(o, bool):
            return
        if isinstance(o, (int, float)):
            if math.isfinite(o):
                nums.append(float(o))
        elif isinstance(o, str):
            for t in _d1_tokens(re.sub(r"\d{4}-\d{2}-\d{2}", " ", o), mask=False):
                if t.kind == "num":
                    nums.append(t.tv)
        elif isinstance(o, (list, tuple)):
            for x in o:
                walk(x)
        elif isinstance(o, dict):
            has_obj.append(1)
            for x in o.values():
                walk(x)
    walk(obj)
    return nums, not has_obj


def _d1_leaves(value):
    """Numeric leaves of a stored value (stored order)."""
    return _d1_leaf_list(value)[0]


def _d1_match1(t, c, widen=1.0):
    """The leaf value (as compared) if token t rounds from leaf c, else None."""
    for cc in ([c, c * 100.0] if t.pct else [c]):
        if t.signed:
            d = abs(t.tv - cc)
        else:
            d = abs(abs(t.tv) - abs(cc))
        # a scientific-notation token compares at a RELATIVE floor: the absolute 1e-15 let any two values below ~1e-15 match (7.2e-93 vs 9.8e-93, S2-D28)
        floor = 1e-9 * max(abs(cc), abs(t.tv)) if t.sci else 1e-15
        if d <= t.tol * widen * (1 + 1e-9) + floor:
            return cc
    return None


def _d1_hit_at(t, leaves, start=0, widen=1.0):
    """(index, compared value) of the first leaf at or after `start` that token t rounds from, else None."""
    for i in range(start, len(leaves)):
        cc = _d1_match1(t, leaves[i], widen)
        if cc is not None:
            return i, cc
    return None


def _d1_hit(t, leaves, widen=1.0):
    """The leaf (as compared) that this token rounds from, else None. widen=2 gives the 'within one unit of the last printed place' test."""
    h = _d1_hit_at(t, leaves, 0, widen)
    return None if h is None else h[1]


_ORDER_CACHE = {}


def d1_order(path=None):
    """Declared display orders {row id: [stored-leaf index read by each display quantity, in reading order]} from R2_d1_display_order.tsv
    (tab-separated: id, order, justification). A missing file is an error, never an empty table."""
    path = path or D1_ORDER_FILE
    st = os.stat(path)
    key = (os.path.abspath(path), st.st_mtime_ns, st.st_size)
    if key not in _ORDER_CACHE:
        out = {}
        with io.open(path, encoding="utf-8", newline="") as fh:
            rd = csv.DictReader((ln for ln in fh if ln.strip() and not ln.startswith("#")), delimiter="\t", quoting=csv.QUOTE_NONE)
            assert rd.fieldnames and rd.fieldnames[:3] == ["id", "order", "justification"], "bad header in %s: %s" % (path, rd.fieldnames)
            for r in rd:
                assert re.fullmatch(r"L\d{3}", r["id"]), "bad id %r in %s" % (r["id"], path)
                assert r["id"] not in out, "duplicate id %s in %s" % (r["id"], path)
                assert re.fullmatch(r"\d+(?: \d+)*", r["order"]), "bad order for %s in %s: %r" % (r["id"], path, r["order"])
                out[r["id"]] = [int(x) for x in r["order"].split()]
        _ORDER_CACHE[key] = out
    return _ORDER_CACHE[key]


# Step 1d (S2-D28), second pass: the files that ship into the manuscript or Additional files through the builders are surfaces too. They are enumerated
# by the same file-name patterns as P31R_tablegate.py's default run (tables/P31R_Table[S]N.md, supplement/P31R_SM_<key>.md and P31R_SNote.md), which
# excludes the *_preS5S1d.* backups by name, plus sections/P31R_Declarations.md (read by P31R_build.py). selftest_s5 enumerates the directories a second,
# different way and asserts the two lists agree.
SURFACE_TABLE_RX = re.compile(r"P31R_TableS?\d+\.md")
SURFACE_SUPP_RX = re.compile(r"P31R_(?:SM_[a-z0-9_-]+|SNote)\.md")


def r2_surface_files():
    """The operative R2 surfaces that can cite a ledger row: the R2 prose / Table 3 / supplement / legend / caveat-map files, the AF3 legends, the figure
    generators, the figure registries, the approved-source manifest and its builder, and the AF1 / AF2 builders; since step 1d (second pass) also the
    shipped tables (tables/), the legacy supplement files and the Supplementary Note (supplement/), and the Declarations (sections/). Every declared
    file must exist (a surface that is not opened is not scanned)."""
    j = os.path.join
    names = list(R2_PROSE) + ["P31R_Table3_R2_draft.md", "P31R_SM_MethodsDetail_R2.md", "P31R_SM_ResultsDetail_R2.md", "P31R_Legend_AF3_R2.md",
                              "P31R_Results_R2_caveat_map.tsv"] + ["P31R_Legend_Fig%d_R2.md" % i for i in range(1, 6)]
    files = [j(R2_SEC, n) for n in names]
    files += [j(ROOT, "figures", "af3", "legends", "P31R_Legend_FigS%d.md" % i) for i in range(1, 7)]
    files += [j(ROOT, "figures", "scripts", "P31R_Fig%d.R" % i) for i in range(1, 6)] + [j(ROOT, "figures", "scripts", "P31R_FigS%d.R" % i) for i in range(1, 7)]
    reg = j(ROOT, "figures", "registry")
    files += [j(reg, "APPROVED_SOURCE_MANIFEST.tsv"), j(reg, "build_approved_source_manifest.py")]
    files += [j(reg, f) for f in sorted(os.listdir(reg)) if re.fullmatch(r"Figure(?:_\d|S\d+(?:_p\d+)?)\.json", f)]
    files += [j(ROOT, "scripts", "P31R_build_af1.py"), j(ROOT, "scripts", "P31R_build_af2.py")]
    for sub, rx in (("tables", SURFACE_TABLE_RX), ("supplement", SURFACE_SUPP_RX)):
        dd = j(ROOT, sub)
        assert os.path.isdir(dd), "citation surface directory missing: %s" % dd
        files += [j(dd, f) for f in sorted(os.listdir(dd)) if rx.fullmatch(f)]
    files.append(j(ROOT, "sections", "P31R_Declarations.md"))
    # Stage 6.2 (S2-D32 repoint a): also scan every canonical file of the component source map, so the front matter the builder reads is gated.
    mp = j(ROOT, "build", "P31R_SOURCE_MAP.tsv")
    for r_ in csv.DictReader(io.open(mp, encoding="utf-8", newline=""), delimiter="\t"):
        if r_["component"] == "references":      # generated list, no source file
            continue
        for part_ in r_["canonical_file"].split(";"):
            files += [h for h in sorted(glob.glob(j(ROOT, part_.strip()))) if not re.search(r"_pre[A-Za-z0-9_]*\.md$", h)]
    files = list(dict.fromkeys(os.path.normpath(os.path.abspath(f)) for f in files))   # one spelling per file: the map's '/' and the hard-coded '\\' paths must dedupe
    missing = [f for f in files if not os.path.isfile(f)]
    assert not missing, "declared R2 citation surface(s) missing: %s" % missing
    return files


_SCAN_CACHE = {}


def uncited_sites(files=None):
    """{exempt row id: ['relative/path:line', ...]} for every citation of a D1_UNCITED id in the surfaces (files, else --surfaces, else r2_surface_files())."""
    files = files or SURFACES_OVERRIDE or r2_surface_files()
    for f in files:
        assert os.path.isfile(f), "citation surface missing: %s" % f
    key = tuple((f, os.stat(f).st_mtime_ns, os.stat(f).st_size) for f in files)
    if key not in _SCAN_CACHE:
        rx = {i: re.compile(r"(?<![A-Za-z0-9_])%s(?![A-Za-z0-9_])" % i) for i in D1_UNCITED}
        out = {i: [] for i in D1_UNCITED}
        for f in files:
            try:
                rel = os.path.relpath(f, ROOT).replace("\\", "/")
            except ValueError:   # another drive
                rel = ".."
            if rel.startswith(".."):
                rel = f.replace("\\", "/")
            for ln, line in enumerate(io.open(f, encoding="utf-8", errors="replace"), 1):
                for i, r in rx.items():
                    if r.search(line):
                        out[i].append("%s:%d" % (rel, ln))
        _SCAN_CACHE.clear()
        _SCAN_CACHE[key] = out
    return _SCAN_CACHE[key]


def exempt_status(rid, cited=None):
    """(is_cited, failure lines) for a D1_UNCITED row. cited: sites the caller already knows (a tag in a checked file); else the surfaces are scanned.
    A cited row is not exempt (full D1 follows); a stale pending allowance is reported but the uncited row stays exempt."""
    sites = list(cited or [])
    allowed = D1_UNCITED_PENDING.get(rid, [])
    stale = False
    if cited is None:
        found = uncited_sites().get(rid, [])
        sites = [s for s in found if s.rsplit(":", 1)[0] not in allowed]
        stale = bool(allowed) and not any(s.rsplit(":", 1)[0] in allowed for s in found)
    out = []
    if sites:
        out.append("D1 %s [EXEMPT-CITED] row is exempt from D1 as 'not cited by R2' but is cited in %s: give it a successor row whose display its value derives (S2-D27)"
                   % (rid, ", ".join(sites[:4]) + (" (+%d more)" % (len(sites) - 4) if len(sites) > 4 else "")))
    if stale:
        out.append("D1 %s [EXEMPT-STALE] the pending-citation allowance (%s) is no longer needed: delete %s from D1_UNCITED_PENDING" % (rid, ", ".join(allowed), rid))
    return bool(sites), out


def d1_row(row, cited=None):
    """D1 for one ledger row. Returns a list of failure strings (empty = consistent). cited: where the caller found the row cited (see exempt_status)."""
    disp, val, rid = row["display"], row["value"], row["id"]
    if not disp.strip():
        return []
    fails = []
    if rid in D1_UNCITED:
        is_cited, fl = exempt_status(rid, cited)
        fails += fl
        if not is_cited:
            return fails
    toks = _d1_tokens(disp)
    leaves, ordered = _d1_leaf_list(val)
    qn = [t for t in toks if t.kind == "num"]
    decl = d1_order().get(rid) if ordered else None
    if not ordered and rid in d1_order():
        fails.append("D1 %s [POSITION] a display order is declared for the row but its value holds a JSON object (not positional): remove the declaration" % rid)
    if decl is not None and (len(decl) != len(qn) or any(i >= len(leaves) for i in decl)):
        fails.append("D1 %s [POSITION] the declared display order (%d entries, leaf indices up to %d) does not fit the row (%d display quantities, %d stored leaves) | display=%r"
                     % (rid, len(decl), max(decl), len(qn), len(leaves), disp[:110]))
        decl = None
    ptr, k = 0, -1
    for t in toks:
        if t.kind == "date":
            miss = [d for d in t.iso if d not in val]
            if miss:
                fails.append("D1 %s [DATE] display date '%s' (= %s) is not in the stored value: %s" % (rid, t.text, ", ".join(miss), val[:100]))
            continue
        k += 1
        anywhere = _d1_hit_at(t, leaves)
        if ordered and decl is not None:
            ok = _d1_match1(t, leaves[decl[k]]) is not None
            where = "its declared leaf #%d (%.10g)" % (decl[k], leaves[decl[k]])
        elif ordered:
            h = _d1_hit_at(t, leaves, ptr)
            ok = h is not None
            if ok:
                ptr = h[0] + 1   # strictly past the matched leaf: a leaf is read by ONE display quantity (S2-D28), so a changed leaf cannot hide behind the next equal one
            where = "an unread leaf at or after #%d in stored order (declare the display order in R2_d1_display_order.tsv if the display legitimately reads one leaf twice or reads the leaves in another order)" % ptr
        else:
            ok, where = anywhere is not None, ""
        if ok:
            continue
        if anywhere is not None:
            fails.append("D1 %s [POSITION] display '%s' (quantity %d) matches stored leaf #%d but not the leaf in its own position: %s | display=%r | value=%r"
                         % (rid, t.text.strip(), k + 1, anywhere[0], where, disp[:110], val[:110]))
            continue
        if re.fullmatch(r"(?:19|20)\d\d", t.text) and re.search(r"[A-Za-z]\s$", disp[:t.start]) and re.search(r"[A-Z]", disp[max(0, t.start - 12):t.start]):
            continue   # a publication / resource year after a name, matching no stored number: a label
        near = _d1_hit(t, leaves, widen=2.0)
        if near is not None:
            fails.append("D1 %s [MISMATCH] display '%s' is not the stored value rounded at its precision (stored value near it: %.10g) | display=%r | value=%r"
                         % (rid, t.text.strip(), near, disp[:110], val[:110]))
        else:
            fails.append("D1 %s [UNDERIVED] display '%s' has no counterpart in the stored value | display=%r | value=%r"
                         % (rid, t.text.strip(), disp[:110], val[:110]))
    return fails


def d1_propose_order(row):
    """Declaration line for a row that passes the any-leaf rule but not the positional one: the leaf each display quantity reads (an unused leaf first, else
    a reused one). None if some quantity matches no leaf, or the row needs no declaration."""
    leaves, ordered = _d1_leaf_list(row["value"])
    if not ordered or not row["display"].strip():
        return None
    qn = [t for t in _d1_tokens(row["display"]) if t.kind == "num"]
    used, idx = set(), []
    for t in qn:
        cand = [i for i in range(len(leaves)) if _d1_match1(t, leaves[i]) is not None]
        if not cand:
            return None
        pick = next((i for i in cand if i not in used), cand[0])
        used.add(pick)
        idx.append(pick)
    return idx


def d1_tags(L, tags, where="a checked file"):
    """D1 over the rows named by `tags` (unknown tags are N3's business). One failure list, ledger order of first mention. A tag in a checked file IS a citation."""
    out = []
    for tg in tags:
        if tg in L:
            out += d1_row(L[tg], cited=["tag in %s" % where] if tg in D1_UNCITED else None)
    return out


def d1_file_tags(path):
    body = strip_md(io.open(path, encoding="utf-8").read())
    seen = []
    for m in TAG.finditer(body):
        if m.group(1) not in seen:
            seen.append(m.group(1))
    return seen


def d1_all(L):
    """D1 over every ledger row; returns (failure lines, rows checked, rows with a display, tokens checked)."""
    fails, ntok, nd = [], 0, 0
    for r in L.values():
        if r["display"].strip():
            nd += 1
            ntok += len(_d1_tokens(r["display"]))
        fails += d1_row(r)
    return fails, len(L), nd, ntok


def strip_md(t):
    t = re.sub(r"^#+ .*$", "", t, flags=re.M)
    return t


def sentences(t):
    t = TAG.sub("", t)
    t = re.sub(r"\s+", " ", t)
    return [s for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", t) if s.strip()]


PAIR_RX = re.compile(r"(?<![\w.\-])(\d{1,3}(?:,\d{3})*(?:\.\d+)?) (\[\[(L\d{3,4})\]\])")


def scan_pairs(text):
    """[(number, '[[Lxxx]]', 'Lxxx')] for every '<number> [[Lxxx]]' in text: the ONE mutation-target scanner (prose and table self-tests). A four-digit tag is
    refused loudly: the id space is L000-L999, and a scanner that skipped [[L1000]] would shrink the mutation pool without any error."""
    out = []
    for m in PAIR_RX.finditer(text):
        assert len(m.group(3)) == 4, "four-digit ledger tag %s in a mutation scan: the id space is L000-L999 (widening is a tested migration, not an accident)" % m.group(2)
        out.append((m.group(1), m.group(2), m.group(3)))
    return out


def check(path, L):
    name = os.path.basename(path)
    raw = io.open(path, encoding="utf-8").read()
    body = strip_md(raw)
    fails, notes = [], []
    # N3: every tag exists
    for m in TAG.finditer(body):
        if re.match(r"L\d{4}", m.group(1)):
            fails.append("N3 four-digit ledger tag %s: the id space is L000-L999 (widening is a tested migration, not an accident)" % m.group(1))
        elif m.group(1) not in L:
            fails.append("N3 unknown tag %s" % m.group(1))
    # mask exempt spans
    tag_spans = [m.span() for m in TAG.finditer(body)]

    def in_tag(a, b):
        return any(a < te and b > ts for ts, te in tag_spans)

    masked = body
    for pat in EXEMPT + EXTRA_EXEMPT:
        masked = re.sub(pat, lambda m: m.group(0) if in_tag(*m.span()) else "§" * len(m.group(0)), masked)
    # N1/N2
    for m in NUM.finditer(masked):
        tok = m.group(0)
        start, end = m.span()
        # skip numbers inside a tag
        if re.search(r"\[\[L\d*$", masked[max(0, start - 6):start]):
            continue
        tail = masked[end:end + 40]
        tm = re.match(r"(\s?(?:%|kb|Mb|bp|pp|-fold|×)?\)?\]?,?)?\s?\[\[(L\d{3}[a-z]?[A-Za-z0-9]*)\]\]", tail)
        if not tm:
            # numbers followed by another number of a list before the tag: "3 / 45 / 33 [[L014]]" -> look ahead within 60 chars without sentence end
            look = masked[end:end + 90]
            tm2 = re.match(r"([^.;\[]{0,80}?)\[\[(L\d{3}[a-z]?[A-Za-z0-9]*)\]\]", look)
            if tm2 and not re.search(r"[A-Za-z]{4,}", tm2.group(1)):
                tag = tm2.group(2)
            else:
                ctx = body[max(0, start - 30):end + 30].replace("\n", " ")
                fails.append("N1 untagged number '%s' ... %s ..." % (tok.strip(), ctx))
                continue
        else:
            tag = tm.group(2)
        if tag not in L:
            continue
        try:
            tv, tol, sci = parse_text_num(tok)
        except ValueError:
            fails.append("N2 unparsable number '%s'" % tok)
            continue
        pct = "%" in masked[end:end + 2]
        cands = nums_in(L[tag])
        mm = matches(tv, tol, sci, cands, pct)
        ctx = body[max(0, start - 30):end + 30].replace("\n", " ")
        if mm == "BOUNDARY":
            fails.append("N2b '%s' sits on a rounding boundary of %s (a rounded copy?): print from full precision ... %s ..." % (tok.strip(), tag, ctx))
        elif not mm:
            fails.append("N2 '%s' not in %s (%s) ... %s ..." % (tok.strip(), tag, L[tag]["display"][:60], ctx))
    # N4 spelled numbers (one..twenty) followed by a tag are checked against the row; untagged ones are listed as a note
    W2N = {w: i for i, w in enumerate("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen "
                                      "sixteen seventeen eighteen nineteen twenty".split())}
    for m in re.finditer(r"\b(%s)\b(\s?\[\[(L\d{3}[a-z]?[A-Za-z0-9]*)\]\])?" % "|".join(W2N), body, flags=re.I):
        if m.group(3):
            if m.group(3) in L and not any(abs(c - W2N[m.group(1).lower()]) < 1e-9 for c in nums_in(L[m.group(3)])):
                fails.append("N4 '%s' (= %d) not in %s (%s)" % (m.group(1), W2N[m.group(1).lower()], m.group(3), L[m.group(3)]["display"][:60]))
    # S1 manuscript prose: no markdown bullet lists, numbered lists or tables in a section file
    # Stage 6.2: Research Insights is a journal-mandated highlight list ("- " lines; format asserted by P31R_build.py), so only it is exempt
    ri_file = os.path.basename(path) == "P31R_ResearchInsights.md"
    for ln, line in enumerate(raw.splitlines(), 1):
        if ri_file and re.match(r"^- \S", line):
            continue
        if re.match(r"^\s*(?:[-*+]\s|\d+[.)]\s|\|)", line):
            fails.append("S1 list/table line %d: %s" % (ln, line.strip()[:70]))
    # S2 supplementary-methods numbers come from {SM:key} keys only (numbered at build in citation order)
    for m in re.finditer(r"\bSM\d+\b", TAG.sub("", body)):
        fails.append("S2 literal '%s': cite supplementary methods as {SM:key}" % m.group(0))
    # V1 banned vocabulary
    plain = TAG.sub("", body)
    for pat in BANNED:
        for m in re.finditer(pat, plain, flags=re.I):
            fails.append("V1 banned '%s' ... %s ..." % (m.group(0), plain[max(0, m.start() - 40):m.end() + 40].replace("\n", " ")))
    for m in re.finditer(r"replicat\w*", plain, flags=re.I):
        win = plain[max(0, m.start() - 40):m.end() + 5]
        if not any(re.search(p, win, flags=re.I) for p in REPLICAT_OK):
            fails.append("V1 'replicat*' outside an allowed negation ... %s ..." % plain[max(0, m.start() - 50):m.end() + 30].replace("\n", " "))
    # W1/W2
    words = len(re.findall(r"\S+", re.sub(r"\s+", " ", plain)))
    b = BUDGET.get(name)
    if b and name in HARD and words > b:
        fails.append("W1 %d words over the hard ceiling %d" % (words, b))
    elif b and words > 1.1 * b:
        fails.append("W1 %d words over budget %d (+10%%)" % (words, b))
    elif b and words < 0.9 * b:
        notes.append("  W1 note: %d words is under budget %d (a ceiling, not a target; check content is complete)" % (words, b))
    ss = sentences(body)
    long_ = [s for s in ss if len(s.split()) > 35]
    share = len(long_) / max(1, len(ss))
    notes.append("%s: %d words (budget %s); %d sentences; >35 words %d (%.0f%%); mean %.1f"
                 % (name, words, b, len(ss), len(long_), 100 * share, sum(len(s.split()) for s in ss) / max(1, len(ss))))
    if share > 0.20:
        notes.append("  W2 note: long-sentence share above 20%")
    # subsection word counts
    for sec in re.split(r"^### ", raw, flags=re.M)[1:]:
        head, _, txt = sec.partition("\n")
        notes.append("  %-70s %4d words" % (head.strip()[:70], len(re.findall(r"\S+", TAG.sub("", txt)))))
    return fails, notes


def load_exemptions(path):
    """R2 mode: extra exemption patterns as DATA (tab-separated: class, regex, justification, example). Refused when a regex does not compile,
    can match the empty string, matches a bare quantity in EXEMPT_PROBES, or does not match its own example."""
    out = []
    with io.open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        assert rd.fieldnames and rd.fieldnames[:4] == ["class", "regex", "justification", "example"], "bad header in %s: %s" % (path, rd.fieldnames)
        for r in rd:
            if not r["class"] or r["class"].startswith("#"):
                continue
            rx = re.compile(r["regex"])
            assert not rx.fullmatch(""), "exemption %s can match the empty string" % r["class"]
            for pr in EXEMPT_PROBES:
                assert not rx.search(pr), "exemption %s masks the bare quantity %r: an exemption must carry its label words" % (r["class"], pr)
            assert rx.search(r["example"]), "exemption %s does not match its own example %r" % (r["class"], r["example"])
            out.append(r["regex"])
    assert out, "no exemption patterns in %s" % path
    return out


def _parse(argv):
    a = {"sections": None, "files": None, "r2": False, "exemptions": None, "ledger": None, "d1_all": False, "d1_order": None, "surfaces": None, "propose": False}
    i = 0
    while i < len(argv):
        x = argv[i]
        if x in ("--sections", "--files", "--exemptions", "--ledger", "--d1-order", "--surfaces"):
            assert i + 1 < len(argv), "%s needs a value" % x
            a[x[2:].replace("-", "_")] = [f for f in argv[i + 1].split(",") if f] if x in ("--files", "--surfaces") else os.path.abspath(argv[i + 1])
            i += 2
        elif x == "--r2":
            a["r2"] = True
            i += 1
        elif x == "--d1-all":
            a["d1_all"] = True
            i += 1
        elif x == "--d1-propose-order":
            a["propose"] = True
            i += 1
        elif x == "--selftest":
            i += 1
        else:
            raise SystemExit("unknown option %s (see the docstring)" % x)
    return a


def main(argv=None):
    global LEDGER, D1_ORDER_FILE
    a = _parse(sys.argv[1:] if argv is None else argv)
    r2 = a["r2"]
    old_ledger, old_order, old_surf = LEDGER, D1_ORDER_FILE, list(SURFACES_OVERRIDE)
    if a["ledger"]:
        LEDGER = a["ledger"]
    if a["d1_order"]:
        D1_ORDER_FILE = a["d1_order"]
    if a["surfaces"]:
        SURFACES_OVERRIDE[:] = [os.path.abspath(f) for f in a["surfaces"]]
    sec = a["sections"] or (R2_SEC if r2 else SEC)
    files = a["files"] or (R2_PROSE if r2 else None)
    ex = a["exemptions"] or (R2_EXEMPT_FILE if r2 else None)
    try:   # D1 lines quote ledger displays (x, superscripts, en dashes): always print as UTF-8, in the default mode too
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    EXTRA_EXEMPT[:] = load_exemptions(ex) if ex else []
    try:
        L = ledger()
        if a["propose"]:
            import tempfile
            fd, tmp_order = tempfile.mkstemp(suffix=".tsv", prefix="d1order_")
            os.close(fd)
            io.open(tmp_order, "w", encoding="utf-8", newline="").write("id\torder\tjustification\n")
            D1_ORDER_FILE = tmp_order
            for r in L.values():
                if any(x.startswith("D1 %s [POSITION]" % r["id"]) for x in d1_row(r)) and not any(("D1 %s [MISMATCH]" % r["id"]) in x or ("D1 %s [UNDERIVED]" % r["id"]) in x for x in d1_row(r)):
                    idx = d1_propose_order(r)
                    if idx is not None:
                        print("%s\t%s\t%s" % (r["id"], " ".join(map(str, idx)), r["display"][:60].replace("\t", " ")))
            return 0
        if a["d1_all"]:
            fl, nrows, nd, ntok = d1_all(L)
            for x in fl:
                print("  FAIL " + x)
            bad_rows = sorted({x.split()[1] for x in fl})
            print("D1 over the whole ledger %s: %d rows, %d with a display, %d display quantities checked; %d inconsistent row(s), %d failing token(s) "
                  "(MISMATCH %d, UNDERIVED %d, DATE %d)" % (LEDGER, nrows, nd, ntok, len(bad_rows), len(fl), sum("[MISMATCH]" in x for x in fl),
                                                           sum("[UNDERIVED]" in x for x in fl), sum("[DATE]" in x for x in fl)))
            ordd = d1_order()
            cited_ids = sorted(i for i in D1_UNCITED if any(x.split()[1] == i and "[EXEMPT-CITED]" in x for x in fl))
            print("D1 positional: %d row(s) declare a display order (%s); exempt as not cited by R2: %s; listed but cited, so NOT exempt: %s; pending-citation allowance: %s"
                  % (len(ordd), os.path.basename(D1_ORDER_FILE), ", ".join(i for i in D1_UNCITED if i not in cited_ids) or "none", ", ".join(cited_ids) or "none", ", ".join("%s only in %s" % (k, "/".join(os.path.basename(f) for f in v)) for k, v in D1_UNCITED_PENDING.items()) or "none"))
            return 1 if fl else 0
        if files is None:
            files = sorted(f for f in os.listdir(sec) if f.startswith("P31R_") and f.endswith(".md") and f not in SKIP)
            if not r2 and not a["sections"]:
                # Stage 6.2: the default mode gates the sections/ files the component source map ships (front matter, Declarations).
                # The pre-R2 Methods/Results/Discussion/legends left in sections/ are superseded by sections_R2/ (gated by --r2) and are
                # named here, not silently skipped.
                mp = os.path.join(ROOT, "build", "P31R_SOURCE_MAP.tsv")
                shipped = set()
                for r_ in csv.DictReader(io.open(mp, encoding="utf-8", newline=""), delimiter="\t"):
                    for part_ in r_["canonical_file"].split(";"):
                        if os.path.dirname(part_.strip().replace("\\", "/")) == "sections":
                            shipped.add(os.path.basename(part_.strip()))
                assert shipped, "source map lists no sections/ file"
                superseded = [f for f in files if f not in shipped]
                files = [f for f in files if f in shipped]
                print("MODE: default = the %d sections/ file(s) of the component source map; superseded, not gated: %s" % (len(files), ", ".join(superseded) or "none"))
        assert files, "no section files"
        if r2 or a["sections"] or a["files"] or ex:
            print("MODE: sections=%s files=%d exemption patterns=%d" % (sec, len(files), len(EXTRA_EXEMPT)))
        total = 0
        for f in files:
            fails, notes = check(os.path.join(sec, f), L)
            fails = fails + d1_tags(L, d1_file_tags(os.path.join(sec, f)), where=f)
            for n in notes:
                print(n)
            for x in fails:
                print("  FAIL " + x)
            if r2:
                npl = len(re.findall(r"\{\{[A-Za-z0-9_.]+\}\}", io.open(os.path.join(sec, f), encoding="utf-8").read()))
                if npl:
                    print("  NOTE %d unresolved {{id}} placeholder(s) in %s (not a gate failure; convert them before this file counts as gated)" % (npl, f))
            total += len(fails)
        print("GATE %s: %d failure(s) over %d file(s)" % ("RED" if total else "GREEN", total, len(files)))
        return 1 if total else 0
    finally:
        EXTRA_EXEMPT[:] = []
        LEDGER, D1_ORDER_FILE = old_ledger, old_order
        SURFACES_OVERRIDE[:] = old_surf


def _cli(argv):
    """Runs main() with command-line options and returns (rc, [failure lines as printed], full stdout)."""
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = main(argv)
    out = buf.getvalue()
    return rc, [l[len("  FAIL "):] for l in out.splitlines() if l.startswith("  FAIL ")], out


def selftest_r2(L):
    """R2 mode through the real command-line path. Pointing at another directory must read THAT directory, and a mutation in a copy of an
    R2 draft must trip N2, N1 and N3. Targets are found in the draft itself ('<number> [[Lxxx]]' occurring exactly once and not already
    failing), so the test keeps working as the draft changes; it FAILS (does not skip) if no target exists. Exemptions must hide labels
    (Rule 3(ii), Stage 1) and never quantities."""
    import collections, shutil, tempfile
    bad = 0

    def rep(lab, ok, detail=""):
        nonlocal bad
        print("selftest R2 %-34s %s %s" % (lab, "OK" if ok else "MISS", detail))
        bad += not ok

    name = "P31R_Results_R2_draft.md"
    src = os.path.join(R2_SEC, name)
    if not os.path.isfile(src):
        rep("R2 draft present", False, src)
        return 1
    d = tempfile.mkdtemp()
    cp = os.path.join(d, name)
    shutil.copy(src, cp)

    def fails_of(sec):
        return _cli(["--sections", sec, "--files", name, "--exemptions", R2_EXEMPT_FILE])

    rc0, base, _ = fails_of(d)
    rc1, real, _ = fails_of(R2_SEC)
    rep("copy == real file (--sections)", base == real and rc0 == rc1, "(%d baseline failure lines)" % len(base))
    rcr, _, outr = _cli(["--r2"])
    rep("--r2 opens the three R2 prose drafts", all(("%s:" % f) in outr for f in R2_PROSE) and "MODE:" in outr and rcr in (0, 1))
    text = io.open(cp, encoding="utf-8", newline="").read()
    cands = []
    for tok, tag, lid in scan_pairs(text):
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
    rep("real R2 draft unchanged", fails_of(R2_SEC)[1] == real)
    # exemptions hide labels, never quantities
    probe = os.path.join(d, "probe.md")
    io.open(probe, "w", encoding="utf-8").write("Rule 3(ii) uses 5,000 draws in Stage 1 and Amendment 1b.\n")
    with_ex = [x for x in _cli(["--sections", d, "--files", "probe.md", "--exemptions", R2_EXEMPT_FILE])[1] if x.startswith("N1")]
    without = [x for x in _cli(["--sections", d, "--files", "probe.md"])[1] if x.startswith("N1")]
    toks = lambda fs: sorted(re.findall(r"untagged number '([^']+)'", " ".join(fs)))
    rep("exemption hides labels only", toks(with_ex) == ["5,000"] and len(toks(without)) > 1, "(with: %s; without: %s)" % (toks(with_ex), toks(without)))
    # the loader refuses a pattern that would hide a quantity, and one that does not match its own example
    badf = os.path.join(d, "bad_exemptions.tsv")
    refused = {}
    for lab, rx, ex_ in (("bare number", r"\b\d+(?:\.\d+)?\b", "0.10"), ("panel size", r"\b(?:503|489)\b", "503"), ("own example", r"\bStage\s+\d+", "Phase 1")):
        io.open(badf, "w", encoding="utf-8").write("class\tregex\tjustification\texample\n%s\t%s\tx\t%s\n" % (lab.replace(" ", "_"), rx, ex_))
        try:
            load_exemptions(badf)
            refused[lab] = False
        except AssertionError:
            refused[lab] = True
    rep("loader refuses quantity-hiding patterns", all(refused.values()), "(%s)" % refused)
    return 1 if bad else 0


def selftest_d1(L):
    """D1 through the real command-line path on LEDGER COPIES (never the real ledger). A consistent row must pass; the same row with ONLY its
    value changed (or ONLY its display changed) must trip D1 on that row's id; a date with a changed day must trip; a row that no checked file
    tags is not D1-checked in file mode but is caught by --d1-all. Synthetic rows cover each display form D1 parses. Mutation targets are
    asserted to exist (the test FAILS rather than skips)."""
    import tempfile
    bad = 0

    def rep(lab, ok, detail=""):
        nonlocal bad
        print("selftest D1 %-44s %s %s" % (lab, "OK" if ok else "MISS", detail))
        bad += not ok

    d = tempfile.mkdtemp()
    with io.open(LEDGER, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh)
        fields = rd.fieldnames
        rows = [dict(r) for r in rd]
    ids = {r["id"] for r in rows}
    # (id, label, stored value, consistent display, mutated stored value)
    SYN = [("L971", "percent", "[0.00376, 133200]", "0.38% of 133,200", "[0.00476, 133200]"),
           ("L972", "x 10^k", "3.2084188911704316e-06", "p = 3.2 × 10^-6", "4.2084188911704316e-06"),
           ("L973", "1.2e-6 form", "3.15151369829389e-06", "p = 3.15e-06", "3.25151369829389e-06"),
           ("L974", "bare 10^-k", "1e-06", "p < 10^-6", "1e-05"),
           ("L975", "word count", "6", "six checks", "7"),
           ("L976", "range", "[0.124920680954014, 0.31995455552256, 6]", "0.12-0.32 (6 pairs)", "[0.124920680954014, 0.35995455552256, 6]"),
           ("L977", "sign", "-0.1910133770663204", "-0.191", "0.1910133770663204"),
           ("L978", "list", '{"Supported": 45, "Flagged": 33, "Convergent": 3}', "3 / 45 / 33", '{"Supported": 46, "Flagged": 33, "Convergent": 3}'),
           ("L979", "multi-number", "[-0.136440539213011, 0.0172824683003147, 2.9092843151557904e-15, \"concordant\"]",
            "-0.136 (SE 0.017), p = 2.9 × 10^-15; concordant", "[-0.136440539213011, 0.0272824683003147, 2.9092843151557904e-15, \"concordant\"]"),
           ("L981", "textual date", '"2026-08-08"', "8 August 2026", '"2026-08-09"'),
           ("L982", "string leaf", '["24.3-29.5"]', "F 24.3-29.5", '["24.3-30.5"]'),
           ("L983", "range percent", "[0.355, 0.66, 6]", "36-66% (6 pairs)", "[0.355, 0.76, 6]"),
           ("L984", "thousands + labels", "[14521, 38]", "14,521 PAM-T2D pairs (rs58290281, GRCh38); 38", "[14522, 38]")]
    clash = [i for i, *_ in SYN if i in ids]
    rep("synthetic ids are free in the real ledger", not clash, "(%s)" % clash)
    rep("synthetic ids are reserved test ids", {i for i, *_ in SYN} <= RESERVED_TEST_IDS)

    def write_ledger(name, edits=None, extra=()):
        p = os.path.join(d, name)
        with io.open(p, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
            w.writeheader()
            for r in rows:
                r2 = dict(r)
                r2.update((edits or {}).get(r["id"], {}))
                w.writerow(r2)
            for xid, val, disp in extra:
                r2 = {k: "" for k in fields}
                r2.update(id=xid, value=val, display=disp)
                w.writerow(r2)
        return p

    def run(ledger_path, text, extra_args=()):
        pr = os.path.join(d, "probe.md")
        io.open(pr, "w", encoding="utf-8").write(text)
        rc, fl, _ = _cli(["--sections", d, "--files", "probe.md", "--ledger", ledger_path] + list(extra_args))
        return rc, [x for x in fl if x.startswith("D1 ")]

    syn_ok = [(i, v, disp) for i, _, v, disp, _m in SYN]
    syn_mut = [(i, m, disp) for i, _, _v, disp, m in SYN]
    text_syn = " ".join("x [[%s]]" % i for i, *_ in SYN)
    text_real = "We tested 14,521 [[L004]] pairs; revised 2026-10-08 [[L918]]; dates 8 August 2026 [[L412]]."
    # consistent copy: nothing may trip
    p_ok = write_ledger("ledger_ok.csv", extra=syn_ok)
    rc, d1 = run(p_ok, text_syn + " " + text_real)
    rep("consistent rows pass (real + synthetic)", not d1 and rc in (0, 1), "(%d D1 lines)" % len(d1))
    # value-only change on real rows (the display is untouched)
    tgt = {"L004": ("value", "14521", "14522"), "L918": ("value", '"2026-10-08"', '"2026-10-09"'),
           "L412": ("value", '["2026-08-08", "2026-08-16", "2026-09-17"]', '["2026-08-08", "2026-08-17", "2026-09-17"]')}
    for rid, (col, old, new) in tgt.items():
        have = next((r[col] for r in rows if r["id"] == rid), None)
        rep("mutation target %s exists unchanged" % rid, have == old, "(%r)" % (have if have != old else ""))
        if have != old:
            continue
        p = write_ledger("ledger_%s.csv" % rid, edits={rid: {col: new}})
        rc, d1 = run(p, text_real)
        want = "[DATE]" if rid != "L004" else "[MISMATCH]"
        rep("value-only change trips D1 on %s" % rid, rc == 1 and any(x.startswith("D1 %s %s" % (rid, want)) for x in d1) and
            all(x.startswith("D1 %s " % rid) for x in d1), "(%s)" % (d1[0][:90] if d1 else "no D1 line"))
    # display-only change (value untouched)
    p = write_ledger("ledger_disp.csv", edits={"L004": {"display": "14,522"}})
    rc, d1 = run(p, text_real)
    rep("display-only change trips D1 on L004", rc == 1 and any(x.startswith("D1 L004 ") for x in d1), "(%s)" % (d1[0][:70] if d1 else "no D1 line"))
    # a row nobody tags is not file-checked, but --d1-all catches it
    p = write_ledger("ledger_untagged.csv", edits={"L004": {"value": "14522"}})
    rc, d1 = run(p, "Nothing numeric here.")
    rep("untagged mutated row: file mode silent", rc == 0 and not d1)
    rc_a, fl_a, _ = _cli(["--d1-all", "--ledger", p])
    rep("untagged mutated row: --d1-all catches it", rc_a == 1 and any(x.startswith("D1 L004 ") for x in fl_a))
    # synthetic forms: each mutation trips ITS OWN id (exact id), and only that id
    p_mut = write_ledger("ledger_mut.csv", extra=syn_mut)
    rc, d1 = run(p_mut, text_syn)
    got = {x.split()[1] for x in d1}
    for i, lab, *_ in SYN:
        rep("synthetic mutation trips %s (%s)" % (i, lab), i in got)
    rep("only the mutated ids trip", got <= {i for i, *_ in SYN} and rc == 1)
    # the real ledger itself is read by the default CLI (no --ledger): same rows as the copy
    rc_r, fl_r, _ = _cli(["--d1-all"])
    rc_c, fl_c, _ = _cli(["--d1-all", "--ledger", write_ledger("ledger_copy.csv")])
    rep("copy == real ledger under --d1-all", fl_r == fl_c and rc_r == rc_c, "(%d failing token lines on the real ledger)" % len(fl_r))
    rep("--ledger override is restored afterwards", LEDGER.endswith("S1R1_CLAIMS_LEDGER.csv"))
    return 1 if bad else 0


def selftest_s5(L):
    """Stage 5 step 1a / 1d: D1 positional matching (strict pointer, hand-declared orders), the scientific-notation floor, the uncited-row exemption (L409
    included) and its assert, the four-digit loud failure, and the reserved test ids. Every check is exercised through the real command-line path on COPIES; each mutation must trip the check it NAMES (a
    specific id and bracket code), and its target is asserted to exist and be unique (the test FAILS rather than skips)."""
    import shutil, tempfile
    bad = 0

    def rep_(lab, ok, detail=""):
        nonlocal bad
        print("selftest S5 %-52s %s %s" % (lab, "OK" if ok else "MISS", detail))
        bad += not ok

    d = tempfile.mkdtemp()
    with io.open(LEDGER, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh)
        fields = rd.fieldnames
        rows = [dict(r) for r in rd]
    byid = {r["id"]: r for r in rows}

    def write_ledger(name, edits=None, extra=()):
        p = os.path.join(d, name)
        with io.open(p, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
            w.writeheader()
            for r in rows:
                r2 = dict(r)
                r2.update((edits or {}).get(r["id"], {}))
                w.writerow(r2)
            for xid, val, disp in extra:
                r2 = {k: "" for k in fields}
                r2.update(id=xid, value=val, display=disp)
                w.writerow(r2)
        return p

    def write_order(name, lines):
        p = os.path.join(d, name)
        io.open(p, "w", encoding="utf-8", newline="").write("id\torder\tjustification\n" + "".join("%s\t%s\ttest\n" % x for x in lines))
        return p

    def d1_of(ledger_path, text, order=None, surf=None):
        pr = os.path.join(d, "probe.md")
        io.open(pr, "w", encoding="utf-8").write(text)
        args = ["--sections", d, "--files", "probe.md", "--ledger", ledger_path, "--d1-order", order or D1_ORDER_FILE]
        if surf:
            args += ["--surfaces", surf]
        rc, fl, _ = _cli(args)
        return rc, [x for x in fl if x.startswith("D1 ")]

    empty_order = write_order("order_empty.tsv", [])
    real_order = D1_ORDER_FILE
    # -- reserved ids -------------------------------------------------------------------------------------------------------
    rep_("reserved test ids: none is a real ledger row", not (RESERVED_TEST_IDS & set(byid)), "(%s)" % sorted(RESERVED_TEST_IDS & set(byid)))
    rep_("reserved test ids: L971-L984 and L999 exactly", RESERVED_TEST_IDS == set("L%d" % n for n in range(971, 985)) | {"L999"})
    # -- positional matching ---------------------------------------------------------------------------------------------------
    # (a) a real multi-leaf row: a value-only change on ONE leaf, and a swap of two leaves, each trip D1 on that row
    for rid, i, j in (("L080", 1, 2), ("L090", 0, 1)):
        row = byid.get(rid)
        vals = json.loads(row["value"]) if row else []
        good = row is not None and not d1_row(row) and isinstance(vals, list) and len(vals) >= 3 and vals[i] != vals[j] and all(isinstance(v, (int, float)) for v in vals[:3])
        rep_("target %s exists, is a clean multi-leaf list" % rid, good, "(%s | %s)" % (row["display"][:40], row["value"][:40]) if row else "missing")
        if not good:
            continue
        txt = "x [[%s]]" % rid
        v1 = list(vals)
        v1[i] = vals[i] + 1
        rc, fl = d1_of(write_ledger("ledger_%s_val.csv" % rid, {rid: {"value": json.dumps(v1)}}), txt, empty_order)
        rep_("%s: value-only change on leaf %d trips D1 [MISMATCH] on %s" % (rid, i, rid), rc == 1 and any(x.startswith("D1 %s [MISMATCH]" % rid) for x in fl) and all(x.startswith("D1 %s " % rid) for x in fl), "(%s)" % (fl[0][:80] if fl else "no D1 line"))
        v2 = list(vals)
        v2[i], v2[j] = v2[j], v2[i]
        rc, fl = d1_of(write_ledger("ledger_%s_swap.csv" % rid, {rid: {"value": json.dumps(v2)}}), txt, empty_order)
        rep_("%s: swapping leaves %d and %d trips D1 [POSITION] only" % (rid, i, j), rc == 1 and bool(fl) and all(x.startswith("D1 %s [POSITION]" % rid) for x in fl), "(%d line(s))" % len(fl))
        rc, fl = d1_of(write_ledger("ledger_%s_ok.csv" % rid), txt, empty_order)
        rep_("%s: the unchanged copy passes" % rid, rc == 0 and not fl)
    # (b) synthetic L980: "9 of 10" stored [10, 9] -- fails without a declaration, passes with it, and the declaration pins the position
    syn = [("L980", "[10, 9]", "9 of 10")]
    p_ok = write_ledger("ledger_980.csv", extra=syn)
    rc, fl = d1_of(p_ok, "x [[L980]]", empty_order)
    rep_("L980 '9 of 10' from [10, 9], no declaration: [POSITION]", rc == 1 and len(fl) == 1 and fl[0].startswith("D1 L980 [POSITION]"), "(%s)" % (fl[0][:70] if fl else "no D1 line"))
    ord_ok = write_order("order_980.tsv", [("L980", "1 0")])
    rc, fl = d1_of(p_ok, "x [[L980]]", ord_ok)
    rep_("L980 declared order '1 0': passes", rc == 0 and not fl)
    p_sw = write_ledger("ledger_980_swap.csv", extra=[("L980", "[9, 10]", "9 of 10")])
    rc, fl = d1_of(p_sw, "x [[L980]]", ord_ok)
    rep_("L980 declared order, leaves swapped: [POSITION] and nothing else", rc == 1 and bool(fl) and all(x.startswith("D1 L980 [POSITION]") for x in fl), "(%d line(s))" % len(fl))
    ord_bad = write_order("order_980_bad.tsv", [("L980", "1 0 0")])
    rc, fl = d1_of(p_ok, "x [[L980]]", ord_bad)
    rep_("L980 declaration with the wrong length: [POSITION]", rc == 1 and any(x.startswith("D1 L980 [POSITION]") and "declared display order" in x for x in fl))
    # (c) the swap mutation on the existing synthetic patterns (L979 multi-number, L983 range percent)
    for sid, v_ok, v_sw, disp in (("L979", '[-0.136440539213011, 0.0172824683003147, 2.9092843151557904e-15, "concordant"]', '[-0.136440539213011, 2.9092843151557904e-15, 0.0172824683003147, "concordant"]',
                                   "-0.136 (SE 0.017), p = 2.9 × 10^-15; concordant"),
                                  ("L983", "[0.355, 0.66, 6]", "[0.66, 0.355, 6]", "36-66% (6 pairs)")):
        rc, fl = d1_of(write_ledger("ledger_%s_a.csv" % sid, extra=[(sid, v_ok, disp)]), "x [[%s]]" % sid, empty_order)
        rep_("%s consistent copy passes" % sid, rc == 0 and not fl)
        rc, fl = d1_of(write_ledger("ledger_%s_b.csv" % sid, extra=[(sid, v_sw, disp)]), "x [[%s]]" % sid, empty_order)
        rep_("%s swapped leaves trip [POSITION] and nothing else" % sid, rc == 1 and bool(fl) and all(x.startswith("D1 %s [POSITION]" % sid) for x in fl), "(%s)" % (fl[0][:60] if fl else "no D1 line"))
    # (d) the declared orders in the shipped file all fit the shipped ledger (a stale declaration is a failure, and the table is not empty)
    declared = d1_order(real_order)
    stale = [i for i in declared if i not in byid or d1_row(byid[i]) != []]
    rep_("every declared display order fits the shipped ledger", bool(declared) and not stale, "(%d declared; stale: %s)" % (len(declared), stale))
    # (e) real rows that declare an order: the declaration pins each quantity to its leaf, so a value-only change and a leaf swap both trip
    for rid in ("L250", "L076"):
        row = byid.get(rid)
        vals = json.loads(row["value"]) if row else []
        good = row is not None and rid in declared and not d1_row(row) and isinstance(vals, list) and len(vals) >= 3 and vals[0] != vals[1]
        rep_("target %s exists, declares a display order and is clean" % rid, good)
        if not good:
            continue
        v1 = list(vals)
        v1[0] = vals[0] + 1
        rc, fl = d1_of(write_ledger("ledger_%s_dv.csv" % rid, {rid: {"value": json.dumps(v1)}}), "x [[%s]]" % rid)
        rep_("%s (declared order): value-only change on leaf 0 trips [MISMATCH]" % rid, rc == 1 and any(x.startswith("D1 %s [MISMATCH]" % rid) for x in fl) and all(x.startswith("D1 %s " % rid) for x in fl), "(%s)" % (fl[0][:70] if fl else "no D1 line"))
        v2 = list(vals)
        v2[0], v2[1] = v2[1], v2[0]
        rc, fl = d1_of(write_ledger("ledger_%s_ds.csv" % rid, {rid: {"value": json.dumps(v2)}}), "x [[%s]]" % rid)
        rep_("%s (declared order): swapping leaves 0 and 1 trips [POSITION] only" % rid, rc == 1 and bool(fl) and all(x.startswith("D1 %s [POSITION]" % rid) for x in fl), "(%d line(s))" % len(fl))
    # (f) Step 1d (S2-D28): the pointer is STRICT. L430 repeats the value 6 in leaves 2, 3 and 5; a change of leaf 2 (6 -> 7) used to be absorbed by leaf 3.
    r430 = byid.get("L430")
    v430 = json.loads(r430["value"]) if r430 else []
    good = r430 is not None and "L430" not in declared and not d1_row(r430) and isinstance(v430, list) and len(v430) == 8 and v430[2] == 6 and v430[3] == 6
    rep_("target L430 exists, is clean, undeclared, and repeats 6 in leaves 2 and 3", good, "(%s)" % (r430["value"] if r430 else "missing"))
    if good:
        v = list(v430)
        v[2] = 7
        rc, fl = d1_of(write_ledger("ledger_430_dv.csv", {"L430": {"value": json.dumps(v)}}), "x [[L430]]")
        rep_("L430: leaf 2 changed 6 -> 7 trips D1 on L430 (the next equal leaf no longer absorbs it)",
             rc == 1 and bool(fl) and all(x.startswith("D1 L430 ") for x in fl) and any(("[POSITION]" in x or "[MISMATCH]" in x) for x in fl), "(%s)" % (fl[0][:60] if fl else "no D1 line"))
        miss = []
        for k in range(len(v430)):
            v = list(v430)
            v[k] = v430[k] + 1000
            rc, fl = d1_of(write_ledger("ledger_430_leaf%d.csv" % k, {"L430": {"value": json.dumps(v)}}), "x [[L430]]")
            if not (rc == 1 and fl and all(x.startswith("D1 L430 ") for x in fl)):
                miss.append(k)
        rep_("L430: a change of ANY one of its 8 leaves trips D1", not miss, "(silent for leaves %s)" % miss)
    # (g) Step 1d: a scientific-notation token compares at a relative floor. L071 stores 7.19e-93 and 2.06e-76; the absolute 1e-15 floor let 9.8e-93 and 2.9e-76 match.
    r071 = byid.get("L071")
    v071 = json.loads(r071["value"]) if r071 else []
    good = r071 is not None and not d1_row(r071) and isinstance(v071, list) and len(v071) == 3 and abs(v071[0]) < 1e-90 and abs(v071[1]) < 1e-70
    rep_("target L071 exists, is clean, and stores p values below 1e-70", good, "(%s)" % (r071["value"][:60] if r071 else "missing"))
    if good:
        for k, newv in ((0, 9.8e-93), (1, 2.9e-76)):
            v = list(v071)
            v[k] = newv
            rc, fl = d1_of(write_ledger("ledger_071_%d.csv" % k, {"L071": {"value": json.dumps(v)}}), "x [[L071]]")
            rep_("L071: leaf %d changed to %g trips D1 on L071" % (k, newv), rc == 1 and bool(fl) and all(x.startswith("D1 L071 ") for x in fl) and any(("[MISMATCH]" in x or "[UNDERIVED]" in x) for x in fl), "(%s)" % (fl[0][:60] if fl else "no D1 line"))
    # (h) Step 1d: the five hand-declared orders. L204, L448 and L961 print one leaf twice (they fail [POSITION] without a declaration); L213 and L240 print one number for
    # two leaves, so the declaration is what pins leaf 0 (a change of leaf 0 is silent without it: the next equal leaf absorbs it).
    HAND = {"L204": lambda v: [v[0] * 3.0] + v[1:], "L448": lambda v: 0.06, "L961": lambda v: [0.06] + v[1:], "L213": lambda v: [0.5, v[1]], "L240": lambda v: [0.0002] + v[1:]}
    for rid, mut in HAND.items():
        row = byid.get(rid)
        good = row is not None and rid in declared and not d1_row(row)
        rep_("target %s exists, is hand-declared and clean" % rid, good)
        if not good:
            continue
        vals = json.loads(row["value"])
        rc, fl = d1_of(write_ledger("ledger_%s_hand.csv" % rid, {rid: {"value": json.dumps(mut(vals))}}), "x [[%s]]" % rid)
        rep_("%s: a change of leaf 0 trips D1 on %s with the declaration" % (rid, rid), rc == 1 and bool(fl) and all(x.startswith("D1 %s " % rid) for x in fl), "(%s)" % (fl[0][:60] if fl else "no D1 line"))
        rc, fl = d1_of(write_ledger("ledger_%s_hand0.csv" % rid), "x [[%s]]" % rid, empty_order)
        if rid in ("L204", "L448", "L961"):
            rep_("%s: without the declaration the unchanged row fails [POSITION]" % rid, rc == 1 and bool(fl) and all(x.startswith("D1 %s [POSITION]" % rid) for x in fl))
        else:
            rc, fl = d1_of(write_ledger("ledger_%s_hand1.csv" % rid, {rid: {"value": json.dumps(mut(vals))}}), "x [[%s]]" % rid, empty_order)
            rep_("%s: without the declaration the same leaf-0 change is SILENT (the declaration is what pins it)" % rid, rc == 0 and not fl)
    v961 = json.loads(byid["L961"]["value"]) if "L961" in byid else []
    if len(v961) == 4:
        rc, fl = d1_of(write_ledger("ledger_961_l3.csv", {"L961": {"value": json.dumps(v961[:3] + [0.06])}}), "x [[L961]]")
        rep_("L961: a change of leaf 3 (the size-difference level) trips D1", rc == 1 and bool(fl) and all(x.startswith("D1 L961 ") for x in fl))
    # -- L409: exempt as uncited (step 1d, S2-D28); every citation moved to the successor L962 ---------------------------------------------
    r409, r962 = byid.get("L409"), byid.get("L962")
    rep_("target L409 exists and is on the uncited list, with no label mechanism left", r409 is not None and "L409" in D1_UNCITED and "D1_ROW_LABELS" not in globals())
    v962 = json.loads(r962["value"]) if r962 else {}
    rep_("successor L962 exists, stores n_core_criteria = 4, and passes D1", r962 is not None and isinstance(v962, dict) and v962.get("n_core_criteria") == 4 and d1_row(r962) == [])
    rep_("L409 is cited on no operative surface (its citations are L962)", not uncited_sites().get("L409"), "(%s)" % uncited_sites().get("L409"))
    if r409 is not None:
        rc, fl = d1_of(write_ledger("ledger_409_x.csv", extra=[("L980", r409["value"], r409["display"])]), "x [[L980]]")
        rep_("L409's display on another row (L980) fails D1 on 'four': nothing but the uncited exemption passes it", rc == 1 and any(x.startswith("D1 L980 ") and "'four'" in x for x in fl), "(%s)" % (fl[0][:60] if fl else "no D1 line"))
    # -- uncited rows: exempt while uncited, [EXEMPT-CITED] once cited ----------------------------------------------------------
    rc_a, fl_a, out_a = _cli(["--d1-all"])
    # step 1d second pass: the scan reaches the tables, supplement and Declarations, so some listed rows ARE cited. The verdict must follow the scan:
    # [EXEMPT-CITED] on exactly the listed rows a non-allowed site cites, and no D1 line at all on the others (it asserts the consistency, not a pinned list).
    st_now = uncited_sites()
    cited_now = sorted(i for i in D1_UNCITED if [s for s in st_now.get(i, []) if s.rsplit(":", 1)[0] not in D1_UNCITED_PENDING.get(i, [])])
    flagged = sorted({x.split()[1] for x in fl_a if "[EXEMPT-CITED]" in x})
    rep_("--d1-all flags [EXEMPT-CITED] on exactly the listed rows the surfaces cite", flagged == cited_now, "(flagged %s)" % ",".join(flagged))
    exempt_lines = [x for x in fl_a if x.split()[1] in D1_UNCITED and x.split()[1] not in cited_now]
    rep_("a listed row nothing cites raises no D1 line on the shipped surfaces", not exempt_lines, "(%d line(s); %d of 9 cited)" % (len(exempt_lines), len(cited_now)))
    rep_("exemption list is exactly the approved 13 rows (8 of S2-D27 + L409 of S2-D28 + the 1c predecessors L256, L261, L262, L454)",
         sorted(D1_UNCITED) == sorted(["L032", "L094", "L403", "L415", "L438", "L439", "L446", "L210", "L409", "L256", "L261", "L262", "L454"]))
    plain = os.path.join(d, "surface_clean.md")
    io.open(plain, "w", encoding="utf-8").write("A surface that cites no exempt row. 14,521 [[L004]].\n")
    rc_c, fl_c, _ = _cli(["--d1-all", "--surfaces", plain])
    rep_("--surfaces: a clean surface cites nothing, and no listed row raises a D1 line (the pending allowance is empty since 1c)",
         not [x for x in fl_c if "[EXEMPT-CITED]" in x or "[EXEMPT-STALE]" in x] and not [x for x in fl_c if x.split()[1] in D1_UNCITED])
    d1_exempt_ids = sorted(D1_UNCITED)
    for rid in d1_exempt_ids:
        src = os.path.join(R2_SEC, "P31R_Results_R2_draft.md")
        text = io.open(src, encoding="utf-8", newline="").read()
        assert ("[[%s]]" % rid) not in text, "mutation target is not unique: %s already cited in the R2 Results draft" % rid
        cp = os.path.join(d, "cited_%s.md" % rid)
        io.open(cp, "w", encoding="utf-8", newline="").write(text + "\nOne more sentence with a counted item 14,521 [[L004]] and the exempt row [[%s]].\n" % rid)
        # file mode: a tag in a checked file is a citation
        rc, fl = d1_of(write_ledger("ledger_same.csv"), io.open(cp, encoding="utf-8").read(), real_order)
        rep_("%s tagged in a COPY of the R2 Results draft: [EXEMPT-CITED]" % rid, rc == 1 and any(x.startswith("D1 %s [EXEMPT-CITED]" % rid) for x in fl), "(%s)" % ("; ".join(x[:40] for x in fl if rid in x)[:80]))
        # surface mode: the bare id in an operative surface
        rc_s, fl_s, _ = _cli(["--d1-all", "--surfaces", cp])
        rep_("%s cited in a surface COPY: --d1-all [EXEMPT-CITED]" % rid, rc_s == 1 and any(x.startswith("D1 %s [EXEMPT-CITED]" % rid) for x in fl_s))
    rc_u, fl_u, _ = _cli(["--d1-all", "--surfaces", os.path.join(d, "cited_L032.md")])
    others = [x for x in fl_u if "[EXEMPT-CITED]" in x and not x.startswith("D1 L032 ")]
    rep_("citing L032 raises [EXEMPT-CITED] for L032 only", not others and any(x.startswith("D1 L032 [EXEMPT-CITED]") for x in fl_u))
    # L210's pending allowance ended at 1c (the manifest names L963): L210 is cited nowhere and the allowance is empty. The [EXEMPT-STALE] mechanism stays
    # tested: a synthetic allowance for a row that no allowed file names must report itself stale (mutation on the in-process check).
    st = uncited_sites()
    sites210 = st.get("L210", [])
    rep_("L210 is named on no operative surface and the pending allowance is empty (1c)", not sites210 and D1_UNCITED_PENDING == {}, "(%s)" % ", ".join(sites210))
    D1_UNCITED_PENDING["L032"] = ["figures/registry/no_such_file.tsv"]
    try:
        _c32, _ln32 = exempt_status("L032")
    finally:
        D1_UNCITED_PENDING.pop("L032", None)
    rep_("a synthetic allowance that no file needs fails as [EXEMPT-STALE] (mutation)", not _c32 and any("[EXEMPT-STALE]" in x and x.startswith("D1 L032 ") for x in _ln32), "(%s)" % "; ".join(_ln32)[:80])
    rep_("surface list opens every declared file", len(r2_surface_files()) >= 40 and all(os.path.isfile(f) for f in r2_surface_files()), "(%d files)" % len(r2_surface_files()))
    # -- scope of the citation scan (step 1d second pass): both axes ------------------------------------------------------------
    # axis 1, files: enumerate the shipped directories a DIFFERENT way (prefix / suffix, not the gate's regexes) and require the scan list to hold every one.
    surf = set(os.path.abspath(f) for f in r2_surface_files())
    want = set()
    for sub, pref in (("tables", "P31R_Table"), ("supplement", "P31R_S")):
        for f in os.listdir(os.path.join(ROOT, sub)):
            if f.startswith(pref) and f.endswith(".md") and "_pre" not in f and f[len(pref)] in "S0123456789MN":
                want.add(os.path.abspath(os.path.join(ROOT, sub, f)))
    want.add(os.path.abspath(os.path.join(ROOT, "sections", "P31R_Declarations.md")))
    rep_("the scan lists every shipped table, supplement file, SNote and the Declarations (independent enumeration)", bool(want) and want <= surf and len(want) >= 20,
         "(%d expected, %d missing: %s)" % (len(want), len(want - surf), sorted(os.path.basename(x) for x in want - surf)[:4]))
    rep_("the scan lists no backup (*_pre*) copy", not [f for f in surf if "_pre" in os.path.basename(f)])
    # Stage 6.2 (repoint a): the component source map's front matter is scanned
    fm = [os.path.abspath(os.path.join(ROOT, "sections", f)) for f in ("P31R_Abstract.md", "P31R_ResearchInsights.md", "P31R_Background.md", "P31R_Conclusions.md", "P31R_TitlePage.md")]
    rep_("the scan lists the source map's front-matter files (Stage 6.2 repoint a)", all(f in surf for f in fm), "(missing: %s)" % [os.path.basename(f) for f in fm if f not in surf])
    # Stage 6.2: the S1 list-line exemption is the Research Insights file only
    for nm, want_s1 in (("P31R_ResearchInsights.md", False), ("P31R_Abstract.md", True)):
        pp = os.path.join(d, nm)
        io.open(pp, "w", encoding="utf-8", newline="").write("## X\n\n**Q?**\n- a highlight line\n")
        has_s1 = any(x.startswith("S1 list/table") for x in check(pp, L)[0])
        rep_("S1: a '- ' line in %s %s" % (nm, "fails" if want_s1 else "is exempt"), has_s1 == want_s1)
        os.remove(pp)
    # axis 2, patterns: a citation of an exempt row planted in a COPY of each newly scanned file class must be found through r2_surface_files() itself.
    keep_root, keep_sec = ROOT, R2_SEC
    tmp_root = os.path.join(d, "root_copy")
    try:
        for f in r2_surface_files():
            dst = os.path.join(tmp_root, os.path.relpath(f, keep_root))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(f, dst)
        os.makedirs(os.path.join(tmp_root, "build"), exist_ok=True)   # Stage 6.2: the source map now defines part of the surface list
        shutil.copy2(os.path.join(keep_root, "build", "P31R_SOURCE_MAP.tsv"), os.path.join(tmp_root, "build", "P31R_SOURCE_MAP.tsv"))
        globals()["ROOT"], globals()["R2_SEC"] = tmp_root, os.path.join(tmp_root, "sections_R2")
        base_sites = uncited_sites().get("L032", [])
        rep_("scan copy: L032 is named in no surface to begin with (mutation target is new)", not base_sites, "(%s)" % base_sites)
        for rel in ("tables/P31R_Table1.md", "tables/P31R_Table2.md", "supplement/P31R_SM_kidney.md", "supplement/P31R_SNote.md", "sections/P31R_Declarations.md"):
            fp = os.path.join(tmp_root, rel)
            orig = io.open(fp, encoding="utf-8", newline="").read()
            io.open(fp, "w", encoding="utf-8", newline="").write(orig + "\nPlanted citation [[L032]].\n")
            hits = uncited_sites().get("L032", [])
            rep_("scan copy: L032 planted in %s is found" % rel, len(hits) == 1 and hits[0].startswith(rel + ":"), "(%s)" % hits)
            io.open(fp, "w", encoding="utf-8", newline="").write(orig)
        bk = os.path.join(tmp_root, "supplement", "P31R_SM_kidney_preS5S1d.md")
        io.open(bk, "w", encoding="utf-8", newline="").write("Backup copy [[L032]].\n")
        rep_("scan copy: a *_preS5S1d.md backup beside a supplement file is not a surface", not uncited_sites().get("L032"))
    finally:
        globals()["ROOT"], globals()["R2_SEC"] = keep_root, keep_sec
    # -- four-digit tags are loud ---------------------------------------------------------------------------------------------
    pr = os.path.join(d, "probe4.md")
    io.open(pr, "w", encoding="utf-8").write("We tested 14,521 [[L1000]] pairs.\n")
    fails4, _ = check(pr, L)
    rep_("a [[L1000]] tag fails loudly as N3 (four-digit)", any(x.startswith("N3 four-digit ledger tag L1000") for x in fails4) and not any("unknown tag L1000" in x for x in fails4), "(%s)" % ("; ".join(x[:50] for x in fails4)[:90]))
    try:
        scan_pairs("so 14,521 [[L004]] and 7 [[L1000]]")
        loud = False
    except AssertionError as e:
        loud = "four-digit" in str(e)
    rep_("scan_pairs refuses a four-digit tag loudly", loud)
    rep_("scan_pairs reads three-digit pairs", scan_pairs("so 14,521 [[L004]] and 7 [[L123]]") == [("14,521", "[[L004]]", "L004"), ("7", "[[L123]]", "L123")])
    return 1 if bad else 0


def selftest():
    """Mutations must trip the check they name."""
    L = ledger()
    import tempfile
    cases = [("untagged", "We tested 14,521 pairs.", "N1"), ("wrong value", "We tested 14,522 [[L004]] pairs.", "N2"),
             ("unknown tag", "We tested 14,521 [[L999]] pairs.", "N3"), ("banned", "The colocalization gate was strict 14,521 [[L004]].", "V1"),
             ("replicat", "deCODE replicated 23 [[L031]] pairs.", "V1"), ("ok", "We tested 14,521 [[L004]] pairs.", None),
             ("ok sci", "Instruments reached p < 10⁻⁶ [[L230]].", None), ("ok list", "Tiers were 3 / 45 / 33 [[L014]].", None),
             ("sign flip", "The ratio was +0.191 [[L050]].", "N2"), ("double round", "B gave +0.040 [[L055]].", "N2"),
             ("bullet", "Text.\n- a bullet item", "S1")]
    bad = 0
    for lab, txt, want in cases:
        d = tempfile.mkdtemp()
        p = os.path.join(d, "X.md")
        io.open(p, "w", encoding="utf-8").write(txt)
        fails, _ = check(p, L)
        codes = {x.split()[0] for x in fails}
        ok = (want is None and not codes) or (want is not None and want in codes)
        print("selftest %-12s want %-5s got %-14s %s" % (lab, want, ",".join(sorted(codes)) or "-", "OK" if ok else "MISS"))
        bad += not ok
    r2bad = selftest_r2(L)
    d1bad = selftest_d1(L)
    s5bad = selftest_s5(L)
    return 1 if (bad or r2bad or d1bad or s5bad) else 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else main())
