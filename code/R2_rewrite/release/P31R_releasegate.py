#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 R2, Stage 5 steps 3-4 and RC1: the release gate (R0-R11) and its self-test.

  python -I -X utf8 P31R_releasegate.py [--release DIR] [--v3 DIR] [--manuscript PATH] [--selftest]

Reads only; writes nothing in the release (the self-test works in a scratch COPY).  Exit 0 only when every check is green.
Check ids (matched on letter + number; R3 is never satisfied by R3b):
  R0  scope and reader: the gate opened every file in out/ (count printed); the text reader recovers a known string from a real shipped
      file and returns a plausible number of characters per byte (two-sided band: neither blind nor flooded).
  R1  every file in out/ is in RELEASE_MANIFEST.tsv and every manifest row exists (the manifest, README.md, REPRODUCE.md are exempt, and
      must not have rows of their own).
  R2  every md5 / byte count matches, every source md5 matches its source file now, rows/columns of csv/tsv agree, a producer script that
      resolves matches its producer_md5.
  R3  no md5 listed in inputs/NOT_SHIPPED.tsv is present in out/ or withheld/ (path + md5 list, never a folder rule).
  R3b each unit REPRODUCE.md states has the class derived in build/regenerability_by_unit.tsv; every command's inputs are in out/ or are a
      declared external download / non-redistributable input.
  R4  no individual-level genotype file (.i8 .i2 .bed .bim .fam .pgen .pvar .psam .bgen .vcf) and no samples x variants shaped table.
  R5  secret scan (BLOCKING: JWT, API key / token / password literals, private keys, Bearer headers; FLAGGED, non-blocking: absolute personal
      or drive paths, e-mail addresses).  Covers every file in out/ and every ship=1 file of inputs/CODE_INCLUDE.tsv.
  R6  every file name the manuscript names exists in out/ (or is declared external in _gate_fixtures/R6_EXTERNAL_NAMES.tsv).  The default reads the built
      R2_rewrite/build/P31R_Manuscript_R2_preview.docx (reader liveness: a known string comes back and text chars per XML byte lies in DOCX_BAND).
      --manuscript takes a .docx or a folder/file of .md.
  R7  every number printed in README.md / REPRODUCE.md sits in a paragraph (or table row) that carries a ledger row id, and the number is
      in that row's value or display.  Nothing is typed.
  R8  the screen: rows == L003, testable == L004, BH column recomputed from p_family == fdr_sig_bh_family == L008 admitted, legacy flag == L009;
      a conservative-prior column ships only in a declared file and column (CONS_PERMITTED) and only on rows whose (gene, outcome) is one of the 81 pairs of the
      shipped canonical tiered table (= L010, = its source file), plus, for Table 2, the pairs the manuscript's own Table 2 prints (bound by IDENTITY, not row count).
  R9  allowlist: every file in out/ maps to an approved inclusion row (inputs/RELEASE_INPUTS.tsv ship_status ship*, which includes the declared module outputs, or
      inputs/CODE_INCLUDE.tsv ship=1); every approved row has its file; unknown entries fail, including every member of every archive (zip / tar, nested to
      R9_ARCHIVE_DEPTH, sniffed by content) unless listed in _gate_fixtures/R9_ARCHIVE_MEMBERS.tsv (member = 'archive::inner::file').
  R10 ledger recount: ledger rows L660-L676 (read from release inputs) equal a recount of their input now (own parser, not ext_r2_release.py).  R10_PENDING declares the
      one known difference (L676: 36 in the ledger, 35 in PROVIDERS.tsv, 'fixed in the Stage 6 ledger batch'); it fails when the recount matches the ledger (stale)
      or when it or any other row differs from the declaration (drift).
  R12 portability, executables (BLOCKING, RC1 fix K01): no shipped *.R / *.py / *.sh ANYWHERE in out/ names an absolute <C>/Users, D:/ or E:/ path in a non-comment line
  R13 portability, every shipped text file (BLOCKING, RC1 fix K02): no shipped text file contains the author's home path or a drive-letter path
  R11 tiers: build/tiers_by_file.tsv (file tier providers decision reason) covers every shipped file; nothing ships against inputs/TIER_RULES.tsv
      (rule_id kind target tier providers handling selector basis reason).  Content-judged, see r11_tiers().  Tier 3 never ships while switch TPMI_DERIVED is 'withheld'.
  R14 withheld-value leak scan (BLOCKING, RC1b fix Q01): every withheld cell value at >= 4 significant digits (a decimal or exponent token; cells of withheld/*.csv|*.tsv that differ from
      the shipped counterpart, or all cells when no counterpart or row count differs; plus the values of the masked ledger rows) is searched, as a number (so 0.0258 and 2.58e-2 are one key),
      in EVERY shipped text file of out/ (code/ included).  A hit fails unless inputs/LEAKSCAN_COINCIDENCES.tsv declares it (out_path, value, withheld_source, reason).
      A declared row that no longer occurs is flagged (non-blocking: a rebuilt release may drop it).
      RC2 (S2-D35(3)): the scan runs on EVERY field of every shipped text file (notes columns and the public ledger included) after numeric normalisation: a token of >= 6 significant digits matches a
      withheld value of at least as many digits, rounded to the token's own precision (the same p printed at 6 and at 7 digits is one value); and ROW-LOCALLY, a withheld cell (or a masked ledger row's private value) is compared at the shorter token's precision (floor 4 digits)
      with every other cell of the same shipped row.  Compared at 4 digits across the whole release the scan meets about 750 chance coincidences (5 digits: 49; 6 digits: 6), which is why the global floor is 6.

  R16 derived columns (BLOCKING, S2-D40 W1): a shipped column must not be a deterministic function of withheld columns.  For every withheld table with a shipped counterpart of the same row
      count, every beta-like x se-like pair of the WITHHELD file gives z = beta/se; no shipped column of the counterpart may equal |z|, z^2 (= F, the chi-square) or the two-sided normal p of z
      (relative error 1e-4 on at least half of its numeric rows, at least three rows).  The check reads the data, not TIER_RULES.tsv: a rule table that omits a derived column (F01 omitted F_stat) cannot hide it.
  R17 masked rows (BLOCKING, S2-D40 W2): in every shipped table that has a notes / note / comment column, a row whose value / se / p / ci_lo / ci_hi cell is a withheld marker has a notes cell that is
      empty, NA or the same marker (a notes cell quotes the values it sits beside).
  R16 (S2-D41 W3) also reads every withheld F-like column (F, F_stat, lead_F, decode_F_stat): z^2 = F, so no shipped column may equal F, |z| = sqrt(F) or the normal p of z, and F must not be recoverable from
      shipped columns either: with b, se1 and se2 shipped, F = b^2 / (se2^2 - se1^2) is the delta-method identity (deCODE: delta_se^2 = decode_se^2 + decode_b^2 / F).
  R18 provider values in free text (BLOCKING, S2-D41 W4): in every shipped table with a notes / note / comment column, (a) no numeric token of >= 4 significant digits equals a withheld value unless the declaration is
      S2-D35 or published_quote (a value public in a cited paper's supplement) (a quoted_value declaration does not excuse a notes cell: an unmasked row's notes may not quote a withheld provider value), and (b) no token of >= 5 significant digits equals, at its own precision,
      the |z| = |value / se| (or |beta / se|) of a row whose value / se cells are withheld in a table that ships a counterpart.
  R19 the journal files (BLOCKING, S2-D43): no number in the built manuscript docx or Additional files 1 and 2 equals a value the release withholds and does not ship
      elsewhere, unless _gate_fixtures/R19_JOURNAL_DECLARED.tsv declares it (PENDING-USER declarations pass with a warning); every hold-H9 pair sits in Table S7a and
      prints no number there except its p.

Contracts the writers of RELEASE_MANIFEST.tsv / REPRODUCE.md / README.md are held to (documented here because the gate parses them):
  * manifest columns exactly: path bytes md5 rows columns source_path source_md5 producer_script producer_md5 software_label cited_by licence_class;
    path is relative to out/; a multi-source row separates source_path and source_md5 with ';'; source_path is relative to V3 (or absolute).
  * REPRODUCE.md: one unit and one class id per line (table row or 'unit: class'); the class is printed as the id from
    regenerability_by_unit.tsv (regenerates_from_release / needs_public_download / needs_nonredistributable_input /
    does_not_regenerate_engine_output; the spaced plan wording is also read).  Commands sit in ``` fences, one per line, starting
    'python -I' or 'Rscript --vanilla', run from the repository root (= out/); output paths follow --out/--output/-o/--outdir or '>'.
  * README/REPRODUCE numbers: each numeral outside code spans sits in a paragraph with a ledger id (L###); dates, versions and ids exempt.
"""
import argparse, base64, collections, csv, gzip, hashlib, io, json, math, os, re, shlex, shutil, sys, tarfile, tempfile, zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST_COLS = ["path", "bytes", "md5", "rows", "columns", "source_path", "source_md5", "producer_script", "producer_md5",
                 "software_label", "cited_by", "licence_class", "provenance_label", "transform"]
SELF_EXEMPT = {"RELEASE_MANIFEST.tsv", "README.md", "REPRODUCE.md", "LICENSE", "CITATION.cff", ".gitattributes"}
MANUSCRIPT_DOCX = "R2_rewrite/build/P31R_Manuscript_R2_preview.docx"
DOCX_BAND = (0.30, 0.90)      # text chars per byte of word/document.xml: the preview docx reads 0.59 (markup is ~40%% of the XML); a blind reader gives ~0, a flooded one > 1
LEDGER_REL = "R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv"
CLASSES = ["regenerates_from_release", "needs_public_download", "needs_nonredistributable_input", "does_not_regenerate_engine_output"]
CLASS_ALIASES = {"needs_a_public_download": "needs_public_download", "needs_a_non_redistributable_input": "needs_nonredistributable_input",
                 "needs_non_redistributable_input": "needs_nonredistributable_input", "does_not_regenerate": "does_not_regenerate_engine_output"}
INDIV_EXT = (".i8", ".i2", ".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam", ".bgen", ".vcf")
MAX_COLS = 200            # R4 declared band: the widest shipped table has 77 columns; a samples x variants matrix has hundreds to thousands
MAX_VARIANT_HEADERS = 20  # R4: more than this many rsID / chr:pos column names in one header = genotype matrix
CODE_EXT = {".py", ".r", ".sh", ".bat", ".md", ".txt", ".json", ".tsv", ".csv", ".yml", ".yaml", ".cfg", ".ini", ".log"}
FILE_TOKEN = re.compile(r"(?<![\w/.\-])[A-Za-z0-9_][A-Za-z0-9_./\-]*\.(?:csv|tsv|json|xlsx|txt|log|py|R|md|docx)(?:\.gz)?(?![\w])")

SECRET_BLOCK = [
    ("JWT", re.compile(r"eyJ[A-Za-z0-9_\-]{4,}\.[A-Za-z0-9_\-]*\.[A-Za-z0-9_\-]*")),      # any header followed by two dot-separated segments of any length (V02)
    ("private key", re.compile(r"-----BEGIN (?:RSA |EC |DSA |OPENSSH |PGP |ENCRYPTED )?PRIVATE KEY")),
    ("Bearer header", re.compile(r"\bBearer\s+[A-Za-z0-9._~+/=\-]{16,}")),
    ("AWS/GitHub key", re.compile(r"\b(?:AKIA[0-9A-Z]{16}|gh[pousr]_[A-Za-z0-9]{30,})\b")),
    ("key/token/password literal", re.compile(
        r"(?i)\b(?:api[_\-]?key|apikey|auth[_\-]?token|access[_\-]?token|synapse[_\w\-]*(?:token|auth)|ldlink[_\w\-]*token|ncbi[_\w\-]*key|"
        r"token|secret|passwd|password|pwd)\b[\"']?\s*(?:[:=]|<-)\s*[\"']?"
        r"(?!(?:Sys\.getenv|os\.environ|os\.getenv|getenv|process\.env|environ)\b)(?=[^\s\"',;)(]{8,}(?:[\s\"',;)]|$))[^\s\"'$,;)(]{8,}")),    # quotes optional (V02)
]
SECRET_FLAG = [
    ("absolute path", re.compile(r"(?:\b[A-Za-z]:[\\/](?=\S)|/Users/\w|/home/\w)")),
    ("e-mail", re.compile(r"[A-Za-z0-9._%+\-]+@[A-Za-z0-9\-]+\.[A-Za-z]{2,}")),
]


LOCAL_PATH_RX = re.compile(r"(?<![A-Za-z0-9_])[A-Z]:[/\\]{1,2}(?=[A-Za-z0-9_$.~])(?!\\?[ntr](?![A-Za-z0-9_$.~/\\]))|Users[/\\]+" + re.escape(os.path.basename(os.path.expanduser("~"))))      # S6.3b2a: a message escape such as an R-colon-backslash-n is not a drive path


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class Res:
    def __init__(self):
        self.fails, self.notes, self.warns = {}, {}, []
        self.order = ["R0", "R1", "R2", "R3", "R3b", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16", "R17", "R18", "R19"]

    def fail(self, cid, msg):
        self.fails.setdefault(cid, []).append(msg)

    def note(self, cid, msg):
        self.notes.setdefault(cid, []).append(msg)

    def tripped(self):
        return set(self.fails)

    def report(self, out=print):
        for cid in self.order:
            f = self.fails.get(cid, [])
            n = "; ".join(self.notes.get(cid, []))
            out("%-3s %s  %s" % (cid, "FAIL (%d)" % len(f) if f else "PASS", n))
            for m in f[:8]:
                out("      - " + m)
            if len(f) > 8:
                out("      - ... %d more" % (len(f) - 8))
        for w in self.warns[:40]:
            out("FLAG " + w)
        if len(self.warns) > 40:
            out("FLAG ... %d more flagged (non-blocking)" % (len(self.warns) - 40))
        return not self.fails


# ---------------------------------------------------------------- readers
def read_text(p):
    """Decode a shipped file to text; None if binary (NUL bytes).  Never silently empty: raises on unreadable files."""
    raw = Path(p).read_bytes()
    if str(p).endswith(".gz"):
        raw = gzip.decompress(raw)
    if b"\x00" in raw[:8192]:
        return None
    return raw.decode("utf-8", errors="replace")


def docx_text(p):
    with zipfile.ZipFile(p) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    paras = re.findall(r"<w:p[ >].*?</w:p>", xml, flags=re.S)
    out = ["".join(re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", x)) for x in paras]
    return "\n".join(out), len(xml)


def manuscript_default(v3):
    r2 = v3 / "R2_rewrite"
    fs = []
    for pat, d in (("P31R_Methods_R2_draft.md", "sections_R2"), ("P31R_Results_R2_draft.md", "sections_R2"), ("P31R_Discussion_R2_draft.md", "sections_R2")):
        fs.append(r2 / d / pat)
    fs += sorted((r2 / "sections_R2").glob("P31R_Legend_*_R2.md"))
    fs += sorted((r2 / "tables").glob("P31R_Table*.md"))
    fs += sorted((r2 / "supplement").glob("P31R_S*.md"))
    fs.append(r2 / "sections" / "P31R_Declarations.md")
    return [f for f in fs if f.exists()]


def r6_docx_sources(v3):
    """The files the built manuscript docx is made from: the canonical files of the component source map, except rows of state 'generated'.
    Those rows are another file's builder or output (references, the Additional file 1 builder, the Additional file 3 legends), and the
    manuscript builder reads none of them, so a newer one does not make the docx stale (Stage 6.3b part 2b: an edit to P31R_build_af1.py
    had turned R6 red)."""
    fs = []
    for r_ in read_tsv(v3 / "R2_rewrite" / "build" / "P31R_SOURCE_MAP.tsv"):
        if r_["state"] == "generated":
            continue
        for part_ in r_["canonical_file"].split(";"):
            fs += [h for h in (v3 / "R2_rewrite").glob(part_.strip()) if not re.search(r"_pre[A-Za-z0-9_]*\.md$", h.name)]
    return fs


def r6_stale(docx, sources):
    """The newest of the docx's sources when it is newer than the docx (the docx is stale), else None."""
    newest = max(sources, key=lambda h: h.stat().st_mtime)
    return newest if docx.stat().st_mtime < newest.stat().st_mtime else None


# ---------------------------------------------------------------- ledger
def load_ledger(v3, ledger_path=None):
    p = Path(ledger_path) if ledger_path else v3 / LEDGER_REL
    with open(p, newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def nums_in(s):
    return {m.replace(",", "") for m in re.findall(r"\d[\d,]*(?:\.\d+)?", s or "")}


LABEL_RX = re.compile(r"(?i)\b(?:Tables?|Figures?|Fig\.|Stage|Step|Paper|Section|Decision|Additional file)s?\s+(?<!\d)\d+(?:\.\d+)?(?![\d])"
                      r"(?!\s+(?:pairs?|rows?|files?|columns?|cells?|scripts?|genes?|variants?|signals?|loci|outcomes?|proteins?|tests?|hits?|entries|records?))")


def _blank(rx, text, flags=0):
    return re.sub(rx, lambda m: " " * len(m.group(0)), text, flags=flags)


def r7_para(R, fn, para, ledger):
    """R7 for one paragraph / table row (V03). Every printed numeral must lie INSIDE the verbatim display of a ledger row cited in the same paragraph.
    Only where NO cited row's display occurs in the paragraph at all is the number-set reading used, and then the number must be a whole token of
    a cited row's value / display in the same position context (the word before it). Labels (Table 2, Stage 5, ...) are exempt only when no count
    noun follows. Returns the number of numerals examined."""
    ids = set(re.findall(r"\bL\d{3}\b", para))
    cited = []
    for i in sorted(ids):
        if i in ledger:
            cited.append(i)
        else:
            R.fail("R7", "%s cites %s, which is not a ledger row" % (fn, i))
    t = _blank(r"<!--.*?-->", para, re.S)
    t = _blank(r"`[^`]*`|\]\([^)]*\)|https?://\S+", t)
    t = _blank(r"\b\d{4}-\d{2}-\d{2}\b|\b\d+(?:\.\d+){2,}\b|\bv\d[\w.]*", t)
    t = _blank(r"(?m)^\s*(?:#+\s*)?\d+[.)]\s", t)
    t = re.sub(LABEL_RX, lambda m: " " * len(m.group(0)), t)
    t = _blank(r"\bship\s*=\s*\d", t)
    t = _blank(r"\bL\d{3}(?:\s*[-–]\s*L?\d{3})?\b", t)
    spans = []
    for i in cited:
        disp = (ledger[i]["display"] or "").strip()
        if not disp:
            continue
        rx = re.compile(r"\s+".join(re.escape(w) for w in disp.split()))
        spans += [m.span() for m in rx.finditer(para)]
    n = 0
    for m in re.finditer(r"(?<![\w./:\-])\d[\d,]*(?:\.\d+)?(?![\w])", t):
        n += 1
        tok = m.group(0).rstrip(",")
        num = tok.replace(",", "")
        if any(a <= m.start() and m.end() <= b for a, b in spans):
            continue
        ctx_txt = t[max(0, m.start() - 25):m.end() + 15].replace("\n", " ")
        if spans or not cited:
            R.fail("R7", "%s: number %s is not inside the display of a cited ledger row (ids here: %s): ...%s..." % (fn, tok, sorted(ids) or "none", ctx_txt))
            continue
        prev = re.search(r"([A-Za-z][\w+\-]*)\W{0,3}$", t[:m.start()])
        good = False
        for i in cited:
            v = (ledger[i]["value"] or "") + " " + (ledger[i]["display"] or "")
            if num not in nums_in(v):
                continue
            if prev is None or re.search(r"(?<![\w])" + re.escape(prev.group(1)) + r"\W{0,3}" + re.escape(tok) + r"(?![\w])", v):
                good = True
        if not good:
            R.fail("R7", "%s: number %s has no ledger row beside it in the same position (ids here: %s): ...%s..." % (fn, tok, sorted(ids) or "none", ctx_txt))
    return n


# ---------------------------------------------------------------- RC1 additions: R8 identity, R9 allowlist, R10 recount, R11 tiers
CONS_COL = re.compile(r"(?i)conservative|posterior_cons|cons_prior")
# R8: the only columns that may carry a conservative-prior value, per shipped file (declared; decision 6 as relaxed by the manuscript's own tables)
CONS_PERMITTED = {
    "source_data/figures/Figure_2/fig2a_plotted.tsv": {"PP.H4_conservative"},
    "tables/Table2_evidence_cards.tsv": {"pp_h4_conservative", "reg_pp_h4_conservative", "eas_pp_h4_conservative"},
}
CANON_TIERED = "colocalization/table_s_pair_evidence_tiered.csv"
TABLE2_MD = "R2_rewrite/tables/P31R_Table2.md"
# S2-D34(5): the three Table 2 cards that lie outside the 81 pairs (FGFR4-T2D, UMOD-eGFR, DPEP1-eGFR) may carry a conservative-prior value. This is a DECLARED list, by pair and column, not
# a rule that reads the pairs from the manuscript: file -> {(gene, outcome): (columns that may hold a value, reason)}. Every other pair keeps the canonical-81 identity rule; an undeclared
# column of a declared pair, or a value in a declared pair with an empty column set, fails. Each declared pair must still be one that the manuscript's Table 2 prints (when the file is reachable).
CONS_DECLARED = {
    "tables/Table2_evidence_cards.tsv": {
        ("FGFR4", "T2D"): ({"pp_h4_conservative"}, "printed in the manuscript's Table 2 (conservative PP.H4 of the East Asian discovery colocalization, ledger L215)"),
        ("UMOD", "eGFR"): (set(), "declared by S2-D34(5); the pair is printed in Table 2, but the current Table 2 text does not print the East Asian conservative PP.H4 (ledger L075), so no conservative-prior column holds a value"),
        ("DPEP1", "eGFR"): (set(), "declared by S2-D34(5); the pair is printed in Table 2 and no conservative-prior column holds a value"),
    },
}
CONS_NO_VALUE = ("", "not reported", "NA")
R9_ARCHIVE_DEPTH = 4
OFFICE_EXT = (".docx", ".xlsx", ".pptx", ".xlsm", ".odt", ".ods", ".odp")
TIER_DECISION_SHIP = re.compile(r"(?i)^(ship\w*|binned\w*|reduced\w*)$")
TIER_DECISION_WITHHOLD = re.compile(r"(?i)^withh\w*$")
TIERS_COLS = ["file", "tier", "providers", "decision", "reason"]
# R11 content detector: a table whose header has a variant-identifier column AND a value column carries per-variant values (Tier 2 at least)
VAR_ID = re.compile(r"(?i)^(rs_?id\w*|snp|snp_?id|variant|variant_?id|signal_variant|marker(?:name)?|chr_?pos|chrpos|pos|position|bp|ukb_variant_id)$")
VAR_VAL = re.compile(r"(?i)(^|_)(beta|se|stderr|pval\w*|p|p_value|maf|eaf|af|freq\w*|z|zscore|or|logp|log10p|nlog\w*|lp)(_|$)")
MASK_TOKEN = "withheld_pending_terms"      # the cell value rel_tables / rel_licence write for a masked value
EAS_EUR_COL = re.compile(r"(?i)^(eur_.*|het_p|.*strict.*)$")


def split_ids(s):
    return [x.strip() for x in re.split(r"[+;,|]", s or "") if x.strip()]


def is_per_variant(header):
    hs = [h.strip().strip('"') for h in header]
    return any(VAR_ID.match(h) for h in hs) and any(VAR_VAL.search(h) for h in hs if not VAR_ID.match(h))


def cc_ok_providers(release):
    """Providers whose terms were read from an official page AND name CC0 or CC BY (not -SA/-NC/-ND): the only providers whose per-variant values may ship."""
    ok = set()
    for r in read_tsv(Path(release) / "inputs" / "PROVIDER_TERMS.tsv"):
        if r.get("status", "").strip() == "confirmed_from_official_page" and re.match(r"(?i)^(CC0|CC[ -]BY)(?![- ]?(?:SA|NC|ND))\b", r.get("licence_name", "").strip()):
            ok.add(r["provider_id"].strip())
    return ok


def canonical_pairs(v3, shipped_text, man_rows):
    """The 81 pairs of the canonical tiered table, from the shipped copy; cross-checked against its source file and the ledger count (L010) when reachable."""
    rows = list(csv.DictReader(io.StringIO(shipped_text, newline="")))
    pairs = {(r["gene"].strip(), r["outcome"].strip()) for r in rows}
    return pairs, len(rows)


def table2_md_pairs(v3):
    p = Path(v3) / TABLE2_MD
    if not p.is_file():
        return None
    return {(m.group(1), m.group(2)) for m in re.finditer(r"(?m)^\|\s*\*\*([A-Za-z0-9]+)\s*[–—-]\s*([A-Za-z0-9]+)\*\*", p.read_text(encoding="utf-8"))}


def _archive_kind(name, data):
    """'zip' / 'tar' / '7z-rar' / None for the bytes of a file (content sniff, not extension); OOXML / ODF documents are not archives."""
    bio = io.BytesIO(data)
    if data[:6] == b"7z\xbc\xaf\x27\x1c" or data[:4] == b"Rar!":
        return "7z-rar"
    if zipfile.is_zipfile(bio):
        z = zipfile.ZipFile(io.BytesIO(data))
        names = set(z.namelist())
        if name.lower().endswith(OFFICE_EXT) and ("[Content_Types].xml" in names or "mimetype" in names):
            return None
        return "zip"
    try:
        tarfile.open(fileobj=io.BytesIO(data), mode="r:*").close()
        return "tar"
    except (tarfile.TarError, EOFError, OSError):
        return None


def archive_members(name, data, depth=0):
    """Yield (chained member path, error) for every member of an archive, recursing into nested archives; an archive that cannot be listed or is nested too deep yields an error."""
    kind = _archive_kind(name, data)
    if kind is None:
        return
    if kind == "7z-rar":
        yield name, "archive format that the gate cannot list"
        return
    if depth >= R9_ARCHIVE_DEPTH:
        yield name, "archive nested deeper than %d levels" % R9_ARCHIVE_DEPTH
        return
    members = []
    if kind == "zip":
        z = zipfile.ZipFile(io.BytesIO(data))
        for zi in z.infolist():
            if not zi.is_dir():
                members.append((zi.filename, lambda zi=zi, z=z: z.read(zi)))
    else:
        t = tarfile.open(fileobj=io.BytesIO(data), mode="r:*")
        for ti in t.getmembers():
            if ti.isfile():
                members.append((ti.name, lambda ti=ti, t=t: t.extractfile(ti).read()))
            elif not ti.isdir():
                yield name + "::" + ti.name, "non-regular archive entry (link / device)"
    for mname, rd in members:
        chained = name + "::" + mname
        yield chained, None
        raw = rd()
        yield from archive_members(chained, raw, depth + 1)


def r9_allowlist(R, release, rel):
    ri = read_tsv(Path(release) / "inputs" / "RELEASE_INPUTS.tsv")
    ok, no = {}, {}
    for r in ri:
        rp = (r.get("release_path") or "").strip()
        if not rp:
            continue
        (ok if (r.get("ship_status") or "").strip().lower().startswith("ship") else no).setdefault(rp, "RELEASE_INPUTS %s" % r.get("input_id", "?"))
    ci = Path(release) / "inputs" / "CODE_INCLUDE.tsv"
    ncode = 0
    for r in (read_tsv(ci) if ci.exists() else []):
        if str(r.get("ship", "")).strip().lower() in ("1", "true", "yes") and (r.get("ship_path") or "").strip():
            ok.setdefault(r["ship_path"].strip(), "CODE_INCLUDE")
            ncode += 1
    decl = {}
    fx = Path(release) / "_gate_fixtures" / "R9_ARCHIVE_MEMBERS.tsv"
    if fx.exists():
        for r in read_tsv(fx):
            decl[r["member"].strip()] = r.get("reason", "")
    nunk = narch = nmem = 0
    for k, p in sorted(rel.items()):
        if k in SELF_EXEMPT:
            continue
        if k not in ok:
            nunk += 1
            R.fail("R9", "%s in out/ maps to no approved inclusion row (RELEASE_INPUTS ship*, CODE_INCLUDE ship=1)%s" % (
                k, "; its row says " + no[k] + " (not shipped)" if k in no else ""))
        data = p.read_bytes()
        found = False
        for mname, err in archive_members(k, data):
            found = True
            nmem += 1
            if err:
                R.fail("R9", "%s: %s" % (mname, err))
            elif mname not in decl:
                R.fail("R9", "archive member not declared in _gate_fixtures/R9_ARCHIVE_MEMBERS.tsv: %s" % mname)
        narch += found
    for k in sorted(set(ok) - set(rel)):
        R.fail("R9", "approved inclusion row has no file in out/: %s (%s)" % (k, ok[k]))
    R.note("R9", "%d files in out/ against %d approved paths (%d RELEASE_INPUTS ship*, %d CODE_INCLUDE ship=1); %d archive(s) listed, %d member(s)" % (
        len(rel) - len(SELF_EXEMPT & set(rel)), len(ok), len(ok) - ncode, ncode, narch, nmem))


# ---- R10: recount of every ledger row that is read from a release input (L660-L675, L533, L676). Independent of ext_r2_release.py: stdlib only, own parsing.
# S6.3b2a: L533 [35, 14] is recounted from the LIVE inputs/PROVIDERS.tsv, exactly; L676 [36, 14] from its FROZEN, md5-pinned input (the 1c list). The pending exception for L676 is retired.
R10_IDS = ["L%d" % i for i in range(660, 676)] + ["L533", "L676"]
L676_FROZEN_MD5 = "e05944719a56be320717c6fa67f4f43d"      # R1_stage1/S1R1_06_claims_ledger/frozen_inputs/PROVIDERS_at_1c_L676.tsv (the file the ledger row names)
# Declared pending exceptions: ledger value that is KNOWN to differ from the recount until the named fix lands. The exception holds only while
# ledger == the stated ledger value AND recount == the stated recount value; a recount that matches the ledger (stale exception) or any other
# difference (drift) fails.
R10_PENDING = {}      # retired (S6.3b2a): no ledger row differs from its recount. The mechanism stays: a re-declared exception that the recount already matches fails as stale.


def _csv_rows(path, delim=","):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter=delim))


def _body_rows(path):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    sep = [i for i, l in enumerate(lines) if re.fullmatch(r"\|(?:\s*:?-{3,}:?\s*\|)+", l.strip())]
    if len(sep) != 1:
        raise ValueError("%s: expected one table separator, found %d" % (path, len(sep)))
    body = []
    for l in lines[sep[0] + 1:]:
        if not l.startswith("|"):
            break
        body.append(l)
    return body


def r10_recount(v3, release):
    """-> {ledger id: (recounted value | None, error | None)}. Reads the files now; V3-wide sources from v3, release-declared lists from `release`."""
    v3, release = Path(v3), Path(release)
    csv.field_size_limit(10 ** 9)
    res = {}

    def put(i, fn):
        try:
            res[i] = (fn(), None)
        except Exception as e:      # a missing / reshaped source is a failure of the row, reported, never skipped
            res[i] = (None, "%s: %s" % (type(e).__name__, e))
    tfr = v3 / "V4" / "source_data" / "tables_from_rds"
    mr_p = tfr / "phase2_mr_screen__mr_all_results.csv"
    dep = v3 / "paper3_deposit_request_v4"
    PRIM = ("Wald ratio", "Inverse variance weighted")

    def mr():
        return _csv_rows(mr_p)
    put("L660", lambda: len(mr()))

    def l661():
        np_ = [r for r in mr() if r["method"] not in PRIM]
        return [len(np_), len({(r["gene"], r["outcome_name"]) for r in np_})]
    put("L661", l661)

    def l662():
        px = [r for r in mr() if r["any_proxy"].strip().upper() == "TRUE"]
        return [len(px), len({(r["gene"], r["outcome_name"]) for r in px})]
    put("L662", l662)
    put("L663", lambda: len(_csv_rows(tfr / "phase3_susie__susie_grid_results.csv")))

    def l664():
        m = [r["relative_path"] for r in _csv_rows(dep / "MANIFEST.tsv", "\t")]
        a = sum(x.startswith("scripts/scripts/") for x in m)
        b = sum(x.startswith("scripts/R/") for x in m)
        return [a + b, a, b]
    put("L664", l664)

    def l665():
        n = 0
        for _, ds, fs in os.walk(v3):
            ds[:] = [d for d in ds if d not in (".git", "__pycache__")]
            n += sum(f.endswith(".i8") for f in fs)
        return n
    put("L665", l665)

    def l666():
        rep = (dep / "REPORT.md").read_text(encoding="utf-8", errors="replace").splitlines()
        words = {w: i for i, w in enumerate("one two three four five six seven eight nine".split(), 1)}
        l6 = [l for l in rep if l.startswith("6. **Three Part B facts")]
        m6 = re.search(r"\(a\) (\w+) writer scripts were edited after their output was written", l6[0])
        l1 = [l for l in rep if l.startswith("1. **Writer scripts whose current bytes are not the version that produced the delivered file**")]
        scripts = sorted(set(re.findall(r"`([\w./-]+\.(?:R|py))`", l1[0].split("For A9")[0])))
        relab = sorted({re.match(r"\| (A\d+) \|", l).group(1) for l in rep if re.match(r"\| A\d+ \|", l) and "hand relabel" in l and "no writer script" in l})
        if len(scripts) != words[m6.group(1).lower()]:
            raise ValueError("scripts listed (%d) != count spelled (%s)" % (len(scripts), m6.group(1)))
        return [len(scripts), len(relab)]
    put("L666", l666)

    def l667():
        fig = json.loads((v3 / "R2_rewrite" / "figures" / "af3" / "P31R_ACTIVE_FIGURES.json").read_text(encoding="utf-8"))["figures"]
        return [sum(int(f["pages"]) for f in fig), len(fig)]
    put("L667", l667)
    put("L668", lambda: len(_csv_rows(v3 / "R2_rewrite" / "review_R2" / "seen" / "wpb" / "proxy_ledger_draft.csv")))
    put("L669", lambda: sum(1 for l in (v3 / "R1_stage1" / "S1R1_03_decoy_method" / "S1R1_decoy_rates.log").read_text(encoding="utf-8", errors="replace").splitlines()
                            if re.match(r"^md5 (\S+)\s+([0-9a-f]{32})\s+(\S+)$", l)))

    def l670():
        ws = {}
        for r in read_tsv(release / "inputs" / "RELEASE_INPUTS.tsv"):
            if r["role"] == "copy" and r["ship_status"] == "ship" and r["release_path"].startswith("workstreams/"):
                f = r["release_path"].split("/")[1]
                ws[f] = ws.get(f, 0) + 1
        return [sum(ws.values())] + list(ws.values())
    put("L670", l670)
    put("L671", lambda: len(_body_rows(v3 / "R2_rewrite" / "tables" / "P31R_Table2.md")))
    put("L673", lambda: len(_body_rows(v3 / "R2_rewrite" / "tables" / "P31R_Table1.md")))

    def l672():
        bi = _csv_rows(v3 / "R1_stage1" / "S1R1_R2_10_integration" / "branches_integrated.tsv", "\t")
        bs = {}
        for r in bi:
            bs[r["source"]] = bs.get(r["source"], 0) + 1
        return [len(bi)] + list(bs.values())
    put("L672", l672)

    def l674():
        cl = [l for l in (dep / "scripts" / "reports" / "Paper3_Program" / "run_logs" / "REGEN_20260817" / "01_mr_screen.log").read_text(encoding="utf-8", errors="replace").splitlines()
              if l.startswith("MR cells:")]
        m = re.fullmatch(r"MR cells: (\d+) with-result, (\d+) attempted-no-result, (\d+) no-outcome-data, (\d+) null-exposure-gene \| full-grid (\d+)/(\d+) COMPLETE", cl[0])
        return [int(m.group(3)), int(m.group(2))]
    put("L674", l674)

    def l675():
        ins = {r["gene"] for r in _csv_rows(tfr / "phase1_instruments__ukb_cis_instruments.csv")}
        genes = {r["gene"] for r in mr()}
        return [len(ins - genes), len(ins), len(genes)]
    put("L675", l675)

    def l533():
        pv = read_tsv(release / "inputs" / "PROVIDERS.tsv")
        return [len(pv), sum(r["in_declarations"] == "yes" for r in pv)]
    put("L533", l533)

    def l676():
        fz = v3 / "R1_stage1" / "S1R1_06_claims_ledger" / "frozen_inputs" / "PROVIDERS_at_1c_L676.tsv"
        if md5_file(fz) != L676_FROZEN_MD5:
            raise ValueError("the frozen 1c provider list is not the md5-pinned file (%s)" % md5_file(fz)[:8])
        pv = read_tsv(fz)
        return [len(pv), sum(r["in_declarations"] == "yes" for r in pv)]
    put("L676", l676)
    return res


def r10_ledger(R, ledger, v3, release):
    rc = r10_recount(v3, release)
    nok = nexc = 0
    for i in R10_IDS:
        if i not in ledger:
            R.fail("R10", "%s is not a ledger row" % i)
            continue
        val, err = rc[i]
        if err:
            R.fail("R10", "%s: recount failed (%s)" % (i, err))
            continue
        try:
            led = json.loads(ledger[i]["value"])
        except ValueError:
            R.fail("R10", "%s: ledger value is not JSON: %r" % (i, ledger[i]["value"][:40]))
            continue
        ex = R10_PENDING.get(i)
        if ex:
            if led == val:
                R.fail("R10", "%s: declared pending exception (%s) is STALE: the recount %s now equals the ledger" % (i, ex["reason"], val))
            elif led != ex["ledger"] or val != ex["recount"]:
                R.fail("R10", "%s: ledger %s / recount %s differ from the declared exception (ledger %s, recount %s; %s): drift" % (i, led, val, ex["ledger"], ex["recount"], ex["reason"]))
            else:
                nexc += 1
        elif led != val:
            R.fail("R10", "%s: ledger %s != recount %s now (%s)" % (i, led, val, ledger[i]["claim"][:60]))
        else:
            nok += 1
    for i in set(R10_PENDING) - set(R10_IDS):
        R.fail("R10", "pending exception %s is not a recounted row" % i)
    R.note("R10", "%d ledger rows (L660-L675, L533 from the live provider list, L676 from its frozen input) recounted from their inputs now: %d equal, %d declared pending exception(s)" % (len(R10_IDS), nok, nexc))


# ---- R15 (S6.3b2a, S2-D37 D4): provenance labels, portable-copy transform, third-party notices, and the provider count in the documents. Independent of the builder: own reading.
PROV_LABELS = ("full run-time hash match", "partial launch-record match", "timestamp evidence only")
DEPOSIT_PREFIXES = ("paper3_deposit_request_v4/", "paper3_deposit_request_v5/")
THIRD_PARTY_RX = re.compile("|".join([r"(?i:copy" + r"right)\s*(?:\(c\)|©)", r"\(c\)\s*(?:19|20)\d\d", "(?i:all rights" + " reserved)", "SPDX-License" + "-Identifier",      # pieces: this file is shipped and scanned by this very check
    "GNU (?:Lesser |Affero )?General Public" + " License", "Apache" + " License", "Permission is hereby" + " granted", "Redistribution and use in" + " source and binary forms"]))


def r15_provenance(R, release, v3, rel, text_of, man):
    ci = {r["ship_path"]: r for r in read_tsv(Path(release) / "inputs" / "CODE_INCLUDE.tsv") if r["ship"] == "1"}
    mrow = {r["path"]: r for r in man}
    n = {lab: 0 for lab in PROV_LABELS}
    ntr = 0
    for sp, c in sorted(ci.items()):
        m = mrow.get(sp)
        if m is None:
            R.fail("R15", "%s: shipped code row has no manifest row" % sp)
            continue
        lab, ev, evid, cr = c.get("provenance_label", ""), c.get("evidence_md5", ""), c.get("provenance_evidence", ""), c.get("code_rights", "")
        if m.get("provenance_label", "") != lab:
            R.fail("R15", "%s: manifest provenance_label %r differs from CODE_INCLUDE %r" % (sp, m.get("provenance_label", ""), lab))
        orig = Path(v3) / c["path"]
        have = md5_file(orig) if orig.is_file() else ""
        if m["source_md5"] != have:
            R.fail("R15", "%s: manifest source_md5 %s is not the md5 of the original %s (%s)" % (sp, m["source_md5"][:8], c["path"], have[:8]))
        tr = m.get("transform", "")
        same = m["md5"] == m["source_md5"]
        if same != (tr == "byte-identical to the original") or (not same and not tr.startswith("derived from the original")):
            R.fail("R15", "%s: transform %r does not match md5 %s vs source_md5 %s (a changed copy must say it is derived, not byte-identical)" % (sp, tr[:50], m["md5"][:8], m["source_md5"][:8]))
        if c["path"].startswith(DEPOSIT_PREFIXES) and not lab:
            R.fail("R15", "%s: a shipped deposit-package script carries no provenance label" % sp)
        if lab:
            if lab not in PROV_LABELS:
                R.fail("R15", "%s: label %r is not one of the three" % (sp, lab))
                continue
            n[lab] += 1
            if not evid or not (cr == "authors" or cr.startswith("third_party:")):
                R.fail("R15", "%s: a label needs stated evidence and code_rights (authors or third_party:<what>)" % sp)
            if lab == PROV_LABELS[0] and not (len(ev) == 32 and ev == have):
                R.fail("R15", "%s: %r but the recorded md5 %r is not the original's md5" % (sp, lab, ev))
            if lab == PROV_LABELS[1] and not (len(ev) == 8 and have.startswith(ev)):
                R.fail("R15", "%s: %r but the logged prefix %r does not start the original's md5" % (sp, lab, ev))
            if lab == PROV_LABELS[2] and ev:
                R.fail("R15", "%s: %r claims no hash but evidence_md5 is set" % (sp, lab))
        if not cr.startswith("third_party:") and sp in text_of and text_of[sp] is not None:
            for i, ln in enumerate(text_of[sp].splitlines(), 1):
                if THIRD_PARTY_RX.search(ln):
                    ntr += 1
                    R.fail("R15", "%s:%d  a third-party copyright or licence notice in a file whose code_rights is not third_party:<what> (flag it, do not relabel it): %s" % (sp, i, ln.strip()[:90]))
                    break
    for k, m in mrow.items():
        if m.get("provenance_label", "") and k not in ci:
            R.fail("R15", "%s: manifest row carries a provenance label but is not a shipped CODE_INCLUDE row" % k)
        if m.get("transform", "") and k not in ci:
            R.fail("R15", "%s: transform set on a row that is not a shipped CODE_INCLUDE row" % k)
    rd = text_of.get("README.md") or ""
    for lab in PROV_LABELS:
        if "`%s`" % lab not in rd:
            R.fail("R15", "README.md does not define the label %r" % lab)
    for doc in ("README.md", "REPRODUCE.md"):
        t = text_of.get(doc) or ""
        if "[L676]" in t:
            R.fail("R15", "%s cites the frozen row L676; the provider count is the live row L533" % doc)
    if "[L533]" not in rd:
        R.fail("R15", "README.md does not print the provider count with its ledger row L533")
    R.note("R15", "%d labelled code rows (%s); manifest labels = CODE_INCLUDE; every code row's transform agrees with its md5 against the original; %d third-party notice(s) in unflagged code; documents cite L533, never L676" % (
        sum(n.values()), ", ".join("%s %d" % (k, v) for k, v in n.items()), ntr))


def _rule_rows(release):
    rules = Path(release) / "inputs" / "TIER_RULES.tsv"
    return read_tsv(rules) if rules.exists() else []


def r11_tiers(R, release, rel, text_of, man):
    """R11: build/tiers_by_file.tsv (file, tier, providers, decision, reason) covers every shipped file, and nothing ships against its rule.
    Rules are read from inputs/TIER_RULES.tsv (rule_id kind target tier providers handling selector basis reason): switch rows SW*, file rows (handling
    columns:a|b | rows:<col>~<regex>:<value cols>[@resource] | binned:<col>:<step> | derived | declared:...).  The gate judges CONTENT, not labels:
      * a file whose header has a variant-id column and a value column is Tier 2 unless a file rule declares it Tier 1 (derived); with no rule it fails (UNRULED);
      * Tier 3 never ships while switch SW1 is 'withheld' (default); Tier 2 with decision 'ship' needs every provider confirmed CC0 / CC BY;
      * reduced_* / ship_reduced_* / binned_* decisions need a file rule and the rule is checked on the shipped file (withheld columns absent, binned values on the grid,
        masked row values present); S2-D34(6): eur_* / het_p / strict columns of east_asian/ files ship only as the file rule declares: selector comparator:<accession> (the accession
        rule says cc0) or derived:<columns> (Tier 1 principle SINGLE_INSTRUMENT_WALD); with no declaration the S2-D33(1) test applies (every provider on the manifest licence_class open)."""
    release = Path(release)
    tp = release / "build" / "tiers_by_file.tsv"
    rules = _rule_rows(release)
    if not rules:
        R.fail("R11", "inputs/TIER_RULES.tsv is absent or empty")
    if not tp.exists():
        R.fail("R11", "build/tiers_by_file.tsv is absent")
        return
    hdr = tp.read_text(encoding="utf-8").splitlines()[0].split("\t")
    if hdr != TIERS_COLS:
        R.fail("R11", "tiers_by_file.tsv columns differ from %s: %s" % (TIERS_COLS, hdr))
        return
    rows = read_tsv(tp)
    ok_prov = cc_ok_providers(release)
    known = {r["provider_id"].strip() for r in read_tsv(release / "inputs" / "PROVIDER_TERMS.tsv")}
    sw = {r["target"]: r for r in rules if r.get("kind") == "switch"}
    frules = collections.defaultdict(list)
    for r in rules:
        if r.get("kind") == "file":
            frules[r["target"].replace("\\", "/")].append(r)
    if "TPMI_DERIVED" not in sw or sw["TPMI_DERIVED"]["handling"].strip() not in ("withheld", "ship"):
        R.fail("R11", "the TPMI switch (kind switch, target TPMI_DERIVED, handling withheld|ship) is absent from TIER_RULES.tsv")
    tpmi_ship = sw.get("TPMI_DERIVED", {}).get("handling", "").strip() == "ship"
    acc_rules = {r["target"].strip(): r for r in rules if r.get("kind") == "accession"}
    princ = {r["target"].strip(): r for r in rules if r.get("kind") == "principle"}
    by = collections.defaultdict(list)
    for r in rows:
        by[re.sub(r"^out/", "", r["file"].strip().replace("\\", "/"))].append(r)
    mlc = {r["path"]: r.get("licence_class", "") for r in man}
    for k, v in by.items():
        if len(v) > 1:
            R.fail("R11", "%s has %d rows in tiers_by_file.tsv" % (k, len(v)))
    n2 = n3 = nder = nred = 0
    for k in sorted(rel):
        if k in SELF_EXEMPT:
            continue
        if k not in by:
            R.fail("R11", "shipped file without a tiers_by_file.tsv row: %s" % k)
            continue
        r = by[k][0]
        tier = re.sub(r"(?i)^\s*(tier|t)?\s*", "", r["tier"]).strip()
        dec = r["decision"].strip()
        provs = split_ids(r["providers"])
        for pv in provs:
            if pv not in known:
                R.fail("R11", "%s: provider %r is not in PROVIDER_TERMS.tsv" % (k, pv))
        if tier not in ("1", "2", "3"):
            R.fail("R11", "%s: tier %r is not 1, 2 or 3" % (k, r["tier"]))
            continue
        if TIER_DECISION_WITHHOLD.match(dec):
            R.fail("R11", "%s ships in out/ but tiers_by_file.tsv says %s" % (k, dec))
            continue
        if not re.match(r"(?i)^(ship|reduced|binned)\w*$", dec):
            R.fail("R11", "%s: decision %r is not ship* / reduced* / binned* / withheld*" % (k, dec))
            continue
        if not r["reason"].strip():
            R.fail("R11", "%s: empty reason" % k)
        t = text_of.get(k)
        hdrs, drows = [], []
        if k.endswith((".csv", ".tsv")) and t:
            delim = "\t" if k.endswith(".tsv") else ","
            hdrs = [h.strip().strip('"') for h in t.splitlines()[0].split(delim)]
            drows = list(csv.DictReader(io.StringIO(t, newline=""), delimiter=delim))
        pervar = bool(hdrs) and is_per_variant(hdrs)
        fr = frules.get(k, [])
        if len(fr) > 1:
            R.fail("R11", "%s has %d file rules in TIER_RULES.tsv" % (k, len(fr)))
        rule = fr[0] if fr else None
        reduced = bool(re.match(r"(?i)^(reduced|binned|ship_reduced)", dec))
        if tier == "3":
            n3 += 1
            if not tpmi_ship:
                R.fail("R11", "%s is Tier 3 (no licence) and ships (decision %s) while the TPMI switch is 'withheld'" % (k, dec))
        if rule is not None and rule["tier"].strip() and rule["tier"].strip() != tier:
            R.fail("R11", "%s: tiers_by_file says Tier %s, rule %s says Tier %s" % (k, tier, rule["rule_id"], rule["tier"]))
        if reduced and rule is None:
            R.fail("R11", "%s: decision %s but no file rule in TIER_RULES.tsv" % (k, dec))
        if tier == "2" and not reduced and tier != "3":
            n2 += 1
            declared_ship = rule is not None and rule["handling"].strip() == "derived" and "S2-D35" in rule["basis"]     # RC2: a provider value the manuscript prints at the analysed variant ships (S2-D35), declared by file rule
            bad = [] if declared_ship else [pv for pv in provs if pv not in ok_prov]
            if not provs:
                R.fail("R11", "%s: Tier 2 with no providers listed" % k)
            if bad:
                R.fail("R11", "%s: Tier 2 ships whole (decision %s) but provider(s) %s are not confirmed CC0 / CC BY in PROVIDER_TERMS.tsv" % (k, dec, bad))
        if pervar and not (rule is not None or (tier == "2" and provs and all(p in ok_prov for p in provs))):
            R.fail("R11", "%s: UNRULED: header has a variant id (%s) and a value column (%s) but no TIER_RULES file row covers it" % (
                k, [h for h in hdrs if VAR_ID.match(h)][:2], [h for h in hdrs if VAR_VAL.search(h) and not VAR_ID.match(h)][:2]))
        if pervar and rule is not None and tier == "1":
            nder += 1
        # the rule, checked on the shipped file
        if rule is not None:
            h = rule["handling"].strip()
            if h.startswith("columns:"):
                nred += 1
                cols = h[len("columns:"):].split("|")
                still = [c for c in cols if c in hdrs]
                if still:
                    R.fail("R11", "%s: rule %s withholds column(s) %s, which ship" % (k, rule["rule_id"], still))
            elif h.startswith("binned:"):
                nred += 1
                _, col, step = h.split(":")
                step = float(step)
                if col not in hdrs:
                    R.fail("R11", "%s: rule %s bins column %r, which is absent" % (k, rule["rule_id"], col))
                else:
                    off = 0
                    for rr in drows:
                        try:
                            v = float(rr[col])
                        except ValueError:
                            continue
                        if abs(v / step - round(v / step)) > 1e-6:
                            off += 1
                    if off:
                        R.fail("R11", "%s: %d value(s) of %s are not on the %s grid (rule %s)" % (k, off, col, step, rule["rule_id"]))
            elif h.startswith("rows:"):
                nred += 1
                body = h[len("rows:"):]
                m = re.match(r"^(\w+)~(.*):([\w|]+)(?:@\w+)?$", body)
                if not m:
                    R.fail("R11", "%s: rule %s has an unreadable rows handling %r" % (k, rule["rule_id"], h[:60]))
                else:
                    sel, rx, vcols = m.group(1), re.compile(m.group(2)), [c for c in m.group(3).split("|") if c in hdrs]
                    hit = [rr for rr in drows if rx.search(rr.get(sel, ""))]
                    masked = [rr for rr in hit if any(rr[c].strip() == MASK_TOKEN for c in vcols)]
                    if hit and not masked:
                        R.fail("R11", "%s: rule %s selects %d row(s) and none is masked (%s)" % (k, rule["rule_id"], len(hit), MASK_TOKEN))
                    # S2-D41: a rows rule whose declared providers are ALL unconfirmed masks every selected row (the provider of a row cannot open it); a mixed list keeps the any-row test above
                    rprov = split_ids(rule.get("providers", ""))
                    if rprov and not any(pv in ok_prov for pv in rprov):
                        unm = [rr for rr in hit if any(rr[c].strip() not in (MASK_TOKEN, "", "NA") for c in vcols)]
                        if unm:
                            R.fail("R11", "%s: rule %s declares only unconfirmed provider(s) %s and %d of its %d selected row(s) still carry a value in %s" % (k, rule["rule_id"], rprov, len(unm), len(hit), vcols))
        if k.startswith("east_asian/") and any(EAS_EUR_COL.match(x) for x in hdrs):
            eur_cols = [x for x in hdrs if EAS_EUR_COL.match(x)]
            sel = (rule.get("selector") or "").strip() if rule is not None else ""
            if sel.startswith("comparator:"):
                # S2-D34(6): the eur_* / het_p / strict columns follow the COMPARATOR, declared in the file rule as accession(s) whose accession rule says cc0
                accs = split_ids(sel[len("comparator:"):])
                bad = [a for a in accs if (acc_rules.get(a) or {}).get("handling", "").strip() != "cc0"]
                if not accs or bad:
                    R.fail("R11", "%s ships eur_* / het_p / strict columns and rule %s declares comparator(s) %s, which no accession rule in TIER_RULES.tsv marks cc0" % (k, rule["rule_id"], bad or "(none)"))
            elif sel.startswith("derived:"):
                # S2-D34(2),(6): single-instrument Wald results are Tier 1 (principle P01); only the declared derived columns may ship
                allowed = set(split_ids(sel[len("derived:"):]))
                extra_c = [x for x in eur_cols if x not in allowed]
                if "SINGLE_INSTRUMENT_WALD" not in princ or princ["SINGLE_INSTRUMENT_WALD"].get("tier", "").strip() != "1":
                    R.fail("R11", "%s declares derived eur_* columns but TIER_RULES.tsv has no Tier 1 principle SINGLE_INSTRUMENT_WALD" % k)
                if extra_c:
                    R.fail("R11", "%s ships comparator column(s) %s that rule %s does not declare derived (permitted: %s)" % (k, extra_c, rule["rule_id"], sorted(allowed)))
            else:
                mp = split_ids(mlc.get(k, "").partition(":")[2])
                bad = [pv for pv in mp if pv not in ok_prov]
                if bad or not mp:
                    R.fail("R11", "%s ships eur_* / het_p / strict columns (S2-D33(1)) while comparator provider(s) %s are not confirmed CC0 / CC BY" % (k, bad or "(none listed)"))
    for k in sorted(set(by) - set(rel)):
        if not TIER_DECISION_WITHHOLD.match(by[k][0]["decision"].strip()):
            R.fail("R11", "tiers_by_file.tsv row for a file that is neither in out/ nor declared withheld: %s (decision %s)" % (k, by[k][0]["decision"]))
    R.note("R11", "%d tier rows cover %d shipped files (%d Tier 2 whole, %d with a rule checked on the file, %d declared Tier 1 derived despite a variant-id header, %d Tier 3 shipped); TPMI switch: %s; CC0/CC BY providers: %d" % (
        len(rows), len(rel) - len(SELF_EXEMPT & set(rel)), n2, nred, nder, n3, "ship" if tpmi_ship else "withheld", len(ok_prov)))


def r8_conservative(R, v3, rel, text_of, man, ledger):
    """R8, second half: a conservative-prior column ships only (a) in a declared file and column, and (b) on rows whose (gene, outcome) is one of the
    canonical 81 pairs (read from the shipped tiered table; its pair count must equal L010 when the ledger is given, and its source file's pair set when reachable);
    Table 2 may additionally carry the declared pairs and columns of CONS_DECLARED (S2-D34(5))."""
    canon_txt = text_of.get(CANON_TIERED)
    if not canon_txt:
        R.fail("R8", "canonical tiered table %s absent" % CANON_TIERED)
        return 0
    pairs, n = canonical_pairs(v3, canon_txt, man)
    if len(pairs) != n:
        R.fail("R8", "canonical tiered table has duplicate (gene, outcome) rows: %d rows, %d pairs" % (n, len(pairs)))
    if len(pairs) != int(json.loads(ledger["L010"]["value"])):
        R.fail("R8", "the canonical tiered table has %d pairs, L010 is %s" % (len(pairs), ledger["L010"]["value"]))
    src = next((r["source_path"] for r in man if r["path"] == CANON_TIERED), "")
    sp = Path(v3) / src.split(";")[0] if src else None
    if sp is not None and sp.is_file():
        spairs = {(r["gene"].strip(), r["outcome"].strip()) for r in _csv_rows(sp)}
        if spairs != pairs:
            R.fail("R8", "the shipped tiered table's pairs differ from its source file (%d vs %d pairs)" % (len(pairs), len(spairs)))
    t2md = table2_md_pairs(v3)
    for fk, decl in CONS_DECLARED.items():
        for pr, (cols_ok, why) in decl.items():
            if not why.strip():
                R.fail("R8", "CONS_DECLARED %s %s has no reason" % (fk, pr))
            if t2md is not None and pr not in t2md:
                R.fail("R8", "CONS_DECLARED %s names %s, which the manuscript's Table 2 does not print" % (fk, pr))
    nscan = 0
    for k, t in text_of.items():
        if not (k.endswith((".csv", ".tsv")) and t):
            continue
        delim = "\t" if k.endswith(".tsv") else ","
        hdr = [h.strip('"') for h in t.splitlines()[0].split(delim)]
        cons = [h for h in hdr if CONS_COL.search(h)]
        nscan += 1
        if not cons:
            continue
        bad_cols = [h for h in cons if h not in CONS_PERMITTED.get(k, set())]
        if bad_cols:
            R.fail("R8", "conservative-prior column(s) %s shipped in %s, which is not a declared carrier (decision 6; permitted: %s)" % (bad_cols, k, sorted(CONS_PERMITTED) or "none"))
            continue
        if "gene" not in hdr or "outcome" not in hdr:
            R.fail("R8", "%s carries a conservative-prior column but has no gene / outcome columns to identify its pairs" % k)
            continue
        decl = CONS_DECLARED.get(k, {})
        for i, r in enumerate(csv.DictReader(io.StringIO(t, newline=""), delimiter=delim), 1):
            pr = (r["gene"].strip(), r["outcome"].strip())
            if pr in pairs:
                continue
            if pr not in decl:
                R.fail("R8", "%s row %d: (%s, %s) is not one of the canonical %d pairs%s, yet the file carries a conservative-prior column" % (
                    k, i, r["gene"], r["outcome"], len(pairs), " or a declared pair (CONS_DECLARED, S2-D34(5))" if decl else ""))
                continue
            held = [c for c in cons if (r.get(c) or "").strip() not in CONS_NO_VALUE + (MASK_TOKEN,)]
            undeclared = [c for c in held if c not in decl[pr][0]]
            md_p = Path(v3) / TABLE2_MD
            if md_p.is_file():      # RC1b: a value the declaration permits must be one that the CURRENT Table 2 text prints (an exact numeral, not a substring of a longer one)
                md_t = md_p.read_text(encoding="utf-8")
                for c in held:
                    if c in decl[pr][0] and not re.search(r"(?<![\d.])%s(?![\d])" % re.escape(r[c].strip()), md_t):
                        R.fail("R8", "%s row %d: declared pair (%s, %s) holds %s = %s, which the current Table 2 text does not print" % (k, i, pr[0], pr[1], c, r[c]))
            if undeclared:
                R.fail("R8", "%s row %d: declared pair (%s, %s) holds a conservative-prior value in %s; the declaration permits %s" % (
                    k, i, pr[0], pr[1], undeclared, sorted(decl[pr][0]) or "none"))
    return nscan


# ---------------------------------------------------------------- the gate
# ---------------------------------------------------------------- R14: withheld values must not reappear in out/ (RC1b fix Q01)
LEAK_FL = re.compile(r"(?<![\w.])(?:\d+\.\d+|\.\d+|\d+(?:\.\d+)?[eE][-+]?\d+)(?![\w])")
LEAK_MIN_SIG = 4
LEAK_ALLOW = "inputs/LEAKSCAN_COINCIDENCES.tsv"
LEAK_COLS = ["out_path", "value", "withheld_source", "reason"]


def leak_sig(tok):
    m = re.split(r"[eE]", tok)[0]
    if "." in m:
        m = m.rstrip("0")
    return len(m.replace(".", "").lstrip("0"))


def leak_key(tok):
    return "%.10g" % float(tok)


def _grid(p):
    t = Path(p).read_text(encoding="utf-8", errors="replace")
    return list(csv.reader(io.StringIO(t, newline=""), delimiter="\t" if str(p).endswith(".tsv") else ","))


def leak_values(release, ledger):
    """{canonical key: (set of withheld sources, a printed token)} for every withheld cell value of >= LEAK_MIN_SIG significant digits.
    Tabular withheld files only (the .md / .json / .log files re-quote values that sit in a withheld table): a cell counts when the shipped counterpart
    (same relative path) has no such column, a different value at the same row and column, or a different row count; no counterpart = every cell."""
    release = Path(release)
    out, wh = release / "out", release / "withheld"
    vals = {}

    def add(text, src, kept=frozenset()):
        for t in LEAK_FL.findall(text):
            if leak_sig(t) >= LEAK_MIN_SIG and leak_key(t) not in kept:
                vals.setdefault(leak_key(t), (set(), t))[0].add(src)
    for w in sorted(wh.rglob("*")):
        if not w.is_file() or w.suffix.lower() not in (".csv", ".tsv"):
            continue
        k = w.relative_to(wh).as_posix()
        W = _grid(w)
        O = _grid(out / k) if (out / k).is_file() else None
        for i, row in enumerate(W[1:], 1):
            for j, c in enumerate(row):
                kept = frozenset()
                if O is not None and len(O) == len(W) and j < len(W[0]) and W[0][j] in O[0] and i < len(O):
                    jj = O[0].index(W[0][j])
                    if jj < len(O[i]) and O[i][jj] == c:
                        continue
                    if W[0][j] in R17_NOTE and jj < len(O[i]) and O[i][jj] in R17_MARKS:
                        continue        # S2-D40 (W2): the free text of a masked row's notes cell is not a table value (R17 guards it); the row's own value / se / p cells are counted
                    if jj < len(O[i]):
                        kept = frozenset(leak_key(x_) for x_ in LEAK_FL.findall(O[i][jj]))      # S2-D41: a token that still stands in the SHIPPED cell (the conditional z_cond beside a redacted marginal beta) is shipped, not withheld
                add(c, k, kept)
    pub = out / "ledger" / "claims_ledger_public.tsv"
    if pub.is_file():
        for r in read_tsv(pub):
            if '"masked":true' in r.get("value", "") and r["id"] in ledger:
                add(ledger[r["id"]].get("value", "") + " " + ledger[r["id"]].get("display", ""), "ledger:" + r["id"])
    return vals


LEAK_SHIP_MIN_SIG = 6        # RC2 (S2-D35(3)): a shipped token of >= 6 significant digits is compared with the withheld values at 6 digits, so a value re-printed at another precision is one value
LEAK_ROW_MIN_SIG = 4         # RC2: inside ONE row (a withheld cell and the other cells of the same shipped row, a masked ledger row and its other fields) the comparison runs at the precision of the shorter token, floor 4 digits


def leak_round(tok, n):
    return float("%.*g" % (n, float(tok)))


def leak_precision_hit(tok, idx):
    """A shipped token of n >= LEAK_SHIP_MIN_SIG significant digits is a re-print of a withheld value of >= n digits when the withheld value, rounded to n digits, is the token (idx[n] = those rounded values)."""
    n = leak_sig(tok)
    return LEAK_SHIP_MIN_SIG <= n <= 15 and leak_round(tok, n) in idx.get(n, ())


def leak_declared(allow, a_by_file, k, tok):
    """The (out_path, key) of the declared coincidence that covers this token, or None: the same value (10 digits), or a declared value of >= n digits that rounds to the token at its n >= LEAK_SHIP_MIN_SIG digits."""
    key = leak_key(tok)
    if (k, key) in allow:
        return (k, key)
    n = leak_sig(tok)
    if LEAK_SHIP_MIN_SIG <= n <= 15:
        for ak, v in a_by_file.get(k, ()):
            if leak_sig(v) >= n and leak_round(v, n) == leak_round(tok, n):
                return ak
    return None


def leak_rowlocal(release, ledger):
    """[(shipped path, row, column, token, withheld token, withheld column or source)] where a shipped cell prints, at the precision of the shorter token (floor LEAK_ROW_MIN_SIG), a withheld value of the SAME row:
    (a) a cell of withheld/<f> that differs from out/<f> (or sits in a column out/<f> lacks; same row count) against every cell of that shipped row; (b) a masked row of the public ledger (the private value and display)
    against the other fields of its public row. Notes columns and claim text are cells like any other."""
    release = Path(release)
    out, wh = release / "out", release / "withheld"
    hits = []

    def tokens(s):
        return [(t, leak_sig(t)) for t in LEAK_FL.findall(s or "") if leak_sig(t) >= LEAK_ROW_MIN_SIG]

    def cmp(wtoks, cells, where, skip=None):
        for wt, wn in wtoks:
            for cname, cell in cells:
                if cname == skip:
                    continue      # the shipped counterpart of the withheld cell itself (masked, binned or rounded on purpose)
                for ot, on in tokens(cell):
                    m = min(wn, on)
                    if leak_round(wt, m) == leak_round(ot, m):
                        hits.append((where[0], where[1], cname, ot, wt, where[2]))
    for w in sorted(wh.rglob("*")):
        if not w.is_file() or w.suffix.lower() not in (".csv", ".tsv"):
            continue
        k = w.relative_to(wh).as_posix()
        if not (out / k).is_file():
            continue
        W, O = _grid(w), _grid(out / k)
        if len(W) != len(O):
            continue
        for i in range(1, len(W)):
            ocells = [(O[0][j], c) for j, c in enumerate(O[i]) if j < len(O[0])]
            for j, c in enumerate(W[i]):
                if j >= len(W[0]):
                    continue
                if W[0][j] in O[0]:
                    jj = O[0].index(W[0][j])
                    if jj < len(O[i]) and O[i][jj] == c:
                        continue
                cmp(tokens(c), ocells, (k, i, W[0][j]), skip=W[0][j])
    pub = out / "ledger" / "claims_ledger_public.tsv"
    if pub.is_file():
        for r in read_tsv(pub):
            if r["id"] in ledger and ('"masked":true' in r.get("value", "")):
                priv = ledger[r["id"]]
                cmp(tokens(priv.get("value", "") + " " + priv.get("display", "")), [(c, r[c]) for c in r if c not in ("value", "display")],
                    ("ledger/claims_ledger_public.tsv", r["id"], "private value"))
    return hits


def r14_leakscan(R, release, rel, text_of, ledger):
    release = Path(release)
    allow_p = release / LEAK_ALLOW
    if not allow_p.is_file():
        R.fail("R14", "%s is absent (the declared coincidences)" % LEAK_ALLOW)
        return
    hdr = allow_p.read_text(encoding="utf-8").splitlines()[0].split("\t")
    if hdr != LEAK_COLS:
        R.fail("R14", "%s columns differ: %s" % (LEAK_ALLOW, hdr))
        return
    allow = {}
    for r in read_tsv(allow_p):
        if not r["reason"].strip():
            R.fail("R14", "declared coincidence without a reason: %s %s" % (r["out_path"], r["value"]))
        allow[(r["out_path"], leak_key(r["value"]))] = r
    vals = leak_values(release, ledger)
    a_by_file = {}
    for (ak, akey), r in allow.items():
        a_by_file.setdefault(ak, []).append(((ak, akey), r["value"]))
    idx = {}
    for (srcs, tok) in vals.values():
        for n in range(LEAK_SHIP_MIN_SIG, min(leak_sig(tok), 15) + 1):
            idx.setdefault(n, {}).setdefault(leak_round(tok, n), (srcs, tok))
    seen, nscan, bad = set(), 0, 0
    for k, t in rel.items():
        txt = text_of.get(k)
        if txt is None:
            continue
        nscan += 1
        for tok in LEAK_FL.findall(txt):
            n = leak_sig(tok)
            if n < LEAK_MIN_SIG:
                continue
            key = leak_key(tok)
            hit = vals.get(key) or (idx.get(n, {}).get(leak_round(tok, n)) if LEAK_SHIP_MIN_SIG <= n <= 15 else None)
            if hit is None:
                continue
            dk = leak_declared(allow, a_by_file, k, tok)
            if dk is not None:
                seen.add(dk)
            else:
                bad += 1
                R.fail("R14", "%s prints %s, a withheld value of %s (compared at 10 digits, and at the token's own precision for a token of >= %d digits), and it is not a declared coincidence" % (
                    k, tok, "; ".join(sorted(hit[0])), LEAK_SHIP_MIN_SIG))
    rl = leak_rowlocal(release, ledger)
    nrl = 0
    for (k, i, cname, ot, wt, wsrc) in rl:
        dk = leak_declared(allow, a_by_file, k, ot)
        if dk is not None:
            seen.add(dk)
            continue
        nrl += 1
        R.fail("R14", "%s row %s, column %s prints %s, which agrees at %d or more digits with the withheld value %s of the same row (%s), and it is not a declared coincidence" % (
            k, i, cname, ot, LEAK_ROW_MIN_SIG, wt, wsrc))
    bad += nrl
    seen_n = len(seen)
    pend = sorted(k_ for k_ in seen if (k_ in allow) and allow[k_]["reason"].startswith(("quoted_value", "REAL-PENDING")))
    if pend:
        R.warns.append("R14 %d declared occurrence(s) are marked quoted_value (a provider per-variant value quoted in shipped code or a free-text cell; declared in README and REPRODUCE, item (d)): %s" % (
            len(pend), "; ".join("%s %s" % x for x in pend[:4]) + (" ..." if len(pend) > 4 else "")))
    stale = sorted(set(allow) - seen)
    for k, key in stale[:20]:
        R.warns.append("R14 declared coincidence no longer occurs: %s %s" % (k, key))
    R.note("R14", "%d withheld values (>= %d significant digits) searched in %d shipped text files (every field, notes columns and the public ledger included; 10 digits, and at the token's own precision for tokens of >= %d digits) and, row by row, at >= %d digits against the other cells of the same row (%d row-local); %d declared coincidences used, %d undeclared; %d declared rows stale" % (
        len(vals), LEAK_MIN_SIG, nscan, LEAK_SHIP_MIN_SIG, LEAK_ROW_MIN_SIG, len(rl), seen_n, bad, len(stale)))


# ---------------------------------------------------------------- R16: a shipped column must not be a function of withheld columns (S2-D40, W1)
R16_BETA = re.compile(r"(?i)(^|_)(beta|b|bx|by|logor|logodds)(_|$)")
R16_SE = re.compile(r"(?i)(^|_)(se|sebx|seby|stderr)(_|$)")
R16_REL_TOL = 1e-4
R16_MIN_ROWS = 3
R16_F = re.compile(r"(?i)(^|_)F(_stat)?(_|$)")      # S2-D41: F, F_stat, lead_F, decode_F_stat
# Declared derived columns that ship although they are a function of withheld columns (file, column) -> reason.  A single-instrument Wald ratio has the z of the outcome association
# (its first-order SE is se_out / |b_exp|), so its p IS the outcome p at the instrument; principle P01 (S2-D34(2), USER 2026-10-18) ships it as the authors' derived per-pair result
# for the pairs the manuscript tests.  Anything not declared here fails.
R16_DECLARED = {("external_outcomes/S1R1_outcome_checks.tsv", "p"): "P01 (S2-D34(2)): the single-instrument Wald p of the seven external-check rows equals the outcome p at the instrument by construction"}


def _fcol(rows, j):
    out_ = []
    for r in rows:
        try:
            out_.append(float(r[j]) if j < len(r) and r[j] not in ("", "NA") else None)
        except ValueError:
            out_.append(None)
    return out_


def r16_derived_columns(R, release):
    release = Path(release)
    out, wh = release / "out", release / "withheld"
    npairs = nshipped = nfcols = nfrec = 0
    declared_hit = set()
    for w in sorted(wh.rglob("*")):
        if not w.is_file() or w.suffix.lower() not in (".csv", ".tsv"):
            continue
        k = w.relative_to(wh).as_posix()
        if not (out / k).is_file():
            continue
        W, O = _grid(w), _grid(out / k)
        if len(W) != len(O) or len(W) < R16_MIN_ROWS + 1 or not W[0] or not O[0]:
            continue
        def withheld_col(j):
            if W[0][j] not in O[0]:
                return True
            jo_ = O[0].index(W[0][j])
            return any(jo_ < len(O[i_]) and j < len(W[i_]) and O[i_][jo_] != W[i_][j] for i_ in range(1, len(W)))
        # S2-D41 (W3): a withheld F-like column is z^2; no shipped column may be F, |z| = sqrt(F) or the p of z, and F may not be recoverable from shipped b / se1 / se2
        for f_ in [j for j, c in enumerate(W[0]) if R16_F.search(c) and withheld_col(j)]:
            Fv = _fcol(W[1:], f_)
            zf = [(math.sqrt(v) if (v is not None and v >= 0) else None) for v in Fv]
            if sum(v is not None for v in zf) < R16_MIN_ROWS:
                continue
            nfcols += 1
            trf = {"|z| = sqrt(F)": zf, "F": Fv, "2-sided normal p of z": [None if v is None or v > 37 else math.erfc(v / math.sqrt(2)) for v in zf]}
            for jo, cname in enumerate(O[0]):
                if cname == W[0][f_]:
                    continue
                ov = _fcol(O[1:], jo)
                nv = sum(v is not None for v in ov)
                if nv < R16_MIN_ROWS:
                    continue
                for tname, tv in trf.items():
                    m = sum(1 for a_, b2 in zip(ov, tv) if a_ is not None and b2 not in (None, 0) and abs(a_ - b2) <= R16_REL_TOL * abs(b2))
                    if m >= R16_MIN_ROWS and m >= 0.5 * nv:
                        R.fail("R16", "%s: shipped column %r equals %s of the withheld F column %r on %d of %d numeric rows (a deterministic function of a withheld column must be withheld)" % (k, cname, tname, W[0][f_], m, nv))
            bo = [jo for jo, c in enumerate(O[0]) if R16_BETA.search(c)]
            so = [jo for jo, c in enumerate(O[0]) if R16_SE.search(c)]
            for b_ in bo:
                for s1 in so:
                    for s2 in so:
                        if s1 == s2:
                            continue
                        bv, v1, v2 = _fcol(O[1:], b_), _fcol(O[1:], s1), _fcol(O[1:], s2)
                        rec = [(x * x / (y * y - w * w)) if (x is not None and y is not None and w is not None and y * y > w * w) else None for x, w, y in zip(bv, v1, v2)]
                        nr = sum(1 for a_, b2 in zip(rec, Fv) if a_ is not None and b2 not in (None, 0))
                        m = sum(1 for a_, b2 in zip(rec, Fv) if a_ is not None and b2 not in (None, 0) and abs(a_ - b2) <= R16_REL_TOL * abs(b2))
                        nfrec += 1
                        if m >= R16_MIN_ROWS and m >= 0.5 * nr:
                            R.fail("R16", "%s: the withheld F column %r is recoverable from shipped columns: %r^2 / (%r^2 - %r^2) equals it on %d of %d rows (withhold the delta-method SE and p with F)" % (k, W[0][f_], O[0][b_], O[0][s2], O[0][s1], m, nr))
        bj = [j for j, c in enumerate(W[0]) if R16_BETA.search(c) and withheld_col(j)]
        sj = [j for j, c in enumerate(W[0]) if R16_SE.search(c) and withheld_col(j)]
        if not bj or not sj:
            continue
        skip = {W[0][j] for j in bj + sj}
        for b_ in bj:
            for s_ in sj:
                bv, sv = _fcol(W[1:], b_), _fcol(W[1:], s_)
                z = [(abs(x / y) if (x is not None and y not in (None, 0.0)) else None) for x, y in zip(bv, sv)]
                if sum(v is not None for v in z) < R16_MIN_ROWS:
                    continue
                npairs += 1
                tr = {"|z|": z, "z^2 (F)": [None if v is None else v * v for v in z],
                      "2-sided normal p of z": [None if v is None or v > 37 else math.erfc(v / math.sqrt(2)) for v in z]}
                for jo, cname in enumerate(O[0]):
                    if cname in skip:
                        continue
                    ov = _fcol(O[1:], jo)
                    nv = sum(v is not None for v in ov)
                    if nv < R16_MIN_ROWS:
                        continue
                    nshipped += 1
                    for tname, tv in tr.items():
                        m = 0
                        for a_, b2 in zip(ov, tv):
                            if a_ is None or b2 is None or b2 == 0:
                                continue
                            if abs(a_ - b2) <= R16_REL_TOL * abs(b2):
                                m += 1
                        if m >= R16_MIN_ROWS and m >= 0.5 * nv:
                            if (k, cname) in R16_DECLARED:
                                declared_hit.add((k, cname))
                                continue
                            R.fail("R16", "%s: shipped column %r equals %s of withheld %r / %r on %d of %d numeric rows (a deterministic function of withheld columns must be withheld)" % (
                                k, cname, tname, W[0][b_], W[0][s_], m, nv))
    for dk in sorted(set(R16_DECLARED) - declared_hit):
        R.warns.append("R16 declared derived column no longer matches a withheld function: %s %s" % dk)
    R.note("R16", "%d withheld beta x se pair(s) in tables that ship a counterpart; %d shipped numeric column(s) tested against |z|, z^2 and the normal p of z; %d declared (P01) column(s) used; %d withheld F column(s) tested the same way and for recovery from %d shipped b / se1 / se2 triple(s)" % (npairs, nshipped, len(declared_hit), nfcols, nfrec))


# ---------------------------------------------------------------- R17: the notes of a masked row are masked (S2-D40, W2)
R17_EST = ("value", "se", "p", "ci_lo", "ci_hi")
R17_NOTE = ("notes", "note", "comment", "comments")
R17_MARKS = ("withheld_pending_terms", "withheld_hold_H9")


def r17_masked_row_notes(R, rel, text_of):
    nfiles = nrows = nnotes = 0
    for k in sorted(rel):
        if not k.endswith((".tsv", ".csv")) or text_of.get(k) is None:
            continue
        grid = list(csv.reader(io.StringIO(text_of[k], newline=""), delimiter="\t" if k.endswith(".tsv") else ","))
        if not grid:
            continue
        hdr = grid[0]
        ei = [j for j, c in enumerate(hdr) if c in R17_EST]
        ni = [j for j, c in enumerate(hdr) if c in R17_NOTE]
        if not ei or not ni:
            continue
        nfiles += 1
        for ln, row in enumerate(grid[1:], 2):
            if any(j < len(row) and row[j] in R17_MARKS for j in ei):
                nrows += 1
                for j in ni:
                    c = row[j] if j < len(row) else ""
                    if c in ("", "NA") or c in R17_MARKS:
                        nnotes += 1
                        continue
                    R.fail("R17", "%s line %d: a value / se / p cell is withheld and %s carries %r (the notes of a masked row must be masked too)" % (k, ln, hdr[j], c[:70]))
    R.note("R17", "%d shipped table(s) with value cells and a notes column; %d masked row(s), %d of their notes cell(s) empty or masked" % (nfiles, nrows, nnotes))


# ---------------------------------------------------------------- R18: provider values in free text (S2-D41, W4)
R18_MIN_SIG = 5


def r18_notes_provider_values(R, release, rel, text_of, ledger):
    release = Path(release)
    out, wh = release / "out", release / "withheld"
    vals = leak_values(release, ledger)
    allow_p = release / LEAK_ALLOW
    s2d35, quoted = set(), set()
    if allow_p.is_file():
        for r in read_tsv(allow_p):
            kk = (r["out_path"], leak_key(r["value"]))
            if r["reason"].strip().startswith(("S2-D35", "published_quote")):       # printed by the manuscript, or public in a cited paper's own supplement
                s2d35.add(kk)
            elif r["reason"].strip().startswith("quoted_value"):
                quoted.add(kk)
    zr, nz = {}, 0
    for w in sorted(wh.rglob("*")):
        if not w.is_file() or w.suffix.lower() not in (".csv", ".tsv"):
            continue
        k = w.relative_to(wh).as_posix()
        if not (out / k).is_file():
            continue
        W, O = _grid(w), _grid(out / k)
        if len(W) != len(O) or len(W) < 2 or not W[0] or not O[0]:
            continue
        bj = [j for j, c in enumerate(W[0]) if R16_BETA.search(c) or c == "value"]
        sj = [j for j, c in enumerate(W[0]) if R16_SE.search(c)]
        for b_ in bj:
            for s_ in sj:
                for i in range(1, len(W)):
                    try:
                        b2, s2 = float(W[i][b_]), float(W[i][s_])
                    except (ValueError, IndexError):
                        continue
                    if s2 == 0:
                        continue
                    jb, js = (O[0].index(W[0][b_]) if W[0][b_] in O[0] else None), (O[0].index(W[0][s_]) if W[0][s_] in O[0] else None)
                    masked = (jb is None or i >= len(O) or jb >= len(O[i]) or O[i][jb] != W[i][b_]) or (js is None or js >= len(O[i]) or O[i][js] != W[i][s_])
                    if not masked:
                        continue
                    z = abs(b2 / s2)
                    nz += 1
                    for n in range(R18_MIN_SIG, 16):
                        zr.setdefault(n, set()).add(float("%.*g" % (n, z)))
    nfiles = ntok = nbad = 0
    for k in sorted(rel):
        if not k.endswith((".tsv", ".csv")) or text_of.get(k) is None:
            continue
        grid = list(csv.reader(io.StringIO(text_of[k], newline=""), delimiter="\t" if k.endswith(".tsv") else ","))
        if not grid:
            continue
        ni = [j for j, c in enumerate(grid[0]) if c in R17_NOTE]
        if not ni:
            continue
        nfiles += 1
        for ln, row in enumerate(grid[1:], 2):
            for j in ni:
                cell = row[j] if j < len(row) else ""
                for tok in LEAK_FL.findall(cell):
                    n = leak_sig(tok)
                    if n < LEAK_MIN_SIG:
                        continue
                    ntok += 1
                    key = leak_key(tok)
                    if key in vals and (k, key) in quoted:
                        nbad += 1
                        R.fail("R18", "%s line %d: the notes of an unmasked row quote %s, a withheld value of %s (a quoted_value declaration does not excuse a notes cell)" % (k, ln, tok, "; ".join(sorted(vals[key][0]))))
                    elif n >= R18_MIN_SIG and n <= 15 and float("%.*g" % (n, float(tok))) in zr.get(n, ()) and (k, key) not in s2d35:
                        nbad += 1
                        R.fail("R18", "%s line %d: the notes print %s, the |z| = |value / se| of a row whose value / se cells are withheld (a provider marginal z at a variant)" % (k, ln, tok))
    R.note("R18", "%d shipped table(s) with a notes column, %d numeric token(s) of >= %d significant digits read against %d withheld values and %d withheld-pair |z| values; %d offending" % (nfiles, ntok, LEAK_MIN_SIG, len(vals), nz, nbad))


# ---------------------------------------------------------------- R19: the journal files print no value the release withholds (S2-D43)
R19_DECL = "_gate_fixtures/R19_JOURNAL_DECLARED.tsv"
R19_COLS = ["file", "sheet", "column", "value", "basis", "reason"]
JOURNAL = {"AF1": "R2_rewrite/build/P31R_Additional_file_1.xlsx", "AF2": "R2_rewrite/build/P31R_Additional_file_2.docx", "MS": MANUSCRIPT_DOCX}
R19_H9 = ("S7a", "colocalization/table_s_pair_evidence_tiered.csv", "eas_beta_mr", "withheld_hold_H9")


def r19_af1_cells(p):
    """Every number in every cell of the workbook: (sheet, column header, first two cells of the row, token, cell is numeric).
    A numeric cell is read as Python's shortest repr; a text cell by LEAK_FL. The header is the first row with at least three filled cells."""
    import openpyxl
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    out_, n, titles = [], 0, {}
    for ws in wb.worksheets:
        rows = [list(r) for r in ws.iter_rows(values_only=True)]
        titles[ws.title] = str(rows[0][0]) if rows and rows[0] and rows[0][0] is not None else ""
        hi = next((i for i, r in enumerate(rows) if sum(c is not None for c in r) >= 3), None)
        hdr = rows[hi] if hi is not None else []
        for i, r in enumerate(rows):
            for j, c in enumerate(r):
                if c is None:
                    continue
                n += 1
                num = isinstance(c, (int, float)) and not isinstance(c, bool)
                col = str(hdr[j]) if hi is not None and i > hi and j < len(hdr) else "(title, legend or header)"
                for t in ([repr(float(c))] if num else LEAK_FL.findall(str(c))):
                    out_.append((ws.title, col, tuple(str(x) for x in r[:2]), t, num))
    return out_, n, titles


def r19_journal(R, release, v3, text_of, ledger, journal=None):
    """Stage 6.3b part 2b fix round (S2-D43). R14 scans the release only; Additional file 1 had printed the hold-H9 East Asian eGFR magnitudes and the UKB-PPP instrument
    values that the release withholds, and no check compared the two. Here every number in the manuscript docx and Additional files 1 and 2 is matched against the release's
    withheld values (R14's list and rules). A value the release also ships elsewhere is public and passes; any other hit needs a declaration in R19_DECL (file, sheet,
    column, value or '*'). A declaration whose reason starts with PENDING-USER passes with a warning. Separately, every hold-H9 pair of the release's tiered table must sit
    in S7a and print no number there except its p."""
    release, v3 = Path(release), Path(v3)
    paths = {k: (Path(journal[k]) if journal and k in journal else v3 / p_) for k, p_ in JOURNAL.items()}
    for k, p_ in paths.items():
        if not p_.is_file():
            R.fail("R19", "journal file %s is absent: %s" % (k, p_))
            return
    dp = release / R19_DECL
    if not dp.is_file():
        R.fail("R19", "%s is absent (the declared journal exceptions)" % R19_DECL)
        return
    if dp.read_text(encoding="utf-8").splitlines()[0].split("\t") != R19_COLS:
        R.fail("R19", "%s columns differ from %s" % (R19_DECL, R19_COLS))
        return
    decl = read_tsv(dp)
    for d_ in decl:
        if not d_["basis"].strip() or not d_["reason"].strip():
            R.fail("R19", "declaration without a basis or a reason: %s %s %s" % (d_["file"], d_["sheet"], d_["column"]))
    vals = leak_values(release, ledger)
    idx = {}
    for (srcs, tok) in vals.values():
        for n in range(LEAK_SHIP_MIN_SIG, min(leak_sig(tok), 15) + 1):
            idx.setdefault(n, {}).setdefault(leak_round(tok, n), (srcs, tok))
    ship_k, ship_r = set(), {}
    for txt in text_of.values():
        for t in LEAK_FL.findall(txt or ""):
            ship_k.add(leak_key(t))
            for n in range(LEAK_SHIP_MIN_SIG, min(leak_sig(t), 15) + 1):
                ship_r.setdefault(n, set()).add(leak_round(t, n))
    if len(ship_k) < 1000:
        R.fail("R19", "the shipped text yielded %d numeric tokens: the public-value set is implausibly small" % len(ship_k))
        return

    def withheld(tok):
        n = leak_sig(tok)
        if n < LEAK_MIN_SIG:
            return None
        if leak_key(tok) in ship_k or (LEAK_SHIP_MIN_SIG <= n <= 15 and leak_round(tok, n) in ship_r.get(n, ())):
            return None
        return vals.get(leak_key(tok)) or (idx.get(n, {}).get(leak_round(tok, n)) if LEAK_SHIP_MIN_SIG <= n <= 15 else None)

    used, pend, nhit, ncov = set(), {}, 0, 0

    def judge(file_, sheet, col, where, tok):
        nonlocal nhit, ncov
        h = withheld(tok)
        if h is None:
            return
        nhit += 1
        d_ = next((i for i, d_ in enumerate(decl) if d_["file"] == file_ and d_["sheet"] == sheet and d_["column"] in (col, "*") and d_["value"] in (tok, "*")), None)
        if d_ is None:
            R.fail("R19", "%s %s, column %r (%s) prints %s, a value the release withholds (%s), and no declaration covers it" % (
                file_, sheet, col, where, tok, "; ".join(sorted(h[0]))[:160]))
            return
        used.add(d_)
        ncov += 1
        if decl[d_]["reason"].startswith("PENDING-USER"):
            pend[d_] = pend.get(d_, 0) + 1
    cells, ncell, titles = r19_af1_cells(paths["AF1"])
    if ncell < 5000 or not titles.get("S1", "").startswith("Table S1."):
        R.fail("R19", "the Additional file 1 reader recovered %d cells and S1 title %r: blind" % (ncell, titles.get("S1")))
        return
    for (sh, col, key2, tok, num) in cells:
        judge("AF1", sh, col, "row %s" % "/".join(key2), tok)
    nchar = {}
    for k in ("AF2", "MS"):
        txt, _ = docx_text(paths[k])
        nchar[k] = len(txt)
        if len(txt) < 50000 or "Additional file" not in txt:
            R.fail("R19", "the %s reader recovered %d characters: blind" % (k, len(txt)))
            continue
        for tok in LEAK_FL.findall(txt):
            judge(k, "text", "*", "running text", tok)
    # hold H9: the pairs whose East Asian magnitudes the release withholds print no number in S7a except their p
    sh9, tab9, col9, mark9 = R19_H9
    T = _grid(release / "out" / tab9)
    gi, oi, bi = T[0].index("gene"), T[0].index("outcome"), T[0].index(col9)
    h9 = {(r[gi], r[oi]) for r in T[1:] if len(r) > bi and r[bi] == mark9}
    if not h9:
        R.fail("R19", "no %s row in %s: the hold-H9 rule cannot run" % (mark9, tab9))
    seen9 = {key2 for (sh, col, key2, tok, num) in cells if sh == sh9 and key2 in h9}
    if seen9 != h9:
        R.fail("R19", "%s does not hold every hold-H9 pair; missing %s" % (sh9, sorted(h9 - seen9)[:5]))
    n9 = 0
    for (sh, col, key2, tok, num) in cells:
        if sh == sh9 and key2 in h9 and num and col != "p":
            n9 += 1
            R.fail("R19", "%s row %s prints %s in column %r: a hold-H9 pair may print only its p" % (sh9, "/".join(key2), tok, col))
    for i in sorted(pend):
        d_ = decl[i]
        R.warns.append("R19 %d value(s) pass on a PENDING-USER declaration: %s %s %r (%s)" % (pend[i], d_["file"], d_["sheet"], d_["column"], d_["basis"]))
    for i, d_ in enumerate(decl):
        if i not in used:
            R.warns.append("R19 declaration no longer used: %s %s %r value %s" % (d_["file"], d_["sheet"], d_["column"], d_["value"]))
    R.note("R19", "%d Additional file 1 cells, %s characters of AF2 / manuscript read; %d withheld values; %d journal tokens equal a withheld value the release does not ship, "
                  "%d of them covered by %d declaration(s) (%d PENDING-USER); hold-H9 pairs %d, all in %s, %d printing a number other than p" % (
                      ncell, "/".join(str(nchar.get(k, 0)) for k in ("AF2", "MS")), len(vals), nhit, ncov, len(used), len(pend), len(h9), sh9, n9))


def run_gate(release, v3, manuscript=None, ledger_path=None, journal=None):
    R = Res()
    release, v3 = Path(release), Path(v3)
    out, wh = release / "out", release / "withheld"
    files = sorted(p for p in out.rglob("*") if p.is_file())
    rel = {p.relative_to(out).as_posix(): p for p in files}
    ledger = load_ledger(v3, ledger_path)

    # ---- R0 scope + reader
    opened, text_of, ratios = set(), {}, []
    for k, p in rel.items():
        t = read_text(p)
        opened.add(k)
        text_of[k] = t
        if t is not None and p.stat().st_size > 1024 and not k.endswith(".gz"):
            ratios.append((k, len(t) / p.stat().st_size))
    if opened != set(rel):
        R.fail("R0", "gate did not open %d file(s)" % len(set(rel) - opened))
    for k, r in ratios:
        if not 0.5 <= r <= 1.0:
            R.fail("R0", "reader chars/byte %.2f outside [0.5, 1.0] for %s (blind or flooded)" % (r, k))
    probe = "screen/screen_ledger.csv"
    if probe in text_of and text_of[probe] is not None:
        if "fdr_sig_bh_family" not in text_of[probe].splitlines()[0]:
            R.fail("R0", "reader did not recover 'fdr_sig_bh_family' from the header of %s" % probe)
    else:
        R.fail("R0", "probe file %s absent or binary" % probe)
    R.note("R0", "opened %d/%d files in out/ (%d text, %d binary); chars/byte band checked on %d files" % (
        len(opened), len(rel), sum(t is not None for t in text_of.values()), sum(t is None for t in text_of.values()), len(ratios)))

    # ---- R1 / R2 manifest
    mpath = out / "RELEASE_MANIFEST.tsv"
    man = []
    if not mpath.exists():
        R.fail("R1", "RELEASE_MANIFEST.tsv is absent")
    else:
        hdr = mpath.read_text(encoding="utf-8").splitlines()[0].split("\t")
        if hdr != MANIFEST_COLS:
            R.fail("R1", "manifest columns differ: %s" % hdr)
        man = read_tsv(mpath)
    mp = [r.get("path", "") for r in man]
    if len(mp) != len(set(mp)):
        R.fail("R1", "duplicate manifest paths")
    for q in mp:
        if q in SELF_EXEMPT:
            R.fail("R1", "manifest has a row for %s (exempt, cannot carry its own md5)" % q)
    shipped = {k for k in rel if k not in SELF_EXEMPT}
    for k in sorted(shipped - set(mp)):
        R.fail("R1", "file in out/ not in manifest: %s" % k)
    for k in sorted(set(mp) - set(rel)):
        R.fail("R1", "manifest row without a file: %s" % k)
    for k in SELF_EXEMPT:
        if k not in rel:
            R.fail("R1", "%s is absent from out/" % k)
    R.note("R1", "%d files in out/ (%d exempt), %d manifest rows" % (len(rel), len(SELF_EXEMPT & set(rel)), len(man)))

    nsrc = nprod = nrc = 0
    for r in man:
        k = r.get("path", "")
        if k not in rel:
            continue
        p = rel[k]
        if str(p.stat().st_size) != r["bytes"]:
            R.fail("R2", "%s bytes %s != %d" % (k, r["bytes"], p.stat().st_size))
        if md5_file(p) != r["md5"]:
            R.fail("R2", "%s md5 differs from the manifest" % k)
        sp, sm = r["source_path"].split(";"), r["source_md5"].split(";")
        if len(sp) != len(sm):
            R.fail("R2", "%s: %d source paths vs %d source md5s" % (k, len(sp), len(sm)))
        for a, b in zip(sp, sm):
            a = a.strip()
            if not a:
                continue
            f = Path(a) if re.match(r"^[A-Za-z]:|^/", a) else (v3.parent / a[len("CKM papers/"):] if a.startswith("CKM papers/") else v3 / a)
            nsrc += 1
            if not f.is_file():
                R.fail("R2", "%s: source missing now: %s" % (k, a))
            elif md5_file(f) != b.strip():
                R.fail("R2", "%s: source md5 changed: %s" % (k, a))
        ps = r["producer_script"].split(":")[0].strip()
        if ps and r["producer_md5"].strip():
            f = next((c for c in (v3 / ps, release / ps, release / "out" / ps) if c.is_file()), None)
            if f is not None:
                nprod += 1
                if md5_file(f) != r["producer_md5"].strip():
                    R.fail("R2", "%s: producer %s md5 differs from producer_md5" % (k, ps))
        if k.endswith((".csv", ".tsv")) and re.fullmatch(r"\d+", r["rows"] or "") and re.fullmatch(r"\d+", r["columns"] or ""):
            rows = list(csv.reader(io.StringIO(text_of[k], newline=""), delimiter="\t" if k.endswith(".tsv") else ","))
            nrc += 1
            if (len(rows) - 1, len(rows[0])) != (int(r["rows"]), int(r["columns"])):
                R.fail("R2", "%s: rows/columns %s/%s in manifest, %d/%d in file" % (k, r["rows"], r["columns"], len(rows) - 1, len(rows[0])))
    R.note("R2", "%d files md5/bytes, %d source md5s, %d producer md5s, %d row/column counts re-derived" % (len(man), nsrc, nprod, nrc))

    # ---- R3 not shipped
    ns = read_tsv(release / "inputs" / "NOT_SHIPPED.tsv")
    nsmd5 = {}
    for r in ns:
        nsmd5.setdefault(r["md5"], []).append(r["path"])
    seen = 0
    for base in (out, wh):
        for p in sorted(x for x in base.rglob("*") if x.is_file()):
            seen += 1
            m = md5_file(p)
            if m in nsmd5:
                R.fail("R3", "%s/%s has the md5 of a not-shipped file (%s)" % (base.name, p.relative_to(base).as_posix(), nsmd5[m][0]))
    R.note("R3", "%d files in out/ + withheld/ against %d not-shipped md5s" % (seen, len(nsmd5)))

    # ---- R3b REPRODUCE
    rep = out / "REPRODUCE.md"
    units = {r["unit"]: r["class"] for r in read_tsv(release / "build" / "regenerability_by_unit.tsv")}
    srcs = read_tsv(release / "build" / "regenerability_by_source.tsv")
    external = {r["source_path"] for r in srcs if r["class"] in ("needs_public_download", "needs_nonredistributable_input")}
    external |= set(nsmd5_paths(ns))
    ncmd = 0
    if not rep.exists():
        R.fail("R3b", "REPRODUCE.md is absent")
    else:
        lines = text_of["REPRODUCE.md"].splitlines()
        infence, stated = False, {}
        for ln in lines:
            if ln.strip().startswith("```"):
                infence = not infence
                continue
            if infence:
                continue
            norm = re.sub(r"[^a-z0-9]+", "_", ln.lower())
            cls = {c for c in CLASSES if re.search(r"(?<![a-z0-9])%s(?![a-z0-9])" % c, norm)}
            cls |= {CLASS_ALIASES[a] for a in CLASS_ALIASES if re.search(r"(?<![a-z0-9])%s(?![a-z0-9])" % a, norm)}
            if not cls:
                continue
            lnu = re.sub(r"`[^`]*`", " ", ln)       # a path in a code span (".../Figure_3/...") is an input, not a unit name
            us = [u for u in units if re.search(r"(?<![\w])%s(?![\w])" % re.escape(u), lnu)]
            if len(us) >= 1 and len(cls) == 1:      # one class on the line: every unit named on it is stated with that class (each is checked below)
                for u_ in us:
                    stated.setdefault(u_, set()).update(cls)
            elif len(us) > 1:
                R.fail("R3b", "line names several units and several classes (ambiguous): %s" % ln[:90])
            elif len(us) == 1:
                stated.setdefault(us[0], set()).update(cls)
        # V04: the unit table is parsed AS a table: every row's class cell is exactly one known class; every derived unit has exactly one row
        tab_rows, in_tab, ccol = collections.defaultdict(list), False, None
        known = set(CLASSES) | set(CLASS_ALIASES)
        for ln in lines:
            cells = [c.strip() for c in ln.strip().strip("|").split("|")] if ln.strip().startswith("|") else None
            if cells is None:
                in_tab, ccol = False, None
                continue
            low = [c.lower() for c in cells]
            if "unit" in low and "class" in low:
                in_tab, ccol = True, low.index("class")
                continue
            if not in_tab or all(re.fullmatch(r":?-{2,}:?", c) for c in cells):
                continue
            cell = re.sub(r"[^a-z0-9]+", "_", cells[ccol].replace("`", "").lower()).strip("_") if ccol < len(cells) else ""
            if cell not in known:
                R.fail("R3b", "unit table row %r: class cell %r is not exactly one known class" % (cells[0], cells[ccol] if ccol < len(cells) else ""))
            else:
                tab_rows[cells[0]].append(CLASS_ALIASES.get(cell, cell))
        for u, c in units.items():
            if len(tab_rows.get(u, [])) != 1:
                R.fail("R3b", "unit %s has %d row(s) in the REPRODUCE unit table (exactly one required)" % (u, len(tab_rows.get(u, []))))
            elif tab_rows[u][0] != c:
                R.fail("R3b", "unit table row %s: class %s, derived class is %s" % (u, tab_rows[u][0], c))
        for u, c in units.items():
            if u not in stated:
                R.fail("R3b", "unit %s has no stated class in REPRODUCE.md" % u)
            elif stated[u] != {c}:
                R.fail("R3b", "unit %s: REPRODUCE says %s, derived class is %s" % (u, sorted(stated[u]), c))
        infence = False
        for ln in lines:
            s = ln.strip()
            if s.startswith("```"):
                infence = not infence
                continue
            s = s.lstrip("$ ").strip()
            if not (infence and (s.startswith("python -I") or s.startswith("Rscript --vanilla"))):
                continue
            ncmd += 1
            try:
                toks = shlex.split(s)
            except ValueError:
                R.fail("R3b", "unparseable command: %s" % s[:90])
                continue
            script_seen, skip = False, False
            for i, t in enumerate(toks[1:], 1):
                if skip:
                    skip = False
                    continue
                if t in ("--out", "--output", "-o", "--outdir", ">"):
                    skip = True
                    continue
                if t == ">>" or t.startswith("-") or t in ("python", "Rscript"):
                    continue
                if not (("/" in t) or FILE_TOKEN.fullmatch(t) or t.endswith((".gz", ".tbi"))):
                    continue
                t = t.replace("\\", "/")
                if t in rel or t in external:
                    script_seen = True
                    continue
                R.fail("R3b", "command input not in out/ and not a declared external input: %s  (%s)" % (t, s[:70]))
        if ncmd == 0 and any(c == "regenerates_from_release" for c in units.values()):
            R.fail("R3b", "a unit regenerates from the release but REPRODUCE.md has no command")
    R.note("R3b", "%d units against regenerability_by_unit.tsv, %d commands" % (len(units), ncmd))

    # ---- R4 individual-level
    nwide = 0
    for k, p in rel.items():
        nm = k.lower()[:-3] if k.lower().endswith(".gz") else k.lower()
        if nm.endswith(INDIV_EXT):
            R.fail("R4", "individual-level genotype file extension: %s" % k)
            continue
        if nm.endswith((".csv", ".tsv")) and text_of[k]:
            hdr = text_of[k].splitlines()[0].split("\t" if nm.endswith(".tsv") else ",")
            nwide += 1
            vh = sum(bool(re.fullmatch(r"(?:rs\d+|(?:chr)?\d+[:_]\d+(?:[:_][ACGT]+)*)", h.strip('"'))) for h in hdr)
            if len(hdr) > MAX_COLS or vh > MAX_VARIANT_HEADERS:
                R.fail("R4", "%s looks like samples x variants (%d columns, %d variant-like headers)" % (k, len(hdr), vh))
    R.note("R4", "extension rule on %d files; column band <= %d, variant headers <= %d on %d tables" % (len(rel), MAX_COLS, MAX_VARIANT_HEADERS, nwide))

    # ---- R5 secrets
    scan = dict(text_of)
    ci = release / "inputs" / "CODE_INCLUDE.tsv"
    nextra = 0
    if ci.exists():
        for r in read_tsv(ci):
            sh = str(r.get("ship", "")).strip().lower() in ("1", "true", "yes")
            pth = r.get("source_path") or r.get("path") or r.get("file") or ""
            f = (v3 / pth) if pth else None
            if sh and f is not None and f.is_file():
                t = read_text(f)
                if t is not None:
                    scan["[code] " + pth] = t
                    nextra += 1
    nblock = nflag = 0
    agg = {}
    for k, t in scan.items():
        if t is None:
            continue
        for i, ln in enumerate(t.splitlines(), 1):
            for nm, rx in SECRET_BLOCK:
                if rx.search(ln):
                    nblock += 1
                    R.fail("R5", "%s:%d  %s" % (k, i, nm))
            for nm, rx in SECRET_FLAG:
                if rx.search(ln):
                    nflag += 1
                    agg.setdefault((k, nm), [i, 0])[1] += 1
    for (k, nm), (first, n) in sorted(agg.items()):
        R.warns.append("R5 %s:%d  %s (%d line(s))" % (k, first, nm, n))
    R.note("R5", "%d text files scanned (%d from the code list outside out/), %d blocking, %d flagged (non-blocking)" % (
        sum(t is not None for t in scan.values()), nextra, nblock, nflag))

    # ---- R6 manuscript names
    ext = {}
    fx = release / "_gate_fixtures" / "R6_EXTERNAL_NAMES.tsv"      # relative to --release (V05)
    if fx.exists():
        ext = {r["name"]: r["reason"] for r in read_tsv(fx)}
    if manuscript is None:
        dd = v3 / MANUSCRIPT_DOCX
        if not dd.is_file():
            R.fail("R6", "default manuscript %s is absent" % MANUSCRIPT_DOCX)
        else:   # Stage 6.2 (S2-D32 repoint b): refuse a built docx older than any of its sources in the component source map
            mp_files = r6_docx_sources(v3)
            if not mp_files:
                R.fail("R6", "the component source map resolved to no files (staleness check cannot run)")
            else:
                newest = r6_stale(dd, mp_files)
                if newest is not None:
                    R.fail("R6", "built docx %s is older than its source %s: rebuild it (python -I scripts/P31R_build.py)" % (MANUSCRIPT_DOCX, newest.name))
        mfiles, mode = [dd] if dd.is_file() else [], "built docx %s" % MANUSCRIPT_DOCX
    else:
        mp_ = Path(manuscript)
        mfiles = [mp_] if mp_.is_file() else sorted(mp_.glob("*.md"))
        mode = "docx %s" % mp_ if mp_.suffix == ".docx" else "markdown %s" % mp_
    mtext, mbytes = [], 0
    for f in mfiles:
        if f.suffix == ".docx":
            t, nb = docx_text(f)
            mbytes += nb
        else:
            t = f.read_bytes().decode("utf-8")
            mbytes += f.stat().st_size
        mtext.append(t)
    alltext = "\n".join(mtext)
    if not mfiles or not alltext.strip():
        R.fail("R6", "manuscript reader returned nothing")
    elif manuscript is None or Path(manuscript).suffix == ".docx":      # reader liveness: a known string comes back, and the chars/byte ratio is two-sided
        ps = read_tsv(release / "inputs" / "PROVIDERS.tsv")[0]["provider_id"]
        if ps not in alltext:
            R.fail("R6", "reader did not recover the known string %r" % ps)
        lo, hi = (DOCX_BAND if mfiles[0].suffix == ".docx" else (0.5, 1.0))
        if not lo <= len(alltext) / max(mbytes, 1) <= hi:
            R.fail("R6", "manuscript reader chars/byte %.2f outside [%.2f, %.2f]" % (len(alltext) / max(mbytes, 1), lo, hi))
    names = sorted(set(FILE_TOKEN.findall(alltext)))
    base_ix = {Path(k).name for k in rel}
    for n in names:
        n2 = n.replace("\\", "/")
        if n2 in rel or Path(n2).name in base_ix or n2 in ext or Path(n2).name in ext:
            continue
        R.fail("R6", "manuscript names %s, which is not in out/" % n)
    R.note("R6", "mode: %s; %d file(s) read, %d names" % (mode, len(mfiles), len(names)))

    # ---- R7 counts from the ledger
    nnum = 0
    for fn in ("README.md", "REPRODUCE.md"):
        if fn not in rel:
            R.fail("R7", "%s is absent" % fn)
            continue
        raw = text_of[fn]
        raw = re.sub(r"```.*?```", lambda m: "\n" * m.group(0).count("\n"), raw, flags=re.S)
        blocks, cur = [], []
        for ln in raw.splitlines():
            if ln.strip().startswith("|"):
                blocks.append([ln])
            elif ln.strip() == "":
                if cur:
                    blocks.append(cur)
                cur = []
            else:
                cur.append(ln)
        if cur:
            blocks.append(cur)
        for b in blocks:
            para = "\n".join(b)
            nnum += r7_para(R, fn, para, ledger)
    R.note("R7", "%d numerals in README/REPRODUCE prose, each matched to a ledger row beside it" % nnum)

    # ---- R8 screen
    sp = "screen/screen_ledger.csv"
    if sp not in rel:
        R.fail("R8", "screen_ledger.csv absent")
    else:
        rows = list(csv.DictReader(io.StringIO(text_of[sp], newline="")))
        lv = lambda i: json.loads(ledger[i]["value"])
        if len(rows) != int(lv("L003")):
            R.fail("R8", "rows %d != L003 %s" % (len(rows), lv("L003")))
        if sum(r["testable"] == "TRUE" for r in rows) != int(lv("L004")):
            R.fail("R8", "testable rows != L004 %s" % lv("L004"))
        p = [float(r["p_family"]) for r in rows]
        m, alpha = len(p), float(lv("L239"))
        order = sorted(range(m), key=lambda i: p[i])
        q, mn = [0.0] * m, 1.0
        for rank in range(m, 0, -1):
            i = order[rank - 1]
            mn = min(mn, p[i] * m / rank)
            q[i] = mn
        flag = [qi < alpha for qi in q]
        if sum(flag) != int(lv("L008")):
            R.fail("R8", "recomputed BH admits %d, L008 is %s" % (sum(flag), lv("L008")))
        if [r["fdr_sig_bh_family"] == "TRUE" for r in rows] != flag:
            R.fail("R8", "fdr_sig_bh_family differs from BH recomputed over p_family at L239")
        if not all(abs(float(r["fdr_q_bh_family"]) - qi) < 1e-9 for r, qi in zip(rows, q)):
            R.fail("R8", "fdr_q_bh_family differs from the recomputed BH q")
        if "fdr_sig_legacy_per_method" not in rows[0]:
            R.fail("R8", "column fdr_sig_legacy_per_method absent (legacy fdr_sig must be renamed)")
        elif sum(r["fdr_sig_legacy_per_method"] == "TRUE" for r in rows) != int(lv("L009")):
            R.fail("R8", "legacy flag count != L009 %s" % lv("L009"))
    ncol = r8_conservative(R, v3, rel, text_of, man, ledger)
    R.note("R8", "screen rows/testable/BH/legacy against L003/L004/L008/L009; %d table headers scanned; a conservative-prior column ships only in a declared file and column, on rows of the canonical pairs (+ Table 2's own)" % ncol)

    # ---- R12 / R13 portability: the scope is every shipped file the reader opened (R0), not code/ alone
    nexec = nlines = 0
    for k, t_ in sorted(text_of.items()):
        if t_ is None:
            continue
        is_exec = k.rsplit(".", 1)[-1].lower() in ("r", "py", "sh")
        nexec += is_exec
        for i, ln in enumerate(t_.splitlines(), 1):
            if not LOCAL_PATH_RX.search(ln):
                continue
            nlines += 1
            if is_exec and not ln.lstrip().startswith("#"):
                R.fail("R12", "%s:%d  absolute local path in a non-comment line of a shipped executable: %s" % (k, i, ln.strip()[:100]))
            R.fail("R13", "%s:%d  absolute local path: %s" % (k, i, ln.strip()[:100]))
    R.note("R12", "%d shipped executables (*.R, *.py, *.sh) anywhere in out/ read; no absolute path in a non-comment line" % nexec)
    R.note("R13", "%d shipped text files read; %d line(s) name a drive-letter or home path" % (sum(v is not None for v in text_of.values()), nlines))

    # ---- R9 allowlist, R10 ledger recount, R11 tiers
    r9_allowlist(R, release, rel)
    r10_ledger(R, ledger, v3, release)
    r11_tiers(R, release, rel, text_of, man)
    r14_leakscan(R, release, rel, text_of, ledger)
    r15_provenance(R, release, v3, rel, text_of, man)
    r16_derived_columns(R, release)
    r17_masked_row_notes(R, rel, text_of)
    r18_notes_provider_values(R, release, rel, text_of, ledger)
    r19_journal(R, release, v3, text_of, ledger, journal)
    return R


def nsmd5_paths(ns):
    return [r["path"] for r in ns]


# ---------------------------------------------------------------- self-test
def standins(rel_copy, v3, ledger):
    """Always write stand-in manifest / REPRODUCE / README into the scratch copy, so the mutation suite tests the gate and does not depend on
    what Builder D has delivered at that moment.  The manifest is derived from the files actually in the copy; source columns come from
    build/BUILT_FILES.tsv where the path matches."""
    out = rel_copy / "out"
    for fn in ("RELEASE_MANIFEST.tsv", "README.md", "REPRODUCE.md"):      # the three written as stand-ins below; LICENSE, CITATION.cff and .gitattributes stay as built (Stage 7 A2.6, Stage 10 step 3a)
        if (out / fn).exists():
            (out / fn).unlink()
    built = {b["release_path"]: b for b in read_tsv(rel_copy / "build" / "BUILT_FILES.tsv")}
    cimap = {r["ship_path"]: r for r in read_tsv(rel_copy / "inputs" / "CODE_INCLUDE.tsv") if r["ship"] == "1"}      # S6.3b2a: code rows carry their original, label and transform
    with open(out / "RELEASE_MANIFEST.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(MANIFEST_COLS)
        for p in sorted(x for x in out.rglob("*") if x.is_file() and x.relative_to(out).as_posix() not in SELF_EXEMPT):      # exempt documents carry no manifest row
            k = p.relative_to(out).as_posix()
            b = built.get(k, {})
            rows = cols = ""
            if k.endswith((".csv", ".tsv")):
                rr = list(csv.reader(io.StringIO(p.read_text(encoding="utf-8"), newline=""), delimiter="\t" if k.endswith(".tsv") else ","))
                rows, cols = str(len(rr) - 1), str(len(rr[0]))
            sp, sm = b.get("source_path", ""), b.get("source_md5", "")
            if sp and not all((v3 / x).is_file() and md5_file(v3 / x) == y for x, y in zip(sp.split(";"), sm.split(";"))):
                sp = sm = ""      # the real build is mid-flight: keep the stand-in self-consistent
            lab_ = tr_ = ""
            if k in cimap and (v3 / cimap[k]["path"]).is_file():
                sp, sm, lab_ = cimap[k]["path"], md5_file(v3 / cimap[k]["path"]), cimap[k].get("provenance_label", "")
                tr_ = "byte-identical to the original" if md5_file(p) == sm else "derived from the original (stand-in)"
            w.writerow([k, p.stat().st_size, md5_file(p), rows, cols, sp, sm, "", "", "post hoc", "", "own_derived", lab_, tr_])
    units = read_tsv(rel_copy / "build" / "regenerability_by_unit.tsv")
    ext = next(r["path"] for r in read_tsv(rel_copy / "inputs" / "NOT_SHIPPED.tsv") if not LOCAL_PATH_RX.search(r["path"]))      # a project-relative name: R13 scans REPRODUCE.md too
    L = ["# REPRODUCE (stand-in)", "", "| unit | class |", "|---|---|"]
    L += ["| %s | %s |" % (u["unit"], u["class"]) for u in units]
    L += ["", "Commands, run from the repository root.", "", "```", "Rscript --vanilla decoy/S1R1_decoy_rates.R decoy/S1R1_decoy_rates.tsv --out decoy/new.tsv",
          "python -I -X utf8 decoy/S1R1_decoy_rates.R %s" % ext, "```", "",
          "The family has %s tests (L003)." % ledger["L003"]["display"], ""]
    (out / "REPRODUCE.md").write_text("\n".join(L), encoding="utf-8")
    (out / "README.md").write_text("# README (stand-in)\n\nThe family has %s tests (L003); %s are testable (L004).\n\nThe provider list holds %s [L533]. Labels: `full run-time hash match`, `partial launch-record match`, `timestamp evidence only`.\n" % (
        ledger["L003"]["display"], ledger["L004"]["display"], ledger["L533"]["display"]), encoding="utf-8")


def man_rewrite(path, fn):
    rows = read_tsv(path)
    fn(rows)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def refresh_row(out, rows, k):
    p = out / k
    for r in rows:
        if r["path"] == k:
            r["bytes"], r["md5"] = str(p.stat().st_size), md5_file(p)
            if k.endswith((".csv", ".tsv")) and r["rows"]:
                rr = list(csv.reader(io.StringIO(p.read_text(encoding="utf-8"), newline=""), delimiter=chr(9) if k.endswith(".tsv") else ","))
                r["rows"], r["columns"] = str(len(rr) - 1), str(len(rr[0]))
            return
    raise AssertionError("no manifest row " + k)


def fake_jwt(header=b'{"alg":"HS256"}', sig="xxx", payload=b"{}"):
    """A JWT-shaped test token built at RUN TIME, so that no eyJ... literal sits in this file (a secret scanner on the deposit would flag one) (V09)."""
    b64 = lambda b: base64.urlsafe_b64encode(b).decode().rstrip("=")
    return b64(header) + "." + b64(payload) + "." + sig


def write_tsv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def strip_eas_comparators(d):
    """Stand-in for the tier withhold of a release not yet rebuilt. S2-D34(6): a file whose TIER_RULES file row declares its comparator columns (selector comparator: / derived:) keeps
    them, and a columns: rule drops exactly its own columns; only a file with no file row gets the S2-D33(1) strip (drop eur_* / het_p / *strict* columns)."""
    out = d / "out"
    frules = {r["target"].replace("\\", "/"): r for r in read_tsv(d / "inputs" / "TIER_RULES.tsv") if r.get("kind") == "file"} if (d / "inputs" / "TIER_RULES.tsv").exists() else {}
    for p in sorted((out / "east_asian").glob("*.tsv")):
        rows = list(csv.reader(io.StringIO(p.read_text(encoding="utf-8"), newline=""), delimiter="\t"))
        fr = frules.get("east_asian/" + p.name)
        if fr is not None and fr["handling"].strip().startswith("columns:"):
            gone = set(fr["handling"].strip()[len("columns:"):].split("|"))
            keep = [i for i, h in enumerate(rows[0]) if h not in gone]
        elif fr is not None:
            continue
        else:
            keep = [i for i, h in enumerate(rows[0]) if not EAS_EUR_COL.match(h)]
        if len(keep) == len(rows[0]):
            continue
        with open(p, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            for r in rows:
                w.writerow([r[i] for i in keep])


def restamp_manifest(d, v3):
    """Real-release copy: refresh bytes / md5 / row counts of changed files and the source md5 of rows whose source moved since the build (stand-in for the rebuild)."""
    out = d / "out"

    have = rel_of(out)

    def f(rows):
        for r in rows:
            if r["path"] in have:
                refresh_row(out, rows, r["path"])
            sp, sm = r["source_path"].split(";"), r["source_md5"].split(";")
            if len(sp) == len(sm) and all(x.strip() for x in sp):
                fs = [(v3 / x.strip()) for x in sp]
                if all(x.is_file() for x in fs):
                    r["source_md5"] = ";".join(md5_file(x) for x in fs)
    man_rewrite(out / "RELEASE_MANIFEST.tsv", f)


def rel_of(out):
    return {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()}


def standin_tiers(d):
    """Stand-in tiers_by_file.tsv / TIER_RULES.tsv (Agent T's files may not exist yet): per-variant tables are Tier 2 with a CC BY provider (CKB), the rest Tier 1."""
    out = d / "out"
    known = {r["provider_id"] for r in read_tsv(d / "inputs" / "PROVIDER_TERMS.tsv")}
    mlc = {r["path"]: r["licence_class"] for r in read_tsv(out / "RELEASE_MANIFEST.tsv")}
    rows = []
    for k in sorted(rel_of(out) - SELF_EXEMPT):
        hdr = []
        if k.endswith((".csv", ".tsv")):
            hdr = [h.strip('"') for h in (out / k).read_text(encoding="utf-8").splitlines()[0].split("\t" if k.endswith(".tsv") else ",")]
        if hdr and is_per_variant(hdr):
            rows.append({"file": k, "tier": "2", "providers": "CKB", "decision": "ship", "reason": "stand-in: per-variant values, CC BY provider"})
        else:
            pv = [x for x in split_ids(mlc.get(k, "").partition(":")[2]) if x in known]
            rows.append({"file": k, "tier": "1", "providers": "+".join(pv), "decision": "ship", "reason": "stand-in: derived results"})
    write_tsv(d / "build" / "tiers_by_file.tsv", rows, TIERS_COLS)
    write_tsv(d / "inputs" / "TIER_RULES.tsv", [{"rule_id": "SW1", "kind": "switch", "target": "TPMI_DERIVED", "tier": "3", "providers": "TPMI", "handling": "withheld",
                                                   "selector": "", "basis": "stand-in", "reason": "stand-in switch"}], ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])


def selftest(release, v3):
    tmp = Path(tempfile.mkdtemp(prefix="p31r_gate_"))
    ok_all = True
    try:
        _real_ok = []

        def use_real_tiers():
            if not _real_ok:
                ok_ = False
                if (Path(release) / "build" / "tiers_by_file.tsv").exists() and (Path(release) / "inputs" / "TIER_RULES.tsv").exists():
                    pr = tmp / "probe_real"
                    for sub in ("out", "withheld", "build", "inputs", "_gate_fixtures"):
                        shutil.copytree(Path(release) / sub, pr / sub)
                    strip_eas_comparators(pr)
                    restamp_manifest(pr, v3)
                    ok_ = "R11" not in run_gate(pr, v3, None).tripped()
                    shutil.rmtree(pr, ignore_errors=True)
                print("tiers_by_file.tsv / TIER_RULES.tsv for the tier-2 copies: %s" % (
                    "the REAL files (they match the real out/)" if ok_ else "STAND-INS (the real files are absent or do not yet match the real out/; rebuild first for a real-tier run)"))
                _real_ok.append(ok_)
            return _real_ok[0]

        def fresh(name, real=False):
            d = tmp / re.sub(r"[^A-Za-z0-9_]+", "_", name)
            for sub in ("out", "withheld", "build", "inputs", "_gate_fixtures"):
                shutil.copytree(Path(release) / sub, d / sub)
            led = load_ledger(v3)
            ci = d / "inputs" / "CODE_INCLUDE.tsv"      # stand-in: a code row whose file the (not yet rebuilt) release does not hold is not an approved row in the scratch copy
            have = rel_of(d / "out")
            cr = read_tsv(ci)
            for r in cr:
                if r["ship"].strip() == "1" and r["ship_path"] not in have:
                    r["ship"] = "0"
            write_tsv(ci, cr, list(cr[0].keys()))
            strip_eas_comparators(d)
            if not real:             # tier 1: generated stand-ins; tier 2 (real=True) keeps the REAL built manifest, README and REPRODUCE (V05)
                standins(d, v3, led)
            else:
                restamp_manifest(d, v3)
            if not use_real_tiers():      # RC1b: the baseline too reads the real tiers once they match the real out/ (the stand-in carried no comparator rule for the East Asian scans, S2-D34(6))
                standin_tiers(d)
            ms = tmp / (d.name + "_ms")
            ms.mkdir()
            for f in manuscript_default(v3):
                (ms / f.name).write_text(f.read_text(encoding="utf-8"), encoding="utf-8")
            (ms / "zz_cited.md").write_text("The frame is Figure_3__3a_forest.csv in the repository.\n", encoding="utf-8")
            return d, ms, led

        base, ms, led = fresh("base")
        R0 = run_gate(base, v3, ms)
        print("baseline on the scratch copy (stand-ins written for any file Builder D has not delivered):")
        ok = R0.report(lambda s: print("   " + s))
        if not ok:
            print("SELFTEST FAIL: baseline is not green, mutations cannot be interpreted")
            return 1
        tests = []

        def mutate(name, expect, fn, real=False):
            d, m, l = fresh(name, real)
            kw = fn(d, m, l) or {}
            r = run_gate(d, v3, kw.get("manuscript", m), kw.get("ledger_path"), kw.get("journal"))
            tr = r.tripped()
            good = expect in tr
            nonself = sorted(tr - {expect})
            print("  %-34s expect %-4s -> tripped %s %s" % (name, expect, sorted(tr), "OK" if good else "MISS"))
            tests.append(good)
            if nonself:
                print("      (also tripped %s: the mutation touches more than one invariant)" % nonself)

        def unique(path_glob_name, root):
            hits = [p for p in root.rglob("*") if p.is_file() and p.name == path_glob_name]
            assert len(hits) == 1, "mutation target not unique: %s" % path_glob_name
            return hits[0]

        def m_extra(d, m, l):
            (d / "out" / "extra_file.txt").write_text("not in the manifest\n", encoding="utf-8")

        def m_flip(d, m, l):
            p = unique("priorart_results.tsv", d / "out")
            b = bytearray(p.read_bytes())
            i = next(j for j in range(len(b) // 2, len(b)) if chr(b[j]).isalpha())
            b[i] ^= 0x20
            p.write_bytes(bytes(b))

        def m_ns(d, m, l):
            ns = read_tsv(d / "inputs" / "NOT_SHIPPED.tsv")
            out_md5 = {md5_file(p) for p in (d / "out").rglob("*") if p.is_file()}
            cand = next(r for r in ns if (v3 / r["path"]).is_file() and int(r["bytes"]) < 300000 and r["md5"] not in out_md5
                        and md5_file(v3 / r["path"]) == r["md5"])
            tgt = d / "out" / "_ns_copy.txt"
            shutil.copyfile(v3 / cand["path"], tgt)

            def f(rows):
                rows.append({"path": "_ns_copy.txt", "bytes": str(tgt.stat().st_size), "md5": cand["md5"], "rows": "", "columns": "",
                             "source_path": cand["path"], "source_md5": cand["md5"], "producer_script": "", "producer_md5": "",
                             "software_label": "post hoc", "cited_by": "", "licence_class": "own_derived"})
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", f)

        def m_i8(d, m, l):
            src = unique("a3_umod_q.tsv", d / "out")
            new = src.with_suffix(".i8")
            src.rename(new)
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: [r.__setitem__("path", "kidney/a3_umod_q.i8") for r in rows if r["path"] == "kidney/a3_umod_q.tsv"])

        def m_secret(d, m, l):
            p = unique("VERDICT.md", d / "out")
            p.write_text(p.read_text(encoding="utf-8") + '\napi' + '_key = "ab12cd34ef56ab12cd34ef56ab12cd34"\n', encoding="utf-8")
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: refresh_row(d / "out", rows, "workstream_A/VERDICT.md"))

        def m_jwt(d, m, l):
            p = unique("VERDICT.md", d / "out")
            p.write_text(p.read_text(encoding="utf-8") + "\n" + fake_jwt(sig="abcDEF123_-xyz", payload=b'{"sub":"1234567890"}') + "\n", encoding="utf-8")
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: refresh_row(d / "out", rows, "workstream_A/VERDICT.md"))

        def append_pin(d, text):
            p = unique("VERDICT.md", d / "out")
            p.write_text(p.read_text(encoding="utf-8") + "\n" + text + "\n", encoding="utf-8")
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: refresh_row(d / "out", rows, "workstream_A/VERDICT.md"))

        def one_replace(path, old, new):
            t = path.read_text(encoding="utf-8")
            assert t.count(old) == 1, "mutation target not unique (%d): %s" % (t.count(old), old[:50])
            path.write_text(t.replace(old, new, 1), encoding="utf-8")

        def r_jwt_short(d, m, l):         # tier 2: a JWT whose payload segment is 3 characters long (the verifier's case)
            append_pin(d, "Authorization: " + fake_jwt(sig="xxx", payload=b"{}"))

        def r_pw_bare(d, m, l):           # tier 2: an unquoted password
            append_pin(d, "pass" + "word=" + "hunter2" * 2)

        def r_wrongcount(d, m, l):        # tier 2: a count that is wrong but whose digit occurs elsewhere in the same ledger display (W3 11 -> W3 9)
            one_replace(d / "out" / "README.md", "W3 11", "W3 9")

        def r_label(d, m, l):             # tier 2: a count hidden behind a label word
            one_replace(d / "out" / "README.md", "Rows of Table 2: 8 [L671].", "Rows of Table 2: 8 [L671]. Table 2 covers Stage 41 pairs.")

        def r_class_spaced(d, m, l):      # tier 2: a unit-table class cell in free wording
            one_replace(d / "out" / "REPRODUCE.md", "| kidney | deliverable | does_not_regenerate_engine_output |", "| kidney | deliverable | regenerates from the release |")

        def r_row_missing(d, m, l):       # tier 2: a unit that has no row in the unit table
            t = (d / "out" / "REPRODUCE.md").read_text(encoding="utf-8").splitlines()
            ix = [i for i, ln in enumerate(t) if ln.startswith("| kidney | deliverable |")]
            assert len(ix) == 1
            del t[ix[0]]
            (d / "out" / "REPRODUCE.md").write_text("\n".join(t) + "\n", encoding="utf-8")

        def m_cited(d, m, l):
            p = unique("Figure_3__3a_forest.csv", d / "out")
            p.unlink()
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: rows.__delitem__(next(i for i, r in enumerate(rows) if r["path"].endswith("Figure_3__3a_forest.csv"))))

        def m_typed(d, m, l):
            p = d / "out" / "README.md"
            t = p.read_text(encoding="utf-8")
            assert "7777 tables" not in t
            p.write_text(t + "\nThe release holds 7777 tables.\n", encoding="utf-8")

        def m_class(d, m, l):
            p = d / "out" / "REPRODUCE.md"
            t = p.read_text(encoding="utf-8")
            lines = t.splitlines()
            idx = [i for i, ln in enumerate(lines) if re.search(r"(?<![\w])screen(?![\w])", ln) and "does_not_regenerate_engine_output" in ln.replace(" ", "_")]
            assert len(idx) == 1, "mutation target line not unique (%d)" % len(idx)
            lines[idx[0]] = re.sub(r"does[_ ]not[_ ]regenerate[_ ]engine[_ ]output", "regenerates_from_release", lines[idx[0]])
            p.write_text("\n".join(lines) + "\n", encoding="utf-8")

        def m_bh(d, m, l):
            p = unique("screen_ledger.csv", d / "out")
            t = p.read_text(encoding="utf-8")
            hdr, body = t.split("\n", 1)
            assert body.count(",FALSE,FALSE,", 1) >= 1
            body = body.replace(",FALSE,FALSE,", ",FALSE,FALSE,", 1)
            i = body.index("1.0,FALSE\n")
            body = body[:i] + "1.0,TRUE\n" + body[i + len("1.0,FALSE\n"):]
            p.write_text(hdr + "\n" + body, encoding="utf-8")
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: refresh_row(d / "out", rows, "screen/screen_ledger.csv"))

        def m_cons(d, m, l):
            p = unique("phase3_coloc__coloc_results_final.csv", d / "out")
            t = p.read_text(encoding="utf-8")
            hdr, body = t.split("\n", 1)
            p.write_text(hdr + ",conservative_prior_posterior\n" + body, encoding="utf-8")
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: refresh_row(d / "out", rows, "colocalization/phase3_coloc__coloc_results_final.csv"))

        # ---- RC1 mutation helpers (R6 docx, R8 identity, R9 allowlist, R10 recount, R11 tiers); each targets one invariant of the real release
        def refresh(d, k):
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: refresh_row(d / "out", rows, k))

        def add_to_release(d, k, data, inputs_row=True):
            p = d / "out" / k
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(data)

            def f(rows):
                rows.append({"path": k, "bytes": str(len(data)), "md5": hashlib.md5(data).hexdigest(), "rows": "", "columns": "", "source_path": "", "source_md5": "",
                             "producer_script": "", "producer_md5": "", "software_label": "post hoc", "cited_by": "", "licence_class": "own_derived"})
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", f)
            if inputs_row:
                ri = d / "inputs" / "RELEASE_INPUTS.tsv"
                rows = read_tsv(ri)
                new = {c: "" for c in rows[0]}
                new.update({"input_id": "RI999", "origin": "curated", "role": "copy", "release_path": k, "ship_status": "ship"})
                write_tsv(ri, rows + [new], list(rows[0].keys()))
            tp = d / "build" / "tiers_by_file.tsv"       # a Tier 1 row, so that only the invariant under test trips
            tr = read_tsv(tp)
            write_tsv(tp, tr + [{"file": k, "tier": "1", "providers": "", "decision": "ship", "reason": "test file"}], TIERS_COLS)

        def edit_tsv_first_row(path, col, value):
            rows = list(csv.reader(io.StringIO(path.read_text(encoding="utf-8"), newline=""), delimiter="\t"))
            rows[1][rows[0].index(col)] = value
            with open(path, "w", newline="", encoding="utf-8") as f:
                csv.writer(f, delimiter="\t", lineterminator="\n").writerows(rows)

        def add_column(path, name, value, delim):
            rows = list(csv.reader(io.StringIO(path.read_text(encoding="utf-8"), newline=""), delimiter=delim))
            rows[0].append(name)
            for r in rows[1:]:
                r.append(value)
            with open(path, "w", newline="", encoding="utf-8") as f:
                csv.writer(f, delimiter=delim, lineterminator="\n").writerows(rows)

        def tiers_edit(d, fn):
            tp = d / "build" / "tiers_by_file.tsv"
            rows = read_tsv(tp)
            fn(rows)
            write_tsv(tp, rows, TIERS_COLS)

        def pervar_file(d):
            tr = read_tsv(d / "build" / "tiers_by_file.tsv")
            return next(r["file"] for r in tr if r["tier"] == "2" and r["decision"] == "ship" and r["file"].endswith((".csv", ".tsv")))

        def numeric_col(d, k):
            p_ = d / "out" / k
            rows_ = list(csv.DictReader(io.StringIO(p_.read_text(encoding="utf-8"), newline=""), delimiter="\t" if k.endswith(".tsv") else ","))
            for c in rows_[0]:
                vals = []
                for r_ in rows_:
                    try:
                        vals.append(float(r_[c]))
                    except ValueError:
                        pass
                if len(vals) == len(rows_) and any(abs(v * 2 - round(v * 2)) > 1e-6 for v in vals):
                    return c
            raise AssertionError("no non-grid numeric column in " + k)

        def add_rule(d, f_, handling):
            rr = read_tsv(d / "inputs" / "TIER_RULES.tsv")
            rr = [x for x in rr if not (x["kind"] == "file" and x["target"] == f_)]
            rr.append({"rule_id": "ZZ1", "kind": "file", "target": f_, "tier": "2", "providers": "CKB", "handling": handling, "selector": "", "basis": "test", "reason": "test"})
            write_tsv(d / "inputs" / "TIER_RULES.tsv", rr, ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])

        def docx_variant(tmpdir, fn, name):
            src = v3 / MANUSCRIPT_DOCX
            dst = tmpdir / name
            with zipfile.ZipFile(src) as zi, zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zo:
                for it in zi.infolist():
                    b = zi.read(it.filename)
                    if it.filename == "word/document.xml":
                        b = fn(b.decode("utf-8")).encode("utf-8")
                    zo.writestr(it, b)
            return dst

        def r6_docx_ghost(d, m, l):
            def fn(x):
                i = x.rfind("<w:sectPr")
                assert i > 0
                return x[:i] + "<w:p><w:r><w:t>Values are in ghost_zz_file.tsv.</w:t></w:r></w:p>" + x[i:]
            return {"manuscript": docx_variant(tmp, fn, d.name + "_ghost.docx")}

        def r6_docx_blind(d, m, l):
            def fn(x):
                assert x.count("UKB-PPP") >= 1
                return x.replace("UKB-PPP", "UKB‑PPP")
            return {"manuscript": docx_variant(tmp, fn, d.name + "_blind.docx")}

        def r8_pair(d, m, l):             # an in-band row count, but a pair that is not one of the canonical 81
            edit_tsv_first_row(unique("fig2a_plotted.tsv", d / "out"), "gene", "ZZ9")
            refresh(d, "source_data/figures/Figure_2/fig2a_plotted.tsv")

        def r8_nondecl(d, m, l):          # 15 rows (under the old 81-row bound) in a file that is not a declared carrier
            add_column(unique("Table1_data_sources.tsv", d / "out"), "pp_h4_conservative_x", "0.5", "\t")
            refresh(d, "tables/Table1_data_sources.tsv")

        def r8_t2pair(d, m, l):           # Table 2 names a pair that is neither in the 81 nor printed in the manuscript's Table 2
            one_replace(unique("Table2_evidence_cards.tsv", d / "out"), "PCSK9-CAD\tPCSK9\tCAD\t", "PCSK9-CAD\tZZ9\tCAD\t")
            refresh(d, "tables/Table2_evidence_cards.tsv")

        def edit_table2_cell(d, pair_key, col, value):     # one cell of the shipped Table 2 (by pair_key), manifest refreshed
            path = unique("Table2_evidence_cards.tsv", d / "out")
            rows = list(csv.reader(io.StringIO(path.read_text(encoding="utf-8"), newline=""), delimiter="\t"))
            hits = [r for r in rows[1:] if r[rows[0].index("pair_key")] == pair_key]
            assert len(hits) == 1, "mutation target not unique: %s" % pair_key
            hits[0][rows[0].index(col)] = value
            with open(path, "w", newline="", encoding="utf-8") as f:
                csv.writer(f, delimiter="\t", lineterminator="\n").writerows(rows)
            refresh(d, "tables/Table2_evidence_cards.tsv")

        def r8_declcol(d, m, l):          # a declared pair (FGFR4-T2D) holds a conservative-prior value in a column the declaration does not permit
            edit_table2_cell(d, "FGFR4-T2D (EAS)", "eas_pp_h4_conservative", "0.5")

        def r8_decl_empty(d, m, l):       # a declared pair whose declaration permits no column (DPEP1-eGFR) holds a value
            edit_table2_cell(d, "DPEP1-eGFR", "pp_h4_conservative", "0.5")

        def r9_unknown(d, m, l):
            add_to_release(d, "zz_unknown/extra.txt", b"not an approved inclusion\n", inputs_row=False)

        def zip_bytes(members):
            bio = io.BytesIO()
            with zipfile.ZipFile(bio, "w") as z:
                for n, b in members:
                    z.writestr(n, b)
            return bio.getvalue()

        def r9_zip(d, m, l):              # the archive itself IS approved; its member is not declared
            add_to_release(d, "zz_arch/data.zip", zip_bytes([("payload.csv", b"a,b\n1,2\n")]))

        def r9_nested(d, m, l):           # approved archive -> declared member (an inner archive) -> undeclared file inside it
            inner = zip_bytes([("deep.txt", b"hidden\n")])
            add_to_release(d, "zz_arch/outer.zip", zip_bytes([("inner.zip", inner)]))
            with open(d / "_gate_fixtures" / "R9_ARCHIVE_MEMBERS.tsv", "a", encoding="utf-8", newline="") as f:
                f.write("zz_arch/outer.zip::inner.zip\tdeclared for the test\n")

        def r9_withheld(d, m, l):         # the file is in out/, but its only inclusion row says withheld
            ri = d / "inputs" / "RELEASE_INPUTS.tsv"
            rows = read_tsv(ri)
            row = next(r for r in rows if r["release_path"] == "literature/pairs81_priorart.tsv" and r["ship_status"] == "ship")
            row["ship_status"] = "withheld_pending_terms"
            write_tsv(ri, rows, list(rows[0].keys()))

        def r10_drop_provider(d, m, l):   # PROVIDERS.tsv loses a row: 35 -> 34, neither the ledger (36) nor the declared recount (35)
            pv = d / "inputs" / "PROVIDERS.tsv"
            rows = read_tsv(pv)
            write_tsv(pv, rows[:-1], list(rows[0].keys()))

        def r10_stale_exception(d, m, l): # S6.3b2a: PROVIDERS.tsv gains a row: the live list (36) no longer equals L533 (35); with the exception retired this is plain drift
            pv = d / "inputs" / "PROVIDERS.tsv"
            rows = read_tsv(pv)
            new = dict(rows[-1])
            new["provider_id"] = "ZZ-STALE-TEST"
            new["in_declarations"] = "no"
            write_tsv(pv, rows + [new], list(rows[0].keys()))

        def ledger_edit(row_id, was, now):    # a copy of the ledger with one row value changed (the real row must carry the value it is said to have)
            src = v3 / LEDGER_REL
            rows = list(csv.reader(io.StringIO(src.read_text(encoding="utf-8"), newline="")))
            hdr = rows[0]
            hit = [r for r in rows if r[0] == row_id]
            assert len(hit) == 1 and hit[0][hdr.index("value")] == was, "R10 mutation target %s is not the expected value %s" % (row_id, was)
            hit[0][hdr.index("value")] = now
            p_ = tmp / ("%s_%s_ledger.csv" % (row_id, now.replace("[", "").replace("]", "").replace(", ", "_")))
            with open(p_, "w", newline="", encoding="utf-8") as f:
                csv.writer(f, lineterminator="\n").writerows(rows)
            return {"ledger_path": p_}

        def r10_l676_off(d, m, l):        # L676 no longer equals its frozen input
            return ledger_edit("L676", "[36, 14]", "[37, 14]")

        def r10_l533_off(d, m, l):        # L533 no longer equals the live provider list
            return ledger_edit("L533", "[35, 14]", "[36, 14]")

        def r10_ledger_row(d, m, l):      # a ledger row (L663 SuSiE grid rows) that no longer equals its input
            src = v3 / LEDGER_REL
            rows = list(csv.reader(io.StringIO(src.read_text(encoding="utf-8"), newline="")))
            hdr = rows[0]
            hit = [r for r in rows if r[0] == "L663"]
            assert len(hit) == 1 and hit[0][hdr.index("value")] == "256"
            hit[0][hdr.index("value")], hit[0][hdr.index("display")] = "257", "257"
            p_ = tmp / (d.name + "_ledger.csv")
            with open(p_, "w", newline="", encoding="utf-8") as f:
                csv.writer(f, lineterminator="\n").writerows(rows)
            return {"ledger_path": p_}

        def r15_label_swap(d, m, l):        # CODE_INCLUDE says "full run-time hash match"; the manifest row of that file says "timestamp evidence only"
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: next(r for r in rows if r.get("provenance_label") == "full run-time hash match").__setitem__("provenance_label", "timestamp evidence only"))

        def r15_evidence_md5(d, m, l):      # a "full run-time hash match" whose recorded md5 is not the original's
            ci = d / "inputs" / "CODE_INCLUDE.tsv"
            rows = read_tsv(ci)
            hit = [r for r in rows if r.get("provenance_label") == "full run-time hash match"]
            assert hit, "R15 mutation target: no full-match row in the real CODE_INCLUDE"
            hit[0]["evidence_md5"] = ("0" if hit[0]["evidence_md5"][0] != "0" else "1") + hit[0]["evidence_md5"][1:]
            write_tsv(ci, rows, list(rows[0].keys()))

        def r15_transform_claim(d, m, l):   # a portable copy that differs from its original but is declared byte-identical
            shipped_ci = {r["ship_path"] for r in read_tsv(d / "inputs" / "CODE_INCLUDE.tsv") if r["ship"] == "1"}      # a real bundle row (the records under code/ are not CODE_INCLUDE rows)
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: next(r for r in rows if r["path"] in shipped_ci and r["md5"] != r["source_md5"]).__setitem__("transform", "byte-identical to the original"))

        def r15_notice(d, m, l):            # a third-party notice in a shipped script whose code_rights says "authors"
            p = unique("33_sensitivity_provenance.R", d / "out")
            p.write_text(p.read_text(encoding="utf-8") + "\n# Copy" + "right (" + "c) 20" + "19 Some Author. All rights re" + "served.\n", encoding="utf-8")
            man_rewrite(d / "out" / "RELEASE_MANIFEST.tsv", lambda rows: refresh_row(d / "out", rows, "code/paper3_deposit_request_v5/scripts/scripts/paper3.1_ckm/phase3/33_sensitivity_provenance.R"))

        def r15_readme_l676(d, m, l):       # the README prints the provider count from the frozen row
            p = d / "out" / "README.md"
            p.write_text(p.read_text(encoding="utf-8") + "\nThe list has 36 providers [L676].\n", encoding="utf-8")

        def r15_unlabelled(d, m, l):        # a deposit-package script loses its label in CODE_INCLUDE
            ci = d / "inputs" / "CODE_INCLUDE.tsv"
            rows = read_tsv(ci)
            hit = [r for r in rows if r["path"].startswith(DEPOSIT_PREFIXES) and r["ship"] == "1" and r.get("provenance_label")]
            assert hit, "R15 mutation target: no labelled deposit-package row in the real CODE_INCLUDE"
            hit[0]["provenance_label"] = ""
            write_tsv(ci, rows, list(rows[0].keys()))

        def r11_norow(d, m, l):
            f_ = pervar_file(d)
            tiers_edit(d, lambda rows: rows.remove(next(r for r in rows if r["file"] == f_)))

        def r11_tier1_pervar(d, m, l):
            f_ = pervar_file(d)
            tiers_edit(d, lambda rows: next(r for r in rows if r["file"] == f_).update({"tier": "1"}))

        def r11_badprov(d, m, l):         # Tier 2 ships for a provider that no official page confirmed CC0 / CC BY
            f_ = pervar_file(d)
            tiers_edit(d, lambda rows: next(r for r in rows if r["file"] == f_).update({"providers": "FinnGen"}))

        def r11_tier3(d, m, l):           # S2-D34(1): the TPMI switch is 'ship' in the real release, so Tier 3 may ship; with the switch put back to 'withheld' it must not
            f_ = pervar_file(d)
            tiers_edit(d, lambda rows: next(r for r in rows if r["file"] == f_).update({"tier": "3", "providers": "TPMI"}))
            rr = read_tsv(d / "inputs" / "TIER_RULES.tsv")
            sw_ = [x for x in rr if x["kind"] == "switch" and x["target"] == "TPMI_DERIVED"]
            assert len(sw_) == 1, "mutation target not unique: TPMI switch"
            sw_[0]["handling"] = "withheld"
            write_tsv(d / "inputs" / "TIER_RULES.tsv", rr, ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])

        def r11_withheld_in_out(d, m, l):
            f_ = pervar_file(d)
            tiers_edit(d, lambda rows: next(r for r in rows if r["file"] == f_).update({"decision": "withheld"}))

        def r11_binned_pervar(d, m, l):
            # a 'binned' decision with no binned rule behind it: the file is a shipped table that no TIER_RULES file row covers (every Tier 2 file of the real release now has a rule)
            ruled = {x["target"] for x in read_tsv(d / "inputs" / "TIER_RULES.tsv") if x["kind"] == "file"}
            f_ = next(r["file"] for r in read_tsv(d / "build" / "tiers_by_file.tsv")
                      if r["decision"] == "ship" and r["file"].endswith((".csv", ".tsv")) and r["file"] not in ruled and r["file"].startswith(("screen/", "source_data/")))
            tiers_edit(d, lambda rows: next(r for r in rows if r["file"] == f_).update({"decision": "binned"}))

        def r11_eas_cols(d, m, l):        # S2-D34(6): the DIAMANTE frequency eur_af, which the T2D rule withholds and does not declare derived, back in the T2D scan
            add_column(unique("phase2_T2D_scan.tsv", d / "out"), "eur_af", "0.5", "\t")
            refresh(d, "east_asian/phase2_T2D_scan.tsv")

        def r11_eas_comp(d, m, l):        # S2-D34(6): the CAD comparator accession is swapped for one that is not cc0 (GCST010118, EMBL-EBI terms)
            rr = read_tsv(d / "inputs" / "TIER_RULES.tsv")
            hit = [x for x in rr if x["kind"] == "file" and x["target"] == "east_asian/phase2_CAD_scan.tsv"]
            assert len(hit) == 1 and hit[0]["selector"] == "comparator:GCST90132314", "mutation target not unique / changed"
            hit[0]["selector"] = "comparator:GCST010118"
            write_tsv(d / "inputs" / "TIER_RULES.tsv", rr, ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])

        def r11_eas_nowald(d, m, l):      # S2-D34(2): the Tier 1 principle behind the derived T2D columns is removed
            rr = read_tsv(d / "inputs" / "TIER_RULES.tsv")
            assert sum(1 for x in rr if x["target"] == "SINGLE_INSTRUMENT_WALD") == 1, "mutation target not unique"
            write_tsv(d / "inputs" / "TIER_RULES.tsv", [x for x in rr if x["target"] != "SINGLE_INSTRUMENT_WALD"], ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])

        def r11_rule_col(d, m, l):        # a rule withholds a column that still ships
            f_ = pervar_file(d)
            add_rule(d, f_, "columns:" + numeric_col(d, f_))

        def r11_rule_bin(d, m, l):        # a rule bins a column whose values are not on the grid
            f_ = pervar_file(d)
            c_ = numeric_col(d, f_)
            add_rule(d, f_, "binned:%s:0.5" % c_)
            tiers_edit(d, lambda rows: next(r for r in rows if r["file"] == f_).update({"decision": "binned_x"}))

        def r11_rule_rows(d, m, l):       # a rows rule selects rows and none is masked
            f_ = pervar_file(d)
            first = (d / "out" / f_).read_text(encoding="utf-8").splitlines()[0].split("\t" if f_.endswith(".tsv") else ",")[0].strip('"')
            add_rule(d, f_, "rows:%s~.:%s" % (first, numeric_col(d, f_)))
            tiers_edit(d, lambda rows: next(r for r in rows if r["file"] == f_).update({"decision": "reduced_rows"}))

        def r11_noswitch(d, m, l):
            rr = [x for x in read_tsv(d / "inputs" / "TIER_RULES.tsv") if x["target"] != "TPMI_DERIVED"] + [
                {"rule_id": "ZZ0", "kind": "switch", "target": "OTHER", "tier": "", "providers": "", "handling": "withheld", "selector": "", "basis": "", "reason": ""}]
            write_tsv(d / "inputs" / "TIER_RULES.tsv", rr, ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])

        def r11_norules(d, m, l):
            (d / "inputs" / "TIER_RULES.tsv").unlink()

        print("mutations (each in its own scratch copy):")
        mutate("R1 extra file", "R1", m_extra)
        mutate("R2 byte flip", "R2", m_flip)
        mutate("R3 not-shipped copy", "R3", m_ns)
        mutate("R3b class mismatch", "R3b", m_class)
        mutate("R4 renamed .i8", "R4", m_i8)
        def r14_pick(d, l):                # a withheld value that appears nowhere in out/ and is not a declared coincidence; read from withheld/ at run time
            vals = leak_values(d, l)
            allow_k = {leak_key(r["value"]) for r in read_tsv(d / LEAK_ALLOW)}
            seen_k = set()
            for p_ in (d / "out").rglob("*"):
                t_ = read_text(p_) if p_.is_file() else None
                if t_:
                    seen_k |= {leak_key(x) for x in LEAK_FL.findall(t_) if leak_sig(x) >= LEAK_MIN_SIG}
            return next(v[1] for k_, v in sorted(vals.items()) if k_ not in seen_k and k_ not in allow_k)

        def r14_code_comment(d, m, l):     # RC1b fix Q01: a withheld value in a comment of a shipped code file (the leak that sat in rel_tables.py)
            f_ = unique("rel_tables.py", d / "out")
            f_.write_text(f_.read_text(encoding="utf-8") + "\n# SW1 planted " + r14_pick(d, l) + "\n", encoding="utf-8")
            refresh(d, f_.relative_to(d / "out").as_posix())

        def r14_data_cell(d, m, l):        # an undeclared withheld value in a data cell of a shipped table
            f_ = unique("S1R1_decoy_rates.tsv", d / "out")
            f_.write_text(f_.read_text(encoding="utf-8") + "zz\t" + r14_pick(d, l) + "\n", encoding="utf-8")
            refresh(d, f_.relative_to(d / "out").as_posix())

        def r14_pick6(d, l):               # RC2: a withheld value of >= 7 digits that no shipped token matches at 6 digits and no declared coincidence covers
            vals = leak_values(d, l)
            allow6 = {(r["out_path"], leak_round(r["value"], LEAK_SHIP_MIN_SIG)) for r in read_tsv(d / LEAK_ALLOW) if leak_sig(r["value"]) >= LEAK_SHIP_MIN_SIG}
            seen6 = set()
            for p_ in (d / "out").rglob("*"):
                t_ = read_text(p_) if p_.is_file() else None
                if t_:
                    seen6 |= {leak_round(x, LEAK_SHIP_MIN_SIG) for x in LEAK_FL.findall(t_) if leak_sig(x) >= LEAK_SHIP_MIN_SIG}
            allow_any = {a_[1] for a_ in allow6}
            return next("%.6g" % float(v[1]) for k_, v in sorted(vals.items()) if leak_sig(v[1]) >= 7 and leak_round(v[1], LEAK_SHIP_MIN_SIG) not in seen6
                        and leak_round(v[1], LEAK_SHIP_MIN_SIG) not in allow_any and "%.6g" % float(v[1]) != v[1] and float("%.6g" % float(v[1])) != float(v[1]))

        def r14_notes_cell_6(d, m, l):      # RC2 (S2-D35(3)): a withheld value re-printed at 6 significant digits in a notes cell of a shipped table
            f_ = unique("S1R1_decoy_rates.tsv", d / "out")
            f_.write_text(f_.read_text(encoding="utf-8") + "zz\tfrom b/se; file p_value " + r14_pick6(d, l) + "\n", encoding="utf-8")
            refresh(d, f_.relative_to(d / "out").as_posix())

        def r14_ledger_value_6(d, m, l):    # RC2: the same re-print in a value cell of the public ledger
            f_ = d / "out" / "ledger" / "claims_ledger_public.tsv"
            x_ = r14_pick6(d, l)
            f_.write_text(f_.read_text(encoding="utf-8") + "\t".join(["L999", "s", "planted", x_, x_, "src", "not separately checked", ""]) + "\n", encoding="utf-8")
            refresh(d, "ledger/claims_ledger_public.tsv")

        def r14_row_local_4(d, m, l):       # RC2: inside one row, a withheld cell's value re-printed at 4 digits in a notes cell of the shipped row
            out_, wh_ = d / "out", d / "withheld"
            for w in sorted(wh_.rglob("*.tsv")):
                k_ = w.relative_to(wh_).as_posix()
                if not (out_ / k_).is_file():
                    continue
                W_, O_ = _grid(w), _grid(out_ / k_)
                if len(W_) != len(O_) or len(W_) < 3:
                    continue
                for i_ in range(1, len(W_)):
                    for j_, c_ in enumerate(W_[i_]):
                        if W_[0][j_] in O_[0] and O_[0].index(W_[0][j_]) < len(O_[i_]) and O_[i_][O_[0].index(W_[0][j_])] == c_:
                            continue
                        toks_ = [x for x in LEAK_FL.findall(c_) if leak_sig(x) >= 7]
                        if toks_:
                            O_[i_][-1] = O_[i_][-1] + " (file value " + "%.4g" % float(toks_[0]) + ")"
                            buf_ = io.StringIO()
                            csv.writer(buf_, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None, quotechar=None).writerows(O_)
                            (out_ / k_).write_text(buf_.getvalue(), encoding="utf-8")
                            refresh(d, k_)
                            return
            raise AssertionError("no row-local mutation target")

        mutate("R14 withheld value in a code comment", "R14", r14_code_comment, real=True)
        mutate("R14 6-digit re-print in a notes cell", "R14", r14_notes_cell_6, real=True)
        mutate("R14 6-digit re-print in a ledger value", "R14", r14_ledger_value_6, real=True)
        mutate("R14 4-digit re-print in a row", "R14", r14_row_local_4, real=True)
        mutate("R14 withheld value in a data cell", "R14", r14_data_cell, real=True)
        mutate("R5 planted api key", "R5", m_secret)
        mutate("R5 planted JWT", "R5", m_jwt)
        mutate("R6 removed cited file", "R6", m_cited)
        mutate("R7 typed number in README", "R7", m_typed)
        mutate("R8 flipped BH flag", "R8", m_bh)
        mutate("R8 conservative column", "R8", m_cons)
        # tier 2: the same gate against a copy of the REAL built release (real manifest, README, REPRODUCE)
        bad_real, mR, _ = fresh("base_real", real=True)
        RR = run_gate(bad_real, v3, mR)
        print("tier 2 baseline on a copy of the real built release:")
        if not RR.report(lambda s_: print("   " + s_)):
            print("SELFTEST FAIL: the real-release baseline is not green (rebuild the release first), so tier-2 mutations cannot be interpreted")
            return 1
        print("tier 2 mutations (copy of the real built release):")
        mutate("T2 R5 short-payload JWT", "R5", r_jwt_short, real=True)
        mutate("T2 R5 unquoted password", "R5", r_pw_bare, real=True)
        mutate("T2 R7 wrong-but-present count", "R7", r_wrongcount, real=True)
        mutate("T2 R7 count behind a label", "R7", r_label, real=True)
        mutate("T2 R3b spaced class wording", "R3b", r_class_spaced, real=True)
        mutate("T2 R3b unit row missing", "R3b", r_row_missing, real=True)
        def r12_exec(d, m, l):            # an absolute path in a non-comment line of a shipped R script outside code/
            f_ = unique("S1R1_decoy_rates.R", d / "out")
            f_.write_text(f_.read_text(encoding="utf-8") + 'V3 <- "' + "C" + ':/Users/' + 'example/Downloads/x"\n', encoding="utf-8")
            refresh(d, f_.relative_to(d / "out").as_posix())

        def r12_comment_only(d, m, l):    # control for the branch: a path in a COMMENT line of an executable must not trip R12 (it trips R13 only)
            f_ = unique("S1R1_decoy_rates.R", d / "out")
            f_.write_text(f_.read_text(encoding="utf-8") + '# ' + "D" + ':/data/x\n', encoding="utf-8")
            refresh(d, f_.relative_to(d / "out").as_posix())

        def r13_cell(d, m, l):            # a personal path in a data cell of a shipped table
            f_ = unique("results.tsv", d / "out" / "workstreams" / "W2")
            f_.write_text(f_.read_text(encoding="utf-8") + "S2.zz\t" + "C" + ":/Users/" + "example/Downloads/x.tsv\n", encoding="utf-8")
            refresh(d, f_.relative_to(d / "out").as_posix())

        def r13_md(d, m, l):              # a drive path in a shipped document (REPRODUCE.md): the check is not limited to tables or code
            f_ = d / "out" / "REPRODUCE.md"
            f_.write_text(f_.read_text(encoding="utf-8") + "\n  - input: " + "D" + ":/data/x.gz\n", encoding="utf-8")

        print("tier 2 mutations for the RC1 checks (R6 docx, R8 identity, R9, R10, R11):")
        mutate("T2 R6 docx names a ghost file", "R6", r6_docx_ghost, real=True)
        mutate("T2 R6 docx reader blind", "R6", r6_docx_blind, real=True)

        def r6_stale_unit():
            # Stage 6.3b part 2b. R6's staleness branch runs only on the real default docx, which mutate() never reaches, so it is tested here
            # on a scratch tree: a source of the docx newer than the docx must be reported; a newer 'generated' row (another file's builder) must not.
            t = tmp / "r6_stale_unit" / "R2_rewrite"
            for sub in ("build", "sections_R2", "scripts"):
                (t / sub).mkdir(parents=True)
            src, gen, doc = t / "sections_R2" / "m.md", t / "scripts" / "b.py", t / "build" / "ms.docx"
            for f_ in (src, gen, doc):
                f_.write_text("x\n", encoding="utf-8")
            (t / "build" / "P31R_SOURCE_MAP.tsv").write_text(
                "component\tcanonical_file\tstate\tconsumers\tevidence\n"
                "methods\tsections_R2/m.md\tR2_operative\t\t\n"
                "af1\tscripts/b.py\tgenerated\t\t\n", encoding="utf-8")
            fs_ = r6_docx_sources(t.parent)
            os.utime(src, (1000, 1000)); os.utime(doc, (2000, 2000)); os.utime(gen, (3000, 3000))
            ok_gen = [f_.name for f_ in fs_] == ["m.md"] and r6_stale(doc, fs_) is None
            os.utime(src, (4000, 4000))
            ok_src = r6_stale(doc, fs_) == src
            print("  %-34s expect none -> %s" % ("R6 newer generated row (control)", "OK" if ok_gen else "MISS"))
            print("  %-34s expect R6   -> %s" % ("R6 newer docx source (stale)", "OK" if ok_src else "MISS"))
            tests.extend([ok_gen, ok_src])
        r6_stale_unit()

        # R19 (S2-D43): each rung edits a COPY of the real built Additional file 1 (or the scratch declarations), never the build folder
        def r19_af1_copy(name, edit):
            import openpyxl
            src_ = v3 / JOURNAL["AF1"]
            dst_ = tmp / ("r19_" + name + ".xlsx")
            wb_ = openpyxl.load_workbook(src_)
            edit(wb_)
            wb_.save(dst_)
            return {"journal": {"AF1": str(dst_)}}

        def r19_plant(d, m, l):             # a withheld value the release does not ship, typed into a legend cell of S2
            tok = r14_pick(d, l)

            def ed(wb_):
                ws_ = wb_["S2"]
                ws_.cell(row=2, column=1).value = str(ws_.cell(row=2, column=1).value) + " Planted " + tok + "."
            return r19_af1_copy("plant", ed)

        def r19_h9(d, m, l):                # a hold-H9 pair's estimate restored as a number
            def ed(wb_):
                ws_ = wb_["S7a"]
                hit_ = [c_ for row_ in ws_.iter_rows() for c_ in row_ if c_.value == "withheld (Additional file 2, Section S-R13)"]
                assert hit_, "mutation target absent: no hold-H9 marker in S7a"
                hit_[0].value = -0.0144
            return r19_af1_copy("h9", ed)

        def r19_undeclare(d, m, l):         # the S3 'Instrument F' declaration removed: the 81 UKB-PPP F values must then fail
            p_ = d / R19_DECL
            rows_ = p_.read_text(encoding="utf-8").splitlines()
            keep_ = [x for x in rows_ if "\tInstrument F\t" not in x]
            assert len(keep_) == len(rows_) - 1, "mutation target not unique: the Instrument F declaration"
            p_.write_text("\n".join(keep_) + "\n", encoding="utf-8")

        def r19_blind(d, m, l):             # a workbook with nothing in it: the reader floor must refuse it
            import openpyxl
            dst_ = tmp / "r19_blind.xlsx"
            wb_ = openpyxl.Workbook()
            wb_.active.title = "S1"
            wb_.save(dst_)
            return {"journal": {"AF1": str(dst_)}}
        print("tier 2 mutations for R19 (the journal files):")
        mutate("T2 R19 withheld value in AF1", "R19", r19_plant, real=True)
        mutate("T2 R19 hold-H9 magnitude in S7a", "R19", r19_h9, real=True)
        mutate("T2 R19 S3 declaration removed", "R19", r19_undeclare, real=True)
        mutate("T2 R19 blind workbook reader", "R19", r19_blind, real=True)
        mutate("T2 R8 non-canonical pair (fig2a)", "R8", r8_pair, real=True)
        mutate("T2 R8 undeclared carrier, 15 rows", "R8", r8_nondecl, real=True)
        mutate("T2 R8 Table 2 foreign pair", "R8", r8_t2pair, real=True)
        mutate("T2 R8 declared pair, undeclared column", "R8", r8_declcol, real=True)
        mutate("T2 R8 declared pair with no permitted column", "R8", r8_decl_empty, real=True)
        mutate("T2 R9 unknown file", "R9", r9_unknown, real=True)
        mutate("T2 R9 undeclared archive member", "R9", r9_zip, real=True)
        mutate("T2 R9 member of nested archive", "R9", r9_nested, real=True)
        mutate("T2 R9 only a withheld row", "R9", r9_withheld, real=True)
        mutate("T2 R10 recount drift (provider)", "R10", r10_drop_provider, real=True)
        mutate("T2 R10 live list drifts from L533", "R10", r10_stale_exception, real=True)
        mutate("T2 R10 L676 off its frozen input", "R10", r10_l676_off, real=True)
        mutate("T2 R10 L533 off the live list", "R10", r10_l533_off, real=True)
        mutate("T2 R10 ledger row off its input", "R10", r10_ledger_row, real=True)
        mutate("T2 R15 manifest label differs from CODE_INCLUDE", "R15", r15_label_swap, real=True)
        mutate("T2 R15 full match with a wrong recorded md5", "R15", r15_evidence_md5, real=True)
        mutate("T2 R15 changed copy declared byte-identical", "R15", r15_transform_claim, real=True)
        mutate("T2 R15 third-party notice in an 'authors' script", "R15", r15_notice, real=True)
        mutate("T2 R15 README cites the frozen L676", "R15", r15_readme_l676, real=True)
        mutate("T2 R15 deposit script without a label", "R15", r15_unlabelled, real=True)
        mutate("T2 R11 shipped file, no tier row", "R11", r11_norow, real=True)
        mutate("T2 R11 per-variant file as Tier 1", "R11", r11_tier1_pervar, real=True)
        mutate("T2 R11 Tier 2 non-CC provider", "R11", r11_badprov, real=True)
        def r11_s2d35(d, m, l):            # RC2: the CLSA W4a copy ships whole only under the declared S2-D35 file rule; without the basis it is Tier 2 from an unconfirmed provider
            rr = read_tsv(d / "inputs" / "TIER_RULES.tsv")
            tg_ = [x for x in rr if x["kind"] == "file" and x["target"] == "workstreams/W4a/results.tsv"]
            assert len(tg_) == 1, "mutation target not unique: W4a rule"
            assert "S2-D35" in tg_[0]["basis"], "mutation target does not carry the S2-D35 basis"
            tg_[0]["basis"] = tg_[0]["basis"].replace("S2-D35", "S2-D00")
            write_tsv(d / "inputs" / "TIER_RULES.tsv", rr, ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])

        mutate("T2 R11 CLSA ships without the S2-D35 basis", "R11", r11_s2d35, real=True)
        mutate("T2 R11 Tier 3 ships", "R11", r11_tier3, real=True)
        mutate("T2 R11 withheld decision in out/", "R11", r11_withheld_in_out, real=True)
        mutate("T2 R11 binned file keeps rsid", "R11", r11_binned_pervar, real=True)
        mutate("T2 R11 eur_ column, comparators", "R11", r11_eas_cols, real=True)
        mutate("T2 R11 CAD comparator not cc0", "R11", r11_eas_comp, real=True)
        mutate("T2 R11 derived columns without the Wald principle", "R11", r11_eas_nowald, real=True)
        mutate("T2 R11 TIER_RULES absent", "R11", r11_norules, real=True)
        mutate("T2 R11 rule: column still ships", "R11", r11_rule_col, real=True)
        mutate("T2 R11 rule: values off the grid", "R11", r11_rule_bin, real=True)
        mutate("T2 R11 rule: selected rows unmasked", "R11", r11_rule_rows, real=True)
        mutate("T2 R11 TPMI switch undeclared", "R11", r11_noswitch, real=True)

        def w1_restore_f(d, rule_too):      # S2-D40 W1: the instruments file with F_stat put back from the withheld copy; rule_too=False also drops F_stat from the F01 handling (the defect as it was: a rule table that omits it)
            f_ = d / "out" / "screen" / "cis_pqtl_instruments.csv"
            W_, O_ = _grid(d / "withheld" / "screen" / "cis_pqtl_instruments.csv"), _grid(f_)
            assert len(W_) == len(O_) and "F_stat" in W_[0] and "F_stat" not in O_[0], "W1 mutation target is not the expected pair"
            jf = W_[0].index("F_stat")
            O_[0].append("F_stat")
            for i_ in range(1, len(O_)):
                O_[i_].append(W_[i_][jf])
            buf_ = io.StringIO()
            csv.writer(buf_, lineterminator="\n").writerows(O_)
            f_.write_text(buf_.getvalue(), encoding="utf-8")
            refresh(d, "screen/cis_pqtl_instruments.csv")
            if not rule_too:
                rr_ = read_tsv(d / "inputs" / "TIER_RULES.tsv")
                t_ = [x for x in rr_ if x["rule_id"] == "F01"]
                assert len(t_) == 1 and t_[0]["handling"].count("|F_stat|") == 1, "mutation target: F01 does not list F_stat"
                t_[0]["handling"] = t_[0]["handling"].replace("|F_stat|", "|")
                write_tsv(d / "inputs" / "TIER_RULES.tsv", rr_, ['rule_id', 'kind', 'target', 'tier', 'providers', 'handling', 'selector', 'basis', 'reason'])

        def r16_rule_omits_f(d, m, l):      # the real defect: F_stat ships AND the rule table does not list it, so R11 has nothing to say
            w1_restore_f(d, False)

        def r11_f_still_ships(d, m, l):     # F_stat ships although F01 withholds it
            w1_restore_f(d, True)

        def r17_notes_number(d, m, l):      # a masked row whose notes cell carries a number
            f_ = d / "out" / "external_outcomes" / "extchecks" / "results.tsv"
            g_ = _grid(f_)
            hit_ = [i_ for i_, r_ in enumerate(g_) if i_ and "value" in g_[0] and r_[g_[0].index("value")] == MASK_TOKEN and r_[g_[0].index("notes")] == MASK_TOKEN]
            assert hit_, "mutation target: no masked row with masked notes in extchecks/results.tsv"
            g_[hit_[0]][g_[0].index("notes")] = "EA C NEA G EAF 0.95 Beta -0.15 SE 0.019 P 4.6e-15"
            buf_ = io.StringIO()
            csv.writer(buf_, delimiter="\t", lineterminator="\n").writerows(g_)
            f_.write_text(buf_.getvalue(), encoding="utf-8")
            refresh(d, "external_outcomes/extchecks/results.tsv")

        mutate("T2 R16 F_stat ships, rule table omits it", "R16", r16_rule_omits_f, real=True)
        mutate("T2 R11 F_stat ships against F01", "R11", r11_f_still_ships, real=True)
        mutate("T2 R17 masked row, notes quote values", "R17", r17_notes_number, real=True)

        # ---- S2-D41 (W3 deCODE F, W4 provider marginal values in notes): each case works on a copy of the real artifact
        def dec_write(d, mut):
            f_ = d / "out" / "replication" / "decode_replication__decode_mr_results.csv"
            W_, O_ = _grid(d / "withheld" / "replication" / "decode_replication__decode_mr_results.csv"), _grid(f_)
            assert len(W_) == len(O_) and W_[0] == O_[0] and "decode_F_stat" in W_[0], "S2-D41 mutation target is not the expected pair"
            assert any(O_[i_][O_[0].index("decode_F_stat")] == MASK_TOKEN for i_ in range(1, len(O_))), "mutation target: decode_F_stat is not masked in the shipped copy"
            mut(W_, O_)
            buf_ = io.StringIO()
            csv.writer(buf_, lineterminator="\n").writerows(O_)
            f_.write_text(buf_.getvalue(), encoding="utf-8")
            refresh(d, "replication/decode_replication__decode_mr_results.csv")

        def s41_f_back(d, m, l):            # decode_F_stat put back from the withheld copy (rule F15 still lists it)
            def mut(W_, O_):
                j_ = W_[0].index("decode_F_stat")
                for i_ in range(1, len(O_)):
                    O_[i_][j_] = W_[i_][j_]
            dec_write(d, mut)

        def s41_z_col(d, m, l):             # an exposure |z| = sqrt(F) column computed from the withheld F
            def mut(W_, O_):
                j_ = W_[0].index("decode_F_stat")
                O_[0].append("decode_exposure_z")
                for i_ in range(1, len(O_)):
                    v_ = W_[i_][j_]
                    O_[i_].append(("%.7g" % math.sqrt(float(v_))) if v_ not in ("", "NA") else "")
            dec_write(d, mut)

        def s41_delta_back(d, m, l):        # the delta-method SE and p put back: F = decode_b^2 / (delta_se^2 - decode_se^2) is recoverable again
            def mut(W_, O_):
                for c_ in ("delta_se", "delta_pval"):
                    j_ = W_[0].index(c_)
                    for i_ in range(1, len(O_)):
                        O_[i_][j_] = W_[i_][j_]
            dec_write(d, mut)

        def w2_note_set(d, pick, text_fn):
            f_ = d / "out" / "workstreams" / "W2" / "results.tsv"
            g_, w_ = _grid(f_), _grid(d / "withheld" / "workstreams" / "W2" / "results.tsv")
            assert len(g_) == len(w_) and g_[0] == w_[0], "S2-D41 mutation target: the W2 shipped / withheld pair"
            ii_, ni_ = g_[0].index("quantity_id"), g_[0].index("notes")
            hit_ = [i_ for i_ in range(1, len(g_)) if pick(g_[i_][ii_], g_[i_][ni_])]
            assert hit_, "S2-D41 mutation target: no redacted S2.marg.ratio note in the W2 copy"
            g_[hit_[0]][ni_] = text_fn(w_, hit_[0])
            buf_ = io.StringIO()
            csv.writer(buf_, delimiter="\t", lineterminator="\n").writerows(g_)
            f_.write_text(buf_.getvalue(), encoding="utf-8")
            refresh(d, "workstreams/W2/results.tsv")

        def s41_note_quote(d, m, l):        # the redacted marginal z of an S2.marg.ratio.* row put back from the withheld copy: an UNMASKED row quoting provider values
            w2_note_set(d, lambda id_, now_: id_.startswith("S2.marg.ratio.") and "withheld_pending_terms" in now_, lambda w_, r_: w_[r_][w_[0].index("notes")])

        def s41_note_z(d, m, l):            # a note that prints, at 6 digits, the |z| = |value / se| of a masked (value, se) pair: not itself a withheld cell value
            W_, O_ = _grid(d / "withheld" / "workstreams" / "W2" / "results.tsv"), _grid(d / "out" / "workstreams" / "W2" / "results.tsv")
            vi_, si_ = W_[0].index("value"), W_[0].index("se")
            zz_ = None
            for i_ in range(1, len(W_)):
                if O_[i_][vi_] == MASK_TOKEN:
                    try:
                        b_, s_ = float(W_[i_][vi_]), float(W_[i_][si_])
                    except ValueError:
                        continue
                    if s_:
                        zz_ = "%.6g" % abs(b_ / s_)
                        break
            assert zz_, "S2-D41 mutation target: no masked (value, se) pair in W2"
            w2_note_set(d, lambda id_, now_: id_.startswith("S2.marg.ratio.") and "withheld_pending_terms" in now_, lambda w_, r_: "protein z %s" % zz_)

        def s42_rad23b_f_back(d, m, l):     # S2-D42: only the RAD23B row's F put back (the row the S2-D41 exception had left plain)
            def mut(W_, O_):
                j_, g_ = W_[0].index("decode_F_stat"), W_[0].index("gene")
                rr_ = [i_ for i_ in range(1, len(O_)) if W_[i_][g_] == "RAD23B"]
                assert len(rr_) == 1 and O_[rr_[0]][j_] == MASK_TOKEN, "S2-D42 mutation target: the RAD23B row is not masked"
                O_[rr_[0]][j_] = W_[rr_[0]][j_]
            dec_write(d, mut)

        def s42_ledger_leaf_back(d, m, l):  # S2-D42: the RAD23B F printed again in a value cell of the public ledger (the L246 leaf put back)
            W_ = _grid(d / "withheld" / "replication" / "decode_replication__decode_mr_results.csv")
            f_ = [r_[W_[0].index("decode_F_stat")] for r_ in W_[1:] if r_[W_[0].index("gene")] == "RAD23B"][0]
            p_ = d / "out" / "ledger" / "claims_ledger_public.tsv"
            p_.write_text(p_.read_text(encoding="utf-8") + "\t".join(["L999", "s", "planted", "[%s, 0.04, 0.55, true]" % f_, "deCODE F", "src", "not separately checked", ""]) + "\n", encoding="utf-8")
            refresh(d, "ledger/claims_ledger_public.tsv")

        mutate("T2 R11 deCODE F put back (rule F15)", "R11", s41_f_back, real=True)
        mutate("T2 R11 RAD23B deCODE F put back (S2-D42)", "R11", s42_rad23b_f_back, real=True)
        mutate("T2 R14 L246 F leaf printed again (S2-D42)", "R14", s42_ledger_leaf_back, real=True)
        mutate("T2 R16 deCODE |z| column from withheld F", "R16", s41_z_col, real=True)
        mutate("T2 R16 deCODE delta SE/p back, F recoverable", "R16", s41_delta_back, real=True)
        mutate("T2 R18 marginal z put back in a note", "R18", s41_note_quote, real=True)
        mutate("T2 R18 note prints a masked pair's |z|", "R18", s41_note_z, real=True)
        print("tier 2 mutations for the portability checks (R12 executables, R13 every shipped text file):")
        mutate("T2 R12 path in a shipped script", "R12", r12_exec, real=True)
        mutate("T2 R13 personal path in a data cell", "R13", r13_cell, real=True)
        mutate("T2 R13 drive path in REPRODUCE.md", "R13", r13_md, real=True)
        dcm, mcm, _ = fresh("ctl_comment", real=True)
        r12_comment_only(dcm, mcm, None)
        rcm = run_gate(dcm, v3, mcm)
        okcm = "R12" not in rcm.tripped() and "R13" in rcm.tripped()
        print("  %-34s expect R13 only -> tripped %s %s" % ("control: path in an executable comment", sorted(rcm.tripped()), "OK" if okcm else "MISS"))
        tests.append(okcm)
        # positive controls: a DECLARED archive member and the default docx must stay green (a check never seen green proves nothing)
        dctl, mctl, _ = fresh("ctl_archive", real=True)
        add_to_release(dctl, "zz_arch/ok.zip", zip_bytes([("payload.csv", b"a,b\n1,2\n")]))
        with open(dctl / "_gate_fixtures" / "R9_ARCHIVE_MEMBERS.tsv", "a", encoding="utf-8", newline="") as f_:
            f_.write("zz_arch/ok.zip::payload.csv\tdeclared for the test\n")
        rc = run_gate(dctl, v3, mctl)
        okc = not rc.tripped()
        print("  %-34s expect green -> tripped %s %s" % ("control: declared archive member", sorted(rc.tripped()), "OK" if okc else "MISS"))
        tests.append(okc)
        rd = run_gate(bad_real, v3, None)
        okd = not rd.tripped()
        print("  %-34s expect green -> tripped %s %s" % ("control: default docx manuscript", sorted(rd.tripped()), "OK" if okd else "MISS"))
        tests.append(okd)
        if not all(tests):
            ok_all = False
        print("SELFTEST %s: %d/%d mutations tripped the check they name" % ("PASS" if ok_all else "FAIL", sum(tests), len(tests)))
        return 0 if ok_all else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", default=str(HERE))
    ap.add_argument("--v3", default=None)
    ap.add_argument("--manuscript", default=None, help="built .docx, or a folder / file of markdown (default: the operative R2 markdown)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    rel = Path(a.release).resolve()
    v3 = Path(a.v3).resolve() if a.v3 else rel.parent.parent
    if a.selftest:
        return selftest(rel, v3)
    R = run_gate(rel, v3, a.manuscript)
    ok = R.report()
    print("RELEASE GATE: %s" % ("GREEN" if ok else "RED (%s)" % ", ".join(c for c in R.order if c in R.fails)))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
