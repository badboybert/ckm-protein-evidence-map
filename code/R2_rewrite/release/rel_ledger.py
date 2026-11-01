"""rel_ledger.py -- the PUBLIC EXTRACT of the claims ledger (release candidate RC1, S2-D33(3)/(4); brief 'Agent C').

Writes out/ledger/claims_ledger_public.tsv: the ledger the shipped figure generators read their printed numbers from
(P31R_fig_common.R, shipped copy), so that a figure can be regenerated from released results + this file + the code
(reproducibility level 2). It is a DECLARED TRANSFORM of the private ledger; nothing is typed here.

Columns: id, section, claim, value, display, source, evidence_status.
  * dropped: key, tag, load_bearing, note, independent_check (internal notes and verification history);
  * source: made archive-relative (the project prefix and drive roots removed; a path of a file that ships is replaced by its
    path inside the archive); an external dataset keeps its name;
  * evidence_status (EVIDENCE_STATUS below; the README defines it): derived from independent_check and load_bearing;
  * masking: the value and display of a row are replaced by a REASONED mask when the row
      - is TPMI-derived (TPMI is one declared switch, TIER_RULES.tsv row SW1; S2-D34(1): 'ship'). With SW1 'ship' a row is masked only if it holds a RAW
        per-variant TPMI value (beta, SE or allele frequency at the variant), found by content; derived per-pair values and study-level counts ship.
        With SW1 'withheld' every TPMI-named or TPMI-valued row is masked, as before;
      - is ChEMBL-derived (the source names ChEMBL); or
      - cites a source file that is withheld as a whole (RELEASE_INPUTS ship_status withheld*, or a file under withheld/ with no
        counterpart under out/, or a 'withheld' row of build/tiers_by_file.tsv), EXCEPT a row that is Open Targets-native (S2-D33(4)).
    A masked value is the JSON object {"masked": true, "reason": "<code>"} (never 0, never blank), so a generator that needs it
    stops with the reason instead of drawing a wrong number; the display reads "[masked: <reason text>]".
  * partial masking (RC2, S2-D35(1)): a row whose value is a JSON list of word statements keeps its row, and only the ChEMBL-vocabulary development-stage clauses are replaced, element by element,
    by the token "[masked: CHEMBL]". A clause is a ChEMBL development-stage statement when it carries ChEMBL stage vocabulary (APPROVAL, PHASE n, a maximum-phase phrase, or a leading
    'approved:' / 'in clinical development:'); in an element that holds one, the clauses before the first such clause (Open Targets-native, e.g. tractability) stay and that clause and the
    clauses after it (they qualify the named agents) are masked. The rule is general (STAGE_RX); no row id is listed. Provenance (RC2): the md5 a source cell cites is the md5 of the SHIPPED file; the
    pre-transform md5 sits in the column source_md5_before_release_transform.
  * asserted: the extract has the private ledger's ids in the private ledger's order; every unmasked value and display is equal
    to the private one; every masked value differs from it and carries a reason.

Module interface (as rel_screen / rel_tables): INPUTS, build(v3, out, withheld, ledger) -> records.
"""
import csv
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

LEDGER_REL = "R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv"
INPUTS = [LEDGER_REL]
OUT_REL = "ledger/claims_ledger_public.tsv"
PUBLIC_COLS = ["id", "section", "claim", "value", "display", "source", "evidence_status", "source_md5_before_release_transform"]
INPUTS_DIR = "R2_rewrite/release/inputs"

TIER_RULES_REL = INPUTS_DIR + "/TIER_RULES.tsv"


def tpmi_switch(v3):
    """The one declared TPMI switch: TIER_RULES.tsv row SW1 / target TPMI_DERIVED, handling 'withheld' or 'ship' (S2-D34(1): 'ship'). Read, never typed here."""
    rows = [r for r in _read_tsv(Path(v3) / TIER_RULES_REL) if r.get("kind") == "switch" and r.get("target") == "TPMI_DERIVED"]
    _req(len(rows) == 1 and rows[0]["handling"].strip() in ("withheld", "ship"), "TIER_RULES.tsv has no single TPMI_DERIVED switch with handling withheld|ship")
    return rows[0]["handling"].strip()


EVIDENCE_STATUS = {
    "re-derived independently": "a second, separately written computation (or a blinded verifier) reproduced the value from the raw sources",
    "source-reconciled": "the value was read from, or asserted equal to, a named source at build time; one reader, not independently re-derived",
    "not separately checked": "computed by the analysis pipeline and carried as computed; no separate check is recorded",
}

TPMI_RX = re.compile(r"(?<![A-Za-z0-9])TPMI(?![A-Za-z])")
CHEMBL_RX = re.compile(r"chembl", re.I)
OT_RX = re.compile(r"open ?targets|opentargets", re.I)
REASON = {
    "TPMI": "TPMI-derived: the TPMI PheWeb carries no licence, so derived values stay withheld (switch SW1 in TIER_RULES.tsv)",
    "TPMI_RAW": "holds a raw per-variant TPMI value (beta, SE or allele frequency at the variant); only derived per-pair values and study-level counts ship (S2-D34(1))",
    "CHEMBL": "ChEMBL-derived drug or development-stage field: withheld from the data file (reason CHEMBL; see the README, Tables and Licence tiers)",
    "SOURCE_WITHHELD": "the cited source file is withheld from this release",
}

STAGE_RX = re.compile(r"\bAPPROVAL\b|\bPHASE\s*[1-4]\b|\b(?:max(?:imum)?|clinical)[ _]phase\b|^\s*(?:approved|in clinical development)\s*:", re.I)
MASK_TOKEN = "[masked: CHEMBL]"


def split_clauses(s):
    """Split a statement at ';' outside parentheses."""
    out, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        if ch == ";" and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur).strip())
    return out


class _Num(str):
    """A JSON number kept as its source text, so that a rebuilt list prints the unmasked elements exactly as the private ledger does."""


def leaf_rules(v3):
    """S2-D42: TIER_RULES.tsv rows of kind 'ledgerleaf' (target = a ledger row id; handling 'leaf:<withheld file>:<column>:<CODE>'; selector = the display pattern whose group 1 is replaced)."""
    return [r for r in _read_tsv(Path(v3) / TIER_RULES_REL) if r.get("kind") == "ledgerleaf"]


def leaf_keys(withheld, spec):
    """(set of 10-digit keys of the numeric cells of the withheld column, CODE). Fails closed: a missing file or column means the rule is stale."""
    parts = spec.split(":")
    _req(len(parts) == 4 and parts[0] == "leaf", "ledgerleaf handling %r is not leaf:<withheld file>:<column>:<CODE>" % spec)
    _, f, col, code = parts
    p = Path(withheld) / f
    _req(p.is_file(), "ledgerleaf: %s is not in withheld/ (the leaf rule is stale)" % f)
    rows = list(csv.DictReader(open(p, encoding="utf-8", newline="")))
    _req(rows and col in rows[0], "ledgerleaf: %s has no column %s" % (f, col))
    keys = set()
    for r in rows:
        try:
            keys.add("%.10g" % float(r[col]))
        except ValueError:
            pass
    _req(keys, "ledgerleaf: the withheld column %s of %s holds no number" % (col, f))
    return keys, code


def leaf_mask(r, rules, withheld):
    """None, or (masked value text, masked display text) for a private ledger row named by a ledgerleaf rule: every numeric element of its JSON list that equals a withheld value becomes
    the token [masked: CODE]; the number the rule's display pattern captures becomes the token too. Asserts at least one element and exactly one display match."""
    for rule in rules:
        if rule["target"] != r["id"]:
            continue
        keys, code = leaf_keys(withheld, rule["handling"])
        v = r["value"].strip()
        _req(v.startswith("["), "%s: a ledgerleaf row must have a JSON list value" % r["id"])
        els = json.loads(v, parse_float=_Num, parse_int=_Num)
        tok = "[masked: %s]" % code
        new, n = [], 0
        for e in els:
            if isinstance(e, _Num) and ("%.10g" % float(e)) in keys:
                new.append(tok)
                n += 1
            else:
                new.append(e)
        _req(n >= 1, "%s: no list element equals a withheld value of %s (the leaf rule is stale)" % (r["id"], rule["handling"]))
        rx = re.compile(rule["selector"])
        ms = list(rx.finditer(r["display"]))
        _req(len(ms) == 1, "%s: the display pattern %r matches %d time(s), want 1" % (r["id"], rule["selector"], len(ms)))
        m = ms[0]
        dis = r["display"][:m.start(1)] + tok + r["display"][m.end(1):]
        val = "[" + ", ".join(e if isinstance(e, _Num) else json.dumps(e, ensure_ascii=False) for e in new) + "]"
        return val, dis
    return None


def mask_stage_element(el):
    """(masked element, changed). Clauses before the first ChEMBL stage clause stay; that clause and the ones after it become one MASK_TOKEN."""
    cl = split_clauses(el)
    k = next((i for i, c in enumerate(cl) if STAGE_RX.search(c)), None)
    if k is None:
        return el, False
    return "; ".join(cl[:k] + [MASK_TOKEN]), True


def partial_mask(value, display):
    """None, or (masked value text, masked display text) for a row whose value is a JSON list of word statements that carries ChEMBL stage clauses. Asserts display == '; '.join(elements)."""
    v = value.strip()
    if not v.startswith("["):
        return None
    try:
        els = json.loads(v)
    except ValueError:
        return None
    if not isinstance(els, list) or not all(isinstance(e, str) for e in els):
        return None
    res = [mask_stage_element(e) for e in els]
    if not any(c for _, c in res):
        return None
    _req(display == "; ".join(els), "partial mask: the display is not the '; '-join of the value elements")
    new = [m for m, _ in res]
    return json.dumps(new, ensure_ascii=False), "; ".join(new)


# The home path is written as concatenated fragments on purpose: this module is itself shipped, and the release gate's absolute-path scans (R12, R13) read every shipped file,
# including this one, so the tool's own constant must not match the pattern it enforces.
_CKM = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")).replace("\\", "/")   # the project root above V3, from this file's place
_DRIVE_RX = re.compile(r"(?<![A-Za-z0-9_])[CDE]:/")
_V3_RX = re.compile(r"(?<![A-Za-z0-9_.-])V3/")


class LedgerError(AssertionError):
    pass


def _req(cond, msg):
    if not cond:
        raise LedgerError(msg)


def _read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def evidence_status(independent_check, load_bearing):
    """re-derived independently / source-reconciled / not separately checked. Deterministic; the order of the tests is the rule."""
    ic = (independent_check or "").strip()
    if not ic:
        return "not separately checked"
    if ic.startswith(("PIPE", "KEY-STALE")):                  # KEY-STALE: the verifier certified a superseded value
        return "not separately checked"
    if ic.startswith("KEY-VISIBLE") or "derivable from visible" in ic:
        return "source-reconciled"                            # the answer was visible to the verifiers: not a blind re-derivation
    if re.search(r"\bnot independently re-derived\b", ic):
        return "source-reconciled"
    if (ic.startswith(("MATCH", "KEY ")) or re.search(r"independent re-derivation|re-derived|reproduced|: MATCH", ic)) and not ic.startswith("MISMATCH"):
        return "re-derived independently"
    if ic.startswith(("SRC-QUOTE", "DERIVED", "MISMATCH", "raised by")):
        return "source-reconciled"
    return "source-reconciled"


def _norm(p):
    return str(p).replace("\\", "/")


def _shipped_maps(v3):
    """source path (project-relative) -> path inside the archive, from RELEASE_INPUTS (copies that ship) and CODE_INCLUDE (ship=1)."""
    m = {}
    d = Path(v3) / INPUTS_DIR
    for r in _read_tsv(d / "RELEASE_INPUTS.tsv"):
        if r["ship_status"] == "ship" and r["release_path"] and r["source_path"]:
            m[_norm(r["source_path"])] = _norm(r["release_path"])
    for r in _read_tsv(d / "CODE_INCLUDE.tsv"):
        if r["ship"] == "1":
            m[_norm(r["path"])] = _norm(r["ship_path"])
    m[LEDGER_REL] = OUT_REL
    return m


def _withheld_sources(v3, out, withheld):
    """Project-relative source paths of files withheld AS A WHOLE (no counterpart under out/). Derived, not listed."""
    d = Path(v3) / INPUTS_DIR
    ri = _read_tsv(d / "RELEASE_INPUTS.tsv")
    by_release = {_norm(r["release_path"]): _norm(r["source_path"]) for r in ri if r["release_path"] and r["source_path"]}
    res = {}
    for r in ri:
        if r["ship_status"].startswith("withheld") and r["source_path"]:
            res[_norm(r["source_path"])] = "RELEASE_INPUTS ship_status %s; %s" % (r["ship_status"], r["note"])
    wh, o = Path(withheld), Path(out)
    if wh.is_dir():
        for p in wh.rglob("*"):
            if p.is_file():
                rel = p.relative_to(wh).as_posix()
                if not (o / rel).exists() and rel in by_release:
                    res[by_release[rel]] = "withheld/ has the file and out/ has no counterpart"
    tf = Path(out).parent / "build" / "tiers_by_file.tsv"
    if tf.is_file():
        for r in _read_tsv(tf):
            rel = _norm(r.get("file", ""))
            if r.get("decision", "").startswith("withheld") and rel in by_release and not (o / rel).exists():
                res[by_release[rel]] = "tiers_by_file.tsv: " + r.get("reason", "")[:80]
    return res


def _sources_in(text, known):
    """Which of the known project-relative paths does a ledger source cell name (longest match, bounded)?"""
    found = []
    for p in sorted(known, key=len, reverse=True):
        if re.search(r"(?<![A-Za-z0-9_./-])" + re.escape(p) + r"(?![A-Za-z0-9_/-]|\.[A-Za-z0-9])", text):
            found.append(p)
    return found


def public_source(src, shipped):
    s = src.replace(_CKM + "/V3/", "").replace(_CKM + "/", "")
    s = _V3_RX.sub("", s)
    # a path of a file that ships becomes its path inside the archive (longest first, bounded)
    for p in sorted(shipped, key=len, reverse=True):
        s = re.sub(r"(?<![A-Za-z0-9_./-])" + re.escape(p) + r"(?![A-Za-z0-9_/-]|\.[A-Za-z0-9])", lambda m, q=shipped[p]: q, s)
    s = _DRIVE_RX.sub("<external>/", s)
    return re.sub(r"[\t\r\n]+", " ", s).strip()


def _cell(v):
    return re.sub(r"[\t\r\n]+", " ", "" if v is None else str(v))


NUM_TOK_RX = re.compile(r"(?<![\w.])\d[\d.]*(?:[eE][+-]?\d+)?(?![\w.]*\d)")


def _key(tok):
    """A number as 12 significant digits: the ledger may print a value shortened to fewer digits than the source file holds (the same value)."""
    try:
        v = float(tok.lstrip("-+"))
    except ValueError:
        return None
    return format(v, ".12g") if v != 0 else None


RAW_COL_RX = re.compile(r"(?i)^(b|beta)_out$|^se_out$|^eaf(_|$)")                    # columns of a check table that hold the resource's own per-variant effect, SE or frequency
RAW_QID_RX = re.compile(r"(?i)\.(beta|se|eaf)_(exp|out)$")                            # rows of a results table that hold them (quantity_id)
RAW_MIN_LEN = 6                                                                      # a raw token is a distinctive string, not a rounded number


def tpmi_raw_values(withheld, out):
    """12-significant-digit keys of the RAW per-variant TPMI values (effect, SE, allele frequency at the variant), read by content from every table that has a
    `resource` column: the columns named b_out / se_out / eaf* of a TPMI row, and the value / se cells of a TPMI row whose quantity_id is a .beta_out / .se_out / .eaf_out
    quantity. Read from the build's own withheld/ copies and the shipped tables, so nothing is typed. The derived values (Wald ratio, SE, p, class, counts) are not in the set."""
    vals = set()
    for root in (Path(withheld), Path(out)):
        if not root.is_dir():
            continue
        for p in sorted(root.rglob("*")):
            if not p.is_file() or p.suffix not in (".tsv", ".csv"):
                continue
            try:
                rows = _read_tsv(p) if p.suffix == ".tsv" else list(csv.DictReader(open(p, encoding="utf-8", newline="")))
            except (UnicodeDecodeError, csv.Error):
                continue
            for r in rows:
                if "TPMI" not in (r.get("resource") or ""):
                    continue
                cells = [v for c, v in r.items() if RAW_COL_RX.match(c or "")]
                if RAW_QID_RX.search(r.get("quantity_id") or ""):
                    cells += [r.get("value"), r.get("se")]
                for v in cells:
                    for t in NUM_TOK_RX.findall(v or ""):
                        if len(t) >= RAW_MIN_LEN and _key(t):
                            vals.add(_key(t))
    return vals


def tpmi_values(withheld, out):
    """Distinctive numeric strings (>= 8 characters, sign dropped) of every TPMI row of a withheld full copy or a shipped table that has a `resource` column. A ledger row
    that prints one of them is TPMI-derived even when its text does not name TPMI (L223 prints the six external checks as one list). Read from the build's own
    withheld/ copies (byte-identical to the declared sources), so nothing is typed."""
    vals = set()
    for root in (Path(withheld), Path(out)):
        if not root.is_dir():
            continue
        for p in sorted(root.rglob("*")):
            if not p.is_file() or p.suffix not in (".tsv", ".csv"):
                continue
            try:
                rows = _read_tsv(p) if p.suffix == ".tsv" else list(csv.DictReader(open(p, encoding="utf-8", newline="")))
            except (UnicodeDecodeError, csv.Error):
                continue
            for r in rows:
                if "TPMI" in (r.get("resource") or ""):
                    for v in r.values():
                        for t in NUM_TOK_RX.findall(v or ""):
                            if len(t) >= 8 and _key(t):
                                vals.add(_key(t))
    return vals


_AGENT_PAREN_RX = re.compile(r" \([A-Za-z][A-Za-z0-9 ,'/-]*\)$")


def _claim_public(claim, mk):
    """G5 (2026-11-01, step 5f list): the claim text of a CHEMBL-masked row may end in a parenthetical naming the agent ('development-stage note (<agent>)',
    L110-L113); the mask covers that parenthetical too, so the public claim column names no ChEMBL-derived agent."""
    if mk and mk[0] == "CHEMBL" and _AGENT_PAREN_RX.search(claim):
        return _AGENT_PAREN_RX.sub(" (agent name masked: CHEMBL)", claim)
    return claim


def mask_reason(r, withheld_src, tpmi_toks=frozenset(), tpmi_mode="withheld", raw_toks=frozenset()):
    """None, or (code, detail) for a private ledger row. tpmi_mode is the SW1 switch ('withheld' | 'ship'); raw_toks are the raw per-variant TPMI keys (ship mode)."""
    txt = " ".join(r[c] for c in ("claim", "source", "display", "value"))
    row_toks = {_key(t) for t in NUM_TOK_RX.findall(r["value"] + " " + r["display"]) if len(t) >= 6}
    if tpmi_mode == "ship" and raw_toks and raw_toks & row_toks:
        return "TPMI_RAW", REASON["TPMI_RAW"]
    if tpmi_mode == "withheld" and tpmi_toks and tpmi_toks & {_key(t) for t in NUM_TOK_RX.findall(r["value"] + " " + r["display"]) if len(t) >= 8}:
        return "TPMI", REASON["TPMI"]
    ot_native = bool(OT_RX.search(txt)) and not CHEMBL_RX.search(r["source"]) and not (tpmi_mode == "withheld" and TPMI_RX.search(txt))
    if tpmi_mode == "withheld" and TPMI_RX.search(txt):
        return "TPMI", REASON["TPMI"]
    if CHEMBL_RX.search(r["source"]) or CHEMBL_RX.search(r["value"]) or CHEMBL_RX.search(r["display"]):
        return "CHEMBL", REASON["CHEMBL"]
    if not ot_native:
        named = _sources_in(_V3_RX.sub("", r["source"].replace(_CKM + "/V3/", "").replace(_CKM + "/", "")), withheld_src)
        # a file withheld ONLY because it prints a TPMI line does not mask its other rows: the TPMI rule masks the TPMI-derived ones
        if named and not TPMI_RX.search(withheld_src[named[0]]):
            return "SOURCE_WITHHELD", REASON["SOURCE_WITHHELD"] + ": " + named[0]
    return None


MD5_RX = re.compile(r"(?P<path>[^\s;+]+)(?P<mid>(?:; (?:lines?|columns?|sheet|rows?|key)\b[^;]*)*); md5 (?P<md5>[0-9a-f]{32})")


def _file_md5(path):
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def shipped_md5_source(src_public, out):
    """(source with every cited md5 replaced by the md5 of the SHIPPED file when the cited path is a file of out/, the list of the md5s as cited). A path that is not a file of out/ keeps its md5."""
    before = []

    def sub(m):
        before.append(m.group("md5"))
        f = Path(out) / m.group("path")
        new = _file_md5(f) if f.is_file() else m.group("md5")
        return m.group("path") + m.group("mid") + "; md5 " + new
    return MD5_RX.sub(sub, src_public), before


def make_public(v3, out, withheld, ledger_rows):
    shipped = _shipped_maps(v3)
    wsrc = _withheld_sources(v3, out, withheld)
    mode = tpmi_switch(v3)
    ttoks = tpmi_values(withheld, out) if mode == "withheld" else frozenset()
    rtoks = tpmi_raw_values(withheld, out) if mode == "ship" else frozenset()
    pub, stats = [], {"TPMI": 0, "TPMI_RAW": 0, "CHEMBL": 0, "SOURCE_WITHHELD": 0, "CHEMBL_CLAUSES": 0, "LEAF": 0}
    leafr = leaf_rules(v3)
    for r in ledger_rows:
        mk = mask_reason(r, wsrc, ttoks, mode, rtoks)
        val, dis = r["value"], r["display"]
        part = None if mk else partial_mask(val, dis)
        if part:
            val, dis = part
            stats["CHEMBL_CLAUSES"] += 1
        leaf = None if (mk or part) else leaf_mask(r, leafr, withheld)
        if leaf:
            val, dis = leaf
            stats["LEAF"] += 1
        if mk:
            code, why = mk
            val = json.dumps({"masked": True, "reason": code}, separators=(",", ":"))
            dis = "[masked: %s]" % why
            stats[code] += 1
        src, before = shipped_md5_source(public_source(r["source"], shipped), out)
        pub.append(dict(id=r["id"], section=_cell(_DRIVE_RX.sub("<external>/", r["section"])), claim=_cell(_DRIVE_RX.sub("<external>/", _claim_public(r["claim"], mk))),
                        value=_cell(val), display=_cell(dis),
                        source=src, evidence_status=evidence_status(r["independent_check"], r["load_bearing"]),
                        source_md5_before_release_transform=";".join(before),
                        _masked=bool(mk), _partial=bool(part), _leaf=bool(leaf)))
    return pub, stats


def write_public(pub, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None, escapechar=None)
        w.writerow(PUBLIC_COLS)
        for r in pub:
            w.writerow([r[c] for c in PUBLIC_COLS])


def check_public(path, ledger_rows, pub):
    """Row identity with the private ledger: same ids in the same order; unmasked values and displays equal; masked ones differ and carry a reason."""
    got = _read_tsv(path)
    _req([g["id"] for g in got] == [r["id"] for r in ledger_rows], "public ledger ids differ from the private ledger's")
    _req(list(got[0].keys()) == PUBLIC_COLS, "public ledger columns are %s" % list(got[0].keys()))
    nm = 0
    for g, r, p in zip(got, ledger_rows, pub):
        if p["_masked"]:
            nm += 1
            j = json.loads(g["value"])
            _req(isinstance(j, dict) and j.get("masked") is True and j.get("reason"), "%s: masked value without a reason: %s" % (g["id"], g["value"]))
            _req(g["value"] != r["value"] and g["display"].startswith("[masked: ") and len(g["display"]) > 12, "%s: mask does not replace the value or has no reason" % g["id"])
        elif p.get("_partial"):
            nm_ = 0
            gv, rv = json.loads(g["value"]), json.loads(r["value"])
            _req(len(gv) == len(rv) and g["display"] == "; ".join(gv), "%s: partially masked row: elements or display differ in shape" % g["id"])
            for ge, re_ in zip(gv, rv):
                if ge != re_:
                    nm_ += 1
                    _req(ge.endswith(MASK_TOKEN) and re_.startswith(ge[:-len(MASK_TOKEN)]) and any(STAGE_RX.search(c) for c in split_clauses(re_)),
                         "%s: a changed element is not a ChEMBL stage clause replaced by the mask token" % g["id"])
                    _req(not any(STAGE_RX.search(c) for c in split_clauses(ge)), "%s: a masked element still carries ChEMBL stage vocabulary" % g["id"])
            _req(nm_ > 0, "%s: partially masked row with no masked element" % g["id"])
        elif p.get("_leaf"):
            gv, rv = json.loads(g["value"]), json.loads(r["value"])
            _req(len(gv) == len(rv), "%s: leaf-masked row: the list changed length" % g["id"])
            ch = [(a_, b_) for a_, b_ in zip(gv, rv) if a_ != b_]
            _req(ch and all(isinstance(a_, str) and a_.startswith("[masked: ") and isinstance(b_, (int, float)) and not isinstance(b_, bool) for a_, b_ in ch),
                 "%s: a changed element is not a number replaced by a [masked: ...] token" % g["id"])
            _req(g["display"] != r["display"] and "[masked: " in g["display"], "%s: leaf-masked row whose display is not masked" % g["id"])
        else:
            _req(g["value"] == _cell(r["value"]) and g["display"] == _cell(r["display"]), "%s: unmasked value/display differs from the private ledger" % g["id"])
            # the general rule (RC2): no unmasked list element of the extract carries ChEMBL development-stage vocabulary
            if g["value"].startswith("["):
                try:
                    els = json.loads(g["value"])
                except ValueError:
                    els = []
                _req(not any(isinstance(e, str) and any(STAGE_RX.search(c) for c in split_clauses(e)) for e in els),
                     "%s: an unmasked list element carries ChEMBL development-stage vocabulary" % g["id"])
        _req(g["evidence_status"] in EVIDENCE_STATUS, "%s: evidence_status %r" % (g["id"], g["evidence_status"]))
        _req(not re.search(r"(?<![A-Za-z0-9_])[CDE]:/|" + re.escape(os.path.expanduser("~").replace("\\", "/")) + r"|\b" + re.escape(os.path.basename(os.path.expanduser("~"))) + r"\b"," ".join(g.values())), "%s: an absolute local path survives in the public row" % g["id"])
        _req(g["value"] not in ("", "0"), "%s: blank or zero value" % g["id"]) if p["_masked"] else None
    return len(got), nm


def _licence_class(v3):
    for r in _read_tsv(Path(v3) / INPUTS_DIR / "RELEASE_INPUTS.tsv"):
        if _norm(r["source_path"]) == LEDGER_REL and r["licence_class"].startswith("third_party_derived:"):
            toks = r["licence_class"].split(":", 1)[1].split("+")
            if tpmi_switch(v3) == "withheld":
                toks = [t for t in toks if t != "TPMI"]
            return "third_party_derived:" + "+".join(toks) if toks else "own_derived"
    return "own_derived"


def _rows_of(ledger):
    cols = ["id", "section", "claim", "value", "display", "source", "key", "tag", "load_bearing", "note", "independent_check"]
    if hasattr(ledger, "to_dict"):
        rows = ledger.to_dict("records")
    else:
        rows = [dict(r) for r in ledger]
    for c in cols:
        _req(all(c in r for r in rows), "ledger has no column %s" % c)
    rows = [{k: ("" if v is None else str(v)) for k, v in r.items()} for r in rows]
    _req(len({r["id"] for r in rows}) == len(rows) and len(rows) >= 100, "ledger ids not unique / implausible row count %d" % len(rows))
    return rows


def build(v3, out, withheld, ledger):
    """Spec interface: write out/ledger/claims_ledger_public.tsv; return one record."""
    rows = _rows_of(ledger)
    pub, stats = make_public(v3, out, withheld, rows)
    dst = Path(out) / OUT_REL
    write_public(pub, dst)
    n, nm = check_public(dst, rows, pub)
    ev = {k: sum(1 for p in pub if p["evidence_status"] == k) for k in EVIDENCE_STATUS}
    note = ("public extract of the claims ledger (%d rows, %d masked: TPMI %d, raw TPMI %d, ChEMBL %d, withheld source %d; %d more rows with ChEMBL stage clauses masked element by element; %d row(s) with a withheld-value leaf masked element by element; TPMI switch '%s'); evidence_status: %s"
            % (n, nm, stats["TPMI"], stats["TPMI_RAW"], stats["CHEMBL"], stats["SOURCE_WITHHELD"], stats["CHEMBL_CLAUSES"], stats["LEAF"], tpmi_switch(v3), ", ".join("%s %d" % (k, v) for k, v in ev.items())))
    return [dict(release_path=OUT_REL, ship_status="ship", source_paths=[LEDGER_REL], producer="rel_ledger.py", licence_class=_licence_class(v3),
                 rows=n, columns=len(PUBLIC_COLS), note=note)]


# ---------------------------------------------------------------------------------------------------------------- self-test
def selftest(v3):
    """Each mutation must trip the check it names."""
    import pandas as pd
    v3 = Path(v3)
    led = pd.read_csv(v3 / LEDGER_REL, dtype=str, keep_default_na=False)
    rows = _rows_of(led)
    tmp = Path(tempfile.mkdtemp(prefix="rel_ledger_selftest_"))
    bad = []

    def rep(name, ok, detail=""):
        print("  %-66s %s %s" % (name, "ok" if ok else "MISS", detail))
        if not ok:
            bad.append(name)

    here = Path(__file__).resolve().parent
    out, wh = tmp / "out", here / "withheld"          # the real withheld/ copies carry the raw TPMI row the content test reads; a scratch copy of out/ is what they are compared with
    if (here / "out").is_dir():
        shutil.copytree(here / "out", out)
    rec = build(v3, out, wh, led)
    rep("build writes the extract", (out / OUT_REL).is_file(), rec[0]["note"][:90])
    pub, _ = make_public(v3, out, wh, rows)
    mode = tpmi_switch(v3)
    byid = {p["id"]: p for p in pub}
    rtoks = tpmi_raw_values(wh, out)
    rep("raw TPMI tokens are found by content (beta, SE, frequency of the TPMI row)", len(rtoks) >= 3, "%d keys" % len(rtoks))
    if mode == "ship":
        derived = ["L045", "L217", "L223", "L261", "L375", "L905", "L960"]
        rep("S2-D34(1): the derived per-pair and count rows ship unmasked", all(not byid[i]["_masked"] for i in derived), "rows %s" % derived)
        base = next(r for r in rows if r["id"] == "L045")
        raw_tok = sorted(rtoks, key=len)[-1]
        planted = dict(base, display=base["display"] + " (beta " + raw_tok + ")")
        rep("M5 a row that prints a raw TPMI value is masked (TPMI_RAW)", (mask_reason(planted, {}, frozenset(), "ship", rtoks) or ("", ""))[0] == "TPMI_RAW")
        rep("M6 the same row with the switch 'withheld' is masked as TPMI", (mask_reason(base, {}, frozenset(), "withheld") or ("", ""))[0] == "TPMI")
        rep("no ledger row holds a raw per-variant TPMI value", not any(p["value"].startswith('{"masked":true,"reason":"TPMI_RAW"') for p in pub))
    else:
        tp = [p for p in pub if p["_masked"] and "TPMI" in p["value"]]
        rep("a TPMI-derived row is masked with a reason", len(tp) > 0 and all("TPMI" in p["display"] for p in tp), "%d rows" % len(tp))
    rep("no masked value is blank or zero", all(p["value"] not in ("", "0") for p in pub if p["_masked"]))

    def trips(name, fn, token):
        try:
            fn()
            rep(name, False, "no check fired")
        except LedgerError as e:
            rep(name, token in str(e), str(e)[:70])

    p2 = [dict(p) for p in pub]
    i = next(k for k, p in enumerate(p2) if not p["_masked"])
    p2[i]["value"] = "999"
    write_public(p2, tmp / "m1.tsv")
    trips("M1 an altered unmasked value", lambda: check_public(tmp / "m1.tsv", rows, p2), "differs from the private ledger")
    p3 = [dict(p) for p in pub]
    j = next(k for k, p in enumerate(p3) if p["_masked"])
    p3[j]["value"] = rows[j]["value"]
    write_public(p3, tmp / "m2.tsv")
    trips("M2 an unmasked value in a row that must be masked", lambda: check_public(tmp / "m2.tsv", rows, p3), "")
    p4 = [dict(p) for p in pub][::-1]
    write_public(p4, tmp / "m3.tsv")
    trips("M3 reordered ids", lambda: check_public(tmp / "m3.tsv", rows, p4), "ids differ")
    p5 = [dict(p) for p in pub]
    p5[i]["source"] = "C" + ":/Us" + "ers/example/x.tsv"
    write_public(p5, tmp / "m4.tsv")
    trips("M4 an absolute local path in a row", lambda: check_public(tmp / "m4.tsv", rows, p5), "absolute local path")
    # RC2 (S2-D35(1)): the ChEMBL stage-clause rule is general and asserted
    for lid in ("L370", "L371", "L374", "L375"):
        rep("P1 %s: ChEMBL stage clauses masked element by element, row kept" % lid, byid[lid]["_partial"] and not byid[lid]["_masked"] and MASK_TOKEN in byid[lid]["value"]
            and not STAGE_RX.search(byid[lid]["value"]), byid[lid]["display"][:60])
    rep("P1 L374 keeps the Open Targets-native 'tractable' part", "tractable, no ARG1 mechanism; " + MASK_TOKEN in byid["L374"]["value"] and "pegzilarginase" not in byid["L374"]["display"])
    gen = partial_mask('["a", "approved: SOMEDRUG inhibits X; agrees", "tractable, no clinical candidate", "x (PHASE 2) y"]', "a; approved: SOMEDRUG inhibits X; agrees; tractable, no clinical candidate; x (PHASE 2) y")
    rep("P2 the rule is general (a new drug clause is masked, a tractability element is not)", gen is not None and "SOMEDRUG" not in gen[0] and "tractable, no clinical candidate" in gen[0] and "PHASE" not in gen[0])
    p6 = [dict(p) for p in pub]
    k6 = next(k for k, p in enumerate(p6) if p["id"] == "L374")
    p6[k6]["value"] = rows[k6]["value"]
    p6[k6]["display"] = rows[k6]["display"]
    p6[k6]["_partial"] = False
    write_public(p6, tmp / "m6.tsv")
    trips("P3 a stage clause left unmasked in a list row", lambda: check_public(tmp / "m6.tsv", rows, p6), "ChEMBL development-stage vocabulary")
    p7 = [dict(p) for p in pub]
    k7 = next(k for k, p in enumerate(p7) if p["id"] == "L375")
    p7[k7]["value"] = p7[k7]["value"].replace(MASK_TOKEN, "in clinical development: FGF401", 1)
    p7[k7]["display"] = "; ".join(json.loads(p7[k7]["value"]))
    write_public(p7, tmp / "m7.tsv")
    trips("P4 a partially masked row whose masked clause is put back as plain stage text", lambda: check_public(tmp / "m7.tsv", rows, p7), "not a ChEMBL stage clause replaced")
    d = tmp / "md5out"
    (d / "w").mkdir(parents=True)
    (d / "w" / "r.tsv").write_bytes(b"shipped\n")
    s_, b_ = shipped_md5_source("w/r.tsv; line 3; column value; md5 " + "0" * 32 + " + other/x.tsv; line 1; md5 " + "1" * 32, d)
    rep("Q1 the cited md5 becomes the shipped file's; a non-shipped path keeps its md5; the cited ones are kept", _file_md5(d / "w" / "r.tsv") in s_ and "1" * 32 in s_ and "0" * 32 not in s_ and b_ == ["0" * 32, "1" * 32])
    rows_pub = _read_tsv(out / OUT_REL)
    chk_ = [r for r in rows_pub if r["source_md5_before_release_transform"]]
    rep("Q2 the new column holds the pre-transform md5 of every row that cites one", len(chk_) > 100 and all(re.fullmatch(r"[0-9a-f]{32}(;[0-9a-f]{32})*", r["source_md5_before_release_transform"]) for r in chk_), "%d rows" % len(chk_))
    print("rel_ledger selftest: %s" % ("PASS" if not bad else "FAIL: %s" % bad))
    return not bad


if __name__ == "__main__":
    import sys
    v3_ = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[2])
    sys.exit(0 if selftest(v3_) else 1)
