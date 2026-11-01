#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 R2, Stage 5 step 2: the release builder (core).

One entry point. It reads ONLY the inputs declared in inputs/RELEASE_INPUTS.tsv (each pinned by md5), writes what ships to out/, what is held
back to withheld/, and its records to build/.  Nothing is typed here that a file can supply: counts come from the files, the held-back column
lists come from figures/af3/P31R_ACTIVE_FIGURES.json (hold H9), the frame lists come from the registries and from the writer calls in
P31R_Fig2.R.  It never renders a figure, never edits a ledger / gate / manifest / manuscript file, never downloads.

Modes (exactly one):
  (default)            strict build: refuse on any md5 mismatch, any unclaimed or leaking file, any frame that the render has not yet written.
  --stage2             as default, but frames that only the step 1c render will create are listed in build/pending_render.tsv instead of refusing.
  --pin                record the current md5 / bytes / cited_by of every declared input in RELEASE_INPUTS.tsv, re-derive the frame rows from the
                       registries, print every changed pin.  Never runs implicitly.
  --derive-not-shipped write inputs/NOT_SHIPPED.tsv by script (overrides are hand-declared in NOT_SHIPPED_OVERRIDES.tsv and never rewritten).
  --providers          write inputs/PROVIDERS.tsv (provider list; a normal build regenerates it in memory and asserts equality, P4).
  --check              rebuild into a temp folder and assert byte identity with out/, withheld/ and build/.
  --selftest           small mutation suite: each defect must trip the check it NAMES (id match on letter+number).
Options: --v3 DIR (default: the V3 that contains this release folder), --release DIR (default: this folder; scratch copies pass their own).

RELEASE_INPUTS.tsv columns that --pin rewrites: source_md5, source_bytes, cited_by, licence_class (derived from PROVIDER_DECLARATIONS.tsv; asserted by L3) and the derived frame rows. Anything that reads this file for a
count (the step 1c ledger rows) must read only the other columns, or the ledger and the pins never settle on a fixed point (the ledger CSV is itself pinned).
licence_class "module_reported" (module_output rows only) means: the class is whatever the module returns for that file (it knows its providers).
Hand-declared tables under inputs/ (DECLARED; no mode rewrites them): PROVIDER_ALIASES, PROVIDER_DECLARATIONS, PROVIDERS_PENDING_DECLARATIONS, CANONICAL_SOURCES,
ORPHAN_FRAMES, NOT_SHIPPED_DECLARED, NOT_SHIPPED_OVERRIDES, DERIVED_PAIR_KEYS, DERIVED_ROWS. SOURCE_CLASS_RULES is hand-declared too (class + basis).
Check ids (letter + number, matched exactly): P pins (P4 PROVIDERS.tsv stale), D declarations (D3 provider declaration, D4 typed count, D5 note path, D6 cite token,
D7 provider list, D8 provider coverage of Table S2 / SM), L licence (L3 derived class, L4 module provider), M modules, C claims, S leak scan (S8-S11 pair-keyed derived scan), H hold H9 (H3 row content),
N not-shipped (N4 undeclared override), F frames (F2 undeclared orphan), G canonical copies, R regenerability, Q copyright, K holds in text, E scope, X cited_by drift.
Judgement calls and the checks they lean on are listed in build/BUILD_LOG.txt after every build.  Python: python -I -X utf8.
"""
import sys
sys.dont_write_bytecode = True
import argparse, collections, csv, fnmatch, hashlib, importlib, io, json, math, os, re, shutil, tempfile
from pathlib import Path

csv.field_size_limit(2 ** 30)

HERE = Path(__file__).resolve().parent
DEFAULT_V3 = HERE.parent.parent

SENT_PENDING = "withheld_pending_terms"
SENT_H9 = "withheld_hold_H9"
SENTINELS = {"tpmi": SENT_PENDING, "ot": SENT_PENDING, "h9": SENT_H9}
ROLES = {"copy", "input", "module_output"}
ORIGINS = {"curated", "derived_frames"}
SHIPS = {"ship", "withheld_pending_terms"}
SHIP_H9 = "withheld_source_hold_H9"     # built rows only: the full counterpart of a file masked by hold H9 alone (source hold, not provider terms)
LICENCE_RX = re.compile(r"own_derived|not_redistributable|third_party_derived:[A-Za-z0-9_+.-]+")    # no space in any token (F18)


PRE_DOC_MODULES = ("rel_portability", "rel_ledger", "rel_records")     # RC1: run after the code bundle is copied and the tier rules are applied, before the regenerability record
LATE_MODULES = PRE_DOC_MODULES + ("rel_docs",)     # rel_docs runs after every other record exists (it documents the release); everything else in MODULES runs first
DOC_FILES = ("RELEASE_MANIFEST.tsv", "README.md", "REPRODUCE.md", "LICENSE", "CITATION.cff", ".gitattributes")     # written by rel_docs; not leak-scanned by the core (the release gate's R5 scans them and code/)


LICENCE_GRAMMAR_RX = re.compile(r"own_derived|third_party_derived:[A-Za-z0-9_.-]+(?:\+[A-Za-z0-9_.-]+)*")


def licence_grammar_ok(s):
    """V13: the ONE grammar of a shipped row's class: 'own_derived' or 'third_party_derived:A+B+...' (no '|', no space)."""
    return bool(LICENCE_GRAMMAR_RX.fullmatch(s or ""))


def licence_ok(s):
    """One class, or several joined by '|' (the table modules report a file built from more than one provider)."""
    return bool(s) and all(LICENCE_RX.fullmatch(x) for x in s.split("|"))
MODULES = {"rel_screen": "screen/", "rel_tables": "tables/", "rel_portability": "code/", "rel_ledger": "ledger/", "rel_records": "release_records/", "rel_docs": ""}   # module -> the only out/ prefix it may write under
RI_COLS = ["input_id", "origin", "role", "deliverable", "release_path", "source_path", "source_md5", "source_bytes", "licence_class",
           "ship_status", "withhold_rule", "producer", "cite_keys", "cited_by", "note"]
BUILT_COLS = ["release_path", "bytes", "md5", "rows", "columns", "source_path", "source_md5", "producer", "licence_class", "ship_status",
              "cited_by", "note"]
NS_COLS = ["path", "md5", "bytes", "reason_code", "evidence", "md5_basis"]
NSO_COLS = ["path", "md5", "manifest_line", "reason"]
RULE_COLS = ["rule_id", "pattern", "class", "resource", "reason", "list_not_shipped", "basis"]
# The plan's regenerability vocabulary (STAGE5_PLAN l.121-127): three classes for the downstream units, plus the stored engine outputs (item 2).
# `basis` (SOURCE_CLASS_RULES.tsv) says WHY an input is not in the release; it is a column, never a fourth class.
CLASS_ORDER = ["regenerates_from_release", "needs_public_download", "needs_nonredistributable_input", "does_not_regenerate_engine_output"]
RULE_CLASSES = set(CLASS_ORDER[1:])
RULE_BASES = {"not_released_pending_decision", "public", "third_party_terms", "individual_level", "terms_unread_or_hold", "not_on_approved_list", "engine_output", "code_bundle_step3"}
# hand-declared tables under inputs/ (never rewritten by --pin or --derive-not-shipped): file -> columns
DECLARED = {
    "PROVIDER_ALIASES.tsv": ["provider_id", "display_name", "aliases", "role", "release_handling"],
    "PROVIDER_DECLARATIONS.tsv": ["decl_id", "pattern", "providers", "inherits", "evidence"],
    "PROVIDERS_PENDING_DECLARATIONS.tsv": ["provider_id", "reason"],
    "CANONICAL_SOURCES.tsv": ["canon_id", "file_name", "source_path", "md5", "evidence"],
    "ORPHAN_FRAMES.tsv": ["path", "md5", "reason"],
    "NOT_SHIPPED_DECLARED.tsv": ["path", "reason_code", "evidence"],
    "NOT_SHIPPED_OVERRIDES.tsv": NSO_COLS,
    "DERIVED_PAIR_KEYS.tsv": ["key_id", "class", "target", "key_regex", "reason"],
    "LEDGER_DYNAMIC_DECLARED.tsv": ["unit", "line_text", "reason"],
    "DERIVED_ROWS.tsv": ["row_id", "class", "release_glob", "where", "status", "derived_cols", "label_cols", "quantities", "reason"],
}
PROVIDER_COLS = ["provider_id", "display_name", "aliases", "role", "named_in", "in_declarations", "licence_class", "terms_status", "release_handling"]
OPERATIVE_DIRS = ("sections_R2", "figures/af3/legends", "tables", "supplement")    # spec section 4: where cited_by is grepped
VERBATIM_COLS = {"context_1", "context_2", "context_3", "sentence", "verbatim", "first_context", "probe_context"}
HASH_CAP = 256 * 1024 * 1024    # larger files: md5 comes from the pin that records it (md5_basis says so)
PAIR_FRAME_EXT = (".tsv", ".csv")

NUM_RX = re.compile(r"(?<![\w.])[+\-−]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+\-−]?\d+)?(?![\w])")
QUOTE_RX = re.compile("[\"“]([^\"”]{40,}?)[\"”]")
PATH_RX = re.compile(r"(?<![A-Za-z0-9])([A-Za-z]):[/\\]|/Users/[A-Za-z]|Downloads[/\\]CKM")    # any drive-letter or home path (F21; non-blocking, plan l.155)


class Refusal(Exception):
    pass


# ------------------------------------------------------------------------------------------------------------------ helpers
def norm(p):
    return str(p).replace("\\", "/")


def md5_file(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def md5_bytes(b):
    return hashlib.md5(b).hexdigest()


def read_table(path, delim="\t"):
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter=delim))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def read_dicts(path, delim="\t"):
    hdr, rows = read_table(path, delim)
    return hdr, [dict(zip(hdr, r + [""] * (len(hdr) - len(r)))) for r in rows]


def write_table(path, hdr, rows, delim="\t"):
    """Deterministic: UTF-8, LF, QUOTE_MINIMAL. rows are lists (or dicts keyed by hdr)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO(newline="")
    w = csv.writer(buf, delimiter=delim, lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    w.writerow(hdr)
    for r in rows:
        w.writerow([r[c] for c in hdr] if isinstance(r, dict) else r)
    path.write_bytes(buf.getvalue().encode("utf-8"))


def delim_of(path):
    s = str(path).lower()
    return "\t" if s.endswith(".tsv") else "," if s.endswith(".csv") else None


def shape(path):
    """(rows, columns) for a delimited file, ('', '') otherwise. rows = data rows."""
    d = delim_of(path)
    if d is None:
        return "", ""
    hdr, rows = read_table(path, d)
    return str(len(rows)), str(len(hdr))


def sig_digits(tok):
    m = re.sub(r"[eE].*$", "", tok).lstrip("+-−")
    digits = m.replace(".", "").lstrip("0")
    digits = digits.rstrip("0") if "." not in m else digits.rstrip("0")
    return len(digits)


def tofloat(tok):
    return float(tok.replace("−", "-"))


def canon(tok, k=None):
    v = abs(tofloat(tok))
    return format(v, ".%dg" % (k if k else max(sig_digits(tok), 1)))


def fp_forms(tok):
    """The strings under which a number may reappear: the number itself and its roundings at 4..12 significant digits."""
    try:
        v = abs(tofloat(tok))
    except ValueError:
        return set()
    s = sig_digits(tok)
    if s < 4 or v == 0:
        return set()
    return {format(v, ".%dg" % k) for k in range(4, min(s, 12) + 1)}


def list_files(root):
    root = Path(root)
    return sorted(norm(p.relative_to(root)) for p in root.rglob("*") if p.is_file())


# ------------------------------------------------------------------------------------------------------------------ context
class Checks:
    def __init__(self):
        self.items = []

    def fail(self, cid, msg):
        self.items.append((cid, "FAIL", msg))

    def warn(self, cid, msg):
        self.items.append((cid, "WARN", msg))

    def info(self, msg):
        self.items.append(("I0", "INFO", msg))

    def fail_ids(self):
        return sorted({c for c, s, _ in self.items if s == "FAIL"})

    def failures(self):
        return [(c, m) for c, s, m in self.items if s == "FAIL"]

    def raise_if_failed(self, where):
        f = self.failures()
        if f:
            lines = ["%s: %d failure(s)" % (where, len(f))] + ["  [%s] %s" % (c, m) for c, m in f[:60]]
            if len(f) > 60:
                lines.append("  ... %d more" % (len(f) - 60))
            raise Refusal("\n".join(lines))


_MD5_CACHE = {}     # (path, size, mtime_ns) -> md5; shared by every Ctx of the process (the key changes when the file does)


class Ctx:
    def __init__(self, v3, release, base=None):
        self.v3 = Path(v3).resolve()
        self.ckm = self.v3.parent
        self.r2 = self.v3 / "R2_rewrite"
        self.release = Path(release).resolve()
        self.inputs = self.release / "inputs"
        self.base = Path(base).resolve() if base else self.release
        self.out = self.base / "out"
        self.withheld = self.base / "withheld"
        self.build = self.base / "build"
        self._md5 = _MD5_CACHE
        self.chk = Checks()

    def resolve(self, src):
        src = norm(src)
        if src.startswith("CKM papers/"):
            return self.ckm / src[len("CKM papers/"):]
        if re.match(r"^[A-Za-z]:/", src):
            return Path(src)
        return self.v3 / src

    def md5(self, path):
        path = Path(path)
        st = path.stat()
        key = (str(path), st.st_size, st.st_mtime_ns)
        if key not in self._md5:
            self._md5[key] = md5_file(path)
        return self._md5[key]

    def rel(self, path):
        p = Path(path).resolve()
        try:
            return norm(p.relative_to(self.v3))
        except ValueError:
            try:
                return "CKM papers/" + norm(p.relative_to(self.ckm))
            except ValueError:
                return norm(p)


# ------------------------------------------------------------------------------------------------------------------ declared tables
def load_ri(ctx):
    p = ctx.inputs / "RELEASE_INPUTS.tsv"
    if not p.is_file():
        raise Refusal("inputs/RELEASE_INPUTS.tsv is missing (%s)" % p)
    hdr, rows = read_dicts(p)
    if hdr != RI_COLS:
        raise Refusal("RELEASE_INPUTS.tsv header differs from the declared columns: %s" % hdr)
    return rows


def write_ri(ctx, rows):
    rows = sorted(rows, key=lambda r: (r["origin"] != "curated", r["input_id"]))
    write_table(ctx.inputs / "RELEASE_INPUTS.tsv", RI_COLS, rows)


def load_rules(ctx):
    p = ctx.inputs / "SOURCE_CLASS_RULES.tsv"
    if not p.is_file():
        raise Refusal("inputs/SOURCE_CLASS_RULES.tsv is missing")
    hdr, rows = read_dicts(p)
    if hdr != RULE_COLS:
        raise Refusal("SOURCE_CLASS_RULES.tsv header differs: %s" % hdr)
    for r in rows:
        if r["class"] not in RULE_CLASSES:
            raise Refusal("rule %s: class %r is not one of %s" % (r["rule_id"], r["class"], sorted(RULE_CLASSES)))
        if not r["reason"].strip():
            raise Refusal("rule %s has no reason (a rule without a reason is a silent default)" % r["rule_id"])
        if r["basis"] not in RULE_BASES:
            raise Refusal("rule %s: basis %r is not one of %s" % (r["rule_id"], r["basis"], sorted(RULE_BASES)))
    return rows


def classify_source(rules, path):
    """First rule whose glob matches the normalised source path; None if no rule matches."""
    p = norm(path)
    for r in rules:
        if fnmatch.fnmatchcase(p, r["pattern"]):
            return r
    return None


def load_declared(ctx, name, required=True):
    """A hand-declared table under inputs/ (DECLARED): header checked, never rewritten by any mode."""
    p = ctx.inputs / name
    if not p.is_file():
        if required:
            raise Refusal("inputs/%s is missing (a hand-declared table; no build runs without it)" % name)
        return []
    hdr, rows = read_dicts(p)
    if hdr != DECLARED[name]:
        raise Refusal("%s header differs from the declared columns: %s" % (name, hdr))
    return rows


# ------------------------------------------------------------------------------------------------------------------ providers and licence classes (F08, F23, F31)
def provider_tables(ctx):
    """(aliases, ids in declared order, declarations, declarations by id). Declarations map a source-path glob to the providers whose data its values are computed from."""
    aliases = load_declared(ctx, "PROVIDER_ALIASES.tsv")
    ids = [a["provider_id"] for a in aliases]
    if len(set(ids)) != len(ids) or not all(re.fullmatch(r"[A-Za-z0-9._-]+", i) for i in ids):
        raise Refusal("PROVIDER_ALIASES.tsv: provider ids must be unique and free of '+' and spaces")
    for a in aliases:
        try:
            re.compile(a["aliases"])
        except re.error as e:
            raise Refusal("PROVIDER_ALIASES.tsv %s: bad alias regex (%s)" % (a["provider_id"], e))
    decls = load_declared(ctx, "PROVIDER_DECLARATIONS.tsv")
    byid = {d["decl_id"]: d for d in decls}
    if len(byid) != len(decls):
        raise Refusal("PROVIDER_DECLARATIONS.tsv: duplicate decl_id")
    for d in decls:
        toks = [] if d["providers"] in ("", "none") else d["providers"].split("+")
        bad = [t for t in toks if t not in ids]
        if bad:
            raise Refusal("PROVIDER_DECLARATIONS.tsv %s names provider(s) %s that PROVIDER_ALIASES.tsv does not declare" % (d["decl_id"], bad))
        for i in [x.strip() for x in d["inherits"].split(",") if x.strip()]:
            if i not in byid:
                raise Refusal("PROVIDER_DECLARATIONS.tsv %s inherits unknown decl_id %s" % (d["decl_id"], i))
        if not d["evidence"].strip():
            raise Refusal("PROVIDER_DECLARATIONS.tsv %s has no evidence" % d["decl_id"])
    return aliases, ids, decls, byid


def decl_closure(d, byid, seen=()):
    if d["decl_id"] in seen:
        raise Refusal("PROVIDER_DECLARATIONS.tsv: inherits cycle at %s" % d["decl_id"])
    out = set() if d["providers"] in ("", "none") else set(d["providers"].split("+"))
    for i in [x.strip() for x in d["inherits"].split(",") if x.strip()]:
        out |= decl_closure(byid[i], byid, seen + (d["decl_id"],))
    return out


def decl_for_path(decls, path):
    p = norm(path)
    for d in decls:
        if fnmatch.fnmatchcase(p, d["pattern"]):
            return d
    return None


def class_for(provs, ids):
    """own_derived only when no provider supplied data the values are computed from; providers joined by '+' in PROVIDER_ALIASES order."""
    return "own_derived" if not provs else "third_party_derived:" + "+".join(i for i in ids if i in provs)


def provider_set_of_class(cls):
    return {t for c in cls.split("|") if c.startswith("third_party_derived:") for t in c.split(":", 1)[1].split("+")}


def unit_of_frame(release_path):
    m = re.match(r"source_data/figures/([^/]+)/", release_path)
    return m.group(1) if m else None


def unit_provider_sets(ctx, decls, byid, chk):
    """Providers per manifest figure id: the union over its non-SUPERSEDED manifest source files (F31: derived from the manifest, never typed per unit)."""
    out = {}
    for m in manifest_rows(ctx):
        if m["disposition"] == "SUPERSEDED":
            continue
        d = decl_for_path(decls, m["source_file"])
        if d is None:
            chk.fail("D3", "manifest line %d: source %s of %s matches no PROVIDER_DECLARATIONS row (no default)" % (m["_line"], m["source_file"], m["figure"]))
            continue
        out.setdefault(m["figure"], set()).update(decl_closure(d, byid))
    return out


def expected_licences(ctx, rows, chk):
    """{input_id: licence_class} derived from the provider declarations: a copy / input from its source path, a frame from its unit's manifest sources."""
    aliases, ids, decls, byid = provider_tables(ctx)
    ups = unit_provider_sets(ctx, decls, byid, chk)
    mids = dict(active_units(ctx))
    out = {}
    for r in rows:
        if r["role"] == "module_output":
            continue
        if r["origin"] == "derived_frames":
            u = unit_of_frame(r["release_path"])
            if u not in mids:
                chk.fail("D3", "%s: frame unit %r is not an active unit" % (r["input_id"], u))
                continue
            if mids[u] not in ups:
                chk.fail("D3", "%s: unit %s (manifest id %s) has no declared provider set (no manifest source rows)" % (r["input_id"], u, mids[u]))
                continue
            out[r["input_id"]] = class_for(ups[mids[u]], ids)
        else:
            d = decl_for_path(decls, r["source_path"])
            if d is None:
                chk.fail("D3", "%s: source %s matches no PROVIDER_DECLARATIONS row (no default)" % (r["input_id"], r["source_path"]))
                continue
            out[r["input_id"]] = class_for(decl_closure(d, byid), ids)
    return out


def check_derived_licences(ctx, rows, chk):
    for iid, want in expected_licences(ctx, rows, chk).items():
        have = [r["licence_class"] for r in rows if r["input_id"] == iid][0]
        if have != want:
            chk.fail("L3", "%s: licence_class %r differs from the class derived from PROVIDER_DECLARATIONS (%r); run --pin" % (iid, have, want))


def check_module_classes(ctx, res, chk):
    """A module-reported class names only declared providers (no placeholder such as unspecified_East_Asian_resource)."""
    _, ids, _, _ = provider_tables(ctx)
    for d in res:
        for c in d["licence_class"].split("|"):
            if c.startswith("third_party_derived:"):
                bad = [t for t in c.split(":", 1)[1].split("+") if t not in ids]
                if bad:
                    chk.fail("L4", "%s: licence_class names %s, which PROVIDER_ALIASES.tsv does not declare (resolve the provider from its files)" % (d["release_path"], bad))


def operative_surfaces(ctx):
    """The text surfaces cited_by is grepped over (spec section 4): the R2 prose / Table 3 / SM / legends under sections_R2, the AF3 legends, tables/ and supplement/.
    Code, registries and the manifest are not surfaces. Returns {path relative to R2_rewrite: text}."""
    sys.path.insert(0, str(ctx.r2 / "scripts"))
    try:
        sys.modules.pop("P31R_gate", None)
        import P31R_gate as G
        files = G.r2_surface_files()
    finally:
        sys.path.pop(0)
    out = {}
    for f in files:
        rel = norm(Path(f).resolve().relative_to(ctx.r2.resolve()))
        if rel.endswith(".md") and any(rel.startswith(d + "/") for d in OPERATIVE_DIRS):
            out[rel] = Path(f).read_text(encoding="utf-8", errors="replace")
    return out


def table3_ids(ctx):
    """Stage 7 A2.4c (S2-D51, review R18): Table 3's repository column moved verbatim out of the manuscript table into
    sections_R2/P31R_Table3_repository_rows.tsv (claim, repository_rows), which the release README prints; the source ids are read from it."""
    t = (ctx.r2 / "sections_R2" / "P31R_Table3_repository_rows.tsv").read_text(encoding="utf-8")
    lines = t.splitlines()
    assert lines and lines[0].split("\t") == ["claim", "repository_rows"], "Table 3 repository map: unexpected header"
    ids = set()
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) == 2:
            ids |= set(re.findall(r"\b(?:W\d+[ab]?|INT|BI)\b|extchecks", cells[1]))
    assert ids, "Table 3 repository map names no source id: the D6 check would run against an empty set"
    return ids


def check_cite_tokens(ctx, rows, chk):
    """cite_keys are plain tokens: the file's own name or a Table 3 source id; no number or phrase keys (F19)."""
    t3 = table3_ids(ctx)
    for r in rows:
        if r["role"] not in ("copy", "module_output") or not r["cite_keys"]:
            continue
        names = {Path(r["release_path"]).name}
        for tok in r["cite_keys"].split("|"):
            if tok not in names and tok not in t3 and not re.fullmatch(r"[A-Za-z0-9_.-]+\.(?:tsv|csv|json|md|R|log)", tok):
                chk.fail("D6", "%s: cite_keys token %r is neither a file name nor a Table 3 source id %s" % (r["input_id"], tok, sorted(t3)))


COUNT_RX = re.compile(r"\b\d[\d,]*(?:\s+[A-Za-z-]+){0,2}\s+(?:rows?|columns?|pairs?|files?|cells?|characters?|units?|genes?|figures?)\b|\b\d+-columns?\b|\b\d{1,3},\d{3}\b")


def check_notes(ctx, rows, chk):
    """F14 / F26: a curated note carries no typed count and names only release paths that exist."""
    known = {r["release_path"] for r in rows if r["release_path"]}
    for r in rows:
        n = r["note"]
        m = COUNT_RX.search(n) if r["origin"] == "curated" else None
        if m:
            chk.fail("D4", "%s: note carries a typed count %r (generate it from the file or drop it)" % (r["input_id"], m.group(0)))
        for m in re.finditer(r"(?:out/|withheld/|ships as )((?:[A-Za-z0-9_]+/)+[A-Za-z0-9_.-]+\.(?:tsv|csv|json))", n):
            if m.group(1) not in known:
                chk.fail("D5", "%s: note names the release path %s, which is not a declared release path" % (r["input_id"], m.group(1)))


def md_body_rows(path):
    """Cells of the body rows of the one markdown table in a file (first pipe table, after the |---| separator)."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    sep = [i for i, l in enumerate(lines) if re.fullmatch(r"\|(?:\s*:?-{3,}:?\s*\|)+", l.strip())]
    if len(sep) != 1:
        raise Refusal("%s: expected one table separator, found %d" % (path, len(sep)))
    body = []
    for l in lines[sep[0] + 1:]:
        if not l.startswith("|"):
            break
        body.append([c.strip() for c in l.strip().strip("|").split("|")])
    return body


def build_providers(ctx, ri, chk):
    """The text of inputs/PROVIDERS.tsv (spec: step 4 provider list, user decision 2). Sources: Table 1's resources, the Declarations acknowledgement list, the
    providers named in the operative R2 surfaces, and the providers in the licence classes. A provider named or used but in neither Table 1 nor the Declarations must
    sit on PROVIDERS_PENDING_DECLARATIONS.tsv (to be added to the Declarations at step 5) or the build fails (D7). terms_status (RC1) is the provider's status in
    inputs/PROVIDER_TERMS.tsv: confirmed_from_official_page | inferred | not_found (a provider with no row there fails: D9)."""
    aliases, ids, decls, byid = provider_tables(ctx)
    pend = {p["provider_id"]: p for p in load_declared(ctx, "PROVIDERS_PENDING_DECLARATIONS.tsv")}
    for pid in pend:
        if pid not in ids:
            raise Refusal("PROVIDERS_PENDING_DECLARATIONS.tsv names %s, which PROVIDER_ALIASES.tsv does not declare" % pid)
    surf = operative_surfaces(ctx)
    t1 = " || ".join(r[1] for r in md_body_rows(ctx.r2 / "tables" / "P31R_Table1.md"))
    dec = (ctx.r2 / "sections" / "P31R_Declarations.md").read_text(encoding="utf-8")
    m = re.search(r"### Acknowledgements\n(.*?)(?=\n###|\Z)", dec, re.S)
    if not m:
        raise Refusal("P31R_Declarations.md has no '### Acknowledgements' section")
    ack = m.group(1)
    check_provider_coverage(ctx, aliases, surf, chk)
    exp = expected_licences(ctx, ri, chk)
    used = set()
    for c in exp.values():
        used |= provider_set_of_class(c)
    terms = {t["provider_id"].strip(): t["status"].strip() for t in read_dicts(ctx.inputs / "PROVIDER_TERMS.tsv")[1]}
    rows = []
    for a in aliases:
        pid, rx = a["provider_id"], re.compile(a["aliases"])
        if terms.get(pid) not in ("confirmed_from_official_page", "inferred", "not_found"):
            chk.fail("D9", "provider %s has no readable status in inputs/PROVIDER_TERMS.tsv (%r)" % (pid, terms.get(pid)))
        in_t1, in_dec = bool(rx.search(t1)), bool(rx.search(ack))
        named = sorted({Path(n).name for n, t in surf.items() if rx.search(t)})
        if (named or pid in used) and not (in_t1 or in_dec) and pid not in pend:
            chk.fail("D7", "provider %s is %s but is in neither Table 1 nor the Declarations acknowledgements: add it to PROVIDERS_PENDING_DECLARATIONS.tsv (for step 5) "
                     "or to the Declarations" % (pid, "named in " + named[0] if named else "in a licence class"))
        if pid in pend and (in_t1 or in_dec):
            chk.fail("D7", "provider %s is on PROVIDERS_PENDING_DECLARATIONS.tsv but Table 1 / the Declarations already name it: remove the pending row" % pid)
        rows.append([pid, a["display_name"], a["aliases"], a["role"], ";".join((["Table1"] if in_t1 else []) + (["Declarations"] if in_dec else []) + named),
                     "yes" if in_dec else ("pending_step5" if pid in pend else "no"), "third_party_derived" if pid in used else "none", terms.get(pid, "not_read"), a["release_handling"]])
    buf = io.StringIO(newline="")
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(PROVIDER_COLS)
    w.writerows(rows)
    return buf.getvalue()


# Terms that name a database, distributor or annotation source. Every occurrence on an operative surface must itself match a declared alias (D8): a source that the
# surfaces name but PROVIDER_ALIASES.tsv lacks fails the build. Extend this list when a surface starts naming a new kind of source.
SOURCE_LEXICON = re.compile(r"ChEMBL|DGIdb|UniProt|Ensembl|Open ?Targets|OpenGWAS|IEU OpenGWAS|GWAS Catalog|GWAS Atlas|pheweb\.jp|GTEx|gnomAD|DrugBank|Reactome|Human Protein Atlas|dbSNP|ClinVar|"
                            r"PharmGKB|eQTLGen|PhenoScanner|DisGeNET|Pan-UK Biobank|Million Veteran", re.I)
ACCESSION_RX = re.compile(r"GCST\d{6,}|ieu-[ab]-\d+|hum\d{4}\.v\d+")


def check_provider_coverage(ctx, aliases, surf, chk):
    """D8: a provider named in Table S2 or the supplementary material, and missing from PROVIDER_ALIASES.tsv, fails the build.
    (a) every Table S2 Resource (or File or accession) cell matches an alias; (b) every database / distributor / annotation term on an operative surface matches an alias;
    (c) every line that cites a GWAS accession (GCST, ieu-, hum) on Table S2 or an SM surface also names a declared provider."""
    rxs = [re.compile(a["aliases"]) for a in aliases]
    cite = re.compile(r"\[@[^\]]+\]|\[\[L\d+[a-z]?\]\]")
    for cells in md_body_rows(ctx.r2 / "tables" / "P31R_TableS2.md"):
        res = cite.sub("", cells[1])
        if not any(r.search(res + " " + cells[2]) for r in rxs):     # the Resource cell, or the File or accession cell where the resource is only a description (EGFR_EAS_2P2M)
            chk.fail("D8", "Table S2 resource %r names no provider declared in PROVIDER_ALIASES.tsv: declare it (and list it in PROVIDERS_PENDING_DECLARATIONS.tsv if the Declarations lack it)" % res)
    for name, text in surf.items():
        base = Path(name).name
        for i, line in enumerate(text.splitlines(), 1):
            for m in SOURCE_LEXICON.finditer(line):
                if not any(r.search(m.group(0)) for r in rxs):
                    chk.fail("D8", "%s l.%d names the source %r, which PROVIDER_ALIASES.tsv does not declare" % (base, i, m.group(0)))
            if ("_SM_" in base or base.startswith("P31R_SNote") or base == "P31R_TableS2.md") and ACCESSION_RX.search(line) and not any(r.search(line) for r in rxs):
                chk.fail("D8", "%s l.%d cites an accession but names no declared provider" % (base, i))


def check_providers(ctx, ri, chk):
    want = build_providers(ctx, ri, chk)
    p = ctx.inputs / "PROVIDERS.tsv"
    if not p.is_file() or p.read_bytes() != want.encode("utf-8"):
        chk.fail("P4", "inputs/PROVIDERS.tsv is missing or stale against the declarations and surfaces: run --providers")
    return want


def validate_ri(ctx, rows, chk):
    ids = collections.Counter(r["input_id"] for r in rows)
    for i, n in ids.items():
        if n > 1:
            chk.fail("D2", "duplicate input_id %s" % i)
    seen = {}
    for r in rows:
        i = r["input_id"]
        if r["role"] not in ROLES:
            chk.fail("D1", "%s: role %r not in %s" % (i, r["role"], sorted(ROLES)))
        if r["origin"] not in ORIGINS:
            chk.fail("D1", "%s: origin %r" % (i, r["origin"]))
        if r["ship_status"] not in SHIPS and r["role"] != "input":
            chk.fail("D1", "%s: ship_status %r not in %s" % (i, r["ship_status"], sorted(SHIPS)))
        lc = r["licence_class"]
        if lc == "UNCLASSIFIED" or not (licence_ok(lc) or (r["role"] == "module_output" and lc == "module_reported")):
            chk.fail("L1", "%s: licence_class %r is not own_derived / third_party_derived:<provider> / not_redistributable" % (i, lc))
        if r["role"] in ("copy", "module_output"):
            if not r["release_path"] or r["release_path"].startswith("/") or ".." in r["release_path"].split("/"):
                chk.fail("D1", "%s: release_path %r" % (i, r["release_path"]))
            key = (r["release_path"], r["ship_status"])
            if key in seen:
                chk.fail("D2", "release_path %s (%s) declared by %s and %s" % (r["release_path"], r["ship_status"], seen[key], i))
            seen[key] = i
        if r["role"] in ("copy", "input") and not r["source_path"]:
            chk.fail("D1", "%s: a %s row needs a source_path" % (i, r["role"]))
        if r["role"] == "copy" and r["ship_status"] == "withheld_pending_terms" and r["withhold_rule"].strip():
            chk.fail("D1", "%s: a whole-file withheld row carries a rule" % i)
        try:
            parse_rule(r["withhold_rule"])
        except Refusal as e:
            chk.fail("D1", "%s: %s" % (i, e))


# ------------------------------------------------------------------------------------------------------------------ withhold-rule DSL
# rule  := clause (' ;; ' clause)*
# clause:= 'tpmi|h9|ot' ':' where '=>' action     (masking / assertion)
#        | 'h9:@H9'                               (columns and source value come from the hold in P31R_ACTIVE_FIGURES.json)
#        | 'labels:' cols | 'exempt:' cols | 'allow:' tokens | 'guard:quote12'
# where := COL '~' REGEX | COL '=' V1 '|' V2 | '*~' REGEX    (ANY cell matches)
# action:= mask(c1,c2) | keep(c1,c2)  (mask every other column) | require(c1,c2)  (assertion: the cells must already carry the sentinel)
CLAUSE_RX = re.compile(r"^(tpmi|h9|ot):(?P<where>[^=]+?(?:=[^>]*?)?)=>(?P<act>mask|keep|require)\((?P<cols>[^)]*)\)$")
WHERE_RX = re.compile(r"^(?P<col>\*|[A-Za-z0-9_.]+)(?P<op>~|=)(?P<arg>.*)$")


class Clause:
    def __init__(self, cls, col=None, op=None, arg=None, act=None, cols=None, special=None):
        self.cls, self.col, self.op, self.arg, self.act, self.cols, self.special = cls, col, op, arg, act, cols, special
        self.rx = re.compile(arg, re.I) if op == "~" else None
        self.hits = 0

    def matches(self, rec):
        if self.col == "*":
            return any(self.rx.search(v) for v in rec.values())
        v = rec.get(self.col)
        if v is None:
            return False
        return bool(self.rx.search(v)) if self.op == "~" else v in self.arg.split("|")


def parse_rule(s):
    s = (s or "").strip()
    if not s:
        return []
    out = []
    for part in [x.strip() for x in s.split(" ;; ")]:
        if part == "h9:@H9":
            out.append(Clause("h9", special="H9"))
        elif part == "h9:absent":
            out.append(Clause("h9", special="absent"))
        elif part.startswith(("labels:", "exempt:", "allow:")):
            k, v = part.split(":", 1)
            out.append(Clause(k, cols=[c.strip() for c in v.split(",") if c.strip()]))
        elif part == "guard:quote12":
            out.append(Clause("guard", special="quote12"))
        else:
            m = CLAUSE_RX.match(part)
            if not m:
                raise Refusal("unparseable withhold clause %r" % part)
            w = WHERE_RX.match(m.group("where"))
            if not w:
                raise Refusal("unparseable where %r" % m.group("where"))
            out.append(Clause(m.group(1), w.group("col"), w.group("op"), w.group("arg"), m.group("act"),
                              [c.strip() for c in m.group("cols").split(",") if c.strip()]))
    return out


def load_holds(ctx):
    d = json.loads((ctx.r2 / "figures" / "af3" / "P31R_ACTIVE_FIGURES.json").read_text(encoding="utf-8"))
    return {h["id"]: h for h in d.get("holds", [])}


def h9_spec(ctx):
    h = load_holds(ctx).get("H9")
    if not h:
        raise Refusal("P31R_ACTIVE_FIGURES.json carries no hold H9")
    return h


def masking_clauses(clauses):
    return [c for c in clauses if c.cls in ("tpmi", "h9", "ot") and (c.special == "H9" or c.act in ("mask", "keep"))]


def h9_declared_absent(clauses):
    return any(c.special == "absent" for c in clauses)


def apply_rules(ctx, hdr, rows, clauses, label, chk):
    """Mask rows in place. Returns the original cells that were masked as [(cls, col, original)] (for the fingerprint)."""
    masked = []
    for cl in masking_clauses(clauses):
        recs_hit = 0
        if cl.special == "H9":
            h = h9_spec(ctx)
            src_cols = [c for c in h["source_columns"] if c in hdr]
            mag = [c for c in h["magnitude_columns"] if c in hdr]
            if mag and not src_cols:
                chk.fail("H2", "%s: carries H9 magnitude columns %s but no source column %s" % (label, mag, h["source_columns"]))
                continue
            for row in rows:
                rec = dict(zip(hdr, row))
                if src_cols and rec[src_cols[0]] == h["source_value"]:
                    recs_hit += 1
                    for c in mag:
                        i = hdr.index(c)
                        if row[i] not in ("", "NA", SENT_H9):
                            masked.append(("h9", c, row[i]))
                            row[i] = SENT_H9
            cl.hits = recs_hit
        else:
            bad = [c for c in cl.cols if c not in hdr]
            if bad:
                chk.fail("S5", "%s: rule names column(s) %s that the file does not have" % (label, bad))
                continue
            for row in rows:
                rec = dict(zip(hdr, row))
                if cl.matches(rec):
                    recs_hit += 1
                    tgt = cl.cols if cl.act == "mask" else [c for c in hdr if c not in cl.cols]
                    for c in tgt:
                        i = hdr.index(c)
                        if row[i] != SENTINELS[cl.cls]:
                            masked.append((cl.cls, c, row[i]))
                            row[i] = SENTINELS[cl.cls]
            cl.hits = recs_hit
        if recs_hit == 0:
            chk.fail("S4", "%s: masking clause %s:%s matched no row (a rule that never fires is not a control)" %
                     (label, cl.cls, cl.special or "%s%s%s" % (cl.col, cl.op, cl.arg)))
    return masked


def require_check(hdr, rows, clauses, label, chk):
    """require(...) clauses: the matching rows must carry the sentinel of their class in every listed column."""
    for cl in clauses:
        if cl.cls in ("tpmi", "h9", "ot") and cl.act == "require":
            bad = [c for c in cl.cols if c not in hdr]
            if bad:
                chk.fail("S3", "%s: require() names column(s) %s that the file does not have" % (label, bad))
                continue
            hit = 0
            for i, row in enumerate(rows, 2):
                rec = dict(zip(hdr, row))
                if cl.matches(rec):
                    hit += 1
                    for c in cl.cols:
                        if rec[c] != SENTINELS[cl.cls]:
                            chk.fail("S3", "%s line %d column %s: expected %s (%s-derived), found %r" %
                                     (label, i, c, SENTINELS[cl.cls], cl.cls, rec[c][:40]))
            if hit == 0:
                chk.fail("S4", "%s: require clause %s:%s%s%s matched no row" % (label, cl.cls, cl.col, cl.op, cl.arg))


# ------------------------------------------------------------------------------------------------------------------ pins
def ri_rows_with_source(rows):
    return [r for r in rows if r["source_path"]]


def verify_pins(ctx, rows, chk):
    for r in ri_rows_with_source(rows):
        p = ctx.resolve(r["source_path"])
        if not p.is_file():
            chk.fail("P1", "%s: source %s is absent" % (r["input_id"], r["source_path"]))
            continue
        if not re.fullmatch(r"[0-9a-f]{32}", r["source_md5"] or ""):
            chk.fail("P1", "%s: %s is not pinned (run --pin)" % (r["input_id"], r["source_path"]))
            continue
        if p.stat().st_size > HASH_CAP:
            chk.fail("P1", "%s: %s exceeds the hash cap; declare it as a manifest-pinned external source instead" % (r["input_id"], r["source_path"]))
            continue
        now = ctx.md5(p)
        if now != r["source_md5"]:
            chk.fail("P2", "%s: %s md5 %s != pinned %s (run --pin only after reading why it changed)" %
                     (r["input_id"], r["source_path"], now[:8], r["source_md5"][:8]))


# ------------------------------------------------------------------------------------------------------------------ units and frames
def active_units(ctx):
    d = json.loads((ctx.r2 / "figures" / "af3" / "P31R_ACTIVE_FIGURES.json").read_text(encoding="utf-8"))
    units = []
    for f in d["figures"]:
        fid = f["id"]
        mid = fid if fid.startswith("Figure_") else "AF3_S" + fid[len("FigureS"):]
        units.append((fid, mid))
        for p in range(2, int(f.get("pages", 1)) + 1):
            units.append(("%s_p%d" % (fid, p), mid))
    return units


def fig2_writer_files(ctx):
    """Frames of Figure 2 come from the writer calls in P31R_Fig2.R, not from a glob: dfd <- file.path(P31R_FIGS, "data", "Figure_2") and
    every file.path(dfd, "<name>")."""
    src = (ctx.r2 / "figures" / "scripts" / "P31R_Fig2.R").read_text(encoding="utf-8")
    m = re.search(r'dfd\s*<-\s*file\.path\(P31R_FIGS((?:,\s*"[^"]+")+)\)', src)
    if not m:
        raise Refusal("P31R_Fig2.R: cannot find the dfd <- file.path(P31R_FIGS, ...) line that names the frame folder")
    parts = re.findall(r'"([^"]+)"', m.group(1))
    names = sorted(set(re.findall(r'file\.path\(dfd,\s*"([^"]+)"\)', src)))
    if not names:
        raise Refusal("P31R_Fig2.R: no file.path(dfd, ...) writer call found")
    return "R2_rewrite/figures/" + "/".join(parts), names


def pending_expected(ctx, unit):
    """Frame names a pending unit will write: from the generator if its edit is applied, else from the staged-edit record. None if neither."""
    gen = ctx.r2 / "figures" / "scripts" / ("P31R_" + unit.replace("Figure_", "Fig") + ".R")
    if gen.is_file():
        src = gen.read_text(encoding="utf-8")
        m = re.search(r"frames\s*=\s*(?:I\()?list\(", src)
        if m:
            tail, depth, seg, segs = src[m.end():], 1, "", []
            for ch in tail:
                if ch in "([{":
                    depth += 1
                elif ch in ")]}":
                    depth -= 1
                    if depth == 0:
                        break
                if ch == "," and depth == 1:
                    segs.append(seg)
                    seg = ""
                else:
                    seg += ch
            segs.append(seg)
            names = [x.group(1).strip("\"'") for x in (re.match(r"\s*([A-Za-z0-9_\"']+)\s*=", t) for t in segs) if x]
            if names:
                return names, "generator P31R_%s.R" % unit.replace("Figure_", "Fig")
    rec = ctx.r2 / "figures" / "registry" / "STAGE5_STEP2_FRAMES_STAGED.md"
    if rec.is_file():
        names = sorted(set(re.findall(r"`%s__([A-Za-z0-9_]+)\.csv`" % re.escape(unit), rec.read_text(encoding="utf-8"))))
        if names:
            return names, "staged record STAGE5_STEP2_FRAMES_STAGED.md (edit not applied, not rendered)"
    return None, ""


def derive_frames(ctx):
    """([ {unit, name, source_path, derivation, release_path} ], [ {unit, expected, why} ], notes)"""
    reg_dir = ctx.r2 / "figures" / "registry"
    frames, pending, notes = [], [], []
    for unit, _ in active_units(ctx):
        rj = reg_dir / (unit + ".json")
        if not rj.is_file():
            raise Refusal("registry %s is missing for unit %s" % (rj, unit))
        reg = json.loads(rj.read_text(encoding="utf-8"))
        if unit == "Figure_1":
            frames.append(dict(unit=unit, name=rj.name, source_path=ctx.rel(rj), derivation="registry JSON (Figure 1 prints ledger counts only)"))
        elif unit == "Figure_2":
            folder, names = fig2_writer_files(ctx)
            for n in names:
                frames.append(dict(unit=unit, name=n, source_path=folder + "/" + n, derivation="writer calls in P31R_Fig2.R"))
        elif reg.get("frames"):
            for fr in reg["frames"]:
                frames.append(dict(unit=unit, name="%s__%s.csv" % (unit, fr), source_path="R2_rewrite/figures/af3/source_data/%s__%s.csv" % (unit, fr),
                                   derivation="registry %s.json frames list" % unit))
        else:
            ex, how = pending_expected(ctx, unit)
            pending.append(dict(unit=unit, expected=(";".join("%s__%s.csv" % (unit, n) for n in ex) + " [from " + how + "]") if ex else "unknown (registry has no frames list)",
                                why="the step 1c render writes the plotted frame; registry %s.json carries no frames list yet" % unit))
    for f in frames:
        f["release_path"] = "source_data/figures/%s/%s" % (f["unit"], f["name"])
    return frames, pending, notes


def sync_frame_rows(ctx, rows, frames):
    """Replace the derived rows by the current derivation, keeping hand-edited fields of a row whose source_path survives."""
    old = {r["source_path"]: r for r in rows if r["origin"] == "derived_frames"}
    keep = [r for r in rows if r["origin"] != "derived_frames"]
    new = []
    for i, f in enumerate(sorted(frames, key=lambda x: x["release_path"]), 1):
        o = old.get(f["source_path"])
        new.append(dict(input_id="FR%03d" % i, origin="derived_frames", role="copy", deliverable="source_data_figures",
                        release_path=f["release_path"], source_path=f["source_path"], source_md5=o["source_md5"] if o else "",
                        source_bytes=o["source_bytes"] if o else "",
                        licence_class=o["licence_class"] if o else "UNCLASSIFIED",
                        ship_status=o["ship_status"] if o else "ship", withhold_rule=o["withhold_rule"] if o else "",
                        producer=f["derivation"], cite_keys=o["cite_keys"] if o else f["name"],
                        cited_by=o["cited_by"] if o else "", note=o["note"] if o else ""))
    return keep + new, sorted(set(old) - {f["source_path"] for f in frames}), sorted({f["source_path"] for f in frames} - set(old))


# ------------------------------------------------------------------------------------------------------------------ cited_by
def surface_texts(ctx):
    """{file name: text} of the operative R2 text surfaces (operative_surfaces); the gate resolves its own root, so a scratch release dir still reads the real tree."""
    return {Path(rel).name: t for rel, t in operative_surfaces(ctx).items()}


def cited_by(texts, keys):
    """Operative surfaces that name the file or its Table 3 source id (spec section 4). keys = plain tokens joined by '|': a file name is matched literally, a
    source id as a whole word. No number or phrase keys exist (check_cite_tokens)."""
    if not keys:
        return ""
    toks = [k for k in keys.split("|") if k]
    rx = re.compile("|".join((re.escape(k) if "." in k else r"\b%s\b" % re.escape(k)) for k in toks))
    return ";".join(sorted(n for n, t in texts.items() if rx.search(t)))


# ------------------------------------------------------------------------------------------------------------------ leak scans
class Fingerprint:
    """Numbers that were masked as TPMI- or H9-derived, and the strings they may reappear as. Whole-number cells only: free text in a masked
    cell is masked by its rule but never fingerprinted (a threshold printed in a note is not a TPMI value).
    strong: roundings at 6..12 significant digits, flagged wherever they reappear (a coincidence at six digits is a ~1e-5 event).
    weak:   roundings at 3..5 digits, flagged only in a cell / line that also names the source (TPMI, FGFR4, 2P2M).
    A token is never flagged in the file and column it was masked from: that file's own rule is the control there, and q-values tie."""

    def __init__(self):
        self.strong = {"tpmi": {}, "h9": {}}
        self.weak = {"tpmi": {}, "h9": {}}

    def add(self, cls, cell, where):
        cell = cell.strip()
        if cls not in self.strong or not NUM_RX.fullmatch(cell):
            return
        try:
            v = abs(tofloat(cell))
        except ValueError:
            return
        s = sig_digits(cell)
        if v == 0 or s < 3:
            return
        for k in range(6, min(s, 12) + 1):
            self.strong[cls].setdefault(format(v, ".%dg" % k), where)
        for k in range(3, min(s, 5) + 1):
            self.weak[cls].setdefault(format(v, ".%dg" % k), where)

    def n(self):
        return {k: len(self.strong[k]) for k in self.strong}


_TPMI_SHIP_FORCE = None      # selftest only: True / False overrides the TIER_RULES switch for one control, None = read the file


def tpmi_ship_switch(release=None):
    if _TPMI_SHIP_FORCE is not None:
        return _TPMI_SHIP_FORCE
    """S2-D34(1): TIER_RULES.tsv row SW1 (kind switch, target TPMI_DERIVED) at 'ship' means a TPMI mention is no leak; the VALUE scan (fingerprint) still applies to the masked raw cells."""
    with open(Path(release or HERE) / "inputs" / "TIER_RULES.tsv", encoding="utf-8", newline="") as f:
        sw = [r for r in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE) if r["kind"] == "switch" and r["target"] == "TPMI_DERIVED"]
    if len(sw) != 1 or sw[0]["handling"].strip() not in ("ship", "withheld"):
        raise Refusal("TIER_RULES.tsv has no single TPMI_DERIVED switch with handling ship|withheld")
    return sw[0]["handling"].strip() == "ship"


def scan_tabular_cells(hdr, rows):
    for ln, row in enumerate(rows, 2):
        for ci, cell in enumerate(row):
            yield ln, (hdr[ci] if ci < len(hdr) else "?"), cell


def check_cell_tokens(fp, cell, label, ln, col, chk, allow):
    named = re.search(r"TPMI|FGFR4|2P2M", cell) is not None
    for m in NUM_RX.finditer(cell):
        tok = m.group(0)
        s = sig_digits(tok)
        if s < 3:
            continue
        try:
            v = abs(tofloat(tok))
        except ValueError:
            continue
        for cls in fp.strong:
            if s >= 6:
                c = format(v, ".%dg" % min(s, 12))
                hit = fp.strong[cls].get(c)
                tier = "strong"
            elif named:
                c = format(v, ".%dg" % s)
                hit = fp.weak[cls].get(c)
                tier = "weak (beside its source name)"
            else:
                continue
            if hit and hit != (label, col) and c not in allow:
                chk.fail("S1", "%s line %d column %s: %s-derived value %s, %s tier (masked at %s column %s)" %
                         (label, ln, col, cls.upper(), tok, tier, hit[0], hit[1]))


def scan_shipped(ctx, root, rule_of, fp, chk, stage2):
    """Scan every text file under root (= out/). rule_of: release_path -> clauses. Returns informational counts."""
    stats = collections.Counter()
    holds = h9_spec(ctx)
    mag = set(holds["magnitude_columns"])
    for rp in list_files(root):
        if rp in DOC_FILES or rp.startswith("code/"):
            continue
        path = Path(root) / rp
        clauses = rule_of.get(rp, [])
        labels = {c for cl in clauses if cl.cls == "labels" for c in cl.cols}
        exempt = {c for cl in clauses if cl.cls == "exempt" for c in cl.cols}
        allow = {t for cl in clauses if cl.cls == "allow" for t in cl.cols}
        guard = any(cl.cls == "guard" for cl in clauses)
        d = delim_of(rp)
        stats["files_scanned"] += 1
        try:
            raw = path.read_bytes()
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            stats["binary_skipped"] += 1
            continue
        pm = [m.group(0) for m in PATH_RX.finditer(text)]
        if pm:
            stats["path_files"] += 1
            cnt = collections.Counter(x.rstrip("/\\") for x in pm)
            chk.warn("Q3", "%s holds %d absolute-path hit(s) %s (step 4 decides whether to scrub; not blocking)" % (rp, len(pm), dict(sorted(cnt.items()))))
        if "[HOLD" in text:
            stats["hold_files"] += 1
            chk.items.append(("K1", "HOLD", rp))
        if d:
            hdr, rows = read_table(path, d)
            vc = VERBATIM_COLS & set(hdr)
            if vc:
                chk.fail("Q2", "%s has verbatim-text column(s) %s (copyright)" % (rp, sorted(vc)))
            for ln, col, cell in scan_tabular_cells(hdr, rows):
                if col in exempt or cell in (SENT_PENDING, SENT_H9):
                    continue
                if re.search(r"tpmi", cell, re.I) and col not in labels and "*" not in labels and not tpmi_ship_switch():
                    chk.fail("S2", "%s line %d column %s mentions TPMI outside a declared label column (%s)" % (rp, ln, col, cell[:50]))
                check_cell_tokens(fp, cell, rp, ln, col, chk, allow)
                if guard:
                    for m in QUOTE_RX.finditer(cell):
                        if len(m.group(1).split()) >= 12:
                            chk.fail("Q1", "%s line %d column %s: quoted run of %d words" % (rp, ln, col, len(m.group(1).split())))
            hset = [c for c in hdr if c in mag]
            if not hset:
                # F32: the hold is matched by ROW CONTENT too, not only by column names. A row that names the held resource (its source value, or '2P2M' in an
                # accession or a quantity id) may hold counts and directions only: a decimal in any other cell is a magnitude (H9 note).
                trx = re.compile(re.escape(holds["source_value"]) + "|2P2M")
                for ln, row in enumerate(rows, 2):
                    if any(trx.search(c) for c in row):
                        for ci, c in enumerate(row):
                            if hdr[ci] not in exempt and NUM_RX.fullmatch(c.strip()) and re.search(r"[.eE]", c):
                                chk.fail("H3", "%s line %d column %s: decimal %r on a row that names the held resource %s (hold H9: counts and directions only)" %
                                         (rp, ln, hdr[ci], c, holds["source_value"]))
            if hset:
                src = [c for c in holds["source_columns"] if c in hdr]
                if not src and h9_declared_absent(clauses):
                    if holds["source_value"] in text:
                        chk.fail("H1", "%s is declared free of %s (h9:absent) but the string occurs in it" % (rp, holds["source_value"]))
                elif not src:
                    chk.fail("H2", "%s carries H9 magnitude columns %s but no source column: hold cannot be checked" % (rp, hset))
                else:
                    si = hdr.index(src[0])
                    for ln, row in enumerate(rows, 2):
                        if row[si] == holds["source_value"]:
                            for c in hset:
                                v = row[hdr.index(c)]
                                if v not in ("", "NA", SENT_H9, SENT_PENDING):
                                    chk.fail("H1", "%s line %d: %s=%r on a %s row" % (rp, ln, c, v, holds["source_value"]))
        else:
            if re.search(r"tpmi", text, re.I) and "*" not in labels and not tpmi_ship_switch():
                chk.fail("S2", "%s mentions TPMI and is not a declared label file" % rp)
            for ln, line in enumerate(text.splitlines(), 1):
                check_cell_tokens(fp, line, rp, ln, "-", chk, allow)
    return stats


# ------------------------------------------------------------------------------------------------------------------ pair-keyed scan of TPMI- and Open Targets-derived cells (F04)
# The string scan above (S2) only sees cells that SAY "TPMI". The quantities below are derived from TPMI (or Open Targets) whatever the cell says: the FGFR4-T2D
# (EAS) external check's estimate, CI, p, n, class, condition-1 verdict and corroboration verdict, and the Open Targets development-stage columns. They are found by
# KEY (the pair label as each file writes it, DERIVED_PAIR_KEYS.tsv) and every row a key finds must be classified by a hand-declared DERIVED_ROWS.tsv entry; a
# 'derived' entry names the columns that must read the sentinel (derived_cols) or the only columns that may hold a value (label_cols). An unclassified pair row,
# an unmasked derived cell, and an entry that matches nothing each fail.
def _where(where):
    if where == "ALL":
        return lambda rec: True
    m = WHERE_RX.match(where)
    if not m:
        raise Refusal("DERIVED_ROWS.tsv: unparseable where %r" % where)
    col, op, arg = m.group("col"), m.group("op"), m.group("arg")
    if op == "~":
        rx = re.compile(arg, re.I)
        return (lambda rec: any(rx.search(v) for v in rec.values())) if col == "*" else (lambda rec: bool(rx.search(rec.get(col, ""))))
    vals = arg.split("|")
    return lambda rec: rec.get(col) in vals


def expand_cols(spec, hdr):
    """derived_cols / label_cols: comma list; a token may be an fnmatch pattern (ext_*) and a token starting '-' excludes. Returns (columns, tokens that matched nothing)."""
    toks = [t.strip() for t in spec.split(",") if t.strip()]
    inc, exc = [t for t in toks if not t.startswith("-")], [t[1:] for t in toks if t.startswith("-")]
    cols = [c for c in hdr if any(fnmatch.fnmatchcase(c, t) for t in inc) and not any(fnmatch.fnmatchcase(c, t) for t in exc)]
    return cols, [t for t in inc if not any(fnmatch.fnmatchcase(c, t) for c in hdr)]


def scan_derived(ctx, root, chk):
    keys = load_declared(ctx, "DERIVED_PAIR_KEYS.tsv")
    ents = load_declared(ctx, "DERIVED_ROWS.tsv")
    for k in keys:
        if k["class"] not in SENTINELS or k["target"] not in ("cell", "header"):
            raise Refusal("DERIVED_PAIR_KEYS.tsv %s: class/target" % k["key_id"])
        k["rx"] = re.compile(k["key_regex"], re.I)
    ids = collections.Counter(e["row_id"] for e in ents)
    for e in ents:
        if ids[e["row_id"]] > 1 or e["status"] not in ("derived", "not_derived") or e["class"] not in SENTINELS or not e["reason"].strip():
            raise Refusal("DERIVED_ROWS.tsv %s: duplicate id, status, class or empty reason" % e["row_id"])
        if e["status"] == "derived" and bool(e["derived_cols"].strip()) == bool(e["label_cols"].strip()):
            raise Refusal("DERIVED_ROWS.tsv %s: a derived entry needs exactly one of derived_cols / label_cols" % e["row_id"])
        e["match"] = (lambda rec: True) if e["where"] == "FILE" else _where(e["where"])
        e["hits"] = 0
    n_rows = 0
    for rp in list_files(root):
        mine = [e for e in ents if fnmatch.fnmatchcase(rp, e["release_glob"])]
        d = delim_of(rp)
        path = Path(root) / rp
        if d is None:
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
            for k in keys:
                if k["target"] == "cell" and k["rx"].search(text):
                    hit = [e for e in mine if e["class"] == k["class"] and e["where"] == "FILE"]
                    if not hit:
                        chk.fail("S8", "%s (not tabular) carries a %s key (%s) and no DERIVED_ROWS 'FILE' entry classifies it" % (rp, k["class"], k["key_id"]))
                    for e in hit:
                        e["hits"] += 1
            continue
        hdr, rows = read_table(path, d)
        for ln, row in enumerate(rows, 2):
            rec = dict(zip(hdr, row))
            for cls in sorted({k["class"] for k in keys}):
                kk = [k for k in keys if k["class"] == cls and ((k["target"] == "cell" and any(k["rx"].search(c) for c in row)) or
                                                                 (k["target"] == "header" and any(k["rx"].fullmatch(h) for h in hdr)))]
                if not kk:
                    continue
                n_rows += 1
                hit = [e for e in mine if e["class"] == cls and e["where"] != "FILE" and e["match"](rec)]
                if not hit:
                    chk.fail("S8", "%s line %d: row carries a %s-derived key (%s) and no DERIVED_ROWS entry classifies it" % (rp, ln, cls.upper(), kk[0]["key_id"]))
                for e in hit:
                    e["hits"] += 1
                    if e["status"] != "derived":
                        continue
                    sent = SENTINELS[cls]
                    if e["derived_cols"].strip():
                        cols, bad = expand_cols(e["derived_cols"], hdr)
                        if bad:
                            chk.fail("S9", "%s: entry %s names derived column(s) %s that the file does not have" % (rp, e["row_id"], bad))
                        e.setdefault("cols_seen", set()).update(cols)
                        for c in cols:
                            if rec[c] not in (sent, "", "NA"):
                                chk.fail("S9", "%s line %d column %s: %s-derived (%s) but not masked: %r" % (rp, ln, c, cls.upper(), e["row_id"], rec[c][:40]))
                    else:
                        labels, bad = expand_cols(e["label_cols"], hdr)
                        labels = set(labels)
                        if bad:
                            chk.fail("S9", "%s: entry %s names label column(s) %s that the file does not have" % (rp, e["row_id"], bad))
                        for c, v in rec.items():
                            if c not in labels and v not in (sent, "", "NA"):
                                chk.fail("S9", "%s line %d column %s: not a declared label column of %s and not masked: %r" % (rp, ln, c, e["row_id"], v[:40]))
    for e in ents:
        if e["hits"] == 0:
            chk.fail("S10", "DERIVED_ROWS.tsv %s (%s, %s) matched nothing: an entry that never fires is not a control" % (e["row_id"], e["release_glob"], e["where"]))
    # S11: a column the card definitions call Open Targets-derived must be declared derived (the declaration may not go stale against the shipped definitions)
    cdef = Path(root) / "tables" / "Table2_evidence_cards.columns.tsv"
    if cdef.is_file():
        dh, drows = read_table(cdef, "\t")
        decl_ot = set()
        for e in ents:
            if e["class"] == "ot" and e["status"] == "derived" and fnmatch.fnmatchcase("tables/Table2_evidence_cards.tsv", e["release_glob"]):
                decl_ot |= set(e.get("cols_seen", ()))
        for r in drows:
            rec = dict(zip(dh, r))
            if rec.get("definition", "").startswith("Open Targets-derived") and rec.get("column") not in decl_ot:
                chk.fail("S11", "Table2_evidence_cards.columns.tsv calls column %s Open Targets-derived but no DERIVED_ROWS ot entry masks it" % rec.get("column"))
    return n_rows


# ------------------------------------------------------------------------------------------------------------------ module protocol
def validate_module_result(name, res, out, withheld, before, chk):
    prefix = MODULES[name]
    req = ["release_path", "ship_status", "source_paths", "producer", "licence_class", "rows", "columns", "note"]
    claimed = set()
    for d in res:
        miss = [k for k in req if k not in d]
        if miss:
            chk.fail("M3", "%s returned a record without %s" % (name, miss))
            continue
        rp, ss = norm(d["release_path"]), d["ship_status"]
        if not rp.startswith(prefix) or ".." in rp.split("/"):
            chk.fail("M2", "%s wrote %s outside its namespace %s" % (name, rp, prefix))
        if ss not in SHIPS:
            chk.fail("D1", "%s: ship_status %r" % (name, ss))
        if not licence_ok(d["licence_class"]):
            chk.fail("L1", "%s: %s licence_class %r" % (name, rp, d["licence_class"]))
        root = out if ss == "ship" else withheld
        p = Path(root) / rp
        if not p.is_file():
            chk.fail("M3", "%s: %s not found under %s" % (name, rp, Path(root).name))
            continue
        claimed.add((Path(root).name, rp))
        r_, c_ = shape(p)
        for key, v in (("rows", r_), ("columns", c_)):
            if v != "" and str(d[key]) != v:
                chk.fail("M4", "%s: %s reports %s=%s, the file has %s" % (name, rp, key, d[key], v))
    after = {("out", f) for f in list_files(out)} | {("withheld", f) for f in list_files(withheld)}
    new = after - before
    if new != claimed:
        for x in sorted(new - claimed):
            chk.fail("M3", "%s left an unreported file %s/%s" % (name, x[0], x[1]))
        for x in sorted(claimed - new):
            chk.fail("M3", "%s reported %s/%s that is not new" % (name, x[0], x[1]))


def call_module(ctx, name, ledger, rows, chk, stage2=False, attrs=None):
    sys.path.insert(0, str(ctx.release))
    try:
        sys.modules.pop(name, None)
        importlib.invalidate_caches()
        try:
            mod = importlib.import_module(name)
        except ImportError as e:
            chk.fail("M5", "module %s cannot be imported: %s" % (name, e))
            return []
        for k_, v_ in (attrs or {}).items():
            setattr(mod, k_, v_)
        declared = {norm(r["source_path"]) for r in rows if r["source_path"] and name in r["producer"].split(";")}
        for p in mod.INPUTS:
            if norm(p) not in declared:
                chk.fail("M1", "%s reads %s, which RELEASE_INPUTS.tsv does not declare for it (add a row with producer %s and --pin)" % (name, p, name))
        used = {norm(p) for p in mod.INPUTS}
        for p in sorted(declared - used):
            if not p.endswith("S1R1_CLAIMS_LEDGER.csv"):     # the core adds the ledger to every record's source_paths (F28)
                chk.warn("M1", "%s declares %s but the module does not list it in INPUTS" % (name, p))
        if chk.failures():
            return []
        before = {("out", f) for f in list_files(ctx.out)} | {("withheld", f) for f in list_files(ctx.withheld)}
        if hasattr(mod, "STRICT"):
            mod.STRICT = not stage2      # rel_tables: strict refuses an orphan [HOLD ...]; --stage2 passes every hold through
        try:
            res = mod.build(ctx.v3, ctx.out, ctx.withheld, ledger)
        except Exception as e:      # the modules raise their own check ids (S01.., T..); surface them, not a traceback
            chk.fail("M6", "%s.build raised %s: %s" % (name, type(e).__name__, str(e)[:600]))
            return []
        validate_module_result(name, res, ctx.out, ctx.withheld, before, chk)
        return res
    finally:
        sys.path.pop(0)


# ------------------------------------------------------------------------------------------------------------------ assemble (copy / mask)
def reset_dir(p, base):
    p = Path(p)
    if p.parent != Path(base) or p.name not in ("out", "withheld", "build"):
        raise Refusal("refusing to clear %s" % p)
    if p.exists():
        shutil.rmtree(p)
    p.mkdir(parents=True)


NL_RX = re.compile(r"\s*(?:\r\n|\r|\n)\s*")
NL_NOTE = "declared transform: %d cell(s) held a line break (CR/LF) inside a quoted cell; replaced by ' / ' so that a TSV reader that does not unquote reads one line per row (V10); no other cell changed"


def normalise_newlines(src, dst):
    """Write dst = src with every CR/LF inside a cell replaced by ' / '. Returns the number of cells changed (0 -> dst is not written, the caller copies verbatim)."""
    text = open(src, encoding="utf-8", newline="").read()
    rows = list(csv.reader(io.StringIO(text, newline=""), delimiter="\t"))
    n = 0
    out = []
    for row in rows:
        new = []
        for c in row:
            c2 = NL_RX.sub(" / ", c)
            n += c2 != c
            new.append(c2)
        out.append(new)
    if n:
        with open(dst, "w", encoding="utf-8", newline="") as f:
            csv.writer(f, delimiter="\t", lineterminator="\n").writerows(out)
    return n


def check_tsv_lines(ctx, chk):
    """V10: every shipped .tsv has physical lines == logical rows (no quoted cell spans lines)."""
    n = 0
    for f in sorted(Path(ctx.out).rglob("*.tsv")):
        text = open(f, encoding="utf-8", errors="replace", newline="").read()
        phys = text.count("\n") + (0 if text.endswith("\n") or not text else 1)
        logical = len(list(csv.reader(io.StringIO(text, newline=""), delimiter="\t")))
        n += 1
        if phys != logical:
            chk.fail("S7", "%s: %d physical lines for %d logical rows (a quoted cell holds a line break)" % (norm(f.relative_to(ctx.out)), phys, logical))
    chk.info("S7 tsv lines == rows in %d shipped .tsv files" % n)


def assemble_copies(ctx, rows, chk, fp):
    """Copy / mask every role=copy row. Returns BUILT dict rows and {release_path: clauses} for out/ files."""
    built, rule_of = [], {}
    for r in rows:
        if r["role"] != "copy":
            continue
        src = ctx.resolve(r["source_path"])
        clauses = parse_rule(r["withhold_rule"])
        rp = r["release_path"]
        smd5 = r["source_md5"]
        if r["ship_status"] == "withheld_pending_terms":
            dst = ctx.withheld / rp
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            built.append(dict(release_path=rp, md5=md5_file(dst), source_path=r["source_path"], source_md5=smd5, producer="P31R_build_release.py:copy",
                              licence_class=r["licence_class"], ship_status="withheld_pending_terms", input_id=r["input_id"], store="withheld",
                              note="whole file withheld until the provider's terms are read" + ("; " + r["note"] if r["note"] else "")))
            continue
        dst = ctx.out / rp
        dst.parent.mkdir(parents=True, exist_ok=True)
        mc = masking_clauses(clauses)
        if not mc:
            nl = normalise_newlines(src, dst) if rp.endswith(".tsv") else 0
            if not nl:
                shutil.copyfile(src, dst)
            built.append(dict(release_path=rp, md5=md5_file(dst), source_path=r["source_path"], source_md5=smd5,
                              producer="P31R_build_release.py:copy" if not nl else "P31R_build_release.py:normalise_newlines",
                              licence_class=r["licence_class"], ship_status="ship", input_id=r["input_id"], store="out",
                              note=(r["note"] + "; " if r["note"] and nl else r["note"]) + (NL_NOTE % nl if nl else "")))
            rule_of[rp] = clauses
            continue
        d = delim_of(rp)
        if d is None:
            chk.fail("S6", "%s: a masking rule needs a delimited file" % rp)
            continue
        hdr, body = read_table(src, d)
        body = [list(x) + [""] * (len(hdr) - len(x)) for x in body]
        # round trip: rewriting without masking must reproduce the parsed cells (guards the csv dialect)
        probe = io.StringIO(newline="")
        csv.writer(probe, delimiter=d, lineterminator="\n").writerows([hdr] + body)
        back = list(csv.reader(io.StringIO(probe.getvalue(), newline=""), delimiter=d))
        if back != [hdr] + body:
            chk.fail("S6", "%s: csv round trip changes cells; refusing to rewrite it" % rp)
            continue
        masked = apply_rules(ctx, hdr, body, clauses, rp, chk)
        for cls, col, orig in masked:
            fp.add(cls, orig, (rp, col))
        require_check(hdr, body, clauses, rp, chk)
        nl = 0
        if d == "\t":          # V10: declared transform on the shipped (masked) copy only; the withheld counterpart stays byte-identical to the source
            for row in body:
                for i, c in enumerate(row):
                    c2 = NL_RX.sub(" / ", c)
                    nl += c2 != c
                    row[i] = c2
        write_table(dst, hdr, body, d)
        wh = ctx.withheld / rp
        wh.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, wh)
        classes = sorted({c for c, _, _ in masked})
        h9only = classes == ["h9"]
        status = SHIP_H9 if h9only else "+".join(x for x in ("withheld_pending_terms", SHIP_H9 if "h9" in classes else "") if x)
        why = ("held for source hold H9 (the identity and N of %s are unconfirmed; P31R_ACTIVE_FIGURES.json holds[H9]); it does not wait for a provider's terms" % h9_spec(ctx)["source_value"]) if h9only else               ("held back until the provider's terms are read" + ("; plus source hold H9 on its %s magnitudes" % h9_spec(ctx)["source_value"] if "h9" in classes else ""))
        built.append(dict(release_path=rp, md5=md5_file(dst), source_path=r["source_path"], source_md5=smd5,
                          producer="P31R_build_release.py:mask(%s)" % ",".join(classes), licence_class=r["licence_class"], ship_status="ship",
                          input_id=r["input_id"], store="out",
                          note="%d cell(s) masked by the declared rule; the full file is in withheld/%s%s%s" % (len(masked), rp, ("; " + r["note"]) if r["note"] else "", ("; " + NL_NOTE % nl) if nl else "")))
        built.append(dict(release_path=rp, md5=md5_file(wh), source_path=r["source_path"], source_md5=smd5, producer="P31R_build_release.py:copy",
                          licence_class=r["licence_class"], ship_status=status, input_id=r["input_id"], store="withheld",
                          note="unmasked counterpart of out/%s; %s" % (rp, why)))
        rule_of[rp] = clauses
    return built, rule_of


def module_rule_clauses(rows):
    out = {}
    for r in rows:
        if r["role"] == "module_output" and r["ship_status"] == "ship":   # the scan reads out/ only
            out[r["release_path"]] = parse_rule(r["withhold_rule"])
    return out


# ------------------------------------------------------------------------------------------------------------------ claims
def check_claims(ctx, root_out, root_wh, built, chk):
    """Every file in out/ and withheld/ has a BUILT row, every BUILT row has its file."""
    have = {("out", f) for f in list_files(root_out)} | {("withheld", f) for f in list_files(root_wh)}
    want = {(b["store"], b["release_path"]) for b in built}
    for s, f in sorted(have - want):
        chk.fail("C1", "unclaimed file %s/%s (no BUILT_FILES row)" % (s, f))
    for s, f in sorted(want - have):
        chk.fail("C2", "BUILT_FILES row without a file: %s/%s" % (s, f))


# ------------------------------------------------------------------------------------------------------------------ not-shipped derivation
def manifest_rows(ctx):
    p = ctx.r2 / "figures" / "registry" / "APPROVED_SOURCE_MANIFEST.tsv"
    hdr, rows = read_dicts(p)
    for i, r in enumerate(rows, 2):
        r["_line"] = i
    return rows


def file_md5_or_pin(ctx, path, pinned=None):
    p = Path(path)
    if not p.is_file():
        return pinned, None, "absent"
    sz = p.stat().st_size
    if sz > HASH_CAP:
        return (pinned, sz, "manifest_pin") if pinned else (None, sz, "unhashed_over_cap")
    return ctx.md5(p), sz, "computed"


def prune_dirs(ctx, root, dirs):
    """Directories a walk over V3 must not enter: caches, and the release builder's OWN outputs (out/, withheld/, build/ of this release and of the
    folder the build is writing to, and the canonical R2_rewrite/release/ ones, which a scratch build or the selftest must not read either). A shipped copy that was edited by a withhold rule (copy B with the H9 cells blanked) has another md5 than its
    source and would otherwise be listed as a non-canonical copy of itself, so a second --derive-not-shipped after a build would refuse the next build (N1)."""
    own = {b / d for b in (ctx.release, ctx.r2 / "release") for d in ("out", "withheld", "build")} | {ctx.out, ctx.withheld, ctx.build}
    return [d for d in dirs if d not in (".git", "__pycache__") and (Path(root) / d).resolve() not in own]


GENO_EXT = (".bed", ".bim", ".fam", ".raw", ".vcf", ".vcf.gz")
WINDOW_DIRS = {"extracts": "THIRD_PARTY_WINDOW", "ensembl_rest": "THIRD_PARTY_WINDOW", "plink": "GENOTYPE_DERIVED", "ref": "GENOTYPE_DERIVED", "ld": "GENOTYPE_DERIVED"}
PRIORART_DL = "R1_stage1/S1R1_R2_09_priorart_downloads_20261010/"


def check_canonical(ctx, ri, chk):
    """F11: the canonical copies are pinned in inputs/CANONICAL_SOURCES.tsv, which --pin never rewrites. A RELEASE_INPUTS row of that file name must point at that path
    with that md5, and the file on disk must have it; otherwise a re-pin could silently ship the wrong copy. Returns {file_name: md5}."""
    canon = load_declared(ctx, "CANONICAL_SOURCES.tsv")
    out = {}
    for c in canon:
        out[c["file_name"]] = c["md5"]
        p = ctx.resolve(c["source_path"])
        if not p.is_file():
            chk.fail("G1", "canonical %s: %s is absent" % (c["canon_id"], c["source_path"]))
        elif ctx.md5(p) != c["md5"]:
            chk.fail("G1", "canonical %s: %s has md5 %s, CANONICAL_SOURCES.tsv declares %s" % (c["canon_id"], c["source_path"], ctx.md5(p)[:8], c["md5"][:8]))
        rows = [r for r in ri if r["role"] == "copy" and Path(r["source_path"]).name == c["file_name"]]
        if len(rows) != 1:
            chk.fail("G1", "canonical %s: %d RELEASE_INPUTS copy row(s) carry the file name %s (exactly one expected)" % (c["canon_id"], len(rows), c["file_name"]))
        for r in rows:
            if norm(r["source_path"]) != norm(c["source_path"]) or r["source_md5"] != c["md5"]:
                chk.fail("G1", "%s ships %s (md5 %s), not the canonical copy %s (md5 %s): refused even after --pin" %
                         (r["input_id"], r["source_path"], (r["source_md5"] or "-")[:8], c["source_path"], c["md5"][:8]))
    return out


def derive_not_shipped(ctx, ri, rules, chk):
    """Rows (path, md5, bytes, reason_code, evidence, md5_basis), path-and-md5 only, never a folder rule. Returns (rows, overrides actually used).
    Overrides are HAND-DECLARED (inputs/NOT_SHIPPED_OVERRIDES.tsv); the derivation refuses any it would have to apply on its own (N4)."""
    out = {}
    decl_ov = {(o["path"], o["md5"]): o for o in load_declared(ctx, "NOT_SHIPPED_OVERRIDES.tsv")}
    used_ov = {}
    ship_src = {norm(r["source_path"]) for r in ri if r["role"] == "copy" and r["ship_status"] in SHIPS}
    ship_md5 = {r["source_md5"] for r in ri if r["role"] == "copy"}
    ri_paths = {norm(r["source_path"]) for r in ri if r["source_path"]}
    decl_md5 = {r["source_md5"] for r in ri if r["source_path"] and r["source_md5"]}

    def override(path, md5, line, why):
        key = (path, md5)
        if key in decl_ov:
            used_ov[key] = decl_ov[key]
        else:
            chk.fail("N4", "%s (md5 %s): %s -- this needs an override, and overrides are hand-declared: add it to NOT_SHIPPED_OVERRIDES.tsv after reading why, or fix the "
                     "declaration" % (path, (md5 or "-")[:8], why))

    def put(path, reason, evidence, pinned=None):
        path = norm(path)
        if path in ri_paths:
            return      # a declared input or deliverable is never listed as unshipped
        p = ctx.resolve(path)
        md5, sz, basis = file_md5_or_pin(ctx, p, pinned)
        if basis == "absent":
            chk.info("not-shipped: %s is absent on this machine (%s)" % (path, reason))
            out.setdefault(path, dict(path=path, md5=pinned or "", bytes="", reason_code=reason, evidence=evidence + "; ABSENT on this machine",
                                      md5_basis="manifest_pin" if pinned else "none"))
            return
        if md5 and md5 in decl_md5:
            override(path, md5, "%s, but byte-identical to a declared input or deliverable of the release" % reason, None)
            return
        if path in out:
            if reason not in out[path]["reason_code"].split("+"):
                out[path]["reason_code"] += "+" + reason
            return
        out[path] = dict(path=path, md5=md5 or "", bytes=str(sz), reason_code=reason, evidence=evidence, md5_basis=basis)

    man = manifest_rows(ctx)
    canon_md5 = check_canonical(ctx, ri, chk)
    # (a) manifest rows marked SUPERSEDED. Not listed when the file is a declared source: that is an override, and only a HAND-DECLARED one is accepted.
    for m in man:
        if m["disposition"] == "SUPERSEDED":
            sp = norm(m["source_file"])
            ev = "APPROVED_SOURCE_MANIFEST.tsv:%d figure=%s panel=%s" % (m["_line"], m["figure"], m["panel"])
            if sp in ship_src or m["md5"] in ship_md5:
                override(sp, m["md5"], "SUPERSEDED manifest row %d (%s %s) names a declared deliverable" % (m["_line"], m["figure"], m["panel"]), None)
            else:
                put(sp, "SUPERSEDED_MANIFEST", ev, m["md5"])
    # (b) spec 7(b), literally: consumer set = the manifest's non-SUPERSEDED source_file values + the RELEASE_INPUTS sources; every V4 / V4_local figure-folder
    # file outside it is listed (a raw third-party window or per-SNP coloc input is a THIRD_PARTY_WINDOW).
    man_src = {norm(m["source_file"]) for m in man if m["disposition"] != "SUPERSEDED"}
    consumer = man_src | ri_paths
    tot = n_cons = n_win = n_other = 0
    for base in ("V4/source_data/figures", "V4_local/source_data/figures"):
        root = ctx.v3 / base
        for f in sorted(root.rglob("*")):
            if not f.is_file():
                continue
            tot += 1
            rel = ctx.rel(f)
            if rel in consumer:
                n_cons += 1
                continue
            win = re.search(r"/(?:slices|coloc_inputs)/", rel) is not None
            before = len(out)
            put(rel, "THIRD_PARTY_WINDOW" if win else "V4_FRAME_NOT_A_SOURCE",
                ("per-window third-party summary statistics or per-SNP coloc inputs (plan l.157-160)" if win else
                 "outside the consumer set (non-SUPERSEDED manifest source_file values + RELEASE_INPUTS sources), spec section 7(b)"))
            if len(out) > before:
                n_win += win
                n_other += (not win)
    chk.info("not-shipped (b): V4/V4_local figure folders hold %d file(s); %d are in the consumer set (%d manifest sources + %d RELEASE_INPUTS sources, union); "
             "%d are listed (%d THIRD_PARTY_WINDOW, %d V4_FRAME_NOT_A_SOURCE); %d remaining are byte-identical to a declared source (declared overrides)" %
             (tot, n_cons, len(man_src), len(ri_paths), n_win + n_other, n_win, n_other, tot - n_cons - n_win - n_other))
    # (c) non-canonical copies, found by content md5, against the PINNED canonical md5 (CANONICAL_SOURCES.tsv), never against RELEASE_INPUTS itself
    for name, cm in sorted(canon_md5.items()):
        for root, dirs, files in os.walk(ctx.v3):
            dirs[:] = prune_dirs(ctx, root, dirs)
            if name in files:
                p = Path(root) / name
                md5 = ctx.md5(p)
                if md5 != cm:
                    put(ctx.rel(p), "NONCANONICAL_COPY", "same file name as the canonical %s (md5 %s, CANONICAL_SOURCES.tsv) and different content" % (name, cm[:8]))
    # (d) decoy-null inputs named in the decoy log (md5 from the log; the on-disk md5 must agree)
    log = ctx.v3 / "R1_stage1" / "S1R1_03_decoy_method" / "S1R1_decoy_rates.log"
    n_log = 0
    for ln, line in enumerate(log.read_text(encoding="utf-8").splitlines(), 1):
        m = re.match(r"^md5 (\S+)\s+([0-9a-f]{32})\s+(\S+)$", line)
        if m:
            n_log += 1
            path, logged = m.group(3), m.group(2)
            put(path, "DECOY_NULL_INPUT", "S1R1_decoy_rates.log:%d (%s, logged md5 %s)" % (ln, m.group(1), logged[:8]), logged)
            if path in out and out[path]["md5"] and out[path]["md5"] != logged:
                chk.warn("N2", "decoy input %s: on-disk md5 %s differs from the logged %s" % (path, out[path]["md5"][:8], logged[:8]))
    if n_log == 0:
        chk.fail("N2", "decoy log lists no md5 line: the six decoy inputs cannot be derived")
    # (e) one walk of V3: individual-level dosages and genotype files; the workstream window extracts and genotype-derived folders (F24); the prior-art downloads
    for root, dirs, files in os.walk(ctx.v3):
        dirs[:] = prune_dirs(ctx, root, dirs)
        for fn in files:
            p = Path(root) / fn
            rel = ctx.rel(p)
            parts = rel.split("/")
            if fn.endswith(".i8"):
                put(rel, "I8_DOSAGE", "individual-level 1000 Genomes dosage matrix (plan l.157-160)")
            elif fn.endswith(".i2"):
                put(rel, "I2_BOOT_INDEX", "bootstrap resample index over the reference-panel members")
            elif fn.endswith(GENO_EXT):
                put(rel, "GENOTYPE_DERIVED", "genotype-format file (%s): individual-level or reference-panel genotypes, never released (plan l.157)" % Path(fn).suffix)
            elif rel.startswith(PRIORART_DL):
                put(rel, "PRIORART_DOWNLOAD", "prior-art download folder (publisher / repository files fetched for the literature search): third-party text, never released")
            elif rel.startswith("R1_stage1/S1R1_R2_") and any(d in WINDOW_DIRS for d in parts[:-1]):
                d = [x for x in parts[:-1] if x in WINDOW_DIRS][0]
                put(rel, WINDOW_DIRS[d], "workstream %s/ folder: %s (plan l.157-160, spec section 3)" %
                    (d, "cis-window extract of third-party summary statistics" if WINDOW_DIRS[d] == "THIRD_PARTY_WINDOW" else "genotype- or reference-panel-derived file"))
            elif rel.startswith("R1_stage1/S1R1_R2_") and fn.startswith("extract_"):
                put(rel, "THIRD_PARTY_WINDOW", "window extract written by a workstream (file name extract_*): a cis-window of third-party summary statistics (the TPMI extract included)")
    mpins = {norm(m["source_file"]): m["md5"] for m in man}
    for rr in rules:
        if rr["list_not_shipped"] != "yes":
            continue
        pat = rr["pattern"]
        if any(ch in pat for ch in "*?["):
            base = pat.split("*")[0].split("?")[0].split("[")[0]
            root = ctx.resolve(base.rsplit("/", 1)[0]) if "/" in base else ctx.v3
            for f in sorted(root.rglob("*")) if root.is_dir() else []:
                if f.is_file() and fnmatch.fnmatchcase(ctx.rel(f), pat):
                    put(ctx.rel(f), "THIRD_PARTY_WINDOW", "SOURCE_CLASS_RULES %s (%s)" % (rr["rule_id"], rr["resource"]), mpins.get(ctx.rel(f)))
        else:
            put(pat, "THIRD_PARTY_WINDOW", "SOURCE_CLASS_RULES %s (%s)" % (rr["rule_id"], rr["resource"]), mpins.get(pat))
    # (f) AF3 plotted frames no registry claims: only the HAND-DECLARED orphans (ORPHAN_FRAMES.tsv, path + md5) are accepted; anything else is F2
    claimed = set()
    for unit, _ in active_units(ctx):
        reg = json.loads((ctx.r2 / "figures" / "registry" / (unit + ".json")).read_text(encoding="utf-8"))
        for fr in reg.get("frames") or []:
            claimed.add("%s__%s.csv" % (unit, fr))
    orphans = {o["path"]: o for o in load_declared(ctx, "ORPHAN_FRAMES.tsv")}
    sd = ctx.r2 / "figures" / "af3" / "source_data"
    seen_orph = set()
    for f in sorted(sd.glob("*")):
        if f.is_file() and f.name not in claimed:
            rel = ctx.rel(f)
            o = orphans.get(rel)
            if o is None or ctx.md5(f) != o["md5"]:
                chk.fail("F2", "af3/source_data/%s is claimed by no registry and is not on the hand-declared ORPHAN_FRAMES.tsv with this md5 (md5 %s)" % (f.name, ctx.md5(f)[:8]))
            else:
                seen_orph.add(rel)
                put(rel, "STALE_ORPHAN_FRAME", "named in no registry frames list; declared in ORPHAN_FRAMES.tsv (%s)" % o["reason"])
    for rel in sorted(set(orphans) - seen_orph):
        chk.warn("F5", "ORPHAN_FRAMES.tsv declares %s, which is no longer an unclaimed frame (stale declaration)" % rel)
    # (g) priorart files that carry verbatim publication text (header scan)
    pa = ctx.v3 / "R1_stage1" / "S1R1_R2_09_priorart"
    for f in sorted(pa.glob("*.tsv")):
        hdr = read_table(f, "\t")[0]
        hit = sorted(VERBATIM_COLS & set(hdr))
        if hit:
            put(ctx.rel(f), "VERBATIM_PUBLICATION_TEXT", "column(s) %s hold passages extracted from publications" % hit)
    for f in sorted(pa.glob("*.txt")) + sorted(pa.glob("findings*.md")):
        put(ctx.rel(f), "VERBATIM_PUBLICATION_TEXT", "prior-art text extract / findings file quoting publications")
    # (h) files behind a number an operative surface prints that are neither shipped nor withheld: a hand-declared disposition with its reason (F09)
    for dd in load_declared(ctx, "NOT_SHIPPED_DECLARED.tsv"):
        put(dd["path"], dd["reason_code"], dd["evidence"])
    for k in sorted(set(decl_ov) - set(used_ov)):
        chk.warn("N5", "NOT_SHIPPED_OVERRIDES.tsv declares an override for %s (md5 %s) that nothing needs any more" % (k[0], k[1][:8]))
    rows = sorted(out.values(), key=lambda r: r["path"])
    return rows, sorted(used_ov.values(), key=lambda r: r["path"])


def load_not_shipped(ctx):
    p = ctx.inputs / "NOT_SHIPPED.tsv"
    if not p.is_file():
        return None
    hdr, rows = read_dicts(p)
    if hdr != NS_COLS:
        raise Refusal("NOT_SHIPPED.tsv header differs: %s" % hdr)
    return rows


def check_not_shipped(ctx, ns, ri, root_out, root_wh, chk):
    """N1: no listed md5 in out/ or withheld/.  N3: no declared deliverable has a listed md5."""
    bad = {r["md5"]: r for r in ns if r["md5"]}
    for root in (root_out, root_wh):
        for f in list_files(root):
            m = md5_file(Path(root) / f)
            if m in bad:
                chk.fail("N1", "%s/%s has md5 %s of a NOT_SHIPPED file (%s: %s)" % (Path(root).name, f, m[:8], bad[m]["path"], bad[m]["reason_code"]))
    for r in ri:
        if r["role"] == "copy" and r["source_md5"] in bad:
            chk.fail("N3", "%s: deliverable source %s has the md5 of NOT_SHIPPED %s" % (r["input_id"], r["source_path"], bad[r["source_md5"]]["path"]))



LEDGER_REL = "R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv"
LEDGER_CALL_RX = re.compile(r"\bp31r_ledger(?:_row)?\s*\(")
LEDGER_DEF_RX = re.compile(r"\bp31r_ledger(?:_row)?\s*<-\s*function")
SOURCE_RX = re.compile(r"\bsource\(\s*[\"']([^\"']+)[\"']")
DIRSRC_RX = re.compile(r"file\.path\(dirname\(\.(?:me|this)\),\s*[\"']([^\"']+)[\"']\)")


def generator_reads(ctx, unit):
    """V01: every input a unit's generator reads that the approved source manifest does not list: the claims ledger (any p31r_ledger() / p31r_ledger_row() call
    in the generator or in a file it source()s, the file that DEFINES p31r_ledger excepted) and each source()d file named by a literal path.
    Returns [(path relative to V3, why)]."""
    d = json.loads((ctx.r2 / "figures" / "af3" / "P31R_ACTIVE_FIGURES.json").read_text(encoding="utf-8"))
    base = re.sub(r"_p\d+$", "", unit)
    gen = next((f["generator"] for f in d["figures"] if f["id"] == base), None)
    if not gen:
        return []
    res, seen, todo, ledger_by = [], set(), [gen], None
    while todo:
        rel = todo.pop(0)
        if rel in seen:
            continue
        seen.add(rel)
        p = ctx.resolve(rel)
        if not p.is_file():
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        code = "\n".join(l for l in txt.splitlines() if not l.lstrip().startswith("#"))
        if not LEDGER_DEF_RX.search(code) and LEDGER_CALL_RX.search(code) and ledger_by is None:
            ledger_by = rel
        for m in SOURCE_RX.finditer(code):
            q = norm(m.group(1))
            if "CKM papers/V3/" in q:
                q = q.split("CKM papers/V3/", 1)[1]
            elif Path(q).is_absolute():
                continue
            else:
                q = norm(Path(rel).parent / q)
            todo.append(q)
        for m in DIRSRC_RX.finditer(code):
            todo.append(norm(Path(rel).parent / m.group(1)))
    if ledger_by:
        res.append((LEDGER_REL, "read by %s (p31r_ledger() / num(...)) through P31R_fig_common.R; not released" % ledger_by))
    shipped = {r["path"] for r in read_dicts(ctx.inputs / "CODE_INCLUDE.tsv")[1] if r["ship"] == "1"}      # source()d code that the bundle ships counts as released
    res += [(r, "source()d by the generator chain") for r in sorted(seen) if r != gen and r not in shipped and ctx.resolve(r).is_file()]
    return res


LEDGER_PUBLIC_REL = "ledger/claims_ledger_public.tsv"
LEDGER_ID_RX = re.compile(r"[\"'](L\d{3,4})[\"']")
LEDGER_NONLIT_RX = re.compile(r"\b(?:p31r_ledger_row|p31r_ledger|num)\(\s*(?![\"'\s]|id\b)")      # a wrapper passing its own parameter `id` through (num <- function(id, idx) p31r_ledger(id, idx)) is not a consumer


def unit_ledger_ids(ctx, unit):
    """(ledger ids a unit's generator chain names as string literals, lines that call the ledger with a non-literal id). Same chain walk as generator_reads."""
    d = json.loads((ctx.r2 / "figures" / "af3" / "P31R_ACTIVE_FIGURES.json").read_text(encoding="utf-8"))
    base = re.sub(r"_p\d+$", "", unit)
    gen = next((f["generator"] for f in d["figures"] if f["id"] == base), None)
    ids, dyn, seen, todo = set(), [], set(), [gen] if gen else []
    while todo:
        rel = todo.pop(0)
        if rel in seen:
            continue
        seen.add(rel)
        p = ctx.resolve(rel)
        if not p.is_file():
            continue
        code = "\n".join(l for l in p.read_text(encoding="utf-8", errors="replace").splitlines() if not l.lstrip().startswith("#"))
        if not LEDGER_DEF_RX.search(code):       # the file that DEFINES p31r_ledger calls its own wrappers with the parameter; it is not a consumer
            ids |= set(LEDGER_ID_RX.findall(code))
            dyn += [l.strip() for l in code.splitlines() if LEDGER_NONLIT_RX.search(l) and "<- function" not in l]
        for m in SOURCE_RX.finditer(code):
            q = norm(m.group(1))
            if "CKM papers/V3/" in q:
                q = q.split("CKM papers/V3/", 1)[1]
            elif Path(q).is_absolute():
                continue
            else:
                q = norm(Path(rel).parent / q)
            todo.append(q)
        for m in DIRSRC_RX.finditer(code):
            todo.append(norm(Path(rel).parent / m.group(1)))
    return ids, dyn


def public_ledger_masks(ctx):
    """{ledger id: reason} of the rows the public extract masks; None when the extract is not in out/."""
    p = ctx.out / LEDGER_PUBLIC_REL
    if not p.is_file():
        return None
    hdr, rows = read_table(p)
    iv, ii = hdr.index("value"), hdr.index("id")
    res = {}
    for r in rows:
        if r[iv].startswith("{"):
            try:
                j = json.loads(r[iv])
            except ValueError:
                continue
            if j.get("masked"):
                res[r[ii]] = j.get("reason", "")
    return res


# ------------------------------------------------------------------------------------------------------------------ regenerability
def level2_map(ctx):
    """RC1b (S2-D34(3)): {source path: LEVEL2_INPUTS row} of the sources a figure generator reads that the archive answers with a declared masked or reduced copy
    (inputs/LEVEL2_INPUTS.tsv; rel_portability writes the extracts and marks the map rows lift = 1)."""
    p = Path(ctx.inputs) / "LEVEL2_INPUTS.tsv"
    if not p.is_file():
        return {}
    hdr, rows = read_table(p)
    return {norm(r[hdr.index("source_path")]): dict(zip(hdr, r)) for r in rows}


def regenerability(ctx, ri, built, rules, chk):
    """Per unit (13 render units) and per non-figure deliverable. Returns (by_unit rows, by_source rows).
    class = the conservative reading: a source counts as released only if a released file has its md5 (verbatim copy).
    class_spec_literal = the plan's wording: also released when the md5 is the SOURCE md5 of a released derived file. That reading counts a ledger
    that a table merely drew a few values from as 'released', so it is reported beside, never instead."""
    man = [m for m in manifest_rows(ctx)]
    out_md5 = {b["md5"] for b in built if b["store"] == "out"}
    src_md5 = set()
    for b in built:
        if b["store"] == "out":
            src_md5.update(x for x in (b["source_md5"] or "").split(";") if x)
    by_source, by_unit = [], []
    deliv_of = {r["release_path"]: set((r["deliverable"] or "").split(";")) - {""} for r in ri if r["release_path"]}

    def masked_of(unit_or_group):
        """Released files of the unit / deliverable group whose out/ copy was edited by a withhold rule (so it is NOT the source file): they are marked masked."""
        res = set()
        for b in built:
            if b["store"] == "out" and "mask(" in b["producer"]:
                if unit_of_frame(b["release_path"]) == unit_or_group or unit_or_group in deliv_of.get(b["release_path"], ()):
                    res.add(b["release_path"])
        return res

    def classify(path, pinned):
        p = ctx.resolve(path)
        md5, sz, basis = file_md5_or_pin(ctx, p, pinned)
        if basis == "absent":
            basis = "manifest_pin(absent)"
        if md5 and md5 in out_md5:
            return dict(strict="regenerates_from_release", spec="regenerates_from_release", cov="verbatim", why="a released file has this md5", md5=md5, basis=basis, nrb="")
        lv = l2.get(norm(path))
        if lv is not None and lv["release_path"] in shipped_paths:
            return dict(strict="regenerates_from_release", spec="regenerates_from_release", cov="level2_lift",
                        why="%s: released as %s (%s); the clean-room regeneration compares the plotted frames with the live frames" % (lv["l2_id"], lv["release_path"], lv["kind"]),
                        md5=md5, basis=basis, nrb="")
        r = classify_source(rules, path)
        if r is None:
            chk.fail("R1", "source %s matches no SOURCE_CLASS_RULES row (no default)" % path)
            return dict(strict="UNCLASSIFIED", spec="UNCLASSIFIED", cov="", why="", md5=md5, basis=basis, nrb="")
        derived = bool(md5 and md5 in src_md5)
        return dict(strict=r["class"], spec="regenerates_from_release" if derived else r["class"], cov="derived" if derived else "",
                    why="%s: %s" % (r["rule_id"], r["reason"]), md5=md5, basis=basis, nrb=r["basis"])

    masks = public_ledger_masks(ctx)
    l2 = level2_map(ctx)
    shipped_paths = {b["release_path"] for b in built if b["store"] == "out"}

    def ledger_class(unit, c):
        """RC1: the claims ledger a generator reads is released as the public extract. The unit stands on it only if every ledger id its chain names is released
        unmasked and the chain makes no non-literal ledger call (a masked or unseen id stops the shipped generator at run time, and is then not claimed)."""
        if masks is None:
            return c
        ids, dyn = unit_ledger_ids(ctx, unit)
        okdyn = {(d["unit"], d["line_text"]) for d in load_declared(ctx, "LEDGER_DYNAMIC_DECLARED.tsv")}      # RC1b: a non-literal call whose ids come from literal vectors of the same chain, declared by hand
        dyn = [l for l in dyn if (unit, l) not in okdyn]
        hit = sorted(i for i in ids if i in masks)
        if hit:
            reasons = sorted({masks[i] for i in hit})
            return dict(c, strict="needs_nonredistributable_input", spec="needs_nonredistributable_input", cov="",
                        why="public extract masks ledger row(s) %s (%s) that the generator reads" % (",".join(hit), ",".join(reasons)), nrb="third_party_terms")
        if dyn:
            return dict(c, strict="needs_nonredistributable_input", spec="needs_nonredistributable_input", cov="",
                        why="the generator chain calls the ledger with a non-literal id (%d line(s)); its ids cannot be checked against the masked rows" % len(dyn), nrb="not_on_approved_list")
        return dict(c, strict="regenerates_from_release", spec="regenerates_from_release", cov="public_extract",
                    why="released as %s; none of the %d ledger id(s) the generator names is masked" % (LEDGER_PUBLIC_REL, len(ids)), nrb="")

    def worst(classes):
        return max([c for c in classes if c in CLASS_ORDER] or ["regenerates_from_release"], key=CLASS_ORDER.index)

    def add_source(unit, s, kind, disp, c):
        by_source.append(dict(unit=unit, source_path=s, kind=kind, manifest_disposition=disp, class_=c["strict"], class_spec_literal=c["spec"],
                              coverage=c["cov"], md5=c["md5"] or "", md5_basis=c["basis"], reason=c["why"], not_released_because=c["nrb"]))

    for unit, mid in active_units(ctx):
        seen = {}
        for m in man:
            if m["figure"] == mid and m["disposition"] != "SUPERSEDED":
                seen.setdefault(norm(m["source_file"]), (m["md5"], m["disposition"], m["kind"]))
        if not seen:
            chk.fail("R2", "unit %s (manifest id %s) has no manifest source rows" % (unit, mid))
        for gp, why in generator_reads(ctx, unit):       # V01: inputs the generator reads that the manifest rows do not list
            if norm(gp) not in seen:
                seen[norm(gp)] = (None, "", "generator read: " + why)
        cs, cl = [], []
        blocking, blocked = [], []
        cnt = collections.Counter()
        for s, (md5, disp, kind) in sorted(seen.items()):
            c = classify(s, md5)
            if norm(s) == LEDGER_REL:       # RC1b: a manifest row for the ledger names the same file the generator reads through the public extract
                c = ledger_class(unit, c)
            if kind.startswith("generator read"):
                c["spec"] = c["strict"]       # a generator's own read is never 'covered' by a derived file that merely drew values from the input
            add_source(unit, s, kind, disp, c)
            if disp == "BLOCKED":
                blocked.append(s)
                continue
            cs.append(c["strict"]); cl.append(c["spec"]); cnt[c["strict"]] += 1
            if c["strict"] != "regenerates_from_release":
                blocking.append(s)
        by_unit.append(dict(unit=unit, kind="figure", manifest_id=mid, n_sources=len(cs), released=cnt["regenerates_from_release"],
                            needs_public_download=cnt["needs_public_download"], does_not_regenerate_engine_output=cnt["does_not_regenerate_engine_output"],
                            needs_nonredistributable_input=cnt["needs_nonredistributable_input"],
                            covered_only_by_derived_file=sum(1 for x in by_source if x["unit"] == unit and x["coverage"] == "derived"),
                            blocked_not_counted=";".join(blocked), class_=worst(cs), class_spec_literal=worst(cl), blocking_sources=";".join(blocking),
                            masked_files=";".join(sorted(masked_of(unit)))))
    # non-figure deliverables: copy rows are released verbatim; the inputs of module outputs are classified like figure sources
    groups = collections.defaultdict(lambda: dict(copies=0, inputs={}))
    for r in ri:
        for g in (r["deliverable"] or "").split(";"):
            if not g:
                continue
            if r["role"] == "copy" and r["origin"] == "curated":
                groups[g]["copies"] += 1
            elif r["role"] == "input" and r["source_path"]:
                groups[g]["inputs"][r["source_path"]] = r["source_md5"]
    for g in sorted(groups):
        cs, cl, cnt, blocking = [], [], collections.Counter(), []
        for s, md5 in sorted(groups[g]["inputs"].items()):
            c = classify(s, md5)
            add_source(g, s, "deliverable input", "", c)
            cs.append(c["strict"]); cl.append(c["spec"]); cnt[c["strict"]] += 1
            if c["strict"] != "regenerates_from_release":
                blocking.append(s)
        if not groups[g]["inputs"]:
            # result tables of the analysis engines / workstreams, shipped as stored (plan l.127, item 2): nothing regenerates them from this release.
            # A group with masked files says so in masked_files; it is never called "verbatim".
            cls = spec = "does_not_regenerate_engine_output"
        else:
            cls, spec = worst(cs), worst(cl)
        by_unit.append(dict(unit=g, kind="deliverable", manifest_id="", n_sources=sum(cnt.values()) + groups[g]["copies"],
                            released=cnt["regenerates_from_release"] + groups[g]["copies"], needs_public_download=cnt["needs_public_download"],
                            does_not_regenerate_engine_output=cnt["does_not_regenerate_engine_output"], needs_nonredistributable_input=cnt["needs_nonredistributable_input"],
                            covered_only_by_derived_file=sum(1 for x in by_source if x["unit"] == g and x["coverage"] == "derived"),
                            blocked_not_counted="", class_=cls, class_spec_literal=spec, blocking_sources=";".join(blocking),
                            masked_files=";".join(sorted(masked_of(g)))))
    return by_unit, by_source


UNIT_COLS = ["unit", "kind", "manifest_id", "n_sources", "released", "needs_public_download", "needs_nonredistributable_input",
             "does_not_regenerate_engine_output", "covered_only_by_derived_file", "blocked_not_counted", "class", "class_spec_literal", "blocking_sources", "masked_files"]
SRC_COLS = ["unit", "source_path", "kind", "manifest_disposition", "class", "class_spec_literal", "coverage", "md5", "md5_basis", "reason", "not_released_because"]


def _unit_row(d):
    return [d["class_"] if c == "class" else d[c] for c in UNIT_COLS]


def _src_row(d):
    return [d["class_"] if c == "class" else d[c] for c in SRC_COLS]


# ------------------------------------------------------------------------------------------------------------------ build
def fill_row(ctx, b, texts, outkeys):
    """bytes / rows / columns / cited_by of a built row added after the first fill (code bundle, release documents)."""
    k = outkeys.get((b["release_path"], b["ship_status"])) or (Path(b["release_path"]).name if not b["input_id"] else "")
    p = Path(ctx.out if b["store"] == "out" else ctx.withheld, b["release_path"])
    b["cited_by"] = cited_by(texts, k)
    b["bytes"] = str(p.stat().st_size)
    b["rows"], b["columns"] = shape(p)


def classify_built(ctx, built, chk):
    """licence_class of every built row = rel_licence.classify(file, row): the providers whose data the file's VALUES derive from (read from its own columns,
    else from its declaration). Module outputs start from the class the module reports."""
    sys.path.insert(0, str(ctx.release))
    try:
        sys.modules.pop("rel_licence", None)
        importlib.invalidate_caches()
        try:
            rl = importlib.import_module("rel_licence")
        except ImportError as e:
            chk.fail("M5", "module rel_licence cannot be imported: %s" % e)
            return
    finally:
        sys.path.pop(0)
    rl.INPUTS = ctx.inputs
    rl._CACHE.clear()
    for b in built:
        f = Path(ctx.out if b["store"] == "out" else ctx.withheld) / b["release_path"]
        try:
            new, ev = rl.classify(f, b)
        except ValueError as e:
            chk.fail("L4", "%s: %s" % (b["release_path"], e))
            continue
        if not licence_ok(new) or not licence_grammar_ok(new):
            chk.fail("L1", "%s: classify returned %r (one grammar: own_derived | third_party_derived:A+B+...)" % (b["release_path"], new))
        if {"TWB", "TWB2-BBJ"} <= set(re.split(r"[:+|]", new)):       # V06: TWB is the cohort of the TWB2-BBJ meta-analysis, never a provider of its own
            chk.fail("L4", "%s: class holds both TWB and TWB2-BBJ" % b["release_path"])
        b["licence_class"], b["licence_evidence"] = new, ev
    check_module_classes(ctx, built, chk)

PROV_LABELS = ("full run-time hash match", "partial launch-record match", "timestamp evidence only")      # S2-D37 D4: exactly these three
DEPOSIT_PREFIXES = ("paper3_deposit_request_v4/", "paper3_deposit_request_v5/")


def check_code_labels(ctx, rows, chk):
    """K4 (S2-D37 D4, S6.3b2a): the provenance label of a CODE_INCLUDE row is one of the three labels and is checked against the file it describes.
    full run-time hash match: evidence_md5 is a full md5 and equals the md5 of the original file now (the run recorded it; the REPORT's Rule 6 and the delivered
    checkpoints are the record, read by the row's author); partial launch-record match: evidence_md5 is the logged prefix and the original starts with it;
    timestamp evidence only: no hash is claimed. Every shipped row from a deposit package carries a label and a code_rights value; a row that is not shipped carries none."""
    n = collections.Counter()
    for r in rows:
        p, lab, ev, evid, cr = r["path"], r.get("provenance_label", ""), r.get("evidence_md5", ""), r.get("provenance_evidence", ""), r.get("code_rights", "")
        shipped = r["ship"] == "1"
        if not shipped:
            if lab or ev or evid or cr:
                chk.fail("K4", "CODE_INCLUDE %s: not shipped, but it carries a provenance label, evidence or code_rights" % p)
            continue
        if p.startswith(DEPOSIT_PREFIXES) and not lab:
            chk.fail("K4", "CODE_INCLUDE %s: a shipped deposit-package script has no provenance label" % p)
        if not lab:
            if ev or evid or cr:
                chk.fail("K4", "CODE_INCLUDE %s: evidence or code_rights without a provenance label" % p)
            continue
        if lab not in PROV_LABELS:
            chk.fail("K4", "CODE_INCLUDE %s: provenance label %r is not one of %s" % (p, lab, list(PROV_LABELS)))
            continue
        if not evid:
            chk.fail("K4", "CODE_INCLUDE %s: a provenance label with no stated evidence" % p)
        if not (cr == "authors" or cr.startswith("third_party:")):
            chk.fail("K4", "CODE_INCLUDE %s: code_rights must be 'authors' or 'third_party:<what>' (MIT only where the authors hold the rights), got %r" % (p, cr))
        src = ctx.resolve(p)
        if not src.is_file():
            continue                                      # the missing file is K2's finding
        have = md5_file(src)
        if lab == PROV_LABELS[0]:
            if not (re.fullmatch(r"[0-9a-f]{32}", ev) and ev == have):
                chk.fail("K4", "CODE_INCLUDE %s: labelled %r, but the recorded full md5 %r is not the md5 of the file (%s)" % (p, lab, ev, have))
        elif lab == PROV_LABELS[1]:
            if not (re.fullmatch(r"[0-9a-f]{8}", ev) and have.startswith(ev)):
                chk.fail("K4", "CODE_INCLUDE %s: labelled %r, but the logged prefix %r does not start the md5 of the file (%s)" % (p, lab, ev, have))
        elif ev:
            chk.fail("K4", "CODE_INCLUDE %s: labelled %r, which claims no hash, but evidence_md5 is %r" % (p, lab, ev))
        n[lab] += 1
    chk.info("code provenance labels: %s" % ", ".join("%s %d" % (k, n[k]) for k in PROV_LABELS))


def copy_code_bundle(ctx, built, chk):
    """Copy the ship=1 rows of inputs/CODE_INCLUDE.tsv (named files, never globs) to out/<ship_path>; each is a built row."""
    taken = {b["release_path"] for b in built}
    n = 0
    check_code_labels(ctx, read_dicts(ctx.inputs / "CODE_INCLUDE.tsv")[1], chk)
    for r in read_dicts(ctx.inputs / "CODE_INCLUDE.tsv")[1]:
        if r["ship"] != "1":
            continue
        src, rp = ctx.resolve(r["path"]), norm(r["ship_path"])
        if not src.is_file():
            chk.fail("K2", "CODE_INCLUDE: %s does not exist" % r["path"])
            continue
        if rp in taken or not rp.startswith("code/"):
            chk.fail("K2", "CODE_INCLUDE: ship_path %s collides or is outside code/" % rp)
            continue
        taken.add(rp)
        dst = ctx.out / rp
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        built.append(dict(release_path=rp, md5=md5_file(dst), source_path=norm(r["path"]), source_md5=ctx.md5(src), producer="P31R_build_release.py:copy_code",
                          licence_class="own_derived", ship_status="ship", input_id="", store="out",
                          note="code bundle (inputs/CODE_INCLUDE.tsv)" + ("; provenance: %s" % r["provenance_label"] if r.get("provenance_label") else "")))
        n += 1
    if n == 0:
        chk.fail("K2", "CODE_INCLUDE.tsv ships no file")


def _load_sibling(ctx, name):
    sys.path.insert(0, str(ctx.release))
    try:
        sys.modules.pop(name, None)
        importlib.invalidate_caches()
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def apply_tiers(ctx, built, chk):
    """RC1 (S2-D33): classify every file of out/ by content (rel_tiers), apply the declared withholding (inputs/TIER_RULES.tsv), and bring the build records in line:
    a file the rules moved whole becomes a withheld row, a reduced file gets its new md5, and the full counterpart in withheld/ gets a row of its own."""
    try:
        rt = _load_sibling(ctx, "rel_tiers")
    except ImportError as e:
        chk.fail("T0", "module rel_tiers cannot be imported: %s" % e)
        return
    tctx = rt.Ctx(ctx.out, ctx.withheld, ctx.build, True)
    rows, unruled = rt.classify(tctx)
    bad = rt.verify_state(tctx, rows)
    hits = rt.leak_scan(tctx, rows, use_removed=True)
    rt.write_outputs(tctx, rows, hits)
    for u in unruled:
        chk.fail("T1", "%s carries a variant id and a value column and no TIER_RULES row covers it" % u)
    for b_ in bad:
        chk.fail("T2", "tier rule not applied: %s" % b_)
    for h in hits:
        if "TPMI source row" in h[2]:
            chk.fail("T3", "%s still holds %r (%s)" % tuple(h))
    # build records
    have = {(b["store"], b["release_path"]) for b in built}
    for b in list(built):
        if b["store"] != "out":
            continue
        pth = ctx.out / b["release_path"]
        if not pth.is_file():
            wh = ctx.withheld / b["release_path"]
            if not wh.is_file():
                chk.fail("T4", "%s left out/ and is not in withheld/" % b["release_path"])
                continue
            if ("withheld", b["release_path"]) in have:
                built.remove(b)         # the withheld counterpart already has its own row
            else:
                b.update(store="withheld", ship_status="withheld_pending_terms", md5=md5_file(wh),
                         note=(b["note"] + "; " if b["note"] else "") + "withheld whole by the licence tier rules (rel_tiers, inputs/TIER_RULES.tsv)")
                have.add(("withheld", b["release_path"]))
            continue
        m = md5_file(pth)
        if m != b["md5"]:
            b["md5"] = m
            b["note"] = (b["note"] + "; " if b["note"] else "") + "reduced by the licence tier rules (rel_tiers, inputs/TIER_RULES.tsv); the full copy is in withheld/"
            b.pop("bytes", None)
    outrows = {b["release_path"]: b for b in built if b["store"] == "out"}
    for f in list_files(ctx.withheld):
        if ("withheld", f) in have:
            continue
        src = outrows.get(f)
        if src is None:
            chk.fail("T4", "withheld/%s has no build record and no shipped counterpart" % f)
            continue
        built.append(dict(release_path=f, md5=md5_file(ctx.withheld / f), source_path=src["source_path"], source_md5=src["source_md5"], producer="rel_tiers",
                          licence_class=src["licence_class"], ship_status="withheld_pending_terms", input_id=src["input_id"], store="withheld",
                          note="full counterpart of out/%s before the licence tier rules removed or binned its third-party values (rel_tiers)" % f))
        have.add(("withheld", f))


def sync_portability(ctx, built, chk):
    """rel_portability rewrote the shipped code copies: the build record of each carries the shipped md5 (the source md5 stays the source's)."""
    mod = sys.modules.get("rel_portability")      # the instance call_module just ran (a fresh import would not carry LAST_UPDATES)
    if mod is None:
        chk.fail("K3", "rel_portability is not loaded")
        return
    upd = {u["release_path"]: u for u in getattr(mod, "LAST_UPDATES", [])}
    n = 0
    for b in built:
        u = upd.get(b["release_path"])
        if u is None or b["store"] != "out":
            continue
        if b["source_md5"] != u["source_md5"] and not u.get("extra"):      # an extra (a copy outside code/) may have been reduced by the tier step before the transform
            chk.fail("K3", "%s: the portability step read source md5 %s, the copy step recorded %s" % (b["release_path"], u["source_md5"][:8], b["source_md5"][:8]))
        b["md5"] = u["md5"]
        if u["rules"] != "none (byte-identical)":
            b["note"] = b["note"] + "; portability transform applied to this copy (rel_portability; code/P31R_PORTABILITY_LOG.tsv maps source md5 to shipped md5)"
        b.pop("bytes", None)
        n += 1
    if n != len(upd):
        chk.fail("K3", "the portability step transformed %d code file(s), %d have a build record" % (len(upd), n))


def run_build(ctx, stage2, ledger_df=None):
    """Run the whole pipeline into ctx.base. Returns the state dict. Raises Refusal on any failure."""
    chk = ctx.chk
    ri = load_ri(ctx)
    rules = load_rules(ctx)
    validate_ri(ctx, ri, chk)
    verify_pins(ctx, ri, chk)
    check_derived_licences(ctx, ri, chk)
    check_notes(ctx, ri, chk)
    check_cite_tokens(ctx, ri, chk)
    providers_text = check_providers(ctx, ri, chk)
    frames, pending, _ = derive_frames(ctx)
    synced, dropped, added = sync_frame_rows(ctx, ri, frames)
    if dropped or added:
        chk.fail("P3", "the frame rows in RELEASE_INPUTS.tsv differ from the registry derivation (extra %s, missing %s): run --pin" % (dropped, added))
    if pending and not stage2:
        for p in pending:
            chk.fail("F3", "unit %s has no plotted frame yet (%s); strict mode refuses, --stage2 lists it" % (p["unit"], p["expected"]))
    ns = load_not_shipped(ctx)
    if ns is None:
        chk.fail("N2", "inputs/NOT_SHIPPED.tsv is missing: run --derive-not-shipped")
        ns = []
    else:
        fresh, ovr = derive_not_shipped(ctx, ri, rules, chk)
        if [dict(r) for r in fresh] != [dict(r) for r in ns]:
            a, b = {r["path"]: r for r in ns}, {r["path"]: r for r in fresh}
            chk.fail("N2", "NOT_SHIPPED.tsv is stale (%d listed, %d derived now; only listed %s; only derived %s): re-run --derive-not-shipped" %
                     (len(a), len(b), sorted(set(a) - set(b))[:3], sorted(set(b) - set(a))[:3]))
    cnt_src = collections.Counter(f["source_path"] for f in frames)
    for sp, n in cnt_src.items():
        if n > 1:
            chk.fail("F2", "frame %s is claimed %d times" % (sp, n))
    chk.raise_if_failed("pre-build")
    if ledger_df is None:
        import pandas as pd
        lrow = [r for r in ri if r["role"] == "input" and r["source_path"].endswith("S1R1_CLAIMS_LEDGER.csv")]
        if len(lrow) != 1:
            raise Refusal("exactly one input row must pin S1R1_CLAIMS_LEDGER.csv")
        ledger_df = pd.read_csv(ctx.resolve(lrow[0]["source_path"]), dtype=str, keep_default_na=False)
    for p in (ctx.out, ctx.withheld, ctx.build):
        ctx.base.mkdir(parents=True, exist_ok=True)
        reset_dir(p, ctx.base)
    fp = Fingerprint()
    built, rule_of = assemble_copies(ctx, ri, chk, fp)
    chk.raise_if_failed("assemble")
    modrows = {}

    def register_module(name, res):
        check_module_classes(ctx, res, chk)
        # F28: the ledger CSV is a source of every file of a module that RELEASE_INPUTS declares it for (the module asserts against / reads it)
        led = {norm(r["source_path"]) for r in ri if r["role"] == "input" and r["source_path"].endswith("S1R1_CLAIMS_LEDGER.csv") and name in r["producer"].split(";")}
        for d in res:
            d["source_paths"] = sorted({norm(x) for x in d["source_paths"]} | led)
            store = "out" if d["ship_status"] == "ship" else "withheld"
            built.append(dict(release_path=norm(d["release_path"]), md5=md5_file(Path(ctx.out if store == "out" else ctx.withheld) / norm(d["release_path"])),
                              source_path=";".join(sorted(norm(x) for x in d["source_paths"])),
                              source_md5=";".join(src_pin(ri, x) for x in sorted(norm(x) for x in d["source_paths"])),
                              producer=d["producer"], licence_class=d["licence_class"], ship_status=d["ship_status"], input_id="", store=store,
                              note=d["note"]))
            for x in d["source_paths"]:
                if not src_pin(ri, x):
                    chk.fail("M1", "%s: source_paths names %s, which is not a declared input" % (d["release_path"], x))
    for name in MODULES:
        if name in LATE_MODULES:
            continue
        res = call_module(ctx, name, ledger_df, ri, chk, stage2)
        chk.raise_if_failed("module " + name)
        modrows[name] = res
        if not res:
            chk.fail("M5", "module %s returned no file" % name)
        register_module(name, res)
    # module outputs the declaration expects (names the spec fixes)
    def expected_outputs(late):
        for r in ri:
            if r["role"] == "module_output" and (bool(set(r["producer"].split(";")) & set(LATE_MODULES)) == late):
                if not any(b["release_path"] == r["release_path"] and b["ship_status"] == r["ship_status"] for b in built):
                    chk.fail("M5", "expected module output %s (%s) was not produced" % (r["release_path"], r["ship_status"]))
    expected_outputs(False)
    for r in ri:
        if r["role"] == "module_output" and not (set(r["producer"].split(";")) & set(LATE_MODULES)):
            got = [b["licence_class"] for b in built if b["release_path"] == r["release_path"] and b["ship_status"] == r["ship_status"]]
            if got and r["licence_class"] != "module_reported" and got[0] != r["licence_class"]:
                chk.warn("L2", "%s (%s): declared licence_class %r, module reports %r" % (r["release_path"], r["ship_status"], r["licence_class"], got[0]))
    rule_of.update(module_rule_clauses(ri))
    # module outputs: require() clauses and fingerprint of the cells they masked cannot be seen, so the assertions are on the shipped file
    for r in ri:
        if r["role"] == "module_output" and r["ship_status"] == "ship" and not (set(r["producer"].split(";")) & set(LATE_MODULES)):
            p = ctx.out / r["release_path"]
            if p.is_file() and delim_of(p):
                hdr, body = read_table(p, delim_of(p))
                require_check(hdr, body, parse_rule(r["withhold_rule"]), r["release_path"], chk)
    chk.raise_if_failed("module outputs")
    check_claims(ctx, ctx.out, ctx.withheld, built, chk)
    stats = scan_shipped(ctx, ctx.out, rule_of, fp, chk, stage2)
    stats["pair_keyed_rows"] = scan_derived(ctx, ctx.out, chk)
    n_out = sum(1 for b in built if b["store"] == "out")
    if stats["files_scanned"] != n_out:
        chk.fail("E1", "the leak scan opened %d file(s) but %d are shipped (a gate's scope is the gate)" % (stats["files_scanned"], n_out))
    if fp.n()["tpmi"] == 0 or fp.n()["h9"] == 0:
        chk.fail("S7", "the leak scan has an empty fingerprint %s: it would pass anything" % fp.n())
    held = sorted({r for c, s, r in chk.items if c == "K1" and s == "HOLD"})
    if held:
        res_text = (ctx.r2 / "sections_R2" / "P31R_Results_R2_draft.md").read_text(encoding="utf-8")
        modfiles = {norm(d["release_path"]) for res in modrows.values() for d in res}
        for rp in held:
            txt = (ctx.out / rp).read_text(encoding="utf-8")
            tokens = sorted(set(re.findall(r"\[HOLD [A-Za-z0-9]+(?: [A-Za-z0-9]+)?", txt)))
            for t in tokens:
                in_res = t in res_text
                if rp in modfiles:
                    chk.warn("K1", "%s carries %s (%s in Results); the module enforces its own strict hold rule" % (rp, t, "also held" if in_res else "NOT in Results by this token"))
                elif not in_res and not stage2:
                    chk.fail("K1", "%s carries %s, which Results does not carry: strict mode refuses" % (rp, t))
                else:
                    chk.warn("K1", "%s carries %s" % (rp, t))
    check_not_shipped(ctx, ns, ri, ctx.out, ctx.withheld, chk)
    texts = surface_texts(ctx)
    keys = {r["input_id"]: r["cite_keys"] for r in ri}
    outkeys = {(r["release_path"], r["ship_status"]): r["cite_keys"] for r in ri if r["role"] in ("copy", "module_output")}
    for b in built:
        k = outkeys.get((b["release_path"], b["ship_status"])) or (Path(b["release_path"]).name if not b["input_id"] else "")
        b["cited_by"] = cited_by(texts, k)
        b["bytes"] = str(Path(ctx.out if b["store"] == "out" else ctx.withheld, b["release_path"]).stat().st_size)
        b["rows"], b["columns"] = shape(Path(ctx.out if b["store"] == "out" else ctx.withheld, b["release_path"]))
    shipped_cb = {b["release_path"]: b["cited_by"] for b in built if b["store"] == "out"}
    for b in built:
        if b["store"] == "withheld" and not b["cited_by"]:
            b["cited_by"] = shipped_cb.get(b["release_path"], "")
    drift = [r["input_id"] for r in ri if r["role"] == "copy" and r["cited_by"] != cited_by(texts, r["cite_keys"])]
    if drift:
        chk.warn("X1", "cited_by in RELEASE_INPUTS differs from the grep now for %d row(s): %s (re-pin to record)" % (len(drift), drift[:6]))
    check_tsv_lines(ctx, chk)
    # RC1 order: copy -> code bundle -> tiers / withhold -> licence classes -> portability transform -> ledger extract -> regenerability -> docs last
    copy_code_bundle(ctx, built, chk)
    chk.raise_if_failed("code bundle")
    apply_tiers(ctx, built, chk)
    chk.raise_if_failed("licence tiers")
    classify_built(ctx, built, chk)
    chk.raise_if_failed("licence classes")
    for b in built:
        if "bytes" not in b:
            fill_row(ctx, b, texts, outkeys)

    def run_late(name, attrs=None):
        res = call_module(ctx, name, ledger_df, ri, chk, stage2, attrs=attrs)
        chk.raise_if_failed("module " + name)
        modrows[name] = res
        if not res:
            chk.fail("M5", "module %s returned no file" % name)
        register_module(name, res)
        for b in built:
            if "bytes" not in b:
                fill_row(ctx, b, texts, outkeys)
    for name in PRE_DOC_MODULES:
        run_late(name)
        if name == "rel_portability":
            sync_portability(ctx, built, chk)
    chk.raise_if_failed("portability")
    scan_shipped(ctx, ctx.out, rule_of, fp, chk, stage2)      # the leak scan again, over the final out/ (the public ledger and the portability records are new since the first scan)
    chk.raise_if_failed("leak scan (final)")
    apply_tiers(ctx, built, chk)      # again: idempotent, and now tiers_by_file.tsv also covers the files the two steps above added
    chk.raise_if_failed("licence tiers (final)")
    expected_outputs(True)
    chk.raise_if_failed("late module outputs")
    by_unit, by_source = regenerability(ctx, ri, built, rules, chk)
    unit_rows = [dict(zip(UNIT_COLS, _unit_row(d))) for d in by_unit]
    src_rows = [dict(zip(SRC_COLS, _src_row(d))) for d in by_source]
    l2_words = {"released_copy": "masked copy", "extract": "reduced copy"}      # REPRODUCE labels each lifted input as declared in LEVEL2_INPUTS.tsv
    l2_labels = {}
    for r_ in level2_map(ctx).values():
        if r_["kind"] not in l2_words:
            chk.fail("Q3", "LEVEL2_INPUTS.tsv %s: kind %r is neither released_copy nor extract" % (r_.get("l2_id"), r_["kind"]))
        else:
            l2_labels[norm(r_["release_path"])] = l2_words[r_["kind"]]
    chk.raise_if_failed("level-2 labels")
    run_late("rel_docs", attrs=dict(PRIOR_ROWS=[dict(b) for b in built], REGEN_UNIT_ROWS=unit_rows, REGEN_SOURCE_ROWS=src_rows, LEVEL2_LABELS=l2_labels))
    chk.raise_if_failed("build")
    write_build_records(ctx, ri, built, pending, by_unit, by_source, stats, fp, ns, stage2)
    return dict(ri=ri, built=built, pending=pending, by_unit=by_unit, by_source=by_source, stats=stats, fp=fp, ns=ns, rule_of=rule_of)


def src_pin(ri, path):
    for r in ri:
        if r["source_path"] and norm(r["source_path"]) == norm(path):
            return r["source_md5"]
    return ""


def write_build_records(ctx, ri, built, pending, by_unit, by_source, stats, fp, ns, stage2):
    order = {"out": 0, "withheld": 1}
    rows = sorted(built, key=lambda b: (order[b["store"]], b["release_path"]))
    write_table(ctx.build / "BUILT_FILES.tsv", BUILT_COLS, [[b[c] if c in b else "" for c in BUILT_COLS] for b in rows])
    write_table(ctx.build / "regenerability_by_unit.tsv", UNIT_COLS, [_unit_row(d) for d in by_unit])
    write_table(ctx.build / "regenerability_by_source.tsv", SRC_COLS, [_src_row(d) for d in by_source])
    write_table(ctx.build / "pending_render.tsv", ["unit", "expected_frames", "why"], [[p["unit"], p["expected"], p["why"]] for p in pending])
    shipped = [b for b in built if b["store"] == "out"]
    held = [b for b in built if b["store"] == "withheld"]
    lines = ["P31R_build_release.py build log (deterministic: no clock, no absolute path)",
             "mode: %s" % ("stage2 (frames the render has not written are listed, not refused)" if stage2 else "strict"),
             "declared input rows: %d (curated %d, derived frames %d)" % (len(ri), sum(1 for r in ri if r["origin"] == "curated"),
                                                                       sum(1 for r in ri if r["origin"] == "derived_frames")),
             "files shipped (out/): %d; withheld (withheld/): %d; NOT_SHIPPED entries: %d; pending_render entries: %d" % (len(shipped), len(held), len(ns), len(pending)),
             "scan: %s" % dict(sorted(stats.items())),
             "fingerprint forms (strong tokens): %s" % fp.n(),
             "producers (md5): " + "; ".join("%s %s" % (f.name, md5_file(f)[:12]) for f in [Path(__file__)] + [ctx.release / (m + ".py") for m in sorted(MODULES)] if f.is_file()),
             "cited_by grep (spec section 4): a token is a file name or a Table 3 source id; %d operative text surfaces: %s" %
             (len(operative_surfaces(ctx)), ", ".join(sorted(operative_surfaces(ctx)))),
             "providers (inputs/PROVIDERS.tsv, regenerated in memory and asserted equal): %d; terms_status from inputs/PROVIDER_TERMS.tsv" % (len((ctx.inputs / "PROVIDERS.tsv").read_text(encoding="utf-8").splitlines()) - 1),
             "licence classes of shipped files: %s" % dict(sorted(collections.Counter(c.split(":")[0] for b in shipped for c in b["licence_class"].split("|")).items())),
             "regenerability (unit class counts): %s" % dict(sorted(collections.Counter("%s:%s" % (d["kind"], d["class_"]) for d in by_unit).items())),
             ""]
    for cid, sev, msg in ctx.chk.items:
        if sev in ("WARN", "INFO", "HOLD"):
            lines.append("%s %s %s" % (sev, cid, msg))
    (ctx.build / "BUILD_LOG.txt").write_bytes(("\n".join(lines) + "\n").encode("utf-8"))


# ------------------------------------------------------------------------------------------------------------------ commands
def cmd_pin(ctx):
    ri = load_ri(ctx)
    frames, pending, _ = derive_frames(ctx)
    ri, dropped, added = sync_frame_rows(ctx, ri, frames)
    # module inputs that RELEASE_INPUTS does not declare yet: appended UNCLASSIFIED so a human classifies them
    for name in MODULES:
        sys.path.insert(0, str(ctx.release))
        try:
            sys.modules.pop(name, None)
            mod = importlib.import_module(name)
        except ImportError as e:
            print("  note: module %s not importable (%s); its inputs are not checked" % (name, e))
            continue
        finally:
            sys.path.pop(0)
        have = {norm(r["source_path"]): r for r in ri if r["source_path"]}
        for p in mod.INPUTS:
            p = norm(p)
            if p in have:
                if name not in have[p]["producer"].split(";"):
                    have[p]["producer"] = ";".join(sorted(set(have[p]["producer"].split(";")) - {""} | {name}))
                    print("  producer of %s += %s" % (p, name))
            else:
                n = 1 + max([int(r["input_id"][2:]) for r in ri if re.fullmatch(r"RI\d+", r["input_id"])] or [0])
                ri.append(dict(input_id="RI%03d" % n, origin="curated", role="input", deliverable="", release_path="", source_path=p, source_md5="",
                               source_bytes="", licence_class="UNCLASSIFIED", ship_status="", withhold_rule="", producer=name, cite_keys="", cited_by="",
                               note="added by --pin from %s.INPUTS; classify it" % name))
                print("  ADDED input %s for %s (licence UNCLASSIFIED: classify)" % (p, name))
    # licence classes are DERIVED from PROVIDER_DECLARATIONS.tsv (a file inherits the union of its inputs' providers); never set by hand
    pchk = Checks()
    exp = expected_licences(ctx, ri, pchk)
    for r in ri:
        if r["input_id"] in exp and r["licence_class"] != exp[r["input_id"]]:
            print("  licence %s: %s -> %s (derived)" % (r["input_id"], r["licence_class"], exp[r["input_id"]]))
            r["licence_class"] = exp[r["input_id"]]
    for c, m in pchk.failures():
        print("  FAIL [%s] %s" % (c, m))
    texts = surface_texts(ctx)
    changed = 0
    for r in ri:
        if not r["source_path"]:
            continue
        p = ctx.resolve(r["source_path"])
        if not p.is_file():
            print("  ABSENT %s %s" % (r["input_id"], r["source_path"]))
            continue
        md5, sz = ctx.md5(p), p.stat().st_size
        cb = cited_by(texts, r["cite_keys"])
        if (r["source_md5"], r["source_bytes"], r["cited_by"]) != (md5, str(sz), cb):
            changed += 1
            print("  pin %s %s: md5 %s -> %s, bytes %s -> %s%s" % (r["input_id"], r["source_path"], (r["source_md5"] or "-")[:8], md5[:8],
                                                                  r["source_bytes"] or "-", sz, "" if r["cited_by"] == cb else ", cited_by changed"))
        r["source_md5"], r["source_bytes"], r["cited_by"] = md5, str(sz), cb
    for p in dropped:
        print("  frame row removed: %s" % p)
    for p in added:
        print("  frame row added:   %s" % p)
    write_ri(ctx, ri)
    print("pinned %d row(s); %d changed; %d unit(s) pending a render: %s" % (sum(1 for r in ri if r["source_path"]), changed, len(pending), [p["unit"] for p in pending]))
    # --pin never rewrites CANONICAL_SOURCES.tsv; say now if a pin moved a canonical copy (the build refuses it)
    cchk = Checks()
    check_canonical(ctx, ri, cchk)
    for c, m in cchk.failures():
        print("  WARNING [%s] %s" % (c, m))


def cmd_providers(ctx):
    ri = load_ri(ctx)
    chk = ctx.chk
    text = build_providers(ctx, ri, chk)
    chk.raise_if_failed("providers")
    (ctx.inputs / "PROVIDERS.tsv").write_bytes(text.encode("utf-8"))
    print("PROVIDERS.tsv: %d providers; terms_status from PROVIDER_TERMS.tsv" % (len(text.splitlines()) - 1))


def cmd_derive_ns(ctx):
    ri, rules = load_ri(ctx), load_rules(ctx)
    chk = ctx.chk
    verify_pins(ctx, ri, chk)
    chk.raise_if_failed("derive-not-shipped (pins)")
    rows, ovr = derive_not_shipped(ctx, ri, rules, chk)
    chk.raise_if_failed("derive-not-shipped")
    write_table(ctx.inputs / "NOT_SHIPPED.tsv", NS_COLS, [[r[c] for c in NS_COLS] for r in rows])
    by = collections.Counter(r["reason_code"] for r in rows)
    print("NOT_SHIPPED.tsv: %d path+md5 rows %s; hand-declared overrides in use (NOT_SHIPPED_OVERRIDES.tsv, never rewritten here): %d" % (len(rows), dict(sorted(by.items())), len(ovr)))
    for cid, sev, msg in chk.items:
        if sev != "FAIL":
            print("  %s %s %s" % (sev, cid, msg))


def build_atomic(ctx, stage2):
    """Build into a temp folder and move out/ withheld/ build/ into place only if every check passed: a refusal leaves the previous release untouched
    (a half-built out/ must never sit where someone could ship it)."""
    tmp = Path(tempfile.mkdtemp(prefix="p31r_rel_build_"))
    try:
        st = run_build(Ctx(ctx.v3, ctx.release, tmp), stage2)
        for d in ("out", "withheld", "build"):
            target = ctx.base / d
            if target.exists():
                try:
                    shutil.rmtree(target)
                except PermissionError:      # the folder itself is held open (e.g. as a shell's working directory): its CONTENTS must be gone, then the new ones move in
                    left = [x for x in target.iterdir()] if target.exists() else []
                    if left:
                        raise
            if target.exists():
                for x in sorted((tmp / d).iterdir()):
                    shutil.move(str(x), str(target / x.name))
            else:
                shutil.move(str(tmp / d), str(target))
        return st
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def tree_digest(root):
    return {f: md5_file(Path(root) / f) for f in list_files(root)} if Path(root).exists() else {}


def cmd_check(ctx, stage2):
    tmp = Path(tempfile.mkdtemp(prefix="p31r_rel_check_"))
    try:
        c2 = Ctx(ctx.v3, ctx.release, tmp)
        run_build(c2, stage2)
        bad = 0
        for d in ("out", "withheld", "build"):
            a, b = tree_digest(ctx.base / d), tree_digest(tmp / d)
            if a != b:
                bad += 1
                print("DIFF in %s: only here %s; only rebuilt %s; changed %s" % (d, sorted(set(a) - set(b))[:4], sorted(set(b) - set(a))[:4],
                                                                                 [f for f in a if f in b and a[f] != b[f]][:4]))
        if not (ctx.base / "out").exists():
            raise Refusal("nothing to check: %s does not exist (build first)" % (ctx.base / "out"))
        if bad:
            raise Refusal("--check: the rebuild is not byte-identical to out/ withheld/ build/")
        print("--check: out/ (%d files), withheld/ (%d), build/ (%d) are byte-identical to a rebuild" %
              (len(tree_digest(ctx.out)), len(tree_digest(ctx.withheld)), len(tree_digest(ctx.build))))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------------------------------------------------------------------------------ selftest
def _expect(name, fn, cid, results):
    """Run fn() -> Checks; the check ids that FAIL must contain cid (letter+number, exact) and the baseline must not."""
    chk = fn()
    ids = chk.fail_ids()
    ok = cid in ids
    results.append((name, cid, ids, ok))
    return ok


def cmd_selftest(ctx, stage2):
    tmp = Path(tempfile.mkdtemp(prefix="p31r_rel_selftest_"))
    results = []
    try:
        c2 = Ctx(ctx.v3, ctx.release, tmp)
        st = run_build(c2, stage2)
        ri, built, rules = st["ri"], st["built"], load_rules(c2)
        # baseline: every check function is green on the clean build
        chk = Checks(); verify_pins(c2, ri, chk); check_claims(c2, c2.out, c2.withheld, built, chk); check_not_shipped(c2, st["ns"], ri, c2.out, c2.withheld, chk)
        base = Checks(); scan_shipped(c2, c2.out, st["rule_of"], st["fp"], base, stage2)
        results.append(("baseline", "-", chk.fail_ids() + base.fail_ids(), not chk.fail_ids() and not base.fail_ids()))
        # P2: a mutated pin
        def m_pin():
            rows = [dict(r) for r in ri]
            tgt = [r for r in rows if r["role"] == "copy" and r["origin"] == "curated"][0]
            tgt["source_md5"] = "0" * 32
            c = Checks(); verify_pins(c2, rows, c); return c
        _expect("mutated pin", m_pin, "P2", results)
        # C1: an unclaimed file planted in a copy of out/
        def m_claim():
            d = tmp / "mut_claim"; shutil.copytree(c2.out, d); (d / "zz_planted.txt").write_text("x\n")
            c = Checks(); check_claims(c2, d, c2.withheld, built, c); return c
        _expect("unclaimed file", m_claim, "C1", results)
        # S1: a TPMI-derived value planted into a copy of a shipped table
        def m_leak():
            d = tmp / "mut_leak"; shutil.copytree(c2.out, d)
            tok, where = sorted(st["fp"].strong["tpmi"].items())[0]
            tgt = [b["release_path"] for b in built if b["store"] == "out" and delim_of(b["release_path"]) and b["rows"] not in ("", "0")
                   and b["release_path"] != where[0]][0]
            f = d / tgt; txt = f.read_text(encoding="utf-8"); assert txt.count(chr(10)) >= 2
            f.write_text(txt.rstrip(chr(10)) + chr(10) + tok + chr(10), encoding="utf-8")
            c = Checks(); scan_shipped(c2, d, st["rule_of"], st["fp"], c, stage2); return c
        _expect("TPMI value planted", m_leak, "S1", results)
        # N1: a file whose md5 is on the NOT_SHIPPED list
        def m_ns():
            d = tmp / "mut_ns"; shutil.copytree(c2.out, d)
            small = sorted([r for r in st["ns"] if r["md5"] and r["bytes"].isdigit() and int(r["bytes"]) < 2_000_000 and ctx.resolve(r["path"]).is_file()],
                           key=lambda r: int(r["bytes"]))[0]
            shutil.copyfile(ctx.resolve(small["path"]), d / "planted_not_shipped.bin")
            c = Checks(); check_not_shipped(c2, st["ns"], ri, d, c2.withheld, c); return c
        _expect("NOT_SHIPPED md5 in out/", m_ns, "N1", results)
        # H1: an H9 magnitude un-blanked in a copy
        def m_h9():
            d = tmp / "mut_h9"; shutil.copytree(c2.out, d)
            h = h9_spec(c2)
            cand = [b["release_path"] for b in built if b["store"] == "out" and delim_of(b["release_path"])]
            for rp in cand:
                hdr, rows = read_table(d / rp, delim_of(rp))
                src = [c for c in h["source_columns"] if c in hdr]; mag = [c for c in h["magnitude_columns"] if c in hdr]
                if src and mag:
                    hit = [r for r in rows if r[hdr.index(src[0])] == h["source_value"]]
                    if hit:
                        assert sum(1 for r in rows if r is hit[0]) == 1
                        hit[0][hdr.index(mag[0])] = "0.1234"
                        write_table(d / rp, hdr, rows, delim_of(rp)); break
            else:
                raise Refusal("selftest: no shipped table carries an H9 row to mutate")
            c = Checks(); scan_shipped(c2, d, st["rule_of"], st["fp"], c, stage2); return c
        _expect("H9 magnitude un-blanked", m_h9, "H1", results)
        # S4: a masking rule that matches no row
        def m_rule():
            hdr = ["pair", "p"]; rows = [["A", "1"], ["B", "2"]]
            cl = parse_rule("tpmi:pair=NOSUCH=>mask(p)"); c = Checks(); apply_rules(c2, hdr, rows, cl, "mut", c); return c
        _expect("rule matches nothing", m_rule, "S4", results)
        # S3: a require() clause finds an unmasked cell
        def m_req():
            hdr = ["pair", "dev"]; rows = [["A", "phase 3"]]
            c = Checks(); require_check(hdr, rows, parse_rule("ot:*~.=>require(dev)"), "mut", c); return c
        _expect("require() unmasked", m_req, "S3", results)
        # M2: a module that writes outside its namespace
        def m_ns2():
            d = tmp / "mut_mod"; (d / "out" / "evil").mkdir(parents=True); (d / "withheld").mkdir(); (d / "out" / "evil" / "x.tsv").write_text("a\tb\n1\t2\n")
            rec = dict(release_path="evil/x.tsv", ship_status="ship", source_paths=[], producer="rel_screen.py", licence_class="own_derived", rows="1", columns="2", note="")
            c = Checks(); validate_module_result("rel_screen", [rec], d / "out", d / "withheld", set(), c); return c
        _expect("module outside namespace", m_ns2, "M2", results)
        # R1: a source that no rule covers
        def m_rule_none():
            c = Checks(); r = classify_source(rules, "zz_no_such_tree/none.csv");
            if r is None: c.fail("R1", "unmatched source")
            return c
        _expect("source with no rule", m_rule_none, "R1", results)
        def plant(rp, text, rename_header=None):
            d = tmp / ("mut_%02d_" % len(results) + re.sub(r"\W", "_", rp)); shutil.copytree(c2.out, d)
            f = d / rp; txt = f.read_text(encoding="utf-8")
            if rename_header:
                old_, new_ = rename_header
                first, rest = txt.split(chr(10), 1)
                assert first.count(old_) == 1, "selftest: header target not unique"
                txt = first.replace(old_, new_) + chr(10) + rest
            f.write_text(txt.rstrip(chr(10)) + chr(10) + text + chr(10), encoding="utf-8")
            c = Checks(); scan_shipped(c2, d, st["rule_of"], st["fp"], c, stage2); return c
        # RC1b (S2-D34(1)): with SW1 = ship the string scan S2 is off by decision, so the control runs with the switch forced to 'withheld' (the state in which a TPMI mention outside a label is a leak)
        def m_s2_withheld():
            global _TPMI_SHIP_FORCE
            _TPMI_SHIP_FORCE = False
            try:
                return plant("decoy/S1R1_decoy_rates.tsv", "TPMI")
            finally:
                _TPMI_SHIP_FORCE = None
        _expect("TPMI named outside a label (switch withheld)", m_s2_withheld, "S2", results)
        _expect("long quoted passage", lambda: plant("literature/priorart_results.tsv", "see " + chr(34) + " ".join("word%d" % i for i in range(16)) + chr(34)), "Q1", results)
        _expect("verbatim-text column", lambda: plant("decoy/S1R1_decoy_rates.tsv", "x", rename_header=("stratum", "sentence")), "Q2", results)
        # full rebuilds with the defect planted in the DECLARATIONS (the layer where a forgotten rule lives)
        def rebuild(name, mut_ri=None, mut_ns=None, mut_decl=None, mut_text=None):
            rel = tmp / ("rel_" + name)
            shutil.copytree(ctx.release / "inputs", rel / "inputs")
            for m in MODULES:
                shutil.copyfile(ctx.release / (m + ".py"), rel / (m + ".py"))
            cx = Ctx(ctx.v3, rel, tmp / ("base_" + name))
            if mut_ri:
                rows = load_ri(cx); mut_ri(rows); write_ri(cx, rows)
            if mut_ns:
                hdr, rws = read_table(cx.inputs / "NOT_SHIPPED.tsv"); mut_ns(rws); write_table(cx.inputs / "NOT_SHIPPED.tsv", hdr, rws)
            for fn_, f_ in (mut_decl or {}).items():        # a hand-declared table, rows as lists
                hdr, rws = read_table(cx.inputs / fn_); f_(rws); write_table(cx.inputs / fn_, hdr, rws)
            for fn_, f_ in (mut_text or {}).items():
                q_ = cx.inputs / fn_; q_.write_bytes(f_(q_.read_text(encoding="utf-8")).encode("utf-8"))
            try:
                run_build(cx, stage2)
            except Refusal:
                pass
            return cx.chk

        def drop_rule(suffix):
            def f(rows):
                hit = [r for r in rows if r["source_path"].endswith(suffix) and r["role"] == "copy"]
                assert len(hit) == 1 and hit[0]["withhold_rule"], "selftest: mutation target %s is not unique / has no rule" % suffix
                hit[0]["withhold_rule"] = ""
            return f
        _expect("forgotten raw-TPMI mask (outcome checks, RI017)", lambda: rebuild("tpmi", mut_ri=drop_rule("S1R1_outcome_checks.tsv")), "S1", results)   # RC1b: fig2b carries only derived values now (SW1 = ship); the raw per-variant cells of RI017 are the ones a forgotten mask leaks
        _expect("forgotten H9 rule (copy B)", lambda: rebuild("h9", mut_ri=drop_rule("table_s_pair_evidence_tiered.csv")), "H1", results)
        def drop_ns(rws):
            assert len(rws) > 3
            del rws[1]
        _expect("stale NOT_SHIPPED list", lambda: rebuild("ns", mut_ns=drop_ns), "N2", results)
        # ---- fix round 2026-10-18: each new control must trip the check it NAMES (the mutation target is asserted unique first)
        def one(rows, pred, what):
            hit = [r for r in rows if pred(r)]
            assert len(hit) == 1, "selftest: mutation target %s is not unique (%d)" % (what, len(hit))
            return hit[0]

        def m_canon(rows):          # a re-pinned repoint to a non-canonical copy (md5 follows the new file, so P2 stays quiet)
            tgt = one(rows, lambda r: r["role"] == "copy" and Path(r["source_path"]).name == "table_s_pair_evidence_tiered.csv", "tiered copy")
            alt = [x for x in st["ns"] if "NONCANONICAL_COPY" in x["reason_code"] and Path(x["path"]).name == "table_s_pair_evidence_tiered.csv" and x["md5"]][0]
            tgt["source_path"], tgt["source_md5"] = alt["path"], alt["md5"]
        _expect("canonical copy repointed and re-pinned", lambda: rebuild("canon", mut_ri=m_canon), "G1", results)

        def m_ovr(rows):            # the SUPERSEDED override is hand-declared: without its row the derivation refuses to apply it on its own
            rows.remove(one(rows, lambda r: r[2] == "88", "override of manifest line 88"))
        _expect("SUPERSEDED override not declared", lambda: rebuild("ovr", mut_decl={"NOT_SHIPPED_OVERRIDES.tsv": m_ovr}), "N4", results)

        def m_orph(rows):           # an unclaimed frame that is not on the hand-declared orphan list
            del rows[0]
        _expect("orphan frame not declared", lambda: rebuild("orph", mut_decl={"ORPHAN_FRAMES.tsv": m_orph}), "F2", results)

        def m_lic(rows):
            one(rows, lambda r: r["input_id"] == "RI013", "RI013")["licence_class"] = "own_derived"
        _expect("licence class set by hand", lambda: rebuild("lic", mut_ri=m_lic), "L3", results)

        def m_decl(rows):           # the unit's manifest sources lose their provider declaration: no silent default
            rows.remove(one(rows, lambda r: r[1] == "R2_rewrite/figures/data/Figure_5/*", "Figure_5 data declaration"))
        _expect("source without provider declaration", lambda: rebuild("decl", mut_decl={"PROVIDER_DECLARATIONS.tsv": m_decl}), "D3", results)

        def m_pend(rows):
            rows.remove(one(rows, lambda r: r[0] == "OpenTargets", "pending OpenTargets"))
        _expect("provider in neither Table 1 nor Declarations", lambda: rebuild("pend", mut_decl={"PROVIDERS_PENDING_DECLARATIONS.tsv": m_pend}), "D7", results)

        def m_giant_a(rows):        # Table S2 names GIANT: its provider declaration (and pending row) removed, so the resource row names no declared provider
            rows.remove(one(rows, lambda r: r[0] == "GIANT", "GIANT alias"))

        def m_giant_p(rows):
            rows.remove(one(rows, lambda r: r[0] == "GIANT", "GIANT pending"))
        _expect("Table S2 provider missing from the declarations", lambda: rebuild("s2prov", mut_decl={"PROVIDER_ALIASES.tsv": m_giant_a, "PROVIDERS_PENDING_DECLARATIONS.tsv": m_giant_p}), "D8", results)

        def m_prov(t):
            assert t.count("not_found") > 1
            return t.replace("not_found", "read", 1)
        _expect("PROVIDERS.tsv stale", lambda: rebuild("prov", mut_text={"PROVIDERS.tsv": m_prov}), "P4", results)

        def m_note(rows):
            one(rows, lambda r: r["input_id"] == "RI013", "RI013")["note"] += " 1,701 pairs"
        _expect("typed count in a note", lambda: rebuild("note", mut_ri=m_note), "D4", results)

        def m_cite(rows):
            one(rows, lambda r: r["input_id"] == "RI013", "RI013")["cite_keys"] = "15,584"
        _expect("number as a cite key", lambda: rebuild("cite", mut_ri=m_cite), "D6", results)

        def weaken_clause(suffix, clause, weaker):    # the masking rule that the F04 leak lacked, weakened in the declaration (the layer where the defect lives)
            def f(rows):
                hit = [r for r in rows if r["role"] in ("copy", "module_output") and (r["source_path"] or r["release_path"]).endswith(suffix) and clause in r["withhold_rule"]]      # RC1b: RI070 (Table 2) is a module output
                assert len(hit) == 1 and hit[0]["withhold_rule"].count(clause) == 1, "selftest: clause target %s is not unique" % suffix
                hit[0]["withhold_rule"] = hit[0]["withhold_rule"].replace(clause, weaker)
            return f
        # RC1b (S2-D34(1)): the derived FGFR4-T2D (EAS) values ship (TIER_RULES SW1 = ship). The two controls that weakened a mask on a DERIVED TPMI cell (INT pair-list row, extchecks class-rule row)
        # therefore asserted a retired requirement; each is re-aimed, at the same layer (the declared withhold rule), at a mask that still exists: the RAW per-variant cells of RI018, and the
        # ChEMBL-derived Table 2 cells of RI070.
        _expect("raw TPMI per-variant mask weakened (RI018)", lambda: rebuild("int81", mut_ri=weaken_clause(
            "S1R1_R2_11_extchecks/results.tsv", "=>mask(value,se,notes)", "=>mask(notes)")), "S1", results)
        def m_t08(rows):      # the declaration claims the INT class / verdict cells are masked (T08 as it stood before S2-D34(1)) while they ship: the shipped file does not match its own declaration
            r_ = one(rows, lambda r: r[0] == "T08", "DERIVED_ROWS T08")
            r_[4], r_[5] = "derived", "value,p,notes"
        _expect("DERIVED_ROWS declares a shipped TPMI-derived cell as masked (T08)", lambda: rebuild("br15", mut_decl={"DERIVED_ROWS.tsv": m_t08}), "S9", results)

        def plant_derived(rp, edit):
            d = tmp / ("mut_%02d_" % len(results) + re.sub(r"\W", "_", rp)); shutil.copytree(c2.out, d)
            f = d / rp; hdr, rows = read_table(f, delim_of(rp)); edit(hdr, rows); write_table(f, hdr, rows, delim_of(rp))
            c = Checks(); scan_derived(c2, d, c); return c

        # RC1b (S2-D34(1)): the verdict of FGFR4-T2D (EAS) is a derived value and ships, so "verdict un-masked" asserted a retired requirement. Re-aimed at the same row family on the RAW side:
        # the shipped RI018 per-variant row C1.FGFR4_T2D.beta_out carries the resource's own token again.
        def m_raw_cell():
            d = tmp / "mut_rawcell"; shutil.copytree(c2.out, d)
            rp = "external_outcomes/extchecks/results.tsv"
            hdr, rows = read_table(d / rp, delim_of(rp))
            tok = sorted(st["fp"].strong["tpmi"])[0]
            tgt = one([dict(zip(hdr, r), _r=r) for r in rows], lambda r: r["quantity_id"] == "C1.FGFR4_T2D.beta_out", "FGFR4 beta_out row")
            tgt["_r"][hdr.index("value")] = tok
            write_table(d / rp, hdr, rows, delim_of(rp))
            c = Checks(); scan_shipped(c2, d, st["rule_of"], st["fp"], c, stage2); return c
        _expect("raw TPMI per-variant cell un-masked (extchecks beta_out)", m_raw_cell, "S1", results)

        def e_newrow(hdr, rows):    # a pair row no entry classifies (the string 'TPMI' is not in it)
            rows.append(["zz.FGFR4_T2D.new_quantity"] + ["0.5"] * (len(hdr) - 1))
        _expect("pair-keyed row not classified", lambda: plant_derived("workstreams/INT/results.tsv", e_newrow), "S8", results)

        def e_ot(hdr, rows):
            rows[0][hdr.index("development_stage")] = "approved"
        _expect("Open Targets column un-masked", lambda: plant_derived("tables/Table2_evidence_cards.tsv", e_ot), "S8", results)   # RC1b: the OT entries are not_derived now (cells are masked by rel_tables OT1), so the planted drug-stage token is caught as an unclassified OT key

        _expect("H9 decimal on a row naming the held resource", lambda: plant("literature/priorart_results.tsv", "q" + chr(9) + "EGFR_EAS_2P2M" + chr(9) + "0.1234"), "H3", results)
        # K4 (S6.3b2a, D4): the provenance labels are checked against the originals; each mutation edits a copy of the REAL CODE_INCLUDE rows and must trip K4
        ci_rows = read_dicts(c2.inputs / "CODE_INCLUDE.tsv")[1]

        def k4(mutate_rows, label_of):
            rws = [dict(r) for r in ci_rows]
            hit = [r for r in rws if r["ship"] == "1" and r.get("provenance_label") == label_of]
            if not hit:
                raise Refusal("selftest K4: the real CODE_INCLUDE has no shipped %r row to mutate" % label_of)
            mutate_rows(hit[0])
            c = Checks(); check_code_labels(c2, rws, c); return c
        base_k4 = Checks(); check_code_labels(c2, [dict(r) for r in ci_rows], base_k4)
        results.append(("code labels baseline", "-", base_k4.fail_ids(), not base_k4.fail_ids()))
        _expect("label outside the three", lambda: k4(lambda r: r.update(provenance_label="run hash match"), PROV_LABELS[0]), "K4", results)
        _expect("full match, md5 not the file's", lambda: k4(lambda r: r.update(evidence_md5=("0" if r["evidence_md5"][0] != "0" else "1") + r["evidence_md5"][1:]), PROV_LABELS[0]), "K4", results)
        _expect("partial match, prefix not the file's", lambda: k4(lambda r: r.update(evidence_md5=("0" if r["evidence_md5"][0] != "0" else "1") + r["evidence_md5"][1:]), PROV_LABELS[1]), "K4", results)
        _expect("timestamp row claims a hash", lambda: k4(lambda r: r.update(evidence_md5="0" * 32), PROV_LABELS[2]), "K4", results)
        _expect("deposit script without a label", lambda: k4(lambda r: r.update(provenance_label="", evidence_md5="", provenance_evidence="", code_rights=""), PROV_LABELS[2]), "K4", results)
        _expect("label without code_rights", lambda: k4(lambda r: r.update(code_rights=""), PROV_LABELS[0]), "K4", results)
        _expect("label on a row that is not shipped", lambda: k4(lambda r: r.update(ship="0"), PROV_LABELS[1]), "K4", results)

        def m_pending_labelled():
            rws = [dict(r) for r in ci_rows]
            hit = [r for r in rws if r["ship"] == "0" and r["role"] == "deposit_v4_script"]
            if not hit:
                raise Refusal("selftest K4: no pending deposit-v4 row to mutate")
            hit[0]["provenance_label"] = PROV_LABELS[2]
            c = Checks(); check_code_labels(c2, rws, c); return c
        _expect("pending row carries a label", m_pending_labelled, "K4", results)
        # K-id exactness: S1 must not satisfy an S19 expectation
        assert re.fullmatch("S1", "S1") and not re.fullmatch("S1", "S19")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    bad = [r for r in results if not r[3]]
    for name, cid, ids, ok in results:
        print("  %-28s expects %-3s got %-22s %s" % (name, cid, ids, "ok" if ok else "FAILED"))
    print("selftest: %d/%d" % (len(results) - len(bad), len(results)))
    if bad:
        raise Refusal("selftest failed")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--pin", action="store_true")
    g.add_argument("--derive-not-shipped", action="store_true")
    g.add_argument("--check", action="store_true")
    g.add_argument("--selftest", action="store_true")
    g.add_argument("--providers", action="store_true")
    ap.add_argument("--stage2", action="store_true")
    ap.add_argument("--v3", default=str(DEFAULT_V3))
    ap.add_argument("--release", default=str(HERE))
    a = ap.parse_args(argv)
    ctx = Ctx(a.v3, a.release)
    try:
        if a.pin:
            cmd_pin(ctx)
        elif a.derive_not_shipped:
            cmd_derive_ns(ctx)
        elif a.providers:
            cmd_providers(ctx)
        elif a.check:
            cmd_check(ctx, a.stage2)
        elif a.selftest:
            cmd_selftest(ctx, a.stage2)
        else:
            st = build_atomic(ctx, a.stage2)
            print("build ok: %d files shipped, %d withheld, %d not-shipped entries, %d pending render" %
                  (sum(1 for b in st["built"] if b["store"] == "out"), sum(1 for b in st["built"] if b["store"] == "withheld"), len(st["ns"]), len(st["pending"])))
            for d in st["by_unit"]:
                print("  %-20s %-11s %-32s (spec-literal reading: %s)" % (d["unit"], d["kind"], d["class_"], d["class_spec_literal"]))
    except Refusal as e:
        print("REFUSED: %s" % e)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
