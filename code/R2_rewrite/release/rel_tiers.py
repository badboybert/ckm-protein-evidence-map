"""rel_tiers.py: licence tier (1 / 2 / 3) of every file in release/out/, decided by CONTENT, and the withholding that follows (Agent T, RC1, 2026-10-18).

Working principle (USER 2026-10-18): minimal but sufficient. The release ships our derived results and our code; third-party per-variant values are pointed to, not redistributed.

  Tier 1  derived per-pair results (MR estimates, colocalization posteriors, counts, classes). Ships with a citation.
  Tier 2  per-variant values copied from a provider: rsID/position rows with the provider's own effect, SE, p or frequency; regional tracks; instrument tables;
          estimates at variants. Ships only for providers whose terms are CONFIRMED CC0 / CC BY in inputs/PROVIDER_TERMS.tsv (status confirmed_from_official_page and a
          licence_name that starts CC0 or CC BY, the same test as the release gate). Otherwise the column / row / file is withheld (full copy in withheld/) or reduced to binned
          plotted values.
  Tier 3  no licence (TPMI): withheld. ONE switch, TIER_RULES.tsv row SW1, default withheld.

What is Tier 1 and what is Tier 2 is a content test, not a folder or file-name test: the authors' own estimate (a Wald ratio and its SE at an instrument, a joint ratio, a posterior)
is Tier 1 even when it sits on a variant row; a provider's association statistic at a variant (beta, SE, p, allele frequency, -log10 p) is Tier 2.

Inputs  : inputs/TIER_RULES.tsv (switch / accession / file rows), inputs/PROVIDER_TERMS.tsv, rel_licence.VALUE_TOKENS (provider tokens in a row's resource cell).
Outputs : build/tiers_by_file.tsv (file, tier, providers, decision, reason), build/tiers_leakscan.tsv, and, when --apply is on, reduced files in out/ with the full copies in withheld/.
Design  : fail closed. A file whose header carries a variant id and a value column (the release gate's R11 detector) and that no rule covers is reported UNRULED and the run exits 2.
          A row-level rule masks a row when the provider named in the row is not open, and also when no provider can be read from the row.
Usage   : python -I -X utf8 rel_tiers.py [--out DIR] [--withheld DIR] [--build DIR] [--apply] [--selftest]
          Without --apply nothing in out/ or withheld/ is touched (classification and the leak scan only; the decisions describe the state AFTER the rules are applied).
"""
import argparse
import csv
import fnmatch
import importlib.util
import io
import re
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
csv.field_size_limit(2 ** 31 - 1)
MASK = "withheld_pending_terms"          # the token rel_tables / rel_licence already use for a withheld cell
NOTE_COLS = ("notes", "note", "comment", "comments")     # S2-D40 (W2): a row whose value / SE / p cell is masked has its notes or comment cell masked too
TIERS_COLS = ["file", "tier", "providers", "decision", "reason"]
# the release gate's R11 header detector (kept identical on purpose; the gate is the authority)
VAR_ID = re.compile(r"(?i)^(rs_?id\w*|snp|snp_?id|variant|variant_?id|signal_variant|marker(?:name)?|chr_?pos|chrpos|pos|position|bp|ukb_variant_id)$")
VAR_VAL = re.compile(r"(?i)(^|_)(beta|se|stderr|pval\w*|p|p_value|maf|eaf|af|freq\w*|z|zscore|or|logp|log10p|nlog\w*|lp)(_|$)")
TABULAR = (".csv", ".tsv")
TEXTUAL = (".csv", ".tsv", ".json", ".md", ".txt", ".log", ".R", ".py")
SELF_EXEMPT = {"RELEASE_MANIFEST.tsv", "README.md", "REPRODUCE.md"}


def _paths(a):
    out = Path(a.out) if a.out else HERE / "out"
    wh = Path(a.withheld) if a.withheld else (out.parent / "withheld" if a.out else HERE / "withheld")
    bd = Path(a.build) if a.build else (out.parent / "build" if a.out else HERE / "build")
    return out, wh, bd


def _tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _licence_tokens():
    spec = importlib.util.spec_from_file_location("rel_licence_for_tiers", HERE / "rel_licence.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.VALUE_RX


def open_providers():
    """Providers whose per-variant values may ship: terms read from an official page AND a licence name that starts CC0 / CC BY (same test as the release gate)."""
    ok = set()
    for r in _tsv(HERE / "inputs" / "PROVIDER_TERMS.tsv"):
        if r["status"].strip() == "confirmed_from_official_page" and re.match(r"(?i)^(CC0|CC[ -]BY)(?![- ]?(?:SA|NC|ND))\b", r["licence_name"].strip()):
            ok.add(r["provider_id"].strip())
    return ok


def load_rules():
    rows = _tsv(HERE / "inputs" / "TIER_RULES.tsv")
    sw = {r["target"]: r["handling"].strip() for r in rows if r["kind"] == "switch"}
    files = [r for r in rows if r["kind"] == "file"]
    acc = {r["target"]: r for r in rows if r["kind"] == "accession"}
    for s in ("TPMI_DERIVED", "CHEMBL_DERIVED", "EAS_COMPARATOR_COLUMNS"):
        if s not in sw:
            raise SystemExit("TIER_RULES.tsv: switch %s is missing" % s)
    for r in files:
        h = r["handling"]
        if not re.match(r"^(columns:[^:]+|rows:[^:]+:[^@]+@\w+|binned:[^:]+:[0-9.]+|file|derived|declared:\w+)$", h):
            raise SystemExit("TIER_RULES.tsv %s: unreadable handling %r" % (r["rule_id"], h))
    return sw, files, acc


def split_ids(s):
    return [x.strip() for x in re.split(r"[+;,|]", s or "") if x.strip()]


def _read_table(path):
    raw = Path(path).read_bytes()
    text = raw.decode("utf-8-sig")
    nl = "\r\n" if "\r\n" in text[:20000] else "\n"
    d = "\t" if str(path).endswith(".tsv") else ","
    rows = list(csv.reader(io.StringIO(text, newline=""), delimiter=d))
    quote_all = text.startswith('"')
    return rows, d, nl, quote_all, raw.startswith(b"\xef\xbb\xbf")


def _write_table(path, rows, d, nl, quote_all, bom):
    buf = io.StringIO(newline="")
    w = csv.writer(buf, delimiter=d, lineterminator=nl, quoting=csv.QUOTE_ALL if quote_all else csv.QUOTE_MINIMAL)
    w.writerows(rows)
    Path(path).write_bytes((b"\xef\xbb\xbf" if bom else b"") + buf.getvalue().encode("utf-8"))


def header_pervar(header):
    hs = [h.strip().strip('"') for h in header]
    return any(VAR_ID.match(h) for h in hs) and any(VAR_VAL.search(h) for h in hs if not VAR_ID.match(h))


def row_providers(text, value_rx):
    ps = []
    for rx, v in value_rx:
        if rx.search(text or ""):
            ps += v
    return ps


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


class Ctx:
    def __init__(self, out, wh, bd, apply):
        self.out, self.wh, self.bd, self.apply = out, wh, bd, apply
        self.open = open_providers()
        self.sw, self.rules, self.acc = load_rules()
        self.vrx = _licence_tokens()
        self.removed = {}       # file -> set of raw strings removed or masked (for the leak scan)
        self.notes_masked = {}  # file -> number of notes / comment cells masked by the W2 rule (counted on the apply run; a verify run adds nothing because nothing is left to mask)
        self.redacts = [r for r in _tsv(HERE / "inputs" / "TIER_RULES.tsv") if r["kind"] == "noteredact"]      # S2-D41 (W4): rule kind 'noteredact'
        self.redacted = {}      # file -> number of notes cells in which a provider value was replaced by the marker
        self.notes = []

    def rule_for(self, rel):
        return [r for r in self.rules if fnmatch.fnmatch(rel, r["target"])]


def apply_one(ctx, rel, rule):
    """Apply one file rule to out/<rel>. Returns (n_cells_or_cols_changed, header_after). Idempotent: a second run changes nothing."""
    path = ctx.out / rel
    h = rule["handling"]
    if h in ("derived",) or h.startswith("declared:"):
        return 0, None
    if h == "file":
        if ctx.apply and path.exists():
            dst = ctx.wh / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(path), str(dst))
        return 1, None
    rows, d, nl, qa, bom = _read_table(path)
    hdr = rows[0]
    gone = ctx.removed.setdefault(rel, set())
    n = 0
    kind, _, rest = h.partition(":")
    if kind == "columns":
        cols = rest.split("|")
        keep = [i for i, c in enumerate(hdr) if c not in cols]
        drop = [i for i, c in enumerate(hdr) if c in cols]
        for r in rows[1:]:
            for i in drop:
                if i < len(r):
                    gone.add(r[i])
        n = len(drop)
        new = [[r[i] if i < len(r) else "" for i in keep] for r in rows]
        hdr_after = new[0]
        if drop and ctx.apply:
            _backup_full(ctx, rel, path)
            _write_table(path, new, d, nl, qa, bom)
        return n, hdr_after
    if kind == "binned":
        cols_s, _, step_s = rest.rpartition(":")
        cols, step = cols_s.split("|"), float(step_s)
        idx = [i for i, c in enumerate(hdr) if c in cols]
        if len(idx) != len(cols):
            raise SystemExit("%s: binned column(s) %s not in header %s" % (rel, cols, hdr))
        new = [hdr]
        for r in rows[1:]:
            r = list(r)
            for i in idx:
                v = _num(r[i])
                if v is not None:
                    b = round(round(v / step) * step, 6)
                    s = repr(b)
                    if s != r[i]:
                        gone.add(r[i])
                        n += 1
                    r[i] = s
            new.append(r)
        if n and ctx.apply:
            _backup_full(ctx, rel, path)
            _write_table(path, new, d, nl, qa, bom)
        return n, hdr
    if kind == "rows":
        sel, _, tail = rest.partition(":")
        idcol, _, rx = sel.partition("~")
        cols_s, _, pcol = tail.partition("@")
        cols = cols_s.split("|")
        rx = re.compile(rx)
        ii = hdr.index(idcol)
        pi = hdr.index(pcol)
        ci = [hdr.index(c) for c in cols if c in hdr]
        ni = [i for i, c in enumerate(hdr) if c in NOTE_COLS]
        new = [hdr]
        for r in rows[1:]:
            r = list(r)
            if rx.search(r[ii]):
                ps = row_providers(r[pi], ctx.vrx)
                if (not ps) or any(p not in ctx.open for p in ps):
                    for i in ci:
                        if r[i] not in ("", "NA", MASK):
                            gone.add(r[i])
                            r[i] = MASK
                            n += 1
                    # S2-D40 (W2), a general rule of every `rows:` handling: when any of the row's value / SE / p cells is masked, its notes or comment cell is masked
                    # too (a notes cell quotes the values it sits beside: beta, SE, p, frequency, z). The provider was read from the row above, before any mask.
                    if any(r[i] == MASK for i in ci):
                        for i in ni:
                            if r[i] not in ("", "NA", MASK):
                                r[i] = MASK
                                n += 1
                                ctx.notes_masked[rel] = ctx.notes_masked.get(rel, 0) + 1
            new.append(r)
        if n and ctx.apply:
            _backup_full(ctx, rel, path)
            _write_table(path, new, d, nl, qa, bom)
        return n, hdr
    raise SystemExit("unreachable handling %r" % h)


def apply_redact(ctx, rel):
    """S2-D41 (W4): rule kind 'noteredact' (target = file glob; handling '<id column>~<id regex>@<notes column>'; selector = a pattern whose capture groups are provider values at a variant).
    In the notes cell of every selected row each capture group is replaced by the withheld marker; the rest of the note stays. Idempotent. Returns the number of cells changed
    (written only when ctx.apply, so a verify run counts what is still unredacted)."""
    rules = [r for r in ctx.redacts if fnmatch.fnmatch(rel, r["target"])]
    path = ctx.out / rel
    if not rules or not path.exists():
        return 0
    rows, d, nl, qa, bom = _read_table(path)
    hdr = rows[0]
    gone = ctx.removed.setdefault(rel, set())
    n = 0
    for rule in rules:
        sel, _, ncol = rule["handling"].partition("@")
        idcol, _, idrx = sel.partition("~")
        _ok = idcol in hdr and ncol in hdr
        if not _ok:
            raise SystemExit("TIER_RULES.tsv %s: column %r or %r is not in the header of %s" % (rule["rule_id"], idcol, ncol, rel))
        ii, ni, irx, vrx = hdr.index(idcol), hdr.index(ncol), re.compile(idrx), re.compile(rule["selector"])

        def sub(m):
            s, off, parts, pos = m.group(0), m.start(), [], 0
            for g in range(1, vrx.groups + 1):
                if m.group(g) is None:
                    continue
                a, b = m.span(g)
                parts.append(s[pos:a - off])
                parts.append(MASK)
                gone.add(m.group(g))
                pos = b - off
            parts.append(s[pos:])
            return "".join(parts)
        for r in rows[1:]:
            if ni < len(r) and ii < len(r) and irx.search(r[ii]) and r[ni] not in ("", "NA", MASK):
                new = vrx.sub(sub, r[ni])
                if new != r[ni]:
                    r[ni] = new
                    n += 1
    if n:
        ctx.redacted[rel] = ctx.redacted.get(rel, 0) + n
        if ctx.apply:
            _backup_full(ctx, rel, path)
            _write_table(path, rows, d, nl, qa, bom)
    return n


def _backup_full(ctx, rel, path):
    dst = ctx.wh / rel
    if not dst.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dst)


def classify(ctx, apply_rules=True):
    """One row per file in out/ (plus one per withheld-only file). The decision describes the file as it will be shipped after the rules."""
    rows, unruled = [], []
    out_files = sorted(p.relative_to(ctx.out).as_posix() for p in ctx.out.rglob("*") if p.is_file()) if ctx.out.exists() else []
    if not ctx.apply:           # classify a would-be state without touching out/: work on the rule table, then read the (unmodified) header
        pass
    for rel in out_files:
        if rel in SELF_EXEMPT:
            continue
        rules = ctx.rule_for(rel)
        pth = ctx.out / rel
        if len(rules) > 1:
            raise SystemExit("%s matches more than one TIER_RULES file row: %s" % (rel, [r["rule_id"] for r in rules]))
        if apply_rules:
            apply_redact(ctx, rel)      # S2-D41 (W4): a file with no file rule can carry a noteredact rule too (workstreams/W3/results.tsv)
        hdr_now = None
        if rel.endswith(TABULAR) and pth.exists():
            hdr_now = _read_table(pth)[0][0]
        if rules:
            r = rules[0]
            n, hdr_after = (0, None)
            if apply_rules:
                n, hdr_after = apply_one(ctx, rel, r)
            h = r["handling"]
            tier, provs = r["tier"].strip(), r["providers"]
            if h == "file":
                dec = "withheld"
            elif h == "derived":
                dec = "ship"
            elif h.startswith("declared:"):
                dec = "ship"
            else:
                kind = h.split(":")[0]
                after = hdr_after if hdr_after is not None else hdr_now
                if ctx.apply is False and kind == "columns":
                    gone = h.partition(":")[2].split("|")
                    after = [c for c in (hdr_now or []) if c not in gone]
                still = bool(after) and header_pervar(after)
                if kind == "binned":
                    dec = "binned_" + re.sub(r"[^0-9a-z]+", "_", h.partition(":")[2].replace(".", "p"))
                elif kind == "rows":
                    dec = "reduced_rows"
                else:
                    dec = ("ship_reduced_columns" if still else "reduced_columns")
            rows.append({"file": rel, "tier": tier, "providers": provs, "decision": dec,
                         "reason": "%s: %s" % (r["rule_id"], r["reason"])})
            continue
        if hdr_now is not None and header_pervar(hdr_now):
            unruled.append(rel)
            rows.append({"file": rel, "tier": "2", "providers": "", "decision": "withheld",
                         "reason": "UNRULED: header carries a variant id and a value column (%s) and no TIER_RULES row covers it" % [h for h in hdr_now if VAR_VAL.search(h)][:3]})
            continue
        own = rel.startswith("code/") or rel.endswith((".R", ".py"))
        rows.append({"file": rel, "tier": "1", "providers": "",
                     "decision": "ship",
                     "reason": "own code" if own else "derived results and counts; header and rules show no provider per-variant value"})
    # withheld-only files
    seen = set(out_files)
    if ctx.wh.exists():
        for p in sorted(ctx.wh.rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(ctx.wh).as_posix()
            if rel in seen:
                continue
            txt = p.read_text(encoding="utf-8", errors="replace") if p.suffix in (".tsv", ".csv", ".tsv", ".log", ".json") else ""
            tier = "3" if re.search(r"TPMI", txt) else "2"
            prov = "TPMI" if tier == "3" else ("OpenTargets;ChEMBL;DGIdb" if re.search(r"opentargets|devstage|annotation", rel, re.I) else "")
            rows.append({"file": rel, "tier": tier, "providers": prov, "decision": "withheld",
                         "reason": "kept in withheld/: full third-party or TPMI-derived values (rel_tables / rel_licence / rel_tiers); the shipped counterpart, if any, is masked or reduced"})
    return rows, unruled


def tpmi_source_values(ctx):
    """Strings of the unmasked TPMI row of the check table (value columns only), read from the stage-1 source file named in the manifest. Empty when it is unreachable."""
    v3 = HERE.parent.parent
    p = v3 / "R1_stage1" / "S1R1_05b_outcome_checks" / "S1R1_outcome_checks.tsv"
    if not p.exists():
        ctx.notes.append("TPMI source row unreachable (%s): TPMI leak scan limited to masked-cell check" % p.name)
        return set()
    vals = set()
    for r in _tsv(p):
        if "TPMI" in (r.get("resource") or ""):
            for k in ("eaf_exposure", "eaf_resource_aligned", "b_out", "se_out", "p_out_file", "wald", "wald_se", "wald_se_delta", "p", "discovery_wald", "discovery_se", "het_z", "het_p"):
                if r.get(k):
                    vals.add(r[k])
    return vals


def _distinctive(x):
    """A string worth searching for: a number of >= 8 characters (sign excluded) that is not the p-value underflow constant. Short values (0.0482) collide by chance."""
    y = x.strip().lstrip("-+")
    v = _num(y)
    return bool(re.fullmatch(r"\d[\d.]*(?:[eE][+-]?\d+)?", y)) and len(y) >= 8 and v is not None and (v == 0 or abs(v) > 1e-300)


def leak_scan(ctx, rows, use_removed=True):
    """Search every shipped text file for the exact strings that were removed or masked (and the TPMI source row). A hit is a provider value that survives in another file."""
    toks = {}
    for rel, s in (ctx.removed.items() if use_removed else ()):
        for x in s:
            if _distinctive(x):
                toks.setdefault(x, set()).add(rel)
    if ctx.sw["TPMI_DERIVED"] == "withheld":
        for x in tpmi_source_values(ctx):
            if _distinctive(x):
                toks.setdefault(x, set()).add("TPMI source row")
    hits = []
    shipped = {r["file"] for r in rows if not r["decision"].startswith("withh")}
    num = re.compile(r"(?<![\w.])\d[\d.]*(?:[eE][+-]?\d+)?(?![\w.]*\d)")      # number-like tokens, sign dropped
    tk = {x.lstrip("-+"): w for x, w in toks.items()}
    for rel in sorted(shipped):
        if rel in SELF_EXEMPT or not rel.endswith(TEXTUAL):
            continue
        p = ctx.out / rel
        if not p.exists():
            continue
        found = set(num.findall(p.read_text(encoding="utf-8", errors="replace")))
        for x in sorted(found & set(tk)):
            if rel in tk[x]:
                continue        # same file: a repeated value (q-value plateau, duplicated instrument), not a copy from another file
            hits.append((rel, x, "|".join(sorted(tk[x]))))
    return hits


def write_outputs(ctx, rows, hits):
    ctx.bd.mkdir(parents=True, exist_ok=True)
    with open(ctx.bd / "tiers_by_file.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        w.writerow(TIERS_COLS)
        for r in rows:
            w.writerow([r[c] for c in TIERS_COLS])
    with open(ctx.bd / "tiers_leakscan.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["shipped_file", "surviving_string", "removed_from"])
        w.writerows(hits)


def verify_state(ctx, rows):
    """After --apply: nothing the rules say to withhold is still present in out/. Returns a list of problems."""
    bad = []
    saved, ctx.apply = ctx.apply, False
    try:
        return _verify(ctx, bad)
    finally:
        ctx.apply = saved


def _verify(ctx, bad):
    for r in ctx.rules:
        h = r["handling"]
        for pth in sorted(ctx.out.glob(r["target"])) if not any(c in r["target"] for c in "[]") else []:
            rel = pth.relative_to(ctx.out).as_posix()
            if h == "file":
                bad.append("%s still in out/ (rule %s withholds it)" % (rel, r["rule_id"]))
            elif h.startswith("columns:"):
                hdr = _read_table(pth)[0][0]
                left = [c for c in h.partition(":")[2].split("|") if c in hdr]
                if left:
                    bad.append("%s still has %s (rule %s)" % (rel, left, r["rule_id"]))
            elif h.startswith("rows:"):
                n, _ = apply_one(ctx, rel, r)
                if n:
                    bad.append("%s: %d cell(s) still unmasked (rule %s)" % (rel, n, r["rule_id"]))
            elif h.startswith("binned:"):
                n, _ = apply_one(ctx, rel, r)
                if n:
                    bad.append("%s: %d value(s) not binned (rule %s)" % (rel, n, r["rule_id"]))
    for rr in ctx.redacts:       # S2-D41 (W4): a notes cell that still prints a provider value at a variant
        for pth in sorted(ctx.out.glob(rr["target"])):
            rel = pth.relative_to(ctx.out).as_posix()
            n = apply_redact(ctx, rel)
            if n:
                bad.append("%s: %d notes cell(s) still print a provider value (rule %s)" % (rel, n, rr["rule_id"]))
    return bad


def run(a):
    out, wh, bd = _paths(a)
    ctx = Ctx(out, wh, bd, a.apply)
    rows, unruled = classify(ctx)
    bad = verify_state(ctx, rows) if a.apply else []
    hits = leak_scan(ctx, rows, use_removed=a.apply)
    write_outputs(ctx, rows, hits)
    from collections import Counter
    c = Counter((r["tier"], r["decision"].split("_")[0]) for r in rows)
    print("rel_tiers: %d files; tier/decision %s; open providers %s" % (len(rows), dict(sorted(c.items())), sorted(ctx.open)))
    print("switches:", ctx.sw)
    for n in ctx.notes:
        print("note:", n)
    for k_, v_ in sorted(ctx.notes_masked.items()):
        print("W2 rule (S2-D40): %d notes cell(s) masked in %s" % (v_, k_))
    for k_, v_ in sorted(ctx.redacted.items()):
        print("W4 rule (S2-D41): %d notes cell(s) with a provider marginal value redacted in %s" % (v_, k_))
    for u in unruled:
        print("UNRULED:", u)
    for b in bad:
        print("NOT APPLIED:", b)
    for h in [h for h in hits if "TPMI source row" in h[2]][:25]:
        print("TPMI LEAK: %s still holds %r (%s)" % h)
    by_file = {}
    for h in hits:
        by_file.setdefault(h[0], 0)
        by_file[h[0]] += 1
    for k, v in sorted(by_file.items(), key=lambda kv: -kv[1])[:8]:
        print("value hit: %s: %d string(s) equal to a withheld value" % (k, v))
    tp = [h for h in hits if "TPMI source row" in h[2]]
    if hits and not tp:
        print("note: %d value hit(s) are informational. A single-instrument MR p-value IS the provider's SNP p-value, so a derived table can repeat a withheld string legitimately; "
              "review build/tiers_leakscan.tsv for hits outside MR-result p columns." % len(hits))
    return 2 if (unruled or bad) else (3 if tp else 0)


def selftest():
    """Copy out/ to a scratch folder, apply, and prove each class of rule bites: a mutation that re-inserts a withheld value must be caught by the leak scan and by verify_state."""
    src = HERE / "out"
    with tempfile.TemporaryDirectory(prefix="rel_tiers_selftest_") as td:
        td = Path(td)
        shutil.copytree(src, td / "out")
        (td / "withheld").mkdir()
        a = argparse.Namespace(out=str(td / "out"), withheld=str(td / "withheld"), build=str(td / "build"), apply=True)
        # the copy has not been tiered yet if out/ is a fresh build; a second run must change nothing (idempotence)
        rc1 = run(a)
        before = {p.relative_to(td / "out").as_posix(): p.read_bytes() for p in (td / "out").rglob("*") if p.is_file()}
        rc2 = run(a)
        after = {p.relative_to(td / "out").as_posix(): p.read_bytes() for p in (td / "out").rglob("*") if p.is_file()}
        fails = []
        if before != after:
            fails.append("second apply changed files: %s" % sorted(k for k in before if before[k] != after.get(k)))
        if rc1 not in (0, 3) or rc2 not in (0, 3):
            fails.append("run exit codes %s / %s" % (rc1, rc2))
        # mutation 1: put a raw instrument column back -> verify_state must say so
        p = td / "out" / "screen" / "cis_pqtl_instruments.csv"
        if p.exists():
            rows, d, nl, qa, bom = _read_table(p)
            rows[0].append("beta")
            for r in rows[1:]:
                r.append("0.123456")
            _write_table(p, rows, d, nl, qa, bom)
            ctx = Ctx(td / "out", td / "withheld", td / "build", True)
            r2, _ = classify(ctx, apply_rules=False)
            if not any("still has" in b and "beta" in b for b in verify_state(ctx, r2)):
                fails.append("mutation 'beta column restored' was not caught by verify_state")
        # mutation 2: an unruled per-variant table must be reported UNRULED
        q = td / "out" / "zz_unruled_probe.tsv"
        q.write_text("rsid\tbeta\nrs1\t0.5\n", encoding="utf-8")
        ctx = Ctx(td / "out", td / "withheld", td / "build", False)
        _, unr = classify(ctx, apply_rules=False)
        if "zz_unruled_probe.tsv" not in unr:
            fails.append("mutation 'unruled per-variant table' was not caught")
        q.unlink()
        # mutation 3: a removed raw string pasted into a shipped markdown file must be a leak hit
        ctx = Ctx(td / "out", td / "withheld", td / "build", True)
        wf = td / "withheld" / "screen" / "cis_pqtl_instruments.csv"
        if wf.exists():
            rows, *_ = _read_table(wf)
            h = rows[0]
            tok = next((r[h.index("pval")] for r in rows[1:] if len(r[h.index("pval")]) >= 8), None)
            if tok is None:
                fails.append("selftest could not pick a p-value probe from the withheld copy")
            else:
                ctx.removed["screen/cis_pqtl_instruments.csv"] = {tok}
                (td / "out" / "workstream_A" / "VERDICT.md").open("a", encoding="utf-8").write("\nprobe %s\n" % tok)
                r3, _ = classify(ctx, apply_rules=False)
                if not any(h_[1] == tok for h_ in leak_scan(ctx, r3, use_removed=True)):
                    fails.append("mutation 'removed value pasted into VERDICT.md' was not caught by the leak scan")
        # mutation 4 (S2-D41 W4): a redacted provider value put back into a notes cell must be caught by verify_state
        w2 = td / "out" / "workstreams" / "W2" / "results.tsv"
        if w2.exists():
            rows, d, nl, qa, bom = _read_table(w2)
            hd = rows[0]
            hit = [r for r in rows[1:] if r[hd.index("quantity_id")].startswith("S2.marg.ratio.") and "marginal z protein" in r[hd.index("notes")]]
            if not hit:
                fails.append("selftest could not find a redacted S2.marg.ratio note in the W2 copy")
            else:
                hit[0][hd.index("notes")] = "marginal z protein 98.7654, outcome 3.21098"
                _write_table(w2, rows, d, nl, qa, bom)
                ctx = Ctx(td / "out", td / "withheld", td / "build", True)
                r4, _ = classify(ctx, apply_rules=False)
                if not any("still print a provider value" in b for b in verify_state(ctx, r4)):
                    fails.append("mutation 'redacted marginal z put back into a W2 note' was not caught by verify_state")
        for f in fails:
            print("SELFTEST FAIL:", f)
        print("rel_tiers selftest: %s" % ("FAIL" if fails else "ok (apply idempotent; 4 mutations each tripped the check they name)"))
        return 1 if fails else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    ap.add_argument("--withheld")
    ap.add_argument("--build")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    sys.exit(selftest() if a.selftest else run(a))
