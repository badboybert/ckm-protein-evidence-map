"""rel_licence.py: per-file licence class from a released file's CONTENT (Builder L, Stage 5 step 3-4).

classify(path_of_released_file, built_row) -> (licence_class, evidence)

Rule: a released file's class names exactly the providers whose data its VALUES derive from. Evidence order:
  1. F-rows of inputs/PROVIDER_DECLARATIONS.tsv: a file-level declaration (decl_id starts with 'F', pattern = the file's source_path) with its own evidence text.
     Used only where the file has no provider column to read (the evidence column says what was read instead).
  2. The file's own content: column names (decode_*, aric_*, ...), the values of its resource / source / accession / outcome-id / panel columns, and, when a file
     names its outcomes only as CAD/T2D/eGFR..., the screen's declared outcome -> provider table (OUTCOME_PROVIDER; Table 1). Exposure-side providers
     (UKB-PPP, deCODE, ARIC, CKB) are kept from the file's declaration, because the instruments are not a column of the outcome tables.
  3. No resolvable provider column: the declaration (D-rows by source path; for module outputs the module's own class) and the evidence says so.
Distributors (OpenGWAS, GWASCatalog, pheweb.jp) and annotation sources (ChEMBL, DGIdb, OpenTargets, Ensembl) enter only where a value names them.

Nothing here reads a folder name. Providers are the ids of inputs/PROVIDER_ALIASES.tsv; tokens that are not ids are refused.
"""
import csv
import fnmatch
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
OUT = HERE / "out"
csv.field_size_limit(2 ** 31 - 1)

# value tokens (regex over the cell text of a resource-like column) -> providers
VALUE_TOKENS = [
    (r"FinnGen|\bOUT_FG_|finngen", ["FinnGen"]),
    (r"MEGASTROKE", ["MEGASTROKE"]),
    (r"GIGASTROKE", ["GIGASTROKE"]),
    (r"Wuttke|CKDGen|CKDGEN", ["CKDGen"]),
    (r"DIAMANTE", ["DIAMANTE"]),
    (r"CARDIOGRAM", ["CARDIoGRAMplusC4D"]),
    (r"ieu[-_]a[-_]7|OpenGWAS", ["CARDIoGRAMplusC4D", "OpenGWAS"]),
    (r"TPMI", ["TPMI"]),
    (r"\bCLSA\b|GCST9020036\d|GCST9020037\d", ["CLSA"]),
    (r"UKB-?PPP|UKBPPP|UKB_PPP|\bOID\d{4,5}", ["UKB-PPP"]),
    (r"deCODE", ["deCODE"]),
    (r"\bARIC\b", ["ARIC"]),
    (r"\bCKB\b|CKB_", ["CKB"]),
    (r"\bAGEN\b", ["AGEN"]),
    (r"BBJ|Biobank Japan|hum0197|hum0014", ["BBJ"]),
    (r"T2DGGI", ["T2DGGI"]),
    (r"TWB2|Taiwan[–-]Japan|Taiwan Biobank|\bTWB\b", ["TWB2-BBJ"]),
    (r"pheweb|PHEWEB", ["pheweb.jp"]),
    (r"2P2M", ["EGFR_EAS_2P2M_LiuH2025"]),
    (r"GCST90668009", ["JapanHF_GCST90668009"]),
    (r"GCST90018602", ["BBJ"]),
    (r"1000G|1000 Genomes|\bLD50\d|\bLD48\d|g1000_eur|\b1KG\b", ["1000Genomes"]),
    (r"Ensembl", ["Ensembl"]),
    (r"Open ?Targets", ["OpenTargets"]),
    (r"ChEMBL", ["ChEMBL"]),
    (r"DGIdb", ["DGIdb"]),
    (r"GWAS Catalog", ["GWASCatalog"]),
]
VALUE_RX = [(re.compile(p, re.I if p.startswith("CARDIOGRAM") else 0), v) for p, v in VALUE_TOKENS]
HEADER_RX = [(re.compile(p, re.I), v) for p, v in [
    (r"decode", ["deCODE"]), (r"(^|_)aric(_|$)", ["ARIC"]), (r"crossplatform", ["deCODE", "ARIC"]),
    (r"tpmi", ["TPMI"]), (r"(^|_)clsa(_|$)", ["CLSA"]), (r"(^|_)ckb(_|$)", ["CKB"]),
    (r"^(evidence_)?tier$", ["deCODE", "ARIC"])]]   # the evidence tier is graded with the deCODE and ARIC exposure-side arms (Methods l.24)
# columns whose values name a resource, source, accession or outcome file
RES_COLS = re.compile(r"^(resource|res|source|source_file|accession|ea_outcome_accession|eas_outcome_source|outcome_id|protein_arm|protein_source|"
                      r"provider|ext_resource|pair|outcome|outcome_name)$", re.I)
PANEL_COLS = re.compile(r"^(panel|ld_panel)$", re.I)
PANEL_RX = re.compile(r"^(503|489|both|503\+489|503[ (].*|pipeline 1KG.*|1000G.*)$")
# screen design (Table 1; outcome_id / outcome_name columns of screen/screen_ledger.csv): outcome label -> the discovery outcome GWAS
OUTCOME_PROVIDER = {"CAD": "FinnGen", "HF": "FinnGen", "AF": "FinnGen", "T2D": "FinnGen", "NAFLD": "FinnGen", "Stroke": "MEGASTROKE", "CKD": "CKDGen", "eGFR": "CKDGen"}
OUTCOME_TYPE = {"FinnGen", "MEGASTROKE", "GIGASTROKE", "CKDGen", "DIAMANTE", "CARDIoGRAMplusC4D", "BBJ", "AGEN", "T2DGGI", "TPMI", "TWB2-BBJ",
                "pheweb.jp", "EGFR_EAS_2P2M_LiuH2025", "JapanHF_GCST90668009"}
OPAQUE_RX = re.compile(r"^(TAB_|LOG_|ENG_|derived$|definition$|W\d|S6 adjudicator|rebuilt|mr_all_results)")   # a value that names an internal table, whose providers the cell does not show
EXTERNAL_CHECK = {"DIAMANTE", "CARDIoGRAMplusC4D", "TPMI"}   # a file that names these resources reads its outcome labels through them, not through the screen table
EXPOSURE_KEEP = {"UKB-PPP", "deCODE", "ARIC", "CKB"}
DEVSTAGE_RX = re.compile(r"\b(approved|in clinical development|tractab|no agent in the clinic|phase [0-9I]+|clinical stage|discontinued|withdrawn)\b", re.I)


def _tsv(name):
    with open(INPUTS / name, encoding="utf8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _ctx():
    if not _CACHE:
        al = _tsv("PROVIDER_ALIASES.tsv")
        _CACHE["ids"] = [a["provider_id"] for a in al]
        _CACHE["decls"] = _tsv("PROVIDER_DECLARATIONS.tsv")
        _CACHE["byid"] = {d["decl_id"]: d for d in _CACHE["decls"]}
    return _CACHE


_CACHE = {}


def _provs(d):
    return set() if d["providers"] in ("", "none") else set(d["providers"].split("+"))


def _closure(d, byid, seen=()):
    out = _provs(d)
    for i in [x.strip() for x in d["inherits"].split(",") if x.strip()]:
        if i not in seen:
            out |= _closure(byid[i], byid, seen + (d["decl_id"],))
    return out


def _decl(source_path):
    """(decl_id, providers) of the first matching D-row, or the F-row; (None, None) when nothing matches."""
    c = _ctx()
    p = source_path.replace("\\", "/")
    for d in c["decls"]:
        if d["decl_id"].startswith("F") and fnmatch.fnmatchcase(p, d["pattern"]):
            return d["decl_id"], _closure(d, c["byid"]), True
    for d in c["decls"]:
        if not d["decl_id"].startswith("F") and fnmatch.fnmatchcase(p, d["pattern"]):
            return d["decl_id"], _closure(d, c["byid"]), False
    return None, None, False


def _fmt(provs, ids):
    bad = provs - set(ids)
    if bad:
        raise ValueError("classify: provider(s) %s not in PROVIDER_ALIASES.tsv" % sorted(bad))
    return "own_derived" if not provs else "third_party_derived:" + "+".join(i for i in ids if i in provs)


def _of_class(cls):
    return {t for c in cls.split("|") if c.startswith("third_party_derived:") for t in c.split(":", 1)[1].split("+")}


def _read(path):
    p = str(path)
    if not p.lower().endswith((".tsv", ".csv")):
        return None, None
    d = "\t" if p.lower().endswith(".tsv") else ","
    with open(p, encoding="utf8", newline="") as f:
        rows = list(csv.reader(f, delimiter=d))
    return (rows[0], rows[1:]) if rows else (None, None)


def scan(path):
    """(providers, evidence parts) read from the file's own columns; empty when nothing resolves."""
    hdr, rows = _read(path)
    if hdr is None:
        return set(), [], False
    found, ev = set(), []
    opaque = False
    for h in hdr:
        for rx, v in HEADER_RX:
            if rx.search(h):
                found.update(v)
                ev.append("column name %s -> %s" % (h, "+".join(v)))
    strong = set()
    outcome_vals = set()
    for i, h in enumerate(hdr):
        if h.lower() in ("outcome", "outcome_name"):
            outcome_vals |= {r[i] for r in rows if i < len(r)}
        if PANEL_COLS.match(h):
            vals = {r[i] for r in rows if i < len(r)}
            if any(PANEL_RX.match(v) for v in vals):
                strong.add("1000Genomes")
                ev.append("column %s -> 1000Genomes" % h)
        if RES_COLS.match(h):
            vals = {r[i] for r in rows if i < len(r) and r[i]}
            if h.lower() not in ("outcome", "outcome_name", "pair") and any(OPAQUE_RX.search(v) for v in vals):
                opaque = True
            hit = set()
            for v in vals:
                for rx, prov in VALUE_RX:
                    if rx.search(v):
                        hit.update(prov)
            if h.lower() == "provider":
                hit |= {x for v in vals for x in re.split(r"[+;|, ]+", v) if x in _ctx()["ids"]}
            if hit:
                strong |= hit
                ev.append("column %s values -> %s" % (h, "+".join(sorted(hit))))
    found |= strong
    if not (strong & EXTERNAL_CHECK) and outcome_vals & set(OUTCOME_PROVIDER):
        conv = {OUTCOME_PROVIDER[v] for v in outcome_vals if v in OUTCOME_PROVIDER}
        found |= conv
        ev.append("outcome labels %s -> %s (screen outcome table, Table 1)" % ("/".join(sorted(outcome_vals & set(OUTCOME_PROVIDER))), "+".join(sorted(conv))))
    return found, ev, opaque


def _classify_raw(path_of_released_file, built_row):
    """-> (licence_class, evidence). built_row: a row of build/BUILT_FILES.tsv (release_path, source_path, producer, licence_class ...)."""
    ids = _ctx()["ids"]
    p = Path(path_of_released_file)
    if not p.exists():
        p = OUT / str(path_of_released_file)
    old = built_row.get("licence_class", "")
    producer = built_row.get("producer", "")
    srcs = [s for s in built_row.get("source_path", "").split(";") if s]
    module = producer.startswith("rel_")
    # 1. file-level declaration
    fid = fdesc = None
    if len(srcs) == 1:
        did, dp, is_f = _decl(srcs[0])
        if is_f:
            return _fmt(dp, ids), "file declaration %s: %s" % (did, _ctx()["byid"][did]["evidence"])
        declared = dp if did else None
        did_s = did
    else:
        declared, did_s = None, None
    if module:
        declared, did_s = _of_class(old), "module-reported class"
    # tables: development-stage prose names its annotation sources
    rel = str(built_row.get("release_path", p.name))
    if rel.startswith("tables/Table2_evidence_cards") and rel.endswith("Table2_evidence_cards.tsv"):
        base = _of_class(old)
        extra = {"OpenTargets", "ChEMBL", "DGIdb"}
        hdr, rows = _read(p)
        k = hdr.index("established_role_reported_before") if hdr and "established_role_reported_before" in hdr else None
        if k is not None and any(DEVSTAGE_RX.search(r[k]) for r in rows):
            base |= extra
            return "|".join(sorted(["own_derived"] + ["third_party_derived:" + x for x in base])), \
                "module-reported providers plus OpenTargets, ChEMBL, DGIdb: established_role_reported_before carries development-stage statements (%d rows match)" % sum(bool(DEVSTAGE_RX.search(r[k])) for r in rows)
    found, ev, opaque = scan(p)
    if declared is not None and not declared:
        return "own_derived", "declaration %s names no provider (data-free output)" % did_s
    resolved = found - {"1000Genomes"}
    if resolved:
        if opaque:
            add, why = declared or set(), "values name internal tables, so the declaration (%s) is added" % did_s
        else:
            add, why = (declared or set()) & EXPOSURE_KEEP, "exposure side from %s" % did_s
        if add - found:
            ev.append("%s: %s" % (why, "+".join(sorted(add - found))))
        return _fmt(found | add, ids), "content: " + "; ".join(ev)
    if declared is not None:
        return _fmt(declared | found, ids), "declaration (%s): no provider column or value to read in this file%s" % (did_s, "; " + "; ".join(ev) if ev else "")
    return old, "built-row class kept: no content evidence, no file declaration, no module class"


def normalise(cls):
    """V13: ONE licence-class grammar for every row: 'own_derived' (no third-party provider) or 'third_party_derived:A+B+...' (providers in PROVIDER_ALIASES order,
    unknown ids last; no space, no '|'). A '|'-joined class (the Table 2 form) is folded into it."""
    provs = set()
    for c in cls.split("|"):
        if c.startswith("third_party_derived:"):
            provs |= {t.replace("1000 Genomes", "1000Genomes") for t in c.split(":", 1)[1].split("+") if t}
    ids = _ctx()["ids"]
    ordered = [i for i in ids if i in provs] + sorted(provs - set(ids))
    return "third_party_derived:" + "+".join(ordered) if ordered else "own_derived"


def row_class(provs):
    """Class of ONE row from the provider ids of its values (R05 row half): the grammar of normalise()."""
    return normalise("|".join("third_party_derived:" + p.replace("third_party_derived:", "") for p in provs))


def classify(path_of_released_file, built_row):
    """-> (licence_class, evidence); the class always in the one grammar (V13)."""
    cls, ev = _classify_raw(path_of_released_file, built_row)
    return normalise(cls), ev


def table2_row_providers(prior_evidence):      # the row's established_role_reported_before cell (Stage 6.2 name)
    """Annotation sources a Table 2 card row's class must add (R05): the development-stage prose comes from Open Targets, ChEMBL and DGIdb.
    rel_tables.py calls this per row with the row's prior_evidence cell; returns provider ids ('third_party_derived:<id>' tokens are the caller's)."""
    return {"OpenTargets", "ChEMBL", "DGIdb"} if DEVSTAGE_RX.search(prior_evidence or "") else set()


if __name__ == "__main__":   # python -I -X utf8 rel_licence.py : old -> new class of every shipped file in build/BUILT_FILES.tsv, with evidence
    n = ch = 0
    with open(HERE / "build" / "BUILT_FILES.tsv", encoding="utf8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["ship_status"] != "ship" or not (OUT / r["release_path"]).exists():
                continue
            new, ev = classify(OUT / r["release_path"], r)
            n += 1
            if new != r["licence_class"]:
                ch += 1
                print("%s\n  old %s\n  new %s\n  evidence %s" % (r["release_path"], r["licence_class"], new, ev))
    print("%d files classified, %d changed" % (n, ch))
