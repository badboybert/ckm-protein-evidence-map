# -*- coding: utf-8 -*-
"""Stage 5 steps 3-4 -- release module for RELEASE_MANIFEST.tsv, REPRODUCE.md and README.md (Builder D; plan l.115-129).

Writes, at the root of out/ (the core's MODULES entry for this module needs the empty prefix ""):
  RELEASE_MANIFEST.tsv   one row per released file: path, bytes, md5, rows, columns, source_path, source_md5, producer_script, producer_md5,
                         software_label, cited_by, licence_class. It lists every file under out/ EXCEPT itself, README.md and REPRODUCE.md
                         (a file cannot carry its own md5; the README says so).
  REPRODUCE.md, README.md

Interface: INPUTS and build(v3, out, withheld, ledger) -> list[dict] (the rel_screen / rel_tables protocol).
Run it LAST: out/ must already hold every other released file, including the files of the code bundle (out/code/..., copied by the core from the
ship=1 rows of inputs/CODE_INCLUDE.tsv). A missing bundle file, an unrecorded non-code file or a README folder without a description is an error.

What the core must give it (the core resets build/ before the modules run, so the records are passed, not read):
  rel_docs.PRIOR_ROWS       list of dicts of every record built so far (BUILT rows and the module records): release_path, ship_status, md5 (optional),
                            source_path (str) or source_paths (list), source_md5 (optional), producer, licence_class, rows, columns, cited_by (optional).
                            None -> read <v3>/R2_rewrite/release/build/BUILT_FILES.tsv (standalone).
  rel_docs.REGEN_UNIT_ROWS  / REGEN_SOURCE_ROWS   the rows of build/regenerability_by_unit.tsv / _by_source.tsv; None -> read them from build/.

Rules: no new statistic, no typed value. Every number printed in README.md or REPRODUCE.md is a ledger row's own display token, followed by the row id in
square brackets, e.g. "48 files (44 in scripts/scripts, 4 in scripts/R) [L664]" (LEDGER_PRINT_RX). A count that has no ledger row is left out.
Machine-readable lines for the release gate:
  REPRODUCE.md unit table:   "| unit | kind | class | ..." rows, the class token being the one derived in build/regenerability_by_unit.tsv;
  each command:              a fenced line "Rscript --vanilla code/..." or "python -I -X utf8 code/...", followed by bullets
                             "  - input: <V3-relative source path> -> <out:RELEASE_PATH | public download | non-redistributable | author file not in this release | engine output | ...>".
Standalone:  python -I -X utf8 rel_docs.py --build OUT_DIR WITHHELD_DIR       (use a folder outside the project; out/ must already be populated)
             python -I -X utf8 rel_docs.py --selftest
"""
import csv
import hashlib
import io
import json
import os
import re
import sys
from pathlib import Path

PRODUCER = "rel_docs.py"
PRIOR_ROWS = None
REGEN_UNIT_ROWS = None
REGEN_SOURCE_ROWS = None
ALLOW_MISSING_CODE = False       # tests only: a ship=1 bundle file absent from out/ is an error

LEDGER_REL = "R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv"
CODE_REL = "R2_rewrite/release/inputs/CODE_INCLUDE.tsv"
PROV_REL = "R2_rewrite/release/inputs/PROVIDERS.tsv"
TERMS_REL = "R2_rewrite/release/inputs/PROVIDER_TERMS.tsv"
T3MAP_REL = "R2_rewrite/sections_R2/P31R_Table3_repository_rows.tsv"   # Stage 7 A2.4c (S2-D51, review R18): Table 3's repository column
INPUTS = [LEDGER_REL, CODE_REL, PROV_REL, TERMS_REL, T3MAP_REL]
REL_DIR = "R2_rewrite/release"
ACTIVE_REL = "R2_rewrite/figures/af3/P31R_ACTIVE_FIGURES.json"

MANIFEST_COLS = ["path", "bytes", "md5", "rows", "columns", "source_path", "source_md5", "producer_script", "producer_md5", "software_label", "cited_by", "licence_class",
                 "provenance_label", "transform"]      # S6.3b2a (D4): the last two are filled for the code bundle only (source_md5 = the original's md5, md5 = the shipped copy's)
PROV_TABLE = [      # S2-D37 D4: the three labels and what each does not establish (the reviewer note's table, wording kept)
    ("full run-time hash match", "a run recorded the full md5 of the file and it equals the delivered original", "proof that all helpers, settings or inputs are identical"),
    ("partial launch-record match", "only a logged launch-time md5 prefix matches the delivered original", "full-file run-version verification"),
    ("timestamp evidence only", "the file's modification time precedes the outputs of the run, and no run recorded a hash of it (a candidate historical version)", "a demonstrated match to the executed version"),
]
SELF_FILES = ("RELEASE_MANIFEST.tsv", "README.md", "REPRODUCE.md", "LICENSE", "CITATION.cff", ".gitattributes")
VERSION = "v2.0.5-R2"   # Stage 7 A2.6 (USER 2026-10-23); v2.0.0-R2 bumped to v2.0.1-R2 at A2.9c (USER 2026-10-25), then to v2.0.2-R2 at Stage 8 step 4 (USER-approved plan 2026-10-26), then to v2.0.3-R2 at Stage 9 (USER-approved 2026-10-27: the release content changed again), then to v2.0.4-R2 at Stage 10 step 1 (USER 2026-10-28: the shipped Additional file 1 builder changed), then to v2.0.5-R2 at G5 (2026-11-01: the repository URL and the Zenodo concept DOI inserted; the AF1, gate, ledger-mask and tier-rule minors of the step 5f list)
# LICENSE and CITATION.cff text lives in release/templates/, outside the shipped code: release gate R15 reads shipped scripts for licence and
# copyright notices, and the authors' own licence text must not sit in a script. {VERSION} is filled in here and asserted.
# G5 (2026-11-01): the deposit's identifiers. Written by review_R2/s10h23_make_pairs.py from the Zenodo API record (review_R2/zenodo_record_v2.0.4-R2/), never
# typed, and asserted against that record every time the README is built (_readme). CITATION.cff carries no doi key: Zenodo mints a version DOI per GitHub release.
CONCEPT_DOI = "10.5281/zenodo.23855497"
REPO_URL = "https://github.com/badboybert/ckm-protein-evidence-map"
ZENODO_RECORD_REL = "R2_rewrite/review_R2/zenodo_record_v2.0.4-R2/g4_record.json"
ZENODO_RAW_REL = "R2_rewrite/review_R2/zenodo_record_v2.0.4-R2/zenodo_record_23855498_edit1.json"
_TPL = Path(__file__).resolve().parent / "templates"
LICENSE_TEXT = (_TPL / "LICENSE.txt").read_text(encoding="utf-8").replace("{VERSION}", VERSION)
CITATION_TEXT = (_TPL / "CITATION.cff").read_text(encoding="utf-8").replace("{VERSION}", VERSION)
assert "{VERSION}" not in LICENSE_TEXT + CITATION_TEXT and VERSION in LICENSE_TEXT and ('version: "%s"' % VERSION) in CITATION_TEXT
assert "{TITLE}" in LICENSE_TEXT      # Stage 10 step 3 fix round: the LICENSE names the manuscript title, filled at write time from the source map (_ms_title), as the README does
GITATTRIBUTES_TEXT = "* -text\n"     # Stage 10 step 3a (R7): git never rewrites line endings, so a clone equals the zip byte for byte

LEDGER_PRINT_RX = re.compile(r"\[L\d{3}\]")

# software_label: run-time = the version string was written by the run itself; post hoc = recorded afterwards in an environment record (L447; SM Software section).
# Assigned by exact path or directory prefix; every other file is `none_recorded` (no version of its producing run is recorded in L447). Unverified beyond these three.
SOFTWARE_RULES = [("screen/", "run-time"), ("colocalization/phase3_coloc__coloc_results_final.csv", "run-time"),
                  ("finemapping/phase3_susie__susie_grid_results.csv", "post hoc")]
NO_LABEL = "none_recorded"
ACC_RX = re.compile(r"GCST\d+|hum\d{4}(?:\.v\d+)?|ieu-[a-z]-\d+|PheCode [\d.]+")

FOLDERS = {
    "annotation": "annotation outputs (development-stage and target annotation)",
    "code": "the code bundle: analysis-tail and figure code, at its project-relative path (release_records/CODE_INCLUDE.tsv, ship = 1)",
    "release_records": "the authors' own release records: the tier rules, the provider terms, the code-bundle list, and the declared leak-scan coincidences without their values",
    "colocalization": "stored colocalization results and the tiered pair-evidence table",
    "decoy": "the decoy-rate script, its log and its result table",
    "east_asian": "East Asian ancestry-matched results (workstream W7 and the stored East Asian scans)",
    "external_outcomes": "external outcome checks, overlap and LDSC noise correlations; extchecks/ holds the pre-stated external-check workstream",
    "finemapping": "stored SuSiE grid results and the candidate edge audit",
    "ledger": "the public extract of the claims ledger (every number the shipped figure code prints; a row is masked, with its reason, when it holds a raw per-variant TPMI value or a ChEMBL-derived field or cites a source file withheld for another reason; a row that cites a file withheld only because it prints a TPMI line is not masked, because the derived TPMI values ship; evidence_status is defined below)",
    "kidney": "renal classification and the UMOD quantity table",
    "literature": "prior-art search results and the search log",
    "metabolite": "metabolite-locus outputs",
    "prep_stamps": "figure input preparation stamps",
    "proxy": "the partial proxy record (the full cache is on the analysis machine)",
    "replication": "deCODE and ARIC replication results",
    "screen": "the MR screen: one row per protein x outcome pair, non-primary estimates, reason codes and the instrument list",
    "sensitivity_uncited": "a sensitivity table the manuscript does not cite",
    "source_data": "per-figure plotted frames (source_data/figures/<figure>/), the metabolite bridge screen and figure preparation stamps",
    "tables": "Table 1 and Table 2 as data, the column definitions and the Table 2 literature search log",
    "workstream_A": "the instrument verification workstream (reconciliation, verdict and the second-reader outputs)",
    "workstreams": "Stage 1 workstream result and branch tables, one folder per source id (see the map below)",
}
# Table 3 source id -> release folder (record F27). Every folder named is asserted to exist in out/.
SOURCE_FOLDER = [("W1", "workstreams/W1/"), ("W2", "workstreams/W2/"), ("W3", "workstreams/W3/"), ("W4a", "workstreams/W4a/"), ("W4b", "workstreams/W4b/"),
                 ("W5", "workstreams/W5/"), ("W6a", "workstreams/W6a/"), ("W6b", "workstreams/W6b/"), ("W7", "east_asian/"),
                 ("INT", "workstreams/INT/"), ("BI", "workstreams/BI/"), ("extchecks", "external_outcomes/extchecks/")]
CLASS_LEGEND = {
    "regenerates_from_release": "every input of the unit is in this release",
    "needs_public_download": "an input is a public download; the script is shipped, the input is not",
    "needs_nonredistributable_input": "an input is a third-party or individual-level file that is not redistributed, or an author file that is not in this release (each input line says which); the unit is not claimed as regenerable",
    "does_not_regenerate_engine_output": "an input is a stored engine output whose engine ran on the analysis machine; the unit is not claimed as regenerable",
}
INPUT_TAG = {"needs_public_download": "public download", "needs_nonredistributable_input": "non-redistributable",
             "does_not_regenerate_engine_output": "engine output"}
AUTHOR_HELD = {"not_on_approved_list", "not_released_pending_decision", "code_bundle_step3"}     # reasons that name an author file or an unreleased own script, not a provider's terms
AUTHOR_TAG = "author file not in this release"


PER_VARIANT_NOTE = ("The release does not redistribute any provider's summary-statistics file or any individual-level data. "
                    "Individual per-variant provider values do appear, in four ways: (a) where the provider's licence permits; "
                    "(b) where a single variant defines a per-pair estimate or the manuscript itself prints the provider's value at the analysed variant (for the FGFR4-T2D pair, TPMI's effect and standard error at the instrument follow exactly from the shipped Wald ratio and the shipped CKB exposure effect, and TPMI's own p at the instrument ships, as printed in Table 2 and Figure 2b; the CLSA arginine marginal p at the two ARG1 variants ships, as printed in the Results); "
                    "(c) as values plotted in a figure (some binned or rounded); and (d) quoted in shipped code or audit notes (listed in `release_records/LEAKSCAN_COINCIDENCES.tsv`, by file, line and provider, with the reason `quoted_value`; the values are not repeated there). "
                    "Separately, Additional file 1 of the article (Table S3) prints, as instrument reporting, the UKB-PPP lead-variant effect, standard error, -log10(p), F statistic and allele frequency of each colocalization-supported pair; "
                    "this release withholds the same columns in `screen/cis_pqtl_instruments.csv`, which covers every instrument and not only those pairs.")


class DocsError(Exception):
    pass


def _req(cond, msg):
    if not cond:
        raise DocsError(msg)


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def norm(p):
    p = str(p).replace("\\", "/")
    return p[3:] if p.startswith("V3/") else p


def read_tsv(p):
    with open(p, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def shape(p):
    """(rows, columns) of a delimited text file (header excluded), else ("", "")."""
    d = {".csv": ",", ".tsv": "\t"}.get(Path(p).suffix.lower())
    if d is None:
        return "", ""
    with open(p, encoding="utf-8", newline="") as fh:
        rd = list(csv.reader(fh, delimiter=d))
    return (str(max(len(rd) - 1, 0)), str(len(rd[0]))) if rd else ("0", "0")


class Ledger:
    def __init__(self, df):
        self.by = {str(r["id"]): r for r in df.to_dict("records")}

    def P(self, rid):
        """Display token of a ledger row followed by its id: how every count reaches README.md / REPRODUCE.md."""
        r = self.by.get(rid)
        _req(r is not None, "ledger row %s does not exist" % rid)
        d = str(r["display"]).strip()
        _req(d and d != "nan", "ledger row %s has no display" % rid)
        return "%s [%s]" % (d, rid)

    def value(self, rid):
        return json.loads(self.by[rid]["value"])


def _prior(v3):
    if PRIOR_ROWS is not None:
        return list(PRIOR_ROWS)
    return read_tsv(Path(v3) / REL_DIR / "build" / "BUILT_FILES.tsv")


def _regen(v3):
    b = Path(v3) / REL_DIR / "build"
    u = REGEN_UNIT_ROWS if REGEN_UNIT_ROWS is not None else read_tsv(b / "regenerability_by_unit.tsv")
    s = REGEN_SOURCE_ROWS if REGEN_SOURCE_ROWS is not None else read_tsv(b / "regenerability_by_source.tsv")
    return list(u), list(s)


def _code_rows(v3):
    rows = read_tsv(Path(v3) / CODE_REL)
    for r in rows:
        _req(r["ship"] in ("0", "1"), "CODE_INCLUDE ship must be 0 or 1: %s" % r["path"])
        _req((r["ship"] == "1") == bool(r["ship_path"]), "CODE_INCLUDE ship_path iff ship=1: %s" % r["path"])
        _req(r["ship"] == "0" or r["group"] in ("analysis_tail_figure", "deposit_engine"), "a shipped code file must be group analysis_tail_figure or deposit_engine: %s" % r["path"])
        _req(r.get("provenance_label", "") in [""] + [x[0] for x in PROV_TABLE], "CODE_INCLUDE provenance_label is not one of the three labels: %s" % r["path"])
        _req("*" not in r["path"] and "?" not in r["path"], "CODE_INCLUDE declares files by name, never by glob: %s" % r["path"])
    return rows


def _software(path):
    for pre, lab in SOFTWARE_RULES:
        if path.startswith(pre):
            return lab
    return NO_LABEL


def _srcs(rec):
    if rec.get("source_paths") is not None and not isinstance(rec.get("source_paths"), str):
        return [norm(x) for x in rec["source_paths"]]
    s = rec.get("source_path", "")
    return [norm(x) for x in s.split(";") if x] if s else []


def _manifest_rows(v3, out, code):
    v3, out = Path(v3), Path(out)
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())
    files = [f for f in files if f not in SELF_FILES]
    ship = {}
    for r in _prior(v3):
        if r.get("ship_status") == "ship" and r.get("store", "out") == "out":
            ship[norm(r["release_path"])] = r
    cmap = {r["ship_path"]: r for r in code if r["ship"] == "1"}
    for sp, r in cmap.items():
        if not (out / sp).is_file() and not ALLOW_MISSING_CODE:
            raise DocsError("code bundle file %s (CODE_INCLUDE) is not in out/" % sp)
    prod_cache = {}

    def producer(prod):
        name = prod.split(":")[0] or "P31R_build_release.py"
        if name not in prod_cache:
            p = v3 / REL_DIR / name
            _req(p.is_file(), "producer script %s not found" % name)
            prod_cache[name] = md5_file(p)
        return name, prod_cache[name]

    rows = []
    for f in files:
        p = out / f
        lab, tr = "", ""
        if f in ship:
            r = ship[f]
            srcs = _srcs(r)
            smd5 = r.get("source_md5") or ";".join(md5_file(v3 / s) if (v3 / s).is_file() else "" for s in srcs)
            pname, pmd5 = producer(r["producer"])
            lic = r["licence_class"]
            cited = r.get("cited_by", "")
            _req(r.get("md5", "") in ("", md5_file(p)), "%s: record md5 differs from the file" % f)
        elif f in cmap:
            r = cmap[f]
            srcs = [r["path"]]
            smd5 = md5_file(v3 / r["path"])
            plog = out / "code" / "P31R_PORTABILITY_LOG.tsv"
            lg = {x["path"]: x for x in read_tsv(plog)} if plog.is_file() else {}
            want = lg[f]["shipped_md5"] if f in lg and lg[f]["source_md5"] == smd5 else smd5      # a transformed copy is checked against the portability log
            _req(want == md5_file(p), "%s: bundle copy differs from its source %s (and from the portability log)" % (f, r["path"]))
            pname, pmd5 = producer("P31R_build_release.py")
            lic, cited = "own_derived", ""
            lab = r.get("provenance_label", "")
            tr = ("byte-identical to the original" if md5_file(p) == smd5 else
                  "derived from the original (source_md5), not byte-identical: %s" % (lg[f]["rules"] if f in lg else "transform not recorded"))
        else:
            raise DocsError("%s is in out/ but has no build record and is not a CODE_INCLUDE ship_path" % f)
        if f in cmap:      # a code-bundle file (its build record is the first branch): label and transform from CODE_INCLUDE and the portability log
            plog0 = out / "code" / "P31R_PORTABILITY_LOG.tsv"
            lg0 = {x["path"]: x for x in read_tsv(plog0)} if plog0.is_file() else {}
            lab = cmap[f].get("provenance_label", "")
            tr = ("byte-identical to the original" if md5_file(p) == smd5 else
                  "derived from the original (source_md5), not byte-identical: %s" % (lg0[f]["rules"] if f in lg0 else "transform not recorded"))
        rr, cc = shape(p)
        rows.append(dict(path=f, bytes=str(p.stat().st_size), md5=md5_file(p), rows=rr, columns=cc, source_path=";".join(srcs), source_md5=smd5,
                         producer_script=pname, producer_md5=pmd5, software_label=_software(f), cited_by=cited, licence_class=lic, provenance_label=lab, transform=tr))
    return rows


def _write_tsv(path, cols, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])


def _write_text(path, text):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)


def _release_of(src, md5, ship_rows):
    """The out/ path a regeneration input corresponds to (by md5 of the source, then by source path), else None."""
    for r in ship_rows:
        if md5 and (r.get("source_md5") == md5 or r.get("md5") == md5):
            return norm(r["release_path"])
    for r in ship_rows:
        if norm(src) in _srcs(r):
            return norm(r["release_path"])
    return None


_VERB = {}      # Stage 10 step 3a (R4): release path -> True (verbatim), False (masked) or None (no source md5 recorded) for the released copy of a source found by its project path
_MD5 = {}       # Stage 10 step 3 fix round (ii): release path -> (shipped md5, source md5) of every shipped manifest row; a copy label is derived from the md5s, never from a declared kind
_LIFT = {}      # Stage 10 step 3 fix round (iv): release path -> the lift flag of code/P31R_INPUT_MAP.tsv ('1' when the figure code reads the released copy)
STOP_CLAUSE = "; the figure code does not read it"      # printed by _tag on a released copy the generator stops on; section 1 and the V11 check key on it


def _copy_of(src, md5, ship_rows):
    """The out/ file that HOLDS the rows of a source: a verbatim copy (same md5) or the masked copy of that single source. A derived file (a table that drew values
    from several sources, e.g. a column-definition file) is not a copy and is never pointed to (V11)."""
    for r in ship_rows:
        if md5 and r.get("md5") == md5:
            return norm(r["release_path"])
    for r in ship_rows:
        if md5 and r.get("source_md5") == md5 and "mask(" in (r.get("producer") or ""):
            return norm(r["release_path"])
    for r in ship_rows:       # R4: a manifest row whose source_path IS this project path holds the source (verbatim when its md5 is the source's), whatever the transform was named
        if _srcs(r) == [norm(src)]:
            rp_ = norm(r["release_path"])
            _VERB[rp_] = (r.get("source_md5") == r.get("md5")) if r.get("source_md5") else None      # Stage 10 step 3 fix round: never 'verbatim' without a recorded source md5
            return rp_
    return None


BASIS_TEXT = {}
LEDGER_PUBLIC_REL = "ledger/claims_ledger_public.tsv"
LEDGER_PUBLIC_TAG = "out:" + LEDGER_PUBLIC_REL + " (public extract; masked rows carry their reason)"
FORBIDDEN_DOC_PHRASES = ("pending the user's decision",)          # a statement the shipped files contradict (the ledger extract ships); RC1 fix K03


def _pub(p):
    """A source path as it is printed in the release: the redacted relative form of the public ledger (a drive path becomes <external>/..., a project path project-relative)."""
    import rel_ledger
    return norm(rel_ledger.public_source(norm(p), {}))


LEVEL2_LABELS = None             # {release_path: 'masked copy' | 'reduced copy'}; None -> read inputs/LEVEL2_INPUTS.tsv (kind released_copy | extract)
LEVEL2_WORDS = {"released_copy": "masked copy", "extract": "reduced copy"}
_L2 = {}


def _level2_labels(v3):
    if LEVEL2_LABELS is not None:
        return dict(LEVEL2_LABELS)
    p = Path(v3) / REL_DIR / "inputs" / "LEVEL2_INPUTS.tsv"
    if not p.is_file():
        return {}
    rows = read_tsv(p)
    for r in rows:
        _req(r["kind"] in LEVEL2_WORDS, "LEVEL2_INPUTS.tsv %s: kind %r is neither released_copy nor extract" % (r["l2_id"], r["kind"]))
    return {norm(r["release_path"]): LEVEL2_WORDS[r["kind"]] for r in rows}


def _figure_env(out):
    """(R packages, base packages) that the shipped figure code loads, read from the code itself: library(x), x:: and the names listed in a requireNamespace loop."""
    base = {"grid", "graphics", "grDevices", "tools", "utils", "stats", "parallel", "methods"}
    pk = set()
    for f in sorted((Path(out) / "code" / "R2_rewrite" / "figures" / "scripts").glob("*.R")):
        txt = f.read_text(encoding="utf-8", errors="replace")
        pk |= set(re.findall(r"library\((\w[\w.]*)\)", txt))
        pk |= set(re.findall(r"\b([A-Za-z][\w.]*)::", txt))
        for ln in txt.splitlines():
            if "requireNamespace" in ln and "in c(" in ln:
                pk |= set(re.findall(r'"([A-Za-z][\w.]*)"', ln.split("requireNamespace")[0]))
    pk -= {"x", "p"}
    return sorted(pk - base), sorted(pk & base)


def _env_versions(pkgs):
    """(R version string, {package: version}, Python version, PyMuPDF version) read from the machine that builds the release; the build stops if one cannot be read."""
    import subprocess
    code = 'cat(R.version.string, "\\n", sep=""); for (p in commandArgs(TRUE)) cat(p, "\\t", as.character(packageVersion(p)), "\\n", sep="")'
    r = subprocess.run(["Rscript", "--vanilla", "-e", code] + list(pkgs), capture_output=True, text=True, encoding="utf-8")
    _req(r.returncode == 0, "Rscript could not report the package versions: %s" % r.stderr[-300:])
    lines = [x for x in r.stdout.splitlines() if x.strip()]
    ver = {x.split("\t")[0]: x.split("\t")[1] for x in lines[1:]}
    py = subprocess.run(["python", "-c", "import sys, fitz; print(sys.version.split()[0]); print(fitz.VersionBind)"], capture_output=True, text=True, encoding="utf-8")
    _req(py.returncode == 0, "python (with PyMuPDF) could not report its version: %s" % py.stderr[-300:])
    pl = py.stdout.split()
    return lines[0], ver, pl[0], pl[1]


def _multiline_files(out):
    """Shipped .csv/.tsv files (outside code/) in which a cell holds a line break, read from out/ itself."""
    res = []
    for p in sorted(Path(out).rglob("*")):
        rel = p.relative_to(out).as_posix()
        if not p.is_file() or p.suffix not in (".csv", ".tsv") or rel.startswith("code/"):
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        rows = list(csv.reader(io.StringIO(txt, newline=""), delimiter="\t" if p.suffix == ".tsv" else ","))
        if any(("\n" in c or "\r" in c) for r in rows for c in r):
            res.append(rel)
    return res


CODE_SHIPPED = {}      # {project path: ship_path} of the shipped CODE_INCLUDE rows, set by _build (an input that is such a file is in the release as its portable copy)


def _input_line(s, rp):
    """One 'input:' line of REPRODUCE.md. The private ledger is released as the public extract; every other path is shown redacted.
    A lifted input (inputs/LEVEL2_INPUTS.tsv) is labelled with the kind of copy the archive holds (RC1b fix, item e)."""
    if norm(s["source_path"]) == norm(LEDGER_REL):
        return "  - input: %s -> %s" % (norm(s["source_path"]), LEDGER_PUBLIC_TAG)
    cs = CODE_SHIPPED.get(norm(s["source_path"]))
    if cs:        # Stage 10 step 3b (S2-D68): the code-bundle label is derived from the manifest md5s like every other copy label, never printed from a fixed string
        _req(norm(cs) in _MD5, "code-bundle input %s has no shipped manifest row" % cs)
        _cm_, _csm_ = _MD5[norm(cs)]
        _cb_ = ("byte-identical to the original" if (_csm_ and _cm_ == _csm_) else
                "derived from the original by the declared portability transform" if _csm_ else
                "a copy; the manifest records no source md5, so whether it is byte-identical is not established")
        return "  - input: %s -> out:%s (code bundle; the shipped copy is %s)" % (_pub(s["source_path"]), cs, _cb_)
    lab = _L2.get(rp) if rp else None
    if lab == "masked copy" and rp in _MD5 and _MD5[rp][1] and _MD5[rp][0] == _MD5[rp][1]:      # Stage 10 step 3 fix round (ii): a level-two released_copy whose shipped md5 IS the source md5 is verbatim, whatever the registry calls it
        lab = "verbatim copy"
    return "  - input: %s -> %s%s" % (_pub(s["source_path"]), _tag(s, rp), " (%s)" % lab if lab else "")


def _tag(s, rp):
    """The release-standing text of one source row (one wording for section 2 and section 3)."""
    if s["class"] == "regenerates_from_release":
        return "out:%s" % rp if rp else "release (derived file)"
    why = s.get("not_released_because") or "see regenerability record"
    if rp and rp in _VERB:       # R4: the input is a manifest source_path, so it is written as shipped, whatever the author-held reason says
        if rp in _L2:
            return "out:%s" % rp
        if _VERB[rp] is None:
            return "out:%s (copy; the manifest records no source md5, so whether it is verbatim is not established)" % rp
        if _VERB[rp]:
            return "out:%s (verbatim copy)" % rp
        _req(_LIFT, "code/P31R_INPUT_MAP.tsv was not loaded before a masked copy was described")      # Stage 10 step 3 fix round (iv)
        if _LIFT.get(rp) == "1":
            return "out:%s (masked copy, lifted for regeneration: `code/P31R_INPUT_MAP.tsv` marks it lifted)" % rp
        return "out:%s (masked copy%s: `code/P31R_INPUT_MAP.tsv` %s, so the generator stops on this input)" % (rp, STOP_CLAUSE, "marks the copy neither verbatim nor lifted" if rp in _LIFT else "does not list it")
    tag = AUTHOR_TAG if s["class"] == "needs_nonredistributable_input" and why in AUTHOR_HELD else INPUT_TAG[s["class"]]
    return "%s (%s)%s" % (tag, BASIS_TEXT.get(why, why), "; a masked copy of this source is at out:%s" % rp if rp else "")


def _units(v3, code):
    act = json.loads((Path(v3) / ACTIVE_REL).read_text(encoding="utf-8"))
    gen = {f["id"]: f["generator"] for f in act["figures"]}
    ship_code = {r["path"]: r["ship_path"] for r in code if r["ship"] == "1"}
    for g in gen.values():
        _req(g in ship_code, "generator %s is not a shipped CODE_INCLUDE row" % g)
    return gen, ship_code


def _doc_assertions(name, txt, out):
    """RC1 fix K03/K02: a statement of a shipped document that a shipped file contradicts, or an absolute local path, stops the build."""
    import rel_portability
    for ph in FORBIDDEN_DOC_PHRASES:
        _req(ph not in txt, "%s says %r, which the shipped files contradict" % (name, ph))
    lk = rel_portability.leaks(txt)
    _req(not lk, "%s names an absolute local path at line %d: %s" % ((name,) + (lk[0] if lk else (0, ""))))
    for m in re.finditer(r"^  - input: (\S.*?) -> (?!out:|release \(derived)", txt, re.M):
        pth = m.group(1)
        _req(not (Path(out) / pth).is_file() and not (Path(out) / "code" / pth).is_file(),
             "%s lists %s as not in the release, but the release holds it" % (name, pth))


def _reproduce(v3, out, led, code, tables):
    units, sources = tables
    ship_rows = [r for r in _prior(v3) if r.get("ship_status") == "ship" and r.get("store", "out") == "out"]
    gen, ship_code = _units(v3, code)
    _L2.clear()
    _L2.update(_level2_labels(v3))
    _MD5.clear()      # Stage 10 step 3 fix round (ii)/(iv): the md5s of every shipped row and the lift flags of the input map, read before any input line is written
    _MD5.update({norm(r["release_path"]): (r.get("md5") or "", r.get("source_md5") or "") for r in ship_rows})
    _mp = Path(out) / "code" / "P31R_INPUT_MAP.tsv"
    _req(_mp.is_file(), "code/P31R_INPUT_MAP.tsv is not in out/ (rel_portability must run before rel_docs)")
    _LIFT.clear()
    _LIFT.update({norm(r["release_path"]): (r.get("lift") or "").strip() for r in read_tsv(_mp)})
    _req(_LIFT, "code/P31R_INPUT_MAP.tsv has no rows")
    P = led.P
    L = []
    a = L.append
    a("# REPRODUCE\n")
    a("What of this release regenerates, what does not, and why. Every count is a claims-ledger row printed from its display, with the row id in square brackets.\n")
    a("## 0. Reproducibility levels\n")
    a("- Level one, inspect results: every released table, plotted frame and stored result is in this release with its size, md5, source file and producing script (`RELEASE_MANIFEST.tsv`).")
    a("- Level two, regenerate a figure: the figure code in `code/` runs from this release alone (released results, `ledger/claims_ledger_public.tsv`, `code/P31R_INPUT_MAP.tsv`). Only a unit of class `regenerates_from_release` is claimed at this level.")
    a("- Level three, recompute from upstream inputs: the analysis scripts are shipped, but the third-party inputs they read (provider summary statistics, reference-panel genotypes) are not redistributed. Each is pointed to by provider in the README; no unit is claimed at this level.\n")
    lv2 = [u["unit"] for u in units if u["class"] == "regenerates_from_release" and u["kind"] == "figure"]
    lv1 = [u["unit"] for u in units if u["class"] != "regenerates_from_release"]
    a("Level two (figure units): %s." % (", ".join("`%s`" % x for x in lv2) if lv2 else "none"))
    a("Level one only: %s. Each is named with the missing input in section 2.\n" % ", ".join("`%s`" % x for x in lv1))
    a("## 1. Regenerability, per unit\n")
    a("The class of each unit is derived by script from the approved source manifest (the table in section 1; the inputs are named per unit in section 2). Classes:\n")
    for k, v in CLASS_LEGEND.items():
        a("- `%s`: %s." % (k, v))
    by_unit = {}
    for s in sources:
        by_unit.setdefault(s["unit"], []).append(s)
    # V11: ONE list per unit feeds both the unit table cell and section 2 (every source not of class regenerates_from_release, blocked or not)
    def _copy_row(sp):      # the manifest row whose source_path is exactly this project path (a derived file of several sources is not a copy)
        return next((r_ for r_ in ship_rows if _srcs(r_) == [norm(sp)]), None)

    def _readable(r_):      # Stage 10 step 3 fix round (iv): the code reads the released copy only when it is verbatim (md5 == source md5) or lifted (code/P31R_INPUT_MAP.tsv)
        return (bool(r_.get("source_md5")) and r_.get("source_md5") == r_.get("md5")) or _LIFT.get(norm(r_["release_path"])) == "1"
    _cand = lambda x: x["class"] != "regenerates_from_release" and norm(x["source_path"]) != norm(LEDGER_REL) and norm(x["source_path"]) not in CODE_SHIPPED
    missing = {u["unit"]: sorted({_pub(x["source_path"]) for x in by_unit.get(u["unit"], []) if _cand(x) and not (_copy_row(x["source_path"]) and _readable(_copy_row(x["source_path"])))}) for u in units}
    stops = {u["unit"]: {_pub(x["source_path"]) for x in by_unit.get(u["unit"], []) if _cand(x) and _copy_row(x["source_path"]) and not _readable(_copy_row(x["source_path"]))} for u in units}
    a("\n| unit | kind | class | inputs not in the release, or released only as a copy the code does not read |")
    a("|---|---|---|---|")
    for u in units:
        blk = missing[u["unit"]]
        a("| %s | %s | %s | %s |" % (u["unit"], u["kind"], u["class"], "; ".join("`%s`%s" % (b, " (released as a masked copy that the code does not read)" if b in stops[u["unit"]] else "") for b in blk) if blk else "none"))
    a("\nA unit is claimed as regenerable only in the class `regenerates_from_release`. The other classes name the input that is missing; the unit is not claimed. For a deliverable unit (a stored result table), the class names the engine output that the released file is a copy of.\n")
    a("## 2. Inputs per unit\n")
    # Stage 10 step 3 fix round (R6): four regenerations of FigureS1 compared byte by byte with the render (review record review_R2/s10e_opus_check_wf_ae157b0b.json, lens H-r5-r6-scripts, 2026-10-29):
    # all four differed in the creation date and the trailer file identifier; two of them also in the value of an unreferenced integer object (object 31: 278 or 279). The sentence states what that comparison showed.
    a("A regenerated PDF is compared with the render it reproduces (the figure files of the article and the pages of Additional file 3; this release ships no render) by its text and raster, not by md5: "
      "in four regenerations compared byte by byte (`P31R_FigS1.R`; the comparison is recorded in the project's review record, which is not shipped), the files had the same size as the render, and "
      "the bytes that differed were the creation date in the document information, the file identifier in the trailer and, in two of the four, the value of one integer object that nothing references.\n")
    a("Each source is shown with where it stands: `out:<path>` if the release holds it (or a file derived from it), else the reason it is not released.\n")
    # Stage 10 step 4 (S2-D69): the four label families defined for the reader (the 3b refute found them undefined); derived from _copy_of / _tag / LEVEL2_WORDS and the LEVEL2_INPUTS reasons
    a("Labels on the input lines: `verbatim copy`, the release holds the input byte for byte (its manifest md5 is the source md5); `masked copy`, the released file holds that input's rows "
      "with the cells that restate a provider's per-variant value, or that a data hold withholds, masked and the other cells as in the source (its md5 differs from the source md5); "
      "`reduced copy`, the release holds an extract of the input, the rows and columns the generator reads, with provider per-variant values left out; "
      "`code bundle; the shipped copy is ...`, a script of the code bundle, byte-identical to the original or derived from it by the portability transform the manifest's transform column names.\n")
    for u in units:
        a("### %s (%s)\n" % (u["unit"], u["class"]))
        for s in sorted(by_unit.get(u["unit"], []), key=lambda x: x["source_path"]):
            rp = _release_of(s["source_path"], s.get("md5", ""), ship_rows) if s["class"] == "regenerates_from_release" else _copy_of(s["source_path"], s.get("md5", ""), ship_rows)
            a(_input_line(s, rp))
        a("")
    a("## 3. Commands for the code the bundle ships\n")
    a("Run every command from the archive root (the folder that holds `code/`, `ledger/` and the released tables). The shipped copies of the scripts are portable: the scripts that read the archive configuration resolve their paths from one configuration, "
      "the environment variable `P31R_ARCHIVE` (default: the archive root, derived from the script's own place) and `P31R_EXTERNAL` (default: `external_inputs/` beside it, for third-party inputs that are not in the release). "
      "`code/P31R_INPUT_MAP.tsv` maps each input's project path to its path in the archive and says whether the released copy is verbatim. "
      "`code/P31R_PORTABILITY_LOG.tsv` maps the source md5 of each shipped script, which is what the md5 pins inside the scripts name, to the shipped md5. "
      "The figure scripts stop with `input not in the release: <path>` and the reason when asked for an input that is not in the release. The analysis scripts shipped for tracing (the workstream and decoy scripts, the ledger builder) have the same path roots but fail with an ordinary file-not-found error on an input that is not in the release.\n")
    _fs = sorted({ship_code[gen[re.sub(r"_p\d+$", "", u["unit"])]] for u in units if u["kind"] == "figure"})
    _noarch = []
    for _p in _fs:
        _txt = (Path(out) / _p).read_text(encoding="utf-8", errors="replace")
        if "P31R_ARCHIVE" not in _txt:
            _req(re.search(r'source\(file\.path\(dirname\([^)]*\), "P31R_fig_common\.R"\)\)', _txt) is not None, "%s neither reads P31R_ARCHIVE nor sources P31R_fig_common.R from its own folder" % _p)
            _noarch.append(_p.split("/")[-1])
    if _noarch:
        a("%s do not read `P31R_ARCHIVE`: each finds its own folder from the command line and sources `P31R_fig_common.R` from that folder, which derives the figures folder from its own place.\n"
          % ", ".join("`%s`" % x for x in _noarch))
    # Stage 10 step 3 fix round: a path printed <external>/<path> (rel_ledger drops the drive letter) is opened by the portable scripts under P31R_EXTERNAL/<drive letter>/<path>;
    # the pairs are derived from the .P31R_EXT literals of the shipped figure scripts, and every external source must have one (a .tbi index is found by its data file).
    _ext_lit = {}
    for _p in _fs:
        for _m in re.finditer(r'\.P31R_EXT,\s*"/([A-Z])/([^"]+)"', (Path(out) / _p).read_text(encoding="utf-8", errors="replace")):
            _ext_lit[_m.group(2)] = _m.group(1)
    _ext_src = sorted({_pub(s["source_path"])[len("<external>/"):] for s in sources if _pub(s["source_path"]).startswith("<external>/")})
    _ext_of = lambda x: _ext_lit.get(x) or _ext_lit.get(x[:-4] if x.endswith(".tbi") else x)
    for _x in _ext_src:
        _req(_ext_of(_x), "external input %s has no .P31R_EXT literal in a shipped figure script" % _x)
    if _ext_src:
        a("A path printed `<external>/<path>` is a path on a drive of the analysis machine with the drive letter removed. The portable scripts open it below `P31R_EXTERNAL` under a folder named for that drive letter, as their text shows: %s.\n"
          % "; ".join("`<external>/%s` is read from `external_inputs/%s/%s`" % (x, _ext_of(x), x) for x in _ext_src))
    nlf = sorted(norm(r["release_path"]) for r in _prior(v3) if "normalise_newlines" in (r.get("producer") or ""))
    if nlf:
        a("A tab-separated frame that holds a line break inside a cell is shipped rewritten (declared transform of the release builder, `normalise_newlines`): each break inside a cell is replaced by ` / `, "
          "the quotes that this makes unnecessary are dropped, and the rows end with a line feed alone. Frames so rewritten, whose every row sits on one line: %s. "
          "A fresh render of such a frame differs from the shipped file in those three ways, so it differs in md5; compare the two by cells with a CSV reader, after replacing each break inside a cell by ` / `. "
          "The other frames are shipped as rendered: most comma-separated frames end their rows with a carriage return and a line feed, and the ledger and the rewritten frames with a line feed." % ", ".join("`%s`" % x for x in nlf))
        ml = _multiline_files(out)
        if ml:
            a("Frames shipped as rendered that hold line breaks inside quoted cells (so a row can span lines; read them with a CSV reader): %s." % ", ".join("`%s`" % x for x in ml))
        a("")
    a("The commands run `python -I` and `Rscript --vanilla`. A command is listed for every active figure unit; only a unit of class `regenerates_from_release` is claimed as regenerable. The figure code reads its printed numbers from `ledger/claims_ledger_public.tsv`; a figure whose numbers include a masked row stops with the reason.\n")
    seen = {}
    for u in units:
        if u["kind"] != "figure":
            continue
        base = re.sub(r"_p\d+$", "", u["unit"])
        seen.setdefault(base, []).append(u)
    for base, us in seen.items():
        path = ship_code[gen[base]]
        a("```\nRscript --vanilla %s\n```" % path)
        a("Unit(s): %s. Class: %s." % (", ".join(x["unit"] for x in us), ", ".join(sorted({x["class"] for x in us}))))
        lines = []          # RC1b fix Q11: a script that makes two units lists each input once
        for u in us:
            for s in sorted(by_unit.get(u["unit"], []), key=lambda x: x["source_path"]):
                rp = _release_of(s["source_path"], s.get("md5", ""), ship_rows) if s["class"] == "regenerates_from_release" else _copy_of(s["source_path"], s.get("md5", ""), ship_rows)
                ln = _input_line(s, rp)
                if ln not in lines:
                    lines.append(ln)
        for ln in lines:
            a(ln)
        a("")
    a("Aggregators (shipped for inspection, not classified by the regenerability table and not claimed as regenerable; they read the claims ledger, of which the release holds the public extract only, and the manuscript sections, which this release does not contain; the Additional file 1 builder also imports the manuscript builder `P31R_build.py`, which is not shipped):\n")
    for r in code:
        if r["ship"] == "1" and r["role"] in ("render_driver", "af1_builder", "assembler", "legend_builder"):
            cmd = "python -I -X utf8 %s" % r["ship_path"]
            a("- `%s`: %s." % (cmd, r["note"]))
    _agg_noarch = []      # Stage 10 step 3 fix round (R5 scope): the aggregators that do not read P31R_ARCHIVE either, read from their shipped text
    for r_a in code:
        if r_a["ship"] == "1" and r_a["role"] in ("render_driver", "af1_builder", "assembler", "legend_builder"):
            _t_a = (Path(out) / r_a["ship_path"]).read_text(encoding="utf-8", errors="replace")
            if "P31R_ARCHIVE" not in _t_a:
                _req("__file__" in _t_a, "%s neither reads P31R_ARCHIVE nor locates itself by __file__" % r_a["ship_path"])
                _agg_noarch.append(r_a["ship_path"].split("/")[-1])
    if _agg_noarch:
        a("\nOf the aggregators, %s do not read `P31R_ARCHIVE` either: each derives the project root from its own place in the tree.\n" % ", ".join("`%s`" % x for x in _agg_noarch))
    a("\nThe release builder (`code/R2_rewrite/release/`) reads the project tree and the withheld files; it is shipped for inspection and is not rerunnable from this repository alone.\n")
    a("## 4. What does not regenerate\n")
    a("- The screen, colocalization, SuSiE and candidate engines ran on the analysis machine. Their stored outputs are released; their scripts are only partly deposited. Stored screen rows: %s." % P("L660"))
    v4 = [r for r in code if r["role"] == "deposit_v4_script"]
    v4_ship = [r for r in v4 if r["ship"] == "1"]
    v4_not = [r for r in v4 if r["ship"] == "0"]
    renal = [r for r in code if r["path"].endswith("paper3_ckd_rerun_v3/scripts_sent/01_renal_triad_mr.R")]
    _req(len(renal) == 1 and renal[0]["ship"] == "1" and renal[0]["role"] == "reference_script" and renal[0]["provenance_label"], "CODE_INCLUDE: 01_renal_triad_mr.R is not a single shipped, labelled reference_script row")
    _req(len(v4_ship) + len(v4_not) == len(v4) and len(v4) + len(renal) == led.value("L664")[0], "the v4 rows of CODE_INCLUDE plus 01_renal_triad_mr.R differ from the files of L664")
    a("- Deposit v4 delivered %s. `release_records/CODE_INCLUDE.tsv` accounts for those files as follows (the counts are computed from that file when this document is built): `%d` rows of role `deposit_v4_script`, "
      "of which `%d` ship (the v4 scripts that the v5 delivery re-confirmed against the analysis machine) and `%d` do not (ship = 0, each with a reason, pending a decision); and `%s`, which ships as reference code under its "
      "`paper3_ckd_rerun_v3` path without a v5 re-confirmation, with the label `%s`." % (P("L664"), len(v4), len(v4_ship), len(v4_not), renal[0]["path"].split("/")[-1], renal[0]["provenance_label"]))
    a("- The v5 delivery (Part A of the request that asked for these scripts, `REQUEST_PROMPT_v5.md`, which is not shipped) supplied the SuSiE, cascade, candidate-chain and East Asian axis scripts. They are in the code bundle with the provenance labels of the README; the labels say how strongly each file is tied to the run that wrote the results, and what they do not establish.")
    a("- `l3fm__conditional.py`, the conditioning engine named in Supplementary Sections S-M6 and S-M7, is not in deposit v4. The copy in the project (`V4_local/builders/`) is in the code bundle with its helpers; that it is the version that produced the results is unverified.")
    a("- The Stage 1 workstream scripts and the ledger builder with its extension modules are in the code bundle, so that each result can be traced to its script. The ledger builder reads the whole project tree and is not rerunnable from this release.")
    a("- The per-cell screen audit files and the full proxy cache stay on the analysis machine; they are not in this release. Rows of the partial proxy record shipped here: %s." % P("L668"))
    a("- The decoy-rate script (`decoy/S1R1_decoy_rates.R`) is released; its inputs are not. Inputs kept out: %s." % P("L669"))
    a("- Individual-level dosage matrices (`.i8`) found: %s; none is released." % P("L665"))
    n_led = _nrows(out / "screen" / "screen_ledger.csv")
    n_np = _nrows(out / "screen" / "screen_nonprimary_estimates.csv")
    _req(n_led == led.value("L003") == led.value("L004") + sum(led.value("L674")), "screen_ledger.csv rows differ from L003, or from L004 plus the untestable combinations of L674")
    _req(led.value("L660") == led.value("L004") + n_np, "the stored screen rows (L660) differ from the tests with an estimate (L004) plus the rows of screen_nonprimary_estimates.csv")
    a("- Screen counts. `screen/screen_ledger.csv` has one row for each protein-outcome combination of the multiple-testing family, %s: the tests with an estimate, %s, and the untestable combinations (%s) "
      "(the estimates and the untestable combinations add up to the family; checked when this document is built). The stored screen output has %s: every estimator's rows for the tested pairs, which is the primary row of each tested pair plus the "
      "non-primary rows shipped as `screen/screen_nonprimary_estimates.csv` (checked: the rows of that file equal the stored rows minus the tests with an estimate). "
      "The coverage statement (%s) refers to that stored output, not to `screen/screen_ledger.csv`.\n" % (P("L003"), P("L004"), P("L674"), P("L660"), P("L675")))
    a("## 5. Provenance caveats\n")
    a("- %s (deposit v4 report, Part B; the report is not part of this release)." % P("L666"))
    a("- Licence classes are decided per file by content (README, \"Licence tiers\"). " + PER_VARIANT_NOTE + "\n")
    a("## 6. Environment\n")
    a("Versions as recorded: %s." % P("L447"))
    pk, basepk = _figure_env(out)
    rv, pv_, pyv, fitzv = _env_versions(pk)
    a("Level-two figure rendering (the regenerating units of section 0) needs, as installed on the machine that built this release: %s; the R packages %s (read from the shipped figure code); the base packages %s; "
      "and Python %s on the path as `python` with PyMuPDF (`fitz`) %s. `code/R2_rewrite/figures/scripts/P31R_fig_common.R` calls `python` (without `-I`) to write the author metadata into each PDF "
      "before the PNG is rasterized from it; without it a generator stops after it has written the frame and the PDF." % (
          rv, ", ".join("`%s` %s" % (k, pv_[k]) for k in pk), ", ".join("`%s`" % x for x in basepk), pyv, fitzv))
    a("The panels of Figure 4 and FigureS5 that read genotype windows also load `Rsamtools`, `GenomicRanges` and `IRanges` and need inputs that are not in this release (section 2).")
    a("`software_label` in `RELEASE_MANIFEST.tsv`: `run-time` (the version was written by the run itself), `post hoc` (recorded afterwards in an environment record), `%s` (no version of the producing run is recorded in the claims ledger). "
      "A version recorded after the run is not a run-time version.\n" % NO_LABEL)
    txt = "\n".join(L)
    # Stage 10 step 3 fix round (ii): every copy label printed after "-> out:" agrees with the manifest md5s (verbatim <=> shipped md5 == source md5; masked / reduced <=> they differ)
    _nlab, _ncb = 0, 0
    for _m in re.finditer(r"-> out:(\S+) \((verbatim copy|masked copy|reduced copy|copy;|code bundle; the shipped copy is byte-identical|code bundle; the shipped copy is derived|code bundle; the shipped copy is a copy;)", txt):
        _rp, _lab = norm(_m.group(1)), _m.group(2)
        if _lab.startswith("code bundle"):        # Stage 10 step 3b (S2-D68): the code-bundle label family is swept like the other three
            _ncb += 1
        _req(_rp in _MD5, "copy label on %s, which has no shipped manifest row" % _rp)
        _md5_, _smd5_ = _MD5[_rp]
        if _lab == "verbatim copy" or _lab.endswith("byte-identical"):
            _req(_smd5_ and _md5_ == _smd5_, "%s is labelled verbatim or byte-identical but its md5 differs from the source md5, or none is recorded" % _rp)
        elif _lab in ("masked copy", "reduced copy") or _lab.endswith("is derived"):
            _req(_smd5_ and _md5_ != _smd5_, "%s is labelled %s but its md5 equals the source md5, or none is recorded" % (_rp, _lab))
        else:
            _req(not _smd5_, "%s is labelled as having no source md5, but the manifest records one" % _rp)
        _nlab += 1
    _req(_nlab >= len(_L2), "the copy-label sweep matched %d line(s), fewer than the %d level-two inputs: the reader is blind" % (_nlab, len(_L2)))
    _req(_ncb >= 1, "the copy-label sweep matched no code-bundle label: the reader is blind to that family")
    # V11 assertion: the unit table's input cell equals the section-2 list of non-released inputs, unit by unit
    tab = {m.group(1): sorted(re.findall(r"`([^`]+)`", m.group(2))) for m in re.finditer(r"^\| (\S+) \| \S+ \| \S+ \| (.*) \|$", txt, re.M)}
    sec2 = txt.split("## 2. Inputs per unit", 1)[1].split("## 3.", 1)[0]
    parts = re.split(r"^### ", sec2, flags=re.M)
    for u in units:
        part = next((b for b in parts if b.startswith(u["unit"] + " (")), "")
        got = sorted(norm(m.group(1)) for m in re.finditer(r"^  - input: (.+?) -> (?!out:|release \(derived)", part, re.M))
        got = sorted(got + [norm(m.group(1)) for m in re.finditer(r"^  - input: (.+?) -> out:\S+ \(masked copy" + re.escape(STOP_CLAUSE), part, re.M)])      # Stage 10 step 3 fix round (iv): a released copy the generator stops on is in the section-1 cell too
        _req(tab.get(u["unit"]) == got, "REPRODUCE unit table cell != section 2 list for %s: %s vs %s" % (u["unit"], tab.get(u["unit"]), got))
        if u["kind"] == "figure" and u["class"] != "regenerates_from_release":
            _req(tab.get(u["unit"]), "%s keeps class %s but, with the ledger extract released, its list of inputs not in the release is empty" % (u["unit"], u["class"]))
    _doc_assertions("REPRODUCE.md", txt, out)
    return txt


def _masked_statement(out):
    """What the public ledger masks in THIS release, read from the shipped extract (RC1b fix Q09). A row that cites a file withheld only because it prints a TPMI line is not masked."""
    rows = read_tsv(Path(out) / LEDGER_PUBLIC_REL)
    reasons = sorted({json.loads(r["value"])["reason"] for r in rows if r["value"].startswith('{"masked":')})
    _req(reasons, "the public ledger masks no row, but the README describes masked rows")
    txt = "Masked rows in this release carry the reason %s. A row that cites a source file withheld only because it prints a TPMI line is not masked, because the derived TPMI values ship." % ", ".join("`%s`" % x for x in reasons)
    txt += " No row holds a raw per-variant TPMI value, so none is masked for `TPMI_RAW`." if "TPMI_RAW" not in reasons else " Rows with reason `TPMI_RAW` hold a raw per-variant TPMI value."
    return txt + "\n"


def _provenance_columns(out):
    """Non-plotted provenance columns, read from the shipped frames' own headers."""
    spec = [("source_data/figures/Figure_4/Figure_4__4a_forest.csv",
             [("ledger_p", "ledger row id of the p value"), ("printed_p", "the p text printed on the figure"),
              ("ref_ld_sensitive", "flag: drawn dotted because the estimate depends on the reference LD panel")]),
            ("source_data/figures/Figure_3/Figure_3__3b_forest.csv",
             [("quantity_ids", "ids of the quantities behind the row"), ("ledger_ids", "ledger row ids behind the row"), ("source_file", "file the row was read from"),
              ("printed_text", "the text printed on the figure")]),
            ("source_data/figures/Figure_3/Figure_3__3a_stored_results.csv", [("ledger_id", "ledger row id of the stored value"), ("source", "file the value was read from")]),
            ("source_data/figures/Figure_4/Figure_4__4b_credible_set_rows.csv", [("printed_text", "the text printed on the figure")])]
    res = []
    for f, cols in spec:
        p = Path(out) / f
        _req(p.is_file(), "provenance frame %s not in out/" % f)
        hdr = next(csv.reader(open(p, encoding="utf-8", newline="")))
        for c, _ in cols:
            _req(c in hdr, "%s has no column %s" % (f, c))
        res.append((f, cols))
    return res


# S2-D40 (defect 4): the successor pairs that the Stage 6 ledger batches created. ASSERTION ONLY: the note in the README is derived from the successors' claim text ("[successor of Lxxx: ...]").
SUPERSEDED_STAGE6 = {"L048": "L532", "L208": "L534", "L329": "L535", "L967": "L536", "L419": "L537", "L676": "L533"}


def _nrows(path):
    return int(shape(path)[0])


SCAN_SUFFIX = (".tsv", ".csv", ".md", ".txt", ".json", ".log", ".R", ".py", ".sh")
UNCONFIRMED_RX = re.compile("N " + "unconfirmed")                                   # pieces: the shipped copy of this file must not carry the phrase it scans for
LEGACY_LABEL = "BBJ_CAD_hum0197" + ".v5"                                          # pieces: the shipped copy of this file must not carry the label it scans for
LEGACY_LABEL_SPACED = "BBJ CAD (hum0197" + ".v5)"
LEGACY_LABEL_RX = re.compile(re.escape(LEGACY_LABEL) + "|" + re.escape(LEGACY_LABEL_SPACED))


def _files_with(out, rx):
    """S2-D41: every shipped text file of out/ whose text matches rx, sorted, found by scanning (never typed). The three self-written files are excluded (they describe the release)."""
    out = Path(out)
    hit = []
    for q in sorted(out.rglob("*")):
        rel = q.relative_to(out).as_posix()
        if q.is_file() and rel not in SELF_FILES and rel.endswith(SCAN_SUFFIX) and rx.search(q.read_text(encoding="utf-8", errors="replace")):
            hit.append(rel)
    return hit


def _legacy_labels_note(out):
    """README section 'Legacy labels' (S2-D41, text a): generated, not typed. The list of files is the result of scanning the shipped tree for the label."""
    files = _files_with(out, LEGACY_LABEL_RX)
    _req(files, "no shipped file carries the legacy BBJ CAD label any more: the README note is stale")
    tables = [f for f in files if f.endswith((".tsv", ".csv"))]
    code = [f for f in files if f not in tables]
    _req(tables, "the legacy BBJ CAD label is no longer in a shipped result table or frame")
    return ("## Legacy labels\n\n"
            "The label `%s` (also written `%s`) in the shipped result tables and frames is the run-time name of the PheWeb.jp CAD file "
            "(the `pheweb.jp` resource above: Kanai et al.; it is not a `hum0197` release). The release keeps the name so that a shipped table can be traced to the run output that carries it. "
            "Result tables and frames that carry it, found by scanning the shipped files: %s. Code that writes or tests the label: %s. "
            "The `hum0197.v3` labels (Biobank Japan T2D, ischemic stroke, CKD and HF) are correct and are unchanged.\n" % (
                LEGACY_LABEL, LEGACY_LABEL_SPACED, ", ".join("`%s`" % f for f in tables), ", ".join("`%s`" % f for f in code) if code else "none"))


def _superseded_note(out, pv):
    """README section 'Superseded ledger rows', derived from the public ledger (nothing typed): every row whose claim says '[successor of Lxxx: ...]' names its predecessor."""
    out = Path(out)
    with open(out / "ledger" / "claims_ledger_public.tsv", encoding="utf-8", newline="") as fh:
        pub = list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    by = {r["id"]: r for r in pub}
    pairs = {}
    for r in pub:
        m = re.search(r"\[successor of (L\d{3}):", r["claim"])
        if m:
            _req(m.group(1) in by, "ledger row %s names the predecessor %s, which is not a ledger row" % (r["id"], m.group(1)))
            _req(m.group(1) not in pairs, "ledger predecessor %s has two successors" % m.group(1))
            pairs[m.group(1)] = r["id"]
    for pred, succ in SUPERSEDED_STAGE6.items():
        _req(pairs.get(pred) == succ, "superseded pair %s -> %s is not what the ledger claims say (%s)" % (pred, succ, pairs.get(pred)))
    um = re.search(r"\(([^()]*unconfirmed[^()]*)\)", by["L967"]["claim"])
    _req(um is not None, "the claim of L967 no longer carries an 'unconfirmed' parenthetical")
    ident = [r["display_name"] for r in pv if r["provider_id"] == "EGFR_EAS_2P2M_LiuH2025"]
    _req(len(ident) == 1, "PROVIDERS.tsv has no single EGFR_EAS_2P2M_LiuH2025 row")
    unconf = _files_with(out, UNCONFIRMED_RX)
    _req("east_asian/ancestry_results.tsv" in unconf, "east_asian/ancestry_results.tsv no longer has a note that says the N is unconfirmed")
    pairs_txt = ", ".join("%s by %s" % (p_, pairs[p_]) for p_ in sorted(pairs))
    return ("## Superseded ledger rows\n\n"
            "`ledger/claims_ledger_public.tsv` keeps the predecessor of every row that a later ledger batch replaced. The claim of the successor names its predecessor as `[successor of Lxxx: ...]`; read the successor. "
            "Predecessor by successor: %s. "
            "Two kinds of text keep wording from before the East Asian eGFR file was identified: the claim of the predecessor %s (\"%s\") and the shipped files that carry the note that its N is unconfirmed, found by scanning the shipped tree: %s. "
            "All of them predate the identification of that file as `%s`; the successors %s and %s carry the identified file. The release does not rewrite any of these texts.\n" % (
                pairs_txt, "L967", um.group(1), ", ".join("`%s`" % f for f in unconf), ident[0], pairs["L329"], pairs["L967"]))


def _ms_title(v3):
    """The manuscript title, read the way scripts/P31R_build.py reads it: the source map's title_page component, and the same pattern on that file (Stage 10 step 3a, R1)."""
    mp = Path(v3) / "R2_rewrite" / "build" / "P31R_SOURCE_MAP.tsv"
    rows = [r for r in csv.DictReader(io.open(mp, encoding="utf-8", newline=""), delimiter="\t") if r["component"].strip() == "title_page"]
    _req(len(rows) == 1, "P31R_SOURCE_MAP.tsv has no single title_page row")
    tp = io.open(Path(v3) / "R2_rewrite" / rows[0]["canonical_file"].strip(), encoding="utf-8").read()
    m = re.search(r"\*\*Title[^*]*\*\*\s*(.+)", tp)
    _req(m is not None, "the title page has no **Title line")
    return m.group(1).strip()


def _readme(v3, out, led, code, prov):
    shipped_src = {s for b in _prior(v3) if b.get("ship_status") == "ship" and b.get("store", "out") == "out" for s in _srcs(b)}
    out = Path(out)
    P = led.P
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file() and p.relative_to(out).as_posix() not in SELF_FILES)
    tops = sorted({f.split("/")[0] for f in files})
    for t in tops:
        _req(t in FOLDERS, "folder %s has no description in rel_docs.FOLDERS" % t)
    for _, d in SOURCE_FOLDER:
        _req(any(f.startswith(d) for f in files), "Table 3 folder %s is not in out/" % d)
    pv = read_tsv(Path(v3) / PROV_REL)
    _req(pv and all(r["terms_status"] in ("confirmed_from_official_page", "inferred", "not_found") for r in pv), "PROVIDERS.tsv: a terms_status outside confirmed_from_official_page / inferred / not_found")
    _zr = json.loads((Path(v3) / ZENODO_RECORD_REL).read_text(encoding="utf-8"))      # G5 (2026-11-01): the DOI constants are checked against the Zenodo API record on every build
    _req(_zr.get("concept_doi") == CONCEPT_DOI, "rel_docs.CONCEPT_DOI differs from the Zenodo record %s" % ZENODO_RECORD_REL)
    _zraw = json.loads((Path(v3) / ZENODO_RAW_REL).read_text(encoding="utf-8"))
    _req(any(str(x.get("identifier", "")).startswith(REPO_URL + "/") for x in _zraw["metadata"].get("related_identifiers", [])), "rel_docs.REPO_URL is not the repository the Zenodo record names")
    L = []
    a = L.append
    a("# Release repository for \"%s\", version `%s`\n" % (_ms_title(v3), VERSION))
    a("Result tables, per-figure plotted frames, tables as data and analysis-tail code for the manuscript. `REPRODUCE.md` says what regenerates and what does not. "
      "`RELEASE_MANIFEST.tsv` lists every released file with its size, md5, row and column counts, source file and md5, producing script and md5, software label, citing document and licence class. "
      "The manifest cannot hold its own md5, so it and the files `README.md`, `REPRODUCE.md`, `LICENSE`, `CITATION.cff` and `.gitattributes` have no manifest row. `.gitattributes` holds the one line `* -text`, so that git never rewrites line endings and a clone equals the zip byte for byte.\n")
    a("## Version, licence and citation\n")
    a("This is release `%s`. `LICENSE` gives the terms: the MIT License for code the authors wrote, and `CC BY 4.0` for tables the authors derived "
      "(`own_derived` in `RELEASE_MANIFEST.tsv`); files derived from third-party data keep their providers' terms (below). "
      "`CITATION.cff` says how to cite the release. The repository is `%s`; its archived versions are on Zenodo under the concept DOI `https://doi.org/%s`, which resolves to the latest version.\n" % (VERSION, REPO_URL, CONCEPT_DOI))
    _zroot = "P31R_R2_release_%s/" % VERSION
    _paths = [str(f).replace("\\", "/") for f in files]      # derived from the built tree, never typed; printed as code spans (R7 blanks code spans)
    _longest = sorted(f for f in _paths if len(f) == max(len(x) for x in _paths))      # Stage 10 step 3 fix round: every path of the maximal length, so a tie is reported, not hidden by max()
    _req(_longest, "the release holds no files")
    a("On Windows, extracting the zip or cloning the repository into a deep folder can fail once a full path exceeds the default Windows path limit (`MAX_PATH`). "
      "The longest path%s inside this release %s %s (for a clone, the repository folder name replaces the zip root `%s`), so extract or clone close to the drive root, or run `git config core.longpaths true`.\n"
      % ("s" if len(_longest) > 1 else "", "are" if len(_longest) > 1 else "is", " and ".join("`%s%s`" % (_zroot, f) for f in _longest), _zroot))
    a("## Licence tiers, third-party data and the TPMI switch\n")
    a("The release ships the authors' derived results and code. Third-party data are pointed to (provider, accession or URL, citation), not redistributed. Each file has a tier, decided by its content:\n")
    a("- First tier, derived per-pair results (MR estimates, colocalization posteriors, counts, classes): ship, with a citation.")
    a("- Second tier, per-variant values copied from a provider: ship only where redistribution is confirmed (CC0 or CC BY); otherwise the column or file is withheld, or reduced to the values plotted (regional tracks are binned). The exceptions that remain are listed in the paragraph after these tiers.")
    import rel_ledger as _rl
    _tpmi = _rl.tpmi_switch(v3)
    a("- Third tier, no licence (TPMI): the derived per-pair values ship with a citation; the TPMI effect, standard error and allele frequency columns do not.\n" if _tpmi == "ship" else "- Third tier, no licence (TPMI): withheld.\n")
    a(PER_VARIANT_NOTE + "\n")
    if _tpmi == "ship":
        a("TPMI is one declared switch (`release_records/TIER_RULES.tsv`, row `SW1`), set to `ship` for this release. The derived per-pair TPMI values ship with a citation of the TPMI PheWeb (`PheCode 250.2`): the Wald ratio, its standard error, p, class and verdict, and the study-level case and control counts. "
          "The TPMI effect, standard error and allele frequency columns of the check table are withheld. A single-instrument Wald ratio with its standard error is a derived per-pair result (first tier) and ships. "
          "ChEMBL- and DGIdb-derived drug and stage data fields are withheld (ChEMBL is CC BY-SA); the mask covers these data fields only, and drug names in manuscript prose come from the literature. Shipped scripts that quote ChEMBL facts as string literals are declared under Code licence. Open Targets (CC0) fields are native and ship.\n")
    else:
        a("TPMI is one declared switch, here `withheld`: the TPMI-derived ledger rows are masked and the TPMI-derived files are kept out. "
          "ChEMBL- and DGIdb-derived drug data fields are withheld (ChEMBL is CC BY-SA); the mask covers these data fields only, and drug names in manuscript prose come from the literature. Shipped scripts that quote ChEMBL facts as string literals are declared under Code licence. Open Targets (CC0) fields are native.\n")
    a("Where each third-party source is obtained, and how to cite it. A source whose redistribution is not confirmed contributes derived results; the cases in which an individual per-variant value appears are listed above. "
      "A licence is given only where the provider states one. The page that was read for each provider, the wording it states and its status are recorded in `release_records/PROVIDER_TERMS.tsv`. "
      "Names, routes and citations are identifiers and are set in code font.\n")
    a("The provider list holds %s.\n" % P("L533"))      # the live list (L533); the frozen 1c row L676 is never cited
    a("| provider | accession or URL | citation |")
    a("|---|---|---|")
    pt = {t["provider_id"].strip(): t for t in read_tsv(Path(v3) / TERMS_REL)}
    cell = lambda x: "`%s`" % ("not stated" if (x or "").strip() in ("", "Unknown", "Not fetched") else x).replace("|", "/").replace("`", "'").replace(chr(10), " ").strip()
    stated = re.compile(r"(?i)^(CC0|CC[ -]BY)")        # a licence name an official page states; 'not named', 'none named' and 'not_found' are not a licence
    for r in pv:
        t = pt.get(r["provider_id"], {})
        route = (t.get("distributor_or_access_route") or "").strip()
        accs = sorted(set(ACC_RX.findall(t.get("exact_release_or_accession_used") or "")))      # accession ids only; sample sizes and file lists stay in PROVIDER_TERMS.tsv
        where = "; ".join(x for x in ([route] if route not in ("", "Unknown", "Not fetched") else []) + [x for x in accs if x not in route])
        lic = (t.get("licence_name") or "").strip()
        cite = (t.get("original_study_or_provider") or "").strip()
        if t.get("status", "").strip() == "confirmed_from_official_page" and stated.match(lic):
            cite += "; licence " + lic.split(" (")[0]
        a("| %s | %s | %s |" % (cell(r["display_name"]), cell(where), cell(cite)))
    a("")
    a("## Folders\n")
    for t in tops:
        fl = [f for f in files if f.split("/")[0] == t]
        a("- `%s/`: %s." % (t, FOLDERS[t]))
        if t != "code":
            a("  Files: %s." % ", ".join("`%s`" % f[len(t) + 1:] for f in fl))
    a("\nWorkstream folders, files released: %s. Branch table: %s.\n" % (P("L670"), P("L672")))
    a("## Table 3 source ids and release folders\n")
    a("The repository-row map of Table 3 (below) and `branches_integrated.tsv` name rows by source id. The ids are not all folder names:\n")
    a("| source id | release folder |")
    a("|---|---|")
    for s, d in SOURCE_FOLDER:
        a("| %s | `%s` |" % (s, d))
    a("\n## Repository rows behind the claims of Table 3\n")
    a("Table 3 of the manuscript prints the values; this list names the result rows that hold them, as code, by folder and row-identifier prefix "
      "(rows of workstreams/BI by workstream and row number).\n")
    a("| claim | repository rows |")
    a("|---|---|")
    _t3 = list(csv.DictReader(io.open(Path(v3) / T3MAP_REL, encoding="utf-8", newline=""), delimiter="\t"))
    _req(len(_t3) >= 1 and set(_t3[0]) == {"claim", "repository_rows"}, "Table 3 repository map is empty or malformed")
    for _r in _t3:
        a("| `%s` | `%s` |" % (_r["claim"].replace("|", "/").replace("`", "'"), _r["repository_rows"].replace("|", "/").replace("`", "'")))
    a("\n## Tables\n")
    a("Rows of Table 1: %s. Rows of Table 2: %s. Table 2 ships whole except the ChEMBL-derived drug and development-stage cells (a stage set by a drug record, the agent list, and statements or verdicts that name a drug), which read `withheld_pending_terms`; the Open Targets tractability class and the derived TPMI values of the external check ship. `tables/Table2_evidence_cards.columns.tsv` defines every column. The ChEMBL facts behind the masked cells are quoted, under the licence of ChEMBL, in shipped scripts, as declared under Code licence.\n" % (P("L673"), P("L671")))
    a("## The public claims ledger\n")
    a("`ledger/claims_ledger_public.tsv` carries the id, section, claim, value, display, source, evidence_status and source_md5_before_release_transform of every ledger row. Internal notes and the verification history are dropped. "
      "A masked value is the JSON object `{\"masked\":true,\"reason\":...}`, never zero or blank, and its display reads `[masked: <reason>]`. "
      "A row whose value is a list of word statements (the Table 2 cards) stays unmasked except for the development-stage clauses that carry ChEMBL vocabulary (a stage word, such as approved or in clinical development, with the agents it names): "
      "each such clause is replaced inside its list element by `[masked: CHEMBL]`, the clauses before it (Open Targets-native, for example tractability) stay, and the display is the same list with the same replacement. "
      "`source` names the file, line and column a value was read from and the md5 of that file AS SHIPPED in this release, with the line numbers of the shipped file. "
      "Where the release rewrote the file (masking, reduction, portability rewrite or line-break replacement), `source_md5_before_release_transform` holds the md5 of the file before the rewrite, separated by `;` when a row cites several files; "
      "for a cited file that is not in the release both are the original md5, and the column is empty for a row that cites no md5. `evidence_status` is one of:\n")
    import rel_ledger
    for k_, v_ in rel_ledger.EVIDENCE_STATUS.items():
        a("- `%s`: %s." % (k_, v_))
    a("")
    a(_masked_statement(out))
    a(_superseded_note(out, pv))
    a(_legacy_labels_note(out))
    a("## Column notes: the p values of the renal files\n")
    a("- `kidney/egfr_mirage_classification.csv`, columns `crea_p`, `cys_p` and `bun_p`: the two-sided normal p of a single-instrument Wald ratio, `2 * pnorm(-abs(b / se))`, recomputed from the effect and standard error that the GWAS prints to four decimal places (script `code/paper3_ckd_rerun_v3/scripts_sent/01_renal_triad_mr.R`). "
      "`class` is concordant, mirage or other, from the cystatin-C result and the sign (the script's header defines the three).")
    a("- `source_data/figures/FigureS3/FigureS3__renal_proteins.csv`, columns `cys_p` and `bun_p`: the CKDGen GWAS's own P value for the instrument SNP; for a single SNP it is the first-order Wald p. They set the filled and hollow markers of FigureS3 panels b and c (filled: nominally significant). "
      "`egfrcrea_p` is the shipped single-instrument eGFR Wald p. The two definitions agree on every call of nominal significance and differ in magnitude (the recomputation from four-decimal values is coarser); Supplementary Section S-R12 documents both.\n")
    a("## Provenance columns in the plotted frames\n")
    a("Some frame columns are not plotted. They tie each plotted value to its ledger row or printed text:\n")
    for f, cols in prov:
        for c, d in cols:
            a("- `%s`, column `%s`: %s." % (f, c, d))
    a("")
    a("## Withheld and not shipped\n")
    a("Files whose provider terms are unread, and files under source hold, are kept out of this repository. Individual-level dosage files and third-party summary-statistic windows are not released; `REPRODUCE.md` names them per unit.\n")
    a("## Code bundle\n")
    a("`code/` holds the shipped files of `release_records/CODE_INCLUDE.tsv` at their project-relative paths, and the records the release builder wrote into it (the input map, the portability log and the reduced inputs of the level-two figures). Declared but not shipped (pending a decision): rows of that file with ship = 0, of role %s.%s The bundle files:\n" % (", ".join("`%s`" % x for x in sorted({r["role"] for r in code if r["ship"] == "0" and norm(r["path"]) not in shipped_src})),
                                                 "".join(" Declared ship = 0 but not pending, because the file ships elsewhere in the release: role `%s` (%s)." % (r["role"], r["note"].strip().rstrip(".")) for r in code if r["ship"] == "0" and norm(r["path"]) in shipped_src)))
    for r in code:
        if r["ship"] == "1":
            a("- `%s`" % r["ship_path"])
    shipped_paths = {r["ship_path"] for r in code if r["ship"] == "1"}
    code_files = [f for f in files if f.startswith("code/")]
    extra = sorted(f for f in code_files if f not in shipped_paths)
    a("")
    _req(sorted(code_files) == sorted(shipped_paths | set(extra)), "README code listing differs from the files of code/")      # the listing is derived from code/ itself, so no count is typed
    a("The release builder also wrote these files into `code/` (the input map, the portability log and the reduced inputs of the level-two figures):\n")
    for f in extra:
        a("- `%s`" % f)
    a("")
    a("## Provenance labels of the engine and deposit scripts\n")
    a("A script that came with a deposit package from the analysis machine (the v4 and v5 deliveries) carries one of three labels, in the `provenance_label` column of `RELEASE_MANIFEST.tsv` and of `release_records/CODE_INCLUDE.tsv`. "
      "A label describes the ORIGINAL file as the deposit holds it. It does not say that the copy in `code/` is that file byte for byte. Where the portability transform changed the copy, "
      "`source_md5` in the manifest is the original's md5, `md5` is the shipped copy's, and `transform` lists the declared changes; the copy is derived from the matched original and is not byte-identical to it. "
      "`evidence_md5` and `provenance_evidence` in `CODE_INCLUDE.tsv` give the recorded hash (or prefix) and where it was recorded, and the build compares them with the original.\n")
    a("| label | what the evidence is | what it does not establish |")
    a("|---|---|---|")
    for lab_, ev_, no_ in PROV_TABLE:
        a("| `%s` | %s | %s |" % (lab_, ev_, no_))
    a("")
    a("A script that evaluates sections of other scripts (the candidate SuSiE script evaluates the function blocks of the core and eGFR SuSiE engines; the sensitivity script patches their text in memory) ships with those parents. "
      "A match for one file does not verify the files it evaluates, or the settings and inputs of the run.\n")
    for lab_, _, _ in PROV_TABLE:
        grp = [r for r in code if r["ship"] == "1" and r.get("provenance_label") == lab_]
        a("Files labelled `%s`:\n" % lab_)
        for r in grp:
            a("- `%s`" % r["ship_path"])
        a("")
    a("## Code licence\n")
    a("The MIT licence (`LICENSE`) applies to all shipped code, which the authors wrote, except a file flagged `third_party:<what>`. The `code_rights` column of `release_records/CODE_INCLUDE.tsv` says `authors` for the engine and deposit scripts that carry a provenance label (an empty `code_rights` does not mark a file as third-party), "
      "on the basis of a scan of each file for third-party copyright or licence notices (none was found). "
      "A file that contains third-party code is flagged `third_party:<what>`, keeps its own licence and is not relabelled. "
      "The R packages these scripts call and external programs such as PLINK are not shipped. A code licence does not cover provider values quoted inside scripts; those follow the cases above and `release_records/LEAKSCAN_COINCIDENCES.tsv`. "
      "Some shipped scripts quote ChEMBL facts as string literals: agent names, development-stage statements, approval status and ChEMBL identifiers, chiefly those behind the masked cells of `tables/Table2_evidence_cards.tsv`. Those facts are reproduced under the Creative Commons Attribution-ShareAlike licence of ChEMBL, with attribution to ChEMBL (EMBL-EBI) and the licence version as recorded in `release_records/PROVIDER_TERMS.tsv`; the MIT grant does not cover them, and a work that redistributes them inherits the share-alike term.\n")
    return "\n".join(L)


def _build(v3, out, withheld, ledger):
    v3, out = Path(v3), Path(out)
    led = Ledger(ledger)
    code = _code_rows(v3)
    CODE_SHIPPED.clear()
    CODE_SHIPPED.update({norm(r["path"]): r["ship_path"] for r in code if r["ship"] == "1"})
    tables = _regen(v3)
    prov = _provenance_columns(out)
    rows = _manifest_rows(v3, out, code)
    _write_tsv(out / "RELEASE_MANIFEST.tsv", MANIFEST_COLS, rows)
    _write_text(out / "REPRODUCE.md", _reproduce(v3, out, led, code, tables))
    _write_text(out / "README.md", _readme(v3, out, led, code, prov))
    _doc_assertions("README.md", (out / "README.md").read_text(encoding="utf-8"), out)
    _lic = LICENSE_TEXT.replace("{TITLE}", _ms_title(v3))      # Stage 10 step 3 fix round: the internal name no longer ships; the title is the one the README prints
    _req("{" not in _lic.splitlines()[0] and _ms_title(v3) in _lic, "LICENSE title line not filled")
    _write_text(out / "LICENSE", _lic)
    _write_text(out / "CITATION.cff", CITATION_TEXT)
    _write_text(out / ".gitattributes", GITATTRIBUTES_TEXT)
    for r in code:         # a row declared as not shipped must not be in the release (README says so; RC1 fix K03)
        if r["ship"] == "0":
            _req(not (out / "code" / r["path"]).is_file(), "CODE_INCLUDE says %s is not shipped, but code/%s exists" % (r["path"], r["path"]))
    # every file under out/ except the three has exactly one manifest row; every row's file exists
    on_disk = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()} - set(SELF_FILES)
    _req(on_disk == {r["path"] for r in rows}, "manifest and out/ disagree")
    common = [LEDGER_REL, CODE_REL, PROV_REL, TERMS_REL]
    return [
        {"release_path": "RELEASE_MANIFEST.tsv", "ship_status": "ship", "source_paths": common, "producer": PRODUCER, "licence_class": "own_derived",
         "rows": len(rows), "columns": len(MANIFEST_COLS),
         "note": "One row per released file except itself, README.md and REPRODUCE.md (a file cannot carry its own md5)."},
        {"release_path": "REPRODUCE.md", "ship_status": "ship", "source_paths": common, "producer": PRODUCER, "licence_class": "own_derived", "rows": "", "columns": "",
         "note": "Regenerability per unit, what does not regenerate, provenance caveats, environment, commands. Every count is a ledger row display with its id."},
        {"release_path": "README.md", "ship_status": "ship", "source_paths": common, "producer": PRODUCER, "licence_class": "own_derived", "rows": "", "columns": "",
         "note": "Folder map, Table 3 source id -> folder map, licence statement, provenance columns of the plotted frames, code bundle."},
        {"release_path": "LICENSE", "ship_status": "ship", "source_paths": common, "producer": PRODUCER, "licence_class": "own_derived", "rows": "", "columns": "",
         "note": "MIT for the authors' code, CC BY 4.0 for the authors' own-derived tables (S2-D27(1)); third-party-derived files keep their providers' terms."},
        {"release_path": "CITATION.cff", "ship_status": "ship", "source_paths": common, "producer": PRODUCER, "licence_class": "own_derived", "rows": "", "columns": "",
         "note": "Citation File Format 1.2.0; version %s; all authors listed; no doi key (Zenodo mints a version DOI per GitHub release; the concept DOI is in the README)." % VERSION},
        {"release_path": ".gitattributes", "ship_status": "ship", "source_paths": common, "producer": PRODUCER, "licence_class": "own_derived", "rows": "", "columns": "",
         "note": "The single line `* -text`: git never rewrites line endings, so a clone equals the zip byte for byte."},
    ]


def build(v3, out, withheld, ledger):
    """Spec section 2 interface. Writes RELEASE_MANIFEST.tsv, REPRODUCE.md and README.md at the root of out/ (out/ must already hold every other released file)."""
    return _build(v3, out, withheld, ledger)


# ----------------------------------------------------------------------------------------------------------------------------------------------
def _default_ledger(v3):
    import pandas as pd
    return pd.read_csv(Path(v3) / LEDGER_REL, dtype=str)


def _populate(v3, out):
    """Standalone test helper: copy the real out/ (read-only) and the ship=1 code rows into a scratch out/."""
    import shutil
    real = Path(v3) / REL_DIR / "out"
    for p in real.rglob("*"):
        if p.is_file() and p.relative_to(real).as_posix() not in SELF_FILES:
            d = Path(out) / p.relative_to(real)
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p, d)
    for r in _code_rows(v3):
        d = Path(out) / r["ship_path"]
        if r["ship"] == "1" and (Path(v3) / r["path"]).is_file() and not d.exists():     # the real out/ already holds the (portability-transformed) bundle copy; the source is only a fallback
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(Path(v3) / r["path"], d)


def _v3():
    return Path(__file__).resolve().parents[2]


def selftest():
    import shutil
    import tempfile
    v3 = _v3()
    led = _default_ledger(v3)
    root = Path(tempfile.mkdtemp(prefix="rel_docs_selftest_"))
    bad = []

    def rep(name, ok, detail=""):
        print("  %-60s %s %s" % (name, "ok" if ok else "MISS", detail))
        if not ok:
            bad.append(name)

    out = root / "out"
    _populate(v3, out)
    recs = build(v3, out, root / "wh", led)
    first = {f: (out / f).read_bytes() for f in SELF_FILES}
    build(v3, out, root / "wh", led)
    rep("two builds are byte-identical", all((out / f).read_bytes() == first[f] for f in SELF_FILES))
    m = read_tsv(out / "RELEASE_MANIFEST.tsv")
    rep("manifest md5 and bytes match the files", all(md5_file(out / r["path"]) == r["md5"] and str(os.path.getsize(out / r["path"])) == r["bytes"] for r in m))
    rd = (out / "README.md").read_text(encoding="utf-8")
    rp = (out / "REPRODUCE.md").read_text(encoding="utf-8")
    rep("ledger tokens are followed by their row id", all(LEDGER_PRINT_RX.search(x) for x in (rd, rp)))

    def trips(name, fn):
        try:
            fn()
        except DocsError as e:
            rep("mutation: %s" % name, True, "(%s)" % str(e)[:70])
            return
        rep("mutation: %s" % name, False, "(no error)")

    d1 = root / "m1"
    shutil.copytree(out, d1)
    (d1 / "stray.txt").write_text("x")
    trips("unrecorded file in out/", lambda: build(v3, d1, root / "wh", led))
    d2 = root / "m2"
    shutil.copytree(out, d2)
    t = next(r for r in _code_rows(v3) if r["ship"] == "1")
    (d2 / t["ship_path"]).write_bytes(b"changed")
    trips("code copy differs from its source", lambda: build(v3, d2, root / "wh", led))
    d3 = root / "m3"
    shutil.copytree(out, d3)
    os.remove(d3 / t["ship_path"])
    trips("bundle file missing from out/", lambda: build(v3, d3, root / "wh", led))
    import pandas as pd
    led2 = led.copy()
    i = led2.index[led2["id"] == "L664"]
    _req(len(i) == 1, "selftest target L664 not unique")
    led2.loc[i, "display"] = ""
    trips("ledger row without display", lambda: build(v3, out, root / "wh", led2))
    rep("selftest scratch kept at", True, str(root))
    print("rel_docs selftest: %s" % ("FAIL " + str(bad) if bad else "all ok"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))     # python -I leaves the script folder off sys.path; the sibling modules (rel_ledger, rel_portability) are imported lazily
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    if "--build" in sys.argv:
        i = sys.argv.index("--build")
        o, w = sys.argv[i + 1], sys.argv[i + 2]
        r = build(_v3(), o, w, _default_ledger(_v3()))
        for x in r:
            print(x["release_path"], x["rows"], x["columns"])
        sys.exit(0)
    print(__doc__)
