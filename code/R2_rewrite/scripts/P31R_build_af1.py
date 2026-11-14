#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 rewrite (Stage 3b, track D) -- Additional file 1, the supplementary-table workbook.

Every sheet is DERIVED here from its source file; nothing is typed. Each sheet's counts are asserted against the claims
ledger (V3/R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv), and every numeral in a title or legend is read from
the ledger row it names. The vocabulary sheet (writing frame section 3) is enforced on every text cell with the prose
gate's own banned list (P31R_gate.BANNED, imported, not copied). Only the README columns that hold file paths are exempt.

Sheets: README; S1 STROBE-MR checklist; S2 data sources (from tables/P31R_TableS2.md); S3 the colocalization-supported
pairs; S4 development stage; S5 external outcome checks; S6 multi-signal candidates; S7a/S7b East Asian per-pair results;
S8a-S8c Supplementary Note tables.

Table S1 reads the manuscript through the component source map (build/P31R_SOURCE_MAP.tsv, via P31R_build.MAP / mread), never from a
directory listing. Every location it prints is a section, the subsection heading as printed, a paragraph number within it and an
anchor phrase that must occur exactly once in that paragraph, in the markdown split AND in the built .docx (Stage 6.3b part 2b).

Run:  python P31R_build_af1.py            -> $P31R_BUILD_OUT (default R2_rewrite/build)/P31R_Additional_file_1.xlsx
      P31R_AF1_DOCX=<built manuscript>    -> the .docx whose paragraphs Table S1 is checked against
                                             (default <OUT>/P31R_Manuscript_R2_preview.docx; used when newer than every map component)
      P31R_AF1_FINAL=1 python ...         -> packaging mode: refuses to write while any sheet or the README shows PLACEHOLDER or a
                                             [... TO CONFIRM] marker, while no current built manuscript checks Table S1, or while
                                             the stale Stage 2c skeleton outputs sit beside the workbook
Exit 1 on any failed assertion. Mutation ladder: P31R_af1_mutation_ladder.py; independent recompute of L450: P31R_af1_verify_devstage.py.
"""
import hashlib, io, json, os, re, sys, datetime
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V3 = os.path.dirname(ROOT)
CKM = os.path.dirname(V3)
S1 = os.path.join(V3, "R1_stage1")
OUT = os.environ.get("P31R_BUILD_OUT") or os.path.join(ROOT, "build")
os.makedirs(OUT, exist_ok=True)
AUTHOR = "Bertrand Chin-Ming Tan"
sys.path.insert(0, HERE)
import P31R_gate as PG  # noqa: E402  (ledger reader, tag pattern, banned vocabulary: one implementation)
import P31R_af1_devstage as DEV  # noqa: E402  (development stage: shared with the ledger extension ext_af1.py)
import P31R_build as MB  # noqa: E402  (the component source map: MB.MAP, MB.mread; importing it has no side effects)

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ------------------------------------------------------------------------------------------------ sources (one place)
SRC = dict(
    tiered=os.path.join(V3, "v3.2/_engine_root/results/phase4_annotation/table_s_pair_evidence_tiered.csv"),
    mr=os.path.join(V3, "V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv"),
    ins=os.path.join(V3, "V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv"),
    vep=os.path.join(V3, "V4_local/source_data/annotation/lead_snp_consequences.csv"),
    master=DEV.MASTER, audit=DEV.AUDIT, oldt4=DEV.OLDT4,
    strobe=os.path.join(ROOT, "data", "P31R_strobe_mr_items.tsv"),
    proxy=os.path.join(S1, "S1R1_R2_08_proxy_record", "proxy_record.tsv"),   # S2-D47 (Stage 7 A2.2): the five proxy pairs; written by w8_proxy_record.py
    tables2=os.path.join(ROOT, "tables", "P31R_TableS2.md"),
    oc=os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_outcome_checks.tsv"),
    ph=os.path.join(ROOT, "posthoc/P31R_posthoc_results.tsv"),
    ckbf=os.path.join(CKM, "paper 3.eas/analysis/data/eas_exposure/FGFR4_cis_eas.tsv"),   # CKB cis-pQTL at FGFR4: source of the FGFR4 discovery allele and frequency in S5 (md5 pinned below)
    harm=os.path.join(V3, "V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz"),     # ARG1 6q23 harmonised panel: source of the signal B allele and frequency in S5 (md5 pinned below)
    cand=os.path.join(V3, "V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv"),
    eas=os.path.join(V3, "V4/source_data/supp_tables/SuppTable_EAS_axis.csv"),
    cen=os.path.join(ROOT, "data", "P31R_eas_census_results.tsv"),      # relocated from V4_local/staging/d2-eas-census (md5 pinned below)
    ckf=os.path.join(ROOT, "data", "P31R_eas_census_ckbfree.tsv"),
    snap=DEV.SNAPSHOT,
    topa=os.path.join(S1, "S1R1_02_topology_note/S1R1_topo_supp_table_a_intercepts.tsv"),
    topb=os.path.join(S1, "S1R1_02_topology_note/S1R1_topo_supp_table_b_pairs.tsv"),
    topc=os.path.join(S1, "S1R1_02_topology_note/S1R1_topo_supp_table_c_excess.tsv"),
    ledger=PG.LEDGER)
# Table S1 locates its items in the manuscript components of the source map (and in the R2 Results Detail, for Additional file 2 sections);
# those files are its sources and are listed in the README with their md5
S1_COMP = [("abstract", "Abstract"), ("background", "Background"), ("methods", "Methods"), ("results", "Results"), ("discussion", "Discussion"),
           ("conclusions", "Conclusions"), ("declarations", "Declarations")]
SRC.update({"ms_" + c: MB.msrc(c) for c, _ in S1_COMP})
SRC["ms_af2_results_detail"] = MB.msrc("af2_results_detail")
S1_SRCS = ["strobe"] + ["ms_" + c for c, _ in S1_COMP] + ["ms_af2_results_detail"]


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


# the ancestry-matched census results were relocated from a run-staging folder into R2_rewrite/data (P31R_DATA_PROVENANCE.md):
# pin their identity, and while the staging originals exist assert they have not changed since (a re-run must be re-relocated)
CENSUS_MD5 = {"cen": "9a62d288ee4919cc6434eaaa1c62a7bc", "ckf": "2ae5f2f93d94467881c417fd716e174a"}
for _k, _orig in (("cen", "results_census.tsv"), ("ckf", "results_ckbfree.tsv")):
    assert md5(SRC[_k]) == CENSUS_MD5[_k], ("relocated census file changed", SRC[_k])
    _o = os.path.join(V3, "V4_local/staging/d2-eas-census", _orig)
    assert not os.path.exists(_o) or md5(_o) == CENSUS_MD5[_k], ("the staging original differs from the relocated copy: re-relocate", _o)
# the proxy record (Stage 7 A2.2, S2-D47) is a source: its identity is pinned here and the ledger row that counts it (L543) must name the same md5 (checked once the ledger is read)
PROXY_MD5 = "ef2af6e512db8981924edd72d1efe724"
assert md5(SRC["proxy"]) == PROXY_MD5, "proxy_record.tsv changed since this builder was written: regenerate with w8_proxy_record.py, rebuild the ledger, then update the pin (do not re-pin silently)"


# README "Source file(s)": the REPOSITORY path of each source, never a local project path (S2-D62). It is read from the release's own records:
# release/out/RELEASE_MANIFEST.tsv (each shipped file, the project paths it was made from, its md5) and release/inputs/NOT_SHIPPED.tsv (files kept out,
# with the reason). A source the repository does not hold is described in words.
def _release_table(p):
    assert os.path.isfile(p), ("release record missing: the README prints repository paths from it (restore or rebuild the release first)", p)
    return pd.read_csv(p, sep="\t", dtype=str).fillna("")


REL_MAN = _release_table(os.path.join(ROOT, "release", "out", "RELEASE_MANIFEST.tsv"))
REL_NOTSHIP = _release_table(os.path.join(ROOT, "release", "inputs", "NOT_SHIPPED.tsv"))
REL_STALE = []   # (source, shipped path): the manifest pins another md5 for the single source than the file read here, so the release must be rebuilt first


def repo_src(p):
    v3rel = os.path.relpath(p, V3).replace("\\", "/")
    here = md5(p)
    # G5 (2026-11-01, release-F1): a manifest row built from a relocated copy of this file names the copy's path in source_path (V4_local/staging/...), so a
    # file with no path match is matched by its source_md5 (path first, so byte-identical vendored twins keep their own row); the two Table S7b census
    # sources (results_census.tsv, results_ckbfree.tsv) are the case in point
    hits = REL_MAN[REL_MAN["source_path"].apply(lambda sp: v3rel in sp.split(";"))]
    if not len(hits):
        hits = REL_MAN[REL_MAN["source_md5"].apply(lambda sm: here in sm.split(";"))]
    if len(hits):
        out = []
        for _, h in hits.iterrows():
            if ";" not in h["source_path"] and h["source_md5"] != here:
                REL_STALE.append((v3rel, h["path"]))
            out.append(h["path"] if h["md5"] == here else h["path"] + " (derived from this file)")
        return "; ".join(out)
    ns = REL_NOTSHIP[REL_NOTSHIP["path"] == v3rel]
    if len(ns) and ns["reason_code"].iloc[0].startswith("THIRD_PARTY"):
        return "third-party file, not redistributed (providers and terms: repository README)"
    if v3rel.startswith("R2_rewrite/sections"):
        return "manuscript source text of this package, not in the repository"
    return "project working file, not in the repository (no entry in RELEASE_MANIFEST.tsv)"


L = PG.ledger()


def LV(i):
    return json.loads(L[i]["value"])


def spell(n):
    """One to nine in words (house rule), the numeral otherwise; the number itself is always read from the ledger."""
    return "zero one two three four five six seven eight nine".split()[n] if 0 <= n <= 9 else str(n)


def eq(a, b, what):
    assert a == b, "%s: workbook %r != ledger %r" % (what, a, b)


def dash(s):
    return s.replace("-", "–")


def tf(x):
    return str(x).strip().upper() == "TRUE"


def yn(x):
    return "yes" if x else "no"


def mapv(series, d, what):
    miss = set(series.dropna().unique()) - set(d)
    assert not miss, "%s: unmapped values %s" % (what, miss)
    return series.map(lambda v: d[v] if v in d else v)


SHEETS = []   # (name, title, legend, DataFrame, {column: description}, [source keys], {column: number format})


def sheet(name, title, legend, df, cols, srcs, fmt=None):
    assert list(df.columns) == list(cols), (name, [c for c in df.columns if c not in cols], [c for c in cols if c not in df.columns])
    assert not df.isna().any().any(), "%s: empty cell(s) in %s" % (name, list(df.columns[df.isna().any()]))
    SHEETS.append((name, title, legend, df, cols, srcs, fmt or {}))


N81, N62 = LV("L010"), LV("L011")
NPX = LV("L543")[0]                                              # pairs among the 81 whose outcome association was read at an LD proxy (L543, Stage 7 A2.2)
assert PROXY_MD5 in L["L543"]["source"], "ledger row L543 does not name the pinned proxy_record.tsv (stale ledger or stale record)"
P_SCREEN = LV("L238")
PRIOR, PRIOR6, RULE_H4, RULE_RATIO = LV("L240")[2], LV("L241"), LV("L236"), LV("L237")


def e_(x):
    m_, ex = ("%.0e" % x).split("e")
    return "%se%d" % (m_, int(ex))


# every threshold printed in a label is read from the ledger row that holds it (numbers-vs-source), never typed
T80 = "%.2f" % RULE_H4                                           # colocalization rule PP.H4 (L236)
HALF = LV("L962")["prior_weak_below"] if "L962" in L else 0.5   # weakness threshold (supplement row L962); asserted against the tier table's flag
F10 = "%g" % LV("L235")                                          # minimum F (L235)
FDRP = "%g%%" % (100 * LV("L239"))                               # FDR level (L239)
PNOM = "%g" % LV("L238")                                         # screen threshold (L238)
E5, E6, GW = e_(PRIOR), e_(PRIOR6), e_(LV("L231"))               # priors (L240, L241); genome-wide threshold (L231)
NFAM = "{:,}".format(LV("L003"))                                 # multiple-testing family (L003)
BON = "%g/%d" % (LV("L961")[0], LV("L268"))                      # Bonferroni for the checks selected in advance (L961, L268)
PUR = "%g" % LV("L247")                                          # SuSiE purity (L247)
C_H4D, C_H4C = "PP.H4, default prior (p12 = %s)" % E5, "PP.H4, conservative prior (p12 = %s)" % E6
C_GW, C_BON, C_PUR = "Instrument p < %s" % GW, "Passes %s" % BON, "SuSiE PP.H4, purity %s" % PUR
AUDIT_DATE = L["L500"]["display"]                                # date of the Open Targets audit (L500)


def sci(x):
    m_, e_ = ("%.0e" % x).split("e")
    return "10⁻%s" % "".join("⁰¹²³⁴⁵⁶⁷⁸⁹"[int(c)] for c in str(abs(int(e_)))) if m_ == "1" else "%s × 10⁻%s" % (m_, "".join("⁰¹²³⁴⁵⁶⁷⁸⁹"[int(c)] for c in str(abs(int(e_)))))


# ================================================================================================ S1 STROBE-MR
def section_paragraphs():
    """{(file stem, subsection or ''): [paragraph text, ...]} with ledger tags stripped, split exactly as the manuscript
    builder's emit() does (P31R_build.py; nested there, so mirrored here): blocks separated by a blank line; lines inside a
    block are joined into ONE paragraph; a '## ' block is a heading whose remaining lines are re-emitted; a '### ' block is a
    subheading whose remaining lines form one paragraph. Each component is read through the source map (MB.mread), not from a
    directory listing. Also returns {stem: level-1 heading} for the docx cross-check."""
    out, h1s = {}, {}
    for comp, stem in S1_COMP:
        sub = ""
        todo = [b for b in re.split(r"\n\s*\n", PG.TAG.sub("", MB.mread(comp)).strip())]
        while todo:
            b = todo.pop(0).strip()
            if not b:
                continue
            if b.startswith("## "):
                h1s.setdefault(stem, b[3:].split("\n", 1)[0].strip())
                sub = ""
                rest = b.split("\n", 1)[1] if "\n" in b else ""
                if rest.strip():
                    todo.insert(0, rest)
            elif b.startswith("### "):
                head, _, rest = b.partition("\n")
                sub = head[4:].strip().replace("*", "")      # Stage 7 A2.4d: gene italics in a heading are markup; the docx heading carries none
                if rest.strip():
                    out.setdefault((stem, sub), []).append(" ".join(rest.split("\n")).strip())
            else:
                out.setdefault((stem, sub), []).append(" ".join(b.split("\n")).strip())
    return out, h1s


PARA, H1 = section_paragraphs()


def docx_paragraphs(path):
    """{(level-1 heading, level-2 heading or ''): [text of each non-empty body paragraph]} in a built manuscript .docx."""
    import docx as _docx
    n, cur, h1 = {}, None, None
    for p in _docx.Document(path).paragraphs:
        st = p.style.name
        if st == "Heading 1":
            h1 = p.text.strip(); cur = (h1, "")
        elif st == "Heading 2":
            cur = (h1, p.text.strip())
        elif cur and p.text.strip():
            n.setdefault(cur, []).append(p.text.strip())
    assert n, ("no paragraphs read from", path)
    return n


# cross-representation check: the paragraphs that Table S1 locates must be those of the manuscript AS BUILT. Compare the mirrored split
# with a built .docx (P31R_AF1_DOCX, else <OUT>/P31R_Manuscript_R2_preview.docx) when one exists that is newer than every map component
# the sheet reads (including the title page, for the Keywords paragraph); packaging mode (P31R_AF1_FINAL=1) refuses without one.
_docx = os.environ.get("P31R_AF1_DOCX") or os.path.join(OUT, "P31R_Manuscript_R2_preview.docx")
_newest = max(os.path.getmtime(f) for c_ in [c for c, _ in S1_COMP] + ["title_page"] for f in MB.MAP[c_]["files"])
DOCXP = None
if os.path.exists(_docx) and os.path.getmtime(_docx) >= _newest:
    DOCXP = docx_paragraphs(_docx)
    # the builder adds the title page's Keywords line as one more paragraph under the Abstract heading
    _kw = int(bool(re.search(r"\*\*Keywords:\*\*", MB.mread("title_page"))))
    _cmp = [((s, sub), len(ps) + (_kw if s == "Abstract" else 0), len(DOCXP.get((H1[s], sub), []))) for (s, sub), ps in PARA.items()]
    _bad = [c for c in _cmp if c[1] != c[2]]
    assert not _bad, ("S1 paragraph split differs from the built manuscript (markdown, docx):", _bad)
    DOCX_CHECK = "S1 paragraph numbering and anchors checked against %s: %d (sub)sections identical" % (_docx, len(_cmp))
else:
    DOCX_CHECK = "!!! S1 paragraph numbering NOT checked against a built manuscript (no .docx newer than the map components at %s)" % _docx
    if os.environ.get("P31R_AF1_FINAL") == "1":
        sys.exit("REFUSING TO WRITE (final mode): " + DOCX_CHECK)
RES1 = "81 colocalization-supported pairs, graded by internal support and exposure-side corroboration"
RESX = "External data corroborate selected T2D associations"
RESA = "ARG1 substrate, product and diabetes associations across two signals"
RESE = "Direction transferred to East Asian populations more widely than significance"
RESK = "Kidney pairs among the exploratory multi-signal candidates, with UMOD as a control"

S1LOG = []   # (item location kind, stem, subsection, paragraph, anchor, occurrences in the markdown, occurrences in the docx or None)


def loc(stem, sub, k, anchor):
    """A manuscript location, asserted: the paragraph exists and holds the anchor phrase EXACTLY ONCE, in the markdown split and, when a
    current built manuscript exists, in the same paragraph of the built .docx. Anchors carry no markdown, citation, ledger-tag or {SM:} text,
    so the two renderings of a paragraph contain them identically."""
    ps = PARA.get((stem, sub))
    assert ps and len(ps) >= k, ("S1 location missing", stem, sub, k)
    n_md = ps[k - 1].count(anchor)
    assert n_md != 0, ("S1 anchor not in markdown paragraph", stem, sub, k, anchor)
    assert n_md == 1, ("S1 anchor repeated in markdown paragraph", stem, sub, k, anchor, n_md)
    n_dx = None
    if DOCXP is not None:
        dps = DOCXP.get((H1[stem], sub), [])
        assert len(dps) >= k, ("S1 location missing in the built docx", stem, sub, k)
        n_dx = dps[k - 1].count(anchor)
        assert n_dx != 0, ("S1 anchor not in built-docx paragraph", stem, sub, k, anchor)
        assert n_dx == 1, ("S1 anchor repeated in built-docx paragraph", stem, sub, k, anchor, n_dx)
    S1LOG.append(("manuscript", stem, sub, k, anchor, n_md, n_dx))
    if stem == "Abstract":
        return "Abstract, %s" % re.match(r"\*\*(\w+):", ps[k - 1]).group(1)
    if stem == "Declarations":
        return "Declarations, %s" % sub
    return ("%s%s" % ({"ResearchInsights": "Research Insights"}.get(stem, stem), (": " + sub) if sub else ""), k)


def sm_sections(path):
    """{'S-R10': section text, ...}: the '### Section S-Rn.' sections of the R2 Results Detail (the source of Additional file 2), tags stripped."""
    parts = re.split(r"(?m)^### Section (S-R\d+)\.[^\n]*\n", PG.TAG.sub("", io.open(path, encoding="utf-8").read()))
    return {parts[i]: parts[i + 1] for i in range(1, len(parts), 2)}


SMR = sm_sections(SRC["ms_af2_results_detail"])


def smloc(sec, anchor):
    """An Additional file 2 location for content that the main text moved there: the section exists in the R2 Results Detail source and
    holds the anchor exactly once. Checked against the markdown source only (no Additional file 2 is built from the R2 text yet)."""
    assert sec in SMR, ("S1 Additional file 2 section missing", sec)
    n = SMR[sec].count(anchor)
    assert n != 0, ("S1 anchor not in Additional file 2 section", sec, anchor)
    assert n == 1, ("S1 anchor repeated in Additional file 2 section", sec, anchor, n)
    S1LOG.append(("additional file 2", "Additional file 2", sec, 0, anchor, n, None))
    return "Additional file 2, Section %s" % sec


def where(items):
    """Join locations; paragraphs of one (sub)section are merged: 'Methods: Colocalization, paragraphs 1 and 3'."""
    out, grp = [], {}
    for it in items:
        if isinstance(it, tuple):
            if it[0] not in grp:
                grp[it[0]] = []
                out.append(it[0])
            if it[1] not in grp[it[0]]:
                grp[it[0]].append(it[1])
        elif it not in out:
            out.append(it)
    txt = []
    for o in out:
        if o in grp:
            ks = sorted(grp[o])
            txt.append("%s, paragraph%s %s" % (o, "s" if len(ks) > 1 else "", ", ".join(map(str, ks[:-1])) + (" and " if len(ks) > 1 else "") + str(ks[-1])))
        else:
            txt.append(o)
    return "; ".join(txt)


def has(stem, sub, anchor):
    return any(anchor in p for p in PARA.get((stem, sub), []))


IVA = "assumes that the instrument"   # anchor of the requested Methods sentence that states the three IV assumptions
iv_stated = any(IVA in p for (s, _), ps in PARA.items() if s == "Methods" for p in ps)
M = "Methods"
S1MAP = {
    ("1", ""): ([loc("Abstract", "", 2, "two-sample Mendelian randomization")], "reported",
                "The title names the graded resource; the Abstract names the design."),
    ("2", ""): ([loc("Background", "", 1, "continuum"), loc("Background", "", 2, "cis protein quantitative trait loci"),
                 loc("Background", "", 3, "Interpreting these maps requires distinguishing support for the cis instrument")], "reported", "--"),
    ("3", ""): ([loc("Background", "", 4, "We therefore built a genetic evidence map of circulating proteins"),
                 loc(M, "External outcome checks", 1, "Before any external outcome result was read"),
                 loc(M, "Kidney and metabolite loci", 3, "a rule written before the signal-specific T2D estimates were computed")],
                "reported" if iv_stated else "partly reported",
                "Objectives are stated. The external outcome checks and the arginase-1 decision rule were fixed before their results were read; "
                "the six externally checked pairs were chosen after discovery." + ("" if iv_stated else
                " What MR estimates under its assumptions is not stated as such.")),
    ("4", ""): ([loc(M, "Study design and data sources", 1, "The association analyses used published or publicly released summary statistics"), "Table 1; Additional file 1: Table S2"], "reported", "--"),
    ("4", "a"): ([loc(M, "Study design and data sources", 1, "participants of European ancestry"),
                  loc(M, "Study design and data sources", 2, "None of these resources includes UK Biobank"), "Table 1; Additional file 1: Table S2"],
                 "partly reported", "Design and populations are given; locations and recruitment dates are those of the source studies and are not restated."),
    ("4", "b"): (["Table 1; Additional file 1: Table S2", loc(M, "External outcome checks", 2, "were checked in DIAMANTE 2018 without UK Biobank")], "partly reported",
                 "Sample sizes are given; a power calculation before the main analysis is not reported."),
    ("4", "c"): ([loc(M, "Instruments and harmonization", 1, "imputation INFO")], "reported", "--"),
    ("4", "d"): ([loc(M, "Study design and data sources", 2, "The eight outcomes were"),
                  loc(M, "Study design and data sources", 3, "ln(eGFR) per standard deviation"), "Table 1"], "reported", "--"),
    ("4", "e"): ([loc("Declarations", "Ethics approval and consent to participate", 1, "previously published or publicly released, de-identified summary-level genetic-association statistics")], "reported", "--"),
    ("5", ""): (([next(("%s: %s" % (M, sub), i + 1) for (s, sub), ps in PARA.items() if s == M for i, p in enumerate(ps) if IVA in p)]
                 if iv_stated else
                 [loc(M, "Colocalization", 1, "shared by protein and outcome"),
                  loc(M, "Colocalization", 2, "Single-variant colocalization assumes one causal variant per trait in the region"),
                  loc(M, "Directionality, intermediate traits and annotation", 1, "The Steiger test")]),
                "reported" if iv_stated else "partly reported",
                "--" if iv_stated else "The three core assumptions are not stated as such; colocalization, the Steiger test and the intermediate-trait scan address parts of them."),
    ("6", ""): ([loc(M, "Screen and multiple testing", 1, "two-sample Mendelian randomization (MR) with TwoSampleMR")], "reported", "--"),
    ("6", "a"): ([loc(M, "Statistical reporting", 1, "Effects are log-odds of the outcome per SD of protein"),
                  loc(M, "Study design and data sources", 3, "on a risk-aligned scale")], "reported", "--"),
    ("6", "b"): ([loc(M, "Instruments and harmonization", 1, "Alleles were aligned to the protein-raising allele")], "reported", "--"),
    ("6", "c"): ([loc(M, "Screen and multiple testing", 1, "We used the Wald ratio")], "partly reported",
                 "The estimator is stated; covariate adjustment in the exposure and outcome GWAS is not described."),
    ("6", "d"): ([loc(M, "Instruments and harmonization", 1, "a proxy in high LD"), loc(M, "Screen and multiple testing", 1, "Untestable combinations counted as non-significant")],
                 "reported", "--"),
    ("6", "e"): ([loc(M, "Screen and multiple testing", 1, "Benjamini–Hochberg")], "reported", "--"),
    ("7", ""): ([loc(M, "Colocalization", 1, "coloc.abf"), loc(M, "Colocalization", 2, "colocalization under the sum of single effects (SuSiE) model"),
                 loc(M, "Directionality, intermediate traits and annotation", 1, "The Steiger test"),
                 loc(M, "Exposure-side corroboration", 1, "deCODE"), loc("Discussion", "", 11, "one cis instrument per pair")], "partly reported",
                "Colocalization, the Steiger test and exposure-side corroboration address relevance and confounding by linkage. With one cis instrument per pair, "
                "pleiotropy-robust estimators were unavailable, so the exclusion restriction is not assessed directly."),
    ("8", ""): ([loc(M, "Instruments and harmonization", 1, "the analyses without proxies and a genome-wide-threshold sensitivity"),
                 loc(M, "Colocalization", 1, "conservative prior"),
                 loc(M, "Colocalization", 3, "decoy-locus stress test"), loc(M, "External outcome checks", 1, "an external outcome check"),
                 loc(M, "External regional colocalization", 1, "We repeated colocalization"),
                 loc(M, "Signal-level analyses", 1, "Signal-level analyses asked whether"),
                 loc(M, "East Asian analyses", 1, "strict sensitivity analysis"), loc(M, "Kidney and metabolite loci", 1, "cystatin-C")], "reported", "--"),
    ("9", ""): (["see 9a and 9b"], "see subitems", "--"),
    ("9", "a"): ([loc(M, "Statistical reporting", 1, "Analyses used R (TwoSampleMR, coloc, susieR)"), "Additional file 2"], "partly reported",
                 "Software is named; versions are given in Additional file 2 where they were recorded, and several runs did not record their package versions."),
    ("9", "b"): ([loc(M, "External outcome checks", 1, "Before any external outcome result was read"),
                  loc(M, "Kidney and metabolite loci", 3, "a rule written before the signal-specific T2D estimates were computed"),
                  loc(M, "Pre-declared follow-up analyses", 1, "in writing before their final runs")], "partly reported",
                 "No public protocol pre-registration is reported. The follow-up analyses were declared in writing before their final runs, "
                 "but not before all of their values were known; the external outcome checks and the arginase-1 decision rule were written down before their results."),
    ("10", ""): (["see 10a to 10d"], "see subitems", "--"),
    ("10", "a"): ([loc("Results", RES1, 1, "protein–outcome tests"), "Figure 1"], "reported", "--"),
    ("10", "b"): (["Table 1; Additional file 1: Table S2"], "partly reported",
                  "Summary statistics only: sample sizes are given; phenotype distributions are those of the source studies."),
    ("10", "c"): (["--"], "not reported", "The outcome GWAS were used as released."),
    ("10", "d"): ([loc(M, "Study design and data sources", 2, "exposure and outcome samples do not overlap"),
                   loc(M, "External outcome checks", 3, "cross-trait LD-score-regression intercept"),
                   loc("Results", RESX, 2, "sampling-noise correlation")], "reported",
                  "Exposure and discovery outcome samples are both of European ancestry."),
    ("11", ""): (["see 11a to 11d"], "see subitems", "--"),
    ("11", "a"): (["Additional file 1: Table S3"], "partly reported",
                  "Variant-exposure associations (effect and SE) are in Table S3; variant-outcome associations are not printed and equal the instrument effect times the MR estimate. "
                  "For the %s proxy pairs, the proxy, its alleles and its r2 with the lead are in Table S3." % spell(NPX)),
    ("11", "b"): ([loc("Results", RESX, 1, "effect estimates and p values in Figure 2b"), loc("Results", RESA, 1, "log-odds per SD"),
                   "Figure 2; Additional file 1: Tables S3 and S5"], "reported", "--"),
    ("11", "c"): (["--"], "not applicable", "Estimates are per SD of protein and are not translated into absolute risk."),
    ("11", "d"): (["Figures 2 to 5"], "reported", "--"),
    ("12", ""): (["see 12a and 12b"], "see subitems", "--"),
    ("12", "a"): ([loc("Results", RES1, 4, "pairs with a supportive multi-signal class"),
                   loc("Results", RES1, 2, "Exposure-side corroboration in deCODE or ARIC was evaluable")], "reported", "--"),
    ("12", "b"): ([loc(M, "Instruments and harmonization", 1, "Every colocalization-supported pair rests on one lead cis variant"),
                   loc("Results", RESX, 4, "gave negative conditional protein-scaled CAD estimates"), loc("Results", RESK, 3, "behaved as a control"),
                   loc("Results", RESA, 2, "Nearly independent signals"), "Table 3"], "partly reported",
                  "Heterogeneity across independent signals is tested at the multi-signal loci (TGFB1, UMOD, ARG1 and DPEP1; Table 3). The 81 pairs rest on one "
                  "instrument each, so heterogeneity across their variants cannot be assessed."),
    ("13", ""): (["see 13a to 13e"], "see subitems", "--"),
    ("13", "a"): ([loc("Results", RES1, 2, "Under the conservative prior"), smloc("S-R10", "pairs kept an instrument at p")], "reported", "--"),
    ("13", "b"): ([loc("Results", RESX, 2, "A post hoc extension of the same check"),
                   smloc("S-R13", "Keeping only outcome files that report per-variant sample sizes"),
                   loc("Results", RESE, 2, "including a resource swap that excludes CKB"),
                   smloc("S-R13", "Replacing outcome resources that include CKB with Biobank Japan"),
                   loc("Results", RESX, 3, "In external regional colocalization"), loc("Results", RESX, 4, "gave negative conditional protein-scaled CAD estimates"),
                   loc("Results", RESK, 3, "behaved as a control"), loc("Results", RESA, 2, "Nearly independent signals"), "Table 3"], "reported",
                  "The strict East Asian subset and the CKB-free swap are reported in Additional file 2, Section S-R13, and per pair in Additional file 1: Tables S7a and S7b. "
                  "The external regional colocalizations and the signal-level analyses (TGFB1, UMOD, ARG1, DPEP1) are in the Results and Table 3."),
    ("13", "c"): ([loc(M, "Directionality, intermediate traits and annotation", 1, "The Steiger test"),
                   smloc("S-R10", "Direction and instrument strength"), "Additional file 1: Table S3"], "reported",
                  "The Steiger result is reported in Additional file 2 and, per pair, in Table S3; the main-text Results do not repeat it."),
    ("13", "d"): ([loc("Discussion", "", 4, "This direction is as expected if loss of PAM function raises risk"),
                   loc("Discussion", "", 7, "Mouse models suggest that TGF-β signaling protects against atherosclerosis")], "partly reported",
                  "Direction is compared with experimental and earlier genetic reports for selected pairs (PAM, TGFB1); no systematic comparison with non-MR effect estimates is made."),
    ("13", "e"): (["--"], "not applicable", "Single instruments: leave-one-out analyses are not possible."),
    ("14", ""): ([loc("Discussion", "", 1, "This study organizes the genetic evidence for"), loc("Discussion", "", 8, "Outside data corroborated"),
                  loc("Discussion", "", 10, "The kidney and ancestry analyses show two limits")], "reported",
                 "The Discussion restates the design and the corroborated pairs and where the kidney and ancestry analyses reach; the grades themselves are summarised "
                 "in the Results (map subsection) and Figure 1."),
    ("15", ""): ([loc("Discussion", "", 11, "This study has several limitations")], "reported", "--"),
    ("16", ""): (["see 16a to 16c"], "see subitems", "--"),
    ("16", "a"): ([loc("Discussion", "", 1, "It identifies which cis instrument supports each association"), loc("Discussion", "", 2, "cover several kinds of biology"),
                   loc("Discussion", "", 12, "would not by themselves establish mediation")], "reported", "--"),
    ("16", "b"): ([loc("Discussion", "", 4, "PAM amidates the C-terminus of many peptide hormones"), loc("Discussion", "", 2, "the dipeptidase DPEP1"),
                   loc("Discussion", "", 5, "ARG1 hydrolyzes arginine to urea and ornithine"),
                   loc("Discussion", "", 3, "a fibroblast growth factor receptor involved in bile acid and lipid metabolism")], "partly reported",
                  "Mechanisms are discussed; gene-environment equivalence is not discussed as such."),
    ("16", "c"): ([loc("Discussion", "", 1, "point to the molecular phenotype and disease outcome to test next"), loc("Conclusions", "", 1, "The map helps researchers prioritize protein–outcome pairs")],
                  "partly reported", "Clinical relevance is addressed through research follow-up: the map prioritizes protein–outcome pairs for follow-up and does not by itself identify a protein-mediated mechanism. "
                  "Effect sizes of possible interventions are not discussed."),
    ("17", ""): ([loc("Discussion", "", 10, "direction transferred more widely than significance"),
                  loc("Discussion", "", 11, "which may not match Finnish or mixed-ancestry outcome data")], "partly reported",
                 "Generalisability across ancestries is discussed; across exposure periods and levels it is not."),
    ("18", ""): ([loc("Declarations", "Funding", 1, "This work was supported")], "reported", "Funding and the funders' role are stated in the Declarations."),
    ("19", ""): ([loc("Declarations", "Availability of data and materials", 1, "published or publicly released genetic-association summary statistics"),
                  loc("Declarations", "Availability of data and materials", 2, "The repository contains the analysis and figure scripts")], "reported",
                 "Data sources, the repository contents, the repository URL and the archive DOI (the Zenodo concept DOI) are stated in the Declarations."),
    ("20", ""): ([loc("Declarations", "Competing interests", 1, "The authors declare that")], "reported", "Competing interests of all authors are stated in the Declarations."),
}
_it = pd.read_csv(SRC["strobe"], sep="\t", dtype=str).fillna("")
eq([int((_it["sub"] == "").sum()), int((_it["sub"] != "").sum())], LV("L453"), "STROBE-MR items / subitems (L453)")
assert set(zip(_it["item"], _it["sub"])) == set(S1MAP), "S1 map does not cover exactly the checklist items"
s1 = pd.DataFrame([dict(Item=r["item"] + r["sub"], Part=r["part"], Heading=r["heading"], **{"Checklist item (short form)": r["label"]},
                        **{"Where addressed": where(S1MAP[(r["item"], r["sub"])][0]), "Status": S1MAP[(r["item"], r["sub"])][1],
                           "Note": S1MAP[(r["item"], r["sub"])][2]}) for _, r in _it.iterrows()])
# G5b (2026-11-06): the Status legend lists the values the S1 table uses, in a fixed order, instead of a typed list ('author-supplied' stayed
# listed after no row carried it); a value outside the order fails the build.
_S1_STATUS_ORDER = ("reported", "partly reported", "not reported", "not applicable", "see subitems")
_s1_used = set(s1["Status"])
assert _s1_used <= set(_S1_STATUS_ORDER), ("S1 Status value outside the legend order", sorted(_s1_used - set(_S1_STATUS_ORDER)))
S1_STATUS_LEGEND = " / ".join(x for x in _S1_STATUS_ORDER if x in _s1_used)
ni, ns = LV("L453")
sheet("S1", "Table S1. STROBE-MR checklist.",
      "The %d items and %d subitems of the STROBE-MR statement (Skrivankova et al., JAMA 2021; PMID 34698778), in short form; the full wording is in the statement and its explanation paper (BMJ 2021; PMID 34702754). For each item: where the manuscript addresses it, by section and paragraph (paragraphs are counted within each subsection), and whether it is reported. Content that the main text places in Additional file 2 is located by section there." % (ni, ns),
      s1, {"Item": "checklist item and subitem letter", "Part": "part of the checklist", "Heading": "item heading in the statement",
           "Checklist item (short form)": "what the item asks for, shortened", "Where addressed": "manuscript section and paragraph (paragraphs counted within each subsection), table or figure",
           "Status": S1_STATUS_LEGEND, "Note": "what is or is not covered"}, S1_SRCS)

# ================================================================================================ S2 data sources
s2_note = ""
if os.path.exists(SRC["tables2"]):
    import P31R_tablegate as TG
    _n, _title, _hdr, _body, _notes = TG.parse_table(io.open(SRC["tables2"], encoding="utf-8").read())
    assert _n == "S2", _n
    CITE = re.compile(r"((?:\s?\[@[^\]]+\])+)")
    keys = list(dict.fromkeys(k for t in [_title, _notes] + _hdr + [c for r in _body for c in r] for k in re.findall(r"\[@([^\]]+)\]", t)))
    # citations become self-contained "(Author et al. Year, PMID n)": Table S2 cites sources that the main text does not, so main-text
    # reference numbers cannot serve; entries come from the manuscript builder's own resolver (checked list by title, else PubMed)
    numbering = {}
    if keys:
        import csv as _csv
        import P31R_build as MB
        ents, _prov = MB.resolve(keys, MB.current_refs()[0], MB.old_titles())
        oldpm = {"old:" + r["n"]: r["pmid"] for r in _csv.DictReader(io.open(MB.OLD_TSV, encoding="utf-8"), delimiter="	")}

        def doi_pmid(doi, preprint=False):
            """PMID for a DOI, from PubMed itself (esearch on [doi]), confirmed by fetching the record and comparing its DOI. A preprint
            that PubMed does not hold (no hit on its DOI) has none: '' (the caller then cites its DOI)."""
            import urllib.parse
            q = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "term": doi + "[doi]"})
            # NCBI answers a rate limit with HTTP 200 and an error body (no idlist), which refcheck._get does not retry: retry a bounded
            # number of times here, then fail loudly (Stage 7 A3-fix: the mutation ladder's burst of AF1 builds hit it twice)
            import time as _time
            for _try in range(4):
                _r = json.loads(MB.refcheck._get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + q))
                if "idlist" in _r.get("esearchresult", {}):
                    break
                _time.sleep(2.0 * (_try + 1))
            else:
                raise AssertionError(("PubMed esearch returned no idlist after 4 tries", doi, _r))
            ids = _r["esearchresult"]["idlist"]
            if preprint and not ids:
                return ""
            assert len(ids) == 1, ("PubMed does not resolve the DOI to exactly one record", doi, ids)
            rec = MB.refcheck.fetch_pubmed(ids)[ids[0]]
            assert (rec.get("doi") or "").lower() == doi.lower(), ("PubMed record carries another DOI", doi, ids[0], rec.get("doi"))
            return ids[0]
        S2_DOI_PMID, S2_PREPRINT = {}, set()
        for k, e in zip(keys, ents):
            first = e.split(",")[0].strip()
            mm = re.match(r"^(.+?)\s[A-Z]{1,3}$", first)
            au = mm.group(1) if mm and not re.search(r"Consortium|Group|Collaboration|Investigators", first) else first
            yy = re.search(r"\.\s(?:medRxiv\s|bioRxiv\s)?((?:19|20)\d\d)[;.\s]", e)
            assert yy, ("cannot read the year", k, e[:100])
            kk = MB.ALIAS.get(k, k)
            pm = MB.NEW_PMID.get(k) or oldpm.get(kk, "")
            doi = re.search(r"doi:\s?(\S+?)\.?$", e)
            if not pm and doi:                       # the reference set holds no PMID: derive it from the DOI in PubMed
                preprint = bool(re.search(r"\b(?:medRxiv|bioRxiv)\b", e))
                pm = S2_DOI_PMID[k] = doi_pmid(doi.group(1), preprint)
                if preprint and not pm:              # a preprint PubMed does not hold is cited by its DOI
                    S2_PREPRINT.add(k)
            ident = ("PMID %s" % pm) if pm else ("doi:%s" % doi.group(1) if doi else "")
            assert ident, ("no PMID or DOI for", k)
            group = au == first                                     # a group author keeps its full name and no "et al."
            many = "et al." in e or len(e.split(". ")[0].split(", ")) > 1
            numbering[k] = "%s%s %s, %s" % (au, " et al." if many and not group else "", yy.group(1), ident)
        assert all(re.search(r", PMID [0-9]+$", v) for k, v in numbering.items() if k not in S2_PREPRINT), (
            [k for k, v in numbering.items() if k not in S2_PREPRINT and not re.search(r", PMID [0-9]+$", v)])
        assert all(re.search(r", doi:10\.\S+$", numbering[k]) for k in S2_PREPRINT), [numbering[k] for k in S2_PREPRINT]
        s2_note = " Citations give the first author, year and PMID of each source" + ("; a preprint that PubMed does not hold is cited by its DOI." if S2_PREPRINT else ".")

    def clean(t):
        t = PG.TAG.sub("", t).replace("**", "")
        if numbering:
            t = CITE.sub(lambda m: " (" + "; ".join(numbering[k] for k in re.findall(r"\[@([^\]]+)\]", m.group(1))) + ")", t)
        return t.strip()
    s2 = pd.DataFrame([[clean(c) for c in r] for r in _body], columns=[clean(h) for h in _hdr])
    s2_title = re.sub(r"^Table S2\.\s*", "", clean(_title)).strip().rstrip(".") + "."
    s2_legend = (" ".join(clean(x) for x in re.split(r"\n\s*\n", _notes) if x.strip()) or "--") + s2_note
    S2COLS = {"Analysis": "analysis the resource serves, in the order of the Methods",
              "Resource": "study or consortium, with its citation (first author, year, PMID) and the outcome or trait where one resource serves several",
              "File or accession": "summary-statistics file, GWAS Catalog accession or release identifier used",
              "Ancestry": "ancestry of the participants in the summary statistics",
              "Trait and scale": "trait definition and the scale of its effects (log-odds for binary traits; SD units or ln-transformed units for quantitative traits)",
              "Cases / controls or N": "cases / controls for binary traits, N for quantitative traits (counts)",
              "Overlap with the exposure sample": "participants shared with the protein (exposure) sample: UKB-PPP for the European instrument, CKB for the ancestry-matched analysis",
              "Sample size read from": "where the sample size was read"}
    assert list(s2.columns) == list(S2COLS), ("Table S2 columns changed: describe them in S2COLS", list(s2.columns))
    sheet("S2", "Table S2. " + s2_title, s2_legend, s2, S2COLS, ["tables2"])
    S2_PLACEHOLDER = False
    # East Asian resources are named in S7a/S7b exactly as Table S2 names them (citations stripped), matched by file or accession
    S2RES = [(PG.TAG.sub("", r[0]).strip(), re.sub(r"\s*\[@[^\]]+\]", "", PG.TAG.sub("", r[1])).replace("**", "").strip(), PG.TAG.sub("", r[2]).strip())
             for r in _body]
else:
    s2 = pd.DataFrame({"Status": ["PLACEHOLDER: tables/P31R_TableS2.md (all data sources) is written by another Stage 3b track; rebuild this workbook once it exists."]})
    sheet("S2", "Table S2. Data sources [PLACEHOLDER].", "Not yet built.", s2, {"Status": "placeholder"}, [])
    S2_PLACEHOLDER = True
    S2RES = []


def _nrm(t):
    return re.sub(r"[^a-z0-9]", "", t.lower())


def s2_name(label, prefer):
    """The Table S2 'Resource' name of an outcome file, matched on its file or accession within the first analysis (by prefix)
    that matches; exactly one row must match."""
    n = _nrm(label)
    for pre in prefer:
        hits = sorted({res for an, res, fi in S2RES if an.startswith(pre) and len(_nrm(fi)) >= 6 and (n in _nrm(fi) or _nrm(fi) in n)})
        if hits:
            assert len(hits) == 1, ("ambiguous Table S2 match", label, pre, hits)
            return hits[0]
    raise AssertionError("no Table S2 row for the outcome file %r (S7 names its resources from Table S2)" % label)

# ================================================================================================ S3 the 81 pairs
pe = pd.read_csv(SRC["tiered"])
mr = pd.read_csv(SRC["mr"])
prim = mr[mr.method.isin(["Wald ratio", "Inverse variance weighted"])]
ins = pd.read_csv(SRC["ins"])
vep = pd.read_csv(SRC["vep"])
pxr = pd.read_csv(SRC["proxy"], sep="\t", dtype=str)
eq(len(pe), N81, "pairs (L010)")
eq(pe.gene.nunique(), N62, "proteins (L011)")
eq(pe.outcome.value_counts().to_dict(), LV("L013"), "pairs per outcome (L013)")
eq(pe.evidence_tier.value_counts().to_dict(), LV("L014"), "tiers (L014)")
x = pe.merge(prim[["gene", "outcome_name", "b", "se", "pval", "any_proxy", "n_snps", "method"]], left_on=["gene", "outcome"],
             right_on=["gene", "outcome_name"], how="left", validate="1:1")
assert x.b.notna().all(), "a pair has no primary MR row"
flip = np.where(x.outcome == "eGFR", -1.0, 1.0)
# the tier table stores mr_beta rounded to 6 decimals: print the CANONICAL full-precision MR-table value, and assert the rounded copy agrees
assert ((x.mr_beta - flip * x.b).abs() <= 5e-7 + 1e-12).all(), "tier-table estimate is not the MR estimate (eGFR sign-flipped, L078)"
assert np.allclose(x.mr_pval, x.pval, rtol=1e-5, atol=0), "tier-table p differs from the MR table"
x["b_star"] = flip * x.b
eq([int((x.outcome == "eGFR").sum())] * 2, LV("L078")[:2], "eGFR pairs sign-flipped (L078)")
assert (x.method == "Wald ratio").all() and (x.n_snps == 1).all(), "a colocalization-supported pair has more than one instrument"
# conservative-prior band from the VALUE (the tier table's prior_tier is a different construct and is not used)
h6 = x["PP.H4_p12_1e6"]
band = np.where(h6 >= RULE_H4, "robust (PP.H4 ≥ %s)" % T80, np.where(h6 >= HALF, "PP.H4 %.2f to %s" % (HALF, T80), "PP.H4 below %.2f" % HALF))
eq({"robust": int((h6 >= RULE_H4).sum()), "atten": int(((h6 >= HALF) & (h6 < RULE_H4)).sum()), "fragile": int((h6 < HALF).sum())}, LV("L016"), "conservative-prior bands (L016)")
assert (x.is_severe_prior_instab.astype(bool) == (h6 < HALF)).all(), "weakness threshold differs from the tier table flag"
eq(x.susie_class.value_counts().to_dict(), LV("L017"), "SuSiE classes (L017)")
eq(x.steiger_result.value_counts().to_dict(), LV("L019"), "Steiger (L019)")
msc = x.mr_sig_class.value_counts().to_dict()
eq([msc["Bonferroni"], msc["FDR"], msc["Suggestive"]], LV("L959")[:3], "screen significance (L959)")
eq(x.loc[x.mr_sig_class == "Suggestive", ["gene", "outcome"]].apply("-".join, axis=1).tolist(), [LV("L959")[3]], "the nominal-only pair (L959)")
# weaknesses (Methods: the five weaknesses that make a pair Flagged)
WK = [("is_severe_prior_instab", "PP.H4 below one half under the conservative prior"), ("is_susie_discordant", "discordant SuSiE-based colocalization"),
      ("is_replication_counter", "exposure-side arm with the opposite sign"), ("is_assay_discordant", "discordant assays within one platform"),
      ("is_steiger_borderline", "Steiger direction not forward")]
wk = x.apply(lambda r: "; ".join(t for c, t in WK if bool(r[c])) or "--", axis=1)
# the one reason the evidence table names per Flagged pair: the first weakness in the order SuSiE, opposing arm, assays, Steiger,
# conservative prior (as Additional file 2 {SM:grading} states), asserted against the stored reason
WKORD = ["is_susie_discordant", "is_replication_counter", "is_assay_discordant", "is_steiger_borderline", "is_severe_prior_instab"]
RSN = {"SuSiE discordant": "is_susie_discordant", "Replication counter-evidence": "is_replication_counter", "Assay discordant": "is_assay_discordant",
       "Steiger borderline": "is_steiger_borderline", "Prior-fragile": "is_severe_prior_instab"}
first = x.apply(lambda r: next((c for c in WKORD if bool(r[c])), None), axis=1)
assert (first.fillna("--") == x.exploratory_reason.map(lambda v: RSN[v] if isinstance(v, str) else "--")).all(), "named reason is not the first weakness in the stated order"
named = first.map(lambda c: dict(WK)[c] if isinstance(c, str) else "--")
eq(bool(((wk != "--") == (x.evidence_tier == "Flagged")).all()), True, "a weakness exactly when Flagged")
eq(int(((wk != "--") & (x.evidence_tier == "Flagged")).sum()), LV("L014")["Flagged"], "Flagged (L014)")
eq(int((x.is_severe_prior_instab.astype(bool) & (x[[c for c, _ in WK]].astype(bool).sum(axis=1) == 1)).sum()), LV("L021")["Prior-fragile"],
   "below one half as the only weakness (L021)")
crit = x[["crit_bonferroni", "crit_prior_stable", "crit_susie_supportive", "crit_abf_dual"]].astype(bool)
assert (crit.all(axis=1) == x.core_complete.astype(bool)).all() and x.crit_abf_dual.astype(bool).all()
assert (x.crit_prior_stable.astype(bool) == (h6 >= RULE_H4)).all() and (x.crit_bonferroni.astype(bool) == (x.mr_sig_class == "Bonferroni")).all()
supP = x[x.evidence_tier == "Supported"]
eq([int((supP.replication_depth == "R0").sum()), len(supP), int(((supP.replication_depth == "R0") & supP.core_complete.astype(bool)).sum())],
   LV("L226"), "Supported without an evaluable arm (L226)")
# exposure-side corroboration (labels as in the Figure 2 legend and Results)
DEC = {"Concordant (deCODE)": "same direction", "Discordant (deCODE)": "opposite direction", "Not evaluable (weak instrument)": "weak instrument (F < %s)" % F10,
       "Not in deCODE release": "not assayed", "No usable variant (deCODE)": "no usable variant"}
ARI = {"Concordant (ARIC)": "same direction", "Discordant (ARIC)": "opposite direction", "Not evaluable (weak instrument)": "weak instrument (F < %s)" % F10,
       "Not on ARIC panel": "not assayed", "No usable variant (ARIC)": "no usable variant"}
XP = {"concordant": "same direction", "discordant": "opposite direction", "conflicting": "cohorts disagree", "not assayed (panel gap)": "not assayed",
      "not evaluable (weak instrument)": "weak instrument", "no usable variant": "no usable variant", "assay discordant (within platform)": "assays disagree"}
eq(int((x.decode_status != "NOT_IN_RELEASE").sum() - (x.decode_status == "NO_EXACT_RSID_MATCH").sum()), LV("L030"), "deCODE-testable (L030)")
eq(int((x.decode_replication_class == "Concordant (deCODE)").sum()), LV("L031"), "deCODE criteria met (L031)")
eq([int(x.aric_evaluable_pair.map(tf).sum()), int((x.aric_replication_class == "Concordant (ARIC)").sum()),
    int((x.aric_replication_class == "Discordant (ARIC)").sum())], LV("L033"), "ARIC evaluable / same / opposite (L033)")
xpc = x.crossplatform_class.value_counts().to_dict()
eq([xpc["concordant"], xpc["discordant"], xpc["conflicting"]], LV("L034"), "combined exposure-side classes (L034)")
eq([xpc["concordant"] + xpc["discordant"] + xpc["conflicting"], N81 - (xpc["concordant"] + xpc["discordant"] + xpc["conflicting"]),
    xpc["not assayed (panel gap)"], xpc["not evaluable (weak instrument)"], xpc["no usable variant"], xpc["assay discordant (within platform)"]],
   LV("L259"), "exposure-side evaluable / not (L259)")
# instruments (UKB-PPP), aligned to the protein-raising allele
x["SNP_ok"] = [((ins.gene == g) & (ins.SNP == s)).sum() == 1 for g, s in zip(x.gene, x.lead_snps)]
assert x.SNP_ok.all(), "a lead variant is not the protein's UKB-PPP instrument"
iv = x[["gene", "lead_snps"]].merge(ins, left_on=["gene", "lead_snps"], right_on=["gene", "SNP"], how="left", validate="m:1")
up = iv.beta > 0
ra = np.where(up, iv.effect_allele, iv.other_allele); oa = np.where(up, iv.other_allele, iv.effect_allele)
rf = np.where(up, iv.eaf, 1 - iv.eaf)
pal = [{a, b} in ({"A", "T"}, {"C", "G"}) for a, b in zip(iv.effect_allele, iv.other_allele)]
# UKB-PPP p values below the double-precision range are stored at the floor (pval_underflow) with the true value in log10p:
# print -log10(p) from log10p for every instrument, after asserting the two columns agree wherever p did not underflow
uf = iv.pval_underflow.astype(bool)
assert (iv.loc[uf, "pval"] < 1e-300).all() and (iv.loc[~uf, "pval"] > 0).all(), "underflow flag does not match the stored p"
assert np.allclose(-np.log10(iv.loc[~uf, "pval"]), iv.loc[~uf, "log10p"], rtol=0, atol=1e-6), "log10p is not -log10(pval)"
assert (iv.log10p > 0).all()
gmin = ins.groupby("gene").pval.min()
gws = x.gene.map(gmin) < LV("L231")
eq(int((~gws).sum()), LV("L023"), "pairs without a genome-wide-significant instrument (L023)")
eq([int(gws.sum()), int(x[gws].gene.nunique())], LV("L025"), "pairs / proteins kept at 5e-8 (L025)")
eq(sorted((x[~gws].gene + "-" + x[~gws].outcome).tolist()), LV("L024")["pairs"], "the pairs without a genome-wide-significant instrument (L024)")
v = x[["gene", "outcome"]].merge(vep, on=["gene", "outcome"], how="left", validate="1:1")
assert (v.lead_snp == x.lead_snps).all()
mis = v.missense.map(tf)
eq([int(mis.sum()), int(v[mis].gene.nunique())], [LV("L020")["pairs"], LV("L020")["genes"]], "missense-led pairs / genes (L020)")
def _cons(mc, aa, m_, c, ms):
    if tf(mc):
        m2 = re.fullmatch(r"([A-Z*])/([A-Z*])(\d+)", aa)
        assert m2, aa
        return "missense %s%s%s (canonical transcript)" % (m2.group(1), m2.group(3), m2.group(2))
    if tf(m_):
        return "missense in a non-canonical transcript of the gene; canonical: %s" % c.replace("_", " ")
    if isinstance(c, str):
        return c.replace("_", " ").replace(";", "; ")
    return "outside the gene's transcripts; most severe consequence: %s" % ms.replace("_", " ")


cons = [_cons(*a) for a in zip(v.missense_canonical, v.amino_acids_canonical, v.missense, v.consequence_canonical, v.most_severe_site)]
SUSIE = {"supportive": "supportive", "inconclusive": "inconclusive", "discordant": "discordant", "not_evaluable": "not evaluable"}
SNA = {"susie_nonconvergence": " (no converged fit)", "MHC_excluded": " (MHC excluded)"}
sus = [SUSIE[c] + (SNA[r] if isinstance(r, str) else "") for c, r in zip(x.susie_class, x.susie_na_reason)]
SIG = {"Bonferroni": "study-wide (Bonferroni)", "FDR": FDRP + " FDR", "Suggestive": "nominal (p < %s)" % PNOM}
# S2-D47 (Stage 7 A2.2): the five proxy pairs. The record's pairs must be exactly the pairs the MR table flags (any_proxy), the count is the ledger's (L543), and the lead and the protein-raising allele it
# names must be the ones this sheet prints.
_key = x.gene + "-" + x.outcome
assert pxr.pair.is_unique
_px = pxr.set_index("pair")
eq(sorted(_key[x.any_proxy.map(tf)]), sorted(_px.index), "proxy pairs: record vs any_proxy of the MR table")
eq(len(_px), NPX, "proxy pairs (L543)")
for _i, _k in enumerate(_key):
    if _k in _px.index:
        assert x.lead_snps.iloc[_i] == _px.loc[_k].lead_rsid and ra[_i] == _px.loc[_k].protein_raising_allele and iv.chrpos.iloc[_i] == _px.loc[_k].lead_pos_grch38, ("proxy record disagrees with S3 on the lead", _k)


def _pcol(f):
    return [f(_px.loc[k]) if k in _px.index else "--" for k in _key]


s3 = pd.DataFrame({
    "Protein": x.gene, "Outcome": x.outcome, "Evidence tier": x.evidence_tier, "Weaknesses (Flagged pairs)": wk, "Reason named (one per Flagged pair)": named,
    "Internal criteria met (of %d)" % crit.shape[1]: crit.sum(axis=1).astype(int), "Screen significance": x.mr_sig_class.map(SIG),
    "MR estimate (log-odds per SD of protein; eGFR: eGFR*)": x.b_star, "MR SE": x.se, "MR p": x.pval,
    C_H4D: x.pp_h4, "PP.H3, default prior": x.pp_h3, C_H4C: h6,
    "Conservative-prior band": band, "SuSiE-based colocalization": sus,
    "Steiger direction": x.steiger_result.map({"FORWARD_CONFIRMED": "forward", "BORDERLINE": "non-forward"}),
    "deCODE": mapv(x.decode_replication_class, DEC, "deCODE class"), "ARIC": mapv(x.aric_replication_class, ARI, "ARIC class"),
    "Exposure-side corroboration": mapv(x.crossplatform_class, XP, "combined class"),
    "Lead variant": x.lead_snps, "Position (GRCh38)": iv.chrpos, "Protein-raising allele": ra, "Other allele": oa,
    "Protein-raising allele frequency": rf, "Instrument effect (SD of protein per protein-raising allele)": iv.beta.abs(),
    "Instrument SE": iv.se, "Instrument -log10(p)": iv.log10p, "Instrument F": iv.F_stat,
    C_GW: [yn(g) for g in gws], "Palindromic instrument": [yn(p) for p in pal],
    "LD proxy in the outcome GWAS": [yn(tf(a)) for a in x.any_proxy],
    "Outcome-side proxy": _pcol(lambda r: r.proxy_rsid), "Proxy position (GRCh38)": _pcol(lambda r: r.proxy_pos_grch38),
    "Proxy alleles in the outcome file (effect/other)": _pcol(lambda r: "%s/%s" % (r.proxy_effect_allele_in_outcome_file, r.proxy_other_allele_in_outcome_file)),
    "Proxy allele on the protein-raising haplotype": _pcol(lambda r: r.proxy_allele_on_protein_raising_haplotype),
    "Proxy r2 with the lead (1000 Genomes European)": _pcol(lambda r: float(r.r2_with_lead)),
    "Lead present in the outcome file": _pcol(lambda r: r.lead_in_outcome_file),
    "Reason for the proxy": _pcol(lambda r: r.proxy_reason),
    "Lead-variant consequence (Ensembl 116)": cons})
s3 = s3.sort_values(["Evidence tier", "Outcome", "Protein"], key=lambda c: c.map({"Convergent": 0, "Supported": 1, "Flagged": 2}) if c.name == "Evidence tier" else c).reset_index(drop=True)
sheet("S3", "Table S3. The %d colocalization-supported pairs." % N81,
      ("One row per protein–outcome pair (%d proteins), with the evidence tier and its parts. Internal support: screen significance, "
       "colocalization under the default and conservative priors, and SuSiE-based colocalization; every pair met the colocalization rule "
       "(PP.H4 ≥ %.2f and PP.H4/PP.H3 ≥ %g at p12 = %s). Exposure-side corroboration: the same lead variant in deCODE and ARIC (SomaScan), "
       "which reuse the discovery outcome data. Effects are per SD of protein: log-odds for binary outcomes; for eGFR, eGFR* = −ln(eGFR) per SD "
       "of protein, so positive values mean lower kidney function. Instrument columns describe the UKB-PPP lead cis variant, aligned to the "
       "protein-raising allele. For %d pairs the outcome association was read at the lead, so the instrument's association with the outcome equals its protein effect times the MR estimate. "
       "For the other %s, it was read at a high-LD proxy (proxy columns, with the reason the lead was not used) and aligned to the protein-raising allele through phased haplotypes, "
       "so the protein effect times the MR estimate is reconstructed, not observed at the lead.")
      % (N62, RULE_H4, RULE_RATIO, sci(PRIOR), N81 - NPX, spell(NPX)),
      s3, {
          "Protein": "gene symbol of the protein", "Outcome": "CAD, HF, AF, T2D, NAFLD, Stroke, eGFR, CKD", "Evidence tier": "Convergent / Supported / Flagged",
          "Weaknesses (Flagged pairs)": "every weakness of the pair; -- for none", "Reason named (one per Flagged pair)": "the first weakness in the order: discordant SuSiE, opposing arm, discordant assays, Steiger, conservative prior",
          "Internal criteria met (of %d)" % crit.shape[1]: "study-wide significance, robustness to the conservative prior, supportive SuSiE-based colocalization, colocalization rule",
          "Screen significance": "study-wide (Bonferroni over %s tests) / %s FDR / nominal" % (NFAM, FDRP),
          "MR estimate (log-odds per SD of protein; eGFR: eGFR*)": "Wald ratio per SD of protein; eGFR* = -ln(eGFR) per SD, positive = lower eGFR",
          "MR SE": "standard error of the MR estimate", "MR p": "p of the MR estimate",
          C_H4D: "coloc.abf posterior for a shared variant", "PP.H3, default prior": "coloc.abf posterior for distinct variants",
          C_H4C: "coloc.abf posterior for a shared variant under the conservative prior",
          "Conservative-prior band": "robust (>= %s) / %.2f to %s / below %.2f" % (T80, HALF, T80, HALF), "SuSiE-based colocalization": "supportive / inconclusive / not evaluable / discordant",
          "Steiger direction": "forward / non-forward", "deCODE": "result in deCODE (F >= %s required)" % F10, "ARIC": "result in ARIC (F >= %s required)" % F10,
          "Exposure-side corroboration": "deCODE and ARIC combined (labels as in Figure 2a)", "Lead variant": "rsID of the UKB-PPP lead cis variant",
          "Position (GRCh38)": "chromosome:position", "Protein-raising allele": "allele that raises the protein", "Other allele": "the other allele",
          "Protein-raising allele frequency": "frequency in UKB-PPP", "Instrument effect (SD of protein per protein-raising allele)": "UKB-PPP effect of the protein-raising allele, therefore positive",
          "Instrument SE": "standard error", "Instrument -log10(p)": "UKB-PPP p as -log10(p), read from the source's log10 p (for the strongest instruments p lies below the double-precision range and cannot be printed as p)", "Instrument F": "F-statistic",
          C_GW: "whether the protein's instrument reaches genome-wide significance", "Palindromic instrument": "A/T or C/G variant",
          "LD proxy in the outcome GWAS": "whether the discovery outcome GWAS needed a high-LD proxy",
          "Outcome-side proxy": "rsID of the proxy variant whose outcome association was used; -- if none", "Proxy position (GRCh38)": "chromosome:position of the proxy; -- if none",
          "Proxy alleles in the outcome file (effect/other)": "the proxy's effect and other allele as the outcome file lists them; -- if none",
          "Proxy allele on the protein-raising haplotype": "the proxy allele that travels with the protein-raising allele of the lead (phased 1000 Genomes European haplotypes); -- if none",
          "Proxy r2 with the lead (1000 Genomes European)": "r2 between lead and proxy; -- if none",
          "Lead present in the outcome file": "whether the outcome file lists the lead (with the label it carries there); -- if no proxy was used",
          "Reason for the proxy": "why the outcome association was not read at the lead: a palindromic lead whose allele frequency could not resolve its strand, a lead listed under another identifier, or a lead absent from the outcome file; -- if no proxy was used",
          "Lead-variant consequence (Ensembl 116)": "VEP consequence on the canonical transcript; missense in any transcript noted"},
      ["tiered", "mr", "ins", "vep", "proxy"],
      {"MR estimate (log-odds per SD of protein; eGFR: eGFR*)": "0.0000", "MR SE": "0.0000", "MR p": "0.00E+00", C_H4D: "0.000",
       "PP.H3, default prior": "0.000", C_H4C: "0.000", "Protein-raising allele frequency": "0.000",
       "Instrument effect (SD of protein per protein-raising allele)": "0.000", "Instrument SE": "0.0000", "Instrument -log10(p)": "0.0", "Instrument F": "0.0",
       "Proxy r2 with the lead (1000 Genomes European)": "0.000"})

# ================================================================================================ S4 development stage
dv, dsum = DEV.table()
eq(len(dv), N62, "S4 proteins (L011)")
L450 = LV("L450")
eq([dsum["approved"], dsum["clinical"], dsum["tractable"], dsum["none"], len(dsum["approved_opposite"]), len(dsum["approved_concordant"]),
    len(dsum["approved_unknown"]), dsum["approved_opposite"], dsum["approved_concordant"], dsum["approved_unknown"], dsum["restaged"],
    dsum["approved_withdrawn"]], L450,
   "development stage (L450)")
phw = DEV.phewas()
eq([phw["consensus_flagged"], phw["genes"], phw["ot_flagged"], phw["ot_queryable"]], LV("L451"), "PheWAS flags (L451)")
PHW = {"clean": "no flag", "pleiotropy_flagged": "pleiotropy flag", "safety_flagged": "safety flag", "not_queryable": "not queryable", "not_queried": "not queried"}
PHC = {"flag_multi_source": "flagged by two or more sources", "flag_single_source": "flagged by one source"}
pdx = DEV.pair_directions()


def dirtext(o, d):
    if o == "eGFR":
        return "higher protein, lower eGFR" if d == "risk" else "higher protein, higher eGFR"
    return "higher protein, higher risk" if d == "risk" else "higher protein, lower risk"


s4 = pd.DataFrame({
    "Protein": dv.gene,
    "Colocalization-supported pairs: outcome, genetic direction (tier)": [
        "; ".join("%s: %s (%s)" % (o, dirtext(o, d), t) for o, d, t in sorted(zip(pdx[pdx.gene == g].outcome, pdx[pdx.gene == g].mr_direction,
                                                                                   pdx[pdx.gene == g].evidence_tier))) for g in dv.gene],
    "Highest evidence tier": dv.best_tier, "Development stage (direct agents)": dv.stage, "Stage changed by the audit": dv.stage_changed,
    "Direct agent(s) at that stage": dv.agents, "ChEMBL action type(s)": dv.action_types.str.lower(), "Agent effect on the protein": dv.effect,
    "Mode of action vs genetic direction": dv.verdict, "UniProt subcellular location (audited proteins)": dv.uniprot_location,
    "Small-molecule tractable": dv.tract_sm, "Antibody tractable": dv.tract_ab,
    "PheWAS flag: Open Targets": mapv(dv.phewas_ot, PHW, "PheWAS OT"), "PheWAS flag: IEU OpenGWAS": mapv(dv.phewas_ieu, PHW, "PheWAS IEU"),
    "PheWAS flag: GWAS Atlas": mapv(dv.phewas_atlas, PHW, "PheWAS Atlas"), "PheWAS flag: any-source consensus": mapv(dv.phewas_consensus, PHC, "PheWAS consensus"),
    "Audit source": dv.audit_source, "Audit note": dv.audit_note})
s4 = s4.assign(_n=dv.stage_n.values).sort_values(["_n", "Protein"]).drop(columns="_n").reset_index(drop=True)
sheet("S4", "Table S4. Development stage and mode of action of the %d proteins." % N62,
      ("One row per protein. Development stage counts only agents that act directly on the single protein (a ChEMBL mechanism on a single-protein "
       "target): approved, in clinical development, tractable (a direct drug-gene interaction in DGIdb, or small-molecule or antibody tractability "
       "in Open Targets), or no known agent or tractability. Every protein with a direct approved or clinical-stage agent was audited in Open Targets on %s "
       "(Audit source), with one rule for every protein: an agent counts at the highest stage it reached in any jurisdiction, as Open Targets records "
       "it, and an approved agent whose entry for the protein Open Targets lists at the withdrawal stage is marked (FN1: ocriplasmin). COL15A1's agents "
       "act on a group of collagens, so its stage was re-assigned from tractability. Mode of action vs genetic "
       "direction compares the agent's effect on the protein with the sign of the genetic estimate per SD of circulating protein: concordant when the "
       "agent moves the protein in the direction the genetic estimate associates with lower risk, opposite when it moves it towards higher risk, and "
       "unknown when the agents act in both directions or the protein is a membrane receptor, whose plasma assay measures a circulating form that need "
       "not track the receptor signalling the agent acts on. A rule specified in advance demoted phenome-wide association flags to a "
       "supplementary annotation if at least %d%% of genes were flagged. It did not fire on the Open Targets flags (%d of %d genes); the any-source "
       "consensus, which flagged every protein (%d of %d), was chosen after both counts were known (Additional file 2), so the flags do not "
       "discriminate between proteins.") % (AUDIT_DATE, LV("L411")["stoploss_pct"], phw["ot_flagged"], phw["ot_queryable"], phw["consensus_flagged"], phw["genes"]),
      s4, {"Protein": "gene symbol", "Colocalization-supported pairs: outcome, genetic direction (tier)": "direction per SD of circulating protein (eGFR: direction of eGFR)",
           "Highest evidence tier": "strongest tier among the protein's pairs",
           "Development stage (direct agents)": "approved / in clinical development / tractable / no known agent or tractability; '(agent listed as withdrawn)' marks an approved agent whose entry for this protein Open Targets lists at the withdrawal stage",
           "Stage changed by the audit": "no / yes: pipeline stage -> audited stage",
           "Direct agent(s) at that stage": "agents acting on this protein alone (ChEMBL via Open Targets); -- for none",
           "ChEMBL action type(s)": "action type of the agent(s) on the protein", "Agent effect on the protein": "lowers / raises / mixed activity of the protein",
           "Mode of action vs genetic direction": "concordant / opposite / unknown / no direct agent / no agent",
           "UniProt subcellular location (audited proteins)": "as returned by Open Targets; a single-pass type I membrane protein is read as a membrane receptor",
           "Small-molecule tractable": "Open Targets tractability", "Antibody tractable": "Open Targets tractability",
           "PheWAS flag: Open Targets": "no flag / pleiotropy flag / safety flag / not queryable", "PheWAS flag: IEU OpenGWAS": "as for Open Targets",
           "PheWAS flag: GWAS Atlas": "as for Open Targets; queried for a subset", "PheWAS flag: any-source consensus": "flagged by one / two or more sources",
           "Audit source": "where and when the agent and its mechanism were checked", "Audit note": "what the audit found"},
      ["master", "tiered", "audit", "snap", "oldt4"])

# ================================================================================================ S5 external outcome checks
oc = pd.read_csv(SRC["oc"], sep="\t")
ph = pd.read_csv(SRC["ph"], sep="\t")
ORD6 = ["PAM-T2D", "RAD23B-T2D", "ARG1-T2D", "PCSK9-CAD", "TGFB1-CAD", "FGFR4-T2D (EAS)"]
eq(oc.pair.tolist()[:6], ORD6, "order of the pre-specified checks")
eq(oc.pair.tolist()[6], "ARG1-T2D signal B", "seventh row")
six = oc.iloc[:6]
eq(len(six), LV("L268"), "pairs selected in advance (L268)")
eq([int((six["class"] == "concordant").sum()), int(six.bonferroni_6.astype(bool).sum())], LV("L048"), "concordant / Bonferroni (L048)")
# Stage 6.3b: L532 is the successor of L048 ([concordant, pairs checked, Bonferroni] = [6, 6, 4]). This sheet prints neither the
# pairs-checked leaf nor anything else the successor changed (the pairs checked are asserted against L268 above), and L532's other two
# leaves equal L048's, so the predecessor keeps being read; the successor is asserted equal here so that a divergence fails the build.
eq([LV("L532")[0], LV("L532")[2]], LV("L048"), "L532 (successor of L048): concordant / Bonferroni leaves")
eq(LV("L532")[1], len(six), "L532 (successor of L048): pairs checked")
eq([float(p) for p in six.p_out_file] + [int((six.p_out_file < LV("L961")[2]).sum())], LV("L960")[:7], "resource-reported p (L960)")
for (pair, o_r, lo_r, hi_r, p_r), (_, r) in zip(LV("L223"), six.iterrows()):
    assert pair == r.pair and abs(o_r - np.exp(r.wald)) < 1e-12 and abs(lo_r - np.exp(r.wald - 1.96 * r.wald_se)) < 1e-12, ("L223", pair)
P05, BON6 = LV("L961")[0], LV("L961")[2]


def classify(est, se_, p, disc):
    same = np.sign(est) == np.sign(disc)
    return "concordant" if same and p < P05 else ("directionally concordant" if same else ("discordant" if p < P05 else "inconclusive"))


# the rule, re-applied here to the reported p, must give the stored class (Methods: classification unchanged by the p source)
for _, r in oc.iterrows():
    eq(classify(r.wald, r.wald_se, r.p_out_file, r.discovery_wald), r["class"], "class of %s" % r.pair)
for _, r in ph.iterrows():
    eq(classify(r.wald, r.se, r.p, r.disc_b), r.cls if not r.cls.startswith("not evaluable") else "?", "class of %s-%s" % (r.gene, r.outcome))
RES = {"DIAMANTE": "DIAMANTE 2018, European, without UK Biobank", "CARDIOGRAM": "CARDIoGRAMplusC4D 2015", "TPMI": "TPMI PheWeb, PheCode 250.2"}
tierof = dict(zip(pe.gene + "-" + pe.outcome, pe.evidence_tier))
DIA_NOTE = "DIAMANTE reports beta and SE to two significant figures, so the Wald p and the heterogeneity p are approximate."
# the allele letter the discovery frequency refers to: the instrument's effect allele (the checks table carries it as `ea`; the post hoc table does not, so it is read from the instruments file and the
# post hoc table's discovery frequency is asserted to be that allele's frequency). Resource frequencies are aligned to the same allele in both scripts (S1R1_02_outcome_checks.R, P31R_posthoc_check.py).
_ins_ea = {(g, s): (e, f) for g, s, e, f in zip(ins.gene, ins.SNP, ins.effect_allele, ins.eaf)}
assert all(re.fullmatch(r"[ACGT]", a) for a in oc.ea), "a discovery effect allele is not a single base"
for _, r in oc.iterrows():
    _k = (r.pair.split("-")[0], r.instrument)
    if _k in _ins_ea:
        assert _ins_ea[_k][0] == r.ea and abs(_ins_ea[_k][1] - r.eaf_exposure) < 6e-7, ("S5 discovery allele or frequency differs from the instruments file", r.pair)
# Stage 7 A3-fix F8: the two rows the instruments-file guard above never reaches (FGFR4 East Asian; ARG1 signal B) are asserted against their source records, read from file
assert md5(SRC["ckbf"]) == "15d2d38c1a34b90b7708239f197b0a55" and md5(SRC["harm"]) == "c394b6a1c0acbe32bfabc65d136ac974", "an S5 discovery source file changed"
_ck = pd.read_csv(SRC["ckbf"], sep="\t", dtype=str)
_ck = _ck[_ck.rsid == "rs451643"]
assert len(_ck) == 1, "rs451643 is not exactly one row of the CKB FGFR4 file"
_hm = pd.read_csv(SRC["harm"], sep="\t", usecols=["rsid_ref", "a1", "a2", "af_ref"])
_hm = _hm[_hm.rsid_ref == "rs17788484"]
assert len(_hm) == 1, "rs17788484 is not exactly one row of the ARG1 harmonised panel"
_f4 = oc[oc.pair == "FGFR4-T2D (EAS)"].iloc[0]
_sb = oc[oc.pair == "ARG1-T2D signal B"].iloc[0]
assert (_f4.instrument, _f4.ea) == ("rs451643", _ck.alt.iloc[0]) and _ck.ref.iloc[0] == "G", "S5 FGFR4 allele differs from the CKB file"
assert abs(_f4.eaf_exposure - float(_ck.af.iloc[0])) < 1e-9, "S5 FGFR4 discovery frequency differs from the CKB file"
assert (_sb.instrument, _sb.ea) == ("rs17788484", _hm.a1.iloc[0]) and _hm.a2.iloc[0] == "C", "S5 signal B allele differs from the harmonised panel"
_hA = pd.read_csv(SRC["harm"], sep="\t", usecols=["rsid_ref", "a1", "af_ref"])
_hA = _hA[_hA.rsid_ref == "rs2781668"]
# orientation anchor: at signal A the panel's a1 and af_ref agree with the UKB-PPP instruments file, so af_ref is the frequency of a1
assert len(_hA) == 1 and _hA.a1.iloc[0] == _ins_ea[("ARG1", "rs2781668")][0] and abs(float(_hA.af_ref.iloc[0]) - _ins_ea[("ARG1", "rs2781668")][1]) < 0.01, "af_ref is not the frequency of a1 at ARG1 signal A"
assert abs(_sb.eaf_exposure - float(_hm.af_ref.iloc[0])) < 5.1e-5, "S5 signal B frequency differs from the harmonised panel's reference frequency"
for _, r in ph.iterrows():
    assert (r.gene, r.instrument) in _ins_ea and abs(_ins_ea[(r.gene, r.instrument)][1] - r.eaf_exp) < 1e-9, ("post hoc discovery frequency is not the instrument's effect-allele frequency", r.gene)
rows = []
for _, r in oc.iterrows():
    key = r.pair.replace(" (EAS)", "").replace(" signal B", "")
    note = DIA_NOTE if r.resource == "DIAMANTE" else ("TPMI's p is SAIGE's score-test p. The CKB protein effect is per SD of rank-inverse-normal protein level and cross-ancestry magnitudes are approximate, so this pair is compared by sign and p."
                                                       if r.resource == "TPMI" else "--")
    if "signal B" in r.pair:
        note = ("Signal B (rs17788484) of the arginase-1 analysis, not one of the six original external checks. For this signal, the exposure-side "
                "allele frequency comes from the European reference panel because the exposure-sample frequency was unavailable for this SNP. ") + note
    rows.append({"Analysis": "original external check" if "signal B" not in r.pair else "arginase-1 signal B",
                 "Pair": dash(key) + (" (East Asian)" if "EAS" in r.pair else ""), "Original external check": yn("signal B" not in r.pair),
                 "Evidence tier": tierof.get(key, "East Asian only"), "External resource": RES[r.resource], "Instrument": r.instrument,
                 "Palindromic instrument": r.palindrome if isinstance(r.palindrome, str) else "no",
                 "Effect allele, discovery": r.ea,
                 "Effect-allele frequency, discovery": r.eaf_exposure, "Effect-allele frequency, resource": r.eaf_resource_aligned,
                 "External estimate (log-odds per SD of protein)": r.wald, "External SE": r.wald_se, "Odds ratio": np.exp(r.wald),
                 "OR 95% CI lower": np.exp(r.wald - 1.96 * r.wald_se), "OR 95% CI upper": np.exp(r.wald + 1.96 * r.wald_se),
                 "p, as reported by the resource": r.p_out_file, "p, Wald ratio from beta and SE": r.p, "Classification": r["class"],
                 C_BON: yn(bool(r.bonferroni_6)) if "signal B" not in r.pair else "not applicable",
                 "Discovery estimate (log-odds per SD of protein)": r.discovery_wald, "Discovery SE": r.discovery_se,
                 # S2-D43 (M5): the TPMI pair is compared by sign and p only (its note, Table 2 and the Results), so no magnitude test is printed for it
                 "Heterogeneity p, external vs discovery": "not reported (sign and p only)" if r.resource == "TPMI" else r.het_p, "Note": note})
assert sum(1 for x in rows if x["Heterogeneity p, external vs discovery"] == "not reported (sign and p only)") == 1, "the sign-and-p-only pair is not exactly one row"
eq(len(ph), LV("L250")[0], "post hoc pairs (L250)")
phe = ph[~ph.cls.str.startswith("not evaluable")]
eq([len(phe), int(phe.same_sign.sum()), int((phe.cls == "concordant").sum())], LV("L250")[:3], "post hoc same sign / concordant (L250)")
phn = phe[~phe.predeclared]
eq([len(phn), int(phn.same_sign.sum()), int((phn.cls == "concordant").sum())], LV("L251")[:3], "post hoc, not selected in advance (L251)")
eq([[a, b, c, float(d)] for a, b, c, d in phe[~phe.same_sign][["gene", "outcome", "tier", "p"]].values], LV("L253"), "opposite sign (L253)")
eq([int((phe.het_p < LV("L961")[3]).sum()), len(phe)], LV("L254"), "post hoc size differences (L254)")
# cross-file identity: the five European pairs selected in advance carry the same numbers in both files
for _, r in ph[ph.predeclared].iterrows():
    o = oc.set_index("pair").loc["%s-%s" % (r.gene, r.outcome)]
    assert abs(o.wald - r.wald) < 1e-12 and abs(o.wald_se - r.se) < 1e-12 and abs(o.p_out_file - r.p) < 1e-15 and abs(o.p - r.p_wald) < 1e-15, (r.gene, "post hoc vs pre-specified")
eq(int(ph.predeclared.sum()), LV("L268") - 1, "European pairs selected in advance inside the post hoc check")
for _, r in ph.iterrows():
    res = "DIAMANTE" if r.outcome == "T2D" else "CARDIOGRAM"
    rows.append({"Analysis": "post hoc (all T2D and CAD pairs)", "Pair": "%s–%s" % (r.gene, r.outcome), "Original external check": yn(bool(r.predeclared)),
                 "Evidence tier": r.tier, "External resource": RES[res], "Instrument": r.instrument,
                 "Palindromic instrument": r.palindrome if isinstance(r.palindrome, str) else "no",
                 "Effect allele, discovery": _ins_ea[(r.gene, r.instrument)][0],
                 "Effect-allele frequency, discovery": r.eaf_exp, "Effect-allele frequency, resource": r.eaf_res,
                 "External estimate (log-odds per SD of protein)": r.wald, "External SE": r.se, "Odds ratio": r.OR, "OR 95% CI lower": r.OR_lo,
                 "OR 95% CI upper": r.OR_hi, "p, as reported by the resource": r.p, "p, Wald ratio from beta and SE": r.p_wald, "Classification": r.cls,
                 C_BON: "not applicable", "Discovery estimate (log-odds per SD of protein)": r.disc_b, "Discovery SE": r.disc_se,
                 "Heterogeneity p, external vs discovery": r.het_p, "Note": DIA_NOTE if res == "DIAMANTE" else "--"})
s5 = pd.DataFrame(rows)
assert set(s5.Classification) <= {"concordant", "directionally concordant", "inconclusive", "discordant"}
nn = {6: "six"}
sheet("S5", "Table S5. External outcome checks.",
      ("The %s pairs of the original external check, with signal B of the arginase-1 analysis, and the post hoc check of all %d colocalization-supported T2D "
       "and CAD pairs, which repeats the European pairs of the original check. Each check used the discovery instrument without proxies; estimates "
       "are log-odds per SD of protein. p is each resource's own p value (for TPMI, SAIGE's score-test p), with the Wald-ratio p from beta and SE "
       "beside it; the change of p source was made after the checks were run, and the resource p values and the first-order Wald p values gave the same concordance classifications. For RAD23B–T2D, rounding of the published DIAMANTE estimates could move the Wald p value across the concordance threshold (main text, Table 3). Concordant: discovery sign with p < %g; "
       "directionally concordant: discovery sign with p ≥ %g; discordant: opposite sign with p < %g; inconclusive: otherwise. Bonferroni for the "
       "original external checks: p < %s. Heterogeneity p compares the external with the discovery estimate (z-test). Frequencies are those of the allele in the column Effect allele, discovery; "
       "the Wald ratios do not depend on which allele is used. The discovery frequency is that of the discovery sample (UK Biobank Pharma Proteomics Project, or the China Kadoorie Biobank for FGFR4). For arginase-1 signal B it comes from the European reference panel, because the exposure-sample frequency was unavailable. The European resources "
       "include Finnish studies and may share participants with FinnGen, so these are external outcome checks rather than independent replication. "
       "Their sampling-noise correlation with FinnGen (the correlation of null z-scores implied by the cross-trait LD-score-regression intercept) was %.3f (SE %.3f) for T2D and %.3f (SE %.3f) "
       "for CAD (Additional file 2).")
      % ((nn[LV("L268")], LV("L250")[0], P05, P05, P05, BON) + tuple(LV("L449"))),
      s5, {"Analysis": "original external check / arginase-1 signal B / post hoc", "Pair": "protein-outcome", "Original external check": "yes / no",
           "Evidence tier": "tier of the pair (FGFR4-T2D: found in the East Asian screen only)", "External resource": "external outcome GWAS",
           "Instrument": "discovery lead variant (no proxies)", "Palindromic instrument": "A/T or C/G, with whether allele frequency orients it",
           "Effect allele, discovery": "the allele that the two frequency columns refer to (the instrument's effect allele in the exposure data)",
           "Effect-allele frequency, discovery": "frequency of the allele in 'Effect allele, discovery', in the exposure data (arginase-1 signal B: European reference panel)",
           "Effect-allele frequency, resource": "frequency of the same allele in the external resource",
           "External estimate (log-odds per SD of protein)": "Wald ratio in the external resource", "External SE": "first-order SE",
           "Odds ratio": "exp(estimate)", "OR 95% CI lower": "exp(estimate - 1.96 SE)", "OR 95% CI upper": "exp(estimate + 1.96 SE)",
           "p, as reported by the resource": "the resource's own p (TPMI: SAIGE score test)", "p, Wald ratio from beta and SE": "for comparison",
           "Classification": "concordant / directionally concordant / inconclusive / discordant", C_BON: "original external checks only",
           "Discovery estimate (log-odds per SD of protein)": "FinnGen (FGFR4-T2D: East Asian discovery estimate)", "Discovery SE": "SE of the discovery estimate",
           "Heterogeneity p, external vs discovery": "z-test of the difference", "Note": "precision and scale caveats"},
      ["oc", "ph", "ckbf", "harm"],
      {c: "0.0000" for c in ["Effect-allele frequency, discovery", "Effect-allele frequency, resource", "External estimate (log-odds per SD of protein)",
                             "External SE", "Discovery estimate (log-odds per SD of protein)", "Discovery SE"]} |
      {c: "0.000" for c in ["Odds ratio", "OR 95% CI lower", "OR 95% CI upper"]} |
      {c: "0.00E+00" for c in ["p, as reported by the resource", "p, Wald ratio from beta and SE", "Heterogeneity p, external vs discovery"]})

# ================================================================================================ S6 multi-signal candidates
ca = pd.read_csv(SRC["cand"])
L452 = LV("L452")
eq([len(ca)] + [int((ca.Status == s).sum()) for s in ("Drawn", "Held (not drawn)", "Dropped")], L452[1:5], "candidate audit (L452)")
eq(L452[2:5], LV("L222"), "drawn / held / dropped (L222)")
eq(L452[0], LV("L263")[3], "SuSiE-testable pool (L263)")
eq([int(ca[ca.Status == "Drawn"].Endpoint.isin(["eGFR", "CKD"]).sum()), int((ca.Status == "Drawn").sum())], LV("L225")[:2], "kidney candidates (L225)")
# the direction label follows the MR estimate: risk = higher protein with higher risk (binary) or lower ln(eGFR)
cm = ca.merge(prim[["gene", "outcome_name", "b"]], left_on=["Gene", "Endpoint"], right_on=["gene", "outcome_name"], how="left", validate="1:1")
assert cm.b.notna().all()
assert ((cm["MR direction"] == "risk") == np.where(cm.Endpoint == "eGFR", cm.b < 0, cm.b > 0)).all(), "candidate MR direction convention"
_it = sorted({int(v) for v in ca.Reason.str.extract(r"max_iter=(\d+)")[0].dropna()})
_pc = sorted({float(v) for v in ca.Reason.str.extract(r"min_abs_corr=([0-9.]+)")[0].dropna()})
assert len(_it) == 1 and _pc == [LV("L247")], (_it, _pc)
REASON = {"PP.H4_SuSiE >= %s with a converged fit" % T80: "SuSiE PP.H4 ≥ %s with a converged fit" % T80,
          "PP.H4_SuSiE >= %s only with a non-convergent SuSiE fit (IBSS did not converge even at max_iter=%d)" % (T80, _it[0]):
              "SuSiE PP.H4 ≥ %s only with a fit that did not converge (IBSS, at most %s iterations)" % (T80, "{:,}".format(_it[0])),
          "PP.H4_SuSiE < %s under purity filter" % T80: "SuSiE PP.H4 < %s after purity filtering" % T80,
          "Collapsed to no pure credible set under purity filter (min_abs_corr=%s)" % PUR: "no credible set passed the purity filter (minimum absolute correlation %s)" % PUR}
STAT = {"Drawn": "kept (multi-signal candidate)", "Held (not drawn)": "not kept: the fit did not converge", "Dropped": "dropped"}   # house term 'kept' (3b critic F16)
pur = ca["SuSiE PP.H4 (purity %s)" % PUR].map(lambda s: "no pure credible set" if str(s).strip() == "—" else float(s))
s6 = pd.DataFrame({"Protein": ca.Gene, "Outcome": ca.Endpoint, "Kidney outcome": [yn(o in ("eGFR", "CKD")) for o in ca.Endpoint],
                   "MR direction": [dirtext(o, d) for o, d in zip(ca.Endpoint, ca["MR direction"])],
                   "coloc.abf PP.H4": ca["coloc.abf PP.H4"], "coloc.abf PP.H3": ca["coloc.abf PP.H3"],
                   "SuSiE PP.H4 before purity filtering": ca["SuSiE PP.H4 (naive)"], C_PUR: pur,
                   "SuSiE fit converged": ca["SuSiE converged (no warning)"], "FDR q < %g" % LV("L239"): ca["FDR q<%g" % LV("L239")],
                   "Status": mapv(ca.Status, STAT, "candidate status"), "Reason": mapv(ca.Reason, REASON, "candidate reason")})
sheet("S6", "Table S6. Multi-signal candidate associations.",
      ("Associations that passed the %s FDR on second-order Wald-ratio standard errors but failed coloc.abf, outside the major histocompatibility "
       "complex, were tested with SuSiE-based colocalization (%d SuSiE-testable associations; LD from the 1000 Genomes European panel). The %d that "
       "reached SuSiE PP.H4 ≥ %.2f before purity filtering are listed. They were re-run with credible-set purity %g: %d were kept as multi-signal "
       "candidates (PP.H4 ≥ %.2f with a converged fit; %d of them kidney pairs), %d was not kept because its fit did not converge, and %d were dropped. "
       "Candidates are shown separately from the colocalization-supported pairs and excluded from all counts. Because the LD comes from a reference "
       "panel, SuSiE results are read qualitatively.")
      % (FDRP, L452[0], L452[1], RULE_H4, LV("L247"), L452[2], RULE_H4, L452[5], L452[3], L452[4]),
      s6, {"Protein": "gene symbol", "Outcome": "outcome", "Kidney outcome": "eGFR or CKD",
           "MR direction": "direction per SD of protein (eGFR: direction of eGFR)", "coloc.abf PP.H4": "single-variant posterior, shared variant",
           "coloc.abf PP.H3": "single-variant posterior, distinct variants", "SuSiE PP.H4 before purity filtering": "SuSiE-based colocalization posterior",
           C_PUR: "after credible-set purity filtering", "SuSiE fit converged": "yes / no",
           "FDR q < %g" % LV("L239"): "on second-order standard errors", "Status": "kept (multi-signal candidate) / not kept: the fit did not converge / dropped", "Reason": "why"},
      ["cand", "mr"], {c: "0.0000" for c in ["SuSiE PP.H4 before purity filtering"]} | {c: "0.00E+00" for c in ["coloc.abf PP.H4", "coloc.abf PP.H3"]})

# ================================================================================================ S7a East Asian outcome-side
ea = pd.read_csv(SRC["eas"])
eq(len(ea), N81, "S7a rows (L010)")
ev = ea[ea.qc_evaluable.map(tf)]
st = ea.eas_outcome_state.value_counts().to_dict()
eq([len(ev), st.get("EAS_COMPATIBLE", 0), st.get("EAS_INCONCLUSIVE", 0), st.get("EAS_DISCORDANT", 0), st.get("EAS_NOT_EVALUABLE", 0)], LV("L080"),
   "outcome-side states (L080)")
sup_, fla_ = ev[ev.evidence_tier == "Supported"], ev[ev.evidence_tier == "Flagged"]
eq([int((ev.direction_class == "SAME").sum()), len(ev), int((sup_.direction_class == "SAME").sum()), len(sup_), int((fla_.direction_class == "SAME").sum()), len(fla_)],
   LV("L224"), "kept the European sign (L224)")
sts = ea.eas_outcome_state_strict.value_counts().to_dict()
eq([sts.get("EAS_COMPATIBLE", 0) + sts.get("EAS_INCONCLUSIVE", 0), sts.get("EAS_COMPATIBLE", 0), sts.get("STRICT_EXCLUDED", 0)], LV("L081"), "strict subset (L081)")
assert (ea[ea.evidence_tier == "Convergent"].eas_outcome_state == "EAS_NOT_EVALUABLE").all(), "a Convergent pair is evaluable outcome-side"
assert ((ev.eas_outcome_state == "EAS_COMPATIBLE") == ((ev.q_eas < 0.05) & (ev.direction_class == "SAME"))).all(), "compatible != FDR-significant with the European sign"
# sign orientation: European estimates are eGFR*, and the East Asian eGFR estimate is eGFR*-oriented (kidney sentinel: eGFR* and CKD agree)
xe = ev.merge(pe[["gene", "outcome", "mr_beta"]], on=["gene", "outcome"])
assert ((np.sign(xe.beta_mr_eas) == np.sign(xe.mr_beta)) == (xe.direction_class == "SAME")).all(), "direction class disagrees with the signs"
kd = ev[ev.gene.isin(["FGF5", "SDCCAG8"]) & ev.outcome.isin(["eGFR", "CKD"])].pivot(index="gene", columns="outcome", values="beta_mr_eas")
assert len(kd) == 2 and (np.sign(kd.eGFR) == np.sign(kd.CKD)).all(), "kidney sentinel failed: East Asian eGFR is not eGFR*-oriented"
EASRES = {lab: s2_name(lab, ["East Asian outcome-side"]) for lab in sorted(ea.eas_source.dropna().unique())}
QCR = {"VARIANT_ABSENT": "variant absent", "LOW_INFO": "low imputation information", "PALINDROMIC_AMBIGUOUS": "palindromic, strand ambiguous",
       "OUTCOME_UNAVAILABLE": "no East Asian outcome GWAS"}
SCB = {"VARIANT_LEVEL": "per variant", "STUDY_LEVEL_APPROX": "study level (approximate)", "UNRESOLVED": "unresolved"}
STT = {"EAS_COMPATIBLE": "compatible", "EAS_INCONCLUSIVE": "not significant", "EAS_NOT_EVALUABLE": "not evaluable", "EAS_DISCORDANT": "discordant",
       "STRICT_EXCLUDED": "excluded (study-level sample sizes)"}
DIRC = {"SAME": "same as European", "OPPOSITE": "opposite to European", "NULL_OR_UNDEFINED": "not evaluable"}
na_ = lambda s: s.map(lambda v: v if pd.notna(v) else "not evaluable")
s7a = pd.DataFrame({"Protein": ea.gene, "Outcome": ea.outcome, "Evidence tier": ea.evidence_tier, "European instrument": ea.instrument_rsID,
                    "East Asian outcome resource": ea.eas_source.map(lambda s: EASRES[s] if pd.notna(s) else "no East Asian outcome GWAS"),
                    "Evaluable": [yn(tf(v)) for v in ea.qc_evaluable],
                    "Reason not evaluable": ea.qc_reason.map(lambda v: QCR[v] if pd.notna(v) else "--"),
                    "Sample-size basis": ea.sample_count_basis.map(SCB),
                    "East Asian estimate (log-odds per SD of protein; eGFR: eGFR*-oriented)": na_(ea.beta_mr_eas), "SE": na_(ea.se_mr_eas),
                    "95% CI lower": na_(ea.ci95_low), "95% CI upper": na_(ea.ci95_high), "p": na_(ea.p_eas), "FDR q": na_(ea.q_eas),
                    "Direction": mapv(ea.direction_class, DIRC, "direction"), "Outcome-side state": mapv(ea.eas_outcome_state, STT, "state"),
                    "State in the strict sensitivity analysis": mapv(ea.eas_outcome_state_strict, STT, "strict state")})
assert set(ea.qc_reason.dropna()) <= set(QCR) and set(ea.sample_count_basis) <= set(SCB)
# S2-D43 (X1): hold H9 withholds the East Asian eGFR magnitudes of the EGFR_EAS_2P2M file (Additional file 2, Section S-R13 says so, and the release masks
# eas_beta_mr, eas_se_mr and eas_q of these rows). Their estimate, SE, CI and FDR q are withheld here too; p, direction and states stay, because the text
# compares this file with the European estimates by sign and p.
H9SRC, H9MARK = "EGFR_EAS_2P2M", "withheld (Additional file 2, Section S-R13)"
H9COLS = ["East Asian estimate (log-odds per SD of protein; eGFR: eGFR*-oriented)", "SE", "95% CI lower", "95% CI upper", "FDR q"]
h9 = ((ea.eas_source == H9SRC) & ea.qc_evaluable.map(tf)).values
_l536 = LV("L536")
eq([int(h9.sum()), int((ea[h9].eas_outcome_state == "EAS_COMPATIBLE").sum()), int((ea[h9].eas_outcome_state == "EAS_INCONCLUSIVE").sum())],
   [_l536["evaluable"], _l536["compatible"], _l536["inconclusive"]], "H9 rows of S7a (L536)")
assert set(ea[h9].outcome) == {"eGFR"}, "an H9 row is not an eGFR row"
for c in H9COLS:
    s7a[c] = s7a[c].astype(object)
    s7a.loc[h9, c] = H9MARK
assert all((s7a.loc[h9, c] == H9MARK).all() for c in H9COLS) and not any((s7a.loc[~h9, c] == H9MARK).any() for c in H9COLS), "H9 masking is not exact"
assert pd.to_numeric(s7a.loc[h9, "p"], errors="coerce").notna().all(), "the p of an H9 row is not printed"
nev = LV("L080")
sheet("S7a", "Table S7a. East Asian outcome-side analysis: European instrument, East Asian outcome GWAS.",
      ("All %d colocalization-supported pairs; %d were evaluable. A pair was compatible when its FDR-adjusted East Asian association was significant "
       "with the European sign, discordant when significant with the opposite sign, and not significant otherwise. Estimates are log-odds per SD of "
       "protein for binary outcomes; eGFR estimates are oriented as eGFR* (positive = lower eGFR). The East Asian eGFR file comes from a "
       "meta-analysis of z-scores (Table S2), so its effects are not ln(eGFR) effects and eGFR is compared with the European estimate by sign. "
       "For the %d pairs tested in that file, the estimate, SE, interval and q are withheld (Additional file 2, Section S-R13); their p, direction and state are shown. "
       "The strict sensitivity analysis keeps the outcome files that report per-variant sample "
       "sizes (T2D and eGFR).") % (N81, nev[0], int(h9.sum())),
      s7a, {"Protein": "gene symbol", "Outcome": "outcome", "Evidence tier": "European tier", "European instrument": "UKB-PPP lead variant",
            "East Asian outcome resource": "outcome GWAS used, named as in Table S2",
            "Evaluable": "yes / no", "Reason not evaluable": "why not", "Sample-size basis": "how the file reports sample size",
            "East Asian estimate (log-odds per SD of protein; eGFR: eGFR*-oriented)": "Wald ratio with the European instrument",
            "SE": "standard error", "95% CI lower": "lower limit", "95% CI upper": "upper limit", "p": "p", "FDR q": "Benjamini-Hochberg q",
            "Direction": "sign relative to the European estimate", "Outcome-side state": "compatible / not significant / discordant / not evaluable",
            "State in the strict sensitivity analysis": "as Outcome-side state, T2D and eGFR files only"},
      ["eas", "tiered"], {"East Asian estimate (log-odds per SD of protein; eGFR: eGFR*-oriented)": "0.0000", "SE": "0.0000", "95% CI lower": "0.0000",
                          "95% CI upper": "0.0000", "p": "0.00E+00", "FDR q": "0.00E+00"})

# ================================================================================================ S7b ancestry-matched
cen = pd.read_csv(SRC["cen"], sep="\t")
ckf = pd.read_csv(SRC["ckf"], sep="\t")
eq([len(cen), cen.gene.nunique()], LV("L216"), "ancestry-matched pairs / proteins (L216)")
tst = cen[cen.coloc_status == "ok"]
sup = tst.support.map(tf); rob = sup & tst["support_p12_1e-6"].map(tf)
eq([len(cen), len(tst), int(sup.sum()), int(rob.sum())], LV("L082"), "ancestry-matched (L082)")
eq([int((tst.call == "NOT_SUPPORTED_H1").sum()), int((tst.call == "NOT_SUPPORTED_H3").sum())], LV("L211"), "no association / distinct variants (L211)")
eq(round(float(tst.loc[(tst.gene == "DNER") & (tst.outcome == "eGFR"), "PP.H4"].iloc[0]), 3), LV("L255"), "DNER-eGFR (L255)")
eq(round(float(tst.F.min()), 1), LV("L221"), "weakest CKB instrument F (L221)")
# eGFR orientation: the census builder reverses the sign of the rank-inverse-normal eGFR outcome (flip=(outcome == "eGFR")), so a
# protein with a nominal eGFR estimate must carry the same sign on CKD (kidney sentinel, as in S7a)
_kd = cen[cen.outcome.isin(["eGFR", "CKD"])].pivot_table(index="gene", columns="outcome", values=["mr_beta", "mr_p"]).dropna()
_kd = _kd[_kd[("mr_p", "eGFR")] < LV("L961")[0]]
assert len(_kd) >= 1 and (np.sign(_kd[("mr_beta", "eGFR")]) == np.sign(_kd[("mr_beta", "CKD")])).all(), "kidney sentinel failed in S7b: eGFR is not risk-oriented"
# Stage 6 (closure row 9): the PheWeb.jp CAD file is not a hum0197 release, and Table S2 now names it PHEWEB_BBJ_CAD. The pinned census results
# keep its run-time label BBJ_CAD_hum0197.v5 (kept for traceability, Stage 6.3b fix record); S7b names the resource as Table S2 does.
S7B_LEGACY = {"BBJ_CAD_hum0197.v5": "PHEWEB_BBJ_CAD"}
assert all(set(pd.concat([cen, ckf]).query("accession == @a").outcome) == {"CAD"} for a in S7B_LEGACY), "a legacy S7b label names more than the CAD file"
S7BRES = {a: s2_name(S7B_LEGACY.get(a, a), ["Ancestry-matched", "East Asian outcome-side"]) for a in sorted(set(cen.accession.dropna()) | set(ckf.accession.dropna()))}
# CKB-free: replace every outcome that includes CKB with its Biobank Japan file
base = cen[cen.ckb_overlap.map(lambda v: str(v).strip().upper() == "FALSE")]
assert set(zip(ckf.gene, ckf.outcome)) == set(zip(cen[cen.ckb_overlap.map(lambda v: str(v).strip().upper() == "TRUE")].gene,
                                                  cen[cen.ckb_overlap.map(lambda v: str(v).strip().upper() == "TRUE")].outcome)), "CKB-free rows != CKB-including rows"
cf = pd.concat([base, ckf])
eq([int(cf.support.map(tf).sum()), int((cf.support.map(tf) & cf["support_p12_1e-6"].map(tf)).sum())], LV("L084"), "CKB-free colocalizations (L084)")
CALL = {"SUPPORTED": "colocalized", "NOT_SUPPORTED_H1": "protein-only model favored (H1 most probable)", "NOT_SUPPORTED_H3": "distinct variants favored (H3 most probable)",
        "NOT_SUPPORTED_H4_below_rule": "below the colocalization rule", "NOT_TESTABLE_NO_EAS_GWAS": "not testable (no East Asian outcome GWAS)",
        "NOT_SUPPORTED_H0": "no association with either trait (H0 most probable)", "NOT_SUPPORTED_H2": "no protein association (H2 most probable)"}
tierof2 = dict(zip(pe.gene + "|" + pe.outcome, pe.evidence_tier))


def s7b_rows(d, label):
    out = pd.DataFrame({"Analysis": label, "Protein": d.gene.values, "Outcome": d.outcome.values,
                        "Evidence tier (European)": [tierof2[g + "|" + o] for g, o in zip(d.gene, d.outcome)],
                        "East Asian outcome resource": d.accession.map(lambda a: S7BRES[a] if pd.notna(a) else "no East Asian outcome GWAS").values,
                        "Outcome resource includes CKB": [yn(str(v).strip().upper() == "TRUE") if pd.notna(v) else "not applicable" for v in d.ckb_overlap],
                        "East Asian instrument": d.instr_label.fillna("--").values, "Instrument is the lead variant": [yn(tf(v)) if pd.notna(v) else "--" for v in d.is_abs_lead],
                        "Instrument is a proxy in the outcome file": [yn(tf(v)) if pd.notna(v) else "--" for v in d.instr_is_proxy],
                        "Variants tied at the lead p": d.n_tied_top_p.fillna(-1).astype(int).map(lambda v: "--" if v < 0 else str(v)).values,
                        "Instrument F": d.F.values, "MR estimate (per CKB protein unit; binary outcomes: log-odds; eGFR: rank-inverse-normal eGFR with its sign reversed)": d.mr_beta.values, "MR SE": d.mr_se.values,
                        "MR p": d.mr_p.values, "PP.H0": d["PP.H0"].values, "PP.H1": d["PP.H1"].values, "PP.H2": d["PP.H2"].values, "PP.H3": d["PP.H3"].values,
                        "PP.H4": d["PP.H4"].values, "Colocalization call, default prior": d.call.map(CALL).values,
                        "PP.H3, conservative prior": d["PP.H3_p12_1e-6"].values, "PP.H4, conservative prior": d["PP.H4_p12_1e-6"].values,
                        "Robust to the conservative prior": [yn(tf(a) and tf(b)) if pd.notna(a) else "not testable" for a, b in zip(d.support, d["support_p12_1e-6"])]})
    assert d.call.isin(list(CALL)).all(), set(d.call) - set(CALL)
    return out


s7b = pd.concat([s7b_rows(cen, "ancestry-matched"), s7b_rows(ckf, "CKB-free sensitivity (Biobank Japan)")], ignore_index=True)
for c in ["Instrument F", "MR estimate (per CKB protein unit; binary outcomes: log-odds; eGFR: rank-inverse-normal eGFR with its sign reversed)", "MR SE", "MR p", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
          "PP.H3, conservative prior", "PP.H4, conservative prior"]:
    s7b[c] = s7b[c].map(lambda v: v if pd.notna(v) else "not testable")
sheet("S7b", "Table S7b. Ancestry-matched East Asian analysis: CKB instrument, East Asian outcome GWAS.",
      ("The %d eligible pairs over %d proteins (CKB Olink cis-pQTL at p < %s within ±%g kb, minor-allele frequency ≥ %g), "
       "followed by the CKB-free sensitivity analysis, which replaces outcome resources that include CKB with Biobank Japan. Colocalization used the "
       "same priors and rule as the European analysis; robust means the pair still met the full rule at p12 = %s. MR estimates use the CKB protein effect, "
       "which is per SD of rank-inverse-normal protein level (Table S2), so cross-ancestry magnitudes are approximate; for eGFR the outcome is "
       "rank-inverse-normal eGFR, oriented as eGFR* (positive = lower eGFR).") % (LV("L216")[0], LV("L216")[1], sci(LV("L231")), LV("L232"), LV("L234"), sci(PRIOR6)),
      s7b, {"Analysis": "ancestry-matched / CKB-free sensitivity", "Protein": "gene symbol", "Outcome": "outcome",
            "Evidence tier (European)": "tier of the European pair", "East Asian outcome resource": "outcome GWAS used, named as in Table S2",
            "Outcome resource includes CKB": "yes / no", "East Asian instrument": "CKB instrument (rsID)", "Instrument is the lead variant": "yes / no",
            "Instrument is a proxy in the outcome file": "yes / no", "Variants tied at the lead p": "variants sharing the file's minimum p",
            "Instrument F": "F-statistic in CKB", "MR estimate (per CKB protein unit; binary outcomes: log-odds; eGFR: rank-inverse-normal eGFR with its sign reversed)": "Wald ratio per SD of rank-inverse-normal CKB protein level; log-odds for binary outcomes; for eGFR, rank-inverse-normal eGFR with the sign reversed, so positive = lower eGFR", "MR SE": "standard error",
            "MR p": "p", "PP.H0": "coloc.abf posterior, no association", "PP.H1": "protein only", "PP.H2": "outcome only", "PP.H3": "distinct variants",
            "PP.H4": "shared variant", "Colocalization call, default prior": "colocalized or the most probable alternative",
            "PP.H3, conservative prior": "at p12 = %s" % E6, "PP.H4, conservative prior": "at p12 = %s" % E6, "Robust to the conservative prior": "full rule at p12 = %s" % E6},
      ["cen", "ckf", "tiered"], {c: "0.000" for c in ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H3, conservative prior", "PP.H4, conservative prior"]} |
      {"Instrument F": "0.0", "MR estimate (per CKB protein unit; binary outcomes: log-odds; eGFR: rank-inverse-normal eGFR with its sign reversed)": "0.0000", "MR SE": "0.0000", "MR p": "0.00E+00"})

# ================================================================================================ S8 Supplementary Note tables
ta = pd.read_csv(SRC["topa"], sep="\t").drop(columns=["section"])
ZTRUNC = re.search(r"variants with \|Z\| < (\d+) in both outcomes",
                   io.open(os.path.join(S1, "S1R1_02_topology_note/S1R1_topo_methods_note.md"), encoding="utf-8").read()).group(1)
tb = pd.read_csv(SRC["topb"], sep="\t").drop(columns=["section"])
tc = pd.read_csv(SRC["topc"], sep="\t").drop(columns=["section"])
ff = ta[(ta.source_1 == "FinnGen R10") & (ta.source_2 == "FinnGen R10") & ta.outcome_1.isin(["CAD", "HF", "AF", "T2D"]) & ta.outcome_2.isin(["CAD", "HF", "AF", "T2D"]) & (ta.outcome_1 != ta.outcome_2)]
eq([float(ff.null_z_correlation.min()), float(ff.null_z_correlation.max()), len(ff)], LV("L206"), "FinnGen null-z correlations (L206)")
eq(int(tb.n_random_snps.dropna().iloc[0]), LV("L100")[2], "random variants (L100)")
eq(int(tc.observed_count.iloc[0]), LV("L102")[0], "multi-outcome count (L102)")
ZT = "|z| < %s" % ZTRUNC
# column -> (header, README description with units and sign conventions)
A = {"outcome_1": ("Outcome 1", "first outcome GWAS of the pair"), "outcome_2": ("Outcome 2", "second outcome GWAS (the same as Outcome 1 on diagonal rows)"),
     "source_1": ("Source 1", "resource of Outcome 1"), "source_2": ("Source 2", "resource of Outcome 2"),
     "n_snps": ("Variants used", "variants entering LD-score regression (count)"),
     "mean_chi_square": ("Mean chi-square (diagonal rows)", "mean association chi-square of the GWAS (unitless); diagonal rows only"),
     "ldsc_intercept": ("LD-score-regression intercept", "cross-trait intercept (off-diagonal rows) or the GWAS's own intercept (diagonal rows); unitless"),
     "ldsc_intercept_se": ("Intercept SE", "block-jackknife standard error of the intercept"),
     "null_z_correlation": ("Null-z correlation", "cross-trait intercept divided by the square root of the two own intercepts: the correlation of null z-scores from shared participants and stratification (unitless, -1 to 1; 1 on diagonal rows)"),
     "null_z_correlation_se": ("Null-z correlation SE", "standard error of the null-z correlation"),
     "chr2_6_truncated": ("Chromosomes 2 and 6, %s: correlation" % ZT, "correlation of the two outcomes' z-scores over variants on chromosomes 2 and 6 (MHC excluded) with %s in both; FinnGen pairs only (unitless)" % ZT),
     "chr2_6_truncated_se": ("Chromosomes 2 and 6: SE", "standard error of that correlation"),
     "chr2_6_truncation_corrected": ("Chromosomes 2 and 6: truncation-corrected", "the same correlation corrected for the truncation, assuming unit-variance null z-scores (unitless)"),
     "chr2_6_corrected_ci_low": ("Truncation-corrected, 95% interval lower", "lower limit of the 95% interval of the corrected correlation"),
     "chr2_6_corrected_ci_high": ("Truncation-corrected, 95% interval upper", "upper limit of the 95% interval of the corrected correlation"),
     "chr2_6_untruncated": ("Chromosomes 2 and 6: all variants", "correlation over all variants on chromosomes 2 and 6 (MHC excluded), without truncation (unitless)")}
B = {"outcome_1": ("Outcome 1", "first outcome of the pair"), "outcome_2": ("Outcome 2", "second outcome of the pair"),
     "both_finngen": ("Both outcomes from FinnGen", "yes / no"), "n_proteins": ("Proteins", "instrumented proteins with both outcomes testable (count)"),
     "table3_manuscript": ("Correlation across instrumented proteins (as first reported)", "Spearman correlation of MR z-scores across proteins as printed in the earlier version of the analysis (eGFR as eGFR*)"),
     "table3_recomputed": ("Correlation across instrumented proteins (recomputed)", "the same Spearman correlation recomputed here (unitless, -1 to 1)"),
     "table3_ci_low": ("95% interval lower", "lower limit of the bootstrap 95%% interval over %d-Mb genomic regions" % LV("L444")[2]),
     "table3_ci_high": ("95% interval upper", "upper limit of the bootstrap 95%% interval over %d-Mb genomic regions" % LV("L444")[2]),
     "overlap_only_expected": ("Expected from shared participants alone", "mean Spearman correlation over simulations in which true effects are independent across outcomes and only the noise is correlated (LD-score-regression correlations)"),
     "overlap_only_simulation_sd": ("Simulation SD", "standard deviation of the expected correlation over simulations"),
     "overlap_only_share_of_table3": ("Expected share of the protein correlation", "expected correlation divided by the recomputed protein correlation (ratio)"),
     "overlap_only_expected_full_ldsc_covariance": ("Expected, full intercept covariance", "as Expected from shared participants alone, with the full LD-score-regression intercept covariance as the noise"),
     "overlap_only_share_full_ldsc_covariance": ("Expected share, full intercept covariance", "that expectation divided by the recomputed protein correlation (ratio)"),
     "random_snp_spearman": ("Correlation at random variants", "Spearman correlation of the two outcomes' z-scores over random genome-wide variants, with the sign of eGFR reversed"),
     "random_snp_ci_low": ("Random variants, 95% interval lower", "lower limit of the 95% interval from resampling contiguous genomic blocks"),
     "random_snp_ci_high": ("Random variants, 95% interval upper", "upper limit of that interval"),
     "n_random_snps": ("Random variants", "random variants used (count)"),
     "random_snp_share_of_table3": ("Random-variant correlation as a share of the protein correlation", "correlation at random variants divided by the recomputed protein correlation (ratio)")}
C = {"analysis_set": ("Analysis set", "the units counted: proteins, proteins whose lead instrument has an LD score, or random variant sets"),
     "unit": ("Unit", "what one count is: a protein, a locus (proteins chained where adjacent gene starts lie within 1 Mb) or a variant"),
     "category": ("Category", "all units, those whose hits are all in FinnGen outcomes, or those with at least one stroke or kidney hit"),
     "null_model": ("Null model", "within-unit permutation (each outcome permuted independently), Gaussian-copula permutation or signal plus correlated noise model (FinnGen outcomes share noise)"),
     "finngen_noise_correlation": ("Noise correlation between FinnGen outcomes", "the noise correlation the null model uses"),
     "renal_unit": ("Kidney outcome group", "how eGFR and CKD form one outcome group"),
     "threshold": ("p threshold", "p below which an outcome counts as a hit"),
     "n_units": ("Units", "units tested (count)"), "n_sets": ("Random sets", "random variant sets (count); random-set rows only"),
     "hit_rate": ("Hit rate", "share of unit-outcome results below the threshold (proportion)"),
     "observed_count": ("Observed count", "units with hits for at least two outcome groups; for random sets, the mean over the sets"),
     "null_mean": ("Null mean", "mean count under the null model; for random sets, the mean over the sets"),
     "null_sd": ("Null SD", "standard deviation of the null count; for random sets, the mean over the sets"),
     "null_q025": ("Null 2.5th percentile", "2.5th percentile of the null count"), "null_q975": ("Null 97.5th percentile", "97.5th percentile of the null count"),
     "z": ("z", "(observed count - null mean) / null SD; for random sets, the mean over the sets"),
     "fold_excess": ("Fold excess", "observed count / null mean; for random sets, the mean of the per-set folds (not the mean observed count over the mean null)"),
     "fold_min": ("Fold, minimum over sets", "smallest per-set fold; random-set rows only"), "fold_max": ("Fold, maximum over sets", "largest per-set fold; random-set rows only"),
     "empirical_p": ("Empirical p", "share of permutations or simulation draws with a count at least as large as observed"),
     "n_replicates": ("Permutations or simulations", "number of permutations or simulation draws"),
     "note": ("Note", "how to read the row")}
assert set(ta.columns) == set(A) and set(tb.columns) == set(B) and set(tc.columns) == set(C)
# every text value of S8c is mapped to reader terms (no source column names, no internal shorthand): an unmapped value fails the build
CTXT = {
    "analysis_set": {"proteins (cis-pQTL MR)": "proteins (cis-pQTL MR)",
                     "random SNP sets matched on LD score and MAF": "random variant sets matched on LD score and minor-allele frequency",
                     "proteins with an LD-scored lead instrument": "proteins whose lead instrument has an LD score",
                     "random SNP sets, unmatched (every 2,000th SNP from a random offset)": "random variant sets, unmatched (every 2,000th variant from a random offset)",
                     "proteins (cis-pQTL MR) (renal = eGFR only)": "proteins (cis-pQTL MR), with eGFR as the only kidney outcome"},
    "unit": {"protein": "protein", "SNP": "variant",
             "locus: one representative protein (re-selected per permutation)": "locus: one representative protein (re-selected per permutation)",
             "locus: minimum p per endpoint-unit across member proteins": "locus: minimum p per outcome group across member proteins"},
    "category": {"all": "all", "all hits in FinnGen outcomes": "all hits in FinnGen outcomes",
                 "at least one hit in Stroke or Renal": "at least one hit in stroke or a kidney outcome"},
    "null_model": {"copula permutation": "Gaussian-copula permutation", "within-unit permutation": "within-unit permutation",
                   "signal + correlated noise model": "signal plus correlated noise model",
                   "within-unit permutation (re-implementation of the paper's locus sweep)": "within-unit permutation (re-implementation of the original locus analysis)"},
    "finngen_noise_correlation": {"none (independent columns)": "none (outcomes independent)",
                                  "chr2/chr6 untruncated correlations": "chromosomes 2 and 6, all variants",
                                  "LDSC null-Z correlations": "LD-score-regression null-z correlations",
                                  "LDSC intercepts, diagonal only (no correlation)": "LD-score-regression own intercepts only (no correlation)",
                                  "LDSC intercept covariance (full)": "LD-score-regression intercept covariance (full)",
                                  "chr2/chr6 |Z|<%s correlations" % ZTRUNC: "chromosomes 2 and 6, %s" % ZT,
                                  "chr2/chr6 correlations, truncation-corrected": "chromosomes 2 and 6, truncation-corrected",
                                  "LDSC intercept covariances, unit variances": "LD-score-regression intercept covariances, unit variances"},
    "renal_unit": {"min p of eGFR, CKD": "minimum p of eGFR and CKD", "eGFR only": "eGFR only"},
    "note": {"the paper's pre-specified primary locus-level control": "the pre-specified primary locus-level control",
             "observed_count, null_mean, null_sd, z and fold_excess are means over the sets (fold_excess is the mean of per-set folds, not observed_count / null_mean); fold_min and fold_max give the range over sets":
                 "observed count, null mean, null SD and z are means over the sets; fold excess is the mean of the per-set folds, not the mean observed count over the mean null, and the smallest and largest folds give the range over the sets",
             "each endpoint-unit takes the minimum p across the locus members, so a locus counts when different members hit different endpoints; the paper's sweep script pre-specifies this unit as expected to saturate, not as evidence against sharing":
                 "each outcome group takes the minimum p across the locus members, so a locus counts when different members reach significance for different outcomes; the original locus analysis specified this unit in advance as expected to saturate, not as evidence against sharing",
             "independence reference for this category: with independent columns the copula permutation is the within-unit permutation":
                 "independence reference for this category: with independent outcomes the copula permutation equals the within-unit permutation",
             "reference row: the model with uncorrelated noise": "reference row: the model with uncorrelated noise"}}
for c, mp in CTXT.items():
    tc[c] = mapv(tc[c], mp, "S8c " + c)
# the random-set rows are the ones whose fold excess is a mean of per-set folds; the legend states both definitions, so pin which rows they are
_rs = tc.analysis_set.str.startswith("random variant sets")
assert (tc.loc[~_rs, "fold_excess"] - tc.loc[~_rs, "observed_count"] / tc.loc[~_rs, "null_mean"]).abs().max() < 1e-9, "fold excess != observed / null mean outside the random sets"
assert tc.loc[_rs, "note"].str.startswith("observed count, null mean").all() and tc.loc[~_rs, "fold_min"].isna().all()


def fill(df, what="not applicable"):
    return df.apply(lambda col: col.map(lambda v: v if pd.notna(v) else what))


s8a, s8b, s8c = (fill(ta.rename(columns={k: v[0] for k, v in A.items()})), fill(tb.rename(columns={k: v[0] for k, v in B.items()})),
                 fill(tc.rename(columns={k: v[0] for k, v in C.items()})))
s8b["Both outcomes from FinnGen"] = s8b["Both outcomes from FinnGen"].map(lambda v: yn(tf(v)))
sheet("S8a", "Table S8a. Shared sampling noise between the outcome GWAS.",
      ("Genome-wide cross-trait LD-score-regression intercepts and the implied null-z correlation (sampling noise shared through overlapping "
       "participants and stratification) for each pair of outcome GWAS; diagonal rows give each GWAS's own intercept. For pairs of FinnGen outcomes, "
       "a check on chromosomes 2 and 6 (MHC excluded) used variants with %s in both outcomes, with and without correction for truncation. "
       "Methods are in the Supplementary Note (Additional file 2).") % ZT,
      s8a, {v[0]: v[1] for v in A.values()}, ["topa"])
sheet("S8b", "Table S8b. Correlation between outcomes across instrumented proteins, from shared participants alone, and at random variants.",
      ("Per outcome pair: the Spearman correlation of MR z-scores across instrumented proteins (eGFR as eGFR*, so that positive means higher "
       "risk), with bootstrap 95%% intervals over %d-Mb genomic regions; the value expected from shared sampling noise alone, by simulation with the LD-score-"
       "regression correlations; and the same correlation at random genome-wide variants. Methods are in the Supplementary Note (Additional file 2).") % LV("L444")[2],
      s8b, {v[0]: v[1] for v in B.values()}, ["topb"])
sheet("S8c", "Table S8c. Proteins associated with more than one outcome: observed counts against null models.",
      ("Counts of units (proteins, loci or random variants) below the p threshold for at least two outcome groups (eGFR and CKD counted as one), "
       "against within-unit permutation, null models that keep the shared noise between FinnGen outcomes, and random variant sets. Fold excess is "
       "the observed count over the null mean. For the random variant sets, the observed count, null mean, null SD and z are means over the sets, "
       "and fold excess is the mean of the per-set folds, with the smallest and largest fold beside it. Methods are in the Supplementary Note "
       "(Additional file 2)."),
      s8c, {v[0]: v[1] for v in C.values()}, ["topc"])

# ================================================================================================ README + vocabulary + write
MS_TITLE = re.search(r"\*\*Title[^*]*\*\*\s*(.+)", MB.mread("title_page")).group(1).strip()   # the title exactly as the manuscript build reads it
BAN = [re.compile(p, re.I) for p in PG.BANNED]


def wrapped_lines(text, cap):
    """Lines a bold header needs when wrapped by words into `cap` characters (Excel wraps at spaces and hyphens)."""
    n, cur = 1, 0
    for wd in re.split(r"(?<=[ -])", text):
        if cur and cur + len(wd) > cap:
            n, cur = n + 1, 0
        cur += len(wd)
    return n


def vocab(text, where):
    bad = [m.group(0) for p in BAN for m in p.finditer(text)]
    for m in re.finditer(r"replicat\w*", text, flags=re.I):
        win = text[max(0, m.start() - 40):m.end() + 5]
        if not any(re.search(p, win, flags=re.I) for p in PG.REPLICAT_OK):
            bad.append(m.group(0))
    assert not bad, "vocabulary (%s): %s in %r" % (where, bad, text[:120])


# where the manuscript cites each sheet, derived from the source map's components (not typed): in the running text, only in a table note or
# figure legend, or nowhere. A citation is "Table S4", "Tables S8a-S8c" (a range) or a list after "Tables" ("Tables S3 (pairs), S4 (..) and S7a").
SHEET_ORDER = [s[0] for s in SHEETS]
_ITEM = r"S\d[a-c]?(?:\s*[–-]\s*S\d[a-c]?)?(?:\s*\([^)]*\))?"
_CITE = re.compile(r"Tables?\s+(%s(?:(?:,\s*|\s+and\s+)%s)*)" % (_ITEM, _ITEM))


def cited_sheets(text):
    out = set()
    for m in _CITE.finditer(PG.TAG.sub("", text)):
        for a_, b_ in re.findall(r"(S\d[a-c]?)(?:\s*[–-]\s*(S\d[a-c]?))?", re.sub(r"\([^)]*\)", "", m.group(1))):
            assert a_ in SHEET_ORDER and (not b_ or b_ in SHEET_ORDER), ("the manuscript cites a sheet the workbook does not have", a_, b_)
            lo, hi = SHEET_ORDER.index(a_), SHEET_ORDER.index(b_ or a_)
            out.update(SHEET_ORDER[lo:hi + 1])
    return out


_run_text = set().union(*[cited_sheets(MB.mread(c)) for c in ("abstract", "research_insights", "background", "methods", "results", "discussion",
                                                                "conclusions", "declarations")])
_note_text = set().union(*[cited_sheets(MB.mread(c)) for c in MB.MAP if re.fullmatch(r"table_\d+|legend_fig_\d+", c)])
CITED = {s: ("running text" if s in _run_text else ("a main-text table or figure legend only, not the running text" if s in _note_text else "not cited"))
         for s in SHEET_ORDER}
readme = []
for name, title, legend, df, cols, srcs, _f in SHEETS:
    readme.append({"Sheet": name, "Title": title, "Contents": legend,
                   "Columns (units and sign conventions)": "\n".join("%s: %s" % (c, cols[c]) for c in df.columns),
                   "Source file(s)": "\n".join(repo_src(SRC[s]) for s in srcs) or "--", "Source md5": "\n".join(md5(SRC[s]) for s in srcs) or "--",
                   "Cited in the manuscript": CITED[name]})
readme = pd.DataFrame(readme)
nchk = 0
for name, title, legend, df, cols, srcs, _f in SHEETS:
    for t in [title, legend] + list(df.columns) + ["%s %s" % (c, cols[c]) for c in df.columns]:
        vocab(str(t), name); nchk += 1
    for c in df.columns:
        for vv in df[c]:
            if isinstance(vv, str):
                vocab(vv, "%s / %s" % (name, c)); nchk += 1
for c in ["Sheet", "Title", "Contents", "Columns (units and sign conventions)", "Cited in the manuscript"]:   # source paths are identifiers and are exempt
    for vv in readme[c]:
        vocab(vv, "README / " + c); nchk += 1
# no printed number may be the double-precision floor that a source uses for an underflowed p (it is not the value)
for name, title, legend, df, cols, srcs, _f in SHEETS:
    for c in df.columns:
        tiny = [v for v in df[c] if isinstance(v, (float, np.floating)) and 0 < abs(v) < 1e-300]
        assert not tiny, ("%s / %s: a value at the underflow floor is printed as a number" % (name, c), tiny[:3])
# open items: a PLACEHOLDER sheet or any visible [... TO CONFIRM] marker (README included) blocks packaging mode
OPEN = re.compile(r"PLACEHOLDER|\[[A-Z ]+ TO CONFIRM\]")
placeholders = {}
for n, t, l, d, cols, *_ in SHEETS:
    k = sum(len(OPEN.findall(str(x))) for x in [t, l] + list(d.columns) + list(cols.values()) + list(d.values.ravel()))
    if k:
        placeholders[n] = k
k = sum(len(OPEN.findall(str(x))) for x in readme.values.ravel())
if k:
    placeholders["README"] = k
# stale Stage 2c skeleton outputs (written by the retired _retired/P31R_build_addfiles.py) must not ship beside the real workbook
STALE = [f for f in ("P31R_Additional_file_1_DRAFT.xlsx", "P31R_Table_S12.csv") if os.path.exists(os.path.join(OUT, f))]
if os.environ.get("P31R_AF1_FINAL") == "1" and (placeholders or STALE or REL_STALE):
    sys.exit("REFUSING TO WRITE (final mode): open markers per sheet %s; stale skeleton outputs in %s: %s; README sources whose release record predates them: %s" % (placeholders, OUT, STALE, REL_STALE))

out = os.path.join(OUT, "P31R_Additional_file_1.xlsx")
tmp = out + ".tmp.xlsx"
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
with pd.ExcelWriter(tmp, engine="openpyxl") as xw:
    readme.to_excel(xw, sheet_name="README", index=False, startrow=4)
    for name, title, legend, df, cols, srcs, fmt in SHEETS:
        df.to_excel(xw, sheet_name=name, index=False, startrow=3)
        ws = xw.sheets[name]
        ws["A1"] = title; ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = legend; ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
        ncol = len(df.columns)
        ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=max(ncol, 8))
        ws.row_dimensions[2].height = min(409, 15 * (1 + len(legend) // 140))
        hdr_lines = 1
        for j, c in enumerate(df.columns, 1):
            ws.cell(row=4, column=j).font = Font(bold=True)
            ws.cell(row=4, column=j).alignment = Alignment(wrap_text=True, vertical="top")
            w = max(10, min(60, 2 + max([len(str(c)) // 2] + [len(str(vv)) for vv in df[c].head(200)])))
            ws.column_dimensions[get_column_letter(j)].width = w
            hdr_lines = max(hdr_lines if j > 1 else 1, wrapped_lines(str(c), int(w) - 3) + 1)
            if c in fmt:
                for i in range(len(df)):
                    cell = ws.cell(row=5 + i, column=j)
                    if isinstance(cell.value, (int, float)):
                        cell.number_format = fmt[c]
        ws.row_dimensions[4].height = min(120, 15 * hdr_lines)
        ws.freeze_panes = "B5"
    ws = xw.sheets["README"]
    ws["A1"] = "Additional file 1"
    ws["A2"] = "Supplementary tables S1 to S8 for: %s" % MS_TITLE
    ws["A3"] = "Each sheet starts with its title (row 1) and legend (row 2); the table header is on row 4. This README describes every sheet's columns and sources."
    ws["A1"].font = Font(bold=True, size=12)
    for a_ in ("A1", "A2", "A3"):
        ws[a_].alignment = Alignment(horizontal="left", vertical="top", wrap_text=False)
    for j, wdt in enumerate([8, 40, 90, 90, 60, 36, 34], 1):
        ws.column_dimensions[get_column_letter(j)].width = wdt
    for row in ws.iter_rows(min_row=6):
        for cell in row:
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    for j in range(1, 8):
        ws.cell(row=5, column=j).font = Font(bold=True)
import openpyxl
wb = openpyxl.load_workbook(tmp)
now = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0, tzinfo=None)   # openpyxl writes UTC and stamps `modified` itself on save
wb.properties.creator = AUTHOR; wb.properties.lastModifiedBy = AUTHOR; wb.properties.title = "Additional file 1"
wb.properties.description = None; wb.properties.subject = None; wb.properties.created = now; wb.properties.modified = now
wb.save(tmp)
# the package's app.xml names the generator ("Microsoft Excel Compatible / Openpyxl"): drop it, as the manuscript build does for its .docx
import zipfile
with zipfile.ZipFile(tmp) as zin, zipfile.ZipFile(tmp + ".zip", "w", zipfile.ZIP_DEFLATED) as zout:
    for it in zin.infolist():
        data = zin.read(it.filename)
        if it.filename == "docProps/app.xml":
            data = re.sub(r"<AppVersion>.*?</AppVersion>", "", re.sub(r"<Application>.*?</Application>", "", data.decode("utf-8"))).encode("utf-8")
        zout.writestr(it, data)
os.replace(tmp + ".zip", tmp)
os.replace(tmp, out)

# ---- read the written file back (a different reader than the writer) and assert every sheet's body survived intact
rb = openpyxl.load_workbook(out, read_only=True)
assert rb.sheetnames == ["README"] + [s[0] for s in SHEETS], rb.sheetnames
# metadata, read back from the written package: the author is the literal default of the global authorship rule (not the AUTHOR variable,
# which a mutation would also change), the generator names are gone, and the dates are in order
_pp = openpyxl.load_workbook(out).properties
assert _pp.creator == "Bertrand Chin-Ming Tan" and _pp.lastModifiedBy == "Bertrand Chin-Ming Tan", ("workbook author", _pp.creator, _pp.lastModifiedBy)
assert not (_pp.description or _pp.subject or _pp.keywords or _pp.category), "workbook properties carry a description, subject or keywords"
assert _pp.created <= _pp.modified, ("workbook created after modified", _pp.created, _pp.modified)
with zipfile.ZipFile(out) as _z:
    _meta = _z.read("docProps/core.xml").decode("utf-8") + _z.read("docProps/app.xml").decode("utf-8")
assert not re.search(r"openpyxl|Microsoft Excel", _meta, re.I), ("generator default in the workbook properties", _meta[:300])
for name, title, legend, df, cols, srcs, fmt in SHEETS:
    back = pd.read_excel(out, sheet_name=name, header=3, keep_default_na=False, na_values=[])
    assert list(back.columns) == list(df.columns) and len(back) == len(df), (name, "read-back shape")
    assert pd.read_excel(out, sheet_name=name, header=None, nrows=1).iloc[0, 0] == title, (name, "title")
    for c in df.columns:
        a, b = df[c].tolist(), back[c].tolist()
        for va, vb in zip(a, b):
            if isinstance(va, (int, float, np.integer, np.floating)) and not isinstance(va, bool):
                assert isinstance(vb, (int, float, np.integer, np.floating)) and (va == vb or abs(va - vb) <= 1e-12 * max(1, abs(va))), (name, c, va, vb)
            else:
                assert str(va) == str(vb), (name, c, va, vb)
print("wrote %s (md5 %s)" % (out, md5(out)))
for name, title, legend, df, cols, srcs, fmt in SHEETS:
    print("  %-4s %3d rows x %2d cols  %s" % (name, len(df), len(df.columns), title))
print("vocabulary checked on %d strings; read-back identical; author %s" % (nchk, AUTHOR))
print(DOCX_CHECK)
if placeholders:
    print("!!! open markers (PLACEHOLDER / [... TO CONFIRM]) per sheet: %s (P31R_AF1_FINAL=1 refuses to write)" % placeholders)
if STALE:
    print("!!! stale Stage 2c skeleton output(s) beside the workbook in %s: %s (delete at packaging; P31R_AF1_FINAL=1 refuses)" % (OUT, STALE))
if REL_STALE:
    print("!!! README source files changed since the release manifest was built (the repository path is of an older version; rebuild the release; P31R_AF1_FINAL=1 refuses): %s" % REL_STALE)
