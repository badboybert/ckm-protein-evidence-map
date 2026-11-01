# -*- coding: utf-8 -*-
"""Release module: the screen ledger (Paper 3.1 R2, Stage 5 step 2; spec section 5, interface section 2).

Writes four files under out/screen/ (nothing under withheld/: no TPMI- or Open Targets-derived value is read here):
  screen_ledger.csv                 one row per protein x outcome in the multiple-testing family (rows = L003)
  screen_nonprimary_estimates.csv   the surplus Egger / weighted-median rows of the multi-instrument pairs (not in the family)
  screen_reason_codes.csv           codebook for the reason_code column, counts derived from the pipeline run log
  cis_pqtl_instruments.csv          the instrument list, byte-for-byte

Rules this module keeps:
  * No value is typed. Every count a check compares against is a ledger row (json.loads of its `value`) or is read from a file.
    The stored columns of the screen are carried as TEXT, so the shipped digits are the screen's own digits.
  * Derived values: `p_family`, `fdr_q_bh_family`, `fdr_sig_bh_family` (decision 6) and `reason_code`. Nothing else.
  * The canonical BH column is the ledger's definition (ext_supp.py l.53-57: untestable combinations enter at p = pad, the
    threshold is the largest sorted p with p <= rank/m*alpha). The pipeline's own form is primary_fdr_set() in 00_config.R
    (p.adjust(BH) over the padded grid, q < alpha); the q column below is that form and must give the same set. The ledger code
    is exec'd inside the builder and cannot be imported, so it is implemented once here and must reproduce L003-L009, L242,
    L400 and L401 or the build refuses. The `withheld` argument is accepted and not used.
  * Reason codes come from a file. The pipeline wrote per-cell audit files (mr_no_result_cells_*.csv, mr_no_outcome_data_cells_*.csv;
    01_mr_screen.R l.1069-1071) that are not on this machine. Its run log (01_mr_screen.log) records the class totals per outcome and,
    for every attempted cell, a harmonisation line, followed by an analysis line when an estimate was reached. A cell with a
    harmonisation line and no analysis line is attempted-no-result; an untestable cell with no harmonisation line is no-outcome-data.
    Both readings are cross-checked against the log's own totals and against the set of testable cells in the screen (S20-S23).
  * Checks raise ScreenCheckError(check_id, message). Ids are S + two digits and are matched exactly.
Run `python -I -X utf8 rel_screen.py --selftest` for the mutation ladder (each mutation must trip the check it names).
"""
import csv
import hashlib
import io
import json
import math
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

PRODUCER = "rel_screen.py"

SCREEN = "V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv"
INSTR = "V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv"
CONV = "V4/source_data/tables_from_rds/RDS_CONVERSIONS.csv"
LOG = "paper3_deposit_request_v4/scripts/reports/Paper3_Program/run_logs/REGEN_20260817/01_mr_screen.log"
CONFIG = "paper3_deposit_request_v4/scripts/R/00_config.R"
CASCADE01 = "paper3_deposit_request_v4/scripts/scripts/paper3.1_ckm/cascade/01_mr_screen.R"
INPUTS = [SCREEN, INSTR, CONV, LOG, CONFIG, CASCADE01]
LEDGER_CSV = "R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv"   # used by the CLI and selftest only; build() gets the frame

OUT_LEDGER = "screen/screen_ledger.csv"
OUT_NONPRIMARY = "screen/screen_nonprimary_estimates.csv"
OUT_CODES = "screen/screen_reason_codes.csv"
OUT_INSTR = "screen/cis_pqtl_instruments.csv"

# Declared, not defaulted: the screen's stored columns, in stored order (S01 refuses any other header).
STORED_COLUMNS = ["gene", "outcome_id", "n_snps", "method", "b", "se", "pval", "OR", "OR_lo95", "OR_up95", "n_proxy",
                  "any_proxy", "Q", "Q_df", "Q_pval", "outcome_name", "bonferroni_sig", "fdr_sig"]
ID_COLUMNS = ["gene", "outcome_id", "outcome_name"]
LEGACY = "fdr_sig_legacy_per_method"
ESTIMATE_COLUMNS = [LEGACY if c == "fdr_sig" else c for c in STORED_COLUMNS if c not in ID_COLUMNS]
LEDGER_HEADER = ["gene", "outcome_name", "outcome_id", "testable", "reason_code"] + ESTIMATE_COLUMNS + \
    ["p_family", "fdr_q_bh_family", "fdr_sig_bh_family"]
NONPRIMARY_HEADER = ["gene", "outcome_name", "outcome_id", "role"] + ESTIMATE_COLUMNS
NONPRIMARY_ROLE = "non_primary_estimator_not_in_family"
CODES_HEADER = ["reason_code", "cells", "outcomes", "definition", "evidence"]

WALD, IVW = "Wald ratio", "Inverse variance weighted"      # the two primary methods; S01 checks them against 00_config.R
CODE_NOD = "no_outcome_data"
CODE_ANR = "attempted_no_result"

# Outcome-GWAS provider by outcome_id prefix (declared; an outcome_id that matches none refuses the build).
PROVIDER_BY_PREFIX = [("FinnGen_", "FinnGen"), ("Wuttke2019_", "CKDGen"), ("MEGASTROKE_", "MEGASTROKE")]
EXPOSURE_PROVIDER = "UKB-PPP"

LEDGER_IDS = ["L001", "L002", "L003", "L004", "L005", "L006", "L007", "L008", "L009", "L239", "L242", "L400", "L401", "L928"]


class ScreenCheckError(AssertionError):
    def __init__(self, check_id, msg):
        super().__init__("%s: %s" % (check_id, msg))
        self.check_id = check_id


def _check(check_id, cond, msg):
    if not cond:
        raise ScreenCheckError(check_id, msg)


def _close(a, b, rel=1e-12):
    return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=0.0)


def _text(p):
    return Path(p).read_text(encoding="utf-8", errors="replace")


def _md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def _lv(ledger, lid):
    r = ledger.loc[ledger["id"] == lid, "value"]
    _check("S00", len(r) == 1, "ledger row %s: expected exactly one row, found %d" % (lid, len(r)))
    return json.loads(r.iloc[0])


# ------------------------------------------------------------------------------------------------ Benjamini-Hochberg
def _bh_flags(p, m, alpha):
    """Ledger definition (ext_supp.py l.53-57, S1R1_build_ledger.py bh_count): threshold = largest sorted p with
    p <= rank/m*alpha, over the family padded to m; flag = p <= threshold. Returns (threshold or None, bool array)."""
    p = np.asarray(p, float)
    _check("S11", len(p) == m, "BH vector has %d entries, family size is %d" % (len(p), m))
    ps = p[np.argsort(p, kind="stable")]
    ok = ps <= np.arange(1, m + 1) / m * alpha
    if not ok.any():
        return None, np.zeros(m, bool)
    thr = float(ps[ok.nonzero()[0].max()])
    return thr, p <= thr


def _bh_q(p):
    """R p.adjust(p, "BH"): min over larger ranks of p*m/rank, capped at 1; tied p share one value."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p, kind="stable")
    qs = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(qs, 1.0)
    return q


def _bh_count_ledger(p, m, q):
    """S1R1_build_ledger.py bh_count, verbatim logic: testable p only, ranks against the full family size m."""
    p = np.sort(np.asarray(p))
    k = np.arange(1, len(p) + 1)
    ok = p <= k / m * q
    return int(ok.nonzero()[0].max() + 1) if ok.any() else 0


# ------------------------------------------------------------------------------------------------ run log
RE_CELLS = re.compile(r"^MR cells: (\d+) with-result, (\d+) attempted-no-result, (\d+) no-outcome-data, "
                      r"(\d+) null-exposure-gene \| full-grid (\d+)/(\d+) (\S+)\s*$")
RE_OUTCOME_ROW = re.compile(r"^\s*\d+:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$")
RE_RUN1 = re.compile(r"^N_PROTEINS=(\d+), Bonf=(\S+) \| rows=(\d+) genes=(\d+) \| primary=(\d+)\s*$")
RE_RUN2 = re.compile(r"^Bonferroni=(\d+) \| FDR\(by=method\)=(\d+) \| suggestive=(\d+)\s*$")
LOG_LOSS_HEADER = "per-outcome MR cell loss (attempted / with-result / no-result):"
RE_HARM = re.compile(r"^Harmonising (\S+) \((\S+)\) \((\S+)\) and (\S+) \((\S+)\)$")
RE_ANALYSING = re.compile(r"^Analysing '(\S+)' on '(\S+)'$")


def _scan_cells(lines):
    """One entry per (gene, outcome_id) the pipeline attempted: its harmonisation lines, its analysis lines, its first line number.
    An analysis line belongs to the harmonised cell directly above it (same exposure label and outcome); any other is an orphan."""
    cells, cur, orphan = {}, None, 0
    for i, l in enumerate(lines, 1):
        m = RE_HARM.match(l)
        if m:
            key, exp = (m.group(1), m.group(4)), m.group(3)
            if cur is not None and cur["key"] == key and cur["exp"] == exp:
                cur["h"] += 1
            else:
                _check("S23", key not in cells, "run log: cell %s is harmonised in two separate places" % (key,))
                cur = dict(key=key, exp=exp, h=1, a=0, line=i)
                cells[key] = cur
            continue
        m = RE_ANALYSING.match(l)
        if m:
            if cur is not None and (cur["exp"], cur["key"][1]) == (m.group(1), m.group(2)):
                cur["a"] += 1
            else:
                orphan += 1
    return cells, orphan


def _read_log(v3):
    lines = _text(Path(v3) / LOG).splitlines()
    out = {}
    out["cells"], out["orphan_analysing"] = _scan_cells(lines)
    hit = [(i, RE_CELLS.match(l)) for i, l in enumerate(lines) if RE_CELLS.match(l)]
    _check("S20", len(hit) == 1, "run log: expected exactly one 'MR cells:' line, found %d" % len(hit))
    i, m = hit[0]
    out["cells_line"] = i + 1
    (out["with_result"], out["attempted_no_result"], out["no_outcome_data"], out["null_exposure"],
     out["grid_num"], out["grid_den"]) = [int(m.group(k)) for k in (1, 2, 3, 4, 5, 6)]
    out["grid_status"] = m.group(7)
    hdr = [i for i, l in enumerate(lines) if l.strip() == LOG_LOSS_HEADER]
    _check("S21", len(hdr) == 1, "run log: expected exactly one per-outcome cell-loss header, found %d" % len(hdr))
    rows, first, skipped, j = [], None, 0, hdr[0] + 1
    while j < len(lines):
        mm = RE_OUTCOME_ROW.match(lines[j])
        if mm:
            first = j + 1 if first is None else first
            rows.append((mm.group(1), int(mm.group(2)), int(mm.group(3)), int(mm.group(4))))
        elif rows:
            break
        else:
            skipped += 1
            _check("S21", skipped <= 3, "run log: per-outcome table not found after its header")
        j += 1
    out["outcome_rows"] = rows
    out["outcome_lines"] = (first, first + len(rows) - 1) if rows else (None, None)
    r1 = [(i, RE_RUN1.match(l)) for i, l in enumerate(lines) if RE_RUN1.match(l)]
    r2 = [(i, RE_RUN2.match(l)) for i, l in enumerate(lines) if RE_RUN2.match(l)]
    _check("S22", len(r1) == 1 and len(r2) == 1, "run log: expected one run-summary pair, found %d and %d" % (len(r1), len(r2)))
    g = r1[0][1]
    out["n_proteins"], out["bonf_text"] = int(g.group(1)), g.group(2)
    out["rows"], out["genes"], out["primary"] = int(g.group(3)), int(g.group(4)), int(g.group(5))
    h = r2[0][1]
    out["bonferroni"], out["fdr_by_method"], out["suggestive"] = int(h.group(1)), int(h.group(2)), int(h.group(3))
    out["summary_lines"] = (r1[0][0] + 1, r2[0][0] + 1)
    return out


def _script_lines(v3):
    """Line numbers (1-based) of the two cell-class statements in 01_mr_screen.R, found by text, plus the audit-file anchors."""
    t = _text(Path(v3) / CASCADE01)
    lines = t.splitlines()
    nod = [i + 1 for i, l in enumerate(lines) if "no_data_cells[[length(no_data_cells)+1]]" in l]
    att = [i + 1 for i, l in enumerate(lines) if "attempted_cells[[length(attempted_cells)+1]]" in l]
    _check("S01", len(nod) == 1 and len(att) == 1, "01_mr_screen.R: cell-class statements not found exactly once")
    for needle in ('mr_artifact("no_result_cells", "csv")', 'mr_artifact("no_outcome_data_cells", "csv")',
                   'p.adjust(pval, method="BH") < FDR_ALPHA, by=method'):
        _check("S01", t.count(needle) >= 1, "01_mr_screen.R: expected text not found: %s" % needle)
    return nod[0], att[0]


# ------------------------------------------------------------------------------------------------ config
def _read_config(v3):
    t = _text(Path(v3) / CONFIG)
    m = re.search(r'^FDR_ALPHA_PRIMARY_METHODS\s*<-\s*c\(([^)]*)\)', t, re.M)
    a = re.search(r'^FDR_ALPHA\s*<-\s*([0-9.eE+-]+)', t, re.M)
    _check("S01", m is not None and a is not None, "00_config.R: FDR_ALPHA / FDR_ALPHA_PRIMARY_METHODS not found")
    _check("S01", "untested cells padded at" in t and "p = 1" in t, "00_config.R: the p = 1 padding definition is not stated")
    return tuple(re.findall(r'"([^"]+)"', m.group(1))), float(a.group(1))


# ------------------------------------------------------------------------------------------------ small writers
def _write_csv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _read_csv_rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        r = csv.reader(f)
        header = next(r)
        return header, list(r)


def _tf(x):
    return "TRUE" if x else "FALSE"


def _scan_tpmi(paths):
    """S31: no shipped cell mentions TPMI (spec section 3). Returns nothing; raises."""
    for p in paths:
        t = Path(p).read_bytes().decode("utf-8", errors="replace")
        _check("S31", "tpmi" not in t.lower(), "%s mentions TPMI" % Path(p).name)


# ------------------------------------------------------------------------------------------------ the build
def _reason_code(attempted_without_result):
    return CODE_ANR if attempted_without_result else CODE_NOD


def _build(v3, out, ledger):
    v3, out = Path(v3), Path(out)
    LV = {i: _lv(ledger, i) for i in LEDGER_IDS}
    primary_methods, alpha_cfg = _read_config(v3)
    alpha = float(LV["L239"])
    _check("S01", _close(alpha_cfg, alpha), "00_config.R FDR_ALPHA %r != ledger L239 %r" % (alpha_cfg, alpha))
    pad = LV["L400"]["pad_p"]

    # ---- inputs
    txt = pd.read_csv(v3 / SCREEN, dtype=str, keep_default_na=False)     # the screen as text (what ships)
    num = pd.read_csv(v3 / SCREEN)                                       # parsed as S1R1_build_ledger.py l.55 parses it
    ins = pd.read_csv(v3 / INSTR, dtype=str, keep_default_na=False)
    _check("S01", list(txt.columns) == STORED_COLUMNS, "screen header differs from the declared stored columns: %s" % list(txt.columns))
    _check("S01", set(primary_methods) == {WALD, IVW} <= set(txt.method), "00_config.R primary methods %s differ from the declared pair" % (primary_methods,))
    _check("S01", "gene" in ins.columns, "instrument list has no gene column")
    conv = pd.read_csv(v3 / CONV, dtype=str, keep_default_na=False).set_index("csv")
    for rel, df in ((SCREEN, txt), (INSTR, ins)):
        key = "tables_from_rds/" + Path(rel).name
        _check("S02", key in conv.index, "RDS_CONVERSIONS.csv has no record for %s" % key)
        c = conv.loc[key]
        _check("S02", c["csv_md5"] == _md5(v3 / rel), "%s: md5 differs from the md5 recorded at RDS conversion" % key)
        _check("S02", int(c["n_rows"]) == len(df) and int(c["n_cols"]) == df.shape[1], "%s: shape differs from the conversion record" % key)

    prim_mask = txt.method.isin(primary_methods).values
    prim = txt[prim_mask]
    pairs = list(zip(prim.gene, prim.outcome_name))
    _check("S03", len(set(pairs)) == len(pairs), "more than one primary-method row for a (gene, outcome) pair")
    outcomes = sorted(txt.outcome_name.unique())
    id_of = dict(zip(txt.outcome_name, txt.outcome_id))
    _check("S04", len({(a, b) for a, b in zip(txt.outcome_name, txt.outcome_id)}) == len(outcomes) == len(set(id_of.values())),
           "outcome_name <-> outcome_id is not one-to-one")
    _check("S04", len(outcomes) == int(LV["L002"]), "outcomes in the screen %d != L002 %s" % (len(outcomes), LV["L002"]))
    _check("S04", all(sum(oid.startswith(pre) for pre, _ in PROVIDER_BY_PREFIX) == 1 for oid in id_of.values()),
           "an outcome_id matches none or several of the declared provider prefixes")
    lic = "third_party_derived:" + "+".join([EXPOSURE_PROVIDER] + sorted({p for oid in id_of.values() for pre, p in PROVIDER_BY_PREFIX if oid.startswith(pre)}))
    inst_genes = sorted(set(ins.gene))
    _check("S05", len(inst_genes) == int(LV["L001"]), "instrument-list genes %d != L001 %s" % (len(inst_genes), LV["L001"]))
    _check("S05", set(txt.gene) <= set(inst_genes), "screen genes outside the instrument list: %s" % sorted(set(txt.gene) - set(inst_genes))[:5])
    fam = len(inst_genes) * len(outcomes)
    _check("S06", fam == int(LV["L003"]) == int(LV["L001"]) * int(LV["L002"]), "grid %d != L003 %s" % (fam, LV["L003"]))
    n_test = int(prim_mask.sum())
    n_wald = int((prim.method == WALD).sum())
    n_ivw = int((prim.method == IVW).sum())
    n_unt = fam - n_test
    _check("S07", n_test == int(LV["L004"]), "testable %d != L004 %s" % (n_test, LV["L004"]))
    _check("S07", n_wald == int(LV["L005"]) and n_ivw == int(LV["L400"]["ivw"]) and n_wald + n_ivw == n_test,
           "Wald %d / IVW %d vs L005 %s / L400 %s" % (n_wald, n_ivw, LV["L005"], LV["L400"]["ivw"]))
    _check("S07", n_unt == int(LV["L400"]["untestable"]), "untestable %d != L400 %s" % (n_unt, LV["L400"]["untestable"]))
    nsn = pd.to_numeric(prim.n_snps)
    _check("S08", (nsn[prim.method == WALD] == 1).all() and (nsn[prim.method == IVW] >= 2).all(),
           "n_snps is not 1 for every Wald ratio and >= 2 for every IVW estimate")

    pv_cr = np.array([float(x) for x in txt.pval])          # correctly rounded parse of the shipped text
    pv_pd = num.pval.values                                  # the ledger's parse
    p_prim_cr, p_prim_pd = pv_cr[prim_mask], pv_pd[prim_mask]
    nominal = float(LV["L928"])                    # the stored 0.05 that the ledger calls the nominal cut; L006 is "p < nominal" over the primary rows
    _check("S09", int((p_prim_pd < nominal).sum()) == int(LV["L006"]), "nominal count %d != L006 %s" % (int((p_prim_pd < nominal).sum()), LV["L006"]))
    bonf_thr = alpha / fam
    _check("S10", _close(bonf_thr, LV["L242"]), "alpha/L003 %r != L242 %r" % (bonf_thr, LV["L242"]))
    stored_bonf = (txt.bonferroni_sig == "TRUE").values
    _check("S10", set(txt.bonferroni_sig) <= {"TRUE", "FALSE"}, "bonferroni_sig is not TRUE/FALSE text")
    _check("S10", int((p_prim_pd < bonf_thr).sum()) == int(LV["L007"]) == int(stored_bonf[prim_mask].sum()),
           "Bonferroni count %d / stored %d != L007 %s" % (int((p_prim_pd < bonf_thr).sum()), int(stored_bonf[prim_mask].sum()), LV["L007"]))
    _check("S10", (stored_bonf[prim_mask] == (p_prim_cr < bonf_thr)).all(), "stored bonferroni_sig differs from p < alpha/L003 on a primary row")

    # ---- the family grid, in a fixed order
    row_of = {k: i for i, k in enumerate(zip(txt.gene[prim_mask], txt.outcome_name[prim_mask]))}
    prim_idx = np.flatnonzero(prim_mask)
    grid = [(g, o) for g in inst_genes for o in outcomes]
    grid_pos = {k: n for n, k in enumerate(grid)}
    is_test = np.array([(g, o) in row_of for g, o in grid])
    _check("S06", len(grid) == fam and int(is_test.sum()) == n_test, "family grid does not contain every testable pair")
    p_fam_cr = np.ones(fam) * float(pad)
    p_fam_pd = np.ones(fam) * float(pad)
    for n, (g, o) in enumerate(grid):
        if is_test[n]:
            j = prim_idx[row_of[(g, o)]]
            p_fam_cr[n], p_fam_pd[n] = pv_cr[j], pv_pd[j]

    # ---- canonical BH (S11)
    thr, flag = _bh_flags(p_fam_cr, fam, alpha)
    thr_pd, flag_pd = _bh_flags(p_fam_pd, fam, alpha)
    q = _bh_q(p_fam_cr)
    _check("S11", int(flag.sum()) == int(LV["L008"]), "canonical BH admits %d, ledger L008 is %s" % (int(flag.sum()), LV["L008"]))
    _check("S11", _bh_count_ledger(p_fam_pd[is_test], fam, alpha) == int(LV["L008"]), "ledger-style bh_count differs from L008")
    _check("S11", ((q < alpha) == flag).all(), "flag differs from BH q < alpha (R p.adjust form)")
    _check("S11", (flag == flag_pd).all() and _close(thr, thr_pd), "flags depend on how the p text is parsed")
    _check("S11", not flag[~is_test].any(), "an untestable combination is flagged")

    # ---- stored (legacy) flag (S12) and L401 (S13)
    _check("S12", set(txt.fdr_sig) <= {"TRUE", "FALSE"}, "fdr_sig is not TRUE/FALSE text")
    legacy = (txt.fdr_sig == "TRUE").values
    calc = np.zeros(len(txt), bool)
    for meth in txt.method.unique():
        ix = (txt.method == meth).values
        calc[ix] = _bh_q(pv_cr[ix]) < alpha
    _check("S12", (calc == legacy).all(), "stored fdr_sig is not BH within each estimator (01_mr_screen.R fdr_sig line)")
    _check("S12", int(legacy[prim_mask].sum()) == int(LV["L009"]), "stored fdr_sig over primary rows %d != L009 %s" % (int(legacy[prim_mask].sum()), LV["L009"]))
    leg_prim = legacy[prim_mask]
    canon_prim = flag[[grid_pos[k] for k in pairs]]                    # canonical flag of each primary row, in screen row order
    L401 = LV["L401"]
    extra = leg_prim & ~canon_prim
    _check("S13", _close(thr, L401["max_p_admitted"]), "largest admitted p %r != L401 %r" % (thr, L401["max_p_admitted"]))
    _check("S13", int(leg_prim.sum()) == int(L401["stored"]) and int(extra.sum()) == int(L401["extra"]) and not (canon_prim & ~leg_prim).any(),
           "stored / extra counts differ from L401")
    _check("S13", _close(p_prim_cr[extra].min(), L401["extra_p_min"]) and _close(p_prim_cr[extra].max(), L401["extra_p_max"]),
           "p range of the extra tests differs from L401")
    _check("S14", (q[~is_test] == 1.0).all() and not flag[~is_test].any(), "an untestable combination has q != 1 or a flag")

    # ---- non-primary rows (S15)
    sur = txt[~prim_mask]
    _check("S15", len(sur) == len(txt) - n_test, "surplus row count")
    _check("S15", all(((g, o) in row_of and txt.method.values[prim_idx[row_of[(g, o)]]] == IVW) for g, o in zip(sur.gene, sur.outcome_name)),
           "a non-primary row belongs to a pair without a primary IVW row")
    _check("S15", set(sur.method) <= {"MR Egger", "Weighted median"} and (sur.fdr_sig == "FALSE").all(), "unexpected non-primary rows")

    # ---- run log: reasons for the untestable combinations (S20-S23)
    lg = _read_log(v3)
    _check("S20", lg["grid_status"] == "COMPLETE" and lg["grid_num"] == lg["grid_den"] == fam, "log grid accounting %s %d/%d vs L003 %d" % (lg["grid_status"], lg["grid_num"], lg["grid_den"], fam))
    _check("S20", lg["with_result"] == n_test, "log with-result %d != testable %d" % (lg["with_result"], n_test))
    _check("S20", lg["with_result"] + lg["attempted_no_result"] + lg["no_outcome_data"] + lg["null_exposure"] == fam, "log cell classes do not sum to L003")
    _check("S20", lg["attempted_no_result"] + lg["no_outcome_data"] == int(LV["L400"]["untestable"]) == n_unt, "log untestable cells != L400")
    _check("S20", lg["null_exposure"] == 0, "the log records null-exposure cells; this module's reason scheme does not cover them")
    orow = {r[0]: r for r in lg["outcome_rows"]}
    _check("S21", len(orow) == len(lg["outcome_rows"]) == len(outcomes) and set(orow) == set(id_of.values()), "log per-outcome table does not list the screen's outcomes")
    test_by_o = {o: sum(1 for g in inst_genes if (g, o) in row_of) for o in outcomes}
    for o in outcomes:
        _, attempted, with_result, no_result = orow[id_of[o]]
        _check("S21", with_result == test_by_o[o] and attempted - with_result == no_result and attempted + (len(inst_genes) - attempted) == len(inst_genes),
               "log row for %s: attempted %d / with-result %d / no-result %d vs %d testable in the screen" % (o, attempted, with_result, no_result, test_by_o[o]))
    _check("S21", sum(r[3] for r in lg["outcome_rows"]) == lg["attempted_no_result"] and
           sum(len(inst_genes) - r[1] for r in lg["outcome_rows"]) == lg["no_outcome_data"], "log per-outcome sums differ from the log's cell totals")
    _check("S22", lg["n_proteins"] == int(LV["L001"]) and lg["rows"] == len(txt) and lg["primary"] == n_test and lg["genes"] == txt.gene.nunique(),
           "log run summary differs from the screen")
    _check("S22", lg["bonferroni"] == int(LV["L007"]) and lg["fdr_by_method"] == int(LV["L009"]) and lg["suggestive"] == int(LV["L006"]),
           "log Bonferroni / FDR(by=method) / suggestive differ from L007 / L009 / L006")
    _check("S22", _close(float(lg["bonf_text"]), bonf_thr, rel=1e-3), "log Bonferroni threshold differs from alpha/L003")
    attempted = set(lg["cells"])
    analysed = {k for k, c in lg["cells"].items() if c["a"] > 0}
    primary_cells = {(g, id_of[o]) for g, o in pairs}
    _check("S23", lg["orphan_analysing"] == 0, "run log: %d analysis line(s) follow no harmonisation line of their cell" % lg["orphan_analysing"])
    _check("S23", analysed == primary_cells,
           "run log: cells with an analysis line differ from the screen's testable cells (%d only in the log, %d only in the screen)"
           % (len(analysed - primary_cells), len(primary_cells - analysed)))
    anr = attempted - analysed                                   # attempted, harmonised, no analysis line: attempted-no-result
    grid_cells = {(g, id_of[o]) for g in inst_genes for o in outcomes}
    _check("S23", attempted <= grid_cells and len(attempted) == lg["with_result"] + lg["attempted_no_result"] and len(anr) == lg["attempted_no_result"],
           "run log: attempted cells (%d) / attempted-no-result cells (%d) differ from the log's own totals" % (len(attempted), len(anr)))
    for o in outcomes:
        _, att_o, _, nores_o = orow[id_of[o]]
        _check("S23", sum(1 for k in attempted if k[1] == id_of[o]) == att_o and sum(1 for k in anr if k[1] == id_of[o]) == nores_o,
               "run log: harmonised / unanalysed cells for %s differ from the per-outcome table" % o)
    nod_line, att_line = _script_lines(v3)

    # ---- write
    est_src = [c for c in STORED_COLUMNS if c not in ID_COLUMNS]
    pad_text = ("%d" % pad) if float(pad) == int(pad) else repr(float(pad))
    rows = []
    for n, (g, o) in enumerate(grid):
        if is_test[n]:
            r = txt.iloc[prim_idx[row_of[(g, o)]]]
            rows.append([g, o, id_of[o], "TRUE", ""] + [r[c] for c in est_src] + [r["pval"], repr(float(q[n])), _tf(flag[n])])
        else:
            rows.append([g, o, id_of[o], "FALSE", _reason_code((g, id_of[o]) in anr)] + [""] * len(est_src) + [pad_text, repr(float(q[n])), _tf(flag[n])])
    _write_csv(out / OUT_LEDGER, LEDGER_HEADER, rows)

    srows = sorted(([r["gene"], r["outcome_name"], r["outcome_id"], NONPRIMARY_ROLE] + [r[c] for c in est_src] for _, r in sur.iterrows()),
                   key=lambda x: (x[0], x[1], x[ESTIMATE_COLUMNS.index("method") + 4]))
    _write_csv(out / OUT_NONPRIMARY, NONPRIMARY_HEADER, srows)

    ev_log = Path(LOG).name
    pad_def = "stays in the multiple-testing family at p = %s" % pad_text
    anr_cells = sorted(anr)
    anr_ev = ", ".join("%s:%d (%s %s)" % (ev_log, lg["cells"][k]["line"], k[0], next(o for o in outcomes if id_of[o] == k[1])) for k in anr_cells)
    base_ev = "%s:%d (class totals); %s:%d-%d (per-outcome table); 01_mr_screen.R:%d-%d (cell classes)" % (
        ev_log, lg["cells_line"], ev_log, lg["outcome_lines"][0], lg["outcome_lines"][1], nod_line, att_line)
    by_code = {CODE_NOD: [r for r in rows if r[3] == "FALSE" and r[4] == CODE_NOD], CODE_ANR: [r for r in rows if r[3] == "FALSE" and r[4] == CODE_ANR]}
    codes = [
        [CODE_NOD, str(len(by_code[CODE_NOD])), ";".join(sorted({r[1] for r in by_code[CODE_NOD]})),
         "The outcome extract held no row for any of the protein's instrument variants, directly or through an LD proxy, so the pipeline did not "
         "attempt the combination (pipeline cell class no-outcome-data; no harmonisation step for it in the run log). The combination " + pad_def + ".",
         base_ev],
        [CODE_ANR, str(len(by_code[CODE_ANR])), ";".join(sorted({r[1] for r in by_code[CODE_ANR]})),
         "The pipeline attempted the combination (at least one instrument variant was in the outcome extract) but no estimate resulted: the run log "
         "shows the harmonisation step for the combination and no analysis step after it (pipeline cell class attempted-no-result). The combination "
         + pad_def + ".",
         base_ev + "; harmonisation without analysis at " + anr_ev],
    ]
    _write_csv(out / OUT_CODES, CODES_HEADER, codes)

    shutil.copyfile(v3 / INSTR, out / OUT_INSTR)

    # ---- read back the shipped representation (S30-S33)
    hdr, back = _read_csv_rows(out / OUT_LEDGER)
    _check("S32", hdr == LEDGER_HEADER and len(back) == fam, "ledger file: header / row count differ after writing")
    ix = {c: i for i, c in enumerate(hdr)}
    keys = [(r[ix["gene"]], r[ix["outcome_name"]]) for r in back]
    _check("S32", keys == sorted(keys) and len(set(keys)) == fam, "ledger file: rows are not unique and sorted")
    pf = np.array([float(r[ix["p_family"]]) for r in back])
    qf = np.array([float(r[ix["fdr_q_bh_family"]]) for r in back])
    thr_b, flag_b = _bh_flags(pf, fam, alpha)
    for r, g_o in zip(back, keys):
        if r[ix["testable"]] == "TRUE":
            s = txt.iloc[prim_idx[row_of[g_o]]]
            _check("S32", all(r[ix[LEGACY if c == "fdr_sig" else c]] == s[c] for c in est_src) and r[ix["outcome_id"]] == s["outcome_id"]
                   and r[ix["p_family"]] == s["pval"] and r[ix["reason_code"]] == "", "ledger file: a testable row differs from the screen text")
        else:
            _check("S32", all(r[ix[c]] == "" for c in ESTIMATE_COLUMNS) and r[ix["reason_code"]] != "" and r[ix["p_family"]] == pad_text,
                   "ledger file: an untestable row carries an estimate, no reason code or a wrong pad")
    _check("S32", (qf == _bh_q(pf)).all() and ((np.array([r[ix["fdr_sig_bh_family"]] == "TRUE" for r in back])) == flag_b).all(),
           "ledger file: BH columns differ when recomputed from the written p_family")
    _check("S32", int(flag_b.sum()) == int(LV["L008"]) and sum(r[ix["testable"]] == "TRUE" for r in back) == int(LV["L004"]) and
           sum(r[ix[LEGACY]] == "TRUE" for r in back) == int(LV["L009"]) and sum(r[ix["testable"]] == "FALSE" for r in back) == int(LV["L400"]["untestable"]),
           "ledger file: counts differ from L008 / L004 / L009 / L400")
    for o in outcomes:                                   # expectation re-derived from the log's per-outcome table
        _, attempted_o, _, no_result_o = orow[id_of[o]]
        have = [r[ix["reason_code"]] for r, k in zip(back, keys) if k[1] == o and r[ix["testable"]] == "FALSE"]
        _check("S30", have.count(CODE_ANR) == no_result_o and have.count(CODE_NOD) == len(inst_genes) - attempted_o and len(have) == len(inst_genes) - orow[id_of[o]][2],
               "reason codes for %s differ from the run log's per-outcome table" % o)
    chdr, crows = _read_csv_rows(out / OUT_CODES)
    _check("S30", chdr == CODES_HEADER and sum(int(c[1]) for c in crows) == int(LV["L400"]["untestable"]), "codebook cells do not sum to L400 untestable")
    nhdr, nrows = _read_csv_rows(out / OUT_NONPRIMARY)
    _check("S15", nhdr == NONPRIMARY_HEADER and len(nrows) == len(txt) - n_test, "non-primary file row count")
    _scan_tpmi([out / OUT_LEDGER, out / OUT_NONPRIMARY, out / OUT_CODES, out / OUT_INSTR])
    _check("S33", _md5(out / OUT_INSTR) == _md5(v3 / INSTR), "instrument copy is not byte-identical to its source")
    ihdr, irows = _read_csv_rows(out / OUT_INSTR)

    return [
        dict(release_path=OUT_LEDGER, ship_status="ship", source_paths=[SCREEN, INSTR, LOG], producer=PRODUCER, licence_class=lic,
             rows=len(back), columns=len(hdr),
             note=("Family of %d protein x outcome combinations (L003): %d testable (stored primary estimate; Wald ratio or IVW), %d untestable at p = %s with a reason_code. "
                   "fdr_sig renamed %s (stored, %d). fdr_q_bh_family / fdr_sig_bh_family are derived for the release: Benjamini-Hochberg over the whole family "
                   "at alpha = ledger L239; they reproduce L008 (%d) and L401 (largest admitted p). Stored columns are the screen's own text. "
                   "Reason codes are read from the pipeline run log (cell-level), see screen_reason_codes.csv. Parameters alpha and pad read from ledger rows L239, L400."
                   % (fam, n_test, n_unt, pad_text, LEGACY, int(LV["L009"]), int(LV["L008"])))),
        dict(release_path=OUT_NONPRIMARY, ship_status="ship", source_paths=[SCREEN], producer=PRODUCER, licence_class=lic,
             rows=len(nrows), columns=len(nhdr),
             note=("Labelled, not part of the family: the %d Egger / weighted-median rows of the multi-instrument pair(s) %s (stored columns, same renames). "
                   "Each belongs to a pair whose primary IVW estimate is in screen_ledger.csv." % (len(nrows), ",".join(sorted(set(sur.gene)))))),
        dict(release_path=OUT_CODES, ship_status="ship", source_paths=[LOG, SCREEN, INSTR, CASCADE01], producer=PRODUCER, licence_class="own_derived",
             rows=len(crows), columns=len(chdr),
             note=("Codebook for reason_code. Cell classes are read from %s: an untestable combination with no harmonisation line is no_outcome_data, one with a "
                   "harmonisation line and no analysis line is attempted_no_result. Both agree with the log's own totals and per-outcome table, and the cells with an "
                   "analysis line are exactly the screen's testable cells. The pipeline's per-cell audit files are not part of this release." % Path(LOG).name)),
        dict(release_path=OUT_INSTR, ship_status="ship", source_paths=[INSTR], producer=PRODUCER, licence_class="third_party_derived:" + EXPOSURE_PROVIDER,
             rows=len(irows), columns=len(ihdr),
             note=("Instrument list as stored (byte-identical copy): per-variant beta, se, p and eaf from UKB-PPP, so the licence class is third-party-derived. "
                   "The terms of UKB-PPP are to be read in step 4; shipping follows spec section 4.")),
    ]


def build(v3, out, withheld, ledger):
    """Interface of spec section 2. Builds in a private staging folder and copies the four files into out/ only when every check has
    passed, so a refused build leaves nothing behind. `withheld` is accepted and not used: nothing here is held back."""
    stage = Path(tempfile.mkdtemp(prefix="rel_screen_"))
    try:
        recs = _build(v3, stage, ledger)
        for rec in recs:
            dst = Path(out) / rec["release_path"]
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(stage / rec["release_path"], dst)
        return sorted(recs, key=lambda r: r["release_path"])
    finally:
        shutil.rmtree(stage, ignore_errors=True)


# ------------------------------------------------------------------------------------------------ selftest
def _csv_rewrite(b, fn):
    """Apply fn(header, rows) -> rows to a CSV held as bytes (CRLF, minimal quoting, as the source tables are).
    Asserts first that the csv round trip is the identity, so the mutation is the only difference."""
    rd = list(csv.reader(io.StringIO(b.decode("utf-8"), newline="")))

    def dump(rows):
        sio = io.StringIO(newline="")
        csv.writer(sio, lineterminator="\r\n").writerows(rows)
        return sio.getvalue().encode("utf-8")

    assert dump(rd) == b, "csv round trip is not the identity; the mutation would confound"
    return dump([rd[0]] + fn(rd[0], rd[1:]))


def _csv_set(b, where, col, new):
    """Set one cell of the single row selected by `where` (a dict predicate); the target must be unique."""
    def fn(h, rows):
        hits = [r for r in rows if where(dict(zip(h, r)))]
        assert len(hits) == 1, "mutation target is not unique (%d rows)" % len(hits)
        hits[0][h.index(col)] = new
        return rows
    return _csv_rewrite(b, fn)


def _lines_drop(b, col, value):
    """Drop the CRLF lines whose comma-split field `col` equals value (for a CSV with quoted empty cells, where a csv round trip is not exact)."""
    crlf = bytes([13, 10])
    lines = b.split(crlf)
    keep = [l for l in lines if l.split(b",")[col:col + 1] != [value]]
    assert len(keep) < len(lines), "mutation dropped nothing"
    return crlf.join(keep)


def _drop_line(b, text):
    """Remove the one line equal to `text` (any line ending)."""
    pat = re.escape(text.encode("utf-8")) + rb"\r?\n"
    assert len(re.findall(pat, b)) == 1, "mutation target is not unique"
    return re.sub(pat, b"", b, count=1)


def _replace_once(b, old, new):
    assert b.count(old) == 1, "mutation target is not unique (%d)" % b.count(old)
    return b.replace(old, new)


def _scratch(v3, tmp, edits, refresh=True):
    """Copy INPUTS to tmp (same relative paths) and apply byte edits. With refresh, RDS_CONVERSIONS is brought in line with an edited
    table (a re-converted file); without it, the edit is a silent change to a converted table, which S02 must catch."""
    tmp = Path(tmp)
    for rel in INPUTS:
        (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(Path(v3) / rel, tmp / rel)
    for rel, fn in edits.items():
        b = (tmp / rel).read_bytes()
        nb = fn(b)
        assert nb != b, "mutation changed nothing"
        (tmp / rel).write_bytes(nb)
    if refresh:
        cv = pd.read_csv(tmp / CONV, dtype=str, keep_default_na=False)
        for rel in (SCREEN, INSTR):
            if rel in edits:
                k = cv.csv == "tables_from_rds/" + Path(rel).name
                cv.loc[k, "csv_md5"] = _md5(tmp / rel)
                cv.loc[k, "n_rows"] = str(len(pd.read_csv(tmp / rel, dtype=str, keep_default_na=False)))
        cv.to_csv(tmp / CONV, index=False, lineterminator="\r\n")
    return tmp


def _ledger_with(ledger, lid, fn):
    l2 = ledger.copy()
    k = l2["id"] == lid
    assert k.sum() == 1
    v = fn(json.loads(l2.loc[k, "value"].iloc[0]))
    l2.loc[k, "value"] = json.dumps(v)
    return l2


def _trip(v3, ledger, patch=None):
    """Build into a throw-away folder; return the id of the check that fired, or 'no failure'."""
    g = globals()
    saved = {k: g[k] for k in (patch or {})}
    try:
        g.update(patch or {})
        with tempfile.TemporaryDirectory() as t:
            build(v3, t, Path(t) / "withheld", ledger)
        return "no failure"
    except ScreenCheckError as e:
        return e.check_id
    finally:
        g.update(saved)


def _wrong_bh_flags(p, m, alpha):
    """Defect planted at the layer where it would live: BH over the analysable rows only (the legacy denominator)."""
    p = np.asarray(p, float)
    keep = np.sort(p[p < 1])
    thr, _ = _ORIG_BH_FLAGS(keep, len(keep), alpha)
    return thr, p <= thr


def _wrong_reason(attempted_without_result):
    """Defect at the labelling layer: every untestable cell is called no_outcome_data."""
    return CODE_NOD


def _write_rounded_q(path, header, rows):
    """Defect at the writer: the q column is written to six significant digits."""
    if header == LEDGER_HEADER:
        k = header.index("fdr_q_bh_family")
        rows = [r[:k] + ["%.6g" % float(r[k])] + r[k + 1:] for r in rows]
    _REAL_WRITE_CSV(path, header, rows)


def _write_changed_cell(path, header, rows):
    """Defect at the writer: one stored cell of the first testable row is changed."""
    if header == LEDGER_HEADER:
        k, t = header.index("b"), header.index("testable")
        i = next(i for i, r in enumerate(rows) if r[t] == "TRUE")
        rows = [list(r) for r in rows]
        rows[i][k] = "9"
    _REAL_WRITE_CSV(path, header, rows)


def _capped_q(p):
    return np.minimum(_ORIG_BH_Q(p), 0.999)


def selftest(v3):
    """Baseline build twice (byte-identical), then one mutation per check. Returns [(name, outcome, ok)]; every mutation must raise
    exactly the check it names, with a unique target. A mutation that trips a different check is reported as a failure."""
    v3 = Path(v3)
    ledger = pd.read_csv(v3 / LEDGER_CSV, dtype=str)
    rep = []
    with tempfile.TemporaryDirectory() as t1, tempfile.TemporaryDirectory() as t2, tempfile.TemporaryDirectory() as ts:
        r1 = build(v3, t1, Path(t1) / "withheld", ledger)
        build(v3, t2, Path(t2) / "withheld", ledger)
        for a in r1:
            same = _md5(Path(t1) / a["release_path"]) == _md5(Path(t2) / a["release_path"])
            rep.append(("two builds, " + a["release_path"], "byte-identical" if same else "DIFFERENT", same))
        txt = pd.read_csv(v3 / SCREEN, dtype=str, keep_default_na=False)
        num = pd.read_csv(v3 / SCREEN)
        prim = num[num.method.isin([WALD, IVW])]
        mid = prim[(prim.pval > 0.0012) & (prim.pval < 0.04)].iloc[0]         # not admitted; nominal and Bonferroni counts stay put
        hi = prim[(prim.method == WALD) & (prim.pval > 0.5)].iloc[0]
        sur = num[~num.method.isin([WALD, IVW])].iloc[0]
        ins = pd.read_csv(v3 / INSTR, dtype=str, keep_default_na=False)
        gone = sorted(set(ins.gene) - set(txt.gene))[0]
        log = _text(v3 / LOG)
        aq = re.search(r"^ *\d+: +FinnGen_R10_I9_AF +\d+ +\d+ +\d+$", log, re.M).group(0)
        an_line = next(l for l in log.splitlines() if RE_ANALYSING.match(l) and log.count(l) == 1)
        hz_line = next(l for l in log.splitlines() if RE_HARM.match(l) and "AMIGO2" in l and "Wuttke2019_CKD_EA" in l and log.count(l) == 1)
        aq_bad = re.sub(r"(FinnGen_R10_I9_AF +)\d+", lambda m: m.group(1) + "9999", aq, count=1)
        at = lambda row: (lambda d: d["gene"] == row.gene and d["outcome_id"] == row.outcome_id and d["method"] == row.method)
        n = [0]

        def mut(name, expect, edits=None, ledger_edit=None, patch=None, refresh=True):
            n[0] += 1
            tmp = Path(ts) / ("m%d" % n[0])
            tmp.mkdir()
            sc = _scratch(v3, tmp, edits or {}, refresh)
            got = _trip(sc, ledger if ledger_edit is None else ledger_edit, patch)
            rep.append((name, "expected %s, got %s" % (expect, got), got == expect))

        def bump(m):
            return m.group(1) + str(int(m.group(2)) + 1).encode() + m.group(3)

        mut("ledger L400 row removed", "S00", ledger_edit=ledger[ledger["id"] != "L400"])
        mut("ledger L003 + 1", "S06", ledger_edit=_ledger_with(ledger, "L003", lambda v: v + 1))
        mut("ledger L400 untestable + 1", "S07", ledger_edit=_ledger_with(ledger, "L400", lambda v: dict(v, untestable=v["untestable"] + 1)))
        mut("ledger L008 - 1", "S11", ledger_edit=_ledger_with(ledger, "L008", lambda v: v - 1))
        mut("ledger L401 max_p_admitted x 1.01", "S13", ledger_edit=_ledger_with(ledger, "L401", lambda v: dict(v, max_p_admitted=v["max_p_admitted"] * 1.01)))
        mut("config: padding sentence removed", "S01", edits={CONFIG: lambda b: _replace_once(b, b"untested cells padded at", b"untested cells left at")})
        mut("screen edited, conversion record not refreshed", "S02", edits={SCREEN: lambda b: _csv_set(b, at(mid), "se", "0.1234")}, refresh=False)
        mut("screen: a non-primary row relabelled as a primary one", "S03", edits={SCREEN: lambda b: _csv_set(b, at(sur), "method", WALD)})
        mut("screen: a gene renamed outside the instrument list", "S05", edits={SCREEN: lambda b: _csv_set(b, at(mid), "gene", "ZZZ9")})
        mut("instruments: a gene of the family removed", "S05", edits={INSTR: lambda b: _lines_drop(b, list(ins.columns).index("gene"), gone.encode())})
        mut("screen: a Wald row with n_snps 2", "S08", edits={SCREEN: lambda b: _csv_set(b, at(hi), "n_snps", "2")})
        mut("screen: a p moved across 0.05", "S09", edits={SCREEN: lambda b: _csv_set(b, at(mid), "pval", "0.06")})
        mut("screen: stored bonferroni_sig flipped", "S10", edits={SCREEN: lambda b: _csv_set(b, at(hi), "bonferroni_sig", "TRUE")})
        mut("screen: one p moved under the BH threshold only", "S11", edits={SCREEN: lambda b: _csv_set(b, at(mid), "pval", "0.0005")})
        mut("BH over the analysable rows only", "S11", patch={"_bh_flags": _wrong_bh_flags})
        mut("screen: stored fdr_sig flipped", "S12", edits={SCREEN: lambda b: _csv_set(b, at(hi), "fdr_sig", "TRUE")})
        mut("screen: a non-primary row moved to a pair without IVW", "S15", edits={SCREEN: lambda b: _csv_set(b, at(sur), "gene", "A1BG")})
        mut("log: no-outcome-data total + 1", "S20", edits={LOG: lambda b: re.sub(rb"(\d+ attempted-no-result, )(\d+)( no-outcome-data)", bump, b, count=1)})
        mut("log: one outcome's attempted count changed", "S21", edits={LOG: lambda b: _replace_once(b, aq.encode(), aq_bad.encode())})
        mut("log: FDR(by=method) + 1", "S22", edits={LOG: lambda b: re.sub(rb"(FDR\(by=method\)=)(\d+)()", bump, b, count=1)})
        mut("reason code: the attempted-no-result cells labelled no_outcome_data", "S30", patch={"_reason_code": _wrong_reason})
        mut("log: an analysis line of a testable cell removed", "S23", edits={LOG: lambda b: _drop_line(b, an_line)})
        mut("log: the harmonisation line of an unanalysed cell removed", "S23", edits={LOG: lambda b: _drop_line(b, hz_line)})
        mut("screen: a column renamed", "S01", edits={SCREEN: lambda b: _replace_once(b, b",Q_df,", b",Qdf,")})
        mut("screen: an outcome_id borrowed from another outcome", "S04", edits={SCREEN: lambda b: _csv_set(b, at(hi), "outcome_id", sorted(set(txt.outcome_id) - {hi.outcome_id})[0])})
        mut("ledger L401 extra + 1", "S13", ledger_edit=_ledger_with(ledger, "L401", lambda v: dict(v, extra=v["extra"] + 1)))
        mut("untestable q capped below 1", "S14", patch={"_bh_q": _capped_q})
        mut("log: genes=... changed", "S22", edits={LOG: lambda b: re.sub(rb"(genes=)(\d+)()", bump, b, count=1)})
        mut("writer: q column rounded", "S32", patch={"_write_csv": _write_rounded_q})
        mut("writer: a stored cell changed", "S32", patch={"_write_csv": _write_changed_cell})
        bad = Path(ts) / "bad.csv"
        bad.write_text("a,b\nx,planted TPMI value\n", encoding="utf-8")
        try:
            _scan_tpmi([bad])
            got = "no failure"
        except ScreenCheckError as e:
            got = e.check_id
        rep.append(("TPMI string planted in a shipped file", "expected S31, got %s" % got, got == "S31"))
        real_copy = shutil.copyfile
        try:
            shutil.copyfile = lambda s, d: Path(d).write_bytes(Path(s).read_bytes()[:-1])
            got = _trip(v3, ledger)
        finally:
            shutil.copyfile = real_copy
        rep.append(("instrument copy truncated by one byte", "expected S33, got %s" % got, got == "S33"))
    return rep


_ORIG_BH_FLAGS = _bh_flags
_ORIG_BH_Q = _bh_q
_REAL_WRITE_CSV = _write_csv


def _main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--v3", default=str(Path(__file__).resolve().parents[2]))
    ap.add_argument("--build", nargs=2, metavar=("OUT", "WITHHELD"), help="build into OUT / WITHHELD (refuses this project's release/out)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    v3 = Path(a.v3)
    if a.build:
        out = Path(a.build[0]).resolve()
        protected = (Path(__file__).resolve().parent / "out").resolve()
        if out == protected or protected in out.parents:
            raise SystemExit("refusing to write into %s: the core builder owns it" % protected)
        for rec in build(v3, out, a.build[1], pd.read_csv(v3 / LEDGER_CSV, dtype=str)):
            print(json.dumps(rec))
    if a.selftest:
        bad = 0
        for name, msg, ok in selftest(v3):
            print(("ok   " if ok else "FAIL ") + name + " -> " + msg)
            bad += (not ok)
        raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    _main()
