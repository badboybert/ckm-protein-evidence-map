# -*- coding: utf-8 -*-
# R2 Stage 3 join: ledger rows L700-L998 (ceiling lowered from L999 at Stage 5 step 1a, S2-D27: L999 and L971-L984 are reserved for the self-tests) for every {{id}} placeholder used in the R2 Results / Table 3 / SM / Discussion drafts.
# Executed by S1R1_build_ledger.py in ITS namespace (add, rel, ROWS, EXT_CHECKS, HERE, V3, CKM and the stdlib modules it imported).
# Every name defined here starts with _r2s_ so nothing in the builder's namespace is overwritten.
#
# INPUT   R2_ext_ids.txt (same folder): "id <TAB> package <TAB> used_in", '#' comments. APPEND-ONLY: add new ids at the end.
# OUTPUT  ledger rows (one row per distinct source CELL = package + quantity_id + column; several ids may share a cell) and
#         R2_id_to_L_map.tsv (id -> Lxxx -> value as printed -> source). The previous map is read back so that an id keeps its
#         Lxxx for ever; a new id takes the next free number.  Existing rows are never renumbered or reshaped.
# VALUES  read from the Stage 1 package results.tsv (or, for a Stage 1b key id, from the KEY_MAP quantity and the KEY_IDS
#         compare_col). Nothing is typed here: every number is parsed from a source file and every printed form is derived
#         from that source string by the print rules below (decimal arithmetic on the source text, half-up).
# PRINT RULES (recommended precision; the gate N2 accepts any rounding of the full-precision value):
#   p      2 significant figures (p < 0.001 as m x 10^e); exactly 1 prints "1"; three significant figures for the RAD23B
#          study-wide-cut family (S1H.lead_p, S1H.proxy_p, S1H.p_second_order) because 0.05/15,584 = 3.21 x 10^-6 is compared with them
#   pp     posterior probabilities 2 decimals (2 significant figures below 0.01)
#   dp3    s, r2, log-odds ratios and differences, ratio of ratios: 3 decimals (2 significant figures below 0.0095)
#   dec2   chi-square and z: 2 decimals
#   share  bootstrap shares: 2 decimals
#   or     odds ratios and their bounds: 3 significant figures (4 when >= 100 or when 3 would print 1.00)
#   int    counts
# ADDED 2026-10-14 (infrastructure round, append-only; the previous Lxxx of every id are unchanged):
#   * R2_ext_ids.txt column 4 (optional) = ALIAS 'quantity_id|column': the id is another name for a numeric cell of the package results.tsv.
#   * R2_protocol_constants.tsv = protocol constants and engine settings that no results.tsv holds (thresholds, dates, seeds, panel sizes, engine
#     settings). One row each, value parsed from a token quoted from a PINNED document/script; document md5, line and quote asserted at every build.
#     Their cells ('CONST', id, 'value') take the next free Lxxx after the results cells, in file order.
#   * A KEY row whose Stage 1b status was scored against a superseded value is labelled KEY-STALE / KEY-VISIBLE (not plain KEY); see _r2s_key_status.
import csv as _r2s_csv, io as _r2s_io, os as _r2s_os, re as _r2s_re, sys as _r2s_sys, hashlib as _r2s_hashlib
from decimal import Decimal as _r2s_D, ROUND_HALF_UP as _r2s_HU, ROUND_HALF_EVEN as _r2s_HE

_r2s_csv.field_size_limit(10 ** 9)
_r2s_R1 = _r2s_os.path.join(V3, "R1_stage1")
_r2s_IDS = _r2s_os.path.join(HERE, "R2_ext_ids.txt")
_r2s_MAP = _r2s_os.path.join(HERE, "R2_id_to_L_map.tsv")
_r2s_LO, _r2s_HI = 700, 998   # L999 and L971-L984 are reserved for the self-tests (S1R1_id_alloc.py); the allocator never returns them
_r2s_sys.path.insert(0, HERE)   # python -I does not put the script's folder on sys.path
import S1R1_id_alloc as _r2s_ida  # noqa: E402
_r2s_CUT3 = {"S1H.lead_p", "S1H.proxy_p", "S1H.p_second_order"}   # three significant figures: compared with the Bonferroni cut
_r2s_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
_r2s_NA = ("", "NA", "NaN", "nan", "NULL")


def _r2s_md5(path):
    h = _r2s_hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _r2s_read_tsv(path):
    """-> (header, [(line_no, dict)]). csv module with tab delimiter: handles the quoted files (04b, 05, 10) and the unquoted ones."""
    with _r2s_io.open(path, encoding="utf-8", newline="") as fh:
        rd = _r2s_csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        out = []
        for r in rd:
            if not r:
                continue
            assert len(r) == len(hdr), "%s line %d: %d fields, header has %d" % (path, rd.line_num, len(r), len(hdr))
            out.append((rd.line_num, dict(zip(hdr, r))))
    return hdr, out


# ------------------------------------------------------------------------------------------------------------------ print rules
def _r2s_round_sig(d, sf):
    """Decimal -> Decimal rounded half-up to sf significant digits; second return = True if half-up and half-even disagree (a tie)."""
    if d == 0:
        return d, False
    e = d.adjusted()
    q = _r2s_D(1).scaleb(e - sf + 1)
    up, ev = d.quantize(q, rounding=_r2s_HU), d.quantize(q, rounding=_r2s_HE)
    tie = up != ev
    if up.adjusted() != e:   # carry (9.96 -> 10.0): re-quantize so the digit count stays sf
        up = up.quantize(_r2s_D(1).scaleb(up.adjusted() - sf + 1), rounding=_r2s_HU)
    return up, tie


def _r2s_fixed(d, dp):
    q = _r2s_D(1).scaleb(-dp)
    return d.quantize(q, rounding=_r2s_HU), d.quantize(q, rounding=_r2s_HE) != d.quantize(q, rounding=_r2s_HU)


def _r2s_sci(d, sf):
    """-> (mantissa Decimal, exponent int, tie) for a positive Decimal."""
    e = d.adjusted()
    m, tie = _r2s_round_sig(d.scaleb(-e), sf)
    if m >= 10:
        m, e = m / 10, e + 1
        m = m.quantize(_r2s_D(1).scaleb(-(sf - 1)))
    return m, e, tie


_r2s_THR = {"p": [_r2s_D(x) for x in ("0.1", "0.05", "0.01", "0.001")], "pp": [_r2s_D(x) for x in ("0.5", "0.8", "0.9", "0.95")],
            "share": [_r2s_D(x) for x in ("0.5", "0.8", "0.9", "0.95")]}


def _r2s_print(kind, s, sf_override=None):
    """Rounds as the kind says; if the rounded text lands exactly on a decision threshold the value itself does not equal (p 0.0497 -> 0.050,
    PP 0.7996 -> 0.80), one more digit is printed so the rounding cannot move the value across the line."""
    a, pr, tie = _r2s_print0(kind, s, sf_override)
    d = _r2s_D(s.strip())
    if kind in _r2s_THR and any(_r2s_D(a) == th and d != th for th in _r2s_THR[kind] if not a.startswith("1 ×") and "×" not in a):
        a, pr, tie = _r2s_print0(kind, s, (sf_override or 2) + 1, extra_dp=1)
    return a, pr, tie


def _r2s_print0(kind, s, sf_override=None, extra_dp=0):
    """Source string s -> (ascii_display, prose_display, tie_flag). ascii uses '10^-6' and '-'; prose uses superscripts and U+2212."""
    d = _r2s_D(s.strip())
    tie = False
    sup = None
    if kind == "p":
        sf = sf_override or 2
        assert _r2s_D(0) <= d <= _r2s_D(1), "p outside [0, 1]: %s" % s
        if d == 1:
            txt = "1"
        elif d == 0:
            txt = "0"
        elif d < _r2s_D("0.001"):
            m, e, tie = _r2s_sci(d, sf)
            txt, sup = "%s × 10^%d" % (m, e), (str(m), e)
        else:
            r, tie = _r2s_round_sig(d, sf)
            assert r < 1, "p %s rounds to 1" % s
            txt = format(r, "f")
    elif kind == "pp" or kind == "share":
        assert _r2s_D(0) <= d <= _r2s_D(1), "%s outside [0, 1]: %s" % (kind, s)
        if d != 0 and d < _r2s_D("0.01"):
            r, tie = _r2s_round_sig(d, 2)
        else:
            r, tie = _r2s_fixed(d, 2 + extra_dp)
        txt = format(r, "f")
    elif kind == "dp3":
        if d != 0 and abs(d) < _r2s_D("0.0095"):
            r, tie = _r2s_round_sig(d, 2)
        else:
            r, tie = _r2s_fixed(d, 3)
        txt = format(r, "f")
    elif kind == "dec2":
        r, tie = _r2s_fixed(d, 2)
        txt = format(r, "f")
    elif kind == "or":
        assert d > 0, "odds ratio <= 0: %s" % s
        sf = 3
        r3, tie = _r2s_round_sig(d, 3)
        if d >= 100 or abs(r3 - 1) < _r2s_D("0.005"):
            sf = 4
        r, tie = _r2s_round_sig(d, sf)
        txt = format(r, "f")
    elif kind == "int":
        assert d == d.to_integral_value(), "count is not an integer: %s" % s
        txt = "{:,}".format(int(d))
    elif kind == "sig3":   # fallback for a unit with no class rule: 3 significant figures, scientific below 0.001
        if d != 0 and abs(d) < _r2s_D("0.001"):
            m, e, tie = _r2s_sci(abs(d), 3)
            sg = "-" if d < 0 else ""
            txt, sup = "%s%s × 10^%d" % (sg, m, e), (sg + str(m), e)
        else:
            r, tie = _r2s_round_sig(d, 3)
            txt = format(r, "f")
    else:
        raise ValueError(kind)
    if txt.startswith("-") and _r2s_D(txt) == 0:
        raise AssertionError("print rule %s gives a signed zero for %s" % (kind, s))
    ascii_ = txt
    if sup is not None:
        prose = "%s × 10%s" % (sup[0], str(sup[1]).translate(_r2s_SUP))
    else:
        prose = txt
    return ascii_, prose.replace("-", "−"), tie


def _r2s_num(s):
    s = s.strip()
    return int(s) if _r2s_re.fullmatch(r"-?\d+", s) else float(s)


# ------------------------------------------------------------------------------------------------------------------ source index
_r2s_PK = {}      # package dir -> (path, header, rows, md5)
_r2s_BYQ = {}     # package dir -> {quantity_id: [(line, row)]}


def _r2s_pkg(pkg):
    if pkg not in _r2s_PK:
        path = _r2s_os.path.join(_r2s_R1, pkg, "results.tsv")
        assert _r2s_os.path.isfile(path), "no results.tsv for package %s" % pkg
        hdr, rows = _r2s_read_tsv(path)
        key = hdr[0]
        assert key in ("quantity_id", "id"), (pkg, hdr[:3])
        by = {}
        for ln, r in rows:
            by.setdefault(r[key], []).append((ln, r))
        _r2s_PK[pkg] = (path, hdr, rows, _r2s_md5(path))
        _r2s_BYQ[pkg] = by
    return _r2s_PK[pkg]


def _r2s_BYQ_get(pkg, qid):
    _r2s_pkg(pkg)
    return _r2s_BYQ[pkg].get(qid, [])


# Stage 1b keys: KEY_IDS (definition, compare_col, kind), KEY_MAP (results file + quantity), stage1b_key_labels (survival status)
_r2s_KEYS = {}
for _s in ("C1", "S1", "S2", "S3", "S4", "S5", "S6", "S7"):
    _p = _r2s_os.path.join(_r2s_R1, "S1R1_R2_1b_verify", _s, "KEY_IDS.tsv")
    with _r2s_io.open(_p, encoding="utf-8", newline="") as _fh:
        for _r in _r2s_csv.DictReader(_fh, delimiter="\t"):
            _r["_sec"] = _s
            _r["_ids_file"] = _p
            _r["_ids_md5"] = None
            _r2s_KEYS[_r["key_id"]] = _r
    _p = _r2s_os.path.join(_r2s_R1, "S1R1_R2_00_resources", "_restricted", "keymap", _s, "KEY_MAP.tsv")
    with _r2s_io.open(_p, encoding="utf-8", newline="") as _fh:
        for _r in _r2s_csv.DictReader(_fh, delimiter="\t"):
            _r2s_KEYS[_r["key_id"]].update({"map_file": _r["results_file"].replace("\\", "/"), "map_md5": _r["results_md5"],
                                            "map_q": _r.get("quantity_ids") or _r.get("quantity_id(s)")})
_r2s_LAB = {}
with _r2s_io.open(_r2s_os.path.join(_r2s_R1, "S1R1_R2_10_integration", "stage1b_key_labels.tsv"), encoding="utf-8", newline="") as _fh:
    for _r in _r2s_csv.DictReader(_fh, delimiter="\t"):
        _r2s_LAB[_r["key_id"]] = _r
assert set(_r2s_LAB) == set(_r2s_KEYS) and all("map_file" in v for v in _r2s_KEYS.values()), "KEY_IDS / KEY_MAP / labels disagree"


def _r2s_isnum(s):
    try:
        float(s)
        return s.strip() not in _r2s_NA
    except ValueError:
        return False


def _r2s_canon_col(row, want):
    """A cell is (package, quantity_id, column). When the value and p columns hold the same number they are one cell ('value'),
    so a key whose compare_col is p and a direct id that reads the value column share one ledger row and one check status."""
    v, p = row.get("value", ""), row.get("p", "")
    if want in ("value", "p") and _r2s_isnum(v) and _r2s_isnum(p) and float(v) == float(p):
        return "value"
    return want


_r2s_KEY_BY_CELL = {}     # (package, quantity_id, canonical column) -> [key_id]
for _k, _v in _r2s_KEYS.items():
    if _v["map_file"].endswith("/results.tsv") and _r2s_re.fullmatch(r"[A-Za-z0-9_.\-]+", _v["map_q"]):
        _pkg_ = _v["map_file"].split("/")[-2]
        _hits_ = _r2s_BYQ_get(_pkg_, _v["map_q"])
        assert len(_hits_) == 1, "key %s: quantity %s found %d time(s) in %s" % (_k, _v["map_q"], len(_hits_), _pkg_)
        _r2s_KEY_BY_CELL.setdefault((_pkg_, _v["map_q"], _r2s_canon_col(_hits_[0][1], _v["compare_col"])), []).append(_k)


# ------------------------------------------------------------------------------------------------------------------ id list
_r2s_lines = []
for _ln in _r2s_io.open(_r2s_IDS, encoding="utf-8").read().split("\n"):
    if not _ln.strip() or _ln.startswith("#"):
        continue
    _c = _ln.split("\t")
    assert 2 <= len(_c) <= 4, "bad line in R2_ext_ids.txt: %r" % _ln
    _r2s_lines.append((_c[0].strip(), _c[1].strip(), _c[2].strip() if len(_c) > 2 else "", _c[3].strip() if len(_c) > 3 else ""))
_r2s_idl = [x[0] for x in _r2s_lines]
assert len(_r2s_idl) == len(set(_r2s_idl)), "duplicate id in R2_ext_ids.txt: %s" % sorted({i for i in _r2s_idl if _r2s_idl.count(i) > 1})

# ---- protocol constants and engine settings (2026-10-14): R2_protocol_constants.tsv, one ledger row per id, value PARSED from the token
# the file quotes from a pinned governing document / deposited script. At every build the document md5 and the verbatim quote on the
# named line are asserted and the token is re-parsed; a changed source stops the build instead of leaving a stale constant.
_r2s_CONST = _r2s_os.path.join(HERE, "R2_protocol_constants.tsv")
_r2s_WORDN = {w: i for i, w in enumerate("zero one two three four five six seven eight nine ten".split())}


def _r2s_tok(tok):
    """Token as written in the source -> int | float | ISO date string. 90% -> 0.9; 5,000 -> 5000; five -> 5; 1e-5 -> 1e-05."""
    t_ = tok.strip()
    if _r2s_re.fullmatch(r"\d{4}-\d{2}-\d{2}", t_):
        return t_
    if t_.lower() in _r2s_WORDN:
        return _r2s_WORDN[t_.lower()]
    pct = t_.endswith("%")
    t_ = (t_[:-1] if pct else t_).replace(",", "")
    d_ = _r2s_D(t_) / (100 if pct else 1)
    return int(d_) if (_r2s_re.fullmatch(r"-?\d+", t_) and not pct) else float(d_)


_r2s_consts = []
if _r2s_os.path.isfile(_r2s_CONST):
    with _r2s_io.open(_r2s_CONST, encoding="utf-8", newline="") as _fh:
        _rd = _r2s_csv.DictReader(_fh, delimiter="\t")
        assert _rd.fieldnames == ["id", "value", "source document", "line", "verbatim quote", "token", "display", "meaning", "doc_md5", "class"], _rd.fieldnames
        _r2s_consts = list(_rd)
    _cid = [r["id"] for r in _r2s_consts]
    assert len(_cid) == len(set(_cid)), "duplicate id in R2_protocol_constants.tsv"
    assert not (set(_cid) & set(_r2s_idl)), "id in both R2_ext_ids.txt and R2_protocol_constants.tsv: %s" % sorted(set(_cid) & set(_r2s_idl))
    _r2s_idl = _r2s_idl + _cid

# ---- derived rows (2026-10-14, D1 round): R2_derived_rows.py fills _r2d_rows (successors of ledger rows whose display printed numbers their value did not
# store, and counts computed from existing local tables). One ledger row each, appended after the constants; ids are never re-pointed.
_r2s_DERIV = _r2s_os.path.join(HERE, "R2_derived_rows.py")
_r2s_derived = []
if _r2s_os.path.isfile(_r2s_DERIV):
    _r2d_rows = []
    exec(compile(_r2s_io.open(_r2s_DERIV, encoding="utf-8").read(), _r2s_DERIV, "exec"))
    _r2s_derived = list(_r2d_rows)
    _did = [r["id"] for r in _r2s_derived]
    assert len(_did) == len(set(_did)), "duplicate id in R2_derived_rows.py"
    assert not (set(_did) & set(_r2s_idl)), "id in both R2_derived_rows.py and R2_ext_ids.txt / R2_protocol_constants.tsv: %s" % sorted(set(_did) & set(_r2s_idl))
    _r2s_idl = _r2s_idl + _did

# previous map: ids keep their Lxxx for ever
_r2s_prev = {}
if _r2s_os.path.isfile(_r2s_MAP):
    with _r2s_io.open(_r2s_MAP, encoding="utf-8", newline="") as _fh:
        for _r in _r2s_csv.DictReader(_fh, delimiter="\t"):
            _r2s_prev[_r["id"]] = _r
    _gone = sorted(set(_r2s_prev) - set(_r2s_idl))
    assert not _gone, "ids in the previous map are missing from R2_ext_ids.txt (append-only): %s" % _gone[:10]


def _r2s_resolve(pid, pkg_in, alias=""):
    """-> dict(pkg, qid, col, line, row, key_id or None). Fails loudly on anything ambiguous.
    alias = 'quantity_id|column' (4th column of R2_ext_ids.txt): the placeholder id is another name for that cell of the package results.tsv."""
    key = _r2s_KEYS.get(pid)
    if alias:
        assert key is None and not pid.startswith("S12."), "%s: an alias id cannot also be a key / S12 id" % pid
        assert pkg_in not in ("", "auto", "KEYMAP"), "%s: an alias needs an explicit package" % pid
        pkg, qid = pkg_in, alias.split("|")[0]
    elif pid.startswith("S12."):
        pkg, qid = "S1R1_R2_12_chen_instrument", pid[4:]
    elif key is not None:
        pkg, qid = key["map_file"].split("/")[-2], key["map_q"]
        assert key["map_file"].endswith("/results.tsv"), "%s maps to %s, not a results.tsv" % (pid, key["map_file"])
    else:
        qid = pid
        pkg = pkg_in
    if pkg_in not in ("", "auto", "KEYMAP") and not alias:
        assert pkg_in == pkg or (key is None and not pid.startswith("S12.")), "%s: package column %s disagrees with resolved %s" % (pid, pkg_in, pkg)
    path, hdr, rows, md5 = _r2s_pkg(pkg)
    hits = _r2s_BYQ[pkg].get(qid, [])
    assert len(hits) == 1, "%s: quantity %s found %d time(s) in %s" % (pid, qid, len(hits), pkg)
    ln, row = hits[0]
    return dict(pkg=pkg, qid=qid, line=ln, row=row, key=key, path=path, md5=md5, hdr=hdr)


_r2s_PNAME = _r2s_re.compile(r"^(p|p[AB]|p_[A-Za-z0-9_]+|wald[12]_(?:p|lo|hi)|p_wald1)$")


def _r2s_column(pid, hit, alias=""):
    row, key = hit["row"], hit["key"]
    vcol = "value"
    if alias:
        col = alias.split("|")[1]
        assert col in row and _r2s_isnum(row[col]), "%s: alias column %s of %s is not a numeric cell" % (pid, col, hit["qid"])
        return col
    if key is not None:
        col = key["compare_col"]
        assert col in row, "%s: compare_col %s is not a column of %s" % (pid, col, hit["pkg"])
        return _r2s_canon_col(row, col)
    v, p = row.get("value", ""), row.get("p", "")
    if _r2s_isnum(v) and (not _r2s_isnum(p) or float(v) == float(p)):
        return vcol
    if _r2s_isnum(v) and _r2s_isnum(p):   # value is an estimate and p its test (ror, or_mle, diff): a p-named id would be ambiguous
        assert not _r2s_PNAME.match(hit["qid"].split(".")[-1]), "%s: value and p differ (%s vs %s) on a p-named quantity %s" % (pid, v, p, hit["qid"])
        return vcol
    if _r2s_isnum(p):
        return "p"
    return vcol   # non-numeric value (a range string such as 12.1376-12.7365) is handled by the caller




def _r2s_kind(pid, hit, col, sval):
    """Print class for the cell. Key ids use their Stage 1b kind; others the unit / name rules (name first: some p rows carry a ratio unit)."""
    row, key, qid = hit["row"], hit["key"], hit["qid"]
    u = (row.get("unit") or "").strip().lower()
    last = qid.split(".")[-1]
    if pid in _r2s_CUT3:
        return "p", 3
    if col == "n":
        return "int", None
    if key is not None:
        k = key["tolerance"]
        if k == "P":
            return "p", None
        if k == "COUNT":
            return "int", None
        if k == "PP":
            return "pp", None
        if k == "MC_SHARE":
            return "share", None
        if k == "LD":
            return "dp3", None
        if k == "EST":
            ku = (key.get("unit") or "").strip().lower()
            if ku == "or":
                return "or", None
            if ku == "z":
                return "dec2", None
            return "dp3", None
        raise AssertionError("%s: key tolerance class %s has no print rule" % (pid, k))
    if _r2s_PNAME.match(last) and col in ("p", "value"):
        return "p", None
    if qid.startswith("r2."):
        return "dp3", None
    if u.startswith("p (") or u == "p":
        return "p", None
    if "posterior probability" in u or u == "pp.h4" or last.endswith("PPH4"):
        return "pp", None
    if u == "share" or last.endswith("share") or last.endswith("share_same"):
        return "share", None
    if u in ("s", "r2", "ratio") or u.startswith("log-odds per") or u.startswith("ornithine units"):
        return "dp3", None
    if u in ("chi-square", "z"):
        return "dec2", None
    if u.startswith("or"):
        return "or", None
    if u in ("df", "pairs", "variants", "settings", "list", "count") or _r2s_re.fullmatch(r"-?\d+", sval.strip()):
        return "int", None
    return "sig3", None


# ------------------------------------------------------------------------------------------------------------------ cells
_r2s_cells = {}       # (pkg, qid, col) -> dict
_r2s_order = []
for _pid, _pk, _use, _alias in _r2s_lines:
    _hit = _r2s_resolve(_pid, _pk, _alias)
    _col = _r2s_column(_pid, _hit, _alias)
    _ck = (_hit["pkg"], _hit["qid"], _col)
    if _ck not in _r2s_cells:
        _r2s_cells[_ck] = dict(hit=_hit, col=_col, ids=[], used=set(), hits={}, alias=[])
        _r2s_order.append(_ck)
    if _alias:
        _r2s_cells[_ck]["alias"].append(_pid)
    _r2s_cells[_ck]["ids"].append(_pid)
    _r2s_cells[_ck]["hits"][_pid] = _hit
    _r2s_cells[_ck]["used"] |= {u.strip() for u in _use.split(",") if u.strip()}

for _r in _r2s_consts:   # constants: one cell each, appended after the results cells (L numbers of earlier ids never move)
    _ck = ("CONST", _r["id"], "value")
    _r2s_cells[_ck] = dict(const=_r, ids=[_r["id"]], used=set(), hits={})
    _r2s_order.append(_ck)

for _r in _r2s_derived:   # derived rows (R2_derived_rows.py): one cell each, appended after the constants
    _ck = ("DERIVED", _r["id"], "value")
    _r2s_cells[_ck] = dict(derived=_r, ids=[_r["id"]], used=set(), hits={})
    _r2s_order.append(_ck)

# Lxxx assignment: previous map first (an id keeps its number; a cell keeps the number of its first id), then new cells in file order
_r2s_cellL = {}
_used_nums = set()
for _pid, _r in _r2s_prev.items():
    _ck = tuple(_r["cell"].split("|"))
    _num = int(_r["L"][1:])
    if _ck in _r2s_cellL:
        assert _r2s_cellL[_ck] == _num, "cell %s carries two Lxxx in the previous map" % (_ck,)
    _r2s_cellL[_ck] = _num
    _used_nums.add(_num)
for _pid, _r in _r2s_prev.items():
    _ck = tuple(_r["cell"].split("|"))
    _mine = [k for k, c in _r2s_cells.items() if _pid in c["ids"]]
    assert _mine == [_ck], "id %s was re-pointed from cell %s to %s: ids are never re-pointed" % (_pid, _ck, _mine)
_r2s_new = [_ck for _ck in _r2s_order if _ck not in _r2s_cellL]
for _ck, _nxt in zip(_r2s_new, _r2s_ida.allocate(_used_nums, len(_r2s_new), _r2s_LO, _r2s_HI)):   # next free numbers, skipping the reserved ids, ceiling L998
    _r2s_cellL[_ck] = _nxt
    _used_nums.add(_nxt)
assert all(_r2s_LO <= n <= _r2s_HI for n in _r2s_cellL.values()), "ext_r2_stage1 needs Lxxx outside L%d-L%d (%d cells)" % (_r2s_LO, _r2s_HI, len(_r2s_cellL))
assert len(set(_r2s_cellL.values())) == len(_r2s_cellL)

# ------------------------------------------------------------------------------------------------------------------ rows
_r2s_sys.path.insert(0, _r2s_os.path.join(V3, "R2_rewrite", "scripts"))
import P31R_gate as _r2s_gate   # one implementation of the N2 rule (parse_text_num, matches, nums_in): the printed form is checked with it

_r2s_SECTION = {"Results": "R2 Results", "Table3": "R2 Table 3", "SM": "R2 SM results detail", "Discussion": "R2 Discussion", "Methods": "R2 Methods"}
_r2s_maprows = []
_r2s_n0 = len(ROWS)   # first row added by this module; its rows are put in ascending Lxxx order at the end (ids are numbered by cell, rows are added in file order)
_r2s_stats = dict(key=0, key_stale=0, key_visible=0, pipe=0, tie=0, fallback=0, stale_changed=0, const=0, derived=0)
_r2s_archive_cache = {}


def _r2s_archived_row(pkg, md5):
    """The results.tsv copy (under the package folder, _pre_* snapshots) whose md5 is the one a key was verified against."""
    if (pkg, md5) not in _r2s_archive_cache:
        found = None
        base = _r2s_os.path.join(_r2s_R1, pkg)
        for d in sorted(_r2s_os.listdir(base)):
            p = _r2s_os.path.join(base, d, "results.tsv")
            if d.startswith("_pre_") and _r2s_os.path.isfile(p) and _r2s_md5(p) == md5:
                found = p
                break
        _r2s_archive_cache[(pkg, md5)] = found
    return _r2s_archive_cache[(pkg, md5)]


# ---- 2026-10-14 (infrastructure round, after the last-round skeptic): a KEY row whose Stage 1b status was scored against a SUPERSEDED value
# is no longer labelled plain 'KEY'. The status is carried into the label and the current value is said to be PIPE, with the cross-checks
# that exist at the current value (verifier lenses L1/L2 in S6_SURVIVAL.tsv; A_third where it computed the quantity), all read from files.
# Relabelled when: the mapped results row CHANGED since the key was verified (archived copy compared), or the label is L1L2_AGREE_IMPL_DIFFERS,
# or the label is SURVIVES_VISIBLE_TARGET. Every other key keeps its previous text byte for byte.
_r2s_XTOL = 1e-6   # tolerance of the L1/L2/A_third cross-checks: STEP2_7FIT_REPORT.md section 5, "1e-6 on every comparison (W9 assert)"
_r2s_S6_SURV = _r2s_os.path.join(_r2s_R1, "S1R1_R2_1b_verify", "S6", "adjudication", "S6_SURVIVAL.tsv")
_r2s_A3_DIR = _r2s_os.path.join(_r2s_R1, "S1R1_R2_1b_verify", "S6", "adjudication", "A_third")
_r2s_xc_cache = {}


def _r2s_xc_tsv(path):
    if path not in _r2s_xc_cache:
        with _r2s_io.open(path, encoding="utf-8", newline="") as fh:
            _r2s_xc_cache[path] = list(_r2s_csv.DictReader(fh, delimiter="\t"))
    return _r2s_xc_cache[path]


def _r2s_g(x):
    return "%.10g" % x


def _r2s_cross_check(kid, cell):
    """-> (text, impl value or None, visible_in). Cross-checks of the CURRENT value of a Stage 1b S6 key by L1, L2 (S6_SURVIVAL.tsv) and A_third."""
    k = _r2s_KEYS[kid]
    cur = float(cell["hit"]["row"][cell["col"]])
    sv = [r for r in _r2s_xc_tsv(_r2s_S6_SURV) if r["key_id"] == kid]
    if len(sv) != 1:
        return "no L1/L2 record for this key in S6_SURVIVAL.tsv", None, ""
    sv = sv[0]
    out, bad = [], False
    for lens in ("L1", "L2"):
        v = float(sv[lens])
        d = abs(v - cur)
        bad = bad or d > _r2s_XTOL
        out.append("%s %s (|diff| %.1e)" % (lens, _r2s_g(v), d))
    a3 = None
    qid = cell["hit"]["qid"]
    mm = _r2s_re.fullmatch(r"S6\.d\.(\w+)\.s_max", qid)
    if mm:
        rows = [r for r in _r2s_xc_tsv(_r2s_os.path.join(_r2s_A3_DIR, "a3_allcand_fitset.tsv")) if r["cid"] == mm.group(1)]
        vals = [float(r[c]) for r in rows for c in ("sP", "sO") if _r2s_isnum(r[c])]
        if vals:
            a3 = max(vals)
    elif qid in ("S6.d.lo", "S6.d.reclassified"):
        rows = [r for r in _r2s_xc_tsv(_r2s_os.path.join(_r2s_A3_DIR, "a3_d_table.tsv"))
                if (r["fit_set"], r["shift_rows"], r["version"]) == ("eight", "keep", "reclassify_ING1")]
        if len(rows) == 1:
            a3 = float(rows[0]["lo" if qid == "S6.d.lo" else "n_reclassified"])
    if a3 is not None:
        d = abs(a3 - cur)
        bad = bad or d > _r2s_XTOL
        out.append("A_third (own .bed decoder) %s (|diff| %.1e)" % (_r2s_g(a3), d))
    else:
        out.append("A_third did not compute this quantity")
    txt = "cross-checked at the current value by verifier lenses " + ", ".join(out) + (
        "" if not bad else " -- OUTSIDE the 1e-6 tolerance; the current value is NOT cross-checked") + (
        "; these are checks against the verifiers' and adjudicator's own pipeline, not independent proof about the fit-set rule" if a3 is not None else "")
    imp = sv.get("impl", "")
    return txt, (float(imp) if _r2s_isnum(imp) else None), sv.get("visible_in", "")


def _r2s_key_status(cell):
    parts, flagged = [], []
    for kid in cell["keys"]:
        k, lab = _r2s_KEYS[kid], _r2s_LAB[kid]
        s = "KEY %s: %s" % (kid, lab["status_corrected"])
        if lab["qualifier"]:
            s += " [%s]" % lab["qualifier"].strip().rstrip(";")
        cur = cell["hit"]["md5"]
        changed, oldval = False, None
        if k["map_md5"] != cur and k["map_file"].endswith("results.tsv"):
            oldp = _r2s_archived_row(cell["hit"]["pkg"], k["map_md5"])
            if oldp is None:
                s += " [verified against results.tsv md5 %s; current %s; verified copy not found, row not compared]" % (k["map_md5"][:8], cur[:8])
            else:
                _, orows = _r2s_read_tsv(oldp)
                orow = [r for _, r in orows if r[cell["hit"]["hdr"][0]] == cell["hit"]["qid"]]
                same = len(orow) == 1 and all(orow[0].get(c) == cell["hit"]["row"].get(c) for c in ("value", "se", "ci_lo", "ci_hi", "p", "n") if c in orow[0])
                s += " [verified against results.tsv md5 %s (%s); current %s; the mapped row is %s]" % (
                    k["map_md5"][:8], _r2s_os.path.relpath(oldp, _r2s_R1).replace("\\", "/"), cur[:8], "IDENTICAL in value/se/ci/p/n" if same else "CHANGED since that verification")
                if not same:
                    _r2s_stats["stale_changed"] += 1
                    changed = True
                    oldval = orow[0].get(cell["col"]) if len(orow) == 1 else None
        st = lab["status_corrected"]
        if changed or st.startswith("L1L2_AGREE_IMPL_DIFFERS") or st.startswith("SURVIVES_VISIBLE_TARGET"):
            xc, imp, vis = _r2s_cross_check(kid, cell) if k["_sec"] == "S6" else ("", None, "")
            if changed or st.startswith("L1L2_AGREE_IMPL_DIFFERS"):
                cert = oldval if oldval is not None and _r2s_isnum(oldval) else (_r2s_g(imp) if imp is not None else "unrecorded")
                if oldval is not None and imp is not None and _r2s_isnum(oldval):
                    assert abs(float(oldval) - imp) < 1e-9, "%s: archived value %s and S6_SURVIVAL impl %s disagree" % (kid, oldval, imp)
                what = "four-fit value" if cell["hit"]["qid"].endswith(".s_max") else "four-fit-derived value"
                verb = ("certified the superseded %s" % what) if st.startswith("SURVIVES") else ("was scored against the superseded %s" % what)
                s = ("KEY-STALE %s: Stage 1b status %s %s %s (results.tsv md5 %s, %s), so it does not speak to the current seven-fit value %s (md5 %s), "
                     "which is PIPE: analyst-computed, not a Stage 1b verdict; %s [label and qualifier: stage1b_key_labels.tsv]" % (
                         kid, st.split(" (")[0], verb, cert, k["map_md5"][:8], "archived pre-step-2 row" if oldval is not None else "S6_SURVIVAL.tsv impl column",
                         cell["hit"]["row"][cell["col"]], cur[:8], xc or "no cross-check recorded for this key's section"))
                flagged.append((0, "stale"))
            else:   # SURVIVES_VISIBLE_TARGET with an unchanged row
                s = ("KEY-VISIBLE %s: Stage 1b status %s, so the verifiers' agreement is not a blind re-derivation (S6_SURVIVAL.tsv visible_in: %s); "
                     "the mapped row is unchanged since verification (md5 %s); %s [label: stage1b_key_labels.tsv]" % (
                         kid, st, vis or "not recorded", cur[:8], xc or "no cross-check recorded for this key's section"))
                flagged.append((1, "visible"))
        else:
            flagged.append((2, "clean"))
        parts.append(s)
    order = sorted(range(len(parts)), key=lambda i: flagged[i][0])   # a relabelled part leads, so the first token of the check names the weakest status
    _r2s_stats["relabelled"] = _r2s_stats.get("relabelled", []) + [(cell["keys"][i], flagged[i][1]) for i in order if flagged[i][1] != "clean"]
    return " | ".join(parts[i] for i in order)


_r2s_NUMTOK = _r2s_re.compile(r"-?\d[\d,]*(?:\.\d+)?%?")


def _r2s_const_row(ck, r):
    """One ledger row + one map row for a protocol constant / engine setting. Everything is re-read from the source at this build."""
    L_ = "L%d" % _r2s_cellL[ck]
    doc, cid = r["source document"], r["id"]
    path = _r2s_os.path.join(V3, doc)
    assert _r2s_os.path.isfile(path), "%s: source document %s does not exist" % (cid, doc)
    md5_ = _r2s_md5(path)
    assert md5_ == r["doc_md5"], "%s: source %s changed (md5 %s, pinned %s): re-read the quote and regenerate R2_protocol_constants.tsv" % (cid, doc, md5_, r["doc_md5"])
    with _r2s_io.open(path, encoding="utf-8", newline="") as fh_:
        lines_ = [x.rstrip("\r") for x in fh_.read().split("\n")]
    ln_ = int(r["line"])
    assert 1 <= ln_ <= len(lines_) and r["verbatim quote"] in lines_[ln_ - 1], "%s: the verbatim quote is not on line %d of %s" % (cid, ln_, doc)
    assert r["token"] in r["verbatim quote"], "%s: token %r is not inside its quote" % (cid, r["token"])
    val_ = _r2s_tok(r["token"])
    assert (r["value"] == val_) if isinstance(val_, str) else (float(r["value"]) == float(val_)), "%s: value column %r disagrees with the token %r (parsed %r)" % (cid, r["value"], r["token"], val_)
    disp_ = r["display"]
    prose_ = _r2s_re.sub(r"10\^(-?\d+)", lambda m_: "10" + m_.group(1).translate(_r2s_SUP), disp_)
    if not isinstance(val_, str):   # round-trip with the gate's own N2 functions, on every simple numeric form the prose can print
        for tk_ in (r["token"], prose_):
            if _r2s_NUMTOK.fullmatch(tk_) or _r2s_re.fullmatch(r"10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+", tk_):
                tv_, tol_, sci_ = _r2s_gate.parse_text_num(tk_.rstrip("%"))
                m_ = _r2s_gate.matches(tv_, tol_, sci_, _r2s_gate.nums_in(dict(value=json.dumps(val_), display=disp_)), tk_.endswith("%"))
                assert m_ is True, "%s: printed %r does not round-trip to the stored value (%r)" % (cid, tk_, m_)
    add(L_, "R2 protocol constant (governing document)" if r["class"] == "protocol" else "R2 engine setting (deposited script)", r["meaning"], val_, disp_,
        "%s; line %d; md5 %s" % (doc, ln_, md5_), "constant=%s; token=%s" % (cid, r["token"]), "[R]", "",
        "Placeholder id(s): %s. Verbatim source line %d: \"%s\". The value is parsed from the token as written in the source; the document md5, the line and the quote are asserted at every build." % (cid, ln_, r["verbatim quote"]))
    EXT_CHECKS[L_] = ("SRC-QUOTE: read from %s line %d; the verbatim quote, the line and the document md5 are asserted at every build; one reader, "
                      "not independently re-derived (a second reader should re-read that line)" % (_r2s_os.path.basename(doc), ln_))
    _r2s_stats["const"] += 1
    _r2s_maprows.append(dict(id=cid, L=L_, printed=prose_, ledger_display=disp_, source_value=r["token"], package="CONST:" + r["class"], quantity_id=cid, column="value",
                             line=str(ln_), results_md5=md5_, check="SRC-QUOTE", print_rule="constant", cell="|".join(ck), shared_with=""))


def _r2s_derived_row(ck, r):
    """One ledger row + one map row for a derived entry of R2_derived_rows.py. The gate's own D1 check runs on the new row at build time:
    a derived row whose display prints a number its value does not store stops the build."""
    L_ = "L%d" % _r2s_cellL[ck]
    add(L_, r["section"], r["claim"], r["value"], r["display"], r["source"], r["key"], "[R]", "", r["note"])
    d1_ = _r2s_gate.d1_row(ROWS[-1])
    assert not d1_, "%s (%s): the new row fails D1: %s" % (L_, r["id"], d1_)
    EXT_CHECKS[L_] = r["check"]
    _r2s_stats["derived"] += 1
    _r2s_maprows.append(dict(id=r["id"], L=L_, printed=r["display"], ledger_display=r["display"], source_value=ROWS[-1]["value"][:200], package="DERIVED", quantity_id=r["id"],
                             column="value", line="", results_md5="", check="DERIVED", print_rule="derived", cell="|".join(ck), shared_with=""))


for _ck in _r2s_order:
    _c = _r2s_cells[_ck]
    if "const" in _c:
        _r2s_const_row(_ck, _c["const"])
        continue
    if "derived" in _c:
        _r2s_derived_row(_ck, _c["derived"])
        continue
    _hit, _col = _c["hit"], _c["col"]
    _c["keys"] = _r2s_KEY_BY_CELL.get(_ck, [])
    _c["keys"] = [k for k in _c["keys"]]
    _sval = _hit["row"][_col].strip()
    _pid0 = _c["ids"][0]
    _rng = _r2s_re.fullmatch(r"(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)", _sval)
    if _rng:   # a range string such as 12.1376-12.7365 (seed range of a bootstrap bound)
        _a, _b = _rng.group(1), _rng.group(2)
        _kind, _sf = "or", None
        _pa, _pr1, _t1 = _r2s_print("or", _a)
        _pb, _pr2, _t2 = _r2s_print("or", _b)
        _ascii, _prose, _tie = "%s to %s" % (_pa, _pb), "%s to %s" % (_pr1, _pr2), _t1 or _t2
        _val = [_r2s_num(_a), _r2s_num(_b)]
        _kind_label = "or (range)"
    else:
        assert _r2s_isnum(_sval), "%s: cell %s holds a non-numeric value %r" % (_pid0, _ck, _sval)
        _others = {_r2s_kind(i, _c["hits"][i], _col, _sval) for i in _c["ids"]}   # ids sharing a cell must agree on the print class
        if len(_others) > 1:   # a direct id with no class rule of its own (default sig3) yields to the rule of the other id (the Stage 1b key's)
            _others = {k for k in _others if k[0] != "sig3"}
        assert len(_others) == 1, "ids %s share cell %s but disagree on the print rule: %s" % (_c["ids"], _ck, _others)
        _kind, _sf = next(iter(_others))
        _ascii, _prose, _tie = _r2s_print(_kind, _sval, _sf)
        _val = _r2s_num(_sval)
        _kind_label = _kind + ("" if _sf is None else "/%d sf" % _sf)
        if _kind == "sig3":
            _r2s_stats["fallback"] += 1
    if _tie:
        _r2s_stats["tie"] += 1
    _L = "L%d" % _r2s_cellL[_ck]
    _row = _hit["row"]
    _est = (_row.get("estimand") or "").strip()
    if not _est:   # package 12 has no estimand column: the description before the first ';' in notes
        _est = (_row.get("notes") or "").split(";")[0].strip()
    _kest = _r2s_KEYS[_c["keys"][0]]["estimand"].strip() if _c["keys"] else ""
    _kdiff = ""
    if _c["keys"] and _col != "value":   # the row's estimand describes its value column; a key on p or n is described by its own definition
        _est = _kest
    elif _c["keys"] and _kest != _est:   # the row text is current; the key text is what was verified
        _kdiff = " Stage 1b key definition at verification (%s): %s" % (_c["keys"][0], _kest)
    _ctx = [x for x in (_row.get("resource"), ("panel " + _row["panel"]) if _row.get("panel") not in (None, "", "NA") else "",
                        _row.get("cov"), (_row.get("se_order") + "-order SE") if _row.get("se_order") not in (None, "", "NA") else "") if x and x not in _r2s_NA]
    _claim = _est + (" [" + "; ".join(_ctx) + "]" if _ctx else "")
    if _c["alias"]:
        _claim = "count (column n) of: " + _claim
    _used = sorted(_c["used"], key=["Results", "Table3", "SM", "Discussion", "Methods"].index) if _c["used"] <= set(_r2s_SECTION) else sorted(_c["used"])
    _sect = "R2 Stage 1 package " + _hit["pkg"].replace("S1R1_R2_", "") + (" (used in " + ", ".join(_used) + ")" if _used else "")
    _src = "%s; line %d; column %s; md5 %s" % (rel(_hit["path"]), _hit["line"], _col, _hit["md5"])
    _note = "Print rule %s; printed %s. Placeholder id(s): %s." % (_kind_label, _prose, ", ".join(_c["ids"]))
    if _col == "p" and _r2s_isnum(_row.get("value", "")):
        _note += " The row's value column holds an estimate (%s); the printed number is the p column." % _row["value"]
    if _c["alias"]:
        _note += " ALIAS: %s names this count under another name: it is column %s of quantity %s (the row's own value column holds %s)." % (
            ", ".join(_c["alias"]), _col, _hit["qid"], (_row.get("value") or "")[:60])
    elif _col == "n" and _row.get("value"):
        _note += " The row's value column holds the list: %s." % _row["value"]
    _src_note = (_row.get("notes") or "").strip()
    if _src_note:
        _note += " Source notes: %s" % _src_note[:400]
    _note += _kdiff
    if _kind == "sig3":
        _note += " DEFAULT print rule (3 significant figures): no class rule exists for unit %r; check the precision before printing." % _row.get("unit", "")
    if _rng:
        _note += " Range row: print each end with its own tag (12.1 [[Lxxx]] to 12.7 [[Lxxx]]); the gate does not look ahead across a decimal point."
    if _tie:
        _note += " WARNING: the source string sits on a rounding tie at this precision."
    _other = [k for (p_, q_, c_), kk in _r2s_KEY_BY_CELL.items() if (p_, q_) == (_hit["pkg"], _hit["qid"]) and c_ != _col for k in kk]
    if _other and not _c["keys"]:
        _note += " The same results row is Stage 1b key %s on column %s; this cell is a different column of it." % (", ".join(_other), _r2s_KEYS[_other[0]]["compare_col"])
    add(_L, _sect, _claim, _val, _ascii, _src, "quantity_id=%s; column=%s" % (_hit["qid"], _col), "[R]", "Y" if _c["keys"] else "", _note)
    # N2 self-check with the gate's own functions: the printed form must round-trip to the stored full-precision value
    _chk = dict(value=ROWS[-1]["value"], display=ROWS[-1]["display"])
    for _tok in ((_ascii.split(" to ") if " to " in _ascii else [_ascii])):
        _tv, _tol, _sci = _r2s_gate.parse_text_num(_tok.replace("10^", "10").replace("−", "-"))
        _m = _r2s_gate.matches(_tv, _tol, _sci, _r2s_gate.nums_in(_chk), False)
        assert _m is True, "%s %s: printed %r does not round-trip to the stored value (%r)" % (_L, _c["ids"], _tok, _m)
    # independent check
    if _c["keys"]:
        EXT_CHECKS[_L] = _r2s_key_status(_c)
        _r2s_stats["key_stale" if EXT_CHECKS[_L].startswith("KEY-STALE") else "key_visible" if EXT_CHECKS[_L].startswith("KEY-VISIBLE") else "key"] += 1
    elif _hit["pkg"] == "S1R1_R2_12_chen_instrument" and _hit["qid"] in ("r2.C_A.503", "r2.C_B.503"):
        # second reader: Python with its own .bed decoder (S1R1_R2_12_chen_instrument_verify/mine.tsv); recorded in RECONCILE.md section 2a
        _vp = _r2s_os.path.join(_r2s_R1, "S1R1_R2_12_chen_instrument_verify", "mine.tsv")
        _lab = {"r2.C_A.503": "rs2781668 A", "r2.C_B.503": "rs17788484 B"}[_hit["qid"]]
        _mine = [float(l.split("\t")[3]) for l in _r2s_io.open(_vp, encoding="utf-8").read().split("\n")
                 if l.startswith("LD\t") and ("rs2246012, " in l) and (_lab in l) and "[r2]" in l and l.split("\t")[2] == "503"]
        assert len(_mine) == 1, "second-reader row for %s not unique in mine.tsv" % _hit["qid"]
        _diff = abs(_mine[0] - float(_sval))
        assert _diff <= 0.005, "second reader disagrees on %s by %g (LD tolerance 0.005)" % (_hit["qid"], _diff)
        EXT_CHECKS[_L] = ("MATCH: independent second reader (Python, own .bed decoder, S1R1_R2_12_chen_instrument_verify/mine.tsv) differs by %.1e (LD tolerance 0.005; "
                          "RECONCILE.md section 2a); the checker was not blind to the AMEND2 v3 record of r2 0.993; not a Stage 1b KEY" % _diff)
        _r2s_stats["pipe"] += 0
    else:
        EXT_CHECKS[_L] = "PIPE: analyst-computed, not independently re-derived"
        _r2s_stats["pipe"] += 1
    for _pid in _c["ids"]:
        _r2s_maprows.append(dict(id=_pid, L=_L, printed=_prose, ledger_display=_ascii, source_value=_sval, package=_hit["pkg"], quantity_id=_hit["qid"],
                                 column=_col, line=str(_hit["line"]), results_md5=_hit["md5"], check=EXT_CHECKS[_L].split(":")[0].split(" ")[0] + (" " + ",".join(_c["keys"]) if _c["keys"] else ""),
                                 print_rule=_kind_label, cell="|".join(_ck), shared_with=",".join(i for i in _c["ids"] if i != _pid)))

ROWS[_r2s_n0:] = _r2s_ida.sorted_by_id(ROWS[_r2s_n0:])   # loud on an id that is not L + 3 digits (int(id[1:4]) would read L1000 as 100)

# map file: written only when every row above was built (an exception leaves the previous map untouched)
_r2s_fields = ["id", "L", "printed", "ledger_display", "source_value", "package", "quantity_id", "column", "line", "results_md5", "check", "print_rule", "shared_with", "cell"]
_r2s_tmp = _r2s_MAP + ".%d.tmp" % _r2s_os.getpid()
with _r2s_io.open(_r2s_tmp, "w", encoding="utf-8", newline="") as _fh:
    _w = _r2s_csv.DictWriter(_fh, fieldnames=_r2s_fields, delimiter="\t", lineterminator="\n")
    _w.writeheader()
    _w.writerows(sorted(_r2s_maprows, key=lambda r: (int(r["L"][1:]), r["id"])))
_r2s_os.replace(_r2s_tmp, _r2s_MAP)
print("ext_r2_stage1: %d ids -> %d cells (L%d-L%d); KEY cells %d, KEY-STALE cells %d, KEY-VISIBLE cells %d, PIPE cells %d, stale-key rows changed %d, rounding ties %d, default print rule %d -> %s" % (
    len(_r2s_idl), len(_r2s_cells), min(_r2s_cellL.values()), max(_r2s_cellL.values()), _r2s_stats["key"], _r2s_stats["key_stale"], _r2s_stats["key_visible"],
    _r2s_stats["pipe"], _r2s_stats["stale_changed"], _r2s_stats["tie"], _r2s_stats["fallback"], _r2s_MAP))
print("ext_r2_stage1: protocol constants / engine settings %d (R2_protocol_constants.tsv); alias ids %d" % (_r2s_stats["const"], sum(len(c["alias"]) for c in _r2s_cells.values() if "alias" in c)))
print("ext_r2_stage1: derived rows %d (R2_derived_rows.py)" % _r2s_stats["derived"])
print("ext_r2_stage1: KEY labels relabelled (not plain KEY): %s" % ", ".join("%s=%s" % x for x in _r2s_stats.get("relabelled", [])))
