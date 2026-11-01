"""rel_portability.py -- DECLARED TRANSFORM of the SHIPPED COPIES of the code bundle (release candidate RC1, brief 'Agent C', item 2).

The live generators and scripts are never edited (an edit forces a re-render and moves pinned md5s). The copies under out/code/ are
rewritten here, from the source file every time (idempotent), so that the code runs from a clean extraction:

  * an absolute path in an EXECUTABLE line (a string literal in code) becomes an archive-relative expression resolved from ONE
    configuration:   environment P31R_ARCHIVE  (default: the archive root, derived from the script's own place, <archive>/code/<path>)
                     environment P31R_EXTERNAL (default: <archive>/external_inputs; third-party inputs that are not in the release)
        project V3 root      -> <archive>/code     (the project tree is mirrored under code/)
        other project folder -> <external>/CKM papers/...
        drive path (D:/ E:/) -> <external>/D/... , <external>/E/...
        Rscript / python / the user's tool folder / the Windows font folder -> found on PATH, or an environment override;
  * an absolute path in a COMMENT, docstring or data cell becomes a redacted relative form (<archive>/code/..., <CKM>/..., <D>/...);
  * P31R_fig_common.R (shipped copy) reads the claims ledger from the PUBLIC extract (ledger/claims_ledger_public.tsv), resolves its
    inputs through code/P31R_INPUT_MAP.tsv (a verbatim released copy is used; anything else stops with 'input not in the release: ...'),
    and stops with the reason when it is asked for a masked ledger value;
  * the same transform covers EVERY other shipped executable (*.R, *.py, *.sh outside code/) and every shipped text table or document (RC1 fix K01/K02): an executable gets
    the expressions and prologue, a data cell the public relative form of the claims ledger (project-relative, or <external>/...);
  * every file's source md5, shipped md5 and applied rules go to build/portability_log.tsv;
  * RC1b (S2-D34(3)), level-2 lift: inputs/LEVEL2_INPUTS.tsv declares, per source a figure generator reads, the released file that stands in for it in the archive
    (kind released_copy: a masked or reduced copy another step already ships) or a reduced extract this module writes under code/lift_inputs/ (kind extract:
    named columns kept or dropped, or a named section of a text). Each such map row carries lift = 1 and the units it serves; the shipped resolver then accepts
    the copy although it is not verbatim, and the figure-source check compares the SOURCE md5 recorded in the map (pinned in RELEASE_INPUTS, verified at build)
    with the approved manifest. Whether a reduced copy is enough is proved by the clean-room regeneration (plotted frames md5-equal to the live frames), not here.

Nothing is typed: every replaced path is matched by pattern in the file itself and verified by a round trip (the expression, with each
root put back, reproduces the original literal). The module is shipped too, so the patterns are built from pieces.
"""
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import ast
import csv
import hashlib
import io
import re
import tempfile
import tokenize
from pathlib import Path

CODE_REL = "R2_rewrite/release/inputs/CODE_INCLUDE.tsv"
RI_REL = "R2_rewrite/release/inputs/RELEASE_INPUTS.tsv"
MAP_REL = "code/P31R_INPUT_MAP.tsv"
LOG_REL = "code/P31R_PORTABILITY_LOG.tsv"
MAP_COLS = ["source_path", "release_path", "md5", "source_md5", "verbatim", "lift", "lift_units"]
L2_REL = "R2_rewrite/release/inputs/LEVEL2_INPUTS.tsv"
L2_COLS = ["l2_id", "units", "kind", "source_path", "release_path", "spec", "needs", "reason"]
L2_KINDS = ("released_copy", "extract")
L2_PREFIX = "code/lift_inputs/"          # the only place this module writes an extract (MODULES lets rel_portability write under code/ only)
_REL_DIR = Path(__file__).resolve().parent         # the release folder this module was imported from (the builder's --release)
_L2_FILE = _REL_DIR / "inputs" / "LEVEL2_INPUTS.tsv"
_RI_FILE = _REL_DIR / "inputs" / "RELEASE_INPUTS.tsv"
LOG_COLS = ["path", "source_md5", "shipped_md5", "rules"]
LEDGER_PUBLIC_REL = "ledger/claims_ledger_public.tsv"


class PortabilityError(AssertionError):
    pass


def _req(cond, msg):
    if not cond:
        raise PortabilityError(msg)


def md5_bytes(b):
    return hashlib.md5(b).hexdigest()


# ---------------------------------------------------------------------------------------------------------- the patterns
# The home path (and the drive roots below) are written as concatenated fragments on purpose: this module is itself shipped, and the release gate's absolute-path scans (R12, R13)
# read every shipped file, including this one, so the tool's own constants must not match the pattern it enforces.
_HOME = Path.home().as_posix()                              # the running user's home folder (no personal path is written in this shipped file)
_CKM = _HOME + "/Downloads/CKM papers"
_PY_EXE = _HOME + "/AppData/Local/Programs/Python/Python312/python.exe"
_RS_EXE = "C" + ":/Program Files/R/R-4.4.1/bin/Rscript"
_FONTS = "C" + ":/Windows/Fonts"
_TOOLS = (_P31R_EXT + "/tools/")
# S6.3b2a (D4): an Rscript path of any installed version and either separator (the v5 script 61 names another version than the figure code does)
_RS_RX = "C" + ":" + r"[/\\]Program Files[/\\]R[/\\]R-[0-9][0-9.]*[/\\]bin[/\\]Rscript(?:\.exe)?"
# one alternation, longest first; group names say what the match is
def _anysep(s):
    """a Windows path may be written with / or with a backslash (raw strings)"""
    return re.escape(s).replace("/", r"[/\\]")


PATH_RX = re.compile(
    r"(?P<scratch>" + re.escape(_HOME + "/AppData/Local/Temp/claude/") + r"[^\s\"']*)"
    r"|(?P<pyexe>" + re.escape(_PY_EXE) + r")"
    r"|(?P<rscript>" + _RS_RX + r")"
    r"|(?P<fonts>" + re.escape(_FONTS) + r")"
    r"|(?P<v3>" + _anysep(_CKM) + r"[/\\]V3(?![A-Za-z0-9_.]))"
    r"|(?P<ckm>" + _anysep(_CKM) + r"(?![A-Za-z0-9_.]))"
    r"|(?P<home>" + _anysep(_HOME) + r"(?![A-Za-z0-9_.]))"
    r"|(?P<tools>" + re.escape(_TOOLS) + r")"
    r"|(?<![A-Za-z0-9_])(?P<drive>[CDE]):[/\\](?=[A-Za-z0-9_$.~])")      # S6.3b2a (D4): a backslash root too (a Python raw string); the text after it is rewritten with '/'
LEAK_RX = re.compile(r"(?<![A-Za-z0-9_])[A-Za-z]:(?:/|\\{1,2})(?=[A-Za-z0-9_$.~])(?!\\?[ntr](?![A-Za-z0-9_$.~/\\]))|" + re.escape(_HOME) + r"|Users(?:/|\\{1,2})" + re.escape(Path.home().name))


def leaks(text):
    """(line number, line) of every line that still names an absolute local path; 'https://' and friends are not drive paths (one letter, upper case)."""
    out = []
    for i, ln in enumerate(text.split("\n"), 1):
        for m in LEAK_RX.finditer(ln):
            if m.group(0)[:1].isupper() or m.group(0).startswith("Users"):
                out.append((i, ln.strip()[:140]))
                break
    return out


# redaction (comments, docstrings, data cells)
def redact(text):
    def sub(m):
        k = m.lastgroup
        if k in ("pyexe",):
            return "python"
        if k == "rscript":
            return "Rscript"
        if k == "scratch":
            return "<scratch>"
        if k == "fonts":
            return "<system-fonts>"
        if k == "v3":
            return "<archive>/code"
        if k == "ckm":
            return "<CKM>"
        if k == "home":
            return "<home>"
        if k == "tools":
            return "<tools>/"
        if k == "drive":
            return "<%s>/" % m.group("drive")
        return m.group(0)
    return PATH_RX.sub(sub, text)


def redact_cell(text):
    """Data cells (a manifest, a registry): a path that names an input keeps a RESOLVABLE relative form (from <archive>/code)."""
    def sub(m):
        k = m.lastgroup
        if k == "ckm":
            return "../external_inputs/CKM papers"
        if k == "drive":
            return "../external_inputs/%s/" % m.group("drive")
        return redact(m.group(0))
    text = re.sub(re.escape(_CKM + "/V3") + r"(?![A-Za-z0-9_.])/?", "", text)      # V3-relative: the path as the manifest keys it
    return PATH_RX.sub(sub, text)


# ---------------------------------------------------------------------------------------------------------- expressions
R_EXPR = {"v3": ".P31R_V3", "ckm": ".P31R_CKM", "ext": ".P31R_EXT", "rscript": 'file.path(R.home("bin"), "Rscript")', "pyexe": '"python"',
          "fonts": '"."', "home": ".P31R_EXT", "tools": None}
PY_EXPR = {"v3": "_P31R_V3", "ckm": "_P31R_CKM", "ext": "_P31R_EXT", "rscript": "_P31R_RSCRIPT", "pyexe": "_p31r_sys.executable",
           "fonts": "_P31R_FONTDIR", "home": "_P31R_EXT", "tools": None}


def split_body(body, lang):
    """Split the text of one string literal into [(kind, text)]: kind 'lit' (kept) or an expression key; also returns the originals for the round trip."""
    parts, pos, bs = [], 0, False        # bs: the last drive root was written with a backslash; the path text after it is rewritten with '/' (the original text is kept for the round trip)
    for m in PATH_RX.finditer(body):
        k = m.lastgroup
        if m.start() > pos:
            seg = body[pos:m.start()]
            parts.append(("lit", seg.replace("\\", "/") if bs else seg, seg))
        bs = False
        if k == "drive":
            bs = m.group(0)[2:3] == "\\"
            parts.append(("ext", "", m.group(0)))
            parts.append(("lit", "/%s/" % m.group("drive"), ""))
        elif k == "tools":
            parts.append(("ext", "", m.group(0)))
            parts.append(("lit", "/tools/", ""))
        elif k == "scratch":
            parts.append(("ext", "", m.group(0)))
            parts.append(("lit", "/scratch", ""))
        else:
            parts.append((k, "", m.group(0)))
        pos = m.end()
    if pos < len(body):
        seg = body[pos:]
        parts.append(("lit", seg.replace("\\", "/") if bs else seg, seg))
    merged = []
    for k, tx, o in parts:                  # adjacent literals are one literal
        if merged and k == "lit" and merged[-1][0] == "lit":
            merged[-1] = ("lit", merged[-1][1] + tx, merged[-1][2] + o)
        else:
            merged.append((k, tx, o))
    return merged


def round_trip(parts, body):
    """Each expression put back as the path it replaced must reproduce the original literal."""
    # a synthetic literal (the '/D/' after a drive root) has orig == "" and is part of the match before it; real literals have orig == txt
    return "".join(orig for _, _, orig in parts) == body


def make_expr(parts, q, lang, fstring=False, qclose=None):
    exprs = R_EXPR if lang == "R" else PY_EXPR
    if fstring:
        return "".join(("{%s}" % exprs[k] if k != "lit" else t) for k, t, _ in parts)
    items = []
    for k, t, _ in parts:
        if k == "lit":
            if t:
                items.append(q + t + (qclose or q))
        else:
            _req(exprs[k] is not None, "no expression for %s" % k)
            items.append(exprs[k])
    if len(items) == 1:
        return items[0]
    return ("paste0(%s)" % ", ".join(items)) if lang == "R" else ("(%s)" % " + ".join(items))


# ---------------------------------------------------------------------------------------------------------- R
def r_segments(t):
    """[(kind, start, end)] with kind code | str | comment. Raw strings are refused (none exist in the bundle)."""
    n, i, segs, cs = len(t), 0, [], 0
    while i < n:
        c = t[i]
        if c == "#":
            if cs < i:
                segs.append(("code", cs, i))
            j = t.find("\n", i)
            j = n if j < 0 else j
            segs.append(("comment", i, j))
            i = cs = j
        elif c in "\"'":
            _req(not (i > 0 and t[i - 1] in "rR" and i + 1 < n and t[i + 1] in "([{-"), "R raw string at offset %d" % i)
            if cs < i:
                segs.append(("code", cs, i))
            j = i + 1
            while j < n and t[j] != c:
                j += 2 if t[j] == "\\" else 1
            _req(j < n, "unterminated R string at offset %d" % i)
            segs.append(("str", i, j + 1))
            i = cs = j + 1
        elif c == "`":
            j = t.find("`", i + 1)
            i = n if j < 0 else j + 1
        else:
            i += 1
    if cs < n:
        segs.append(("code", cs, n))
    return segs


def _kinds(text):
    k = ["c"] * len(text)
    for kind, a, b in r_segments(text):
        if kind != "code":
            for i in range(a, b):
                k[i] = "s"
    return k


def _enclosing_open(text, kinds, pos):
    depth, i = 0, pos - 1
    while i >= 0:
        if kinds[i] == "c":
            if text[i] in ")]}":
                depth += 1
            elif text[i] in "([{":
                if depth == 0:
                    return i
                depth -= 1
        i -= 1
    return -1


def _matching_close(text, kinds, o):
    depth = 0
    for i in range(o, len(text)):
        if kinds[i] == "c":
            if text[i] in "([{":
                depth += 1
            elif text[i] in ")]}":
                depth -= 1
                if depth == 0:
                    return i
    return -1


def _split_top(text, kinds, a, b):
    parts, depth, s = [], 0, a
    for i in range(a, b):
        if kinds[i] == "c":
            if text[i] in "([{":
                depth += 1
            elif text[i] in ")]}":
                depth -= 1
            elif text[i] == "," and depth == 0:
                parts.append(text[s:i]); s = i + 1
    parts.append(text[s:b])
    return parts


def rewrite_named_calls(text):
    """R cannot take an expression as an argument NAME. A call c("<path>" = v, ...) / list(...) whose names are path literals becomes
    setNames(c(v, ...), c("<path>", ...)), the same named vector; the literals are then rewritten by the ordinary pass."""
    n = 0
    while True:
        hit = None
        for kind, a, b in r_segments(text):
            if kind == "str" and PATH_RX.search(text[a:b]) and re.match(r"\s*=(?!=)", text[b:b + 12]):
                hit = a
                break
        if hit is None:
            return text, n
        kinds = _kinds(text)
        o = _enclosing_open(text, kinds, hit)
        _req(o > 0, "a path-named argument outside a call")
        m = re.search(r"(?<![A-Za-z0-9_.])(c|list)\s*$", text[:o])
        _req(m, "a path-named argument in a call that is not c() or list(): %s" % text[max(0, o - 20):o + 20])
        c = _matching_close(text, kinds, o)
        args = _split_top(text, kinds, o + 1, c)
        names, vals = [], []
        for ar in args:
            mm = re.match(r"^\s*(\"(?:[^\"\\]|\\.)*\"|'(?:[^'\\]|\\.)*')\s*=(?!=)(.*)$", ar, re.S)
            _req(mm, "a call mixes path-named and other arguments: %s" % ar.strip()[:60])
            names.append(mm.group(1))
            vals.append(mm.group(2).strip())
        body = ", ".join(vals)
        if any(k == "comment" for k, _, _ in r_segments(body)):
            body += chr(10)
        text = text[:m.start()] + "setNames(%s(%s), c(%s))" % (m.group(1), body, ", ".join(names)) + text[c + 1:]
        n += 1


def transform_r(text, depth):
    text, nnamed = rewrite_named_calls(text)
    rules, out, used = ({"P3 named-argument call": nnamed} if nnamed else {}), [], set()
    for kind, a, b in r_segments(text):
        seg = text[a:b]
        if kind == "comment":
            new = redact(seg)
            if new != seg:
                rules["P6 comment"] = rules.get("P6 comment", 0) + 1
        elif kind == "str" and PATH_RX.search(seg):
            q, body = seg[0], seg[1:-1]
            parts = split_body(body, "R")
            _req(round_trip(parts, body), "round trip failed for an R literal: %s" % seg[:80])
            new = make_expr(parts, q, "R")
            for k, _, _ in (x for x in parts if x[0] != "lit"):
                used.add(k)
                key = "P1 %s" % k if k in ("v3", "ckm", "ext") else "P2 %s" % k
                rules[key] = rules.get(key, 0) + 1
        else:
            new = seg
        out.append(new)
    new = "".join(out)
    if used:
        new = _insert_r_prologue(new, depth)
        rules["P7 prologue"] = 1
    return new, rules


R_PROLOGUE = (
    "# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is\n"
    "# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.\n"
    "if (!exists(\".P31R_ARCHIVE\")) {\n"
    "  .p31r_f <- sub(\"^--file=\", \"\", grep(\"^--file=\", commandArgs(FALSE), value = TRUE))\n"
    "  .P31R_ARCHIVE <- if (nzchar(Sys.getenv(\"P31R_ARCHIVE\"))) Sys.getenv(\"P31R_ARCHIVE\") else {\n"
    "    if (!length(.p31r_f)) stop(\"run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root\")\n"
    "    normalizePath(file.path(dirname(.p31r_f[1]), %(up)s), winslash = \"/\", mustWork = FALSE)\n"
    "  }\n"
    "  .P31R_V3 <- file.path(.P31R_ARCHIVE, \"code\")\n"
    "  .P31R_EXT <- Sys.getenv(\"P31R_EXTERNAL\", file.path(.P31R_ARCHIVE, \"external_inputs\"))\n"
    "  .P31R_CKM <- file.path(.P31R_EXT, \"CKM papers\")\n"
    "}\n")


def _insert_r_prologue(text, depth):
    lines = text.split("\n")
    k = 0
    while k < len(lines) and (not lines[k].strip() or lines[k].lstrip().startswith("#")):
        k += 1
    _req(k < len(lines), "an R file with paths and no code line")
    up = ", ".join('".."' for _ in range(depth))
    return "\n".join(lines[:k]) + ("\n" if k else "") + (R_PROLOGUE % {"up": up}) + "\n".join(lines[k:])


# ---------------------------------------------------------------------------------------------------------- Python
PY_PROLOGUE = (
    "# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is\n"
    "# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.\n"
    "import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil\n"
    "_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get(\"P31R_ARCHIVE\") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), %(up)s)).replace(chr(92), \"/\")\n"
    "_P31R_V3 = _P31R_ARCHIVE + \"/code\"\n"
    "_P31R_EXT = (_p31r_os.environ.get(\"P31R_EXTERNAL\") or _P31R_ARCHIVE + \"/external_inputs\").replace(chr(92), \"/\")\n"
    "_P31R_CKM = _P31R_EXT + \"/CKM papers\"\n"
    "_P31R_RSCRIPT = _p31r_shutil.which(\"Rscript\") or \"Rscript\"\n"
    "_P31R_FONTDIR = _p31r_os.environ.get(\"P31R_FONTDIR\") or _p31r_os.path.join(_p31r_os.environ.get(\"SystemRoot\") or \"/usr/share\", \"Fonts\")\n")


def _py_doc_ranges(tree):
    r = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant) and isinstance(n.value.value, str):
            r.append((n.lineno, n.end_lineno))
    return r


def transform_py(text, depth):
    tree = ast.parse(text)
    docs = _py_doc_ranges(tree)
    lines = text.splitlines(keepends=True)
    starts = [0]
    for ln in lines:
        starts.append(starts[-1] + len(ln))

    def off(rc):
        return starts[rc[0] - 1] + rc[1]

    edits, rules, used = [], {}, set()

    def rule(k, n=1):
        rules[k] = rules.get(k, 0) + n

    toks = list(tokenize.generate_tokens(io.StringIO(text).readline))
    in_f = 0
    # implicit concatenation ("a" "b"): a run of STRING tokens separated only by NL / comments is ONE literal; it is rewritten as one expression
    # (the run's lines are kept, so line numbers do not move)
    handled = set()
    i = 0
    while i < len(toks):
        if toks[i].type == tokenize.STRING:
            run, k = [i], i + 1
            while True:
                k2 = k
                while k2 < len(toks) and toks[k2].type in (tokenize.NL, tokenize.COMMENT):
                    k2 += 1
                if k2 < len(toks) and toks[k2].type == tokenize.STRING:
                    run.append(k2); k = k2 + 1
                else:
                    break
            if len(run) > 1 and not any(lo <= toks[run[0]].start[0] and toks[run[-1]].end[0] <= hi for lo, hi in docs):
                ms = [re.match(r"^([A-Za-z]*)(\"\"\"|'''|\"|')(.*)\2$", toks[r].string, re.S) for r in run]
                _req(all(ms), "cannot split an implicitly concatenated literal")
                body = "".join(m.group(3) for m in ms)
                if PATH_RX.search(body):
                    _req(len({(m.group(1), m.group(2)) for m in ms}) == 1 and ms[0].group(1) in ("", "r", "u") and ms[0].group(2) in ("\"", "'"),
                         "an implicitly concatenated literal with mixed or exotic quoting names a path: %s" % body[:60])
                    pre, q = ms[0].group(1), ms[0].group(2)
                    parts = split_body(body, "py")
                    _req(round_trip(parts, body), "round trip failed for an implicitly concatenated literal: %s" % body[:80])
                    expr = make_expr(parts, pre + q, "py", qclose=q)
                    a0, b0 = off(toks[run[0]].start), off(toks[run[-1]].end)
                    nl_in = text[a0:b0].count("\n")
                    if not expr.startswith("("):
                        expr = "(" + expr + ")"
                    edits.append((a0, b0, expr[:-1] + "\n" * nl_in + ")"))
                    for k_, _, _ in (x for x in parts if x[0] != "lit"):
                        used.add(k_); rule("P1 %s" % k_ if k_ in ("v3", "ckm", "ext") else "P2 %s" % k_)
                    handled.update(range(run[0], run[-1] + 1))
            i = run[-1] + 1 if len(run) > 1 else i + 1
        else:
            i += 1
    for i, tk in enumerate(toks):
        if i in handled:
            continue
        a, b = off(tk.start), off(tk.end)
        s = text[a:b]
        if tk.type == tokenize.COMMENT:
            new = redact(s)
            if new != s:
                edits.append((a, b, new)); rule("P6 comment")
        elif tk.type == tokenize.FSTRING_START:
            in_f += 1
        elif tk.type == tokenize.FSTRING_END:
            in_f -= 1
        elif tk.type == tokenize.FSTRING_MIDDLE and PATH_RX.search(s):
            if any(lo <= tk.start[0] <= hi for lo, hi in docs):
                new = redact(s); edits.append((a, b, new)); rule("P6 comment"); continue
            parts = split_body(s, "py")
            _req(round_trip(parts, s), "round trip failed for an f-string part: %s" % s[:80])
            edits.append((a, b, make_expr(parts, '"', "py", fstring=True)))
            for k, _, _ in (x for x in parts if x[0] != "lit"):
                used.add(k); rule("P1 %s" % k if k in ("v3", "ckm", "ext") else "P2 %s" % k)
        elif tk.type == tokenize.STRING and PATH_RX.search(s):
            if any(lo <= tk.start[0] and tk.end[0] <= hi for lo, hi in docs):
                new = redact(s)
                if new != s:
                    edits.append((a, b, new)); rule("P6 comment")
                continue
            m = re.match(r"^([A-Za-z]*)(\"\"\"|'''|\"|')(.*)\2$", s, re.S)
            _req(m, "cannot split the Python literal %s" % s[:60])
            pre, q, body = m.group(1), m.group(2), m.group(3)
            _req("b" not in pre.lower(), "a bytes literal names a path: %s" % s[:60])
            parts = split_body(body, "py")
            _req(round_trip(parts, body), "round trip failed for a Python literal: %s" % s[:80])
            edits.append((a, b, make_expr(parts, pre + q, "py", qclose=q)))
            for k, _, _ in (x for x in parts if x[0] != "lit"):
                used.add(k); rule("P1 %s" % k if k in ("v3", "ckm", "ext") else "P2 %s" % k)
    new = text
    for a, b, r in sorted(edits, reverse=True):
        new = new[:a] + r + new[b:]
    if used:
        node = next((n for n in tree.body if not (isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant) and isinstance(n.value.value, str))
                     and not (isinstance(n, ast.ImportFrom) and n.module == "__future__")), None)
        _req(node is not None, "a Python file with paths and no statement")
        first = min([node.lineno] + [d.lineno for d in getattr(node, "decorator_list", [])])
        nl = new.splitlines(keepends=True)
        # the edits above never add or remove a line, so the line numbers of the original still hold
        nl.insert(first - 1, PY_PROLOGUE % {"up": ", ".join('".."' for _ in range(depth))})
        new = "".join(nl)
        rules["P7 prologue"] = 1
    ast.parse(new)
    return new, rules


# ---------------------------------------------------------------------------------------------------------- special edits
def _replace_once(text, old, new, what):
    _req(text.count(old) == 1, "special edit %s: the target text occurs %d times (want exactly 1); the source changed" % (what, text.count(old)))
    return text.replace(old, new)


def special_fig_common(text):
    """P31R_fig_common.R shipped copy: public ledger, input resolver through the release map, masked ledger values refuse."""
    text = _replace_once(
        text,
        'P31R_LEDGER <- file.path(V3_ROOT, "R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv")',
        '.p31r_arch  <- normalizePath(file.path(V3_ROOT, ".."), winslash = "/", mustWork = FALSE)          # the archive root (V3_ROOT is <archive>/code)\n'
        'P31R_LEDGER <- file.path(.p31r_arch, "ledger/claims_ledger_public.tsv")                           # the PUBLIC extract (masked rows carry their reason)',
        "S1 ledger path")
    text = _replace_once(
        text,
        'p31r_use <- function(path) {\n  p <- normalizePath(path, winslash = "/", mustWork = FALSE)\n  if (!file.exists(p)) stop("data file missing: ", p)\n',
        '# an input that is not at its original place is looked up in the release map (code/P31R_INPUT_MAP.tsv): a verbatim released copy is used,\n'
        '# anything else stops. The map is keyed by the path relative to V3_ROOT (= <archive>/code).\n'
        '.p31r_map <- NULL\n'
        '.p31r_resolve <- function(p) {\n'
        '  if (file.exists(p)) return(p)\n'
        '  pre <- paste0(V3_ROOT, "/")\n'
        '  key <- if (startsWith(p, pre)) substr(p, nchar(pre) + 1, nchar(p)) else p\n'
        '  if (is.null(.p31r_map)) {\n'
        '    mf <- file.path(V3_ROOT, "P31R_INPUT_MAP.tsv")\n'
        '    if (!file.exists(mf)) stop("input not in the release: ", key, " (the release map code/P31R_INPUT_MAP.tsv is missing)", call. = FALSE)\n'
        '    .p31r_map <<- utils::read.delim(mf, sep = "\\t", quote = "", colClasses = "character", na.strings = character(0), check.names = FALSE, encoding = "UTF-8")\n'
        '  }\n'
        '  r <- .p31r_map[.p31r_map$source_path == key, , drop = FALSE]\n'
        '  if (nrow(r) == 0) stop("input not in the release: ", key, " (third-party, withheld or not a released file; see REPRODUCE.md)", call. = FALSE)\n'
        '  r <- r[1, ]\n'
        '  if (r$verbatim != "1" && r$lift != "1") stop("input not in the release: ", key, " (the release holds only a reduced or masked copy, ", r$release_path, ")", call. = FALSE)\n'
        '  q <- file.path(.p31r_arch, r$release_path)\n'
        '  if (!file.exists(q)) stop("input not in the release: ", key, " (", r$release_path, " is missing from this extraction)", call. = FALSE)\n'
        '  if (!identical(unname(tools::md5sum(q)), r$md5)) stop("released input ", r$release_path, " does not match its md5 in the release map", call. = FALSE)\n'
        '  q\n'
        '}\n'
        '# level 2 (RC1b): a lifted copy (map column lift = 1) is a declared masked or reduced stand-in for a source the generator reads. The figure-source check compares\n'
        '# the MD5 OF THE SOURCE with the approved manifest, so it asks for the source md5 recorded in the map (pinned in RELEASE_INPUTS and verified when the release was built).\n'
        '.p31r_declared_md5 <- function(p) {\n'
        '  if (file.exists(p) || is.null(.p31r_map)) return(NULL)\n'
        '  pre <- paste0(V3_ROOT, "/")\n'
        '  key <- if (startsWith(p, pre)) substr(p, nchar(pre) + 1, nchar(p)) else p\n'
        '  r <- .p31r_map[.p31r_map$source_path == key, , drop = FALSE]\n'
        '  if (nrow(r) == 0 || r$lift[1] != "1" || !nzchar(r$source_md5[1])) return(NULL)\n'
        '  r$source_md5[1]\n'
        '}\n'
        'p31r_use <- function(path) {\n  p <- .p31r_resolve(normalizePath(path, winslash = "/", mustWork = FALSE))\n',
        "S2 input resolver")
    text = _replace_once(
        text,
        '    led <- utils::read.csv(P31R_LEDGER, colClasses = "character", encoding = "UTF-8", check.names = FALSE, na.strings = character(0))\n',
        '    led <- utils::read.delim(P31R_LEDGER, sep = "\\t", quote = "", colClasses = "character", encoding = "UTF-8", check.names = FALSE, na.strings = character(0))   # public extract: TSV, JSON values unquoted\n',
        "S3 ledger reader")
    text = _replace_once(
        text,
        '  # idx: the SAME convention as the figure gate\'s registry:',
        '  if (is.list(x) && !is.null(names(x)) && identical(names(x), c("masked", "reason")))\n'
        '    stop("ledger ", id, " is masked in the public extract (", x$reason, "): ", p31r_ledger_row(id)$display, call. = FALSE)\n'
        '  # idx: the SAME convention as the figure gate\'s registry:',
        "S4 masked ledger value")
    return text


def special_af3_common(text):
    return _replace_once(
        text,
        '  have <- if (!file.exists(want)) stop("source file missing: ", want) else .md5(want)\n',
        '  have <- { q_ <- .p31r_resolve(want); sm_ <- .p31r_declared_md5(want); if (is.null(sm_)) .md5(q_) else sm_ }   # released copy (verbatim: its md5; lifted: the md5 of its source), or "input not in the release"\n',
        "S5 source md5 through the release map")


def special_figs6(text):
    """P31R_FigS6.R shipped copy (level 2): the released tiered table masks the hold-H9 cells (eas_beta_mr, eas_se_mr, eas_q of the EGFR_EAS_2P2M rows) with a text
    marker, so fread reads those columns as text. The generator drops those rows (bin <- ev[!HELD]) before it computes with them; the shipped copy converts the
    columns to numbers right after the read (the marker becomes NA; every other cell is the same number, as the frame comparison of the clean-room run shows)."""
    return _replace_once(
        text,
        'mr   <- af3_read("V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv", "a")\n',
        '# level 2 (RC1b): the released tiered table masks the hold-H9 cells with a text marker; those rows are dropped below (H9), so read the three columns as numbers (marker -> NA)\n'
        'for (.cc in c("eas_beta_mr", "eas_se_mr", "eas_q")) if (is.character(tier[[.cc]])) set(tier, j = .cc, value = suppressWarnings(as.numeric(tier[[.cc]])))\n'
        'mr   <- af3_read("V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv", "a")\n',
        "S6 numeric H9 columns")


def special_decode_v5_script(text):
    """S2-D41 W3: the comment of 01_decode_replication_mr.R quotes two deCODE exposure F values (TINAGL1 ~ CAD, rs945214: the F of the discarded and of the kept contrast).
    They are the withheld decode_F_stat of that pair (TIER_RULES F15) and the manuscript does not print either, so the shipped copy names the quantity and not its values."""
    return _replace_once(
        text,
        "57" + ".465 (discarded T/A) in place of 163" + ".857 (kept G/A)",        # pieces: this file is shipped and its text is leak-scanned (R14)
        "[value withheld, release rule F15] (discarded T/A) in place of [value withheld, release rule F15] (kept G/A)",
        "S8 deCODE F values in a comment")


def special_b01_ld(text):
    """S2-D41 W4: b01_ld.R writes the note of workstream_A/results.tsv row variant.A.rsid with the UKB-PPP marginal z of the ARG1 protein at signal 1 typed into the string. The shipped table
    withholds that number (TIER_RULES N01); the shipped copy of the script names the quantity and not the value, so it writes the note the shipped table carries."""
    return _replace_once(
        text,
        "(z_marginal -24" + ".066; pos37 131897278)",
        "(z_marginal withheld_pending_terms; pos37 131897278)",
        "S9 marginal z in a note string")


def special_w1_compare_seen(text):
    """S2-D42: w1_compare_seen.R carries the RAD23B deCODE F as a typed 'seen' literal. The F is withheld (TIER_RULES F15); the shipped copy names the quantity, not the value."""
    return _replace_once(
        text,
        'add("S1.rad23b.decode_F", "10' + '.14", L,',
        'add("S1.rad23b.decode_F", "[value withheld, release rule F15]", L,',
        "S10 deCODE F literal in a comparison script")


SPECIAL = {"R2_rewrite/figures/scripts/P31R_fig_common.R": special_fig_common,
           "R1_stage1/S1R1_R2_01_rad23b/w1_compare_seen.R": special_w1_compare_seen,
           "R1_stage1/S1R1_R2_12_chen_instrument/scripts/b01_ld.R": special_b01_ld,
           "R2_rewrite/figures/af3/P31R_af3_common.R": special_af3_common,
           "R2_rewrite/figures/scripts/P31R_FigS6.R": special_figs6,
           "paper3_deposit_request_v5/scripts/scripts/_shared/decode/01_decode_replication_mr.R": special_decode_v5_script}


# ---------------------------------------------------------------------------------------------------------- one file
def transform_file(path, data, ship_path):
    """path: project-relative; data: the source bytes. Returns (new bytes, {rule: count}). Raises when an absolute local path survives."""
    text = data.decode("utf-8")
    nl = "\r\n" if "\r\n" in text else "\n"
    work = text.replace("\r\n", "\n")
    rules = {}
    sp = SPECIAL.get(path)
    if sp:
        work = sp(work)
        rules["S special edit"] = 1
    depth = ship_path.count("/")            # code/a/b/f -> archive is 'depth' levels above the file's folder
    ext = path.rsplit(".", 1)[-1].lower()
    if ext == "r":
        work, r = transform_r(work, depth)
    elif ext == "py":
        work, r = transform_py(work, depth)
    elif ext in ("tsv", "csv", "json", "md", "txt"):
        new = redact_cell(work)
        r = {"P10 data cell": 1} if new != work else {}
        work = new
    else:
        raise PortabilityError("no transform for %s" % path)
    for k, v in r.items():
        rules[k] = rules.get(k, 0) + v
    lk = leaks(work)
    _req(not lk, "%s: an absolute local path survives the transform at line %d: %s" % (path, lk[0][0], lk[0][1]) if lk else "")
    return work.replace("\n", nl).encode("utf-8"), rules


# ---------------------------------------------------------------------------------------------------------- every other shipped file (RC1 fix K01/K02)
EXEC_EXT = ("r", "py", "sh")                                    # executables anywhere in out/ come under the same transform as code/
TEXT_EXT = ("tsv", "csv", "json", "md", "txt", "log")           # shipped text tables and documents: absolute local paths in cells are redacted
EXTRA_SKIP_PREFIX = ("code/", "ledger/")                        # code/ is the loop above; ledger/ is written later, already in the redacted form
EXTRA_SKIP_FILES = ("RELEASE_MANIFEST.tsv", "README.md", "REPRODUCE.md")        # written by rel_docs after this step, from redacted sources


def redact_public(text):
    """Data cells and documents outside code/: the redacted relative form of the public ledger (rel_ledger.public_source): a project path becomes project-relative,
    a drive path becomes <external>/..., anything else of the forms in PATH_RX goes through redact()."""
    s = text.replace(_CKM + "/V3/", "").replace(_CKM + "/", "")
    s = re.sub(r"(?<![A-Za-z0-9_])[CDE]:/(?=[A-Za-z0-9_$.~])", "<external>/", s)
    return redact(s)


def extra_paths(out):
    out = Path(out)
    res = []
    for p in sorted(out.rglob("*")):
        if not p.is_file():
            continue
        k = p.relative_to(out).as_posix()
        if k.startswith(EXTRA_SKIP_PREFIX) or k in EXTRA_SKIP_FILES:
            continue
        if k.rsplit(".", 1)[-1].lower() in EXEC_EXT + TEXT_EXT:
            res.append(k)
    return res


def transform_extra(k, data):
    """k: path inside out/. Returns (new bytes, {rule: count}); the bytes are unchanged when the file names no absolute local path. Raises when one survives."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return data, {}
    if not leaks(text):
        return data, {}
    nl = "\r\n" if "\r\n" in text else "\n"
    work = text.replace("\r\n", "\n")
    ext = k.rsplit(".", 1)[-1].lower()
    depth = k.count("/")
    if ext == "r":
        work, rules = transform_r(work, depth)
    elif ext == "py":
        work, rules = transform_py(work, depth)
    elif ext == "sh":
        raise PortabilityError("%s: no transform for a shell script that names an absolute local path" % k)
    else:
        new = redact_public(work)
        rules = {"P11 data cell / document (public relative form)": 1} if new != work else {}
        work = new
    lk = leaks(work)
    _req(not lk, "%s: an absolute local path survives the transform at line %d: %s" % (k, lk[0][0], lk[0][1]) if lk else "")
    return work.replace("\n", nl).encode("utf-8"), rules


# ---------------------------------------------------------------------------------------------------------- the build step
def _read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def code_rows(v3):
    return [r for r in _read_tsv(Path(v3) / CODE_REL) if r["ship"] == "1"]


def input_map_rows(v3, out):
    """One row per released NON-code file that has a source path: RELEASE_INPUTS copies plus module outputs it declares. verbatim = the released file's md5 equals the
    source file's md5 NOW (a masked or reduced copy is listed, with verbatim 0, so a generator can say why it cannot use it)."""
    v3, out = Path(v3), Path(out)
    rows, seen = [], set()
    for r in _read_tsv(v3 / RI_REL):
        rp, sp = r["release_path"], r["source_path"]
        if r["ship_status"] != "ship" or not rp or not sp or rp.startswith("code/") or rp in seen:
            continue
        f, src = out / rp, v3 / sp
        if not f.is_file():
            continue
        seen.add(rp)
        smd5 = md5_bytes(src.read_bytes()) if src.is_file() and src.stat().st_size < 200_000_000 else ""
        m = md5_bytes(f.read_bytes())
        rows.append(dict(source_path=sp.replace("\\", "/"), release_path=rp, md5=m, source_md5=smd5, verbatim="1" if smd5 and smd5 == m else "0", lift="0", lift_units=""))
    return sorted(rows, key=lambda x: x["source_path"])


# ---------------------------------------------------------------------------------------------------------- level 2 (RC1b)
def _delim(path):
    return "\t" if str(path).lower().endswith(".tsv") else ","


def l2_rows(path=None):
    """The declared lifts (inputs/LEVEL2_INPUTS.tsv of the release folder this module was imported from), validated."""
    p = Path(path) if path else _L2_FILE
    if not p.is_file():
        return []
    rows = _read_tsv(p)
    _req(not rows or list(rows[0].keys()) == L2_COLS, "L2-V1 %s: columns are %s, want %s" % (p.name, list(rows[0].keys()) if rows else "?", L2_COLS))
    ids, srcs, rels = set(), set(), set()
    for r in rows:
        _req(r["l2_id"] and r["l2_id"] not in ids, "L2-V2 duplicate or empty l2_id %r" % r["l2_id"])
        ids.add(r["l2_id"])
        _req(r["kind"] in L2_KINDS, "L2-V3 %s: kind %r is not one of %s" % (r["l2_id"], r["kind"], L2_KINDS))
        for k in ("source_path", "release_path", "units", "reason"):
            _req(r[k].strip(), "L2-V4 %s: empty %s" % (r["l2_id"], k))
        for k in ("source_path", "release_path"):
            _req(".." not in r[k].split("/") and not re.match(r"[A-Za-z]:|/", r[k]), "L2-V5 %s: %s %r is not a relative path" % (r["l2_id"], k, r[k]))
        _req(r["source_path"] not in srcs, "L2-V6 %s: source %s is declared twice" % (r["l2_id"], r["source_path"]))
        _req(r["release_path"] not in rels, "L2-V7 %s: release path %s is declared twice" % (r["l2_id"], r["release_path"]))
        srcs.add(r["source_path"])
        rels.add(r["release_path"])
        if r["kind"] == "extract":
            _req(r["release_path"].startswith(L2_PREFIX), "L2-V8 %s: an extract is written under %s only" % (r["l2_id"], L2_PREFIX))
            _req(r["spec"].partition(":")[0] in ("keep", "drop", "section"), "L2-V9 %s: spec %r is not keep:/drop:/section:" % (r["l2_id"], r["spec"]))
        else:
            _req(not r["release_path"].startswith(L2_PREFIX) and r["needs"].strip(), "L2-V10 %s: a released_copy names a file another step ships and the columns the generator needs" % r["l2_id"])
    return rows


def extract_bytes(spec, data, source_path):
    """spec: keep:a|b|c (those columns, in that order) | drop:a|b (all but those) | section:START~~END (the lines from the one line that starts with START to the one that starts with END, both included).
    Returns (bytes, data rows or lines, columns). Every named column or anchor must exist exactly as written; the output is re-read and compared with the projection of the source."""
    kind, _, arg = spec.partition(":")
    text = data.decode("utf-8").replace("\r\n", "\n")
    if kind == "section":
        a, _, b = arg.partition("~~")
        ln = text.split("\n")
        ia = [i for i, x in enumerate(ln) if x.startswith(a)]
        _req(a and b and len(ia) == 1, "L2-X1 section start %r occurs %d times in %s (want exactly one)" % (a, len(ia), source_path))
        ib = [i for i, x in enumerate(ln) if x.startswith(b) and i > ia[0]]
        _req(len(ib) >= 1 and sum(1 for x in ln if x.startswith(b)) == 1, "L2-X2 section end %r is not unique after the start in %s" % (b, source_path))
        out = ln[ia[0]:ib[0] + 1]
        return ("\n".join(out) + "\n").encode("utf-8"), len(out), 1
    _req(kind in ("keep", "drop"), "L2-X3 unknown spec %r" % spec)
    d = _delim(source_path)
    rows = [r for r in csv.reader(io.StringIO(text), delimiter=d) if r]
    hdr = rows[0]
    cols = arg.split("|")
    _req(all(cols) and len(set(cols)) == len(cols), "L2-X4 spec %r names an empty or repeated column" % spec)
    _req(len(set(hdr)) == len(hdr), "L2-X5 %s has a repeated header name" % source_path)
    miss = [c for c in cols if c not in hdr]
    _req(not miss, "L2-X6 %s has no column %s" % (source_path, miss))
    _req(all(len(r) == len(hdr) for r in rows), "L2-X7 %s has a ragged row" % source_path)
    idx = [hdr.index(c) for c in cols] if kind == "keep" else [i for i, c in enumerate(hdr) if c not in cols]
    _req(idx, "L2-X8 the extract of %s would hold no column" % source_path)
    proj = [[r[i] for i in idx] for r in rows]
    buf = io.StringIO()
    csv.writer(buf, delimiter=d, lineterminator="\n").writerows(proj)
    blob = buf.getvalue().encode("utf-8")
    back = [r for r in csv.reader(io.StringIO(blob.decode("utf-8")), delimiter=d) if r]
    _req(back == proj, "L2-X9 the extract of %s does not read back to the projection of its source" % source_path)
    return blob, len(proj) - 1, len(idx)


def _ri_field(field, source_path):
    """A field of the RELEASE_INPUTS row (of the release folder this module came from) that declares source_path; '' when none."""
    if not _RI_FILE.is_file():
        return ""
    hit = [r for r in _read_tsv(_RI_FILE) if r["source_path"].replace("\\", "/") == source_path and r["role"] in ("input", "copy")]
    return hit[0][field] if hit else ""


def _header(f):
    with open(f, encoding="utf-8", newline="") as fh:
        first = fh.readline().rstrip("\r\n")
    return next(csv.reader([first], delimiter=_delim(f)), [])     # RC1b: a masked copy may carry quoted header cells (R fread reads them as the bare names)


def level2(v3, out, mp, rows=None):
    """Apply the declared lifts to the map rows mp (list of dicts). Writes the extracts under out/code/lift_inputs/. Returns (mp, one record per extract)."""
    v3, out = Path(v3), Path(out)
    by_src = {r["source_path"]: r for r in mp}
    new = []
    for r in (rows if rows is not None else l2_rows()):
        src, rp = r["source_path"], r["release_path"]
        needs = [c for c in r["needs"].split("|") if c]
        f = out / rp
        if r["kind"] == "released_copy":
            _req(f.is_file(), "L2-R1 %s: %s is not in out/ (another step must ship it before this one runs)" % (r["l2_id"], rp))
            miss = [c for c in needs if c not in _header(f)]
            _req(not miss, "L2-R2 %s: %s lacks the column(s) %s that the generator reads" % (r["l2_id"], rp, miss))
            sf = v3 / src
            _req(sf.is_file(), "L2-R3 %s: the source %s is missing" % (r["l2_id"], src))
            smd5 = md5_bytes(sf.read_bytes())
        else:
            sf = v3 / src
            _req(sf.is_file(), "L2-E1 %s: the source %s is missing" % (r["l2_id"], src))
            data = sf.read_bytes()
            smd5 = md5_bytes(data)
            pin = _ri_field("source_md5", src)
            _req(not pin or pin == smd5, "L2-E2 %s: %s has md5 %s, RELEASE_INPUTS pins %s (run --pin)" % (r["l2_id"], src, smd5[:8], pin[:8]))
            blob, nr, nc = extract_bytes(r["spec"], data, src)
            miss = [c for c in needs if c not in _header_of(blob, src, r["spec"])]
            _req(not miss, "L2-E3 %s: the extract lacks the column(s) %s that the generator reads" % (r["l2_id"], miss))
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_bytes(blob)
            lic = _ri_field("licence_class", src)
            _req(lic, "L2-E4 %s: %s has no licence_class in RELEASE_INPUTS (declare the source as an input of rel_portability)" % (r["l2_id"], src))
            new.append(dict(release_path=rp, ship_status="ship", source_paths=[src, L2_REL], producer="rel_portability.py", licence_class=lic, rows=nr, columns=nc,
                            note="level-2 lift %s (%s): reduced extract of %s [%s]; for %s. %s" % (r["l2_id"], r["kind"], src, r["spec"], r["units"].replace(";", ", "), r["reason"])))
        m = md5_bytes(f.read_bytes())
        row = dict(source_path=src, release_path=rp, md5=m, source_md5=smd5, verbatim="1" if smd5 == m else "0", lift="1", lift_units=r["units"])
        if src in by_src:
            _req(by_src[src]["release_path"] == rp, "L2-M1 %s: the map sends %s to %s, LEVEL2_INPUTS to %s" % (r["l2_id"], src, by_src[src]["release_path"], rp))
            by_src[src].update(row)
        else:
            mp.append(row)
            by_src[src] = row
    return sorted(mp, key=lambda x: x["source_path"]), new


def _header_of(blob, source_path, spec):
    if spec.startswith("section:"):
        return []
    return blob.decode("utf-8").split("\n", 1)[0].split(_delim(source_path))


def write_tsv(path, cols, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])


def run(v3, out, build=None):
    """Transform every shipped code file in out/code/ from its source; write code/P31R_INPUT_MAP.tsv and build/portability_log.tsv.
    Returns {'updates': [{release_path, md5, source_md5, rules}], 'new': [module-style records for the files this step ADDED]}.
    Order in the build: copy -> tiers/withhold -> THIS -> ledger extract -> docs."""
    v3, out = Path(v3), Path(out)
    updates, log = [], []
    for r in code_rows(v3):
        src, dst = v3 / r["path"], out / r["ship_path"]
        _req(src.is_file(), "source %s is missing" % r["path"])
        _req(dst.is_file(), "the code copy %s is not in out/ (run the copy step first)" % r["ship_path"])
        data = src.read_bytes()
        new, rules = transform_file(r["path"].replace("\\", "/"), data, r["ship_path"])
        dst.write_bytes(new)
        sm, nm = md5_bytes(data), md5_bytes(new)
        log.append(dict(path=r["ship_path"], source_md5=sm, shipped_md5=nm, rules="; ".join("%s x%d" % (k, v) for k, v in sorted(rules.items())) or "none (byte-identical)"))
        updates.append(dict(release_path=r["ship_path"], md5=nm, source_md5=sm, rules=log[-1]["rules"]))
    decl = {r["release_path"]: r["source_md5"] for r in _read_tsv(v3 / RI_REL) if r["role"] == "copy" and r["release_path"]}
    for k in extra_paths(out):                 # every other shipped executable and text table (K01/K02): the copy as built so far is the input of this step
        f = out / k
        data = f.read_bytes()
        new, rules = transform_extra(k, data)
        if new == data:
            continue
        f.write_bytes(new)
        sm, nm = md5_bytes(data), md5_bytes(new)
        rl = "; ".join("%s x%d" % (a, b) for a, b in sorted(rules.items()))
        dm = decl.get(k, "")
        log.append(dict(path=k, source_md5=sm, shipped_md5=nm, rules="X0 outside code/ (source_md5 = the copy before this step%s); " % ("; the declared source file has md5 %s" % dm[:8] if dm else "") + rl))
        updates.append(dict(release_path=k, md5=nm, source_md5=sm, rules=log[-1]["rules"], extra=True))
    mp = input_map_rows(v3, out)
    mp, l2new = level2(v3, out, mp)               # RC1b: declared lifts (extracts written first, so that the map md5 is the shipped md5)
    write_tsv(out / MAP_REL, MAP_COLS, mp)
    b = Path(build) if build else out.parent / "build"
    write_tsv(b / "portability_log.tsv", LOG_COLS, log)
    write_tsv(out / LOG_REL, LOG_COLS, log)          # shipped too: the md5 pins inside the scripts are SOURCE md5s; this maps them to the shipped copies
    new = [dict(release_path=LOG_REL, ship_status="ship", source_paths=[CODE_REL], producer="rel_portability.py", licence_class="own_derived", rows=len(log), columns=len(LOG_COLS),
                note="per shipped code file: source md5, shipped md5 and the portability rules applied (md5 pins inside the scripts name the SOURCE md5)"),
           dict(release_path=MAP_REL, ship_status="ship", source_paths=[CODE_REL], producer="rel_portability.py", licence_class="own_derived", rows=len(mp), columns=len(MAP_COLS),
                note="release map of the shipped code: project path of an input -> its path in this archive, md5, and whether the released copy is verbatim (rel_portability.py)")]
    return dict(updates=updates, new=new + l2new)


# ---------------------------------------------------------------------------------------------------------- module interface (as rel_screen / rel_tables / rel_ledger)
INPUTS = [CODE_REL] + ([L2_REL] + sorted({r["source_path"] for r in l2_rows() if r["kind"] == "extract"}) if _L2_FILE.is_file() else [])
LAST_UPDATES = []        # set by build(): the shipped md5 of every transformed code copy, for the core to put into its build records


def build(v3, out, withheld, ledger):
    """Spec interface: transform the shipped copies in out/code/, write code/P31R_INPUT_MAP.tsv and code/P31R_PORTABILITY_LOG.tsv; return one record per NEW file."""
    global LAST_UPDATES
    res = run(v3, out)
    LAST_UPDATES = res["updates"]
    return res["new"]


# ---------------------------------------------------------------------------------------------------------- self-test
def selftest(v3):
    v3 = Path(v3)
    bad = []

    def rep(name, ok, detail=""):
        print("  %-66s %s %s" % (name, "ok" if ok else "MISS", detail))
        if not ok:
            bad.append(name)

    def trips(name, fn, token):
        try:
            fn()
            rep(name, False, "no check fired")
        except PortabilityError as e:
            rep(name, token in str(e), str(e)[:80])

    q = '"' + _CKM + "/V3/R1_stage1/x/%s.tsv" + '"'
    r_src = "# Run: Rscript " + _RS_EXE + ".exe " + q + "\nV3 <- " + '"' + _CKM + "/V3" + '"' + "\nf <- sprintf(" + q + ", 1)\nd <- " + '"D' + ':/data/a b/c.gz"' + "\n"
    new, rules = transform_r(r_src, 1)
    rep("R: executable literals become expressions, comment redacted", ".P31R_V3" in new and "<archive>/code/R1_stage1" in new and "Rscript" in new and not leaks(new), str(rules))
    rep("R: a drive path goes under the external root", 'paste0(.P31R_EXT, "/D/data/a b/c.gz")' in new)
    rep("R: sprintf template survives (paste0 + format)", 'sprintf(paste0(.P31R_V3, "/R1_stage1/x/%s.tsv"), 1)' in new)
    p_src = '"""doc ' + _CKM + '/V3/x"""\nimport os\nP = "' + _CKM + '/V3/a/b.tsv"\nQ = f"' + _CKM + '/V3/{P}/c"\nR = "' + "D" + ':/data/x"  # ' + _CKM + "\n"
    new, rules = transform_py(p_src, 1)
    rep("Python: literals, f-string, comment and docstring", "_P31R_V3" in new and not leaks(new) and "{_P31R_V3}" in new, str(rules))
    bs = "\\"
    trips("M1 a Windows backslash path (a form no rule matches) is refused, not shipped",
          lambda: transform_file("R2_rewrite/figures/scripts/z.R", ("x <- 'C" + ":" + bs + bs + "Users" + bs + bs + "example" + bs + bs + "q'\n").encode(), "code/R2_rewrite/figures/scripts/z.R"), "survives the transform")
    # S6.3b2a (D4): the two declared extensions, each with a mutation that must trip the check it names
    pb = ("import os\nB = r\"E" + ":" + bs + "20260331-CKM" + bs + "Paper 3\"\nOUT = os.path.join(B, \"results\")\n")
    new, rules = transform_py(pb, 1)
    rep("D4 Python: a raw string with a backslash drive root becomes the external root with '/' separators",
        '(_P31R_EXT + r"/E/20260331-CKM/Paper 3")' in new and not leaks(new) and rules.get("P1 ext") == 1, str(rules))
    pr = ("import subprocess\n_r = subprocess.run([r\"C" + ":" + bs + "Program Files" + bs + "R" + bs + "R-4.5.1" + bs + "bin" + bs + "Rscript.exe\", \"-e\"])\n")
    new, rules = transform_py(pr, 1)
    rep("D4 Python: an Rscript path of another version and backslashes becomes the Rscript expression", "subprocess.run([_P31R_RSCRIPT, " in new and not leaks(new), str(rules))
    trips("M7 a drive other than C, D, E in a backslash raw string is refused, not shipped",
          lambda: transform_file("R2_rewrite/figures/scripts/z.py", ("import os\nB = r\"F" + ":" + bs + "x\"\n").encode(), "code/R2_rewrite/figures/scripts/z.py"), "survives the transform")
    rep("D4 leaks(): an escape sequence in a message string is not a drive path, a real backslash path still is",
        not leaks("sys.exit(\"61: could not read R" + ":" + bs + "n\" + x)") and bool(leaks("p = r\"C" + ":" + bs + "temp" + bs + "x\"")) and bool(leaks("p = \"D" + ":" + bs + "n" + bs + "x\"")))
    trips("M2 a round trip that does not reproduce the literal",
          lambda: _req(round_trip([("v3", "", _CKM + "/V3"), ("lit", "/x", "/y")], _CKM + "/V3/x"), "round trip failed"), "round trip")
    trips("M3 special edit whose anchor changed", lambda: special_af3_common("nothing here\n"), "occurs 0 times")
    real = (v3 / "R2_rewrite/figures/scripts/P31R_fig_common.R").read_text(encoding="utf-8")
    done = special_fig_common(real.replace("\r\n", "\n"))
    rep("fig_common special edit: public ledger, resolver, TSV reader, masked-value stop",
        all(s in done for s in ("claims_ledger_public.tsv", "input not in the release", "read.delim(P31R_LEDGER", "is masked in the public extract")))
    trips("M4 fig_common source drifts (anchor edited)", lambda: special_fig_common(real.replace("p31r_use <- function(path) {", "p31r_use <- function(path, x) {")), "occurs 0 times")
    d_src = ("a\tb\n1\t" + _CKM + "/V3/V4_local/raw/x.tsv\n2\tE" + ":/20260331-CKM/Paper 3/y.rds\n").encode()
    nb, rl = transform_extra("workstreams/W2/results.tsv", d_src)
    rep("extra: a table cell goes to the public relative form", b"\tV4_local/raw/x.tsv" in nb and b"<external>/20260331-CKM/Paper 3/y.rds" in nb and not leaks(nb.decode()), str(rl))
    nb, rl = transform_extra("decoy/z.R", ('V3 <- "' + _CKM + '/V3"\nx <- 1\n').encode())
    rep("extra: a shipped R script outside code/ gets the same expressions and prologue", b".P31R_V3" in nb and b"P31R_ARCHIVE" in nb and not leaks(nb.decode()), str(rl))
    nb, rl = transform_extra("decoy/ok.tsv", b"a\tb\n1\t2\n")
    rep("extra: a file without a local path is byte-identical", nb == b"a\tb\n1\t2\n" and not rl)
    trips("M5 an absolute path of a form no rule matches, in a table outside code/", lambda: transform_extra("x/t.tsv", b"a\nZ" + b":/odd/path\n"), "survives the transform")
    # ---- level 2 (RC1b): the extract grammar, the table validator, the shipped resolver and the S6 edit; each mutation must trip the check it names
    blob, nr, nc = extract_bytes("keep:b|a", b"a,b,c\n1,2,3\n4,5,6\n", "x.csv")
    rep("L2 extract keep: the columns in the order named, every row", blob == b"b,a\n2,1\n5,4\n" and (nr, nc) == (2, 2), repr(blob))
    blob, nr, nc = extract_bytes("drop:b", b"a\tb\tc\r\n1\t2\t3\r\n", "x.tsv")
    rep("L2 extract drop: all but the columns named (tab-separated, CRLF in)", blob == b"a\tc\n1\t3\n", repr(blob))
    blob, nr, nc = extract_bytes("section:## A~~## C", b"x\n## A\nk\n## B\n## C\ny\n", "x.md")
    rep("L2 extract section: from the start line to the end line, both included", blob == b"## A\nk\n## B\n## C\n" and nr == 4, repr(blob))
    trips("L2-X6 a named column that the source lacks", lambda: extract_bytes("keep:zz", b"a,b\n1,2\n", "x.csv"), "L2-X6")
    trips("L2-X1 a section start that occurs twice", lambda: extract_bytes("section:## A~~## C", b"## A\n## A\n## C\n", "x.md"), "L2-X1")
    trips("L2-X2 a section end that is not unique", lambda: extract_bytes("section:## A~~## C", b"## A\n## C\n## C\n", "x.md"), "L2-X2")
    trips("L2-X8 a drop that leaves no column", lambda: extract_bytes("drop:a|b", b"a,b\n1,2\n", "x.csv"), "L2-X8")
    trips("L2-X5 a repeated header name", lambda: extract_bytes("keep:a", b"a,a\n1,2\n", "x.csv"), "L2-X5")
    trips("L2-X7 a ragged row", lambda: extract_bytes("keep:a", b"a,b\n1\n", "x.csv"), "L2-X7")
    trips("L2-X3 an unknown spec", lambda: extract_bytes("take:a", b"a\n1\n", "x.csv"), "L2-X3")
    good = ["l2_id\tunits\tkind\tsource_path\trelease_path\tspec\tneeds\treason",
            "L2-T1\tFigure_1\textract\tR1/a.tsv\tcode/lift_inputs/a.tsv\tkeep:a\ta\twhy",
            "L2-T2\tFigure_1\treleased_copy\tR1/b.tsv\tx/b.tsv\t\ta\twhy"]

    def l2_check(lines):
        with tempfile.TemporaryDirectory() as td:
            fp = Path(td) / "LEVEL2_INPUTS.tsv"
            fp.write_text("\n".join(lines) + "\n", encoding="utf-8")
            return l2_rows(fp)
    rep("L2 table: a valid table reads", len(l2_check(good)) == 2)
    trips("L2-V2 a repeated l2_id", lambda: l2_check(good + [good[1].replace("R1/a.tsv", "R1/c.tsv").replace("lift_inputs/a", "lift_inputs/c")]), "L2-V2")
    trips("L2-V3 an unknown kind", lambda: l2_check([good[0], good[1].replace("extract", "copy")]), "L2-V3")
    trips("L2-V5 an absolute source path", lambda: l2_check([good[0], good[1].replace("R1/a.tsv", (_P31R_EXT + "/C/x/a.tsv"))]), "L2-V5")
    trips("L2-V6 the same source twice", lambda: l2_check(good + [good[1].replace("L2-T1", "L2-T3").replace("lift_inputs/a", "lift_inputs/c")]), "L2-V6")
    trips("L2-V8 an extract outside code/lift_inputs/", lambda: l2_check([good[0], good[1].replace("code/lift_inputs/a.tsv", "x/a.tsv")]), "L2-V8")
    trips("L2-V10 a released_copy without the columns the generator needs", lambda: l2_check([good[0], good[2].replace("\ta\twhy", "\t\twhy")]), "L2-V10")
    rep("level 2 resolver: a lifted copy is accepted, its source md5 is answered, a non-lifted reduced copy still stops",
        all(s in done for s in ('r$verbatim != "1" && r$lift != "1"', ".p31r_declared_md5 <- function(p)", "the release holds only a reduced or masked copy")))
    af3 = special_af3_common(
        "x\n  have <- if (!file.exists(want)) stop(\"source file missing: \", want) else .md5(want)\ny\n")
    rep("level 2 figure-source check: a lifted copy answers with its source md5", ".p31r_declared_md5(want)" in af3 and ".md5(q_)" in af3)
    trips("M6 the FigS6 edit whose anchor changed", lambda: special_figs6("nothing here\n"), "occurs 0 times")
    s6 = special_figs6('a\nmr   <- af3_read("V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv", "a")\nb\n')
    rep("FigS6 edit: the three H9 columns are read as numbers, inserted before the screen read", s6.count("eas_beta_mr") == 1 and s6.index("as.numeric") < s6.index("mr   <- af3_read"))
    with tempfile.TemporaryDirectory() as td:
        tv, to = Path(td) / "v3", Path(td) / "out"
        (tv / "R1").mkdir(parents=True)
        (to / "x").mkdir(parents=True)
        (tv / "R1" / "b.tsv").write_text("a\tb\n1\t2\n", encoding="utf-8")
        (to / "x" / "b.tsv").write_text("a\tb\n1\t2\n", encoding="utf-8")
        row = dict(l2_id="L2-T2", units="Figure_1", kind="released_copy", source_path="R1/b.tsv", release_path="x/b.tsv", spec="", needs="a|b", reason="why")
        m, _ = level2(tv, to, [], rows=[row])
        rep("L2 released_copy: a lifted map row, verbatim when the md5s agree", len(m) == 1 and m[0]["lift"] == "1" and m[0]["verbatim"] == "1" and m[0]["lift_units"] == "Figure_1", str(m))
        trips("L2-R1 the released file is not in out/", lambda: level2(tv, to, [], rows=[dict(row, release_path="x/none.tsv")]), "L2-R1")
        trips("L2-R2 the released file lacks a column the generator needs", lambda: level2(tv, to, [], rows=[dict(row, needs="a|zz")]), "L2-R2")
        trips("L2-M1 the map sends the source somewhere else", lambda: level2(tv, to, [dict(source_path="R1/b.tsv", release_path="y/other.tsv", md5="", source_md5="", verbatim="0", lift="0", lift_units="")], rows=[row]), "L2-M1")
    print("rel_portability selftest: %s" % ("PASS" if not bad else "FAIL: %s" % bad))
    return not bad


if __name__ == "__main__":
    import sys
    v3_ = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[2])
    sys.exit(0 if selftest(v3_) else 1)
