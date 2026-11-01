"""Stage 2 step 2: RE-APPLICATION of the Stage 1b KEY-ID NORMALISER (S6 N_tools/rewrite_06b_branch_inputs.py) to the regenerated
branches.tsv, as a step of run_all.R (10_branches.R writes value-bearing `inputs` cells, so every rerun needs it again).
Differences from the Stage 1b script, and only these: (i) no pinned pre-rewrite/results md5 (the files are regenerated; the
pre-rewrite file is kept in fix_N_prefix_step2/ with its md5); (ii) the 6(d) and Rule 3(i) rows also list the four-fit sensitivity
ids (S6.d.<cid>.s_max_4fit; s_protein_4fit, s_outcome_4fit); (iii) it refuses to run when the file is already normalised.
Original docstring follows.

Stage 1b KEY-ID NORMALISER, section S6, step 1 (STAGE1B_BRIEFS_v3.md, "Order of work" 1).

W6b's branches.tsv `inputs` column holds values and descriptions, not quantity_ids. This script rewrites that column, and
only that column, as quantity_ids. Ids from W6b's own results.tsv are written bare; ids from W6a's results.tsv carry the
prefix "06a:". Nothing else in the file changes (asserted byte for byte on every other field, and by a round-trip check of
the writer on the original file). The original file is kept in fix_N_prefix/ with its md5.

Checks (each with a mutation that must trip it):
  C1 every id written exists exactly once in the results.tsv it names;
  C2 every field outside `inputs` is unchanged, and the row count is unchanged;
  C3 the rewritten `inputs` column is value-free: each cell is a '; '-separated list of ids and nothing else;
  C4 the writer reproduces the original file byte for byte (so quoting and line endings are preserved).
Usage: python -I rewrite_06b_branch_inputs.py            (refuses to run if branches.tsv no longer has the pre-rewrite md5)
"""
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import csv, hashlib, io, os, re, shutil, sys

R1 = (_P31R_V3 + "/R1_stage1")
B6B = R1 + "/S1R1_R2_06b_renal_d_umod/branches.tsv"
RES6B = R1 + "/S1R1_R2_06b_renal_d_umod/results.tsv"
RES6A = R1 + "/S1R1_R2_06a_renal/results.tsv"
PREFIX_DIR = R1 + "/S1R1_R2_06b_renal_d_umod/fix_N_prefix_step2"
LOG = R1 + "/S1R1_R2_06b_renal_d_umod/work/10b_normalise_inputs.log"
PRE_MD5 = None      # regenerated at every run; recorded in fix_N_prefix_step2/MD5SUMS.txt
RES6A_MD5 = "2a85c748f20495b74227cb4780a5f759"
RES6B_MD5 = None

out = []
def say(s):
    print(s); out.append(s)

def md5b(b): return hashlib.md5(b).hexdigest()
def md5f(p):
    with open(p, "rb") as f: return md5b(f.read())

def parse(b):
    return list(csv.reader(io.StringIO(b.decode("utf-8"), newline=""), delimiter="\t", quotechar='"'))

def serialise(rows):
    s = io.StringIO(newline="")
    w = csv.writer(s, delimiter="\t", quotechar='"', quoting=csv.QUOTE_ALL, lineterminator="\r\n")
    w.writerows(rows)
    return s.getvalue().encode("utf-8")

def ids_of(path):
    with open(path, encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        ids = [row["quantity_id"] for row in r]
    return ids

def mapping(rows):
    """Return the new inputs cell for every data row. Each row is identified by its rule_ref and the leading token of the
    old inputs cell (a candidate id or a signal label), never by a value."""
    new = []
    for i, row in enumerate(rows[1:], start=1):
        rule, inp = row[0], row[1]
        if rule == "6(d); Rule 3(i)":
            m = re.match(r"^([A-Za-z0-9]+_[A-Za-z0-9]+): s_protein ", inp)
            if not m: raise SystemExit("row %d: unrecognised 6(d) inputs" % i)
            c = m.group(1)
            new.append("S6.d.%s.s_protein; S6.d.%s.s_outcome; S6.d.%s.s_max; S6.d.%s.s_max_4fit" % (c, c, c, c))
        elif rule == "6(d)":
            new.append("; ".join(["S6.d.%s" % k for k in "abcd"] + ["S6.d.kept.%s" % k for k in "abcd"]))
        elif rule == "11 D6":
            new.append("06a:S6.primary.or_lo; S6.d.lo")
        elif rule == "6 renal OR wording; 11 D6":
            new.append("; ".join(["06a:S6.primary.or_mle", "06a:S6.primary.or_lo", "06a:S6.primary.or_hi", "06a:S6.primary.p",
                                  "06a:S6.b.locus_outcome.or_lo", "06a:S6.b.locus_outcome.or_hi", "06a:S6.c.pooled_lo", "06a:S6.c.pooled_hi",
                                  "06a:S6.b.a_pair.or_lo", "06a:S6.b.a_pair.or_hi", "S6.d.lo", "S6.d.hi"]))
        elif rule == "6 UMOD control":
            new.append("; ".join(["S6.umod.H4.A.503", "S6.umod.H4.B.503", "S6.umod.H4.A.489", "S6.umod.H4.B.489",
                                  "S6.umod.Q.aware_2nd.Q", "S6.umod.Q.aware_2nd.df", "S6.umod.Q.aware_2nd.p",
                                  "S6.umod.Q.p489.aware_2nd.p", "S6.umod.Q.boot.share"]))
        elif rule == "Rule 3(ii) (UMOD primary cell: all-signal GLS, 503)":
            new.append("S6.umod.Q.aware_2nd.p; S6.umod.Q.p489.aware_2nd.p; S6.umod.Q.boot.share")
        elif rule == "Rule 3(iv) (UMOD two-signal contrast)":
            m = re.match(r"^signal ([AB]):", inp)
            if not m: raise SystemExit("row %d: unrecognised Rule 3(iv) inputs" % i)
            X = m.group(1)
            new.append("S6.umod.z.%s.503; S6.umod.z.%s.489; S6.umod.H4.%s.503; S6.umod.H4.%s.489" % (X, X, X, X))
        elif rule == "6 UMOD secondary two-signal test; Rule 4b; Rule 5":
            new.append("; ".join(["S6.umod.two.503.kA", "S6.umod.two.503.kB", "S6.umod.two.503.p_aware_2nd", "S6.umod.two.503.p_ign",
                                  "S6.umod.two.503.p_aware_1st", "S6.umod.two.489.p_aware_2nd"]))
        elif rule == "Rule 3(i) (SuSiE with reference LD)":
            m = re.match(r"^([A-Za-z0-9]+_[A-Za-z0-9]+): s_protein ", inp)
            if not m: raise SystemExit("row %d: unrecognised Rule 3(i) inputs" % i)
            c = m.group(1)
            new.append("S6.d.%s.s_protein; S6.d.%s.s_outcome; S6.d.%s.s_protein_4fit; S6.d.%s.s_outcome_4fit" % (c, c, c, c))
        elif rule == "Rule 3(i) scope":
            new.append("S6.d.rule3i_scope")
        else:
            raise SystemExit("row %d: rule_ref not mapped: %r" % (i, rule))
    return new

ID_CELL = re.compile(r"^(?:06a:)?S6\.[A-Za-z0-9_.]+(?:; (?:06a:)?S6\.[A-Za-z0-9_.]+)*$")

def check_C1(cells, ids6a, ids6b):
    bad = []
    for c in cells:
        for tok in c.split("; "):
            pool, q = (ids6a, tok[4:]) if tok.startswith("06a:") else (ids6b, tok)
            if pool.count(q) != 1: bad.append(tok)
    return bad

def check_C2(old, new):
    if len(old) != len(new): return ["row count %d -> %d" % (len(old), len(new))]
    bad = []
    for i, (a, b) in enumerate(zip(old, new)):
        if len(a) != len(b): bad.append("row %d field count" % i); continue
        for j in range(len(a)):
            if j == 1 and i > 0: continue
            if a[j] != b[j]: bad.append("row %d col %d" % (i, j))
    return bad

def check_C3(cells):
    return [c for c in cells if not ID_CELL.match(c)]

def main():
    if md5f(RES6A) != RES6A_MD5: raise SystemExit("STOP: md5 of %s is not %s" % (RES6A, RES6A_MD5))
    raw = open(B6B, "rb").read()
    global PRE_MD5
    PRE_MD5 = md5b(raw)
    if all(ID_CELL.match(r[1]) for r in parse(raw)[1:]):
        raise SystemExit("STOP: branches.tsv inputs are already id lists (already normalised); nothing to do")
    say("pre-rewrite branches.tsv md5 %s (asserted)" % PRE_MD5)
    old = parse(raw)
    if old[0] != ["rule_ref", "inputs", "branch", "filled_template", "qualifiers"]: raise SystemExit("STOP: header")
    # C4 round trip: the writer must reproduce the original bytes
    rt = serialise(old)
    say("C4 writer round-trip on the original: %s" % ("PASS" if rt == raw else "FAIL"))
    if rt != raw: raise SystemExit("STOP: C4 failed")
    ids6a, ids6b = ids_of(RES6A), ids_of(RES6B)
    cells = mapping(old)
    new = [old[0]] + [[r[0], cells[i]] + r[2:] for i, r in enumerate(old[1:])]
    c1, c2, c3 = check_C1(cells, ids6a, ids6b), check_C2(old, new), check_C3(cells)
    say("C1 every id exists once in the results.tsv it names: %s %s" % ("PASS" if not c1 else "FAIL", c1))
    say("C2 only the inputs column changed, row count unchanged: %s %s" % ("PASS" if not c2 else "FAIL", c2))
    say("C3 inputs cells are id lists only: %s %s" % ("PASS" if not c3 else "FAIL", c3))
    if c1 or c2 or c3: raise SystemExit("STOP: a check failed; nothing written")
    # mutations: each must trip the check it names
    m1 = check_C1([cells[0] + "_MUT"] + cells[1:], ids6a, ids6b)
    newm = [r[:] for r in new]; newm[5][2] = newm[5][2] + " "
    m2 = check_C2(old, newm)
    m3 = check_C3([cells[0] + "; 0.5"] + cells[1:])
    m4 = (serialise([r[:] for r in old]).replace(b"\r\n", b"\n", 1) == raw)
    say("M1 (id suffixed _MUT) trips C1: %s" % ("TRIPPED" if m1 else "NOT TRIPPED"))
    say("M2 (a branch cell altered) trips C2: %s" % ("TRIPPED" if m2 else "NOT TRIPPED"))
    say("M3 (a number appended to an inputs cell) trips C3: %s" % ("TRIPPED" if m3 else "NOT TRIPPED"))
    say("M4 (one line ending altered) trips C4: %s" % ("TRIPPED" if not m4 else "NOT TRIPPED"))
    if not (m1 and m2 and m3 and not m4): raise SystemExit("STOP: a mutation did not trip; nothing written")
    # keep the original, then write
    os.makedirs(PREFIX_DIR, exist_ok=True)
    keep = os.path.join(PREFIX_DIR, "branches.tsv")
    if os.path.exists(keep):
        pass
    shutil.copyfile(B6B, keep)
    with open(os.path.join(PREFIX_DIR, "MD5SUMS.txt"), "w", newline="\n") as f:
        f.write("%s *./branches.tsv\n" % md5f(keep))
    say("original kept: %s (md5 %s)" % (keep, md5f(keep)))
    body = serialise(new)
    with open(B6B, "wb") as f: f.write(body)
    post = md5f(B6B)
    # re-read what was written and re-run every check on the file itself
    back = parse(open(B6B, "rb").read())
    ok = (back == new) and not check_C2(old, back) and not check_C3([r[1] for r in back[1:]]) and not check_C1([r[1] for r in back[1:]], ids6a, ids6b)
    say("post-rewrite re-read, C1-C3 on the written file: %s" % ("PASS" if ok else "FAIL"))
    say("post-rewrite branches.tsv md5 %s; rows %d (header + %d)" % (post, len(back), len(back) - 1))
    with open(LOG, "w", encoding="utf-8", newline="\n") as f: f.write("\n".join(out) + "\n")
    if not ok: raise SystemExit("FAIL after write: restore from %s" % keep)

if __name__ == "__main__":
    main()
