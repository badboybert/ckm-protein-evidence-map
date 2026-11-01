"""w2_step2_gate_lib.py -- ID-LIMITED, content-pinned expected-change gate for W2 (fix round 2026-10-10; same algorithm as the R library
step2_gate_lib.R in W6b and W4b). W2 had NO equivalent check before: the Stage 2 step 2 outputs were compared with the pre-step-2 copies by hand
only, so a wrong value in any cell (an unrelated one, a per-fit s, an s_max) left run_all.py exiting 0.

A regenerated results.tsv / branches.tsv may differ from the baseline (_pre_step2_20261010/) ONLY in the cells and added rows listed in
step2_expected_changes.tsv (columns file, key, column, new_md5, kind, reason), and each listed cell must carry exactly the recorded content
(md5 of the cell text; md5 of the tab-joined row for an added row). Anything else is a PROBLEM line naming the id and the column.
key = quantity_id (results.tsv) or 'row N (rule_ref)' (branches.tsv)."""
import csv
import hashlib

COLS = ["file", "key", "column", "new_md5", "kind", "reason"]


def md5(s):
    return hashlib.md5(s.encode("utf-8")).hexdigest()


def rd(path):
    with open(path, encoding="utf-8", newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        names = next(r)
        return names, [dict(zip(names, row)) for row in r]


def keys(fname, rows):
    return [x["quantity_id"] for x in rows] if fname == "results.tsv" else ["row %d (%s)" % (i + 1, x["rule_ref"]) for i, x in enumerate(rows)]


def rowmd5(names, row):
    return md5("\t".join(row[c] for c in names))


def clip(x):
    x = x.replace("\r", " ").replace("\n", " ").replace("\t", " ")
    return x if len(x) <= 70 else x[:70] + "..."


def check(fname, new, old, exp):
    nn, nr = new
    on, orr = old
    if nn != on:
        return ["%s: column names differ from the baseline" % fname]
    P = []
    kn, ko = keys(fname, nr), keys(fname, orr)
    if len(set(kn)) != len(kn):
        P.append("%s: duplicate key in the new file" % fname)
    if len(set(ko)) != len(ko):
        P.append("%s: duplicate key in the baseline" % fname)
    E = {(e["key"], e["column"]): e["new_md5"] for e in exp if e["file"] == fname}
    seen = set()
    for k in ko:
        if k not in set(kn):
            P.append("%s id=%s column=*: row REMOVED" % (fname, k))
    ixn, ixo = {k: i for i, k in enumerate(kn)}, {k: i for i, k in enumerate(ko)}
    for k in kn:
        if k in ixo:
            continue
        seen.add((k, "*"))
        got = rowmd5(nn, nr[ixn[k]])
        if (k, "*") not in E:
            P.append("%s id=%s column=*: UNEXPECTED NEW ROW (not in the expected-changes table)" % (fname, k))
        elif E[(k, "*")] != got:
            P.append("%s id=%s column=*: new row differs from the recorded expected row (md5 %s, recorded %s)" % (fname, k, got, E[(k, "*")]))
    for k in kn:
        if k not in ixo:
            continue
        a, b = orr[ixo[k]], nr[ixn[k]]
        for c in nn:
            if fname == "results.tsv" and c == "quantity_id":
                continue
            if a[c] == b[c]:
                continue
            seen.add((k, c))
            if (k, c) not in E:
                P.append("%s id=%s column=%s: UNEXPECTED CHANGE (old '%s' -> new '%s')" % (fname, k, c, clip(a[c]), clip(b[c])))
            elif E[(k, c)] != md5(b[c]):
                P.append("%s id=%s column=%s: changed cell differs from the recorded expected content (new '%s', md5 %s, recorded %s)" % (
                    fname, k, c, clip(b[c]), md5(b[c]), E[(k, c)]))
    for (k, c) in E:
        if (k, c) not in seen:
            P.append("%s id=%s column=%s: EXPECTED CHANGE ABSENT (listed in the table, not present in the output)" % (fname, k, c))
    return P


def load_exp(path):
    names, rows = rd(path)
    assert names == COLS, names
    assert len({(r["file"], r["key"], r["column"]) for r in rows}) == len(rows), "duplicate table entry"
    assert all(r["reason"].strip() for r in rows), "an entry without a reason"
    return rows


def run(new_dir, pre_dir, exp_path):
    exp = load_exp(exp_path)
    P = []
    for f in ("results.tsv", "branches.tsv"):
        P += check(f, rd(new_dir + "/" + f), rd(pre_dir + "/" + f), exp)
    return P
