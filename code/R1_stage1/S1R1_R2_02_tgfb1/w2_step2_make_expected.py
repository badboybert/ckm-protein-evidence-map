# HONESTY NOTE (2026-10-10): These tables pin the md5 of what the pipeline emitted at step 2; they are REGRESSION guards (reject any later
# change) and do NOT establish that the pinned values are correct. Independent anchors: W6b s_max - A_third + L1/L2 within 2.1e-8; W2/W4b
# s_max, per-fit s and text cells - none.
"""w2_step2_make_expected.py -- ONE-OFF generator of step2_expected_changes.tsv (2026-10-10). NOT part of run_all.py: the table is the recorded,
reviewed expected-change list; regenerating it from a new output would defeat the gate. Records the post-step-2 output that FIXLOG.md describes
(9 fit-s rows + 3 s_max_4fit rows added; the three S2.susie.<trait>.s_max rows changed; branches row 7 inputs/qualifiers). Every cell needs a
reason or the script stops."""
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
sys.path.insert(0, HERE)
import w2_step2_gate_lib as G  # noqa: E402

PRE = HERE + "/_pre_step2_20261010"


def reason_res(k, col, kind):
    if kind == "added":
        if re.fullmatch(r"S2\.susie\.[a-z0-9_]+\.s_max_4fit", k):
            return "new row: the four-fit maximum kept as the labelled sensitivity"
        if re.fullmatch(r"S2\.susie\.[a-z0-9_]+\.(p2np_503panel|susie_489panel|p1fr_489panel)\.s", k):
            return "new row: s of one of the three added SuSiE fits (s-only, no credible set)"
        return None
    if re.fullmatch(r"S2\.susie\.[a-z0-9_]+\.s_max", k):
        return {"value": "s_max is now the maximum over all seven distinct fits (was four)", "estimand": "estimand text names the seven-fit set"}.get(col)
    return None


def reason_br(col):
    return {"inputs": "Rule 3(i) row lists the new seven-fit / sensitivity ids", "qualifiers": "Rule 3(i) row prints both maxima (seven-fit, four-fit)"}.get(col)


rows = []
for f in ("results.tsv", "branches.tsv"):
    nn, nr = G.rd(HERE + "/" + f)
    on, orr = G.rd(PRE + "/" + f)
    kn, ko = G.keys(f, nr), G.keys(f, orr)
    ixo = {k: i for i, k in enumerate(ko)}
    for i, k in enumerate(kn):
        if k not in ixo:
            r = reason_res(k, "*", "added") if f == "results.tsv" else None
            if r is None:
                sys.exit("unreasoned added row %s %s" % (f, k))
            rows.append({"file": f, "key": k, "column": "*", "new_md5": G.rowmd5(nn, nr[i]), "kind": "added", "reason": r})
            continue
        a, b = orr[ixo[k]], nr[i]
        for c in nn:
            if f == "results.tsv" and c == "quantity_id":
                continue
            if a[c] != b[c]:
                r = reason_res(k, c, "changed") if f == "results.tsv" else reason_br(c)
                if r is None:
                    sys.exit("unreasoned changed cell %s %s %s" % (f, k, c))
                rows.append({"file": f, "key": k, "column": c, "new_md5": G.md5(b[c]), "kind": "changed", "reason": r})
with open(HERE + "/step2_expected_changes.tsv", "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=G.COLS, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("step2_expected_changes.tsv: %d entries (%d added rows, %d changed cells)" % (len(rows), sum(r["kind"] == "added" for r in rows), sum(r["kind"] == "changed" for r in rows)))
