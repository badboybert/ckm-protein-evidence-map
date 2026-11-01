"""w2_rerun_diff.py -- compare results.tsv / branches.tsv of this folder with a from-scratch rerun (rerun_check/).
Numeric fields compared at relative tolerance 1e-9; text fields exactly (notes that embed cache-file md5s are listed)."""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402


def load(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def same(a, b):
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    return abs(x - y) <= 1e-9 * max(1.0, abs(x), abs(y))


def main():
    other = sys.argv[1] if len(sys.argv) > 1 else W.OUT + "/rerun_check"
    a = {r["quantity_id"]: r for r in load(W.OUT + "/results.tsv")}
    b = {r["quantity_id"]: r for r in load(other + "/results.tsv")}
    diffs = []
    for q in sorted(set(a) | set(b)):
        if q not in a or q not in b:
            diffs.append((q, "missing in " + ("primary" if q not in a else "rerun"), ""))
            continue
        for c in a[q]:
            if not same(a[q][c], b[q][c]):
                diffs.append((q, c, "%s | %s" % (a[q][c][:120], b[q][c][:120])))
    ba, bb = load(W.OUT + "/branches.tsv"), load(other + "/branches.tsv")
    bdiff = [(x["rule_ref"], k) for x, y in zip(ba, bb) for k in x if x[k] != y[k]] + ([("n_rows", len(ba), len(bb))] if len(ba) != len(bb) else [])
    W.log("RERUN DIFF vs %s: %d results rows primary / %d rerun; %d differing fields; branches differing fields %d" % (
        other, len(a), len(b), len(diffs), len(bdiff)))
    for d in diffs:
        W.log("  rerun diff: %s [%s] %s" % d)
    for d in bdiff:
        W.log("  rerun branch diff: %s" % (d,))


if __name__ == "__main__":
    main()
