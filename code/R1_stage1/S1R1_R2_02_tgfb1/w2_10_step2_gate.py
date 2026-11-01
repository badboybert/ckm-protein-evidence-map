"""w2_10_step2_gate.py -- run by run_all.py after w2_08_assemble.py (fix round 2026-10-10). ID-LIMITED, content-pinned gate: results.tsv /
branches.tsv may differ from the pre-step-2 copies in _pre_step2_20261010/ only in the cells and added rows listed in
step2_expected_changes.tsv, each with exactly the recorded content (see w2_step2_gate_lib.py). Any other difference prints
'PROBLEM <file> id=.. column=..' and the script exits 1, which stops run_all.py.
The gate is ACTIVE for every output folder (W2_OUT, as in w2_common.py), including from-scratch runs into another folder. ONE named exemption:
the single cell (S2.annot.rs1800470, notes) -- the md5 of the Ensembl REST JSON cache, which legitimately differs when the cache is re-fetched. Its
use is logged loudly ('EXEMPTION USED'); any other cell, including another column of the same row, still fails.
Test overrides (the mutation tests run on COPIES): W2_GATE_NEW (dir with results.tsv/branches.tsv), W2_GATE_EXPECTED (table)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
sys.path.insert(0, HERE)
import w2_step2_gate_lib as G  # noqa: E402


EXEMPT = ("results.tsv", "S2.annot.rs1800470", "notes")   # the ONE exempt (file, id, column)
EXEMPT_PREFIX = "%s id=%s column=%s: UNEXPECTED CHANGE" % EXEMPT


def main():
    new = os.environ.get("W2_GATE_NEW")
    exp = os.environ.get("W2_GATE_EXPECTED", HERE + "/step2_expected_changes.tsv")
    if new is None:
        import w2_common as W
        new = W.OUT
    P = G.run(new, HERE + "/_pre_step2_20261010", exp)
    used = [p for p in P if p.startswith(EXEMPT_PREFIX)]
    P = [p for p in P if p not in used]
    msg = "STEP2 GATE (id-limited, content-pinned; table %s): %d expected-change entries; %d problems" % (os.path.basename(exp), len(G.load_exp(exp)), len(P))
    emsg = ("*** EXEMPTION USED: %s id=%s column=%s differs from the pre-step-2 copy and was NOT checked (Ensembl REST cache md5 note); "
            "all other cells were checked ***" % EXEMPT) if used else None
    if os.environ.get("W2_GATE_NEW") is None:
        import w2_common as W
        W.log(msg)
        if emsg:
            W.log(emsg)
    else:
        print(msg)
        if emsg:
            print(emsg)
    for p in P:
        print("PROBLEM " + p, file=sys.stderr)
    return 1 if P else 0


if __name__ == "__main__":
    sys.exit(main())
