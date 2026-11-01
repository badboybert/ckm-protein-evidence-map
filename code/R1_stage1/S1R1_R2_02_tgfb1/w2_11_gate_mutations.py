"""w2_11_gate_mutations.py -- mutation tests of the id-limited gate w2_10_step2_gate.py (fix round 2026-10-10). Every mutation is applied to a
COPY of results.tsv / branches.tsv (and of the expected-changes table for the legitimate-change case) in a temp folder; the real gate runs on the
copy as a separate process. A mutation is TRIPPED only if the exit status is non-zero AND the output names the mutated id and column (and the
kind of problem). Controls (unmutated; a legitimate change recorded in the table) are TRIPPED when they exit 0 as required.
Writes mutation_tests.tsv; exits 1 if any case does not behave."""
import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
sys.path.insert(0, HERE)
import w2_step2_gate_lib as G  # noqa: E402

GATE = HERE + "/w2_10_step2_gate.py"
EXPT = HERE + "/step2_expected_changes.tsv"
PRE = HERE + "/_pre_step2_20261010"
OUTD = os.environ.get("W2_OUT", HERE).replace("\\", "/")   # output folder under test (as w2_common.OUT): mutation copies are made from it, mutation_tests.tsv is written to it
CK = "w2_10_step2_gate.py id-limited gate (exit non-zero, message names id and column)"
OUT = []


def write(path, names, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=names, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def gate_copy(mut=None):
    """mut(R, B, E) mutates the row lists in place (R results, B branches, E expected table rows)."""
    rn, R = G.rd(OUTD + "/results.tsv")
    bn, B = G.rd(OUTD + "/branches.tsv")
    E = G.load_exp(EXPT)
    if mut:
        mut(R, B, E)
    td = tempfile.mkdtemp(prefix="w2gate_")
    write(td + "/results.tsv", rn, R)
    write(td + "/branches.tsv", bn, B)
    write(td + "/exp.tsv", G.COLS, E)
    env = dict(os.environ, W2_GATE_NEW=td, W2_GATE_EXPECTED=td + "/exp.tsv")
    r = subprocess.run([sys.executable, "-X", "utf8", GATE], env=env, capture_output=True, text=True, encoding="utf-8")
    return r.returncode, r.stdout + r.stderr


def one(rows, qid):
    hits = [r for r in rows if r["quantity_id"] == qid]
    assert len(hits) == 1, (qid, len(hits))          # a mutation target must be unique
    return hits[0]


def plus(x):
    return "%.10g" % (float(x) + 1e-3)


def rec(mid, check, mutation, ok, detail):
    OUT.append({"mutation_id": mid, "check_named": check, "mutation": mutation, "result": "TRIPPED" if ok else "NOT TRIPPED", "detail": detail})
    print("MUTATION %s -> %s (%s)" % (mid, "TRIPPED" if ok else "NOT TRIPPED", detail), flush=True)


def case(mid, mutation, mut, must_have):
    rc, out = gate_copy(mut)
    ok = rc != 0 and all(m in out for m in must_have)
    rec(mid, CK, mutation, ok, "exit %d; message must contain: %s" % (rc, " | ".join(must_have)))


def main():
    rc, out = gate_copy()
    ctl = G.run(OUTD, PRE, EXPT)
    rec("G0", "w2_10 gate (control: unmutated output must pass)", "none (copy of the real results.tsv/branches.tsv/table)", rc == 0 and not ctl and "0 problems" in out,
        "exit %d; %d problems (a control is TRIPPED when it passes as required)" % (rc, len(ctl)))

    def m1(R, B, E):
        one(R, "S2.r.503")["value"] = plus(one(R, "S2.r.503")["value"])
    case("G1", "+1e-3 on the unrelated cell S2.r.503 value", m1, ["id=S2.r.503 column=value", "UNEXPECTED CHANGE"])

    def m2a(R, B, E):
        x = one(R, "S2.susie.protein.s_max_4fit")                       # a sensitivity row is pinned as a whole row
        x["value"] = plus(x["value"])
    case("G2a", "+1e-3 on the four-fit sensitivity row S2.susie.protein.s_max_4fit (value)", m2a, ["id=S2.susie.protein.s_max_4fit column=*", "differs from the recorded expected row"])

    def m2b(R, B, E):
        x = one(R, "S2.susie.protein.s_max")
        x["value"] = plus(x["value"])
    case("G2b", "+1e-3 on a listed s_max (S2.susie.protein.s_max value): listed id, wrong content", m2b, ["id=S2.susie.protein.s_max column=value", "differs from the recorded expected content"])

    def m2c(R, B, E):
        x = one(R, "S2.susie.protein.s_max")
        x["se"] = "0.001"
    case("G2c", "an s_max row changed in a column not on the list (S2.susie.protein.s_max se := 0.001)", m2c, ["id=S2.susie.protein.s_max column=se", "UNEXPECTED CHANGE"])

    def m3(R, B, E):
        x = one(R, "S2.susie.protein.susie_489panel.s")
        x["value"] = plus(x["value"])
    case("G3", "+1e-3 on a new-fit s (S2.susie.protein.susie_489panel.s value)", m3, ["id=S2.susie.protein.susie_489panel.s column=*", "differs from the recorded expected row"])

    def m3b(R, B, E):
        ids = [r["quantity_id"] for r in R if r["quantity_id"].startswith("S2.susie.protein.") and r["quantity_id"].endswith(".s")]
        old = [i for i in ids if "panel" not in i]
        assert old, ids
        x = one(R, old[0])
        x["value"] = plus(x["value"])
        m3b.target = old[0]
    rc, out = gate_copy(m3b)
    t = m3b.target
    rec("G3b", CK, "+1e-3 on a pre-existing per-fit s (%s value)" % t, rc != 0 and ("id=%s column=value" % t) in out and "UNEXPECTED CHANGE" in out,
        "exit %d; message must contain: id=%s column=value | UNEXPECTED CHANGE" % (rc, t))

    def m5(R, B, E):
        x = dict(R[0])
        x["quantity_id"] = "S2.zz.unlisted_new_row"
        R.append(x)
    case("G5", "an unlisted new row appended", m5, ["id=S2.zz.unlisted_new_row column=*", "UNEXPECTED NEW ROW"])

    def m6(R, B, E):
        assert not any(e["file"] == "branches.tsv" and e["key"].startswith("row 3 (") for e in E)
        B[2]["qualifiers"] += " x"
    case("G6", "an unrelated branches.tsv cell changed (row 3, qualifiers)", m6, ["branches.tsv id=row 3", "column=qualifiers", "UNEXPECTED CHANGE"])

    def m7(R, B, E):
        x = one(R, "S2.r.503")
        R.remove(x)
    case("G7", "a results.tsv row removed (S2.r.503)", m7, ["id=S2.r.503 column=*", "row REMOVED"])

    def m8(R, B, E):
        on, orr = G.rd(PRE + "/results.tsv")
        one(R, "S2.susie.protein.s_max")["value"] = one(orr, "S2.susie.protein.s_max")["value"]
    case("G8", "a listed change reverted to the pre-step-2 value (S2.susie.protein.s_max value): EXPECTED CHANGE ABSENT", m8, ["id=S2.susie.protein.s_max column=value", "EXPECTED CHANGE ABSENT"])

    def legit(upd):
        def f(R, B, E):
            x = one(R, "S2.susie.protein.s_max")
            x["estimand"] += " [legit]"
            j = [e for e in E if e["file"] == "results.tsv" and e["key"] == "S2.susie.protein.s_max" and e["column"] == "estimand"]
            assert len(j) == 1
            if upd:
                j[0]["new_md5"] = G.md5(x["estimand"])
        return f
    rc, out = gate_copy(legit(True))
    rec("G4", "w2_10 gate (control: a legitimate allowed change recorded in the table must pass)", "S2.susie.protein.s_max estimand text changed AND its table md5 updated", rc == 0, "exit %d (must be 0)" % rc)
    case("G4b", "the same legitimate-looking change WITHOUT updating the table", legit(False), ["id=S2.susie.protein.s_max column=estimand", "differs from the recorded expected content"])

    def e1(R, B, E):
        one(R, "S2.annot.rs1800470")["notes"] = "variation + vep JSON md5 " + "0" * 32 + " / " + "1" * 32      # re-fetched cache: different md5s
    rc, out = gate_copy(e1)
    rec("E1", "w2_10 gate (control: the ONE named exemption cell may differ, logged loudly)", "S2.annot.rs1800470 notes replaced by other md5s", rc == 0 and "EXEMPTION USED" in out and "S2.annot.rs1800470" in out,
        "exit %d (must be 0) and the output must contain EXEMPTION USED" % rc)

    def e2(R, B, E):
        x = one(R, "S2.annot.rs1800470")
        x["value"] = x["value"] + " x"
    case("E2", "another column of the exempt row changed (S2.annot.rs1800470 value + ' x')", e2, ["id=S2.annot.rs1800470 column=value", "UNEXPECTED CHANGE"])

    def e3(R, B, E):
        t = [r for r in R if r["quantity_id"].startswith("S2.annot.") and r["quantity_id"] != "S2.annot.rs1800470" and r["notes"]]
        assert t
        t[0]["notes"] += " x"
        e3.target = t[0]["quantity_id"]
    rc, out = gate_copy(e3)
    rec("E3", CK, "the notes cell of a different annotation row changed (%s)" % e3.target, rc != 0 and ("id=%s column=notes" % e3.target) in out and "UNEXPECTED CHANGE" in out,
        "exit %d; message must contain: id=%s column=notes | UNEXPECTED CHANGE" % (rc, e3.target))

    def e4(R, B, E):
        one(R, "S2.annot.rs1800470")["notes"] = "changed"
        one(R, "S2.r.503")["value"] = plus(one(R, "S2.r.503")["value"])
    case("E4", "exempt cell changed AND an unrelated cell changed (exemption must not mask S2.r.503 value)", e4, ["id=S2.r.503 column=value", "UNEXPECTED CHANGE", "EXEMPTION USED"])

    write(OUTD + "/mutation_tests.tsv", ["mutation_id", "check_named", "mutation", "result", "detail"], OUT)
    n = sum(o["result"] == "TRIPPED" for o in OUT)
    print("W2 gate mutation tests: %d of %d tripped (G0 and G4 are controls: TRIPPED = passed as required)" % (n, len(OUT)))
    return 0 if n == len(OUT) else 1


if __name__ == "__main__":
    sys.exit(main())
