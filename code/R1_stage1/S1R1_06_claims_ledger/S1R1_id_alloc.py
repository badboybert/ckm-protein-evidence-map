#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Ledger id-space rules (Stage 5 step 1a, S2-D27, 2026-10-16). One implementation, imported by S1R1_build_ledger.py and ext_r2_stage1.py.

  * Ids are L + 3 digits (L000-L999). Every tool that parses ids reads three digits; a four-digit id is refused loudly (sorted_by_id), never sorted wrong.
  * RESERVED_TEST_IDS (L971-L984 and L999) belong to the self-tests of P31R_gate.py / P31R_figgate.py: no real row may take them, and allocate() never
    returns one. The set must equal P31R_gate.RESERVED_TEST_IDS (selftest).
  * allocate() hands out the next free numbers after the highest used one, skipping the reserved ids, up to the ceiling `hi` (L998 for ext_r2_stage1.py).
    Running out is an AssertionError, never a wrap or a reuse.
  * assert_ext_row_id (step 1d, S2-D28): every row an extension module adds has an id that is L + 3 digits (fullmatch) inside that module's range.
  * assert_no_overlap / assert_unowned_ranges: the extension-module registry (EXT_RANGES) has no overlap, and no row sits inside the range of a module that
    has not been created yet (a registered range is reserved for its module).
Run:  python -I -X utf8 S1R1_id_alloc.py --selftest      (exit 1 on any MISS)
"""
import os
import re
import sys

RESERVED_TEST_IDS = frozenset(["L%d" % n for n in range(971, 985)] + ["L999"])
_RESERVED_NUMS = frozenset(int(i[1:]) for i in RESERVED_TEST_IDS)
ID_RX = re.compile(r"L\d{3}")


def allocate(used_nums, n, lo, hi):
    """The next n free id numbers above max(used_nums | {lo - 1}), skipping RESERVED_TEST_IDS; AssertionError if one would exceed hi."""
    used = set(used_nums)
    clash = sorted(used & _RESERVED_NUMS)
    assert not clash, "id number(s) %s in use but reserved for the self-tests" % clash
    cur = max(used | {lo - 1}) + 1
    out = []
    for _ in range(n):
        while cur in _RESERVED_NUMS:
            cur += 1
        assert cur <= hi, "id space exhausted: %d more id(s) needed after L%d (ceiling L%d; L971-L984 and L999 are reserved for the self-tests)" % (n - len(out), out[-1] if out else max(used | {lo - 1}), hi)
        out.append(cur)
        cur += 1
    return out


def sorted_by_id(rows):
    """rows sorted by their id number; AssertionError (never a silent misorder) if an id is not L + 3 digits."""
    bad = [r["id"] for r in rows if not ID_RX.fullmatch(r["id"])]
    assert not bad, "id(s) %s are not L + 3 digits: the sort and the range checks read three digits (widening the id space is a tested migration, not an accident)" % bad[:5]
    return sorted(rows, key=lambda r: int(r["id"][1:]))


def assert_no_overlap(ranges):
    """ranges: {module: (lo, hi)}. No two registered ranges overlap."""
    rng = sorted(ranges.values())
    assert all(a[1] < b[0] for a, b in zip(rng, rng[1:])), "EXT_RANGES overlap: %s" % rng


def assert_unowned_ranges(ranges, rows):
    """Rows present before any extension module ran (`rows`) take no id inside a registered range or a reserved id; and for a registered module that
    does not exist on disk, no row at all lies in its range."""
    for r in rows:
        assert r["id"] not in RESERVED_TEST_IDS, "row %s takes an id reserved for the self-tests" % r["id"]
        if ID_RX.fullmatch(r["id"]):
            num = int(r["id"][1:])
            own = [m for m, (lo, hi) in ranges.items() if lo <= num <= hi]
            assert not own, "row %s sits inside the registered range of %s" % (r["id"], own)


def assert_unowned_after(ranges, rows, here):
    """After the extensions ran: a registered module not created yet has no row in its range (a registered range is reserved for that module)."""
    for m, (lo, hi) in ranges.items():
        if not os.path.exists(os.path.join(here, m)):
            sq = [r["id"] for r in rows if ID_RX.fullmatch(r["id"]) and lo <= int(r["id"][1:]) <= hi]
            assert not sq, "%s is registered (L%d-L%d) but not created, yet row(s) exist in its range: %s" % (m, lo, hi, sq[:5])


def assert_ext_row_id(module, rid, lo, hi):
    """Every row an extension module adds: an id that is L + 3 digits (fullmatch) and inside the module's registered range (Stage 5 step 1d, S2-D28). The earlier
    check read int(id[1:4]), which took 'L5000' and 'L500x' for L500 and let both through a strict build."""
    assert ID_RX.fullmatch(rid), "%s added %r: an extension row id must be L + 3 digits (widening the id space is a tested migration, not an accident)" % (module, rid)
    num = int(rid[1:])
    assert lo <= num <= hi, "%s added %s outside its range L%d-L%d" % (module, rid, lo, hi)


def assert_no_reserved(rows):
    """No real row takes an id reserved for the self-tests."""
    bad = sorted(r["id"] for r in rows if r["id"] in RESERVED_TEST_IDS)
    assert not bad, "row(s) %s take ids reserved for the self-tests (L971-L984, L999)" % bad


# ------------------------------------------------------------------------------------------------------------------------------------------------
def selftest():
    bad = 0

    def rep(lab, ok, detail=""):
        nonlocal bad
        print("selftest ids %-72s %s %s" % (lab, "OK" if ok else "MISS", detail))
        bad += not ok

    def trips(fn, frag):
        try:
            fn()
        except AssertionError as e:
            return frag in str(e)
        return False

    used = set(range(700, 971)) | set(range(985, 989))   # the shipped state: L700-L970 and L985-L988 used; L971-L984 are the self-test block
    got = allocate(used, 10, 700, 998)
    rep("allocate: 10 ids after L988 are L989-L998", got == list(range(989, 999)), "(%s..%s)" % (got[0], got[-1]))
    rep("allocate: an 11th id (one more than the free set) fails loudly", trips(lambda: allocate(used, 11, 700, 998), "id space exhausted"))
    rep("allocate: L999 is never returned even with a ceiling of 999", trips(lambda: allocate(used, 11, 700, 999), "id space exhausted") and 999 not in allocate(used, 10, 700, 999))
    low = set(range(700, 970))                         # a fixture whose free set crosses the reserved block
    got = allocate(low, 15, 700, 998)
    rep("allocate: skips L971-L984 (L970 then L985-L998)", got == [970] + list(range(985, 999)) and not (set(got) & _RESERVED_NUMS), "(%s)" % got[:3])
    rep("allocate: 16 ids from that fixture exceed the ceiling", trips(lambda: allocate(low, 16, 700, 998), "id space exhausted"))
    rep("allocate: a used reserved number is refused", trips(lambda: allocate(used | {975}, 1, 700, 998), "reserved for the self-tests"))
    rep("allocate: nothing used starts at the floor", allocate(set(), 2, 700, 998) == [700, 701])
    starts = [set(range(700, s)) for s in range(700, 972)] + [set(range(700, 971)) | set(range(985, s)) for s in range(985, 990)]
    rep("allocate: no result ever intersects the reserved block (every legal used set, 5 ids each)", all(not (set(allocate(u, 5, 700, 998)) & _RESERVED_NUMS) for u in starts), "(%d fixtures)" % len(starts))
    rows = [{"id": "L702"}, {"id": "L701"}, {"id": "L700"}]
    rep("sorted_by_id: three-digit ids sort", [r["id"] for r in sorted_by_id(rows)] == ["L700", "L701", "L702"])
    rep("sorted_by_id: L1000 fails loudly (it would sort as 100)", trips(lambda: sorted_by_id(rows + [{"id": "L1000"}]), "not L + 3 digits"))
    rep("sorted_by_id: a suffixed id fails loudly", trips(lambda: sorted_by_id(rows + [{"id": "L049a"}]), "not L + 3 digits"))
    rep("assert_no_overlap: disjoint ranges pass", not trips(lambda: assert_no_overlap({"a": (500, 599), "b": (610, 639), "c": (700, 998)}), "overlap"))
    rep("assert_no_overlap: an overlapping fixture entry trips the assert", trips(lambda: assert_no_overlap({"a": (500, 599), "b": (590, 609)}), "EXT_RANGES overlap"))
    ranges = {"ext_r2_stage5.py": (500, 599), "ext_r2_stage1.py": (700, 998)}
    rep("assert_unowned_ranges: a base row outside every range passes", not trips(lambda: assert_unowned_ranges(ranges, [{"id": "L250"}, {"id": "L049a"}]), "sits inside"))
    rep("assert_unowned_ranges: a base row inside the stage 5 range trips", trips(lambda: assert_unowned_ranges(ranges, [{"id": "L550"}]), "sits inside the registered range of ['ext_r2_stage5.py']"))
    rep("assert_unowned_ranges: a base row on a reserved id trips", trips(lambda: assert_unowned_ranges({}, [{"id": "L975"}]), "reserved for the self-tests"))
    here = os.path.dirname(os.path.abspath(__file__))
    rep("assert_unowned_after: a row in the range of a module not on disk trips", trips(lambda: assert_unowned_after({"ext_nonexistent_s5.py": (500, 599)}, [{"id": "L510"}], here), "registered (L500-L599) but not created"))
    rep("assert_unowned_after: no row in that range passes", not trips(lambda: assert_unowned_after({"ext_nonexistent_s5.py": (500, 599)}, [{"id": "L250"}], here), "registered"))
    rep("assert_unowned_after: a module that exists may own its range", not trips(lambda: assert_unowned_after({os.path.basename(__file__): (500, 599)}, [{"id": "L510"}], here), "registered"))
    rep("assert_no_reserved: L975 trips, L989 passes", trips(lambda: assert_no_reserved([{"id": "L975"}]), "reserved for the self-tests") and not trips(lambda: assert_no_reserved([{"id": "L989"}]), "reserved"))
    rep("assert_ext_row_id: L500 inside L500-L599 passes", not trips(lambda: assert_ext_row_id("ext_x.py", "L500", 500, 599), "ext_x.py"))
    rep("assert_ext_row_id: L5000 (four digits, read as L500 by id[1:4]) trips", trips(lambda: assert_ext_row_id("ext_x.py", "L5000", 500, 599), "must be L + 3 digits"))
    rep("assert_ext_row_id: L500x (a suffix, read as L500 by id[1:4]) trips", trips(lambda: assert_ext_row_id("ext_x.py", "L500x", 500, 599), "must be L + 3 digits"))
    rep("assert_ext_row_id: a suffixed id L049a, a short id L5 and a lower-case id l500 trip", all(trips(lambda i=i: assert_ext_row_id("ext_x.py", i, 0, 998), "must be L + 3 digits") for i in ("L049a", "L5", "l500")))
    rep("assert_ext_row_id: L600 and the reserved L975 outside L500-L599 trip the range", all(trips(lambda i=i: assert_ext_row_id("ext_x.py", i, 500, 599), "outside its range L500-L599") for i in ("L600", "L975")))
    # the set equals the gate's own copy
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(here)), "R2_rewrite", "scripts"))
    try:
        import P31R_gate as G
        rep("RESERVED_TEST_IDS equals P31R_gate.RESERVED_TEST_IDS", frozenset(G.RESERVED_TEST_IDS) == RESERVED_TEST_IDS, "(%d ids)" % len(RESERVED_TEST_IDS))
    except ImportError as e:
        rep("RESERVED_TEST_IDS equals P31R_gate.RESERVED_TEST_IDS", False, "(P31R_gate not importable: %s)" % e)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv[1:] else 0)
