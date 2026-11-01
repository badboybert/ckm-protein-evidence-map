"""Render every active figure unit once, resumably (Stage 5 step 1c; plan: one render of all 13 units after the ledger batch).

Reads the active figure list (af3/P31R_ACTIVE_FIGURES.json), runs each generator once with Rscript, sequentially, and logs each run to
review_R2/render_1c_logs/<generator>.log. A campaign file (render_log/_campaign.json, written on the first run) fixes the start time and the
ledger md5: a unit counts as done when its render_log row names the current generator md5 and was rendered after the campaign start, so a
re-invocation after a crash re-runs only the generators whose units are not done. Prints one flushed progress line per generator.

    python -I -X utf8 render_all_units.py            # start or resume the campaign
    python -I -X utf8 render_all_units.py --status   # liveness / progress only (renders nothing)
    python -I -X utf8 render_all_units.py --new      # start a new campaign (only after a further ledger change)
"""
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

FIG = Path(__file__).resolve().parent
R2 = FIG.parent
V3 = R2.parent
LEDGER = V3 / "R1_stage1" / "S1R1_06_claims_ledger" / "S1R1_CLAIMS_LEDGER.csv"
CAMPAIGN = FIG / "render_log" / "_campaign.json"
LOGS = R2 / "review_R2" / "render_1c_logs"
HEART = LOGS / "_heartbeat.json"


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def units_by_generator():
    act = json.loads((FIG / "af3" / "P31R_ACTIVE_FIGURES.json").read_text(encoding="utf-8"))
    out = {}
    for f in act["figures"]:
        ids = [f["id"]] + ["%s_p%d" % (f["id"], k) for k in range(2, int(f.get("pages", 1)) + 1)]
        out.setdefault(f["generator"], []).extend(ids)
    return out


def unit_done(unit, gen_md5, start):
    p = FIG / "render_log" / ("%s.tsv" % unit)
    if not p.is_file():
        return False
    rows = list(csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"))
    return bool(rows) and rows[-1]["generator_md5"] == gen_md5 and rows[-1]["rendered"] >= start


def main(argv):
    LOGS.mkdir(parents=True, exist_ok=True)
    if "--new" in argv and CAMPAIGN.exists():
        CAMPAIGN.unlink()
    if not CAMPAIGN.exists():
        if "--status" in argv:
            print("no campaign started"); return 0
        CAMPAIGN.write_text(json.dumps({"start": time.strftime("%Y-%m-%d %H:%M:%S"), "ledger_md5": md5(LEDGER)}, indent=1), encoding="utf-8")
    camp = json.loads(CAMPAIGN.read_text(encoding="utf-8"))
    if md5(LEDGER) != camp["ledger_md5"]:
        print("REFUSED: the ledger changed since the campaign started (%s -> %s); start a new campaign with --new" % (camp["ledger_md5"][:8], md5(LEDGER)[:8]))
        return 2
    gens = units_by_generator()
    todo = []
    for g, units in gens.items():
        gm = md5(V3 / g)
        done = [u for u in units if unit_done(u, gm, camp["start"])]
        if len(done) < len(units):
            todo.append((g, units))
    n_units = sum(len(u) for u in gens.values())
    print("campaign %s (ledger %s): %d generators, %d units; %d generator(s) still to render" % (camp["start"], camp["ledger_md5"][:8], len(gens), n_units, len(todo)), flush=True)
    if "--status" in argv:
        if HEART.exists():
            print("heartbeat:", HEART.read_text(encoding="utf-8"))
        for g, units in todo:
            print("  to do:", g, units)
        return 0
    for i, (g, units) in enumerate(todo, 1):
        log = LOGS / (Path(g).stem + ".log")
        HEART.write_text(json.dumps({"generator": g, "pid": os.getpid(), "since": time.strftime("%Y-%m-%d %H:%M:%S")}), encoding="utf-8")
        t0 = time.time()
        with log.open("w", encoding="utf-8") as fh:
            rc = subprocess.run(["Rscript", str(V3 / g)], cwd=str(FIG), stdout=fh, stderr=subprocess.STDOUT, timeout=3600).returncode
        gm = md5(V3 / g)
        ok = rc == 0 and all(unit_done(u, gm, camp["start"]) for u in units)
        print("[%d/%d] %s -> %s in %.0f s (exit %d; log %s)" % (i, len(todo), g, "OK " + ",".join(units) if ok else "FAILED", time.time() - t0, rc, log.name), flush=True)
        if not ok:
            print("STOPPED at the first failure; fix it, then re-run this script to resume.", flush=True)
            return 1
    HEART.unlink(missing_ok=True)
    print("ALL %d units rendered in this campaign." % n_units, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
