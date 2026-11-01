"""Regenerate sections_R2/P31R_Legend_AF3_R2.md from the per-figure AF3 legend sources (registry order).

Sources: the 'legend' of every kind=="af3" figure in P31R_ACTIVE_FIGURES.json (paths relative to the V3 root).
Output is LF; each source is read as UTF-8 with CRLF normalised. Identity is checked on the copy AS WRITTEN
TO DISK: after writing, and in --check, the file is re-read and every '## Figure S<n>' section is compared
with its source re-read from disk. --check is the gate; run it in the step 1c gate suite.
Usage: python -I -X utf8 build_legend_AF3_R2.py [--out PATH] [--check]
  --check: assert the copy at --out (default: the operative copy) is current and section-identical.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(V3, "R2_rewrite", "sections_R2", "P31R_Legend_AF3_R2.md")
HEADER = ("# Additional file 3: figure legends (R2). Generated from figures/af3/legends/P31R_Legend_FigS1-S6.md, "
          "which are the files the gate reads; do not edit this copy.")
SECTION = re.compile(r"(?m)^(?=## Figure S\d)")


def rd(p):
    return open(p, "rb").read().decode("utf-8").replace("\r\n", "\n").strip()


def af3_figures():
    reg = json.load(open(os.path.join(HERE, "P31R_ACTIVE_FIGURES.json"), encoding="utf-8"))
    figs = sorted((f for f in reg["figures"] if f["kind"] == "af3"), key=lambda f: f["order"])
    assert figs, "no af3 figures in registry"
    return figs


def build(figs):
    secs = []
    for f in figs:
        s = rd(os.path.join(V3, f["legend"]))
        assert s.startswith("## %s" % f["id"].replace("Figure", "Figure ", 1)), (f["id"], s[:30])
        secs.append(s)
    return HEADER + "\n\n" + "\n\n".join(secs) + "\n"


def verify_on_disk(path, figs):
    """Re-read the copy from disk; each section must equal its source re-read from disk."""
    cur = open(path, "rb").read().decode("utf-8").replace("\r\n", "\n")
    parts = SECTION.split(cur)[1:]
    assert len(parts) == len(figs), ("section count", len(parts), len(figs))
    for p, f in zip(parts, figs):
        src = rd(os.path.join(V3, f["legend"]))
        assert p.strip() == src, "copy section differs from its source: " + f["id"]
    return len(parts)


if __name__ == "__main__":
    figs = af3_figures()
    text = build(figs)
    out = OUT
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    if "--check" in sys.argv:
        cur = open(out, "rb").read().decode("utf-8")
        assert cur == text, "AF3_R2 copy is stale: rerun build_legend_AF3_R2.py"
        n = verify_on_disk(out, figs)
        print("AF3_R2 copy current;", n, "sections equal their sources")
    else:
        open(out, "wb").write(text.encode("utf-8"))
        n = verify_on_disk(out, figs)
        print("wrote", out, len(text), "chars,", n, "sections, each equal to its source")
