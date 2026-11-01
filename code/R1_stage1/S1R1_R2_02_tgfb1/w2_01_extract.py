"""w2_01_extract.py -- stream-cut the locus windows W2 needs from every large outcome file (never loaded whole;
brief section 0 rule 7). One pass per source file; resumable: a source whose stamp records the same md5 and windows
is skipped. Extracts are summary-statistic rows copied verbatim (header + matching lines) into OUT/extracts/.

Windows: GRCh38 gene +/- 1 Mb from the UKB-PPP instrument table (cis_start - 500 kb, cis_end + 500 kb; cis = gene
+/- 500 kb), widened to cover the whole UKB-PPP TGFB1 cis file. GRCh37 files: both bounds lifted hg38 -> hg19 with
REF_LIFTOVER_hg38ToHg19 (md5 asserted), then widened by 100 kb (bounds only; no row is joined on a lifted position here).
"""
import csv
import gzip
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

LOCI = ["TGFB1", "CELSR2", "IL6R", "ABO", "GALNT3", "SPON1", "FGF5", "FGFR4", "ARG1"]
# source id -> (build, loci, chrom_idx, pos_idx, sep(None = whitespace), n_header_skip_prefix)
SOURCES = {
    "OUT_FG_R10_I9_CHD": ("38", ["TGFB1", "CELSR2", "IL6R"], 0, 1, b"\t", b"#chrom"),
    "OUT_FG_R10_T2D": ("38", ["ABO", "GALNT3", "ARG1", "FGFR4"], 0, 1, b"\t", b"#chrom"),
    "OUT_FG_R10_I9_AF": ("38", ["SPON1"], 0, 1, b"\t", b"#chrom"),
    "OUT_CARDIOGRAMplusC4D_ieu-a-7": ("37", ["TGFB1"], 0, 1, b"\t", b"#CHROM"),
    "OUT_DIAMANTE2018_noUKB_EUR_chrpos": ("37", ["FGFR4"], 1, 2, b"\t", b"SNP\t"),
    "OUT_CKDGEN2019_eGFR_EA": ("37", ["FGF5"], 0, 1, None, b"Chr "),
    "EAS_BBJ_CAD_hum0197v5": ("37", ["TGFB1", "CELSR2", "IL6R"], 1, 2, b"\t", b"v\t"),
    "EAS_AGEN_T2D_Spracklen2020": ("37", ["ABO", "GALNT3", "ARG1", "FGFR4"], 2, 3, b"\t", b"MarkerName"),
    "EAS_BBJ_AF_2022": ("37", ["SPON1"], 0, 1, b"\t", b"CHR\t"),
    "EAS_TWB2_BBJ_eGFR": ("37", ["FGF5"], 8, 9, None, b"MarkerName"),
    "EAS_BBJ_T2D_hum0197v3": ("37", ["FGFR4"], 1, 2, b"\t", b"v\t"),
}
MARGIN37 = 100000


def windows38():
    ins = list(csv.DictReader(open(W.res("EXP_UKBPPP_INSTRUMENT_TABLE"), encoding="utf-8")))
    w = {}
    for g in LOCI:
        r = [x for x in ins if x["gene"] == g]
        assert len(r) == 1, (g, len(r))
        r = r[0]
        w[g] = (r["gene_chrom"], int(r["cis_start"]) - 500000, int(r["cis_end"]) + 500000)
    # TGFB1: cover the whole UKB-PPP cis file too
    c, lo, hi = w["TGFB1"]
    w["TGFB1"] = (c, min(lo, 40319325 - 1000), max(hi, 42319263 + 1000))
    return w


def lift_bounds(w38):
    from pyliftover import LiftOver
    lo = LiftOver(W.res("REF_LIFTOVER_hg38ToHg19"))
    out = {}
    for g, (c, a, b) in w38.items():
        def lift(p, step):
            for k in range(0, 200000, 1000):
                r = lo.convert_coordinate("chr" + c, p + step * k - 1)
                if r and r[0][0] == "chr" + c:
                    return r[0][1] + 1
            raise SystemExit("cannot lift bound %s %s" % (g, p))
        a37, b37 = lift(a, 1), lift(b, -1)
        if a37 > b37:
            a37, b37 = b37, a37
        out[g] = (c, a37 - MARGIN37, b37 + MARGIN37)
    return out


def open_any(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")


def cut(rid, spec, w38, w37):
    build, loci, ci, pi, sep, hdr = spec
    wins = {g: (w38[g] if build == "38" else w37[g]) for g in loci}
    bychr = {}
    for g, (c, a, b) in wins.items():
        bychr.setdefault(c.encode(), []).append((g, a, b))
    src = W.res(rid)
    srcmd5 = W.md5_file(src)
    stamp_p = "%s/extracts/%s.stamp.json" % (W.OUT, rid)
    if os.path.exists(stamp_p):
        st = json.load(open(stamp_p, encoding="utf-8"))
        if st["source_md5"] == srcmd5 and st["windows"] == {g: list(v) for g, v in wins.items()} and \
                all(os.path.exists(W.OUT + "/extracts/" + f) for f in st["files"]):
            W.log("[extract] %s up to date (skip)" % rid)
            return st
    t0 = time.time()
    outs = {g: open("%s/extracts/%s__%s.tsv.part" % (W.OUT, rid, g), "wb") for g in loci}
    nkeep = {g: 0 for g in loci}
    header = None
    n = 0
    nsci = [0, 0, 0]       # [sci-notation coordinates on window chromosomes, unresolvable AND possibly in a window,
    #                         unresolvable but provably outside every window (rounding interval misses them all)]
    maxsplit = max(ci, pi) + 1
    with open_any(src) as fh:
        for line in fh:
            if header is None:
                if line.startswith(b"##"):
                    continue
                if line.startswith(hdr):
                    header = line
                    for f in outs.values():
                        f.write(line)
                    continue
                raise SystemExit("%s: unexpected first line %r" % (rid, line[:80]))
            n += 1
            if n % 5000000 == 0:
                W.log("[extract] %s %d rows read, kept %s, %.0fs" % (rid, n, nkeep, time.time() - t0))
            f = line.split(sep, maxsplit) if sep else line.split(None, maxsplit)
            c = f[ci]
            if c.startswith(b"chr"):
                c = c[3:]
            ws = bychr.get(c)
            if not ws:
                continue
            try:
                p = int(f[pi])
            except ValueError:
                # scientific-notation coordinate (hygiene section 7): resolve from the text coordinate embedded in the
                # first field (BBJ 'v' = chr:pos:ref:alt; AGEN/TWB2 MarkerName = chr:pos_...) and check the float agrees
                pf = float(f[pi])
                mant, ex = f[pi].lower().split(b"e")
                ndig = len(mant.replace(b".", b"").replace(b"-", b"").lstrip(b"0")) or 1
                tol = 0.5 * 10 ** (int(ex) - (ndig - 1))
                txt = f[0].split(b":")
                pt = int(txt[1].split(b"_")[0]) if len(txt) > 1 and txt[1].split(b"_")[0].isdigit() else None
                nsci[0] += 1
                if pt is None or abs(pt - pf) > tol:
                    # no usable text coordinate: provably outside every window if the rounding interval misses them all
                    if any(a <= pf + tol and pf - tol <= b for _, a, b in ws):
                        nsci[1] += 1
                        W.log("ASSERT_FAILED %s: unresolvable coordinate %r may fall inside a window" % (rid, line[:80]))
                    else:
                        nsci[2] += 1
                    continue
                p = pt
            for g, a, b in ws:
                if a <= p <= b:
                    outs[g].write(line)
                    nkeep[g] += 1
    files = []
    for g, f in outs.items():
        f.close()
        fn = "%s__%s.tsv" % (rid, g)
        os.replace(W.OUT + "/extracts/" + fn + ".part", W.OUT + "/extracts/" + fn)
        files.append(fn)
        if nkeep[g] == 0:
            W.log("ASSERT_FAILED extract %s %s kept 0 rows" % (rid, g))
    st = {"source": rid, "path": src, "source_md5": srcmd5, "source_bytes": os.path.getsize(src), "build": build,
          "rows_read": n, "rows_kept": nkeep, "sci_notation_coords_on_window_chroms": nsci[0],
          "sci_notation_unresolved_possibly_in_window": nsci[1], "sci_notation_unresolved_provably_outside": nsci[2], "windows": {g: list(v) for g, v in wins.items()}, "files": files,
          "file_md5": {fn: W.md5_file(W.OUT + "/extracts/" + fn) for fn in files},
          "seconds": round(time.time() - t0, 1), "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
    with open(stamp_p, "w", encoding="utf-8") as fh:
        json.dump(st, fh, indent=1)
    W.log("[extract] %s done: %d rows read, kept %s, %.0fs; scientific-notation coordinates on window chromosomes %d "
          "(resolved from the text coordinate %d; no text coordinate but provably outside every window %d; possibly inside %d)" % (
              rid, n, nkeep, time.time() - t0, nsci[0], nsci[0] - nsci[1] - nsci[2], nsci[2], nsci[1]))
    return st


def main():
    os.makedirs(W.OUT + "/extracts", exist_ok=True)
    W.log("==== w2_01_extract ====")
    w38 = windows38()
    w37 = lift_bounds(w38)
    for g in LOCI:
        W.log("[extract] window %s GRCh38 %s ; GRCh37 (lifted bounds + %d) %s" % (g, w38[g], MARGIN37, w37[g]))
    with open(W.OUT + "/extracts/windows.json", "w", encoding="utf-8") as fh:
        json.dump({"grch38": w38, "grch37": w37}, fh, indent=1)
    todo = sys.argv[1:] or list(SOURCES)
    for rid in todo:
        cut(rid, SOURCES[rid], w38, w37)


if __name__ == "__main__":
    main()
