#!/usr/bin/env python
"""W3 step 0 (pre-declaration §3): governance md5 asserts, exposure build check, outcome stream-cut.

Reads every input path from RESOURCES.json by id and asserts its md5 at read time.
- Exposures (UKB-PPP EUR cis windows): parse the GRCh37 chr:pos:ref:alt embedded in ID, check ID ref/alt
  against ALLELE0/ALLELE1, lift the GRCh37 position to GRCh38 with REF_LIFTOVER_hg19ToHg38 and compare it
  with the file's GRCh38 `pos` column; the GRCh37 window bounds are min/max of the ID positions.
- Outcomes: DIAMANTE 2018 noUKB (plain text, GRCh37, 1.1 GB) and ieu-a-7 (BGZF VCF, GRCh37, no .tbi) are
  streamed line by line (never loaded whole); rows inside each pair's GRCh37 window are written to
  extract_<pair>.tsv in this folder. The full stream is counted to prove nothing was truncated.
Checkpoint: if extract_manifest.json exists with the same source md5s and the extract md5s still match,
the stream-cut is skipped (logged).
"""
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import gzip, hashlib, json, os, sys, time

OUT = (_P31R_V3 + "/R1_stage1/S1R1_R2_03_extcoloc")
RESDIR = (_P31R_V3 + "/R1_stage1/S1R1_R2_00_resources")
GOV = {
    RESDIR + "/STAGE1_BRIEFS.md": "40d0d52f7055542d774c03f1edf04e00",
    RESDIR + "/P31R_R2_STAGE1_PRESPEC_REDACTED.md": "a281c3464619b1b67c28e13e7392c41c",
    RESDIR + "/P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md": "f265d600973c43239bb8377e0cd39d33",
    (_P31R_V3 + "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"): "3319dc39740d9f5b3b8822133d9108b7",
    (_P31R_V3 + "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"): "377508832ea9bd384133e869c456c909",
    RESDIR + "/gls_q.R": "b7390fdf5b6f2f9c1d30c298d6afdafc",
    RESDIR + "/RESOURCES.json": "6d21238155a62fb58f71e05a01638dd0",
}
PAIRS = [  # pair, exposure id, outcome id
    ("PAM_T2D", "EXP_UKBPPP_EUR_PAM_OID20291", "OUT_DIAMANTE2018_noUKB_EUR_chrpos"),
    ("RAD23B_T2D", "EXP_UKBPPP_EUR_RAD23B_OID21166", "OUT_DIAMANTE2018_noUKB_EUR_chrpos"),
    ("TGFB1_CAD", "EXP_UKBPPP_EUR_TGFB1_OID20621", "OUT_CARDIOGRAMplusC4D_ieu-a-7"),
    ("PCSK9_CAD", "EXP_UKBPPP_EUR_PCSK9_OID20235", "OUT_CARDIOGRAMplusC4D_ieu-a-7"),
]
PLACEHOLDERS = {"", "NA", "NaN", "nan", ".", "-", "NULL", "#N/A", "TBD", "Inf", "-Inf"}


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S") + " [extract] " + msg
    print(line, flush=True)
    with open(os.path.join(OUT, "run.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_md5(path, expected, label):
    got = md5(path)
    if got != expected:
        log("MD5 MISMATCH %s: %s expected %s got %s -- STOP" % (label, path, expected, got))
        sys.exit(3)
    log("md5 OK %s %s %s" % (label, expected, path))


def fnum(tok, what):
    if tok in PLACEHOLDERS:
        raise ValueError("placeholder token %r in numeric column %s" % (tok, what))
    return float(tok)


def main():
    t0 = time.time()
    log("python %s; argv %s" % (sys.version.split()[0], sys.argv))
    for p, m in GOV.items():
        assert_md5(p, m, "governing")
    res = {e["id"]: e for e in json.load(open(RESDIR + "/RESOURCES.json", encoding="utf-8"))}
    log("RESOURCES.json: %d entries" % len(res))

    # ---- exposures: parse, allele check, liftover check, GRCh37 window ----
    from pyliftover import LiftOver
    chain = res["REF_LIFTOVER_hg19ToHg38"]
    assert_md5(chain["path"], chain["md5"], "REF_LIFTOVER_hg19ToHg38")
    import pyliftover
    log("pyliftover %s" % getattr(pyliftover, "__version__", "unknown"))
    lo = LiftOver(chain["path"])
    windows = {}
    for pair, eid, oid in PAIRS:
        e = res[eid]
        assert e["role"] == "raw_input", (eid, e["role"])
        assert_md5(e["path"], e["md5"], eid)
        n = n_allele_ok = n_lift_ok = n_unlift = n_mismatch = 0
        chrs, pos37s, pos38_of = set(), [], {}
        with open(e["path"], encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            assert hdr == ["ID", "chrom", "pos", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"], hdr
            for line in fh:
                f = line.rstrip("\n").split("\t")
                n += 1
                for j in (5, 6, 7, 8):
                    fnum(f[j], hdr[j])
                idp = f[0].split(":")
                assert len(idp) == 6 and idp[4] == "imp" and idp[5] == "v1", f[0]
                c37, p37, ref, alt = idp[0], int(idp[1]), idp[2], idp[3]
                assert c37 == f[1], (f[0], f[1])
                chrs.add(c37)
                pos37s.append(p37)
                pos38_of[p37] = int(f[2])
                if ref == f[3] and alt == f[4]:
                    n_allele_ok += 1
                r = lo.convert_coordinate("chr" + c37, p37 - 1)
                if not r:
                    n_unlift += 1
                elif r[0][0] == "chr" + c37 and r[0][1] + 1 == int(f[2]):
                    n_lift_ok += 1
                else:
                    n_mismatch += 1
        assert len(chrs) == 1, chrs
        chrom = chrs.pop()
        # The window is cut in GRCh38 (+/-1 Mb); a variant can sit far away in GRCh37 (build relocation).
        # Bulk = GRCh37 positions within 1.5 Mb of the median; outliers are kept as exact positions.
        med = sorted(pos37s)[len(pos37s) // 2]
        bulk = [q for q in pos37s if abs(q - med) <= 1_500_000]
        outl = sorted(set(q for q in pos37s if abs(q - med) > 1_500_000))
        windows[pair] = dict(exposure_id=eid, outcome_id=oid, chrom=chrom, start37=min(bulk), end37=max(bulk),
                             extra_pos37=outl, extra_pos38=[pos38_of[q] for q in outl],
                             n_exposure_rows=n, n_id_alleles_eq_A0A1=n_allele_ok, n_lift_ok=n_lift_ok,
                             n_unliftable=n_unlift, n_lift_mismatch=n_mismatch)
        log("exposure %s %s: %d rows; chr%s GRCh37 bulk window %d-%d (%.0f kb); GRCh37 outliers (>1.5 Mb from the "
            "median, kept as exact positions): %s with GRCh38 pos %s; ID ref:alt == ALLELE0:ALLELE1 on %d/%d; "
            "hg19->hg38 lift of ID pos == pos column on %d/%d (unliftable %d, mismatched %d)"
            % (pair, eid, n, chrom, min(bulk), max(bulk), (max(bulk) - min(bulk)) / 1e3, outl,
               [pos38_of[q] for q in outl], n_allele_ok, n, n_lift_ok, n, n_unlift, n_mismatch))
        assert n_allele_ok == n, "ID alleles disagree with ALLELE0/ALLELE1"
        assert n_lift_ok == n, "GRCh37 ID position does not lift to the GRCh38 pos column on every row"

    # ---- outcomes: md5 + stream-cut ----
    out_ids = sorted({oid for _, _, oid in PAIRS})
    src = {}
    for oid in out_ids:
        e = res[oid]
        assert e["role"] == "raw_input", (oid, e["role"])
        assert os.path.getsize(e["path"]) == e["bytes"], (oid, os.path.getsize(e["path"]), e["bytes"])
        assert_md5(e["path"], e["md5"], oid)
        src[oid] = e["md5"]
    for pair, eid, _ in PAIRS:
        src[eid] = res[eid]["md5"]

    man_path = os.path.join(OUT, "extract_manifest.json")
    if os.path.exists(man_path):
        man = json.load(open(man_path, encoding="utf-8"))
        ok = man.get("source_md5") == src and all(
            os.path.exists(os.path.join(OUT, x["file"])) and md5(os.path.join(OUT, x["file"])) == x["md5"]
            for x in man["extracts"].values())
        if ok:
            log("checkpoint: extract_manifest.json matches source md5s and extract md5s; stream-cut skipped")
            log("extract step done in %.1f s" % (time.time() - t0))
            return
        log("checkpoint present but stale/mismatched; re-running stream-cut")

    extracts = {}
    # DIAMANTE (plain text, GRCh37)
    oid = "OUT_DIAMANTE2018_noUKB_EUR_chrpos"
    dpairs = [p for p, _, o in PAIRS if o == oid]
    path = res[oid]["path"]
    fhs, counts, bad = {}, {p: 0 for p in dpairs}, {p: 0 for p in dpairs}
    prefixes = tuple(sorted({windows[p]["chrom"] + ":" for p in dpairs}))
    t1 = time.time()
    nrows = 0
    with open(path, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        assert hdr == ["SNP", "Chr", "Pos", "EA", "NEA", "EAF", "Beta", "SE", "Pvalue"], hdr
        for p in dpairs:
            fhs[p] = open(os.path.join(OUT, "extract_%s.tsv" % p), "w", encoding="utf-8", newline="\n")
            fhs[p].write("\t".join(hdr) + "\n")
        for line in fh:
            nrows += 1
            if nrows % 5_000_000 == 0:
                log("DIAMANTE stream: %d rows read (%.0f s); kept %s" % (nrows, time.time() - t1, counts))
            if not line.startswith(prefixes):
                continue
            f = line.rstrip("\n").split("\t")
            for p in dpairs:
                w = windows[p]
                q = int(f[2])
                if f[1] == w["chrom"] and (w["start37"] <= q <= w["end37"] or q in w["extra_pos37"]):
                    assert f[0] == "%s:%s" % (f[1], f[2]), f[0]
                    try:
                        for j in (5, 6, 7, 8):
                            fnum(f[j], hdr[j])
                    except ValueError:
                        bad[p] += 1
                    fhs[p].write(line if line.endswith("\n") else line + "\n")
                    counts[p] += 1
    for p in dpairs:
        fhs[p].close()
    log("DIAMANTE stream complete: %d data rows (RESOURCES notes say 21,508,698 rows by wc -l); %.0f s; kept %s; "
        "rows with a placeholder/non-numeric token in EAF/Beta/SE/Pvalue: %s" % (nrows, time.time() - t1, counts, bad))
    src_rows = {oid: nrows}
    for p in dpairs:
        fn = "extract_%s.tsv" % p
        extracts[p] = dict(file=fn, md5=md5(os.path.join(OUT, fn)), rows=counts[p], nonnumeric_rows=bad[p],
                           source=oid, source_rows_streamed=nrows)

    # ieu-a-7 (BGZF VCF, GRCh37)
    oid = "OUT_CARDIOGRAMplusC4D_ieu-a-7"
    vpairs = [p for p, _, o in PAIRS if o == oid]
    path = res[oid]["path"]
    out_hdr = ["CHROM", "POS", "ID", "REF", "ALT", "ES", "SE", "LP", "AF", "SS", "multiallelic"]
    fhs, counts, bad, multi = {}, {p: 0 for p in vpairs}, {p: 0 for p in vpairs}, {p: 0 for p in vpairs}
    prefixes = tuple(sorted({windows[p]["chrom"] + "\t" for p in vpairs}))
    for p in vpairs:
        fhs[p] = open(os.path.join(OUT, "extract_%s.tsv" % p), "w", encoding="utf-8", newline="\n")
        fhs[p].write("\t".join(out_hdr) + "\n")
    t1 = time.time()
    nrec, total_declared, fmt_seen = 0, None, set()
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##SAMPLE="):
                    for kv in line.strip()[len("##SAMPLE=<"):-1].split(","):
                        if kv.startswith("TotalVariants="):
                            total_declared = int(kv.split("=")[1])
                continue
            if line.startswith("#CHROM"):
                vh = line.rstrip("\n").split("\t")
                assert vh[:9] == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"], vh
                continue
            nrec += 1
            if nrec % 2_000_000 == 0:
                log("ieu-a-7 stream: %d records read (%.0f s); kept %s" % (nrec, time.time() - t1, counts))
            if not line.startswith(prefixes):
                continue
            f = line.rstrip("\n").split("\t")
            for p in vpairs:
                w = windows[p]
                q = int(f[1])
                if f[0] == w["chrom"] and (w["start37"] <= q <= w["end37"] or q in w["extra_pos37"]):
                    keys = f[8].split(":")
                    fmt_seen.add(f[8])
                    vals = dict(zip(keys, f[9].split(":")))
                    is_multi = "," in f[4]
                    if is_multi:
                        multi[p] += 1
                    row = [f[0], f[1], f[2], f[3], f[4]] + [vals.get(k, "NA") for k in ("ES", "SE", "LP", "AF", "SS")] + [str(int(is_multi))]
                    try:
                        for k in ("ES", "SE", "LP", "AF", "SS"):
                            fnum(vals.get(k, "NA"), k)
                    except ValueError:
                        bad[p] += 1
                    fhs[p].write("\t".join(row) + "\n")
                    counts[p] += 1
    for p in vpairs:
        fhs[p].close()
    log("ieu-a-7 stream complete: %d records; ##SAMPLE TotalVariants=%s; %.0f s; kept %s; FORMAT strings in windows %s; "
        "rows with placeholder/non-numeric ES/SE/LP/AF/SS: %s; multi-allelic ALT rows: %s"
        % (nrec, total_declared, time.time() - t1, counts, sorted(fmt_seen), bad, multi))
    assert total_declared is not None and nrec == total_declared, "ieu-a-7 record count != header TotalVariants (truncated stream?)"
    src_rows[oid] = nrec
    for p in vpairs:
        fn = "extract_%s.tsv" % p
        extracts[p] = dict(file=fn, md5=md5(os.path.join(OUT, fn)), rows=counts[p], nonnumeric_rows=bad[p],
                           multiallelic_rows=multi[p], source=oid, source_rows_streamed=nrec)

    for p, x in extracts.items():
        assert x["rows"] > 0, "empty outcome window for %s: blind reader or wrong build" % p
    man = dict(source_md5=src, source_rows=src_rows, windows=windows, extracts=extracts,
               created=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(man_path, "w", encoding="utf-8") as fh:
        json.dump(man, fh, indent=1)
    for p, x in extracts.items():
        log("extract %s: %s rows=%d md5=%s" % (p, x["file"], x["rows"], x["md5"]))
    log("extract step done in %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
