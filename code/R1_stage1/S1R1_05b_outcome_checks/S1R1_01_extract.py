"""Stream each check resource ONCE (read-only) and write (a) the instrument rows and (b) rows at the LD-score SNP set.
Resumable: a resource with a .done file is skipped. Prints counts only, never an instrument value.
Run: python -u S1R1_01_extract.py   (progress every 2M lines)"""
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import gzip, hashlib, os, sys, time

OUT = (_P31R_EXT + "/D/data/paper31_r1_stage1/extracts")
LDRS = (_P31R_EXT + "/D/data/paper31_r1_stage1/finngen_ldsc/ldsc_rs.txt")
os.makedirs(OUT, exist_ok=True)
INSTR = {"rs149802978", "rs7035725", "rs2781668", "rs17788484", "rs11591147", "rs73045269"}
FGFR4 = {"38": ("5", "177084431"), "37": ("5", "176511432")}   # rs451643; build decided at a neutral SNP first
NEUTRAL = {"38": ("19", "44908822"), "37": ("19", "45412079")}  # rs7412 (APOE) -- presence only

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 22), b""):
            h.update(ch)
    return h.hexdigest()

def done(tag):
    return os.path.exists(os.path.join(OUT, tag + ".done"))

def finish(tag, msg):
    with open(os.path.join(OUT, tag + ".done"), "w") as f:
        f.write(msg + "\n")
    print(tag, "DONE", msg, flush=True)

ldrs = None
def load_ld():
    global ldrs
    if ldrs is None:
        with open(LDRS) as f:
            ldrs = {l.strip() for l in f if l.strip()}
        print("LD-score rsIDs:", len(ldrs), flush=True)
    return ldrs

# ---- DIAMANTE 2018 noUKB (rsID version): rsID BP SNP CHR EA NEA EAF Beta SE Pvalue ----
def diamante():
    tag, src = "DIAMANTE", (_P31R_EXT + "/D/data/EUR_new_outcome/Mahajan_T2D/Mahajan_T2D_rsid.txt")
    if done(tag): print(tag, "skip (done)"); return
    rs = load_ld(); t0 = time.time(); n = ni = nl = 0
    with opener(src) as f, open(f"{OUT}/{tag}_instr.tsv", "w") as fi, open(f"{OUT}/{tag}_ldsc.tsv", "w") as fl:
        hdr = f.readline().rstrip("\n").split("\t"); assert hdr[:3] == ["rsID", "BP", "SNP"], hdr
        fi.write("\t".join(hdr) + "\n"); fl.write("\t".join(hdr) + "\n")
        for line in f:
            n += 1
            r = line.split("\t", 1)[0]
            if r in INSTR: fi.write(line); ni += 1
            if r in rs: fl.write(line); nl += 1
            if n % 2_000_000 == 0: print(tag, n, "lines", round(time.time() - t0), "s", flush=True)
    finish(tag, f"src_md5={md5(src)} rows={n} instr_rows={ni} ldsc_rows={nl}")

# ---- CARDIoGRAMplusC4D 1000G, OpenGWAS VCF ieu-a-7: FORMAT ES:SE:LP:AF:SS:... ----
def cardiogram():
    tag, src = "CARDIOGRAM", (_P31R_EXT + "/D/data/EUR_new_outcome/ieu-a-7.vcf.gz")
    if done(tag): print(tag, "skip (done)"); return
    rs = load_ld(); t0 = time.time(); n = ni = nl = 0
    cols = "rsid\tchr\tpos\tref\talt\tbeta_alt\tse\tlp\taf_alt\tss\n"
    with opener(src) as f, open(f"{OUT}/{tag}_instr.tsv", "w") as fi, open(f"{OUT}/{tag}_ldsc.tsv", "w") as fl:
        fi.write(cols); fl.write(cols)
        for line in f:
            if line.startswith("#"): continue
            n += 1
            p = line.rstrip("\n").split("\t")
            r = p[2]
            if r not in INSTR and r not in rs: continue
            keys = p[8].split(":"); vals = dict(zip(keys, p[9].split(":")))
            row = "\t".join([r, p[0], p[1], p[3], p[4], vals.get("ES", "NA"), vals.get("SE", "NA"), vals.get("LP", "NA"),
                             vals.get("AF", "NA"), vals.get("SS", "NA")]) + "\n"
            if r in INSTR: fi.write(row); ni += 1
            if r in rs: fl.write(row); nl += 1
            if n % 2_000_000 == 0: print(tag, n, "records", round(time.time() - t0), "s", flush=True)
    finish(tag, f"src_md5={md5(src)} records={n} instr_rows={ni} ldsc_rows={nl}")

# ---- TPMI SAIGE PheCode 250.2: CHR POS MarkerID Allele1 Allele2 ... BETA (Allele2) SE ... N_case N_ctrl QC ----
def tpmi():
    tag, src = "TPMI_250.2", (_P31R_EXT + "/D/data/twb/saige/250.2.saige.txt.gz")
    if done(tag): print(tag, "skip (done)"); return
    t0 = time.time(); n = 0; hits = {"38": 0, "37": 0}; keep = []
    with opener(src) as f:
        hdr = f.readline().rstrip("\n").split("\t"); assert hdr[0] == "CHR" and hdr[1] == "POS", hdr
        for line in f:
            n += 1
            c, p = line.split("\t", 2)[:2]
            for b, (cc, pp) in NEUTRAL.items():
                if c == cc and p == pp: hits[b] += 1
            for b, (cc, pp) in FGFR4.items():
                if c == cc and p == pp: keep.append((b, line))
            if n % 2_000_000 == 0: print(tag, n, "lines", round(time.time() - t0), "s", flush=True)
    build = "38" if hits["38"] and not hits["37"] else ("37" if hits["37"] and not hits["38"] else "AMBIGUOUS")
    with open(f"{OUT}/{tag}_instr.tsv", "w") as fi:
        fi.write("build_candidate\t" + "\t".join(hdr) + "\n")
        for b, line in keep:
            if b == build: fi.write(b + "\t" + line)
    finish(tag, f"src_md5={md5(src)} rows={n} neutral_hits={hits} build={build} instr_rows_written={sum(b == build for b, _ in keep)}")

if __name__ == "__main__":
    which = sys.argv[1:] or ["diamante", "cardiogram", "tpmi"]
    for w in which:
        {"diamante": diamante, "cardiogram": cardiogram, "tpmi": tpmi}[w]()
