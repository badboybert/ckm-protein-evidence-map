# Read-only extraction of Chen 2023 (Nat Genet) supplementary rows for the three ARG1-locus variants.
# Run: python -I scripts/a01_chen_extract.py   (xlsx is untrusted data: read only, md5 asserted)
# Outputs: chen_rows.tsv (raw rows, sheet/row addressed) and chen_estimates.tsv (named cells, sheet!row!col addressed)
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import openpyxl, hashlib, csv, re
X = (_P31R_V3 + "/R1_stage1/S1R1_R2_09_priorart_downloads_20261010/chen2023_41588_2022_1270_MOESM4_ESM.xlsx")
D = (_P31R_V3 + "/R1_stage1/S1R1_R2_12_chen_instrument/")
h = hashlib.md5(open(X, 'rb').read()).hexdigest()
assert h == "9a245ed156dd900b4e53c9ff4c5a7586", h
wb = openpyxl.load_workbook(X, read_only=True)
VARS = {"rs2246012", "rs17788484", "rs2781668"}
S = {n: list(wb[n].iter_rows(values_only=True)) for n in ["ST5", "ST14", "ST17", "ST19", "ST24"]}
# ---- raw dump: every row of ST5/ST14/ST17/ST19/ST24 that names one of the three rsIDs ----
with open(D + "chen_rows.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sheet", "xlsx_row", "cells (None = empty)"])
    for n, rows in S.items():
        for i, r in enumerate(rows, start=1):
            txt = " ".join("" if c is None else str(c) for c in r)
            if set(re.findall(r"rs\d+", txt)) & VARS:
                w.writerow([n, i] + ["" if c is None else str(c) for c in r])
# ---- named cells ----
def find(sheet, pred):
    hits = [(i, r) for i, r in enumerate(S[sheet], start=1) if pred(r)]
    assert len(hits) == 1, (sheet, len(hits))
    return hits[0]
out = []
def add(variant, sheet, row, what, vals):
    out.append([variant, sheet, row] + [what] + vals)
# ST5 header cols (row 3): 4 SNP,5 CHR,6 Pos38,7 EA,8 NEA,9 FREQ,10 BETA,11 SE,12 P,13 N,14 Metabolite
assert S["ST5"][2][6] == "Position(hg38)" and S["ST5"][2][9] == "FREQ" and S["ST5"][2][10] == "BETA"
for v in ("rs2246012", "rs17788484"):
    i, r = find("ST5", lambda r, v=v: r[4] == v and r[14] == "arginine")
    add(v, "ST5", i, "arginine GWAS (Chen CLSA)", [r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[12], r[13]])
# rs2781668 is not an ST5 row
assert not any(r[4] == "rs2781668" for r in S["ST5"] if r and len(r) > 4)
# ST17 header (row 4): 0 SNP,1 CHR,2 POS_hg38,3 POS_hg37,4 metabolite,5 EA,6 NEA,7 EAF,8 beta,9 se,10 P,11 effector,12 N ; trait blocks start at (label col+4)
assert S["ST17"][3][3] == "POS_hg37" and S["ST17"][3][67] == "Effect Allele"
i, r = find("ST17", lambda r: r[0] == "rs2246012" and r[4] == "arginine")
add("rs2246012", "ST17", i, "instrument row: chr,pos38,pos37,EA,NEA,EAF,beta,se,P,effector,N", [r[1], r[2], r[3], r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[12]])
add("rs2246012", "ST17", i, "T2D block (cols 68-73, 1-based): EA,NEA,EAF,beta,se,P", [r[67], r[68], r[69], r[70], r[71], r[72]])
assert not any(r[0] in ("rs17788484", "rs2781668") for r in S["ST17"][4:] if r and r[4] == "arginine")  # not an ST17 instrument for arginine
# ST19 Wald ratio arginine->T2D ; ST24 coloc arginine/T2D
assert S["ST19"][2][6] == "Beta / Odds ratio"
i, r = find("ST19", lambda r: r[0] == "arginine" and r[1] == "T2D")
add("rs2246012", "ST19", i, "Wald ratio arginine->T2D: SNP,method,nSNP,OR,lo,hi,P,Padj", [r[2], r[4], r[5], r[6], r[7], r[8], r[9], r[10]])
assert S["ST24"][2][10] == "PP.H4.abf"
i, r = find("ST24", lambda r: r[1] == "arginine" and r[2] == "T2D")
add("rs2246012", "ST24", i, "coloc arginine/T2D: SNP,effector,nSNP,PP.H0..PP.H4", [r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]])
# ST14 effector gene
i, r = find("ST14", lambda r: r[0] == "rs2246012" and r[1] == "arginine")
add("rs2246012", "ST14", i, "effector gene,protein type", [r[3], r[4]])
with open(D + "chen_estimates.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["variant", "sheet", "xlsx_row", "what", "values..."])
    for o in out: w.writerow([str(x) for x in o])
print("xlsx md5", h, "| named cells", len(out))
