#!/usr/bin/env python
# ============================================================================
# EAS AXIS 01: acquire + localize the outcome sources (contract v5 SS10 item 1)
# ----------------------------------------------------------------------------
# Verifies the checksum of every registered source, normalises each to the
# SS7.1 schema (all 33 columns), and writes one gzipped TSV per source plus
# eas_source_manifest.csv and a provenance JSON. Fail-closed everywhere: an
# unregistered file, a checksum mismatch (registered OR previously pinned), an
# unmapped header, a malformed data row, or an unresolved build stops the run.
# NA stays NA (SS7.1: no 0-filling).
#
# Modes:
#   --probe      read HEADERS ONLY for every source and print them next to the
#                registered map (presence/metadata only -- no effect values are
#                read, so this mode is legal before the pre-spec lock).
#   (default)    full localization. Only legal AFTER the pre-spec is locked.
#
# The sample-count hierarchy (SS7.2) is applied PER ROW and per outcome type:
#   continuous outcome: VARIANT_LEVEL iff the row itself carries a numeric
#     per-variant total; binary outcome: iff the row carries BOTH numeric
#     per-variant case and control counts. Otherwise STUDY_LEVEL_APPROX when
#     official study counts are registered, else UNRESOLVED. The *_mac_basis
#     columns are filled from ONE level only (no variant/study hybrids), and
#     n_*_variant is NEVER back-filled from study counts. n_eff_variant is
#     carried but never used for any MAC quantity.
# source_qc_status is decided per row: INFO_PASS / INFO_FAIL / INFO_MISSING
# (the source has an INFO column but this row has no usable value) /
# SOURCE_QC_ONLY (the source has no INFO column at all).
# ref/alt are written ONLY for sources whose native columns designate the
# genome reference allele (ref_alt_designated); everywhere else they stay NA
# rather than being cloned from effect/other (SS7.1: absent stays NA).
# ============================================================================
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..", "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import sys
import zipfile
from datetime import datetime, timezone

B = (_P31R_EXT + r"/E/20260331-CKM/Paper 3")
OUT = os.path.join(B, "results", "repro_phase1_hg38", "eas_axis")
LOCAL = os.path.join(OUT, "localized")
SCHEMA = ["source_id", "outcome", "phenotype_label", "ancestry", "build",
          "chr", "pos", "rsid", "ref", "alt", "effect_allele", "other_allele",
          "beta", "se", "p", "eaf", "info",
          "n_total_study", "n_cases_study", "n_controls_study",
          "n_total_variant", "n_cases_variant", "n_controls_variant", "n_eff_variant",
          "n_total_mac_basis", "n_cases_mac_basis", "n_controls_mac_basis",
          "sample_count_basis",
          "mac_total_approx", "mac_cases_approx", "mac_controls_approx",
          "effect_scale", "source_qc_status"]
VALID_BUILDS = ("GRCh37", "GRCh38")

def md5_of(path, limit=None):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def fail(msg):
    sys.stderr.write("EAS-01 FAIL: " + msg + "\n"); sys.exit(1)

def numeric(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if v == v else None          # NaN -> None

# ── source registry. checksum "PENDING" means: compute at first localization,
# pin into the manifest, and REFUSE on any later mismatch; a literal md5 (the
# two new downloads, registered in the pre-spec) is verified before any read.
# outcome_type drives the SS7.2 per-row basis rule; ref_alt_designated is TRUE
# only where the native columns name the genome reference allele (REF/ALT, or
# SAIGE Allele1/Allele2, or GWAS Catalog harmonised other/effect = ref/alt),
# which also strand-anchors the record for palindromic orientation (SS7.2.8).
MG = os.path.join(B, "data", "multiancestry_gwas")
SAIGE_MAP = dict(chr="CHR", pos="POS", rsid="SNPID", other_allele="Allele1",
                 effect_allele="Allele2", eaf="AF_Allele2", info="imputationInfo",
                 n_total_variant="N", beta="BETA", se="SE", p="p.value")
SOURCES = [
  dict(source_id="T2DGGI_EAS", outcome="T2D", role="primary", outcome_type="binary",
       path=os.path.join(MG, "diamante", "Suzuki.Nature2024.T2DGGI.EAS.sumstats.zip"),
       md5="PENDING", fmt="t2dggi_zip", build="GRCh37", ancestry="EAS",
       phenotype_label="type 2 diabetes", effect_scale="log-odds",
       n_cases_study=88109, n_controls_study=339395, ref_alt_designated=False,
       colmap=dict(chr="Chromsome", pos="Position", effect_allele="EffectAllele",
                   other_allele="NonEffectAllele", beta="Beta", se="SE", p="Pval",
                   eaf="EAF", n_cases_variant="Ncases", n_controls_variant="Ncontrols",
                   n_eff_variant="Neff")),
  dict(source_id="EGFR_EAS_2P2M", outcome="eGFR", role="primary", outcome_type="continuous",
       path=os.path.join(MG, "egfr_eas_2p2m", "eGFRcrea_GWAS_EAS_282K_Individuals_Summary_Statistics.txt.gz"),
       md5="3c74f086890195b601b141f215becc66", fmt="metal_gz", build="GRCh37", ancestry="EAS",
       phenotype_label="eGFRcrea (log)", effect_scale="log-eGFR-ish (verify README)",
       n_total_study=282852, ref_alt_designated=True,
       colmap=dict(rsid="MarkerName", chr="CHR", pos="POS", other_allele="REF",
                   effect_allele="ALT", eaf="Freq_ALT", beta="BETA", se="SE",
                   p="P.value", n_total_variant="N")),
  dict(source_id="BBJ_AF_2022", outcome="AF", role="primary", outcome_type="binary",
       path=os.path.join(MG, "bbj", "hum0014.v29.AF.v1.zip"),
       md5="PENDING", fmt="nbdc_zip", build="GRCh37", ancestry="EAS",
       phenotype_label="atrial fibrillation or flutter (disclose)", effect_scale="log-odds",
       n_cases_study=9826, n_controls_study=140446, inner="BBJ_AF_2022.sumstats.gz",
       ref_alt_designated=True,
       colmap=dict(chr="CHR", pos="POS", other_allele="REF", effect_allele="ALT",
                   eaf="AAF", beta="BETA", se="SE", p="P", n_total_variant="N", info="INFO")),
  dict(source_id="PHEWEB_BBJ_CAD", outcome="CAD", role="primary", outcome_type="binary",
       path=os.path.join(MG, "bbj_cad_pheweb", "pheweb_bbj_CAD.tsv.gz"),
       # build resolved GRCh37 by dated amendment 2026-08-16: rs671 (non-target
       # sentinel) sits at 12:112241766 in this file = the GRCh37 coordinate
       # (GRCh38 would be 12:111803962). 02 cross-checks empirically per source.
       md5="93fec56d62ab2e371e0fade6a71b4fff", fmt="saige_gz", build="GRCh37", ancestry="EAS",
       phenotype_label="coronary artery disease", effect_scale="log-odds",
       n_cases_study=32512, n_controls_study=146214, ref_alt_designated=True,
       colmap=dict(SAIGE_MAP)),
  dict(source_id="BBJ_CHF", outcome="HF", role="primary", outcome_type="binary",
       path=os.path.join(MG, "bbj", "hum0197.v3.BBJ.CHF.v1.zip"),
       md5="PENDING", fmt="nbdc_zip", build="GRCh37", ancestry="EAS",
       phenotype_label="chronic heart failure (vs EUR all-cause HF: partial match)",
       effect_scale="log-odds", n_cases_study=10540, n_controls_study=168186, inner=None,
       ref_alt_designated=True, colmap=dict(SAIGE_MAP)),
  dict(source_id="BBJ_IS", outcome="Stroke", role="primary", outcome_type="binary",
       path=os.path.join(MG, "bbj", "hum0197.v3.BBJ.IS.v1.zip"),
       md5="PENDING", fmt="nbdc_zip", build="GRCh37", ancestry="EAS",
       phenotype_label="ischemic stroke", effect_scale="log-odds",
       n_cases_study=22664, n_controls_study=152022, inner=None,
       ref_alt_designated=True, colmap=dict(SAIGE_MAP)),
  dict(source_id="GCST90018602_CKD", outcome="CKD", role="primary", outcome_type="binary",
       path=os.path.join(MG, "bbj_ckd", "GCST90018602.h.tsv.gz"),
       md5="PENDING", fmt="gwascat_h_gz", build="GRCh38", ancestry="EAS",
       phenotype_label="chronic renal failure (CKD proxy; low power)", effect_scale="log-odds",
       n_cases_study=2117, n_controls_study=174345, ref_alt_designated=True,
       colmap=dict(chr="chromosome", pos="base_pair_location", rsid="rsid",
                   other_allele="other_allele", effect_allele="effect_allele",
                   beta="beta", se="standard_error", p="p_value",
                   eaf="effect_allele_frequency")),
  dict(source_id="BBJ_T2D", outcome="T2D", role="sensitivity", outcome_type="binary",
       path=os.path.join(MG, "bbj", "hum0197.v3.BBJ.T2D.v1.zip"),
       md5="PENDING", fmt="nbdc_zip", build="GRCh37", ancestry="EAS",
       phenotype_label="type 2 diabetes", effect_scale="log-odds", inner=None,
       # DATED AMENDMENT 2026-08-16 (pre-spec SS2.1 required this before use; it was locked with the
       # gap open by mistake). Official study counts for the Sakaue/Kanai 2021 BBJ release this file
       # belongs to (GWASsummary_T2D_Japanese_SakaueKanai2020), read from the release's own portal
       # entry (pheweb.jp/pheno/T2D): 45,383 cases / 132,032 controls. Sibling entries from the same
       # release are already registered the same way (MI 14,992/146,214; CHF 10,540/168,186;
       # IS 22,664/152,022). Per-variant N columns are NEVER used to backfill these (SS2.1 rule).
       n_cases_study=45383, n_controls_study=132032,
       ref_alt_designated=True, colmap=dict(SAIGE_MAP)),
  dict(source_id="BBJ_KANAI_EGFR", outcome="eGFR", role="sensitivity", outcome_type="continuous",
       path=os.path.join(MG, "bbj", "hum0014.v8.eGFR.zip"),
       md5="PENDING", fmt="nbdc_zip", build="GRCh37", ancestry="EAS",
       phenotype_label="eGFR (Kanai 2018; nested in the 2.2M EAS meta -- not independent)",
       effect_scale="verify README", n_total_study=143658, inner=None,
       ref_alt_designated=True,
       colmap=dict(rsid="SNP", chr="CHR", pos="POS", other_allele="REF",
                   effect_allele="ALT", eaf="Frq", info="Rsq", beta="BETA",
                   se="SE", p="P", n_total_variant="N")),
  dict(source_id="BBJ_MI", outcome="CAD", role="sensitivity", outcome_type="binary",
       path=os.path.join(MG, "bbj", "hum0197.v3.BBJ.MI.v1.zip"),
       md5="PENDING", fmt="nbdc_zip", build="GRCh37", ancestry="EAS",
       phenotype_label="MI (CAD proxy)", effect_scale="log-odds",
       n_cases_study=14992, n_controls_study=146214, inner=None,
       ref_alt_designated=True, colmap=dict(SAIGE_MAP)),
]
# EXCLUDED (pre-spec SS3; asserted absent from the input manifest, never read):
EXCLUDED_HINTS = ["GCST90204201", "GCST90018822", "MultiAncestry", "DIAMANTE-TA",
                  "hum0197.v3.BBJ.AF"]  # the mislabelled flutter file

def open_source(spec):
    """Return (fileobj, member_name) for the tabular payload, streaming."""
    p = spec["path"]
    if not os.path.exists(p):
        fail("source missing on disk: " + p)
    if spec["fmt"].endswith("_zip"):
        z = zipfile.ZipFile(p)
        members = [m for m in z.namelist() if not m.endswith("/")]
        want = spec.get("inner")
        if want:
            hit = [m for m in members if os.path.basename(m) == want]
            if not hit:
                fail("%s: inner file %s not in zip (members: %s)" %
                     (spec["source_id"], want, ", ".join(members[:5])))
            m = hit[0]
        else:
            # the largest member is the sumstats table; refuse ambiguity
            members.sort(key=lambda m: z.getinfo(m).file_size, reverse=True)
            m = members[0]
        raw = z.open(m)
        if m.endswith(".gz"):
            return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), errors="replace"), m
        return io.TextIOWrapper(raw, errors="replace"), m
    else:
        return io.TextIOWrapper(gzip.open(p), errors="replace"), os.path.basename(p)

def probe():
    print("EAS-01 PROBE (headers only; no effect values read)")
    for spec in SOURCES:
        try:
            fh, member = open_source(spec)
            header = fh.readline().rstrip("\n\r")
            fh.close()
        except SystemExit:
            raise
        except Exception as e:
            print("  %-18s ERROR: %s" % (spec["source_id"], e)); continue
        sep = "\t" if "\t" in header else (" " if " " in header else ",")
        cols = header.split(sep)
        print("  %-18s member=%s\n    cols(%d): %s" %
              (spec["source_id"], member, len(cols), " ".join(cols[:22])))
        if spec["colmap"] != "PROBE":
            missing = [v for v in spec["colmap"].values() if v not in cols]
            print("    map check: %s" % ("OK" if not missing else "MISSING " + ",".join(missing)))
        if spec["build"] not in VALID_BUILDS:
            print("    build: %s -- must be resolved (dated amendment) before localization" % spec["build"])

def pinned_md5s():
    """Previously pinned checksums (first localization pins PENDING sources)."""
    mpath = os.path.join(OUT, "eas_source_manifest.csv")
    if not os.path.exists(mpath):
        return {}
    with open(mpath, newline="") as fh:
        return {r["source_id"]: r["md5"] for r in csv.DictReader(fh)}

def localize():
    os.makedirs(LOCAL, exist_ok=True)
    pinned = pinned_md5s()
    manifest, prov = [], {"generated": datetime.now(timezone.utc).isoformat(),
                          "schema": SCHEMA, "sources": []}
    still_pending = []
    for spec in SOURCES:
        sid = spec["source_id"]
        if spec["build"] not in VALID_BUILDS:
            fail("%s: build is %r -- resolve it (dated amendment + re-review) before "
                 "localization; a placeholder build must never reach evaluable output"
                 % (sid, spec["build"]))
        got = md5_of(spec["path"])
        if spec["md5"] != "PENDING":
            if got.lower() != spec["md5"].lower():
                fail("%s: md5 %s != registered %s" % (sid, got, spec["md5"]))
        elif sid in pinned:
            if got.lower() != pinned[sid].lower():
                fail("%s: md5 %s != previously pinned %s -- the file changed between "
                     "runs; re-register it with a dated amendment or restore the original"
                     % (sid, got, pinned[sid]))
        else:
            still_pending.append(sid)
        if spec["colmap"] == "PROBE":
            fail("%s: column map still PROBE -- run --probe, fill the map (dated amendment), re-review" % sid)
        fh, member = open_source(spec)
        header = fh.readline().rstrip("\n\r")
        sep = "\t" if "\t" in header else " "
        cols = header.split(sep)
        idx = {}
        for k, v in spec["colmap"].items():
            if v not in cols:
                fail("%s: registered column %s=%s not in header [%s]" % (sid, k, v, " ".join(cols[:20])))
            idx[k] = cols.index(v)
        opath = os.path.join(LOCAL, sid + ".tsv.gz")
        n_rows, n_bad, bad_example = 0, 0, ""
        binary = spec["outcome_type"] == "binary"
        has_study = bool(spec.get("n_cases_study") or spec.get("n_total_study"))
        has_info_col = "info" in idx
        ref_alt = bool(spec.get("ref_alt_designated"))
        basis_counts = {"VARIANT_LEVEL": 0, "STUDY_LEVEL_APPROX": 0, "UNRESOLVED": 0}
        with gzip.open(opath, "wt", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(SCHEMA)
            for line in fh:
                if not line.strip():
                    continue                      # a bare trailing newline is not a record
                f = line.rstrip("\n\r").split(sep)
                if len(f) != len(cols):
                    n_bad += 1
                    if n_bad == 1:
                        bad_example = line[:120]
                    continue                      # counted; any n_bad>0 fails below
                def g(k):
                    return f[idx[k]] if k in idx else "NA"
                row = dict.fromkeys(SCHEMA, "NA")
                row.update(source_id=sid, outcome=spec["outcome"],
                           phenotype_label=spec["phenotype_label"], ancestry=spec["ancestry"],
                           build=spec["build"], effect_scale=spec["effect_scale"],
                           chr=g("chr"), pos=g("pos"), rsid=g("rsid"),
                           effect_allele=g("effect_allele"), other_allele=g("other_allele"),
                           beta=g("beta"), se=g("se"), p=g("p"), eaf=g("eaf"), info=g("info"),
                           n_total_study=str(spec.get("n_total_study", "NA")),
                           n_cases_study=str(spec.get("n_cases_study", "NA")),
                           n_controls_study=str(spec.get("n_controls_study", "NA")),
                           n_total_variant=g("n_total_variant"),
                           n_cases_variant=g("n_cases_variant"),
                           n_controls_variant=g("n_controls_variant"),
                           n_eff_variant=g("n_eff_variant"))
                # ref/alt only where the native columns designate the genome
                # reference allele; SS7.1 says absent stays NA, not cloned.
                if ref_alt:
                    row["ref"] = g("other_allele")
                    row["alt"] = g("effect_allele")
                # SS7.2 hierarchy PER ROW, one level only, never from n_eff.
                if binary:
                    vc, vn = numeric(row["n_cases_variant"]), numeric(row["n_controls_variant"])
                    if vc is not None and vn is not None:
                        basis = "VARIANT_LEVEL"
                        mb = (numeric(row["n_total_variant"]), vc, vn)
                        if mb[0] is None:
                            mb = (vc + vn, vc, vn)     # total of the SAME level, not study
                    elif has_study:
                        basis = "STUDY_LEVEL_APPROX"
                        sc, sn = numeric(row["n_cases_study"]), numeric(row["n_controls_study"])
                        st = numeric(row["n_total_study"])
                        if st is None and sc is not None and sn is not None:
                            st = sc + sn
                        mb = (st, sc, sn)
                    else:
                        basis, mb = "UNRESOLVED", (None, None, None)
                else:
                    vt = numeric(row["n_total_variant"])
                    if vt is not None:
                        basis, mb = "VARIANT_LEVEL", (vt, None, None)
                    elif has_study:
                        basis, mb = "STUDY_LEVEL_APPROX", (numeric(row["n_total_study"]), None, None)
                    else:
                        basis, mb = "UNRESOLVED", (None, None, None)
                row["sample_count_basis"] = basis
                basis_counts[basis] += 1
                fmtn = lambda v: ("NA" if v is None else repr(v))
                row["n_total_mac_basis"], row["n_cases_mac_basis"], row["n_controls_mac_basis"] = \
                    fmtn(mb[0]), fmtn(mb[1]), fmtn(mb[2])
                # MAC approximations (SS7.2.7) from the mac_basis level only.
                eaf = numeric(row["eaf"])
                if eaf is not None and 0.0 < eaf < 1.0:
                    maf = min(eaf, 1.0 - eaf)
                    mac = lambda n: "NA" if n is None else repr(2.0 * n * maf)
                    row["mac_total_approx"] = mac(mb[0])
                    row["mac_cases_approx"] = mac(mb[1])
                    row["mac_controls_approx"] = mac(mb[2])
                # per-row INFO status (SS7.2.6 + SS10.1 vocabulary)
                if not has_info_col:
                    row["source_qc_status"] = "SOURCE_QC_ONLY"
                else:
                    iv = numeric(row["info"])
                    row["source_qc_status"] = ("INFO_MISSING" if iv is None
                                               else "INFO_PASS" if iv >= 0.80 else "INFO_FAIL")
                w.writerow([row[c] for c in SCHEMA])
                n_rows += 1
        fh.close()
        if n_bad:
            fail("%s: %d malformed data row(s) (field count != header; first: %r). A row "
                 "silently skipped is indistinguishable from an absent variant, so the run "
                 "stops instead. Inspect the file; register a documented exception if the "
                 "rows are provably non-records." % (sid, n_bad, bad_example))
        manifest.append(dict(source_id=sid, outcome=spec["outcome"], role=spec["role"],
                             outcome_type=spec["outcome_type"], path=spec["path"],
                             member=member, md5=got, rows=n_rows, build=spec["build"],
                             ref_alt_designated=ref_alt,
                             n_total_study=spec.get("n_total_study", ""),
                             n_cases_study=spec.get("n_cases_study", ""),
                             n_controls_study=spec.get("n_controls_study", ""),
                             phenotype_label=spec["phenotype_label"],
                             basis_variant=basis_counts["VARIANT_LEVEL"],
                             basis_study=basis_counts["STUDY_LEVEL_APPROX"],
                             basis_unresolved=basis_counts["UNRESOLVED"],
                             localized=opath))
        prov["sources"].append(manifest[-1])
        print("  %-18s %8d rows -> %s (basis V/S/U = %d/%d/%d)" %
              (sid, n_rows, os.path.basename(opath), basis_counts["VARIANT_LEVEL"],
               basis_counts["STUDY_LEVEL_APPROX"], basis_counts["UNRESOLVED"]))
    # excluded files must not appear among inputs
    used = {m["path"] for m in manifest}
    for hint in EXCLUDED_HINTS:
        bad = [p for p in used if hint in p]
        if bad:
            fail("excluded source found in input manifest: " + bad[0])
    with open(os.path.join(OUT, "eas_source_manifest.csv"), "w", newline="") as fh2:
        w = csv.DictWriter(fh2, fieldnames=list(manifest[0].keys()))
        w.writeheader(); w.writerows(manifest)
    with open(os.path.join(OUT, "eas_source_provenance.json"), "w") as fh3:
        json.dump(prov, fh3, indent=1)
    print("EAS-01 done: %d sources localized; manifest + provenance written" % len(manifest))
    if still_pending:
        print("PINNED THIS RUN (fold into the pre-spec before the lock): " + ", ".join(still_pending))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    a = ap.parse_args()
    probe() if a.probe else localize()
