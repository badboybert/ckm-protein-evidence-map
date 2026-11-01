# -*- coding: utf-8 -*-
# Stage 3b, track T1 (Table 1 and Additional file 1: Table S2, data sources): ledger rows L300-L339.
# Executed by S1R1_build_ledger.py in ITS namespace (add, rel, pd, np, json, V3, CKM, S1, ROWS, EXT_CHECKS).
# Every name defined here starts with _t1_ so nothing in the builder's namespace is overwritten.
# Each value is read from its source file and asserted; a value stated only in a publication carries
# source 'PMID n (where stated)' and tag [P], and was read in that publication (abstract or PMC full text) on 2026-10-06.
# EXT_CHECKS records a cross-source match made by the T1 builder (a second, different source), not a separate verifier.
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import ast as _t1_ast, gzip as _t1_gzip, io as _t1_io, os as _t1_os, re as _t1_re

_t1_D = (_P31R_EXT + "/D/data")


def _t1_row(rid):
    _r = [r_ for r_ in ROWS if r_["id"] == rid]
    assert len(_r) == 1, rid
    return json.loads(_r[0]["value"])


def _t1_cc(a, b):
    return "{:,} / {:,}".format(int(a), int(b))


# ---------------------------------------------------------------- the analysis configuration (declared outcome N)
_t1_cfgp = _t1_os.path.join(V3, "paper3_deposit_request_v4/scripts/R/00_config.R")
_t1_cfg = _t1_io.open(_t1_cfgp, encoding="utf-8").read()


def _t1_rvec(name):
    _m = _t1_re.search(name + r"\s*<-\s*c\(([^)]*)\)", _t1_cfg, flags=_t1_re.S)
    assert _m, name
    _o = {}
    for _k, _v in _t1_re.findall(r"(\w+)\s*=\s*(NA_integer_|\d+)L?", _m.group(1)):
        _o[_k] = None if _v.startswith("NA") else int(_v)
    return _o


_t1_ntot, _t1_ncase = _t1_rvec("OUTCOME_N_TOTAL"), _t1_rvec("OUTCOME_N_CASE")
_t1_ids = dict(_t1_re.findall(r'"(\w+)"\s*=\s*"((?:FinnGen|MEGASTROKE|Wuttke)\w+)"', _t1_cfg))
assert _t1_ids == {"CAD": "FinnGen_R10_I9_CHD", "Stroke": "MEGASTROKE_EUR_GCST006906", "HF": "FinnGen_R10_I9_HEARTFAIL",
                   "T2D": "FinnGen_R10_T2D", "eGFR": "Wuttke2019_eGFR_EA", "CKD": "Wuttke2019_CKD_EA",
                   "NAFLD": "FinnGen_R10_NAFLD", "AF": "FinnGen_R10_I9_AF"}, _t1_ids

# ---------------------------------------------------------------- L300-L301: exposure (UKB-PPP)
add("L300", "Table 1 / Table S2: exposure", "UKB-PPP European-ancestry discovery cohort for pQTL mapping (Olink Explore 3072; protein levels inverse-rank normalised)",
    34557, "34,557", "PMID 37794186 (where stated): Sun et al. 2023, PMC10567551 full text, 'Discovery pQTL analyses were performed in participants of European ancestry ... (n = 34,557)'",
    "Results: pQTL discovery cohort", "[P]")
_t1_ins = pd.read_csv(_t1_os.path.join(V3, "V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv"), usecols=["gene", "samplesize"])
_t1_ss = _t1_ins.samplesize.astype(float)
assert len(_t1_ins) == 1964 and _t1_row("L200")[0] == len(_t1_ins), "instrument count differs from L200"
assert _t1_ss.max() <= 34557 and float(_t1_ss.median()).is_integer()
add("L301", "Table S2: exposure", "UKB-PPP per-protein sample size at the instruments (samplesize column of the instrument table): instruments / minimum / median / maximum",
    [len(_t1_ins), int(_t1_ss.min()), int(_t1_ss.median()), int(_t1_ss.max())],
    "%d instruments; N %s to %s (median %s)" % (len(_t1_ins), "{:,}".format(int(_t1_ss.min())), "{:,}".format(int(_t1_ss.max())), "{:,}".format(int(_t1_ss.median()))),
    "V3/V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv", "samplesize", "[D]")

# ---------------------------------------------------------------- L302-L303: exposure-side corroboration cohorts
add("L302", "Table 1 / Table S2: corroboration", "deCODE plasma proteomics: Icelanders / SomaScan aptamers", [35559, 4907], "35,559 Icelanders; 4,907 aptamers",
    "PMID 34857953 (where stated): Ferkingstad et al. 2021, abstract", "abstract", "[P]")
_t1_aric = _t1_io.open(_t1_os.path.join(V3, "paper3_deposit_request_v4/scripts/reports/Paper3_Program/run_logs/RERUN_20260915/12b_aric.log"), encoding="utf-8", errors="replace").read()
assert "data/aric/EA.zip" in _t1_aric, "the ARIC run did not read the European-American file"
add("L303", "Table 1 / Table S2: corroboration", "ARIC plasma proteomics used: European-American participants (the EA.zip file read by the ARIC run)", 7213, "7,213 European Americans",
    "PMID 35501419 (where stated): Zhang et al. 2022, abstract; file identity from run_logs/RERUN_20260915/12b_aric.log", "abstract; zip: EA.zip", "[P]")

# ---------------------------------------------------------------- L304-L308: FinnGen R10 outcomes (release manifest)
_t1_man = pd.read_csv(_t1_os.path.join(_t1_D, "paper31_v4_raw/ukbfree_outcomes/R10_manifest.tsv"), sep="\t")
_t1_fg = [("L304", "CAD", "I9_CHD"), ("L305", "HF", "I9_HEARTFAIL"), ("L306", "AF", "I9_AF"), ("L307", "T2D", "T2D"), ("L308", "NAFLD", "NAFLD")]
for _t1_id, _t1_o, _t1_code in _t1_fg:
    _t1_r = _t1_man[_t1_man.phenocode == _t1_code]
    assert len(_t1_r) == 1, _t1_code
    _t1_ca, _t1_co = int(_t1_r.num_cases.iloc[0]), int(_t1_r.num_controls.iloc[0])
    assert _t1_ids[_t1_o] == "FinnGen_R10_" + _t1_code
    assert _t1_ncase[_t1_o] == _t1_ca and _t1_ntot[_t1_o] == _t1_ca + _t1_co, (_t1_o, "manifest differs from the analysis configuration")
    add(_t1_id, "Table 1 / Table S2: outcomes", "FinnGen R10 %s (%s, '%s'): cases / controls / total" % (_t1_o, _t1_code, _t1_r.phenotype.iloc[0]),
        [_t1_ca, _t1_co, _t1_ca + _t1_co], _t1_cc(_t1_ca, _t1_co),
        (_P31R_EXT + "/D/data/paper31_v4_raw/ukbfree_outcomes/R10_manifest.tsv (release md5 9a6b2ab7, FETCH_LOG); asserted equal to OUTCOME_N_CASE/OUTCOME_N_TOTAL in paper3_deposit_request_v4/scripts/R/00_config.R"),
        "phenocode %s: num_cases, num_controls, phenotype" % _t1_code, "[D]")
    EXT_CHECKS[_t1_id] = "MATCH (cross-source, T1 builder 2026-10-06): FinnGen R10 public API r10.finngen.fi/api/pheno/%s gives %d cases / %d controls" % (_t1_code, _t1_ca, _t1_co)

# ---------------------------------------------------------------- L309-L312: MEGASTROKE and CKDGen
_t1_sc, _t1_sn = _t1_ncase["Stroke"], _t1_ntot["Stroke"]
assert (_t1_sc, _t1_sn - _t1_sc) == (40585, 406111), "MEGASTROKE European any stroke differs from Malik 2018 (40,585 cases; 406,111 controls)"
add("L309", "Table 1 / Table S2: outcomes", "MEGASTROKE European any stroke (GCST006906): cases / controls / total", [_t1_sc, _t1_sn - _t1_sc, _t1_sn], _t1_cc(_t1_sc, _t1_sn - _t1_sc),
    "paper3_deposit_request_v4/scripts/R/00_config.R (OUTCOME_N_CASE/OUTCOME_N_TOTAL Stroke); PMID 29531354 (where stated)", "Stroke", "[D][P]")
EXT_CHECKS["L309"] = "MATCH (cross-source, T1 builder 2026-10-06): Malik 2018 full text (PMC5968830) 'European participants only (40,585 cases; 406,111 controls)'; GWAS Catalog GCST006906 initial sample 40,585 / 406,111"
_t1_eg = _t1_ntot["eGFR"]
assert _t1_eg == 567460 and _t1_ncase["eGFR"] is None
add("L310", "Table 1 / Table S2: outcomes", "CKDGen European-ancestry creatinine-based eGFR (Wuttke 2019, no UK Biobank): N", _t1_eg, "{:,}".format(_t1_eg),
    "paper3_deposit_request_v4/scripts/R/00_config.R (OUTCOME_N_TOTAL eGFR); PMID 31152163 (where stated)", "eGFR", "[D][P]")
EXT_CHECKS["L310"] = ("MATCH (cross-source, T1 builder 2026-10-06): Wuttke 2019 full text (PMC6698888) 'European (n = 567,460)'; maximum n_total_sum over 8,885,712 rows of " + _P31R_EXT + "/D/data/ckdgen/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz = 567,460")
_t1_kc, _t1_kn = _t1_ncase["CKD"], _t1_ntot["CKD"]
assert (_t1_kc, _t1_kn - _t1_kc, _t1_kn) == (41395, 439303, 480698)
add("L311", "Table 1 / Table S2: outcomes", "CKDGen European-ancestry CKD (Wuttke 2019, no UK Biobank): cases / controls / total", [_t1_kc, _t1_kn - _t1_kc, _t1_kn], _t1_cc(_t1_kc, _t1_kn - _t1_kc),
    "paper3_deposit_request_v4/scripts/R/00_config.R (OUTCOME_N_CASE/OUTCOME_N_TOTAL CKD)", "CKD", "[D]")
EXT_CHECKS["L311"] = "MATCH (cross-source, T1 builder 2026-10-06): CKDGen download page for Wuttke 2019 ('CKD EA: n=480,698 (41,395 cases and 439,303 controls)', file CKD_overall_EA_JW_20180223_nstud23)"
add("L312", "Table 1: outcomes", "CKDGen CKD definition: eGFR below 60 mL/min per 1.73 m2", [60, 1.73], "eGFR < 60 mL/min per 1.73 m2",
    "PMID 31152163 (where stated): Wuttke 2019 Methods, PMC6698888 ('CKD was defined as an eGFR below 60 ml min-1 per 1.73 m2')", "Methods", "[P]")

# ---------------------------------------------------------------- L313-L314: European external outcome checks
_t1_dc, _t1_dn = 74124 - 19119, 824006 - 423698
assert [_t1_dc, _t1_dn] == _t1_row("L218"), "DIAMANTE without UK Biobank differs from L218"
add("L313", "Table 1 / Table S2: external checks", "DIAMANTE 2018 European T2D without UK Biobank: cases / controls (all-sample totals minus UK Biobank)",
    [_t1_dc, _t1_dn, 74124, 824006, 19119, 423698], _t1_cc(_t1_dc, _t1_dn),
    "PMID 30297969 (where stated): Mahajan 2018 full text, PMC6287706: 74,124 cases and 824,006 controls in all; UK Biobank 19,119 cases and 423,698 controls",
    "Results 'Study overview'; Methods 'Estimating phenotypic variance explained by SNPs'", "[P]",
    note="Re-derives L218 (single-check row) from the publication; the alternative_outcomes_summary.csv on D: carries the WITH-UKB totals for the noUKB file.")
EXT_CHECKS["L313"] = "MATCH (cross-source, T1 builder 2026-10-06): equals L218, the Paper 4 adjudication recorded in S1R1_PRESPEC_outcome_checks.md (N 455,313)"
with _t1_gzip.open(_t1_os.path.join(_t1_D, "EUR_new_outcome/ieu-a-7.vcf.gz"), "rt") as _t1_fh:
    for _t1_ln in _t1_fh:
        if not _t1_ln.startswith("#"):
            break
_t1_fmt = _t1_ln.rstrip("\n").split("\t")
_t1_ssv = int(float(dict(zip(_t1_fmt[8].split(":"), _t1_fmt[9].split(":")))["SS"]))
assert _t1_ssv == 60801 + 123504 and _t1_row("L219")[:2] == [60801, 123504]
add("L314", "Table 1 / Table S2: external checks", "CARDIoGRAMplusC4D 2015 (ieu-a-7): cases / controls; ancestry shares European / South Asian / East Asian (%); VCF sample size per variant",
    [60801, 123504, 77, 13, 6, _t1_ssv], "60,801 / 123,504; 77% European, 13% South Asian, 6% East Asian",
    ("PMID 26343387 (where stated): Nikpay 2015 full text, PMC4589895 ('60,801 cases and 123,504 controls from 48 studies'; 'The majority (77%) ... European ancestry; 13% and 6% ... south and east Asian'); " + _P31R_EXT + "/D/data/EUR_new_outcome/ieu-a-7.vcf.gz SS field"),
    "SS = 184,305 = cases + controls", "[P][D]")
EXT_CHECKS["L314"] = "MATCH (cross-source, T1 builder 2026-10-06): equals L219 (alternative_outcomes_summary.csv) and the ieu-a-7 VCF per-variant SS (184,305)"

# ---------------------------------------------------------------- L315-L317: measured sampling-noise correlation with FinnGen
_t1_tp = pd.read_csv(_t1_os.path.join(S1, "S1R1_02_topology_note/S1R1_topo_supp_table_a_intercepts.tsv"), sep="\t")
_t1_tp = _t1_tp[_t1_tp.outcome_1 != _t1_tp.outcome_2]
_t1_fgset = {"CAD", "HF", "AF", "T2D", "NAFLD"}


def _t1_pairs(other):
    _m = _t1_tp[((_t1_tp.outcome_1 == other) & _t1_tp.outcome_2.isin(_t1_fgset)) | ((_t1_tp.outcome_2 == other) & _t1_tp.outcome_1.isin(_t1_fgset))]
    assert len(_m) == 5, (other, len(_m))
    return _m.null_z_correlation.astype(float)


_t1_ms, _t1_ck = _t1_pairs("Stroke"), _t1_pairs("eGFR")
_t1_ff = _t1_tp[_t1_tp.outcome_1.isin(_t1_fgset) & _t1_tp.outcome_2.isin(_t1_fgset)].null_z_correlation.astype(float)
_t1_f4 = _t1_tp[_t1_tp.outcome_1.isin(_t1_fgset - {"NAFLD"}) & _t1_tp.outcome_2.isin(_t1_fgset - {"NAFLD"})].null_z_correlation.astype(float)
assert len(_t1_ff) == 10 and len(_t1_f4) == 6
assert abs(_t1_f4.min() - _t1_row("L206")[0]) < 1e-12 and abs(_t1_f4.max() - _t1_row("L206")[1]) < 1e-12, "the four-outcome subset no longer reproduces L206"
add("L315", "Table 1: overlap", "Null-z correlation (cross-trait LDSC intercept) between MEGASTROKE any stroke and each FinnGen outcome: minimum / maximum over 5 pairs",
    [float(_t1_ms.min()), float(_t1_ms.max())], "%.3f to %.3f" % (_t1_ms.min(), _t1_ms.max()),
    "V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_supp_table_a_intercepts.tsv", "null_z_correlation, pairs Stroke x {CAD,HF,AF,T2D,NAFLD}", "[R]")
add("L316", "Table 1: overlap", "Null-z correlation (cross-trait LDSC intercept) between CKDGen eGFR and each FinnGen outcome: minimum / maximum over 5 pairs (CKD not estimated)",
    [float(_t1_ck.min()), float(_t1_ck.max())], "%.4f to %.4f" % (_t1_ck.min(), _t1_ck.max()),
    "V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_supp_table_a_intercepts.tsv", "null_z_correlation, pairs eGFR x {CAD,HF,AF,T2D,NAFLD}", "[R]")
# Table 1 states the four-outcome range (L206, the set the Results and the Supplementary Note name) and, separately, the
# four pairs that include NAFLD; together they cover all ten FinnGen pairs (asserted).
_t1_fn = _t1_tp[((_t1_tp.outcome_1 == "NAFLD") & _t1_tp.outcome_2.isin(_t1_fgset - {"NAFLD"})) |
                ((_t1_tp.outcome_2 == "NAFLD") & _t1_tp.outcome_1.isin(_t1_fgset - {"NAFLD"}))].null_z_correlation.astype(float)
assert len(_t1_fn) == 4 and len(_t1_ff) == len(_t1_f4) + len(_t1_fn)
assert _t1_fn.min() == _t1_ff.min() and _t1_f4.max() == _t1_ff.max() and _t1_fn.max() < _t1_f4.min()
add("L317", "Table 1: overlap", "Null-z correlation (cross-trait LDSC intercept) between NAFLD and each of CAD, HF, AF and T2D in FinnGen: minimum / maximum over 4 pairs (the other 6 FinnGen pairs are L206)",
    [float(_t1_fn.min()), float(_t1_fn.max())], "%.3f to %.3f" % (_t1_fn.min(), _t1_fn.max()),
    "V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_supp_table_a_intercepts.tsv", "null_z_correlation, pairs NAFLD x {CAD,HF,AF,T2D}", "[R]")

# ---------------------------------------------------------------- L318-L319: LD references
_t1_fams = [_t1_os.path.join(CKM, "paper 4/independent_build/data/g1000_eas_merged.fam"),
            _t1_os.path.join(CKM, "lhcmr_eas_refs/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.6.fam")]
_t1_nf = [sum(1 for l_ in _t1_io.open(p_) if l_.strip()) for p_ in _t1_fams]
assert _t1_nf[0] == _t1_nf[1], _t1_nf
add("L318", "Table S2: LD references", "1000 Genomes East Asian LD panel: individuals (.fam rows; the panel of the ancestry-matched analysis and of the r2 in L269)", _t1_nf[0], "%d" % _t1_nf[0],
    "paper 4/independent_build/data/g1000_eas_merged.fam (eas__census_common.G1K_FAM); lhcmr_eas_refs/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.6.fam", "row count", "[D]")
_t1_rd = _t1_io.open(_t1_os.path.join(_t1_D, "ldscore/lhcMR_README.md"), encoding="utf-8").read()
_t1_m = _t1_re.search(r"computed based on ([\d,]+) individuals from the UK10K study", _t1_rd)
assert _t1_m
_t1_uk = int(_t1_m.group(1).replace(",", ""))
add("L319", "Table S2: LD references", ("LD scores for the cross-trait LDSC intercepts (" + _P31R_EXT + "/D/data/ldscore/LDscores_filtered.csv): UK10K individuals they were computed from"), _t1_uk, "{:,}".format(_t1_uk),
    (_P31R_EXT + "/D/data/ldscore/lhcMR_README.md (lhcMR documentation, Darrous 2021 PMID 34907193)"), "README: 'LD scores and regression weights were computed based on 3,781 individuals from the UK10K study'", "[D]")

# ---------------------------------------------------------------- L320-L324: intermediate traits (frozen gwasinfo of the analysis)
_t1_gi = pd.read_csv(_t1_os.path.join(V3, "paper3_deposit_request_v4/files/data/intermediate_gwas/frozen/intermediate_gwasinfo.csv")).set_index("id")
# BMI: the frozen gwasinfo records 339,224 (Locke 2015 all-ancestry total), but the run that built the intermediate extracts read the
# European-descent file (41_3_freeze.log: 'GIANT_Locke2015_EUR322k'); Locke 2015 full text: 322,154 European + 17,072 non-European = 339,224.
_t1_frz = _t1_io.open(_t1_os.path.join(V3, "paper3_deposit_request_v4/scripts/reports/Paper3_Program/run_logs/RERUN_20260915/41_3_freeze.log"), encoding="utf-8", errors="replace").read()
assert _t1_re.search(r"BMI\s+GIANT_Locke2015_EUR322k", _t1_frz), "the BMI extract is no longer the European 322k file"
assert int(_t1_gi.loc["GIANT_Locke2015_BMI", "sample_size"]) == 339224   # Locke 2015: 'up to 339,224 individuals' (all ancestries)
add("L320", "Table S2: intermediate traits", "Intermediate trait BMI (GIANT, Locke 2015): European-descent N (the file the analysis read) / non-European / all-ancestry total (the value in the frozen gwasinfo)",
    [322154, 17072, 339224], "322,154 European (all ancestries up to 339,224)",
    "PMID 25673413 (where stated): Locke 2015 full text, PMC4382211 ('322,154 individuals of European descent and 17,072 individuals of non-European descent'); file identity run_logs/RERUN_20260915/41_3_freeze.log",
    "41_3_freeze.log 'GIANT_Locke2015_EUR322k'; intermediate_gwasinfo.csv sample_size 339224", "[P][D]",
    note="The frozen gwasinfo labels the BMI source European with N 339,224, the all-ancestry total; the extract log names the European 322k file. Locke's two subgroup counts sum to 339,226, two more than its stated total.")
EXT_CHECKS["L320"] = "MATCH (cross-source, T1 builder 2026-10-06): Locke 2015 Extended Data Table 1 (European joint GWAS+Metabochip 322,154; all ancestries 339,224) and the freeze-log file label EUR322k"
_t1_im = [("L322", "GLGC2013_LDL", 188577, "PMID 24097068 (abstract: '188,577 individuals')"),
          ("L323", "MAGIC2017_HbA1c", 123665, None),
          ("L324", "ieu-b-35", 204402, "PMID 30388399 (abstract: '204,402 European individuals')")]
for _t1_id, _t1_k, _t1_n, _t1_x in _t1_im:
    _t1_v = int(_t1_gi.loc[_t1_k, "sample_size"])
    assert _t1_v == _t1_n, (_t1_k, _t1_v)
    add(_t1_id, "Table S2: intermediate traits", "Intermediate trait %s (%s): N" % (_t1_k, _t1_gi.loc[_t1_k, "trait"]), _t1_v, "{:,}".format(_t1_v),
        "paper3_deposit_request_v4/files/data/intermediate_gwas/frozen/intermediate_gwasinfo.csv", "sample_size", "[D]",
        note=("" if _t1_x else "GWAS Catalog for PMID 28898252 lists the European HbA1c analysis as 88,355 (GCST007954; all four ancestries GCST004903); the paper's text states only 'up to 159,940' across four ancestries. Neither the 123,665 nor the European label of the analysis record (intermediate_gwasinfo.csv, PMID NA) is stated in the publication, and the file is not on this machine. Printed with [N TO CONFIRM] and [ANCESTRY TO CONFIRM]."))
    if _t1_x:
        EXT_CHECKS[_t1_id] = "MATCH (cross-source, T1 builder 2026-10-06): " + _t1_x
assert int(_t1_gi.loc["GLGC2013_TG", "sample_size"]) == 188577
_t1_sb = int(_t1_gi.loc["SBP_Evangelou2018_local", "sample_size"])
assert _t1_sb == 458577 + 299024 and "UKB-INCLUSIVE" in _t1_gi.loc["SBP_Evangelou2018_local", "note"]
add("L321", "Table S2: intermediate traits", "Systolic blood pressure (Evangelou 2018 discovery: UK Biobank + ICBP): N / UK Biobank / ICBP", [_t1_sb, 458577, 299024],
    "757,601 (UK Biobank 458,577; ICBP 299,024)", "intermediate_gwasinfo.csv sample_size; PMID 30224653 (where stated): Evangelou 2018 full text, PMC6284793 (UKB 458,577; ICBP 299,024; discovery up to 757,601)",
    "sample_size; Online Methods", "[D][P]")
EXT_CHECKS["L321"] = "MATCH (cross-source, T1 builder 2026-10-06): GWAS Catalog GCST006624 / GCST006628 initial sample '757,601 European ancestry individuals'"

# ---------------------------------------------------------------- L325-L327: kidney marker, metabolites, CKB pQTL (publications)
add("L325", "Table S2: kidney markers", "Cystatin-C-based eGFR (Stanzick 2021): N / UK Biobank / CKDGen", [460826, 436765, 24061], "460,826 (UK Biobank 436,765; CKDGen 24,061)",
    "PMID 34272381 (where stated): Stanzick 2021, abstract and PMC8285412 full text", "Results: eGFRcys in UKB (n = 436,765) meta-analysed with CKDGen (n = 24,061)", "[P]")
assert 436765 + 24061 == 460826
add("L326", "Table S2: metabolites", "CLSA plasma metabolite GWAS (Chen 2023): unrelated European individuals", 8299, "8,299",
    "PMID 36635386 (where stated): Chen 2023, PMC7614162 full text", "Methods: Study cohort", "[P]")
add("L327", "Table S2: ancestry-matched exposure", "China Kadoorie Biobank participants with Olink Explore (and SomaScan) plasma proteomics (Wang 2025)", 3976, "3,976",
    "PMID 39984443 (where stated): Wang 2025, abstract", "abstract", "[P]",
    note="Cohort size of the cited study, not the N of the summary statistics analysed (that is L340). CKB PheWeb (pheweb.ckbiobank.org/about) attributes its Olink protein GWAS to Pozarickij et al. 2025 (medRxiv); Mazidi 2026 (PMID 42524807) describes the CKB Olink case-cohort as 1,976 cases and 2,001 subcohort controls.")
# L340: the N of the CKB Olink summary statistics the ancestry-matched analysis read (CKB PheWeb phenotype manifest, md5-pinned in
# the census builder), over the proteins of that analysis (results_census.tsv, status ok = the 37 proteins of L216).
_t1_ph = _t1_os.path.join(CKM, "paper 3.eas/analysis/data/ckb_meta/ckb_phenotypes.json")
import hashlib as _t1_hl
assert _t1_hl.md5(open(_t1_ph, "rb").read()).hexdigest() == "254bfd16ae7e0ab6617f5582d79f5e8a", "CKB phenotype manifest differs from the pinned file (eas__census_common.PHEN_MD5)"
_t1_phd = {d_["phenocode"]: d_ for d_ in json.load(_t1_io.open(_t1_ph, encoding="utf-8"))}
_t1_rc = pd.read_csv(_t1_os.path.join(V3, "V4_local/staging/d2-eas-census/results_census.tsv"), sep="\t")
_t1_cg = sorted(set(_t1_rc[_t1_rc.status == "ok"].gene))
assert len(_t1_cg) == _t1_row("L216")[1] and len(_t1_rc) == _t1_row("L216")[0], "protein / pair count differs from L216"
_t1_cn = pd.Series([_t1_phd[g_]["num_samples"] for g_ in _t1_cg]).astype(int)
assert set(_t1_cn) == {3968, 3974} and _t1_cn.max() < 3976
add("L340", "Table S2: ancestry-matched exposure", "CKB PheWeb Olink summary statistics: per-protein N (manifest num_samples) over the 37 proteins of the ancestry-matched analysis: the two values / proteins with each",
    [3968, 3974, int((_t1_cn == 3968).sum()), int((_t1_cn == 3974).sum())], "3,968 or 3,974 per protein",
    "paper 3.eas/analysis/data/ckb_meta/ckb_phenotypes.json (md5 254bfd16, eas__census_common.PHEN) x V3/V4_local/staging/d2-eas-census/results_census.tsv (status ok)",
    "num_samples", "[D]", note="No protein in the manifest has N 3,976; the manifest gives 3,968 (1,455 traits) or 3,974 (1,463 traits).")

# ---------------------------------------------------------------- L328-L334: East Asian outcome-side files (European instrument)
_t1_ex = pd.read_csv(_t1_os.path.join(V3, "V4/source_data/supp_tables/SuppTable_EAS_axis.csv"))


def _t1_study(src):
    _s = _t1_ex[_t1_ex.eas_source == src]
    assert len(_s), src
    _o = []
    for c_ in ("n_cases_study", "n_controls_study", "n_total_study"):
        _u = sorted(set(_s[c_].dropna().astype(int)))
        assert len(_u) <= 1, (src, c_, _u)
        _o.append(_u[0] if _u else None)
    return _o


_t1_es = [("L328", "T2D", "T2DGGI_EAS", "GWAS Catalog GCST90492734 (the multi-ancestry T2DGGI study, PMID 38374256): East Asian stratum 88,109 cases, 339,395 controls"),
          ("L330", "AF", "BBJ_AF_2022", "Miyazawa 2023 abstract (PMID 36653681): 9,826 cases among 150,272 individuals"),
          ("L331", "CAD", "PHEWEB_BBJ_CAD", "pheweb.jp CAD: 32,512 cases, 146,214 controls (pheweb.jp attributes CAD to Kanai et al., medRxiv 2021; PMID 34594039 has no CAD study in the GWAS Catalog)"),
          ("L332", "CKD", "GCST90018602_CKD", "GWAS Catalog GCST90018602 (chronic renal failure, PMID 34594039): 2,117 / 174,345"),
          ("L333", "Stroke", "BBJ_IS", "GWAS Catalog GCST90018644 (ischemic stroke, East Asian, PMID 34594039): 22,664 / 152,022; the local hum0197.v3 Sakaue-Kanai IS file reproduces p, |beta/bx| and SE/bx of the analysis table at all 5 resolved instruments and lacks the 2 unresolved ones (scratchpad t1_bbj_match.py, t1_bbj_pos.py)"),
          ("L334", "HF", "BBJ_CHF", "GWAS Catalog GCST90018586 (chronic heart failure, East Asian, PMID 34594039): 10,540 / 168,186; the local hum0197.v3 Sakaue-Kanai CHF file reproduces p (3 s.f.) and per-variant N 178,726 at all 6 resolved instruments and lacks the 2 unresolved ones (scratchpad t1_bbj_match.py)")]
_t1_esv = {}
for _t1_id, _t1_o, _t1_src, _t1_x in _t1_es:
    _t1_ca, _t1_co, _t1_tt = _t1_study(_t1_src)
    assert _t1_tt is None and _t1_ca and _t1_co
    _t1_esv[_t1_src] = (_t1_ca, _t1_co)
    add(_t1_id, "Table S2: East Asian outcome-side", "East Asian outcome-side file %s (%s): cases / controls / total" % (_t1_src, _t1_o), [_t1_ca, _t1_co, _t1_ca + _t1_co], _t1_cc(_t1_ca, _t1_co),
        "V3/V4/source_data/supp_tables/SuppTable_EAS_axis.csv", "n_cases_study, n_controls_study (one value per source)", "[D]",
        note=("File identified as the hum0197.v3 Sakaue-Kanai release (Sakaue 2021, old:23) by matching the analysis table against the local file at every resolved instrument; Methods still carries [SOURCE TO CONFIRM] (edit requested)." if _t1_o in ("Stroke", "HF")
              else "Sakaue 2021 (old:23) has no CAD study; the pheweb.jp CAD release is attributed to Kanai et al. (medRxiv 2021, not in PubMed). Printed with [SOURCE TO CONFIRM]." if _t1_o == "CAD" else ""))
    EXT_CHECKS[_t1_id] = "MATCH (cross-source, T1 builder 2026-10-06): " + _t1_x
_t1_eN = _t1_study("EGFR_EAS_2P2M")
assert _t1_eN[0] is None and _t1_eN[1] is None and _t1_eN[2]
_t1_pv = _t1_ex[_t1_ex.eas_source == "EGFR_EAS_2P2M"].n_total_variant.dropna().astype(int)
assert _t1_pv.max() == _t1_eN[2]
add("L329", "Table S2: East Asian outcome-side", "East Asian outcome-side eGFR file EGFR_EAS_2P2M: study-level N (equal to the largest per-variant N at the instruments)", _t1_eN[2], "{:,}".format(_t1_eN[2]),
    "V3/V4/source_data/supp_tables/SuppTable_EAS_axis.csv", "n_total_study; n_total_variant", "[D]", note="File identity and scale are [SOURCE TO CONFIRM] / [SCALE TO CONFIRM] (source machine).")

# ---------------------------------------------------------------- L335-L339: ancestry-matched outcomes (resource registry of that analysis)
_t1_cp = _t1_os.path.join(V3, "V4_local/builders/eas__census_common.py")
_t1_tree = _t1_ast.parse(_t1_io.open(_t1_cp, encoding="utf-8").read())
_t1_reg = None
for _t1_nd in _t1_tree.body:
    if isinstance(_t1_nd, _t1_ast.Assign) and any(getattr(t_, "id", None) == "OUTCOMES" for t_ in _t1_nd.targets):
        _t1_reg = _t1_nd.value
assert isinstance(_t1_reg, _t1_ast.Dict)
_t1_res = {}
for _t1_k, _t1_v in zip(_t1_reg.keys, _t1_reg.values):
    _t1_kw = {k_.arg: k_.value for k_ in _t1_v.keywords}
    _t1_res[_t1_k.value] = {a_: (_t1_ast.literal_eval(_t1_kw[a_]) if a_ in _t1_kw else None) for a_ in ("label", "accession", "ncase", "nctrl", "ntot", "ckb_overlap")}
assert set(_t1_res) == {"T2D", "CAD", "eGFR", "CKD", "AF", "HF", "Stroke", "T2D_BBJ", "Stroke_BBJIS"}, set(_t1_res)
for _t1_a, _t1_b in (("AF", "BBJ_AF_2022"), ("CAD", "PHEWEB_BBJ_CAD"), ("CKD", "GCST90018602_CKD"), ("Stroke_BBJIS", "BBJ_IS")):
    assert (_t1_res[_t1_a]["ncase"], _t1_res[_t1_a]["nctrl"]) == _t1_esv[_t1_b], (_t1_a, _t1_b)
assert [k_ for k_, v_ in _t1_res.items() if v_["ckb_overlap"]] == ["T2D", "Stroke"]
_t1_am = [("L335", "T2D", "GWAS Catalog GCST010118: 77,418 / 356,122; Spracklen 2020 abstract (PMID 32499647)"),
          ("L336", "Stroke", "GWAS Catalog GCST90104544: 27,413 / 237,242; Mishra 2022 full text (PMC9524349) 27,413 East Asian any-stroke cases"),
          ("L338", "HF", "GWAS Catalog GCST90668009: 16,251 / 197,577; Enzan 2025 abstract (PMID 41184235) 16,251 cases among 213,828"),
          ("L339", "T2D_BBJ", "pheweb.jp T2D (Sakaue & Kanai 2021): 45,383 / 132,032; GWAS Catalog GCST90018706 (type 2 diabetes, East Asian, PMID 34594039): 45,383 / 132,032"),
          ("L342", "Stroke_BBJIS", "GWAS Catalog GCST90018644 (ischemic stroke, East Asian, PMID 34594039): 22,664 / 152,022")]
assert "hum0197.v3" in _t1_res["Stroke_BBJIS"]["accession"] and "hum0197.v3" in _t1_res["T2D_BBJ"]["accession"]
for _t1_id, _t1_k, _t1_x in _t1_am:
    _t1_r = _t1_res[_t1_k]
    add(_t1_id, "Table S2: ancestry-matched", "Ancestry-matched outcome %s, %s (%s; CKB overlap %s): cases / controls / total" % (_t1_k, _t1_r["label"], _t1_r["accession"], _t1_r["ckb_overlap"]),
        [_t1_r["ncase"], _t1_r["nctrl"], _t1_r["ncase"] + _t1_r["nctrl"]], _t1_cc(_t1_r["ncase"], _t1_r["nctrl"]),
        "V3/V4_local/builders/eas__census_common.py (OUTCOMES registry)", "%s: ncase, nctrl, ckb_overlap" % _t1_k, "[D]")
    EXT_CHECKS[_t1_id] = "MATCH (cross-source, T1 builder 2026-10-06): " + _t1_x
_t1_r = _t1_res["eGFR"]
assert _t1_r["ncase"] is None and _t1_r["accession"] == "TWB2_BBJ_eGFR_PMID39472450"
add("L337", "Table S2: ancestry-matched", "Ancestry-matched eGFR, %s (%s): N" % (_t1_r["label"], _t1_r["accession"]), _t1_r["ntot"], "{:,}".format(_t1_r["ntot"]),
    "V3/V4_local/builders/eas__census_common.py (OUTCOMES registry)", "eGFR: ntot", "[D]")
EXT_CHECKS["L337"] = "MATCH (cross-source, T1 builder 2026-10-06): Chen 2024 abstract (PMID 39472450) 'meta-analysis of multiple GWASs (n = 244,952)'"

# ---------------------------------------------------------------- L341: BUN (CKDGen European-ancestry release, Wuttke 2019)
# The file (CKDGen_BUN_EA_Wuttke2019.txt.gz) is on the analysis machine's E: drive; its read log (paper3_ckd_rerun_v3/partA/a1_BUN.log)
# records the header (n_total_sum present) and 8,358,026 data rows but no N summary. Wuttke 2019 states the European-ancestry BUN
# meta-analysis N and 8,358,347 European-ancestry BUN variants. The file maximum of n_total_sum has not been checked.
_t1_bl = json.load(_t1_io.open(_t1_os.path.join(V3, "paper3_ckd_rerun_v3/partA/a1_BUN.json"), encoding="utf-8"))
assert _t1_bl["path"].endswith("CKDGen_BUN_EA_Wuttke2019.txt.gz") and "n_total_sum" in _t1_bl["header"] and _t1_bl["data_rows"] == 8358026
add("L341", "Table S2: kidney markers", "Blood urea nitrogen, CKDGen European-ancestry meta-analysis (Wuttke 2019): N stated in the publication",
    243029, "243,029", "PMID 31152163 (where stated): Wuttke 2019 full text, PMC6698888 ('BUN GWAS meta-analysis of CKDGen European-ancestry studies (n = 243,029)'); file identity paper3_ckd_rerun_v3/partA/a1_BUN.json",
    "Results: Association of eGFR loci with BUN", "[P]",
    note="Printed with [N TO CONFIRM]: the maximum n_total_sum of the analysed file (E: drive) has not been checked. Wuttke's Methods read 'linear regression models to log(eGFR) and BUN'; the ln(BUN) scale in Methods rests on Stanzick 2021 (old:34), which combined CKDGen BUN results with UK Biobank log(BUN).")
