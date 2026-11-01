# -*- coding: utf-8 -*-
# Stage 3b, track C2 (Supplementary Methods part 2 and the Supplementary Note): ledger rows L425-L449.
# Executed by S1R1_build_ledger.py in ITS namespace (add, rel, pd, np, V3, CKM, S1, LF, LSRC, prim, EXT_CHECKS ...).
# Every name defined here starts with _s2_ so nothing in the builder's namespace is overwritten. Every value is read from
# its source file here and asserted; nothing is typed. EXT_CHECKS is set (end of file) only for rows that the Stage 3b C2
# verifier's separately written R scripts cover in full; the other rows stay unset.
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
import os as _p31r_os, sys as _p31r_sys, shutil as _p31r_shutil
_P31R_ARCHIVE = _p31r_os.path.normpath(_p31r_os.environ.get("P31R_ARCHIVE") or _p31r_os.path.join(_p31r_os.path.dirname(_p31r_os.path.abspath(__file__)), "..", "..", "..")).replace(chr(92), "/")
_P31R_V3 = _P31R_ARCHIVE + "/code"
_P31R_EXT = (_p31r_os.environ.get("P31R_EXTERNAL") or _P31R_ARCHIVE + "/external_inputs").replace(chr(92), "/")
_P31R_CKM = _P31R_EXT + "/CKM papers"
_P31R_RSCRIPT = _p31r_shutil.which("Rscript") or "Rscript"
_P31R_FONTDIR = _p31r_os.environ.get("P31R_FONTDIR") or _p31r_os.path.join(_p31r_os.environ.get("SystemRoot") or "/usr/share", "Fonts")
import io as _s2_io, os as _s2_os, re as _s2_re, json as _s2_json, subprocess as _s2_sp
from scipy.stats import norm as _s2_norm

_s2_P = lambda *a: _s2_os.path.join(V3, *a)
_s2_txt = lambda p: _s2_io.open(p, encoding="utf-8").read()
_s2_OC = _s2_P("R1_stage1", "S1R1_05b_outcome_checks")
_s2_TOPO = _s2_P("R1_stage1", "S1R1_02_topology_note")
_s2_DEP = _s2_P("paper3_deposit_request_v4")
_s2_EAS = _s2_os.path.join(CKM, "paper 3.eas", "analysis")
_s2_CEN = _s2_P("V4_local", "staging", "d2-eas-census")
_s2_FM = _s2_P("V4_local", "staging", "l3-finemap")

# ---------------------------------------------------------------------------------------------------- external outcome checks
# L425: cross-trait LD-score-regression settings, stated identically in the external-check run log and the Note's methods
_s2_olog = _s2_txt(_s2_os.path.join(_s2_OC, "S1R1_outcome_checks.log"))
_s2_tnote = _s2_txt(_s2_os.path.join(_s2_TOPO, "S1R1_topo_methods_note.md"))
assert "MAF >= 1%, chi2 < 80, 200-block jackknife" in _s2_olog
assert "MAF >= 1%; chi-square < 80" in _s2_tnote and "200-block jackknife" in _s2_tnote
add("L425", "SM: external outcome checks / Supplementary Note", "Cross-trait LD-score regression settings: minimum MAF; chi-square ceiling; jackknife blocks",
    [0.01, 80, 200], "MAF >= 1%; chi-square < 80; 200-block jackknife",
    "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_outcome_checks.log; V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_methods_note.md",
    "'MAF >= 1%, chi2 < 80, 200-block jackknife' (both files)", "[R]")

# L426: palindrome rule of the external checks (pre-specified; the post hoc script applies the same numbers)
_s2_spec = _s2_txt(_s2_os.path.join(_s2_OC, "S1R1_PRESPEC_outcome_checks.md"))
_s2_phpy = _s2_txt(_s2_P("R2_rewrite", "posthoc", "P31R_posthoc_check.py"))
assert "If MAF ≤ 0.42 and the frequencies disagree by more than 0.10" in _s2_spec and "RAD23B (MAF 0.45)" in _s2_spec
assert "maf > 0.42" in _s2_phpy and "abs(eafy - r.eaf) > 0.10" in _s2_phpy
add("L426", "SM: external outcome checks", "Palindromic instruments: frequency-checked when MAF <= 0.42, flagged if the effect-allele frequencies differ by more than 0.10; above 0.42 labelled forward-strand assumption",
    [0.42, 0.10], "MAF <= 0.42; |delta EAF| > 0.10", "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md section 4; V3/R2_rewrite/posthoc/P31R_posthoc_check.py",
    "harmonisation rule text", "[R]")

# L427: post hoc design counts (spec, erratum, run log)
_s2_plog = _s2_txt(_s2_P("R2_rewrite", "posthoc", "P31R_posthoc.log"))
_s2_pspec = _s2_txt(_s2_P("R2_rewrite", "posthoc", "P31R_POSTHOC_SPEC.md"))
_s2_perr = _s2_txt(_s2_P("R2_rewrite", "posthoc", "P31R_POSTHOC_SPEC_ERRATUM.md"))
_s2_m = _s2_re.search(r"pairs (\d+) distinct instruments (\d+)", _s2_plog)
_s2_n28, _s2_nv = int(_s2_m.group(1)), int(_s2_m.group(2))
_s2_m = _s2_re.search(r"outcome T2D \((\d+)\) or CAD \((\d+)\)", _s2_pspec)
_s2_nt, _s2_nc = int(_s2_m.group(1)), int(_s2_m.group(2))
_s2_m = _s2_re.search(r"NOT PRE-DECLARED: .*same sign (\d+) of (\d+)", _s2_plog)
_s2_n23 = int(_s2_m.group(2))
assert "the other 22 pairs" in _s2_pspec and "**What is correct: 23.**" in _s2_perr and _s2_n23 == 23
assert _s2_nt + _s2_nc == _s2_n28 == len(ph) and (ph.outcome == "T2D").sum() == _s2_nt and ph.instrument.nunique() == _s2_nv
add("L427", "SM: external outcome checks (post hoc)", "Post hoc check design: pairs (T2D / CAD); distinct lead variants; count of pairs not selected in advance as planned in the spec / as corrected by the erratum",
    [_s2_n28, _s2_nt, _s2_nc, _s2_nv, 22, _s2_n23], "%d pairs (%d T2D, %d CAD); %d variants; planned 22, corrected %d" % (_s2_n28, _s2_nt, _s2_nc, _s2_nv, _s2_n23),
    "V3/R2_rewrite/posthoc/P31R_POSTHOC_SPEC.md + _ERRATUM.md + P31R_posthoc.log", "spec text; 'pairs N distinct instruments M'; erratum", "[R]")

# L428: heterogeneity p values of the T2D pairs under DIAMANTE's two-significant-figure printing (half a unit of the last digit)
_s2_ins = pd.read_csv(_s2_P("V4", "source_data", "tables_from_rds", "phase1_instruments__ukb_cis_instruments.csv"))
_s2_dia = pd.read_csv((_P31R_EXT + "/D/data/paper31_r1_stage1/extracts/posthoc/DIAMANTE_posthoc.tsv"), sep="\t")
_s2_rng = {}
for _s2_i, _s2_r in ph[ph.outcome == "T2D"].iterrows():
    _s2_x = _s2_ins[(_s2_ins.gene == _s2_r.gene) & (_s2_ins.SNP == _s2_r.instrument)].iloc[0]
    _s2_y = _s2_dia[_s2_dia.rsid == _s2_r.instrument].iloc[0]
    _s2_sg = 1 if (_s2_y.ea == _s2_x.effect_allele and _s2_y.oa == _s2_x.other_allele) else -1
    assert _s2_sg == 1 or (_s2_y.ea == _s2_x.other_allele and _s2_y.oa == _s2_x.effect_allele)
    def _s2_hp(b, se, _x=_s2_x, _r=_s2_r, _sg=_s2_sg):
        return float(2 * _s2_norm.sf(abs((_sg * b / _x.beta - _r.disc_b) / np.sqrt((se / abs(_x.beta)) ** 2 + _r.disc_se ** 2))))
    _s2_h = lambda v: 0.5 * 10 ** (int(np.floor(np.log10(abs(v)))) - 1)
    assert float("%.2g" % _s2_y.b) == _s2_y.b and float("%.2g" % _s2_y.se) == _s2_y.se, (_s2_r.gene, _s2_y.b, _s2_y.se)   # printed to <= 2 s.f.
    _s2_grid = [_s2_hp(_s2_y.b + db, _s2_y.se + ds) for db in np.linspace(-_s2_h(_s2_y.b), _s2_h(_s2_y.b), 21)
                for ds in np.linspace(-_s2_h(_s2_y.se), _s2_h(_s2_y.se), 21)]
    assert abs(_s2_hp(_s2_y.b, _s2_y.se) - _s2_r.het_p) < 1e-12, _s2_r.gene
    _s2_rng[_s2_r.gene] = (float(_s2_r.het_p), min(_s2_grid), max(_s2_grid))
_s2_cad = ph[ph.outcome == "CAD"]
_s2_cmin = int((_s2_cad.het_p < 0.05).sum()) + sum(1 for v in _s2_rng.values() if v[2] < 0.05)
_s2_cmax = int((_s2_cad.het_p < 0.05).sum()) + sum(1 for v in _s2_rng.values() if v[1] < 0.05)
_s2_strad = sorted(g for g, v in _s2_rng.items() if v[1] < 0.05 <= v[2])
assert _s2_strad == ["IL17RB", "PAM"] and int((ph.het_p < 0.05).sum()) == 5 and (_s2_cmin, _s2_cmax) == (4, 6), (_s2_strad, _s2_cmin, _s2_cmax)
_s2_I, _s2_A = _s2_rng["IL17RB"], _s2_rng["PAM"]
add("L428", "SM: external outcome checks (post hoc)", "T2D heterogeneity p values whose range under DIAMANTE's two-significant-figure beta/SE crosses 0.05: IL17RB-T2D (p, min, max); PAM-T2D (p, min, max); resulting range of the count with heterogeneity p < 0.05 (of 28)",
    [_s2_I[0], _s2_I[1], _s2_I[2], _s2_A[0], _s2_A[1], _s2_A[2], _s2_cmin, _s2_cmax],
    "first pair %.3f (%.3f-%.3f); second pair %.3f (%.3f-%.3f); count %d-%d" % (_s2_I + _s2_A + (_s2_cmin, _s2_cmax)),
    (_P31R_EXT + "/D/data/paper31_r1_stage1/extracts/posthoc/DIAMANTE_posthoc.tsv (printed beta, SE) + V3/R2_rewrite/posthoc/P31R_posthoc_results.tsv + UKB-PPP instruments"),
    "beta and SE varied over +/- half a unit of the last printed digit (21 x 21 grid); same z-test as het_p", "[R]",
    note="Pairs: first = IL17RB-T2D, second = PAM-T2D (names kept out of value/display: the gate reads every digit there, including those inside gene names). "
         "Stage 3b (C2): the 3a status recorded IL17RB; PAM, a pair selected in advance, also crosses 0.05 within the rounding.")

# ---------------------------------------------------------------------------------------------------- East Asian quality control
_s2_ax = pd.read_csv(_s2_P("V4", "source_data", "supp_tables", "SuppTable_EAS_axis.csv"))
_s2_qc = _s2_ax.qc_reason.value_counts().to_dict()
_s2_pal = sorted(_s2_ax[_s2_ax.qc_reason == "PALINDROMIC_AMBIGUOUS"].gene + "-" + _s2_ax[_s2_ax.qc_reason == "PALINDROMIC_AMBIGUOUS"].outcome)
_s2_li = _s2_ax[_s2_ax.qc_reason == "LOW_INFO"]
assert (_s2_ax.qc_evaluable == False).sum() == 29 and _s2_pal == ["MINDY1-T2D", "PAM-T2D", "RAD23B-T2D"]
assert set(_s2_li.eas_source) == {"BBJ_AF_2022"} and set(_s2_li.source_qc_status) == {"INFO_MISSING"}
assert list(_s2_ax[_s2_ax.qc_reason == "OUTCOME_UNAVAILABLE"].outcome) == ["NAFLD"]
_s2_v429 = [_s2_qc["VARIANT_ABSENT"], _s2_qc["LOW_INFO"], _s2_qc["PALINDROMIC_AMBIGUOUS"], _s2_qc["OUTCOME_UNAVAILABLE"]]
assert sum(_s2_v429) == 29
add("L429", "SM: East Asian quality control", "Outcome-side pairs not evaluable, by reason: instrument absent / imputation information missing (all BBJ atrial fibrillation) / palindromic and strand-ambiguous / no East Asian outcome GWAS (APOE-NAFLD)",
    _s2_v429, "%d absent; %d INFO missing (BBJ AF); %d palindromic ambiguous; %d no outcome GWAS" % (_s2_v429[0], _s2_v429[1], _s2_v429[2], _s2_v429[3]),
    rel(_s2_P("V4", "source_data", "supp_tables", "SuppTable_EAS_axis.csv")), "qc_reason, eas_source, source_qc_status", "[D]",
    note="Palindromic ambiguous: %s (names in the note: the gate reads every digit in value/display)." % ", ".join(_s2_pal))

_s2_ev = _s2_ax[_s2_ax.qc_evaluable == True]
_s2_ex = _s2_ev[_s2_ev.eas_outcome_state_strict == "STRICT_EXCLUDED"]
_s2_exs = _s2_ex.eas_source.value_counts().to_dict()
_s2_lost = sorted(_s2_ev[(_s2_ev.eas_outcome_state == "EAS_COMPATIBLE") & (_s2_ev.eas_outcome_state_strict == "STRICT_EXCLUDED")].apply(lambda r: "%s-%s" % (r.gene, r.outcome), axis=1))
_s2_fam = _s2_ev[_s2_ev.sample_count_basis == "VARIANT_LEVEL"].eas_source.value_counts().to_dict()
assert set(_s2_ex.sample_count_basis) == {"STUDY_LEVEL_APPROX"} and _s2_fam == {"EGFR_EAS_2P2M": 14, "T2DGGI_EAS": 13}
_s2_v430 = [len(_s2_ex), _s2_exs["PHEWEB_BBJ_CAD"], _s2_exs["BBJ_CHF"], _s2_exs["GCST90018602_CKD"], _s2_exs["BBJ_IS"], len(_s2_lost), _s2_fam["T2DGGI_EAS"], _s2_fam["EGFR_EAS_2P2M"]]
assert _s2_v430[0] == sum(_s2_v430[1:5]) == 25 and _s2_v430[6] + _s2_v430[7] == 27
add("L430", "SM: East Asian quality control", "Strict sensitivity: evaluable pairs excluded for study-level sample counts (total; BBJ CAD / heart failure / CKD / ischaemic stroke); compatible pairs lost only by exclusion; strict family by outcome (T2D / eGFR)",
    _s2_v430, "%d excluded (CAD %d, HF %d, CKD %d, ischaemic stroke %d); %d compatible lost; family T2D %d + eGFR %d" % tuple(_s2_v430),
    rel(_s2_P("V4", "source_data", "supp_tables", "SuppTable_EAS_axis.csv")), "eas_outcome_state_strict, sample_count_basis, eas_source", "[D]",
    note="Compatible pairs lost by exclusion: %s (names in the note: the gate reads every digit in value/display)." % ", ".join(_s2_lost))

_s2_rc = pd.read_csv(_s2_os.path.join(_s2_CEN, "results_census.tsv"), sep="\t")
_s2_px = sorted(_s2_rc[_s2_rc.instr_is_proxy == True].apply(lambda r: "%s-%s" % (r.gene, r.outcome), axis=1))
_s2_tz = sorted(_s2_rc[_s2_rc.n_tied_top_z > 1].apply(lambda r: "%s-%s" % (r.gene, r.outcome), axis=1))
_s2_ccm = _s2_txt(_s2_P("V4_local", "builders", "eas__census_common.py"))
assert '"""User-decided tie rule: reported p ascending, |beta/sebeta| descending, position ascending."""' in _s2_ccm
assert _s2_px == ["ABO-Stroke", "ADAM15-eGFR", "GSTA3-eGFR"] and _s2_tz == ["ABO-Stroke", "ABO-T2D", "IL6R-CAD"]
add("L431", "SM: East Asian quality control", "Ancestry-matched analysis: instruments other than the CKB lead (lead not usable in the outcome file); instruments chosen among variants tied on p and |z|",
    [len(_s2_px), len(_s2_tz)], "%d not the CKB lead; %d tied" % (len(_s2_px), len(_s2_tz)),
    "V3/V4_local/staging/d2-eas-census/results_census.tsv; tie rule in V3/V4_local/builders/eas__census_common.py key_new", "instr_is_proxy, n_tied_top_z", "[D]",
    note="Not the CKB lead: %s; tied: %s (names in the note: the gate reads every digit in value/display)." % (", ".join(_s2_px), ", ".join(_s2_tz)))

_s2_v07 = _s2_txt(_s2_os.path.join(_s2_DEP, "scripts", "scripts", "_shared", "eas_axis", "07_verify_eas_axis.py"))
assert _s2_v07.count('q < 0.05 and dc == "SAME"') == 2 and _s2_v07.count('r["sample_count_basis"] == "VARIANT_LEVEL"') >= 2
add("L432", "SM: East Asian quality control", "Outcome-side classification threshold: Benjamini-Hochberg q < 0.05 within the evaluable family (primary) and within the per-variant-count family (strict)",
    0.05, "q < 0.05", "V3/paper3_deposit_request_v4/scripts/scripts/_shared/eas_axis/07_verify_eas_axis.py (state recompute)", "q < 0.05 and dc == SAME", "[D]")

# ---------------------------------------------------------------------------------------------------- the East Asian screen behind FGFR4-T2D
_s2_ci = pd.read_csv(_s2_os.path.join(_s2_EAS, "data", "ckb_meta", "ckb_cis_instruments.tsv"), sep="\t", comment="#")
assert _s2_ci.protein.is_unique
_s2_fs = LF["fgfr4_scan"]
# C2-V06: phase2_binary_scan.py keeps an instrument only if its rsid starts with 'rs' and beta/sebeta are present, then skips
# palindromic alleles; the four retained arms' result files give the proteins actually tested and their weakest F.
_s2_rsid = _s2_ci[_s2_ci.rsid.astype(str).str.startswith("rs") & _s2_ci.beta.notna() & _s2_ci.sebeta.notna()]
_s2_PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
_s2_np = _s2_rsid[[(str(a).upper(), str(r).upper()) not in _s2_PAL for a, r in zip(_s2_rsid.alt, _s2_rsid.ref)]]
_s2_arm = {o: pd.read_csv(_s2_os.path.join(_s2_EAS, "results", f), sep="\t") for o, f in
           [("T2D", "phase2_T2D_scan.tsv"), ("CAD", "phase2_CAD_scan.tsv"), ("Stroke", "phase2_Stroke_scan.tsv"), ("eGFR", "phase2_egfr_scan.tsv")]}
assert {o: len(t) for o, t in _s2_arm.items()} == _s2_fs["arms"], ({o: len(t) for o, t in _s2_arm.items()}, _s2_fs["arms"])
_s2_tst = pd.concat(_s2_arm.values()).drop_duplicates("protein")
assert len(_s2_tst) == 373 and set(_s2_tst.protein) <= set(_s2_np.protein)
_s2_Ft = float(_s2_ci[_s2_ci.protein.isin(_s2_tst.protein)].F.min())
assert abs(float(_s2_ci.F.min()) - 88.6) < 0.05 and _s2_Ft > float(_s2_ci.F.min())
add("L433", "SM: FGFR4 screen", "CKB cis-pQTL instruments of the East Asian screen (one per protein): number; with an rsID (the screen's requirement); of these, not palindromic; weakest F among the proteins tested in a retained arm",
    [len(_s2_ci), len(_s2_rsid), len(_s2_np), round(_s2_Ft, 1)],
    "%d instruments; %d with an rsID; %d not palindromic; weakest F among tested proteins %.1f" % (len(_s2_ci), len(_s2_rsid), len(_s2_np), _s2_Ft),
    "paper 3.eas/analysis/data/ckb_meta/ckb_cis_instruments.tsv; paper 3.eas/analysis/results/phase2_{T2D,CAD,Stroke,egfr}_scan.tsv; rule in src/phase2_binary_scan.py",
    "rows; rsid startswith 'rs'; PAL; F of proteins in any arm", "[D]",
    note="Stage 3b C2-V06: the weakest F of all 588 (88.6) belongs to an instrument that never entered and is deliberately not a value of this row, so a reversion to it fails the gate.")
assert sum(_s2_fs["arms"].values()) == _s2_fs["tests"] == 1128 and _s2_fs["cand"] == _s2_fs["strict"] + LF["scan_egfr_carried"]
_s2_v434 = [_s2_fs["arms"]["T2D"], _s2_fs["arms"]["CAD"], _s2_fs["arms"]["Stroke"], _s2_fs["arms"]["eGFR"], _s2_fs["strict"], LF["scan_egfr_carried"]]
add("L434", "SM: FGFR4 screen", "East Asian screen: proteins tested per outcome (T2D / CAD / stroke / eGFR); candidates meeting every filter; eGFR candidates carried to colocalization although they failed the filters",
    _s2_v434, "T2D %d, CAD %d, stroke %d, eGFR %d; %d met every filter; %d eGFR carried" % tuple(_s2_v434), "l2l4_facts: " + LSRC["fgfr4_scan"] + "; " + LSRC["scan_egfr_carried"],
    "fgfr4_scan arms / strict; scan_egfr_carried", "[G]")
_s2_bs = _s2_txt(_s2_os.path.join(_s2_EAS, "src", "phase2_binary_scan.py"))
assert 'r["strict"]=(r["eas_q"]<0.05 and r["eur_p"]>0.05 and common and r["het_p"]<0.05 and not r["pleio"])' in _s2_bs
assert 'common=(r["eas_af"] and 0.05<r["eas_af"]<0.95) and (r["eur_af"] is None or 0.05<r["eur_af"]<0.95)' in _s2_bs
add("L435", "SM: FGFR4 screen", "Screen filters: East Asian BH q < 0.05; European p > 0.05; heterogeneity p < 0.05; effect-allele frequency between 0.05 and 0.95 in both ancestries; outside the APOE, ABO and MHC regions",
    [0.05, 0.05, 0.05, 0.05, 0.95], "q < 0.05; European p > 0.05; heterogeneity p < 0.05; 0.05 < EAF < 0.95", "paper 3.eas/analysis/src/phase2_binary_scan.py", "strict definition", "[D]")
_s2_t2 = pd.read_csv(_s2_os.path.join(_s2_EAS, "results", "phase2_T2D_scan.tsv"), sep="\t").set_index("protein")
_s2_f4 = _s2_t2.loc["FGFR4"]
assert _s2_f4.rsid == "rs451643" and bool(_s2_f4.strict) and abs(float(_s2_f4.eur_p) - LF["fgfr4_eur_src"]["scan_eurp"]) < 1e-12
add("L436", "SM: FGFR4 screen", "FGFR4-T2D in the screen: East Asian q (AGEN); European p (DIAMANTE, the screen's comparator); East Asian vs European heterogeneity p",
    [float(_s2_f4.eas_q), float(_s2_f4.eur_p), float(_s2_f4.het_p)], "q = %.1e; European p = %.2f; heterogeneity p = %.1e" % (_s2_f4.eas_q, _s2_f4.eur_p, _s2_f4.het_p),
    "paper 3.eas/analysis/results/phase2_T2D_scan.tsv", "FGFR4 row: eas_q, eur_p, het_p", "[D]",
    note="The main text prints the European estimate on FinnGen R10 (L215, p = 0.27); the screen compared against DIAMANTE (p = 0.21).")
# L437: the FGFR4 tie at the CKB p floor and the AGEN T2D Wald ratios of the tied variants
_s2_fx = pd.read_csv(_s2_os.path.join(_s2_EAS, "data", "eas_exposure", "FGFR4_cis_eas.tsv"), sep="\t")
_s2_fx = _s2_fx[_s2_fx.rsid.astype(str).str.startswith("rs")]
_s2_pf = float(_s2_fx.pvalue.min())
_s2_td = _s2_fx[_s2_fx.pvalue == _s2_pf]
assert _s2_pf <= 1e-300 and len(_s2_td) == LF["fgfr4_tied"] == 4 and _s2_td.rsid.iloc[0] == "rs451643" and "rs351855" in set(_s2_td.rsid)
_s2_g1 = pd.read_csv(_s2_os.path.join(_s2_CEN, "slices", "G1K__FGFR4__pilot.tsv"), sep="\t").set_index("rsid")
_s2_ag = pd.read_csv(_s2_os.path.join(_s2_CEN, "slices", "T2D__FGFR4__pilot.tsv"), sep="\t")
_s2_wr = []
for _s2_i, _s2_r in _s2_td.iterrows():
    _s2_y = _s2_ag[_s2_ag.pos == int(_s2_g1.loc[_s2_r.rsid, "pos"])]
    assert len(_s2_y) == 1 and {_s2_y.ea.iloc[0], _s2_y.oa.iloc[0]} == {_s2_r.alt, _s2_r.ref}
    _s2_wr.append((1 if _s2_y.ea.iloc[0] == _s2_r.alt else -1) * float(_s2_y.beta.iloc[0]) / float(_s2_r.beta))
assert abs(_s2_wr[0] - LF["fgfr4"]["eas"]) < 5e-4 and max(_s2_wr) - min(_s2_wr) < 0.002
add("L437", "SM: FGFR4 screen / Methods", "FGFR4 East Asian instrument: rsID-bearing CKB variants tied at the file's p-value floor (incl. Gly388Arg rs351855); the floor; AGEN T2D Wald ratios of the tied variants (min, max)",
    [len(_s2_td), _s2_pf, min(_s2_wr), max(_s2_wr)], "%d tied at p = %.0e; Wald ratio %.4f to %.4f per SD" % (len(_s2_td), _s2_pf, min(_s2_wr), max(_s2_wr)),
    "paper 3.eas/analysis/data/eas_exposure/FGFR4_cis_eas.tsv; V3/V4_local/staging/d2-eas-census/slices/{G1K,T2D}__FGFR4__pilot.tsv", "pvalue floor; beta(AGEN)/beta(CKB)", "[R]",
    note="Supports the Methods sentence 'the tied variants gave almost identical estimates'.")

# ---------------------------------------------------------------------------------------------------- kidney benchmark
_s2_mc = pd.read_csv(_s2_P("v3.2", "_engine_root", "results", "repro_phase1_hg38", "pam_conditional", "egfr_mirage_classification.csv"))
_s2_cr = pd.read_csv(_s2_P("V4", "source_data", "tables_from_rds", "phase3_coloc__coloc_results_final.csv"))
_s2_eg = _s2_cr[_s2_cr.outcome_name == "eGFR"].copy()
_s2_eg["rule"] = (_s2_eg["PP.H4"] >= 0.80) & (_s2_eg["PP.H4"] >= 5 * _s2_eg["PP.H3"])
assert (_s2_eg.rule == _s2_eg.coloc_support).all() and _s2_eg.gene.is_unique and set(_s2_mc.gene) <= set(_s2_eg.gene)
_s2_pass = set(_s2_eg[_s2_eg.rule].gene)
_s2_conf = _s2_mc[_s2_mc["class"] == "concordant"]
_s2_notc = _s2_mc[_s2_mc["class"] != "concordant"]
_s2_ret_c = int(_s2_conf.gene.isin(_s2_pass).sum()); _s2_ret_n = sorted(set(_s2_notc.gene) & _s2_pass)
_s2_opp = sorted(_s2_mc[_s2_mc["class"] == "other"].gene)
assert _s2_ret_c == SF["renal"]["retained"] == 12 and len(_s2_notc) == 37
assert ((_s2_mc[_s2_mc["class"] == "mirage"].cys_p >= 0.05).all() and (_s2_mc[_s2_mc["class"] == "other"].cys_p < 0.05).all()
        and (np.sign(_s2_mc[_s2_mc["class"] == "other"].cys_b) != np.sign(_s2_mc[_s2_mc["class"] == "other"].crea_b)).all())
_s2_rt = _s2_txt(_s2_os.path.join(_s2_DEP, "scripts", "scripts", "_shared", "renal_triad", "01_renal_triad_mr.R"))
assert "(eaf_exp > 0.42 & eaf_exp < 0.58)" in _s2_rt and "(ofrq    > 0.42 & ofrq    < 0.58)" in _s2_rt and 'CKDGen_BUN_EA_Wuttke2019.txt.gz' in _s2_rt
# C2-V01: the benchmark's cystatin-C p is recomputed from the GWAS's Effect/StdErr (m[, p := 2 * pnorm(-abs(b / se))]), which the
# file prints to 4 decimal places. Back-compute both from the classification x the instrument beta, assert they sit on the 4-dp grid,
# and find the proteins whose p range (each value +/- half a unit of the 4th decimal) spans 0.05.
assert "m[, p := 2 * pnorm(-abs(b / se))]" in _s2_rt and '"Effect", "StdErr", "RSID")' in _s2_rt
_s2_ib = _s2_ins.drop_duplicates("gene").set_index("gene")
_s2_k = _s2_mc.set_index("gene").join(_s2_ib[["beta"]], how="left")
assert _s2_k.beta.notna().all()
_s2_k["eff"] = _s2_k.cys_b * _s2_k.beta
_s2_k["z"] = _s2_norm.isf(_s2_k.cys_p / 2)
assert (abs(_s2_k.eff * 1e4 - np.round(_s2_k.eff * 1e4)) < 1e-6).all()            # Effect printed to 4 dp
_s2_k = _s2_k[_s2_k.z > 0]                                                              # Effect 0.0000 -> p = 1 (no SE recoverable)
_s2_k["sef"] = np.round(abs(_s2_k.eff) / _s2_k.z, 4)
assert (abs(abs(_s2_k.eff) / _s2_k.z - _s2_k.sef) < 1e-9).all()                       # StdErr printed to 4 dp
_s2_hu = 5e-5
_s2_k["pmin"] = 2 * _s2_norm.sf((abs(_s2_k.eff) + _s2_hu) / (_s2_k.sef - _s2_hu))
_s2_k["pmax"] = 2 * _s2_norm.sf(np.maximum(abs(_s2_k.eff) - _s2_hu, 0) / (_s2_k.sef + _s2_hu))
_s2_x = _s2_k[(_s2_k.pmin < 0.05) & (_s2_k.pmax >= 0.05)]
_s2_same = np.sign(_s2_x.cys_b) == np.sign(_s2_x.crea_b)
_s2_cmn = len(_s2_conf) - int((_s2_x["class"] == "concordant").sum())
_s2_cmx = len(_s2_conf) + int(((_s2_x["class"] != "concordant") & _s2_same).sum())
_s2_omn = len(_s2_opp) - int((_s2_x["class"] == "other").sum())
_s2_z2 = sorted(_s2_k[abs(_s2_k.z - 2) < 1e-9].index)
assert sorted(_s2_x.index) == ["AMIGO1", "APOBR", "CASP9", "ENPP5", "MEGF9", "MSRA", "NFKB1", "OPLAH", "TNFAIP8L2"], sorted(_s2_x.index)
assert (_s2_cmn, _s2_cmx, _s2_omn) == (47, 54, 0) and set(_s2_opp) <= set(_s2_z2) and not (set(_s2_x.index) & _s2_pass)
add("L438", "SM: kidney benchmark", "Benchmark proteins not confirmed on cystatin C: cystatin p >= 0.05 / opposite sign with p < 0.05 (names) / total; of these, passing the colocalization rule (names); confirmed proteins passing the rule (recount of L070); palindrome band dropped (frequency 0.42-0.58); decimal places of the cystatin-C Effect/StdErr; proteins whose recomputed p spans 0.05 within that rounding; confirmed count range; opposite-sign count range; z of the two opposite-sign proteins",
    [int((_s2_mc["class"] == "mirage").sum()), len(_s2_opp), len(_s2_notc), len(_s2_ret_n), _s2_ret_c, len(_s2_conf), 0.42, 0.58,
     4, len(_s2_x), _s2_cmn, _s2_cmx, _s2_omn, len(_s2_opp), 2],
    "%d not significant; %d opposite; %d in all; %d of them pass the rule; confirmed: %d of %d pass; palindromes with frequency 0.42-0.58 dropped; Effect/SE at 4 dp: %d proteins span p = 0.05 (none passes the rule), confirmed %d-%d, opposite %d-%d (both at z = 2)" % (
        (_s2_mc["class"] == "mirage").sum(), len(_s2_opp), len(_s2_notc), len(_s2_ret_n), _s2_ret_c, len(_s2_conf),
        len(_s2_x), _s2_cmn, _s2_cmx, _s2_omn, len(_s2_opp)),
    "v3.2/_engine_root/results/repro_phase1_hg38/pam_conditional/egfr_mirage_classification.csv (+ UKB-PPP instrument beta); V4 phase3_coloc__coloc_results_final.csv (eGFR); renal_triad/01_renal_triad_mr.R",
    "class; PP.H4 >= 0.80 and >= 5 x PP.H3; back-computed Effect/StdErr +/- 0.00005", "[D]",
    note="Names (kept out of value/display, where the gate reads every digit, including those inside gene names): opposite %s; not confirmed but passing the rule %s; spanning p = 0.05 %s. " % (
        ", ".join(_s2_opp), ", ".join(_s2_ret_n), ", ".join(sorted(_s2_x.index))) +
         "The 12 of 52 is re-counted from the colocalization table directly (L070 counts it through the decoy-observed file, which omits three MHC genes that also fail). "
         "Stage 3b C2-V01: the benchmark recomputes p from 4-dp Effect/StdErr, whereas the triangulation (L208/L245) uses the GWAS's own P.value column; re-running with P.value needs the source machine.")
_s2_sw = pd.read_csv(_s2_P("v3.2", "_engine_root", "results", "repro_phase1_hg38", "pam_conditional", "egfr_mirage_threshold_sweep.csv"))
_s2_sw = _s2_sw.set_index("threshold_label")
_s2_b = [k for k in _s2_sw.index if k.startswith("Bonferroni")][0]
assert int(_s2_sw.loc["p<1e-3", "n_classifiable"]) == len(_s2_mc) and int(_s2_sw.loc["p<1e-3", "n_mirage"]) == int((_s2_mc["class"] == "mirage").sum())
_s2_v439 = [float(_s2_sw.loc[_s2_b, "threshold"])] + [int(_s2_sw.loc[k, c]) for k in [_s2_b, "p<1e-4", "p<1e-3", "p<0.05"] for c in ("n_mirage", "n_classifiable")]
add("L439", "SM: kidney benchmark", "Benchmark threshold sweep: proteins with cystatin-C p >= 0.05 of those classifiable at creatinine-eGFR p < Bonferroni / 1e-4 / 1e-3 / 0.05",
    _s2_v439, "Bonferroni (%.1e) %d of %d; p < 1e-4 %d of %d; p < 1e-3 %d of %d; p < 0.05 %d of %d" % tuple(_s2_v439),
    "v3.2/_engine_root/results/repro_phase1_hg38/pam_conditional/egfr_mirage_threshold_sweep.csv", "n_mirage, n_classifiable", "[D]")

# ---------------------------------------------------------------------------------------------------- gene attribution (48 settings)
_s2_rob = _s2_txt(_s2_P("V4_local", "builders", "l3fm__robustness.py"))
for _s2_s in ["GRID: 2 x 3 x 2 x 2 x 2 = 48 cells", "palindromic harmonised SNP with reference MAF <= 0.40", "GWAS MAF <= 0.40 and |GWAS frequency of the reference A1 allele",
              "match) - reference A1 frequency| <= 0.10", "estimate_s_rss s < 0.10 for BOTH fits", ">= 180 of 200 (90%) LD-reference bootstrap",
              "PROP iff |dev| <= 3.5 (v1)", "PROP iff |dev| / sd <= 1.96", "blocks only if |z_Q(A|B)| >= 5.4513",
              "DECLARATION (written 2026-09-27 BEFORE any class was computed under any cell of this grid"]:
    assert _s2_s in _s2_rob, _s2_s
assert LF["grid_cells"] == 48
add("L440", "SM: gene attribution", "48-setting grid: levels per dimension (palindromes 2, SuSiE evidence 3, LD stability 2, proportionality test 2, rival rule 2); palindrome MAF <= 0.40 and frequency difference <= 0.10; SuSiE LD-mismatch guard s < 0.10; bootstrap stability in >= 180 of 200 draws; |deviation| <= 3.5 or standardised <= 1.96; rival genome-wide |z| >= 5.4513",
    [2, 3, 2, 2, 2, 48, 0.40, 0.10, 0.10, 180, 200, 3.5, 1.96, 5.4513], "2 x 3 x 2 x 2 x 2 = 48; 0.40 / 0.10; s < 0.10; 180 of 200; 3.5 / 1.96; |z| >= 5.4513",
    "V3/V4_local/builders/l3fm__robustness.py (declaration docstring, md5-stamped before any class)", "GRID / level definitions", "[D]")
_s2_bc = pd.read_csv(_s2_P("V4_local", "source_data", "l3_metabolite", "bridge_class_summary.csv")).set_index(["protein", "metabolite"])
_s2_a1 = _s2_bc.loc[("ARG1", "arginine")]; _s2_d1 = _s2_bc.loc[("DPEP1", "cysteinylglycine")]
assert set(_s2_bc.annotation) == {"Ensembl 116 protein-coding genes (GRCh38)"} and (_s2_bc.n_cells == 48).all()
_s2_v441 = [int(_s2_a1.n_gene_attributed), int(_s2_a1.glist_n_gene_attributed), int(_s2_a1.ens116all_n_gene_attributed),
            int(_s2_d1.n_gene_attributed), int(_s2_d1.glist_n_gene_attributed), int(_s2_d1.ens116all_n_gene_attributed),
            LF["dpep1_halves"]["n"], int(_s2_re.search(r"Ensembl (\d+) protein-coding", _s2_bc.annotation.iloc[0]).group(1))]
assert _s2_v441[:6] == [48, 48, 48, 0, 16, 0] and _s2_v441[4] == LF["dpep1_attr"]["glist"] and LF["dpep1_halves"]["drop_gene"] == "VPS9D1"
assert _s2_d1.modal_class == "locus-level colocalization, gene uncertain" and _s2_a1.modal_class == "gene-attributed bridge"
add("L441", "SM: gene attribution", "Settings (of 48) giving gene attribution: ARG1-arginine under Ensembl 116 protein-coding / glist-hg19 / every-biotype annotation; DPEP1-cysteinylglycine under the same three; settings in which DPEP1's second signal is a VPS9D1 variant (palindromes dropped); Ensembl release",
    _s2_v441, "ARG1 %d / %d / %d; DPEP1 %d / %d / %d; VPS9D1 second signal in %d; Ensembl release %d" % tuple(_s2_v441),
    "V3/V4_local/source_data/l3_metabolite/bridge_class_summary.csv; l2l4_facts dpep1_halves", "n_gene_attributed, glist_n_gene_attributed, ens116all_n_gene_attributed", "[D][G]")
_s2_ps = {}
for _s2_k in ("susie_p1np", "susie_p2np"):
    _s2_t = pd.read_csv(_s2_os.path.join(_s2_FM, "DPEP1_16q24", _s2_k, "pair_summary.tsv"), sep="\t")
    _s2_t = _s2_t[(_s2_t.protein == "P__DPEP1") & (_s2_t.other == "M__GCST90200364")]
    assert len(_s2_t) == 1 and bool(_s2_t.susie_evaluable.iloc[0])
    _s2_ps[_s2_k] = (float(_s2_t.susie_max_H4.iloc[0]), float(_s2_t.abf_H4_R.iloc[0]),
                     _s2_json.load(_s2_io.open(_s2_os.path.join(_s2_FM, "DPEP1_16q24", _s2_k, "RUN_STAMP.json"), encoding="utf-8"))["ref_n"])
# C2-V08: s for the SNP set the paragraph describes (palindromes dropped, p1np), read from each locus's susie_p1np fit; the
# declaration's d2 quotes the same DPEP1 value (0.385). The all-SNP values (finemap_pairs s_hat_P) are 0.403 and 0.144.
_s2_st1 = lambda loc, tid: float(pd.read_csv(_s2_os.path.join(_s2_FM, loc, "susie_p1np", "susie_traits.tsv"), sep="\t").set_index("tid").loc[tid, "s_hat"])
_s2_sD, _s2_sA = _s2_st1("DPEP1_16q24", "P__DPEP1"), _s2_st1("ARG1_6q23", "P__ARG1")
assert "DPEP1 estimate_s_rss s = 0.385 (p1np)" in _s2_rob and abs(_s2_sD - 0.385) < 5e-4
assert _s2_ps["susie_p1np"][2] == 503 and _s2_ps["susie_p2np"][2] == 489 and _s2_ps["susie_p1np"][0] > 0.99 and _s2_ps["susie_p2np"][0] < 0.05
add("L442", "SM: gene attribution", "DPEP1-cysteinylglycine SuSiE-based PP.H4 with the 503-person panel and its 489-person subset (palindromes dropped); coloc.abf PP.H4 (LD-free); estimated LD mismatch s (palindromes dropped) for the DPEP1 protein fit and, for contrast, the ARG1 protein fit",
    [_s2_ps["susie_p1np"][0], _s2_ps["susie_p2np"][0], 503, 489, _s2_ps["susie_p1np"][1], _s2_sD, _s2_sA],
    "SuSiE PP.H4 %.3f (503) vs %.3f (489); coloc.abf %.3f; s %.3f (DPEP1) vs %.3f (ARG1), palindromes dropped" % (_s2_ps["susie_p1np"][0], _s2_ps["susie_p2np"][0], _s2_ps["susie_p1np"][1], _s2_sD, _s2_sA),
    "V3/V4_local/staging/l3-finemap/DPEP1_16q24/susie_p1np|susie_p2np/pair_summary.tsv + RUN_STAMP.json; {DPEP1_16q24,ARG1_6q23}/susie_p1np/susie_traits.tsv", "susie_max_H4, abf_H4_R, ref_n, s_hat", "[D]",
    note="Stage 3b C2-V08: s was the all-SNP value (0.403); the paragraph is about the palindromes-dropped set, whose s is 0.385 (DPEP1) and 0.139 (ARG1). Fix-round note 2026-10-10: the ARG1 0.139 is the p1np/503 per-fit s-hat, NOT the ARG1 s_max; the Rule 3(i) s_max for ARG1 is 0.148 (maximum over all seven fits, p1fr on the 489 panel; W4b S4.s.max.ARG1_*; the earlier four-fit maximum 0.144 is the labelled sensitivity).")

# ---------------------------------------------------------------------------------------------------- ARG1 signal-specific analysis
_s2_ar = _s2_txt(_s2_P("R1_stage1", "S1R1_01_arg1_triangle", "S1R1_arg1_triangle.R"))
_s2_m = _s2_re.search(r"NCAS <- (\d+); NCON <- (\d+)", _s2_ar)
_s2_ncas, _s2_ncon = int(_s2_m.group(1)), int(_s2_m.group(2))
assert _s2_ncas + _s2_ncon == int(SF["t1"]["n"]["T2D"]) == 400197
assert "keep <- r^2 < 0.8 & seq_along(r) != c" in _s2_ar and "pmin(af_ref, 1 - af_ref) <= 0.40 & pmin(afD, 1 - afD) <= 0.40 & abs(afD - af_ref) <= 0.10" in _s2_ar
assert "stopifnot(nrow(h) == 3067)" in _s2_ar and 'H    <- file.path(ROOT, "V3/V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz")' in _s2_ar
_s2_m = _s2_re.search(r"joined to harmonised panel: (\d+) of (\d+); after palindrome rule: (\d+)", tlog)
_s2_v443 = [_s2_ncas, _s2_ncon, int(_s2_m.group(2)), int(_s2_m.group(1)), int(_s2_m.group(3)), 0.8]
add("L443", "SM: ARG1", "ARG1 signal analysis inputs: FinnGen R10 T2D cases / controls; variants in the locus's harmonised fine-mapping panel (protein, metabolites, 1000 Genomes reference; not every variant in the window) / joined to FinnGen / kept after the palindrome rule; r^2 cut for conditional colocalization",
    _s2_v443, "%s cases / %s controls; %s variants, %s joined, %s kept; r2 < 0.8" % ("{:,}".format(_s2_ncas), "{:,}".format(_s2_ncon), "{:,}".format(_s2_v443[2]), "{:,}".format(_s2_v443[3]), "{:,}".format(_s2_v443[4])),
    "V3/R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_triangle.R (NCAS, NCON from R10_manifest) + .log", "NCAS/NCON; join line; cond()", "[R]",
    note="Stage 3b C2-V12: 3,067 is the harmonised panel (l3-finemap/ARG1_6q23/harmonised.tsv.gz, stopifnot nrow == 3067); the FinnGen window itself has 8,212 rows.")

# ---------------------------------------------------------------------------------------------------- Supplementary Note
_s2_tb = pd.read_csv(_s2_os.path.join(_s2_TOPO, "S1R1_topo_supp_table_b_pairs.tsv"), sep="\t")
_s2_ta = pd.read_csv(_s2_os.path.join(_s2_TOPO, "S1R1_topo_supp_table_a_intercepts.tsv"), sep="\t")
_s2_tc = pd.read_csv(_s2_os.path.join(_s2_TOPO, "S1R1_topo_supp_table_c_excess.tsv"), sep="\t")
_s2_zl = _s2_txt(_s2_os.path.join(_s2_TOPO, "S1R1_topo_step_outputs", "S1R1_topo_s01_zmatrix.log"))
_s2_m = _s2_re.search(r"Z matrix: (\d+) genes x 8 outcomes; non-missing cells: (\d+)", _s2_zl)
_s2_ng = int(_s2_m.group(1))
assert _s2_ng == prim.gene.nunique() and int(_s2_m.group(2)) == len(prim)
_s2_s01 = _s2_txt(_s2_os.path.join(_s2_TOPO, "S1R1_topo_s01_zmatrix.R")); _s2_s08 = _s2_txt(_s2_os.path.join(_s2_TOPO, "S1R1_topo_s08_random_table3.R"))
_s2_s07 = _s2_txt(_s2_os.path.join(_s2_TOPO, "S1R1_topo_s07_expected_table3.R"))
assert "boot_blocks <- function(block_id, B = 2000)" in _s2_s01 and "floor(G$start / 10e6)" in _s2_s01
assert "seq(100, .N, by = 100)" in _s2_s08 and "vapply(1:1000, function(b)" in _s2_s08 and "200 contiguous genomic blocks" in _s2_s08
assert "vapply(1:1000, function(b)" in _s2_s07 and "(EM, 0.05 grid)" in _s2_tnote
add("L444", "SM: Supplementary Note", "Correlation analysis settings: proteins with at least one testable outcome; protein-level bootstrap (resamples of 10-Mb regions); random variants every 100th; random-variant bootstrap (resamples of 200 contiguous blocks); simulations of the shared-noise model; prior grid step",
    [_s2_ng, 2000, 10, 100, 1000, 200, 1000, 0.05], "%s proteins; 2,000 resamples of 10-Mb regions; every 100th variant; 1,000 resamples of 200 blocks; 1,000 simulations; grid 0.05" % "{:,}".format(_s2_ng),
    "V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_s01_zmatrix.R/.log, s07, s08, S1R1_topo_methods_note.md", "Z matrix line; B; by = 100; 1:1000", "[R-topo]")
_s2_r = _s2_tb.dropna(subset=["random_snp_spearman"]).copy()
assert len(_s2_r) == 21
_s2_fg = {"CAD", "HF", "AF", "T2D", "NAFLD"}
_s2_ff = _s2_r[_s2_r.outcome_1.isin(_s2_fg) & _s2_r.outcome_2.isin(_s2_fg)]
_s2_fs_ = _s2_r[((_s2_r.outcome_1 == "Stroke") & _s2_r.outcome_2.isin(_s2_fg)) | ((_s2_r.outcome_2 == "Stroke") & _s2_r.outcome_1.isin(_s2_fg))]
_s2_above = lambda d: int((d.table3_recomputed > d.random_snp_spearman).sum())
_s2_out = lambda d: int((d.table3_recomputed > d.random_snp_ci_high).sum())
_s2_inside = int(((_s2_r.random_snp_spearman >= _s2_r.table3_ci_low) & (_s2_r.random_snp_spearman <= _s2_r.table3_ci_high)).sum())
_s2_xr = _s2_ta[(_s2_ta.source_1 != _s2_ta.source_2)]
_s2_x4 = _s2_tb[_s2_tb.outcome_1.isin(["CAD", "HF", "AF", "T2D"]) & _s2_tb.outcome_2.isin(["CAD", "HF", "AF", "T2D"])]
_s2_v445 = [len(_s2_ff), _s2_above(_s2_ff), len(_s2_fs_), _s2_above(_s2_fs_), _s2_out(_s2_fs_), _s2_inside, len(_s2_r),
            float(_s2_xr.null_z_correlation.abs().max()), float(_s2_x4.overlap_only_share_of_table3.min()), float(_s2_x4.overlap_only_share_of_table3.max()),
            float(_s2_x4.overlap_only_share_full_ldsc_covariance.min()), float(_s2_x4.overlap_only_share_full_ldsc_covariance.max())]
_s2_below = lambda d: int((d.table3_recomputed < d.random_snp_ci_low).sum())
_s2_v445 += [_s2_out(_s2_ff), _s2_below(_s2_ff)]                                   # C2-V04: FinnGen pairs on the interval measure
assert _s2_v445[:7] == [10, 9, 5, 4, 4, 20, 21] and _s2_v445[12:] == [5, 1] and len(_s2_x4) == 6, _s2_v445
assert list(_s2_ff[_s2_ff.table3_recomputed < _s2_ff.random_snp_ci_low][["outcome_1", "outcome_2"]].itertuples(index=False, name=None)) == [("CAD", "HF")]
add("L445", "SM: Supplementary Note", "Protein vs random-variant correlations: FinnGen-FinnGen pairs above the random value (of 10); FinnGen-stroke pairs above (of 5) and above the random-variant interval; random value inside the protein interval (of 21); largest |null-z correlation| between different resources; shared-noise-only share of the CAD/HF/AF/T2D correlations (normalised; full intercept covariance); FinnGen-FinnGen pairs above / below the random-variant interval",
    _s2_v445, "FinnGen %d of %d above; FinnGen-stroke %d of %d above (%d outside the random interval); random inside protein interval %d of %d; cross-resource |r0| <= %.3f; overlap share %.0f-%.0f%% (full covariance %.0f-%.0f%%); FinnGen pairs above the random interval %d, below it %d (CAD-HF)" % (
        _s2_v445[1], _s2_v445[0], _s2_v445[3], _s2_v445[2], _s2_v445[4], _s2_v445[5], _s2_v445[6], _s2_v445[7], 100 * _s2_v445[8], 100 * _s2_v445[9], 100 * _s2_v445[10], 100 * _s2_v445[11], _s2_v445[12], _s2_v445[13]),
    "V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_supp_table_a_intercepts.tsv, _b_pairs.tsv", "table3_recomputed vs random_snp_spearman and CIs; null_z_correlation; overlap_only_share*", "[R-topo]",
    note="Full-precision shares (L101 holds 3-dp copies that sit on a rounding boundary for '36%').")
_s2_c = _s2_tc
_s2_mrow = lambda thr: _s2_c[(_s2_c.analysis_set == "random SNP sets matched on LD score and MAF") & (_s2_c.threshold == thr)].iloc[0]
_s2_p0 = _s2_c[(_s2_c.analysis_set == "proteins (cis-pQTL MR)") & (_s2_c.unit == "protein") & (_s2_c.threshold == 0.01) & (_s2_c.null_model == "within-unit permutation")
               & (_s2_c.category == "all")].iloc[0]
_s2_lp = _s2_c[_s2_c.analysis_set == "proteins with an LD-scored lead instrument"]
assert int(_s2_p0.n_replicates) == 5000 and int(_s2_p0.observed_count) == 128 and _s2_lp.n_units.nunique() == 1
_s2_pp = _s2_c[(_s2_c.analysis_set == "proteins (cis-pQTL MR)") & (_s2_c.unit == "protein") & (_s2_c.threshold == 0.01) & (_s2_c.category == "all")
              & (_s2_c.observed_count == 128)]
_s2_corr = _s2_pp[~_s2_pp.finngen_noise_correlation.str.startswith(("none", "LDSC intercepts, diagonal only"))]
_s2_ind = _s2_pp[(_s2_pp.null_model == "signal + correlated noise model") & _s2_pp.finngen_noise_correlation.str.startswith("none")]
assert len(_s2_corr) == 11 and len(_s2_ind) == 1 and (_s2_corr.n_replicates == 20000).all()
_s2_lr = lambda thr: _s2_lp[_s2_lp.threshold == thr].iloc[0]
_s2_v446 = [0.01, 7, int(_s2_p0.n_replicates), int(_s2_lp.n_units.iloc[0]), int(_s2_mrow(0.01).n_sets),
            float(_s2_mrow(0.01).fold_min), float(_s2_mrow(0.01).fold_max), float(_s2_mrow(0.005).fold_min), float(_s2_mrow(0.005).fold_max),
            float(_s2_p0.null_mean), float(_s2_corr.null_mean.min()), float(_s2_corr.null_mean.max()), float(_s2_ind.null_mean.iloc[0]),
            float(_s2_lr(0.005).fold_excess), float(_s2_lr(0.005).hit_rate), float(_s2_lr(0.001).fold_excess), float(_s2_lr(0.001).hit_rate),
            float(_s2_mrow(0.01).fold_excess), float(_s2_mrow(0.01).hit_rate), float(_s2_mrow(0.005).fold_excess), float(_s2_mrow(0.005).hit_rate), 128]
add("L446", "SM: Supplementary Note", "Multi-outcome count (proteins with p < 0.01 for >= 2 of 7 outcome groups, eGFR and CKD as one): threshold; groups; within-unit permutations; LD-scored proteins; matched random-variant sets; matched-set fold range at p < 0.01 and 0.005; null means (within-unit permutation; models with correlated FinnGen noise, min and max; signal model with independent noise); proteins' fold and hit rate at p < 0.005 and 0.001; matched random sets' mean fold and hit rate at p < 0.01 and 0.005; observed count",
    _s2_v446, ("p < 0.01; 7 groups; %d permutations; %d proteins; %d sets; fold range %.2f-%.2f (p < 0.01), %.2f-%.2f (p < 0.005); null mean %.1f (permutation), %.1f-%.1f (correlated noise), %.1f (independent noise); "
               "proteins fold %.2f (hit %.1f%%), %.2f (%.1f%%); random fold %.2f (hit %.1f%%), %.2f (%.1f%%); observed %d") % tuple(
        _s2_v446[2:13] + [_s2_v446[13], 100 * _s2_v446[14], _s2_v446[15], 100 * _s2_v446[16], _s2_v446[17], 100 * _s2_v446[18], _s2_v446[19], 100 * _s2_v446[20], _s2_v446[21]]),
    "V3/R1_stage1/S1R1_02_topology_note/S1R1_topo_supp_table_c_excess.tsv", "n_replicates, n_units, n_sets, fold_min, fold_max, null_mean, fold_excess, hit_rate", "[R-topo]",
    note="Full precision. L102's 'overlap-aware' minimum (91.3) is the signal model with INDEPENDENT noise; the models with correlated FinnGen noise give 91.9-95.6.")

# L448: the ARG1 pre-stated decision rule's thresholds
_s2_ru = _s2_txt(_s2_P("R1_stage1", "S1R1_01_arg1_triangle", "S1R1_ARG1_PRESTATED_RULE.md"))
for _s2_s in ["ratio at B has the same sign as at A, p_B < 0.05, and A-vs-B heterogeneity p ≥ 0.05", "same sign, heterogeneity p ≥ 0.05, but p_B ≥ 0.05",
              "heterogeneity p < 0.05, or the B ratio has the opposite sign with p_B < 0.05", "Written 2026-10-02, before any FinnGen R10 T2D row in the ARG1 window was read"]:
    assert _s2_s in _s2_ru, _s2_s
add("L448", "SM: ARG1", "ARG1 pre-stated decision rule thresholds (p at signal B; heterogeneity p)", 0.05, "p < 0.05 / p >= 0.05",
    "V3/R1_stage1/S1R1_01_arg1_triangle/S1R1_ARG1_PRESTATED_RULE.md (decision rule)", "rule text", "[R]")

# ---------------------------------------------------------------------------------------------------- software versions (as recorded at run time)
_s2_rs = _P31R_RSCRIPT
_s2_env = _s2_sp.run([_s2_rs, "-e", 'cat(attr(readRDS(commandArgs(TRUE)[1]), "ckm_declared_partial_env"))',
                      _s2_os.path.join(_s2_DEP, "files", "results", "repro_phase1_hg38", "mr_all_results.rds")],
                     capture_output=True, text=True, timeout=120).stdout.strip()
assert _s2_env == "R=4.5.1; TwoSampleMR=0.7.1; data.table=1.17.8; coloc=5.2.3", _s2_env
_s2_rl = _s2_os.path.join(_s2_DEP, "scripts", "reports", "Paper3_Program", "run_logs", "RERUN_20260915")
assert _s2_txt(_s2_os.path.join(_s2_rl, "04_coloc.log")).startswith("This is coloc version 5.2.3")
assert "packages (same caveat -- present at sidecar generation): data.table 1.17.8 | susieR 0.14.2 | coloc 5.2.3" in _s2_txt(_s2_os.path.join(_s2_rl, "33_sensitivity_provenance.log"))
_s2_st = _s2_json.load(_s2_io.open(_s2_os.path.join(_s2_FM, "DPEP1_16q24", "susie_p1np", "RUN_STAMP.json"), encoding="utf-8"))
assert (_s2_st["susieR"], _s2_st["coloc"], _s2_st["R"][:15]) == ("0.14.2", "5.2.3", "R version 4.4.1")
assert "[LD local] EUR plink" in _s2_txt(_s2_os.path.join(_s2_rl, "07_susie.log")) and "PLINK v1.90b7.2 64-bit" in _s2_txt(_s2_os.path.join(_s2_rl, "07_susie.log"))
assert "**Reproducibility.** R 4.4.1 (data.table" in _s2_tnote
assert "R coloc 5.2.3 European-consistent priors: 0 call changes" in _s2_txt(_s2_os.path.join(_s2_CEN, "RUN_LOG.txt"))
assert _s2_plog.startswith("spec md5 573fd781") and "prespec md5 recorded c9079f2de1ffe75ce31f8a1d597d1fda; now c9079f2de1ffe75ce31f8a1d597d1fda" in _s2_olog
add("L447", "SM: software", "Software versions as recorded: MR screen (attribute written by the screen at run time); colocalization run log; fine-mapping run stamps; environment record (susieR for the main SuSiE runs, recorded after the runs); PLINK",
    {"screen": _s2_env, "coloc_log": "5.2.3", "finemap": "R 4.4.1, susieR 0.14.2, coloc 5.2.3", "env_record_susieR": "0.14.2", "plink_susie_ld": "v1.90b7.2", "note": "R 4.4.1", "ancestry_matched_sensitivity_coloc": "5.2.3"},
    "screen: R 4.5.1, TwoSampleMR 0.7.1, data.table 1.17.8, coloc 5.2.3; coloc log 5.2.3; fine-mapping R 4.4.1, susieR 0.14.2, coloc 5.2.3; susieR 0.14.2 (environment record, after the runs); PLINK v1.90b7.2 (SuSiE LD); Supplementary Note R 4.4.1; ancestry-matched sensitivity coloc 5.2.3",
    "paper3_deposit_request_v4/files/results/repro_phase1_hg38/mr_all_results.rds (attr ckm_declared_partial_env); RERUN_20260915/04_coloc.log, 33_sensitivity_provenance.log; 07_susie.log; l3-finemap RUN_STAMP.json; S1R1_topo_methods_note.md",
    "version strings", "[D]", note="Version strings are not numbers the prose gate can check (printed with a 'v' prefix or in its exempt 'R x.y.z' / 'coloc 5.2.3' forms); this row asserts them against the sources.")

# L449: the measured FinnGen overlap at full precision (L049a/L049b hold 4-dp copies; the SE 0.0155 sits on a rounding boundary)
_s2_ov = pd.read_csv(_s2_os.path.join(_s2_OC, "S1R1_overlap_ldsc.tsv"), sep="\t")
_s2_ovT = _s2_ov[_s2_ov.pair.str.contains("T2D")].iloc[0]; _s2_ovC = _s2_ov[_s2_ov.pair.str.contains("CHD")].iloc[0]
assert abs(_s2_ovT.null_z_corr - 0.0612) < 1e-4 and abs(_s2_ovC.null_z_corr - 0.0542) < 1e-4
add("L449", "SM: external outcome checks", "Measured FinnGen overlap at full precision: null-z correlation and SE, T2D (DIAMANTE) and CAD (CARDIoGRAMplusC4D)",
    [float(_s2_ovT.null_z_corr), float(_s2_ovT.null_z_corr_se), float(_s2_ovC.null_z_corr), float(_s2_ovC.null_z_corr_se)],
    "T2D %.3f (SE %.3f); CAD %.3f (SE %.3f)" % (_s2_ovT.null_z_corr, _s2_ovT.null_z_corr_se, _s2_ovC.null_z_corr, _s2_ovC.null_z_corr_se),
    "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_overlap_ldsc.tsv", "null_z_corr, null_z_corr_se", "[R]",
    note="Same values as L049a/L049b; printed SEs should be tagged here (full precision).")

# ---------------------------------------------------------------------------------------------------- independent checks (Stage 3b C2 verifier; different agent, R 4.4.1, separately written)
# Scripts copied to the V3 session scratchpad c2v/ (het.R, m.R, kid.R, kid2.R, k4.R, eas.R) and re-run by the C2 builder on 2026-10-06;
# their printed output (c2v/rerun_20261006.txt) was read against each row. Rows the scripts cover only in part (L437, L445 and the
# rest) stay unset.
_s2_XV = "Stage 3b C2 verifier, scratchpad c2v/%s (R 4.4.1; separately written; re-run 2026-10-06): %s"
EXT_CHECKS["L427"] = _s2_XV % ("het.R", "MATCH: 28 pairs, 16 T2D, 12 CAD, 27 variants; planned 22 / corrected 23 read from P31R_POSTHOC_SPEC.md and its erratum")
EXT_CHECKS["L428"] = _s2_XV % ("het.R, m.R", "MATCH: het p < 0.05 in 5 of 28; PAM 0.0297 (0.0152-0.0539), IL17RB 0.0551 (0.0488-0.0620); every other T2D pair on one side of 0.05. "
                               "het.R reads the half-unit from the printed string, which widens MINDY1 ('0.2'); m.R at two significant figures keeps it at 0.066-0.110, and the builder found no "
                               "trailing zero and no third significant figure in 3,000,000 rows of the DIAMANTE file, so '0.2' is 0.20")
EXT_CHECKS["L430"] = _s2_XV % ("eas.R", "MATCH: strict-excluded BBJ CAD 8, CHF 6, CKD 6, IS 5 (25); variant-level T2D 13, eGFR 14; compatible 12 primary vs 6 strict, the lost six CELSR2, FLT1, FURIN, IL6R, PECAM1 (CAD), TJAP1 (stroke); BH q reproduced within 9e-4")
EXT_CHECKS["L438"] = _s2_XV % ("kid.R, kid2.R, k4.R", "MATCH: 52 / 35 / 2 (APOBR, OPLAH); 12 of 52 confirmed and AK1, DNER, MINK1, SMPD3 of the 37 pass the rule; Effect/StdErr reconstructed on the 4-dp grid (max error 3e-18); "
                               "9 proteins span 0.05 (5 confirmed, 2 opposite, MEGF9 and TNFAIP8L2 same sign), none passes the rule, hence 47-54 and 0-2. The 0.42-0.58 band is asserted from 01_renal_triad_mr.R by this file")
EXT_CHECKS["L439"] = _s2_XV % ("kid.R", "MATCH: Bonferroni (3.21e-06) 7 of 39; p < 1e-4 17 of 58; p < 1e-3 35 of 89; p < 0.05 198 of 308")
