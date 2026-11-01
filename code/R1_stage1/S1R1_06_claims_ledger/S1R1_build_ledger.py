#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paper 3.1 rewrite -- CLAIMS LEDGER (Stage 1, item 6). Every value is DERIVED here from a named source; none is typed.

Sources:
  - the gated fact modules of the current build (v3.2/prose/scripts: v4_facts, l2l4_facts, sidecar_claims) -- tag [G];
  - Stage 1 outputs in V3/R1_stage1/ (own code this session) -- tag [R];
  - R1 skeptic re-derivations archived in V3/v3.2/HANDOFF_20261002_review_R1_evidence.zip -- tag [S];
  - V4 source tables read directly here -- tag [D].
Run:  python S1R1_build_ledger.py      -> S1R1_CLAIMS_LEDGER.csv (+ .md), and asserts the cross-row invariants.
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
import csv, io, json, os, sys, hashlib
import pandas as pd, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.dirname(os.path.dirname(HERE))
CKM = os.path.dirname(V3)
S1 = os.path.join(V3, "R1_stage1")
sys.path.insert(0, os.path.join(V3, "v3.2", "prose", "scripts"))
import v4_facts, l2l4_facts, sidecar_claims  # noqa: E402

VF, VSRC = v4_facts.facts()
LF, LSRC = l2l4_facts.facts()
SF = sidecar_claims.facts()
ROWS = []


def rel(p):
    return os.path.relpath(p, CKM).replace("\\", "/")


def add(cid, section, claim, value, display, source, key, tag, lb="", note=""):
    if isinstance(value, (np.integer,)):
        value = int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    ROWS.append(dict(id=cid, section=section, claim=claim, value=json.dumps(value, ensure_ascii=False, default=str),
                     display=display, source=source, key=key, tag=tag, load_bearing=lb, note=note))
    return value


def tsv(path):
    return pd.read_csv(path, sep="\t")


def f3(x):
    return "%.3f" % x


def sci(p):
    m, e = ("%.1e" % p).split("e")
    return "%s × 10^%d" % (m, int(e))


# ======================================================================================================== screen
MR = os.path.join(V3, "V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv")
m = pd.read_csv(MR)
prim = m[m.method.isin(["Wald ratio", "Inverse variance weighted"])]
n_prot = add("L001", "Methods/Results: screen", "Proteins with a cis-pQTL instrument (UKB-PPP, p < 10^-6)", SF["is"]["prot"],
             "{:,}".format(SF["is"]["prot"]), "sidecar_claims.facts()['is']", "prot", "[G]", "Y")
n_out = add("L002", "Methods", "Outcomes", 8, "8", "sidecar_claims.facts()['t1']['ids']", "len(ids)", "[G]", "Y")
assert len(SF["t1"]["ids"]) == n_out
fam = add("L003", "Methods/Results: screen", "Protein-outcome tests in the multiple-testing family (proteins x outcomes)",
          n_prot * n_out, "{:,}".format(n_prot * n_out), "L001 x L002", "product", "[D]", "Y")
tests = add("L004", "Results: screen", "Tests with an instrument available in the outcome GWAS (primary method)",
            int(prim[["gene", "outcome_name"]].drop_duplicates().shape[0]), "", rel(MR), "unique gene x outcome_name, method in {Wald, IVW}", "[D]", "Y")
assert tests == LF["f9"][("L1", "1")]
ROWS[-1]["display"] = "{:,}".format(tests)
single = add("L005", "Methods", "Primary tests that are single-instrument Wald ratios", int((prim.method == "Wald ratio").sum()),
             "", rel(MR), "method == 'Wald ratio'", "[D]")
ROWS[-1]["display"] = "{:,}".format(single)
nom = add("L006", "Results: screen", "Nominal associations (p < 0.05), all taken to colocalization", int((prim.pval < 0.05).sum()), "",
          rel(MR), "pval < 0.05 (primary)", "[D]", "Y")
ROWS[-1]["display"] = "{:,}".format(nom)
assert nom == VF["coloc_rows"]


def bh_count(p, M, q=0.05):
    p = np.sort(np.asarray(p)); k = np.arange(1, len(p) + 1); ok = p <= k / M * q
    return int(ok.nonzero()[0].max() + 1) if ok.any() else 0


bon = add("L007", "Results: screen", "Study-wide significant (Bonferroni, p < 0.05/15,584 = 3.2 × 10^-6)", int((prim.pval < 0.05 / fam).sum()), "",
          rel(MR), "pval < 0.05/L003", "[D]", "Y")
ROWS[-1]["display"] = str(bon)
assert bon == int(prim.bonferroni_sig.sum())
fdr = add("L008", "Results: screen", "FDR < 5% (Benjamini-Hochberg over the full 15,584-test family; the canonical flag)",
          bh_count(prim.pval.values, fam), "", rel(MR), "BH(pval, m = L003)", "[D]", "Y",
          note="The stored fdr_sig column counts 304 (BH over stored rows only); the 12 extra pairs are inert (R1 assessment rec. 8). Use 292 everywhere.")
ROWS[-1]["display"] = str(fdr)
add("L009", "Supplement", "Stored fdr_sig flag (BH over analysable rows) -- NOT for print", int(prim.fdr_sig.sum()), str(int(prim.fdr_sig.sum())),
    rel(MR), "fdr_sig", "[D]", note="Explains the 292 vs 304 difference if a reader recomputes.")

# ======================================================================================================== colocalization-supported pairs
add("L010", "Results: map", "Colocalization-supported pairs (dual rule PP.H4 >= 0.80 and PP.H4/PP.H3 >= 5, p12 = 1e-5)", VF["n_pairs"], "81",
    "v4_facts:" + VSRC["n_pairs"], "n_pairs", "[G]", "Y")
add("L011", "Results: map", "Proteins among them", VF["n_genes"], "62", "v4_facts:" + VSRC["n_genes"], "n_genes", "[G]", "Y")
add("L012", "Results: map", "Study-wide significant among the 81", VF["bonferroni_pairs"], "41", "v4_facts:" + VSRC["bonferroni_pairs"], "bonferroni_pairs", "[G]", "Y")
add("L013", "Results: map", "Pairs per outcome", VF["pairs_per_outcome"], "", "v4_facts:" + VSRC["pairs_per_outcome"], "pairs_per_outcome", "[G]")
t = VF["tiers"]
add("L014", "Results: map", "Evidence tiers Convergent / Supported / Flagged", t, "%d / %d / %d" % (t["Convergent"], t["Supported"], t["Flagged"]),
    "v4_facts:" + VSRC["tiers"], "tiers", "[G]", "Y")
assert sum(t.values()) == VF["n_pairs"]
add("L015", "Results: map", "Convergent pairs", ["PAM-T2D", "PCSK9-CAD", "RAD23B-T2D"], "PAM-T2D, PCSK9-CAD, RAD23B-T2D",
    "sidecar/v4 A11 evidence_tier == Convergent", "evidence_tier", "[G]", "Y")
pb = VF["prior_bands"]
add("L016", "Results: map", "Robust to the conservative prior (PP.H4 >= 0.80 at p12 = 1e-6) / attenuated / fragile", pb,
    "%d / %d / %d" % (pb["robust"], pb["atten"], pb["fragile"]), "v4_facts:" + VSRC["prior_bands"], "prior_bands", "[G]", "Y")
add("L017", "Results: map", "SuSiE-based colocalization: supportive / inconclusive / not evaluable / discordant", VF["susie_counts"],
    "47 / 30 / 3 / 1", "v4_facts:" + VSRC["susie_counts"], "susie_counts", "[G]")
add("L018", "Results: map", "Core-robust pairs (coloc + SuSiE + conservative prior jointly)", VF["core_robust"], "18", "v4_facts:" + VSRC["core_robust"], "core_robust", "[G]")
add("L019", "Results: map", "Steiger: forward / borderline", VF["steiger"], "80 / 1", "v4_facts:" + VSRC["steiger"], "steiger", "[G]",
    note="R1 MF7: drop 'borderline' wording for CCND2 if it is the borderline pair; check which.")
add("L020", "Results: map", "Pairs led by a missense variant (pairs / genes)", VF["missense_led"], "12 / 9", "v4_facts:" + VSRC["missense_led"], "missense_led", "[G]")
add("L021", "Results: map", "Flagged reasons", VF["flag_reasons"], "", "v4_facts:" + VSRC["flag_reasons"], "flag_reasons", "[G]")
add("L022", "Results: map", "Core-complete pairs supported only for lack of an evaluable external arm", VF["core_complete_no_evaluable_arm"], "12",
    "v4_facts:" + VSRC["core_complete_no_evaluable_arm"], "core_complete_no_evaluable_arm", "[G]", note="R1 rec. 2: 11 of them for lack of an arm (CCND2-T2D has a weakness).")

# ---- instrument strength (MF7) [S]
S4 = os.path.join(V3, "v3.2/HANDOFF_20261002_review_R1_evidence.zip") + "::exec/rev1-skeptic2/S4_output.txt"
ins = pd.read_csv(os.path.join(V3, "V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv"))
pe = pd.read_csv(os.path.join(V3, "v3.2/_engine_root/results/phase4_annotation/table_s_pair_evidence_tiered.csv"))
gmin = ins.groupby("gene").pval.min()
pe["minp"] = pe.gene.map(gmin)
lost = pe[pe.minp >= 5e-8][["gene", "outcome", "evidence_tier"]]
add("L023", "Methods/SM: instruments", "Pairs whose instrument does not reach p < 5 × 10^-8", len(lost), str(len(lost)),
    "V4 tables_from_rds/phase1_instruments__ukb_cis_instruments.csv x phase4 tiered table", "gene min pval >= 5e-8", "[D][S]", "Y")
assert len(lost) == VF["pairs_instrument_lost_5e8"]
_lp = sorted((lost.gene + "-" + lost.outcome).tolist())
add("L024", "SM: instruments", "The 11 pairs (gene-outcome) without a genome-wide-significant instrument; pairs per gene for BRAP and CCND2",
    {"pairs": _lp, "BRAP": int((lost.gene == "BRAP").sum()), "CCND2": int((lost.gene == "CCND2").sum())}, "AK1-eGFR; BRAP x5; CCND2 x3; FLT1-CAD; KCTD5-T2D",
    "as L023", "list; count by gene", "[D][S]", note="Counts per gene added 3a-ii so a figure can register them (idx 'BRAP' / 'CCND2').")
kept = pe[pe.minp < 5e-8]
add("L025", "SM: instruments", "Pairs / proteins remaining at p < 5 × 10^-8", [len(kept), kept.gene.nunique()], "%d / %d" % (len(kept), kept.gene.nunique()),
    "as L023", "complement", "[D][S]")
kt = kept.evidence_tier.value_counts().to_dict()
add("L026", "SM: instruments", "Tiers among the 70 (Convergent / Supported / Flagged)", kt, "%d / %d / %d" % (kt.get("Convergent", 0), kt.get("Supported", 0), kt.get("Flagged", 0)),
    "as L023", "evidence_tier", "[D][S]")
add("L027", "SM: instruments", "Strongest hub lost: BRAP lead p (borderline miss); F range of the lost instruments",
    [float(gmin["BRAP"]), VF["lost_5e8_F_range"]], "p 5.7 × 10^-8; F 24.3-29.5", "as L023 / v4_facts lost_5e8_F_range", "BRAP min pval", "[D][S]")
add("L028", "Results: map", "Proteins linked to >= 3 outcome units (eGFR+CKD = one unit), all 81 pairs", VF["genes_three_plus_units"], "APOE, BRAP (4), CCND2, FURIN",
    "v4_facts:" + VSRC["genes_three_plus_units"], "genes_three_plus_units", "[G]", note="BRAP and CCND2 drop at 5e-8; APOE and FURIN keep three units.")

# ======================================================================================================== exposure-side corroboration
tt = SF["tt"]
add("L030", "Results: corroboration", "deCODE-testable pairs", tt["dtest"], "48", "sidecar_claims.facts()['tt']", "dtest", "[G]", "Y")
add("L031", "Results: corroboration", "Met exposure-side corroboration criteria in deCODE (F >= 10 and concordant sign)", tt["dhigh"], "23",
    "sidecar_claims.facts()['tt']", "dhigh", "[G]", "Y", note="The p < 0.05 criterion never binds (R1 rev1-bio: all 34 F >= 10 assay rows have p <= 0.0137). deCODE reuses the discovery outcome, so it corroborates the exposure side only.")
add("L032", "Results: corroboration", "Direction-concordant in deCODE", tt["dconc"], "39 of 48", "sidecar_claims.facts()['tt']", "dconc", "[G]")
add("L033", "Results: corroboration", "ARIC-evaluable / concordant / discordant", [VF["aric_evaluable"], VF["aric_evaluable_concordant"], VF["aric_evaluable_discordant"]],
    "22 / 16 / 6", "v4_facts aric_*", "aric_evaluable*", "[G]")
xp = tt["xp"]
add("L034", "Results: corroboration", "Combined deCODE/ARIC class: concordant / discordant / conflicting", [xp["concordant"], xp["discordant"], xp["conflicting"]],
    "21 / 4 / 3", "sidecar_claims.facts()['tt']['xp']", "xp", "[G]", "Y")
EVZ = os.path.join(V3, "v3.2/HANDOFF_20261002_review_R1_evidence.zip")
add("L035", "Results: corroboration", "Of the 21 concordant: two-arm / one-arm", [15, 6], "15 / 6",
    "evidence zip exec/rev1-bio/adhoc_checks.txt [xplat] crosstab (concordant x n_evaluable)", "n_evaluable 2 / 1", "[1]",
    note="TYPED from the archived crosstab; re-derive from the tiered table (crossplatform_class x n_evaluable) before printing.")
add("L036", "Methods", "Both external pQTL resources are SomaScan (deCODE; ARIC)", "SomaScan", "SomaScan", "R1 assessment MF5 [1]", "platform", "[1]",
    note="Name the platform; cite the ARIC pQTL source.")

# ======================================================================================================== external outcome checks (NEW)
oc = tsv(os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_outcome_checks.tsv"))
ov = tsv(os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_overlap_ldsc.tsv"))
OCS = "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_outcome_checks.tsv"
for i, r in oc.iterrows():
    add("L04%d" % i, "Results: external outcome checks", "%s: Wald ratio in %s (per SD protein; log-odds); class" % (r.pair, r.resource),
        [float(r.wald), float(r.wald_se), float(r.p), r["class"], float(r.discovery_wald), float(r.discovery_se), float(r.het_p)],
        "%+.3f (SE %.3f), p = %s; %s" % (r.wald, r.wald_se, sci(r.p), r["class"]), OCS, "pair == %s" % r.pair, "[R]",
        "Y" if r.pair in ("ARG1-T2D", "PCSK9-CAD", "PAM-T2D", "FGFR4-T2D (EAS)") else "",
        note="Pre-specified (S1R1_PRESPEC_outcome_checks.md, md5 c9079f2d). European checks are UKB-free external consortium checks, not independent replication.")
add("L048", "Results: external outcome checks", "Pre-declared pairs concordant in the external resource (p < 0.05) / at p < 0.05/6",
    [int((oc["class"].iloc[:6] == "concordant").sum()), int(oc.bonferroni_6.iloc[:6].sum())], "6 of 6 / 4 of 6", OCS, "first six rows", "[R]", "Y")
for j, r in ov.iterrows():
    add("L049%s" % "ab"[j], "Methods/SM: external outcome checks", "Measured overlap with FinnGen: %s, null-z correlation (LDSC intercept)" % r.pair,
        [round(r.null_z_corr, 4), round(r.null_z_corr_se, 4)], "%.3f (SE %.3f)" % (r.null_z_corr, r.null_z_corr_se),
        "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_overlap_ldsc.tsv", r.pair, "[R]", "Y")

# ======================================================================================================== ARG1 case
rt = tsv(os.path.join(S1, "S1R1_01_arg1_triangle/S1R1_arg1_ratios.tsv"))
co = tsv(os.path.join(S1, "S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv"))
ARGR = "V3/R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_ratios.tsv"
ARGC = "V3/R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv"
a1 = LF["arg1_t2d"]; a1m = LF["arg1_t2d_mr"]
add("L050", "Results: ARG1", "ARG1-T2D MR (FinnGen): log-odds per SD ARG1, p", [a1m["b"], a1m["p"]], "%.3f, p = %s" % (a1m["b"], sci(a1m["p"])),
    "l2l4_facts:" + LSRC["arg1_t2d_mr"], "arg1_t2d_mr", "[G]", "Y")
add("L051", "Results: ARG1", "ARG1-T2D PP.H4 at p12 = 1e-5 / 1e-6", [a1["h4"], a1["h4_6"]], "%.3f / %.3f" % (a1["h4"], a1["h4_6"]),
    "l2l4_facts:" + LSRC["arg1_t2d"], "arg1_t2d", "[G]", "Y")
add("L052", "Results: ARG1", "ARG1-arginine PP.H4 (enzyme-substrate colocalization)", LF["bridge_top2"][1][2], "0.996",
    "l2l4_facts:" + LSRC["bridge_top2"], "bridge_top2", "[G]", "Y")
r = rt[(rt["mode"] == "marginal") & (rt.ratio == "M per P")].iloc[0]
add("L053", "Results: ARG1", "Arginine per SD ARG1 at signal A / signal B (beta-ratio), heterogeneity p", [float(r.A_ratio), float(r.B_ratio), float(r.het_p)],
    "%.2f and %.2f; p = %.2f" % (r.A_ratio, r.B_ratio, r.het_p), ARGR, "marginal, M per P", "[R][1]", "Y", note="Matches R1 rev1-bio arg1_beta_ratio.txt (-0.7255/-0.7158, p 0.947).")
r = rt[(rt["mode"] == "conditional") & (rt.ratio == "D per P")].iloc[0]
add("L054", "Results: ARG1", "T2D per SD ARG1 (FinnGen) at signal A|B / B|A; heterogeneity p (PRE-STATED PRIMARY)",
    [float(r.A_ratio), float(r.A_se), float(r.B_ratio), float(r.B_se), float(r.het_p)],
    "%.3f (SE %.3f) vs %.3f (SE %.3f); p = %.3f" % (r.A_ratio, r.A_se, r.B_ratio, r.B_se, r.het_p), ARGR, "conditional, D per P", "[R]", "Y",
    note="Pre-stated rule (S1R1_ARG1_PRESTATED_RULE.md; primary = conditional) -> WEAKENS. Marginal p = 0.053 would read 'consistent, not independently confirmed': the FinnGen verdict depends on conditioning (vR). DIAMANTE (L055) weakens on marginal estimates with B of opposite sign.")
oa = oc.set_index("pair")
ra, rb = oa.loc["ARG1-T2D"], oa.loc["ARG1-T2D signal B"]
hp_dia = float(open(os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_outcome_checks.log"), encoding="utf-8").read().split("heterogeneity z")[1].split(" p ")[1].split(";")[0])
add("L055", "Results: ARG1", "T2D per SD ARG1 in DIAMANTE (no UKB) at signal A / B; heterogeneity p",
    [float(ra.wald), float(rb.wald), hp_dia],
    "%.3f vs %+.3f; p = %.3f" % (ra.wald, rb.wald, hp_dia), "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_outcome_checks.log", "ARG1 in DIAMANTE line", "[R]", "Y",
    note="Same rule applied in DIAMANTE -> WEAKENS. Expected z at B under A's ratio -2.16 (observed +0.75).")
add("L056", "Results: ARG1", "Signal B frequency: Finnish vs non-Finnish European (EUR ref) vs DIAMANTE", [0.0073, 0.0219, round(rb.eaf_resource_aligned, 3)],
    "0.7% / 2.2% / 2.5%", "V3/R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_triangle.log + S1R1_outcome_checks.tsv", "Finnish A1 freq line; eaf_resource_aligned", "[R]",
    note="0.0073 and 0.0219 are read from the triangle log line; re-derive if printed.")
cA = co[(co.pair == "P-D") & (co.p12 == 1e-5)]
add("L057", "Results: ARG1", "Signal-specific coloc ARG1-T2D: A|B / B|A PP.H4 (p12 = 1e-5)",
    [round(cA[cA.analysis.str.startswith("signal A")].pp_h4.iloc[0], 3), round(cA[cA.analysis.str.startswith("signal B")].pp_h4.iloc[0], 3)],
    "0.954 / 0.009", ARGC, "pair P-D, p12 1e-5", "[R]")
add("L058", "Results: ARG1", "r between the two ARG1 signals (1000 Genomes EUR, n = 503)", -0.0445, "-0.04", "S1R1_arg1_triangle.log; = registered full_cond.tsv r_AB",
    "r(A,B)", "[R][G]")
add("L059", "Results: ARG1", "ARG1 -> BMI (SD of BMI per SD ARG1), p", [LF["arg1_bmi"]["b"], LF["arg1_bmi"]["p"]], "-0.090, p = 4.2 × 10^-5",
    "l2l4_facts:" + LSRC["arg1_bmi"], "arg1_bmi", "[G]")
add("L060", "Results: ARG1/ancestry", "ARG1-T2D in East Asians: ancestry-matched coloc call H3 (distinct signals); outcome-side estimate",
    [("ARG1–T2D" in LF["cen_h3_pairs"]), LF["arg1_outcomeside"]["beta"]], "distinct variants (PP.H3); outcome-side -0.118 (compatible)",
    "l2l4_facts cen_h3_pairs / arg1_outcomeside", "S23N", "[G]", note="The triangle did not recur in East Asians.")
add("L061", "Discussion: ARG1", "Prior ARG1-T2D risk reports (novelty status: strengthened prior hypothesis)", ["PMID 39870808", "PMID 30233283"],
    "Liu 2025 (Sci Rep; ARG1 among proteomic-MR mediators of BDH1 -> T2D); Fawad Ali Shah S 2018 (EXCLI J; rs2781666/rs2781665 and higher serum arginase-1 in T2D)",
    "PubMed metadata (2026-10-02)", "PMIDs", "[1][vR]", "Y",
    note="PMID 39885775 (Chen 2025, Biomark Med) is ARG1 variation and vascular complications WITHIN T2D: cite as related, not as a T2D-risk report (vR flag). PubMed indexes the 2018 first author's surname as 'Fawad Ali Shah' (forename Syed). Observational higher arginase-1 in T2D runs against the genetic direction.")

# ======================================================================================================== kidney loci
rn = SF["renal"]
add("L070", "Results: kidney", "Cystatin-C-concordant renal associations retained by single-signal colocalization", [rn["retained"], rn["concordant"]], "12 of 52",
    "sidecar_claims.facts()['renal']", "retained / concordant", "[G]", "Y", note="Benchmark retention, NOT a false-negative rate (MF3).")
add("L071", "Results: kidney", "UMOD association p with creatinine eGFR / cystatin eGFR; DPEP1 with creatinine eGFR", [rn["umod_crea"], rn["umod_cys"], rn["dpep1_crea"]],
    "7.2 × 10^-93 / 2.1 × 10^-76; 3.3 × 10^-10", "sidecar_claims.facts()['renal']", "umod_crea / umod_cys / dpep1_crea", "[G]")
add("L072", "Results: kidney", "DPEP1-cysteinylglycine PP.H4", LF["bridge_top2"][0][2], "0.998", "l2l4_facts:" + LSRC["bridge_top2"], "bridge_top2", "[G]", "Y")
add("L073", "Results: kidney", "DPEP1 locus: lead rsID, second-signal missense", [LF["dpep1_ann"]["lead"], LF["dpep1_ann"]["aa"]], "rs58290281; Glu351Gln",
    "l2l4_facts:" + LSRC["dpep1_ann"], "dpep1_ann", "[G]")
add("L074", "Results: kidney", "DPEP1 competing gene assignment (frequency-rule convention)", LF["dpep1_halves"]["freq_gene"], "CHMP1A",
    "l2l4_facts:" + LSRC["dpep1_halves"], "dpep1_halves", "[G]")
add("L075", "Results: kidney", "UMOD-eGFR East Asian colocalization PP.H4 (p12 1e-5 / 1e-6)", [LF["umod"]["h4"], LF["umod"]["h4_6"]], "0.992 / 0.929",
    "l2l4_facts:" + LSRC["umod"], "umod", "[G]")
rs_ = VF["renal_scrub"]
add("L076", "SM: kidney", "Renal pairs with cystatin-C / BUN concordance (marker triangulation)", [rs_["n"], rs_["cys_concordant"], rs_["with_bun"], rs_["bun_concordant"]],
    "cystatin 22 of 23; BUN 19 of 22", "v4_facts:" + VSRC["renal_scrub"], "renal_scrub", "[G]")
add("L077", "Methods/SM: eGFR scale", "eGFR effects are ln(eGFR) per SD protein; coloc used sdY = 0.1831", 0.1831, "ln(eGFR); sdY 0.1831",
    "evidence zip exec/rev1-skeptic1/work/f6_scale.log + rev1-calib/egfr_scale_check.json", "sdY", "[S]", "Y", note="Real SD units would be 5.46x larger. Relabel, do not rescale (USER default).")

mrE = prim[prim.outcome_name == "eGFR"][["gene", "b"]]
eE = pe[pe.outcome == "eGFR"].merge(mrE, on="gene")
flip = [len(eE), int((np.sign(eE.mr_beta) != np.sign(eE.b)).sum()), float((eE.mr_beta + eE.b).abs().max())]
assert flip[0] == flip[1]
add("L078", "Methods/legends: eGFR scale", "eGFR effects printed from the 81-pair table (S12, forest) are eGFR* = -ln(eGFR) per SD protein (sign-flipped in every eGFR pair)",
    flip, "%d of %d eGFR pairs sign-flipped (max |mr_beta + b| %.1e)" % (flip[1], flip[0], flip[2]),
    "phase4 tiered table mr_beta vs " + rel(MR) + " b", "outcome == eGFR", "[D][vR]", "Y",
    note="Label every printed eGFR* effect '-ln(eGFR) per SD of protein (positive = lower eGFR)'; the raw MR table carries +ln(eGFR) (vR flag).")

# ======================================================================================================== East Asian transfer
es = VF["eas_states"]
add("L080", "Results: ancestry", "East Asian outcome-side: evaluable / compatible / inconclusive / discordant / not evaluable",
    [es["EAS_COMPATIBLE"] + es["EAS_INCONCLUSIVE"], es["EAS_COMPATIBLE"], es["EAS_INCONCLUSIVE"], 0, es["EAS_NOT_EVALUABLE"]], "52 / 12 / 40 / 0 / 29",
    "v4_facts:" + VSRC["eas_states"], "eas_states", "[G]", "Y")
_ex = pd.read_csv(os.path.join(V3, "V4/source_data/supp_tables/SuppTable_EAS_axis.csv"))
_st = _ex.eas_outcome_state_strict.value_counts()
_sev = int(_st.get("EAS_COMPATIBLE", 0) + _st.get("EAS_INCONCLUSIVE", 0)); _sco = int(_st.get("EAS_COMPATIBLE", 0)); _sxc = int(_st.get("STRICT_EXCLUDED", 0))
_sv = _ex[_ex.eas_outcome_state_strict.isin(["EAS_COMPATIBLE", "EAS_INCONCLUSIVE"])]
assert (_sv.sample_count_basis == "VARIANT_LEVEL").all() and set(_sv.outcome) == {"T2D", "eGFR"}
assert (_ex[_ex.eas_outcome_state_strict == "STRICT_EXCLUDED"].sample_count_basis == "STUDY_LEVEL_APPROX").all()
assert [_sev, _sco, _sxc] == [27, 6, 25], [_sev, _sco, _sxc]
add("L081", "Results: ancestry/SM", "Strict sensitivity (outcome files reporting per-variant sample sizes: the T2D and eGFR files only): evaluable / compatible; excluded rows (all study-level counts)",
    [_sev, _sco, _sxc], "%d / %d" % (_sev, _sco), "V3/V4/source_data/supp_tables/SuppTable_EAS_axis.csv", "eas_outcome_state_strict; sample_count_basis", "[D]", "Y",
    note="Re-derived 3a-iii from the EAS axis table (was typed from a skeptic log); T2D %d + eGFR %d, all VARIANT_LEVEL. Lost compatibles: CELSR2, FLT1, FURIN, IL6R, PECAM1 (CAD) and TJAP1 (stroke)." % ((_sv.outcome == "T2D").sum(), (_sv.outcome == "eGFR").sum()))
add("L082", "Results: ancestry", "Ancestry-matched East Asian analysis: pairs / testable / colocalize / robust at 1e-6",
    [LF["cen_n"], LF["cen_test"], len(LF["cen_coloc"]), len(LF["cen_robust"])], "48 / 47 / 7 / 3", "l2l4_facts cen_*", "S23N", "[G]", "Y")
add("L083", "Results: ancestry", "Ancestry-matched colocalizations; robust subset", [LF["cen_coloc"], LF["cen_robust"]], "", "l2l4_facts cen_coloc / cen_robust", "S23N", "[G]")
add("L084", "Results: ancestry", "CKB-free colocalizations (default / conservative prior)", list(LF["ckbfree_totals"]), "8 / 4", "l2l4_facts:" + LSRC["ckbfree_totals"], "ckbfree_totals", "[G]", "Y",
    note="CKB-free swap also changes phenotype (any -> ischaemic stroke) and BBJ sits inside AGEN.")
tg = LF["tgfb1"]
add("L085", "Results: ancestry", "TGFB1-CAD: European / East Asian log-odds per SD; heterogeneity z; instruments 32 kb apart", [tg["eur"], tg["eas"], tg["hetz"], tg["dist_kb"]],
    "+0.164 / -0.371; z 8.08; 32 kb", "l2l4_facts:" + LSRC["tgfb1"], "tgfb1", "[G]", "Y", note="'opposite signs on instruments 32 kb apart', not 'sign-reversed'. European variant ~1.5% in CKB, absent from BBJ CAD [1].")
fg = LF["fgfr4"]
add("L086", "Results: ancestry", "FGFR4-T2D (East Asian only): log-odds per SD, p; PP.H4; CKB-free BBJ PP.H4", [fg["eas"], fg["easp"], fg["h4"], fg["ckbfree"]],
    "+0.048, p = 3.9 × 10^-12; 0.989; 0.991", "l2l4_facts:" + LSRC["fgfr4"], "fgfr4", "[G]")
add("L087", "Results: ancestry", "FGFR4 instrument in LD with Gly388Arg (r^2, 1000G EAS)", LF["fgfr4_g388r_r2"], "0.97", "l2l4_facts:" + LSRC["fgfr4_g388r_r2"], "fgfr4_g388r_r2", "[G]")
add("L088", "SM: ancestry", "Palindromic ancestry-matched colocalizations and their orientation check",
    LF["cen_pal_coloc"], "FGF5-eGFR, GALNT3-T2D: consistent with forward strand",
    "l2l4_facts cen_pal_coloc; FGF5 via BBJ creatinine (this session: rs16998073-A +0.026 creatinine vs -0.021 eGFR; T 0.30 BBJ vs 0.41 CKB); GALNT3 A 0.38 CKB / 0.34 AGEN / 0.31 BBJ (rev1-bio adhoc)",
    "cen_pal_coloc", "[G][R]", note="Fix the Figure 3 legend's 'six colocalizations in the same-sign quadrants' (it counts the palindromic pairs).")
s26 = LF["s26_q"]
S26F = os.path.join(V3, "V4_local/source_data/figures/S26_metabolite_coloc_deflation/S26_metabolite_coloc_deflation__quadrant__source.csv")
q26 = pd.read_csv(S26F)


def bhq(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); r = p[o] * n / np.arange(1, n + 1)
    q = np.empty(n); q[o] = np.minimum.accumulate(r[::-1])[::-1]; return np.minimum(q, 1)


q26["q_all"] = bhq(q26.p_het.values)
dv = q26[q26.quadrant.isin(["EUR-only", "EAS-only", "discordant"])].copy(); dv["q_in"] = bhq(dv.p_het.values)
v89 = [len(dv), int((dv.p_het < 0.05).sum()), int((dv.q_all < 0.05).sum()), int((dv.q_in < 0.05).sum()), len(q26)]
assert v89[0] == s26["EUR-only"] + s26["EAS-only"] + s26["discordant"]
add("L089", "SM: metabolites x ancestry", "Metabolite-ancestry cells detected in one ancestry only or discordant; nominal heterogeneity; BH over all 253 cells; BH within the 70",
    v89, "%d cells; %d nominal; %d after BH over all %d cells (%d if BH within the %d)" % (v89[0], v89[1], v89[2], v89[4], v89[3], v89[0]),
    rel(S26F), "quadrant, p_het", "[D][vR]", note="'different ancestry-specific detection patterns', not 'ancestry-divergent'. Name the BH family when printing (vR flag).")

# ======================================================================================================== decoy-locus stress test
dr = tsv(os.path.join(S1, "S1R1_03_decoy_method/S1R1_decoy_rates.tsv")).set_index("stratum")
DRS = "V3/R1_stage1/S1R1_03_decoy_method/S1R1_decoy_rates.tsv"
for cid, k, lab in [("L090", "non-eGFR (seven outcomes, CKD included)", "seven non-eGFR outcomes (CKD included)"), ("L091", "eGFR", "eGFR"), ("L092", "CKD (all decoy draws)", "CKD, all draws"),
                    ("L093", "CKD, selection-matched, p < 0.05 at the located instrument (sel_p; 154 contributing targets)", "CKD, selection-matched (154 targets contribute)")]:
    r = dr.loc[k]
    add(cid, "Methods/SM: decoy-locus stress test", "Decoy pass rate, %s (two-way cluster bootstrap 95%% interval)" % lab,
        [int(r.decoy_passes), int(r.draws), round(r.rate, 5), round(r.twoway_boot_lo, 5), round(r.twoway_boot_hi, 5)],
        "%d of %s (%.2f%%; %.2f-%.2f%%)" % (r.decoy_passes, "{:,}".format(int(r.draws)), 100 * r.rate, 100 * r.twoway_boot_lo, 100 * r.twoway_boot_hi),
        DRS, k, "[R][S]", "Y" if cid in ("L090", "L093") else "", note="Pass rates are NOT error rates; no specificity or FDR is printed (MF2).")
dl = open(os.path.join(S1, "S1R1_03_decoy_method/S1R1_decoy_rates.log"), encoding="utf-8").read()
add("L094", "SM: decoy-locus stress test", "CKD draws truncating the target / cutting a located instrument row; passes with it cut", [8160, 1891, 16],
    "50.4%; 11.7%; 16 of 50", "V3/R1_stage1/S1R1_03_decoy_method/S1R1_decoy_rates.log", "truncation lines", "[R][S]")
assert "8160 (50.4%)" in dl and "(idx > k): 1891 (11.7%)" in dl and "cut away: 16 of 50" in dl

# ======================================================================================================== sharing (one sentence) -- topology note (agent; pending)
TOPO = os.path.join(S1, "S1R1_02_topology_note")
tb = pd.read_csv(os.path.join(TOPO, "S1R1_topo_supp_table_b_pairs.tsv"), sep="	")
tc = pd.read_csv(os.path.join(TOPO, "S1R1_topo_supp_table_c_excess.tsv"), sep="	")
rp = pd.read_csv(os.path.join(TOPO, "S1R1_topo_reproduction.csv"))
assert len(rp) == 201 and set(rp.status) <= {"identical", "agrees within stated tolerance"}
TOPOS = "V3/R1_stage1/S1R1_02_topology_note/"
tr = tb.dropna(subset=["random_snp_spearman"])
r100 = float(np.corrcoef(tr.table3_manuscript, tr.random_snp_spearman)[0, 1])
assert tr.n_random_snps.nunique() == 1
add("L100", "Results: one sentence / SM topology note", "Cross-outcome correlations at the protein instruments (old Table 3) vs at random SNPs: r across outcome pairs; number of random SNPs",
    [round(r100, 4), len(tr), int(tr.n_random_snps.iloc[0])], "r = %.3f across %d outcome pairs; %s random SNPs" % (r100, len(tr), "{:,}".format(int(tr.n_random_snps.iloc[0]))), TOPOS + "S1R1_topo_supp_table_b_pairs.tsv", "table3_manuscript vs random_snp_spearman",
    "[S][R-topo]", "Y", note="Wording: 'largely mirrored'; FinnGen-stroke pairs sit slightly above random SNPs, so a small protein-specific excess cannot be excluded (topology agent; visible in its Figure a/b).")
cm4 = tb[tb.outcome_1.isin(["CAD", "HF", "AF", "T2D"]) & tb.outcome_2.isin(["CAD", "HF", "AF", "T2D"])]
add("L101", "SM topology note", "Overlap-only (shared-participant) share of each CAD-HF-AF-T2D correlation",
    [round(cm4.overlap_only_share_of_table3.min(), 3), round(cm4.overlap_only_share_of_table3.max(), 3), len(cm4)],
    "%.0f-%.0f%% (%d pairs)" % (100 * cm4.overlap_only_share_of_table3.min(), 100 * cm4.overlap_only_share_of_table3.max(), len(cm4)),
    TOPOS + "S1R1_topo_supp_table_b_pairs.tsv", "overlap_only_share_of_table3", "[S][R-topo]")
pr = tc[(tc.analysis_set == "proteins (cis-pQTL MR)") & (tc.unit == "protein") & (tc.threshold == 0.01) & (tc.observed_count == 128)]
wu = pr[pr.null_model == "within-unit permutation"].null_mean.iloc[0]
oa = pr[pr.null_model != "within-unit permutation"].null_mean
def one(aset, thr, col):
    v = tc[(tc.analysis_set == aset) & (tc.threshold == thr)][col]
    assert len(v) == 1, (aset, thr, len(v))
    return float(v.iloc[0])


PRT, RND = "proteins with an LD-scored lead instrument", "random SNP sets matched on LD score and MAF"
v102 = [128, round(float(wu), 2), round(float(oa.min()), 1), round(float(oa.max()), 1),
        round(one(PRT, 0.005, "fold_excess"), 2), round(one(RND, 0.01, "fold_excess"), 2),
        round(one(PRT, 0.001, "fold_excess"), 2), round(one(RND, 0.005, "fold_excess"), 2)]
hr = [100 * one(PRT, 0.005, "hit_rate"), 100 * one(RND, 0.01, "hit_rate"), 100 * one(PRT, 0.001, "hit_rate"), 100 * one(RND, 0.005, "hit_rate")]
add("L102", "SM topology note", "Multi-outcome excess: observed 128 vs within-unit null; vs overlap-aware nulls; proteins vs LD/MAF-matched random SNPs at similar hit rates",
    v102 + [round(h, 1) for h in hr], "128 vs %.1f; %.1f-%.1f overlap-aware; fold %.2f vs %.2f (hit %.1f%% vs %.1f%%), %.2f vs %.2f (%.1f%% vs %.1f%%)"
    % (v102[1], v102[2], v102[3], v102[4], v102[5], hr[0], hr[1], v102[6], v102[7], hr[2], hr[3]),
    TOPOS + "S1R1_topo_supp_table_c_excess.tsv", "observed_count, null_mean, fold_excess", "[S][R-topo]",
    note="The skeptic's '1.61 vs 1.46' compared different hit rates; at similar hit rates proteins sit slightly above the random mean but inside its range. The locus-minimum-p unit (z 0.15) was pre-specified to saturate: not evidence either way, do not cite it.")

# ======================================================================================================== development stage (Table 4 -> S15) audit
ta = pe.set_index(["gene", "outcome"])
for cid, g, o, agent, why in [("L110", "COL15A1", "AF", "ocriplasmin", "target is a collagen PROTEIN COMPLEX GROUP (CHEMBL2364188), not a direct single-protein target -> re-tier to 'no direct agent'"),
                              ("L111", "FN1", "CAD", "ocriplasmin", "degrades fibronectin; MR direction protective -> agent acts against the implied therapeutic direction"),
                              ("L112", "PLAU", "CKD", "urokinase", "the agent IS the PLAU protein (exogenous); MR direction risk -> adding PLAU is the harmful direction"),
                              ("L113", "DPEP1", "AF", "cilastatin", "DPEP1 inhibitor; MR direction protective -> inhibition is the harmful direction")]:
    rr = ta.loc[(g, o)]
    add(cid, "Table S15 / evidence cards", "%s-%s development-stage note (%s)" % (g, o, agent), [float(rr.mr_beta), rr.mr_direction], "%s (MR %+.3f)" % (rr.mr_direction, rr.mr_beta),
        "v3.2/_engine_root/results/phase4_annotation/table_s_pair_evidence_tiered.csv + ChEMBL (R1 skeptic1 F4)", "mr_beta, mr_direction", "[D][S]", note=why)

# ======================================================================================================== Monte Carlo / analysis settings
fw = SF["fwer"]
add("L120", "SM", "Global maximum-statistic test: exceedances / permutations / cells", [0, fw["nperm"], fw["family"], fw["p"]], "0 of 20,000; p = 5.0 × 10^-5 (resolution floor); 30 cells",
    "sidecar_claims.facts()['fwer']", "fwer", "[G]", note="Only if the sharing analysis stays in the supplement.")
add("L121", "Methods", "Analytic settings in the ARG1/DPEP1 attribution grid", LF["grid_cells"], "48", "l2l4_facts:" + LSRC["grid_cells"], "grid_cells", "[G]")
add("L122", "Methods", "LD reference panel (1000 Genomes EUR) size", LF["ld_ref_n"], "503", "l2l4_facts:" + LSRC["ld_ref_n"], "ld_ref_n", "[G]")
add("L123", "Methods", "SuSiE maximum number of effects (L)", LF["susie_L"], "10", "l2l4_facts:" + LSRC["susie_L"], "susie_L", "[G]")

# ======================================================================================================== Table 1 sizes
for o, n in SF["t1"]["n"].items():
    add("L13%s" % o, "Table 1", "Outcome total N: %s (%s)" % (o, SF["t1"]["ids"][o]), int(n), "{:,}".format(int(n)), "sidecar_claims.facts()['t1']", "n[%s]" % o, "[G]",
        note="Table 1 gains cases/controls, transformation (ln eGFR), source cohort and overlap fields (MF4, §5.5).")

# ======================================================================================================== STAGE 2a additions: parameters + Results details
import re as _re
CFG = os.path.join(V3, "paper3_deposit_request_v4/scripts/R/00_config.R")
cfg = io.open(CFG, encoding="utf-8").read()


def cval(name):
    m = _re.search(r"^%s\s*<-\s*([0-9.eE+-]+)" % name, cfg, _re.M)
    assert m, name
    return float(m.group(1))


for cid, nm, lab, disp in [("L230", "P_INSTRUMENT", "Instrument p-value threshold", "p < 10^-6"), ("L231", "P_GWAS", "Strict instrument threshold (sensitivity)", "p < 5 × 10^-8"),
                           ("L232", "CIS_WINDOW_KB", "Cis window (kb either side of the gene)", "±500 kb"), ("L233", "INFO_MIN", "Minimum imputation INFO", "0.80"),
                           ("L234", "MAF_MIN", "Minimum minor-allele frequency", "0.01"), ("L235", "F_STAT_MIN", "Minimum F-statistic", "10"),
                           ("L236", "COLOC_PP4_MIN", "Colocalization rule: minimum PP.H4", "0.80"), ("L237", "COLOC_H4H3_MIN", "Colocalization rule: minimum PP.H4/PP.H3", "5"),
                           ("L238", "SUGGESTIVE_P", "Screen threshold for colocalization", "p < 0.05"), ("L239", "FDR_ALPHA", "Benjamini-Hochberg FDR level", "5%")]:
    add(cid, "Methods: parameters", lab, cval(nm), disp, rel(CFG), nm, "[D]")
NUL = os.path.join(V3, "paper3_deposit_request_v4/scripts/scripts/_shared/coloc_null/05_run_null.R")
pr_ = _re.search(r"p1=([0-9e.-]+),p2=([0-9e.-]+),p12=([0-9e.-]+)", io.open(NUL, encoding="utf-8").read())
add("L240", "Methods: parameters", "coloc.abf priors p1, p2, p12 (default)", [float(x) for x in pr_.groups()], "p1 = p2 = 10^-4, p12 = 10^-5", rel(NUL), "coloc.abf call", "[D]")
dist_ = _re.search(r"abs\(cj\$s-gci\$s\)>([0-9e.]+)", io.open(NUL, encoding="utf-8").read())
add("L243", "Methods: decoy", "Decoy donor distance from the target gene (or another chromosome)", float(dist_.group(1)) / 1e6, "5 Mb", rel(NUL), "abs(cj$s-gci$s)>5e6", "[D]")
assert "PP.H4_p12_1e6" in pe.columns
add("L241", "Methods: parameters", "Conservative prior p12", 1e-6, "p12 = 10^-6", "tiered table column PP.H4_p12_1e6", "column name", "[D]")
add("L242", "Methods: parameters", "Bonferroni threshold over the 15,584-test family", 0.05 / fam, "3.2 × 10^-6", "0.05 / L003", "derived", "[D]")
INS = os.path.join(V3, "V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv")
q = ins.F_stat.quantile([0.25, 0.5, 0.75]).values
nsing = int((ins.groupby("gene").size() == 1).sum())
add("L200", "Results/Methods: instruments", "Instruments / proteins; proteins with a single instrument", [len(ins), ins.gene.nunique(), nsing],
    "%s instruments for %s proteins; %s single-instrument" % ("{:,}".format(len(ins)), "{:,}".format(ins.gene.nunique()), "{:,}".format(nsing)),
    rel(INS), "rows; gene", "[D]")
add("L201", "Results: instruments", "Instrument F: median (IQR); min-max", [round(q[1], 1), round(q[0]), round(q[2]), round(ins.F_stat.min()), round(ins.F_stat.max())],
    "median %.1f (IQR %s-%s); %d-%s" % (q[1], "{:,}".format(int(round(q[0]))), "{:,}".format(int(round(q[2]))), round(ins.F_stat.min()), "{:,}".format(int(round(ins.F_stat.max())))),
    rel(INS), "F_stat", "[D]")
pt = pe.set_index(["gene", "outcome"])
add("L202", "Results: map", "Established targets recovered: PCSK9-CAD and IL6R-CAD PP.H4", [float(pt.loc[("PCSK9", "CAD"), "pp_h4"]), float(pt.loc[("IL6R", "CAD"), "pp_h4"])],
    "%.3f and %.3f" % (pt.loc[("PCSK9", "CAD"), "pp_h4"], pt.loc[("IL6R", "CAD"), "pp_h4"]), "phase4 tiered table", "pp_h4", "[D]")
add("L203", "Results: map", "Pairs not study-wide significant among the 81", VF["not_bonferroni"], "40", "v4_facts:" + VSRC["not_bonferroni"], "not_bonferroni", "[G]")
add("L204", "Results: map", "RAD23B-T2D at the edge: MR p vs Bonferroni; PP.H4 at p12 1e-6; IL6R-CAD PP.H4 at p12 1e-6",
    [float(pt.loc[("RAD23B", "T2D"), "mr_pval"]), round(float(pt.loc[("RAD23B", "T2D"), "PP.H4_p12_1e6"]), 3), round(float(pt.loc[("IL6R", "CAD"), "PP.H4_p12_1e6"]), 3)],
    "p = %.2e vs 3.2e-6; %.3f; IL6R %.3f" % (pt.loc[("RAD23B", "T2D"), "mr_pval"], pt.loc[("RAD23B", "T2D"), "PP.H4_p12_1e6"], pt.loc[("IL6R", "CAD"), "PP.H4_p12_1e6"]),
    "phase4 tiered table", "mr_pval, PP.H4_p12_1e6", "[D]")
add("L205", "Results: map", "SuSiE-discordant pair", VF["susie_discordant_pairs"], "FES-AF", "v4_facts:" + VSRC["susie_discordant_pairs"], "susie_discordant_pairs", "[G]")
ta_ = pd.read_csv(os.path.join(TOPO, "S1R1_topo_supp_table_a_intercepts.tsv"), sep="\t")
ff = ta_[(ta_.source_1 == "FinnGen R10") & (ta_.source_2 == "FinnGen R10") & (ta_.outcome_1 != ta_.outcome_2)
         & ta_.outcome_1.isin(["CAD", "HF", "AF", "T2D"]) & ta_.outcome_2.isin(["CAD", "HF", "AF", "T2D"])]
add("L206", "Results: external checks / SM", "Null-z correlation between FinnGen outcomes (CAD/HF/AF/T2D; LDSC intercept)",
    [float(ff.null_z_correlation.min()), float(ff.null_z_correlation.max()), len(ff)],
    "%.2f-%.2f" % (ff.null_z_correlation.min(), ff.null_z_correlation.max()), TOPOS + "S1R1_topo_supp_table_a_intercepts.tsv", "null_z_correlation", "[S][R-topo]")
CA = os.path.join(V3, "V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv")
ca = pd.read_csv(CA).set_index(["Gene", "Endpoint"])
u_s, d_s = float(ca.loc[("UMOD", "eGFR"), "SuSiE PP.H4 (purity 0.5)"]), float(ca.loc[("DPEP1", "eGFR"), "SuSiE PP.H4 (purity 0.5)"])
u_a, d_a = float(ca.loc[("UMOD", "eGFR"), "coloc.abf PP.H4"]), float(ca.loc[("DPEP1", "eGFR"), "coloc.abf PP.H4"])
add("L207", "Results: kidney", "Multi-signal (SuSiE, purity 0.5) PP.H4 for UMOD-eGFR and DPEP1-eGFR; coloc.abf PP.H4 for each", [u_s, d_s, u_a, d_a],
    "UMOD %.2f vs %.1e; DPEP1 %.2f vs %.1e" % (u_s, u_a, d_s, d_a), rel(CA), "SuSiE PP.H4 (purity 0.5); coloc.abf PP.H4", "[D]",
    note="External 1000 Genomes LD: read qualitatively.")
rp_ = VF["renal_scrub_proteins"]
add("L208", "Results: kidney", "Renal proteins: cystatin-C concordant / with cystatin-C; BUN concordant / with BUN; the non-concordant protein",
    [rp_["cys_concordant"], rp_["with_cys"], rp_["bun_concordant"], rp_["with_bun"], rp_["cys_not_concordant"]],
    "%d of %d; %d of %d; %s" % (rp_["cys_concordant"], rp_["with_cys"], rp_["bun_concordant"], rp_["with_bun"], ", ".join(rp_["cys_not_concordant"])),
    "v4_facts:" + VSRC["renal_scrub_proteins"], "renal_scrub_proteins", "[G]")
add("L209", "Results: kidney / metabolites", "Protein-metabolite screen: pairs / proteins / metabolites; colocalizations / at proteins",
    [LF["bsr_n"], LF["bsr_prot"], LF["bsr_metab"], LF["bridges"], LF["bridge_prot"]],
    "%d pairs, %d proteins, %d metabolites; %d colocalizations at %d proteins" % (LF["bsr_n"], LF["bsr_prot"], LF["bsr_metab"], LF["bridges"], LF["bridge_prot"]),
    "l2l4_facts bsr_* / bridges", "bridge_screen_results.csv", "[G]")
add("L210", "Results: kidney", "DPEP1-cysteinylglycine attributed to DPEP1 in at most N of 48 settings", LF["dpep1_attr"]["glist"], "16 of 48",
    "l2l4_facts:" + LSRC["dpep1_attr"], "dpep1_attr glist", "[G]")
add("L211", "Results: ancestry", "Ancestry-matched non-colocalizing pairs: no East Asian outcome association (PP.H1) / distinct variants (PP.H3)",
    [LF["cen_h1"], len(LF["cen_h3_pairs"])], "%d / %d" % (LF["cen_h1"], len(LF["cen_h3_pairs"])), "l2l4_facts cen_h1 / cen_h3_pairs", "S23N", "[G]")
c2 = LF["celsr2"]
add("L212", "Results: ancestry", "CELSR2-CAD ancestry-matched: PP.H4 (1e-5 / 1e-6); East Asian vs European log-odds per SD; heterogeneity p",
    [c2["h4"], c2["h4_6"], c2["eas"], c2["eur"], c2["hetp"]], "%.3f / %.3f; %.3f vs %.3f; p = %.2f" % (c2["h4"], c2["h4_6"], c2["eas"], c2["eur"], c2["hetp"]),
    "l2l4_facts:" + LSRC["celsr2"], "celsr2", "[G]")
add("L213", "Results: ancestry", "FGF5-eGFR ancestry-matched PP.H4 (both priors)", [LF["fgf5"]["h4"], LF["fgf5"]["h4_6"]], "%.3f under both priors" % LF["fgf5"]["h4"],
    "l2l4_facts:" + LSRC["fgf5"], "fgf5", "[G]")
e16 = LF["eur_1e6"][("TGFB1", "CAD")]
add("L214", "Results: ancestry", "TGFB1-CAD ancestry-matched PP.H4 (1e-5 / 1e-6); European PP.H4 at 1e-6", [tg["h4"], tg["h4_6"], round(e16, 3)],
    "%.3f / %.3f; European %.3f" % (tg["h4"], tg["h4_6"], e16), "l2l4_facts tgfb1 / eur_1e6", "S23N", "[G]")
add("L215", "Results: ancestry", "FGFR4-T2D: PP.H4 at 1e-6; European p", [fg["h4_6"], round(fg["eurp"], 2)], "%.3f; European p = %.2f" % (fg["h4_6"], fg["eurp"]),
    "l2l4_facts:" + LSRC["fgfr4"], "fgfr4", "[G]")
add("L216", "Methods: ancestry", "Ancestry-matched analysis: eligible pairs / proteins (CKB cis-pQTL p < 5e-8)", [LF["cen_n"], LF["cen_prot"]],
    "%d pairs over %d proteins" % (LF["cen_n"], LF["cen_prot"]), "l2l4_facts cen_n / cen_prot", "S23N", "[G]")
tp_ = pd.read_csv((_P31R_EXT + "/D/data/paper31_r1_stage1/extracts/TPMI_250.2_instr.tsv"), sep="\t")
add("L217", "Methods: external checks", "TPMI PheCode 250.2 cases / controls", [int(tp_.N_case.iloc[0]), int(tp_.N_ctrl.iloc[0])],
    "{:,} / {:,}".format(int(tp_.N_case.iloc[0]), int(tp_.N_ctrl.iloc[0])), (_P31R_EXT + "/D/data/paper31_r1_stage1/extracts/TPMI_250.2_instr.tsv"), "N_case, N_ctrl", "[D]")
add("L218", "Methods: external checks", "DIAMANTE 2018 European T2D without UKB: cases / controls", [55005, 400308], "55,005 / 400,308",
    "S1R1_PRESPEC_outcome_checks.md (Paper 4 adjudication; N 455,313)", "typed", "[1]", note="Re-derive from the DIAGRAM download page before printing in Table 1.")
add("L219", "Methods: external checks", "CARDIoGRAMplusC4D 2015: cases / controls; European share", [60801, 123504, 0.77], "60,801 / 123,504; 77% European",
    (_P31R_EXT + "/D/data/EUR_new_outcome/alternative_outcomes_summary.csv + Nikpay 2015 (PMID 26343387)"), "n_cases, n_controls", "[1]", note="Re-derive before printing in Table 1.")
add("L220", "Results: ancestry", "ABO-T2D ancestry-matched at the conservative prior: PP.H4 and PP.H4/PP.H3", [LF["abo_t2d"]["h4_6"], LF["abo_t2d"]["ratio_6"]],
    "PP.H4 %.3f, ratio %.2f" % (LF["abo_t2d"]["h4_6"], LF["abo_t2d"]["ratio_6"]), "l2l4_facts:" + LSRC["abo_t2d"], "abo_t2d", "[G]")
cs_ = VF["candidate_status"]
add("L222", "Results: map", "Multi-signal candidate pairs (failed coloc.abf, SuSiE PP.H4 >= 0.80): drawn / held / dropped", [cs_["Drawn"], cs_["Held (not drawn)"], cs_["Dropped"]],
    "%d drawn, %d held, %d dropped" % (cs_["Drawn"], cs_["Held (not drawn)"], cs_["Dropped"]), "v4_facts:" + VSRC["candidate_status"], "candidate_status", "[G]")
# ---- 2b review additions (each derived; finding id in the note)
import math as _m
orw = []
for i_, r_ in oc.iloc[:6].iterrows():
    orw.append([r_.pair, _m.exp(r_.wald), _m.exp(r_.wald - 1.96 * r_.wald_se), _m.exp(r_.wald + 1.96 * r_.wald_se), float(r_.p)])
add("L223", "Results: external outcome checks", "Odds ratio per SD of protein (95% CI) and p for the six external checks", orw,
    "; ".join("%s %.2f (%.2f-%.2f)" % (a, b, c, d) for a, b, c, d, e in orw), OCS, "exp(wald +/- 1.96 SE), first-order SE", "[R]", note="F16 (2b review).")
EAX = os.path.join(V3, "V4/source_data/supp_tables/SuppTable_EAS_axis.csv")
ea = pd.read_csv(EAX)
ev = ea[ea.qc_evaluable.astype(str).str.upper() == "TRUE"]
assert len(ev) == 52
same = int((ev.direction_class == "SAME").sum())
tcol = [c for c in ev.columns if c.lower() in ("tier", "evidence_tier")][0]
sup = ev[ev[tcol] == "Supported"]; fla = ev[ev[tcol] == "Flagged"]
add("L224", "Results: ancestry", "East Asian outcome-side: evaluable pairs with the European sign; Supported vs Flagged keeping the sign (post hoc)",
    [same, len(ev), int((sup.direction_class == "SAME").sum()), len(sup), int((fla.direction_class == "SAME").sum()), len(fla)],
    "%d of %d same sign; Supported %d of %d vs Flagged %d of %d" % (same, len(ev), (sup.direction_class == "SAME").sum(), len(sup), (fla.direction_class == "SAME").sum(), len(fla)),
    rel(EAX), "qc_evaluable, direction_class, %s" % tcol, "[D]", note="F22 (2b review). Tier comparison is post hoc and exploratory.")
cak = pd.read_csv(CA)
dr_ = cak[cak.Status == "Drawn"]
kid = int(dr_.Endpoint.isin(["eGFR", "CKD"]).sum())
add("L225", "Results: kidney", "Multi-signal-only candidate pairs that are kidney pairs (of the drawn); kidney pairs among the 81",
    [kid, len(dr_), VF["renal_pairs"], VF["n_pairs"]], "%d of %d candidates; %d of %d colocalization-supported" % (kid, len(dr_), VF["renal_pairs"], VF["n_pairs"]),
    rel(CA) + " + v4_facts renal_pairs", "Status == Drawn, Endpoint in eGFR/CKD", "[D]", note="F24 (2b review). The audited pool was itself selected (26 candidates).")
supP = pe[pe.evidence_tier == "Supported"]
noarm = int((supP.replication_depth == "R0").sum())
cc_noarm = int(((supP.replication_depth == "R0") & (supP.core_complete.astype(str).str.upper() == "TRUE")).sum())
add("L226", "Results: map", "Supported pairs with no evaluable exposure-side arm; of these, meeting every internal criterion",
    [noarm, len(supP), cc_noarm], "%d of %d; %d core-complete" % (noarm, len(supP), cc_noarm), "phase4 tiered table", "replication_depth, core_complete", "[D]", note="F13 (2b review).")
tlog = io.open(os.path.join(S1, "S1R1_01_arg1_triangle/S1R1_arg1_triangle.log"), encoding="utf-8").read()
olog = io.open(os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_outcome_checks.log"), encoding="utf-8").read()
mF = _re.search(r"expected z ([-0-9.]+) \(observed z ([-0-9.]+)\)", tlog)
mD = _re.search(r"T2D z marginal A ([-0-9.]+) B ([-0-9.]+).*expected z at B under A's ratio ([-0-9.]+)", olog)
add("L227", "Results: ARG1", "Expected vs observed T2D z at signal B under signal A's ratio: FinnGen (conditional); DIAMANTE (marginal)",
    [float(mF.group(1)), float(mF.group(2)), float(mD.group(3)), float(mD.group(2))], "FinnGen %s vs %s; DIAMANTE %s vs %s" % (mF.group(1), mF.group(2), mD.group(3), mD.group(2)),
    "S1R1_arg1_triangle.log; S1R1_outcome_checks.log", "expected z lines", "[R]", note="F23 (2b review); required by the pre-stated rule.")
rm_ = rt[(rt["mode"] == "marginal") & (rt.ratio == "M per P")].iloc[0]
lr = _m.log(rm_.A_ratio / rm_.B_ratio); slr = _m.sqrt((rm_.A_se / rm_.A_ratio) ** 2 + (rm_.B_se / rm_.B_ratio) ** 2)
add("L228", "Results: ARG1", "Arginine-per-ARG1 ratio of ratios (signal A / B), 95% CI", [rm_.A_ratio / rm_.B_ratio, _m.exp(lr - 1.96 * slr), _m.exp(lr + 1.96 * slr)],
    "%.2f (%.2f-%.2f)" % (rm_.A_ratio / rm_.B_ratio, _m.exp(lr - 1.96 * slr), _m.exp(lr + 1.96 * slr)), ARGR, "marginal M per P", "[R]", note="F08 (2b review).")
rmd = rt[(rt["mode"] == "marginal") & (rt.ratio == "D per P")].iloc[0]
add("L229", "Results: ARG1", "FinnGen T2D-per-ARG1 heterogeneity on MARGINAL estimates", float(rmd.het_p), "%.3f" % rmd.het_p, ARGR, "marginal D per P", "[R]", note="F08 (2b review): verdict depends on conditioning.")
fs = LF["fgfr4_scan"]
add("L244", "Methods: ancestry", "East Asian screen that found FGFR4-T2D: tests / proteins; candidates to colocalization; at default prior; at conservative prior",
    [fs["tests"], fs["prot"], fs["cand"], fs["dflt"], 1], "%s tests over %d proteins; %d candidates; %d at default, 1 at conservative prior" % ("{:,}".format(fs["tests"]), fs["prot"], fs["cand"], fs["dflt"]),
    "l2l4_facts:" + LSRC["fgfr4_scan"], "fgfr4_scan", "[G]", note="F15 (2b review).")
RS = os.path.join(V3, "V4/source_data/supp_tables/SuppTable_renal_cystatin_BUN_scrub.csv")
rsx = pd.read_csv(RS)
colc = [c for c in rsx.columns if "concordant_and_sig" in c and "cys" in c][0]
gcol = [c for c in rsx.columns if c.lower() in ("gene", "protein")][0]
ug = rsx.drop_duplicates(gcol)
add("L245", "Results: kidney", "Kidney proteins concordant AND nominally significant on cystatin-C eGFR", [int(ug[colc].astype(str).str.upper().eq("TRUE").sum()), len(ug)],
    "%d of %d" % (ug[colc].astype(str).str.upper().eq("TRUE").sum(), len(ug)), rel(RS), colc, "[D]", note="F21 (2b review).")
DEC = os.path.join(V3, "V4/source_data/tables_from_rds/decode_replication__decode_mr_results.csv")
dc_ = pd.read_csv(DEC)
rdF = float(dc_[(dc_.gene == "RAD23B") & (dc_.outcome_name == "T2D")].decode_F_stat.iloc[0])
rdo = oc[oc.pair == "RAD23B-T2D"].iloc[0]
rdx = pd.read_csv(os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_eur_discovery_mr.tsv"), sep="\t")
rdp = bool(rdx[rdx.gene == "RAD23B"].any_proxy.iloc[0])
add("L246", "Results: map", "RAD23B-T2D near the cuts: deCODE F; external check p; instrument effect-allele frequency (palindromic G/C); discovery used an LD proxy",
    [rdF, float(rdo.p), float(rdo.eaf_exposure), rdp], "deCODE F %.2f; external p = %.3f; EAF %.2f (palindromic); proxy in discovery: %s" % (rdF, rdo.p, rdo.eaf_exposure, rdp),
    rel(DEC) + " + S1R1_outcome_checks.tsv + S1R1_eur_discovery_mr.tsv", "decode_F_stat; p; eaf_exposure; any_proxy", "[D][R]", note="F04 (2b review).")
pur_ = _re.search(r"purity ([0-9.]+)", " ".join(pd.read_csv(CA).columns))
add("L247", "Methods: parameters", "SuSiE credible-set purity filter for multi-signal candidates", float(pur_.group(1)), "purity %s" % pur_.group(1), rel(CA), "column name", "[D]")
# ---- post hoc external check of all 28 T2D/CAD pairs (spec R2_rewrite/posthoc/P31R_POSTHOC_SPEC.md, md5 573fd781)
from scipy.stats import binomtest as _bt
PH = os.path.join(V3, "R2_rewrite/posthoc/P31R_posthoc_results.tsv")
ph = pd.read_csv(PH, sep="	")
phe = ph[~ph.cls.str.startswith("not evaluable")]
phn = phe[~phe.predeclared]
PHS = "V3/R2_rewrite/posthoc/P31R_posthoc_results.tsv"
for cid, d, lab in [("L250", phe, "all T2D/CAD pairs"), ("L251", phn, "T2D/CAD pairs not pre-declared")]:
    k = int(d.same_sign.sum()); n = len(d); c = int((d.cls == "concordant").sum())
    add(cid, "Results: external outcome checks (post hoc)", "Post hoc check, %s: evaluable; same sign as discovery; concordant at p < 0.05; binomial p (descriptive)" % lab,
        [n, k, c, _bt(k, n, 0.5).pvalue], "%d of %d same sign; %d concordant; binomial p %.1e" % (k, n, c, _bt(k, n, 0.5).pvalue), PHS, "same_sign, cls", "[R]", "Y",
        note="POST HOC (user 2026-10-02). External resources overlap FinnGen slightly; not independent replication.")
tt_ = []
for t in ["Convergent", "Supported", "Flagged"]:
    d = phe[phe.tier == t]; tt_ += [t, len(d), int(d.same_sign.sum()), int((d.cls == "concordant").sum())]
add("L252", "Results: external outcome checks (post hoc)", "Post hoc check by tier: n, same sign, concordant", tt_,
    "Convergent %d/%d (%d); Supported %d/%d (%d); Flagged %d/%d (%d)" % (tt_[2], tt_[1], tt_[3], tt_[6], tt_[5], tt_[7], tt_[10], tt_[9], tt_[11]), PHS, "tier", "[R]")
op_ = phe[~phe.same_sign]
add("L253", "Results: external outcome checks (post hoc)", "Opposite-sign pair(s) in the post hoc check: gene-outcome, tier, p",
    [[a, b, c, float(d)] for a, b, c, d in op_[["gene", "outcome", "tier", "p"]].values], "; ".join("%s-%s (%s; p = %.2f)" % (a, b, c, d) for a, b, c, d in op_[["gene", "outcome", "tier", "p"]].values),
    PHS, "same_sign == False", "[R]")
add("L254", "Results: external outcome checks (post hoc)", "Post hoc: pairs whose external estimate differed from discovery (heterogeneity p < 0.05)", [int((phe.het_p < 0.05).sum()), len(phe)],
    "%d of %d" % ((phe.het_p < 0.05).sum(), len(phe)), PHS, "het_p", "[R]")
# ---- 2c review additions
add("L255", "Results: ancestry", "Ancestry-matched pair just below the colocalization rule: DNER-eGFR PP.H4", LF["dner_h4"], "%.3f" % LF["dner_h4"],
    "l2l4_facts:" + LSRC["dner_h4"], "dner_h4", "[G]", note="VN09 (2c): 32 + 7 + 1 = 40 non-colocalizing.")
msc = pe.mr_sig_class.value_counts().to_dict()
add("L256", "Results: map", "Screen significance among the 81: Bonferroni / FDR only / nominal only (the one nominal pair)",
    [msc.get("Bonferroni", 0), msc.get("FDR", 0), msc.get("Suggestive", 0), "-".join(pe[pe.mr_sig_class == "Suggestive"][["gene", "outcome"]].values[0]), msc.get("Bonferroni", 0) + msc.get("FDR", 0)],
    "%d / %d / %d (%s); %d of 81 at 5%% FDR" % (msc.get("Bonferroni", 0), msc.get("FDR", 0), msc.get("Suggestive", 0), "-".join(pe[pe.mr_sig_class == "Suggestive"][["gene", "outcome"]].values[0]),
                                               msc.get("Bonferroni", 0) + msc.get("FDR", 0)),
    "phase4 tiered table", "mr_sig_class", "[D]", note="V04 (2c).")
_p = np.sort(prim.pval.values); _k = np.arange(1, len(_p) + 1); _thr = _p[(_p <= _k / fam * 0.05).nonzero()[0].max()]
_fd = prim[prim.pval <= _thr]
_col = set(zip(pe.gene, pe.outcome))
_nc = _fd[[(g, o) not in _col for g, o in zip(_fd.gene, _fd.outcome_name)]]
add("L257", "Results/Methods: kidney", "FDR-significant associations not colocalization-supported (the candidate pool) and the kidney pairs among them",
    [len(_nc), int(_nc.outcome_name.isin(["eGFR", "CKD"]).sum())], "%d; %d kidney" % (len(_nc), _nc.outcome_name.isin(["eGFR", "CKD"]).sum()),
    rel(MR) + " + tiered table", "BH(m = 15,584) minus the 81", "[D]", note="V01 (2c): fair base for the 11-of-16 candidate comparison.")
_h = phe[phe.het_p < 0.05]
_sh = _h[_h.wald.abs() < _h.disc_b.abs()]; _lg = _h[_h.wald.abs() >= _h.disc_b.abs()]
add("L258", "Results: external outcome checks (post hoc)", "Post hoc pairs whose external estimate differed in size: shrank / larger",
    [len(_sh), len(_lg), sorted(_sh.gene), sorted(_lg.gene)], "%d shrank (%s); %d larger (%s)" % (len(_sh), ", ".join(sorted(_sh.gene)), len(_lg), ", ".join(sorted(_lg.gene))),
    PHS, "het_p < 0.05; |wald| vs |disc_b|", "[R]", note="V03 (2c).")
xpv = VF["crossplatform"]
ev_ = xpv.get("concordant", 0) + xpv.get("discordant", 0) + xpv.get("conflicting", 0)
add("L259", "Results: map", "Exposure-side corroboration: pairs with an evaluable arm / without (not assayed, weak instrument, no usable variant, assay-discordant)",
    [ev_, VF["n_pairs"] - ev_, xpv.get("not assayed (panel gap)", 0), xpv.get("not evaluable (weak instrument)", 0), xpv.get("no usable variant", 0), xpv.get("assay discordant (within platform)", 0)],
    "%d with an evaluable arm; %d without" % (ev_, VF["n_pairs"] - ev_), "v4_facts:" + VSRC["crossplatform"], "crossplatform", "[G]", note="V11 (2c).")
MC = os.path.join(V3, "v3.2/_engine_root/results/repro_phase1_hg38/pam_conditional/egfr_mirage_classification.csv")
mcx = pd.read_csv(MC)
cc_ = mcx[mcx["class"] == "concordant"]
add("L260", "Methods: kidney", "Kidney benchmark: proteins with creatinine-eGFR p < 1e-3; concordant on cystatin-C eGFR (same sign, p < 0.05)",
    [len(mcx), len(cc_), 1e-3, 0.05, float(mcx.crea_p.max()), float(cc_.cys_p.max())], "%d proteins at creatinine-eGFR p < 10^-3; %d concordant (cystatin-C same sign, p < 0.05)" % (len(mcx), len(cc_)), rel(MC), "class, crea_p, cys_p", "[D]",
    note="Thresholds confirmed in the generating script (_shared/renal_triad/01_renal_triad_mr.R: GATE1_P = 1e-3; concordant = cystatin p < 0.05 and same sign); the 89 are those testable on cystatin C. Data range agrees (max creatinine p %.4f; max cystatin p %.3f)." % (mcx.crea_p.max(), cc_.cys_p.max()))
assert len(cc_) == SF["renal"]["concordant"] and mcx.crea_p.max() < 1e-3 and cc_.cys_p.max() < 0.05
add("L221", "Results: ancestry", "Weakest CKB instrument F", round(LF["cen_Fmin"], 1), "%.1f" % LF["cen_Fmin"], "l2l4_facts:" + LSRC["cen_Fmin"], "cen_Fmin", "[G]")

# ---- 3a-ii additions (2026-10-06): corrections found by the Stage 3a figure checks (figures/review_3a/P31R_3a_STATUS.md)
from scipy.stats import norm as _norm, fisher_exact as _fx
# L261: external-check p as reported by each resource (USER 2026-10-03). For a single-variant Wald ratio with first-order SE the
# test is the variant-outcome test, so the resource's own p is the exact form of the pre-specified test: DIAMANTE prints
# beta/SE to two significant figures (beta/SE p is a lossy copy), and TPMI's p is SAIGE's score-test p.
_ord6 = ["PAM-T2D", "RAD23B-T2D", "ARG1-T2D", "PCSK9-CAD", "TGFB1-CAD", "FGFR4-T2D (EAS)"]
_oci = oc.set_index("pair")   # NB: `oa` is reassigned further up (line ~328), so index afresh
_pf = [float(_oci.loc[k, "p_out_file"]) for k in _ord6]
_nb6 = int(sum(p < 0.05 / 6 for p in _pf))
assert _nb6 == int(oc.bonferroni_6.iloc[:6].sum()) == 4, "Bonferroni count changes with the resource p"
assert all((p < 0.05) == (q < 0.05) for p, q in zip(_pf, _oci.loc[_ord6, "p"])), "classification changes with the resource p"
add("L261", "Results: external outcome checks", "External check p as reported by each resource: PAM, RAD23B, ARG1 (DIAMANTE); PCSK9, TGFB1 (CARDIoGRAMplusC4D); FGFR4 (TPMI); number passing 0.05/6",
    _pf + [_nb6], "; ".join("%s p = %s" % (k.split("-")[0], sci(p)) for k, p in zip(_ord6, _pf)) + "; %d pass 0.05/6" % _nb6, OCS, "p_out_file (first six rows)", "[D]", "Y",
    note="USER 2026-10-03: the resource's own p everywhere (decided after the checks were run). Wald p from beta/SE: PAM 2.9e-15, RAD23B 0.043, ARG1 0.014, FGFR4 0.0024; TPMI's p is SAIGE's score-test p (SPA and non-SPA agree at this variant; 3a-ii R check). Classification and the Bonferroni count are unchanged (asserted).")
# L262: thresholds of the external outcome checks (pre-specification section 5 and the post hoc spec)
_ps = io.open(os.path.join(S1, "S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), encoding="utf-8").read()
assert "The classification itself uses 0.05" in _ps and "0.0083 (0.05/6)" in _ps
add("L262", "Methods/Results: external outcome checks", "External check thresholds: classification p < 0.05; Bonferroni for six checks 0.05/6; post hoc size difference heterogeneity p < 0.05",
    [0.05, 0.05 / 6, 0.05], "classification p < 0.05; Bonferroni 0.05/6 = 0.0083; size difference p < 0.05", "V3/R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md section 5; R2_rewrite/posthoc/P31R_POSTHOC_SPEC.md",
    "pre-specified text", "[D]", note="Replaces the L238 tag (the colocalization screen threshold, a different concept with the same value) wherever the concordance rule is stated.")
# L263/L264: the multi-signal candidate pool as the pipeline built it (run log 28_candidate_pool.log), re-derived here
_RL = os.path.join(V3, "paper3_deposit_request_v4/scripts/reports/Paper3_Program/run_logs/RERUN_20260915/28_candidate_pool.log")
_rl = io.open(_RL, encoding="utf-8").read()
_m248 = int(_re.search(r"MR FDR-significant pairs = (\d+)", _rl).group(1))
_mp = _re.search(r"candidate pool = (\d+) coloc-fail pairs \| testable\(SuSiE\)=(\d+) \[(\d+) eGFR, (\d+) non-eGFR\] \| MHC=(\d+)", _rl)
_ins = pd.read_csv(INS)
_ni = _ins.groupby("gene").size()
_sg = _ins[_ins.gene.map(_ni) == 1].set_index("gene")
_pr = prim.copy()
_w = (_pr.method == "Wald ratio") & _pr.gene.isin(_sg.index)
_bx = _pr.loc[_w, "gene"].map(_sg["beta"]); _sx = _pr.loc[_w, "gene"].map(_sg["se"])
_pr["p2"] = _pr.pval
_pr.loc[_w, "p2"] = 2 * _norm.sf(np.abs(_pr.loc[_w, "b"] / np.sqrt(_pr.loc[_w, "se"] ** 2 + _pr.loc[_w, "b"] ** 2 * (_sx / _bx) ** 2)))
_p2 = np.sort(np.r_[_pr.p2.values, np.ones(fam - len(_pr))]); _k2 = np.arange(1, fam + 1)
_t2 = _p2[(_p2 <= _k2 / fam * 0.05).nonzero()[0].max()]
_f2 = _pr[_pr.p2 <= _t2]
_co = pd.read_csv(os.path.join(V3, "V4/source_data/tables_from_rds/phase3_coloc__coloc_results_final.csv"))
_cp = set(zip(_co[(_co["PP.H4"] >= 0.8) & (_co["PP.H4"] / _co["PP.H3"] >= 5)].gene, _co[(_co["PP.H4"] >= 0.8) & (_co["PP.H4"] / _co["PP.H3"] >= 5)].outcome_name))
_pool = _f2[[(g, o) not in _cp for g, o in zip(_f2.gene, _f2.outcome_name)]].copy()
_g1 = _ins.drop_duplicates("gene").set_index("gene")
_pool["mhc"] = [str(_g1.loc[x, "gene_chrom"]).replace("chr", "") == "6" and _g1.loc[x, "gene_start"] <= 34000000 and _g1.loc[x, "gene_end"] >= 25000000 for x in _pool.gene]
_tst = _pool[~_pool.mhc]
_kid = ["eGFR", "CKD"]
assert (len(_f2), len(_pool), len(_tst), int((_tst.outcome_name == "eGFR").sum()), int((_tst.outcome_name != "eGFR").sum()), int(_pool.mhc.sum())) == \
       (_m248, int(_mp.group(1)), int(_mp.group(2)), int(_mp.group(3)), int(_mp.group(4)), int(_mp.group(5))), "candidate pool does not reproduce the run log"
_cand = pd.read_csv(CA)
_drawn = _cand[_cand.Status == "Drawn"]
_dp = set(zip(_drawn.Gene, _drawn.Endpoint))
assert len(_dp) == 16 and _dp <= set(zip(_tst.gene, _tst.outcome_name)), "the 16 drawn candidates are not all in the SuSiE-testable pool"
_kp, _kt = int(_pool.outcome_name.isin(_kid).sum()), int(_tst.outcome_name.isin(_kid).sum())
add("L263", "Methods/Results: kidney", "Multi-signal candidate pool (second-order-SE FDR): FDR-significant / failed coloc.abf / in the MHC / SuSiE-testable; kidney pairs among the failures / among the testable",
    [len(_f2), len(_pool), int(_pool.mhc.sum()), len(_tst), _kp, _kt], "%d FDR; %d failed coloc.abf (%d kidney); %d MHC; %d testable (%d kidney)" % (len(_f2), len(_pool), _kp, int(_pool.mhc.sum()), len(_tst), _kt),
    rel(MR) + " + instruments + coloc table; asserted equal to " + rel(_RL), "BH(m = 15,584) on second-order Wald-ratio SEs; minus coloc.abf passes; MHC chr6:25-34 Mb", "[D]", "Y",
    note="3a figure checks (Fig 1 NV-F1-01): the pool SuSiE ran on, replacing L257's first-order 212/90 as the comparator. Reproduced by two agents and asserted against the pipeline run log.")
_dk = int(_drawn.Endpoint.isin(_kid).sum())
assert _dk == 11
_oth = (_kt - _dk, (len(_tst) - 16) - (_kt - _dk))
_pp = _fx([[_dk, 16 - _dk], list(_oth)])[1]
_um = _drawn[_drawn.Gene == "UMOD"]
assert len(_um) == 2 and _um.Endpoint.isin(_kid).all()
# locus-level sensitivity, SYMMETRIC (3a-ii verification: merging UMOD alone, on the candidate side only, was one-sided).
# Unit = a locus within an outcome group (eGFR and CKD merged as 'kidney'): genes on one chromosome whose cis windows
# start within 1 Mb of each other are chained (single linkage), in candidates and comparator alike.
_u = _tst[["gene", "outcome_name"]].copy()
_u["grp"] = np.where(_u.outcome_name.isin(_kid), "kidney", _u.outcome_name)
_u["chr"] = [str(_g1.loc[x, "gene_chrom"]).replace("chr", "") for x in _u.gene]
_u["pos"] = [float(_g1.loc[x, "gene_start"]) for x in _u.gene]
_u["cand"] = [(g, o) in _dp for g, o in zip(_u.gene, _u.outcome_name)]
_u = _u.sort_values(["grp", "chr", "pos"]).reset_index(drop=True)
_lid, _k = [], 0
for i in range(len(_u)):
    if i == 0 or _u.grp[i] != _u.grp[i - 1] or _u.chr[i] != _u.chr[i - 1] or _u.pos[i] - _u.pos[i - 1] > 1e6:
        _k += 1
    _lid.append(_k)
_u["locus"] = _lid
_lc = _u.groupby("locus").agg(kid=("grp", lambda x: x.iloc[0] == "kidney"), cand=("cand", "any"))
_ck, _co_ = int((_lc.cand & _lc.kid).sum()), int((_lc.cand & ~_lc.kid).sum())
_nk, _no = int((~_lc.cand & _lc.kid).sum()), int((~_lc.cand & ~_lc.kid).sum())
_pl = _fx([[_ck, _co_], [_nk, _no]])[1]
add("L264", "Results: kidney", "Kidney over-representation among the 16 multi-signal candidates vs the other 144 associations of the SuSiE-testable pool (pairs); the same at locus level (both arms, 1 Mb, eGFR+CKD merged)",
    [_dk, 16, _kt, len(_tst), _pp, _oth[0], _oth[0] + _oth[1], _ck, _ck + _co_, _nk, _nk + _no, _pl],
    "%d of 16 vs %d of the other %d (Fisher p = %.3f); loci: %d of %d vs %d of %d (p = %.3f)" % (_dk, _oth[0], _oth[0] + _oth[1], _pp, _ck, _ck + _co_, _nk, _nk + _no, _pl),
    "L263 + " + rel(CA) + " + instrument gene positions", "fisher_exact([[11,5],[57,87]]); loci by single-linkage <= 1 Mb within outcome group", "[D]", "Y",
    note="USER 2026-10-03: keep 'over-represented' (pair level). 3a-ii verification (R): a symmetric locus-level test does not reach 0.05, so the claim is stated for pairs, with the locus result given.")
# L265: signal-specific ARG1-T2D colocalization under the conservative prior
cA6 = co[(co.pair == "P-D") & (co.p12 == 1e-6)]
add("L265", "Results: ARG1", "Signal-specific coloc ARG1-T2D under the conservative prior: A|B / B|A PP.H4 (p12 = 1e-6)",
    [float(cA6[cA6.analysis.str.startswith("signal A")].pp_h4.iloc[0]), float(cA6[cA6.analysis.str.startswith("signal B")].pp_h4.iloc[0])],
    "%.3f / %.4f" % (cA6[cA6.analysis.str.startswith("signal A")].pp_h4.iloc[0], cA6[cA6.analysis.str.startswith("signal B")].pp_h4.iloc[0]), ARGC, "pair P-D, p12 1e-6", "[R]",
    note="3a Fig 4 (VIS-F4-03): signal A falls below the 0.80 rule under the conservative prior.")
# L266: ARG1-T2D in the ancestry-matched East Asian analysis, with and without CKB in the outcome data
_cen = pd.read_csv(os.path.join(V3, "V4_local/staging/d2-eas-census/results_census.tsv"), sep="\t")
_ckf = pd.read_csv(os.path.join(V3, "V4_local/staging/d2-eas-census/results_ckbfree.tsv"), sep="\t")
_a = _cen[(_cen.gene == "ARG1") & (_cen.outcome == "T2D")].iloc[0]; _b = _ckf[(_ckf.gene == "ARG1") & (_ckf.outcome == "T2D")].iloc[0]
assert bool(_a.ckb_overlap) and not bool(_b.ckb_overlap)
add("L266", "Results: ARG1/ancestry", "ARG1-T2D ancestry-matched: AGEN (includes CKB) PP.H3 / PP.H4; Biobank Japan (CKB-free) PP.H1 / PP.H3 / PP.H4",
    [float(_a["PP.H3"]), float(_a["PP.H4"]), float(_b["PP.H1"]), float(_b["PP.H3"]), float(_b["PP.H4"])],
    "AGEN H3 %.3f, H4 %.3f; BBJ H1 %.2f, H3 %.2f, H4 %.3f" % (_a["PP.H3"], _a["PP.H4"], _b["PP.H1"], _b["PP.H3"], _b["PP.H4"]),
    "V3/V4_local/staging/d2-eas-census/results_census.tsv + results_ckbfree.tsv", "gene ARG1, outcome T2D", "[D]",
    note="3a Fig 4/5 (NUM-F4-06, NV F5-05): the distinct-variant call holds only on CKB-containing AGEN.")
# L267: expected z at signal B under signal A's ratio: 95% interval from signal A's ratio alone, and power at the point expectation
_fe = _re.search(r"expected T2D beta at B under A's ratio: ([-0-9.]+); FinnGen SE at B ([0-9.]+)", tlog)
_zF = float(_fe.group(1)) / float(_fe.group(2))
_rF = rt[(rt["mode"] == "conditional") & (rt.ratio == "D per P")].iloc[0]
_bxB = float(rb.b_out) / float(rb.wald)
_zD = float(ra.wald) * _bxB / float(rb.se_out)
assert abs(_zF - float(mF.group(1))) < 0.006 and abs(_zD - float(mD.group(3))) < 0.006
def _iv(z, r, se):
    f = 1.96 * se / abs(r); return [z * (1 - f), z * (1 + f)]
def _pw(z):
    return float(_norm.cdf(abs(z) - 1.96) + _norm.cdf(-abs(z) - 1.96))
_iF, _iD = _iv(_zF, float(_rF.A_ratio), float(_rF.A_se)), _iv(_zD, float(ra.wald), float(ra.wald_se_delta))   # delta-method SE in both (as Figure 4c)
add("L267", "Results/Discussion: ARG1", "Expected T2D z at signal B under signal A's ratio: point, 95% interval from signal A's ratio (delta-method SE), power at the point expectation; FinnGen (conditional), DIAMANTE (marginal)",
    [_zF, _iF[0], _iF[1], _pw(_zF), _zD, _iD[0], _iD[1], _pw(_zD)],
    "FinnGen %.2f (%.2f to %.2f), power %.2f; DIAMANTE %.2f (%.2f to %.2f), power %.2f" % (_zF, _iF[1], _iF[0], _pw(_zF), _zD, _iD[1], _iD[0], _pw(_zD)),
    "S1R1_arg1_triangle.log; S1R1_arg1_ratios.tsv; S1R1_outcome_checks.tsv", "z_exp * (1 -/+ 1.96 SE_A/|ratio_A|); Phi(|z|-1.96) + Phi(-|z|-1.96)", "[R]",
    note="3a Fig 4 (NUM-F4-03). The interval carries only the uncertainty in signal A's ratio, not the observed z's own sampling SD of 1; the heterogeneity tests (L054, L055) are the tests.")

# ---- 3a-iii additions (2026-10-06): rows the figure round needed (figures/review_3a/P31R_3a_STATUS.md)
_sel = [k for k in oc.pair if not k.endswith("signal B")]
assert len(_sel) == 6 and _sel == _ord6
add("L268", "Methods/Results: external outcome checks", "Pairs selected in advance for an external outcome check (the denominator of 'six checks'; L048 holds the concordant counts)",
    len(_sel), "%d pairs" % len(_sel), OCS, "rows excluding the ARG1 signal-B locus test", "[D]", note="3a-iii (Fig 1 NUM-F1-02, Fig 2 F2-N3): 'six' was tagged to L048, the concordant count.")
_h = pd.read_csv(os.path.join(V3, "V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz"), sep="	", usecols=["rsid_ref", "pos38"])
_pA, _pB = int(_h[_h.rsid_ref == "rs2781668"].pos38.iloc[0]), int(_h[_h.rsid_ref == "rs17788484"].pos38.iloc[0])
_ce = _cen[(_cen.gene == "ARG1") & (_cen.outcome == "T2D")].iloc[0]
_pE = int(_ce.instr_pos_hg38)
assert str(_ce.instr_label) == "rs150766204" and str(_b.instr_label) == "rs150766204"
assert not (_h.pos38 == _pE).any(), "the East Asian instrument is in the European LD panel after all"
# LD of the East Asian instrument with the two signals in the 1000 Genomes East Asian panel (n = 481; plink .bed read directly)
_EP = os.path.join(CKM, "lhcmr_eas_refs/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.6")
_bim = [l.split() for l in io.open(_EP + ".bim")]
_nE = sum(1 for _ in io.open(_EP + ".fam"))
_bi = {r_[1]: i for i, r_ in enumerate(_bim)}
def _geno(rs):
    bpr = (_nE + 3) // 4
    with open(_EP + ".bed", "rb") as fh:
        fh.seek(3 + _bi[rs] * bpr); b_ = np.frombuffer(fh.read(bpr), dtype=np.uint8)
    g_ = np.zeros(bpr * 4)
    for k in range(4):
        c_ = (b_ >> (2 * k)) & 3
        g_[k::4] = np.select([c_ == 0, c_ == 2, c_ == 3, c_ == 1], [2, 1, 0, np.nan])
    return g_[:_nE]
_gI, _gA = _geno("rs150766204"), _geno("rs2781668")
_mk = ~np.isnan(_gI) & ~np.isnan(_gA)
_r2A = float(np.corrcoef(_gI[_mk], _gA[_mk])[0, 1] ** 2)
_Babs = "rs17788484" not in _bi
assert _Babs and _r2A < 0.05
add("L269", "Results/Discussion: ARG1 ancestry", "East Asian (CKB) ARG1 instrument rs150766204: GRCh38 position; distance to signal B (rs17788484); distance to signal A (rs2781668); r2 with signal A in 1000 Genomes EAS (signal B absent from that panel and rs150766204 absent from the EUR panel)",
    [_pE, abs(_pE - _pB), abs(_pE - _pA), _r2A], "rs150766204 at 6:%d; %d bp from signal B, %d bp from signal A; r2 with A in 1000G EAS %.4f; B not in the EAS panel" % (_pE, abs(_pE - _pB), abs(_pE - _pA), _r2A),
    "V3/V4_local/staging/d2-eas-census/results_census.tsv + results_ckbfree.tsv (instr); ARG1_6q23/harmonised.tsv.gz (pos38)", "instr_pos_hg38; pos38 of rs17788484 / rs2781668", "[D]",
    note="3a-iii (Fig 4 F4-NUM-01): the East Asian analyses instrument the protein beside signal B, so they do not test signal A.")
add("L270", "Figures: kidney", "Nominal threshold for the marker-triangulation fills (cystatin-C eGFR and BUN; filled = p < 0.05)", 0.05, "p < 0.05",
    "V3/R2_rewrite/figures/scripts/P31R_Fig3.R (P_NOM); 01_renal_triad_mr.R (concordant = cystatin p < 0.05)", "convention", "[D]",
    note="3a-iii (Fig 3): the BUN fill had been registered against L260, the benchmark's cystatin-C concordance threshold.")

# ---- Stage 3b extension modules: each track adds its own rows in its own file and ID range, so parallel builders never
# edit this file. ext_table1.py L300-L349, ext_table2.py L350-L399, ext_supp.py L400-L449, ext_af1.py L450-L499. Each file
# is executed in THIS namespace (it may use add(), rel(), pd, np, V3, CKM, S1 and any table read above) and must only add
# rows inside its range. A row's "independent_check" is set by the file itself (EXT_CHECKS dict) after verification.
EXT_RANGES = {"ext_table1.py": (300, 349), "ext_table2.py": (350, 399), "ext_supp.py": (400, 424), "ext_supp2.py": (425, 449), "ext_af1.py": (450, 499)}
# R2 review (Stage 0 freeze, R2_rewrite/review_R2/P31R_R2_STAGE0_FREEZE.md): one module per Stage 1 topic, created as the
# topic runs; a registered module that does not exist yet is skipped (and reported).
# Stage 5 step 1a (2026-10-16, S2-D27): the five topic modules the freeze registered at L500-L609 (ext_r2_rad23b, _ancestry, _extcoloc, _arg1, _dpep1)
# were never created (Stage 1 registered its rows through ext_r2_stage1.py at L700-L999), so they are deregistered; L500-L599 are re-allocated to
# ext_r2_stage5.py (successor and release rows). REGISTERED ONLY: the module is written in step 2 / 1c, and no row may sit in its range until it
# exists (S1R1_id_alloc.assert_unowned_ranges / assert_unowned_after). L600-L609 are free.
EXT_RANGES.update({"ext_r2_stage5.py": (500, 520), "ext_r2_renal.py": (610, 639),
                   "ext_r2_numbers.py": (640, 659), "ext_r2_release.py": (660, 689), "ext_r2_text.py": (690, 699)})
# R2 Stage 3 join (2026-10-14): one row per placeholder cell of the R2 Results / Table 3 / SM / Discussion drafts. The 211 ids exceed the
# 200 rows allocated in the Stage 0 freeze (L500-L699, split by topic: ARG1 alone has ~80 ids against 30), so the join module takes the
# first unallocated block. It reads R2_ext_ids.txt and keeps ids stable through R2_id_to_L_map.tsv (append-only). Three-digit ids only
# (the tag regex and the range check parse L + 3 digits), hence L999 is the ceiling.
# Stage 5 step 1a (S2-D27): the ceiling is L998. L999 and L971-L984 are reserved for the self-tests (S1R1_id_alloc.RESERVED_TEST_IDS): no row may take them.
EXT_RANGES.update({"ext_r2_stage1.py": (700, 998)})
# Stage 6.3b (2026-10-20, S2-D37..S2-D39): ext_r2_stage5.py is narrowed to its used span L500-L520 (the ledger holds exactly those rows); L521-L599 go to ext_r2_stage6.py
# (11 new rows + 6 successor rows). It is registered LAST on purpose: the registry runs the modules in this order, and ext_r2_stage6.py reads rows that ext_r2_release.py
# (L663, L676) and ext_r2_stage1.py (L967) build, so it must follow both.
EXT_RANGES.update({"ext_r2_stage6.py": (521, 542)})
# Stage 7 A2.2 (2026-10-21, S2-D47/48/50/54/55): ext_r2_stage6.py is narrowed to its used span L521-L542; L543-L599 go to ext_r2_stage7.py (proxy record, CCND2-T2D Steiger r2, post hoc Fisher tests,
# Wang JH 2025 overlap, CLSA per-file N). It is registered AFTER stage6 on purpose: the registry runs the modules in this order, and ext_r2_stage7.py reads rows built earlier
# (L010, L011, L252, L326 base; L668 ext_r2_release.py; L970 ext_r2_stage1.py).
EXT_RANGES.update({"ext_r2_stage7.py": (543, 599)})
sys.path.insert(0, HERE)   # python -I does not put the script's folder on sys.path
import S1R1_id_alloc as IDA  # noqa: E402
_unreg = sorted(f for f in os.listdir(HERE) if _re.fullmatch(r"ext_\w+\.py", f) and f not in EXT_RANGES)
assert not _unreg, "extension module(s) on disk but not registered in EXT_RANGES (their rows would be silently skipped): %s" % _unreg
IDA.assert_no_overlap(EXT_RANGES)
IDA.assert_unowned_ranges(EXT_RANGES, ROWS)   # ROWS holds only the base rows here: none sits in a registered range or takes a reserved test id
EXT_CHECKS = {}
for _ef, (_lo, _hi) in EXT_RANGES.items():
    _ep = os.path.join(HERE, _ef)
    if not os.path.exists(_ep):
        if _ef.startswith("ext_r2_"):
            print("extension %s: not created yet (registered L%d-L%d)" % (_ef, _lo, _hi))
        continue
    _n0 = len(ROWS)
    try:
        exec(compile(io.open(_ep, encoding="utf-8").read(), _ep, "exec"))
    except Exception as _e:   # a track's half-written file must not block the others; strict mode (final build) refuses
        del ROWS[_n0:]
        if os.environ.get("P31R_LEDGER_STRICT") == "1":
            raise
        print("!!! EXTENSION %s FAILED AND WAS SKIPPED: %s: %s" % (_ef, type(_e).__name__, _e))
        continue
    for _r in ROWS[_n0:]:
        IDA.assert_ext_row_id(_ef, _r["id"], _lo, _hi)   # L + 3 digits (fullmatch) and inside the module's range; int(id[1:4]) had let 'L5000' and 'L500x' through
    print("extension %s: %d row(s)" % (_ef, len(ROWS) - _n0))
IDA.assert_unowned_after(EXT_RANGES, ROWS, HERE)   # a registered module that does not exist has no row in its range
IDA.assert_no_reserved(ROWS)                       # L971-L984 and L999 stay free for the self-tests

# ---------------------------------------------------------------------------------------------------- independent check (R, blind; verify_R/)
EXTRA = {"L045": "This row keeps the pre-specified Wald p from beta/SE (0.0024); the text prints TPMI's own SAIGE score-test p, 1.2e-3 (L261; USER 2026-10-03). CKB pQTL beta scale (per SD?) not confirmed from local files (vR).",
         "L040": "DIAMANTE beta/SE printed to 2-3 decimals: Wald ratio good to ~2-4% (vR).",
         "L041": "DIAMANTE precision as L040 (vR).", "L042": "DIAMANTE precision as L040 (vR).",
         "L085": "CKB TGFB1 beta printed to 2 significant digits (vR)."}
VR = os.path.join(HERE, "verify_R", "S1R1_vR_compare.csv")
vr = pd.read_csv(VR).set_index("id") if os.path.isfile(VR) else None
RESOLVED = {"L054": "MATCH after reconciliation: vR's beta-space conditioning with EUR-panel genotype variances mis-scaled FinnGen at B; its z-space re-run reproduced -0.0184, p 0.033",
            "L061": "MISMATCH -> FIXED: PMID 39885775 is not a T2D-risk report (removed from the claim; cited as related)"}
for r in ROWS:
    if r["id"] in EXTRA:
        r["note"] = (r["note"] + " " if r["note"] else "") + EXTRA[r["id"]]
    if r["id"] in RESOLVED:
        r["independent_check"] = RESOLVED[r["id"]]
    elif vr is not None and r["id"] in vr.index:
        r["independent_check"] = str(vr.loc[r["id"], "status"])
    elif r["id"] in ("L250", "L251", "L252", "L253", "L254"):
        r["independent_check"] = "MATCH: re-derived in R from raw DIAMANTE / ieu-a-7 by the 2c numbers agent before reading the results (review_2c/numbers/P31R_vn_posthoc_rederived_counts.txt)"
    elif r["id"] in ("L261", "L263", "L265", "L266"):
        r["independent_check"] = "MATCH: re-derived in R from the raw files, blind to the values, by the 3a-ii verifier (wf_d6c1d649-0f5)"
    elif r["id"] in EXT_CHECKS:
        r["independent_check"] = EXT_CHECKS[r["id"]]
    elif r["id"] in ("L268", "L269"):
        r["independent_check"] = "MATCH: derived independently in R by the Figure 1/2 (L268) and Figure 4 (L269) owners in 3a-iii before the row existed (asserted in P31R_Fig1/2/4.R)"
    elif r["id"] in ("L264", "L267"):
        r["independent_check"] = "MATCH: revised after the 3a-ii verifier (symmetric loci; delta-method SE), then re-derived independently in R/Python by a SECOND verifier from the raw sources (wf_8b690056-01f: 9/14 vs 44/120, p 0.0798; both intervals)"
    elif r["id"] in ("L100", "L101", "L102"):
        r["independent_check"] = "201/201 skeptic values reproduced (70 identical, 131 within tolerance); refuting agent in Python/PyMuPDF found no blocker"
    elif r["id"] in ("L078", "L089"):
        r["independent_check"] = "raised by vR; derived here and re-derived by me in Python (2026-10-02)"
    else:
        r["independent_check"] = ""
nb = sum(1 for r in ROWS if r["load_bearing"] == "Y" and not r["independent_check"] and r["value"] != '"PENDING"')
print("load-bearing rows without an independent check (excluding PENDING):", nb, [r["id"] for r in ROWS if r["load_bearing"] == "Y" and not r["independent_check"] and r["value"] != '"PENDING"'])

# ---------------------------------------------------------------------------------------------------- write
out = os.path.join(HERE, "S1R1_CLAIMS_LEDGER.csv")
_tmp = out + ".%d.tmp" % os.getpid()   # atomic: parallel builders rebuild the ledger; a reader never sees a half file
with io.open(_tmp, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ROWS[0].keys())); w.writeheader(); w.writerows(ROWS)
os.replace(_tmp, out)
ids = [r["id"] for r in ROWS]
assert len(ids) == len(set(ids)), "duplicate ids"
md = ["# Paper 3.1 rewrite -- claims ledger (generated by S1R1_build_ledger.py; do not edit by hand)", "",
      "Independent check: verify_R/ (R, blind to these values; 44 in-scope rows + 3 extra). Column 'check'.", "",
      "Tags: [G] gated build fact module; [D] read here from a V4 table; [R] Stage 1 own code; [S] reproduced by an R1 independent skeptic; [1] one check only (re-derive before print).", "",
      "| id | section | claim | display | tag | LB | check | source / key | note |", "|---|---|---|---|---|---|---|---|---|"]
for r in ROWS:
    md.append("| %s | %s | %s | %s | %s | %s | %s | `%s` · %s | %s |" % (r["id"], r["section"], r["claim"], r["display"] or r["value"], r["tag"], r["load_bearing"],
                                                                     r["independent_check"], r["source"], r["key"], r["note"]))
io.open(os.path.join(HERE, "S1R1_CLAIMS_LEDGER.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")
with io.open(os.path.join(HERE, "S1R1_CLAIMS_LEDGER_BLIND.csv"), "w", encoding="utf-8", newline="") as fh:   # for the independent check: no values
    keep = ["id", "section", "claim", "source", "key", "tag", "load_bearing"]
    w = csv.DictWriter(fh, fieldnames=keep); w.writeheader(); w.writerows([{k: r[k] for k in keep} for r in ROWS])
print("ledger rows: %d (load-bearing %d) -> %s" % (len(ROWS), sum(1 for r in ROWS if r["load_bearing"]), out))
