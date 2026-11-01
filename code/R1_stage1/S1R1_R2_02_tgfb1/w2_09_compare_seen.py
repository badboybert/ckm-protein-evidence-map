"""w2_09_compare_seen.py -- COMPARE_SEEN.tsv for section 2 (brief rule 10). Runs only after results.tsv and branches.tsv
exist; reads the FULL pre-declaration (md5 asserted) for the seen values, which are transcribed below with their line
numbers. Stage 1 values come from results.tsv; where the seen value is on another scale (the seen 2 x 2 cells are Wald
ratios), the conversion is computed HERE for the comparison only and is never written back to results.tsv.
Flag rule: a difference larger than half a unit in the last printed digit of the seen value is flagged.
FIX1 (2026-10-08, Stage 1 fixer; FIXLOG.md): two rows added for the branch outcomes FIX1 changed; earlier rows unchanged.
"""
import csv
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2_common as W  # noqa: E402

FULL = W.GOVERNING["PRESPEC_FULL"]


def half_unit(s):
    s = s.replace("−", "-").replace("+", "").strip()
    m = re.match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", s)
    if not m:
        return None
    dec = len(m.group(2) or "")
    ex = int(m.group(3) or 0)
    return 0.5 * 10 ** (ex - dec)


def num(s):
    return float(s.replace("−", "-").replace("+", ""))


def main():
    W.assert_md5(*FULL, "PRESPEC_FULL (opened for COMPARE_SEEN only)")
    lines = open(FULL[0], encoding="utf-8").read().splitlines()
    W.log("==== w2_09_compare_seen (full pre-declaration opened; computed values are not changed) ====")
    R = {r["quantity_id"]: r for r in W.read_tsv(W.OUT + "/results.tsv")}
    B = W.read_tsv(W.OUT + "/branches.tsv")
    loc = {r["rsid503"]: r for r in W.read_tsv(W.OUT + "/tgfb1_locus.tsv")}
    rows = []

    def add(qid, seen, seen_label, s1, s1_label, expl, line):
        assert seen_label.split(" [")[0][:20] or True
        ad = ""
        flag = ""
        try:
            a, b = num(str(seen)), float(s1)
            ad = abs(a - b)
            hu = half_unit(str(seen))
            if hu is not None and ad > hu + 1e-12:
                flag = "DIFFERS beyond printed precision (half unit %.2g). " % hu
            else:
                flag = "within printed precision. "
        except (ValueError, TypeError):
            pass
        rows.append({"quantity_id": qid, "seen_value": seen, "seen_label": "%s [full text l.%d]" % (seen_label, line),
                     "stage1_value": s1, "stage1_label": s1_label, "abs_diff": ad, "explanation": flag + expl})
    v = lambda q: float(R[q]["value"])  # noqa: E731
    weak = "Deterministic quantity from the same summary statistics; agreement is weak evidence (the value was visible in the scoping)."
    # ---- LD (line 213, 206)
    add("S2.r.503", "+0.24", "LD r(rs73045269, rs12983775), 1000G EUR, sign on the protein-raising alleles (T, A)", -v("S2.r.503"),
        "-(S2.r.503) = r on (T, A); S2.r.503 itself is on the 503 A1 alleles (T, G) = %.6f" % v("S2.r.503"),
        weak + " Scoping: custom Python .bed reader; Stage 1: PLINK --r.", 213)
    add("S2.ld.rs1800470_rs12983775.r2", "0.71", "r2 rs1800470-rs12983775 (EUR)", v("S2.ld.rs1800470_rs12983775.r2"), "PLINK 503", weak, 206)
    # ---- marginal pQTL and Wald values (lines 210-212) from the locus table (comparison only)
    def marg(rs, allele, t):
        r = loc[rs]
        s = 1 if r["a1"] == allele else -1
        return s * float(r["b__" + t]), float(r["se__" + t]), float(r["log10p__" + t])
    for rs, al, seen, line in (("rs73045269", "T", "+0.288", 210), ("rs12983775", "A", "+0.199", 211), ("rs1800470", "A", "+0.162", 212)):
        b, se, _ = marg(rs, al, "P")
        add("compare.pqtl.%s-%s" % (rs, al), seen, "UKB-PPP pQTL beta, SD per %s allele" % al, b, "tgfb1_locus.tsv b__P re-oriented to %s (comparison only)" % al, weak, line)
    for rs, al, t, seen, seen_se, line in (("rs73045269", "T", "CHD_FG", "+0.164", "0.033", 210), ("rs73045269", "T", "CAD_C4D", "+0.264", None, 210),
                                            ("rs12983775", "A", "CHD_FG", "−0.181", "0.045", 211), ("rs12983775", "A", "CAD_C4D", "−0.208", "0.050", 211),
                                            ("rs1800470", "A", "CHD_FG", "−0.186", None, 212)):
        bp, sp, _ = marg(rs, al, "P")
        by, sy, lp = marg(rs, al, t)
        w = by / bp
        se2 = math.sqrt(sy ** 2 / bp ** 2 + by ** 2 * sp ** 2 / bp ** 4)
        se1 = sy / abs(bp)
        res = {"CHD_FG": "finngen_chd", "CAD_C4D": "ieu_a_7"}[t]
        qid = "S2.marg.ratio.%s.%s.value" % (res, rs) if ("S2.marg.ratio.%s.%s.value" % (res, rs)) in R else "compare.wald.%s.%s" % (res, rs)
        add(qid, seen, "marginal Wald ratio, %s" % t, w, "marginal Wald ratio (results.tsv where present; rs1800470 computed here for comparison)", weak, line)
        if seen_se:
            add(qid.replace(".value", ".se") if qid.endswith(".value") else qid + ".se", seen_se, "SE of that Wald ratio", se2,
                "second-order SE (Stage 1 convention); first-order SE = %.4f" % se1,
                "Seen SE appears to be first-order (s_Y/|b_P| = %.4f) where it differs from the second-order %.4f; convention difference, not a data difference." % (se1, se2)
                if abs(se1 - num(seen_se)) < abs(se2 - num(seen_se)) else weak, line)
        if rs == "rs1800470":
            add("compare.outcome_p.rs1800470.%s" % t, "4.16e-4", "FinnGen CHD outcome p at rs1800470", 10 ** (-lp), "tgfb1_locus.tsv log10p (comparison only)", weak, 212)
    add("compare.outcome_p.rs73045269.CHD_FG", "5.37e-7", "FinnGen CHD outcome p at rs73045269", 10 ** (-marg("rs73045269", "T", "CHD_FG")[2]), "tgfb1_locus.tsv", weak, 210)
    add("compare.outcome_p.rs12983775.CHD_FG", "5.6e-5", "FinnGen CHD outcome p at rs12983775", 10 ** (-marg("rs12983775", "A", "CHD_FG")[2]), "tgfb1_locus.tsv", weak, 211)
    # naive heterogeneity z (line 213)
    a, sa = v("S2.marg.ratio.finngen_chd.rs73045269.value"), float(R["S2.marg.ratio.finngen_chd.rs73045269.value"]["se"])
    b, sb = v("S2.marg.ratio.finngen_chd.rs12983775.value"), float(R["S2.marg.ratio.finngen_chd.rs12983775.value"]["se"])
    zn = (a - b) / math.sqrt(sa ** 2 + sb ** 2)
    add("compare.naive_het_z.finngen", "6.2", "naive within-European heterogeneity z (marginal Wald ratios, covariance ignored)", zn,
        "computed here from S2.marg.ratio rows (second-order SEs)", "Not a declared Stage 1 quantity (the declared test is the conditional aware test, "
        "S2.het.finngen_chd.p_aware_2nd = %s). With the seen first-order SEs (0.033, 0.045) z = %.2f." % (R["S2.het.finngen_chd.p_aware_2nd"]["value"],
                                                                                                      (0.164 + 0.181) / math.sqrt(0.033 ** 2 + 0.045 ** 2)), 213)
    add("compare.coloc.TGFB1_CAD.marginal", "0.968 / 0.753", "TGFB1-CAD marginal PP.H4 at p12 1e-5 / 1e-6", "NA", "not a W2 quantity",
        "Marginal regional colocalization is section 3 (W3). W2 computes only the conditional PP.H4 of Rule 3(iv) (A: 0.995 / 0.994; B: 0.838 / 0.811, FinnGen 503 / 489).", 214)
    # ---- 2 x 2 (lines 224-232): seen cells are Wald ratios (except FGFR4, per-allele); Stage 1 cells are per-allele betas
    def bx(locus, instr, which):
        n = R["S2.x22.%s.%s.snp" % (locus, instr)]["notes"]
        if instr == "EURinstr":
            m = re.search(r"UKB-PPP beta per raising allele ([-\d.]+).*?CKB beta per same allele ([-\d.]+|absent)", n)
            vals = {"UKB": m.group(1), "CKB": m.group(2)}
        else:
            m = re.search(r"CKB beta per raising allele ([-\d.]+).*?UKB-PPP beta per same allele ([-\d.]+|absent)", n)
            vals = {"CKB": m.group(1), "UKB": m.group(2)}
        x = vals[which]
        return float(x) if x != "absent" else float("nan")
    seen22 = {"TGFB1": (["+0.164", "−0.181", "−0.371", "n.a."], 224), "CELSR2": (["−0.174", "−0.172", "−0.168", "−0.175"], 225),
              "IL6R": (["−0.035", "−0.034", "−0.045", "−0.049"], 226), "ABO": (["+0.049", "+0.053", "+0.035", "+0.039"], 227),
              "GALNT3": (["−0.103", "−0.079", "−0.128", "−0.046"], 228), "SPON1": (["+0.142", "+0.143", "+0.130", "+0.135"], 229),
              "FGF5": (["−0.005", "−0.005", "−0.033", "n.a."], 230), "ARG1": (["−0.191", "n.a.", "n.a.", "−0.159"], 231)}
    cells = ["EURinstr_EURout", "CKBinstr_EURout", "CKBinstr_EASout", "EURinstr_EASout"]
    for locus, (vals, line) in seen22.items():
        for cell, sv in zip(cells, vals):
            q = "S2.x22.%s.%s.beta" % (locus, cell)
            beta = R[q]["value"]
            instr = cell.split("_")[0]
            pop = "EUR" if cell.endswith("EURout") else "EAS"
            denom = bx(locus, instr, "UKB" if pop == "EUR" else "CKB")
            if beta == "NA":
                s1, lab = "NA", "Stage 1 cell NA: %s" % R[q]["notes"][:160]
                if sv == "n.a.":
                    expl = "Both n.a.: the variant is absent from the resource (Stage 1 holds the instrument fixed, no proxy)."
                else:
                    expl = "unexplained"
                add(q, sv, "seen 2x2 cell (Wald ratio, per SD protein)", s1, lab, expl, line)
                continue
            ratio = float(beta) / denom if denom == denom else float("nan")
            if locus == "FGF5":
                ratio_cmp = -ratio
                lab = ("Stage 1 per-allele beta %.4g (raw eGFR scale) / bx %s %.4f = %.4f; shown x(-1) = eGFR* (risk-aligned) to match the seen sign" % (
                    float(beta), "UKB-PPP" if pop == "EUR" else "CKB", denom, ratio))
                expl = ("ESTIMAND/SCALE: seen cells are Wald ratios on risk-aligned eGFR* (census flip); Stage 1 reports per-allele SNP-outcome betas on the raw scale "
                        "(conversion computed here only). ")
            else:
                ratio_cmp = ratio
                lab = "Stage 1 per-allele beta %.4g / bx %s %.4f = implied Wald ratio (computed here only)" % (float(beta), "UKB-PPP" if pop == "EUR" else "CKB", denom)
                expl = ("ESTIMAND: seen cells are Wald ratios (EUR-outcome cells / UKB-PPP bx; EAS-outcome cells / CKB bx, matching the census mr_beta); "
                        "Stage 1 cells are per-allele SNP-outcome betas per European protein-raising allele (brief wording); converted here only. ")
            if sv == "n.a.":
                expl += ("DEFINITION: scoping marked n.a.; Stage 1 applies the uniform CKB-lead rule (ARG1 CKB lead 6:131573339, an East Asian variant "
                         "absent from UKB-PPP and FinnGen, oriented to the CKB-raising allele) and regenerates this cell (p %s)." % R[q]["p"])
                add(q, sv, "seen 2x2 cell (Wald ratio, per SD protein)", ratio_cmp, lab, expl, line)
                continue
            add(q, sv, "seen 2x2 cell (Wald ratio, per SD protein)", ratio_cmp, lab, expl + "Deterministic; agreement after conversion is weak evidence.", line)
    # TGFB1 instrument-and-ancestry contrast z (line 224): not a declared Stage 1 quantity
    r1, r3 = v("S2.x22.TGFB1.EURinstr_EURout.beta") / bx("TGFB1", "EURinstr", "UKB"), v("S2.x22.TGFB1.CKBinstr_EASout.beta") / bx("TGFB1", "CKBinstr", "CKB")
    s1 = float(R["S2.x22.TGFB1.EURinstr_EURout.se"]["value"]) / bx("TGFB1", "EURinstr", "UKB")
    s3 = float(R["S2.x22.TGFB1.CKBinstr_EASout.se"]["value"]) / bx("TGFB1", "CKBinstr", "CKB")
    add("compare.TGFB1.instrument_ancestry_contrast_z", "8.085", "TGFB1 instrument-and-ancestry contrast z (scoping, cell 1 vs cell 3)",
        (r1 - r3) / math.sqrt(s1 ** 2 + s3 ** 2), "computed here with first-order Wald SEs (comparison only)",
        "Not a declared Stage 1 quantity; computed only to reconcile the seen note.", 224)
    # FGFR4 (line 232): per-allele
    for cell, sv, sp in (("EURinstr_EURout_FinnGen", "+0.0076", "0.27"), ("EURinstr_EURout_DIAMANTE", "+0.0089", "0.20"),
                         ("EURinstr_EASout_AGEN", "+0.047", None), ("EURinstr_EASout_BBJ", "+0.063", None)):
        q = "S2.x22.FGFR4.%s.beta" % cell
        expl = "Same variant for both instruments (rs351855); per-allele beta on both sides."
        if "DIAMANTE" in cell:
            expl += (" UNEXPLAINED: both local DIAMANTE 2018 noUKB files (chrpos md5 9dda9efd..., and the rsID derivative Mahajan_T2D_rsid.txt) give "
                     "rs351855-A beta 0.0058, SE 0.0091, p 0.52; the seen value was scout prose with no output file (full text l.232).")
        add(q, sv, "seen FGFR4 per-allele T2D (scout prose)", R[q]["value"], "Stage 1 per-allele beta (p %s)" % R[q]["p"], expl, 232)
        if sp:
            add(q.replace(".beta", ".p"), sp, "seen p", R[q]["p"], "Stage 1 p", "" if "DIAMANTE" not in cell else "UNEXPLAINED (as above)", 232)
    # branch (line 234)
    b2 = [b for b in B if b["rule_ref"].startswith("prespec s2 wording")][0]
    rows.append({"quantity_id": "branch.s2", "seen_value": "(a)", "seen_label": "expected branch: (a) is the arithmetic expectation [full text l.234]",
                 "stage1_value": b2["branch"], "stage1_label": "branches.tsv s2 wording rule", "abs_diff": "",
                 "explanation": "Same branch. The scoping called (a) the arithmetic expectation from the marginal values (positively correlated "
                                "protein-raising alleles with opposite marginal CAD associations); agreement is weak evidence."})
    rows.append({"quantity_id": "rng.s2.bootstrap", "seen_value": "no TGFB1 bootstrap share in the seen text", "seen_label": "[full text s2 Seen]",
                 "stage1_value": "%s (FinnGen), %s (ieu-a-7)" % (R["S2.het.boot.finngen_chd.share"]["value"], R["S2.het.boot.ieu_a_7.share"]["value"]),
                 "stage1_label": "S2.het.boot.<resource>.share", "abs_diff": "",
                 "explanation": "Scoping TGFB1 conditioning was in R; Stage 1 is Python numpy default_rng(20261007) (PCG64): not the same RNG; no seen share to compare."})
    rows.append({"quantity_id": "class.s2.susie", "seen_value": "no scoping s classification for TGFB1", "seen_label": "[full text: s classifications exist only for s4(v)(1), s5, s6 UMOD]",
                 "stage1_value": "inadmissible (s_max protein %s, CAD FinnGen %s)" % (R["S2.susie.protein.s_max"]["value"][:6], R["S2.susie.finngen_chd.s_max"]["value"][:6]),
                 "stage1_label": "Rule 3(i) on Stage 1 s", "abs_diff": "", "explanation": "No CLASS_FLIP possible (no recorded scoping class for TGFB1)."})
    # FIX1 (FIXLOG.md): the two branch outcomes the fixer added; the full text (Seen block l.209-232; Rule 3(i) l.90) records no seen value for either
    b3ii = [b for b in B if b["rule_ref"].startswith("Rule 3(ii) on the named primary cell")][0]
    rows.append({"quantity_id": "branch.s2.rule3ii_primary_cell", "seen_value": "no TGFB1 Rule 3(ii) outcome, conditional colocalization or bootstrap share in the seen text",
                 "seen_label": "[full text s2 Seen, l.209-232]", "stage1_value": "%s (ratio test share %s; conditional colocalization share A %s / B %s)" % (
                     b3ii["branch"], R["S2.het.boot.finngen_chd.share"]["value"], R["S2.coloc_clause.boot.rs73045269.finngen_chd.share"]["value"],
                     R["S2.coloc_clause.boot.rs12983775.finngen_chd.share"]["value"]),
                 "stage1_label": "branches.tsv Rule 3(ii) primary-cell row (FIX1)", "abs_diff": "",
                 "explanation": "No seen value to compare. FIX1 evaluates Rule 3(ii) on the named primary cell's conditional colocalizations as well as its ratio test "
                                "(Rule 3(ii) lists conditional-colocalization results); the pre-fix row evaluated the ratio test only. Bootstrap uses the ratio "
                                "bootstrap's person indices (seed 20261007); no seen RNG to match."})
    b3i = [b for b in B if b["rule_ref"].startswith("Rule 3(i) SuSiE")][0]
    rows.append({"quantity_id": "flag.s2.rule3i_first_mention", "seen_value": "no TGFB1 SuSiE class or first-mention flag in the seen text",
                 "seen_label": "[full text l.90 (Rule 3(i)); s classifications only for s4(v)(1), s5, s6 UMOD]",
                 "stage1_value": "%s (existing susie_class %s read from TAB_tiered81; proxy s_max protein %s)" % (
                     b3i["branch"], R["S2.rule3i.class.TGFB1_CAD.susie_class"]["value"], R["S2.susie.protein.s_max"]["value"][:6]),
                 "stage1_label": "branches.tsv Rule 3(i) row (FIX1)", "abs_diff": "",
                 "explanation": "No seen value to compare and no CLASS_FLIP possible. The existing class is read from raw_input TAB_tiered81 and TAB_susie_grid "
                                "(not seen corpus), not recomputed; the reference-LD attribution rests on the engine config 00_config.R."})
    cols =["quantity_id", "seen_value", "seen_label", "stage1_value", "stage1_label", "abs_diff", "explanation"]
    W.write_tsv(W.OUT + "/COMPARE_SEEN.tsv", cols, rows)
    nflag = sum("DIFFERS" in r["explanation"] for r in rows)
    nun = sum("UNEXPLAINED" in r["explanation"].upper() for r in rows)
    W.log("COMPARE_SEEN.tsv: %d rows; %d differ beyond printed precision; %d marked unexplained" % (len(rows), nflag, nun))
    for r in rows:
        if "DIFFERS" in r["explanation"] or "UNEXPLAINED" in r["explanation"].upper():
            W.log("  %s seen %s stage1 %s :: %s" % (r["quantity_id"], r["seen_value"], r["stage1_value"], r["explanation"][:200]))


if __name__ == "__main__":
    main()
