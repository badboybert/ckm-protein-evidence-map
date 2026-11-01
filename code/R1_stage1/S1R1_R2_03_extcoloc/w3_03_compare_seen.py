#!/usr/bin/env python
"""W3 step 3 (brief §0 rule 10, compute first, compare second): writes COMPARE_SEEN.tsv.

Run only AFTER results.tsv and branches.tsv exist. It reads results.tsv (values are copied, never recomputed),
the full pre-declaration (md5 3319dc39...; §3 Seen block and the §1/§2 Seen lines that touch these four pairs),
and TAB_tiered81 (raw_input) to name the analysis behind two seen colocalization values. It asserts that
results.tsv and branches.tsv still have the md5s recorded in run.log before the full text was opened.
Seen values below are transcribed from the full pre-declaration text (quoted line by line in seen_label).
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
import csv, hashlib, json, os, re, sys, time

OUT = (_P31R_V3 + "/R1_stage1/S1R1_R2_03_extcoloc")
FULL = (_P31R_V3 + "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md")
FULL_MD5 = "3319dc39740d9f5b3b8822133d9108b7"
AMEND = (_P31R_V3 + "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md")
AMEND_MD5 = "377508832ea9bd384133e869c456c909"
RESJ = (_P31R_V3 + "/R1_stage1/S1R1_R2_00_resources/RESOURCES.json")


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S") + " [compare_seen] " + msg
    print(line, flush=True)
    with open(os.path.join(OUT, "run.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def main():
    runlog = open(os.path.join(OUT, "run.log"), encoding="utf-8").read()
    m = re.search(r"results.tsv md5 before opening the full pre-declaration: ([0-9a-f]{32}); branches ([0-9a-f]{32})", runlog)
    assert m, "run.log lacks the pre-compare md5 line"
    assert md5(os.path.join(OUT, "results.tsv")) == m.group(1), "results.tsv changed after the full text was opened"
    assert md5(os.path.join(OUT, "branches.tsv")) == m.group(2), "branches.tsv changed after the full text was opened"
    log("results.tsv %s and branches.tsv %s unchanged since before the full text was opened" % (m.group(1), m.group(2)))
    assert md5(FULL) == FULL_MD5 and md5(AMEND) == AMEND_MD5
    log("full pre-declaration %s and full amendment %s opened (md5 OK) AFTER results.tsv/branches.tsv were written" % (FULL_MD5, AMEND_MD5))
    full = open(FULL, encoding="utf-8").read()
    s3 = full[full.index("## 3. External regional colocalization"):full.index("## 4. ARG1")]
    for needle in ["- Nothing for these four regional tests.", "PCSK9–CAD PP.H3 0.99998", "PAM–T2D MR −0.112 (p 0.0026, same sign), PP.H4 0.26 (H1 most probable)"]:
        assert needle in s3, needle
    for needle in ["**Colocalization:** TGFB1–CAD PP.H4 0.968 at the default prior and 0.753 at p12 10⁻⁶.",
                   "conservative PP.H4 0.8273", "resource p 0.038 at the lead (the declared check's instrument) and 0.042 at the proxy",
                   "**UKB-PPP at rs7856112:** T β −0.0637"]:
        assert needle in full, needle
    amend = open(AMEND, encoding="utf-8").read()
    assert not re.search(r"PAM|RAD23B|PCSK9|TGFB1|§3", amend), "amendment mentions a §3 pair: review it"
    res = {r["quantity_id"]: r for r in csv.DictReader(open(os.path.join(OUT, "results.tsv"), encoding="utf-8"), delimiter="\t", quoting=csv.QUOTE_NONE)}
    br = list(csv.DictReader(open(os.path.join(OUT, "branches.tsv"), encoding="utf-8"), delimiter="\t", quoting=csv.QUOTE_NONE))

    # TAB_tiered81 (raw_input) names the analysis behind the §1/§2 seen colocalization values
    R = {e["id"]: e for e in json.load(open(RESJ, encoding="utf-8"))}
    t = R["TAB_tiered81"]
    assert t["role"] == "raw_input" and md5(t["path"]) == t["md5"]
    tier = {(r["gene"], r["outcome"]): r for r in csv.DictReader(open(t["path"], encoding="utf-8"))}
    log("TAB_tiered81 %s read (raw_input; md5 OK) only to name the analysis behind the seen §1/§2 colocalization values" % t["md5"])

    rows = []
    def add(q, sv, sl, s1v, s1l, diff, expl):
        rows.append(dict(quantity_id=q, seen_value=sv, seen_label=sl, stage1_value=s1v, stage1_label=s1l, abs_diff=diff, explanation=expl))

    NOT_SEEN = "full pre-declaration §3 Seen: 'Nothing for these four regional tests.'"
    for pair in ("PAM_T2D", "RAD23B_T2D", "TGFB1_CAD", "PCSK9_CAD"):
        for p12 in ("1e-05", "1e-06"):
            for h in ("H4", "H3"):
                q = "S3.%s.pmaf.%s.%s" % (pair, p12, h)
                add(q, "none", NOT_SEEN, res[q]["value"], "Stage 1 primary input (p + MAF + N)", "NA",
                    "No seen value exists for the §3 regional tests; the Stage 1 value is new, not a reproduction.")
        b = [x for x in br if x["rule_ref"].startswith("§3 wording rule (branches 1-4)") and x["inputs"].startswith(pair)]
        assert len(b) == 1
        add("S3.%s.branch" % pair, "none", NOT_SEEN, b[0]["branch"], "Stage 1 §3 branch (primary input)", "NA",
            "Branch first determined by Stage 1.")
    ctl = [x for x in br if x["rule_ref"].startswith("§3 PCSK9 control criterion")][0]
    add("S3.PCSK9_CAD.control", "none", NOT_SEEN, ctl["branch"], "Stage 1 PCSK9 control criterion", "NA",
        "First determined by Stage 1. The criterion's consequences were written after the East Asian PP.H3 (next row) was seen.")

    # Seen values from other analyses that touch these pairs (no agreement expected)
    add("EAS.PCSK9_CAD.PP.H3 (native East Asian arm)", "0.99998", "§3 Seen, native East Asian arm (a different analysis)",
        "%s (1e-5) / %s (1e-6)" % (res["S3.PCSK9_CAD.pmaf.1e-05.H3"]["value"], res["S3.PCSK9_CAD.pmaf.1e-06.H3"]["value"]),
        "S3.PCSK9_CAD.pmaf.{1e-05,1e-06}.H3 (European external, CARDIoGRAMplusC4D)", "NA",
        "Different analysis (CKB East Asian cis-pQTL x East Asian CAD vs UKB-PPP European x CARDIoGRAMplusC4D); printed for context only, not a comparison. Both favour distinct variants more than the discovery screen did (next rows).")
    add("EAS.PAM_T2D.PP.H4 (native East Asian arm)", "0.26 (H1 most probable)", "§3 Seen, native East Asian arm (a different analysis)",
        "%s (1e-5)" % res["S3.PAM_T2D.pmaf.1e-05.H4"]["value"], "S3.PAM_T2D.pmaf.1e-05.H4 (European external, DIAMANTE noUKB)", "NA",
        "Different analysis and ancestry; context only.")
    add("EAS.PAM_T2D.MR (native East Asian arm)", "-0.112 (p 0.0026, same sign)", "§3 Seen, native East Asian arm", "NOT RUN by W3",
        "W3 runs no MR", "NA", "Not a §3 quantity.")
    tg = tier[("TGFB1", "CAD")]
    add("DISC.TGFB1_CAD.PP.H4 (discovery screen)", "0.968 (1e-5) / 0.753 (1e-6)", "§2 Seen: 'TGFB1–CAD PP.H4 0.968 at the default prior and 0.753 at p12 10⁻⁶'",
        "%s / %s" % (res["S3.TGFB1_CAD.pmaf.1e-05.H4"]["value"], res["S3.TGFB1_CAD.pmaf.1e-06.H4"]["value"]),
        "S3.TGFB1_CAD.pmaf.{1e-05,1e-06}.H4 (CARDIoGRAMplusC4D)", "%.4f / %.4f" % (abs(float(res["S3.TGFB1_CAD.pmaf.1e-05.H4"]["value"]) - 0.968), abs(float(res["S3.TGFB1_CAD.pmaf.1e-06.H4"]["value"]) - 0.753)),
        "Different estimand: the seen values equal the discovery-screen colocalization in TAB_tiered81 (pp_h4 %s, PP.H4_p12_1e6 %s), i.e. the screen's own outcome GWAS, not CARDIoGRAMplusC4D. No agreement expected; the external test is new." % (tg["pp_h4"], tg["PP.H4_p12_1e6"]))
    rd = tier[("RAD23B", "T2D")]
    add("DISC.RAD23B_T2D.PP.H4_p12_1e6 (discovery screen)", "0.8273", "§1 Seen: 'RAD23B near-threshold criteria: ... conservative PP.H4 0.8273'",
        res["S3.RAD23B_T2D.pmaf.1e-06.H4"]["value"], "S3.RAD23B_T2D.pmaf.1e-06.H4 (DIAMANTE noUKB)", "%.4f" % abs(float(res["S3.RAD23B_T2D.pmaf.1e-06.H4"]["value"]) - 0.8273),
        "Different estimand: equals TAB_tiered81 RAD23B–T2D PP.H4_p12_1e6 %s (discovery screen; pp_h4 %s at 1e-5). The external DIAMANTE regional test (branch 4) is new; the contrast with the screen is a Stage 1 result, not a discrepancy." % (rd["PP.H4_p12_1e6"], rd["pp_h4"]))
    for gene, oc in (("PAM", "T2D"), ("PCSK9", "CAD")):
        r = tier[(gene, oc)]
        q5, q6 = "S3.%s_%s.pmaf.1e-05.H4" % (gene, oc), "S3.%s_%s.pmaf.1e-06.H4" % (gene, oc)
        add("DISC.%s_%s.PP.H4 (discovery screen, TAB_tiered81; not transcribed in the pre-declaration)" % (gene, oc), "%s / %s" % (r["pp_h4"], r["PP.H4_p12_1e6"]),
            "TAB_tiered81 raw_input (read after results were written)", "%s / %s" % (res[q5]["value"], res[q6]["value"]),
            "%s / %s (external)" % (q5, q6), "NA", "Different estimand (discovery outcome GWAS vs the external resource); context for the §3 branch only.")
    # deterministic values from the same files (agreement is weak evidence)
    lead_notes = res["S3.RAD23B_T2D.lead_exposure"]["notes"]
    m2 = re.search(r"outcome p at this SNP ([0-9.e-]+)", lead_notes)
    add("S1-seen.DIAMANTE_p_at_RAD23B_proxy_rs7856112", "0.042", "§1 Seen: DIAMANTE resource p at the proxy",
        m2.group(1), "S3.RAD23B_T2D.lead_exposure notes (outcome p printed by DIAMANTE at rs7856112)", "%.3g" % abs(float(m2.group(1)) - 0.042),
        "Same printed value from the same file; deterministic, so agreement is weak evidence.")
    ext = {}
    with open(os.path.join(OUT, "extract_RAD23B_T2D.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["SNP"] == "9:110049914":
                ext = r
    add("S1-seen.DIAMANTE_p_at_RAD23B_lead_rs7035725", "0.038", "§1 Seen: DIAMANTE resource p at the lead",
        ext["Pvalue"], "looked up in this folder's extract_RAD23B_T2D.tsv after the full text was opened (not a results.tsv quantity: the lead is palindromic, DIAMANTE EAF %s / UKB-PPP 0.550, so §3 drops it)" % ext["EAF"],
        "%.3g" % abs(float(ext["Pvalue"]) - 0.038), "Same printed value from the same file; deterministic, weak evidence.")
    cin = {r["ID"]: r for r in csv.DictReader(open(os.path.join(OUT, "coloc_input_RAD23B_T2D.tsv"), encoding="utf-8"), delimiter="\t", quoting=csv.QUOTE_NONE)}
    b = cin["9:110072945:G:T:imp:v1"]
    assert b["ALLELE1"] == "T"
    add("S1-seen.UKBPPP_beta_rs7856112_T", "-0.0637", "§1 Seen: UKB-PPP at rs7856112 T beta", b["BETA"],
        "coloc_input_RAD23B_T2D.tsv BETA (ALLELE1 = T)", "%.2g" % abs(float(b["BETA"]) + 0.0637), "Same file value; deterministic, weak evidence.")

    cols = ["quantity_id", "seen_value", "seen_label", "stage1_value", "stage1_label", "abs_diff", "explanation"]
    with open(os.path.join(OUT, "COMPARE_SEEN.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            assert not any(("\t" in str(r[c]) or "\n" in str(r[c])) for c in cols)
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    log("wrote COMPARE_SEEN.tsv: %d rows (%d §3 quantities with no seen value; the rest are seen values from other analyses or deterministic file values)"
        % (len(rows), sum(1 for r in rows if r["seen_value"] == "none")))


if __name__ == "__main__":
    main()
