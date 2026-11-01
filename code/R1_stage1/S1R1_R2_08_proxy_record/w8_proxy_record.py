# -*- coding: utf-8 -*-
"""w8_proxy_record.py (Stage 7 A2.2, S2-D47): the proxy record for the 81-pair map. python -I -X utf8 w8_proxy_record.py

Writes   proxy_record.tsv           (5 rows, PUBLIC: no per-variant outcome beta, SE or p; provider terms do not establish redistribution, release/inputs/PROVIDER_TERMS.tsv FinnGen and MEGASTROKE rows)
         proxy_record_internal.tsv  (checks only: every input of the asserts, provider values included; NOT for the release)
Reads    only saved files, each pinned by md5 below (a changed file stops the script: do not re-pin, report). NOTHING IS TYPED: leads, proxies, alleles, positions, effects and r2 are read.
Order    1 pins and anchors -> 2 the 81 pairs and their proxy flag -> 3 LD rows (w8_export_ld.R) -> 4 outcome rows (w8_raw_scan.py, cross-checked against the saved extracts) ->
         5 Wald ratios at proxy and lead -> 6 asserts (i)-(vii) -> 7 outputs.
Orientation (assert v): the protein-raising (PR) allele of the lead is the instrument's effect allele if its UKB-PPP beta is > 0, else the other allele. The PLINK in-phase pair table lists the two
phased haplotypes (a_hap1,b_hap1) and (a_hap2,b_hap2); the proxy allele on the PR haplotype is the b_hap whose a_hap equals the PR allele. This is the rule of orient_proxy() in 01_mr_screen.R
(paper3_deposit_request_v4/.../cascade/01_mr_screen.R l.606: p_inphase <- if (ie==pr$a_hap1) pr$b_hap1 else if (ie==pr$a_hap2) pr$b_hap2). The outcome beta is the beta of the effect allele of the
outcome file (FinnGen: alt; MEGASTROKE: hm_effect_allele), sign-flipped when the proxy allele on the PR haplotype is the outcome file's other allele. The opposite haplotype is tested too and must NOT reproduce the screen.
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
import csv, hashlib, io, json, math, os, re, sys

V3 = _P31R_V3
HERE = V3 + "/R1_stage1/S1R1_R2_08_proxy_record"
WPB = V3 + "/R2_rewrite/review_R2/seen/wpb"
csv.field_size_limit(10 ** 9)

PINS = {
    "tier":    (V3 + "/R2_rewrite/release/out/colocalization/table_s_pair_evidence_tiered.csv", "f6ce42dae779294a041d848fc0822531"),
    "mr":      (V3 + "/V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv", "2d54dadda4da121a55bef600d6b00095"),
    "panel":   (V3 + "/R2_rewrite/release/out/replication/replication_sensitivity__panel1_2.csv", "c31efc6bceabccf7831bf369561f8043"),
    "ins":     (V3 + "/V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv", "7170afca7086d63918093b445a8fe6ff"),
    "ld_rds":  (V3 + "/paper3_deposit_request_v5/files/data/ld_reference/1kg_eur/instr_ld_proxies_phase1_hg38.rds", "9057d5006b1a1817fa223dd04ea88acc"),
    "ld_log":  (V3 + "/paper3_deposit_request_v5/files/data/ld_reference/1kg_eur/instr_ld_proxies_phase1_hg38.log", "ccf3757e15763408e5ef04309edd5ed0"),
    "ld_rows": (HERE + "/w8_ld_rows.tsv", "c1e586351627eee80b09a5a782dc933a"),
    "ld_pan":  (HERE + "/w8_ld_panel.txt", "efa2a9206f706c31428be63262b2a24f"),
    "raw":     (HERE + "/w8_outcome_rows_raw.tsv", "41118b6815671db664b5297a8ecb4a7b"),
    "fg_T2D":  (WPB + "/fg_T2D.txt", "7708ca9b4f3f892dea5714fd2f005b08"),
    "fg_CHD":  (WPB + "/fg_CHD.txt", "862e02cb1c435e7a14b6e1d2f9aea6c5"),
    "fg_AF":   (WPB + "/fg_AF.txt", "591ecdabdc9289031c6ced7584c37624"),
    "mega":    (WPB + "/mega.txt", "241f9d873bf5c9798ecef40d71b3ccb7"),
    "draft":   (V3 + "/R2_rewrite/review_R2/seen/wpb/proxy_ledger_draft.csv", "f6edc3bdb193ff592b8454d8d3217da8"),   # the saved draft; the release now ships proxy_record.tsv
    "ledger":  (V3 + "/R1_stage1/S1R1_06_claims_ledger/S1R1_CLAIMS_LEDGER.csv", None),   # rebuilt at every stage: NOT pinned, rows asserted by content (L668, L713, L714)
}
RAW_DONE_MD5 = {"FinnGen_R10_T2D": "00f930552bdff89d40a921377d944f68", "FinnGen_R10_I9_CHD": "8a5ea2ee65abd375d532b7cc79b35d7b",
                "FinnGen_R10_I9_AF": "aea6fc3288b9eba7a9a0200e471283ad", "MEGASTROKE_EUR_GCST006906": "e182e23019ee482cfdf266ad2b6c8751"}   # md5 of the raw outcome files, as recorded by w8_raw_scan.py
LOG = []


def say(msg):
    LOG.append(msg)
    print(msg, flush=True)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


for k, (p, m) in PINS.items():
    if m is not None:
        assert md5(p) == m, "source md5 differs from the pinned one (%s): %s; do not re-pin, report" % (p, md5(p))
say("pins: %d files match their md5" % sum(1 for v in PINS.values() if v[1]))


def rd(key, delim=","):
    with io.open(PINS[key][0], encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def close(a, b, rel):
    return abs(a - b) <= rel * max(abs(a), abs(b))


# ---- 2. the 81 pairs and the proxy flag
tier = rd("tier")
assert len(tier) == 81 and len({(r["gene"], r["outcome"]) for r in tier}) == 81
mr_multi = {}
for r in rd("mr"):
    mr_multi.setdefault((r["gene"], r["outcome_name"]), []).append(r)
assert all(len(mr_multi.get((r["gene"], r["outcome"]), [])) == 1 for r in tier), "an 81-map pair has no unique (gene, outcome_name) row in mr_all_results"   # other genes repeat (e.g. CAD in two resources); the 81 do not
mr = {k: v[0] for k, v in mr_multi.items() if len(v) == 1}
prox_pairs = [r for r in tier if mr[(r["gene"], r["outcome"])]["any_proxy"] == "TRUE"]
assert len(prox_pairs) == 5, len(prox_pairs)                                                  # (vii)
assert all(mr[(r["gene"], r["outcome"])]["method"] == "Wald ratio" and mr[(r["gene"], r["outcome"])]["n_snps"] == "1" and mr[(r["gene"], r["outcome"])]["n_proxy"] == "1" for r in prox_pairs)
say("(vii) exactly 5 proxy pairs among the 81 (any_proxy TRUE in mr_all_results; all single-instrument Wald ratios with n_proxy 1): " + ", ".join("%s-%s" % (r["gene"], r["outcome"]) for r in prox_pairs))
order = [(r["gene"], r["outcome"]) for r in prox_pairs]

# ---- 3. LD rows
ld = {r["pair"]: r for r in rd("ld_rows", "\t")}
assert sorted(ld) == sorted("%s-%s" % k for k in order), (sorted(ld), order)
pan = {(r["gene"], r["outcome_name"]): r for r in rd("panel")}
ins_all = rd("ins")
ins = {}
for r in ins_all:
    assert (r["SNP"], r["gene"]) not in ins
    ins[(r["SNP"], r["gene"])] = r
ld_panel_line = dict(l.split("\t", 1) for l in io.open(PINS["ld_pan"][0], encoding="utf-8").read().splitlines())["log_line"]
n_ld = int(re.match(r"(\d+) people", ld_panel_line).group(1))
assert re.search(r"--bfile .*1kg_eur/EUR", io.open(PINS["ld_log"][0], encoding="utf-8").read()) and "--r2 in-phase with-freqs" in io.open(PINS["ld_log"][0], encoding="utf-8").read()
ld_panel_label = "1000 Genomes European (n = %d), PLINK --r2 in-phase" % n_ld

# ---- 4. outcome rows (raw scan, by position), cross-checked against the saved extracts
raw = {}
done = {}
for line in io.open(PINS["raw"][0], encoding="utf-8").read().splitlines()[1:]:
    a = line.split("\t")
    if a[0] == "#DONE":
        done[a[1]] = dict(kv.split("=", 1) for kv in a[4:6])
        continue
    assert (a[0], a[1]) not in raw
    raw[(a[0], a[1])] = dict(zip(["outcome_id", "label", "chrom", "pos", "ref", "alt", "rsid", "p", "beta", "se", "eaf", "extra"], a))
assert set(done) == set(RAW_DONE_MD5) and all(done[o]["md5"] == RAW_DONE_MD5[o] for o in done), done
say("raw outcome files read by position: md5 recorded by w8_raw_scan.py equals the pinned md5 for all 4 files")
EXTRACT = {"FinnGen_R10_T2D": "fg_T2D", "FinnGen_R10_I9_CHD": "fg_CHD", "FinnGen_R10_I9_AF": "fg_AF"}
n_ext = 0
for oid, key in EXTRACT.items():
    for line in io.open(PINS[key][0], encoding="utf-8").read().splitlines():
        a = line.split("\t")
        # extract columns (FinnGen): chrom pos ref alt rsid gene p mlogp beta sebeta af_alt ...; the header is not in the extract, so the column roles are asserted by consistency
        z = float(a[8]) / float(a[9])
        assert abs(float(a[7]) + math.log10(float(a[6]))) < 1e-3, ("mlogp != -log10(p)", a)
        assert any(r["outcome_id"] == oid and r["pos"] == a[1] and r["chrom"] == a[0] and r["rsid"] == a[4] and r["ref"] == a[2] and r["alt"] == a[3] and r["beta"] == a[8] and r["se"] == a[9] and r["p"] == a[6]
                   for r in raw.values()), ("extract row not equal to the raw row", oid, a[:5])
        n_ext += 1
mh = io.open(PINS["mega"][0], encoding="utf-8").read().splitlines()
mrow = dict(zip(mh[0].split("\t"), mh[1].split("\t")))
assert any(r["outcome_id"] == "MEGASTROKE_EUR_GCST006906" and r["pos"] == mrow["hm_pos"] and r["chrom"] == mrow["hm_chrom"] and r["rsid"] == mrow["hm_rsid"] and r["ref"] == mrow["hm_other_allele"]
           and r["alt"] == mrow["hm_effect_allele"] and r["beta"] == mrow["hm_beta"] and r["se"] == mrow["standard_error"] for r in raw.values()), "mega.txt row not equal to the raw row"
n_ext += 1
say("saved extracts (fg_T2D, fg_CHD, fg_AF, mega) = raw-file rows for all %d extract rows (ref, alt, rsID, p, beta, SE)" % n_ext)


def aligned(ref, alt, beta, allele):
    """beta for `allele` given the beta of `alt` against `ref` (no strand flipping: an allele outside {ref, alt} stops)."""
    if allele == alt:
        return beta
    if allele == ref:
        return -beta
    raise AssertionError("allele %s is not in the outcome pair %s/%s" % (allele, ref, alt))


def draft_row(pair):
    for r in rd("draft"):
        if r["pair"] == pair:
            return r
    raise KeyError(pair)


# ledger anchors for the control (vi) and the count
led = {r["id"]: r for r in rd("ledger")}
assert json.loads(led["L668"]["value"]) == 5

pub, internal = [], []
for (gene, outc) in order:
    pair = "%s-%s" % (gene, outc)
    t = next(r for r in tier if r["gene"] == gene and r["outcome"] == outc)
    m = mr[(gene, outc)]
    p = pan[(gene, outc)]
    L = ld[pair]
    lead = t["lead_snps"]
    assert re.fullmatch(r"rs\d+", lead)
    ii = ins[(lead, gene)]
    # (iii) proxy ids equal proxy_snps of the panel file and the LD table
    mm = re.fullmatch(r"(rs\d+)<-(rs\d+)\(r2=([0-9.]+)\)", p["proxy_snps"])
    assert mm and mm.group(1) == lead and mm.group(2) == L["proxy_rsid"], (p["proxy_snps"], lead, L["proxy_rsid"])
    assert L["lead_rsid"] == lead and L["panel_instrument_id"] == ii["rsid_ld_panel"] and L["outcome_id"] == m["outcome_id"] == p["outcome_id"]
    prox = L["proxy_rsid"]
    # (iv) r2 equals the LD table: the panel file carries it at 4 dp in the text and at 6 dp in proxy_r2_min; the draft at 4 dp
    r2 = float(L["r2"])
    assert abs(r2 - float(p["proxy_r2_min"])) < 5e-7 and abs(round(r2, 4) - float(mm.group(3))) < 1e-12 and abs(round(r2, 4) - float(draft_row(pair)["r2_503"])) < 1e-12, (pair, r2)
    # instrument and PR allele
    bX = float(ii["beta"])
    ea, oa = ii["effect_allele"], ii["other_allele"]
    pr = ea if bX > 0 else oa
    bX_pr = abs(bX)
    chrom, lpos = ii["chrpos"].split(":")
    haps = [(L["a_hap1"], L["b_hap1"]), (L["a_hap2"], L["b_hap2"])]
    assert pr in {h[0] for h in haps} and len({h[0] for h in haps}) == 2 and len({h[1] for h in haps}) == 2, (pair, pr, haps)
    prox_pr = next(h[1] for h in haps if h[0] == pr)
    prox_other = next(h[1] for h in haps if h[0] != pr)
    # outcome rows
    oid = m["outcome_id"]
    pr_row = next(r for r in raw.values() if r["outcome_id"] == oid and r["label"] == "proxy:" + prox)
    ld_row = next((r for r in raw.values() if r["outcome_id"] == oid and r["label"] == "lead:" + lead), None)
    assert pr_row["chrom"] == chrom
    # the lead row is taken by POSITION (the instruments file's chrpos); a row at the position must carry the instrument's allele pair
    if ld_row is not None:
        assert ld_row["pos"] == lpos and {ld_row["ref"], ld_row["alt"]} == {ea, oa}, (pair, ld_row)
    assert {pr_row["ref"], pr_row["alt"]} == {prox_pr, prox_other}, (pair, pr_row, prox_pr, prox_other)
    bYp = aligned(pr_row["ref"], pr_row["alt"], float(pr_row["beta"]), prox_pr)
    wald_p = bYp / bX_pr
    bYp_alt = aligned(pr_row["ref"], pr_row["alt"], float(pr_row["beta"]), prox_other)
    wald_p_alt = bYp_alt / bX_pr
    # (i) equals the executed MR estimate within 1e-6 relative
    b_mr = float(m["b"])
    assert close(wald_p, b_mr, 1e-6), ("(i)", pair, wald_p, b_mr)
    # the opposite haplotype must not reproduce it
    assert not close(wald_p_alt, b_mr, 1e-3) and close(wald_p_alt, -b_mr, 1e-6), ("(v) opposite mapping", pair, wald_p_alt, b_mr)
    # (ii) equals the draft CSV
    dr = draft_row(pair)
    assert close(wald_p, float(dr["wald_proxy"]), 1e-6), ("(ii)", pair, wald_p, dr["wald_proxy"])
    # lead
    if ld_row is not None:
        bYl = aligned(ld_row["ref"], ld_row["alt"], float(ld_row["beta"]), pr)
        wald_l = bYl / bX_pr
        assert dr["wald_lead"] != "" and close(wald_l, float(dr["wald_lead"]), 1e-6), ("lead vs draft", pair, wald_l, dr["wald_lead"])
        lead_state = "yes, as " + ld_row["rsid"]
        # frequency sanity of the alignment at the lead: the instrument's effect-allele EAF and the outcome file's EAF for the same allele
        eaf_o = float(ld_row["eaf"]) if ld_row["alt"] == ea else 1 - float(ld_row["eaf"])
        eaf_gap = abs(float(ii["eaf"]) - eaf_o)
        if {ea, oa} in ({"A", "T"}, {"C", "G"}):                                             # a palindromic lead: the strand is only checkable by frequency; the gap must be small
            assert eaf_gap < 0.1, (pair, ii["eaf"], eaf_o)                                  # a non-palindromic lead is oriented by its alleles; the gap is recorded only (ABO: UKB vs Finnish frequency differ)
    else:
        bYl = wald_l = None
        lead_state = "no"
        assert dr["wald_lead"] == "", (pair, dr["wald_lead"])
        eaf_gap = None
    palin = {ea, oa} in ({"A", "T"}, {"C", "G"})
    amb = palin and (0.42 < float(ii["eaf"]) < 0.58)
    label = oid.replace("_", " ")
    # Stage 7 A2.2 (main session): why the screen read the outcome at a proxy. 01_mr_screen.R L670-684 routes a palindromic direct row whose EAF
    # lies inside 0.42-0.58 (exposure or outcome) to the proxy path, because harmonise_data(action=2) would discard it as ambiguous; the screen
    # matched outcome rows by rsID, so a lead filed under another identifier was not found directly. Exactly one reason per pair, or refuse.
    _alias = lead_state.startswith("yes, as ") and lead_state[len("yes, as "):] != lead
    _why = [w for w, c in (("palindromic lead whose allele frequency could not resolve its strand", amb),
                           ("lead listed under another identifier in the outcome file; the screen matched by rsID", (not amb) and _alias),
                           ("lead absent from the outcome file", lead_state == "no")) if c]
    assert len(_why) == 1, (pair, lead_state, amb, _why)
    pub.append(dict(pair=pair, outcome_resource=label, lead_rsid=lead, lead_pos_grch38=ii["chrpos"], protein_raising_allele=pr, proxy_rsid=prox, proxy_pos_grch38="%s:%s" % (pr_row["chrom"], pr_row["pos"]),
                    proxy_effect_allele_in_outcome_file=pr_row["alt"], proxy_other_allele_in_outcome_file=pr_row["ref"], proxy_allele_on_protein_raising_haplotype=prox_pr,
                    r2_with_lead="%.6f" % r2, ld_panel=ld_panel_label, denominator="UKB-PPP effect of the lead", lead_in_outcome_file=lead_state,
                    proxy_reason=_why[0], wald_at_proxy="%.10g" % wald_p, wald_at_lead=("" if wald_l is None else "%.10g" % wald_l)))
    internal.append(dict(pair=pair, outcome_id=oid, gene=gene, lead=lead, panel_instrument_id=L["panel_instrument_id"], lead_chrpos=ii["chrpos"], ukb_effect_allele=ea, ukb_other_allele=oa, ukb_beta=ii["beta"], ukb_se=ii["se"],
                         ukb_eaf=ii["eaf"], ukb_samplesize=ii["samplesize"], ukb_F=ii["F_stat"], protein_raising_allele=pr, bX_pr=repr(bX_pr), proxy=prox, proxy_chrpos="%s:%s" % (pr_row["chrom"], pr_row["pos"]),
                         haplotypes="%s-%s|%s-%s" % (haps[0] + haps[1]), proxy_allele_pr=prox_pr, outcome_proxy_ref=pr_row["ref"], outcome_proxy_alt=pr_row["alt"], outcome_proxy_beta_alt=pr_row["beta"],
                         outcome_proxy_se=pr_row["se"], outcome_proxy_p=pr_row["p"], outcome_proxy_eaf_alt=pr_row["eaf"], bY_proxy_aligned=repr(bYp),
                         outcome_lead_label=("" if ld_row is None else ld_row["rsid"]), outcome_lead_ref=("" if ld_row is None else ld_row["ref"]), outcome_lead_alt=("" if ld_row is None else ld_row["alt"]),
                         outcome_lead_beta_alt=("" if ld_row is None else ld_row["beta"]), outcome_lead_se=("" if ld_row is None else ld_row["se"]), outcome_lead_p=("" if ld_row is None else ld_row["p"]),
                         outcome_lead_eaf_alt=("" if ld_row is None else ld_row["eaf"]), bY_lead_aligned=("" if bYl is None else repr(bYl)),
                         wald_proxy=repr(wald_p), wald_lead=("" if wald_l is None else repr(wald_l)), mr_all_results_b=m["b"], rel_diff_vs_mr=repr(abs(wald_p - b_mr) / abs(b_mr)),
                         draft_wald_proxy=dr["wald_proxy"], draft_wald_lead=dr["wald_lead"], wald_proxy_opposite_haplotype=repr(wald_p_alt), r2_ld_table=L["r2"], r2_panel_proxy_r2_min=p["proxy_r2_min"],
                         lead_palindromic=str(palin), lead_exposure_eaf_in_0p42_0p58=str(0.42 < float(ii["eaf"]) < 0.58), screen_would_drop_direct_row_as_ambiguous=str(amb),
                         lead_eaf_gap_ukb_vs_outcome=("" if eaf_gap is None else "%.4f" % eaf_gap)))
    say("%-12s lead %s (%s) PR %s | proxy %s PR-haplotype allele %s | r2 %.6f | wald at proxy %.6g = mr_all b %.6g (rel %.1e); opposite haplotype %.6g | lead in file: %s%s" %
        (pair, lead, ii["chrpos"], pr, prox, prox_pr, r2, wald_p, b_mr, abs(wald_p - b_mr) / abs(b_mr), wald_p_alt, lead_state, "" if wald_l is None else "; wald at lead %.6g" % wald_l))

say("(i) wald_at_proxy = aligned bY_proxy / bX_lead equals mr_all_results b within 1e-6 relative: 5/5")
say("(ii) equals the draft CSV wald_proxy within 1e-6 relative: 5/5; wald_at_lead equals the draft wald_lead where the lead is present: %d/%d (ABO-Stroke blank in both)" % (sum(1 for r in pub if r["wald_at_lead"]), sum(1 for r in pub if r["wald_at_lead"])))
say("(iii) proxy ids equal proxy_snps of replication_sensitivity__panel1_2.csv and the LD table: 5/5")
say("(iv) r2: LD-table value equals proxy_r2_min (6 dp, abs < 5e-7), the panel text and the draft (4 dp): 5/5")
say("(v) orientation: the phased haplotype of the lead's protein-raising allele fixes the proxy allele; the opposite haplotype gives the negated ratio and fails (i): 5/5")

# (vi) control: RAD23B against ledger L713 / L714
c = next(r for r in pub if r["pair"] == "RAD23B-T2D")
l713, l714 = led["L713"], led["L714"]
assert c["proxy_rsid"] in l713["claim"] and c["proxy_rsid"] in l714["claim"] and c["lead_rsid"] in l714["claim"], (l713["claim"], l714["claim"])
r2_led = json.loads(l714["value"])
assert close(float(c["r2_with_lead"]), r2_led, 1e-6) and close(float(next(i for i in internal if i["pair"] == "RAD23B-T2D")["r2_ld_table"]), r2_led, 1e-6), (c["r2_with_lead"], r2_led)
say("(vi) control RAD23B-T2D: proxy %s and lead %s are named in ledger L713/L714; LD-table r2 %s vs L714 value %s: relative difference %.1e (PLINK prints 6 digits, L714 is the 8-digit recomputation in S1R1_R2_01_rad23b)" %
    (c["proxy_rsid"], c["lead_rsid"], next(i for i in internal if i["pair"] == "RAD23B-T2D")["r2_ld_table"], l714["value"], abs(float(next(i for i in internal if i["pair"] == "RAD23B-T2D")["r2_ld_table"]) - r2_led) / r2_led))
say("count: 5 = ledger L668 value (draft row count) = 5 any_proxy pairs among the 81")
say("lead present in the outcome file: " + "; ".join("%s: %s" % (r["pair"], r["lead_in_outcome_file"]) for r in pub))
say("direct row dropped as ambiguous by the screen's rule (palindromic lead, exposure EAF in 0.42-0.58): " + "; ".join("%s: %s" % (i["pair"], i["screen_would_drop_direct_row_as_ambiguous"]) for i in internal))


def write(path, rows):
    cols = list(rows[0].keys())
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            assert all("\t" not in str(r[c]) and "\n" not in str(r[c]) for c in cols)
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


write(HERE + "/proxy_record.tsv", pub)
write(HERE + "/proxy_record_internal.tsv", internal)
banned = ("beta", "se", "p", "bY", "outcome_proxy", "outcome_lead")
assert not [c for c in pub[0] if any(c == b or c.startswith(b + "_") for b in banned)], "a provider-value column reached the public file"
io.open(HERE + "/w8_proxy_record.log", "w", encoding="utf-8", newline="\n").write("\n".join(LOG) + "\n")
say("wrote proxy_record.tsv (%d rows, %d columns) and proxy_record_internal.tsv" % (len(pub), len(pub[0])))
