# -*- coding: utf-8 -*-
"""Stage 5 step 2 (written 2026-10-17; BUILT at step 1c): ledger rows L660-L689 for the numbers that the release README, REPRODUCE and the
release manifest summary will print and that had no ledger row (1b finding 7, STAGE5_STEP1AB_wf_1ac9282b.json, plus the plan's step 1d list,
STAGE5_PLAN_R2_2026-10-16.md, "Numbers the new README and REPRODUCE lines will print").

Executed by S1R1_build_ledger.py in ITS namespace (add, rel, ROWS, pd, np, V3, CKM, HERE, io, os, json are available). Every name defined here starts
with xr_ so nothing of the builder is overwritten. NOTHING IS TYPED: every value is read from the file named in the row, and the structural
anchors of each source (a column, a header sentence, a log line) are asserted, so a changed source fails the build instead of being miscounted.
Counts that the ledger already holds are cross-checked against those rows (L001, L004, L005, L400), never retyped.

Release counts are read from the DECLARED lists (R2_rewrite/release/inputs/RELEASE_INPUTS.tsv: columns role, release_path and ship_status only), never
from release/out/. The columns that `P31R_build_release.py --pin` rewrites (source_md5, source_bytes, cited_by) are not read: RELEASE_INPUTS.tsv pins the
ledger CSV, so a ledger row that read them would never settle on a fixed point.

Rows are packed (several leaves per row, stored in display order) so that no hand-declared order is needed in R2_d1_display_order.tsv:
  L660 stored screen rows | L661 non-primary rows, pairs | L662 proxied rows, pairs | L663 SuSiE grid rows | L664 deposited script files
  L665 .i8 files | L666 post-run edits in the deposit | L667 render units, figures | L668 proxy draft rows | L669 decoy-null inputs
  L670 files per workstream folder | L671 Table 2 rows | L672 branches_integrated rows per source | L673 Table 1 rows
  L674 untestable combinations by recorded reason | L675 instrument genes without a screen row | L676 providers in the release provider list
L670-L675 beyond the two lists named above are candidates the step 2 builders listed (Table 1 rows; the reason-code split and the absent genes);
delete a row here if the README will not print it. The PROVIDER COUNT (1b finding 7) is L676: it is read from release/inputs/PROVIDERS.tsv, which
`P31R_build_release.py --providers` writes from the declared provider tables and which a normal build regenerates in memory and asserts equal (fix round
2026-10-18, user decision 2). The Declarations sentence that holds 14 changes at step 5 (release/STEP2_RECORD.md, open question 1).
"""
import csv as xr_csv, io as xr_io, json as xr_json, os as xr_os, re as xr_re, hashlib as xr_hl

xr_csv.field_size_limit(10 ** 9)
xr_LED = {r["id"]: xr_json.loads(r["value"]) for r in ROWS}      # rows built so far (base rows and earlier extension modules)
xr_WORD = {w: i for i, w in enumerate("one two three four five six seven eight nine".split(), 1)}


def xr_tsv(path):
    with xr_io.open(path, encoding="utf-8", newline="") as fh:
        return list(xr_csv.DictReader(fh, delimiter="\t"))


def xr_txt(path):
    return xr_io.open(path, encoding="utf-8", errors="replace").read()


def xr_n(x):
    return "{:,}".format(x)


xr_RW = xr_os.path.join(V3, "R2_rewrite")
xr_DEP = xr_os.path.join(V3, "paper3_deposit_request_v4")
xr_TFR = xr_os.path.join(V3, "V4", "source_data", "tables_from_rds")
xr_INT = xr_os.path.join(V3, "R1_stage1", "S1R1_R2_10_integration")

# ---------------------------------------------------------------------------------------------------- the stored screen output
xr_MR = xr_os.path.join(xr_TFR, "phase2_mr_screen__mr_all_results.csv")
xr_mr = pd.read_csv(xr_MR)
assert xr_mr["any_proxy"].dtype == bool, "any_proxy is not parsed as a boolean column: %s" % xr_mr["any_proxy"].dtype
xr_KEY = ["gene", "outcome_name"]
xr_prim = xr_mr[xr_mr.method.isin(["Wald ratio", "Inverse variance weighted"])]   # the ledger's own definition (S1R1_build_ledger.py, L004)
xr_nonprim = xr_mr.drop(xr_prim.index)
xr_pairs = xr_mr[xr_KEY].drop_duplicates()
assert len(xr_prim) == len(xr_pairs) == xr_LED["L004"], "primary rows / unique pairs differ from the ledger's testable count L004"
assert int((xr_prim.method == "Wald ratio").sum()) == xr_LED["L005"] and int((xr_prim.method == "Inverse variance weighted").sum()) == xr_LED["L400"]["ivw"]
xr_np_pairs = xr_nonprim[xr_KEY].drop_duplicates()
assert len(xr_nonprim) == len(xr_mr) - xr_LED["L004"], "non-primary rows are not the surplus over the unique pairs"
assert set(xr_nonprim.method) == {"MR Egger", "Weighted median"}, "the claim of L661 names MR Egger and weighted median: %s" % sorted(set(xr_nonprim.method))
assert len(xr_np_pairs.merge(xr_prim[xr_KEY], on=xr_KEY)) == len(xr_np_pairs), "a non-primary row belongs to a pair that has no primary row"
xr_px = xr_mr[xr_mr.any_proxy]
assert len(xr_px) and len(xr_px[xr_KEY].drop_duplicates()) <= xr_LED["L004"]

add("L660", "Release: screen", "Rows in the stored MR screen output (every estimator; the release screen file has one row per protein x outcome pair, not per stored row)",
    int(len(xr_mr)), xr_n(len(xr_mr)), rel(xr_MR), "row count of phase2_mr_screen__mr_all_results.csv", "[D]",
    note="Stage 5 step 2 (1b finding 7): the README prints the stored row count beside the family size (L003) and the testable count (L004 = the unique pairs).")
add("L661", "Release: screen", "Stored screen rows that are not the primary estimate of their pair (MR Egger and weighted median), and the protein x outcome pairs they belong to",
    [int(len(xr_nonprim)), int(len(xr_np_pairs))], "%d non-primary rows from %d pairs" % (len(xr_nonprim), len(xr_np_pairs)), rel(xr_MR),
    "rows with method not in {Wald ratio, IVW}; unique (gene, outcome_name) among them", "[D]",
    note="Stage 5 step 2 (1b finding 6): surplus rows of a few multi-instrument pairs; their primary is the IVW estimate counted in L400. Shipped in a separate file.")
add("L662", "Release: screen", "Stored screen rows that carry a proxy instrument variant (any_proxy), and the protein x outcome pairs they belong to",
    [int(len(xr_px)), int(len(xr_px[xr_KEY].drop_duplicates()))], "%s proxied rows in %s proxied pairs" % (xr_n(len(xr_px)), xr_n(len(xr_px[xr_KEY].drop_duplicates()))), rel(xr_MR),
    "rows with any_proxy TRUE; unique (gene, outcome_name) among them", "[D]",
    note="Stage 5 step 2 (1b finding 6): the plan's '745 proxied rows' counts pairs; the rows are the larger number (the non-primary rows are all proxied).")

xr_SUSIE = xr_os.path.join(xr_TFR, "phase3_susie__susie_grid_results.csv")
xr_su = pd.read_csv(xr_SUSIE)
add("L663", "Release: finemapping", "Rows of the stored SuSiE grid results", int(len(xr_su)), xr_n(len(xr_su)), rel(xr_SUSIE),
    "row count of phase3_susie__susie_grid_results.csv", "[D]", note="Stage 5 step 2 (1b finding 7).")

# ---------------------------------------------------------------------------------------------------- the deposited code (MANIFEST.tsv of deposit v4)
xr_MAN = xr_os.path.join(xr_DEP, "MANIFEST.tsv")
xr_man = xr_tsv(xr_MAN)
assert set(xr_man[0]) >= {"relative_path", "md5_source"}
xr_pa = [r["relative_path"] for r in xr_man if r["relative_path"].startswith("scripts/scripts/")]
xr_pb = [r["relative_path"] for r in xr_man if r["relative_path"].startswith("scripts/R/")]
xr_ext = {}
for xr_p in xr_pa + xr_pb:
    xr_e = xr_os.path.splitext(xr_p)[1]
    xr_ext[xr_e] = xr_ext.get(xr_e, 0) + 1
# F13: the manifest rows OUTSIDE the two prefixes are described from the manifest itself (code-extension files, their folders, byte copies of in-scope files)
xr_in = set(xr_pa + xr_pb)
xr_CODEEXT = (".R", ".py", ".sh", ".bat")
xr_outc = [r for r in xr_man if r["relative_path"] not in xr_in and xr_os.path.splitext(r["relative_path"])[1] in xr_CODEEXT]
xr_oext, xr_odir = {}, {}
for xr_r in xr_outc:
    xr_e = xr_os.path.splitext(xr_r["relative_path"])[1]
    xr_oext[xr_e] = xr_oext.get(xr_e, 0) + 1
    xr_d = "/".join(xr_r["relative_path"].split("/")[:2])
    xr_odir[xr_d] = xr_odir.get(xr_d, 0) + 1
xr_inmd5 = {r["md5_source"] for r in xr_man if r["relative_path"] in xr_in}
xr_dupn = sum(1 for r in xr_outc if r["md5_source"] in xr_inmd5)
add("L664", "Release: code", "Files of the deposited code listed in the deposit v4 MANIFEST under scripts/scripts and scripts/R (total, scripts/scripts, scripts/R)",
    [len(xr_pa) + len(xr_pb), len(xr_pa), len(xr_pb)],
    "%d files (%d in scripts/scripts, %d in scripts/R)" % (len(xr_pa) + len(xr_pb), len(xr_pa), len(xr_pb)), rel(xr_MAN),
    "MANIFEST.tsv relative_path prefixes scripts/scripts/ and scripts/R/", "[D]",
    note=("Stage 5 step 2 (1b finding 7). In scope, by extension: %s. The manifest lists %d further code-extension files outside the two prefixes (by extension: %s; by folder: %s); "
          "byte copies (md5_source) of in-scope files among them: %d. Those further files are NOT counted in this row. The remaining manifest rows are not code-extension files." %
          (", ".join("%s %d" % (k, v) for k, v in sorted(xr_ext.items())), len(xr_outc), ", ".join("%s %d" % (k, v) for k, v in sorted(xr_oext.items())),
           ", ".join("%s %d" % (k, v) for k, v in sorted(xr_odir.items())), xr_dupn)))

# ---------------------------------------------------------------------------------------------------- individual-level dosage files under V3
xr_i8 = 0
for xr_root, xr_dirs, xr_files in xr_os.walk(V3):
    xr_dirs[:] = [d for d in xr_dirs if d not in (".git", "__pycache__")]
    xr_i8 += sum(1 for f in xr_files if f.endswith(".i8"))
assert xr_i8 > 0
add("L665", "Release: not shipped", "Individual-level dosage matrices (.i8 files) found under V3; none is released",
    int(xr_i8), "%d" % xr_i8, "V3/ (walk; file name suffix .i8)", "count of files ending .i8 (the same walk as the release not-shipped list, reason I8_DOSAGE)", "[D]",
    note="Stage 5 step 2 (1b finding 7, plan l.157-160).")

# ---------------------------------------------------------------------------------------------------- post-run edits recorded by the deposit REPORT
xr_REP = xr_os.path.join(xr_DEP, "REPORT.md")
xr_rep = xr_txt(xr_REP).splitlines()
xr_l6 = [l for l in xr_rep if l.startswith("6. **Three Part B facts")]
assert len(xr_l6) == 1
xr_m6 = xr_re.search(r"\(a\) (\w+) writer scripts were edited after their output was written", xr_l6[0])
assert xr_m6 and xr_m6.group(1).lower() in xr_WORD
xr_l1 = [l for l in xr_rep if l.startswith("1. **Writer scripts whose current bytes are not the version that produced the delivered file**")]
assert len(xr_l1) == 1
xr_scripts = sorted(set(xr_re.findall(r"`([\w./-]+\.(?:R|py))`", xr_l1[0].split("For A9")[0])))
assert len(xr_scripts) == xr_WORD[xr_m6.group(1).lower()], "the scripts listed in Part B (%s) do not match the count spelled in the first section (%s)" % (xr_scripts, xr_m6.group(1))
xr_relab = sorted({xr_re.match(r"\| (A\d+) \|", l).group(1) for l in xr_rep if xr_re.match(r"\| A\d+ \|", l) and "hand relabel" in l and "no writer script" in l})
assert len(xr_relab) == 1 and "(b)" in xr_l6[0] and ("%s intermediate_gwasinfo.csv" % xr_relab[0]) in xr_l6[0] and "hand relabel" in xr_l6[0]
add("L666", "Release: reproduce", "Deposit v4 writer scripts edited after their output was written, and delivered inputs that are a hand relabel with no writer script",
    [len(xr_scripts), len(xr_relab)], "%d writer scripts edited after their output was written; %d hand-relabelled input (%s)" % (len(xr_scripts), len(xr_relab), xr_relab[0]),
    rel(xr_REP), "REPORT.md Part B (writer scripts: %s); Part B table rows with 'hand relabel' and 'no writer script' (%s)" % (", ".join(xr_scripts), ", ".join(xr_relab)), "[D]",
    note="Stage 5 step 2 (1b finding 7): REPRODUCE states that the current md5 of these scripts is not the version that ran. The spelled count in section 0 item 6 agrees with the scripts listed in Part B.")

# ---------------------------------------------------------------------------------------------------- render units
xr_ACT = xr_os.path.join(xr_RW, "figures", "af3", "P31R_ACTIVE_FIGURES.json")
xr_act = xr_json.load(xr_io.open(xr_ACT, encoding="utf-8"))["figures"]
xr_units = sum(int(f["pages"]) for f in xr_act)
assert xr_units >= len(xr_act) and len({f["id"] for f in xr_act}) == len(xr_act)
add("L667", "Release: figures", "Render units of the active figure list (a figure of two pages is two units), and the figures in the list",
    [xr_units, len(xr_act)], "%d render units from %d figures" % (xr_units, len(xr_act)), rel(xr_ACT), "sum of pages over figures; number of figures", "[D]",
    note="Stage 5 step 2 (1b finding 7): the unit count of the single step 1c render.")

# ---------------------------------------------------------------------------------------------------- proxy draft
xr_PX = xr_os.path.normpath(xr_os.path.join(xr_RW, "..", "R1_stage1", "S1R1_R2_08_proxy_record", "proxy_record.tsv"))   # Stage 7 A2.6: the shipped proxy record
xr_pxd = pd.read_csv(xr_PX, sep="\t")   # the proxy record is tab-separated
add("L668", "Release: proxy", "Rows of the proxy record shipped in the release (built from the saved sources by w8_proxy_record.py, S2-D47)", int(len(xr_pxd)), "%d" % len(xr_pxd), rel(xr_PX),
    "row count of proxy_record.tsv", "[D]", note="Stage 5 step 2 (1b finding 7).")

# ---------------------------------------------------------------------------------------------------- decoy-null inputs named in the decoy log
xr_LOG = xr_os.path.join(V3, "R1_stage1", "S1R1_03_decoy_method", "S1R1_decoy_rates.log")
xr_dec = [xr_re.match(r"^md5 (\S+)\s+([0-9a-f]{32})\s+(\S+)$", l) for l in xr_txt(xr_LOG).splitlines() if l.startswith("md5 ")]
assert xr_dec and all(xr_dec), "a decoy log md5 line does not match the expected form"
add("L669", "Release: decoy", "Inputs of the decoy-rate script that stay out of the release, with the md5 each carries in the script's log",
    len(xr_dec), "%d" % len(xr_dec), rel(xr_LOG), "lines 'md5 <name> <md5> <path>' of S1R1_decoy_rates.log", "[D]",
    note="Stage 5 step 2 (plan l.263, l.269). Logged md5 (name, first 8 digits): %s." % "; ".join("%s %s" % (m.group(1), m.group(2)[:8]) for m in xr_dec))

# ---------------------------------------------------------------------------------------------------- files per workstream folder (declared lists)
xr_RI = xr_os.path.join(xr_RW, "release", "inputs", "RELEASE_INPUTS.tsv")
xr_ri = xr_tsv(xr_RI)
assert set(xr_ri[0]) >= {"role", "release_path", "ship_status"}
xr_ws = {}
for xr_r in xr_ri:
    if xr_r["role"] == "copy" and xr_r["ship_status"] == "ship" and xr_r["release_path"].startswith("workstreams/"):
        xr_f = xr_r["release_path"].split("/")[1]
        xr_ws[xr_f] = xr_ws.get(xr_f, 0) + 1
assert xr_ws
add("L670", "Release: workstreams", "Files released per workstream folder of the repository (total, then one count per folder in declared order)",
    [sum(xr_ws.values())] + list(xr_ws.values()),
    "%d files across the workstream folders: %s" % (sum(xr_ws.values()), ", ".join("%s %d" % (k, v) for k, v in xr_ws.items())), rel(xr_RI),
    "RELEASE_INPUTS.tsv rows with role copy, ship_status ship and release_path under workstreams/ (columns role, release_path, ship_status only)", "[R]",
    note="Stage 5 step 2 (plan l.269): read from the declared include list, which is also the list of files the builder copies. Folder ids as Table 3's Source-id column.")

# ---------------------------------------------------------------------------------------------------- Table 1 and Table 2 rows
xr_T1 = xr_os.path.join(xr_RW, "tables", "P31R_Table1.md")
xr_T2 = xr_os.path.join(xr_RW, "tables", "P31R_Table2.md")


def xr_body_rows(path):
    """Body rows of the one markdown table in a table file: the pipe lines after the |---| separator."""
    lines = xr_txt(path).splitlines()
    sep = [i for i, l in enumerate(lines) if xr_re.fullmatch(r"\|(?:\s*:?-{3,}:?\s*\|)+", l.strip())]
    assert len(sep) == 1, "%s: expected one table separator, found %d" % (path, len(sep))
    cols = lines[sep[0] - 1].count("|") - 1
    body = []
    for l in lines[sep[0] + 1:]:
        if not l.startswith("|"):
            break
        body.append(l)
    assert body and all(b.count("|") >= cols + 1 for b in body)
    return body


xr_b1, xr_b2 = xr_body_rows(xr_T1), xr_body_rows(xr_T2)
add("L671", "Release: tables", "Rows of Table 2 (evidence cards), hence of Table2_evidence_cards.tsv", len(xr_b2), "%d" % len(xr_b2), rel(xr_T2),
    "body rows of the markdown table in P31R_Table2.md", "[R]", note="Stage 5 step 2 (plan l.269, S2-D29 gate (d)): the release builder asserts the TSV has this many rows.")

# ---------------------------------------------------------------------------------------------------- branches_integrated rows per source
xr_BI = xr_os.path.join(xr_INT, "branches_integrated.tsv")
xr_bi = xr_tsv(xr_BI)
assert "source" in xr_bi[0]
xr_bs = {}
for xr_r in xr_bi:
    xr_bs[xr_r["source"]] = xr_bs.get(xr_r["source"], 0) + 1
add("L672", "Release: workstreams", "Rows of branches_integrated.tsv (total, then one count per source workstream in file order)", [len(xr_bi)] + list(xr_bs.values()),
    "%d rows: %s" % (len(xr_bi), ", ".join("%s %d" % (k, v) for k, v in xr_bs.items())), rel(xr_BI),
    "row count; rows per value of the source column (order of first appearance)", "[D]",
    note="Stage 5 step 2 (plan l.269): Table 3 cites 'BI <source> r<source_row>', so a reader needs the rows per source.")

add("L673", "Release: tables", "Rows of Table 1 (data sources), hence of Table1_data_sources.tsv", len(xr_b1), "%d" % len(xr_b1), rel(xr_T1),
    "body rows of the markdown table in P31R_Table1.md", "[R]",
    note="Stage 5 step 2: ADDED beyond the two lists (the README names the file beside Table2_evidence_cards.tsv); delete this row if the README does not print it.")

# ---------------------------------------------------------------------------------------------------- untestable combinations by recorded reason; instrument genes without a screen row
xr_RL = xr_os.path.join(xr_DEP, "scripts", "reports", "Paper3_Program", "run_logs", "REGEN_20260817", "01_mr_screen.log")
xr_cl = [l for l in xr_txt(xr_RL).splitlines() if l.startswith("MR cells:")]
assert len(xr_cl) == 1
xr_mc = xr_re.fullmatch(r"MR cells: (\d+) with-result, (\d+) attempted-no-result, (\d+) no-outcome-data, (\d+) null-exposure-gene \| full-grid (\d+)/(\d+) COMPLETE", xr_cl[0])
assert xr_mc, xr_cl[0]
xr_wr, xr_anr, xr_nod, xr_nul, xr_g1, xr_g2 = [int(x) for x in xr_mc.groups()]
assert xr_wr == xr_LED["L004"] and xr_g1 == xr_g2 == xr_LED["L003"], "the run log's cell totals differ from the ledger (L003, L004)"
assert xr_anr + xr_nod + xr_nul == xr_LED["L400"]["untestable"], "attempted + no-outcome + null-exposure cells differ from the untestable count in L400"
add("L674", "Release: screen", "Untestable combinations by the reason the run log records: no outcome data for any instrument variant, attempted with no estimate",
    [xr_nod, xr_anr], "%s with no outcome data; %d attempted with no estimate" % (xr_n(xr_nod), xr_anr), rel(xr_RL),
    "log line 'MR cells: ...' (with-result, attempted-no-result, no-outcome-data, null-exposure-gene); their sum is asserted equal to L400's untestable count",
    "[D]", note="Stage 5 step 2: ADDED beyond the two lists (the screen codebook prints it); null-exposure-gene cells: %d." % xr_nul)

xr_INS = xr_os.path.join(xr_TFR, "phase1_instruments__ukb_cis_instruments.csv")
xr_ins = pd.read_csv(xr_INS, usecols=["gene"])
xr_absent = set(xr_ins.gene) - set(xr_mr.gene)
assert len(set(xr_ins.gene)) == xr_LED["L001"] and set(xr_mr.gene) <= set(xr_ins.gene)
add("L675", "Release: screen", "Instrument genes with no row in the stored screen output, of all instrument genes, and the genes the stored screen holds",
    [len(xr_absent), len(set(xr_ins.gene)), int(xr_mr.gene.nunique())],
    "%d of %s instrument genes have no stored screen row (the stored screen holds %s genes)" % (len(xr_absent), xr_n(len(set(xr_ins.gene))), xr_n(xr_mr.gene.nunique())),
    rel(xr_INS) + " + " + rel(xr_MR), "set difference of gene between the instrument list and the stored screen; unique genes", "[D]",
    note="Stage 5 step 2: ADDED beyond the two lists. The family (L003) is instrument genes x outcomes, so these genes appear in the release screen file as untestable rows entered at p = 1 (L400).")

# ---------------------------------------------------------------------------------------------------- providers (fix round 2026-10-18, user decision 2)
# Stage 6.3b (2026-10-20): L676 is a predecessor, so it reads the provider list as it stood at the 1c ledger batch (36 rows, terms not yet read),
# frozen from release/_pre_S5S34_20261018/inputs/PROVIDERS.tsv and md5-pinned; its successor (ext_r2_stage6.py) reads the live list.
xr_PRV = xr_os.path.join(HERE, "frozen_inputs", "PROVIDERS_at_1c_L676.tsv")
assert xr_hl.md5(open(xr_PRV, "rb").read()).hexdigest() == "e05944719a56be320717c6fa67f4f43d", "the frozen 1c provider list changed"
xr_prv = xr_tsv(xr_PRV)
assert set(xr_prv[0]) >= {"provider_id", "in_declarations", "terms_status"} and len({r["provider_id"] for r in xr_prv}) == len(xr_prv)
assert all(r["terms_status"] == "not_read" for r in xr_prv), "the provider list records no provider terms as read"
xr_indecl = sum(1 for r in xr_prv if r["in_declarations"] == "yes")
add("L676", "Release: providers", "Data providers in the release provider list (Table 1 resources, Declarations acknowledgements, providers named in the operative surfaces and in the licence classes), and how many of them the Declarations acknowledgements already name",
    [len(xr_prv), xr_indecl], "%d providers, %d named in the Declarations acknowledgements; terms read for none" % (len(xr_prv), xr_indecl), rel(xr_PRV),
    "rows of PROVIDERS.tsv; rows with in_declarations yes (columns provider_id, in_declarations, terms_status only)", "[R]",
    note="Fix round 2026-10-18 (1b finding 7; user decision 2): replaces the Declarations count of 14 once the pending providers are added to the Declarations at step 5. PROVIDERS.tsv is written by P31R_build_release.py --providers.")
