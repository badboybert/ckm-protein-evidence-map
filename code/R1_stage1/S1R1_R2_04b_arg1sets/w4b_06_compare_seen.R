# w4b_06_compare_seen.R -- Stage 1 W4b, step 6: COMPARE_SEEN.tsv. Runs ONLY after results.tsv / branches.tsv were frozen
# (run.log FREEZE line). Stage 1 values are read from results.tsv (never retyped); seen values are transcribed from the full
# pre-declaration P31R_R2_STAGE1_PRESPEC.md (md5 3319dc39...) section 4 'Seen' and Rule 3(iii)/4(v) text, and from the
# registered grid columns (definition files read before Stage 1 computed s). No computed value is changed here.
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets/w4b_common.R"))
logf("==== W4b step 6 (COMPARE_SEEN) start")
assert_md5(GOV$PRESPEC_FULL$path, GOV$PRESPEC_FULL$md5, "PRESPEC_FULL (opened for COMPARE_SEEN)")
assert_md5(GOV$AMEND1_FULL$path, GOV$AMEND1_FULL$md5, "AMEND1_FULL (opened for COMPARE_SEEN)")
RES <- fread(file.path(OUT, "results.tsv"), colClasses = "character")
D <- readRDS(file.path(OUT, "data", "w4b_data.rds")); CSR <- readRDS(file.path(OUT, "data", "w4b_credsets.rds"))
span_kb <- function(k) { p <- D$H[match(CSR$FULL[[k]]$members, snp)]$pos38; (max(p) - min(p)) / 1e3 }
BT <- fread(file.path(OUT, "boot", "boot_credsets_long.tsv.gz"))[method == "m3" & signal == "A"]
rq <- quantile(BT$r_AB, c(0.025, 0.975))
v <- function(id) { x <- RES[quantity_id == id]$value; stopifnot(length(x) == 1); x }
num <- function(id) as.numeric(v(id))
rows <- list()
cmp <- function(qid, seen_value, seen_label, stage1_value, stage1_label, abs_diff, explanation)
  rows[[length(rows) + 1]] <<- data.table(quantity_id = qid, seen_value = seen_value, seen_label = seen_label, stage1_value = stage1_value,
                                          stage1_label = stage1_label, abs_diff = abs_diff, explanation = explanation)
SAME <- "same method and inputs as the scoping (gls_q.R re-implements seen/wpD/wpD_03_multisignal.py; same harmonised table; 503 LD from PLINK equals the R1 i8 dosage r to 4e-16): agreement expected by construction and is weak evidence"
d <- function(a, b) formatC(abs(a - b), digits = 4, format = "g")
# ---- 4(iv) five-signal GLS (full text section 4 Seen: 'arginine Q 0.33 (4 df; upper p 0.99; lower-tail p 0.012); T2D Q 6.51 (p 0.164)')
cmp("S4.omni.arginine.5sig.aware.Q", "0.33", "4 Seen: five-signal GLS, arginine Q (cov/se_order not stated)", v("S4.omni.arginine.5sig.aware.Q"), "aware, 2nd order, 503",
    d(0.33, num("S4.omni.arginine.5sig.aware.Q")), paste("agrees at the printed precision (2 dp);", SAME))
cmp("S4.omni.arginine.5sig.aware.df", "4", "4 Seen: df", v("S4.omni.arginine.5sig.aware.df"), "df", "0", "identical")
cmp("S4.omni.arginine.5sig.p_aware_2nd", "0.99", "4 Seen: arginine upper p", v("S4.omni.arginine.5sig.p_aware_2nd"), "aware 2nd-order upper p (primary)",
    d(0.99, num("S4.omni.arginine.5sig.p_aware_2nd")), paste("agrees at the printed precision;", SAME))
cmp("S4.omni.arginine.5sig.aware.p_lower", "0.012", "4 Seen: arginine lower-tail p", v("S4.omni.arginine.5sig.aware.p_lower"), "aware 2nd-order lower-tail p",
    d(0.012, num("S4.omni.arginine.5sig.aware.p_lower")), paste("agrees at the printed precision (2 s.f.);", SAME))
cmp("S4.omni.T2D_FG.5sig.aware.Q", "6.51", "4 Seen: five-signal GLS, T2D Q (resource not stated; FinnGen R10 is the declared T2D resource)", v("S4.omni.T2D_FG.5sig.aware.Q"),
    "FinnGen R10, aware, 2nd order, 503", d(6.51, num("S4.omni.T2D_FG.5sig.aware.Q")), paste("agrees at the printed precision;", SAME))
cmp("S4.omni.T2D_FG.5sig.p_aware_2nd", "0.164", "4 Seen: T2D upper p", v("S4.omni.T2D_FG.5sig.p_aware_2nd"), "FinnGen R10 aware 2nd-order upper p (primary)",
    d(0.164, num("S4.omni.T2D_FG.5sig.p_aware_2nd")), paste("agrees at the printed precision (3 dp);", SAME))
cmp("S4.omni.T2D_DM.5sig.*; S4.omni.*.6sig*", "none", "no seen value", "see results.tsv", "first computed in Stage 1", "", "no seen comparator (DIAMANTE omnibus and the six-signal sensitivity were not in the Seen block)")
# ---- s (full text: 's: ARG1 0.139-0.144 across the four 1000G-reference fits (503 and 489 panels; the maximum, 0.144, is from p1fr)')
sP <- c(susie = num("S4.s.P__ARG1.susie.503"), p1fr = num("S4.s.P__ARG1.p1fr.503"), p1np = num("S4.s.P__ARG1.p1np.503"), p2np = num("S4.s.P__ARG1.p2np.489"))
cmp("S4.s.P__ARG1.<set>.<panel> (range)", "0.139-0.144", "4 Seen: ARG1 s across the four fits", sprintf("%.4f-%.4f", min(sP), max(sP)), "Stage 1 range over susie/p1fr/p1np (503) and p2np (489)",
    sprintf("%.4f / %.4f", abs(0.139 - min(sP)), abs(0.144 - max(sP))), "agrees at the printed precision (3 dp); same estimator (susieR estimate_s_rss) and n_hat as the engine, so agreement is expected and weak evidence; Stage 1 LD from PLINK, SNP sets re-built in R")
cmp("S4.s.max.ARG1_T2D (argmax set)", "0.144 (p1fr)", "4 Seen: maximum and its set", sprintf("%s (p1fr)", v("S4.s.P__ARG1.p1fr.503")), "Stage 1 maximum is P__ARG1 on p1fr", d(0.144, num("S4.s.P__ARG1.p1fr.503")), "agrees; same set attains the maximum")
grid <- list(c("P__ARG1", "p1np", "0.1391"), c("P__ARG1", "p1fr", "0.1438"), c("arginine", "p1np", "0.0453"), c("arginine", "p1fr", "0.0474"),
             c("ornithine", "p1np", "0.0419"), c("ornithine", "p1fr", "0.0441"))
for (g in grid) {
  st <- num(sprintf("S4.s.%s.%s.503", g[1], g[2])); dd <- abs(as.numeric(g[3]) - st)
  cmp(sprintf("S4.s.%s.%s.503", g[1], g[2]), g[3], sprintf("registered grid column s_%s (pal=%s; engine run, 4 dp)", if (g[1] == "P__ARG1") "protein" else "metabolite", if (g[2] == "p1np") "drop" else "freq"),
      formatC(st, digits = 4, format = "f"), "Stage 1 s on the rebuilt set", formatC(dd, digits = 2, format = "g"),
      if (round(st, 4) == as.numeric(g[3])) "identical at 4 dp (p1np is fully determined by the harmonised table; engine settings re-implemented)"
      else "differs beyond the printed precision: Stage 1's p1fr set is rebuilt with the palindrome frequency rule over the traits fitted here (ARG1, arginine, ornithine and FinnGen T2D allele frequencies) and without CCN2; the engine's p1fr used every trait at the locus (including CCN2, not FinnGen). Different palindromes enter the set (Stage 1 3,045 SNPs); explained by the SNP-set definition")
}
cmp("4(v)(1) admissibility (Rule 3(i))", "inadmissible, s 0.14", "full text 4(v)(1) scoping classification", sprintf("inadmissible, s_max %s", formatC(num("S4.s.max.ARG1_T2D"), digits = 4, format = "f")),
    "Stage 1 Rule 3(i) on own s", "", "no CLASS_FLIP")
cmp("Rule 3(i) ARG1-arginine / ARG1-ornithine", "not stated", "the full text classifies only method (1)", sprintf("inadmissible (s_max %s / %s)", formatC(num("S4.s.max.ARG1_arginine"), digits = 4, format = "f"), formatC(num("S4.s.max.ARG1_ornithine"), digits = 4, format = "f")),
    "Stage 1", "", "no scoping classification to flip; s_max is the protein's, as for method (1)")
# ---- credible sets (full text: '(1) and (2): A = {rs2781668, rs2246012} in the ARG1 transcript; B = {rs17788484}; (3): 6 variants in MED23, excluding rs2781668; (4): 33 variants over 157 kb (top PIP 0.085)')
cmp("S4.cs.2.A.members", "rs2781668;rs2246012", "4 Seen: credible set (2) at A", v("S4.cs.2.A.members"), "ST16 CS rsIDs", "0", "identical; read from the same deposited table (deterministic)")
cmp("S4.cs.1_protein.A.members", "rs2781668;rs2246012", "4 Seen: credible set (1) at A (the seen (1) is the protein SuSiE set)", v("S4.cs.1_protein.A.members"), "ARG1 protein SuSiE CS at A, p1np 503", "0",
    "identical members; same estimator and settings as the engine run, so agreement is expected")
cmp("S4.cs.1_T2D.A.members", "not in the seen record", "the seen (1) reads only the protein credible set", sprintf("%s variants", RES[quantity_id == "S4.cs.1_T2D.A.members"]$n), "T2D FinnGen SuSiE CS at A, p1np 503",
    "", "INTERPRETATION DIFFERENCE: Stage 1 also reads the T2D SuSiE credible set as part of method (1), because the effector rule concerns signal A's T2D association; the effector is 'unresolved' under both readings (logged in results notes)")
cmp("S4.cs.1_protein.B.members; S4.cs.2.B.members", "rs17788484", "4 Seen: (1) and (2) at B", paste(v("S4.cs.1_protein.B.members"), v("S4.cs.2.B.members"), sep = " / "), "Stage 1", "0", "identical")
cmp("S4.cs.2.A.genes_<annot>", "in the ARG1 transcript", "4 Seen: location of the (1)/(2) set at A", sprintf("pc %s; all %s; glist %s", v("S4.cs.2.A.genes_ens116pc"), v("S4.cs.2.A.genes_ens116all"), v("S4.cs.2.A.genes_glist")),
    "genes whose body holds a member", "", "DESCRIPTION DIFFERENCE: both members are in ARG1, and also inside MED23's gene body under all three annotations (Ensembl 116 MED23 6:131,572,244-131,628,313 overlaps ARG1 131,470,832-131,584,332 in GRCh38; glist-hg19 MED23 starts at 131,895,105, before rs2781668 at 131,897,278). The seen wording names ARG1 only")
cmp("S4.cs.3.A.members", "6 variants in MED23, excluding rs2781668", "4 Seen: credible set (3)", sprintf("%s variants in %s, rs2781668 %s", RES[quantity_id == "S4.cs.3.A.members"]$n, v("S4.cs.3.A.genes_ens116pc"),
    if (grepl("rs2781668", v("S4.cs.3.A.members"))) "included" else "excluded"), "Stage 1 p1np 503", "2 variants",
    "gene conclusion and the exclusion of rs2781668 agree; size differs. Explained by the SNP set: the seen script (prespec_review3/cs.R) kept frequency-resolved palindromes, Stage 1 uses p1np (the pre-declaration's '503' label). Post-full-text diagnostic (diagnostics/w4b_diag_cs_after_fulltext.R; no value changed): with frequency-resolved palindromes the same code gives 6 variants; the collinearity cut (0.8 vs 0.9) and the protein prior (sdY 1 vs estimated) change nothing")
cmp("S4.cs.4.A.members", "33 variants over 157 kb (top PIP 0.085)", "4 Seen: credible set (4)", sprintf("%s variants over %.1f kb (top %.3f)", RES[quantity_id == "S4.cs.4.A.members"]$n, span_kb("m4 A 503"), CSR$FULL[["m4 A 503"]]$top_pp), "Stage 1 p1np 503", "8 variants",
    "explained by the SNP set as for (3): with frequency-resolved palindromes the diagnostic gives 33 variants over 157.1 kb, top 0.085 (the seen values exactly)")
cmp("S4.cs.5.A.*; S4.cs.<m>.robust.*", "none", "no seen value", "see results.tsv", "first computed in Stage 1", "", "no seen comparator for method (5) or for any credible-set robustness share")
cmp("S4.effector", "unresolved", "full text 4 Effector: 'On the seen sets this yields unresolved'", v("S4.effector"), "Stage 1 effector rule", "", "agrees; under Stage 1 the sets of methods (1)-T2D and (4) span ARG1, ENPP3 and MED23, and (3)/(5) lie in MED23 only")
# ---- Rule 3(iii) grid (redacted text keeps 'ARG1-arginine 48/48')
cmp("S4.genename.ARG1_arginine.<annot> (unrestricted)", "48/48", "Rule 3(iii) text: SM attribution grid ARG1-arginine 48/48", "48/48 gene-attributed in each annotation file (unrestricted; see notes)", "grid files (definition)", "0",
    "identical (deterministic read of the registered grid). Restricted to settings meeting Rule 3(i) on Stage 1 s: 32/32 in each annotation (first computed in Stage 1)")
cmp("Rule 3(iii) ARG1-ornithine", "no grid exists for the pair", "Rule 3(iii) premise (full and redacted texts)",
    sprintf("the registered grids hold 48 ARG1-ornithine cells per annotation (0/48 gene-attributed; restricted %s/%s, %s/%s, %s/%s)",
            v("S4.genename.ARG1_ornithine.ens116pc.assigned"), v("S4.genename.ARG1_ornithine.ens116pc.total"), v("S4.genename.ARG1_ornithine.ens116all.assigned"),
            v("S4.genename.ARG1_ornithine.ens116all.total"), v("S4.genename.ARG1_ornithine.glist.assigned"), v("S4.genename.ARG1_ornithine.glist.total")), "grid files (definition)", "",
    "PREMISE CONFLICT (the user decides at Stage 2): the premise is false (visible in the redacted text, found before the FREEZE). FIXLOG F1: the branch applies the general gene-naming rule (ARG1 not named, every cell classed 'shared variant, mediation uncertain'); 'attribution not tested' is not used because the attribution was tested; either reading writes no attribution sentence")
# ---- by-product printed in run.log (W4a owns the quantity)
cmp("r(A,B) bootstrap 95% (W4a quantity S4.ld.boot.*)", "-0.107 to +0.029", "4 Seen: LD bootstrap 95% interval", sprintf("%.4f to %.4f", rq[1], rq[2]), "W4b by-product (run.log), R, set.seed(20261007), sample.int draws",
    sprintf("%.4f / %.4f", abs(-0.107 - rq[1]), abs(0.029 - rq[2])), "not a W4b quantity; printed only because the W4b credible-set bootstrap records r(A,B) per draw. The seen interval's generator is not stated (numpy and R scoping scripts both exist); Monte Carlo and generator differences are the likely cause; W4a reports the declared value")
cmp("S4.persig.*", "none comparable", "the Seen block's Rule 3(iv) values are two-signal (A|B) and belong to W4a", "see results.tsv", "Stage 1 per-signal values conditioned on the other four signals", "", "different estimand; no seen comparator")
CS <- rbindlist(rows)
fwrite(CS, file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", quote = TRUE)
logf("COMPARE_SEEN.tsv: %d rows (md5 %s)", nrow(CS), md5(file.path(OUT, "COMPARE_SEEN.tsv")))
logf("REPRODUCIBILITY CHECK (Stage 1 history): results.tsv md5 %s; the FREEZE line recorded 72c729526796514c5d975856989cab4f -> %s (differs by design after Stage 2 step 2, which changed s_max; see the scope check below)", md5(file.path(OUT, "results.tsv")),
     if (md5(file.path(OUT, "results.tsv")) == "72c729526796514c5d975856989cab4f") "IDENTICAL" else "DIFFERS")
# FIX ROUND 2026-10-10: the frozen-copy check that stood here (regex scopes; it only logged ASSERT_FAILED, exit 0, no cell named; a changed
# s_max, a new-fit s or a genename note passed) is replaced by the ID-LIMITED, content-pinned gate w4b_06b_step2_gate.R (next step in
# w4b_run_all.R; table step2_expected_changes.tsv; exit non-zero naming the id and column). Old block: _pre_fix_20261010/w4b_06_compare_seen.R.
logf("==== W4b step 6 done")
