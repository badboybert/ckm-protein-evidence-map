# w4b_05_rules.R -- Stage 1 W4b, step 5: s results and Rule 3(i) admissibility, Rule 3(iii) gene naming (restricted grid),
# the effector rule, the omnibus wording branches; assembles results.tsv and branches.tsv from the step parts.
# Uses only the redacted pre-declaration's rules (the full text is opened only afterwards, for COMPARE_SEEN.tsv).
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
logf("==== W4b step 5 (rules, assembly) start")
D <- readRDS(file.path(OUT, "data", "w4b_data.rds")); H <- D$H
CSR <- readRDS(file.path(OUT, "data", "w4b_credsets.rds")); OMR <- readRDS(file.path(OUT, "data", "w4b_omnibus.rds"))
acc <- new_rows(); BR <- list()
addbr <- function(rule_ref, inputs, branch, filled_template, qualifiers = "") {
  BR[[length(BR) + 1]] <<- data.table(rule_ref = rule_ref, inputs = inputs, branch = branch, filled_template = filled_template, qualifiers = qualifiers)
}
f3 <- function(x) formatC(x, digits = 3, format = "f")
fp <- function(p) if (p < 0.001) formatC(p, digits = 1, format = "e") else formatC(signif(p, 2), digits = 2, format = "fg", flag = "#")

# ---------------------------------------------------------------- s (Rule 3(i))
S <- fread(file.path(OUT, "susie", "susie_s_table.tsv"))
stopifnot(nrow(S) == 16, !anyNA(S$s))
# Stage 2 step 2 (user decision 2026-10-09): s_max over ALL SEVEN distinct fits per trait (all variant sets, both panels); the
# four-fit value (susie, p1fr, p1np on 503; p2np on 489) is kept as a labelled sensitivity. The 12 extra fits are s-only
# (susie/susie_s_extra.tsv from w4b_02b_fitset7.R): p2np on 503, susie on 489, p1fr on 489.
SX <- fread(file.path(OUT, "susie", "susie_s_extra.tsv"))
stopifnot(nrow(SX) == 12, !anyNA(SX$s))
S4 <- copy(S)
S <- rbind(S, SX[, names(S4)[names(S4) %in% names(SX)], with = FALSE], fill = TRUE)
stopifnot(nrow(S) == 28, !anyDuplicated(S[, paste(trait, set, panel)]))
TLAB <- c(P__ARG1 = "UKB-PPP OID21445 (ARG1)", arginine = "CLSA GCST90200372", ornithine = "CLSA GCST90200404", T2D_FG = "FinnGen R10 T2D")
for (i in seq_len(nrow(S))) acc <- add_row(acc, sprintf("S4.s.%s.%s.%s", S$trait[i], S$set[i], S$panel[i]), "4(v)(1); Rule 3(i)",
   sprintf("susieR::estimate_s_rss (null-mle) for %s on SNP set %s, panel %s (PLINK r)", S$trait[i], S$set[i], S$panel[i]), TLAB[[S$trait[i]]],
   S$panel[i], value = S$s[i], n = S$n[i], unit = "s",
   notes = if (i <= 16) sprintf("%d SNPs; SuSiE L=10 converged %s, %d CS, evaluable %s; n = %s; DEVIATION: engine path-locked (QC and SuSiE settings re-implemented in R)",
                   S$nsnp[i], S$converged[i], S$n_cs[i], S$evaluable[i], if (S$trait[i] == "T2D_FG") "N_eff from case/control counts" else "n_hat (engine rule)") else
           sprintf("%d SNPs; s only (no SuSiE fit); added at Stage 2 for the seven-fit s_max (%s on panel %s); n = %s", S$nsnp[i], S$set[i], S$panel[i], if (S$trait[i] == "T2D_FG") "N_eff from case/control counts" else "n_hat (engine rule)"))
PAIRS <- list(arginine = c("P__ARG1", "arginine"), ornithine = c("P__ARG1", "ornithine"), T2D = c("P__ARG1", "T2D_FG"))
SMAX <- list(); SARG <- list()
for (pr in names(PAIRS)) {
  sp <- S[trait %in% PAIRS[[pr]]]; mx <- max(sp$s); w <- sp[which.max(sp$s)]
  sp4 <- S4[trait %in% PAIRS[[pr]]]; mx4 <- max(sp4$s); w4 <- sp4[which.max(sp4$s)]
  SMAX[[pr]] <- mx; SARG[[pr]] <- sprintf("%s on the %s panel", w$set, w$panel)
  adm <- mx < 0.10
  acc <- add_row(acc, sprintf("S4.s.max.ARG1_%s", pr), "Rule 3(i)", sprintf("max s over every fit for the pair ARG1-%s (both traits; all seven distinct fits: all variant sets, both panels)", pr),
                 paste(TLAB[PAIRS[[pr]]], collapse = " + "), "503+489", value = mx, unit = "s",
                 notes = sprintf("max at %s/%s/%s; per trait max: %s; admissible (s < 0.10): %s", w$trait, w$set, w$panel,
                                 paste(sprintf("%s %.4f", PAIRS[[pr]], vapply(PAIRS[[pr]], function(t) max(S[trait == t]$s), 0)), collapse = ", "), adm))
  acc <- add_row(acc, sprintf("S4.s.max_4fit.ARG1_%s", pr), "Rule 3(i) sensitivity", sprintf("sensitivity: max s for the pair ARG1-%s over the earlier four-fit set (susie, p1fr, p1np on 503; p2np on 489)", pr),
                 paste(TLAB[PAIRS[[pr]]], collapse = " + "), "503+489", value = mx4, unit = "s",
                 notes = sprintf("max at %s/%s/%s; admissible (s < 0.10): %s", w4$trait, w4$set, w4$panel, mx4 < 0.10))
  addbr("Rule 3(i)", sprintf("S4.s.max.ARG1_%s; S4.s.max_4fit.ARG1_%s; S4.s.<trait>.<set>.<panel> (7 distinct fits per trait)", pr, pr), if (adm) "admissible" else "inadmissible",
        sprintf("ARG1-%s: SuSiE with reference LD (PP.H4, credible sets) is %s (s_max %s, %s/%s/%s)%s", pr, if (adm) "admissible" else "inadmissible",
                f3(mx), w$trait, w$set, w$panel, if (adm) "" else "; reported in the supplement, labelled, never cited as support in new text"),
        sprintf("locally computed s for a pipeline fit is a proxy and is labelled as one; four-fit sensitivity s_max %s (%s/%s/%s)", f3(mx4), w4$trait, w4$set, w4$panel))
}

# ---------------------------------------------------------------- Rule 3(iii) gene naming (restricted grid)
GRIDS <- list(ens116pc = rpath("GRID_bridge_class_robustness__ens116pc"), ens116all = rpath("GRID_bridge_class_robustness__ens116all"),
              glist = rpath("GRID_bridge_class_robustness"))
GA <- "gene-attributed bridge"
NAMED <- list()
# FIXLOG F1: the Rule 3(iii) clauses for one pair, evaluated on that pair's grid in the declared order (used for ARG1-ornithine):
#   no grid rows for the pair -> no attribution sentence, the SM says 'attribution not tested';
#   restricted grid empty -> unresolved;
#   ARG1 assigned in a majority of restricted settings under each annotation -> named;
#   otherwise -> ARG1 not named, no attribution sentence.
# n_rows, assigned, total, ga_all: integer vectors named by annotation; cls: description of the classes over the pair's cells.
genename_branch <- function(pair, n_rows, assigned, total, ga_all, cls) {
  cnt <- paste(sprintf("%s %d/%d", names(total), assigned, total), collapse = ", ")
  if (all(n_rows == 0)) return(list(branch = "no grid for the pair", text = "no attribution sentence is written; the SM says 'attribution not tested'"))
  if (any(total == 0)) return(list(branch = "unresolved", text = sprintf("the attribution of the %s association to a gene is unresolved (restricted grid empty)", pair)))
  if (all(assigned > total / 2)) return(list(branch = "named", text = sprintf("ARG1 is named for the %s association: the restricted grid assigns it in a majority of settings under each of the three annotations (restricted %s)", pair, cnt)))
  list(branch = "not named (grid exists, no majority)",
       text = sprintf("no attribution sentence is written; the SM prints the existing ARG1-%s grid under its own standard: ARG1 assigned in restricted settings %s (unrestricted %s); %s",
                      pair, cnt, paste(sprintf("%s %d/%d", names(n_rows), ga_all, n_rows), collapse = ", "), cls))
}
GN <- list(); GCLS <- list()
for (met in c("arginine", "ornithine")) {
  smx <- SMAX[[met]]; NAMED[[met]] <- list(); GN[[met]] <- list(); GCLS[[met]] <- character(0)
  for (a in names(GRIDS)) {
    g <- fread(GRIDS[[a]])[protein == "ARG1" & metabolite == met]
    stopifnot(nrow(g) == 48L, !anyDuplicated(g$variant))
    # a setting meets Rule 3(i) unless its class rests on SuSiE-with-reference-LD evidence that is inadmissible on our own s:
    # E1 (the SuSiE colocalization flag) TRUE in that cell while s_max >= 0.10
    e1 <- toupper(as.character(g$E1)); stopifnot(all(e1 %in% c("TRUE", "FALSE")))      # fread may read True/False as logical
    meets <- if (smx < 0.10) rep(TRUE, nrow(g)) else e1 != "TRUE"
    stopifnot(sum(meets) == if (smx < 0.10) nrow(g) else sum(e1 == "FALSE"))
    alt_axis <- if (smx < 0.10) rep(TRUE, nrow(g)) else g$susie %in% c("none", "sguard")
    ass <- sum(meets & g$class == GA); tot <- sum(meets)
    NAMED[[met]][[a]] <- tot > 0 && ass > tot / 2
    GN[[met]][[a]] <- c(n_rows = nrow(g), assigned = ass, total = tot, ga_all = sum(g$class == GA))
    GCLS[[met]] <- union(GCLS[[met]], unique(g$class))
    nt <- sprintf("restricted = cells meeting Rule 3(i) on our own s (s_max %.4f): %s; classes in the restricted cells: %s; unrestricted gene-attributed %d/48; alternative restriction by the susie axis only (none, sguard): %d/%d",
                  smx, if (smx < 0.10) "all 48 (SuSiE admissible)" else "cells whose E1 (SuSiE coloc flag) is FALSE",
                  paste(sprintf("%s=%d", names(table(g$class[meets])), table(g$class[meets])), collapse = "; "), sum(g$class == GA),
                  sum(alt_axis & g$class == GA), sum(alt_axis))
    if (met == "ornithine") nt <- paste0(nt, "; PREMISE CONFLICT: the pre-declaration and brief state no grid exists for ARG1-ornithine, but the grid file holds these 48 cells; the Rule 3(iii) branch applies the general gene-naming rule to them (FIXLOG F1)")
    acc <- add_row(acc, sprintf("S4.genename.ARG1_%s.%s.assigned", met, a), "Rule 3(iii)", sprintf("restricted grid cells assigning ARG1 (class '%s'), annotation %s", GA, a),
                   basename(GRIDS[[a]]), "grid", value = ass, n = tot, unit = "cells", notes = nt)
    acc <- add_row(acc, sprintf("S4.genename.ARG1_%s.%s.total", met, a), "Rule 3(iii)", sprintf("restricted grid cells (meeting Rule 3(i)), annotation %s", a),
                   basename(GRIDS[[a]]), "grid", value = tot, unit = "cells", notes = nt)
    logf("GENENAME ARG1-%s %s: assigned %d / restricted %d (unrestricted GA %d/48; axis-only restriction %d/%d)", met, a, ass, tot, sum(g$class == GA),
         sum(alt_axis & g$class == GA), sum(alt_axis))
  }
}
ga_in <- function(met) paste(sprintf("S4.genename.ARG1_%s.%s.{assigned,total}", met, names(GRIDS)), collapse = "; ")
named_arg <- all(unlist(NAMED$arginine))
addbr("Rule 3(iii) gene naming (ARG1-arginine)", paste0(ga_in("arginine"), "; S4.s.max.ARG1_arginine"), if (named_arg) "named" else "unresolved",
      if (named_arg) "ARG1 is named for the arginine association: the restricted grid assigns it in a majority of settings under each of the three annotations (counts in the SM)"
      else "the attribution of the arginine association to a gene is unresolved (restricted grid empty or no majority under every annotation)",
      sprintf("restricted counts (%s) printed under the grid's own standard in the SM beside the unrestricted 48-cell counts; the grid's A2 test reads 'ARG1 in the gene list', and signal A (rs2781668) also lies in MED23's gene body under all three annotations (signal B under the two Ensembl annotations); ratio sentences per SD of the measured protein need no grid; the registered class does not rest on the inadmissible SuSiE evidence (the same class in every restricted cell)",
              paste(sprintf("%s %s", names(GRIDS), vapply(names(GRIDS), function(a) { x <- acc[[which(vapply(acc, function(r) r$quantity_id == sprintf("S4.genename.ARG1_arginine.%s.assigned", a), TRUE))]]; paste0(x$value, "/", x$n) }, "")), collapse = ", ")))
gv <- function(met, k) vapply(names(GRIDS), function(a) as.integer(GN[[met]][[a]][[k]]), 0L)
orn_cls <- if (length(GCLS$ornithine) == 1) sprintf("every cell classed '%s' under each annotation", GCLS$ornithine) else sprintf("classes over the cells: %s", paste(GCLS$ornithine, collapse = "; "))
ob <- genename_branch("ornithine", gv("ornithine", "n_rows"), gv("ornithine", "assigned"), gv("ornithine", "total"), gv("ornithine", "ga_all"), orn_cls)
logf("GENENAME ARG1-ornithine branch (FIXLOG F1): %s", ob$branch)
addbr("Rule 3(iii) gene naming (ARG1-ornithine)", paste0(ga_in("ornithine"), "; S4.s.max.ARG1_ornithine"), ob$branch, ob$text,
      sprintf("PREMISE CONFLICT, the user decides at Stage 2 (FIXLOG F1): Rule 3(iii) names ARG1-ornithine as the pair with no grid and gives it the SM phrase 'attribution not tested', but the registered grid files hold %s ARG1-ornithine cells under each annotation, so the attribution was tested and that phrase would misstate it; the clause's condition (no grid for the pair) is evaluated on the registered files and is false, so the general gene-naming rule of Rule 3(iii) is applied; no estimand, threshold or rule is changed; under either reading no attribution sentence is written; the ornithine wording rule (ornithine per SD of ARG1) attributes nothing and needs no grid",
              paste(unique(gv("ornithine", "n_rows")), collapse = "/")))

# ---------------------------------------------------------------- effector rule (signal A's T2D association)
TAB <- CSR$TAB[signal == "A"]
ANN <- c("ens116pc", "ens116all", "glist")
adm_m <- c("1_T2D" = SMAX$T2D < 0.10, "1_protein" = SMAX$T2D < 0.10, "2" = TRUE, "3" = TRUE, "4" = TRUE, "5" = TRUE)
stopifnot(setequal(TAB$method, names(adm_m)))
EFF <- list()
for (a in ANN) {
  w <- setNames(strsplit(TAB[[paste0("within_", a)]], ";"), TAB$method)
  cand <- Reduce(intersect, w)
  cand <- cand[vapply(cand, function(g) any(vapply(names(w), function(m) adm_m[[m]] && g %in% w[[m]], TRUE)), TRUE)]
  # several genes meeting the criteria (overlapping gene bodies) do not identify one effector -> unresolved, with the names
  EFF[[a]] <- if (length(cand) == 1) cand else if (length(cand) > 1) sprintf("unresolved (criteria met by %s: overlapping gene bodies)", paste(cand, collapse = " and ")) else "unresolved"
  # readings of method (1) (protein-only / T2D-only) printed beside
  alt <- vapply(c("1_T2D", "1_protein"), function(drop) { ww <- w[setdiff(names(w), drop)]; cc <- Reduce(intersect, ww); if (length(cc)) paste(cc, collapse = "+") else "unresolved" }, "")
  acc <- add_row(acc, sprintf("S4.effector.%s", a), "4 Effector; Rule 3(iii) annotations", "gene placed by EVERY credible-set method (1)-(5) at signal A (all members inside it), with >= 1 admissible method",
                 "credible_sets_by_method.tsv", "503 (methods 1,3-5); in-sample (2)", value = EFF[[a]], unit = "gene",
                 notes = sprintf("within per method: %s; admissible: %s; reading (1) as protein-only: %s; as T2D-only: %s",
                                 paste(sprintf("%s=%s", names(w), vapply(w, function(x) if (length(x)) paste(x, collapse = "+") else "-", "")), collapse = ", "),
                                 paste(names(adm_m)[adm_m], collapse = ","), alt[["1_T2D"]], alt[["1_protein"]]))
  logf("EFFECTOR %s: %s (within: %s)", a, EFF[[a]], paste(sprintf("%s=%s", names(w), vapply(w, function(x) paste(x, collapse = "+"), "")), collapse = ", "))
}
eff_all <- if (length(unique(unlist(EFF))) == 1 && !startsWith(EFF[[1]], "unresolved")) EFF[[1]] else "unresolved"
acc <- add_row(acc, "S4.effector", "4 Effector", "effector of signal A's T2D association (same gene under all three annotations)", "", "", value = eff_all, unit = "gene")
addbr("4 Effector", "S4.effector.{ens116pc,ens116all,glist}; S4.cs.<method>.A.{members,genes_<annot>}; S4.s.max.ARG1_T2D", if (eff_all == "unresolved") "unresolved" else "named",
      if (eff_all == "unresolved") "the attribution of signal A's T2D association to a gene is unresolved"
      else sprintf("%s is named as the effector of signal A's T2D association: every credible-set method places the set within %s under each annotation, and at least one admissible method does", eff_all, eff_all),
      sprintf("every 'unresolved' gene attribution travels (Principle 0); method admissibility and robustness are reported per method; no gene holds every method's set at signal A (ens116pc within-gene per method: %s)",
              paste(sprintf("(%s) %s", TAB$method, ifelse(nzchar(TAB$within_ens116pc), TAB$within_ens116pc, "none")), collapse = "; ")))

# credible-set methods: admissibility / robustness status rows
for (k in names(CSR$ROB)) {
  m <- sub(" .*", "", k); X <- sub(".* ", "", k); r <- CSR$ROB[[k]]
  idb <- if (X == "A") sprintf("S4.cs.%s.robust", sub("m", "", m)) else sprintf("S4.cs.%s.B.robust", sub("m", "", m))
  addbr("4(v) method status; Rule 3(ii)", sprintf("%s.panel489_same; %s.boot_share_same_genes; S4.cs.%s.%s.genes_<annot>", idb, idb, sub("m", "", m), X),
        if (r$robust) "robust" else "not robust",
        sprintf("method (%s) credible set at signal %s: %s", sub("m", "", m), X, if (r$robust) "the gene(s) it spans are the same on 489 and in >= 90% of 5,000 503 draws" else "sensitive to the reference LD"),
        if (r$robust) "" else "add 'sensitive to the reference LD' wherever this set is cited (Principle 0)")
}
addbr("4(v) method status; Rule 3(i)", "S4.s.max.ARG1_T2D; S4.cs.1_T2D.A.members; S4.cs.1_protein.A.members", if (SMAX$T2D < 0.10) "admissible" else "inadmissible",
      sprintf("method (1) SuSiE with the 503 reference (credible sets, p1np): %s (s_max %s, the maximum over all seven fits on the 503 and 489 panels; attained by %s)", if (SMAX$T2D < 0.10) "admissible" else "inadmissible; supplement only, labelled, never cited as support", f3(SMAX$T2D), SARG$T2D), "")
addbr("4(v) method status", "S4.cs.2.A.members", "admissible", "method (2) UKB-PPP in-sample SuSiE: admissible (in-sample LD)", "")

# ---------------------------------------------------------------- omnibus wording (secondary; outside the Rule 4 family)
R0 <- rbindlist(lapply(list.files(file.path(OUT, "parts"), "^results_0[34].*\\.tsv$", full.names = TRUE), fread, colClasses = "character"))
getv <- function(id) { v <- R0[quantity_id == id]$value; stopifnot(length(v) == 1); as.numeric(v) }
omni_sentence <- function(tr, lab) {
  Q <- getv(sprintf("S4.omni.%s.5sig.aware.Q", tr)); df <- getv(sprintf("S4.omni.%s.5sig.aware.df", tr))
  p <- getv(sprintf("S4.omni.%s.5sig.p_aware_2nd", tr)); pl <- getv(sprintf("S4.omni.%s.5sig.aware.p_lower", tr))
  core <- if (p < 0.05) sprintf("%s: heterogeneous across the five signals (Q %s, df %d, p %s)", lab, f3(Q), df, fp(p))
          else sprintf("%s: no detectable heterogeneity across the five signals (Q %s, df %d, p %s)", lab, f3(Q), df, fp(p))
  if (pl < 0.05) core <- paste0(core, sprintf("; Q was below its expected value (lower-tail p %s)", fp(pl)))
  list(text = core, branch = if (p < 0.05) "p < 0.05" else "p >= 0.05", lower = pl < 0.05,
       beside = sprintf("ign p %s; 1st-order p %s; COJO NOT RUN for GLS Q", fp(getv(sprintf("S4.omni.%s.5sig.p_ign_2nd", tr))), fp(getv(sprintf("S4.omni.%s.5sig.p_aware_1st", tr)))))
}
oa <- omni_sentence("arginine", "arginine per SD of ARG1"); ot <- omni_sentence("T2D_FG", "T2D (FinnGen R10) per SD of ARG1"); od <- omni_sentence("T2D_DM", "T2D (DIAMANTE, beside)")
ids <- function(tr) sprintf("S4.omni.%s.5sig.{aware.Q,aware.df,p_aware_2nd,aware.p_lower,p_ign_2nd,p_aware_1st}", tr)
q_omni <- "secondary: outside the Rule 4 family and never in the Abstract; the text states that the omnibus and the A-vs-B test answer different questions; no calibration or power statement (section 12); per-signal z and PP.H4 in the supplementary table (Rule 3(iv)), no clause added; Rule 4b wording; Rule 3(ii) is not evaluated (the ARG1 omnibus is not a named primary cell)"
addbr("4 Omnibus (arginine)", ids("arginine"), oa$branch, oa$text, paste0(q_omni, "; beside: ", oa$beside))
addbr("4 Omnibus (T2D)", paste(ids("T2D_FG"), ids("T2D_DM"), sep = "; "), ot$branch,
      paste0(ot$text, "; DIAMANTE beside: ", sub("^T2D \\(DIAMANTE, beside\\): ", "", od$text)), paste0(q_omni, "; beside: FinnGen ", ot$beside, "; DIAMANTE ", od$beside))
addbr("4 Omnibus (six-signal sensitivity)", "S4.omni.arginine.6sig.*; S4.omni.T2D_DM.6sig.*; S4.omni.T2D_FG.6sigproxy.*; S4.omni.T2D_FG.6sigdrop5.*", "no wording rule (declared sensitivity, printed beside the omnibus)",
      sprintf("six-signal UKB-PPP in-sample set: arginine p %s; T2D FinnGen (proxy for rs2698107) p %s, (five available signals) p %s; T2D DIAMANTE p %s",
              fp(getv("S4.omni.arginine.6sig.p_aware_2nd")), fp(getv("S4.omni.T2D_FG.6sigproxy.p_aware_2nd")), fp(getv("S4.omni.T2D_FG.6sigdrop5.p_aware_2nd")),
              fp(getv("S4.omni.T2D_DM.6sig.p_aware_2nd"))),
      "post hoc set choice disclosed: the five-signal set was chosen after seeing the five-signal Q values; DEVIATION: rs2698107 absent from FinnGen R10")
write_part(acc, "05_rules")

# ---------------------------------------------------------------- assemble results.tsv and branches.tsv
RES <- rbindlist(lapply(sort(list.files(file.path(OUT, "parts"), "^results_.*\\.tsv$", full.names = TRUE)), fread, colClasses = "character"))
stopifnot(!anyDuplicated(RES$quantity_id), all(names(RES) == RCOLS))
# no NA where a value is expected
bad <- RES[value == "" | is.na(value)]
if (nrow(bad)) { logf("EMPTY VALUES in %d rows: %s", nrow(bad), paste(head(bad$quantity_id, 10), collapse = ", ")); stop("empty values in results") }
fwrite(RES, file.path(OUT, "results.tsv"), sep = "\t", quote = TRUE)
BRT <- rbindlist(BR); fwrite(BRT, file.path(OUT, "branches.tsv"), sep = "\t", quote = TRUE)
logf("results.tsv: %d rows (md5 %s); branches.tsv: %d rows (md5 %s)", nrow(RES), md5(file.path(OUT, "results.tsv")), nrow(BRT), md5(file.path(OUT, "branches.tsv")))
logf("==== W4b step 5 done")
