# 07_fisher_d.R -- section 6(d): s per candidate (max over ALL SEVEN distinct fits: the four l3-finemap fits susie, p1fr, p1np on 503 and p2np
# on 489, plus p2np on 503, susie on 489 and p1fr on 489 from 05b_fitset7.R, Stage 2 step 2; the four-fit max is kept as s_max_4fit), the
# Rule 3(i) reclassification (s >= 0.10 for either trait -> not a candidate; s not computable -> reclassified, 'kept' version beside),
# and the pair-level Fisher test / conditional-MLE OR with exact 95% CI, both ways. DPEP1-eGFR run rows (S6.dpep1egfr.*).
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "00_common.R"))
say("==== W6b 07_fisher_d ====")
cand <- fread(file.path(WORK, "candidates_windows.tsv"))
SETS <- c(susie = "503", p1fr = "503", p1np = "503", p2np = "489")
rows <- list(); add <- function(...) rows[[length(rows) + 1]] <<- rrow(...)
PROXY <- "s = susieR::estimate_s_rss on a local SuSiE-RSS fit; a proxy for the pipeline fit (pipeline fit on E:, not local)"
outres <- c(eGFR = "CKDGen 2019 EA eGFR (OUT_CKDGEN2019_eGFR_EA)", CKD = "CKDGen 2019 EA CKD (OUT_CKDGEN2019_CKD_EA)", T2D = "FinnGen R10 T2D (OUT_FG_R10_T2D)",
            CAD = "FinnGen R10 I9_CHD (OUT_FG_R10_I9_CHD)", HF = "FinnGen R10 I9_HEARTFAIL (OUT_FG_R10_I9_HEARTFAIL)")
S <- list()
for (i in seq_len(nrow(cand))) {
  cc <- cand[i]; cid <- cc$cid
  st <- fromJSON(file.path(WORK, "harm", paste0(cid, "__stamp.json")))
  u <- lapply(names(SETS), function(s) { f <- file.path(WORK, "susie", cid, paste0(s, "_", SETS[[s]], ".rds")); if (!file.exists(f)) stop("missing checkpoint ", f); readRDS(f) })
  names(u) <- names(SETS)
  f7 <- readRDS(file.path(WORK, "susie", cid, "fitset7.rds")); stopifnot(identical(f7$cid, cid))
  sv <- rbindlist(lapply(names(u), function(s) { x <- u[[s]]
    if (!is.null(x$note)) return(data.table(set = s, panel = SETS[[s]], trait = c("P", "O"), s_hat = NA_real_, computable = FALSE, nsnp = 0L, conv = NA, n_cs = NA_integer_, n_warn = NA_integer_, H4 = NA_real_))
    tr <- x$traits
    data.table(set = s, panel = SETS[[s]], trait = tr$trait, s_hat = if ("s_hat" %in% names(tr)) tr$s_hat else NA_real_, computable = tr$computable, nsnp = x$nsnp,
               conv = if ("converged" %in% names(tr)) tr$converged else NA, n_cs = if ("n_cs" %in% names(tr)) tr$n_cs else NA_integer_,
               n_warn = if ("n_warn" %in% names(tr)) tr$n_warn else NA_integer_, H4 = x$pair$susie_max_H4) }), fill = TRUE)
  if (any(sv$computable & is.na(sv$s_hat))) stop(cid, ": estimate_s_rss returned NA for a computable trait")
  sv4 <- copy(sv)                                                    # the four-fit set (sensitivity)
  sx <- f7$fits[, .(set, panel, trait, s_hat, computable, nsnp, conv = NA, n_cs = NA_integer_, n_warn = NA_integer_, H4 = NA_real_)]
  sv <- rbind(sv, sx, fill = TRUE)                                   # all seven distinct fits
  stopifnot(nrow(sv) == 14L, !anyDuplicated(sv[, paste(set, panel, trait)]))
  s_of <- function(v, tr) if (any(v[trait == tr]$computable)) max(v[trait == tr & computable == TRUE]$s_hat) else NA_real_
  sP4 <- s_of(sv4, "P"); sO4 <- s_of(sv4, "O"); smax4 <- max(c(sP4, sO4), na.rm = TRUE); at4 <- sv4[computable == TRUE][which.max(s_hat)]
  sP <- if (any(sv[trait == "P"]$computable)) max(sv[trait == "P" & computable == TRUE]$s_hat) else NA_real_
  sO <- if (any(sv[trait == "O"]$computable)) max(sv[trait == "O" & computable == TRUE]$s_hat) else NA_real_
  smax <- max(c(sP, sO), na.rm = TRUE)
  at <- sv[computable == TRUE][which.max(s_hat)]
  P_comp <- any(sv[trait == "P"]$computable)
  reclass <- !P_comp || is.na(sO) || smax >= 0.10
  reclass_kept <- P_comp && !is.na(sO) && smax >= 0.10       # 'the version that keeps it': a non-computable candidate is kept
  why <- if (!P_comp) sprintf("s_protein not computable (protein window %s not local)%s", cc$oid, if (!is.na(sO) && sO >= 0.10) sprintf("; s_outcome %.4f >= 0.10", sO) else "") else
    if (smax >= 0.10) sprintf("s_max %.4f >= 0.10 (%s %s/%s)", smax, at$trait, at$set, at$panel) else sprintf("s_max %.4f < 0.10", smax)
  S[[cid]] <- data.table(cid = cid, gene = cc$gene, outcome = cc$outcome, kidney = cc$kidney, s_protein = sP, s_outcome = sO, s_max = smax, P_comp = P_comp,
                         reclass = reclass, reclass_kept = reclass_kept, why = why, s_protein_4fit = sP4, s_outcome_4fit = sO4, s_max_4fit = smax4)
  pre <- paste0("S6.d.", cid)
  est0 <- sprintf("Rule 3(i) proxy s for %s-%s, max over all seven distinct fits (all variant sets, both panels: susie/p1fr/p1np/p2np on 503 and susie/p1fr/p1np=p2np on 489); %s", cc$gene, cc$outcome, PROXY)
  srcP <- if (P_comp) sprintf("UKB-PPP EUR %s (%s), %s", cc$gene, cc$oid, basename(cc$protein_path)) else "UKB-PPP window not local"
  add(paste0(pre, ".s_protein"), "6(d); Rule 3(i)", est0, srcP, "503+489", value = if (P_comp) sP else "s_protein not computable", n = st$n_P, unit = "s",
      notes = if (P_comp) "" else sprintf("protein window %s_chr%s_cis.tsv absent locally", cc$oid, cc$chr))
  add(paste0(pre, ".s_outcome"), "6(d); Rule 3(i)", est0, outres[[cc$outcome]], "503+489", value = sO, n = st$n_O, unit = "s",
      notes = if (P_comp) st$n_O_source else paste("outcome-only SNP set (outcome x 503 reference by rsID; no protein window): not the pipeline's SNP set;", st$n_O_source))
  add(paste0(pre, ".s_max"), "6(d); Rule 3(i)", est0, paste(srcP, "|", outres[[cc$outcome]]), "503+489", value = if (P_comp) smax else sprintf("not computable (s_outcome %.6g)", sO), unit = "s",
      notes = sprintf("argmax: trait %s, set %s, panel %s; %s", at$trait, at$set, at$panel, why))
  est4 <- sprintf("sensitivity: Rule 3(i) proxy s for %s-%s, max over the earlier four-fit set (susie/p1fr/p1np on 503, p2np on 489); %s", cc$gene, cc$outcome, PROXY)
  add(paste0(pre, ".s_protein_4fit"), "6(d); Rule 3(i) sensitivity", est4, srcP, "503+489", value = if (P_comp) sP4 else "s_protein not computable", n = st$n_P, unit = "s")
  add(paste0(pre, ".s_outcome_4fit"), "6(d); Rule 3(i) sensitivity", est4, outres[[cc$outcome]], "503+489", value = sO4, n = st$n_O, unit = "s")
  add(paste0(pre, ".s_max_4fit"), "6(d); Rule 3(i) sensitivity", est4, paste(srcP, "|", outres[[cc$outcome]]), "503+489", value = if (P_comp) smax4 else sprintf("not computable (s_outcome %.6g)", sO4), unit = "s",
      notes = sprintf("argmax: trait %s, set %s, panel %s; four-fit class_d would be %s", at4$trait, at4$set, at4$panel, if (!P_comp || is.na(sO4) || smax4 >= 0.10) "reclassified" else "kept"))
  add(paste0(pre, ".window"), "6(d)", "SuSiE window (GRCh38) as recorded in LOG_29 (29_cand_egfr.log / 29_cand_nonegfr.log)", "LOG_29", "",
      value = sprintf("chr%s:%d-%d", cc$chr, cc$lo38, cc$hi38), notes = sprintf("GRCh37 reference region %d-%d (%s)", cc$lo37, cc$hi37, cc$bounds37_source))
  add(paste0(pre, ".source"), "6(d)", "inputs of the local proxy fit", "", "", value = paste(srcP, "|", outres[[cc$outcome]], "| LD PLINK --r square on LD503 (g1000_eur) / LD489 (1000G.EUR.QC)"),
      notes = sprintf("SNPs susie %d, p1fr %d, p1np %d, p2np %d; n_P %s; n_O %.0f", st$n_susie, st$n_p1fr, st$n_p1np, st$n_p2np, if (P_comp) sprintf("%.0f", st$n_P) else "NA", st$n_O))
  add(paste0(pre, ".class_d"), "6(d)", "candidate status under the Rule 3(i) sensitivity (s >= 0.10 for either trait, or s not computable -> not a candidate)", "", "",
      value = if (reclass) "reclassified (not a candidate)" else "kept", notes = why)
  for (k in seq_len(nrow(sv))) if (sv$computable[k]) {
    add(sprintf("%s.s.%s.%s.%s", pre, if (sv$trait[k] == "P") "protein" else "outcome", sv$set[k], sv$panel[k]), "6(d); Rule 3(i)", paste("estimate_s_rss per fit;", PROXY),
        if (sv$trait[k] == "P") srcP else outres[[cc$outcome]], sv$panel[k], value = sv$s_hat[k], n = sv$nsnp[k], unit = "s",
        notes = if (is.na(sv$conv[k]) && is.na(sv$n_cs[k])) "s only (no SuSiE fit); added at Stage 2 for the seven-fit s_max (05b_fitset7.R)" else sprintf("SuSiE L=10 converged %s, n_cs %s, warnings %s", sv$conv[k], sv$n_cs[k], sv$n_warn[k]))
  }
  for (s in names(u)) if (is.null(u[[s]]$note) && P_comp) {
    h4 <- u[[s]]$pair$susie_max_H4; trs <- u[[s]]$traits
    nr <- if (is.na(h4)) sprintf("NOT RUN: coloc.susie needs >= 1 credible set per trait (n_cs protein %s, outcome %s)", trs[trait == "P"]$n_cs, trs[trait == "O"]$n_cs) else ""
    add(sprintf("%s.coloc_susie_H4.%s.%s", pre, s, SETS[[s]]), "6(d) supplement (local fit; Rule 3(i) status per s)",
      "coloc.susie max PP.H4 over credible-set pairs, p12 = 5e-6, on the local proxy fit", paste(srcP, "|", outres[[cc$outcome]]), SETS[[s]],
      value = if (is.na(h4)) "not computed" else h4, unit = "PP.H4",
      notes = if (nzchar(nr)) nr else sprintf("PP.H3 at max %.4g; p12 1e-5 max %.4g; cs pairs %d", u[[s]]$pair$susie_H3_at_max, u[[s]]$pair$susie_max_H4_p12_1e5, u[[s]]$pair$n_cs_pairs))
  }
  # DPEP1-eGFR run deliverable
  if (cid == "DPEP1_eGFR") for (s in names(u)) {
    tr <- u[[s]]$traits; pp <- sprintf("S6.dpep1egfr.%s.%s", s, SETS[[s]])
    e1 <- "DPEP1-eGFR l3-finemap-style run (UKB-PPP OID20897 x CKDGen 2019 EA eGFR, own DPEP1-eGFR harmonised set); SuSiE-RSS L = 10"
    add(paste0(pp, ".s"), "6 'also built'; 5 (W5 reads eGFR s)", paste(e1, "- eGFR estimate_s_rss"), outres[["eGFR"]], SETS[[s]], value = tr[trait == "O"]$s_hat, n = u[[s]]$nsnp, unit = "s")
    add(paste0(pp, ".s_protein"), "6 'also built'", paste(e1, "- DPEP1 estimate_s_rss"), "UKB-PPP OID20897", SETS[[s]], value = tr[trait == "P"]$s_hat, n = u[[s]]$nsnp, unit = "s")
    add(paste0(pp, ".cs"), "6 'also built'", paste(e1, "- eGFR credible sets (purity 0.5)"), outres[["eGFR"]], SETS[[s]], value = tr[trait == "O"]$n_cs, unit = "count",
        notes = sprintf("converged %s, niter %s, warnings %s", tr[trait == "O"]$converged, tr[trait == "O"]$niter, tr[trait == "O"]$n_warn))
    add(paste0(pp, ".cs_protein"), "6 'also built'", paste(e1, "- DPEP1 credible sets (purity 0.5)"), "UKB-PPP OID20897", SETS[[s]], value = tr[trait == "P"]$n_cs, unit = "count",
        notes = sprintf("converged %s, niter %s, warnings %s", tr[trait == "P"]$converged, tr[trait == "P"]$niter, tr[trait == "P"]$n_warn))
    add(paste0(pp, ".coloc_H4"), "6 'also built'", paste(e1, "- coloc.susie DPEP1 x eGFR, max PP.H4 over credible-set pairs, p12 = 5e-6"), "UKB-PPP OID20897 x CKDGen eGFR", SETS[[s]],
        value = u[[s]]$pair$susie_max_H4, unit = "PP.H4", notes = sprintf("PP.H3 at max %.4g; p12 = 1e-5 max PP.H4 %.4g; Rule 3(i) admissibility per S6.d.DPEP1_eGFR.s_max",
                                                                          u[[s]]$pair$susie_H3_at_max, u[[s]]$pair$susie_max_H4_p12_1e5))
  }
}
S <- rbindlist(S)
print(S[, .(cid, kidney, s_protein = signif(s_protein, 4), s_outcome = signif(s_outcome, 4), s_max = signif(s_max, 4), reclass, reclass_kept)])
fwrite(S, file.path(WORK, "d_classification.tsv"), sep = "\t")

base <- readRDS(file.path(WORK, "primary_table_from_logs.rds"))
# FIXER F1 (FIXLOG.md in this folder and in W6a): the conditional-MLE OR and its exact CI come from fisher.test's own equations
# solved at tol 1e-13 (one shared implementation, W6a/exact_cond_or.R, md5 asserted); fisher.test's default uniroot tolerance
# moved the upper bounds in the 3rd decimal. The p value is still fisher.test's (an exact sum; no root finding).
ECOR <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_06a_renal/exact_cond_or.R")
assert_md5(ECOR, "1217326d2a6880b174f91d377c5b21b9", "SHARED:exact_cond_or.R (W6a folder; FIXER F1)")
source(ECOR)
# Stage 2 step 2: with seven fits the (d) table is 2/0/66/92 (OR infinite), where Inf - Inf is NaN; an infinite bound equal in sign on both sides
# counts as a relative difference of 0 (any sign or finiteness mismatch stays NaN and fails the check).
relmax <- function(a, b) { r <- ifelse(is.infinite(a) & is.infinite(b) & sign(a) == sign(b), 0, abs(a - b) / abs(b)); max(r) }
fisher_row <- function(drop_k, drop_nk, tag, label) {
  a <- base$a - drop_k; b <- base$b - drop_nk; c_ <- base$c + drop_k; d <- base$d + drop_nk
  f0_ <- fisher.test(matrix(c(a, c_, b, d), nrow = 2))
  e <- exact_cond_or(matrix(c(a, c_, b, d), nrow = 2), conf.level = 0.95, tol = 1e-13)
  ftv <- c(unname(f0_$estimate), f0_$conf.int)
  ft <- list(p.value = f0_$p.value, estimate = e$estimate, conf.int = e$conf.int)
  pq <- if (tag == "") "S6.d" else paste0("S6.d.", tag)
  log_assert(paste0(pq, ".exact_roots"), "<1e-10", sprintf("%.2g", max(abs(e$resid))), max(abs(e$resid)) < 1e-10,
             "exact_cond_or: max |residual| of E[X]=x, P(X>=x)=0.025, P(X<=x)=0.025 at the reported MLE and bounds (FIXLOG F1)")
  log_assert(paste0(pq, ".matches_fishertest"), "<0.01", sprintf("%.2g", relmax(c(e$estimate, e$conf.int), ftv)), relmax(c(e$estimate, e$conf.int), ftv) < 0.01,
             "same estimand as fisher.test: MLE and bounds within 1% of fisher.test's default-tolerance values (FIXLOG F1)")
  ftn <- function(v) sprintf("exact CI: fisher.test's equations solved at tol 1e-13 (exact_cond_or.R); R fisher.test (uniroot at its default tolerance) gave %.10g; superseded (FIXLOG F1)", v)
  est <- sprintf("%s: pair-level Fisher exact test and conditional-MLE OR (exact 95%% CI), rows candidate yes/no x columns kidney yes/no; a = kidney candidates, b = other candidates, c = kidney non-candidates, d = other non-candidates", label)
  pp <- if (tag == "") "S6.d" else paste0("S6.d.", tag)
  add(paste0(pp, ".a"), "6(d)", est, "160-pair pool (29-log headers) x section 6(d) classes", value = a, unit = "count")
  add(paste0(pp, ".b"), "6(d)", est, "", value = b, unit = "count")
  add(paste0(pp, ".c"), "6(d)", est, "", value = c_, unit = "count")
  add(paste0(pp, ".d"), "6(d)", est, "", value = d, unit = "count")
  add(paste0(pp, ".p"), "6(d)", est, "", value = ft$p.value, p = ft$p.value, unit = "p (two-sided)")
  add(paste0(pp, ".or"), "6(d)", est, "", value = unname(ft$estimate), ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2], p = ft$p.value, unit = "OR (conditional MLE)", notes = ftn(ftv[1]))
  add(paste0(pp, ".lo"), "6(d)", est, "", value = ft$conf.int[1], unit = "OR lower 95% exact", notes = ftn(ftv[2]))
  add(paste0(pp, ".hi"), "6(d)", est, "", value = ft$conf.int[2], unit = "OR upper 95% exact", notes = ftn(ftv[3]))
  say("%s: %d/%d vs %d/%d -> p %.4g OR %.4f (%.4f-%.4f)", label, a, a + b, c_, c_ + d, ft$p.value, ft$estimate, ft$conf.int[1], ft$conf.int[2])
  list(a = a, b = b, c = c_, d = d, p = ft$p.value, or = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2])
}
fd <- fisher_row(S[reclass == TRUE & kidney == TRUE, .N], S[reclass == TRUE & kidney == FALSE, .N], "", "6(d) Rule 3(i) sensitivity (non-computable candidates reclassified)")
fk <- fisher_row(S[reclass_kept == TRUE & kidney == TRUE, .N], S[reclass_kept == TRUE & kidney == FALSE, .N], "kept", "6(d) 'kept' version (the non-computable candidate kept; computable candidates as in (d))")
f0 <- fisher_row(0, 0, "check_primary", "check: W6b recomputation of the primary table from the 29-log pool (W6a holds the primary)")
add("S6.d.reclassified", "6(d)", "candidates reclassified in the (d) estimate", "", value = paste(S[reclass == TRUE]$cid, collapse = ";"), n = S[reclass == TRUE, .N], unit = "list")
add("S6.d.kept.reclassified", "6(d)", "candidates reclassified in the 'kept' version", "", value = paste(S[reclass_kept == TRUE]$cid, collapse = ";"), n = S[reclass_kept == TRUE, .N], unit = "list")
add("S6.d.rule3i_scope", "Rule 3(i); 6(d)", "first-mention flag scope", "", value = "s not computed (pipeline fit on E:)",
    notes = "DEVIATION: Rule 3(i) first-mention flag computed only for the slot loci and the 16 candidates; the declared rule covers every pair whose class rests on reference-LD SuSiE; user decision pending at Stage 2")
for (sl in c("UMOD_eGFR", "UMOD_CKD", "DPEP1_eGFR")) add(sprintf("S6.slot_s.%s", sl), "Rule 3(i) (slot loci touched by W6b)", "s_max over both traits and all seven distinct fits (see S6.d.<cid>.s_*; four-fit sensitivity in S6.d.<cid>.s_max_4fit)",
    "", "503+489", value = S[cid == sl]$s_max, unit = "s", notes = sprintf("s_protein %.6g; s_outcome %.6g; admissible iff both < 0.10; four-fit sensitivity s_max %.6g", S[cid == sl]$s_protein, S[cid == sl]$s_outcome, S[cid == sl]$s_max_4fit))
saveRDS(list(S = S, fd = fd, fk = fk, f0 = f0), file.path(WORK, "d_summary.rds"))
write_rows(rbindlist(rows), "07_fisher_d")
say("07_fisher_d done")
