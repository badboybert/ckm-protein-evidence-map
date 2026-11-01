# 43b (2026-06-30, gate fix #2): BEFORE Phase 5 (44), assert the Phase-4 annotation master covers
# EVERY corrected strong gene. phase4/12 + 09_merge_master.py only iterate the cached gene universe
# (cache/phase4_gene_ids.csv, the prior 47-gene set), so any NEW strong gene produced by the
# UKB-free rerun would be SILENTLY un-annotated. This fails loud instead. If it aborts, annotate the
# listed gene(s) (add to gene_ids + run the 42* annotation scripts + re-merge), then re-run this.
# Run: Rscript scripts/paper3.1_ckm/cascade/16_assert_annotation.R
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
suppressPackageStartupMessages(library(data.table))
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_guard.R"))
RC <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
v6f <- file.path(RC, "coloc_results_v7.rds")
if (!file.exists(v6f)) rerun_abort(sprintf("43b: %s missing — run Phase 3 (05/06/08) first", v6f))
v6 <- as.data.table(readRDS(v6f))
strong_genes <- sort(unique(v6[is_coloc_pass(PP.H4, PP.H3), gene]))   # dual coloc rule (D1 2026-07-04)
guard_nonempty(data.frame(g = strong_genes), "corrected strong genes (PP.H4>=0.80)", min_rows = 1L)
# 15_stage_annotation.R points this at a CANDIDATE file so the master is validated before it is
# published to the canonical name -- publish-then-retract loses the retraction if the run is
# interrupted. Unset (the standalone case) this is the canonical path, exactly as before.
mf <- if (nzchar(Sys.getenv("PE_MASTER_IN"))) Sys.getenv("PE_MASTER_IN") else file.path(RC, "phase4_annotation_master.csv")
if (!file.exists(mf)) rerun_abort(sprintf("43b: annotation master absent: %s — stage it (15_stage_annotation.R, which runs 09_merge_master.py with the output pointed at %s) before 44", mf, RC))
master <- fread(mf)
if (!"gene" %in% names(master)) rerun_abort("43b: phase4_annotation_master.csv has no 'gene' column")
# 2026-06-30 (editor enhancement #3): a master with gene rows but NO annotation VALUE columns would
# silently NA every annotation downstream (left_join adds nothing). Require >=1 expected value column.
if (!length(intersect(c("drug_tier_direct","drug_tier_any_evidence","tissue_class","tissue_breadth_hpa","phewas_consensus"), names(master))))
  rerun_abort("43b: annotation master has gene rows but NONE of the expected value columns (drug_tier_direct/drug_tier_any_evidence/tissue_class/tissue_breadth_hpa/phewas_consensus) — master malformed; would silently NA all annotations")
# The breadth axis specifically. tissue_class alone would satisfy the check above, and a master built
# before 2026-08-04 has it -- but no tissue_breadth_hpa, so every breadth statement in the Results and
# in Supplementary Figure S8 would silently read as empty rather than failing.
if (!"tissue_breadth_hpa" %in% names(master))
  rerun_abort("43b: annotation master has no tissue_breadth_hpa column -> it predates the HPA breadth axis; re-run phase4/09_merge_master.py (via cascade/15_stage_annotation.R) so the breadth counts in the Results and Figure S8 are not silently empty")
miss <- setdiff(strong_genes, unique(as.character(master$gene)))
if (length(miss))
  rerun_abort(sprintf("43b: annotation master MISSING %d/%d strong gene(s): %s -- annotate before 44; do NOT proceed with silent gaps",
                      length(miss), length(strong_genes), paste(miss, collapse = ", ")))
# 2026-06-30 (gate obj #3): report stale-extra master genes (not in current strong set) — benign but
# made visible (the master is built from the cached gene universe, not exactly the current strong set).
.extra <- setdiff(unique(as.character(master$gene)), strong_genes)
if (length(.extra)) cat(sprintf("43b note: master has %d gene(s) NOT in current strong set (stale-extra, benign annotation rows): %s\n",
                                length(.extra), paste(head(.extra, 20), collapse = ", ")))
# 2026-07-01 (gate finding 3): a strong gene present in the master but ABSENT from the SOURCE annotation
# caches would carry unrefreshed DEFAULT values (drug_tier=4 / not_queryable), i.e. missing external
# assessment silently encoded as "no evidence". Assert per-source coverage of the strong genes.
.cachedir <- file.path(if (exists("BASE_DIR")) BASE_DIR else dirname(RESULTS_DIR), "cache")
# The identifier map every source was supposed to be queried with. 01_gene_id_map.R rebuilds it from
# the certified Phase 1 table each run, so it is the current generation's answer by construction;
# the per-source check below compares each cache's own identifier columns against it.
.idmapf <- file.path(.cachedir, "phase4_gene_ids.csv")
if (!file.exists(.idmapf)) rerun_abort(sprintf("43b: identifier map absent: %s -> run phase4/01_gene_id_map.R first", .idmapf))
.idmap <- fread(.idmapf)
if (!"gene" %in% names(.idmap)) rerun_abort("43b: phase4_gene_ids.csv has no 'gene' column -> malformed identifier map")
# 2026-07-02: gwasatlas is EXEMPT from full-coverage (it is a TARGETED secondary source -- 08c only queries the
# SNPs that ALREADY had an upstream PheWAS hit, not all instruments; a strong gene whose SNP had NO upstream hit
# is legitimately absent = no non-CKM safety flags, NOT an unrefreshed-cache gap). It is reported below instead.
for (.src in c("phase4_gene_ids.csv", "phase4_opentargets.csv", "phase4_gtex.csv", "phase4_chembl.csv",
               "phase4_dgidb.csv", "phase4_hpa.csv", "phase4_phewas.csv", "phase4_phewas_ieu.csv")) {  # gwasatlas exempt (targeted secondary)
  .sp <- file.path(.cachedir, .src)
  if (!file.exists(.sp)) rerun_abort(sprintf("43b: REQUIRED annotation source cache absent: %s -> re-run the source before Phase 5 (a missing cache would silently annotate every strong gene as 'no evidence')", .src))  # 2026-07-02 (Codex d): fail, don't skip a required source
  .sd <- fread(.sp); if (!"gene" %in% names(.sd)) rerun_abort(sprintf("43b: annotation source %s has no 'gene' column -> malformed cache", .src))
  # 2026-07-02 (Codex phase4 #4b): SNP-level caches carry multi-gene cells like "KCTD5;PRSS27" -> split on ;/,
  # before matching the strong set (else the individual PRSS27/KCTD5 wrongly look ABSENT).
  .cache_genes <- unique(trimws(unlist(strsplit(as.character(.sd$gene), "[;,]"))))
  .sm <- setdiff(strong_genes, .cache_genes)
  if (length(.sm))
    rerun_abort(sprintf("43b: %d strong gene(s) ABSENT from annotation source %s -> unrefreshed cache = default/uninformative annotation read as 'no evidence'; re-query the sources before Phase 5: %s",
                        length(.sm), .src, paste(.sm, collapse = ", ")))
  # 2026-07-02 (Codex phase4 #4/#4a): a strong gene present but with ONLY a FAILED query status (api_error) is a
  # TECHNICAL failure that 09_merge_master silently encodes as "no evidence"/"not_queried". Assert each strong
  # gene has >=1 successfully-queried row. Includes ieu_query_status (IEU PheWAS). (base-R vectors -> no DT NSE)
  .stcol <- intersect(c("query_status", "atlas_query_status", "ieu_query_status", "status"), names(.sd))
  if (length(.stcol)) {
    .stat     <- as.character(.sd[[.stcol[1]]])
    .ok_genes <- unique(trimws(unlist(strsplit(as.character(.sd$gene[!grepl("error|fail", .stat, ignore.case = TRUE)]), "[;,]"))))
    .sfail    <- setdiff(intersect(strong_genes, .cache_genes), .ok_genes)
    if (length(.sfail))
      rerun_abort(sprintf("43b: %d strong gene(s) have ONLY FAILED %s in %s (technical failure silently read as 'no evidence') -> re-query before Phase 5: %s",
                          length(.sfail), .stcol[1], .src, paste(.sfail, collapse = ", ")))
  }
  # AND THAT IT WAS QUERIED FOR THE RIGHT GENE. Coverage and query status both pass on a cache built
  # for a DIFFERENT identifier: 09_merge_master.py indexes these files by gene SYMBOL alone, so a row
  # labelled with the right symbol but fetched against the wrong accession merges as that gene's
  # evidence. Not hypothetical -- the July cache carried ACP1 with ENSG00000004779, which is NDUFAB1,
  # and reported drug tier 1 / 3 drugs / phase 4 (METFORMIN) under the ACP1 label. The identifier map
  # is rebuilt from Phase 1 each run; a source cache is only current if it agrees with it.
  #
  # ONLY where the identifier is what the source was actually QUERIED with. Open Targets sends the
  # Ensembl ID (03_opentargets.py:48-49) and ChEMBL the UniProt accession (04_chembl.py:20-23), so a
  # disagreement there means the row was fetched for another gene. GTEx resolves by SYMBOL and only
  # falls back to the Ensembl ID (06_gtex.py:123,135), and HPA queries gene_name (07_hpa.py:30) --
  # in those two the identifier column is metadata copied from the map, so comparing it would abort
  # on a correct cache and force an unnecessary refetch. ACP1's GTEx row is exactly that case: it
  # carries the stale ensembl_id but was fetched on gencode ENSG00000143727.15, which is right.
  .QKEY <- c("phase4_opentargets.csv" = "ensembl_id", "phase4_chembl.csv" = "uniprot_id")
  for (.idc in intersect(if (.src %in% names(.QKEY)) .QKEY[[.src]] else character(0), names(.sd))) {
    if (!.idc %in% names(.idmap)) next
    .j <- merge(unique(data.table::data.table(gene = trimws(as.character(.sd$gene)),
                                              cache_id = trimws(as.character(.sd[[.idc]])))),
                unique(data.table::data.table(gene = trimws(as.character(.idmap$gene)),
                                              map_id = trimws(as.character(.idmap[[.idc]])))), by = "gene")
    .bad <- .j[gene %in% strong_genes & nzchar(cache_id) & !is.na(cache_id) &
               nzchar(map_id) & !is.na(map_id) & cache_id != map_id]
    if (nrow(.bad))
      rerun_abort(sprintf("43b: %d strong gene(s) in %s carry a %s that disagrees with the current identifier map (annotation fetched for a different gene reads as this gene's evidence) -> re-query the source: %s",
                          nrow(.bad), .src, .idc,
                          paste(sprintf("%s cache=%s map=%s", .bad$gene, .bad$cache_id, .bad$map_id), collapse = "; ")))
  }
}
# AND THAT THE PheWAS CACHE IS KEYED ON THIS GENERATION'S LEAD SNPs. Coverage and query status both
# pass on a cache built for a DIFFERENT variant: 08_phewas.py takes its SNPs solely from lead_snps in
# table_s_pair_evidence.csv (08_phewas.py:128-149), and 09_merge_master.py then folds EVERY cache row
# carrying a gene's symbol into one gene-level phewas_class (09_merge_master.py:111-118, 252). So a row
# left from an earlier generation contributes that variant's safety flags to this gene, and a lead SNP
# that changed is never queried at all -- both silently. Present on disk today: the July cache carries
# NOTCH2/rs587635864, a SNP in no pair-evidence row.
.pef <- file.path(RC, "table_s_pair_evidence.csv")
if (!file.exists(.pef)) rerun_abort(sprintf("43b: %s absent -> phase4/12_rebuild_pair_evidence.R must run before the PheWAS cache can be checked against the current lead SNPs", .pef))
.pe <- fread(.pef)
if (!all(c("gene", "lead_snps") %in% names(.pe)))
  rerun_abort("43b: table_s_pair_evidence.csv lacks gene/lead_snps -> cannot verify the PheWAS caches are keyed on the current lead SNPs")
.L <- unique(.pe[, .(snp = trimws(unlist(strsplit(as.character(lead_snps), ",")))), by = gene][nzchar(snp) & gene %in% strong_genes, .(gene, snp)])
# BOTH PheWAS caches, not just the Open Targets one. 09_merge_master.py folds the IEU rows in by gene
# on its own path (09_merge_master.py:82-89, 283), so an IEU cache left on an earlier generation's
# variants reaches the master even when the OT cache is current. Concrete on disk today: PAM's current
# instrument is rs149802978, while the July IEU cache holds PAM/rs80134320 and PAM/rs77460585, both
# flagged, 20 non-CKM associations between them, none of which belong to this generation's variant.
for (.pwf in c("phase4_phewas.csv", "phase4_phewas_ieu.csv")) {
  .pw <- fread(file.path(.cachedir, .pwf))
  if (!"snp" %in% names(.pw)) rerun_abort(sprintf("43b: %s lacks a snp column -> cannot verify it is keyed on the current lead SNPs", .pwf))
  # Duplicate rows are NOT collapsed before comparing. 09_merge_master.py:252 sums and classifies every
  # row it finds for a gene, so a cache holding an old and a new row for the same SNP doubles that
  # gene's association and safety counts; unique() here would have let exactly that through.
  .dup <- .pw[, .N, by = .(snp = trimws(as.character(snp)))][N > 1L]
  if (nrow(.dup))
    rerun_abort(sprintf("43b: %s has %d SNP(s) on more than one row (%s) -- the merger sums every row, so the gene-level association and safety counts would be inflated; rebuild the cache with phase4/08_phewas.py",
                        .pwf, nrow(.dup), paste(sprintf("%s x%d", .dup$snp, .dup$N), collapse = ", ")))
  .C <- unique(.pw[, .(gene = trimws(unlist(strsplit(as.character(gene), "[;,]")))), by = .(snp = trimws(as.character(snp)))][nzchar(gene) & gene %in% strong_genes, .(gene, snp)])
  .mi <- fsetdiff(.L, .C); .ex <- fsetdiff(.C, .L)
  if (nrow(.mi) || nrow(.ex))
    rerun_abort(sprintf("43b: %s is not keyed on this generation's lead SNPs -> %d current lead SNP(s) never queried (%s); %d cached SNP(s) left from an earlier generation, whose flags merge as this gene's evidence (%s) -- re-run phase4/08_phewas.py (and 08b_phewas_ieu.R, which takes its SNP list from the OT cache)",
                        .pwf,
                        nrow(.mi), if (nrow(.mi)) paste(sprintf("%s/%s", .mi$gene, .mi$snp), collapse = ", ") else "none",
                        nrow(.ex), if (nrow(.ex)) paste(sprintf("%s/%s", .ex$gene, .ex$snp), collapse = ", ") else "none"))
}
# 2026-07-02: gwasatlas coverage is REPORT-ONLY (targeted secondary source, flagged-SNP subset only)
.gwa <- file.path(.cachedir, "phase4_phewas_gwasatlas.csv")
if (file.exists(.gwa)) { .gd <- fread(.gwa)
  if ("gene" %in% names(.gd)) {
    .gwg <- length(intersect(strong_genes, unique(trimws(unlist(strsplit(as.character(.gd$gene), "[;,]"))))))
    cat(sprintf("43b [info] gwasatlas (targeted; only SNPs with an upstream PheWAS hit) covers %d/%d strong genes -- partial by design\n", .gwg, length(strong_genes)))
    # 2026-07-02 (Codex b): coverage is exempt, but a gwasatlas api_error for a SNP that WAS queried is a
    # TECHNICAL failure (not a legit no-hit) -> still fail on it.
    if ("atlas_query_status" %in% names(.gd)) {
      .gwerr <- intersect(strong_genes, unique(trimws(unlist(strsplit(as.character(.gd$gene[grepl("error|fail", as.character(.gd$atlas_query_status), ignore.case = TRUE)]), "[;,]")))))
      if (length(.gwerr)) rerun_abort(sprintf("43b: %d strong gene(s) have a FAILED atlas_query_status in gwasatlas (technical failure, not a no-hit) -- re-query: %s", length(.gwerr), paste(.gwerr, collapse = ", ")))
    }
  }
}
cat("43b: all strong genes present + successfully-queried in the 8 required annotation source caches (gene_ids/opentargets/gtex/chembl/dgidb/hpa/phewas/phewas_ieu); gwasatlas report-only\n")
cat(sprintf("43b OK: annotation master covers all %d corrected strong genes.\n", length(strong_genes)))
