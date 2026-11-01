# Paper 3.1 rewrite (R2) -- Additional file 3, Figure S4: DPEP1 locus and attribution (moved here from the old Figure 3d; S2-D24). Illustrates SM S-R8.
#
# STATE (S2-D25 fix round): DRAWN. The inputs that were not approved for AF3_S4 now have approved rows in registry/APPROVED_SOURCE_MANIFEST.tsv
# (traits.tsv, PREP_STAMP.json, conditional/signals.tsv, the 503-person dosage matrix, the attribution grid, the locus annotation, the Ensembl window
# response and its fetch record, the two Ensembl info responses and the eGFR window extract). Every file is opened through af3_read()/af3_source()
# only, so gate G12 re-checks each one by path and md5. Two cross-checks of the old Figure 3d generator are NOT repeated because the files they read
# are not approved for this figure: the 192-MB CKDGen source (the window extract's md5 is checked against the fetch record instead) and the V4 DPEP1
# Ensembl cache (strands and release come from the approved Ensembl response and its fetch record).
#
# Content (association tracks, -log10 p, as Figure 4b of the main set):
#   DPEP1 protein (UKB-PPP), cysteinylglycine (CLSA metabolite GWAS), creatinine eGFR (CKDGen 2019, European ancestry); points coloured by LD
#   (1000 Genomes EUR, n = 503) with signal A or B (r2 >= 0.5 dark, 0.2-0.5 light); signal A and B of the attribution grid's primary setting
#   (palindromic variants kept) marked with diamonds and dashed lines, the top shared variant of the conditional colocalization at A with a square and a
#   dotted line in every track; gene models = protein-coding transcript extents (grey), canonical transcript (black), and, dashed, the part of the
#   CHMP1A gene record held only by its transcript without a coding sequence (which contains the top shared variant).
#   Evidence labels (S2-D25 visual check): the header states WHICH method gave WHICH result. Under the single-causal-variant model (coloc.abf)
#   DPEP1-cysteinylglycine colocalized (PP.H4 0.998, ledger L072, declared prior); the fit that is inadmissible is the reference-LD multi-signal
#   fine-mapping, which is not cited as support. The attribution between DPEP1 and CHMP1A is unresolved. (The previous legend reversed the two methods.)
# Run:  Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_FigS4.R"

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
source(paste0(.P31R_V3, "/R2_rewrite/figures/scripts/P31R_fig_common.R"))
source(paste0(.P31R_V3, "/R2_rewrite/figures/af3/P31R_af3_common.R"))
suppressPackageStartupMessages({ library(ggtext) })
UNIT <- "FigureS4"
af3_begin(UNIT)
EN <- "\u2013"
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))
INK <- "grey15"
fmtd <- function(d) function(v) formatC(v, format = "f", digits = d)
# a literal that is not a tick (scale constants of the key, the log base)
lit <- function(text, why, anchor = NULL, count = 1L) {
  REG$literals[[length(REG$literals) + 1]] <<- Filter(Negate(is.null), list(text = text, why = why, count = as.integer(count), anchor = anchor))
  text
}
# axis ticks, one registry entry per printed break; stacked tracks may join tick lines, so the anchor allows other numbers on the line
TICK_ANCHOR <- "^(?:[\u2212+-]?[0-9.,]+\\s)*{}(?:\\s[\u2212+-]?[0-9.,]+)*$"
ticks_set <- function(axis, breaks, lims, why, labels = p31r_lab(breaks), transform = "linear") {
  stopifnot(all(breaks >= lims[1] & breaks <= lims[2]), !anyDuplicated(breaks), length(labels) == length(breaks))
  for (t in labels) REG$literals[[length(REG$literals) + 1]] <<- list(text = t, why = paste0("axis_tick: ", why), axis = axis, transform = transform, count = 1L, anchor = TICK_ANCHOR)
  labels
}

# ================================================================== data (approved sources only)
FM <- "V4_local/staging/l3-finemap/DPEP1_16q24"
DOS_REL <- "V4_local/raw/l3-finemap/DPEP1_16q24__ref_A1dosage_503xP.i8"
DD <- "R2_rewrite/figures/data/Figure_3"
h <- af3_read(file.path(FM, "harmonised.tsv.gz"), "", sep = "\t")
H_F <- .abs_path(file.path(FM, "harmonised.tsv.gz"))
stamp <- jsonlite::fromJSON(af3_source(file.path(FM, "PREP_STAMP.json"), ""))
stopifnot(unname(tools::md5sum(H_F)) == stamp$outputs[["harmonised.tsv.gz"]], nrow(h) == stamp$n_snps, !anyDuplicated(h$snp))
tr <- af3_read(file.path(FM, "traits.tsv"), "", sep = "\t")
stopifnot(unname(tools::md5sum(.abs_path(file.path(FM, "traits.tsv")))) == stamp$outputs[["traits.tsv"]])
MT <- tr[kind == "metabolite" & label == "cysteinylglycine", tid]
stopifnot(length(MT) == 1L, tr[tid == MT, d1_H4] == num("L072"), "P__DPEP1" %in% tr$tid)
grid <- af3_read("V4_local/source_data/l3_metabolite/bridge_class_robustness__ens116pc.csv", "", encoding = "UTF-8")[protein == "DPEP1" & metabolite == "cysteinylglycine"]
stopifnot(nrow(grid) == num("L121"), sum(grid$is_primary == TRUE) == 1L, grid[is_primary == TRUE, pal] == "freq")
gf <- grid[pal == "freq"]; gd <- grid[pal == "drop"]
stopifnot(uniqueN(gf$signal_A) == 1L, uniqueN(gf$signal_B) == 1L, uniqueN(gf$A2_lead_snp) == 1L, uniqueN(gf$A2_signal_snp) == 1L, gf$A2_signal_snp[1] == gf$signal_A[1])
TOP <- h[snp == gf$A2_lead_snp[1], rsid_ref]; SIGA <- gf$signal_A[1]; SIGB <- gf$signal_B[1]
stopifnot(length(TOP) == 1L)
ann <- af3_read("V4_local/source_data/l3_metabolite/dpep1_locus_annotation.csv", "")
stopifnot(all(c(TOP, SIGA, SIGB) %in% ann$rsid),
          grepl("signal_A\\[[a-z,]*freq", ann[rsid == SIGA, grid_roles]), grepl("signal_B\\[freq\\]", ann[rsid == SIGB, grid_roles]),
          grepl("A2_lead_snp\\[freq\\]", ann[rsid == TOP, grid_roles]), isTRUE(as.logical(ann[rsid == SIGB, study_palindromic])))
aa <- sub(".*p\\.([A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}).*", "\\1", ann[rsid == SIGB, vep_canonical])
stopifnot(grepl("missense_variant", ann[rsid == SIGB, vep_canonical]), grepl("^[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}$", aa))
stopifnot(identical(c(TOP, aa), unlist(p31r_ledger("L073"))))                      # the printed variant names are the ledger's
sigs <- af3_read(file.path(FM, "conditional/signals.tsv"), "", sep = "\t")[tid == "P__DPEP1"]
stopifnot(identical(sigs[order(signal)][1:3, rsid], c(SIGA, SIGB, TOP)))
MK <- data.table(rsid = c(SIGA, SIGB, TOP), role = c("A", "B", "T"), lab = c(paste("A", SIGA), sprintf("B %s (%s)", SIGB, aa), TOP))
MK <- merge(MK, ann[, .(rsid, pos38)], by = "rsid", sort = FALSE)
MK <- merge(MK, h[, .(rsid = rsid_ref, pos38h = pos38)], by = "rsid", sort = FALSE)
stopifnot(nrow(MK) == 3L, all(MK$pos38 == MK$pos38h))
# Ensembl window response: md5s equal the fetch record; release and assembly from the approved info responses
fr <- jsonlite::fromJSON(af3_source(file.path(DD, "FETCH_RECORD.json"), ""))
for (nm in c("overlap", "info_data", "info_software", "info_assembly")) {
  f <- af3_source(file.path(DD, paste0("ensembl_", nm, ".json")), "")
  stopifnot(unname(tools::md5sum(f)) == fr$ensembl$responses[[nm]]$md5)
}
OV_F <- .abs_path(file.path(DD, "ensembl_overlap.json"))
stopifnot(jsonlite::fromJSON(.abs_path(file.path(DD, "ensembl_info_software.json")))$release == fr$ensembl$release, fr$ensembl$release == 116,
          jsonlite::fromJSON(.abs_path(file.path(DD, "ensembl_info_assembly.json")))$default_coord_system_version == "GRCh38", fr$ensembl$assembly == "GRCh38")
WIN <- as.numeric(fr$ensembl$window[2:3])
stopifnot(all(MK$pos38 > WIN[1] & MK$pos38 < WIN[2]))
# with palindromic variants dropped signal B moves OUTSIDE the plotted window; the further conditional signals lie outside too
stopifnot(uniqueN(gd$signal_A) == 1L, gd$signal_A[1] == SIGA, uniqueN(gd$signal_B) == 1L, gd$signal_B[1] != SIGB,
          !(ann[rsid == gd$signal_B[1], pos38] %between% WIN), h[rsid_ref == SIGB, palindromic] == 1,
          all(!(h[rsid_ref %in% sigs[signal > 3, rsid], pos38] %between% WIN)))
# the eGFR window extract: the file and row count equal the fetch record (the 192-MB source is not opened here: it has no approved row for this figure)
eg <- af3_read(file.path(DD, "ckdgen_eGFR_EA_chr16_window_grch37.tsv"), "", sep = "\t", colClasses = list(character = c("Allele1", "Allele2")))
stopifnot(unname(tools::md5sum(.abs_path(file.path(DD, "ckdgen_eGFR_EA_chr16_window_grch37.tsv")))) == fr$egfr$extract_md5, nrow(eg) == fr$egfr$rows,
          identical(fr$egfr$window_from, file.path(FM, "PREP_STAMP.json")))
eg[, `:=`(A1 = toupper(Allele1), A2 = toupper(Allele2))]
j <- merge(h[, .(snp, pos37, pos38, rsid = rsid_ref, a1, a2, af_ref, zP = get("z__P__DPEP1"), zM = get(paste0("z__", MT)))],
           eg[, .(pos37 = Pos_b37, rs_e = RSID, A1, A2, Freq1, Effect, StdErr, P = `P-value`)], by = "pos37", all.x = TRUE)
j[, s := fifelse(A1 == a1 & A2 == a2, 1, fifelse(A1 == a2 & A2 == a1, -1, NA_real_))]
j[, f1 := fifelse(s == 1, Freq1, 1 - Freq1)]
j[!is.na(s) & abs(f1 - af_ref) > stamp$af_diff_max, s := NA_real_]          # same allele-frequency rule as the fine-mapping panel
j <- j[order(pos38, is.na(s))][!duplicated(pos38)]
stopifnot(nrow(j) == nrow(h), all(j[!is.na(s), rs_e == rsid]), all(!is.na(j[rsid %in% MK$rsid, s])))
# -log10 p: protein and metabolite from z (two-sided normal, log scale: DPEP1 |z| reaches ~105), creatinine eGFR from CKDGen's own p
mlp <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
j[, zE := s * Effect / StdErr]
j[, sP := sign(zP)]
j[, `:=`(sM = zM * sP, sE = -zE * sP)]                                      # signed, DPEP1-raising allele; eGFR* = -ln(eGFR)
j[, `:=`(yP = mlp(zP), yM = mlp(zM), yE = -log10(P))]
stopifnot(max(abs(j$yE - mlp(j$zE)), na.rm = TRUE) < 0.05, all(j$P > 0, na.rm = TRUE))
# the direction the panel note states: the DPEP1-raising allele LOWERS cysteinylglycine at A, B and the top shared variant, and lowers eGFR at A and B
stopifnot(all(j[rsid %in% MK$rsid, sM] < 0), all(j[rsid %in% c(SIGA, SIGB), sE] > 0))
# LD with the two signal variants: the fine-mapping panel's own 1000 Genomes EUR dosage (n = L122), SNP order = harmonised
DOS_F <- af3_source(DOS_REL, "")
stopifnot(unname(tools::md5sum(DOS_F)) == stamp$outputs[[paste0("V3/", DOS_REL)]], stamp$ref_n == num("L122"))
NREF <- stamp$ref_n
X <- matrix(as.numeric(readBin(DOS_F, "integer", n = NREF * nrow(h), size = 1, signed = TRUE)), nrow = NREF, ncol = nrow(h), byrow = TRUE)
stopifnot(file.info(DOS_F)$size == NREF * nrow(h), max(abs(colMeans(X) / 2 - h$af_ref)) < 1e-5)
jw <- j[pos38 >= WIN[1] & pos38 <= WIN[2]]
iw <- match(jw$snp, h$snp); iA <- match(SIGA, h$rsid_ref); iB <- match(SIGB, h$rsid_ref)
stopifnot(!anyNA(iw), !is.na(iA), !is.na(iB))
jw[, `:=`(rA = suppressWarnings(cor(X[, iw], X[, iA]))[, 1], rB = suppressWarnings(cor(X[, iw], X[, iB]))[, 1])]
jw[is.na(rA), rA := 0]; jw[is.na(rB), rB := 0]
R2CUT <- c(hi = 0.5, lo = 0.2)
jw[, ld := fifelse(rA^2 >= rB^2,
                   fifelse(rA^2 >= R2CUT[["hi"]], "A_hi", fifelse(rA^2 >= R2CUT[["lo"]], "A_lo", "other")),
                   fifelse(rB^2 >= R2CUT[["hi"]], "B_hi", fifelse(rB^2 >= R2CUT[["lo"]], "B_lo", "other")))]
stopifnot(jw[rA^2 >= R2CUT[["hi"]] & rB^2 >= R2CUT[["hi"]], .N] == 0, jw[rsid == SIGA, ld] == "A_hi", jw[rsid == SIGB, ld] == "B_hi",
          jw[rsid == TOP, rA^2 < R2CUT[["lo"]] & rA^2 > rB^2], jw[which.max(yE), ld] == "A_lo")
stopifnot(all(is.na(jw[is.na(yE), rs_e])))
TRK_LAB <- c(P = "DPEP1 protein\n(UKB-PPP)", M = "Cysteinylglycine\n(CLSA)", E = "Creatinine eGFR\n(CKDGen)")
TRK <- rbind(jw[, .(trk = "P", pos = pos38, y = yP, rsid, ld)], jw[, .(trk = "M", pos = pos38, y = yM, rsid, ld)], jw[!is.na(yE), .(trk = "E", pos = pos38, y = yE, rsid, ld)])
LDLEV <- c("A_hi", "A_lo", "B_hi", "B_lo", "other")
TRK[, ld := factor(ld, levels = LDLEV)]
setorder(TRK, -ld)
# gene models: protein-coding genes of the window; grey line = extent of the gene's PROTEIN-CODING transcripts, black = canonical transcript
ov <- as.data.table(jsonlite::fromJSON(OV_F)[, c("feature_type", "id", "Parent", "external_name", "biotype", "start", "end", "strand", "is_canonical")])
gn <- ov[feature_type == "gene" & biotype == "protein_coding", .(gene = id, g = external_name, strand)]
stopifnot(setequal(gn$g, c("DPEP1", "CHMP1A", "SPATA33")))
pcx <- ov[feature_type == "transcript" & biotype == "protein_coding" & Parent %in% gn$gene, .(x0 = min(start), x1 = max(end)), by = .(gene = Parent)]
txc <- ov[feature_type == "transcript" & is_canonical %in% 1 & Parent %in% gn$gene, .(tx = id, gene = Parent, t0 = start, t1 = end, biotype)]
stopifnot(nrow(txc) == nrow(gn), !anyDuplicated(txc$gene), all(txc$biotype == "protein_coding"), nrow(pcx) == nrow(gn))
gn <- merge(merge(gn, pcx, by = "gene"), txc, by = "gene")
stopifnot(all(gn$x0 <= gn$t0 & gn$t1 <= gn$x1))
exn <- ov[feature_type == "exon" & Parent %in% gn$tx, .(tx = Parent, e0 = start, e1 = end)]
stopifnot(all(gn$tx %in% exn$tx))
stopifnot(gn[g == "DPEP1", strand] %in% c(-1, 1), gn[g == "CHMP1A", strand] %in% c(-1, 1))
.tp <- MK[rsid == TOP, pos38]
stopifnot(gn[g == "DPEP1", x1] < .tp, .tp < gn[g == "CHMP1A", x0],
          .tp - gn[g == "DPEP1", x1] == as.numeric(ann[rsid == TOP, DPEP1_pc_tx_bp]),
          gn[g == "CHMP1A", x0] - .tp == as.numeric(ann[rsid == TOP, CHMP1A_pc_tx_bp]),
          !any(gn$t0 <= .tp & gn$t1 >= .tp))
.ct <- ov[feature_type == "transcript" & Parent == gn[g == "CHMP1A", gene] & start <= .tp & end >= .tp, .(id, biotype, start, end)]
stopifnot(nrow(.ct) == 1L, .ct$biotype == "protein_coding_CDS_not_defined", as.numeric(ann[rsid == TOP, CHMP1A_gene_bp]) == 0,
          startsWith(ann[rsid == TOP, CHMP1A_containing_tx], .ct$id))
grec <- ov[feature_type == "gene" & id %in% gn$gene, .(gene = id, r0 = start, r1 = end)]
gn <- merge(gn, grec, by = "gene")
ext <- rbind(gn[r0 < x0 - 100, .(g, a = r0, b = x0, strand)], gn[r1 > x1 + 100, .(g, a = x1, b = r1, strand)])
stopifnot(nrow(ext) == 1L, ext$g == "CHMP1A", ext$a == .ct$start, ext$a < .tp, .tp < ext$b)

# ================================================================== drawing
MB <- 1e6
X_D <- seq(89.60, 89.68, by = 0.02); stopifnot(isTRUE(all.equal(range(X_D), WIN / MB)))
xs_d <- function(lab) scale_x_continuous(limits = WIN / MB, breaks = X_D, labels = lab, expand = c(0, 0))
VAB <- MK[role != "T"]; VT <- MK[role == "T"]
COL_SIG <- PAL_SIGNAL
vl_ab <- geom_vline(data = VAB, aes(xintercept = pos38 / MB, colour = role), linetype = "dashed", linewidth = 0.35, show.legend = FALSE)
vl_t  <- geom_vline(data = VT, aes(xintercept = pos38 / MB), linetype = "dotted", colour = "black", linewidth = 0.45)
R_HI <- lit(format(R2CUT[["hi"]]), "scale: LD colour bins (r2 with a DPEP1 signal variant)", anchor = "(?:\u2265 |\u2013){}", count = 4L)
R_LO <- lit(format(R2CUT[["lo"]]), "scale: LD colour bins (r2 with a DPEP1 signal variant)", anchor = "{}\u2013", count = 2L)
LDLAB <- c(A_hi = sprintf("A, r\u00b2 \u2265 %s", R_HI), A_lo = sprintf("A, r\u00b2 %s\u2013%s", R_LO, R_HI), B_hi = sprintf("B, r\u00b2 \u2265 %s", R_HI),
           B_lo = sprintf("B, r\u00b2 %s\u2013%s", R_LO, R_HI), other = "Other variants")
stopifnot(all(LDLEV %in% TRK$ld))
COL_LD <- c(A_hi = PAL_SIGNAL[["A"]], A_lo = PAL_SIGNAL_LIGHT[["A"]], B_hi = PAL_SIGNAL[["B"]], B_lo = PAL_SIGNAL_LIGHT[["B"]], other = PAL_NEUTRAL[["other"]])
LOG_TTL <- "\u2212log<sub><span style='font-size:7pt'>10</span></sub> <i>p</i>"
invisible(lit("10", "scale: base of the logarithm in the y title (-log10 p)", anchor = "log {}"))
DIR_NOTE <- "DPEP1-raising allele:
lower cysteinylglycine at A and B"
track <- function(tk, first = FALSE, ytitle = FALSE) {
  d <- TRK[trk == tk]; d[, strip_lab := factor(TRK_LAB[[tk]])]
  hiab <- d[rsid %in% VAB$rsid]; hiab[, role := VAB$role[match(rsid, VAB$rsid)]]
  hit <- d[rsid %in% VT$rsid]
  lim <- range(c(0, d$y)); EXP <- c(0.08, if (tk == "M") 0.62 else 0.08); edge <- lim + c(-EXP[1], EXP[2]) * diff(lim)
  brk <- scales::breaks_extended(n = if (tk == "M") 3 else 4)(lim); brk <- brk[brk >= edge[1] & brk <= edge[2] - 0.08 * diff(edge)]
  p <- ggplot(d, aes(pos / MB, y)) + vl_ab + vl_t +
    geom_point(aes(colour = ld), size = 0.6, stroke = 0) +
    geom_point(data = hit, shape = 15, colour = "black", size = 1.5) +
    geom_point(data = hiab, aes(fill = role), shape = 23, colour = "black", size = 2.1, stroke = 0.4, show.legend = FALSE) +
    facet_grid(strip_lab ~ ., switch = "y") +
    scale_colour_manual(values = c(COL_LD, COL_SIG), breaks = LDLEV, labels = LDLAB[LDLEV], name = "DPEP1 signals", drop = FALSE,
                        guide = guide_legend(nrow = 1, override.aes = list(size = 1.8))) +
    scale_fill_manual(values = COL_SIG, guide = "none") +
    xs_d(NULL) +
    scale_y_continuous(limits = lim, breaks = brk, labels = ticks_set(paste0("S4_y_", tk), brk, lim, paste0("S4 y, track ", tk)), expand = expansion(mult = EXP)) +
    labs(x = NULL, y = if (ytitle) LOG_TTL else NULL) + theme_p31r() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.title.y = element_markdown(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, lineheight = 0.95),
          plot.margin = margin(2.5, 12, 2.5, 4),
          legend.position = if (first) "top" else "none", legend.justification = "right", legend.title = element_text(size = BODY_PT, face = "bold"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0), legend.key.size = unit(8, "pt"))
  if (tk == "M") p <- p + annotate("text", x = VT$pos38 / MB + 0.0006, y = Inf, label = DIR_NOTE, hjust = 0, vjust = 1.15, size = ANNOT, family = BASE_FONT, lineheight = 0.95)
  p
}
pd1 <- track("P", first = TRUE)
pd2 <- track("M", ytitle = TRUE)
pd3 <- track("E")
GY <- 2.25
gn[, y := fifelse(strand > 0, 1, -1)]
gseg <- gn[, .(x = pmax(x0, WIN[1]) / MB, xend = pmin(x1, WIN[2]) / MB, y)]
cseg <- gn[, .(x = pmax(t0, WIN[1]) / MB, xend = pmin(t1, WIN[2]) / MB, y)]
exr <- merge(exn, gn[, .(tx, y)], by = "tx")[, .(xmin = pmax(e0, WIN[1]) / MB, xmax = pmin(e1, WIN[2]) / MB, y)][xmax > xmin]
glab <- gn[, .(x = (pmax(x0, WIN[1]) + pmin(x1, WIN[2])) / 2 / MB, y = y * GY, lab = fifelse(strand > 0, paste(g, "\u2192"), paste("\u2190", g)))]
MK[, `:=`(hj = fifelse(role == "T", 0, 1), dx = fifelse(role == "T", 0.0006, -0.0006))]
X_D_LAB <- ticks_set("S4_x_chr16", X_D, WIN / MB, "S4 x, chromosome 16 position (Mb, GRCh38)", labels = scales::label_number(accuracy = 0.01, style_negative = "minus")(X_D))
eseg <- ext[, .(x = a / MB, xend = b / MB, y = fifelse(strand > 0, 1, -1))]
pdg <- ggplot() + vl_ab + vl_t +
  geom_segment(data = eseg, aes(x = x, xend = xend, y = y, yend = y), colour = "grey55", linewidth = 0.35, linetype = "22") +
  geom_segment(data = gseg, aes(x = x, xend = xend, y = y, yend = y), colour = "grey55", linewidth = 0.35) +
  geom_segment(data = cseg, aes(x = x, xend = xend, y = y, yend = y), colour = INK, linewidth = 0.35) +
  geom_rect(data = exr[(xmax - xmin) * MB >= 150], aes(xmin = xmin, xmax = xmax, ymin = y - 0.4, ymax = y + 0.4), fill = INK, colour = NA) +
  geom_segment(data = exr[(xmax - xmin) * MB < 150], aes(x = (xmin + xmax) / 2, xend = (xmin + xmax) / 2, y = y - 0.4, yend = y + 0.4), colour = INK, linewidth = 0.3) +   # exons under 150 bp (below 0.25 pt at this scale) as ticks
  geom_text(data = glab, aes(x = x, y = y, label = lab), size = ANNOT, family = BASE_FONT, fontface = "italic") +
  geom_text(data = MK, aes(x = pos38 / MB + dx, y = 4.45, label = lab, hjust = hj), colour = INK, size = ANNOT, family = BASE_FONT) +
  scale_colour_manual(values = COL_SIG, guide = "none") +
  xs_d(X_D_LAB) + scale_y_continuous(limits = c(-3.1, 5.2), expand = c(0, 0)) +
  labs(x = "Chromosome 16 position (Mb, GRCh38)", y = NULL) + theme_p31r() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = margin(2, 12, 4, 4))

# header (own strip, 7 pt): which method gave which result (S2-D25: the reference-LD fit is the inadmissible one, not the single-variant model)
PP72 <- led("L072", NULL, fmtd(3), where = "S4 header: DPEP1-cysteinylglycine coloc.abf PP.H4", anchor = "coloc\\.abf, PP\\.H4 {} at the declared prior")
ph <- ggplot() +
  annotate("text", x = 0, y = 2, hjust = 0, vjust = 1, size = ANNOT, family = BASE_FONT, lineheight = 1.15, fontface = "plain",
           label = paste0("DPEP1", EN, "cysteinylglycine colocalized under the single-causal-variant model (coloc.abf, PP.H4 ", PP72, " at the declared prior).\n",
                          "The reference-LD fine-mapping fit is inadmissible and is not cited as support. Attribution between DPEP1 and CHMP1A is unresolved.")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2.2), expand = FALSE) + theme_void()
row3 <- (pd1 / pd2 / pd3 / pdg) + plot_layout(heights = c(1, 1, 1, 0.95))
P <- (wrap_elements(full = ph) / wrap_elements(full = row3)) + plot_layout(heights = c(8, 92))
af3_save(P, h_mm = 104, panels = character(0), panel_map = character(0),
         frames = list(tracks = TRK[, .(trk, pos, y, rsid, ld)], marked_variants = MK[, .(rsid, role, lab, pos38)],
                       genes = gn[, .(g, strand, x0, x1, t0, t1)]))
cat(sprintf("S4: window %s-%s; markers %s; LD bins in window: %s; r(A,B) %.4f; window variants %d, with eGFR %d; top shared variant r2 A %.3f B %.3f\n",
            WIN[1], WIN[2], paste(MK$rsid, MK$lab, MK$pos38, collapse = "; "),
            paste(sprintf("%s %d", LDLEV, sapply(LDLEV, function(z) jw[ld == z, .N])), collapse = ", "), cor(X[, iA], X[, iB]),
            nrow(jw), jw[!is.na(yE), .N], jw[rsid == TOP, rA^2], jw[rsid == TOP, rB^2]))
