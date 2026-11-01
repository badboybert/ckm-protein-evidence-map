# Paper 3.1 rewrite (R2) -- Additional file 3, Figure S3: kidney benchmark retention and the cystatin-C and BUN concordance of the renal proteins
# (moved here from the old Figure 3c; S2-D24). Illustrates SM section S-R12.
#   a  coloc.abf retention of the cystatin-C-confirmed kidney benchmark (L070); 95% Wilson interval on the share (plotted, not printed)
#   b  creatinine eGFR* against cystatin-C eGFR* per SD of protein for the 17 colocalization-supported kidney proteins (asinh scales)
#   c  creatinine eGFR* against BUN (natural-log BUN) per SD of protein (16 proteins with a BUN estimate)
#   Filled = p < 0.05 on the vertical marker; shaded = opposite signs; dashed lines = zero (and y = x in b).
# Sources: ONLY through the approved manifest (AF3_S3): the eGFR mirage classification (benchmark, 52 concordant of 89) and the full-precision
# kidney cystatin-C / BUN concordance table. The Stage 3a-iii generator also cross-checked these against the coloc, MR and instrument tables and
# the V4 scrubbed supplementary table; none of those is approved for AF3_S3, so those cross-checks are not repeated here (retention, the 12 of 17
# and the direction counts are asserted against the ledger instead).
# Cystatin-C p values: the full-precision file carries every p as a number (no placeholder is plotted or printed).
# Run:  Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_FigS3.R"

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
suppressPackageStartupMessages({ library(ggrepel) })
UNIT <- "FigureS3"
af3_begin(UNIT)
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))
INK <- "grey15"; LINE_GREY <- "grey45"; COL_BEN <- PAL_NEUTRAL[["mid"]]
TXT <- function(...) geom_text(..., size = ANNOT, family = BASE_FONT)
# ticks with some labels left blank (an asinh scale packs close breaks; the blank ones are not printed, so they are not registered)
ticks_q <- function(axis, breaks, lims, why, quiet = numeric(0)) {
  stopifnot(all(breaks > lims[1] & breaks < lims[2]), !anyDuplicated(breaks))
  labs <- p31r_lab(breaks); labs[breaks %in% quiet] <- ""
  for (t in labs[nzchar(labs)]) REG$literals[[length(REG$literals) + 1]] <<- list(text = t, why = paste0("axis_tick: ", why), axis = axis,
                                                                                   transform = "asinh", count = 1L, anchor = "^{}$")
  labs
}
wilson <- function(k, n, z = qnorm(0.975)) {
  p <- k / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d; hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = 100 * (c0 - hw), hi = 100 * (c0 + hw))
}

ben <- af3_read("v3.2/_engine_root/results/repro_phase1_hg38/pam_conditional/egfr_mirage_classification.csv", "a")
kc  <- af3_read("paper3_ckd_rerun_v3/out/D1/kidney_cystatin_bun_concordance_E_current_20260922.csv", "b")

# ---- a: benchmark retention ----------------------------------------------------------------------------------------------------------
P_NOM <- num("L260", 4)
stopifnot(nrow(ben) == num("L260", 1), !anyDuplicated(ben$gene), all(ben$crea_p < num("L260", 3)))
conc <- ben[class == "concordant"]
stopifnot(nrow(conc) == num("L260", 2), nrow(conc) == num("L070", 2), all(conc$cys_p < P_NOM), all(sign(conc$cys_b) == sign(conc$crea_b)))

# ---- b, c: the renal proteins ---------------------------------------------------------------------------------------------------------
PC <- c("SNP", "egfrcrea_b", "egfrcrea_p", "cys_b", "cys_p", "bun_b", "bun_p")
stopifnot(kc[, .(n = uniqueN(.SD)), by = gene, .SDcols = PC][, all(n == 1L)])        # one instrument, one estimate set per protein
k <- unique(kc, by = "gene")
stopifnot(nrow(k) == num("L208", 2), nrow(k) == num("L245", 2), !anyNA(k$cys_p), all(k$cys_p > 0))
# eGFR* = -ln(eGFR); BUN on its source GWAS scale, ln(BUN) (positive = higher BUN = lower kidney function)
k[, `:=`(x = -egfrcrea_b, y_cys = -cys_b, y_bun = bun_b)]
k[, `:=`(agree_cys = sign(x) == sign(y_cys), agree_bun = sign(x) == sign(y_bun))]
stopifnot(sum(k$agree_cys) == num("L208", 1), sum(k$agree_cys & k$cys_p < P_NOM) == num("L245", 1),
          identical(k[agree_cys == FALSE, gene], p31r_ledger("L208", 5)), k[gene == "DNER", cys_b] == 0,
          sum(!is.na(k$bun_b)) == num("L208", 4), sum(k$agree_bun, na.rm = TRUE) == num("L208", 3), identical(k[is.na(bun_b), gene], "PLAU"))
# the retained 12 of the benchmark are the 12 renal proteins concordant and nominally significant on cystatin C (L070 = L245)
stopifnot(num("L070", 1) == num("L245", 1), all(k[agree_cys & cys_p < P_NOM, gene] %in% conc$gene))
B0 <- data.table(x = 1, k = num("L070", 1), n = num("L070", 2))
B0[, `:=`(pct = 100 * k / n, lo = wilson(k, n)$lo, hi = wilson(k, n)$hi)]
B0[, lab := paste(led("L070", 1, function(v) formatC(v, format = "f", digits = 0), where = "S3a: retained", anchor = "^{} of 52$"), "of",
                  led("L070", 2, function(v) formatC(v, format = "f", digits = 0), where = "S3a: benchmark size", anchor = "^12 of {}$"))]
Y_A <- c(0, 50, 100)
pa <- ggplot(B0, aes(x, pct)) + geom_col(width = 0.6, fill = COL_BEN) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.35, colour = INK) +
  TXT(aes(y = hi, label = lab), vjust = -0.7) +
  scale_x_continuous(breaks = 1, labels = "Kidney\nbenchmark", limits = c(0.4, 1.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = Y_A, labels = ticks("S3_a_y", Y_A, c(-1, 101), "linear", "S3 a, retained by coloc.abf (%)"),
                     expand = expansion(mult = c(0, 0.14)), guide = guide_axis(cap = "upper")) +
  labs(x = NULL, y = "Retained by\ncoloc.abf (%)") + theme_p31r() +
  theme(axis.text.x = element_text(lineheight = 0.9), axis.title.y = element_text(lineheight = 0.9))

P_TXT <- led("L270", NULL, function(v) formatC(v, format = "f", digits = 2), where = "S3 key: p threshold (drawn twice)", anchor = "^p [<\u2265]\\s?{}$", count = 2L)
SIG_LV <- c(paste("p <", P_TXT), paste("p \u2265", P_TXT))
k[, `:=`(sig_cys = ifelse(cys_p < P_NOM, SIG_LV[1], SIG_LV[2]), sig_bun = ifelse(bun_p < P_NOM, SIG_LV[1], SIG_LV[2]))]
TR_C <- scales::pseudo_log_trans(sigma = 0.002, base = 10)                      # asinh(x / 0.004) / ln 10: an asinh scale
BR_C <- c(-0.03, -0.01, 0, 0.01, 0.03, 0.1, 0.3)
LIM_X <- c(-0.05, 0.05); stopifnot(all(abs(k$x) < 0.05))
LIM_YC <- c(-0.05, 0.3); stopifnot(all(k$y_cys > LIM_YC[1] & k$y_cys < LIM_YC[2]))
LIM_YB <- c(-0.05, 0.05); stopifnot(all(abs(k$y_bun) < 0.05, na.rm = TRUE))
quad <- data.table(xmin = c(-Inf, 0), xmax = c(0, Inf), ymin = c(0, -Inf), ymax = c(Inf, 0))
SIGFILL <- setNames(c(INK, "white"), SIG_LV)
X_TTL <- "Creatinine eGFR*\n(\u2212natural-log eGFR)\nper SD of protein"
pc_base <- function(dat, yv, sv, ylim, axis_id) {
  ybr <- BR_C[BR_C > ylim[1] & BR_C < ylim[2]]; xbr <- BR_C[BR_C > LIM_X[1] & BR_C < LIM_X[2]]
  ggplot(dat, aes(x, .data[[yv]])) +
    geom_rect(data = quad, inherit.aes = FALSE, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = SHADE_OPPOSITE) +
    geom_hline(yintercept = 0, colour = LINE_GREY, linetype = "dashed", linewidth = 0.3) + geom_vline(xintercept = 0, colour = LINE_GREY, linetype = "dashed", linewidth = 0.3) +
    geom_point(aes(fill = .data[[sv]]), shape = 21, size = 1.9, colour = INK, stroke = 0.45) +
    scale_x_continuous(trans = TR_C, limits = LIM_X, breaks = xbr, labels = ticks_q(paste0(axis_id, "_x"), xbr, LIM_X, paste0("S3 ", axis_id, " x, asinh scale"), quiet = c(-0.01, 0.01)),
                       expand = expansion(mult = 0.03)) +
    scale_y_continuous(trans = TR_C, limits = ylim, breaks = ybr, labels = ticks_q(paste0(axis_id, "_y"), ybr, ylim, paste0("S3 ", axis_id, " y, asinh scale")), expand = expansion(mult = 0.03)) +
    scale_fill_manual(values = SIGFILL, breaks = names(SIGFILL), name = NULL) +
    annotate("text", x = LIM_X[1], y = ylim[2], label = "asinh\nscales", hjust = 0, vjust = 1, size = ANNOT, family = BASE_FONT, lineheight = 0.95) +
    labs(x = X_TTL) + theme_p31r() +
    theme(legend.position = "none", axis.title.x = element_text(lineheight = 0.9), axis.title.y = element_text(lineheight = 0.9),
          plot.title = element_text(size = ANNOT_PT, face = "plain", hjust = 0, margin = margin(b = 2)), plot.title.position = "panel")
}
REPEL <- list(size = ANNOT, family = BASE_FONT, seed = 3, box.padding = 0.45, min.segment.length = 0, segment.size = 0.25, segment.colour = "grey40")
lab_c1 <- k[gene %in% c(p31r_ledger("L208", 5), "BRAP")]
pb <- pc_base(k, "y_cys", "sig_cys", LIM_YC, "b") +
  annotate("segment", x = LIM_X[1], xend = LIM_X[2], y = LIM_X[1], yend = LIM_X[2], linetype = "dashed", colour = LINE_GREY, linewidth = 0.3) +
  do.call(geom_text_repel, c(list(data = lab_c1, mapping = aes(label = gene), nudge_x = ifelse(lab_c1$gene == "DNER", 0.35, 0)), REPEL)) +
  labs(title = sprintf("%s of %s same direction",
                       led("L208", 1, function(v) formatC(v, format = "f", digits = 0), where = "S3b title: same direction", anchor = "^{} of 17 same direction$"),
                       led("L208", 2, function(v) formatC(v, format = "f", digits = 0), where = "S3b title: proteins", anchor = "^16 of {} same direction$")),
       y = "Cystatin-C eGFR*\n(\u2212natural-log eGFR)\nper SD of protein")
lab_c2 <- k[!is.na(bun_b) & agree_bun == FALSE]
stopifnot(nrow(lab_c2) == num("L208", 4) - num("L208", 3), "DNER" %in% lab_c2$gene)
pc <- pc_base(k[!is.na(bun_b)], "y_bun", "sig_bun", LIM_YB, "c") +
  do.call(geom_text_repel, c(list(data = lab_c2, mapping = aes(label = gene), nudge_x = ifelse(lab_c2$gene == "DNER", 0.35, 0)), REPEL)) +
  labs(title = sprintf("%s of %s same direction",
                       led("L208", 3, function(v) formatC(v, format = "f", digits = 0), where = "S3c title: same direction", anchor = "^{} of 16 same direction$"),
                       led("L208", 4, function(v) formatC(v, format = "f", digits = 0), where = "S3c title: proteins with BUN", anchor = "^14 of {} same direction$")),
       y = "BUN\n(natural-log BUN)\nper SD of protein") +
  theme(legend.position = "right", legend.key.size = unit(9, "pt"))

top <- (pa + labs(tag = "a")) + (pb + labs(tag = "b")) + (pc + labs(tag = "c")) + plot_layout(widths = c(0.6, 1.15, 1.45))
note <- ggplot() + annotate("text", x = 0, y = 1, hjust = 0, vjust = 0.5, size = ANNOT, family = BASE_FONT,
                            label = "Dashed lines mark zero; the dashed diagonal in panel b marks equal estimates on the two markers.") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2), expand = FALSE) + theme_void()
fig <- top / wrap_elements(full = note) + plot_layout(heights = c(75, 5))
af3_save(fig, h_mm = 80, panels = c("a", "b", "c"),
         panel_map = c(a = "coloc.abf retention of the kidney benchmark", b = "creatinine versus cystatin-C eGFR* per SD of protein, 17 renal proteins",
                       c = "creatinine eGFR* versus BUN per SD of protein, 16 renal proteins"),
         frames = list(renal_proteins = k[, .(gene, outcome, evidence_tier, SNP, egfrcrea_b, egfrcrea_p, cys_b, cys_p, bun_b, bun_p, x, y_cys, y_bun, agree_cys, agree_bun)],
                       benchmark = B0[, .(retained = k, benchmark = n, pct, wilson_lo = lo, wilson_hi = hi)]))
