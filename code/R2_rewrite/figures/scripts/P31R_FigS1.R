# Paper 3.1 rewrite (R2) -- Additional file 3, Figure S1: post hoc external extension, all 28 T2D and CAD pairs (moved here from the old
# Figure 2c; decision F3, S2-D24). Illustrates SM section S-R11. ONE panel, no tag.
#   x  discovery estimate (FinnGen), y  external-check estimate (DIAMANTE for T2D, CARDIoGRAMplusC4D for CAD); log-odds per SD of protein,
#   asinh scale on both axes (one scale, so y = x stays the diagonal). Fill = evidence tier (PAL_TIER); filled = external p < 0.05 (the resource's
#   own p), open = p >= 0.05; diamonds = the five pairs selected in advance; the opposite-sign quadrants are shaded.
# Sources: ONLY through the approved manifest (AF3_S1): the post hoc results table and the tiered pair table (copy B, reconciled 2026-10-16).
# The Stage 3a-iii generator also cross-checked the discovery estimates against two R1 outcome-check files; both are SUPERSEDED in the manifest, so
# those two cross-checks are not repeated here (the discovery estimate is checked against the tiered table's mr_beta, as before).
# Run:  Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_FigS1.R"

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
suppressPackageStartupMessages(library(ggtext))
UNIT <- "FigureS1"
af3_begin(UNIT)
EN <- "\u2013"
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))

# ---- sources ---------------------------------------------------------------------------------------------------------------------
ph <- af3_read("R2_rewrite/posthoc/P31R_posthoc_results.tsv", "", sep = "\t")
pe <- af3_read("v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv", "")

# ---- invariants: what is plotted must reproduce the ledger ---------------------------------------------------------------------------
stopifnot(nrow(pe) == num("L010"), !anyDuplicated(pe[, .(gene, outcome)]))
stopifnot(nrow(ph) == num("L254", 2), nrow(ph) == num("L250", 1), !anyDuplicated(ph$gene), is.logical(ph$predeclared),
          setequal(paste(ph$gene, ph$outcome), paste(pe[outcome %in% c("T2D", "CAD")]$gene, pe[outcome %in% c("T2D", "CAD")]$outcome)))
phm <- merge(ph, pe[, .(gene, outcome, mr_beta, evidence_tier)], by = c("gene", "outcome"))
stopifnot(nrow(phm) == nrow(ph), max(abs(phm$disc_b - phm$mr_beta)) < 1e-5, all(phm$tier == phm$evidence_tier))
P05 <- num("L961", 1); PHET <- num("L961", 4)   # L961 (S2-D25): [p < 0.05, 6 checks, 0.05/6, heterogeneity p < 0.05]
stopifnot(P05 == 0.05, PHET == 0.05)
ph[, same := sign(wald) == sign(disc_b)]
nps <- ph[predeclared == FALSE]
stopifnot(nrow(nps) == num("L251", 1), sum(nps$same) == num("L251", 2), sum(nps$same & nps$p < P05) == num("L251", 3),
          sum(ph$same) == num("L250", 2), sum(ph$same & ph$p < P05) == num("L250", 3))
opp <- ph[same == FALSE]
L253 <- p31r_ledger("L253", 1)
stopifnot(nrow(opp) == 1, opp$gene == L253[[1]], opp$outcome == L253[[2]], opp$tier == L253[[3]], abs(signif(opp$p, 2) - L253[[4]]) < 1e-12)
het <- ph[het_p < PHET]
stopifnot(nrow(het) == num("L254", 1), identical(sort(het[abs(wald) < abs(disc_b)]$gene), unlist(p31r_ledger("L258", 3))),
          identical(sort(het[abs(wald) >= abs(disc_b)]$gene), unlist(p31r_ledger("L258", 4))))
stopifnot(all(ph[same & p < P05]$cls == "concordant"), all(ph[same & p >= P05]$cls == "directionally concordant"),
          all(ph[!same & p >= P05]$cls == "inconclusive"), identical(ph$p < P05, ph$cls == "concordant"))
stopifnot(all(ph[predeclared == TRUE]$p < P05), sum(ph$predeclared) == 5L)

# ---- encodings ------------------------------------------------------------------------------------------------------------------
ph[, kind := factor(fifelse(predeclared == TRUE, "sel", fifelse(p < P05, "sig", "ns")), levels = c("sig", "ns", "sel"))]
ph[, `:=`(fillc = fifelse(p < P05, unname(PAL_TIER[tier]), "white"),
          colc = fifelse(kind == "sel", "white", fifelse(p < P05, "grey15", unname(PAL_TIER[tier]))))]
stopifnot(!anyNA(ph$fillc), !anyNA(ph$colc), all(ph$tier %in% names(PAL_TIER)))
SIG <- 0.1
Tf <- function(x) asinh(x / SIG); Ti <- function(t) SIG * sinh(t)
tr <- scales::new_transform("asinh01", Tf, Ti)
CBX <- c(-3, -1, -0.3, 0, 0.3, 1); CBY <- CBX
LIMX_T <- c(-5.05, 3.96); LIMY_T <- c(-4.55, 3.75)
stopifnot(all(Tf(CBX) > LIMX_T[1] & Tf(CBX) < LIMX_T[2]), all(Tf(CBY) > LIMY_T[1] & Tf(CBY) < LIMY_T[2]),
          min(Tf(ph$disc_b)) > LIMX_T[1] + 0.2, max(Tf(ph$disc_b)) < LIMX_T[2] - 0.2,
          min(Tf(ph$wald)) > LIMY_T[1] + 0.2, max(Tf(ph$wald)) < LIMY_T[2] - 0.2)
CBX_LAB <- ticks("S1_x_discovery", CBX, Ti(LIMX_T), "asinh", "S1 x, log-odds per SD of protein, asinh scale")
CBY_LAB <- ticks("S1_y_external", CBY, Ti(LIMY_T), "asinh", "S1 y, log-odds per SD of protein, asinh scale")
# disclosed: KCTD5-T2D and PRSS27-T2D share one instrument; FN1-CAD and GRP-T2D coincide
shared <- ph[, .N, by = instrument][N > 1]
stopifnot(nrow(shared) == 1, setequal(ph[instrument %in% shared$instrument]$gene, c("KCTD5", "PRSS27")))
dd <- as.matrix(dist(ph[, .(Tf(disc_b), Tf(wald))])); diag(dd) <- Inf
stopifnot(setequal(ph$gene[unique(as.vector(which(dd < 0.1, arr.ind = TRUE)))], c("FN1", "GRP")),
          ph[gene == "FN1"]$kind == ph[gene == "GRP"]$kind, ph[gene == "FN1"]$tier == ph[gene == "GRP"]$tier)
# labels: the five pairs selected in advance and every pair whose external estimate differed in size, placed by hand in asinh units
# (tx, ty) = label centre; (lx, ly) = the point its leader reaches (NA: label level with its marker, no leader). No repel.
LAB_C <- data.table(gene = c("SLC16A1", "ARG1", "RAD23B", "PAM", "CELSR2", "TGFB1", "PCSK9", "NCAN", "LY96"),
                    tx = c(-3.89, -2.05, -3.47, -0.795, -1.28, 1.05, 2.959, 2.97, 3.04),
                    ty = c(-0.72, -0.28, -1.70, -2.12, -3.06, 2.60, 1.24, 0.44, -0.30),
                    lx = c(-2.819, -1.40, NA, -0.81, -2.26, 0.29, 2.042, 2.12, NA),
                    ly = c(-0.72, -0.46, NA, -1.94, -2.86, 2.423, 1.42, 0.30, NA))
stopifnot(setequal(LAB_C$gene, ph[predeclared == TRUE | het_p < PHET]$gene),
          all(LAB_C$tx > LIMX_T[1] + 1 & LAB_C$tx < LIMX_T[2] - 0.5), all(LAB_C$ty > LIMY_T[1] + 0.5 & LAB_C$ty < LIMY_T[2] - 0.5))
cl <- merge(LAB_C, ph[, .(gene, outcome, disc_b, wald, kind)], by = "gene")
cl[, txt := paste0(gene, EN, outcome)]
cl[, `:=`(px = Tf(disc_b), py = Tf(wald))]
R_EXT <- c(sig = 0.103, ns = 0.135, sel = 0.155)
ld <- cl[!is.na(lx)]
ld[, d := sqrt((lx - px)^2 + (ly - py)^2)][, `:=`(ux = (lx - px) / d, uy = (ly - py) / d)]
ld[, r0 := fifelse(kind == "sel", R_EXT[["sel"]] / (abs(ux) + abs(uy)), R_EXT[as.character(kind)]) + 0.025]
ld[, `:=`(sx = px + ux * r0, sy = py + uy * r0, ex = lx - ux * 0.05, ey = ly - uy * 0.05)]
stopifnot(!anyNA(ld$r0), all(ld$d - ld$r0 - 0.05 >= 0.15))
oth <- ph[, .(gene, ox = Tf(disc_b), oy = Tf(wald), orad = R_EXT[as.character(kind)])]
for (i in seq_len(nrow(ld))) { o <- oth[gene != ld$gene[i]]; stopifnot(min(sqrt((o$ox - ld$sx[i])^2 + (o$oy - ld$sy[i])^2) - o$orad) >= 0.10) }

quad <- data.table(xmin = c(-Inf, 0), xmax = c(0, Inf), ymin = c(0, -Inf), ymax = c(Inf, 0))
PTS <- function(d, stroke, size = 1.8) geom_point(data = d, aes(shape = kind, fill = fillc, colour = colc), size = size, stroke = stroke)
SEL_ORDER <- c("ARG1", "RAD23B", "PAM", "TGFB1", "PCSK9")
stopifnot(setequal(SEL_ORDER, ph[kind == "sel"]$gene))
P_TXT <- led("L961", 1, function(v) formatC(v, format = "f", digits = 2), where = "S1 key: external p < / >= threshold (drawn twice)",
             anchor = "^External p [<\u2265]\\s?{}$", count = 2L)

p1 <- ggplot(ph, aes(disc_b, wald)) +
  geom_rect(data = quad, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = SHADE_OPPOSITE, colour = NA) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.35, colour = "grey45") +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35, colour = "grey45") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.35, colour = "grey45") +
  geom_segment(data = ld, aes(x = Ti(sx), y = Ti(sy), xend = Ti(ex), yend = Ti(ey)), inherit.aes = FALSE, linewidth = 0.25, colour = "grey30") +
  PTS(ph[kind == "ns"], 0.8) + PTS(ph[kind == "sig"], 0.3) +
  PTS(ph[kind == "sel"][match(SEL_ORDER, gene)], 0.5, size = 2.3) +
  geom_text(data = cl, aes(x = Ti(tx), y = Ti(ty), label = txt), inherit.aes = FALSE, hjust = 0.5, vjust = 0.5, size = ANNOT, family = BASE_FONT) +
  annotate("text", x = Ti(LIMX_T[1] + 0.15), y = Ti(LIMY_T[2] - 0.15), label = "asinh scale", hjust = 0, vjust = 1, size = ANNOT, family = BASE_FONT, colour = "grey25") +
  scale_shape_manual(values = c(sig = 21, ns = 21, sel = 23), drop = FALSE, name = NULL,
                     labels = c(sig = paste0("External <i>p</i> < ", P_TXT), ns = paste0("External <i>p</i> \u2265 ", P_TXT), sel = "Original external check")) +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_continuous(trans = tr, breaks = CBX, labels = CBX_LAB) +
  scale_y_continuous(trans = tr, breaks = CBY, labels = CBY_LAB) +
  coord_cartesian(xlim = Ti(LIMX_T), ylim = Ti(LIMY_T), expand = FALSE) +
  labs(x = "Discovery (FinnGen)\nlog-odds per SD of protein", y = "External check\nlog-odds per SD of protein") +
  theme_p31r() +
  theme(aspect.ratio = diff(LIMY_T) / diff(LIMX_T), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0), legend.key.spacing.y = unit(0, "pt"), legend.key.spacing.x = unit(8, "pt"),
        legend.text = ggtext::element_markdown(size = BODY_PT, family = BASE_FONT)) +
  guides(shape = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(fill = c(unname(PAL_NEUTRAL[["dark"]]), "white", unname(PAL_NEUTRAL[["dark"]])),
                                                  colour = c("grey15", unname(PAL_NEUTRAL[["dark"]]), "white"), stroke = c(0.3, 0.8, 0.5), size = c(1.8, 1.8, 2.3))))

# tier key (own strip: the fill of the markers carries the evidence tier, as in Figure 2a; a tier hue means the tier and nothing else)
tk <- data.table(tier = factor(names(PAL_TIER), levels = names(PAL_TIER)), x = c(0.00, 0.17, 0.34))
pk <- ggplot(tk) +
  geom_point(aes(x = x, y = 1, fill = tier), shape = 21, size = 1.8, stroke = 0.3, colour = "grey15") +
  geom_text(aes(x = x + 0.03, y = 1, label = tier), hjust = 0, size = ANNOT, family = BASE_FONT) +
  annotate("text", x = -0.02, y = 1, label = "Evidence tier (as in Figure 2a)", hjust = 1, size = ANNOT, family = BASE_FONT) +
  scale_fill_manual(values = PAL_TIER, guide = "none") +
  coord_cartesian(xlim = c(-0.62, 0.62), ylim = c(0.5, 1.5), expand = FALSE) + theme_void()

fig <- p1 / pk + plot_layout(heights = c(1, 0.06))
af3_save(fig, h_mm = 150, panels = character(0), panel_map = character(0),
         frames = list(posthoc_28 = ph[, .(gene, outcome, tier, instrument, predeclared, disc_b, disc_se, wald, se, p, het_p, same, kind, cls)]))
