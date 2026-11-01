# =============================================================================
# theme_publication_ckm.R  --  Paper 3.EAS project figure theme
# = the publication-figures-tables skill theme (theme_and_setup.R) with the
# USER RULE applied: HARMONIZED 2-TIER hierarchy at +1 pt (decided 2026-06-18).
# Body type (axis text = axis title = legend) is a single size; panel labels
# are larger + bold. This supersedes the earlier +2 pt 4-level spec: the 4-level
# steps (legend<axis-text, title only +1) read as inconsistency, not hierarchy,
# and +2 clipped dense panels. Dimensions + gridlines (off) UNCHANGED.
# Source this instead of the skill theme for CKM figures.
# =============================================================================
suppressMessages({ library(ggplot2); library(scales)
  if (requireNamespace("patchwork", quietly = TRUE)) try(library(patchwork), silent = TRUE) })

# --- fonts: harmonized 2-tier, +1 pt ----------------------------------------
# body = 9 pt everywhere in the plot field; panel label 11 pt bold. (skill
# default body = 8 pt; "+1" = 9 pt. Sparse panels also tolerate +2, but +1 is
# the locked anchor because it fits dense panels without clipping.)
BASE_FONT            <- "Arial"
BASE_FONT_SIZE       <- 9    # body: axis text
BASE_FONT_SIZE_TITLE <- 9    # axis title  = body (harmonized)
BASE_FONT_SIZE_LEGEND<- 9    # legend      = body (harmonized)
PANEL_LABEL_SIZE     <- 11   # panel label: body + 2, bold
GEOM_TEXT_MM         <- 3.2  # in-plot annotations ~9 pt (= body)

theme_publication_ckm <- function(base_size = BASE_FONT_SIZE,
                                   base_family = BASE_FONT,
                                   legend_position = "right") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text        = element_text(family = base_family, color = "black"),
      plot.title  = element_text(size = BASE_FONT_SIZE_TITLE, face = "plain", hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 1, color = "grey30", hjust = 0, margin = margin(b = 6)),
      axis.title  = element_text(size = BASE_FONT_SIZE_TITLE, color = "black"),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      axis.text   = element_text(size = base_size, color = "black"),
      axis.line   = element_line(linewidth = 0.4, color = "black"),
      axis.ticks  = element_line(linewidth = 0.4, color = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid  = element_blank(),         # gridlines off (unchanged)
      panel.border = element_blank(),
      legend.position = legend_position,
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key  = element_rect(fill = "transparent", color = NA),
      legend.text = element_text(size = BASE_FONT_SIZE_LEGEND, color = "black"),
      legend.title = element_text(size = BASE_FONT_SIZE_LEGEND, face = "plain", color = "black"),
      legend.key.size = unit(9, "pt"),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text  = element_text(size = BASE_FONT_SIZE_TITLE, face = "plain", margin = margin(t = 2, b = 2)),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

# dimensions UNCHANGED from the skill (89 single / 183 double mm)
DIM <- list(single_col_mm = 89, double_col_mm = 183, max_height_mm = 240)

save_figure_ckm <- function(plot, filename_base, width_mm = 183, height_mm = 120,
                            dpi = 300, output_dir = "../results/figures/") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggsave(file.path(output_dir, paste0(filename_base, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(output_dir, paste0(filename_base, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = dpi)
  message(sprintf("Saved %s.{pdf,png} (%g x %g mm)", filename_base, width_mm, height_mm))
}

set_publication_defaults_ckm <- function() {
  theme_set(theme_publication_ckm())
  update_geom_defaults("point",   list(size = 1.7))
  update_geom_defaults("text",    list(size = GEOM_TEXT_MM, family = BASE_FONT))
  update_geom_defaults("errorbar",list(linewidth = 0.4, width = 0.15))
}

# =============================================================================
# VARIANT "flat inside" (committed 2026-06-28; axes param added 2026-06-28) ---
# Uniform body font at the LOCKED size (9 pt) + inside bordered legend + boxed
# in-plot annotations, on the CLEAN base (gridlines off) — NOT the framed/gridlined
# full-flat (8 pt) variant. Callable for any CKM figure that wants the self-
# annotated, "inside-boxed" Nat-Commun-ish look while keeping the project's 9 pt.
#
# ⚠ AXIS-LINE RULE (do not hand-blank lines per panel — use the `axes` arg):
#   the base is theme_classic, which draws BOTH the x (bottom) and y (left) axis
#   lines (an L). On panels where an axis carries no real/labelled scale (heatmaps;
#   forests whose categories live in an inside legend instead of y-tick labels),
#   that bare line reads as a STRAY partial frame. `axes` blanks the unused line(s):
#     axes="both" (DEFAULT) — standard x/y plot: keep the L (e.g. regional Manhattan, scatter).
#     axes="x"    — keep bottom only      (e.g. forest with inside legend / no y labels: F4d, F4f).
#     axes="y"    — keep left only.
#     axes="none" — no axis lines/ticks   (e.g. tile heatmap: F3b; or pair with no-axis text for funnels).
#
#   usage:  p + theme_flat_inside(legend_position = c(0.80, 0.22), axes = "x")
#   boxed callouts (facet-aware): geom_box_label(data = df, aes(x, y, label), ...)
theme_flat_inside <- function(base_size = BASE_FONT_SIZE, base_family = BASE_FONT,
                              legend_position = c(0.80, 0.22),
                              axes = c("both", "x", "y", "none")) {
  axes <- match.arg(axes)
  th <- theme_publication_ckm(base_size = base_size, base_family = base_family) +
    theme(
      plot.subtitle    = element_text(size = base_size),   # force uniform (flat)
      strip.text       = element_text(size = base_size),
      legend.position  = legend_position,
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.margin    = margin(2, 4, 2, 4),
      legend.key       = element_rect(fill = "white", color = NA),
      legend.key.size  = unit(11, "pt")
    )
  th + switch(axes,
    both = theme(),
    x    = theme(axis.line.y = element_blank(), axis.ticks.y = element_blank()),
    y    = theme(axis.line.x = element_blank(), axis.ticks.x = element_blank()),
    none = theme(axis.line = element_blank(), axis.ticks = element_blank()))
}

# White-filled, thin-bordered text callout at the flat size (ggplot2 4.0: border
# is `linewidth`, not the deprecated `label.size`). Works inside facets (pass data).
geom_box_label <- function(..., size_pt = BASE_FONT_SIZE, color = "grey20") {
  ggplot2::geom_label(..., fill = "white", color = color, linewidth = 0.25,
                      size = size_pt / .pt, label.padding = unit(1.4, "pt"))
}
