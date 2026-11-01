# Paper 3.1 rewrite (R2) -- AF3 status glyph: the HATCHED marker (visual-status rule 3, "inadmissible fit"; plan v3 section 0, S2-D24).
# source() AFTER P31R_fig_common.R (and after P31R_af3_common.R). Its md5 is appended to P31R_DEPS, so the render log pins it in `deps`.
#
# Why a geom of its own: the five statuses need five distinct encodings, and the crossed circle (shape 13) already means "no colocalization" in
# Figure 4a(vi), AF3 S5a/b and AF3 S6c. An inadmissible fit is therefore drawn as a circle with three parallel diagonal chords, in every AF3 figure
# that has one (S2, S5 page 2). The marker is built in millimetres around each data point (grid units), so it does not depend on the axis scales.
#   geom_hatchdot(aes(x, y), data = ..., r_mm = 1.2, n = 3, ink = "grey15", linewidth = 0.3)

stopifnot(exists("P31R_DEPS"), exists(".p31r_this"))
P31R_DEPS <- c(P31R_DEPS, .p31r_this())

GeomHatchDot <- ggplot2::ggproto("GeomHatchDot", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(alpha = NA),
  draw_key = ggplot2::draw_key_point,
  draw_panel = function(data, panel_params, coord, r_mm = 1.2, n = 3, ink = "grey15", fill = "white", linewidth = 0.3) {
    cd <- coord$transform(data, panel_params)
    N <- nrow(cd); if (N == 0) return(grid::nullGrob())
    lwd <- linewidth * ggplot2::.pt
    circ <- grid::circleGrob(x = grid::unit(cd$x, "npc"), y = grid::unit(cd$y, "npc"), r = grid::unit(r_mm, "mm"),
                             gp = grid::gpar(fill = fill, col = ink, lwd = lwd))
    # chords at 45 degrees (rising to the right), offsets -(n-1)/2 ... (n-1)/2 times the spacing along the normal; endpoints on the circle
    sp <- 2 * r_mm / (n + 1); off <- (seq_len(n) - (n + 1) / 2) * sp
    h <- sqrt(r_mm^2 - off^2); u <- c(1, 1) / sqrt(2); v <- c(-1, 1) / sqrt(2)
    ex0 <- off * v[1] - h * u[1]; ey0 <- off * v[2] - h * u[2]; ex1 <- off * v[1] + h * u[1]; ey1 <- off * v[2] + h * u[2]
    ix <- rep(seq_len(N), each = n); ik <- rep(seq_len(n), times = N)
    seg <- grid::segmentsGrob(x0 = grid::unit(cd$x[ix], "npc") + grid::unit(ex0[ik], "mm"), y0 = grid::unit(cd$y[ix], "npc") + grid::unit(ey0[ik], "mm"),
                              x1 = grid::unit(cd$x[ix], "npc") + grid::unit(ex1[ik], "mm"), y1 = grid::unit(cd$y[ix], "npc") + grid::unit(ey1[ik], "mm"),
                              gp = grid::gpar(col = ink, lwd = lwd, lineend = "butt"))
    grid::gTree(children = grid::gList(circ, seg))
  })

geom_hatchdot <- function(mapping = NULL, data = NULL, r_mm = 1.2, n = 3, ink = "grey15", fill = "white", linewidth = 0.3, inherit.aes = TRUE, ...) {
  ggplot2::layer(geom = GeomHatchDot, mapping = mapping, data = data, stat = "identity", position = "identity", show.legend = FALSE,
                 inherit.aes = inherit.aes, params = list(r_mm = r_mm, n = n, ink = ink, fill = fill, linewidth = linewidth, ...))
}

# DOTTED RING (visual-status rule 5, "not robust to the reference LD"): an open circle drawn with a dotted line, radius in millimetres (scale-free).
GeomDotRing <- ggplot2::ggproto("GeomDotRing", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(alpha = NA),
  draw_key = ggplot2::draw_key_point,
  draw_panel = function(data, panel_params, coord, r_mm = 1.0, ink = "grey15", linewidth = 0.4, dash = "11") {
    cd <- coord$transform(data, panel_params)
    if (nrow(cd) == 0) return(grid::nullGrob())
    grid::circleGrob(x = grid::unit(cd$x, "npc"), y = grid::unit(cd$y, "npc"), r = grid::unit(r_mm, "mm"),
                     gp = grid::gpar(fill = NA, col = ink, lwd = linewidth * ggplot2::.pt, lty = dash))
  })

geom_dotring <- function(mapping = NULL, data = NULL, r_mm = 1.0, ink = "grey15", linewidth = 0.4, dash = "11", inherit.aes = TRUE, ...) {
  ggplot2::layer(geom = GeomDotRing, mapping = mapping, data = data, stat = "identity", position = "identity", show.legend = FALSE,
                 inherit.aes = inherit.aes, params = list(r_mm = r_mm, ink = ink, linewidth = linewidth, dash = dash, ...))
}
