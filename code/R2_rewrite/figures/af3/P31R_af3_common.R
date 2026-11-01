# Paper 3.1 rewrite (R2) -- shared helpers for the AF3 generators (and any generator that follows the S2-D24 conventions).
# source() AFTER P31R_fig_common.R, at top level:
#     source(".../figures/scripts/P31R_fig_common.R"); source(".../figures/af3/P31R_af3_common.R")
#
# What it adds to the main-figure conventions (it does NOT change P31R_fig_common.R, whose md5 every main render log pins):
#   * this file's md5 is appended to P31R_DEPS, so the render log pins it in `deps` beside the theme and palette (gate G8 re-derives it);
#   * af3_begin(unit): refuses a generator that is not the ACTIVE generator of that unit in af3/P31R_ACTIVE_FIGURES.json;
#   * af3_source(source_file, panel) / af3_read(): the ONLY way to open plotted data. The file must have an APPROVED row for this figure (and panel)
#     in registry/APPROVED_SOURCE_MANIFEST.tsv at the manifest's md5 (fail closed: a superseded, blocked or edited file stops the render). The
#     declared sources are written to the registry for gate G12, which re-checks them and the render-log pins against the manifest;
#   * led() / ticks(): the ledger-registration and axis-tick helpers every main generator carries a private copy of, defined ONCE here
#     (same semantics as spec section 3: register at the point of printing; idx 1-based for lists, a key for objects);
#   * af3_save(): p31r_save() (cairo PDF, PNG rasterised from the PDF, render-log row with generator + common + deps + data md5s) PLUS the
#     registry JSON and the plotted frames, written at save time. 170 mm wide, at most 225 mm high PER PAGE. A figure that needs more room
#     continues on unit <id>_p2, <id>_p3 (set `pages` in the active registry); it is never shrunk.
# No G1b for AF3 (decision F9); G1, G2 (7/9/11 pt), G11 and the ledger-number gate G7 apply.

stopifnot(exists("p31r_save"), exists("P31R_FIGS"), exists("P31R_DEPS"), exists(".p31r_this"))
P31R_AF3_COMMON <- .p31r_this()
P31R_AF3 <- file.path(P31R_FIGS, "af3")
P31R_DEPS <- c(P31R_DEPS, P31R_AF3_COMMON)          # pinned in the render log (`deps`); p31r_save() reads P31R_DEPS when it is called

.AF3 <- new.env()
REG <- list(ledger = list(), literals = list(), sources = list())
af3_reset <- function() { REG <<- list(ledger = list(), literals = list(), sources = list()); invisible(NULL) }

# ---- the active registry ----------------------------------------------------------------------------------------
af3_active <- function() {
  p <- file.path(P31R_AF3, "P31R_ACTIVE_FIGURES.json")
  if (!file.exists(p)) stop("active figure registry missing: ", p)
  jsonlite::fromJSON(p, simplifyVector = FALSE)
}
# af3_begin("FigureS5") or af3_begin("FigureS5_p2"): the unit this run renders. Call once, first.
af3_begin <- function(unit) {
  stopifnot(is.character(unit), length(unit) == 1, grepl("^Figure(_[0-9]+|S[0-9]+)(_p[0-9]+)?$", unit))
  fid <- sub("_p[0-9]+$", "", unit)
  fe <- Filter(function(f) identical(f$id, fid), af3_active()$figures)
  if (length(fe) != 1) stop(unit, ": not an active figure (see af3/P31R_ACTIVE_FIGURES.json)")
  fe <- fe[[1]]
  page <- if (grepl("_p[0-9]+$", unit)) as.integer(sub(".*_p", "", unit)) else 1L
  if (page > fe$pages) stop(unit, ": the active registry gives ", fe$id, " ", fe$pages, " page(s)")
  gen <- sub(paste0("^", V3_ROOT, "/"), "", .p31r_generator())
  if (!identical(fe$generator, gen)) stop("this script (", gen, ") is not the active generator of ", fid, " (", fe$generator, "): a retired generator cannot render")
  .AF3$unit <- unit; .AF3$figure <- fid; .AF3$entry <- fe
  af3_reset()
  invisible(fe)
}

# ---- sources: only through the approved manifest -----------------------------------------------------------------
# registry/APPROVED_SOURCE_MANIFEST.tsv (S2-D24 task B): figure, panel, series, source_file, md5, columns_or_ids, status, disposition, kind, notes.
.MANIFEST_COLS <- c("figure", "panel", "series", "source_file", "md5", "columns_or_ids", "status", "disposition", "kind", "notes")
af3_manifest <- function() {
  p <- Sys.getenv("P31R_SOURCE_MANIFEST", file.path(P31R_FIGS, "registry", "APPROVED_SOURCE_MANIFEST.tsv"))
  if (!file.exists(p)) stop("approved source manifest missing: ", p)
  m <- utils::read.delim(p, sep = "\t", colClasses = "character", quote = "", na.strings = character(0), check.names = FALSE, encoding = "UTF-8")
  if (!identical(names(m), .MANIFEST_COLS)) stop("manifest columns are ", paste(names(m), collapse = ","), ", want ", paste(.MANIFEST_COLS, collapse = ","))
  m
}
.man_names <- function(fid) if (grepl("^FigureS[0-9]+$", fid)) c(fid, sub("^FigureS", "AF3_S", fid), sub("^FigureS", "Figure_S", fid)) else fid
.abs_path <- function(p) normalizePath(if (grepl("^([A-Za-z]:|/)", p)) p else file.path(V3_ROOT, p), winslash = "/", mustWork = FALSE)
.panel_ok <- function(row_panel, want) { if (is.null(want) || !nzchar(want)) return(rep(TRUE, length(row_panel)))
  rp <- vapply(strsplit(trimws(row_panel), "\\s+"), function(x) if (length(x)) x[1] else "", ""); rp %in% c(want, "-", "flow") }
# the path of an APPROVED source for this figure (and panel), after every check; recorded in the registry's `sources` (gate G12)
af3_source <- function(source_file, panel = NULL) {
  if (is.null(.AF3$unit)) stop("call af3_begin(unit) first")
  m <- af3_manifest(); want <- .abs_path(source_file)
  same <- vapply(m$source_file, function(x) identical(.abs_path(x), want), TRUE)
  mine <- m$figure %in% .man_names(.AF3$figure)
  have <- { q_ <- .p31r_resolve(want); sm_ <- .p31r_declared_md5(want); if (is.null(sm_)) .md5(q_) else sm_ }   # released copy (verbatim: its md5; lifted: the md5 of its source), or "input not in the release"
  sup <- same & m$disposition == "SUPERSEDED" & m$md5 == have
  if (any(sup)) stop("source ", source_file, " is a SUPERSEDED version (", m$notes[which(sup)[1]], "): not drawn from")
  appr <- same & mine & m$disposition == "APPROVED"
  if (!any(appr)) {
    if (any(same & mine & m$disposition == "BLOCKED")) stop("source ", source_file, " is BLOCKED for ", .AF3$figure, ": ", m$notes[which(same & mine & m$disposition == "BLOCKED")[1]])
    stop("source ", source_file, " has no APPROVED manifest row for ", .AF3$figure)
  }
  if (!any(appr & .panel_ok(m$panel, panel))) stop("source ", source_file, " is approved for ", .AF3$figure, " but not for panel '", panel, "'")
  ok <- appr & .panel_ok(m$panel, panel)
  if (!any(ok & (m$md5 == have | m$kind == "ledger"))) stop("source ", source_file, " is at md5 ", have, "; the manifest approves ", paste(unique(m$md5[ok]), collapse = ", "))
  REG$sources[[length(REG$sources) + 1]] <<- list(path = sub(paste0("^", V3_ROOT, "/"), "", want), panel = if (is.null(panel)) "" else panel)
  p31r_use(want)
}
af3_read <- function(source_file, panel = NULL, ...) p31r_read(af3_source(source_file, panel), ...)

# ---- registration at the point of printing (spec section 3) -----------------------------------------------------
# led(): fetch ONE ledger element, print it with fmt, record {text, id, idx, count, anchor, where} from THIS call. anchor: a regex with {}
# for the number, matched against the printed line ("^{}$" = a line holding only the number). count = times the label is DRAWN.
led <- function(id, idx = NULL, fmt = function(v) as.character(v), where, anchor, count = 1L) {
  v <- p31r_ledger(id, idx)
  stopifnot(length(v) == 1, is.numeric(v), is.character(where), is.character(anchor))
  text <- fmt(v)
  REG$ledger[[length(REG$ledger) + 1]] <<- c(list(text = text, id = id), if (!is.null(idx)) list(idx = idx),
                                              list(count = as.integer(count), anchor = anchor, where = where))
  text
}
# ticks(): the labels a scale draws, from its own break vector (p31r_lab: U+2212 minus), one registry entry per printed break
ticks <- function(axis, breaks, lims, transform, why) {
  stopifnot(all(breaks > lims[1] & breaks < lims[2]), !anyDuplicated(breaks), transform %in% c("linear", "log10", "sqrt", "asinh", "custom"))
  labs <- p31r_lab(breaks)
  for (t in labs) REG$literals[[length(REG$literals) + 1]] <<- list(text = t, why = paste0("axis_tick: ", why), axis = axis,
                                                                    transform = transform, count = 1L, anchor = "^{}$")
  labs
}

# the registry as the gate reads it (figgate g7_numbers / G0 panel tags / G12 sources); a function of its own so the round trip is testable
af3_registry_json <- function(panels, panel_map = character(0), frame_names = character(0)) {
  jsonlite::toJSON(list(panels = I(panels), ledger = REG$ledger, derived = list(), literals = REG$literals, sources = REG$sources,
                        panel_map = as.list(panel_map), frames = I(frame_names)), auto_unbox = TRUE, pretty = TRUE)
}

# ---- save: figure + registry + plotted frames ---------------------------------------------------------------------
# p: the patchwork/ggplot; h_mm <= 225 for this PAGE; panels: the panel letters this page carries (character(0) for none);
# panel_map: named character vector letter -> what it shows; frames: named list of data.tables = what was plotted (written to
# af3/source_data/<unit>__<name>.csv so a number on the page can be traced to a row).
af3_save <- function(p, h_mm, panels, panel_map = character(0), frames = list()) {
  if (is.null(.AF3$unit)) stop("call af3_begin(unit) first")
  stopifnot(is.character(panels), setequal(names(panel_map), panels), length(REG$sources) >= 1, is.list(frames), length(frames) >= 1,
            all(vapply(frames, function(d) is.data.frame(d) && nrow(d) > 0, TRUE)), !is.null(names(frames)), !anyDuplicated(names(frames)))
  unit <- .AF3$unit
  dir.create(file.path(P31R_AF3, "source_data"), showWarnings = FALSE)
  for (nm in names(frames)) data.table::fwrite(frames[[nm]], file.path(P31R_AF3, "source_data", sprintf("%s__%s.csv", unit, nm)))
  p31r_save(p, unit, w_mm = FIG_W_MM, h_mm = h_mm)
  regf <- file.path(P31R_FIGS, "registry", paste0(unit, ".json"))
  dir.create(dirname(regf), showWarnings = FALSE)
  writeLines(af3_registry_json(panels, panel_map, names(frames)), regf, useBytes = TRUE)
  message("registry: ", length(REG$ledger), " ledger entries, ", length(REG$literals), " literals, ", length(REG$sources), " source(s) -> ", regf)
  af3_reset()
  invisible(regf)
}
