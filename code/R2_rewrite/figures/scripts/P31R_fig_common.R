# Paper 3.1 rewrite (R2) -- shared figure helpers. source() this first from every figure generator.
#
# What it gives a generator:
#   * the LOCKED house theme (v3.2/scripts/theme_publication_ckm.R) and palette (ckm_palette.R), read-only;
#   * size rules: body 9 pt (axis text = axis title = legend), panel tags 11 pt bold, EVERY in-panel text layer 7 pt
#     (ANNOT, the USER rule of 2026-09-25); page 170 mm (or 85 mm) wide, at most 225 mm high (BMC);
#   * p31r_read(): a FAIL-CLOSED reader that records the file for the render log (fread(<missing path with a space>)
#     runs a shell command and returns an empty table under "CKM papers"; this refuses instead);
#   * p31r_ledger(): numbers printed on a figure come from the claims ledger, never typed;
#   * p31r_save(): cairo PDF (the shipped representation) + a PNG RASTERISED FROM THAT PDF (so the two cannot
#     disagree on glyphs or layout), and a render-log row pinning generator, this file, every data file and both outputs.
# The PDF gate (P31R_figgate.py) re-derives every pin from disk.

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(jsonlite); library(patchwork); library(scales)
})

.p31r_this <- function() {
  f <- sys.frames(); for (i in rev(seq_along(f))) { o <- f[[i]]$ofile; if (!is.null(o)) return(normalizePath(o, winslash = "/")) }
  stop("P31R_fig_common.R must be source()d")
}
P31R_COMMON <- .p31r_this()
P31R_FIGS   <- normalizePath(file.path(dirname(P31R_COMMON), ".."), winslash = "/")       # R2_rewrite/figures
P31R_ROOT   <- normalizePath(file.path(P31R_FIGS, ".."), winslash = "/")                  # R2_rewrite
V3_ROOT     <- normalizePath(file.path(P31R_ROOT, ".."), winslash = "/")                  # V3
.p31r_arch  <- normalizePath(file.path(V3_ROOT, ".."), winslash = "/", mustWork = FALSE)          # the archive root (V3_ROOT is <archive>/code)
P31R_LEDGER <- file.path(.p31r_arch, "ledger/claims_ledger_public.tsv")                           # the PUBLIC extract (masked rows carry their reason)
.V32 <- file.path(V3_ROOT, "v3.2/scripts")
P31R_DEPS <- c(file.path(.V32, "theme_publication_ckm.R"), file.path(.V32, "ckm_palette.R"))
for (.d in P31R_DEPS) { if (!file.exists(.d)) stop("missing house file: ", .d); source(.d) }

# ---- size rules ------------------------------------------------------------------------------------------------
BODY_PT  <- 9                     # axis text = axis title = legend (locked 9/11)
TAG_PT   <- 11                    # panel letters, bold
ANNOT_PT <- 7                     # every in-panel text layer
ANNOT    <- ANNOT_PT / ggplot2::.pt   # geom_text/geom_label size (mm)
FIG_W_MM  <- 170                  # BMC double column
FIG_W1_MM <- 85                   # BMC single column
FIG_HMAX_MM <- 225

theme_p31r <- function(base_size = BODY_PT, ...) {
  theme_publication_ckm(base_size = base_size, base_family = BASE_FONT, ...) +
    theme(plot.tag = element_text(size = TAG_PT, face = "bold", family = BASE_FONT),
          plot.title = element_blank(), plot.subtitle = element_blank())
}
# composite: tags a, b, c ... at 11 pt bold
p31r_tags <- function() plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(size = TAG_PT, face = "bold")))

# ---- palettes (ONE key per concept across Figures 1-5; same hexes as the v3.2 package) ---------------------------
OI <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442",
        blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7", grey = "#999999")   # Okabe-Ito
PAL_OUTCOME <- CKM_OUTCOME[c("CAD", "HF", "AF", "Stroke", "T2D", "NAFLD", "eGFR", "CKD")]
PAL_TIER    <- c(Convergent = CKM_DEFLATE[["anchor"]], Supported = CKM_DEFLATE[["colocalized"]], Flagged = CKM_DEFLATE[["tested"]])
PAL_DIR     <- c(risk = CKM_DIR[["risk"]], protective = CKM_DIR[["protective"]], null = CKM_DIR[["null"]])
PAL_ANC_COL <- CKM_ANC_COL; PAL_ANC_SHAPE <- CKM_ANC_SHAPE
# ONE MEANING PER HUE across Figures 1-5 (3a cross-figure critic X1-X3; spec figures/P31R_FIGURE_SPEC.md). Tier hues
# (green / blue / grey above) mean the evidence tier and nothing else; everything below avoids them.
PAL_SIGNAL       <- c(A = "#CC79A7", B = "#E69F00")    # the two signals at one locus (ARG1 in Figure 4, DPEP1 in Figure 3d)
PAL_SIGNAL_LIGHT <- c(A = "#E8B9D3", B = "#F3CF80")    # their weaker LD bin (r2 0.2-0.5)
PAL_GROUP        <- c(kidney = "#333333", other = "#BDBDBD")   # kidney vs other outcomes where tiers are not shown
PAL_NEUTRAL      <- c(dark = "#333333", mid = "#888888", light = "#D9D9D9", other = "#C8C8C8")
SHADE_OPPOSITE   <- "#EEEEEE"   # the two quadrants where paired estimates have OPPOSITE signs (Figures 2c, 3c, 5b, 5c)
# axis labels: U+2212 minus on every numeric scale (gate G7m); use as scale_x_continuous(labels = p31r_lab)
p31r_lab <- scales::label_number(style_negative = "minus", big.mark = ",", drop0trailing = TRUE)

# ---- fail-closed data access, recorded for the render log ------------------------------------------------------
.p31r_read_set <- new.env()
# an input that is not at its original place is looked up in the release map (code/P31R_INPUT_MAP.tsv): a verbatim released copy is used,
# anything else stops. The map is keyed by the path relative to V3_ROOT (= <archive>/code).
.p31r_map <- NULL
.p31r_resolve <- function(p) {
  if (file.exists(p)) return(p)
  pre <- paste0(V3_ROOT, "/")
  key <- if (startsWith(p, pre)) substr(p, nchar(pre) + 1, nchar(p)) else p
  if (is.null(.p31r_map)) {
    mf <- file.path(V3_ROOT, "P31R_INPUT_MAP.tsv")
    if (!file.exists(mf)) stop("input not in the release: ", key, " (the release map code/P31R_INPUT_MAP.tsv is missing)", call. = FALSE)
    .p31r_map <<- utils::read.delim(mf, sep = "\t", quote = "", colClasses = "character", na.strings = character(0), check.names = FALSE, encoding = "UTF-8")
  }
  r <- .p31r_map[.p31r_map$source_path == key, , drop = FALSE]
  if (nrow(r) == 0) stop("input not in the release: ", key, " (third-party, withheld or not a released file; see REPRODUCE.md)", call. = FALSE)
  r <- r[1, ]
  if (r$verbatim != "1" && r$lift != "1") stop("input not in the release: ", key, " (the release holds only a reduced or masked copy, ", r$release_path, ")", call. = FALSE)
  q <- file.path(.p31r_arch, r$release_path)
  if (!file.exists(q)) stop("input not in the release: ", key, " (", r$release_path, " is missing from this extraction)", call. = FALSE)
  if (!identical(unname(tools::md5sum(q)), r$md5)) stop("released input ", r$release_path, " does not match its md5 in the release map", call. = FALSE)
  q
}
# level 2 (RC1b): a lifted copy (map column lift = 1) is a declared masked or reduced stand-in for a source the generator reads. The figure-source check compares
# the MD5 OF THE SOURCE with the approved manifest, so it asks for the source md5 recorded in the map (pinned in RELEASE_INPUTS and verified when the release was built).
.p31r_declared_md5 <- function(p) {
  if (file.exists(p) || is.null(.p31r_map)) return(NULL)
  pre <- paste0(V3_ROOT, "/")
  key <- if (startsWith(p, pre)) substr(p, nchar(pre) + 1, nchar(p)) else p
  r <- .p31r_map[.p31r_map$source_path == key, , drop = FALSE]
  if (nrow(r) == 0 || r$lift[1] != "1" || !nzchar(r$source_md5[1])) return(NULL)
  r$source_md5[1]
}
p31r_use <- function(path) {
  p <- .p31r_resolve(normalizePath(path, winslash = "/", mustWork = FALSE))
  if (isTRUE(file.info(p)$size == 0)) stop("data file empty: ", p)
  assign(p, TRUE, envir = .p31r_read_set)
  p
}
p31r_read <- function(path, ...) {
  p <- p31r_use(path)
  # fread reports a truncated or mangled read ("Stopped early on line N", bad quoting) only as a WARNING and returns the
  # rows it got; a figure drawn from part of a file is the CKDGen-truncation defect class, so any warning is fatal here.
  # Warnings are collected and muffled, then raised AFTER fread returns: stopping inside the handler aborted fread before
  # its cleanup, and the next read in the session then failed with "previous fread() session was not cleaned up".
  ws <- character()
  d <- withCallingHandlers(data.table::fread(file = p, ...),
                           warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (length(ws)) stop("fread warning on ", p, ": ", paste(ws, collapse = " | "), call. = FALSE)
  if (nrow(d) == 0) stop("read 0 rows from ", p)
  d
}

# ---- claims ledger: the only source of a number PRINTED on a figure -------------------------------------------
.p31r_led <- NULL
p31r_ledger_row <- function(id) {
  if (is.null(.p31r_led)) {
    # read.csv, not fread: fread (quote rule 0) leaves the CSV's doubled quotes inside JSON values ({""a"": 1}), which
    # then failed fromJSON silently for 36 rows (gate adversary, 3a). Every JSON-looking value must parse, or stop.
    p31r_use(P31R_LEDGER)
    led <- utils::read.delim(P31R_LEDGER, sep = "\t", quote = "", colClasses = "character", encoding = "UTF-8", check.names = FALSE, na.strings = character(0))   # public extract: TSV, JSON values unquoted
    js <- grepl("^[[:space:]]*[[{]", led$value)
    bad <- led$id[js][!vapply(led$value[js], function(v) isTRUE(tryCatch({ jsonlite::fromJSON(v); TRUE }, error = function(e) FALSE)), logical(1))]
    if (length(bad)) stop("ledger values that do not parse as JSON: ", paste(bad, collapse = ", "))
    if (nrow(led) < 100 || anyDuplicated(led$id)) stop("ledger read implausible: ", nrow(led), " rows")
    .p31r_led <<- led
  }
  .want <- id   # data.table scoping: inside [ ], a bare `id` is the COLUMN, so `id == id` selects every row
  r <- as.data.frame(.p31r_led)[which(.p31r_led[["id"]] == .want), , drop = FALSE]
  if (nrow(r) != 1) stop("ledger id not found exactly once: ", id)
  r
}
# value: a scalar or a list element; idx is 1-based for list/array values (JSON in the ledger's value column)
p31r_ledger <- function(id, idx = NULL) {
  v <- p31r_ledger_row(id)$value
  x <- if (grepl("^[[:space:]]*([[{\"]|-?[0-9])", v) || v %in% c("true", "false", "null")) jsonlite::fromJSON(v, simplifyVector = FALSE) else stop("ledger ", id, ": value is not JSON: ", v)
  if (is.list(x) && !is.null(names(x)) && identical(names(x), c("masked", "reason")))
    stop("ledger ", id, " is masked in the public extract (", x$reason, "): ", p31r_ledger_row(id)$display, call. = FALSE)
  # idx: the SAME convention as the figure gate's registry: an integer is a 1-based position in a JSON list, a string is
  # a key of a JSON object, and a vector (or list) of these walks a nested value. Integers on an object are refused.
  if (!is.null(idx)) {
    for (k in as.list(idx)) {
      if (!is.list(x)) stop("ledger ", id, ": idx step ", format(k), " on a scalar")
      if (is.character(k)) {
        if (is.null(names(x)) || !(k %in% names(x))) stop("ledger ", id, ": no key ", k)
      } else if (is.numeric(k)) {
        if (!is.null(names(x))) stop("ledger ", id, ": integer idx ", k, " on an object: use the key")
        if (k != round(k) || k < 1 || k > length(x)) stop("ledger ", id, ": idx ", k, " must be an integer in 1..", length(x))
      } else stop("ledger ", id, ": idx steps are integers or strings")
      x <- x[[k]]
    }
  }
  if (is.list(x) && length(x) == 1 && !is.list(x[[1]])) x <- x[[1]]
  x
}
# printable forms: U+2212 minus; thousands separators for counts
p31r_num <- function(x, digits = 0, big = ",") {
  if (is.list(x) || length(x) != 1 || !is.numeric(suppressWarnings(as.numeric(x))) || !is.finite(suppressWarnings(as.numeric(x))))
    stop("p31r_num: not one finite number: ", paste(format(x), collapse = " "))
  s <- formatC(as.numeric(x), format = "f", digits = digits, big.mark = big)
  sub("^-", "−", s)
}

# ---- save: cairo PDF + PNG rasterised from the PDF + render-log row --------------------------------------------
.p31r_generator <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) != 1) stop("run generators with Rscript <file> so the render log can pin the generator")
  normalizePath(f, winslash = "/")
}
.md5 <- function(p) unname(tools::md5sum(p))
p31r_save <- function(p, name, w_mm = FIG_W_MM, h_mm) {
  stopifnot(w_mm %in% c(FIG_W_MM, FIG_W1_MM), h_mm > 0, h_mm <= FIG_HMAX_MM)
  out <- file.path(P31R_FIGS, "out"); dir.create(out, showWarnings = FALSE)
  pdf_f <- file.path(out, paste0(name, ".pdf")); png_f <- file.path(out, paste0(name, ".png"))
  ggsave(pdf_f, p, width = w_mm / 25.4, height = h_mm / 25.4, units = "in", device = cairo_pdf)
  # author/creator metadata (global CLAUDE.md section 8): cairo_pdf takes none, so PyMuPDF sets it BEFORE the md5 is pinned
  rc <- system2("python", c("-c", shQuote(paste0("import fitz,sys,os; p=sys.argv[1]; d=fitz.open(p); ",
         "d.set_metadata({'author':'Bertrand Chin-Ming Tan','creator':'Bertrand Chin-Ming Tan','producer':d.metadata.get('producer',''),",
         "'title':os.path.splitext(os.path.basename(p))[0].replace('_',' ')}); d.save(p+'.tmp'); d.close(); os.replace(p+'.tmp', p)")), shQuote(pdf_f)))
  if (rc != 0) stop("could not set PDF metadata on ", pdf_f)
  bm <- pdftools::pdf_render_page(pdf_f, page = 1, dpi = 300, numeric = FALSE)
  png::writePNG(bm, png_f, dpi = 300, text = c(Author = "Bertrand Chin-Ming Tan"))   # pHYs = 300 dpi, so viewers report the true size   # raster OF the shipped PDF
  gen <- .p31r_generator()
  data <- sort(ls(.p31r_read_set))
  row <- data.table(figure = name, w_mm = w_mm, h_mm = h_mm,
                    generator = sub(paste0("^", V3_ROOT, "/"), "", gen), generator_md5 = .md5(gen),
                    common_md5 = .md5(P31R_COMMON), deps = paste0(sub(paste0("^", V3_ROOT, "/"), "", P31R_DEPS), "=", .md5(P31R_DEPS), collapse = ";"),
                    # pins under V3 are written relative to it (as generator and deps), so a copied tree checks its own files
                    data = paste0(sub(paste0("^", V3_ROOT, "/"), "", data), "=", .md5(data), collapse = ";"),
                    pdf_md5 = .md5(pdf_f), png_md5 = .md5(png_f), rendered = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  # one log file PER FIGURE: generators run concurrently, and a shared read-modify-write log would drop rows
  logd <- file.path(P31R_FIGS, "render_log"); dir.create(logd, showWarnings = FALSE)
  data.table::fwrite(row[, lapply(.SD, as.character)], file.path(logd, paste0(name, ".tsv")), sep = "\t")
  message("saved ", pdf_f, " (", w_mm, " x ", h_mm, " mm) + raster PNG; render log updated")
  invisible(pdf_f)
}
