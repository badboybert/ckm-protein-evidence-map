# 02_candidates.R -- the 16 multi-signal candidates and their SuSiE windows, from the pool rules (DEF_SM_CANDIDATES,
# LOG_28) and the run logs 29/30/31 (W6a-independent read). Pool composition (kidney vs other) from the 29 log headers.
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
say("==== W6b 02_candidates ====")
res("DEF_SM_CANDIDATES"); res("LOG_28_candidate_pool")
L29e <- readLines(res("LOG_29_cand_egfr"), warn = FALSE)
L29n <- readLines(res("LOG_29_cand_nonegfr"), warn = FALSE)
L30  <- readLines(res("LOG_30_candidate_audit"), warn = FALSE)
L31e <- readLines(res("LOG_31_candidate_egfr"), warn = FALSE)
L31n <- readLines(res("LOG_31_candidate_nonegfr"), warn = FALSE)

# --- 31 strict-SuSiE summary tables (gene, outcome, orig, strict, n, n_warn, run_status, survives) ---
parse31 <- function(L, mode) {
  i0 <- grep("^\\s+gene outcome_name\\s+orig\\s+strict", L)
  if (length(i0) != 1) stop("31 ", mode, ": summary header not unique")
  out <- list(); i <- i0 + 2
  while (i <= length(L) && grepl("^\\s*[0-9]+:\\s", L[i])) {
    f <- strsplit(trimws(sub("^\\s*[0-9]+:\\s*", "", L[i])), "\\s+")[[1]]
    out[[length(out) + 1]] <- data.table(gene = f[1], outcome = f[2], orig = as.numeric(f[3]), strict = suppressWarnings(as.numeric(f[4])),
                                         n_cs = f[5], n_warn = as.integer(f[6]), run_status = f[7], mode = mode)
    i <- i + 1
  }
  dt <- rbindlist(out)
  # survives column is printed on the continuation block (egfr) or the second block (nonegfr): re-read it
  j0 <- grep("^\\s+(collapse_reason|survives)", L)
  sv <- c()
  for (j in j0) { k <- j + 2; while (k <= length(L) && grepl("^\\s*[0-9]+:\\s", L[k])) { f <- strsplit(trimws(sub("^\\s*[0-9]+:\\s*", "", L[k])), "\\s+")[[1]]; sv <- c(sv, f[1]); k <- k + 1 } }
  dt[, survives_raw := NA_character_]
  # egfr: the 'survives' value is the last field on the first block; nonegfr: first field of the second block
  first_block_last <- vapply(seq_len(nrow(dt)), function(r) { f <- strsplit(trimws(sub("^\\s*[0-9]+:\\s*", "", L[i0 + 1 + r])), "\\s+")[[1]]; f[length(f)] }, "")
  dt[, survives_raw := ifelse(first_block_last %in% c("TRUE", "FALSE"), first_block_last, sv[seq_len(nrow(dt))])]
  dt[, survives := survives_raw == "TRUE"]
  m <- regmatches(L, regexpr("SURVIVED\\(strict PP.H4>=0.80\\)=[0-9]+", L))
  say("31 %s: %d edges parsed; log line: %s", mode, nrow(dt), m)
  stopifnot(sum(dt$survives) == as.integer(sub(".*=", "", m)))
  dt
}
c31 <- rbind(parse31(L31e, "egfr"), parse31(L31n, "nonegfr"))
# candidate rule (SM_candidates + log 28/31): strict PP.H4 >= 0.80 with converged fits (no warning)
c31[, candidate := survives & strict >= 0.80 & n_warn == 0L]
cand <- c31[candidate == TRUE]
held <- c31[survives & n_warn > 0]
say("candidates: %d (%s); held for a warning (non-converged fit): %s", nrow(cand), paste(cand$gene, cand$outcome, sep = "~", collapse = ", "),
    paste(held$gene, held$outcome, sep = "~", collapse = ", "))

# 30: 26 dashed edges; every candidate must be one of them
i30 <- grep("== DASHED EDGES", L30)
d30 <- c(); k <- i30 + 1
while (!grepl("^\\s*[0-9]+:\\s", L30[k])) k <- k + 1        # skip the column-name and <type> rows
while (grepl("^\\s*[0-9]+:\\s", L30[k])) { f <- strsplit(trimws(sub("^\\s*[0-9]+:\\s*", "", L30[k])), "\\s+")[[1]]; d30 <- c(d30, paste(f[1], f[2], sep = "~")); k <- k + 1 }
say("30: %d dashed edges parsed", length(d30))
stopifnot(length(d30) == 26, all(paste(cand$gene, cand$outcome, sep = "~") %in% d30))

# cross-check against the V4 edge audit table (raw_input): Status == Drawn
ea <- fread(res("TAB_edge_audit"), encoding = "UTF-8")
drawn <- ea[Status == "Drawn", paste(Gene, Endpoint, sep = "~")]
same <- setequal(drawn, paste(cand$gene, cand$outcome, sep = "~"))
log_assert("S6.d.candidates_vs_edge_audit", "identical sets", same, same, "31-log rule (strict>=0.80, n_warn 0) vs TAB_edge_audit Status=Drawn")

# --- windows from the 29 logs: '  susie: GENE ~ LABEL | chrC:lo-hi' ---
win <- rbindlist(lapply(c(L29e, L29n)[grepl("^\\s+susie: \\S+ ~ \\S+ \\| chr", c(L29e, L29n))], function(l) {
  m <- regmatches(l, regexec("susie: (\\S+) ~ (\\S+) \\| chr([0-9XY]+):([0-9]+)-([0-9]+)", l))[[1]]
  data.table(gene = m[2], label29 = m[3], chr = m[4], lo38 = as.integer(m[5]), hi38 = as.integer(m[6]))
}))
lab2out <- c(eGFR = "eGFR", Wuttke2019_CKD_EA = "CKD", FinnGen_R10_T2D = "T2D", FinnGen_R10_I9_CHD = "CAD", FinnGen_R10_I9_HEARTFAIL = "HF",
             FinnGen_R10_I9_AF = "AF", MEGASTROKE_EUR = "Stroke")
win[, outcome := unname(lab2out[label29])]
win <- unique(win)
cand <- merge(cand, win, by = c("gene", "outcome"), all.x = TRUE, sort = FALSE)
if (anyNA(cand$lo38)) stop("candidate without a 29-log window: ", paste(cand[is.na(lo38)]$gene, collapse = ","))
if (anyDuplicated(cand[, .(gene, outcome)])) stop("a candidate has more than one 29-log window")
cand[, kidney := outcome %in% c("eGFR", "CKD")]
cand[, cid := paste(gene, outcome, sep = "_")]
setorder(cand, -kidney, gene)
print(cand[, .(cid, gene, outcome, chr, lo38, hi38, label29, strict, n_warn, kidney)])
fwrite(cand[, .(cid, gene, outcome, label29, chr, lo38, hi38, orig, strict, n_cs, n_warn, kidney)], file.path(WORK, "candidates.tsv"), sep = "\t")

# --- pool composition from the 29 headers: [k/N] GENE ~ OUTCOME ---
hdr <- function(L, N) {
  h <- regmatches(L, regexpr(sprintf("^\\[[0-9]+/%d\\] \\S+ ~ \\S+", N), L))
  dt <- rbindlist(lapply(h, function(x) { m <- regmatches(x, regexec("^\\[([0-9]+)/[0-9]+\\] (\\S+) ~ (\\S+)", x))[[1]]; data.table(k = as.integer(m[2]), gene = m[3], outcome = m[4]) }))
  unique(dt)
}
pe <- hdr(L29e, 61); pn <- hdr(L29n, 99)
stopifnot(nrow(pe) == 61, setequal(pe$k, 1:61), nrow(pn) == 99, setequal(pn$k, 1:99))
pool <- rbind(pe, pn)
stopifnot(!anyDuplicated(pool[, .(gene, outcome)]))
pool[, kidney := outcome %in% c("eGFR", "CKD")]
pool[, candidate := paste(gene, outcome, sep = "~") %in% paste(cand$gene, cand$outcome, sep = "~")]
tab <- pool[, .N, by = .(candidate, kidney)]
a <- pool[candidate == TRUE & kidney == TRUE, .N]; b <- pool[candidate == TRUE & kidney == FALSE, .N]; c_ <- pool[candidate == FALSE & kidney == TRUE, .N]; d <- pool[candidate == FALSE & kidney == FALSE, .N]
say("pool from 29 headers: %d pairs (%d eGFR, %d CKD, %d other); candidates %d (kidney %d)", nrow(pool), pool[outcome == "eGFR", .N],
    pool[outcome == "CKD", .N], pool[kidney == FALSE, .N], a + b, a)
fwrite(pool, file.path(WORK, "pool160_from_29_headers.tsv"), sep = "\t")
# assertion targets from the pre-declaration (sec 6 primary: 11/16 vs 57/144) -- computed first, compared second (rule 9)
log_assert("S6.d.pool.rows", 160, nrow(pool), nrow(pool) == 160, "unique [k/61] + [k/99] headers of LOG_29_cand_egfr / _nonegfr")
log_assert("S6.d.candidates", 16, a + b, a + b == 16, "LOG_31 strict PP.H4 >= 0.80 and n_warn == 0 (egfr + nonegfr)")
log_assert("S6.d.primary.table", "11/16 vs 57/144", sprintf("%d/%d vs %d/%d", a, a + b, c_, c_ + d), a == 11 && b == 5 && c_ == 57 && d == 87,
           "candidates x kidney (eGFR or CKD) over the 29-log pool")
saveRDS(list(a = a, b = b, c = c_, d = d), file.path(WORK, "primary_table_from_logs.rds"))
say("02_candidates done")
