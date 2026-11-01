# The Phase 1 instrument set every downstream script reads. ONE definition.
#
# This was 66 hardcoded "repro_phase1_1682" strings across 59 files. That is not merely repetitive:
# it makes a PARTIAL switch possible, where part of the cascade reads one instrument set and part
# reads another, and every table produced in between looks perfectly normal. The directory name also
# encodes a protein count, which stops being true the moment the set is rebuilt.
#
# Overridable by pre-setting, matching the idiom R/00_config.R already uses for RESULTS_DIR.
# THE DEFAULT IS THE REBUILT SET, as of 2026-07-31 when 09 produced it. It used to be
# repro_phase1_1682, which meant that after the rebuild landed the clean default path -- source
# config with no environment variables -- hit the guard below and aborted before a single threshold
# was defined. A default that always fails is not a safe default, it is a broken one; the safe
# default is the generation the analysis is actually meant to use.
#
# "0"/"false" IS NOT PERMISSION. nzchar() is TRUE for both, so setting a flag to a FALSE value used
# to license exactly what it looks like it forbids; and as.logical() alone is not the fix either,
# because as.logical("1") is NA in R, so the affirmative spelling would have been refused. The
# affirmative values are named, ONCE -- this file and 00_config.R both gate on environment flags,
# and a second copy of the list is a second chance to disagree with the first.
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
ckm_env_yes <- function(name) tolower(Sys.getenv(name)) %in% c("1", "true", "t", "yes")

# THE PREVIOUS GENERATIONS, BY CONTENT. CKM_ALLOW_STALE_PHASE1 says it selects "the previous
# generation" and both messages that mention it say so, but it waived the rebuild check and the
# missing-certificate check for ANY readable table: a probe with an invented set name printed
# "deliberately using the PREVIOUS generation" and would have derived manuscript thresholds from
# whatever that path held. An override is only as honest as the thing it authorises, so it
# authorises these specific tables. Add a line here, deliberately, to admit another.
CKM_PREVIOUS_GENERATION_MD5 <- c(
  repro_phase1_1682  = "6f14eb19759daa8ca6370076597dcdae",   # 1,723 instruments / 1,682 genes
  phase1_instruments = "72af3e4027462a562f901c7d2e19048d"    # the 1,635-gene table before it
)
ckm_is_known_previous <- function(path) {
  .m <- tryCatch(tolower(unname(tools::md5sum(path))), error = function(e) NA_character_)
  !is.na(.m) && .m %in% tolower(CKM_PREVIOUS_GENERATION_MD5)
}

# nzchar(), not the two-argument Sys.getenv(): its default applies only when the variable is
# ABSENT, so a launcher that represents "unset" as "" -- which is how an unset optional variable
# usually reaches a subprocess -- selected results//ukb_cis_instruments.rds and aborted with no
# thresholds at all. An empty value means "no preference", which is what the default is for.
if (!exists("PHASE1_SET"))
  PHASE1_SET <- if (nzchar(Sys.getenv("CKM_PHASE1_SET"))) Sys.getenv("CKM_PHASE1_SET") else
                "phase1_hg38"
# ...AND RE-SOURCING WITH A CHANGED SELECTION IS A CONTRADICTION, not a no-op. The exists() guard
# above retains whatever the session already chose, so changing CKM_PHASE1_SET and sourcing again
# silently kept the first choice: measured, phase1_hg38 stayed selected after switching the variable
# to repro_phase1_1682, and a script that believed it had switched read 1,948 genes. The threshold
# stays consistent with the retained file, so nothing downstream notices.
if (nzchar(Sys.getenv("CKM_PHASE1_SET")) &&
    !identical(PHASE1_SET, Sys.getenv("CKM_PHASE1_SET")))
  stop(sprintf(paste0("_phase1_set.R: CKM_PHASE1_SET asks for %s but this session already selected ",
                      "%s. Re-sourcing does not switch generations -- start a fresh R process, or ",
                      "unset the variable to keep the current selection."),
               Sys.getenv("CKM_PHASE1_SET"), PHASE1_SET))
.PHASE1_FROM_SET <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"), "results", PHASE1_SET,
                              "ukb_cis_instruments.rds")
if (!exists("PHASE1_INSTRUMENTS")) PHASE1_INSTRUMENTS <- .PHASE1_FROM_SET
# ...AND THE TWO MUST AGREE. Both are retained by their own exists() guard, so pre-setting
# PHASE1_INSTRUMENTS alone left PHASE1_SET naming one generation and the path naming another: a
# consumer that builds results/<PHASE1_SET>/ukb_cis_instruments.rds then read 1,948 genes while
# everything using PHASE1_INSTRUMENTS read 1,682, and the Bonferroni threshold followed the latter
# (3.7158e-06 against a 1,948-gene analysis, 15.8% too loose -- three stored MR rows, CTRB1-T2D,
# ITIH1-eGFR and CEP170-eGFR, sit between the two thresholds). Re-sourcing after changing only
# CKM_PHASE1_SET was silently ignored for the same reason.
if (!identical(normalizePath(PHASE1_INSTRUMENTS, mustWork = FALSE),
               normalizePath(.PHASE1_FROM_SET, mustWork = FALSE)))
  stop(sprintf(paste0("_phase1_set.R: PHASE1_SET is %s, which names %s, but PHASE1_INSTRUMENTS is ",
                      "%s. A consumer that follows the set name and one that follows the path would ",
                      "read different instrument generations. Set CKM_PHASE1_SET, or pre-set both ",
                      "consistently."),
               PHASE1_SET, .PHASE1_FROM_SET, PHASE1_INSTRUMENTS))

# ── A REBUILD IN PROGRESS STOPS THE LAYER, rather than letting it finish on the old set ───────────
# Phase 1 was rebuilt because the cis windows were computed with GRCh37 gene boundaries against
# GRCh38 variant positions. While that rebuild is only half-landed -- the extraction has produced
# results/phase1_hg38/, the harmonised RDS it feeds does not exist yet -- every script here still
# resolves to a PREVIOUS-generation table and runs to completion looking entirely successful.
# Measured against the rebuilt extraction (2,053 instruments / 2,036 genes): the 1,635-gene table
# omits 427 rebuilt genes and keeps 26 the rebuild does not have; the 1,682-gene table omits 385 and
# keeps 31. Results produced in that state are the old generation wearing the new one's timestamps,
# and they land in the same shared directories as genuinely rebuilt artifacts.
#
# So: if the rebuilt extraction exists and the selected set is not derived from it, refuse. This is
# the only place that knows both facts. Set CKM_ALLOW_STALE_PHASE1=1 to run the previous generation
# deliberately -- for a comparison, or to reproduce an old number -- which is a different act from
# reaching it by default.
local({
  .rebuilt <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"), "results", "phase1_hg38",
                        "ukb_cis_instruments.csv")
  if (!file.exists(.rebuilt)) return(invisible())              # no rebuild under way: nothing to say
  # THE FILE'S CONTENTS, not the spelling of its path. `grepl("phase1_hg38", ...)` passed for any RDS
  # under a directory with that name whatever it held. A table derived from the rebuild can only
  # contain genes the rebuild found; the previous generations each carry genes it does not have (26
  # in the 1,635-gene table, 31 in the 1,682-gene one), so a single set difference separates them
  # without needing anything new to be written anywhere.
  # INSTRUMENT IDENTITY, not gene membership. Comparing gene sets alone let a mosaic through: the two
  # generations share 1,651 genes, and among those, 351 of the old table's 1,692 (gene, chrpos) pairs
  # do not occur in the rebuild -- for 320 common genes the two sets of instrument positions do not
  # overlap at all. So a table carrying new genes beside stale rows for the common ones would pass,
  # and its own md5 would then be stamped into every extract and MR result, making every downstream
  # equality check bless the mosaic instead of detecting it. The (gene, chrpos) pair is the identity
  # that matters, and it costs one more column to compare.
  # The two carry the position under different names -- the instrument RDS as `chrpos` ("19:58347310")
  # and the extraction CSV as CHROM + GENPOS -- so the key is built rather than read.
  #
  # CONTAINMENT IS NOT ENOUGH, and neither is the key alone. "No key the rebuild lacks" is satisfied
  # by any SUBSET: an RDS holding only the 1,341 instruments the two generations happen to share
  # passes while dropping 712 rebuilt instruments and 705 proteins. And a shared key does not mean a
  # shared row -- TNF has the same (gene, chrpos) in both, but the old row is OID20848 with beta
  # -0.102519 while the rebuild has OID20074 and -0.091522, so a mosaic carrying the stale row would
  # pass on keys alone and put a previous generation's exposure value into a new MR estimate.
  #
  # So: every key the selected table holds must be in the rebuild (no stale variants), and every row
  # they share must carry the same beta (no stale values). It is deliberately NOT required to be
  # complete -- 09 removes instruments it cannot name, so a legitimate harmonised table is a strict
  # subset -- but the size of that gap is itself recorded by 09 and is checked where it is known.
  # ...AND THE INTERVAL, not only the estimate. Key plus beta still admits the substitution this
  # generation exists to prevent: a row whose gene and chrpos and beta all agree with the rebuild,
  # but whose gene_start/gene_end/cis_start/cis_end came from GRCh37. Its own md5 and row count can
  # be updated in the manifest, and every check above would pass while colocalisation reconstructed
  # a region the instruments were never selected in. Both files carry these five columns under the
  # same names, so the signature costs one paste and no extra read.
  .wsig <- function(d) paste(d$gene_chrom, d$gene_start, d$gene_end, d$cis_start, d$cis_end, sep = ":")
  .sel <- tryCatch({ .d <- as.data.frame(readRDS(PHASE1_INSTRUMENTS))
                     data.frame(k = paste(.d$gene, .d$chrpos), b = as.numeric(.d$beta),
                                w = .wsig(.d), stringsAsFactors = FALSE) }, error = function(e) NULL)
  .reb <- tryCatch({ .d <- utils::read.csv(.rebuilt, stringsAsFactors = FALSE)
                     data.frame(k = paste(.d$gene, paste0(.d$CHROM, ":", .d$GENPOS)),
                                b = as.numeric(.d$BETA),
                                w = .wsig(.d), stringsAsFactors = FALSE) },
                   error = function(e) NULL)
  if (!is.null(.sel) && !is.null(.reb) && nrow(.sel) && nrow(.reb) &&
      !length(setdiff(.sel$k, .reb$k))) {
    .m <- merge(unique(.sel), unique(.reb), by = "k", suffixes = c(".sel", ".reb"))
    .drift <- .m[(is.finite(.m$b.sel) & is.finite(.m$b.reb) &
                  abs(.m$b.sel - .m$b.reb) > 1e-9) | .m$w.sel != .m$w.reb, , drop = FALSE]
    if (!nrow(.drift)) return(invisible())
    .isw <- .drift$w.sel[1] != .drift$w.reb[1]
    stop(sprintf(paste0("_phase1_set.R: %s holds %d instrument(s) whose (gene, position) the rebuild ",
                        "also has but with a DIFFERENT %s (e.g. %s: %s here, %s in %s). The keys ",
                        "match, so a key-only check would pass this table and a previous generation's ",
                        "%s would enter the analysis."),
                 PHASE1_INSTRUMENTS, nrow(.drift),
                 if (.isw) "gene/cis interval" else "exposure estimate", .drift$k[1],
                 if (.isw) .drift$w.sel[1] else format(.drift$b.sel[1]),
                 if (.isw) .drift$w.reb[1] else format(.drift$b.reb[1]), .rebuilt,
                 if (.isw) "coordinates" else "beta"))
  }
  # ...and "0"/"false" is not permission -- see ckm_env_yes() at the top of this file. Nor is the
  # flag a licence for an arbitrary table: it names a previous generation, so it grants one.
  if (ckm_env_yes("CKM_ALLOW_STALE_PHASE1") && !ckm_is_known_previous(PHASE1_INSTRUMENTS))
    stop(sprintf(paste0("_phase1_set.R: CKM_ALLOW_STALE_PHASE1 authorises a PREVIOUS GENERATION, ",
                        "but %s is not one of them (%s). It is neither the rebuild nor a recorded ",
                        "previous generation, so nothing states what it is. Add its md5 to ",
                        "CKM_PREVIOUS_GENERATION_MD5 if it really is one."),
                 PHASE1_INSTRUMENTS, paste(names(CKM_PREVIOUS_GENERATION_MD5), collapse = ", ")))
  if (ckm_env_yes("CKM_ALLOW_STALE_PHASE1")) {
    cat(sprintf("[_phase1_set] CKM_ALLOW_STALE_PHASE1 is set -- deliberately using the PREVIOUS generation: %s\n",
                PHASE1_INSTRUMENTS))
    return(invisible())
  }
  stop(sprintf(paste0(
    "_phase1_set.R: a GRCh38 rebuild of Phase 1 exists (%s) but this run would use %s, which is a ",
    "previous generation. Every table produced from it would look successful while being built on ",
    "the cis windows the rebuild exists to replace. Finish the rebuild and point CKM_PHASE1_SET at ",
    "it, or set CKM_ALLOW_STALE_PHASE1=1 to use the previous generation on purpose."),
    .rebuilt, PHASE1_INSTRUMENTS))
})
