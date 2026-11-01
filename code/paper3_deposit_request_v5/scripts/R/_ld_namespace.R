# ============================================================================
# Translating between the two rsID namespaces an instrument can live in.
#
# WHY THIS EXISTS
# ---------------
# A variant can carry one name in the outcome GWAS and a different one in the
# 1000G EUR reference panel -- merged or retired identifiers for the same
# variant, matched on position AND on the unordered allele pair, so they are two
# names for one variant rather than two variants.
#
# 09_phase1_harmonise_hg38.R therefore emits both:
#   SNP            the name the OUTCOME GWAS resolve. Every outcome lookup and
#                  every join keyed on a variant uses this.
#   rsid_ld_panel  the name the REFERENCE PANEL knows, and NA where the panel
#                  does not contain the variant at all (about 15% of rows).
#
# One column could not serve both: with SNP holding the panel's name the outcome
# lookup missed variants that were present in the file; with SNP holding the
# outcome's name PLINK cannot find them in the panel.
#
# THE TWO DIRECTIONS ARE NOT SYMMETRIC
# ------------------------------------
#   SNP -> panel   ld_panel_key(). "Not in the map" means the panel does not
#                  have this variant, so the answer is NA: do not query it.
#                  Falling back to SNP here would ask PLINK for a name the panel
#                  was already known not to hold, and if that string happens to
#                  exist for a record the position/allele matching REJECTED, LD
#                  from the wrong record comes back and is used.
#   panel -> SNP   ld_names_to_snp(). "Not in the map" means this is an ordinary
#                  panel variant that is not one of our instruments -- a proxy --
#                  and it must pass through unchanged.
# ============================================================================
suppressPackageStartupMessages(library(data.table))

# The identifier to hand to PLINK / ieugwasr for each row, or NA where the panel does not have the
# variant. Callers must drop the NAs; they are not a failure, they are the panel's coverage.
ld_panel_key <- function(inst) {
  stopifnot(is.data.frame(inst), "SNP" %in% names(inst))
  if (!"rsid_ld_panel" %in% names(inst)) {
    # LEGACY tables, which is every instrument table built before this change. Their SNP column IS
    # the panel's name -- the old rule resolved disagreements in the panel's favour -- so returning
    # it reproduces exactly what these consumers did before, with no row lost.
    #
    # THE CALLER'S CONTRACT, because this function cannot enforce it -- and the WORST CASE, measured
    # rather than argued (2026-08-01, against the pinned panel):
    #
    #   1,964 instruments = 1,757 whose SNP and rsid_ld_panel agree
    #                     +   202 with no panel counterpart (rsid_ld_panel NA)
    #                     +     5 aliases (ADM, CD84, CEP170, LRRC37A2, PGD)
    #   Of those 207, how many of their SNP names exist anywhere in EUR.bim (8,550,156 variants)?
    #   ZERO.
    #
    # So if a CURRENT table lost this column through a projection and landed here, the header's worst
    # case -- "LD from the wrong record comes back and is used" -- has no instance on this panel:
    # every affected row simply returns nothing. That is availability, not correctness, and it is not
    # silent (01_mr_screen warns and writes *_UNSEEN_KEYS.txt). The guarantee is tied to the pinned
    # panel: 01 compares the panel files' md5 against the Phase-1 certificate before PLINK runs.
    #
    # The fallback stays. A data frame carries no signal separating a genuine legacy table from a
    # current projection -- 1,757 identifiers are identical in both -- so `data.frame(SNP = "rs...")`
    # is consistent with either, and making the column mandatory would break legitimate reruns of the
    # 20-column legacy table. So: ANY PROJECTION PASSED HERE MUST KEEP rsid_ld_panel. Every current
    # caller passes the whole read_phase1_instruments() result or a full readRDS; the `inst_al`
    # projections elsewhere are not passed to these helpers. The same "no column == legacy" test
    # appears at seven other decision points (01:245, 01:307, 02:61, 03:94, 03:166, 53:199, 53:228),
    # so this note documents a convention, not a single site.
    return(as.character(inst$SNP))
  }
  as.character(inst$rsid_ld_panel)
}

# panel name -> outcome-resolvable SNP, for translating a panel-space artifact (a PLINK proxy cache,
# an LD-matrix dimname) back into the namespace every join uses.
ld_to_snp_map <- function(inst) {
  k <- ld_panel_key(inst); s <- as.character(inst$SNP)
  ok <- !is.na(k) & nzchar(k) & !is.na(s) & nzchar(s)
  # CONTRADICTORY INPUT IS REFUSED, not silently resolved. Dropping the NA rows first meant an
  # instrument SNP appearing twice -- legitimately, when one variant instruments two genes -- with a
  # panel name on one row and NA on the other survived as "this SNP maps to that panel name". The
  # reverse map then answers a query for the row explicitly marked "not in the panel", and 53 goes on
  # to use panel LD for a variant 09 said the panel does not hold. The two rows disagree about a fact;
  # neither is a safe default.
  .bad <- unique(s[!is.na(s) & nzchar(s) & s %in% s[ok] & !ok])
  if (length(.bad))
    stop("_ld_namespace.R: ", length(.bad), " instrument identifier(s) carry a reference-panel name ",
         "on one row and none on another (e.g. ", paste(utils::head(.bad, 3), collapse = ", "),
         "). The table contradicts itself about whether the panel holds these variants.")
  k <- k[ok]; s <- s[ok]
  pair <- !duplicated(paste(k, s))
  k <- k[pair]; s <- s[pair]
  # One panel name must not stand for two instruments, or a translated proxy cache would attach one
  # variant's proxies to another's.
  if (anyDuplicated(k))
    stop("_ld_namespace.R: ", sum(duplicated(k)), " reference-panel identifier(s) map to more than ",
         "one instrument (e.g. ", paste(utils::head(unique(k[duplicated(k)]), 3), collapse = ", "),
         "). Translating panel-space results back would attach them to the wrong variant.")
  # ...and the REVERSE must be a function too, because 53 inverts this map to go SNP -> panel. One
  # SNP with two different panel names would make that inversion pick whichever came first: 03 would
  # query both panel variants while 53 used only one.
  if (anyDuplicated(s))
    stop("_ld_namespace.R: ", sum(duplicated(s)), " instrument identifier(s) carry more than one ",
         "reference-panel name (e.g. ", paste(utils::head(unique(s[duplicated(s)]), 3),
                                              collapse = ", "),
         "). The SNP -> panel direction would not be well defined.")
  m <- setNames(s, k)
  # The FULL outcome-identifier set travels with the map. ld_names_to_snp() passes an unmapped panel
  # name through unchanged -- correct for an ordinary proxy -- but if that string is also the
  # outcome-space SNP of a DIFFERENT instrument, the pass-through resolves the proxy to that
  # instrument's outcome association: its beta is then attached to the gene the proxy was found for,
  # with LD and phase belonging to an unrelated panel record. The names are needed to tell the two
  # cases apart, so they are carried rather than looked up again by every caller.
  attr(m, "outcome_ids") <- unique(as.character(inst$SNP))
  m
}

# SNP -> panel, for the callers that must QUERY the panel. NA where the panel has no such variant;
# this is deliberately NOT a pass-through, for the reason given at the top.
ld_snp_to_panel_map <- function(inst) {
  m <- ld_to_snp_map(inst)
  setNames(names(m), unname(m))
}

# Translate panel-space names into SNP-space. Names absent from the map pass through UNCHANGED: a
# proxy variant is not an instrument, so it never appears in the map and must survive. Where a proxy
# IS also one of our instruments -- 27 such variants and 34 proxy edges in the current cache -- the
# mapping is known and is applied, which is exactly what this pass-through gives.
ld_names_to_snp <- function(x, map) {
  x <- as.character(x)
  hit <- !is.na(x) & x %in% names(map)
  x[hit] <- unname(map[x[hit]])
  # An unmapped name that is nonetheless one of our outcome identifiers is NOT an ordinary proxy: it
  # is a panel record the position/allele matching refused for that instrument, and letting it
  # through would attach that instrument's outcome association to whichever gene the proxy was found
  # for. NA means "unusable", which every consumer already drops.
  oi <- attr(map, "outcome_ids")
  if (!is.null(oi)) {
    bad <- !hit & !is.na(x) & x %in% oi
    if (any(bad)) x[bad] <- NA_character_
  }
  x
}

# SNP-space names into panel space. Unmapped becomes NA -- "do not query" -- rather than passing
# through, which is the difference between this and ld_names_to_snp().
ld_names_to_panel <- function(x, map) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  hit <- !is.na(x) & x %in% names(map)
  out[hit] <- unname(map[x[hit]])
  out
}

# The EXACT set of identifiers handed to PLINK. 03, 16 and cascade/01 build their snplist from it,
# and 02 fingerprints the proxy cache against it -- and they must agree to the character, or the
# cache 03 writes is one 02 refuses. They disagreed once: 03 filtered to rs-format while 02 hashed
# every non-NA key including the 79 chr:pos entries in the legacy table, so the legacy path had no
# working route at all. One definition, used by both.
ld_query_keys <- function(inst) {
  k <- ld_panel_key(inst)
  # A NON-EMPTY key that is not rs-format is a CONTRADICTION, not a row to drop. `rsid_ld_panel`
  # asserts "the panel holds this variant under this name"; a malformed or ss-prefixed value asserts
  # it while naming something PLINK cannot be asked for. Filtered silently, the instrument got no
  # proxies, was absent from the .unseen report that exists to name exactly this loss, and its column
  # still claimed panel coverage. NA is different and stays a legitimate "not in the panel".
  .bad <- unique(k[!is.na(k) & nzchar(k) & !grepl("^rs[0-9]+$", k)])
  # The legacy path returns SNP itself as the key, and the shipped legacy table really does carry 79
  # chr:pos identifiers there; those keep their long-standing exclusion rather than aborting a
  # reproduction of the previous analysis.
  if (length(.bad) && !is.null(inst$rsid_ld_panel))
    stop("_ld_namespace.R: ", length(.bad), " instrument(s) carry a non-rs reference-panel ",
         "identifier (e.g. ", paste(utils::head(.bad, 3), collapse = ", "), "). PLINK cannot be ",
         "asked for it, so the variant would silently have no proxies while its column says the ",
         "panel holds it. Set rsid_ld_panel to NA if the panel does not name it.")
  sort(unique(k[!is.na(k) & grepl("^rs[0-9]+$", k)]))
}

# The fingerprint of that set, as recorded by 03 and checked by 02.
ld_query_index_md5 <- function(inst) {
  digest::digest(paste(ld_query_keys(inst), collapse = "\n"), algo = "md5", serialize = FALSE)
}
