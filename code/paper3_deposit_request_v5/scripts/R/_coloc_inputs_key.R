# ============================================================================
# One definition of "which inputs this colocalisation run depends on".
#
# The split (21), the workers (22) and the merge (23) each need to state what a
# result was computed from, and the three statements are only comparable if they
# come from the same code -- so they come from here.
#
# Why it exists: the split manifest used to pin the pair set, the eGFR id and the
# coordinate table, and nothing else. Everything else that decides a coloc result
# -- the engine in R/coloc.R, the configuration, the instrument table, the rsID
# map, the liftover chain, the outcome GWAS files, the coloc priors -- could be
# replaced between two chunks, or between the chunks and the merge, and every
# coverage check would still pass. A failed chunk re-run after any of those
# changes merged silently into results from before them.
# ============================================================================

# Content hashes for everything small enough to hash in a second or two.
.coloc_key_hashed <- function() {
  b <- BASE_DIR
  f <- c(
    engine        = file.path(b, "R/coloc.R"),
    config        = file.path(b, "R/00_config.R"),
    phase1_set    = file.path(b, "R/_phase1_set.R"),
    worker        = file.path(b, "scripts/_shared/localize/22_coloc_chunk_worker.R"),
    splitter      = file.path(b, "scripts/_shared/localize/21_coloc_split.R"),
    # THIS file. The protocol compares fingerprints, so all three participants must be computing
    # them the same way; the one script guaranteeing that was the one it did not cover.
    inputs_key    = file.path(b, "R/_coloc_inputs_key.R"),
    # The MERGER too. It is the script that assembles the final scientific table -- it decides row
    # order, applies the coverage assertions and writes the published files -- and it was the one
    # participant in this protocol absent from the fingerprint the protocol compares.
    merger        = file.path(b, "scripts/_shared/localize/23_merge_parallel_coloc.R"),
    coords        = GENE_HG38_COORDS,
    coords_cert   = paste0(GENE_HG38_COORDS, ".manifest.tsv"),
    instruments   = PHASE1_INSTRUMENTS,
    # mr_suggestive is fingerprinted on the FIELDS coloc reads, below, not on its bytes: re-saving
    # the identical object with different RDS compression, or changing Q or OR, altered the byte
    # hash and made every worker reject unchanged scientific work.
    rsid_map      = file.path(b, "cache/cis_rsid_map.tsv"),
    # The chain the engine ACTUALLY reads. liftover_to_hg37() passes LIFTOVER_CHAIN_38_TO_37
    # (R/coloc.R:318); this hashed LIFTOVER_CHAIN, the hg19->hg38 chain, which no coloc code path
    # reads -- so replacing the downward chain left the fingerprint unchanged, which is the one
    # substitution this key exists to catch.
    chain_38_37   = if (exists("LIFTOVER_CHAIN_38_TO_37", envir = globalenv()))
                      get("LIFTOVER_CHAIN_38_TO_37", envir = globalenv())
                    else file.path(b, "data/ld_reference/hg38ToHg19.over.chain.gz")
  )
  out <- vapply(f, function(p) if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_,
                character(1))
  names(out) <- names(f)
  # The pair universe as coloc sees it: the columns it reads, in the order it reads them.
  .sugf <- file.path(b, file.path("results", paste0("repro_", PHASE1_TAG), "mr_suggestive.rds"))
  out <- c(out, mr_suggestive = if (!file.exists(.sugf)) NA_character_ else {
    .s <- data.table::as.data.table(readRDS(.sugf))
    # `with = FALSE`, not the `..` prefix: data.table's `..keep` looks for a variable called `keep`,
    # so naming it `.keep` made the lookup fail at run time rather than at parse time.
    keepcols <- intersect(c("gene", "outcome_id", "outcome_name", "method"), names(.s))
    digest::digest(paste(do.call(paste, c(as.list(.s[, keepcols, with = FALSE]), sep = "\x1f")),
                         collapse = "\n"),
                   algo = "md5", serialize = FALSE)
  })
  # Every piece of software whose version can change a number. coloc computes the posteriors;
  # data.table does the joins; the engine shells out to Python for region extraction and to
  # pyliftover for coordinate conversion, so pandas and pyliftover decide which variants come back.
  # Recording R and coloc alone said "the software is pinned" while three of the five were not.
  .pv <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  py <- tryCatch({
    pybin <- if (exists("PYTHON_BIN", envir = globalenv())) get("PYTHON_BIN", envir = globalenv()) else "python"
    o <- suppressWarnings(system2(pybin, c("-c", shQuote(paste(
      "import sys;",
      "v=lambda m:__import__('importlib.metadata',fromlist=['x']).version(m);",
      "print(sys.version.split()[0], end=' ');",
      "print(v('pandas'), end=' ');",
      "print(v('pyliftover'))")))
      , stdout = TRUE, stderr = FALSE))
    if (length(o)) paste(o, collapse = " ") else NA_character_
  }, error = function(e) NA_character_)
  out <- c(out,
           coloc_pkg  = .pv("coloc"),
           data_table = .pv("data.table"),
           r_version  = paste(R.version$major, R.version$minor, sep = "."),
           python_stack = py)          # "<python> <pandas> <pyliftover>"
  out
}

# R/coloc.R holds both GWAS_META and resolve_ukbppp_tar(), so it is loaded ONCE, here, before any
# component that needs either. .coloc_key_pqtl() used to run first and call resolve_ukbppp_tar()
# before anything had sourced it: the tryCatch around that call turned the resulting error into NA
# for every gene, so the archives silently contributed nothing to the key while the field still
# read as though they had been fingerprinted.
.coloc_key_engine <- function() {
  if (exists("GWAS_META", envir = globalenv()) &&
      exists("resolve_ukbppp_tar", envir = globalenv())) return(invisible(TRUE))
  # COLOC_SOURCE_ONLY must be set in the GLOBAL environment, because that is where R/coloc.R:1127
  # looks for it. Set as a local, the engine could not see it and sourcing it to read one lookup
  # table launched a complete colocalisation run -- observed, not theorised.
  .had <- exists("COLOC_SOURCE_ONLY", envir = globalenv())
  .old <- if (.had) get("COLOC_SOURCE_ONLY", envir = globalenv()) else NULL
  assign("COLOC_SOURCE_ONLY", TRUE, envir = globalenv())
  on.exit({
    if (.had) assign("COLOC_SOURCE_ONLY", .old, envir = globalenv())
    else if (exists("COLOC_SOURCE_ONLY", envir = globalenv()))
      rm("COLOC_SOURCE_ONLY", envir = globalenv())
  }, add = TRUE)
  # Did OUR source() create it, not does it exist. An interactive session that already holds an
  # older coloc_results was diagnosed as having just executed main().
  .had_res <- exists("coloc_results", envir = globalenv())
  sys.source(file.path(BASE_DIR, "R/coloc.R"), envir = globalenv())
  if (!.had_res && exists("coloc_results", envir = globalenv()))
    stop("_coloc_inputs_key.R: sourcing R/coloc.R executed main(). The source-only guard did not ",
         "hold, so building this key would run the whole analysis. Source the engine yourself ",
         "with COLOC_SOURCE_ONLY <- TRUE before calling coloc_inputs_key().")
  invisible(TRUE)
}

# A sampling digest (.edge_digest, sixteen 1 MB windows) stood here and was DELETED 2026-08-01. It
# had no callers -- .coloc_key_stamped() uses a full tools::md5sum and .coloc_key_pqtl() hashes the
# whole member -- but it carried a live comment saying "it is still not a full content hash, and it
# is not presented as one". In a file whose entire purpose is to assert that the key covers the input
# bytes, the first thing a reader met was a self-declared weaker guarantee belonging to code nothing
# called. Dead code that documents a superseded promise is worse than dead code.

# ONLY THE ARCHIVES THIS RUN READS. Stamping all 2,940 meant a tar for a protein absent from the
# pair set invalidated a legitimate run, and hashing 1.6 TB was never an option; the analysis reads
# one archive per gene in its own pair set.
.coloc_key_pqtl <- function(genes = NULL) {
  d <- file.path(BASE_DIR, "data/ukb_ppp/european")
  if (!dir.exists(d)) return(c(pqtl_archives = NA_character_))
  if (is.null(genes)) {
    # No scope given: name the universe rather than pretend to identify it, and say so in the value.
    fs <- list.files(d, pattern = "\\.tar$")
    return(c(pqtl_archives = sprintf("UNSCOPED:%d files present", length(fs))))
  }
  genes <- sort(unique(as.character(genes)))
  if (!length(genes)) return(c(pqtl_archives = "<no genes in scope>"))
  # ARCHIVE IDENTITY, NOT ARCHIVE PAYLOAD -- and the blind spot is stated, not implied.
  #
  # This used to seek into each gene's tar, locate the member for that gene's chromosome, and hash
  # roughly 15 MB of payload. Measured: 195 s for a 332-gene scope, 585 s for 927 genes, of which 91%
  # was locating and reading members (47.8 s navigation, 109.0 s payload, 14.5 s digest). Across the
  # 23 call sites in a full parallel run -- 21 once, each of six workers twice at both ends, 23 once,
  # 05 seven times, 06 twice -- that is about 2 h 52 m and roughly 235 GB of reading, and 68% of it
  # came from the twelve worker calls, each over the FULL manifest scope rather than its own chunk.
  #
  # What that bought was detection of ONE failure the fields below miss: a member payload replaced by
  # same-size content with the mtime preserved. That is out of scope for this pipeline's operating
  # model -- one operator, one machine, steps launched deliberately, archives downloaded once from
  # Synapse and not written again. IT IS NOT ENFORCED, though: all 2,940 tars have the Windows
  # read-only flag unset, so their immutability is a convention. THIS KEY DOES NOT DETECT A
  # SAME-SIZE, SAME-MTIME PAYLOAD SUBSTITUTION. Anyone who needs that must hash the bytes.
  #
  # Nor could payload hashing ever establish that the download matches Synapse: the stored release
  # listing carries name, id and version but NO checksum and NO size, so there is no external record
  # to compare against. It only ever noticed LOCAL mutation after the fact.
  #
  # What is kept identifies the release and notices every ordinary way an archive changes -- a
  # re-download, a restore, a rebuild, a different panel or OID, a truncated copy:
  #   - the resolved basename, which carries the gene symbol, UniProt, OID, version and panel;
  #   - size and mtime, which move on any ordinary rewrite;
  #   - the hash of eur_folder_contents.json, which pins the Synapse id and version of all 2,940.
  # The historical defect this generation exists to remove -- GRCh37 boundaries on GRCh38 positions
  # -- would not have been caught by payload hashing either. The code, coordinate and instrument
  # hashes elsewhere in this key are what catch that class.
  .listing <- file.path(BASE_DIR, "data/ukb_ppp/eur_folder_contents.json")
  .lh <- if (file.exists(.listing)) unname(tools::md5sum(.listing)) else "<absent>"
  fs <- vapply(genes, function(g) {
    p <- tryCatch(resolve_ukbppp_tar(g), error = function(e) NULL)
    if (is.null(p) || !length(p)) NA_character_ else as.character(p)[1]
  }, character(1))
  # A gene whose archive cannot be resolved is a gene whose exposure data this key does not
  # describe. Left as NA it simply dropped out of the digest, so the key claimed to cover archives
  # it had never seen.
  if (anyNA(fs))
    stop("_coloc_inputs_key.R: no UKB-PPP archive resolves for ", sum(is.na(fs)), " gene(s) in the ",
         "pair set (", paste(utils::head(names(fs)[is.na(fs)], 5), collapse = ", "),
         "). Their exposure data cannot be fingerprinted, so the key would describe less than the ",
         "run reads.")
  .id <- sprintf("%s:%.0f:%s", basename(fs), as.numeric(file.size(fs)),
                 format(file.mtime(fs), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"))
  c(pqtl_archives = sprintf("%d archive(s) for %d gene(s), listing %s, %s",
                            sum(!is.na(fs)), length(genes), substr(.lh, 1, 12),
                            digest::digest(paste(sort(.id), collapse = "\n"),
                                           algo = "md5", serialize = FALSE)))
}

# The outcome GWAS are 0.2-0.8 GB each, so they are not hashed whole. They are identified by size
# plus a digest of their first and last 1 MB, which is strictly stronger than the size-and-mtime
# this used to record: a restored copy, a re-download or a rebuild all change those bytes, and a
# preserved timestamp no longer hides it. It is still not a full hash, and an edit confined to the
# middle of a file would pass -- that is stated rather than implied.
#
# And only the files this analysis READS. Stamping the whole directory recursively meant 40 files --
# frozen extracts, local MR extracts, .lock files, unrelated studies, opengwas_jwt.txt -- so an
# unrelated file appearing or a lock being written invalidated a legitimate run.
.coloc_key_stamped <- function(outcome_ids = NULL) {
  .coloc_key_engine()
  gm <- get("GWAS_META", envir = globalenv())
  # Restricted to the outcomes in scope. All eight were stamped although the non-eGFR analysis uses
  # seven, so replacing the one it does not read invalidated the run.
  ids <- if (length(outcome_ids)) sort(unique(as.character(outcome_ids))) else unname(unlist(OUTCOMES))
  fs <- vapply(ids, function(i) {
    m <- gm[[i]]
    if (is.null(m) || is.null(m$file)) NA_character_ else file.path(OUTCOME_DIR, m$file)
  }, character(1))
  names(fs) <- ids
  # FULL md5. The sampler exists for the 1.6 TB of pQTL archives; the seven outcome files total
  # 6.0 GB, which hashes in about 30 s -- so sampling 2-9% of them, and then treating equality of
  # that sample as proof two workers read the same GWAS, was a cost saving that bought nothing.
  # Five of them are BGZF, where a same-length block outside the sampled windows can be swapped
  # with every sampled byte and the total size unchanged.
  out <- vapply(fs, function(q) if (is.na(q) || !file.exists(q)) "<absent>" else
                  sprintf("%.0f:%s", file.size(q), unname(tools::md5sum(q))), character(1))
  names(out) <- paste0(names(fs), " (", basename(fs), ")")
  out[order(names(out))]
}

# The coloc priors are part of the answer, not part of the environment: the same data with a
# different p12 gives a different PP.H4.
#
# They are NOT configuration variables -- they are literal arguments in the coloc.abf() call in
# R/coloc.R -- so they are read back out of the engine's source. Looking them up as COLOC_P1 etc.
# would have recorded NA for all three, which reads like "the priors are pinned" while pinning
# nothing. A missing call is fatal here rather than NA for the same reason.
.coloc_key_priors <- function() {
  src <- paste(readLines(file.path(BASE_DIR, "R/coloc.R"), warn = FALSE), collapse = "\n")
  m <- regmatches(src, regexpr("coloc::coloc\\.abf\\([^)]*\\)", src))
  if (!length(m))
    stop("_coloc_inputs_key.R: no coloc::coloc.abf(...) call found in R/coloc.R, so the priors ",
         "this run uses cannot be established. The engine's shape changed; update this reader ",
         "rather than recording nothing.")
  pr <- vapply(c("p1", "p2", "p12"), function(a) {
    v <- regmatches(m, regexpr(paste0("\\b", a, "\\s*=\\s*[0-9eE.+-]+"), m))
    if (!length(v)) stop("_coloc_inputs_key.R: coloc.abf() call does not pass ", a)
    gsub("\\s+", "", v)
  }, character(1))
  # These two ARE configuration, and they decide which pairs are called strong.
  g <- function(n) if (exists(n, envir = globalenv())) {
    v <- get(n, envir = globalenv())
    if (is.numeric(v) && length(v) == 1L) format(v, digits = 17) else NA_character_
  } else NA_character_
  c(pr, pp4 = g("COLOC_PP4_MIN"), h4h3 = g("COLOC_H4H3_MIN"))
}

# `scope` is the work this key describes: the genes whose pQTL archives will be read and the
# outcome ids that will be looked up. 21 computes it from the pair set and records it in the
# manifest; 22 and 23 take it from there, so all three build the SAME key over the SAME data.
coloc_inputs_key <- function(scope = NULL) {
  .coloc_key_engine()
  h <- c(.coloc_key_hashed(), .coloc_key_pqtl(scope$genes))
  s <- .coloc_key_stamped(scope$outcome_ids); p <- .coloc_key_priors()
  detail <- list(hashed = h, stamped = s, priors = p)
  # Sorted names on every component, so the key describes the inputs and not the order a
  # particular R session happened to list them in.
  flat <- c(paste0("H:", names(h)[order(names(h))], "=", h[order(names(h))]),
            paste0("S:", s),
            paste0("P:", names(p)[order(names(p))], "=", p[order(names(p))]))
  list(key = digest::digest(paste(flat, collapse = "\n"), algo = "md5", serialize = FALSE),
       detail = detail)
}

# Names the first component that differs, rather than reporting only that two 32-character
# strings are not equal -- which is what the caller actually needs in order to act.
coloc_inputs_diff <- function(a, b) {
  out <- character(0)
  for (part in c("hashed", "stamped", "priors")) {
    x <- a[[part]]; y <- b[[part]]
    for (n in union(names(x), names(y))) {
      xv <- if (n %in% names(x)) x[[n]] else "<absent>"
      yv <- if (n %in% names(y)) y[[n]] else "<absent>"
      if (!identical(xv, yv))
        out <- c(out, sprintf("  %s/%s: %s -> %s", part, n,
                              substr(xv, 1, 44), substr(yv, 1, 44)))
    }
  }
  if (!length(out)) "  (no component differs; the two keys were built by different code)" else out
}
