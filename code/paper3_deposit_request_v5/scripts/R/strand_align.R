# ============================================================================
# strand_align_R() -- STRAND-AWARE LD SIGN ALIGNMENT for the SuSiE engines (15 / 26 / 16; 27 reuses 26).
#
# THE BUG THIS FIXES (journal-editor ruling 2026-07-17, root-cause fix 2(b)). The engines build z from the
# pQTL sumstats oriented to `ea_pqtl` (z = pqtl_beta/se) and R (the LD matrix) oriented to the local 1000G EUR
# `.bim` A1 allele (plink --keep-allele-order). They then ASSUME `.bim` A1 == ea_pqtl for every SNP and feed
# (z, R) to susie_rss with NO sign harmonisation. That assumption was only spot-checked in 3 of 48 regions and
# fails sparsely: 37/38 found e.g. ATP5IF1~eGFR rs74064859 has `.bim` A1 = C = oa_pqtl while ea_pqtl = T, so z
# and R are anti-aligned at that SNP -> the fit's LD is sign-wrong there.
#
# WHAT THIS DOES. Given the harmonized allele table (ea_pqtl/oa_pqtl) and the LD matrix's ALLELE-SUFFIXED
# rownames (ieugwasr returns `rsID_A1_A2`; the engines strip the suffix -- this uses it BEFORE it is stripped),
# it re-orients R to the pQTL effect allele so z and R agree, per SNP, STRAND-AWARE:
#   - aligned  (.bim A1 is ea_pqtl, directly or by reverse-complement)         -> no change.
#   - SWAPPED  (.bim A1 is oa_pqtl, directly or by reverse-complement)         -> flip that SNP's R row & column
#                                                                                 (double flip keeps diag = 1);
#                                                                                 equivalently re-orients R to ea_pqtl.
#   - palindromic ({ea,oa} is A/T or C/G)                                       -> LEFT UNCHANGED. Allele matching
#       cannot decide strand for palindromes; these are verified SEPARATELY by allele frequency (pQTL EAF vs
#       1000G EUR MAF; both UKB-PPP and 1000G are forward-strand, and the flagged high-|z| palindromes were
#       confirmed same-strand). Resolving them here would need frequency plumbing; out of scope for this fix.
#   - allele_mismatch (.bim pair matches the pQTL pair neither directly nor by reverse-complement)
#                                                                              -> DROP the SNP (a different
#       biallelic coding cannot be reliably oriented). In 37/38 this count was 0; handled defensively.
#
# NO-OP WHEN THE DATA IS ALREADY ALIGNED: with 12,377 of 12,378 drawn-edge SNPs already aligned (37), this
# changes only the sparse mismatched SNPs. It is the standard sumstat-to-LD sign harmonisation.
#
# INVARIANTS ASSERTED (this touches frozen core engines, so it fails loud, never silently mis-aligns):
#   - the LD raw rownames parse to the SAME rsIDs as common_ld$snp_id, in the SAME order;
#   - R is square, symmetric, unit-diagonal on input;
#   - after flipping, R stays symmetric with unit diagonal.
# ============================================================================

.STRAND_COMP <- c(A = "T", T = "A", C = "G", G = "C")

.strand_revcomp <- function(a) {
  a <- toupper(a)
  vapply(a, function(x) {
    if (is.na(x) || !grepl("^[ACGT]+$", x)) return(NA_character_)
    paste(rev(unname(.STRAND_COMP[strsplit(x, "")[[1]]])), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

.strand_is_palindromic <- function(ea, oa) {
  ea <- toupper(ea); oa <- toupper(oa)
  is_snp <- ea %in% names(.STRAND_COMP) & oa %in% names(.STRAND_COMP)
  is_snp & (!is.na(.STRAND_COMP[ea]) & .STRAND_COMP[ea] == oa)
}

# Decide per SNP: "aligned" | "SWAPPED" | "palindromic" | "allele_mismatch". Vectorised.
# DATA-INTEGRITY FIRST (matches 37's order, cold review High): a pair that matches the pQTL alleles NEITHER
# directly NOR by reverse-complement is `allele_mismatch` even if it happens to be A/T or C/G -- palindromic is
# only assigned when the alleles actually correspond but the strand is undecidable.
.strand_verdict <- function(ea, oa, a1, a2) {
  ea <- toupper(ea); oa <- toupper(oa); a1 <- toupper(a1); a2 <- toupper(a2)
  ce <- .strand_revcomp(ea); co <- .strand_revcomp(oa)
  pal <- .strand_is_palindromic(ea, oa)
  same_pair <- mapply(function(x, y, p, q) setequal(c(x, y), c(p, q)), a1, a2, ea, oa)
  comp_pair <- !is.na(ce) & !is.na(co) & mapply(function(x, y, p, q) setequal(c(x, y), c(p, q)), a1, a2, ce, co)
  matched   <- same_pair | comp_pair
  out <- rep("allele_mismatch", length(ea))                   # default: bim pair corresponds to neither
  out[matched & pal] <- "palindromic"                         # alleles correspond but A/T,C/G -> strand undecidable, left as-is
  np <- matched & !pal
  out[np & same_pair & (a1 == ea)] <- "aligned"
  out[np & same_pair & (a1 == oa)] <- "SWAPPED"
  out[np & !same_pair & comp_pair & (a1 == ce)] <- "aligned"
  out[np & !same_pair & comp_pair & (a1 == co)] <- "SWAPPED"
  out
}

# common_ld : data.table/frame with snp_id, ea_pqtl, oa_pqtl, in the SAME row order as R's rows.
# ld_snps_raw_kept : character vector of the ALLELE-SUFFIXED LD rownames (`rsID_A1_A2`) for the kept SNPs,
#                    SAME order as common_ld / R.
# R : the LD correlation matrix (ld_mat_f), rownames = plain rsID, same order.
# label : a short string for log messages (edge/pair identity).
# Returns list(common_ld, R, log) with SWAPPED SNPs' R rows/cols flipped and allele_mismatch SNPs dropped.
strand_align_R <- function(common_ld, ld_snps_raw_kept, R, label = "") {
  stopifnot(is.matrix(R) || inherits(R, "Matrix"))
  if (!data.table::is.data.table(common_ld))   # `common_ld[keep]` below is data.table ROW-subset; a data.frame would column-subset (cold re-review Low)
    stop(sprintf("strand_align_R[%s]: common_ld must be a data.table (row-subset semantics assumed at the drop step)", label))
  R <- as.matrix(R)
  n <- nrow(R)
  if (ncol(R) != n) stop(sprintf("strand_align_R[%s]: R is not square (%dx%d)", label, n, ncol(R)))
  if (length(ld_snps_raw_kept) != n || nrow(common_ld) != n)
    stop(sprintf("strand_align_R[%s]: length mismatch -- R=%d, raw=%d, common_ld=%d",
                 label, n, length(ld_snps_raw_kept), nrow(common_ld)))
  # input R must already be a symmetric unit-diagonal correlation matrix
  .asym <- max(abs(R - t(R))); .dg <- max(abs(diag(R) - 1))
  if (!is.finite(.asym) || !is.finite(.dg) || .asym > 1e-8 || .dg > 1e-8)
    stop(sprintf("strand_align_R[%s]: input R is not symmetric unit-diagonal (max|R-Rt|=%.2e, max|diag-1|=%.2e)", label, .asym, .dg))

  # parse .bim A1/A2 from the allele-suffixed raw rownames `rsID_A1_A2`; verify the rsID matches common_ld order
  a1  <- toupper(sub("^.*_([A-Za-z]+)_([A-Za-z]+)$", "\\1", ld_snps_raw_kept))
  a2  <- toupper(sub("^.*_([A-Za-z]+)_([A-Za-z]+)$", "\\2", ld_snps_raw_kept))
  rs  <- sub("_[A-Za-z]+_[A-Za-z]+$", "", ld_snps_raw_kept)
  # Non-parsing rownames are a contract break and must stop. (Multi-char/indel A1/A2 DO parse here but are
  # dropped at the verdict step below as non-SNP, since this is a SNP-only sign harmonisation.)
  if (any(!grepl("_[A-Za-z]+_[A-Za-z]+$", ld_snps_raw_kept)))
    stop(sprintf("strand_align_R[%s]: %d LD rowname(s) lack the expected `rsID_A1_A2` allele suffix -- cannot recover .bim A1",
                 label, sum(!grepl("_[A-Za-z]+_[A-Za-z]+$", ld_snps_raw_kept))))
  # duplicate rsIDs would let match()/identical() pass while z and R carry a duplicated variable (cold re-review Medium)
  if (anyDuplicated(rs) || anyDuplicated(as.character(common_ld$snp_id)))
    stop(sprintf("strand_align_R[%s]: duplicate rsID(s) in LD/common_ld -- order/identity assertions cannot guarantee a 1:1 SNP mapping", label))
  if (!identical(rs, as.character(common_ld$snp_id)))
    stop(sprintf("strand_align_R[%s]: parsed LD rsIDs do not match common_ld$snp_id in order -- refusing to align on a mis-ordered table", label))
  if (!identical(rownames(R), rs) || !identical(colnames(R), rs))   # the row/col we flip must BE these SNPs (cold re-review Low)
    stop(sprintf("strand_align_R[%s]: R row/col names do not equal the parsed rsIDs in order -- refusing to flip a permuted matrix", label))

  v <- .strand_verdict(common_ld$ea_pqtl, common_ld$oa_pqtl, a1, a2)

  # SNP-ONLY harmonisation: a SNP whose .bim OR pQTL alleles are not a single A/C/G/T base cannot be strand-
  # decided -> force to drop (closes the unused-`bad` gap; cold re-review Medium). NA/multi-char -> not single.
  .single <- function(x) toupper(x) %in% c("A", "C", "G", "T")
  nonsnp  <- !(.single(a1) & .single(a2) & .single(common_ld$ea_pqtl) & .single(common_ld$oa_pqtl))

  # PALINDROMES ARE LEFT UNCHANGED here (PI decision 2026-07-17, simpler path). Allele matching cannot decide
  # strand for A/T,C/G SNPs; UKB-PPP and 1000G are both forward-strand (so same-strand in practice, confirmed
  # for the flagged high-|z| palindromes), and every palindromic fit-SNP is verified by a SEPARATE mandatory,
  # comprehensive allele-frequency GATE AFTER the re-run (39_palindromic_freq_gate.R), rather than by threading
  # frequency data through the core engines.
  flip <- which(v == "SWAPPED" & !nonsnp)          # non-palindromic single-base sign flips only
  drop <- which(v == "allele_mismatch" | nonsnp)   # allele_mismatch + non-SNP dropped; palindromic + aligned kept

  # SENSITIVITY drop-hook (guarded, INERT by default). When options(ckm.palindrome_drop = data.frame(label,rsid))
  # is set, ALSO drop the listed rsIDs FOR THIS label -- the drop-palindrome robustness re-fit (editor DECISION 2)
  # removes frequency-UNRESOLVABLE palindromes to show the SuSiE-tier fits are unchanged. No effect in normal runs.
  .dropopt <- getOption("ckm.palindrome_drop", NULL)
  if (is.data.frame(.dropopt) && nrow(.dropopt) && all(c("label", "rsid") %in% names(.dropopt))) {
    .dr_rs <- as.character(.dropopt$rsid[as.character(.dropopt$label) == label])
    if (length(.dr_rs)) drop <- union(drop, which(rs %in% .dr_rs))
  }

  if (length(flip)) {                       # re-orient R to ea_pqtl: flip the SWAPPED SNPs' rows AND columns
    R[flip, ] <- -R[flip, ]
    R[, flip] <- -R[, flip]                 # double flip => R[i,i] returns to +1; R[i,j] both-swapped unchanged
  }
  keep <- setdiff(seq_len(n), drop)
  if (length(drop)) { R <- R[keep, keep, drop = FALSE]; common_ld <- common_ld[keep] }

  if (nrow(R) == 0L)   # every SNP dropped -> return empty so the caller's `<10` gate skips, not a 0x0 postcondition crash (cold re-review Low)
    return(list(common_ld = common_ld, R = R,
                log = list(n = n, aligned = 0L, swapped = length(flip), palindromic = 0L,
                           allele_mismatch = sum(v == "allele_mismatch" & !nonsnp), non_snp = sum(nonsnp),
                           dropped = length(drop), swapped_rsids = rs[flip], dropped_rsids = rs[drop],
                           palindromic_rsids = character(0))))

  # post-condition: still a valid correlation matrix
  .asym2 <- max(abs(R - t(R))); .dg2 <- max(abs(diag(R) - 1))
  if (!is.finite(.asym2) || !is.finite(.dg2) || .asym2 > 1e-8 || .dg2 > 1e-8)
    stop(sprintf("strand_align_R[%s]: after flipping, R is not symmetric unit-diagonal (max|R-Rt|=%.2e, max|diag-1|=%.2e)", label, .asym2, .dg2))

  # AUDIT SINK for 39_palindromic_freq_gate.R. Written from this single choke point (no per-engine plumbing) ONLY
  # when this locus will actually be FIT: >= 10 SNPs after drops, matching the callers' `< 10` skip, so the sink
  # reflects exactly the fitted set -- a palindrome on a locus the caller then skips is NOT recorded (cold re-review
  # Medium). ENV VAR ONLY (single mechanism, so engine force-fresh, PSOCK propagation and 39 all key on the same
  # signal; cold re-review Medium). Unique file per call so multi-core workers (26/31) never race. Uses FINAL rows.
  .sink <- Sys.getenv("CKM_PALINDROME_SINK", "")
  if (nzchar(.sink) && nrow(R) >= 10L) {
    dir.create(.sink, showWarnings = FALSE, recursive = TRUE)
    .rid  <- Sys.getenv("CKM_RUN_ID", ""); .safe <- gsub("[^A-Za-z0-9]+", "-", label)
    .vk <- v[keep]; .nk <- nonsnp[keep]; .rk <- rs[keep]; .ak <- a1[keep]   # align verdict/allele vectors to the FINAL (post-drop) common_ld
    pidx <- which(.vk == "palindromic" & !.nk)
    # SEEN marker on EVERY fitted call (even 0 palindromes): proof capture ran, so 39 can distinguish a true
    # zero-palindrome run from capture-off and reconcile counts. n_palindromic must equal the palin rows below.
    data.table::fwrite(data.table::data.table(label = label, run_id = .rid, n_palindromic = length(pidx)),
                       tempfile(pattern = sprintf("seen_%s_%d_", .safe, Sys.getpid()), tmpdir = .sink, fileext = ".csv"))
    if (length(pidx)) {
      # BOTH z's: a palindrome enters the exposure AND outcome SuSiE fit on the SAME (possibly mis-signed) R, so 39
      # weighs severity by max(|z_exp|,|z_out|). out_beta is ea_pqtl-oriented after harmonize, so z_out is its outcome-fit z.
      .zx <- if (all(c("pqtl_beta", "pqtl_varbeta") %in% names(common_ld))) common_ld$pqtl_beta[pidx] / sqrt(common_ld$pqtl_varbeta[pidx]) else NA_real_
      .zo <- if (all(c("out_beta", "out_varbeta") %in% names(common_ld))) common_ld$out_beta[pidx] / sqrt(common_ld$out_varbeta[pidx]) else NA_real_
      data.table::fwrite(data.table::data.table(
        label = label, rsid = .rk[pidx],
        ea_pqtl = toupper(common_ld$ea_pqtl)[pidx], oa_pqtl = toupper(common_ld$oa_pqtl)[pidx],
        bim_a1 = .ak[pidx], z_pqtl = .zx, z_out = .zo, run_id = .rid),
        tempfile(pattern = sprintf("palin_%s_%d_", .safe, Sys.getpid()), tmpdir = .sink, fileext = ".csv"))
    }
  }

  log <- list(n = n, aligned = sum(v == "aligned" & !nonsnp), swapped = length(flip),
              palindromic = sum(v == "palindromic" & !nonsnp),
              allele_mismatch = sum(v == "allele_mismatch" & !nonsnp), non_snp = sum(nonsnp),
              dropped = length(drop), swapped_rsids = rs[flip], dropped_rsids = rs[drop],
              palindromic_rsids = rs[v == "palindromic" & !nonsnp])
  if (length(flip) || length(drop))
    cat(sprintf("    [strand_align %s] n=%d aligned=%d SWAPPED(flipped)=%d palindromic(kept)=%d allele_mismatch=%d non_snp=%d dropped=%d\n",
                label, log$n, log$aligned, log$swapped, log$palindromic, log$allele_mismatch, log$non_snp, log$dropped))
  list(common_ld = common_ld, R = R, log = log)
}
