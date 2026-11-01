# w4b_04_credsets.R -- Stage 1 W4b, step 4: section 4(v) credible sets by methods (1)-(5), gene membership under the three
# Rule 3(iii) annotations, and Rule 3(ii) robustness of methods (3)-(5) (489 agreement; 5,000-draw 503 bootstrap, seed 20261007).
# Method definitions (redacted pre-declaration 4(v); method reading of seen script prespec_review3/cs.R logged in run.log):
#  (1) SuSiE-RSS with the 503 reference (SNP set p1np): the T2D (FinnGen) credible set at signal X and the ARG1 protein
#      credible set at signal X (both are read; the effector rule needs both inside the gene).
#  (2) UKB-PPP in-sample SuSiE (Sun 2023 ST16): the protein credible set of the signal X.
#  (3) shared P-T2D(FinnGen) set at X: protein and T2D each z-space-conditioned on the other signal (single SNP), Wakefield
#      ABF (W 0.15 quant sdY 1; W 0.2 cc), 95% set of SNP.PP.H4 (prop. to exp(lABF_P + lABF_D)) = coloc.abf's SNP.PP.H4.
#  (4) T2D(FinnGen) X|other set: T2D conditioned on the other signal, 95% set of exp(lABF_D).
#  (5) as (3) with DIAMANTE 2018 EUR noUKB (printed beta/SE).
# Conditioning drops SNPs with r^2 > 0.9 to the conditioning SNP (l3fm__conditional COLLIN) and the SNP itself.
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets/w4b_common.R"))
suppressMessages({library(coloc); library(readxl)})
logf("==== W4b step 4 (credible sets) start")
D <- readRDS(file.path(OUT, "data", "w4b_data.rds")); H <- D$H
acc <- new_rows()
SIG <- c(A = "6:131897278:T:C", B = "6:131894358:T:C")
OTHER <- c(A = "B", B = "A")

# ---------------------------------------------------------------- annotations (glist semantics; three annotations)
gt <- fread(rpath("GRID_ens116_genes"), colClasses = list(character = c("gene_id", "symbol", "biotype", "chr")))
gt <- unique(gt[chr == "6", .(gene_id, symbol, biotype, start = start38, end = end38)])
stopifnot(!anyDuplicated(gt$gene_id))
gt[symbol == "" | is.na(symbol), symbol := gene_id]
GT <- list(ens116pc = gt[biotype == "protein_coding"], ens116all = gt)
gl <- fread(xpath("GLIST_HG19"), header = FALSE, col.names = c("chr", "start", "end", "symbol"), colClasses = c("character", "integer", "integer", "character"))
GT$glist <- gl[chr == "6"]
logf("annotations: ens116pc %d genes, ens116all %d, glist chr6 %d", nrow(GT$ens116pc), nrow(GT$ens116all), nrow(GT$glist))
ANN <- names(GT)
annot_snp <- function(a, snps) {
  h <- H[match(snps, snp)]; pos <- if (a == "glist") h$pos37 else h$pos38
  lapply(pos, function(p) genes_at(GT[[a]], p))
}
SNPANN <- setNames(lapply(ANN, function(a) setNames(annot_snp(a, H$snp), H$snp)), ANN)
# a set's genes: union of in-body genes, intergenic members as ~nearest; 'within' = genes whose body holds EVERY member
set_genes <- function(members, a) {
  g <- SNPANN[[a]][members]
  lab <- unique(unlist(lapply(g, function(x) if (x$dist == 0) x$genes else paste0("~", x$genes))))
  inb <- lapply(g, function(x) if (x$dist == 0) x$genes else character(0))
  within <- Reduce(intersect, inb)
  list(genes = paste(sort(lab), collapse = ";"), within = paste(sort(within), collapse = ";"))
}
# check the annotation against the engine's own lookup on a known case: rs2781668 lies in ARG1 under all three
for (a in ANN) logf("annotation check %s: rs2781668 -> %s (dist %d)", a, paste(SNPANN[[a]][[SIG[["A"]]]]$genes, collapse = ","), SNPANN[[a]][[SIG[["A"]]]]$dist)

# ---------------------------------------------------------------- ABF helpers (checked against coloc.abf below)
labf <- function(z, V, W) { r <- W^2 / (W^2 + V); 0.5 * (log(1 - r) + r * z^2) }
cs95 <- function(l, snps) {
  pp <- exp(l - max(l)); pp <- pp / sum(pp); o <- order(-pp)
  n <- which(cumsum(pp[o]) >= 0.95)[1]
  list(members = snps[o[seq_len(n)]], top = snps[o[1]], top_pp = pp[o[1]])
}
cond1 <- function(z, se, r, j) {  # single-SNP z-space conditioning on SNP j given r = cor(all SNPs, SNP j)
  q <- r^2; zc <- (z - r * z[j]) / sqrt(pmax(1 - q, 1e-12)); sc <- se / sqrt(pmax(1 - q, 1e-12))
  bad <- q > 0.9; zc[bad] <- NA; sc[bad] <- NA
  list(z = zc, se = sc)
}
METH <- list(m3 = list(trait = "FG", shared = TRUE), m4 = list(trait = "FG", shared = FALSE), m5 = list(trait = "DM", shared = TRUE))
WCC <- 0.2; WQ <- 0.15
# one method at signal X, given the vector r (all set SNPs vs the conditioning SNP)
run_cs <- function(m, X, snps, r_cond, condsnp) {
  h <- H[match(snps, snp)]
  zD <- h[[paste0("z__", METH[[m]]$trait)]]; sD <- h[[paste0("se__", METH[[m]]$trait)]]
  rr <- r_cond; j <- match(condsnp, snps); rr[j] <- 1
  cD <- cond1(zD, sD, rr, j)
  if (METH[[m]]$shared) {
    cP <- cond1(h$z__P__ARG1, h$se__P__ARG1, rr, j)
    ok <- !is.na(cD$z) & !is.na(cP$z)
    l <- labf(cD$z[ok], cD$se[ok]^2, WCC) + labf(cP$z[ok], cP$se[ok]^2, WQ)
  } else {
    ok <- !is.na(cD$z)
    l <- labf(cD$z[ok], cD$se[ok]^2, WCC)
  }
  c(cs95(l, snps[ok]), list(nsnp = sum(ok), X_in = X %in% snps[ok]))
}
SETSN <- list(m3 = "FG", m4 = "FG", m5 = "DM")
set_for <- function(m, setname) { s <- D$SETS[[setname]]; s[!is.na(H[match(s, snp)][[paste0("z__", SETSN[[m]])]])] }

# ---------------------------------------------------------------- full panels (503 = p1np, 489 = p2np)
FULL <- list()
for (pn in c("503", "489")) {
  setname <- if (pn == "503") "p1np" else "p2np"; Rf <- if (pn == "503") D$R503 else D$R489
  for (m in names(METH)) for (X in names(SIG)) {
    snps <- set_for(m, setname); cs_ <- SIG[[OTHER[[X]]]]
    stopifnot(SIG[[X]] %in% snps, cs_ %in% snps)
    r <- run_cs(m, SIG[[X]], snps, Rf[snps, cs_], cs_)
    FULL[[paste(m, X, pn)]] <- r
  }
}
# check my ABF against R coloc.abf (method 3, signal A, 503): SNP.PP.H4 vectors must agree
{
  snps <- set_for("m3", "p1np"); rr <- D$R503[snps, SIG[["B"]]]; h <- H[match(snps, snp)]
  jB <- match(SIG[["B"]], snps); cD <- cond1(h$z__FG, h$se__FG, rr, jB); cP <- cond1(h$z__P__ARG1, h$se__P__ARG1, rr, jB); ok <- !is.na(cD$z) & !is.na(cP$z)
  ab <- suppressWarnings(suppressMessages(coloc.abf(list(beta = cP$z[ok] * cP$se[ok], varbeta = cP$se[ok]^2, type = "quant", sdY = 1, snp = snps[ok]),
                                                    list(beta = cD$z[ok] * cD$se[ok], varbeta = cD$se[ok]^2, type = "cc", s = D$CASEFRAC$FG, snp = snps[ok]),
                                                    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)))
  l <- labf(cD$z[ok], cD$se[ok]^2, WCC) + labf(cP$z[ok], cP$se[ok]^2, WQ); pp <- exp(l - max(l)); pp <- pp / sum(pp)
  dd <- max(abs(pp - ab$results$SNP.PP.H4[match(snps[ok], ab$results$snp)]))
  logf("CHECK own ABF vs coloc.abf SNP.PP.H4 (method 3, A, 503): max |diff| %.3g; conditioned PP.H4 %.4f", dd, ab$summary[["PP.H4.abf"]])
  stopifnot(dd < 1e-8)
}

# ---------------------------------------------------------------- bootstrap (5,000 draws of the 503 people; seed 20261007)
bfile <- file.path(OUT, "data", "boot_index_seed20261007.rds")
if (!file.exists(bfile)) { set.seed(SEED); IDX <- vapply(seq_len(NBOOT), function(i) sample.int(503L, 503L, replace = TRUE), integer(503)); saveRDS(IDX, bfile) }
IDX <- readRDS(bfile); stopifnot(dim(IDX) == c(503L, NBOOT))
logf("bootstrap index: set.seed(%d); draw d = sample.int(503, 503, replace = TRUE), d = 1..%d (R %s); md5 %s", SEED, NBOOT, R.version.string, md5(bfile))
gp <- rpath("R1_ARG1_REF503_I8")      # read in place, never copied
G <- matrix(readBin(gp, "integer", n = 503L * 3067L + 1L, size = 1L, signed = TRUE), nrow = 503, byrow = TRUE)
stopifnot(ncol(G) == nrow(H), all(abs(colMeans(G) / 2 - H$af_ref) < 1e-5)); colnames(G) <- H$snp
U <- unique(unlist(lapply(names(METH), function(m) set_for(m, "p1np"))))
G <- G[, U]; storage.mode(G) <- "double"
ck <- file.path(OUT, "boot"); dir.create(ck, showWarnings = FALSE)
CH <- 250L
full_lab <- function(m, X, pn, a) set_genes(FULL[[paste(m, X, pn)]]$members, a)$genes
for (c0 in seq(1L, NBOOT, by = CH)) {
  cf <- file.path(ck, sprintf("boot_%04d_%04d.rds", c0, c0 + CH - 1L))
  if (file.exists(cf)) next
  t0 <- Sys.time(); rows <- list(); nmono <- 0L
  for (d in c0:(c0 + CH - 1L)) {
    Gd <- G[IDX[, d], , drop = FALSE]
    vv <- colSums(Gd^2) - colSums(Gd)^2 / nrow(Gd); nmono <- nmono + sum(vv <= 1e-9)
    rAB <- suppressWarnings(cor(Gd, Gd[, SIG]))
    rAB[is.na(rAB)] <- 0                                   # monomorphic in the draw -> r = 0
    for (m in names(METH)) for (X in names(SIG)) {
      snps <- set_for(m, "p1np"); cs_ <- SIG[[OTHER[[X]]]]
      r <- run_cs(m, SIG[[X]], snps, rAB[snps, cs_], cs_)
      gl_ <- vapply(ANN, function(a) set_genes(r$members, a)$genes, "")
      rows[[length(rows) + 1]] <- data.table(draw = d, method = m, signal = X, size = length(r$members), top = r$top,
                                             same_members = setequal(r$members, FULL[[paste(m, X, "503")]]$members),
                                             genes_ens116pc = gl_[["ens116pc"]], genes_ens116all = gl_[["ens116all"]], genes_glist = gl_[["glist"]],
                                             r_AB = rAB[SIG[["A"]], "6:131894358:T:C"])
    }
  }
  bt <- rbindlist(rows); saveRDS(bt, paste0(cf, ".part")); file.rename(paste0(cf, ".part"), cf)
  logf("boot draws %d-%d done in %.0fs (monomorphic SNP-draws %d); done %d of %d", c0, c0 + CH - 1L, as.numeric(difftime(Sys.time(), t0, units = "secs")),
       nmono, c0 + CH - 1L, NBOOT)
}
BT <- rbindlist(lapply(sort(list.files(ck, "^boot_.*\\.rds$", full.names = TRUE)), readRDS))
stopifnot(nrow(BT) == NBOOT * length(METH) * length(SIG), !anyDuplicated(BT[, .(draw, method, signal)]), setequal(unique(BT$draw), seq_len(NBOOT)))
fwrite(BT, file.path(OUT, "boot", "boot_credsets_long.tsv.gz"), sep = "\t")
logf("bootstrap table: %d rows (%d draws x %d methods x %d signals); r(A,B) over draws: 2.5%% %.4f, 97.5%% %.4f", nrow(BT), NBOOT, length(METH), length(SIG),
     quantile(BT[method == "m3" & signal == "A"]$r_AB, 0.025), quantile(BT[method == "m3" & signal == "A"]$r_AB, 0.975))

# ---------------------------------------------------------------- method (1): SuSiE 503 (p1np) credible sets at A and B
cs_at <- function(tr, setname, X) {
  o <- readRDS(file.path(OUT, "susie", sprintf("%s__%s.rds", tr, setname)))
  if (is.null(o$cs) || !nrow(o$cs)) return(list(members = character(0), r2 = NA_real_, cs = "none", o = o))
  Rf <- if (D$PANEL[[setname]] == "503") D$R503 else D$R489
  r2 <- vapply(strsplit(o$cs$members, ";"), function(mm) max(Rf[mm, X]^2), 0)
  j <- which.max(r2)
  list(members = strsplit(o$cs$members[j], ";")[[1]], r2 = r2[j], cs = o$cs$cs[j], o = o, lead = o$cs$lead_snp[j])
}
# ---------------------------------------------------------------- method (2): UKB-PPP in-sample SuSiE (ST16)
st16 <- as.data.table(read_excel(rpath("EXP_SUN2023_SUPP_TABLES_ST16"), sheet = "ST16", skip = 4)); names(st16)[1] <- "protein_id"
x6 <- st16[protein_id == "ARG1:P05089:OID21445:v1" & `Cis/trans` == "cis"]
ukb_cs <- function(X) {
  p37 <- H[snp == X]$pos37
  r <- x6[grepl(paste0("^6:", p37, ":"), `Variant ID`)]; stopifnot(nrow(r) == 1)
  ids <- strsplit(r$`CS variant IDs`, ", ")[[1]]; rs <- strsplit(r$`CS rsIDs`, ", ")[[1]]
  p <- as.integer(vapply(strsplit(ids, ":"), `[`, "", 2)); al <- lapply(strsplit(ids, ":"), `[`, 3:4)
  # map to the harmonised table (GRCh37 pos + alleles); a member outside it is located through the UKB-PPP window (GRCh37 -> GRCh38)
  hh <- H[match(p, pos37)]
  stopifnot(all(!is.na(hh$snp) & mapply(function(a, x, y) setequal(a, c(x, y)), al, hh$a1, hh$a2)))
  list(members = hh$snp, rsids = rs)
}

# ---------------------------------------------------------------- assemble per method x signal
TAB <- list()
mk <- function(method, X, members, admis, robust_note, extra = "") {
  ga <- lapply(ANN, function(a) set_genes(members, a)); names(ga) <- ANN
  TAB[[length(TAB) + 1]] <<- data.table(method = method, signal = X, size = length(members),
                                        members = paste(H[match(members, snp)]$rsid_ref, collapse = ";"),
                                        genes_ens116pc = ga$ens116pc$genes, genes_ens116all = ga$ens116all$genes, genes_glist = ga$glist$genes,
                                        within_ens116pc = ga$ens116pc$within, within_ens116all = ga$ens116all$within, within_glist = ga$glist$within,
                                        admissibility = admis, robustness = robust_note, extra = extra)
}
S <- fread(file.path(OUT, "susie", "susie_s_table.tsv"))
# Stage 2 step 2: s_max over all seven distinct fits (the four fits + the three of susie/w4b_02b_fitset7.R)
S <- rbind(S[, .(trait, set, panel, s)], fread(file.path(OUT, "susie", "susie_s_extra.tsv"))[, .(trait, set, panel, s)])
s_pair <- max(S[trait %in% c("P__ARG1", "T2D_FG")]$s)
adm1 <- if (s_pair < 0.10) "admissible (Rule 3(i): max s < 0.10)" else sprintf("inadmissible (Rule 3(i): max s over every ARG1 and T2D fit = %.4f >= 0.10)", s_pair)
ROB <- list()
for (X in names(SIG)) {
  dcs <- cs_at("T2D_FG", "p1np", SIG[[X]]); pcs <- cs_at("P__ARG1", "p1np", SIG[[X]])
  mk("1_T2D", X, if (!is.na(dcs$r2) && dcs$r2 >= 0.8) dcs$members else character(0), adm1, "Rule 3(i) governs; no bootstrap",
     sprintf("T2D FinnGen SuSiE CS %s (lead %s) max r2 with signal %.3f; 'at signal' requires r2 >= 0.8", dcs$cs, H[match(dcs$lead, snp)]$rsid_ref, dcs$r2))
  mk("1_protein", X, if (!is.na(pcs$r2) && pcs$r2 >= 0.8) pcs$members else character(0), adm1, "Rule 3(i) governs; no bootstrap",
     sprintf("ARG1 SuSiE CS %s (lead %s) max r2 with signal %.3f", pcs$cs, H[match(pcs$lead, snp)]$rsid_ref, pcs$r2))
  u <- ukb_cs(SIG[[X]]); mk("2", X, u$members, "admissible (in-sample LD)", "not applicable (in-sample)", paste("ST16 CS rsIDs:", paste(u$rsids, collapse = ",")))
  for (m in names(METH)) {
    f5 <- FULL[[paste(m, X, "503")]]; f4 <- FULL[[paste(m, X, "489")]]
    b <- BT[method == m & signal == X]
    same489 <- vapply(ANN, function(a) set_genes(f5$members, a)$genes == set_genes(f4$members, a)$genes, TRUE)
    bshare <- vapply(ANN, function(a) mean(b[[paste0("genes_", a)]] == set_genes(f5$members, a)$genes), 0)
    ball <- mean(b$genes_ens116pc == set_genes(f5$members, "ens116pc")$genes & b$genes_ens116all == set_genes(f5$members, "ens116all")$genes &
                 b$genes_glist == set_genes(f5$members, "glist")$genes)
    robust <- all(same489) && ball >= 0.90
    ROB[[paste(m, X)]] <- list(same489 = same489, bshare = bshare, ball = ball, robust = robust, smem = mean(b$same_members),
                               f4 = f4, f5 = f5)
    mk(sub("m", "", m), X, f5$members, "Rule 3(ii) applies", sprintf("%s: 489 same genes (all 3 annotations) %s; boot share same genes (all 3) %.4f",
                                                                     if (robust) "robust" else "NOT robust (sensitive to the reference LD)", all(same489), ball),
       sprintf("top %s (pp %.3f); 489 set: %s (%s)", H[match(f5$top, snp)]$rsid_ref, f5$top_pp, paste(H[match(f4$members, snp)]$rsid_ref, collapse = ";"),
               set_genes(f4$members, "ens116pc")$genes))
  }
}
TAB <- rbindlist(TAB); fwrite(TAB, file.path(OUT, "credible_sets_by_method.tsv"), sep = "\t")
for (i in seq_len(nrow(TAB))) logf("CS %-9s %s size %3d genes pc=%s all=%s glist=%s within pc=%s all=%s glist=%s | %s | %s", TAB$method[i], TAB$signal[i], TAB$size[i],
                                   TAB$genes_ens116pc[i], TAB$genes_ens116all[i], TAB$genes_glist[i], TAB$within_ens116pc[i], TAB$within_ens116all[i],
                                   TAB$within_glist[i], TAB$robustness[i], TAB$extra[i])

# ---------------------------------------------------------------- results rows
MLAB <- c("1_T2D" = "(1) SuSiE-RSS, 503 reference (p1np), T2D FinnGen credible set at the signal",
          "1_protein" = "(1) SuSiE-RSS, 503 reference (p1np), ARG1 protein credible set at the signal",
          "2" = "(2) UKB-PPP in-sample SuSiE (Sun 2023 ST16), protein credible set of the signal",
          "3" = "(3) shared ARG1-T2D(FinnGen) 95% set, both conditioned on the other signal (z-space, 503 LD)",
          "4" = "(4) T2D(FinnGen) 95% ABF set conditioned on the other signal (z-space, 503 LD)",
          "5" = "(5) shared ARG1-T2D(DIAMANTE) 95% set, both conditioned on the other signal (z-space, 503 LD)")
MRES <- c("1_T2D" = "FinnGen R10 T2D", "1_protein" = "UKB-PPP OID21445", "2" = "Sun 2023 ST16", "3" = "UKB-PPP OID21445 x FinnGen R10 T2D",
          "4" = "FinnGen R10 T2D", "5" = "UKB-PPP OID21445 x DIAMANTE noUKB")
for (i in seq_len(nrow(TAB))) {
  mth <- TAB$method[i]; X <- TAB$signal[i]
  base <- sprintf("S4.cs.%s.%s", mth, X)
  pnl <- if (mth == "2") "UKB-PPP in-sample" else "503"
  acc <- add_row(acc, paste0(base, ".members"), "4(v)", MLAB[[mth]], MRES[[mth]], pnl, value = if (TAB$size[i]) TAB$members[i] else "NO CREDIBLE SET AT THE SIGNAL",
                 n = TAB$size[i], unit = "rsIDs", notes = paste(TAB$admissibility[i], "|", TAB$robustness[i], "|", TAB$extra[i]))
  for (a in ANN) {
    acc <- add_row(acc, paste0(base, ".genes_", a), "4(v); Rule 3(iii) annotations", paste(MLAB[[mth]], "- genes holding a member (~ = nearest gene of an intergenic member)"),
                   MRES[[mth]], pnl, value = if (TAB$size[i]) TAB[[paste0("genes_", a)]][i] else "NO CREDIBLE SET AT THE SIGNAL", unit = "gene symbols", notes = sprintf("within (every member in the gene body): %s", TAB[[paste0("within_", a)]][i]))
  }
}
for (k in names(ROB)) {
  m <- sub(" .*", "", k); X <- sub(".* ", "", k); r <- ROB[[k]]
  base <- if (X == "A") sprintf("S4.cs.%s.robust", sub("m", "", m)) else sprintf("S4.cs.%s.B.robust", sub("m", "", m))
  acc <- add_row(acc, paste0(base, ".panel489_same"), "4(v); Rule 3(ii)", sprintf("%s at signal %s: 489 (p2np) set spans the same genes under all three annotations", MLAB[[sub("m", "", m)]], X),
                 MRES[[sub("m", "", m)]], "489", value = all(r$same489), unit = "logical",
                 notes = sprintf("per annotation pc/all/glist: %s; 489 members %s", paste(r$same489, collapse = "/"), paste(H[match(r$f4$members, snp)]$rsid_ref, collapse = ";")))
  acc <- add_row(acc, paste0(base, ".boot_share_same_genes"), "4(v); Rule 3(ii)", sprintf("%s at signal %s: share of 5,000 503 bootstrap draws whose set spans the same genes under all three annotations", MLAB[[sub("m", "", m)]], X),
                 MRES[[sub("m", "", m)]], "503 bootstrap", value = r$ball, n = NBOOT, unit = "share",
                 notes = sprintf("per annotation pc/all/glist: %s; identical membership share %.4f; seed %d (R); robust (489 same AND share >= 0.90): %s",
                                 paste(sprintf("%.4f", r$bshare), collapse = "/"), r$smem, SEED, r$robust))
}
write_part(acc, "04_credsets")
saveRDS(list(TAB = TAB, ROB = ROB, FULL = FULL, s_pair = s_pair), file.path(OUT, "data", "w4b_credsets.rds"))
logf("==== W4b step 4 done")
