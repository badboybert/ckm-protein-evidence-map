# Population prevalence of each BINARY outcome, for cascade/09_steiger.R only.
#
# NOT in R/00_config.R deliberately. That file is sourced by every layer and its md5 is hashed into
# the L3 provenance sidecar (33_sensitivity_provenance.R lists it in EXTRA_LOGIC and enforces that
# every certified artifact POSTDATES its listed code). Adding an L4-only table there re-dated the
# shared config mid-run and would have failed 33's staleness gate against artifacts written earlier
# the same evening -- an L4 addition invalidating the L1-L3 provenance chain, for nothing.
#
# WHAT THIS IS. get_r_from_lor()'s `prevalence` argument is the prevalence of the disease IN THE
# POPULATION. The study's case:control ratio is conveyed separately through ncase/ncontrol. Until
# 2026-08-04, 09_steiger passed n_case/n_total -- the GWAS sample case fraction -- so every binary
# outcome's r2_out, steiger_z and steiger_p was built on the wrong liability transformation.
# Case-control GWAS are case-enriched by design, so the two numbers are not interchangeable.
#
# DELIBERATELY EMPTY. Each entry is a per-disease epidemiological input and needs a citation, in the
# population the GWAS was drawn from: Finnish for the FinnGen outcomes (CAD, HF, T2D, NAFLD, AF),
# European for MEGASTROKE (Stroke) and Wuttke (CKD). 09_steiger stops with the outcome named when an
# entry is missing or outside (0,1), rather than proceeding on a substituted value -- a
# plausible-looking default here would reproduce the original defect while looking deliberate.
#
# Fill in as:  "FinnGen_R10_T2D" = 0.xx,   # source, population, year
OUTCOME_POP_PREVALENCE <- list(
  # ── FinnGen R10 endpoints ──────────────────────────────────────────────────────────────────────
  # Source: Risteys FinnGen R10 (https://r10.risteys.finngen.fi/endpoints/<ENDPOINT>), "unadjusted
  # prevalence (%)" over the whole FinnGen cohort. THE R10 PAGES, not the current R14 ones: the GWAS
  # summary statistics used here are R10, and R14 has four further releases of follow-up. Checked
  # rather than assumed -- CHD is 21.68% at R14 against 11.40% at R10, a 1.9-fold difference.
  #
  # This is the prevalence of the phenotype AS FinnGen ASCERTAINS IT (the ICD-coded endpoint), which
  # is the quantity get_population_allele_frequency() needs: it weights the case and control allele
  # frequencies to rebuild the source population, so K must mean the same thing as the case group it
  # is weighting. Community estimates of the disease are a different quantity -- Finnish surveys put
  # imaging-defined fatty liver near 21%, thirty times the register-diagnosed NAFLD endpoint.
  "FinnGen_R10_I9_CHD"       = 0.1140,   # Risteys R10 I9_CHD,       46,959 cases; our extract 11.39%
  "FinnGen_R10_I9_HEARTFAIL" = 0.0720,   # Risteys R10 I9_HEARTFAIL, 29,672 cases; our extract 7.20%
  "FinnGen_R10_T2D"          = 0.1580,   # Risteys R10 T2D,          65,083 cases; our extract 16.26%
  "FinnGen_R10_NAFLD"        = 0.0062,   # Risteys R10 NAFLD,         2,568 cases; our extract 0.62%
  # AF is where the distinction earns its keep. Risteys and our extract agree on the case count
  # (50,742 / 50,743) but not the denominator: Risteys divides by the whole cohort (412,181 -> 12.3%)
  # while the GWAS analysed 261,395 (-> 19.4%). The AF analysis is case-enriched relative to the
  # cohort, so here the sample fraction genuinely is not the population prevalence.
  "FinnGen_R10_I9_AF"        = 0.1232    # Risteys R10 I9_AF,        50,742 cases; our extract 19.41%
)

# ── Stroke and CKD: NO de-ascertainment, deliberately, and kept in a SEPARATE object ─────────────
# These two are NOT population prevalences and must not sit in the list above, where a reader or a
# later rebuild would be unable to tell which entries are registry figures and which are placeholders.
#
# Both are supplied as the analysed sample's own case fraction. That is not a fabricated number: when
# K equals ncase/(ncase+ncontrol), get_population_allele_frequency() returns the reported allele
# frequency unchanged (verified to 1.1e-16), so this is the algebraically exact statement "apply no
# de-ascertainment" -- a defined option, not a guess dressed as a citation.
#
# WHY NEITHER HAS A POPULATION PREVALENCE:
#   MEGASTROKE (MEGASTROKE_EUR_GCST006906) meta-analyses dozens of studies of heterogeneous design -- population
#   cohorts, hospital-based case-control, matched controls -- so the correct de-ascertainment weight is
#   per-study and the meta-analysed frequency is a weighted mixture of different case:control ratios.
#   A single scalar cannot unmix that; this is a definitional limit, not a gap in the literature. The
#   candidate values are also the wrong quantity: European point prevalence (~5%) counts surviving
#   stroke history, German 12-month prevalence (1.6%) counts a period, and lifetime risk (~1 in 4 over
#   25) is cumulative incidence -- none is the cross-sectional composition this weight requires.
#
#   CKDGen (Wuttke2019_CKD_EA) applies an eGFR < 60 threshold WITHIN each contributing cohort. Nothing
#   over-samples cases, so the sample case fraction already IS a direct estimate of the prevalence of
#   eGFR < 60 in the source population, at its actual age structure and actual creatinine equation.
#   The published European range (1% central Italy to 5.9% northeast Germany) would be worse here: it
#   describes a general population, not these cohorts.
#
# DERIVED, NEVER LITERAL. The identity above holds only when K equals the sample fraction EXACTLY; a
# hardcoded 0.0908 would silently become "a mysterious number near the sample fraction" the moment an
# extract's N changed, at which point neither reading would be true.
OUTCOME_NO_DEASCERTAINMENT <- c("MEGASTROKE_EUR_GCST006906", "Wuttke2019_CKD_EA")

# ── Study-wide case fraction, for the Fisher-z denominator only ─────────────────────────────────
# A THIRD quantity, kept apart from the two above because conflating it with either is exactly the
# mistake this file has now made twice. It is the study's OWN published cases / total participants.
#
# Not the extract-implied fraction. The extract's median per-variant N is smaller than the study
# total (CKD: 410,882 against 480,698), so dividing the study-wide case count by it mixes two
# granularities -- study-level numerator, variant-level denominator -- and inflates the fraction
# (0.1007 against the published 0.0861). That inflated value is correct for .prev, where it must
# match ncase/(ncase+ncontrol) as passed to get_r_from_lor for the de-ascertainment identity, and
# wrong for phi, where it is meant to describe the study.
#
# Used as: a variant analysed in N subjects contributes n_eff = 4*phi*(1-phi)*N. That is an
# explicitly approximate missingness model -- in a meta-analysis, which cohorts contribute to a
# given variant can change the case:control mix -- so the resulting z is approximate. The Steiger
# DIRECTION does not depend on it.
# DERIVED from the analysed GWAS, not transcribed. The five FinnGen entries used to be Risteys R10
# registry prevalences under a comment claiming they "agree with Risteys R10 to the third decimal".
# For two of them that was false: T2D is 65,085/400,197 = 0.1626 against the 0.1580 recorded, and AF
# is 50,743/261,395 = 0.1941 against 0.1232 -- AF's registry figure divides by the whole cohort
# (412,181) while the GWAS analysed 261,395. A registry prevalence answers "how common is this in
# Finland"; phi has to answer "what was the case:control mix in the analysis these betas came from",
# and those differ whenever the endpoint's control set is not the whole cohort.
#
# OUTCOME_N_CASE / OUTCOME_N_TOTAL are the right source because config checks the five FinnGen values
# against finngen_R10_manifest.tsv at load time, so a transcription error cannot survive there either.
OUTCOME_STUDY_CASE_FRACTION <- local({
  if (!exists("OUTCOMES") || !exists("OUTCOME_N_CASE") || !exists("OUTCOME_N_TOTAL"))
    stop("R/_outcome_prevalence.R: source R/00_config.R first. The study case fractions are DERIVED ",
         "from OUTCOME_N_CASE / OUTCOME_N_TOTAL; without them this file would have to hardcode them ",
         "again, which is the defect it was rewritten to remove.")
  .out <- list()
  for (.k in names(OUTCOMES)) {
    .nc <- OUTCOME_N_CASE[[.k]]; .nt <- OUTCOME_N_TOTAL[[.k]]
    if (is.null(.nc) || is.na(.nc)) next          # quantitative outcome: no case fraction to define
    if (is.null(.nt) || is.na(.nt) || .nt <= 0 || .nc <= 0 || .nc >= .nt)
      stop(sprintf("R/_outcome_prevalence.R: %s has ncase=%s of N=%s, which cannot be a case fraction",
                   .k, .nc, .nt))
    .out[[unname(OUTCOMES[[.k]])]] <- .nc / .nt
  }
  .out
})
