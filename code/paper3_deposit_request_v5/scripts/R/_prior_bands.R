# ============================================================================
# PRIOR-SENSITIVITY BANDS: ONE DEFINITION (editor ruling 2026-08-01).
#
# `stable` used to mean three different things: phase3/10 and phase3/25 required PP.H4 >= 0.80 at
# BOTH p12 = 1e-6 and 1e-5; reporting/prior_sweep.R used 1e-6 only; the Methods prose said three
# priors. They agreed only by construction -- p12 = 1e-5 IS coloc.abf's default, so that arm re-runs
# the main analysis and is >= 0.80 for every coloc-supported pair (verified on this generation:
# 76/76 at 1e-5, and 23 either way at 1e-6). A definition that cannot disagree cannot be audited,
# and Figure 4's caption presented a 3-band partition built from two of those rules.
#
# WHY THIS IS ITS OWN FILE AND NOT R/00_config.R. That is where is_coloc_pass() lives and where
# these belonged by kinship -- but 00_config.R is hashed into the COLOC INPUTS KEY
# (R/_coloc_inputs_key.R: `config`), which every coloc artifact is stamped with and which
# 47_a6_stage_fig_inputs.R verifies before staging. Adding a classifier there changed the key from
# 84be48207e18 to 0d4d27c1b88c and made this generation's staged artifacts unverifiable -- a
# post-hoc CLASSIFIER of sweep output cannot be allowed to invalidate the colocalisation it
# describes. Nothing in here participates in any computation the key exists to fingerprint.
#
# The band is decided at the MOST CONSERVATIVE prior, the only reading under which the three bands
# partition the coloc-supported set. The 0.50 line has no other role anywhere in the study.
# ============================================================================
PRIOR_STABLE_MIN  <- 0.80   # PP.H4 at p12 = 1e-6 at or above this -> robust
PRIOR_FRAGILE_MAX <- 0.50   # ... below this -> fragile; in between -> attenuated

is_prior_stable <- function(pp_h4_1e6) {
  !is.na(pp_h4_1e6) & is.finite(pp_h4_1e6) & pp_h4_1e6 >= PRIOR_STABLE_MIN
}

# Band names per the same ruling: "Prior-sensitive" collided with a differently-defined `sensitive`
# flag elsewhere (42 vs 23 on the previous generation), so the middle band is "Prior-attenuated".
prior_band <- function(pp_h4_1e6) {
  data.table::fcase(
    is.na(pp_h4_1e6) | !is.finite(pp_h4_1e6), NA_character_,
    pp_h4_1e6 >= PRIOR_STABLE_MIN,            "Prior-robust",
    pp_h4_1e6 >= PRIOR_FRAGILE_MAX,           "Prior-attenuated",
    default =                                 "Prior-fragile")
}
