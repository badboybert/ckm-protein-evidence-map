# ckm_palette.R -- the locked Okabe-Ito-based CKM house palette (from the cross-paper style survey,
# `paper 6/analysis/src/fig_common.R` + `_skills_extracted/publication-figures-tables/color_palettes.R`).
# Source this alongside theme_publication_ckm.R in every v3.2 figure generator.

# colocalization / evidence-tier key (reused across ALL figures)
CKM_COLOC <- c(causal      = "#0072B2",   # shared-variant / coloc-causal (H4)
               distinct    = "#D55E00",   # H3-favored / distinct-variant / LD-confounded
               ambiguous   = "#9AA0A6",   # ambiguous / tested-not-supported
               prior_stable= "#04386b")   # prior-robust (dark navy)

# direction of effect
CKM_DIR <- c(risk = "#B2182B", protective = "#2166AC", null = "grey60")

# CKM outcomes (organ-system keyed)
CKM_OUTCOME <- c(T2D="#E63946", CAD="#A8201A", HF="#9C27B0", Stroke="#1E88E5",
                 AF="#00ACC1", CKD="#1B7837", eGFR="#66BB6A", NAFLD="#E67E22")

# cross-ancestry: EUR = grey circle, EAS = blue triangle (redundant shape encoding)
CKM_ANC_COL   <- c(EUR = "grey55", EAS = "#0072B2")
CKM_ANC_SHAPE <- c(EUR = 16,       EAS = 17)

# significance shape (redundant encoding; filled = supported, open = not)
CKM_SIG_SHAPE <- c(`TRUE` = 16, `FALSE` = 1)

# deflation gradient (tested -> colocalized -> anchor), Okabe-Ito hues
CKM_DEFLATE <- c(tested = "#9AA0A6", colocalized = "#0072B2", anchor = "#009E73")
