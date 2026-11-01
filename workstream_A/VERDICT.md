# VERDICT: Workstream A (Chen 2023 instrument rs2246012 vs ARG1 protein signals A and B)

Pre-declaration: `P31R_R2_STAGE2_AMEND2_v3.md`, md5 7daeec61dde5bd9a98956498c44a9f34 (asserted by `scripts/b01_ld.R` before any computation).
Status: **builder result, not yet independently checked.** The Opus checker (Python, raw genotypes) has not been run; see "Not done".

## Category: STRONGLY TAGS A

| Panel | n | r2 rs2246012 vs A (rs2781668) | r2 rs2246012 vs B (rs17788484) | Category |
|---|---|---|---|---|
| 503 (primary) | 503 | 0.99274860605258064 | 0.0020139120068619937 | STRONGLY TAGS A |
| 489 (nested subset) | 489 | 0.99257549411518131 | 0.0013035738289601778 | STRONGLY TAGS A |

Final category (503 primary, 489 agrees): **STRONGLY TAGS A**. Rule applied to unrounded values: r2 with A at least 0.80 and r2 with B below 0.20, in both panels.
Context pair, not part of the verdict: r2 A vs B = 0.0019781 (503), 0.0012679 (489); signed r = -0.0445 and -0.0356, which equal the r(A,B) already recorded in the Stage 1 pre-declaration (-0.0445 on 503, -0.036 on 489).

Wording the pre-declaration permits: rs2246012 strongly tags the lead of signal A (rs2781668) and has no pairwise LD with the lead of signal B, in a 503-person and a nested 489-person European reference panel. This is a description of pairwise tagging. It is not proof of an identical causal signal, of an independent causal perturbation, or of protein-mediated disease causation, and it changes no existing branch.

## Reading it correctly

1. **This is a confirmation, not a discovery.** `S1R1_ARG1_PRESTATED_RULE.md` line 9 (2026-10-02; disclosed in AMEND2 v3 section 4) already recorded r2 0.993 between rs2246012 and rs2781668. The genotype calculation reproduces it (0.99275 on 503) and adds the 489 panel and the comparison with B. I did not open that file myself. I only know its content from AMEND2 v3 section 4.
2. **The two panels are not two cohorts.** The 489 panel is the 1000G Phase 3 EUR QC set, whose 489 IIDs are all in the 503 panel (489/489). As a diagnostic, the 503 panel restricted to those 489 IIDs gives r2 identical to the 489 panel to full precision (rows `r2.*.503in489ids`). The agreement between panels is a check on sample composition and file handling, not independent replication.
3. **Orientation.** The three variants are non-palindromic. The signed r between rs2246012-C and rs2781668-T is +0.996, so the two alleles whose effects are reported (C and T) sit on the same haplotype. Effects at the two variants are numerically close in every source that reports both (UKB-PPP ARG1 protein -0.236 vs -0.237; CLSA arginine +0.173 vs +0.172; FinnGen T2D +0.0454 vs +0.0453). With r2 of 0.99 this is expected and is not additional evidence of identity.
4. **Credible sets (existing outputs, read only).** rs2246012 and rs2781668 are both in the protein sets of signal A under the in-sample UKB-PPP fit (admissible) and the 503-reference SuSiE fit (inadmissible under Rule 3(i), s max 0.1482). rs2246012 is also in the FinnGen T2D SuSiE set (inadmissible) and the T2D ABF set (method 4) at A. Methods 3 and 5 (503) list four and five variants that do not include either variant; their 489 versions list exactly rs2781668 and rs2246012 (flagged "NOT robust to the reference LD" in W4b). rs17788484 (B) is in the B sets of methods 1 (protein), 2, 3 and 5 and in none of the A sets. Full 54-row table: `credible_set_status.tsv`. Membership in an inadmissible fit is an existing observation, not support for a main claim. Overlap of credible sets does not establish a common causal variant.
5. **Build.** AMEND2 section 4 noted that the build of Chen's chr6:131,577,068 "was not established". It is established: Chen ST5 labels that column `Position(hg38)`, and Chen ST17 gives POS_hg38 131,577,068 with POS_hg37 131,898,208. The latter equals the pos37 in the harmonised table and the panel .bim. No lift-over was performed and no chain file was used.

## Not done, by design or by sequence

- The independent checker (Opus, Python from raw genotype files, results saved before seeing these) is not run. This verdict stays "builder result" until it is, and until the numerical-agreement convention is recorded and applied.
- No new fine-mapping, association estimate, MR or download. Chen's OR 1.238 and PP.H4 0.877 are quoted as published (ST19 row 130, ST24 row 25), not re-derived.
- Not written: the Additional file 1 sheet and the decisions-log entry (section 6 puts them downstream of the check).

## Limits

- 489 and 503 share the same people for the 489 (point 2); n for rs17788484 rests on 22 and 20 copies of the minor allele (503 and 489). The B comparison is far from both cut-offs (r2 about 0.002), so the category does not depend on that sparse count; the precision of r2 with B is low.
- Method 2 (Sun 2023 ST16 in-sample set): the file lists the set only, so for a variant not in the set I cannot say whether it was in the fit's input; recorded as such, not forced into "tested" or "not in input".
- Methods 3 to 5 are ABF/colocalization on z-space-conditioned data, not SuSiE. The Rule 3(i) wording ("SuSiE with reference LD") does not cover them, so their admissibility is recorded as "unresolved under Rule 3(i) wording", with the Rule 3(ii) robustness label alongside. The 489 versions of these sets are parsed from the `extra` field of `credible_sets_by_method.tsv`; W4b stores no 489 version for methods 1 and 2.
- Rule 3(i) scope was checked against `P31R_R2_STAGE1_PRESPEC.md` (md5 3319dc39...) line 85: it governs SuSiE evidence only. The descriptive genotype r2 here is not SuSiE evidence, so no conflict with AMEND2 section 3 arose and no exception was invented.
- Frequencies differ by source (rs2246012 effect allele C: 0.160 Chen/CLSA, 0.167 UKB-PPP, 0.171 and 0.171 panels, 0.231 FinnGen). This does not affect r2, which is computed within each reference panel.
- Prior exposure of the builder: AMEND2 v3 section 4 (credible-set memberships, beta 0.1733 vs 0.172, the 0.10 bound, the 0.993 record) and the W4b credible-set outputs were read before computing. The LD values and the category were first seen at 2026-10-13 13:25. This analysis is not blind.
