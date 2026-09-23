# Archive: superseded first-round analyses

Files in this folder belong to the first-round (pre-revision) version of the analysis. None of them is read by
`Analysis_pipeline.qmd` or by any script in `Scripts/`. They are kept so that the numbers reviewers saw in the
first round can still be checked against the revised analysis. Why each was replaced, and what replaced it, is
given below; the replacements are described in the Supplementary Appendix (sections given in brackets).

## Scripts, figures and objects

| File | Why it was superseded | What replaced it |
|---|---|---|
| `Scripts/gllvm_analysis.R` | GLLVM of OTU-level habitat responses (server script). Genus-level summaries were maxima over per-OTU coefficients, many of which were not identifiable because the OTU never occurred in the forest reference habitat (\|β\| up to ~300 on the log scale). The random effects were community-level row effects only, and the published model comparison cannot be reproduced from this script | Per-genus negative binomial GLMMs in `Scripts/Figure4_pooled_glmm.R` (Appendix 6.3, Figure 4C, Data S2) |
| `objects/p_heat.rds` | GLLVM heatmap (published Figure 4C) | `objects/p_heat_glmm.rds` |
| `Scripts/Figure4_option1_site_slopes.R`, `figures/Figure_4_option1_site_slopes.png` | Alternative Figure 4A–B with site-specific elevation slopes. Site-specific slopes are not supported (Appendix Table S18b), so the across-site version was adopted | `Scripts/Figure4_pooled_glmm.R` → `figures/Figure_4_option2_pooled.png` (main-text Figure 4) |

The archived Figure 4 script still runs from the repository root (`Rscript archive/Scripts/Figure4_option1_site_slopes.R`,
after `Scripts/Figure4_pooled_glmm.R`). `gllvm_analysis.R` uses server paths and objects not in this repository and
is kept as a record only.

## Tables

| File | Why it was superseded | What replaced it |
|---|---|---|
| File | Why it was superseded | What replaced it |
|---|---|---|
| `tables/Table_S16_iNEXT3D_TD.csv` | Habitat-level iNEXT3D run on 120 root samples, so the two roots of a plant counted as independent units and the CIs were too narrow | `tables/Table_S16_iNEXT3D_TD_plant.csv` (appendix Table S12) |
| `tables/Table_S17_iNEXT3D_PD.csv` | Same problem, phylogenetic diversity | `tables/Table_S17_iNEXT3D_PD_plant.csv` (appendix Table S13) |
| `tables/Table_S18_TD_parallelism_anova.csv` | Sequential (Type I) ANOVA at q = 1 from the interaction model, so habitat was not adjusted for site | `perplant_hill_size_std_omnibus.csv`, TD rows (appendix Table S14) |
| `tables/Table_S19_PD_parallelism_anova.csv` | Same problem, phylogenetic diversity | `perplant_hill_size_std_omnibus.csv`, PD rows (appendix Table S15) |
| `tables/Table_S15b_habitat_posthoc.csv` | Tukey contrasts at q = 1 on unrarefied hillR values, taken from a different model than the omnibus test above | `perplant_hill_size_std_pairwise.csv` (appendix Table S15b) |
| `tables/File_S2_GLLVM_table.csv` | GLLVM genus summaries whose coefficients had diverged under separation (\|beta\| up to 260 on the log scale) | `genus_glmm_results.csv` + `genus_glmm_presence_by_site.csv` (appendix Data S2) |
| `tables/Data_S2_GLLVM_table.csv` | The `Data_files/` copy of the GLLVM genus table, byte-identical to `File_S2_GLLVM_table.csv` above | `Data_files/Data_S2_genus_GLMM_table.csv` (a copy of `tables/genus_glmm_results.csv`) |
| `tables/Table_S20_NRI_NTI_median_IQR.csv` | Per-plant medians recomputed with the corrected inputs | `nri_nti_median_iqr_per_plant.csv` (appendix Table S19) |
| `tables/Table_S21_NRI_NTI_ANOVA.csv` | Type I tests, so habitat was not adjusted for site | `adjusted_tests_centroid_nri.csv`, NRI_aw / NTI_aw rows (appendix Table S20) |

Tables S5–S10 (composition) are **not** superseded: they reproduce exactly under the explicit zeros-retained robust
CLR and remain in `tables/` (see Appendix Section 1, "Software versions", for the vegan ≥ 2.7-1 change).
