# Superseded tables

These CSVs were produced before the revision re-analyses (see `Scripts/REVISION_NOTES.md`
and `AUDIT_NOTES.md`). They are kept only as a record of the published version and are no
longer read by `Analysis_pipeline.qmd` or by any script in `Scripts/`.

| File | Why it was superseded | What replaced it |
|---|---|---|
| `Table_S16_iNEXT3D_TD.csv` | Habitat-level iNEXT3D run on 120 root samples, so the two roots of a plant counted as independent units and the CIs were too narrow | `Table_S16_iNEXT3D_TD_plant.csv` (appendix Table S12) |
| `Table_S17_iNEXT3D_PD.csv` | Same problem, phylogenetic diversity | `Table_S17_iNEXT3D_PD_plant.csv` (appendix Table S13) |
| `Table_S18_TD_parallelism_anova.csv` | Sequential (Type I) ANOVA at q = 1 from the interaction model, so habitat was not adjusted for site | `perplant_hill_size_std_omnibus.csv`, TD rows (appendix Table S14) |
| `Table_S19_PD_parallelism_anova.csv` | Same problem, phylogenetic diversity | `perplant_hill_size_std_omnibus.csv`, PD rows (appendix Table S15) |
| `Table_S15b_habitat_posthoc.csv` | Tukey contrasts at q = 1 on unrarefied hillR values, taken from a different model than the omnibus test above | `perplant_hill_size_std_pairwise.csv` (appendix Table S15b) |
| `File_S2_GLLVM_table.csv` | GLLVM genus summaries whose coefficients had diverged under separation (\|beta\| up to 260 on the log scale) | `genus_glmm_results.csv` + `genus_glmm_presence_by_site.csv` (appendix Data S2) |
| `Data_S2_GLLVM_table.csv` | The `Data_files/` copy of the GLLVM genus table, byte-identical to `File_S2_GLLVM_table.csv` above | `Data_files/Data_S2_genus_GLMM_table.csv` (a copy of `tables/genus_glmm_results.csv`) |
| `Table_S20_NRI_NTI_median_IQR.csv` | Per-plant medians recomputed with the corrected inputs | `nri_nti_median_iqr_per_plant.csv` (appendix Table S19) |
| `Table_S21_NRI_NTI_ANOVA.csv` | Type I tests, so habitat was not adjusted for site | `adjusted_tests_centroid_nri.csv`, NRI_aw / NTI_aw rows (appendix Table S20) |

Tables S5–S10 (composition) are **not** superseded: they reproduce exactly and remain in
`tables/`. See issue 8 of the revision notes for the vegan rCLR change behind that check.
