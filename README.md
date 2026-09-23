# Paramo_SkyIslands_RootFungi

This repository contains the data, scripts, and analysis documents supporting the manuscript:

> **Parallel elevation filtering of Ericaceae root-associated fungal communities in Andean páramo ecosystems** — *[Authors]* — *[Journal, Year]*

Sequence data are deposited at ENA under accession **PRJEB107725**.

## Supplementary Appendix

- Web version (with all code in drop-downs): https://ssangulo.github.io/Paramo_SkyIslands_RootFungi/appendix.html
- Repository copy of the rendered HTML: `Appendix_HTML.html` (identical to `docs/appendix.html`, which GitHub Pages serves)
- Code-free PDF version for the journal: `Appendix_PDF.pdf`
- Quarto source: `Analysis_pipeline.qmd`

The appendix chunks are `eval: false`: rendering reads the saved tables in `tables/` and figures in `figures/`, and does not re-run analyses. Render from the repository root:

```
quarto render Analysis_pipeline.qmd --to html    # -> Appendix_HTML.html (copy to docs/appendix.html)
quarto render Analysis_pipeline.qmd --to pdf     # -> Appendix_PDF.pdf (needs xelatex)
```

## Repository contents

| Folder / file | Contents |
| ---- | ----------- |
| `Analysis_pipeline.qmd` | Quarto source of the Supplementary Appendix |
| `Data_files/` | Data S1 (sample metadata), Data S2 (genus-level GLMM results), Data S3 (soil physicochemistry) |
| `Scripts/` | Bioinformatic pipeline and all analysis scripts (below) |
| `objects/` | Saved R objects: `ps_individual.rds` (per-plant phyloseq object, 60 plants), `tree_ps.rds` (phyloseq object with the ML tree, 120 root samples), and intermediate outputs of the scripts |
| `tables/` | All result tables (CSV) read by the appendix |
| `figures/` | Main-text figures (PNG/JPEG/EPS/TIFF) and supplementary figures (PNG) |
| `logs/` | Console output of the last full run of `Scripts/run_revision.sh` |
| `docs/` | GitHub Pages site (rendered appendix) |
| `archive/` | Superseded first-round analyses (GLLVM, original Figure 4 alternative, pre-revision tables), kept for comparison; see `archive/README.md` |

## Scripts

### Bioinformatics and figures (not part of the scripted run)

| Script | What it does |
|---|---|
| `full_pipeline_DADA2.R` | Read processing (demultiplexing, trimming, DADA2 denoising, chimera removal, LULU), control and soil subtraction, 2-of-4 PCR replicate filter, taxonomy (UNITE), per-plant collapsing; writes `ps_individual.rds`. Uses server paths; kept as a record of the pipeline |
| `funguild_script.R` | FUNGuild guild assignment and Figure S6 |
| `map_code.R` (+ `map_code_panelB.py`, called by it) | Figure 1 (study-area map and elevation profiles). Needs an R with the spatial stack (see the script header); downloads a DEM to `objects/` on first run |

### Analyses (run in this order by `bash Scripts/run_revision.sh`, from the repository root)

| Step | Script | What it does | Appendix |
|---|---|---|---|
| 1 | `iNEXT3D_plant_level.R` | Habitat-level Hill numbers (TD, meanPD; q = 0–2) with plants as incidence units, coverage-standardised estimates and pairwise habitat comparisons; Figure 3A. Slowest step (nboot = 500) | 5.2 |
| 2 | `perplant_hill_size_standardised.R` | Per-plant Hill numbers rarefied to 28,122 reads; gamma GLMs (habitat + site; interaction by model comparison; Tukey contrasts); Figure 3B and the assembled Figure 3 | 5.3 |
| 3 | `parallelism_diagnostics.R` | Which site drives the PD interaction; parallelism tests and site slopes for Helotiales and Sebacinales | 5.3, 6.1 |
| 4 | `Figure4_pooled_glmm.R` | Per-genus negative binomial GLMMs (Figure 4C, Data S2) and main-text Figure 4 (panels A–D) | 6.2–6.3 |
| 5 | `helotiales_fraction_check.R` | Helotiales elevation trend split by genus assignment and family | 6.3 |
| 6 | `adjusted_tests_centroid_nri.R` | Site-adjusted tests for the centroid model and NRI/NTI, Tukey contrasts, location-mean sensitivity. Reuses the cached null-model runs in `objects/nri_nti_per_plant_999.rds` | 4.3, 8 |
| 7 | `reproduce_composition_tables.R` | Verifies that Tables S5–S10 reproduce under vegan ≥ 2.7-1 with the explicit zeros-retained robust CLR | 1, 4 |
| 8 | `model_selection_aicc.R` | AIC, AICc, BIC and Akaike weights for every additive vs interaction decision; Jacobian-corrected family screen | 10 |
| 9 | `family_sensitivity.R` | Per-plant habitat results and model structure under gaussian, log-normal and gamma | 5.3, 10 |

Steps 2–3 depend on step 1 (for the combined Figure 3) and step 2; step 8 needs steps 2, 4 and 6; step 9 needs step 2. The per-plant models use a gamma GLM with a log link (`FAMILY` at the top of steps 2 and 3).

Software: revision analyses were run in R 4.5.2 (vegan 2.7.2, phyloseq 1.54.0, iNEXT.3D 1.0.12, glmmTMB 1.1.14) and reproduced in R 4.3.3 (vegan 2.7.2, phyloseq 1.46.0, iNEXT.3D 1.0.12, glmmTMB 1.1.9, betareg 3.2.4, emmeans 2.0.4, picante 1.8.2) with identical tables.

## File names and appendix numbers

Some file names predate the final numbering of the appendix. The appendix number of each file is:

**Figures**

| Appendix | File |
|---|---|
| Figure 1 (main text) | `figures/Fig_1_map_regenerated.*` |
| Figure 2 (main text) | `figures/Figure_2_dbRDA_centroid.*` |
| Figure 3 (main text) | `figures/Figure_3_corrected.*` (panels: `Figure_3A_iNEXT_plant_units.png`, `Figure_3B_perplant_size_std.png`) |
| Figure 4 (main text) | `figures/Figure_4_option2_pooled.*` (panel C alone: `genus_glmm_heatmap.png`) |
| Figure S1 | `figures/Fig_S1_Venn.png` |
| Figure S2 | `figures/Fig_S2.png` |
| Figure S3 | `figures/Fig_S3.png` |
| Figure S3b | `figures/diag_parallelism_sites.png` |
| Figure S4 | `figures/Fig_S5_betareg_diagnostics_helotiales.png` |
| Figure S5 | `figures/Fig_S6_betareg_diagnostics_sebacinales.png` |
| Figure S6 | `figures/Fig_S7_FunGuild.png` |
| Figure S7 | `figures/Fig_S10_soil.png` |

**Tables** (Tables S1, S2 and S4 are written directly in the appendix)

| Appendix | File (`tables/` unless stated) |
|---|---|
| Table S3 | `Table_S3_OTU_filtering_summary.csv` |
| Tables S5–S6 | `Table_S5_dbRDA_terms.csv`, `Table_S6_dbRDA_axes.csv` |
| Table S7 | `Table_S7_PERMANOVA.csv` |
| Table S8 | `Table_S8_dispersion_habitat.csv`, `Table_S8_betadisper_test_habitat.csv` |
| Table S9 | `Table_S9_betadisper_test_site.csv` (group distances: `Table_S9_dispersion_site.csv`) |
| Table S10 | `Table_S10_centroid_model.csv` |
| Tables S10b–S10c | `adjusted_tests_centroid_nri.csv`, `adjusted_pairwise_centroid_nri.csv` (`dist_centroid` rows) |
| Table S11 | `Table_S11_Sorensen_partitioning.csv` (abundance-based Bray–Curtis partitioning, not shown: `Table_S11_Bray_partitioning.csv`) |
| Table S12 | `Table_S16_iNEXT3D_TD_plant.csv` |
| Table S13 | `Table_S17_iNEXT3D_PD_plant.csv` |
| Table S13b | `iNEXT3D_estimates_plant_vs_root.csv` (plant units, coverage-standardised rows) |
| Table S13c | `iNEXT3D_pairwise_habitat_plant_vs_root.csv` (plant units, coverage-standardised rows) |
| Tables S14–S15 | `perplant_hill_size_std_omnibus.csv` (TD / PD rows) |
| Table S15b | `perplant_hill_size_std_pairwise.csv` |
| Table S15c | `diag_leave_one_site_out.csv` |
| Table S15d | `diag_within_site_contrasts.csv` (PD, q = 1–2) |
| Table S15e | `diag_PD_interaction_vs_lineages.csv` |
| Table S16 | `Table_S12_betareg_helotiales.csv` |
| Table S17 | `Table_S13_betareg_sebacinales.csv` |
| Table S18 | `Table_S14_betareg_other_orders.csv` |
| Table S18b | `diag_lineage_parallelism_tests.csv` |
| Table S18c | `diag_lineage_site_slopes.csv` |
| Table S18d | `helotiales_fraction_check.csv` |
| Data S2 | `Data_files/Data_S2_genus_GLMM_table.csv` (= `genus_glmm_results.csv`; presence per site: `genus_glmm_presence_by_site.csv`) |
| Table S19 | `nri_nti_median_iqr_per_plant.csv` |
| Table S20 | `adjusted_tests_centroid_nri.csv` (`NRI_aw` / `NTI_aw` rows) |
| Table S21 | `family_sensitivity_structure.csv` |
| Table S22 | `model_selection_aicc.csv` (rows other than per-plant diversity) |
| Table S23 | `model_selection_aicc_genus.csv` |

Other CSVs in `tables/` are supporting outputs of the scripts (per-plant values, emmeans, trends, sensitivity refits, reproduction checks).
