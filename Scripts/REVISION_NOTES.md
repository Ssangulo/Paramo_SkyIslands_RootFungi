# Revision notes: re-analyses and what they change

These notes cover the re-analyses that came out of answering Reviewer 2's comments on collapsing and post-hoc habitat tests. Once we looked closely, a few problems turned up beyond what the reviewer asked about, some small and some that change what we can claim. Everything here was run on the current `objects/ps_individual.rds` and `objects/tree_ps.rds` (R 4.5.2, vegan 2.7.2, phyloseq 1.54.0, iNEXT.3D 1.0.12, glmmTMB 1.1.14).

Issues are ordered roughly by how much they change the paper. For each one there is a short description of the problem, what we did, what it showed, and a table of what the manuscript (or appendix) says now against what the new results show. Manuscript locations refer to the track-changes draft *E&E revisions – Root-fungi paper draft – V2*.

## What checks out

- The pipeline does what the Methods say: soil subtraction on raw counts, then the 2-of-4 PCR replicate filter (within each root sample), then PCR replicates summed, then the two roots summed per plant.
- The beta regressions (Tables S16–S18) and the original per-plant Hill ANOVA reproduce exactly from `ps_individual`.
- The composition tables (S5–S10: dbRDA, PERMANOVA, betadisper, centroid model) reproduce exactly, as long as the robust CLR is computed without zero imputation (see issue 8).
- NRI/NTI (Tables S19–S20) were computed on the 60 plants, as they should be.
- betapart (Table S11) was already redone on the 60 plants.

## Scripts

| Script | What it does |
|---|---|
| `iNEXT3D_plant_level.R` | Habitat-level Hill numbers with plants (not root samples) as incidence units, pairwise habitat comparisons, corrected Figure 3A |
| `perplant_hill_size_standardised.R` | Per-plant Hill numbers rarefied to 28,122 reads, omnibus and post-hoc tests from one model, corrected Figure 3B |
| `parallelism_diagnostics.R` | Which site drives the PD interaction; tests whether Helotiales and Sebacinales trends are parallel across sites |
| `Figure4_option1_site_slopes.R`, `Figure4_option2_pooled.R` | Figure 4A–B redrawn without forcing parallel site lines (two options) |
| `genus_glmm_habitat.R` | Genus-level habitat responses replacing the GLLVM summaries (Figure 4C, Data S2) |
| `helotiales_fraction_check.R` | Splits the Helotiales elevation trend into genus-assigned and unassigned OTUs, and by family |
| `adjusted_tests_centroid_nri.R` | Site-adjusted tests for the centroid model and NRI/NTI, plus location-mean sensitivity |
| `reproduce_composition_tables.R` | Reproduces Tables S5–S10 with the robust CLR used for the published results |
| `model_selection_aicc.R` | AIC, AICc, BIC and Akaike weights for every additive vs interaction decision, plus the family screening redone |
| `family_sensitivity.R` | Per-plant habitat results refitted under gaussian, lognormal and gamma |

### Run order

Everything runs from the repository root. `bash Scripts/run_revision.sh` runs all of it in this order and saves each script's console output to `logs/`. The outputs from our runs are already committed in `tables/`, `figures/` and `objects/`, so you only need to re-run to reproduce them.

| Step | What to run and what it needs |
|---|---|
| 1 | `Rscript Scripts/iNEXT3D_plant_level.R`. Needs `ps_individual.rds` and `tree_ps.rds`. Slowest step (iNEXT3D bootstraps, nboot = 500; add `20` as an argument for a quick test). Writes the plant-unit tables, `out_TD_plant.rds`, `out_PD_plant.rds` and the Figure 3A panel |
| 2 | `Rscript Scripts/perplant_hill_size_standardised.R`. Needs step 1 for the combined Figure 3 (it still runs without it, but only saves Figure 3B). Writes `perplant_hill_size_std.rds` and the per-plant tables |
| 3 | `Rscript Scripts/parallelism_diagnostics.R`. Needs `perplant_hill_size_std.rds` from step 2 |
| 4 | `Rscript Scripts/genus_glmm_habitat.R`. Writes `genus_glmm_results.csv` and the heatmap `p_heat_glmm.rds` |
| 5 | `Rscript Scripts/helotiales_fraction_check.R` |
| 6 | `Rscript Scripts/Figure4_option2_pooled.R` (main-text Figure 4) and, for reference only, `Figure4_option1_site_slopes.R`. Both need `p_heat_glmm.rds` from step 4 |
| 7 | `Rscript Scripts/adjusted_tests_centroid_nri.R`. Slow the first time (999 null-model runs); reuses `objects/nri_nti_per_plant_999.rds` if it exists |
| 8 | `Rscript Scripts/reproduce_composition_tables.R` |
| 9 | `Rscript Scripts/model_selection_aicc.R`. Needs steps 2, 4 and 7 |
| 10 | `Rscript Scripts/family_sensitivity.R`. Needs step 2 |

Steps 2 and 3 use a gamma GLM with a log link (`FAMILY` at the top of each; set both to `"lognormal"` to reproduce the earlier version). Step 9 recomputes the structure comparison on the log-normal fits, and step 10 repeats it within each family.

Steps 4, 5, 7 and 8 don't depend on steps 1–3 and can be run in any order. `full_pipeline_DADA2.R`, `funguild_script.R` and `gllvm_analysis.R` are the original server scripts and aren't part of the revision run (the GLLVM is superseded by step 4).

---

## How model structure was decided

Every additive vs interaction decision below was made with a likelihood-ratio test or a model-comparison F-test, and then checked against AIC, AICc, BIC and Akaike weights (`model_selection_aicc.R`). For context, the published analysis used AIC only to pick the error family for the per-plant Hill models (appendix 5.3) and the distribution and latent-variable structure for the GLLVM. AICc was not used anywhere, and interaction structure was decided on p-values.

The rule applied throughout: where the interaction is not supported, the additive model is the model, and the additive model is what gets reported. Site-specific estimates from an unsupported interaction are not reported as results.

The two frameworks agree everywhere except two places, both of which make the revision more conservative rather than less:

| Comparison | Test | ΔAICc (interaction − additive) | Agreement |
|---|---|---|---|
| Per-plant TD, all q | p = 0.21–0.25 | +4.5 to +5.0 | additive, both |
| Per-plant PD q = 0 | p = 0.38 | +6.4 | additive, both |
| Per-plant PD q = 1 | p = 0.008 | −4.8 (weight 0.89) | interaction, both |
| Per-plant PD q = 2 | p < 0.001 | −15.7 (weight 0.9995) | interaction, both |
| Centroid distance | p = 0.29 | +5.5 | additive, both |
| NTI_aw | p = 0.11 | +2.5 | additive, both |
| NRI_aw | p = 0.053 | +0.35 | **disagree**: AICc prefers site only over both; additive retained |
| Helotiales, habitat × site | p = 0.042 | +1.5 (weight 0.28) | **disagree**: interaction not supported; additive retained |
| Helotiales, elevation × site | p = 0.060 | +0.6 (weight 0.38) | interaction not supported; additive retained |
| Sebacinales, both | p = 0.19, 0.33 | +4.6, +5.1 | additive, both |

On that basis the interaction is retained only for abundance-weighted PD at q = 1 and q = 2 (issue 3). Everywhere else the additive model stands, including both dominant lineages (issue 4).

AICc also has something to say about whether habitat belongs in the model at all, which the p-values were never asked. Habitat is supported for per-plant TD q = 0 (ΔAICc −9.9 vs site only), PD q = 0 (−10.7), centroid distance (−7.2) and both PD q = 1–2 models; it is not supported for TD q = 1 (+0.9), TD q = 2 (+2.8), NRI (+3.4) or NTI (+3.2). That matches the significance results exactly.

---

## 1. Per-plant diversity tests (the reviewer's post-hoc comment)

**The problem.** Two things. First, the habitat test in Tables S14–S15 came from `anova()` on `log(diversity) ~ habitat * site`, which gives sequential (Type I) tests, so habitat was tested without adjusting for site. The Tukey contrasts in Table S15b came from the additive model, where habitat is adjusted for site. Because BEL has no páramo and MA no forest, habitat and site are partly confounded, so the two tests were asking different questions. That is why we had a significant habitat effect but no pair of habitats that differed. Second, per-plant Hill numbers were calculated with `hillR` on unrarefied counts (libraries range from 28,122 to 882,416 reads), and `hillR` PD isn't on a common reference time, although the manuscript says PD was.

**What we did.** Rarefied every plant to 28,122 reads with `iNEXT.3D::estimate3D` (TD and meanPD, q = 0, 1, 2), then took all tests from one model: habitat adjusted for site, the habitat × site interaction by model comparison, and Tukey contrasts from the same additive model. Coverage-based standardisation isn't an option here: no plant has a single singleton OTU (a side effect of the 2-of-4 replicate filter), so coverage is ≥ 0.999 at any depth.

**What it showed.**

- The published numbers reproduce exactly (TD F₂,₅₀ = 4.39, p = 0.017). Adjusted for site, the same values give F₂,₅₄ = 1.92, p = 0.16 (PD: 3.03 becomes 0.72, p = 0.49). The published habitat effect was site confounding.
- On the rarefied values (gamma GLM with a log link; see the family results below):

| Metric | Habitat (adjusted for site) | Habitat × site | Forest vs subpáramo | Forest vs páramo | Subpáramo vs páramo |
|---|---|---|---|---|---|
| TD q = 0 | F₂,₅₄ = 9.21, p = 0.0004 | p = 0.22 | 1.52× (1.18–1.96), p = 0.0006 | 1.69× (1.26–2.26), p = 0.0002 | 1.11×, p = 0.64 |
| TD q = 1 | F₂,₅₄ = 3.21, p = 0.048 | p = 0.24 | 1.48× (0.99–2.20), p = 0.055 | 1.53×, p = 0.078 | 1.03×, p = 0.98 |
| TD q = 2 | p = 0.15 | p = 0.25 | 1.42×, p = 0.14 | 1.33×, p = 0.38 | 0.93×, p = 0.93 |
| PD q = 0 | F₂,₅₄ = 9.61, p = 0.0003 | p = 0.42 | 1.39× (1.13–1.72), p = 0.001 | 1.58× (1.24–2.02), p = 0.0001 | 1.13×, p = 0.40 |
| PD q = 1 | p = 0.023 | **p = 0.009** | 1.18×, p = 0.12 | 1.27×, p = 0.047 | 1.07×, p = 0.73 |
| PD q = 2 | p = 0.015 | **p = 0.0001** | 1.16×, p = 0.023 | 1.15×, p = 0.094 | 0.99×, p = 0.97 |

Linear polynomial contrasts across the ordered habitats are significant for TD q = 0 (p < 0.0001) and q = 1 (p = 0.031), and for PD at all three orders (p = 0.00004, 0.018, 0.038).

- So per plant, forest roots hold more taxa than subpáramo or páramo roots, which don't differ from each other. The effect is clearest in richness (q = 0), still present at q = 1 but with no pairwise contrast resolved (best p = 0.055), and gone by q = 2. That points to rare taxa being lost. PD q = 1–2 has a supported interaction, which is issue 3.

**Family sensitivity** (`family_sensitivity.R`). The habitat effect and the three Tukey contrasts, refitted under all three families:

| Metric | Habitat p (gaussian / lognormal / gamma) | Forest vs subpáramo | Forest vs páramo | Subpáramo vs páramo |
|---|---|---|---|---|
| TD q = 0 | 0.00004 / 0.0012 / 0.0004 | 0.0003 / 0.005 / 0.0006 | 0.0002 / 0.003 / 0.0002 | ns in all |
| TD q = 1 | 0.037 / 0.157 / 0.048 | 0.062 / 0.18 / 0.055 | 0.071 / 0.28 / 0.078 | ns in all |
| TD q = 2 | 0.15 / 0.36 / 0.15 | ns in all | ns in all | ns in all |
| PD q = 0 | 0.00004 / 0.0008 / 0.0003 | 0.0006 / 0.005 / 0.001 | 0.0001 / 0.002 / 0.0001 | ns in all |
| PD q = 1 | 0.040 / 0.049 / 0.023 | ns in all | 0.047 / 0.057 / 0.047 | ns in all |
| PD q = 2 | 0.012 / 0.017 / 0.015 | 0.015 / 0.020 / 0.023 | 0.056 / 0.074 / 0.094 | ns in all |

- **The headline richness result is robust.** For TD q = 0 and PD q = 0 the habitat effect and both forest contrasts hold under every family, with fold changes barely moving (TD forest/subpáramo 1.52 under both lognormal and gamma; forest/páramo 1.65 vs 1.69). Subpáramo and páramo never differ.
- **One statement is family-sensitive.** For TD q = 1 the habitat effect is p = 0.157 under lognormal but 0.037 (gaussian) and 0.048 (gamma). No pairwise contrast resolves under any family (best p = 0.055), so "no pairwise habitat differences in abundance-weighted taxonomic diversity" holds throughout, whereas "no habitat effect at q = 1" does not. Worth wording as the former.
- **Family used: gamma with a log link.** Lowest summed AICc across all six, preferred or tied for four of them, gives ratios as reported, and avoids different families at different orders of q. `perplant_hill_size_standardised.R` and `parallelism_diagnostics.R` both have `FAMILY <- "gamma_log"` at the top (set both to `"lognormal"` to reproduce the earlier version). The tables above are the gamma results; the family table below shows the log-normal and gaussian equivalents.
- **Model structure is the same under every family, with one exception** (`family_sensitivity_structure.csv`): habitat + site for TD q = 0 and PD q = 0, site only for TD q = 2, and habitat × site for PD q = 1 and q = 2. The exception is TD q = 1, where habitat is retained under gaussian (ΔAICc −2.3) and gamma (−2.0) but not under lognormal (+0.9). The interaction is never supported for TD at any order, under any family.
- **What that changes at TD q = 1.** Under gamma, habitat stays in the model (p = 0.048) but no pairwise contrast resolves (forest/subpáramo 1.48, p = 0.055; forest/páramo 1.53, p = 0.078; subpáramo/páramo 1.03, p = 0.98), though the linear trend across habitats is significant (p = 0.031). So the per-plant decline from forest is clearest in richness, still present but unresolved pairwise at q = 1, and gone by q = 2, rather than being confined to richness alone.

| In the paper now | What the new results show |
|---|---|
| Methods: per-sample Hill diversity at q = 1 from `hillR` on each collapsed sample; interaction model; contrasts from the additive model | Per-plant values rarefied to 28,122 reads with iNEXT.3D at q = 0–2; every test comes from the additive model, with the interaction tested by model comparison. `hillR` no longer used |
| Results, second diversity paragraph: habitat F₂,₅₀ = 4.39 (TD) and 3.03 (PD), "confirming directional diversity decline"; none of the three Tukey contrasts resolved | Under the retained model (gamma, habitat + site), habitat affects per-plant richness (TD q = 0 F₂,₅₄ = 9.21, p = 0.0004; PD q = 0 F₂,₅₄ = 9.61, p = 0.0003): forest > subpáramo (1.52×, 1.39×) and forest > páramo (1.69×, 1.58×), subpáramo = páramo. At q = 1 habitat is retained but no pairwise contrast resolves (best p = 0.055); at q = 2 habitat is dropped |
| Tables S14–S15 (Type I ANOVA, q = 1) and S15b (contrasts, q = 1) | Replaced by `perplant_hill_size_std_omnibus.csv` and `perplant_hill_size_std_pairwise.csv`, all three orders of q |
| Figure 3B: unrarefied `hillR` q = 1 values by site | Rarefied values. The richness result (q = 0) is where the habitat effect is, so worth considering showing q = 0 as well |
| Appendix 5.3: narrative about the omnibus test and "no additional support for a stepwise decline" | Superseded; code and text need to match the new analysis |
| Appendix 5.3: family screening picks lognormal, "best model is lognormal_lm by far" | The screening compared AIC across response scales without a Jacobian correction, so the log-scale model was not comparable. Corrected and redone on the rarefied values, no single family wins across the board: gaussian for TD q = 0 and PD q = 0, gamma for TD q = 1, lognormal for TD q = 2, PD q = 1 and PD q = 2. Gamma with a log link is the best overall compromise (see the sensitivity results below) and keeps ratios, as currently reported |
| References: Li (2018) for `hillR`; iNEXT.3D cited as v1.0.6 | `hillR` no longer used; analyses run with iNEXT.3D 1.0.12 |

---

## 2. Habitat-level Hill numbers used root samples as independent units

**The problem.** The iNEXT3D habitat curves (Figure 3A, Tables S12–S13) were run on `tree_ps`, which keeps both roots per plant, so there were 50 / 40 / 30 "independent" sampling units. Two roots from one plant aren't independent, so the confidence intervals were too narrow. It also meant the Methods statement that roots were collapsed for Hill numbers wasn't true for this analysis. On top of that, the Results numbers match no table in the repo (e.g. Shannon 2,821 / 2,134 / 1,234 in the text against 2,927 / 2,342 / 1,316 in `Table_S16_iNEXT3D_TD.csv`), and the PD values are described as "Faith's PD in branch length units" when the analysis used meanPD (an effective number of lineages).

**What we did.** Re-ran iNEXT3D with plants as incidence units (an OTU is present in a plant if it's in either root), 25 / 20 / 15 plants, same tree, reference time, seed and 500 bootstraps. Added coverage-standardised estimates and pairwise habitat comparisons.

**What it showed.**

- Confidence intervals widen by up to about 2×.
- Páramo is roughly half as diverse as forest and subpáramo for every metric and order (ratios 1.46–2.56), with no CI overlap under either coverage-standardised or asymptotic estimation.
- Forest and subpáramo no longer differ for TD (coverage-standardised ratios 1.00–1.07, overlapping CIs) and differ only a little for PD (1.04–1.32). With root samples as units, every forest–subpáramo comparison had looked separated.
- Asymptotic richness from 15–25 plants runs away (forest 6,172, more than the 4,386 OTUs in the whole dataset), so the coverage-standardised values are the ones to trust.

Coverage-standardised estimates (60% coverage):

| Metric | Forest | Subpáramo | Páramo |
|---|---|---|---|
| TD q = 0 | 2,710 | 2,533 | 1,229 |
| TD q = 1 | 2,106 | 2,105 | 999 |
| TD q = 2 | 1,475 | 1,421 | 709 |
| PD (meanPD) q = 0 | 196 | 189 | 102 |
| PD q = 1 | 115 | 110 | 64 |
| PD q = 2 | 71 | 59 | 40 |

Putting issues 1 and 2 together: per plant, the drop is forest to subpáramo (forest > subpáramo = páramo); per habitat, it's subpáramo to páramo (forest ≈ subpáramo > páramo). That fits subpáramo plants being individually poorer but more different from each other, and páramo plants being as rich as subpáramo plants but more alike. This is qualitative only, since the per-plant values are abundance-based and the habitat values incidence-based.

| In the paper now | What the new results show |
|---|---|
| Methods: incidence-based Hill numbers across habitats; coverage-based rarefaction/extrapolation | Same approach but with plants as incidence units (25 / 20 / 15); coverage-standardised estimates at 60% coverage; habitats compared pairwise by CI overlap |
| Methods, collapsing sentence: roots collapsed for "Hill-number computation" | Now true for both the habitat-level and per-plant analyses |
| Results, first diversity paragraph: "consistent decline"; asymptotic richness 4,587 / 3,066 / 1,951; Shannon and Simpson values; "Faith's PD ~162 branch length units"; CIs "did not overlap among habitats for any diversity order or metric" | Current numbers don't match any table. Plant-level values in the table above. Páramo about half of forest and subpáramo (non-overlapping CIs); forest and subpáramo overlap for TD at all orders. PD is meanPD (effective lineages), not Faith's PD in branch length units. Asymptotic estimates unreliable at this sample size |
| Figure 3A and caption (title "Parallel decline…"; habitat level "across all sites") | Curves with plants as units: wider bands, forest and subpáramo overlapping, páramo clearly lower |
| Tables S12–S13 (root samples as units) | Plant-unit versions (`Table_S16/S17_*_plant.csv`), plus the coverage-standardised estimates |
| Discussion: habitat-level diversity "declined monotonically from forest to subpáramo to páramo" | Not monotonic at habitat level (forest ≈ subpáramo > páramo). Per plant it's forest > subpáramo = páramo |
| Appendix 5.2 (root-sample code; "Faith's PD" in places) | Needs to match the plant-unit analysis; meanPD throughout |

---

## 3. Abundance-weighted PD trajectories aren't parallel

**The problem.** The published per-plant models found no habitat × site interaction for TD or PD, and the Discussion uses that as the second of three lines of evidence for parallel restructuring. With PD on a common reference time (issue 1), the interaction is significant for abundance-weighted PD (q = 1: p = 0.008; q = 2: p < 0.001).

**What we did.** Leave-one-site-out tests, NV + DOM only (the two complete gradients), site-by-site trajectories, influential plants, and PD re-tested with Helotiales and Sebacinales dominance as covariates (`parallelism_diagnostics.R`).

**What it showed.**

- It's DOM. The interaction disappears only when DOM is dropped (q = 1: p = 0.18; q = 2: p = 0.23) and stays significant or borderline (p ≤ 0.054) when any other site is dropped. It's strong within NV + DOM alone (p = 0.012 and p = 0.0001).
- DOM forest roots have abundance-weighted PD of 2.21 (q = 1) against 1.43–1.44 in the other forests, then drop to subpáramo and páramo levels (DOM forest/subpáramo 1.69×, p = 0.0016; forest/páramo 1.81×, p = 0.0004). Every other site is flat.
- DOM forest is a single location (DOM_3, 3,224 m), so this could be one stand (issue 9).
- It isn't one odd plant: removing influential plants leaves it significant (p = 0.029 and p = 0.0004).
- It's partly Helotiales: Helotiales dominance strongly predicts lower abundance-weighted PD at q = 2 (p < 0.001, though not at q = 1, p = 0.18), and adding it as a covariate shrinks the q = 2 interaction by about 40% (F 7.27 to 4.51, still p = 0.0035).
- TD is parallel between the two complete gradients (NV + DOM, all p ≥ 0.43). Without NV the TD interaction reaches p = 0.02–0.04, which is weak given how many subsets were tested.
- AICc agrees strongly. For PD q = 1 the interaction model carries an Akaike weight of 0.89 (ΔAICc −4.8) and for q = 2 it carries 0.9995 (ΔAICc −15.7), while for TD at every order and PD q = 0 the additive model is preferred by 4.5 to 6.4 AICc units. This is the one interaction in the paper that both frameworks support outright.

| In the paper now | What the new results show |
|---|---|
| Abstract: "the magnitude and direction of community assembly were statistically consistent across mountaintops in composition, diversity trajectories and lineage-specific turnover" | Holds for TD (all orders) and PD richness (q = 0). Not for abundance-weighted PD (q = 1–2), which is site dependent because of DOM. Composition and lineages: issues 4 and 6 |
| Results, "Parallel restructuring" paragraph: no habitat × site interaction for TD (F₄,₅₀ = 1.52) or PD (F₄,₅₀ = 1.22); rate of decline didn't differ among sites | TD: no interaction at any q, including NV + DOM alone. PD q = 0: no interaction. PD q = 1: F₄,₅₀ = 3.81, p = 0.009; q = 2: F₄,₅₀ = 7.27, p = 0.0001, driven by DOM |
| Discussion, "Three independent lines of evidence", second line: no site dependence for TD or PD | True for TD and PD richness; abundance-weighted PD depends on site (DOM), and this is linked to Helotiales dominance |
| Appendix Figure 3 caption: "supporting parallel diversity erosion across regions" | Not supported for abundance-weighted PD |
| Conclusion: convergence "recovered across composition, diversity, lineage turnover, and phylogenetic structure" | Supported for magnitude of compositional turnover, TD, PD richness and NTI. Not for abundance-weighted PD, the identity of the taxa involved (issue 6) or lineage abundance shifts (issue 4) |

---

## 4. Dominant-lineage trends: parallelism was assumed, now tested

**The problem.** Three things. The Helotiales and Sebacinales beta regressions are `elevation + site`, which forces every site to share the same elevation slope, so "parallel lineage turnover" was built into the model rather than tested. For the same reason, the site lines in Figure 4A–B are parallel by construction. The coefficient is on `betareg`'s default logit link, so 0.101 per 100 m is not a 10% change in relative abundance. Exponentiated it is a 10.6% increase in μ/(1−μ), where μ is the mean proportion; on the proportion scale that is roughly +2.5 percentage points per 100 m at the mean Helotiales proportion of 0.43. Reporting the coefficient with its p-value plus a proportion-scale statement avoids the issue. And the Figure 4 caption says "individual root samples" when the models were fitted to the 60 plants.

**What we did.** Likelihood-ratio tests of elevation × site and habitat × site, site-specific slopes (`parallelism_diagnostics.R`, part B), and AICc model comparison (`model_selection_aicc.R`). Figure 4A–B redrawn two ways: site-specific slopes (option 1) or one site-averaged line (option 2).

**What it showed.**

| Order | Comparison | LR p | ΔAICc (interaction − additive) | Akaike weight of the interaction |
|---|---|---|---|---|
| Helotiales | elevation × site | 0.060 | +0.6 | 0.38 |
| Helotiales | habitat × site | 0.042 | +1.5 | 0.28 |
| Sebacinales | elevation × site | 0.33 | +4.6 | 0.06 |
| Sebacinales | habitat × site | 0.19 | +5.1 | 0.02 |

Site-specific elevation slopes (percentage points per 100 m):

| Order | BEL | DOM | MA | NV |
|---|---|---|---|---|
| Helotiales | +2.8 (ns) | +5.7 (2.9 to 8.5) | +1.0 (ns) | 0.0 (−2.6 to 2.6) |
| Sebacinales | −0.5 (p = 0.042) | −2.2 (p = 0.055) | 0.0 | −0.2 |

- The across-site trends hold and are supported by AICc: elevation earns its place for both orders (Helotiales ΔAICc −3.0 against site only, Sebacinales −1.5).
- **Site-specific slopes are not supported for either order** (ΔAICc +0.6 to +5.1), so the additive model is retained and the across-site trends are what gets reported. The parallelism that was previously built into the model has now been tested, and the additive model is the supported one. The original claim therefore stands, but on tested rather than assumed grounds.
- For Sebacinales, AICc goes further: habitat as a predictor is not supported at all (site only preferred, ΔAICc +2.5), so there is no support for a subpáramo peak.
- Elevation varies among locations, not plants, so the pooled p-values are on the optimistic side (issue 9).
- What still changes in the paper: the coefficients are odds, not relative abundance; the Figure 4 caption says root samples when the models are per plant; and the interaction test now needs reporting so that parallelism is stated as a result rather than an assumption.

| In the paper now | What the new results show |
|---|---|
| Abstract: "This parallel restructuring was driven by contrasting turnover of the two dominant lineages" | Stands. Parallelism was previously assumed by the model structure; tested now, site-specific slopes are not supported (ΔAICc +1.5 / +0.6) and the additive model is retained |
| Methods: beta regression with elevation and site | No test of whether slopes differ among sites; now tested (elevation × site LR test) |
| Results, lineage paragraph: Helotiales "increased in relative abundance by ~10.1% per 100 m"; Sebacinales "decreasing by ~10% per 100 m"; "parallel restructuring … concentrated within the two dominant lineages" | These are logit-scale coefficients, not changes in relative abundance: +0.101 and −0.100 per 100 m (p = 0.017, 0.033), which on the proportion scale is about +2.5 percentage points per 100 m for Helotiales at its mean of 0.43. The elevation × site and habitat × site comparisons are new and should be reported: neither is supported (ΔAICc +0.6 to +5.1), so the additive model stands for both orders |
| Figure 4A–B: one line per site from the additive model; "+10.1% per 100 m, p = 0.016" annotation; caption says "individual root samples" | The additive model is the supported one, so the figure structure is right. Use option 2 (one across-site line, points coloured by site), which shows the fitted model without implying four independently estimated site trends. Data are per plant (n = 60), not root samples; the annotation is a logit-scale coefficient, not a percentage change in abundance |
| Discussion opening: parallel restructuring "in composition, diversity trajectories, and lineage-specific turnover"; restructuring concentrated in Helotiales and Sebacinales with "a consistent set of environmental filters" | Stands for lineage turnover, now on tested grounds. Two qualifications: the taxa involved differ somewhat between the two complete gradients (issue 6) and abundance-weighted PD responds differently at DOM (issue 3) |
| Discussion, lineage paragraph: "≈10% per 100 m elevation"; "observed parallel trajectories are largely mediated through changes in these two groups" | Logit-scale coefficient, not a change in relative abundance. The parallel-trajectory statement stands under the retained additive model |
| Appendix 6.1–6.2: "+10.1% per 100 m", Figure 4 caption | Same issues as above |
| (no table) | New supplementary table of site slopes and interaction tests: `diag_lineage_site_slopes.csv`, `diag_lineage_parallelism_tests.csv`, plus the AICc comparison in `model_selection_aicc.csv` |

---

## 5. Genus-level responses (the GLLVM, Figure 4C, Data S2)

**The problem.**

- The GLLVM coefficients have blown up. Genus summaries reach |β| = 260 on the log scale (*Coniochaeta*), *Acephala* (one OTU) has β = 78, and Data S2 reports "% enrichment" up to 2 × 10¹¹⁵. That's what you get when a taxon is absent from the reference habitat (forest): the coefficient can't be estimated and just runs off. The heatmap squishes anything beyond ±35, so most tiles only show sign.
- The random effects (`row.eff = ~(1|site) + (1|Unique_ID)`) are community-level row effects, shared by every OTU. There are no OTU-specific site effects, so a taxon confined to MA (no forest) looks páramo-enriched whatever it actually does.
- The model comparison in the appendix (AIC 80,295 vs 186,832) can't be reproduced from `gllvm_analysis.R` (`fit_nb` is never defined; `fit_nb_1` uses a term not in `studyDesign`).
- *Meliniomyces* and *Hyaloscypha* are kept separate, with opposite signs, although we told Reviewer 1 they're synonyms.

**What we did.** Pooled reads to genus per plant (*Meliniomyces* merged into *Hyaloscypha*), kept genera in ≥ 10 plants with ≥ 1,000 reads (56 genera), and fitted a negative binomial GLMM per genus: `count ~ habitat + site + offset(log(library size)) + (1 | location)`. Site as a genus-specific fixed effect handles the confounding; the location random effect makes the 12 locations the replicates for habitat. Fold changes with CIs, BH q-values across genera, and genera absent from a habitat flagged rather than given a fold change (`genus_glmm_habitat.R`). Then split the Helotiales trend into genus-assigned and unassigned OTUs and by family (`helotiales_fraction_check.R`).

**What it showed.**

- Five genera have an omnibus habitat effect at q < 0.05 (*Acephala*, *Conlarium*, *Capnobotryella*, *Oidiodendron*, *Trichoderma*), and AICc supports habitat in nine (ΔAICc < −2). Separately, four genera have a significant contrast against forest, and **every one of them is a decline**. No genus is significantly more abundant in subpáramo or páramo than in forest.

Genera with a significant vs-forest contrast:

| Genus | Result | Plants present (F / S / P) | ΔAICc for habitat | NV, DOM (páramo vs forest) |
|---|---|---|---|---|
| *Oidiodendron* | Down in páramo, q = 0.014 | 17 / 8 / 1 | −6.3 | forest only, forest only |
| *Capnobotryella* | Down in subpáramo (q = 0.009) and páramo (q = 0.014) | see presence table | −6.4 | forest only, higher in forest |
| *Cephalotheca* | Down in subpáramo, q = 0.009 | see presence table | −3.2 | higher in páramo, forest only |
| *Gyoerffyella* | Down in páramo, q = 0.013 | 14 / 4 / 4 | −0.2 | higher in forest, higher in páramo |

*Gyoerffyella* is the one to treat carefully: the contrast is significant but neither the omnibus test (q = 0.26) nor AICc supports including habitat at all, which happens when a contrast is driven by a near-empty cell.

Genera where AICc supports habitat but no vs-forest contrast is estimable or significant: *Acephala* (ΔAICc −9.8, the strongest habitat support of any genus), *Conlarium* (−8.5), *Trichoderma* (−5.7), *Pezicula* (−4.5), *Cyphellophora* (−2.9), *Cladosporium* (−2.7).

- The genera we highlighted as páramo-enriched mostly don't hold up. *Hyaloscypha* (ΔAICc +4.2) and *Leohumicola* (+4.5) show no habitat effect on any measure, *Lachnum* is only borderline higher in subpáramo (contrast q = 0.059, ΔAICc −0.4), and *Serendipita* shows no significant contrast (subpáramo/forest 0.62, q = 0.86; páramo/forest 0.12, q = 0.18) and no AICc support (+0.8).
- *Acephala* is the exception and needs careful wording. It has the strongest habitat support of any genus (omnibus q = 0.029, ΔAICc −9.8), and no fold change is estimable only because it is never detected in forest. But it occurs at just two locations, DOM_1 (páramo) and DOM_2 (subpáramo), and is absent from DOM_3 (forest) and from NV entirely. So the habitat signal is real within one mountain and says nothing about the others. *Pezicula* is similar in spirit (ΔAICc −4.5) but goes in opposite directions at NV and DOM.
- The presence/absence patterns that are harder to put down to chance go the other way: genera found at forest locations on all three forested mountains but at no páramo location.

| Genus | Absent from | Where it occurs | Chance of that pattern |
|---|---|---|---|
| *Acephala* | forest | DOM_1, DOM_2 | 0.32 |
| *Pochonia* | páramo | all 5 forest locations, 2 subpáramo | 0.05 |
| *Tetracladium* | páramo | 4 forest, 2 subpáramo (BEL, DOM, NV) | 0.09 |
| *Cyphellophora* | páramo | 4 forest, 1 subpáramo (BEL, DOM, NV) | 0.16 |
| *Conlarium* | subpáramo | 4 forest, NV páramo | 0.07 |

(Rough guide only: treats locations as exchangeable and isn't adjusted for looking at several genera.)

- The order-level Helotiales increase is real, and it's carried by the OTUs with no genus assignment:

| Helotiales fraction | Share of Helotiales reads | exp(β) per 100 m (logit scale) | p | Site × elevation p |
|---|---|---|---|---|
| All (436 OTUs) | 100% | +10.6% | 0.017 | 0.060 |
| Genus unassigned (297) | 72% | +11.1% | 0.019 | 0.11 |
| Genus assigned (139) | 28% | +4.6% | 0.36 | 0.060 |
| Hyaloscyphaceae (43) | 10% | +11.1% | 0.021 | 0.018 (BEL only) |

- Site dependence is not supported in the unassigned fraction either (elevation × site p = 0.11). Hyaloscyphaceae is the one fraction with a significant interaction (p = 0.018; it increases at BEL only), and it holds 10% of Helotiales reads. No other family trends.
- These are relative abundances. Helotiales is ~43% of reads per plant, so its share goes up when Sebacinales and the forest genera go down, and the data can't separate a Helotiales gain from losses elsewhere.

So: Helotiales becomes more dominant upslope, carried by lineages we can't resolve to genus (no site-specific slopes are supported, issue 4). At genus level the signal is forest-associated taxa dropping out upslope. The one genus with strong statistical support for higher-elevation occurrence, *Acephala*, occurs on a single mountain, so it cannot carry a general specialist claim.

| In the paper now | What the new results show |
|---|---|
| Abstract: Helotiales increase "driven by a subset of high-elevation specialist taxa" | No genus is significantly enriched above the treeline. The increase is carried by genus-unassigned Helotiales (72% of Helotiales reads, exp(β) = 1.11 per 100 m, p = 0.019); genus-assigned Helotiales show no trend |
| Methods: GLLVM with habitat fixed, random intercepts for site and plant, one latent variable | Replaced by per-genus negative binomial GLMMs (habitat + site, location random effect, BH across 56 genera); taxa absent from a habitat reported by presence |
| Results, genus paragraph: increase "driven by a subset of specialized taxa"; strongest páramo enrichment in *Acephala*, *Hyaloscypha* and *Pezicula* ("large positive β"); *Lachnum* and *Leohumicola* strongest in subpáramo; *Serendipita* enriched in subpáramo, declined sharply in páramo | Every significant vs-forest contrast is a decline (*Oidiodendron*, *Capnobotryella*, *Cephalotheca*, *Gyoerffyella*). *Acephala* has the strongest habitat support of any genus but occurs at two DOM locations only; *Hyaloscypha* and *Leohumicola* show no effect at all; *Pezicula* has AICc support but opposite directions at NV and DOM; *Lachnum* borderline; *Serendipita* no significant contrast. Forest genera (*Pochonia*, *Tetracladium*, *Cyphellophora*) absent from all páramo locations |
| Discussion, Helotiales paragraph: *Hyaloscypha* and *Acephala* "consistently enriched in páramo roots"; *Lachnum* and *Leohumicola* strongest in subpáramo; transition "selects for a small number of high-elevation specialists within Helotiales" | "Consistently" is not supported: *Hyaloscypha* and *Leohumicola* show no habitat effect, and *Acephala*, which does, is confined to one mountain. Helotiales dominance rises on average via unresolved lineages; the repeatable genus-level signal is loss of forest taxa upslope, including the ErM genus *Oidiodendron* |
| Discussion, Sebacinales paragraph: relative abundance "peaked in subparamos … declining sharply in páramo"; subpáramo is "the upper elevational limit of Sebacinales dominance, above which Helotiales take over" | No evidence of a subpáramo peak: *Serendipita* subpáramo/forest 0.62 (ns), and at order level AICc does not support habitat as a predictor of Sebacinales at all (site only preferred). Sebacinales declines with elevation on average (p = 0.033), but the decline is small and only detectable at BEL and DOM |
| Figure 4C and caption: GLLVM β heatmap, 15 hand-picked genera, colour limits ±35, *Meliniomyces* and *Hyaloscypha* separate | New heatmap from the GLMMs: log₂ fold change vs forest, stars for q, "abs. F" where a genus isn't in forest, *Meliniomyces* merged (`objects/p_heat_glmm.rds`) |
| Data S2: GLLVM genus table ("% enrichment" up to 2 × 10¹¹⁵) | Replaced by `genus_glmm_results.csv` and `genus_glmm_presence_by_site.csv` |
| Appendix 6.3: GLLVM model selection (AIC comparison), diagnostics | AIC comparison can't be reproduced; section needs to describe the GLMMs |
| References: Niku et al. (2019) ×2 | GLLVM no longer used; glmmTMB (Brooks et al., 2017) is |
| (no table) | New supplementary table of the Helotiales fractions (`helotiales_fraction_check.csv`) |

---

## 6. Composition: unreported dispersion result, unadjusted habitat test, and a site-specific component

**The problem.** Three things, none of which overturn the composition results. The centroid model's habitat and site tests are Type I, so habitat is not adjusted for site (the interaction test, entered last, is fine). The betadisper result is not reported anywhere, although habitats differ significantly in dispersion (F = 11.3, p = 0.001), which matters for how the PERMANOVA habitat effect is read. And the PERMANOVA on the two complete gradients has a significant habitat × site interaction (R² = 0.09 against a habitat main effect of 0.10, p = 0.001), which is worth reporting alongside the centroid result because it picks up a site-specific component the centroid model cannot see.

**What we did.** Site-adjusted tests for the centroid model, Tukey contrasts, and a location-means version (n = 12) (`adjusted_tests_centroid_nri.R`, using the standard robust CLR).

**What it showed.**

- Habitat adjusted for site: F₂,₅₄ = 6.11, p = 0.004 (published, unadjusted: F₂,₅₀ = 14.7).
- Forest plants sit further from their site centroid than subpáramo (difference 3.6, p = 0.008) or páramo plants (3.9, p = 0.015); subpáramo and páramo don't differ.
- Site adjusted for habitat: F₃,₅₄ = 3.52, p = 0.021. The interaction doesn't change (p = 0.29).
- On location means (n = 12, 6 residual df) neither habitat (p = 0.12) nor site (p = 0.26) is detectable.
- Forest is also the most dispersed habitat (betadisper mean distance 23.6 vs 17.8 subpáramo and 17.1 páramo), so the centroid "turnover" metric is partly picking up heterogeneity within forest.
- The composition evidence points both ways and both parts are real. The dbRDA conditioned on site finds a habitat effect after site differences are removed (F = 1.41, p = 0.001), and the centroid model finds that the size of the habitat-associated shift does not differ among sites (interaction p = 0.29). Against that, the PERMANOVA interaction says the habitat effect on composition is not identical between NV and DOM. A shared compositional response with a site-specific component on top is the fair summary, rather than consistency or divergence alone.

| In the paper now | What the new results show |
|---|---|
| Abstract: "the magnitude and direction of community assembly were statistically consistent across mountaintops in composition, diversity trajectories and lineage-specific turnover" | Supported for diversity (the decline with elevation is the same size and direction at every site for TD at all orders and for PD richness; abundance-weighted PD is the exception, issue 3) and for the magnitude of compositional turnover (centroid interaction p = 0.29). Worth adding one qualification: the habitat effect on composition is not identical between the two complete gradients (PERMANOVA interaction R² = 0.09, p = 0.001) |
| Discussion, "Three lines of evidence": "Isolated páramo regions therefore retained distinct community baselines and species pools, but the magnitude and direction of community restructuring, in both compositional identity and diversity, were statistically consistent across sky islands" | Largely stands. The distinct-baselines clause is fine, and the diversity half holds for TD at all orders and PD richness (abundance-weighted PD is the exception, issue 3). The one phrase worth softening is "in both compositional identity and diversity": the taxa involved are not identical between NV and DOM (PERMANOVA interaction R² = 0.09, p = 0.001), even though the magnitude of turnover is consistent |
| Results, centroid paragraph: habitat F₂,₅₀ = 14.7, p < 0.001; site F₃,₅₀ = 3.6, p = 0.020 | Adjusted for site: habitat F₂,₅₄ = 6.11, p = 0.004 (forest further from centroid than subpáramo and páramo); site F₃,₅₄ = 3.52, p = 0.021. Interaction unchanged |
| Methods: "Homogeneity of multivariate dispersion was assessed with betadisper" (no result reported) | Dispersion differs among habitats (F = 11.3, p = 0.001), highest in forest, so part of the PERMANOVA habitat effect may be spread rather than location |
| Table S10 (Type I) | Add the adjusted tests and Tukey contrasts (`adjusted_tests_centroid_nri.csv`, `adjusted_pairwise_centroid_nri.csv`) |

---

## 7. NRI/NTI

**The problem.** The appendix code computes NRI/NTI on `tree_ps` (120 root samples), but Tables S19–S20 have n = 25 / 20 / 15 and 54 residual df, so they were made from the plants and the code doesn't match the output. The habitat test is Type I (unadjusted); the site test, entered second, is adjusted.

**What we did.** Recomputed NRI/NTI on `ps_individual` with the same tree and settings, then site-adjusted tests, interaction, and location means.

**What it showed.** Medians match Table S19 within null-model noise (e.g. forest NRI_aw 1.06 vs 0.99). No habitat effect on either index adjusted for site (NRI p = 0.48; NTI p = 0.43), and AICc agrees: for both indices the site-only model is preferred over habitat + site (ΔAICc +3.4 and +3.2). The borderline NRI habitat × site term (p = 0.053) gets no AICc support either (ΔAICc +0.35 against the additive model, and both lose to site only). NRI site effect holds at plant level (F₃,₅₄ = 3.29, p = 0.027) and NTI still has none (p = 0.19), so the deep vs tip-level decoupling argument survives. On location means the NRI site effect isn't detectable (p = 0.20).

| In the paper now | What the new results show |
|---|---|
| Results, phylogenetic structure: habitat NTI p = 0.60, NRI p = 0.29; site NTI p = 0.16, NRI p = 0.023 | Adjusted: habitat NTI p = 0.43, NRI p = 0.48; site NTI p = 0.19, NRI p = 0.027. Same conclusions |
| Discussion, NRI/NTI paragraph: same p-values | As above; NRI site effect not detectable on location means (p = 0.20) |
| Table S20 (Type I) | Adjusted tests |
| Appendix 8 code uses `tree_ps` | Tables were made from plants; code needs to use `ps_individual` with the tree |

---

## 8. vegan changed the robust CLR (reproducibility)

**The problem.** From vegan 2.7-1, `decostand(method = "rclr")` and `vegdist(method = "robust.aitchison")` fill in zeros by matrix completion (optspace) by default. That isn't the robust CLR we describe (log of the non-zero counts centred on their mean, zeros left as zeros), and it isn't what produced our tables. Under vegan 2.7.2 the appendix code gives distances that are ~45% different on average and completely different results (e.g. dbRDA F = 6.15 instead of 1.41, centroid habitat F = 0.82 instead of 14.7). The same release drops row names from the rclr output, which is why the appendix's `stopifnot` on row names fails. The dbRDA chunk also passes `otu_table(ps)` to `vegdist()`, which now errors. Because the appendix chunks are `eval: false`, the rendered appendix just reads the old CSVs and none of this shows up.

**What we did.** Computed the robust CLR explicitly without imputation and re-ran the composition analyses (`reproduce_composition_tables.R`).

**What it showed.** Tables S5–S10 reproduce exactly. The published composition results and Figure 2 are fine; the problem is only that the code no longer reproduces them under current vegan. Any other project using `rclr` or `robust.aitchison` with default settings will shift in the same way if re-run.

| In the paper now | What the new results show |
|---|---|
| Methods: Aitchison distances from rCLR-transformed abundances (vegan) | Published results used the zeros-retained rCLR; vegan ≥ 2.7-1 imputes zeros by default, so the vegan version matters |
| Appendix 4.1–4.3 code (`decostand(..., "rclr")`, `vegdist(..., "robust.aitchison")`, `vegdist(otu_table(ps))`) | Doesn't reproduce the tables under current vegan, and the dbRDA chunk errors |
| Appendix Section 0/1: no package versions | Versions needed for anyone to reproduce the tables |

---

## 9. Design caveat: 12 locations, mostly one per habitat per site

Only BEL forest and NV forest have two locations; every other habitat × site cell is one location of five plants. That means habitat × site interactions can't be separated from differences between individual locations (the DOM PD pattern is one forest stand, DOM_3), and plant-level analyses treat plants within a location as independent replicates of habitat, so p-values for habitat and elevation are on the optimistic side. When run on location means (n = 12), the centroid habitat effect (p = 0.12) and the NRI site effect (p = 0.20) aren't detectable. The genus GLMMs already include location as a random effect.

| In the paper now | What the new results show |
|---|---|
| Discussion, limitations: only NV and DOM complete; interaction terms from an incomplete design with reduced power | Also: most cells are single locations, so interactions and location effects are confounded; the DOM patterns (PD, Helotiales) may be one stand; plant-level p-values optimistic, and some effects disappear on location means |

---

## 10. Appendix housekeeping

- Everything is `eval: false`, so rendering never re-checks the tables. Worth saying so in Section 1 along with the package versions.
- Code that doesn't match output or doesn't run: dbRDA (`vegdist(otu_table(ps))`), rclr/robust Aitchison (issue 8), NRI/NTI (`tree_ps`, issue 7), `gllvm_analysis.R` (undefined objects).
- The iNEXT3D chunk still contains a redundant tree-inference block (its own comment says so).
- GLLVM diagnostics captioned Figure S6 but saved as `Fig_S8_…`; NTI diagnostics captioned S8 but saved as `Fig_S9_…`.
- CSV names don't match displayed table numbers (e.g. `Table_S16_iNEXT3D_TD.csv` is shown as Table S12).
- Methods say soil counts were subtracted from "every root sample"; the code does it per root PCR replicate (the appendix gets this right).
- "GTR+G+I" and "GTR+Γ+I" both appear in the Methods.

---

## 11. Response to reviewers

| What the letter says now | What the new results show |
|---|---|
| Opening summary, point (iv): post-hoc contrasts added (Table S15b) and the habitat-level vs per-sample distinction carried through | The post-hoc work is now bigger: habitat-level diversity redone with plants as units, per-plant diversity rarefied, all habitat tests site-adjusted from one model, lineage parallelism tested, genus-level analysis replaced |
| R2, minor point 1 (Abstract): results report "no divergence in the magnitude or direction of community assembly was detectable across mountaintops" | Mostly holds. Two qualifications to add: abundance-weighted PD responds to habitat differently at DOM (issue 3), and the habitat effect on composition is not identical between the two complete gradients (issue 6) |
| R2, key point 2 (collapsing): GLLVMs used the full sampling structure, "correctly specifying random effects" | GLLVM random effects were community-level only; the replacement GLMMs work on plants with location as a random effect |
| R2, key point 2: list of analyses run on collapsed data ("…Hill-number computation and NRI/NTI") | Accurate now (issues 2 and 7) |
| R2, p. 10 l. 24: per-sample Hill numbers "operate on summed raw counts with library size handled within each analysis" | Per-plant Hill numbers weren't depth-standardised; now rarefied to 28,122 reads |
| R2, p. 15 l. 35 (post-hoc tests): no response yet | Issue 1 has the answer: the old omnibus and post-hoc tests came from different models; everything now comes from one site-adjusted model, and pairwise tests exist at both the per-plant and habitat level |
| R2, p. 17 l. 26–34: quotes the Discussion sentence on selection for "a small number of high-elevation specialists within Helotiales" | That claim isn't supported (issue 5) |
| R2, minor point 3 (limitations): quotes the limitations paragraph | The limitations need to cover the single-location issue (issue 9) |
| R1, "Within Sebacinales, Serendipita": adopted sentence says *Serendipita* "was enriched in subpáramo relative to forest, but declined sharply in páramo roots" | Not supported by the genus GLMM (issue 5) |
| R1, genus list: *Meliniomyces* treated as a synonym of *Hyaloscypha* | True in the text, but the GLLVM heatmap kept them separate; now merged |
| R1 and R2, Figure 4A–B: taxon labels added to the y-axes | Lines and the "+10.1% per 100 m" annotation change too (issue 4) |
| R1, Figure 3: Results paragraph now cites Fig. 3A, "consistent decline" wording | That paragraph changes substantially (issue 2) |
| R1, soil subtraction: "subtracted from the raw counts of every root sample" | It's every root PCR replicate |

---

## Decisions still to make

| Open question | What the results say |
|---|---|
| Figure 4A–B: option 1 (site-specific slopes) or option 2 (one across-site line)? | Settled: option 2, since the additive model is the retained one. Option 1 is no longer needed |
| Whether to report gamma or log-normal per-plant models | Settled on gamma (lowest summed AICc, keeps ratios); the script has a `FAMILY` switch. The only consequence is that habitat is retained at TD q = 1, where no pairwise contrast resolves |
| Figure 3B: keep q = 1 only (as published) or add q = 0? | The per-plant habitat effect is in richness (q = 0); q = 1 TD shows no habitat effect, and q = 1 PD shows the DOM interaction |
| Which new supplementary tables, and their numbers | Candidates: coverage-standardised habitat estimates (with S12–S13), adjusted centroid tests (with S10), lineage site slopes and interaction tests, Helotiales fractions, genus GLMM results (Data S2) |
| Report location-mean analyses, or only mention them in the limitations? | They matter for the centroid habitat effect (p = 0.12) and the NRI site effect (p = 0.20); nothing else was run at location level |
| How far to rework the "parallel" framing | Parallel holds for turnover magnitude, TD, PD richness and NTI; not for abundance-weighted PD, taxa identity or lineage abundance. Title can stay |
| How much of the genus-level detail to report | Four genera have significant vs-forest contrasts, five have omnibus effects, nine have AICc support, and the three sets only partly overlap. Worth settling on one criterion for the Results |

---

## Appendix: what replaces what

| Appendix section now | What replaces it |
|---|---|
| Section 1 (object loading; analyses described as run on the server) | Add package versions, and say that analysis chunks are `eval: false` and the tables are read from saved outputs |
| 4.1–4.3 composition code (`decostand(..., "rclr")`, `vegdist(..., "robust.aitchison")`, `vegdist(otu_table(ps))`) | Explicit robust CLR without imputation, as in `reproduce_composition_tables.R`; count matrix passed to distance functions; adjusted centroid tests and Tukey contrasts from `adjusted_tests_centroid_nri.R` |
| 5.2 iNEXT3D (root samples as units) | `iNEXT3D_plant_level.R`: `Table_S16_iNEXT3D_TD_plant.csv`, `Table_S17_iNEXT3D_PD_plant.csv`, coverage-standardised values from `iNEXT3D_estimates_plant_vs_root.csv`, pairwise comparisons from `iNEXT3D_pairwise_habitat_plant_vs_root.csv` |
| 5.3 per-sample Hill diversity (`hillR`, Type I ANOVA, S15b) | `perplant_hill_size_standardised.R`: `perplant_hill_size_std_omnibus.csv`, `_pairwise.csv`, `_emmeans.csv`, `_trend.csv`; the Type I vs adjusted diagnostic (`perplant_hill_published_typeI_vs_adjusted.csv`) is worth including as it explains the change |
| 5.4 Figure 3 | `figures/Figure_3_corrected.png` (or panels `Figure_3A_iNEXT_plant_units.png`, `Figure_3B_perplant_size_std.png`) |
| (new) parallelism of PD and lineages | `parallelism_diagnostics.R`: `diag_leave_one_site_out.csv`, `diag_site_trajectories.csv`, `diag_within_site_contrasts.csv`, `diag_lineage_parallelism_tests.csv`, `diag_lineage_site_slopes.csv`, `diag_PD_interaction_vs_lineages.csv` |
| 6.1 beta regressions | Keep, plus the elevation × site tests and site slopes above; fix the "+10.1% per 100 m" wording |
| 6.2 Figure 4 | `Figure4_option1_site_slopes.R` or `Figure4_option2_pooled.R` |
| 5.3 family screening (AIC across response scales) | `model_selection_aicc.R` section E (Jacobian-corrected) and `family_sensitivity.R` |
| (new) model structure | `model_selection_aicc.R`: `model_selection_aicc.csv`, `model_selection_aicc_genus.csv`, `model_selection_family_screen.csv` |
| 6.3 GLLVM (model selection, diagnostics, genus table) | `genus_glmm_habitat.R` (`genus_glmm_results.csv`, `genus_glmm_presence_by_site.csv`, heatmap) and `helotiales_fraction_check.R` (`helotiales_fraction_check.csv`) |
| 8 NRI/NTI (`tree_ps`, Type I) | `adjusted_tests_centroid_nri.R`: `nri_nti_per_plant.csv`, `nri_nti_median_iqr_per_plant.csv`, `adjusted_tests_centroid_nri.csv`, `adjusted_pairwise_centroid_nri.csv` |

---

## Where "parallel" appears in the manuscript

| In the paper now | What the new results show |
|---|---|
| Title and running title | Fine: composition magnitude and TD are parallel |
| Abstract: prediction sentence | Fine |
| Abstract: results sentences on consistent magnitude and direction, and lineage turnover driving parallel restructuring | Broadly supported; qualifications in issues 3 and 6 |
| Introduction: research question (2) | Fine |
| Methods: "to test whether diversity trajectories were parallel" | Fine as a question; the method underneath changes (issue 1) |
| Figure 3 caption title: "Parallel decline…" | Not supported for abundance-weighted PD; habitat-level decline isn't a steady decline (issue 2) |
| Results: "Abundance-sensitive metrics showed parallel reductions" | Paragraph changes (issue 2) |
| Results: "Parallel restructuring across isolated sky islands" section | Holds for turnover magnitude and TD; not PD q = 1–2 (issue 3) |
| Results: lineage paragraph opening "The parallel restructuring…" | Issue 4 |
| Discussion: opening paragraph and "Three independent lines of evidence" | Issues 3, 4, 6 |
| Discussion: "parallel assembly still emerges" (host filter paragraph) | Fine with the qualifications above |
| Discussion: limitations | Issue 9 |
| Discussion: lineage paragraph "parallel trajectories are largely mediated through these two groups" | Issue 4 |
| Conclusion | Issue 3 |
