#!/usr/bin/env Rscript
# =============================================================================
# model_selection_aicc.R
#
# Information-criterion comparison for every model-structure decision made in
# the revision, so that the choice between additive and interaction models does
# not rest on p-values alone.
#
# What the published analysis actually used AIC for:
#   - per-plant Hill models: family screening only (gaussian vs lognormal vs
#     gamma), plain AIC, appendix Section 5.3;
#   - GLLVM: distribution and latent-variable structure, plain AIC.
#   Interaction structure was decided with anova() p-values, not AIC, and AICc
#   was not used anywhere.
#
# This script reports AIC, AICc, BIC, delta AICc and Akaike weights for:
#   A. per-plant size-standardised Hill numbers (TD and PD, q = 0, 1, 2):
#        site-only vs habitat + site vs habitat * site
#   B. centroid distance, NRI_aw, NTI_aw: same three-model set
#   C. Helotiales and Sebacinales beta regressions:
#        site-only vs elevation + site vs elevation * site, and
#        site-only vs habitat + site vs habitat x site (cell means)
#   D. genus GLMMs: delta AICc for including habitat, alongside the LR q-values
#   E. family screening for the per-plant models, redone on the rarefied values
#      with AICc (gaussian vs lognormal vs gamma)
#
# k is taken from logLik (so it includes sigma / the precision parameter /
# variance components), and aliased coefficients in the interaction models are
# not counted. AICc = AIC + 2k(k+1)/(n-k-1).
#
# Run from the repository root, after perplant_hill_size_standardised.R and
# adjusted_tests_centroid_nri.R:
#   Rscript Scripts/model_selection_aicc.R
#
# Inputs : objects/perplant_hill_size_std.rds, objects/ps_individual.rds,
#          tables/nri_nti_per_plant.csv, tables/genus_glmm_results.csv
# Outputs: tables/model_selection_aicc.csv        (A-C, one row per model)
#          tables/model_selection_aicc_genus.csv  (D)
#          tables/model_selection_family_screen.csv (E)
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(betareg); library(glmmTMB); library(dplyr); library(tidyr)
})

out_tab_dir <- "tables"; out_obj_dir <- "objects"
dir.create(out_tab_dir, showWarnings = FALSE)
hab_levels <- c("forest", "subparamo", "paramo")
fmt <- function(x) signif(x, 4)

# ---- Information criteria ---------------------------------------------------
ic <- function(m, n) {
  ll <- logLik(m)
  k  <- attr(ll, "df")
  aic <- -2 * as.numeric(ll) + 2 * k
  aicc <- if (n - k - 1 > 0) aic + 2 * k * (k + 1) / (n - k - 1) else NA_real_
  data.frame(logLik = as.numeric(ll), k = k, n = n, AIC = aic, AICc = aicc,
             BIC = -2 * as.numeric(ll) + log(n) * k)
}

compare <- function(models, n, comparison) {
  tab <- bind_rows(lapply(names(models), function(nm) cbind(model = nm, ic(models[[nm]], n))))
  tab$dAICc <- tab$AICc - min(tab$AICc, na.rm = TRUE)
  w <- exp(-0.5 * tab$dAICc); tab$weight_AICc <- w / sum(w)
  tab$dAIC  <- tab$AIC  - min(tab$AIC)
  tab$dBIC  <- tab$BIC  - min(tab$BIC)
  cbind(comparison = comparison, tab)
}

# Verdict for a nested pair, using the usual delta < 2 rule of thumb
verdict <- function(tab, simple, complex) {
  d <- tab$AICc[tab$model == complex] - tab$AICc[tab$model == simple]
  if (is.na(d)) return("not computable")
  if (d <= -2) "interaction clearly better"
  else if (d < 0) "interaction better by < 2 (equivalent support; keep simpler)"
  else if (d < 2) "additive better by < 2 (equivalent support; keep simpler)"
  else "additive clearly better"
}

all_tabs <- list(); verdicts <- list()
add_set <- function(models, n, label, simple = "habitat + site", complex = "habitat * site") {
  tab <- compare(models, n, label)
  all_tabs[[length(all_tabs) + 1]] <<- tab
  verdicts[[length(verdicts) + 1]] <<- data.frame(
    comparison = label,
    dAICc_interaction_vs_additive = tab$AICc[tab$model == complex] - tab$AICc[tab$model == simple],
    dAICc_additive_vs_null = tab$AICc[tab$model == simple] - tab$AICc[tab$model == names(models)[1]],
    best_by_AICc = tab$model[which.min(tab$AICc)],
    verdict = verdict(tab, simple, complex))
  invisible(tab)
}

# =============================================================================
# A. Per-plant size-standardised Hill numbers
# =============================================================================
vals <- readRDS(file.path(out_obj_dir, "perplant_hill_size_std.rds")) %>%
  mutate(habitat = factor(as.character(habitat), levels = hab_levels),
         site    = factor(as.character(site)))

grid <- distinct(vals, diversity, q) %>% arrange(desc(diversity), q)
for (i in seq_len(nrow(grid))) {
  d <- filter(vals, diversity == grid$diversity[i], q == grid$q[i])
  add_set(list("site only"       = lm(log(estimate) ~ site, data = d),
               "habitat + site"  = lm(log(estimate) ~ habitat + site, data = d),
               "habitat * site"  = lm(log(estimate) ~ habitat * site, data = d)),
          n = nrow(d), label = paste0("per-plant ", grid$diversity[i], " q=", grid$q[i]))
}

# =============================================================================
# B. Centroid distance and phylogenetic structure
# =============================================================================
nri_path <- file.path(out_tab_dir, "nri_nti_per_plant.csv")
if (file.exists(nri_path)) {
  nri <- read.csv(nri_path, stringsAsFactors = FALSE) %>%
    mutate(habitat = factor(habitat, levels = hab_levels), site = factor(site))
  for (resp in c("dist_centroid", "NRI_aw", "NTI_aw")) {
    d <- nri[!is.na(nri[[resp]]), ]
    f <- function(rhs) as.formula(paste(resp, "~", rhs))
    add_set(list("site only"      = lm(f("site"), data = d),
                 "habitat + site" = lm(f("habitat + site"), data = d),
                 "habitat * site" = lm(f("habitat * site"), data = d)),
            n = nrow(d), label = resp)
  }
} else {
  message("tables/nri_nti_per_plant.csv not found; run adjusted_tests_centroid_nri.R first. Skipping section B.")
}

# =============================================================================
# C. Dominant-lineage beta regressions
# =============================================================================
ps  <- readRDS(file.path(out_obj_dir, "ps_individual.rds"))
otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)[colnames(otu), ]
ordv <- tolower(sub("^[a-z]__", "", as.character(tax$Order))); ordv[is.na(ordv)] <- ""
sd0 <- as(sample_data(ps), "data.frame")[rownames(otu), ]

lin <- data.frame(site = factor(as.character(sd0$site)),
                  habitat = factor(as.character(sd0$habitat), levels = hab_levels),
                  elev_100 = as.numeric(sd0$elevation) / 100)
lin$cell <- interaction(lin$site, lin$habitat, drop = TRUE)
nL <- nrow(lin); total <- rowSums(otu)
for (o in c("helotiales", "sebacinales")) {
  p <- colSums(t(otu[, grepl(o, ordv, fixed = TRUE), drop = FALSE])) / total
  lin[[o]] <- (p * (nL - 1) + 0.5) / nL
}

for (o in c("helotiales", "sebacinales")) {
  f <- function(rhs) as.formula(paste(o, "~", rhs))
  add_set(list("site only"        = betareg(f("site"), data = lin),
               "elevation + site" = betareg(f("elev_100 + site"), data = lin),
               "elevation * site" = betareg(f("elev_100 * site"), data = lin)),
          n = nL, label = paste(o, "(elevation)"),
          simple = "elevation + site", complex = "elevation * site")
  add_set(list("site only"      = betareg(f("site"), data = lin),
               "habitat + site" = betareg(f("habitat + site"), data = lin),
               "habitat * site" = betareg(f("cell"), data = lin)),
          n = nL, label = paste(o, "(habitat)"))
}

model_tab <- bind_rows(all_tabs)
verdict_tab <- bind_rows(verdicts)
write.csv(model_tab, file.path(out_tab_dir, "model_selection_aicc.csv"), row.names = FALSE)

# =============================================================================
# D. Genus GLMMs: is habitat worth including?
# =============================================================================
gres_path <- file.path(out_tab_dir, "genus_glmm_results.csv")
if (file.exists(gres_path)) {
  gres <- read.csv(gres_path, stringsAsFactors = FALSE)
  gen <- sub("^[a-z]__", "", as.character(tax$Genus)); gen[is.na(gen)] <- ""
  gen[gen == "Meliniomyces"] <- "Hyaloscypha"
  keep_otu <- gen != "" & !grepl("incertae_sedis", gen, ignore.case = TRUE) &
    !grepl(" (Kingdom|Phylum|Class|Order|Family)$", gen) & !grepl("_gen_", gen, fixed = TRUE)
  gcounts <- t(rowsum(t(otu[, keep_otu, drop = FALSE]), group = gen[keep_otu]))

  meta <- data.frame(site = factor(as.character(sd0$site)),
                     habitat = factor(as.character(sd0$habitat), levels = hab_levels),
                     location = factor(as.character(sd0$site_elevation)),
                     lib = total)
  genus_ic <- bind_rows(lapply(gres$genus, function(g) {
    if (!g %in% colnames(gcounts)) return(NULL)
    d <- meta; d$y <- gcounts[, g]
    q <- function(f) suppressWarnings(tryCatch(
      glmmTMB(f, family = nbinom2, data = d), error = function(e) NULL))
    m1 <- q(y ~ habitat + site + offset(log(lib)) + (1 | location))
    m0 <- q(y ~ site + offset(log(lib)) + (1 | location))
    if (is.null(m1) || is.null(m0)) return(data.frame(genus = g, dAICc_habitat = NA_real_))
    data.frame(genus = g,
               dAICc_habitat = ic(m1, nrow(d))$AICc - ic(m0, nrow(d))$AICc)
  })) %>%
    left_join(select(gres, genus, order, omnibus_p, omnibus_q, published_fig4C,
                     paramo_vs_forest_q, subparamo_vs_forest_q), by = "genus") %>%
    mutate(habitat_supported_AICc = dAICc_habitat < -2) %>%
    arrange(dAICc_habitat)
  write.csv(genus_ic, file.path(out_tab_dir, "model_selection_aicc_genus.csv"), row.names = FALSE)
} else {
  genus_ic <- NULL
  message("tables/genus_glmm_results.csv not found; run Figure4_pooled_glmm.R first. Skipping section D.")
}

# =============================================================================
# E. Family screening on the rarefied per-plant values (appendix 5.3 redone)
# =============================================================================
fam_screen <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  d <- filter(vals, diversity == grid$diversity[i], q == grid$q[i])
  fits <- list(gaussian     = lm(estimate ~ habitat + site, data = d),
               lognormal_lm = lm(log(estimate) ~ habitat + site, data = d),
               gamma_log    = glm(estimate ~ habitat + site, data = d, family = Gamma(link = "log")))
  tab <- bind_rows(lapply(names(fits), function(nm) cbind(family = nm, ic(fits[[nm]], nrow(d)))))
  # lognormal is on a different response scale: add the Jacobian so AIC is comparable
  jac <- sum(log(d$estimate))
  tab$AIC[tab$family == "lognormal_lm"]  <- tab$AIC[tab$family == "lognormal_lm"]  + 2 * jac
  tab$AICc[tab$family == "lognormal_lm"] <- tab$AICc[tab$family == "lognormal_lm"] + 2 * jac
  tab$dAICc <- tab$AICc - min(tab$AICc)
  cbind(metric = paste0(grid$diversity[i], " q=", grid$q[i]), tab)
}))
write.csv(fam_screen, file.path(out_tab_dir, "model_selection_family_screen.csv"), row.names = FALSE)

# ---- Console summary --------------------------------------------------------
cat("\n==== A-C. Interaction vs additive (negative dAICc favours the interaction) ====\n")
print(as.data.frame(verdict_tab %>%
  transmute(comparison,
            dAICc_int_vs_add = fmt(dAICc_interaction_vs_additive),
            dAICc_add_vs_null = fmt(dAICc_additive_vs_null),
            best_by_AICc, verdict)), row.names = FALSE)

cat("\n==== Full model table ====\n")
print(as.data.frame(model_tab %>%
  transmute(comparison, model, k, AICc = fmt(AICc), dAICc = fmt(dAICc),
            weight_AICc = fmt(weight_AICc), dAIC = fmt(dAIC), dBIC = fmt(dBIC))), row.names = FALSE)

if (!is.null(genus_ic)) {
  cat("\n==== D. Genus GLMMs: habitat supported by AICc (dAICc < -2) ====\n")
  cat("Genera with dAICc < -2:", sum(genus_ic$habitat_supported_AICc, na.rm = TRUE),
      "| genera with LR q < 0.05:", sum(genus_ic$omnibus_q < 0.05, na.rm = TRUE), "\n")
  print(as.data.frame(genus_ic %>%
    filter(habitat_supported_AICc | omnibus_q < 0.05 | published_fig4C) %>%
    transmute(genus, order, dAICc_habitat = fmt(dAICc_habitat),
              omnibus_q = fmt(omnibus_q), PvF_q = fmt(paramo_vs_forest_q),
              published_fig4C)), row.names = FALSE)
}

cat("\n==== E. Family screening on rarefied values (AIC comparable across scales) ====\n")
print(as.data.frame(fam_screen %>%
  transmute(metric, family, AICc = fmt(AICc), dAICc = fmt(dAICc))), row.names = FALSE)

cat("\nDone.\n")
