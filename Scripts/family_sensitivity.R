#!/usr/bin/env Rscript
# =============================================================================
# family_sensitivity.R
#
# Are the per-plant habitat conclusions robust to the error family?
#
# Why: the appendix screened gaussian vs lognormal vs gamma by AIC and kept
# lognormal ("best by far"), but AIC was compared across response scales
# without a Jacobian correction, so the log-scale model was not comparable.
# With the correction (model_selection_aicc.R, section E), the preferred family
# on the rarefied values is gaussian for TD q = 0 and PD q = 0, gamma for
# TD q = 1, and lognormal for TD q = 2, PD q = 1 and PD q = 2. The main
# richness result is currently fitted on the log scale.
#
# For each metric and diversity order this refits habitat + site under all
# three families and reports the habitat test (marginal F / LR, adjusted for
# site) and the three Tukey-adjusted habitat contrasts, so the conclusions can
# be compared directly across families.
#
# Run from the repository root, after perplant_hill_size_standardised.R:
#   Rscript Scripts/family_sensitivity.R
#
# Input  : objects/perplant_hill_size_std.rds
# Output : tables/family_sensitivity_habitat.csv
#          tables/family_sensitivity_contrasts.csv
#          tables/family_sensitivity_structure.csv
# =============================================================================

suppressPackageStartupMessages({ library(emmeans); library(dplyr) })

hab_levels <- c("forest", "subparamo", "paramo")
fmt <- function(x) signif(x, 3)
dir.create("tables", showWarnings = FALSE)

vals <- readRDS("objects/perplant_hill_size_std.rds") %>%
  mutate(habitat = factor(as.character(habitat), levels = hab_levels),
         site    = factor(as.character(site)))
grid <- distinct(vals, diversity, q) %>% arrange(desc(diversity), q)

fits_for <- function(d) list(
  gaussian  = lm(estimate ~ habitat + site, data = d),
  lognormal = lm(log(estimate) ~ habitat + site, data = d),
  gamma_log = glm(estimate ~ habitat + site, data = d, family = Gamma(link = "log")))

hab_test <- function(m) {
  if (inherits(m, "glm")) {
    a <- drop1(m, test = "F")           # F test, dispersion estimated
    c(stat = a["habitat", "F value"], p = a["habitat", "Pr(>F)"])
  } else {
    a <- drop1(m, test = "F")
    c(stat = a["habitat", "F value"], p = a["habitat", "Pr(>F)"])
  }
}

hab_rows <- list(); con_rows <- list()
for (i in seq_len(nrow(grid))) {
  d <- filter(vals, diversity == grid$diversity[i], q == grid$q[i])
  lab <- paste0(grid$diversity[i], " q=", grid$q[i])
  for (fam in names(fits_for(d))) {
    m <- fits_for(d)[[fam]]
    ht <- hab_test(m)
    hab_rows[[length(hab_rows) + 1]] <- data.frame(
      metric = lab, family = fam, F = unname(ht["stat"]), p = unname(ht["p"]))

    emm <- emmeans(m, ~ habitat)
    # ratios on the log scales, differences on the identity scale
    pw <- as.data.frame(summary(pairs(emm, adjust = "tukey"),
                                type = if (fam == "gaussian") "link" else "response",
                                infer = c(TRUE, TRUE)))
    est_col <- intersect(c("ratio", "estimate", "odds.ratio"), names(pw))[1]
    con_rows[[length(con_rows) + 1]] <- data.frame(
      metric = lab, family = fam, contrast = as.character(pw$contrast),
      scale = if (fam == "gaussian") "difference" else "ratio",
      estimate = pw[[est_col]], p = pw$p.value)
  }
}

hab_tab <- bind_rows(hab_rows)
con_tab <- bind_rows(con_rows)
write.csv(hab_tab, "tables/family_sensitivity_habitat.csv", row.names = FALSE)
write.csv(con_tab, "tables/family_sensitivity_contrasts.csv", row.names = FALSE)

cat("\n==== Habitat effect (adjusted for site) by error family ====\n")
print(as.data.frame(hab_tab %>% mutate(F = fmt(F), p = fmt(p)) %>%
  tidyr::pivot_wider(names_from = family, values_from = c(F, p))), row.names = FALSE)

cat("\n==== Tukey habitat contrasts by error family (p-values) ====\n")
print(as.data.frame(con_tab %>% select(metric, family, contrast, p) %>%
  mutate(p = fmt(p)) %>%
  tidyr::pivot_wider(names_from = family, values_from = p)), row.names = FALSE)

# ---- Model structure under each family --------------------------------------
# The structure comparison in model_selection_aicc.R was run on the log-normal
# models. This repeats it within each family, so family and structure are not
# chosen from different fits.
ic <- function(m, n) {
  ll <- logLik(m); k <- attr(ll, "df")
  aic <- -2 * as.numeric(ll) + 2 * k
  aic + 2 * k * (k + 1) / (n - k - 1)
}
struct <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  d <- filter(vals, diversity == grid$diversity[i], q == grid$q[i])
  lab <- paste0(grid$diversity[i], " q=", grid$q[i])
  jac <- sum(log(d$estimate))          # log-scale models put on the response scale
  bind_rows(lapply(c("gaussian", "lognormal", "gamma_log"), function(fam) {
    f <- function(rhs) {
      if (fam == "lognormal") as.formula(paste("log(estimate) ~", rhs))
      else as.formula(paste("estimate ~", rhs))
    }
    fit <- function(rhs) {
      if (fam == "gamma_log") glm(f(rhs), data = d, family = Gamma(link = "log"))
      else lm(f(rhs), data = d)
    }
    adj <- if (fam == "lognormal") 2 * jac else 0
    v <- vapply(c("site", "habitat + site", "habitat * site"),
                function(rhs) ic(fit(rhs), nrow(d)) + adj, numeric(1))
    data.frame(metric = lab, family = fam,
               AICc_site = v[1], AICc_additive = v[2], AICc_interaction = v[3],
               dAICc_habitat = v[2] - v[1], dAICc_interaction = v[3] - v[2],
               best = c("site only", "habitat + site", "habitat * site")[which.min(v)])
  }))
}))
write.csv(struct, "tables/family_sensitivity_structure.csv", row.names = FALSE)

cat("\n==== Model structure within each family (AICc) ====\n")
print(as.data.frame(struct %>%
  transmute(metric, family, dAICc_habitat = fmt(dAICc_habitat),
            dAICc_interaction = fmt(dAICc_interaction), best)), row.names = FALSE)

cat("\n==== Total AICc across all six metric/order combinations, by family ====\n")
print(as.data.frame(struct %>% group_by(family) %>%
  summarise(total_dAICc_vs_best = fmt(sum(AICc_additive)), .groups = "drop") %>%
  arrange(total_dAICc_vs_best)), row.names = FALSE)

cat("\n==== Tukey habitat contrasts: effect sizes ====\n")
print(as.data.frame(con_tab %>% mutate(estimate = fmt(estimate)) %>%
  select(metric, contrast, family, scale, estimate)), row.names = FALSE)

cat("\nDone.\n")
