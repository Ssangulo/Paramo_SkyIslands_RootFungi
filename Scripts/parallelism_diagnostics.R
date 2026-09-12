#!/usr/bin/env Rscript
# =============================================================================
# parallelism_diagnostics.R
#
# Two checks before finalising the "parallel restructuring" interpretation.
#
# A. Which site(s) drive the significant habitat x site interaction for
#    abundance-weighted per-plant PD (q = 1, 2)?
#      - site-specific habitat trajectories and within-site contrasts
#      - leave-one-site-out interaction tests (+ NV & DOM only, the two
#        complete gradients, as for the PERMANOVA)
#      - influential plants (Cook's distance) and whether the interaction
#        survives without them
#
# B. Does it sit uneasily with the dominant-lineage result?
#      - The published beta regressions (Tables S16-S17) fit
#        Proportion ~ elevation + site, which ASSUMES equal elevation slopes
#        across sites; parallelism of lineage turnover was never tested.
#        Here elevation x site and habitat x site interactions are tested
#        (likelihood-ratio tests) and site-specific slopes are reported.
#      - Per-plant PD (q = 1, 2) is related to Helotiales/Sebacinales
#        dominance, and the PD interaction is re-tested with lineage
#        dominance as a covariate.
#
# Run from the repository root, after perplant_hill_size_standardised.R:
#   Rscript Scripts/parallelism_diagnostics.R
#
# Inputs : objects/perplant_hill_size_std.rds, objects/ps_individual.rds
# Outputs: tables/diag_*.csv, figures/diag_parallelism_sites.png
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(emmeans)
  library(betareg)
  library(ggplot2)
})

out_tab_dir <- "tables"
out_fig_dir <- "figures"
out_obj_dir <- "objects"
dir.create(out_tab_dir, showWarnings = FALSE)
dir.create(out_fig_dir, showWarnings = FALSE)

hab_levels  <- c("forest", "subparamo", "paramo")
labs_map    <- c(forest = "Forest", subparamo = "Subpáramo", paramo = "Páramo")
site_colors <- c(BEL = "#E69F00", DOM = "#009E73", MA = "#56B4E9", NV = "#D55E00")
site_labels <- c(BEL = "Belmira", DOM = "Las Domínguez", MA = "Matarredonda", NV = "La Nevera")

# Error family for the per-plant models; must match perplant_hill_size_standardised.R
FAMILY <- "gamma_log"          # "gamma_log" or "lognormal"

fit_model <- function(rhs, d) {
  if (FAMILY == "gamma_log") {
    glm(as.formula(paste("estimate ~", rhs)), data = d, family = Gamma(link = "log"))
  } else {
    lm(as.formula(paste("log(estimate) ~", rhs)), data = d)
  }
}

fmt <- function(x) signif(x, 3)
show <- function(title, df) {
  cat("\n==== ", title, " ====\n", sep = "")
  print(as.data.frame(df), row.names = FALSE)
}

# ---- Data -------------------------------------------------------------------
vals <- readRDS(file.path(out_obj_dir, "perplant_hill_size_std.rds")) %>%
  mutate(habitat = factor(as.character(habitat), levels = hab_levels),
         site    = factor(as.character(site)))
sites <- levels(vals$site)

ps_individual <- readRDS(file.path(out_obj_dir, "ps_individual.rds"))

# F-test of habitat x site: additive vs interaction model (optional covariates)
int_test <- function(d, extra = NULL) {
  rhs_add <- "habitat + site"
  rhs_int <- "habitat * site"
  if (!is.null(extra)) {
    rhs_add <- paste(rhs_add, "+", extra)
    rhs_int <- paste(rhs_int, "+", extra)
  }
  m_add <- fit_model(rhs_add, d)
  m_int <- fit_model(rhs_int, d)
  a <- if (FAMILY == "gamma_log") anova(m_add, m_int, test = "F") else anova(m_add, m_int)
  if (nrow(a) < 2 || is.na(a$Df[2]) || a$Df[2] == 0) {
    return(data.frame(n_plants = nrow(d), df1 = 0, df2 = NA, F = NA, p = NA))
  }
  a_rdf <- grep("^Res", names(a), value = TRUE)[1]   # "Res.Df" (lm) / "Resid. Df" (glm)
  a_p   <- grep("^Pr",  names(a), value = TRUE)[1]
  data.frame(n_plants = nrow(d), df1 = abs(a$Df[2]), df2 = a[[a_rdf]][2], F = a$F[2],
             p = a[[a_p]][2])
}

grid <- distinct(vals, diversity, q) %>% arrange(desc(diversity), q)

# =============================================================================
# A. Which site drives the interaction?
# =============================================================================

# ---- A1. Leave-one-site-out interaction tests (all metrics, for comparison) --
loso <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  d <- filter(vals, diversity == grid$diversity[i], q == grid$q[i])
  rows <- list(cbind(subset = "all sites", int_test(d)))
  for (s in sites) {
    rows[[length(rows) + 1]] <- cbind(subset = paste("without", s),
                                      int_test(droplevels(filter(d, site != s))))
  }
  if (all(c("NV", "DOM") %in% sites)) {
    rows[[length(rows) + 1]] <- cbind(subset = "NV + DOM only (complete gradients)",
                                      int_test(droplevels(filter(d, site %in% c("NV", "DOM")))))
  }
  cbind(diversity = grid$diversity[i], q = grid$q[i], bind_rows(rows))
}))
write.csv(loso, file.path(out_tab_dir, "diag_leave_one_site_out.csv"), row.names = FALSE)

# ---- A2. Site-specific trajectories and within-site contrasts ----------------
traj <- list()
contr <- list()
for (i in seq_len(nrow(grid))) {
  d <- filter(vals, diversity == grid$diversity[i], q == grid$q[i])
  m <- fit_model("habitat * site", d)
  emm <- emmeans(m, ~ habitat | site)
  traj[[i]] <- as.data.frame(summary(emm, type = "response")) %>%
    filter(!is.na(response)) %>%
    mutate(diversity = grid$diversity[i], q = grid$q[i], .before = 1)
  contr[[i]] <- as.data.frame(summary(pairs(emm), type = "response", infer = c(TRUE, TRUE))) %>%
    filter(!is.na(ratio)) %>%
    mutate(diversity = grid$diversity[i], q = grid$q[i], .before = 1)
}
traj  <- bind_rows(traj)
contr <- bind_rows(contr)
write.csv(traj,  file.path(out_tab_dir, "diag_site_trajectories.csv"),   row.names = FALSE)
write.csv(contr, file.path(out_tab_dir, "diag_within_site_contrasts.csv"), row.names = FALSE)

# ---- A3. Influential plants (PD q = 1, 2) ------------------------------------
infl <- bind_rows(lapply(c(1, 2), function(qq) {
  d <- filter(vals, diversity == "PD", q == qq)
  m <- fit_model("habitat * site", d)
  cd <- cooks.distance(m)
  flagged <- d$plant[is.finite(cd) & cd > 4 / nrow(d)]
  without <- int_test(droplevels(filter(d, !plant %in% flagged)))
  data.frame(diversity = "PD", q = qq,
             plant = d$plant, site = d$site, habitat = d$habitat,
             estimate = d$estimate, cooks_d = cd) %>%
    filter(plant %in% flagged) %>%
    arrange(desc(cooks_d)) %>%
    mutate(interaction_p_without_flagged = without$p,
           interaction_F_without_flagged = without$F)
}))
write.csv(infl, file.path(out_tab_dir, "diag_influential_plants_PD.csv"), row.names = FALSE)

# =============================================================================
# B. Dominant lineages
# =============================================================================

# ---- B1. Per-plant lineage proportions (same preparation as Tables S16-S17) --
otu <- as(otu_table(ps_individual), "matrix")
if (!taxa_are_rows(ps_individual)) otu <- t(otu)
tax <- as.data.frame(tax_table(ps_individual), stringsAsFactors = FALSE)
stopifnot("Order" %in% names(tax))
ord <- tolower(as.character(tax[rownames(otu), "Order"]))
ord[is.na(ord)] <- ""

sd <- as(sample_data(ps_individual), "data.frame")[colnames(otu), , drop = FALSE]
if (!"elevation" %in% names(sd)) stop("No 'elevation' column in ps_individual sample data.")

total <- colSums(otu)
lin <- data.frame(plant     = colnames(otu),
                  site      = factor(as.character(sd$site)),
                  habitat   = factor(as.character(sd$habitat), levels = hab_levels),
                  elevation = as.numeric(sd$elevation),
                  stringsAsFactors = FALSE)
lin$elev_100 <- lin$elevation / 100
lin$cell     <- interaction(lin$site, lin$habitat, drop = TRUE)  # observed site x habitat cells

n_pl <- nrow(lin)
for (o in c("helotiales", "sebacinales")) {
  p <- colSums(otu[grepl(o, ord, fixed = TRUE), , drop = FALSE]) / total
  lin[[paste0(o, "_prop")]] <- p
  lin[[paste0(o, "_sv")]]   <- (p * (n_pl - 1) + 0.5) / n_pl  # Smithson-Verkuilen, as published
}

# ---- B2. Are the lineage trends parallel across sites? -----------------------
lr_test <- function(m0, m1) {
  ll0 <- logLik(m0); ll1 <- logLik(m1)
  stat <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  df   <- attr(ll1, "df") - attr(ll0, "df")
  data.frame(chisq = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
}

lineage_tests <- list()
lineage_slopes <- list()
for (o in c("helotiales", "sebacinales")) {
  y <- paste0(o, "_sv")

  # Published structure (raw elevation) — compare with Tables S16-S17
  m_pub <- betareg(as.formula(paste(y, "~ elevation + site")), data = lin)
  pub_coef <- summary(m_pub)$coefficients$mean["elevation", ]

  # Elevation slopes: common vs site-specific
  m_el_add <- betareg(as.formula(paste(y, "~ elev_100 + site")), data = lin)
  m_el_int <- betareg(as.formula(paste(y, "~ elev_100 * site")), data = lin)

  # Habitat effects: additive vs full site x habitat (cell-means; equivalent to
  # habitat * site with empty cells dropped)
  m_h_add <- betareg(as.formula(paste(y, "~ habitat + site")), data = lin)
  m_h_int <- betareg(as.formula(paste(y, "~ cell")), data = lin)

  lineage_tests[[o]] <- bind_rows(
    cbind(order = o, test = "elevation x site (LR)", lr_test(m_el_add, m_el_int)),
    cbind(order = o, test = "habitat x site (LR)",   lr_test(m_h_add,  m_h_int))
  ) %>%
    mutate(published_elev_coef = pub_coef["Estimate"],
           published_elev_p    = pub_coef["Pr(>|z|)"])

  lineage_slopes[[o]] <- tryCatch(
    as.data.frame(summary(emtrends(m_el_int, ~ site, var = "elev_100"),
                          infer = c(TRUE, TRUE))) %>%
      mutate(order = o, .before = 1),
    error = function(e) {
      message("emtrends failed for ", o, ": ", conditionMessage(e))
      NULL
    })
}
lineage_tests  <- bind_rows(lineage_tests)
lineage_slopes <- bind_rows(lineage_slopes)
write.csv(lineage_tests,  file.path(out_tab_dir, "diag_lineage_parallelism_tests.csv"), row.names = FALSE)
write.csv(lineage_slopes, file.path(out_tab_dir, "diag_lineage_site_slopes.csv"),       row.names = FALSE)

# ---- B3. Does lineage dominance account for the PD interaction? --------------
pd_lin <- vals %>%
  filter(diversity == "PD", q %in% c(1, 2)) %>%
  left_join(select(lin, plant, helotiales_sv, sebacinales_sv), by = "plant") %>%
  mutate(hel_logit = qlogis(helotiales_sv),
         seb_logit = qlogis(sebacinales_sv))
stopifnot(!anyNA(pd_lin$hel_logit))

pd_vs_lineage <- bind_rows(lapply(c(1, 2), function(qq) {
  d <- filter(pd_lin, q == qq)
  assoc <- function(v) {
    m  <- fit_model(paste(v, "+ habitat + site"), d)
    co <- summary(m)$coefficients
    pc <- grep("^Pr", colnames(co), value = TRUE)[1]
    c(estimate = co[v, "Estimate"], p = co[v, pc])
  }
  a_h <- assoc("hel_logit")
  a_s <- assoc("seb_logit")
  bind_rows(
    cbind(covariate = "none",                    int_test(d)),
    cbind(covariate = "Helotiales (logit)",      int_test(d, "hel_logit")),
    cbind(covariate = "Sebacinales (logit)",     int_test(d, "seb_logit")),
    cbind(covariate = "Helotiales + Sebacinales", int_test(d, "hel_logit + seb_logit"))
  ) %>%
    mutate(diversity = "PD", q = qq, .before = 1,
           hel_assoc_slope = a_h[["estimate"]], hel_assoc_p = a_h[["p"]],
           seb_assoc_slope = a_s[["estimate"]], seb_assoc_p = a_s[["p"]])
}))
write.csv(pd_vs_lineage, file.path(out_tab_dir, "diag_PD_interaction_vs_lineages.csv"), row.names = FALSE)

# ---- Diagnostic figure: all four sites, including MA -------------------------
plot_df <- bind_rows(
  vals %>% filter(diversity == "PD", q %in% c(1, 2)) %>%
    transmute(plant, site, habitat, panel = paste0("PD (q=", q, ")"), value = estimate),
  lin %>% transmute(plant, site, habitat, panel = "Helotiales (proportion)",  value = helotiales_prop),
  lin %>% transmute(plant, site, habitat, panel = "Sebacinales (proportion)", value = sebacinales_prop)
) %>%
  mutate(habitat_lab = factor(unname(labs_map[as.character(habitat)]), levels = unname(labs_map)),
         site  = factor(as.character(site), levels = names(site_colors)),
         panel = factor(panel, levels = c("PD (q=1)", "PD (q=2)",
                                          "Helotiales (proportion)", "Sebacinales (proportion)")))

dodge <- position_dodge(width = 0.30)
p_diag <- ggplot(plot_df, aes(x = habitat_lab, y = value, colour = site, group = site)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.30),
             alpha = 0.4, size = 1) +
  stat_summary(fun = mean, geom = "line",  position = dodge, linewidth = 0.8) +
  stat_summary(fun = mean, geom = "point", position = dodge, size = 2) +
  facet_wrap(~ panel, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = site_colors, labels = site_labels, name = "Site") +
  labs(x = "Habitat", y = NULL,
       title = "Diagnostic: site trajectories (all sites; lines join site means)") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out_fig_dir, "diag_parallelism_sites.png"),
       plot = p_diag, width = 18, height = 14, units = "cm", dpi = 300, bg = "white")

# ---- Console summary --------------------------------------------------------
show("A1. Habitat x site interaction, leaving out one site at a time",
     loso %>% mutate(F = fmt(F), p = fmt(p)))

show("A2. Site-specific means (PD q = 1, 2; back-transformed)",
     traj %>% filter(diversity == "PD", q %in% c(1, 2)) %>%
       transmute(q, site, habitat, response = fmt(response),
                 lower = fmt(lower.CL), upper = fmt(upper.CL)))

show("A2. Within-site habitat contrasts (PD q = 1, 2; ratios)",
     contr %>% filter(diversity == "PD", q %in% c(1, 2)) %>%
       transmute(q, site, contrast, ratio = fmt(ratio), p = fmt(p.value)))

show("A3. Influential plants (Cook's D > 4/n) and interaction without them",
     infl %>% mutate(across(c(estimate, cooks_d, interaction_p_without_flagged,
                              interaction_F_without_flagged), fmt)))

show("B2. Parallelism of lineage trends (likelihood-ratio tests)",
     lineage_tests %>% mutate(across(c(chisq, p, published_elev_coef, published_elev_p), fmt)))

if (nrow(lineage_slopes) > 0) {
  slope_col <- intersect(c("elev_100.trend", "elev_100"), names(lineage_slopes))[1]
  show("B2. Site-specific elevation slopes (proportion scale, per 100 m)",
       lineage_slopes %>%
         transmute(order, site, slope = fmt(.data[[slope_col]]),
                   lower = fmt(asymp.LCL), upper = fmt(asymp.UCL), p = fmt(p.value)))
}

show("B3. PD interaction with lineage dominance as covariate",
     pd_vs_lineage %>% mutate(across(c(F, p, hel_assoc_slope, hel_assoc_p,
                                       seb_assoc_slope, seb_assoc_p), fmt)))

cat("\nFigure: ", file.path(out_fig_dir, "diag_parallelism_sites.png"), "\nDone.\n", sep = "")
