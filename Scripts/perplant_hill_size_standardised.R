#!/usr/bin/env Rscript
# =============================================================================
# perplant_hill_size_standardised.R
#
# Per-plant Hill diversity (TD and meanPD; q = 0, 1, 2) standardised to a
# common sequencing effort with iNEXT.3D (size-based rarefaction), followed by
# omnibus and post-hoc habitat tests taken from the SAME model.
#
# Replaces the Section 5.3 inputs (hillR on unrarefied summed counts) and fixes
# the mismatch between tests:
#   - Tables S14-S15 used anova() on log(y) ~ habitat * site. anova() on an lm
#     is sequential (Type I): habitat was tested first, NOT adjusted for site,
#     in a design where habitat and site are partly confounded (MA has no
#     forest, BEL no paramo).
#   - Table S15b Tukey contrasts came from log(y) ~ habitat + site, i.e.
#     habitat adjusted for site.
#   The omnibus and pairwise tests therefore asked different questions.
#
# Here, for each metric:
#   1. Parallelism: habitat x site interaction via anova(additive, interaction).
#   2. Omnibus habitat effect, adjusted for site: drop1() on the additive model.
#   3. Tukey pairwise contrasts (ratios) from the same additive model.
#   4. Polynomial (linear/quadratic) contrasts on ordered habitat: direction/shape.
# A diagnostic section re-fits the published unstandardised hillR q = 1 values
# with both Type I and site-adjusted tests to show the source of the mismatch.
#
# Run from the repository root:
#   Rscript Scripts/perplant_hill_size_standardised.R          # depth = min library
#   Rscript Scripts/perplant_hill_size_standardised.R 25000    # custom depth
#
# Inputs : objects/ps_individual.rds (per plant; root samples summed)
#          objects/tree_ps.rds       (ML tree)
# Outputs: objects/perplant_hill_size_std.rds
#          tables/perplant_hill_size_std_values.csv
#          tables/perplant_hill_size_std_omnibus.csv
#          tables/perplant_hill_size_std_pairwise.csv
#          tables/perplant_hill_size_std_trend.csv
#          tables/perplant_hill_size_std_emmeans.csv
#          tables/perplant_hill_published_typeI_vs_adjusted.csv (diagnostic)
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(iNEXT.3D)
  library(ape)
  library(phangorn)
  library(dplyr)
  library(emmeans)
})

# Error family for the per-plant models. The appendix screened gaussian,
# lognormal and gamma by AIC and kept lognormal, but that comparison was made
# across response scales without a Jacobian correction. Corrected
# (model_selection_aicc.R section E; family_sensitivity.R), gamma with a log
# link has the lowest summed AICc across the six metric/order combinations and
# is preferred or tied for four of them. It also gives ratios, as reported.
# Set to "lognormal" to reproduce the previous version of this script.
FAMILY <- "gamma_log"          # "gamma_log" or "lognormal"

args        <- commandArgs(trailingOnly = TRUE)
out_tab_dir <- "tables"
out_obj_dir <- "objects"
dir.create(out_tab_dir, showWarnings = FALSE)

hab_levels <- c("forest", "subparamo", "paramo")
labs_map   <- c(forest = "Forest", subparamo = "Subpáramo", paramo = "Páramo")

message("iNEXT.3D version: ", as.character(packageVersion("iNEXT.3D")))

# ---- Load objects, shared tree ----------------------------------------------
ps_individual <- readRDS(file.path(out_obj_dir, "ps_individual.rds"))
tree_ps       <- readRDS(file.path(out_obj_dir, "tree_ps.rds"))

tr_full     <- phy_tree(tree_ps)
common_taxa <- intersect(taxa_names(ps_individual), tr_full$tip.label)
tr <- ape::keep.tip(tr_full, common_taxa)
tr$node.label <- NULL
if (!ape::is.rooted(tr)) tr <- phangorn::midpoint(tr)
reftime <- max(ape::node.depth.edgelength(tr))
message("Taxa in tree: ", length(common_taxa), " | PD reference time: ", signif(reftime, 5))

# ---- Community matrix (taxa x plants, summed read counts) --------------------
ps   <- prune_taxa(common_taxa, ps_individual)
comm <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) comm <- t(comm)
comm <- comm[tr$tip.label, , drop = FALSE]
storage.mode(comm) <- "numeric"

meta <- as(sample_data(ps), "data.frame")[colnames(comm), , drop = FALSE]
meta$plant   <- colnames(comm)
meta$habitat <- factor(as.character(meta$habitat), levels = hab_levels)
meta$site    <- factor(as.character(meta$site))
stopifnot(!anyNA(meta$habitat), !anyNA(meta$site))

lib   <- colSums(comm)
DEPTH <- if (length(args) >= 1) as.numeric(args[1]) else min(lib)
if (DEPTH > min(lib)) warning("DEPTH exceeds the smallest library; some plants will be extrapolated.")

message("Plants: ", ncol(comm), " (", paste(names(table(meta$habitat)), table(meta$habitat),
                                          sep = "=", collapse = ", "), ")")
message("Library size range: ", min(lib), "-", max(lib), " | standardisation depth: ", DEPTH)
message("Singletons per plant (median [range]): ", median(colSums(comm == 1)),
        " [", min(colSums(comm == 1)), "-", max(colSums(comm == 1)), "]")

# ---- Size-standardised Hill numbers per plant -------------------------------
message("Estimating size-standardised TD ...")
est_td <- estimate3D(data = comm, diversity = "TD", q = c(0, 1, 2),
                     datatype = "abundance", base = "size", level = DEPTH, nboot = 0)
message("Estimating size-standardised PD (meanPD) ...")
est_pd <- estimate3D(data = comm, diversity = "PD", q = c(0, 1, 2),
                     datatype = "abundance", base = "size", level = DEPTH, nboot = 0,
                     PDtree = tr, PDreftime = reftime, PDtype = "meanPD")

pick_col <- function(df, exact, what) {
  hit <- exact[exact %in% names(df)]
  if (length(hit) == 0) stop("Cannot find ", what, " column. Columns: ",
                             paste(names(df), collapse = ", "))
  hit[1]
}

tidy_est <- function(df, div) {
  df <- as.data.frame(df)
  asm <- pick_col(df, c("Assemblage", "Community"), "assemblage")
  ord <- pick_col(df, c("Order.q", "Order"), "order")
  est <- pick_col(df, c(paste0("q", div), "qD"), "estimate")
  data.frame(plant     = as.character(df[[asm]]),
             diversity = div,
             q         = as.numeric(df[[ord]]),
             estimate  = df[[est]],
             coverage  = if ("SC" %in% names(df)) df[["SC"]] else NA_real_,
             method    = if ("Method" %in% names(df)) as.character(df[["Method"]]) else NA_character_,
             stringsAsFactors = FALSE)
}

vals <- bind_rows(tidy_est(est_td, "TD"), tidy_est(est_pd, "PD")) %>%
  left_join(meta[, c("plant", "habitat", "site")], by = "plant")
stopifnot(!anyNA(vals$habitat), nrow(vals) == ncol(comm) * 2 * 3)

message("Error family for the habitat models: ", FAMILY)
message("Methods used: ", paste(names(table(vals$method)), table(vals$method),
                                sep = "=", collapse = ", "))
message("Sample coverage at standardisation depth (TD): ",
        paste(signif(range(vals$coverage[vals$diversity == "TD"], na.rm = TRUE), 4),
              collapse = "-"))

saveRDS(vals, file.path(out_obj_dir, "perplant_hill_size_std.rds"))
write.csv(vals, file.path(out_tab_dir, "perplant_hill_size_std_values.csv"), row.names = FALSE)

# ---- Models: omnibus and post-hoc from the same (additive) model -------------
fit_model <- function(rhs, d) {
  if (FAMILY == "gamma_log") {
    glm(as.formula(paste("estimate ~", rhs)), data = d, family = Gamma(link = "log"))
  } else {
    lm(as.formula(paste("log(estimate) ~", rhs)), data = d)
  }
}

fit_one <- function(d) {
  d     <- droplevels(d)
  m_int <- fit_model("habitat * site", d)
  m_add <- fit_model("habitat + site", d)

  a_int <- if (FAMILY == "gamma_log") anova(m_add, m_int, test = "F") else anova(m_add, m_int)
  a_rdf <- grep("^Res", names(a_int), value = TRUE)[1]   # "Res.Df" (lm) / "Resid. Df" (glm)
  a_p   <- grep("^Pr",  names(a_int), value = TRUE)[1]
  d1    <- drop1(m_add, test = "F")
  df2   <- df.residual(m_add)

  omnibus <- data.frame(
    term = c("habitat x site (parallelism)", "habitat (adjusted for site)", "site (adjusted for habitat)"),
    df1  = c(abs(a_int$Df[2]), d1["habitat", "Df"], d1["site", "Df"]),
    df2  = c(a_int[[a_rdf]][2], df2, df2),
    F    = c(a_int$F[2], d1["habitat", "F value"], d1["site", "F value"]),
    p    = c(a_int[[a_p]][2], d1["habitat", "Pr(>F)"], d1["site", "Pr(>F)"]))

  emm <- emmeans(m_add, ~ habitat)
  pw  <- as.data.frame(summary(pairs(emm, adjust = "tukey"), type = "response",
                               infer = c(TRUE, TRUE)))
  trd <- as.data.frame(summary(contrast(emm, "poly"), infer = c(TRUE, TRUE)))
  mm  <- as.data.frame(summary(emm, type = "response"))

  list(omnibus = omnibus, pairwise = pw, trend = trd, emmeans = mm)
}

grid <- unique(vals[, c("diversity", "q")])
res  <- lapply(seq_len(nrow(grid)), function(i) {
  d <- vals[vals$diversity == grid$diversity[i] & vals$q == grid$q[i], ]
  r <- fit_one(d)
  lapply(r, function(x) cbind(diversity = grid$diversity[i], q = grid$q[i], x))
})

omnibus  <- bind_rows(lapply(res, `[[`, "omnibus"))
pairwise <- bind_rows(lapply(res, `[[`, "pairwise"))
trend    <- bind_rows(lapply(res, `[[`, "trend"))
emm_tab  <- bind_rows(lapply(res, `[[`, "emmeans"))

write.csv(omnibus,  file.path(out_tab_dir, "perplant_hill_size_std_omnibus.csv"),  row.names = FALSE)
write.csv(pairwise, file.path(out_tab_dir, "perplant_hill_size_std_pairwise.csv"), row.names = FALSE)
write.csv(trend,    file.path(out_tab_dir, "perplant_hill_size_std_trend.csv"),    row.names = FALSE)
write.csv(emm_tab,  file.path(out_tab_dir, "perplant_hill_size_std_emmeans.csv"),  row.names = FALSE)

# ---- Diagnostic: published unstandardised q = 1 values -----------------------
# Reproduces Tables S14-S15 (Type I, interaction model) and compares with the
# site-adjusted habitat test from the additive model.
diag_tab <- NULL
if (requireNamespace("hillR", quietly = TRUE)) {
  get_hill_phylo <- function(x, tree, qval = 1) {
    out <- tryCatch(hillR::hill_phylo(x, tree = tree, q = qval), error = function(e) NULL)
    if (is.null(out)) out <- hillR::hill_phylo(x, tree, qvalue = qval)
    out
  }
  get_hill_taxa <- function(x, qval = 1) {
    out <- tryCatch(hillR::hill_taxa(x, q = qval), error = function(e) NULL)
    if (is.null(out)) out <- hillR::hill_taxa(x, qvalue = qval)
    out
  }
  cm  <- t(comm)
  raw <- list(TD = get_hill_taxa(cm, 1), PD = get_hill_phylo(cm, tr, 1))
  diag_tab <- bind_rows(lapply(names(raw), function(nm) {
    d <- meta
    d$estimate <- as.numeric(raw[[nm]][d$plant])
    tI <- anova(lm(log(estimate) ~ habitat * site, data = d))
    d1 <- drop1(lm(log(estimate) ~ habitat + site, data = d), test = "F")
    data.frame(diversity = nm, q = 1,
               test = c("Type I habitat, interaction model (published S14-S15)",
                        "habitat adjusted for site, additive model"),
               df1 = c(tI["habitat", "Df"], d1["habitat", "Df"]),
               df2 = c(tI["Residuals", "Df"], df.residual(lm(log(estimate) ~ habitat + site, data = d))),
               F   = c(tI["habitat", "F value"], d1["habitat", "F value"]),
               p   = c(tI["habitat", "Pr(>F)"], d1["habitat", "Pr(>F)"]))
  }))
  write.csv(diag_tab, file.path(out_tab_dir, "perplant_hill_published_typeI_vs_adjusted.csv"),
            row.names = FALSE)
} else {
  message("hillR not installed; published-model diagnostic skipped.")
}

# ---- Console summary --------------------------------------------------------
fmt <- function(x) signif(x, 3)

cat("\n==== Omnibus tests (size-standardised; habitat adjusted for site) ====\n")
print(omnibus %>% mutate(F = fmt(F), p = fmt(p)), row.names = FALSE)

cat("\n==== Tukey pairwise contrasts (ratio of effective diversity) ====\n")
ci_lo <- intersect(c("lower.CL", "asymp.LCL"), names(pairwise))[1]
ci_hi <- intersect(c("upper.CL", "asymp.UCL"), names(pairwise))[1]
print(pairwise %>%
        transmute(diversity, q, contrast, ratio = fmt(ratio),
                  lower = fmt(.data[[ci_lo]]), upper = fmt(.data[[ci_hi]]),
                  p_tukey = fmt(p.value)),
      row.names = FALSE)

cat("\n==== Polynomial contrasts on ordered habitat (log scale) ====\n")
print(trend %>% transmute(diversity, q, contrast, estimate = fmt(estimate), p = fmt(p.value)),
      row.names = FALSE)

if (!is.null(diag_tab)) {
  cat("\n==== Diagnostic: published unstandardised q = 1 values ====\n")
  print(diag_tab %>% mutate(F = fmt(F), p = fmt(p)), row.names = FALSE)
}

cat("\nDone.\n")

# =============================================================================
# Figure 3B (corrected): per-plant, size-standardised Hill diversity by site
# Same aesthetic as the published Figure 3B: jittered per-plant values with
# site-specific model means (+/- 95% CI) from log(y) ~ habitat * site; MA
# excluded from the plot (no forest baseline). Uses `vals` from above if in
# memory; otherwise loads objects/perplant_hill_size_std.rds, so this block
# can be run on its own.
# If objects/Figure_3A_panel_plant_units.rds exists (from iNEXT3D_plant_level.R),
# the full corrected Figure 3 is also assembled.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(emmeans)
})

Q_PLOT <- c(0, 1, 2)   # diversity orders to show; c(1) reproduces the published layout
FIG_W  <- 30           # full Figure 3 width (cm); published was 28 with q = 1 only
FIG_H  <- 12           # full Figure 3 height (cm)

out_fig_dir <- "figures"
dir.create(out_fig_dir, showWarnings = FALSE)

if (!exists("vals")) vals <- readRDS(file.path("objects", "perplant_hill_size_std.rds"))
depth_lab <- if (exists("DEPTH")) format(DEPTH, big.mark = ",") else "common depth"

# Palette and theme copied from the published Figure 3 code
site_colors <- c(BEL = "#E69F00", DOM = "#009E73", MA = "#56B4E9", NV = "#D55E00")
site_labels <- c(BEL = "Belmira", DOM = "Las Domínguez", MA = "Matarredonda", NV = "La Nevera")
hab_lab     <- c(forest = "Forest", subparamo = "Subpáramo", paramo = "Páramo")

base_theme <- theme_classic(base_size = 11) +
  theme(
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 7, color = "black"),
    legend.title     = element_text(size = 10.2, face = "bold"),
    legend.text      = element_text(size = 10.2),
    legend.key.size  = unit(0.89, "lines"),
    legend.position  = "right",
    plot.tag         = element_text(size = 13, face = "bold"),
    strip.background = element_blank()
  )

metric_levels <- c(paste0("TD (q=", Q_PLOT, ")"), paste0("PD (q=", Q_PLOT, ")"))

vals_plot <- vals %>%
  filter(q %in% Q_PLOT) %>%
  mutate(
    habitat = factor(as.character(habitat), levels = names(hab_lab)),
    site    = factor(as.character(site), levels = names(site_colors)),
    metric  = factor(paste0(diversity, " (q=", q, ")"), levels = metric_levels)
  )

# Site-specific means from the interaction model (as in the published panel);
# non-estimable habitat x site cells (MA forest, BEL paramo) are dropped.
emm_plot <- bind_rows(lapply(split(vals_plot, vals_plot$metric, drop = TRUE), function(d) {
  m <- fit_model("habitat * site", d)
  as.data.frame(summary(emmeans(m, ~ habitat | site), infer = c(TRUE, TRUE))) %>%
    filter(!is.na(emmean)) %>%
    mutate(response = exp(emmean),
           lower.CL = exp(lower.CL),
           upper.CL = exp(upper.CL),
           metric   = as.character(unique(d$metric)))
})) %>%
  mutate(metric = factor(metric, levels = metric_levels))

relabel <- function(df) {
  df %>% mutate(
    habitat_lab = factor(unname(hab_lab[as.character(habitat)]), levels = unname(hab_lab)),
    site_fac    = factor(as.character(site), levels = names(site_colors))
  )
}

vals_B <- relabel(vals_plot) %>% filter(site != "MA")
emm_B  <- relabel(emm_plot)  %>% filter(site != "MA")

model_dodge <- position_dodge(width = 0.30)

panel_B <- ggplot(vals_B, aes(x = habitat_lab, y = estimate, color = site_fac)) +
  geom_point(aes(group = site_fac),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.30),
             alpha = 0.45, size = 1.2) +
  geom_errorbar(data = emm_B,
                aes(x = habitat_lab, ymin = lower.CL, ymax = upper.CL,
                    color = site_fac, group = site_fac),
                width = 0.10, linewidth = 0.5, position = model_dodge, inherit.aes = FALSE) +
  geom_line(data = emm_B,
            aes(x = habitat_lab, y = response, color = site_fac, group = site_fac),
            linewidth = 0.9, position = model_dodge, inherit.aes = FALSE) +
  geom_point(data = emm_B,
             aes(x = habitat_lab, y = response, color = site_fac, group = site_fac),
             size = 2.2, position = model_dodge, inherit.aes = FALSE) +
  facet_wrap(~ metric, scales = "free_y", nrow = if (length(Q_PLOT) == 1) 1 else 2) +
  scale_color_manual(values = site_colors, labels = site_labels, name = "Site") +
  labs(x = "Habitat",
       y = paste0("Per-plant Hill diversity (", depth_lab, " reads)"),
       tag = "B") +
  base_theme +
  theme(legend.position = "right",
        axis.title  = element_text(size = 9),
        axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1),
        strip.text  = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 2, linewidth = 0.8)))

saveRDS(panel_B, file.path("objects", "Figure_3B_panel_size_std.rds"))
ggsave(file.path(out_fig_dir, "Figure_3B_perplant_size_std.png"),
       plot = panel_B, width = 15, height = 11, units = "cm", dpi = 900, bg = "white")
message("Figure 3B written to ", file.path(out_fig_dir, "Figure_3B_perplant_size_std.png"))

# ---- Assemble full corrected Figure 3 if Panel A is available ---------------
panel_A_path <- file.path("objects", "Figure_3A_panel_plant_units.rds")
if (file.exists(panel_A_path)) {
  panel_A_plant <- readRDS(panel_A_path)
  fig3 <- (panel_A_plant | panel_B) + plot_layout(widths = c(1.1, 1), guides = "keep")
  ggsave(file.path(out_fig_dir, "Figure_3_corrected.png"),
         plot = fig3, width = FIG_W, height = FIG_H, units = "cm", dpi = 900, bg = "white")
  message("Full Figure 3 written to ", file.path(out_fig_dir, "Figure_3_corrected.png"))
} else {
  message("Panel A object not found; run the Figure 3A block of iNEXT3D_plant_level.R to assemble the full figure.")
}
