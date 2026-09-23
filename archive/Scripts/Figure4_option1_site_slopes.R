#!/usr/bin/env Rscript
# =============================================================================
# Figure4_option1_site_slopes.R
#
# Figure 4 with panels A-B redrawn from beta regressions with SITE-SPECIFIC
# elevation slopes (Proportion ~ elevation * site). The published panels used
# Proportion ~ elevation + site, which forces a common slope, so the site lines
# were parallel by construction. Panel C is the genus-level GLMM heatmap.
#
# Annotations give the across-site elevation effect from the published
# additive model (as a change in odds per 100 m) and the likelihood-ratio test
# of elevation x site.
#
# SUPERSEDED (reference only; archived). Run from the repository root:
#   Rscript archive/Scripts/Figure4_option1_site_slopes.R
# Inputs : objects/ps_individual.rds, objects/p_heat_glmm.rds (from Figure4_pooled_glmm.R)
# Output : figures/Figure_4_option1_site_slopes.png
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(betareg)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

out_fig_dir <- "figures"
out_obj_dir <- "objects"
dir.create(out_fig_dir, showWarnings = FALSE)

# ---- Shared aesthetics (as published Figure 4) ------------------------------
site_colors <- c(BEL = "#E69F00", DOM = "#009E73", MA = "#56B4E9", NV = "#D55E00")
site_labels <- c(BEL = "Belmira", DOM = "Las Domínguez", MA = "Matarredonda", NV = "La Nevera")

base_theme <- theme_classic(base_size = 11) +
  theme(
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 7, color = "black"),
    legend.title     = element_text(size = 7, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.8, "lines"),
    legend.position  = "right",
    plot.tag         = element_text(size = 13, face = "bold"),
    strip.background = element_blank()
  )

inv_logit <- function(x) 1 / (1 + exp(-x))

# ---- Per-plant order proportions (same preparation as Tables S16-S17) -------
ps_individual <- readRDS(file.path(out_obj_dir, "ps_individual.rds"))

make_order_df <- function(ps, order_name) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  sd  <- as(sample_data(ps), "data.frame")[colnames(otu), , drop = FALSE]
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  ord <- tolower(as.character(tax[rownames(otu), "Order"]))
  ord[is.na(ord)] <- ""
  counts <- colSums(otu[grepl(tolower(order_name), ord, fixed = TRUE), , drop = FALSE])
  total  <- colSums(otu)
  df <- data.frame(plant     = colnames(otu),
                   prop_raw  = counts / total,
                   site      = factor(as.character(sd$site), levels = names(site_colors)),
                   elevation = as.numeric(sd$elevation),
                   stringsAsFactors = FALSE)
  n <- nrow(df)
  df$Proportion <- (df$prop_raw * (n - 1) + 0.5) / n   # Smithson-Verkuilen, as published
  droplevels(df)
}

df_hel <- make_order_df(ps_individual, "helotiales")
df_seb <- make_order_df(ps_individual, "sebacinales")
message("Plants: ", nrow(df_hel))

lr_test <- function(m0, m1) {
  ll0 <- logLik(m0); ll1 <- logLik(m1)
  stat <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  df   <- attr(ll1, "df") - attr(ll0, "df")
  list(chisq = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
}

# Across-site elevation effect from the published additive model, as odds change per 100 m
pooled_label <- function(m_add) {
  co <- summary(m_add)$coefficients$mean["elevation", ]
  sprintf("Across sites: beta = %+.3f per 100 m (logit scale), p = %.3f",
          100 * co[["Estimate"]], co[["Pr(>|z|)"]])
}

x_shared <- range(df_hel$elevation, na.rm = TRUE)

# =============================================================================
# Option 1: site-specific slopes (Proportion ~ elevation * site)
# =============================================================================
fit_models <- function(df) {
  list(add = betareg(Proportion ~ elevation + site, data = df),   # published model
       int = betareg(Proportion ~ elevation * site, data = df))   # site-specific slopes
}
m_hel <- fit_models(df_hel)
m_seb <- fit_models(df_seb)
lr_hel <- lr_test(m_hel$add, m_hel$int)
lr_seb <- lr_test(m_seb$add, m_seb$int)

make_label <- function(m, lr) {
  paste0(pooled_label(m$add), "\n",
         sprintf("Elevation \u00D7 site: \u03C7\u00B2(%d) = %.2f, p = %.3f", lr$df, lr$chisq, lr$p))
}
lab_hel <- make_label(m_hel, lr_hel)
lab_seb <- make_label(m_seb, lr_seb)
message("Helotiales: ", gsub("\n", " | ", lab_hel))
message("Sebacinales: ", gsub("\n", " | ", lab_seb))

plot_site_slopes <- function(model, data, tag, y_limits, y_label, effect_label) {
  newdat <- data %>%
    group_by(site) %>%
    summarize(elev_min = min(elevation), elev_max = max(elevation), .groups = "drop") %>%
    rowwise() %>%
    mutate(elevation = list(seq(elev_min, elev_max, length.out = 100))) %>%
    unnest(elevation) %>%
    ungroup() %>%
    select(site, elevation)
  newdat$site <- factor(as.character(newdat$site), levels = levels(data$site))

  mean_terms <- delete.response(terms(model, model = "mean"))
  X    <- model.matrix(mean_terms, newdat)
  beta <- coef(model, model = "mean")
  V    <- vcov(model, model = "mean")
  eta  <- drop(X %*% beta)
  se   <- sqrt(rowSums((X %*% V) * X))
  newdat$fit <- inv_logit(eta)
  newdat$lwr <- inv_logit(eta - 1.96 * se)
  newdat$upr <- inv_logit(eta + 1.96 * se)

  ggplot() +
    geom_point(data = data, aes(x = elevation, y = prop_raw, colour = site),
               size = 1.8, alpha = 0.7) +
    geom_ribbon(data = newdat, aes(x = elevation, ymin = lwr, ymax = upr, fill = site, group = site),
                alpha = 0.15, colour = NA) +
    geom_line(data = newdat, aes(x = elevation, y = fit, colour = site, group = site),
              linewidth = 0.75) +
    scale_color_manual(values = site_colors, labels = site_labels, name = "Site") +
    scale_fill_manual(values = site_colors, guide = "none") +
    scale_y_continuous(y_label) +
    coord_cartesian(ylim = y_limits) +
    annotate("text", x = Inf, y = Inf, label = effect_label,
             hjust = 1.02, vjust = 1.2, size = 2.6, lineheight = 0.95) +
    labs(tag = tag) +
    base_theme +
    theme(plot.tag = element_text(size = 12, face = "bold", hjust = -0.05))
}

panel_A <- plot_site_slopes(m_hel$int, df_hel, "A", c(0, 1),    "Relative abundance of Helotiales",  lab_hel)
panel_B <- plot_site_slopes(m_seb$int, df_seb, "B", c(0, 0.75), "Relative abundance of Sebacinales", lab_seb)

# ---- Panel C: genus-level GLMM heatmap ---------------------------------------
# Genus-level GLMM heatmap from Figure4_pooled_glmm.R (replaces the GLLVM p_heat)
p_heat_path <- file.path(out_obj_dir, "p_heat_glmm.rds")
if (!file.exists(p_heat_path)) stop("Run Figure4_pooled_glmm.R first: ", p_heat_path, " not found")
p_heat <- readRDS(p_heat_path)
stopifnot(inherits(p_heat, "ggplot"))

# Genus label colours mark the two focal orders only (Panel C's selection rule);
# every other order, Chaetothyriales included, stays neutral grey.
genus_levels  <- levels(p_heat$data$Genus)
if (is.null(genus_levels)) genus_levels <- unique(as.character(p_heat$data$Genus))
ord_by_genus  <- tapply(as.character(p_heat$data$Order), as.character(p_heat$data$Genus), `[`, 1)
genus_orders  <- unname(ord_by_genus[genus_levels])
order_colors  <- c(Helotiales = "#7B4F9E", Sebacinales = "#E69F00")
y_axis_colors <- unname(ifelse(genus_orders %in% names(order_colors), order_colors[genus_orders], "grey30"))

panel_C <- p_heat +
  labs(tag = "C") +
  scale_x_discrete(labels = c("Páramo vs forest"    = "Páramo\nvs forest",
                              "Subpáramo vs forest" = "Subpáramo\nvs forest")) +
  theme(
    axis.text.y       = element_text(color = y_axis_colors, size = 8.5, face = "italic", hjust = 1),
    axis.text.x       = element_text(size = 8.5, angle = 0, hjust = 0.5, vjust = 1),
    plot.margin       = margin(t = 5, r = 5, b = 2, l = 5),
    plot.tag          = element_text(size = 13, face = "bold"),
    legend.title      = element_text(size = 8.5),
    legend.text       = element_text(size = 7),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.4, "cm")
  )

# ---- Panel A/B layout tweaks (as published) ---------------------------------
panel_A <- panel_A +
  scale_x_continuous("", limits = x_shared) +
  theme(legend.position = "none")

panel_B <- panel_B +
  scale_x_continuous("Elevation (m)", limits = x_shared) +
  theme(
    legend.position       = "bottom",
    legend.direction      = "horizontal",
    legend.box            = "horizontal",
    legend.justification  = "center",
    legend.title          = element_text(size = 8, face = "bold"),
    legend.text           = element_text(size = 8),
    legend.key.size       = unit(0.7, "lines"),
    axis.title.x          = element_text(size = 9, margin = margin(t = 1, b = 0)),
    legend.box.margin     = margin(t = -6, b = 0, l = 0, r = 0),
    plot.margin           = margin(t = 2, r = 5, b = 1, l = 5)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = "none")

fig4 <- ((panel_A / panel_B) | panel_C) + plot_layout(widths = c(1.5, 0.75))

# ---- Save -------------------------------------------------------------------
ggsave(file.path(out_fig_dir, "Figure_4_option1_site_slopes.png"),
       plot = fig4, width = 22, height = 14, units = "cm", dpi = 900, bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
message("Written: ", file.path(out_fig_dir, "Figure_4_option1_site_slopes.png"))
