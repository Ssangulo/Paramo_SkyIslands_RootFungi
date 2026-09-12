#!/usr/bin/env Rscript
# =============================================================================
# Figure4_option2_pooled.R
#
# Figure 4 with panels A-B redrawn as ONE across-site trend: the published
# beta regression (Proportion ~ elevation + site), shown as a single
# site-averaged line with 95% CI, and per-plant points coloured by site.
# This displays the common elevation effect without drawing per-site lines
# that are parallel by construction. Panel C is the genus-level GLMM heatmap.
#
# Run from the repository root:  Rscript Scripts/Figure4_option2_pooled.R
# Inputs : objects/ps_individual.rds, objects/p_heat_glmm.rds (from genus_glmm_habitat.R)
# Output : figures/Figure_4_option2_pooled.png
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
# Option 2: one across-site trend (published additive model, site-averaged)
# =============================================================================
m_hel <- betareg(Proportion ~ elevation + site, data = df_hel)   # published model
m_seb <- betareg(Proportion ~ elevation + site, data = df_seb)

lab_hel <- pooled_label(m_hel)
lab_seb <- pooled_label(m_seb)
message("Helotiales: ", lab_hel)
message("Sebacinales: ", lab_seb)

# Site-averaged prediction: the design row is averaged over sites (equal
# weights) at each elevation, so the line carries the common elevation slope
# at the mean site intercept, with a delta-method 95% CI on the logit scale.
plot_pooled <- function(model, data, tag, y_limits, y_label, effect_label) {
  sites     <- levels(data$site)
  elev_grid <- seq(x_shared[1], x_shared[2], length.out = 200)
  grid      <- expand.grid(elevation = elev_grid, site = factor(sites, levels = sites))

  mean_terms <- delete.response(terms(model, model = "mean"))
  X    <- model.matrix(mean_terms, grid)
  Xbar <- rowsum(X, match(grid$elevation, elev_grid)) / length(sites)
  beta <- coef(model, model = "mean")
  V    <- vcov(model, model = "mean")
  eta  <- drop(Xbar %*% beta)
  se   <- sqrt(rowSums((Xbar %*% V) * Xbar))
  line_df <- data.frame(elevation = elev_grid,
                        fit = inv_logit(eta),
                        lwr = inv_logit(eta - 1.96 * se),
                        upr = inv_logit(eta + 1.96 * se))

  ggplot() +
    geom_ribbon(data = line_df, aes(x = elevation, ymin = lwr, ymax = upr),
                fill = "grey40", alpha = 0.18, colour = NA) +
    geom_point(data = data, aes(x = elevation, y = prop_raw, colour = site),
               size = 1.8, alpha = 0.7) +
    geom_line(data = line_df, aes(x = elevation, y = fit),
              colour = "black", linewidth = 0.8) +
    scale_color_manual(values = site_colors, labels = site_labels, name = "Site") +
    scale_y_continuous(y_label) +
    coord_cartesian(ylim = y_limits) +
    annotate("text", x = Inf, y = Inf, label = effect_label,
             hjust = 1.02, vjust = 1.2, size = 2.6) +
    labs(tag = tag) +
    base_theme +
    theme(plot.tag = element_text(size = 12, face = "bold", hjust = -0.05))
}

panel_A <- plot_pooled(m_hel, df_hel, "A", c(0, 1),    "Relative abundance of Helotiales",  lab_hel)
panel_B <- plot_pooled(m_seb, df_seb, "B", c(0, 0.75), "Relative abundance of Sebacinales", lab_seb)

# ---- Panel C: genus-level GLMM heatmap ---------------------------------------
# Genus-level GLMM heatmap from genus_glmm_habitat.R (replaces the GLLVM p_heat)
p_heat_path <- file.path(out_obj_dir, "p_heat_glmm.rds")
if (!file.exists(p_heat_path)) stop("Run genus_glmm_habitat.R first: ", p_heat_path, " not found")
p_heat <- readRDS(p_heat_path)
stopifnot(inherits(p_heat, "ggplot"))

# Genus label colours by order, taken from the heatmap data (other orders grey)
genus_levels  <- levels(p_heat$data$Genus)
if (is.null(genus_levels)) genus_levels <- unique(as.character(p_heat$data$Genus))
ord_by_genus  <- tapply(as.character(p_heat$data$Order), as.character(p_heat$data$Genus), `[`, 1)
genus_orders  <- unname(ord_by_genus[genus_levels])
order_colors  <- c(Helotiales = "#7B4F9E", Sebacinales = "#E69F00", Chaetothyriales = "#FF0000")
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
ggsave(file.path(out_fig_dir, "Figure_4_option2_pooled.png"),
       plot = fig4, width = 22, height = 14, units = "cm", dpi = 900, bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
message("Written: ", file.path(out_fig_dir, "Figure_4_option2_pooled.png"))
