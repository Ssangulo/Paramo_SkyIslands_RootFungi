# =============================================================================
# iNEXT3D_plant_level.R
#
# Re-runs the habitat-level iNEXT3D analyses (Appendix Section 5.2; Fig. 3A;
# Tables S12-S13) with PLANTS as the incidence sampling units instead of
# root samples, and compares the two.
#
# This was done because the the original run used tree_ps, which keeps both root samples per plant
# as separate incidence units (datatype = "incidence_raw"). iNEXT3D treats sampling 
# units as independent draws; two roots from the same plant are not
# independent, which inflates effective sample size and narrows the bootstrap
# CIs. Here a plant is scored as an incidence for an OTU if the OTU was
# detected in either of its root samples (ps_individual, root samples summed),
# matching the collapsing used for all other analyses requiring independence.
#
# Run from the repository root:
#   Rscript Scripts/iNEXT3D_plant_level.R          # nboot = 500 (as published)
#   Rscript Scripts/iNEXT3D_plant_level.R 20       # quick test run
#
# Inputs : objects/ps_individual.rds  (per plant; root samples summed)
#          objects/tree_ps.rds        (per root sample; carries the ML tree)
#
# Outputs (existing Table_S16/S17 CSVs are NOT overwritten):
#   objects/out_TD_plant.rds, objects/out_PD_plant.rds      iNEXT3D objects
#   objects/out_TD_root_rerun.rds, objects/out_PD_root_rerun.rds
#   tables/Table_S16_iNEXT3D_TD_plant.csv   drop-in replacement for Table S12
#   tables/Table_S17_iNEXT3D_PD_plant.csv   drop-in replacement for Table S13
#   tables/iNEXT3D_estimates_plant_vs_root.csv   asymptotic + coverage-std
#   tables/iNEXT3D_pairwise_habitat_plant_vs_root.csv   pairwise habitat tests
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(iNEXT.3D)
  library(ape)
  library(phangorn)
  library(dplyr)
})

# ---- Settings ---------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
NBOOT <- if (length(args) >= 1) as.integer(args[1]) else 500L
SEED  <- 1L
RUN_ROOT_LEVEL <- TRUE  # rerun root-sample units under identical settings

out_tab_dir <- "tables"
out_obj_dir <- "objects"
dir.create(out_tab_dir, showWarnings = FALSE)
dir.create(out_obj_dir, showWarnings = FALSE)

hab_levels <- c("forest", "subparamo", "paramo")
labs_map   <- c(forest = "Forest", subparamo = "Subpáramo", paramo = "Páramo")
lvl_order  <- unname(labs_map)

message("iNEXT.3D version: ", as.character(packageVersion("iNEXT.3D")),
        " | nboot = ", NBOOT)

# ---- Load objects -----------------------------------------------------------
ps_individual <- readRDS(file.path(out_obj_dir, "ps_individual.rds"))
tree_ps       <- readRDS(file.path(out_obj_dir, "tree_ps.rds"))

# ---- Shared taxa and tree (identical for both levels) -----------------------
tr_full     <- phy_tree(tree_ps)
common_taxa <- Reduce(intersect, list(taxa_names(ps_individual),
                                      taxa_names(tree_ps),
                                      tr_full$tip.label))
message("Taxa: ps_individual = ", ntaxa(ps_individual),
        ", tree_ps = ", ntaxa(tree_ps),
        ", shared and in tree = ", length(common_taxa))

tr <- ape::keep.tip(tr_full, common_taxa)
tr$node.label <- NULL
if (!ape::is.rooted(tr)) tr <- phangorn::midpoint(tr)
reftime <- max(ape::node.depth.edgelength(tr))  # iNEXT3D default reference time
message("PD reference time (tree height): ", signif(reftime, 5))

# ---- Incidence matrices by habitat ------------------------------------------
incidence_by_habitat <- function(ps, taxa) {
  ps  <- prune_taxa(taxa, ps)
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)          # taxa x units
  sd  <- as(sample_data(ps), "data.frame")
  hab <- as.character(sd[colnames(otu), "habitat"])
  if (anyNA(hab) || !all(hab %in% hab_levels)) {
    stop("Unexpected habitat labels: ",
         paste(unique(hab[!hab %in% hab_levels]), collapse = ", "))
  }
  out <- lapply(hab_levels, function(h) {
    M <- (otu[, hab == h, drop = FALSE] > 0) * 1
    storage.mode(M) <- "numeric"
    M
  })
  names(out) <- hab_levels
  out
}

inc_plant <- incidence_by_habitat(ps_individual, common_taxa)
inc_root  <- incidence_by_habitat(tree_ps, common_taxa)

message("Sampling units per habitat (plants): ",
        paste(names(inc_plant), vapply(inc_plant, ncol, 1L), sep = "=", collapse = ", "))
message("Sampling units per habitat (roots):  ",
        paste(names(inc_root), vapply(inc_root, ncol, 1L), sep = "=", collapse = ", "))

# ---- Cross-check: plant incidence == "present in either root" ---------------
# Compares per-habitat OTU incidence frequencies, so it does not depend on
# sample naming. Requires the plant ID column (Unique_ID) in tree_ps.
sd_root <- as(sample_data(tree_ps), "data.frame")
if ("Unique_ID" %in% names(sd_root)) {
  ok <- vapply(hab_levels, function(h) {
    M   <- inc_root[[h]]
    ids <- as.character(sd_root[colnames(M), "Unique_ID"])
    Mp  <- do.call(cbind, lapply(split(seq_len(ncol(M)), ids),
                                 function(j) as.numeric(rowSums(M[, j, drop = FALSE]) > 0)))
    ncol(Mp) == ncol(inc_plant[[h]]) &&
      isTRUE(all.equal(unname(rowSums(Mp)), unname(rowSums(inc_plant[[h]]))))
  }, logical(1))
  if (all(ok)) {
    message("Cross-check passed: ps_individual incidence matches root samples merged by Unique_ID.")
  } else {
    warning("Cross-check FAILED for: ", paste(hab_levels[!ok], collapse = ", "),
            ". ps_individual may not correspond to tree_ps merged by Unique_ID; inspect before use.")
  }
} else {
  message("Unique_ID not found in tree_ps sample data; cross-check skipped.")
}

# ---- Run iNEXT3D (asymptotic + coverage-standardised) -----------------------
# Coverage-standardised estimates for BOTH unit sets are taken at the plant-unit
# default level (Cmin of the plant incidence data, 0.600), so plant and root
# rows in the comparison tables are at the same coverage.
COV_LEVEL <- iNEXT.3D:::check.level(inc_plant, "incidence_raw", "coverage", NULL)
message("Coverage level for standardised estimates: ", signif(COV_LEVEL, 6))

run_inext <- function(inc, label) {
  message("Running iNEXT3D TD (", label, ") ...")
  set.seed(SEED)
  td <- iNEXT3D(data = inc, diversity = "TD", q = c(0, 1, 2),
                datatype = "incidence_raw", nboot = NBOOT)
  message("Running iNEXT3D PD (", label, ") ...")
  set.seed(SEED)
  pd <- iNEXT3D(data = inc, diversity = "PD", q = c(0, 1, 2),
                datatype = "incidence_raw", nboot = NBOOT,
                PDtree = tr, PDreftime = reftime, PDtype = "meanPD")
  message("Running coverage-standardised estimates (", label, ") ...")
  set.seed(SEED)
  td_cov <- estimate3D(data = inc, diversity = "TD", q = c(0, 1, 2),
                       datatype = "incidence_raw", base = "coverage", level = COV_LEVEL, nboot = NBOOT)
  set.seed(SEED)
  pd_cov <- estimate3D(data = inc, diversity = "PD", q = c(0, 1, 2),
                       datatype = "incidence_raw", base = "coverage", level = COV_LEVEL, nboot = NBOOT,
                       PDtree = tr, PDreftime = reftime, PDtype = "meanPD")
  list(TD = td, PD = pd, TD_cov = td_cov, PD_cov = pd_cov)
}

res_plant <- run_inext(inc_plant, "plant units")
saveRDS(res_plant$TD, file.path(out_obj_dir, "out_TD_plant.rds"))
saveRDS(res_plant$PD, file.path(out_obj_dir, "out_PD_plant.rds"))

if (RUN_ROOT_LEVEL) {
  res_root <- run_inext(inc_root, "root-sample units")
  saveRDS(res_root$TD, file.path(out_obj_dir, "out_TD_root_rerun.rds"))
  saveRDS(res_root$PD, file.path(out_obj_dir, "out_PD_root_rerun.rds"))
}

# ---- Drop-in replacements for Tables S12-S13 (same format as Section 5.2) ---
relabel_hab <- function(x) factor(labs_map[as.character(x)], levels = lvl_order)

tab_td <- as.data.frame(res_plant$TD$TDAsyEst)
tab_td$Assemblage <- relabel_hab(tab_td$Assemblage)
tab_td <- tab_td %>%
  mutate(Diversity_order = case_when(
    qTD == 0 ~ "Species richness", qTD == 1 ~ "Shannon diversity",
    qTD == 2 ~ "Simpson diversity", TRUE ~ as.character(qTD))) %>%
  rename(Habitat = Assemblage)
write.csv(tab_td, file.path(out_tab_dir, "Table_S16_iNEXT3D_TD_plant.csv"), row.names = FALSE)

tab_pd <- as.data.frame(res_plant$PD$PDAsyEst) %>% rename(Habitat = Assemblage)
tab_pd$Habitat <- relabel_hab(tab_pd$Habitat)
write.csv(tab_pd, file.path(out_tab_dir, "Table_S17_iNEXT3D_PD_plant.csv"), row.names = FALSE)

# ---- Standardise iNEXT output columns across package versions ---------------
pick_col <- function(df, exact, regex = NULL, what) {
  hit <- exact[exact %in% names(df)]
  if (length(hit) == 0 && !is.null(regex)) hit <- grep(regex, names(df), value = TRUE)
  if (length(hit) == 0) {
    stop("Cannot find ", what, " column. Columns are: ", paste(names(df), collapse = ", "))
  }
  hit[1]
}

q_as_num <- function(x) {
  if (is.numeric(x)) return(x)
  x   <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  dig <- is.na(out) & grepl("[0-9]", x)
  out[dig] <- as.numeric(sub("^\\D*(\\d+).*$", "\\1", x[dig]))
  lx <- tolower(x)
  out[is.na(out) & grepl("richness", lx)] <- 0
  out[is.na(out) & grepl("shannon",  lx)] <- 1
  out[is.na(out) & grepl("simpson",  lx)] <- 2
  if (anyNA(out)) stop("Cannot parse diversity order from: ",
                       paste(unique(x[is.na(out)]), collapse = ", "))
  out
}

# Asymptotic tables: order in qTD/qPD, estimate in TD_asy/PD_asy
std_asy <- function(df, div, level) {
  df  <- as.data.frame(df)
  asm <- pick_col(df, c("Assemblage", "Community"), what = "assemblage")
  ord <- pick_col(df, c("Order.q", paste0("q", div)), what = "order")
  est <- pick_col(df, c(paste0(div, "_asy"), "Estimator", "Estimate"), what = "estimate")
  se  <- pick_col(df, c("s.e.", "Est_s.e.", "se"), what = "s.e.")
  lcl <- pick_col(df, c(paste0("q", div, ".LCL"), "95% Lower", "LCL"), "LCL|Lower", "LCL")
  ucl <- pick_col(df, c(paste0("q", div, ".UCL"), "95% Upper", "UCL"), "UCL|Upper", "UCL")
  data.frame(level = level, basis = "asymptotic", diversity = div,
             q = q_as_num(df[[ord]]), habitat = as.character(df[[asm]]),
             estimate = df[[est]], se = df[[se]], lcl = df[[lcl]], ucl = df[[ucl]],
             coverage = NA_real_, stringsAsFactors = FALSE)
}

# Coverage-standardised tables: order in Order.q, estimate in qTD/qPD
std_cov <- function(df, div, level) {
  df  <- as.data.frame(df)
  asm <- pick_col(df, c("Assemblage", "Community"), what = "assemblage")
  ord <- pick_col(df, c("Order.q", "Order"), what = "order")
  est <- pick_col(df, c(paste0("q", div), "Estimate", "qD"), what = "estimate")
  se  <- pick_col(df, c("s.e.", "se"), what = "s.e.")
  lcl <- pick_col(df, c(paste0("q", div, ".LCL"), "LCL"), "LCL|Lower", "LCL")
  ucl <- pick_col(df, c(paste0("q", div, ".UCL"), "UCL"), "UCL|Upper", "UCL")
  sc  <- pick_col(df, c("SC", "goalSC"), what = "coverage")
  data.frame(level = level, basis = "coverage_standardised", diversity = div,
             q = q_as_num(df[[ord]]), habitat = as.character(df[[asm]]),
             estimate = df[[est]], se = df[[se]], lcl = df[[lcl]], ucl = df[[ucl]],
             coverage = df[[sc]], stringsAsFactors = FALSE)
}

collect <- function(res, level) {
  bind_rows(std_asy(res$TD$TDAsyEst, "TD", level),
            std_asy(res$PD$PDAsyEst, "PD", level),
            std_cov(res$TD_cov, "TD", level),
            std_cov(res$PD_cov, "PD", level))
}

est_all <- collect(res_plant, "plant")
if (RUN_ROOT_LEVEL) est_all <- bind_rows(est_all, collect(res_root, "root_sample"))
est_all <- est_all %>% mutate(ci_width = ucl - lcl)
write.csv(est_all, file.path(out_tab_dir, "iNEXT3D_estimates_plant_vs_root.csv"), row.names = FALSE)

# ---- Pairwise habitat comparisons -------------------------------------------
# ci_overlap: whether 95% bootstrap CIs overlap (conservative criterion).
# z / p_z: Wald test on the difference using bootstrap s.e. (assemblages are
# bootstrapped independently); p_holm: Holm-adjusted across the three pairs.
pairs_list <- combn(hab_levels, 2, simplify = FALSE)

pairwise <- est_all %>%
  group_by(level, basis, diversity, q) %>%
  group_modify(function(g, key) {
    out <- bind_rows(lapply(pairs_list, function(p) {
      a <- g[g$habitat == p[1], ]
      b <- g[g$habitat == p[2], ]
      if (nrow(a) != 1 || nrow(b) != 1) return(NULL)
      z <- (a$estimate - b$estimate) / sqrt(a$se^2 + b$se^2)
      tibble(pair = paste(labs_map[p[1]], "vs", labs_map[p[2]]),
             est_1 = a$estimate, lcl_1 = a$lcl, ucl_1 = a$ucl,
             est_2 = b$estimate, lcl_2 = b$lcl, ucl_2 = b$ucl,
             ratio = a$estimate / b$estimate,
             ci_overlap = a$lcl <= b$ucl & b$lcl <= a$ucl,
             z = z, p_z = 2 * pnorm(-abs(z)))
    }))
    out$p_holm <- p.adjust(out$p_z, method = "holm")
    out
  }) %>%
  ungroup() %>%
  arrange(basis, diversity, q, pair, level)

write.csv(pairwise, file.path(out_tab_dir, "iNEXT3D_pairwise_habitat_plant_vs_root.csv"),
          row.names = FALSE)

# ---- Console summary --------------------------------------------------------
cat("\n==== Pairwise habitat comparisons (CI overlap; Holm-adjusted Wald p) ====\n")
print(as.data.frame(pairwise %>%
        transmute(level, basis, diversity, q, pair,
                  ratio = round(ratio, 2), ci_overlap,
                  p_holm = signif(p_holm, 3))),
      row.names = FALSE)

if (RUN_ROOT_LEVEL) {
  cat("\n==== CI width: plant units relative to root-sample units ====\n")
  width_cmp <- est_all %>%
    select(level, basis, diversity, q, habitat, ci_width) %>%
    tidyr::pivot_wider(names_from = level, values_from = ci_width) %>%
    mutate(width_ratio_plant_over_root = round(plant / root_sample, 2))
  print(as.data.frame(width_cmp), row.names = FALSE)
}

cat("\nCoverage level used for coverage-standardised estimates:\n")
print(unique(est_all[est_all$basis == "coverage_standardised", c("level", "diversity", "coverage")]),
      row.names = FALSE)

cat("\nDone. Outputs written to ", out_tab_dir, "/ and ", out_obj_dir, "/\n", sep = "")

# =============================================================================
# Figure 3A (corrected): iNEXT3D coverage-based curves, plants as incidence units
# Same aesthetic as the published Figure 3A. Uses res_plant from above if it is
# in memory; otherwise loads objects/out_TD_plant.rds and out_PD_plant.rds, so
# this block can also be run on its own without re-running iNEXT3D.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

out_fig_dir <- "figures"
dir.create(out_fig_dir, showWarnings = FALSE)

labs_map  <- c(forest = "Forest", subparamo = "Subpáramo", paramo = "Páramo")
lvl_order <- unname(labs_map)

out_TD_plant <- if (exists("res_plant")) res_plant$TD else readRDS(file.path("objects", "out_TD_plant.rds"))
out_PD_plant <- if (exists("res_plant")) res_plant$PD else readRDS(file.path("objects", "out_PD_plant.rds"))

# Palette and theme copied from the published Figure 3 code
habitat_colors <- c("Forest" = "#332288", "Subpáramo" = "#AA4499", "Páramo" = "#661100")

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

prep_curves <- function(df, metric_name) {
  df <- as.data.frame(df)
  div_col <- intersect(c("qTD", "qPD", "qD"), names(df))[1]
  lcl_col <- intersect(c("qTD.LCL", "qPD.LCL", "qD.LCL"), names(df))[1]
  ucl_col <- intersect(c("qTD.UCL", "qPD.UCL", "qD.UCL"), names(df))[1]
  if (anyNA(c(div_col, lcl_col, ucl_col)) || !"SC" %in% names(df)) {
    stop("Unexpected iNEXT3D columns: ", paste(names(df), collapse = ", "))
  }
  data.frame(
    habitat     = factor(unname(labs_map[as.character(df$Assemblage)]), levels = lvl_order),
    metric      = metric_name,
    q_label     = paste0("q=", df$Order.q),
    coverage    = df$SC,
    diversity   = df[[div_col]],
    ci_low      = df[[lcl_col]],
    ci_up       = df[[ucl_col]],
    Method_plot = case_when(
      grepl("rare",  df$Method, ignore.case = TRUE) ~ "Rarefaction",
      grepl("extra", df$Method, ignore.case = TRUE) ~ "Extrapolation",
      grepl("obs",   df$Method, ignore.case = TRUE) ~ "Observed",
      TRUE ~ as.character(df$Method)
    ),
    stringsAsFactors = FALSE
  )
}

combined_data <- bind_rows(
  prep_curves(out_TD_plant$TDiNextEst$coverage_based, "TD"),
  prep_curves(out_PD_plant$PDiNextEst$coverage_based, "PD")
)

observed_pts <- combined_data %>% filter(Method_plot == "Observed")

# Observed point is added to both segments so the solid and dashed lines join
line_data <- bind_rows(
  combined_data %>% filter(Method_plot != "Observed"),
  observed_pts %>% mutate(Method_plot = "Rarefaction"),
  observed_pts %>% mutate(Method_plot = "Extrapolation")
) %>%
  arrange(metric, q_label, habitat, Method_plot, coverage)

base_iNEXT_plot <- function(data, obs_pts, lines, y_lab, tag_label = NULL) {
  ggplot(data, aes(x = coverage, y = diversity, color = habitat)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_up, fill = habitat),
                alpha = 0.15, color = NA, show.legend = FALSE) +
    geom_line(data = lines,
              aes(linetype = Method_plot, group = interaction(habitat, Method_plot, q_label)),
              linewidth = 0.85, alpha = 0.85) +
    geom_point(data = obs_pts, size = 1.8, alpha = 0.9) +
    facet_wrap(~ q_label, nrow = 1, scales = "free_y") +
    scale_color_manual(values = habitat_colors, name = NULL) +
    scale_fill_manual(values = habitat_colors, guide = "none") +
    scale_linetype_manual(values = c("Rarefaction" = "solid", "Extrapolation" = "dashed"),
                          name = NULL) +
    labs(x = "Coverage", y = y_lab, tag = tag_label) +
    base_theme +
    theme(axis.title = element_text(size = 9), strip.text = element_text(size = 8))
}

panel_A_TD <- base_iNEXT_plot(
  filter(combined_data, metric == "TD"),
  filter(observed_pts,  metric == "TD"),
  filter(line_data,     metric == "TD"),
  y_lab = "Taxonomic diversity", tag_label = "A"
) +
  theme(legend.position = "none")

panel_A_PD <- base_iNEXT_plot(
  filter(combined_data, metric == "PD"),
  filter(observed_pts,  metric == "PD"),
  filter(line_data,     metric == "PD"),
  y_lab = "Phylogenetic diversity (meanPD)"
) +
  guides(linetype = guide_legend(order = 1, nrow = 1),
         color    = guide_legend(order = 2, nrow = 1)) +
  theme(
    legend.position       = "bottom",
    legend.direction      = "horizontal",
    legend.box            = "horizontal",
    legend.justification  = "center",
    legend.title          = element_blank(),
    legend.text           = element_text(size = 10.2, margin = margin(l = 3)),
    legend.key.width      = unit(0.95, "lines"),
    legend.spacing.x      = unit(0.15, "lines"),
    legend.key.spacing.x  = unit(3, "pt"),
    legend.margin         = margin(t = 1, r = 0, b = 0, l = 0)
  )

panel_A_plant <- panel_A_TD / panel_A_PD

# Panel object saved so it can be combined with a revised Panel B later
saveRDS(panel_A_plant, file.path("objects", "Figure_3A_panel_plant_units.rds"))

# Width matches Panel A's share of the published 28 cm Figure 3 (widths 1.1:1)
ggsave(file.path(out_fig_dir, "Figure_3A_iNEXT_plant_units.png"),
       plot = panel_A_plant, width = 15, height = 11, units = "cm", dpi = 900, bg = "white")

message("Figure 3A (plant units) written to ", file.path(out_fig_dir, "Figure_3A_iNEXT_plant_units.png"))
