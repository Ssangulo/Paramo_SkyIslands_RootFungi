#!/usr/bin/env Rscript
# =============================================================================
# Figure4_pooled_glmm.R
#
# Main-text Figure 4, end to end. Merges the former genus_glmm_habitat.R
# (Part 1) and Figure4_option2_pooled.R (Part 2) into one script, so the
# genus models and the figure they feed cannot drift apart.
#
# Part 1 - Genus-level habitat responses (replaces the GLLVM summaries)
#   Why the GLLVM summaries needed replacing
#     1. Coefficient divergence: many OTU-level habitat coefficients are
#        extremely large (genus summaries up to |beta| = 260 on the log scale;
#        Data S2 "% enrichment" up to 2e115). This is the signature of taxa
#        absent from the reference habitat (forest), for which the coefficient
#        is not identifiable. Abundance-weighted means of these values are not
#        effect sizes.
#     2. Random effects: row.eff = ~(1|site) + (1|Unique_ID) are community-level
#        row effects shared by all OTUs. There are no taxon-specific site
#        effects, so, because habitat and site are partly confounded (MA has no
#        forest, BEL no paramo), taxa specific to a site can appear as habitat
#        responses.
#     3. Meliniomyces was kept separate from its synonym Hyaloscypha.
#   What Part 1 does
#     - Aggregates reads to genus per plant (ps_individual; root samples summed),
#       merging Meliniomyces into Hyaloscypha.
#     - For each genus meeting a prevalence threshold, fits a negative binomial
#       GLMM (glmmTMB, nbinom2):
#           count ~ habitat + site + offset(log(library size)) + (1 | location)
#       Site is a genus-specific fixed effect (handles the confounding); the
#       location random intercept makes the 12 sampling locations the effective
#       replicates for habitat.
#     - Omnibus habitat test (likelihood ratio vs model without habitat) and
#       three habitat contrasts as fold changes in relative abundance
#       (subparamo/forest, paramo/forest, paramo/subparamo), with 95% CI;
#       Benjamini-Hochberg q-values across genera for each test.
#     - Flags genera absent from a habitat (contrast not estimable; reported as
#       presence/absence instead) and non-converged fits.
#     - Reports presence by habitat and the direction of the paramo-vs-forest
#       difference within NV and DOM (the two complete gradients).
#     - Builds the Panel C heatmap (log2 fold change vs forest).
#
# Part 2 - Figure 4 itself, a 2 x 2 panel
#   A  Helotiales relative abundance per plant vs elevation
#   B  Sebacinales relative abundance per plant vs elevation
#      Both from the published additive beta regression (Proportion ~ elevation
#      + site), drawn as ONE site-averaged line with a delta-method 95% CI and
#      per-plant points coloured by site, so the common elevation effect is
#      shown without implying four independently estimated site trends.
#   C  The Part 1 genus heatmap.
#   D  The Panel A trend split into genus-assigned and genus-unassigned
#      Helotiales OTUs, each fitted with the same additive beta regression.
#      Coefficients in A, B and D are all on the logit scale.
#
# Run from the repository root:  Rscript Scripts/Figure4_pooled_glmm.R
#
# Input  : objects/ps_individual.rds   (tracked in the repo; written upstream
#          by full_pipeline_DADA2.R)
# Outputs: figures/Figure_4_option2_pooled.png     (main-text Figure 4)
#          tables/genus_glmm_results.csv           (full results; new Data S2)
#          tables/genus_glmm_presence_by_site.csv  (presence per site x habitat)
#          objects/p_heat_glmm.rds                 (Panel C, also read by
#                                                   archive/Scripts/Figure4_option1_site_slopes.R)
#          figures/genus_glmm_heatmap.png          (Panel C on its own)
#
# Panel D's coefficients reproduce tables/helotiales_fraction_check.csv; that
# file is a cross-check, not an input.
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(glmmTMB)
  library(emmeans)
  library(betareg)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
# PART 1 -- per-genus negative binomial GLMMs (was genus_glmm_habitat.R)
# =============================================================================

# ---- Settings ---------------------------------------------------------------
MIN_PLANTS    <- 10     # genus must be present in at least this many plants
MIN_READS     <- 1000   # and have at least this many reads in total
Q_THRESHOLD   <- 0.05   # BH threshold for highlighting
LOG2_CAP      <- 6      # colour scale limit for log2 fold change

# Rows shown in Figure 4C (see the selection block below)
FOCAL_ORDERS    <- c("Helotiales", "Sebacinales")
MIN_FOCAL_READS <- 40000
KEEP_GENERA     <- "Capronia"

# Genera shown in the published Figure 4C; always reported if they pass the
# prevalence filter, and flagged if they do not.
published_genera <- c("Acephala", "Hyaloscypha", "Pezicula", "Pezoloma", "Gyoerffyella",
                      "Coniochaeta", "Oidiodendron", "Sclerococcum", "Lachnum", "Leohumicola",
                      "Capronia", "Pseudoplectania", "Xenochalara", "Serendipita")
# (Meliniomyces is merged into Hyaloscypha below.)

out_tab_dir <- "tables"; out_fig_dir <- "figures"; out_obj_dir <- "objects"
for (d in c(out_tab_dir, out_fig_dir, out_obj_dir)) dir.create(d, showWarnings = FALSE)

hab_levels <- c("forest", "subparamo", "paramo")
labs_map   <- c(forest = "Forest", subparamo = "Subpáramo", paramo = "Páramo")

message("glmmTMB ", as.character(packageVersion("glmmTMB")),
        " | emmeans ", as.character(packageVersion("emmeans")))

# ---- Data -------------------------------------------------------------------
ps <- readRDS(file.path(out_obj_dir, "ps_individual.rds"))

otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)                 # plants x OTUs
storage.mode(otu) <- "numeric"

meta <- as(sample_data(ps), "data.frame")[rownames(otu), , drop = FALSE]
stopifnot(all(c("site", "habitat", "site_elevation") %in% names(meta)))
meta <- data.frame(plant    = rownames(otu),
                   site     = factor(as.character(meta$site)),
                   habitat  = factor(as.character(meta$habitat), levels = hab_levels),
                   location = factor(as.character(meta$site_elevation)),
                   lib      = rowSums(otu),
                   stringsAsFactors = FALSE)
stopifnot(!anyNA(meta$habitat))
message("Plants: ", nrow(meta), " | locations: ", nlevels(meta$location),
        " | library size range: ", min(meta$lib), "-", max(meta$lib))

# ---- Genus labels -----------------------------------------------------------
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)[colnames(otu), , drop = FALSE]
clean <- function(x) {
  x <- as.character(x)
  x <- sub("^[a-z]__", "", x)
  x[is.na(x)] <- ""
  trimws(x)
}
genus_raw <- clean(tax$Genus)
order_raw <- clean(tax$Order)
order_raw <- sub(" (Kingdom|Phylum|Class|Order)$", "", order_raw)

# A genus is "assigned" if it is not empty, not incertae sedis, and not a
# placeholder filled from a higher rank by tax_fix (e.g. "Helotiales Order").
assigned <- genus_raw != "" &
  !grepl("incertae_sedis", genus_raw, ignore.case = TRUE) &
  !grepl(" (Kingdom|Phylum|Class|Order|Family)$", genus_raw) &
  !grepl("_gen_", genus_raw, fixed = TRUE)

genus <- genus_raw
genus[genus == "Meliniomyces"] <- "Hyaloscypha"      # synonym (Fehrer et al., 2019)
message("OTUs with genus assigned: ", sum(assigned), " of ", length(assigned),
        " (", sum(genus_raw == "Meliniomyces"), " Meliniomyces OTUs merged into Hyaloscypha)")

g_otu   <- otu[, assigned, drop = FALSE]
g_names <- genus[assigned]
g_order <- tapply(order_raw[assigned], g_names, function(v) names(sort(table(v), decreasing = TRUE))[1])
g_notu  <- table(g_names)

gcounts <- t(rowsum(t(g_otu), group = g_names))      # plants x genera
stopifnot(identical(rownames(gcounts), meta$plant))

# ---- Prevalence filter ------------------------------------------------------
n_plants <- colSums(gcounts > 0)
n_reads  <- colSums(gcounts)
keep     <- n_plants >= MIN_PLANTS & n_reads >= MIN_READS
message("Genera: ", ncol(gcounts), " assigned | ", sum(keep), " pass the prevalence filter (>= ",
        MIN_PLANTS, " plants, >= ", MIN_READS, " reads)")

not_modelled <- setdiff(published_genera, colnames(gcounts)[keep])
if (length(not_modelled) > 0) {
  message("Published Figure 4C genera NOT modelled (absent or below prevalence filter): ",
          paste(not_modelled, collapse = ", "))
}

# ---- Presence by habitat and by site x habitat ------------------------------
presence_hab <- sapply(hab_levels, function(h) colSums(gcounts[meta$habitat == h, , drop = FALSE] > 0))
colnames(presence_hab) <- paste0("plants_present_", hab_levels)
n_plants_hab <- table(meta$habitat)

presence_site <- expand.grid(site = levels(meta$site), habitat = hab_levels, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(n_plants = sum(meta$site == site & meta$habitat == habitat)) %>%
  ungroup() %>%
  filter(n_plants > 0)
pres_site_tab <- bind_rows(lapply(seq_len(nrow(presence_site)), function(i) {
  idx <- meta$site == presence_site$site[i] & meta$habitat == presence_site$habitat[i]
  data.frame(genus = colnames(gcounts),
             site = presence_site$site[i], habitat = presence_site$habitat[i],
             n_plants = presence_site$n_plants[i],
             plants_present = colSums(gcounts[idx, , drop = FALSE] > 0),
             mean_rel_abund = colMeans(gcounts[idx, , drop = FALSE] / meta$lib[idx]),
             stringsAsFactors = FALSE)
}))
write.csv(pres_site_tab, file.path(out_tab_dir, "genus_glmm_presence_by_site.csv"), row.names = FALSE)

# Direction of paramo vs forest within the two complete gradients
within_site_dir <- function(g, s) {
  f <- pres_site_tab$mean_rel_abund[pres_site_tab$genus == g & pres_site_tab$site == s & pres_site_tab$habitat == "forest"]
  p <- pres_site_tab$mean_rel_abund[pres_site_tab$genus == g & pres_site_tab$site == s & pres_site_tab$habitat == "paramo"]
  if (length(f) == 0 || length(p) == 0) return(NA_character_)
  if (f == 0 && p == 0) return("absent")
  if (f == 0) return("paramo only")
  if (p == 0) return("forest only")
  if (p > f) "higher in paramo" else "higher in forest"
}

# ---- Per-genus GLMM ---------------------------------------------------------
contrast_list <- list("subparamo / forest" = c(-1, 1, 0),
                      "paramo / forest"    = c(-1, 0, 1),
                      "paramo / subparamo" = c(0, -1, 1))

fit_one <- function(g) {
  d <- meta
  d$y <- gcounts[, g]
  zero_hab <- hab_levels[tapply(d$y, d$habitat, sum)[hab_levels] == 0]

  warn <- character(0)
  fit_quiet <- function(f) withCallingHandlers(
    tryCatch(glmmTMB(f, family = nbinom2, data = d), error = function(e) NULL),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })

  m1 <- fit_quiet(y ~ habitat + site + offset(log(lib)) + (1 | location))
  m0 <- fit_quiet(y ~ site + offset(log(lib)) + (1 | location))

  base <- data.frame(genus = g, order = unname(g_order[g]), n_OTUs = as.integer(g_notu[g]),
                     total_reads = n_reads[g], plants_present = n_plants[g],
                     t(presence_hab[g, ]), absent_from = paste(zero_hab, collapse = ";"),
                     NV_paramo_vs_forest = within_site_dir(g, "NV"),
                     DOM_paramo_vs_forest = within_site_dir(g, "DOM"),
                     stringsAsFactors = FALSE, row.names = NULL)

  if (is.null(m1) || is.null(m0)) {
    return(cbind(base, converged = FALSE, note = "fit failed"))
  }
  pd_ok <- isTRUE(m1$sdr$pdHess)

  lrt <- tryCatch(anova(m0, m1), error = function(e) NULL)
  omni <- if (is.null(lrt)) c(chisq = NA, df = NA, p = NA) else
    c(chisq = lrt$Chisq[2], df = lrt[["Chi Df"]][2], p = lrt[["Pr(>Chisq)"]][2])

  cons <- tryCatch({
    emm <- emmeans(m1, ~ habitat)
    as.data.frame(summary(contrast(emm, method = contrast_list), type = "response",
                          infer = c(TRUE, TRUE), adjust = "none"))
  }, error = function(e) NULL)

  out <- cbind(base, converged = pd_ok,
               note = if (length(warn)) paste(unique(warn), collapse = " | ") else "",
               omnibus_chisq = omni[["chisq"]], omnibus_df = omni[["df"]], omnibus_p = omni[["p"]])

  for (cn in names(contrast_list)) {
    key <- gsub(" / ", "_vs_", cn)
    hab_in <- strsplit(cn, " / ")[[1]]
    estimable <- !any(hab_in %in% zero_hab)
    r <- if (!is.null(cons)) cons[cons$contrast == cn, , drop = FALSE] else NULL
    lo <- intersect(c("lower.CL", "asymp.LCL"), names(r))[1]
    hi <- intersect(c("upper.CL", "asymp.UCL"), names(r))[1]
    ok <- estimable && !is.null(r) && nrow(r) == 1
    out[[paste0(key, "_ratio")]] <- if (ok) r$ratio   else NA_real_
    out[[paste0(key, "_lcl")]]   <- if (ok) r[[lo]]   else NA_real_
    out[[paste0(key, "_ucl")]]   <- if (ok) r[[hi]]   else NA_real_
    out[[paste0(key, "_p")]]     <- if (ok) r$p.value else NA_real_
    out[[paste0(key, "_status")]] <- if (!estimable) paste("not estimable: absent from",
                                                          paste(intersect(hab_in, zero_hab), collapse = "/"))
                                     else if (!pd_ok) "unreliable: non-positive-definite Hessian"
                                     else "estimated"
  }
  out
}

genera_fit <- colnames(gcounts)[keep]
message("Fitting ", length(genera_fit), " genus models ...")
res <- bind_rows(lapply(genera_fit, fit_one))

# Unreliable fits: keep estimates but exclude from multiplicity correction
for (key in c("omnibus", gsub(" / ", "_vs_", names(contrast_list)))) {
  pcol <- if (key == "omnibus") "omnibus_p" else paste0(key, "_p")
  usable <- res$converged & !is.na(res[[pcol]])
  q <- rep(NA_real_, nrow(res))
  q[usable] <- p.adjust(res[[pcol]][usable], method = "BH")
  res[[sub("_p$", "_q", pcol)]] <- q
}

res <- res %>%
  mutate(published_fig4C = genus %in% published_genera) %>%
  arrange(paramo_vs_forest_q, desc(total_reads))
write.csv(res, file.path(out_tab_dir, "genus_glmm_results.csv"), row.names = FALSE)

# ---- Heatmap (drop-in p_heat for Figure 4C) ---------------------------------
# Rows shown. Two a priori clauses and one named exception:
#   (1) genera of the two orders the paper is about (FOCAL_ORDERS) with at least
#       MIN_FOCAL_READS reads. The cut falls between Gyoerffyella (46,851) and
#       Arachnopeziza (32,036), so it does not split a tie. Acephala enters here,
#       not on its fold change, which is not estimable.
#   (2) any genus with a resolved contrast against forest (BH q < Q_THRESHOLD),
#       whatever its order or abundance, so that every significant result in
#       Data S2 is visible in the figure.
#   (3) KEEP_GENERA: Capronia, retained from the published Figure 4C as its most
#       abundant Chaetothyriales genus. This is an EXCEPTION, not a rule -- no
#       abundance threshold admits Capronia (102,051 reads) without also admitting
#       Sclerococcum (293,718), Athelopsis (173,848), Cryptodiscus (160,772) and
#       Teratosphaeria (109,824), none of which has any habitat support.
# The old rule kept all 14 genera of the published GLLVM figure regardless of
# support, which put the largest unresolved fold changes at the top of the panel.
show <- res %>%
  filter(converged | absent_from != "",
         (order %in% FOCAL_ORDERS & total_reads >= MIN_FOCAL_READS) |
           (!is.na(subparamo_vs_forest_q) & subparamo_vs_forest_q < Q_THRESHOLD) |
           (!is.na(paramo_vs_forest_q) & paramo_vs_forest_q < Q_THRESHOLD) |
           genus %in% KEEP_GENERA) %>%
  arrange(desc(total_reads))
message("Figure 4C rows: ", nrow(show), " (", paste(show$genus, collapse = ", "), ")")

heat_df <- bind_rows(lapply(c("subparamo_vs_forest", "paramo_vs_forest"), function(key) {
  data.frame(Genus  = show$genus, Order = show$order,
             contrast = if (key == "subparamo_vs_forest") "Subpáramo vs forest" else "Páramo vs forest",
             log2fc = log2(show[[paste0(key, "_ratio")]]),
             q      = show[[paste0(key, "_q")]],
             status = show[[paste0(key, "_status")]],
             stringsAsFactors = FALSE)
})) %>%
  mutate(
    label = case_when(
      grepl("absent from forest", status) ~ "abs. F",
      grepl("absent from subparamo", status) ~ "abs. S",
      grepl("absent from paramo", status) ~ "abs. P",
      !is.na(q) & q < 0.001 ~ "***",
      !is.na(q) & q < 0.01  ~ "**",
      !is.na(q) & q < 0.05  ~ "*",
      TRUE ~ ""),
    fill_val = pmax(pmin(log2fc, LOG2_CAP), -LOG2_CAP),
    contrast = factor(contrast, levels = c("Subpáramo vs forest", "Páramo vs forest"))
  )

genus_levels <- heat_df %>%
  filter(contrast == "Páramo vs forest") %>%
  arrange(fill_val) %>%
  pull(Genus)
heat_df$Genus <- factor(heat_df$Genus, levels = unique(c(genus_levels, heat_df$Genus)))

p_heat <- ggplot(heat_df, aes(x = contrast, y = Genus, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 2.6) +
  scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#EF4444", midpoint = 0,
                       limits = c(-LOG2_CAP, LOG2_CAP), na.value = "grey85",
                       name = expression(log[2]~"fold change")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"))

saveRDS(p_heat, file.path(out_obj_dir, "p_heat_glmm.rds"))
ggsave(file.path(out_fig_dir, "genus_glmm_heatmap.png"), p_heat,
       width = 12, height = 3 + 0.35 * nrow(show), units = "cm", dpi = 600, bg = "white")

# ---- Console summary --------------------------------------------------------
fmt <- function(x) signif(x, 3)
cat("\n==== Genus GLMM: models fitted =", nrow(res), "| converged =", sum(res$converged), "====\n")
cat("Omnibus habitat effect, q <", Q_THRESHOLD, ":", sum(res$omnibus_q < Q_THRESHOLD, na.rm = TRUE), "genera\n")

cat("\n==== Published Figure 4C genera ====\n")
print(as.data.frame(res %>%
  filter(published_fig4C) %>%
  transmute(genus, order, n_OTUs, plants_present,
            present_F = plants_present_forest, present_S = plants_present_subparamo,
            present_P = plants_present_paramo,
            SvF = fmt(subparamo_vs_forest_ratio), SvF_q = fmt(subparamo_vs_forest_q),
            PvF = fmt(paramo_vs_forest_ratio), PvF_q = fmt(paramo_vs_forest_q),
            NV = NV_paramo_vs_forest, DOM = DOM_paramo_vs_forest, converged)), row.names = FALSE)
if (length(not_modelled) > 0) cat("Not modelled:", paste(not_modelled, collapse = ", "), "\n")

cat("\n==== All genera with q <", Q_THRESHOLD, "for paramo vs forest or subparamo vs forest ====\n")
print(as.data.frame(res %>%
  filter((subparamo_vs_forest_q < Q_THRESHOLD) | (paramo_vs_forest_q < Q_THRESHOLD)) %>%
  transmute(genus, order,
            SvF = fmt(subparamo_vs_forest_ratio), SvF_q = fmt(subparamo_vs_forest_q),
            PvF = fmt(paramo_vs_forest_ratio), PvF_q = fmt(paramo_vs_forest_q),
            PvS = fmt(paramo_vs_subparamo_ratio), PvS_q = fmt(paramo_vs_subparamo_q),
            NV = NV_paramo_vs_forest, DOM = DOM_paramo_vs_forest)), row.names = FALSE)

cat("\n==== Genera absent from a habitat (contrasts not estimable) ====\n")
print(as.data.frame(res %>% filter(absent_from != "") %>%
  transmute(genus, order, absent_from, plants_present,
            plants_present_forest, plants_present_subparamo, plants_present_paramo)), row.names = FALSE)

cat("\nOutputs: tables/genus_glmm_results.csv, tables/genus_glmm_presence_by_site.csv,",
    "objects/p_heat_glmm.rds, figures/genus_glmm_heatmap.png\nDone.\n")


# =============================================================================
# PART 2 -- Figure 4 (was Figure4_option2_pooled.R)
# =============================================================================

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
    plot.tag          = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag.location = "plot",
    plot.tag.position = c(0.005, 0.99),
    strip.background  = element_blank()
  )

inv_logit <- function(x) 1 / (1 + exp(-x))

# ---- Per-plant order proportions (same preparation as Tables S16-S17) -------
ps_individual <- ps   # read once in Part 1

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
  sprintf("Across sites: beta = %+.3f per 100 m\n(logit scale), p = %.3f",
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
pooled_band <- function(model, sites) {
  elev_grid <- seq(x_shared[1], x_shared[2], length.out = 200)
  grid      <- expand.grid(elevation = elev_grid, site = factor(sites, levels = sites))

  mean_terms <- delete.response(terms(model, model = "mean"))
  X    <- model.matrix(mean_terms, grid)
  Xbar <- rowsum(X, match(grid$elevation, elev_grid)) / length(sites)
  beta <- coef(model, model = "mean")
  V    <- vcov(model, model = "mean")
  eta  <- drop(Xbar %*% beta)
  se   <- sqrt(rowSums((Xbar %*% V) * Xbar))
  data.frame(elevation = elev_grid,
             fit = inv_logit(eta),
             lwr = inv_logit(eta - 1.96 * se),
             upr = inv_logit(eta + 1.96 * se))
}

plot_pooled <- function(model, data, tag, y_limits, y_label, effect_label) {
  line_df <- pooled_band(model, levels(data$site))

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
    base_theme
}

panel_A <- plot_pooled(m_hel, df_hel, "A", c(0, 1),    "Relative abundance of Helotiales",  lab_hel)
panel_B <- plot_pooled(m_seb, df_seb, "B", c(0, 0.75), "Relative abundance of Sebacinales", lab_seb)

# ---- Panel D: is the Helotiales trend carried by named genera? --------------
# Splits the per-plant Helotiales read share into genus-assigned and genus-
# unassigned OTUs and fits the same additive beta regression to each. The
# placeholder rule matches helotiales_fraction_check.R, so the coefficients
# reproduce tables/helotiales_fraction_check.csv.
placeholder <- function(x) {
  x <- sub("^[a-z]__", "", as.character(x)); x[is.na(x)] <- ""; x <- trimws(x)
  x == "" | grepl("incertae_sedis", x, ignore.case = TRUE) |
    grepl(" (Kingdom|Phylum|Class|Order|Family)$", x) | grepl("_gen_|_fam_", x)
}

make_split_df <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  sd  <- as(sample_data(ps), "data.frame")[colnames(otu), , drop = FALSE]
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)[rownames(otu), , drop = FALSE]
  hel <- grepl("helotiales", tolower(as.character(tax$Order)), fixed = TRUE)
  hel[is.na(hel)] <- FALSE
  asg <- !placeholder(tax$Genus)
  total <- colSums(otu)
  n <- ncol(otu)
  out <- lapply(c(unassigned = FALSE, assigned = TRUE), function(a) {
    raw <- colSums(otu[hel & (asg == a), , drop = FALSE]) / total
    data.frame(prop_raw = raw,
               Proportion = (raw * (n - 1) + 0.5) / n,
               site = factor(as.character(sd$site), levels = names(site_colors)),
               elevation = as.numeric(sd$elevation),
               n_otu = sum(hel & (asg == a)),
               share_hel = sum(otu[hel & (asg == a), ]) / sum(otu[hel, ]),
               stringsAsFactors = FALSE)
  })
  lapply(out, droplevels)
}

split_df <- make_split_df(ps_individual)
m_unassigned <- betareg(Proportion ~ elevation + site, data = split_df$unassigned)
m_assigned   <- betareg(Proportion ~ elevation + site, data = split_df$assigned)

# Reported on the logit scale, as in panels A-B: betareg's coefficient is a change
# in log-odds, NOT a percentage change in relative abundance.
split_label <- function(model, what, n_otu, share_hel) {
  co <- summary(model)$coefficients$mean["elevation", ]
  sprintf("%s: %d OTUs, %.0f%% of Helotiales reads\nbeta = %+.3f per 100 m (logit scale), p = %.3f",
          what, n_otu, 100 * share_hel, 100 * co[["Estimate"]], co[["Pr(>|z|)"]])
}
lab_unassigned <- split_label(m_unassigned, "Genus unassigned",
                              split_df$unassigned$n_otu[1], split_df$unassigned$share_hel[1])
lab_assigned   <- split_label(m_assigned, "Genus assigned",
                              split_df$assigned$n_otu[1], split_df$assigned$share_hel[1])
message("Helotiales, ", gsub("\n", " | ", lab_unassigned))
message("Helotiales, ", gsub("\n", " | ", lab_assigned))

split_colors <- c(unassigned = "#7B4F9E", assigned = "#BDB0CC")
band_D <- rbind(
  data.frame(pooled_band(m_unassigned, levels(split_df$unassigned$site)), frac = "unassigned"),
  data.frame(pooled_band(m_assigned,   levels(split_df$assigned$site)),   frac = "assigned"))
pts_D <- rbind(
  data.frame(elevation = split_df$unassigned$elevation,
             y = split_df$unassigned$prop_raw, frac = "unassigned", row.names = NULL),
  data.frame(elevation = split_df$assigned$elevation,
             y = split_df$assigned$prop_raw,   frac = "assigned",   row.names = NULL))

# The two annotations are colour-keyed to the lines and carry the series
# descriptions, so panel D needs no separate legend.
panel_D <- ggplot() +
  geom_ribbon(data = band_D, aes(elevation, ymin = lwr, ymax = upr, fill = frac),
              alpha = 0.18, colour = NA) +
  geom_point(data = pts_D, aes(elevation, y, colour = frac), size = 1.5, alpha = 0.65) +
  geom_line(data = band_D, aes(elevation, fit, colour = frac), linewidth = 0.8) +
  scale_colour_manual(values = split_colors, guide = "none") +
  scale_fill_manual(values = split_colors, guide = "none") +
  scale_x_continuous("Elevation (m)", limits = x_shared) +
  scale_y_continuous("Helotiales OTUs,\nshare of total reads per plant",
                     labels = function(x) paste0(100 * x, "%")) +
  coord_cartesian(ylim = c(0, 1.02)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.25, size = 2.2,
           lineheight = 0.95, colour = split_colors[["unassigned"]], label = lab_unassigned) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 3.1, size = 2.2,
           lineheight = 0.95, colour = "grey45", label = lab_assigned) +
  labs(tag = "D") +
  base_theme +
  theme(axis.title.x = element_text(size = 9, margin = margin(t = 1, b = 0)),
        plot.margin = margin(t = 2, r = 5, b = 1, l = 2))

# ---- Panel C: genus-level GLMM heatmap ---------------------------------------
# p_heat is the heatmap built in Part 1 (also saved to objects/p_heat_glmm.rds
# for archive/Scripts/Figure4_option1_site_slopes.R).
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
    plot.margin       = margin(t = 2, r = 2, b = 0, l = 2),
    plot.tag          = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag.location = "plot",
    plot.tag.position = c(0.005, 0.99),
    legend.title      = element_text(size = 8.5),
    legend.text       = element_text(size = 7),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.4, "cm")
  )

# ---- Panel A/B layout tweaks (as published) ---------------------------------
panel_A <- panel_A +
  scale_x_continuous("", limits = x_shared) +      # B carries the shared x title
  theme(legend.position = "none",
        plot.margin = margin(t = 2, r = 5, b = 0, l = 2))

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
    plot.margin           = margin(t = 0, r = 5, b = 1, l = 2)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = "none")

fig4 <- ((panel_A / panel_B) |
         ((panel_C / panel_D) + plot_layout(heights = c(1, 1.5)))) +
  plot_layout(widths = c(1, 1.15))

# ---- Save -------------------------------------------------------------------
# ---- Main-text figure export: PNG + JPEG + EPS ------------------------------
# EPS is written through cairo, which keeps text and lines as vectors; layers
# with alpha are rasterised at fallback_resolution, since the EPS format has no
# transparency.
save_main_fig <- function(stem, plot, width, height, units = "cm", dpi = 900) {
  has_ragg <- requireNamespace("ragg", quietly = TRUE)
  ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height,
         units = units, dpi = dpi, bg = "white",
         device = if (has_ragg) ragg::agg_png else "png")
  ggsave(paste0(stem, ".jpeg"), plot = plot, width = width, height = height,
         units = units, dpi = dpi, bg = "white", quality = 95,
         device = if (has_ragg) ragg::agg_jpeg else "jpeg")
  ggsave(paste0(stem, ".eps"), plot = plot, width = width, height = height,
         units = units, bg = "white",
         device = grDevices::cairo_ps, fallback_resolution = 600)
  message("Written: ", stem, ".png / .jpeg / .eps")
}

save_main_fig(file.path(out_fig_dir, "Figure_4_option2_pooled"), fig4,
              width = 22, height = 17)
