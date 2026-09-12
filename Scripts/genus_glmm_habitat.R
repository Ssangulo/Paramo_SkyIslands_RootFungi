#!/usr/bin/env Rscript
# =============================================================================
# genus_glmm_habitat.R
#
# Genus-level habitat responses: replacement for the GLLVM summaries behind
# Figure 4C and Data S2.
#
# Why the GLLVM summaries need replacing
#   1. Coefficient divergence: many OTU-level habitat coefficients are
#      extremely large (genus summaries up to |beta| = 260 on the log scale;
#      Data S2 "% enrichment" up to 2e115). This is the signature of taxa that
#      are absent from the reference habitat (forest), for which the
#      coefficient is not identifiable. Abundance-weighted means of these
#      values are not effect sizes.
#   2. Random effects: row.eff = ~(1|site) + (1|Unique_ID) are community-level
#      row effects shared by all OTUs. There are no taxon-specific site
#      effects, so, because habitat and site are partly confounded (MA has no
#      forest, BEL no paramo), taxa specific to a site can appear as habitat
#      responses.
#   3. Meliniomyces was kept separate from its synonym Hyaloscypha.
#
# What this script does
#   - Aggregates reads to genus per plant (ps_individual; root samples summed),
#     merging Meliniomyces into Hyaloscypha.
#   - For each genus meeting a prevalence threshold, fits a negative binomial
#     GLMM (glmmTMB, nbinom2):
#         count ~ habitat + site + offset(log(library size)) + (1 | location)
#     Site is a genus-specific fixed effect (handles the confounding); the
#     location random intercept makes the 12 sampling locations the effective
#     replicates for habitat.
#   - Omnibus habitat test (likelihood ratio vs model without habitat) and
#     three habitat contrasts as fold changes in relative abundance
#     (subparamo/forest, paramo/forest, paramo/subparamo), with 95% CI;
#     Benjamini-Hochberg q-values across genera for each test.
#   - Flags genera absent from a habitat (contrast not estimable; reported
#     as presence/absence instead) and non-converged fits.
#   - Reports, for each genus, presence by habitat and the direction of the
#     paramo-vs-forest difference within NV and DOM (the two complete
#     gradients), to check "consistent" enrichment claims.
#   - Builds a heatmap (log2 fold change vs forest) as a drop-in p_heat for
#     the Figure 4 scripts.
#
# Run from the repository root:  Rscript Scripts/genus_glmm_habitat.R
# Inputs : objects/ps_individual.rds
# Outputs: tables/genus_glmm_results.csv          (full results; new Data S2)
#          tables/genus_glmm_presence_by_site.csv (presence per site x habitat)
#          objects/p_heat_glmm.rds                (heatmap for Figure 4C)
#          figures/genus_glmm_heatmap.png
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(glmmTMB)
  library(emmeans)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---- Settings ---------------------------------------------------------------
MIN_PLANTS    <- 10     # genus must be present in at least this many plants
MIN_READS     <- 1000   # and have at least this many reads in total
Q_THRESHOLD   <- 0.05   # BH threshold for highlighting
HEATMAP_MAX_N <- 30     # max genera shown in the heatmap
LOG2_CAP      <- 6      # colour scale limit for log2 fold change

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
# Genera shown: any with q < threshold for a vs-forest contrast, plus published
# Figure 4C genera that were modelled; capped by total reads.
show <- res %>%
  filter(converged | absent_from != "",
         published_fig4C |
           (!is.na(subparamo_vs_forest_q) & subparamo_vs_forest_q < Q_THRESHOLD) |
           (!is.na(paramo_vs_forest_q) & paramo_vs_forest_q < Q_THRESHOLD) |
           grepl("not estimable", paramo_vs_forest_status) |
           grepl("not estimable", subparamo_vs_forest_status)) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = HEATMAP_MAX_N)

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
