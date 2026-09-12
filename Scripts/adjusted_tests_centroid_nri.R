#!/usr/bin/env Rscript
# =============================================================================
# adjusted_tests_centroid_nri.R
#
# Re-tests the centroid-based turnover model (Table S10) and phylogenetic
# community structure (NRI/NTI; Tables S19-S20) with tests that match the
# design.
#
# Why
#   - Both published analyses used anova() on lm(), i.e. sequential (Type I)
#     tests with habitat entered first. The habitat test was therefore NOT
#     adjusted for site, in a design where habitat and site are partly
#     confounded (MA has no forest, BEL no paramo). (Site, entered second, was
#     adjusted for habitat; the centroid interaction test, entered last, is
#     unaffected.)
#   - The appendix NRI/NTI code uses tree_ps (120 root samples), but the
#     published tables have n = 60 plants, so the code does not reproduce the
#     output. NRI/NTI are recomputed here on ps_individual (plants) with the
#     same tree and null-model settings.
#   - Habitat varies among the 12 sampling locations, not among plants, so a
#     location-level sensitivity analysis (location means; n = 12) is added.
#
# For each response:
#   1. Published test (Type I, as in the appendix) — reproduction check
#   2. Habitat adjusted for site and site adjusted for habitat (drop1 on the
#      additive model)
#   3. Habitat x site interaction (additive vs interaction model)
#   4. Tukey pairwise habitat contrasts from the additive model
#   5. Location-level sensitivity: same additive tests on location means
#
# Run from the repository root:
#   Rscript Scripts/adjusted_tests_centroid_nri.R          # 999 null runs (as published)
#   Rscript Scripts/adjusted_tests_centroid_nri.R 99       # quick test
# The robust CLR is computed explicitly (see section A), not with vegan's decostand().
# Per-plant NRI/NTI are cached in objects/nri_nti_per_plant_<runs>.rds and
# reused on later runs.
#
# Inputs : objects/ps_individual.rds, objects/tree_ps.rds
# Outputs: tables/adjusted_tests_centroid_nri.csv
#          tables/adjusted_pairwise_centroid_nri.csv
#          tables/nri_nti_per_plant.csv, tables/nri_nti_median_iqr_per_plant.csv
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(picante)
  library(ape)
  library(phangorn)
  library(emmeans)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
RUNS <- if (length(args) >= 1) as.integer(args[1]) else 999L
MIN_TOTAL_TAXON_READS <- 10   # as in the appendix

out_tab_dir <- "tables"; out_obj_dir <- "objects"
dir.create(out_tab_dir, showWarnings = FALSE)
hab_levels <- c("forest", "subparamo", "paramo")

ps_individual <- readRDS(file.path(out_obj_dir, "ps_individual.rds"))
tree_ps       <- readRDS(file.path(out_obj_dir, "tree_ps.rds"))

X <- as(otu_table(ps_individual), "matrix")
if (taxa_are_rows(ps_individual)) X <- t(X)        # plants x OTUs
storage.mode(X) <- "numeric"

md <- as(sample_data(ps_individual), "data.frame")[rownames(X), , drop = FALSE]
stopifnot(all(c("site", "habitat", "site_elevation") %in% names(md)))
md <- data.frame(plant    = rownames(X),
                 site     = factor(as.character(md$site)),
                 habitat  = factor(as.character(md$habitat), levels = hab_levels),
                 location = factor(as.character(md$site_elevation)),
                 stringsAsFactors = FALSE)
message("Plants: ", nrow(md), " | locations: ", nlevels(md$location))

# ---- Test helpers -----------------------------------------------------------
tests_for <- function(d, response) {
  f <- function(rhs) as.formula(paste(response, "~", rhs))
  m_add <- lm(f("habitat + site"), data = d)
  m_int <- lm(f("habitat * site"), data = d)
  t1    <- anova(m_int)                               # published (Type I)
  d1    <- drop1(m_add, test = "F")                   # adjusted
  ai    <- anova(m_add, m_int)                        # interaction

  # Location-level sensitivity (location means)
  loc <- d %>% group_by(location, site, habitat) %>%
    summarise(y = mean(.data[[response]]), .groups = "drop")
  d1_loc <- drop1(lm(y ~ habitat + site, data = loc), test = "F")

  bind_rows(
    data.frame(test = "habitat, Type I (published)", df1 = t1["habitat", "Df"],
               df2 = t1["Residuals", "Df"], F = t1["habitat", "F value"], p = t1["habitat", "Pr(>F)"]),
    data.frame(test = "site, Type I (published)", df1 = t1["site", "Df"],
               df2 = t1["Residuals", "Df"], F = t1["site", "F value"], p = t1["site", "Pr(>F)"]),
    data.frame(test = "habitat adjusted for site", df1 = d1["habitat", "Df"],
               df2 = df.residual(m_add), F = d1["habitat", "F value"], p = d1["habitat", "Pr(>F)"]),
    data.frame(test = "site adjusted for habitat", df1 = d1["site", "Df"],
               df2 = df.residual(m_add), F = d1["site", "F value"], p = d1["site", "Pr(>F)"]),
    data.frame(test = "habitat x site interaction", df1 = ai$Df[2],
               df2 = ai$Res.Df[2], F = ai$F[2], p = ai$`Pr(>F)`[2]),
    data.frame(test = "habitat adjusted for site, location means (n = 12)", df1 = d1_loc["habitat", "Df"],
               df2 = nrow(loc) - 1 - 2 - (nlevels(droplevels(loc$site)) - 1),
               F = d1_loc["habitat", "F value"], p = d1_loc["habitat", "Pr(>F)"]),
    data.frame(test = "site adjusted for habitat, location means (n = 12)", df1 = d1_loc["site", "Df"],
               df2 = nrow(loc) - 1 - 2 - (nlevels(droplevels(loc$site)) - 1),
               F = d1_loc["site", "F value"], p = d1_loc["site", "Pr(>F)"])
  ) %>% mutate(response = response, .before = 1)
}

pairwise_for <- function(d, response) {
  m_add <- lm(as.formula(paste(response, "~ habitat + site")), data = d)
  as.data.frame(summary(pairs(emmeans(m_add, ~ habitat), adjust = "tukey"), infer = c(TRUE, TRUE))) %>%
    mutate(response = response, .before = 1)
}

# =============================================================================
# A. Centroid-based turnover model (replicates the appendix code)
# =============================================================================
# Standard robust CLR (log of non-zero counts centred on their row mean; zeros
# stay zero), as described in the Methods. Computed explicitly because
# vegan 2.7.x decostand(method = "rclr") no longer reproduces it, whereas this
# definition reproduces the published Tables S5-S10 exactly.
rclr_std <- function(x) {
  lx <- log(x); lx[!is.finite(lx)] <- NA
  out <- sweep(lx, 1, rowMeans(lx, na.rm = TRUE)); out[is.na(out)] <- 0
  out
}
rfy <- rclr_std(X)
centroids <- rowsum(rfy, group = as.character(md$site)) / as.vector(table(as.character(md$site)))
md$dist_centroid <- vapply(seq_len(nrow(rfy)), function(i)
  sqrt(sum((rfy[i, ] - centroids[as.character(md$site[i]), ])^2)), numeric(1))

# =============================================================================
# B. NRI / NTI per plant
# =============================================================================
cache <- file.path(out_obj_dir, paste0("nri_nti_per_plant_", RUNS, ".rds"))
if (file.exists(cache)) {
  message("Loading cached NRI/NTI: ", cache)
  nri <- readRDS(cache)
} else {
  tr <- phy_tree(tree_ps)
  if (!ape::is.rooted(tr)) tr <- phangorn::midpoint(tr)
  comm <- X[, intersect(colnames(X), tr$tip.label), drop = FALSE]
  comm <- comm[, colSums(comm) >= MIN_TOTAL_TAXON_READS, drop = FALSE]
  tr   <- ape::keep.tip(tr, colnames(comm))
  dist_phy <- cophenetic(tr)
  message("NRI/NTI: ", nrow(comm), " plants x ", ncol(comm), " taxa; ", RUNS,
          " null runs (this can take a while) ...")
  set.seed(1)
  mpd_aw  <- ses.mpd (comm, dist_phy, null.model = "taxa.labels", runs = RUNS, abundance.weighted = TRUE)
  mntd_aw <- ses.mntd(comm, dist_phy, null.model = "taxa.labels", runs = RUNS, abundance.weighted = TRUE)
  comm_pa <- 1 * (comm > 0)
  mpd_pa  <- ses.mpd (comm_pa, dist_phy, null.model = "taxa.labels", runs = RUNS, abundance.weighted = FALSE)
  mntd_pa <- ses.mntd(comm_pa, dist_phy, null.model = "taxa.labels", runs = RUNS, abundance.weighted = FALSE)
  nri <- data.frame(plant  = rownames(comm),
                    NRI_aw = -mpd_aw$mpd.obs.z,  NTI_aw = -mntd_aw$mntd.obs.z,
                    NRI_pa = -mpd_pa$mpd.obs.z,  NTI_pa = -mntd_pa$mntd.obs.z)
  saveRDS(nri, cache)
}
md <- left_join(md, nri, by = "plant")
write.csv(md[, c("plant", "site", "habitat", "location", "dist_centroid",
                 "NRI_aw", "NTI_aw", "NRI_pa", "NTI_pa")],
          file.path(out_tab_dir, "nri_nti_per_plant.csv"), row.names = FALSE)

fmt_med_iqr <- function(x) sprintf("%.3f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
med_tab <- md %>% group_by(habitat) %>%
  summarise(n = n(), NRI_aw = fmt_med_iqr(NRI_aw), NTI_aw = fmt_med_iqr(NTI_aw),
            NRI_pa = fmt_med_iqr(NRI_pa), NTI_pa = fmt_med_iqr(NTI_pa), .groups = "drop")
write.csv(med_tab, file.path(out_tab_dir, "nri_nti_median_iqr_per_plant.csv"), row.names = FALSE)

# =============================================================================
# Tests
# =============================================================================
d_ok <- md[complete.cases(md[, c("dist_centroid", "NRI_aw", "NTI_aw")]), ]
if (nrow(d_ok) < nrow(md)) message("Dropped ", nrow(md) - nrow(d_ok), " plants with missing NRI/NTI")

responses <- c("dist_centroid", "NRI_aw", "NTI_aw")
tests <- bind_rows(lapply(responses, function(r) tests_for(d_ok, r)))
pw    <- bind_rows(lapply(responses, function(r) pairwise_for(d_ok, r)))
write.csv(tests, file.path(out_tab_dir, "adjusted_tests_centroid_nri.csv"), row.names = FALSE)
write.csv(pw,    file.path(out_tab_dir, "adjusted_pairwise_centroid_nri.csv"), row.names = FALSE)

# ---- Console summary --------------------------------------------------------
fmt <- function(x) signif(x, 3)
cat("\n==== Reproduction check: NRI/NTI median (IQR) per habitat ====\n")
cat("Published Table S19: forest NRI_aw 0.991 (3.02), NTI_aw 3.053 (1.70); paramo 1.435 (1.85), 3.277 (1.62);\n",
    "                     subparamo 1.675 (1.77), 2.831 (1.40)\n")
print(as.data.frame(med_tab), row.names = FALSE)
cat("Published Table S10 (centroid, Type I): habitat F(2,50) = 14.7; site F(3,50) = 3.60; habitat x site F(4,50) = 1.28, p = 0.29\n")
cat("Published Table S20 (Type I): NTI_aw habitat F = 0.51, p = 0.60; site F = 1.82, p = 0.16;",
    "NRI_aw habitat F = 1.27, p = 0.29; site F = 3.42, p = 0.024\n")

cat("\n==== Tests ====\n")
print(as.data.frame(tests %>% mutate(F = fmt(F), p = fmt(p))), row.names = FALSE)

cat("\n==== Tukey pairwise habitat contrasts (additive model) ====\n")
lo <- intersect(c("lower.CL", "asymp.LCL"), names(pw))[1]
hi <- intersect(c("upper.CL", "asymp.UCL"), names(pw))[1]
print(as.data.frame(pw %>% transmute(response, contrast, estimate = fmt(estimate),
                                      lower = fmt(.data[[lo]]), upper = fmt(.data[[hi]]),
                                      p = fmt(p.value))), row.names = FALSE)
cat("\nDone.\n")
