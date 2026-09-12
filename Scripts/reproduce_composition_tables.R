#!/usr/bin/env Rscript
# =============================================================================
# reproduce_composition_tables.R
#
# Reproduction check of the composition analyses (Tables S5-S10) from the
# current objects/ps_individual.rds, using the appendix code.
#
# Background: under vegan >= 2.7-1, decostand(method = "rclr") and
# vegdist(method = "robust.aitchison") impute zeros by matrix completion
# (optspace) by default, so they no longer give the robust CLR used for the
# published tables (log of non-zero counts centred on their mean, zeros kept
# at zero). This script computes that robust CLR explicitly (rclr_std) and
# reproduces Tables S5-S10 with it. The vegan defaults are also reported for
# comparison.
#
# Run from the repository root:  Rscript Scripts/reproduce_composition_tables.R
# Inputs : objects/ps_individual.rds, tables/Table_S5-S10 CSVs
# Output : console comparison; tables/reproduction_check_composition.csv
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
})

message("vegan ", as.character(packageVersion("vegan")),
        " | phyloseq ", as.character(packageVersion("phyloseq")), " | R ", getRversion())

ps <- readRDS("objects/ps_individual.rds")
fmt <- function(x) signif(x, 4)
rows <- list()
add <- function(table, quantity, published, recomputed) {
  rows[[length(rows) + 1]] <<- data.frame(table = table, quantity = quantity,
                                          published = published, recomputed = recomputed)
}
pub <- function(f) read.csv(file.path("tables", f), check.names = FALSE)

rclr_std <- function(x) {                      # robust CLR, zeros kept at zero
  lx <- log(x); lx[!is.finite(lx)] <- NA
  out <- sweep(lx, 1, rowMeans(lx, na.rm = TRUE)); out[is.na(out)] <- 0
  out
}

X <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) X <- t(X)
sam <- data.frame(sample_data(ps))
sam$habitat <- factor(sam$habitat)
sam$site    <- factor(sam$site)
message("ps_individual: ", nrow(X), " plants x ", ncol(X), " OTUs; total reads ", sum(X))

# ---- Row-order check ---------------------------------------------------------
# The appendix code pairs otu_table rows with sample_data rows by POSITION.
# If the two are stored in different orders, habitat/site labels are attached
# to the wrong plants.
order_ok <- identical(rownames(X), rownames(sam))
cat("\n==== Row order: otu_table vs sample_data identical:", order_ok, "====\n")
if (!order_ok) {
  cat("Same set of names:", setequal(rownames(X), rownames(sam)), "\n")
  cat("First 8 otu_table rows :", head(rownames(X), 8), "\n")
  cat("First 8 sample_data rows:", head(rownames(sam), 8), "\n")
  cat("Rows whose habitat label would be wrong if paired by position:",
      sum(as.character(sam$habitat) != as.character(sam[rownames(X), "habitat"])), "of", nrow(sam), "\n")
}
sam_al <- sam[rownames(X), , drop = FALSE]   # aligned by name

# ---- Tables S5-S6: dbRDA (appendix Section 4.1) -----------------------------
# NOTE: the appendix passes otu_table(ps) directly; under vegan 2.7.x /
# phyloseq 1.54 that errors ("logical subscript too long"), so the plain
# matrix X is used here (same values).
d_ait   <- dist(rclr_std(X))                  # robust Aitchison, no imputation
d_vegan <- vegdist(X, method = "robust.aitchison")   # vegan default (imputes zeros)
cat("vegan default vs standard robust Aitchison, mean relative difference:",
    signif(mean(abs(as.vector(d_vegan) - as.vector(d_ait)) / as.vector(d_ait)), 3), "\n")
mod_ait <- capscale(d_ait ~ habitat + Condition(site), data = sam)
set.seed(1)
a_terms <- anova.cca(mod_ait, by = "terms", permutations = 999)
s5 <- pub("Table_S5_dbRDA_terms.csv")
add("S5 dbRDA", "habitat F", s5$F[1], a_terms$F[1])
add("S5 dbRDA", "habitat variance", s5$Variance[1], a_terms$Variance[1])
add("S5 dbRDA", "residual variance", s5$Variance[2], a_terms$Variance[2])
if (!order_ok) {
  mod_al <- capscale(d_ait ~ habitat + Condition(site), data = sam_al)
  set.seed(1)
  a_al <- anova.cca(mod_al, by = "terms", permutations = 999)
  add("S5 dbRDA (aligned by name)", "habitat F", s5$F[1], a_al$F[1])
  add("S5 dbRDA (aligned by name)", "habitat p", 0.001, a_al$`Pr(>F)`[1])
}

# ---- Tables S7-S9: PERMANOVA and betadisper (Section 4.2) --------------------
bps <- subset_samples(ps, site %in% c("DOM", "NV"))
bps <- subset_samples(bps, site_elevation != "NV_4")
bps <- prune_taxa(taxa_sums(bps) > 0, bps)
Xb  <- as(otu_table(bps), "matrix"); if (taxa_are_rows(bps)) Xb <- t(Xb)
sdf <- data.frame(sample_data(bps)); sdf$site <- factor(sdf$site); sdf$habitat <- factor(sdf$habitat)
order_ok_b <- identical(rownames(Xb), rownames(sdf))
cat("Row order (NV + DOM subset) identical:", order_ok_b, "\n")
sdf_al <- sdf[rownames(Xb), , drop = FALSE]
dists <- dist(rclr_std(Xb))
set.seed(1)
perma <- adonis2(dists ~ site * habitat, by = "terms", data = sdf, permutations = 999)
s7 <- pub("Table_S7_PERMANOVA.csv")
for (i in 1:3) add("S7 PERMANOVA", paste(s7$Source[i], "R2"), s7$R2[i], perma$R2[i])
for (i in 1:3) add("S7 PERMANOVA", paste(s7$Source[i], "F"), s7[["Pseudo.F"]][i], perma$F[i])
if (!order_ok_b) {
  set.seed(1)
  perma_al <- adonis2(dists ~ site * habitat, by = "terms", data = sdf_al, permutations = 999)
  for (i in 1:3) add("S7 PERMANOVA (aligned by name)", paste(s7$Source[i], "R2"), s7$R2[i], perma_al$R2[i])
  for (i in 1:3) add("S7 PERMANOVA (aligned by name)", paste(s7$Source[i], "p"), s7[["Pr..F."]][i], perma_al$`Pr(>F)`[i])
  bd_hab_al <- betadisper(dists, sdf_al$habitat)
  set.seed(1); pt_hab_al <- permutest(bd_hab_al, permutations = 999)
  add("S8 betadisper (aligned by name)", "habitat F", pub("Table_S8_betadisper_test_habitat.csv")$F, pt_hab_al$tab[1, "F"])
}

bd_hab  <- betadisper(dists, sdf$habitat)
bd_site <- betadisper(dists, sdf$site)
set.seed(1); pt_hab  <- permutest(bd_hab,  permutations = 999)
set.seed(1); pt_site <- permutest(bd_site, permutations = 999)
add("S8 betadisper", "habitat F", pub("Table_S8_betadisper_test_habitat.csv")$F, pt_hab$tab[1, "F"])
add("S9 betadisper", "site F",    pub("Table_S9_betadisper_test_site.csv")$F,    pt_site$tab[1, "F"])
s8d <- pub("Table_S8_dispersion_habitat.csv")
for (h in s8d$habitat) add("S8 dispersion", paste(h, "mean distance"),
                          s8d$avg_distance_to_centroid[s8d$habitat == h], bd_hab$group.distances[[h]])

# ---- Table S10: centroid model (Section 4.3), with diagnostics ---------------
meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
rfy  <- rclr_std(X)
cat("decostand kept rownames:", !is.null(rownames(rfy)), "| rfy vs X order identical:",
    identical(rownames(rfy), rownames(X)), "| rfy vs sample_data order identical:",
    identical(rownames(rfy), rownames(meta)), "\n")
if (is.null(rownames(rfy))) rownames(rfy) <- rownames(X)
meta_pos <- meta                                   # as the appendix pairs them (by position)
meta     <- meta[rownames(rfy), , drop = FALSE]    # aligned by name
centroids <- rowsum(rfy, group = meta$site) / as.vector(table(meta$site))
dist_a <- vapply(seq_len(nrow(rfy)), function(i)
  sqrt(sum((rfy[i, ] - centroids[meta$site[i], ])^2)), numeric(1))

# Alternative: distances to site centroid from betadisper on robust Aitchison
# distances (should equal dist_a if both use the same rclr)
bd_all <- betadisper(d_ait, factor(meta$site))
dist_b <- bd_all$distances

s10 <- pub("Table_S10_centroid_model.csv")
fit_s10 <- function(dv) {
  df <- meta; df$y <- dv
  df$site <- factor(df$site); df$habitat <- factor(df$habitat, levels = c("forest", "subparamo", "paramo"))
  anova(lm(y ~ habitat * site, data = df))
}
a_a <- fit_s10(dist_a); a_b <- fit_s10(dist_b)

# Positional pairing, as the appendix code would do if the stopifnot were absent
centroids_pos <- rowsum(rfy, group = meta_pos$site) / as.vector(table(meta_pos$site))
dist_pos <- vapply(seq_len(nrow(rfy)), function(i)
  sqrt(sum((rfy[i, ] - centroids_pos[meta_pos$site[i], ])^2)), numeric(1))
df_pos <- meta_pos; df_pos$y <- dist_pos
df_pos$site <- factor(df_pos$site); df_pos$habitat <- factor(df_pos$habitat, levels = c("forest", "subparamo", "paramo"))
a_pos <- anova(lm(y ~ habitat * site, data = df_pos))
add("S10 centroid (positional pairing)", "habitat F", s10$F.value[1], a_pos["habitat", "F value"])
add("S10 centroid (positional pairing)", "residual mean square", s10$Mean.Square[4], a_pos["Residuals", "Mean Sq"])
add("S10 centroid", "habitat F", s10$F.value[1], a_a["habitat", "F value"])
add("S10 centroid", "site F",    s10$F.value[2], a_a["site", "F value"])
add("S10 centroid", "habitat x site F", s10$F.value[3], a_a["habitat:site", "F value"])
add("S10 centroid", "residual mean square", s10$Mean.Square[4], a_a["Residuals", "Mean Sq"])
add("S10 centroid", "habitat sum of squares", s10$Sum.of.Squares[1], a_a["habitat", "Sum Sq"])

out <- bind_rows(rows) %>%
  mutate(rel_diff = ifelse(is.finite(published) & published != 0,
                           (recomputed - published) / abs(published), NA_real_),
         matches = abs(rel_diff) < 0.05)
write.csv(out, "tables/reproduction_check_composition.csv", row.names = FALSE)

cat("\n==== Published vs recomputed (matches = within 5%) ====\n")
print(as.data.frame(out %>% mutate(across(c(published, recomputed, rel_diff), fmt))), row.names = FALSE)

cat("\n==== Centroid diagnostics ====\n")
cat("Distance to site centroid, appendix method: mean", fmt(mean(dist_a)), "sd", fmt(sd(dist_a)),
    "range", fmt(min(dist_a)), "-", fmt(max(dist_a)), "\n")
cat("Distance to site centroid, betadisper:      mean", fmt(mean(dist_b)), "sd", fmt(sd(dist_b)), "\n")
cat("Correlation between the two:", fmt(cor(dist_a, dist_b)), "\n")
cat("Published residual mean square:", fmt(s10$Mean.Square[4]),
    "(residual SD", fmt(sqrt(s10$Mean.Square[4])), ")\n")
cat("\nMean distance by site x habitat (appendix method):\n")
print(as.data.frame(data.frame(site = meta$site, habitat = meta$habitat, d = dist_a) %>%
  group_by(site, habitat) %>% summarise(n = n(), mean = fmt(mean(d)), sd = fmt(sd(d)), .groups = "drop")),
  row.names = FALSE)
cat("\nType I table, appendix method:\n"); print(a_a)
cat("\nType I table, betadisper distances:\n"); print(a_b)
cat("\nDone.\n")
