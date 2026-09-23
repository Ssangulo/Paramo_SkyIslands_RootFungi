#### GLLVM analysis 
library(gllvm)

ps <- alldat.N[[2]]

comm <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) comm <- t(comm)   # samples x taxa

# Metadata
md <- as(sample_data(ps), "data.frame")

# total reads per sample
libsize <- rowSums(comm)

# store in metadata
md$libsize <- libsize
md$log_libsize <- log(libsize)

range(libsize)
all(rownames(comm) == rownames(md))

# Quick peek
head(libsize)
table(md$habitat)
table(md$site)

# keep taxa with >= 50 reads across all samples
keep <- colSums(comm) >= 50
comm_filt <- comm[, keep]

dim(comm)       # before
dim(comm_filt)  # after

#.....



#There is clear overdispersion and unbalanced design - first making a clean data prep from scratch

# ---- settings ----
min_total_reads <- 1500
presence_threshold <- 200
min_prevalence_samples <- 2
habitat_ref <- "forest"
out_prefix <- "gllvm_results"

Y <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) Y <- t(Y)
md <- as(sample_data(ps), "data.frame")

stopifnot(nrow(Y) == nrow(md))
stopifnot(all(rownames(Y) == rownames(md)))

stopifnot(all(c("habitat","site","Unique_ID") %in% names(md)))

# factors + reference
md$habitat   <- droplevels(factor(md$habitat))
md$site      <- droplevels(factor(md$site))
md$Unique_ID <- droplevels(factor(md$Unique_ID))

if (!(habitat_ref %in% levels(md$habitat))) {
  stop(sprintf("habitat_ref='%s' not found in md$habitat levels: %s",
               habitat_ref, paste(levels(md$habitat), collapse=", ")))
}
md$habitat <- relevel(md$habitat, ref = habitat_ref)

# offsets
md$libsize     <- rowSums(Y)
md$log_libsize <- log(md$libsize)

message("Samples per habitat:"); print(table(md$habitat))
message("Samples per site:");    print(table(md$site))
message("Library size summary:"); print(summary(md$libsize))

# ---- filter ultra-rare taxa ----
totals <- colSums(Y)
prev   <- colSums(Y >= presence_threshold)
keep1 <- totals >= min_total_reads
keep2 <- prev   >= min_prevalence_samples
keep  <- keep1 & keep2
Yf <- Y[, keep, drop = FALSE]
message(sprintf("Taxa before: %d | after filtering: %d", ncol(Y), ncol(Yf)))
message(sprintf("Matrix sparsity after filtering: %.2f%% zeros", mean(Yf == 0) * 100))
if (ncol(Yf) == 0) stop("All taxa were filtered out. Loosen thresholds.")

# ---- design (only site + Unique_ID) ----
Xhab        <- data.frame(habitat = md$habitat)
studyDesign <- md[, c("site","Unique_ID")]

## MODEL A: Abundance (counts)
##    NB family + offset(log libsize) AND poisson
## ---------------------------------------


#Individual random effect with unique ID
fit_nb_2 <- gllvm(
  y           = Yf,
  X           = Xhab,  
  formula     = ~ habitat,
  family      = "negative.binomial",
  offset      = md$log_libsize,
  row.eff     = ~(1|site) + (1|Unique_ID),  # site + individual 
  studyDesign = studyDesign,
  num.lv      = 1,      
  sd.errors   = TRUE,
  method      = "VA"
)

saveRDS(fit_nb_2, file = "/data/lastexpansion/_ang/models/fit_nb_2.rds")
#fit_nb_2 <- readRDS("/data/lastexpansion/_ang/models/fit_nb_2.rds")

#Trying 2 latent variables -- worse AIC fit than lv = 1 -> keeping fit_nb_2
fit_nb_2_2 <- gllvm(
  y           = Yf,
  X           = Xhab,  
  formula     = ~ habitat,
  family      = "negative.binomial",
  offset      = md$log_libsize,
  row.eff     = ~(1|site) + (1|Unique_ID),  # site + individual 
  studyDesign = studyDesign,
  num.lv      = 2,      
  sd.errors   = TRUE,
  method      = "VA"
)



fit_nb_2 <- gllvm(
  y           = Yf,
  X           = Xhab,  
  formula     = ~ habitat,
  family      = "negative.binomial",
  offset      = md$log_libsize,
  row.eff     = ~(1|site) + (1|Unique_ID),  # site + individual 
  studyDesign = studyDesign,
  num.lv      = 1,      
  sd.errors   = TRUE,
  method      = "VA"
)


#Repeating model with individual (nested within site_elevation) as random effect 

fit_nb_1 <- gllvm(
  y           = Yf,
  X           = Xhab,  
  formula     = ~ habitat,
  family      = "negative.binomial",
  offset      = md$log_libsize,
  row.eff     = ~(1|site) + (1 | site_elevation:Individual),  # site + individual 
  studyDesign = studyDesign,
  num.lv      = 1,      
  sd.errors   = FALSE,
  method      = "VA"
)


fit_poisson <- gllvm(
  y           = Yf,
  X           = Xhab,
  formula     = ~ habitat,
  family      = poisson(),
  offset      = md$log_libsize,
  row.eff     = ~(1|site),
  studyDesign = studyDesign,
  num.lv      = 0,
  sd.errors   = FALSE,
  method      = "VA"
)

# AIC and log-likelihood
AIC(fit_poisson); logLik(fit_poisson)
AIC(fit_nb); logLik(fit_nb)
AIC(fit_nb_1); logLik(fit_nb_1)


# Compare models directly
anova(fit_poisson, fit_nb)

# Save the default residual diagnostics
png(filename = "/data/lastexpansion/_ang/Plots2/gllvm_fit_nb_2_residuals.png", width = 2400, height = 1200, res = 500)
par(mfrow = c(1,1))
par(mar = c(5,5,2,2))
plot(fit_nb_2, var.colors = 1)   
dev.off()


graphics.off()
png("/data/lastexpansion/_ang/Plots2/gllvm_fit_nb_2_residuals.png",
    width = 3600, height = 2400, res = 300)

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1), ask = FALSE)
for (i in 1:5) {
  plot(fit_nb_2, which = i, var.colors = 1)
}
plot.new()  # 6th empty cell to complete 3x2 layout

dev.off()


#### NB outperforms poisson in AIC and loglikehood

######## 

# PRIMARY ANALYSIS

library(dplyr)
library(stringr)


## 0) full coef table 
coef_nb <- as.data.frame(coef(fit_nb_2)$Xcoef) %>%
  tibble::rownames_to_column("taxon")

## 1) raw taxonomy
tax <- as.data.frame(tax_table(ps)) %>%
  tibble::rownames_to_column("taxon")

## 1b) APPLY FILTER: keep only taxa assigned at Genus level (g__), exclude Incertae_sedis
taxa_genus_assigned <- grepl("^g__", tax$Genus)
exclude_incertae    <- !grepl("incertae_sedis", tax$Genus, ignore.case = TRUE)
keep_taxa           <- taxa_genus_assigned & exclude_incertae

# filter tax and coef_nb to those taxa only
tax_filt    <- tax[keep_taxa, ]
coef_nb_filt <- coef_nb %>% dplyr::filter(taxon %in% tax_filt$taxon)

## 2) join taxonomy (full), keep only needed ranks
annot_raw <- coef_nb_filt %>%
  left_join(tax_filt %>% select(taxon, Genus, Family, Order, Class, Phylum),
            by = "taxon")

## 3) clean prefixes/suffixes like "o__" and " Family"/" Order"
clean_tax <- function(x) {
  x <- as.character(x)
  x <- gsub("^\\w__","", x)                         # drop rank prefixes (g__, o__, etc.)
  x <- gsub(" (Family|Order|Class|Phylum)$","", x)  # drop trailing rank words if present
  ifelse(is.na(x) | x == "", "Unassigned", x)
}

annot_raw <- annot_raw %>%
  mutate(
    Genus = clean_tax(Genus),
    Order = clean_tax(Order)
  )

## 3) sanity: ensure coef columns are numeric
stopifnot(is.numeric(annot_raw$habitatparamo))
if ("habitatsubparamo" %in% names(annot_raw)) stopifnot(is.numeric(annot_raw$habitatsubparamo))

## 4) per-OTU abundance (weights) — use colSums, not rowSums
otu_abund <- colSums(Yf)  # Yf: samples x taxa used in the model
annot_raw$w_abund <- pmax(1, otu_abund[match(annot_raw$taxon, names(otu_abund))])

## 5) GENUS-LEVEL summary across ALL OTUs
min_n <- 5
genus_summary <- annot_raw %>%
  group_by(Order, Genus) %>%
  summarise(
    n_OTUs   = n(),
    mean_beta_paramo = mean(habitatparamo, na.rm = TRUE),
    q90_beta_paramo  = quantile(habitatparamo, 0.90, na.rm = TRUE),
    max_beta_paramo  = max(habitatparamo, na.rm = TRUE),
    wmean_beta_paramo = weighted.mean(habitatparamo, w = w_abund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_OTUs >= min_n) %>%
  arrange(desc(q90_beta_paramo), desc(wmean_beta_paramo))

## 6) ORDER-LEVEL summary across ALL OTUs
order_summary <- annot_raw %>%
  group_by(Order) %>%
  summarise(
    n_OTUs   = n(),
    median_beta_paramo = median(habitatparamo, na.rm = TRUE),
    q90_beta_paramo    = quantile(habitatparamo, 0.90, na.rm = TRUE),
    max_beta_paramo    = max(habitatparamo, na.rm = TRUE),
    wmean_beta_paramo  = weighted.mean(habitatparamo, w = w_abund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_OTUs >= min_n) %>%
  arrange(desc(q90_beta_paramo), desc(wmean_beta_paramo))

# peek
head(genus_summary, 15)
head(order_summary, 15)

#######

### CLEAN GENUS TABLE
has_subparamo <- "habitatsubparamo" %in% names(annot_raw)
## ---- filtering parameters ----
min_abs_beta_paramo <- 2        # keep genera with |median β_paramo| >= 2
min_total_w_abund   <- 20000    # hard cutoff
min_n_OTUs          <- 3        # keep genera with at least 3 OTUs


## ---- helper for weighted median ----
w_median <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  keep <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(NA_real_)
  o <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}


## ---- genus-level summary table ----
genus_table <- annot_raw %>%
  dplyr::group_by(Order, Genus) %>%
  dplyr::summarise(
    n_OTUs        = dplyr::n(),
    
    # abundance metrics
    total_w_abund = sum(w_abund, na.rm = TRUE),
    mean_w_abund  = mean(w_abund, na.rm = TRUE),
    
    ## --- PÁRAMO: unweighted stats ---
    median_beta_paramo = median(habitatparamo, na.rm = TRUE),
    mean_beta_paramo   = mean(habitatparamo,   na.rm = TRUE),
    q90_beta_paramo    = quantile(habitatparamo, 0.90, na.rm = TRUE),
    max_beta_paramo    = max(habitatparamo,    na.rm = TRUE),
    
    ## --- PÁRAMO: abundance-weighted stats ---
    w_median_beta_paramo = w_median(habitatparamo, w_abund),
    w_mean_beta_paramo   = stats::weighted.mean(habitatparamo, w_abund, na.rm = TRUE),
    
    ## --- SUBPÁRAMO: unweighted stats ---
    median_beta_subparamo =
      if (has_subparamo) median(habitatsubparamo, na.rm = TRUE) else NA_real_,
    mean_beta_subparamo   =
      if (has_subparamo) mean(habitatsubparamo,   na.rm = TRUE) else NA_real_,
    q90_beta_subparamo    =
      if (has_subparamo) quantile(habitatsubparamo, 0.90, na.rm = TRUE) else NA_real_,
    max_beta_subparamo    =
      if (has_subparamo) max(habitatsubparamo,    na.rm = TRUE) else NA_real_,
    
    ## --- SUBPÁRAMO: abundance-weighted stats ---
    w_median_beta_subparamo =
      if (has_subparamo) w_median(habitatsubparamo, w_abund) else NA_real_,
    w_mean_beta_subparamo   =
      if (has_subparamo) stats::weighted.mean(habitatsubparamo, w_abund, na.rm = TRUE) else NA_real_,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ## --- % enrichment (unweighted, as in the old table) ---
    perc_enriched_paramo =
      (exp(median_beta_paramo) - 1) * 100,
    
    perc_enriched_subparamo =
      if (has_subparamo) (exp(median_beta_subparamo) - 1) * 100 else NA_real_,
    
    ## --- % enrichment (weighted medians) ---
    perc_enriched_paramo_w =
      (exp(w_median_beta_paramo) - 1) * 100,
    
    perc_enriched_subparamo_w =
      if (has_subparamo) (exp(w_median_beta_subparamo) - 1) * 100 else NA_real_
  ) %>%
  ## ---- HARD ecological filter + secondary rules (same as old table) ----
dplyr::filter(
  total_w_abund >= min_total_w_abund &
    (
      abs(median_beta_paramo) >= min_abs_beta_paramo |
        n_OTUs >= min_n_OTUs |
        Order == "Helotiales"
    )
) %>%
  dplyr::arrange(dplyr::desc(w_median_beta_paramo))

# quick look
head(genus_order_table, 20)
head(genus_table, 20)

write_xlsx(
  genus_table,
  "/data/lastexpansion/_ang/models/genus_table_NEW.xlsx"
)

write.csv(
  genus_table,
  "/data/lastexpansion/_ang/models/genusr_table_NEW.csv"
)

######## HEATMAP

##  genera to show in heatmap
genera_keep <- c(
  "Acephala",
  "Hyaloscypha",
  "Pezicula",
  "Pezoloma",
  "Gyoerffyella",
  "Meliniomyces",
  "Coniochaeta",
  "Oidiodendron",
  "Sclerococcum",
  "Lachnum",
  "Leohumicola",
  "Capronia",
  "Pseudoplectania",
  "Xenochalara",
  "Serendipita"
)

# order genera by páramo weighted mean beta
genus_order_for_plot <- genus_table %>%
  dplyr::filter(Genus %in% genera_keep) %>%
  dplyr::arrange(dplyr::desc(w_mean_beta_paramo)) %>%
  dplyr::pull(Genus)

# NO abundance labels
genus_label_levels <- genus_order_for_plot

heatmap_df <- genus_table %>%
  dplyr::filter(Genus %in% genera_keep) %>%
  dplyr::select(
    Order, Genus, total_w_abund,
    w_mean_beta_paramo, w_mean_beta_subparamo
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("w_mean_beta_"),
    names_to  = "habitat",
    values_to = "mean_beta"
  ) %>%
  dplyr::mutate(
    habitat_label = dplyr::case_when(
      habitat == "w_mean_beta_paramo"    ~ "Páramo vs forest",
      habitat == "w_mean_beta_subparamo" ~ "Subpáramo vs forest",
      TRUE ~ habitat
    ),
    Genus = factor(Genus, levels = genus_label_levels)
  )

p_heat <- ggplot(heatmap_df,
                 aes(x = habitat_label,
                     y = Genus,
                     fill = mean_beta)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#3B82F6",
    mid      = "white",
    high     = "#EF4444",
    midpoint = 0,
    limits   = c(-35, 35),
    oob      = scales::squish,
    name = expression(paste(beta, " (weighted mean)"))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_blank(),
    axis.text.y  = element_text(size = 9),
    axis.text.x  = element_text(size = 11, angle = 25, hjust = 1),
    legend.position = "right"
  )

ggsave(
  filename = "/data/lastexpansion/_ang/models/heatmap_15_WMEAN.png",
  plot     = p_heat,
  width    = 8,
  height   = 10,
  units    = "in",
  dpi = 900
)

# Save p_heat as RDS for loading in Quarto document
dir.create("objects", showWarnings = FALSE, recursive = TRUE)
saveRDS(p_heat, file = "objects/p_heat_gllvm.rds")


#Making bar plot alongside heatmap with 95% CI 

## labels with n and total abundance (k reads)
genus_labels <- genus_table %>%
  filter(Genus %in% genera_keep) %>%
  mutate(
    Genus_label = sprintf(
      "%s (n=%d, %dk)",
      Genus, n_OTUs, round(total_w_abund / 1000)
    )
  ) %>%
  select(Genus, Genus_label)


## weighted bootstrap CI for β (páramo)
boot_wmean <- function(beta, w, R = 500) {
  beta <- as.numeric(beta)
  w    <- as.numeric(w)
  
  keep <- is.finite(beta) & is.finite(w) & w > 0
  beta <- beta[keep]
  w    <- w[keep]
  
  if (length(beta) < 2) return(c(lo = NA_real_, hi = NA_real_))
  
  prob <- w / sum(w)
  
  boots <- replicate(R, {
    idx <- sample(seq_along(beta), replace = TRUE, prob = prob)
    ww  <- w[idx] / sum(w[idx])
    weighted.mean(beta[idx], ww)
  })
  
  quantile(boots, c(0.025, 0.975), na.rm = TRUE)
}

### --- 3) Compute weighted means + CI ---
genus_ci <- annot_raw %>%
  filter(Genus %in% genera_keep) %>%
  group_by(Genus, Order) %>%
  summarise(
    w_mean_beta = weighted.mean(habitatparamo, w_abund, na.rm = TRUE),
    {
      ci <- boot_wmean(habitatparamo, w_abund, R = 500)
      tibble(ci_lo = ci[1], ci_hi = ci[2])
    },
    .groups = "drop"
  ) %>%
  left_join(genus_labels, by = "Genus")


x_lim <- 100  # visualization window

genus_ci <- genus_ci %>%
  mutate(
    # truncated values for plotting
    plot_x  = pmax(pmin(w_mean_beta, x_lim), -x_lim),
    plot_lo = pmax(pmin(ci_lo,       x_lim), -x_lim),
    plot_hi = pmax(pmin(ci_hi,       x_lim), -x_lim),
    
    # flag genera whose CI extends beyond the window
    extreme_pos = ci_hi > x_lim,
    extreme_neg = ci_lo < -x_lim,
  )

p_bar_w <- ggplot(
  genus_ci,
  aes(
    x = plot_x,
    y = reorder(Genus_label, -plot_x),
    color = Order
  )
) +
  # CI bars (truncated)
  geom_errorbarh(
    aes(xmin = plot_lo, xmax = plot_hi),
    height = 0.3,
    linewidth = 0.8
  ) +
  
  # main points
  geom_point(size = 3) +
  
  # --- Horizontal arrow for extreme positive genus (Coniochaeta) ---
  geom_segment(
    data = subset(genus_ci, extreme_pos),
    aes(
      x = x_lim - 15, 
      xend = x_lim + 5, 
      y = Genus_label,
      yend = Genus_label
    ),
    arrow = arrow(
      length = unit(0.25, "cm"),
      type   = "closed",
      ends   = "last"
    ),
    linewidth = 1,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # vertical zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  
  # styling
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(-x_lim, x_lim + 10)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.y       = element_blank(),
    legend.position    = "right",
    axis.title.x       = element_text(margin = margin(t = 8))
  ) +
  labs(
    x = expression(paste("Weighted mean ", beta, " (Páramo vs forest)")),
    color = "Order"
  )


ggsave(
  filename = "/data/lastexpansion/_ang/models/barplot_15_NEW.png",
  plot = p_bar_w,
  width = 8,
  height = 10,
  units = "in", 
  dpi = 900
)




######## 
# PANEL PLOT #

library(patchwork)

## 
p_heat_tight <- p_heat +
  theme(
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    legend.margin     = margin(t = 0, r = 2, b = 0, l = 0)
  )

p_bar_tight <- p_bar_w +
  theme(
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.margin     = margin(t = 0, r = 2, b = 0, l = 0),
    axis.text.y       = element_text(size = 8),  # keep y labels readable
    # bring x-axis title closer to the axis
    axis.title.x = element_text(
      size = 10,
      margin = margin(t = 0)   # smaller top margin = closer to axis
    ))

p_heat_clean <- p_heat_tight +
  labs(title = NULL, subtitle = NULL) +
  theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank()
  )

p_bar_clean <- p_bar_tight +
  labs(title = NULL, subtitle = NULL) +
  theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank()
  )

panel_AB <- p_heat_clean + p_bar_clean +
  plot_layout(widths = c(0.9, 1.3)) +   # <-- B is wider than A
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ""   # gives A, B
  )


ggsave(
  filename = "/data/lastexpansion/_ang/models/panel_AB_heatmap_barplot_NEW.png",
  plot     = panel_AB,
  width    = 11,
  height   = 6.5,
  units    = "in",
  dpi      = 900
)




