#!/usr/bin/env Rscript
# =============================================================================
# helotiales_fraction_check.R
#
# Is the order-level Helotiales increase with elevation (Table S16) carried by
# OTUs that are unassigned at genus level, and by which families?
#
# Splits per-plant Helotiales reads into (i) genus-assigned and (ii) genus-
# unassigned OTUs, and (iii) the main Helotiales families, then fits the
# published beta regression (proportion ~ elevation + site) to each fraction,
# with a likelihood-ratio test of elevation x site and site-specific slopes.
#
# Run from the repository root:  Rscript Scripts/helotiales_fraction_check.R
# Input : objects/ps_individual.rds
# Output: tables/helotiales_fraction_check.csv
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(betareg); library(emmeans); library(dplyr)
})

MIN_FAMILY_SHARE <- 0.01   # families with >= 1% of Helotiales reads
MIN_PLANTS       <- 10

ps  <- readRDS("objects/ps_individual.rds")
otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
storage.mode(otu) <- "numeric"
md  <- as(sample_data(ps), "data.frame")[rownames(otu), ]
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)[colnames(otu), ]

clean <- function(x) { x <- sub("^[a-z]__", "", as.character(x)); x[is.na(x)] <- ""; trimws(x) }
ord <- tolower(clean(tax$Order)); gen <- clean(tax$Genus); fam <- clean(tax$Family)
placeholder <- function(x) x == "" | grepl("incertae_sedis", x, ignore.case = TRUE) |
  grepl(" (Kingdom|Phylum|Class|Order|Family)$", x) | grepl("_gen_|_fam_", x)

hel      <- grepl("helotiales", ord, fixed = TRUE)
assigned <- !placeholder(gen)
fam_ok   <- !placeholder(fam)
total    <- rowSums(otu)

fractions <- list(
  "Helotiales (all)"               = hel,
  "Helotiales, genus assigned"     = hel & assigned,
  "Helotiales, genus unassigned"   = hel & !assigned
)
hel_reads <- sum(otu[, hel])
fam_share <- tapply(colSums(otu[, hel, drop = FALSE]), ifelse(fam_ok[hel], fam[hel], "(family unassigned)"), sum) / hel_reads
for (f in names(fam_share)[fam_share >= MIN_FAMILY_SHARE]) {
  fractions[[paste0("  family: ", f)]] <- hel & (if (f == "(family unassigned)") !fam_ok else fam == f)
}

base <- data.frame(site = factor(as.character(md$site)), elevation = as.numeric(md$elevation))
base$elev_100 <- base$elevation / 100
n <- nrow(base)

lr <- function(m0, m1) {
  s <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
  d <- attr(logLik(m1), "df") - attr(logLik(m0), "df")
  pchisq(s, d, lower.tail = FALSE)
}

res <- bind_rows(lapply(names(fractions), function(nm) {
  cols <- fractions[[nm]]
  p <- rowSums(otu[, cols, drop = FALSE]) / total
  d <- base; d$y <- (p * (n - 1) + 0.5) / n
  out <- data.frame(fraction = nm, n_OTUs = sum(cols),
                    share_of_Helotiales_reads = sum(otu[, cols]) / hel_reads,
                    plants_present = sum(p > 0), mean_prop = mean(p))
  if (sum(p > 0) < MIN_PLANTS) return(cbind(out, note = "too sparse"))
  m_add <- tryCatch(betareg(y ~ elev_100 + site, data = d), error = function(e) NULL)
  m_int <- tryCatch(betareg(y ~ elev_100 * site, data = d), error = function(e) NULL)
  if (is.null(m_add)) return(cbind(out, note = "fit failed"))
  co <- summary(m_add)$coefficients$mean["elev_100", ]
  sl <- tryCatch(as.data.frame(summary(emtrends(m_int, ~ site, var = "elev_100"))), error = function(e) NULL)
  get_sl <- function(s) if (is.null(sl)) NA_real_ else sl[sl$site == s, grep("trend", names(sl))[1]]
  cbind(out,
        odds_change_pct_per_100m = 100 * (exp(co[["Estimate"]]) - 1),
        p_elevation = co[["Pr(>|z|)"]],
        p_elevation_x_site = if (is.null(m_int)) NA_real_ else lr(m_add, m_int),
        slope_BEL = get_sl("BEL"), slope_DOM = get_sl("DOM"),
        slope_MA = get_sl("MA"), slope_NV = get_sl("NV"),
        note = "")
}))

write.csv(res, "tables/helotiales_fraction_check.csv", row.names = FALSE)
cat("\nSite slopes are on the proportion scale (change in relative abundance per 100 m).\n\n")
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
