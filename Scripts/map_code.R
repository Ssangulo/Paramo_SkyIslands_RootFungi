# =============================================================================
# Figure 1 -- study-area map (panel A) + elevation profiles (panel B)
#
# Manuscript: "Parallel elevation filtering of Ericaceae root-associated fungal
#              communities in Andean paramo ecosystems" (Ecology and Evolution)
#
# PROVENANCE ------------------------------------------------------------------
# The original Figure 1 code was lost.  This script was reconstructed from
#   (i)  a surviving fragment of the panel-A zoom map (the `mapa_zoom` ggplot),
#   (ii) Table S1 of Analysis_pipeline.qmd (the only source of site coordinates
#        in the repository -- none of the Data_S*.csv files carry lat/long), and
#   (iii) direct measurement of the published figures/Fig_1_map.png: axis
#        ranges, panel pixel geometry, point and label positions, the elevation
#        colour ramp, and the DEM cell size (~0.011 deg -> elevatr z = 6; the
#        surviving fragment said z = 7, but that fragment is demonstrably an
#        older variant -- it also labels the legend "Main Cities"/"Study Sites"
#        where the published figure reads "Large cities"/"Study sites").
#
# Panel B lives in map_code_panelB.py (matplotlib), which this script calls.
#
# ENVIRONMENT -----------------------------------------------------------------
# Needs an R with the spatial stack.  It was built as a dedicated conda env so
# that r_env (used by the rest of the analysis) is left untouched:
#
#   conda create -y -n r_geo -c conda-forge r-base=4.3.3 r-ggplot2 r-patchwork \
#       r-scales r-png r-sf r-terra r-raster r-elevatr r-rnaturalearth \
#       r-rnaturalearthdata r-progress
#
# Run with:
#   cd <repo root> && conda run -n r_geo Rscript Scripts/map_code.R
#
# The DEM is downloaded once from the AWS terrain tiles and cached in
# objects/dem_colombia_z6.tif (gitignored); later runs are offline and instant.
#
# OUTPUTS ---------------------------------------------------------------------
#   figures/Fig_1A_maps.png          panel A only        (3600 x 1802 px)
#   figures/Fig_1B_profiles.png      panel B only        (3600 x 2400 px)
#   figures/Fig_1_map_regenerated.png  composite         (3600 x 4202 px)
# The published figures/Fig_1_map.png is never modified.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(raster); library(elevatr); library(rnaturalearth)
  library(ggplot2); library(patchwork); library(scales); library(png)
})

HERE     <- tryCatch(dirname(normalizePath(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
              error = function(e) file.path(getwd(), "Scripts"))
if (is.na(HERE) || !nzchar(HERE)) HERE <- file.path(getwd(), "Scripts")
SCRIPT_DIR <- HERE
HERE     <- dirname(SCRIPT_DIR)   # repository root (this script lives in Scripts/)
FIG_DIR  <- file.path(HERE, "figures")
OBJ_DIR  <- file.path(HERE, "objects")

# ---- render geometry --------------------------------------------------------
# The layout is designed on a 3600 px-wide canvas: panel B is 3600 x 2400 px and
# panel A was 1802 px in the published figure, 1830 here so the enlarged axis
# text fits without shrinking the maps.
DPI       <- 300
PANEL_A_W <- 3600; PANEL_A_H <- 1830
PANEL_B_H <- 2400

# ---- print resolution -------------------------------------------------------
# The figure is exported at OUT_DPI for a printed width of OUT_WIDTH_CM, matching
# the dpi = 900 convention used by Figures 2-4 of this repository.  RENDER_SCALE
# multiplies the pixel count only: layout fractions and point sizes are untouched,
# so the figure looks identical, just at print resolution.
#   17.4 cm at 900 dpi -> 6165 px wide (17.4 cm is the double-column width used
#   by Figure 2; at that size the axis tick labels print at about 6 pt).
OUT_DPI      <- 900
OUT_WIDTH_CM <- 17.4
OUT_W_PX     <- round(OUT_WIDTH_CM / 2.54 * OUT_DPI)
RENDER_SCALE <- OUT_W_PX / PANEL_A_W

# ---- sizes, tuned against the published figure ------------------------------
# TEXT_SCALE multiplies every piece of text on panel A at once.  1.00 reproduces
# the published figure; it was raised to 1.20 on 2026-08-21 at the authors'
# request, for legibility at print size.  PANEL_A_H was grown to match so that
# the two map panels keep their published pixel size rather than shrinking.
# It is exported to map_code_panelB.py as FIG1_TEXT_SCALE, so panel B's text
# always matches panel A's -- the two panels' base sizes are pixel-matched at
# scale 1, so scaling both by the same factor keeps them matched.
TEXT_SCALE <- 1.55

BASE_SIZE  <- 9      # theme_minimal(base_size = ) for both maps
SITE_SIZE  <- 2.4    # study-site point size
CITY_SIZE  <- 2.0    # city point size
SITE_LAB   <- 2.77 * TEXT_SCALE   # study-site label size
CITY_LAB   <- SITE_LAB            # city labels drawn at the site-label size
TAG_SIZE   <- 19.3   # the "A" tag (pt)
TAG_XY     <- c(57, 75)   # tag centre, px from the top-left of panel A

# Panel-A layout.  Tuned so that the two map panels land on the same pixels as
# the published figure: left panel x 386..1557, zoom panel x 2048..3119, both
# spanning y 39..1651 in the 3600 x 1802 panel-A canvas.
PLOT_MARGIN   <- c(t = 3.90, r = 5.5, b = 7.70, l = 5.5)   # pt
PAD_L1        <- 34.95    # extra left padding, country map (pt)
PAD_L2        <- 45.97    # extra left padding, zoom map (pt)
BAR_W_PX      <- 62; BAR_H_PX <- 430                      # colourbar, px
# (bar lengthened from the published 308 px: the enlarged legend text made its
#  0..5000 labels crowd each other on a bar that short)
# Text sizes (pt) and legend spacing, all measured off the published figure --
# they are not simply theme_minimal(base_size)'s proportions.
FS_AXIS_TITLE <- 9.44 * TEXT_SCALE
FS_AXIS_TEXT  <- 6.55 * TEXT_SCALE
FS_LEG_TITLE  <- 9.55 * TEXT_SCALE
FS_LEG_TEXT   <- 9.00 * TEXT_SCALE   # raised from 7.44: the published size read
                                     # too small against the enlarged map text
AXIS_TITLE_PAD <- 6.6     # axis.title margin towards the panel (pt)
AXIS_TEXT_PAD  <- 3.2     # axis.text margin towards the panel (pt)
LEG_SPACING   <- 0        # between the colourbar and the point legend (pt)
LEG_KEY_H_PX  <- 38       # point-legend key height (px)
LEG_KEY_GAP   <- 5.5      # vertical gap between point-legend keys (pt)
LAYOUT_WIDTHS <- c(1.00, 1.00)                            # patchwork widths

ZOOM_BBOX <- c(xmin = -77, xmax = -73, ymin = 2, ymax = 8)
DEM_Z     <- 6

# ---- the elevation ramp (verbatim from the surviving fragment) --------------
ELEV_COLS <- c("#006400", "#00CC00", "#FFFF66", "#FF9933", "#CD5C5C", "#FFFFFF")
ELEV_VALS <- c(0, 500, 1000, 2500, 4000, 5000)

# =============================================================================
# 1. Spatial data
# =============================================================================

colombia <- ne_countries(country = "Colombia", scale = "medium", returnclass = "sf")

dem_file <- file.path(OBJ_DIR, sprintf("dem_colombia_z%d.tif", DEM_Z))
if (!file.exists(dem_file)) {
  message("downloading DEM (elevatr z = ", DEM_Z, ") -- this happens once")
  dir.create(OBJ_DIR, showWarnings = FALSE)
  writeRaster(get_elev_raster(locations = colombia, z = DEM_Z, clip = "locations"),
              dem_file, overwrite = TRUE)
}
elev_full <- raster(dem_file)

elev_df_full <- as.data.frame(rasterToPoints(elev_full))
colnames(elev_df_full) <- c("long", "lat", "elev")

elev_df_zoom <- as.data.frame(rasterToPoints(
  mask(crop(elev_full, extent(ZOOM_BBOX)), colombia)))
colnames(elev_df_zoom) <- c("long", "lat", "elev")

# Both maps share one colour mapping (identical limits) so the single collected
# legend is valid for both.  `values` is deliberately scales::rescale(ELEV_VALS)
# with no `from =`: that is what the surviving fragment does, and it is what the
# published colourbar shows -- the ramp anchors sit at 0, 0.1, 0.2, 0.5, 0.8 and
# 1 of the DATA range, not at the literal elevations in ELEV_VALS.
ELEV_LIMITS <- range(elev_df_full$elev, na.rm = TRUE)

# ---- study regions ----------------------------------------------------------
# One point per region = the mean of that region's sampling locations in
# Table S1 (verified against the published figure to within ~1 px for MA, DOM
# and NV; BEL sits ~3 px higher there, within marker-centroid error).
paramos <- data.frame(
  etiqueta = c("BEL", "MA", "DOM", "NV"),
  long     = c(-75.666824, -74.010986, -76.100389, -76.064840),
  lat      = c(  6.633500,   4.555806,   3.734824,   3.531368),
  stringsAsFactors = FALSE
)

# ---- cities -----------------------------------------------------------------
ciudades <- data.frame(
  nombre = c("Bogota", "Medellin", "Cali"),
  long   = c(-74.0721, -75.5812, -76.5320),
  lat    = c(  4.7110,   6.2442,   3.4516),
  stringsAsFactors = FALSE
)

# ---- label offsets (the ifelse ladders of the surviving fragment) -----------
site_hjust <- function(l) ifelse(l == "NV", -0.8, ifelse(l == "DOM", 0,
                         ifelse(l == "BEL", 0.5, ifelse(l == "MA", -0.5, 1))))
site_vjust <- function(l) ifelse(l == "NV",  1.7, ifelse(l == "DOM", -1,
                         ifelse(l == "BEL",  -1, ifelse(l == "MA",  0, 0.5))))
city_hjust <- function(n) ifelse(n == "Cali", 0, -0.2)
# Cali sits just below-left of NV.  At the enlarged label size its label ran
# into NV's on the country map, so it is dropped a further line.
city_vjust <- function(n) ifelse(n == "Cali", 2.6, -1)

# =============================================================================
# 2. Shared map layers
# =============================================================================

map_layers <- function(df) {
  list(
    geom_raster(data = df, aes(x = long, y = lat, fill = elev)),
    scale_fill_gradientn(
      colours   = ELEV_COLS,
      values    = scales::rescale(ELEV_VALS),
      limits    = ELEV_LIMITS,
      name      = "Elevation (m)",
      na.value  = "transparent"
    ),
    geom_sf(data = colombia, fill = NA, colour = "black", linewidth = 0.15),
    geom_point(data = paramos, aes(x = long, y = lat, colour = "Study sites"),
               size = SITE_SIZE, shape = 21, fill = "red"),
    geom_text(data = paramos, aes(x = long, y = lat, label = etiqueta),
              hjust = site_hjust(paramos$etiqueta),
              vjust = site_vjust(paramos$etiqueta),
              fontface = "bold", size = SITE_LAB, colour = "black"),
    geom_point(data = ciudades, aes(x = long, y = lat, colour = "Large cities"),
               size = CITY_SIZE, shape = 21, fill = "blue"),
    geom_text(data = ciudades, aes(x = long, y = lat, label = nombre),
              hjust = city_hjust(ciudades$nombre),
              vjust = city_vjust(ciudades$nombre),
              fontface = "bold", size = CITY_LAB, colour = "black"),
    scale_colour_manual(name = "",
      values = c("Large cities" = "blue", "Study sites" = "red")),
    labs(x = "Longitude", y = "Latitude"),
    guides(fill = guide_colourbar(order = 1, theme = theme(
             legend.key.width  = unit(BAR_W_PX / DPI, "in"),
             legend.key.height = unit(BAR_H_PX / DPI, "in"))),
           colour = guide_legend(order = 2, theme = theme(
             legend.key.height = unit(LEG_KEY_H_PX / DPI, "in")))),
    theme_minimal(base_size = BASE_SIZE),
    theme(legend.position = "right", plot.title = element_blank(),
          axis.title.x = element_text(size = FS_AXIS_TITLE,
                                      margin = margin(t = AXIS_TITLE_PAD)),
          axis.title.y = element_text(size = FS_AXIS_TITLE, angle = 90,
                                      margin = margin(r = AXIS_TITLE_PAD)),
          axis.text.x      = element_text(size = FS_AXIS_TEXT,
                                          margin = margin(t = AXIS_TEXT_PAD)),
          axis.text.y      = element_text(size = FS_AXIS_TEXT,
                                          margin = margin(r = AXIS_TEXT_PAD)),
          legend.title     = element_text(size = FS_LEG_TITLE),
          legend.text      = element_text(size = FS_LEG_TEXT),
          legend.spacing.y = unit(LEG_SPACING, "pt"),
          legend.key.spacing.y = unit(LEG_KEY_GAP, "pt"),
          plot.margin = margin(PLOT_MARGIN["t"], PLOT_MARGIN["r"],
                               PLOT_MARGIN["b"], PLOT_MARGIN["l"], "pt"))
  )
}

# =============================================================================
# 3. The two maps
# =============================================================================

bb <- st_bbox(colombia)

mapa_pais <- ggplot() + map_layers(elev_df_full) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  guides(fill = "none", colour = "none") +
  theme(plot.margin = margin(PLOT_MARGIN["t"], PLOT_MARGIN["r"],
                             PLOT_MARGIN["b"], PLOT_MARGIN["l"] + PAD_L1, "pt"))

mapa_zoom <- ggplot() + map_layers(elev_df_zoom) +
  scale_x_continuous(breaks = seq(ZOOM_BBOX["xmin"], ZOOM_BBOX["xmax"], 1)) +
  scale_y_continuous(breaks = seq(ZOOM_BBOX["ymin"], ZOOM_BBOX["ymax"], 1)) +
  coord_sf(xlim = c(ZOOM_BBOX["xmin"], ZOOM_BBOX["xmax"]),
           ylim = c(ZOOM_BBOX["ymin"], ZOOM_BBOX["ymax"]), expand = FALSE) +
  theme(plot.margin = margin(PLOT_MARGIN["t"], PLOT_MARGIN["r"],
                             PLOT_MARGIN["b"], PLOT_MARGIN["l"] + PAD_L2, "pt"))

panel_A <- (mapa_pais | mapa_zoom) + plot_layout(widths = LAYOUT_WIDTHS)

# The "A" tag is drawn straight onto the device rather than via labs(tag=): a
# ggplot/patchwork tag reserves a row above the panels, which shrinks them.
png(file.path(FIG_DIR, "Fig_1A_maps.png"),
    width  = round(PANEL_A_W * RENDER_SCALE),
    height = round(PANEL_A_H * RENDER_SCALE),
    res    = DPI * RENDER_SCALE, bg = "white")
print(panel_A)
grid::grid.text("A", x = grid::unit(TAG_XY[1] / PANEL_A_W, "npc"),
                     y = grid::unit(1 - TAG_XY[2] / PANEL_A_H, "npc"),
                gp = grid::gpar(fontsize = TAG_SIZE))
invisible(dev.off())
message("wrote panel A")
if (nzchar(Sys.getenv("FIG1_PANEL_A_ONLY"))) quit(save = "no")

# =============================================================================
# 4. Panel B (matplotlib) and the composite
# =============================================================================

message("building panel B ...")
stopifnot(system2("python3", shQuote(file.path(SCRIPT_DIR, "map_code_panelB.py")),
                  env = c(sprintf("FIG1_RENDER_SCALE=%.10f", RENDER_SCALE),
                          sprintf("FIG1_TEXT_SCALE=%.6f", TEXT_SCALE))) == 0L)

A <- readPNG(file.path(FIG_DIR, "Fig_1A_maps.png"))
B <- readPNG(file.path(FIG_DIR, "Fig_1B_profiles.png"))
if (dim(A)[3] != dim(B)[3]) {                      # drop alpha if only one has it
  A <- A[, , 1:3, drop = FALSE]; B <- B[, , 1:3, drop = FALSE]
}
stopifnot(dim(A)[2] == dim(B)[2])

target_h <- round((PANEL_A_H + PANEL_B_H) * RENDER_SCALE)
comp <- array(1, dim = c(target_h, dim(A)[2], dim(A)[3]))
comp[seq_len(dim(A)[1]), , ] <- A                                  # panel A on top
comp[target_h - dim(B)[1] + seq_len(dim(B)[1]), , ] <- B           # panel B below

out <- file.path(FIG_DIR, "Fig_1_map_regenerated.png")
writePNG(comp, out, dpi = OUT_DPI)
message(sprintf("wrote %s  (%d x %d px = %.1f x %.1f cm at %d dpi)",
                out, dim(comp)[2], dim(comp)[1],
                dim(comp)[2] / OUT_DPI * 2.54, dim(comp)[1] / OUT_DPI * 2.54, OUT_DPI))
