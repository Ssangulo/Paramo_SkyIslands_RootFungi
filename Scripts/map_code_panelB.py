#!/usr/bin/env python3
# =============================================================================
# Figure 1, panel B -- schematic elevation profiles of the four paramo regions
#
# Manuscript: "Parallel elevation filtering of Ericaceae root-associated fungal
#              communities in Andean paramo ecosystems" (Ecology and Evolution)
#
# PROVENANCE ------------------------------------------------------------------
# The original panel-B script was lost.  The render settings below were
# recovered by measuring the published figures/Fig_1_map.png directly
# (sub-pixel axis calibration against the y-tick centroids, then reading fill
# colours, belt boundaries, dashed-line elevations and the terrain polyline out
# of the raster).  See map_code.R for the panel-A counterpart and the
# compositing step.
#
# REVISION (2026-08-21), at the authors' request:
#   * dashed sampling lines drawn thicker;
#   * Matarredonda redrawn as a plateau, with no Andean-forest belt (the site
#     has none) and its summit set to the upper sampling elevation, 3678 m;
#   * Belmira's paramo belt removed (the site has none) and its summit set to
#     the upper sampling elevation, 3254 m;
#   * La Nevera and Las Dominguez unchanged;
#   * optional COMMON_Y mode putting all four panels on one elevation axis.
#
# Runs on the plain system python3 (matplotlib + numpy only; no scipy -- the
# natural cubic spline used for Belmira is implemented inline below).
#
#     python3 Scripts/map_code_panelB.py
#
# -> figures/Fig_1B_profiles.png   (3600 x 2400 px)
# =============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root (script is in Scripts/)
OUT = os.path.join(HERE, "figures", "Fig_1B_profiles.png")

# =============================================================================
# ======================  EDIT BELOW TO RESHAPE THE MOUNTAINS  ================
# =============================================================================

# ---- one elevation axis for all four panels, or one per panel? --------------
# True  -> every panel spans YLIM_COMMON, so the four mountains and their
#          vegetation belts can be compared directly (reviewer 1's suggestion).
# False -> each panel gets its own limits from YLIM_INDEPENDENT, so each
#          mountain fills its panel but absolute elevations are not comparable.
COMMON_Y   = True
YLIM_COMMON = (2300.0, 4400.0)

# Each profile is a list of (x, elevation) control points.  x runs 0 -> 1 across
# the panel and has no units (the panel deliberately has no x axis); elevation
# is metres a.s.l.
#
#   form = "linear"  -> straight segments between the control points
#   form = "spline"  -> natural cubic spline through the control points
#
# Add, drop or move control points freely; the belts, the dashed sampling lines
# and (in COMMON_Y = False mode) the y limits all follow automatically.

PROFILES = {
    # La Nevera -- unchanged from the published figure
    "NEV": dict(
        form="linear",
        nodes=[(0.0000, 2498.9),
               (0.3025, 2766.1),
               (0.7075, 3906.4),
               (1.0000, 4212.3)],
    ),
    # Matarredonda -- REDRAWN as a plateau.  The site has no Andean forest, so
    # the profile starts at the foot of the subparamo belt rather than in
    # forest.  The flat top sits at 3720 m, a little ABOVE the upper sampling
    # line (3678 m, MA_1), so that line meets the slope just below the plateau
    # rather than running along the summit itself.
    "MAT": dict(
        form="linear",
        nodes=[(0.0000, 3200.0),
               (0.3200, 3310.0),
               (0.6000, 3570.0),
               (0.7000, 3720.0),
               (1.0000, 3720.0)],
    ),
    # Las Dominguez -- unchanged from the published figure
    "DOM": dict(
        form="linear",
        nodes=[(0.0000, 2797.2),
               (0.2500, 2947.1),
               (0.7600, 3855.1),
               (1.0000, 4189.2)],
    ),
    # Belmira -- reshaped.  The site reaches no paramo, so the summit is set a
    # little above the upper sampling line (3254 m, BEL_1) at 3320 m, and the
    # upper half of the curve is flattened so that the profile crosses the
    # 3150 m forest/subparamo boundary at about x = 0.75.  That gives the
    # subparamo belt a visible width instead of the thin corner triangle a
    # steep summit produced.
    "BEL": dict(
        form="spline",
        nodes=[(0.0000, 2491.8),
               (0.1667, 2560.0),
               (0.3333, 2690.0),
               (0.5000, 2860.0),
               (0.6667, 3055.0),
               (0.8333, 3225.0),
               (1.0000, 3320.0)],
    ),
}

# Vegetation belts, lowest first, as (label, upper bound in m).  The lowest
# belt starts at the foot of the profile; `None` means "up to the summit".
# Only the belts listed for a site are drawn and only they appear in that
# panel's legend -- Matarredonda has no forest, Belmira reaches no paramo.
BELTS = {
    "NEV": [("Andean Forest", 3400.0), ("Subparamo", 3600.0), ("Páramo", None)],
    "MAT": [("Subparamo", 3550.0), ("Páramo", None)],
    "DOM": [("Andean Forest", 3400.0), ("Subparamo", 3600.0), ("Páramo", None)],
    "BEL": [("Andean Forest", 3150.0), ("Subparamo", None)],
}

# Sampling locations, from Table S1 of the supplementary appendix
# (Analysis_pipeline.qmd).  Drawn as dashed lines running from x = 0.1 to the
# first point where they meet the mountain profile.
#
# NOTE: the published figure draws Matarredonda's upper line at ~3708 m, which
# disagrees with Table S1 (3678 m) and with Data_S1_metadata.csv (3677 m).  All
# eleven other lines match Table S1 exactly, so this was a slip in the lost
# original.  It is corrected here.
SAMPLING = {
    "NEV": [3828.0, 3550.0, 3333.0, 3078.0],
    "MAT": [3678.0, 3381.0],
    "DOM": [3815.0, 3534.0, 3234.0],
    "BEL": [3254.0, 2950.0, 2715.0],
}

# Dashed-line colour per region (measured off the published figure)
DASH_COLOR = {"NEV": "#FF0000", "MAT": "#FF0000", "DOM": "#0000FF", "BEL": "#800080"}

# Per-panel limits used only when COMMON_Y is False
YLIM_INDEPENDENT = {
    "NEV": (2300.0, 4400.0),
    "MAT": (3000.0, 3900.0),
    "DOM": (2600.0, 4400.0),
    "BEL": (2300.0, 3500.0),
}

TITLES = {"NEV": "La Nevera", "MAT": "Matarredonda",
          "DOM": "Las Domínguez", "BEL": "Belmira"}

# Panel order: top-left, top-right, bottom-left, bottom-right
PANEL_ORDER = ["NEV", "MAT", "DOM", "BEL"]

# =============================================================================
# ======================  RENDER SETTINGS  ====================================
# =============================================================================

FIG_W_IN, FIG_H_IN, DPI = 18.0, 12.0, 200      # -> 3600 x 2400 px at scale 1

# map_code.R sets FIG1_RENDER_SCALE so the whole figure can be exported at the
# print resolution requested there (OUT_DPI / OUT_WIDTH_CM).  It multiplies the
# pixel count only -- the layout fractions and the point sizes are untouched,
# so the panel looks identical, just at higher resolution.
RENDER_SCALE = float(os.environ.get("FIG1_RENDER_SCALE", "1"))

# Axes geometry, as fractions of the figure.  These are NOT matplotlib defaults
# and are not what tight_layout() produces -- they were measured off the spine
# positions of the published Fig_1_map.png (left spines at x = 454.5 / 1951.5,
# right spines at 1783 / 3280.5, axes rows spanning y = 1871.5..2970 and
# 3066..4165 within the 3600 x 2400 panel-B half of the composite).
ADJUST = dict(left=0.126250, right=0.911250, bottom=0.015417, top=0.970521,
              wspace=0.126837, hspace=0.087395)

BELT_COLORS = {                                 # full-strength; fills use ALPHA
    "Andean Forest": "#0B6623",
    "Subparamo":     "#A2D729",
    "Páramo":        "#C9DF8A",
}
BELT_ALPHA = 0.7

TERRAIN_LW   = 3.5          # pt
DASH_LW      = 1.3          # pt  (0.8 in the published figure)
DASH_X0      = 0.10         # dashed lines start here

# TEXT_SCALE is set by map_code.R (see below).  The base sizes here are
# pixel-matched to panel A's text at scale 1 (panel B is drawn on an 18 in
# canvas at 200 dpi, panel A on a 12 in canvas at 300 dpi, so a point is worth
# 2.78 px here and 4.17 px there -- the published figure already balanced the
# two, and scaling both by the same factor keeps them balanced).
# map_code.R exports FIG1_TEXT_SCALE so the two halves cannot drift apart; the
# literal below is the fallback when this script is run on its own.
TEXT_SCALE = float(os.environ.get("FIG1_TEXT_SCALE", "1.55"))
FS_TITLE   = 16.0 * TEXT_SCALE   # subplot title (no panel-A counterpart)
FS_LABEL   = 14.1 * TEXT_SCALE   # "Elevation (m)"  <-> panel A's axis titles
FS_TICK    = 10.0 * TEXT_SCALE   # tick labels      <-> panel A's tick labels
FS_LEGEND  = 13.5 * TEXT_SCALE   # legend text      <-> panel A's legend text
# labelspacing/borderpad are tighter than matplotlib's defaults: at the enlarged
# text size a default-spaced three-row legend reached down to La Nevera's 3828 m
# sampling line and clipped its left end.
LEGEND_KW = dict(loc="upper left", framealpha=1.0, handlelength=1.87,
                 handleheight=0.58, handletextpad=0.9,
                 labelspacing=0.25, borderpad=0.35)
GRID_KW = dict(axis="y", color="#e6e6e6", linewidth=0.8, zorder=0)

TAG_B = dict(x=0.0912, y=0.9815, s="B", fontsize=28, va="top", ha="left")

# =============================================================================
# ======================  MACHINERY  ==========================================
# =============================================================================

def natural_cubic(xk, yk, xq):
    """Natural cubic spline interpolation (stand-in for scipy's CubicSpline)."""
    xk, yk, xq = map(np.asarray, (xk, yk, xq))
    n = len(xk)
    if n < 3:
        return np.interp(xq, xk, yk)
    h = np.diff(xk)
    A = np.zeros((n, n))
    r = np.zeros(n)
    A[0, 0] = A[-1, -1] = 1.0                    # natural: y'' = 0 at both ends
    for i in range(1, n - 1):
        A[i, i - 1] = h[i - 1]
        A[i, i] = 2 * (h[i - 1] + h[i])
        A[i, i + 1] = h[i]
        r[i] = 6 * ((yk[i + 1] - yk[i]) / h[i] - (yk[i] - yk[i - 1]) / h[i - 1])
    m = np.linalg.solve(A, r)                    # second derivatives
    j = np.clip(np.searchsorted(xk, xq) - 1, 0, n - 2)
    a = (xk[j + 1] - xq) / h[j]
    b = 1.0 - a
    return (a * yk[j] + b * yk[j + 1]
            + ((a ** 3 - a) * m[j] + (b ** 3 - b) * m[j + 1]) * h[j] ** 2 / 6.0)


def profile_curve(site, n=1201):
    """Dense (x, elevation) polyline for one mountain."""
    spec = PROFILES[site]
    xk = np.array([p[0] for p in spec["nodes"]], float)
    yk = np.array([p[1] for p in spec["nodes"]], float)
    x = np.linspace(xk[0], xk[-1], n)
    y = natural_cubic(xk, yk, x) if spec["form"] == "spline" else np.interp(x, xk, yk)
    return x, y


def x_first_reaching(x, y, elev):
    """First x at which the profile reaches `elev` (handles flat summits)."""
    if elev <= y[0]:
        return x[0]
    if elev > y[-1]:
        return x[-1]
    i = int(np.argmax(y >= elev))                # first index at or above
    if i == 0 or y[i] == y[i - 1]:
        return float(x[i])
    t = (elev - y[i - 1]) / (y[i] - y[i - 1])
    return float(x[i - 1] + t * (x[i] - x[i - 1]))


def belt_spans(site, ymin, ymax):
    """(label, lower, upper) for each belt actually present at this site."""
    out, lo = [], ymin
    for label, top in BELTS[site]:
        hi = ymax if top is None else min(top, ymax)
        if hi > lo:
            out.append((label, lo, hi))
        lo = hi
    return out


def draw_panel(ax, site, show_ylabel=True):
    x, y = profile_curve(site)
    ymin, ymax = float(y.min()), float(y.max())

    ax.set_axisbelow(True)
    ax.grid(**GRID_KW)

    # ---- vegetation belts: the mountain interior, split at the belt bounds ---
    spans = belt_spans(site, ymin, ymax)
    for label, b0, b1 in spans:
        ax.fill_between(x, np.clip(y, b0, b1), b0,
                        where=(y > b0), interpolate=True,
                        facecolor=BELT_COLORS[label], alpha=BELT_ALPHA,
                        linewidth=0, zorder=1)

    # ---- dashed sampling elevations -----------------------------------------
    for elev in SAMPLING[site]:
        ax.plot([DASH_X0, x_first_reaching(x, y, elev)], [elev, elev],
                linestyle="--", color=DASH_COLOR[site], linewidth=DASH_LW,
                zorder=3, solid_capstyle="butt")

    # ---- mountain silhouette -------------------------------------------------
    ax.plot(x, y, color="black", linewidth=TERRAIN_LW, zorder=4,
            solid_capstyle="butt", solid_joinstyle="miter")

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(*(YLIM_COMMON if COMMON_Y else YLIM_INDEPENDENT[site]))
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=FS_TICK)
    # Only the left column carries the axis title: all four panels share the
    # same elevation axis, and at the enlarged text size the right column's
    # label overran the gap and was drawn across the left column's plot.
    if show_ylabel:
        ax.set_ylabel("Elevation (m)", fontsize=FS_LABEL)
    ax.set_title(TITLES[site], fontsize=FS_TITLE, fontweight="bold")
    handles = [Patch(facecolor=BELT_COLORS[lab], label=lab)
               for lab, _, _ in reversed(spans)]        # highest belt first
    ax.legend(handles=handles, fontsize=FS_LEGEND, **LEGEND_KW)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(FIG_W_IN, FIG_H_IN), dpi=DPI)
    for i, (ax, site) in enumerate(zip(axes.ravel(), PANEL_ORDER)):
        draw_panel(ax, site, show_ylabel=(i % 2 == 0))   # left column only
    fig.subplots_adjust(**ADJUST)
    fig.text(color="black", **TAG_B)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=DPI * RENDER_SCALE, facecolor="white")
    plt.close(fig)
    print("wrote %s  (COMMON_Y = %s, render scale %.4f)"
          % (OUT, COMMON_Y, RENDER_SCALE))


if __name__ == "__main__":
    main()
