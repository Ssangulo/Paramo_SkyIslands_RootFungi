#!/usr/bin/env bash
# =============================================================================
# run_revision.sh   (final)
#
# Runs every revision analysis in dependency order. From the repository root:
#   bash Scripts/run_revision.sh
#
# Console output for each script is saved to logs/<script>.log.
#
# Settings worth knowing:
#   - Per-plant diversity models use a gamma GLM with a log link (FAMILY at the
#     top of perplant_hill_size_standardised.R and parallelism_diagnostics.R;
#     set both to "lognormal" to reproduce the earlier version).
#   - Model structure: interaction retained only where AICc supports it, which
#     is abundance-weighted PD at q = 1 and q = 2. Everywhere else additive.
#   - Figure 4 uses one across-site line (additive model); the superseded
#     site-slopes version is in archive/.
# =============================================================================
set -euo pipefail

if [ ! -f objects/ps_individual.rds ] || [ ! -f objects/tree_ps.rds ]; then
  echo "Run this from the repository root (objects/ps_individual.rds and objects/tree_ps.rds not found)."
  exit 1
fi
mkdir -p logs tables figures objects

run () {
  echo "=== $(date '+%H:%M:%S')  $1"
  Rscript "Scripts/$1" 2>&1 | tee "logs/${1%.R}.log"
}

# --- Diversity ----------------------------------------------------------------
run iNEXT3D_plant_level.R              # slowest (iNEXT3D bootstraps); writes the Figure 3A panel
run perplant_hill_size_standardised.R  # rarefied per-plant Hill numbers, tests, Figure 3B and full Figure 3
run parallelism_diagnostics.R          # needs objects/perplant_hill_size_std.rds from the previous step

# --- Taxa ---------------------------------------------------------------------
run Figure4_pooled_glmm.R              # genus GLMMs + main-text Figure 4 (one script)
run helotiales_fraction_check.R

# --- Composition and phylogenetic structure -----------------------------------
run adjusted_tests_centroid_nri.R      # slow unless objects/nri_nti_per_plant_999.rds exists
run reproduce_composition_tables.R

# --- Model selection and sensitivity ------------------------------------------
run model_selection_aicc.R
run family_sensitivity.R

echo "=== $(date '+%H:%M:%S')  All done. Logs are in logs/."
