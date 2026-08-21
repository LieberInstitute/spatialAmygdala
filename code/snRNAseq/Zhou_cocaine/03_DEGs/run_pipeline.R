#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R — Execute the full enrichment analysis pipeline
# =============================================================================
#
# BEFORE RUNNING:
#   1. Edit 00_config.R with your file paths and domain column name
#   2. Place the DEG file and spatial atlas .rds in the working directory
#      (or update paths in 00_config.R)
#
# USAGE:
#   Rscript run_pipeline.R          # run all steps
#   Rscript run_pipeline.R 3 4      # run only steps 3 and 4
#
# PIPELINE:
#   01  Process DEGs + ortholog mapping     (needs: DEG file)
#   02  Load atlas + score spots            (needs: atlas .rds, step 01)
#   03  Approach 1: linear models           (needs: step 02)
#   04  Approach 2: Fisher tests            (needs: step 02)
#   05  Spatial visualizations              (needs: step 02)
#   06  Convergence analysis                (needs: steps 03 + 04)
# =============================================================================

scripts <- c(
  "01_degs_and_orthologs.R",
  "02_score_atlas.R",
  "03_approach1_linear_model.R",
  "04_approach2_fisher.R",
  "05_spatial_plots.R",
  "06_convergence.R"
)

# Parse command line: which steps to run
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  steps <- as.integer(args)
  scripts <- scripts[steps]
  message("Running selected steps: ", paste(steps, collapse = ", "))
} else {
  message("Running full pipeline (steps 1-6)")
}

# Execute
for (s in scripts) {
  message("\n", strrep("═", 60))
  message("  Running: ", s)
  message(strrep("═", 60), "\n")
  source(s)
}

message("\n", strrep("═", 60))
message("  PIPELINE COMPLETE")
message(strrep("═", 60))
message("Results directory: results/")
message("Key outputs:")
message("  results/approach1_lm_results.csv     — composition-corrected domain tests")
message("  results/approach1_heatmap.pdf         — adjusted means heatmap")
message("  results/approach2_fisher_results.csv  — marker overlap tests")
message("  results/approach2_heatmap.pdf         — Fisher test heatmap")
message("  results/convergence_summary.csv       — domains significant in both")
message("  results/spatial_deg_scores.pdf        — spatial score maps")
