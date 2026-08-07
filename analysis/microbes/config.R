# Shared choices for the 2016 microbial-diversity workflow.

microbe_config <- list(
  input = list(
    dada2 = "results/dada2_2016",
    analysis_inputs = "data/derived/analysis_inputs"
  ),
  output = "results/diversity_2016",
  figures = "figures",
  rarefaction_depth = 10000L,
  hill_orders = c(0, 1, 2),
  hill_bootstraps = 50L,
  permutations = 999L,
  random_seed = 2016L,
  habitat_order = c("planktonic", "hyporheic", "sediment", "soil"),
  habitat_labels = c(
    planktonic = "Planktonic streamwater",
    hyporheic = "Hyporheic porewater",
    sediment = "Stream sediment",
    soil = "Terrestrial soil"
  ),
  habitat_colors = c(
    planktonic = "#2C7FB8",
    hyporheic = "#41B6C4",
    sediment = "#8C6D31",
    soil = "#4D9221"
  )
)
