# Result guide

Compact analysis tables, input audits, and provenance are retained for
reproducibility. The root `.gitignore` names specific caches, pilot outputs,
and duplicate package exports; it does not exclude this directory wholesale.

| Directory | Role and current status |
| --- | --- |
| `diversity_2016/soil_porewater_sediment/` | Current reviewed 76-library comparative analysis: habitat diversity/composition, soil inventories, paired comparisons, and porewater checks; [methods](../analysis/microbes/README_SOIL_POREWATER_SEDIMENT.md) |
| `fticr_2016/matched_sediment_dispersion/` | Matched 33-site microbial-molecular dispersion, shared selections, and distance matrices; [methods](../analysis/fticr/README_MATCHED_SEDIMENT_DISPERSION.md) |
| `fticr_2016/sediment_characteristics/` | Eligible-sample sediment-condition models and paired property adjustment; [methods](../analysis/fticr/README_SEDIMENT_CHARACTERISTICS.md) |
| `fticr_2016/bipartite_networks/`, `fticr_2016/stable_group_distribution/` | Exploratory associations, explicit nulls, selected-group identities, and maps; [methods](../analysis/fticr/README_BIPARTITE_NETWORKS.md) |
| Other `fticr_2016/` subdirectories | Molecular and measured-chemistry analyses; [FT-ICR guide](../analysis/fticr/README.md). Water-community models predate the cohort exclusions. |
| `sediment_manuscript_figures_2016/` | Plot data and checksums for reformatted within-sediment figures; [figure report](../analysis/figures/README_SEDIMENT_FIGURES.md) |
| Other `diversity_2016/` subdirectories | Baseline diversity, UniFrac, network, and soil-routing/sharing results. Four-habitat and water outputs are historical; sediment-only values retain their sample set. |
| `dada2_2016/` | Completed inference, counts, taxonomy, QC, and provenance; preserve as source products |
| `phylogeny_2016/` | Completed master tree and provenance; downstream scripts prune it in memory |
| `icamp_2016/catchment_castor/` | Completed four-habitat analysis including suspect libraries; historical for the current manuscript; [results](../analysis/microbes/README_ICAMP_RESULTS.md) |
| `script_clarity_2016/` | Historical readability-rewrite validation and timings |

Sediment iCAMP outputs belong in `icamp_2016/sediment_castor/` when obtained
from submitted job 432228; its current state has not been checked. The reviewed
three-habitat iCAMP pool has not been run. Do not infer either status from the
existence of a directory or distance descriptor.

The new comparative directory saves model designs, seeds, sample exclusions,
null statistics, sensitivity selections, hashes, and session information.
Use its tables with its own cohort, rather than combining them with unfiltered
100-library summaries. No results were moved, deleted, or recalculated during
the October 10 repository tidy-up.
