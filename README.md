# HJA synoptic: 2016 microbial and molecular ecology

Soil, hyporheic porewater, and sediment communities across the HJ Andrews
stream network. Start with the six figures in [figures/manuscript](figures/manuscript).

## Figure workflow

The numbered scripts in `analysis/` contain the analyses and their plotting
code. Open them directly in Positron; paths use `here()` and need no arguments.
Run `01_prepare_cohort.R` first, then the figure scripts in order.

| Script | Figure | Main output |
| --- | --- | --- |
| [01_prepare_cohort.R](analysis/01_prepare_cohort.R) | Reviewed sample set and exclusion audit | `data/derived/soil_porewater_sediment_2016/` |
| [02_habitat_diversity.R](analysis/02_habitat_diversity.R) | 1. Hill diversity (q = 1); 2. weighted UniFrac with habitat ellipses | [Alpha](figures/manuscript/2016_manuscript_alpha_diversity.pdf), [ordination](figures/manuscript/2016_manuscript_weighted_unifrac.pdf) |
| [03_first_detection.R](analysis/03_first_detection.R) | 3. First-detected inventory along the network | [Flow figure](figures/manuscript/2016_first_detection_flowpaths.pdf) |
| [04_sediment_convergence.R](analysis/04_sediment_convergence.R) | 4. Microbial and molecular dispersion on matched sites | [Convergence](figures/manuscript/2016_matched_sediment_dispersion_comparison.pdf) |
| [05_molecular_classes.R](analysis/05_molecular_classes.R) | 5. Site-level molecular classes from headwaters to mainstem | [Class composition](figures/manuscript/2016_fticr_network_group_class_comparison.pdf) |
| [06_network_localization.R](analysis/06_network_localization.R) | 6. Exploratory microbial–molecular group distribution | [Catchment map](figures/manuscript/2016_sediment_stable_group_colocation_map.pdf) |

The checkout contains the prepared inputs and completed results used by these
scripts. Script 02 uses saved UniFrac distances; 04 refits dispersion and its
balanced-sampling sensitivity; 05–06 read completed molecular/network results.
They do not rerun DADA2, phylogeny construction, or network null simulations.
Each script declares its packages and analytical settings at the top.

## Where things belong

- `analysis/preparation/`: sample mapping, completed sequence processing,
  microbial tables, phylogeny, and UniFrac. See its [preparation guide](analysis/preparation/README.md)
  only when inputs need rebuilding.
- `analysis/supporting/`: numbered molecular, chemistry, network, and sediment
  analyses. These supply the main workflow and additional figures.
- `figures/manuscript/`: the six selected figures. `figures/supplement/`:
  individual panels, diagnostics, and supporting analyses.
- `data/`: original inputs; `data/derived/`: constructed inputs.
  `results/`: analysis tables, saved selections, checksums, and checkpoints.
- `archive/`: legacy files and their recovery manifest. Older workflows and
  removed duplicate figures are also recoverable from Git history.

## Sample set and interpretation

The reviewed fixed-10K cohort has **76 samples: 15 soil, 27 porewater, and
34 sediment**. Script 01 records the 15 suspect water-library exclusions
individually, preserves the original rarefaction, and separates QC exclusions
from the decision to omit plankton. The cause of the anomalies is unresolved.
There are 18 sediment–porewater pairs and 33 sediment microbial–molecular
matches; sediment site 43 has no molecular profile. Soils form a regional pool.

[Methods](analysis/METHODS.md) gives the definitions and interpretation limits
needed for the figures. Numerical results and diagnostics belong in the saved
tables rather than duplicated reports.

## Cluster work

Sediment iCAMP job **432228** was submitted October 10, 2026; its current state
has not been checked. Its deployed script remains at
`analysis/microbes/06_icamp_sediment.R`, with launcher
`analysis/cluster/run_icamp_sediment.sh`. These are intentionally unchanged.
The completed four-habitat catchment run is historical and includes excluded
libraries. It does not represent the reviewed cohort.

Do not submit, check, transfer, or change remote work without explicit approval
of that action. Obtain approval of exact paths and message before committing;
pushes require separate authorization.
