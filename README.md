# HJA synoptic: 2016 microbial and molecular ecology

Soil, hyporheic porewater, and sediment communities across the HJ Andrews
stream network. Start with the four figures in [figures/manuscript](figures/manuscript).

## Figure workflow

The numbered scripts in `analysis/` contain the analyses and their plotting
code. Open them directly in Positron; paths use `here()` and need no arguments.
Run `01_prepare_cohort.R` first, then the figure scripts in order.

| Script | Figure | Main output |
| --- | --- | --- |
| [01_prepare_cohort.R](analysis/01_prepare_cohort.R) | Reviewed sample set and exclusion audit | `data/derived/soil_porewater_sediment_2016/` |
| [02_habitat_diversity.R](analysis/02_habitat_diversity.R) | 1. Hill diversity (q = 1) and weighted UniFrac | [Habitat diversity](figures/manuscript/2016_habitat_diversity.pdf) |
| [03_first_detection.R](analysis/03_first_detection.R) | 2. First detection, paired soil sharing, and taxonomic resolution | [Paired comparisons](figures/manuscript/2016_first_detection_flowpaths.pdf) |
| [04_sediment_convergence.R](analysis/04_sediment_convergence.R) | 3. Dispersion and molecular classes on the same 33 sediment sites | [Dispersion and classes](figures/manuscript/2016_matched_sediment_dispersion_comparison.pdf) |
| [05_network_localization.R](analysis/05_network_localization.R) | 4. Exploratory microbial–molecular group distribution | [Catchment map](figures/manuscript/2016_sediment_stable_group_colocation_map.pdf) |

The checkout contains the prepared inputs and completed results used by these
scripts. Script 02 uses saved UniFrac distances and supplies the paired tests
for 03; 04 refits dispersion and its balanced-sampling sensitivity, then reads
completed molecular-class results; 05 reads completed network results.
They do not rerun DADA2, phylogeny construction, or network null simulations.
Each script declares its packages and analytical settings at the top.
Microbial abundance comparisons use Hellinger distance. Script 03 compares
same-site resemblance across taxonomic ranks and audits classification coverage.
The [sampling map](figures/supplement/2016_soil_porewater_sediment_sampling_map.pdf),
individual panels, the fixed genus-classified sensitivity, additional same-site
metrics, and the full 44-site molecular class figure remain in `figures/supplement/`.

[Supporting script 12](analysis/supporting/12_sediment_icamp_figures.R) maps the
completed sediment iCAMP assignments and compares headwater, intermediate, and
mainstem turnover. Start with the [network composition](figures/supplement/2016_sediment_icamp_network_composition.pdf)
and [catchment maps](figures/supplement/2016_sediment_icamp_catchment_maps.pdf);
the [balanced comparison](figures/supplement/2016_sediment_icamp_network_comparison.pdf)
shows sampling sensitivity. It reads saved results without rerunning iCAMP.
The [all-pair distributions](figures/supplement/2016_sediment_icamp_all_within_pairs.pdf)
show every within-position sample pair, including same-segment pairs, with
equal weight; their [composition view](figures/supplement/2016_sediment_icamp_all_within_pair_composition.pdf)
shows arithmetic means.

## Where things belong

- `analysis/preparation/`: sample mapping, completed sequence processing,
  microbial tables, phylogeny, and UniFrac. See its [preparation guide](analysis/preparation/README.md)
  only when inputs need rebuilding.
- `analysis/supporting/`: numbered molecular, chemistry, network, and sediment
  analyses. These supply the main workflow and additional figures.
- `figures/manuscript/`: the four selected figures. `figures/supplement/`:
  individual panels, diagnostics, and supporting analyses.
- `data/`: original inputs; `data/derived/`: constructed inputs.
  `results/`: analysis tables, saved selections, checksums, and checkpoints.
- `archive/`: legacy files and their recovery manifest. Older workflows and
  removed duplicate figures are also recoverable from Git history.

Results retain current analyses and their required inputs. Retired workflows
and redundant plot-data, settings, and summary exports have been removed;
their earlier versions remain in Git history.

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

Sediment iCAMP job **432228** completed October 11, 2026, with exit code `0:0`.
Its small results are saved in `results/icamp_2016/sediment_castor/`; the
phylogenetic distance cache remains on the cluster. Its deployed script remains at
`analysis/microbes/06_icamp_sediment.R`, with launcher
`analysis/cluster/run_icamp_sediment.sh`. These are intentionally unchanged.
The completed four-habitat catchment run is historical and includes excluded
libraries. It does not represent the reviewed cohort.

[Supporting script 13](analysis/supporting/13_reviewed_pool_icamp.R) defines a
fresh iCAMP sensitivity using the reviewed soil, porewater, and sediment pool.
It is built but has not been run. The [cluster launcher](analysis/cluster/run_icamp_reviewed_pool.sh)
requests 16 CPUs, 192 GB, and a seven-day walltime limit; runtime is unknown.
Results will go in `results/icamp_2016/reviewed_pool_castor/`. The same script
will make separate sediment and porewater maps, network-position comparisons,
and all-within-position pair distributions, plus an 18-shared-site comparison.
Reopening a completed run reads its saved result and regenerates the figures.

Do not submit, check, transfer, or change remote work without explicit approval
of that action. Obtain approval of exact paths and message before committing;
pushes require separate authorization.
