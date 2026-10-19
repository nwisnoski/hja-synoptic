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
| [03_first_detection.R](analysis/03_first_detection.R) | 2. First detection and matched-habitat iCAMP contrasts | [Detection and assembly](figures/manuscript/2016_first_detection_flowpaths.pdf) |
| [04_sediment_convergence.R](analysis/04_sediment_convergence.R) | 3. Molecular classes and dispersion on the same 33 sediment sites | [Classes and dispersion](figures/manuscript/2016_matched_sediment_dispersion_comparison.pdf) |
| [05_network_localization.R](analysis/05_network_localization.R) | 4. Exploratory microbial–molecular group distribution | [Catchment map](figures/manuscript/2016_sediment_stable_group_colocation_map.pdf) |

The checkout contains the prepared inputs and completed results used by these
scripts. Script 02 uses saved UniFrac distances and supplies the paired tests
for 03; 04 refits dispersion and its balanced-sampling sensitivity, then reads
completed molecular-class results; 05 reads completed network results.
They do not rerun DADA2, phylogeny construction, or network null simulations.
Each script declares its packages and analytical settings at the top.
Microbial abundance comparisons use Hellinger distance. Figure 2 B–C compare
dispersal limitation and homogeneous selection for sediment and porewater on
the same 18 sites, reusing completed reviewed-pool iCAMP assignments.
Paired soil sharing, same-site taxonomic resolution, and classification coverage
remain supporting analyses in script 03.
The [supporting figure guide](figures/supplement/README.md) points to the sampling
map, individual figure panels, and the retained bulk transformation summaries.
The [methods guide](analysis/METHODS.md) defines the measurements and inference
limits; numerical results remain in the saved tables.

## Bulk molecular transformations

Run [supporting 08](analysis/supporting/08_molecular_transformations.R), then
[supporting 09](analysis/supporting/09_transformation_summary.R), using the pinned
[Danczak reference](reference/danczak_2020/README.md). These describe the catchment
inventory of putative chemical relationships. Start with the
[readable assignment summary](figures/supplement/transformations/2016_transformation_assignment_summary.pdf)
and [label-rank distribution](figures/supplement/transformations/2016_transformation_label_rank_distribution.pdf).
Pair counts, formula checks, and site coverage are in
`results/fticr_2016/transformations/`; they do not measure realized reaction rates.

## Where things belong

- `analysis/preparation/`: sample mapping, completed sequence processing,
  microbial tables, phylogeny, and UniFrac. See its [preparation guide](analysis/preparation/README.md)
  only when inputs need rebuilding.
- `analysis/supporting/`: molecular inputs and dispersion (01–03), network
  inference and stability (04–07), and bulk transformations (08–09). The `icamp/`
  subfolder retains its separate 01–03 sequence for completed-run provenance and
  plotting.
- `analysis/cluster/`: batch launchers for sequence processing, phylogeny, and iCAMP.
- `figures/manuscript/`: the four selected figures. `figures/supplement/`:
  panels and supporting analyses in named subfolders; its guide lists entry points.
- `data/`: original inputs; `data/derived/`: constructed inputs.
  `results/`: analysis tables, saved selections, checksums, and checkpoints.
- `archive/`: legacy files and their recovery manifest. Older workflows and
  removed duplicate figures are also recoverable from Git history.

Retired chemistry comparisons, sediment side analyses, enzyme associations, and
transformation gradients are packaged in [the cleanup archive](archive/README.md),
including their results and failed-run diagnostics. Raw inputs, completed
sequence processing, the master phylogeny, and inference checkpoints are intact.

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

## Completed iCAMP results

[The iCAMP result guide](results/icamp_2016/README.md) distinguishes the completed
reviewed-pool run used by Figure 2, the sediment-only sensitivity, and the
historical four-habitat run that includes excluded libraries. Reuse completed
assignments when plotting. The reviewed-pool script can start expensive inference
if its completed result is absent; do not rerun inference for this cleanup.

Remote work requires approval of each exact action; this cleanup is local.
Obtain approval of exact paths and message before committing; pushes require
separate authorization.
