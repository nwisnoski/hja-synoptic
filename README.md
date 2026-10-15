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
The [sampling map](figures/supplement/habitat_diversity/2016_soil_porewater_sediment_sampling_map.pdf),
individual panels, the fixed genus-classified sensitivity, additional same-site
metrics, and the full 44-site molecular class figure are grouped by analysis in
[figures/supplement](figures/supplement/README.md).

For the full exploratory iCAMP context, start with the five overview figures in
[reviewed_pool](figures/supplement/icamp/README.md), produced by
[the reviewed-pool analysis](analysis/supporting/icamp/01_reviewed_pool.R): sediment
and porewater composition, their maps, and the shared-site habitat comparison.
The sediment-only run is retained as a regional-pool sensitivity in
`figures/supplement/icamp/sediment_sensitivity/`. Each run has separate
`panels/` and `diagnostics/` folders. Completed assignments are reused for plotting.

## Where things belong

- `analysis/preparation/`: sample mapping, completed sequence processing,
  microbial tables, phylogeny, and UniFrac. See its [preparation guide](analysis/preparation/README.md)
  only when inputs need rebuilding.
- `analysis/supporting/`: numbered molecular, chemistry, network, and sediment
  analyses. These supply the main workflow and additional figures. Its `icamp/`
  subfolder groups the reviewed-pool analysis and sediment-only sensitivity.
- `analysis/cluster/`: batch launchers for sequence processing, phylogeny, and iCAMP.
- `figures/manuscript/`: the four selected figures. `figures/supplement/`:
  panels and supporting analyses in named subfolders; its guide lists entry points.
- `data/`: original inputs; `data/derived/`: constructed inputs.
  `results/`: analysis tables, saved selections, checksums, and checkpoints.
- `archive/`: legacy files and their recovery manifest. Older workflows and
  removed duplicate figures are also recoverable from Git history.

Results retain completed analyses and their required inputs. Retired code and
figures are recoverable from Git history; [cleanup manifests](archive/README.md)
record moved and removed paths. Raw inputs, sequence outputs, and inference
checkpoints are preserved.

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
phylogenetic distance cache remains on the cluster. The local inference script is
`analysis/supporting/icamp/02_sediment_sensitivity.R`, with launcher
`analysis/cluster/run_icamp_sediment.sh`. The script content is unchanged from
the completed deployment; the local launcher uses its new path. This cleanup
does not synchronize the remote checkout.
The completed four-habitat catchment run is historical and includes excluded
libraries. It does not represent the reviewed cohort.

The reviewed soil–porewater–sediment iCAMP run is also complete; its result,
input hashes, and run manifest are in `results/icamp_2016/reviewed_pool_castor/`.
[The iCAMP result guide](results/icamp_2016/README.md) distinguishes this run,
the sediment-only sensitivity, and the historical four-habitat run.
The reviewed-pool script validates the saved run identity before regenerating plots.
If its completed result is absent, it can start expensive inference; the
[cluster launcher](analysis/cluster/run_icamp_reviewed_pool.sh) requests 16 CPUs,
192 GB, and a seven-day walltime limit. Future-run runtime is unknown.

Do not submit, check, transfer, or change remote work without explicit approval
of that action. Obtain approval of exact paths and message before committing;
pushes require separate authorization.
