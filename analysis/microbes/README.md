# 2016 microbial diversity

Soil-to-stream sharing and drainage localization are implemented in
`10_soil_drainage_routing.py`, `11_soil_stream_localization.R`, and
`12_soil_branch_signatures.R`. See
[the soil drainage report](README_SOIL_STREAM_LOCALIZATION.md) for definitions,
GPS exclusions, public terrain provenance, results, and figure captions.
`13_soil_stream_distance.R` provides an independent distance-only nearness
analysis requiring no terrain or drainage-membership inputs.
`14_soil_asv_catchment_figure.R` makes one combined-habitat plot of drainage
area against soil-detected ASV read percentage, with an optional richness
fraction setting. Separate descriptive GAM curves show habitat-specific
patterns across drainage area. `15_first_detection_flowpaths.R` adapts the
Ruiz-Gonzalez upstream-first detection inventory to Soil, Sediment, Aquatic
(combined hyporheic/planktonic detections), and Sediment + Aquatic ties across
headwater, intermediate, and mainstem groups. Receiving habitat panels remain
separate, with a 1,000-draw equal-sample sensitivity. The flowpath figure is
recorded as a growing descriptive result; retain both rows while the
manuscript choice of richness or read percentages remains open.

ASV example and distance-correlation plots are disabled and archived under
`figures/archive/soil_first_pass/`. The initial branch plots are also archived
pending refinement; the Python drainage backbone and analysis tables remain.

Open the repository as the working folder in Positron. The scripts require no
command-line arguments.

Run the relevant files in order; the completed phylogeny need not be rebuilt.
The revised iCAMP scripts are ready for a cluster pilot; obtain approval
before transferring files or submitting jobs. See the rerun plan below.

```r
source("analysis/microbes/01_prepare_diversity.R")
source("analysis/microbes/02_basic_diversity.R")
source("analysis/microbes/03_build_phylogeny.R")
source("analysis/microbes/04_network_environment.R")
source("analysis/microbes/05_unifrac.R")
# Submit analysis/cluster/run_icamp_sediment.sh for the long iCAMP run.
# Submit analysis/cluster/run_icamp_catchment.sh for the catchment run.
```

`01_prepare_diversity.R` joins the DADA2 samples to the environmental metadata,
retains assigned bacterial and archaeal ASVs, removes chloroplast and
mitochondrial sequences, records basic sequencing summaries, excludes libraries
below 10,000 target reads, and performs one reproducible rarefaction. It writes
four ordinary CSV files to
`data/derived/microbial_diversity_2016/`:

- `sample_metadata.csv` contains all 120 libraries and an `included_10k` column;
- `asv_metadata_screened.csv` contains the full bacterial/archaeal ASV catalog,
  its sequences, abundances, lengths, and taxonomy;
- `asv_counts_unrarefied.csv` contains the eligible libraries and is used
  only to calculate expected iNEXT diversity at 10,000 reads;
- `asv_counts_10k.csv` is the primary community table for all ecological
  analyses.

The count tables place ASVs in rows and samples in columns so they can be
opened in spreadsheet software. No RDS files are created.

`02_basic_diversity.R` reads those CSVs and calculates alpha and gamma Hill
numbers, Hellinger PCA, and a simple habitat RDA using `vegan` and `iNEXT`.
All ecological analyses use the 10K table. Tables are written to
`results/diversity_2016/tables/`, and figures are written to `figures/`.

`03_build_phylogeny.R` aligns all ASVs in the full screened catalog with
`DECIPHER::AlignSeqs` and estimates one unrooted nucleotide tree with FastTree.
Downstream analyses should prune this master tree to the ASVs in the community
table being analyzed; a different rarefaction therefore does not require a new
alignment or tree. FASTA, Newick, log, and summary files are written to
`results/phylogeny_2016/`.

Build the full tree on the cluster:

```sh
mkdir -p logs
sbatch analysis/cluster/run_phylogeny.sh
```

The job loads `openmpi4` and R 4.4.0, uses the OpenMP FastTree installation at
`~/scratch/shared/bioinformatics/FastTree/FastTree`, and requests 8 CPUs,
128 GB RAM, and 96 hours.

`04_network_environment.R` is deliberately exploratory. It joins Hill diversity
to network position, makes separate Hellinger PCAs within each aquatic habitat,
and plots alpha diversity and community scores along drainage area, distance to
the outlet, planar coordinates, and sediment EEA. It fits no permutation tests
or environmental models. Inspect these figures before choosing a small number
of models in the next analysis step.

`05_unifrac.R` prunes the completed master tree to the ASVs in the fixed 10K
table, roots it at its midpoint in memory, and calculates unweighted and
weighted UniFrac with `phyloseq`. It writes the distance matrices and PCoA
scores to `results/diversity_2016/tables/` and the ordination figure to
`figures/`. The full tree on disk is not changed and is retained for later
iCAMP analyses. Phylogenetic Hill diversity is not calculated because the
current implementation is impractically slow for this 75K-tip microbial tree.

The iCAMP analysis uses two scales: all four habitats together to test
catchment-wide habitat assembly, and sediment alone for the detailed microbial,
FT-ICR-MS, and sediment-environment integration. An aquatic-only intermediate
run is intentionally omitted.

`06_icamp_sediment.R` runs the sediment analysis on the 34 sediment libraries
retained at 10K reads. It uses the pruned sediment tree, 1,000 null
randomizations, and writes pairwise and bin-level process tables to
`results/icamp_2016/sediment_castor/`. Run it on the cluster with:

```sh
mkdir -p logs
sbatch analysis/cluster/run_icamp_sediment.sh
```

`07_icamp_catchment.R` uses all 100 libraries in the fixed 10K table and one
observed catchment-wide metacommunity pool. It uses the same iCAMP settings as
the sediment run, prunes the tree to the 49,260 observed ASVs, and writes
results to `results/icamp_2016/catchment_castor/`. The pairwise CSV includes habitat
labels; `process_importance_by_habitat_pair.csv` gives descriptive mean process
fractions and pair counts for each within- or across-habitat comparison.

For the final run, submit from the repository root after the approved pilot:

```sh
mkdir -p logs
sbatch analysis/cluster/run_icamp_catchment.sh
```

The existing catchment job requests 32 CPUs, 256 GB RAM, and 96 hours; the
replacement walltime will be selected after the distance/pilot benchmarks.
Both analyses
retain a detailed RDS checkpoint alongside inspectable CSV results. The
phylogenetic distance matrices are large computational intermediates. Both
scripts calculate tip distances with castor, validate their disk-backed copy,
and pass the completed matrix explicitly to iCAMP. A completion record and
input hashes prevent silent reuse of partial or outdated matrices. No ASVs
or samples are removed to speed up the distance calculation. Set `test_run`
near the top of the same R script for a 100-draw pilot with separate outputs;
leave it `FALSE` for the final 1,000-draw analysis. No command-line arguments
or separate resume script are required.

See `analysis/README_SEDIMENT_INTEGRATION.md` for the sediment FT-ICR-MS,
environmental, enzyme, and source-context analysis plan.

`08_sediment_figures.R` summarizes the 34 fixed-10K sediment samples with
Hellinger PCA and Hill q = 1 figures.

`09_network_dispersion.R` compares headwater, intermediate, and mainstem
microbial beta diversity separately for sediment, hyporheic, and planktonic
communities. It uses the fixed 10K table, Bray-Curtis dispersion, a Jaccard
incidence sensitivity, and segment-balanced resampling. It runs directly
with `source(here::here("analysis", "microbes", "09_network_dispersion.R"))`.
See [README_NETWORK_DISPERSION.md](README_NETWORK_DISPERSION.md) for methods,
sample counts, results, figure captions, and the deferred beta-NTI design.

Catchment iCAMP job 431635 timed out during distance construction. See
[the rerun plan](../README_ICAMP_RERUN.md) before submitting another job;
the existing distance matrix is incomplete and must not be reused. The revised
scripts use fresh `_castor` output folders and leave the original intact.
