# 2016 microbial diversity

Open the repository as the working folder in Positron. The scripts require no
command-line arguments.

Run these files in order:

```r
source("analysis/microbes/01_prepare_diversity.R")
source("analysis/microbes/02_basic_diversity.R")
source("analysis/microbes/03_build_phylogeny.R")
source("analysis/microbes/04_network_environment.R")
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

`09_network_dispersion.R` compares headwater, intermediate, and mainstem
microbial beta diversity separately for sediment, hyporheic, and planktonic
communities. It uses the fixed 10K table, Bray-Curtis dispersion, a Jaccard
incidence sensitivity, and segment-balanced resampling. It runs directly
with `source(here::here("analysis", "microbes", "09_network_dispersion.R"))`.
See [README_NETWORK_DISPERSION.md](README_NETWORK_DISPERSION.md) for methods,
sample counts, results, figure captions, and the deferred beta-NTI design.
