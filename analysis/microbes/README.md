# 2016 microbial diversity

Open the repository as the working folder in Positron. The scripts require no
command-line arguments.

Run these files in order:

```r
source("analysis/microbes/01_prepare_diversity.R")
source("analysis/microbes/02_basic_diversity.R")
```

`01_prepare_diversity.R` joins the DADA2 samples to the environmental metadata,
records basic sequencing summaries, excludes libraries below 10,000 reads, and
performs one reproducible rarefaction. It writes three ordinary CSV files to
`data/derived/microbial_diversity_2016/`:

- `sample_metadata.csv` contains all 120 libraries and an `included_10k` column;
- `asv_counts_unrarefied.csv` contains the 102 eligible libraries and is used
  only to calculate expected iNEXT diversity at 10,000 reads;
- `asv_counts_10k.csv` is the primary community table for all ecological
  analyses.

The count tables place ASVs in rows and samples in columns so they can be
opened in spreadsheet software. No RDS files are created.

`02_basic_diversity.R` reads those CSVs and calculates alpha and gamma Hill
numbers, Hellinger PCA, and a simple habitat RDA using `vegan` and `iNEXT`.
All ecological analyses use the 10K table. Tables are written to
`results/diversity_2016/tables/`, and figures are written to `figures/`.
