# 2016 microbial diversity

Run the data alignment, then the basic diversity analysis:

```sh
Rscript analysis/microbes/01_prepare_diversity.R
Rscript analysis/microbes/02_basic_diversity.R
```

`02_basic_diversity.R` has one explicit branch. The raw ASV table is used only
for sequencing-depth and detection summaries, plus a breakaway sensitivity
table. Samples below 10,000 reads are then excluded, and one fixed-seed 10K
count table is saved as `results/diversity_2016/derived/diversity_10k.rds`.
All ecological analyses use that saved table.

The script uses established package functions:

- `vegan::specnumber` and `vegan::diversity` for alpha Hill numbers q = 0, 1,
  and 2;
- `iNEXT::estimateD` for expected alpha Hill numbers at exactly 10,000 reads,
  avoiding variation from a random rarefaction draw;
- `iNEXT::estimateD` for incidence-based gamma Hill numbers at common sample
  coverage;
- `vegan::decostand(..., "hellinger")` and `vegan::rda` for Hellinger PCA and
  the initial habitat RDA.

All three diversity components use one count table rarefied to 10,000 reads.
This retains 102 of 120 libraries. Phylogenetic diversity and iCAMP are later
analyses and are intentionally not included here.

Tables go to `results/diversity_2016/tables/`; figures go to `figures/`.
Breakaway remains a raw-data sensitivity analysis and is not used as the
downstream ecological response.
