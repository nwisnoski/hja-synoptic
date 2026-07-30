# Loading analysis data

Use one entry point for later analysis scripts:

```r
source("analysis/00_load_analysis_data.R")
hja <- load_hja_analysis_data()
hja
```

The loader is read-only. It verifies raw-source MD5 checksums, unique
identifiers, the 60-site environment/FT-ICR row order, the primary feature
filter, and—after DADA2 finishes—the 44-site multiblock alignment.

## Main objects

- `hja$environment$fticr_60_sites`
- `hja$environment$sediment_44_sites`
- `hja$environment$sediment_44_complete_primary_predictors`
- `hja$environment$sediment_44_incomplete_primary_predictors`
- `hja$environment$sediment_44_sensitivity_predictors`
- `hja$environment$predictor_blocks`
- `hja$fticr$intensity_60_primary_features`
- `hja$fticr$peak_metadata`
- `hja$fticr$laboratory_standards`

When `hja$dada2_ready` is `TRUE`, it also contains:

- `hja$microbes$asv_counts_all`
- `hja$microbes$sample_metadata`
- `hja$sediment_multiblock`
- `hja$source_pools`

## Explicit in-memory views

With the default `prepare_views = TRUE`, each appropriate data block has
clearly named:

- `raw`
- `row_relative`
- `hellinger`
- `presence_absence`

These do not replace or modify raw matrices. The formulas and thresholds are
listed in `hja$transformation_log`. Any zero-total microbial rows remain loaded,
are listed in the view's `zero_total_rows`, and receive `NA` relative/Hellinger
values rather than being dropped. No environmental log transformation,
standardization, imputation, rarefaction, feature pruning, or statistical
filtering is performed by the loader.

Set `prepare_views = FALSE` to load only raw data. Set
`require_dada2 = TRUE` when running a microbial analysis that should fail
instead of returning a partial FT-ICR/environment-only object.
