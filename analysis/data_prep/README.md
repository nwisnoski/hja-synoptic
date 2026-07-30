# Reproducible preparation of HJA analysis tables

This is the R-only preparation layer between immutable source files and later
ecological analyses. It never edits, renames, reorders, or saves changes into a
source CSV or workbook.

## Run

From the repository root:

```r
system2(
  file.path(R.home("bin"), "Rscript"),
  "analysis/data_prep/run_data_preparation.R"
)
```

The runner executes inventory/crosswalk, environmental preparation, and
FT-ICR-MS preparation. It executes the sediment multiblock and source-pool
steps only after `results/dada2_2016/session_info.txt` confirms that DADA2
finished. Rerun the same command then; partial DADA2 output is never consumed.

The only preparation dependency beyond base R is `readxl`.

## Written record

- `config.R` declares every path, source column, derived name, source sample
  type, unit, suggested future transformation, and FT-ICR-MS filter.
- `01_inventory_and_crosswalk.R` records source MD5 checksums, all 76 workbook
  column positions, exact site joins, and unmatched records.
- `02_prepare_environment.R` reads declared variables from Stream, Hyporheic,
  or Sediment master-table rows and reports missingness.
- `03_prepare_fticr.R` keeps raw intensities, separates the two lab standards,
  and records every feature-level inclusion decision.
- `04_build_sediment_multiblock.R` aligns the 44 sediment ASV, FT-ICR-MS, and
  environmental rows without transforming them.
- `05_build_source_pool_tables.R` prepares same-site water/hyporheic comparisons
  and a clearly labeled regional soil pool without inferring colonization.

Outputs go to `data/derived/analysis_inputs/`. Small audit CSVs are suitable for
version control; large `.rds` and `.csv.gz` matrices are reproducible build
products and are ignored by Git.

## Primary FT-ICR-MS table

The primary filter is fully declared and audited: measured mass 200–900,
`C13 <= 0`, detection at intensity greater than zero in at least two of the 60
site columns, and a formula assignment represented by `carbon_count > 0`.
No normalization or transformation occurs in data preparation.
