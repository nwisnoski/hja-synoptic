# Analysis sequence

This roadmap separates reproducible data construction from ecological analysis.
Scripts should consume derived tables; they should never edit the master CSV,
the FT-ICR-MS workbook, FASTQ files, or DADA2 results in place.

## 0. Build ASVs

Run the [DADA2 workflow](README_DADA2.md). Its completion marker is
`results/dada2_2016/session_info.txt`. Until that file exists, downstream
microbial scripts must not consume the output directory.

## 1. Construct analysis inputs

Run [the data-preparation workflow](data_prep/README.md):

```r
system2(
  file.path(R.home("bin"), "Rscript"),
  "analysis/data_prep/run_data_preparation.R"
)
```

The first three steps can run while DADA2 is computing. The same command later
adds the DADA2-dependent 44-site sediment multiblock and source-pool tables.

Load and validate all available objects with:

```r
source("analysis/00_load_analysis_data.R")
hja <- load_hja_analysis_data()
```

See [the loader documentation](README_LOADING_DATA.md) for object names and
the explicit in-memory analysis views.

## 2. Microbial habitat and spatial patterns

Planned script group: `analysis/microbes/`.

1. Report sample counts and site overlap before hypothesis tests.
2. Summarize sequencing depth, ASV richness, and composition by planktonic,
   hyporheic, sediment, and regional soil habitat.
3. Test habitat differences with methods that respect unequal sample sizes.
4. Within aquatic habitats, relate composition to stream order, drainage area,
   valley geomorphology, and spatial coordinates.
5. Treat repeated site codes as blocks where habitats are compared at a site;
   do not imply a complete three-habitat design.

The terrestrial soils are a regional comparison set, not paired observations
from the aquatic sites.

## 3. Sediment FT-ICR-MS patterns

Planned script group: `analysis/fticr/`.

1. Describe retained-feature counts and total signal, with laboratory standards
   reported separately.
2. Examine molecular composition using both primary intensities and
   presence/absence as sensitivity views.
3. Relate molecular patterns to landscape position, hydrology, sediment
   organic content and enzyme activities, surface/hyporheic DOM optics, and
   nutrient chemistry.
4. Use the environmental missingness audit to define a complete primary
   predictor set; reserve discharge and sediment texture for smaller sensitivity
   subsets.

The feature filter is fixed in `data_prep/config.R` and is not re-created inside
analysis scripts.

## 4. Paired sediment integration

Planned script group: `analysis/integration/`. The primary sample universe is
the 44 sites with sediment ASVs, an FT-ICR-MS profile, and master-table rows.

1. Compare whole-community microbial and molecular dissimilarity patterns.
2. Ask whether shared environmental gradients explain both data blocks.
3. Relate predeclared taxonomic or putative functional groups to molecular
   classes and measured process proxies such as NAG, LAP, GLU, AP, organic
   content, nutrients, and optical DOM indices.
4. Label taxon–molecule associations as hypotheses or signatures, not direct
   evidence of microbial function or metabolite production.
5. Keep high-dimensional pairwise association searches explicitly exploratory
   and use multiple-testing control plus stability/sensitivity checks.

## 5. Potential source context

Same-site source comparisons are available only where the corresponding
libraries exist: 15 paired sediment–planktonic sites, 26
sediment–hyporheic sites, and 10 sites with both water sources. Regional soil
samples can describe a terrestrial source pool but cannot support same-site
source attribution. These analyses can compare ASV sharing or compositional
similarity; they should not infer colonization direction from this cross-section
alone.

## Reproducibility rule

Every analysis script should:

- read only declared inputs from `data/derived/analysis_inputs/` or completed
  DADA2 outputs;
- record all exclusions, transformations, and model formulas in code;
- write tables and figures to a dedicated generated-results directory;
- save a small run summary and `sessionInfo()`;
- fail on duplicate or misaligned sample identifiers.
