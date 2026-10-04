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

Implemented script group: [`analysis/microbes/`](microbes/README.md).

1. Report sample counts and site overlap before hypothesis tests.
2. Summarize sequencing depth, ASV richness, and composition by planktonic,
   hyporheic, sediment, and regional soil habitat.
3. Test habitat differences with methods that respect unequal sample sizes.
4. Within aquatic habitats, relate composition to stream order, drainage area,
   valley geomorphology, and spatial coordinates.
5. Treat repeated site codes as blocks where habitats are compared at a site;
   do not imply a complete three-habitat design.

Run the complete baseline with:

```sh
Rscript analysis/microbes/01_prepare_diversity.R
Rscript analysis/microbes/02_basic_diversity.R
Rscript analysis/microbes/03_build_phylogeny.R
Rscript analysis/microbes/04_network_environment.R
```

These scripts assume the repository root is the working folder and take no
command-line arguments. The first writes inspectable flat tables to
`data/derived/microbial_diversity_2016/`; the second reads those tables and
writes the diversity results and figures.

The network/environment script first makes descriptive figures within habitats;
it does not launch a grid of models or permutation tests. Formal tests should be
chosen only after inspecting those patterns. Catchment-gradient plots use
drainage area and distance to the outlet; formal stream-network autocorrelation
will require a connected reach topology or pairwise along-network distances
that are not currently in the derived data.

The terrestrial soils are a regional comparison set, not paired observations
from the aquatic sites.

All microbial analyses retain assigned Bacteria and Archaea and remove
chloroplast and mitochondrial sequences before depth filtering, rarefaction,
diversity estimation, or phylogeny construction.

The phylogeny is constructed once from the full screened ASV catalog, not from
a rarefied table. Each phylogenetic analysis then prunes that master tree to
the ASVs present in its selected community table.

## 3. Sediment FT-ICR-MS patterns

Initial descriptive workflow implemented in [`analysis/fticr/`](fticr/README.md).
`01_explore_sediment.R` uses presence/absence for molecular analyses and saves
Jaccard PCoA, van Krevelen, and feature-count figures. Raw intensities are used
only for detection and QC diagnostics. `02_river_composition.R` now explores
chemical signatures, drainage associations, replacement/nestedness-resultant
components, and repeated equal-feature subsampling. See the
[river-composition report](fticr/README_RIVER_COMPOSITION.md) for results and
published PNNL/GitHub precedents. Sediment-property models remain planned.

`03_dom_optical_integration.R` now compares EEM peak summaries, fluorescence
index, and SUVA254 across water compartments and against molecular composition.
`04_integrated_chemistry_pca.R` supplies FT-ICR property and joint chemistry
axes with explicit block weighting. See [DOM integration](fticr/README_DOM_INTEGRATION.md).

`05_network_dispersion.R` now uses all 44 sites in three-group dispersion and
richness-constrained Raup–Crick comparisons (headwater, intermediate, mainstem).
The equal-feature sensitivity shows an ordered decline in dispersion, but
complete-profile, chemical-property, replacement, and RC omnibus tests are
not significant. Balanced draws retain six sites per group on distinct mapped
segments. See [network dispersion](fticr/README_NETWORK_DISPERSION.md), including
RC saturation and geographic-extent limitations.

`06_aquatic_microbe_chemistry.R` adds habitat-matched aquatic comparisons:
20 planktonic communities/surface optics and 31 hyporheic communities/hyporheic
optics. Contextual sediment FT-ICR matches are 19/30, using fixed property-PC
projections from the earlier reference. None of four partial RDA models survives
correction. See [aquatic microbial chemistry](fticr/README_AQUATIC_MICROBE_CHEMISTRY.md).

1. Describe retained-feature counts and total signal, with laboratory standards
   reported separately.
2. Examine molecular composition using presence/absence only, following the
   user's decision on 2026-10-03. Retain intensities for detection/QC only.
3. Relate molecular patterns to landscape position, hydrology, sediment
   organic content and enzyme activities, surface/hyporheic DOM optics, and
   nutrient chemistry.
4. Use the environmental missingness audit to define a complete primary
   predictor set; reserve discharge and sediment texture for smaller sensitivity
   subsets.

The feature filter is fixed in `data_prep/config.R` and is not re-created inside
analysis scripts.

## 4. Paired sediment integration

Planned script group: `analysis/integration/`. The initial biogeochemistry
universe is 44 sites with sediment libraries, molecular profiles, and
master-table rows. The verified fixed-10K microbial intersection is 33 sites:
site `43` has a retained microbial sample but no molecular profile. Paired
analyses must use this intersection rather than assume all 34 retained sediment
samples have molecular data.

Initial paired dissimilarity, two-axis Procrustes, and one chemical-signature
model are implemented in `fticr/02_river_composition.R`. They provide no robust
overall molecular–microbial concordance under the examined comparisons; specific
taxon/formula associations and variation partitioning remain planned.

The subsequent chemistry-axis partial RDA finds a modest FT-ICR property
association after conditioning on drainage and detection count;
adding water optics does not strengthen it. Broad environmental variation
partitioning and taxon/formula associations remain prospective.

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

- read only declared inputs using paths relative to the repository root;
- record all exclusions, transformations, and model formulas in code;
- write tables and figures to a dedicated generated-results directory;
- write inspectable flat tables when an intermediate data product is needed;
- fail when sample identifiers or count tables are misaligned.
