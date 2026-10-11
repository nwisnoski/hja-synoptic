# Input guide

This repository covers the 2016 survey. Raw environmental CSVs, the FT-ICR-MS
workbook, sequencing/sample records, reads, and completed DADA2 products are
source records. Never edit them in place to implement sample exclusions.

| Input group | Use |
| --- | --- |
| `hja-synoptic_sequence-sample-list.csv` | Original library-to-sample identities; mapping exceptions are audited in code |
| `hja-env_data_clean.csv`, `hja-synoptic_env-data-soils.csv` | Aquatic and soil environmental sources |
| `hja-FTICRMS.xlsx` | Original molecular feature workbook; raw signal is detection/QC information |
| `StreamCenterline.csv` | Recorded stream centerlines, network order, and segment geography |
| `derived/analysis_inputs/` | Constructed environmental, FT-ICR, sediment, and source-pool inputs; [preparation](../analysis/data_prep/README.md) |
| `derived/microbial_diversity_2016/` | Original fixed 10K draw and full 120-library metadata; source for the reviewed cohort and deployed sediment iCAMP |
| `derived/soil_porewater_sediment_2016/` | Current 76-library manuscript cohort, exact exclusion ledger, and 18-site sediment-porewater pairing audit; [methods](../analysis/microbes/README_SOIL_POREWATER_SEDIMENT.md) |
| `derived/fticr_2016/` | Binary molecular ecological inputs and property summaries; [methods](../analysis/fticr/README.md) |
| `derived/sediment_characteristics_2016/` | Sediment-condition tables with model-specific availability; [methods](../analysis/fticr/README_SEDIMENT_CHARACTERISTICS.md) |
| `spatial/soil_routing/`, `derived/soil_drainage_2016/` | Frozen public terrain, routing products, and membership uncertainty; [spatial provenance](spatial/soil_routing/README.md) |
| Legacy mothur files | Historical OTU analysis; not inputs to the current ASV cohort |

The reviewed cohort removes the 15 exact suspect libraries and all planktonic
communities from manuscript analyses while preserving the retained libraries'
original rarefaction counts. It is distinct from the source 100-library table.
All 15 soils enter general sharing; only 13 have resolved mapped locations.
Missing profiles remain distinct from nondetection.

Constructed data belong here; analysis summaries belong in `results/` and
figures in `figures/`. All transformations and exclusions must remain in the
reproducible scripts. No source or constructed input files were moved or
changed during the October 10 repository tidy-up.
