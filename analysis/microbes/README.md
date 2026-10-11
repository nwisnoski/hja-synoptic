# 2016 microbial analysis

The current manuscript compares **15 soils, 27 hyporheic porewater, and 34
sediment communities**, retaining the existing 10,000-read draw. Start with
[the six planned figures](../MANUSCRIPT_FIGURES.md) and
[reviewed cohort methods](README_SOIL_POREWATER_SEDIMENT.md).

## Current local workflow

```r
source(here::here("analysis", "microbes", "18_prepare_manuscript_cohort.R"))
source(here::here("analysis", "microbes", "19_soil_porewater_sediment.R"))
source(here::here("analysis", "microbes", "15_first_detection_flowpaths.R"))
source(here::here("analysis", "figures", "02_current_manuscript.R"))
```

Script 18 records exclusions for all 120 original libraries and constructs
`data/derived/soil_porewater_sediment_2016/`. It excludes 15 exact suspect water
libraries, then omits all remaining planktonic libraries from manuscript scope.
It does not rarefy again or diagnose the cause of the anomalies.

Script 19 refits habitat diversity, Hellinger PCA, soil-inventory overlap,
coverage sensitivities, same-site comparisons at 18 sediment–porewater pairs,
and porewater environmental checks. Model-specific missingness remains explicit.
Script 15 generates the [reviewed first-detection inventory](README_FIRST_DETECTION.md)
with sediment and hyporheic detections as parallel pools. Figure script 02
combines Hill diversity panels and refits weighted UniFrac PCoA on the reviewed
76-library submatrix of saved distances. No phylogeny rebuild is required.

## Retained dependencies and sediment analyses

- `01_prepare_diversity.R`: original screened catalog, metadata, and fixed
  100-library rarefaction; preserve the completed draw for provenance.
- `03_build_phylogeny.R`: completed master tree; do not routinely rebuild.
- `05_unifrac.R`: shared original distance-matrix construction and source PCoA
  tables. It no longer creates four-habitat figures. Use saved distances.
- `06_icamp_sediment.R`: deployed sediment-only inference. Job 432228 was
  submitted October 10; its current state has not been checked. Preserve its
  code and inputs. See [deployment](../cluster/README_SEDIMENT_DEPLOYMENT.md).
- `08_sediment_figures.R`: 34-sample sediment Hellinger PCA and q = 1 panels.
- `10_soil_drainage_routing.py`: frozen-terrain routing and identity/QC backbone;
  see [spatial provenance](README_SOIL_STREAM_LOCALIZATION.md).

[The analysis roadmap](../README_ANALYSIS_PLAN.md) links the FT-ICR, sediment
condition, matched dispersion, and exploratory association workflows.

## Removed historical workflows

At Nathan's request, plankton-inclusive analyses and unreviewed water figures
were removed from the working tree on October 10. The
[removal manifest](../history/removed_plankton_workflows_20261010.csv) records
13 scripts and 47 PDFs, their checksums, and recovery commit `908d201`.
Retrieve a historical file with `git show 908d201:path/to/file` when needed.
The historical reports and compact result tables remain for provenance.
Raw sources, completed DADA2, the master tree, and inference outputs remain
intact. The completed four-habitat catchment iCAMP includes excluded libraries
and is historical; a reviewed three-habitat run requires a new null pool.
