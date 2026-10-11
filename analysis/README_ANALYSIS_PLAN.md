# Analysis roadmap

[Six planned manuscript figures](MANUSCRIPT_FIGURES.md) lists current PDFs, scripts, and caption notes.

Updated October 10, 2026. The active manuscript compares soils, hyporheic
porewater, and sediments and integrates microbial and molecular composition
within sediment. Start with the [manuscript outline](SOIL_POREWATER_SEDIMENT_MANUSCRIPT.md),
[current findings](KEY_FINDINGS_2016.md), and [figure guide](../figures/README.md).

## 1. Completed preparation and source products

DADA2 inference, the screened ASV catalog, the master tree, and fixed 10K
rarefaction are complete. Routine manuscript work consumes these products;
rebuilding them requires a specific reason. Raw reads, environmental sources,
the FT-ICR workbook, and completed inference remain intact.

| Dependency | Documentation | Use |
| --- | --- | --- |
| Library identity and mapping | `build_sample_manifest.R`; [preparation guide](data_prep/README.md) | Preserve unresolved soil identity and explicit site translations |
| Completed DADA2 | [workflow](README_DADA2.md) | Counts, taxonomy, QC, and provenance |
| Original microbial preparation | `microbes/01_prepare_diversity.R` | Existing 100-library fixed 10K draw, not the reviewed cohort |
| Environmental and molecular inputs | [preparation guide](data_prep/README.md) | Audited joins, binary molecular inputs, and missingness |
| Original input loader | [loader guide](README_LOADING_DATA.md) | Read-only access to original preparation objects; no automatic manuscript exclusions |
| Master phylogeny and UniFrac | `microbes/03_build_phylogeny.R`, `05_unifrac.R` | Reuse the completed tree/distances; prune or subset for each analysis |

## 2. Current comparative workflow

Run these directly in Positron or source them from any directory inside the
project. Settings and seeds are defined in the scripts; no arguments are needed.

```r
source(here::here("analysis", "microbes", "18_prepare_manuscript_cohort.R"))
source(here::here("analysis", "microbes", "19_soil_porewater_sediment.R"))
```

Script 18 subsets the existing draw and records all 120 libraries in an
inclusion audit. It excludes the 15 exact suspect water libraries and omits
all planktonic communities from this manuscript. It retains 76 libraries:
15 soil, 34 sediment, and 27 porewater, with 18 sediment-porewater pairs.
The original draw has 20 sediment-porewater pairs; the older count of 26 is
not supported by that input. Suspected contamination or sample-identity
problems remain unresolved causes, not confirmed diagnoses.

Script 19 refits Hellinger PCA and Hill diversity, habitat distribution tests,
soil-inventory sharing and coverage sensitivity, same-site resemblance tests,
and porewater network/environment checks. Model-specific missingness never
imposes one complete-case filter on unrelated analyses. Soil communities form
a regional inventory, not paired aquatic observations. Sequencing plate is
technical metadata and is excluded from ecological predictors and conditioning.

Use [the cohort report](microbes/README_SOIL_POREWATER_SEDIMENT.md) for methods,
nulls, multiplicity families, captions, and results. Constructed inputs belong
in `data/derived/soil_porewater_sediment_2016/`, results in
`results/diversity_2016/soil_porewater_sediment/`, and figures in `figures/`.

## 3. Within-sediment integration

The sediment analyses retain all 34 microbial samples and 44 molecular
profiles when describing each dataset separately. Direct comparisons use the
identical 33-site overlap; site 43 remains in microbial-only work and in
eligible sediment-porewater comparisons. Missing conductivity or texture
restricts only the model requiring that measurement.

| Question | Scripts and documentation |
| --- | --- |
| Baseline composition and detection | `fticr/01_explore_sediment.R`, `02_river_composition.R`, `microbes/08_sediment_figures.R`; [molecular methods](fticr/README_RIVER_COMPOSITION.md) |
| Water optics and property axes | `fticr/03_dom_optical_integration.R`, `04_integrated_chemistry_pca.R`; [integration report](fticr/README_DOM_INTEGRATION.md) |
| Network dispersion | `fticr/05_network_dispersion.R`, `07_matched_sediment_dispersion.R`; [matched methods](fticr/README_MATCHED_SEDIMENT_DISPERSION.md) |
| Sediment conditions and paired adjustment | `fticr/15_sediment_characteristics.R`; [methods](fticr/README_SEDIMENT_CHARACTERISTICS.md) |
| Exploratory associations and selected groups | `fticr/08_bipartite_feasibility.R` through `11_stable_group_distribution.R`; [nulls and limits](fticr/README_BIPARTITE_NETWORKS.md) |
| Manuscript-ready review panels | `figures/01_sediment_manuscript.R`; [panel guide](figures/README_SEDIMENT_FIGURES.md) |

Molecular ecology uses incidence. Raw signal and detection count are QC
information; equal-feature thinning remains a sensitivity analysis. Whole-
composition concordance is weak, and the molecular-property association
weakens after sediment-condition adjustment. Formula signatures and selected
ASV associations do not establish metabolites, metabolic guilds, production,
or consumption. The core-2 branch hypothesis still requires audited geography.

## 4. Assembly analysis and deployment state

The four-habitat catchment iCAMP run completed October 9 as job 432103.
It includes suspect and planktonic libraries and is historical for this
manuscript. Its bins and null pool do not represent the reviewed 76-sample
cohort. See [completed results](microbes/README_ICAMP_RESULTS.md).

Sediment-only job 432228 was submitted with approval October 10 for all 34
sediment samples and 17,946 ASVs, using 1,000 randomizations. Its current
queue/run state has not been checked. Keep the deployed scripts and inputs
unchanged. The [deployment record](cluster/README_SEDIMENT_DEPLOYMENT.md)
contains its resources, log paths, and approved submission history; the
[rerun record](README_ICAMP_RERUN.md) documents the earlier failed and pilot jobs.

A reviewed three-habitat iCAMP analysis remains planned. It needs its own
community pool, binning, outputs, validation, and approved remote actions.
Do not infer its fractions by removing pairs from the old catchment results.
Different null pools answer different questions; their process percentages
are not interchangeable. Review binning and within-pool phylogenetic signal
before mechanistic interpretation.

## 5. Spatial follow-up and historical workflows

The existing soil routing backbone, frozen public terrain, coordinate audits,
and membership sensitivities remain available. See [spatial methods](microbes/README_SOIL_STREAM_LOCALIZATION.md).
Water and pooled-habitat localization results predate the exclusions and need
cohort-specific refitting before current manuscript use. Retain all 15 soils
in general sharing; use only resolved coordinates and stable memberships for
spatial claims. Cross-sectional similarity does not identify movement direction.

Original four-habitat PCA/diversity, water-community chemistry, microbial
water dispersion, and catchment iCAMP retain historical result tables.
Their obsolete scripts and figures were removed; see the
[removal manifest](history/removed_plankton_workflows_20261010.csv).
Script 15 now regenerates reviewed soil–hyporheic–sediment detection inventories. The [sediment-only outline](SEDIMENT_ONLY_MANUSCRIPT.md)
is an earlier alternative; its sediment results remain applicable.
Dated decisions are preserved in [project history](history/PROJECT_HISTORY.md).
Legacy mothur source files remain intact; the obsolete mixed-habitat R Markdown
analysis is recorded in the removal manifest.

## 6. Reproducibility and housekeeping

Use `here::here()` paths, explicit settings/seeds, inspectable inclusion audits,
and model-specific sample counts. Keep tables in `results/`, constructed
inputs in `data/derived/`, and figures in `figures/`. Preserve raw sources and
completed inference, and fail on identifier/count mismatches. Do not edit a
script while it is running.

The [input guide](../data/README.md), [result guide](../results/README.md),
and [figure guide](../figures/README.md) distinguish current, supporting,
and historical products. Main figure outputs go to `figures/manuscript/`; supporting panels, alternatives,
and diagnostics go to `figures/supplement/`. Commits, pushes, and each
remote action require their own authorization; tidy-up does not include them.
