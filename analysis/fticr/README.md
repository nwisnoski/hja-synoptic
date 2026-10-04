# Sediment FT-ICR-MS exploration

## Run

If the ignored molecular matrices are missing, rebuild them from the immutable
workbook using the established filter:

```sh
Rscript analysis/data_prep/03_prepare_fticr.R
```

Run the analysis scripts directly in R or Positron, without arguments:

```r
source(here::here("analysis", "fticr", "01_explore_sediment.R"))
source(here::here("analysis", "fticr", "02_river_composition.R"))
source(here::here("analysis", "fticr", "03_dom_optical_integration.R"))
source(here::here("analysis", "fticr", "04_integrated_chemistry_pca.R"))
source(here::here("analysis", "fticr", "05_network_dispersion.R"))
source(here::here("analysis", "fticr", "06_aquatic_microbe_chemistry.R"))
source(here::here("analysis", "microbes", "08_sediment_figures.R"))
```

All use `here()` paths and work from directories within the repository. They
preserve existing baseline figures and do not rerun DADA2, rarefaction, or
phylogeny. Molecular analyses use **presence/absence only**, as specified on
2026-10-03. Raw intensities are retained only for detection and QC summaries.

## Samples and features

The fixed filter in `analysis/data_prep/config.R` retains measured mass 200-900,
C13 indicator <= 0, detection at positive intensity in at least two of 60 site
profiles, and carbon count > 0. It retains 4,760 peak features. Laboratory
standards are excluded from ecological analyses and summarized separately.
Their site-filtered QC summaries use the same site-defined feature catalog.

The initial analysis uses the existing 44-site sediment multiblock crosswalk.
Of the 4,760 primary features, 4,718 occur in this subset; zero columns remain
in the aligned matrix and do not affect Jaccard distances. Prevalence is
calculated across these 44 profiles without reapplying the two-site filter.
The remaining 16 molecular profiles remain in preparation outputs and the
membership audit but are outside this initial analysis universe.

The fixed 10K microbial table contains 34 sediment samples and 17,946 ASVs
detected in sediment. Only **33** sediment samples overlap the 44-site molecular
set: site `43` has a retained microbial sample but no FT-ICR-MS profile. Eleven
sites in the 44-site set have no sediment sample retained at 10K. The microbial
figures use all 34 sediment samples; paired integration must use the audited
33-site intersection and refit ordinations on that subset.

`results/fticr_2016/tables/site_overlap_audit.csv` records membership and analysis
roles for all molecular profiles plus the unmatched sediment site. Raw-source
checksums are verified against the preparation inventory before analysis.
No source worksheets or sample identifiers are changed.

## Methods and interpretation

- A peak is detected when its prepared intensity is > 0. The exported matrix
  `data/derived/fticr_2016/sediment_44_presence_absence.csv` contains site rows
  and binary peak columns. No new intensity or signal-to-noise cutoff is imposed.
- Molecular ordination uses binary Jaccard dissimilarity and PCoA. It requests
  a Lingoes correction and records both its constant and the uncorrected
  negative eigenvalue count. In this run, there are no appreciable negative
  eigenvalues and the constant is approximately 1.1e-16. Axis percentages
  refer to positive corrected eigenvalues.
- Van Krevelen coordinates are atomic H:C and O:C ratios from prepared formula
  metadata. Points represent retained peak features detected in the 44-site
  subset and color shows site prevalence. Shared ratios can overlap. Formula
  ratios and source class labels do not establish structure or metabolite identity.
- Site-level mean H:C and O:C ratios give equal weight to every detected feature.
  No relative-intensity normalization, intensity-weighted property summaries,
  or intensity-based molecular ordination is used.
- Microbial PCA uses Hellinger-transformed counts from the fixed 10K table,
  column centering, and no column scaling. Hill diversity at q = 0, 1, and 2
  is checked against the existing baseline. These microbial choices are
  separate from the molecular presence/absence decision.
- The initial drainage plots use logarithmic x axes and show observations.
  `02_river_composition.R` adds exploratory chemical gradients, a descriptive
  H:C fit, molecular distance sensitivities, and paired microbial comparisons.
  See [the river-composition analysis](README_RIVER_COMPOSITION.md) for methods,
  results, published PNNL workflows, inspected GitHub scripts, and limitations.

Presence/absence avoids treating uncalibrated relative signal as molecular
abundance. Ionization efficiency and matrix effects can alter signal, as
demonstrated by [Lechtenfeld et al. (2024)](https://pubs.acs.org/doi/10.1021/acs.est.3c07219).
Detection itself remains affected by instrument response and sample preparation;
an absent peak means nondetection in this measurement rather than established
absence from the sample.

Filtering, incidence summaries, and formula-ratio visualization follow the
exploratory framework in [Bramer et al. (2020)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007654).
Microbial Hellinger ordination follows [Legendre and Gallagher (2001)](https://doi.org/10.1007/s004420100716).
These scripts use vegan and ggplot2 rather than ftmsRanalysis.

## Figures and draft captions

All figures are in `figures/`. Four sparse panels use vector PDF; the dense
van Krevelen panel uses a white-background, 600-dpi PNG.

| File | Draft caption |
| --- | --- |
| `2016_fticr_sediment_jaccard_pcoa.pdf` | PCoA of binary Jaccard dissimilarity in detected primary molecular features at 44 sediment sites. Color indicates log10 drainage area in hectares. The first two axes explain 27.3% and 10.5% of positive corrected eigenvalues. |
| `2016_fticr_sediment_richness_drainage.pdf` | Detected primary molecular feature count against drainage area at 44 sediment sites. Counts describe the filtered profiles and remain sensitive to differences in detection. |
| `2016_fticr_sediment_van_krevelen.png` | Van Krevelen distribution of 4,718 primary peak features detected in the 44-site sediment subset. Color indicates the proportion of sites detecting each feature. |
| `2016_sediment_microbes_hellinger_pca.pdf` | Hellinger PCA of sediment microbial composition for 34 samples at the fixed 10,000-read depth. Color indicates log10 drainage area in hectares. PC1 and PC2 explain 16.9% and 13.6% of variance, respectively. |
| `2016_sediment_microbes_q1_drainage.pdf` | Sediment microbial Hill diversity at q = 1 against drainage area for 34 samples at the fixed 10,000-read depth. Each point represents one sample. |

## Diagnostics and next analysis

Detected molecular feature counts range from 391 to 2,638 per site and correlate
strongly with total primary signal (site-level Spearman rho = 0.969). This is a
QC diagnostic, not an intensity-based ecological analysis. Detection differences
need investigation before interpreting feature counts as ecological richness.

River-gradient, source-class, DOM optical-index, and 33-site paired Procrustes
comparisons are now implemented in `02_river_composition.R`. PCoA1 strongly tracks
detection count (rho = −0.984). The small drainage association is sensitive to
richness treatment; the examined chemical signatures do not show a convincing
monotonic downstream shift. Overall molecular–microbial concordance is weak.
The new report supplies the seven additional figure captions and result tables.

`03_dom_optical_integration.R` adds 59-site surface and 56-site hyporheic water
optics, with 44/42 molecular matches. `04_integrated_chemistry_pca.R` fits a
42-site chemistry PCA and tests its first two axes against the 33 matched
microbial communities. See [DOM and chemistry integration](README_DOM_INTEGRATION.md)
for the measurement definitions, block weighting, additional figures, and results.
The FT-ICR property PCs have a modest adjusted microbial association (4.04%,
P = 0.0103; FDR-adjusted P = 0.0309); joint chemistry PCs do not strengthen it.
These property summaries differ from the earlier overall formula-distance tests.

`05_network_dispersion.R` uses all 44 sites in three-group tests: 24 first/second
order headwaters, 14 third/fourth order intermediate sites, and six fifth-order
mainstem sites. It tests dispersion and occurrence-weighted, richness-constrained
Raup–Crick, with 1,000 balanced draws of six distinct segments per group.
The equal-feature sensitivity shows decreasing dispersion across those groups
(omnibus adjusted P = 0.003); complete profiles, replacement, chemical properties,
and RC do not show significant omnibus differences. RC is strongly saturated.
See [network dispersion](README_NETWORK_DISPERSION.md) for all comparisons,
spatial limitations, provenance, and six figure captions.

`06_aquatic_microbe_chemistry.R` matches all 20 planktonic communities to
surface-water optics and all 31 hyporheic communities to hyporheic-water optics.
It also compares 19/30 communities with contextual sediment FT-ICR signatures,
extending the fixed property-PC projection to all 60 profiles. None of the four
partial RDA models survives BH correction; matched water optical associations
remain unconfirmed after conditioning on drainage. See
[aquatic microbial chemistry](README_AQUATIC_MICROBE_CHEMISTRY.md) for sample
exclusions, the water/sediment distinction, sensitivity checks, and three figures.

Next, select a small sediment-property or enzyme comparison using the missingness
audit. Partial RDA now tests chemistry after
conditioning on drainage and detection count; broader environmental
variation partitioning remains planned.
Retain the complete binary profiles as primary data; the 1,000 equal-feature
draws are a detection-count sensitivity, not standardized analytical sampling.

Sequencing plate is excluded from ecological predictors and conditioning
covariates and retained only as technical QC metadata.

## Verification on 2026-10-03

All eight scripts completed locally (FT-ICR preparation, six FT-ICR analyses,
and sediment microbial figures). All four inventoried source checksums matched,
feature identifiers matched the preparation filter, primary intensities were
finite and nonnegative, and every profile had positive total signal. Microbial
Hill summaries agreed with the baseline to tolerance 1e-8. The PDFs were
rendered with Poppler through R's pdftools and visually inspected; the PNG was
inspected separately. Outputs include inspectable tables and session records.
