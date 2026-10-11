# Reviewed soil-porewater-sediment analysis

October 10, 2026. This is the current comparative manuscript workflow. The
four-habitat analyses elsewhere in this repository remain historical outputs
and include libraries excluded here. Use the
[current manuscript outline](../SOIL_POREWATER_SEDIMENT_MANUSCRIPT.md) for framing.

## Run locally

Both scripts run directly in Positron from any directory within the repository.
They reuse the completed fixed 10K draw and saved UniFrac distances. They do
not rerun DADA2, rarefaction, phylogeny, or iCAMP.

```r
source(here::here("analysis", "microbes", "18_prepare_manuscript_cohort.R"))
source(here::here("analysis", "microbes", "19_soil_porewater_sediment.R"))
```

Constructed inputs are in `data/derived/soil_porewater_sediment_2016/`.
Results are in `results/diversity_2016/soil_porewater_sediment/`, and all 12
vector figures have prefix `figures/2016_soil_porewater_sediment_`.

Start figure review with the [sampling map](../../figures/2016_soil_porewater_sediment_sampling_map.pdf),
[habitat context](../../figures/2016_soil_porewater_sediment_habitat_context.pdf),
[soil sharing](../../figures/2016_soil_porewater_sediment_soil_sharing.pdf), and
[same-site comparison](../../figures/2016_soil_porewater_sediment_same_site_comparison.pdf).

The baseline `microbial_diversity_2016` directory remains the immutable source
of the existing rarefaction and deployed sediment-iCAMP inputs. Its 100-sample
table is not the reviewed manuscript cohort. Future analyses of the current
manuscript should consume the new cohort, not the unfiltered four-habitat table.

## Exact exclusions and cohort

Nathan authorized excluding the anomalous water block on October 10. Its cause
is unresolved: contamination and sample misidentification remain possibilities.
The exclusions are identity-based and checked against site, habitat, plate,
and lane-sample identifiers before any output is written. `S` labels below
are recorded sequencing identifiers, not verified physical extraction wells.

| Library | Site | Habitat | Plate 2 lane sample |
| --- | --- | --- | --- |
| hja2016_099 | 162 | Planktonic | S1 |
| hja2016_100 | 64 | Planktonic | S2 |
| hja2016_101 | 58 | Planktonic | S3 |
| hja2016_102 | 73 | Hyporheic | S4 |
| hja2016_103 | 260 | Planktonic | S5 |
| hja2016_105 | 178 | Planktonic | S6 |
| hja2016_106 | 170 | Planktonic | S7 |
| hja2016_107 | 57 | Hyporheic | S8 |
| hja2016_109 | 169 | Planktonic | S9 |
| hja2016_110 | 76 | Hyporheic | S10 |
| hja2016_111 | 160 | Planktonic | S11 |
| hja2016_114 | WS1-5 | Planktonic | S60 |
| hja2016_116 | 53 | Planktonic | S13 |
| hja2016_119 | 56 | Planktonic | S14 |
| hja2016_121 | 84 | Hyporheic | S15 |

Do not replace this table with a numeric ID range or exclusion of all Plate 2
libraries. Other Plate 2 samples remain eligible. The source conversation's
description of consecutive wells is not a verified plate map; the actual
recorded identifiers skip S12 and include S60. No ASVs are removed based on
the contaminant-associated genera mentioned in that conversation.

The preparation script writes a 120-library inclusion audit that separates
the 15 suspected-library exclusions, the established 10K threshold, and the
manuscript's exclusion of all planktonic communities. Nine unflagged planktonic
libraries are outside scope; they are not labeled contaminated. The new table
contains 76 libraries: 15 soils, 34 sediments, and 27 porewater samples. Its
33,974 ASVs are those detected in this cohort. Each library retains its exact
original 10,000 counts; only ASVs absent from every retained library are dropped.

The pairing audit independently joins sample identities by site. Twenty
sediment-porewater pairs occur in the original fixed 10K table; sites 57 and 84
lose their porewater library, leaving 18 pairs. Flagged porewater sites 73 and
76 do not have fixed-10K sediment counterparts. The previously reported 26
pairs are not supported by this fixed 10K input. Site 43 remains in the paired
comparison even though it lacks FT-ICR. Soil 04 lacks GPS and hja2016_200 has
unresolved site identity; both remain in the regional soil inventory, with
neither assigned a location or an aquatic-site pairing.

## Methods and null hypotheses

**Habitat context.** Hellinger PCA is refitted on all 76 samples. Hill q = 0,
1, and 2 use the fixed counts. Bray-Curtis PERMANOVA compares all three
habitats and each of three pairs, with 9,999 unrestricted permutations.
These four exploratory P values form one BH family. Each comparison also
has a spatial-median, bias-adjusted PERMDISP test; its four P values form a
separate BH family. Unrestricted exchanges do not control regional spatial
dependence or repeated aquatic sites. Strong dispersion differences prevent
interpreting PERMANOVA as pure centroid separation or a specific filtering
mechanism. A separate 36-library paired PERMANOVA conditions on the 18 sites
and exchanges habitat labels only within site; the site nuisance term has no
reported test. See the primary vegan documentation for
[PERMANOVA](https://vegandevs.github.io/vegan/reference/adonis.html),
[dispersion](https://vegandevs.github.io/vegan/reference/betadisper.html), and
[exchangeability](https://vegandevs.github.io/vegan/reference/permutations.html).

**Soil inventory.** An ASV is soil-detected if it occurs at least once in any
of the 15 fixed-depth soil libraries. For each sediment or porewater library,
report both the percentage of reads and percentage of detected ASVs belonging
to this inventory. Soil sharing is observed overlap, not source attribution.
One thousand independent, seeded draws of three distinct soil libraries hold
all receiving communities fixed and rebuild the inventory. Saved selections
make each draw reproducible. Quantiles describe variation among inventory
selections, not confidence intervals. This deliberately isolates soil-inventory
coverage; it differs from the historical workflow that also subsampled aquatic
habitat-stage cells. Do not substitute the old balanced-inventory percentages.

Drainage slopes use arithmetic means of sharing and mean log10 drainage within
each stream segment, separately by habitat. There are 27 sediment and 23
porewater segments. The two metrics in two habitats form a four-test BH family.
Lines are shown only for adjusted P <= 0.05. Segment aggregation avoids treating
multiple sites in one segment as independent replicates; residual spatial
dependence between segments remains possible. A separate paired soil-sharing
contrast exchanges the signs of the 18 within-site sediment-minus-porewater
read-percentage differences and uses a two-sided mean-difference statistic.

**Same-site resemblance.** Use exactly the same 18 paired sites for Bray-Curtis,
weighted UniFrac, unweighted UniFrac, and the percentage of porewater reads in
ASVs detected in sediment. UniFrac distances are subset from the completed
distance tables; neither tree nor distances are recalculated. For each metric,
the statistic is the mean of the 18 same-site comparisons. Under the null,
reassign intact porewater identities one-to-one among the sediment profiles.
This preserves each community's abundance profile and within-habitat covariance
while breaking site identity. Self matches can occur in randomized assignments.
The unrestricted null exchanges all 18 sites; the sensitivity exchanges sites
only within their exact stream order. Neither null fully controls geography
or shared segments.

For distances, the one-sided alternative is lower same-site distance; for read
sharing, it is greater same-site overlap. Monte Carlo P = (extreme + 1) /
(permutations + 1). Four metrics form one BH family for each null separately.
All null statistics, per-run seeds, and per-site observed comparisons are saved.
The plotted other-site value averages over the 17 nonmatching sediment sites
for each porewater sample. This differs slightly from the randomized baseline,
which permits self matches. Nonsignificance is not a test of equivalence.

**Porewater network and chemistry.** Refit Bray-Curtis network-group PERMANOVA
and median-based PERMDISP on the 27 porewater samples. Eight predeclared
continuous variables enter separate Spearman Mantel checks: drainage, porewater
DOC (the NPOC field), SUVA254, fluorescence index, dissolved oxygen, temperature,
TDN, and sediment hydraulic conductivity. Drainage, DOC, TDN, and conductivity
are log10-transformed; other variables retain their measured scale. Environmental
distance is absolute difference in each transformed variable. The positive-
association alternative asks whether greater environmental difference accompanies
greater community dissimilarity. Each model uses its own available sites and
records every inclusion/exclusion. Unrestricted and within-exact-order exchanges
form separate eight-test BH families. Following
[vegan's Mantel method](https://vegandevs.github.io/vegan/reference/mantel.html),
permutations exchange whole site labels, not individual pairwise distances.
These checks cannot exclude every form of environmental response.

## Completed results

| Comparison | Reviewed result |
| --- | --- |
| Three-habitat Bray-Curtis PERMANOVA | R² = 0.2018, P = 0.0001; dispersion P = 0.0001 |
| Paired sediment-porewater habitat contrast | Habitat R² = 0.1647 of total variation, conditional pseudo-F = 6.716, P = 0.0003 |
| Median richness, soil / sediment / porewater | 1,835 / 2,845.5 / 726 ASVs |
| Median soil-detected reads, sediment / porewater | 70.26% / 34.38% |
| Median soil-detected ASVs, sediment / porewater | 56.74% / 34.72% |
| Three-soil sensitivity, median habitat read percentages | 41.58% / 19.43%; draw quantiles 22.43-64.13% / 9.50-28.98% |
| Paired soil-sharing difference | Mean sediment excess 36.15 percentage points; two-sided P = 0.0001 |
| Soil-sharing drainage slopes | All four BH P >= 0.825 |
| Same-site Bray-Curtis | Mean 0.91898 vs randomized 0.91953; raw P = 0.4131 |
| Same-site weighted / unweighted UniFrac | Raw P = 0.8595 / 0.3664 |
| Sediment-detected porewater reads | Mean same-site 25.47%, offsite 24.40%; raw P = 0.1456 |
| Within-order pairing sensitivity | No metric BH-supported; minimum raw P = 0.1121 |
| Porewater network composition / dispersion | P = 0.2932 / 0.9266 |
| Porewater environmental associations | No BH-supported positive Mantel association; n = 16-27 |

The global habitat effect differs from the preliminary conversation's estimate;
the table above reports the reproducible R result. Habitat dispersion differences
are material. Same-site tests agree qualitatively with the preliminary check,
while the distinctions between means, medians, and randomized/offsite baselines
are now explicit. All four same-site metrics have BH P >= 0.5508 under the
unrestricted null and >= 0.3582 within order.

## Figure captions

All filenames below follow the prefix `2016_soil_porewater_sediment_` in the
top-level `figures/` directory. Figures are vector PDFs with white backgrounds.

| Suffix | Caption |
| --- | --- |
| `sampling_map.pdf` | Reviewed community sampling and sediment molecular coverage in the stream network. Separate panels show 13 mapped soils of 15 analyzed, 27 porewater communities, 34 sediment communities, and 44 molecular profiles. Black rings in the sediment panels mark the 33 microbial-molecular paired sites. Unmapped soils remain in nonspatial analyses. Coordinates are in UTM zone 10N. |
| `hellinger_pca.pdf` | PCA of Hellinger-transformed fixed-10K counts for all 76 reviewed communities. Color and shape identify habitat. Axes report explained total variance. |
| `diversity_q0.pdf`, `diversity_q1.pdf`, `diversity_q2.pdf` | Observed Hill diversity at 10,000 reads, separately for soil, sediment, and porewater. Points are individual libraries, with horizontal jitter only; boxes show medians and interquartile ranges. |
| `habitat_context.pdf` | Combined Hellinger PCA and q = 1 diversity panels; both use the reviewed three-habitat cohort. |
| `soil_read_percent.pdf`, `soil_asv_percent.pdf` | Soil-inventory overlap in each sediment and porewater community across drainage area. Read and detected-ASV fractions are shown separately. The inventory pools all 15 soils. No fitted slope survives the declared BH correction, so no regression line is drawn. |
| `soil_sharing.pdf` | The read-percentage and ASV-percentage soil-overlap panels combined. Absolute fractions depend on soil sampling coverage and do not estimate terrestrial-source percentages. |
| `same_site_distances.pdf` | Same-site dissimilarity versus mean dissimilarity to the 17 other paired sediment profiles for each of 18 porewater samples. Lines connect comparisons for the same porewater sample. Bray-Curtis and the two UniFrac metrics have separate vertical scales. Tests use whole-assignment randomization, not independent tests of the displayed pairs. |
| `same_site_read_sharing.pdf` | Percentage of each paired porewater library's reads in ASVs detected in its same-site sediment versus the mean across the other 17 paired sediment profiles. Lines retain porewater identity. Shared detections do not establish exchange direction. |
| `same_site_comparison.pdf` | The dissimilarity and read-sharing panels combined. These analyses find no supported preferential same-site resemblance under either declared null. |

## Provenance and remaining limitations

The preparation and analysis scripts check identifiers, counts, sample sizes,
pairing, and input hashes. All 12 PDFs were rendered with `pdftools` and visually
inspected. A separate numerical check compared the new cohort with the original
draw, sediment diversity with the existing summaries, pairing with source
metadata, and reported soil/same-site statistics with direct calculations.

The exclusions do not establish which handling step caused the anomaly. Resolve
that through plate maps and extraction records before describing confirmed
contamination. Low-biomass artifacts in retained porewater remain possible;
no ASV-level decontamination is claimed. Nathan's transient-porewater context
motivates interpretation, but is not inferred from these null comparisons.

All sediment microbial-molecular results remain available and use unchanged
samples. The completed four-habitat catchment iCAMP includes excluded libraries
and is historical for this manuscript. A reviewed three-habitat iCAMP run has
not been performed. Sediment-only job 432228 retains its deployed inputs and
scripts; its current state has not been checked. No remote commands, job
submissions, commits, or pushes accompanied this local analysis.
