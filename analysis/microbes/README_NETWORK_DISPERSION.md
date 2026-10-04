# Aquatic microbial beta diversity across the 2016 stream network

Run `source(here::here("analysis", "microbes", "09_network_dispersion.R"))`
in R or Positron. The script uses the existing fixed 10K ASV table and sample
metadata. It does not rerun DADA2, rarefaction, or phylogeny construction.

## Question and design

The question is whether among-site microbial compositional heterogeneity
declines from headwater tributaries to the mainstem, as in the molecular
equal-feature dispersion figure. Headwaters are first/second order,
intermediate sites are third/fourth order, and mainstem sites are fifth order.
These are the same definitions used in `analysis/fticr/05_network_dispersion.R`.
Order describes network branching, not simply distance from the basin outlet.
All retained aquatic samples are analyzed, separately by habitat; regional soils
are excluded because they do not represent paired aquatic sites.

| Habitat | Headwater samples (segments) | Intermediate samples (segments) | Mainstem samples (segments) |
| --- | ---: | ---: | ---: |
| Sediment | 19 (13) | 10 (9) | 5 (5) |
| Hyporheic | 15 (12) | 11 (10) | 5 (5) |
| Planktonic | 10 (7) | 7 (7) | 3 (3) |

Sample IDs remain unchanged, and every aquatic habitat/site combination is
unique. Coordinates and stream order are checked against the recorded segment
in `data/StreamCenterline.csv`. Basin-outlet stream distance is derived from
this centerline; the valley-local metadata distance is not used. The 34
sediment samples include site 43, which has no molecular profile. Molecular
dispersion uses 44 sites and has 33 microbial matches, so these full-site
microbial and molecular results are not a comparison on identical site sets.
For a direct figure comparison, use the separate
[matched 33-site analysis](../fticr/README_MATCHED_SEDIMENT_DISPERSION.md), which
refits both composition views and uses identical sites in every balanced draw.

## Distance measures and tests

Primary abundance composition uses Bray-Curtis dissimilarity on the fixed
10,000-read communities. Equal totals make this equivalent to Bray-Curtis on
relative abundances. Binary Jaccard is the incidence sensitivity, using the
same table. It gives equal weight to each detected ASV and remains sensitive
to observed richness differences. Microbial read rarefaction and molecular
equal-feature thinning standardize different quantities and do not make their
dispersion scales directly comparable. See the
[vegan distance documentation](https://vegandevs.github.io/vegan/reference/vegdist.html).

For each habitat/metric, calculate distance to the group spatial median with
`vegan::betadisper(type = "median", bias.adjust = TRUE)`. The adjustment
multiplies distances by sqrt(n/(n-1)), addressing downward bias with unequal
group sizes. Spatial medians follow the package's default, avoiding its
documented anti-conservative centroid permutation procedure. This use of
dispersion as beta diversity follows
[Anderson et al. (2006)](https://doi.org/10.1111/j.1461-0248.2006.00926.x) and the
[vegan documentation](https://vegandevs.github.io/vegan/reference/betadisper.html).
The hypothesis concerns heterogeneity around group centers, not a shift in mean
composition. Box plots show site distributions, not uncertainty in group means.

Estimate the Lingoes correction separately on each full habitat distance
matrix and retain this geometry in resampled subsets. The recorded corrections
are zero or numerical roundoff (maximum distance change below 1.2e-15), so
they do not materially alter these distance scales.

Each three-group PERMDISP uses 9,999 residual permutations. Each test has its
own seed recorded in `test_design.csv`. Apply BH correction across the three
primary Bray-Curtis tests, and separately across the three Jaccard sensitivity
tests. Secondary pairwise tests use vegan's pairwise permutation output; apply
BH separately across all nine contrasts for each metric. These are omnibus
group comparisons, not formal tests of a monotonic stream-order trend.

Pairwise raw dissimilarity and geographic separation are exported as
descriptive diagnostics. Shared samples make pairs dependent; no pairwise row
is treated as an independent replicate in a test. Unrestricted residual
permutations are exploratory because exchangeability is not established in
this spatially connected network. Site groups differ in tributary identity,
sampling extent, and catchment conditions; stream order does not isolate a
causal connectivity effect.

## Segment-balanced sensitivity

Draw 1,000 balanced site sets per habitat. Each set selects distinct mapped
headwater and intermediate segments, then one sample per selected segment,
while retaining every mainstem sample. There are five sites per group for
sediment and hyporheic communities and three per group for planktonic
communities. All groups participate in every draw. Per-draw seeds and exact
sample selections are retained. Refit spatial medians in each subset using the
full habitat's distance correction and retain failures or warnings as flagged
rows. Report dispersion ratios and the fraction of successful draws with a
larger value in the first named group. These distributions are sensitivity
summaries, not confidence intervals, permutation P values, or new independent
observations. Balancing does not remove geographic confounding, and planktonic
mainstem dispersion remains estimated from only three sites.

## Initial results verified on 2026-10-04

| Habitat | Bray-Curtis headwater mean | Intermediate mean | Mainstem mean | Headwater/mainstem ratio | Three-group P | BH-adjusted P |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Sediment | 0.4480 | 0.4368 | 0.4832 | 0.927 | 0.5007 | 0.7511 |
| Hyporheic | 0.6634 | 0.6714 | 0.6537 | 1.015 | 0.4226 | 0.7511 |
| Planktonic | 0.6470 | 0.6043 | 0.5790 | 1.118 | 0.8221 | 0.8221 |

Planktonic mean dispersion has the hypothesized ordering, but the three-group
test does not support a difference. Sediment mainstem dispersion is larger
than headwater dispersion, and hyporheic group means differ little. Binary
Jaccard yields the same qualitative distinctions: headwater/mainstem ratios
are 0.964, 1.003, and 1.067, respectively; all adjusted omnibus P values are
0.8609. Thus, these comparisons do not establish a general downstream decline
in microbial beta diversity. They also do not demonstrate equivalent
heterogeneity, given the small mainstem samples and broad planktonic variation.

All nine Bray-Curtis and nine Jaccard secondary contrasts are nonsignificant
after their respective corrections (minimum adjusted P = 0.7345 and 0.8952).
Segment-balanced summaries retain habitat differences:

| Habitat | Median balanced Bray-Curtis headwater/mainstem ratio | Draws with headwater > mainstem | Median balanced Jaccard ratio | Draws with headwater > mainstem |
| --- | ---: | ---: | ---: | ---: |
| Sediment | 0.918 | 2.9% | 0.963 | 5.6% |
| Hyporheic | 1.017 | 83.5% | 1.005 | 66.2% |
| Planktonic | 1.110 | 83.4% | 1.066 | 82.5% |

Sediment therefore retains the opposite direction from the molecular
equal-feature result in most balanced draws. Hyporheic differences are small
despite the frequent positive ratios. Planktonic Bray-Curtis ratios span
0.714-1.206 across draws, showing sensitivity to which tributaries are sampled.
These draw fractions do not establish statistical significance. Mean
headwater/mainstem geographic separations are 5.54/3.21 km for sediment,
5.14/3.54 km for hyporheic samples, and 6.72/5.11 km for planktonic samples;
the sampled spatial extent remains unequal after matching group sizes.

The final analysis completed without analytical warnings; all 6,000 balanced
fits succeeded without diagnostic flags. All draws use distinct segments
within each group. All three analyzed input checksums remained unchanged.
Every PDF was rendered using Poppler through `pdftools` and inspected visually
for labels, sample counts, clipping, and point/box alignment. No remote action,
commit, or push was performed, and prior local changes were preserved.

## Figures and reproducibility outputs

All figures are vector PDFs in `figures/`. Each individual panel uses the same
group colors and box-and-point layout as `2016_fticr_dispersion_equal_features.pdf`.
The x-axis reports group sample counts. The habitat is identified on the y-axis;
no overall title, subtitle, or explanatory footnote is placed on the figure.

- `2016_sediment_microbes_dispersion_bray.pdf`
- `2016_hyporheic_microbes_dispersion_bray.pdf`
- `2016_planktonic_microbes_dispersion_bray.pdf`
- `2016_sediment_microbes_dispersion_jaccard.pdf`
- `2016_hyporheic_microbes_dispersion_jaccard.pdf`
- `2016_planktonic_microbes_dispersion_jaccard.pdf`
- `2016_aquatic_microbes_dispersion_bray.pdf` combines the three primary panels.

Draft caption for the primary figure: Among-site microbial compositional
dispersion in first/second-order headwaters, third/fourth-order intermediate
sites, and fifth-order mainstem sites. Points represent bias-adjusted distances
to habitat-specific group spatial medians in Bray-Curtis space, calculated from
the fixed 10,000-read ASV table. Boxes describe site distributions, with all
sites plotted. Habitats are fitted independently and have separate y-axis
scales. The Jaccard panels use ASV incidence with otherwise identical methods.

Tables are in `results/diversity_2016/network_dispersion/`, including primary
tests, secondary contrasts, group sizes, sample/network audits, pairwise
diagnostics, distance corrections, balanced draws and sample selections,
input checksums, settings, and session information. No sequencing plate is
included as an ecological predictor or conditioning covariate.

## Deferred beta-NTI analysis

Beta-NTI has not been calculated. The master tree named in the project
continuity notes is absent in this checkout; no remote check or transfer has
been performed. Unlike dispersion, beta-NTI is a signed pairwise standardized
departure of beta mean nearest-taxon distance from a phylogenetic null
expectation. See
[Stegen et al. (2013)](https://doi.org/10.1038/ismej.2013.93).
It must not be treated as an ordinary Euclidean distance matrix in PERMDISP.

Before calculation, verify the existing tree's tips and branch lengths, match
the community ASVs, evaluate relevant phylogenetic signal, specify the null
pool and abundance weighting, and benchmark a representative local run.
Within-habitat pools and a catchment pool answer different questions; use the
same declared pool across network groups within a comparison, rather than
changing the pool for headwaters and mainstem. Null runs must retain degenerate
null standard deviations and other failures as diagnostic rows. Future group
comparisons can summarize within-group signed beta-NTI and threshold fractions
with site-level permutations, retaining shared-sample dependence. The existing
iCAMP bMPD/Confidence settings do not supply whole-community beta-NTI. This
extension remains prospective and does not alter the current conclusions.
