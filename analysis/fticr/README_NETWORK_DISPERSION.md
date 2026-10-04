# Molecular heterogeneity across the stream network

Run `source(here::here("analysis", "fticr", "05_network_dispersion.R"))` in R or
Positron. This script uses existing preparation and river-composition outputs;
it does not change raw data or weight ecological comparisons by FT-ICR signal.

## Question and sampling design

The hypothesis is that molecular composition varies more among headwater sites
than among mainstem sites, even if mean composition changes little with drainage
area. Dispersion tests address variation around group centers; the earlier
drainage models addressed compositional location. These are different responses.

Headwaters are defined prospectively as first- and second-order streams: 24 sites
on 16 distinct mapped stream segments. Mainstem sites are the six fifth-order
sites, each on a separate segment, along the highest-order trunk visible in
`data/Segments686_streams.jpg`. The 14 third- and fourth-order sites on 13 segments
form the intermediate group in the **primary three-group tests using all 44
sites**, in every figure, and in the balanced-segment sensitivity. All three
pairwise group contrasts are secondary comparisons. Only three sites are first order, so
an isolated first-order contrast is not the primary analysis.
Branch order defines headwater status here: small tributaries can occur near
the basin outlet, so headwater does not mean greatest basin-outlet distance.

Coordinates and segment identifiers are matched to the nearest centerline vertex
within the recorded segment. All 44 sites match within 1.1e-8 m and every stream
order agrees. This is numerical agreement with the mapped coordinates, not
independent evidence of field positioning accuracy. The geometry documentation
(`data/_ReadMe_Geometry.docx`) identifies
centerline D as distance along the stream to the **basin** outlet. The previously
prepared `site_distance_to_outlet_m` originates from valley-segment distance and
must not be interpreted as basin-outlet distance. The audit preserves both fields.
The original centerline has one nonfinite outlet elevation; coordinates, stream
distance, drainage area, order, and segment remain finite. Elevation is not used.
For the map only, centerline vertices are thinned to approximately 10-m spacing.

## Distance views and inference

1. **Binary Jaccard:** all 4,718 features detected in the 44-site sediment subset.
2. **Jaccard replacement:** the Baselga component previously calculated and
   checked in `02_river_composition.R`. It separates replacement from the
   nestedness-resultant component but does not identify the cause of nestedness.
3. **Equal-feature Jaccard:** the mean distance from the previously completed
   1,000 uniform draws of 391 detected features per site. This is a sensitivity
   to detection-count differences, not laboratory-depth or coverage rarefaction.
4. **Formula properties:** Euclidean distance across all five standardized site
   summaries: mean detected H:C, mean detected O:C, N-bearing fraction,
   protein-like fraction, and tannin-like fraction. Each variable is standardized
   across the full 44-site set. No PCA axes are discarded. These partly correlated
   variables define the same FT-ICR property view used in the chemistry integration;
   they are not five independent chemical dimensions. Their labels describe
   formula-based signatures rather than confirmed carbon structures or sources.

For each view, use distance to the group spatial median, with small-sample bias
adjustment sqrt(n/(n-1)) because group sizes differ. PERMDISP uses 9,999 residual
permutations across all 44 sites with seed 20261008. Three-group tests use the
PERMDISP F statistic; the three secondary contrasts use vegan's pairwise
permutation procedure from that same fit. The spatial-median default avoids the
anti-conservative centroid permutation test documented by
[vegan](https://vegandevs.github.io/vegan/reference/betadisper.html).
The underlying method is
[Anderson et al. (2006)](https://doi.org/10.1111/j.1461-0248.2006.00926.x).
Lingoes constants are estimated from each full 44-site distance matrix and
retained when examining subsets. This puts subsets in a common geometry.
Correction size is recorded, particularly because replacement has substantial
non-Euclidean structure. Distances after that correction are not raw replacement.

## Richness-constrained Raup–Crick

The shared regional pool contains all 4,718 detected features from all 44 sites.
The `vegan` r1 null keeps each site's observed detection count and draws distinct
features with probabilities proportional to their frequencies across this pool.
Three preliminary simulations audit binary entries and exact site detection
counts. Intermediate sites contribute to both the common pool and all three-group
tests. Separate pools by group would change
the reference expectation and are not used.

There are 1,999 simulated 44-site communities, seed 20261007, generated in batches
of 20 on one process; the full null run took 215.6 s locally. A 99-simulation
benchmark took 10.9 s. Completed null matrices are reused only when the input
MD5 checksum, seed, simulation count, null method, and tie rule match. Setting
`reuse_matching_null <- FALSE` explicitly recomputes the simulations. For each
pair, the signed score is:

`RC = 2 * (fraction(null shared > observed shared) + 0.5 * fraction(ties)) - 1`

Thus +1 indicates fewer shared features than expected (greater dissimilarity),
−1 indicates more shared features (greater similarity), and zero indicates the
middle of the null expectation. Midpoint tie handling follows the
[PNNL-associated binary FT-ICR script](https://github.com/danczakre/Meta-Metabolome_Ecology/blob/master/FTICR_RCBC_PresAbs.R)
and the split-tie option of
[Chase et al. (2011)](https://doi.org/10.1890/ES10-00117.1).
[The vegan documentation](https://vegandevs.github.io/vegan/reference/raupcrick.html)
describes its occurrence-weighted null and alternative tie conventions. This
analysis obtains null shared counts with `oecosimu` and explicitly splits ties;
it does not use vegan's default unsplit ties. Binary Bray and Jaccard both vary
monotonically with shared-feature count conditional on the two fixed richness
values, so the shared-count comparison captures the same ordering of null pairs.

Within-group RC means are calculated separately for headwater, intermediate,
and mainstem sites. The omnibus statistic is the sum of squared deviations of
these three means from their unweighted mean. Shuffle **site labels** 9,999
times across all 44 sites, preserving group sizes 24/14/6, and compare this
statistic to its simulated reference. Secondary RC contrasts compare each
pair of group means against that same three-group label-permutation reference.
Report both absolute-contrast two-sided probabilities and one-sided probabilities
for the first group being more dissimilar. These RC contrasts shuffle labels
across all three groups rather than restricting each contrast to two groups.
Pairs are never treated as independent replicates.

Apply BH correction to the five primary omnibus tests (four dispersion views
and RC), then separately to the 15 secondary two-sided contrasts. Directional
RC probabilities are reported descriptively without a separate significance
claim. These views are correlated; they are prespecified comparisons rather
than independent confirmations. Do not apply Euclidean dispersion analysis
directly to the signed RC matrix.

Pairwise diagnostics retain null expected overlap, tie frequency, saturation,
and Monte Carlo standard error. With 1,999 simulations, signed scores have
Monte Carlo uncertainty up to approximately 0.022 near the middle of the null
distribution. Extreme ±1 scores mean no opposite outcome was observed in these
simulations, not an infinitely precise probability.

## Reach-density and sample-size sensitivity

Each of 1,000 draws selects six of the 16 headwater segments and then one site
within each selected segment. Likewise select six of the 13 intermediate
segments and one site per selected segment. Retain the same six distinct
mainstem segments. Every draw therefore contains 18 sites, six per group.
Per-run seeds are 20261009 + run number; the selected sites are exported.
Recalculate spatial-median dispersion with six sites per group and summarize
all three pairwise dispersion ratios for each distance view. Calculate within-group mean
RC differences using the original 44-site null reference. These draws address
unequal group size and repeated sampling within mapped segments. Their ranges
and fractions above zero are sensitivity summaries, not confidence intervals,
permutation P values, or new independent observations.

This cross-sectional comparison cannot separate network position from tributary
identity, geography, catchment conditions, or detection thresholds. Mapped
segments are not independent watersheds, and connected mainstem segments remain
dependent. Pairwise geographic distances are exported to characterize sampling
extent. Unrestricted site-label/residual permutations provide an exploratory
reference; they do not establish exchangeability across a spatial network.
Raup–Crick controls observed feature count under its stated null model; it does
not recover missed detections or establish an assembly mechanism.

## Results verified on 2026-10-03

| Distance view | Headwater mean | Intermediate mean | Mainstem mean | Three-group P | BH-adjusted P |
| --- | ---: | ---: | ---: | ---: | ---: |
| Complete binary Jaccard dispersion | 0.3723 | 0.3286 | 0.3625 | 0.3626 | 0.6043 |
| Corrected replacement dispersion | 0.7673 | 0.7719 | 0.7618 | 0.5803 | 0.7254 |
| Equal-feature mean Jaccard dispersion | 0.6405 | 0.6234 | 0.5803 | 0.0006 | 0.0030 |
| Five-property dispersion | 1.9863 | 1.9897 | 1.6320 | 0.7575 | 0.7575 |
| Within-group mean signed RC | −0.8762 | −1.0000 | −1.0000 | 0.3197 | 0.6043 |

Dispersion values are mean bias-adjusted distances to the spatial median;
RC values are within-group pairwise means. Their scales and tests differ.
The equal-feature sensitivity is the one view with convincing differences
under these exploratory permutation references. Its headwater/mainstem ratio
is 1.104; the secondary headwater–mainstem contrast has P = 0.0011 and
BH-adjusted P = 0.0165 across 15 contrasts. Intermediate dispersion falls
between the endpoints, but neither adjacent contrast survives that correction
(headwater–intermediate adjusted P = 0.1575;
intermediate–mainstem adjusted P = 0.1575).

Balanced six-site draws preserve the equal-feature pattern: headwater dispersion
exceeds mainstem dispersion in all 1,000 draws (median ratio 1.104), intermediate
exceeds mainstem in all draws (median ratio 1.074), and headwater exceeds
intermediate in 92.4% (median ratio 1.026). Chemical-property dispersion exceeds
mainstem in 71.9% of headwater draws and 72.8% of intermediate draws; this is not
strong evidence of an ordered chemical-property gradient. Complete-profile
headwater/mainstem dispersion is greater in 68.6% of balanced draws, and
replacement ratios stay close to one. Draw fractions are descriptive.

The RC hypothesis has the expected direction but weak inferential support:
headwater–mainstem mean difference 0.1238, two-sided P = 0.2143,
one-sided P = 0.0579. Intermediate and mainstem RC means are identical.
Of 946 total site pairs, **867 saturate at −1 and 14 at +1**. All 91
intermediate-within-group pairs and all 15 mainstem-within-group pairs are −1.
Such saturation means the null provides little resolution in those groups;
−1 does not mean identical observed molecular profiles. Headwater RC exceeds
mainstem in 84.3% of balanced draws and ties it in the remainder.

Replacement's Lingoes constant is 0.5391; the other constants are zero or
numerical roundoff. The correction substantially changes replacement distance
geometry, so its near-one dispersion ratios deserve limited chemical
interpretation. Mean within-group geographic separation is 5.74 km for
headwaters, 3.14 km for intermediate sites, and 3.13 km for mainstem sites.
Higher headwater dispersion can therefore reflect wider spatial extent as well
as tributary heterogeneity. Segment balancing does not remove that confounding.

The result supports further investigation of heterogeneity after equal-feature
thinning, rather than a general claim that all molecular or structural measures
become less variable downstream. This signal is clearer than the earlier
thinned-profile mean-composition association with drainage (P = 0.1339), but the
two tests have different hypotheses and P values do not measure comparable
biological effect strength. Neither the complete profiles nor five-property
variation nor RC provides significant omnibus support. Detection-count
standardization remains a sensitivity because measured intensities and extraction
settings cannot establish equal analytical sampling depth here.

## Outputs and draft captions

The final script completed without analytical warnings. All six PDFs were
rendered with Poppler through `pdftools` and visually inspected. The 1,000
balanced draws each contain six distinct segments in each of the three groups;
all five omnibus tests use 44 sites. All four inventoried raw-input checksums
still match. No raw inputs, remote files, commits, or pushes were changed.

Tables are under `results/fticr_2016/network_dispersion/`. Important files are
`three_group_tests.csv`, `pairwise_group_contrasts.csv`,
`headwater_mainstem_tests.csv`, `balanced_segment_summary.csv`,
`site_dispersion.csv`, `site_network_position_audit.csv`,
`raup_crick_pairwise_diagnostics.csv`, `null_model_audit.csv`,
`distance_correction_audit.csv`, `molecular_network_pairs.csv`,
and the complete balanced-draw site selections. All figures are vector PDFs.

| Figure in `figures/` | Draft caption |
| --- | --- |
| `2016_fticr_network_groups.pdf` | Locations of the 44 sediment molecular profiles on the LiDAR-derived stream network. Green: first/second order; orange: third/fourth order; blue: fifth-order mainstem. Centerlines are simplified for plotting only. |
| `2016_fticr_dispersion_jaccard.pdf` | Bias-adjusted distances to group spatial medians for binary Jaccard composition. Points represent sites: 24 headwater, 14 intermediate, six mainstem. Box plots describe site distributions, not uncertainty in group means. |
| `2016_fticr_dispersion_equal_features.pdf` | Bias-adjusted distances to group spatial medians for the mean Jaccard matrix from 1,000 equal-feature draws at 391 detections per site. This is a detection-count sensitivity, not measurement-depth standardization. |
| `2016_fticr_dispersion_chemical_properties.pdf` | Bias-adjusted distances to group spatial medians in the full five-variable standardized FT-ICR property space. Formula summaries give equal weight to each detected feature within a site. |
| `2016_fticr_network_raup_crick.pdf` | Within-group richness-constrained RC under a common occurrence-weighted regional pool. Each point is a pair, so points share sites and are dependent. Group inference permutes site labels. Zero is the center of the null expectation, not zero observed dissimilarity. |
| `2016_fticr_balanced_dispersion.pdf` | All three group dispersion ratios from 1,000 draws of six headwater and six intermediate sites on distinct segments, against all six mainstem sites. Values above one indicate greater dispersion in the first named group. These are sensitivity distributions, not confidence intervals. Replacement uses a Lingoes-corrected distance space. |
