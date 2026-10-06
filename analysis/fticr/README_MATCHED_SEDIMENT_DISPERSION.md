# Molecular and microbial dispersion on matched sediment sites

Run `source(here::here("analysis", "fticr", "07_matched_sediment_dispersion.R"))`
in R or Positron. This analysis refits both molecular and microbial dispersion
on their exact shared site set, preserving all full-site results and figures.

## Matched design

The intersection contains 33 sites: 18 first/second-order headwaters on 12
mapped segments, 10 third/fourth-order intermediate sites on nine segments,
and five fifth-order mainstem sites on five segments. Both composition tables
have identical site identities, ordering, group labels, and segment identifiers.
The audit compares coordinates, segments, and orders between microbial metadata
and the previously audited molecular network positions without renaming sites.

Site 43 is excluded because its retained 10K sediment microbial sample has no
molecular profile. Eleven molecular sites have no sediment sample retained at
10K: WS1-2, WS1-8, WS3-4, WS3-5, 47, 54, 58, 65, 68, 170, and 360.
The full union and every inclusion/exclusion are exported in
`site_inclusion_audit.csv`. No raw inputs are altered.

## Four complementary composition views

Primary views are (1) molecular mean equal-feature Jaccard dissimilarity from
1,000 completed draws at 391 detected features per site and (2) microbial
Bray-Curtis dissimilarity from the fixed 10,000-read ASV table. Complete binary
molecular Jaccard and binary microbial Jaccard are sensitivity views. Removing
all-zero columns within the matched subset does not change these distances.

The minimum molecular feature count remains 391 on the matched sites.
Equal-feature thinning is performed independently within each site; excluding
other sites does not change a retained pair's completed mean distance. Therefore
the established mean-distance matrix is subset by exact site IDs, without
repeating the feature draws or changing the thinning target. Microbial
rarefaction is also retained, not repeated. Equalizing molecular feature count
and microbial read count standardizes different quantities; the two primary
dispersion scales are not directly interchangeable.

For every view, reconstruct principal coordinates on the full 33-site matrix,
estimate its Lingoes correction, and refit group spatial medians with
`vegan::betadisper(type = "median", bias.adjust = TRUE)`. These are new
matched-site fits, not previously calculated distances to the 44-site or
34-sample group medians. Methods follow
[Anderson et al. (2006)](https://doi.org/10.1111/j.1461-0248.2006.00926.x) and
the [vegan documentation](https://vegandevs.github.io/vegan/reference/betadisper.html).

All four tests use the same 9,999 residual permutations, with seed 20261015.
Apply BH correction across the two primary omnibus tests, separately from
the two sensitivity omnibus tests. Apply BH across the six secondary pairwise
contrasts for the primary views and separately across the six sensitivity
contrasts. This correction is specific to the matched-site comparison, not
pooled with the earlier full-site analysis. Group tests do not formally test
an ordered trend or a difference between molecular and microbial trends.
Unrestricted residual permutations remain exploratory because spatial-network
exchangeability is not established.

## Shared segment-balanced draws

Each of 1,000 draws selects five distinct headwater segments and five distinct
intermediate segments, followed by one site per segment. All five mainstem
sites are retained. The resulting 15 sites are used for every molecular and
microbial view within that draw. This keeps both sample sets identical in
every sensitivity fit, not only in the full matched analysis. Per-draw seeds
are 20261016 + run number; selected sites and seeds are exported once for all
four views. Retain the full 33-site correction in each subset, refit group
medians, and retain any failures or warnings as diagnostic rows.

Ratios and fractions of draws above one describe sensitivity to sampled
segments; they are not confidence intervals, permutation P values, or
independent replicates. Matching and balancing remove site-membership and
sample-count discrepancies, but do not remove geographic or catchment
confounding, mainstem connectivity, or detection effects.

## Matched results verified on 2026-10-04

The shared subset contains 17,764 microbial ASVs and 4,496 detected molecular
features. All four distance matrices are 33 by 33 with identical site order.
Lingoes corrections are zero or numerical roundoff (maximum distance change
5.6e-17), so the correction does not materially change the geometry.

| View | Headwater mean dispersion | Intermediate mean | Mainstem mean | Headwater/mainstem ratio | Omnibus P | BH-adjusted P |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Molecular equal-feature Jaccard | 0.6422 | 0.6217 | 0.5701 | 1.126 | 0.0004 | 0.0008 |
| Microbial Bray-Curtis | 0.4481 | 0.4368 | 0.4832 | 0.927 | 0.5174 | 0.5174 |
| Complete molecular Jaccard | 0.3733 | 0.3408 | 0.3565 | 1.047 | 0.7312 | 0.7312 |
| Binary microbial Jaccard | 0.5268 | 0.5286 | 0.5459 | 0.965 | 0.6946 | 0.7312 |

The molecular equal-feature headwater/mainstem contrast has raw P = 0.0024
and adjusted P = 0.0144 across six primary-view contrasts. Neither adjacent
molecular contrast survives that correction (headwater/intermediate adjusted
P = 0.0759; intermediate/mainstem adjusted P = 0.1080). All microbial
contrasts are nonsignificant. Mean molecular dispersion is ordered across
the three groups, but the pairwise results do not establish a decline at
each successive transition.

| View | Median balanced headwater/mainstem ratio | Draws with headwater > mainstem |
| --- | ---: | ---: |
| Molecular equal-feature Jaccard | 1.122 | 100.0% |
| Microbial Bray-Curtis | 0.919 | 2.1% |
| Complete molecular Jaccard | 1.055 | 65.9% |
| Binary microbial Jaccard | 0.962 | 5.2% |

All equal-feature draws retain higher headwater and intermediate dispersion
than mainstem; headwater exceeds intermediate in 92.3% of draws. Microbial
Bray-Curtis mainstem dispersion exceeds headwater dispersion in 97.9% of
shared draws, but its omnibus test does not support a group difference.
Thus, aligning sites preserves the molecular equal-feature pattern without
revealing a corresponding microbial decline. Complete molecular profiles do
not show significant group dispersion differences. A significant molecular
test and nonsignificant microbial test do not themselves establish a
statistically significant difference between the two responses. Interpret
this as contrasting observed patterns, not a tested response-by-network
interaction or evidence of a causal consolidation mechanism.

All 4,000 balanced fits succeeded without diagnostic flags. The resampling
loop completed locally in 21.6 seconds; the representative four-view draw
took approximately 0.04 seconds. Every draw has five distinct segments per
group, and the same sites enter all four views. Three new vector PDFs were
rendered with Poppler through `pdftools` and visually inspected. All recorded
input and original network-dispersion output checksums remained unchanged.
No remote action, commit, push, or beta-NTI calculation was performed.

## Outputs and draft caption

Tables are in `results/fticr_2016/matched_sediment_dispersion/`: site inclusion
audit, matched manifest, group counts, community dimensions, four aligned
distance matrices, site dispersion, omnibus tests, secondary contrasts,
distance corrections, test and resampling designs, balanced fits/selections,
input/output checksums, settings, and session information.

Three vector PDFs are in `figures/`:

- `2016_matched_sediment_dispersion_molecular_equal_features.pdf`
- `2016_matched_sediment_dispersion_microbial_bray.pdf`
- `2016_matched_sediment_dispersion_comparison.pdf` combines the two panels.

Draft caption: Molecular and microbial compositional dispersion at the same
33 sediment sites across first/second-order headwaters, third/fourth-order
intermediate sites, and fifth-order mainstem sites (18, 10, and five sites).
Points are bias-adjusted distances to group spatial medians, independently
refitted for molecular mean equal-feature Jaccard dissimilarity and microbial
Bray-Curtis dissimilarity. Molecular distances average 1,000 draws of 391
detected features per site; microbial distances use the fixed 10,000-read ASV
table. Boxes describe site distributions, not uncertainty in means. The two
panels have separate y-axis scales. Molecular signatures are formula-level
features, not confirmed metabolites or structures; a downstream dispersion
pattern does not establish transport or consolidation as its cause.

For manuscript comparisons, use these matched panels together. The original
44-site molecular and 34-sample sediment microbial panels remain useful as
full-sampling summaries, with their respective original tests preserved.

## Two-panel contrast figure added on 2026-10-06

Run `source(here::here("analysis", "fticr", "12_plot_matched_dispersion.R"))`
to draw the completed primary results without repeating inference. The script
checks the original input/output hashes, shared site identities and group labels,
and agreement of plotted means with the saved test table. It retains the older
figures and creates:

- `figures/2016_sediment_microbes_fticr_dispersion_comparison.pdf`
- `figures/2016_sediment_microbes_fticr_dispersion_microbial.pdf`
- `figures/2016_sediment_microbes_fticr_dispersion_molecular.pdf`

Plot data, group means, and source hashes are in the `figure_comparison/`
subdirectory of the matched results. The script runs from within the repository
without arguments. All three new PDFs were rendered and visually inspected;
the existing analysis tables and recorded inputs remained unchanged.

Draft caption: Contrasting microbial and molecular compositional dispersion
across network groups at the same 33 sediment sites. Panel A shows microbial
Bray-Curtis dispersion from fixed-10,000-read communities; panel B shows FT-ICR
molecular dispersion from mean Jaccard dissimilarity across 1,000 draws of 391
detected features per site. Points are bias-adjusted distances to each group's
spatial median, following [Anderson et al. (2006)](https://doi.org/10.1111/j.1461-0248.2006.00926.x)
and the [vegan implementation](https://vegandevs.github.io/vegan/reference/betadisper.html).
Boxes show the median and interquartile range, whiskers extend to the most
extreme observation within 1.5 interquartile ranges, and black diamonds show
group means. All sites are plotted, with horizontal jitter only. Headwaters,
intermediate reaches, and mainstem reaches comprise stream orders 1-2, 3-4,
and 5, respectively, with 18, 10, and five sites. The panels retain separate
y-axis scales because the two distance metrics are not directly interchangeable.
Mean molecular dispersion declines across groups (0.642, 0.622, 0.570;
omnibus BH-adjusted P = 0.0008), whereas microbial means are 0.448, 0.437,
and 0.483 (adjusted P = 0.5174). The contrast describes the observed patterns;
it does not test a difference between microbial and molecular responses, nor
does it establish a significant change at every successive network transition.
