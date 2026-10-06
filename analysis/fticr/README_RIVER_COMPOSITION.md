# River composition and microbial integration

Exploratory local analysis completed on 2026-10-03. Molecular analyses use
presence/absence exclusively; intensities enter detection and QC summaries only.
The executable workflow is `02_river_composition.R`, with all analytical choices
and seeds in its first numbered section. Run `01_explore_sediment.R` first, then:

```r
source(here::here("analysis", "fticr", "02_river_composition.R"))
```

Tables are in `results/fticr_2016/river_composition/`. All seven new figures are
single-page vector PDFs in `figures/`. The script ran successfully from its own
directory, demonstrating repository-relative path handling.

## Published methods and source code

The user's description identifies the basic PNNL workflow. Related publications
support presence/absence, but they do not identify the exact HJA extraction or
acquisition settings. Two relevant examples illustrate why that distinction matters:

- [Graham et al. (2017), with James Stegen](https://www.emsl.pnnl.gov/sites/default/files/2021-06/2.1.Stegen_SuppMaterials_Graham2017.pdf)
  describe sequential water, methanol, and chloroform sediment extraction,
  negative-mode electrospray, and conversion to presence/absence. Ion accumulation
  was optimized by sample, so equal acquisition effort cannot be assumed from a
  processed intensity table. Their solvent sequence is context, not a verified
  description of these HJA samples.
- [Danczak et al. (2023), with James Stegen](https://www.frontiersin.org/journals/water/articles/10.3389/frwa.2023.1087108/full)
  describe aqueous DOM from the Yakima basin, PPL extraction, carbon preparation,
  negative-mode electrospray, Formularity assignments, and binary detections.
  Their aqueous preparation should not be substituted for a sediment protocol.

The following public repositories were inspected on 2026-10-03:

- [Meta-Metabolome_Ecology: FTMS_Analysis.R](https://github.com/danczakre/Meta-Metabolome_Ecology/blob/master/FTMS_Analysis.R)
  processes Formularity columns corresponding to our workbook schema. Its site
  summaries select positive detections and calculate class counts and unweighted
  chemical-property summaries. It demonstrates an incidence-based approach;
  its several classification options do not establish which boundaries generated
  the HJA workbook's existing labels.
- [Meta-Metabolome_Ecology: FTICR_RCBC_PresAbs.R](https://github.com/danczakre/Meta-Metabolome_Ecology/blob/master/FTICR_RCBC_PresAbs.R)
  implements a presence/absence Raup–Crick comparison. Random communities retain
  each site's observed feature count and sample the regional pool with weights
  proportional to feature occurrence. This is a richness-constrained null model,
  rather than thinning all samples to the same count. Its binary Bray–Curtis
  comparison is equivalent to Sørensen dissimilarity. We inspected this code;
  we have not run this null analysis on HJA.
- [Yakima-River-Basin-Functional-Diversity](https://github.com/danczakre/Yakima-River-Basin-Functional-Diversity)
  supplies the manuscript's data and identifies `RC2_Analysis_Final.R` as its
  Rao functional-diversity workflow. The repository README was accessible;
  retrieval of that individual script failed, so its implementation was not verified.

The HJA workbook does not supply the needed extraction, carbon loading, injection,
ion-accumulation, scan, blank, replicate, and detection-threshold records. We retain
the established preparation filter and do not impose a literature S/N threshold
on already processed intensities. Formula ambiguity is also retained rather than
silently imposing a new candidate-count filter.

## Rarefaction: precedent and limits

Molecular feature subsampling has published precedent. For example,
[Freeman et al. (2024)](https://www.nature.com/articles/s41467-023-44431-4)
subsampled study-level compound pools to compare overlap; their
[author manuscript](https://api.repository.cam.ac.uk/server/api/core/bitstreams/f04e362b-fa78-41c2-93a5-b6b5b171132e/content)
describes 1,000 draws of 6,000 compounds per study. That is a pooled comparison,
not evidence that per-site molecular rarefaction is a universal FT-ICR-MS standard.

Sequencing counts represent repeated observations of sequence types. A binary
molecular profile records one detection state per retained feature, with unequal
ionization and detection probabilities. Uniformly drawing detected features
standardizes the number of retained entries, not carbon loading, ion sampling,
or detection completeness. It cannot recover compounds that were not detected.
Total signal is not a calibrated sampling-depth denominator.

We therefore retain the complete binary profiles as the primary descriptive
data and use **equal-feature subsampling as a sensitivity analysis**:

1. Draw 391 observed features without replacement at each of the 44 sites.
2. Repeat 1,000 times using the exported run design and individual seeds.
3. Recalculate Jaccard dissimilarity, drainage R², paired microbial-distance
   correlation, and mean H:C. Retain run failures and error messages if they occur.
4. Average distance matrices across successful runs for one additional PCoA
   and site-permutation comparison.

All 1,000 runs succeeded. Exported 2.5–97.5% ranges describe subsampling variation,
not confidence intervals for an ecological population. Thinning retains only
14.8% of detections at the richest site and therefore discards substantial
information. Uniform thinning preserves the expectation of a site's mean
chemical property or class fraction; it primarily tests distance sensitivity.
We do not estimate standardized richness from these draws, because every draw
has 391 features by construction.

Within-site incidence rarefaction or coverage estimation would require suitable
independent sampling units, such as documented replicate measurements. Treating
different river sites as those replicates would change the question to pooled
regional richness. A richness-constrained null model is another defensible future
comparison, but neither it nor thinning establishes the cause of nondetection.

## Drainage-area PCoA and chemical composition

There are 44 molecular sites spanning 10.2–6,127 ha of drainage area. Primary
PCoA1 correlates almost perfectly with detection count (Spearman rho = −0.984).
This identifies a strong richness component, but does not determine whether
richness differences arise from ecology, analytical detection, or both.

PCoA2 has a clearer chemical interpretation within these data: positive scores
track protein-like detections (rho = 0.908), nitrogen-bearing formulas (0.857),
and mean H:C (0.783). Negative scores track mean O:C (−0.826) and tannin-like
detections (−0.844). These are descriptive associations among quantities derived
from the same molecular matrix, not independent validation of a source gradient.

Site contrasts demonstrate this chemical variation:

| Site | Drainage area (ha) | Mean H:C | Mean O:C | Protein-like detections | Tannin-like detections |
| --- | ---: | ---: | ---: | ---: | ---: |
| WS1-6 | 82.7 | 1.19 | 0.523 | 13.6% | 15.0% |
| CC-4 | 58.5 | 1.20 | 0.510 | 12.9% | 15.6% |
| KC-3 | 114.0 | 1.54 | 0.426 | 32.3% | 3.3% |
| 56 | 99.5 | 1.52 | 0.423 | 30.6% | 3.6% |

The complete keyed site tables allow comparison of individual bars, ordination
coordinates, chemistry, drainage area, and stream order.

There is no convincing monotonic shift in these signatures with drainage area.
None of eight chemical-property correlations with log drainage area survives
within-family BH adjustment. Lignin-like fractions average 36.9–41.9% across
stream orders 1–5, and protein-like fractions average 21.3–26.5%, without a
consistent downstream progression. Only three sites have stream order 1;
21 have order 2. These are descriptive averages in an unbalanced design.

| Molecular comparison | Drainage-area R² | Site-permutation P |
| --- | ---: | ---: |
| Complete-profile Jaccard | 0.0515 | 0.0173 |
| Jaccard replacement component | 0.0230 | 0.6217 |
| Mean Jaccard after equal-feature subsampling | 0.0245 | 0.1339 |

The complete-profile association is small and sensitive to the richness treatment.
Across individual subsamples, median drainage R² is 0.0244, with a subsampling
range of 0.0229–0.0260. Loss of an association after thinning is not proof of
an analytical artifact, because thinning also removes genuine information.

The replacement/nestedness-resultant partition follows
[Baselga (2012)](https://doi.org/10.1111/j.1466-8238.2011.00756.x), implemented
with `adespatial::beta.div.comp(coef = "BJ")`. Replacement accounts for 56.9%
of summed pairwise Jaccard dissimilarity and nestedness-resultant differences for
43.1%. Neither component isolates analytical nondetection. Pairwise shares are
descriptive; the 946 pairs are not independent observations. Replacement PCoA
requires a Lingoes correction and displays only 8.2% on its first two axes.
Correction and axis diagnostics are exported in `sensitivity_pcoa_summary.csv`.

## River interpretation

A decline in aromatic terrestrial signatures and increase in aliphatic signatures
along rivers is a relevant hypothesis, discussed by
[Creed et al. (2015)](https://www.usgs.gov/publications/river-a-chemostat-fresh-perspectives-dissolved-organic-matter-flowing-down-river).
Our sediment data do not currently show a clear version of that gradient.
Higher H:C and protein-like fractions are compatible with comparatively
hydrogen-rich organic matter; oxygen-rich and tannin/lignin-like signatures are
compatible with some terrestrially influenced material. These overlaps do not
provide terrestrial or aquatic source percentages, confirmed carbon structures,
or measurements of microbial availability. Workbook class labels are retained
as “-like” signatures, with their original boundaries unverified.

Drainage area summarizes position across a branched network. It does not trace
one parcel of water from upstream to downstream, and sediment extracts differ
from water-column DOM. Local retention, sediment properties, tributary inputs,
and processing could generate site contrasts, but their contributions remain
untested. Surface-water SUVA254 and fluorescence index do not provide strong
corroboration of these sediment signatures after the exploratory multiple-test
adjustment; the two measurements represent different compartments.

## Molecular–microbial composition

The paired analysis uses the verified 33-site intersection and refits both
ordinations on that intersection. Microbial composition is Hellinger-transformed
from the fixed 10K table. Sample labels and read totals are checked explicitly.

- Complete molecular Jaccard versus microbial Hellinger distance: Spearman
  Mantel rho = 0.112, P = 0.1019 with unrestricted site permutations.
- Two-axis Procrustes correlation = 0.130, P = 0.8281. The two dimensions retain
  40.7% of paired molecular variation and 30.7% of microbial variation, so this
  is a limited projection. The display centers and normalizes both configurations
  to unit inertia and rotates them without subsequently shrinking one.
- Replacement versus microbial distance gives rho = 0.119. P changes from 0.0600
  under unrestricted permutations to 0.0128 within stream order; BH-adjusted
  P is 0.0512 for the latter among the four planned Mantel comparisons. This
  small, permutation-sensitive result does not establish robust concordance.
- Mean equal-feature molecular distance versus microbial distance: rho = 0.0053,
  P = 0.4670. Across individual draws, median rho = −0.0011, with a subsampling
  range of −0.0342–0.0320.
- Mean H:C explains 3.36% of microbial variation in a small marginal PERMANOVA
  accounting for log drainage area and detection count; P = 0.3279. Sequencing
  plate is excluded from ecological models and retained only as QC metadata.

These results provide no robust overall molecular–microbial composition link
under the comparisons examined. They do not demonstrate absence of chemical
effects or rule out specific taxa or formula subsets. Permutations exchange
sites, not individual distance pairs. Within-order permutations probe one
alternative exchangeability scheme; spatial/network dependence remains unresolved.
No association establishes microbial production, consumption, or transport.

The subsequent [DOM and chemistry integration](README_DOM_INTEGRATION.md) tests
chemical-property PCs rather than overall formula distances. It finds a modest
FT-ICR property association with microbial composition after conditioning on
drainage and detection count (adjusted unique fraction 4.04%, P = 0.0103),
while adding surface/hyporheic optics does not strengthen the association.

## Figures and draft captions

| File | Draft caption |
| --- | --- |
| `2016_fticr_site_class_composition.pdf` | Fraction of detected primary features assigned to each source-workbook class at 44 sediment sites, ordered by increasing drainage area. Classes indicate formula-based signatures; percentages are feature-count fractions. |
| `2016_fticr_pcoa_formula_overlay.pdf` | Existing 44-site binary Jaccard PCoA with drainage area shown by color and detection count by point size. Arrows show scaled correlations of unweighted mean H:C and O:C with the displayed axes; they do not represent flow direction or source proportions. |
| `2016_fticr_h_c_drainage.pdf` | Unweighted mean H:C against drainage area at 44 sites, colored by stream order. The descriptive linear fit uses log10 drainage area; shading is its 95% mean-response confidence band. |
| `2016_fticr_protein_like_drainage.pdf` | Protein-like fraction of detected primary features against drainage area at 44 sites, colored by stream order. Workbook labels do not confirm proteins or an aquatic source. |
| `2016_paired_molecular_microbes_procrustes.pdf` | Rotated two-axis molecular Jaccard PCoA and microbial Hellinger PCA at 33 matched sites. Both configurations have unit inertia; circles and triangles indicate molecular and microbial positions, respectively. Lines connect measurements from the same site, and color indicates drainage area. |
| `2016_fticr_replacement_pcoa.pdf` | Lingoes-corrected PCoA of the Jaccard replacement component at 44 molecular sites, colored by drainage area. The displayed axes retain 4.8% and 3.4% of corrected positive eigenvalues. |
| `2016_fticr_equal_feature_pcoa.pdf` | PCoA of mean Jaccard dissimilarity across 1,000 uniform draws of 391 detected features per site. Color indicates drainage area; axes retain 5.2% and 3.3% of positive eigenvalues. This is a sensitivity analysis of observed features, not standardized analytical sampling depth. |

## Network-group chemical-class figures added on 2026-10-06

Run `source(here::here("analysis", "fticr", "14_plot_fticr_network_classes.R"))`
to draw chemical-class comparisons using the established stream-order groups:
headwater orders 1-2 (24 sites), intermediate orders 3-4 (14 sites), and mainstem
order 5 (six sites). These figures use all 44 sediment FT-ICR profiles, rather
than the 33-site microbial intersection. The original drainage-ordered figure
is retained unchanged.

- `2016_fticr_network_group_class_comparison.pdf` shows three stacked bars of
  mean site-level fractions, with each site weighted equally within its group.
- `2016_fticr_site_class_composition_by_network_group.pdf` retains all individual
  site bars in three panels, ordered by increasing drainage area within each
  group. Group labels appear below the x-axis site labels.
- `2016_fticr_site_class_composition_headwater.pdf`,
  `2016_fticr_site_class_composition_intermediate.pdf`, and
  `2016_fticr_site_class_composition_mainstem.pdf` retain the individual panels.

Site-level plot data, group means and descriptive standard deviations, and
source hashes are in `results/fticr_2016/river_composition/network_class_figures/`.
The script checks the audited group assignments, verifies the saved class
fractions against current binary detections and source-workbook labels, and
confirms that the original inputs and figure remain unchanged. All site and
group fraction sums equal one within numerical tolerance. No new significance
tests or regression lines are introduced. All five PDFs were rendered and
visually inspected.

Draft caption for the compact comparison: Mean chemical-class composition of
detected FT-ICR molecular features across headwater, intermediate, and mainstem
sediment sites (24, 14, and six sites). Each site's class counts are divided by
its total detected primary features, then site-level fractions are averaged
within network groups. This gives each site equal weight and avoids weighting
feature-rich sites more heavily. Colors and class assignments follow the
original workbook and earlier site-class figure. Fractions reflect detected
feature counts, not intensity, concentration, or verified compounds. Bars are
descriptive group means; they do not show sampling uncertainty or establish
significant differences between groups. The individual-site counterpart shows
within-group variation using the same fractions. Published methods and the
limitations of the workbook's class assignments are discussed above.

## Verification and next work

Binary entries, row/column alignment, class-fraction sums, matched microbial
labels, and fixed 10K read totals passed checks. The beta partition reproduced
vegan Jaccard and summed correctly to tolerance 1e-12. All subsampling runs
succeeded. PDFs were rendered through R's Poppler-backed `pdftools` and visually
inspected. Seeds, run diagnostics, statistical summaries, site coordinates,
and session information are exported. Raw inputs remain unchanged.

Useful next analyses are a small sediment-property comparison, explicit network
structure where reliable coordinates/connectivity permit it, and a molecular
richness-constrained null sensitivity with a clearly defined regional pool.
Exact extraction and acquisition records remain necessary to distinguish
analytical detection differences from ecological variation. No cluster work
or source-pool/functional attribution is implied by these exploratory results.
