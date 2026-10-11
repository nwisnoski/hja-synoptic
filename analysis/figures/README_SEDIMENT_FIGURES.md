# Sediment-focused manuscript figures

Generated October 10, 2026, by
[`01_sediment_manuscript.R`](01_sediment_manuscript.R). Run directly with:

```sh
Rscript --vanilla analysis/figures/01_sediment_manuscript.R
```

The script reformats completed results and performs no new hypothesis tests.
Original figures and analysis tables are preserved. All new PDFs have the
prefix `2016_sediment_ms_`. The selected matched-dispersion composite is in
`figures/manuscript/`; the other twenty PDFs are in `figures/supplement/`. Twenty-one
PDFs include four candidate main figures, their individual panels, and
supporting views. These are a manuscript review set, not a journal submission.

## Candidate main figures and captions

### Figure 1. Sediment sampling across the stream network

[`2016_sediment_ms_sampling_map.pdf`](../../figures/supplement/2016_sediment_ms_sampling_map.pdf)

Sediment sampling locations across the H.J. Andrews stream network. Colors
indicate first/second-order headwaters, third/fourth-order intermediate reaches,
and fifth-order mainstem reaches. Circles identify 33 sites with both retained
microbial and molecular profiles, the triangle identifies microbial-only site
43, and squares identify 11 molecular-only sites. The union contains 45 sites,
including all 34 microbial samples retained at 10,000 reads and all 44 sediment
molecular profiles. Coordinates and stream orders are reused from completed
centerline audits; site 43 is retained in microbial analyses despite lacking a
molecular profile. Nearby sites may overlap at catchment scale. Profile
availability differs from absence of individual ASVs or molecular features.

### Figure 2. Sediment composition and molecular signatures

[`2016_sediment_ms_composition.pdf`](../../figures/supplement/2016_sediment_ms_composition.pdf)

(A) Hellinger PCA of all 34 fixed-10K sediment microbial samples. Axes explain
16.9% and 13.6% of microbial composition variation and are the saved sediment
PCA coordinates. Colors identify network positions. (B) Mean fractions of
detected molecular features in workbook chemical classes across all 44
sediment profiles, with 24 headwater, 14 intermediate, and six mainstem sites.
Each site receives equal weight in its group mean. Formula-based class labels
are chemical signatures rather than confirmed structures or biological
sources. The two panels describe their full available sampling pools; they
do not test cross-dataset correspondence or differences among network groups.

### Figure 3. Paired microbial and molecular dispersion

[`2016_sediment_ms_dispersion_paired.pdf`](../../figures/manuscript/2016_sediment_ms_dispersion_paired.pdf)

Microbial and molecular dispersion across the same 33 sediment sites:
18 headwater, 10 intermediate, and five mainstem. (A) Microbial Bray–Curtis
distance to the network-group spatial median. (B) Molecular dispersion using
the mean Jaccard distances from 1,000 uniform draws of 391 detected features
per profile. Points are individual sites, diamonds are group means, and boxes
show medians and interquartile ranges, with whiskers extending to observations
within 1.5 interquartile ranges. Existing bias-adjusted dispersion and
permutation results are reused without refitting. Microbial omnibus P = 0.5174;
the equal-feature molecular test has raw P = 0.0004 and BH P = 0.0008 across
the two primary matched views. Complete molecular Jaccard and binary microbial
Jaccard sensitivities remain nonsignificant. Contrasting significance does not
test a difference between responses. Equal-feature thinning is a detection
sensitivity, not calibrated analytical-depth rarefaction.

### Figure 4. Sediment conditions and microbial–molecular correspondence

[`2016_sediment_ms_characteristics.pdf`](../../figures/supplement/2016_sediment_ms_characteristics.pdf)

(A) Adjusted composition variation associated separately with organic content,
the two-axis enzyme profile, or hydraulic conductivity. Microbial models
condition on log10 drainage area; molecular models additionally condition on
detected-feature count. All eligible samples are used for each question:
34 microbial samples for organic content and enzymes, 33 for conductivity;
44 molecular sites for organic content and enzymes, 41 for conductivity.
Fractions are not mutually exclusive, and sample pools, response geometries,
and conditioning variables differ between response types. Negative adjusted
fractions are retained. Organic content explains 3.99% adjusted microbial
variation (primary-family BH P = 0.0102); the other primary tests are
nonsignificant. (B) Adjusted microbial variation associated with the two fixed
molecular-property PCs at all 33 paired sites, first conditioning on drainage
and detection count, then additionally on organic content and enzyme PCs.
The final conductivity sensitivity uses 32 paired sites. Fractions decline
from 4.04% to 2.44% and 2.34%; the sediment-adjusted models do not survive the
declared BH correction. Points are model estimates; no confidence intervals
are implied. These cross-sectional associations do not establish causes,
mediation, molecular production/consumption, or in situ process rates.
Full tests and sensitivity results are documented in the
[characteristics report](../fticr/README_SEDIMENT_CHARACTERISTICS.md).

## Supporting views

| Figure | Purpose and retained samples |
| --- | --- |
| [`2016_sediment_ms_dispersion_all_samples.pdf`](../../figures/supplement/2016_sediment_ms_dispersion_all_samples.pdf) | Full microbial Bray–Curtis (34), complete molecular Jaccard (44), and equal-feature molecular sensitivity (44). Retains every eligible site; panels are descriptive parallel views with unequal site sets. Use the paired main figure for direct comparisons. |
| [`2016_sediment_ms_diversity.pdf`](../../figures/supplement/2016_sediment_ms_diversity.pdf) | Hill q = 1 across drainage area for all 34 sediment microbial samples, colored by network group. No regression or smoothing is added. |
| [`2016_sediment_ms_unifrac.pdf`](../../figures/supplement/2016_sediment_ms_unifrac.pdf) | Weighted and unweighted UniFrac PCoA refitted on the same 34 sediment samples using saved pairwise distances. Axes describe within-sediment variation; no other habitats enter the ordination. |
| [`2016_sediment_ms_molecular_pcoa.pdf`](../../figures/supplement/2016_sediment_ms_molecular_pcoa.pdf) | Complete incidence-Jaccard composition of all 44 molecular profiles, colored by network group. Interpret with the established strong PCo1/detection-count association. |
| [`2016_sediment_ms_molecular_pcoa_equal.pdf`](../../figures/supplement/2016_sediment_ms_molecular_pcoa_equal.pdf) | Existing equal-feature sensitivity coordinates for all 44 molecular profiles. The first two axes explain only 5.25% and 3.26%; the display is a limited summary of the distance geometry. |
| [`2016_sediment_ms_microbial_pca_organic.pdf`](../../figures/supplement/2016_sediment_ms_microbial_pca_organic.pdf) | The unchanged all-34 microbial PCA colored by measured sediment organic content; the multivariate test is supplied separately, not inferred from this color overlay. |

Individual panels are retained with the same prefix. The
[`figure_manifest.csv`](../../results/sediment_manuscript_figures_2016/figure_manifest.csv)
lists all 21 filenames, dimensions, sampling counts, and intended roles.
Existing stable-group maps and networks remain available as an optional
feature-level figure; their selected spatial interpretation remains exploratory.
The sediment iCAMP figure awaits results from job 432228. Catchment-null
within-sediment pairs are not substituted for the sediment-only null pool.

Existing P values and multiplicity families are unchanged. In particular,
the full-sample microbial tests retain their original across-habitat
corrections; this plotting script does not silently define new test families.
No fitted trends, group ellipses, or new causal interpretations are added.

## QA and reproducibility

All 21 PDFs were rendered and visually inspected. Labels, legends, sample
counts, panel alignment, and original dispersion group means were checked.
Both UniFrac response matrices retain all 34 sediment samples. Weighted
UniFrac uses Lingoes constant 0.00063539 and explains 32.43%/17.06% on its
first two axes; unweighted UniFrac needs only a numerically negligible
correction and explains 13.56%/10.64%. Corrected full-dimensional distances
reconstruct to within 2.7e-15. No pairwise UniFrac distances were recalculated.

Checksums confirmed that all 136 inventoried sources and pre-existing figures
were unchanged. Exported plot data, within-sediment UniFrac coordinates,
ordination diagnostics, and integrity checks are in
[`results/sediment_manuscript_figures_2016/`](../../results/sediment_manuscript_figures_2016/).
New figures are vector PDFs with white backgrounds and consistent network
colors. No remote commands, commits, or pushes were performed during this
figure work; deployed iCAMP scripts were not modified.
