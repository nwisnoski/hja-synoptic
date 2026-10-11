# Sediment composition and characteristics

Completed locally October 10, 2026, with
[`15_sediment_characteristics.R`](15_sediment_characteristics.R).
Run directly with `Rscript --vanilla analysis/fticr/15_sediment_characteristics.R`.

## Sample retention

Each question uses every eligible sample. Missing texture or conductivity
never excludes a sample from unrelated analyses. Raw inputs and the fixed
10K microbial table remain unchanged. Sequencing plate is not a predictor.

| Comparison | Microbial samples | Molecular sites |
| --- | ---: | ---: |
| Organic content | 34 | 44 |
| Four-enzyme profile | 34 | 44 |
| Hydraulic conductivity | 33 | 41 |
| Median grain size | 14 | 17 |
| Organic content conditional on enzymes, or the reverse | 34 | 44 |

Microbial site 43 remains in all eligible microbial analyses, including iCAMP;
it is absent only from comparisons requiring a molecular profile. Site 56
is missing conductivity. Direct microbial–molecular models use all 33 paired
sites, with a separate 32-site conductivity sensitivity. The 45-site union
and every model's inclusion flags are exported. No imputation is used.

## Methods and null hypotheses

Microbial composition is Hellinger-transformed fixed-10K abundance. Molecular
composition is complete presence/absence Jaccard dissimilarity, with the
previously completed mean distances from 1,000 uniform draws of 391 detected
features as a sensitivity. The latter is not calibrated analytical-depth
rarefaction. Full-dimensional PCoA coordinates preserve distances, using
Lingoes correction when needed. The largest reconstruction error is below
1.4e-14. RDA of these coordinates preserves the selected Euclidean geometry;
no two-axis response truncation is used.

Organic content, conductivity, median grain size, and drainage area are
log10-transformed. Each enzyme activity is transformed as
`asinh(activity / median_positive_activity)` and standardized before PCA.
This retains measured zeros, including NAG zeros. Two enzyme PCs explain
93.63% and 93.19% of enzyme variance in the full microbial and molecular
sampling pools. Axes are estimated without microbial or molecular responses;
scales, loadings, variance, and scores are exported. The pools have separate
enzyme references and are not treated as identical replicated experiments.

The three primary predictors are organic content, the two-axis enzyme
profile, and conductivity, fitted separately. Microbial models condition on
log10 drainage area; molecular models additionally condition on **untransformed
detected-feature count**, consistent with the earlier molecular-property
partial RDA. Neither missing molecular data nor molecular detection count
restricts microbial-only models. The null for each test is no additional
composition association with the tested predictor after its stated conditions.
These are marginal associations beyond those conditions, not independent
effects after every sediment variable has been controlled.

There are 9,999 reduced-model residual permutations per test. BH correction
applies jointly to the six primary tests (three predictors, two response types).
The three equal-feature sensitivities form their own correction family.
The two grain-size tests and four organic/enzyme conditional tests have
separate secondary families, with separate equal-feature sensitivities.
No variable selection or individual-ASV significance tests are performed.
Negative adjusted fractions are retained.

Within-stream-order residual permutations are reported as a sensitivity, with
BH correction within the same families. Leave-one-segment-out fits keep the
parent response geometry and full-pool enzyme reference fixed. All 512 fits
succeeded. These checks assess stream-order dependence and segment influence;
they do not establish spatial exchangeability or eliminate geographic
autocorrelation. P values remain exploratory for this observational design.

The implementation follows [vegan's permutation-test documentation](https://vegandevs.github.io/vegan/reference/anova.cca.html)
and [variation-partitioning guidance](https://vegandevs.github.io/vegan/reference/varpart.html).
Hellinger transformation follows
[Legendre and Gallagher (2001)](https://doi.org/10.1007/s004420100716).

## Results

| Predictor | Microbial adjusted variation | BH P | Molecular adjusted variation | BH P |
| --- | ---: | ---: | ---: | ---: |
| Organic content | 3.99% | 0.0102 | 0.52% | 0.2858 |
| Enzyme profile | 0.41% | 0.3835 | 0.98% | 0.2858 |
| Conductivity | 0.27% | 0.3835 | -0.61% | 0.8948 |

The organic-content microbial association remains supported with within-order
permutations (BH P = 0.0084). Segment-deletion adjusted fractions range from
2.68% to 4.80%. Conditioning additionally on the enzyme profile gives 4.50%
(secondary-family BH P = 0.0044). These are modest composition associations,
not evidence that organic content causes the observed community differences.

Grain-size tests are nonsignificant (BH P = 0.4807 for both response types).
All three main equal-feature molecular tests are nonsignificant after BH
correction, including organic content (adjusted fraction 0.16%; BH P = 0.0729;
within-order BH P = 0.2658). The present models do not support a general
conductivity, texture, or enzyme association with both datasets. Enzyme
measurements describe biochemical activity potential and cannot identify
in situ reaction rates, microbial producers, or transformations of particular
molecular features.

## Direct microbial–molecular comparison

The earlier two molecular-property PCs retain their fixed 42-site reference.
The new baseline reproduces its adjusted 4.04% microbial fraction to within
2.4e-15. Monte Carlo P values differ slightly because permutation seeds and
the current test family differ from the earlier analysis; old outputs are
preserved.

| Conditions | Paired sites | Adjusted molecular-property fraction | Raw P | BH P | Within-order raw P |
| --- | ---: | ---: | ---: | ---: | ---: |
| Drainage and molecular detection count | 33 | 4.04% | 0.0120 | 0.0360 | 0.0344 |
| Above plus organic content and enzyme PCs | 33 | 2.44% | 0.0447 | 0.0622 | 0.1000 |
| Above plus conductivity | 32 | 2.34% | 0.0622 | 0.0622 | 0.1477 |

The three nested comparisons form one BH family. The molecular-property
association weakens after adjustment for measured sediment characteristics;
its adjusted test does not pass this family's BH threshold or the within-order
sensitivity. A fraction that declines after conditioning does not identify
mediation, transport, or a causal pathway. Direct profile matching is the only
reason these models omit site 43; the full microbial analyses retain it.

## Outputs and figure captions

Tables and diagnostics are in
[`results/fticr_2016/sediment_characteristics/`](../../results/fticr_2016/sediment_characteristics/).
Constructed site and model inputs are in
[`data/derived/sediment_characteristics_2016/`](../../data/derived/sediment_characteristics_2016/).
All seven input checksums remained unchanged. The new PDFs were rendered and
visually inspected after the final run.

- `figures/supplement/2016_sediment_characteristics_coverage.pdf`: Availability of microbial
  libraries, molecular profiles, and measured sediment characteristics across
  the 45-site union. Missing profiles are distinguished from missing predictors.
- `figures/supplement/2016_sediment_characteristics_unique_fractions.pdf`: Adjusted
  composition variation beyond drainage area (microbes) or drainage area and
  molecular detection count (molecules). Each characteristic uses all available
  samples, with counts printed beside points. Predictors are fitted separately;
  fractions are not mutually exclusive, and response geometries, sample sets,
  and conditions differ. This figure does not test whether effect sizes differ
  between microbes and molecules. Equal-feature molecular results are sensitivity
  estimates. All P values, including nonsignificant tests, are in the tables.
