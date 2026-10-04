# Aquatic microbial composition and chemical signatures

Run `source(here::here("analysis", "fticr", "06_aquatic_microbe_chemistry.R"))`
in R or Positron after completing the earlier optical and chemistry preparation.
The script uses the fixed 10K microbial table and preserves raw files.

## Comparisons and sample overlap

The aquatic comparison matches **20 planktonic microbial samples with surface
water optics** and **31 hyporheic microbial samples with hyporheic water optics**.
Every retained sample has complete SUVA254, fluorescence index, and EEM peak T:C.
These water properties describe DOM optical signatures. They are not molecular
formula inventories, source percentages, or identified algal metabolites.

The existing FT-ICR dataset is treated as sediment/site biogeochemistry under the
project design. Its workbook is indexed by site; separate surface-water and
hyporheic-water formula inventories have not been identified. Consequently,
FT-ICR comparisons provide **cross-compartment site context**, while water
optical comparisons match the microbial habitat. Neither comparison establishes
microbial production or consumption of the detected chemical signatures.

Use the full 60-profile FT-ICR table rather than restrict this task to the earlier
44-site sediment multiblock. There are **19 planktonic and 30 hyporheic FT-ICR
matches**. Planktonic site 67 and hyporheic site 43 lack FT-ICR profiles; they
remain in the water-optics analyses. Each habitat is modeled separately, so a
site shared by the two habitats does not become two independent replicates in
a pooled model. All stream orders, including intermediate sites, are retained.

The new constructed inputs are `data/derived/fticr_2016/site_60_presence_absence.csv`
and `site_60_formula_properties.csv`. The binary table has 60 rows and 4,760
prepared primary features, preserving the established filter from the source
workbook. Forty-two of these features were undetected in the earlier 44-site
subset; they are retained for the expanded site set. Signal only defines
detection at greater than zero; it does not weight ecological comparisons.

## Methods

Microbial responses are Hellinger-transformed abundances from fixed 10,000-read
samples. This transformation supports Euclidean PCA and RDA
([Legendre and Gallagher 2001](https://doi.org/10.1007/s004420100716)). Unconstrained
PCA figures describe the habitat-specific data before conditioning. They do not
display the variation uniquely associated with chemistry, and the first two
axes do not represent all community variation.

Four small models are declared before the workflow:

1. Planktonic composition against the three surface-water optical variables,
   conditioning on log10 drainage area.
2. Hyporheic composition against the three hyporheic-water optical variables,
   with the same conditions.
3. Planktonic composition against two fixed FT-ICR property PCs, additionally
   conditioning on detected-feature count.
4. Hyporheic composition against those same FT-ICR PCs, with the same conditions.

The optical predictors are SUVA254, fluorescence index, and log10 peak T:C,
standardized using the earlier 59-site surface-water reference. All three enter
the models; optical PC1 is used only for figure coloring. That PC is projected
onto hyporheic optics using the same reference, and both figures share the
same color scale. See `README_DOM_INTEGRATION.md` for the source index definitions
and the limits of terrestrial/humic versus aquatic/protein-like interpretations.

The five FT-ICR properties are mean detected H:C, mean detected O:C, N-bearing
fraction, protein-like fraction, and tannin-like fraction. Every detection
contributes equally within a site. Reuse the **42-site reference centering,
scaling, and two FT-ICR-only loading vectors** from `04_integrated_chemistry_pca.R`.
Project all 60 profiles onto those vectors rather than fitting new axes using
the aquatic microbial outcomes. Recomputed traits match the earlier 44-site
summaries within 2.22e-16; projected scores match the 42-site reference within
2.84e-15. The first two reference axes represent 95.3% of the five-property
variance in their training set; that percentage is not recalculated for each
new microbial subset. Formula classes remain signatures rather than confirmed
structures, sources, or metabolites.

Partial RDA uses all microbial PCA coordinates for efficient permutation tests.
No response axes are discarded: pairwise Hellinger distances, constrained and
conditioned inertia, and adjusted fractions are checked against the complete
ASV response. Maximum inertia difference is 4.12e-12; all checks pass a 1e-10
tolerance. Models retain 15/26 residual degrees of freedom for planktonic/
hyporheic optics, and 14/25 for their respective FT-ICR comparisons.

There are 9,999 reduced-model residual permutations for each partial RDA, with
seed 20261011 + model ID. Report raw added fractions and adjusted unique fractions
of **total** community variation. Negative adjusted fractions are retained:
they indicate an effect estimate no larger than expected after the complexity
adjustment, not negative biological variation. Correct the four global model
P values together using BH. The partial-RDA implementation follows
[vegan's documentation](https://vegandevs.github.io/vegan/reference/cca.html).

As a network sensitivity, permute sites only within their exact stream orders
(seeds 20261011 + 100 + model ID); report a separate four-model BH family.
Restriction does not remove spatial dependence among connected sites.
Sequencing plate is excluded from ecological predictors and conditioning
covariates per the user's 2026-10-03 clarification; it remains QC metadata.
Habitat sample sets and
sample sizes differ, so significance in one habitat and not another would not
by itself establish a habitat difference in chemical coupling.

Whole-distance concordance compares microbial Hellinger distance with either
three-variable matched-water optical distance or contextual sediment **binary
Jaccard**. Spearman Mantel tests permute sites, not pairs, using the same site
permutation matrices. These are vegan's one-sided tests for positive concordance.
The four distance tests have their own BH correction; stream-order-restricted tests are an additional four-test sensitivity family.
[Vegan describes this site-permutation test](https://vegandevs.github.io/vegan/reference/mantel.html).
Jaccard concordance is descriptive and remains susceptible to molecular
detection-count variation; it is not the count-conditioned property model.
No new equal-feature thinning is claimed for the expanded 60-site set.

Two further diagnostics are declared: rank-transformed water-optics models and
leave-one-microbial-site-out adjusted fractions for all four models. Ranking
occurs within each full available water compartment before microbial joining,
reducing sensitivity to unusually large or small optical values. Rank tests
have a separate two-test BH correction. Leave-one-site-out fractions assess
influence, not predictive validation or confidence intervals. No individual
taxon/formula screening, model selection, or causal process attribution is done.

## Results on 2026-10-03

| Microbial habitat | Chemistry view | Sites | Adjusted unique variation | Raw P | BH P |
| --- | --- | ---: | ---: | ---: | ---: |
| Planktonic | Matched surface optics | 20 | 4.30% | 0.1101 | 0.3192 |
| Hyporheic | Matched hyporheic optics | 31 | 0.49% | 0.1810 | 0.3192 |
| Planktonic | Contextual sediment FT-ICR PCs | 19 | 1.74% | 0.2394 | 0.3192 |
| Hyporheic | Contextual sediment FT-ICR PCs | 30 | −0.14% | 0.5642 | 0.5642 |

The matched-water optical properties have a modest estimated association in
planktonic samples and a small one in hyporheic samples; neither global test
survives correction. The rank sensitivity gives 4.42% in planktonic samples
(BH P = 0.1093) and 0.76% in hyporheic samples (BH P = 0.1093), likewise not
confirming an association under these tests.
None of the four order-restricted RDA tests survives correction either.

| Microbial habitat | Signature distance | Mantel rho | Raw P | BH P |
| --- | --- | ---: | ---: | ---: |
| Planktonic | Matched surface optics | 0.127 | 0.1012 | 0.2024 |
| Hyporheic | Matched hyporheic optics | −0.030 | 0.5845 | 0.5845 |
| Planktonic | Contextual sediment Jaccard | 0.279 | 0.0223 | 0.0892 |
| Hyporheic | Contextual sediment Jaccard | 0.059 | 0.2805 | 0.3740 |

The strongest distance association is planktonic composition versus sediment
formula incidence, but its nominal P = 0.0223 does not survive correction and
does not control detection count. The corresponding count-conditioned property
model is also not significant. Its leave-one-site-out adjusted fractions range
from **−0.90% to 3.55%**, compared with a full estimate of 1.74%; this is a
tentative site-level association rather than a stable result. The other model
ranges are 2.79% to 9.08% for planktonic water optics, −0.12% to 1.69% for
hyporheic water optics, and −0.39% to 0.07% for hyporheic FT-ICR properties.
No order-restricted Mantel test survives its four-test correction.

These comparisons do not establish robust coupling between aquatic microbial
composition and the measured water DOM properties. They also do not contradict
the earlier modest sediment-microbial/FT-ICR-property association: that analysis
has a different microbial habitat and sample set. Sample sizes are small and
the chemical measures cover only selected properties of DOM.

## Tables, figures, and verification

Inspect tables under `results/fticr_2016/aquatic_microbe_chemistry/`, especially
`habitat_overlap_audit.csv`, `aquatic_site_chemistry_matches.csv`,
`aquatic_chemistry_partial_rda.csv`, `aquatic_chemistry_mantel.csv`,
`water_optics_rank_sensitivity.csv`, and `leave_one_site_out_summary.csv`.
Detailed influence rows and coordinate/reference checks are retained separately.
All three individual panels use vector PDF and a white background.

| Figure in `figures/` | Draft caption |
| --- | --- |
| `2016_planktonic_microbes_water_optics.pdf` | Hellinger PCA of 20 planktonic communities at fixed 10K depth, colored by matched surface-water optical PC1. Shapes retain headwater, intermediate, and mainstem groups. Microbial PC1 and PC2 explain 24.4% and 7.4% of variation. The unconditioned display is descriptive; chemistry models condition on drainage. |
| `2016_hyporheic_microbes_water_optics.pdf` | Hellinger PCA of 31 hyporheic communities at fixed 10K depth, colored by matched hyporheic-water optical PC1 projected onto the surface reference. The color scale is shared with the planktonic figure. Microbial PC1 and PC2 explain 8.9% and 4.9%; models use the full community response. |
| `2016_aquatic_microbes_chemistry_associations.pdf` | Adjusted unique fractions of total microbial variation from four partial RDA models. Matching water optics use all three optical predictors; contextual sediment FT-ICR uses two fixed property axes and adds detection count as a covariate. Negative adjusted values are retained. Points are effect estimates, not uncertainty intervals. No global model survives four-test BH correction. |

The final script completed without analytical warnings. All three PDFs were
rendered with Poppler through `pdftools` and visually inspected after the
final run, retaining the shared color scale. Exact habitat/site joins,
10K response depths, projection
and response-geometry checks, and all four raw-source checksums passed.
`git diff --check` passed. No raw inputs, commits, pushes, or remote files changed.
