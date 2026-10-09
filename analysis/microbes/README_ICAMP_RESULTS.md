# Catchment iCAMP results and diagnostics

The full catchment job 432103 completed October 9, 2026, after 2 days,
18 hours, and 9 minutes. It used iCAMP 1.9.1, 1,000 randomizations, the fixed
10K community (100 samples, 49,260 ASVs), and 893 phylogenetic bins.
The sediment-only run remains separate and has not been submitted here.

Run `16_icamp_figures.R` directly in Positron to read the completed result,
summarize diagnostics, and generate four PDFs in `figures/`. No inference,
distance construction, or randomization is repeated. Diagnostic CSVs go under
`results/icamp_2016/catchment_castor/diagnostics/`; original exports remain intact.

## What is saved

`icamp_result.rds` contains the final pairwise process fractions and `detail`:

- `SigbMPDi`: 4,950 sample pairs by 893 bins, plus sample identifiers; signed
  confidence scores comparing abundance-weighted betaMPD with its within-bin
  phylogenetic tip-shuffle null.
- `SigBCa`: the same pairs and bins; signed confidence scores comparing bin
  Bray-Curtis turnover with its across-all-taxa randomization null.
- `bin.weight`: mean relative abundance of each bin across the two samples.
- `taxabin`, `setting`, `comm`, and the final pairwise process table.

For ordinary scores, a positive value is the fraction of null distances below
the observed distance; a negative value is minus the fraction above it, using
whichever directional fraction is larger. Exact ties count in neither tail.
The cutoff is strictly greater than +0.975 or less than -0.975. Positive
significant phylogenetic scores imply heterogeneous selection; negative scores
imply homogeneous selection. Compositional scores classify dispersal limitation
or homogenizing dispersal only where phylogenetic selection is not significant.
These scores are not raw betaMPD/Bray distances, z scores, betaNRI/betaNTI, or
Raup-Crick indices. The +/-1.1 phylogenetic scores are special-case corrections,
not probabilities; the plots display them separately from the ordinary range.

The run used `detail.null = FALSE`. Full null-distance distributions and the
ordinary observed bin-level distances were not retained. Observed Bray-Curtis
distances at ASV and pooled-bin resolution are cheaply recalculated from the
saved community for the diagnostic figure; they are not additional null scores.

## First interpretation

Overall process fractions are 42.6% drift/other, 32.1% dispersal limitation,
18.5% homogeneous selection, 5.0% heterogeneous selection, and 1.8% homogenizing
dispersal. Within sediment, homogeneous selection is 31.8%, dispersal limitation
17.2%, and drift/other 47.1%. Homogeneous selection is 37.5% within soil versus
17.1% within planktonic and 10.2% within hyporheic communities. Water habitats
have larger dispersal-limitation fractions (40.7/41.4%). These are mean
abundance-weighted bin classifications, not fractions of taxa or samples.

Between-habitat heterogeneous selection averages 5.9%, versus 2.6% within
habitats; this does not reproduce the strong compositional habitat contrast
as a dominant heterogeneous-selection fraction. Soil-planktonic mean observed
Bray-Curtis is 0.978 at ASV resolution and 0.833 after pooling into bins.
Approximately 48.1% of their bin abundance weight involves a bin absent from
one community; the corresponding soil-hyporheic values are 0.966, 0.788, and
43.8%. Broad turnover among bins is therefore substantial.

The phylogenetic null shuffles tip identities within each bin, keeping each
sample's bin membership and abundance fixed. Broad differences in which bins
occur do not automatically produce heterogeneous selection. If a bin is absent
from one sample, its ordinary betaMPD is zero for every tip shuffle; iCAMP's
special-case procedure can modify some resulting scores. Ecological differences
among close relatives can also escape a phylogenetic turnover test. These are
methodological explanations to investigate, not proof that they explain all of
the discrepancy. Low heterogeneous selection is not evidence against habitat
filtering, and dispersal limitation is not a direct measurement of transport.

Within-bin habitat-affiliation signal is assessed below; signal along continuous
environmental gradients remains unchecked.
All fractions depend on the observed regional pool, binning, tree, and null
model. The drift/other category includes unresolved processes. Pairs sharing
samples are not independent replicates; no new significance tests are presented.

## Figure captions

1. `2016_icamp_catchment_process_fractions.pdf`: mean abundance-weighted process
   fractions for ten within-/between-habitat comparisons. Each sample pair has
   equal weight within its comparison; habitat groups have unequal pair counts.
2. `2016_icamp_catchment_phylogenetic_scores.pdf`: distribution of signed
   within-bin betaMPD null-comparison scores, weighted by pair-specific bin
   abundance. Ordinary scores use 0.05-wide intervals; +/-1.1 special-case
   corrections are shown separately. Dashed lines mark +/-0.975. The color
   scale uses a square-root transformation and reports mean abundance weight.
3. `2016_icamp_catchment_compositional_scores.pdf`: equivalent distribution
   for the Bray-Curtis null-comparison scores. This includes all bins, including
   those already classified as selection; it is not a standalone dispersal
   classification. The process figure applies the sequential classification.
   Compositional scores remain within [-1, 1], so no special-case columns are shown.
4. `2016_icamp_catchment_turnover_resolution.pdf`: mean observed Bray-Curtis
   dissimilarity at ASV and pooled-bin resolution for each habitat comparison.
   Gray connectors join the two resolutions of the same comparison; they are
   not regression fits. Aggregation necessarily removes within-bin turnover.

## Validation and references

### Within-bin habitat-affiliation signal

Run `17_icamp_phylogenetic_signal.R` locally. It uses the saved 893 united
bins and fixed-10K community; it does not change bins or repeat iCAMP inference.
Soils lack the water-chemistry measurements used for aquatic samples, so this
first check tests observed habitat affiliation, not continuous chemistry.

For each ASV, abundance-weighted means of four habitat indicators give its
affiliation profile. The primary profile weights each sample inversely by
the number sampled in its habitat (20 planktonic, 31 hyporheic, 34 sediment,
15 soil). An unbalanced profile is retained as a sampling sensitivity.
Euclidean distances between the four-component profiles represent habitat
affiliation differences. Within each existing bin, `vegan::mantel` correlates
these differences with branch-length phylogenetic distances, following the
Mantel-test core of `iCAMP::ps.bin`. Small pruned trees avoid a full distance
matrix. This is an explicit categorical-habitat extension of the package's
abundance-weighted niche-value approach, not its default continuous-gradient
analysis or an unchanged call to `ps.bin`.

ASVs must occur in at least three samples; bins need at least six retained
ASVs. Every bin remains in the output with its status, retained abundance,
correlation, and one-sided positive-association P value. Tests request 999
ASV-label permutations, with complete enumeration where fewer permutations
exist. The seed is 20161009 plus the bin ID, shared across the two profiles.
BH correction is applied across bins separately for each profile. Signal
summaries use r >= 0.1 and either nominal or BH-adjusted P <= 0.05, not the
package's permissive default P <= 0.2. Reported abundance fractions use the
entire catchment community as denominator; excluded/untestable ASVs do not
disappear through renormalization. Raw r and P values remain available.

Tables are in `results/icamp_2016/catchment_castor/phylogenetic_signal/`.
Realized habitat associations reflect sampling, geography, and detection as
well as ecological preferences. This diagnostic neither measures physiological
traits nor tests all selection-relevant environmental axes. Future chemistry
checks should use comparable measurements and explicitly restricted sample
sets; do not silently assign water chemistry to terrestrial soils.

The October 9 local run tested 775 of 893 bins, with 17,380 eligible ASVs
representing 80.86% of total community abundance. In the habitat-balanced
analysis, 221 bins meet the nominal criterion (30.49% of total abundance),
and 95 meet the BH-adjusted criterion (14.48%). The abundance-weighted mean
Mantel r across testable bins is 0.122. The unbalanced sensitivity gives 208
nominal and 97 BH-supported bins, representing 26.44% and 13.85% of abundance,
with mean r = 0.118. Signal is therefore supported in a subset of the community,
not uniformly across the bins. Nonsignificant bins are not proof of absent
signal, especially with few retained ASVs.

These profiles are estimated from the same observations as the assembly
analysis, not independent niche measurements. Positive habitat signal supports
one association axis but does not identify selection causes or validate the
process percentages. Continue to treat assembly assignments as exploratory.
The run took approximately three minutes locally. Independent checks reproduced
a seeded bin-level Mantel result, matched the weighted-mean niche calculation
in `dniche` (excluding its undefined normalization of a constant indicator),
and confirmed pruning preserves the required phylogenetic distances.

### Git retention

Track `icamp_result.rds` as the canonical detailed result: it contains the
bin-level scores needed to regenerate the diagnostic figures without expensive
null randomization. Track the CSV summaries, ASV/bin taxonomy, sample/settings
provenance, diagnostic tables, habitat-signal outputs, and four iCAMP PDFs.
The source scripts and cluster settings document how these outputs were made.

Explicit `.gitignore` entries exclude distance caches, the catchment pilot,
the recovery checkpoint, and the duplicate package detail/process exports.
These files are retained locally or on the cluster; none were deleted.

### Numerical validation

A read-only review matched the saved community to the current fixed-10K input,
confirmed all 4,950 unique pairs and finite unit-sum process fractions, and
reconstructed pairwise/bin/overall summaries from the saved confidence scores
and weights (maximum pairwise error 3.22e-15). Habitat labels, pair counts, all
49,260 ASV assignments, and minimum united bin size 24 passed checks. Both
saved full-result formats agree, and the checkpoint is readable with 893 bMPD
bin results. These numerical checks do not establish ecological assumptions.

- [Ning et al. 2020, iCAMP framework](https://doi.org/10.1038/s41467-020-18560-z).
- [Author documentation](https://github.com/DaliangNing/iCAMP1).
- [Within-bin null and special cases](https://github.com/cran/iCAMP/blob/master/R/bNRIn.p.r).
- [betaMPD calculation](https://github.com/cran/iCAMP/blob/master/R/bmpd.r).
- [Sequential process classification](https://github.com/cran/iCAMP/blob/master/R/qp.bin.js.r).
- [Within-bin phylogenetic signal](https://github.com/cran/iCAMP/blob/master/R/ps.bin.r).
- [Abundance-weighted niche values](https://github.com/cran/iCAMP/blob/master/man/dniche.Rd).
