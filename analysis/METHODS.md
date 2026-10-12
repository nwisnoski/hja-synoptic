# Methods needed to read the figures

Settings, seeds, transformations, sample selections, and diagnostics are kept
in the scripts and result tables. This guide defines what the figures measure.

## Samples and preparation

Figures distinguish soil, sediment, and hyporheic porewater (hereafter porewater).
Source metadata retain the original habitat code `hyporheic` for porewater samples.

`01_prepare_cohort.R` validates each suspect library against its site, habitat,
lane ID, and technical plate metadata before excluding it. The 76 retained
libraries preserve the existing 10,000-read draw; only ASVs absent from this
cohort are removed. The full 120-library audit distinguishes low read depth,
suspect-library QC, and plankton outside the manuscript scope. Plate is never
an ecological predictor. Contamination remains a hypothesis, not a diagnosis.

Sediment analyses retain all 34 eligible samples. Molecular analyses use
44 sediment profiles; direct comparisons refit on the identical 33-site
intersection. Site 43 lacks a molecular profile. Missing profiles remain
distinct from nondetections. Thirteen soils have resolved locations; soil 04
lacks GPS and `hja2016_200` has unresolved identity. All 15 enter regional sharing.

FT-ICR ecology uses binary detections. The preparation filter retains assigned
formulas with measured mass 200–900, C13 <= 0, and positive intensity in at
least two of 60 profiles. Raw intensity is retained for QC. Detected-feature
counts depend on detection and are not calibrated molecular richness.

## Supplementary sampling coverage

Four map panels show the 13 located soils, 27 porewater libraries, 34 sediment
libraries, and 44 sediment molecular profiles. Black rings identify the 33
sediment sites with both data types. The two unmapped soils remain in the
regional inventory. Maps use equal coordinate scales.

## Figure 1: habitat diversity

The main alpha-diversity panel shows Hill q = 1, the exponential of Shannon
entropy; q = 0 and 2 remain supplementary. Weighted UniFrac PCoA subsets the completed normalized
distance matrix to the 76 samples and refits with a Lingoes correction. The
underlying master tree was pruned and midpoint-rooted in memory for UniFrac.
Supplementary Hellinger PCA is also refitted on the reviewed cohort. Panel A
shows Hill diversity and panel B shows individual libraries in the weighted
UniFrac ordination, without centroid ellipses.

Companion Hellinger habitat PERMANOVA and median-based, bias-adjusted
PERMDISP use separate BH families for the overall test and three pairwise contrasts. Unrestricted
permutations do not control geography or repeated sites; dispersion differences
limit a pure centroid interpretation, especially with unequal group sizes
([Anderson and Walsh, 2013](https://doi.org/10.1890/12-2010.1)). These tests use
Hellinger rather than the plotted weighted UniFrac distances. Hellinger is
Euclidean distance between square-root relative-abundance profiles, with range
0 to sqrt(2) ([Legendre and Gallagher, 2001](https://doi.org/10.1007/s004420100716)). The 18-site
sediment–porewater test conditions on site and exchanges habitat labels only within site.

Script 02 also retains porewater environmental Mantel checks, which exchange
whole site labels and use separate eight-test BH families for unrestricted and
within-order permutations.

## Figure 2: first detection and paired habitat comparisons

A single read establishes detection. Soil detections have precedence; remaining
ASVs are assigned at the first surveyed stream-order group containing them:
orders 1–2, 3–4, then 5. Sediment and porewater pools are parallel; ties remain
shared. Panels A–B show sediment and porewater, respectively; bars summarize
percentages of ASV reads per library. Percentages of detected ASVs remain
supplementary.
The 1,000 balanced draws select three libraries in each of six stream
habitat–network cells and three soils, rebuilding the inventory each time.
Draw quantiles describe sampling sensitivity, not confidence intervals.
This inventory does not establish contributing drainage, transport direction,
source proportions, viability, or recruitment. It adapts the inventory framework
of [Ruiz-Gonzalez et al. (2015)](https://doi.org/10.1111/ele.12499).

Panel C connects soil-sharing percentages for sediment and porewater at each
of the 18 paired sites. Sharing is the percentage of reads belonging to ASVs
detected in any of the 15 soils. Its two-sided test exchanges habitat labels
within sites and tests the mean paired difference. This panel tests the habitat
contrast underlying the soil component of panels A–B; it is not independent
corroboration of that inventory. Script 02 retains 1,000 three-soil draws with
receiving communities fixed, plus drainage regressions on stream-segment means
with a four-test BH correction; fitted lines require adjusted P <= 0.05.

Panel D shows individual same-site Hellinger distances as thin gray lines and
their mean as a single black line, from ASV through genus, family, order, class,
and phylum across the 18 pairs.
The mean curve for comparisons to the other 17 sediments nearly overlaps the
same-site curve at every rank and is omitted for clarity;
both comparisons remain in the saved tables and supplementary sensitivity plot.
Aggregation sums counts before the Hellinger transformation; its expected
reduction in distance does not itself
demonstrate ecological convergence. Taxa are grouped by complete lineage
prefixes to avoid conflating homonyms. ASVs with missing or generic unresolved
labels at any rank in the prefix remain ASV-specific, preserving every read
without treating unclassified taxa as shared. Sample-level classification
coverage is saved. A supplementary sensitivity retains the same ASVs with
complete genus lineages at every rank, renormalizing each sample.

Tests permute intact porewater profiles among sediment sites, globally and
within exact stream order, using the same 9,999 permutations at every rank and
in both read scopes. Each scope and null has its own six-rank BH family.
The ASV-level Hellinger, UniFrac, and read-sharing comparison in script 02 has
four metrics corrected separately per null. No detectable enhancement of
same-site resemblance does not establish equivalence, temporal transience,
or sediment retention of a local signal.

## Figure 3: sediment dispersion and molecular classes

Panels A–B show molecular equal-feature Jaccard and microbial Hellinger
dispersion, respectively, refitted on the same 33 sites. The molecular
sensitivity averages 1,000 draws of 391 detections per site. Complete molecular and binary microbial Jaccard are
additional views. Network groups are orders 1–2, 3–4, and 5. Median-based
dispersion tests use 9,999 permutations; the two primary views form one BH
family. The 1,000 segment-balanced draws use identical site selections across
all views. Their frequencies are descriptive, not independent tests. A difference
in significance does not itself test a difference between microbial and molecular
responses. Lower molecular dispersion is supported in the equal-feature view;
complete molecular Jaccard does not show a significant network-group dispersion
difference. Network position is confounded with geography in this survey.

Panel C shows molecular class fractions at these identical 33 sites, ordered
by drainage within network groups. Bars use workbook classifications and
fractions of binary detections. Broad class fractions provide chemical context;
they do not directly measure convergence in individual formula detections.
The full 44-site class figure and equally weighted group means remain
supplementary. Formula classes are broad signatures, not identified structures
or confirmed metabolites.

## Figure 4: exploratory association-group distribution

The sediment bipartite network retains features present and absent in at least
10 of 33 sites, collapses identical detection patterns, and selects phi >= 0.6
positive links. Negative links use phi <= -0.6 and do not define modules.
BRIM fits [Barber bipartite modularity](https://doi.org/10.1103/PhysRevE.76.066102)
with eight starts; the optimum is heuristic. Individual links have no calibrated
P values. Each of four nulls has 999 draws:

- Whole-profile shuffles preserve molecular covariance and break site identity.
- Within-order profile shuffles also preserve broad stream-order structure.
- Molecular fixed-margin shuffles preserve site totals and feature prevalence,
  while disrupting molecular covariance; ASV profiles stay fixed.
- Fixed-degree graph shuffles preserve node degrees and test modularity.

Association-count and modularity tests use separate BH families. Segment-removal
refits assess stability. The 80% link/co-membership screen, with positive partial
phi after detection-count/drainage adjustment, identifies candidate connected
groups rather than uniformly stable blocks. The mapped group has 38 ASVs and
eight molecular features. Its apparent upper Lookout/Mack concentration is an
exploratory map-based hypothesis; no selection-adjusted branch test establishes
localization, direct interaction, production, consumption, or a metabolic guild.

## Supplementary sediment iCAMP maps and network comparisons

Supporting script 12 summarizes the completed sediment-only iCAMP run
(34 fixed-10K samples, 1,000 null randomizations). The 561 sample-pair assignments
are relative process importance under the [iCAMP framework](https://doi.org/10.1038/s41467-020-18560-z).
The residual category includes drift and other processes; these assignments
are not measured local rates or evidence of transport direction.

Each mapped site averages its comparisons with sites in other stream segments,
first within each partner segment and then equally across the 26 partner
segments. Same-segment partners are excluded. Map color scales differ among
processes to show spatial variation; their legends give percentages. The
network-position map identifies orders 1-2, 3-4, and 5. These are relational
site summaries, not process estimates for the land area draining to each site.

Network comparisons exclude same-segment pairs and first average sample pairs
within each distinct, unordered segment pair. Segment pairs then receive equal
weight within each of six within-position or between-position comparisons.
Composition bars show the three within-position means on a common percentage
scale. In 1,000 balanced draws, each position contributes five randomly selected
segments and one randomly selected site per segment; all processes and pair
types use the same selection. Points and ranges show draw medians and middle
95% ranges; hollow diamonds show full-data segment-pair means. Panel y-axis
ranges differ. All five mainstem sites enter every draw, so their draw range is
zero. These ranges describe sampling sensitivity, not confidence intervals.
Selections and draw-level contrasts are retained in
`results/icamp_2016/sediment_castor/network_figures/`.

Segment weighting and balanced selection address unequal representation and
repeated sampling within segments. They do not model flow connectivity,
shared upstream drainage, or spatial dependence among branches. Pairwise rows
share samples, and network position is confounded with geography; no independent
pairwise tests or causal network-position claims are made.

The all-within-position figures additionally show every headwater-headwater,
intermediate-intermediate, and mainstem-mainstem sample pair with equal weight
(171, 45, and 10 pairs). Same-segment pairs are included and marked with crosses;
white diamonds show arithmetic means. Boxes show medians and interquartile
ranges with standard 1.5-IQR whiskers; every pair is plotted, including values
beyond the whiskers. The stacked figure shows the same arithmetic means.
Panel y-axis ranges differ. Pair counts are not independent sample sizes, and
these figures carry no significance tests or confidence intervals.

## Reviewed-cohort regional-pool iCAMP sensitivity

Supporting script 13 defines a fresh pooled analysis of the reviewed 76-sample,
33,974-ASV cohort: soil, hyporheic porewater, and sediment. It preserves the
original fixed-10K counts and uses 1,000 randomizations (seed 2016), midpoint
rooting after pruning, bMPD/Confidence, within-bin phylogenetic randomization,
across-all taxonomic randomization, minimum bin size 24, and ds = 0.2. It
requires iCAMP 1.9.1 and a distinct, validated distance cache. This analysis
has been prepared but not executed; it does not replace the sediment-only run.

The regional abundance weights are mean relative abundances across all 76
samples, following the iCAMP default; habitats receive weights proportional
to their sample counts. Adding soils and porewater changes the ASV pool, bins,
and null reference. This is a sensitivity to the regional-pool definition,
not source attribution or a test of directional soil-to-stream dispersal.
All soils remain in the regional pool, including unresolved identity/GPS flags.

Figures use only sediment-sediment or porewater-porewater assignments from
the fresh pooled run. Within each habitat, maps, segment-pair means, five-segment
balanced draws, and all-within-position pair distributions follow script 12.
Map color scales vary by process and habitat; compare numerical legends rather
than color alone. Sediment and porewater have different spatial coverage. A
separate sensitivity retains their 18 shared sites (11 headwater, three
intermediate, four mainstem) and compares the same unordered within-position
site pairs in both habitats, with pair-level porewater-minus-sediment contrasts
retained. Lines connect descriptive group means; no significance tests or
confidence intervals are attached. The shared-site intermediate group is sparse.
