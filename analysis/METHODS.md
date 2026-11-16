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

The [2016 preparation methods](https://doi.org/10.5194/essd-11-1567-2019)
describe three sequential water extractions of 500 mg sediment, pooled to 5 mL,
followed by PPL extraction and fixed acquisition settings. They do not report
matching extracts to a common carbon concentration.

FT-ICR ecology uses binary detections. The preparation filter retains assigned
formulas with measured mass 200–900, C13 <= 0, and positive intensity in at
least two of 60 profiles. Raw intensity is retained for QC. Detected-feature
counts depend on detection and are not calibrated molecular richness.

## Putative molecular transformations

`supporting/08_molecular_transformations.R` applies the per-sample method and
July 2020 database from [Danczak et al. (2020)](https://doi.org/10.1038/s41467-020-19989-y)
to the existing sediment feature set. The [reference record](../reference/danczak_2020/README.md)
pins the original source and database. Absolute differences between measured m/z
values are matched within an inclusive 0.000010 mass-unit tolerance; both peaks
must be detected in the sample. Each unordered pair counts once and receives the
last matching database label in source order, reproducing the published code.
All candidate labels are retained in the pair-assignment table. The optimized
matching is checked against the unmodified upstream script on a complete sample.

Transformation profiles count assigned peak pairs by label. Their descriptive
PCA uses square roots of within-sample pair-count proportions (Hellinger), with
no intensity weighting. Connected-peak degrees include zeros; match counts and
possible-pair counts are retained as detection diagnostics. These connections
do not establish molecular structures, reaction direction, realized activity,
or biogeochemical rates. The database, inherited feature filter, measurement
error, and label-order ambiguity condition the profiles. Environmental hypothesis
tests and causal process assignments are not part of this descriptive analysis.

`supporting/09_transformation_summary.R` ranks reference labels by their mean
within-sample fraction and reports site coverage, distinct pairs, alternative
labels, and mass error. Its curated subset of simple motifs additionally requires
exact net elemental differences and the same measured-mass tolerance. Both
endpoints having one reported formula candidate distinguishes less ambiguous
assignments from pairs with multiple candidates; neither category identifies a
molecular structure or establishes a reaction. Motif checks are calculated
within this script from the completed pair assignments and molecular inventory.

The bulk summaries describe chemical relationships across sampled sediments.
Hydrogen, CH2, oxygen, water, acetyl-residue, and carboxyl-group equivalents lead
the curated subset; hexose and amino acid residue equivalents are less frequent
but broadly distributed. Sparse nutrient-bearing differences do not diagnose
nutrient limitation or the absence of nutrient cycling. These inventories do not
establish functional convergence or assign metabolism to the mapped ASV group.
Enzyme associations, detection-conditioned gradients, and transformation
network-context tests were retired; their methods and results are recoverable
from [the cleanup archive](../archive/README.md).

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

## Figure 2: first detection and matched-habitat iCAMP contrasts

A single read establishes detection. Soil detections have precedence; remaining
ASVs are assigned at the first surveyed stream-order group containing them:
orders 1–2, 3–4, then 5. Sediment and porewater pools are parallel; ties remain
shared. Panel A has sediment and porewater facets; bars summarize
percentages of ASV reads per library. Percentages of detected ASVs remain
supplementary.
The 1,000 balanced draws select three libraries in each of six stream
habitat–network cells and three soils, rebuilding the inventory each time.
Draw quantiles describe sampling sensitivity, not confidence intervals.
This inventory does not establish contributing drainage, transport direction,
source proportions, viability, or recruitment. It adapts the inventory framework
of [Ruiz-Gonzalez et al. (2015)](https://doi.org/10.1111/ele.12499).

Panels B–C show inferred dispersal limitation and homogeneous selection,
respectively, from the completed reviewed soil–porewater–sediment iCAMP pool.
Both habitats use the same 18 sites: 11 headwater, three intermediate, and four
mainstem. Small translucent points show all unordered within-position site
pairs in each habitat (55, three, and six pairs), including same-segment pairs.
Each site pair receives the same horizontal jitter in both habitats and both
process panels; there is no vertical jitter. Large opaque points and connecting
lines show arithmetic means. Both panels use the same
0–60% scale; no significance tests or confidence intervals are attached.
Pairwise rows share sites, and the intermediate group is sparse. These
assignments depend on the null model and regional-pool definition; they do not
measure dispersal or selection rates or establish soil-to-stream transport.
The first-detection inventory and iCAMP panels use different summaries and
sample coverage; their juxtaposition does not establish a causal link.
The full process composition, balancing sensitivities, and sediment-only run
remain supplementary; the iCAMP methods are given below.

The supplementary paired soil-sharing plot connects percentages for sediment
and porewater at each of the 18 paired sites. Sharing is the percentage of reads belonging to ASVs
detected in any of the 15 soils. Its two-sided test exchanges habitat labels
within sites and tests the mean paired difference. This plot tests the habitat
contrast underlying the soil component of panel A; it is not independent
corroboration of that inventory. Script 02 retains 1,000 three-soil draws with
receiving communities fixed, plus drainage regressions on stream-segment means
with a four-test BH correction; fitted lines require adjusted P <= 0.05.

The supplementary taxonomic-resolution plot shows individual same-site
Hellinger distances as thin gray lines and their mean as a single black line,
from ASV through genus, family, order, class,
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

## Figure 3: molecular richness and richness differences along flow paths

All panels use the 33 matched sediment sites. Each panel has a fixed printed size
and the canvas is measured around the panels. Panels A and B relate sediment
organic content, detected feature count, and drainage area by least-squares
regression of log10 variables; slopes are scaling exponents. Three simple
regressions (feature count on drainage area, organic content on drainage area,
feature count on organic content) form one BH family; the first is saved as a
supplementary panel. A joint model of feature count on organic content and
drainage area asks whether network position explains richness beyond organic
content. The same models on all 44 molecular sites are a check, and ASV richness
at equal sequencing depth is the microbial contrast. Sites are treated as
independent, although some share stream segments.

Extracts used a fixed sediment mass, not a common carbon concentration, and
feature count tracks total signal intensity. Feature count is therefore detected
features per unit sediment, not calibrated molecular richness. More organic
sediment may yield more features through detection alone.

Panels C and D use presence and absence for both data types. Binary Jaccard
dissimilarity between every pair of sites (528 pairs) is the sum of replacement,
where each site holds features or ASVs the other lacks, and richness difference,
where one site holds more than the other (Podani and Schmera, 2011;
[Carvalho et al., 2012](https://doi.org/10.1111/j.1466-8238.2011.00694.x)).
Total beta diversity and the share from each component summarize each data type
over all pairs ([Legendre, 2014](https://doi.org/10.1111/geb.12207)). The
abundance-weighted form is not used: it measures differences in total abundance,
which are zero for samples rarefied to equal depth. Panels C and D plot each
pair's richness-difference share, its richness difference divided by its Jaccard
dissimilarity, against watercourse distance.

Pairs are flow-connected when their watercourse distance equals the difference
of their outlet distances, so water passes from one site to the other (120 of
528 pairs). The stream network is built from the surveyed centerline, joining
consecutive points within segments and segments at shared confluence points.
Graph distances to the outlet reproduce the centerline's own values; one
unsampled first-order segment (135) does not reach the outlet. Replacement,
richness difference, and the richness-difference share are each correlated with
watercourse distance and with the log drainage-area ratio, separately for
flow-connected and unconnected pairs, and the two slopes are compared. Null
distributions permute whole sites in the component matrix (9,999 permutations)
while separations and flow connections stay in place; tests are two-sided. For
each separation, slopes within pair sets and slope differences form separate BH
families. In panels A and B, fitted lines are drawn only for adjusted
P <= 0.05. In panels C and D, every flow-connected and unconnected slope is drawn
for comparison, solid for adjusted P <= 0.05 and dashed otherwise; lines span the
distances of the pairs they were fitted to.

Replacement cannot exceed one minus richness difference, so the two components
are not independent: where richness difference rises, replacement must fall.
Richness difference does not imply that the poorer site holds a subset of the
richer one. Flow-connected pairs are not independent either: each mainstem site
is connected to 13–29 others. The same detection limit applies as for feature
count, since lower-carbon extracts could lose rare features without any change in
the underlying molecular pool. The pattern is consistent with downstream loss of
site-specific terrestrial features
([Creed et al., 2015](https://doi.org/10.1139/cjfas-2014-0400)) but does not
establish it.

Dispersion among network groups is supplementary. Median-based tests use 9,999
permutations on Lingoes-corrected distances. Primary views are microbial
Hellinger and complete molecular Jaccard; the Jaccard replacement views form a
second family. Network groups are orders 1–2, 3–4, and 5. The 1,000
segment-balanced draws use identical site selections across all views; their
frequencies are descriptive. Replacement is not Euclidean, so its Lingoes
constant is large; the saved raw matrices allow correction-free checks. The equal-feature
view, which thins every site to 391 detections, is retained only as a
diagnostic: thinning binary detections removes more of a richer profile, so a
site's distance rises with its feature count (`richness_dependence.csv`).
Network position is confounded with geography in this survey.

Molecular class fractions at the matched sites and at all 44 sites are
supplementary. Bars use workbook classifications and fractions of binary
detections. Formula classes are broad signatures, not identified structures or
confirmed metabolites.

## Supporting analysis: distance decay on the matched sediment sites

`supporting/10_distance_decay.R` relates microbial Hellinger and molecular
Jaccard dissimilarity between every pair of the 33 matched sites (528 pairs) to
watercourse and straight-line distance. It reads the complete-data matrices and
the pair separations and flow connections saved by script 04. Their replacement
and richness-difference components are tested in script 04.

Pearson Mantel tests permute whole sites, share 9,999 permutations, and are
two-sided. The microbial–molecular contrast standardizes each dissimilarity
across pairs and relates the pairwise difference to each distance. Flow-connected
and unconnected pairs are tested as for Figure 3. Figures show watercourse
distance only, with a line for a pair set only where adjusted P <= 0.05.

## Figure 4: soil-shared taxa in sediment across the stream network

`05_soil_affinity.R` calls an ASV soil-core when it occurs in at least five of
the 15 reviewed soils. The class uses the soil samples only, so the sediment
data it is tested on play no part in defining it. Its share of a sample's
10,000 reads is measured for all 34 sediment and 27 porewater samples, counting
every soil-core ASV. Occurrence in soil does not establish origin or movement;
soil-core ASVs may be habitat generalists.

Each soil-core ASV is also classed by whether its mean abundance is higher in
sediment than in soil, and the share of soil-core sediment reads from ASVs more
abundant in sediment is recorded. This split does not depend on which lineages
panel C shows; selecting lineages by soil abundance makes most of them appear
diluted in sediment by construction.

Sediment and porewater shares are compared at the 18 paired sites by exchanging
habitats within sites (9,999 sign flips). Spatial structure in the sediment
share uses Moran's I with inverse watercourse-distance weights. Each sample is
assigned to the north branch (above segment 134), the Lookout branch (above
segment 39), or the channels below their fork, from the same stream graph as
Figure 3. Drainage-area trends use Spearman correlations; the six site tests
form one BH family.

Panel C shows the 25 soil-core lineages with the highest mean percentage of
soil reads, so it shows the lineages that dominate soils rather than taxa rare
in soil; together they hold 69% of soil-core reads in soil. Each lineage sums
its ASVs' mean percentages in soil, sediment, and porewater, drawn on a log
axis; a zero is not drawn. Rows are ordered by the ratio of soil to sediment
abundance, from soil-dominant lineages at the top to lineages more abundant in
sediment than in soil below the dashed line. A lineage is the family, or
"Unclassified" at the lowest assigned rank; SILVA placeholder families are
named by their order.

## Supporting analysis: sediment bipartite network and the stable group

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
eight molecular features. Most of its ASVs never occur in soil. It is
concentrated on channels draining Lookout Mountain (the Lookout branch and the
KC tributary) and nearly absent from the north branch. The group was selected
from a network whose association counts do not exceed the profile-shuffle nulls,
so its microbial–molecular concordance is not an independent result, and no
selection-adjusted test establishes direct interaction, production,
consumption, or a metabolic guild.

## Supplementary sediment iCAMP maps and network comparisons

`analysis/supporting/icamp/03_sediment_sensitivity_figures.R` summarizes the
completed sediment-only iCAMP sensitivity
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

`analysis/supporting/icamp/01_reviewed_pool.R` contains the completed pooled
analysis of the reviewed 76-sample,
33,974-ASV cohort: soil, hyporheic porewater, and sediment. It preserves the
original fixed-10K counts and uses 1,000 randomizations (seed 2016), midpoint
rooting after pruning, bMPD/Confidence, within-bin phylogenetic randomization,
across-all taxonomic randomization, minimum bin size 24, and ds = 0.2. It
requires iCAMP 1.9.1 and a distinct, validated distance cache. This analysis
is the primary exploratory iCAMP view; the sediment-only run remains a sensitivity
to the regional-pool definition. Saved assignments are reused when plotting.

The regional abundance weights are mean relative abundances across all 76
samples, following the iCAMP default; habitats receive weights proportional
to their sample counts. Adding soils and porewater changes the ASV pool, bins,
and null reference. This is a sensitivity to the regional-pool definition,
not source attribution or a test of directional soil-to-stream dispersal.
All soils remain in the regional pool, including unresolved identity/GPS flags.

Figures use only sediment-sediment or porewater-porewater assignments from
the completed pooled run. Within each habitat, maps, segment-pair means, five-segment
balanced draws, and all-within-position pair distributions follow the
sediment-only figure script.
Map color scales vary by process and habitat; compare numerical legends rather
than color alone. Sediment and porewater have different spatial coverage. A
separate sensitivity retains their 18 shared sites (11 headwater, three
intermediate, four mainstem) and compares the same unordered within-position
site pairs in both habitats, with pair-level porewater-minus-sediment contrasts
retained. Lines connect descriptive group means; no significance tests or
confidence intervals are attached. The shared-site intermediate group is sparse.
