# Sediment ASV–molecular bipartite network design

Status on 2026-10-04: the full bipartite workflow is implemented. The completed
results and verification record are reported below. The original design is
retained afterward to distinguish implemented choices from possible extensions.

## Run the implemented workflow

Run directly in R or Positron, without arguments:

```r
source(here::here("analysis", "fticr", "08_bipartite_feasibility.R"))
source(here::here("analysis", "fticr", "09_bipartite_networks.R"))
source(here::here("analysis", "fticr", "10_bipartite_consensus.R"))
source(here::here("analysis", "fticr", "11_stable_group_distribution.R"))
source(here::here("analysis", "fticr", "check_bipartite_networks.R"))
```

Scripts use `here()` paths. The main inference and consensus scripts have
`test_run <- FALSE`; the main script's test mode runs nine draws per null on
the full-size input matrices in a separate pilot directory. Pilot outputs are
excluded explicitly from Git; the main analysis tables and figures are retained.
The representative full-size pilot took about 6.5 seconds for 36 null networks.

## Implemented methods

The primary analysis retains ASVs and molecular features present and absent in
at least 10 of the 33 paired sediment sites. This excludes both rare features
and nearly ubiquitous features from the incidence network; ubiquitous ASVs can
still vary in abundance, so these results do not cover all possible associations.
The filter retains 2,405 ASVs and 902 molecular features. Identical detection
vectors are collapsed to 2,392 ASV and 850 molecular patterns, with every original
identifier retained in `feature_pattern_map.csv`. This yields 2,033,200 tested
cross-type correlation coefficients rather than treating duplicate patterns as
independent evidence.

Edges are descriptive effect-size selections: phi >= 0.6 for positive links and
phi <= -0.6 for negative links. There are no individual edge P values or claims
of edge-level false-discovery control. Negative associations are exported but
do not define positive modules. Each positive graph is fitted by a transparent
BRIM implementation in `bipartite_helpers.R`, using Barber's bipartite modularity.
It alternates assignments of the two node sets, comparing eight starts: connected
components, one initial group per ASV pattern, and random initial group budgets
of 4, 8, 16, 32, 64, and 128 within each connected component. Disconnected
components remain separate. The maximum-sweep budget is 100 and the objective
tolerance is 1e-10. Reported modularity is a heuristic lower bound, not a certified
global maximum. Observed, null, adjusted, and resampled graphs use the same budget.

There are 999 draws for each of four distinct nulls: whole molecular-profile
permutations; whole profiles permuted within recorded stream order; molecular
incidence randomization with fixed row and column sums; and positive graph
randomization with fixed degrees in both node sets. The ASV matrix remains fixed
in every data-level null. Fixed-margin and fixed-degree nulls use vegan's
nonsequential binary `quasiswap`, with `thin = 100`. Graph simulation omits
zero-degree rows and columns, which have no effect on this objective and make
the simulation unnecessarily slow. Every fixed-margin draw is checked for binary
values and exact margin preservation. Each draw has an exported simulation seed
and a separate optimizer seed; all failure and convergence fields are retained.

For each statistic, the upper-tail P value is (1 + number of null values at least
as large as observed) / 1,000. Six association-count tests form one BH correction
family, and four modularity tests form another. Graph randomization cannot test
edge counts because those counts are fixed. The profile permutation results are
exploratory: stream-order restriction preserves broad order structure but does
not ensure geographic or stream-segment exchangeability.

Module and edge stability are assessed by refitting after omitting each of the
26 stream segments. These refits retain the original feature universe and use
per-refit seeds. Co-membership comparisons use memberships rather than numerical
module labels, and an absent graph node counts as unsupported co-membership.
The complete pair-stability table covers every cross-type pair in an observed
module, rather than only its selected edges.

The consensus screen retains observed positive links with at least 80% support
for both the phi >= 0.6 threshold and co-membership across segment-removal fits,
plus positive partial phi after adjustment for total molecular detection count
and log10 drainage area. Connected components of these links are candidate
stable groups. Connectedness does not guarantee every cross-type pair has 80%
co-membership; `stable_core_summary.csv` therefore reports the mean, minimum,
and qualifying fraction across *all* pairs within each candidate group. The
three-pattern minimum in each node set is a descriptive reporting rule.
Neither individual modules nor consensus groups receive a selection-adjusted
P value. Partial correlations and the adjusted graph are descriptive, with no
conditional permutation inference.

Additional sensitivity graphs use raw phi thresholds of 0.5 and 0.7 and the
adjusted phi threshold of 0.6. They do not receive separate null tests. Molecular
data remain binary throughout; sequencing plate is excluded from ecological
adjustment. An ASV-abundance sensitivity and a fully spatially conditional null
remain possible extensions, rather than completed results.

## Output interpretation and figure captions

Main tables are in `results/fticr_2016/bipartite_networks/`. Start with
`null_tests.csv`, `module_summary.csv`, `stable_core_summary.csv`, and
`candidate_edges.csv`. Original identifiers, taxonomy, and molecular properties
are in the feature-membership and corresponding taxonomy/property tables.
Chemical-property means describe original features, whereas graph calculations
use unique detection patterns.

## Completed results on 2026-10-04

All 3,996 null draws and all 26 segment-removal fits completed and converged.
The null calculations took 687 seconds locally. The positive graph contains
554 links among 242 ASV and 245 molecular patterns; there are also 293 negative
links. Its fitted modularity is 0.7534, with 74 modules, of which 14 contain at
least three patterns of each node type. Most small modules consist of few links;
module counts alone are not evidence for distinct biological guilds.

| Comparison | Statistic | Observed | Null mean | P | BH-adjusted P |
| --- | --- | ---: | ---: | ---: | ---: |
| Whole-profile shuffling | Positive links | 554 | 434.6 | 0.200 | 0.400 |
| Whole profiles within stream order | Positive links | 554 | 508.7 | 0.313 | 0.470 |
| Molecular incidence with fixed margins | Positive links | 554 | 367.3 | 0.001 | 0.003 |
| Graph with fixed degrees | Modularity | 0.7534 | 0.6350 | 0.001 | 0.004 |

Negative link counts also exceed the fixed-margin expectation: 293 observed
versus 181.0 on average (P = 0.001, BH-adjusted P = 0.003). No whole-profile
negative-count or data-level modularity test is significant. In particular,
raw optimized modularity is *lower* than the fixed-margin data-null mean of
0.8810, illustrating why modularity cannot be interpreted without the density
and topology constraints of the specific comparison.

Thus, strong associations exceed incidence-margin constraints, and the selected
positive graph has more modular topology than expected from node degrees. The
whole-profile nulls provide no convincing evidence that the actual cross-domain
site matching has more links or modularity than expected while preserving each
dataset's internal covariance. The null that preserves detection counts disrupts
molecular–molecular covariance, so its rejection does not isolate microbe-specific
chemistry coupling from that covariance or shared environmental structure.

Broad module boundaries are unstable: mean all-cross-pair co-membership in the
14 larger observed modules ranges from 0.133 to 0.581 across segment removals.
All 554 observed positive links retain positive partial phi after adjusting
for molecular detection count and drainage, but this is a descriptive direction
check, not a conditional significance test. Threshold sensitivity is substantial:
phi >= 0.5 gives 4,842 positive links and phi >= 0.7 gives 52, compared with 554
at the primary 0.6 threshold. Conclusions about a specific partition require
more support than the apparent presence of modules.

The stable-link screen retains 261 positive links and produces 93 connected
groups, mostly pairs or small stars. Only one contains at least three patterns
of each node type. This candidate group, core 2, contains 37 ASV patterns
representing **38 ASVs**, eight molecular features, and **99 stable links**.
Their median raw phi is 0.645 and mean partial phi is 0.691. Selected-link
threshold support averages 98.9% of segment-removal refits. Across all 296
cross-type pairs in this group, mean co-membership is 83.1%; the minimum is
53.8%, and 63.2% of pairs reach 80% support. This is a connected set of stable
links, rather than a uniformly stable all-pair block.

The group spans 12 microbial phyla. The largest assignments are Acidobacteriota
(11 ASVs), Pseudomonadota (10), and Gemmatimonadota (5); most ASVs lack a genus
assignment. Six of the eight molecular formulas contain nitrogen and four contain
sulfur. Mean H/C and O/C are 1.598 and 0.348. These descriptions identify candidate
taxonomic and chemical signatures; they do not establish substrate use or
microbial production. The group has no individually calibrated selection-adjusted
P value. Its stable links and the global null results motivate further examination,
while the intact-profile null results limit claims about specific cross-domain
coupling.

The validation script passed analytical block/complete/empty graph examples,
verified the binary phi formula against all Pearson coefficients (maximum
difference 7.2e-16), independently recalculated modularity (difference 1.2e-14),
and confirmed no improving coordinate assignment at the observed solution. All
36 full-size pilot draws reproduced at the same seeds in the final run. All five
inference input checksums and all consensus input checksums were unchanged.
Four PDFs were rendered through macOS CoreGraphics and visually inspected,
along with the 600-dpi association-matrix PNG. No remote action, commit, or push
was performed. These results remain exploratory.

- `2016_sediment_bipartite_association_nulls.pdf`: distributions of positive
  cross-type link counts under the three data-level nulls. The orange line marks
  the observed count; each panel uses its own axis range. Links are selected by
  phi >= 0.6, not individual edge significance.
- `2016_sediment_bipartite_module_null.pdf`: optimized Barber modularity under
  fixed-degree bipartite graph randomization. The orange line marks observed
  modularity. This is a topology comparison conditional on the selected graph.
- `2016_sediment_bipartite_groups.pdf`: observed modules containing at least
  three patterns of each node type. Labels give module number and ASV/molecular
  pattern counts; position gives all-pair co-membership across segment removals,
  point size gives within-module links, and color gives mean adjusted phi among
  those links. Module numbers are arbitrary identifiers.
- `2016_sediment_bipartite_association_matrix.png`: raw phi coefficients for
  positive-graph active nodes, ordered by observed module. Blue denotes positive
  and orange denotes negative correlation. Rows/columns exclude graph isolates.
  This dense matrix is a white-background raster at 600 dpi; the other figures
  are vector PDFs.
- `2016_sediment_bipartite_stable_group.pdf`: the largest candidate stable-link
  group (core 2). Circles are ASV detection patterns and squares are molecular
  features. Only the 99 retained stable positive links are drawn. Labels show
  original representative feature IDs; `(+1)` denotes one additional ASV with
  the same detection pattern, recorded in the feature-pattern mapping. ASV node
  colors show the representative's phylum, with less frequent phyla pooled for
  display. This is a selected exploratory group, without an individual group
  significance test.

## Core 2 identities and catchment distribution

`11_stable_group_distribution.R` describes the existing selected group without
refitting the network or testing a newly chosen spatial hypothesis. Tables are
in `results/fticr_2016/stable_group_distribution/`. It includes all 34 fixed-10K
sediment microbial sites and all 44 sediment molecular sites: 45 distinct sites,
33 with both profiles. Missing profiles remain missing, not nondetections.
Original ASV counts and read percentages are exported; molecular summaries remain
binary. The 37 collapsed ASV patterns are summarized only on the 33 inference
sites, since equivalence on those sites need not extend to microbial-only site 43.

Taxonomy resolves 11 ASVs to a genus: Bryobacter (3), Candidatus Solibacter,
Hirschia, JTB255 marine benthic group, Nitrospira, Nordella, Phaeodactylibacter,
Sphingorhabdus, and Terrimonas (one each). Twenty-seven ASVs lack a genus and
none has a species assignment. The group contains Acidobacteriota (11),
Pseudomonadota (10), Gemmatimonadota (5), Bacteroidota (3), Planctomycetota (2),
and one each from seven other phyla. `core_asv_identities.csv` preserves all
38 IDs, their taxonomic ranks, detection counts, and their shares of group reads.

| Feature | Workbook formula | Broad workbook signature | Candidates |
| --- | --- | --- | --- |
| FTICR_11736 | C10H20N2O3 | Protein-like | 1 |
| FTICR_14319 | C13H25N3O4 | Protein-like | 1 |
| FTICR_14857 | C14H25N3O4 | Lipid-like | 1 |
| FTICR_15869 | C13H25N3O6 | Protein-like | 1 |
| FTICR_18143 | C17H29NO3S2 | Lipid-like | 1 |
| FTICR_20959 | C18H12O9S | Condensed-hydrocarbon-like | 2 |
| FTICR_20963 | C12H24N2O3S5 | Lipid-like | 4 |
| FTICR_21109 | C18H14O9S | Condensed-hydrocarbon-like | 3 |

The original `Candidates` column is retained as `candidate_count`; three peaks
have multiple candidates. These formulas and elemental-ratio class labels do
not identify structures, proteins, lipids, or metabolic substrates. A unique
formula candidate also leaves structural isomers unresolved. The nitrogen-rich,
high-H/C formulas and two low-H/C, oxygen/sulfur-containing formulas describe
different chemical signatures within the group. Structural characterization
would require complementary evidence, such as fragmentation, separation, and
reference standards; see [Leyva et al., structural classification from fragmentation
pathways](https://pmc.ncbi.nlm.nih.gov/articles/PMC11293370/). The five-sulfur
assignment has four candidates and should remain provisional. No exact compound
names or structures were assigned in this analysis.

On the **33 paired sites**, median group read percentages are 0.08% in headwaters
(18 sites), 1.075% in intermediate reaches (10), and 0.06% in the mainstem (5).
Median ASV detections are 3.5, 27.5, and 4 out of 38; median molecular detections
are 3, 8, and 7 out of eight. These are descriptive summaries of a selected group,
not independent spatial significance tests. High representation is patchy:

- Site 190 contains 36 ASVs and all eight features (1.82% of microbial reads).
- CC-4 contains 36 ASVs and six features (1.69%); WS1-5 has 25 and six (1.41%).
- KC-1/KC-3, 53, 66, and 72 contain all eight features and 27–35 ASVs.
- Mainstem 49 and 156 contain 30/35 ASVs and eight/seven features, respectively.
- Sites 161 and 163 have neither component detected. The three paired WS3 sites
  have only 0–2 group ASVs; their molecular detections are 3, 2, and 0.
- Site 50 contains all eight features but only four group ASVs (0.06% of reads).
  Thus, the chemical signature is not a deterministic marker of the microbial set.

At least one member occurs in 28/34 microbial and 34/44 molecular profiles;
all eight molecular features occur at 12 sites, 11 of them paired. There is no
single binary definition of a whole cluster being present. The full profiles
in `core_site_profiles.csv` allow each component to be inspected separately.

All 45 coordinates and orders were rechecked against their recorded centerline
segments; maximum coordinate discrepancy is below 1 m and all orders match.
Stream distance uses the centerline distance to the basin outlet, not valley-local
metadata distance. All seven input checksums remained unchanged. Three vector
PDFs were rendered and visually inspected:

- `2016_sediment_stable_group_asv_map.pdf`: fraction of the original 38 ASVs
  detected at each of 34 sites, with site identifiers. Crosses mark 11 additional
  molecular sites with no retained microbial profile; they do not indicate absence.
- `2016_sediment_stable_group_molecular_map.pdf`: fraction of eight molecular
  features detected at each of 44 sites. A cross marks microbial-only site 43.
- `2016_sediment_stable_group_catchment_maps.pdf`: the same components side by
  side, ASVs in panel A and molecular features in B. Dark colors indicate more
  group members detected, not higher molecular concentration or microbial biomass.

## Spatial interpretation and priority follow-up

On October 4, 2026, Nathan interpreted the catchment maps as placing much of
core 2 in upper Lookout Creek and the Mack Creek branch draining Lookout Mountain,
commonly upstream of their confluence. This is a map-based geographic
interpretation; branch assignments and an above/below-confluence comparison
have not yet been audited or quantified. It supplies a specific hypothesis
for a drainage-associated microbial–molecular signature.

This selected group is a candidate integration finding despite weak
whole-composition correspondence (33-site Mantel rho = 0.112, P = 0.1019;
two-axis Procrustes correlation = 0.130, P = 0.8281). The modest molecular-property
partial-RDA association remains a separate positive result. A subset can have
coherent distributions while aggregate composition shows limited correspondence,
but the group's stability does not establish production, consumption, or coupling
independent of shared environment. Whole-profile nulls remain nonsignificant,
and core 2 has no individual selection-adjusted significance test.

Priority follow-up is to audit branch/confluence assignments, quantify both
components on the paired sites, assess confounding with order and drainage,
and examine detection count, shared segments, sediment properties, and enzyme
activities. Shared geology and organic-matter inputs are proposed explanations;
they are not established mechanisms. Preserve the post hoc status of this
spatial hypothesis in any subsequent comparison. The current synthesis and
detailed follow-up list are in [2016 key findings](../KEY_FINDINGS_2016.md).

## Original design and interpretation boundaries

## Scientific questions

A bipartite association network can identify ASVs and molecular signatures whose
distributions covary across sediment sites. The two node sets are ASVs and
FT-ICR-MS features; every edge connects the two sets. Sites supply the replicate
observations used to estimate edges, rather than forming a third node set.
Retain feature identifiers because formula assignments do not uniquely identify
chemical compounds or metabolites.

Separate two questions: (1) whether cross-site ASV–molecular alignment exceeds
that expected when the two datasets are unrelated, and (2) whether inferred
associations form modules beyond what node degree alone would predict. Module
stability is a third criterion: a statistically unusual graph need not contain
reproducible memberships.

An overall association network includes shared responses to environmental
gradients. A conditional network targets associations remaining after specified
gradients and detection effects are accounted for. These are different ecological
questions; report both explicitly if both are implemented. Neither establishes
microbial production, consumption, or direct interaction.

## Verified inputs and feasibility

Run `analysis/fticr/08_bipartite_feasibility.R` directly in R or with `Rscript`.
It uses the fixed `asv_counts_10k.csv` and
`data/derived/fticr_2016/sediment_44_presence_absence.csv`. Raw inputs and previous
analyses are preserved. Audits are in
`results/fticr_2016/bipartite_feasibility/`.

The intersection contains 33 sites, 26 distinct stream segments, 17,764 detected
ASVs, and 4,496 detected molecular features. Sediment site 43 lacks a molecular
profile. All microbial samples contain 10,000 reads. Molecular detection counts
range from 391 to 2,638 per paired site. Stream orders 1–5 have 2, 16, 5, 5, and 5
sites, respectively. Four segments have repeated sites; resampling should account
for those groups.

The complete detected catalog permits 79,866,944 cross-domain candidate edges.
For a binary–binary incidence analysis, illustrative filters requiring at least
the following numbers of present and absent sites give:

| Minimum present and absent | ASVs | Molecular features | Candidate edges |
| --- | ---: | ---: | ---: |
| 5 | 5,862 | 1,872 | 10,973,664 |
| 7 | 4,118 | 1,415 | 5,826,970 |
| 10 | 2,405 | 902 | 2,169,310 |
| 11 | 1,952 | 778 | 1,518,656 |

These were feasibility checks before selecting the primary 10-site rule. Presence-only
filters would be appropriate for an abundance-based ASV sensitivity because
ubiquitous ASVs can still vary in abundance. In an incidence analysis, identical
site-detection patterns cannot be distinguished: at the 10-site filter, the
2,405 ASVs represent 2,392 patterns and the 902 molecular features represent 850
patterns. Retain identifiers and document these equivalence groups rather than
interpreting their duplicate edges as independent evidence.

## Proposed first analysis

1. Begin with a transparent incidence-based network: molecular presence/absence
   and ASV presence/absence at fixed sequencing depth. Quantify positive and
   negative associations separately; fit modules to positive associations so
   that avoidance does not become evidence of grouping. A phi coefficient is
   Pearson correlation between the two binary detection vectors. An ASV
   Hellinger-abundance sensitivity can assess dependence on incidence, while
   molecular data remain binary.
2. Choose a prevalence rule and an edge effect-size rule before inspecting
   network results. Test an overall statistic before pursuing millions of
   individual edge tests. For example, compare the number of sufficiently strong
   positive associations to the same statistic in reconstructed null networks.
   Effect-size-selected edges are descriptive unless they also pass an explicit
   inferential procedure. Do not label a thresholded graph a significant-edge
   network without multiplicity control across all eligible pairs.
3. Use bipartite-specific modularity, such as Barber's modularity, with the same
   optimization procedure, starting-point budget, and seeds for observed and null
   graphs. Finding modules or obtaining positive modularity is insufficient on
   its own: optimized random graphs can also contain apparent modules.
4. Evaluate stability with leave-one-site-out influence checks and leave-one-
   segment-out or segment-level resampling. Reestimate associations and refit
   modules each time. Track sign/effect consistency and co-membership frequency,
   which does not depend on arbitrary module labels. Flag nodes with insufficient
   variation after resampling rather than counting them as failed associations.
5. Summarize stable groups by ASV taxonomy and molecular elemental/property
   signatures, retaining the original feature mapping. Relate their distributions
   to measured sediment properties only after establishing stability. Compare
   threshold and detection sensitivities before interpreting modules biologically.

## Null hypotheses and what each preserves

| Null | Preserved structure | Question answered |
| --- | --- | --- |
| Permute whole molecular site profiles relative to ASV profiles | All within-molecular associations, all within-microbial associations, feature prevalences, and each dataset's distribution of site richness | Is cross-domain site matching unusually structured? |
| Randomize binary incidence with fixed site richness and feature prevalence | Row and column sums of each randomized incidence matrix | Does structure exceed detection/richness and prevalence constraints? |
| Rewire inferred positive edges within the bipartite graph, preserving both sets of node degrees | Node identities, bipartite restriction, edge count, and degree sequence | Are modules stronger than expected from graph degree alone? |

These nulls are complementary and are not interchangeable. Whole-profile
permutation relocates molecular richness among sites; it does not preserve each
site's original detection count. Fixed-margin incidence randomization disrupts
within-molecular or within-microbial associations. Graph rewiring tests topology
conditional on inferred edges and cannot validate the original associations or
test an excess number of edges.

Rebuild networks from randomized data with the same analysis settings for
data-level nulls. If edge selection uses estimated significance or resampling,
that selection must also be reproduced under the null. For degree-preserving
rewiring, document swap acceptance, degree preservation, and adequate mixing;
restricted degree sequences can sharply limit possible graphs. Start with a
representative timed run locally before choosing a randomization count or
considering cluster execution.

For an overall alignment comparison, restricted whole-profile permutations within
stream order or justified spatial groups are candidates, but exchangeability must
be assessed. With two first-order sites and repeated segments, there is no
automatic valid spatial restriction. Restricting by order does not remove
within-order geographic dependence. A conditional comparison should address
molecular detection count and a small specified network-position adjustment,
such as log drainage area, and use a justified conditional randomization or
reduced-model residual-permutation procedure. Simply residualizing both matrices
and freely shuffling rows is not automatically valid. Sequencing plate is
excluded from ecological adjustment.

## Interpretation and planned outputs

The intended result is a set of reproducible cross-site associations and, if
supported, groups of ASVs and molecular signatures with consistent distributions.
A sparse or unstable result is informative and should be retained. Formula-level
chemistry and 16S associations support ecological hypotheses, not metabolic
assignments. Environmental filtering and detection can generate the same patterns
as resource relationships.

An implemented analysis should export the site and feature inclusion audits,
edge definitions and weights, full null-run diagnostics with per-run seeds,
observed/null statistics, node memberships and co-membership stability, and
taxonomy/property mappings as inspectable tables. Useful figures are the null
distribution with the observed statistic and a module-ordered cross-domain
association heatmap. A dense network drawing can be secondary. Put all figures
in `figures/`; use vector PDF for compact panels and a white 600-dpi raster for
a dense heatmap if needed.

## Methodological grounding

- [Barber (2007), Modularity and community detection in bipartite networks](https://doi.org/10.1103/PhysRevE.76.066102):
  bipartite modularity and module optimization.
- [Connor, Barberán, and Clauset (2017), Using null models to infer microbial co-occurrence networks](https://doi.org/10.1371/journal.pone.0176751):
  distinguish data-level association nulls from graph-level structural nulls.
- [Carr et al. (2019), Use and abuse of correlation analyses in microbial ecology](https://doi.org/10.1038/s41396-019-0459-z):
  compositionality, environmental filtering, and limits of interaction inference.
- [Winkler et al. (2020), Permutation inference for canonical correlation analysis](https://arxiv.org/abs/2002.10046):
  residualization can affect exchangeability in cross-domain multivariate tests.
- [vegan binary null-model documentation](https://vegandevs.github.io/vegan/reference/commsim.html):
  margin constraints and nonsequential quasiswap simulation.

These references support methodological choices. The implementation and results
above identify the completed work; the original design includes additional
possibilities that were not all implemented.
