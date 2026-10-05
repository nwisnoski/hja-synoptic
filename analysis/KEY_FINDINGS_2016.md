# 2016 synoptic key findings and follow-up

Updated October 5, 2026. This synthesis distinguishes completed analyses from
map-based interpretation and proposed follow-up. The sediment integration uses
33 paired sites; descriptive microbial and molecular coverage is 34 and 44
sites, respectively. Molecular ecological analyses use incidence, and microbial
analyses use the fixed 10K ASV table.

## Working integration finding

Microbial–molecular correspondence may be concentrated in geographically
structured subsets of features. Whole-composition comparisons show weak
correspondence, but the bipartite analysis identifies a candidate stable-link
group of 38 ASVs and eight molecular features. Nathan's interpretation of the
catchment maps on October 4 places much of this group in upper Lookout Creek
and the Mack Creek branch draining Lookout Mountain, commonly upstream of their
confluence. The geographic observation has not yet been quantified with an
audited branch assignment. This group provides a focused integration result
and a hypothesis about a shared drainage-associated microbial–molecular
signature, even where aggregate compositional correspondence is weak.

## Key findings so far

1. **A candidate microbial–molecular group provides a spatial integration lead.**
   Core 2 comprises 37 ASV occurrence patterns representing 38 ASVs, eight
   molecular features, and 99 stable positive links. Across segment-removal
   refits, retained-link threshold support averages 98.9%, and all-cross-pair
   co-membership averages 83.1%. The minimum all-pair support is 53.8%, so the
   group is a connected stable-link set with variable internal support.
   Acidobacteriota, Pseudomonadota, and Gemmatimonadota account for 11, 10, and
   five ASVs. The formulas include six nitrogen-containing and four
   sulfur-containing assignments. The Lookout Mountain drainage interpretation
   gives this selected group a concrete spatial hypothesis for follow-up.
   [Network analysis](fticr/README_BIPARTITE_NETWORKS.md#core-2-identities-and-catchment-distribution)
   and [group summary](../results/fticr_2016/bipartite_networks/stable_core_summary.csv).

2. **Weak aggregate correspondence coexists with a modest molecular-property
   association.** At 33 paired sites, complete molecular Jaccard and microbial
   Hellinger distances have Mantel rho = 0.112 (P = 0.1019); two-axis Procrustes
   correlation is 0.130 (P = 0.8281). Equal-feature mean molecular distances
   show still weaker correspondence (rho = 0.0053, P = 0.4670). However, two
   molecular-property PCs explain 4.04% adjusted unique microbial variation
   after conditioning on drainage and detection count (P = 0.0103,
   BH-adjusted P = 0.0309). These outcomes motivate feature- and property-level
   integration without implying strong whole-assemblage concordance.
   [Composition report](fticr/README_RIVER_COMPOSITION.md#molecularmicrobial-composition)
   and [partial-RDA results](../results/fticr_2016/integrated_chemistry/microbial_chemistry_partial_rda.csv).

3. **Network structure depends on what the null preserves.** Strong positive
   and negative link counts exceed molecular fixed-margin expectations
   (both BH-adjusted P = 0.003), and modularity exceeds the fixed-degree graph
   null (BH-adjusted P = 0.004). Whole-profile shuffling, including within
   stream order, does not support excess positive links or modularity. The
   fixed-margin null disrupts molecular covariance; the group has no individual
   selection-adjusted significance test. The completed network results support
   an exploratory association signature, while specific microbial–chemical
   coupling remains unresolved.
   [Null comparisons](../results/fticr_2016/bipartite_networks/null_tests.csv).

4. **Molecular and microbial dispersion have different descriptive network
   patterns.** On the same 33 sites, equal-feature molecular headwater/mainstem
   dispersion ratio is 1.126 (omnibus BH-adjusted P = 0.0008); microbial
   Bray–Curtis ratio is 0.927 (P = 0.5174). Complete molecular and binary
   microbial Jaccard sensitivities are nonsignificant. Geography and detection
   remain relevant, and contrasting significance does not itself test a
   difference between the two responses.
   [Matched-site dispersion](fticr/README_MATCHED_SEDIMENT_DISPERSION.md)
   and [test results](../results/fticr_2016/matched_sediment_dispersion/three_group_tests.csv).

5. **Detection effects and chemical identification constrain interpretation.**
   Molecular detection count correlates with total primary signal
   (rho = 0.969), and complete-profile PCoA1 tracks detection count
   (rho = -0.984). Equal-feature thinning is a sensitivity analysis rather than
   calibrated analytical-depth rarefaction. A monotonic downstream source shift
   is unresolved. Formula assignments and protein-like/lipid-like labels do
   not identify molecular structures, metabolites, or production/consumption;
   three of core 2's eight features have multiple formula candidates.
   These boundaries apply to the candidate integration finding.
   [Detection and composition report](fticr/README_RIVER_COMPOSITION.md)
   and [molecular identities](../results/fticr_2016/stable_group_distribution/core_molecular_identities.csv).

6. **Water-compartment integration remains exploratory.** Matched surface and
   hyporheic optics are available for all 20 planktonic and 31 hyporheic
   communities, respectively; contextual sediment FT-ICR matches are 19 and
   30. None of four partial-RDA comparisons survives BH correction. The sediment
   molecular inventory is contextual site chemistry for these aquatic
   comparisons, so compartment-specific molecular associations remain untested.
   [Aquatic integration report](fticr/README_AQUATIC_MICROBE_CHEMISTRY.md).

7. **Soil-associated sharing distinguishes sediment from water across the
   network.** In the grouped detection inventory, soil-detected ASVs comprise
   mean sediment read percentages of 72.3/67.7/69.9% across headwater,
   intermediate, and mainstem samples. Water communities contain a larger
   fraction of ASVs assigned to the combined hyporheic/planktonic detection
   pool. Soil-detected percentages do not decline monotonically in every
   habitat. The four inventory labels are Soil, Sediment, Aquatic, and
   Sediment + Aquatic; they describe earliest observed detections, not proven
   sources or movement. This is a growing descriptive result and a candidate
   manuscript figure, adapted from
   [Ruiz-Gonzalez et al. (2015)](https://doi.org/10.1111/ele.12499).
   Retain both richness and read rows for now; the paper will likely use one
   row, with the choice deferred.
   [Methods and captions](microbes/README_SOIL_STREAM_LOCALIZATION.md#earliest-detection-along-network-stages),
   [richness row](../figures/2016_first_detection_richness_percent.pdf),
   [read row](../figures/2016_first_detection_read_percent.pdf), and
   [plotted summaries](../results/diversity_2016/first_detection_flowpaths/full_inventory_network_summary.csv).

## Priority follow-up: Lookout Mountain drainage and the confluence

**Status: noted, not implemented.** The spatial interpretation comes from
Nathan's reading of the maps on October 4, 2026. Shared geology, sediment
properties, and organic-matter inputs are possible explanations, not measured
mechanisms established by the group.

1. Audit each site's branch and position relative to the upper Lookout–Mack
   Creek confluence using the mapped stream topology and Nathan's geographic
   interpretation. Record a reproducible site-to-branch table, source, and
   unresolved assignments; do not infer membership from stream order alone.
2. Quantify the group separately in the Lookout Mountain drainage, the other
   branch, and downstream reaches where the design permits. Use the 33 paired
   sites for joint comparisons; retain the additional microbial/molecular sites
   as explicitly unpaired descriptive coverage. Report ASV detections, group
   read percentages, and molecular detections without inventing a binary
   threshold for the whole group being present.
3. Assess whether branch can be distinguished from stream order and drainage
   in this sample design. Account for molecular detection count and shared
   segments, and inspect branch-specific sampling balance and influential
   sites before choosing a small model or an appropriate spatial null.
4. Compare the group with a small set of measured sediment properties,
   organic content, enzyme activities, and optical signatures. Examine the
   individual members as well as the aggregate group, including mismatches
   such as site 50, which has all eight features but only four group ASVs.
5. Treat same-dataset branch comparisons as exploratory because the group and
   spatial hypothesis were selected after inspecting these data. Any formal
   test must address selection and spatial dependence; independent validation
   would strengthen the proposed drainage-associated signature.

The aim is to determine whether the candidate group reflects a shared
catchment-associated environmental signature and whether cross-type associations
persist after that spatial context is considered.

## Soil ASVs and contributing drainage areas: October 5, 2026

The new soil analysis uses the fixed-10K communities and public-terrain drainage
areas for individual aquatic sites. Thirteen of 15 soils have resolved identities
and GPS coordinates; the other two remain in general sharing. Forty-six of 54
aquatic site catchments pass all outlet-snapping QC scenarios. Unstable soil
memberships and failed catchments are flagged and excluded from spatial
comparisons, while all microbial profiles remain in general sharing.

Of 13,969 soil-detected ASVs, 6,593 occur aquatically and 4,033 occur in mainstem
samples. A descriptive model accounting for aquatic-sample identity, soil-profile
identity, and geographic proximity estimates modest positive localization in
sediment: +1.53 percentage points in the fraction of soil ASVs detected, or
+1.68 points for the operational soil-enriched subset. Water-column effects are
small, and individual-ASV responses include positive, absent, and reversed
associations. No calibrated spatial P values or soil-source percentages are
claimed.

Two disjoint tributary drainage areas above sites 47 and 66 contain two and
three confidently assigned soil profiles. A soil-only abundance/prevalence
screen retains 225/291 branch-concentrated ASVs. Aggregate sediment signatures
are weakly differentiated between the sampled branches and occur downstream
in fifth-order mainstem samples (median read percentages 4.98%/5.62%).
Individual examples show stronger differentiation: Nocardioides ASV_000708
is detected in 3/5 sediment samples in its sampled soil branch, 0/4 in the
other, and 2/3 downstream; it also occurs at seven other sediment sites.
This supports spatially variable sharing rather than exclusive branch sources.
The small water-column branch coverage and shared sediment segments constrain
inference. These branch signatures are separate from the previously selected
bipartite core 2; its branch hypothesis remains untested.

See [soil drainage methods, results, and captions](microbes/README_SOIL_STREAM_LOCALIZATION.md).
Source identities, individual-ASV contrasts, drainage uncertainties, and
mainstem profiles are exported to
`results/diversity_2016/soil_stream_localization/`. All 11 vector figures were
rendered and visually inspected. Scientific framing follows
[Crump et al. (2012)](https://doi.org/10.1038/ismej.2012.9) and
[Read et al. (2015)](https://doi.org/10.1038/ismej.2014.166).

The distance-only follow-up in `analysis/microbes/13_soil_stream_distance.R`
uses straight-line soil-to-aquatic-site distance as an explicitly approved
nearness proxy, independent of drainage assignments. It retains all 13 mapped
soils and 85 aquatic samples. Descriptive soil/site-adjusted distance slopes
are small and mixed in sign, with no consistent aggregate decline in sharing.
Individual soil and ASV correlations vary; 513 of 877 eligible soil-enriched
sediment ASVs have negative nearest-soil distance–abundance correlations
(median rho -0.0485). These are descriptive associations, not calibrated
locality tests. Distance inputs and results are in
`data/derived/soil_stream_distance_2016/` and
`results/diversity_2016/soil_stream_distance/`; the additional vector figure
passed rendered visual inspection.


### Soil figure refinement and earliest detection: October 5, 2026

Following figure review, ASV-example and distance-correlation plots are
removed from active outputs; initial branch plots are archived pending a
revised comparison. The three ASV-sharing maps remain active, and Python
drainage preprocessing, uncertainty audits, and all analysis tables remain
intact. Archived PDFs are in `figures/archive/soil_first_pass/`; plot generation
is disabled by default in the corresponding scripts.

The catchment soil-ASV figure now has separate descriptive GAM curves in
log10 drainage area for each habitat (REML, basis dimension four). Sediment
is approximately linear; water-compartment curves allow modest curvature.
No spatially calibrated tests or confidence bands are asserted.

`analysis/microbes/15_first_detection_flowpaths.R` adapts
[Ruiz-Gonzalez et al. (2015)](https://doi.org/10.1111/ele.12499)
using soils first, then the earliest aquatic network class detecting each
remaining ASV. The current source inventory combines hyporheic and
planktonic water detections as Aquatic, with four categories: Soil, Sediment,
Aquatic, and Sediment + Aquatic. The shared category retains sediment-water
ties at the earliest network stage. Receiving habitat panels remain separate.
Mean within-sample richness and read percentages are displayed
across headwater (orders 1-2), intermediate (3-4), and mainstem (5) groups.
Spatial pairing is not required. This is an observed detection inventory,
not source attribution or a verified connected sequence of flowpaths.

Soil-detected richness percentages are 60.0/53.5/57.4% in sediment,
32.2/27.1/34.2% in hyporheic samples, and 19.2/21.8/13.9% in planktonic
samples across the three stages. There is no universal monotonic decline.
A 1,000-draw sensitivity selecting three profiles per habitat-stage cell
and three soils lowers absolute soil-detected percentages substantially,
while retaining these broad shapes. Equal sample numbers do not establish
equal detection completeness.

Tables, ASV assignments, sampling selections, and validation records are in
`results/diversity_2016/first_detection_flowpaths/`. The combined figure is
`figures/2016_first_detection_flowpaths.pdf`, with both individual rows also
saved. All revised/new PDFs passed rendered visual inspection; input hashes
remained unchanged. No remote action, commit, or push was performed.

Water-pool simplification on October 5 preserves all soil percentages and
earliest network stages. The former shared category splits into 2,744
water-only ties now assigned to Aquatic and 2,134 sediment-water ties retained
as Sediment + Aquatic. Original separate-water tables are archived under
`results/diversity_2016/first_detection_flowpaths/archive_separate_water_pools/`.
The sampling sensitivity still balances original habitats (three profiles
per cell), so the merged water inventory has six profiles per stage versus
three for sediment; it does not balance the merged source pools.
