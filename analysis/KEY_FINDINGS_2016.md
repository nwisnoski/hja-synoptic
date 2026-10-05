# 2016 synoptic key findings and follow-up

Updated October 4, 2026. This synthesis distinguishes completed analyses from
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
