# Planned manuscript figures

October 10, 2026. Nathan's six-figure sequence uses the reviewed soil,
hyporheic porewater, and sediment cohort. Figures are candidates for manuscript
review; supporting panels remain available in `figures/`.

| Figure | Current PDF | Data and script |
| --- | --- | --- |
| 1. Alpha diversity across habitats | [Hill diversity](../figures/2016_manuscript_alpha_diversity.pdf) | 15 soil, 27 hyporheic, 34 sediment; `figures/02_current_manuscript.R` |
| 2. Beta diversity ordination | [Weighted UniFrac PCoA](../figures/2016_manuscript_weighted_unifrac.pdf) | Same 76 libraries, refitted axes; `figures/02_current_manuscript.R` |
| 3. First-detected flowpath inventory | [ASV richness and reads](../figures/2016_first_detection_flowpaths.pdf) | Reviewed cohort only; `microbes/15_first_detection_flowpaths.R` |
| 4. Sediment microbial and molecular convergence | [Matched dispersion](../figures/2016_sediment_ms_dispersion_paired.pdf) | Identical 33 sites; `fticr/07_matched_sediment_dispersion.R`, reformatted by `figures/01_sediment_manuscript.R` |
| 5. Molecular classes along the network | [Class composition](../figures/2016_fticr_network_group_class_comparison.pdf) | 44 molecular sites, grouped 24/14/6; `fticr/14_plot_fticr_network_classes.R` |
| 6. Localized microbial–molecular group | [Core-2 colocation map](../figures/2016_sediment_stable_group_colocation_map.pdf) | 34 microbial and 44 molecular profiles, 45 distinct sites; `fticr/13_plot_stable_group_colocation.R` |

## Caption and interpretation notes

**1.** Hill numbers at q = 0, 1, and 2 describe richness, exponential Shannon
diversity, and inverse Simpson diversity from the existing 10,000-read draw.
Points are individual libraries; boxes show medians and interquartile ranges,
with whiskers extending to observations within 1.5 interquartile ranges.
Soils are a regional comparison pool, not paired stream-site observations.
No significance brackets have been added. Individual panels are retained as
`2016_soil_porewater_sediment_diversity_q0/q1/q2.pdf`.

**2.** Normalized weighted UniFrac emphasizes abundant phylogenetic lineages.
The saved distance matrix uses the master tree pruned to the original fixed-10K
catalog and midpoint-rooted in memory. The plotted PCoA is newly fitted to the
76-library submatrix, with a Lingoes correction. Axis percentages use positive
corrected eigenvalues. This is a descriptive ordination; habitat tests from
script 19 use Bray-Curtis and cannot be reported as UniFrac tests. Habitat
dispersion differs in Bray-Curtis, limiting a centroid-only interpretation.
Hellinger PCA remains a supporting view.

**3.** ASVs detected in the soil inventory receive soil precedence. Remaining
ASVs are assigned to the first network class with a detection, retaining
sediment–hyporheic ties. Receiving sediment and hyporheic panels remain
separate. Bars show mean within-library percentages, not percentages of
pooled reads. Both richness and read panels are retained until the manuscript
choice is settled. Sampling sensitivity rebuilds the inventory in 1,000 draws
of three libraries per habitat–network cell plus three soils (21 libraries per
draw). First detection does not establish geographic origin or movement.
See [the inventory methods](microbes/README_FIRST_DETECTION.md).

**4.** Distances to group spatial medians compare headwater (orders 1–2), intermediate
(3–4), and mainstem (5) sediment profiles. Microbial Bray-Curtis and molecular
equal-feature Jaccard are refitted on the same 33 sites, using shared
segment-balanced selections. The molecular panel shows the equal-feature
sensitivity, not the complete-profile primary representation; complete-profile
and full-site supporting panels must accompany interpretation. Molecular
dispersion declines in that sensitivity, whereas the microbial test is not
significant. This contrast does not itself test a difference between responses
or establish downstream mixing. Methods and nulls are in
[the matched-site report](fticr/README_MATCHED_SEDIMENT_DISPERSION.md).

**5.** Bars average within-site proportions of detected formulas assigned to
broad molecular classes, weighting each site equally. Headwater/intermediate/
mainstem groups contain 24/14/6 sites. These are formula-based chemical classes,
not identified molecular structures or confirmed metabolites. Peak ambiguity,
detection effects, and class boundaries limit interpretation. Individual-site
panels remain available. See [the molecular methods](fticr/README_RIVER_COMPOSITION.md).

**6.** The map describes the selected stable-link group of 38 ASVs and eight
formula signatures. Missing profiles remain distinct from nondetection.
Nathan's map-based observation is localization in upper Lookout Creek and the
Mack Creek branch draining Lookout Mountain. Branch membership and an
above/below-confluence comparison remain unaudited; describe this as an
exploratory localization hypothesis, not a confirmed branch effect or a
metabolic guild. See [group identities and distribution](fticr/README_BIPARTITE_NETWORKS.md).

## Regeneration

```r
source(here::here("analysis", "microbes", "18_prepare_manuscript_cohort.R"))
source(here::here("analysis", "microbes", "19_soil_porewater_sediment.R"))
source(here::here("analysis", "microbes", "15_first_detection_flowpaths.R"))
source(here::here("analysis", "figures", "02_current_manuscript.R"))
```

The last plotting script reuses completed distances and writes ordination
scores, variance, input checksums, and session information to
`results/manuscript_figures_2016/`. Sediment figure scripts above reuse their
existing inputs; no DADA2, phylogeny, iCAMP, or network-null rerun is required.
