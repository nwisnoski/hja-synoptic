# Sediment-only manuscript exploration

**Current framing:** the October 10
[soil-porewater-sediment outline](SOIL_POREWATER_SEDIMENT_MANUSCRIPT.md)
supersedes this alternative for manuscript planning. The sediment analyses
and figures below remain available as within-sediment results.

October 10, 2026. This is an alternative manuscript structure; existing
catchment-wide analyses and figures are preserved. Sediment iCAMP job 432228
was submitted with approval on October 10; queue/run state is not yet checked.

## Central question and framing

Determine how sediment microbial communities and molecular signatures vary
across a dendritic river network, and whether their correspondence remains
after accounting for measured sediment conditions. Network position provides
the spatial setting; organic content, texture, hydraulic conductivity, and
enzyme activities provide distinct measurements of local sediment conditions.
The study can identify spatial associations and null-model assembly signatures,
but cannot measure dispersal direction or assign molecular production and
consumption to taxa.

Working title: **Sediment microbial communities and molecular signatures
across a dendritic river network**.

Use three linked questions:

1. How do microbial composition and molecular signatures vary across network
   positions? A shared downstream homogenization pattern is a hypothesis;
   existing results support detection-equalized molecular convergence without
   a supported microbial dispersion decline.
2. Which measured sediment characteristics covary with these patterns, and
   does microbial–molecular correspondence remain after adjustment? Organic
   content has a modest microbial association. The fixed molecular-property
   association weakens after organic-content/enzyme adjustment. Neither result
   establishes a causal mechanism.
3. What assembly signatures occur within the sediment metacommunity, and are
   any microbial–molecular associations localized? Sediment-only iCAMP is the
   pending assembly analysis. The stable bipartite group is an exploratory
   feature-level result whose geographic interpretation remains post hoc.

The introduction can move from dendritic connectivity and local sediment
conditions to the unresolved correspondence between community composition and
organic-matter composition. The study contribution is the joint comparison
within one habitat across the network. Avoid introducing soil recruitment or
water-habitat filtering as central mechanisms in this version.

## Sampling and analysis structure

Retain all 34 legitimate fixed-10K sediment microbial samples in composition,
diversity, dispersion, and iCAMP. Molecular descriptions use all 44 sediment
sites. Only analyses requiring both datasets use their 33-site intersection.
Site 43 remains in microbial analyses. Missing conductivity or texture affects
only the analysis requiring that measurement: microbial counts are 33 and 14,
respectively; molecular counts are 41 and 17. Organic-content/enzyme comparisons
retain all 34 microbial samples and all 44 molecular sites.

Do not use a single environmental complete-case filter for the manuscript.
Show the union and coverage audit, distinguish unavailable profiles from
nondetection, and report each model's own sample size. Matched 33-site
comparisons remain useful for direct cross-dataset figures, accompanied by
the full 34/44-site analyses so the paper retains all legitimate information.

## Results sequence and figures

The first sediment-focused review set is now generated. See the
[figure report and captions](figures/README_SEDIMENT_FIGURES.md) for the four
candidate main figures, retained individual panels, full-sample supporting
views, and sediment-only UniFrac ordinations. All original figures remain
intact. The table below records the resulting figure sequence.

| Figure | Purpose | Existing material or work required |
| --- | --- | --- |
| 1 | Network geography and sampling | `2016_sediment_ms_sampling_map.pdf` shows the full 45-site union: 34 microbial, 44 molecular, 33 paired profiles. Site 43 is labeled. Sediment coverage remains available as a supplement. |
| 2 | Within-sediment composition and chemical signatures | `2016_sediment_ms_composition.pdf` combines all-34 microbial Hellinger PCA and all-44 molecular class summaries. Individual panels, molecular PCoA, diversity, and within-sediment UniFrac are retained. |
| 3 | Network variation in dispersion | `2016_sediment_ms_dispersion_paired.pdf` uses the same 33 sites in both panels. `2016_sediment_ms_dispersion_all_samples.pdf` retains all 34/44 sites and complete molecular profiles as supporting views. Different significance is not a test of response differences. |
| 4 | Sediment conditions and cross-dataset correspondence | `2016_sediment_ms_characteristics.pdf` combines all-eligible-sample characteristic models with paired molecular-property adjustment. Individual panels and all-34 microbial PCA colored by organic content are retained. |
| 5 | Sediment assembly | Pending sediment-only iCAMP: process fractions and appropriate bin/null-score diagnostics. Validate within-sediment environmental phylogenetic signal before strong process interpretation. |
| 6, optional | Localized microbial–molecular signatures | Existing `2016_sediment_stable_group_colocation_map.pdf` and stable-link network. Use only if the feature-level result contributes beyond Figures 2–4; otherwise place in supporting material. |

This is a candidate figure set. No publication composites have been declared
final. Water optics, soil sharing/routing, cross-habitat diversity, and catchment
iCAMP remain available as broader context or a separate manuscript; they do
not need to carry the main argument of this version. Water measurements must
retain their actual compartment labels if included as sediment context.

## New local results

The [characteristics report](fticr/README_SEDIMENT_CHARACTERISTICS.md) documents
the implemented analyses and limitations. Organic content explains 3.99%
adjusted microbial variation beyond drainage (BH P = 0.0102), retaining all
34 samples and remaining positive across segment-deletion checks. Primary
enzyme-profile, conductivity, and grain-size tests are nonsignificant.
Complete-profile molecular associations are small and nonsignificant after
the declared corrections.

At all 33 paired sites, the fixed molecular-property association declines
from 4.04% to 2.44% adjusted variation after controlling organic content and
enzymes (raw P = 0.0447, BH P = 0.0622; within-order raw P = 0.1000).
This weakens a claim of microbial–molecular correspondence independent of
local sediment conditions, while leaving selected feature-level associations
as explicit exploratory hypotheses. The study cannot infer that sediment
conditions mediate an underlying causal relationship.

## Remaining work

1. Obtain approval to check all-34 sediment iCAMP job 432228 and its logs, as
   recorded in [the deployment note](cluster/README_SEDIMENT_DEPLOYMENT.md). Do not replace
   it with within-sediment pairs from the catchment run: the null pools differ.
2. Review sediment iCAMP completion, binning, null distributions, sample/input
   alignment, and within-sediment environmental phylogenetic signal. Existing
   catchment habitat-affiliation signal tests cannot validate this null pool.
3. Review the generated sediment-focused figure set and captions for manuscript
   selection. Individual panels and all-sample counterparts are retained.
4. If the localized group remains central, audit branch assignments and assess
   covariation with sediment conditions and detection. Do not frame its
   selected spatial pattern as a confirmed branch effect or metabolic guild.

Results and deployment preparation are local and uncommitted. An approved
October 10 remote preflight found matching inputs, no queued jobs, and no
sediment output paths. The separately approved two-script transfer completed
October 10; the separately approved submission was accepted as job 432228.
Its queue/run state has not been checked. No commit or push was performed.
