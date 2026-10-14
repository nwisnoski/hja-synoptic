# iCAMP figures

Start with the completed reviewed soil–porewater–sediment pool. Its five overview
figures are in `reviewed_pool/`:

| View | Sediment | Porewater |
| --- | --- | --- |
| Segment-weighted process composition | [Composition](reviewed_pool/2016_reviewed_pool_icamp_sediment_network_composition.pdf) | [Composition](reviewed_pool/2016_reviewed_pool_icamp_porewater_network_composition.pdf) |
| Site summaries across partner segments | [Maps](reviewed_pool/2016_reviewed_pool_icamp_sediment_catchment_maps.pdf) | [Maps](reviewed_pool/2016_reviewed_pool_icamp_porewater_catchment_maps.pdf) |

The [shared-site habitat comparison](reviewed_pool/2016_reviewed_pool_icamp_matched_habitats.pdf)
uses the same 18 sites and unordered site pairs in both habitats.

`sediment_sensitivity/` holds the sediment-only
[composition](sediment_sensitivity/2016_sediment_icamp_network_composition.pdf)
and [maps](sediment_sensitivity/2016_sediment_icamp_catchment_maps.pdf).
Changing the regional pool changes the bins and null reference; these are
separate inference runs.

Within each run, `panels/` retains individual process panels. `diagnostics/`
contains balanced-sampling comparisons, within/between-position heatmaps, and
all-pair distributions and composition. Balanced ranges describe sampling
sensitivity. All-pair views give every sample pair equal weight, including
same-segment pairs; overview composition weights distinct segment pairs equally
and excludes same-segment pairs. These views answer different descriptive
questions and do not provide independent pairwise tests.

Producers: [supporting 13](../../../analysis/supporting/13_reviewed_pool_icamp.R)
for the reviewed pool and [supporting 12](../../../analysis/supporting/12_sediment_icamp_figures.R)
for the sediment-only sensitivity. Neither needs fresh inference with completed
results present. See [Methods](../../../analysis/METHODS.md) for interpretation
and [saved runs](../../../results/icamp_2016/README.md) for provenance.
