# Supporting figures

Start with the four selected figures in [../manuscript](../manuscript).
The retained supporting material provides sampling context, individual panels,
and checks needed to interpret those results. Presence here does not select a
figure for the manuscript supplement.

| Main result | Supporting material | Producer in `analysis/` |
| --- | --- | --- |
| Figure 1: habitat diversity and assembly | [Sampling map, habitat panels, and iCAMP panels](habitat_diversity) | `02_habitat_diversity.R` |
| Figure 2: sediment molecular richness and richness difference | [Individual panels, dispersion views, molecular classes, and balancing checks](sediment_comparison) | `03_molecular_gradients.R` |
| Distance decay on the matched sediment sites | [Primary views](sediment_comparison/2016_matched_sediment_distance_decay_primary.pdf) against watercourse distance, by flow connection | Supporting 10 |
| Figure 3: first detection | [Individual panels and paired-habitat checks](first_detection) | `04_soil_sharing.R` |
| Figure 4: soil-shared taxa in sediment | [Individual panels](soil_affinity) | `04_soil_sharing.R` |
| Stable association group | [Group membership, stability, null diagnostics, and catchment maps](networks) | Supporting 04–07 |
| Catchment transformation inventory | [Readable assignment summary](transformations/2016_transformation_assignment_summary.pdf) and [label-rank distribution](transformations/2016_transformation_label_rank_distribution.pdf) | Supporting 08 and 09 |

[molecular/](molecular) retains full-site composition and detection diagnostics
used by the sediment comparison. [icamp/](icamp/README.md) retains completed-run
views and diagnostics; the selected assembly contrasts are already in Figure 1.
The additional transformation PCA is descriptive; its detection gradient does
not identify a metabolic pathway.

Retired chemistry, enzyme, sediment-side, and transformation-context figures are
recoverable from [the archive](../../archive/README.md). Definitions and
interpretation limits are in [METHODS.md](../../analysis/METHODS.md).
