# Reviewed first-detection inventory

Run `15_first_detection_flowpaths.R` after scripts 18 and 19. The input is the
reviewed fixed-10K cohort: 15 soil, 34 sediment, and 27 hyporheic porewater
libraries. No planktonic or flagged water libraries enter the inventory.

The descriptive inventory adapts the upstream-first framework of
Ruiz-Gonzalez et al. (2015), doi:10.1111/ele.12499. A single read counts as a
detection. Soil detections have precedence. For remaining ASVs, search stream
orders 1–2, then 3–4, then 5; assign the first detected pool as sediment,
hyporheic, or a sediment–hyporheic tie. The two stream compartments are parallel
pools; there is no assumed flow order between them. Network classes describe
the survey inventory, not an audited upstream path for every receiving site.
Soils are regional rather than individually assigned contributing sources.

The four categories exhaust each receiving library. Bars average individual
library richness or read percentages within network groups. Sediment counts
are 19/10/5 and hyporheic counts are 14/8/5 in headwater/intermediate/mainstem
groups. All 15 soils enter the full inventory. All 33,974 retained ASVs receive
a full-inventory assignment, including soil-only ASVs.

The 1,000-draw equal-sample sensitivity selects three libraries per each of six
stream habitat–network cells and three soils: 21 libraries per draw. The
inventory is rebuilt and evaluated on selected receiving samples each time.
The seed is 20261005 plus the draw ID. Equal counts do not guarantee equal
detection completeness. Draw quantiles are sampling sensitivities, not
confidence intervals or independent replication.

Outputs in `results/diversity_2016/first_detection_flowpaths/` retain sample
identities, ASV assignments, per-sample profiles, full and balanced summaries,
all 21,000 selections, diagnostics, settings, and unchanged-input checksums.
The three current vector PDFs are:

- [Combined richness and reads](../../figures/manuscript/2016_first_detection_flowpaths.pdf)
- [Richness](../../figures/supplement/2016_first_detection_richness_percent.pdf)
- [Reads](../../figures/supplement/2016_first_detection_read_percent.pdf)

These displays describe detection overlap. They do not estimate source
contributions, recruitment, viability, transport rates, or direction. The
previous plankton-inclusive inventories and plots are recoverable from Git;
any archived result subdirectories retain their historical meanings.
