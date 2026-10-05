# Soil ASVs and contributing drainage areas

Updated locally on October 5, 2026. The active figures are the ASV-sharing
maps, the catchment plot with habitat-specific regression curves, and the
first-detection flowpath figure. ASV examples and distance-correlation figures
were dropped from active outputs; the initial branch figures are archived
pending a revised comparison. Python drainage preprocessing and all tables
are retained.

This first analysis finds extensive
soil–aquatic ASV sharing, a modest sediment localization signal, and variable
individual-ASV responses. Soil signatures from two tributary drainage areas
also occur outside their sampled branch and downstream in the mainstem.
These descriptive results support a mixture of localization and broad sharing;
they do not establish soil-to-stream transport, growth, or source percentages.

## Questions and methods

1. Which ASVs detected in sampled soils occur in planktonic, hyporheic, and
   sediment communities, including the fifth-order mainstem?
2. Is a soil profile's ASV signature more prevalent or abundant where that soil
   lies within the aquatic site's contributing drainage area?
3. Do signatures concentrated in one of two separately sampled tributary soil
   pools remain concentrated in that branch, or occur across the network?

These questions follow terrestrial inoculation and subsequent aquatic sorting
along hydrological continua ([Crump et al. 2012](https://doi.org/10.1038/ismej.2012.9))
and changes in riverine bacterioplankton along network gradients
([Read et al. 2015](https://doi.org/10.1038/ismej.2014.166)). Shared 16S ASVs
provide spatial source context; resemblance alone does not establish origin.

### Community inputs and operational definitions

All ecological results use the existing fixed-10K table: 49,260 ASVs in 100
samples, including 15 soils, 20 planktonic, 31 hyporheic, and 34 sediment
samples. The existing bacterial/archaeal screening and rarefaction are retained.
No DADA2 or phylogeny work is rerun. Habitats are analyzed separately.

- **Soil-detected:** at least one read in at least one of the 15 soils.
- **Soil-enriched sensitivity:** at least five total soil reads and mean
  relative abundance across soils at least five times the mean across all
  85 aquatic samples. These thresholds are operational choices. The subset
  uses aquatic data for screening and is interpreted descriptively.
- **Mapped source:** one of 13 soils with a resolved identity and GPS points.
  Soil 04 lacks GPS; `hja2016_200` has an unresolved soil identity. Both remain
  in overall sharing but are excluded from spatial source assignments.
- **Contributing:** the soil point lies within the DEM-derived drainage area
  upstream of the aquatic sampling outlet. This is not a nearest-stream rule.
- **Branch-concentrated:** detected in at least two confidently assigned soils
  of one branch, at least five total reads in that source pool, and mean
  relative abundance at least three times the mean in the other branch's
  sampled soils. This screening uses soils alone. It does not establish
  exclusivity to a branch or comparison with every other sampled soil.

Primary detection requires one read. A three-read sensitivity is applied to
both soil source detections and aquatic detections in the full localization
analysis. Relative abundance is reported as percentage of the 10,000 reads;
it is not absolute cell abundance or a percentage attributed to a source.

### Drainage membership and uncertainty

`10_soil_drainage_routing.py` uses a frozen public USGS 3DEP elevation extract
on a 10-m grid. See [reference-layer provenance](../../data/spatial/soil_routing/README.md)
for extent, source URLs, horizontal-datum assumptions, and download parameters.
This terrain reference differs from the original 2008 LiDAR geometry documented
by [Ward et al. (2019)](https://doi.org/10.5194/essd-11-1567-2019).

Pysheds 0.5 fills pits and depressions, resolves flats, calculates D8 flow
direction and accumulated area, and delineates catchments. The general
DEM-conditioning and drainage-delineation approach follows
[Jenson and Domingue (1988)](https://pubs.usgs.gov/publication/70142175).
Outlets snap to the highest-accumulation pixel within 30 m, with 15- and 50-m
sensitivities. Each mapped soil is also shifted by -20, 0, and +20 m on both
axes, giving 27 membership scenarios per soil–site pair. These offsets assess
local instability; they are not a calibrated GPS-error distribution.

A catchment passes QC if it does not touch the terrain-extract edge and its
area is between 0.5 and 2 times the independently recorded site drainage area.
Forty-six of 54 distinct aquatic sites pass at all three snapping radii.
At the primary radius, retained area ratios range from 0.641 to 1.266, with
median 1.010. Eight sites fail at least one radius and retain flagged audit
rows: 43, 60, 64, 155, 165, 170, 190, and 260. For example, site 155 snaps to a
larger neighboring channel. Their microbial samples remain in general sharing.

Of 702 mapped soil–site pairs, 23 have membership that changes across the
coordinate/snapping scenarios. The spatial comparisons retain 581 pairs that
are stable and pass every outlet-QC scenario. Because habitats share sites,
these become 202 planktonic, 315 hyporheic, and 341 sediment sample–soil pairs.
There are 16, 25, and 27 aquatic samples with at least one usable soil pairing,
respectively. Missing and uncertain relationships are excluded explicitly.

The public gaged-watershed polygons independently place soils 01, 03, and 18
within Watershed 1. The full aquatic-site drainage areas correctly distinguish
sites above and below these near-outlet soils; watershed membership does not
mean a soil contributes to every upstream sample in that watershed.

### Comparisons and limits of inference

`11_soil_stream_localization.R` reports, for each soil and aquatic sample,
the fraction of the soil's ASVs detected aquatically and the aquatic read
percentage belonging to those ASVs. It compares contributing and other sites,
including contrasts restricted to the same stream order. A descriptive linear
model includes aquatic-sample identity, soil-profile identity, log1p Euclidean
distance, and contributing membership. Sample identity absorbs site-level
differences such as drainage area and order. No sequencing-plate predictor is
used. The model coefficients are effect estimates; no independent-pair P values
or confidence intervals are presented.

Individual-ASV contrasts compare aquatic samples with a contributing mapped
soil detection against samples without one. At least three samples per category
are required within a habitat. If any detected candidate soil has uncertain
membership at a site, that ASV–aquatic comparison is excluded. All nondetections
and ties are retained. Positive and negative counts are descriptive, not a
calibrated prevalence of statistically localized taxa.

Pairs share soil profiles, aquatic samples, ASVs, and stream segments. Geographic
proximity and shared habitat can produce resemblance without transport.
Unsampled soils and imperfect detection limit source inventories, and even
exact ASV matches do not establish strain identity or viable organisms. The
current analysis therefore estimates spatial associations without formal
source-attribution or spatially calibrated hypothesis tests.

## Results

### General sharing and mainstem detections

Of 13,969 soil-detected ASVs, 6,593 (47.2%) occur in at least one aquatic sample
and 4,033 (28.9%) occur in a fifth-order mainstem sample. Soil-detected ASVs
comprise a median 17.40%, 33.07%, and 70.26% of reads in planktonic, hyporheic,
and sediment samples. Much of this sharing includes widespread taxa.

The soil-enriched subset contains 3,479 ASVs, of which 1,634 occur aquatically.
Its median aquatic read percentages are 0.70%, 1.31%, and 5.36%, respectively.
The difference between these operational pools shows why all soil-detected
reads cannot be interpreted as a soil contribution percentage.

### Drainage localization is modest and habitat dependent

After the descriptive sample/source/proximity adjustment, contributing soils
have +0.32, +0.18, and +1.53 percentage points higher fractions of their ASVs
detected in planktonic, hyporheic, and sediment samples. Corresponding aquatic
read-percentage effects are +0.10, +0.02, and +1.39 points.

For the soil-enriched sensitivity, detection-fraction effects are +0.065,
+0.072, and +1.68 points, and read-percentage effects are +0.018, -0.058, and
+0.434 points. The clearest aggregate descriptive association is in sediment.

Within-order soil-enriched contrasts likewise favor sediment locality:
mean detection differences are +3.08 points at the one-read threshold and
+2.65 points at three reads. Hyporheic differences are +0.76/+0.33 points;
planktonic differences are -0.86/-0.58 points. Available strata and soils differ
among habitats, so these are complementary descriptions rather than replicated
tests of habitat differences.

Individual responses vary. Among the 1,210 eligible soil-enriched ASVs detected
in retained sediment comparisons, 510 have higher detection in contributing
catchments and 680 have higher detection without a detected contributing
source; 20 tie. Eligibility excludes ASVs without both spatial categories, so
this count is not representative of the entire soil ASV pool. Aggregate sharing
and individual localization do not imply the same response.

### Two tributary soil signatures also occur beyond their branch

`12_soil_branch_signatures.R` compares two disjoint drainage areas above sites
47 and 66 (1,553.61 and 1,362.42 ha in the reference DEM). Site identifiers
label verified drainage areas; mountain or named-creek assignments are not
asserted. Soils 12 and 19 confidently belong to branch 47; soils 07, 08, and 09
belong to branch 66. Soil 10 is excluded from this source contrast because its
membership near a divide is unstable. Branch choice is exploratory and motivated
by contrasting mapped soil coverage.

The soil-only screen retains 225 branch-47-concentrated and 291
branch-66-concentrated ASVs. Retained sediment coverage is five samples in
branch 47, four in branch 66 (only two stream segments), three downstream of
both, and 14 elsewhere. Planktonic coverage is only two/one/two/11, and
hyporheic coverage three/two/two/17; strong water-column branch inference is
not supported by this small coverage.

Sediment median read percentages for branch-47-concentrated ASVs are 3.58%
in branch 47 and 3.16% in branch 66. Branch-66-concentrated ASVs are 1.48% in
branch 66 and 1.21% in branch 47. Both signatures occur downstream, with medians
of 4.98% and 5.62% at the three mainstem sediment samples. The signatures are
therefore weakly differentiated between their sampled branches and broadly
detectable beyond them. Higher downstream percentages may also reflect
sorting, sediment habitat, and cumulative drainage differences.

There are localized individual examples. *Nocardioides* ASV_000708 occurs in
both branch-47 soils and none of the three branch-66 soils. It is detected in
3/5 sediment samples in branch 47, 0/4 in branch 66, 2/3 downstream mainstem
samples, and 7/14 other sediment samples. Its mean relative abundances are
0.012%, 0%, 0.0167%, and 0.010%, respectively. It provides an example of branch
differentiation with broader downstream/elsewhere detection, not exclusive
localization. Acidobacteriota ASV_000752 shows the opposite branch association:
2/4 branch-66 sediment samples, 0/5 branch-47 samples, and 2/3 downstream.
These examples were selected after inspecting contrasts and have no
selection-adjusted significance claims.

## Outputs, captions, and execution

Constructed drainage inputs and all membership sensitivities are in
`data/derived/soil_drainage_2016/`. Results are in
`results/diversity_2016/soil_stream_localization/`:

- `soil_asv_identities.csv`, `sharing_overview.csv`, and
  `habitat_sharing_summary.csv`: full ASV identities, definitions, and sharing.
- `soil_aquatic_pair_sharing.csv`, `individual_soil_contrasts.csv`,
  `within_order_soil_contrasts.csv`, and `adjusted_descriptive_effects.csv`:
  source-specific spatial results and exploratory adjustment.
- `asv_localization_contrasts.csv` and `asv_localization_summary.csv`:
  per-ASV contributing/noncontributing results at both detection thresholds.
- `branch_source_soils.csv`, `branch_site_context_audit.csv`, and
  `branch_sample_counts.csv`: the two-branch spatial and sampling design.
- `branch_soil_asv_signatures.csv`, `branch_signature_site_profiles.csv`,
  `branch_signature_asv_contexts.csv`, and `branch_signature_asv_contrasts.csv`:
  soil-defined signatures and all aquatic contexts, including nondetections.
- Input checksums, DEM-conditioning diagnostics, Python timing/version, and
  R session information document provenance and execution.

The initial 11 figures were rendered and visually inspected. Eight initial
example, branch, and distance figures now reside in
`figures/archive/soil_first_pass/`; their scripts disable these plots by default.
Sharing maps and individual-soil contrasts remain active. Original captions
for the first-pass figures are retained below:

1. `2016_soil_asv_sharing_map_<habitat>.pdf` (three files): Aquatic read percentage
   belonging to ASVs detected in any sampled soil. Orange triangles mark mapped
   soils. All aquatic samples remain included; fill scales differ by habitat.
   Read percentages do not estimate soil-source contributions.
2. `2016_individual_soil_localization.pdf`: Each soil's contributing-minus-other
   detection-fraction contrast, separately by habitat. Only stable catchment
   assignments enter. Sources lacking either category are absent; no uncertainty
   intervals or independent-source inference are asserted.
3. `2016_soil_asv_examples_<habitat>.pdf` (three files): Four soil-enriched ASVs
   selected by absolute contributing-minus-other abundance contrast. Relative
   abundance is plotted against drainage area, with source context and stream
   order encoded. Panels have separate y scales; selected taxa are descriptive.
4. `2016_soil_branch_drainage_map.pdf`: Disjoint drainage areas above sites 47
   and 66, confidently assigned source soils, other/uncertain mapped soils,
   aquatic sites, and the supplied stream centerlines. Soil 10 is visually
   near branch 47 but excluded from that source pool by sensitivity checks.
5. `2016_soil_branch_signatures_<habitat>.pdf` (three files): Individual aquatic
   sample read percentages belonging to the two soil-defined branch-concentrated
   ASV sets, across branch, downstream, and other contexts. Downstream samples
   are fifth-order mainstem. Groups have unequal samples and shared segments;
   source pools have two versus three soil samples and 225 versus 291 ASVs.

Run the R scripts directly in Positron or from any directory within the project.
Routing is a separate preprocessing step; a pinned Python requirements file
is supplied because pysheds 0.5 is incompatible with NumPy 2.4 and later.

```sh
python3 analysis/microbes/10_soil_drainage_routing.py
Rscript analysis/microbes/11_soil_stream_localization.R
Rscript analysis/microbes/12_soil_branch_signatures.R
```

The full routing took 8.1 seconds after dependencies were loaded, using four
Numba threads. All models were estimable at full rank, and input checksums
remained unchanged. No remote action, commit, or push was performed.

For further inference, prioritize a spatially constrained whole-profile null
and leave-one-soil/segment influence analysis, including sensitivity to the
soil-only concentration thresholds. Recovering the two unresolved spatial soil
records and improving water-column coverage would expand the drainage comparison.

## Distance-only nearness proxy: October 5 follow-up

At Nathan's suggestion, `13_soil_stream_distance.R` also analyzes geographic
nearness without requiring a drainage assignment. It reads the fixed-10K table
and GPS/site metadata directly and includes all 13 mapped soils and all 85
aquatic samples: 260 planktonic, 403 hyporheic, and 442 sediment source–sample
pairs. Soils 04 and the unresolved library remain excluded for missing locations;
catchment-QC failures and drainage-divide uncertainty do not exclude samples
from this distance analysis.

The predictor is straight-line distance from each soil GPS point to each aquatic
sampling site in UTM Zone 10N. Distance to the nearest mapped channel point is
exported separately as a soil-level audit; it is not interchangeable with
soil-to-sampling-site distance. The original-datum assumptions are the same as
above. Euclidean distance is explicitly a proxy for nearness, rather than a
verified contributing watershed or transport route. Read et al. (2015), cited
above, distinguish planar distance, drainage area, and network distance as
different spatial descriptors in river microbial biogeography.

For each soil profile and habitat, the script reports Spearman correlations
between distance and ASV detection fraction or aquatic read percentage.
Descriptive models with soil and aquatic-sample intercepts report the change
per doubling of distance. No routing variables or sequencing plate enter.
The source pools and one-/three-read detection sensitivities follow the earlier
definitions, preserving direct comparison with the drainage analysis.

For the soil-enriched pool, adjusted detection-fraction effects per distance
doubling are +0.0023, -0.0184, and +0.0498 percentage points for planktonic,
hyporheic, and sediment communities. Read-percentage effects are -0.0059,
-0.0106, and +0.0146 points, respectively. These estimates do not show a
consistent aggregate decline in sharing with geographic distance. Individual
soil-profile correlations also vary in direction. This does not negate the
drainage results: the descriptors and retained sampling sets differ.

Individual-ASV descriptions use distance to the nearest mapped soil in which
that ASV is detected, requiring at least three aquatic detections within a
habitat. Aquatic nondetections remain in the correlations. Among 877 eligible
soil-enriched sediment ASVs at the one-read threshold, 513 have negative and
362 positive distance–abundance correlations, with two ties and median rho
-0.0485. This eligibility screen and repeated ASV/source/sample use mean that
directional counts are descriptive, not calibrated tests of localized taxa.
Nearest sampled source does not establish actual source identity.

Constructed coordinate/distance inputs are in
`data/derived/soil_stream_distance_2016/`. Tables are in
`results/diversity_2016/soil_stream_distance/`, including
`soil_aquatic_distance_sharing.csv`, `individual_soil_distance_correlations.csv`,
`adjusted_distance_effects.csv`, and `asv_nearest_soil_distance_correlations.csv`.
The latter retains ASV identifiers and genus/phylum assignments. Models are
full rank, input checksums are unchanged, and no remote action or commit was
performed.

`figures/archive/soil_first_pass/2016_soil_stream_distance_correlations.pdf`
was rendered and visually inspected before being removed from active outputs. Caption: Spearman correlation between geographic distance and the
fraction of each soil's enriched ASVs detected in aquatic samples, separately
by habitat. Negative values indicate greater sharing near the sampled soil.
All mapped soils and aquatic profiles enter, irrespective of drainage QC.
Points have no independent-source confidence intervals or significance claims.

## Catchment-wide soil ASV figure

`14_soil_asv_catchment_figure.R` saves one combined habitat figure to
`figures/2016_soil_asv_catchment_pattern.pdf`. All 85 aquatic samples are shown:
34 sediment (green), 20 planktonic (blue), and 31 hyporheic (orange). Drainage
area is plotted in hectares on a logarithmic axis so that headwaters remain
visible alongside mainstem sites. The y-axis is the percentage of each aquatic
sample's fixed-10K reads assigned to ASVs detected in at least one of the 15
sampled soils. Points are individual samples. Separate Gaussian GAM curves
use log10 drainage area, a thin-plate basis of dimension four, and REML
smoothness selection. Each line stays within its habitat's observed area range.
The small basis permits gentle curvature and approximately linear fits.
See the [mgcv method documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.html).
No independent-sample tests or confidence bands are reported because
samples share segments and the mainstem is sparsely sampled.
This is catchment-wide sharing, irrespective of soil GPS or drainage-membership
QC, and does not estimate a percentage of reads originating in soils.

The script's `y_metric` setting can alternatively show
`soil_pool_richness_fraction`: the number of soil-detected ASVs present in each
aquatic sample divided by that sample's total observed ASV richness at 10K.
Both response columns are retained in `catchment_soil_asv_plot_data.csv` in the
soil-localization results directory. Read totals and detections were checked
against the community matrix, input hashes remained unchanged, and the vector
PDF passed rendered visual inspection.


## Earliest detection along network stages

`15_first_detection_flowpaths.R` adapts the upstream-first inventory in
[Ruiz-Gonzalez, Nino-Garcia, and del Giorgio (2015), Ecology Letters
18:1198-1206](https://doi.org/10.1111/ele.12499), particularly Figure 2.
Their terrestrial-aquatic comparison follows a sequence of ecosystem types.
The source inventory now combines hyporheic and planktonic water detections
as **Aquatic**, while keeping **Sediment** separate. The three receiving
habitats remain separate figure panels. Pooling follows Nathan's interpretation
of strong water-compartment coupling in this system; surface-subsurface water
exchange supplies general hydrological context
([Boano et al. 2014](https://doi.org/10.1002/2012RG000417)). Hydrological
coupling does not by itself demonstrate interchangeable microbial communities.
Sediment and the combined water pool are parallel at each network position,
with simultaneous detections retained rather than imposing a flow order.
Novelty relative to all previous studies has not been established.

The rules use the fixed-10K table and a one-read detection threshold:

1. An ASV detected in any of the 15 soils enters the **Soil** category. All soils
   remain included, including the two without usable mapped locations. This
   is a regional inventory, with no assumed soil-site pairing.
2. For ASVs not detected in soil, find the earliest surveyed aquatic stage:
   **Headwater** (orders 1-2), **Intermediate** (orders 3-4), or **Mainstem**
   (order 5). These are network classes, not a verified connected flowpath
   or observations at successive times.
3. Within that earliest stage, classify detections as **Sediment**, **Aquatic**,
   or **Sediment + Aquatic**. Aquatic means detection in hyporheic water,
   planktonic water, or both, without a sediment detection at that stage.
   Sediment + Aquatic retains a tie between sediment and either water
   compartment. A taxon assigned to Sediment can appear in water downstream
   without changing its category. Habitat labels identify earliest observed
   detection pools, not habitat specificity or demonstrated source identity.

For each sample, richness percentages divide detected ASVs in each category
by that sample's total detected ASVs; read percentages divide category reads
by 10,000. The bars show means of those within-sample percentages rather than
fractions of a pooled group inventory. Each sample has four mutually exclusive
categories summing to 100%. Counts above the bars show the unequal sample
coverage. These percentages represent membership in detection inventories,
not estimated habitat contributions, proven recruitment, movement, or growth.

Unequal coverage can change where a taxon is first observed, even at equal
read depth. A sensitivity analysis repeats the entire classification in
1,000 independently seeded selections of three samples from each of the nine
aquatic habitat-by-stage cells and three soils (30 profiles per draw). It
recomputes percentages only for that draw's selected communities. This
balances the original habitat profile counts; the merged Aquatic inventory
contains six profiles per stage versus three sediment profiles. It does not
balance merged source-pool counts, spatial extent, soil heterogeneity, sampling
completeness, or independent stream segments. Its 2.5th-97.5th percentiles
are descriptive sampling ranges, not confidence intervals. Equal sampling
also reduces detection of rare taxa; the full inventory remains the primary
description. No additional read rarefaction or taxon filtering is performed.

Outputs are in `results/diversity_2016/first_detection_flowpaths/`:

- `sample_network_audit.csv` and `sampling_cells.csv`: identity and coverage.
- `asv_first_detection_assignments.csv`: every ASV's category and first stage.
- `sample_first_detection_profiles.csv` and
  `full_inventory_network_summary.csv`: sample denominators and plotted values.
- `sample_first_stage_profiles.csv`: soil and stage-specific new detections.
- `balanced_inventory_draws.csv`, `balanced_inventory_summary.csv`, and
  `balanced_sample_selections.csv`: all draws, ranges, and selected sample IDs.
- `inventory_diagnostics.csv`, `analysis_settings.csv`, `input_checksums.csv`,
  and `session_info.txt`: computational and input provenance.

Caption for `figures/2016_first_detection_flowpaths.pdf`: Mean within-sample
ASV richness (upper row) and read percentages (lower row) assigned to their
earliest surveyed detection pool, by network stage and aquatic habitat.
Soil-detected ASVs have precedence; other ASVs are assigned to Sediment,
Aquatic (hyporheic and/or planktonic water), or Sediment + Aquatic according
to detections at their earliest sampled network stage. Sediment + Aquatic
retains a tie between sediment and water. Habitat labels do not denote
habitat-restricted taxa. The receiving panels retain the three original
habitats. Counts indicate the
number of sampled profiles. The figure uses all 85 aquatic and 15 soil samples,
without requiring spatial pairing. Detection categories do not establish
source identity or movement. Individual rows are also saved as
`2016_first_detection_richness_percent.pdf` and
`2016_first_detection_read_percent.pdf`. All are vector PDFs.

Run the two current figure scripts directly in Positron or with:

```sh
Rscript analysis/microbes/14_soil_asv_catchment_figure.R
Rscript analysis/microbes/15_first_detection_flowpaths.R
```

Neither requires rerunning the drainage preprocessing. The initial branch
analysis remains a candidate for refinement using soil coverage, verified
tributary topology, shared-segment sensitivity, and habitat-specific responses.

### Current results and validation

The sediment read-percentage smooth is approximately linear (effective
degrees of freedom 1.00); planktonic and hyporheic curves have modest curvature
(edf 1.84 each). Fitted read percentages at 100/1,000 ha are 70.94/69.39%
(sediment), 22.10/20.98% (planktonic), and 31.48/26.98% (hyporheic).
These reference-area differences describe the curves, not formal contrasts
among habitats. Predictions and fit summaries are exported alongside the
catchment plot data.

Mean soil-detected richness fractions for headwater/intermediate/mainstem
communities are 60.02/53.47/57.38% in sediment, 32.17/27.11/34.17% in
hyporheic communities, and 19.15/21.85/13.88% in planktonic communities.
The pattern therefore does not support a universal monotonic downstream
decline. Corresponding soil-detected read percentages are
72.26/67.74/69.90%, 34.69/26.00/33.68%, and 20.03/20.81/13.51%.
The per-sample inventories exactly match the original soil-ASV definition.

Equal-sample inventories retain the same broad nonmonotonic habitat patterns
but substantially reduce observed soil membership. Mean balanced soil-detected
richness percentages are 34.0/29.3/33.0% (sediment), 18.0/14.7/19.1%
(hyporheic), and 10.4/12.0/7.9% (planktonic). Reducing the soil inventory
from 15 to three profiles is one reason that absolute labels change; this
sensitivity balances counts while deliberately reducing inventory coverage.
Thus the absolute percentages depend strongly on which soils were sampled.

The initial five-category analysis completed all 1,001 inventories without
diagnostic flags in 43.5 seconds locally. The pooled-water revision reruns
the same 1,000 selections and retains the original soil assignments and
first-stage ordering. Previous separate-water tables are preserved in
`archive_separate_water_pools/` beneath the results directory. Every sample and group sums to 100% for both
metrics. Independent checks confirmed all 13,969 soil assignments, all 49,260
ASV identities, absence of future-stage assignments in upstream samples, and
three selected samples per cell in each of the 1,000 draws. Input hashes were
unchanged. The catchment figure and all three flowpath PDFs were rendered and
visually inspected. Disabled first-pass figure scripts passed syntax checks.

A useful next branch comparison would use each aquatic site's confidently
contributing soil inventory, compared with equally sampled noncontributing
soil pools, separately by habitat. Source-pool size and observed richness
would need explicit sensitivity checks; sparse soil coverage and uncertain
divide/outlet assignments should remain visible. Drainage area and repeated
segments would require attention before attributing an apparent difference
to locality. This would reuse the Python membership backbone without treating
the two initial branch-concentrated ASV sets as established source signatures.

The current pooled-water categories contain 13,969 Soil, 9,345 Sediment,
23,812 Aquatic, and 2,134 Sediment + Aquatic ASVs. Of the former shared
aquatic category, 2,744 ASVs have water-compartment ties without sediment
and therefore move into Aquatic. Soil-detected percentages are unchanged.

Validation of the pooled-water revision confirmed all 49,260 category
assignments against an independently reconstructed mapping, unchanged
soil percentages and earliest stages, and unchanged selections in all
1,000 sampling draws. Every sample still sums to 100% for both metrics.
All 1,001 inventories completed without flags, input hashes remained
unchanged, and all three revised PDFs passed rendered visual inspection.

## Manuscript status

Nathan identified the pooled-water flowpath figure as an interesting growing
result on October 5, 2026 and requested committing the completed soil work.
Both richness and read-percentage rows remain available; the paper will
likely use one row, with the choice deferred. The catchment plot and sharing
maps provide complementary spatial descriptions. Treat the detection-pool
contrast as descriptive habitat structure rather than source attribution.
