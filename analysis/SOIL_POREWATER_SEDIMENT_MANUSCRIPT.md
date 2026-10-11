# Soil, porewater, and sediment manuscript

October 10, 2026. This is the current manuscript framing, following Nathan's
review of the anomalous Plate 2 water libraries. It supersedes the
[sediment-only exploration](SEDIMENT_ONLY_MANUSCRIPT.md) as the working outline.
The reviewed workflow is [documented here](microbes/README_SOIL_POREWATER_SEDIMENT.md).

## Central question

How do microbial communities in soils, stream sediments, and hyporheic
porewater differ in composition and sharing across a dendritic network, and
how do sediment microbial and molecular signatures correspond within that
habitat? Soils provide regional terrestrial context; sediment and porewater
provide co-located compartments with different retention and transport settings.
Sediment FT-ICR profiles provide a second compositional measurement within
sediment, rather than a molecular inventory measured separately in porewater.

Working title: **Microbial communities in soils, porewater, and sediments
across a dendritic stream network**.

Nathan's hydrologic interpretation is that porewater is transient. Use that
context to motivate the contrast, while distinguishing it from a residence-time
effect established by these community analyses. Soil sharing does not measure
recruitment, viability, or movement direction. The cross-sectional data do not
identify sediment retention rates or microbial production and consumption of
individual molecular features.

## Reviewed sample set

Use 15 soil, 34 sediment, and 27 hyporheic porewater communities, retaining
76 libraries and 33,974 detected ASVs from the existing fixed 10K draw. Exclude
the 15 exact flagged Plate 2 libraries on suspected contamination or
sample-identity grounds; their cause remains unresolved. The four excluded
hyporheic sites are 73, 57, 76, and 84. All planktonic communities are outside
this manuscript, including the nine unflagged libraries. Sequencing plate
remains technical metadata, never an ecological covariate.

There are 18 sediment-porewater pairs. The original fixed 10K table has 20,
and excluding porewater sites 57 and 84 reduces this to 18; the other two
flagged porewater sites lack retained sediment counterparts. The older count
of 26 is not reproduced in the fixed 10K table. The full pairing audit records
absent libraries, depth exclusions, and review exclusions separately.

Retain all 34 sediment microbial communities and 44 molecular profiles for
their respective analyses, and their identical 33-site intersection for direct
microbial-molecular comparisons. Site 43 remains in microbial-only analyses
and in the sediment-porewater comparison. All 15 soils enter general sharing;
only 13 have resolved mapped locations. Do not impose a single environmental
complete-case filter across analyses.

## Results sequence and candidate figures

1. **Sampling and habitat context.** The new sampling map shows separate soil,
   porewater, sediment microbial, and sediment molecular panels, with the
   33 paired microbial-molecular sites marked. The new Hellinger PCA and Hill
   diversity panels use the reviewed 76-sample cohort. Bray-Curtis habitat
   PERMANOVA R² is 0.202, but dispersion differs strongly; report habitat-
   associated distributions without assigning a specific filtering mechanism.
   Median richness is 2,846 ASVs in sediment, 1,835 in soil, and 726 in porewater.
2. **Soil inventory sharing across the network.** Sediment and porewater have
   median soil-detected read percentages of 70.26% and 34.38%. The paired-site
   difference averages 36.15 percentage points. Neither habitat has a supported
   drainage slope in the segment-mean analysis. Both absolute percentages
   depend on inventory coverage: three-soil draws yield median habitat read
   percentages of 41.58% and 19.43%. These are detection overlaps, not source
   contributions. Retain both read and ASV fractions for figure review.
3. **Same-site sediment-porewater resemblance.** Mean same-site Bray-Curtis
   distance is 0.9190, compared with a mean randomized assignment of 0.9195
   (one-sided P = 0.4131). Weighted and unweighted UniFrac and sediment-detected
   porewater read share also show no supported same-site excess resemblance.
   Within-order reassignments give the same qualitative conclusion. This is
   a failure to detect preferential pairing, not an equivalence test or proof
   that the compartments are disconnected.
4. **Within-sediment network and chemical structure.** Retain the existing
   matched 33-site dispersion figure, showing equal-feature molecular thinning
   explicitly as a sensitivity analysis alongside complete profiles. Retain
   the all-34 microbial/all-44 molecular supporting views. The organic-content
   association and molecular-property adjustment use their own eligible samples.
   The property fraction falls from 4.04% to 2.44% after organic-content/enzyme
   adjustment (BH P = 0.0622); do not repeat the unadjusted association as an
   independent chemical mechanism.
5. **Localized microbial-molecular signatures.** Core 2 remains exploratory.
   Place it after the aggregate sediment analyses, or in supporting material,
   pending the branch-assignment audit. Its selected map pattern does not
   establish a Lookout Mountain branch effect or metabolic guild.

The new [figure captions and result tables](microbes/README_SOIL_POREWATER_SEDIMENT.md)
identify the individual panels and combinations. The existing
[sediment figure report](figures/README_SEDIMENT_FIGURES.md) supplies the
within-sediment panels. These are candidate figures for review, not a finalized
publication layout.

## Porewater interpretation and assembly analysis

The reviewed porewater set has no supported network-group association
(Bray-Curtis PERMANOVA P = 0.2932) or dispersion difference (P = 0.9266).
Eight declared environmental-distance checks show no BH-supported positive
association. Their available sample sizes range from 16 to 27; lack of support
does not establish environmental invariance or a transient-water mechanism.
Potential low-biomass/reagent effects remain unresolved without sequenced
blanks. Do not remove ASVs based solely on genus names from the conversation.

The completed catchment iCAMP includes flagged and planktonic libraries and
does not characterize this reviewed three-habitat null pool. Keep it as a
historical analysis. A new three-habitat iCAMP run would require a distinct pool,
binning, outputs, diagnostics, and separately approved remote actions.

Sediment-only iCAMP job 432228 was submitted October 10 for the unchanged
34-sample sediment set. Its current state has not been checked. Await its
results and diagnostics before deciding whether sediment assembly adds to
the manuscript. Do not change or restart the deployed job for this cohort update.

## Remaining work

Audit the selected core-2 branch hypothesis, review the new comparative figures,
and reconcile the framing with the actual manuscript draft when that file is
available. Verify extraction/sample handling records before assigning a cause
to the excluded libraries. A three-habitat assembly analysis remains planned;
no remote commands, job submissions, commits, or pushes were performed here.
