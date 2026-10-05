# Catchment-wide assembly and sediment microbial–FT-ICR-MS analysis

This analysis treats stream sediment as the focal metacommunity. The primary
microbial dataset contains the 34 sediment libraries retained at 10,000 reads;
the broader sediment biogeochemistry dataset contains 44 sites with FT-ICR-MS,
environmental, and enzyme-activity measurements. Microbe–molecule integration
uses the verified 33-site overlap: site 43 has a retained sediment library but
no molecular profile. The full 34-sample sediment set remains the iCAMP input.

## 1. Microbial assembly at two ecological scales

Use two complementary iCAMP analyses rather than adding an intermediate
aquatic-only analysis.

### 1.1 Catchment-wide habitat assembly

The catchment-wide analysis uses all 100 libraries retained at 10,000 reads:
planktonic streamwater, hyporheic water, stream sediment, and terrestrial soil.
All observed ASVs form one catchment-scale metacommunity pool. Summarize the
pairwise process assignments within each habitat and across each habitat pair.

Run `analysis/microbes/07_icamp_catchment.R` with
`analysis/cluster/run_icamp_catchment.sh` (32 CPUs, 256 GB, 96 hours). It uses
the same 1,000-randomization settings as the sediment run and writes CSV
results to `results/icamp_2016/catchment_castor/`, including descriptive means and
pair counts for each habitat comparison. The observed pool is weighted by
the sampled communities; unequal habitat sample counts therefore influence
its composition. Habitat summaries are reported separately rather than
letting the most numerous comparison dominate the ecological interpretation.

Job 431635 timed out during phylogenetic-distance construction; the null
randomizations have not run. Follow [the rerun plan](README_ICAMP_RERUN.md)
before resubmitting. The settings above describe the intended full analysis.

This analysis tests whether turnover across habitats is more often assigned to
heterogeneous selection, while comparisons within a habitat are more often
assigned to homogeneous selection or dispersal-related processes. These are
hypotheses, not expected outcomes. The choice of regional pool changes the
question posed by a null model, so interpret this run specifically as habitat
filtering from the observed catchment-wide metacommunity
([Chase 2011](https://doi.org/10.1890/ES10-00117.1)). It does not establish
directional movement from soil or water into sediment.

### 1.2 Within-sediment assembly

Run `analysis/microbes/06_icamp_sediment.R` on the cluster with
`analysis/cluster/run_icamp_sediment.sh`. The primary run keeps all 17,946 ASVs
detected in the 34 sediment libraries and uses 1,000 null randomizations,
abundance-weighted bMPD, Bray-Curtis turnover, a minimum bin size of 24, and the
`Confidence` significance index. These settings follow the primary iCAMP
framework, which partitions turnover into heterogeneous selection, homogeneous
selection, dispersal limitation, homogenizing dispersal, and drift/other
processes ([Ning et al. 2020](https://doi.org/10.1038/s41467-020-18560-z)).

Interpret the categories as null-model evidence, not direct observations of
selection, dispersal, or demographic drift. A later sensitivity run should
exclude ASVs occurring in only one sediment sample and test whether the broad
process proportions are stable. The prevalence filter must be declared rather
than silently replacing the primary analysis.

Before treating bin-level selection as ecological mechanism, inspect bin sizes
and test whether niche differences show phylogenetic signal within the distance
and binning scale used by iCAMP. Use the package's `dniche` and `ps.bin`
functions with a small, standardized environmental matrix, then repeat the
analysis with an alternative bin-size or bacterial-only tree if conclusions
depend strongly on the short V4 phylogeny.

After iCAMP finishes, relate pairwise process fractions and bin-level results
to drainage area, distance to the outlet, stream order, geomorphology,
sediment texture, hydraulic conductivity, organic content, and chemistry.
Network position matters because dendritic connectivity constrains dispersal
and changes the balance between local and regional controls
([Brown et al. 2011](https://doi.org/10.1899/10-129.1)).

The sediment-only run is the primary process analysis for subsequent links to
FT-ICR-MS, sediment properties, enzyme activities, and local environmental
gradients. The catchment-wide run supplies broader habitat context; its raw
process percentages should not be treated as directly interchangeable with the
sediment-only results because the two analyses use different metacommunity
pools.

## 2. Sediment FT-ICR-MS patterns

Begin with all 44 sediment FT-ICR sites before integrating microbial data.

1. Summarize detected formula richness, detection-based molecular properties,
   and formula classes. Retain total signal for detection/QC diagnostics only.
   ESI FT-ICR-MS response varies among molecules and requires careful QC
   ([Liu et al. 2020](https://doi.org/10.1021/acsomega.0c01055)).
2. Use presence/absence for molecular ecological analyses, following the
   decision on 2026-10-03. Use Jaccard distances and the documented
   equal-feature sensitivity; raw intensities remain detection/QC inputs.
3. Examine molecular composition against a small, predeclared predictor set:
   network position/geomorphology, grain size and organic content, hydraulic
   conductivity, nutrients, DOM optics/EEM indices, and enzyme activities.
   Add variables with substantial missingness only as sensitivity analyses.
4. Report formula-level summaries by broad molecular class, H/C, O/C, DBE,
   modified aromaticity index, and nominal oxidation state where those fields
   are available. These are molecular signatures, not confirmed compound
   identities or metabolic pathways.

## 3. Paired microbial–molecular analysis

Use the 33 sediment sites shared by the 10K microbial table and FT-ICR data.

1. Compare microbial Hellinger PCA with molecular Jaccard PCoA using
   symmetric Procrustes analysis and `vegan::protest`. Procrustes is an
   interpretable first test of whole-community concordance, but does not by
   itself identify causal taxon–molecule links
   ([McHardy et al. 2013](https://doi.org/10.1186/gb-2013-14-1-r1)).
2. Use RDA and variation partitioning to ask how much microbial and molecular
   composition is associated with shared environmental and spatial/network
   predictors. Keep model size small relative to 33 sites.
3. Join iCAMP pairwise process fractions to FT-ICR dissimilarity and
   environmental/network distances. First visualize whether molecularly or
   environmentally dissimilar site pairs show more heterogeneous selection or
   dispersal limitation; formal pairwise models must use permutations or
   resampling that respect non-independence among site pairs.
4. Only after the global tests, explore stable taxon–molecular-class
   associations. Aggregate ASVs to iCAMP bins or well-supported taxonomic
   groups and FT-ICR features to molecular classes before testing associations.
   Apply false-discovery control to formal tests and label results as signatures
   rather than evidence that a taxon produced or consumed a particular molecule.

Initial composition, chemistry, matched dispersion, and bipartite analyses are
now implemented under `analysis/fticr/`. Their reports and
[key findings](KEY_FINDINGS_2016.md) supersede the earlier prospective portions
of this plan. In particular, the selected stable-link group is exploratory;
individual edges do not have calibrated P values or edge-level FDR claims.

MIMOSA2 is not the first-line tool here. It relies on mapped microbial metabolic
potential and identifiable metabolites, whereas these data are 16S ASVs paired
with thousands of largely formula-level FT-ICR features. Global concordance,
environmental partitioning, and class-level associations are better matched to
the information content of this dataset.

## 4. Potential colonization sources

Use planktonic and hyporheic libraries as same-site source context where they
exist (15 sediment–planktonic pairs, 26 sediment–hyporheic pairs, and 10 sites
with both). Treat terrestrial soils as a regional comparison pool because they
are not paired to the sediment sites.

Start with ASV sharing and Bray-Curtis/UniFrac similarity between each sediment
sample and its available water sources. FEAST can later provide a sensitivity
analysis of source mixtures ([Shenhav et al. 2019](https://doi.org/10.1038/s41592-019-0431-x)),
but incomplete source coverage and the cross-sectional design prevent strong
claims about colonization direction. Source similarity should complement—not
replace—the iCAMP and environmental evidence. In particular, the
catchment-wide iCAMP analysis describes non-random turnover among habitats;
FEAST addresses resemblance to candidate sources. These are distinct
inferences and should be reported separately.
