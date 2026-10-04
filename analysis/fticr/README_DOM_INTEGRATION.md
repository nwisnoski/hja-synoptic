# DOM optics and integrated chemistry

Exploratory analyses completed locally on 2026-10-03. These workflows extend
the presence/absence molecular analysis with independent water-column DOM
measurements and tests of microbial composition. Run the scripts in order:

```r
source(here::here("analysis", "fticr", "01_explore_sediment.R"))
source(here::here("analysis", "fticr", "02_river_composition.R"))
source(here::here("analysis", "fticr", "03_dom_optical_integration.R"))
source(here::here("analysis", "fticr", "04_integrated_chemistry_pca.R"))
```

Both new scripts run directly without arguments and have been executed from
`analysis/fticr/`. Settings, seeds, and variable choices appear at the top.
Tables are under `results/fticr_2016/dom_optical_integration/` and
`results/fticr_2016/integrated_chemistry/`. Constructed CSV inputs are under
`data/derived/fticr_2016/`; eleven new single-panel vector PDFs are in `figures/`.

## Measurements and interpretation

The local master data contain 59 surface-water and 56 hyporheic-water profiles
with NPOC, SUVA254, fluorescence index, and EEM peaks A, C, and T. All three
primary optical variables are complete within those available profiles. The
site-inclusion audit retains the three missing hyporheic profiles explicitly.
There are 44 surface-water and 42 hyporheic-water matches to the sediment
molecular set. WS3-4 and WS3-5 lack hyporheic optics; neither is in the retained
33-site microbial intersection. Thus, integration of all three chemistry
blocks retains all 33 microbial–molecular sites.

The preparation dictionary links each optical value to its raw master column
and water compartment. The new script verifies missingness and numeric values
against `data/hja-env_data_clean.csv`; no source values are changed.

| Measure | Operational interpretation | Limit |
| --- | --- | --- |
| SUVA254 | Carbon-normalized absorbance related to aromatic character | Aromaticity proxy, not a terrestrial carbon fraction or reactivity measurement |
| EEM A and C | Humic-like fluorescence at the source-reported excitation/emission pairs | Humic-like fluorescence has multiple sources and can be modified by processing |
| EEM T | Tryptophan/protein-like fluorescence | Does not establish algal origin, a protein concentration, or molecular identity |
| T:C and T:A | Relative protein-like versus humic-like fluorescence contrast | Optical ratios, not percentages of molecules or carbon |
| Fluorescence index | Fluorescence spectral-shape indicator with historical terrestrial/microbial source associations | Depends on wavelength definition, corrections, and the fluorescent material measured |
| Source spectral slope ratios | Secondary absorbance-shape measures | One definition is unspecified; the explicit 275–295/300–350 ratio differs from the canonical 275–295/350–400 ratio |

The interpretations follow [Coble (1996)](https://doi.org/10.1016/0304-4203(95)00062-3),
[Weishaar et al. (2003)](https://pubs.usgs.gov/publication/70025134),
[McKnight et al. (2001)](https://doi.org/10.4319/lo.2001.46.1.0038),
[Cory et al. (2010)](https://doi.org/10.4319/lom.2010.8.67), and
[Helms et al. (2008)](https://doi.org/10.4319/lo.2008.53.3.0955).
McKnight's original approximate 1.4/1.9 end members used a 450/500-nm emission
ratio. The HJA column names Cory 2010; we do not apply those original values
as classification cutoffs or calculate source fractions from FI.

We located EEM peak summaries rather than full excitation–emission matrices.
HIX, BIX/freshness, and validated PARAFAC components cannot be reconstructed
from three peak measurements. Correction and instrument records are not
established from the summary table. Both terrestrial and autochthonous organic
matter can occur in freshwater; the available measurements do not uniquely
separate algal material from bacterial products or processed terrestrial inputs.

## Optical contrasts and the river gradient

The primary optical summary uses SUVA254, fluorescence index, and log10(T:C).
NPOC-normalized A, C, and T, log10(T:A), and the source slope ratios are secondary
diagnostics. Peaks A and C correlate closely across sites (Spearman rho = 0.986
in surface water, 0.983 in hyporheic water), so including both raw humic peaks
as separate primary PCA variables would duplicate much of their information.

There is a **modest aromatic/humic-like versus protein-like contrast**:
SUVA254 versus log10(T:C) has rho = −0.406 in surface water (P = 0.0014) and
−0.306 in hyporheic water (P = 0.0216). Fluorescence index does not corroborate
that same contrast: its correlations with log10(T:C) are 0.011 and 0.003.
Consequently, we do not label this a verified terrestrial-to-algal source axis.

The 59-site surface-water PCA standardizes the three primary variables and
explains 49.4% and 31.6% on its first two axes. Surface-water centering and
scaling also project hyporheic observations into the same reference coordinates.
Unusual FI at CC-1 affects the ordinary surface PCA: omitting CC-1 changes the
FI loading on PC1 from positive to negative. The rank-based PC1 correlates
with the ordinary PC1 at rho = 0.801. All observations remain in the primary
analysis; these diagnostics identify sensitivity rather than justify exclusion.

| Primary optical measure | Surface drainage rho; P | Hyporheic drainage rho; P |
| --- | --- | --- |
| SUVA254 | −0.117; 0.376 | −0.055; 0.688 |
| Fluorescence index | 0.093; 0.485 | −0.183; 0.177 |
| log10(T:C) | 0.205; 0.119 | −0.001; 0.997 |

None of the six primary drainage correlations survives BH adjustment; there is
no clear monotonic optical source shift from small to large drainage areas.
The explicit source slope ratio shows an unadjusted surface-water correlation
of 0.293 (P = 0.0246), but BH-adjusted P is 0.319 in the secondary family.
Its wavelength definition also requires caution before comparison with Helms.

Across 56 paired water observations, hyporheic log10(T:C) is higher by a median
0.0617 (approximately a 15% paired ratio increase), Wilcoxon P = 0.0059 and
BH-adjusted P = 0.0177. Neither paired SUVA nor FI differences are compelling.
This is a compartment contrast and does not establish algal production or
downstream change.

## Cross-method agreement and richness sensitivity

The existing molecular Jaccard partition separates replacement from
nestedness-resultant differences following
[Baselga (2012)](https://doi.org/10.1111/j.1466-8238.2011.00756.x).
Across 44 sediment sites, replacement contributes 56.9% of summed Jaccard
dissimilarity and nestedness-resultant differences 43.1%. Neither component
identifies the cause of missing detections. This decomposition complements
the 1,000 equal-feature draws; it does not calibrate analytical sampling depth.

The optical distance is Euclidean across the three primary variables,
standardized using the 59-site surface-water reference. Molecular comparisons
retain presence/absence exclusively. Site labels are matched before all tests;
surface and hyporheic measurements are tested separately rather than counted
as independent replicates from the same site.

| Water compartment | Matched sites | Complete Jaccard rho; P | Replacement rho; P | Mean equal-feature rho; P |
| --- | ---: | --- | --- | --- |
| Surface | 44 | 0.136; 0.0662 | −0.028; 0.6519 | 0.050; 0.3124 |
| Hyporheic | 42 | 0.089; 0.1641 | −0.103; 0.8964 | 0.139; 0.1130 |

These are Spearman Mantel correlations with 9,999 site permutations; none
survives BH adjustment across the six comparisons. The three optical variables
jointly explain an additional 6.0% of total molecular Jaccard inertia in surface
water after conditioning on drainage and detection count (P = 0.1687), and
5.3% for hyporheic water (P = 0.5871). The conditional dbRDA tests use reduced
models. None of fifteen optical-versus-formula-property correlations per
compartment survives its BH adjustment.

Thus, optical DOM measurements describe useful additional water-column
variation, but they provide little independent confirmation of the sediment
molecular pattern. Water optics and sediment extracts measure different pools;
agreement is a testable possibility rather than an expectation imposed on the data.

## Integrated chemistry PCA

The joint analysis uses 42 sites with complete measurements in three blocks:

1. Five presence-based FT-ICR summaries: unweighted mean H:C and O:C, and
   nitrogen-bearing, protein-like, and tannin-like detected-feature fractions.
2. Surface-water SUVA254, fluorescence index, and log10(T:C).
3. Hyporheic-water SUVA254, fluorescence index, and log10(T:C).

Each variable is centered and scaled across these 42 sites. Each block is then
multiplied by 1/sqrt(number of variables), giving each block total variance one.
PCA uses these weighted columns without further scaling. Exported scaling and
block-variance audits make the weighting explicit. Correlated chemical summaries
remain correlated; equal block variance does not make them independent assays.
Feature counts, FT-ICR intensities, drainage, and microbial outcomes are not PCA inputs.

The first two axes explain **30.8% and 22.6%**, or 53.5% together:

- PC1 primarily separates hydrogen-rich/protein-like/N-bearing sediment
  signatures from oxygen-rich/tannin-like signatures. Variable correlations
  include H:C = 0.915, protein-like fraction = 0.967, O:C = −0.943, and
  tannin-like fraction = −0.946. Most optical correlations with PC1 are modest.
- PC2 primarily captures water fluorescence contrasts: surface and hyporheic
  log10(T:C) correlate at −0.906 and −0.868, respectively. SUVA correlations
  are positive but smaller. The FI directions do not consistently fit one
  terrestrial/autochthonous ordering.

The result preserves two complementary dimensions rather than identifying
one source gradient shared by all methods. PCA signs orient H:C toward positive
scores and have no inherent ecological direction. A rank-based joint PCA is
retained as a sensitivity to unusual measurements.

For comparison, a PCA of the five FT-ICR properties alone on the same 42 sites
explains 85.9% and 9.5% on its first two axes. PC1 captures the major
hydrogen-rich versus oxygen-rich contrast. PC2 chiefly separates N-bearing
feature fraction from that contrast (rotation loading = −0.905 for N fraction).
These chemical-property axes answer a different question from a PCoA of all
formula detections.

## Chemistry axes and microbial variation

Microbial responses use Hellinger-transformed screened ASVs from the fixed
10,000-read table at all 33 matched sites. The first two chemistry axes were
specified before microbial testing and learned without microbial outcomes.
Partial RDA conditions on log10 drainage area and molecular detection count.
Sequencing plate is excluded from ecological models per the user's 2026-10-03
clarification; it remains technical metadata for QC. Each model tests two chemistry axes jointly using the same
9,999 reduced-model site permutations.

| Chemistry view | Additional total microbial variation | Adjusted unique fraction | Permutation P | BH-adjusted P |
| --- | ---: | ---: | ---: | ---: |
| FT-ICR property PC1 + PC2 | 9.84% | 4.04% | 0.0103 | 0.0309 |
| Integrated chemistry PC1 + PC2 | 8.09% | 2.04% | 0.0912 | 0.1096 |
| Rank-based integrated PC1 + PC2 | 7.88% | 1.80% | 0.1096 | 0.1096 |

The FT-ICR property model supplies a **modest exploratory association** that
was not apparent in the earlier overall formula-distance comparison. The
individual FT-ICR PC2 test has unadjusted P = 0.0060 and adjusted P = 0.0360
across six marginal axis tests. This supports an association with that fixed
property axis; it does not establish a chemical or microbial mechanism.

Leaving out each microbial site in turn keeps the FT-ICR model's adjusted
unique fraction positive, ranging from 3.12% to 4.81%; all 33 influence fits
succeed. This is effect-size sensitivity, not out-of-sample prediction or
independent replication. The joint model does not strengthen this association
when water-column optics are included. Two chemistry dimensions also leave
46.5% of joint chemistry variation outside the tested model.

The adjusted fraction corrects for model degrees of freedom and conditioning;
it is not the unadjusted variance fraction divided by the residual variance.
Permutation significance and adjusted effect size are reported separately.
Spatial/network dependence remains unresolved, and drainage adjustment does
not itself ensure site exchangeability. The analysis does not establish
microbial production, consumption, source fractions, or algal contributions.

## Figure captions

| File | Draft caption |
| --- | --- |
| `2016_dom_surface_optical_pca.pdf` | PCA of standardized SUVA254, fluorescence index, and log10 EEM peak T:C at 59 surface-water sites. Color indicates drainage area; arrows show scaled variable loadings. CC-1 is labeled because of its unusual optical profile and retained in the analysis. |
| `2016_dom_surface_suva_peak_ratio.pdf` | Surface-water SUVA254 versus EEM peak T:C at 59 sites, colored by drainage area. Higher T:C indicates relatively stronger protein-like than humic-like fluorescence, not a source percentage. |
| `2016_dom_suva254_drainage.pdf` | SUVA254 against drainage area for 59 surface and 56 hyporheic water observations. Shapes and colors identify compartment; overlapping measurements at a site are paired. |
| `2016_dom_fluorescence_index_drainage.pdf` | Source-reported Cory-2010 fluorescence index against drainage area in surface and hyporheic water. No universal source-classification threshold is imposed. |
| `2016_dom_t_c_ratio_drainage.pdf` | EEM peak T:C against drainage area in surface and hyporheic water, with both axes logarithmic. Points represent observations rather than fractions of dissolved carbon. |
| `2016_dom_optical_molecular_protein_comparison.pdf` | Water-column EEM peak T:C versus sediment protein-like detected-feature fraction at 44 surface-water and 42 hyporheic-water matches. Water measurements and sediment extracts represent different organic matter pools. |
| `2016_dom_paired_water_optical_positions.pdf` | Surface and hyporheic optical profiles projected into the surface-water PCA reference at 56 paired sites. Lines connect water compartments at the same site. |
| `2016_fticr_property_pca.pdf` | PCA of five standardized presence-based chemical summaries at 42 sediment sites with complete water optics. Axes retain 85.9% and 9.5% of summary-property variance; color indicates drainage area. |
| `2016_integrated_chemistry_pca.pdf` | PCA of eleven standardized chemistry variables at 42 sites, weighted to give equal total variance to FT-ICR, surface-water optics, and hyporheic-water optics. Axes retain 30.8% and 22.6% of weighted variation. |
| `2016_integrated_chemistry_loadings.pdf` | Correlations of eleven chemistry variables with the first two joint PCA axes. Values describe loadings-derived relationships within these data, not independent source validation. |
| `2016_microbial_pca_integrated_chemistry.pdf` | Independent Hellinger PCA of 33 matched sediment microbial communities at the fixed 10K depth, colored by integrated chemistry PC1. Microbial axes retain 17.3% and 13.3%; the display is separate from the partial RDA test. |

## Verification

Source optical values, sample-key uniqueness, explicit compartment subsets,
44/42 molecular matches, 42 complete chemistry sites, 33 microbial matches,
fixed read totals, and equal block variances passed checks. Numeric intermediates,
inclusion decisions, loadings, coordinates, seeds, influence fits, and session
information are exported. Raw inputs remain intact.

Permutation calculations use **all** microbial PCA coordinates, which preserve
full Hellinger geometry while reducing redundant response columns from 17,764
to 33. Distances agree to 1e-10, and each model is also fitted to the original
Hellinger matrix to verify total, conditioned, and constrained inertia and
adjusted fractions. Maximum inertia discrepancy is 1.4e-12. No microbial
variance is truncated for the test. The original slow local permutation run
was interrupted before changing its script; the completed coordinate-based
run reproduces model quantities and completes each model's tests in about one
second. The
[vegan adjusted-R² documentation](https://vegandevs.github.io/vegan/reference/RsquareAdj.html)
describes the conditional-model adjustment used here.

All eleven PDFs were rendered through Poppler-backed `pdftools` and visually
inspected. These are exploratory associations in one cross-sectional survey;
no cluster actions, commit, or push were performed.
