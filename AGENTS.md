# HJA synoptic project continuity

This repository analyzes the 2016 HJ Andrews LTER microbial and sediment
biogeochemistry survey. Read this file when starting on a new computer, then
inspect the actual checkout and inputs. Paths below are relative to the
repository root unless explicitly identified as machine or cluster paths.

## Working preferences

Preference sources, checked on 2026-10-02, are the synced Dropbox files under
`work/style_notes/`: `programming_style.md`,
`nathan_writing_style_guide.md`, and `claude_working_notes.md`. On the current
computer, locate them under the Dropbox root; use the synced equivalent on
other computers. Read their current versions when available.
This file summarizes the applicable rules, so it remains useful without
Dropbox access. New explicit user instructions take precedence.

- Write clear, mostly linear R scripts using established packages such as
  vegan, iNEXT, phyloseq, and iCAMP. Avoid one-use helper functions, elaborate
  wrappers, and large grids of models before inspecting the data.
- Prefer R with tidyverse, `here`, and patchwork where useful. For new or
  revised scripts, route paths through `here::here()` so they run from any
  working directory within the project. Existing scripts assume the repository
  root; that convention describes the current code, not the preferred future
  path handling. Scripts must run directly in Positron without arguments.
- Load each package on its own line at the top. Use one numbered settings
  block with all analysis choices and seeds, then numbered workflow sections.
  Environment variables may describe machine resources (for example
  `SLURM_CPUS_PER_TASK`), but must not silently change analytical settings.
  Use one statement per line and comments explaining what the code does and why.
- Prefer plain loops when they make the analysis easier to follow. Use design
  tables and per-run seeds for multiple runs; keep failed runs and diagnostics
  as flagged rows rather than silently dropping them.
- Validate readability rewrites against the original outputs with the same
  settings and seeds. Never edit a script while it is running. Benchmark a
  representative run before deciding that new work needs the cluster.
- Put all figure outputs in the top-level `figures/` directory. Write analysis
  tables to `results/` and constructed inputs to `data/derived/`.
- Prefer inspectable CSV/TSV intermediates. Existing large preparation objects
  and detailed iCAMP checkpoints use RDS; these are documented exceptions.
- Keep raw CSVs, workbooks, reads, and completed DADA2 outputs intact. Record
  every sample translation, exclusion, and transformation in reproducible R
  code. Never silently rename samples or restructure source worksheets.
- Ground scientific methods and interpretations in published literature.
  Keep inference proportional to this unbalanced, cross-sectional design.
- Exclude sequencing plate from ecological predictors and conditioning
  covariates, per the user's clarification on 2026-10-03. Retain plate only
  as technical metadata for lab QC and plate-specific DADA2 error learning.
- Use American English in code, variables, documentation, and commit text.
  For manuscript prose, use the writing guide: precise technical vocabulary,
  calibrated claims, and transitions that express actual logic. Avoid generic
  AI phrasing, stacked hedges, filler, and exaggerated claims. The manuscript
  voice is formal; routine progress updates should remain concise.
- Preserve unrelated changes. Before committing, obtain approval of the exact
  files and commit message. Use simple descriptive subjects without
  Conventional Commit prefixes. Push only with user authorization.
- Do not blanket-ignore `results/`. Add only explicit exclusions for unwanted
  files. ASV tables and compact results are needed to regenerate analyses.
- This repository covers 2016. Later surveys and temporal comparisons belong
  in their own repositories; do not add 2025 inputs or workflows here.

## Figure preferences

- Define the plotting theme near the top of the script. Use a white background,
  colorblind-safe palettes, useful gridlines only, and no gray facet-label boxes.
- Keep titles, subtitles, and explanatory footnotes in captions or associated
  notes rather than on the figure itself.
- By default, show regression lines only when the corresponding slope passes
  the stated significance threshold, using multiplicity adjustment when
  appropriate. Retain all observations and report all tested slopes in tables,
  including nonsignificant slopes. Do not add untested trend lines.
- Save single panels; use patchwork for combinations when useful and retain
  the individual panels.
- Default to vector PDF. Use raster for dense plots that would make PDF files
  unwieldy, with 600 dpi and a white background. State the format choice when
  it is not obvious. All figures go in `figures/`.

## Inputs and sample identity

- Raw paired reads: `sequences/` (not tracked by Git).
- Sequence/sample index: `data/hja-synoptic_sequence-sample-list.csv`.
- Aquatic master metadata: `data/hja-env_data_clean.csv`.
- Soil metadata: `data/hja-synoptic_env-data-soils.csv`.
- Raw FT-ICR-MS workbook: `data/hja-FTICRMS.xlsx`.
- Mapping code: `analysis/build_sample_manifest.R`; inspect its outputs in
  `data/derived/`, especially `hja_2016_mapping_issues.csv`.

The original inventory has 120 paired libraries: 24 planktonic, 36 hyporheic,
45 sediment, and 15 soil. There are only 10 complete aquatic triplets; never
assume each site has all habitats. `hja2016_174` translates `WS1_S8` to
`WS1-8`. Soil library `hja2016_200` lacks a resolved site in the sample list;
retain that uncertainty and consult the mapping audit before using its location.
Soils are a regional comparison pool, not paired aquatic-site observations.

## Active workflow and decisions

1. **ASV inference:** `analysis/dada2_pipeline.R` and
   `analysis/cluster/run_dada2.sh`; see `analysis/README_DADA2.md`.
   DADA2 has completed. Separate sequencing plates learn separate error
   models; their sequence tables are merged before chimera removal. Plate is
   a technical distinction, not an ecological grouping. Do not rerun expensive
   steps or overwrite completed outputs without a specific reason and authority.
2. **Environment and FT-ICR preparation:** `analysis/data_prep/README.md` and
   `analysis/data_prep/config.R` document source columns, feature filtering,
   sample joins, and audits. `run_data_preparation.R` writes
   `data/derived/analysis_inputs/`; inspect overlap and missingness audits.
3. **Primary microbial tables:** `analysis/microbes/01_prepare_diversity.R`
   writes `data/derived/microbial_diversity_2016/`. Retain assigned Bacteria
   and Archaea and remove chloroplast/mitochondrial ASVs before analysis.
   All main ecological analyses use the fixed `asv_counts_10k.csv` table;
   use unrarefied data only for declared sequencing/coverage diagnostics.
4. **Diversity and environmental exploration:** `02_basic_diversity.R`
   uses Hill numbers and Hellinger PCA; `04_network_environment.R` plots
   within-habitat patterns. Use Hellinger data for Euclidean PCA/RDA. Start
   descriptively and select a small set of ecological models afterward.
5. **Phylogeny:** `03_build_phylogeny.R` aligns the full screened catalog
   with DECIPHER and builds a FastTree nucleotide tree. The completed master
   tree is `results/phylogeny_2016/asv_tree_screened_fasttree.nwk` (74,782 tips).
   Prune it in memory for each analysis. Avoid repeatedly rebuilding the tree
   or running the prohibitively slow phylogenetic Hill workflow.
6. **UniFrac:** `05_unifrac.R` computes weighted/unweighted UniFrac and PCoA.
   Tables are in `results/diversity_2016/tables/`; figures are in `figures/`.
7. **iCAMP:** `06_icamp_sediment.R` and `07_icamp_catchment.R` analyze two
   pools: sediment only and all four habitats. No aquatic-only run is planned.
   Both use 1,000 randomizations, bMPD, Bray turnover, `Confidence`, bin size
   minimum 24, `ds = 0.2`, phylogenetic randomization within bins, and
   taxonomic randomization across all taxa. With these settings, the pairwise
   result is `icamp_result$CbMPDiCBraya`. Check package arguments and output
   names before changing versions. Detailed checkpoints and CSV exports go
   to `results/icamp_2016/sediment_castor/` or `results/icamp_2016/catchment_castor/`.
   Both scripts now calculate tip distances with castor, validate the disk-backed
   matrix, and pass its descriptor/tip order explicitly to iCAMP. Never reuse a
   distance descriptor without its completion record and matching input hashes.

The fixed 10K table currently has 100 samples and 49,260 ASVs: 20 planktonic,
31 hyporheic, 34 sediment, and 15 soil samples. The sediment subset has 17,946
ASVs. Recheck these counts if the upstream tables change.

## Integration questions

Use `analysis/README_SEDIMENT_INTEGRATION.md` as the detailed design. Analyze
the 44-site sediment biogeochemistry set first, then use the verified 33-site
intersection with the 34 sediment samples retained at 10K. Site 43 lacks a
molecular profile; retain it in sediment-only microbial analyses.
Relate microbial and molecular composition to network position, geomorphology,
sediment properties, DOM optics/EEMs, nutrients, and enzyme activities using
small models, microbial Hellinger ordinations, molecular presence/absence
distances, Procrustes, and variation partitioning. Use molecular incidence only
for ecological analyses; retain raw signal for detection and QC diagnostics.
Composition, chemistry, network dispersion, and exploratory bipartite analyses
are implemented; consult `analysis/KEY_FINDINGS_2016.md` and the dated updates
below. Broader environmental partitioning and iCAMP integration remain planned.

Catchment iCAMP supplies habitat-assembly context; sediment iCAMP supports the
FT-ICR comparison. Different null pools answer different questions, so their
process percentages are not directly interchangeable. Validate phylogenetic
signal and binning before strong mechanistic interpretation. Pairwise process
summaries are descriptive; pairs sharing samples are not independent replicates.

ASV sharing/UniFrac and optional later FEAST provide source context. Available
10K overlaps are 15 sediment–planktonic pairs, 26 sediment–hyporheic pairs,
and 10 sites with both. Cross-sectional resemblance cannot establish movement
direction. Formula-level FT-ICR features are molecular signatures, not confirmed
metabolites, and taxon associations do not demonstrate production or consumption.
MIMOSA2 is deferred given the current 16S/formula-level data.

Method references: [iCAMP, Ning et al. 2020](https://doi.org/10.1038/s41467-020-18560-z),
[null pools, Chase 2011](https://doi.org/10.1890/ES10-00117.1),
[Hellinger, Legendre & Gallagher 2001](https://doi.org/10.1007/s004420100716),
and [FEAST, Shenhav et al. 2019](https://doi.org/10.1038/s41592-019-0431-x).

## Cluster and machine continuity

- Before every command touching a remote machine, show the exact command,
  host, and intended effect and wait for explicit approval. This includes
  read-only checks, transfers, package installation, cleanup, and submission.
- For a job, also show the full script or command, partition, CPUs, memory,
  walltime, and expected runtime (or say that runtime is unknown). Obtain
  approval before submission or cancellation. Prefer `sbatch` for work longer
  than a few seconds; cap any future job array at five concurrent tasks (`%5`).
- Keep each remote action in its own approved call. In particular, separate
  updating code, submitting jobs, and checking status. A push requires its own
  authorization; rsync also requires approval. Use the recorded remote path
  below rather than guessing or creating an unverified location.
- SSH: `niw7@lugh`.
- Cluster checkout: `/mnt/home/niw7/scratch/GitHub/hja-synoptic`;
  `scratch` resolves to `/mnt/scratch/wisnoskilab` on this cluster.
- Local checkout locations vary by machine; use repository-relative paths.
  Check remote/local Git status before synchronizing; raw reads and
  ignored intermediates need separate transfer. Do not use destructive sync.
- Cluster scripts load `R/4.4.0`. User packages are under
  `~/R/x86_64-pc-linux-gnu-library/4.4`. Bioconductor packages such as DECIPHER
  and Biostrings are installed through BiocManager rather than CRAN.
- FastTree: `~/scratch/shared/bioinformatics/FastTree/FastTree`, version 2.2.0,
  OpenMP build; the tree job also loads `openmpi4`.
- Catchment job: `analysis/cluster/run_icamp_catchment.sh`, one node, 32 CPUs,
  256 GB, 96 hours. Sediment job: `run_icamp_sediment.sh`, one node, 16 CPUs,
  192 GB, 96 hours. Create `logs/` before submitting because Slurm opens logs
  before the script runs.

### Local analysis update: 2026-10-03

The local checkout was clean at the start of this session. It lacked the ignored
molecular matrices; `analysis/data_prep/03_prepare_fticr.R` rebuilt them from the
raw workbook and reproduced 60 profiles x 4,760 primary features.

The user specified presence/absence only for FT-ICR-MS ecological analyses.
`analysis/fticr/01_explore_sediment.R` now audits the 44-site subset and saves
Jaccard PCoA, van Krevelen, and detected-feature-count figures. Raw intensities
are used only to define detections and for QC diagnostics. The constructed
binary table is `data/derived/fticr_2016/sediment_44_presence_absence.csv`.
`analysis/microbes/08_sediment_figures.R` makes sediment Hellinger PCA and Hill
q = 1 figures from the fixed 10K microbial table. See `analysis/fticr/README.md`
for methods, captions, diagnostics, and outputs. Tables are in
`results/fticr_2016/tables/` and `results/diversity_2016/tables/`.

Verified paired overlap is **33**, not 34: sediment site `43` is retained at
10K but has no molecular profile. Microbial figures use all 34 sediment samples;
paired analyses refit on the 33-site intersection. Detected molecular
feature counts correlate with total primary signal (Spearman rho = 0.969), so
inspect detection effects before interpreting richness.

`analysis/fticr/02_river_composition.R` completed the initial river-composition
and paired microbial exploration, plus 1,000 uniform draws of 391 detected
features per site as a sensitivity analysis. Tables are in
`results/fticr_2016/river_composition/`; seven additional PDFs are in `figures/`.
See `analysis/fticr/README_RIVER_COMPOSITION.md` for published PNNL methods,
inspected public analysis code, captions, and limitations. Exact HJA extraction
and acquisition settings are not established from the related papers.

PCoA1 tracks detection count (rho = −0.984). Drainage explains 5.15% of complete
Jaccard variation (P = 0.0173), but 2.45% of mean equal-feature dissimilarity
(P = 0.1339). Chemical signatures do not show a convincing monotonic downstream
shift; 33-site molecular–microbial Mantel rho is 0.112 (P = 0.1019), and two-axis
Procrustes correlation is 0.130 (P = 0.8281). These are exploratory results,
not source or process attribution. All draws succeeded and figures passed
rendered visual inspection. Sediment-property models and variation partitioning
remain prospective.

`03_dom_optical_integration.R` now audits surface/hyporheic EEM peaks A/C/T,
fluorescence index, SUVA254, and secondary slope measures (59/56 water profiles;
44/42 molecular matches). `04_integrated_chemistry_pca.R` combines three chemistry
blocks at 42 sites, giving each block equal total variance. All 33 microbial
matches retain both optical compartments. Reports and captions are in
`analysis/fticr/README_DOM_INTEGRATION.md`; new tables are in
`results/fticr_2016/dom_optical_integration/` and `integrated_chemistry/`.

Water optics show a modest aromatic/humic-like versus protein-like contrast,
but no clear drainage gradient or consistent FI corroboration. The first two
FT-ICR property PCs explain 9.84% additional microbial variation, or 4.04% after
adjustment, conditioning on drainage and detection count (P = 0.0103,
FDR-adjusted P = 0.0309). Joint chemistry PCs explain an adjusted 2.04%
(P = 0.0912). These are exploratory associations; source percentages and algal
contributions are not established. Leave-one-microbial-site-out adjusted
fractions for the FT-ICR model range 3.12–4.81%. Eleven new PDFs were rendered
and inspected. Full EEM matrices were not located; HIX, BIX, and PARAFAC are
not reconstructed from peak summaries.

The laptop checkout at that time lacked `05_unifrac.R`, both iCAMP scripts, the master tree,
and `analysis/README_SEDIMENT_INTEGRATION.md` named above. Their historical
completion/status has not been verified from this machine. No remote action,
commit, or push was performed. The new local work is uncommitted; cluster status
below remains historical and requires an approved check.

Local network heterogeneity work is complete in
`analysis/fticr/05_network_dispersion.R`; methods and results are in
`analysis/fticr/README_NETWORK_DISPERSION.md`. Primary tests include **all 44
sites**, split into 24 first/second-order headwaters, 14 third/fourth-order
intermediate sites, and six fifth-order mainstem sites. Six vector figures
are in `figures/`; tables are in `results/fticr_2016/network_dispersion/`.
Equal-feature mean Jaccard dispersion declines across groups (omnibus
BH-adjusted P = 0.003), with headwater/mainstem ratio 1.104. All 1,000 draws
of six distinct segments per group retain higher headwater and intermediate
dispersion than mainstem. Complete profiles, replacement, chemical properties,
and RC do not have significant omnibus tests. RC saturates at -1 in 867/946
pairs; interpret its limited resolution explicitly. Headwaters span a wider
geographic extent, so network position is not isolated causally.
Centerline D is basin-outlet stream distance; prepared
`site_distance_to_outlet_m` is valley-local. All site coordinates and stream
orders match the recorded centerline segments. Null outputs reuse only an
input-checksum/seed/simulation-count/method/tie-rule match. Raw inputs remain
intact; no commit, push, or remote action was performed.

Aquatic microbial/chemical comparisons are complete in
`analysis/fticr/06_aquatic_microbe_chemistry.R`; see
`analysis/fticr/README_AQUATIC_MICROBE_CHEMISTRY.md`. All 20 planktonic samples
match surface-water optics; all 31 hyporheic samples match hyporheic-water optics.
The full 60-profile FT-ICR set matches 19/30 of these communities (planktonic
site 67 and hyporheic site 43 lack FT-ICR). FT-ICR is contextual sediment/site
chemistry, not a separately measured water-compartment formula inventory.
The new 60-site binary and property inputs are in `data/derived/fticr_2016/`;
the original 42-site property-PC reference is preserved and checked numerically.
Adjusted microbial fractions for matched water optics are 4.30%/0.49%, and
for contextual FT-ICR PCs 1.74%/-0.14%; none of four partial RDA tests survives
BH correction. The planktonic FT-ICR hint is sensitive to site removal.
Results and audits are in `results/fticr_2016/aquatic_microbe_chemistry/`;
three vector PDFs are in `figures/`. Habitats are analyzed separately, all
network orders are retained, and rank/within-order/influence checks are exported.

### Local microbial network dispersion update: 2026-10-04

`analysis/microbes/09_network_dispersion.R` now compares network-group
dispersion separately for all 34 sediment, 31 hyporheic, and 20 planktonic
fixed-10K samples. Headwater/intermediate/mainstem counts are 19/10/5,
15/11/5, and 10/7/3, respectively. Methods, captions, results, and deferred
beta-NTI design are in `analysis/microbes/README_NETWORK_DISPERSION.md`.
Seven vector PDFs are in `figures/`; tables and audits are in
`results/diversity_2016/network_dispersion/`.

Bray-Curtis headwater/mainstem dispersion ratios are 0.927 (sediment),
1.015 (hyporheic), and 1.118 (planktonic); all three-group tests are
nonsignificant (BH-adjusted P = 0.751/0.751/0.822). Binary Jaccard
sensitivities and all secondary contrasts are also nonsignificant.
One thousand segment-balanced draws per habitat retain higher sediment
mainstem than headwater dispersion in 97.1% of Bray-Curtis draws, while
planktonic headwater dispersion is higher in 83.4%. Draw fractions are
descriptive sensitivities, not P values; planktonic mainstem has only three
sites. All 6,000 fits succeeded without diagnostic flags. All seven PDFs
were rendered and visually inspected; input checksums remained unchanged.
All microbial network positions and orders match their recorded centerlines.
Beta-NTI remains uncalculated; the master tree is absent locally. No remote
action or push was performed. The full-site microbial and matched sediment
dispersion analyses are committed together with their figures and results.

### Matched sediment dispersion update: 2026-10-04

For direct microbial/molecular network comparisons, use
`analysis/fticr/07_matched_sediment_dispersion.R` and
`analysis/fticr/README_MATCHED_SEDIMENT_DISPERSION.md`. Both datasets use
the identical 33-site intersection: 18 headwater, 10 intermediate, and five
mainstem sites. Site 43 and 11 molecular sites without retained 10K sediment
samples are excluded explicitly. Four composition views are refitted on this
subset, with 1,000 shared segment-balanced site selections for every view.
The established 1,000 molecular feature draws at 391 detections are reused;
the minimum remains 391 in the matched set. Fixed microbial rarefaction is
retained. Original full-site outputs are preserved and checksum-verified.

Matched equal-feature molecular headwater/mainstem dispersion ratio is 1.126
(omnibus raw P = 0.0004, BH-adjusted P = 0.0008 across two primary views).
Sediment microbial Bray-Curtis ratio is 0.927 (P = 0.5174); complete molecular
and binary microbial Jaccard sensitivities are nonsignificant. Molecular
headwater dispersion exceeds mainstem in all shared balanced draws, while
microbial mainstem exceeds headwater in 97.9%. These are sensitivity fractions,
not independent tests; contrasting significance does not itself test a
difference between responses. All 4,000 fits succeeded without diagnostics.
Three new matched-site PDFs are in `figures/`; tables, aligned distance
matrices, selections, and audits are in
`results/fticr_2016/matched_sediment_dispersion/`. All PDFs were rendered and
visually inspected. No remote action, push, or beta-NTI calculation was
performed. This matched analysis is committed together with the full-site
microbial dispersion analysis, figures, results, and documentation.

### Local bipartite network update: 2026-10-04

`analysis/fticr/08_bipartite_feasibility.R`, `09_bipartite_networks.R`, and
`10_bipartite_consensus.R` now implement the sediment ASV-molecular network
workflow. Methods, captions, results, and limits are in
`analysis/fticr/README_BIPARTITE_NETWORKS.md`. Tables are in
`results/fticr_2016/bipartite_networks/`; four vector PDFs and a 600-dpi
association-matrix PNG are in `figures/`. All use the fixed 10K ASV table
and binary molecular data, with 33 paired sites and 26 stream segments.

Requiring presence and absence in at least 10 sites retains 2,405 ASVs and
902 molecular features, collapsed to 2,392/850 unique detection patterns with
full identifier mappings. The primary threshold is phi >= 0.6 for positive
links; negative links use phi <= -0.6 and do not define positive modules.
There are 554 positive and 293 negative links. BRIM fits Barber bipartite
modularity with eight starts; the observed heuristic modularity is 0.7534.
Individual links have no calibrated P values or edge-level FDR claims.

All 999 draws for each of four nulls succeeded (3,996 total). Positive and
negative counts exceed molecular fixed-margin expectations of 367.3/181.0
(both BH-adjusted P = 0.003). Modularity exceeds a fixed-degree bipartite
graph null (mean 0.6350; BH-adjusted P = 0.004). Whole-profile shuffling,
including shuffling within stream order, does not support excess positive
links (P = 0.200/0.313) or modularity. The fixed-margin null disrupts
molecular covariance; its rejection does not establish specific microbial
chemistry coupling. Within-order shuffling does not fully control geography.

Segment-removal refits show unstable broad module boundaries. The exploratory
80% link/co-membership screen retains 261 positive links and one larger
stable-link group containing 38 ASVs, eight molecular features, and 99 links.
Mean all-cross-pair co-membership is 83.1%, but the minimum is 53.8%, so this
is a connected stable-link set rather than a uniformly stable block. The
group spans multiple phyla and is not an established metabolic guild.
Detection-count/drainage partial correlations are descriptive; individual
selected groups have no selection-adjusted significance tests. ASV-abundance
and fully spatially conditional null sensitivities remain possible extensions.

`check_bipartite_networks.R` passed analytical objective and phi checks;
all 36 full-size pilot null draws reproduced in the final run at the same
seeds. All inference and consensus inputs remained unchanged. All five
figures were visually inspected after rendering. No remote action or push was
performed. This network analysis is committed together with stable-group
identities, catchment maps, and key findings. Pilot outputs are explicitly
ignored, while main analysis tables and figures are retained.

### Local stable-group identity and distribution update: 2026-10-04

`analysis/fticr/11_stable_group_distribution.R` describes core 2 from the completed
network analysis, without rerunning inference. Tables are in
`results/fticr_2016/stable_group_distribution/`; three vector map PDFs are in
`figures/`. Coverage is 34 microbial and 44 molecular sites, 45 unique sites and
33 paired. Missing profiles remain distinct from nondetection; collapsed ASV
pattern summaries are restricted to the original paired sites.

Of 38 ASVs, 11 have genus assignments and none has a species assignment.
Bryobacter has three ASVs; eight other named genera/groups have one each.
Acidobacteriota/Pseudomonadota/Gemmatimonadota contain 11/10/5 ASVs. Eight assigned
formulas include three protein-like, three lipid-like, and two
condensed-hydrocarbon-like workbook signatures. Three peaks report 2/4/3 formula
candidates, respectively; no structures or exact metabolites are identified.

Paired-site median group read percentages are 0.08/1.075/0.06% in
headwater/intermediate/mainstem reaches. Representation is patchy, with strong
signals at 190, CC-4, KC-1/KC-3 and several intermediate sites, plus mainstem
49/156. Site 50 has all eight formulas but only four group ASVs. Distribution
summaries are descriptive for a selected group, without spatial significance
tests. All 45 site coordinates and orders match their centerline segments; all
seven input checksums are unchanged. All three maps were rendered and visually
inspected. No remote action or push was performed. These outputs and the key
findings are committed together with the bipartite network analysis.

On October 4, Nathan interpreted the group as common in upper Lookout Creek and
the Mack Creek branch draining Lookout Mountain, upstream of their confluence.
`analysis/KEY_FINDINGS_2016.md` records this map-based interpretation as a priority
follow-up and integrates core 2 as a candidate localized microbial–molecular
result alongside weak whole-composition concordance and the modest property-PC
association. Branch membership and above/below-confluence contrasts are not yet
audited or tested. Follow-up should consider detection count, order/drainage,
shared segments, sediment properties, and group/spatial-hypothesis selection.
Do not describe this observation as a confirmed branch effect or metabolic guild.

### Historical cluster handoff: 2026-10-02

Catchment job **431635** was submitted on 2026-09-30 and verified running on
`lughc2` that day. Installed cluster packages were ape 5.8.1, phangorn 2.12.1,
iCAMP 1.9.1, bigmemory 4.6.6, and vegan 2.7.3. The sediment job has not been
submitted in this session. Current job status has not been checked on October 2.

Next, propose a check of the existing job before submitting a duplicate.
The following commands are a reference for the user; an agent must obtain
approval before executing each remote action:

```sh
ssh niw7@lugh
cd /mnt/home/niw7/scratch/GitHub/hja-synoptic
squeue -j 431635
sacct -j 431635 --format=JobID,State,Elapsed,ExitCode
tail -n 40 logs/hja2016_icamp_catchment-431635.out
tail -n 40 logs/hja2016_icamp_catchment-431635.err
```

Review completed CSV exports and `session_info.txt` before calling the analysis
finished. Retain useful results locally; large phylogenetic-distance files and
verbose logs need deliberate transfer/version-control decisions. The iCAMP,
UniFrac, phylogeny, and integration documentation changes were still uncommitted
at this handoff. Update this dated section after meaningful progress so future
agents can distinguish completed work from plans and historical job status.

### Reconciled checkout and iCAMP status: 2026-10-05

This checkout now includes the nine laptop commits through `f11cdc6`, with
the local UniFrac/iCAMP README additions merged and their scripts/tree retained.
The merge preserves the new presence/absence FT-ICR decision and verified
33-site paired overlap. Older statements about missing local files describe
the laptop at the time of its analysis, not this reconciled checkout.

An approved cluster check found job 431635 timed out on October 4 after 96
hours during phylogenetic-distance construction, before null randomization.
The saved distance matrix is incomplete and must not be reused merely because
its descriptor exists. No final 1,000-draw or sediment job has been submitted.
`analysis/README_ICAMP_RERUN.md` records the implemented faster tip-distance
calculation and full-community, 100-draw cluster pilot. Local castor
benchmarks on 500/2,000/5,000-tip subsets matched ape and pruning checks to
numerical precision. Both scripts passed small local end-to-end iCAMP/export
diagnostics, including cache reuse and incomplete-cache rejection (iCAMP 1.8.6).
The cluster pilot still needs to validate iCAMP 1.9.1 and measure full-size
runtime. Fresh `_castor`
directories leave the failed outputs intact; `test_run <- FALSE` is the default
in both local repository R scripts. The final 1,000-draw job is still pending.
An approved remote check and installation on October 5 confirmed castor 1.8.7,
here 1.0.2, and iCAMP 1.9.1 under cluster R 4.4. Catchment pilot copies were
transferred with approval on October 5 (100 null draws, normal partition,
32 CPUs, 256 GB, 24 hours). Rsync completed successfully; original cluster
scripts have suffix `.before-castor-pilot-20261005-8zLKhk`. The local repository
scripts retain the final-run defaults. Sediment scripts have not been transferred.
Pilot job **432023** was submitted with approval on October 5. Its current
queue/run state has not been checked; do not claim it is running or complete
without an approved status check. Logs are `logs/hja2016_icamp_pilot-432023.out`
and `.err`; outputs belong in `results/icamp_2016/catchment_castor/pilot/`.
Review pilot timing, memory, distance checks, binning, and exports before
proposing final-run resources. Do not modify the active cluster scripts.
Obtain approval before each remote action and before any commit or push of
these local changes.

### Local soil drainage localization update: 2026-10-05

`analysis/microbes/10_soil_drainage_routing.py`,
`11_soil_stream_localization.R`, and `12_soil_branch_signatures.R` analyze
soil ASV sharing and contributing drainage membership. Methods, results,
definitions, and captions are in `README_SOIL_STREAM_LOCALIZATION.md` in that
directory. Frozen public USGS 3DEP terrain and partial gaged-watershed polygons
are in `data/spatial/soil_routing/`, with source/datum assumptions documented.
This current terrain reference is not the original 2008 LiDAR model. The
Python preprocessing uses pysheds 0.5 and NumPy 2.3.5; a pinned requirements
file is provided. Constructed catchments and membership sensitivities are in
`data/derived/soil_drainage_2016/`.

All fixed-10K samples remain in general sharing. Thirteen soils have resolved
identities and GPS; soil 04 lacks GPS and hja2016_200 remains unresolved.
D8 catchments use 30-m outlet snapping, checked at 15/50 m and nine soil
coordinate offsets. Forty-six of 54 distinct aquatic sites pass all area/edge
QC scenarios. Eight failing catchments and 23 unstable soil-site memberships
retain flagged audit rows; 581 stable soil-site pairs enter spatial comparisons.
Never infer contribution by Euclidean nearest-stream distance alone.

Of 13,969 soil-detected ASVs, 6,593 occur aquatically and 4,033 in fifth-order
mainstem samples. A sample/source/proximity-adjusted descriptive comparison
shows modest sediment localization (+1.53 percentage points in soil-ASV
detection fraction; +1.68 for the operational soil-enriched sensitivity),
with small water-column effects and variable individual responses. Two disjoint
tributary source pools above sites 47/66 contain soils 12/19 versus 07/08/09;
soil 10 is excluded from that source contrast for divide sensitivity. Their
soil-only, branch-concentrated signatures contain 225/291 ASVs and are broadly
detected, including downstream mainstem sediment. Site identifiers label the
verified drainage areas; named-creek/mountain attribution is not established
here. These signatures are separate from bipartite core 2.

Tables, source identities, per-ASV contrasts, and checksums are in
`results/diversity_2016/soil_stream_localization/`; 11 vector PDFs are in
`figures/` and were rendered and visually inspected. Input checksums remained
unchanged. This is exploratory spatial sharing, not calibrated transport/source
attribution: pairs share profiles, ASVs, and stream segments. No remote action,
commit, push, DADA2 rerun, or phylogeny rebuild was performed. All work from
this update remains uncommitted.

Nathan subsequently approved geographic soil-to-stream-sampling-site distance
as a proxy for nearness even when drainage membership is uncertain.
`analysis/microbes/13_soil_stream_distance.R` implements this independently of
terrain/routing inputs, retaining all 13 mapped soils and 85 aquatic samples.
Use Euclidean distance for descriptive nearness; do not label it verified
contribution. Soil-to-nearest-channel distance is a separate exported audit.
The descriptive source/site-adjusted distance slopes are small and mixed in
sign; individual profiles and ASVs vary. Constructed inputs are in
`data/derived/soil_stream_distance_2016/`, and results are in
`results/diversity_2016/soil_stream_distance/`. Methods and the additional
vector figure caption are appended to the soil drainage report. The figure
was rendered and inspected, all models are full rank, and input checksums
remained unchanged. This follow-up is also uncommitted; no remote action or
push occurred.

`analysis/microbes/14_soil_asv_catchment_figure.R` adds the requested single-panel
catchment-wide figure: drainage area (ha, log scale) versus soil-detected ASV
reads (%), colored by the three aquatic habitats. It retains all 85 aquatic
samples and the ASV detection pool from all 15 soils. Its settings also allow
the soil-detected fraction of aquatic observed ASV richness. The vector PDF
is `figures/2016_soil_asv_catchment_pattern.pdf`; plot data and input hashes
are in the soil-localization results directory. Counts were verified against
the fixed-10K matrix and the PDF was rendered and inspected. No commit or
push was performed.


### Soil figure review and grouped detection update: 2026-10-05

Following Nathan's review, ASV-example and distance-correlation plots are
disabled by default and archived under `figures/archive/soil_first_pass/`.
The initial branch map/signature plots are also archived pending a better
comparison; retain the Python drainage backbone and all source/membership
tables. Do not reinstate these figures on ordinary reruns. ASV-sharing maps
and the catchment-wide figure remain active. Script 14 now adds separate
Gaussian GAM curves against log10 drainage area (REML, basis dimension four).
They are descriptive, with no independent-sample confidence bands or tests.

`analysis/microbes/15_first_detection_flowpaths.R` adapts Ruiz-Gonzalez et al.
(2015; doi:10.1111/ele.12499). Soil detections have precedence; remaining
ASVs use their earliest observed headwater/intermediate/mainstem stage,
then the compartment(s) detecting them in that stage. The three aquatic
habitats are parallel, with ties retained. Never impose an unverified
sediment-to-hyporheic-to-planktonic flow order. All 15 soils enter this regional
inventory; GPS exclusions apply only to spatial routing. Means of within-sample
richness/read percentages avoid unequal pooled group denominators, but inventory
coverage still differs. A 1,000-draw sensitivity selects three profiles in each
of ten cells (nine aquatic habitat-stage cells plus soils) and reconstructs
each inventory. These are sampling sensitivities, not confidence intervals.

Soil-detected richness percentages across headwater/intermediate/mainstem
are 60.0/53.5/57.4% in sediment, 32.2/27.1/34.2% in hyporheic samples,
and 19.2/21.8/13.9% in planktonic samples. There is no universal monotonic
decline. Three-soil inventories substantially lower absolute percentages.
New tables/assignments/selections are in
`results/diversity_2016/first_detection_flowpaths/`. The combined vector figure
is `figures/2016_first_detection_flowpaths.pdf`; each row is also retained.
Methods, captions, limitations, and checks are in the soil report. All four
revised/new PDFs passed rendered inspection; all inputs remained unchanged.
No remote action, commit, or push was performed; these changes are uncommitted.

Nathan clarified the flowpath labels: use "Sediment first", "Hyporheic first",
and "Planktonic first" rather than "only". These indicate earliest observed
aquatic detection, not habitat specificity. The classification itself is unchanged.


### Pooled-water detection labels: 2026-10-05

The current flowpath legend uses Soil, Sediment, Aquatic, and Sediment + Aquatic.
Nathan requested dropping "first" from habitat labels and combining hyporheic
and planktonic detections into Aquatic, then explicitly chose to retain
sediment-water ties as a fourth category. This supersedes the labels above.
The earliest-stage rule and soil precedence remain unchanged; the three
receiving habitats remain separate panels. Do not interpret Aquatic as
including sediment in this operational legend. Water-only ties now enter
Aquatic; original separate-water tables are archived in the results directory.
Sampling sensitivities continue to balance original habitat-stage cells,
not merged source pools (six water profiles versus three sediment per stage).


### Soil result commit and manuscript candidate: 2026-10-05

Nathan authorized committing the completed soil analysis and identified the
pooled-water flowpath figure as a growing result. Both richness and read
rows are retained; the paper will likely use one row, with the choice deferred.
Key finding 7 records the descriptive habitat contrast and links both panels.
The soil commit excludes the unrelated local iCAMP, UniFrac, and phylogeny
changes, including their portions of shared documentation. No push is authorized.
