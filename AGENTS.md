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
   to `results/icamp_2016/sediment/` or `results/icamp_2016/catchment/`.

The fixed 10K table currently has 100 samples and 49,260 ASVs: 20 planktonic,
31 hyporheic, 34 sediment, and 15 soil samples. The sediment subset has 17,946
ASVs. Recheck these counts if the upstream tables change.

## Integration questions

Use `analysis/README_SEDIMENT_INTEGRATION.md` as the detailed design. Analyze
the 44-site sediment biogeochemistry set first, then join the 34 sediment
samples retained at 10K; verify the exact overlap rather than assuming it.
Relate microbial and molecular composition to network position, geomorphology,
sediment properties, DOM optics/EEMs, nutrients, and enzyme activities using
small models, microbial Hellinger ordinations, molecular presence/absence
distances, Procrustes, and variation partitioning. Use molecular incidence only
for ecological analyses; retain raw signal for detection and QC diagnostics.
Much of this integration is planned rather than implemented.

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

This local checkout lacks `05_unifrac.R`, both iCAMP scripts, the master tree,
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
