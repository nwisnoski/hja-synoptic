# HJA synoptic project continuity

## Current manuscript cohort: October 10, 2026

Nathan authorized removing 15 anomalous Plate 2 water libraries for suspected
contamination or sample-identity problems and focusing on soils, hyporheic
porewater, and sediments. Their cause is unresolved; do not report confirmed
contamination or exclude all Plate 2 samples. Exact identities and metadata
assertions are in `analysis/microbes/18_prepare_manuscript_cohort.R`.
It subsets the existing fixed 10K draw without new rarefaction, writing
`data/derived/soil_porewater_sediment_2016/`: 76 libraries, 33,974 ASVs,
15 soil / 27 porewater / 34 sediment. The 120-library inclusion audit separates
suspect-library, read-depth, and manuscript-scope exclusions. Nine unflagged
planktonic libraries are outside manuscript scope, not labeled contaminated.
All 15 soils remain in general sharing, including the two without resolved GPS.

Run `analysis/microbes/19_soil_porewater_sediment.R` for reviewed habitat PCA,
Hill diversity, soil sharing and inventory sensitivity, same-site comparisons,
and porewater environmental checks. Methods and captions are in
`analysis/microbes/README_SOIL_POREWATER_SEDIMENT.md`; current framing is in
`analysis/SOIL_POREWATER_SEDIMENT_MANUSCRIPT.md`. Results are under
`results/diversity_2016/soil_porewater_sediment/`; 12 vector PDFs use prefix
`2016_soil_porewater_sediment_` in `figures/`. The prior four-habitat outputs
remain historical and must not supply current pooled/water manuscript claims.
The old source draw is retained to preserve exact rarefaction and deployed
iCAMP hashes; new manuscript analyses should use the reviewed cohort.

The verified fixed-10K sediment-porewater overlap is 20 before review and 18
after excluding porewater sites 57 and 84. The older count of 26 is stale for
this input. All 34 sediment communities, 44 molecular profiles, and their
33-site intersection remain unchanged. Site 43 remains eligible for microbial
and paired porewater analyses. The completed four-habitat catchment iCAMP
includes the suspect block and planktonic samples and does not represent this
cohort. A new three-habitat iCAMP run remains planned, with separate null pool,
binning, outputs, and remote authorization. Sediment job 432228 and its deployed
scripts/inputs are untouched; its current state has not been checked.

The six planned manuscript figures and captions are in
`analysis/MANUSCRIPT_FIGURES.md`. Script 15 now uses only reviewed soil,
hyporheic, and sediment inputs; script `analysis/figures/02_current_manuscript.R`
plots habitat Hill diversity and weighted UniFrac refitted on all 76 libraries.
Obsolete plankton-inclusive scripts and 47 water/mixed figures were removed
from the working tree at Nathan's request. Their checksums and Git recovery
commit `908d201` are in `analysis/history/removed_plankton_workflows_20261010.csv`.
Historical reports retain original results; do not run removed scripts from
historical command examples. Common source preparation and saved distances,
raw inputs, completed inference, and deployed sediment scripts remain intact.

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
3. **Original microbial tables:** `analysis/microbes/01_prepare_diversity.R`
   writes `data/derived/microbial_diversity_2016/`. Retain assigned Bacteria
   and Archaea and remove chloroplast/mitochondrial ASVs before analysis.
   This fixed `asv_counts_10k.csv` table is the source draw. Current manuscript
   analyses use its reviewed subset in `soil_porewater_sediment_2016/`, built
   by script 18 above. Use unrarefied data only for declared diagnostics.
4. **Diversity and environmental exploration:** `19_soil_porewater_sediment.R`
   uses the reviewed cohort for Hill numbers, Hellinger PCA, and environmental
   checks. The original mixed-habitat scripts 02/04 were removed. Use Hellinger data for Euclidean PCA/RDA. Start
   descriptively and select a small set of ecological models afterward.
5. **Phylogeny:** `03_build_phylogeny.R` aligns the full screened catalog
   with DECIPHER and builds a FastTree nucleotide tree. The completed master
   tree is `results/phylogeny_2016/asv_tree_screened_fasttree.nwk` (74,782 tips).
   Prune it in memory for each analysis. Avoid repeatedly rebuilding the tree
   or running the prohibitively slow phylogenetic Hill workflow.
6. **UniFrac:** `05_unifrac.R` computes weighted/unweighted UniFrac and PCoA.
   Source tables are in `results/diversity_2016/tables/`; reviewed figures are
   made by `analysis/figures/02_current_manuscript.R` in `figures/`.
7. **iCAMP:** `06_icamp_sediment.R` implements the sediment pool. The original
   four-habitat script 07 is historical (removed; Git `908d201`); it completed as job 432103 on
   October 9 and is historical for the reviewed manuscript. Sediment job
   432228 was submitted October 10; its state has not been checked. A reviewed
   three-habitat pool is planned but has not been run.
   Both use 1,000 randomizations, bMPD, Bray turnover, `Confidence`, bin size
   minimum 24, `ds = 0.2`, phylogenetic randomization within bins, and
   taxonomic randomization across all taxa. With these settings, the pairwise
   result is `icamp_result$CbMPDiCBraya`. Check package arguments and output
   names before changing versions. Detailed checkpoints and CSV exports go
   to `results/icamp_2016/sediment_castor/` or `results/icamp_2016/catchment_castor/`.
   Both scripts now calculate tip distances with castor, validate the disk-backed
   matrix, and pass its descriptor/tip order explicitly to iCAMP. Never reuse a
   distance descriptor without its completion record and matching input hashes.

The original fixed 10K table has 100 samples and 49,260 ASVs: 20 planktonic,
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

The completed original catchment iCAMP describes the unreviewed four-habitat
pool and cannot supply current cohort assembly fractions. Sediment iCAMP
is the unchanged FT-ICR-linked analysis awaiting results. Different null pools answer different questions, so their
process percentages are not directly interchangeable. Validate phylogenetic
signal and binning before strong mechanistic interpretation. Pairwise process
summaries are descriptive; pairs sharing samples are not independent replicates.

ASV sharing/UniFrac and optional later FEAST provide source context. Verified
original 10K overlaps are nine sediment-planktonic pairs, 20 sediment-porewater
pairs, and six sites with both. The reviewed manuscript uses 18 sediment-
porewater pairs and no planktonic communities. Cross-sectional resemblance cannot establish movement
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
- Historical catchment job: one node, 32 CPUs, 256 GB, 14 days; its launcher
  was removed from the local checkout. Sediment job: `run_icamp_sediment.sh`, one node, 16 CPUs,
  192 GB, 96 hours. Create `logs/` before submitting because Slurm opens logs
  before the script runs.

## Dated project history

Detailed October 2-10 preparation, network, soil, and cluster handoffs are
preserved in [the project history](analysis/history/PROJECT_HISTORY.md). Those
records retain historical states; use the current cohort and deployment
summaries above for new work.
