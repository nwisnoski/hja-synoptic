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
small models, Hellinger ordinations, Procrustes, and variation partitioning.
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

### Last recorded handoff: 2026-10-02

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
