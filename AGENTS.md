# HJA synoptic: working instructions

This repository covers the 2016 HJ Andrews microbial and sediment survey.
Start with `README.md` for the figure workflow and `analysis/METHODS.md` for
interpretation. Consult the actual checkout before relying on historical status.

## Working with Nathan

Read the current synced Dropbox files under `work/style_notes/` when available:
`programming_style.md`, `nathan_writing_style_guide.md`,
`claude_working_notes.md`, and `agents_notes.md`. New explicit user instructions
win. Check in briefly that these project instructions are visible and understood.

- Write clear, mostly linear R with established packages, one package per line,
  a numbered settings block, numbered workflow sections, and explanatory comments.
- Use `here::here()` paths; analysis scripts run directly in Positron without
  arguments. Resource variables may set CPUs, never hidden analytical choices.
- Put analysis and its plots together. Keep the numbered main workflow in
  `analysis/`, input construction in `analysis/preparation/`, and additional
  analyses in `analysis/supporting/`. Avoid separate formatting scripts or wrappers.
- Figures go in `figures/manuscript/` or `figures/supplement/`; tables in
  `results/`; constructed inputs in `data/derived/`. Prefer inspectable CSV/TSV.
- Use white backgrounds, colorblind-safe palettes, no gray facet boxes, and
  vector PDF by default. Dense raster figures use 600 dpi and white backgrounds.
  Keep explanatory titles/footnotes in captions. Save individual panels.
- Give comparable panels consistent, pleasing aspect ratios and enough space;
  size the canvas around panels instead of squashing them. Preserve equal
  coordinate scales for ordinations and maps. Check the rendered figure set.
- Draw regression lines only for slopes passing the stated significance
  threshold, with multiplicity adjustment as appropriate; retain all tested
  slopes and observations in tables.
- Keep documentation short. Do not duplicate numerical output across reports.
  Keep only workflow, input provenance, methodological definitions, and necessary
  interpretation limits. Write American English with precise, calibrated claims.
- Preserve unrelated changes and all raw inputs. Encode translations, exclusions,
  and transformations reproducibly. Keep completed DADA2 outputs and the master
  phylogeny intact. Do not rerun expensive processing without reason and authority.
- Validate readability changes against the same inputs, settings, and seeds.
  Never edit running scripts. Benchmark before proposing cluster computation.
- Retain failed runs and flagged diagnostics. Do not blanket-ignore `results/`.
- Obtain approval of exact files and commit message before committing. Use simple
  descriptive subjects without Conventional Commit prefixes or `Co-Authored-By`.
  Push only with separate authorization.
- Later surveys belong in separate repositories; do not add 2025 workflows here.

## Current analytical decisions

- `analysis/01_prepare_cohort.R` encodes the exact 15 suspect water libraries.
  Their cause remains unresolved; do not call contamination confirmed or exclude
  all Plate 2 libraries. Plate is technical QC, never an ecological covariate.
- Preserve the original 10K draw. The reviewed cohort is 76 samples and 33,974
  ASVs: 15 soil, 27 porewater, 34 sediment. Audits distinguish QC, depth, and scope.
- The original 100-sample source table remains intact for reproducibility and
  deployed iCAMP hashes. New pooled analyses use the reviewed cohort.
- There are 18 sediment–porewater pairs and 33 microbial–molecular sediment
  matches. Site 43 lacks FT-ICR; retain it in sediment-only microbial analyses.
- FT-ICR ecological analyses use presence/absence. Intensity and feature counts
  are detection/QC diagnostics, not calibrated richness. Formulas are not confirmed
  metabolites; network associations do not establish direct interactions.
- Soils are a regional pool. Soil 04 lacks GPS and `hja2016_200` is unresolved.
  First detection does not establish movement or source fractions. Align exact
  sites and balancing selections before comparing microbial and molecular results.
- The mapped stable group's Lookout/Mack pattern is exploratory. Branch effects
  and a metabolic guild have not been established.

## Cluster continuity and authorization

Host: `niw7@lugh`. Checkout: `/mnt/home/niw7/scratch/GitHub/hja-synoptic`
(`scratch` resolves to `/mnt/scratch/wisnoskilab`). Cluster R is 4.4.0; user
packages are under `~/R/x86_64-pc-linux-gnu-library/4.4`. FastTree 2.2.0 is at
`~/scratch/shared/bioinformatics/FastTree/FastTree`; its job loads `openmpi4`.

Before **every** remote action, show the exact command, host, and effect and
wait for explicit approval, including read-only checks, transfers, cleanup,
installation, submissions, and cancellations. Keep each approved action separate.
Before a job, also show its full script, partition, CPUs, memory, walltime, and
expected runtime or state that it is unknown. Prefer `sbatch` for substantial
work; cap arrays at five concurrent tasks. Inspect local/remote status before
approved synchronization, preserve unrelated files, and never use destructive sync.

Sediment job **432228** was submitted October 10, 2026; its current state has
not been checked. Keep `analysis/microbes/06_icamp_sediment.R` and
`analysis/cluster/run_icamp_sediment.sh` unchanged while that deployment stands.
Resources: normal, 16 CPUs, 192 GB, 96 hours; runtime unknown. Logs:
`logs/hja2016_icamp_sediment-432228.out` and `.err`; results:
`results/icamp_2016/sediment_castor/`. Create `logs/` before future submissions.

Catchment job **432103** completed October 9 using the original four-habitat
pool, including now-excluded libraries. It is historical, not a reviewed-cohort
analysis. Earlier incomplete distance matrices must not be reused. Require
matching hashes/tip order, a completion record, and numerical distance checks;
resume checkpoints only from the same run with unchanged settings and versions.
