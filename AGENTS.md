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
- Save substantive results, required inputs, and meaningful diagnostics. Avoid
  separate CSVs for duplicate plot data, settings already declared in scripts,
  or simple summaries of saved tables. Remove outputs of retired workflows.
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

## Review of November 15, 2026

Checks against the saved tables. Figure 3 has been rebuilt; the other figure
changes are pending.

- Equal-feature Jaccard dispersion tracks detected feature count (Spearman
  rho = 0.95 across the 33 sites). Thinning binary detections removes more of a
  richer profile. It is retained in script 04 only as a diagnostic. No
  complete-data view shows a network-group dispersion difference for molecules
  or microbes. Raup–Crick was tried, saturated at zero, and removed.
- Figure 3 now shows the supported pattern on the 33 matched sites: sediment
  organic content declines with drainage area, feature count follows organic
  content, and drainage area adds nothing once organic content is in the model.
  Fixed sediment mass does not standardize carbon, so treat feature count as
  detected features per unit sediment, not chemodiversity.
- The direct feature-count decline with drainage area is the weakest link: it
  passes on the 33 matched sites (P = 0.033) and not on all 44 (P = 0.063).
- Figure 3 (rebuilt November 16) partitions binary Jaccard into replacement and
  richness difference (Podani–Carvalho, `adespatial::beta.div.comp`) for both
  data types. Richness difference is 60% of molecular and 10% of microbial beta
  diversity. Among the 120 flow-connected pairs, the molecular richness-difference
  share rises with watercourse distance (r = 0.48) and is flat among unconnected
  pairs; microbes show neither. Replacement is bounded by 1 minus richness
  difference, so its fall along flow paths is not a separate result. Panels:
  organic content and drainage area, feature count and organic content, and the
  richness-difference share against watercourse distance for each data type.
- Avoid "nested subset" as a finding. Pairwise Baselga nestedness tracked the
  difference in log feature count (r = 0.97), and site NODF ordered by drainage
  area is slightly below a curveball null that fixes site richness and feature
  prevalence (53.3 against 53.7–53.9; terminal check, not in a script). Along
  flow paths the downstream site is the poorer one in 97 of 120 pairs.
- `supporting/10_distance_decay.R` now tests only microbial Hellinger and
  molecular Jaccard, reading pair separations from script 04. Over all pairs,
  decay is weak (Mantel r = 0.12–0.17); among flow-connected pairs molecular
  Jaccard rises with channel distance (r = 0.58). Flow-connected pairs
  concentrate on the five mainstem sites.
- Every figure is built at printed size from `analysis/figure_style.R` (November
  16): black plain text, 9-point titles and 8-point tick and legend labels,
  2.1 x 1.6 in standard panels, 2.6 x 2.0 in panels for overlays, map and
  ordination heights from the data ranges, canvas measured around the panels.
  Composites use flat `wrap_plots(design = ...)` layouts; nested patchworks hide
  fixed panel sizes from the measurement.
- Open question: whether the four-panel Figure 3 is final. The direct
  feature-count decline with drainage area is a supplementary panel.
- Planned, not started: rename `04_sediment_convergence.R`, and restyle the
  scripts for readability.
- The Figure 2B–C porewater–sediment iCAMP contrast does not vary with porewater
  Hill diversity, pre-rarefaction depth, or plate.
- The mapped group is present at 17 sediment sites and nearly absent at 17. This
  is unrelated to plate, depth, feature count, or sampling date. The network's
  association counts do not exceed either profile-shuffle null, so the figure is
  a candidate for the supplement.
- Network position returns no detectable effect in first detection, iCAMP,
  microbial dispersion, porewater composition, or soil sharing. The supported
  results are habitat contrasts on the 18 paired sites. Five mainstem sites
  limit power; a null here does not establish equivalence.
- Porewater now uses the same square marker and blue in Figures 1 and 2 and the
  matched-habitat iCAMP panels. Figure fixes pending: the
  supplementary dispersion panels reuse the habitat colors for network groups;
  Figure 1 plots weighted UniFrac but tests Hellinger; Figure 2A shows the full
  inventory, whose soil share is effort-dependent.

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

Sediment job **432228** has completed results in the local checkout; its saved
environment record ends October 11, 2026. iCAMP scripts are now grouped under
`analysis/supporting/icamp/`: `01_reviewed_pool.R`, `02_sediment_sensitivity.R`,
and `03_sediment_sensitivity_figures.R`. Their contents were preserved during
relocation; local launchers use the new paths. The remote checkout was not
checked or synchronized during this cleanup. Never edit a script while it is running.
Resources: normal, 16 CPUs, 192 GB, 96 hours; runtime unknown. Logs:
`logs/hja2016_icamp_sediment-432228.out` and `.err`; results:
`results/icamp_2016/sediment_castor/`. Create `logs/` before future submissions.

Catchment job **432103** completed October 9 using the original four-habitat
pool, including now-excluded libraries. It is historical, not a reviewed-cohort
analysis. Earlier incomplete distance matrices must not be reused. Require
matching hashes/tip order, a completion record, and numerical distance checks;
resume checkpoints only from the same run with unchanged settings and versions.
