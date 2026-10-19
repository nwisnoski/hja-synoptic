# Rebuilding inputs

The completed inputs are already available. Use the main scripts in `analysis/`
for figures; rebuild preparation stages only when their inputs actually change.
Raw CSVs, workbooks, reads, completed DADA2 outputs, and the master tree stay intact.

## Sequence and microbial inputs

Run these stages in numerical order when regenerating from raw reads:

| Script | Purpose and output |
| --- | --- |
| `01_sample_manifest.R` | Join the sequence index, environmental metadata, and soils; audit identity in `data/derived/hja_2016_mapping_issues.csv`. |
| `02_dada2.R` | Infer ASVs using `config/dada2_run_config.csv`; completed outputs in `results/dada2_2016/`. |
| `03_dada2_summary.R` | Read retention and taxonomy QC. |
| `04_microbial_tables.R` | Screen Bacteria/Archaea, remove organelles, and construct fixed-10K tables in `data/derived/microbial_diversity_2016/`. |
| `05_phylogeny.R` | DECIPHER alignment and FastTree of the screened catalog; completed tree in `results/phylogeny_2016/`. |
| `06_unifrac.R` | Prune/root that tree in memory; save normalized weighted and unweighted distance tables in `results/diversity_2016/tables/`. |

DADA2 learns errors by technical sequencing plate and merges sequence tables
before chimera removal. `check_dada2_setup.R` and `plot_quality_profiles.R`
are QC utilities. Large reads under `sequences/` and SILVA reference files under
`reference/` are untracked; see [taxonomy references](../../reference/README.md).
The master tree has 74,782 tips. Do not reconstruct it for individual figures.
The original 100-sample 10K table remains the frozen source for the reviewed
cohort and deployed iCAMP inputs; it is not the current manuscript cohort.

## Environment and molecular inputs

```r
source(here::here("analysis", "preparation", "environment", "run_data_preparation.R"))
```

The five numbered steps inventory sources, prepare environment, filter FT-ICR,
align sediment blocks, and build source-pool tables. `environment/config.R`
defines source columns, units, mappings, and filters; `helpers.R` contains shared
preparation checks. Outputs and audits go to `data/derived/analysis_inputs/`.
The last two steps require the completed DADA2 marker. The FT-ICR primary table
has 60 profiles and 4,760 features; its inclusion decisions are saved explicitly.
`load_analysis_data.R` is an optional reader of these completed blocks.

## Supporting analyses

The numbered scripts in `analysis/supporting/` have this dependency order:

- 01–02 prepare sediment molecular profiles, class fractions, and the saved
  equal-feature Jaccard reference. These supply main figures 3–4.
- 03 examines full-site molecular dispersion and supplies the network-position
  audit.
- 04–07 prepare, fit, assess stability, and describe the bipartite network.
  They supply main figure 4. `check_bipartite_networks.R` checks the methods;
  `bipartite_helpers.R` is shared by inference and its checks.
- 08 constructs putative molecular transformation profiles; 09 independently
  checks the curated motifs and creates the bulk summaries.
- `soil_drainage_routing.py` reproduces optional drainage membership; its
  requirements are in `analysis/supporting/soil_routing_requirements.txt` and terrain
  provenance is [documented with the inputs](../../data/spatial/soil_routing/README.md).
- `icamp/` groups the reviewed-pool analysis and sediment-only sensitivity.
  See the [figure guide](../../figures/supplement/icamp/README.md).

| Script in `analysis/supporting/icamp/` | Role |
| --- | --- |
| [01_reviewed_pool.R](../supporting/icamp/01_reviewed_pool.R) | Primary exploratory view: validates and reuses the completed reviewed-pool result, then plots; can start inference if the result is absent. |
| [02_sediment_sensitivity.R](../supporting/icamp/02_sediment_sensitivity.R) | Sediment-only inference; retained byte-for-byte from the completed deployment and stops if its completed result exists. |
| [03_sediment_sensitivity_figures.R](../supporting/icamp/03_sediment_sensitivity_figures.R) | Maps and comparisons from completed sediment-only assignments; does not rerun inference. |

Reuse completed results when plotting. Network inference and sequence processing
are substantial computations. `analysis/cluster/` contains batch launchers;
all remote actions require explicit approval. The local iCAMP launchers use the
paths above; script moves do not synchronize the remote checkout.
A distance-cache descriptor is reusable only with matching input hashes, tip
order, a completion record, and numerical validation. Resume checkpoints only
within the same run with unchanged inputs, settings, and package version.
The completed original catchment run used all four habitats and is historical;
the completed reviewed three-habitat run is the primary exploratory iCAMP view.
