# Early script readability and plotting validation

Validated locally on October 6, 2026, against the clean checkout at
`cca972bb9f2fde7e32cadbb5deb2546c65ef8e2d`. This pass applies the current
Dropbox programming guide and the repository rules to the rerunnable early
preparation and microbial diversity workflows.

## Changes

Thirteen R files were revised: `build_sample_manifest.R`; the five numbered
`data_prep/` steps, their configuration, helpers, and runner; and microbial
scripts `01_prepare_diversity.R`, `02_basic_diversity.R`,
`04_network_environment.R`, and `05_unifrac.R`.

Packages are loaded separately at the top, settings are declared before the
workflow, and numbered sections follow the preparation or analysis in order.
Paths resolve through `here`, allowing direct sourcing in Positron from a
nested project directory. Exact-key matches replace per-site lookup wrappers;
plain loops expose missingness summaries, feature exclusion reasons, and iNEXT
inputs. Shared helpers remain where they consistently preserve source headers,
write audit tables, or check identities. The preparation runner sources its
steps in separate environments and stops on errors.

The manifest retains its legacy path flags because the existing DADA2 cluster
launcher uses `--project-root`. Its relative FASTQ paths, explicit `WS1_S8`
translation, unresolved soil identity, and incomplete aquatic triplets are
preserved. Preparation paths and molecular filters remain in `data_prep/config.R`;
the preparation scripts no longer require command-line setup.

The three early plotting scripts follow the later microbial palette and theme:
planktonic blue, hyporheic orange, sediment green, and a distinct vermilion soil
color; 12-point text; white backgrounds; no gray facet boxes; and white vector
PDF exports. Only gradient plots retain light horizontal guides. Ordination
axes use equal physical scales. Figure jitter has its own fixed seed without
changing count-table rarefaction or diversity estimates. Short alpha labels
avoid the clipped labels observed during rendering. Explanations and method
names are recorded in the microbial README and filenames, rather than overall
titles or footnotes. Individual panels accompany the existing overview
filenames; Hill and UniFrac combinations use patchwork.

## Validation

The original and revised scripts ran on separate local copies of the actual
inputs. Raw CSVs/workbook, reads, completed DADA2 outputs, and the master tree
were shared read-only inputs; preparation and analysis outputs went to separate
temporary directories. The original scripts ran from the project root. The
revised scripts were sourced with `source(here::here(...))` from
`analysis/microbes/`, exercising the interactive execution path.

All **59 generated CSV, compressed CSV, and RDS outputs matched exactly**.
This includes sample translations, missingness and exclusion audits, molecular
matrices, fixed rarefied and unrarefied ASV tables, sample metadata, Hill/iNEXT
summaries, Hellinger scores/RDA summaries, and weighted/unweighted UniFrac
distances and PCoA scores. Compressed CSVs were compared after decompression;
RDS files matched byte for byte. The source inventory comparison excluded only
`absolute_path` and `modified_time`, which differ when files are regenerated
in separate directories. Its relative paths, sizes, and input MD5 hashes matched.
Machine-specific session logs and intentionally restyled PDF bytes were not
included in the numerical equivalence comparison.

The existing stored primary ASV count tables and all eleven early diversity
result tables also matched the rerun of the original scripts. Current analysis
tables were left intact. All 175 protected raw-source/workbook and completed
DADA2 files retained their initial MD5 hashes.

Additional checks:

- The legacy manifest path flags reproduced all six mapping CSVs.
- With the DADA2 completion marker absent, the runner completed steps 01–03,
  recorded steps 04–05 as `waiting_for_dada2`, and did not read a deliberately
  invalid partial ASV file.
- The unchanged public `00_load_analysis_data.R` loader accepted the rebuilt
  inputs with source checksum verification, DADA2 requirements, and its default
  in-memory views enabled.
- All modified R sources parsed, and `git diff --check` passed.
- All 42 vector PDFs were rendered and visually inspected, including the
  separate panels and existing overview filenames.

The original full workflow took about 224 seconds locally. Revised timings
include the additional individual-panel exports; these are single runs, not
replicated performance benchmarks.

Inspect the retained validation tables:

- [Every output comparison](../results/script_clarity_2016/output_comparison.csv)
- [Original/revised source hashes](../results/script_clarity_2016/source_comparison.csv)
- [Full-size runtimes and exit codes](../results/script_clarity_2016/runtime_comparison.csv)

R was 4.6.1; package versions were here 1.0.2, readxl 1.5.0, vegan 2.7.6,
iNEXT 3.0.2, ggplot2 4.0.3, ape 5.8.1, phangorn 2.12.1, phyloseq 1.56.0, and patchwork 1.3.2.

## Added q = 2 network figures

The October 6 follow-up adds `2016_alpha_q2_network_gradients.pdf` and six
individual habitat-by-gradient panels. These reuse the existing q = 1 plotting
layout with the fixed-10K table's already calculated q = 2 response. The rerun
preserved the alpha input, alpha-network table, and within-habitat PCA score
table checksums. The panel data update uses ggplot's supported `+ data.frame`
syntax rather than its deprecated `%+%` operator. The source hash record includes
this addition; the runtime table describes the preceding readability run.

## Drainage-focused plotting follow-up

The later October 6 request changes the network figures to drainage area only
and introduces explicit slope tests to govern line display. This is an added
analysis rather than a readability-only change. The settings, segment-clustered
inference, three-habitat BH test families, complete estimates, and limitations
are documented in the microbial README. The original equivalence report above
remains a record of the readability pass; the new two model-audit CSVs extend
that output set. The fixed count tables, metadata, alpha estimates, original
alpha-network table, and PCA score table retained their checksums on the new
run. All nine models succeeded, with no exclusions and no displayed lines.
The three drainage overviews and nine individual drainage panels were rendered
and inspected. The source-hash record includes this later plotting change;
the runtime table remains the preceding readability benchmark.

The subsequent figure preference is recorded in `AGENTS.md`: show regression
lines only after the corresponding slope passes the stated significance
threshold. Untested sediment EEA fits and their bands were removed from the
overview and four individual panels; their observations are retained.
All five EEA PDFs were regenerated, rendered, and visually inspected.

## Scope retained

DADA2 inference and QC utilities, phylogeny construction, cluster launchers,
iCAMP, completed later ecological analyses, and the historical mothur/R Markdown
workflow were not rewritten or rerun. The public data loader remains a reusable
function interface; its helpers serve repeated loading and validation. No remote
action, commit, or push was performed.

Analytical settings and ecological interpretation were unchanged by the
readability pass; the later drainage test/display settings are declared above. Method
references and figure notes are in
[the microbial README](microbes/README.md#early-diversity-figure-notes-and-method-references),
including [Hill-number rarefaction and extrapolation](https://doi.org/10.1890/13-0133.1),
[Hellinger ordination](https://doi.org/10.1007/s004420100716), and
[UniFrac](https://doi.org/10.1128/AEM.71.12.8228-8235.2005).
