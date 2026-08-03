# 2016 DADA2 ASV workflow

This is the active analysis path for the 2016 HJA synoptic 16S V4 reads. Run
each numbered checkpoint from the repository root. Do not begin with the old
mothur shared table: DADA2 starts from the paired FASTQ files in `sequences/`.

The two sequencing plates are treated as separate runs for filtering, error
learning, and denoising. Their exact sequence tables are merged only afterward,
before consensus chimera removal.

## One-time setup

Open the `hja-synoptic.Rproj` project in RStudio, or set the working directory
to the repository root. Confirm that DADA2 and ggplot2 are installed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("dada2")
install.packages("ggplot2")
```

To run a script from the R console in a clean R process:

```r
run_script <- function(path, args = character()) {
  status <- system2(file.path(R.home("bin"), "Rscript"), c(path, args))
  if (!identical(status, 0L)) stop("Step failed: ", path)
}
```

The equivalent `Rscript ...` commands can be pasted directly into the RStudio
Terminal.

## Fresh cluster run

The Slurm entry point is:

`analysis/cluster/run_dada2.sh`

It is configured for the remote repository at
`/mnt/home/niw7/scratch/GitHub/hja-synoptic` and requests one node, 8 CPUs,
128 GB RAM, and 48 hours. From that repository:

```sh
mkdir -p logs
sbatch analysis/cluster/run_dada2.sh
```

The job assumes a fresh output path. It rebuilds the 2016 manifest, checks that
all 240 FASTQs and required R packages are present, runs DADA2, and creates the
QC summaries. It refuses to write into an existing output directory unless
replacement of the default output is explicitly authorized:

```sh
sbatch --export=ALL,DADA2_OVERWRITE=1 analysis/cluster/run_dada2.sh
```

That command permanently removes only
`results/dada2_2016` and then starts the complete workflow from the raw FASTQs.
It does not reuse any intermediate files from the interrupted run. Overwrite is
deliberately restricted to that exact default output path.

By default the output is `results/dada2_2016`. To preserve a separate attempt:

```sh
sbatch --export=ALL,DADA2_OUTPUT_DIR=results/dada2_2016_attempt2 \
  analysis/cluster/run_dada2.sh
```

If the SILVA genus training file is present under `reference/`, the job assigns
taxonomy; otherwise it completes ASV inference without taxonomy. Slurm output
and error logs are written under `logs/`.

## 1. Rebuild and review the sample manifest

```r
run_script("analysis/build_sample_manifest.R")
```

Review:

- `data/derived/hja_2016_mapping_summary.csv`
- `data/derived/hja_2016_mapping_issues.csv`
- `data/derived/hja_2016_sequence_crosswalk.xlsx`

Expected audit:

- 120 paired libraries: 24 planktonic streamwater, 36 hyporheic water,
  45 stream sediment, and 15 terrestrial soil.
- 118 direct metadata matches.
- `hja2016_174`: `WS1_S8` is normalized to `WS1-8`.
- `hja2016_200`: retained in DADA2 as a soil library, but its site remains
  unresolved.

Only 10 of the 59 represented aquatic sites have the complete
streamwater/hyporheic/sediment triplet. This is an unbalanced design, not a
mapping failure.

## 2. Check the local setup

```r
run_script("analysis/check_dada2_setup.R")
```

This is a fast, read-only check of packages, FASTQ pairing, run configuration,
available disk space, and the manifest. It writes:

`results/dada2_2016/setup_check.csv`

Do not start DADA2 if the final `ready_for_quality_review` row is `FALSE`.

## 3. Generate and inspect quality profiles

```r
run_script("analysis/plot_quality_profiles.R")
```

Open both PDFs:

- `results/dada2_2016/runs/HJA2016_Plate1/quality_profiles_input.pdf`
- `results/dada2_2016/runs/HJA2016_Plate2/quality_profiles_input.pdf`

The script selects samples across each run rather than only the first files.
The inspected raw reads are 250 bp and the 515F/806R primers are already absent,
so the current configuration uses `trim_left_f=0` and `trim_left_r=0`.

Review the quality tails, then edit only:

`config/dada2_run_config.csv`

The supplied `trunc_len_f=240` and `trunc_len_r=200` are starting values.
Because V4 is roughly 250 bp, those values normally retain ample overlap for
paired-read merging. Keep the two plate rows separate even if their settings
end up identical.

This is the deliberate pause in the workflow. The next step is computationally
expensive and should begin only after the PDFs have been reviewed.

## 4. Infer ASVs with DADA2

Without taxonomy:

```r
run_script("analysis/dada2_pipeline.R")
```

The default command uses:

- manifest: `data/derived/hja_2016_sample_manifest.csv`
- run settings: `config/dada2_run_config.csv`
- output: `results/dada2_2016`
- random seed: `2016`

To assign SILVA 138.2 taxonomy in the same run:

```r
run_script(
  "analysis/dada2_pipeline.R",
  c(
    "--taxonomy", "reference/silva_nr99_v138.2_toGenus_trainset.fa.gz",
    "--species", "reference/silva_v138.2_assignSpecies.fa.gz"
  )
)
```

Taxonomy is optional at this stage. The ASV sequences can be classified later
without rerunning filtering and denoising.

Each plate is independently filtered, used to learn an error model, and
denoised with pseudo-pooling. Exact sequence tables are then merged, followed
by consensus chimera removal. The workflow records copies of the manifest,
configuration, command, and R session information with the results.

## 5. Review read retention and ASV lengths

```r
run_script("analysis/summarize_dada2_results.R")
```

Review:

- `results/dada2_2016/qc/run_qc_summary.csv`
- `results/dada2_2016/qc/sample_qc_flags.csv`
- `results/dada2_2016/qc/asv_length_summary.csv`

Investigate samples with low final retention or no non-chimeric reads rather
than automatically deleting them. The old mothur analysis used a depth filter;
that threshold should not be silently carried into the new ASV analysis.

## 6. Use the ASV outputs downstream

Primary files:

- `asv_count_table.csv` — samples by ASV counts.
- `asv_sequences.csv` and `asv_sequences.fasta` — exact sequences and stable
  biological keys.
- `asv_taxonomy.csv` — present only when taxonomy was requested.
- `sample_metadata_and_read_tracking.csv` — sample manifest joined to read QC.
- `read_tracking.csv` — input, filtered, denoised, merged, and non-chimeric
  read counts.
- `sequence_table_asv.rds` — native DADA2 sample-by-sequence table.

Labels such as `ASV_000001` are assigned within this completed run. Retain the
exact nucleotide sequence when a stable biological identity is needed; do not
interpret the numeric label itself as a taxonomic or functional identifier.

## Troubleshooting and safe reruns

- The raw FASTQs are never modified.
- `sequences/`, `results/`, and taxonomy FASTA files are ignored by Git.
- Rerunning step 1 replaces only `data/derived/hja_2016_*` mapping tables.
- Rerunning step 3 replaces only quality PDFs and their sample index.
- Step 4 stops if `sequence_table_asv.rds` already exists. On the cluster, set
  `DADA2_OVERWRITE=1` to remove the incomplete default output and rerun the
  entire workflow, or use a new output directory for a deliberately separate
  analysis attempt.
- If filtering or merging is poor, change the configuration and use a new
  output directory so the attempts remain auditable.

## References

- DADA2 paired-end tutorial: https://benjjneb.github.io/dada2/tutorial
- DADA2 multi-run guidance: https://benjjneb.github.io/dada2/bigdata.html
- DADA2 taxonomy training files: https://benjjneb.github.io/dada2/training.html
- Public metadata for unresolved `hja2016_200`:
  https://www.ncbi.nlm.nih.gov/biosample/SAMN12129905
