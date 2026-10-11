# Sediment-only iCAMP deployment

Prepared October 10, 2026. The approved read-only preflight completed on
October 10. All four remote input hashes match locally, `logs/` exists,
`squeue -u niw7` lists no jobs, and no sediment output paths were found.
Remote Git status includes untracked catchment scripts, backups, and results;
preserve them. The separately approved two-script transfer completed
successfully on October 10. The separately approved submission was accepted
as job **432228**. Queue/run state has not been checked after submission.

## Reviewable run

- Full R analysis: [`../microbes/06_icamp_sediment.R`](../microbes/06_icamp_sediment.R).
- Full batch script: [`run_icamp_sediment.sh`](run_icamp_sediment.sh).
- Host: `niw7@lugh`.
- Checkout: `/mnt/home/niw7/scratch/GitHub/hja-synoptic`.
- Request: `normal`, one node, one task, 16 CPUs, 192 GB, 96 hours.
- Expected runtime: unknown. The larger catchment run took 66 hours 9 minutes
  on 32 CPUs; it is a reference, not a sediment-runtime prediction.
- Inputs verified locally: all 34 sediment samples at exactly 10,000 reads,
  17,946 detected ASVs; the full master tree is pruned/rooted in memory.
- Settings retained: 1,000 randomizations; seed 2016; bMPD; Bray; Confidence;
  within-bin phylogenetic and across-all taxonomic randomization; bin minimum
  24; `ds = 0.2`. No habitat or environmental completeness filtering is applied.
- Outputs: `results/icamp_2016/sediment_castor/`; 561 expected sample pairs.
- Local changes: explicit `normal` partition; `temp_saved_file <- NULL` exposed
  and passed to iCAMP, matching the completed catchment workflow. Resume only
  from the same sediment run with unchanged inputs/settings.
- Distance cache: build a sediment-specific matrix and use it only after hash,
  tip-order, completion-record, and numerical checks. Preserve incomplete
  matrices; never substitute the catchment descriptor.

Local R parsing and batch-shell syntax checks passed. The local R installation
does not have iCAMP, so no new local end-to-end iCAMP run was performed. The
castor workflow has earlier local diagnostics and a completed cluster catchment
run, documented in [`../README_ICAMP_RERUN.md`](../README_ICAMP_RERUN.md).

## Completed read-only preflight

The approved command below checked remote changes, queue entries, input
fingerprints, existing sediment outputs, and Slurm's log directory. It
modified no project files and submitted no job.

```sh
ssh niw7@lugh 'cd /mnt/home/niw7/scratch/GitHub/hja-synoptic && git status --short && squeue -u niw7 && ls -ld logs && md5sum data/derived/microbial_diversity_2016/asv_counts_10k.csv data/derived/microbial_diversity_2016/asv_metadata_screened.csv data/derived/microbial_diversity_2016/sample_metadata.csv results/phylogeny_2016/asv_tree_screened_fasttree.nwk && find results/icamp_2016 -maxdepth 2 -iname "*sediment*" -print'
```

Locally verified fingerprints for comparison:

```text
asv_counts_10k.csv: d9121c9599889ba92cc330e258a4c569
asv_metadata_screened.csv: 7fc8eb797b749ee90fd4bab5816f68f6
sample_metadata.csv: 0aa9d2533bf5910bd965b64e4a11c928
asv_tree_screened_fasttree.nwk: e29b5c119cf476c07529eab3c6c50aa4
```

## Completed two-script transfer

The following transfer to `niw7@lugh` completed with exit code 0 after explicit
approval. Rsync reported both requested scripts transferred. It preserved
relative directories and enabled backups for existing destination scripts.
No deletion, commit, push, or submission was included.

```sh
rsync -av --backup --suffix=.before-sediment-deploy-20261010 --relative analysis/microbes/06_icamp_sediment.R analysis/cluster/run_icamp_sediment.sh niw7@lugh:/mnt/home/niw7/scratch/GitHub/hja-synoptic/
```

Obtain separate approval for a read-only cluster package check if needed.
Do not overwrite active scripts or synchronize unrelated local/remote changes.
No commit or push is needed to transfer these scripts. Input mismatches must
be resolved explicitly before submission.

## Accepted submission and next checks

The following command was executed after separate approval and returned
`Submitted batch job 432228` (exit code 0):

```sh
ssh niw7@lugh 'cd /mnt/home/niw7/scratch/GitHub/hja-synoptic && sbatch analysis/cluster/run_icamp_sediment.sh'
```

Logs are `logs/hja2016_icamp_sediment-432228.out` and `.err` in the remote
checkout. Outputs belong in `results/icamp_2016/sediment_castor/`. Submission
acceptance does not establish that the job is running or complete. Do not
modify its deployed scripts. Later status/log checks and result transfers
each require their own authorization. Review completion, exports, input
alignment, distance validation, binning, and phylogenetic signal before
interpreting results.
