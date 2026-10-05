# iCAMP rerun plan

Updated on 2026-10-05. Both R scripts now calculate distances with castor.
Local validation is described below. The catchment pilot R and batch scripts
were transferred to the cluster with approval on October 5. Pilot job
**432023** was submitted with approval that day. Its current queue/run state
has not been checked. No final 1,000-draw or sediment job has been submitted.

## 1. What failed

Catchment job 431635 reached its 96-hour limit on 2026-10-04. It was still
inside iCAMP's phylogenetic-distance construction for 49,260 ASVs. No bin-level
null randomizations or final process results were produced. Peak resident
memory was approximately 18.4 GiB, well below the 256 GB allocation.

The saved `results/icamp_2016/catchment/phylogenetic_distance/pd.bin` has a
logical size of about 19 GB, but this size does not establish completion. A
12-tip check found 49 of 66 off-diagonal entries still zero; the same pairs
all have positive distances in the input tree. Treat the matrix as incomplete.
iCAMP can reuse an existing descriptor automatically, so resubmitting the
unchanged script is unsafe. The saved tree paths alone are not a completed
distance checkpoint.

## 2. Implemented distance calculation

Keep the existing screened phylogeny, fixed 10K community table, all 100
samples, and all 49,260 detected ASVs. Improve distance computation without
changing the ecological dataset or iCAMP settings.

Both scripts use `castor::get_all_pairwise_distances()` with
`only_clades = seq_along(tree$tip.label)` and `as_edge_counts = FALSE`.
This calculates branch-length distances between tips and avoids requesting
the much larger matrix that includes internal nodes. Castor is designed for
large phylogenies ([Louca & Doebeli 2018](https://doi.org/10.1093/bioinformatics/btx701));
its [function documentation](https://search.r-project.org/CRAN/refmans/castor/html/get_all_pairwise_distances.html)
specifies the tip order and branch-length metric. Castor 1.8.7 is installed
locally and passed its load check in the cluster R 4.4 user library on
2026-10-05. The approved remote check also found here 1.0.2 and iCAMP 1.9.1.

Local benchmarks on spread-out subsets of the actual pruned, midpoint-rooted
catchment tree took 0.002, 0.018, and 0.207 seconds for 500, 2,000, and 5,000
tips, respectively. These timings cover distance calculation only, not input
loading or disk writing. Each calculation matched `ape::cophenetic.phylo()`
on 64 spread-out tips; the largest absolute discrepancy was below 5e-15.
Distances for these same tips in the original tree also matched those after
pruning. These are small-tree checks, not a measured full-size speedup.

Both scripts also passed an end-to-end local diagnostic using eight real
samples, 600 abundant ASVs, ten null draws, and two workers (iCAMP 1.8.6).
This exercised fresh distance construction, safe cache reuse, rejection of
an incomplete cache, and all CSV exports. All 28 sample pairs had the expected
five process columns and fractions summing to one. Diagnostic outputs stayed
outside the repository and are not ecological results. The cluster's iCAMP
1.9.1 must still be validated by the full-community pilot.

The complete tip matrix itself needs approximately 18.1 GiB. Additional working
memory and a disk-backed copy are required. Write the completed matrix to
iCAMP's `bigmemory` format in blocks, preserving the exact tip-name order,
then pass `pd.desc`, `pd.spname`, and `pd.wd` explicitly to `icamp.big()`.
The [iCAMP interface](https://search.r-project.org/CRAN/refmans/iCAMP/html/icamp.big.html)
accepts a precomputed distance matrix.

The distance stage writes columns in 256-tip blocks, then checks the entire
disk-backed matrix for finite/nonnegative entries, symmetry, and a zero
diagonal. It independently checks 64 spread-out tips against ape. A
`distance_complete.csv` input-checksum record is written last. On reuse, the
script requires matching input hashes and exact tip order, and repeats the
matrix checks. A descriptor without that record is rejected.

New output directories are `results/icamp_2016/catchment_castor/` and
`results/icamp_2016/sediment_castor/`. The failed `catchment/` directory stays
in place, untouched. No cleanup or folder move is needed for the rerun.

## 3. Cluster execution order

1. Package/version check and castor installation completed with approval on
   2026-10-05. Here 1.0.2, iCAMP 1.9.1, and castor 1.8.7 are available.
2. Catchment pilot R and batch scripts were transferred with approval on
   2026-10-05. Rsync completed successfully and retained previous script copies
   with suffix `.before-castor-pilot-20261005-8zLKhk`. No data or result files
   were transferred or removed. Sediment scripts have not been transferred.
3. For the catchment pilot, set `test_run <- TRUE` near the top of the existing
   R script before transferring that version. That pilot copy is now on the
   cluster; the local repository version remains `FALSE`. Its batch copy declares
   partition `normal`, a 24-hour limit, a distinct pilot job name, and a pilot
   environment-record path. Both copies replaced the corresponding cluster
   script names, with the old scripts backed up. Pilot job 432023 was submitted
   with approval; the submission command below is a record, not a request to
   submit a duplicate:

   ```sh
   cd /mnt/home/niw7/scratch/GitHub/hja-synoptic
   mkdir -p logs
   sbatch --partition=normal --time=24:00:00 --job-name=hja2016_icamp_pilot \
     analysis/cluster/run_icamp_catchment.sh
   ```

   This is one node, one task, 32 CPUs, and 256 GB. Runtime is unknown; 24 hours
   is the pilot limit, not an estimate. The pilot keeps all 100 communities
   and 49,260 ASVs but uses 100 null draws. Outputs go under
   `catchment_castor/pilot/`; distances go under the main directory so the
   final run can reuse the validated matrix. The committed/default setting
   is `test_run <- FALSE`. Do not interpret pilot process fractions as final.
4. Inspect pilot completion, elapsed time, memory, binning, and distance checks.
   Select final resources from these measurements, not from an extrapolation
   of the small local benchmarks. Obtain approval before each remote action.
5. Set `test_run <- FALSE`, transfer the approved final version, and run
   the full catchment analysis with the established 1,000 randomizations,
   bMPD/Bray/Confidence settings and habitat-pair summaries. Then run the
   sediment-only analysis with the validated distance method.

The existing catchment request is one node, 32 CPUs, 256 GB, and 96 hours.
The actual runtime of the revised method is unknown. Confirm the partition and
allocation against the pilot before proposing a replacement submission; simply
extending the original slow calculation is not the proposed first step.
Every remote check, transfer, installation, pilot, and submission needs
its own explicit approval. Only pilot job 432023 has been submitted as part
of this plan; the final catchment and sediment runs remain pending.

## 4. Integration after completion

Sediment iCAMP still uses all 34 retained sediment communities. Subset its
pairwise outputs to the verified 33-site molecular intersection for FT-ICR
integration; exclude site 43 from that paired comparison only. Molecular
ecological analyses use presence/absence; intensities are reserved for QC.
Use the completed network analyses and `analysis/KEY_FINDINGS_2016.md` to select
focused comparisons rather than repeating their existing permutation grids.
