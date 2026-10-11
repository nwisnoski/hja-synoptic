# Cluster logs

Slurm writes standard output and error here for the batch scripts in
`analysis/cluster/`, including DADA2, phylogeny, and iCAMP. Generated `.out`
and `.err` files are ignored by Git. Batch jobs require this directory to exist
before submission because Slurm opens the log files before the script runs.

Sediment job 432228 uses `hja2016_icamp_sediment-432228.out` and `.err` on the
cluster. Its current state has not been checked. See the
[deployment record](../analysis/cluster/README_SEDIMENT_DEPLOYMENT.md).
Reading remote logs, transferring them, or checking job state requires approval.
