#!/bin/bash -l

# Run the reviewed 2016 regional-pool iCAMP analysis on the Slurm cluster.
# Submit from the repository root with:
#   sbatch analysis/cluster/run_icamp_reviewed_pool.sh

# Runtime is unknown. Sediment-only benchmark: 6:16:45 on 16 CPUs.
# The pooled run has more ASVs and pairs; this is a walltime limit, not an estimate.
#SBATCH --time=7-00:00:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192GB
#SBATCH --job-name=hja2016_icamp_reviewed_pool
#SBATCH --chdir=/mnt/home/niw7/scratch/GitHub/hja-synoptic
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

HJA_PROJECT_ROOT="/mnt/home/niw7/scratch/GitHub/hja-synoptic"
HJA_OUTPUT_DIR="${HJA_PROJECT_ROOT}/results/icamp_2016/reviewed_pool_castor"

cd "${HJA_PROJECT_ROOT}"
mkdir -p logs "${HJA_OUTPUT_DIR}"

required_inputs=(
  "analysis/supporting/13_reviewed_pool_icamp.R"
  "data/derived/soil_porewater_sediment_2016/asv_counts_10k.csv"
  "data/derived/microbial_diversity_2016/asv_metadata_screened.csv"
  "data/derived/soil_porewater_sediment_2016/sample_metadata.csv"
  "results/phylogeny_2016/asv_tree_screened_fasttree.nwk"
  "data/StreamCenterline.csv"
)
for required_input in "${required_inputs[@]}"; do
  if [[ ! -f "${required_input}" ]]; then
    echo "Missing required input: ${HJA_PROJECT_ROOT}/${required_input}"
    exit 3
  fi
done

if [[ -f "${HJA_OUTPUT_DIR}/icamp_result.rds" ]]; then
  echo "The completed iCAMP result already exists: ${HJA_OUTPUT_DIR}/icamp_result.rds"
  exit 2
fi

module load R/4.4.0

Rscript --vanilla -e '
required <- c("ape", "phangorn", "iCAMP", "bigmemory", "vegan", "castor", "here",
  "readr", "dplyr", "tidyr", "ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
stopifnot(as.character(packageVersion("iCAMP")) == "1.9.1")
'

{
  echo "job_id=${SLURM_JOB_ID}"
  echo "host=$(hostname)"
  echo "start_time=$(date --iso-8601=seconds)"
  echo "allocated_cpus=${SLURM_CPUS_PER_TASK}"
  echo "allocated_memory=192GB"
  module list
  R --version
} > "logs/hja2016_icamp_reviewed_pool-${SLURM_JOB_ID}.environment.txt" 2>&1

echo "[$(date --iso-8601=seconds)] Starting reviewed-pool iCAMP analysis."
Rscript --vanilla analysis/supporting/13_reviewed_pool_icamp.R
echo "end_time=$(date --iso-8601=seconds)" \
  >> "logs/hja2016_icamp_reviewed_pool-${SLURM_JOB_ID}.environment.txt"
echo "[$(date --iso-8601=seconds)] Reviewed-pool iCAMP analysis complete."
