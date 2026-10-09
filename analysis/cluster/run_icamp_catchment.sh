#!/bin/bash -l

# Run the 2016 catchment iCAMP analysis on the Slurm cluster.
# Submit from the repository root with:
#   sbatch analysis/cluster/run_icamp_catchment.sh

#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB
#SBATCH --job-name=hja2016_icamp_catchment
#SBATCH --chdir=/mnt/home/niw7/scratch/GitHub/hja-synoptic
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

HJA_PROJECT_ROOT="/mnt/home/niw7/scratch/GitHub/hja-synoptic"
HJA_OUTPUT_DIR="${HJA_PROJECT_ROOT}/results/icamp_2016/catchment_castor"

cd "${HJA_PROJECT_ROOT}"
mkdir -p logs "${HJA_OUTPUT_DIR}"

required_inputs=(
  "analysis/microbes/07_icamp_catchment.R"
  "data/derived/microbial_diversity_2016/asv_counts_10k.csv"
  "data/derived/microbial_diversity_2016/asv_metadata_screened.csv"
  "data/derived/microbial_diversity_2016/sample_metadata.csv"
  "results/phylogeny_2016/asv_tree_screened_fasttree.nwk"
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
required <- c("ape", "phangorn", "iCAMP", "bigmemory", "vegan", "castor", "here")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
'

{
  echo "job_id=${SLURM_JOB_ID}"
  echo "host=$(hostname)"
  echo "start_time=$(date --iso-8601=seconds)"
  echo "allocated_cpus=${SLURM_CPUS_PER_TASK}"
  echo "allocated_memory=256GB"
  module list
  R --version
} > "${HJA_OUTPUT_DIR}/cluster_environment.txt" 2>&1

echo "[$(date --iso-8601=seconds)] Starting catchment iCAMP analysis."
Rscript --vanilla analysis/microbes/07_icamp_catchment.R
echo "end_time=$(date --iso-8601=seconds)" \
  >> "${HJA_OUTPUT_DIR}/cluster_environment.txt"
echo "[$(date --iso-8601=seconds)] Catchment iCAMP analysis complete."
