#!/bin/bash -l

# Build the full screened 2016 HJA ASV phylogeny on the Slurm cluster.
# Submit from the repository root with:
#   sbatch analysis/cluster/run_phylogeny.sh

#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --job-name=hja2016_tree
#SBATCH --chdir=/mnt/home/niw7/scratch/GitHub/hja-synoptic
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

HJA_PROJECT_ROOT="/mnt/home/niw7/scratch/GitHub/hja-synoptic"
FASTTREE_BIN="${HOME}/scratch/shared/bioinformatics/FastTree/FastTree"

cd "${HJA_PROJECT_ROOT}"
mkdir -p logs results/phylogeny_2016

if [[ ! -f "data/derived/microbial_diversity_2016/asv_metadata_screened.csv" ]]; then
  echo "Missing data/derived/microbial_diversity_2016/asv_metadata_screened.csv"
  echo "Run analysis/microbes/01_prepare_diversity.R before submitting this job."
  exit 3
fi
if [[ ! -x "${FASTTREE_BIN}" ]]; then
  echo "FastTree is missing or is not executable: ${FASTTREE_BIN}"
  exit 3
fi

module load openmpi4
module load R/4.4.0

export FASTTREE_BIN
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

Rscript --vanilla -e '
required <- c("Biostrings", "DECIPHER", "ape")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
'

{
  echo "job_id=${SLURM_JOB_ID}"
  echo "host=$(hostname)"
  echo "start_time=$(date --iso-8601=seconds)"
  echo "allocated_cpus=${SLURM_CPUS_PER_TASK}"
  echo "fasttree=${FASTTREE_BIN}"
  module list
  R --version
} > results/phylogeny_2016/cluster_environment.txt 2>&1

echo "[$(date --iso-8601=seconds)] Starting full screened-ASV phylogeny."
Rscript --vanilla analysis/microbes/03_build_phylogeny.R
echo "end_time=$(date --iso-8601=seconds)" \
  >> results/phylogeny_2016/cluster_environment.txt
echo "[$(date --iso-8601=seconds)] Phylogeny workflow complete."
