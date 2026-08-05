#!/bin/bash -l

# Temporarily continue the 2016 HJA DADA2 run from the completed per-run
# pre-chimera sequence tables. Restore the fresh-output guard after the final
# outputs are generated and the full R processing block is re-enabled.
#
# Submit from the repository root with:
#   sbatch analysis/cluster/run_dada2.sh

#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --job-name=hja2016_dada2
#SBATCH --chdir=/mnt/home/niw7/scratch/GitHub/hja-synoptic
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

HJA_PROJECT_ROOT="/mnt/home/niw7/scratch/GitHub/hja-synoptic"
HJA_OUTPUT_SETTING="${DADA2_OUTPUT_DIR:-results/dada2_2016}"
HJA_GENUS_TRAINING="reference/silva_nr99_v138.2_toGenus_trainset.fa.gz"
HJA_SPECIES_TRAINING="reference/silva_v138.2_assignSpecies.fa.gz"

cd "${HJA_PROJECT_ROOT}"

if [[ "${HJA_OUTPUT_SETTING}" = /* ]]; then
  HJA_OUTPUT_PATH="${HJA_OUTPUT_SETTING}"
else
  HJA_OUTPUT_PATH="${HJA_PROJECT_ROOT}/${HJA_OUTPUT_SETTING}"
fi

if [[ ! -d "${HJA_OUTPUT_PATH}" ]]; then
  echo "Missing the existing DADA2 output directory containing checkpoints:"
  echo "  ${HJA_OUTPUT_PATH}"
  exit 2
fi
for run_id in HJA2016_Plate1 HJA2016_Plate2; do
  for checkpoint in \
    sequence_table_prechimera.rds \
    read_tracking_prechimera.csv; do
    checkpoint_path="${HJA_OUTPUT_PATH}/runs/${run_id}/${checkpoint}"
    if [[ ! -f "${checkpoint_path}" ]]; then
      echo "Missing required checkpoint: ${checkpoint_path}"
      exit 2
    fi
  done
done

required_inputs=(
  "analysis/build_sample_manifest.R"
  "analysis/check_dada2_setup.R"
  "analysis/dada2_pipeline.R"
  "analysis/summarize_dada2_results.R"
  "config/dada2_run_config.csv"
  "data/hja-synoptic_sequence-sample-list.csv"
  "data/hja-env_data_clean.csv"
  "data/hja-synoptic_env-data-soils.csv"
)
for required_input in "${required_inputs[@]}"; do
  if [[ ! -f "${required_input}" ]]; then
    echo "Missing required input: ${HJA_PROJECT_ROOT}/${required_input}"
    exit 3
  fi
done
if [[ ! -d "sequences" ]]; then
  echo "Missing sequence directory: ${HJA_PROJECT_ROOT}/sequences"
  exit 3
fi

fastq_count="$(
  find sequences -maxdepth 1 -type f \
    -name 'GSF2238-Lennon-Plate*-hja2016_*_R[12]_001.fastq.gz' |
    wc -l |
    tr -d '[:space:]'
)"
if [[ "${fastq_count}" != "240" ]]; then
  echo "Expected 240 paired-read FASTQ files; found ${fastq_count}."
  exit 3
fi

module load R/4.4.0

HJA_ENVIRONMENT_FILE="${HJA_OUTPUT_PATH}/cluster_environment_continuation.txt"
{
  echo "job_id=${SLURM_JOB_ID:-not_available}"
  echo "job_name=${SLURM_JOB_NAME:-not_available}"
  echo "host=$(hostname)"
  echo "start_time=$(date --iso-8601=seconds)"
  echo "project_root=${HJA_PROJECT_ROOT}"
  echo "output_path=${HJA_OUTPUT_PATH}"
  echo "allocated_cpus=${SLURM_CPUS_PER_TASK:-not_available}"
  echo "fastq_count=${fastq_count}"
  module list
  R --version
} > "${HJA_ENVIRONMENT_FILE}" 2>&1

echo "[$(date --iso-8601=seconds)] Rebuilding the 2016 sample manifest."
Rscript --vanilla analysis/build_sample_manifest.R \
  --project-root "${HJA_PROJECT_ROOT}"

echo "[$(date --iso-8601=seconds)] Running the DADA2 preflight checks."
Rscript --vanilla analysis/check_dada2_setup.R \
  --project-root "${HJA_PROJECT_ROOT}" \
  --output "${HJA_OUTPUT_PATH}"

taxonomy_args=()
if [[ -f "${HJA_GENUS_TRAINING}" ]]; then
  taxonomy_args+=(--taxonomy "${HJA_GENUS_TRAINING}")
  if [[ -f "${HJA_SPECIES_TRAINING}" ]]; then
    taxonomy_args+=(--species "${HJA_SPECIES_TRAINING}")
  else
    echo "Species reference is absent; genus-level taxonomy will still run."
  fi
else
  echo "SILVA genus training file is absent; ASVs will be inferred without taxonomy."
  echo "Taxonomy can be assigned later from asv_sequences.fasta."
fi

echo "[$(date --iso-8601=seconds)] Starting DADA2."
Rscript --vanilla analysis/dada2_pipeline.R \
  --project-root "${HJA_PROJECT_ROOT}" \
  --output "${HJA_OUTPUT_PATH}" \
  "${taxonomy_args[@]}"

echo "[$(date --iso-8601=seconds)] Building DADA2 QC summaries."
Rscript --vanilla analysis/summarize_dada2_results.R \
  --project-root "${HJA_PROJECT_ROOT}" \
  --input "${HJA_OUTPUT_PATH}"

echo "end_time=$(date --iso-8601=seconds)" \
  >> "${HJA_ENVIRONMENT_FILE}"
echo "[$(date --iso-8601=seconds)] DADA2 cluster workflow complete."
