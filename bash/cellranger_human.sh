#!/bin/bash
#SBATCH --job-name=human_Cellranger
#SBATCH --partition=pe2
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=220G 
#SBATCH --time=120:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.err

# module purge

# Define the path to Cell Ranger
CELLRANGER_PATH="/gpfs/commons/home/acorman/software/cellranger-8.0.1/cellranger"

# Define the path to the transcriptome
TRANSCRIPTOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A"

# Output directory where you want to store results
OUTPUT_BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human"

# Path to the library CSV file
LIBRARY_CSV="/gpfs/commons/home/jameslee/Cas13/liver_atlas/library_1.csv"

# Create logs directory if it doesn't exist
mkdir -p ${OUTPUT_BASE_DIR}/logs

# Get the sample list for Homo sapiens
SAMPLE_CSV="/gpfs/commons/home/jameslee/Cas13/liver_atlas/human_SRR_list_1.txt"
samples=($(awk '{print $1}' $SAMPLE_CSV))

# Directory containing the fastq files
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq"

# List available samples in the fastq directory
available_samples=($(ls ${FASTQ_DIR}))

# Print available samples
echo "Available samples in ${FASTQ_DIR}:"
for available_sample in "${available_samples[@]}"; do
  echo "${available_sample}"
done

for sample in "${samples[@]}"; do
  # Define the output directory for this sample (shortened for --id)
  OUTPUT_DIR="${OUTPUT_BASE_DIR}/cellranger_${sample}"
  
  # Use a shorter id (only the sample name)
  ID="${sample}"

  # Run Cell Ranger for the sample using the library CSV
  ${CELLRANGER_PATH} count \
      --id=${ID} \
      --description="Output for ${sample}" \
      --transcriptome=${TRANSCRIPTOME} \
      --libraries=${LIBRARY_CSV} \
      --create-bam=true \
      --output-dir=${OUTPUT_DIR}
done