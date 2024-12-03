#!/bin/bash
#SBATCH --job-name=mouse_Cellranger
#SBATCH --partition=pe2
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=220G 
#SBATCH --time=120:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.err

module purge

# Define the path to Cell Ranger
CELLRANGER_PATH="/gpfs/commons/home/acorman/software/cellranger-8.0.1/cellranger"

# Define the path to the transcriptome
TRANSCRIPTOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-mm10-2024-A"

# Output directory where you want to store results
OUTPUT_BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/mouse"

# Path to the FASTQ files directory
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq"

# Path to the library CSV file
LIBRARY_CSV="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/library_mouse_1.csv"

# Create logs directory if it doesn't exist
mkdir -p ${OUTPUT_BASE_DIR}/logs

# Get the sample list for Mus musculus
SAMPLE_CSV="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/human_mouse_SRR_list.csv"
samples=($(awk -F, '$2 == "Mus musculus" {print $1}' $SAMPLE_CSV | tail -n +2))

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
      --fastqs=${FASTQ_DIR} \
      --libraries=${LIBRARY_CSV} \
      --create-bam=true \
      --output-dir=${OUTPUT_DIR}
done


