#!/bin/bash
#SBATCH --job-name=human_Cellranger
#SBATCH --partition=pe2
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=120:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human/logs/%x_%j.err

CELLRANGER_PATH="/gpfs/commons/home/acorman/software/cellranger-8.0.1/cellranger"
TRANSCRIPTOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A"
OUTPUT_BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger/human"
CHEMISTRY_FILE="/gpfs/commons/home/jameslee/Cas13/liver_atlas/chemistry_mapping.csv"
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/fastq"

mkdir -p ${OUTPUT_BASE_DIR}/logs

# Read the mapping file and process each sample
while IFS=, read -r sample chemistry; do
    # Skip the header line
    if [[ "$sample" == "sample_id" ]]; then
        continue
    fi

    # Define output directory and ID
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/cellranger_${sample}"
    ID="${sample}"

    echo "Processing sample ${sample} with chemistry ${chemistry}..."

    # Run Cell Ranger with specified chemistry
    ${CELLRANGER_PATH} count \
        --id=${ID} \
        --description="Output for ${sample}" \
        --transcriptome=${TRANSCRIPTOME} \
        --libraries=${LIBRARY_CSV} \
        --chemistry=${chemistry} \
        --create-bam=true \
        --output-dir=${OUTPUT_DIR}
done < ${CHEMISTRY_FILE}