#!/bin/bash
#SBATCH -J cb_AF_array
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH --array=1-7    
#SBATCH -o /mnt/archive/farhadie/public_studies/Bondareva_AF/log/cb_AF_%A_%a.out  
#SBATCH -e /mnt/archive/farhadie/public_studies/Bondareva_AF/log/cb_AF_%A_%a.err

module load cuda/11.6

# Define sample array
SAMPLES=(AF1 AF2 AF3 AF4 AF5 AF6 AF7)

# Get current sample based on array task ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: $SAMPLE"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

SIF="/mnt/archive/farhadie/container/cellbender_0.3.0.sif"
INPUT="/mnt/archive/farhadie/public_studies/Bondareva_AF/cellranger_runs/${SAMPLE}/${SAMPLE}/outs/raw_feature_bc_matrix.h5"
OUTPUT="/mnt/archive/farhadie/public_studies/Bondareva_AF/cellranger_runs/${SAMPLE}/${SAMPLE}/outs/${SAMPLE}_cellbender.h5"

singularity exec --nv \
  --bind /mnt/archive:/mnt/archive \
  $SIF cellbender remove-background \
    --input $INPUT \
    --output $OUTPUT \
    --expected-cells 7500 \
    --total-droplets-included 30000 \
    --fpr 0.01 \
    --epochs 150 \
    --cuda

echo "Output files for $SAMPLE:"
ls -lh /mnt/archive/farhadie/public_studies/Bondareva_AF/cellranger_runs/${SAMPLE}/${SAMPLE}/outs/${SAMPLE}_cellbender*
