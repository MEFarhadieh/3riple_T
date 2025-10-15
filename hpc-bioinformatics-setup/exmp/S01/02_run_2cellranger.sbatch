#!/bin/bash
#SBATCH -J cr_count
#SBATCH -p hpc
#SBATCH -t 18:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --array=1-34
#SBATCH -o $LAB_SCRATCH/hill_sc/logs/cellranger_%A_%a.out
#SBATCH -e $LAB_SCRATCH/hill_sc/logs/cellranger_%A_%a.err

set -euo pipefail
set -x

# --- Paths from your env ---
PROJ=hill_sc
BASE=$LAB_ARCHIVE/public_studies/Hill
SCR=$LAB_SCRATCH/$PROJ
REN=$SCR/fastqs_renamed
QC=$SCR/qc
LOG=$SCR/logs

SAMPLES_FILE=$SCR/samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

# FASTQs per-sample (contains multiple lanes)
FASTQS="$SCR/fastqs_by_sample/$SAMPLE"

# Transcriptome reference (adjust if stored elsewhere)
TRANSCRIPTOME=/mnt/archive/farhadie/ref/hg38/cellranger/refdata-gex-GRCh38-2024-A

# Run directory (keep logs/intermediates) and final destination
RUNROOT=$BASE/cellranger_runs/$SAMPLE
DEST=$BASE/cellranger_out/$SAMPLE
mkdir -p "$RUNROOT" "$DEST"

# Optional: safer TMP on archive if scratch is small
export TMPDIR=$BASE/tmp
export TEMP=$BASE/tmp
export SINGULARITY_TMPDIR=$BASE/tmp
export SINGULARITY_CACHEDIR=$BASE/singularity_cache
mkdir -p "$TMPDIR" "$SINGULARITY_TMPDIR" "$SINGULARITY_CACHEDIR"

# snRNA toggle (set true for nuclei)
INTRONS=true
INTRON_FLAG=""
if [ "$INTRONS" = "true" ]; then
  INTRON_FLAG="--include-introns=true"
fi

# Move to run root so cellranger writes here
cd "$RUNROOT"

# Cell Ranger (BAM will be created by default; do NOT pass --no-bam)
cellranger count \
  --id="$SAMPLE" \
  --transcriptome="$TRANSCRIPTOME" \
  --fastqs="$FASTQS" \
  --sample="$SAMPLE" \
  --localcores=8 \
  --localmem=90 \
  --create-bam true \
  $INTRON_FLAG

# Collect outputs (matrices + web summary + metrics + BAM)
OUT="$SAMPLE/outs"
mkdir -p "$DEST"
