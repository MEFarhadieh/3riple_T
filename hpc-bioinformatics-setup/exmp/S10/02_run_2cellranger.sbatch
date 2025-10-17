#!/bin/bash
#SBATCH -J cr_count_AF
#SBATCH -p hpc
#SBATCH -t 18:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --array=1-8
#SBATCH -o $LAB_SCRATCH/Bondareva_AF_sc/logs/cellranger_%A_%a.out
#SBATCH -e $LAB_SCRATCH/Bondareva_AF_sc/logs/cellranger_%A_%a.err

set -euo pipefail
set -x

# --- Paths for Bondareva AF project ---
PROJ=Bondareva_AF_sc
BASE=$LAB_ARCHIVE/public_studies/Bondareva_AF
SCR=$LAB_SCRATCH/$PROJ
REN=$SCR/fastqs_renamed
QC=$SCR/qc
LOG=$SCR/logs

SAMPLES_FILE=$SCR/samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

# FASTQs per-sample (contains multiple lanes merged)
FASTQS="$SCR/fastqs_by_sample/$SAMPLE"

# Transcriptome reference
TRANSCRIPTOME=/mnt/archive/farhadie/ref/hg38/cellranger/refdata-gex-GRCh38-2024-A

# Run directory and final destination
RUNROOT=$BASE/cellranger_runs/$SAMPLE
DEST=$BASE/cellranger_out/$SAMPLE
mkdir -p "$RUNROOT" "$DEST"

# Temp directories
export TMPDIR=$BASE/tmp
export TEMP=$BASE/tmp
export SINGULARITY_TMPDIR=$BASE/tmp
export SINGULARITY_CACHEDIR=$BASE/singularity_cache
mkdir -p "$TMPDIR" "$SINGULARITY_TMPDIR" "$SINGULARITY_CACHEDIR"

# snRNA-seq or scRNA-seq toggle
# TODO: CHECK if this is nuclei (snRNA) or whole cells (scRNA)
# If nuclei -> INTRONS=true
# If whole cells -> INTRONS=false
INTRONS=true  # snRNA-seq
INTRON_FLAG=""
if [ "$INTRONS" = "true" ]; then
  INTRON_FLAG="--include-introns=true"
fi

# Move to run root
cd "$RUNROOT"

# Run CellRanger count
cellranger count \
  --id="$SAMPLE" \
  --transcriptome="$TRANSCRIPTOME" \
  --fastqs="$FASTQS" \
  --sample="$SAMPLE" \
  --localcores=8 \
  --localmem=90 \
  --create-bam=true \
  $INTRON_FLAG

# Copy outputs to final destination
OUT="$SAMPLE/outs"
if [ -d "$OUT" ]; then
  cp -r "$OUT"/* "$DEST/"
  echo "=== SUCCESS: Outputs copied to $DEST ==="
else
  echo "=== ERROR: Output directory not found: $OUT ===" >&2
  exit 1
fi
