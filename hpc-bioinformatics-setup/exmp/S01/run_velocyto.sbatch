#!/bin/bash
#SBATCH -J vly_run10x
#SBATCH -p hpc
#SBATCH -t 12:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --array=1-34
#SBATCH -o $LAB_ARCHIVE/public_studies/Hill/cellranger_runs/_logs/velocyto_%A_%a.out
#SBATCH -e $LAB_ARCHIVE/public_studies/Hill/cellranger_runs/_logs/velocyto_%A_%a.err

set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"   # ensure velocyto is visible

BASE=$LAB_ARCHIVE/public_studies/Hill
RUNS=$BASE/cellranger_runs
SAMPLES_FILE=$RUNS/samples.txt
GTF="/mnt/archive/farhadie/ref/hg38/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf"

SAMPLE_DIR="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")"  # parent of 'outs'
SAMPLE_ID="$(basename "$SAMPLE_DIR")"

BAM="$SAMPLE_DIR/outs/possorted_genome_bam.bam"
BAI="$SAMPLE_DIR/outs/possorted_genome_bam.bam.bai"
BCS_DIR="$SAMPLE_DIR/outs/filtered_feature_bc_matrix"

[[ -s "$BAM" ]] || { echo "[ERROR] Missing BAM: $BAM"; exit 1; }
[[ -s "$BAI" ]] || { echo "[ERROR] Missing BAI: $BAI"; exit 1; }
[[ -d "$BCS_DIR" ]] || { echo "[ERROR] Missing matrix dir: $BCS_DIR"; exit 1; }
[[ -s "$GTF" || -s "${GTF}.gz" ]] || { echo "[ERROR] Missing GTF: $GTF"; exit 1; }

TMP_GTF=""
if [[ ! -s "$GTF" && -s "${GTF}.gz" ]]; then
  TMP_GTF="$SAMPLE_DIR/genes.tmp.gtf"
  gunzip -c "${GTF}.gz" > "$TMP_GTF"
  GTF="$TMP_GTF"
fi

echo "[INFO] Running velocyto on $SAMPLE_ID"
velocyto run10x "$SAMPLE_DIR" "$GTF"

[[ -n "$TMP_GTF" ]] && rm -f "$TMP_GTF"
