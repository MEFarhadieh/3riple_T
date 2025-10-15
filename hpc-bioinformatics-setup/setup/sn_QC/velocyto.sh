# Build a sample list from directories matching (in my case all samples id start with AF)
BASE="/path/to/cellranger_run"
mkdir -p "$BASE"/_logs

find "$BASE" -mindepth 1 -maxdepth 1 -type d -name "AF*" \
 | sort > "$BASE/samples.txt"

wc -l "$BASE/samples.txt"
cat "$BASE/samples.txt"

# Write SLURM Array job for velocyto + fraction_unspliced extraction
# Write SLURM Array job: run velocyto + extract fraction_unspliced
# Write a minimal SLURM array job for velocyto run10x
cat > "$BASE/run_velocyto_run10x.sbatch" <<'SBATCH'
#!/bin/bash
#SBATCH -J vly_run10x
#SBATCH -p hpc
#SBATCH -t 08:00:00
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --array=1-__NSAMPLES__%3
#SBATCH -o __BASE__/_logs/velocyto_%A_%a.out
#SBATCH -e __BASE__/_logs/velocyto_%A_%a.err

set -euo pipefail

# Ensure velocyto is found (user pip installs)
export PATH="$HOME/.local/bin:$PATH"

BASE_DIR="__BASE__"
SAMPLES_FILE="$BASE_DIR/samples.txt"
GTF="/mnt/archive/farhadie/ref/hg38/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf"

SAMPLE_DIR="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")"
SAMPLE_ID="$(basename "$SAMPLE_DIR")"
echo "[INFO] Sample: $SAMPLE_ID"

# Expect these to exist already (Cell Ranger outputs):
BAM="$SAMPLE_DIR/possorted_genome_bam.bam"
BAI="$SAMPLE_DIR/possorted_genome_bam.bam.bai"
BCS_DIR="$SAMPLE_DIR/filtered_feature_bc_matrix"

[[ -s "$BAM" ]] || { echo "[ERROR] Missing BAM: $BAM"; exit 1; }
[[ -s "$BAI" ]] || { echo "[ERROR] Missing BAI: $BAI"; exit 1; }
[[ -d "$BCS_DIR" ]] || { echo "[ERROR] Missing matrix dir: $BCS_DIR"; exit 1; }
[[ -s "$GTF" || -s "${GTF}.gz" ]] || { echo "[ERROR] Missing GTF"; exit 1; }

# Create a minimal 10x-like layout for run10x (outs/), if not present
OUTS="$SAMPLE_DIR/outs"
mkdir -p "$OUTS"
# Symlink only if targets not already under outs/
[[ -e "$OUTS/possorted_genome_bam.bam" ]] || ln -sf ../possorted_genome_bam.bam "$OUTS/possorted_genome_bam.bam"
[[ -e "$OUTS/possorted_genome_bam.bam.bai" ]] || ln -sf ../possorted_genome_bam.bam.bai "$OUTS/possorted_genome_bam.bam.bai"
[[ -e "$OUTS/filtered_feature_bc_matrix" ]] || ln -sfn ../filtered_feature_bc_matrix "$OUTS/filtered_feature_bc_matrix"

# Uncompress GTF temporarily if only .gz exists
TMP_GTF=""
if [[ ! -s "$GTF" && -s "${GTF}.gz" ]]; then
  TMP_GTF="$SAMPLE_DIR/genes.tmp.gtf"
  gunzip -c "${GTF}.gz" > "$TMP_GTF"
  GTF="$TMP_GTF"
fi

# Run velocyto (pure and simple)
echo "[INFO] Running: velocyto run10x $SAMPLE_DIR $GTF"
velocyto run10x "$SAMPLE_DIR" "$GTF"

# Cleanup temp GTF if created
[[ -n "$TMP_GTF" ]] && rm -f "$TMP_GTF"

echo "[DONE] $SAMPLE_ID"
SBATCH


# Replace placeholders
NSAMPLES=$(wc -l < "$BASE/samples.txt")
sed -i "s|__NSAMPLES__|$NSAMPLES|g" "$BASE/run_velocyto_run10x.sbatch"
sed -i "s|__BASE__|$BASE|g"        "$BASE/run_velocyto_run10x.sbatch"

sbatch "$BASE/run_velocyto_run10x.sbatch"

# check job status
squeue -u $USER -o "%.18i %.9P %.10T %.20R %.8M %.2t %.10M %.30j"
