#!/bin/bash
#SBATCH -J vly_run10x_raw
#SBATCH -p hpc
#SBATCH -t 12:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --array=1-34
#SBATCH -o $LAB_ARCHIVE/public_studies/Hill/cellranger_runs/_logs/velocyto_raw_%A_%a.out
#SBATCH -e $LAB_ARCHIVE/public_studies/Hill/cellranger_runs/_logs/velocyto_raw_%A_%a.err

set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"

BASE=$LAB_ARCHIVE/public_studies/Hill
RUNS=$BASE/cellranger_runs
SAMPLES_FILE=$RUNS/samples.txt
GTF="/mnt/archive/farhadie/ref/hg38/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf"

SAMPLE_DIR="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")"
SAMPLE_ID="$(basename "$SAMPLE_DIR")"

BAM="$SAMPLE_DIR/outs/possorted_genome_bam.bam"
BAI="$SAMPLE_DIR/outs/possorted_genome_bam.bam.bai"
RAW_BCS="$SAMPLE_DIR/outs/raw_feature_bc_matrix/barcodes.tsv.gz"
OUTDIR="$SAMPLE_DIR/velocyto_raw"

[[ -s "$BAM" ]] || { echo "[ERROR] Missing BAM: $BAM"; exit 1; }
[[ -s "$BAI" ]] || { echo "[ERROR] Missing BAI: $BAI"; exit 1; }
[[ -s "$RAW_BCS" ]] || { echo "[ERROR] Missing RAW barcodes: $RAW_BCS"; exit 1; }
[[ -s "$GTF" || -s "${GTF}.gz" ]] || { echo "[ERROR] Missing GTF: $GTF"; exit 1; }

TMP_GTF=""
if [[ ! -s "$GTF" && -s "${GTF}.gz" ]]; then
  TMP_GTF="$SAMPLE_DIR/genes.tmp.gtf"
  gunzip -c "${GTF}.gz" > "$TMP_GTF"
  GTF="$TMP_GTF"
fi

mkdir -p "$OUTDIR"

velocyto run -b "$RAW_BCS" -o "$OUTDIR" "$BAM" "$GTF"

LOOM=$(ls -1 "$OUTDIR"/*.loom | head -n 1) || true
[[ -s "$LOOM" ]] || { echo "[ERROR] Loom not found in $OUTDIR"; exit 1; }

python - <<'PY'
import scvelo as scv, pandas as pd, numpy as np, os
loom_path = os.environ.get("LOOM")
ldata = scv.read_loom(loom_path, sparse=True)
u = ldata.layers['unspliced'].sum(1).A1
s = ldata.layers['spliced'].sum(1).A1
frac = u / (u + s + 1e-9)
df = pd.DataFrame({'barcode': ldata.obs_names, 'fraction_unspliced': frac})
out_csv = os.path.join(os.path.dirname(loom_path), "fraction_unspliced.csv")
df.to_csv(out_csv, index=False)
print(f"[OK] Wrote {out_csv} with {df.shape[0]} rows")
PY
export LOOM="$LOOM"

[[ -n "$TMP_GTF" ]] && rm -f "$TMP_GTF"
