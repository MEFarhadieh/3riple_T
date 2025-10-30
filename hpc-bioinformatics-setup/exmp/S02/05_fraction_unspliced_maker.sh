# Make fraction_unspliced.csv file for QClus
BASE=$LAB_ARCHIVE/public_studies/Leblanc
RUNS=$BASE/cellranger_runs
SAMPLES_FILE=$RUNS/samples.txt

while read -r SAMPLE_DIR; do
  [ -z "$SAMPLE_DIR" ] && continue
  SAMPLE_ID=$(basename "$SAMPLE_DIR")
  
  # Try to find loom file in velocyto subdirectory first
  LOOM=("$SAMPLE_DIR/velocyto"/*.loom)
  [ ! -f "${LOOM[0]}" ] && LOOM=("$SAMPLE_DIR"/*.loom)
  [ ! -f "${LOOM[0]}" ] && { echo "[SKIP] $SAMPLE_ID: no loom file"; continue; }
  
  # Use the first matching file
  LOOM="${LOOM[0]}"
  
  CSV="$SAMPLE_DIR/velocyto/fraction_unspliced.csv"
  echo "[INFO] $SAMPLE_ID → $(basename "$CSV")"

  export LOOM CSV
  python - <<'PY'
import os
import qclus as qc
loom = os.environ["LOOM"]
csv_out = os.environ["CSV"]
fu = qc.utils.fraction_unspliced_from_loom(loom)
os.makedirs(os.path.dirname(csv_out), exist_ok=True)
fu.to_csv(csv_out)
print(f"[DONE] {csv_out}")
PY
done < "$SAMPLES_FILE"
