# Make fraction_unspliced.csv file for QClus
BASE=$LAB_ARCHIVE/public_studies/Hill
RUNS=$BASE/cellranger_runs
SAMPLES_FILE=$RUNS/samples.txt

while read -r SAMPLE_DIR; do
  [ -z "$SAMPLE_DIR" ] && continue
  SAMPLE_ID=$(basename "$SAMPLE_DIR")
  LOOM="$SAMPLE_DIR/velocyto/${SAMPLE_ID}.loom"
  [ ! -f "$LOOM" ] && LOOM="$SAMPLE_DIR/${SAMPLE_ID}.loom"
  [ ! -f "$LOOM" ] && { echo "[SKIP] $SAMPLE_ID: no loom file"; continue; }

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
