while read -r SAMPLE; do
  [ -z "$SAMPLE" ] && continue
  [ "${SAMPLE:0:1}" = "#" ] && continue

  SAMPLE_ID="$(basename "$SAMPLE")"
  COUNTS="$SAMPLE/outs/filtered_feature_bc_matrix.h5"
  CSV="$SAMPLE/velocyto/fraction_unspliced.csv"
  OUT="$SAMPLE/outs/${SAMPLE_ID}_QClus.h5ad"

  if [ ! -f "$COUNTS" ] || [ ! -f "$CSV" ]; then
    echo "[SKIP] $SAMPLE: missing counts or csv"
    continue
  fi

  # Make variables visible to the Python child process
  export SAMPLE COUNTS CSV OUT

  python - <<'PYCODE'
# Read environment variables exported from Bash
import os
import pandas as pd
from qclus import qclus as qc

sample   = os.environ["SAMPLE"]
counts   = os.environ["COUNTS"]
csv_path = os.environ["CSV"]
out_path = os.environ["OUT"]

# 1) read as Series (index=barcodes, values=fraction)
fu = pd.read_csv(csv_path, index_col=0)["fraction_unspliced"]
fu = pd.to_numeric(fu, errors="coerce").dropna()

# 2) make both index variants: stripped and with "-1"
base_idx = fu.index.astype(str).str.strip().str.replace(r"-\d+$", "", regex=True)
fu_both = pd.concat([
    pd.Series(fu.values, index=base_idx + "-1"),
    pd.Series(fu.values, index=base_idx),
])
fu_both = fu_both[~fu_both.index.isin(["", "nan", "None"])]

# 3) run QClus and write output
adata = qc.run_qclus(counts, fu_both)
adata.write_h5ad(out_path, compression="gzip")
print(f"Wrote {out_path}")
PYCODE
done < /mnt/archive/farhadie/public_studies/Leblanc/cellranger_runs/samples.txt
