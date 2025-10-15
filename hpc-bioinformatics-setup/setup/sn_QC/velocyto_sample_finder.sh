# Build samples.txt for velocyto without deep find
# (English comments in code)
set -euo pipefail

PROJ=hill_sc
BASE=$LAB_ARCHIVE/public_studies/Hill
RUNS=$BASE/cellranger_runs             # contains <SAMPLE>/
SAMPLES_TXT="$RUNS/samples.txt"
TMP="$SAMPLES_TXT.tmp"

: > "$TMP"                               # truncate temp file

shopt -s nullglob                        # avoid literal globs when empty
for S in "$RUNS"/*; do
  [ -d "$S" ] || continue
  BN=$(basename "$S")
  PARENT="$S/$BN"                        # expected parent of 'outs'
  OUTS="$PARENT/outs"
  BAM="$OUTS/possorted_genome_bam.bam"
  BAI="$OUTS/possorted_genome_bam.bam.bai"

  # Fast checks: no recursion, no find
  if [ -d "$OUTS" ] && [ -s "$BAM" ] && [ -s "$BAI" ]; then
    echo "$PARENT" >> "$TMP"
  else
    echo "[SKIP] $BN (missing outs/BAM/BAI)" >&2
  fi
done
shopt -u nullglob

mv -f "$TMP" "$SAMPLES_TXT"
echo "[INFO] samples.txt -> $SAMPLES_TXT"
wc -l "$SAMPLES_TXT"
nl -ba "$SAMPLES_TXT" | head
