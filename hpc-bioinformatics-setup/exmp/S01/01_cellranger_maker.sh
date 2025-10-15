############################################################################################
# Make suitable environment and prepare files and directories for cellranger counts
############################################################################################
# Define project variablea
PROJ=hill_sc
BASE=$LAB_ARCHIVE/public_studies/Hill
SCR=$LAB_SCRATCH/$PROJ
REN=$SCR/fastqs_renamed        
QC=$SCR/qc
LOG=$SCR/logs

# Make run(s) to sample files based on SRA run on NCBI
cat > "$SCR/run2sample.tsv" <<'EOF'
SRR27950916     1279_3n
SRR27950917     1279_3n
SRR27950915     1296_1n
SRR27950913     1334_3n
SRR27950914     1334_3n
SRR27950912     1357_2n
SRR27950910     1360_1n
SRR27950911     1360_1n
SRR27950909     1365_1n
SRR27950907     1369_3n
SRR27950908     1369_3n
SRR27950905     1377_3n
SRR27950906     1377_3n
SRR27950895     1440_1n
SRR27950893     1465_1n
SRR27950894     1465_1n
SRR27950892     1488_2n
SRR27950890     1490_3n
SRR27950891     1490_3n
SRR27950889     1498_1n
SRR27950887     1513_1n
SRR27950888     1513_1n
SRR27950886     1515_1n
SRR27950884     1540_3n
SRR27950885     1540_3n
SRR27950903     1543_1n
SRR27950904     1543_1n
SRR27950902     1558_1n
SRR27950901     1561_1n
SRR27950900     1570_1n
SRR27950899     1582_1n
SRR27950898     1588_1n
SRR27950897     1591_1n
SRR27950896     1593_1n
SRR27950881     1603_1n
SRR27950880     1610_1n
EOF27950882     1789_1n
EOF

# Symlink per-sample with lane indexing
SRC="$REN"                          # where SRR-based fastqs live (adjust if under $BASE)
DST="$SCR/fastqs_by_sample"
mkdir -p "$DST"

declare -A LANE
while read -r RUN SAMPLE; do
  [[ -z "$RUN" || -z "$SAMPLE" ]] && continue
  mkdir -p "$DST/$SAMPLE"
  LANE["$SAMPLE"]=$(( ${LANE["$SAMPLE"]:-0} + 1 ))
  LNUM=$(printf "L%03d" ${LANE["$SAMPLE"]})   # L001, L002, ...

  # Find R1/R2 from this run (follow symlinks if needed)
  for r in 1 2; do
    f=$(find -L "$SRC/$RUN" -maxdepth 1 -type f -name "*_R${r}_*.fastq.gz" | head -n1)
    if [ -z "$f" ]; then
      echo "WARN: Missing R${r} for $RUN" >&2
      continue
    fi
    ln -sf "$f" "$DST/$SAMPLE/${SAMPLE}_S1_${LNUM}_R${r}_001.fastq.gz"
  done
done < "$SCR/run2sample.tsv"

# Build samples.txt (unique sample names)
cut -f2 "$SCR/run2sample.tsv" | sort -u > "$SCR/samples.txt"
