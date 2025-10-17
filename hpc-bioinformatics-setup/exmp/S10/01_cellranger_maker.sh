############################################################################################
# Make suitable environment and prepare files and directories for cellranger counts
############################################################################################

# Define project variables
PROJ=Bondareva_AF_sc
BASE=$LAB_ARCHIVE/public_studies/Bondareva_AF
SCR=$LAB_SCRATCH/$PROJ
REN=$SCR/fastqs_renamed        
QC=$SCR/qc
LOG=$SCR/logs

# Create necessary directories
mkdir -p "$REN" "$QC" "$LOG"

# Build sample mapping automatically from directory names
# Extract base sample names (AF1, AF10, etc.) and their variants
cd "$BASE"
> "$SCR/dir2sample.tsv"

for dir in AF*; do
    [[ ! -d "$dir" ]] && continue
    
    # Extract base sample name (remove _deep, _repeat suffixes)
    SAMPLE=$(echo "$dir" | sed -E 's/_(deep|repeat)$//')
    
    # Write mapping: directory -> base_sample
    echo -e "${dir}\t${SAMPLE}" >> "$SCR/dir2sample.tsv"
done

# Sort and preview the mapping
sort -V "$SCR/dir2sample.tsv" | tee "$SCR/dir2sample_sorted.tsv"
mv "$SCR/dir2sample_sorted.tsv" "$SCR/dir2sample.tsv"

echo "=== Sample mapping created ==="
cat "$SCR/dir2sample.tsv"

# Symlink per-sample with lane indexing
DST="$SCR/fastqs_by_sample"
mkdir -p "$DST"

declare -A LANE

while read -r DIR SAMPLE; do
    [[ -z "$DIR" || -z "$SAMPLE" ]] && continue
    
    mkdir -p "$DST/$SAMPLE"
    
    # Increment lane counter for this sample
    LANE["$SAMPLE"]=$(( ${LANE["$SAMPLE"]:-0} + 1 ))
    LNUM=$(printf "L%03d" ${LANE["$SAMPLE"]})   # L001, L002, L003, ...
    
    echo "Processing: $DIR -> $SAMPLE (Lane: $LNUM)"
    
    # Find R1/R2 fastq files in this directory
    for r in 1 2; do
        # Look for fastq.gz files with R1 or R2
        f=$(find -L "$BASE/$DIR" -maxdepth 1 -type f -name "*_R${r}_*.fastq.gz" -o -name "*_R${r}.fastq.gz" | head -n1)
        
        if [ -z "$f" ]; then
            echo "WARN: Missing R${r} for $DIR" >&2
            continue
        fi
        
        # Create symlink with cellranger-compatible naming
        ln -sf "$f" "$DST/$SAMPLE/${SAMPLE}_S1_${LNUM}_R${r}_001.fastq.gz"
    done
done < "$SCR/dir2sample.tsv"

# Build samples.txt (unique sample names)
cut -f2 "$SCR/dir2sample.tsv" | sort -u > "$SCR/samples.txt"

echo "=== Final sample list ==="
cat "$SCR/samples.txt"

echo "=== File structure preview ==="
ls -lh "$DST"/*/*.fastq.gz | head -20
