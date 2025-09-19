# How to Run
---
### 1- Define the vraiables:
Change the below block base on the projects.
```
SAMPLE=pbmc_unsorted_10k
LIBS=/PATH_to_the_file/libraries.csv
REF=/PATH_to_the_file/refdata-cellranger-arc-GRCh38-2024-A
```
---
### 2- Run for one sample:
```
nextflow run main.nf -profile slurm \
  --sample_id $SAMPLE \
  --libraries $LIBS \
  --reference $REF
```
---
### 3- Run for multiple sample:
```
  nextflow run multi-sample_main.nf -profile slurm \
  --samplesheet /PATH_to_the_file/samplesheet.csv \
  --reference $REF \
  --outdir_root /PATH_to_the_output/$SAMPLE \
  --max_parallel 20
```






