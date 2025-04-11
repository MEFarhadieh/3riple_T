#!/bin/bash

##HyperParameters
snp_miss_th1=0.2
snp_miss_th2=0.02
ind_miss_th1=0.2
ind_miss_th2=0.02
maf_th=0.05
hwe_th1=1e-6
hwe_th2=1e-10
pihat_th=0.2

##Paths
root_path="/Users/erfan/Desktop/fututre_projects/QTL_MR/FDG_PET_GWAS/"
work_path=$root_path"QualityControl/"
mkdir -p $work_path
data_path=$root_path"ADNI_test"
cov_path=$root_path"cov_pheno.txt"
final_path=$work_path
code_path=$root_path"scripts/"
#Utility Scripts path
QC_path=$code_path"QC_GWAS/"

#############################################################################################
# 0 step
#############################################################################################
#merge data
plink --bfile ADNI_1_GWAS_Plink/ADNI_cluster_01_forward_757LONI \
      --bmerge ADNI_GO_2_OmniExpress/ADNI_GO_2_Forward_Bin.bed \
               ADNI_GO_2_OmniExpress/ADNI_GO_2_Forward_Bin.bim \
               ADNI_GO_2_OmniExpress/ADNI_GO_2_Forward_Bin.fam \
      --make-bed --out GWAS_1_2

plink --bfile GWAS_1_2 \
      --bmerge ADNI_GO_2_2nd/ADNI_GO2_GWAS_2nd_orig_BIN.bed \
               ADNI_GO_2_2nd/ADNI_GO2_GWAS_2nd_orig_BIN.bim \
               ADNI_GO_2_2nd/ADNI_GO2_GWAS_2nd_orig_BIN.fam \
      --make-bed --out GWAS_1_2_2

plink --bfile GWAS_1_2_2 \
      --bmerge ADNI3_PLINKFinal/ADNI3_PLINK_Final.bed \
               ADNI3_PLINKFinal/ADNI3_PLINK_Final.bim \
               ADNI3_PLINKFinal/ADNI3_PLINK_Final.fam \
      --make-bed --out GWAS_1_2_3

# Remove mis-matching SNPs and merge again
#rs10250779
#rs16910526
#rs17107315
#rs17602729
#rs2274083
#rs35067814
#rs3825942
plink --bfile GWAS_1_2_2 --exclude GWAS_1_2_3-merge.missnp --make-bed --out GWAS_1_2_2

plink --bfile ADNI3_PLINKFinal/ADNI3_PLINK_Final \
      --exclude GWAS_1_2_3-merge.missnp \
      --make-bed --out ADNI3_clean

plink --bfile GWAS_1_2_2 \
      --bmerge ADNI3_clean \
      --make-bed --out GWAS_3

plink --bfile GWAS_3 \
      --bmerge ADNI3_PLINKFinal_2nd/ADNI3_PLINK_FINAL_2nd.bed \
               ADNI3_PLINKFinal_2nd/ADNI3_PLINK_FINAL_2nd.bim \
               ADNI3_PLINKFinal_2nd/ADNI3_PLINK_FINAL_2nd.fam \
      --make-bed --out GWAS_total

# Remove mis-matching SNPs and merge again
plink --bfile GWAS_3 --exclude GWAS_total-merge.missnp --make-bed --out GWAS_3

plink --bfile ADNI3_PLINKFinal_2nd/ADNI3_PLINK_FINAL_2nd \
      --exclude GWAS_total-merge.missnp \
      --make-bed --out ADNI3_2_clean

plink --bfile GWAS_3 \
      --bmerge ADNI3_2_clean \
      --make-bed --out GWAS_final

#############################################################################################
# 1st step
#############################################################################################
cd $work_path

# Investigate missingness per individual and per SNP and make histograms
plink --bfile $data_path --pheno $cov_path --pheno-name DIAG --missing --noweb

# Generate plots to visualize the missingness results
Rscript --no-save  $QC_path/hist_miss.R --noweb

# Delete SNPs with missingness >0.02
plink --bfile $data_path --pheno $cov_path --pheno-name DIAG --geno 0.02 --mind 0.05 --make-bed --out "GWAS_clean_2" --noweb

#############################################################################################
# 2nd step
#############################################################################################

# Check for sex discrepancy
plink --bfile "GWAS_clean_2" --pheno $cov_path --pheno-name DIAG --check-sex --noweb

# Generate plots to visualize the sex-check results
Rscript --no-save $QC_path/gender_check.R

# Delete individuals with sex discrepancy
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt

# impute-sex
plink --bfile "GWAS_clean_2" --pheno $cov_path --pheno-name DIAG --impute-sex --make-bed --out "GWAS_clean_3" --noweb

#############################################################################################
# 3rd step
#############################################################################################

# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF)

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' "GWAS_clean_3".bim > snp_1_22.txt
plink --bfile "GWAS_clean_3" --pheno $cov_path --pheno-name DIAG --extract snp_1_22.txt --make-bed --out "GWAS_clean_4" --noweb

# Generate a plot of the MAF distribution.
plink --bfile "GWAS_clean_4" --pheno $cov_path --pheno-name DIAG --freq --out MAF_check --noweb
Rscript --no-save $QC_path/MAF_check.R

# Remove SNPs with a low MAF frequency.
plink --bfile "GWAS_clean_4" --pheno $cov_path --pheno-name DIAG --maf $maf_th --make-bed --out "GWAS_clean_5" --noweb

#############################################################################################
# 4th step
#############################################################################################

# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE)
plink --bfile "GWAS_clean_5" --pheno $cov_path --pheno-name DIAG --hardy --noweb

# Selecting SNPs with HWE p-value below 0.00001, required for one of the two plot generated by the next Rscript, allows to zoom in on strongly deviating SNPs
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
Rscript --no-save $QC_path/hwe.R

plink --bfile "GWAS_clean_5" --pheno $cov_path --pheno-name DIAG --hwe $hwe_th1 --make-bed --out "GWAS_1_2_3_clean_"$analysis"hwe_filter_step1" --noweb
plink --bfile  "GWAS_1_2_3_clean_"$analysis"hwe_filter_step1" --pheno $cov_path --pheno-name DIAG --hwe $hwe_th2 --hwe-all --make-bed --out "GWAS_1_2_3_clean_"$analysis"9" --noweb
