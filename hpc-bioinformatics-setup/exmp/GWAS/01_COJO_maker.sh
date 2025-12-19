# Change UKBB to COJO format
gunzip -c II_GWAS.tsv.bgz | \
awk 'BEGIN{OFS="\t"}
NR==1 {
  print "SNP","A1","A2","freq","b","se","P","N";
  next
}
{
  split($1,a,":");        # a[3]=ref, a[4]=alt
  A1=$2;                  # minor allele = effect allele
  if (A1==a[3]) A2=a[4];
  else if (A1==a[4]) A2=a[3];
  else next;              # safety check

  print $1, A1, A2, $3, $9, $10, $12, $6
}' > II_GCTA_COJO.txt
