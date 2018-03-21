#!/usr/bin/env bash
#

# Args
left_sum=data/20002_1226.nealeUKB_20170915.assoc.clean.tsv.gz
inplink=../../reference_data/1000Genomes_phase3/plink_format/EUR/1kg_p3.20130502.EUR.CHROM
chrom=1
pos=108375289
buffer=50000
threads=2

# Make range string
range=$(($pos-$buffer))-$(($pos+$buffer))

# Make covariance matrix using LDstore
mkdir -p temp
cov_prefix=temp/$chrom"_"$(echo $range | sed -e 's/[^A-Za-z0-9._]/_/g')
if [ ! -f $cov_prefix.mat.tsv ]; then
  # Make bcor and merge
  ldstore --bplink ${inplink/CHROM/$chrom} --incl-range $range --bcor $cov_prefix.bcor --n-threads $threads
  ldstore --bcor $cov_prefix.bcor --merge $threads
  rm $cov_prefix.bcor_*
  # Extract LD matrix
  ldstore --bcor $cov_prefix.bcor --matrix $cov_prefix.mat.tsv
  ldstore --bcor $cov_prefix.bcor --meta $cov_prefix.meta.txt
fi

# Run vloc
mkdir -p output
outpref=output/test
python scripts/vloc.py \
  --left_sumstats $left_sum \
  --right_sumstats $left_sum \
  --left_cov $cov_prefix.mat.tsv \
  --left_covmeta $cov_prefix.meta.txt \
  --left_kmax 1 \
  --range $chrom":"$range \
  --outprefix $outpref

rm -rf temp

# Top loci for 20002_1226:
# 1_114377568_A_G
# 12_111884608_T_C
# 9_100546600_C_T
# 2_204745003_T_C
# 12_112591686_G_A
# 1_108375289_G_A
# 3_188115682_C_A
# 13_28604007_T_C
# 6_90990050_G_A
# 2_191970120_C_G
# 13_24786576_G_A
# 12_9905690_A_G
# 6_167403873_C_A
# 1_113840826_C_T
# 22_37591290_T_G
# 4_149634572_T_C
# 10_63779871_C_T
# 11_95311422_T_C
# 8_8319963_G_C
