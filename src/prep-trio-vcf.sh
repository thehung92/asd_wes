#!/bin/bash
# INPUT="./data/ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz"
INPUT=$1
# FAMILY="ASD002_2_bwa,ASD002_3_bwa,ASD002_1_bwa" ~ mother,father,child
FAMILY=$2
LOOP=$3
mkdir -p ./temp/trio/
TEMP="./temp/trio/asd.trio"${LOOP}".vcf.gz"
OUTPUT="./temp/trio/asd.trio"${LOOP}"_QC.vcf.gz"
bcftools view -s ${FAMILY} ${INPUT} -Ou |
  bcftools view -i 'FILTER="PASS"' -Ou |
  bcftools view -M 2 -Ou |
  bcftools view -e 'CHROM ~ "_"' -Ou |
  bcftools annotate --rename-chrs ./temp/chr_names.txt -Oz -o ${TEMP}
# argument is: mother, father, child
# QC the output
bcftools +mendelian ${TEMP} -m d -t ${FAMILY} | # delete GT with mendelian violation
  bcftools view -g ^miss -Ou | # rm sites with missing genotype
  bcftools view -i 'INFO/AC>0' -Oz -o ${OUTPUT} # rm sites without alternate allele
#
rm ${TEMP}