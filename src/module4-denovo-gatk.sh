#!/bin/bash
# annotate and sort variant for gatk
# create chr renaming file
for CHR in {{1..22},X,Y}; do
    echo ${CHR} chr${CHR}
done >> ./temp/chr_names.txt
# run on multi-sample file
# annotate the multi-sample vcf
INPUT="data/asd.hg38_QC.vcf.gz"
ARG="./temp/chr_names.txt"
OUTPUT="temp/data/asd.hg38_rename.vcf.gz"
bcftools annotate --rename-chrs $ARG ${INPUT} -Oz -o ${OUTPUT}
bcftools index -t $OUTPUT
# create the pedigree file with R for vcf
Rscript src/create-ped-gatk.R
# run gatk VariantAnnotator on multisample
mkdir -p output/denovo
INPUT="temp/data/asd.hg38_rename.vcf.gz"
ARG1="/home/reference/GRCh38/Homo_sapiens_assembly38.fasta"
ARG2="./temp/data/asd.284_ped-gatk.txt"
OUTPUT="./output/denovo/asd.284_denovo-gatk.vcf.gz"
gatk VariantAnnotator -R ${ARG1} \
  -V ${INPUT} \
  -A PossibleDeNovo \
  -ped ${ARG2} \
  -O ${OUTPUT}
# filter high confidence denovo variants in proband
INPUT="./output/denovo/asd.284_denovo-gatk.vcf.gz"
OUTPUT="./output/denovo/asd.284_hiconfdenovo.vcf.gz"
bcftools view -i 'INFO/hiConfDeNovo ~ "_1" || INFO/hiConfDeNovo ~ "_4"' ${INPUT} -Oz -o ${OUTPUT}
bcftools index ${OUTPUT}
# extract info for annotate with VEP online tools
INPUT="./output/denovo/asd.284_hiconfdenovo.vcf.gz"
OUTPUT="./output/denovo/asd.284_denovo-var.vcf"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $INPUT > $OUTPUT
# annotate "./output/denovo/asd.284_denovo-var.vcf" with VEP online tools
# we will get ./output/denovo/asd.284_denovo-vep.txt