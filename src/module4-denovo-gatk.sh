#!/bin/bash
# annotate and sort variant for gatk
# create chr renaming file
for CHR in {{1..22},X,Y}; do
    echo ${CHR} chr${CHR}
done >> chr_names.txt
# run on multi-sample file
# annotate the multi-sample vcf
INPUT="asd.hg38_autosomes.vcf.gz"
OUTPUT="./output_gatk/asd.autosomes_chr.vcf.gz"
bcftools annotate --rename-chrs chr_names.txt ${INPUT} -Oz -o ${OUTPUT}
# create the IBD_cleaned fam file with R
write.table
# run gatk on multisample
cd output_gatk
HG38="/Users/hung/Tools/Library/HumanGenomeBuild38/Homo_sapiens_assembly38.fasta"
INPUT="asd.autosomes_chr.vcf.gz"
OUTPUT="asd.autosomes_denovo_gatk.vcf"
PED="IBD_cleaned.ped"
gatk VariantAnnotator -R ${HG38} \
  -V ${INPUT} \
  -A PossibleDeNovo \
  -ped ${PED} \
  -O ${OUTPUT}
# compress and index
bcftools convert ${OUTPUT} -Oz -o ${OUTPUT}.gz
bcftools index -t ${OUTPUT}.gz
# filter high confidence denovo variants in proband
INPUT="asd.autosomes_denovo_gatk.vcf.gz"
OUTPUT="asd.autosomes_hiconfdenovo_2.vcf.gz"
bcftools view -i 'INFO/hiConfDeNovo ~ "_1" || INFO/hiConfDeNovo ~ "_4"' ${INPUT} -Oz -o ${OUTPUT}
bcftools index ${OUTPUT}
# filter low confidence denovo variants in proband
INPUT="asd.autosomes_denovo_gatk.vcf.gz"
OUTPUT="asd.autosomes_alldenovo.vcf.gz"
bcftools view -i 'INFO/hiConfDeNovo ~ "1" || INFO/loConfDeNovo ~ "1"' ${INPUT} > ${OUTPUT}

#### ####
# filter denovo for chrX
# create chrx file from the original file (annotate chr)
INPUT="/Users/hung/Data/Autism_vinmec_coop/ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz"
CHR_ANNOT="chr_names.txt"
OUTPUT="asd.chrX.vcf.gz"
bcftools view -i 'FILTER="PASS"' ${INPUT} -Ou |\
    bcftools view -i 'CHROM ~ "X"' -Ou |\
    bcftools annotate --rename-chrs ${CHR_ANNOT} -Oz -o ${OUTPUT}
bcftools index -t ${OUTPUT}
# run gatk possibledenovo on chrx
HG38="/Users/hung/Tools/Library/HumanGenomeBuild38/Homo_sapiens_assembly38.fasta"
INPUT="asd.chrX.vcf.gz"
OUTPUT="asd.chrX_denovo_gatk.vcf"
PED="IBD_cleaned.ped"
gatk VariantAnnotator -R ${HG38} \
  -V ${INPUT} \
  -A PossibleDeNovo \
  -ped ${PED} \
  -O ${OUTPUT}
# compress and index
bcftools convert ${OUTPUT} -Oz -o ${OUTPUT}.gz
bcftools index -t ${OUTPUT}.gz

# filter high confidence denovo variants in proband
INPUT="asd.chrX_denovo_gatk.vcf.gz"
OUTPUT="asd.chrX_hiconfdenovo.vcf.gz"
bcftools view -i 'INFO/hiConfDeNovo ~ "_1" || INFO/hiConfDeNovo ~ "_4"' ${INPUT} -Oz -o ${OUTPUT}
# filter low confidence denovo variants in proband
INPUT="asd.chrX_denovo_gatk.vcf.gz"
OUTPUT="asd.chrX_alldenovo.vcf.gz"
bcftools view -i 'INFO/hiConfDeNovo ~ "1" || INFO/loConfDeNovo ~ "1"' ${INPUT} > ${OUTPUT}


#### merge hiconfdenovo ####
OUTPUT="asd.all_hiconfdenovo.vcf.gz"
bcftools concat asd.autosomes_hiconfdenovo.vcf.gz asd.chrX_hiconfdenovo.vcf.gz -Oz -o ${OUTPUT}
bcftools index -t ${OUTPUT}

