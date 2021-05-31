#!/usr/bin/env bash
# cd to woking dir
cd /home/public
# populate dir structure
mkdir -p data docs output src temp/data temp/subset
# filter pass
# remove variants called on contigs
# split multi-allelic
# add tag MAF, F_missing
INPUT="./data/ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz"
OUTPUT="./data/asd.hg38_QC.vcf.gz"
bcftools view -i 'FILTER="PASS"' ${INPUT} -Ou |\
    bcftools view -e 'CHROM ~ "_"' -Ou |\
    bcftools +fill-tags  -Oz -o ${OUTPUT} -- -t MAF,F_MISSING
bcftools index ${OUTPUT}
# work with biallelic variants only in plink format
INPUT="./data/asd.hg38_QC.vcf.gz"
OUTPUT="./temp/data/asd.hg38_QC2"
plink2 --vcf ${INPUT} \
	--id-delim --allow-extra-chr --chr 1-22,X,Y --max-alleles 2 \
	--set-missing-var-ids @:# \
	--make-bed --out ${OUTPUT}
# backup famfile
mv ./temp/data/asd.hg38_QC2.fam ./temp/data/asd.hg38_QC2.fam.bu
# use Rscript to generate fam file with metadata from excel file
EXCEL="./data/Sample-information.xlsx"
FAM="./temp/data/asd.hg38_QC2.fam.bu"
OUTPUT="./temp/data/r-output.fam"
Rscript src/infer-fam.R $EXCEL $FAM $OUTPUT
# replace asd.hg38_QC2.fam
mv ./temp/data/r-output.fam ./temp/data/asd.hg38_QC2.fam

# filter for common SNP
INPUT="./temp/data/asd.hg38_QC2"
OUTPUT="./temp/data/asd.hg38_common"
plink --bfile $INPUT \
    --maf 0.01 \
    --make-bed --out $OUTPUT
# filter missingness at 95% per person and 95% per SNP
INPUT="./temp/data/asd.hg38_common"
OUTPUT="./temp/data/asd.hg38_filter"
plink --bfile $INPUT \
    --geno 0.05 --mind 0.95 \
    --make-bed --out $OUTPUT
# remove SNP deviate from hwe 1*10^-5
INPUT="./temp/data/asd.hg38_filter"
OUTPUT="./temp/data/asd.hg38_hwe"
plink --bfile $INPUT \
    --hwe 0.00001 \
    --make-bed --out $OUTPUT

# prune snp by LD
INPUT="./temp/data/asd.hg38_hwe"
OUTPUT="./temp/subset/asd.hg38_pairwise"
plink --bfile $INPUT --indep-pairwise 1500 150 0.2 --out $OUTPUT
OUTPUT2="./temp/subset/asd.hg38_prune"
plink --bfile $INPUT --extract ${OUTPUT}.prune.in --make-bed --out $OUTPUT2

# exclude SNP in high LD regions and in non-autosomal region
INPUT="./temp/subset/asd.hg38_prune"
TEMPDIR="./temp/subset/"
OUTPUT="./temp/subset/exclude-highld-nonautosome.list"
awk -f src/highLDregions4bim_b38.awk ${INPUT}.bim > ${TEMPDIR}/temp1.txt
awk '($1 < 1) || ($1 > 22) {print $2}' ${INPUT}.bim > ${TEMPDIR}/temp2.txt
cat ${TEMPDIR}/temp1.txt ${TEMPDIR}/temp2.txt > $OUTPUT
OUTPUT2="./temp/subset/asd.hg38_autosome"
plink --bfile $INPUT --exclude $OUTPUT --make-bed --out $OUTPUT2

# check sex with asd.hg38_prune
INPUT="./temp/subset/asd.hg38_prune"
OUTPUT="./temp/subset/asd.hg38_split"
plink2 --bfile $INPUT --split-par b38 --make-bed --out $OUTPUT # remove pseudo-autosome region with plink2
OUTPUT2="./temp/subset/asd.hg38_sexcheck"
plink --bfile $OUTPUT --check-sex ycount 0.3 0.8 5 10 --out $OUTPUT2 # check sex with plink
OUTPUT3="""./temp/subset/discordant_sex.txt"
awk '$5 ~ /PROBLEM/{print $0}' ${OUTPUT2}.sexcheck > $OUTPUT3 # print problem ID

# manually select individual to exclude from `cat ./temp/subset/discordant_sex.txt`
# ASD034, ASD043, ASD101 family will be excluded because F-stats and Y-snpcount do not match in the proband
# ASD014 sex will be updated to 1
# extract FID & IID to be excluded
ARG1="./temp/subset/asd.hg38_autosome"
ARG2="./temp/data/discordant_sex.list"
awk '($1~/ASD034/) || ($1~/ASD043/) || ($1~/ASD101/) {print $1, $2}' ${ARG1}.fam > ${ARG2}
# edited in QC dataset
INPUT="./temp/data/asd.hg38_hwe"
OUTPUT="./temp/data/asd.hg38_sexcheck-clean"
plink --bfile $INPUT --remove $ARG2 --make-bed --out $OUTPUT
mv $OUTPUT.fam $OUTPUT.fam.bu # backup the fam file
awk '$2=="ASD014_1" {gsub(/./, "1", $5)} 1' $OUTPUT.fam.bu > $OUTPUT.fam # update sex of ASD014_1
# edited in pruned dataset
INPUT="./temp/subset/asd.hg38_autosome"
OUTPUT="./temp/subset/asd.hg38_prune-sexcheck-clean"
plink --bfile $INPUT --remove $ARG2 --make-bed --out $OUTPUT # pruned dataset
mv $OUTPUT.fam $OUTPUT.fam.bu # backup the fam file
awk '$2=="ASD014_1" {gsub(/./, "1", $5)} 1' $OUTPUT.fam.bu > $OUTPUT.fam # update sex of ASD014_1

# Pairwise identical-by-descent (IBD) check all the samples from their pruned autosome SNP
INPUT="./temp/subset/asd.hg38_prune-sexcheck-clean"
OUTPUT="./temp/subset/asd.hg38_IBD"
plink --bfile $INPUT --genome --make-bed --out $OUTPUT
# filter the pi_hat value between individuals from the same family with Rscript
Rscript src/ibd-filter.R $OUTPUT.genome $OUTPUT.genome.filter
# extract FID & IID of sample to be excluded
OUTPUT2="./temp/subset/discordant_ibd.list"
awk ' NR>1 {print $3, $4}' $OUTPUT.genome.filter > $OUTPUT2

# remove discordant_ibd sample from dataset and output
INPUT="./temp/data/asd.hg38_sexcheck-clean"
ARG="./temp/subset/discordant_ibd.list"
OUTPUT="./output/asd.hg38_ibd-clean"
plink --bfile $INPUT --remove $ARG --make-bed --out $OUTPUT