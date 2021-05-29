# This project will analyze WES data from ASD cohort

## create data structure
```
# populate directory structure
mkdir -p data docs output src temp/data temp/subset

├── README.md
├── data
│   ├── ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz
│   ├── ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz.csi
│   ├── README.md
│   ├── Sample-information.xlsx
│   ├── asd.hg38_QC.vcf.gz
│   ├── asd.hg38_QC.vcf.gz.csi
│   └── high-ld-region-hg38.txt
├── docs
│   ├── QC_summary.xlsx
│   ├── table1.csv
│   ├── tdt_annotate_vep.R
│   └── tdt_manhattan.R
├── output
│   ├── README.md
│   ├── asd.hg38_ibd-clean.bed
│   ├── asd.hg38_ibd-clean.bim
│   ├── asd.hg38_ibd-clean.fam
│   ├── asd.hg38_ibd-clean.hh
│   ├── asd.hg38_ibd-clean.log
│   └── tdt
│       ├── asd.290_ibd-clean.hh
│       ├── asd.290_ibd-clean.log
│       ├── asd.290_ibd-clean.tdt
│       └── asd.290_ibd-clean.tdt.adjusted
├── src
│   ├── exclude-high-ld-region.R
│   ├── github-src.md
│   ├── highLDregions4bim_b38.awk
│   ├── ibd-filter.R
│   ├── infer-fam.R
│   ├── module1-QC-data.sh
│   ├── module2-tdt.sh
│   └── setup-docker.sh
└── temp
    ├── README.md
    ├── data
    │   ├── README.md
    │   ├── asd.hg38_QC2.bed
    │   ├── asd.hg38_QC2.bim
    │   ├── asd.hg38_QC2.fam
    │   ├── asd.hg38_QC2.fam.bu
    │   ├── asd.hg38_QC2.log
    │   ├── asd.hg38_common.bed
    │   ├── asd.hg38_common.bim
    │   ├── asd.hg38_common.fam
    │   ├── asd.hg38_common.hh
    │   ├── asd.hg38_common.log
    │   ├── asd.hg38_common.nosex
    │   ├── asd.hg38_filter.bed
    │   ├── asd.hg38_filter.bim
    │   ├── asd.hg38_filter.fam
    │   ├── asd.hg38_filter.hh
    │   ├── asd.hg38_filter.log
    │   ├── asd.hg38_filter.nosex
    │   ├── asd.hg38_hwe.bed
    │   ├── asd.hg38_hwe.bim
    │   ├── asd.hg38_hwe.fam
    │   ├── asd.hg38_hwe.hh
    │   ├── asd.hg38_hwe.log
    │   ├── asd.hg38_hwe.nosex
    │   ├── asd.hg38_sexcheck-clean.bed
    │   ├── asd.hg38_sexcheck-clean.bim
    │   ├── asd.hg38_sexcheck-clean.fam
    │   ├── asd.hg38_sexcheck-clean.fam.bu
    │   ├── asd.hg38_sexcheck-clean.hh
    │   ├── asd.hg38_sexcheck-clean.log
    │   ├── asd.hg38_sexcheck-clean.nosex
    │   └── discordant_sex.list
    └── subset
        ├── README.md
        ├── asd.hg38_IBD.bed
        ├── asd.hg38_IBD.bim
        ├── asd.hg38_IBD.fam
        ├── asd.hg38_IBD.genome
        ├── asd.hg38_IBD.genome.filter
        ├── asd.hg38_IBD.log
        ├── asd.hg38_IBD.nosex
        ├── asd.hg38_autosome.bed
        ├── asd.hg38_autosome.bim
        ├── asd.hg38_autosome.fam
        ├── asd.hg38_autosome.log
        ├── asd.hg38_autosome.nosex
        ├── asd.hg38_pairwise.hh
        ├── asd.hg38_pairwise.log
        ├── asd.hg38_pairwise.nosex
        ├── asd.hg38_pairwise.prune.in
        ├── asd.hg38_pairwise.prune.out
        ├── asd.hg38_prune-sexcheck-clean.bed
        ├── asd.hg38_prune-sexcheck-clean.bim
        ├── asd.hg38_prune-sexcheck-clean.fam
        ├── asd.hg38_prune-sexcheck-clean.fam.bu
        ├── asd.hg38_prune-sexcheck-clean.log
        ├── asd.hg38_prune-sexcheck-clean.nosex
        ├── asd.hg38_prune.bed
        ├── asd.hg38_prune.bim
        ├── asd.hg38_prune.fam
        ├── asd.hg38_prune.hh
        ├── asd.hg38_prune.log
        ├── asd.hg38_prune.nosex
        ├── asd.hg38_sexcheck.hh
        ├── asd.hg38_sexcheck.log
        ├── asd.hg38_sexcheck.nosex
        ├── asd.hg38_sexcheck.sexcheck
        ├── asd.hg38_split.bed
        ├── asd.hg38_split.bim
        ├── asd.hg38_split.fam
        ├── asd.hg38_split.log
        ├── discordant_ibd.list
        ├── discordant_sex.txt
        ├── exclude-highld-nonautosome.list
        ├── temp1.txt
        └── temp2.txt

8 directories, 106 files
```