# This project will analyze WES data from ASD cohort

## create data structure
```
# populate directory structure
mkdir -p data docs output src temp/data temp/subset
```

## content of project dir after running on macOS
```
├── README.md
├── asd_wes.Rproj
├── data
│   ├── ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz
│   ├── ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz.csi
│   ├── README.md
│   ├── SFARI-Gene_genes_10-29-2020release_11-04-2020export.csv
│   ├── Sample-information.xlsx
│   ├── asd.hg38_QC.vcf.gz
│   ├── asd.hg38_QC.vcf.gz.csi
│   ├── gene_mutationrate.xls
│   ├── genetic_map_GRCh38_merged.tab.gz
│   ├── matrix.r2.unrelated.variant1.txt
│   ├── matrix.r2.unrelated.variant2.txt
│   ├── matrix.r2.unrelated.variant3.txt
│   └── matrix.r2.unrelated.variant4.txt
├── docker
│   ├── Dockerfile
│   └── install-R-lib.R
├── docs
│   └── QC_summary.xlsx
├── output
│   ├── README.md
│   ├── asd.hg38_ibd-clean.bed
│   ├── asd.hg38_ibd-clean.bim
│   ├── asd.hg38_ibd-clean.fam
│   ├── asd.hg38_ibd-clean.hh
│   ├── asd.hg38_ibd-clean.log
│   ├── denovo
│   │   ├── asd.284_denovo-gatk.vcf.gz
│   │   ├── asd.284_denovo-gatk.vcf.gz.tbi
│   │   ├── asd.284_denovo-var.vcf
│   │   ├── asd.284_denovo-vep.txt
│   │   ├── asd.284_hiconfdenovo.vcf.gz
│   │   └── asd.284_hiconfdenovo.vcf.gz.csi
│   ├── plot-table
│   │   ├── figure2.pdf
│   │   ├── figure3.png
│   │   ├── figureS1.pdf
│   │   ├── figureS2.png
│   │   ├── figureS3.png
│   │   ├── figureS4.png
│   │   ├── table1.csv
│   │   ├── table2.csv
│   │   ├── tables1.csv
│   │   └── tables2.csv
│   └── tdt
│       ├── asd.290_ibd-clean.hh
│       ├── asd.290_ibd-clean.log
│       ├── asd.290_ibd-clean.tdt
│       ├── asd.290_ibd-clean.tdt.adjusted
│       ├── asd.290_tdt-sig-var.vcf
│       ├── asd.290_tdt-vep.txt
│       └── regions.txt
├── src
│   ├── TADA.v.1.2.R
│   ├── create-ped-gatk.R
│   ├── exclude-high-ld-region.R
│   ├── function_plot-gene-region.R
│   ├── github-src.md
│   ├── highLDregions4bim_b38.awk
│   ├── ibd-filter.R
│   ├── infer-fam.R
│   ├── module1-QC-data.sh
│   ├── module2-tdt.sh
│   ├── module3-tdt-analysis.R
│   ├── module4-compute-ld.R
│   ├── module4-denovo-gatk.sh
│   ├── module5-plot-gene-region.R
│   ├── module6-gprofiler.R
│   ├── module7-tada.R
│   ├── prep-trio-vcf.sh
│   └── setup-docker.sh
└── temp
    ├── README.md
    ├── chr_names.txt
    ├── data
    │   ├── asd.284_ped-gatk.txt
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
    │   ├── asd.hg38_rename.vcf.gz
    │   ├── asd.hg38_rename.vcf.gz.tbi
    │   ├── asd.hg38_sexcheck-clean.bed
    │   ├── asd.hg38_sexcheck-clean.bim
    │   ├── asd.hg38_sexcheck-clean.fam
    │   ├── asd.hg38_sexcheck-clean.fam.bu
    │   ├── asd.hg38_sexcheck-clean.hh
    │   ├── asd.hg38_sexcheck-clean.log
    │   ├── asd.hg38_sexcheck-clean.nosex
    │   └── discordant_sex.list
    ├── ld-r2
    │   ├── asd.unrelated.fam
    │   ├── matrix.r2.unrelated.variant1.txt
    │   ├── matrix.r2.unrelated.variant2.txt
    │   ├── matrix.r2.unrelated.variant3.txt
    │   ├── matrix.r2.unrelated.variant4.txt
    │   ├── variant1-region.txt
    │   ├── variant1-region_ldr2.hh
    │   ├── variant1-region_ldr2.ld
    │   ├── variant1-region_ldr2.log
    │   ├── variant1-region_ldr2.nosex
    │   ├── variant2-region.txt
    │   ├── variant2-region_ldr2.ld
    │   ├── variant2-region_ldr2.log
    │   ├── variant2-region_ldr2.nosex
    │   ├── variant3-region.txt
    │   ├── variant3-region_ldr2.hh
    │   ├── variant3-region_ldr2.ld
    │   ├── variant3-region_ldr2.log
    │   ├── variant3-region_ldr2.nosex
    │   ├── variant4-region.txt
    │   ├── variant4-region_ldr2.ld
    │   ├── variant4-region_ldr2.log
    │   └── variant4-region_ldr2.nosex
    └── subset
        ├── asd.hg38_IBD.bed
        ├── asd.hg38_IBD.bim
        ├── asd.hg38_IBD.fam
        ├── asd.hg38_IBD.genome
        ├── asd.hg38_IBD.genome.filter
        ├── asd.hg38_IBD.log
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

12 directories, 165 files

```
