
# coloc_qtl Intro

coloc_qtl is a collection of scripts for easy colocalization analysis automation using GWAS summary statistics and quantitative trait loci (QTL) summary statistics from 3 resources:

- Genotype Tissue Expression Project (GTEx) (https://www.gtexportal.org/home/aboutGTEx) 
- Religious Orders Study and the Rush Memory and Aging Project (ROSMAP) (https://mostafavilab.stat.ubc.ca/xqtl/)
- Trans-Omics for Precision Medicine (TOPMed) Program

Running coloc_qtl necessitates having the summary statistics of these QTL studies downloaded and available for use. GTEx and ROSMAP summary statistics are publicly available for download, while TOPMed summary statistics can be downloaded from dbGaP after obtaining both significant QTLs and full summary statistics. 

This README is a work in progress. If you have specific questions about coloc_qtl, please reach out to me at Lincoln.Shade@uky.edu.

# Installation

Clone the repository. The main interface to the pipeline is `./coloc_qtl.py`, but the scripts in `./src/` can all be run individually. Basically, after dowloading QTL data and a `.bed` file containing SNP information and running setup scripts, you will run `python3 ./coloc_qtl.py` (or make executable and run `./coloc_qtl.py`) with the appropriate option flags and arguments. Run `python3 ./coloc_qtl.py -h` to see descriptions of option flags.

## dbSNP .bed file

A reference SNP file needs to be used to harmonize variants between different data sources that don't format variants the same way. The file used when testing coloc_qtl was the dbSNP153Common.bb file from UCSC genome browser, which can be dowloaded with the below command line instructions. This file is ~1.42 GB in size.

```
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb
```

A text file is needed for coloc_qtl, so the UCSC utility bigBedToBed will need to be downloaded. UCSC utilities can be downloaded from here: https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads. Make sure you're following UCSC licenses and proper procedure when downloading. The Linux version is can be downloaded using the below command line command.

```
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
```

## GTEx

coloc_qtl has been tested only with GTEx v8 European-descended summary statistics, which can be downloaded from here: https://www.gtexportal.org/home/datasets. You will need both significant QTLs and all QTL summary statistics for whatever tissues you're interested in.

## ROSMAP

ROSMAP QTL summary statistics can be downloaded from https://mostafavilab.stat.ubc.ca/xqtl/. Again, both significant and all summary statistics are needed. coloc_qtl can only be used with eQTLs and mQTLs, not haQTLs.

## TOPMed

TOPMed QTL results do not come with significant QTLs. You will need to subset the nominal summary statistics to create a subset of significant QTLs in a separate directory with the same file names. I used a p-value threshold of 1e-5. Additionally, using the whole blood MAF $\ge 0.001$ is not recommended because coloc_qtl currently only works with common (MAF $\ge 0.01$) SNPs. 

```
./src/common/extract_variant_ids.R

./src/common/subset_snp_list.R

# depending on which data sets you want to use
./src/find_qtls_gtex.py
./src
```

# Setup and running

1. Format SNP reference file

Create a SNP reference text file. Example command-line code below.

```
bigBedToBed dbSNP153Common.bb dbSNPCommon.bed

./src/write_dbsnp_common_list.R \
    -f dbSNPCommon.bed \
    -o snp_list.txt
```

2. Run ./coloc_qtl.py with appropriate arguments. Minimal example for using TOPMed QTLs with a continuous outcome below, but please run `./coloc_qtl.py` to see help messages and default values for each argument.

```
./coloc_qtl.py \
    --ref_file snp_list.txt \
    --gwas_sumstats [your GWAS Summary statistics file] \
    --max_pval [highest GWAS pval for a variant to be checked for QTL activity. Recommended 1e-5 or lower] \ 
    --topmed \
    --eqtl_dir_topmed [directory with significant TOPMed eQTLs] \
    --eqtl_all_dir_topmed [directory with all TOPMed eQTL sumstats] \
    --out_prefix [directory for output] \
    --snp_col [SNPID column in --gwas_sumstats file] \
    --pval_col [pval column in --gwas_sumstats file] \
    --maf_col [MAF column in --gwas_sumstats file] \
    --beta_col [beta column in --gwas_sumstats file] \
    --betasd_col [beta SE column in --gwas_sumstats file] \
    --N [GWAS sample size] \
    --phenotype_type 'quant' \
    --phenotype_sd [SD of the GWAS outcome] 
```
