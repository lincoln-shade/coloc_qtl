
# coloc_qtl Intro

coloc_qtl is a collection of scripts for easy colocalization analysis automation using GWAS summary statistics and quantitative trait loci (QTL) summary statistics from 3 resources:

- Genotype Tissue Expression Project (GTEx) (https://www.gtexportal.org/home/aboutGTEx) 
- Religious Orders Study and the Rush Memory and Aging Project (ROSMAP) (https://mostafavilab.stat.ubc.ca/xqtl/)
- Trans-Omics for Precision Medicine (TOPMed) Program

Running coloc_qtl necessitates having the summary statistics of these QTL studies downloaded and available for use. GTEx and ROSMAP summary statistics are publicly available for download, while TOPMed summary statistics can be downloaded from dbGaP after obtaining both significant QTLs and full summary statistics. 

# Installation

Install the repository. The main interface to the pipeline is `./coloc_qtl.py`, but the scripts in `./src/` can all be run individually. Basically, you will run `python3 ./coloc_qtl.py` with the appropriate option flags and arguments. Run `python3 ./coloc_qtl.py -h` to see 

# GTEx

coloc_qtl has been tested only with GTEx v8 European-descended summary statistics. 

# ROSMAP

coloc_qtl can only be used with eQTLs and mQTLs, not haQTLs.

# TOPMed

TOPMed QTL results do not come with significant QTLs. You will need to subset the nominal summary statistics to create a subset of significant QTLs in a separate directory with the same file names. I used a p-value threshold of 1e-5. Additionally, using the whole blood MAF $\ge 0.001$ is not recommended because coloc_qtl currently only works with common (MAF $\ge 0.01$) SNPs. 

```
./src/common/extract_variant_ids.R

./src/common/subset_snp_list.R

# depending on which data sets you want to use
./src/find_qtls_gtex.py
./src
```