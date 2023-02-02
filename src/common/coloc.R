#! /usr/bin/Rscript --vanilla

#===============================================================
# Single-variant colocalization analysis
#===============================================================

library(data.table)
library(magrittr)
library(coloc)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-g", "--gwas", 
                    help="GWAS sumstats .Rds file")
parser$add_argument("-q", "--qtl", 
                    help="QTL sumstats .Rds file")
parser$add_argument("-o", "--out",
                    default = "./coloc_results.txt",
                    help="output file path")
parser$add_argument("--p12", default = "1e-5", 
                    help="prior probability of colocalization")
args <- parser$parse_args()

p12 <- as.numeric(args$p12)
gwas <- readRDS(args$gwas)
qtl <- readRDS(args$qtl)

# # format to avoid errors with coloc.abf function
if (is.null(qtl$metadata$sdY)) {
  qtl$metadata$sdY <- NULL
}

if (is.null(gwas$metadata$sdY)) {
  gwas$metadata$sdY <- NULL
}
# 
if (gwas$metadata$type == "quant") {
  gwas$metadata$s <- NULL
}

# subset so gwas and qtl sumstats have same variants and variant id format
# setnames(gwas$dataset, qtl$metadata$variant_id_format, "snp")
gwas$dataset[, snp := get(qtl$metadata$variant_id_format)]
qtl$dataset <- qtl$dataset[snp %in% gwas$dataset$snp]
gwas$dataset <- gwas$dataset[snp %in% qtl$dataset$snp]

# check if MAF is present in QTL data set.
# If not, use MAF from GWAS data set.
if (is.null(qtl$dataset[["MAF"]])) {
  qtl$dataset <- merge(qtl$dataset, gwas$dataset[, .(snp, MAF)], "snp")
}

# perform colocalization analysis
make_coloc_dataset <- function(l) {
  ds <- c(l$dataset, l$metadata)
}

results <- coloc.abf(
  dataset1 = make_coloc_dataset(gwas), 
  dataset2 = make_coloc_dataset(qtl), 
  p12 = p12
)

# write output
posteriors <- formatC(unname(results$summary))
top_snp <- gwas$dataset[pvalues == min(pvalues), snp]
gwas_top <- gwas$dataset[snp == top_snp]

qtl_top <- qtl$dataset[snp == top_snp]

cols <- c('GWAS_Phenotype', 'QTL_Source', 'QTL_Type', 'QTL_Phenotype', 
          'Tissue', 'Chromosome', 'rsID', 'P_GWAS', 'beta_QTL', 'P_QTL', 
          'NSNP', 'PPH0', 'PPH1', 'PPH2', 'PPH3', 'PPH4')

out <- c(gwas$metadata$phenotype, qtl$metadata$data_source, 
         qtl$metadata$qtl_type, qtl$metadata$phenotype, qtl$metadata$tissue,
         gwas$dataset$chr[1], gwas_top$rsID, gwas_top$pvalues, qtl_top$beta,
         qtl_top$pvalues, posteriors)
if (!(file.exists(args$out))) {
  write(cols,
        file = args$out,
        ncolumns = length(cols))
}
write(out, 
      file = args$out, 
      ncolumns = length(cols), 
      append = TRUE)
