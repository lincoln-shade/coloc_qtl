#! /usr/bin/Rscript --vanilla
#=====================================
# Tidy GTEx QTL files
#=====================================

library(data.table)
library(magrittr)
library(stringi)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--file", "-f", default = "./tmp/gtex_qtl_test.txt")
parser$add_argument("--out", "-o", default = "./tmp/gtex_qtl_tidy.txt")
args <- parser$parse_args()

qtls <- fread(args$file, header = FALSE)
setnames(
  qtls, 
  colnames(qtls), 
  c("filename", "phenotype_id", "variant_id",	"tss_distance",	"maf",
    "ma_samples", "ma_count",	"pval_nominal",	"slope", "slope_se",
    "pval_nominal_threshold", "min_pval_nominal", "pval_beta"
    )
)

# create variable for QTL type
qtls[!(grep("clu", phenotype_id)), qtl_type := "eQTL"]
qtls[grep("clu", phenotype_id), qtl_type := "sQTL"]

# create variable for tissue type
qtls[
  , 
  tissue := stri_replace_last_fixed(filename, ".v8.EUR.signif_pairs.txt", "")
]

# create variable for chromosome
qtls[, chr := stri_replace_first_regex(variant_id, "_.*", "")]

# select only the top QTL for each phenotype to keep in output
top_qtl_index <- 
  qtls[, .I[pval_nominal == min(pval_nominal)], .(phenotype_id, tissue)] %>%
  .[, V1]
top_qtls <- qtls[top_qtl_index]
setcolorder(top_qtls, c("phenotype_id", "tissue", "chr", "qtl_type"))

# some variants in perfect LD have same pval. Remove all but one
dup_rows <- top_qtls[, .I[duplicated(pval_nominal)], .(phenotype_id, tissue)]$V1
top_qtls <- top_qtls[!dup_rows]


# write output
fwrite(top_qtls, file = args$out, sep = " ", col.names = TRUE, quote = FALSE)
