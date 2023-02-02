#! /usr/bin/Rscript --vanilla
library(data.table)
library(magrittr)
library(stringi)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--file", "-f", default = "./rosmap_qtl_test.txt")
parser$add_argument("--out", "-o", default = "./rosmap_qtl_format_test.txt")
args <- parser$parse_args()

qtls <- fread(args$file, header = FALSE)
setnames(
    qtls, 
    colnames(qtls), 
    c(
        "filename", "qtl_type", "rsID", "chr", "pos", "A1", "A2", "phenotype_id", 
        "CpG_pos_or_gene_TSS", "beta", "beta_sd", "t", "pvalues"
    )
)

qtls[, tissue := "Brain_Prefrontal_Cortex"]
qtls[, chr := paste0("chr", chr)]
# select only the top QTL for each phenotype to keep in output
top_qtl_index <- 
    qtls[, .I[pvalues == min(pvalues)], .(phenotype_id, tissue)] %>%
    .[, V1]
top_qtls <- qtls[top_qtl_index]
setcolorder(top_qtls, c("phenotype_id", "tissue", "chr", "qtl_type"))

# some variants in perfect LD have same pval. Remove all but one
dup_rows <- top_qtls[, .I[duplicated(pvalues)], .(phenotype_id, tissue)]$V1
top_qtls <- top_qtls[!dup_rows]

# write output
fwrite(top_qtls, file = args$out, sep = " ", quote = FALSE)
