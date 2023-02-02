#! /usr/bin/Rscript --vanilla

#===============================================================
# subsets variant info file using list of variant ids provided
#===============================================================

library(data.table)
library(argparse)

#---------------------------
# parse command-line args
#---------------------------
parser <- ArgumentParser()
parser$add_argument("--ref_file", "-r", default = "data/snp_list.txt")
parser$add_argument("--query_file", "-q")
parser$add_argument("--snp_format", default = "rsID")
parser$add_argument("--out", "-o", default = "./snp_query.csv")
args <- parser$parse_args()

#--------------------------------
# main script
#--------------------------------

# read in files
ref <- fread(args$ref_file)
query <- fread(args$query_file)

# merge
snps <- merge(ref, query, c("chr", args$snp_format))

# write output
fwrite(snps, file = args$out, quote = FALSE)
