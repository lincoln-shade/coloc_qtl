#! /usr/bin/Rscript --vanilla
#=============================================================
# Input text file with SNP IDs and chromosome as two columns
# and output list of SNP IDs and chromosomes.
# can use an optional p-value threshold
#=============================================================

library(data.table)
library(argparse)

#---------------------------
# parse command-line args
#---------------------------
parser <- ArgumentParser()
parser$add_argument("-f", "--file", help = "text file w SNP and chrom cols")
parser$add_argument("--snp_col", default = "SNPID", help = "SNP column name")
parser$add_argument(
    "--snp_format", 
    default = "rsID", 
    help = paste(
        "format that the SNP IDs are in.",
        "Default is rsID.",
        "Options are rsID, chr_bp_ref_alt, or chr_bp_ref_alt_b38."
    )
)
parser$add_argument("-o", "--out", default = "./snp_list.csv",
                    help = "output filename")
parser$add_argument(
    "--chrom_col", 
    default = "chr", 
    help = paste(
        "chromosome column.",
        "Values should be formatted as \"[chr]\" or \"chr[chr]\""
    )
)
parser$add_argument(
    "--pval_col", 
    help = "p-value column for use with --max_pval"
)
parser$add_argument("--max_pval", help = "optional pval threshold for SNPs")
args <- parser$parse_args()

#--------------------------------
# main script
#--------------------------------

# apply p-valye threshold if provided
snps <- fread(args$file)
if (!is.null(args$max_pval) & !is.null(args$pval_col)) {
    max_pval <- as.numeric(args$max_pval)
    if (0 < max_pval & 1 >= max_pval) {
        snps <- snps[get(args$pval_col) <= max_pval]
    } else {
        stop("--max_pval should be between 0 and 1")
    }
}

# recode chromosome column in "chr[#]" format if needed
chr_char <- "chr"
snps[, (chr_char) := as.character(get(args$chrom_col))]
not_have_chr_char <- !grepl(chr_char, snps[[chr_char]])
snps[
    not_have_chr_char, 
    (chr_char) := paste0(chr_char, eval(parse(text=chr_char)))
]

# rename SNP column to match --snp_format
if (!(args$snp_format %in% c("rsID", "chr_bp_ref_alt", "chr_bp_ref_alt_b38"))) {
    stop("--snp_format not formatted correctly. See --help")
}

setnames(snps, args$snp_col, args$snp_format)


# write output
keep_cols <- c(chr_char, args$snp_format)
fwrite(snps[, ..keep_cols],
       file = args$out,
       quote = FALSE)

