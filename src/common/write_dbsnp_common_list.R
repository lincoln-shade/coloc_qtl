
#========================================================================
# Take dbSNP[version]Common.bed file and
# create file with columns rsID, chr_bp_ref_alt, and chr_bp_ref_alt_b38
# to find QTLs in ROSMAP, TOPMed, and GTEx studies
#========================================================================

library(data.table)
library(argparse)
library(stringi)

parser <- ArgumentParser()
parser$add_argument("-f", "--file", help = "File with SNP information")
parser$add_argument("--rsid_col", default = "V4", help = "rsID column name")
parser$add_argument("--bp_start_col", default = "V2", 
                    help = "variant start column (0 indexed)")
parser$add_argument("--ref_col", default = "V5", 
                    help = "reference allele column")
parser$add_argument("--alt_col", default = "V7", 
                   help = "reference allele column")
parser$add_argument("--num_alt_col", default = "V6", 
                   help = "number of alt alleles column")
parser$add_argument("--chrom_col", default = "V1", 
                    help = "chromosome column")
parser$add_argument("-o", "--out", default = "./data/snp_list.txt",
                    help = "output filename")
args <- parser$parse_args()

snp <- fread(args$file)
# lots of non-standard chromosome labels, keep only chr1-chr22, chrX, chrY
chr_labs <- c(paste0("chr", 1:22), "chrX", "chrY")
snp <- snp[get(args$chrom_col) %in% chr_labs]
# 1 alt allele
snp <- snp[get(args$num_alt_col) == 1]
# remove indels
snp[, (args$alt_col) := 
             stri_replace_last_fixed(get(args$alt_col), ",", "")]
snp <- snp[nchar(get(args$alt_col)) == 1]
snp <- snp[nchar(get(args$ref_col)) == 1]

# format output columns
# chr, bp, ref, alt, rsID, chr_bp_ref_alt, chr_bp_ref_alt_b38
snp[
    , 
    chr_bp_ref_alt := paste0(eval(parse(text=args$chrom_col)),
                             "_", 
                             eval(parse(text=args$bp_start_col)) + 1,
                             "_",
                             eval(parse(text=args$ref_col)),
                             "_",
                             eval(parse(text=args$alt_col)))
]

snp[
    ,
    chr_bp_ref_alt_b38 := paste0(chr_bp_ref_alt, "_b38")
]

setnames(snp, 
         c(args$rsid_col, args$chrom_col, args$ref_col, args$alt_col), 
         c("rsID", "chr", "ref", "alt")
)
snp[, bp := get(args$bp_start_col) + 1]
snp <- snp[, .(chr, bp, ref, alt, rsID, chr_bp_ref_alt, chr_bp_ref_alt_b38)]

# write output file, 9,978,867 variants from dbSNP153Common.bed
fwrite(snp, file = args$out, quote = FALSE)
