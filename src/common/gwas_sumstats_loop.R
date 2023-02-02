#! /usr/bin/Rscript --vanilla
#==================================================================
# Loop through a set of QTLs and subset GWAS sumstats to perform 
# colocalization analysis with each.
#==================================================================

library(data.table)
library(argparse)
library(parallel)
library(stringi)

parser <- ArgumentParser()
parser$add_argument("--maf_col", default = "MAF")
parser$add_argument("--pval_col", default = "pval.spa")
parser$add_argument("--snpid_col", default = "SNPID")
parser$add_argument("--beta_col")
parser$add_argument("--betavar_col")
parser$add_argument("--betasd_col")
parser$add_argument("--radius", default = "2e5")
parser$add_argument("-r", "--results", help="path to GWAS results file")
parser$add_argument(
    "--qtl_format", 
    help = "rsID, chr_bp_ref_alt, or chr_bp_ref_alt_b38 allowed"
)
parser$add_argument(
    "--snp_format", 
    default = "rsID", 
    help = "rsID, chr_bp_ref_alt, or chr_bp_ref_alt_b38 allowed")
parser$add_argument("--out_prefix", "-o")
parser$add_argument(
    "--snp_file", "-s", 
    default = "./data/snp_list.txt",
    help = "SNP file that's output of src/common/write_dbsbp_common_list.R"
)

# arguments for QTL file
parser$add_argument("--qtl_file")
parser$add_argument("--qtl_id_col")
parser$add_argument("--qtl_chr_col")

# arguments for determining features of phenotype for colocalization
parser$add_argument("--phenotype_file")
parser$add_argument(
    "--phenotype_col", default = "GWAS",
    help = "phenotype column if --phenotype_file given, or pheno name if not"
)
parser$add_argument(
    "--phenotype_type", 
    default = "cc", 
    help = "'cc' or 'quant'")
parser$add_argument("--case_prop", help = "for phenotypes of type 'cc'")
parser$add_argument(
    "--phenotype_sd", 
    help = "phenotype SD when --phenotype_type = 'quant'"
)
parser$add_argument("--phenotype_na_string", default = "NA")
parser$add_argument("--N")
args <- parser$parse_args()

#---------------------------
# QTL file
#---------------------------
qtls <- fread(args$qtl_file)

if (is.integer(qtls[[args$qtl_chr_col]])) {
    qtls[, "chr" := paste0("chr", get(args$qtl_chr_col))]
} else {
    setnames(qtls, args$qtl_chr_col, "chr")
}

qtls_tab <- qtls[, .N, c("chr", args$qtl_id_col)]

append_if_not_null <- function(cmd, arg) {
    if (!is.null(args[[arg]])) {
        cmd <- c(cmd, paste0("--", arg), args[[arg]])
    }
    return(cmd)
}

#--------------------------------------------------
# subset GWAS sumstats for each QTL
#--------------------------------------------------
for (i in 1:nrow(qtls_tab)) {
    qtl <- qtls_tab[[args$qtl_id_col]][i]
    cmd_args <- c(
        "--results", args$results,
        "--maf_col", args$maf_col,
        "--pval_col", args$pval_col,
        "--snpid_col", args$snpid_col,
        "--radius", args$radius,
        "--qtl_format", args$qtl_format,
        "--snp_format", args$snp_format,
        "--snp_file", args$snp_file,
        "--chr", qtls_tab$chr[i],
        "--qtl", qtl,
        "--out", paste0(args$out_prefix, "_", qtl, ".Rds")
    )
    
    cmd_args <- append_if_not_null(cmd_args, "beta_col")
    cmd_args <- append_if_not_null(cmd_args, "betavar_col")
    cmd_args <- append_if_not_null(cmd_args, "phenotype_file")
    cmd_args <- append_if_not_null(cmd_args, "phenotype_col")
    cmd_args <- append_if_not_null(cmd_args, "phenotype_type")
    cmd_args <- append_if_not_null(cmd_args, "case_prop")
    cmd_args <- append_if_not_null(cmd_args, "phenotype_sd")
    cmd_args <- append_if_not_null(cmd_args, "phenotype_na_string")
    cmd_args <- append_if_not_null(cmd_args, "N")
    cmd_args <- append_if_not_null(cmd_args, "betasd_col")
    
    system(
        paste(
            "Rscript --vanilla src/common/gwas_sumstats.R", 
            paste(cmd_args, collapse = " "))
    )
}
