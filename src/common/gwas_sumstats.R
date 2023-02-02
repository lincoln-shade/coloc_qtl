#!/usr/bin/env Rscript --vanilla

#=======================================
# subset gwas data for coloc analysis
#=======================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--maf_col", default = "MAF")
parser$add_argument("--pval_col", default = "pval.spa")
parser$add_argument("--snpid_col", default = "SNPID")
parser$add_argument("--beta_col")
parser$add_argument("--betavar_col")
parser$add_argument("--betasd_col")
parser$add_argument("--radius", default = "2e5")
parser$add_argument("-r", "--results", help="path to GWAS results file")
parser$add_argument("--chr", help="chromosome in chr[CHR] format")
parser$add_argument("--qtl")
parser$add_argument(
    "--qtl_format", 
    help = "rsID, chr_bp_ref_alt, or chr_bp_ref_alt_b38 allowed"
)
parser$add_argument(
    "--snp_format", 
    default = "rsID", 
    help = "rsID, chr_bp_ref_alt, or chr_bp_ref_alt_b38 allowed")
parser$add_argument("--out", "-o")
parser$add_argument(
    "--snp_file", "-s", 
    default = "/home/lmsh224/projects/coloc/data/snp_list.txt",
    help = "SNP file that's output of src/common/write_dbsbp_common_list.R"
)

# arguments for determining features of phenotype for colocalization
parser$add_argument("--phenotype_file")
parser$add_argument(
    "--phenotype_col", 
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

#----------------------
# SNP data
#----------------------
read_snp_cmd <- function(snp_file, chr) {
    paste0("cat ", args$snp_file, " | awk -F \",\" 'NR == 1 || $1 == \"", 
           chr, "\"'")
}

snps <- fread(cmd = read_snp_cmd(args$snp_file, args$chr))

#-----------------------
# phenotype meta data
#-----------------------

if (!is.null(args$phenotype_file)) {
    if (is.null(args$phenotype_col)) {
        stop("need --phenotype_col if --phenotype_file given")
    } else {
        pheno <- fread(
            args$phenotype_file, 
            na.strings = args$phenotype_na_string
        )
        pheno <- pheno[[args$phenotype_col]]
    }
}

if (args$phenotype_type == "cc") {
    if (exists("pheno")) {
        pheno <- as.logical(pheno[!is.na(pheno)])
        N <- length(pheno)
        case_prop <- sum(pheno) / N
    } else {
        N <- as.integer(args$N)
        case_prop <- as.numeric(args$case_prop)
    }
} else if (args$phenotype_type == "quant") {
    if (exists("pheno")) {
        pheno <- pheno[!is.na(pheno)]
        N <- length(pheno)
        sdY <- sd(pheno)
    } else {
        N <- as.integer(args$N)
        if (!is.null(args$sdy)) {
            sdY <- as.numeric(args$sdy)
        } else {
            print("sdY will be estimated from beta, varbeta, MAF, and N")
            if (
                is.null(N) | is.null(args$beta_col) | 
                (is.null(args$betavar_col) & is.null(args$betasd_col)) | 
                is.null(args$maf_col)) {
                stop("Cannot estimate sdY. beta, varbeta, MAF, or N missing")
            }
        }
    }
} else {
    stop("--phenotype_type should be either 'cc' or 'quant'")
}

if (!is.null(args$phenotype_col)) {
    phenotype <- args$phenotype_col
} else {
    phenotype <- "GWAS"
}

if (!exists("sdY")) {sdY <- NULL}
if (!exists("case_prop")) {case_prop <- NULL}

metadata <- list(
    phenotype = phenotype,
    N = N,
    type = args$phenotype_type,
    sdY = sdY,
    s = case_prop,
    qtl = args$qtl
)


#-----------------------
# subset GWAS sumstats
#-----------------------
gwas <- fread(args$results)
gwas_keep_cols <- c(
    args$snpid_col, args$maf_col, args$pval_col, # necessary 
    args$beta_col, args$betavar_col, args$betasd_col # optional in some cases
)

keep_col_index <- !c(
    is.null(args$snpid_col), 
    is.null(args$maf_col),
    is.null(args$pval_col),
    is.null(args$beta_col),
    is.null(args$betavar_col),
    is.null(args$betasd_col)
)

coloc_cols <- c(args$snpid_col, "MAF", "pvalues", "beta", "varbeta", "sdbeta")
coloc_keep_cols <- coloc_cols[keep_col_index]
gwas <- gwas[, ..gwas_keep_cols]
setnames(gwas, gwas_keep_cols, coloc_keep_cols)

# calculate varbeta from sdbeta, if present
if (!is.null(gwas[["sdbeta"]])) {
    gwas[, varbeta := sdbeta^2]
    gwas[, sdbeta := NULL]
}

# merge with snp_list and subset
gwas <- merge(snps, gwas, by.x = args$snp_format, by.y = args$snpid_col)

qtl_bp <- snps[get(args$qtl_format) == args$qtl, bp]
radius <- as.integer(args$radius)
gwas <- gwas[(bp >= qtl_bp - radius) & (bp <= qtl_bp + radius)]

# write results
if (!is.null(args$out)) {
    if (grepl("[.]Rds$", args$out, ignore.case = TRUE)) {
        out <- args$out
    } else {
        out <- paste0(args$out, ".Rds")
    }
} else {
    out <- paste0("./", phenotype, "_", args$qtl, ".Rds")
}

saveRDS(list(dataset = gwas, metadata = metadata), file = out)
