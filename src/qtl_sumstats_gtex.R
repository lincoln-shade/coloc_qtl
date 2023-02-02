#! /usr/bin/Rscript --vanilla

#===================================================
# Subset GTEx QTL file for colocalization analysis
#===================================================

library(data.table)
library(argparse)
library(arrow)

parser <- ArgumentParser()
parser$add_argument("--maf_col", default = "maf")
parser$add_argument("--pval_col", default = "pval_nominal")
parser$add_argument("--snpid_col", default = "variant_id")
parser$add_argument("--beta_col", default = "slope")
parser$add_argument("--beta_sd_col", default = "slope_se")
parser$add_argument("--phenotype_col", default = "phenotype_id")
parser$add_argument("--tss_distance_col", default = "tss_distance")
parser$add_argument("--ma_count_col", default = "ma_count")
parser$add_argument("--radius", default = "2e5")
parser$add_argument("--chr", help="chromosome in chr[CHR] format")
parser$add_argument("--qtl")
parser$add_argument("--tissue")
parser$add_argument("--phenotype")
parser$add_argument("--qtl_type")
parser$add_argument(
    "--eqtl_dir", 
    default = paste0(
        "/data_global/GTEx/GTEx_Analysis_v8_QTLs/",
        "GTEx_Analysis_v8_EUR_eQTL_all_associations/"
    )
)
parser$add_argument(
    "--sqtl_dir", 
    default = paste0(
        "/data_global/GTEx/GTEx_Analysis_v8_QTLs/",
        "GTEx_Analysis_v8_EUR_sQTL_all_associations/"
    )
)
parser$add_argument("--out", "-o")
args <- parser$parse_args()

#-------------------------
# read in parquet file
#-------------------------
if (args$qtl_type == "eQTL") {
    qtl_dir <- args$eqtl_dir
    qtl_file <- paste0(args$tissue, ".v8.EUR.allpairs.", args$chr, ".parquet")
} else if (args$qtl_type == "sQTL") {
    qtl_dir <- args$sqtl_dir
    qtl_file <- paste0(
        args$tissue, ".v8.EUR.sqtl_allpairs.", args$chr, ".parquet")
} else {
    stop("--qtl_type must be 'eQTL' or 'sQTL'")
}

results <- setDT(read_parquet(paste0(qtl_dir, qtl_file)))

#----------------------------
# subset qtl results
#----------------------------
results <- results[get(args$phenotype_col) == args$phenotype]

qtl_tss_distance <- results[
    get(args$snpid_col) == args$qtl, 
    get(args$tss_distance_col)
]

radius <- as.integer(args$radius)

results <- results[
    (get(args$tss_distance_col) >= qtl_tss_distance - radius) &
        (get(args$tss_distance_col) <= qtl_tss_distance + radius)
]

#-----------------------------
# prepare output for coloc
#-----------------------------

# sample size estimated from MAF and MA count
N <- results[, round(mean(
    (get(args$ma_count_col) / get(args$maf_col)), 
    na.rm = TRUE
    ) / 2)
]

# beta variance from beta SD
results[, varbeta := get(args$beta_sd_col)^2]

# rename columns for coloc
setnames(
    results,
    c(args$pval_col, args$maf_col, args$beta_col, args$snpid_col),
    c("pvalues", "MAF", "beta", "snp")
)

results <- results[, .(snp, MAF, beta, varbeta, pvalues)]

# write output
if (!is.null(args$out)) {
    if (grepl("[.]Rds$", args$out, ignore.case = TRUE)) {
        out <- args$out
    } else {
        out <- paste0(args$out, ".Rds")
    }
} else {
    out <- paste0(
        "./GTEx_", args$tissue, "_", args$phenotype, "_", args$qtl, ".Rds"
    )
}

output <- list(
    dataset = results,
    metadata = list(
        data_source = "GTEx",
        qtl_type = args$qtl_type,
        variant_id_format = "chr_bp_ref_alt_b38",
        tissue = args$tissue,
        phenotype = args$phenotype,
        lead_gwas_qtl = args$qtl,
        type = "quant",
        N = N
    )
)

saveRDS(output, file = args$out)
