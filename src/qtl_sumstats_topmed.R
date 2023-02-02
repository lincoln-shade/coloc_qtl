#! /usr/bin/Rscript --vanilla

#=====================================================
# Subset TOPMed QTL file for colocalization analysis
#=====================================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--maf_col", default = "af")
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
parser$add_argument("--qtl_type", default = "eQTL")
parser$add_argument(
    "--eqtl_dir", 
    default = paste0(
        "/data_global/TOPMed/Omics_qtl_summaries/rnaseq/freeze-1RNA/cis-eqtl/",
        "nominal/"
    )
)
parser$add_argument("--maf_min", default = "001")
parser$add_argument("--out", "-o")
args <- parser$parse_args()

#-------------------------
# read in parquet file
#-------------------------
if (args$qtl_type == "eQTL") {
    qtl_dir <- args$eqtl_dir
    qtl_file <- paste0(args$tissue, ".maf", args$maf_min, ".tsv.gz")
} else {
    stop("--qtl_type must be 'eQTL'")
}

results <- fread(paste0(qtl_dir, qtl_file))
results <- results[`#chrom` == args$chr]

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
results[get(args$maf_col) > 0.5, (args$maf_col) := 1 - get(args$maf_col)]
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
    if (grepl("\\.Rds$", args$out, ignore.case = TRUE)) {
        out <- args$out
    } else {
        out <- paste0(args$out, ".Rds")
    }
} else {
    out <- paste0(
        "./TOPMed_", args$tissue, "_", args$phenotype, "_", args$qtl, ".Rds"
    )
}

output <- list(
    dataset = results,
    metadata = list(
        data_source = "TOPMed",
        qtl_type = args$qtl_type,
        variant_id_format = "chr_bp_ref_alt",
        tissue = args$tissue,
        phenotype = args$phenotype,
        lead_gwas_qtl = args$qtl,
        type = "quant",
        N = N
    )
)

saveRDS(output, file = args$out)
