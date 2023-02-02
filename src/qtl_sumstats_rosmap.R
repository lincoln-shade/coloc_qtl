#===================================================
# Subset GTEx QTL file for colocalization analysis
#===================================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--pval_col", default = "p")
parser$add_argument("--snpid_col", default = "rsID")
parser$add_argument("--beta_col", default = "beta")
parser$add_argument("--beta_sd_col", default = "se")
parser$add_argument("--phenotype_col")
parser$add_argument("--pos_col", default = "pos")
parser$add_argument("--radius", default = "2e5")
parser$add_argument("--chr", help="chromosome in chr[CHR] format")
parser$add_argument("--qtl")
parser$add_argument("--tissue", default = "Brain_Prefrontal_Cortex")
parser$add_argument("--phenotype")
parser$add_argument("--qtl_type")
parser$add_argument(
    "--eqtl_dir", 
    default = "/data_global/ROSMAP/xQTL_updated_data/eQTLs/"
)
parser$add_argument(
    "--mqtl_dir", 
    default = "/data_global/ROSMAP/xQTL_updated_data/mQTLs/"
)
parser$add_argument("--out", "-o")
args <- parser$parse_args()

#-------------------------
# read in qtl results file
#-------------------------
if (args$qtl_type == "eQTL") {
    qtl_dir <- args$eqtl_dir
    qtl_file <- paste0("eQTL", args$chr, "_1Mb.csv")
    N <- 534L # sample size taken from https://mostafavilab.stat.ubc.ca/xqtl/
} else if (args$qtl_type == "mQTL") {
    qtl_dir <- args$mqtl_dir
    qtl_file <- paste0("mQTL", args$chr, "_50Kb.csv")
    N <- 543L # sample size taken from https://mostafavilab.stat.ubc.ca/xqtl/ 
} else {
    stop("--qtl_type must be 'eQTL' or 'mQTL'")
}

results <- fread(paste0(qtl_dir, qtl_file))

#----------------------------
# subset qtl results
#----------------------------
results <- results[get(args$phenotype_col) == args$phenotype]

qtl_pos <- results[
    get(args$snpid_col) == args$qtl, 
    get(args$pos_col)
]

radius <- as.integer(args$radius)

results <- results[
    (get(args$pos_col) >= qtl_pos - radius) &
        (get(args$pos_col) <= qtl_pos + radius)
]

#-----------------------------
# prepare output for coloc
#-----------------------------

# beta variance from beta SD
results[, varbeta := get(args$beta_sd_col)^2]

# rename columns for coloc
setnames(
    results,
    c(args$pval_col, args$beta_col, args$snpid_col),
    c("pvalues", "beta", "snp")
)

results <- results[, .(snp, beta, varbeta, pvalues)]

# write output
if (!is.null(args$out)) {
    if (grepl("[.]Rds$", args$out, ignore.case = TRUE)) {
        out <- args$out
    } else {
        out <- paste0(args$out, ".Rds")
    }
} else {
    out <- paste0(
        "./ROSMAP_", args$tissue, "_", args$phenotype, "_", args$qtl, ".Rds"
    )
}

output <- list(
    dataset = results,
    metadata = list(
        data_source = "ROSMAP",
        qtl_type = args$qtl_type,
        variant_id_format = "rsID",
        tissue = args$tissue,
        phenotype = args$phenotype,
        lead_gwas_qtl = args$qtl,
        type = "quant",
        N = N
    )
)

saveRDS(output, file = args$out)