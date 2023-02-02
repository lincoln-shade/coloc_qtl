#! /usr/bin/Rscript --vanilla
library(data.table)
library(argparse)
library(magrittr)

parser <- ArgumentParser()
# properties of QTL analysis summary stats files
parser$add_argument("--maf_col", default = "af")
parser$add_argument("--pval_col", default = "pval_nominal")
parser$add_argument("--snpid_col", default = "variant_id")
parser$add_argument("--beta_col", default = "slope")
parser$add_argument("--beta_sd_col", default = "slope_se")
parser$add_argument("--phenotype_col", default = "phenotype_id")
parser$add_argument("--tss_distance_col", default = "tss_distance")
parser$add_argument("--ma_count_col", default = "ma_count")
parser$add_argument("--maf_min", default = "001")


# properties of file with top QTLs for each phenotype/tissue combo
parser$add_argument("--qtl_file", "-q")
parser$add_argument("--qtl_col", default = "variant_id")
parser$add_argument("--chr_col", default = "chr")
parser$add_argument("--tissue_col", default = "tissue")
parser$add_argument("--qtl_type_col", default = "qtl_type")

parser$add_argument(
    "--radius", 
    default = "2e5", 
    help = "radius (basepairs) around top QTL to keep in subset, default = 200k"
)
parser$add_argument(
    "--eqtl_dir", 
    default = paste0(
        "/data_global/TOPMed/Omics_qtl_summaries/rnaseq/freeze-1RNA/cis-eqtl/",
        "nominal/"
    )
)
parser$add_argument("--out_prefix", "-o", default = "tmp/TOPMed/")
args <- parser$parse_args()

#---------------------------
# QTL file
#---------------------------
qtls <- fread(args$qtl_file)

if (is.integer(qtls[[args$chr_col]])) {
    qtls[, "chr" := paste0("chr", get(args$chr_col))]
} else {
    setnames(qtls, args$chr_col, "chr")
}

# loop
append_if_not_null <- function(cmd, arg) {
    if (!is.null(args[[arg]])) {
        cmd <- c(cmd, paste0("--", arg), args[[arg]])
    }
    return(cmd)
}

script_args <- names(args)[
    !(names(args) %in% c(
        "qtl_file", "qtl_col", "chr_col", "tissue_col", 
        "qtl_type_col", "out_prefix"
    ))
]

cmd_args_base <- vector(mode = "character")
for (arg in script_args) {
    cmd_args_base <- append_if_not_null(cmd_args_base, arg)
}

if (!dir.exists(args$out_prefix)) {
    dir.create(args$out_prefix, recursive = TRUE)
}

for (i in 1:nrow(qtls)) {
    qtli <- qtls[i]
    qtl <- qtli[[args$qtl_col]]
    tissue <- qtli[[args$tissue_col]]
    phenotype <- qtli[[args$phenotype_col]]
    cmd_args <- c(
        cmd_args_base,
        "--qtl", qtl,
        "--tissue", tissue,
        "--phenotype", phenotype,
        "--qtl_type", qtli[[args$qtl_type_col]],
        "--chr", qtli[["chr"]],
        "--out",
        paste0(args$out_prefix, tissue, "_", phenotype, "_", qtl, ".Rds")
    )
    
    system(
        paste(
            "Rscript --vanilla src/qtl_sumstats_topmed.R",
            paste(cmd_args, collapse = " ")
        )
    )
}