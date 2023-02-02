#! /usr/bin/Rscript --vanilla
library(data.table)
library(argparse)
library(magrittr)

parser <- ArgumentParser()
# properties of QTL analysis summary stats files
parser$add_argument("--pval_col", default = "p")
parser$add_argument("--snpid_col", default = "rsID")
parser$add_argument("--beta_col", default = "beta")
parser$add_argument("--beta_sd_col", default = "se")
parser$add_argument("--pos_col", default = "pos")
# phenotypes have different col names for mQTL and eQTL files
parser$add_argument("--eqtl_phenotype_col", default = "ENSG")
parser$add_argument("--mqtl_phenotype_col", default = "CpG")


# properties of file with top QTLs for each phenotype/tissue combo
parser$add_argument("--qtl_file", "-q")
parser$add_argument("--qtl_col", default = "rsID")
parser$add_argument("--chr_col", default = "chr")
parser$add_argument("--tissue_col", default = "tissue")
parser$add_argument("--qtl_type_col", default = "qtl_type")
parser$add_argument("--phenotype_col", default = "phenotype_id")

parser$add_argument(
    "--radius", 
    default = "2e5", 
    help = "radius (basepairs) around top QTL to keep in subset, default = 200k"
)
parser$add_argument(
    "--eqtl_dir", 
    default = "/data_global/ROSMAP/xQTL_updated_data/eQTLs/"
)
parser$add_argument(
    "--mqtl_dir", 
    default = "/data_global/ROSMAP/xQTL_updated_data/mQTLs/"
)
parser$add_argument("--out_prefix", "-o", default = "tmp/ROSMAP/")
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

qtls[qtl_type == "eQTL", phenotype_col := args$eqtl_phenotype_col]
qtls[qtl_type == "mQTL", phenotype_col := args$mqtl_phenotype_col]

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
        "qtl_type_col", "out_prefix", "phenotype_col", "eqtl_phenotype_col",
        "mqtl_phenotype_col"
    ))
]

cmd_args_base <- vector(mode = "character")
for (arg in script_args) {
    cmd_args_base <- append_if_not_null(cmd_args_base, arg)
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
        "--phenotype_col", qtli[["phenotype_col"]],
        "--qtl_type", qtli[[args$qtl_type_col]],
        "--chr", qtli[["chr"]],
        "--out",
        paste0(args$out_prefix, tissue, "_", phenotype, "_", qtl, ".Rds")
    )
    
    system(
        paste(
            "Rscript --vanilla src/qtl_sumstats_rosmap.R",
            paste(cmd_args, collapse = " ")
        )
    )
}