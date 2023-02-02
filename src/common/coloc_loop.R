#! /usr/bin/Rscript --vanilla
library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-g", "--gwas_dir", default = "tmp/gwas/")
parser$add_argument("-q", "--qtl_dir", default = "tmp/GTEx/")
parser$add_argument("--qtl_format", default = "chr_bp_ref_alt_b38")
parser$add_argument("--p12", default = "1e-5")
parser$add_argument("-o", "--out_file")
args <- parser$parse_args()

gwas_files <- list.files(args$gwas_dir)
qtl_files <- list.files(args$qtl_dir)

for (gf in gwas_files) {
    path <- paste0(args$gwas_dir, gf)
    gwas <- readRDS(path)
    qtls <- qtl_files[grep(gwas$metadata$qtl, qtl_files)]
    for (qf in qtls) {
        cmd_args <- c(
            "-g", path, 
            "-q", paste0(args$qtl_dir, qf), 
            "-o", args$out_file,
            "--p12", args$p12
        )
        
        system(
            paste0(
                "Rscript --vanilla src/common/coloc.R ", 
                paste(cmd_args, collapse = " ")
            )
        )
    }
}


