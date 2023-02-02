#! /usr/bin/python3

import sys
import argparse
import os
import subprocess

def append_if_not_none(cmd, string, x):
    if x:
        cmd = cmd + [string, x]
    
    return cmd

def run_extract_variant_ids(
        gwas_sumstats, out_file, snp_col=None, chrom_col=None, snp_format=None, 
        pval_col=None, max_pval=None):
    
    cmd = ["Rscript", "--vanilla", "./src/common/extract_variant_ids.R", 
           "-f", gwas_sumstats, "-o", out_file]
    cmd = append_if_not_none(cmd, "--snp_col", snp_col)
    cmd = append_if_not_none(cmd, "--chrom_col", chrom_col)
    cmd = append_if_not_none(cmd, "--snp_format", snp_format)
    cmd = append_if_not_none(cmd, "--pval_col", pval_col)
    cmd = append_if_not_none(cmd, "--max_pval", max_pval)
    subprocess.run(cmd)

def run_subset_snp_list(ref_file, query_file, snp_format, out):
    cmd = ["Rscript", "--vanilla", "./src/common/subset_snp_list.R", "-r", ref_file, "-q", query_file,
           "-o", out]
    cmd = append_if_not_none(cmd, "--snp_format", snp_format)
    subprocess.run(cmd)

def run_gwas_sumstats_loop(
        gwas_results, out_pref, qtl_file, qtl_format, qtl_id_col, qtl_chr_col, 
        snp_file, 
        radius, snp_format, maf_col, pval_col, snpid_col, beta_col, betavar_col, 
        betasd_col, phenotype_file, phenotype_type, phenotype_col, case_prop,
        phenotype_sd, phenotype_na_string, N):
    
    cmd = ["Rscript", "--vanilla", "./src/common/gwas_sumstats_loop.R", 
           "-r", gwas_results, 
           "--qtl_file", qtl_file, "--qtl_format", qtl_format, "--qtl_id_col", 
           qtl_id_col, "--qtl_chr_col", qtl_chr_col, "--snp_file", snp_file,
           "--radius", radius, "--snp_format", snp_format, "-o", out_pref]
    cmd = append_if_not_none(cmd, "--maf_col", maf_col)
    cmd = append_if_not_none(cmd, "--pval_col", pval_col)
    cmd = append_if_not_none(cmd, "--snpid_col", snpid_col)
    cmd = append_if_not_none(cmd, "--beta_col", beta_col)
    cmd = append_if_not_none(cmd, "--betavar_col", betavar_col)
    cmd = append_if_not_none(cmd, "--betasd_col", betasd_col)
    cmd = append_if_not_none(cmd, "--phenotype_file", phenotype_file)
    cmd = append_if_not_none(cmd, "--phenotype_type", phenotype_type)
    cmd = append_if_not_none(cmd, "--phenotype_col", phenotype_col)
    cmd = append_if_not_none(cmd, "--case_prop", case_prop)
    cmd = append_if_not_none(cmd, "--phenotype_sd", phenotype_sd)
    cmd = append_if_not_none(cmd, "--phenotype_na_string", phenotype_na_string)
    cmd = append_if_not_none(cmd, "--N", N)
    subprocess.run(cmd)

def run_coloc_loop(gwas_dir, qtl_dir, qtl_format, p12, outfile):
    cmd = ["Rscript", "--vanilla", "./src/common/coloc_loop.R", "-g", gwas_dir, "-q", qtl_dir, "--p12",
           p12, "-o", outfile]
    subprocess.run(cmd)

# GTEx functions
def run_find_qtls_gtex(infile, outfile, pool, eqtl_dir, sqtl_dir):
    cmd = ["python3", "./src/find_qtls_gtex.py", infile, "-o", outfile, "--eqtl_dir", 
           eqtl_dir, "--sqtl_dir", sqtl_dir, "-p", pool]
    subprocess.run(cmd)

def run_tidy_qtls_gtex(infile, outfile):
    cmd = ["Rscript", "--vanilla", "./src/tidy_qtls_gtex.R", "-f", infile, "-o", outfile]
    subprocess.run(cmd)

def run_qtl_sumstats_gtex_loop(infile, outfile, radius, eqtl_dir, sqtl_dir):
    cmd = ["Rscript", "--vanilla", "./src/qtl_sumstats_gtex_loop.R", "-q", infile, "-o", outfile,
           "--radius", radius, "--eqtl_dir", eqtl_dir, "--sqtl_dir", sqtl_dir]
    subprocess.run(cmd)

# ROSMAP functions
def run_find_qtls_rosmap(infile, outfile, pool, eqtl_dir, mqtl_dir):
    cmd = ["python3", "./src/find_qtls_rosmap.py", infile, "-o", outfile, "--eqtl_dir", 
           eqtl_dir, "--mqtl_dir", mqtl_dir, "-p", pool]
    subprocess.run(cmd)

def run_tidy_qtls_rosmap(infile, outfile):
    cmd = ["Rscript", "--vanilla", "./src/tidy_qtls_rosmap.R", "-f", infile, "-o", outfile]
    subprocess.run(cmd)

def run_qtl_sumstats_rosmap_loop(infile, outfile, radius, eqtl_dir, mqtl_dir):
    cmd = ["Rscript", "--vanilla", "./src/qtl_sumstats_rosmap_loop.R", "-q", infile, "-o", outfile,
           "--radius", radius, "--eqtl_dir", eqtl_dir, "--mqtl_dir", mqtl_dir]
    subprocess.run(cmd)

# TOPMed functions
def run_find_qtls_topmed(infile, outfile, pool, eqtl_dir):
    cmd = ["python3", "./src/find_qtls_topmed.py", infile, "-o", outfile, "--eqtl_dir", 
           eqtl_dir, "-p", pool]
    subprocess.run(cmd)

def run_tidy_qtls_topmed(infile, outfile):
    cmd = ["Rscript", "--vanilla", "./src/tidy_qtls_topmed.R", "-f", infile, "-o", outfile]
    subprocess.run(cmd)

def run_qtl_sumstats_topmed_loop(infile, outfile, radius, eqtl_dir):
    cmd = ["Rscript", "--vanilla", "./src/qtl_sumstats_topmed_loop.R", "-q", infile, "-o", outfile,
           "--radius", radius, "--eqtl_dir", eqtl_dir]
    subprocess.run(cmd)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action = "store_true", help = "run colocalization using GTEx, ROSMAP, and TOPMed QTLs")
    parser.add_argument("--gtex", action = "store_true", help = "run colocalization using GTEx QTLs")
    parser.add_argument("--rosmap", action = "store_true", help = "run colocalization using ROSMAP QTLs")
    parser.add_argument("--topmed", action = "store_true", help = "run colocalization using TOPMed QTLs")
    parser.add_argument("--gwas_sumstats", "-g", help = "gwas sumstats file")
    parser.add_argument("--out_pref", "-o", default = "./tmp/coloc", help = "directory for coloc_qtl output")
    # extract_variant_ids args
    parser.add_argument("--snp_col", default = "SNPID", help = "SNP column name in --gwas_sumstats")
    parser.add_argument("--snp_format", default = "rsID", help=("format that the SNP IDs are in in --gwas_sumstats. Default is rsID. Options are rsID, chr_bp_ref_alt, or chr_bp_ref_alt_b38."))
    parser.add_argument("--chrom_col", default = "chr", help=("chromosome column in --gwas_sumstats. Values should be formatted as \"[chr]\" or \"chr[chr]\""))
    parser.add_argument("--pval_col", help = "p-value column in --gwas_sumstats")
    parser.add_argument("--max_pval", help = "optional maximum pval threshold for GWAS variants")
    # subset_snp_list args
    # parser.add_argument("--query_file", "-q")
    parser.add_argument("--ref_file", "-r", default = "data/snp_list.txt", help = "file with SNP info output from ./src/common/write_dbsnp_common_list.R")
    # parser.add_argument("--out", "-o", default = "./snp_query.csv")
    # find_qtls_gtex args
    parser.add_argument("--eqtl_dir_gtex", default=("/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_Eur/eqtls/"), help = "directory with GTEx significant eQTLs")
    parser.add_argument("--sqtl_dir_gtex", default=("/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_sQTL_Eur/GTEx_Analysis_v8_sQTL_EUR/"), help = "directory with GTEx significant sQTLs")
    parser.add_argument('-p', '--pool', default="1", help = "number of threads to use for looking up significant QTLs")
    # qtl_sumstats_gtex_loop args
    parser.add_argument("--radius", default = "2e5", help = "radius (basepairs) around top QTL to keep in subset, default = 200k")
    parser.add_argument("--eqtl_all_dir_gtex", default = "/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/", help = "directory with GTEx all eQTL results")
    parser.add_argument("--sqtl_all_dir_gtex", default = "/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_sQTL_all_associations/", help = "directory with GTEx all sQTL results")
    # find_qtls_rosmap args
    parser.add_argument("--eqtl_dir_rosmap", default = "/data_global/ROSMAP/xQTL_updated_data/sigeQTLs/", help = "directory with ROSMAP significant eQTLs")
    parser.add_argument("--mqtl_dir_rosmap", default = "/data_global/ROSMAP/xQTL_updated_data/sigmQTLs/", help = "directory with ROSMAP significant mQTLs")
    # qtl_sumstats_rosmap args
    parser.add_argument("--eqtl_all_rosmap_dir", default  = "/data_global/ROSMAP/xQTL_updated_data/eQTLs/", help = "directory with ROSMAP all eQTL results")
    parser.add_argument("--mqtl_all_rosmap_dir", default  = "/data_global/ROSMAP/xQTL_updated_data/mQTLs/", help = "directory with ROSMAP all mQTLs results")
    # find_qtls_topmed args
    parser.add_argument("--eqtl_dir_topmed", default = "./data/topmed_qtl/", help = "directory with TOPMed significant eQTLs")
    # qtls_sumstats_topmed_loop args
    parser.add_argument("--eqtl_all_dir_topmed",default = "/data_global/TOPMed/Omics_qtl_summaries/rnaseq/freeze-1RNA/cis-eqtl/nominal/")
    # gwas_sumstats_loop args
    parser.add_argument("--maf_col", help = "MAF column in --gwas_sumstats")
    parser.add_argument("--beta_col", help = "Beta column in --gwas_sumstats")
    parser.add_argument("--betavar_col", help = "Beta variance column in --gwas_sumstats")
    parser.add_argument("--betasd_col", help = "Beta SE column in --gwas_sumstats")
    parser.add_argument("--phenotype_file", help = "text file with individual phenotype data")
    parser.add_argument("--phenotype_type", help = "type of GWAS phenotype, 'cc' (case/control) or 'quant' (quantitative) allowed")
    parser.add_argument("--phenotype_col", help = "phenotype column (if --phenotype_file used) or phenotype name")
    parser.add_argument("--case_prop", help = "case proportion (if --phenotype_type = 'cc')")
    parser.add_argument("--phenotype_sd", help = "phenotype SD (if --phenotype_type = 'quant')")
    parser.add_argument("--phenotype_na_string", default = "NA", help = "String for missing values in --phenotype_file")
    parser.add_argument("--N", help = "GWAS sample size")
    # coloc_loop args
    parser.add_argument("--p12", default = "1e-5", help = "Prior probability a variant is associated with both GWAS and QTL phenotypes")
    args = parser.parse_args()
    
    if args.all:
        args.gtex = True
        args.rosmap = True
        args.topmed = True
    
    with open(args.gwas_sumstats, 'r') as fp:
        for count, line in enumerate(fp):
            pass
    if count > 10e4 and not args.max_pval:
        sys.exit("Warning, --gwas_sumstats argument has >10,000 variants and no --max_pval argument."
                 "\nThis will lead to every variant in the file being checked for QTL activity, even if it is not associated with the GWAS phenotype."
                 "\nIf this was intentional, re-run with --max_pval set to 1.")
    
    if not os.path.exists(os.path.split(args.out_pref)[0]):
        os.makedirs(os.path.split(args.out_pref)[0])
    
    extract_variant_ids_out = args.out_pref + "variant_list.csv"
    run_extract_variant_ids(
        args.gwas_sumstats, extract_variant_ids_out, args.snp_col,
        args.chrom_col, args.snp_format, args.pval_col, args.max_pval
    )
    
    subset_snp_list_out = args.out_pref + "snp_list.csv"
    run_subset_snp_list(args.ref_file, extract_variant_ids_out, args.snp_format,
                        subset_snp_list_out)
    
    gwas_sumstats_loop_out_dir = args.out_pref + "GWAS/"
    if not os.path.exists(gwas_sumstats_loop_out_dir):
        os.makedirs(gwas_sumstats_loop_out_dir)
    
    if args.gtex:
        find_qtls_gtex_out = args.out_pref + "qtls_gtex.txt"
        run_find_qtls_gtex(subset_snp_list_out, find_qtls_gtex_out, args.pool,
                           args.eqtl_dir_gtex, args.sqtl_dir_gtex)
        
        tidy_qtls_gtex_out = args.out_pref + "qtls_gtex_tidy.txt"
        run_tidy_qtls_gtex(find_qtls_gtex_out, tidy_qtls_gtex_out)
        
        qtl_sumstats_gtex_loop_out = args.out_pref + "GTEx/"
        if not os.path.exists(qtl_sumstats_gtex_loop_out):
            os.makedirs(qtl_sumstats_gtex_loop_out)
        
        run_qtl_sumstats_gtex_loop(
            tidy_qtls_gtex_out, qtl_sumstats_gtex_loop_out,
            args.radius, args.eqtl_all_dir_gtex, args.sqtl_all_dir_gtex
        )
        
        run_gwas_sumstats_loop(
            args.gwas_sumstats, gwas_sumstats_loop_out_dir + "gwas",
            tidy_qtls_gtex_out, "chr_bp_ref_alt_b38",
            "variant_id", "chr", args.ref_file, args.radius, args.snp_format,
            args.maf_col, args.pval_col, args.snp_col, args.beta_col,
            args.betavar_col, args.betasd_col, args.phenotype_file,
            args.phenotype_type, args.phenotype_col, args.case_prop,
            args.phenotype_sd, args.phenotype_na_string, args.N
        )
        
        coloc_loop_out = args.out_pref + "coloc_results.txt"
        run_coloc_loop(gwas_sumstats_loop_out_dir, qtl_sumstats_gtex_loop_out,
                       "chr_bp_ref_alt_b38", args.p12, coloc_loop_out)
    
    if args.rosmap:
        find_qtls_rosmap_out = args.out_pref + "qtls_rosmap.txt"
        run_find_qtls_rosmap(subset_snp_list_out, find_qtls_rosmap_out,
                             args.pool, args.eqtl_dir_rosmap,
                             args.mqtl_dir_rosmap)
        
        tidy_qtls_rosmap_out = args.out_pref + "qtls_rosmap_tidy.txt"
        run_tidy_qtls_rosmap(find_qtls_rosmap_out, tidy_qtls_rosmap_out)
        
        qtl_sumstats_rosmap_loop_out = args.out_pref + "ROSMAP/"
        if not os.path.exists(qtl_sumstats_rosmap_loop_out):
            os.makedirs(qtl_sumstats_rosmap_loop_out)
        
        run_qtl_sumstats_rosmap_loop(
            tidy_qtls_rosmap_out, qtl_sumstats_rosmap_loop_out,
            args.radius, args.eqtl_all_rosmap_dir, args.mqtl_all_rosmap_dir
        )
        
        run_gwas_sumstats_loop(
            args.gwas_sumstats, gwas_sumstats_loop_out_dir + "gwas",
            tidy_qtls_rosmap_out, "rsID",
            "rsID", "chr", args.ref_file, args.radius, args.snp_format,
            args.maf_col, args.pval_col, args.snp_col, args.beta_col,
            args.betavar_col, args.betasd_col, args.phenotype_file,
            args.phenotype_type, args.phenotype_col, args.case_prop,
            args.phenotype_sd, args.phenotype_na_string, args.N
        )
        
        coloc_loop_out = args.out_pref + "coloc_results.txt"
        run_coloc_loop(gwas_sumstats_loop_out_dir, qtl_sumstats_rosmap_loop_out,
                       "rsID", args.p12, coloc_loop_out)
    
    if args.topmed:
        find_qtls_topmed_out = args.out_pref + "qtls_topmed.txt"
        run_find_qtls_topmed(subset_snp_list_out, find_qtls_topmed_out,
                             args.pool, args.eqtl_dir_topmed)
        
        tidy_qtls_topmed_out = args.out_pref + "qtls_topmed_tidy.txt"
        run_tidy_qtls_topmed(find_qtls_topmed_out, tidy_qtls_topmed_out)
        
        qtl_sumstats_topmed_loop_out = args.out_pref + "TOPMed/"
        if not os.path.exists(qtl_sumstats_topmed_loop_out):
            os.makedirs(qtl_sumstats_topmed_loop_out)
        
        run_qtl_sumstats_topmed_loop(
            tidy_qtls_topmed_out, qtl_sumstats_topmed_loop_out,
            args.radius, args.eqtl_all_dir_topmed
        )
        
        run_gwas_sumstats_loop(
            args.gwas_sumstats, gwas_sumstats_loop_out_dir + "gwas",
            tidy_qtls_topmed_out, "chr_bp_ref_alt",
            "variant_id", "chr", args.ref_file, args.radius, args.snp_format,
            args.maf_col, args.pval_col, args.snp_col, args.beta_col,
            args.betavar_col, args.betasd_col, args.phenotype_file,
            args.phenotype_type, args.phenotype_col, args.case_prop,
            args.phenotype_sd, args.phenotype_na_string, args.N
        )
        
        coloc_loop_out = args.out_pref + "coloc_results.txt"
        run_coloc_loop(gwas_sumstats_loop_out_dir, qtl_sumstats_topmed_loop_out,
                       "chr_bp_ref_alt", args.p12, coloc_loop_out)
    
    print("\nDone!")
