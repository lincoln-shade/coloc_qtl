#!/usr/bin/python3

import sys
import argparse
import multiprocessing as mp
import os
from itertools import chain
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "infile", 
    nargs = '?', 
    default="./snp_query.csv"
)
parser.add_argument(
    "--eqtl_dir", 
    default=(
        "/data_global/GTEx/GTEx_Analysis_v8_QTLs/"
        "GTEx_Analysis_v8_eQTL_Eur/eqtls/"
    )
)
parser.add_argument(
    "--sqtl_dir", 
    default=(
        "/data_global/GTEx/GTEx_Analysis_v8_QTLs/"
        "GTEx_Analysis_v8_sQTL_Eur/GTEx_Analysis_v8_sQTL_EUR/"
    )
)
parser.add_argument('-o', '--outfile', default="./gtex_qtl_test.txt")
parser.add_argument('-p', '--pool', type=int, default=22)
args = parser.parse_args()
# snps = [line.strip() for line in args.infile]
snps = pd.read_csv(args.infile)
snp_id = snps["chr_bp_ref_alt_b38"]

if len(snp_id) == 0:
    sys.exit("No variant names were entered.")

snp_set = set(snp_id)
eqtl_dir = args.eqtl_dir
sqtl_dir = args.sqtl_dir
eqtl_dir_files = os.listdir(eqtl_dir)
sqtl_dir_files = os.listdir(sqtl_dir)
signif_pairs = '.signif_pairs.txt'

# remove the 'xgenes' files from the list of files in the directories
def remove_genes_files(files, suffix):
    x = []
    for file in files:
        if suffix in file:
            x.append(file)
    return x

eqtl_files = remove_genes_files(files = eqtl_dir_files, suffix = signif_pairs)
sqtl_files = remove_genes_files(files = sqtl_dir_files, suffix = signif_pairs)

def find_qtl_in_file(i, file, snp_set, path):
    out = []
    file_path = path + file
    with open(file_path) as qtl_file:
        lines = qtl_file.readlines()
        for line in lines:
            # variant IDs in second column
            if line.split('\t')[1] in snp_set: 
                out.append(file + '\t' + line)
    return (i, out)

def collect_result(result):
    global results
    results.append(result)

# create pool and find QTLs in GTEx files
results = []

if args.pool:
    pool_n = args.pool
else:
    pool_n = 1

pool = mp.Pool(pool_n)

for i, file in enumerate(eqtl_files):
    pool.apply_async(find_qtl_in_file, 
                     args=(i, file, snp_set, eqtl_dir), 
                     callback=collect_result
    )
    
for i, file in enumerate(sqtl_files):
    pool.apply_async(find_qtl_in_file, 
                     args=(i, file, snp_set, sqtl_dir), 
                     callback=collect_result
    )

pool.close()

# postpones the execution of next line of code until all 
# processes in the queue are done.
pool.join()  


# sort and format results
results.sort(key=lambda x: x[0])
results_final = [r for i, r in results]
results_final = list(chain.from_iterable(results_final))

# write output
if args.outfile:
    with open(args.outfile, 'x') as out_file:
        for i in results_final:
            out_file.write(i)
else:
    for i in results_final:
        print(i, end="", file=sys.stdout)
