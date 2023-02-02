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
    default=("/data_global/ROSMAP/xQTL_updated_data/sigeQTLs/")
)
parser.add_argument(
    "--mqtl_dir", 
    default=("/data_global/ROSMAP/xQTL_updated_data/sigmQTLs/")
)
parser.add_argument('-o', '--outfile', default="./rosmap_qtl_test.txt")
parser.add_argument('-p', '--pool', type=int, default=22)
args = parser.parse_args()
snps = pd.read_csv(args.infile)
snp_id = snps["rsID"]

if len(snp_id) == 0:
    sys.exit("No variant names were entered.")

snp_set = set(snp_id)
eqtl_dir = args.eqtl_dir
mqtl_dir = args.mqtl_dir
eqtl_files = os.listdir(eqtl_dir)
mqtl_files = os.listdir(mqtl_dir)

def find_mqtl_in_file(i, file, snp_set, path):
    out = []
    file_path = path + file
    with open(file_path) as qtl_file:
        lines = qtl_file.readlines()
        for line in lines:
            # variant IDs in second column
            if line.split(',')[0] in snp_set: 
                out.append(file + ',mQTL,' + line)
    return (i, out)

def find_eqtl_in_file(i, file, snp_set, path):
    out = []
    file_path = path + file
    with open(file_path) as qtl_file:
        lines = qtl_file.readlines()
        for line in lines:
            # variant IDs in second column
            if line.split(',')[0] in snp_set:
                line_split = line.split(',')
                del line_split[6] # remove gene symbol so columns match mQTL
                line = ','.join(line_split)
                out.append(file + ',eQTL,' + line)
    return (i, out)

def collect_result(result):
    global results
    results.append(result)

# create pool and find QTLs in files
results = []

if args.pool:
    pool_n = args.pool
else:
    pool_n = 1

pool = mp.Pool(pool_n)

for i, file in enumerate(eqtl_files):
    pool.apply_async(find_eqtl_in_file, 
                     args=(i, file, snp_set, eqtl_dir), 
                     callback=collect_result
    )
    
for i, file in enumerate(mqtl_files):
    pool.apply_async(find_mqtl_in_file, 
                     args=(i, file, snp_set, mqtl_dir), 
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
