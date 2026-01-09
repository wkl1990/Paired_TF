#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="doublet removal for qc h5ad file")
#parser.add_argument("-i", "--input", type=str, dest="input", help="input h5ad file")
parser.add_argument("-m", "--meta_file", type=str, dest="meta_file", help="meta csv file")
parser.add_argument("-f", "--info_file", type=str, dest="info_file", help="info txt file")
parser.add_argument("-b", "--bam_dir", type=str, dest="bam_dir", help="bam path")
parser.add_argument("-t", "--n_cpu", type=int, default=16, dest="n_cpu", help="cpu number")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

from multiprocessing import Pool
import os
import pandas as pd
import numpy as np
import pysam
import re

# Input files
meta_file = args.meta_file
info_file = args.info_file
bam_dir = args.bam_dir
n_cpu = args.n_cpu
out_prefix = args.output


meta = pd.read_csv(meta_file, sep = ",")

meta.head()

CTCF_info =  pd.read_csv(info_file, sep = "\t")

CTCF_info.head()

CTCF_info["path"] = bam_dir + CTCF_info["RNA"] + "_mm10_sorted_rmdup.bam"

print(all(sample in pd.unique(CTCF_info["RNA"]) for sample in meta["sampleID"]))

meta.head()


NCPU = n_cpu

def run():
    """Entry point of the program"""
    print("Filtering out bam files")
    generate_bams(CTCF_info, meta, out_prefix)

def generate_bams(CTCF_info, meta, outPrefix):
    """Generate separate bam files based on cell types and samples"""
    p = Pool(NCPU)
    for idx in np.unique(meta['celltype'].astype(str)):
        for sample in np.unique(CTCF_info['RNA'].astype(str)):
            bamf = CTCF_info.loc[CTCF_info["RNA"] == sample, "path"].values[0]
            #bamf = "/" + bamf + "/outs/gex_possorted_bam.bam"
            p.apply_async(generate_bam_worker, (bamf, meta, idx, sample, outPrefix))
            #p.apply_async(generate_bam_worker, (bamf, meta, idx, sample, outPrefix), callback=func.call_back, error_callback=func.err_call_back)
    p.close()
    p.join()

def generate_bam_worker(bamf, meta, cluster, sample, prefix):
    """Worker function to generate filtered bam files"""
    print(cluster)
    cluster_name = cluster.replace(" ", "_")
    print(sample)
    # Check if the output BAM file already exists
    bam_fname = f"{prefix}metacell_{cluster_name}.{sample}.bam"
    if os.path.exists(bam_fname):
        print("For metaCell =", cluster, "Sample =", sample, "BAM file already exists. Skipping.")
        bamF.close()
        return
    name = bamf
    bamF = pysam.AlignmentFile(name)
    qnames =  list(meta[(meta['celltype'].astype(str) == cluster) & (meta['sampleID'] == sample)]['barcode'].astype(str))
    qnames_set = set(qnames)
    print(len(qnames_set))
    if len(qnames_set) > 0:
        bam_fname = f"{prefix}metacell_{cluster_name}.{sample}.bam"
        print("For metaCell =", cluster, "The filtered bam is writing to:", bam_fname)
        obam = pysam.AlignmentFile(bam_fname, "wb", template=bamF)
        for b in bamF.fetch(until_eof=True):
            tag = re.search(':\w\w:\w\w:\w\w:', b.query_name).group().strip(":")
            if tag in qnames_set:
                obam.write(b)
        obam.close()
        bamF.close()
        print("metaCell =", cluster, "Sample =", sample, "writing finished.")

if __name__ == "__main__":
    run()


