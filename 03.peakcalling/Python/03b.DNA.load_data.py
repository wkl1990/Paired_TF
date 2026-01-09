#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="load fragment file to snap")
parser.add_argument("-i", "--input", type=str, dest="input", help="input raw file")
parser.add_argument("-g", "--genome", type=str, default="hg38", dest="genome", help="genome")
parser.add_argument("-c", "--cells", type=str, default=None, dest="cells", help="cell barcode file")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ 
from pathlib import Path
import numpy as np
import pandas as pd
import sys

# Input files
frag_path = args.input
cells_file = args.cells
genome = args.genome
out_prefix = args.output
file = Path(frag_path)

if genome == "mm10":
	snap_genome = snap.genome.mm10
elif genome == "hg38":
	snap_genome = snap.genome.hg38
else:
	sys.exit("Genome version is not support!")

# load data
h5ad_filename = out_prefix + ".h5ad"
if cells_file is None:
	data = snap.pp.import_data(file, chrom_sizes=snap_genome, file=h5ad_filename, sorted_by_barcode=False, min_num_fragments=0)
else:
	data = snap.pp.import_data(file, chrom_sizes=snap_genome, file=h5ad_filename, sorted_by_barcode=False, min_num_fragments=0, whitelist=cells_file)

ATAC_cells = np.array(data.obs_names)
cells_filename = out_prefix + ".cells.txt"
np.savetxt(cells_filename, ATAC_cells, fmt="%s", delimiter=',')
data.close()

print("Job is done!")


