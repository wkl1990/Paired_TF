#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="subset h5ad file based on obs pattern")
parser.add_argument("-i", "--input", type=str, dest="input", help="input raw file")
parser.add_argument("-p", "--pattern", type=str, dest="pattern", help="pattern to subset")
parser.add_argument("-v", "--var", type=str, default="obs_names", dest="var", help="which var to search")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ 
from pathlib import Path
import numpy as np
import pandas as pd
import sys
import re

# Input files
data_path = args.input
pattern = args.pattern
var = args.var
out_prefix = args.output
file = Path(data_path)

# load data
h5ad_filename = out_prefix + ".h5ad"

data = snap.read(file, backed='r')

if var == "obs_names":
	obs = data.obs_names
elif var == "celltype":
	obs = data.obs[var]
else:
	sys.exit("Var is not support!")

obs_idx = [i for i, item in enumerate(obs) if re.search(pattern, item)]
#data_obs_subset = np.take(data.obs_names, obs_idx)
#data_subset = data.subset(obs_idx, out=h5ad_filename)
data_subset = data.subset(obs_idx, inplace=False)
data_subset.write(h5ad_filename)

data.close()

ATAC_cells = np.array(data_subset.obs_names)
cells_filename = out_prefix + ".cells.txt"
np.savetxt(cells_filename, ATAC_cells, fmt="%s", delimiter=',')
#data_subset.close()

print("Job is done!")


