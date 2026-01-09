#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="subset h5ads file based on obs pattern")
parser.add_argument("-i", "--input", type=str, dest="input", help="input h5ads file")
parser.add_argument("-m", "--meta", type=str, dest="meta", default=None, help="input meta file")
parser.add_argument("-p", "--pattern", type=str, dest="pattern", help="pattern to subset")
parser.add_argument("-v", "--var", type=str, default="obs_names", dest="var", help="which var to search")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ #'2.5.3'
from pathlib import Path
import numpy as np
import pandas as pd
import sys
import re

# Input files
data_path = args.input
meta_path = args.meta
pattern = args.pattern
var = args.var
out_prefix = args.output
file = Path(data_path)

# load data
#h5ad_filename = out_prefix + ".h5ad"

data = snap.read_dataset(file, mode='r')

if var == "obs_names":
	obs = data.obs_names
elif var in data.obs:
	obs = data.obs[var]
elif meta_path is not None:
	df = pd.read_csv(meta_path)
	if var in df.columns:
		#obs = [df.loc[df['barcode']==cell, var].to_list()[0] for cell in data.obs_names]
		exist_barcode = set(df['barcode'].to_list())
		obs = [df.loc[df['barcode']==cell, var].to_list()[0] if cell in exist_barcode else 'NA' for cell in data.obs_names]
	else:
		sys.exit("Var is not support!")
else:
	sys.exit("Var is not support!")

obs_idx = [i for i, item in enumerate(obs) if re.search(pattern, item)]
#data_obs_subset = np.take(data.obs_names, obs_idx)
#data_subset = data.subset(obs_idx, out=h5ad_filename)
data_subset = data.subset(obs_idx, out=out_prefix)

data.close()

ATAC_cells = np.array(data_subset[0].obs_names)
cells_filename = out_prefix + ".cells.txt"
np.savetxt(cells_filename, ATAC_cells, fmt="%s", delimiter=',')
data_subset[0].close()

print("Job is done!")


