#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="subset h5ad file based on obs pattern")
parser.add_argument("-i", "--input", type=str, dest="input", help="input raw file")
parser.add_argument("-n", "--ngroup", type=str, dest="ngroup", default=None, help="inside group obs")
parser.add_argument("-g", "--group", type=str, dest="group", default=None, help="input group csv file")
parser.add_argument("-s", "--selection", type=str, dest="selection", default=None, help="selection to output")
parser.add_argument("-v", "--var", type=str, default="celltype", dest="var", help="which var to output")
parser.add_argument("-p", "--path", type=str, dest="outpath", help="output path")
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
ngroup = args.ngroup
group = args.group
selection = args.selection
var = args.var
out_path = args.outpath
out_prefix = args.output
file = Path(data_path)

# load data
merged_data = snap.read_dataset(file, mode='r')
#cluster_data = merged_data.copy(h5ads_filename)
#merged_data.close()

if group is not None:
	df = pd.read_csv(group)
	#df['barcode'] = df['barcode'].str.replace('_',':')
	exist_barcode = set(df['barcode'].to_list())
	merged_data_group = [df.loc[df['barcode']==cell, var].to_list()[0] if cell in exist_barcode else 'NA' for cell in merged_data.obs_names]
if ngroup is not None:
	if ngroup in merged_data.obs[:].columns:
		merged_data_group = merged_data.obs[ngroup].to_list()
if merged_data_group is None:
	sys.exit("Group is not NULL!")

# output
snap.ex.export_fragments(merged_data, groupby=merged_data_group, selections=selection, out_dir=out_path, prefix=out_prefix, suffix='.bed.gz')

merged_data.close()

print("Job is done!")


