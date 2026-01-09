#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="clustering for annset h5ads file")
parser.add_argument("-i", "--input", type=str, dest="input", help="input annset h5ads file")
parser.add_argument("-c", "--precelltype", type=str, dest="precelltype", default=None, help="input previous celltype csv file")
parser.add_argument("-v", "--var", type=str, default="celltype", dest="var", help="which var to use")
parser.add_argument("-m", "--mnc", help="mnc correct", action="store_true", dest="mnc")
parser.add_argument("-y", "--harmony", help="harmony correct", action="store_true", dest="harmony")
parser.add_argument("-f", "--n_features", type=int, default=500000, dest="n_features", help="n_features")
parser.add_argument("-l", "--blacklist", type=str, default=None, dest="blacklist", help="blacklist")
parser.add_argument("-w", "--whitelist", type=str, default=None, dest="whitelist", help="whitelist")
parser.add_argument("-d", "--filter_lower_quantile", type=float, default=0.005, dest="filter_lower_quantile", help="filter_lower_quantile")
parser.add_argument("-u", "--filter_upper_quantile", type=float, default=0.005, dest="filter_upper_quantile", help="filter_upper_quantile")
parser.add_argument("-n", "--n_neighbors", type=int, default=50, dest="n_neighbors", help="n_neighbors")
parser.add_argument("-r", "--resolution", type=float, default=1, dest="resolution", help="resolution")
parser.add_argument("-g", "--height", type=int, default=500, dest="height", help="plot height")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ 
from pathlib import Path
import numpy as np
import pandas as pd
import sys

# Input files
qc_path = args.input
precelltype = args.precelltype
var = args.var
mnc_correct = args.mnc
harmony = args.harmony
n_features = args.n_features
blacklist_path = args.blacklist
whitelist_path = args.whitelist
filter_lower_quantile = args.filter_lower_quantile
filter_upper_quantile = args.filter_upper_quantile
n_neighbors = args.n_neighbors
resolution = args.resolution
height = args.height
out_prefix = args.output
file = Path(qc_path)
use_rep = "X_spectral"

# qc 
h5ads_filename = out_prefix + ".h5ads"

umap_sample_filename = out_prefix + ".sample.umap.pdf"
umap_mnc_sample_filename = out_prefix + ".mnc.sample.umap.pdf"
umap_harmony_sample_filename = out_prefix + ".harmony.sample.umap.pdf"

umap_cluster_filename = out_prefix + ".cluster.umap.pdf"
umap_mnc_cluster_filename = out_prefix + ".mnc.cluster.umap.pdf"
umap_harmony_cluster_filename = out_prefix + ".harmony.cluster.umap.pdf"

umap_precelltype_filename = out_prefix + ".precelltype.umap.pdf"
umap_mnc_precelltype_filename = out_prefix + ".mnc.precelltype.umap.pdf"
umap_harmony_precelltype_filename = out_prefix + ".harmony.precelltype.umap.pdf"

merged_data = snap.read_dataset(file, mode='r+')
#cluster_data = merged_data.copy(h5ads_filename)
#merged_data.close()

if precelltype is not None:
	df = pd.read_csv(precelltype)
	#df['barcode'] = df['barcode'].str.replace('_',':')
	if var in df.columns:
		merged_data_celltype = [df.loc[df['barcode']==cell, var].to_list()[0] for cell in merged_data.obs_names]
		merged_data.obs[var] = merged_data_celltype
	else:
		sys.exit("Var is not support!")

if blacklist_path != None and whitelist_path != None:
	snap.pp.select_features(merged_data, n_features=n_features, filter_lower_quantile=filter_lower_quantile, filter_upper_quantile=filter_upper_quantile, blacklist=blacklist_path, whitelist=whitelist_path)
elif blacklist_path == None and whitelist_path != None:
	snap.pp.select_features(merged_data, n_features=n_features, filter_lower_quantile=filter_lower_quantile, filter_upper_quantile=filter_upper_quantile, whitelist=whitelist_path)
elif blacklist_path != None and whitelist_path == None:
	snap.pp.select_features(merged_data, n_features=n_features, blacklist=blacklist_path, filter_lower_quantile=filter_lower_quantile, filter_upper_quantile=filter_upper_quantile)
else:
	snap.pp.select_features(merged_data, n_features=n_features, filter_lower_quantile=filter_lower_quantile, filter_upper_quantile=filter_upper_quantile)


snap.tl.spectral(merged_data)
snap.tl.umap(merged_data, use_rep=use_rep)
snap.pl.umap(merged_data, color="sample", show=False, out_file=umap_sample_filename, height=height)

snap.pp.knn(merged_data, n_neighbors=n_neighbors, use_rep=use_rep)
snap.tl.leiden(merged_data, resolution=resolution)
snap.pl.umap(merged_data, color='leiden', show=False, out_file=umap_cluster_filename, height=height)
if precelltype is not None:
	snap.pl.umap(merged_data, color=var, show=False, out_file=umap_precelltype_filename, height=height)


if mnc_correct:
	snap.pp.mnc_correct(merged_data, batch="sample")
	#use_rep="X_spectral_mnn"
	snap.tl.umap(merged_data, use_rep="X_spectral_mnn")
	snap.pl.umap(merged_data, color="sample", show=False, out_file=umap_mnc_sample_filename, height=height)
	snap.pp.knn(merged_data, n_neighbors=n_neighbors, use_rep="X_spectral_mnn")
	snap.tl.leiden(merged_data, resolution=resolution)
	snap.pl.umap(merged_data, color='leiden', show=False, out_file=umap_mnc_cluster_filename, height=height)
	if precelltype is not None:
		snap.pl.umap(merged_data, color=var, show=False, out_file=umap_mnc_precelltype_filename, height=height)

if harmony:
	snap.pp.harmony(merged_data, batch="sample", max_iter_harmony=20)
	#use_rep="X_spectral_harmony"
	snap.tl.umap(merged_data, use_rep="X_spectral_harmony")
	snap.pl.umap(merged_data, color="sample", show=False, out_file=umap_harmony_sample_filename, height=height)
	snap.pp.knn(merged_data, n_neighbors=n_neighbors, use_rep="X_spectral_harmony")
	snap.tl.leiden(merged_data, resolution=resolution)
	snap.pl.umap(merged_data, color='leiden', show=False, out_file=umap_harmony_cluster_filename, height=height)
	if precelltype is not None:
		snap.pl.umap(merged_data, color=var, show=False, out_file=umap_harmony_precelltype_filename, height=height)


merged_data.close()

