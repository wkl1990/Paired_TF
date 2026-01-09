#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="doublet removal for qc h5ad file")
parser.add_argument("-i", "--input", type=str, dest="input", help="input h5ad file")
parser.add_argument("-b", "--bin_size", type=int, default=500, dest="bin_size", help="bin size")
parser.add_argument("-f", "--n_features", type=int, default=500000, dest="n_features", help="number of features")
parser.add_argument("-l", "--blacklist", type=str, default=None, dest="blacklist", help="blacklist")
parser.add_argument("-w", "--whitelist", type=str, default=None, dest="whitelist", help="whitelist")
#parser.add_argument("-r", "--rate", type=float, default=0.1, dest="rate", help="expected doublet rate")
#parser.add_argument("-p", "--probability", type=float, default=0.5, dest="probability", help="probability threshold")
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
bin_size = args.bin_size
n_features = args.n_features
blacklist_path = args.blacklist
whitelist_path = args.whitelist
#rate = args.rate
#probability = args.probability
out_prefix = args.output
file = Path(qc_path)

# qc 
h5ad_filename = out_prefix + ".h5ad"

qc_data = snap.read(file, backed='r')
doublet_data = qc_data.copy(h5ad_filename)
qc_data.close()

snap.pp.add_tile_matrix(doublet_data, bin_size=bin_size)
if blacklist_path != None and whitelist_path != None:
	snap.pp.select_features(doublet_data, n_features=n_features, blacklist=blacklist_path, whitelist=whitelist_path)
elif blacklist_path == None and whitelist_path != None:
	snap.pp.select_features(doublet_data, n_features=n_features, whitelist=whitelist_path)
elif blacklist_path != None and whitelist_path == None:
	snap.pp.select_features(doublet_data, n_features=n_features, blacklist=blacklist_path)
else:
	snap.pp.select_features(doublet_data, n_features=n_features)

#snap.pp.scrublet(doublet_data, expected_doublet_rate=rate)
#snap.pp.filter_doublets(doublet_data, probability_threshold=probability)

ATAC_cells = np.array(doublet_data.obs_names)
cells_filename = out_prefix + ".cells.txt"
np.savetxt(cells_filename, ATAC_cells, fmt="%s", delimiter=',')
qc_data.close()
doublet_data.close()

