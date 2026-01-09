#!/bin/python

# example: python scripts/Python/03a.DNA.merge_data.py -i samle_info_celltypeage.txt -n PairedTF1.h5ad,PairedTF2.h5ad -o subset/

import argparse

parser = argparse.ArgumentParser(description="extract cell data from each sample")
parser.add_argument("-i", "--input", type=str, dest="input", required=True, help="input sample files")
parser.add_argument("-n", "--name", type=str, dest="name", required=True, help="input sample names")
#parser.add_argument("-c", "--cell", type=str, dest="cell", required=True, help="input cell info. file")
#parser.add_argument("-p", "--processes", type=int, dest="processes", help="number of processes to use")
parser.add_argument("-o", "--output", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
print(snap.__version__)
import numpy as np
#from multiprocessing import Pool
import sys


def main():
	inputs = args.input
	names = args.name
#	cell_file = args.cell
#	processes = args.processes
	out_prefix = args.output

	sample_files = inputs.split(',')
	sample_names = names.split(',')
	if len(sample_names) != len(sample_files):
		print("Error: length of names and files do not match!")
		sys.exit(0)
	dataset = dict(zip(sample_names, sample_files))

	hdfiles_bk = [(sample, snap.read(dataset[sample], backed='r')) for sample in dataset.keys()] 
	output_path = "".join([out_prefix, ".h5ads"])
	merge_dataset = snap.AnnDataSet(
		adatas=hdfiles_bk,
		filename=output_path,
		add_key='sample',
	)
	merge_dataset.obs_names = merge_dataset.obs['sample'].to_numpy() + ":" + merge_dataset.obs_names

	print("Subset dataset written to file:", output_path)

	sys.exit(0)

if __name__ == "__main__":
	main()

