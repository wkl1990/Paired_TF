#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="prepare fragment file")
parser.add_argument("-i", "--input", type=str, dest="input", help="input bam file")
parser.add_argument("-q", "--mapq", type=int, default=30, dest="mapq", help="min mapq value")
parser.add_argument("-p", "--paired", action='store_true', dest="paired", help="paired-end reads")
parser.add_argument("-f", "--field", type=str, dest="field", help="tag or regex to indicate parameter")
parser.add_argument("-b", "--barcode", type=str, dest="barcode", help="barcode tag or regex")
parser.add_argument("-u", "--umi", type=str, dest="umi", default=None, help="umi tag or regex")
parser.add_argument("-s", "--stat", action='store_true', dest="stat", help="output stat")
parser.add_argument("-o", "--out", type=str, dest="output", help="output fragment file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ 
from pathlib import Path
import numpy as np
import pandas as pd
import os

# Input files
bam = args.input
mapq = args.mapq
paired = args.paired
field = args.field
barcode = args.barcode
umi = args.umi
stat = args.stat
output_file = args.output
base_name=os.path.splitext(output_file)[0]
stat_file=base_name+".stat"
bam_file = Path(bam)

# load data
if umi==None:
	if field=="tag":
		if paired:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=True, barcode_tag=barcode)
		else:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=False, barcode_tag=barcode)
	elif field=="regex":
		if paired:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=True, barcode_regex=barcode)
		else:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=False, barcode_regex=barcode)
else:
	if field=="tag":
		if paired:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=True, barcode_tag=barcode, umi_tag=umi)
		else:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=False, barcode_tag=barcode, umi_tag=umi)
	elif field=="regex":
		if paired:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=True, barcode_regex=barcode, umi_regex=umi)
		else:
			fragment_stat=snap.pp.make_fragment_file(bam_file, output_file, min_mapq=mapq, is_paired=False, barcode_regex=barcode, umi_regex=umi)

if stat:
	with open(stat_file, 'w') as f:
		f.writelines(str(fragment_stat))

print("Job is done!")


