#!/bin/python

import argparse

parser = argparse.ArgumentParser(description="call peak for annset h5ads file")
parser.add_argument("-i", "--input", type=str, dest="input", help="input annset h5ads file")
parser.add_argument("-c", "--precelltype", type=str, dest="precelltype", default=None, help="input previous celltype csv file")
parser.add_argument("-v", "--var", type=str, default="celltype", dest="var", help="which var to use")
parser.add_argument("-q", "--qvalue", type=float, default=0.01, dest="qvalue", help="qvalue")
parser.add_argument("-b", "--blacklist", type=str, default=None, dest="blacklist", help="blacklist")
parser.add_argument("-g", "--genome", type=str, default="mm10", dest="genome", help="genome")
parser.add_argument("-u", "--pvalue", type=float, default=0.01, dest="pvalue", help="pvalue")
parser.add_argument("-f", "--bmat_features", type=int, default=500000, dest="bmat_features", help="bmat number of features")
parser.add_argument("-m", "--pmat_features", type=int, default=1000, dest="pmat_features", help="pmat number of features")
parser.add_argument("-k", "--peak", help="call peak", action="store_true", dest="peak")
parser.add_argument("-d", "--roughdiff", help="rough diff peak", action="store_true", dest="roughdiff")
parser.add_argument("-r", "--regressdiff", help="regression diff peak", action="store_true", dest="regressdiff")
parser.add_argument("-l", "--umap", help="plot umap", action="store_true", dest="umap")
parser.add_argument("-o", "--out", type=str, dest="output", help="prefix of output file")

args = parser.parse_args()

import snapatac2 as snap
snap.__version__ 
from pathlib import Path
import numpy as np
import pandas as pd
import sys
import polars as pl
import os
import anndata

# Input files
TF_path = args.input
file = Path(TF_path)
precelltype = args.precelltype
var = args.var
peak_qvalue = args.qvalue
blacklist_path = args.blacklist
genome = args.genome
diff_pvalue = args.pvalue
bmat_features = args.bmat_features
pmat_features = args.pmat_features

peak_call = args.peak
rough_diff = args.roughdiff
regress_diff = args.regressdiff
umap_plot = args.umap
out_prefix = args.output


# out_files
peak_filename = out_prefix + ".peaks.union.csv"
pmat_filename = out_prefix + ".pmat.h5ad"
bmat_umap_filename = out_prefix + ".bmat.umap.pdf"
pmat_umap_filename = out_prefix + ".pmat.umap.pdf"
count_filename = out_prefix + ".count_peakByCelltype.csv"
cpm_filename = out_prefix + ".cpm_peakByCelltype.csv"
heatmap_filename = out_prefix + ".rough.diffpeaks.pdf"


def call_peak(file=file, precelltype=precelltype, var=var, qvalue=peak_qvalue, blacklist_path=blacklist_path, genome=genome, peak_filename=peak_filename, pmat_filename=pmat_filename):
    TF_data = snap.read_dataset(file)

    if precelltype is not None:
        df = pd.read_csv(precelltype)
        #df['barcode'] = df['barcode'].str.replace('_',':')
        if var in df.columns:
            TF_data_celltype = [df.loc[df['barcode']==cell, var].to_list()[0] for cell in TF_data.obs_names]
            TF_data.obs[var] = TF_data_celltype
        else:
            sys.exit("Var is not support!")

    snap.tl.macs3(TF_data, groupby=var, qvalue=qvalue, shift=-75, extsize=150, blacklist=blacklist_path)

    if genome == "mm10":
        chrom = snap.genome.mm10
    elif genome == "hg38":
        chrom = snap.genome.hg38
    else:
        sys.exit("Only mm10 and hg38 are supported for genome!")

    peaks = snap.tl.merge_peaks(TF_data.uns['macs3'], chrom_sizes=chrom)
    peaks_df = peaks.to_pandas()
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True) 
    peaks_df.to_csv(peak_filename, sep=',', header=True, index=False)
    peak_mat = snap.pp.make_peak_matrix(TF_data, use_rep=peaks['Peaks'])
    peak_mat.write_h5ad(pmat_filename)
    TF_data.close()
    #peak_mat.close()

    # get cpm by celltype for nmf
    peak_mat: anndata.AnnData = anndata.read_h5ad(pmat_filename)
    npeak = peak_mat.shape[1]
    ann_meta = peak_mat.obs
    grouped = ann_meta.groupby(var)
    nsc = grouped.ngroups
    cnt_pbysc = pd.DataFrame(
        np.zeros( (npeak, nsc) , dtype = np.float64),
        columns = list(grouped.groups.keys()),
        index = peak_mat.var_names)

    for group, idx in grouped.indices.items():
        cnt_pbysc[group] = np.ravel(
            peak_mat.X[idx, ].sum(axis = 0, dtype = np.float64))

    cnt_pbysc.to_csv(
        count_filename,
        index = True,
        header = True
    )

    sum_pofsc = cnt_pbysc.sum(axis = 0)
    values = cnt_pbysc.values
    factors = sum_pofsc.values

    cpm_pbysc = (values / factors ) * 10e6 # bug: should be 1e6
    cpm_pbysc_df = pd.DataFrame(cpm_pbysc, columns = cnt_pbysc.columns)
    cpm_pbysc_df.index = cnt_pbysc.index
    cpm_pbysc_df.to_csv(
        cpm_filename,
        index = True,
        header = True   
    )



# rough diff peaks
def rough_diffpeak(pmat_filename=pmat_filename, var=var, diff_pvalue=diff_pvalue, heatmap_filename=heatmap_filename):
    TF_peak_mat = snap.read(pmat_filename, backed="r")
    TF_marker_peaks = snap.tl.marker_regions(TF_peak_mat, groupby=var, pvalue=diff_pvalue)
    snap.pl.regions(TF_peak_mat, groupby=var, peaks=TF_marker_peaks, interactive=False, out_file=heatmap_filename)

    diff_path = os.path.join(os.path.dirname(out_prefix), "diffpeaks_rough")
    os.makedirs(diff_path, exist_ok=True) 
    for celltype,peaks in TF_marker_peaks.items():
        diffpeak_file = os.path.join(diff_path, ".".join([celltype, "diffpeak.txt"]))
        peaks.to_frame().to_csv(diffpeak_file, header=False, index=False)
    TF_peak_mat.close()


# regression diff peaks
def regress_diffpeak(pmat_filename=pmat_filename, peak_filename=peak_filename, var=var, diff_pvalue=diff_pvalue):
    TF_peak_mat = snap.read(pmat_filename, backed="r")
    diff_path = os.path.join(os.path.dirname(out_prefix), "diffpeaks_regression")
    os.makedirs(diff_path, exist_ok=True) 
    #pd.value_counts(TF_peak_mat.obs[var])
    TF_peaks = pd.read_csv(peak_filename)
    for celltype in set(TF_peak_mat.obs[var]):
        print(celltype)
        barcodes = np.array(TF_peak_mat.obs_names)
        background = []
        for i in np.unique(TF_peak_mat.obs[var]):
            if i != celltype:
                cells = np.random.choice(barcodes[TF_peak_mat.obs[var] == i], size=100, replace=False)
                background.append(cells)
        celltype_cells = TF_peak_mat.obs[var] == celltype
        background = np.concatenate(background)
        diff_peaks = snap.tl.diff_test(
            TF_peak_mat,
            cell_group1=celltype_cells,
            cell_group2=background,
            features=TF_peaks[celltype].to_numpy(),
            direction="positive",
            min_log_fc = 0.1,
            min_pct = 0.05,
        )
        if diff_peaks.shape[0]>0:
            #diff_peaks = diff_peaks.filter(pl.col('adjusted p-value') < 0.1)
            diff_peaks = diff_peaks.filter(pl.col('p-value') < diff_pvalue)
            diffpeak_file = os.path.join(diff_path, ".".join([celltype, "diffpeak.txt"]))
            diff_peaks.to_pandas().to_csv(diffpeak_file, header=False, index=False)
            heatmap_file = os.path.join(diff_path, ".".join([celltype, "diffpeaks.pdf"]))
            snap.pl.regions(
                TF_peak_mat,
                groupby=var,
                peaks={ celltype: diff_peaks['feature name'].to_numpy() },
                interactive=False,
                out_file=heatmap_file
            )

def plot_umap(file=file, pmat_filename=pmat_filename, precelltype=precelltype, var=var, bmat_umap_filename=bmat_umap_filename, pmat_umap_filename=pmat_umap_filename):
    TF_data = snap.read_dataset(file)

    if precelltype is not None:
        df = pd.read_csv(precelltype)
        #df['barcode'] = df['barcode'].str.replace('_',':')
        if var in df.columns:
            TF_data_celltype = [df.loc[df['barcode']==cell, var].to_list()[0] for cell in TF_data.obs_names]
            TF_data.obs[var] = TF_data_celltype
        else:
            sys.exit("Var is not support!")
    # bmat umap
    snap.pp.select_features(TF_data, n_features=bmat_features)
    snap.tl.spectral(TF_data)
    snap.tl.umap(TF_data)
    snap.pl.umap(TF_data, color=var, show=False, out_file=bmat_umap_filename, height=500)
    TF_data.close()
    # pmat umap
    peak_mat = snap.read(pmat_filename)
    snap.pp.select_features(peak_mat, n_features=pmat_features)
    snap.tl.spectral(peak_mat)
    snap.tl.umap(peak_mat)
    peak_mat_celltype = [df.loc[df['barcode']==cell, var].to_list()[0] for cell in peak_mat.obs_names]
    peak_mat.obs[var] = peak_mat_celltype
    snap.pl.umap(peak_mat, color=var, show=False, out_file=pmat_umap_filename, height=500)
    peak_mat.close()


def main():
    if peak_call:
        call_peak()
    if rough_diff:
        rough_diffpeak()
    if regress_diff:
        regress_diffpeak()
    if umap_plot:
        plot_umap()

if __name__ == "__main__":
    main()
