import numpy as np
import cv2
import seaborn
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import cooltools
import cooler
import cooltools.lib.plotting
from cooltools import numutils
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LogNorm, Normalize, NoNorm, AsinhNorm, FuncNorm, PowerNorm, SymLogNorm, TwoSlopeNorm
from scipy import ndimage as nd
from skimage.transform import rescale, resize
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# funcitons
def plot_hic_real(chr, start, end, resolution, cool_file, out_prefix, vmin=0, vmax=1, clip=False):
	chromp = chr
	anchor_pos1 = [start, end]
	res = resolution
	#nempt = 5
	#vvmax = 0.1
	anchor1 = [int(ii / res) for ii in anchor_pos1]
	Qchr_subset = np.concatenate([np.arange(anchor1[0], anchor1[1])])
	Qannotate = [ii / 1000000 for ii in anchor_pos1]
	Qtitle = chromp + ':' + '-'.join([str(element) for element in anchor_pos1]) 

	clr = cooler.Cooler(cool_file + '::resolutions/' + str(res))
	Q = clr.matrix(balance=False).fetch(chromp)
	#Q = Q - np.diag(np.diag(Q))

	### insert white space
	Qnew = Q[Qchr_subset][:, Qchr_subset]
	vvmax = np.quantile(Qnew, vmax)
	vvmin = np.quantile(Qnew, vmin)
	if vvmin<=0:
		print("Warning: vmin is less than 0")
		vvmin = None

	Qrotate = nd.rotate(Qnew, 45, order=0, reshape=True, prefilter=False, cval=0)
	pos_list_new = [0, (anchor1[1] - anchor1[0])]
	new_pos_list = [round(np.sqrt(2)*ii) for ii in pos_list_new]

	f, axes = plt.subplots(1, 1, figsize=(9, 9))

	ax = axes#[idx]
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	im = ax.matshow(Qrotate, cmap='fall', norm=LogNorm(vmax=vvmax, vmin=vvmin, clip=clip))

	h = len(Qrotate)
	ax.set_ylim([0.5*h, 0])
	ax.set_xlim([0, h])
	ax.set_yticks([])
	ax.set_yticklabels([])
	ax.set_xticks([])
	ax.set_xticks(new_pos_list)
	ax.set_xticklabels(Qannotate)
	ax.xaxis.set_ticks_position('bottom')
	ax.set_title(Qtitle)

	plt.colorbar(im, fraction=0.01, pad=0.04, label='counts')#Balanced counts')
	plt.tight_layout()
	plt.savefig(out_prefix + f"_nodiag_quantile1_{chromp}_{anchor_pos1[0]}_{anchor_pos1[1]}_res{res}.pdf")
	plt.close('all')

def hic_cv2_resize(npy_file, resolution):
	# Load the .npy file
	matrix = np.load(npy_file)
	# Calculate the binning factor
	factor = int(resolution) / 8192  # 10kb/8kb = 1.25, so approximately group every 5 bins into 4
	# Ensure dimensions are divisible by binning factor
	new_size = int(matrix.shape[0] // factor) #* int(factor)  # Adjust size to be divisible
	# Resize the image
	resized_matrix = cv2.resize(matrix, (new_size, new_size))
	return resized_matrix

def plot_hic_pred(chr, start, end, resolution, npy_mtx, out_prefix, vmin=0, vmax=1, clip=False):
	chromp = chr
	anchor_pos1 = [start, end]
	res = resolution
	#nempt = 5
	#vvmax = 0.1
	anchor1 = [int(ii / res) for ii in anchor_pos1]
	Qchr_subset = np.concatenate([np.arange(anchor1[0], anchor1[1])])
	Qannotate = [format(ii / 1000000, '.2f') for ii in anchor_pos1]
	Qtitle = chromp + ':' + '-'.join([str(element) for element in anchor_pos1]) 

	Qnew = npy_mtx
	vvmax = np.quantile(Qnew, vmax)
	vvmin = np.quantile(Qnew, vmin)
	if vvmin<=0:
		print("Warning: vmin is less than 0")
		vvmin = None

	Qrotate = nd.rotate(Qnew, 45, order=0, reshape=True, prefilter=False, cval=0)
	pos_list_new = [0, (anchor1[1] - anchor1[0])]
	new_pos_list = [round(np.sqrt(2)*ii) for ii in pos_list_new]

	f, axes = plt.subplots(1, 1, figsize=(9, 9))

	ax = axes#[idx]
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	im = ax.matshow(Qrotate, cmap='fall', norm=LogNorm(vmax=vvmax, vmin=vvmin, clip=clip))

	h = len(Qrotate)
	ax.set_ylim([0.5*h, 0])
	ax.set_xlim([0, h])
	ax.set_yticks([])
	ax.set_yticklabels([])
	ax.set_xticks([])
	ax.set_xticks(new_pos_list)
	ax.set_xticklabels(Qannotate)
	ax.xaxis.set_ticks_position('bottom')
	ax.set_title(Qtitle)

	plt.colorbar(im, fraction=0.01, pad=0.04, label='counts')#Balanced counts')
	plt.tight_layout()
	plt.savefig(out_prefix + f"_{chromp}_{anchor_pos1[0]}_{anchor_pos1[1]}_res{res}.pdf")
	plt.close('all')

def get_hic_real(cool_file, chr, start, end, resolution):
	chromp = chr
	anchor_pos1 = [start, end]
	res = resolution
	anchor = [int(ii / res) for ii in anchor_pos1]
	Qchr_subset = np.concatenate([np.arange(anchor[0], anchor[1])])

	clr = cooler.Cooler(cool_file + '::resolutions/' + str(res))
	Q = clr.matrix(balance=False).fetch(chromp)
	#Q = Q - np.diag(np.diag(Q))

	### insert white space
	Q_subset = Q[Qchr_subset][:, Qchr_subset]
	return Q_subset

def cal_similarity(matrix1, matrix2, method='pearson'):
	# Flatten the matrices
	matrix1_flat = matrix1.flatten()
	matrix2_flat = matrix2.flatten()
	# Remove NaN values
	mask = ~np.isnan(matrix1_flat) & ~np.isnan(matrix2_flat)
	matrix1_flat = matrix1_flat[mask]
	matrix2_flat = matrix2_flat[mask]
	# Calculate similarity
	if method == 'pearson':
		correlation, _ = pearsonr(matrix1_flat, matrix2_flat)
		return correlation
	elif method == 'spearman':
		correlation, _ = spearmanr(matrix1_flat, matrix2_flat)
		return correlation

# plot real hic
marker_regions_pd = pd.read_csv("/tscc/projects/PairedTF/hic_pred/marker_regions.txt", sep="\t", header=0)
celltype_subclass_pd = pd.read_csv("/tscc/projects/PairedTF/CTCF/celltype_subclass.txt", sep="\t", header=0)

for i in range(marker_regions_pd.shape[0]):
	for j in range(celltype_subclass_pd.shape[0]):
		chr = marker_regions_pd['chr'][i]
		start = marker_regions_pd['region1'][i]
		end = int(start+2100000)
		cool_file = "/tscc/projects/PairedTF/TAD/mouse-m3c/Q.10K/" + celltype_subclass_pd['subclass1'][j] + ".Q.10K.mcool"
		plot_dir = os.path.join("/tscc/projects/PairedTF/hic_pred/outputs/marker_regions", "plot", celltype_subclass_pd['celltype1'][j])
		os.makedirs(plot_dir, exist_ok=True)
		plot_hic_real(chr, start, end, 10000, cool_file, plot_dir + "/Real")

# plot predicted hic
for i in range(marker_regions_pd.shape[0]):
	for j in range(celltype_subclass_pd.shape[0]):
		chr = marker_regions_pd['chr'][i]
		start = marker_regions_pd['region1'][i]
		end = int(start+2097152)
		bam2bw_file = "/tscc/projects/PairedTF/hic_pred/outputs/marker_regions/" + celltype_subclass_pd['celltype1'][j] + "_bam2bw/prediction/npy/" + chr + "_" + str(start) + ".npy"
		npy_bam2bw = np.load(bam2bw_file)
		npy_bam2bw_cv10k = hic_cv2_resize(bam2bw_file, 10000)
		zeroCTCF_file = "/tscc/projects/PairedTF/hic_pred/outputs/marker_regions/" + celltype_subclass_pd['celltype1'][j] + "_zeroCTCF/prediction/npy/" + chr + "_" + str(start) + ".npy"
		npy_zerobw = np.load(zeroCTCF_file)
		npy_zerobw_cv10k = hic_cv2_resize(zeroCTCF_file, 10000)
		plot_bam2bw = os.path.join("/tscc/projects/PairedTF/hic_pred/outputs/marker_regions", "plot", celltype_subclass_pd['celltype1'][j] + "_bam2bw")
		os.makedirs(plot_bam2bw, exist_ok=True)
		plot_hic_pred(chr, start, end, 10000, npy_bam2bw_cv10k, plot_bam2bw + "/bam2bw_cv10K")
		plot_zerobw = os.path.join("/tscc/projects/PairedTF/hic_pred/outputs/marker_regions", "plot", celltype_subclass_pd['celltype1'][j] + "_zerobw")
		os.makedirs(plot_zerobw, exist_ok=True)
		plot_hic_pred(chr, start, end, 10000, npy_zerobw_cv10k, plot_zerobw + "/zerobw_cv10K")

# get predicted hic matrices
for i in range(marker_regions_pd.shape[0]):
	for j in range(celltype_subclass_pd.shape[0]):
		chr = marker_regions_pd['chr'][i]
		start = marker_regions_pd['region1'][i]
		real_end = int(start + 2100000)
		pred_end = int(start + 2097152)
		subclass = celltype_subclass_pd['subclass1'][j]
		celltype = celltype_subclass_pd['celltype1'][j]
		cool_file = "/tscc/projects/PairedTF/TAD/mouse-m3c/Q.10K/" + subclass + ".Q.10K.mcool"
		mtx_dir = os.path.join("/tscc/projects/PairedTF/hic_pred/outputs/marker_regions", "mtx", celltype)
		os.makedirs(mtx_dir, exist_ok=True)
		mtx_real_10K = get_hic_real(cool_file, chr, start, real_end, 10000)
		np.savetxt(os.path.join(mtx_dir, f"Real_{chr}_{start}_{real_end}_res10K.txt"), mtx_real_10K, delimiter="\t")

		bam2bw_file = "/tscc/projects/PairedTF/hic_pred/outputs/marker_regions/" + celltype + "_bam2bw/prediction/npy/" + chr + "_" + str(start) + ".npy"
		npy_bam2bw = np.load(bam2bw_file)
		npy_bam2bw_cv10k = hic_cv2_resize(bam2bw_file, 10000)
		zeroCTCF_file = "/tscc/projects/PairedTF/hic_pred/outputs/marker_regions/" + celltype + "_zeroCTCF/prediction/npy/" + chr + "_" + str(start) + ".npy"
		npy_zerobw = np.load(zeroCTCF_file)
		npy_zerobw_cv10k = hic_cv2_resize(zeroCTCF_file, 10000)

		mtx_bam2bw = os.path.join("/tscc/projects/PairedTF/hic_pred/outputs/marker_regions", "mtx", celltype + "_bam2bw")
		os.makedirs(mtx_bam2bw, exist_ok=True)
		np.savetxt(os.path.join(mtx_bam2bw, f"bam2bw_cv10K_{chr}_{start}_{pred_end}.txt"), npy_bam2bw_cv10k, delimiter="\t")

		mtx_zerobw = os.path.join("/tscc/projects/PairedTF/hic_pred/outputs/marker_regions", "mtx", celltype + "_zerobw")
		os.makedirs(mtx_zerobw, exist_ok=True)
		np.savetxt(os.path.join(mtx_zerobw, f"zerobw_cv10K_{chr}_{start}_{pred_end}.txt"), npy_zerobw_cv10k, delimiter="\t")



