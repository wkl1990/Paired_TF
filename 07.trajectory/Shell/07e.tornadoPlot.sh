path=/tscc/projects/PairedTF

# plot ATAC
computeMatrix reference-point \
 -S /tscc/projects/CEMBA2/326_OPC_NN.bw /tscc/projects/CEMBA2//327_Oligo_NN.bw \
 -R ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_pos.bed ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_neg.bed \
 --referencePoint center \
 -a 250 -b 250 \ 
 -out ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_ATAConly.tab.gz

plotHeatmap -m ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_ATAConly.tab.gz -out ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_ATAConly.pdf --samplesLabel "ATAC_OPC" "ATAC_OGC" \
 --regionsLabel "ATAC_pos" "ATAC_neg" --outFileSortedRegions ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_ATAConly.sort.bed

# plot CTCF
computeMatrix reference-point \
 -S ${path}/OGCvsOPC/bigwig/CTCF_OPC_treat_pileup.srt.bw ${path}/OGCvsOPC/bigwig/CTCF_OGC_treat_pileup.srt.bw \
 -R ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_ATAConly.sort.bed \
 --referencePoint center \
 -a 250 -b 250 \ 
 -out ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_CTCF.tab.gz

plotHeatmap -m ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_CTCF.tab.gz -out ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_CTCF.pdf --samplesLabel "CTCF_OPC" "CTCF_OGC" \
 --regionsLabel "ATAC_pos" "ATAC_neg" --sortRegions no

# plot SOX2
computeMatrix reference-point \
 -S ${path}/biwig/SOX2_OPC.bw ${path}/biwig/SOX2_OGC.bw \
 -R ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_ATAConly.sort.bed \
 --referencePoint center \
 -a 250 -b 250 \ 
 -out ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_SOX2.tab.gz

plotHeatmap -m ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_SOX2.tab.gz -out ${path}/OGCvsOPC/tornado/ATAC_log2FC_diffpeak_SOX2.pdf --samplesLabel "SOX2_OPC" "SOX2_OGC" \
 --regionsLabel "ATAC_pos" "ATAC_neg" --sortRegions no


