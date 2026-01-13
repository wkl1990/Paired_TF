export PERL5LIB=~/bin/lib/site_perl/
export PATH=$PATH:~/bin/bin

path=/tscc/projects/PairedTF
# merge peaks
cat /tscc/projects/CEMBA2/OPC_NN.bed /tscc/projects/CEMBA2/Oligo_NN.bed | sort -k1,1 -k2,2n | uniq > ${path}/OGCvsOPC/ATAC_OPC_Oligo_NN.bed
# get mtx
multiBamSummary BED-file --BED ${path}/OGCvsOPC/ATAC_OPC_Oligo_NN.bed --bamfiles /tscc/projects/CEMBA2/326_OPC_NN.srt.bam /tscc/projects/CEMBA2/327_Oligo_NN.srt.bam -p 2 --smartLabels -o ${path}/OGCvsOPC/ATAC_OPC_OGC.npz --outRawCounts ${path}/OGCvsOPC/ATAC_OPC_OGC.tab --scalingFactors ${path}/OGCvsOPC/ATAC_OPC_OGC.txt &
# diffreps
nohup ~/bin/bin/diffReps.pl -tr /tscc/projects/CEMBA2/327_Oligo_NN.bed -co /tscc/projects/CEMBA2/326_OPC_NN.bed -ch /tscc/Reference/genome/mm10/mm10.chrom.sizes -re ${path}/OGCvsOPC/ATAC_OGC_OPC_diffreps.gt.txt -me gt 2>&1 > ${path}/OGCvsOPC/ATAC_OGC_OPC_diffreps.gt.log &

cat ${path}/OGCvsOPC/ATAC_OGC_OPC_diffreps.gt.txt | grep "Up" > ${path}/OGCvsOPC/ATAC_OGC_OPC_diffreps.gt.pos.bed
cat ${path}/OGCvsOPC/ATAC_OGC_OPC_diffreps.gt.txt | grep "Down" > ${path}/OGCvsOPC/ATAC_OGC_OPC_diffreps.gt.neg.bed

# annotate diff peaks
for diff in pos neg
do  
diff_file=${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.bed
bedtools intersect -wao -f 0.5 -a ${diff_file} -b /tscc/annot/gencode.vM23.gene.tssUpDn1k.bed >${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.tssUpDn1k.bed
awk 'BEGIN{FS=OFS="\t"}{if($6!="."){print $1,$2,$3,$5,"proximal",$12}else{print $1,$2,$3,$5,"distal","nan"}}' ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.tssUpDn1k.bed | sort -k1,1 -k2,2n | uniq > ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.bed
awk 'BEGIN{FS=OFS="\t"}{print $1"_"$2"_"$3,$4,$5,$6}' ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.bed | sort -k1,1 > ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.peaks
cat ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.tssUpDn1k.proximal.distal.peaks | sed '1i peak\tlogFC\ttype\tgene' | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}$3=="proximal"{print $0}' > ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_${diff}.proximal.peaks
done
# motif analysis with homer
ls ${path}/OGCvsOPC/ATAC_peak_cpm_log2FC_*.bed | while read id 
do
file=$(basename $id)
name=${file%.*}
echo $name
nohup findMotifsGenome.pl ${id} mm10 ${path}/OGCvsOPC/motif/${name}_given_ATACbg/ -bg /tscc/projects/CEMBA2/mba.whole.sa2.final.peak.srt.bed -p 8 -size given 2>&1 > ${path}/OGCvsOPC/motif/${name}.homer.motif_given_ATACbg.log &
sleep 3
#loadavg
done

